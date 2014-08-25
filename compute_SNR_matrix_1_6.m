function [SNR_mat,Shift_mat] = compute_SNR_matrix_1_6(InData,OutInfo,LarmorFreq,dwelltime,Phantom_flag,ControlParameters)
%
% compute_SNR_matrix_x_y Compute the SNR for each voxel of CSI data.
%
% This function was written by Bernhard Strasser, 2011 - 2012.
%
%
% Compute the SNR for each voxel of a CSI-matrix. The matrix must have the size (ROW,COL,SLC,vecSize) and it must be phased. One can also pass an Array of processed spectra
% to the function (this must have size [ROW COL SLC]), e.g. the spectra outputted by LCModel in the ".Coord" files. 
% The function will write an error-log-file and plots concerning the Metabo_seeking process and the Noise_seeking process if the struct field OutInfo.out_dir is provided.
% As output one gets an SNR and a Shift array of size [ROW COL SLC].
% 
%
% [SNR_mat,Shift_mat] = compute_SNR_matrix_x_y(InData,OutInfo,LarmorFreq,dwelltime,Phantom_flag,ControlParameters)
%
% Input:
% InData:               Structure, that can have the following fields:
%                       -   .csi: Phased csi-data. size(InData.csi) = [ROW COL SLC vecSize].
%                       -   .LCMspec (optional): Struct containing fields .ppm and .Spectra. Spectra & Chemical Shift of the same data as csi, but already somehow processed, 
%                           e.g. the Baseline was subtracted. Example: The spectra from the LCModel .COORD files. size(InData.LCMspec) = [ROW COL SLC]. If this field is passed over,
%                           the Signal is computed with LCMspec, but the noise computed with csi.
%                       -   .mask (optional): The mask, determining for which voxels in [ROW COL SLC] the SNR should be computed. size(InData.mask) = [ROW COL SLC].
%                       -   .Linewidth (optional): If you have the linewidths of the spectra at hand (e.g. computed by LCModel) this can provide good information for calculating the SeekMetabo.BasisRegionPeakdistance_CS
%                           for each spectrum individually. In this case, you should also provide ControlParameters.SeekMetabo.LinewidthToBasisDistanceFactor and ControlParameters.SeekMetabo.ClipBasisRegionPeakdistance and you
%                           can omit ControlParameters.SeekMetabo.BaisRegionPeakdistance_CS, since it doesn't have any influence.
%
% OutInfo:              Structure with infos about the output data. Containing the following fields:
%                       -   .out_dir: The output directory for writing the error-log file and the plots. If not existant, don't write these files.
%                       -   .quiet: If true, suppress display output (e.g. "Processing Voxel x = 32, y = 32, z = 1").
%                       -   .print_individual_spectra_flag: If true, save the plots documentating the SNR computation, like the metabolite peak, the noise, etc.
% 
% LarmorFreq:           Larmor frequency
% dwelltime:            Dwelltime
% Phantom_flag:         Passed over to seek_peak_x_y function, which uses different settings if a phantom was used. (Use 1,0 or true,false)
%
% ControlParameters:    This variable can either be a structure or a string referring to a file, with info about how to process the data. In the latter case, the file gets executed and overwrites all those variables, 
%                       that are in this file. If it is a structure, it can contain the following fields, and also overwrites all variables defined in cell "0. COVENTIONS, PREPARATIONS, DEFINITIONS":
% MANDATORY:
% Nothing
%
% OPTIONAL, RECCOMENDED:
%                       -   NoiseRegion_CS: Cell, containing 1x2 Arrays which determine the regions in ppm where the Noise should be extracted from. E.g. {[5.7 7.3],[7.9 8.4]}.
%                       -   zerofilling: Factor determining how much the InData.csi should be zerofilled for computing the Signal.
%                       -   SeekMetabo: Structure containing info about Seeking Metabo, containing fields:
%                           *   .PeakRegionWidth_IfRefFound: In which area around the ideal chemical shift - Ref-shift should be searched for the peak if Ref was found
%                           *   .PeakRegionWidth_IfRefNotFound: In which area around the ideal chemical shift - Ref-shift should be searched for the peak if Ref was not found
%                           *   .BasisRegionPeakdistance_CS: The distance of the regions for searching the "Basis points" to the peak in ppm. 
%                                Used to determine where to search for the "basis" of the peak. size(SeekMetabo.BasisRegionPeakdistance_CS) = [1 2] (2 for HiFi and LoFi)
%                           *   .LinewidthToBasisDistanceFactor: If you provide the linewidth of each spectrum, this info is used for defining the distance of the basisregions to the peak. 
%                               This factor determines how the linewidth is translated for calculating the BasisRegionPeakdistance_CS.
%                           *   .ClipBasisRegionPeakdistance: All Peakdistances that are higher than this value get clipped to that value.
%                           *   .InRealOrAbs: If 'real' Search peak in real part of spectrum, if 'abs' in absolute values.
%                           *   .SNRThreshold: Threshold for considering peak as Metabo.
%                       -   IdealRef_CS: The ideal chemical shift of the reference, e.g. if water is the reference, this would be 4.7 (or 4.65 for phantoms?)
%                       -   IdealMetabo_CS: The ideal chemical shift of the peak for which the SNR should be computed, e.g. if NAA should be used, this would be 2.01
% OPTIONAL
%                       - Any other variable defined in cell '0. CONVENTIONS, PREPARATIONS, DEFINITIONS'. The passed-over variable is then used instead of the defined one.
%
%
% Output
% SNR_mat:              Array containing the SNR infos
% Shift_mat:            Array containing the shift in ppm of the whole spectrum. Computed using the Reference Peak. Can be considered as B0-map.
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: seek_peak_x_y, noise_computation_new_x_y 

% Further remarks:





%% 0. COVENTIONS, PREPARATIONS, DEFINITIONS

% CONVENTIONS

% CS = Chemical Shift, so the value val_CS is a Chemical Shift
% SP = Spectral Point(s), so the value val_SP is not a Chemical Shift, but the number of points inside the spectrum, counting: left to right (e.g. if you have a vecSize = 1024, then the last point (LoFi) would have 1024 
% SS = Spectral Signal, so the value val_SS is a Signal of the spectrum; this has nothing to do with nazis!
% HiFi = HighField or Upfield, i.e. the Chemical Shift at high field
% LoFi = Low Field or Downfield
% If vectors of CS or SP are given, the first entry is always the HiFi (low SP), the second the LoFi (high SP)
%
% Compute Spectral Point out of Chemical Shift
% Example:  NoiseRegion_SP(2) = find(min(abs(CS_vec - NoiseRegion_CS(2))) == abs(CS_vec - NoiseRegion_CS(2)));
% What is done here? 
% 1) Compute the differences of the CS_vector (= all possible CS-values) and the given CS-value (CS_vec - NoiseRegion_CS(2))
% 2) Find the value inside the CS_vector that is closest to the given CS_value (min(abs(CS_vec - NoiseRegion_CS(2))))
% 3) Create a array where the position of this closest point is given in the CS_vector (min(abs(...)) == abs(...)) by a logical 1 (all other entries of the array are zero)
% 4) Find the position of the 1




% ROW, COL, SLC, vecSize
if(numel(size(InData.csi)) > 4)
    display([char(10) 'Your input matrix has too many dimensions. The program needs matrix of size ROW x COL x SLC x vecSize. Aborting.' char(10)])
    return
end
[ROW COL SLC vecSize] = size(InData.csi);


% VARIABLE STANDARD ASSIGNMENTS
if(~isfield(InData,'mask'))
    InData.mask = ones([ROW COL SLC]);
end
if(~isfield(OutInfo,'quiet'))
    OutInfo.quiet = false;
end
if(~isfield(OutInfo,'print_individual_spectra_flag'))
    OutInfo.print_individual_spectra_flag = false;
end
if(~exist('Phantom_flag','var'))
    Phantom_flag = false;
end



% DEFINITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                            THESE ARE STANDARD VALUES.                                                               %%%%%
%%%%%   DON'T CHANGE THESE VALUES IF YOU NEED DIFFERENT ADJUSTMENTS, BUT INSTEAD DEFINE THESE VARIABLES IN "ControlParameters" ACCORDING TO YOUR NEEDS.   %%%%%
%%%%%               ONLY CHANGE THE VALUES, IF YOU THINK DIFFERENT STANDARD VALUES ARE APPROPRIATE (BECAUSE THESE ARE TOO SPECIAL E.G.).                  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Variables for seeking Ref
IdealRef_CS = 4.7;                                      % Chemical Shift of Ref in ideal case  
SeekRef_PeakRegion_DistToIdeal_CS = [0.17 0.17];        % Search for the reference peak from Ideal_Ref_CS - SeekRef_PeakRegion_DistToIdeal_CS(1) to Ideal_Ref_CS + SeekRef_PeakRegion_DistToIdeal_CS(1)
SeekRef_PeakRegion_StepSize_CS = 0.04;                  % Stepsize with which the PeakRegion gets increased in each loop, if the reference was not found.
SeekRef_TotalSteps = 1;                                 % Total Steps for increasing the PeakRegion    
SeekRef_SNRthreshold = 15.0;                            % Threshold for SNR in order to consider peak as Ref peak
SeekRef_BasisRegionWidth_CS = 0.08;                     % The Basis values will be searched in a certain region. This value defines how large this region should be in ppm (chemshift)
SeekRef_BasisRegionPeakdistance_HiFi_CS = 0.1;          % Compute the "Basis"-value in the HiFi region compared to peak at a Chemical Shift distance of this value
SeekRef_BasisRegionPeakdistance_LoFi_CS = 0.1;          % Same for LoFi


% Variables for seeking Metabo
IdealMetabo_CS = 2.01;                                  % Chemical Shift of Metabo in ideal case
SeekMetabo.PeakRegionWidth_IfRefFound = 0.2;            % Seek for Metabo in the region (IdealMetabo_CS - Chemshift of RefPeak +- this value) If RefPeak was found
SeekMetabo.PeakRegionWidth_IfRefNotFound = 0.2;         % Same but if RefPeak was not found
SeekMetabo.BasisRegionWidth_CS = 0.05;                  % The Basis values will be searched in a certain region. This value defines how large this region should be in ppm (chemshift)
SeekMetabo.BasisRegionPeakdistance_CS = [0.16 0.16];    % Compute the "Basis"-value in the HiFi region compared to peak at a Chemical Shift distance of this value
SeekMetabo.LinewidthToBasisDistanceFactor = [9/3, 9/3]; % How to derive the distance of the Basis regions to the peak from the Linewidth.
SeekMetabo.ClipBasisRegionPeakdistance = 0.16;
SeekMetabo.SNRthreshold = 3.5;                          % Threshold for SNR in order to consider peak as Metabo peak
SeekMetabo.SimilarBasisSignals_Ratio = 9/5;             % This value says: If the difference between the two basis-region-points in signal is this value times the difference between the higher basis-region-point-signal and the
                                                        % found peak, the peak is not considered as a peak (the basis signals are too far away from each other).
SeekMetabo.InRealOrAbs = 'real';


% Variables for Noise computation
NoiseRegion_CS = {[6 7.5]};                             % Region where to compute the standard deviation of the noise


% Variables for plotting the Metabolite Region
MetaboRegion_HiFi_CS = 4.0;
MetaboRegion_LoFi_CS = 1.0;


% Zerofilling
zerofilling = 2;


% OVERWRITE VARIABLES GIVEN BY ControlParameters
if(ischar(ControlParameters) && exist(ControlParameters,'file'))  
    
    % Run the file with path ControlParameters if it is a file.
    eval(['run ' ControlParameters]);
    
elseif(isstruct(ControlParameters))
    
    % Create the variables by the name of the fieldnames of ControlParameters. E.g.: ControlParameters.IdealRef_CS = 3. Then the variable IdealRef_CS is defined with IdealRef_CS = 3;
    struct_fieldnames = fieldnames(ControlParameters);
    for field_no = 1:numel(struct_fieldnames); 
        eval([struct_fieldnames{field_no} ' = ControlParameters.' struct_fieldnames{field_no} ';']);
    end

end


if(numel(SeekMetabo.LinewidthToBasisDistanceFactor)<2)
    SeekMetabo.LinewidthToBasisDistanceFactor(2) = SeekMetabo.LinewidthToBasisDistanceFactor;
end






%% 1. Compute Preparatory Stuff


% Compute BasisRegionPeakdistance_CS if linewidth is provided
if(isfield(InData, 'Linewidth'))
    SeekMetabo.BasisRegionPeakdistance_CS = zeros([ROW COL SLC 2]);
    SeekMetabo.BasisRegionPeakdistance_CS(:,:,:,1) = SeekMetabo.LinewidthToBasisDistanceFactor(1) * InData.Linewidth;                                          % The left one
    SeekMetabo.BasisRegionPeakdistance_CS(:,:,:,2) = SeekMetabo.LinewidthToBasisDistanceFactor(2) * InData.Linewidth;                                          % The right one 
    % Clip Peakdistance, because otherwise the basis region is much too far away (e.g. at 2.4 ppm) at some voxels.
    SeekMetabo.BasisRegionPeakdistance_CS(SeekMetabo.BasisRegionPeakdistance_CS > SeekMetabo.ClipBasisRegionPeakdistance) = SeekMetabo.ClipBasisRegionPeakdistance;  
end


% Create error log file
if(isfield(OutInfo,'out_dir') && ~strcmpi(OutInfo.out_dir,'No'))
    mkdir(OutInfo.out_dir);
    mkdir([OutInfo.out_dir '/failed']);
    mkdir([OutInfo.out_dir '/succeeded']);    
    error_log_file = sprintf('%s/compute_SNR_error_log.txt',OutInfo.out_dir);
    error_log_fid = fopen(error_log_file, 'w');
    fprintf(error_log_fid,'Error Log File\nThese errors, sorted by x,y,z of matrix were found:\n\n');
end

% initialization
SNR_mat = zeros(ROW,COL,SLC);
Shift_mat = zeros(ROW,COL,SLC);
CurrentVoxel = 1;


% Total Voxels to Process
TotalVoxels = numel(find(InData.mask));


% Reference-Peak Region
SeekRef_PeakRegion_LoFi_CS = IdealRef_CS - SeekRef_PeakRegion_DistToIdeal_CS(1);       % initial LoFi (low chemical shift, low frequency) chemshift region to search for Ref in ppm
SeekRef_PeakRegion_HiFi_CS = IdealRef_CS + SeekRef_PeakRegion_DistToIdeal_CS(2);       % initial HiFi region


% chemshift vector
CS_vec = compute_chemshift_vector_1_1(LarmorFreq,dwelltime/10^9,vecSize);                   % Compute the Chemical Shift: For each point of spectrum this function gives the corresponding chemical shift
CS_vec_zf = compute_chemshift_vector_1_1(LarmorFreq,dwelltime/10^9,vecSize*zerofilling);    % dwelltime gets not increased with zero_filling: zeroes get just added at the END of vector


% ZeroFilling and FFT
csi_zf = zeros(ROW,COL,SLC,zerofilling*vecSize);
csi_zf(:,:,:,1:vecSize) = InData.csi;

csi_mat_freq_zf = fftshift(fft(csi_zf, [],4),4);
csi_mat_freq = fftshift(fft(InData.csi, [],4),4);


% Define the spectra and chemical shift for searching Metabo. 
if(isfield(InData,'LCMspec'))
    SeekMetabo_specmat = InData.LCMspec.Spectra;
    SeekMetabo_specmat_zf = InData.LCMspec.Spectra;
    SeekMetabo_CS = InData.LCMspec.ppm;
    SeekMetabo_CS_zf = InData.LCMspec.ppm;             % No oversampling implemented for LCModel spectra.
else
    SeekMetabo_specmat = csi_mat_freq;
    SeekMetabo_specmat_zf = csi_mat_freq_zf;    
    SeekMetabo_CS = CS_vec;
    SeekMetabo_CS_zf = CS_vec_zf;
end


% Convert SeekRef CS --> SP
% CS_vec_zf(1) - SeekRef_BasisRegionWidth_CS defines the distance to the first CS-point. Find out how many points this CS-distance corresponds to.
% Since the CS-scale is linear, this point distance is the same at any CS. -1 at the end because if SeekRef.BasisRegionPeakdistance_CS = 0, then otherwise we would still get 1
SeekRef_BasisRegionWidth_SP = find(min(abs(CS_vec_zf(1) - SeekRef_BasisRegionWidth_CS - CS_vec_zf)) == abs(CS_vec_zf(1) - SeekRef_BasisRegionWidth_CS - CS_vec_zf));
SeekRef_BasisRegionPeakdistance_HiFi_SP = find(min(abs(CS_vec_zf(1) - SeekRef_BasisRegionPeakdistance_HiFi_CS - CS_vec_zf)) == abs(CS_vec_zf(1) - SeekRef_BasisRegionPeakdistance_HiFi_CS - CS_vec_zf)) - 1;
SeekRef_BasisRegionPeakdistance_LoFi_SP = find(min(abs(CS_vec_zf(1) - SeekRef_BasisRegionPeakdistance_LoFi_CS - CS_vec_zf)) == abs(CS_vec_zf(1) - SeekRef_BasisRegionPeakdistance_LoFi_CS - CS_vec_zf)) - 1;

% Convert SeekMetabo CS --> SP
SeekMetabo_BasisRegionWidth_SP = find(min(abs(SeekMetabo_CS_zf(1) - SeekMetabo.BasisRegionWidth_CS - SeekMetabo_CS_zf)) == abs(SeekMetabo_CS_zf(1) - SeekMetabo.BasisRegionWidth_CS - SeekMetabo_CS_zf));
if(numel(SeekMetabo.BasisRegionPeakdistance_CS) <= 2)
    SeekMetabo_BasisRegionPeakdistance_HiFi_SP = find(min(abs(SeekMetabo_CS_zf(1) - SeekMetabo.BasisRegionPeakdistance_CS(1) - SeekMetabo_CS_zf)) == abs(SeekMetabo_CS_zf(1) - SeekMetabo.BasisRegionPeakdistance_CS(1) - SeekMetabo_CS_zf)) - 1;
    SeekMetabo_BasisRegionPeakdistance_LoFi_SP = find(min(abs(SeekMetabo_CS_zf(1) - SeekMetabo.BasisRegionPeakdistance_CS(2) - SeekMetabo_CS_zf)) == abs(SeekMetabo_CS_zf(1) - SeekMetabo.BasisRegionPeakdistance_CS(2) - SeekMetabo_CS_zf)) - 1;
end


% Convert NoiseRegions CS --> SP
NoOfNoiseRegions = size(NoiseRegion_CS,2);
NoiseRegion_SP = cell(1,NoOfNoiseRegions);
NoiseRegion_SP_vec = cell(1,NoOfNoiseRegions);
for NoiseIndex = 1:NoOfNoiseRegions
    NoiseRegion_CS{NoiseIndex} = sort(NoiseRegion_CS{NoiseIndex},'descend');
    NoiseRegion_SP{NoiseIndex}(1) = find(min(abs(CS_vec - NoiseRegion_CS{NoiseIndex}(1))) == abs(CS_vec - NoiseRegion_CS{NoiseIndex}(1)));
    NoiseRegion_SP{NoiseIndex}(2) = find(min(abs(CS_vec - NoiseRegion_CS{NoiseIndex}(2))) == abs(CS_vec - NoiseRegion_CS{NoiseIndex}(2)));
    NoiseRegion_SP_vec{NoiseIndex} = NoiseRegion_SP{NoiseIndex}(1) : NoiseRegion_SP{NoiseIndex}(2);
end
NoiseSubplot_Row = round(sqrt(NoOfNoiseRegions));
NoiseSubplot_Col = ceil(sqrt(NoOfNoiseRegions));


% Convert MetaboRegion CS --> SP
MetaboRegion_SP = [find(min(abs(SeekMetabo_CS - MetaboRegion_HiFi_CS)) == abs(SeekMetabo_CS - MetaboRegion_HiFi_CS)), find(min(abs(SeekMetabo_CS - MetaboRegion_LoFi_CS)) == abs(SeekMetabo_CS - MetaboRegion_LoFi_CS))];



%% 1.D1 DEBUG MODE: How close are the CS_vec values to the SeekRef_PeakRegion_LoFi_CS and HiFi?; is the zerofilling enough?

% diff_SeekRefLoFi_to_chemshiftvec = CS_vec_zf(min(abs(CS_vec_zf - SeekRef_PeakRegion_LoFi_CS)) == abs(CS_vec_zf - SeekRef_PeakRegion_LoFi_CS)) - SeekRef_PeakRegion_LoFi_CS
% diff_SeekRefHiFi_to_chemshiftvec = CS_vec_zf(min(abs(CS_vec_zf - SeekRef_PeakRegion_HiFi_CS)) == abs(CS_vec_zf - SeekRef_PeakRegion_HiFi_CS)) - SeekRef_PeakRegion_HiFi_CS
% clear diff_SeekRefLoFi_to_chemshiftvec diff_SeekRefHiFi_to_chemshiftvec




%% 1.D2 DEBUG MODE2: Apodize the FID to see what is in the Noise & Metabolite region


% csi_apod = transpose(squeeze(InData.csi(22,16,SLC,:))) .* exp(-(1:vecSize)*0.00);
% csi_mat_freq_apod = fftshift(fft(csi_apod));
% % csi_mat_freq = fftshift(fft(squeeze(InData.csi(14,10,SLC,:))));
% % figure
% % plot(CS_vec(MetaboRegion_SP(1):MetaboRegion_SP(2)),real(csi_mat_freq_apod(MetaboRegion_SP(1):MetaboRegion_SP(2))) ,'r','Linewidth',1.9)
% % hold on
% % plot(CS_vec(MetaboRegion_SP(1):MetaboRegion_SP(2)),real(csi_mat_freq(MetaboRegion_SP(1):MetaboRegion_SP(2))))
% % hold off
% set(gca,'XDir','reverse');
% figure
% plot(CS_vec(NoiseRegion_SP_vec),real(csi_mat_freq_apod(NoiseRegion_SP_vec)))
% set(gca,'XDir','reverse');
% pause



%% 1.D3 DEBUG MODE3: Mask
% InData.mask = zeros([64 64 4]);
% % InData.mask(28,10) = 1;
% % InData.mask(24,14) = 1;
% InData.mask(23,23,1:4) = 1;



%% 2. Process each voxel individually: LOOPS

for x = 1:ROW
	for y = 1:COL
        for z = 1:SLC

           
            
            
            %% 4. Voxel-specific Preparations
            
            if(InData.mask(x,y,z) == 0)
				continue
            end

            if(~OutInfo.quiet)
                display([ char(10) 'Processing Voxel x = ' num2str(x) ', y = ' num2str(y) ', z = ' num2str(z) ', Voxel ' num2str(CurrentVoxel) ' of ' num2str(TotalVoxels) ' . . .'])  
            end
            
            
			spectrum_zf = transpose(squeeze(csi_mat_freq_zf(x,y,z,:)));
            spectrum = transpose(squeeze(csi_mat_freq(x,y,z,:)));
            spectrum_real = real(spectrum);
            spectrum_zf_abs = abs(spectrum_zf); %#ok<NASGU>
            spectrum_zf_real = real(spectrum_zf); %#ok<NASGU>
            SeekMetabo_spec_zf = transpose(squeeze(SeekMetabo_specmat_zf(x,y,z,:)));
            SeekMetabo_spec = transpose(squeeze(SeekMetabo_specmat(x,y,z,:)));            
            
            if(numel(find(isnan(spectrum))) > 0)
                display('Problem occurred: The spectrum contains NaNs.')             
                if(isfield(OutInfo,'out_dir'))  
                    fprintf(error_log_fid,'z_%d_y_%d_x_%d\n',z,y,x);     
                    fprintf(error_log_fid,'E2 Error: The Spectrum contained NaNs.\n');
                end
                CurrentVoxel = CurrentVoxel + 1;
                continue
            end
            
            
            % Define the distance where to search for the Basis of the peak if for each voxel this is given seperately (computed with the linewidth)
            if(numel(SeekMetabo.BasisRegionPeakdistance_CS) > 2)
                
                if(isnan(SeekMetabo.BasisRegionPeakdistance_CS(x,y,z,1)) || SeekMetabo.BasisRegionPeakdistance_CS(x,y,z,1) <= 0)
                    SeekMetabo.BasisRegionPeakdistance_CS(x,y,z,:) = [0.15 0.15];
                end
                
                SeekMetabo_BasisRegionPeakdistance_HiFi_SP = find(min(abs(SeekMetabo_CS_zf(1) - SeekMetabo.BasisRegionPeakdistance_CS(x,y,z,1) - SeekMetabo_CS_zf)) == abs(SeekMetabo_CS_zf(1) - SeekMetabo.BasisRegionPeakdistance_CS(x,y,z,1) - SeekMetabo_CS_zf)) - 1;
                SeekMetabo_BasisRegionPeakdistance_LoFi_SP = find(min(abs(SeekMetabo_CS_zf(1) - SeekMetabo.BasisRegionPeakdistance_CS(x,y,z,2) - SeekMetabo_CS_zf)) == abs(SeekMetabo_CS_zf(1) - SeekMetabo.BasisRegionPeakdistance_CS(x,y,z,2) - SeekMetabo_CS_zf)) - 1;             
                
            end
            
            

            
            %% 5. Compute the standard deviation of the noise for searching the Ref peak

            [Noise,NoiseSignal_Smoothie,NoiseSignal_SubSmooth,Noise_ExcludeSubregions_StartEndPoints] = noise_computation_new_1_3(spectrum,NoiseRegion_SP,'Smooth');            
            
            
            
			%% 6. Search for Ref_peak

			% Define searching region		
			
			% initialization
			NoRefPeak_found = 1;
			SeekRef_Step = 0;
            SeekRef_InRealOrAbs = 'abs';
            SeekRef_PeakRegion_LoFi_CS_dummy = SeekRef_PeakRegion_LoFi_CS;
            SeekRef_PeakRegion_HiFi_CS_dummy = SeekRef_PeakRegion_HiFi_CS;


			% while peak not found increase region to look for peak and check again
			while(NoRefPeak_found && (SeekRef_Step <= SeekRef_TotalSteps-1))

                % convert chemshift to points
				% The exact chemshift is most probable not within the CS_vec, so search for the minimum
				SeekRef_PeakRegion_LoFi_SP = find(min(abs(CS_vec_zf - SeekRef_PeakRegion_LoFi_CS_dummy)) == abs(CS_vec_zf - SeekRef_PeakRegion_LoFi_CS_dummy));
				SeekRef_PeakRegion_HiFi_SP = find(min(abs(CS_vec_zf - SeekRef_PeakRegion_HiFi_CS_dummy)) == abs(CS_vec_zf - SeekRef_PeakRegion_HiFi_CS_dummy));
                SeekRef_PeakRegion_SP = [SeekRef_PeakRegion_HiFi_SP, SeekRef_PeakRegion_LoFi_SP];

         
				% Search for RefPeak
                [NoRefPeak_found,FoundRefPeak_SP,FoundRefPeak_SS,FoundRef_Basis_SP,FoundRef_Basis_SS,FoundRef_MeanBasis_SS,SeekRef_HiFiBasisRegion_SP,SeekRef_LoFiBasisRegion_SP] = ...
                 seek_peak_1_3(spectrum_zf,SeekRef_PeakRegion_SP,SeekRef_InRealOrAbs,0,Phantom_flag,'SNR_Criterion', ...
                 {{Noise,SeekRef_SNRthreshold,SeekMetabo.SimilarBasisSignals_Ratio,SeekRef_BasisRegionWidth_SP,SeekRef_BasisRegionPeakdistance_HiFi_SP,SeekRef_BasisRegionPeakdistance_LoFi_SP}});


				% If peak not found: Set Searching Routine to search in abs value of spectrum, if this has already be done: Enlarge seek-region
                if(NoRefPeak_found)
                    if(strcmpi(SeekRef_InRealOrAbs,'real'))
                        SeekRef_InRealOrAbs = 'abs';
                    else
                        SeekRef_PeakRegion_LoFi_CS_dummy = SeekRef_PeakRegion_LoFi_CS_dummy - SeekRef_PeakRegion_StepSize_CS;
                        SeekRef_PeakRegion_HiFi_CS_dummy = SeekRef_PeakRegion_HiFi_CS_dummy + SeekRef_PeakRegion_StepSize_CS;
                        SeekRef_Step = SeekRef_Step + 1;
                        if(SeekRef_Step <= SeekRef_TotalSteps-1)
                            SeekRef_InRealOrAbs = 'real';              % start to search again in real part.
                        end
                    end
                end

			end %while

            if(~NoRefPeak_found)
                Shift = CS_vec_zf(FoundRefPeak_SP) - IdealRef_CS;
            end

            
            
            
            %% 7. Compute Standard Deviation of Noise again, this time at the right (shifted) CS, computed with the shift of the Ref peak 
            
            if(~NoRefPeak_found)
                NoiseRegion_CS_ShiftCorr = NoiseRegion_CS;
                for NoiseIndex = 1:NoOfNoiseRegions
                    NoiseRegion_CS_ShiftCorr{NoiseIndex} = NoiseRegion_CS{NoiseIndex} + Shift;
                    NoiseRegion_SP{NoiseIndex}(1) = find(min(abs(CS_vec - NoiseRegion_CS_ShiftCorr{NoiseIndex}(1))) == abs(CS_vec - NoiseRegion_CS_ShiftCorr{NoiseIndex}(1)));
                    NoiseRegion_SP{NoiseIndex}(2) = find(min(abs(CS_vec - NoiseRegion_CS_ShiftCorr{NoiseIndex}(2))) == abs(CS_vec - NoiseRegion_CS_ShiftCorr{NoiseIndex}(2)));
                    NoiseRegion_SP_vec{NoiseIndex} = NoiseRegion_SP{NoiseIndex}(1) : NoiseRegion_SP{NoiseIndex}(2);
                end

                [Noise,NoiseSignal_Smoothie,NoiseSignal_SubSmooth,Noise_ExcludeSubregions_StartEndPoints] = noise_computation_new_1_3(spectrum,NoiseRegion_SP,'Smooth');         
            end
            
            
            
            
            
			%% 8. Seek & Destroy Metabo


			% 8.0 Define Searching Region

            SeekMetabo_PeakRegion_LoFi_CS = IdealMetabo_CS;
            
            if(~isfield(InData,'LCMspec') && ~NoRefPeak_found)                              % Shift CS search region if no LCMspectrum is used and the Ref peak found.
                SeekMetabo_PeakRegion_LoFi_CS = SeekMetabo_PeakRegion_LoFi_CS + Shift;                        % IdealRef_CS - CS_vec_zf(FoundRefPeak_SP) gives the Shift of the CS
            end                                                                                                                 % so IdealMetabo_CS - Shift gives the shifted Metabo_peak  

            if(NoRefPeak_found)
                SeekMetabo_PeakRegion_HiFi_CS = SeekMetabo_PeakRegion_LoFi_CS + SeekMetabo.PeakRegionWidth_IfRefNotFound/2;                
                SeekMetabo_PeakRegion_LoFi_CS = SeekMetabo_PeakRegion_LoFi_CS - SeekMetabo.PeakRegionWidth_IfRefNotFound/2;
            else
                SeekMetabo_PeakRegion_HiFi_CS = SeekMetabo_PeakRegion_LoFi_CS + SeekMetabo.PeakRegionWidth_IfRefFound/2;         
                SeekMetabo_PeakRegion_LoFi_CS = SeekMetabo_PeakRegion_LoFi_CS - SeekMetabo.PeakRegionWidth_IfRefFound/2;
            end

            
            % convert chemshift to point
            SeekMetabo_PeakRegion_LoFi_SP = find(min(abs(SeekMetabo_CS_zf - SeekMetabo_PeakRegion_LoFi_CS)) == abs(SeekMetabo_CS_zf - SeekMetabo_PeakRegion_LoFi_CS));
            SeekMetabo_PeakRegion_HiFi_SP = find(min(abs(SeekMetabo_CS_zf - SeekMetabo_PeakRegion_HiFi_CS)) == abs(SeekMetabo_CS_zf - SeekMetabo_PeakRegion_HiFi_CS));             
            SeekMetabo_PeakRegion_SP = [SeekMetabo_PeakRegion_HiFi_SP, SeekMetabo_PeakRegion_LoFi_SP];

        

            % 8.1 Search for Metabo, compute SNR
            [NoMetabopeak_found,FoundMetaboPeak_SP,FoundMetaboPeak_SS,FoundMetabo_Basis_SP,FoundMetabo_Basis_SS,FoundMetabo_MeanBasis_SS,SeekMetabo_HiFiBasisRegion_SP,SeekMetabo_LoFiBasisRegion_SP,SeekMetabo_Polyfit,SeekMetabo_SubBas] = ...
            seek_peak_1_3(SeekMetabo_spec_zf,SeekMetabo_PeakRegion_SP,SeekMetabo.InRealOrAbs,0,Phantom_flag,'SNR_Criterion', ...
            {{Noise,SeekMetabo.SNRthreshold,SeekMetabo.SimilarBasisSignals_Ratio,SeekMetabo_BasisRegionWidth_SP,SeekMetabo_BasisRegionPeakdistance_HiFi_SP, SeekMetabo_BasisRegionPeakdistance_LoFi_SP}});            
                       
         
          

            %% 8. Save SNR & Shift to Matrices
            
            if(~NoMetabopeak_found)
                SNR_mat(x,y,z) = (FoundMetaboPeak_SS - FoundMetabo_MeanBasis_SS)/(2*Noise);
            else
                SNR_mat(x,y,z) = NaN;
            end
            if(~NoRefPeak_found)
            	Shift_mat(x,y,z) = Shift;
            else
            	Shift_mat(x,y,z) = NaN; 
            end



         
            
            %% 9. Display Messages 
            
            if(~OutInfo.quiet)
                % RefPeak
                if(~NoRefPeak_found)
                    display(['Ref found at ', num2str(round(CS_vec_zf(FoundRefPeak_SP)*1000)/1000), ' ppm. Shift = ', num2str(round(Shift*1000)/1000) ' ppm.'])
                else
                    display('No Ref found')
                end

                % Metabopeak
                if(NoMetabopeak_found)
                    display('Seek ''n destroy Metabo   FAILED!')
                else
                    display('Seek ''n destroy Metabo   SUCCEEDED!')
                end
            end
            
            
            
            %% 8. Write error log file
            
            if(isfield(OutInfo,'out_dir') && (NoRefPeak_found || NoMetabopeak_found))  
                % Write Voxel in error log file
                fprintf(error_log_fid,'z_%d_y_%d_x_%d\n',z,y,x);              
                if(NoRefPeak_found)
                	fprintf(error_log_fid,'W1 WARNING: No RefPeak found. Assuming Ref suppression.\n');
                end
                
                if(NoMetabopeak_found)
                	fprintf(error_log_fid,'E1 ERROR: No Metabo-peak found. Move to next voxel.\n');
                    
%                     if(error2_no_Metabobasis_HiFi) 
%                     	fprintf(error_log_fid,'E2 ERROR: Basis of Metabo-peak at high chemical shift side not found. Move to next voxel.\n');
%                     end
%                     
%                     if(error3_no_Metabobasis_LoFi) 
%                         fprintf(error_log_fid,'E3 ERROR: Basis of Metabo-peak at low chemical shift side not found. Move to next voxel.\n');
%                     end
                    
                end   
            end	
            
            
            
            

			%% 9. plot: 
            %            - Noise Region
            %            - Noise Region SubPeak
            %            - Metabolite spectrum
            %            - Metabo Region
            % TURNED OFF - whole spectrum zerofilling
            %            - Ref Region

            % Don't try to understand all this crap if you don't need to. 

            if(isfield(OutInfo,'out_dir') && ~strcmpi(OutInfo.out_dir,'No') && OutInfo.print_individual_spectra_flag)


                
                % 9.0 Preparations
                
                % Define out directory for plots
                if(~NoMetabopeak_found)
                    out_dir_FailedOrWon = [OutInfo.out_dir '/succeeded'];
                else
                    out_dir_FailedOrWon = [OutInfo.out_dir '/failed'];
                end
                
                % Make Spectrum Real
                SeekMetabo_spec_zf = eval([SeekMetabo.InRealOrAbs '(SeekMetabo_spec_zf);']);
                SeekMetabo_spec = eval([SeekMetabo.InRealOrAbs '(SeekMetabo_spec);']);

                    
                %% 9.1 Noise Region
                
                NoisePeak_fig = figure('visible','off');             
                
                for NoiseIndex = 1:NoOfNoiseRegions
                    subplot(NoiseSubplot_Col,NoiseSubplot_Row,NoiseIndex);
                    plot(CS_vec(NoiseRegion_SP_vec{NoiseIndex}),spectrum_real(NoiseRegion_SP_vec{NoiseIndex}))
                    hold on
                    %plot(Plot_NoiseRegion_x,polyval(Polyfit_Noise_Real,NoiseRegion_SP_vec),'r')
                    plot(CS_vec(NoiseRegion_SP_vec{NoiseIndex}),real(NoiseSignal_Smoothie{NoiseIndex}(NoiseRegion_SP_vec{NoiseIndex})),'r','Linewidth',0.3)   
                    set(gca,'XDir','reverse');
                    hold off
                    if(NoiseIndex ==ceil(NoOfNoiseRegions/2))
                        title(sprintf('Noise Peak subtracted, voxel x_%d y_%d z_%d', x,y,z),'Interpreter','none')
                    end
                    xlabel('Chemical Shift')
                    %ylabel('Signal')        
                end
                
                % Save & Close
                %saveas(NoisePeak_fig,sprintf('%s/x_%s_y_%s_z_%s_NoiseAndPeak', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)),'fig')
                saveas(NoisePeak_fig,sprintf('%s/x_%s_y_%s_z_%s_NoiseAndPeak', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)),'epsc2')
                close(NoisePeak_fig)

                
                

                %% 9.2 NoiseSubPeak
                 
                NoiseSubPeak_fig = figure('visible','off');                
                
                for NoiseIndex = 1:NoOfNoiseRegions

                    subplot(NoiseSubplot_Col,NoiseSubplot_Row,NoiseIndex);
                    plot(CS_vec(NoiseRegion_SP_vec{NoiseIndex}),real(NoiseSignal_SubSmooth{NoiseIndex}(NoiseRegion_SP_vec{NoiseIndex})))
                    hold on
                    % plot X-es in Regions that are excluded
                    for ExcludeSubregions_LoopIndex = 1:size(Noise_ExcludeSubregions_StartEndPoints{NoiseIndex},2)
                        StartPoint = Noise_ExcludeSubregions_StartEndPoints{NoiseIndex}{ExcludeSubregions_LoopIndex}(1);
                        EndPoint = Noise_ExcludeSubregions_StartEndPoints{NoiseIndex}{ExcludeSubregions_LoopIndex}(2);
                        StartCS = CS_vec(StartPoint);
                        EndCS = CS_vec(EndPoint);
                        plot([StartCS,EndCS],[min(real(NoiseSignal_SubSmooth{NoiseIndex}(StartPoint:EndPoint))),max(real(NoiseSignal_SubSmooth{NoiseIndex}(StartPoint:EndPoint)))],'r','Linewidth',1.4) % plot line increasing diagonal from left to right
                        plot([StartCS,EndCS],[max(real(NoiseSignal_SubSmooth{NoiseIndex}(StartPoint:EndPoint))),min(real(NoiseSignal_SubSmooth{NoiseIndex}(StartPoint:EndPoint)))],'r','Linewidth',1.4) % plot line decreasing diagonal from left to right
                    end  
                    hold off
                    set(gca,'XDir','reverse');
                    if(NoiseIndex == floor(NoiseSubplot_Col/2))
                        title(sprintf('Noise Peak subtracted, voxel x_%d y_%d z_%d', x,y,z),'Interpreter','none')
                    end
                        
                    xlabel('Chemical Shift')
                    %ylabel('Signal')
                end

                % Save & Close
                %saveas(NoiseSubPeak_fig,sprintf('%s/x_%s_y_%s_z_%s_NoiseSubPeak', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)),'fig')
                saveas(NoiseSubPeak_fig,sprintf('%s/x_%s_y_%s_z_%s_NoiseSubPeak', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)),'epsc2') 
                close(NoiseSubPeak_fig)


                
                
                
                %% 9.3 Metabolite Spectrum
                
                
                Metabo_spectrum_fig = figure('visible','off');
               
                plot(SeekMetabo_CS(MetaboRegion_SP(1):MetaboRegion_SP(2)),SeekMetabo_spec(MetaboRegion_SP(1):MetaboRegion_SP(2)))         
                hold on
                plot(SeekMetabo_CS_zf(FoundMetaboPeak_SP),FoundMetaboPeak_SS,'o','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','k');               % Where peak was computed
                plot(SeekMetabo_CS_zf(FoundMetabo_Basis_SP(1)),FoundMetabo_Basis_SS(1),'o','MarkerSize',6,'MarkerFaceColor','c','MarkerEdgeColor','k');     % Where HiFi_Basis was computed 
                plot(SeekMetabo_CS_zf(FoundMetabo_Basis_SP(2)),FoundMetabo_Basis_SS(2),'o','MarkerSize',6,'MarkerFaceColor','b','MarkerEdgeColor','k');     % Where LoFi_Basis was computed
                hold off
                legend('Spectrum','Metabolite Peak','Metabolite Basis Le','Metabolite Basis Ri')

                set(gca,'XDir','reverse');
                title(sprintf('Metabolite Spectrum, voxel x_%d y_%d z_%d', x,y,z), 'Interpreter', 'none')
                xlabel('Chemical Shift')
                ylabel('Signal')


                % Save & Close
                %saveas(Metabo_spectrum_fig,sprintf('%s/x_%s_y_%s_z_%s_Metabo_spectrum.fig', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))
                saveas(Metabo_spectrum_fig,sprintf('%s/x_%s_y_%s_z_%s_Metabo_spectrum', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)),'epsc2')                 
                close(Metabo_spectrum_fig)
                
                
                
                
                
                %% 9.4 Metabo REGION
                
                
                % Definitions
                
                % Plotting Region
                Plot_Metaboregion_HiFi_SP = min([SeekMetabo_PeakRegion_SP(1), SeekMetabo_HiFiBasisRegion_SP(1)]);            % If the PeakRegion is farther left then plot spectrum from PeakRegion
                Plot_Metaboregion_LoFi_SP = max([SeekMetabo_PeakRegion_SP(2), SeekMetabo_LoFiBasisRegion_SP(2)]);            % if Basis Region farther left then Basis Region; Same for right point of spectrum

                
                % SPECTRUM ITSELF
                Plot_SeekMetabo_Spec_x = SeekMetabo_CS_zf(Plot_Metaboregion_HiFi_SP:Plot_Metaboregion_LoFi_SP);
                Plot_SeekMetabo_Spec_y = SeekMetabo_spec_zf(Plot_Metaboregion_HiFi_SP:Plot_Metaboregion_LoFi_SP);


                % Vertical Lines for Peak & Basis Regions
                Plot_SeekMetabo_VertLine_Min = min(Plot_SeekMetabo_Spec_y);
                if(Plot_SeekMetabo_VertLine_Min >= 0)
                   Plot_SeekMetabo_VertLine_Min = -max(Plot_SeekMetabo_Spec_y)/20;
                end            


                % PEAK REGION VERTICAL LINE
                Plot_SeekMetabo_PeakRegion_VertLine_Max = max([SeekMetabo_spec_zf(SeekMetabo_PeakRegion_SP(1)), SeekMetabo_spec_zf(SeekMetabo_PeakRegion_SP(2))]);
                Plot_SeekMetabo_PeakRegion_VertLine_y = [Plot_SeekMetabo_VertLine_Min,Plot_SeekMetabo_PeakRegion_VertLine_Max];
                Plot_SeekMetabo_PeakRegion_HiFi_VertLine_x = repmat(SeekMetabo_CS_zf(SeekMetabo_PeakRegion_SP(1)),[1 2]);                % if you use just plot(x,vector), x scalar, only dots get plotted
                Plot_SeekMetabo_PeakRegion_LoFi_VertLine_x = repmat(SeekMetabo_CS_zf(SeekMetabo_PeakRegion_SP(2)),[1 2]);              % So plot([x,x],[y_start,y_end])    

                % HiFi BASIS REGION VERTICAL LINE
                Plot_SeekMetabo_HiFiBasisRegion_VertLine_Max = max([SeekMetabo_spec_zf(SeekMetabo_HiFiBasisRegion_SP(1)), SeekMetabo_spec_zf(SeekMetabo_HiFiBasisRegion_SP(2))]);
                Plot_SeekMetabo_HiFiBasisRegion_VertLine_y = [Plot_SeekMetabo_VertLine_Min,Plot_SeekMetabo_HiFiBasisRegion_VertLine_Max];
                Plot_SeekMetabo_HiFiBasisRegion_HiFi_VertLine_x = repmat(SeekMetabo_CS_zf(SeekMetabo_HiFiBasisRegion_SP(1)),[1 2]);
                Plot_SeekMetabo_HiFiBasisRegion_LoFi_VertLine_x = repmat(SeekMetabo_CS_zf(SeekMetabo_HiFiBasisRegion_SP(2)), [1 2]);

                % LoFi BASIS REGION VERTICAL LINE
                Plot_SeekMetabo_LoFiBasisRegion_VertLine_Max = max([SeekMetabo_spec_zf(SeekMetabo_LoFiBasisRegion_SP(1)), SeekMetabo_spec_zf(SeekMetabo_LoFiBasisRegion_SP(2))]);            
                Plot_SeekMetabo_LoFiBasisRegion_VertLine_y = [Plot_SeekMetabo_VertLine_Min,Plot_SeekMetabo_LoFiBasisRegion_VertLine_Max];
                Plot_SeekMetabo_LoFiBasisRegion_HiFi_VertLine_x = repmat(SeekMetabo_CS_zf(SeekMetabo_LoFiBasisRegion_SP(1)), [1 2]);
                Plot_SeekMetabo_LoFiBasisRegion_LoFi_VertLine_x = repmat(SeekMetabo_CS_zf(SeekMetabo_LoFiBasisRegion_SP(2)), [1 2]); 
                
                
                
                
                % Plotting
                SeekMetabo_fig = figure('visible','off'); 
                
                plot(Plot_SeekMetabo_Spec_x,Plot_SeekMetabo_Spec_y,'k','Linewidth',1.6);
                hold on

                % Baseline
                %plot(CS_vec_zf(Plot_Metaboregion_HiFi_SP:Plot_Metaboregion_LoFi_SP),polyval(SeekMetabo_Polyfit,Plot_Metaboregion_HiFi_SP:Plot_Metaboregion_LoFi_SP),'r','Linewidth',1.6)

                % Peak & Basis Points
                MetaboLegendHandi = zeros([1 4]);
                plot(SeekMetabo_CS_zf(FoundMetaboPeak_SP),FoundMetaboPeak_SS,'o','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','k');                               % Where peak was computed
                plot(SeekMetabo_CS_zf(FoundMetabo_Basis_SP(1)),FoundMetabo_Basis_SS(1),'o','MarkerSize',6,'MarkerFaceColor','c','MarkerEdgeColor','k');                    % Where HiFi_Basis was computed 
                plot(SeekMetabo_CS_zf(FoundMetabo_Basis_SP(2)),FoundMetabo_Basis_SS(2),'o','MarkerSize',6,'MarkerFaceColor','b','MarkerEdgeColor','k');                         % Where LoFi_Basis was computed
                MetaboLegendHandi(1) = plot(SeekMetabo_CS_zf(FoundMetaboPeak_SP),FoundMetabo_MeanBasis_SS,'o','MarkerSize',6,'MarkerFaceColor','g','MarkerEdgeColor','k');

                % Peak & Basis Regions
                % HiFi & LoFi SeekMetaboPeak_region
                MetaboLegendHandi(2) = plot(Plot_SeekMetabo_PeakRegion_HiFi_VertLine_x,Plot_SeekMetabo_PeakRegion_VertLine_y, '--r', 'Linewidth',1.6);                                                         
                plot(Plot_SeekMetabo_PeakRegion_LoFi_VertLine_x,Plot_SeekMetabo_PeakRegion_VertLine_y,'--r','Linewidth',1.6);
                % HiFi SeekMetaboBasis_HiFiRegion 
                MetaboLegendHandi(3) = plot(Plot_SeekMetabo_HiFiBasisRegion_HiFi_VertLine_x,Plot_SeekMetabo_HiFiBasisRegion_VertLine_y,'--c','Linewidth',1.6);                                          
                plot(Plot_SeekMetabo_HiFiBasisRegion_LoFi_VertLine_x,Plot_SeekMetabo_HiFiBasisRegion_VertLine_y,'--c','Linewidth',1.6); 
                % LoFi SeekMetaboBasis_LoFiRegion
                MetaboLegendHandi(4) = plot(Plot_SeekMetabo_LoFiBasisRegion_HiFi_VertLine_x,Plot_SeekMetabo_LoFiBasisRegion_VertLine_y,'--b','Linewidth',1.6);
                plot(Plot_SeekMetabo_LoFiBasisRegion_LoFi_VertLine_x,Plot_SeekMetabo_LoFiBasisRegion_VertLine_y,'--b','Linewidth',1.6); 


                hold off
                set(gca,'XDir','reverse');

                title(sprintf('Metabopeak %s, Metabopeak found = %d,voxel x_%d y_%d z_%d',SeekMetabo.InRealOrAbs,~NoMetabopeak_found, x,y,z),'Interpreter','none')
                xlabel('Chemical Shift')
                ylabel('Signal')
                leg = legend(MetaboLegendHandi, 'Mean Basis','Peak Region','HiFi Basis Region','LoFi Basis Region');
                set(leg,'FontSize',6);
                set(leg,'Location','Best');
                
                

                % Save & Close
                %saveas(SeekMetabo_fig,sprintf('%s/x_%s_y_%s_z_%s_SeekMetabo.fig', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))
                saveas(SeekMetabo_fig,sprintf('%s/x_%s_y_%s_z_%s_SeekMetabo', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)),'epsc2')  
                close(SeekMetabo_fig)
                
                
                
                %% 9.5 Whole Spectrum zf
%                 whole_spectrum_zf_fig = figure('visible','off');
%                 set(0,'CurrentFigure',whole_spectrum_zf_fig)                                                                   
%                 plot(CS_vec_zf,spectrum_zf_real)
%                 set(gca,'XDir','reverse');
%     		      title(sprintf('Spectrum With zerofilling, voxel x_%d y_%d z_%d', x,y,z),'Interpreter','none')
%                 xlabel('Chemical Shift')
%                 ylabel('Signal')         
%                 saveas(whole_spectrum_zf_fig,sprintf('%s/x_%s_y_%s_z_%s_whole_spectrum_zf.fig', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))
%                 saveas(whole_spectrum_zf_fig,sprintf('%s/x_%s_y_%s_z_%s_whole_spectrum_zf.eps', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))  




                %% 9.6 REFERENCE REGION
                
                
                % Definitions
                
                
                Plot_RefRegion_HiFi_SP = min([SeekRef_PeakRegion_SP(1), SeekRef_HiFiBasisRegion_SP(1)]);      
                Plot_RefRegion_LoFi_SP = max([SeekRef_PeakRegion_SP(2), SeekRef_LoFiBasisRegion_SP(2)]);
                

                % SPECTRUM ITSELF
                Plot_SeekRef_Spec_x = CS_vec_zf(Plot_RefRegion_HiFi_SP:Plot_RefRegion_LoFi_SP);
                Plot_SeekRef_Spec_y = eval([ 'spectrum_zf_' SeekRef_InRealOrAbs '(Plot_RefRegion_HiFi_SP:Plot_RefRegion_LoFi_SP)' ]);       % plot either real or absolute value depending on where Ref was sought


                % Vertical Lines for Peak & Basis Regions
                Plot_SeekRef_VertLine_Min = min(Plot_SeekRef_Spec_y);
                if(Plot_SeekRef_VertLine_Min >= 0)
                   Plot_SeekRef_VertLine_Min = -max(Plot_SeekRef_Spec_y)/20;                                                                % If the absolute value is plotted, then the vertical lines were so small
                end

                % PEAK REGION VERTICAL LINE
                Plot_SeekRef_PeakRegion_VertLine_Max = eval([ 'max(spectrum_zf_' SeekRef_InRealOrAbs '(SeekRef_PeakRegion_SP(1)),spectrum_zf_' SeekRef_InRealOrAbs '(SeekRef_PeakRegion_SP(2)));' ]);
                Plot_SeekRef_PeakRegion_VertLine_y = [Plot_SeekRef_VertLine_Min,Plot_SeekRef_PeakRegion_VertLine_Max];
                Plot_SeekRef_PeakRegion_VertLine_HiFi_x = repmat(CS_vec_zf(SeekRef_PeakRegion_SP(1)),[1 2]);                                % if you use just plot(x,vector), x scalar, only dots get plotted
                Plot_SeekRef_PeakRegion_VertLine_LoFi_x = repmat(CS_vec_zf(SeekRef_PeakRegion_SP(2)),[1 2]);                                % So plot([x,x],[y_start,y_end]) 

                % HiFi BASIS REGION VERTICAL LINE
                Plot_SeekRef_HiFiBasisRegion_VertLine_Max = ...
                eval([ 'max([spectrum_zf_' SeekRef_InRealOrAbs '(SeekRef_HiFiBasisRegion_SP(1)),spectrum_zf_' SeekRef_InRealOrAbs '(SeekRef_HiFiBasisRegion_SP(2))]);' ]);
                Plot_SeekRef_HiFiBasisRegion_VertLine_y = [Plot_SeekRef_VertLine_Min,Plot_SeekRef_HiFiBasisRegion_VertLine_Max];
                Plot_SeekRef_HiFiBasisRegion_HiFi_VertLine_x = repmat(CS_vec_zf(SeekRef_HiFiBasisRegion_SP(1)),[1 2]);
                Plot_SeekRef_HiFiBasisRegion_LoFi_VertLine_x = repmat(CS_vec_zf(SeekRef_HiFiBasisRegion_SP(2)), [1 2]);

                % LoFi BASIS REGION VERTICAL LINE
                Plot_SeekRef_LoFiBasisRegion_VertLine_Max = ...
                eval([ 'max([spectrum_zf_' SeekRef_InRealOrAbs '(SeekRef_LoFiBasisRegion_SP(1)),spectrum_zf_' SeekRef_InRealOrAbs '(SeekRef_LoFiBasisRegion_SP(2))]);' ]);
                Plot_SeekRef_LoFiBasisRegion_VertLine_y = [Plot_SeekRef_VertLine_Min,Plot_SeekRef_LoFiBasisRegion_VertLine_Max];
                Plot_SeekRef_LoFiBasisRegion_HiFi_VertLine_x = repmat(CS_vec_zf(SeekRef_LoFiBasisRegion_SP(1)), [1 2]);
                Plot_SeekRef_LoFiBasisRegion_LoFi_VertLine_x = repmat(CS_vec_zf(SeekRef_LoFiBasisRegion_SP(2)), [1 2]); 


                % Found Ref Peak Point
                Plot_SeekRef_PeakPoint_y = eval([ 'spectrum_zf_' SeekRef_InRealOrAbs '(FoundRefPeak_SP);' ]);
                
                
                
                
                % Plot
                SeekRef_fig = figure('visible','off'); 
                %set(0,'CurrentFigure',SeekRef_fig) 
                plot(Plot_SeekRef_Spec_x,Plot_SeekRef_Spec_y,'k','Linewidth',1.6)
                hold on

                % Peak & Basis Points            
                RefLegendHandi = zeros([1 4]);
                plot(CS_vec_zf(FoundRefPeak_SP),Plot_SeekRef_PeakPoint_y,'o','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','k')                           % Where peak was computed
                plot(CS_vec_zf(FoundRef_Basis_SP(1)),FoundRef_Basis_SS(1),'o','MarkerSize',6,'MarkerFaceColor','c','MarkerEdgeColor','k')                          % Where HiFi_Basis was computed 
                plot(CS_vec_zf(FoundRef_Basis_SP(2)),FoundRef_Basis_SS(2),'o','MarkerSize',6,'MarkerFaceColor','b','MarkerEdgeColor','k')                          % Where LoFi_Basis was computed
                RefLegendHandi(1) = plot(CS_vec_zf(FoundRefPeak_SP),FoundRef_MeanBasis_SS,'o','MarkerSize',6,'MarkerFaceColor','g','MarkerEdgeColor','k');

                % Peak & Basis Regions
                % HiFi & LoFi SeekRefPeak_region
                RefLegendHandi(2) = plot(Plot_SeekRef_PeakRegion_VertLine_HiFi_x,Plot_SeekRef_PeakRegion_VertLine_y, '--r', 'Linewidth',1.6);                       % LegendHandi means plot handle, Set handle for Legend
                plot(Plot_SeekRef_PeakRegion_VertLine_LoFi_x,Plot_SeekRef_PeakRegion_VertLine_y,'--r','Linewidth',1.6);
                % HiFi SeekRefBasis_HiFiRegion 
                RefLegendHandi(3) = plot(Plot_SeekRef_HiFiBasisRegion_HiFi_VertLine_x,Plot_SeekRef_HiFiBasisRegion_VertLine_y,'--c','Linewidth',1.6);               % Only the LegendHandis get displayed in Legend
                plot(Plot_SeekRef_HiFiBasisRegion_LoFi_VertLine_x,Plot_SeekRef_HiFiBasisRegion_VertLine_y,'--c','Linewidth',1.6); 
                % LoFi SeekRefBasis_LoFiRegion
                RefLegendHandi(4) = plot(Plot_SeekRef_LoFiBasisRegion_HiFi_VertLine_x,Plot_SeekRef_LoFiBasisRegion_VertLine_y,'--b','Linewidth',1.6);
                plot(Plot_SeekRef_LoFiBasisRegion_LoFi_VertLine_x,Plot_SeekRef_LoFiBasisRegion_VertLine_y,'--b','Linewidth',1.6); 

                hold off
                set(gca,'XDir','reverse');

                title(sprintf('RefPeak %s, RefPeak found = %d,voxel x_%d y_%d z_%d',SeekRef_InRealOrAbs,~NoRefPeak_found, x,y,z),'Interpreter','none')
                xlabel('Chemical Shift')
                ylabel('Signal')
                leg = legend(RefLegendHandi, 'Mean Basis Point','Peak Region','HiFi Basis Region','LoFi Basis Region');
                set(leg,'FontSize',6);
                set(leg,'Location','Best');


                % Save & Close
                %saveas(SeekRef_fig,sprintf('%s/x_%s_y_%s_z_%s_SeekRef.fig', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)))
                saveas(SeekRef_fig,sprintf('%s/x_%s_y_%s_z_%s_SeekRef', out_dir_FailedOrWon,num2str(x),num2str(y),num2str(z)),'epsc2') 
                close(SeekRef_fig)





            end % if(isfield(OutInfo,'out_dir'))
             
            
            

            %% 10. Loop END Preps
            
            CurrentVoxel = CurrentVoxel + 1;
            
            
     
        end % x-loop
	end % y-loop
end % z-loop





%% 11. Clip & Save SNR map and Shift Map, Close error log file


SNR_mat(SNR_mat == 0) = NaN;
SNR_mat_Reshape = reshape(SNR_mat, [numel(SNR_mat) 1]);            
SNR_mat_MedMad = nanmedian(SNR_mat_Reshape) + 13*mad(SNR_mat_Reshape,1);                  % clip values that are 13 median-absolute-deviations away from the median.
SNR_mat(SNR_mat > repmat(SNR_mat_MedMad,size(SNR_mat))) = NaN;                               % Clip Data.

Shift_mat(Shift_mat == 0) = NaN;
Shift_mat_Reshape = reshape(Shift_mat, [numel(Shift_mat) 1]);            
Shift_mat_MedMad = nanmedian(Shift_mat_Reshape) + 13*mad(Shift_mat_Reshape,1);                  % clip values that are 13 median-absolute-deviations away from the median.
Shift_mat(Shift_mat > repmat(Shift_mat_MedMad,size(Shift_mat))) = NaN;                               % Clip Data.


% Save Maps
if(isfield(OutInfo,'out_dir') && ~strcmpi(OutInfo.out_dir,'No'))
    %save(sprintf('%s/SNR_and_Shift.mat',OutInfo.out_dir), 'Shift_mat', 'SNR_mat')

    for z = 1:SLC
        SNRMap_fig = figure('visible','off');
        imagesc(squeeze(SNR_mat(:,:,z)))
        colorbar
        %saveas(SNRMap_fig,sprintf('%s/SNR_map.fig', OutInfo.out_dir))
        saveas(SNRMap_fig,sprintf('%s/SNR_map_slice%d', OutInfo.out_dir,z),'epsc2')
        close(SNRMap_fig)

        ShiftMap_fig = figure('visible','off');
        imagesc(squeeze(Shift_mat(:,:,z)))
        colorbar
        %saveas(ShiftMap_fig,sprintf('%s/Shift_map.fig', OutInfo.out_dir))
        saveas(ShiftMap_fig,sprintf('%s/Shift_map_slice%d', OutInfo.out_dir,z),'epsc2')
        close(ShiftMap_fig)
    end
    
    fclose(error_log_fid);
end



