function [Output, AdditionalOut] = op_ChangeCellDependencyFor3DCRT(Output,AdditionalIn,Settings)
%
% read_csi_dat Read in raw data from Siemens
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function can read in MRS(I) data in the Siemens raw file format ".DAT" and performs
% some easy Postprocessing steps like zerofilling, Hadamard decoding, Noise Decorrelation etc.
%
%
% [Output, AdditionalOut] = op_ReconstructNonCartMRData(Output,AdditionalIn,Settings)
%
% Input: 
% -         ?                     ...     
% -         ?                     ...     
% -         ?             ...     
%
% Output:
% -         ?                      ...     
% -         ?                        ...     
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: This function uses FFTs to get from k- to image-space. This is mathematically wrong, but Siemens seems to do the same when
% creating DICOMS. The only difference is that the images are flipped left/right and up/down.

% This function expects the input to be of form
% [nCha, nAngInt 


%% 0. Preparations


if(~exist('Settings','var') || ~isfield(Settings,'ReshapeToMatrix_flag'))
    Settings.ReshapeToMatrix_flag = true;
end

if(~isfield(Output,'RecoPar'))
    if(~isfield(Output,'Par'))
        error('Output must have field Par or RecoPar.')
    end
    Output.RecoPar = Output.Par;
end
Output.RecoPar.DataSize = [size_MultiDims(Output.OutTraj.GM,[3 4]) Output.RecoPar.nPartEnc_Meas Output.RecoPar.nSLC ...
                           Output.RecoPar.vecSize Output.RecoPar.total_channel_no_measured Output.RecoPar.nRep];



%% Change cell-dependency of nAngInt --> nPart (each Partition one cell-element)
% from {nAngInt}(nTrajPoints x 1 x nPart x nSlc x nTempInt*vecSize x nCha) --> {nPart}(nAngInt*nTrajPoints x Rest)

Temp2 = cell([1 Output.RecoPar.nPartEnc_Meas]);
Temp3 = cell([1 Output.RecoPar.nPartEnc_Meas]);
for CurPartEnc = 1:Output.RecoPar.nPartEnc_Meas
    Sz = cat(1,Output.Par.DataSize{Output.RecoPar.AngIntsPerPartEnc(CurPartEnc,:)}); 
    Temp2{CurPartEnc} = zeros([sum(Sz(:,1)) Sz(1,2) 1 Sz(1,4:end)],'single');
    CurPt = 1;
    for ii = 1:numel(Output.Data)
        if(Output.RecoPar.AngIntsPerPartEnc(CurPartEnc,ii))
            Temp1 = Output.Data{ii}(:,:, CurPartEnc-sum(~Output.RecoPar.AngIntsPerPartEnc(1:CurPartEnc-1,ii)),:,:,:);
            Temp2{CurPartEnc}(CurPt:CurPt+size(Temp1,1)-1,:,:,:,:,:) = Temp1;
            CurPt = CurPt + size(Temp1,1);
        end
    end
    SizeData_k{CurPartEnc} = size(Temp2{CurPartEnc}); %SizeData_k = cat(2,SizeData_k,ones([1 5-numel(SizeData_k)]));
    if(Settings.ReshapeToMatrix_flag)
        Temp2{CurPartEnc} = Temp2{CurPartEnc}(:,:);
    end
    
    if(isfield(Output,'NoiseData') && numel(Output.NoiseData) > 1)
        Temp3{CurPartEnc} = zeros([sum(Sz(:,1)) Sz(1,2) 1 Sz(1,4:end)],'single');
        CurPt = 1;
        for ii = 1:numel(Output.NoiseData)
            if(Output.RecoPar.AngIntsPerPartEnc(CurPartEnc,ii))
                Temp1 = Output.NoiseData{ii}(:,:, CurPartEnc-(Output.RecoPar.nPartEnc_Meas-sum(Output.RecoPar.AngIntsPerPartEnc(:,ii)))/2,:,:,:);
                Temp3{CurPartEnc}(CurPt:CurPt+size(Temp1,1)-1,:,:,:,:,:) = Temp1;
                CurPt = CurPt + size(Temp1,1);
            end
        end
        if(Settings.ReshapeToMatrix_flag)
            Temp3{CurPartEnc} = Temp3{CurPartEnc}(:,:);
        end
    end
end

Temp4{1} = cat(2,Output.InTraj.GM{:}); Output.InTraj.GM = Temp4;

if(isfield(Output,'NoiseData') && numel(Output.NoiseData) > 1)
    Output.NoiseData = Temp3;
end
Output.Data = Temp2; clear Temp1 Temp2 Temp3 Temp4


