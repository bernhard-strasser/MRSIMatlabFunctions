function MRStruct = supp_FixPars(MRStruct)
%
% supp_FixPars Fix Parameters to for Hacking or Making Consistent
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function fixes the .Par and .RecoPar fields of an MRStruct to hard-code/hack parameters for special cases, 
% or to make Parameters consistent between different read-in-methods.
%
%
% [MRStruct] = supp_FixPars(MRStruct)
%
% Input: 
% -         MRStruct                     ...    MRStruct with fields Par or RecoPar   
%
% Output
% MRStruct:
% -         MRStruct                      ...     MRStruct with fixed fields Par or RecoPar 
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m


%% 0. Preparations

if(~isfield(MRStruct,'Par'))
    MRStruct.Par = struct;
end


%% Loop over Par and RecoPar

for CurParName2 = {'Par','RecoPar'}

    CurParName = CurParName2{1};
    if(~isfield(MRStruct,CurParName))
        continue;
    end

    % Some 'trivial' stuff
    if( ~isfield(MRStruct.(CurParName),'GyroMagnRatioOverTwoPi'))
        MRStruct.(CurParName).GyroMagnRatioOverTwoPi = 42.57747892 * 10^6;
    end
    if( ~isfield(MRStruct.(CurParName),'DataSize') && isfield(MRStruct,'Data'))
        MRStruct.(CurParName).DataSize = cellfun(@size,MRStruct.Data,'uni',false);     % Generalize later for more than 1 partition!
    end
    if( ~isfield(MRStruct.(CurParName),'nFreqEnc'))
        MRStruct.(CurParName).nFreqEnc = 64;
    end
    if( ~isfield(MRStruct.(CurParName),'nPhasEnc'))
        MRStruct.(CurParName).nPhasEnc = 64;
    end
    if( ~isfield(MRStruct.(CurParName),'nPartEnc'))
        MRStruct.(CurParName).nPartEnc = size(MRStruct.Data{1},3);
    end
    if( ~isfield(MRStruct.(CurParName),'nSLC'))
        MRStruct.(CurParName).nSLC = 1;
    end
    if( ~isfield(MRStruct.(CurParName),'vecSize'))
        MRStruct.(CurParName).vecSize = size(MRStruct.Data{1},5);
    end
    if( ~isfield(MRStruct.(CurParName),'FoV_Read'))
        MRStruct.(CurParName).FoV_Read = 220;
    end
    if( ~isfield(MRStruct.(CurParName),'FoV_Phase'))
        MRStruct.(CurParName).FoV_Phase = 220;
    end
    if( ~isfield(MRStruct.(CurParName),'FoV_Partition'))
        MRStruct.(CurParName).FoV_Partition = 12;
    end
    if( ~isfield(MRStruct.(CurParName),'Pos_Cor'))
        MRStruct.(CurParName).Pos_Cor = 0;
    end
    if( ~isfield(MRStruct.(CurParName),'Pos_Sag'))
        MRStruct.(CurParName).Pos_Sag = 0;
    end
    if( ~isfield(MRStruct.(CurParName),'Pos_Tra'))
        MRStruct.(CurParName).Pos_Tra = 0;
    end
    if( ~isfield(MRStruct.(CurParName),'SliceNormalVector_x'))
        MRStruct.(CurParName).SliceNormalVector_x = 0;
    end
    if( ~isfield(MRStruct.(CurParName),'SliceNormalVector_y'))
        MRStruct.(CurParName).SliceNormalVector_y = 0;
    end
    if( ~isfield(MRStruct.(CurParName),'SliceNormalVector_z'))
        MRStruct.(CurParName).SliceNormalVector_z = 1;
    end
    
    
    if( ~isfield(MRStruct.(CurParName),'total_channel_no_measured'))
        MRStruct.(CurParName).total_channel_no_measured = 1;
    end
    if( ~isfield(MRStruct.(CurParName),'ADC_dt'))
        MRStruct.(CurParName).ADC_dt = 5000;
    end
    if( ~isfield(MRStruct.(CurParName),'ADC_OverSamp'))
        MRStruct.(CurParName).ADC_OverSamp = 10000/MRStruct.(CurParName).ADC_dt;
    end
    if( ~isfield(MRStruct.(CurParName),'nAngInts'))
        MRStruct.(CurParName).nAngInts = size(MRStruct.Data,2);
    end
    
    
    %% Create nTempIntsPerAngInt
    if(~isfield(MRStruct.(CurParName),'nTempIntsPerAngInt') && isfield(MRStruct.(CurParName),'nTempInt'))
        MRStruct.(CurParName).nTempIntsPerAngInt = repmat(MRStruct.(CurParName).nTempInt,[1 MRStruct.(CurParName).nAngInts]);
    end

    
    %% ADC_dt
    if(~isfield(MRStruct.(CurParName),'ADCdtPerAngInt_ns'))
        if(isfield(MRStruct.(CurParName),'ADC_dt_ns'))
            if(numel(MRStruct.(CurParName).ADC_dt_ns) == MRStruct.(CurParName).nAngInts)
                MRStruct.(CurParName).ADCdtPerAngInt_ns = MRStruct.(CurParName).ADC_dt_ns;
            else
                 MRStruct.(CurParName).ADCdtPerAngInt_ns = repmat(MRStruct.(CurParName).ADC_dt_ns,[1 MRStruct.(CurParName).nAngInts]);
            end
            MRStruct.(CurParName) = rmfield(MRStruct.(CurParName),'ADC_dt_ns');
        end
        if(isfield(MRStruct.(CurParName),'ADC_dt'))
            if(numel(MRStruct.(CurParName).ADC_dt) == MRStruct.(CurParName).nAngInts)
                MRStruct.(CurParName).ADCdtPerAngInt_ns = MRStruct.(CurParName).ADC_dt;
            else
                 MRStruct.(CurParName).ADCdtPerAngInt_ns = repmat(MRStruct.(CurParName).ADC_dt,[1 MRStruct.(CurParName).nAngInts]);
            end
            MRStruct.(CurParName) = rmfield(MRStruct.(CurParName),'ADC_dt');
        end
    end

    %% Guess dimnames
    if(~isfield(MRStruct.(CurParName),'dimnames'))
        MRStruct.(CurParName).dimnames = {'angint','trajpt','kz','slc','t','cha','contrast'};           % A very rough guess. Make more sophisticated later?          
    end
    
    
    %% Cut dimnames
    if(isfield(MRStruct,'Data'))
        MRStruct.(CurParName).dimnames = MRStruct.(CurParName).dimnames(1:numel(size(MRStruct.Data)));
    end
    
    
    %% Create dims
    if(~isfield(MRStruct.(CurParName),'dims'))
        for ii = 1:numel(MRStruct.(CurParName).dimnames)
            MRStruct.(CurParName).dims.(MRStruct.(CurParName).dimnames{ii}) = ii;
        end
    end    
    
     
end


%% Postparations

MRStruct = supp_UpdateRecoSteps(MRStruct,struct());

