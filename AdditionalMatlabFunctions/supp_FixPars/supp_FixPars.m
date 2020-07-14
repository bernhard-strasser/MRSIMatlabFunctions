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



%% Loop over Par and RecoPar

for CurParName2 = {'Par','RecoPar'}

    CurParName = CurParName2{1};
    if(~isfield(MRStruct,CurParName))
        continue;
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
    MRStruct.(CurParName).dimnames = MRStruct.(CurParName).dimnames(1:numel(size(MRStruct.Data)));

    
    
    %% Create dims
    if(~isfield(MRStruct.(CurParName),'dims'))
        for ii = 1:numel(MRStruct.(CurParName).dimnames)
            MRStruct.(CurParName).dims.(MRStruct.(CurParName).dimnames{ii}) = ii;
        end
    end    
    
     
end


%% Postparations

MRStruct = supp_UpdateRecoSteps(MRStruct,struct());

