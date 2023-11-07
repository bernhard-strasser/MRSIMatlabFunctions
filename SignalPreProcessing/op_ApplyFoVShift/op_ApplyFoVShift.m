function [MRStruct, AdditionalOut] = op_ApplyFoVShift(MRStruct,AdditionalIn,Settings)
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
% [MRStruct, AdditionalOut] = op_ReconstructNonCartMRData(MRStruct,AdditionalIn,Settings)
%
% Input: 
% -         ?                     ...     
% -         ?                     ...     
% -         ?             ...     
%
% MRStruct:
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


if(~exist('AdditionalIn','var'))
    AdditionalIn = struct();
end
if(~exist('Settings','var'))
    Settings = struct();
end
if(~isfield(Settings,'UndoFoVShift_flag'))
    Settings.UndoFoVShift_flag = false;    
end
if(~isfield(MRStruct,'RecoPar'))
    if(~isfield(MRStruct,'Par'))
        error('MRStruct must have field Par or RecoPar.')
    end
    MRStruct.RecoPar = MRStruct.Par;
end
TmpDataSize = [size_MultiDims(MRStruct.OutTraj.GM,[3 4]) MRStruct.RecoPar.nPartEnc MRStruct.RecoPar.nSLC ...
                           MRStruct.RecoPar.vecSize MRStruct.RecoPar.total_channel_no_measured];


                      
% MRStruct = supp_FixPars(MRStruct);  % To hard-code/hack parameters for special cases, or to make Parameters consistent between different read-in-methods.



%% FOV SHIFTs Add correct phaseses to the data and shift the FOV to the image center.

LPH=[MRStruct.RecoPar.Pos_Cor MRStruct.RecoPar.Pos_Sag MRStruct.RecoPar.Pos_Tra];
Normal1=[MRStruct.RecoPar.SliceNormalVector_y MRStruct.RecoPar.SliceNormalVector_x MRStruct.RecoPar.SliceNormalVector_z];
Normal2=[0 0 1];
v=vrrotvec(Normal1,Normal2);
Rot=vrrotvec2mat(v);

% In Plane Rotation?
% Normal2=[1 0 0];
% Normal1 = [cos(InPlaneRot) sin(InPlaneRot) 0];
% v=vrrotvec(Normal1,Normal2);
% Rot2=vrrotvec2mat(v);
% Rot=Rot2*Rot;


AllFields = fieldnamesr(MRStruct.RecoSteps);
ConjFields = AllFields(~cellfun(@isempty,regexp(AllFields,'ConjInkSpace_flag|ConjIniSpace_flag|ConjFlag')));
ConjNumber = 0; for ii=1:numel(ConjFields); ConjNumber = ConjNumber + eval(['MRStruct.RecoSteps.' ConjFields{ii}]);end
ConjFlag = mod(ConjNumber,2);

if(ConjFlag)
   ConjSign = -1;
else
  ConjSign = 1; 
end
if(Settings.UndoFoVShift_flag)
   Sign = -1; 
else
    Sign = 1;
end
Sign = ConjSign * Sign;

PRS=Rot*LPH';
FOVShift = cellfun( @(x) transpose(exp(Sign*1i*x(1,:)/0.5*TmpDataSize(2)*pi*-PRS(2)/MRStruct.RecoPar.FoV_Read)), MRStruct.InTraj.GM , 'uni', false);

FOVShift2 = cellfun( @(x) transpose(exp(Sign*1i*x(2,:)/0.5*TmpDataSize(1)*pi*PRS(1)/MRStruct.RecoPar.FoV_Phase)),MRStruct.InTraj.GM,'uni',false);

MRStruct.Data = cellfun( @(x,y,z) x.*(y.*z),MRStruct.Data, FOVShift,FOVShift2,'uni',false);



%% Postparations

MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);



