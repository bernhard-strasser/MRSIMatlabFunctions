function [SNR,Signal,Noise] = CalcSNRFromLCModel(CoordFolderPath,SNRCalcPar,MatSize)
%
% CalcSNRFromLCModel Calculate SNR Based on LCModel Coord-Files
%
% This function was written by Bernhard Strasser, April 2018.
%
%
% The function reads "Coordinate" LCModel files, and calculates the SNR using a user-specified range . 
%
%
% [SNR,Signal,Noise] = CalcSNRFromLCModel(CoordFolderPath,SNRCalcPar,MatSize)
%
% Input: 
% -         CoordFolderPath             ...   Path of folder containing all the coord-files. 
% -         SNRCalcPar                  ...   Cell of structure containing the parameters for calculating the SNR. Each cell must have fields:
%                      .Name            ...   Name of Peak for which the SNR should be calculated.
%                      .PeakPPMRange    ...   Range in ppm where the peak should be searched.
%                      .NoisePPMRange   ...   Range in ppm where the std for the noise-std should be calculated.

% -         MatSize                     ...   Matrix size of the output. Not all coord-files of all voxels might be available, so the
%                                             resulting SNR-map might have not the expected size.
%
% Output:
% -         SNR                         ...   Matrix containing the calculated SNR
% -         Signal                      ...   Matrix containing the calculated Signal
% -         Noise                       ...   Matrix containing the calculated Noise
%
%
%
% REMARK: IF YOUR COORD-FILES DO NOT HAVE THE VOXEL-POSITION ENCODED AS "x[0-9][0-9]" "y[0-9][0-9]" AND "z[0-9][0-9]",
% YOU WILL GET A VECTOR AS RESULTING SNR
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None






%% 0. Preparations & Housekeeping




%% 1. Get a List of All Coord-Files

CoordFiles = dir(CoordFolderPath);

CoordFiles = {CoordFiles.name};
CoordFiles(1:2) = []; 
testEmpty = regexpi(CoordFiles,'\.coord');
testEmpty = cellfun('isempty',testEmpty);
CoordFiles(testEmpty) = [];



%% Try to Find x,y,z Positions of CoordFiles

XInd = ones([1 numel(CoordFiles)]); YInd = ones([1 numel(CoordFiles)]); ZInd = ones([1 numel(CoordFiles)]);
XYZIndFound = true;
for CurVox = 1:numel(CoordFiles)
    
    Dummy_X = regexpi(CoordFiles{CurVox},'x[0-9][0-9]');
    Dummy_Y = regexpi(CoordFiles{CurVox},'y[0-9][0-9]');
    Dummy_Z = regexpi(CoordFiles{CurVox},'z[0-9][0-9]');

    if( ~isempty(Dummy_X) && ~isempty(Dummy_Y) && ~isempty(Dummy_Z) )
        XInd(CurVox) = str2double(CoordFiles{CurVox}(Dummy_X+1:Dummy_X+2));
        YInd(CurVox) = str2double(CoordFiles{CurVox}(Dummy_Y+1:Dummy_Y+2));
        ZInd(CurVox) = str2double(CoordFiles{CurVox}(Dummy_Z+1:Dummy_Z+2));
    else
        XYZIndFound = false;
        break;
    end
end
if(~XYZIndFound)
    XInd = 1:numel(CoordFiles);
    YInd = ones([1 numel(CoordFiles)]); ZInd = ones([1 numel(CoordFiles)]);
end

if(~exist('MatSize','var'))
    MatSize = [max(XInd) max(YInd) max(ZInd)];
end


%% Initialize SNR

for CurPeak = 1:numel(SNRCalcPar)
    SNR.(SNRCalcPar{CurPeak}.Name) = NaN(MatSize);
    if(nargout > 1)
        Signal.(SNRCalcPar{CurPeak}.Name) = NaN(MatSize);
    end
    if(nargout > 2)
        Noise.(SNRCalcPar{CurPeak}.Name) = NaN(MatSize);
    end
end


%% Read in Coord-files and Calculate SNR for each

for CurVox = 1:numel(CoordFiles)
    
    CurSpec = read_coord_file([CoordFolderPath '/' CoordFiles{CurVox}]);
    if(numel(CurSpec.ppm_points) == 1 && CurSpec.ppm_points == 0) % In case an error occurs when reading the coord-file
        continue
    end
    SignalData = CurSpec.Spec - CurSpec.Baseline;

    
    % Loop over all peaks for which the SNR should be calculated
    for CurPeak = 1:numel(SNRCalcPar)

        Pointys = FindClosestIndex(CurSpec.ppm_points,[SNRCalcPar{CurPeak}.NoisePPMRange SNRCalcPar{CurPeak}.PeakPPMRange]);
        NoisePointRange = [Pointys{[2 1]}];
        PeakPointRange = [Pointys{[4 3]}];

        % Calc Signal
        CurSignal = max(SignalData(PeakPointRange(1):PeakPointRange(2)));

        % Calc Noise
        CurNoise = std(CurSpec.Spec(NoisePointRange(1):NoisePointRange(2)));
        
        SNR.(SNRCalcPar{CurPeak}.Name)(XInd(CurVox),YInd(CurVox),ZInd(CurVox)) = CurSignal/(2*CurNoise);

        if(nargout > 1)
            Signal.(SNRCalcPar{CurPeak}.Name)(XInd(CurVox),YInd(CurVox),ZInd(CurVox)) = CurSignal;
        end
        if(nargout > 2)
            Noise.(SNRCalcPar{CurPeak}.Name)(XInd(CurVox),YInd(CurVox),ZInd(CurVox)) = CurNoise;
        end

    end
    
end





