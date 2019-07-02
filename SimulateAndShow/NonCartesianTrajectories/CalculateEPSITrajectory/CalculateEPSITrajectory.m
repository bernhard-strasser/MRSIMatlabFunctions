function [GradWavForm,SpectrRange_Used,GradInfo] = CalculateEPSITrajectory(Par)
%
% CalculateEPSITrajectory Calculate gradients for standard EPSI trajectory
%
% This function was written by Bernhard Strasser, July 2017.
%
% The function takes some input Par like the maximum slewrate to use and the duration of the ramp time etc. and gives back the gradient waveform along
% with some additional information about at which point which gradient waveform element starts and ends (ramp up, ramp down, flat top positive, flat top negative).
%
%
% [GradWavForm,GradInfo] = CalculateEPSITrajectory(Par.Sys.SlewRate_Max,Par.MeasProt.vecSizePerTI,TRamp,TFlat,Par.Sys.GRAD_RASTER_TIME)
%
% Input: 
% -         Par                               ...    struct with parameteris with fields
% --            .MeasProt                     ...    struct defining the measurement protocol with fields
% ---                    .RampSampling_Flag   ...    Flag which determines if the samples during the ramps will be used. If they will, we don't need to go out that far in
%                                                    k-space with the flat-top gradients, because during de-acceleration we still go a little further in the same direction.
%                                                    Default: true.
% ---                    .SpectrRange_Target  ...    The desired spectral range/bandwidth. The used bandwidth might slightly differ because the gradient points need to lie
%                                                    on the gradient grid defined by Par.Sys.GRAD_RASTER_TIME. Default: 2700 Hz.
% ---                    .nTI                 ...    The number of temporal interleaves. So far, the gradient waveforms are not calculated for the different TIs!!!!!!
%                                                    Default: 1.
% ---                    .Nx                  ...    The matrix size in x-direction. Default: 8.
% ---                    .FOVx                ...    The FOVx in x-direction. Default: 0.22 m
% ---                    .vecSizePerTI        ...    The simulated vector size, i.e. how often the trajectory should go back and forth. Default: 4.
% --            .Sys                          ...    struct defining the system properties with fields
% ---                    .Gamma               ...    The gyromagnetic ratio. Default 42.5756*10^6 Hz/T.
% ---                    .SlewRate_Max         ...    The maximum allowed slew rate. Default: 200 mT/m/ms = T/m/s.
% ---                    .GRAD_RASTER_TIME    ...    The raster time in s, i.e. the time between two calculated gradient time points. Default: 1E-6 s
%
% Output:
% -         GradWavForm                       ...    The gradient waveform in T/m.
% -         SpectrRange_Used                  ...    The given SBW cannot be achieved, the gradient points need to lie on the gradient raster. This is the updated SBW.
% -         GradInfo                          ...    Additional information about the gradient waveform, i.e. when e.g. the ramp up, ramp down, flat top etc. start and end
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: ???

% Further remarks: 

%% Definitions


if(exist('Par','var') && ~isstruct(Par))
    fprintf('\n\nWarning: The input variable ''Par'' needs to be a structure, see help CalculateEPSITrajectory. Assign standard values.\n')
    clear Par;
end
       
    
if(~exist('Par','var') || ~isfield(Par,'Sys') || ~isstruct(Par.Sys))
    Par = rmfield(Par,'Sys');
    Par.Sys.Gamma = 42.5756*10^6;
    Par.Sys.SlewRate_Max = 200;
    Par.Sys.GRAD_RASTER_TIME = 1E-6;
end
if(~isfield(Par.Sys,'Gamma') || Par.Sys.Gamma == 0)
    Par.Sys.Gamma = 42.5756*10^6;
end
if(~isfield(Par.Sys,'SlewRate_Max') || Par.Sys.SlewRate_Max == 0)
    Par.Sys.SlewRate_Max = 200;
end
if(~isfield(Par.Sys,'GRAD_RASTER_TIME') || Par.Sys.GRAD_RASTER_TIME == 0)
    Par.Sys.GRAD_RASTER_TIME = 1E-6;
end

if(~isfield(Par,'MeasProt') || ~isstruct(Par.MeasProt))
    Par = rmfield(Par,'MeasProt');
    Par.MeasProt.RampSampling_Flag = true;    
    Par.MeasProt.SpectrRange_Target = 2700;
    Par.MeasProt.nTI = 1;
    Par.MeasProt.Nx = 8;
    Par.MeasProt.FOVx = 0.22;
    Par.MeasProt.vecSizePerTI = 4;
end
if(~isfield(Par.MeasProt,'RampSampling_Flag'))
    Par.MeasProt.RampSampling_Flag = true;
end
if(~isfield(Par.MeasProt,'SpectrRange_Target') || Par.MeasProt.SpectrRange_Target == 0)
    Par.MeasProt.SpectrRange_Target = 2700;
end
if(~isfield(Par.MeasProt,'nTI') || Par.MeasProt.nTI == 0)
    Par.MeasProt.nTI = 1;
end
if(~isfield(Par.MeasProt,'Nx') || Par.MeasProt.Nx == 0)
    Par.MeasProt.Nx = 8;
end
if(~isfield(Par.MeasProt,'FOVx') || Par.MeasProt.FOVx == 0)
    Par.MeasProt.FOVx = 0.22;
end
if(~isfield(Par.MeasProt,'vecSizePerTI') || Par.MeasProt.vecSizePerTI == 0)
    Par.MeasProt.vecSizePerTI = 4;
end






%% Calculate the Timings and Effective SlewRate


% w/o ramp sampling
if(~Par.MeasProt.RampSampling_Flag)
    TRamp = Par.MeasProt.nTI./(8*Par.MeasProt.SpectrRange_Target) - sqrt((Par.MeasProt.nTI./Par.MeasProt.SpectrRange_Target).^2 / 64 - ...
            (Par.MeasProt.Nx-1)./(2*Par.Sys.Gamma*Par.MeasProt.FOVx*Par.Sys.SlewRate_Max));
    TFlat = Par.MeasProt.nTI ./ Par.MeasProt.SpectrRange_Target/2 - 2*TRamp;    
else
    % For Ramp sampling, the formula change slightly, because the k-space extent we need to cover is now the area of the flat-top + 2*Ramp
    % instead of only the area of the flat-top in the normal sampling
    TFlat = 2*sqrt((Par.MeasProt.nTI./(4*Par.MeasProt.SpectrRange_Target)).^2 - (Par.MeasProt.Nx-1)./(Par.Sys.Gamma*Par.MeasProt.FOVx*Par.Sys.SlewRate_Max));
    TRamp = Par.MeasProt.nTI./(4*Par.MeasProt.SpectrRange_Target) - TFlat/2;
end


if(abs(imag(TRamp)) > 0)
    GradWavForm = NaN; SpectrRange_Used = NaN; GradInfo = NaN;
    return;
end


% Round TRamp and TFlat to Par.Sys.GRAD_RASTER_TIME
TRamp = ceil(TRamp/Par.Sys.GRAD_RASTER_TIME)*Par.Sys.GRAD_RASTER_TIME;
TFlat = ceil(TFlat/Par.Sys.GRAD_RASTER_TIME)*Par.Sys.GRAD_RASTER_TIME;

% Update real spectral bandwidth (the given SBW cannot be achieved, the gradient points need to lie on the gradient raster, and therefore we need a slightly different SBW)
% This parameter is not used, it's only an output-parameter
SpectrRange_Used = Par.MeasProt.nTI./(4*TRamp+2*TFlat);

% Since we needed to round up the TRamp, and round TFlat, we need to change the used slewrate, so that the area remains the same.
if(Par.MeasProt.RampSampling_Flag)
    SlewRate_Used = (Par.MeasProt.Nx-1)./(Par.Sys.Gamma*Par.MeasProt.FOVx) ./ (TRamp .* (TRamp + TFlat));
else
    SlewRate_Used = (Par.MeasProt.Nx-1)./(Par.Sys.Gamma*Par.MeasProt.FOVx) ./ (TRamp .* TFlat);    
end


%% Calculate the Trajectory


GradWavForm = zeros([1 round(Par.MeasProt.vecSizePerTI*(4*TRamp+2*TFlat)/Par.Sys.GRAD_RASTER_TIME)]);
GradInfo{1,1} = 'GradType'; GradInfo{1,2} = 'BeginPoint(incl)'; GradInfo{1,3} = 'EndPoint(incl)';



NoOfRampPts = round(TRamp/Par.Sys.GRAD_RASTER_TIME);
NoOfAcqPts = round(TFlat/Par.Sys.GRAD_RASTER_TIME);

TimePtsRamp = (1:NoOfRampPts);
CurTime = 1;    

% Ramp Up
GradWavForm((CurTime+1):(CurTime + NoOfRampPts)) = SlewRate_Used * TimePtsRamp;
CurTime = CurTime + NoOfRampPts;
GradInfo{2,1} = 'RampUpTime'; GradInfo{2,2} = CurTime-NoOfRampPts; GradInfo{2,3} = CurTime-1;


% Flat Top
GradWavForm((CurTime+1):(CurTime + (NoOfAcqPts))) = GradWavForm(CurTime);
CurTime = CurTime + NoOfAcqPts;
GradInfo{3,1} = 'AcquisitionTime'; GradInfo{3,2} = CurTime-NoOfAcqPts; GradInfo{3,3} = CurTime;


% Ramp Down
GradWavForm((CurTime+1):(CurTime + NoOfRampPts)) = GradWavForm(CurTime) - SlewRate_Used * TimePtsRamp;
CurTime = CurTime + NoOfRampPts - 1;        % -1: Without that we would append the first 0 of the trajectory beginning with the 0 at the end. But we don't want two 0's


% Take the whole first positive part, and negate and append it
GradWavForm((CurTime+1):(CurTime + 2*NoOfRampPts+NoOfAcqPts)) = -GradWavForm(1:CurTime);
CurTime = CurTime + 2*NoOfRampPts+NoOfAcqPts;
GradInfo{4,1} = 'RampDownTime'; GradInfo{4,2} = CurTime-2*NoOfRampPts-NoOfAcqPts-NoOfRampPts+1+1; GradInfo{4,3} = CurTime-NoOfAcqPts-NoOfRampPts;
GradInfo{5,1} = 'AcquisitionTime'; GradInfo{5,2} = CurTime-NoOfAcqPts-NoOfRampPts+1; GradInfo{5,3} = CurTime-NoOfRampPts+1;
GradInfo{6,1} = 'RampUpTime'; GradInfo{6,2} = CurTime-NoOfRampPts+1+1; GradInfo{6,3} = CurTime;


% Now Replicate this whole thing vecSize Times -1 (-1: One we have already...)
GradWavForm((CurTime+1):end) = repmat(GradWavForm(1:CurTime), [1 Par.MeasProt.vecSizePerTI-1]);

for ii = 1:Par.MeasProt.vecSizePerTI-1
    GradInfo(7+(ii-1)*5:6+ii*5,1:3) = GradInfo(2:6,1:3);
    for jj=1:5
        GradInfo{6+(ii-1)*5+jj,2} = GradInfo{6+(ii-1)*5+jj,2} + ii*GradInfo{6,3};
        GradInfo{6+(ii-1)*5+jj,3} = GradInfo{6+(ii-1)*5+jj,3} + ii*GradInfo{6,3};
    end
end


% Rescale GradWavForm to take into account that each point was Par.Sys.GRAD_RASTER_TIME long, unit: T/m
GradWavForm = GradWavForm * Par.Sys.GRAD_RASTER_TIME;


% Now add the time information to it
GradWavForm = cat(1,(1:numel(GradWavForm))*Par.Sys.GRAD_RASTER_TIME,GradWavForm);  








