function [SteadyStateSignal] = sim_SimulateSteadyStateSignal(S0,TR,T1,TE,T2Star,FlipAngle_deg)
%
% SteadyStateSignal Calculate the signal in steady state according to the formula
%
% This function was written by Bernhard Strasser, July 2025.
%
% 
%
%
% [SteadyStateSignal] = CalculateSteadyStateSignal(S0,TR,T1,alpha)
%
% Input: 
% -         S0                     ...     Signal in beginning.
% -         TR                     ...     TR of GRE sequence. Scalar or row vector (lying vector)
% -         T1                     ...     T1 of tissue.
% -         alpha                  ...     flip angle of sequence. Scalar or column vector (standing vector)
% or
%
% Output:
% -         SteadyStateSignal      ...     Signal in steady state
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None


%%

FlipAngle_deg = FlipAngle_deg/180*pi;
SteadyStateSignal = S0*sin(FlipAngle_deg).*(1-exp(-TR/T1))./(1-cos(FlipAngle_deg).*exp(-TR/T1)) .* exp(-TE/T2Star);



