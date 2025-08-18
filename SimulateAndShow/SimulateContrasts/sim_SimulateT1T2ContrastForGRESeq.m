function [Contrast, Signal] = sim_SimulateT1T2ContrastForGRESeq(TR_vec,TE_vec,T1_vec,T2Star_vec,FlipAngle_deg_vec)
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
% -         TR                     ...     TR of GRE sequence.
% -         T1                     ...     T1 of tissue.
% -         alpha                  ...     flip angle of sequence.
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

TR_vec = reshape(TR_vec,[1 numel(TR_vec)]);
FlipAngle_deg_vec = reshape(FlipAngle_deg_vec,[numel(FlipAngle_deg_vec) 1]);
TE_vec = reshape(TE_vec,[1 1 numel(TE_vec)]);



%%

Signal = ones([numel(FlipAngle_deg_vec) numel(TR_vec) numel(TE_vec) numel(T1_vec)]);
for ii = 1:numel(T1_vec)
    Signal(:,:,:,ii) = sim_SimulateSteadyStateSignal(1,TR_vec,T1_vec(ii),TE_vec, T2Star_vec(ii), FlipAngle_deg_vec);
end


%%

Contrast(:,:,:,1:numel(T1_vec)-1) = Signal(:,:,:,2:end) ./ Signal(:,:,:,1);


