function Output = UntiltkSpaceTrajectory(Input,nTI,nrew,zff)
%
% Spice_SynthesizeMeasData 
%
% This function was written by Bernhard Strasser, April 2018.
%
%
% 
%
%
% Output = Spice_SynthesizeMeasData(U,V,B0Shiftmap,SamplingOperator)
%
% Input: 
% -         U                    ...     
% Output:
% -         Output                         ...     



% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!

% Further remarks: 




%% 0. Preparations



%%

% nTI = 1;            
ns = size(Input,1);
% nrew = Par.nRewPtsPerInt;
nc = size(Input,2);
vs = size(Input,3);

Output = Input;

% Normal Phaseroll
vs = zff*vs;
for t=1:ns

    kfaxis=Output(t,:,:);
    kfaxis = cat(3,kfaxis,zeros([size(kfaxis,1) size(kfaxis,2) (zff-1)*size(kfaxis,3)]));
    timeoffset=nTI*(t-1)/(nrew + ns);  % Do I really need to add the rewinder?
    phasecorr = reshape(exp(-2*1i*pi*timeoffset*((0:vs-1)/vs-0.5)),[1 1 vs]);


    kfaxis=fftshift(ifft(kfaxis,[],3),3);
    kfaxis=kfaxis.*repmat((phasecorr),[1 nc 1]);
    kfaxis=fft(fftshift(kfaxis,3),[],3);
    Output(t,:,:)=kfaxis(:,:,1:vs/zff);

end
Output(:,:,1) = Input(:,:,1);       % I dont know if this is helpful...



