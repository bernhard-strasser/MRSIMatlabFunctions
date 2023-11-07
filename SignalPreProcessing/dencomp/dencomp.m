

function [deltak,Kprime,A] = dencomp(n,res,FOV,deltak1,x)
%dencomp calculates the increments deltak for the cirlces used in CONCEPT
%such that the density is hamming distributed.
%   n ... numbers of samples per circle (same for each circle)
%   res....resolution/matrixsize in 1D. (e.g. 64) 
%   FOV...field of view in meter e.g. FOV=0.22m
%   --->kmax=res/FOV ... maximum k-space radius
%   N ... number of circles
%   x... gives the factor x=1,2,3,... which reduces the number of circles
%   such that N->N/x
%   deltak1....first increment=radius of the inner circle=start of
%   iteration.
%   deltakmax=1/FOV.....maximal increment

N=1;

deltak=zeros(1,N);
deltak(1)=deltak1;
count=2;
N=ceil((N)/x);
A=zeros(1,N+1);

kmax=res/FOV/2;

while (sum(deltak)<kmax)

    deltakS=sum(deltak);
    
    Hammy = (0.54-0.46*cos(2*pi*(deltakS-kmax)/2/kmax));
    
    deltak(count)=2*sqrt(n/Hammy/pi+(deltakS-deltak(count-1)/2)^2)-2*deltakS;    
    %deltak(count)=deltak(count);
    
    if (deltak(count)<0 && deltak(count-1)>0)                   % is this alright? from here...
                                                                 
        deltak(count-1)=deltak(count-1)+deltak(count);
        deltak(count)=-deltak(count);
       
    end;                                                        %...to here.
        
    
    if deltak(count)>1/FOV;
        deltak(count)=1/FOV;
    end;
    
    
    
    count=count+1;
    
end;


%reduce number of circles by a factor of 1/x,
% deleting circles is meaning that the increments are adding up
for i=2:x:length(deltak)-x+1
   deltak(i)=sum(deltak(i:i+x-1));
end;

deltak=[deltak(1) deltak(2:x:length(deltak)-x+1)];
deltak=[0 deltak]; % because i want that the first circle always stays there (it shall never be deleted)


% for i=1:length(deltak)-1                       % this is not optimal
%    if deltak(i)+deltak(i+1)<1/FOV
%      deltak(i)=sum(deltak(i:i+1));
%      deltak(i+1)=[];
%      deltak=[deltak 0];
%    end;
% end;

%A=[0 0 A];
%A=A(2:x:length(A));






deltak=[deltak deltak(end)];
%Radius=cumsum(deltak);


%correction for the first three points which are wrong
A=zeros(1,length(deltak)+1);
%deltak(2)=2/sqrt(pi); % irgendwie so!!!!!
%A(1)=pi*deltak(2)^2/4;
A(1)=pi*(deltak(2)+deltak(3))^2;
%correction ends here



for numb=2:(length(deltak)-1)

    %deltak(N+2)=0;
    A(numb)=((sum(deltak(2:numb))+deltak(numb+1)/2)^2-(sum(deltak(2:numb))-deltak(numb)/2)^2)*pi/n;   


end;

A(2)=pi*(deltak(2)+deltak(3))^2;
%A(3)=pi*(deltak(2)+deltak(3))^2;

deltak(1)=[];
%deltak(end)=[];
A(end)=[];
A(end)=[];

% FTs
P=zeros(1,length(deltak)+1);
C=zeros(1,length(deltak)+1);
P(1)=[];
syms r P p z C c

%the polar-psf was derived analytically by hankel-transforming the sampling
%density. p is devided by a constant factor to improve visualization.

for m=2:(length(deltak))
P(m)=besselj(0,2*pi*r*sum(deltak(2:m)))*sum(deltak(2:m));
end;

p=sum(P(2:end));

%cpsf = cartesian psf

for m=2:(length(deltak))
C(m)=2*cos(2*pi*z*sum(deltak(2:m)));
end;

c=1+sum(C(2:end));

%PLOTs

% polar PSF= ppsf 

subplot(3,3,1)
%figure
ezplot(p/1000);

% cpsf


subplot(3,3,2)
%figure
ezplot(c/100);




%line plot of density and hamming filter



subplot(3,3,3)
Y=(0.54-0.46*cos(2*pi*(cumsum(deltak)-max(cumsum(deltak)))/2/max(cumsum(deltak))));
%figure
plot(cumsum(deltak),Y,cumsum(deltak),(1./A)*x); 


kprime=1:length(deltak);
Kprime=(5*kmax*sqrt(pi))./(sqrt(length(deltak))*sqrt(abs(length(deltak)*atan(1/5*sqrt(2)*tan(pi*kprime/length(deltak)/2)))).*(27+23*cos(2*pi*kprime/length(deltak)/2)));


subplot(3,3,7)
plot(1:length(deltak),deltak,kprime,Kprime);


k=1:length(deltak);
K=kmax*sqrt(2/pi)/sqrt(length(deltak))*sqrt(abs(length(deltak)*atan(1/5*sqrt(2)*tan(pi*k/length(deltak)/2))));


subplot(3,3,8)
plot(1:length(deltak),cumsum(deltak),k,K);

subplot(3,3,9)
plot(1:length(deltak),1./(K.*Kprime));

subplot(3,3,4)
scatter(cumsum(deltak),Y)
subplot(3,3,5)
scatter(cumsum(deltak),1./A*x,'g')




%k-space trajectory
theta=0:0.1:2.01*pi;


subplot(3,3,6)
%figure;
for counti=drange(1:length(deltak))
hold all
polar(theta,sum(deltak(1:counti)).*sqrt(sin(theta).^2+cos(theta).^2),'-r');
end;




end

