function [fsl,pd_fsl] = getFSLpoisson(sptimes,delay)
%%
% The poisson surprise method, after Chase & Young, 2007, PNAS

% input:
%   sptimes: a vector with all of the spike times in the recording
%   delay: the delay between start of recording and start of the click
% output:
%   fsl: the first spike latency based on the poisson surprise method.
%   pd_fsl: the poisson probability density function, can be plotted as a 
%       function of the recording time. 
% 
% By: Amarins Heeringa, based on code by Calvin Wu.

%%

sptimesN=sort(sptimes); % all sorted
spac=(length(sptimesN(sptimesN<delay)))/delay;
pd_fsl=[];
for m=0:length(sptimesN)-1
    pd_fsl(m+1)=poisspdf(m,spac*sptimesN(m+1));
end

if ~isempty(find(pd_fsl<10^-6&sptimesN>delay))
    fx=find(pd_fsl<10^-6&sptimesN>delay);
    fsl=sptimesN(fx(1));
else
    fsl=NaN;
end
