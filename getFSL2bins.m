function fsl2bins = getFSL2bins(N,cent,delay)
%%
% The 2 bins method, after Köppl, 1997, J Neurophysiol

% input:
%   N: after a histogram was made from all spiketimes, N is the height of
%       the bins in the histogram
%   cent: cent are the centers of the bins in the histogram, a time vector
%   delay: the delay between start of recording and start of the click
% output:
%   fsl2bins: the first spike latency based on the 2 bins method.
% 
% By: Amarins Heeringa

%%
totDel = delay;
SFRbins = cent<totDel; SFRbins(1:10) = 0;
EFRbins = find(cent>totDel);
threshold = max(N(SFRbins));
if isempty(threshold)
    threshold = 1;
end
highbins = EFRbins(N(EFRbins)>threshold);
pnt = highbins(find(diff(highbins)==1,1));
fsl2bins=cent(pnt)-totDel;
if isempty(fsl2bins)
    fsl2bins=NaN;
end
