function [tlen,StartTime,EndTime,nwSpikCnt,nwSpikRat]=getSpikeCounts(curvedata,curvesettings)

% Tytospan can (unintendedly) change the duration, starttime, and endtime
% of the file if the user is not paying attention. This also changes the
% curvedata.spikecounts and can affect the analysis. 
% This function remakes the spikecounts paradigm to avoid these problems by
% taking the start and endtimes from the more robust curvesettings.stim

% Amarins Heeringa, March 16 2018, Koeppl Lab

% get real timing measures
tlen = curvesettings.stim.Duration;
internaldelay = 2.76192; % fixed internal delay introduced by TDT hardware
StartTime = curvesettings.stim.Delay + internaldelay;
EndTime = StartTime + tlen;

[in,jn]=size(curvedata.spike_times);
nwSpikCnt=nan(in,jn); nwSpikRat=nan(in,jn);
for i=1:in
    for j=1:jn
        spks=curvedata.spike_times{i,j};
        spks2=spks(spks>StartTime & spks<EndTime);
        nwSpikCnt(i,j)=length(spks2);
        nwSpikRat(i,j)=length(spks2)/(tlen/1000);
    end
end
