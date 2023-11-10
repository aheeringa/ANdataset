function [m,s,x,msp,ssp,freq] = getDrivenRatesRLF(curvedata,curvesettings)
%% 
% input:
%   curvedata: a curvedata struct, containing the spike times
%   curvesettings: a curvesettings struct, containing the recording metadata
% output:
%   m: a vector of the mean spike rates during the presentation of all
%       stimulus levels. 
%   s: a vector of the standard deviation of the mean spike rates during
%       the presentation of all stimulus levels. 
%   x: the levels at which the rate-level curve was collected,
%       derives from the curvesettings directly
%   msp: the mean spontaneous rate in spikes/s, based on the firing rate during 
%       the silent trials
%   ssp: the standard deviation of the mean spontaneous rate in spikes/s, 
%       based on the firing rate during the silent trials
%   freq: the frequency at which the rate-level curve was collected,
%       derives from the curvesettings directly

% By: Amarins Heeringa

%%
% delete a repetition from the data file
% ri = 10;
% curvedata.spike_times(:,ri)=[];
% curvedata.isspont(:,ri)=[];
% curvedata.spike_counts(:,ri)=[];

[~,~,~,~,nwSpikRat]=getSpikeCounts(curvedata,curvesettings);

% loop variables
x1 = curvedata.depvars_sort(:,1,1);
x2 = curvedata.depvars_sort(:,1,2);

% calculate mean and std of spike rates 
% y = curvedata.spike_counts * 1000 / tlen; 
y = nwSpikRat;
m1 = mean(y,2);  % average wrt 2nd index
s1 = std(y,1,2); % std wrt 2nd index 

% data for non-spont trials
j = (curvedata.isspont(:,1)==0); % find out non-spont trials
x = x1(j);
m = m1(j);
s = s1(j);

% % data for spont trials - old version
% jsp = (curvedata.isspont(:,1)==1); % find out spont trials
% msp = m1(jsp);
% ssp = s1(jsp);
 
% data for spont trials - new version, takes into account whole acquisition
% duration
jsp = (curvedata.isspont(:,1)==1); % find out spont trials
[~,jn]=size(curvedata.spike_times);
dur = curvesettings.tdt.AcqDuration/1000; sr = [];
for i = 1:jn
    spks = curvedata.spike_times{jsp,i};
    sr = [sr length(spks)/dur];
end
msp = mean(sr);
ssp = std(sr);

% extract frequency from the data structure
tmpfreq = sort( unique( vertcat( curvesettings.stimcache.Freq{:} ) ) );
if curvesettings.curve.Spont
    freq = tmpfreq(2); % tmpfreq(1) should be spont (-99999)
else
    freq = tmpfreq(1);
end
end
