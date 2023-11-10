function [lat_poisson,lat_2bins,fsl_mean,fsl_median,fsl_jit_std,fsl_jit_var,fsl_jit_iqr] = Clickextract_func(curvedata, curvesettings)
%%
% input:
%   curvedata: a curvedata struct, containing the spike times
%   curvesettings: a curvesettings struct, containing the recording metadata
% output:
%   lat_poisson: the click latency calculated based on the poisson surprise
%       method (Chase & Young, PNAS, 2007)
%   lat_2bins: the click latency calculated based on the time when 
%       2 consecutive bins exceed the SR threshold (Köppl, 1997)
%   fsl_mean: the mean of the latencies of all first spike after the click
%       onset.
%   fsl_median: the median of the latencies of all first spike after the 
%       click onset.
%   fsl_jit_std: the jitter calculated by the standard deviation of the 
%       latencies of all first spike after the click onset.
%   fsl_jit_var: the jitter calculated by the variance of the  latencies of 
%       all first spike after the click onset.
%   fsl_jit_iqr: the jitter calculated by the interquartile range of the 
%       latencies of all first spike after the click onset.
% 
% By: Amarins Heeringa

%%

% set some additional variables
nrSpikes = 75; % minimum amount of spikes (75 = one spike in 25% of the trials)
binSize = 0.00005; % set binsize for PSTH in sec
acqLen=curvesettings.tdt.AcqDuration; % length of one trial in ms
dur=acqLen/1000; % duration of one trial in sec
stimDelay = curvesettings.delay; % click delay

% do some checks to see if the files were collected normally
if curvesettings.DAscale ~= 5
    disp('DAscale is not 5 !!')
    return
end
if curvesettings.click.Samples ~= 2
    disp('Number of samples is not 2 !!')
    return
end

% determine number of trials and spikes in each trial
ntrials = length(curvedata.spike_times);
nwSpikCnt=nan(1,ntrials);
for i=1:ntrials
    nwSpikCnt(i)=length(curvedata.spike_times{1,i});
end
trials = 1:ntrials; % vector with all of the trials


% only proceed if the total number of spikes exceeds the minimum amount of spikes 
if sum(nwSpikCnt)>nrSpikes
    
    % plotting some raw data
    figure;
    subplot(2,2,1); plot(nwSpikCnt./dur,'k-o');
    title('Spike rates'); 
    xlabel('Trial number'); 
    ylabel('Rate (spikes/s)');

    % get the rasterplot data and the vector of all the first spikes
    start = 0; stop = 5;
    [rastDataFS, rastDataRest] = getRasterData3(curvedata,trials,stimDelay,start,stop,nwSpikCnt);    
    
    %%%%====---- Calculate FSL ----====%%%%

    % calculate and plot the poisson fsl probability plot
    [latAbs,pd_fsl] = getFSLpoisson(rastDataRest(:,2)',stimDelay);
    lat_poisson = latAbs - stimDelay;
    subplot(2,2,3)
    plot(sort(rastDataRest(:,2)),pd_fsl,'k'); hold on
    set(gca,'yscale','log')
    xlim([0 dur*1000])
    ylim([10^-12 1])
    line(xlim,[10^-6 10^-6],'linestyle',':','color','r')
    xlabel('Time (ms)');
    ylabel('Probability');
    title('Poisson probability plot')
    str={'latency:',sprintf('%0.2f ms',lat_poisson)};
    text(1,10^-10,str,'BackgroundColor','white');

    % make the psth and calculate fsl with 2 bins
    nbin = round(dur/binSize);
    sptimes = horzcat(curvedata.spike_times{trials});
    [N,cent] = hist(sptimes,nbin); N=(N/ntrials)/binSize;
    lat_2bins = getFSL2bins(N,cent,stimDelay); % determine fsl with 2 consecutive bins above max(SFR)
    subplot(2,2,2); bar(cent,N,0.9,'k'); yval=ylim; hold on
    patch([0 stimDelay stimDelay 0],[0 0 yval(2)*0.2 yval(2)*0.2],'blue','FaceAlpha',.3)
    title(sprintf('PSTH, binwidth %4.2f ms',binSize*1000))
    xlabel('Time (ms)'); ylabel('Rate (spikes/s)')
    line([lat_2bins+stimDelay lat_2bins+stimDelay],ylim,'Color','r')
    xlim([0 dur*1000])
    str={'latency:',sprintf('%0.2f ms',lat_2bins)};
    text(1,yval(2)*0.8,str,'BackgroundColor','white');
    text(1,yval(2)*0.25,'Click delay')

    % make a rasterplot and get all first spikes after delay
    fsVec = rastDataFS(:,2)-stimDelay;
    subplot(2,2,4)
    scatter(rastDataRest(:,2),rastDataRest(:,1),5,'k','filled');
    hold on
    scatter(rastDataFS(:,2),rastDataFS(:,1),5,'r','filled');
    xlim([0 dur*1000]); xlabel('time (ms)'); ylabel('trial number');
    title('Rasterplot')
    str = {sprintf('fsl = %0.2f ms',mean(fsVec,'omitnan')),sprintf('jitter = %0.2f ms^{2}',var(fsVec,'omitnan'))};
    text(1,length(trials)*0.2,str,'BackgroundColor','white')%,'EdgeColor','black')
    fsl_mean=mean(fsVec,'omitnan');
    fsl_median=median(fsVec,'omitnan');
    fsl_jit_var=var(fsVec,'omitnan');
    fsl_jit_std=std(fsVec,'omitnan');
    fsl_jit_iqr = iqr(fsVec);
    
else
    disp(filename)
    disp('Not analyzed: spike count below spike number threshold')
    lat_poisson = NaN; lat_2bins = NaN; fsl_mean = NaN; fsl_median = NaN; 
    fsl_jit_var = NaN; fsl_jit_std = NaN; fsl_jit_iqr = NaN;
end
end