function [threshold, sr, rates, stdevs, levels, freq] = RLFextract_func(curvedata,curvesettings)
%% 
% input:
%   curvedata: a curvedata struct, containing the spike times
%   curvesettings: a curvesettings struct, containing the recording metadata
% output:
%   threshold: the estimated threshold based on the following criteria: the
%       first level where the rate is higher than SR + 1.2*std(SR) and 
%       higher than 15 spikes/s (since SR of low-SR fibers is often 0, so 
%       any one spike will trigger threshold)
%   sr: the spontaneous rate in spikes/s, based on the firing rate during 
%       the silent trials
%   rates: a vector of the mean spike rates during the presentation of all
%       stimulus levels. 
%   stdevs: a vector of the standard deviation of the mean spike rates during
%       the presentation of all stimulus levels. 
%   levels: the levels at which the rate-level curve was collected,
%       derives from the curvesettings directly
%   freq: the frequency at which the rate-level curve was collected,
%       derives from the curvesettings directly

% By: Amarins Heeringa

% Altered by Amarins to make it into a function, Sept 5, 2017
% Add getSpikeCounts.m to overcome problems with incorrect window, Amarins,
% April, 2018

%%
% delete a repetition from the data file
% ri = 10;
% curvedata.spike_times(:,ri)=[];
% curvedata.isspont(:,ri)=[];
% curvedata.spike_counts(:,ri)=[];

% do some checks on the file to see if it was collected using the standard
% settings of the Koeppl lab
if curvesettings.stim.Duration ~= 50
    fprintf('Duration is not 50 ms! But %d ms ..\n',curvesettings.stim.Duration)
    return
end
if curvesettings.stim.ISI ~= 150
    fprintf('ISI is not 150 ms! But %d ms ..\n',curvesettings.stim.ISI)
    return
end
if curvesettings.stim.Ramp ~= 5
    fprintf('Ramp is not 5 ms! But %d ms .. \n',curvesettings.stim.Ramp)
    return
end
if curvesettings.tdt.AcqDuration ~= 80
    fprintf('Acquisition duration is not 80 ms! But %d ms ..\n',curvesettings.tdt.AcqDuration)
    %     return
end

[rates, stdevs, levels, msp, ssp, freq] = getDrivenRatesRLF(curvedata,curvesettings);

% when this is not a RLF file, output NaN for each variable
if isempty(rates) || isempty(levels)
    disp('Not a valid RLF file')
    freq=NaN; threshold=NaN; sr=NaN; rates=NaN; stdevs=NaN; levels=NaN;
else

    x = levels;
    m = rates;
    s = stdevs;
    sr = msp; 

    % determine threshold, saturation and dynamic range
    threshold = x(max(min(find(m>(msp+1.2*ssp))),min(find(m>15)))); % new, find out if it works
    if isempty(threshold)
        threshold=NaN;
    end

    go2 = 1;
    while go2 == 1 % loop until the threshold is ok

        % plot the data
        figure;
        errorbar(x,m,s,'-bo','LineWidth',1);
        hold on;
        errorbar(min(x)-10,msp,ssp,'or','MarkerEdgeColor','k','MarkerFaceColor','r')
        yl=ylim;
        xl=xlim; xl(1) = xl(1) - 2; xl(2) = xl(2) + 2;
        line([threshold threshold],yl,'LineStyle','--','Color','k')
        xlabel('Stimulus level (dB SPL)');
        ylabel('Spike rate (spikes/s)');
        set(gca,'FontSize',18)
        legend('raw data','spont rate','threshold','location','southeast')
        ylim(yl)
        xlim(xl)

        newtry2 = input('Do you want to change the threshold? (y/n) ','s');
        if strcmp(newtry2,'y')
            go2=1;
            fprintf('The old threshold was %4.0d dB\n',threshold)
            threshold = input('Enter the new threshold in dB: ');
            close gcf
        else
            go2 = 0;
        end
    end
end
end