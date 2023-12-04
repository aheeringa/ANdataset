% Script to plot average spike waveforms, rate-level functions, and PSTH
% shape. These outcomes are used to determine whether fibers derived from
% the AN bundle. 
% This script uses a single data structure of one animal (exp). 

% By: Amarins Heeringa, October 4, 2023

%% load in the struct

clear; close all;

id = input('Enter the animal ID of the struct you want to see: \n',"s");
load(sprintf('%s.mat',id))

fprintf('This animal (%s) is %.1f months old and has an ABR threshold of %d dB SPL\n',exp.info.sex, exp.info.age/(365/12), exp.info.ABR.threshold)

%% settings for the plots

fs = 22; % fontsize
sym = 'o'; col = 'k';
sz = 70; % markersize

%% make a quick plot of the BF and thresholds

nunit = length(exp.data);
bf = nan(1,nunit);
thr = nan(1,nunit);
srT = nan(4,nunit);

for i = 1:nunit
    if ~isempty(exp.data(i).BF)
    bf(i) = exp.data(i).BF.analysis.bf;
    srT(4,i) = exp.data(i).BF.analysis.sr;
    end
    if ~isempty(exp.data(i).RLF)
    thr(i) = exp.data(i).RLF.analysis.threshold;
    srT(3,i) = exp.data(i).RLF.analysis.sr;
    end
    if ~isempty(exp.data(i).SR)
        srT(1,i) = exp.data(i).SR.analysis.sr;
    end
    if ~isempty(exp.data(i).PH)
        srT(2,i) = exp.data(i).PH.analysis.sr;
    end
end

sr = nan(1,nunit);
for i = 1:nunit
    if ~isnan(srT(1,i))
        sr(1,i) = srT(1,i);
    elseif ~isnan(srT(2,i))
        sr(1,i) = srT(2,i);
    elseif ~isnan(srT(3,i))
        sr(1,i) = srT(3,i);
    elseif ~isnan(srT(4,i))
        sr(1,i) = srT(4,i);
    end
end

% set limitations and ticks of the two axes
xl = [0.3 20]; 
yl = [0 100];
xt = [0.2 0.5 1 2 5 10 20];

lsr = find(sr<18);
hsr = find(sr>= 18);

figure;
scatter(bf(hsr)/1000,thr(hsr),sz,col,sym,'filled')
hold on
scatter(bf(lsr)/1000,thr(lsr),sz,col,sym,'LineWidth',1)
xlabel('Best Frequency (kHz)')
ylabel('Threshold (dB SPL)')
set(gca,'FontSize',fs,'XScale','log','XTick',xt)
xlim(xl); ylim(yl);
legend('High SR','Low SR','Location','northeast')

%% plot the average waveforms of the BF files

% i = 12 is good for AHG68 as an example
i = 12;
% for i = 1:nunit
    checkAN(exp.data(i).BF.curvedata, exp.data(i).BF.curveresp, exp.data(i).BF.curvesettings)
% end

%% plot all the RLFs of this gerbil

xl = [5 85];
yl = [0 400];

figure;

for i = 1:nunit
    if ~isempty(exp.data(i).RLF)
    rate = exp.data(i).RLF.analysis.rates;
    std = exp.data(i).RLF.analysis.stdevs;
    levs = exp.data(i).RLF.analysis.levels;

    % errorbar(levs,rate,std,'LineWidth',1.5)
    plot(levs,rate,'LineWidth',1.5)
    hold on
    end
end

set(gca,'FontSize',fs)
xlabel('Level (dB SPL)')
ylabel('Rate (spikes/s)')
xlim(xl)
% ylim(yl)
box off

%% make psth

close all
TestLevel = 20;
i = 12;

% for i = 1:nunit
    makePSTH(exp,i,TestLevel)
% end


