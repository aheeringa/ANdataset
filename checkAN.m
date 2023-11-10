function checkAN(curvedata,curveresp,curvesettings)

% This script plots the interspike interval information, the first 300
% waveforms and the average +/- std waveform. 
% The ISI information gives an indication whether a single or a multi unit
% was recorded from, there should be close to no ISIs < 1 ms and no ISIs 
% < 0.6 ms.
% The waveform plots reveals if the spikes have a prepotential, which 
% would mean that they were recorded from AVCN neurons, instead of from the 
% auditory nerve. 
% Furthermore, one could look at the click latency to see whether a single
% unit derived from AN or AVCN, see eg Fitzgerald et al., 2001.

% By: Amarins Heeringa, Koeppl Lab, July 2018

Fs = curvesettings.Fs(1,1); 

T = length(curveresp{1,1})/Fs;
t = 1/Fs:1/Fs:T;
yspik = zeros(1,length(curvedata.spike_times{1,1})); 
yspik = yspik + max(curveresp{1,1})*1.1;

% plot the first datatrace and spikes
figure;
subplot(2,2,1)
plot(t*1000,curveresp{1,1},'b')
hold on
scatter(curvedata.spike_times{1,1},yspik,[],'r','*')
xlabel('time (ms)')
ylabel('Amplitude (mV)')
title('First data trace, unfiltered')

% prepare variables to get single waveforms
[a,b,c] = size(curveresp);
ms1 = round(Fs*0.00125); % amount of datapoints for 1.25 ms
lwv = (ms1*2)+1; % total length waveform in datapoints
spks = nan(lwv,400000); % pre-allocate waveform matrix
cnt = 1;
isis = []; 

% get all the waveforms
for an = 1:a
    for bn = 1:b
        for cn = 1:c
            spT = curvedata.spike_times{an,bn,cn}/1000; % spike times in ms
            isis = [isis diff(spT)]; % get the interspike intervals
            spDP = round(spT*Fs); % spike times in data points
            resp = curveresp{an,bn,cn}; % get data trace
            for sp = 1:length(spT) % loop through spikes of that trial
                if spDP(sp)>ms1 && spDP(sp)+ms1<length(resp) % only include spikes not too close to the edges of the response
                    spks(:,cnt) = resp(spDP(sp)-ms1:spDP(sp)+ms1); % get the waveform
                    cnt = cnt+1;
                end
            end
        end
    end
end

% plot ISIs and info
subplot(2,2,2)
edges = [0:1:50];
h=histogram(isis*1000,edges,'FaceAlpha',0.8,'FaceColor',[11 154 51]/256);
xlabel('Interspike intervals (ms)')
ylabel('Count')
title('Interspike interval histogram')
yl=ylim;
text(23,yl(2)*0.9,sprintf('total ISIs = %d',length(isis)));
text(23,yl(2)*0.82,sprintf('ISIs < 1 ms = %d',sum(isis<0.001)));
text(23,yl(2)*0.74,sprintf('ISIs < 0.6 ms = %d',sum(isis<0.0006)));
text(23,yl(2)*0.66,sprintf('ISI median = %0.1f ms',median(isis)*1000));
text(23,yl(2)*0.58,sprintf('min ISI = %0.1f ms',min(isis)*1000));

spks=spks(:,1:cnt-1);
twv = linspace(-1.25,1.25,lwv);

% plot the first 300 waveforms
subplot(2,2,3)
plot(twv,spks(:,1),'k','LineWidth',0.3)
hold on
if cnt>300
    n=300;
else
    n=cnt-1;
end
for i=2:n
    plot(twv,spks(:,i),'k','LineWidth',0.3)
end
title('First 300 spike waveforms')
xlabel('Time rel. spike peak (ms)')
ylabel('Amplitude (mV)')
set(gca,'XTick',[-1 0 1],'XTickLabel',{'-1 ms','0','1 ms'})
xlim([-1.3 1.3]); yl=ylim;

% plot the average waveform
% mWV = mean(spks,2); 
% stdWV = std(spks,[],2);
mWV = median(spks,2); 
ci1 = quantile(spks,0.025,2);
ci2 = quantile(spks,0.975,2);
subplot(2,2,4)
% figure;
plot(twv,mWV,'r','LineWidth',2);
hold on
% plot(twv,mWV+2*stdWV,'r',twv,mWV-2*stdWV,'r','LineWidth',0.3)
plot(twv,ci1,'r',twv,ci2,'r','LineWidth',0.3)
patch([twv fliplr(twv)],[(ci2)' fliplr((ci1)')],'r','FaceAlpha',0.25,'EdgeColor','none')
text(-0.9,max(mWV),sprintf('n = %d',cnt))
xlabel('Time rel. spike peak (ms)')
ylabel('Amplitude (mV)')
title('Waveform median +/- 95% CI')
set(gca,'XTick',[-1 0 1],'XTickLabel',{'-1 ms','0','1 ms'})
xlim([-1.3 1.3]); ylim(yl)
