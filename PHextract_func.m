function [freq,lvlVec,vs,prob,vsVec,pVec,erVec,msp] = PHextract_func(curvedata, curvesettings)
%%
% input:
%   curvedata: a curvedata struct, containing the spike times
%   curvesettings: a curvesettings struct, containing the recording metadata
% output:
%   freq: the frequency at which the rate-level curve was collected,
%       derives from the curvesettings directly
%   lvlVec: the levels at which the rates and vector strengths were collected,
%       derives from the curvesettings directly
%   vs: the vector strength value of all levels combined
%   prob: the p-value for the vs value
%   vsVec: the vector strength values for each level separately in a vector
%       corresponding to the lvlVec
%   pVec: the p-values for each vector strength value at each level
%   erVec: the rate at each stimulus level
%   msp: the spontaneous rate in spikes/s, based on the firing rate during
%       the silent trials

% By: Amarins Heeringa

%%
% delete a repetition from the data file
% ri = 25:40;
% curvedata.spike_times(:,ri)=[];
% curvedata.isspont(:,ri)=[];
% curvedata.spike_counts(:,ri)=[];
% curvedata.depvars_sort(:,ri,:)=[];

% Determine some variables
tlen = curvesettings.stim.Duration;
dur = curvesettings.tdt.AcqDuration/1000;
StartTime = curvesettings.stim.Delay;
EndTime = StartTime + tlen;
binSize = 0.001; % set binsize for PSTH
tmpfreq = sort(unique(vertcat(curvesettings.stimcache.Freq{:})));
if curvesettings.curve.Spont
    freq = tmpfreq(2); % tmpfreq(1) should be spont (-99999)
else
    freq = tmpfreq(1);
end
tmplvl = sort(unique(vertcat(curvesettings.stimcache.ABI(:))));
if curvesettings.curve.Spont
    lvlVec = tmplvl(2:end);
else
    lvlVec = tmplvl;
end
x1 = curvedata.depvars_sort(:,1,1);


if curvesettings.stimcache.nloopvars == 0 % in case of only one level

    % plot rasterplot
    figure;
    subplot(2,2,1)
    [in,jn]=size(curvedata.spike_times);
    cnt=1; spiks=[]; spksEv=[]; NspkEv=[]; spksSr=[]; NspkSr=[];
    for i=1:in
        for j=1:jn
            if ~curvedata.isspont(i,j)
                spks=curvedata.spike_times{i,j};
                spiksEv=spks(spks>StartTime & spks<EndTime);
                spksEv=[spksEv spiksEv];
                NspkEv=[NspkEv length(spiksEv)];
                spiks=[spiks spks];
                scatter(spks,(ones(1,length(spks)))*cnt,'.k')
                hold on
                cnt=cnt+1;
            else
                spks=curvedata.spike_times{i,j};
                %                 spiksSr=spks(spks>StartTime & spks<EndTime);
                spiksSr=spks;
                spksSr=[spksSr spiksSr];
                NspkSr=[NspkSr length(spiksSr)];
            end
        end
    end
    ylabel('Trial #'); xlabel('Time (ms)')
    title(sprintf('Rasterplot, Freq = %d Hz',freq))
    xlim([0 100])
    box on

    % plot psth
    subplot(2,2,2);
    nbin = round(dur/binSize);
    [N,cent] = hist(spiks,nbin); N=(N/cnt)/binSize;
    bar(cent,N,0.9,'k');
    xlabel('Time (ms)'); ylabel('Rate (spikes/s)')
    xlim([0 100])
    title(sprintf('PSTH, binwidth %4.2f ms',binSize*1000))

    % plot rate-level function errorbar
    evR=mean(NspkEv./(tlen/1000));
    stevR=std(NspkEv./(tlen/1000));
    subplot(2,3,4)
    errorbar(lvlVec,evR,stevR,'ko-');
    hold on;
    if ~isempty(NspkSr)
        %         sR=mean(NspkSr./(tlen/1000));
        %         stsR=std(NspkSr./(tlen/1000));
        sR=mean(NspkSr./dur); % adjusted new SR, April 2021, ANH
        stsR=std(NspkSr./dur);
        errorbar(lvlVec-10,sR,stsR,'or','MarkerEdgeColor','k','MarkerFaceColor','r')
    else
        sR=NaN;
    end
    hold off;
    xlabel('Level (dB SPL)');
    ylabel('Rate (spikes/s)');
    title('rate');
    yl=ylim;
    ylim([0 yl(2)])
    % xl=xlim;
    xlim([lvlVec-15 lvlVec+5])

    % calculate vector strength
    nbin = 60;
    [vs,prob,ntot,xbin,nspikes]=TytoView_calcVS(spksEv,freq,nbin);

    % now plot
    subplot(2,3,6)
    bar(xbin,nspikes,1,'k');
    xlim([0,1]);
    xlabel('phase (cycle)'); ylabel('spike count');
    title('phase histogram');
    yl=ylim;
    text(0.6,yl(2)*0.9,sprintf('VS = %.2f',vs))

    % plot both VS
    subplot(2,3,5)
    scatter(lvlVec,vs,'ko');
    hold on;
    if ~isempty(NspkSr) % plot VS of SR
        nbin=60;
        [vsS,probS,ntotS,xbinS,nspikesS]=TytoView_calcVS(spksSr,freq,nbin);
        plot(lvlVec-10,vsS,'ko-','MarkerEdgeColor','k','MarkerFaceColor','r');
        if ~isnan(probS)
            if probS<0.001
                text(lvlVec-10,vsS-0.15,'*','FontSize',20,'HorizontalAlignment','center')
                scatter(lvlVec-10,vsS,'ko','filled')
            end
        end
    end
    if ~isnan(prob)
        if prob<0.001
            text(lvlVec,vs-0.15,'*','FontSize',20,'HorizontalAlignment','center')
            scatter(lvlVec,vs,'ko','MarkerEdgeColor','k','MarkerFaceColor','b')
        end
    end
    xlabel('Level (dB SPL)');
    ylabel('vector strength'); ylim([0 1]); %xlim(xl);
    title('vector strength');
    xlim([lvlVec-15 lvlVec+5])
    box on

    lvlVec=x1;
    vsVec=[vs,vsS];
    pVec=[prob,probS];
    erVec=evR;
    msp=sR;

elseif curvesettings.stimcache.nloopvars ==1 % in case they loop over the levels

    % plot rasterplot
    figure;
    subplot(2,2,1)
    [in,jn]=size(curvedata.spike_times);
    %     jn = 33; % temp line to exclude last 7 repetitions
    cnt=1; lvlV=[]; spiks=[]; NspkEv=[]; NspkSr=[];
    for i=1:in
        for j=1:jn
            if ~curvedata.isspont(i,j)
                spks=curvedata.spike_times{i,j};
                spiksEv=spks(spks>StartTime & spks<EndTime);
                NspkEv=[NspkEv; length(spiksEv)];
                spiks=[spiks spks];
                scatter(spks,(ones(1,length(spks)))*cnt,'.k')
                hold on
                cnt=cnt+1;
                lvlV=[lvlV; curvedata.depvars_sort(i,1) cnt-1];
            else
                spks=curvedata.spike_times{i,j};
                %                 spiksSr=spks(spks>StartTime & spks<EndTime);
                spiksSr=spks;
                NspkSr=[NspkSr; length(spiksSr)];
            end
        end
    end
    lvls=unique(lvlV(:,1));
    newLvl=nan(length(lvls),2);
    for i=1:length(lvls)
        ind=max(find(lvlV(:,1)==lvls(i)));
        newLvl(i,:)=lvlV(ind,:);
    end
    set(gca,'YTick',newLvl(:,2),'YTickLabels',newLvl(:,1))
    ylabel('Level (dB SPL)'); xlabel('Time (ms)')
    title(sprintf('Rasterplot, F = %d Hz',freq))
    xlim([0 100])
    box on

    % plot psth
    subplot(2,2,2);
    nbin = round(dur/binSize);
    [N,cent] = hist(spiks,nbin); N=(N/cnt)/binSize;
    bar(cent,N,0.9,'k');
    xlabel('time (ms)'); ylabel('Rate (spikes/s)')
    xlim([0 100])
    title(sprintf('PSTH, binwidth %4.2f ms',binSize*1000))

    % plot evoked rates
    m=nan(1,length(lvlVec)); s=nan(1,length(lvlVec));
    for i=1:length(lvlVec)
        m(1,i)=mean(NspkEv(lvlV(:,1)==lvlVec(i))./(tlen/1000));
        s(1,i)=std(NspkEv(lvlV(:,1)==lvlVec(i))./(tlen/1000));
    end
    subplot(2,3,4)
    errorbar(lvlVec, m, s, 'ko-');
    hold on;
    if ~isempty(NspkSr)
        %         msp=mean(NspkSr./(tlen/1000));
        %         ssp=std(NspkSr./(tlen/1000));
        msp=mean(NspkSr./dur);
        ssp=std(NspkSr./dur);
        errorbar(min(lvlVec)-10,msp,ssp,'or','MarkerEdgeColor','k','MarkerFaceColor','r')
    else
        msp=NaN;
    end
    hold off;
    xlabel('Level (dB SPL)');
    ylabel('Rate (spikes/s)');
    title('rate');
    erVec = m;
    yl=ylim;
    ylim([0 yl(2)])
    xlim([lvlVec(1)-5 lvlVec(end)+5]);


    % calculate vector strength
    nbin = 60;
    allspiketimes0 = horzcat(curvedata.spike_times{(curvedata.isspont==0)});
    allspiketimes = allspiketimes0( allspiketimes0>StartTime & allspiketimes0<EndTime ); % added Nov 2015
    [vs,prob,ntot,xbin,nspikes]=TytoView_calcVS(allspiketimes,freq,nbin);

    % now plot
    subplot(2,3,6)
    bar(xbin,nspikes,1,'k');
    xlim([0,1]);
    xlabel('phase (cycle)'); ylabel('spike count');
    title('phase histogram');
    yl=ylim;
    text(0.6,yl(2)*0.9,sprintf('VS = %.2f',vs))


    % get all VS for different levels
    if strcmp(curvesettings.stimcache.loopvars{1}, 'ABI')

        % calculate vector strength
        nbin = 60;
        vsVec=nan(1,length(x1));
        pVec=nan(1,length(x1));

        for i1=1:length(x1)

            % calculate VS
            allspiketimes0=horzcat(curvedata.spike_times{(curvedata.depvars_sort(:,:,1)==x1(i1))});
            allspiketimes=allspiketimes0(allspiketimes0>StartTime&allspiketimes0<EndTime); % added Nov 2015
            [vsTemp,pTemp,ntot,xbin,nspikes]=TytoView_calcVS(allspiketimes,freq,nbin);
            if ~isempty(pTemp)
                pVec(1,i1)=pTemp;
            end
            vsVec(1,i1)=vsTemp;
            
        end

        nonSpont = (curvedata.isspont(:,1)==0);
        Spont = (curvedata.isspont(:,1)==1);
        vsSR = vsVec(Spont); pSR = pVec(Spont); levSR=min(x1(nonSpont))-10;
        vsER = vsVec(nonSpont); pER = pVec(nonSpont); levER=x1(nonSpont);
        subplot(2,3,5)
        plot(levER,vsER,'ko-');
        hold on;
        if ~isempty(NspkSr) % plot VS of SR
            plot(levSR,vsSR,'ko-','MarkerEdgeColor','k','MarkerFaceColor','r');
            if ~isnan(pSR)
                if pSR<0.001
                    text(levSR,vsSR-0.15,'*','FontSize',20,'HorizontalAlignment','center')
                    scatter(levSR,vsSR,'ko','filled')
                end
            end
        end
        for iprob=1:length(pER)
            if ~isnan(pER(iprob))
                if pER(iprob)<0.001
                    text(levER(iprob),vsER(iprob)-0.15,'*','FontSize',20,'HorizontalAlignment','center')
                    scatter(levER(iprob),vsER(iprob),'ko','MarkerEdgeColor','k','MarkerFaceColor','b')
                end
            end
        end
        xlabel('Level (dB SPL)');
        ylabel('vector strength'); ylim([0 1]);
        % xlim(xl);
        xlim([lvlVec(1)-5 lvlVec(end)+5]);
        title('vector strength');
    end


    lvlVec=x1;

else
    disp('not a valid PH file')
end
drawnow
