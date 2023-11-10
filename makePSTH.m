function makePSTH(andata,i,TestLevel)
%%
% input:
%   andata: the struct that contains the data that needs to be plotted
%   i: the unit number in the struct that needs to be tested
%   TestLevel: the level in dB above threshold at which the PSTH needs to
%   be plotted
% output:
%   a figure with the PSTHs at TestLevel above threshold

% By: Amarins Heeringa

%%
figure;
binSize = 0.001; % set binsize for PSTH

if ~isempty(andata.data(i).BF)
    bf = andata.data(i).BF.analysis.bf;
else
    bf = [];
end
if ~isempty(andata.data(i).RLF)
    thr = andata.data(i).RLF.analysis.threshold;
else
    thr = [];
end

% plot PSTH of BF file as closest to BF as possible
if ~isempty(andata.data(i).BF)
    filename = andata.data(i).BF.filename;
    
    % BF = load(strcat(Direc,filename));
    BF = andata.data(i).BF;
    BFdur = BF.curvesettings.tdt.AcqDuration/1000;
    BFdelay = BF.curvesettings.stim.Delay;
    idx2 = find(abs(BF.curvedata.depvars_sort(:,1,1) - bf)==min(abs(BF.curvedata.depvars_sort(:,1,1) - bf)));
    if length(idx2) > 1
        idx2 = min(idx2);
    end
    reps = length(BF.curvedata.depvars_sort(1,:,1));
    spiksBF=[];
    for r = 1:reps
        spks = BF.curvedata.spike_times{idx2,r};
        spiksBF=[spiksBF spks];        
    end
    
    subplot(2,2,1)
    nbin = round(BFdur/binSize);
    [N,cent] = hist(spiksBF,nbin); N=(N/reps)/binSize;
    bar(cent,N,0.9,'k');
    xlabel('time (ms)'); ylabel('spike rate (sp/s)')
    xlim([0 100])
    
    level = unique(andata.data(i).BF.curvesettings.stimcache.ABI); level = level(end);
    title(sprintf('BF = %d Hz, PSTH at %d Hz\n Thr = %d dB, Level = %d dB',...
        bf, BF.curvedata.depvars_sort(idx2,1,1), thr, level))
else
    BF = [];
end

% plot PSTH of RLF file at TestLevel above threshold
if ~isempty(andata.data(i).RLF) && ~isnan(thr)
    filename = andata.data(i).RLF.filename;
    
    % ABI = load(strcat(Direc,filename));
    ABI = andata.data(i).RLF;
    ABIdur = ABI.curvesettings.tdt.AcqDuration/1000;
    ABIdelay = ABI.curvesettings.stim.Delay;
    lev = thr + TestLevel;
    idx2 = find(abs(ABI.curvedata.depvars_sort(:,1,1) - lev)==min(abs(ABI.curvedata.depvars_sort(:,1,1) - lev)));
    if length(idx2) > 1
        idx2 = max(idx2);
    end
    
    reps = length(ABI.curvedata.depvars_sort(1,:,1));
    spiksABI=[];
    for r = 1:reps
        spks = ABI.curvedata.spike_times{idx2,r};
        spiksABI=[spiksABI spks];        
    end
    
    subplot(2,2,2)
    nbin = round(ABIdur/binSize);
    [N,cent] = hist(spiksABI,nbin); N=(N/reps)/binSize;
    bar(cent,N,0.9,'k');
    xlabel('time (ms)'); ylabel('spike rate (sp/s)')
    xlim([0 100])
    
    title(sprintf('BF = %d Hz, PSTH at %d Hz\n Thr = %d dB, Level = %d dB',...
        bf,andata.data(i).RLF.analysis.frequency,thr,ABI.curvedata.depvars_sort(idx2,1,1)))
else
    ABI = [];    
end

% plot PSTH of PH file at TestLevel above threshold
if ~isempty(andata.data(i).PH) && ~isempty(thr)
    filename = andata.data(i).PH.filename;
    
    % PH = load(strcat(Direc,filename));
    PH = andata.data(i).PH;
    PHdur = PH.curvesettings.tdt.AcqDuration/1000;
    PHdelay = PH.curvesettings.stim.Delay;
    lev = thr + TestLevel;
    idx2 = find(abs(PH.curvedata.depvars_sort(:,1,1) - lev)==min(abs(PH.curvedata.depvars_sort(:,1,1) - lev)));
    if length(idx2) > 1
        idx2 = max(idx2);
    elseif isempty(idx2)
        idx2 = 1;
    end
    
    reps = length(PH.curvedata.depvars_sort(1,:,1));
    spiksPH=[];
    for r = 1:reps
        spks = PH.curvedata.spike_times{idx2,r};
        spiksPH=[spiksPH spks];        
    end
    
    subplot(2,2,3)
    nbin = round(PHdur/binSize);
    [N,cent] = hist(spiksPH,nbin); N=(N/reps)/binSize;
    bar(cent,N,0.9,'k');
    xlabel('time (ms)'); ylabel('spike rate (sp/s)')
    xlim([0 100])
    
    title(sprintf('BF = %d Hz, PSTH at %d Hz\n Thr = %d dB, Level = %d dB',...
        bf,andata.data(i).PH.analysis.freq,thr,PH.curvedata.depvars_sort(idx2,1,1)))
else
    PH = [];

end

% Plot PSTH with as much as possible data combined
if isempty(ABI) && isempty(PH) || isempty(BF) && isempty(ABI) || isempty(BF) && isempty(PH)
    % do nothing, more than one file missing
    spiks = [];
elseif isempty(BF) || isempty(ABI) || isempty(PH)
    % only file missing, try to combine other two
    if isempty(BF)
        if ABIdelay == PHdelay && ABIdur == PHdur
            spiks = [spiksABI spiksPH];
            dur = ABIdur;
        else
            spiks = [];
        end
    elseif isempty(ABI)
        if BFdelay == PHdelay && BFdur == PHdur
            spiks = [spiksBF spiksPH];
            dur = BFdur;
        else
            spiks = [];
        end
    elseif isempty(PH)
        if BFdelay == ABIdelay && BFdur == ABIdur
            spiks = [spiksBF spiksABI];
            dur = BFdur;
        else
            spiks = [];
        end
    end
else
    % all files there, try to combine
    if ABIdelay == PHdelay && ABIdur == PHdur
        spiks = [spiksABI spiksPH];
        dur = ABIdur;
        if ABIdelay == BFdelay && ABIdur == BFdur
            spiks = [spiks spiksBF];
            dur = BFdur;
        end
    elseif PHdelay == BFdelay && PHdur == BFdur
        spiks = [spiksPH spiksBF];
        dur = BFdur;
    end
end

if ~isempty(spiks)
subplot(2,2,4)
% figure;
nbin = round(dur/binSize);
[N,cent] = hist(spiks,nbin); N=(N/reps)/binSize;
bar(cent,N,0.9,'k');
xlabel('time (ms)'); ylabel('spike rate (sp/s)')
xlim([0 100])
title('All data together')
end

end