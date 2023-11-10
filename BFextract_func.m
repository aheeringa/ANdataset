function [bf,sr,bandwidth,level] = BFextract_func(curvedata,curvesettings)
%% 
% input:
%   curvedata: a curvedata struct, containing the spike times
%   curvesettings: a curvesettings struct, containing the recording metadata
% output:
%   bf: the estimated best frequency in Hz, based on a curve fitted to the data
%   sr: the spontaneous rate in spikes/s, based on the firing rate during the silent trials
%   bandwidth: the bandwidth in Hz at 50% of the maximum firing rate, is
%       NaN if the frequency range was not wide enough
%   level: the level at which the frequency-response curve was collected,
%       derives from the curvesettings directly

% By: Go Ashida, Roberta Aralla, and Amarins Heeringa

% Altered by Amarins to make it into a function, Sept 5, 2017
% Add getSpikeCounts.m to overcome problems with incorrect window, Amarins,
% April, 2018

%%
% delete a repetition from the data file
% ri = 5;
% curvedata.spike_times(:,ri)=[];
% curvedata.isspont(:,ri)=[];
% curvedata.spike_counts(:,ri)=[];

[tlen,StartTime,EndTime,nwSpikCnt,nwSpikRat]=getSpikeCounts(curvedata,curvesettings);

% loop variables
x1 = curvedata.depvars_sort(:,1,1);
x2 = curvedata.depvars_sort(:,1,2);

% calculate mean and std of spike rates 
y = nwSpikRat;
m1 = mean(y,2);  % average wrt 2nd index
s1 = std(y,1,2); % std wrt 2nd index 

% data for non-spont trials
j = (curvedata.isspont(:,1)==0); % find out non-spont trials
x = x1(j);
m = m1(j);
s = s1(j);

% data for spont trials
jsp = (curvedata.isspont(:,1)==1); % find out spont trials
[~,jn]=size(curvedata.spike_times); 
dur = curvesettings.tdt.AcqDuration/1000; sr = [];
for i = 1:jn
    spks = curvedata.spike_times{jsp,i};
    sr = [sr length(spks)/dur];
end
msp = mean(sr);
ssp = std(sr);
% msp = m1(jsp);   % old version by GA, revised Friday, April 9, 2021
% ssp = s1(jsp);

v1 = sort( unique(x1(j) ) );
v2 = max(cell2mat(curvesettings.stimcache.Freq)); % for tone
%v2 = unique(cell2mat(cellfun(@(x){min(x)}, curvesettings.stimcache.Freq))); % for tone and noise (by GA, Jul 2015)

level = max(curvesettings.stimcache.ABI);

if curvesettings.stimcache.nloopvars == 2
    
   % get dependent variables (without spont)
    v2 = sort( unique(x2(j) ) );

    % matrix to store mean rates
    m = zeros(length(v1),length(v2));
    s = zeros(length(v1),length(v2));
    tmpstr = '';
    for i2 = 1:length(v2)
        m(:,i2) = m1( x2==v2(i2) );
        s(:,i2) = s1( x2==v2(i2) );
        tmpstr = [ tmpstr ', ' num2str(v2(i2))];
    end
end

%Create a ISI histogram
vertspikes = curvedata.spike_times(:);
ISI=cell(length(vertspikes),1);
for n=1:length(vertspikes)
    TEMPSPIKES=cell2mat(vertspikes(n));
    TEMPISI=length(TEMPSPIKES-1);
    for i=2:length(TEMPSPIKES)
        TEMPISI(i-1)=TEMPSPIKES(i)-TEMPSPIKES(i-1);
    end
    ISI{n}=TEMPISI';
end
    
ISI=cell2mat(ISI);
if numel(ISI)==0
    disp('This recording has not enough spikes, unable to plot histogram')
%else
    %hist(ISI,max(ISI)*10)
end

% ShortISI=find(ISI<1);
ShortISI = find(ISI<0.6); % new value from Heil paper, ANH 04-2018
SHORTISIS=((length(ShortISI)/length(ISI))*100)


%%%% ======================---------------------===================%%%%%%
% add here a section from BFanalysis of Roberta that calculates the BF and bandwidth
% comes from file: A_BFanalysis_sort_BW.m

lw=1.5; % line width

% Fitting option
fit_option = fitoptions('Method','SmoothingSpline','SmoothingParam',0.00002); %small values reduce the fitting accuracy

% plot
figure;
depvars=v1;
meanresponse=m;
errorbar(depvars,meanresponse,s,'ok-','LineWidth',lw); % RA (4.5.2017)
hold on
% fit a curve and search for peaks
[curve, goodness, output] = fit(depvars,meanresponse,'smoothingspline',fit_option);
% goodness
plot(min(depvars):max(depvars),curve(min(depvars):max(depvars)),'r-','LineWidth',lw);
yFitted = feval(curve,min(depvars):1:max(depvars));
[ypk1,idx] = findpeaks(yFitted);
xpk1 = idx+min(depvars);
if ~isempty(ypk1)
    [resp,bfInd]=max(ypk1);
    bf=xpk1(bfInd);
else
    bf = NaN;
end

sr = msp;

%Calculation for the bandwidth
hold on
x_axis_diff = min(depvars):max(depvars);
y_axis_diff = curve(min(depvars):max(depvars));
if ~isempty(ypk1)
    [bandwidth,oor] = bandwidthBF(x_axis_diff, y_axis_diff, ypk1, msp);
    if oor==1 % out of range
        bandwidth=NaN;
    end
else
    bandwidth=NaN;
end

% plot spont
if ~isempty(msp)
    errorbar([min(x) max(x)], [msp msp], [ssp ssp], 'ko--','LineWidth',lw)
end

% plot calculated bf
if ~isempty(bf) && ~isnan(bf)
    yl=ylim;
    line([bf,bf],[yl(1),resp],'Color','m','LineStyle','--','LineWidth',lw)
    bfTxt=['BF = ',num2str(bf),' Hz'];
    text(bf+70,yl(1)+10,bfTxt)
end
grid on;
set(gca,'fontsize',15)
xlabel('Frequency (Hz)');
ylabel('Response (spikes/s)');

end