% Script to investigate the full, compressed dataset (does not contain raw data). 
% Insert link to dataset: doi:// ...

% This script uses the full dataset to investigate the metadata per animal

% By: Amarins Heeringa, September 5, 2023

%% Load the data and prepare variables

clear; close all

load('all_AN_data_27_04_2023.mat')

% first move the _2 and _3 structs together.
idx = [];
for i = 1:length(all_exp)
    if contains(all_exp(i).animalID,'_1')
        if contains(all_exp(i+1).animalID,'_2')
            all_exp(i).data = [all_exp(i).data all_exp(i+1).data];
            idx = [idx i+1];
        end
        if contains(all_exp(i+2).animalID,'_3')
            all_exp(i).data = [all_exp(i).data all_exp(i+2).data];
            idx = [idx i+2];
        end
        if contains(all_exp(i+2).animalID,'_4')
            all_exp(i).data = [all_exp(i).data all_exp(i+3).data];
            idx = [idx i+3];
        end
    end
end
all_exp(idx) = [];

% pre-allocate variables that will store the metadata per animal
largeNumber = 2000; % make sure this is larger than the sum of all units
ages = nan(1,largeNumber);
sex = nan(1,largeNumber);
abr = nan(1,largeNumber);
bfs = nan(1,largeNumber);
thr = nan(1,largeNumber);
srsT = nan(5,largeNumber);
clk = nan(4,largeNumber);
vss = nan(1,largeNumber);

cnt = 1;

% loop through the animals and store the animal specific metadata
for i = 1:length(all_exp)

    % loop through the units of one animal and store the characteristics of
    % the units
    for j = 1:length(all_exp(i).data)
        
        if ~isempty(all_exp(i).data(j).BF)
            bfs(1,cnt) = all_exp(i).data(j).BF.analysis.bf;
            srsT(4,cnt) = all_exp(i).data(j).BF.analysis.sr;
        end
        if ~isempty(all_exp(i).data(j).RLF)
            thr(1,cnt) = all_exp(i).data(j).RLF.analysis.threshold;
            srsT(3,cnt) = all_exp(i).data(j).RLF.analysis.sr;
        end
        if ~isempty(all_exp(i).data(j).CLICK)
            if length(all_exp(i).data(j).CLICK) == 1
            clk(1,cnt) = all_exp(i).data(j).CLICK.analysis.fsl_mean;
            clk(2,cnt) = all_exp(i).data(j).CLICK.analysis.fsl_median;
            clk(3,cnt) = all_exp(i).data(j).CLICK.analysis.latency_2bins;
            clk(4,cnt) = all_exp(i).data(j).CLICK.analysis.latency_poisson;
            end
        end
        if ~isempty(all_exp(i).data(j).SR)
            srsT(1,cnt) = all_exp(i).data(j).SR.analysis.sr;
        end
        if ~isempty(all_exp(i).data(j).PH)
            srsT(2,cnt) = all_exp(i).data(j).PH.analysis.sr;
            tmp = find(~isnan(all_exp(i).data(j).PH.analysis.prob) & all_exp(i).data(j).PH.analysis.prob<0.001);
            if ~isempty(tmp)
            vss(1,cnt) = max(all_exp(i).data(j).PH.analysis.vs(tmp));
            end
        end

        ages(1,cnt) = all_exp(i).info.age;
        if strcmp(all_exp(i).info.sex,'F')
            sex(1,cnt) = 0;
        elseif strcmp(all_exp(i).info.sex,'M')
            sex(1,cnt) = 1;
        end
        abr(1,cnt) = all_exp(i).info.ABR.threshold;

        cnt = cnt + 1;

    end
end

% find the indices for young, middle-aged and old animals
yo = find(ages <= 365);
ma = find(ages > 365 & ages <= 3*365);
ol = find(ages > 3*365);

% find the indices for male and female animals
F = find(sex == 0);
M = find(sex == 1);

% get the best estimate of spontaneous rate
srs = nan(1,largeNumber);
for i = 1:largeNumber
    if ~isnan(srsT(1,i))
        srs(1,i) = srsT(1,i);
    elseif ~isnan(srsT(2,i))
        srs(1,i) = srsT(2,i);
    elseif ~isnan(srsT(3,i))
        srs(1,i) = srsT(3,i);
    elseif ~isnan(srsT(4,i))
        srs(1,i) = srsT(4,i);
    end
end


%% settings for the plots

fs = 22; % fontsize
yo_sym = 'o'; yo_col = 'b'; % symbol and color to plot data of young-adult gerbils
ma_sym = '^'; ma_col = [251 186 0]./256; % symbol and color to plot data of middle-aged gerbils
ol_sym = 's'; ol_col = 'r'; % symbol and color to plot data of quiet_aged gerbils
sz = 70; % markersize

%% BF vs threshold - Fig 6A

% set limitations and ticks of the two axes
xl = [0.3 20]; 
yl = [0 100];
xt = [0.2 0.5 1 2 5 10 20];

% make running average
logvec = logspace(log10(xl(1)),log10(xl(2)),25); 
YLra = nan(3,length(logvec)-1); 
BFra = nan(1,length(logvec)-1);
for Age = 1:3 % loop through the age groups
    if Age == 1; ag = yo;
    elseif Age == 2; ag = ma;
    else; ag = ol;
    end
    for i = 1:length(logvec)-1 % loop through the logarithmically spaced vector
        idx = find(bfs(ag)/1000>logvec(i) & bfs(ag)/1000<=logvec(i+1));
        data = thr(ag);
        YLra(Age,i) = nanmean(data(idx));
        BFra(1,i) = logspace(log10(logvec(i)),log10(logvec(i+1)),1);
    end
end


figure;
scatter(bfs(yo)/1000,thr(yo),sz,yo_col,yo_sym)
hold on
scatter(bfs(ma)/1000,thr(ma),sz,ma_col,ma_sym)
scatter(bfs(ol)/1000,thr(ol),sz,ol_col,ol_sym)
plot(BFra,YLra(1,:),yo_col,'LineWidth',2)
plot(BFra,YLra(2,:),'color',ma_col,'LineWidth',2)
plot(BFra,YLra(3,:),ol_col,'LineWidth',2)
xlabel('Best Frequency (kHz)')
ylabel('Threshold (dB SPL)')
set(gca,'FontSize',fs,'XScale','log','XTick',xt)
xlim(xl); ylim(yl);
legend('Young adult','Middle aged','Old','Location','northwest')

%% BF distribution - Fig 6B

edges = logspace(log10(200),log10(20000),25);
xl = [300 20000];
xt = [200 500 1000 2000 5000 10000 20000];

figure; 
histogram(bfs(yo),edges,'FaceColor','k','EdgeColor','w')
set(gca,'XScale','log','FontSize',fs,'XTick',xt,'XTickLabel',{'0.2','0.5','1','2','5','10','20'})
xlim(xl)
ylabel('Fiber count')
xlabel('Best Frequency (kHz)')
box off

%% SR distribution

ag = yo;

edges = 0:5:200;
bfl = find(bfs<=3500);
bfh = find(bfs>3500);
ylo = intersect(ag,bfl);
yhi = intersect(ag,bfh);
xl = [0 200];
yl = [0 20];

figure
subplot(1,2,1)
histogram(srs(ylo),edges,'FaceColor','k','EdgeColor','w','Normalization','percentage')
set(gca,'FontSize',fs) %,'XTick',xt,'XTickLabel',{'0.2','0.5','1','2','5','10','20'})
xlim(xl)
ylim(yl)
ylabel('Distribution (%)')
xlabel('SR (spikes/s)')
title('BF < 3.5 kHz')
box off

subplot(1,2,2)
histogram(srs(yhi),edges,'FaceColor','k','EdgeColor','w','Normalization','percentage')
set(gca,'FontSize',fs) %,'XTick',xt,'XTickLabel',{'0.2','0.5','1','2','5','10','20'})
xlim(xl)
ylim(yl)
ylabel('Distribution (%)')
xlabel('SR (spikes/s)')
title('BF > 3.5 kHz')
box off

%% SR vs threshold

figure;
scatter(srs(yo),thr(yo),sz,'k',yo_sym)

%% SR vs bf - Fig 6C

xl = [0.3 20];
yl = [0 200];
xt = [0.2 0.5 1 2 5 10 20];

figure
plot([3.5 3.5],yl,'--k','LineWidth',1)
hold on
patch([3.5 xl(2) xl(2) 3.5],[yl(1) yl(1) yl(2) yl(2)],[0 0 1],'FaceAlpha',0.2,'EdgeColor','none')
scatter(bfs(yo)/1000,srs(yo),sz,'k',yo_sym)
set(gca,'FontSize',fs,'XScale','log','XTick',xt)
xlim(xl)
ylabel('Spontaneous rate (spikes/s)')
xlabel('Best Frequency (kHz)')

%% BF vs phase locking - Fig 6D

xl = [0.3 6]; 
yl = [0 1];
xt = [0.2 0.5 1 2 5 10 20];

% differentiate between low- and high-SR fibers, find the indices
srl = find(srs<18);
srh = find(srs>=18);
ylo = intersect(yo,srl);
yhi = intersect(yo,srh);
mlo = intersect(ma,srl);
mhi = intersect(ma,srh);
olo = intersect(ol,srl);
ohi = intersect(ol,srh);

figure
scatter(bfs(yhi)/1000,vss(yhi),sz,yo_col,yo_sym,'filled')
hold on
scatter(bfs(ylo)/1000,vss(ylo),sz,yo_col,yo_sym,'LineWidth',1)
scatter(bfs(mhi)/1000,vss(mhi),sz,ma_col,ma_sym,'filled')
scatter(bfs(mlo)/1000,vss(mlo),sz,ma_col,ma_sym,'LineWidth',1)
scatter(bfs(ohi)/1000,vss(ohi),sz,ol_col,ol_sym,'filled')
scatter(bfs(olo)/1000,vss(olo),sz,ol_col,ol_sym,'LineWidth',1)
set(gca,'FontSize',fs,'XScale','log','XTick',xt)
xlim(xl)
ylim(yl)
ylabel('Vector strength (max)')
xlabel('Best Frequency (kHz)')
% legend('Young adult','Middle aged','Old','Location','southwest')

% legend
x1 = 2.7; x2 = 3.4; x3 = 4;
y1 = 0.92; y2 = 0.86; y3 = 0.8;
fs2 = fs - 6;
scatter(x1,y1,sz,yo_col,yo_sym,'filled')
scatter(x2,y1,sz,yo_col,yo_sym,'LineWidth',1)
scatter(x1,y2,sz,ma_col,ma_sym,'filled')
scatter(x2,y2,sz,ma_col,ma_sym,'LineWidth',1)
scatter(x1,y3,sz,ol_col,ol_sym,'filled')
scatter(x2,y3,sz,ol_col,ol_sym,'LineWidth',1)
text(x3,y1,'Young adult','FontSize',fs2)
text(x3,y2,'Middle aged','FontSize',fs2)
text(x3,y3,'Old','FontSize',fs2)
text(x1,y1+0.04,'high SR','Rotation',45,'FontSize',fs2)
text(x2,y1+0.04,'low SR','Rotation',45,'FontSize',fs2)

%% BF vs click latency

% set limitations and ticks of the two axes
xl = [0.3 20];
yl = [0 3.5];
xt = [0.2 0.5 1 2 5 10 20];

% select the type of click latency analysis to be plotted
type = 3; % 1 = mean fsl, 2 = median fsl, 3 = 2 bins method, 4 = poisson method

% differentiate between low- and high-SR fibers, find the indices
srl = find(srs<18);
srh = find(srs>=18);
ylo = intersect(yo,srl);
yhi = intersect(yo,srh);
mlo = intersect(ma,srl);
mhi = intersect(ma,srh);
olo = intersect(ol,srl);
ohi = intersect(ol,srh);

% plot the latencies
figure;
scatter(bfs(yhi)/1000,clk(type,yhi),sz,yo_col,yo_sym,'filled')
hold on
scatter(bfs(ylo)/1000,clk(type,ylo),sz,yo_col,yo_sym,'LineWidth',1)
scatter(bfs(mhi)/1000,clk(type,mhi),sz,ma_col,ma_sym,'filled')
scatter(bfs(mlo)/1000,clk(type,mlo),sz,ma_col,ma_sym,'LineWidth',1)
scatter(bfs(ohi)/1000,clk(type,ohi),sz,ol_col,ol_sym,'filled')
scatter(bfs(olo)/1000,clk(type,olo),sz,ol_col,ol_sym,'LineWidth',1)
xlabel('Best Frequency (kHz)')
ylabel('Click latency (ms)')
set(gca,'FontSize',fs,'XScale','log','XTick',xt)
xlim(xl); ylim(yl);

% legend
x1 = 5; x2 = 6.5;
y1 = 3.2; y2 = 2.95; y3 = 2.7;
fs2 = fs - 6;
scatter(x1,y1,sz,yo_col,yo_sym,'filled')
scatter(x2,y1,sz,yo_col,yo_sym,'LineWidth',1)
scatter(x1,y2,sz,ma_col,ma_sym,'filled')
scatter(x2,y2,sz,ma_col,ma_sym,'LineWidth',1)
scatter(x1,y3,sz,ol_col,ol_sym,'filled')
scatter(x2,y3,sz,ol_col,ol_sym,'LineWidth',1)
text(8,y1,'Young adult','FontSize',fs2)
text(8,y2,'Middle aged','FontSize',fs2)
text(8,y3,'Old','FontSize',fs2)
text(x1,y1+0.15,'high SR','Rotation',45,'FontSize',fs2)
text(x2,y1+0.15,'low SR','Rotation',45,'FontSize',fs2)
