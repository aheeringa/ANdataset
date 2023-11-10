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
ages = nan(1,length(all_exp));
sex = nan(1,length(all_exp));
abr = nan(1,length(all_exp));
weights = nan(1,length(all_exp));
nunit = nan(1,length(all_exp));
bfsR = nan(2,length(all_exp));
oxy = nan(1,length(all_exp));
noise = nan(1,length(all_exp));
sps = nan(1,length(all_exp));
cvc = nan(1,length(all_exp));

cnt = 1;

% loop through the animals and store the animal specific metadata
for i = 1:length(all_exp)

    ages(i) = all_exp(i).info.age;
    if strcmp(all_exp(i).info.sex,'F')
        sex(i) = 0;
    elseif strcmp(all_exp(i).info.sex,'M')
        sex(i) = 1;
    end
    abr(i) = all_exp(i).info.ABR.threshold;
    weights(i) = all_exp(i).info.weight;
    nunit(i) = length(all_exp(i).data);
    oxy(i) = all_exp(i).info.anesthesia.oxygen;
    
    bfs = []; noiT = []; spsT = []; cvcT = [];
    for j = 1:length(all_exp(i).data)
        
        if ~isempty(all_exp(i).data(j).BF)
            bfs = [bfs all_exp(i).data(j).BF.analysis.bf];
        elseif ~isempty(all_exp(i).data(j).CF)
            bfs = [bfs all_exp(i).data(j).CF.analysis.cf];
        end
        if ~isempty(all_exp(i).data(j).NOISE)
            noiT = [noiT j];
        end
        if ~isempty(all_exp(i).data(j).SPS)
            spsT = [spsT j];
        end
        if ~isempty(all_exp(i).data(j).CVC)
            cvcT = [cvcT j];
        end
    end
    bfsR(1,i) = min(bfs);
    bfsR(2,i) = max(bfs);
    noise(i) = length(noiT);
    sps(i) = length(spsT);
    cvc(i) = length(cvcT);
end

% find the indices for young, middle-aged and old animals
yo = find(ages <= 365);
ma = find(ages > 365 & ages <= 3*365);
ol = find(ages > 3*365);

% find the indices for male and female animals
F = find(sex == 0);
M = find(sex == 1);

%% Input for Table 1

% number of animals (females)
fprintf('Total number of young-adult animals (females): %d (%d)\n',length(yo),length(yo)-sum(sex(yo)))
fprintf('Total number of middle-aged animals (females): %d (%d)\n',length(ma),length(ma)-sum(sex(ma)))
fprintf('Total number of quiet-aged animals (females): %d (%d)\n',length(ol),length(ol)-sum(sex(ol)))

% age in months
mo = 365/12;
fprintf('Mean +/- std (range) of age in months, young-adults: %.1f +/- %.1f (%.1f - %.1f) mo\n',nanmean(ages(yo))/mo,nanstd(ages(yo))/mo,min(ages(yo))/mo,max(ages(yo))/mo)
fprintf('Mean +/- std (range) of age in months, middle-aged: %.1f +/- %.1f (%.1f - %.1f) mo\n',nanmean(ages(ma))/mo,nanstd(ages(ma))/mo,min(ages(ma))/mo,max(ages(ma))/mo)
fprintf('Mean +/- std (range) of age in months, quiet-aged: %.1f +/- %.1f (%.1f - %.1f) mo\n',nanmean(ages(ol))/mo,nanstd(ages(ol))/mo,min(ages(ol))/mo,max(ages(ol))/mo)

% abr threshold
fprintf('Mean +/- std (range) of abr threshold, young-adults: %.1f +/- %.1f (%.1f - %.1f) mo\n',nanmean(abr(yo)),nanstd(abr(yo)),min(abr(yo)),max(abr(yo)))
fprintf('Mean +/- std (range) of abr threshold, middle-aged: %.1f +/- %.1f (%.1f - %.1f) mo\n',nanmean(abr(ma)),nanstd(abr(ma)),min(abr(ma)),max(abr(ma)))
fprintf('Mean +/- std (range) of abr threshold, quiet-aged: %.1f +/- %.1f (%.1f - %.1f) mo\n',nanmean(abr(ol)),nanstd(abr(ol)),min(abr(ol)),max(abr(ol)))

% weights
fprintf('Mean +/- std (range) of weight, young-adults: %.1f +/- %.1f (%.1f - %.1f) mo\n',nanmean(weights(yo)),nanstd(weights(yo)),min(weights(yo)),max(weights(yo)))
fprintf('Mean +/- std (range) of weight, middle-aged: %.1f +/- %.1f (%.1f - %.1f) mo\n',nanmean(weights(ma)),nanstd(weights(ma)),min(weights(ma)),max(weights(ma)))
fprintf('Mean +/- std (range) of weight, quiet-aged: %.1f +/- %.1f (%.1f - %.1f) mo\n',nanmean(weights(ol)),nanstd(weights(ol)),min(weights(ol)),max(weights(ol)))

% total number of single units
fprintf('Total number of units, young-adults: %d\n',sum(nunit(yo)))
fprintf('Total number of units, middle-aged: %d\n',sum(nunit(ma)))
fprintf('Total number of units, quiet-aged: %d\n',sum(nunit(ol)))

% range of BFs
fprintf('Young-adults, min BF: %d, max BF: %d\n',min(bfsR(1,yo)),max(bfsR(2,yo)))
fprintf('Middle-aged, min BF: %d, max BF: %d\n',min(bfsR(1,ma)),max(bfsR(2,ma)))
fprintf('Quiet-aged, min BF: %d, max BF: %d\n',min(bfsR(1,ol)),max(bfsR(2,ol)))

%% number of experiments with complex stimuli

fprintf('In %d out of %d experiments, noise responses were collected\n',length(find(noise)),length(noise))
fprintf('In %d out of %d experiments, sps responses were collected\n',length(find(sps)),length(sps))
fprintf('In %d out of %d experiments, cvc responses were collected\n',length(find(cvc)),length(cvc))

%% make the metadata sheet

T = table;

for i = 1:length(all_exp)
    
    % remove the '_1' from the structs that need it
    if contains(all_exp(i).animalID,'_1')
        all_exp(i).animalID = all_exp(i).animalID(1:end-2);
    end
    T.animalID{i} = all_exp(i).animalID;
    T.age_days(i) = all_exp(i).info.age;
    T.sex{i} = all_exp(i).info.sex;
    T.ABR_threshold(i) = all_exp(i).info.ABR.threshold;
    T.n_units(i) = length(all_exp(i).data);
    nbf = 0; bfs = []; ncf = 0; nph = 0; nclick = 0; nrlf = 0; nsr = 0; nnoise = 0; nsps = 0; ncvc = 0; ntfs1 = 0; nvcv = 0;
    for j = 1:length(all_exp(i).data)
        if ~isempty(all_exp(i).data(j).BF)
            nbf = nbf + 1;
            bfs = [bfs all_exp(i).data(j).BF.analysis.bf];
        end
        if ~isempty(all_exp(i).data(j).CF)
            ncf = ncf + 1;
        end
        if ~isempty(all_exp(i).data(j).PH)
            nph = nph + 1;
        end
        if ~isempty(all_exp(i).data(j).CLICK)
            nclick = nclick + 1;
        end
        if ~isempty(all_exp(i).data(j).RLF)
            nrlf = nrlf + 1;
        end
        if ~isempty(all_exp(i).data(j).SR)
            nsr = nsr + 1;
        end
        if ~isempty(all_exp(i).data(j).NOISE)
            nnoise = nnoise + 1;
        end
        if ~isempty(all_exp(i).data(j).SPS)
            nsps = nsps + 1;
        end
        if ~isempty(all_exp(i).data(j).CVC)
            ncvc = ncvc + 1;
        end
        if ~isempty(all_exp(i).data(j).TFS1)
            ntfs1 = ntfs1 + 1;
        end
        if ~isempty(all_exp(i).data(j).VCV)
            nvcv = nvcv + 1;
        end
    end
    T.BF_range{i} = [num2str(min(bfs)) ' - ' num2str(max(bfs))];
    T.n_BF(i) = nbf;
    T.n_CF(i) = ncf;
    T.n_PH(i) = nph;
    T.n_CLICK(i) = nclick;
    T.n_RLF(i) = nrlf;
    T.n_SR(i) = nsr;
    T.n_NOISE(i) = nnoise;
    T.n_SPS(i) = nsps;
    T.n_CVC(i) = ncvc;
    T.n_TFS1(i) = ntfs1;
    T.n_VCV(i) = nvcv;
end

writetable(T,'metadata_sheet.csv','Delimiter',',')

%% plot the animal's age against it's ABR threshold

mo = 365/12;
sz = 70;
f_col = [11 154 51]./256; % green
m_col = [111 32 130]./256; % purple
fs = 25;
xl = [0 45];

figure; 
scatter(ages(F)./mo, abr(F), sz, f_col, 'filled')
hold on
scatter(ages(M)./mo, abr(M), sz, m_col, 'filled')
xlabel('Age in months')
ylabel('ABR threshold in dB SPL')
set(gca,'FontSize',fs,'XTick',[0 6 12 18 24 30 36 42],'XTickLabel',{'0','','12','','24','','36',''})
legend('Female','Male', 'location','northwest')
xlim(xl)
