% This script shows examples on how to call the function that were used to
% analyze that data. Each of these functions requires the 'curvedata' and
% 'curvesettings' field. 

% By: Amarins Heeringa, November 6, 2023

%%
clear; close all; clc;

% load the required data struct, 'exp', from one animal
load('/Users/amarins/Documents/PostDoc 2 stuff/Oldenburg/Writings/2023 AN data/dataset/dataset 11.9.23/G220908.mat')

% indicate the unit idx (unit) or unit name
unit = 8;

% === optional code, when the unit_name is entered to find the unit index
% unit_name = '1p_449';
% for i = 1:length(exp.data)
%     if strcmp(exp.data(i).unit,unit_name)
%         unit = i;
%         break
%     end
% end


%% BF
if ~isempty(exp.data(unit).BF)
    curvedata = exp.data(unit).BF.curvedata;
    curvesettings = exp.data(unit).BF.curvesettings;
    [BF_bf, BF_sr, BF_bandwidth, BF_level] = BFextract_func(curvedata, curvesettings)
else
    disp('There is no BF recording for this unit')
end

%% CF
if ~isempty(exp.data(unit).CF)
    curvedata = exp.data(unit).CF.curvedata;
    curvesettings = exp.data(unit).CF.curvesettings;
    [CF_cf, CF_threshold, CF_Q10, CF_sr] = CFextract_func(curvedata, curvesettings)
else
    disp('There is no CF recording for this unit')
end

%% PH
if ~isempty(exp.data(unit).PH)
    curvedata = exp.data(unit).PH.curvedata;
    curvesettings = exp.data(unit).PH.curvesettings;
    [PH_freq, PH_levels, PH_vs, PH_prob, PH_vsVec, PH_pVec, PH_rates, PH_sr] ...
        = PHextract_func(curvedata, curvesettings)
else
    disp('There is no PH recording for this unit')
end

%% CLICK
if ~isempty(exp.data(unit).CLICK)
    curvedata = exp.data(unit).CLICK.curvedata;
    curvesettings = exp.data(unit).CLICK.curvesettings;
    [CLICK_lat_poisson, CLICK_lat_2bins, CLICK_fsl_mean, CLICK_fsl_median, ...
        CLICK_fsl_jit_std, CLICK_fsl_jit_var, CLICK_fsl_jit_iqr] = Clickextract_func(curvedata, curvesettings)
else
    disp('There is no CLICK recording for this unit')
end

%% RLF
if ~isempty(exp.data(unit).RLF)
    curvedata = exp.data(unit).RLF.curvedata;
    curvesettings = exp.data(unit).RLF.curvesettings;
    [RLF_threshold, RLF_sr, RLF_rates, RLF_stdevs, RLF_levels, RLF_freq] = ...
        RLFextract_func(curvedata, curvesettings)
else
    disp('There is no RLF recording for this unit')
end

%% SR
if ~isempty(exp.data(unit).SR)
    curvedata = exp.data(unit).SR.curvedata;
    curvesettings = exp.data(unit).SR.curvesettings;
    stim = exp.data(unit).SR.acoustic_stimulus.waveform;
    [SR_sr] = SRextract_func(curvedata, curvesettings, stim)
else
    disp('There is no SR recording for this unit')
end

