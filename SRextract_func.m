function [sr] = SRextract_func(curvedata, curvesettings, stim)

% Function to analyze the spontaneous rate from a designated spontaneous
% rate recording obtained by the EXTSTIM function of Tytology.

% input:
%   curvedata: a curvedata struct, containing the spike times
%   curvesettings: a curvesettings struct, containing the recording metadata
%   stim: the presented acoustic stimulus, an empty wav-file, is used to
%   determine the length of the recording.
% output:
%   sr: the spontaneous rate in spikes/s

% October, 2019
% By: Amarins Heeringa, AG Koeppl

nwSpikCnt=[]; 
for i=1:length(curvedata.spike_times)
    nwSpikCnt(i)=length(curvedata.spike_times{1,i});
end

T = length(stim)/curvesettings.Fs(1); 
sr = sum(nwSpikCnt)/(T*curvesettings.extstim.Reps);
