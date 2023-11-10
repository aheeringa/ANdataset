function [rastDataFS, rastDataRest] = getRasterData3(curvedata,trials,totDelay,start,stop,SpikCnt)
%%
% input:
%   curvedata: a curvedata struct, containing the spike times
%   trials: a vector with all the trials that need to be included
%   totDelay: the delay between start of recording and start of the click
%   start, stop: indicate the window between which the first spikes can be
%       taken into account
%   SpikCnt: a vector with all the spikes counts in each trial
% output:
%   rastDataFS: a matrix with the first column being the trial number and
%       the second column being the first spike after totDelay + start. If the
%       first spike is > totDelay + stop, no first spike is stored (NaN).
%   rastDataRest: a matrix with the first column being the trial number and
%       the second column being all spikes in the recording. Plotting first
%       column on y-axis and second column on the x-axis will show the
%       rasterplot of the recording. 
% 
% By: Amarins Heeringa

%%

% pre-allocation
rastDataFS = nan(length(trials),2);
rastDataRest = nan(sum(SpikCnt(trials)),2);
cnt=1;

% loop through the trials
for sp=1:length(trials)
    tmpSpik=curvedata.spike_times{1,trials(sp)};
    sp1=find(tmpSpik>totDelay+start & tmpSpik<totDelay+stop);

    % store the first spike
    if ~isempty(sp1)
        rastDataFS(sp,1)=sp; % store the trial number
        rastDataFS(sp,2)=tmpSpik(sp1(1)); % store the spike time of the first spike
        % tmpSpik(sp1(1))=[];
    end

    % store all of the spikes
    rastDataRest(cnt:cnt+length(tmpSpik)-1,1)=ones(1,length(tmpSpik)).*sp; % store the trial number
    rastDataRest(cnt:cnt+length(tmpSpik)-1,2)=tmpSpik; % store the spikes
    cnt=cnt+length(tmpSpik); % counter to know where to start in the next loop
end
