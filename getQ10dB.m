function [thr,CF,q10,indx]=getQ10dB(tc,freqs)
%% 
% input:
%   tc: the tuning curve, a vector of firing rate thresholds at each 
%       stimulus frequency 
%   freqs: the vector with corresponding stimulus frequencies
% output:
%   thr: the threshold at CF
%   CF: the estimated characteristic frequency in Hz, based on the tuning 
%       curve, that is visually checked, i.e. the frequency with the lowest
%       threshold in the tuning curve
%   q10: the Q10dB value, calculated by the CF/bandwidth of the tuning
%       curve at 10 dB above threshold
%   indx: the indices in the frequency vector between which the bandwidth
%       is situated, can be used for plotting the bandwidth 

% By Amarins Heeringa

%%
thr=min(tc);
db10=thr+10;
CFs=freqs(tc==thr);
CF=mean(CFs);

ind=find(tc==thr);

if length(ind)==1
    ind1=ind;
    while tc(ind1)<db10 && ind1>1
        ind1=ind1-1;
    end
    ind2=ind;
    while tc(ind2)<db10 && ind2<length(tc)
        ind2=ind2+1;
    end
    q10=CF/(freqs(ind2)-freqs(ind1));
elseif length(ind)>1 && length(ind)<8
    ind1=ind(1);
    while tc(ind1)<db10 && ind1>1
        ind1=ind1-1;
    end
    ind2=ind(end);
    while tc(ind2)<db10 && ind2<length(tc)
        ind2=ind2+1;
    end
    q10=CF/(freqs(ind2)-freqs(ind1));
else
    disp('not possible to determine q10')
    ind1=NaN; ind2=NaN; q10=NaN;    
end
indx=[ind1 ind2];



% determine out of range
if ind1==1 && tc(ind1)<db10; 
    oor=1;
    q10=NaN;
elseif ind2==length(tc) && tc(ind2)<db10; 
    oor=1;
    q10=NaN;
elseif isnan(ind1)
    oor=1;
    q10=NaN;
else
    oor=0;
end
    