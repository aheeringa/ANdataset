function [thr,q10,indx]=getQ10dB_withCF(tc,freqs,CF)
%% 
% The same function as getQ10dB, except that with this function the CF is
% already known and can be entered in the function.
% 
% input:
%   tc: the tuning curve, a vector of firing rate thresholds at each 
%       stimulus frequency 
%   freqs: the vector with corresponding stimulus frequencies
%   CF: the known characteristic frequency. 
% output:
%   thr: the threshold at CF
%   q10: the Q10dB value, calculated by the CF/bandwidth of the tuning
%       curve at 10 dB above threshold
%   indx: the indices in the frequency vector between which the bandwidth
%       is situated, can be used for plotting the bandwidth 

% By Amarins Heeringa

%%
idx = find(freqs==CF);
thr = tc(idx);
db10=thr+10;

ind1=idx;
while tc(ind1)<db10 && ind1>1
    ind1=ind1-1;
end
ind2=idx;
while tc(ind2)<db10 && ind2<length(tc)
    ind2=ind2+1;
end
q10=CF/(freqs(ind2)-freqs(ind1));
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
    