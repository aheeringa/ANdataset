function [ vs, prob, ntot, xbin, nspikes ] = TytoView_calcVS(sptimes, freq, nbin)
%------------------------------------------------------------------------
% TytoView_calcVS.m
%------------------------------------------------------------------------
%  calculating VS from spike timing data and 
%  calculating phase histogram data 
%------------------------------------------------------------------------
% Input:
%   sptimes         [1xN] array for spike timings (in ms)
%   freq            frequency (Hz)
%   nbin            number of bins per cycle (for phase histograms)
% 
% Output:
%   vs              vector strength
%   prob            significance probability (Rayleigh test)
%   ntot            total bumber of spikes counted
%   xbin            [1xnbin] array for histogram (in cycle)
%   nspikes         [1xnbin] spike count for each bin 
%   ---- Note: bar(xbin,nspikes,1) gives a phase histogram
%------------------------------------------------------------------------
%  Go Ashida 
%   ashida@umd.edu
%------------------------------------------------------------------------
% Created: 18 March, 2012 by GA
%
% Revisions: 
% 
%------------------------------------------------------------------------

% calculate VS
ntot = length(sptimes);
spphase = 2 * pi * freq / 1000 * sptimes; 
c = sum( cos(spphase) ) / ntot;
s = sum( sin(spphase) ) / ntot;
vs = sqrt( c*c + s*s ); 

% significance probability
% P = exp(-N r^2) gives a good approximation for N>50
ptmp = exp(-ntot * vs * vs);

if ntot > 50
    prob = ptmp;
else 
    prob = [];
end

% calculate spike counts in each phase bin
wbin = 1.0/nbin;
xbin = ((1:nbin)-0.5) / nbin; 
modphase = mod(sptimes*freq/1000, 1); % phase in cycles [0,1]
nspikes = zeros(1,nbin);
for i = 1:nbin
    t0 =  (i-1) * wbin; 
    t1 =    i   * wbin; 
    nspikes(i) = sum( (modphase>=t0) & (modphase<t1) );
end


