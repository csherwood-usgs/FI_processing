function [Hm,Tm,t]=fetchlim_deep(Ua,F)
% FETCHLIM_DEEP - Wave height, period, and evolution time, deepwater
% function [Hm,Tm,t]=fetchlim_deep(Ua,F)
%
% Returns fetch-limited wave height, period, and evolution time
% in deep water given wind speed (at 10 m) and fetch [MKS Units]
% [ref: 1984 Shore Protection Manual]

% Equation 3-33a
Hm = 5.112e-4 * Ua .* sqrt(F);
% Equation 3-34a
Tm = 6.238e-2 * (Ua .* F).^(1./3.);
% Equation 3-35a
t = 3.215e+1 * (F.*F./Ua).^(1./3.);
