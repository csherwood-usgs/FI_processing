function ub = ub_f(T,k,H,h)
% UB_F  Calculate near-bottom wave-orbital velocity amplitude
%
% function ub = ub_f(T,k,H,h)
%
% Input:
%   T  = sig. wave period
%   H  = sig. wave height ( = 2 * amplitude of surface wave)
%   h  = water depth
%   k  = wave number (use waven or qkhf)
%  Returns:
%   ub_f = Hs*pi / T*sinh(kh)

% Chris Sherwood, USGS
% March 17, 1999
%
% cf. Dyer(1986) eqn. 3.50, p. 98
% or  Komar(1976) p. 42

twopi = 2*pi;
w= twopi/T;
kh = k*h;
amp = H / (2.*sinh(kh));
ub = w * amp;



