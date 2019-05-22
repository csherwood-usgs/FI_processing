function ws = stokes( rcm, rhof, rhos, mu )
% STOKES - calculate Stokes settling velocity (cgs units)
% ws = stokes( rcm, rhof, rhos, mu )
%
% Input: 
%   rcm - radius (NOT DIAMETER) in cm 
%   rhof - fluid density g/cm^3
%   rhos - sediment density g/cm^3
%   mu - dynamic fluid viscosity in poise
% Returns:
%   ws - settling velocity in cm/s
%
% 
% Rich Signell, USGS

g = 980.665;

ws = 2*g*rcm.^2.*(rhos-rhof)./(9.*mu);














