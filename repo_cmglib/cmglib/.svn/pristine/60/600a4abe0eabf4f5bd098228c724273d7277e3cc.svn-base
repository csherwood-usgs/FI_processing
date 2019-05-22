function ws = gibbs( rcm, rhof, rhos, mu )
% GIBBS - calculate Gibbs settling velocity (cgs units)
% ws = gibbs( rcm, rhof, rhos, mu )
%
% Input: 
%   rcm - radius (NOT DIAMETER) in cm 
%   rhof - fluid density g/cm^3
%   rhos - sediment density g/cm^3
%   mu - dynamic fluid viscosity in poise
% Returns:
%   ws - settling velocity in cm/s
%
% Eqn. on p. 10 of Gibbs, Matthews, and Link,
% (1971). J. Sed. Petrol. 41:7-18.

% Chris Sherwood, USGS
% Last revised Sept 4, 2001 
g = 980.665;

ws=(-3*mu+sqrt(9*mu*mu+g*rcm.*rcm*rhof*(rhos-rhof).*(0.015476+0.19841*rcm)))...
     ./(rhof*(0.011607+0.14881*rcm));













