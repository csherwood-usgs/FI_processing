function rcm = ungibbs( ws, rhof, rhos, mu )
% UNGIBBS - Calculate sphere size from Gibbs settling velocity (CGS units)
%
% Input:
%   ws - settling velocity in cm/s
%   rhof - fluid density g/cm^3
%   rhos - sediment density g/cm^3
%   mu - dynamic fluid viscosity in poise
%
% Returns:
%   rcm - radius (not diameter) in cm
%
% Eqn. on p. 10 of Gibbs, Matthews, and Link,
% (1971). J. Sed. Petrol. 41:7-18.

% Chris Sherwood, USGS
% Last revised Sept. 4, 2001

g = 980.665;

rcm=(0.055804*ws*ws*rhof+...
     sqrt(0.003114*ws.^4*rhof^2+(g*(rhos-rhof))*(4.5*mu*ws+0.008705*ws*ws*rhof)))...
    /(g*(rhos-rhof));

