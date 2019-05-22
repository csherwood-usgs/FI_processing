function g=gravity(lat,z)
% Gravity - Calculate g as function of latitude and elevation
%  g=gravity(lat,z)
%
% Input:
%   lat - Latitude in degrees
%   z - elevation wrt to the geoid (positive up; meters)
% Returned:
%   g - acceleration due to gravity (m s-2)
%
% "Atmosphere-Ocean Dynamics, A. E. Gill, 1982, Appendix Two, p. 597

% csherwood@usgs.gov
% 7 June 2006
a = 6371000; % radius of earth in m
rlat = (pi/180).*abs(lat);
g = (9.78032+0.005172*sin(rlat).^2-0.00006*sin(2*rlat).^2) .* (1+z/a).^-2;