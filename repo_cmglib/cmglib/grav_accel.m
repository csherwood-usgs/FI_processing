function g0 = grav_accel( lat )
% GRAV_ACCEL - IAG Geodetic Reference System 1980 Intl. Gravity Formula
% g0 = grav_accel( lat )
%
% Input: latitude (degrees)
% Returns: gravitational acceleration (m/s)

% http://geophysics.ou.edu/solid_earth/notes/potential/igf.htm
% csherwood@usgs.gov
% last revised August 5, 2006
rlat = pi*abs(lat)/180;
g0 = 9.7803267714*( (1+0.00193185138639*sin(rlat).^2)./...
   sqrt((1-0.00669437999013*sin(rlat).^2)) );