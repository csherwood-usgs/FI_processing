function ws = ws_soulsby( d, rhos, rhow, nu, g )
% WS_SOULSBY - Sand settling velocity according to Soulsby (1997)
% ws = ws_soulsby( d ,[rhos, rhow, nu, g] )
%
% Input:
%   d - grain diameter (m)
%   rhos - sediment density (kg / m^3 ) [default = 2650]
%   rhow - water density (kg / m^3) [default = 1027]
%   nu = kinematic viscosity (m^2 /s ) [ default = 1.36e-6 ]
% Returned:
%   ws - settling velocity (m/s)
%
% Note: will work with CGS units if all five arguments are provided in CGS
%
% Eqn. SC(102) in "Dynamics of Marine Sands", Soulsby (1997)

% csherwood@usgs
% Last revised December 20, 2004

% if no arguments, assume MKS units
if(exist('g')~=1),
    g = 9.81;
end
% default values for T = 10, S = 35
if(exist('nu')~=1),
    nu = 1.36e-6;
end
if(exist('rhow')~=1),
    rhow = 1027;
end
if(exist('rhos')~=1),
    rhos = 2650;
end
s = rhos/rhow;
dstr = d .* ( g*(s-1)/(nu^2) ).^(1/3);
ws = (nu ./ d ) .* ( sqrt( 10.36^2 + 1.049 * dstr.^3) - 10.36 );

