function svout = soulsby_vanrijn( svin )
% soulsby_vanrijn - Compute total load (wave+current) transport
% svout = soulsby_vanrijn( svin )
%
% Input: svin
% svin.d50 = .25e-3 % median grain diameter (m)
% svin.d90 = .5e-3  % 90th percentile grain diameter (m)
% svin.urms = .26 % rms wave-orbital velocity (m/s)
% svin.U = .6  % depth-mean velocity magnitude(m/s)
% svin.tanB = 1/100 % tan(bed slope) in streamwise direction, positive if current is going uphill
% svin.rhos = 2650 % sediment density (kg/m3)
% svin.rhow =  999 % water density (kg/m3)
% svin.nu = 1.14e-6;     % kinematic viscosity of water (m2/s)
% svin.h = 5.            % water depth (m)
% svin.zo = 0.006 % bed roughness length (m)
%
% Output: svout
% svout.Cd = Cd        % drag coefficient for currents along
% svout.Dstar = Dstar  % dimensionless grain size
% svout.qt = qt        % volume transport rate (m2/s)
% svout.Qt = Qt        % mass transport rate (kg m-1 s-1
%
% Calculates total load according Soulsby (1997) Section 10.4
% This is appropriate for rippled bed condtions with Urms < U
%
% csherwood@usgs.gov

% constants
g = 9.81;        % gravitational acceleration (ms2)
vk = 0.4;        % von Karmans const.

% unpack input structure
d50 = svin.d50;
d90 = svin.d90;
U = svin.U;
urms = svin.urms;
tanB = svin.tanB;
rhos = svin.rhos;
rhow = svin.rhow;
nu = svin.nu;
h = svin.h;
zo = svin.zo;

s = rhos/rhow;
Cd = (vk/(log(h/zo)-1.)).^2;
Dstar = d50*((g*(s-1))/nu^2)^(1./3.);

% van Rijn current threshold of motion (Soulsby, 1997, Eqn. 71)
if(d50>=100e-6 && d50<=500e-6)
    Ucr = 0.19*d50^.1 * log10(4*h/d90); % eqn 71a
elseif( d50>500e-6 && d50 <= 2000.e-6) 
    Ucr = 8.5*d50^0.6 * log10(4*h/d90); % eqn 71b
else
    error('d50 out of range - consider adding Soulsby eqns. 72ab')
end
% Equations 136b - d
denom = ((s-1)*g*d50)^1.2;
Asb = (0.005*h*(d50/h)^1.2)/denom;
Ass = (0.012*d50*Dstar^(-0.6))/denom;
As = Asb+Ass;
% Eqn. 136a
qt = 0.;
if( sqrt(U^2 + (0.018*urms^2)/Cd) - Ucr  )>0.
   qt = ( As*U*( sqrt( U^2 + (0.018*urms^2)/Cd)-Ucr )^2.4 )*(1.-1.6*tanB);
end
Qt = qt*rhos;

% pack output structure
svout.Cd = Cd;
svout.Dstar = Dstar;
svout.qt = qt;
svout.Qt = Qt;