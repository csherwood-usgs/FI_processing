function ws = dset( csf, p, Dcm, rhos, rho, nu )
% DSET  Calculate grain settling velocity according to Dietrich (1982)
%
% ws = dset( csf, p, Dcm, rhos, rho, nu )
%
% use cgs units
% ws = settling velocity (cm s-1)
% Dcm = diameter (cm)
% rhos = density solids (g cm-3)
% rho = density water (g cm-3)
% nu = kinematic viscosity (stoke, cm2 s-1)
% **note stoke = poise/fluid density
%
% Dietrich, W.E. (1982) Water Res. Research 18(6):1615-1626
% Original program by Pat Wiberg
% Reasonable (well-rounded) Corey shape factor (csf)
% value is 0.8, range is 0.2 - 1.0
% Typical Powers roundness value is 5, range is 2.0 - 6.0
% Note: number returned is negative, in keeping w/ pos. upward coords
% Only Dcm may be a vector, and rhos must be constant
%
%
% Chris Sherwood, USGS
% March 17, 1999

g = 980.665;
delrho = (rhos-rho)/rho;
zeta = delrho*g*Dcm.*Dcm.*Dcm ./(nu*nu);
zetal = log10( zeta );
zetal2 = zetal.*zetal;
zetal3 = zetal.*zetal2;
zetal4 = zetal2.*zetal2;

%/* eqn 9 */
a = -3.76715 + 1.92944*zetal -0.09815*zetal2 -0.00575*zetal2 +0.00056*zetal4;
omc=1.-csf;

%/* eqn 16 */
b = log10(1.-omc/.85)- omc^2.3 * tanh( zetal-4.6 )+0.3*(0.5-csf)*...
    omc*omc*( zetal-4.6 );
%/* eqn 18 */
c = (.65-(csf/2.83)*tanh(zetal-4.6)).^(1.+(3.5-p)/2.5);
wstr = c .* 10 .^(a+b);
ws = -((delrho*g*nu*wstr).^(1./3.));











