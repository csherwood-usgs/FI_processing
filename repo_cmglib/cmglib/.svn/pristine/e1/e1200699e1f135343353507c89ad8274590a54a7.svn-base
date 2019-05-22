function tcrit = taucrit(rhos,rho,rnu,Dcm)
% TAUCRIT - Calculate critical shear stress for erosion - CGS units
% tcrit = taucrit(rhos,rho,rnu,Dcm)
%
%  Smith's tau crit curve (no shape correction)
%  Fit of Pat Wiberg
g=980.665;
delrho=(rhos-rho)./rho;
zeta=delrho*g.*Dcm.^3 ./rnu.^2;
zetal=log10(zeta);
if (zetal<=3.6),
  tst=.00056*zetal.^3+.00115*zetal.*zetal-0.0316*zetal+.10455;
else
  tst=-.00070*zetal.^3+.01204*zetal.*zetal-0.0584*zetal+.11918;
end
tcrit=tst .* (rhos-rho)*g.*Dcm;


