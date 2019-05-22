function [theta,s02]=circvar(u,v)
% [theta,s02]=circvar(u,v)
% Calculate circular variance
% Davis, 1986, pp 314 - 321
%
% Rbar is the length of the sum of the unit vectors (ranges from 0 to 1)
% s02 is 1-Rbar
u = u(:);
v = v(:);
n = length(u);
[ss,dd]=pcoord(u,v);
[ur,vr]=xycoord(ones(size(dd)),dd);
Cbar = sum(ur)/n;
Sbar = sum(vr)/n;
sqrt(Cbar.^2+Sbar.^2);
[R,theta]=pcoord(sum(ur),sum(vr));
Rbar = R/n;
s02 = 1-Rbar;           % Eqn. 5.44 in Davis   