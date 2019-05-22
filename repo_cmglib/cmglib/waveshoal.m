function [Kr,Ks,H1,theta1] = waveshoal(h0,H0,T,theta0,h1)
% WAVESHOAL - Refraction and shoaling of linear waves
%             assuming straight and parallel depth contours
%             neglects wave dissipation, breaking, and wind growth !
%
% [Kr,Ks,H1,theta1] = waveshoal(h0,H0,T,theta0,h1)
%
% Input:
%    h0     - Water depth at location of known wave condition (m)
%    H0     - Wave height at depth h0 (m)
%    T      - Wave period (s)
%    theta0 - Angle of wave approach at depth h0
%             (degrees relative to shore-normal)
%             (should be in range -180 to 180)
%    h1     - Depth for output wave parameters (m)
%             (may be a vector of depths)
% Output:
%    Kr     - Refraction coefficient
%    Ks     - Shoaling coefficient
%    H1     - Wave height at depth h1 (m)
%    theta1 - Angle of wave approach at depth h1 (degrees)
%
% Chris Sherwood & Giles Lesser, USGS
% Last revised 6/1/07

g = 9.81;
thetar0 = pi*theta0/180; % radians
k0  = waven(T,h0);
for f=1:length(h1);
    k1(f) = waven(T,h1(f));
end
kh0 = k0.*h0;
kh1 = k1.*h1;
L0  = 2*pi/k0;
L1  = 2*pi./k1;
C0  = L0/T;
C1  = L1./T;
n0  = 0.5 + kh0./(sinh(2*kh0));
n1  = 0.5 + kh1./(sinh(2*kh1));
Cg0 = n0.*C0;
Cg1 = n1.*C1;

% wave direction
sinthetar1 = C1*sin(thetar0)./C0;
thetar1 = asin(sinthetar1);
theta1 = (180/pi)*thetar1;

% wave height
Kr = sqrt( cos(thetar0) ./cos(thetar1) );
Ks = sqrt( Cg0 ./ Cg1 );
H1 = H0.*Ks.*Kr;
