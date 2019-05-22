function wavetrace(Ho,T,ho,thetao,dhdx)
% WAVETRACE - Refraction and shoaling of linear waves
% wavetrace(Ho,T,ho,thetao,dhdx)
%
% Input:
%    Ho - Local wave height (m)
%    T  - Period (s)
%    ho - Local water depth (m)
%    thetao - Local angle of wave approach (degrees)
%    dhdx - Bottom slope everywhere
%
% Right now, nothing is returned
% Chris Sherwood, USGS
% Last revised 3/30/04

g = 9.81;
thetaro = pi*thetao/180; % radians
ko = waven(T,ho);
Lo = 2*pi/ko
Ldeep = g*(T.^2)/(2*pi)
Co = g*T/(2*pi)*tanh( 2*pi*ho ./ Lo )
no = 0.5*(1+ 2*ko*ho./(sinh(2*ko*ho)))
Cgo = no*Co
Co*sin(thetaro)


hmin = .5*Ho;
hmax = max(Lo, ho);

dx = 1000
% this builds a depth array that depends on spacing dx,
% centered about the input depth value
h_offshore = flipud( (ho : dhdx*dx :hmax)' );
h_onshore = (ho : -dhdx*dx : hmin )';
h = [h_offshore(1:end-1);h_onshore];
w = 2*pi/T;
% Explicit wavelength approximation from Fenton and McKee (1989)
% as quoted by Fenton, The Sea, vol. 9, A, 1990 and Dalrymple's web page
% (or use qkhf or waven)
L = Ldeep*(tanh( (w.^2 .* h / g).^(3/4) )).^(2/3);
k = 2*pi ./L;
kh = k.*h;
n = 0.5*(1+ 2*kh./(sinh(2*kh)));
C = g*T/(2*pi)*tanh( 2*pi*h ./ L );
Cg = n.*C;
sinthetar = C*sin(thetaro)./Co;
thetar = asin(sinthetar);
theta = (180/pi)*thetar;

Kr = sqrt( cos(thetaro) ./cos(thetar) );
Ks = sqrt( Cgo ./ Cg );
H = Ho.*Ks.*Kr;

% Print depth, local wavelenth, local angle, refraction and shoaling
% coeffs, and local height
disp([h,L,theta,Kr,Ks,H])
