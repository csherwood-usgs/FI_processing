function ub = ubqf( T, Hs, h )
% UBQF - Quick approximation of near-bottom orbital velocity
% ub = ubqf( T, Hs, h )
%
% (eg, Dyer, 1986, Eqn 3.50, p 98).
%
% Uses iterative method of Soulsby (2006) for kh
% Calls qkhfs.m
%
% Input:
%  T Wave period [s]
%  Hs Wave height [m]
%  h Water depth [m]
% Returns:
%  Wave-orbital velocity amplitude ub = Hs*pi / T*sinh(kh) [m/s]

% Chris Sherwood, USGS
% March 17, 1999
% Modified April 27, 2005 for 2D matrices
% Modified Oct 5 to call qkhfs instead of qkhf
% Modified Oct 6 to remove error checks
% stderr = 1;
% % Idiot lights
% if( any(T < 0.0001) )
%   fprintf(stderr,'Warning from ubqf: small T = %g\n',T);
%   ub(T<0) = 0.;
%   return
% end
% if(any(Hs(:) >= 0.78.*h(:))),
%     fprintf(stderr,'Warning from ubqf: wave is breaking; Hs limited.\n');
%     Hs( Hs >= 0.78.*h )=0.78*Hs;
% end
w=2*pi./T;
kh=qkhfs(w,h);
amp = Hs ./ (2.*sinh(kh));
ub = w .* amp;




