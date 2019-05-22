function del = tbdel( Ua, Ub, lat )
% TBDEL - Thickness of tidal boundary layer
% del = tbdel( Ua, Ub, lat )
%
% Soulsby (1997) Eqn. 27
% Input:
%   Ua - alongshore velocity (m/s)
%   Ub - crossshore velocity (m/s)
%   lat - latitude
f = 1.4544e-4*sin( lat*pi/180 )
sig = 1.4052e-4
del = 0.0038*(Ua*sig-Ub*f)/(sig.^2-f.^2)