function Urms = soulsby_ub( Hs, Tz, h )
% Orbital velocity using method of Soulsby and Smallman (1986)
% [Urms_Tn_over_Hs,t] = soulsby_ub( Tz, h )
%
% As presented on page 6 of
% Soulsby, R.L., 2006. "Simplified calculation of wave orbital velocities"
% HR Wallingford TR 155 Release 1.0.

% csherwood@usgs.gov
% 27-Feb-2013
g = 9.81;
Tn = ( h./g ).^(1/2); % Eqn 8
t = Tn./Tz;           % Eqn 27
A = (6500+(0.56 + 15.54.*t).^6).^(1./6.);   % Eqn 26
Urms_Tn_over_Hs = 0.25 ./((1.+A.*t.^2).^3); % Eqn 25
Urms = Urms_Tn_over_Hs*Hs/Tn
return