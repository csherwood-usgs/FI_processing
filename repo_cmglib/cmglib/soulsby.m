function m = soulsby( ubr, wr, ucr, zr, phiwc, kN, iverbose )
% SOULSBY - Wave-current model from Soulsby (1997)
% function m =  soulsby( ubr, wr, ucr, zr, phiwc, kN, iverbose )
%
%   Input:
%      ubr = rep. wave-orbital velocity amplitude outside wbl [m/s]
%      wr = rep. angular wave frequency = 2pi/T [rad/s]
%      ucr = current velocity at height zr [m/s]
%      zr = reference height for current velocity [m]
%      phiwc = angle between currents and waves at zr (radians)
%      kN = bottom roughness height (e.q. Nikuradse k) [m]
%      iverbose = switch; when 1, extra output
%   Returned in structure m:
%      m.ustrc  = current friction velocity         u*c [m/s]
%      m.ustw   = w-c combined friction velocity    u*w [m/s]
%      m.ustrwm = wave max. friction velocity      u*wm [m/s]
%      m.fw = wave friction factor [ ]
%      m.zoa = apparent bottom roughness [m]


% Chris Sherwood, USGS
% Last revised December, 2004
vk = 0.41;
zo = kN/30.;

if(ubr <= 0.01),
   if(ucr <= 0.01),
%     ...no waves or currents
      m.ustrc = 0;
      m.ustrw = 0;
      m.ustrwm = 0;
      m.fw  =   0;
      m.zoa =   zo;
      return
   end
%  ...no waves
   ustrc = ucr * vk / log(zr/zo); 
   m.ustrc = ustrc;
   m.ustrw = 0;
   m.ustrwm = ustrc;
   m.fw = 0;
   m.zoa = zo;
   return
end
ustrc = ucr * vk / log(zr/zo); 
azo = ubr/(wr*zo);
fw = 1.39*(azo)^-0.52;                                         % Eqn 62a
tauw = 0.5*fw*ubr.^2;                                          % Eqn 57
taum = ustrc*ustrc*...
       (1+1.2*(tauw/(ustrc*ustrc+tauw)).^(3/2));               % Eqn 69
taumax = sqrt((taum+tauw*cos(phiwc)).^2+(tauw*sin(phiwc)).^2); % Eqn 70

m.ustrw = sqrt(tauw);
m.ustrc = sqrt(taum);
m.ustrwm= sqrt(taumax);
m.fw = fw;
m.zoa = exp( log(zr)-ucr*vk/m.ustrc )';
return


