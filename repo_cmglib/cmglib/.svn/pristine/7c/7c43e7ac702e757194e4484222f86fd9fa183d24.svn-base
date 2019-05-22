function mu = visc( T )
% VISC - Calculates dynamic viscosity of water from temperature T
% mu = visc( T )
%
% Input:
%   T - Temperature in degrees centigrade
% Returns:
%  mu - Dynamic viscosity in Poise ( g/(cm*s) )
%
% Multiply by 100 to get Centipoise
% Divide by 10 to get MKS units of Pa*s
% Divide by density to get kinematic viscosity in Stokes
%
% Uses linear interpolation between values from Bingham and
% Jackson (BULL. BUR. STDS.,11,75,(1918))

% Chris Sherwood, USGS
% Last revised Sept. 4, 2001

v = [1.7313,1.6728,1.6191,1.5674,1.5188,...
     1.4728,1.4284,1.3860,1.3462,1.3077,...
     1.2713,1.2363,1.2028,1.1709,1.1404,...
     1.1111,1.0828,1.0539,1.0299,1.0030,...
      .9810,.9579,.9358,.9142,.8937,...
      .8737,.8545,.8369,.8180,.8007];

if( T > 30. ),
   fprintf( 'Temperature must be less than 30C, Temp. = %f\n',T );
   mu = .8007/100.;
   return
end
if( T < 0. ),
   fprintf( 'Temperature is less than zero. Call a glaciologist. T= %f\n',T );
   mu = 1.7921/100.;
   return
end

itemp = floor( T );
delt = T - itemp;
if( T < 1. ),
   mu = 0.608 * T + 1.7921;
else
   mu = ( v(itemp-1)-v(itemp) )*delt + v(itemp-1);
end
mu = ( mu/100. );









