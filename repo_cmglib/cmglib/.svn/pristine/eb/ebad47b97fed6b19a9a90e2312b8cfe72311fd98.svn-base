function mu = visc_tsp( T, S, P )
% VISC_TSP - Dynamic (absolute) viscosity of water as function of temp, salinity, pressure
% mu = visc_tsp( T, S, P )
%
% Input:
%   T - Temperature (degrees centigrade)
%   S - Salinity (psu)
%   P - Pressure ( decibars ~= m water depth)
% (any one can be a vector)
%
% Returns:
%  mu - Viscosity in centipoise ( 1e-2 g/(cm*s) )
%       (divide by 1000 to get MKS units of Pa*s)

% Matthaus, cited in Kukulka, D.J., Gebhart, B, and Mollendorf,
% J. C. (1987) Thermodynamic and transport
% properties of pure and saline water. Adv. Heat Transfer 8:325-363.
% Cited in Boudreau, Diagenetic Models and Their Implementation, p. 93-94.
if(exist('T')~=1),T=20;,end
if(exist('S')~=1),S=0;,end
if(exist('P')~=1),P=0;,end
P = P./10.; % convert decibars -> bars
mu = 1.7910-(6.144e-2).*T + (1.451e-3).*(T.^2) - (1.6826e-5).*(T.^3) ...
     -(1.5920e-4).*P + (8.3885e-8).*(P.^2) + (2.4727e-3).*S ...
     + T.* ((6.00574e-6).*P - (2.6760e-9).*(P.^2)) ...
     + S.*((4.8429e-5).*T - (4.7172e-6).*(T.^2) + (7.5986e-8).*(T.^3)); 