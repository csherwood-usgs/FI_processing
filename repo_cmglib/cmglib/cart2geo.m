function gdeg = cart2geo( cdeg );
% Convert from Cartesian angle (degrees) to geographic angle (degrees)
% Cartesian angles are measured CCW from x axis (Matlab)
% Geographic angles are measured CW from y axis (e.g., compasses)
%
% See also geo2cart.m
gdeg = 90-cdeg;
gdeg(gdeg<0)=gdeg(gdeg<0)+360;
gdeg(gdeg>=360)=gdeg(gdeg>=360)-360;