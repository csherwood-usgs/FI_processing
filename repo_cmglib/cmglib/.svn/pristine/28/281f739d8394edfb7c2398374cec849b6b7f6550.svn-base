function Cz = rouse(z,za,h,ws,ustrb)
% ROUSE  Classic Rouse concentration profile
%
% Cz = rouse(z,za,h,ws,ustrb)
%
%   Computes suspended sediment concentration
%   at heights in vector z according to the
%   Rouse equation for constant-stress eddy
%   diffusivity profile.
%
%   z     = vector of elevations above bottom
%   za    = reference elevation
%   h     = water depth
%   ws    = settling velocity (negative for sediment)
%   ustrb = bottom friction velocity
%
%   Multiply results by Cref = Cz(za)

% Chris Sherwood, USGS
% March 17, 1999

vk = 0.41;   % von Karman's constant
smallz = find(z<za);
Cz = (((h-z)./z) * za/(h-za)).^(-ws /(vk*ustrb));
Cz(smallz)=ones(size(smallz));