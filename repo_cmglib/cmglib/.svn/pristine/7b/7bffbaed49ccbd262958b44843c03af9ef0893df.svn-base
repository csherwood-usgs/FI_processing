function [ustar,u1] = USTAR_UBAR( ubar, h, zo )
% USTAR_UBAR - Calculate u* and u at 1 meter from depth-mean velocity
% [ustar,u1] = USTAR_UBAR( ubar, h, zo )
%
% Input:
%   ubar - depth-mean velocity (m/s)
%   h    - water depth (m)
%   zo   - bottom roughness length (m)
%
% Returned:
%   ustar - friction velocity (m/s)
%   u1    - velocity at z = 1 meter (m/s)
%
% Assumes log profile u = ustar/vk * ln( z/zo )

% Chris Sherwood, USGS
% Last revised 15 November 2001

vk = 0.41;
ustar = vk*ubar / (log(h/zo)-1+zo/h);
u1 = ustar/vk * log(1/zo);
