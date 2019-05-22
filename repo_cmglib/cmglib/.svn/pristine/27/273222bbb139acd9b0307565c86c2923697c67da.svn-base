function K=cubick(z,h,ustrs,ustrb)
% CUBICK  Cubic eddy viscosity of Signell et al., 1990
%
% function K=cubick(z,h,ustrs,ustrb)

% Chris Sherwood, USGS
% March 17, 1999

vk = 0.41;
K = vk*ustrs*(h - z) + ...
     (vk*(ustrb-2.*ustrs)/h)*(h - z).^2 - ...
     (vk*(ustrb-ustrs)/(h*h))*(h - z).^3 ;
