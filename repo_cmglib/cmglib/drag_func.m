function [zo, Cd] = drag_func( u, ustr, zr )
% drag_func - Calculate zo and Cd from u, u*, and zr
% [zo, Cd] = drag_func( u, ustr, zr )
vk = 0.41;
u = u(:);
ustr = ustr(:);
vk = vk*ones(size(u));
zr = zr(:);
if(length(zr)==1)
   zr=zr*ones(size(u))
end
zo = exp( -(  (vk.*u./ustr) - log(zr) ) );
Cd = (vk./log(zr./zo)).^2;