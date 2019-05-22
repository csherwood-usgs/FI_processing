function [xutm,yutm]=pvutm(xpv,ypv);
% pv2utm - Convert from PV coords to UTM (zone 11)
% [xutm,yutm]=pvutm(xpv,ypv)
rot = 56;
xoff = 379900;
yoff = 3722150;
zone = 11;
[r,az]=pcoord(xpv,ypv);
az = az-rot;
[xr,yr]=xycoord(r,az);
xutm=xr+xoff;
yutm=yr+yoff;