function [xpv,ypv]=ll2pv(lon,lat)
% Convert from lat/lon to PV coordinates
% [xpv,ypv]=ll2pv(lon,lat)
rot = 56;
xoff = 379900;
yoff = 3722150;
zone=11;
[xutm,yutm]=ll2utm(lon,lat,zone);
[r,az]=pcoord(xutm-xoff,yutm-yoff);
az = az+rot;
[xpv,ypv]=xycoord(r,az);