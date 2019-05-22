function [unew, vnew] = new_north( u, v, az )
% new_north - Rotate vectors so vnew points toward az
% [unew, vnew] = new_north( u, v, az )
[s,d]=pcoord(u,v);
d = d-az;
[unew, vnew] = xycoord( s, d );