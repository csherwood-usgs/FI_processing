function [sd1,az1,sd2,az2,sd3,az3]=pcastats3(u,v,w,s,ipost)
% PCASTATS#  Principal components of 3-d (e.g. current meter) data
%
% function [sd1 az1 sd2 az2]=pcastats(u,v,[[s],ipost])
%
% u and v are column vectors with x and y data
% s is an optional length for the axes
% ipost is an optional flag 0=no action, >1=plot input data
 
% Chris Sherwood, USGS
% Sept 26, 2005

if(nargin <5),ipost=0;,end
if(nargin <4),ipost=0;s=0;,end

mu = mean(u);
mv = mean(v);
mw = mean(w);
C = cov([u' v' w']);
[V,D] = eig(C);

x1 = [.5*sqrt(D(1,1))*V(1,1);-.5*sqrt(D(1,1))*V(1,1)];
y1 = [.5*sqrt(D(1,1))*V(2,1);-.5*sqrt(D(1,1))*V(2,1)];
x2 = [.5*sqrt(D(2,2))*V(1,2);-.5*sqrt(D(2,2))*V(1,2)];
y2 = [.5*sqrt(D(2,2))*V(2,2);-.5*sqrt(D(2,2))*V(2,2)];
