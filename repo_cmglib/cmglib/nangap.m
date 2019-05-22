function [xi,yi]=nangap(x,y,fac)
% NANGAP - Finds gaps in x and insert new xi, yi pairs with yi = NaN
% [xi,yi]=nangap(x,y,[fac])
%
% Values x = x+dt and y = NaN
% This flags gaps and allows Matlab line plots to disconnect.
% Criteria for a gap is median(dx)*fac [default fac = 1.2]

% csherwood@usgs.gov
x = x';
y = y';
if(exist('fac')~=1),fac = 1.2;,end;

dt = diff(x);
mdt = median(dt);
igap=find(dt>mdt*fac);

xigap = x(igap)+mdt;
yigap = NaN*xigap;

[xi,i]=sort([x;xigap]);
yi = [y;yigap];
yi = yi(i);
return
