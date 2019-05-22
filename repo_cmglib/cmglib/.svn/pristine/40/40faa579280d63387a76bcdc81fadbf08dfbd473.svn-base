function [us,zo,r2,r2a,use,zoe,res,P,P2]=mlogfit(u,z,iplot);
% MLOGFIT - Log profile fits and error bars; returns reduced r2 for
% quadratic fit
% [us,zo,r2,r2a,use,zoe,res,P,P2]=mlogfit(u,z,iplot);
% 
% Last revised by Chris Sherwood, USGS, January 20, 2005
if(exist('iplot')~=1),iplot=0;,end;
vk = 0.408;
n = length(z);

y = abs(u);
x = log(z);

% Chris' method (Sachs, 1984, p. 417)
if(0),
sumx = sum(x);
sumy = sum(y);
sumx2= sum(x.^2);
sumy2 = sum(y.^2);
sumxy = sum(x.*y);
Qx = sumx2-sumx^2/n;
Qy = sumy2-sumy^2/n;
Qxy = sumxy-(1/n)*sumx*sumy;
xbar = sumx/n;
ybar = sumy/n;
b = Qxy/Qx;
a = (sumy-b*sumx)/n;
us_crs = vk*b;
zo_crs = exp( -a/b );
r2_crs = (Qxy/sqrt(Qx*Qy)).^2;
Qydotx = Qy-b*Qxy;
sydotx = sqrt( Qydotx/(n-2) );
sb = sydotx/sqrt(Qx);
sa = sb*sqrt(sumx2/n);
end

% linear fit
[P,S]=polyfit(x,y,1);
yhat = polyval(P,x);
SSTO = sum((y-mean(y)).^2);
res = y-yhat;
SSE =  sum((y-yhat).^2);
r2 = 1-( (n-1)/(n-2) )*SSE/SSTO;

% quadratic fit
[P2,S2]=polyfit(x,y,2);
yhat2 = polyval(P2,x);
SSE2 = sum((y-yhat2).^2);
r2a = 1-( (n-1)/(n-3) )*SSE2/SSTO;

us = vk*P(1);
zo = exp( -P(2)/P(1) );

% Table 27 in Sachs (p. 136) for alpha = 0.1, two-sided test
studentt=[ 6.314 2.920 2.353 2.132 2.015,...
        1.943,1.895,1.860,1.833,1.812,1.796,1.782,1.771,...
        1.761,1.753,1.746,1.740,1.734,1.729,1.725,1.721,...
        1.717,1.714,1.711,1.708,1.706,1.703,1.701,1.699];
studentt=studentt(:);
DF = n-2;
t = 0;
if(DF>0),
if(DF<30), 
    t=studentt(DF);
else
    zalpha = 1.644854; % Table 43 in Sachs
    t = zalpha+(zalpha^3+zalpha)/(4*DF);
end
end

r = sqrt(r2); 
rr=sqrt((r^(-2)-1)/(n-2));

use=us*t*rr;
zoe=t*sqrt( (1/n)*sum(log(z).^2) )*rr; % jl
if(0),
use_crs = vk*t*sb;
zoe_crs = abs(-a/b)*sqrt( (t*sa/a).^2 + (t*sb/b).^2 );
fprintf(1,'pol us: %8f use: %8f zo: %8f zoe: %8f\n',us,use,zo,zoe)
fprintf(1,'crs us: %8f use: %8f zo: %8f zoe: %8f\n',us_crs,use_crs,zo_crs,zoe_crs)
end














