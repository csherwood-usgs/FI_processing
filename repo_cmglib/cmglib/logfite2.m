function [us,zo,r2,use,zoe_crs,res]=logfite2(u,z,iplot);
% LOGFITE2 - Log profile fits and error bars
% [us,zo,r2,use,zoe,res]=logfite2(u,z,iplot);
% 
% Last revised by Chris Sherwood, USGS, June 24, 2003

if(exist('iplot')~=1),iplot=0;,end;
vk = 0.41;
n = length(z);

y = abs(u);
x = log(z);

% Chris' method (Sachs, 1984, p. 417)
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
%sa = sydotx*(sqrt(1/n+xbar^2/Qx));

% Jessies method
Z = [ log(z) ones(n,1) ];
A =  Z\u;
us  = vk*A(1);
zo  = exp(-A(2)/A(1));
uhat = (us/vk)*log(z./zo);
C = cov([ u log(z) ]);
r2 = ( C(1,2)./(std(u)*std(log(z))) ).^2;

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

use_crs = vk*t*sb;
% this value has been fixed
zoe_crs = abs(-a/b)*sqrt( (t*sa/a).^2 + (t*sb/b).^2 );

res = u-uhat;
if(iplot>0),
  figure(iplot); clf;
  plot( log10(z),u,'ok',log10(z),uhat,'-r' );
  drawnow
end
if (0),
   fprintf(1,'Jessie:\nust: %g+-(%g)\nzo: %g+-(%g)\nr2: %g\n',...
      us,use,zo,zoe,r2);
   fprintf(1,'Chris:\nust: %g+-(%g)\nzo: %g+-(%g)\nr2: %g\n',...
      us_crs,use_crs,zo_crs,zoe_crs,r2_crs);
end















