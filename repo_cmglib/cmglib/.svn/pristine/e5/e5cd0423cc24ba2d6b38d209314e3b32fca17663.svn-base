function [a,b,r2,sa,sb,n] = nanfit(x,y,iplot,t)
% NANFIT  Least-squares fit for non-NaN values
%
% [a,b, r2, sa, sb ] = lsfit(x,y,[iplot])
%
% Least-squares fit of
%   y = a + b * x
% x and y must be column vectors of same length.
% I think sa and sb are standard error of offset and slope.
%
% Makes a plot if iplot == 1, No plot if iplot is absent
% If provided, sd is plotted.
% r2 calculated according to Sachs(1982) p. 389

%    Chris Sherwood, UW and USGS
%    Last revised March 29,2004
echo off
format compact
if(nargin<3),iplot=0;,end;


x = x(:);
y = y(:);
noriginal = length(x);
ok = find( (~isnan(x(:)))&(~isnan(y(:))));
x = x(ok);
y = y(ok);

% matlab ls method
% x = [x ones(x) ];
% b = x\y

% or, equivalently
% b = (x'*x)\(x'*y)
% b = polyfit( log( z ./ za ), y, 1 )
% b = inv(x'*x)*(x'*y)
% predictions
% yhat = x*b;
% or
% yhat2 = polyval( b2, log( z ./ za ) );

% or, longhand
n = max(size(x));
if(noriginal > n),
  fprintf(1,'Eliminated %d values with x or y = NaN\n',noriginal-n);
end
sumx = sum(x);
sumy = sum(y);
sumx2 = sum(x .* x);
sumy2 = sum(y .* y);
sumxy = sum(x .* y);

% (see Taylor, p. 157)
% del = n*sumx2 - sumx*sumx;
% a = (1 ./del)*(sumx2*sumy - sumx*sumxy)
% b = (1 ./del)*(n*sumxy-sumx*sumy)
% yhat = [x ones(x)]*[ b; a ];
% varY = ( 1 ./(n-2) )*sum( (y-yhat).^2 );
% uncertainty in a = b(2)
% varA = varY * sumx2/del;
% uncertainty in b = b(1)
% varB = n*varY/del;
% sdp = sqrt( varB );
% r2  = sum( (yhat-mean(y)).^2 ) / sum( (y-mean(y)).^2 )

% (Sachs, p. 417)
Qx = sumx2-sumx^2/n;
Qy = sumy2-sumy^2/n;
Qxy = sumxy-(1/n)*sumx*sumy;

xbar = sumx/n;
ybar = sumy/n;

b = Qxy/Qx;
% or b = (sumxy-(sumx*sumy)/n)/(sumx2-sumx^2/n)
a = (sumy-b*sumx)/n;

sx = sqrt(Qx/(n-1));
sy = sqrt(Qy/(n-1));
r = Qxy/sqrt(Qx*Qy);
r2 = r*r;
Qydotx = Qy-b*Qxy;
sydotx = sqrt( Qydotx/(n-2) );
% or sydotx = sy*sqrt( (1-r2)*(n-1)/(n-2) );

% s.d. of intercept and slope
sb = sydotx/sqrt(Qx);
sa = sydotx*(sqrt(1/n+xbar^2/Qx));
%fprintf(1,'Slope: %f +- %f\nIntercept: %f +- %f\n',b,sb,a,sa);

t = 2.132; %n-2=4, two sided, alpha = .9

xline = linspace(min(x),max(x),20);
yline = a+b*xline;
syhatdot = sydotx*sqrt(1+1/n+((xline-xbar).^2)./Qx);
yplus = yline+syhatdot;
yminus = yline-syhatdot;
yhat = a+b*x;
residuals = (y - yhat)';


if(iplot==1),
hold on
%plot( x,y,'ok','markersize',8)
plot(xline,yline,'--k','linewidth',2)
plot(xline,yplus,':r','linewidth',1)
plot(xline,yminus,':r','linewidth',1)
end






