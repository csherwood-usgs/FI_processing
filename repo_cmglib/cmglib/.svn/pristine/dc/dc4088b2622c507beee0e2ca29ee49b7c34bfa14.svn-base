function ynew = fillnan(x,y);
% FILLNAN - Fills NaN gaps in time series with linear interpolation
% ynew = fillnan(x,y);
%

% Chris Sherwood, USGS
if(any(isnan(x))),
   error('Cant cope with NaNs in x')
end
[nt,nc]=size(x);
ynew = y;
if(~any(isnan(y))),return,end

[nr,nc]=size(y);
for j=1:nc,
   list = find(isnan(y(:,j)));
   nbad = length(list);
   for i=1:nbad,
      lastgood = max(1, list(i)-1);
      while((lastgood > 1) & isnan(y(lastgood,j))),
         lastgood = lastgood-1;
      end
      nextgood = min(nt, list(i)+1);
      while((nextgood < nt) & isnan(y(nextgood,j))),
         nextgood = nextgood+1;
      end
      if((lastgood==1 & isnan(y(1,j)))|(nextgood==nt&isnan(y(nt,j)))),
         disp('FILLNAN: Cant fix problem at ends of series')
      else
         dydx = (y(nextgood,j)-y(lastgood,j))/(x(nextgood)-x(lastgood));
         ynew(list(i),j) = y(lastgood,j)+ (x(list(i))-x(lastgood))*dydx;
      end
   end
end
nbadout=sum(isnan(ynew(:)));
fprintf('FILLNAN: NaNs in input: %d; in output: %d\n',nbad,nbadout)
return