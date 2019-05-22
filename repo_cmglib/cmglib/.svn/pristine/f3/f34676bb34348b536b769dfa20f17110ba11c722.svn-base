function ynew = fillnanends(x,y);
% FILLNANENDS - Fills NaN gaps in time series with linear interpolation,
% pads ends with first/last good value
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
         disp('FILLNANENDS found NaNs at ends of series')
      else
         dydx = (y(nextgood,j)-y(lastgood,j))/(x(nextgood)-x(lastgood));
         ynew(list(i),j) = y(lastgood,j)+ (x(list(i))-x(lastgood))*dydx;
      end
   end
   if( isnan( ynew(1,j) ) )
      ifg = find( isfinite(y(:,j)),1,'first');
      ynew(1:ifg-1,j)=y(ifg,j);
      fprintf(1,'FILLNANDENDS: padded beginning with %d values\n',ifg-1);
   end
   if( isnan( ynew(end,j) ) )
      ifg = find( isfinite(y(:,j)),1,'last');
      ynew(ifg+1:end,j)=y(ifg,j);
      fprintf(1,'FILLNANDENDS: padded end with %d values\n',nr-ifg);
   end
end
nbadout=sum(isnan(ynew(:)));
fprintf('FILLNANENDS: NaNs in input: %d; in output: %d\n',nbad,nbadout)
return