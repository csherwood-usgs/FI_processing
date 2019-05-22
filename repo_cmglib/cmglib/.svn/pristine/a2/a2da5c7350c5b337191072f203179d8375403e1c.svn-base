function d = shortdist(x,iplot)
% SHORTDIST - Returns distribution statistics (more bombproof than dist)

% csherwood@usgs.gov
nn = 0;
n = 0;
zmean = NaN;
zstd = NaN;
zmin = NaN;
zmax = NaN;
pct16 = NaN;
pct50 = NaN;
pct84 = NaN;
n = length(x);
nn = sum(~isnan(x));
if(nn>=3),
   z = sort(x);
   zmin = z(1);
   zmax = nanmax(x);
   index = 100*cumsum(ones(n,1))/n;
   %pct1 =  interp1(index,z,1);
   % pct5 =  interp1(index,z,5);
   pct16 =  interp1(index,z,16);
   %pct25 = interp1(index,z,25);
   pct50 = interp1(index,z,50);
   % pct75 = interp1(index,z,75);
   pct84 = interp1(index,z,84);
   % pct95 = interp1(index,z,95);
   %pct99 = interp1(index,z,99);
   zmean = nanmean(x);
   zstd = nanstd(x);

end
if(nn==2),
   zmean = nanmean(x);
   zmin = nanmin(x);
   zmax = nanmax(x);
end

fprintf('\n   N : %g   NonNaN: %d\n',n,nn);
fprintf('   Mean    : %g\n',zmean);
fprintf('   Std Dev : %g\n\n',zstd);
fprintf('   Minimum   : %g\n',zmin);
% %fprintf('    1-pctile : %g\n',pct1);
% fprintf('    5-pctile : %g\n',pct5);
fprintf('   16-pctile : %g\n',pct16);
%fprintf('   25-pctile : %g\n',pct25);
fprintf('   50-pctile : %g\n',pct50);
%fprintf('   75-pctile : %g\n',pct75);
fprintf('   84-pctile : %g\n',pct84);
% fprintf('   95-pctile : %g\n',pct95);
% %fprintf('   99-pctile : %g\n',pct99);
fprintf('   Maximum   : %g\n\n',zmax);
if ( nargin == 2 ),
   plot(z,index);
end

d.n = n;
d.nn = nn;
d.mean = zmean;
d.std  = zstd;
d.min  = zmin;
d.max  = zmax;
%d.pct5 = pct5;
d.pct16 = pct16;
d.pct50 = pct50;
d.pct84 = pct84;
%d.pct95 = pct95;

return

