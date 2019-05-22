function v = prctile(x,p)
% PRCTILE - Return p'th percentile of column x
%
% Quick and dirty version...no checks, multidimensions, or optimization!

% csherwood@usgs.gov Sept 15, 2005
n = length(x);
z = sort(x);
zmin = z(1);
zmax = z(n);
index = 100*cumsum(ones(n,1))/n;
%pct1 =  interp1(index,z,1);
v =  interp1(index,z,p);
return






