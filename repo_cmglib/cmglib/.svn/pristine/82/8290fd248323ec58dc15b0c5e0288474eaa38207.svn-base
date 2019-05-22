function d= matpctile(x,p)
% MATPCTILE - Generates a row vector of the p-th perctile values in
% each column of array x.
% d = matpctile(x,p)
%
% For ABS data where x(range,time), passing it (x',1) will return
% the profiles corresponding to the 1%-ile values for all time at each
% range.
%


% Chris Sherwood, USGS
% June 3, 2004

[nr,nc] = size(x);
z = sort(x); % ascending, by columns
zmin = z(:,1);
zmax = z(:,nc);
index = 100*cumsum(ones(nr,1))/nr;
d =  interp1(index,z,p);
