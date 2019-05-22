function h=circle(r,x,y,pstring)
% CIRCLE - Draws a circle centered at x,y with radius r
%
% function h=(r,x,y,pstring)
%   pstring is optional string to describe line type
%   Returns handle to arguement

% Chris Sherwood, USGS
% Last revised 15 June 99
% comment to test svm
len = 100;
az = [0:5:360]';
r1 = ones(length(az),1)*r;
[x1,y1] = xycoord(r1,az);
x1 = x1+x;
y1 = y1+y;
if (nargin == 3),
   h=plot(x1,y1,'-b');
end
if (nargin > 3),
   h=plot(x1,y1,pstring);
end
