function xf = medfilt( x, m )
% MEDFILT - Moving median filter
% xf = medfilt( x, m )

if( nargin ~=2 ), m=5;, end;
if(mod(m,2)==0),error('m must be odd'),end
m2 = fix(m/2);
n = length(x);
x = x(:);
xf = x(:);

for i=(m2+1):(n-m2-1),
   xf(i) = nanmedian( x((i-m2:i+m2),1) );
end
