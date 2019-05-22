function pf = pfit(C,z,iplot,za)
% pf = pfit(C,z,iplot,za)
%    Fit Rouse parameter to concentration profiles
%
% Fits y = a + bx, where
%   y = ln(C)
%   x = ln(z/z(1))
%   a = ln(C(z(1)))
%
% Input
%   C - Column vector of concentrations at z
%   z - Column vector of elevations (z(1)=lowest)
%   iplot - 1=plot data and fit (optional)
%   za - Reference height for Ca (optional...otherwise lowest height)
%
% Returns:
%    pf.N = length(C);
%    pf.p = b;
%    pf.Ca = exp(a);
%    pf.sa = sa;
%    pf.sb = sb;
%    pf.za = za;
%    pf.r2 = r2;


% See Appendix B in Lynch et al., 1994. Cont. Shelf Res.
% 14(10/11):1139-1165
%
% csherwood@usgs.gov
% last modfied 2-Feb-2012
% comments edited 4-Dec-2014

if(exist('iplot')~=1),iplot=0;,end
if(exist('za','var')~=1),za = z(1);, end;
y = log(C);
x = log(z./za);
[a,b,r2,sa,sb]=lsfit(x,y);
zest = logspace( log10(za),log10(max(z)), 20);
Cest = exp( a + b*log(zest./za) );
if(iplot)
   loglog(C,z,'ob');
   hold on
   loglog(Cest,zest,'-r');
end
pf.N = length(C);
pf.p = b;
pf.Ca = exp(a);
pf.sa = sa;
pf.sb = sb;
pf.za = za;
pf.r2 = r2;



