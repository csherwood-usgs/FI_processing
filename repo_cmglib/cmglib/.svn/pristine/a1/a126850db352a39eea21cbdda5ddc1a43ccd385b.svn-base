function [ubr,Tbr,Hs] = ubspecdat_hs(h,s,f,df)
% UBSPECDAT_HS - Calculate ubr Tbr and Hs from measured spectra
% [ubr,Tbr,Hs] = ubspecdat(h,s,f,df)
%
% Input:
%   h = water depth (m) - scalar or col. vector with length
%   s(nf) or s(nt,nf) = array of spectral densities normalized so that
%   Hs = 4*sqrt(sum(s*df,2))
%   f = row vector with central frequencies (Hz)
%   df = (optional) scalar or row vector with freq. bandwidths (Hz)
% Returns:
%   ubr = representative bottom orbital velocity (m/s)
%   Tbr = representative bottom wave period (s)
%
% Equation numbers refer to Wiberg and Sherwood, 2007, Computers and
% Geosciences, in press.
% The alternative bottom period, Tbz, is also calculated (see text).
%
% Chris Sherwood, USGS
% Modified by Benedicte Ferre to add Hs output
% Last revised May 10, 2007

% g = 9.81;
[nt,nf]=size(s);
% T = 1 ./f;
w = 2*pi*f;
% Determine kh using Soulsby (2006) method (see Appendix 5)
kh = qkhfs(w,h);
w = repmat(w,nt,1);
fm = repmat(f,nt,1);
kh = repmat(kh,nt,1);
if(exist('df')),
   if(length(df)==1),
      df = df*ones(nt,nf);
   elseif(length(df)==nf);
      df = repmat(df,nt,1);
   end
else
   fi = f(1)-(f(2)-f(1));
   fe = f(end)+(f(end)-f(end-1));
   dfb = diff([fi f(:)']);
   dff = diff([f(:)' fe]);
   df = mean([dfb;dff]);
   df = repmat(df,nt,1);
%    dff = diff(f);
%    dff = [dff(1) dff];
%    df = repmat(dff,nt,1);
end
Su = ((w.^2)./((sinh(kh)).^2)).*s; %Eq. 5
ubr = sqrt(2.*sum((Su.*df)')' ); %Eq. 6
fr = (sum((Su.*fm.*df)')') ./ (sum((Su.*df)')'); %Eq. 9
Tbr = 1./fr;
fz = sqrt((sum((Su.*fm.^2.*df)')') ./ (sum((Su.*df)')'));
Tbz = 1./fz;
Hs=4*sqrt(sum(s.*df,2));
