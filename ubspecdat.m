function [ubr,Tbr,Su,df] = ubspecdat(h,s,f,df)
% UBSPECDAT - Calculate ubr and Tbr from measured spectra
% [ubr,Tbr] = ubspecdat(h,s,f,df)
%
% Input:
%   h = water depth (m) - scalar or col. vector with length
%   s(nf) or s(nt,nf) = array of spectral densities normalized so that
%   Hs = 4*sum(s,2)*df
%   f = row vector with central frequencies (Hz)
%   df = (optional) scalar or row vector with freq. bandwidths (Hz)
% Returns:
%   ubr = representative bottom orbital velocity (m/s)
%   Tbr = representative bottom wave period (s)
%
% Equation numbers refer to Wiberg and Sherwood, 2008, Computers and
% Geosciences, 34:1243-1262..
% The alternative bottom period, Tbz, is also calculated (see text).
%
% Chris Sherwood, USGS
% Last revised July 24, 2009 (updated comments only)

% g = 9.81 assumed by qkhfs
[nt,nf]=size(s);
% T = 1 ./f;
w = 2*pi*f;
% Determine kh using Soulsby (2006) method (see Appendix 5)
kh = qkhfs(w,h);
w = repmat(w,nt,1);
fm = repmat(f,nt,1);
%fm=f; 
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
end 


Su = ((w.^2)./((sinh(kh)).^2)).*s; %Eq. 5
ubr = sqrt(2.*sum((Su.*df)')' ); %Eq. 6
fr = (sum((Su.*fm.*df)')') ./ (sum((Su.*df)')')  ;%Eq. 9
Tbr = 1./fr;
fz = sqrt((sum((Su.*fm.^2.*df)')') ./ (sum((Su.*df)')'));
Tbz = 1./fz;