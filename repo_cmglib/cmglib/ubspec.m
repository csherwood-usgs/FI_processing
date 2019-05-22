function ubs = ubspec( h, s, f, method, df )
% UBSPEC - Representative near-bottom wave-orbital velocity and period
% ubs = ubspec( h, s, f, method, df )
%
% Input:
%   h = water depth (m) - scalar or col. vector with length
%   s(nf) or s(nt,nf) = array of spectral densities normalized so that
%                       Hs = 4*sqrt(sum(s*df))
%   f       = row vector with central frequencies (Hz)
%   method  = (optional) 'iterate' (slow but precise wavenumber calc) 
%                     or 'default' (fast and good enough)
%   df      = (optional) scalar or row vector with freq. bandwidths (Hz)
%
% Returns;
%   ubs.Hs - Significant wave height (to check scaling)
%   ubs.ubr - Representative wave-orbital velocity (m/s)
%   ubs.Tr - representative period (s)

% Chris Sherwood, USGS
% Revised to return a structure Nov 2005
% Last revised August 3, 2005

g = 9.81;
[nr,nc]=size(h);

if(exist('method')~=1),method='default',end
[nt,nf]=size(s);
T = 1 ./f;
w = 2*pi*f;
if(strcmp(lower(method),'iterate')),
  for i=1:length(T),
    k(i) = waven( T(i), h );
  end
  kh = (k*h);
else,
  kh = qkhf(w,h);
end
w = repmat(w,nt,1);
T = repmat(T,nt,1);
kh = repmat(kh,nt,1);
h = h*ones(nt,nf);

if(exist('df')),
  if(length(df)==1),
    df = df*ones(nt,nf);
  elseif(length(df)==nf);
    df = repmat(df,nt,1);
  end
else
  df = diff(f(:)');
  df = [df df(end)];
  df = repmat(df,nt,1);
end
Hs = 4*sqrt(sum(s.*df)); % this matches NODC estimates
amp = sqrt(2 .* s);
ubms = (w .* amp) ./ (sinh(kh));
ubr = sqrt( sum((ubms.^2 .*df)')' );
Tr = (sum(((ubms.*df).*T)')') ./ (sum((ubms.*df)')');
ubs.Hs = Hs;
ubs.ubr = ubr;
ubs.Tr = Tr;