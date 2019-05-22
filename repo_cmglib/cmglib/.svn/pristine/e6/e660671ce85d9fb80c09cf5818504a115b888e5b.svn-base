function ripple=hlpred(D,do)
% HLPRED - Predict ripple dimensions according to Wiberg & Harris, 1994.
% ripple=hlpred(D,do)
% Input:
%   do orbital diameter [m]
%   D  grain size [m]
% Outputs expected ripple dimensions, lam [m] and H [m]

% Based on Wiberg, P.L. and C.K. Harris, 1994. Ripple geometry in
% wave-dominated environments. Journal of Geophysical Research, 99(C1),775-789.
%
% Coded by ckharris, UVa
% Converted to Matlab by C. Sherwood, USGS
% Last revised January, 2005
% Assume anrbital type, and calculate expected ripple dimensions
%
ripple.type_no = -9;
ripple.type_text = 'none yet';
ripple.H = -9.;
ripple.lam = -9.;
ripple.stp = -9.;
ripple.D = D;
ripple.do = do;
D = D*1000;
do = do*1000;
%
% Initial guess of lam and steepness
% Assume anorbital
lam=535.*D;
stp = stp_clc(lam,do,D,0.1);
H = stp*lam;
if (H>0.),
  doH=do/H;
else
  doH=0.;
end

% Now, decide ripple type
if (stp<=0.01 | doH <= 0.),
  ripple.type_no = 6;
  ripple.type_text='upper plane bed';
elseif (doH>0 & doH<20.),
  ripple.type_no = 1;
  ripple.type_text='orbital';
elseif ( doH>=20 & doH<=35 ),
  ripple.type_no = 3;
  ripple.type_text='suborbital/orbital';
elseif (doH>35 & doH<=65),
  ripple.type_no = 3;
  ripple.type_text='suborbital';
elseif (doH>65. & doH<100),
  ripple.type_no = 3;
  ripple.type_text='suborbital/anorbital';
elseif (doH>=100.),
  ripple.type_no = 2;
  ripple.type_text='anorbital';
else
  error('Fell through all cases in hlpred.m')
end

if(ripple.type_no == 1),
  % orbital, need new calcs
  ripple.lam=0.62*do;
  ripple.stp = 0.17;
  ripple.H = stp*ripple.lam;
elseif (ripple.type_no == 2),
  % anorbital, as prev. calculated
  ripple.lam=lam;
  ripple.H = H;
  ripple.stp = stp;
elseif (ripple.type_no == 3),
  % suborbital, interpolate based on weighted average		
  frac=(log(doH)-log(100.)) / (log(20.)-log(100.));
  lam=exp(frac*(log(0.62*do)-log(lam))+log(lam));
  stp=stp_clc(lam,do,D,0.1);
  ripple.lam = lam;
  ripple.stp = stp;
  ripple.H=stp*lam;
elseif( ripple.type_no == 6),
  ripple.H = 0.;
  ripple.stp =0.;
else
  error('Fell through second set of cases')
end
ripple.H = ripple.H/1000;
ripple.lam = ripple.lam/1000;
return


function stp_out = stp_clc(lam,do,D,stp)
MAXIT = 30;
stp_old = stp;
for ni=1:MAXIT
  % Iterate to consistent steepness, u*, z0.
  H = stp * lam;
  if (H<=(3.*D)),
    H = 0.0;
    stp_out = 0.0;
    return
  end
  %
  % Recalculate steepness based on a regression relationship
  doH = do/H;
  ln_doH = log(do/H);
  stp_new = exp(-0.0950*ln_doH^2 + ...
		0.4415*ln_doH-2.2827);
  % Make sure that new steepness is reasonable, less than .17
  if (stp_new > 0.17 | doH <= 10.),stp_new=0.17; end
  % See if steepness has converged
  if (((abs(stp_new-stp)/stp)<0.005) | (ni == MAXIT)),
    stp_out = stp_new;
    return
  else
    if(ni==1) stp_old = stp; , end
    % Pick new steepness
    if ((stp_old<stp_new) & (stp<stp_new)),
      stp_old = stp;
      stp = stp_new+(stp_new-stp);
    elseif ((stp_old>stp_new) & (stp>stp_new)),
      stp_old = stp;
      stp = stp_new+(stp_new-stp);
    elseif ((stp_old<stp_new) & (stp>stp_new)),
      stp_old = stp;
      stp = 0.5*(stp_new+stp);
    elseif ((stp_old>stp_new) & (stp<stp_new)),
      stp_old = stp;
      stp = 0.5*(stp_new+stp);
    end
    if (stp>0.17), stp = 0.17; end
  end
end
return
