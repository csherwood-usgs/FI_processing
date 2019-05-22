function [tauw, tauwe, tauwn] = wstress( w_sp, w_dir)
%[tauw, tauwe, tauwn] = wstress(jt, w_sp, w_dir)
%
% tauw = magnitude of wind stress
% tauwe = eastward component of wind stress
% tauwn = northward component of wind stress
%inputs:
%w_sp is a vector of wind speeds, in m/s
%w_dir is  direction wind is coming from in degrees, azimuthal
%north
%
% Cd based on Large and Pond 1981
% assumes wind measurements taken at 10 m elevation
% J. Lacy, USGS, 9-7-01

n= length(w_sp);
rhoa = 1.225;

cd = 0.0012*ones(n,1);
aa = find( w_sp >=11 &  w_sp <=25);
cd(aa)  = (0.49+0.065*w_sp(aa))*0.001;
bb = find( w_sp>25);
cd(bb) = 0.00212;

 
  tauw = rhoa*cd.*w_sp.^2;
  [tauwe, tauwn] = xycoord(tauw, w_dir);
  tauwe = -tauwe;
  tauwn = -tauwn;