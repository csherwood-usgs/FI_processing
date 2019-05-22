% puv_proc_nc.m - Process ADV bursts
% Last revised by CRS on desktop 8/1/2005

kN1 = 30*0.01e-2 % roughness to be used on GM model
kN2 = 30*0.5e-2

% atmos = 10. % atmospheric pressure in decibars
% load atmospheric pressure from combined Ancona and Falconara data
%load atmos_comb;

load atmosp758 % contains ap (millibars) and jt_ap from buoy 46025
% Palos Verdes B6 tripod
dirname = [''] % for RMH use
dirname = ['g:\data\PV_2004\'] % for use on CRS desktop

advbfn = [dirname,'adv7582_2b-cal.nc'];
advsfn = [dirname,'adv7582_2vp-cal-vfix.nc'];
ncload(advsfn);
nc = netcdf(advbfn);

gb = (1:length(burst)); %.nc files have now been cut to gb only  - RMH 7/28/05
% check this CRS - elevation of paros above adv sample elev
pdelz = 1.375; % was .6, maybe that # was from different deployment? RMH

depth = (P_4023-ap)/100;
jt = time+time2/(3600*24*1000);
zr = fillnan(jt, vrange); %vrange has already been ydragged, but there are still NaN's from ecorr
%fs = 10;
%nsamp =7200;
fs = nc.ADVDeploymentSetupSampleRate(:);
nsamp = nc.ADVDeploymentSetupSamplesPerBurst(:);
for n=1:length(burst),
  u = nc{'u_1205'}(n,:)'/100; %? - yes, u was in cm/s RMH
  v = nc{'v_1206'}(n,:)'/100;
  p = nc{'P_4022'}(n,:)'/100; %? - yes, p was in mbar RMH
  PUV(n) = puv( p(19:nsamp), u(19:nsamp), v(19:nsamp), ...
		depth(n), zr(n)+pdelz, zr(n), fs, 1024, 1035 );
  UBS(n) = ubstats( u, v, fs );
  ubr= [PUV(n).ubr];
  wr = [PUV(n).omegar];
  ucr = [UBS(n).mean_spd];
  phiwc = (pi/180)*( [UBS(n).mean_dir] - [UBS(n).maj_az]);
  % make sure all these inputs are reasonable CRS
  M1(n) = m94( ubr, wr, ucr, zr(n), phiwc, kN1, 0 );
  M2(n) = m94( ubr, wr, ucr, zr(n), phiwc, kN2, 0 );
end
ncclose

save puv_proc_08_01.mat ap depth PUV zr UBS M1 M2 jt kN1 kN2