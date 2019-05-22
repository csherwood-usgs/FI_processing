clear all ; close all ; clc ; 
wh=fullfile('/media/taran/DATADRIVE2/Obs_data/data_netcdf/9921whp-cal.nc'); % statistics filename
% 
 netcdf_load(wh)

jt = time+time2/(3600*24*1000);
dn = j2dn(time,time2);

 Hs(:)=squeeze(wh_4061(1,1,:));
 Td(:)=squeeze(wp_peak(1,1,:)); 
 
 nt1=680; nt2=730; 
 
 plot(dn(nt1:nt2),Hs(nt1:nt2),'o')
 datetick('x',2) % 'keeplimits')