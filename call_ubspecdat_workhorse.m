
clear all ; close all ; clc; 
 wh=fullfile('/media/taran/DATADRIVE2/Obs_data/data_netcdf/9921whp-cal.nc'); % statistics filename
% 
  netcdf_load(wh)
  
h_1d(:)=double(hght_18(1,1,:));
s(:,:)=double(sspec(1,1,:,:));  
band_width=0.015625  ;
f=double(frequency(:,1));
 df=band_width ;
% Convert spectra from mm to m 
 
for t=1:length(h_1d); 
 % for ft=1:length(f)
    [ubr(t),Tbr(t)]=ubspecdat(squeeze(h_1d(t)),0.001*s(:,t)',f(:,1)',df);
end 
save('ubr_from_spectra.mat','ubr','Tbr') 