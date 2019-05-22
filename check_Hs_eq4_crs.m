clear all ; close all ; clc ;

% From work horse -------------> 
% verify the significant wave height based on eq 4 of CRS paper
wh=fullfile('/media/taran/DATADRIVE2/Obs_data/data_netcdf/9921whp-cal.nc'); % statistics filename
% 
netcdf_load(wh); 
% equation 4 of finding significant wave height 

band_width=0.015625  ;
f=frequency ; 

for t=1:length(wh_4061)
  summation_Sn_dfi=0.0 ; 
 % for i=2:length(f)
 %   Sn(i)=double(sspec(1,1,i,t)) ; 
 %   df(i)=double(f(i)-f(i-1))    ; 
%  df(i)=double(f(i)); 
%  f=frequency(:,1);
 %   summation_Sn_dfi=Sn(i).*df(i)+summation_Sn_dfi;  % add all the frequency of spectra
  %end 
  m0=trapz(double(f(1:end)),0.001*double(squeeze(sspec(1,1,1:end,t))) ); %% mm to m 
  Hs(t)=4.0*sqrt(m0); 
  Hrms_wh(t)=Hs(t)/sqrt(2.0) ; 
end 



% From ADV directly ------------->
load('crs_adv_9917.mat','dn','jtb_rec','depth','Hrmsu',..........
          'Tr','Ubr','Ur')
      
 nt1=680; nt2=730; 
 
%plot(Hrmsu(nt1:nt2),Hrms_wh(nt1:nt2))

plot(dn(nt1:nt2),Hrmsu(nt1:nt2),'linewidth',2); 
hold on 
plot(dn(nt1:nt2),Hrms_wh(nt1:nt2),'linewidth',2); 

 xlim([dn(nt1) dn(nt2)]);
  datetick('x',2) % 'keeplimits')
legend('Hrms Adv','Hrms workhorse')
%xlim([dn(nt1) dn(nt2)]);
%  datetick('x',2) % 'keeplimits')
