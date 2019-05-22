  clear all ; close all ; clc ; 
% % compare the skewness calculated from Ursell number that adv predicts
% % WITH velocity skewness directly from adv burst data 
 load('crs_adv_9917.mat','dn','jtb_rec','depth','Hrmsu',..........
          'Tr','Ubr','Ur','velu_skew','velv_skew')
dn_adv=dn; 
dnsb_adv = datestr(datenum(dn_adv))   ;
Ur_adv=Ur; 
ok = find(~isnan(depth));
 
nt1=680; nt2=730; 

% rp = ruessink_asymm( Ur );
%Su = rp.Su;
%Au = rp.Au;
%r = rp.r;

wh=fullfile('/media/taran/DATADRIVE2/Obs_data/data_netcdf/9921whp-cal.nc'); % statistics filename
% 
 netcdf_load(wh)

jt = time+time2/(3600*24*1000);
dn = j2dn(time,time2);

Hs(:)=squeeze(wh_4061(1,1,:));
Td(:)=squeeze(wp_peak(1,1,:)); 
hght(:)=squeeze(hght_18(1,1,:)) ; 

nt1=680; nt2=730; 
 
 for i=1:length(Hs)
     if (Hs(i)>100);
        Hs(i)=0.0;
     end
     if (Td(i)>30); 
         Td(i)=0.0;
     end 
  kh(i)   = qkhfs( Td(i), hght(i) );
  k(i)    =kh(i)./hght(i);
  
  Hs_wh(i)=Hs(i); 
  Ur_wh(i)=0.75*0.5*1.4*Hs_wh(i)*k(i)./(kh(i).^3); % RRvR Eqn. 6 end
 end 

% Get Ursell number directly from the calculated Hrms,depth, Tr of the ADV
% % data
 figure(1)
 plot(dn(nt1:nt2),Td(nt1:nt2),'linewidth',2);
 hold on
 plot(dn(nt1:nt2),Tr(nt1:nt2),'linewidth',2); 
 xlim([dn(nt1) dn(nt2)]);
 datetick('x',2) % 'keeplimits')
 legend('adv-Ubr from burst data','workhorse-linear wave theory')
 ylabel('Ur')     
% %       
% 
% figure(1)
% plot(dn(nt1:nt2),Su(nt1:nt2),'linewidth',2);
% hold on   
% plot(dn(nt1:nt2),velu_skew(nt1:nt2),'linewidth',2);
%  xlim([dn(nt1) dn(nt2)]);
%   datetick('x',2) % 'keeplimits')
%   legend('Su','velu skew')