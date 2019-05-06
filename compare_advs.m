ncload('../9916advs-cal.nc'); % load the statistics file
jt = time+time2/(3600*24*1000);
dn = j2dn(time,time2);

figure(1);clf
subplot(411)
h1=plot(dn,Hdg_1215,'.');
hold on
ylabel('Heading [degT]')

subplot(412)
h2=plot(dn,Ptch_1216);
hold on
h3=plot(dn,Roll_1217);
ylabel('Pitch/Roll [deg]')

subplot(413)
h4=plot(dn,P_4023);
hold on
ylabel('Pressure [dBar]')

subplot(414)
h5=plot(dn,brange);
hold on
ylabel('brange [m]')

ncload('../9917advs-cal.nc'); % load the statistics file
jt = time+time2/(3600*24*1000);
dn = j2dn(time,time2);

figure(1);
subplot(411)
h1b=plot(dn,Hdg_1215,'.');
legend([h1;h1b],'9916','9917');

subplot(412)
h2b=plot(dn,Ptch_1216);
h3b=plot(dn,Roll_1217);

subplot(413)
h4b=plot(dn,P_4023);

subplot(414)
h5b=plot(dn,brange);
