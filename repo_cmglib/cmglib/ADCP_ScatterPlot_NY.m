%ADCP_ScatterPlot_NY
%
%Program to do a scatter plot of eastward component vs. northward component
%of current from the New York Bight ADCP data.  Superimposes the ellipse
% with the principal axes.  Needs to be give a specific
%station, and creates a four-paneled plot with hourly-averaged and low-passed
%data from 5 mbs and the bottom-most depth bin.
%
%Soupy Alexander, 1/10/02

%Select a station

station = input('Select a station. ');

%Pull out the current data for the proper station
if station == 'A';
    file_name_h = '5951adc-a1h.nc';
    file_name_lp = '5951adc-alp.nc';
elseif station == 'B';
    file_name_h = '5971adc-a1h.nc';
    file_name_lp = '5971adc-alp.nc';
elseif station == 'C';
    file_name_h = '5991adc-a1h.nc';
    file_name_lp = '5991adc-alp.nc';
elseif station == 'D';
    file_name_h = '6011adc-a1h.nc';
    file_name_lp = '6011adc-alp.nc';
elseif station == 'E';
    file_name_h = '6031adc-a1h.nc';
    file_name_lp = '6031adc-alp.nc';
elseif station == 'F';
    file_name_h = '6041adc-a1h.nc';
    file_name_lp = '6041adc-alp.nc';
end

file_h = netcdf(file_name_h,'nowrite');
file_lp = netcdf(file_name_lp,'nowrite');

%Take care of the hourly averaged data:  pull out current from proper bin
t_h = file_h{'time'}(:);
t2_h = file_h{'time2'}(:);
time_h_orig = singleJD(t_h,t2_h);

depth_h = file_h{'depth'}(:);

index_5 = value2Index(depth_h,5);
index_deep = find(depth_h == max(depth_h));

depth_5 = depth_h(index_5);
depth_deep = depth_h(index_deep);

u5_h_orig = file_h{'u_1205'}(:,index_5);
v5_h_orig = file_h{'v_1206'}(:,index_5);

[u5_h,nh] = ridnan(u5_h_orig);
[v5_h,nh] = ridnan(v5_h_orig);
time_h_5 = time_h_orig(nh);

u5_h = ridfill(u5_h);
v5_h = ridfill(v5_h);

udeep_h_orig = file_h{'u_1205'}(:,index_deep);
vdeep_h_orig = file_h{'v_1206'}(:,index_deep);

[udeep_h,nh_d] = ridnan(udeep_h_orig);
[vdeep_h,nh_d] = ridnan(vdeep_h_orig);
time_h_deep = time_h_orig(nh_d);

udeep_h = ridfill(udeep_h);
vdeep_h = ridfill(vdeep_h);

%Calculate the principal ellipses and means.
theta=2*pi*(1:64)/64;
theta=[0 theta];

[majax5_h,majaz5_h,minax5_h,minaz5_h,ellip5_h] = pcaben(u5_h,v5_h);
[x15_h,y15_h]=cmgspd2uv(majax5_h,majaz5_h);
[x25_h,y25_h]=cmgspd2uv(minax5_h,minaz5_h);
xx5_h=(majax5_h*cos(theta))';
yy5_h=(minax5_h*sin(theta))';
angle5_h=-majaz5_h*pi/180+pi/2;
xxx5_h=xx5_h(:).*cos(angle5_h)-yy5_h.*sin(angle5_h);
yyy5_h=xx5_h(:).*sin(angle5_h)+yy5_h.*cos(angle5_h);

meanu5_h = nanmean(u5_h);
meanv5_h = nanmean(v5_h);
[mean5sp_h,mean5dir_h] = UVtoSpDir(meanu5_h,meanv5_h);

[majaxdeep_h,majazdeep_h,minaxdeep_h,minazdeep_h,ellipdeep_h] = pcaben(udeep_h,vdeep_h);
[x1deep_h,y1deep_h]=cmgspd2uv(majaxdeep_h,majazdeep_h);
[x2deep_h,y2deep_h]=cmgspd2uv(minaxdeep_h,minazdeep_h);
xxdeep_h=(majaxdeep_h*cos(theta))';
yydeep_h=(minaxdeep_h*sin(theta))';
angledeep_h=-majazdeep_h*pi/180+pi/2;
xxxdeep_h=xxdeep_h(:).*cos(angledeep_h)-yydeep_h.*sin(angledeep_h);
yyydeep_h=xxdeep_h(:).*sin(angledeep_h)+yydeep_h.*cos(angledeep_h);

meanudeep_h = nanmean(udeep_h);
meanvdeep_h = nanmean(vdeep_h);
[meandeepsp_h,meandeepdir_h] = UVtoSpDir(meanudeep_h,meanvdeep_h);

%Now take care of the low-passed data
t_lp = file_lp{'time'}(:);
t2_lp = file_lp{'time2'}(:);
time_lp_orig = singleJD(t_lp,t2_lp);

depth_lp = file_lp{'depth'}(:);

index_5 = value2Index(depth_lp,5);
index_deep = find(depth_lp == max(depth_lp));

depth_5 = depth_lp(index_5);
depth_deep = depth_lp(index_deep);

u5_lp_orig = file_lp{'u_1205'}(:,index_5);
v5_lp_orig = file_lp{'v_1206'}(:,index_5);

[u5_lp,nlp] = ridnan(u5_lp_orig);
[v5_lp,nlp] = ridnan(v5_lp_orig);

[u5_lp] = ridfill(u5_lp);
[v5_lp] = ridfill(v5_lp);
time_lp = time_lp_orig(nlp);

udeep_lp_orig = file_lp{'u_1205'}(:,index_deep);
vdeep_lp_orig = file_lp{'v_1206'}(:,index_deep);

[udeep_lp,nlp] = ridnan(udeep_lp_orig);
[vdeep_lp,nlp] = ridnan(vdeep_lp_orig);

udeep_lp = ridfill_nan(udeep_lp);
vdeep_lp = ridfill_nan(vdeep_lp);

%Calculate the principal ellipses and means.
theta=2*pi*(1:64)/64;
theta=[0 theta];

[majax5_lp,majaz5_lp,minax5_lp,minaz5_lp,ellip5_lp] = pcaben(u5_lp,v5_lp);
[x15_lp,y15_lp]=cmgspd2uv(majax5_lp,majaz5_lp);
[x25_lp,y25_lp]=cmgspd2uv(minax5_lp,minaz5_lp);
xx5_lp=(majax5_lp*cos(theta))';
yy5_lp=(minax5_lp*sin(theta))';
angle5_lp=-majaz5_lp*pi/180+pi/2;
xxx5_lp=xx5_lp(:).*cos(angle5_lp)-yy5_lp.*sin(angle5_lp);
yyy5_lp=xx5_lp(:).*sin(angle5_lp)+yy5_lp.*cos(angle5_lp);

meanu5_lp = nanmean(u5_lp);
meanv5_lp = nanmean(v5_lp);
[mean5sp_lp,mean5dir_lp] = UVtoSpDir(meanu5_lp,meanv5_lp);

[majaxdeep_lp,majazdeep_lp,minaxdeep_lp,minazdeep_lp,ellipdeep_lp] = pcaben(udeep_lp,vdeep_lp);
[x1deep_lp,y1deep_lp]=cmgspd2uv(majaxdeep_lp,majazdeep_lp);
[x2deep_lp,y2deep_lp]=cmgspd2uv(minaxdeep_lp,minazdeep_lp);
xxdeep_lp=(majaxdeep_lp*cos(theta))';
yydeep_lp=(minaxdeep_lp*sin(theta))';
angledeep_lp=-majazdeep_lp*pi/180+pi/2;
xxxdeep_lp=xxdeep_lp(:).*cos(angledeep_lp)-yydeep_lp.*sin(angledeep_lp);
yyydeep_lp=xxdeep_lp(:).*sin(angledeep_lp)+yydeep_lp.*cos(angledeep_lp);

meanudeep_lp = nanmean(udeep_lp);
meanvdeep_lp = nanmean(vdeep_lp);
[meandeepsp_lp,meandeepdir_lp] = UVtoSpDir(meanudeep_lp,meanvdeep_lp);

%Now plot the data
%Hourly Averaged, 5 mbs
subplot(2,2,1)
plot(u5_h,v5_h,'.')
hold on
plot([0 meanu5_h],[0 meanv5_h],'r','linewidth',2);
plot(xxx5_h+meanu5_h,yyy5_h+meanv5_h,'r','linewidth',2);
axis equal
axis tight
axis([-50 50 -50 50])
labels = [-50 -40 -30 -20 -10 0 10 20 30 40 50];
text(-45,-23,['Mean U = ' num2str(meanu5_h,'%0.2f') ...
        '; Mean V = ' num2str(meanv5_h,'%0.2f')]);
text(-45,-30,['Vector speed = ' num2str(mean5sp_h,'%0.1f') ... 
        '; Direction = ' num2str(angleXtoCompass(mean5dir_h),'%0.0f') '^o']);
text(-45,-38,['Major Axis = ' num2str(majax5_h,'%0.1f') ...
        '; Minor Axis = ' num2str(minax5_h,'%0.1f')]);
text(-45,-45,['Prin. Orient. = ' num2str(majaz5_h,'%0.0f') '^o']);
l = line([-50 50], [0 0]);
l2 = line([0 0], [-50 50]);
set(l,'color','k', 'linestyle','-.');
set(l2,'color','k', 'linestyle','-.');
set(gca,'xtick',labels);
set(gca,'ytick',labels);
set(gca,'xticklabel','');
set(gca,'yticklabel','');


%Low-Passed, 5 mbs
subplot(2,2,2)
plot(u5_lp,v5_lp,'.')
hold on
plot([0 meanu5_lp],[0 meanv5_lp],'r','linewidth',2);
plot(xxx5_lp+meanu5_lp,yyy5_lp+meanv5_lp,'r','linewidth',2);
axis equal
axis tight
axis([-50 50 -50 50])
labels = [-50 -40 -30 -20 -10 0 10 20 30 40 50];
text(-45,-23,['Mean U = ' num2str(meanu5_lp,'%0.2f') ...
        '; Mean V = ' num2str(meanv5_lp,'%0.2f')]);
text(-45,-30,['Vector speed = ' num2str(mean5sp_lp,'%0.1f') ... 
        '; Direction = ' num2str(angleXtoCompass(mean5dir_lp),'%0.0f') '^o']);
text(-45,-38,['Major Axis = ' num2str(majax5_lp,'%0.1f') ...
        '; Minor Axis = ' num2str(minax5_lp,'%0.1f')]);
text(-45,-45,['Prin. Orient. = ' num2str(majaz5_lp,'%0.0f') '^o']);
l = line([-50 50], [0 0]);
l2 = line([0 0], [-50 50]);
set(l,'color','k', 'linestyle','-.');
set(l2,'color','k', 'linestyle','-.');
set(gca,'xtick',labels);
set(gca,'ytick',labels);
set(gca,'xticklabel','');
set(gca,'yticklabel','');

%Hourly Averaged, Deep
subplot(2,2,3)
plot(udeep_h,vdeep_h,'.')
hold on
plot([0 meanudeep_h],[0 meanvdeep_h],'r','linewidth',2);
plot(xxxdeep_h+meanudeep_h,yyydeep_h+meanvdeep_h,'r','linewidth',2);
axis equal
axis tight
axis([-50 50 -50 50])
labels = [-50 -40 -30 -20 -10 0 10 20 30 40 50];
text(-45,-23,['Mean U = ' num2str(meanudeep_h,'%0.2f') ...
        '; Mean V = ' num2str(meanvdeep_h,'%0.2f')]);
text(-45,-30,['Vector speed = ' num2str(meandeepsp_h,'%0.1f') ... 
        '; Direction = ' num2str(angleXtoCompass(meandeepdir_h),'%0.0f') '^o']);
text(-45,-38,['Major Axis = ' num2str(majaxdeep_h,'%0.1f') ...
        '; Minor Axis = ' num2str(minaxdeep_h,'%0.1f')]);
text(-45,-45,['Prin. Orient. = ' num2str(majazdeep_h,'%0.0f') '^o']);
l = line([-50 50], [0 0]);
l2 = line([0 0], [-50 50]);
set(l,'color','k', 'linestyle','-.');
set(l2,'color','k', 'linestyle','-.');
set(gca,'xtick',labels);
set(gca,'ytick',labels);
set(gca,'xticklabel','');
set(gca,'yticklabel','');

%Low-Passed, Deep
subplot(2,2,4)
plot(udeep_lp,vdeep_lp,'.')
hold on
plot([0 meanudeep_lp],[0 meanvdeep_lp],'r','linewidth',2);
plot(xxxdeep_lp+meanudeep_lp,yyydeep_lp+meanvdeep_lp,'r','linewidth',2);
axis equal
axis tight
axis([-50 50 -50 50])
labels = [-50 -40 -30 -20 -10 0 10 20 30 40 50];
text(-45,-23,['Mean U = ' num2str(meanudeep_lp,'%0.2f') ...
        '; Mean V = ' num2str(meanvdeep_lp,'%0.2f')]);
text(-45,-30,['Vector speed = ' num2str(meandeepsp_lp,'%0.1f') ... 
        '; Direction = ' num2str(angleXtoCompass(meandeepdir_lp),'%0.0f') '^o']);
text(-45,-38,['Major Axis = ' num2str(majaxdeep_lp,'%0.1f') ...
        '; Minor Axis = ' num2str(minaxdeep_lp,'%0.1f')]);
text(-45,-45,['Prin. Orient. = ' num2str(majazdeep_lp,'%0.0f') '^o']);
l = line([-50 50], [0 0]);
l2 = line([0 0], [-50 50]);
set(l,'color','k', 'linestyle','-.');
set(l2,'color','k', 'linestyle','-.');
set(gca,'xtick',labels);
set(gca,'ytick',labels);
set(gca,'xticklabel','');
set(gca,'yticklabel','');
