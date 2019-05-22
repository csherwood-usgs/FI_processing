%ADCP_Plotmaker
%m-file to create plots of ADCP data from each station (A-F) for the NY Bight
%Labels/Location of text in general was modified for appearance
%using Adobe Illustrator
%
%Soupy Alexander, 11/18/2001
%Requires NetCDF toolbox, Soupy's tools "singleJD.m" and "gregaxdNM.m", and Rich's
%tool "gregaxd.m"

%Determine the ADCP file corresponding to the station of interest

station = input('Enter the station of interest, use single quotes ');

if station == 'A';
    ADCP_file = netcdf('5951adc-alp.nc','nowrite');
elseif station == 'B';
    ADCP_file = netcdf('5971adc-alp.nc','nowrite');
elseif station == 'C';
    ADCP_file = netcdf('5991adc-alp.nc','nowrite');
elseif station == 'D';
    ADCP_file = netcdf('6011adc-alp.nc','nowrite');
elseif station == 'E';
    ADCP_file = netcdf('6031adc-alp.nc','nowrite');
elseif station == 'F';
    ADCP_file = netcdf('6041adc-alp.nc','nowrite');
end

ADCP_t1 = ADCP_file{'time'}(:);
ADCP_t2 = ADCP_file{'time2'}(:);
ADCP_singleJD = singleJD(ADCP_t1,ADCP_t2);

ADCP_depth = ADCP_file{'depth'}(:);
ADCP_u_1205 = ADCP_file{'u_1205'}(:);
%ADCP_u_1205 = ridfill_nan(ADCP_u_1205);
ADCP_v_1206 = ADCP_file{'v_1206'}(:);
%ADCP_v_1206 = ridfill_nan(ADCP_v_1206);

%Wind Stress from Ambrose Lighthouse
ambrose = netcdf('alsn6-a.nc','nowrite');
wu_ambrose = ambrose{'WU_422'}(:);
wv_ambrose = ambrose{'WV_423'}(:);
time_ambrose = ambrose{'time'}(:);
time2_ambrose = ambrose{'time2'}(:);

jdtime_ambrose = singleJD(time_ambrose,time2_ambrose);

    data_scale_wind = 2;
    plot_scale_wind = 1;

%set boundaries for interesting data for plot

crange = [-50 50];

startjd = julian(1999,12,01,00);
endjd = julian(2000,4,20,00);

subplot(3,1,1)
title(['Low-passed Wind Stress and Current from Station ' station]);
whiskWindNY(jdtime_ambrose,wu_ambrose,wv_ambrose,data_scale_wind,plot_scale_wind);
ylabel('dynes/cm^2')
text(julian(2000,1,1,00)*plot_scale_wind,0,'Estimated 10 m Wind Stress, Ambrose Light');

subplot(3,1,2)
imagesc(ADCP_singleJD,ADCP_depth,ADCP_u_1205',crange);
axis([startjd endjd min(ADCP_depth) max(ADCP_depth)])
gregaxdNM(ADCP_singleJD,5)
ylabel('depth (m)')
title('Eastward Velocity, Low-Pass Filtered')

subplot(3,1,3)
imagesc(ADCP_singleJD,ADCP_depth,ADCP_v_1206',crange);
axis([startjd endjd min(ADCP_depth) max(ADCP_depth)])
gregaxd(ADCP_singleJD,5)
ylabel('depth (m)')
title('Northward Velocity, Low-Pass Filtered')

orient landscape
