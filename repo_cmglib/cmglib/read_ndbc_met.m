function m = readndbcmet(metfil)
% read_ndbc_met - Read an NDBC buoy met file
% m = readndbcmet(metfil)
% 
% Input
% metfil	- complete pathname of datafile
% Outputs -- N is the number of points
% dn		(1,N) floats
%		Dates in MATLAB datenum format.
% winddir	(1,N) floats
%		direction wind is from rel to true N
% windspd	(1,N) floats
%		wind speed averaged over 8 minutes, m/s
% gust		(1,N) floats
%		peak wind speed observed over 8 minutes, m/s
% sighgt	(1,N) floats
%		significant wave height, m (highest 1/3 of waves observed
%		in 20 min)
% sigper	(1,N) floats
%		period of peak waves in 20-min interval, s
% avgper	(1,N) floats
%		mean period of waves in 20-min interval, s
% avgdir	(1,N) floats
%		mean direction waves at peak period are from, deg true N
% baro		(1,N) floats
%		atmospheric pressure at surface, bar
% airtemp	(1,N) floats
%		air temperature at surface, C
% h20temp	(1,N) floats
%		water temperature at surface, C 
% dewpt		(1,N) floats
%		dewpoint temperature, C
% vis		(1,N) floats
%		station visibility
%
% Based on readndbcmet.m by Kurt Hanson, USGS St. Pete
% Modified by Chris Sherwood, USGS

inu = fopen(metfil, 'r');
head1 = fgetl(inu)
head2 = fgetl(inu)

[kahuna,count] = fscanf(inu, '%f', [18 inf]);
fclose(inu);
nrec = length(kahuna(1,:));

m.N = nrec;
m.dn = datenum(kahuna(1,:), kahuna(2,:), kahuna(3,:), ...
  kahuna(4,:), kahuna(5,:), zeros(1,nrec))';
m.jd = julian([kahuna(1,:)', kahuna(2,:)', kahuna(3,:)', ...
  kahuna(4,:)', kahuna(5,:)', zeros(1,nrec)']);
m.winddir = kahuna(6,:)';
m.windspd = kahuna(7,:)';
m.gust = kahuna(8,:)';
m.sighgt = kahuna(9,:)';
m.sigper = kahuna(10,:)';
m.avgper = kahuna(11,:)';
m.avgdir = kahuna(12,:)';
m.baro = kahuna(13,:)';
m.airtemp = kahuna(14,:)';
m.h20temp = kahuna(15,:)';
m.dewpt = kahuna(16,:)';
m.vis = kahuna(17,:)';

