function PrT_tides(varargin)

% PrT_tides.m
%
% plot time series of tidal elevations for specified site and time interval
%
% Syntax
%
%   PrT_tides
%   PrT_tides(days)
%   PrT_tides(startdate,enddate)
%   PrT_tides(switches)
%   PrT_tides(days,switches)
%   PrT_tides(startdate,enddate,switches)
%
% Description
%
%   days = scalar, number of days to plot (default = 2)
%   startdate, enddate = date strings (mm-dd-yyyy)
%   *** days or startdate,enddate must be listed before switches***
%
%   switches =
%     'moonoff' -- omits moon phase from plots (on by default)
%     'dayntoff' -- omits day/night shading from plot background (on by
%     default)
%     'hthi' -- put high tide text on plot -  puts date, time, elevation of
%       4 or less roughly evenly spaced high tides - labels above tidal
%       elevation plot
%     'htlo' -- same as high but puts labels on bottom of plot
%     'station' -- defines station to plot, followed by station name string
%     (default = Mattapoisett, MA), ex. PrT_tides('station','boston')
%     'lonlat' -- find station nearest specified longitude and
%     latitude, negative for west longitude, ex.
%     PrT_tides('lonlat',-70.6,41.6) returns tides for Woods Hole, MA
% 
% Examples
%
% PrT_tides(20)
% PrT_tides(20,'moonoff','station','Boston')
% PrT_tides('station','los angeles','hthi')
% PrT_tides('3-10-2010','4-10-2010','dayntoff')
%
% Dependency
%
% must have t_tides and air_sea toolboxes installed, available at: 
% http://www.eos.ubc.ca/~rich/
%
% ****** Note******
% t_xtides uses an old xtides database. this function does not work for all
% stations supported in the current version of xtides.
%
% ****** Note *********
% uses local timebase by correcting default (GMT) with output field
% .timezone. On east coast US, time is EST. Does not correct for daylight
% savings.
%
% by Patrick Dickhudt  31-Mar-2010

fsz = 14;  % fontsize
set(0,'defaultaxesfontsize',fsz)
set(0,'defaulttextfontsize',fsz)

figure
fp = get(gcf,'position');
set(gcf,'position',[fp(1) fp(2) 800 400])

% determine number of days to plot
daz = 2;
switch nargin
    case 0
    case 1
        if isnumeric(varargin{1})
        daz = varargin{:};
        end
    otherwise
        if isnumeric(varargin{1})
            daz = varargin{1};
        else 
            try
                daz = datenum(varargin{1}):1/96:datenum(varargin{2}); 
            end
        end
end

% set station
if any(strcmpi('station',varargin))
    staind = find(strcmpi('station',varargin))+1;
    station = varargin{staind};       
elseif any(strcmpi('lonlat',varargin))
    llinds = find(strcmpi('lonlat',varargin))+[1 2];
    stlon = varargin{llinds(1)};
    stlat = varargin{llinds(2)};
    stout = t_xtide(stlon,stlat,'format','info');
    station = stout.station;
    comma = find(station == ',');
    station = station(1:comma);
else
        station = 'Mattapoisett';
end

% set units
units = 'meters';
ylab = 'Water Level';
MLLW = 'MLLW';
yout = t_xtide(station,daz,'units',units,'format','info');

if strcmpi('knots',yout.units)
    units = 'm/s';
    ylab = 'Current Velocity';
    MLLW = '';
end

% make plot with t_xtide and get t_xtide data structs
t_xtide(station,daz,'units',units)
hold on
yout = t_xtide(station,daz,'units',units,'format','info');
yout2 = t_xtide(station,daz,'units',units,'format','times');

lat = yout.latitude;
lon = -yout.longitude;

% get limits from data
xl = [floor(yout2.mtime(1)) ceil(yout2.mtime(end))];
dlim(1) = min(yout2.value);
dlim(2) = max(yout2.value);

if(any(strcmpi('htlo',varargin)))
yl(1) = floor((dlim(1) - .25*diff(dlim))*2)/2;
yl(2) = ceil((dlim(2) + .15*diff(dlim))*2)/2;    
elseif(any(strcmpi('hthi',varargin)))  
yl(1) = floor((dlim(1) - .1*diff(dlim))*2)/2;
yl(2) = ceil((dlim(2) + .25*diff(dlim))*2)/2;    
else
yl(1) = floor((dlim(1) - .1*diff(dlim))*2)/2;
yl(2) = ceil((dlim(2) + .15*diff(dlim))*2)/2;
end
dts = xl(1):1:xl(2);

% calculate time of sunrise and sunset and shade plot gray when dark
if(~any(strcmpi('dayntoff',varargin)))
for i = 1:length(dts)
    
dv = datevec(dts(i));

[rhr,rmin,shr,smin] = sunrise(dv(2),dv(3),dv(1),lat,lon);

rhrst = rhr+yout.timezone;
shrst = shr+yout.timezone;

if rhrst >= 24
rhrst = rhrst - 24;
elseif rhrst <= 0
rhrst = rhrst + 24;
end

if shrst >= 24
shrst = shrst - 24;
elseif shrst <= 0
shrst = shrst + 24;
end

snris = datenum([dv(1:3) rhrst rmin 0]);
snset = datenum([dv(1:3) shrst smin 0]);

rx1 = datenum(dv);
rw1 = (snris-rx1);
rectangle('position',[rx1,yl(1),rw1,diff(yl)],'facecolor',[.85 .85 .85],'edgecolor','none')

rx2 = datenum(dv)+1;
rw2 = (rx2-snset);
rectangle('position',[snset,yl(1),rw2,diff(yl)],'facecolor',[.85 .85 .85],'edgecolor','none')
end
end

% format plot
t_xtide(station,daz,'units',units) % plot again to bring to front
hold on
plot(xl,[0 0],'k-','linewidth',1.25)
ylabel(['           ' ylab ' (' yout.units ')'])
title({yout.station; [datestr(xl(1),1) '  -  ' datestr(xl(2),1)]})
fo1 = findobj(gcf,'type','line');
set(fo1,'linewidth',2)
text(xl(1)-diff(xl).*.04,0,MLLW,'horizontalalignment','right')
xlim(xl)
ylim(yl)

% set x ticks and x tick labels
xt = xl(1):diff(xl)/4:xl(end);
set(gca,'xtick',xt);
set(gca,'tickdir','out');

if diff(xl) <20
   datetick('x','dd-mmm HH:MM:SS','keepticks','keeplimits')
else
   datetick('x',1,'keepticks','keeplimits')
end

% put high tide times on plot above water level line
if(any(strcmpi('hthi',varargin)))   
    hti = find(yout2.type);  % high tide indices
    htm = yout2.mtime(hti);  % time of high tide
    hth = yout2.value(hti);  % height of high tides
    hty = hth + diff(yl).*.15;  % text location = height of high tide + offset
    htx = xl(1):diff(xl)/5:xl(end);
    htx = htx(2:end-1);  % 4 evenly spaced times withing plot window
    hti2 = nearest(htm,htx);

for i3 = 1:length(hti2)
    ind = hti2(i3);
    text(htm(ind),hty(ind),{datestr(htm(ind),1);datestr(htm(ind),13);[num2str(hth(ind),3) ' m']},'horizontalalignment','center')
end
end      

% put high tide times on bottom of plot
if(any(strcmpi('htlo',varargin)))   
    hti = find(yout2.type);  % high tide indices
    htm = yout2.mtime(hti);  % time of high tide
    hth = yout2.value(hti);  % height of high tides
    hty = yl(1) + diff(yl)*.02;  % bottom of plot + offset
    htx = xl(1):diff(xl)/5:xl(end);
    htx = htx(2:end-1);  % 4 evenly spaced times withing plot window
    hti2 = nearest(htm,htx);

for i3 = 1:length(hti2)
    ind = hti2(i3);
    text(htm(ind),hty,{datestr(htm(ind),1);datestr(htm(ind),13);[num2str(hth(ind),3) ' m']},'horizontalalignment','center','verticalalignment','bottom')
end
end      

% calculate moon phase data and some plot parameters for drawing moons
if(~any(strcmpi('moonoff',varargin)))
mp = moonphase(dts);

mpd = [mp.nw mp.fq mp.fl mp.lq];  % make vector of dates with new, 1st 1/4, full, or last 1/4 moon
mps = cellstr([repmat('nw',length(mp.nw),1);repmat('fq',length(mp.fq),1);repmat('fl',length(mp.fl),1);repmat('lq',length(mp.lq),1)]);
gp = get(gcf,'position');
asz = get(gca,'position');
whrat = gp(3)./asz(3)./(gp(4)./asz(4));  % calc ratio of x-axis to y-axis length - use to correct so moons are round

rdx = diff(xl)./25./whrat;
rdy = diff(yl)./25;

t = linspace(0,2*pi);
t2 = linspace(0,pi);
t3 = linspace(pi,2*pi);

% draw moons on plot
for i2 = 1:length(mpd)
    
switch mps{i2}
    case 'nw'
        patch(mpd(i2)+rdx*sin(t),yl(2)-2*rdy+rdy*cos(t),'k-','facecolor',[.2 .2 .2])
    case 'fq'
        patch(mpd(i2)+rdx*sin(t),yl(2)-2*rdy+rdy*cos(t),'k-','facecolor',[.95 .95 .95])
        patch(mpd(i2)+rdx*sin(t3),yl(2)-2*rdy+rdy*cos(t3),'k-','facecolor',[.2 .2 .2])
    case 'fl'
        patch(mpd(i2)+rdx*sin(t),yl(2)-2*rdy+rdy*cos(t),'k-','facecolor',[.95 .95 .95])
    case 'lq'
        patch(mpd(i2)+rdx*sin(t),yl(2)-2*rdy+rdy*cos(t),'k-','facecolor',[.95 .95 .95])
        patch(mpd(i2)+rdx*sin(t2),yl(2)-2*rdy+rdy*cos(t2),'k-','facecolor',[.2 .2 .2])
end
end
end
%% moonphase function

function mp = moonphase(dts)
    
% determine days of new, first quarter, full, and last quarter moon for
% time period specified by dts
%
% Inputs
% dts = scalar or vector of matlab datenums
%
% Outputs
% mp.nw = day(s) with new moon (datenum)
% mp.fq = day(s) with first quarter moon (datenum)
% mp.fl = day(s) with full moon (datenum)
% mp.lq = day(s) with last quarter moon (datenum)

mnp = 29.530588853;  % period of moons orbit in relation to tides (or something like that)

dtsl = dts(1):.01:dts(end);
dtj = julian(datevec(dtsl));
ip =(dtj-2451550.1)./mnp - fix((dtj-2451550.1)./mnp);

nw = abs(ip) < .005/mnp | abs(ip) > 1-.005/mnp;
fq = abs(ip-.25) < .005/mnp;
fl = abs(ip-.5) < .005/mnp;
lq = abs(ip-.75) < .005/mnp;

mp.nw = dtsl(nw);
mp.fq = dtsl(fq);
mp.fl = dtsl(fl);
mp.lq = dtsl(lq);

end

%% sunrise function

function [rhr,rmin,shr,smin]=sunrise(mon,da,yr,lat,long)

% NOTE: this is included in air_sea but didn't work for me in that form.
% This is a slightly modified version that seems to work. Commented below.
% SUNRISE: computes sunrise and sunset times for specified day and
% location. 
% [rhr,rmin,shr,smin] = SUNRISE(mon,da,yr,lat,long) computes the time 
% of sunrise rhr:rmin and sunset shr:smin to the nearest minute in GMT 
% for a calendar day(s) and a specified (scalar) position.   
% 
% INPUT:  mon - month (e.g., Jan is 1)
%         da - day (e.g., Jan 10th is 10)
%         yr - year (e.g., 1995)
%         lat - latitude [deg]
%         long - longitude (west is positive)  [deg] 
%
% OUTPUT: rhr,rmin  - sunrise in GMT hours and minutes
%          shr,smin  - sunset in GMT hours and minutes.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8/28/98: version 1.1 (contributed by RP)
% 8/5/99: version 2.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert calender time to julian yd
j=julianmd(yr,mon,da,0);
j0=julianmd(yr,1,1,0);
yd=j(:)-j0(:);

% compute solar altitude for entire day
dt=1./2880;

% we don't want abs(long)>180...
if long<-180, long=long+360; end;
if long>180,  long=long-360; end;

time=dt.*[0:2879]'+long/360; % have a whole day, beginning at midnight (near enough)
yday=yd(ones(1,2880),:)+time(:,ones(length(yd),1));

if length(yr)>1,
  yr=yr(:,ones(1,2880))';
end;

[z,sorad]=soradna1(yday(:),yr(:),long,lat);

z=reshape(z,2880,length(yd));
sorad=reshape(sorad,2880,length(yd));

[ir,jr]=find(sorad(1:2879,:)==0 & sorad(2:2880,:)>0);
[is,js]=find(sorad(2:2880,:)==0 & sorad(1:2879,:)>0);



srise=zeros(length(yd),1);
sset=any(sorad>0);

srise(jr)=yday(ir+(jr-1)*2880);

% sset(js) =yday(is+(js-1)*2880);   % original code - outputs 0 for shr and
%                                       smin
sset =yday(is+(js-1)*2880);  % modifed by Patrick Dickhudt

rhr=fix(rem(srise,1)*24);
rmin=rem(rem(srise,1)*1440,60);
shr=fix(rem(sset,1)*24);
smin=rem(rem(sset,1)*1440,60);

end

%% nearest function

function n=nearest(jd,jdo)
%function n=nearest(jd,jdo)
%
% purpose: return element numbers in jd of specified times jdo
%
% inputs: jd=julian day row vector at which element number is unknown
%         jdo=julian day row vector with desired times
%
% outputs: n=element numbers of nearest times in array jd
%
% by David Schoellhamer 11/3/95
% modified slightly by Patrick Dickhudt 2-Apr-2010 (removed error warnings)

n=nan*ones(size(jdo));
njd=length(jd);
njdo=length(jdo);
for i=1:njdo
  if i>1
    zi=n(i-1):njd;
  else
    zi=1:njd;
  end
 
  [z,zn]=min(abs(jd(zi)-jdo(i)));
  n(i)=zn+zi(1)-1;
end
end
end