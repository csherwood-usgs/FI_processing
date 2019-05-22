function zone = utmzone(lon)
% UTMZONE - Returns zone for input longitude
% (West longitudes should be negative; does not handle weird regions
% around Norway and Svalbard...see 
% http://www.fmnh.helsinki.fi/english/botany/afe/map/utm.htm)
% UTM zones are 6 degrees wide, numbered from 1 (180W to 174W) to 60
% (174E to 180E)
if( any(lon>180) | any(lon<-180) ),error('lon must be -180 <= lon <= 180'),end
zone = min(60,floor((lon+180+6+eps)/6));
