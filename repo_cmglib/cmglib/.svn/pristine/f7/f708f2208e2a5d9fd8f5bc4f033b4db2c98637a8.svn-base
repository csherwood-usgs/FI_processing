function isok=buoyftp(sta,yy,mm);
% BUOYFTP - Get buoy data for a month in ASCII F291 Format
%
% Input:
%    sta = Alpha station name (e.g., '46025')
%    yy  = Numeric year (e.g., 2004)
%    mm  = Numeric month (e.g., 10)
%
% These data come in compressed (.Z) format.
% http://www.nodc.noaa.gov/BUOY/46025.html
% csherwood@usgs.gov January 31, 2006
isok = 0;
fsn = 'ftp.nodc.noaa.gov';
fdn1 = '/pub/f291/'
fdn2 = sprintf('%04d%02d',yy,mm)
fdn = [fdn1,fdn2];
ffn = sprintf('%s_%04d%02d.Z',sta,yy,mm)
   
f = ftp(fsn,'anonymous','csherwood@usgs.gov');
binary(f);
cd(f,fdn);
mget(f,ffn);
close(f);
isok =1; % actually, no tests done

