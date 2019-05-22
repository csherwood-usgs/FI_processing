function [m,S]=read_nodc291( fn );
% READ_NODC291.M - Read NODC Buoy ASCII data file fn
% [m,S]=read_nodc291( fn );
%
% Metadata are returned in m, time series are returned in S.
% This is an incomplete implementation. Values are returned with units
% described in
% http://www.nodc.noaa.gov/General/NODC-Archive/f291.html

% February 9 - Added record type K
% csherwood@usgs.gov

i=0; % record counter
ii=0; % metadata counter
md_changed=1; % metadata flag

fid = fopen( fn );
while 1,
   s = fgetl(fid);
   if ~ischar(s),break,end
   if s(1:3)~='291',error('291 missing from beginning of line\n'),end
   rec_type = s(10);
   switch lower(rec_type)
      case {'a'}
         i = i+1;
         crn = 0; % reset counter for record type c
         irn = 0; % reset counter for record type i
         krn = 0; % reset counter for record type k
         
         YYYY = str2double(s(4:7));
         MM = str2double(s(19:20));
         DD = str2double(s(21:22));
         HH = str2double(s(23:24));
         MN = str2double(s(25:26));
         SS = 0.;
         S(i).jd = julian([ YYYY MM DD HH MN SS ]); % julian date (GMT)
         S(i).dn = datenum(gregorian(S(i).jd));     % matlab datenum (GMT)

         % temporary copy of metadata
         mt.station =s(11:11+5);
         mt.jd = S(i).jd;
         mt.dn = S(i).dn;
         mt.lat =  (str2double(s(27:28))+str2double(s(29:30))/60 ...
            +str2double(s(31:32))/3600);
         if (s(33)=='S'), mt.lat = -mt.lat;, end
         mt.lon =  (str2double(s(34:36))+str2double(s(37:38))/60  ...
            +str2double(s(39:40))/3600);
         if (s(41)=='W'), mt.lon = -mt.lon;, end
         mt.depth = 0.1*str2double(s(42:46));
         mt.mag_var = str2double(s(47:50));
         mt.buoy_heading = str2double(s(51:53));
         mt.wave_samp_rate = str2double(s(54:57))/10;       % measurements / minute
         mt.wave_samp_duration = str2double(s(58:61))/100;  % minutes
         mt.nf = str2double(s(62:64));
         mt.wind_samp_duration = str2double(s(105:107))/10;
         if( i==1 ),
            md_changed = 1;
         else
            % test for metadata changes here
            md_changed = ~ ( ...
            ( mt.depth == m(ii).depth ) & ...
               ( mt.lon == m(ii).lon ) & ...
               ( mt.lat == m(ii).lat ) ...
            );
         end        
         if( md_changed ),
            ii = ii+1;
            md_changed = 0;
            fprintf(1,'New metadata beginning %s\n',datestr(mt.dn) );
            m(ii)=mt; % copy temp metadata to md
         end
      case {'b'}
         S(i).airt = str2double(s(30:33))/10;  % air temp (deg C)
         S(i).dewp = str2double(s(34:37))/10;  % dew point (deg C)
         S(i).barop = str2double(s(38:42))/10; % atmos pressure at sea level (millibars)
         S(i).wspd = str2double(s(43:46))/100; % wind speed (m/s)
         S(i).wdir = str2double(s(47:50))/10;  % avg wind direction (degrees True)
         S(i).pptn = str2double(s(55:58));     % millimeters
         S(i).solar_rad1 = str2double(s(59:61))/100; % atmospheric radiaition,
         % wavelength < 3.6
         % microns (Langleys / min)
         S(i).solar_rad2 = str2double(s(62:64))/100; % 4 to 50 microns (Langleys / min)
         S(i).hsig = str2double(s(65:67))/10;        % sig wave height (m)
         S(i).tav = str2double(s(68:70))/10;         % avg wave period (s)
         S(i).wvdir = str2double(s(71:73))/10;       % mean wave direc (deg T)
         S(i).sst = str2double(s(80:83))/100;        % sea surface temp (deg C)
         S(i).sal = str2double(s(84:88))/1000;       % salinity (psu)
         S(i).con = str2double(s(89:93))/1000;       % conductivity (millisiemens/cm)
         S(i).tdom = str2double(s(94:96))/10;
         S(i).hmax= str2double(s(97:99))/10;
         S(i).hstp= str2double(s(100:102));
         S(i).gust= str2double(s(103:106))/100;
         S(i).gust_period= str2double(s(107:108));
         S(i).gust2= str2double(s(109:112))/100;
         S(i).gust2_period= str2double(s(113:114));
         S(i).wspd58 = str2double(s(115:117))/10;
         S(i).wdir58 = str2double(s(118:120));
      case {'i'} % directional wave parameter
         ic = str2double(s(27));
         s(strfind(s,' '))='0'; % put zeros in blank spots
         [a]=sscanf(s(28: 27+(ic*(4+4+4+4+4+4+6))),'%04d%04d%04d%04d%04d%04d%06d',[7,ic]);
         S(i).fi(irn+1:irn+ic)=a(1,:)/1e4;  % center of frequency interval (Hz)
         S(i).dfi(irn+1:irn+ic)=a(2,:)/1e4; % resolution of f interval width (Hz)
         S(i).R1(irn+1:irn+ic)=a(3,:)/100; % nondimensional
         S(i).R2(irn+1:irn+ic)=a(4,:)/100; % nondimensional
         S(i).alpha1(irn+1:irn+ic)=a(5,:)/10; % direction in degrees
         S(i).alpha2(irn+1:irn+ic)=a(6,:)/10; % direction in degrees
         S(i).C11i(irn+1:irn+ic)=a(7,:)/1e3; % spectral density (m2/Hz)
         irn=irn+ic;
      case {'k'} % expanded resolution non-directional wave spectra
         ic = str2double(s(34));
         [a]=sscanf(s(35: 34+(ic*(4+4+9))),'%04d%04d%09d',[3,ic]);
         S(i).f(krn+1:krn+ic)=a(1,:)/1e4;  % center frequency (Hz)
         S(i).df(krn+1:krn+ic)=a(2,:)/1e4; % f interval width (Hz)
         S(i).sd(krn+1:krn+ic)=a(3,:)/1e5; % spectral density (m2/Hz)
         krn=krn+ic;
      case {'l'} % expanded resolution co and quad spectra for directional waves data
         % ignore silently
      case {'c'}
         ic = str2double(s(34));
         [a]=sscanf(s(35: 34+(ic*(4+4+6))),'%04d%04d%06d',[3,ic]);
         S(i).f(crn+1:crn+ic)=a(1,:)/1e3;  % center frequency (Hz)
         S(i).df(crn+1:crn+ic)=a(2,:)/1e4; % f interval width (Hz)
         S(i).sd(crn+1:crn+ic)=a(3,:)/1e3; % spectral density (m2/Hz)
         crn=crn+ic;
      case {'d','e','f','h','j'}
         %fprintf(1,'Record type %s not handled yet.\n',rec_type);
      otherwise
         fprintf(1,'Unknown record type: %s\n.',rec_type);
   end
end
