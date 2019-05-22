% UBR_FROM_BUOY - Script to calculate ubr and tr from bouy formats
if(1),
    InDir = 'g:\\data\\46025\\';
    %InDir = 'c:\\crs\\data\\46025\\';
    InBuoy = '46025'
    InYY = '2002'
    InMO = '10_12'
    InExt = '.txt';
    fn = [InDir,InBuoy,'_',InYY,InMO,InExt];
    [m,S]=read_nodc291( fn );
    save nov_2002 S m
else
    load nov_2002
end

ubdepth = 60;
for i=1:length([S.hsig]),
  if(S(i).hsig > 0.1),
    ubs=ubspec( ubdepth, S(i).sd, S(i).f,'default',S(i).df );
    [ubpw,Tpw]=ubspecfun( S(i).hsig, S(i).tdom, ubdepth );
    S(i).ubr = ubs.ubr;
    S(i).tr = ubs.Tr;
    S(i).ubpw = ubpw;
    S(i).tpw = Tpw;
  else
    S(i).ubr = NaN;
    S(i).tr = NaN;
    S(i).ubpw = NaN;
    S(i).tpw = NaN;
  end
end
