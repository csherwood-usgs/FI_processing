function ss = read_seds( fn );
% READ_SEDS - Reads ascii file, returns structure of sediment props

if(nargin < 1),
  fn = 'seds.txt';
end
fid = fopen(fn,'r');
inline = fgetl(fid);
Temp=sscanf(inline,'%f');
inline = fgetl(fid);
Sal=sscanf(inline,'%f');
inline = fgetl(fid);
Press=sscanf(inline,'%f');
inline = fgetl(fid);
rhow=sscanf(inline,'%f');
inline = fgetl(fid);
nu=sscanf(inline,'%f');
inline = fgetl(fid);
ns = sscanf(inline,'%d');
fr=zeros(ns,1);
phi=zeros(ns,1);
Dmm=zeros(ns,1);
rhos=zeros(ns,1);
ws=zeros(ns,1);
tc=zeros(ns,1);
for j=1:ns,
  inline = fgetl(fid);
  a = sscanf(inline,'%f %f %f %f %f %f');
  fr(j)=a(1);
  phi(j)=a(2);
  Dmm(j)=a(3);
  rhos(j)=a(4);
  ws(j)=a(5);
  tc(j)=a(6);
end
fclose(fid);

% get rid of zero-fraction sizes
i = find(fr>eps);
ns = length(i);
fr = fr(i);
phi = phi(i);
Dmm = Dmm(i);
rhos = rhos(i);
ws = ws(i);
tc = tc(i);

% sort by size
if( ns > 1 ),
  [junk,idx]=sort(phi(:)); % idx is index of ascending phi
  idx = flipud(idx); % but we want descending phi
  fr = fr(idx);
  phi = phi(idx);
  Dmm = Dmm(idx);
  rhos = rhos(idx);
  ws = ws(idx);
  tc = tc(idx);
end

fprintf(1,'Temperature         :% 7.2f (deg C)\n',Temp);
fprintf(1,'Salinity            :% 7.2f (psu)\n',Sal);
fprintf(1,'Pressure (depth)    :% 7.2f (decibars)\n',Press);
fprintf(1,'Water density (rhow):% 7.2f (kg/m3)\n',rhow);
fprintf(1,'Kin. viscosity (nu) :% 14.6e (m2/s)\n',nu);
fprintf(1,'ns                  :% 3d sed classes\n',ns);


% check fraction math
if( abs(sum(fr(:))-1)>1e-3 ),
  fprintf(1,'Fractions add up to %g. Renormalizing to 1\n',sum(fr));
  fr = fr ./ sum(fr);
end

fprintf(1,' frac  Size    Size    rhos       ws       tcrit\n') 
fprintf(1,'       (phi)   (mm)   (kg/m3)    (m/s)     (N/m2)\n')
for j=1:ns,
  fprintf(1,'% 5.3f % 5.2f % 8.5f % 6.1f % 10.7f % 9.6f\n',...
	  fr(j),phi(j),Dmm(j),rhos(j),ws(j),tc(j));
end

if(ns>1),
  % the tails are from van Rijn
  % calculate d10, d50, d90
  fr = fr+eps; % avoid probs in interp1
  hulp = cumsum(fr);
  if(hulp(1) > 0.1),
    fine_tail = max( 7, phi(1)+1 )
    d10 = fine_tail-0.1/hulp(1)*phi(1);
  else
    d10 = interp1(hulp,phi,.1,'linear',NaN);
  end
  d50 = interp1(hulp,phi,.5,'linear','extrap');
  d90 = interp1(hulp,phi,.9,'linear',NaN);
  if(isnan(d90)), d90 = phi(ns)+.5;, end 
else
  d10 = phi-.5;
  d50 = phi;
  d90 = phi+.5;
end
fprintf(1,'d10: %f\nd50: %f\nd90: %f\n',d10,d50,d90);

ss.T = Temp;
ss.S = Sal;
ss.P = Press;
ss.rhow = rhow;
ss.nu = nu;
ss.ns = ns;
ss.d10 = d10;
ss.d50 = d50;
ss.d90 = d90;
ss.fr = fr;
ss.phi = phi;
ss.Dmm = Dmm;
ss.rhos = rhos;
ss.ws = ws;
ss.tc = tc;
% Interpolate to find ws and tc for d50 and d90
ss.d50_ws = interp1( phi,ws,d50);
ss.d90_ws = interp1( phi,ws,d90);
ss.d50_tc = interp1( phi,tc,d50);
ss.d90_tc = interp1( phi,tc,d90);
fprintf(1,'d50 ws: %f d90 ws: %f\n',ss.d50_ws,ss.d90_ws)
fprintf(1,'d50 tc: %f d90 tc: %f\n',ss.d50_tc,ss.d90_tc)