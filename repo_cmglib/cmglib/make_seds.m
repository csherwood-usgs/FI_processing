% MAKE_SEDS - Calculate physical parameters for grain size distribution
% This writes a file called seds.txt, which can be edited and read in
% by READ_SEDS.m

% Chris Sherwood, USGS
% Last revised April 27, 2005

% Change these two lines to specify phi sizes and fractions
phi = [0.5 3.5 ];      % sediment sizes (phi units)
fr =  [.85 .15  ];        % fractions
ns = length(phi);
% fr is not important in ROMS...the relative amount of each size is
% determined in make_init
%fr = ones(ns,1) ./ (ns) ;
%fr = ones(ns,1); % placeholder for fractional amount per size class in bed

% Calls dset, dmm, taucrit, visc_tsp (all in CMGLib) and CSIRO sw_dens
Temp = 15;
Sal = 34;
Press = 13;
rhow = sw_dens( Sal, Temp, Press );      % water density (kg/m^3)
mu = visc_tsp( Temp, Sal, Press )/1000;  % kg m-1 s-1
nu = mu/rhow;                            % kinematic viscocity (m^2/s)

Dcm = 0.1*dmm(phi); 
rhos = 2650*ones(ns,1);
ws = zeros(ns,1);
tc = zeros(ns,1);
for j=1:ns,
    % I/O for these functions are in CGS; ws is negative
    ws(j) = dset( 0.8, 5, Dcm(j), rhos(j)/1000, rhow/1000, nu*10000 )./100;
    tc(j) = taucrit( rhos(j)/1000, rhow/1000, nu*10000, Dcm(j) )/10;
end

fid = fopen('seds.txt','w');
  fprintf(fid,'% 7.2f\n',Temp);
  fprintf(fid,'% 7.2f\n',Sal);
  fprintf(fid,'% 7.2f\n',Press);
  fprintf(fid,'% 7.2f\n',rhow);
  fprintf(fid,'% 14.6e\n',nu);
  fprintf(fid,'% 3d\n',ns);
for j=1:ns,
  fprintf(fid,'% 5.3f % 5.2f % 8.5f % 6.1f % 10.7f % 9.6f\n',...
	  fr(j),phi(j),10*Dcm(j),rhos(j),ws(j),tc(j));
end
fclose(fid);
  
  
  
  