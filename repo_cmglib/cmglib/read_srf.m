function srf = read_srf( fn );
% READ_SRF - Read ASCII Surfer file
% srf = read_srf( fn );

% csherwood@usgs.gov 7 July 2006
srf.id = 'XXXX';
srf.nx = -1;
srf.ny = -1;
srf.xlo = NaN;
srf.xhi = NaN;
srf.ylo = NaN;
srf.yhi = NaN;
srf.zlo = NaN;
srf.zhi = NaN;
srf.x = NaN;
srf.y = NaN;
srf.z = NaN;
srf.dx = NaN;
srf.dy = NaN;

fid = fopen( fn, 'r');
if(fid < 1),
   es=sprintf('READ_SRF can''t open %s\n',fn);
   error(es);
end
srf.id = fgetl(fid);
if( ~strcmp( srf.id, 'DSAA') )
   es= sprintf('READ_SRF wants DSAA in first line; got %s\n',srf.id);
   error(es)
end
tline = fgetl(fid);
[A,COUNT] = sscanf(tline,'%d');
srf.nx = A(1);
srf.ny = A(2);
tline = fgetl(fid);
[A,COUNT] = sscanf(tline,'%f');
srf.xlo = A(1);
srf.xhi = A(2);
tline = fgetl(fid);
[A,COUNT] = sscanf(tline,'%f');
srf.ylo = A(1);
srf.yhi = A(2);
tline = fgetl(fid);
[A,COUNT] = sscanf(tline,'%f');
srf.zlo = A(1);
srf.zhi = A(2);
srf.z=fscanf(fid,'%f ',[srf.nx,srf.ny])';
fclose(fid);
srf.dx = (srf.xhi-srf.xlo)/(srf.nx-1);
srf.dy = (srf.yhi-srf.ylo)/(srf.ny-1);
srf.x = srf.xlo:srf.dx:srf.xhi;
srf.y = srf.ylo:srf.dy:srf.yhi;
%pcolor(x,y,srf.z)