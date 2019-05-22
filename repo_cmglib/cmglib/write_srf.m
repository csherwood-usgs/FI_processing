function iret = write_srf(fn,x,y,z);
% WRITE_SRF - Write an ASCII Surfer file
% iret = write_srf(fn,x,y,z );

% csherwood@usgs.gov 6 Oct 2009
id = 'DSAA';
nx = length(x);
ny = length(y);
xlo = x(1);
xhi = x(end);
ylo = y(1);
yhi = y(end);
zlo = nanmin(z(:));
zhi = nanmax(z(:));

fid = fopen( fn, 'w');
if(fid==-1),iret=-1;,return,end
fprintf(fid,'%s\n','DSAA');
fprintf(fid,'%d %d\n',nx,ny);
fprintf(fid,'%f %f\n',xlo,xhi);
fprintf(fid,'%f %f\n',ylo,yhi);
fprintf(fid,'%f %f\n',zlo,zhi);

% This is a kludge:
% replace NaNs and way-big values with special Surfer missing value
% (If this is not formatted correctly, Surfer won't recognize them as such)
ng = find(z>=1.e35 | isnan(z));
z(ng) = 1.70141e+038;
for j=1:ny,
   fprintf(fid,' %12.5e',z(j,:) );
   fprintf(fid,'\n');
end;
iret=fclose(fid);
