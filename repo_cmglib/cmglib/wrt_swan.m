% WRT_SWAN - Reads Matlab grid, writes SWAN grid
% 
% These data are not laid out in conventional order 
slope = 1e-4;
xg = (0:100:10000)';
yg = 100;
zg = xg*slope;

[nx,my]=size(xg)

fid = fopen('fetch.bot','w')
% idla = 1;
% write in default SWAN format (looks like a map when printed out)
for n=(1:my)
  fprintf(fid,' % 7.2f',zg(1:end,n));
  fprintf(fid,'\n ');
end
fclose(fid);

% SWAN input parameters (not conventional MATLAB indexing)
XPC = xg(1,end)
YPC = yg(1,end)
ALPC = 0.0
XLENC = xg(end,1)-xg(1,end)
YLENC = yg(end,1)-yg(1,end)
fprintf(1,'CGRID REGULAR %f %f %f %f %f %d %d\n',...
	XPC,YPC,ALPC,XLENC,YLENC,nx-1,my-1)
if(1),
fid=fopen('fetch.xp','w')
for n=(1:my)
  fprintf(fid,' %g',xg(1:end,n));
  fprintf(fid,'\n ');
end
fclose(fid);

fid=fopen('fetch.yp','w')
for n=(1:my)
  fprintf(fid,'% g',yg(1:end,n));
  fprintf(fid,'\n ');
end
fclose(fid);
end