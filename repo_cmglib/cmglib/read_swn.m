function a = read_swn( nx, ny, fn, idla )
% READ_SWN - Read SURFER SWAN file
if(exist('idla')~=1),idla=1;,end
fid=fopen(fn,'r');
if(~fid)
  fprintf(1,'Problem reading file: %s\n',fn)
  error('in read_swn');
end
a = NaN*ones(1,nx*ny);
if((idla==1)|(idla==2)),
  i = 1;
  while(~feof(fid)),
    s = fgetl(fid);
    [A,c]=sscanf(s,'%f ');
    a(1,i:i+c-1)=A(:)';
    i=i+c;
  end
  a = reshape(a,nx,ny);
else
  error('read_swn: idla = 3 or 4 not implemented')
end