function getWW3SpecPars(dates,bound,utm)

%getWW3SpecPars    Extracts data from WW3 for SWAN input.
%   getWW3SpecPars(DATES,BOUND,UTM) extracts WW3 data from the *.yyyymm.grb
%   (the yyyymm comes from DATES) in pwd for the DATE(S) (yyyymmddHH) along
%   the segment specified by BOUND (lon,lat or x,y). Set UTM=1 if bound is 
%   given in UTM. The BOUNDSPEC lines for the SWAN input file are written 
%   to pars_DATE.bnd or pars_DATES(1)_DATES(2).bnd.  If lenght(DATES)>1 the 
%   TPAR files are written to tpar_num_DATES(1)_DATES(2).dat.
%
%   NOTE: this code is for model_id.par.yyyymm.grb files. So for example, 
%   in the pwd you should have the nah.dp.200409.grb, nah.hs.200409.grb and
%   nah.tp.200409.grb files for this to work.
%
%   getWW3SpecPars(dates,bound,utm)
%
% 2/05 DMT
% 6/05 DMT - updated for getting data for time dependent runs
% 7/05 DMT - cut out all pre-Matlab wgrib stuff and do it here.

if ~exist('utm','var')
   utm = 0;
end

% get the grib file names from the date(s) provided
fnames = dir(['*.',datestr(datenum(num2str(dates(1)),'yyyymmddHH'),'yyyymm'),'.grb']);

% round down/up the dates to the nearest 3 hr mark provided by the WW3 data
dates = dates(:);
if rem(rem(dates(1),100),3)~=0
   dates(1) = dates(1) - rem(rem(dates(1),100),3);
end
if length(dates)>1
   if rem(rem(dates(2),100),3)~=0
      dates(2) = dates(2) + (3-rem(rem(dates(2),100),3));
   end
end
dates = datenum(num2str(dates),'yyyymmddHH');
% get the grib record numbers for the date(s) provided
days = str2num(datestr(dates,'DD'));
hrs = str2num(datestr(dates,'HH'));
rns = (days*8 - 7) + (hrs/3 + 1) - 1;
if length(dates)>1
   rns = rns(1):1:rns(2);
end

% get axes of WW3 data
tmp = evalc(['!wgrib ',fnames(1).name,' -d 1 -V']);
lat = char(regexp(tmp,'lat .* nxny','match'));
lat = str2num(char(regexp(lat,'-*\d+\.\d+','match')));
yax = lat(2):lat(3):lat(1);
lon = char(regexp(tmp,'long .* (','match'));
lon = str2num(char(regexp(lon,'-*\d+\.\d+','match')));
xax = (lon(1):lon(3):lon(2))-360;

D = struct('H',[],'T',[],'D',[]);
more off
disp(' gribbing ww3 files now...')
for zz=1:length(rns)
   % read in the data
   % data = [x,y,par]; par order should be (dp,hs,tp)
   for yy=1:3
      eval(['!wgrib ',fnames(yy).name,' -d ',num2str(rns(zz)),...
         ' -text -o dump.tmp'])
      fid = fopen('dump.tmp','r');
      np = fscanf(fid,'%i',[1 2]);  % this should be the same through the set
      tmp = fscanf(fid,'%f',[np(1) np(2)]);
      fclose(fid);
      !rm dump.tmp
      tmp = flipdim(permute(tmp,[2 1 3]),1);
      tmp(tmp==9.999e20) = nan; % 9.999e20 = WW3 flag for no data
      foo = char(regexp(fnames(yy).name,'\.[htd]','match'));
      if find('htd'==foo(2))==1
         D.H(:,:,zz) = tmp;
      elseif find('htd'==foo(2))==2
         D.T(:,:,zz) = tmp;
      else
         D.D(:,:,zz) = tmp;
      end
   end
end
disp(' done gribbing.')
% WW3 direction is nautical convention: direction FROM measured CW from
%  North
% I'm using SWAN cartesian convention: direction TOWARD measured CCW
%  from East
% So, fix it
D.D = (270-D.D<0)*360 + (270-D.D);

% trim data to region of interest using bound. Need to do this in lon/lat,
% so if you want to work in utm convert bound to lon/lat for this trimming.
% +/-4 degrees seems good...
if utm==1
   zone = input(' So, you want to work in UTM. Which zone? : ','s');
   fid = fopen('utm.tmp','w');
   fprintf(fid,'%.2f %.2f\n',bound');
   fclose(fid);
   eval(['!proj -I -f %.4f +ellps=GRS80 +proj=utm +zone=',zone,...
      '+units=m utm.tmp > ll.tmp'])
   fid = fopen('ll.tmp','r');
   bnd = fscanf(fid,'%f',[2 prod(size(bound))]);
   fclose(fid);
   bnd = bnd';
   ! rm ll.tmp utm.tmp
else
   bnd = bound;
end
xid = find(xax>=min(bnd(:,1))-4 & xax<=max(bnd(:,1))+4);
yid = find(yax>=min(bnd(:,2))-4 & yax<=max(bnd(:,2))+4);
xax = xax(xid);
yax = yax(yid);
D.H = D.H(yid,xid,:);
D.T = D.T(yid,xid,:);
D.D = D.D(yid,xid,:);
[xg,yg] = meshgrid(xax,yax);
   
% transform grid to utm if desired
if utm==1
   fid = fopen('ll.tmp','w');
   fprintf(fid,'%.2f %.2f\n',[xg(:)';yg(:)']);
   fclose(fid);
   eval(['!proj -f %.4f +ellps=GRS80 +proj=utm +zone=',zone,...
      '+units=m ll.tmp > utm.tmp'])
   fid = fopen('utm.tmp','r');
   tmp = fscanf(fid,'%f',[2 prod(size(xg))]);
   fclose(fid);
   ! rm ll.tmp utm.tmp
   xg = reshape(tmp(1,:),size(xg));
   yg = reshape(tmp(2,:),size(yg));
end

% find indices of xg & yg (WW3-points) that are shortest distance from
%  each [xbound ybound] point
xbound = bound(:,1);
ybound = bound(:,2);
for ii = 1:size(bound,1)
   r  = sqrt((xg-xbound(ii)).^2 + (yg-ybound(ii)).^2);
   [row(ii),col(ii)] = find(min(min(r))==r);
end

% plot what we've got so far (just the H from the first date)
figure
pcolor(xg,yg,D.H(:,:,1))
axis equal
axis tight
hold on

% difference row and column indices... will be used for determining whether
%  boundary is north-south or east-west and how many points inbetween
coldif = diff(col);
rowdif = diff(row);

% THE loop
data = struct('seg',[],'par',[]);
cdist = 0;
dist = 0;
for zz=1:length(rns)
   dat = [];
   for ii=1:length(xbound)-1      % step through all [xbound,ybound] points
      if coldif(ii)==0              % moving north or south along boundary
         if rowdif(ii)>0            % moving north
            mm = ii;
            nn = ii+1;
         else                       % moving south
            mm = ii+1;
            nn = ii;
         end
         id = row(mm):row(nn);
         res = mean(diff(yg(id,col(ii))));
         num = ceil((ybound(nn)-ybound(mm))/res);
         res = (ybound(nn)-ybound(mm))/num;
         yb = [ybound(mm):res:ybound(nn)];   % x boundary points
         xb = ones(length(yb),1)*xbound(ii); % y boundary points
         % interp WW3 to the boundary points
         hb = griddata(xg,yg,D.H(:,:,zz),xb,yb);
         hb = hb(:,1);
         tb = griddata(xg,yg,D.T(:,:,zz),xb,yb);
         tb = tb(:,1);
         db = griddata(xg,yg,D.D(:,:,zz),xb,yb);
         db = db(:,1);
         % nans may result if boundary points are close to land and the WW3
         % resolution is course.  Fill these nans with the closest WW3 data
         % point
         for jj=1:length(hb)
            if isnan(hb(jj))
               r = sqrt((xg(id,col(ii))-xb(jj)).^2 +...
                  (yg(id,col(ii))-yb(jj)).^2);
               hb(jj) = D.H(id(find(min(r))),col(ii),zz);
               tb(jj) = D.T(id(find(min(r))),col(ii),zz);
               db(jj) = D.D(id(find(min(r))),col(ii),zz);
            end
         end
         % get the orientation of things correct
         if rowdif(ii)>0
         else
            yb = fliplr(yb);
            hb = flipud(hb); tb = flipud(tb); db = flipud(db);
         end
         yb = yb(:);
      else                 % moving east or west along boundary
         if coldif(ii)>0   % moving east
            mm = ii;
            nn = ii+1;
         else
            mm = ii+1;
            nn = ii;       % moving west
         end
         id = col(mm):col(nn);
         res = mean(diff(xg(row(ii),id)));
         num = ceil((xbound(nn)-xbound(mm))/res);
         res = (xbound(nn)-xbound(mm))/num;
         xb = [xbound(mm):res:xbound(nn)];   % x boundary points
         yb = ones(length(xb),1)*ybound(ii); % y boundary points
         % interp WW3 to the boundary points
         hb = griddata(xg,yg,D.H(:,:,zz),xb,yb);
         hb = hb(1,:);
         tb = griddata(xg,yg,D.T(:,:,zz),xb,yb);
         tb = tb(1,:);
         db = griddata(xg,yg,D.D(:,:,zz),xb,yb);
         db = db(1,:);
         % nans may result if boundary points are close to land and the WW3
         % resolution is course.  Fill these nans with the closest WW3 data
         % point
         for jj=1:length(hb)
            if isnan(hb(jj))
               r = sqrt((xg(row(ii),id)-xb(jj)).^2 +...
                  (yg(row(ii),id)-yb(jj)).^2);
               hb(jj) = D.H(row(ii),id(find(min(r))),zz);
               tb(jj) = D.T(row(ii),id(find(min(r))),zz);
               db(jj) = D.D(row(ii),id(find(min(r))),zz);
            end
         end
         % get the orientation of things correct
         if coldif(ii)>0
         else
            xb = fliplr(xb);
            hb = fliplr(hb); tb = fliplr(tb); db = fliplr(db);
         end
         xb = xb(:); hb = hb(:); tb = tb(:); db = db(:);
      end
      % trim repeat points
      if ii~=1
         xb = xb(2:end); yb = yb(2:end); hb = hb(2:end); tb = tb(2:end);
         db = db(2:end);
         cdist = dist+res;
      end
      % get cumulative distance
      if zz==1
      for jj=1:length(xb)-1
         cdist = [cdist;cdist(length(cdist))+res];
         dist = cdist(length(cdist));
      end
      data.seg = [data.seg;xb,yb,cdist];
      end
      % pack in the data
      dat = [dat;hb,tb,db];
      if ii==length(xbound)-1
         data.par(:,:,zz) = dat;
      end
   end
end

plot(data.seg(:,1),data.seg(:,2),'r.')

% write .bnd file
if length(dates)==1
   bndname = ['pars_',datestr(dates(1),'yyyymmddHH'),'.bnd'];
else
   bndname = ['pars_',datestr(dates(1),'yyyymmddHH'),'_',...
      datestr(dates(2),'yyyymmddHH'),'.bnd'];
   datnames = [];
   for ii=1:size(data.par,1)
      datnames = strvcat(datnames,['tpar_',num2str(ii),'_',...
         datestr(dates(1),'yyyymmddHH'),'_',...
         datestr(dates(1),'yyyymmddHH'),'.dat']);
   end
end

fid = fopen(bndname,'w');
for ii=1:size(data.seg,1)
   if ii==1
      bform = 'BOUN SEGM XY %13.4f %13.4f &\n';
   else
      bform = '             %13.4f %13.4f &\n';
   end
   fprintf(fid,bform,data.seg(ii,1:2));
end

if length(dates)==1
   for ii=1:size(data.par,1)
      if ii==1
         vpl = '     VAR PAR %11.2f %7.2f %7.2f %7.2f 6.00 &\n';
      elseif ii==size(data.par,1)
         vpl = '             %11.2f %7.2f %7.2f %7.2f 6.00';
      else
         vpl = '             %11.2f %7.2f %7.2f %7.2f 6.00 &\n';
      end
      fprintf(fid,vpl,data.seg(ii,3),data.par(ii,:));
   end
else
   for ii=1:size(data.par,1)
      if ii==1
         vpl = ['     VAR FILE %11.2f ''%s'' 1 &\n'];
      elseif ii==size(data.par,1)
         vpl = ['              %11.2f ''%s'' 1'];
      else
         vpl = ['              %11.2f ''%s'' 1 &\n'];
      end
      fprintf(fid,vpl,data.seg(ii,3),strrep(datnames(ii,:),' ',''));
   end
end
fclose(fid);

% write the tpar files if length(dates)>1
if length(dates)>1
   dates = dates(1):3/24:dates(2);
   for ii=1:size(datnames,1)
      fid = fopen(strrep(datnames(ii,:),' ',''),'w');
      for jj=1:length(dates)
         fprintf(fid,'%s %7.2f %7.2f %7.2f 6.00\n',...
            datestr(dates(jj),'yyyymmdd.HHMM'),data.par(ii,:,jj));
      end
      fclose(fid);
   end
end
