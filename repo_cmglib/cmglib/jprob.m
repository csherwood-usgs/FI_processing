function n = jprob(data, xbinlims, ybinlims )
% JPROB - Joint probability plot
% n = jprob(data, xbinlims, ybinlims )
% 
% help by Marinna
% where 
%   data  = [x, y] where x and y are columnar vectors
%   xbinlins = the limits between which you want to detect x, or data(:,1)
%   ybinlins = the limits between which you want to detect y, or data(:,2)
%   n = [ny, nx] number of items found between xlims & ylims
ny = length(ybinlims)-1;
nx = length(xbinlims)-1;
xbinctr = xbinlims(1:nx)+diff(xbinlims)/2;
ybinctr = ybinlims(1:ny)+diff(ybinlims)/2;
n = zeros( ny, nx );
for jx=1:nx,
    for iy=1:ny,
        n(iy,jx)= sum( (data(:,1)>xbinlims(jx)) & (data(:,1)<=xbinlims(jx+1)) &...
                       (data(:,2)>ybinlims(iy)) & (data(:,2)<=ybinlims(iy+1)) );
    end
end
pct = 100*n/(sum(n(:)));
pslice(xbinctr,ybinctr,pct);
%contourf(xbinctr,ybinctr,pct);
%surf(xbinctr,ybinctr,n)