function v = prctile_mat(x,p)
% PRCTILE_MAT - Return p'th percentile of
%
% Calculates percentiles of the 3 rd dimension of 3d array x

% csherwood@usgs.gov 25 Nov 2012
[nx, ny, nz] = size(x);
np = length(p);
v = NaN*ones(nx,ny,np);
for j=1:ny
   for i=1:nx      
      z = sort(squeeze(x(i,j,:)));
      zmin = z(1);
      zmax = z(nz);
      index = 100*cumsum(ones(nz,1))/nz;
      %pct1 =  interp1(index,z,1);
      v(i,j,:) =  interp1(index,z,p);
   end
end
return






