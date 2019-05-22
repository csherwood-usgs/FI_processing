function ncvarlist(fn)
% ncvarlist(fn)
% List all variables in a netCDF file using native Matlab netcdf calls
ncid = netcdf.open(fn,'NC_NOWRITE');
[ndims, nvars, ngatts, unlimdimid] = netcdf.inq(ncid);
% list the variables
fprintf(1,'Found variables in %s\n:',fn)
for m=0:nvars-1
   [varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,m);
   disp(varname)
end
netcdf.close(ncid)
