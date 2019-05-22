function gpx2csv( fnbase )
% gpx2csv - Read Garmin .gpx files, write simple .csv file
% gpx2csv( base_filename) 
% Assumes file has a .gpx extension
% Writes base_filename.csv
% Needs Mapping Toolbox
t = gpxread([fnbase,'.gpx'])
fnout = [fnbase,'.csv']
fid = fopen(fnout,'w');
for i=1:size(t)
    fprintf(fid,'%s, %8.5f, %9.5f, %s\n',t(i).Name,t(i).Latitude,t(i).Longitude,t(i).Time);
end
fclose(fid);