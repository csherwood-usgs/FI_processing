function [drv c] = which_computer
% which_computer - Returns drive letter and computer name
% [drv c] = which_computer

% csherwood@usgs.gov
c = getenv('COMPUTERNAME');
if(strcmpi(c,'IGSAGIEGLTCSH70'))
   drv = 'C:';
   fprintf(1,'CRS Laptop %s\n',c)
   return
end
if(strcmpi(c,'IGSAGIEGLTCSH72'))
   drv = 'C:';
   fprintf(1,'CRS HP Laptop %s\n',c)
   return
end
if(strcmpi(c,'IGSAGIEGWSCSH70'))
   drv = 'D:';
   fprintf(1,'CRS Workstation %s\n',c)
   return
end
drv = '?:';
return
