function [xyzt,R] = r3d(xyz,pry,gimballed)
% r3d - 3d coordinate transformation
% [xyzt,R] = r3d(xyz,pry,[gimballed])
%    xyz = n x 3 matrix of original points
%    pry = 1 x 3 vector of pitch, roll, and yaw (radians)
%    gimballed =  0 (default) or 1 - determines order of rotation
%                 about z axis
%
%    xyzt = n x 3 matrix of transformed points
%    R = transformation matrix (useful to save for inversions)
%
% The axis and angle names are in normal, right-handed conventions.
% The rotation matrix is Rx*Ry*Rz when gimballed = 0, and Rz*Rx*Ry when
% gimballed = 1. Set gimballed = 0 if the instrument axis is fixed in x y
% (roll pitch) space and the instrument rotates about this axis.

% TODO - check to see if this interpretation of gimballed is standard
% TODO - gimballed = 1 should probably be the default

% csherwood@usgs.gov
[nr,nc]=size(xyz);
if(nc~=3),error('xyz must be n x 3');end
if(exist('gimballed','var')~=1),gimballed = 0; end

pitch = pry(1);
roll = pry(2);
yaw = pry(3);

% Rotation around vertical axis (heading)
Rz = [ cos(yaw)   -sin(yaw)   0        0 ;...
   sin(yaw)    cos(yaw)   0        0 ;...
   0           0          1        0 ;...
   0           0          0        1];
% Rotation around starboard/port  axis (roll)
Ry = [ cos(roll)  0          sin(roll) 0;...
   0           1          0          0;...
   -sin(roll)  0          cos(roll) 0;...
   0           0          0          1];
% Rotation around axis fore/aft(pitch)
Rx = [ 1           0          0         0;...
   0           cos(pitch) -sin(pitch) 0;...
   0           sin(pitch)  cos(pitch) 0;...
   0           0          0         1];

if(~gimballed)
   % Order matters...this is appropriate for instruments where the azimuth
   % is defined relative to the instrument axis...not a compass.
   R = Rx*Ry*Rz;
elseif (gimballed)
   % This is appropriate if yaw == heading...for example, a floating
   % compass
   R = Rz*Rx*Ry;
end

xyzt = NaN*ones(nr,3);
for i=1:nr
   p = [ xyz(i,1) ; xyz(i,2); xyz(i,3); 0 ];
   %pt = Rx*(Ry*(Rz*p)); % original crs code
   pt = R*p;
   xyzt(i,:)=pt(1:3,1)';
end