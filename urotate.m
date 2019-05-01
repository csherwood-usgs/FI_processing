function [S] = ustats( u, v, fs, rot )
% UBSTATS - Calculate statistics on primary wave component in u, v time
% series. This will rotate the data into the axis of maximum variance and
% calculate statistics for velocities along that axis only.
%
% Note that the major axis can be +/-180, so the sign of the skewness stats
% are arbitrary, and we have oriented toward the east. We need to fix that.
% TODO - add p so we can determine direction of wave propagation

% [S] = usstats( u, v, fs, [rotate] )
% 
% Input:
%    u, v - arrays of velocities (eg, m/s)
%    fs   - sampling frequency (Hz)
%    rot - flag to indicate treatment of rotation, as follows:
%       1 [default] - Major axis toward east
%      -1  - Major axis to west (flip sign of u)
%       0  - No rotation...treat u as given (v is ignored)
%
% Returned:
%    S    - structure of statistics
%    S.umean - mean u (m/s)and many more
%
%
% Uses geographical coord. system, i.e., assumes u postive
% eastward (azimuth = 90), v positive northward (azimuth = 0 = 360 )

% Chris Sherwood, USGS
% Last revised 5/1/2019
if(~exist('rot',1)
   rot = 1;
end

dt = 1 ./ fs;

S.umean = mean(u);
S.vmean = mean(v);
up = u-S.umean;
vp = v-S.vmean;

% Calculate principal components
C = cov(up,vp);
[V,D] = eig(C);

x1 = [.5*sqrt(D(1,1))*V(1,1);-.5*sqrt(D(1,1))*V(1,1)];
y1 = [.5*sqrt(D(1,1))*V(2,1);-.5*sqrt(D(1,1))*V(2,1)];
x2 = [.5*sqrt(D(2,2))*V(1,2);-.5*sqrt(D(2,2))*V(1,2)];
y2 = [.5*sqrt(D(2,2))*V(2,2);-.5*sqrt(D(2,2))*V(2,2)];
[S.mean_spd S.mean_dir]=pcoord( S.umean, S.vmean );
[ majr, S.maj_az ] = pcoord( x1(1), y1(1) );
[ minr, S.min_az ] = pcoord( x2(1), y2(1) );
if(majr < minr),
  ltemp = majr; aztemp = S.maj_az;
  majr = minr;    S.maj_az = S.min_az;
  minr = ltemp;   S.min_az = aztemp;
end
S.maj_sd = 2*majr;
S.min_sd = 2*minr;

if rot==1
% aim major axis toward East
S.maj_az = S.maj_az-(S.maj_az >= 180)*180;
S.min_az = S.maj_az-90;
end
if rot==-1
S.maj_az = S.maj_az+(S.maj_az <= 180)*180;
S.min_az = S.maj_az+90;
end

% Rotate into principal coordinates
[s,d]=pcoord(up,vp);
d = d-S.min_az;
[ur,vr]=xycoord(s,d);

% We want vr


% Calc stats on one component ofl rotated velocities
S.sd = std(ur);
S.sk = skewness(ur);
[S.ur_e,S.ur_ep]=eskew(ur);
dudt = accel( ur, dt );
S.durdt_mn = mean(dudt);
S.durdt_sd = std( dudt );
S.durdt_sk = skewness( dudt );

S.vr_sd = std(vr);
S.vr_sk = skewness(vr);
[S.vr_e,S.vr_ep]=eskew(vr);
dvdt = accel( vr, dt );
S.dvrdt_mn = mean(dvdt);
S.dvrdt_sd = std( dvdt );
S.dvrdt_sk = skewness( dvdt );

