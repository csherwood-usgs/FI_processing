function [S] = ubstatsr( u, v, fs )
% UBSTATS - Calculate statistics from ADV burst structure
% [S] = ubstatsr( u, v, fs )
% 
% Input:
%    u, v - arrays of velocities (eg, m/s)
%    fs   - sampling frequency (Hz)
% Returned:
%    S    - structure of statistics
%    S.umean - mean u (m/s)and many more
%
%
% Uses geographical coord. system, i.e., assumes u postive
% eastward (azimuth = 90), v positive northward (azimuth = 0 = 360 )

% Chris Sherwood, USGS
% the ...r indicates revised for Fire Island processing
% Last revised Oct 5, 2001
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

% aim major axis toward East
S.maj_az = S.maj_az-(S.maj_az >= 180)*180;
S.min_az = S.maj_az-90;

% Calc unrotated stats
S.u_sd = std(up);
S.u_sk = skewness(up);
[S.ue,S.uep]=eskew(up);
dudt = accel( up, dt );
S.dudt_mn = mean(dudt);
S.dudt_sd = std( dudt );
S.dudt_sk = skewness( dudt );

S.v_sd = std(vp);
S.v_sk = skewness(vp);
[S.ve,S.vep]=eskew(vp);
dvdt = accel( vp, dt );
S.dvdt_mn = mean(dvdt);
S.dvdt_sd = std( dvdt );
S.dvdt_sk = skewness( dvdt );

% Rotate into principal coordinates
[s,d]=pcoord(up,vp);
d = d-S.min_az;
[ur,vr]=xycoord(s,d);

% Calc rotated stats
S.ur_sd = std(ur);
S.ur_sk = skewness(ur);
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

function a = accel(u,dt)
% accel(u,dt) - First-order difference
a = [0; diff(u)]./dt;

function [E, Ep] = eskew(u)
% ESKEW - Skewness and normalized skewness
% [E, Ep ] = eskew( u )
% Eqns. 1 and 2 in Wright et al., 1994 Marine Geology 118:67-77.
up = u-mean(u);
E = mean( up.^2 .* up );
Ep = E /  (mean(up.^2)).^(3/2);