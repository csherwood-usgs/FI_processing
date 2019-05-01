function [S] = ustats( u, fs )
% USTATS - Calculate statistics on primary wave component in u
% [S] = usstats( u, fs )
% 
% Input:
%    u - arrays of velocities (eg, m/s)
%    fs   - sampling frequency (Hz)
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
% Calc stats on one component ofl rotated velocities
S.mean = mean(u);
S.sd = std(u);
S.rms = rms(u);
S.var = var(u);
S.min = min(u);
S.max = max(u);
S.sk = skewness(u);
S.u3 = sum(u.^3);
S.uspike = mean(u.^3)/mean(u.^2); % Eqn. 3a in Malarky and Davies, 2012
S.S = mean(u.^3)/rms(u.^3);       % Eqn. 4a in Malarky and Davies, 2012
% TODO add eqns. 4b requiring Hilbert transform

a = accel(u, 1./fs);
S.Aa = mean(a.^3)./rms(a.^3);     % Eqn. 4c in Malarky and Davies, 2012
S.aspike = mean(a.^3)/mean(a.^2); % Eqn. 3b in Malarky and Davies, 2012

[S.u_e,S.u_ep]=eskew(u);
dudt = accel( ur, dt );
S.dudt_mn = mean(dudt);
S.dudt_sd = std( dudt );
S.dudt_sk = skewness( dudt );

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