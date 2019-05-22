% test_soulsby_vanrijn - Script to test / demo soulsby_vanrijn

% put the input info into a structure
% These numbers are from Souslby (1997), p 191
% My subroutine produces different values for term Ass, even though the equation
% looks right. This leads to slightly different final results, when
% compared with the example calcs on p. 192
svin.d50 = .25e-3; % median grain diameter (m)
svin.d90 = .5e-3;  % 90th percentile grain diameter (m)
svin.urms = .26; % rms wave-orbital velocity (m/s)
svin.U = .6;  % depth-mean velocity magnitude(m/s)
svin.tanB = 1/100; % tan(bed slope) in streamwise direction, positive if current is going uphill
svin.rhos = 2650; % sediment density (kg/m3)
svin.rhow =  999; % water density (kg/m3)
svin.nu = 1.14e-6;     % kinematic viscosity of water (m2/s)
svin.h = 5.;            % water depth (m)
svin.zo = 0.006; % bed roughness length (m)

svout = soulsby_vanrijn( svin )

% check orbital velocity calculations
Hs = 3
Tp = 8*1.28     % to make sure Tz = 8
Tz = Tp/1.28    % not sure where I got this, but think it is correct

h = 10.
Urms  = soulsby_ub(Hs, Tz, h )   % Soulsby, based on Tz
[ub,Tbav] = ubspecfun2(Hs,Tp,h)  % Wiberg-Sherwood, based on Tp
Urms2 = ub/sqrt(2)


%% This section reproduces the results on Fig. 31
svin.d50 = .25e-3; % median grain diameter (m)
svin.d90 = .5e-3;  % 90th percentile grain diameter (m)
svin.urms = .26; % rms wave-orbital velocity (m/s)
svin.U = .6;  % depth-mean velocity magnitude(m/s)
svin.tanB = 1/100; % tan(bed slope) in streamwise direction, positive if current is going uphill
svin.rhos = 2650; % sediment density (kg/m3)
svin.rhow =  999; % water density (kg/m3)
svin.nu = 1.14e-6;     % kinematic viscosity of water (m2/s)
svin.h = 5.;            % water depth (m)
svin.zo = 0.006; % bed roughness length (m)
Hs = [3,1,0.001];
Tp = [8,6,1];
U = [.2:.1:2];
clear qt
for i=1:3
   [ub,Tbav] = ubspecfun2(Hs(i),Tp(i),h);
   fprintf('%f %f %f\n',Hs(i),Tp(i),ub)
   svin.urms = ub/sqrt(2.);
   % Basically the same answer from:
   svin.urms = soulsby_ub(Hs(i), Tp(i)/1.28, h );
   for j = 1:length(U)
      svin.U = U(j);
      svout = soulsby_vanrijn( svin );
      qt(i,j) = svout.qt;
   end
end
%% reproduce Fig. 31
qt(qt<=0.)=NaN;
figure(1); clf
semilogy(U,qt(3,:),'-k','linewidth',2)
hold on
semilogy(U,qt(2,:),'-k','linewidth',2)
semilogy(U,qt(1,:),'-k','linewidth',2)
xlabel('Depth-averaged flow velocity (m/s)')
ylabel('Total load transport in current direction (m^2/s)')
set(gca,'xtick',[0:.2:2])
title('Soulsby (1997) Figure 31')
xlim([0. 2.])
ylim([1e-8,1e-1])

%% This section makes a plot that might apply to Barnegat Bay
svin.d50 = .1e-3; % median grain diameter (m)
svin.d90 = .12e-3;  % 90th percentile grain diameter (m)
svin.urms = NaN; % rms wave-orbital velocity (m/s)
svin.U = NaN;  % depth-mean velocity magnitude(m/s)
svin.tanB = 1/100; % tan(bed slope) in streamwise direction, positive if current is going uphill
svin.rhos = 2650; % sediment density (kg/m3)
svin.rhow =  1025; % water density (kg/m3)
svin.nu = 1.14e-6;     % kinematic viscosity of water (m2/s)
svin.h = 2.;            % water depth (m)
svin.zo = 0.006; % bed roughness length (m)
Hs = [1,.75,.5];
Tp = [4,3,2];
U = [.1:.1:2];
clear qt
for i=1:3
   [ub,Tbav] = ubspecfun2(Hs(i),Tp(i),h);
   fprintf('%f %f %f\n',Hs(i),Tp(i),ub)
   svin.urms = ub/sqrt(2.);
   % Basically the same answer from:
   % svin.urms = soulsby_ub(Hs(i), Tp(i)/1.28, h );
   for j = 1:length(U)
      svin.U = U(j);
      svout = soulsby_vanrijn( svin );
      qt(i,j) = svout.qt;
   end
end
%% reproduce Fig. 31
qt(qt<=0.)=NaN;
figure(2); clf
semilogy(U,qt(3,:),'-k','linewidth',2)
hold on
semilogy(U,qt(2,:),'-k','linewidth',2)
semilogy(U,qt(1,:),'-k','linewidth',2)
xlabel('Depth-averaged flow velocity (m/s)')
ylabel('Total load transport in current direction (m^2/s)')
set(gca,'xtick',[0:.2:2])
title('Hs = 1., .75, and 0.5; T = 4, 3, 2; h = 2')
xlim([0. 2.])
ylim([1e-8,1e-1])
