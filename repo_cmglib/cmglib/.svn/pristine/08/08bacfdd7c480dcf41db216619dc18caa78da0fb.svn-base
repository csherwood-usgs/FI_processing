% Flag bad ADV data
% 
% Assumes data are in one big array a
% a = [u v w agc1 agc2 agc3 cor1 cor2 cor3]
% Returns flags in one array

% with
% 
%  1   Velocity (Beam1|X|East)          (m/s)
%  2   Velocity (Beam2|Y|North)         (m/s)
%  3   Velocity (Beam3|Z|Up)            (m/s)
%  4   Amplitude (Beam1)                (counts)
%  5   Amplitude (Beam2)                (counts)
%  6   Amplitude (Beam3)                (counts)
%  7   Correlation (Beam1)              (%)
%  8   Correlation (Beam2)              (%)
%  9   Correlation (Beam3)              (%)

% Array of flags to return

% 1 = Final criteria (0 = good, 1 = bad)
% 2 = NaN in one or more velocities
% 3 = correlation threshold not met
% 4 = amplitude threshold not met
% 5 = speed threshold exceeded
% 6 = std. deviation of velocity threshold exceeded
% 7 = too close to bottom
% 8 = second (higher) cor. threshold not met
% 9 = std dev. of unflagged data from lp filter exceeded
% n = add others as necessary
fs = 8;         % sampling frequency
nflags = 9;
flags = zeros(length(a),nflags);
corr_thresh = 50;
agc_thresh =90;
corr_thresh = 55;
uthresh = .6; % m/s
stdthresh = 2.4; % std deviations
[nr,nc]=size(a);

% NaN in velocity
flags(:,2) = any(isnan(a(:,1:3)'))';

% Flag 3 - Correlation low
flags(:,3) = sum( a(:,7:9)<repmat(corr_thresh,nr,3), 2);

% Flag 4 - Amplitude low
flags(:,4) = sum( a(:,4:6)<repmat(agc_thresh,nr,3), 2);

% Flag 5 - Speed out of bounds
flags(:,5) =  sum( abs(a(:,1:3))>repmat(uthresh,nr,3), 2);

% Flag 6 - Deviation from filtered time series out of bounds
% Using fir1
fc = 1/6;         % 6-s cutoff
Wn = (2/fs)*fc;
Norder = 20;
b = fir1(Norder,Wn,'low',kaiser(21,3));
Hd = designfilt('lowpassfir','FilterOrder',Norder,'CutoffFrequency',fc, ...
       'DesignMethod','window','Window',{@kaiser,3},'SampleRate',fs);
uf = filtfilt(b,1,a(:,1));
vf = filtfilt(b,1,a(:,2));
wf = filtfilt(b,1,a(:,3));
res = a(:,1:3)-[uf vf wf];
rstd = std(res);
flags(:,6) = sum( abs(res) > stdthresh*repmat(rstd,nr,1),2 );

% Flag 7 - Too close to bottom
flags(:,7) = zi<0.05;

% Flag 8 - alternative corr. threshold
flags(:,8) = sum( a(:,7:9)<repmat(65,nr,3), 2);

% TODO - Other tests
% ? du/dt out of bounds
% ? d dir /dt out of bounds

% Final flags depend on critera.
% 2/3 low cor or any bad velocity or 2/3 outliers
if(0) % first file done this way
   flags(:,1) = any( ...
      any( flags(:,5),2 ) + ...
      any( flags(:,3)>1,2) + ...
      any( flags(:,6)>1,2) + ...
      flags(:,7) ...
      ,2 );
end
if(1) % more strict cor criteria
   flags(:,1) = any( ...
      any( flags(:,5),2 ) + ...
      any( flags(:,3)>0,2) +...
      any( flags(:,8)>1,2) + ...
      any( flags(:,6)>1,2) + ...
      flags(:,7) ...
      ,2 );
end
bad = find(flags(:,1));

% repeat residual test using unflagged data
good = find(~flags(:,1));
% make a full array with high residuals
res = 99.99*ones(nr,3);
res(good,:) = a(good,1:3)-[uf(good) vf(good) wf(good)];
rstd = std(res(good,:));
flags(good,9) = sum( abs(res(good,:)) > stdthresh*repmat(rstd,length(good),1),2 );
flags(:,1) =  any( (any(flags(:,1)>0,2) + any(flags(:,9)>1,2)),2 );

figure(3); clf
subplot(411)
sec = (1:nr)/fs;
hp=plot( sec, a(:,1), '.r');set(hp,'color',[.8 .2 .2]);
hold on
hf=plot( sec, uf,'-r');
plot(sec(bad),a(bad,1),'.k')
ylabel('u (m/s)')
txstr = sprintf('%s %6.3f pct. flagged.',datestr(f.ts(i,k)),100*sum(flags(:,1))/nr);
set(gca,'xticklabel',[])
title(txstr)

subplot(412)
hp=plot( sec, a(:,2), '.b'); set(hp,'color',[.2 .2 .8]);
hold on
hf=plot( sec, vf,'-b'); 
plot(sec(bad),a(bad,2),'.k')
ylabel('v (m/s')
set(gca,'xticklabel',[])

subplot(413)
hp=plot( sec, a(:,3), '.k');set(hp,'color',[.8 .8 .8]);
hold on
hf=plot( sec, wf,'-k');
plot(sec(bad),a(bad,3),'.k')
ylabel('w (m/s)')
set(gca,'xticklabel',[])

subplot(414)
hp=plot(sec,100*zi,'.g'); set(hp,'color',[.6 .3 .6])
hold on
hp=plot( sec, a(:,7:9), '.k');set(hp,'color',[.8 .8 .8])
plot(sec(bad),a(bad,7:9),'.k')
ylim([0 200])
ylabel('Corr., Elev (cm)')
xlabel('Time (s)')
drawnow
shg
if(k==1)
   pfn = sprintf('./armqa/q%03d.png',i)
   eval(['print -dpng ',pfn])
end
