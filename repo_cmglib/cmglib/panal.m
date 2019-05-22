% PANAL - Analyze log profile results

%% load log-fit results
%load profiles7022x; gb = (1:length(jtime))';
% load profiles7122; gb = (1:length(jtime))';
load profiles7033; gb = (1:length(jtime))';
if(1)
   load processed_703
   

%%
pr_jt = jtime(gb);
st = sqrt(ut(gb,:).^2+vt(gb,:).^2);
p = P29(gb,1);
[spd,cdir]=pcoord(ut(gb,4),vt(gb,4));
acc = diff(spd)./diff(pr_jt);
acc = [acc(1); acc; acc(end)];
acc = 0.5*[acc(1:end-1)+acc(2:end)]*24*3600; % units of m s-2
ok = find( (abs(st(:,5))>=0.05)&(us4(gb)>0)&(r24(gb)>.9) );
south = find( (abs(st(:,5))>=0.05)&(us4(gb)>0)&(r24(gb)>.9) & (cdir>90) & (cdir <=270) );
north = find( (abs(st(:,5))>=0.05)&(us4(gb)>0)&(r24(gb)>.9) & (cdir<=90) | (cdir > 270) );

%%
figure(1)
clf
[usu,usv]=xycoord(us4(ok),cdir(ok));
[sd1 az1 sd2 az2]=pcastats(usu,usv,.05,1)

%%
figure(2)
clf
subplot(221)
bins = -.1:.005:.1;
[N,X]=hist(p,bins);
[XX,YY]=stairs(bins,100*N./sum(N));
han = fill(XX,YY,[.7 .7 .7]);
set(han,'edgecolor',[.7 .7 .7]);
hold on
[Nok,Xok]=hist(p(ok),bins);
stairs(bins,100*Nok./sum(N), '-k','linewidth',2);
%fill(XXok,YYok,'k')
xlabel('Curvature Parameter')
ylabel('Percent')

subplot(222)
plot( p(north), acc(north), '.b')
hold on
plot( p(south), acc(south), '.r')
ylabel('Acceleration (m s^{-2})')
xlabel('Curvature Parameter')
grid

subplot(223)
plot( p(north), spd(north), '.b')
hold on
plot( p(south), spd(south), '.r')
ylabel('Speed (m s^{-1})')
xlabel('Curvature Parameter')
grid

subplot(224)
plot( p(north), cdir(north), '.b')
hold on
plot( p(south), cdir(south), '.r')
ylabel('Direction {\circ}T')
xlabel('Curvature Parameter')
grid