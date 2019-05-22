% ACOV_TS - Script for calcuating autocovariance time scale
n = length(x);
n4 = floor(n/4);
[C,lags]=xcov(x,'biased');
lag0 = find(lags==0);
figure(1);
subplot(211)
plot(lags(lag0:lag0+n4),C(lag0:lag0+n4),'-r')
hold on
plot([lags(lag0),lags(lag0+n4)],[0,0],'--k')
T = cumsum(C(lag0:lag0+n4,1))./C(lag0);
subplot(212)
plot(T,'-r')