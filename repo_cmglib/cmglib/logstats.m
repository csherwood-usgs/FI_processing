function [m,v,cv]=logstats(x)
% mean, variance, and coefficient of variation of data with logarithmic distribution
% [m, v, cv]=logstats(x)
x(x<=0)=NaN;
mu = nanmean(log(x));
sig2 = nanvar(log(x));
m = exp(mu+sig2./2);
v = (exp(sig2)-1).*exp(2*mu+sig2);
cv = sqrt(exp(sig2)-1);