function S = stats(x)
% STATS  Routine statistics on vector or matrix (mean, pop std dev, min, max)
%
% function S = stats(x)  

% Chris Sherwood, USGS
% March 23, 1999

[m,n] = size(x);
if m == 1,
   m = n;          %handle isolated row vector
end
S(1,:) = (sum(x)/m);
S(2,:) = (sqrt(sum(x.^2)/m-S(1,:).^2));
S(3,:)=  (min(x));                      
S(4,:)=  (max(x));
fprintf('N = %g\n',m);
fprintf('Mean: ');
disp(S(1,:))
fprintf('S.D.: ');
disp(S(2,:))
fprintf('Min:  ');
disp(S(3,:))
fprintf('Max:  ');
disp(S(4,:))
