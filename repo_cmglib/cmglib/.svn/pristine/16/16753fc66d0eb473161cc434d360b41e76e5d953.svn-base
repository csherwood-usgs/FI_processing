function j = jacknife( d, w )
% jacknife - Compute jacknife estimate of mean and std. dev.
% j = jacknife( d, w )
%
% d is data, w is weights for data (optional)
%
% Returns:
%    j.mean
%    j.std
%
% See Emery & Thompson, 2001, Data Analysis Methods in 
%    Physical Oceanography, p. 301-303

% chserwood@usgs.gov

d = d(:);
mn = NaN*ones(size(d));
if(exist('w')==1),
   w = w(:);
else
   w = ones(size(d));
end

list = (1:length(d))';
for i = 1:length(d),
   subset = find(list~=i);
   mn(i)= sum( d(subset).*w(subset) )/sum(w(subset));
end
j.mean = mean(mn);
j.std = sum( sqrt( (mn(:)-j.mean).^2 ) );
