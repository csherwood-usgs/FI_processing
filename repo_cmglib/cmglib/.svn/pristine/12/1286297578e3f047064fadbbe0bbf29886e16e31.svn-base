function s = skill(data,model,null_model)
% SKILL - Calculate various skill statistics
% s = skill(data,model,null_model)
%
% Input:
%   data, model, and (optional) null_model must be column vectors of equal
%   length (no NaNs)
% Returns structure s with:
%   s.willmott Eqn. 5 in Warner et al. (2005) JGR 110, C05001
%      based on Willmott, C. J. (1981) On the validation of models. Phys.
%      Geogr. 2:184-194.
%   s.brier Based on CRS notes (see also Sutherland, Peet, Soulsby, 2004)
model = model(:);
data = data(:);
m = length(model(:));
d = length(data(:));
if(m~=d),error('Model and data must be same length in skill.'),end
ok = find(isfinite(model+data));
o = length(ok);
fprintf(1,'%d missing.\n',m-o);
model = model(ok);
data = data(ok);
if(exist('null_model','var')~=1),
    null_model = ones(o,1)*mean(data);
else
    null_model = null_model(ok);
end


s.willmott = 1 -  ( sum( (abs( model - data )).^2 ) ./ ...
   sum( ( abs( model-mean(data) )+ abs( data - mean (data ) )).^2 ));

s.brier = 1 - sum( (data-model) .^2 ) ./ sum( (data-null_model).^2 );
s.RMS = rms(data-model);
R = corrcoef([data model]);
s.R = R(1,2);