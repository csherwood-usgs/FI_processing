function fw = fw_soulsbyi( a_over_zo,imeth )
% FW_SOULSBYI - Interpolation from Table 10, Soulsby (1977)
%
% Input:
%   a_over_zo = a/zo = ub*(T/2*pi) / zo (m/s)
%   imeth = Method: 1 = GM79; 2 = Soulsby DATA2 or DATA13 method
fw = NaN*ones(size(a_over_zo));
t10x = log10([1e2;1e3;1e4;1e5]);
t10y = [.1057   .1268;...
        .0316   .0383;...
        .0135   .0116;...
        .0069   .0035];
%for i=1:length(a_over_zo),
%  fw(i) = interp1(t10x(:,1),t10y(:,imeth),a_over_zo(i),'spline');
%end
fw = interp1(t10x(:,1),t10y(:,imeth),log10(a_over_zo),'spline');
