function a = SkAs( eta )
% Calculate skewness and asymmetry
% (Matlab Hilber function returns analytic signal X = Xr + i*Xi
% such that Xi is the Hilbert transform of real vector Xr)
H = imag( hilbert( eta ) );
eta = eta-mean(eta);
denom=1./mean(eta.^2)^(3./2.);
a.Sk = mean(eta.^3)*denom;
a.As = mean(H.^3)*denom;
return
