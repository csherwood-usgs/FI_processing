function ws = wsf( D, nf, rhos, rhow, mu, alpha, beta)
% Fractal settling velocity (Winterwerp, 1988)
% ws = wsf( D, [ nf, rhos, rhow, mu, Re, alpha, beta])
%
% D can be a vector. Other input is optional.
% When nf = 3 (default), returns Stokes settling velocity.
%
% Input:
%    D - particle diameter [m]
%    nf - fractal dimension (optional, range ~1.5 - 3, default = 3) []
%    rhos - particle density (optional, default = 2650) [kg m-3]
%    rhow - water density (optional, default = 1020) [kg m-3]
%    mu - molecular viscosity (optional, default = 1e-3) [Pa s]
%    alpha and beta - shape parameters (optional, default = 1) []
% Returns:
%    ws - settling velocity [m s-1]
%
% Winterwerp and van Kesteren (2004), p 124 - 127.

% csherwood@usgs.gov
% 28 November 2009

g = 9.81;
MAXIT = 10;
acc = 1e-3;

if(exist('nf','var')~=1), nf = 3; end;
if(exist('rhos','var')~=1),rhos = 2650; end
if(exist('rhow','var')~=1),rhow = 1020; end
if(exist('mu','var')~=1),mu=1e-3; end
if(exist('alpha','var')~=1),alpha=1; end
if(exist('beta','var')~=1),beta=1;end
nu = mu/rhow;
D=D(:);
Dp = 4e-6;
[nd,itest] = size(D);
if(itest ~= 1),error('D must be a vector'),end
i=1; maxerr = 999;
% Stokes
wss=NaN*ones(nd,MAXIT);
Re=zeros(nd,MAXIT);
coeff = alpha/(18*beta)*(rhos-rhow)*g/mu*Dp.^(3-nf);
wss(:,i) = coeff .* D.^(nf-1);
while i<MAXIT && maxerr>acc
   i=i+1;
   % Use an average Re to reduce oscillation in ws at each iteration
   Re(:,i) = 0.5*(Re(:,i-1)+(wss(:,i-1).*D)./nu);
   wss(:,i) = coeff .* (D.^(nf-1)./(1+0.15*Re(:,i).^0.687));
   maxerr =  max( abs((wss(:,i-1)-wss(:,i)))./(0.5*(wss(:,i-1)+wss(:,i))));
end
ws = squeeze(wss(:,i));