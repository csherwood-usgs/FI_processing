function sg = sgbbl(Ub,Ab,Ur,Zr,Deg,d_median,znot_in,s,nu,eta_def,lab_def,Alpha,beta)
% SGBBL - Styles and Glenn bottom boundary layer model (CRS version of BBLM02) 
% sg = sgbbl(Ub,Ab,Ur,Zr,Deg,d_median,[znot_in,s,nu,eta_def,lab_def,Alpha,beta])
%
% Note: all variables are CGS units

%
% This version based on the code downloaded from the CBLAAST site 11/8/04
% Modified by cherwood@usgs.gov
% Last revised December, 2004
%
% Following changes have been made:
%   1) now outputs a structure and will not work with multiple input
%      values
%   2) no longer takes array of sediment sizes (d)
%   3) if optional parameter znot is present, uses that instead of
%      calculating znot from sed. transport, ripples, etc.
%   4) 


% INPUT:
% Ub       - bottom wave orbital velocity  (cm/s)
% Ab       - bottom wave excursion amplitude (cm)
% Ur       - mean current at a known height above the bed zr  (cm/s)
% zr       - height above the bed the mean current Ur is measured  (cm)
% DEG      - angle between the wave and current (degrees)
% d        - sediment grain size (can be a scalar or vector) (cm)
% d_median - median (or mean) grain diameter used to compute
%            psicr and Psi (cm)
% s        - relative sediment density
% nu       - kinematic viscosity of the water (default is for seawater @ 15C) (cm^2/s)
% eta_def  - default value for ripple height (cm)
% lab_def  - default value for ripple wavelength (cm)
% alpha    - closure constant
% beta     - closure constant
%
% COMPUTED MODEL PARAMETERS
% znot     - hydraulic roughness  (cm)
% kb       - bottom roughness (=30*znot) (cm)
% ETA      - ripple height (cm)
% LAM      - ripple wavelength (cm)
% Psi      - Shields parameter based on skin friction
% psicr    - critical shear stress for initiation of sediment motion.
% mu       - ratio of the magnitude of the maximum wave shear velocity (ustarwm)
%            to the magnitude of the combined shear velocity (ustarcw).
% epsilon  - ratio of the magnitude of the time averaged shear velocity
%            (ustarc) to the magnitude of the combined shear velocity (ustarcw).
% Ro       - internal friction Rossby number =utarcw/(omega*znot)
% sigma    - ub/ustarcw
%
% OUTPUT PRODUCTS
% SKNPRMS - skin friction parameters.
% BBLMPRMS - selected model parameters needed to compute shear stresses &
%            velocity profiles.

% Original disclaimer:
%         This software at its present stage of development is not intended
%         for commercialization.  This software and any copies or derivatives 
%         is intended to be used for research and evaluation purposes only
%         and is provided as is WITHOUT ANY WARRANTY.
%         WARRANTIES OF MERCHANTABILITY AND OF FITNESS FOR A PARTICULAR
%         PURPOSE ARE EXPRESSLY DISCLAIMED.
%         The authors shall not be liable for any loss or damages arising
%         from any use, defect, omission, failure or the like of said 
%         software, nor shall they have any obligation to make available any
%         corrections, improvements, or other modifications or to provide 
%         any assistance or service of any kind.

global alpha kappa z1p mp
kappa=0.4;
g=981;
mm=50;

% Put default params in place if not passed as arguments
if(exist('nu')~=1),
  nu=0.0119;
end
  % Alpha values based on work presented in Styles and Glenn, JGR, (2002)
  % a value of Alpha = 0.3 and Con = 6.4 is presently recommended in
  % the presence of waves over a sandy bed.
  AlphaMatrix=[0.15 0.3 0.5];
  AlphaCon=[9.8 6.4 4.3]; AlphaIdx=2;
  if(exist('Alpha')~=1),
    Alpha=AlphaMatrix(AlphaIdx);
  end
  Con=AlphaCon(AlphaIdx);
if(exist('beta')~=1),
  beta=0.7;
end
if(exist('d_median')~=1),
  d_median=0.04;
end
if(exist('s')~=1),
  s=2.65;
end
tol=1.0e-4;
if(exist('eta_def')~=1),
  eta_def=1;
end
if(exist('lam_def')~=1),
  lam_def=15;
end
kbr_def=3;
star=d_median/(4*nu).*sqrt(g*d_median*(s-1));
psicrnorm=(s-1)*g*d_median;
% den=0.011607+0.0744*d;
% wf=(-3*nu+sqrt(9*nu^2+g*d.^2*(s-1).*...
%   (0.003869+0.0248*d)))./den;
psinorm=(s-1)*g*d_median;
% calculate critcal shear stress for initiation of sediment motion

psicr=shldc(star);

% LOAD the data file with time series of Ub, Ab, Ur, Zr, Deg
%load test2.mat;
%INPUT=DATA;
%Time=INPUT(:,1);
%Ub=INPUT(:,2);
%Ab=INPUT(:,3);
%Ur=INPUT(:,4);
%Zr=INPUT(:,5);
%Deg=INPUT(:,6);
[m n]=size(Ub);
if( (m~=1)|(n~=1) ),error('CRS version is not vectorized'),end

Omega=Ub./Ab;

% compute skin friction shear stress based on Ole Madsen's formula.
arg_ole=Ab/d_median;
fwcskn=exp(5.61*arg_ole.^(-0.109)-7.30);
ustarwmsknole2=0.5*fwcskn.*(1.42*Ub).^2;
Psi=ustarwmsknole2./psinorm;
SKNPRMS=[Psi fwcskn arg_ole];

ETA=zeros(m,1)+eta_def;
LAM=zeros(m,1)+lam_def;

CHI=4*nu*Ub.^2/(d_median*((s-1)*g*d_median)^(1.5));
Iclt=find(CHI < 2);
Icgt=find(CHI >= 2);
if isempty(Icgt) == 1;
  ETA=Ab*0.30.*CHI.^(-0.39);
  LAM=Ab*1.96.*CHI.^(-0.28);
elseif isempty(Iclt) == 1;
  ETA=Ab*0.45.*CHI.^(-0.99);
  LAM=Ab*2.71.*CHI.^(-0.75);
else;
  ETA(Iclt)=Ab(Iclt)*0.30.*CHI(Iclt).^(-0.39);
  ETA(Icgt)=Ab(Icgt)*0.45.*CHI(Icgt).^(-0.99);
  LAM(Iclt)=Ab(Iclt)*1.96.*CHI(Iclt).^(-0.28);
  LAM(Icgt)=Ab(Icgt)*2.71.*CHI(Icgt).^(-0.75);
end;

for i=1:m;
  ub=Ub(i);
  ab=Ab(i);
  deg=Deg(i);
  zr=Zr(i);
  theta=deg*pi/180;
  ur=Ur(i);
  omega=ub/ab;
  if Psi(i)-psicr <= 0;
    kbr=kbr_def;
    ETA(i)=eta_def;
    LAM(i)=lam_def;
  else;
    kbr=Con*ETA(i);
  end;
  ubour=ub/ur;
  ubokur=ubour/kappa;
  kbs=ab*0.0655*(ub^2/((s-1)*g*ab)).^(1.4);
  kb=d_median+kbr+kbs;
  if(exist('znot_in')==1),
 %   fprintf(1,'User input znot: %f\n',znot_in)
    znot = znot_in;
  else
    znot=kb/30;
  end
  zrozn=zr/znot;
  abozn=ab/znot;
  alpha=Alpha*(1+beta*kb/ab);
  z1p=alpha;
  delta=1/sqrt(2*z1p);
  mp=delta+sqrt(-1)*delta;
  t1=-alpha*kappa*z1p;
  a1=1.0e-6;
   
  INP=[abozn zrozn ubokur theta a1];

  OUT1=bstress2(INP);
  fofa=OUT1(end);
% compute pure wave limit for upper bound
  if abozn < 6.25;
      ubouwmgs=1/abs(t1*mp);
  elseif abozn >= 6.25 & abozn < 10;
      ubouwmgs=exp(1.488)*abozn^(-0.653)*...
             abozn^(0.185*log(abozn));
  elseif abozn >= 10 & abozn < 100;
      ubouwmgs=exp(0.4599)*abozn^(0.1977)*...
               abozn^(0.0085*log(abozn));
  elseif abozn >= 100;
      ubouwmgs=exp(0.13996)*abozn^(0.3539)*...
               abozn^(-0.0106*log(abozn));
  end;
  ubouwm=pwave(abozn,ubouwmgs);
  b1=ubouwm;
  fofb=-fofa;
  c1=0.5*(a1+b1);
  INP(end)=c1;
  OUT1=bstress2(INP);
  fofc=OUT1(end);
  cnt1=0;
  for iii=1:mm;
    cnt1=cnt1+1;
    sgn=fofb*fofc;
    if sgn < 0;
       a1=c1;
    else
       b1=c1;
    end;
    c1=0.5*(a1+b1);
    INP(end)=c1;
    OUT1=bstress2(INP);
    fofc=OUT1(end);
    if b1-c1 < tol; break; end;
  end;
  uboucw=c1;  
% format for OUT1: Ro mu epsilon z1ozn z2ozn zroz1 zroz2 fofx
%  BBLMPRMS(i,:)=[OUT1 kbs kbr znot ub ab ur];
  sg.ro = OUT1(i,1);
  sg.mu = OUT1(i,2);
  sg.epsilon = OUT1(i,3);
  sg.z1ozn = OUT1(i,4);
  sg.z2ozn = OUT1(i,5);
  sg.zroz1 = OUT1(i,6);
  sg.zroz2 = OUT1(i,7);
  sg.fofx = OUT1(i,8);
  sg.kbs = kbs;
  sg.kbr = kbr;
end;

%sg.bblmprms = BBLMPRMS;
%SKNPRMS=[Psi fwcskn arg_ole];
sg.psi = Psi;
sg.psicr = psicr;
sg.fwcskn = fwcskn;
sg.arg_ole = arg_ole;
% sg.wf = wf;
% sg.d = d;
sg.d_median = d_median;
sg.eta = ETA;
sg.lam = LAM;
sg.ub = ub;
sg.ab = ab;
sg.ur = ur;
sg.znot = znot;
sg.omega=ub./ab;
sg.ustarcw=sg.ro.*znot.*omega;
sg.ustarc=sg.epsilon*sg.ustarcw;
sg.ustarw = sg.mu*sg.ustarcw;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ubouwm = pwave(Abozn,ubouwmgs);
% UBOUWM - Calculates nondimensional wave shear, phiw 
% ubouwm = pwave(Abozn,ubouwmgs);

% input arguments: znotp,z1p,z2p
% output arguments: phiw
% this program is for the 2-layer model when z1p/znotp > 1;
global kappa z1p mp

tol=1.0e-4;
ubouwmn=0.5234511947*ubouwmgs;
ubouwm=ubouwmgs;
ubouwmo=ubouwmgs*.293847;
phio=.9; phiw=phio;
cnt=0;
while abs((ubouwmn-ubouwm)/ubouwmn) > tol;
  if cnt > 0;
    ubouwmo=ubouwm;
    ubouwm=ubouwmn;
    phio=phiw;
  end;
  if ubouwmn < 0; ubouwm=1.0e-8; end;
  cnt=cnt+1;
  if cnt > 40;
    cnt
      return;
  end;
%  ubouwm=ubouwmn;
  Ro=Abozn/ubouwm;
  znotp=1/(kappa*Ro);

  if z1p/znotp > 1;
    x=[2.*sqrt(znotp),2.*sqrt(z1p)];
    y=x*exp(3*pi*i/4);
    ber=real(bessel(0,y));
    bei=imag(bessel(0,y));
    ker=real(0.5*pi*i*besselh(0,y));
    kei=imag(0.5*pi*i*besselh(0,y));
    ber1=real(bessel(1,y));
    bei1=imag(bessel(1,y));
    berp=(ber1+bei1)/sqrt(2);
    beip=(-ber1+bei1)/sqrt(2);
    ker1=real(0.5*pi*i*besselh(1,y));
    kei1=imag(0.5*pi*i*besselh(1,y));
    kerp=(ker1+kei1)/sqrt(2);
    keip=(-ker1+kei1)/sqrt(2);
% calculate the coefficients for phiw
    bnot=ber(1)+i*bei(1);
    knot=ker(1)+i*kei(1);
    bnotp=(berp(1)+i*beip(1))/sqrt(znotp);
    knotp=(kerp(1)+i*keip(1))/sqrt(znotp);

    b1=ber(2)+i*bei(2);
    k1=ker(2)+i*kei(2);
    b1p=(berp(2)+i*beip(2))/sqrt(z1p);
    k1p=(kerp(2)+i*keip(2))/sqrt(z1p);

    ll=mp*b1+b1p;
    nn=mp*k1+k1p;
    argi=bnotp*nn/(bnot*nn-knot*ll)+...
    knotp*ll/(knot*ll-bnot*nn);
    gammai=-kappa*znotp*argi;
    phiw=abs(gammai);
  else
   gammai=-kappa*z1p*mp;
   phiw=abs(gammai);
  end
  fofsigma=ubouwm-1/phiw;
  fofsigmao=ubouwmo-1/phio;
  ubouwmn=ubouwm-fofsigma*(ubouwm-ubouwmo)/(fofsigma-fofsigmao);
end;
ubouwm=ubouwmn;
return



function Psicr=shldc(Star);
% CALCULATE CRITICAL SHIELDS PARAMETER
% FOR INITIATION OF SEDIMENT MOTION
% FROM SHIELDS DIAGRAM.
% SCF IS CORRECTION FACTOR.  FORMULAT CAN HANDLE MULTIPLE SIZE CLASSES
scf=1;
for i=1:length(Star);
    star=Star(i);
  if star < 1.5;
    psicr=scf*0.0932*star^(-.707);
  elseif star < 4.0;
    psicr=scf*0.0848*star^(-.473);
  elseif star < 10.0;
    psicr=scf*0.0680*star^(-.314);
  elseif star < 34.0;
    psicr=scf*0.033;
  elseif star < 270;
    psicr=scf*0.0134*star^.255;
  elseif star >= 270;
    psicr=scf*0.056;
  end;
  Psicr(i)=psicr;
end;


function OUT1=bstress2(INP);
global alpha kappa z1p mp
%     program bstress2
% Last update:  19-Oct-99
% This program calculates mu for the 3-layer model
% in matlab, but neglects the region > z2 when calculating
% the wave shear.  It is the same as the fortran program
% fric2.f.

Abozn=INP(1);
zrozn=INP(2);
ubokur=INP(3);
theta=INP(4);
uboucw=INP(5);


Ro=Abozn/uboucw;
znotp=1/(kappa*Ro);
if z1p/znotp > 1;
   phi=phi2_1(znotp);
else
   gammai=-kappa*z1p*mp;
   phi=abs(gammai);
end;

mu=sqrt(uboucw*phi);
eps2=-mu^2*abs(cos(theta))+sqrt(1-mu^4*abs(sin(theta)^2));
epsilon=sqrt(eps2);
% if mu > 1; epsilon=1.0e-1; end;

z2p=z1p/epsilon;
Ror=Ro/zrozn;
zroz1=1./(alpha*kappa*Ror);
zroz2=epsilon*zroz1;
z1ozn=alpha*kappa*Ro;
z2ozn=z1ozn/epsilon;


if zroz2 > 1 & z1ozn > 1;
  fofx=ubokur*epsilon*(log(zroz2)+1-epsilon+...
       epsilon*log(z1ozn))-uboucw;
end

if zroz2 <= 1 & zroz1 > 1 & z1ozn > 1;
  fofx=ubokur*epsilon^2*(zroz1-1+log(z1ozn))-uboucw;
end

if zroz1 <= 1 & z1ozn > 1;
  fofx=ubokur*epsilon^2*log(zrozn)-uboucw;
end

if zroz2 > 1 & z1ozn <= 1 & z2ozn > 1;
  fofx=ubokur*epsilon*(log(zroz2)+1....
       -1/z2ozn)-uboucw;
end

if zroz2 <= 1 & zroz1 > 1 &...
  z1ozn <= 1 & z2ozn > 1;
  fofx=ubokur*epsilon^2*(zroz1-1./z1ozn)-uboucw;
end

if zroz2 > 1 & z2ozn <= 1;
  fofx=ubokur*epsilon*log(zrozn)-uboucw;
end
OUT1=[Ro mu epsilon z1ozn z2ozn zroz1 zroz2 fofx];
return


function phi = phi2_1(znotp);
global kappa z1p mp
% calculates nondimensional wave shear, phi
% input arguments: znotp,z1p,z2p
% oputput arguments: phi
% this program is for the 2-layer model when z1p/znotp > 1;
i=sqrt(-1);
x=[2.*sqrt(znotp),2.*sqrt(z1p)];
y=x*exp(3*pi*i/4);
ber=real(bessel(0,y));
bei=imag(bessel(0,y));
ker=real(0.5*pi*i*besselh(0,y));
kei=imag(0.5*pi*i*besselh(0,y));
ber1=real(bessel(1,y));
bei1=imag(bessel(1,y));
berp=(ber1+bei1)/sqrt(2);
beip=(-ber1+bei1)/sqrt(2);
ker1=real(0.5*pi*i*besselh(1,y));
kei1=imag(0.5*pi*i*besselh(1,y));
kerp=(ker1+kei1)/sqrt(2);
keip=(-ker1+kei1)/sqrt(2);
% calculate the coefficients for phi
bnot=ber(1)+i*bei(1);
knot=ker(1)+i*kei(1);
bnotp=(berp(1)+i*beip(1))/sqrt(znotp);
knotp=(kerp(1)+i*keip(1))/sqrt(znotp);

b1=ber(2)+i*bei(2);
k1=ker(2)+i*kei(2);
b1p=(berp(2)+i*beip(2))/sqrt(z1p);
k1p=(kerp(2)+i*keip(2))/sqrt(z1p);

ll=mp*b1+b1p;
nn=mp*k1+k1p;
argi=bnotp*nn/(bnot*nn-knot*ll)+...
knotp*ll/(knot*ll-bnot*nn);
gammai=-kappa*znotp*argi;
phi=abs(gammai);


