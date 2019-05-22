function rip = hlpred2(D,do);
% HLPRED2 - Ripple dimensions from wave-orbitial diam. and grain size
% rip = hlpred2(D,d0);
%
% Input orbital diameter do and grain diameter D in meters
% Returns rip.eta and rip.lamda in meters and rip.stp

% Written by Art Trembanis
% Last revised 1/6/05 csherwood@usgs.gov

%mfile to calculate non-iterative version of Wiberg and Harris ripple model
%after Malarkey and Davies, 2003.  "A non-iterative procedure for the
%Wiberg and Harris (1994) Oscillatory Sand Ripple Predictor."
%JCR vol. 19 no. 3 pp. 738-739.
%compiled by ACT 10/3/03
%inputs do=orbital diameter; D=grain diameter
%output=ripple height (eta); ripple wavelength(lamda); and
%steepness(eta/lamda)

A1=0.095;%constants
A2=0.442;
A3=2.28;
B1=1./A1;
B2=.5*(1+A2).*B1;
B3=B2.^2-A3*B1;
%first solve eqn 4 for do/lam
dolam1=do./(535*D);%sets first guess to anorbital value
doeta1=exp(B2-sqrt(B3-B1*log(dolam1)));%get doeta based on dolam value from above
lamorb=0.62*do;%get lamorb eq2.a
lamanorb=535*D;%get lamanorb eq2.c
if doeta1<20
    dolam=1./0.62;%get dolam eq4.a
    %      lamorb=0.62*do;%get lamorb eq2.a
elseif doeta1>100
    dolam=dolam1;%get dolam eq4.c
    %    lamanorb=535*D;%get lamanorb eq2.c
else
    fdo_etaano=-log(lamorb./lamanorb)*log(0.01*doeta1)./log(5);
    % lamsub=lamanorb.*exp(-fdo_etaano);%get lamsub eq2.b
    dolam=dolam1.*exp(-fdo_etaano);
end
%next solve for do/eta using equation 3
doeta2=exp(B2-sqrt(B3-B1*log(dolam)));%get doeta based on dolam value from eq 4 routine

lamda=do./dolam;%ripple length
eta=do./doeta2;%ripple height
steepness=eta./lamda;
rip.eta = eta;
rip.lamda = lamda;
rip.stp = steepness;
rip.type = -1; 

