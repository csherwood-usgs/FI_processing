% function ub = ubcalcs( Hs, Tp, h )

Hs = 1;
Tp = 10;
h = 5;
g = 9.81;

hlist = linspace(5,40,20);
ub = zeros(5,length(hlist));

Tz = Tp/1.28
wz = 2*pi./Tz;
wp = 2*pi./Tp;


for i=1:length(hlist)
   h = hlist(i)

% Eqn 28
term = -(((3.65/Tz)*sqrt(h/g)).^2.1);
urms_28 = (Hs/4.)*sqrt(g/h)*exp(term)

Tn = ( h./g ).^(1/2); % Eqn 8
t = Tn./Tz;           % Eqn 27
A = (6500+(0.56 + 15.54.*t).^6).^(1./6.);   % Eqn 26
urms_25 = (0.25*Hs/Tn) ./((1.+A.*t.^2).^3)  % Eqn 25

khz = qkhfs( wz, h )
urms_khTz = ((Hs/2.)*(wz./sinh(khz)) )/sqrt(2.)
khp = qkhfs( wp, h )
urms_khTp = ((Hs/2.)*(wz./sinh(khp)) )/sqrt(2.)

ubs = ubspecfun( Hs, Tp, h )./sqrt(2)

ub(:,i)=[ urms_28, urms_25, urms_khTz, urms_khTp, ubs];
end

figure(1); clf;
han = plot( hlist,ub,'linewidth',2);
legend( han, 'S 28','S 25','lin Tz','lin Tp','WS')

figure(2); clf
plot(hlist, ub(4,:)./ub(5,:))

