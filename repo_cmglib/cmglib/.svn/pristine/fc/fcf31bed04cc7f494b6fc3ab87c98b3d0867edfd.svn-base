% SHALLOW_UB
slope = 1e-4;
Ua = 20;
dx = 100;
xmax = 50000;
x = (0:dx:xmax)';
F_off = x; % offshore winds
F_on = xmax-x+100000;
h = x*slope;

Hm_off = zeros(length(x),1);
Hm_on = zeros(length(x),1);
Tm_off = zeros(length(x),1);
Tm_on = zeros(length(x),1);
ub_off = zeros(length(x),1);
ub_on = zeros(length(x),1);


for i=1:length(x)
    [Hm_off(i),Tm_off(i)]=fetchlim_shallow(Ua,F_off(i),h(i));
    ub_off(i) = ubqf( Tm_off(i), Hm_off(i), h(i) );
    [Hm_on(i),Tm_on(i)]=fetchlim_shallow(Ua,F_on(i),h(i));
    ub_on(i) = ubqf( Tm_on(i), Hm_on(i), h(i) );
end
Hms_on= min([Hm_on h/2],[],2);
ubs_on = 0.5*Hms_on.*sqrt(9.81 ./h);
Hms_off= min([Hm_off h/2],[],2);
ubs_off = 0.5*Hms_off.*sqrt(9.81 ./h);

clf
han(1)=plot(x,ub_off)
hold on
han(2)=plot(x,ub_on,'-r')
han(3)=plot(x,ubs_off,'--b')
han(4)=plot(x,ubs_on,'--r')
xlabel('Fetch (m)')
ylabel('Orbital Velocity (m/s)')
legend(han,'Offshore','Onshore','SW Offshore','SW Onshore')