doH = logspace( 0, 3, 20)
ln_doH = log(doH);

stp_new = exp(-0.0950*ln_doH.^2 + ...
		0.4415*ln_doH-2.2827)

loglog( doH, stp_new)
hold on
loglog( [1 10],[.17 .17],'-k')
