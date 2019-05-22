function s = diffhist( ta, ya, tb, yb, sa, sb);

% ta is primary time base
isa = 1; % assume ta starts second
isb = find(tb>ta(1),1,'first');
if( ta(1) < tb(1) ),
   isa = find( ta>tb(1),1,'first');
   isb = 1;
end

iea = length(ta); % assume ta ends first
ieb = find(tb<ta(end),1,'last');
if( ta(end) > tb(end) ),
   iea = find( ta<tb(end),1,'last');
   ieb = length(tb);
end

fprintf(1,'Start: A(%d) at %s,  B(%d) at %s \n',...
   isa,datestr(ta(isa)),isb,datestr(tb(isb)));
fprintf(1,'End: A(%d) at %s,  B(%d) at %s \n',...
   iea,datestr(ta(iea)),ieb,datestr(tb(ieb)));

alist = (isa:iea);
if(any(isnan(ta(isa:iea)))),
   fprintf(1,'NaNs in time A\n')
   alist = find( ~isnan(ta(isa:iea)));
end

ybi = interp1(tb,yb,ta(alist),'linear');

subplot(311)
title(sa); hold on
X = (-80:5:80);
N=histc( ya(alist),X );
bar(X,100*N./sum(N),1,'k');
axis([-75 75 0 20])
ylabel('Percent')
d=dist(ya(alist) );
ts = sprintf('N = %d\nMedian = %4.0f\n16-ptile = %4.0f\n84-ptile = %4.0f',d.n,d.pct50,d.pct16,d.pct84)
text(35,15,ts);

subplot(312)
title(sb); hold on
N=histc( ybi,X );
bar(X,100*N./sum(N),1,'k');
axis([-75 75 0 20])
ylabel('Percent')
d=dist( ybi );
ts = sprintf('N = %d\nMedian = %4.0f\n16-ptile = %4.0f\n84-ptile = %4.0f',d.n,d.pct50,d.pct16,d.pct84)
text(35,15,ts);

subplot(313); hold on
title(['Difference: ',sb,' - ',sa])
X = (-80:2:80);
N=histc( ybi-ya(alist), X );
bar(X,100*N./sum(N),1,'k');
axis([-75 75 0 15])
ylabel('Percent')
xlabel('Direction waves are from, relative to due south (\circ)')
d=dist( ybi-ya(alist) );
ts = sprintf('N = %d\nMedian = %4.0f\n16-ptile = %4.0f\n84-ptile = %4.0f',d.n,d.pct50,d.pct16,d.pct84)
text(35,10,ts);
   

