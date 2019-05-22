function e = nevents(x,thresh,nf)
% NEVENTS - Event statistics
% e = nevents(x,thresh,nf,nmax)
%
% Input
%   x      - 1D array (no NaNs)
%   thresh - threshold defining event [default = median(x)]
%   nf     - length of median filter; must be odd; 1=no filter [default=3]
% Returned
%   e.length      - length of input array
%   e.thresh      - threshold
%   e.pct_exceed  - pct of time >= threshold (after filter)
%   e.na          - number of >= threshold events
%   e.nb          - number of <  threshold non-events
%   e.a(na)       - array of structures with info for each event 
%      e.a.n      -   event number
%      e.a.len    -   length of event 
%      e.a.ext    -   extreme (max) value of event
%      e.a.mean   -   mean value of event
%   e.b(nb)       - array of structures with info for each non-event
%      e.a.n      -   event number
%      e.b.len    -   length of non-event 
%      e.b.ext    -   extreme (min) value of non-event
%      e.b.mean   -   mean value of non-event
%
%  Uses DIST and MEDFILT

%  Chris Sherwood, USGS
%  Last revised November 19, 2003

% I *think* performance will degrade if the number of events is way
% bigger than nmax...but it still works
if(exist('nmax')~=1),nmax = 50;,end;
if(exist('nf')~=1),nf = 1;,end;
if(exist('thresh')~=1),thresh = median(x);, end
x = x(:);
n = length(x);
xf = x;
if(nf>1),xf=medfilt(x,nf);,end
X = (xf>=thresh);
e.length = n;
e.thresh = thresh;
e.pct_exceed = 100*sum(X)/n;

a = repmat( struct('n',NaN,'len',0,'ext',-1e35,'mean',0),nmax,1);
b = repmat( struct('n',NaN,'len',0,'ext', 1e35,'mean',0),nmax,1);

na = 0;
nb = 0;
aflag=1;
i=1;
if(X(1)),aflag=0;end; % perverse flag at start
while(i<=n),
  while(i<=n & ~X(i)),
    if(aflag),
      nb = nb+1;
      aflag = 0;    
  end
    b(nb).n = nb;
    b(nb).len = b(nb).len + 1;
    if(x(i) < b(nb).ext),b(nb).ext=x(i);,end
    b(nb).mean = b(nb).mean + x(i);
    i=i+1;
  end
  while(i<=n & X(i)),
    if(~aflag),
      na = na+1;
      aflag = 1;
    end
    a(na).n = na;
    a(na).len = a(na).len + 1;
    if(x(i) > a(na).ext), a(na).ext=x(i);,end
    a(na).mean = a(na).mean + x(i);
    i=i+1;
  end
end
% number in events + number non-events should add up
if( sum([a(1:na).len])+sum([b(1:nb).len]) ~= n ),fprintf(1,'problema\n'),end

e.na = na;
e.nb = nb;

% calculate mean value of events
for i=1:na,
    a(i).mean = a(i).mean / a(i).len;
end
for i=1:nb,
    b(i).mean = b(i).mean / b(i).len;
end

e.a = a(1:na);
e.b = b(1:nb);

if(1), % diplay results
  fprintf(1,'Results from nevents:\n')
  fprintf(1,'  Record length: %d\n',[e.length])
  fprintf(1,'  Threshold: %f\n',[e.thresh])
  fprintf(1,'  Number of events: %d;   non-events: %d\n',[e.na],[e.nb])
  fprintf(1,'  Percent time threshold exceeded: %f\n',[e.pct_exceed])
  fprintf(1,'\n  Distribution of Event Duration:\n')
  dist( [e.a(:).len] );
  fprintf(1,'\n  Distribution of Non-Event Duration:\n')
  dist( [e.b(:).len] );

  clf
  plot(xf,'-k');
  hold on
  plot(x)
  hold on
  mx = thresh*(~X)+x.*X;
  plot(mx,'-r')
end