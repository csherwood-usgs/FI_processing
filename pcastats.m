function [sd1,az1,sd2,az2]=pcastats(u,v,s,ipost)
% PCASTATS  Principal components of 2-d (e.g. current meter) data
%
% function [sd1 az1 sd2 az2]=pcastats(u,v,[[s],ipost])
%
% u and v are column vectors with x and y data
% s is an optional length for the axes
% ipost is an optional flag 0=no action, >1=plot input data
 
% Chris Sherwood, USGS
% Sept 26, 2005
% minor edits 30 Dec 2015

if(nargin <4),ipost=0;end
if(nargin <3),ipost=0;s=0;end

mu = mean(u);
mv = mean(v);
C = cov(u,v);
[V,D] = eig(C);

x1 = [.5*sqrt(D(1,1))*V(1,1);-.5*sqrt(D(1,1))*V(1,1)];
y1 = [.5*sqrt(D(1,1))*V(2,1);-.5*sqrt(D(1,1))*V(2,1)];
x2 = [.5*sqrt(D(2,2))*V(1,2);-.5*sqrt(D(2,2))*V(1,2)];
y2 = [.5*sqrt(D(2,2))*V(2,2);-.5*sqrt(D(2,2))*V(2,2)];
[mspd, mdir]=pcoord( mu, mv );
[ l1, az1 ] = pcoord( x1(1), y1(1) );
[ l2, az2 ] = pcoord( x2(1), y2(1) );
if(l1 < l2),
  ltemp = l1; aztemp = az1;
  l1 = l2;    az1 = az2;
  l2 = ltemp; az2 = aztemp;
end
sd1 = 2*l1;
sd2 = 2*l2;

if(ipost),
   hh1=plot(u(:),v(:),'ob');
   set(hh1,'markerfacecolor',[.9 .2 .2],'markeredgecolor',[.2 .2 .2],'markersize',4);
   hold on
   axis('square');
   axis([-s s -s s])
   ts = ['Speed= ',sprintf('%4.2f',mspd),'; Dir.= ',sprintf('%3.0f',mdir) ];
   hh2=text(-.8*s,.8*s,ts);
   set(hh2,'fontsize',14);
   for i=1:2
      eval(['[ leng(i), az(i) ] = pcoord( x' num2str(i) '(1), y' num2str(i) '(1) );']);
   end
   ts = ['Major axis: Mag.= ',sprintf('%4.2f',sd1),'; Az.= ',sprintf('%3.0f',az1) ];
   hh2=text(-.8*s,.6*s,ts);
   set(hh2,'fontsize',14);
   hh3=plot(x1,y1,'-r',x2,y2,'color',[1 .1 .1],'linewidth',2);
   hh4=plot([0; mu],[0; mv],'color',[.1 .1 1],'linewidth',2.5);
   xlabel('East ({\itu}) component (cm/s)','fontsize',14);
   ylabel('North ({\itv}) component (cm/s)','fontsize',14);
   grid
end
