function h = stick2(t,u,v,ax)
% STICK2  Vector plot
%
% h = stick2(time,x,y,ax)
%
% Input:
%   time - array of times for x-axis (size n x 1)
%   x, y - x and y arrays with components of vector (each size n x 1)
%   ax   - array with axis parameters [xlow, xhigh, ylow, yhigh]
%
% Returns:
%   h    - handle from plot call
%
% Typical calling sequence:
%   ax = [timemin timemax ymin ymax];
%   orient('anything you want')
%   subplot(i,j,k)    % anything you want
%   h=stick2(time,x,y,ax)
%   set(h,'properties of plot')

% This one figures out the paper orientation, size
% for the horizontal to width ratio
% Chris Sherwood, March 1999
% Tom Gross, May 1993

axis(ax);
pos = get(gca,'Position');
pap = get(gcf,'PaperPosition');
hwratio = (pos(4)*pap(4))/(pos(3)*pap(3));

% Call to Chris' stick.m routine
% vec = stick(t,x,y,ax, nratio);
dt = ax(2)-ax(1);
dv = ax(4)-ax(3);
sf = hwratio*dt/dv;
s = [0; 1; 0];
n = length(t);
us = u*sf;
vec = zeros( n*3, 2 );
id = (1:3');
for i=1:n
   vec(id,:) = s*[us(i), v(i)];
   vec(id,1) = vec(id,1)+ones(3,1)*t(i);
   id = id+3;
end
h = plot(vec(:,1),vec(:,2));
axis(ax)

