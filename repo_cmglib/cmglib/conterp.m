function [xi, yi, wi] = conterp(x, y, z, level, w)
% CONTERP - Interpolate along a contour.
%  conterp('demo') demonstrates itself.
%  conterp(N) demonstrates itself with an N+1-by-N+1 array.
%  [xi, yi, wi] = conterp(x, y, z, level, w) interpolates
%   3-D data w along the given contour level in the (x, y, z)
%   2-D array data.  Each row of the result represents one
%   layer in the original w array.  If the contour contains
%   more than one sequence, the pieces are separated by NaNs.
%
% Also see: CONTOURC, INTERP2.
 
% Dr. Charles R. Denham 
% Version of 28-Jul-2003 15:18:24.
% Updated    29-Jul-2003 09:52:44.

if nargin < 1, help(mfilename), x = 'demo'; end

if isequal(x, 'demo'), x = 20; end
if ischar(x), x = eval(x); end

% Demonstration.

if length(x) == 1
    N = x;
    [x, y] = meshgrid(1:N, 1:N);
    x = 0:N;
    y = 0:N;
    z = sort(rand(N+1, N+1));
    level = 0.5;
    w = z;   % Easy.
    tic
    [xxi, yyi, wwi] = feval(mfilename, x, y, z, level, w);
    t = round(toc*100)/100;
    disp([' ## ' mfilename ' ' int2str(N) ...
            ' -- Elapsed time: ' num2str(t) ' s'])
    if nargout > 0
        xi = xxi; yi = yyi; wi = wwi;
    end
    surf(x, y, z)
    hold on
    plot3(xxi, yyi, wwi+0.05, 'r-', 'LineWidth', 3)
    hold off
    s = [mfilename ' ' int2str(N)];
    title(s)
    xlabel x, ylabel y, zlabel z
    set(gca, 'CLim', [0 1])
    colorbar
    view(2)
    drawnow
    set(gcf, 'Name', s)
    figure(gcf)
    pause(2)
    view(3)
    return
end

% Get the contour matrix.

if length(level) == 1, level = [level level]; end

c = contourc(x, y, z, level);

% Get indices of actual (x, y) locations in c.

len = size(c, 2);
k = 1:len;

j = 0;
while j < length(k)
    j = j+1;
    count = c(2, j);
    k(j) = 0;
    j = j + count;
end

k = k(k > 0);

% Get the (x, y) values along the contour.

xxi = c(1, k);
yyi = c(2, k);

% Interpolate each layer of the w array along
%  the (x, y) contour.

sz = size(w);
while length(sz) < 3, sz(end+1) = 1; end
nLayers = sz(3);

wwi = zeros([nLayers length(xxi)]);

for i = 1:nLayers
    wtemp = w(:, :, i);
    wwitemp = interp2(x, y, wtemp, xxi, yyi);
    wwi(i, :) = wwitemp;
end

% Construct NaN-separated vectors for each layer.

xi = zeros(size(wwi)) + NaN;
yi = xi;
wi = xi;

for i = 1:nLayers
    xi(i, k) = xxi;
    yi(i, k) = yyi;
    wi(:, k) = wwi;
end

% Discard the first element, always a NaN.

xi(:, 1) = [];
yi(:, 1) = [];
wi(:, 1) = [];
