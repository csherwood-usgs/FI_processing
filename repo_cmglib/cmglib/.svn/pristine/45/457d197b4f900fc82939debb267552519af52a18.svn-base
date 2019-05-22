function [xout, yout] = ydrag(x, y, varargin)

% ydrag -- Drag points in y-direction only.
%  ydrag('demo') demonstrates itself.
%  ydrag(N) demonstrates itself with N points.
%  [xout, yout] = ydrag(x, y, ...) enables y-dragging
%   of (x, y) points.  Additional arguments are sent intact
%   to the embedded "PLOT" command.  To terminate and return
%   the data, press any active key, or delete the window by
%   clicking on its "goaway" box.  If no output arguments
%   are given, "xout" and "yout" are assigned to the caller's
%   workspace.
%  ... = ydrag(aHandle) enables the line with the given
%   handle, such as the "gco".
%  ... = ydrag (no argument) enables an existing line, with
%   priority given to one that was previously targetted
%   by this routine.
%  xyout = ... returns a two-column array of the (x, y) data.

% Copyright (C) 2001 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
 
% Version of 03-Dec-2001 11:35:44.
% Updated    03-Dec-2001 16:49:55.

persistent theLine
persistent thePointIndex
persistent oldFigureName
persistent oldNumberTitle
persistent oldPointer
persistent deleteFlag

CURSOR = 'circle';

if nargout > 0
	xout = [];
	yout = [];
end

% Demonstration.

if nargin == 1 & isequal(x, 'demo')
	help(mfilename)
	x = 10;
end

if nargin == 1 & length(x) == 1 & ~ischar(x)
	xmax = x;
	set(gcf, 'Name', [mfilename ' ' int2str(xmax)])
	xdemo = 0:xmax;
	ydemo = rand(size(xdemo));
	xo = [];
	yo = [];
	[xo, yo] = feval(mfilename, xdemo, ydemo);
	if nargout == 1
		xout = [xo(:) yo(:)];
	elseif nargout == 2
		xout = xo;
		yout = yo;
	else
		assignin('caller', 'xout', xo)
		assignin('caller', 'yout', yo)
		disp(' ## See results in "xout" and "yout".')
	end
	if any(deleteFlag)
		closereq
	end
	return
end

% New or existing line.

doInitialize = ~~0;
if nargin >= 2 & ~ischar(x)
	if length(varargin) > 0
		theLine = plot(x, y, '*', varargin{:});
	else
		theLine = plot(x, y, '-*', 'MarkerEdgeColor', 'r');
	end
	doInitialize = ~~1;
elseif nargin == 1 & ~ischar(x) & ishandle(x)
	theLine = x;
	doInitialize = ~~1;
elseif nargin < 1
	h = findobj(gca, 'Type', 'line', 'Tag', mfilename);
	if ~any(h)
		h = findobj(gca, 'Type', 'line');
	end
	if ~any(h), return, end
	theLine = h(1);
	doInitialize = ~~1;
end

% Initialize, then wait for keypress.
%  We need to figure out how to handle
%  a premature "CloseRequestFcn" callback.

if doInitialize & ishandle(theLine)
	theFigure = gcf;
	set(theLine, ...
			'ButtonDownFcn', [mfilename ' down'], ...
			'Tag', mfilename)
	try
		zoomsafe on
	catch
	end
	set(theFigure, 'UserData', [])
	set(theFigure, 'KeyPressFcn', [mfilename ' keypress'])
	set(theFigure, 'CloseRequestFcn', [mfilename ' closereq'])
	figure(theFigure)
	waitfor(theFigure, 'UserData', 'done')
	if ishandle(theFigure)
		set(theFigure, 'UserData', [])
		set(theFigure, 'KeyPressFcn', '')
		set(theLine, 'ButtonDownFcn', '')
		if ishandle(theLine)
			xo = get(theLine, 'XData');
			yo = get(theLine, 'YData');
			set(theLine, 'ButtonDownFcn', '')
			if nargout == 1
				xout = [x(:) y(:)];
			elseif nargout == 2
				xout = xo;
				yout = yo;
			else
				assignin('caller', 'xout', xo)
				assignin('caller', 'yout', yo)
				disp(' ## See results in "xout" and "yout".')
			end
		end
	end
end

if nargin == 1 & ischar(x)
	theCommand = x;
	switch theCommand
	case 'closereq'
		deleteFlag = ~~1;
		set(gcbf, 'UserData', 'done')
	case 'keypress'
		set(gcbf, 'UserData', 'done')
	case 'down'
		xy = get(gca, 'CurrentPoint');
		oldFigureName = get(gcbf, 'Name');
		oldNumberTitle = get(gcbf, 'NumberTitle');
		oldPointer = get(gcbf, 'Pointer');
		set(gcbf, 'Name', num2str(xy(1, 2)))
		set(gcbf, 'NumberTitle', 'off')
		set(gcbf, 'Pointer', CURSOR)
		xdown = xy(1, 1);
		ydown = xy(1, 2);
		xlim = get(gca, 'XLim');
		ylim = get(gca, 'YLim');
		xdata = get(theLine, 'XData');
		ydata = get(theLine, 'YData');
		
% Use pixel coordinates to detect point.

		oldUnits = get(gca, 'Units');
		set(gca, 'Units', 'pixels')
		pos = get(gca, 'Position');
		set(gca, 'Units', oldUnits)
		
		width = pos(3);
		height = pos(4);

		xdist = width * (xdata - xdown) / diff(xlim);
		ydist = height * (ydata - ydown) / diff(ylim);
		
		dist = abs(xdist + sqrt(-1) * ydist);
		f = find(dist == min(dist));
		thePointIndex = f(1);
		set(gcbf, 'WindowButtonMotionFcn', [mfilename ' motion'])
		set(gcbf, 'WindowButtonUpFcn', [mfilename ' up'])
	case 'motion'
		xy = get(gca, 'CurrentPoint');
		set(gcbf, 'Name', num2str(xy(1, 2)))
	case 'up'
		xy = get(gca, 'CurrentPoint');
		set(gcbf, 'Name', num2str(xy(1, 2)))
		xup = xy(1, 1);
		yup = xy(1, 2);
		y = get(theLine, 'YData');
		y(thePointIndex) = yup;
		set(theLine, 'YData', y);
		set(gcbf, 'Name', oldFigureName)
		set(gcbf, 'NumberTitle', oldNumberTitle)
		set(gcbf, 'Pointer', oldPointer)
		set(gcbf, 'WindowButtonMotionFcn', '')
		set(gcbf, 'WindowButtonUpFcn', '')
	case 'clear'
	otherwise
	end
end
