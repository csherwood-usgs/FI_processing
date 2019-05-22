function fig=cmgtssview(jdaytime,dep,temp)%View profile data of temperature, salinity or density collected from moorings.% % Syntax:  cmgtssview(jdaytime,dep,temp);% % jdaytime = time vector in true julian days% temp = temperature, salinity or density matrix of size (M x N) where% M = length of jdaytime% depth = array of water depths where measurements are made% % jpx, USGS, 03-03-00global yhand tssdata temperif nargin<1 help(mfilename); return; end;	temp(temp>1e5)=nan;tssdata.temper=temp;tssdata.time=jdaytime;tssdata.depth=dep;tssdata.maxtemper=max(max(tssdata.temper));tssdata.mintemper=min(min(tssdata.temper));tssdata.range=tssdata.maxtemper-tssdata.mintemper;colordepth=100;yhand.temper=figure('Integerhandle','off',...	'NumberTitle','off', ...	'name','TSS VIEW',...	'unit','norm',...	'colormap',jet(colordepth),...	'Tag','fig tempview 1');b = uimenu(    'Parent',yhand.temper, ...        'Label','[Save As]', ...         'Tag','File Menu');c = uimenu (    'Parent', b, ...            'Label', 'Save as M-file', ...            'Callback', 'timeplt_command print_to_mfile', ...            'Tag', 'print to mfile menu item' );c = uimenu (    'Parent', b, ...            'Label', 'Save as JPEG', ...            'Callback', 'timeplt_command print_jpeg', ...            'Tag', 'print jpeg menu item' );c = uimenu (    'Parent', b, ...            'Label', 'Save as PS', ...            'Callback', 'timeplt_command print_ps', ...            'Tag', 'print ps menu item' );c = uimenu (    'Parent', b, ...            'Label', 'Save as EPS', ...            'Callback', 'timeplt_command print_eps', ...            'Tag', 'print eps menu item' );c = uimenu (    'Parent', b, ...            'Label', 'Print to Printer', ...            'Callback', 'timeplt_command print_to_printer', ...            'Tag', 'print to printer menu item' );if findstr(inputname(3),'sigma')	set(yhand.temper,'colormap',flipud(jet(colordepth)));	tssdata.label='\sigma-\theta';	tssdata.unit='\sigma';elseif findstr(inputname(3),'temp')	tssdata.label='Temperature';	tssdata.unit='Deg. C';elseif findstr(inputname(3),'sal')	tssdata.label='Salinity';	tssdata.unit='ppt';else	tssdata.label='  ';	tssdata.unit='  ';end;set(yhand.temper,'DoubleBuffer','on');yhand.vprofile=axes('parent',yhand.temper,...	'position',[1/32 19/32 6/32 12/32],...	'xgrid','on');yhand.image=axes('parent',yhand.temper,...	'position',[9/32 19/32 22/32 12/32],...	'xticklabel','');yhand.hprofile=axes('parent',yhand.temper,...	'position',[9/32 5/32 22/32 12/32]);yhand.txtbox=axes('parent',yhand.temper,...	'position',[1/32 10/32 6/32 7/32],'visible','off');yhand.txt=text(0.1,0.6,'');yhand.ensambleslide = uicontrol('Parent',yhand.temper,...	'units','norm',...	'Position',[9/32 2/32 22/32 3/32], ...	'Style','slider',...	'callback','tssviewing(1)'...	);prenxtstr={'<<Prev.','Next>>'};prenxtpos={[8/32 0 3/32 2/32],[12/32 0 3/32 2/32]};for i=1:2	yhand.prenxt(i) = uicontrol('Parent',yhand.temper, ...		'Units','norm', ...		'Callback','tssviewing(4);', ...		'Style','pushbutton',...		'String',prenxtstr{i},...		'position',prenxtpos{i});end;yhand.windl = uicontrol('Parent',yhand.temper, ...	'units','norm',...	'callback','tssviewing(10)',...	'HorizontalAlignment','left', ...	'Position',[27/32 0 4/32 2/32], ...	'String',{'3 days','1 day','7 days','15 dyas'}, ...	'Style','popupmenu', ...	'Value',1,...	'TooltipString','Window Length');wl=24*[3 1 7 15];temper.wl=wl(get(yhand.windl,'value'));tssviewing(1);