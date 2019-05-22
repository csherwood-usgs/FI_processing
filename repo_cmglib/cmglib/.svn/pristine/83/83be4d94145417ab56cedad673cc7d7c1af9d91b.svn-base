function h=fig_doc(text_string);
% FIG_DOC - Put documentation data and filename
% h=fig_doc([text_string]);

% from Rocky Geyer on 4/20/09

hh=get(gcf,'currentaxes');
info=dbstack;
dirname=pwd;
info.name;
nam=ans;
titl=[dirname,'\',nam];
if strcmp(titl,'D:\crs\matlab\cmglib\fig_doc.m')==1;
   titl=pwd;
   titl=['no m-file ',titl];
end;

if nargin>0 
   if length(text_string)>1
        titl=[text_string ' ' titl];
    end;
end;

%gh=axes('position',[.6 .02 .45 .015],'visible','off');
gh=axes('position',[.95 .02 .05 .95],'visible','off');
h =text(.5,1,[date ' ' titl],'interpreter','none');
set(h,'fontsize',8,'horizontalalignment','right',...
		 'verticalalignment','top','rotation',90);
%h =text(.5,.9,titl);
%set(h,'fontsize',8,'horizontalalignment','right',...
%		 'verticalalignment','top','interpreter','none','rotation',90);
set(gcf,'currentaxes',hh);

if nargin>0 
   if length(text_string)==1
       eval(['print -dpng -r300 ',nam,int2str(text_string)]);
    end;
end