function [out]=iicom(CAM,CAM1);

%IICOM	UI commands for FRF data visualization within the IIGUI interface
%
%       synopsis: iicom
%		  iicom CommandString
%                 iicom('CommandString')
%
%       Accepted commands are
%
%	 +,-        next/previous channel/pole/deformation
%	 cax(i), ca+   make axis (i) or the next axis current
%	 ch(i),chc(i)  show channel(s)/pole(s)/deformation(s) number (i)
%		       (i) is a string and shortcuts such as 1:10 are accepted
%		       ch  is applied to all axes with current plot function
%		       chc only applies the command to the current axis
%	 gui [one,two] [fe,gui]  interface configuration
%	 head [Tit (dx dy LCR Size FontName), clear, AdLabel]  page headers
%		       manipulations
%	 Print	       prints the figure of the current axis
%	 PrePrint      sets the current figure size to the printing size
%	 scaxt(i j)    set current axis to PlotFcn (i) and PlotType (j)
%	 show type     selects plot type (abs,nyq ... pole,err ... line ...)
%	 sub (i j n PlotFcn PlotType ParentFig) creates a division of the
%		       drawing figure in (i) vertical and (j) horizontal
%		       divisions (see SUBPLOT). Optional arguments are n (to
%		       create only one specific axis of the subdivision) and
%	               the default plot function (1: IIPLOT, 2: IDIPLOT, 
%		       3: FEPLOT), plot type , and parent figure.
%	 sub[MagPhase,Iso]  standard subdivisions for magnitude/phase FRF plots
%		       and 4 standard views of a deformed structure
%	 sub[leg]      creates a legend axis associated to the current axis
%	 TitOpt(i)     set (i) for title option parameter (see iiplot, idiplot)
%
%	See also: feplot, iimouse

%       Etienne Balmes  02/02/94, 04/15/96
%       Copyright (c) 1990-2004 by SDTools,
%       All Rights Reserved.

global IIgui IIAxes IIMenu

if nargin==0&nargout==0;  commode iicom; return;  end
if nargin<2 CAM1=' '; end

[CAM,Cam]=comstr(CAM,1);
if any(size(IIgui)<[10 8])  IIgui(10,8)=0; end
if IIgui(3,5)>100 disp(['iicom(''' Cam ''')']); end

% Indicate working status
eval('set(IIgui(5,2),''backgroundcolor'',[1 0 0]);','IIgui(5,2)=0;')

out='done';yplot=0;

if comstr(version,'5') v5=1; else v5=0;end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Command determination section (mostly alphabetical order of commands)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% ---------------------------------------------------------------------------
if length(Cam)==1
 if Cam=='+' CAM='ch+';Cam='ch+';end; if Cam=='-' CAM='ch-';Cam='ch-';end
end

if comstr(Cam,';')  commode('iicom',CAM(2:length(CAM)));
elseif strcmp(Cam,'return')
 if IIgui(5,2)~=0 set(IIgui(5,2),'backgroundcolor',[.706 .706 .706]); end
 return

% ---------------------------------------------------------------------------
elseif comstr(Cam,'ad') [CAM,Cam]=comstr(CAM,3);

% AdButton : address button handling - - - - - - - - - - - - - - - - - - - -
if strcmp(Cam,'button')
	aa=strfind(get(IIgui(5,2),'string'),'idfeme');
	if isempty(aa) | aa==1
	   set(IIgui(5,2),'string','fe');
	   set(IIgui(5,3),'Callback',['commode(''fecom,iicom'', ' ...
               'get(gco,''string''));comgui(''Refresh'');']);
          comgui('@Commands sent to fecom,iicom');
        elseif aa==3
	   set(IIgui(5,2),'string','me');
	   set(IIgui(5,3),'Callback',['commode(''femesh,fecom'', ' ...
               'get(gco,''string''));comgui(''Refresh'');']);
          comgui('@Commands sent to femesh,fecom');
        elseif aa==5
	   set(IIgui(5,2),'string','id');
	   set(IIgui(5,3),'Callback',['commode(''idcom,iicom'', ' ...
               'get(gco,''string''));comgui(''Refresh'');']);
          comgui('@Commands sent to idcom,iicom');
	end


else out='unknown'; end % of ad subcommands - - - - - - - - - - - - - - - - -

% ---------------------------------------------------------------------------
elseif comstr(Cam,'ca')

 if comstr(Cam,'cax')		% change the current axis
        opt=comstr(comstr(CAM,4),[-1 IIgui(3,7)]);
	if opt(1)<=size(IIAxes,1) IIgui(3,7)=opt(1); end

 elseif comstr(Cam,'ca+')	% change the current axis
	IIgui(3,7)=remi(IIgui(3,7)+1,size(IIAxes,1));
        if IIAxes(IIgui(3,7),2)==102
          IIgui(3,7)=remi(IIgui(3,7)+1,size(IIAxes,1));
        end
        eval('figure(IIAxes(IIgui(3,7),3))','');
        if IIgui(3,4)==0 eval('figure(abs(IIgui(1,2)));',';'); end
 else out='unknown';iicom('return');return; end

 if IIAxes(IIgui(3,7),1)==3&~comstr(CAM1,';')  feplot('plot');
 else comgui('Refresh'); end
 

% ---------------------------------------------------------------------------
elseif comstr(Cam,'ch')	% goto ith channel

  if isempty(IIAxes) return; end
  yplot=0;if ~any(1:size(IIAxes,1)==IIgui(3,7)) IIgui(3,7)=1;end
  i1=IIAxes(IIgui(3,7),[1 3]);

  if comstr(Cam,'chc')  ind=IIgui(3,7);CAM=CAM(4:length(CAM));
  elseif i1(1)==3 fecom(CAM);iicom('return');return;
  else ind=find(IIAxes(:,1)==i1(1)&IIAxes(:,3)==i1(2));CAM=CAM(3:length(CAM));
  end
  [CAM,Cam]=comstr(CAM,1);
  if comstr(CAM,'+')
    if length(CAM)>1 opt=comstr(Cam(2:length(Cam)),[-1 1]);else opt=1;end
    opt=IIAxes(ind,21:size(IIAxes,2))+opt;
    i2=IIAxes(ind,20);
  elseif comstr(CAM,'-')
    if length(CAM)>1 opt=comstr(Cam(2:length(Cam)),[-1 1]);else opt=1;end
    opt=IIAxes(ind,21:size(IIAxes,2))-opt;
    i2=IIAxes(ind,20);
  else
   opt=round(comstr(CAM,-1));i2=[];
   opt=ones(length(ind),1)*opt(:)';
  end

  if isempty(opt)
  elseif any(any(opt>0)) & IIAxes(ind(1),1)==1
     for j1=1:length(ind)
      i1=opt(j1,:);
      if ~isempty(i2) %+-commands
        i1=i1(1:i2(j1));if any(i1<1|i1>size(IIxf,2)) i1=[];end
      else i1=i1(find(i1>0&i1<=size(IIxf,2)));end
      if isempty(i1) if isempty(i2)  disp('Invalid channel number(s)');end
      else IIAxes(ind(j1),20+[0:length(i1)])= [length(i1) i1];yplot=1;
      end
     end
  elseif  IIAxes(ind(1),1)==2
     for j1=1:length(ind)
      if IIAxes(ind(j1),15)==1 i3=size(IIpo1,1);else i3=size(IIpo,1);end
      i1=opt(j1,:);if ~isempty(i2)
         i1=i1(1:i2(j1));if any(i1<1|i1>i3) i1=[];end
      else i1=i1(find(i1>0&i1<i3+1));end
      if isempty(i1) if isempty(i2) disp('Invalid pole number(s)');end
      else IIAxes(ind(j1),20+[0:length(i1)])= [length(i1) i1];yplot=1;
      end
     end
  else out='unknown';end

% ---------------------------------------------------------------------------
elseif comstr(Cam,'clear')	% reset defaults

IIAxes=[]; IIgui(3,7)=0;

% ---------------------------------------------------------------------------
% UI command diagnostics
elseif strcmp(Cam,'diag') 

  if IIgui(3,5)==0 IIgui(3,5)=101;else IIgui(3,5)=0; end

% ---------------------------------------------------------------------------
elseif comstr(Cam,'gui') comgui(CAM);

% ---------------------------------------------------------------------------
elseif strcmp(Cam,'gcaxt')	% update current plot type from menu

   ind = get(IIgui(4,3),'userdata');
   i2 = setstr(get(IIgui(4,5),'userdata'));
   ind=ind(get(IIgui(4,3),'value'),abs(i2(1,1)));

   if ind==102  disp('Cannot change to a header axis');comgui('refresh');
   elseif ind~=0
     IIAxes(IIgui(3,7),2)=ind;set(IIgui(4,2),'userdata',IIAxes);yplot=1;
     if IIAxes(IIgui(3,7),1)==3 set(IIAxes(IIgui(3,7),4),'userdata',[]);end
   else % default value
     st=get(IIgui(4,3),'string');st = deblank(st(get(IIgui(4,3),'value'),:));
     yplot=1;
     if comstr(st,'iiplot')
       IIAxes(IIgui(3,7),15:size(IIAxes,2))=zeros(1,size(IIAxes,2)-14);
       IIAxes(IIgui(3,7),[1 2 11:18 20 21])=[1 1 0 90 0 0   IIgui(2,1:4) 1 1];
     elseif comstr(st,'idiplot')
       IIAxes(IIgui(3,7),15:size(IIAxes,2))=zeros(1,size(IIAxes,2)-14);
       IIAxes(IIgui(3,7),[1 2 11:18 20 21])=[2 50 0 90 0 0 IIgui(2,11:14) 1 1];
     elseif comstr(st,'feplot')
       if IIgui(8,1)==0 IIgui(8,1)=2; end
       IIAxes(IIgui(3,7),[1 2 5:size(IIAxes,2)])= ...
           [3 IIgui(8,1) zeros(1,size(IIAxes,2)-4)];
       IIAxes(IIgui(3,7),11:13)=IIgui(9,1:3);
       if IIAxes(IIgui(3,7),1)==3 set(IIAxes(IIgui(3,7),4),'userdata',[]);end
     end; 
   end

% ---------------------------------------------------------------------------
elseif comstr(Cam,'head') [CAM,Cam]=comstr(CAM,5);

ca=findobj(IIgui(1,2),'tag','iihead');
if ~isempty(ca)
  i2=IIAxes(find(IIAxes(:,4)==ca),4);
  if isempty(i2) delete(ca);ca=[];end % axis exist but is not in IIAxes 
end

if isempty(ca) % create a 102 axis
  if ~isempty(IIAxes) i2 = find(IIAxes(:,3)~=IIgui(1,2)|IIAxes(:,2)~=102);
  if ~isempty(i2) IIAxes=IIAxes(i2,:);end;end
  IIAxes(size(IIAxes,1)+1,1:12) = [1 102 IIgui(1,2) 0  3 0 0 1 1  0 0 90];
  figure(IIgui(1,2));
  IIAxes(size(IIAxes,1),4)=axes('units','normalized','position',[0 0 1 1], ...
         'visible','off','tag','iihead','nextplot','add'); ca=gca;
else  set(ca,'nextplot','add'); end
if ~v5 axes(ca);end
st = get(IIgui(1,2),'paperunits'); set(IIgui(1,2),'paperunits','centimeters')
i1 = get(IIgui(1,2),'paperposition'); set(IIgui(1,2),'paperunits',st);
i1=i1(3:4);

% default title - - - - - - - - - - - - - - - - - - 
if comstr(Cam,'titdef')
   st = [deblank(setstr(XFopt(1,:))) '  ' ...
           deblank(setstr(XFopt(2,:))) '        ' ...
	   deblank(setstr(XFopt(3,:)))];
   out=findobj(ca,'tag','IITitDef');
   if isempty(out) if v5 out=text('parent',ca); else out=text;end;end
   set(out,'string',st,'units','normalized','position',[.5 .98 0], ...
      'HorizontalAlignment','center','verticalalignment','top', ...
      'tag','IITitDef');

% title insertion - - - - - - - - - - - - - - - - - - 
% dx dy LCR Size FontName
elseif comstr(Cam,'tit') [CAM,Cam]=comstr(CAM,4);

 [i2,i3,i4,i5]=sscanf(CAM,'%i');
 if i3==4 st=sscanf(CAM(i5:length(CAM)),'%s',1); % read fontname
 else
   i3=[0 0 -1 12]; i2(length(i2)+1:4)= i3(length(i2)+1:4);st='times';
 end
 if i2(3)==-1 st1='left';elseif i2(3)==0 st1='center';else st1='right';end
 i3=[0 0 0];
 if i2(1)>0 i3(1)=i2(1)/i1(1); else i3(1)=1+i2(1)/i1(1); end % xposition
 if i2(2)>0 i3(2)=i2(2)/i1(2); else i3(2)=1+i2(2)/i1(2); end % yposition

 if v5
   out=text('parent',ca,'string',CAM1,'units','normalized','position',i3, ...
    'fontsize',i2(4),'HorizontalAlignment',st1);
 else
   out=text('string',CAM1,'units','normalized','position',i3, ...
    'fontsize',i2(4),'HorizontalAlignment',st1);
 end
% clear all information - - - - - - - - - - - - - - - - - - - - - - - - 
elseif comstr(Cam,'clear') [CAM,Cam]=comstr(CAM,4);

  delete(findobj(ca,'type','text','visible','on'));

% titlelabels - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
elseif comstr(Cam,'adlabel') [CAM,Cam]=comstr(CAM,8);

 if nargin==2
 elseif comstr(Cam,'dof') % dof labels
   i1=[];st=['|x|y|z|rx|ry|rz|-x|-y|-z|-rx|-ry|-rz|s|s|?|'];i3=find(st=='|');
   for j1=1:size(XFdof,1)
         i2=XFdof(j1,1); i4=find([1:12 99 0]==round(rem(i2,1)*100));
         if isempty(i4) i4=15;end
         st1 = sprintf('%i%s',fix(i2),st(i3(i4)+1:i3(i4+1)-1));
         i1(j1,1:length(st1))=st1;
   end
 end
 i2=findobj(0,'type','text','string','AdLabel');
 if isempty(i2)
   if v5 i2=text(0,0,'AdLabel','parent',ca);else i2=text(0,0,'AdLabel');end
 end
 set(i2(1),'visible','off','userdata',i1);
 if isempty(CAM) yplot=1; end

% end of head subcommands - - - - - - - - - - - - - - - - - - - - - - - - - - 
else out='unknown'; end
% ---------------------------------------------------------------------------
elseif strcmp(Cam,'help')	help iicom

% ---------------------------------------------------------------------------
elseif comstr(Cam,'iiguiup')		% which window to leave up

   if (IIgui(3,4)==1|~isempty(strfind(Cam,'off'))) & isempty(strfind(Cam,'on'))
     IIgui(3,4)=0;set(IIMenu(6),'checked','off');comgui('@IIguiUp off');
   else
     IIgui(3,4)=1;set(IIMenu(6),'checked','on'); comgui('@IIguiUp on');
   end

 
% ---------------------------------------------------------------------------
elseif comstr(Cam,'info')	% togles display of poles

  eval('i1=IIAxes(IIgui(3,7),1);','i1=1;');
  i2=get(IIgui(5,2),'string');
  if all(i2=='me') femesh('info');
  else fecom('info');
  end

  disp(' ');disp('The current drawing axes are')
  disp('Num PlotFcn PlotType ParentFig')
  disp(sprintf('%3i %5i  %5i %9i\n',[[1:size(IIAxes,1)]' IIAxes(:,[1:3])]'))
  if nargout==1 out=IIAxes(:,4);end


% ---------------------------------------------------------------------------
elseif comstr(Cam,'preprint')	%set to printing size

  %i3=IIgui(1,2); if IIgui(1,1)==2&i3==2 i3=1; end 
  i1 = get(gcf,'paperposition'); st=get(gcf,'PaperOrientation');
  set(gcf,'units',get(gcf,'paperunit'));
  i2 = get(gcf,'position');  i1(1:2)=i2(1:2);
  set(gcf,'position',i1)

% ---------------------------------------------------------------------------
elseif comstr(Cam,'print')	%print the figure of the current axis

  i1=gcf;figure(IIAxes(IIgui(3,7),11));
  disp(sprintf('Printing figure %i',IIAxes(IIgui(3,7),11)))
  eval(Cam);
  figure(i1); yplot=0;

% ---------------------------------------------------------------------------
elseif comstr(Cam,'scaxt')	% set current axis type

   i2=comstr(comstr(CAM,6),[-1 4 10]);
   i1=get(IIgui(4,3),'userdata');
   if i2(1)>size(i1,2) error('iicom scaxt: not a supported PlotFunction');end
   i1=find(i1(:,i2(1))==i2(2));
   if isempty(i1)     error('iicom scaxt: not a supported PlotType');end 
   if i1(2)==102 disp('Cannot change to a header axis');
   elseif IIAxes(IIgui(3,7),2)==102
     disp('Cannot change a header axis');IIAxes(IIgui(3,7),1)=1;
   else
    IIAxes(IIgui(3,7),1:2)=i2;
    comgui(sprintf('CheckAxes %i %i %i',IIAxes(IIgui(3,7),[1 3]),IIgui(3,7))); 
    yplot=1;
   end

% ---------------------------------------------------------------------------
elseif comstr(Cam,'set') [CAM,Cam]=comstr(CAM,4); % set current axis type

if comstr(Cam,'cax') [CAM,Cam]=comstr(CAM,4);
  i1=IIgui(3,7);if i1<1|i1>size(IIAxes,1) IIgui(3,7)=1;i1=1; end
  if isempty(IIAxes) error('iicom setCax: no declared drawing axes');end
  if comstr(Cam,'p') [CAM,Cam]=comstr(CAM,-6);
   opt=comstr(CAM,[-1 IIAxes(i1,6:9)]);
   set(IIAxes(i1,4),'position',opt);IIAxes(i1,6:9)=opt;
  else out='unknown';end
 
else out='unknown';end

% ---------------------------------------------------------------------------
elseif comstr(Cam,'show')  [CAM,Cam]=comstr(CAM,5);

if ~any(1:size(IIAxes,1)==IIgui(3,7)) IIgui(3,7)=1;end
if isempty(IIAxes) error('Use iigui to initialize a IIGUI axis'); end

i1=IIAxes(IIgui(3,7),:);
if i1(2)==102
  comgui('@Warning: cannot change the type of a header axis');
elseif comstr(Cam,'abs')|comstr(Cam,'amp')|comstr(Cam,'mag')
   IIAxes(IIgui(3,7),1:2)=[1 1];
elseif comstr(Cam,'ar') 	IIAxes(IIgui(3,7),1:2)=[3 7];
elseif comstr(Cam,'li') 	IIAxes(IIgui(3,7),1:2)=[3 2];
elseif comstr(Cam,'pa') 	IIAxes(IIgui(3,7),1:2)=[3 1];
elseif comstr(Cam,'se') 	IIAxes(IIgui(3,7),1:2)=[3 3];
else comgui('@IICOM SHOW : Unknown plot type');
end

if i1(1)==3 & IIAxes(IIgui(3,7),1)~=3
  IIAxes(IIgui(3,7),15:size(IIAxes,2))=zeros(1,size(IIAxes,2)-14);
  IIAxes(IIgui(3,7),19)=2;
elseif i1(1)~=3 & IIAxes(IIgui(3,7),1)==3
  IIAxes(IIgui(3,7),15:size(IIAxes,2))=zeros(1,size(IIAxes,2)-14);
end

comgui(sprintf('CheckAxes %i %i %i',IIAxes(IIgui(3,7),[1 3]),IIgui(3,7)));
if IIAxes(IIgui(3,7),1)==3 set(IIAxes(IIgui(3,7),4),'userdata',[]);yplot=1;
else yplot=1; end

% ---------------------------------------------------------------------------
% sub [ND no delete] [Ret no iiplot] [nv nh nc PlotFcn PlotType ParentFig ND]
% sub [magpha]
elseif comstr(Cam,'sub')  [CAM,Cam]=comstr(CAM,4);

% subplot legend
if comstr(Cam,'clear') IIAxes=[]; CAM='';Cam=''; end
if comstr(Cam,'leg')
  figure(IIgui(1,2));
  if ~any(IIgui(3,7)==1:size(IIAxes,1))
   error('No current IIGUI axis to place legend');
  end
  st = '|inches|centimeters|normalized|points|pixels|';
  i3=IIAxes(IIgui(3,7),:);
  i1=find(st=='|');st=st(i1(i3(5))+1:i1(i3(5)+1)-1);
  figure(IIgui(1,2));
  i3=IIAxes(IIgui(3,7),:);
  IIAxes(size(IIAxes,1)+1,1:13)=[1 101 IIgui(1,2) axes('units',st, ...
      'position',i3(6:9)) i3(5) i3(6:9) 0 0 90 0];
  if IIAxes(IIgui(3,7),1)==3
   IIAxes(size(IIAxes,1),30+[0:i3(30)])=i3(30+[0:i3(30)]);
  else
    IIAxes(size(IIAxes,1),20+[0:i3(20)])=i3(20+[0:i3(20)]);
  end
  iiplot; IIgui(3,7)=size(IIAxes,1);

else %  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
opt=1;yplot=1;

% ND no delete and Ret without iiplot
if ~isempty(strfind(Cam,'nd'))  st1='subplotnd';opt(1)=0;i1=strfind(Cam,'nd');
  [CAM,Cam]=comstr(CAM([1:i1-1 i1+2:length(CAM)]),1);
else st1='subplot';end
if ~isempty(strfind(Cam,'ret'))  yplot=0;i1=strfind(Cam,'ret');
  [CAM,Cam]=comstr(CAM([1:i1-1 i1+3:length(CAM)]),1);
end

% definition of options [nv nh PlotFcn PlotType ParentFig ND]
if     comstr(Cam,'mag')   opt=[2 1 0 1 0 0 opt 1];
elseif comstr(Cam,'iso')   opt=[2 2 0 3 0 0 opt 1];
else     opt=comstr(CAM,[-1 0 0 0 IIgui(1,4) 0 IIgui(1,2) opt 0]);
end

if ~any(1:size(IIAxes,1)==IIgui(3,7)) IIgui(3,7)=1;end  % check IIgui(3,7)
if opt(4)==0 opt(4)=1;end				% default PlotFcn

i2=[];if ~isempty(IIAxes) i2=find(IIAxes(:,1)==1);end  % reset iiplot defaults
if ~isempty(i2) IIgui(2,[1:4])=IIAxes(i2(1),[15:18]);
else            IIgui(2,[1:4])=  [1 111 2 2]; end
i2=[];if ~isempty(IIAxes) i2=find(IIAxes(:,1)==2);end  % reset idiplot default
if ~isempty(i2) IIgui(2,[11:14])=IIAxes(i2(1),[15:18]);
else            IIgui(2,[11:14])=[0 111 2 2]; end
i2=[];if ~isempty(IIAxes) i2=find(IIAxes(:,1)==3);end  % reset feplot default

% i3 contains the defaults
  if ~isempty(IIAxes) i3=find(IIAxes(:,1)==opt(4));else i3=[];end
  if ~isempty(i3) i3=IIAxes(i3(1),:); else i3=zeros(1,35); end
  if opt(4)==3
     i3(1:2)=opt(4:5); if opt(5)==0 i3(2)=IIgui(8,1); end %default plot type
     i3(1,[19 21 22 23 24 26])=IIgui(8,[2:7]);
     i3(11:13)=IIgui(9,1:3);
     if ~any(i3(2)==[1:3 7]) i3(2)=2;end
     if any(opt(1:2)==0) opt(1:2)=[1 1];end
  elseif opt(4)==2
     i3(1:2)=opt(4:5);i3(1,[15:18])=IIgui(2,[11:14]);
     i3(11:12)=[0 90];
     if ~any(i3(2)==[50:54]) i3(2)=50;end
     if any(opt(1:2)==0) opt(1:2)=[2 1];end
  else 
     i3(1:2)=[1 opt(5)];i3(1,[15:18])=IIgui(2,[1:4]);
     i3(11:12)=[0 90];
     if ~any(i3(2)==[1:10 100 101 110]) i3(2)=1;end
     if any(opt(1:2)==0) opt(1:2)=[2 1];end
  end

% set the current figure and its options

  if any(get(0,'children')==opt(6)) IIgui(1,2)=opt(6);end
  if ~any(get(0,'children')==IIgui(1,2)) IIgui(1,2)=3-IIgui(1,1);end
  if IIgui(1,1)==2&IIgui(1,2)==2 IIgui(1,2)=1; end % no axes on toolbar only
  figure(IIgui(1,2))

  set(gcf,'WindowButtonDown','iimouse(''ObjectInfo'')', ...
      'KeyPressFcn','iimouse(''key'');','interrup','on')

if comstr(Cam,'mag') % SubMagPha

  h=[comgui(st1,[.13 .31  .775 .615],'normalized');
     comgui(st1,[.13 .11  .775 .200],'normalized')];
  IIAxes=IIAxes(find(IIAxes(:,4)~=h(1)&IIAxes(:,4)~=h(2)),:);
  IIAxes(size(IIAxes,1)+[1:2],[1:14])= ...
          [1 1  gcf h(1) 3 .13 .31  .775 .615 0 0 90 0 0
           1 2  gcf h(2) 3 .13 .11  .775 .200 0 0 90 0 0];
  IIAxes(size(IIAxes,1)+[-1:0],15:18)= ...
          [IIgui(2,1) 145 IIgui(2,3:4);IIgui(2,1) 110 IIgui(2,3:4)];
  IIAxes(size(IIAxes,1)+[-1:0],20:21)=[1 1;1 1];
  IIgui(3,7)=size(IIAxes,1)-1;

elseif comstr(Cam,'iso') % SubIsoViews

  h= [.1 .5 .38 .38 0    90    0 0;
      .5 .5 .38 .38 0     0   90 0;
      .1 .1 .38 .38 0     0    0 0;
      .5 .1 .38 .38 0 -37.5 30.0 0];
  i3(3)=gcf;ind=ones(4,1)*i3;ind(:,6:13)=h;
  ind(:,4)=[comgui(st1,h(1,1:4),'normalized');
            comgui(st1,h(2,1:4),'normalized');
            comgui(st1,h(3,1:4),'normalized');
            comgui(st1,h(4,1:4),'normalized')];
  IIAxes(size(IIAxes,1)+[1:4],1:size(ind,2))=ind;IIgui(3,7)=size(IIAxes,1)-3;

else %typical subplot division : sub nv nh PlotFcn PlotType ParentFig

     if ~any(1:10==opt(1)) | ~any(1:10==opt(2))
       disp('iicom: subplot separation >10 not accepted');
     else
     ind =1:opt(1)*opt(2); if opt(3)~=0 ind=find(ind==opt(3));end
     st = '|inches|centimeters|normalized|points|pixels|'; j2=size(IIAxes,1);
     for j1 = ind

        h=comgui(st1,[opt(1) opt(2) j1]);
        if size(IIAxes,2)>4 j3=find(IIAxes(:,4)==h); else j3=[];end
        if isempty(j3) j3=size(IIAxes,1)+1; end
        if j1==ind(1) IIgui(3,7)=j3;end % set current axis
        IIAxes(j3,1:length(i3)) = i3;
        ua = length(find(st(1:strfind(st,get(gca,'units'))-1)=='|'));
        IIAxes(j3,[3:9])=[double(gcf) double(h) ua get(h,'position')];
        if IIAxes(j3,1)==3
          if i3(1,31)==1 % self incrementation of deformations
            i2=i3(1,30+[1:5:i3(1,30)*5]);
            IIAxes(j3,30+[0 1:5:i3(1,30)*5]) = [length(i2) i2+j1-1];
          end
          set(h,'userdata',[]);cla;    
        else
         if all(opt(1:2)==[2 1]) & any([1 2]==i3(2)) IIAxes(j3,2)=j1;end
        end
      end
     end % of valid type of subplot division	

end % of subcommand selection

end %  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

i2=ones(size(IIAxes,1),1); %eliminated deleted axes of current figure
for j1=find(IIAxes(:,3)==gcf)' 
 if IIAxes(j1,4)~=0 eval('findobj(IIAxes(j1,4));','i2(j1)=0;');end
end
IIAxes=IIAxes(find(i2),:);

% T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T
% ---------------------------------------------------------------------------
elseif comstr(Cam,'titopt') % title option

   if length(Cam)<7
   elseif Cam(7)=='c' i1=str2num(CAM(8:length(Cam))); i2= IIgui(3,7);
   else  i1=str2num(CAM(7:length(Cam))); i2 = 1:size(IIAxes,1);
   end
    
   if ~isempty(i2)
     i3 = i2(find(IIAxes(i2,1)==1|IIAxes(i2,1)==2));
     if ~isempty(i3) IIAxes(i3,16)=i1(1)*ones(length(i3),1); end
     i3 = i2(find(IIAxes(i2,1)==3));
     if ~isempty(i3) IIAxes(i3,27)=i1(1)*ones(length(i3),1); end
   end
   set(IIgui(4,2),'userdata',IIAxes);yplot=1;



% ---------------------------------------------------------------------------
else
	out='matlab';        eval([CAM CAM1])
% ---------------------------------------------------------------------------
end % of command selection

if IIgui(5,2)~=0 set(IIgui(5,2),'backgroundcolor',[.706 .706 .706]); end









