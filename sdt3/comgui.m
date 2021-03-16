function [out,out1] = comgui(CAM,CAM1,CAM2)

%comgui general handling of GUI objects for the Structural Dynamics Toolbox
%       
% ---------------------------------------------------------------------------
%	You should never have to call COMGUI by yourself.
% ---------------------------------------------------------------------------
%
% The following information is provided with no guarantee for advanced users
% IImenu variable contents
%
% (1:4, )     [File       Edit       View     Var.]
% (5:10, )    [PoleLine   IIguiUp    xLin     yLin      D.Axes     C.Figure]
% (11, )      [Print]
% (12:16, )   [IIxf       IIxe       IIxh     IIxi      IIpo]
%
% (30:35, )   [FEPlotType Patch      Line     Sens      Unused    Unused]
% (36:39, )   [Undef      Visible    Dot      Line]
% (40:41, )   [Triax      ColType]
%
% IIgui variable contents
%
% (1,1:5) [OneTwoFig PlotFig BackColor PlotTypeDef MouseCheck]
%
% (2, 1:10) [PlotVarII TitoptII ScaleTypeII PoleLineII  5 6 7 8 9 nw]
% (2,11:14) [PlotVarID TitoptID ScaleTypeID PoleLineID]
% (3,1:7)   [CurrentPole AxesFntSize 3 guiUp Diagnos 6 iCurrentAxes]
% (4,1:4)   [FrameBack	  PushCurrentAxes       PopCurrentAxesType 4 ]
% (4,5:6)   [TextObjectInfo 6]
% (5,1:4)   [PushPrevious   PushCommandFcn        EditCommand        PushNext]
% (5,5:8)   [5              EditWMin		PushW              EditWMax]
%
% (6,1:4) [ULength]
% (6,5:8) [CheckAnimation AxesTriax        7    CurrentStrctPtr]
% (7,1:3) [PopSlider1Type TextSlider1Value Slider1]
% (7,4:6) [PopSlider2Type TextSlider2Value Slider2]
% (7,7:9) [PopSlider3Type TextSlider3Value Slider3]
% (8,1:4) [TypeDflt	   UndefDflt	nCycleDflt      nodeDflt]
% (8,5:8) [ScaleDefDflt    PatchTypDflt ColorTypeDflt   EraseMode]
% (9,1:3) [AzimDflt ElevDflt SelfRDflt]
%
% (10,1) Cu			(10,2)
% (10,3:5)  CurrentOrigin X Y Z
% (10,6:8)  CurrentDirection X Y Z
%
% IIAxes variable content
%
% ( ,1:9)    [PltFcn PltTpe ParentFig AxesPtr AxesUnits AxesPositionRect]
% ( ,10:14)  [TextPtr Azim Elev SelfR Unused]
% For feplot
% ( ,15:20) [CenterX CenterY CenterZ CenterScale Undef UndefPtr]
% ( ,21:26) [nCycle ScaleDof ScaleDef PatchType nSelGroup ColorType]
% ( ,29:30) [LinePtr nDef]
% ( ,31:30+nDef*5) [Def# DefScale DefAn DefAn0 DefPtr ...]
% ( ,30+nDef*5+[1:nSelGroup]) [Gp1 Gp2 ...]
%
% ---------------------------------------------------------------------------
%	You should never have to call COMGUI by yourself.
% ---------------------------------------------------------------------------

%       Etienne Balmes  02/02/94, 10/15/96
%       Copyright (c) 1990-2004 by SDTools,
%       All Rights Reserved.

global IIgui IIAxes IIMenu 

if nargin==0 out=3.0; return; else Cam=lower(CAM); end
if any(size(IIgui)<[10 8])  IIgui(10,8)=0; end

if IIgui(3,5)>100 disp(['comgui(''' Cam ''')']); end

%--------------------------------------------------------------------
% general creation
if comstr(Cam,'@')

   if IIgui(4,5)~=0	set(IIgui(4,5),'string',CAM(2:length(CAM)));
   else disp(CAM(2:length(CAM))); end
   return

%--------------------------------------------------------------------
% general creation
elseif comstr(Cam,'gui')   [CAM,Cam] =comstr(CAM,4);

comgui('CheckMenu');
if ~isempty(strfind(Cam,'back')) st=['back'];else st=''; end

set(0,'units','pixels');eval('i1=get(0,''screensize'');','i1=[0 0 800 600];');

if     ~isempty(strfind(Cam,'di'))
  if IIgui(1,4)==3&IIgui(1,1)==2 else IIgui(1,4)=1; end
elseif ~isempty(strfind(Cam,'ii'))  IIgui(1,4)=1;
elseif ~isempty(strfind(Cam,'dfe')) IIgui(1,4)=3;
elseif ~isempty(strfind(Cam,'fe'))  IIgui(1,4)=3;
elseif IIgui(1,4)~=3 IIgui(1,4)=1;  end


if     ~isempty(strfind(Cam,'one')) IIgui(1,1:2)=[1 2];
elseif ~isempty(strfind(Cam,'two')) IIgui(1,1:2)=[2 1];
elseif any(IIgui(1,1)==[1 2])
elseif i1(3)<799|(i1(3)>=799&IIgui(1,4)==3) IIgui(1,1:2)=[1 2];
else IIgui(1,1:2)=[2 1];
end
if isempty(strfind(Cam,'one'))&isempty(strfind(Cam,'two'))
  if i1(3)<799|(i1(3)>=799&IIgui(1,4)==3) Cam='one';
  else Cam='two'; end
end

if  IIgui(1,1)==1 % one figure

   % check of figure position
   i1 = get(2,'position');
   if all(i1(3:4)==[540 120])|all(i1(3:4)==[330 90])
	     set(2,'units',get(0,'defaultfigureunits'), ...
	        'position',get(0,'defaultfigureposition'),'resize','on');
   end

   if ~isempty(IIAxes) % use figure 1 axes if any
      i1=find(IIAxes(:,3)==1);if ~isempty(i1)
        IIAxes=IIAxes(find(IIAxes(:,3)~=2),:);
        i1=find(IIAxes(:,3)==1); IIgui(3,7)=i1(1);for j1=1:length(i1)
            st=get(IIAxes(i1(j1),4),'units');i2=get(IIAxes(i1(j1),4),'units');
            delete(IIAxes(i1(j1),4));figure(2);
            IIAxes(i1(j1),3:4)=[2 comgui('subplot',i2,st)];
        end
      end
   end

elseif IIgui(1,1)==2 %two figure

   IIgui(1,1:2)=[2 1];
   if ~isempty(IIAxes) % use figure 2 axes if any
      i1=find(IIAxes(:,3)==2);if ~isempty(i1)
       IIAxes=IIAxes(find(IIAxes(:,3)~=1),:);
       i1=find(IIAxes(:,3)==2);IIgui(3,7)=i1(1);for j1=1:length(i1)
            st=get(IIAxes(i1(j1),4),'units');i2=get(IIAxes(i1(j1),4),'units');
            delete(IIAxes(i1(j1),4));figure(1);
            IIAxes(i1(j1),3:4)=[1 comgui('subplot',i2,st)];
       end
      end
   end

end
IIgui(1,2)=2; % force in figure 2
comgui(strcat('initbar',st));

%--------------------------------------------------------------------
elseif comstr(Cam,'check')  [CAM,Cam] =comstr(CAM,6);

% CheckAxes - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if comstr(Cam,'axes')

  % options are PlotType PlotFigure AxesNumber

  opt = comstr(comstr(CAM,5),[-1 0 0 0]);
  st = '|inches|centimeters|normalized|points|pixels|';

  % check of IIgui
  if IIgui(3,2)==0 IIgui(3,2)=12; end

  % set default PlotType
  if any([1 2 3]==opt(1)) IIgui(1,4)=opt(1);
  else eval('IIgui(1,4)=IIAxes(IIgui(3,7),1);','IIgui(1,4)=1;');end

  % determine current drawing figure where 0 implies a search for the good one

  if IIgui(1,1)==2&IIgui(1,2)==2 IIgui(1,2)=1; end % no axes on toolbar only
  % if current plot figure wrong set default
  if ~any(get(0,'children')==IIgui(1,2)) IIgui(1,2)=3-IIgui(1,1);end
  if any(get(0,'children')==opt(2))
   if ~(IIgui(1,1)==2&opt(2)==2) IIgui(1,2)=opt(2);end
  end
  IIgui(1,2)=2;figure(IIgui(1,2))

  % i2 contains axes to check
  i2=opt(3:length(opt)); i2=i2(find(i2>0&i2<=size(IIAxes,1)));
  if isempty(i2)  i2=find(IIAxes(:,3)==gcf); %  check all types current figure
  elseif ~isempty(i2) i2=i2(find(IIAxes(i2,1)==IIgui(1,4)));
  end

  for l1=1:size(IIAxes,1)  if any(i2==l1)

     % axis existence and default position
     eval('findobj(IIAxes(l1,4));i1=get(IIAxes(l1,4),''parent'');', ...
          'IIAxes(l1,4)=0;i1=0;');%i1 parent
     if IIAxes(l1,4)~=0 % if good axis retain current position and units
       IIAxes(l1,6:9)=get(IIAxes(l1,4),'position');
       i3=get(IIAxes(l1,4),'units');i4=find(st=='|');
       IIAxes(l1,5)=find(i4==strfind(st,['|' i3 '|']));
     end
     if i1~=IIAxes(l1,3) eval('delete(IIAxes(l1,4));','');IIAxes(l1,4)=0; end
     if IIAxes(l1,4)==0 % reset axis and position
       if all([1:5]~=IIAxes(l1,5)) IIAxes(l1,5)=3; end % default units
       i1=find(st=='|'); i3 = st(i1(IIAxes(l1,5))+1:i1(IIAxes(l1,5)+1)-1);
       IIgui(1,2)=2; IIAxes(l1,3)=2; figure(IIAxes(l1,3));
       IIAxes(l1,4)=comgui('subplot do',IIAxes(l1,6:9),i3);
	% I really want to delete overlaying axes
       if all(IIAxes(l1,1:2)==[1 102])
         set(IIAxes(l1,4),'tag','iihead','visible','off');
       end
     end

     % set/check other properties
      if IIAxes(l1,1)==1
       if IIAxes(l1,20)==0 IIAxes(l1,20:21)=[1 1]; end
       if IIAxes(l1,15)==0 IIAxes(l1,15)=IIgui(2,1); end
      elseif IIAxes(l1,1)==2
       if IIAxes(l1,20)==0 IIAxes(l1,20:21)=[1 1]; end
       if ~any([0 1]==IIAxes(l1,15)) IIAxes(l1,15)=IIgui(2,11); end
      elseif IIAxes(l1,1)==3
       axes(IIAxes(l1,4));
       set(gca,'xscale','lin','yscale','lin','zscale','lin','visible','off');
       if size(IIAxes,2)<31 IIAxes(l1,31)=0; end
       if IIAxes(l1,30)==0  IIAxes(l1,30:31)=[1 1];IIAxes(l1,19)=2; end  % nDef
       if IIAxes(l1,21)<2   IIAxes(l1,21)=15; end	 % nCycle
     end
  end;end % of loop on axes

% CheckMenu - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'menu') [CAM,Cam] =comstr(CAM,5);

  i1=get(0,'children'); if isempty(find(i1==1)) figure(1);end
  i1=get(0,'children'); if isempty(find(i1==2)) figure(2);end
  i1=findobj(2,'type','uimenu');
  if length(IIMenu)<4 IIMenu(4)=0; end
  if isempty(i1)|isempty(i1==IIMenu(3))|isempty(i1==IIMenu(4))
     comgui('initmenu');
  end

% CheckBar - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'bar') [CAM,Cam] =comstr(CAM,5);

comgui('initbar');

out='unknown';end % of sub command selection

%--------------------------------------------------------------------
% Menu creation and handling
elseif comstr(Cam,'initmenu')   [CAM,Cam] =comstr(CAM,9);

if length(get(0,'children')==0) figure(1);set(1,'visible','off');end
IIgui(1,2)=2;figure(IIgui(1,2));

% erase things
eval('delete(IIMenu(3:4))','');
%delete(findobj(2,'type','uimenu'));
i1=findobj(gcf,'type','axes');
for j1=i1(:)' ua=get(j1,'userdata');
 if size(ua,2)<2 elseif ~all(ua(1,1:2)==[Inf 3.5]) delete(j1);end
end

if ~ishandle(2) figure(2);end
% main menus
IIMenu=IIMenu(:);
IIMenu(3,1)=uimenu(2,'label','FEplot');

% sub menus in iiplot - - - - - - - - - - - - - - - - - - - -

IIMenu(10)=uimenu(IIMenu(3),'label','Command Figure','callback','');
uimenu(IIMenu(10),'label','OneFig','callback','comgui(''guione'')');
uimenu(IIMenu(10),'label','TwoFig','callback','comgui(''guitwo'')');

IIMenu(9)=uimenu(IIMenu(3),'label','Drawing Axes','callback','');
uimenu(IIMenu(9),'label','sub 2,1','callback','iicom(''sub 2 1'');');
uimenu(IIMenu(9),'label','sub 1,1','callback','iicom(''sub 1 1'');');
uimenu(IIMenu(9),'label','sub IsoViews','callback','iicom(''subiso'');');

c = computer;st='off';
%if strcmp(c(1:2),'PC')|strcmp(c(1:2),'MA') IIgui(3,4)=1;st='on';end
IIMenu(6)=uimenu(IIMenu(3),'label','IIguiUp', ...
                 'checked',st,'callback','iicom(''IIguiUp'');');

% feplot menus

IIMenu(30)=uimenu(IIMenu(3),'label','DefType','separator','on');
uimenu(IIMenu(30),'label','Patch','callback','feplot(''PlotCreate 4 1'');');
uimenu(IIMenu(30),'label','Line','callback','feplot(''PlotCreate 4 2'');');
uimenu(IIMenu(30),'label','Sensor','callback','feplot(''PlotCreate 4 3'');');
uimenu(IIMenu(30),'label','Arrow','callback','feplot(''PlotCreate 4 7'');');

IIMenu(41)=uimenu(IIMenu(3),'label','ColType');
uimenu(IIMenu(41),'label','Order','callback','fecom(''ColordataUniform'');');
uimenu(IIMenu(41),'label','Element','callback','fecom(''ColordataElt'');');
uimenu(IIMenu(41),'label','Node','callback','fecom(''ColordataNode'');');
uimenu(IIMenu(41),'label','Group','callback','fecom(''ColordataGroup'');');
uimenu(IIMenu(41),'label','Light','callback','fecom(''ColordataLight'');');

IIMenu(36)=uimenu(IIMenu(3), ...
     'label','Undef','separator','off');
IIMenu(37)=uimenu(IIMenu(36), ...
     'label','Visible','callback','fecom(''undef'');');
IIMenu(38)=uimenu(IIMenu(36), ...
     'label','Dot','callback','fecom(''undefdot'');');
IIMenu(39)=uimenu(IIMenu(36), ...
     'label','Line','callback','fecom(''undefline'');');

%erasemode menu
%h=uimenu(IIgui(2,3),'label','EraseMode','separator','on');
%uimenu(h,'label','normal','callback','fecom(''EraseModeNormal'');')
%uimenu(h,'label','back','callback','fecom(''EraseModeBackGround'');')
%uimenu(h,'label','xor','callback','fecom(''EraseModeXOR'');')

IIMenu(40)=uimenu(IIMenu(3), ...
     'label','Triax','callback','fecom(''triax'');','checked','on');


%--------------------------------------------------------------------
% Toolbar creation and handling
elseif comstr(Cam,'initbar') [CAM,Cam] = comstr(CAM,8);

if any(size(IIgui)<[10 8]) IIgui(10,8)=0; end
IIgui(1,2)=2;figure(IIgui(1,2));
if length(IIMenu)<10|~ishandle(IIMenu(10)); comgui('initmenu');end
i2=get(IIMenu(10),'children');

fprintf('%% ------------------------------------------------------\n')
fprintf('%% ------------------------------------------------------\n')
fprintf('%% This OLD version of FEPLOT is provided AS IS. It is   \n')
fprintf('%% know to have many bugs. You can purchase an up to date\n')
fprintf('%% version from SDTools www.sdtools.com                  \n')
fprintf('%% ------------------------------------------------------\n')
fprintf('%% ------------------------------------------------------\n')

% BarIIONE - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if IIgui(1,1)==1 & IIgui(1,4)==1

     Dims =[
	4 1	-.1 -.1 16.2 1.2	.706 0	% initial frame
	4 4	0 4   1 1	.706 0	% PushContextInfo
	5 1	0 2   1 1	.706 0	% PrevChannelButton
	5 4	0 3   1 1	.706 0	% NextChannelButton
	5 2	0 0   1 1	.706 0	% AdText
	5 3	1 0   4 1	.706 0	% AdEdit
	4 2	5 0   1 1	.706 0	% Current axes menu
	4 3	6 0   3 1	.706 0	% CurrentAxesPlotTypeMenu
	4 5     0 1  11 1	.000 1	% ObjectInfoText
	5 6     9 0   3 1	.700 0	% wMinEdit
	5 7    12 0   1 1	.706 0	% wText
	5 8    13 0   3 1	.700 0	% wMaxEdit
	6 5     0 5   1   1	.706 0	% animButton
	];
  set(i2([1]),'checked','on');set(i2([2]),'checked','off');

% BarIITwo- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif IIgui(1,1)==2 & IIgui(1,4)==1

     Dims =[
	4 1	-.1 -.1 11.2 3.2	.706 0	% initial frame
	5 1	0 1   1 1	.706 0	% PrevChannelButton
	5 4	1 1   1 1	.706 0	% NextChannelButton
	5 2	2 1   1 1	.706 0	% CommandTypeText
	5 3	3 1   4 1	.706 0	% CommandEdit
	4 2	7 1   1 1	.706 0	% Current axes menu
	4 3	8 1   3 1	.706 0	% CurrentAxesPlotTypeMenu
	4 5	0 0  11 1	.706 0	% ObjectInfoText
	5 6	0 2   3 1	.706 0	% wMinEdit
	5 7	3 2   1 1	.706 0	% wText
	5 8	4 2   3 1	.706 0	% wMaxEdit
	4 4	7 2   1 1	.706 0	% PushContextInfo
	6 5	8 2   1 1	.706 0	% animButton
	];

  set(gcf,'units','pixels','NumberTitle','on','resize','off','visible','on');
  i1 = get(gcf,'position');
  if ~all(i1(3:4)==[330 90])|any(i1(1:2)<0)
    set(gcf,'resize','on','position',[0 57 330 90],'resize','off');
  end
  set(i2([1]),'checked','on');set(i2([2]),'checked','off');

% BarFEone - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif IIgui(1,1)==1 & IIgui(1,4)==3

     Dims =[
	4 1	-.1 -.1 30.2 1.2 .706 0	% initial frame
	4 4	0 5   1 1	.706 0	% PushContextInfo
	5 1	0 3   1 1	.706 0	% PrevChannelButton
	5 2	0 1   1 1	.706 0	% CommandTypeText
	5 3	0 0   3 1	.706 0	% CommandEdit
	5 4	0 4   1 1	.706 0	% NextChannelButton
	4 5	1 1   10 1	.000 1	% ObjectInfoText
	4 2	0 6   1 1	.706 0	% Current axes menu

	6 5	  0   2   1   1	.706 0	% animButton
	7 1	  3   0   3   1	.706 0  % SliderMenu1
	7 2	  6   0   3   1	.706 0  % SliderValue1
	7 3	  9   0.2 3 0.6	.706 0  % SliderSlider1
	7 4	 12   0   3   1	.706 0  % SliderMenu2
	7 5	 15   0   3   1	.706 0  % SliderValue2
	7 6	 18   0.2 3 0.6	.706 0  % SliderSlider2
%	7 7	 21   0   3   1	.706 0  % SliderMenu3
%	7 8	 24   0   3   1	.706 0  % SliderValue3
%	7 9	 27   0.2 3 0.6	.706 0  % SliderSlider3
	];
  set(i2([2]),'checked','on');set(i2([1]),'checked','off');

% BarFEtwo- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif IIgui(1,1)==2 & IIgui(1,4)==3

     Dims =[
	4 1	-.1 -.1 18.2 7.2	.706 0	% initial frame
	5 1	0 1   1 1	.706 0	% PrevChannelButton
	5 4	1 1   1 1	.706 0	% NextChannelButton
	5 2	2 1   1 1	.706 0	% AdText
	5 3	3 1   4 1	.706 0	% AdEdit
	4 2	7 1   1 1	.706 0	% Current axes menu
	4 3	8 1   3 1	.706 0	% CurrentAxesPlotTypeMenu
	4 5	0 0  18 1	.706 0	% ObjectInfoText
	5 6	11 1   3 1	.706 0	% wMinEdit
	5 7	14 1   1 1	.706 0	% wText
	5 8	15 1   3 1	.706 0	% wMaxEdit

	7 1	0 3   3 1	.706 0  % SliderMenu1
	7 2	3 3   3 1	.706 0  % SliderValue1
	7 3	6 3.2 3 0.7	.706 0  % SliderSlider1
	7 4	9 3   3 1	.706 0  % SliderMenu2
	7 5	12 3   3 1	.706 0  % SliderValue2
	7 6	15 3.2 3 0.7	.706 0  % SliderSlider2
	7 7	9 2   3 1	.706 0  % SliderMenu3
	7 8	12 2   3 1	.706 0  % SliderValue3
	7 9	15 2.2 3 0.7	.706 0  % SliderSlider3
	6 5	0 2   1 1	.706 0	% animButton
	4 4	1 2   1 1	.706 0	% PushContextInfo
	];

  set(gcf,'units','pixels','NumberTitle','on','resize','off','visible','on');
  i1 = get(gcf,'position');
  if ~all(i1(3:4)==[540 120])|any(i1(1:2)<0)
    set(gcf,'resize','on','position',[0 57 540 120],'resize','off');
  end
  set(i2([2 3]),'checked','on');set(i2([1 4]),'checked','off');

else

  if nargin<2 CAM1=[];end
  if size(CAM1,2)<9 error('Bad specification of desired buttons'); end
  Dims=CAM1;

end % of supported type initialisation

% is initialization needed if so do it
Ptr =[4 1;4 2;4 3;4 4;5 1 ;5 4;5 2;5 3;4 2;4 3;4 5;5 6;5 7;5 8;
      7 1;7 2;7 3;7 4;7 5;7 6; 7 7;7 8;7 9;6 5];
i1=0;
for j1=1:size(Ptr,1);
 if IIgui(Ptr(j1,1),Ptr(j1,2))==0 i1=j1;break;
 else
  try; findobj(IIgui(Ptr(j1,1),Ptr(j1,2)));
  catch; i1=j1;break;end
 end
end

if i1>0 % initialize the controls

% start by deleting
for j1=1:size(Ptr,1); eval('delete(IIgui(Ptr(j1,1),Ptr(j1,2)));','');end

% FrameBack
IIgui(4,1)=uicontrol('Style','frame');

% PushCurrentAxes
IIgui(4,2) = uicontrol('Style','pushbutton','String','1','Callback', ...
	'iicom(''ca+'');','userdata',IIAxes);

% PopCurrentAxesType


[out2,outs]=feplot('supported');
out(1:length(out2(:)),3)=out2(:);
i1 = strfind(outs,'|');
IIgui(4,3) = uicontrol('Style','Popup','String', ...
	[outs 'unknown'], ...
	'userdata',out,'Callback','iicom(''gcaxt'');');

% PushContextInfo
IIgui(4,4) = uicontrol('Style','pushbutton','String','IN','Callback', ...
	'iicom(''info'');','userdata',IIAxes);

% TextObjectInfo
IIgui(4,5) = uicontrol('Style','text','String', ...
	'current object info','userdata',outs,'horizontalalignment','left');

% 4,6:7 [CurrentXYPositionLine CurrentPoleText]

% PushPrevious
IIgui(5,1) = uicontrol('Style','pushbutton','string','-', ...
      'Callback','iicom(''ch-'');');
% PushNext
IIgui(5,4) = uicontrol('Style','pushbutton','string','+','Callback', ...
      'iicom(''ch+'');');

% PushCommandFcn
IIgui(5,2) = uicontrol('Style','pushbutton','string','id', ...
      'Callback','iicom(''adbutton'');','tag','idfeme');
% EditCommand
IIgui(5,3)=uicontrol('Style','edit','string','','Callback', ...
   'commode(''idcom,iicom'',get(gco,''string''));comgui(''Refresh'');');

% EditWMin
IIgui(5,6) = uicontrol('Style','edit','string','','Callback',...
     'iicom([''wmin '' get(gco,''string'')]);');

% PushW
IIgui(5,7) = uicontrol('Style','pushbutton','string','w', ...
   'callback','iicom(''wnext'');');

% EditWMax
IIgui(5,8) = uicontrol('Style','edit','string','','Callback',...
     'iicom([''wmax '' get(gco,''string'')]);');

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% toolbar for the FECOM interface
if comstr(version,'5') on='on'; else on='off';end

opt = [1 0 0   2 90 90   9 0 0];
for j1=1:3:7
  IIgui(7,j1)=uicontrol(gcf,'Style','popup','HorizontalAlignment','left', ...
	 'string',['Azim|Elev|SelfR|CenterX|CenterY|CenterZ|CenterSc|' ...
                   'nCycle|ScaleDef'], ...
	 'value', opt(j1), ...
         'CallBack',['fecom(''SetSlider' num2str((j1+2)/3) ''');']);
  IIgui(7,j1+1)=uicontrol(gcf,'Style','edit', ...
         'CallBack',['fecom(''SetSlider' num2str((j1+2)/3) 'Text'');']);
  IIgui(7,j1+2)= uicontrol(gcf,'Style','slider', ...
         'CallBack',['fecom(''SetSlider' num2str((j1+2)/3) 'Value'');']);

end

% CheckAnimation
IIgui(6,5) = uicontrol(gcf,'Style','pushbutton','value',0,'string','AN', ...
	'Interruptible',on,'Callback','fecom(''anim'');');

end % is initialization needed

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% positionning of the bar elements

for j1=1:size(Ptr,1)
   i1=find(Dims(:,1)==Ptr(j1,1)&Dims(:,2)==Ptr(j1,2));
   if ~isempty(i1)
     set(IIgui(Dims(i1,1),Dims(i1,2)),'units','pixels',...
               'position',Dims(i1,3:6)*30+[3 3 -6 -6], ...
               'backgroundcolor',[1 1 1]*Dims(i1,7),...
               'foregroundcolor',[1 1 1]*Dims(i1,8), ...
               'visible','on');
   else
     set(IIgui(Ptr(j1,1),Ptr(j1,2)),'visible','off');
   end
end
if IIgui(1,1)==1 set(IIgui(4,5),'backgroundcolor',get(gcf,'color')); end
if all(get(IIgui(4,5),'backgroundcolor')>.5)& ...
   all(get(IIgui(4,5),'foregroundcolor')>.9)
   set(IIgui(4,5),'foregroundcolor',[0 0 0]);
end
if any(fix(double(i2)/1000)==7); fecom('SetSlider'); end % xxx wrong i2

% end of toolbar initialization - - - - - - - - - - - - - - - - - - - - -

if comstr(Cam,'back')
elseif IIgui(1,4)==3  iicom('sub 1 1 0 3 2');
else                  iicom('sub 2 1 0 1 1'); end

%------------------------------------------------------------------
elseif  comstr(Cam,'refresh')

% reset the values in the iigui window
global IIw IIAxes IIMenu

  eval('set(IIgui(5,3),''string'','''')',';')
  if isempty(IIAxes) comgui('@No defined IIgui axes');return; end
  if IIgui(4,2)==0 return; end; % button not defined

  if any(size(IIgui)<[10 15]) IIgui(10,15)=0;end

  %current axis default
  if all(1:size(IIAxes,1)~=IIgui(3,7)) IIgui(3,7)=0; 
  elseif IIAxes(IIgui(3,7),1)<0 IIgui(3,7)=0;end
  if IIgui(3,7)==0
   i1=[find(IIAxes(:,3)==IIgui(1,2)&IIAxes(:,1)>0);1];
   if ~isempty(i1) IIgui(3,7)=i1(1); else IIgui(3,7)=1;end
  end
  % current plot figure
  IIgui(1,2)=IIAxes(IIgui(3,7),3);

  % enable feplot axes
  st='on';if isempty(find(IIAxes(:,1)==3)) st='off';end
  set(IIMenu([30 36 40 41]),'enable',st);

  % command figure
  aa=strfind('idfeme',get(findall(2,'tag','idfeme'),'string'));
  if IIAxes(IIgui(3,7),1)==3 & aa(1)==1
	   set(IIgui(5,2),'string','fe');
	   set(IIgui(5,3),'Callback',['commode(''fecom,iicom'', ' ...
               'get(gco,''string''));comgui(''Refresh'');']);
  elseif any([1 2]==IIAxes(IIgui(3,7),1)) & any([3 5]==aa(1));
	   set(IIgui(5,2),'string','id');
	   set(IIgui(5,3),'Callback',['commode(''idcom,iicom'', ' ...
               'get(gco,''string''));comgui(''Refresh'');']);
  end

  % PushCurrentAxes       PopCurrentAxesType

  set(IIgui(4,2),'string',num2str(IIgui(3,7)),'value',IIgui(3,7));
  ind = get(IIgui(4,3),'userdata');
  i1 = IIAxes(IIgui(3,7),1); ind = find(ind(:,i1)==IIAxes(IIgui(3,7),2));
  outm = setstr(get(IIgui(4,5),'userdata'));

  if isempty(ind)	ind=get(IIgui(4,3),'max');  end
  set(IIgui(4,3),'value',ind)

  % option check if iigui started before variables defined
  if IIgui(3,4)~=0 figure(2); end

%------------------------------------------------------------------
% comgui('subplot [nd,do]',[PosRect],'units')
elseif  comstr(Cam,'subplot')

 if all(size(CAM1)==[1 4])
   pos = CAM1; unit0 = CAM2;
 elseif any(CAM1(1:2)<1) error('Illegal plot number.')
 elseif CAM1(3)>CAM1(1)*CAM1(2)
                error('Index exceeds number of subplots.')
 else
	% This is the percent offset from the subplot grid of the plotbox.
        off=[2*0.09  2*0.045 2*0.09  2*0.045];
	if CAM1(1) > 2  off(1:2)=off(1:2)*.9; end
	if CAM1(2) > 2  off(3:4)=off(3:4)*.9; end

	row = [(CAM1(1)-1)-fix((CAM1(3)-1)/CAM1(2)) rem(CAM1(3)-1,CAM1(2))];

    % For this to work the default axes position must be normalized coordinates
	def = get(gcf,'DefaultAxesPosition');
        ofs = [def(3)*(off(3)+off(4))/(CAM1(2)-off(3)-off(4)) ...
               def(4)*(off(1)+off(2))/(CAM1(1)-off(1)-off(2))];
	twid = def(3:4)+ofs;
        wid  = twid./CAM1([2 1])-ofs;
        pos = [def(1:2)+row([2 1]).*twid./CAM1([2 1]) wid];

 	if wid(1) <= 0.5*twid(1)/CAM1(2)
		pos(1) = def(1)+row(2)*(def(3)/CAM1(2));
		pos(3) = 0.7*(def(3)/CAM1(2));
	end
	if wid(2) <= 0.5*twid(2)/CAM1(1)
		pos(2) = def(2)+row(1)*(def(4)/CAM1(1));
		pos(4) = 0.7*(def(4)/CAM1(1));
	end
   unit0 = 'normalized';
 end

chi = get(gcf, 'Children'); out = 0;
for i1 = 1:length(chi) if(strcmp(get(chi(i1),'Type'),'axes'))
  units = get(chi(i1),'Units');   set(chi(i1),'Units',unit0);
  cpos = get(chi(i1),'Position'); set(chi(i1),'Units',units);
  % kill if overlap
  if ~((pos(1) >= cpos(1) + cpos(3)) | (cpos(1) >= pos(1) + pos(3)) | ...
       (pos(2) >= cpos(2) + cpos(4)) | (cpos(2) >= pos(2) + pos(4)))
	if any(cpos ~= pos)
          st = get(chi(i1),'tag');
          if comstr(st,'iihead')|comstr(st,'Struct1')|comstr(st,'Triax')
          elseif ~isempty(strfind(Cam,'nd')) % no delete
          elseif ~isempty(strfind(Cam,'do')) % delete non IIAxes
            if size(IIAxes,2)<4 delete(chi(i1));
            elseif ~any(IIAxes(:,4)==chi(i1)) delete(chi(i1));end
          else delete(chi(i1)); end
	else out=chi(i1);set(gcf,'CurrentAxes',chi(i1)); end
	end
end;end

% create the axis:
if out==0
	out = axes('units',unit0,'Position', pos);
        set(out,'units',get(gcf,'defaultaxesunits'))
end

%--------------------------------------------------------------------
end % of command selection


