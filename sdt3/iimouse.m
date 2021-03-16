function iimouse(CAM)

% IIMOUSE mouse related call-backs used within the graphical user interfaces
%
%	Synopsis: iimouse('ModeName')
%
%	The only time you may want to use IIMOUSE by yourself is to set it on
%	or off using IIMOUSE ON and IIMOUSE OFF.
%	
%	Supported modes are
%
%	OBJECTINFO (the default) which supports
%	  2-D ZOOMING to the area within the rubber box, double click for reset
%         INFO about objects. Selection type normal (data sets, poles, residues
%	    nodes), Selection type extended (deformations),
%	  CURRENT POLE update (updates residue and local Nyquist plots, see
%	    idiplot and iiplot)
%	  KEY modification of the current plot when a key is pressed. Lower
%	    and upper case keys have opposite effects. 'xyz' translations,
%	    'uvw' rotations, 'a' global shrink, 'i' nominal zoom
%	    For FEPLOT axes : 's' deformation scale.
%	    For other axes  : 'rst' xyz limits shrink/expand
%
%	MOVEAXIS can be used to move axes around.
%	POSITION
%
%	See also IICOM,  IIPLOT, IDIPLOT

%       Etienne Balmes  02/02/94, 06/24/96
%       Copyright (c) 1990-2004 by SDTools,
%       All Rights Reserved.

global IIgui IIAxes

if nargin==0 CAM='on';  else     co=gco;if isempty(co) return;end; end
Cam=lower(CAM);
if any(size(IIgui)<[10 7]) IIgui(10,7)=0;end
IIgui(1,5)=IIgui(1,5)+1;
if IIgui(3,5)>100 disp(['iimouse(''' Cam ''')']); end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if comstr(Cam,'on')
  if comstr(version,'5') on='on'; set(gcf,'busyaction','cancel');
  else on='on';end
  set(gcf,'WindowButtonDown','iimouse(''ObjectInfo'')', ...
      'KeyPressFcn','iimouse(''key'');','interrup',on)
  IIgui(1,5)=0;
elseif comstr(Cam,'off')
  set(gcf,'WindowButtonDown','','KeyPressFcn','')
  IIgui(1,5)=0;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'objectinfo')

% go up to current [root figure axes ...]
eval('coty=get(co,''type'');','co=gca;coty=''axes'';')
aa=co;while ~strcmp(get(aa(1),'type'),'root') aa=[get(aa(1),'parent') aa]; end
if size(aa)<3 aa=[0 gcf gca]; end
if size(aa)<4 aa(4)=0; end
if aa(3)==0 aa(3)=gca; end; if aa(2)==0 aa(2)=gcf;end
axes(aa(3));

sel=get(aa(2),'selectiontype');
% safe operation of open selection (for MACs)
if IIgui(1,5)>1 & ~strcmp(sel,'open') sel='open';
elseif IIgui(1,5)>1 IIgui(1,5)=IIgui(1,5)-1;return; end

% determine figure units and selection area

if comstr(version,'5') on='on'; else on='on';end
ab=get(gcf,'units');set(gcf,'units','pixels','interr',on);
pf=get(gcf,'currentpoint');pa=get(gca,'currentpoint');
i1=get(gcf,'position');

rbbox([pf 0 0],pf); drawnow
pf=[pf;get(gcf,'currentpoint')];pa=[pa;get(gca,'currentpoint')];
set(gcf,'units','inches');
i1=[[1;2.54;72]*get(gcf,'position');0 0 1 1;i1];i1=i1([1 2 4 3 5],:);
set(gcf,'units',ab); pf = pf./i1([5 5],[3 4]);

% determine figure units and selection area

st = '|inches|centimeters|normalized|points|pixels|';
if aa(4)~=0 uo = get(aa(4),'userdata');else uo=0;end
if size(uo,2)<8 uo(1,8)=0; end
ua = get(aa(3),'userdata');
if comstr(version,'5') eval('if isstruct(ua) ua=0;end',''); end
if size(ua,2)<9 ua(1,9)=0; end
if size(IIAxes,2)>4 l1=find(IIAxes(:,4)==gca);else l1=[];end

%-----------------------------------------------------------------------
% ObjectInfo section for a triax axis
if     all(ua(1,1:2)==[Inf 4])
   iimouse('moveaxis'); IIgui(1,5)=IIgui(1,5)-1;return;
elseif     all(ua(1)==[1.1])
   iimouse('moveaxis'); IIgui(1,5)=IIgui(1,5)-1;return;

%-----------------------------------------------------------------------
% Axis zooming
elseif max(abs(pf(1,:)-pf(2,:)))>.01 | strcmp(sel,'open')

  if all(ua==0) bb=axis;if size(bb)<6 bb(5:6)=[NaN NaN];end
                ua=[Inf 0 0 bb];set(gca,'userdata',ua);
  end

  %zoom resseting
  if strcmp(sel,'open')
    if all(ua(1,1:2)==[Inf 3])
      if any(ua(3)==[1 2 3 4 7])
         fecom('centerscale');comgui('@Center Scale');
      end
    else
        if ua(1,4)>=ua(1,5) ua(4:5)=[NaN NaN];end
        if ua(1,6)>=ua(1,7) ua(6:7)=[NaN NaN];end
        if ua(1,8)>=ua(1,9) ua(8:9)=[NaN NaN];end
	if ~isnan(ua(1,4))&~isnan(ua(1,5)) set(gca,'xlim',ua(1,4:5));
	else set(gca,'xlimmode','auto');end
	if ~isnan(ua(1,6))&~isnan(ua(1,7)) set(gca,'ylim',ua(1,6:7));
	else set(gca,'ylimmode','auto');end
	if ~isnan(ua(1,8))&~isnan(ua(1,9)) set(gca,'zlim',ua(1,8:9));
	else set(gca,'zlimmode','auto');end
    end
  else
    set(gca,'xlim',[min(pa(:,1)) max(pa(:,1))], ...
            'ylim',[min(pa(:,2)) max(pa(:,2))], ...
	    'zlim',[min(pa(:,3)) max(pa(:,3))])
  end

%-----------------------------------------------------------------------
% ObjectInfo section for a deformed structure axis
elseif all(ua(1,1:2)==[Inf 3])

 if isempty(l1) disp('Not a current feplot axis');
     IIgui(1,5)=IIgui(1,5)-1;return;
 end
 if strcmp(sel,'normal') % find the closest node

        Xform=fecom('viewmtx',IIAxes(l1,:));
        eval('set(IIgui(6,6),''xform'',Xform)','');
	us=get(IIgui(6,8),'userdata'); node=get(us(3),'userdata');
        i4=[node(:,5:7) ones(size(node,1),1)]*Xform';
        i1=sum(((i4(:,1:2)-ones(size(i4,1),1)*pa(1,1:2)).^2)');
        i3=find(isfinite(i1));[i1,i2]=sort(i1(i3));i1=[i1 Inf];i2=i2(1); 
        if i1(1)<i1(length(i1)-1)/50
          set(gca,'zlimmode','auto');i6=get(gca,'zlim');
          if comstr(version,'5') sv = 'marker';else sv='linestyle'; end
          set(IIAxes(l1,29),'xdata',i4(i3(i2),1),'ydata',i4(i3(i2),2), ...
            'zdata',i6(2), ...
            sv,'o','color',[1 0 0],'clipping','off','visible','on');
          comgui(sprintf('@Node %i @ %.3g,%.3g,%.3g',node(i3(i2),[1 5:7])))
        else set(IIAxes(l1,29),'visible','off');end

 elseif isfinite(uo(1))
 elseif any(uo(1,2)==[1:4 7]) % currently no object info

   i1=get(gco,'xdata');i2=get(gco,'ydata');i3=get(gco,'zdata');
   i1=[i1(:) i2(:) i3(:)];

   i3=(i1(:,1)-pa(1,1)).^2+(i1(:,2)-pa(1,2)).^2;i4=find(isfinite(i3));
   [i6,i3]=min(i3(i4)); % find the closest point
   i5=axis;i5=i5(2)-i5(1)+i5(4)-i5(3);

   if i6(1)<.05*i5
     us=get(IIgui(6,8),'userdata'); node=get(us(3),'userdata');
     i4=i4(i3(1));i3=i1(i4,1:3);set(gca,'zlimmode','auto');i6=get(gca,'zlim');
     if comstr(version,'5') sv = 'marker';else sv='linestyle'; end
     set(IIAxes(l1,29),'xdata',i3(1),'ydata',i3(2),'zdata',i6(2), ...
       sv,'o','color',[1 0 0],'clipping','off','visible','on');
     if uo(1,2)==4
       i3=node(uo(1,i4+6),:);
       comgui(sprintf('@Node %i @ %.3g,%.3g,%.3g',i3(1),i3(5:7)))
     elseif any(uo(1,2)==[3 7])
	Sens=get(us(9),'userdata');Sens=Sens(uo(4:length(uo)),:);
        if uo(1,2)==3 Sens=Sens(fix(i4/3)+1,:);
        else Sens=Sens(fix(i4/6)+1,:);end
        comgui(sprintf('@Sensor %g @ node %i',Sens(1),node(Sens(2),1)));
     elseif any(uo(1,2)==[1 2])
	Def=get(us(10),'userdata');i4=uo(1,i4+6);
        i5=size(node,1);i5=[i4 i4+i5 i4+2*i5];
       comgui(sprintf('@Node %i Def %.3g,%.3g,%.3g',node(i4,1),Def(i5,uo(1,3))))

     end 
   else comgui('@Not a Node');
      if IIAxes(l1,29)~=0 eval('set(IIAxes(l1,29),''visible'',''off'')','');end
   end
 elseif uo(1,2)==6 % this is a patch

  comgui(sprintf('@Group %i elt %i',uo(4:5)));
  %first row nPtr nCol size(elt,1)
  %[n1 n2 n3 n4 EltN GroupN pt1 ... ptN nCol1 ... nColn];

 end

% legend axes - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif any([100:110]==ua(1,3))
   iimouse('moveaxis'); IIgui(1,5)=IIgui(1,5)-1;return;
%-----------------------------------------------------------------------
% ObjectInfo section for IIPLOT, IDIPLOT
elseif ~isfinite(ua(1)) & any(ua(2)==[1 2])

if isempty(l1)
  IIgui(1,5)=IIgui(1,5)-1;disp('Not a current iiplot axis');return;
end
if isfinite(uo(1,1))
  i1=find(IIAxes(:,1)==1|IIAxes(:,1)==2);
  set(IIAxes(i1,[10 19]),'visible','off');

elseif uo(1,2)==99 % pole line - - - - - - - - - - - - - - - - - - -

    po=ii_pof(uo(3:size(uo,1),1:2),abs(uo(1,3)));
    set(IIAxes(l1,10),'units','data','fontsize',9,'visible','on');

    if uo(1,3)==2 & ua(2)==1 % pole line from iiplot
     [i1,i2]=min(abs(po(:,2)-pa(2,1)));
     set(IIAxes(l1,10),'position',[po(i2,2) pa(2,2) 0], ...
         'string',sprintf('%10.3e + i %10.3e',po(i2,:)));

    elseif uo(1,3)==2 & ua(2)==2 % real/imag from idiplot

     [i1,i2]=min(abs(po(:,2)-pa(2,2))/max(po(:,2))+abs(po(:,1)-pa(2,1)));
     set(IIAxes(l1,10),'position',[po(i2,1) po(i2,2) 0], ...
         'string',sprintf('%10.3e + i %10.3e',po(i2,:)));

    elseif uo(1,3)==3
     [i1,i2]=min(abs(po(:,1)-pa(2,1))/max(po(:,1))+abs(po(:,2)-pa(2,2)));
     set(IIAxes(l1,10),'position',[po(i2,1) po(i2,2) 0], ...
         'string',sprintf('%10.3e + i %10.3e',po(i2,:)))

    elseif uo(1,3)==-4
     po=ii_pof(uo(3:size(uo,1),1:2),3);
     i1=find(uo(3:size(uo,1),3)==round(pa(2,2)));
     [i3,i2]=min(abs(po(i1,1)-pa(2,1)));i2=i1(i2);
     set(IIAxes(l1,10),'position',[po(i2,1) uo(i2+2,3) 0], ...
         'string',sprintf('%10.3e + i %10.3e',po(i2,:)))

     global IIpo1; IIpo1=po(i2,:);
    else
     [i1,i2]=min(abs(po(:,1)-pa(2,1)));
     set(IIAxes(l1,10),'position',[po(i2,1) pa(2,2) 0], ...
         'string',sprintf('w %10.3e z %10.3e',po(i2,:)))
    end

    i1=deblank(comstr(setstr(uo(2,:)),1));
    if isempty(i1)
     comgui(sprintf('@Pole %i w %g z %g',uo(i2+2,3),ii_pof(po(i2,:),3)));
    else
     comgui(sprintf('@%s(%i,:) at w %.5g z %.4g',i1,uo(i2+2,3), ...
            ii_pof(po(i2,:),3)));

     if strcmp(i1,'IIpo')|strcmp(i1,'IIpo1') % update modeshape & LNY plots

       i1=find(IIAxes(:,1)==2&IIAxes(:,20)==1&IIAxes(:,15)==0&IIAxes(:,2)>51);
       if ~isempty(i1) IIAxes(i1,21)=ones(size(i1))*uo(i2+2,3);
          idiplot('plot',i1);
       end
       i1=find(IIAxes(:,1)==1&IIAxes(:,2)==10);
       if ~isempty(i1) % update lny plots
          IIgui(3,1)=uo(i2+2,3);iiplot('plot',i1);
       end 
     end % updates if IIpo
    end

% data line - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif any([1:10]==uo(1,2))

    x=get(co,'xdata');y=get(co,'ydata');
    i1 = axis;i1=[i1(2)-i1(1) i1(4)-i1(3)];
    [i1,i2]=min(abs(x-pa(1,1))/i1(1)+abs(y-pa(1,2))/i1(2));

    if comstr(version,'5') sv = 'marker';else sv='linestyle'; end
    set(IIAxes(l1,19),'xdata',x(i2),'ydata',y(i2),'zdata',1, ...
        sv,'+','markersize',6,'color',[1 0 0],'visible','on');
    i1=[find(~uo(:));length(uo)+1];st=setstr(uo(3:i1(1)-1));

    comgui(sprintf('@%s(%i,%i) at x %7.3e, y %7.3e ', ...
                   st,i2+uo(i1(1)+1)-1,uo(i1(1)+3),x(i2),y(i2)))

elseif any([0]==uo(1,2))
    comgui(['@' setstr(uo(3:length(uo)))]);
end % of iiplot objects

%-----------------------------------------------------------------------
% not a IIGUI axis
else

%-----------------------------------------------------------------------
end %of type of axis sections

if ~isfinite(ua(1,1)) & any([1 2 3]==ua(1,2))
  i1=find(IIAxes(:,4)==aa(3));
  if ~isempty(i1) IIgui(3,7)=i1(1);comgui('refresh');end
end

%-----------------------------------------------------------------------
elseif comstr(Cam,'moveaxis')

[CAM,Cam] = comstr(CAM,9);
if strcmp(Cam,'ch')

   set(gca,'position',[get(gcf,'currentpoint') 0 0]+ ...
                       get(gca,'position')*diag([0 0 1 1]));
   if strcmp(get(gcf,'SelectionType'),'open') iimouse('moveaxisoff'); end

elseif strcmp(Cam,'in') | isempty(Cam)

   comgui('@moveaxis on');
   set(gcf,'windowbuttonmotionfcn','iimouse(''MoveAxisCh'');')
   set(gcf,'windowbuttondownfcn','iimouse(''MoveAxisoff'');')
   set(gcf,'units',get(gca,'units'));

elseif strcmp(Cam,'off')
   set(gcf,'WindowButtonDown','iimouse(''ObjectInfo'')')
   set(gcf,'windowbuttonmotionfcn','')
   comgui('@moveaxis off');
end

%-----------------------------------------------------------------------
elseif comstr(Cam,'position') [CAM,Cam] = comstr(CAM,9);

eval('figure(get(IIAxes(IIgui(3,7),4),''parent''))','');
if isempty(CAM)

  set(gcf,'WindowButtonDown','iimouse(''PositionSet'')')
  disp('Indicate new position of current axis');

elseif comstr(Cam,'set')

 if comstr(version,'5') on='on'; else on='on';end
 ab=get(gcf,'units');set(gcf,'units','pixels','interr',on);
 pf=get(gcf,'currentpoint');i1=get(gcf,'position');
 rbbox([pf 0 0],pf);
 pf=[pf;get(gcf,'currentpoint')];
 set(gcf,'units',ab);

 if strcmp(get(gcf,'SelectionType'),'extend')
   set(gcf,'WindowButtonDown','iimouse(''ObjectInfo'')');
   i1=find(IIAxes(:,4)==gca);
   if ~isempty(i1) IIAxes(i1(1),6:9)=get(gca,'position'); end
   disp('New position accepted');
 elseif comstr(get(gcf,'SelectionType'),'alt')
   set(gcf,'WindowButtonDown','iimouse(''ObjectInfo'')');
   i1=find(IIAxes(:,4)==gca);
   if ~isempty(i1) set(gca,'position',IIAxes(i1(1),6:9));
           disp('New position rejected');
   else disp('Can''t reject new position of a non IIGUI axis');end

 else
  eval('axes(IIAxes(IIgui(3,7),4))','');ab=get(gca,'units');
  aa=[min(pf) abs(pf(1,:)-pf(2,:))];
  if all(aa(3:4)>0) set(gca,'units','pixels','position',aa); end
  set(gca,'units',ab);
  disp('Extend to accept new position, Alternate to reject');
 end
end

% ----------------------------------------------------------------------------
elseif comstr(Cam,'key')

char = get(gcf, 'CurrentCharacter'); dx = [zeros(1,3);ones(1,3);zeros(1,3)];
if isempty(char); char=' ';end

if ~isempty(IIAxes)
 i1=find(IIAxes(:,3)==gcf);
 if isempty(i1==IIgui(3,7)) IIgui(3,7)=i1(1); end
 l1=IIgui(3,7);axes(IIAxes(l1,4));
else l1=[]; IIAxes=zeros(1,30);end

ua=get(gca,'userdata');
if isempty(ua)|all(ua==0) bb=axis;if size(bb)<6 bb(5:6)=[NaN NaN];end
                ua=[Inf 0 0 bb];set(gca,'userdata',ua);
end
axin = [get(gca,'xlim') get(gca,'ylim') get(gca,'zlim')];
i1 = axin(2:2:6)-axin(1:2:5);
i1(2,1)=norm(i1)/10;

if     any(char == '+=') iicom('ch+');
elseif any(char == '-_') iicom('ch-');
elseif IIAxes(l1,1)==3 % FEPLOT axes

if     char == 'x' IIAxes(l1,15)=IIAxes(l1,15)+i1(2,1); fecom('center');
elseif char == 'X' IIAxes(l1,15)=IIAxes(l1,15)-i1(2,1); fecom('center');
elseif char == 'y' IIAxes(l1,16)=IIAxes(l1,16)+i1(2,1); fecom('center');
elseif char == 'Y' IIAxes(l1,16)=IIAxes(l1,16)-i1(2,1); fecom('center');
elseif char == 'z' IIAxes(l1,17)=IIAxes(l1,17)+i1(2,1); fecom('center');
elseif char == 'Z' IIAxes(l1,17)=IIAxes(l1,17)-i1(2,1); fecom('center');
elseif char == 'l' IIAxes(l1,23)=IIAxes(l1,23)/2;feplot('plot');
elseif char == 'L' IIAxes(l1,23)=IIAxes(l1,23)*2;feplot('plot');
elseif char == 'a' IIAxes(l1,18)=IIAxes(l1,18)/1.2;fecom('center');
elseif char == 'A' IIAxes(l1,18)=IIAxes(l1,18)*1.2;fecom('center');
elseif char == 'i' IIAxes(l1,11:13)=IIgui(9,1:3);
                   fecom(';centerscale;centernode');
elseif char == 'u' IIAxes(l1,11)=rem(IIAxes(l1,11)+10,360);
                   feplot('plotnoscale');
elseif char == 'U' IIAxes(l1,11)=rem(IIAxes(l1,11)-10,360);
                   feplot('plotnoscale');
elseif char == 'v' IIAxes(l1,12)=rem(IIAxes(l1,12)+10,360);
                   feplot('plotnoscale');
elseif char == 'V' IIAxes(l1,12)=rem(IIAxes(l1,12)-10,360);
                   feplot('plotnoscale');
elseif char == 'w' IIAxes(l1,13)=rem(IIAxes(l1,13)+10,360);
                   feplot('plotnoscale');
elseif char == 'W' IIAxes(l1,13)=rem(IIAxes(l1,13)-10,360);
                   feplot('plotnoscale');
end
fecom('setslider')

else % Other axes

if     char == 'u' dx(1,1)=-10;
elseif char == 'U' dx(1,1)=+10;
elseif char == 'v' dx(1,2)=-10;
elseif char == 'V' dx(1,2)=+10;
elseif char == 'A' dx(2,1:3) = 1/1.2*[1 1 1];
elseif char == 'a' dx(2,1:3) = 1.2*[1 1 1];
elseif char == 'x' dx(2,1) = .9;    elseif char == 'X'	dx(2,1) = 1.1; 
elseif char == 'y' dx(2,2) = 0.9;   elseif char == 'Y'	dx(2,2) = 1.1; 
elseif char == 'z' dx(2,3) = 0.9;   elseif char == 'Z'	dx(2,3) = 1.1; 
elseif char == 'r' dx(3,1) = 0.1;   elseif char == 'R'	dx(3,1) = -.1; 
elseif char == 's' dx(3,2) = 0.1;   elseif char == 'S'	dx(3,2) = -.1; 
elseif char == 't' dx(3,3) = 0.1;   elseif char == 'T'	dx(3,3) = -.1; 
end

if char == 'i'
	ua=get(gca,'userdata'); if size(ua,2)<9 ua(1,4:9)=ones(1,6)*NaN;end
        if ua(1,4)>=ua(1,5) ua(4:5)=[NaN NaN];end
        if ua(1,6)>=ua(1,7) ua(6:7)=[NaN NaN];end
        if ua(1,8)>=ua(1,9) ua(8:9)=[NaN NaN];end
	if ~isnan(ua(1,4))&~isnan(ua(1,5)) set(gca,'xlim',ua(1,4:5));
	else set(gca,'xlimmode','auto');end
	if ~isnan(ua(1,6))&~isnan(ua(1,7)) set(gca,'ylim',ua(1,6:7));
	else set(gca,'ylimmode','auto'); end
	if ~isnan(ua(1,8))&~isnan(ua(1,9)) set(gca,'zlim',ua(1,8:9));
	else set(gca,'zlimmode','auto');end
        if ~isempty(l1) IIAxes(l1,11:12)=[0 90];view([0 90]);end

else
    i2 = [(axin(2:2:6)+axin(1:2:5))/2;(axin(2:2:6)-axin(1:2:5))/2];
    i2(2,:)=i2(2,:).*dx(2,:);i2(1,:)=i2(1,:)+dx(3,:).*i2(2,:);
    i3=find(i2(2,:)<eps*abs(i2(1,:)));
    if ~isempty(i3) i2(2,i3)=i2(2,i3)+eps*abs(i2(1,:));end
    lim=reshape([i2(1,:)-i2(2,:);i2(1,:)+i2(2,:)],2,3);
    set(gca,'xlim',lim(1:2),'ylim',lim(3:4),'zlim',lim(5:6))
    if any(dx(1,:))
      if isempty(l1) i2=get(gca,'view');else i2=IIAxes(l1,11:12);end
      view(i2+dx(1,1:2)); if ~isempty(l1) IIAxes(l1,11:12)=i2+dx(1,1:2);end
    end
   set(gca,'xlimmode','manual','ylimmode','manual','zlimmode','manual')
end %~='i'
end %other axes

IIAxes=IIAxes(find(IIAxes(:,1)~=0),:);

%-----------------------------------------------------------------------
elseif comstr(Cam,'feplot') [CAM,Cam] = comstr(CAM,7);

if strcmp(Cam,'ch')

   i1 = get(gcf,'CurrentCharacter');

if isempty(i1)
elseif i1=='q' iimouse('feplot off');
elseif i1=='r'
   i2 = (get(gcf,'CurrentPoint')-IIgui(9,1:2))*180+IIgui(9,3:4);
   fecom(sprintf('view(%f,%f)',i2));
end

elseif strcmp(Cam,'in') | isempty(Cam)

   set(gcf,'units','normalized');
   IIgui(9,1:4)=[get(gcf,'currentpoint') get(gca,'view')];
   comgui('@FEMouse on');
   set(gcf,'windowbuttonmotionfcn','iimouse(''FeplotCh'');', ...)
           'windowbuttonupfcn','iimouse(''FeplotOff'');','interruptible','no')
   i1 = 26+[1:IIAxes(IIgui(3,7),30)]*5;
   set(IIAxes(IIgui(3,7),i1+4),'visible','off');

elseif strcmp(Cam,'off')
   set(gcf,'WindowButtonDown','iimouse(''ObjectInfo'')')
   set(gcf,'windowbuttonupfcn','')
   set(gcf,'windowbuttonmotionfcn','')
   comgui('@FEMouse off');
   i1 = 26+[1:IIAxes(IIgui(3,7),30)]*5;
   set(IIAxes(IIgui(3,7),i1+4),'visible','on');

end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
else  error('Not a valid mode for IIMOUSE')
end

IIgui(1,5)=IIgui(1,5)-1;
