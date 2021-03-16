function [out0] = fecom(CAM,CAM1,CAM2,CAM3,CAM4,CAM5)

%FECOM	UI command functions for deformed structure visualizations
%
%	Synopsis : fecom
%		   fecom command
%		   fecom('command','tail')
% 
%	Commands available through the command mode are
%
%	  Anim [, start, stop]
%	  cax [i], ca+	   : changes the current axis number (see IICOM)
%	  Ch [+,-,(i)]	   : changes the current deformation displayed
%	  Center [Node (i), Pos (x y z)] : will now put node (i) or position
%			    (x y z) in the middle of the figure.
%	  CenterNode	   : gives the default centering
%	  Center Scale (i) : for (i)=1 the structure fits in the axis
%			     other values make it smaller or bigger
%	  Center	   : centers the current axis
%	  ColorData [Order, Element, Node, Group, Lighting, ElNode]
%			   : defines how patch colordata is defined
%	  Color[Face,Edge] [None,Interp,Flat,White,Black,Order]
%			   : sets default patch edge and face color properties
%	  Group [(i), All, Rem (i), Add (i)]
%			   : selects group displayed in the current figure
%	  Head [,freq]     : titles for each deformation (call with arguments)
%	  LineOptim        : optimizes line segmentation for wire-frame lines
%         NodeText[, Off, (i)] : displays all node numbers, no node numers or
%			     node numbers (i). You can use feplot('nodetext',i)
%	  DofText [ ,(i), Off] : displays all DOF info, no DOF info or
%			     info about DOFs (i). Use also feplot('doftext',i)
%	  Scale[ , Def [i], Dof [i], equal, auto]
%			   : set procedure used to scale the deformation
%	  Show[line,patch,sens,arrow] : reinit the current axis with PlotType
%	  Sub [i j]	   : creation of more than one plot axis (see IICOM)
%	  Triax, TriaxOff  : show/not show orientation triax
%	  Titopt	   : figure title type (0 none, 1 title)
%	  View(az el selr) : sets the view angles for the current plot
%	  Undef [ , Dot, Line] : modify the way the undeformed structure is
%			     shown. 
%
%	Commands are insensitive to caps, they can be chained using COMMODE
% 	or simply by starting a command by a semicolumn ';'
%	For the initialisation of a plot use FEPLOT

%	The following are internal commands
%	InitAxes, Key, SetSlider

%	Etienne Balmes 05/11/92, 01/23/96
%       Copyright (c) 1990-2004 by SDTools,
%       All Rights Reserved.

global IIgui IIAxes IIMenu

if nargin==0 commode('fecom'); return;
elseif isstr(CAM)  [CAM,Cam]=comstr(CAM,1); end
CAM0=CAM; out='done';
if nargin<2 CAM1=[]; end

if IIgui(3,5)>100 disp(['fecom(''' Cam ''')']); end

yplot=0;
if isempty(IIAxes) l1=[];
else
  l1 = find(IIAxes(:,1)==3);
  if comstr(Cam,'plot') feplot(CAM,CAM1); if nargout==1 out0='done';end;return;
  elseif isempty(l1) disp('WARNING no existing feplot axis');
  elseif ~isempty(find(l1==IIgui(3,7))) l1=IIgui(3,7);
  else l1=l1(1);IIgui(3,7)=l1; end
end

% ----------------------------------------------------------------------------
if comstr(Cam,';')  commode('fecom',CAM(2:length(CAM)));

% ----------------------------------------------------------------------------
elseif comstr(Cam,'anim') [CAM,Cam]=comstr(CAM,5); % Modeshape animation

 ind=[];j3=1;% set amplification vector ind and index j3

 if comstr(Cam,'stop')
    set(IIgui(6,5),'value',0,'callback','fecom(''animstart'');'); return;
 elseif comstr(Cam,'start')
    set(IIgui(6,5),'value',1,'callback','fecom(''animstop'');');
 elseif get(IIgui(6,5),'value')==1
   if strcmp(get(IIgui(6,5),'callback'),'fecom(''animstop'');')
     set(IIgui(6,5),'value',0,'callback','fecom(''animstart'');'); return;
   else
    set(IIgui(6,5),'callback','fecom(''animstop'');');
    if ~isempty(CAM1) ind=CAM1(:)';
     ind=ind([2:length(ind) 1])./ind;
     if any(~isfinite(ind)) error('Cannot animate with null deformation');end
     IIAxes(l1,21)=length(ind);
    end
   end
 end

 if IIgui(6,8)==0 return; else us=get(IIgui(6,8),'userdata'); end
 node=get(us(3),'userdata');
 Def=get(us(10),'userdata');Sens=get(us(9),'userdata');
 nInf=size(node,1);NNode(node(:,1))=[1:length(node(:,1))]';

 ua=get(IIAxes(l1,4),'userdata');
 set(gca,'XLimMode','manual','YLimMode','manual','ZLimMode','manual',...
     'drawmode','fast')
 set(IIgui(6,5),'backgroundcolor',[1 0 0]);

while get(IIgui(6,5),'value') % this is the infinite loop - - - - - - - - - -

  opt=clock;  
  i1=IIAxes(l1,11:13)*pi/180;i1=[cos(i1);sin(i1)]';
  Xform=[i1(3,1) i1(3,2) 0 0;-i1(3,2) i1(3,1) 0 0;0 0 1 0;0 0 0 1] * ...
      [i1(1,1) i1(1,2) 0 0;-i1(2,2)*i1(1,2) i1(2,2)*i1(1,1) i1(2,1) 0;
       i1(2,1)*i1(1,2) -i1(2,1)*i1(1,1) i1(2,2) 0; 0 0 0 1];
  %Xform=fecom('viewmtx',IIAxes(l1,:));

if any(ua(2:size(ua,1),1)==1)  feplot('PlotNoscale');

else %plots without patches

i1=2:size(ua,1);i1=find((ua(i1,1)==3|ua(i1,1)==2|ua(i1,1)==7)&ua(i1,4)~=0)+1;

% loop on objects of the current plot
for j2=i1(:)' jDef=26+ua(j2,4)*5;

% update a deformed line object - - - - - - - - - - - - - - - - - -
if ua(j2,1)==2

  uo=get(ua(j2,3),'userdata');  uo(1,3)=IIAxes(l1,jDef);
  LD2=uo(7:length(uo));
  aa=[node(LD2,5:7)+real([Def(LD2,uo(1,3)) ...
       Def(LD2+nInf,uo(1,3)) Def(LD2+2*nInf,uo(1,3))]* ...
       IIAxes(l1,jDef+1)*exp(i*IIAxes(l1,jDef+2))) ones(length(LD2),1)] ...
       *Xform';
  set(ua(j2,3),'xdata',aa(:,1),'ydata',aa(:,2), ...
      'zdata',aa(:,3),'userdata',uo);

% Anim a sensor line object - - - - - - - - - - - - - - - - - -
elseif any(ua(j2,1)==[3 7])

  uo=get(ua(j2,3),'userdata'); uo(1,3)=IIAxes(l1,jDef);
  LD2=uo(4:length(uo));
   if ~isempty(LD2)
     i2=Sens(LD2,:);
     LD2=i2(:,2);LD2=[LD2 LD2 nInf*ones(size(LD2,1),1)]'; LD2=LD2(:);
     a0=real(IIAxes(l1,jDef+1)*exp(i*IIAxes(l1,jDef+2))* ...
        (i2(:,3).*Def(i2(:,2),IIAxes(l1,jDef))+ ...
        i2(:,4).*Def(i2(:,2)+nInf,IIAxes(l1,jDef))+ ...
        i2(:,5).*Def(i2(:,2)+2*nInf,IIAxes(l1,jDef))));
     a0=a0(:,[1 1 1]).*i2(:,3:5);
     aa=[node(LD2,5:7) ones(length(LD2),1)];
     aa(2:3:size(aa,1),1:3)=aa(2:3:size(aa,1),1:3)+a0;
     aa=aa*Xform';
     if ua(j2,1)==7 % make the arrow
      i2=2:3:size(aa,1); aa=[aa(i2,1:2) aa(i2,1:2)-aa(i2-1,1:2)];
      aa(:,5:6)=aa(:,[3:4])*[0 -1;1 0]/20;
      aa=[aa(:,1:2)-aa(:,3:4) aa(:,1:2) aa(:,1:2)-aa(:,3:4)/10+aa(:,5:6) ...
        aa(:,1:2)-aa(:,3:4)/10-aa(:,5:6) aa(:,1:2) aa(:,1:2)*Inf];%
      i2=size(aa,1);
      aa=[reshape(aa(:,1:2:size(aa,2))',i2*size(aa,2)/2,1) ...
        reshape(aa(:,2:2:size(aa,2))',i2*size(aa,2)/2,1) ...
        ones(i2*size(aa,2)/2,1)];
     end % make the arrow

     set(ua(j2,3),'xdata',aa(:,1),'ydata',aa(:,2),'zdata',aa(:,3));
   end

end      % of object types - - - - - - - - - - - - - - - - - -
end;  % j2 loop on objects of the current plot

end % of plots with/without patches 

drawnow;
%while etime(clock,opt)<.08 end
%disp(etime(clock,opt))

% change phase/amplitude
j1=26+[1:IIAxes(l1,30)]*5;
if length(ind)~=IIAxes(l1,21)
  ind=exp(i*[0:IIAxes(l1,21)-1]*2*pi/IIAxes(l1,21));
  ind=ind([2:length(ind) 1])./ind;
end
aa=IIAxes(l1,j1+1).*exp(i*IIAxes(l1,j1+2))*ind(j3);
j3=remi(j3+1,IIAxes(l1,21));% increment index count
IIAxes(l1,j1+2)=atan2(imag(aa),real(aa)); IIAxes(l1,j1+1)=abs(aa);

if IIgui(3,5)==99 %save each as a PS file
  eval(sprintf('print -deps def%i',j3)); if j3==1 break;end
end

end % of loop waiting for button to be pressed up
 set(IIgui(6,5),'backgroundcolor',[1 1 1]*.706);

% ----------------------------------------------------------------------------
elseif comstr(Cam,'ca') iicom(CAM,CAM1);
% ----------------------------------------------------------------------------
elseif comstr(Cam,'clear')

 eval('delete(findobj(''tag'',''Struct1''));','');

% ----------------------------------------------------------------------------
elseif comstr(Cam,'sub') iicom(CAM,CAM1);
% ----------------------------------------------------------------------------
elseif comstr(Cam,'center') [CAM,Cam]=comstr(CAM,7);

% get current structure info
if IIgui(6,8)==0 return; else us=get(IIgui(6,8),'userdata'); end
node=get(us(3),'userdata');nInf=size(node,1);
Line=get(us(6),'userdata');NNode(node(:,1))=[1:length(node(:,1))]';

LD2=Line(IIAxes(l1,30+IIAxes(l1,30)*5+[1:IIAxes(l1,25)]),2:size(Line,2))+1;
LD2=LD2(:);LD2=find(sparse(LD2,1,LD2))-1;
if LD2(1)==0 LD2=LD2(2:length(LD2));end
if LD2(length(LD2))==nInf LD2=LD2(1:length(LD2)-1);end

if comstr(Cam,'scale') i1=comstr(comstr(CAM,6),[-1 1]); IIAxes(l1,18)=i1(1);end
if comstr(Cam,'node') % CenterNode

  %one particular node
  i1 = comstr(comstr(CAM,5),[-1 0]);
  i1 = find(node(:,1)==i1(1));
  if ~isempty(i1) if ~all(isfinite(node(i1,5:7))) i1=[]; end;end

  % the mean of selected nodes
  if isempty(i1)
    if isempty(LD2) LD2=1:size(node,1)-1; end
    if ~strcmp(CAM1,';') disp('Centering on mean of selected nodes');end
    i1=node(LD2,5:7);if size(i1,1)>1 i1=mean([min(i1);max(i1)]);end
    IIAxes(l1,15:17) = i1;
  else % a particular node
    disp(sprintf('Centering on node %i',node(i1,1)));
    IIAxes(l1,15:17) = node(i1,5:7);   
  end
  if ~isfinite(IIAxes(l1,18)) IIAxes(l1,18)=1; end

elseif comstr(Cam,'pos')
  opt=comstr(CAM(4:length(CAM)),[-1]);opt=opt(:)';
  if length(opt)~=3 | any(~isfinite(opt))
    error(['fecom length' CAM 'does not specify xyz coordinates']);
  end
  IIAxes(l1,15:17)=opt;
end

% check on the appropriate axis
eval('axes(IIAxes(l1,4));','IIAxes(l1,4)=gca;');
ua = get(gca,'userdata'); if size(ua,2)<2 ua(1,2)=0; end

if all(ua(1,1:2)==[Inf 3]) % this is a feplot axes

  if ~isfinite(IIAxes(l1,18))|IIAxes(l1,18)==0  IIAxes(l1,18)=1; end
  Xform = fecom('viewmtx',IIAxes(l1,:));
  i1 = [IIAxes(l1,15:17) 1]* Xform';
  i2=[node(LD2,5:7) ones(length(LD2),1)]*Xform';
  i2=[min(i2);max(i2)];i2=[mean(i2);(i2(2,:)-i2(1,:))/2];
  lim=[i1-i2(2,:)/IIAxes(l1,18);i1+i2(2,:)/IIAxes(l1,18)];
  i2=find(lim(1,1:3)==lim(2,1:3));
  if ~isempty(i2) lim(2,i2)=lim(2,i2)+1e-5*abs(mean(lim(:,i2))+eps);end

  set(gca,'xlim',lim(:,1)','ylim',lim(:,2)','zlim',lim(:,3)');

else
  disp('fecom center: current axis is not a feplot axis');
end % of all(ua(1,1:2)==[Inf 3])

% ----------------------------------------------------------------------------
elseif comstr(Cam,'ch')

  if isempty(IIAxes) disp('fecom(''ch'') : No axis defined');return; end
  yplot=0;if ~any(1:size(IIAxes,1)==IIgui(3,7)) IIgui(3,7)=1;end
  i1=IIAxes(IIgui(3,7),[1 3]);

  if comstr(Cam,'chc')  ind=IIgui(3,7);CAM=CAM(4:length(CAM));
  else ind=find(IIAxes(:,1)==i1(1)&IIAxes(:,3)==i1(2));CAM=CAM(3:length(CAM));
  end

  if IIgui(6,8)==0 return; else us=get(IIgui(6,8),'userdata'); end
  Def=get(us(10),'userdata');CDefE=[];CDefN=[]; 

 for j1=ind(:)'
  if isempty(CAM) disp('No specified channel');
  elseif CAM(1)=='+'
     if length(CAM)>1 i1=comstr(CAM(2:length(CAM)),[-1 1]);else i1=1;end
     i1=IIAxes(j1,30+[1:5:IIAxes(l1,30)*5])+i1;
     if max(i1)>size(Def,2) i1=[]; end
  elseif CAM(1)=='-'
     if length(CAM)>1 i1=comstr(CAM(2:length(CAM)),[-1 1]);else i1=1;end
     i1=IIAxes(j1,30+[1:5:IIAxes(l1,30)*5])-i1;
     if min(i1)<1 i1=[]; end
  else
     i1 = comstr(CAM,-1); i1 = i1(find(i1<=size(Def,2)));
  end

  i3=IIAxes(j1,30+IIAxes(j1,30)*5+[1:IIAxes(j1,25)]); % group selected
  if ~isempty(i1)
    if length(i1)==IIAxes(j1,30)&yplot~=4 yplot=2; else yplot=4; end
    IIAxes(j1,30+[0 1:5:length(i1)*5]) = [length(i1) i1(:)'];
    IIAxes(j1,30+IIAxes(j1,30)*5+[1:IIAxes(j1,25)])=i3;
    comgui(['@Deformation' sprintf(' %i',IIAxes(j1,26+[1:IIAxes(j1,30)]*5))])
  else if yplot==0 yplot=-1;end; end
 end % loop on axes
 if     yplot==-1 comgui('@Not a valid deformation')
 elseif ~isempty(CAM1)&isstr(CAM1)&comstr(CAM1,';')
 elseif yplot==2&length(ind)==1 feplot('plot');
 elseif yplot==2&length(ind)~=1 feplot('plotall');
 elseif yplot==4&length(ind)==1 feplot('plotcreate');
 elseif yplot==4&length(ind)~=1 feplot('plotallcreate');
 end; yplot=0;

% ----------------------------------------------------------------------------
elseif comstr(Cam,'color') [CAM,Cam]=comstr(CAM,6);

if comstr(Cam,'data')  [CAM,Cam]=comstr(CAM,5);
 %   DefColor: 0 Order 1 Element 2 Node 3 Group  4 Lighting 5 element nodes
 %   EdgeColor 00 none 01 interp 02 flat 03 white 04 black 05 order
 %   FaceColor 00 none 10 interp 20 flat 30 white 40 black 50 order

  if     comstr(Cam,'u')   opt=[0 43];
  elseif comstr(Cam,'eln') opt=[5 11];
  elseif comstr(Cam,'el')  opt=[1 24];
  elseif comstr(Cam,'n')   opt=[2 11];
  elseif comstr(Cam,'gr')  opt=[3 24];
  elseif comstr(Cam,'l')   opt=[4 33];
  else opt=comstr(comstr(CAM,10),[-1 0 0]);end
  if ~any([0:5]==opt(1)) opt(1)=0; end % reset default
  IIAxes(l1,[26 24])=opt;

elseif comstr(Cam,'edge')  [CAM,Cam]=comstr(CAM,5);
  i1=[rem(IIAxes(l1,24),10) fix(IIAxes(l1,24)/10)];
  if     comstr(Cam,'f') i1(1)=2;
  elseif comstr(Cam,'i') i1(1)=1;
  elseif comstr(Cam,'n') i1(1)=0;
  elseif comstr(Cam,'w') i1(1)=3;
  elseif comstr(Cam,'k')|comstr(Cam,'b') i1(1)=4;
  elseif comstr(Cam,'o') i1(1)=5;
  else out='unknown';i1(1)=Inf;end
  if isfinite(i1(1)) IIAxes(l1,24)=i1(1)+10*i1(2); end

elseif comstr(Cam,'face') [CAM,Cam]=comstr(CAM,5);

  i1=[rem(IIAxes(l1,24),10) fix(IIAxes(l1,24)/10)];
  if     comstr(Cam,'f') i1(2)=2;
  elseif comstr(Cam,'i') i1(2)=1;
  elseif comstr(Cam,'n') i1(2)=0;
  elseif comstr(Cam,'w') i1(2)=3;
  elseif comstr(Cam,'k')|comstr(Cam,'b') i1(2)=4;
  elseif comstr(Cam,'o') i1(2)=5;
  else out='unknown';i1(2)=Inf;end
  if isfinite(i1(2)) IIAxes(l1,24)=i1(1)+10*i1(2);end

else out='unknown'; end

yplot=2;

% ----------------------------------------------------------------------------
elseif comstr(Cam,'doftext') [CAM,Cam]=comstr(CAM,8);

 % find the current axis
 if IIAxes(l1,1)~=3 error('fecom nodetext: not a feplot axis');end
 eval('axes(IIAxes(l1,4));','');

 if IIgui(6,8)==0 return; else us=get(IIgui(6,8),'userdata'); end
 % find current DOF text objects and delete them
 ua=get(gca,'userdata'); i1=find(ua(2:size(ua,1),1)==6)+1;
 if isempty(i1)
      ua(size(ua,1)+1,1:3)=[6 1 line('visible','off')];i1=size(ua,1);
 else
      uo=get(ua(i1(1),3),'userdata');
      eval('delete(uo(2:size(uo,1),1));','');
 end
 eval('findobj(ua(i1(1),3));','ua(i1(1),3)=line;')

 node=get(us(3),'userdata');Sens = get(us(9),'userdata');
 NNode(node(:,1))=[1:length(node(:,1))]'; nInf=size(node,1);
 Xform=fecom('viewmtx',IIAxes(l1,:));

 if comstr(Cam,'off') opt=[];
   delete(ua(i1(1),3));ua=ua(find([1:size(ua,1)]~=i1(1)),:);
 % create the DOF text strings 
 elseif isempty(Sens) disp('fecom(''doftext'') initialize sensors first');
 else 
  opt=comstr(CAM,-1);
  if isempty(opt) & ~isempty(CAM1) & ~isstr(CAM1) opt=CAM1(:,1); end
  if isempty(opt) opt=Sens(:,1); end
  uo=[];

  aa=[node(Sens(:,2),5:7)+Sens(:,3:5)*IIAxes(l1,23) ...
        ones(size(Sens,1),1)];uo(:,[2:4])=aa(:,1:3);aa=aa*Xform';
  i2=' x y z-x-y-z s a';
  for j2=1:length(opt)
     j1=find(Sens(:,1)==opt(j2));
     if isempty(j1) disp(sprintf('DofText : %g not found',opt(j2)))
     else
      i3=round(rem(Sens(j1,1),1)*100);st=find(i3==[1:3 7:9 99 98]);
      if isempty(st) st='s';else st=i2(st*2+[-1:0]);end
      uo(j2,1)=text(aa(j1,1),aa(j1,2),aa(j1,3), ...
          [num2str(fix(Sens(j1,1))) st],...
          'FontSize',8,'userdata',[Inf 6]);
     end
  end
  comgui('@Remove using DofTextOff');
  set(ua(i1(1),3),'userdata',[Inf 6 0 0;uo]);
 end
 
 set(gca,'userdata',ua);

% ----------------------------------------------------------------------------
elseif comstr(Cam,'erasemode')

[CAM,Cam]=comstr(CAM,10);
if isempty(CAM)
   if     IIgui(8,8)==1 out0 = 'background';
   elseif IIgui(8,8)==2 out0 = 'xor';
   else                 out0 = 'normal'; end
   return;

elseif strfind(Cam,'back') IIgui(8,8)=1; out = 'background';
elseif strfind(Cam,'xor')  IIgui(8,8)=2; out = 'xor';
else			   IIgui(8,8)=0; out = 'normal'; end
yplot=2;comgui(['@EraseMode set to ' out]);

% ----------------------------------------------------------------------------
% Group [i], GroupAdd [i], GroupRem [i], GroupAll [i]
elseif comstr(Cam,'group') [CAM,Cam] = comstr(CAM,6);

 if IIgui(6,8)==0 return; else us=get(IIgui(6,8),'userdata'); end
 group=get(us(5),'userdata');

 if isempty(group)
   st='No declared group cannot plot';comgui(['@' st]); return;
 end

 i1=length(find(group=='|'))-1;

 if isempty(strfind(Cam,'back')) i1(4)=0; else i1(4)=1;end
 if      comstr(Cam,'c')  ind=l1; [CAM,Cam]=comstr(CAM,2); i1(5)=4;
 else                     ind = find(IIAxes(:,1)==3);i1(5)=5; end
 if      comstr(Cam,'all')  i1(3)=0; CAM=[];Cam=[];
 elseif  comstr(Cam,'rem')  i1(3)=1; [CAM,Cam]=comstr(CAM,4);
 elseif  comstr(Cam,'add')  i1(3)=2; [CAM,Cam]=comstr(CAM,4);
 else i1(3)=3; end

 i2=comstr(CAM,-1); %get group numbers
 if isempty(i2)&~isempty(Cam) % seek by group name
    i2=find(group(1:i1)=='|');
    while ~isempty(Cam) [st,i4,i5,i6]=sscanf(Cam,'%s',1);
     if i4==0 Cam=[]; end; Cam=Cam(i6:length(Cam));
     i4=strfind(st,group);
     if ~isempty(i4) i3=[i3 length(find(i2<i4))]; end
    end; i2=i3;
 end 

 for j1= ind(:)'

   if 30+IIAxes(j1,30)*5+IIAxes(j1,25)<size(IIAxes,2)
     i3=IIAxes(j1,30+IIAxes(j1,30)*5+[1:IIAxes(j1,25)]);
   else i3=[]; end

   if     i1(3)==0	i3=1:i1(1); % all
   elseif i1(3)==1 %rem
     for j1=1:length(i2)
       if ~isempty(find(i3==i2(j1))) i3(find(i3==i2(j1)))=0; end
     end; i3=i3(find(i3));
   elseif i1(3)==2 i3 = [i3 i2(:)']; %add
   elseif ~isempty(i2) i3 = i2;
   end
   i3 = i3(find(i3>0 & i3<=i1(1))); if isempty(i3) i3=1:i1(1); end
   IIAxes(j1,25)=length(i3);
   IIAxes(j1,30+IIAxes(j1,30)*5+[1:IIAxes(j1,25)]) =i3;

 end

 if i1(4)==0 yplot=i1(5); end% isempty(strfind(Cam,'back')) 

% get commands ---------------------------------------------------------------
elseif comstr(Cam,'get')  [CAM,Cam] = comstr(CAM,4);

if comstr(Cam,'cax') % - - - - - - - - - - - - - - - - - - - - - -
 eval('out=get(IIAxes(l1,4),''userdata'');','out=[];');
else out='unknown';end

% help commands ---------------------------------------------------------------
elseif comstr(Cam,'head') [CAM,Cam] = comstr(CAM,5);

 up=findobj(2,'type','axes','tag','Struct1');
if isempty(up) disp('No defined structure');
else
  us=get(up,'userdata');
  Def=get(us(10),'userdata');  Head=get(us(15),'userdata');  
  if comstr(CAM,'freq')&size(CAM1,1)==size(Def,2)
   Head=[];for j1=1:size(CAM1,1)
     Head=str2mat(Head,sprintf('Mode %i at %.4g Hz',j1,CAM1(j1,1)));
   end
   Head=Head(2:size(Head,1),:);
  elseif comstr(CAM,'fs')&size(CAM1,1)==size(Def,2)
   Head=[];for j1=1:size(CAM1,1)
     Head=str2mat(Head,sprintf('%i @ %.4g Hz',j1,CAM1(j1,1)));
   end
   Head=Head(2:size(Head,1),:);
  elseif isstr(CAM1)&size(CAM1,1)==size(Def,2)
     Head=CAM1;
  elseif size(CAM1,1)~=size(Def,2)
   disp('Using default deformation labels');
   i1=fix(log10(size(Def,2))+eps)+1;st=sprintf('Def. %%%ii',i1);
   Head = reshape(sprintf(st,1:size(Def,2)),length(st)+i1-3,size(Def,2))';
  end

  if isempty(Head) disp('Head was not set');
  else  set(us(15),'userdata',Head); end
end

% help commands ---------------------------------------------------------------
elseif comstr(Cam,'help')  [CAM,Cam] = comstr(CAM,5);

help fecom

% info commands ---------------------------------------------------------------
elseif comstr(Cam,'info')
 
up=findobj(2,'type','axes','tag','Struct1');
if isempty(up) disp('No defined structure');
else

  us=get(up,'userdata');
  node=get(us(3),'userdata');nInf=size(node,1);
  elt=get(us(4),'userdata'); Line=get(us(6),'userdata');
  Patch=get(us(8),'userdata');  Sens=get(us(9),'userdata');
  Def=get(us(10),'userdata');   CDefE=[];CDefN=[];
  disp(sprintf('Structure with %i nodes and %i element group(s)', ...
       nInf-1,length(find(~isfinite(elt)))));
  disp(sprintf('%i deformation(s) and %i sensors',size(Def,2),size(Sens,1)));
  ua=get(IIAxes(l1,4),'userdata');
  disp(sprintf('The current axis is %i which contains',l1));
  st1=str2mat('patch','line','sensor','undeformed','nodetext','doftext', ...
      'arrow');
  for j1=2:size(ua,1)
   st=sprintf('Object %i %s def %i group(s)',j1-1,st1(ua(j1,1),:),ua(j1,4));
   st=[st sprintf(' %i',ua(j1,10:ua(j1,7)))];disp(st)
  end
  if nargout==1 out=ua(2:size(ua,1),3); end
end

% init commands ---------------------------------------------------------------
elseif comstr(Cam,'init')

  if nargin==1 feplot(CAM);
  elseif nargin==2 feplot(CAM,CAM1);
  elseif nargin==3 feplot(CAM,CAM1,CAM2);
  elseif nargin==4 feplot(CAM,CAM1,CAM2,CAM3);
  end

% line commands ---------------------------------------------------------------
elseif comstr(Cam,'line')  [CAM,Cam] = comstr(CAM,5);

if comstr(Cam,'init') IIgui(8,1)=2; feplot('PlotCreate 4 2'); return
elseif comstr(Cam,'plot') yplot=2;
elseif strcmp(Cam,'optim')

 if IIgui(6,8)==0 error('feplot InitNode: no current structure');
 else us=get(IIgui(6,8),'userdata'); end
 Line=get(us(6),'userdata');
 node=get(us(3),'userdata');nInf=size(node,1);
 elt=[]; for j1=1:size(Line,1);
  i1=[Line(j1,2:Line(j1,1)-1);Line(j1,3:Line(j1,1))];
  i1=i1(:,find(~any(i1==nInf)&~any(i1==0)))';
  elt=[elt;Inf abs('beam1'); i1 zeros(size(i1,1),4)];
 end
 Line=feutil('get Line',[],elt,nInf);
 set(us(6),'userdata',Line)

else out='unknown';end

% ----------------------------------------------------------------------------
elseif comstr(Cam,'nodetext') [CAM,Cam]=comstr(CAM,9);

 % find the current axis
 if IIAxes(l1,1)~=3 error('fecom nodetext: not a feplot axis');end
 eval('axes(IIAxes(l1,4));','');

 % find the undeformed line or make it
 if IIgui(6,8)==0 return; else us=get(IIgui(6,8),'userdata'); end
 ua=get(gca,'userdata');
 i1=find(ua(2:size(ua,1),1)==5)+1;
 if isempty(i1)
      ua(size(ua,1)+1,1:3)=[5 1 line('visible','off')];i1=size(ua,1);
 else
      uo=get(ua(i1(1),3),'userdata');
      eval('delete(uo(2:size(uo,1),1));','');
 end
 eval('findobj(ua(i1(1),3));','ua(i1(1),3)=line;')

 node=get(us(3),'userdata');
 NNode(node(:,1))=[1:length(node(:,1))]'; nInf=size(node,1);
 Xform=fecom('viewmtx',IIAxes(l1,:));

 % fat to mark the dots
 if comstr(Cam,'fat') i2='fat'; [CAM,Cam]=comstr(CAM,4); else i2=[]; end
 if comstr(Cam,'off') opt=[];
  uo=get(ua(i1(1),3),'userdata');uo=uo(2:size(uo,1),1:size(uo,2));
  eval('delete(uo(:,1))','');
  ua=ua(find([1:size(ua,1)]~=i1(1)),:);
 else
  opt=comstr(CAM,-1);
  if isempty(opt) & ~isempty(CAM1) & ~isstr(CAM1)
    opt=NNode(CAM1(:));
  elseif isempty(opt) opt=get(us(6),'userdata');opt=opt(:,2:size(opt,2));
      opt=opt(:);opt=opt(find(opt>0&opt<nInf)); % nodes of the current model
  else opt=opt(find(opt<=nInf));opt=NNode(opt); opt=opt(find(opt));end
  tnode=[node(:,1:4) [node(:,5:7) ones(size(node,1),1)]*Xform'];
  opt=find(sparse(opt(:)'+1,1,opt(:)'))'-1;

  if comstr(i2,'fat')
   if comstr(version,'5')
     line(tnode(opt(:),5),tnode(opt(:),6),tnode(opt(:),7),'marker','.',...
          'markersize',12,'color','r');
   else
     line(tnode(opt(:),5),tnode(opt(:),6),tnode(opt(:),7),'linestyle','.',...
          'markersize',12,'color','r');
   end 
     tnode(opt,5)=tnode(opt,5)+(max(tnode(opt,5))-min(tnode(opt,5)))/20;
  end

  uo=[];for j1=opt
     uo(j1) = text(tnode(j1,5),tnode(j1,6),tnode(j1,7), ...
          sprintf('%i',node(j1,1)),'FontSize',8, ...
          'userdata',[Inf 5],'visible','on','color',[1 1 1]-IIgui(1,3));
  end
  uo=uo';i3=find(uo); uo=[Inf 5 zeros(1,6);uo(i3) node(i3,:)];
  set(ua(i1(1),3),'userdata',uo);
 end
 
 set(gca,'userdata',ua);

% ----------------------------------------------------------------------------
elseif comstr(Cam,'optcheck') feplot('optcheck');
% ----------------------------------------------------------------------------
elseif strcmp(Cam,'orient') disp('This is done automatically in version 2.1');

% ----------------------------------------------------------------------------
% PatchInit PatchUpdate PatchFast
elseif comstr(Cam,'patch')  [CAM,Cam] = comstr(CAM,6);

  if     comstr(Cam,'init')   IIgui(8,1)=1;  feplot('PlotCreate 4 1');
  elseif comstr(Cam,'plot') IIgui(8,1)=1;  yplot=2;
  else

    i1 = comstr(CAM,[-1 IIAxes(l1,24)]); i1 = i1(1);
    if isempty(find([0:5]==rem(i1,10)))| ...
      isempty(find([0 1 2 3 4]==fix(i1/10)))  % face and edge colors
      disp(['fecom patchcolor value not valid' 7]);
    else IIAxes(l1,24)=i1;yplot=2;
    end
  end

% Modeshape scaling  ----------------------------------------------------------
elseif comstr(Cam,'scale') [CAM,Cam]=comstr(CAM,6);

 if IIgui(6,8)==0 return; else us=get(IIgui(6,8),'userdata'); end
 node=get(us(3),'userdata');
 NNode(node(:,1))=[1:length(node(:,1))]';nInf=size(node,1);

 Line=get(us(6),'userdata'); Def=get(us(10),'userdata');

 if comstr(Cam,'node') 
   i1=comstr(comstr(CAM,5),[-1 0]);Cam=[];
   if any(node(:,1)==i1(1))  IIAxes(l1,22)=i1(1);end
   disp('You should use fecom(''scaledof'')');
 elseif comstr(Cam,'dof')

   i1=comstr(comstr(CAM,4),[-1 0]);
   i1 = [fix(i1(1)) round(rem(i1(1),1)*100)];if i1(2)==0 i1(2)=1;end
   if i1(1)~=0 & (~any(node(:,1)==i1(1)) | ~any(i1(2)==[1 2 3]))
     out='not a valid dof';
   else IIAxes(l1,22)=rem(IIAxes(l1,22),1)+i1(1)+(i1(2)-1)*nInf;end

 elseif comstr(Cam,'equal')|comstr(Cam,'qual')
   i1=find(IIAxes(:,1)==3);IIAxes(i1,22)=fix(IIAxes(i1,22))+.01;
 elseif comstr(Cam,'auto')
   i1=find(IIAxes(:,1)==3);IIAxes(i1,22)=fix(IIAxes(i1,22));
 % ScaleDef changes current scaling length - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'def')

  i1=comstr(comstr(CAM,4),-1);
  if isempty(i1)   i1=IIgui(8,5); i2=0;
  elseif i1(1)==0  i1=1; i2=0;
  else i2 = phaseb(i1)/180*pi;
       IIgui(8,5)=abs(i1);i1 = i1(1)/IIAxes(l1,23);
  end
  IIAxes(l1,23) = IIAxes(l1,23)*i1;
  j1 = 26+[1:IIAxes(l1,30)]*5;          IIAxes(l1,j1+2) = IIAxes(l1,j1+3);
  IIAxes(l1,j1+1) = IIAxes(l1,j1+1)*i1; IIAxes(l1,j1+3) = IIAxes(l1,j1+3)+i2;

 % Other types of scaling - - - - - - - - - - - - - - - - - - - - - - - - -
 else

 % check of deformations numbers
 i1 = IIAxes(l1,26+[1:IIAxes(l1,30)]*5);
 i2 = find(i1>size(Def,2));i1(i2)=ones(size(i2));
 IIAxes(l1,26+[1:IIAxes(l1,30)]*5)=i1;
 if IIAxes(l1,23)==0 IIAxes(l1,23)=IIgui(8,5);end

 i2=fix(IIAxes(l1,22)); i3=fix(rem(IIAxes(l1,22),1)*100);
 if i2>0&rem(i2,nInf)~=0&i2<3*nInf
   i2=Def(i2,i1);i4=fix(IIAxes(l1,22))*ones(size(i2));%specified scaledof
   if i2==0 i2=Def(nInf*3+1,i1);i4=Def(nInf*3+2,i1);
    disp('ScaleDOF with zero response, resetting');
   end
 else  i2=Def(nInf*3+1,i1);i4=Def(nInf*3+2,i1); end % max response

 if i3==0     % all dofs to scaledef
  scales = 1; if length(i1)>1
    ind = [1:nInf-1 nInf+1:2*nInf-1 2*nInf+1:3*nInf-1];
    for j1=2:length(i1) scales(j1) = Def(ind,i1(j1))\ Def(ind,i1(1));
      scales(j1)=scales(j1)/abs(scales(j1));
    end
  end
  scales=scales.*IIAxes(l1,23)./abs(i2);
 elseif i3==1 % first def to ScaleDef and all other equal
  scales=IIAxes(l1,23)./i2(1)*ones(size(i2));
  i2=Def(i4(1),i1);i4=i4(1)*ones(size(i4));
 elseif i3==2 % try to match phases
  error('Not availlable yet');
 end

 j1 = 26+[1:IIAxes(l1,30)]*5;
 IIAxes(l1,j1+1) = abs(scales);
 IIAxes(l1,j1+2) = phaseb(scales)/180*pi;
 IIAxes(l1,j1+3) = phaseb(scales)/180*pi;

 % display of scaling information - - - - - - - - - - - - - - - - - - - -
 if ~isempty(strfind(Cam,'set'))
 elseif IIAxes(l1,23)==eps
 else
    i4=node(rem(i4,nInf),1)'+(fix(i4/nInf)+1)/100;
    disp(sprintf(['D %i. Amp %10.2e , %6.1f (deg) ' ...
         'DOF %.2f : %11.2e, %6.1f(deg)\n'], ...
         [i1;abs(scales);phaseb(scales);i4;abs(i2);phaseb(i2)]))
 end
 Cam=[];
 end % of scale subcommand - - - - - - - - - - - - - - - - - - - -
 if ~isempty(Cam) yplot=2;end

% ----------------------------------------------------------------------------
elseif comstr(Cam,'sens')  [CAM,Cam] = comstr(CAM,5);

  if     comstr(Cam,'init')   IIgui(8,1)=3;  feplot('PlotCreate 4 3');
  elseif comstr(Cam,'plot') IIgui(8,1)=3;  yplot=2;
  end

% ----------------------------------------------------------------------------
elseif comstr(Cam,'set')[CAM,Cam] = comstr(CAM,4);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if comstr(Cam,'slider');[CAM,Cam] = comstr(CAM,7);% slidermenu handling

if isempty(Cam) j1=[1:3]*3-2; else j1=[sscanf(Cam(1),'%i')*3-2 1];j1=j1(1); end
Cam =Cam(2:length(Cam));
if isempty(l1) return; end

for j3 = j1 % loop on sliders to update
  eval('i1=get(IIgui(7,j3),''value'');','i1=[];');i2=[];
  if ~isempty(i1)
   if     strcmp(Cam,'text')   i2 = str2num(get(IIgui(7,j3+1),'string'));
	     if i1==7 | i1==9  i2 = log10(i2); end
   elseif strcmp(Cam,'value')  i2 = get(IIgui(7,j3+2),'value'); end

   i3=[0 0 0]; % for orient, center, update
   if ~isempty(i2)
     if     i1==1 & IIAxes(l1,11)~=i2  IIAxes(l1,11) = i2;i3(1:2)=[1 1];
     elseif i1==2 & IIAxes(l1,12)~=i2  IIAxes(l1,12) = i2;i3(1:2)=[1 1];
     elseif i1==3 & IIAxes(l1,13)~=i2  IIAxes(l1,13) = i2;i3(1:2)=[1 1];
     elseif i1==4 & IIAxes(l1,15)~=i2  IIAxes(l1,15) = i2;i3(2)=1;
     elseif i1==5 & IIAxes(l1,16)~=i2  IIAxes(l1,16) = i2;i3(2)=1;
     elseif i1==6 & IIAxes(l1,17)~=i2  IIAxes(l1,17) = i2;i3(2)=1;
     elseif i1==7 & IIAxes(l1,18)~=i2  IIAxes(l1,18) = 10^i2;i3(2)=1;
     elseif i1==8 & IIAxes(l1,21)~=i2
       if i2>IIAxes(l1,21) IIAxes(l1,21) = max(fix(i2),IIAxes(l1,21)+1);
       else IIAxes(l1,21) = min(fix(i2),IIAxes(l1,21)-1); end
     elseif i1==9 & IIAxes(l1,23)~=i2  IIAxes(l1,23) = 10^(i2);
     end
   end

   if IIAxes(l1,23)==0 IIAxes(l1,23)=eps; end
   if IIAxes(l1,18)==0 IIAxes(l1,18)=1; end

   if     i1==1 i2 = [IIAxes(l1,11) -90 270 IIAxes(l1,11)];	
   elseif i1==2 i2 = [IIAxes(l1,12) -90 270 IIAxes(l1,12)];	
   elseif i1==3 i2 = [IIAxes(l1,13) -90 270 IIAxes(l1,13)];	
   elseif i1==4 i2=1:size(node,1)-1; i2=IIAxes(l1,15)+ ...
            [0 max(max(abs(IIAxes(l1,15:17))),1)*[-1 1] 0];
   elseif i1==5 i2=1:size(node,1)-1; i2=IIAxes(l1,16)+ ...
            [0 max(max(abs(IIAxes(l1,15:17))),1)*[-1 1] 0];
   elseif i1==6 i2=1:size(node,1)-1; i2=IIAxes(l1,17)+ ...
            [0 max(max(abs(IIAxes(l1,15:17))),1)*[-1 1] 0];
   elseif i1==7 i2 = [IIAxes(l1,18) -.301 log10(IIAxes(l1,18))+[2 0]];
   elseif i1==8 i2 = [IIAxes(l1,21) 2 100 IIAxes(l1,21)];	
   elseif i1==9 i2 = [IIAxes(l1,23) log10(IIAxes(l1,23))+[-1 1 0]];
   end

   set(IIgui(7,j3+1),'string',num2str(i2(1)));
   set(IIgui(7,j3+2),'min',i2(2),'max',i2(3),'value',i2(4));
  end
end

% orient and center

if IIgui(6,6)~=0  set(IIMenu(40),'checked','on');
     eval('set(IIgui(6,6),''xform'',fecom(''viewmtx'',IIAxes(l1,:)));', ...
          'IIgui(6,6)=0;');
else     set(IIMenu(40),'checked','off'); end

% reset values for defaults
if IIAxes(l1,1)==3
  IIgui(8,1:7)=IIAxes(l1,[2 19 21 22 23 24 26]);
  IIgui(9,1:3)=IIAxes(l1,[11:13]);
end
if isempty(Cam) comgui('refresh'); else yplot=2; end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'o')

   ua=get(IIAxes(l1,4),'userdata');
   i1=strfind(Cam,'def');if ~isempty(i1) ind(1)=i1;end
   i1=strfind(Cam,'group');if ~isempty(i1) ind(2)=i1;end
   i1=strfind(Cam,'ty');if ~isempty(i1) ind(3)=i1;end
   ind(4)=length(Cam)+1;

   i2=fix(comstr(CAM(2:min(ind(find(ind~=0)))-1),[-1]))';
   i2=i2(find(i2>0));

   if isempty(i2) error('fecom(''set o'') : no corresponding object');end
   if ind(3)~=0 % type
     i3=comstr(CAM(ind(3):min(ind(find(ind>ind(3))))-1),'type','%i');
     if any([1:4 7]==i3(1)) ua(i2+1,1)=i3(1)*ones(size(i2));
     elseif i3(1)==0
       if any(ua(i2+1,1)==[2 3 4 7]) eval('delete(ua(i2+1,3))','');
       elseif any(ua(i2+1,1)==[1])
         eval('delete(get(ua(i2+1,3),''userdata''))','');
       elseif any(ua(i2+1,1)==[5 6])
        eval('uo=get(ua(i2+1,3),''userdata'');delete(uo(2:size(uo,1),1));','');
       end
       ua(i2+1,1)=zeros(size(i2));
     else error('not a valid type to set'); end
   end
   if ind(1)~=0 % deformation
     i3=comstr(CAM(ind(1):min(ind(find(ind>ind(1))))-1),'deformation','%i');
     if isempty(i3)|i3>IIAxes(l1,30)
        error('not a valid deformation number');
     end
     ua(i2+1,4)=i3(1)*ones(size(i2));
   end
   if ind(2)~=0 % group
     i3=fix(comstr(CAM(ind(2)+5:min(ind(find(ind>ind(2))))-1),[-1]));
     ua(i2+1,7)=9+length(i3); ua(i2+1,10:9+length(i3))=i3(:)';
   end
   set(IIAxes(l1,4),'userdata',ua);
   yplot=2;

else out='unknown';end

% ----------------------------------------------------------------------------
elseif comstr(Cam,'show') iicom(CAM,CAM1);

% ----------------------------------------------------------------------------
elseif comstr(Cam,'titopt') iicom(CAM,CAM1);

% ----------------------------------------------------------------------------
elseif comstr(Cam,'triaxmenu')

   h = get(get(IIgui(2,1),'parent'),'currentmenu');
   if strcmp(get(h,'checked'),'on')
      set(h,'checked','off');
      eval('delete(IIgui(6,6));IIgui(6,6)=0;','IIgui(6,6)=0;')
   else fecom('triax',';');set(h,'checked','on'); end

elseif comstr(Cam,'triax')

   if (comstr(get(IIMenu(40),'checked'),'on')&~comstr(Cam,'triaxon')) | ...
      comstr(Cam,'triaxoff')
     eval('delete(IIgui(6,6));IIgui(6,6)=0;','IIgui(6,6)=0;')
     set(IIMenu(40),'checked','off');
     if nargout>0|~comstr(out,'done') out0=out; end
     return;
   end

   IIAxes(l1,3)=2;figure(IIAxes(l1,3)); eval('ax = axes(IIAxes(l1,4))','ax=gca;');
   chi = get(IIAxes(l1,3),'children');

   for j1 = 1:length(chi) ua = get(chi(j1),'userdata');
     if size(ua,2)>1 if all(ua(1,1:2)==[Inf 4]) IIgui(6,6)=chi(j1); end;end
   end
   eval('get(IIgui(6,6),''type'');','IIgui(6,6)=0;')
   Xform=fecom('viewmtx',IIAxes(l1,:));
   if IIgui(6,6)==0 | isempty(get(IIgui(6,6),'children'))
	IIgui(6,6) = axes('units','inches','position',[.4 .7 .25 .25], ...
           'visible','off','xform',Xform,'userdata',[Inf 4]);
	line([1 0 0 0 0],[0 0 1 0 0],[0 0 0 0 1],'color',[1 1 1]-IIgui(1,3),...
             'clipping','off')
	text([1.3 0 0],[0 1.3 0],[0 0 1.3],['x';'y';'z'],'fontsize',9, ...
             'HorizontalAlignment','center','VerticalAlignment','middle');

   else
        set(IIgui(6,6),'xform',Xform)
        set(get(IIgui(6,6),'children'),'visible','on')
   end
   set(IIMenu(40),'checked','on');axes(ax);

% ----------------------------------------------------------------------------
elseif comstr(Cam,'undef') [CAM,Cam]=comstr(CAM,6);

  ua=get(IIAxes(l1,4),'userdata');
  i1=[];if size(ua,1)>2 i1=find(ua(2:size(ua,1),1)==4); end
  if ~isempty(i1) end
  if isempty(i1)
    disp('fecom undef: no undeformed structure object in current plot');
  else
    i1=i1(1);eval('get(ua(i1+1,3),''type'');','i1=0;');
    if i1==0 yplot=2;i1=line;else i1=ua(i1+1,3); end

    if isempty(Cam) & strcmp(get(i1,'visible'),'on')
     set(i1,'visible','off'); if IIAxes(l1,19)==0 IIAxes(l1,19)=2;end
    elseif isempty(Cam)
      set(i1,'visible','on');if IIAxes(l1,19)==0 IIAxes(l1,19)=2; end
    elseif strcmp(Cam,'dot')
      if comstr(version,'5')
       set(i1,'LineStyle','none','marker','.','Markersize',3,'visible','on');
       IIAxes(l1,19)=1;
      else
       set(i1,'LineStyle','.','Markersize',3,'visible','on');IIAxes(l1,19)=1;
      end
    elseif strcmp(Cam,'line')
      if comstr(version,'5')
       set(i1,'LineStyle',':','marker','none','Markersize',3,'visible','on');
       IIAxes(l1,19)=2;
      else
       set(i1,'LineStyle',':','Markersize',6,'visible','on');IIAxes(l1,19)=2;
      end
    end
  end

% ----------------------------------------------------------------------------
elseif comstr(Cam,'viewmtx') [CAM,Cam]=comstr(CAM,5);

  if length(CAM1)>13 i1=CAM1(1,11:13)*pi/180;i1=[cos(i1);sin(i1)]';
  else i1=IIAxes(l1,11:13)*pi/180;i1=[cos(i1);sin(i1)]'; end
  out=[i1(3,1) i1(3,2) 0 0;-i1(3,2) i1(3,1) 0 0;0 0 1 0;0 0 0 1] * ...
      [i1(1,1) i1(1,2) 0 0;-i1(2,2)*i1(1,2) i1(2,2)*i1(1,1) i1(2,1) 0;
       i1(2,1)*i1(1,2) -i1(2,1)*i1(1,1) i1(2,2) 0; 0 0 0 1];

%  out=[cos(i1)		sin(i1)		0	0
%      -sin(i2)*sin(i1)  sin(i2)*cos(i1) cos(i2) 0
%       cos(i2)*sin(i1) -cos(i2)*cos(i1) sin(i2) 0
%       0                0               0       1];
%  if i3~=0 out= ...
%      [cos(i3) sin(i3) 0 0;-sin(i3) cos(i3) 0 0;0 0 1 0;0 0 0 1]*out;
%  end

% ----------------------------------------------------------------------------
elseif comstr(Cam,'view') [CAM,Cam]=comstr(CAM,5);

  eval('axes(IIAxes(l1,4));','figure(IIAxes(l1,3));IIAxes(l1,4)=gca;');
  if CAM(1)=='('; CAM=CAM(2:length(CAM)-1); end
  opt=comstr(CAM,[-1]);
  if     length(opt)==1 & opt(1)==2 IIAxes(l1,11:13)=[0 90 0];
  elseif length(opt)==1 & opt(1)==3 IIAxes(l1,11:13)=[-37.5 30 0];
  else if length(opt)<4 opt(4)=0;end;IIAxes(l1,11:14)=rem(opt(1:4),360);
  end
  yplot=2;

% ----------------------------------------------------------------------------
elseif comstr(Cam,'quit')
% ----------------------------------------------------------------------------
else

  eval([CAM CAM1],'out = ''unknown''; ')

end % of selection of commands

if ~isstr(CAM1)&~isempty(CAM1) CAM1 = ';'; end
if ~isempty(CAM1)&CAM1==';'
elseif yplot==2         feplot('plot');
elseif yplot==3         feplot('plotall');
elseif yplot==4         feplot('plotCreate');
elseif yplot==5         feplot('PlotAllCreate');
elseif yplot==1         drawnow;	end
if strcmp(out,'unknown') comgui(['@ Unknown: fecom(''' CAM0 ''')']);end

if nargout>0|~comstr(out,'done') out0=out; end


