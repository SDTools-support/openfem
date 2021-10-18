function [out,out1] = feplot(CAM,CAM1,CAM2,CAM3,CAM4,CAM5)

%FEPLOT	initialisation of 3-D deformation plots
%
%	Synopsis : feplot(node,elt,  mode,mdof,opt)
%		   feplot(node,elt,  mode,mdof,opt,cmode)
%		   feplot(node,LDraw,mode,mdof,opt)
% 
%	FECOM provides a UI command mode for modifications of the plot.
%	Commands used for the initialisation of a plot are also accessible 
%	with the following formats
%
%	  feplot('Initmodel',node,elt)
%	  feplot('InitDef',mode,mdof)
%         feplot('InitCdef',cmode,mdof)
%         feplot('InitCdef',STRESS)
%         feplot('InitCdef',ENERGY)
%
%	FEPLOT arguments are
%
%	NODE coordinates in the Universal File format
%	     ( ,1)      node number (positive integers)
%	     ( ,2:4)    (unused) PosCoorID DisCoorID NodeGroupID
%	     ( ,5:7)    x-y-z global coordinates of the node 
%
%	ELT    finite element model description matrix (see fe_mk for details)
%
%	MODE,MDOF deformation shape and corresponding DOF definition vector
%	  The first column of MDOF gives DOF identifications of the form 
%	  Node#.Dof# with the following DOFs used by FE_PLOT
%	     node#.01 to node#.03 are the x, y, z translations
%	     node#.07 to node#.09 are the -x, -y, -z translations
%	  Unspecified DOFs are set to zero. XFdof (see XFopt) can be used to
%	  animate the deformation corresponding to a frequency point of a set
%	  of measured FRFs.
%	  To construct interpolations of DOFs not included in your mode
%	  consider using FE_EXP.
%
%	FEopt options for the plot
%	 (1,1) =  plot type  1: patch, 2: lines, 3: sensor sticks
%	 (1,2) =  00: none	 01: UndefDot		02: UndefLine
%        (1,3) =  number of deformations per cycle
%        (1,4) =  number of the node used for modeshape scaling
%        (1,5) =  maximum displacement
%		  (A negative number forces unique scaling for multiple modes)
%
%	CMODE is used to specify values on the colorscale for patch plots
%	   NODE COLORS can be defined directly using a size(node,1) by
%	   size(mode,2) matrix or a size(mdof,1) by size(mode,2) (the color is
%	   then proportional to the sum of the squares of DOF translations
%	   associated to the node)
%	   ELEMENT COLORS can be defined using a size(elt,1) by size(mode,2)
%	   matrix where a row of CMODE gives the color value of the element
%	   described in the same row of ELT.
%
%       THIS IS NOT THE CURRENT FEPLOT FROM SDTOOLS
%
%	See also FECOM, FE_MK, FE_EXP

%	Etienne Balmes
%       Copyright (c) 2001-2003 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       $Revision: 1.8 $  $Date: 2015/01/12 18:44:32 $


global IIgui IIAxes

if nargin==0
 if ~ishandle(2)|isempty(findobj(2,'tag','Struct1')) ;comgui('guifeone');
 else feplot('plotall'); end
 return;
elseif isstr(CAM) [CAM,Cam]=comstr(CAM,1);
elseif isstruct(CAM) 
  feplot('initmodel',CAM);
  if nargin==2 & isstruct(CAM1) feplot('initdef',CAM1);end
  return;
else % initialisation
   if nargin<4 error('4 arguments needed for initialisation'); end
   if nargin==4 CAM4=2;end
   feplot('optcheck',CAM4);
   feplot('initmodelback',CAM,CAM1);
   feplot('initdefBack',CAM2,CAM3);
   if nargin==6 feplot('initcdefback',CAM5,CAM3); end
   fecom('GroupAll');
   return
end

if isempty(IIAxes) l1=[]; else l1 = find(IIAxes(:,1)==3); end
if ~comstr(Cam,'optcheck')& ~comstr(Cam,'supported')
  if isempty(l1) iicom('sub1 1 1 3 2'); l1 = find(IIAxes(:,1)==3); end
  if ~isempty(find(l1==IIgui(3,7))) l1=IIgui(3,7);
  else   l1=l1(1);IIgui(3,7)=l1; end
  if IIgui(3,5)>100 disp(['feplot(''' Cam ''')']); end
end

eval('set(IIgui(5,2),''backgroundcolor'',[1 0 0]);','IIgui(5,2)=0;')
if nargin==1 CAM1=[]; end

% IIgui(6,5) uo contains pointers to structure axes
% For each structure axis object one defines the objects
%  3  node
%  4  elt
%  5  group
%  6  Line
%  7  Line2
%  8  Patch
%  9  Sens
%  10 Def
%  11 mdof
%  12 coor
%  13 CDefE
%  14 CDefN
%  15 Labels
% For each plot axis, one defines the objects
if comstr(Cam,'update')
  [CAM,Cam]=comstr(CAM,7);CAM=['plot' CAM];Cam=['plot' Cam];
  disp('use the plot command instead of update');
end
% ----------------------------------------------------------------------------
if comstr(Cam,'plot') [CAM,Cam]=comstr(CAM,5);

% Make sure there is a current structure

eval('findobj(IIgui(6,8));','IIgui(6,8)=0;')
if IIgui(6,8)==0
  up=findobj(2,'type','axes','tag','Struct1');
  set(IIgui(6,5),'userdata',up);
  if ~isempty(up) IIgui(6,8)=up(1);
  else disp(['feplot : No defined structure' 7]); end
end
if IIgui(6,8)==0  comgui(sprintf('CheckAxes 3 %i %i',IIAxes(l1,3),l1))
     cla;text(.5,.5,'No defined structure','units','normalized', ...
              'HorizontalAlignment','center');
     comgui('@No defined structure');return;
else us=get(IIgui(6,8),'userdata'); end

% PlotAll - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if comstr(Cam,'all') [CAM,Cam]=comstr(CAM,4);

     if comstr(Cam,'a') [CAM,Cam]=comstr(CAM,2);
       i1=find(IIAxes(:,1)==3); i2=IIgui(3,7);
     else
       i1=find(IIAxes(:,1)==3&IIAxes(:,3)==IIgui(1,2)); i2=IIgui(3,7);
     end
     for j1=i1' IIgui(3,7)=j1; feplot(['plot' CAM],';');end
     IIgui(3,7)=i2;comgui('refresh');

% PlotUsual - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
else

  st1='center';

  if comstr(Cam,'create') ua=[];
  elseif comstr(Cam,'add')      error('Not available');ua=[];
  elseif nargin==2 & ~isstr(CAM1) ua=CAM1;
  elseif ~ishandle(IIAxes(l1,4)) ua=[];
  else  ua=get(IIAxes(l1,4),'userdata');end
  if isempty(ua)|comstr(Cam,'create')
     [CAM,Cam]=comstr(CAM,7);
     i1 = comstr(CAM,[-1 IIAxes(l1,2) 4]);
     ua=[Inf 3 i1(1) 0];
     % check that a new axis exists and clear it
     comgui(sprintf('CheckAxes 3 %i %i',IIAxes(l1,3),l1));
     % duplicate objects types that need it (123 and 7)
     i3=[1:IIAxes(l1,30)]';
     if i1(2)==4 i1=[4 i1(1)*ones(1,IIAxes(l1,30));0 1:IIAxes(l1,30)]';
     else i1 = [i1 ; ones(size(i1))]'; end
     ua(2:size(i1,1)+1,[1 4])= i1;
     IIgui(8,1)=ua(3,1);IIAxes(l1,2)=ua(3,1);
     axes(IIAxes(l1,4));cla
     IIAxes(l1,29)=line('visible','off');    iimouse; st1='centernode';
  end
  if size(ua,2)<10 ua(1,10)=0;end

  if isempty(strfind(Cam,'noscale')) fecom('scale'); else st1='';end
  axes(IIAxes(l1,4));

% Local variables are
%   lSel : selected line
%   nSel : nodes of the selected line
%   jDef : index of deformation in IIAxes
%   j2   : loop on objects in ua 2:size(ua,1)
%   uo   : object user data
%   ua   : axis userdata
%   Nr   : projected reference position of nodes
%   Nc   : projected current position of nodes
%    [Inf 3 RefDefType ... 
%    [Type ParentStruct Ptr DefIndex ColorIndex LineStyleIndex nGroup]
%    [10: (groupnumbers)]
%    Recogized types are
%    1: deformed patch  2: deformed line 3: deformed sensor line
%    4: undeformed line 5: node text     6: DOF text
%    7: sensor with arrow
%    negative types are made invisible

  % get associated structure info
  node=get(us(3),'userdata');nInf=size(node,1);
  elt=get(us(4),'userdata'); Line=get(us(6),'userdata');
  Patch=get(us(8),'userdata');  Sens=get(us(9),'userdata');
  Def=get(us(10),'userdata');   CDefE=[];CDefN=[];
  Head = get(us(15),'userdata');
  cmap=get(gca,'colororder');   emode=fecom('EraseMode');

% checks on inits and definition of variables - - - - - - - - - - - - - - - - -

  j2=2:size(ua,1);
  if any(ua(j2,1)==2 | ua(j2,1)==4)&isempty(Line)
      feplot('InitLine');Line=get(us(6),'userdata');
  end
  if any(ua(j2,1)==1)&isempty(Patch)
      feplot('InitPatch');Patch=get(us(8),'userdata');
  end
  if any(ua(j2,1)==3|ua(j2,1)==7)&isempty(Sens)
      feplot('InitSensback');Sens=get(us(9),'userdata');
  end
  if any(ua(j2,1)==1|ua(j2,1)==2|ua(j2,1)==3)
     if isempty(Def)
       error('Use feplot(''initdef'') to declare a deformation');
     end
   if any(IIAxes(l1,26+[1:IIAxes(l1,30)]*5)>size(Def,2))
       disp('Reset to deformation 1');
       i3=IIAxes(l1,30+IIAxes(l1,30)*5+[1:IIAxes(l1,25)]); % group selected
       IIAxes(l1,30:31)=[1 1];
       IIAxes(l1,30+IIAxes(l1,30)*5+[1:IIAxes(l1,25)])=i3;
   end
  end
  if IIAxes(l1,2)==1 & IIAxes(l1,26)==2   % node color
      if isempty(CDefN) CDefN=get(us(14),'userdata'); end
      if isempty(CDefN) feplot('initCDef','a');CDefN=get(us(14),'userdata');end
  end
  if IIAxes(l1,2)==1 & IIAxes(l1,26)==1   % element color
    if isempty(CDefE) CDefE=get(us(13),'userdata');
     if isempty(CDefE) feplot('initCDef','group');CDefE=get(us(13),'userdata');
     end
    end
    if ~all(size(CDefE,1)==size(elt,1))
      disp('warning incoherent element color data');
    end
  end

  % check the number of groups
   i2=[];i1=IIAxes(l1,30+IIAxes(l1,30)*5+[1:IIAxes(l1,25)]);i1=i1(find(i1));
   if ~isempty(Line) i2=i1(find(i1<=size(Line,1)));
   else
   end
   if isempty(i2) fecom('groupallback');
   elseif length(i2)~=length(i1)
     disp(['Selected element group(s) ' sprintf('%i ',i2)]);
     IIAxes(25)=length(i2);
     IIAxes(l1,30+IIAxes(l1,30)*5+[1:IIAxes(l1,25)])=i2(:)';
   end
  i1=IIAxes(l1,11:13)*pi/180;i1=[cos(i1);sin(i1)]';
  Xform=[i1(3,1) i1(3,2) 0 0;-i1(3,2) i1(3,1) 0 0;0 0 1 0;0 0 0 1] * ...
      [i1(1,1) i1(1,2) 0 0;-i1(2,2)*i1(1,2) i1(2,2)*i1(1,1) i1(2,1) 0;
       i1(2,1)*i1(1,2) -i1(2,1)*i1(1,1) i1(2,2) 0; 0 0 0 1];
  %Xform=fecom('viewmtx',IIAxes(l1,:));
  Nr=[node(:,5:7) ones(size(node,1),1)]*Xform';

ua=ua(find(ua(:,1)),:);

% loop on objects to be created/updated -------------------------------------
for j2=2:size(ua,1)  

 eval('findobj(ua(j2,3));','ua(j2,3)=0;');
 % check deformation number
 jDef=ua(j2,4); if jDef~=0&jDef<1|jDef>IIAxes(l1,30)
  jDef=0;disp(sprintf('Object %i using invalid deformation (%i)',j2,jDef));
 elseif jDef>0 jDef=26+jDef*5; end
 % check group number
 if ua(j2,7)==0
   i1=IIAxes(l1,30+IIAxes(l1,30)*5+[1:IIAxes(l1,25)]);i1=i1(find(i1));
   indG=i1(:)';
 else
   indG = ua(j2,10:ua(j2,7));indG=indG(find(indG>0&indG<=size(Line,1)));
 end

 % definition of the selected line
 lSel=[];for j1=indG
        lSel=[lSel Line(j1,2:Line(j1,1)) nInf];
 end
 i1 = find(lSel==0); if ~isempty(i1) lSel(i1)=i1*0+nInf; end
 nSel=[]; if ~isempty(lSel) nSel=find(sparse(lSel+1,1,lSel))-1;
 else disp(sprintf('Warning: no line for object %i',j2-1));end

% create an undeformed line object - - - - - - - - - - - - - - - -
if ua(j2,1)==4

   if isempty(lSel) lSel=1:size(node,1)-1;IIAxes(l1,19)=01; end

   if ua(j2,3)==0
     ua(j2,1:3)=[4 1.01 ...
        double(line(Nr(lSel,1),Nr(lSel,2),Nr(lSel,3), ...
        'color',[1 1 1]-IIgui(1,3), ...
        'userdata',[Inf 4 0 0 0 0 lSel],'clipping','off','markersize',3))];
     IIAxes(l1,20)=ua(j2,3);
   else
    set(ua(j2,3),'xdata',Nr(lSel,1),'ydata',Nr(lSel,2),'zdata',Nr(lSel,3), ...
        'userdata',[Inf 4 0  0 0 0 lSel]);
   end
   if comstr(version,'5')
      if any([2 3]==rem(IIAxes(l1,19),10))
        set(ua(j2,3),'linestyle',':','marker','none');
      else % if any([1 3]==rem(IIAxes(l1,19),10))
        set(ua(j2,3),'linestyle','none','marker','.');
      end
   else
     if any([1 3]==rem(IIAxes(l1,19),10)) set(ua(j2,3),'linestyle','.');
     else set(ua(j2,3),'linestyle',':');
     end
   end
   if rem(IIAxes(l1,19),10)==0  set(ua(j2,3),'visible','off'); end

% update a node text object
elseif ua(j2,1)==5

  uo=get(ua(j2,3),'userdata');uo=uo(2:size(uo,1),1:size(uo,2));
  tnode=[[uo(:,6:8) ones(size(uo,1),1)]*Xform'];
  for j1=1:size(uo,1) set(uo(j1,1),'position',tnode(j1,1:3)); end

% update a dof text object
elseif ua(j2,1)==6 

  uo=get(ua(j2,3),'userdata');uo=uo(2:size(uo,1),1:size(uo,2));
  tnode=[[uo(:,2:4) ones(size(uo,1),1)]*Xform'];
  for j1=1:size(uo,1) set(uo(j1,1),'position',tnode(j1,1:3)); end

% line object creation - - - - - - - - - - - - - - - - - - - - - - - -
elseif ua(j2,1)==2

 if ~comstr(get(ua(j2,3),'type'),'line') %create
   ua(j2,1:2)=[2 1.01];     uo=[Inf 2 0  0 0 0 lSel]; LD2=lSel;
   ua(j2,3)=double(line('color',cmap(remi(ua(j2,4),size(cmap,1)),:), ...
       'clipping','off'));
 end

 % current deformation at needed nodes with projection
 uo=[Inf 2 0  0 0 0 lSel]; LD2=lSel;
 Nc=Nr;
 if jDef>0 Nc(nSel,:)=Nc(nSel,:)+[real([Def(nSel,IIAxes(l1,jDef)) ...
          Def(nSel+nInf,IIAxes(l1,jDef)) Def(nSel+2*nInf,IIAxes(l1,jDef))]* ...
          IIAxes(l1,jDef+1)*exp(1i*IIAxes(l1,jDef+2))) ones(length(nSel),1)] ...
          *Xform';uo(3)=IIAxes(l1,jDef);
 else uo(3)=0; 
 end

 set(ua(j2,3),'xdata',Nc(LD2,1),'ydata',Nc(LD2,2),'zdata',Nc(LD2,3), ...
     'userdata',uo)

% patch object creation/update - - - - - - - - - - - - - - - - - - - - - -
elseif ua(j2,1)==1

 ua(j2,1:2)=[1 1.01];uo=[Inf 1 0  0 0 0];
 % group selection
  i3=indG;
  i3=full(sparse(i3,1,1:length(i3),max([Patch(2:size(Patch,1),6);i3(:)]),1));
  i3=find(i3(Patch(2:size(Patch,1),6)))+1;
  if ~isempty(i3) nSel=Patch(i3,1:4);nSel=find(sparse(nSel(:)+1,1,nSel(:)))-1;
  else nSel=[];end %nodes needed

if comstr(version,'5') & ~isempty(i3)

 Nc=Nr;if jDef>0 Nc(nSel,:)=Nc(nSel,:)+[real([Def(nSel,IIAxes(l1,jDef)) ...
          Def(nSel+nInf,IIAxes(l1,jDef)) Def(nSel+2*nInf,IIAxes(l1,jDef))]* ...
          IIAxes(l1,jDef+1)*exp(i*IIAxes(l1,jDef+2))) ones(length(nSel),1)] ...
          *Xform';
 end
 uo = [uo i3(:)'];
 i6=[];i6(nSel+1)=1:length(nSel); i6(1)=NaN;
 i6 = reshape(i6(Patch(i3,1:4)+1),length(i3),4); %face connectivity matrix

 if ~comstr(get(ua(j2,3),'type'),'patch') % creation
     ua(j2,3) = patch('vertices',Nc(nSel,1:3),'faces',i6,'facecolor','w', ...
      'clipping','off','userdata',uo,'BackFaceLighting','lit');
    % do something for mass and beams
 elseif ua(j2,3)~=0
     set(ua(j2,3), 'vertices',Nc(nSel,1:3),'faces',i6,'userdata',uo);
     % do something for mass and beams
 end

 i4 = [rem(IIAxes(l1,24),10) fix(IIAxes(l1,24)/10)];
 % DefColor: 0 Order 1 Element 2 Node 3 Group  4 Lighting 5 element nodes
 %   EdgeColor 00 none 01 interp 02 flat 03 white 04 black 05 order
 %   FaceColor 00 none 01 interp 02 flat 03 white 04 black 05 order

 if IIAxes(l1,26)==3       % group 

   if     i4(1)==3  st2 = 'w';             elseif i4(1)==4 st2 = 'k';
   elseif i4(1)==5  st2 = cmap(remi(uo(1,3),size(cmap,1)),:);
   else    st2 = 'none'; end
   if     i4(2)==3  st  = 'w';             elseif i4(2)==4 st2 = 'k';
   else st  = 'flat'; end

   set(ua(j2,3),'cdata',Patch(i3,6),'edgecolor',st2,'facecolor',st)

 elseif IIAxes(l1,26)==2   % node 

   if uo(1,3)>size(CDefN,2) error('Not enough CDefN vectors');end
   if     i4(1)==3  st2 = 'w';             elseif i4(1)==4 st2 = 'k';
   elseif i4(1)==4  st2 = cmap(remi(uo(1,3),size(cmap,1)),:);
   else    st2 = 'interp'; end
   if     i4(2)==3  st  = 'w';             elseif i4(2)==4 st2 = 'k';
   elseif i4(2)==4  st  = cmap(remi(uo(1,3),size(cmap,1)),:);
   else st  = 'interp'; end
   set(ua(j2,3),'cdata',CDefN(nSel,IIAxes(l1,jDef)), ...
           'edgecolor',st2,'facecolor',st)

 elseif IIAxes(l1,26)==1   % element

   if     i4(1)==3  st2 = 'w';             elseif i4(1)==4 st2 = 'k';
   elseif i4(1)==4  st2 = cmap(remi(uo(1,3),size(cmap,1)),:);
   else    st2 = 'none'; end

   if     i4(2)==3  st  = 'w';             elseif i4(2)==4 st2 = 'k';
   elseif i4(2)==4  st  = cmap(remi(uo(1,3),size(cmap,1)),:);
   else st  = 'flat'; end

   set(ua(j2,3),'CData',CDefE(Patch(i3,5),IIAxes(l1,jDef)), ...
           'edgecolor',st2,'facecolor',st)

 elseif IIAxes(l1,26)==5   % element node
   warning('ColordataELN :Not availlable with the vectorized patch');

 else % default 
   if IIAxes(l1,26)==4 % lighting
     warning('Using default coloring, use light object for lighting');
   end
   if     i4(1)==3  st2 = 'w';             elseif i4(1)==4 st2 = 'k';
   elseif i4(1)==4  st2 = cmap(remi(uo(1,3),size(cmap,1)),:);
   else    st2 = 'k'; end
   if     i4(2)==3  st  = 'w';             elseif i4(2)==4 st2 = 'k';
   elseif i4(2)==4  st  = cmap(remi(uo(1,3),size(cmap,1)),:);
   else st  = cmap(remi(uo(1,3),size(cmap,1)),:); end
   set(ua(j2,3),'cdata',[],'facecolor',st,'edgecolor',st2)
 end % of color selection

else % now Matlab version 4.0
 % current deformation at needed nodes with projection
 Nc=Nr;if jDef>0 Nc(nSel,:)=Nc(nSel,:)+[real([Def(nSel,IIAxes(l1,jDef)) ...
          Def(nSel+nInf,IIAxes(l1,jDef)) Def(nSel+2*nInf,IIAxes(l1,jDef))]* ...
          IIAxes(l1,jDef+1)*exp(i*IIAxes(l1,jDef+2))) ones(length(nSel),1)] ...
          *Xform';
 end

  % i2 will contain the pointers
  if ua(j2,3)~=0 i2=get(ua(j2,3),'userdata'); end
  % sort i3 by barycenter for better hidden face removal

 i5=[Nc(Patch(i3,1),3)+Nc(Patch(i3,2),3)+Nc(Patch(i3,3),3)];
 i6=find(Patch(i3,4)); if ~isempty(i6) i5(i6)=i5(i6)+Nc(Patch(i3(i6),4),3);end
 [i5,i6]=sort(i5);i3=i3(i6);
 if comstr(version,'5') st3='marker'; else st3='linestyle';end
 for jElt=1:length(i3)

    LD2 = nonzeros(Patch(i3(jElt),1:4))';LD2=[LD2 LD2(1)];
    uo = [Inf 1 0 Patch(i3(jElt),[6 5]) length(LD2)+6 LD2];
    if jDef>0 uo(1,3)=IIAxes(l1,jDef);end
    i6=length(find(sparse(LD2,1,LD2))); % number of distinct nodes

    if ua(j2,3)==0 & i6~=1 % creation
      i2(jElt)=patch('clipping','off', ...
       'xdata', Nc(LD2,1),'ydata', Nc(LD2,2), ...
       'zdata', jElt*ones(size(LD2)), ... % 'zdata', Nc(LD2,3), ...
       'userdata',uo);
    elseif ua(j2,3)==0 & i6==1 % creation of mass
       i2(jElt)=surface('clipping','off','marker','o','linestyle','none', ...
       'xdata',Nc(LD2,1),'ydata',Nc(LD2,2),'zdata',jElt*ones(5),'userdata',uo);
    elseif ua(j2,3)~=0 & i6~=1
      set(i2(jElt),'clipping','off', ...
       'xdata', Nc(LD2,1),'ydata', Nc(LD2,2), ...
       'zdata', jElt*ones(size(LD2)), ... % 'zdata', Nc(LD2,3), ...
       'userdata',uo);
    elseif ua(j2,3)~=0 & i6==1
      set(i2(jElt),'clipping','off',st3,'o', ...
       'xdata',Nc(LD2,1),'ydata',Nc(LD2,2),'zdata',jElt*ones(5),'userdata',uo);
    end

 i6=[i6 i2(jElt) remi(uo(4),size(cmap,1)) uo(1,3)];

 % DefColor: 0 Order 1 Element 2 Node 3 Group  4 Lighting 5 element nodes

 i4=rem(IIAxes(l1,24),10); % edge color
 if i4==0     i4 = 'none'; elseif i4==1 i4 = 'interp';
 elseif i4==2 i4 = 'flat'; elseif i4==3 i4 = 'w'; 
 elseif i4==4 i4 = 'k';    elseif i4==5 i4 = cmap(i6(4),:);
 else  i4=cmap(ua(j2,4),:); end

 st=fix(IIAxes(l1,24)/10); % face color
 if st==0     st = 'none'; elseif st==1 st = 'interp';
 elseif st==2 st = 'flat'; elseif st==3 st = 'w'; 
 elseif st==4 st = 'k';    elseif st==5 st = cmap(i6(4),:);
 else   st=cmap(ua(j2,4),:); end

 if i6(1)==1 % case of a single node
   if IIAxes(l1,26)==3       % group 
     set(i6(2),'cdata',[],'edgecolor',cmap(i6(3),:),'facecolor',cmap(i6(3),:))
   elseif IIAxes(l1,26)==2   % node 
     set(i6(2),'cdata',CDefN(LD2(1),i6(4))*ones(5),...
                      'edgecolor','flat','facecolor','flat');
   elseif IIAxes(l1,26)==1   % element
     set(i6(2),'cdata',CDefE(uo(5),i6(4))*eye(5),...
                      'edgecolor','flat','facecolor','flat');
   elseif IIAxes(l1,26)==4   % lighting
     set(i6(2),'cdata',[],'edgecolor','w','facecolor','w');
   elseif IIAxes(l1,26)==5   % element color
     set(i6(2),'cdata',[],'edgecolor','w','facecolor','w');
   else
     i5=cmap(remi(i6(4),size(cmap,1)),:);set(i6(2),'cdata',[],...
           'facecolor',i5,'edgecolor',i5)
   end

 elseif IIAxes(l1,26)==3       % group 

       i5=i4; if i6(1)<3 i5=cmap(i6(3),:);
       elseif comstr(i5,'interp')|comstr(i5,'flat') i5='k'; end
       set(i6(2),'cdata',[],'edgecolor',i5,'facecolor',cmap(i6(3),:))

 elseif IIAxes(l1,26)==2   % node 
      if i6(4)>size(CDefN,2) error('Not enough CDefN vectors');end

      i5=i4; if i6(1)<3 i5=st;end
      set(i6(2),'cdata',CDefN(LD2,i6(4)),'edgecolor',i5,'facecolor',st)

 elseif IIAxes(l1,26)==1   % element

   i5=i4; if i6(3)<3 i5=st;end
   set(i6(2),'cdata',ones(size(LD2))*CDefE(uo(5),i6(4)), ...
           'edgecolor',i5,'facecolor',st)

 elseif IIAxes(l1,26)==5   % element node
       if isempty(CDefE) CDefE=get(us(13),'userdata');
          if isempty(CDefE)
	   error('Use feplot(''InitCDef'') to define element node colors');
          end
       end
       i8=elt(uo(5),:);i8=i8(find(i8));i8=i8(1:length(i8)-2);
       i7=[];i7(i8)=1:length(i8);
       % find color value for the appropriate nodes
       set(i6(2),'cdata', ...
           CDefE(uo(5),i7(node(uo(7:uo(6)),1))), ...
           'edgecolor',i4,'facecolor','interp')

 elseif IIAxes(l1,26)==4 % lighting

       i7=Nc(LD2,1:3)-ones(length(LD2),1)*Nc(LD2(1),1:3);
       [i7,i8,i9]=svd(i7,0);
       r = diffuse(i9(1,3),i9(2,3),i9(3,3),IIAxes(l1,11:12));
       set(i6(2),'cdata',r*ones(size(LD2)), ...
           'edgecolor',i4,'facecolor','flat')
 else % default 
       if comstr(i4,'k')|comstr(i4,'w')|comstr(i4,'n')
       elseif i6(1)<3 i4=cmap(remi(i6(4),size(cmap,1)),:);
       else i4='k';end
       if comstr(st,'k')|comstr(st,'w')|comstr(st,'n')
       else st=cmap(remi(i6(4),size(cmap,1)),:);end
       set(i6(2),'cdata',[],'facecolor',st,'edgecolor',i4)
 end % of color selection

 end  % of loop jElt
ua(j2,3)=text('visible','off','userdata',i2);

end % matlab version selection


% sensor and arrowed sens object creation - - - - - - - - - - - - - - - - - -
elseif ua(j2,1)==3|ua(j2,1)==7

 i1=[];i1(nSel)=1:length(nSel); %eliminate sensors at unselected nodes
 if length(i1)<max(Sens(:,2)) i1(max(Sens(:,2)))=0;end
 i1=find(i1(Sens(:,2))); Sens=Sens(i1,:);
 LD2=Sens(:,2);LD2=[LD2 LD2 nInf*ones(size(LD2,1),1)]'; LD2=LD2(:);
 ua(j2,2)=[1.01]; uo=[Inf ua(j2,1) 0 i1(:)'];

 % current deformation at needed nodes with projection
 Nc=Nr;if jDef>0
    a0=real(IIAxes(l1,jDef+1)*exp(i*IIAxes(l1,jDef+2))* ...
       (Sens(:,3).*Def(Sens(:,2),IIAxes(l1,jDef))+ ...
        Sens(:,4).*Def(Sens(:,2)+nInf,IIAxes(l1,jDef))+ ...
        Sens(:,5).*Def(Sens(:,2)+2*nInf,IIAxes(l1,jDef))));
    a0=a0(:,[1 1 1]).*Sens(:,3:5);uo(3)=IIAxes(l1,jDef);
 else  a0=IIAxes(l1,23)*Sens(:,3:5);uo(3)=0;
 end

 aa=[node(LD2,5:7) ones(length(LD2),1)];
 aa(2:3:size(aa,1),1:3)=aa(2:3:size(aa,1),1:3)+a0;
 aa=aa*Xform';

 if ua(j2,1)==7
   % make the arrow
   aa=[aa(2:3:size(aa,1),1:2) aa(2:3:size(aa,1),1:2)-aa(1:3:size(aa,1),1:2)];
   aa(:,5:6)=aa(:,[3:4])*[0 -1;1 0]/20;
   aa=[aa(:,1:2)-aa(:,3:4) aa(:,1:2) aa(:,1:2)-aa(:,3:4)/10+aa(:,5:6) ...
     aa(:,1:2)-aa(:,3:4)/10-aa(:,5:6) aa(:,1:2) aa(:,1:2)*Inf];%
   i1=size(aa,1);
   aa=[reshape(aa(:,1:2:size(aa,2))',i1*size(aa,2)/2,1) ...
     reshape(aa(:,2:2:size(aa,2))',i1*size(aa,2)/2,1) ...
     ones(i1*size(aa,2)/2,1)];
 end

 if ua(j2,3)==0  ua(j2,3) =line('erasemode',emode,'clipping','off'); end
 set(ua(j2,3),'xdata',aa(:,1),'ydata',aa(:,2),'zdata',aa(:,3), ...
      'userdata',uo,'color',cmap(remi(ua(j2,4),size(cmap,1)),:));

else % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 disp(sprintf('%i unknown object type',ua(j2,1)));
end % of object type creation/update - - - - - - - - - - - - - - - - - -
end % of j2 loop on objects to be created - - - - - - - - - - - - - - - - -

if comstr(version,'5')
 set(gca,'userdata',ua,'view',[0 90],'DataAspectRatio',[1 1 1], ...
    'visible','off','xdir','normal','ydir','normal','zdir','normal');
else
 set(gca,'userdata',ua,'view',[0 90],'DataAspectRatio',[1 1 1], ...
    'visible','off', ...
    'xdir','normal','ydir','normal','zdir','normal');
end

if IIAxes(l1,27)~=0
  if size(Head,1)~=size(Def,2) fecom('head');Head = get(us(15),'userdata'); end
  i1=IIAxes(l1,26+[1:IIAxes(l1,30)]*5);
  st = Head(i1,:); if size(st,1)>1 st=[st setstr(32*ones(size(st,1),1))];end
  st=st';st=st(:)'; if length(st)>80 st=[st(1:75) ' ...'];end
  title(st);
  set(get(gca,'title'),'visible','on','fontsize',IIgui(3,2))
end

if ~isempty(st1)
 fecom(st1,';');fecom('setslider'); %comgui(sprintf('@Updated axis %i',l1));
end

end % of plot commands - - - - - - - - - - - - - - - - - - - - - - - - -

if nargout>0 out=ua;end

% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------
elseif comstr(Cam,'init') [CAM,Cam]=comstr(CAM,5);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if comstr(Cam,'model') [CAM,Cam]=comstr(CAM,6); % 'Initmodel'


% node initialization

 if nargin>1 node=CAM1;if isnumeric(node) node=full(node);end
 else node=[];end
 if size(node,1)==0
 elseif isstruct(node) 
  feplot('optcheck',2);
  feplot('initmodelback',node.Node,node.Elt);
  feplot('initdefBack',[],[]);
  fecom('GroupAll');
  return
 elseif size(node,2) == 3 %backward compatibility for xyz nodes
     node = [[1:size(node,1)]'*[1 0 0 0] full(node(:,1:3))];
 elseif size(node,2)~=7 error('unknown node format');  end


 if IIgui(6,8)==0  error('feplot InitNode: no current structure');
 else us=get(IIgui(6,8),'userdata'); end


 if ~isempty(node)
   [i1,i2]=sort(node(:,1));i1=[1;find(diff(i1))+1];
   if length(i1)~=size(node,1)
     disp(['feplot initnode: eliminated duplicated nodes' 7]);
     node=node(i2(i1),:);
   end
   if ~all(node(size(node,1),5:7)==[Inf Inf Inf])
    node = [node;max(node(:,1))+1 0 0 0 Inf Inf Inf];
   end
   set(us(3),'userdata',node);
 else node=get(us(3),'userdata'); end
 if isempty(node)
   error('feplot(''initmodel'') No declared nodes and no preexisting nodes');
 end

% element initialization

 if nargin>2 elt=full(CAM2); end
 if isempty(elt)
    i1=[node(:,1)*[1 0]]';i1=i1(:)';
    elt(1,[1 3:7 82+[1:length(i1)]])=[length(i1) abs('nodes') i1];
 end
 set(us(4),'userdata',elt); set(us(6:8),'userdata',[])

% group initialization

 [nElt,i1]=size(elt); group=[];

 if nElt>0&~isfinite(elt(1,1)) % elements have been defined
  EGroup = [find(elt(:,1)==Inf)' size(elt,1)+1];nGroup = length(EGroup)-1;
  for jGroup = 1:nGroup   group=[group '|' sprintf('%i',jGroup)]; end

 elseif nElt>0              % trace line plot
   for j1 = 1:size(elt,1)
        st=comstr(elt(j1,3:82),1);
        if isempty(st) st=sprintf('%i',j1); end
	group=[group '|' abs(st)];
   end
 else nGroup=0;EGroup=[];elt([1 83])=[1 0];group='|none|';  end

 group=[group '|']; set(us(5),'userdata',group);

 if ~comstr(Cam,'ba') fecom('groupall');end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'node') [CAM,Cam]=comstr(CAM,5);

 disp('You should use feplot(''initmodel'')');
 feplot('initmodel',CAM1,[]);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'elt') [CAM,Cam]=comstr(CAM,4);

 disp('You should use feplot(''initmodel'')');
 feplot('initmodel',[],elt);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Initialisation of a line1 property for a given structure
elseif comstr(Cam,'line') [CAM,Cam]=comstr(CAM,5);

 comgui('@Initialization of line plot');drawnow;
 if IIgui(6,8)==0 error('feplot InitElt: no current structure');
 else us=get(IIgui(6,8),'userdata'); end
 elt=get(us(4),'userdata');node=get(us(3),'userdata');
 LD2 = feutil(['Getline' CAM],node,elt); set(us(6),'userdata',LD2);
 i1 = LD2(:,2:size(LD2,2));i1=find(sparse(i1+1,1,i1))-1;
 i2 = [1:size(node,1)]';i2(i1)=i1*0;i2(length(i2))=0;i2=find(i2);
 if ~isempty(i2)
  comgui(sprintf('@Done (%i groups), unused nodes will not show',size(LD2,1)));
  fprintf(1,'\nUnused nodes %s will never show\n', ...
        sprintf(' %i',node(i2,1)))
 else
  comgui(sprintf('@Done InitLine (%i groups)',size(LD2,1)));
 end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Initialisation of a patch property for a given structure
elseif comstr(Cam,'patch') [CAM,Cam]=comstr(CAM,6);

 comgui('@Initialization of patch plot');drawnow;
 if IIgui(6,8)==0 return; else us=get(IIgui(6,8),'userdata'); end
 elt=get(us(4),'userdata');node=get(us(3),'userdata');
 LD2 = feutil('GetPatch',node,elt);
 for j1=2:4  
   i2=find(LD2(:,j1)==0); if ~isempty(i2) LD2(i2,j1)=LD2(i2,j1-1);end
 end

 set(us(8),'userdata',LD2);
 comgui(sprintf('@Done InitPatch (%i patches)',size(LD2,1)-1));

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Initialisation of a group property for a given structure
elseif comstr(Cam,'group') [CAM,Cam]=comstr(CAM,6);

 error('use initmodel');

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'sens')

% get info on the current structure
if IIgui(6,8)==0 return; else us=get(IIgui(6,8),'userdata'); end
elt=get(us(4),'userdata');node=get(us(3),'userdata');

% sens matrix now contains [SensID NodeIndex dx dy dz]

if isempty(CAM1) mdof=get(us(11),'userdata'); else mdof=CAM1; end
nInf=size(node,1);NNode(node(:,1))=[1:length(node(:,1))]';

if size(mdof,2)==5
  if max(abs(rem(mdof(:,2),1)))<1e-6 Sens=CAM1;% sensor matrix is  defined
  else mdof=mdof(:,1);end;
elseif isempty(mdof) error('empty DOF definition vector');
else mdof=mdof(:,1);end

if size(mdof,2)==1 % classical case where DOFs are defined

 % keep nodal translational DOFs only 
 i1=find(mdof>0);i2=round(rem(mdof(i1),1)*100);
 i2=[(i2==1)-(i2==7) (i2==2)-(i2==8) (i2==3)-(i2==9)];
 i3=find(any(i2'));i1=i1(i3);i2=i2(i3,:);
 Sens=[mdof(i1) fix(mdof(i1)) i2];
 if isempty(Sens) error('No translation DOF declared');end
end

% check on existing nodes
if max(Sens(:,2))>length(NNode) NNode(max(Sens(:,2)))=0;end
Sens(:,2)=NNode(Sens(:,2))';
i2 = find(~Sens(:,2)); if ~isempty(i2)
   disp('Removing sensors with undeclared node');
   disp(Sens(i2,:)); Sens=Sens(find(Sens(:,2)),:);
end

Line=get(us(6),'userdata');
i1 = Line(:,2:size(Line,2));i1=find(sparse(i1+1,1,i1))-1;
NNode=[];NNode(i1,1)=[1:length(i1)]';
if length(NNode)<max(Sens(:,2)) NNode(max(Sens(:,2)))=0;end
i2 = find(~NNode(Sens(:,2)));
if ~isempty(i2)
  comgui(sprintf('@Some sensors will never show'));
  fprintf(1, ...
    '\nSensors %s are connected to unused nodes and will never show\n', ...
    sprintf(' %.2f',Sens(i2,1)))
end

if ~isempty(strfind(Cam,'out')) out=Sens;
else set(us(9),'userdata',Sens);
 if isempty(strfind(Cam,'back')) feplot('plot'); end
end

% 'InitDef' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'def')

% get info on the current structure
if IIgui(6,8)==0 return; else us=get(IIgui(6,8),'userdata'); end
elt=get(us(4),'userdata');node=get(us(3),'userdata');

if nargin<2 mode=[]; else mode = CAM1; end
if nargin<3 mdof=[]; else mdof = CAM2; end
if isfield(mode,'def')&isfield(mode,'DOF')
 mdof=mode.DOF;mode=mode.def;
end

st='s';if length(mode)==0 & length(mdof)==0
   mode = ones(size(node,1)*3,1); st='n';
   mdof=[node(:,1)+.01 node(:,1)+.02 node(:,1)+.03]'; mdof=mdof(:);
elseif length(mode)==0   mode = ones(size(mdof,1),1);   end

% case with expansion

if size(mode,1)~=size(mdof,1) error('MODE and MDOF are not compatible'); end
if fix(max(mdof))>max(node(:,1)) error('fix(max(mdof))>max(node(:,1))');end
mdof=mdof(:,1);nInf=size(node,1);NNode(node(:,1))=[1:length(node(:,1))]';

Def  = zeros(3*nInf,size(mode,2));
Def([nInf 2*nInf 3*nInf],:)=Inf*ones(3,size(mode,2));
i1 = round(rem(mdof,1)*100); aa = [1 2 3 7 8 9;0 1 2 0 1 2;1 2 3 -1 -2 -3]';
for j1=1:6
 i2 = find(mdof>0 & i1==aa(j1,1));
 if ~isempty(i2)
   i3=NNode(fix(mdof(i2)));
   if aa(j1,3)>0 Def(i3+nInf*aa(j1,2),:) = mode(i2,:);
   else          Def(i3+nInf*aa(j1,2),:) = -mode(i2,:); end
 end
end

if size(mode,2)<5 & IIgui(8,1)~=1
   i3=IIAxes(l1,30+IIAxes(l1,30)*5+[1:IIAxes(l1,25)]); %selected group numbers
   IIAxes(l1,30+[0 1:5:size(mode,2)*5]) = [size(mode,2) 1:size(mode,2)];
   IIAxes(l1,30+IIAxes(l1,30)*5+[1:IIAxes(l1,25)])=i3;
   comgui(['@Deformation' sprintf(' %i',1:size(mode,2))]);
else
   if size(IIAxes,2)>30
    i1=31:min(size(IIAxes,2),30+IIAxes(l1,30)*5);
    IIAxes(l1,i1)=zeros(length(l1),length(i1));
   end
   IIAxes(l1,30:31) = [1 1]; comgui('@Deformation 1');
end

% set scale parameters end rows with [ValueAtScaleNode MaxDof]

i1=[1:nInf-1 nInf+1:2*nInf-1 2*nInf+1:3*nInf-1]; [i2,i3]=max(abs(Def(i1,:)));
for j1=1:length(i3) i2(j1)=Def(i1(i3(j1)),j1);end
Def(3*nInf+[1:2],:)=[i2;i1(i3)];

% clean up before exit
if ~any([1:nInf-1]==IIgui(8,4)) IIgui(8,4)=i3(1); end % ScaleDOFDflt
if IIgui(8,5)<2*eps % ScaleDefDflt
     IIgui(8,5)=1/max(max(abs(node(1:nInf-1,5:7))))/4;
end

i1=find(IIAxes(:,1)==3);
if ~isempty(i1)
 if st=='s' i2=find(IIAxes(i1,23)<2*eps);
            IIAxes(i1(i2),23)=ones(size(i2))*IIgui(8,5);
 else IIAxes(i1,23)=ones(size(i1))*eps; end
end

set(us(10),'userdata',full(Def)); set(us(11),'userdata',full(mdof))

if isempty(strfind(Cam,'back')) feplot('plotcreate'); end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'cdef') %'InitCdef'
[CAM,Cam]=comstr(CAM,5);

% get info on the current structure
if IIgui(6,8)==0 return; else us=get(IIgui(6,8),'userdata'); end
elt=get(us(4),'userdata');node=get(us(3),'userdata');
Def=get(us(10),'userdata');

cmode = CAM1;
if nargin>2 mdof = CAM2; else mdof=[];end
if ~isstr(cmode) & size(cmode,2)~=size(Def,2) & size(cmode,2)~=size(elt,2)
    disp('feplot InitCDef WARNING: not as many Cmodes as Modes')
end
nInf=size(node,1);  NNode(node(:,1))=[1:length(node(:,1))]';
% DefColor: 0 Order 1 Element 2 Node 3 Group 4 Lighting

% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if isempty(cmode)  disp('feplot InitCDef: empty color mode');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif isstr(cmode)

  CDefN=[];CDefE=[];
  if strcmp(cmode,'x')     CDefN=Def(nInf*0+[1:nInf],:);i1=[14 2];
  elseif strcmp(cmode,'y') CDefN=Def(nInf*1+[1:nInf],:);i1=[14 2];
  elseif strcmp(cmode,'z') CDefN=Def(nInf*2+[1:nInf],:);i1=[14 2];
  elseif strcmp(cmode,'a') CDefN=sqrt(Def(nInf*0+[1:nInf],:).^2+ ...
     Def(nInf*1+[1:nInf],:).^2+Def(nInf*2+[1:nInf],:).^2);i1=[14 2];
  elseif strcmp(cmode,'e') %element color to vertex color

   disp('Defining nodal colors based on element colors');
   CDefE = get(us(13),'userdata');
   %define nodal colors associated to the element colors
   EGroup=[find(~isfinite(elt(:,1)))' size(elt,1)+1];nGroup=length(EGroup)-1;
   CDefN = zeros(size(node,1),size(CDefE,2)); i4=CDefN;
   for jGroup = 1:nGroup %loop on element groups
     ElemF= feutil('getelemf',elt(EGroup(jGroup),:),jGroup);
     cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;
     i1=fe_super('node',ElemF); % node indices
     i2=NNode(elt(cEGI,i1)); for j1=1:size(CDefE,2)
        i3=CDefE(cEGI,j1*ones(size(i1)));
        i5 = full(sparse(i2(:),1,ones(size(i3))));
        i3 = full(sparse(i2(:),1,i3(:)));
        CDefN(1:length(i3),j1)=CDefN(1:length(i3),j1)+i3;
        i4(1:length(i3),j1)=i4(1:length(i3),j1)+i5;
     end
   end % of loop on groups
   i2=find(i4);CDefN(i2)=CDefN(i2)./i4(i2);
   i1= [14 2];
  elseif strcmp(cmode,'group')
     CDefE=[1:size(elt,1)]';i1=[24 4];
  else i1=[11 0];CDefE=0; end
  IIAxes(l1,[24 26])=i1;
  if ~isempty(CDefN) set(us(14),'userdata',CDefN);end
  if ~isempty(CDefE) set(us(13),'userdata',CDefE);end

% NodeColor - - - - - - - - - - - - - - - - - - - - - - - -
elseif size(cmode,1)==size(node,1) | size(cmode,1)==size(node,1)-1 | ...
   (isfield(cmode,'DOF')&round(rem(cmode.DOF(1),1)*100)==99)

  if isfield(cmode,'DOF') mdof=cmode.DOF;
   if isfield(cmode,'def') cmode=cmode.def;else cmode=cmode.data;end
  end

  CDefN(nInf,1)=0; CDefN(1:size(cmode,1),1:size(cmode,2))=cmode;
  set(us(14),'userdata',CDefN); IIAxes(l1,26)=2;

% NodeColor - - - - - - - - - - - - - - - - - - - - - - - -
elseif size(cmode,1)==length(mdof)|isfield(cmode,'DOF')

  if isfield(cmode,'DOF') mdof=cmode.DOF;
   if isfield(cmode,'def') cmode=cmode.def;else cmode=cmode.data;end
  end

  i1=find(mdof>0); if length(i1)<length(mdof)
    disp('feplot InitCdef: eliminated element DOFs');
    mdof=mdof(i1);cmode=cmode(i1,:);
  end
  i1=fix(mdof);[i1,i2]=sort(i1);cmode=cmode(i2,:);
  if any(diff(i1)==0) 
    disp('feplot InitCdef: more than 1 DOF per node, norm will be used');
    cmode=cmode.^2;CDefN=sparse(i1,1:size(cmode,2),cmode,nInf,size(cmode,2));
    CDefN=sqrt(full(CDefN));
  else
    CDefN=zeros(nInf,size(cmode,2)); CDefN(i1,:)=cmode;
  end
  set(us(14),'userdata',CDefN); IIAxes(l1,26)=2;

% EltColor - - - - - - - - - - - - - - - - - - - - - - - -
elseif size(cmode,1)==size(elt,1)|isfield(cmode,'EltId')

 if isfield(cmode,'EltId')
   eltid=feutil('eltidfix',elt);nind(eltid+1)=1:length(eltid);
   r1=zeros(size(elt,1),size(cmode.data,2));
   r1(nind(cmode.EltId+1),:)=cmode.data;
   set(us(13),'userdata',r1); 
 else set(us(13),'userdata',cmode); 
 end

 IIAxes(l1,26)=1;
 if ~isempty(CAM) feplot('initcdef','e');Cam='back';end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
else   disp('feplot InitCDef: type of color mode not supported');

end
if isempty(strfind(Cam,'back')) feplot('PlotCreate 4 1'); end 

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
else out='unknown'; end % of init subcommands - - - - - - - - - - - - - - - -

% ----------------------------------------------------------------------------
elseif comstr(Cam,'optcheck')

 opt=Inf*ones(7,7);
 if nargin>1 opt(1:size(CAM1,1),1:size(CAM1,2))=CAM1; end

 comgui('CheckGUI FE');
 if isfinite(opt(1,1)) IIgui(8,1)=opt(1,1); end
 if ~any([1:3 7]==IIgui(8,1))  % type of plot
  disp(sprintf('PlotType (%i) not known and set to 2',IIgui(8,1)));
  IIgui(8,1)=2;
 end

 % type of undeformed
 if isfinite(opt(1,2)) IIgui(8,2)=opt(1,2); else IIgui(8,2)=2;end
 if ~any([0 1 2 10 11 12 20 21 22]==opt(1,2))  IIgui(8,2)=2; end

 %nCycle
 if isfinite(opt(1,3)) IIgui(8,3)=opt(1,3); end;
 if IIgui(8,3)==0 IIgui(8,3)=25;end
 IIgui(8,3)=max(fix(real(IIgui(8,3))),2);

 if isfinite(opt(1,4)) IIgui(8,4)=opt(1,4); end
 if isfinite(opt(1,5)) IIgui(8,5)=opt(1,5); end

 if isfinite(opt(1,6)) IIgui(8,6)=opt(1,6); end
 if ~isfinite(opt(1,6))&IIgui(8,6)==0 IIgui(8,6)=11; end
 if ~any([0 1 2 3 4]==rem(IIgui(8,6),10))| ...  % face and edge colors
    ~any([0 1 2 3 4]==fix(IIgui(8,6)/10))  IIgui(8,6)=11; end

 if isfinite(opt(1,7)) IIgui(8,7)=opt(1,7); end  % type of EraseMode
 if ~any([0 1 2 ]==IIgui(8,7))  IIgui(8,7)=0; end

 %check of IIgui
 eval('get(IIgui(7,1),''type'');','IIgui(7,1)=0;')
 if all(IIgui(9,1:3)==0) IIgui(9,1:2)=[0 90];end
 if IIgui(7,1)==0
   if IIgui(1,1)==0 IIgui(1,1)=1;end;comgui('guife back');
 end

 %checks in the IIAxes variable

 if isempty(IIAxes) l1=[]; else l1 = find(IIAxes(:,1)==3); end
 if isempty(l1) iicom('sub ret1 1 1 3 2');end
 i1=[]; if size(IIAxes,2)>0 i1 = find(IIAxes(:,1)==3); end
 if ~isempty(i1) if ~any(IIAxes(i1,3)==double(gcf)) 
         IIAxes(i1(1),3)=2;figure(IIAxes(i1(1),3));end;
 end

 comgui(sprintf('CheckAxes 3 %i 0',double(gcf)));
 l1 = find(IIAxes(:,1)==3);
 if ~isempty(l1)
   if any(l1==IIgui(3,7)) l1=IIgui(3,7);else   l1=l1(1);IIgui(3,7)=l1; end
   IIAxes(l1,[2 19 21:24])=IIgui(8,[1:6]); 
 end

 %StructPtr : Object to store pointers to the structures
 up=findobj(2,'type','axes','tag','Struct1');
 if ~isempty(up) if any(IIAxes(:,4)==up) delete(up);up=[];end;end
 if isempty(up)  figure(2);
   up=axes('position',[0 0 eps eps],'visible','off','tag','Struct1');
 else axes(up); end
 IIgui(6,8)=up;
 
 us = get(up,'userdata');
 if length(us)<3   us(1:2) = [Inf 3.5]; end
 if length(us)<15  us(15) = [0]; end

 i1 = findobj(up,'type','text');
 for j1=3:15  i2 = find(i1==us(j1));
   if isempty(i2) us(j1)= text('visible','off'); else i1(i2)=0;end
 end
 i1 = i1(find(i1~=0));if ~isempty(i1) delete(i1);end
 set(up,'userdata',us);

% ----------------------------------------------------------------------------
elseif strfind(Cam,'supported')

out = [1 2 3 7 102];
out1 = ['|feplot ...|patch|line|sensor|arrow|head|'];

% ----------------------------------------------------------------------------
elseif ~isempty(CAM)

  if nargin<2 CAM1=[];end
  fecom(CAM,CAM1);

end % of selection of commands

eval('set(IIgui(5,2),''backgroundcolor'',[.706 .706 .706]);','');
