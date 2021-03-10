function [out,out1]=femesh(varargin)

%FEMESH	GUI commands for finite element mesh generation
%
%	Syntax: femesh
%		femesh CommandString
%		femesh('Command', ... )
%
%	The following a gives a summary of available commands. Refer to the
%       HTML manual for details (type sdtweb('femesh') at the MATLAB prompt).
%
%	; (command chaining), q (exit command mode) (see COMMODE for details)
%
%	Add [FEel(i) FEel(j), Sel]      : combine models
%	Add [Node, NodeNew]             : add new nodes
%	Add Test (NodeShift)    : femesh('addtest (n0)',node,ldraw)
%	Divide [(div1 div2 div3), InGroups, Group (i)]
%       DispFlag [on/off]  : level of infos returned
%	Extrude (nRep tx ty tz)
%	Find Dof (Element Selectors) : see FindElt for selectors
%	FindElt [[with,in] (node), mat (i), pro (i), group (i), egid(i)]
%		+ possibility to use &,|,{} and trailing numeric arguments
%	Find Node [Group (i), Groupa (i), Plane (Orig# nx ny nz), (x y z), ...
%                 rad <= r x y z, [x,y,z][>,<,==,>=,<=](i)] 
%               + possibility to combine with &| and trailing numeric arguments
%	Get Node [...] same as Find but returns nodes rather than node numbers
%	Info [ , FEel(i), node(i)]
%	Join [ElemtName, Group (i)]
%       Model  return model structure
%	ObjectBeamLine (i)
%	ObjectHoleInPlate CtrN# Edge1N# Edge2N# r1 r2 nDiv1 nDiv2 nQuadrant
%	ObjectMass m 
%	Optim [Model, Nodenum]
%	Orient (i) n (nx ny nz) 
%	Plot [Elt,El0]
%	quad2tria, quad42quadb, hexa82hexa20, tria32tria6, penta62penta15
%                hexa2tetra, hexa2penta
%	RemoveEl[t,0] (FindElt Selectors)
%	RepeatSel (nITE tx ty tz Orig# AnleDeg nx ny nz)
%	Refine beam (l)
%	Rev (nRep Orig# AngleDeg nx ny nz tx ty tz)
%	RotateSel (Orig# AnleDeg nx ny nz)
%	Sel Edge [,{El0 selection}] [Line, Patch] [,(element selectors)]
%	SelElt (calls FindElt and puts result in FEel0)
%	SelGroup (i)
%	SetGroup[ ,a] (i) [Mat (i), Pro (i), Name (s)]
%	SymSel (Orig# nx ny nz)
%	TransSel (tx ty tz)
%	UnJoin (i j) or femesh('unjoin','EltSel1','EltSel2')
%
%	Some commands can be used with FEMESH input arguments as detailed in
%	the manual.
%
%	FEMESH uses a number of global variables
%	  FEnode/FEn0/FEn1    main/selected/alternate node set
%	  FEelt/FEel0/FEel1   main/selected/alternate model description matrix
%       FEMESH automatically checks that these variables are declared as
%	global in your base and caller workspaces.
%
%       See sdtweb      fem, femesh
%       See also help   feplot, fecom, fe_mk, feutil. 
%		         demos  gartfe, d_truss, d_ubeam, ...

%       Etienne Balmes
%       Copyright (c) 2001-2020 by INRIA and SDTools
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       All Rights Reserved.

%#ok<*NOSEM,*NASGU>

persistent FE FEInfo DispFlag UseLegacy
% if 1, Use legacy formulation for test. See sdtweb femesh('TestUseLegacy')
if isempty(UseLegacy); UseLegacy=0; end 

global FEnode FEn0 FEn1 FEelt FEel0 FEel1

if nargin==1 && comstr(varargin{1},'cvs')
 out='$Revision: 1.184 $  $Date: 2020/12/16 20:44:41 $'; return;
end

epsl=sdtdef('epsl');
if isempty(DispFlag); DispFlag=1; end

% Provide a model to fill in the global variables 
if nargin>0 && ~ischar(varargin{1})
   
     evalin('caller', ...
        'iigui({''FEnode'',''FEn0'',''FEn1'',''FEelt'',''FEel0'',''FEel1''},''caller'');');
   r1=varargin{1};
   if isfield(r1,'Node')&&isfield(r1,'Elt')
     FE=r1;if isa(FE,'v_handle');FE=FE.GetData;end
     FEnode=FE.Node;
     FEelt=FE.Elt;
   else;warning('Input argument to femesh must be a model');
   end
   FEn0=[];FEn1=[];FEel0=[]; FEel1=[]; carg=2;
   st={'bas','Stack','pl','il'};

% Rather obsolete INIT call (never documented)
elseif nargin==1&&ischar(varargin{1})&&strcmp(varargin{1},'init')

  warning('This call is obsolete and will be removed in the future');
  [CAM,Cam]=comstr(varargin{1},5); if isempty(CAM); CAM='FE';end 
  if evalin('base',sprintf('exist(''%s'',''var'')',CAM))
   FE=evalin('base',CAM);
  elseif ~isfield(FE,'BaseVar')
   FE=struct('Node',[],'Elt',[],'bas',[],'pl',[],'il',[]);
  end
  FE.BaseVar=CAM;
  assignin('base',FE.BaseVar,FE);
  return;

% Just return the model 
elseif nargin==0 && nargout==1

  if isa(FE,'v_handle');FE=FE.GetData;end
  out=FE;  out.Elt=FEelt; 
  if ~isempty(FEel0);out.El0=FEel0;
  else;out=feutil('rmfield',out,'El0');
  end
  out.Node=FEnode; return;

else
   carg=1;
   if isfield(FE,'BaseVar') % new format with persistent variable
     FE=evalin('base',FE.BaseVar);
     FEnode=FE.Node;FEelt=FE.Elt;
     if isfield(FE,'El0'); FEel0=FE.El0; else;FEel0=[]; end
     if isfield(FE,'N0'); FEn0=FE.N0; else;FEn0=[]; end
   else  % Standard attempt to use the local values of caller as global
    if sdtdef('isdeployed'); st='commode';
    else  
      r1=dbstack;if length(r1)>1; [wd,st]=fileparts(r1(2).name); else;st='';end
    end
    if nargin==0 || ~comstr(st,'commode')
     evalin('caller', ...
        'iigui({''FEnode'',''FEn0'',''FEn1'',''FEelt'',''FEel0'',''FEel1''},''caller'');');
    end
   end
end  % How to deal with global(s)
if isa(FE,'v_handle');FE=FE.GetData;end

if nargin<carg && nargout==0 
 return; 
elseif isstruct(varargin{carg})
elseif ischar(varargin{carg}) 
 [CAM,Cam]=comstr(varargin{carg},1);carg=carg+1;
 if ~comstr(Cam,'dispflag') && (DispFlag==1)
  DispFlag=varargin{nargin}; if ischar(DispFlag)&&DispFlag(end)==';';DispFlag=0;
  elseif CAM(end)==';'; DispFlag=0; else;DispFlag=1; end
 end
end

if ~isempty(FEel0)&&isfinite(FEel0(1)); FEel0=[];end

out='done'; yPlot=0; CAM0=CAM;

%% #; commode separator ------------------------------------------------------
if comstr(Cam,';')

  commode('femesh',CAM(2:length(CAM)))

%% #DispFlag -----------------------------------------------------------------
elseif comstr(Cam,'dispflag'); [CAM,Cam] = comstr(CAM,9);

if isempty(Cam); out=DispFlag; return; end
if      comstr(Cam,'on');  DispFlag=1.;
elseif  comstr(Cam,'off'); DispFlag=0.;
else;   error('Command required: outputinfo on/off');
end

%% #Add ----------------------------------------------------------------------
elseif comstr(Cam,'add'); [CAM,Cam] = comstr(CAM,4);

%% #AddNode: ('addnode',OldNode,NewNode) - - - - - - - - - - - - - - - - - - -
if comstr(Cam,'node')

 if carg<=nargin; r1=varargin{carg};carg=carg+1;
  if carg<=nargin; r2=varargin{carg};carg=carg+1; end
 else 
  if DispFlag; disp('Adding FEn0 to FEnode'); end
  r1 = FEnode; r2=FEn0; 
 end

 [out,out1]=feutil(['add' Cam],r1,r2);
 if nargin<3; FEnode = out; out = 'done'; end

%% #AddFeel : Add element models - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'feel')

i1 = max([size(FEelt,2) size(FEel0,2) size(FEel1,2) ]);
if ~isempty(FEelt)&&size(FEelt,2)<i1; FEelt(1,i1)=0; end
if ~isempty(FEel0)&&size(FEel0,2)<i1; FEel0(1,i1)=0; end
if ~isempty(FEel1)&&size(FEel1,2)<i1; FEel1(1,i1)=0; end

i1=1; i3=[];
while ~isempty(i1)
   if Cam(i1+4)=='t'     ; i2='FEelt';
   elseif Cam(i1+4)=='0' ; i2='FEel0';
   elseif Cam(i1+4)=='1' ; i2='FEel1';
   end
   if i1==1; i3=[i2 '= [' i2]; else;i3 = [i3 ';' i2]; end
   i1 = i1+4;
   i1 = i1+strfind(Cam(i1:length(Cam)),'feel')-1;
end
i3 = [i3 '];']; eval(i3);

%% #AddSel : selection to main model - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'sel')

   i1 = max(size(FEelt,2),size(FEel0,2));
   if isempty(FEelt); FEelt=FEel0;
   elseif size(FEel0,1)>1
    FEelt(size(FEelt,1)+[1:size(FEel0,1)],1:size(FEel0,2))=FEel0;
   else % check for single superelement xxx
    disp('AddSel: empty selection');
   end

   if DispFlag; femesh('infoFEelt'); end

% ----------------------------------------------------------------------------
%% #AddTest - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% femesh('addtest (NodeIDShift)',node,ldraw);
elseif comstr(Cam,'test'); [CAM,Cam]=comstr(Cam,5);

r1 = varargin{carg};carg=carg+1;
if ~isfield(r1,'Node'); r1=struct('Node',r1);end
if ~isfield(r1,'Elt');  r1.Elt=varargin{carg};carg=carg+1;end

model=struct('Node',FEnode,'Elt',FEelt);

[model,i1]=feutil(['AddTest' CAM],model,r1);
FEnode=model.Node; FEelt=model.Elt;
if nargout>0; out=i1; end % return new nodes

%%
else;out='unknown';
end % subcommand selection - - - - - - - - - - - - - -

% ----------------------------------------------------------------------------
%% #_Divide : mesh refinement functions DIVIDE -------------------------------
elseif comstr(Cam,'divide') ;  [CAM,Cam]=comstr(Cam,7);

%% #DivideInGroups  - - - - - - - - - - - - - - - - - - - - - - - -
if comstr(Cam,'in')

FEel0=feutil('divideingroups',FEnode,FEel0,[]);
 
%% #DivideBy : 'divide by face  - - - - - - - - - - - - - - NOT DOCUMENTED
% divides by finding sharp edges in a group of quad4 elements
elseif comstr(Cam,'by'); [opt,CAM,Cam]=comstr(CAM,'by','%i');

[EGroup,nGroup]=getegroup(FEel0);
NNode(FEnode(:,1))=[1:length(FEnode(:,1))]'; elt=[];

for jGroup = 1:nGroup %loop on element groups
   [ElemF,i1,ElemP]= feutil('getelemf',FEel0(EGroup(jGroup),:),jGroup);
   if comstr(ElemP,'quad4')
    cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;
    elt(end+1,1:6) = [Inf abs('beam1')];
    i1=[1:4 1];i2=[];i3=[];i4=[];
    for j2 = 1:length(i1)-1
        i2=[i2 [1 length(NNode)+1]*sort(NNode(FEel0(cEGI,i1(j2+[0 1]))'))];
        i3=[i3 1:length(cEGI)];i4(j2*length(cEGI)+[-length(cEGI)+1:0])=j2;
    end
    r3=sparse(i3,i2,i4);
    
    i3=find(i2==1)'; elt(end+[1:length(i3)],1:2)= ...
         [rem(i3,length(NNode)+1) round(i3/length(NNode)+1)];

    for jElt=1:length(cEGI)
      [r2,r4]=basis(FEnode(NNode(FEel0(cEGI(jElt),1:4)),5:7));
      r1(jElt,1:3)=r2(:,3)';
    end
    for j2=find(any(r3))  % loop on edges
       i2=find(r3(:,j2));
       if length(i2)==1 || ( length(i2)>1 && abs(r1(i2(1),:)*r1(i2(2),:)')<.5)
         elt(end+1,1:2)=FEel0(cEGI(i2(1)),i1(r3(i2(1),j2)+[0 1]));
       end
    end
   end % type of element
end % of loop on groups
FEel0=elt;
 

%% #divideGroup :('divide group (i) (element selectors)')  - - - - - - - - - -
elseif comstr(Cam,'group') 

 FEelt=feutil(['divide' CAM],FEnode,FEelt,FEel0,varargin{2:end});

%% #Divide : Divide div1 div2 div3 - - - - - - - - - - - - - - - - - - - - - -
else

  st=['divideelt' CAM]; if DispFlag; st=['DivideElt ' CAM ' -disp '];end
  r1=feutil(st,struct('Node',FEnode,'Elt',FEel0),varargin{carg:end});
  FEnode =r1.Node;	FEel0 =r1.Elt;   yPlot=2;

end % subcommand selection - - - - - - - - - - - - - - - - - - - - - -

%% #Extrude ------------------------------------------------------------------
elseif comstr(Cam,'extrude'); [CAM,Cam]=comstr(CAM,8);

% extrude nRep tx ty tz
opt = comstr(CAM,[-1 1 0 0 1]);
femesh(sprintf('Rev %i 0 0 0 0 0 %.20g %.20g %.20g',opt(1:4)),varargin{2:end});

%% #Find ------------------------------------------------------------------------
elseif comstr(Cam,'find');  [CAM,Cam]=comstr(CAM,5);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% FindElt (selectors)
if comstr(Cam,'el')

if     comstr(Cam,'elt')
  [out,out1]=feutil(['find' CAM],FEnode,FEelt,FEel0,varargin{2:end});
elseif comstr(Cam,'el0')
  [out,out1]=feutil(['find' CAM],FEnode,FEel0,FEelt,varargin{2:end});
else;error('Not a valid element description matrix');
end


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% #FindSym[Plane Orig# nx ny nz] symmetry condition - - - - - - - - - - - - -
elseif comstr(Cam,'sym'); [CAM,Cam]=comstr(CAM,4);
if comstr(Cam,'plane'); [CAM,Cam]=comstr(CAM,6);

  if comstr(Cam,'o')
    opt=comstr(CAM(2:length(CAM)),[-1 0 0 0 0 0 0]);
    i1=opt(1:3); opt=opt(3:6);
  else
    opt =  comstr(CAM,[-1 0 0 0 0]);i1=find(FEnode(:,1)==opt(1));
    if isempty(i1); i1=[0 0 0]';else;i1=FEnode(i1(1),[5:7])'; end
  end
  if norm(opt(2:4))==0; error('FindSymPlane: you must define nx ny nz');end
  i1(:,2)=opt(2:4)'/norm(opt(2:4));
  [EGroup,nGroup]=getegroup(FEelt);
  NNode(FEnode(:,1))=[1:length(FEnode(:,1))]';

  %i2 contains the Nodes of interest, i3 the DOFs of interest
  i2 = femesh(sprintf('FindNode Group %i:%i',1,nGroup'));
  i6=[eye(3) i1(:,1);0 0 0 1]* ...
     [eye(3)-2*i1(:,2)*i1(:,2)' zeros(3,1);zeros(1,3) 1]* ...
     [eye(3) -i1(:,1);0 0 0 1];
  
  i6 = (i6*[FEnode(NNode(i2),5:7) ones(length(i2),1)]')';
  [FEnode,i3]=feutil('AddNode',FEnode,i6(:,1:3));i4=[];i4(i3)=[1:length(i3)]';
  if any(~i4(i2)); error('Some nodes have no symmetric equivalent');end

  i5 = femesh(sprintf('FindDof Group %i:%i',1,nGroup'));
  out=[]; % contains the symmetry condition matrix

  for j1=1:length(i2)
     if i2(j1)==i3(j1)
     elseif i2(i4(i2(j1)))>i2(j1)
        out=[out;fe_c(i5,[i2(j1)])-fe_c(i5,i2(i4(i2(j1))))];
     end
  end
  out=out';
end % of plane

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% #FindDof determines DOFs in a set of element groups
elseif comstr(Cam,'dof'); [CAM,Cam]=comstr(CAM,4);

 if carg<=nargin&&isa(varargin{carg},'double')&&~isfinite(varargin{carg}(1))
   elt=varargin{carg};carg=carg+1;
 else;elt=FEelt;
 end

 if isempty(elt);error('FindDof : Element definition matrix is empty');end
 sdtw('femesh FindDof is obsolete use feutil(''GetDof'',model) instead');
 out=feutil(['getdof ' CAM],FEnode,elt,FEel0,varargin{carg:end});

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% #FindNode'
elseif comstr(Cam,'node'); [CAM,Cam]=comstr(CAM,5);

if comstr(Cam,'get')
 [out1,out]=feutil(['findnode' CAM(4:end)],FEnode,FEelt,FEel0,varargin{2:end});
elseif nargout<2
  out=feutil(['findnode' CAM],FEnode,FEelt,FEel0,varargin{2:end});
else
  [out,out1]=feutil(['findnode' CAM],FEnode,FEelt,FEel0,varargin{2:end});
end

else;out='unknown';
end % subcommand selection - - - - - - - - - - - - - -

%% #Get ----------------------------------------------------------------------
elseif comstr(Cam,'get');  [CAM,Cam]=comstr(CAM,4);

if comstr(Cam,'node') 
  [out1,out]=feutil(['find' CAM],FEnode,FEelt,FEel0,varargin{2:end});
else;out='unknown';
end % subcommand selection - - - - - - - - - - - - - -

%% #Help : help commands -----------------------------------------------------
elseif comstr(Cam,'help');  [CAM,Cam] = comstr(CAM,5);

help femesh

%% Info ----------------------------------------------------------------------
elseif comstr(Cam,'info');  [CAM,Cam]=comstr(CAM,5);


if comstr(Cam,'node')
  feutil(['info' Cam],FEnode,FEelt);
elseif comstr(Cam,'el0'); feutil('infoeltFEel0',FEel0);
elseif comstr(Cam,'el'); [opt,CAM,Cam]=comstr(CAM,'eltid','%i');

  if carg<=nargin; opt=varargin{carg};carg=carg+1;
  else;opt=[]; end
  if isempty(opt); feutil('infoeltFEelt',FEelt)
  else;feutil('eltidinfo',FEelt,opt); end

else

if comstr(Cam,'feel') % info on a particular element definition matrix

  if     Cam(5)=='0'; feutil('infoeltFEel0',FEel0);
  elseif Cam(5)=='1'; feutil('infoeltFEel1',FEel1);
  elseif Cam(5)=='t'; feutil('infoeltFEelt',FEelt);
  else;feutil('infoeltFEelt',FEelt);feutil('infoeltFEel0',FEel0);
  end
else;feutil('infoeltFEelt',FEelt);feutil('infoeltFEel0',FEel0);
end

if DispFlag
 disp(sprintf('\nFEnode contains %i nodes, FEn0 contains %i nodes', ...
     size(FEnode,1),size(FEn0,1)));
end

end

%% #Join ---------------------------------------------------------------------
elseif comstr(Cam,'join'); [CAM,Cam]=comstr(CAM,5);

if comstr(Cam,'el0'); isel0=1; elt=FEel0; [CAM,Cam]=comstr(CAM,4);
else;elt=FEelt; isel0=0; end

elt=feutil(['join' CAM],elt);

if isempty(elt)
else
  if isel0; FEel0=elt; if DispFlag; femesh('infoFEel0'); end
  else;    FEelt=elt; if DispFlag; femesh('infoFEelt'); end;end
end

%% #Model : return model structure -------------------------------------------
elseif comstr(Cam,'model') ; [CAM,Cam] = comstr(CAM,6);

out=femesh;
if comstr(Cam,'0'); out.Elt=FEel0;end
out=feutil('rmfield',out,'El0');

% ----------------------------------------------------------------------------
%% #Object : standard object library -----------------------------------------
elseif comstr(Cam,'object'); [CAM,Cam] = comstr(CAM,7);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 'ObjectHoleInPlate'
% CenterNode# Edge1Node# Edge2Node# r1 r2 nDiv1 nDiv2 nQuadrant 
if comstr(Cam,'holeinplate'); [CAM,Cam]=comstr(CAM,12);

model=feutil(['ObjectHoleinPlate',CAM],struct('Node',FEnode,'Elt',[]));
FEnode=model.Node; FEel0=model.Elt;
yPlot=2;
    
else % objects in feutil

 out=feutil(['Object',CAM],struct('Node',FEnode,'Elt',[]), ...
     varargin{carg:end});
 if isnumeric(out)&&~isfinite(out(1)); FEel0=out;out='done';
 elseif nargout==0||isequal(varargin{end},';')||isequal(varargin{end},' '); 
     FEnode=out.Node; FEel0=out.Elt;
 end

end % subcommand selection - - - - - - - - - - - - - - -

% ----------------------------------------------------------------------------
%% #Optim :optimisation functions --------------------------------------------
elseif comstr(Cam,'optim'); [CAM,Cam] = comstr(CAM,6);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% #OptimModel: removes unused nodes from FEnode
if comstr(Cam,'model');   [FEnode,FEelt]=feutil('optimmodel',FEnode,FEelt);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% #OptimNodeNum: removes unused nodes from FEnode
elseif comstr(Cam,'nodenum')
  [FEnode,FEelt]=feutil('optimnodenum',FEnode,FEelt);

else;out='unknown';
end % subcommand selection - - - - - - - - - - - - - -


% ----------------------------------------------------------------------------
%% #Orient : Surface orientation ---------------------------------------------
elseif comstr(Cam,'orient'); [CAM,Cam]=comstr(CAM,7);

if comstr(Cam,'el0')
  FEel0=feutil(['orient' CAM(4:end)],FEnode,FEel0,[],varargin{carg:end});
else
  FEelt=feutil(['orient' CAM],FEnode,FEelt,FEel0,varargin{carg:end});
end


%% #Plot ---------------------------------------------------------------------
elseif comstr(Cam,'plot');  [CAM,Cam] = comstr(CAM,5);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% #PlotNode - - - -
if comstr(Cam,'node');  [CAM,Cam] = comstr(CAM,5);

 disp('Use the fecom(''nodetext'') command')

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% #PlotElt - - - -
elseif comstr(Cam,'el')

    
if comstr(Cam,'elt')
  model=femesh('model');st=sprintf('Plotting %s','FEelt');
elseif comstr(Cam,'el0')
  model=femesh('model0');st=sprintf('Plotting %s','FEel0');
else;error('Not a valid call to femesh PlotEl');
end
if comgui<105; feplot(model.Node,model.Elt,[],[],2);
else 
    r1=struct('FEnode',FEnode,'FEelt',FEelt,'FEel0',FEel0);
     if ~isfield(FEInfo,'Feplot')||~ishandle(FEInfo.Feplot)||...
      ~strcmp(get(FEInfo.Feplot,'tag'),'feplot')
      cf = sdth.urn('feplot'); FEInfo.Feplot=cf.opt(1);
     else;cf=feplot(FEInfo.Feplot);
     end
    FEnode=r1.FEnode;FEelt=r1.FEelt;FEel0=r1.FEel0;
    cf.model=model; 
end
  
else; error('Plot%s unknown',CAM);
end % subcommand selection - - - - - - - - - - - - - -

%% #Prop ---------------------------------------------------------------------
elseif comstr(Cam,'prop');  [CAM,Cam] = comstr(CAM,5); %#ok<*ASGLU>

 r1=struct('mdl',femesh,'name','femesh');
 eval('feutilg(''initmodel'',r1);');

% ----------------------------------------------------------------------------
elseif comstr(Cam,'quad42quad9')
   [FEnode,FEel0]=feutil('quad42quad9',FEnode,FEel0);

elseif comstr(Cam,'quad42q5p')
   [FEnode,FEel0]=feutil('quad42q5p',FEnode,FEel0);
elseif comstr(Cam,'quad42quadb')
    [FEnode,FEel0]=feutil('quad42quadb',FEnode,FEel0);
elseif comstr(Cam,'tria32tria6')
    [FEnode,FEel0]=feutil('lin2quad',FEnode,FEel0,[],'tria3');
elseif comstr(Cam,'hexa82hexa20')
    [FEnode,FEel0]=feutil('lin2quad',FEnode,FEel0,[],'hexa8');
elseif comstr(Cam,'penta62penta15')
    [FEnode,FEel0]=feutil('lin2quad',FEnode,FEel0,[],'penta6');
elseif comstr(Cam,'tetra42tetra10')
    [FEnode,FEel0]=feutil('lin2quad',FEnode,FEel0,[],'tetra4');
elseif comstr(Cam,'hexa82hexa27'); [CAM,Cam] = comstr(CAM,13);

  [EGroup,nGroup]=getegroup(FEel0);

  i2 =[]; % i2 first group nodes, i3 second group nodes
  i3=[1 2;2 3;3 4;4 1;1 5;2 6;3 7;4 8;5 6;6 7;7 8;8 5];
  i4=[1 2 3 4;1 4 8 5;1 2 6 5;5 6 7 8;2 6 7 3;3 7 8 4];
  for jGroup = 1:nGroup %loop on element groups
   ElemF= feutil('getelemf',FEel0(EGroup(jGroup),:),jGroup);
   if strcmp(ElemF,'hexa8')
     if comstr(Cam,'b'); FEel0(EGroup(jGroup),1:8)=[Inf abs('hexa27b')];
     else;              FEel0(EGroup(jGroup),1:7)=[Inf abs('hexa27')];
     end
     for j1=EGroup(jGroup)+1:EGroup(jGroup+1)-1
       NNode(FEnode(:,1))=[1:length(FEnode(:,1))]';
       i2 = FEnode(NNode(FEel0(j1,1:8)),5:7);
       i2=[(i2(i3(:,1),:)+i2(i3(:,2),:))/2;
           (i2(i4(:,1),:)+i2(i4(:,2),:)+i2(i4(:,3),:)+i2(i4(:,4),:))/4;
           mean(i2(1:8,:))];
       [FEnode,i2]=feutil('AddNode',FEnode,i2);
       FEel0(j1,1:29)=[FEel0(j1,1:8) FEnode(i2(:),1)' FEel0(j1,9:10)];
     end
   end
  end % of jGroup loop

% ------------------------------------------------------------------------
%% #hexa2tetra #hexa2penta #penta2tetra #refine #lin2quad
elseif comstr(Cam,'hexa2tetra')||comstr(Cam,'hexa2penta')||comstr(Cam,'hexa2pyra')||  ...
        comstr(Cam,'penta2tetra')||comstr(Cam,'refine')|| ...
        comstr(Cam,'lin2quad')
 model=feutil(CAM,FEnode,FEel0,[],varargin{carg:end});
 FEnode=model.Node; FEel0=model.Elt;

%% #quad2tria
elseif comstr(Cam,'quad2tria');  FEel0=feutil(CAM,FEel0);


% ----------------------------------------------------------------------------
%% #Remove : REMOVE functions ------------------------------------------------
elseif comstr(Cam,'remove'); [CAM,Cam] = comstr(CAM,7);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% #RemoveElt (findelt selectors)
if comstr(Cam,'elt') 
 FEelt = feutil(['remove' CAM],FEnode,FEelt,[],varargin{2:end});
elseif comstr(Cam,'el0') 
 FEel0 = feutil(['removeelt' CAM(4:end)],FEnode,FEel0,[],varargin{2:end});


else;out='unknown';
end % subcommand selection - - - - - - - - - - - - - -

% ----------------------------------------------------------------------------
%%#Repeat : repeat functions -------------------------------------------------
elseif comstr(Cam,'repeat'); [CAM,Cam] = comstr(CAM,7);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% #RepeatSel       nITE tx ty tz Orig# AnleDeg nx ny nz - - -
if comstr(Cam,'sel')  

model=feutil(['repeat', CAM],FEnode,FEel0,[],varargin{carg:end});
FEnode=model.Node; FEel0=model.Elt; yPlot=2;

%%
else;out='unknown';
end % subcommand selection - - - - - - - - - - - - - -

%% #Reset --------------------------------------------------------------------
elseif comstr(Cam,'reset'); [CAM,Cam] = comstr(CAM,7);
 FEnode=[]; FEelt=[];FEel0=[];
 FE=struct('Node',[],'Elt',[],'pl',[],'il',[]);
 FEInfo.Feplot=-1;
 if sdtdef('isdeployed');evalin('caller','global FEnode FEelt FEn1 FEel0');end
 
% ----------------------------------------------------------------------------
%% #Trans : translation functions --------------------------------------------
elseif comstr(Cam,'trans'); [CAM,Cam] = comstr(CAM,6);

% TransSel dx dy dz - - - - - - - - - - - - - - - - - - - - - - - - - - -
if comstr(Cam,'sel')||comstr(Cam,'el')

 model=feutil(['trans' CAM],FEnode,FEel0,[],varargin{carg:end});
 FEnode=model.Node; FEel0=model.Elt; yPlot=2;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% translation/rotation of selected group of nodes
elseif comstr(Cam,'node'); opt = comstr(CAM(5:length(CAM)),[-1 0 0 0]);

if DispFlag
 disp(['Translating FEn0 ' sprintf('tx %15.3f ty %15.3f tz %15.3f',opt)])
end
FEn0(:,5:7) = FEn0(:,5:7)+opt(ones(size(FEn0,1),1),1:3);

else;out='unknown';
end % subcommand selection - - - - - - - - - - - - - -

%% #Rev : REVOLUTION functions REV -------------------------------------------
elseif comstr(Cam,'rev'); 

 model=feutil(CAM,FEnode,FEel0,[],varargin{carg:end});
 FEnode=model.Node; FEel0=model.Elt; yPlot=2;

%% #RotateSel Orig# AnleDeg nx ny nz -----------------------------------------
elseif comstr(Cam,'rotatesel');

 model=feutil(CAM,FEnode,FEel0,[],varargin{carg:end});
 FEnode=model.Node; FEel0=model.Elt; yPlot=2;

% ----------------------------------------------------------------------------
%% #Sel : selection functions SEL --------------------------------------------
elseif comstr(Cam,'sel'); [CAM,Cam] = comstr(CAM,4);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% #SelElt - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if comstr(Cam,'el')

 [i4,FEel0]=femesh(['Find' CAM],varargin{2:end});

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% #SelGroup - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'group'); [CAM,Cam] = comstr(CAM,6);

 if comstr(Cam,'a')
   [EGroup,nGroup]=getegroup(FEel0);
   [CAM,Cam] = comstr(CAM,2);
   i1 = comstr(CAM,-1); 
   if isempty(i1)&&nargin>=carg; i1=varargin{carg};carg=carg+1; end
   if isempty(i1); error('No group selection specified'); end
   if any(i1<0|i1>nGroup); error('Invalid group selection'); end
   i2=[];for j1 = i1; i2 = [i2 EGroup(j1):EGroup(j1+1)-1]; end
   FEel0 = FEel0(i2,:);
   st=['FEel0 contains FEel0 group(s)' sprintf('  %i',i1)]; 
   if DispFlag; disp(st); end
 else
   [EGroup,nGroup]=getegroup(FEelt);
   i1 = comstr(CAM,-1);
   if isempty(i1)&&nargin>=carg; i1=varargin{carg};carg=carg+1; end
   if isempty(i1); error('No group selection specified'); end
   if any(i1<0|i1>nGroup); error('Invalid group selection'); end
   i2=[];for j1 = i1; i2 = [i2 EGroup(j1):EGroup(j1+1)-1]; end
   FEel0 = FEelt(i2,:);
   st=['FEel0 contains FEelt group(s)' sprintf('  %i',i1)]; 
   if DispFlag; disp(st); end
 end

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% #SelNode : finds the node closest to the current point - - - - - - - - - - 
elseif comstr(Cam,'node')

   NNode(FEnode(:,1))=[1:length(FEnode(:,1))]';
   if ~isempty(Cam) % selection of nodes given by number 
     i1 = comstr(CAM(5:length(CAM)),[-1]); % indices in node numbers
     i1 = i1(i1<=length(NNode)&i1>0);
     i1 = NNode(i1);	  % indices in position number
     i1 = i1(i1~=0);
   else;error('not available');
   end
   if ~isempty(i1)
    FEn0 = FEnode(i1,:);
     st=['Selected nodes (stored in FEn0) ' 10  ...
         sprintf('%6i%6i%6i%6i%6i%6i%6i%6i\n',FEnode(i1,1))];
   else;st='No node selected';
   end
   if DispFlag; disp(st); end

%%
else;out='unknown';
end % subcommand selection - - - - - - - - - - - - - -

% ----------------------------------------------------------------------------
%% #Set functions ------------------------------------------------------------
elseif comstr(Cam,'set');	[CAM,Cam] = comstr(CAM,4);

% #SetGroup [i] Mat[j] Prop [k] Name [s] EGID (i) - - - - - - - - - - - - - -
if comstr(Cam,'group'); [CAM,Cam] = comstr(CAM,6);

if comstr(Cam,'a')&&~comstr(Cam,'all');
 elt=FEel0;st1='FEel0'; [CAM,Cam] = comstr(CAM,2);
else;elt = FEelt;st1='FEelt'; 
end

[elt,st]=feutil(['setgroup' CAM],FEnode,elt,[],st1);
if DispFlag ; disp(st); end;

if comstr(st1,'FEel0'); FEel0=elt; else;FEelt=elt; end

else;out='unknown';
end % subcommand selection - - - - - - - - - - - - - -

% ----------------------------------------------------------------------------
%% #Stick : pushes the nodes of FEn1 onto a surface of quad selected in FEel0
elseif comstr(Cam,'stick')

[EGroup,nGroup]=getegroup(FEel0);
[i1,r1] = feutil(sprintf('findnode group1:%i',nGroup),FEnode,FEel0,[]);
NNode=[];NNode(i1)=1:length(i1);
r1=r1(:,5:7);
if carg<=nargin&&isa(varargin{carg},'double')
  node=varargin{carg};carg=carg+1;
else;node=FEn1;
end
ind = find(isfinite(FEel0(:,1)));

for j1=1:size(node,1)
  r2 = ((r1-node(j1*ones(size(r1,1),1),5:7)).^2)*[1;1;1];
  [r3,i2]=min(r2);
  [i3,i4] = find(FEel0(ind,1:4)==i1(i2));
  [r3,i4]=min(sum(r2(NNode(FEel0(ind(i3),1:4)))'));i3=ind(i3(i4));
  r3 = r1(NNode(FEel0(i3,1:4)),:);
  r3=r3-ones(size(r3,1),1)*mean(r3);
  r4 = basis(r3(2,:)/norm(r3(2,:))+r3(3,:)/norm(r3(3,:)),r3(3,:));
  out1(j1,1:3)=r4(:,3)'*((node(j1,5:7)-r1(i2,:))*r4(:,3));
  node(j1,5:7) = node(j1,5:7)-r4(:,3)'*((node(j1,5:7)-r1(i2,:))*r4(:,3));
end
if nargout==0; FEn1=node; else;out=node; end

% ----------------------------------------------------------------------------
%% #SymSel Orig# nx ny nz ----------------------------------------------------
elseif comstr(Cam,'symsel'); 

model=feutil(CAM,FEnode,FEel0,[],varargin{carg:end});
FEnode=model.Node; FEel0=model.Elt;

%% #TestPatch ----------------------------------------------------------------
elseif comstr(Cam,'testpatch') 

out.Node=[1 0 0 0   0.  10. 0.
          2 0 0 0   0.  0.  0.
          3 0 0 0   10. 0.  0.
          4 0 0 0   10. 10. 0.
          5 0 0 0   2.  2.  0.
          6 0 0 0   8.  3.  0.
          7 0 0 0   8.  6.  0.
          8 0 0 0   4.  7.  0.   ];

out.Elt=[Inf abs('quad4');
       1 2 5 8 1 1
       2 3 6 5 1 1
       3 4 7 6 1 1
       4 1 8 7 1 1
       5 6 7 8 1 1 ];

% ----------------------------------------------------------------------------
%% #Test a series of test objects
elseif comstr(Cam,'test'); [CAM,Cam] = comstr(CAM,5); 
 
 if comstr(Cam,'uselegacy'); % #TestUseLegacy
  UseLegacy=varargin{carg}; carg=carg+1; return
 end
    

RunOpt=struct('Struct',0,'Div',0,'Back',0,'Plot',0);
i1=strfind(Cam,'struct');
if isempty(i1); RunOpt.Struct=0;
else;RunOpt.Struct=1;CAM(i1+[0:5])='';[CAM,Cam]=comstr(CAM,1);
end

[CAM,Cam,RunOpt.Back]=comstr('back',[-25 3],CAM,Cam);
[CAM,Cam,RunOpt.Plot]=comstr('plot',[-25 3],CAM,Cam);
[CAM,Cam,RunOpt.Clear]=comstr('clear',[-25 3],CAM,Cam);

FE=struct('Node',[],'Elt',[],'pl',[],'il',[]);
% 100 steel
% 101 water
% 102 air
FE.pl=m_elastic('dbval 100 Steel','dbval 101 Water','dbval 102 Air');
r1=FE.pl(FE.pl(:,1)==100,:);r1(1)=112;r1([3 5])=2*r1([3 5]);
FE.pl=m_elastic(FE.pl,'dbval 112',r1);

% 110 Shell 1cm
% 111 full solid
% 112 1x1 beam section
% 113 full 2D solid
FE.il=[110 fe_mat('p_shell',1,1) 0 0 0 .01 0]; % default shell thickness
FE.il= p_solid(FE.il,'dbval 113 Plane stress',[111 fe_mat('p_solid','SI',1) 0 0 0]);
FE.il= p_beam(FE.il,'dbval 112 rectangle .02 .02');

st={'2bay','ubeam','1dbar'};ElemF='';
for j1=1:length(st); if comstr(Cam,st{j1}); ElemF=st{j1}; ElemP=st{j1}; end; end

if isempty(ElemF) % match a standard element test
 st={'q4pb',[100 113],'q8pb',[100 113],'t3pb',[100 113],'t6pb',[100 113],...
     'q4p',[100 113],'q8p',[100 113],'t3p',[100 113], ...
     't6p',[100 113],'q5p',[100 113],'q9a',[100 113],...
     'hexa8b',[100 111],'hexa20b',[100 111],'hexa27',[100 111], ...
     'penta6b',[100 111],'penta15b',[100 111],'tetra4b',[100 111], ...
     'tetra10b',[100 111],'pyra5',[100 111],'pyra13',[100 111],...
     'hexa8',[100 111],'hexa20',[100 111],'penta6',[100 111], ...
     'penta15',[100 111],'tetra4',[100 111],'tetra10',[100 111],...
     'tria3',[100 110],'tria6',[100 110],'quad4',[100 110], ...
     'quadb',[100 110],'quadc',[100 110], ...
     'quad9',[100 110],'mitc4',[100 110], ...
     'dktp',[100 110], ...
     'bar1',[100 112],'beam1',[100 112], ...
     'flui4',[102 111],'flui6',[102 111],'flui8',[102 111],...
     'celas',[1 1],'cbush',[1 1],'mck3',[1 1]};
 st=reshape(st,2,length(st)/2)';
 if comstr(Cam,'list'); disp(st(:,1)); return; end
 ElemF=''; i5=[]; j1=1;
 for j1=1:size(st,1); 
   if comstr(Cam,st{j1,1}); ElemF=st{j1,1};i5=st{j1,2}; break; end; 
 end
 if ~isempty(ElemF); 
 else
   ElemF=sscanf(CAM,'%s',1);RunOpt.CAM=CAM(length(ElemF)+1:end);
   if ~exist(ElemF,'file');error('Not a femesh teststruct'); end
 end
 ElemP=feval(ElemF,'parent'); FEel0=[Inf abs(ElemF)];
end
% use integration rule for b elements
if any(strncmp(ElemF,{'hex','pen','tet'},3))&&(ElemF(end)=='b'||~UseLegacy);
   FE.il(FE.il(:,1)==111,4)=-3; 
end
if any(strncmp(ElemF,{'q4p','q8p','t3p','t6p'},3))&&(ElemF(end)=='b'||~UseLegacy);
   FE.il(FE.il(:,1)==113,5)=-3; 
end
%
div=strfind(Cam,'divide'); % Allow for divide definition in command
if ~isempty(div);
 i1=div;[div,i2,i3,i4]=sscanf(Cam(i1+6:end),'%i');
 CAM(i1:i1+4+i4)='';[CAM,Cam]=comstr(CAM,1);
 if length(div)<3;div(end+1:3)=3;end
end

if ~isempty(div);RunOpt.Div=1;
elseif any(strncmpi(ElemP,{'quad4','quadb','tria3','tria6'},5))
 if RunOpt.Struct; div=[6 6];else;div=[1 1];end % division default for 2D
elseif RunOpt.Struct; div=[3 3 3];
else;div=[1 1 1];
end
il=FE.il;
% - - - - - - - - - - - - - - - - - - - - - - - create elements
% - - - - - - - - - - - - - - -  - - - - - - - - - - 2D 
if comstr(ElemP,'quad4') || comstr(ElemP,'quadb')

  FEel0=[Inf abs('quad4')];
  FEnode=[1 0 0 0  0 0 0;2 0 0 0  0 1 0;3 0 0 0 1 1 0;4 0 0 0 1 0 0];
  FEel0(end+1,1:6)=[1 4 3 2 i5];
  femesh(sprintf('divide%i %i',div(1:2)),';');
  if comstr(ElemP,'quadb');  femesh('quad42quadb',';');  end
  femesh(sprintf('set group a1 name %s matid %i proid %i',ElemF,i5));
elseif comstr(ElemP,'tria3') || comstr(ElemP,'tria6')
  FEnode=[1 0 0 0  0 0 0;2 0 0 0  0 1 0;3 0 0 0 1 1 0;4 0 0 0 1 0 0];
  if RunOpt.Struct;  
   femesh('teststructquad4 divide 1 1 -back'); femesh('quad2tria');
   femesh(sprintf('divide%i %i',div(1:2)),';');
  else;     FEel0=[Inf abs('tria3') ;1 2 3 i5 0 ];
  end
  model=femesh('model0'); model.Elt(1,1:6)=[Inf abs('t3p') 0 0];
  FEel0=feutil('orient',model);FEel0(1,1:6)=[Inf abs('tria3')];

  if comstr(ElemP,'tria6');  femesh('tria32tria6'); end
  femesh(sprintf('set group a1 name %s matid %i proid %i',ElemF,i5));
elseif comstr(ElemP,'q5p')
  FEnode=[1 0 0 0  0 0 0;2 0 0 0  0 1 0;3 0 0 0 1 1 0;
          4 0 0 0  1 0 0;5 0 0 0 .5 .5 0];
  FEel0(end+1,1:7)=[1 4 3 2 5 i5];
elseif comstr(ElemP,'quad9')
  FEnode=[1 0 0 0  0 0 0;4 0 0 0  0 1 0;3 0 0 0 1 1 0;2 0 0 0 1 0 0;
  5 0 0 0  .5 .0 .0; 6 0 0 0 1. .5 .0;7 0 0 0 .5 1. 0.;8 0 0 0 .0 .5 0
  9 0 0 0 .5 .5 0];
  FEel0(end+1,1:11)=[1:9 i5];
% - - - - - - - - - - - - - - -  - - - - - - - - - - 1D 
elseif comstr(ElemP,'beam1')
  FEnode=[1 0 0 0  0 0 0;2 0 0 0  1 0 0;3 0 0 0  0 0 1];
  FEel0(end+1,1:5)=[1:2 i5 3];
  femesh(sprintf('divide%i',div(1)),';');
elseif comstr(ElemP,'mck3')
  FEnode=[1 0 0 0  0 0 0;2 0 0 0  1 0 0;3 0 0 0  0.5 0 0];
  FEel0(end+1,1:10)=[1:3 i5 3 10 10 1000 1000 .2];
  femesh(sprintf('divide%i',div(1)),';');
elseif comstr(ElemP,'celas')
  FEnode=[1 0 0 0  0 0 0;2 0 0 0  1 0 0];
  FEel0(end+1,1:4)=[1:2 i5 3];
  femesh(sprintf('divide%i',div(1)),';');
elseif comstr(ElemP,'cbush')
  FEnode=[1 0 0 0  0 0 0;2 0 0 0  1 0 0];
  FEel0(end+1,1:4)=[1:2 i5 3];
  femesh(sprintf('divide%i',div(1)),';');
% - - - - - - - - - - - - - - -  - - - - - - - - - - 3D 
elseif comstr(ElemP,'tetra4')
  if ~RunOpt.Struct
   FEnode=[1 0 0 0  0 0 0;2 0 0 0  1 0 0;3 0 0 0 1 1 0;4 0 0 0 0 0 1];
   FEel0(end+1,1:6)=[1 2 3 4 i5];
  else
   if ~isempty(div);femesh(sprintf('testhexa8 divide %i %i %i -back',div(1:3)));
   else;femesh('teststructhexa8 -back');
   end
   femesh hexa2tetra; 
  end
  FE.il=il;
  femesh(strcat('set group a1 name ',ElemF,sprintf(' mat %i pro %i',i5(1),i5(2))));
elseif comstr(ElemP,'tetra10')
  if ~RunOpt.Struct; femesh('test tetra4');femesh('tetra42tetra10');
  else
   femesh('testhexa8 divide 1 1 1 -back');
   femesh(sprintf('divide%i %i %i',div(1:3)),';');
   femesh hexa2tetra;  femesh('tetra42tetra10');
   FEel0(1,1:length(ElemF)+1)=[Inf abs(ElemF)];
  end
  FE.il=il;
elseif comstr(ElemP,'pyra5')
  if ~RunOpt.Struct
    FEnode=[1 0 0 0  0 0 0;2 0 0 0  1 0 0;3 0 0 0  1 1 0;
            4 0 0 0  0 1 0;5 0 0 0  .5 .5 1;];
    FEel0(end+1,1:7)=[1:5 i5];
  else
   if ~isempty(div);femesh(sprintf('testhexa8 divide %i %i %i -back',div(1:3)));
   else;femesh('teststructhexa8 -back');
   end
   femesh hexa2pyra; 
  end
  FEel0=feutil(sprintf('set group1 name %s mat %i pro %i ', ...
      ElemF,i5(1:2)),FEnode,FEel0);
  FE.il=il;
elseif comstr(ElemP,'pyra13')
 if ~RunOpt.Struct;   femesh('test pyra5');femesh('lin2quad');
  else
   femesh('testhexa8 divide 1 1 1 -back');
   femesh(sprintf('divide%i %i %i',div(1:3)),';');
   femesh('hexa2pyra;');femesh('lin2quad'); 
  end
  FEel0=feutil(sprintf('set group1 name %s mat %i pro %i ', ...
      ElemF,i5(1:2)),FEnode,FEel0);
  FE.il=il;
elseif comstr(ElemP,'penta6')
  if ~RunOpt.Struct
    FEnode=[1 0 0 0  0 0 0;2 0 0 0  1 0 0;3 0 0 0  1 1 0;
            4 0 0 0  0 0 1;5 0 0 0  1 0 1;6 0 0 0  1 1 1];
    FEel0(end+1,1:8)=[1:6 i5];
  else
   femesh('testhexa8 divide 1 1 1 -back');
   femesh(sprintf('divide%i %i %i',div(1:3)),';');
   femesh('hexa2penta;');
  end
  FEel0=feutil(sprintf('set group1 name %s mat %i pro %i ', ...
      ElemF,i5(1:2)),FEnode,FEel0);
  FE.il=il;
elseif comstr(ElemP,'penta15')
  if ~RunOpt.Struct;   femesh('test penta6');femesh('penta62penta15');
  else
   femesh('testhexa8 divide 1 1 1 -back');
   femesh(sprintf('divide%i %i %i',div(1:3)),';');
   femesh('hexa2penta;');femesh penta62penta15; 
  end
  FEel0=feutil(sprintf('set group1 name %s mat %i pro %i ', ...
      ElemF,i5(1:2)),FEnode,FEel0);
  FE.il=il;
elseif comstr(ElemP,'hexa8')
  FEnode=[1 0 0 0  0 0 0;2 0 0 0  1 0 0;3 0 0 0  1 1 0;4 0 0 0  0 1 0;
          5 0 0 0  0 0 1; 6 0 0 0  1 0 1;7 0 0 0  1 1 1;8 0 0 0  0 1 1];
  FEel0(end+1,1:10)=[1:8 i5];
  if RunOpt.Struct||RunOpt.Div; 
   femesh(sprintf('divide%i %i %i',div(1:3)),';');
  end
  FEel0=feutil(sprintf('set group1 name %s mat %i pro %i ', ...
      ElemF,i5(1:2)),FEnode,FEel0);
  FE.il=il;
elseif comstr(ElemP,'hexa20')
  femesh('testhexa8 divide 1 1 1');
  femesh(sprintf('divide%i %i %i',div(1:3)),';');
  femesh('hexa82hexa20');
  FEel0=feutil(sprintf('set group1 name %s mat %i pro %i ', ...
      ElemF,i5(1:2)),FEnode,FEel0);
  FE.il=il;
elseif comstr(ElemP,'hexa27')
  if ~RunOpt.Struct; femesh('test hexa8');femesh('hexa82hexa27');
  else
   femesh('testhexa8 divide 1 1 1');
   femesh(sprintf('divide%i %i %i',div(1:3)),';');
   femesh('hexa82hexa27');
   FEel0=feutil(sprintf('set group1 name %s mat %i pro %i ', ...
      ElemF,i5(1:2)),FEnode,FEel0);
  end
  FE.il=il;
% - - - - - - - - - - - - - - -  - - - - - - - - - - OTHERS
elseif comstr(ElemP,'2bay') % #Test2Bay
  FEnode=[1 0 0 0  0 0 0;2 0 0 0    0 1 0;
          3 0 0 0  1 0 0;4 0 0 0    1 1 0;
          5 0 0 0  2 0 0;6 0 0 0    2 1 0];
  FEel0=[Inf abs('beam1'); 1 3 1 1 0 0;2 4 1 1 0 0;3 4 1 2 0 0;1 4 1 3 0 0;
         Inf abs('beam1'); 3 5 1 1 0 0;4 6 1 1 0 0;5 6 1 2 0 0;3 6 1 3 0 0];
  FEelt=FEel0; 
  FE.pl= m_elastic(FE.pl,'dbval 1 Steel');
  if nargout==1
    i1=fe_mat('p_beam','SI',1);
    FE.il = [1 i1 5.0000e-09 5.0000e-09   5.0000e-09   2.0000e-05   % longerons
     2 i1     5.0000e-09   5.0000e-09   5.0000e-09   2.0000e-05	  % diagonals
     3 i1     5.0000e-09   5.0000e-09   5.0000e-09   2.0000e-05]; % battens
   out=femesh('Model'); out.DOF=feutil('getdof',out);
   out.name='2Bay'; out.unit='SI';
   if RunOpt.Plot; feplot(out);end; return;
  end

elseif comstr(ElemP,'ubeam') % #TestUBeam
 FEnode=[1 0 0 0  -.5 -.5 0;2  0 0 0  -.5+1/6 -.5 0;3 0 0 0  -.5 .5-1/6 0
         4 0 0 0  -.5+1/6 .5-1/6 0;5 0 0 0  -.5 .5 0;6 0 0 0 -.5+1/6 .5 0
         7 0 0 0 .5-1/6 .5 0;8 0 0 0 .5 .5 0;9 0 0 0 .5-1/6 .5-1/6 0
         10 0 0 0 .5 .5-1/6 0;11 0 0 0 .5-1/6 -.5 0;12 0 0 0 .5 -.5 0];
 FEelt=[Inf abs('quad4');4 6 5 3 1 1;9 10 8 7 1 1];
 FEel0=[Inf abs('quad4');1 2 4 3 1 1;11 12 10 9 1 1];
 femesh(';divide 5 1;addsel;');
 FEel0=[Inf abs('quad4');4 6 7 9 1 1];
 femesh(';divide 4 1;addsel;');
 femesh('join group 1:3');
 femesh(';selgroup1;extrude 10 0 0 .25;');
 FEelt=FEel0; 
 FEelt(2:end,1:8)=FEelt(2:end,[5:8 1:4]);
 FEelt=feutil('orient;',FEnode,FEelt);FEel0=FEelt;
 FE.pl=m_elastic('dbval 1 Steel');
 FE.il=p_solid('dbval 1 d3 -3');

 if nargout==1; 
  out=femesh('model'); out.name='UBeam'; out.unit='SI';
  if RunOpt.Plot;feplot(out);end;return; 
 end

elseif comstr(ElemP,'bar')  % #TestBar
error('Change name xxx');
  FEnode=[1 0 0 0    0   0  0;2 0 0 0    50  0  0];
  FEel0=[Inf abs('bar1')];
  FEel0(2,:)=[1 2 1 1 1];
  femesh(sprintf('divide%i',div(1)),';');
  FEelt=FEel0; %FEel0=[];
  FE.pl=m_elastic('dbval 1 Steel');
  FE.il= p_beam('dbval 1 rectangle 1 1');

else; 
 FE=feval(ElemF,['TestMesh' RunOpt.CAM]);FEnode=FE.Node;
 FEel0=FE.Elt; FEelt=FE.Elt;
end

% - - - - - - - - - - - - - - - - - - - - - - - - - end of create elements

if RunOpt.Back;
  out=struct('Node',FEnode,'Elt',FEel0,'pl',FE.pl,'il',FE.il);out1=[];
  out.name='ElemP'; out.unit='SI';
  return;
end

if nargout==1 && ischar(out)
 out=struct('Node',FEnode,'Elt',FEel0,'pl',FE.pl,'il',FE.il);
end
Rand=strfind(Cam,'rand'); % Allow internal node randomization
if ~isempty(Rand);
 i1=Rand;[Rand,i2,i3,i4]=sscanf(Cam(i1+4:end),'%g');
 CAM(i1:i1+2+i4)='';[CAM,Cam]=comstr(CAM,1);
else;Rand=0;
end

if RunOpt.Struct 
%% StructTests - - - - - - - - - - - - - - - - 

model=struct('Node',FEnode,'Elt',FEel0,'pl',FE.pl,'il',FE.il);
def=[];out1='def = fe_load(model);';
switch Cam(1:3)
case {'bea','bar'} % - - - - - - - - - - - - - - - - 1D

  model=femesh('test2bay'); model.DOF=[];
  if ~exist('CAM','var'); Cam='eig'; end
  if comstr(Cam,'bar1')
   model.Elt=feutil('set group 1:2 name bar1',model);
   model=fe_case(model,'fixdof','2D',.03);
  end
  % reassign new material to half the structure
  model.Elt(feutil('findelt withnode {x>1}',model),3)=112;

  model=fe_case(model,'fixdof','Edge','x==0');
  [m,k,mdof]=fe_mk(model);
  [Case,model.DOF]=fe_mknl('init',model);
  k1=fe_mknl('assemble',model,Case,1);  
  if ~isequal(mdof,Case.DOF); i1=fe_c(Case.DOF,mdof,'ind');k1=k1(i1,i1);end  
  if norm(diag(k)-diag(k1))/norm(diag(k))>sqrt(eps) 
   error('Problem in NL assembly');
  end

  if ~isempty(strfind(Cam,'eig')); def=fe_eig(model,[105 6 1e3 11]);
  else % load
    data=struct('sel','groupall','dir',[0 1 0]);
    model=fe_case(model,'FVol','Volume Load',data);
  end

case {'q4p','q8p','t3p','t6p','q5p','tri','qua'} % - - - - - - - - - - - - 2D

  if Rand; % Randomize node positions
   i1=find(FEnode(:,5)<1);
   FEnode(i1,5)=FEnode(i1,5).*(1+rand(size(i1))/div(1)*Rand);
   i1=find(FEnode(:,6)<1);
   FEnode(i1,6)=FEnode(i1,6).*(1+rand(size(i1))/div(2)*Rand);
  end
  model.Node=FEnode;model.Elt=feutil('orient;',model);
 
  switch Cam(1:3)
  case {'tri','qua'}
   if ~isempty(strfind(Cam,'_1'));     model.il(model.il(:,1)==110,3)=1; 
   elseif ~isempty(strfind(Cam,'_2')); model.il(model.il(:,1)==110,3)=2; 
   elseif ~isempty(strfind(Cam,'_5')); model.il(model.il(:,1)==110,3)=5; 
   end
  otherwise
   if ~isempty(strfind(Cam,'_0'));     model.il(model.il(:,1)==113,3)=0; 
   elseif ~isempty(strfind(Cam,'_1')); model.il(model.il(:,1)==113,3)=1; 
   elseif ~isempty(strfind(Cam,'_2')); model.il(model.il(:,1)==113,3)=2; 
   end
   if ElemF(end)=='b';model.il(model.il(:,1)==113,5)=-3;end
  end

  if nargin==2; model.pl=varargin{2}; end % deals with specified pl (eg. anisotropic)

  % attribute material 112 to 1/2 structure
  i2=feval(ElemF,'prop'); ind=feutil('findelt withnode {x>.51}',model);
  model.Elt(ind,i2(1))=112;

  model=fe_case(model,'fixdof','Edge','x==0');

  if comstr(Cam(3),'p') % 2D
   data=struct('sel','groupall','dir',[1 0 0]);
   model=fe_case(model,'FVol','Volume Load',data);
   data=struct('sel','y==1','eltsel','withnode {x>0.51}', ...
               'def',1,'DOF',.02,'type','edge');
   model=fe_case(model,'Fsurf','Surface load',data);
   data=struct('sel','y==1','eltsel','withnode {x>0.51}', ...
               'def',1,'DOF',.19,'type','edge');
   model=fe_case(model,'Fsurf','Pressure load',data);
  else % SHELL
   data=struct('sel','groupall','dir',[0 0 9.81]);
   model=fe_case(model,'FVol','Volume Load',data);
  end
 
  if ~isempty(strfind(Cam,'eig'));def=fe_eig(model,[105 10 1e3 11]);end

case {'tet','pen','hex','pyr'} % - - - - - - - - - - - - - - - - - - - - - 3D 

  model.Elt=feutil('orient;',model);

  if Rand; % Randomize node positions
   i1=find(FEnode(:,5)<1);
   FEnode(i1,5)=FEnode(i1,5).*(1+rand(size(i1))/div(1)*Rand);
   i1=find(FEnode(:,6)<1);
   FEnode(i1,6)=FEnode(i1,6).*(1+rand(size(i1))/div(2)*Rand);
   i1=find(FEnode(:,7)<1);
   FEnode(i1,7)=FEnode(i1,7).*(1+rand(size(i1))/div(3)*Rand);
  end
  FEnode(:,7)=FEnode(:,7)/3/max(FEnode(:,7));model.Node=FEnode;
  if ElemF(end)=='b'||comstr(ElemP,'pyr');
       model.il(model.il(:,1)==111,4)=-3;
  end

  % attribute material 112 to 1/2 structure
  ind=feutil('findelt innode {x>.51}',model);
  i2=feval(ElemF,'prop');
  model.Elt(ind,i2(1))=112;
  model=fe_case(model,'fixdof','Base','x==0');

  if ~isempty(strfind(Cam,'eig')); def=fe_eig(model,[105 10 1e3 11]);
  else % load
    data=struct('sel','groupall','dir',[0 0 9.81]);
    model=fe_case(model,'FVol','Volume Load',data);
    data=struct('sel','z==1/3','eltsel','withnode {x>0.51}', ...
               'def',1,'DOF',.01);
    model=fe_case(model,'Fsurf','Surface load',data);
    data=struct('sel','z==1/3','eltsel','withnode {x>0.51}', ...
                'def',1,'DOF',.19);
    model=fe_case(model,'Fsurf','Pressure load',data);
  end
 
case 'flu'   % - - - - - - - - - - - - - - - - - - - fluids

   model=femesh('model0');
   FEnode(:,7)=FEnode(:,7)/3/max(FEnode(:,7));model.Node=FEnode;

   [m,k,mdof]=fe_mk(model);
   [Case,model.DOF]=fe_mknl('init',model);
   k1=fe_mknl('assemble',model,Case,1);  
   if ~isequal(mdof,Case.DOF); i1=fe_c(Case.DOF,mdof,'ind');k1=k1(i1,i1);end  
   if norm(diag(k)-diag(k1))/norm(diag(k))>sqrt(eps) 
    error('Problem in NL assembly');
   end
   if isempty(strfind(Cam,'eig')) && isempty(strfind(Cam,'load'))
     Cam='eig';
   end
   if ~isempty(strfind(Cam,'eig')) % eig
    model = fe_case(model,'FixDof','clamped dofs','z==0');
    def=fe_eig(model,[105 10 1e3 11]);
   end


otherwise % - - - - - - - - - - - - - - - - - - - Others
end


 % the same for all elements
 if ~RunOpt.Back&&~(RunOpt.Struct&&nargout==2)&&~isempty(def); 
     feplot(model,def);sdtw('_ewt','Obsolete ?');
 end
 if ~isempty(strfind(Cam,'str')) % try display stress
  if isempty(strfind(Cam,'eig'))
   [m,k,mdof]=fe_mknl(model);eval(out1);def.def=ofact(k,def.def);
   feplot(model,def);
  end
  fecom(';colordatastress'); 
 elseif ~isempty(strfind(Cam,'ener')) % try display stress
  if isempty(strfind(Cam,'eig'))
   [m,k,mdof]=fe_mknl(model);def.def=ofact(k,def.def);
   feplot(model,def);
  end
  fecom(';colordataenerk'); 
 end
 out=model;
 if nargout<2; 
 else
     if isempty(def)&&~isempty(out1);eval(out1);end
     out1=def; 
 end% return def
 if RunOpt.Clear
   st='clearvars -global FEel0 FEel1 FEelt FEn0 FEn1 FEnode';eval(st);
   evalin('caller',st);
   return
 end
else % Not a basic elt test

 %if ~isempty(div); femesh(strcat('divide',sprintf(' %i',div)));end
 if isstruct(out); 
  out.Node=FEnode;out.Elt=FEel0;
  out.Elt=feutil('orient;',out);
 [eltid,out.Elt]=feutil('eltidfix;',out);
 end
end % basic test of elements

if RunOpt.Plot
 if isfield(out,'Elt'); mdl=out; 
 else;mdl=femesh; if isempty(FEelt); mdl.Elt=FEel0; end
 end
 feplot(mdl);  
 if sp_util('issdt') 
  cf=feplot;assignin('base','cf',cf);out=cf.mdl;
  disp('Base workspace variable ''cf'' points to the feplot figure')
 end
end

% ------------------------------------------------------------------------
% UnJoin command disconnects nodes that are common to two groups
elseif comstr(Cam,'unjoin')

model=feutil(CAM,femesh('model'),varargin{carg:end});
FEnode=model.Node; FEelt=model.Elt;

% ------------------------------------------------------------------------
else; out='unknown';
end

if nargout==0&&ischar(out)&&comstr(out,'done'); clear out;
elseif ischar(out)&&comstr(out,'unk')
 sdtw('_nb','%s is not a known femesh command',CAM0);clear out
end

if isfield(FE,'Node') % new format with persistent variable
     FE.Node=FEnode;FE.Elt=FEelt;
     if ~isempty(FEel0); FE.El0=FEel0; end
     if ~isempty(FEn0); FE.N0=FEn0; end
     if ~sdtdef('isdeployed'); r1=dbstack;
      if length(r1)>1&& ...
          ~isempty(strfind(r1(1).name,'femesh'))
        assignin('caller','FEnode',FEnode);       
        assignin('caller','FEelt',FEelt);       
        assignin('caller','FEel0',FEel0);       
      end
     end

end
