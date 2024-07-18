function [out,out1,out2]=fe_mknl(varargin)

%FE_MKNL model assembly optimized for non-linear operation
%
% Typically for non-linear operation
%
%   Case=fe_mknl('init',model);                 % performed once
%   mat=fe_mknl('assemble',model,Case,MatType); % performed at every step
%
% For linear matrix assembly. Three calling formats are accepted
%   model=fe_mknl(model)
%   [k,mdof]=fe_mknl(model);
%   [m,k,mdof]=fe_mknl(model);
%
%  The init phase generates Case.GroupInfo which gives
%    GroupInfo = {DofPos Pointers Integ Constit gstate ElMap ...
%       InfoAtNode EltConst}
%    Pointers has one column per element giving %
%     [OutSize1 OutSize2 u3 NdNRule MatDes IntegOffset ConstitOffset ...
%       StateOffset]
%
%  By default, assembling is performed taking into account the boundary
%  conditions (pre/post mult. by Case.T). To obtain the matrices regardless
%  of boundary conditions, use call 'assemble not'. To obtain the matrices
%  without elimination of (accidental) zero coefficients, use call
%  'assemble not nop' (useful for generic symbolic factorization).
%
%  Utilities 
%   NodePos=fe_mknl('NodePos',NNode,model.Elt,cEGI,ElemF);
%   Case.GroupInfo{jGroup,7}=fe_mknl('orientmap',model,cEGI,RunOpt,eltid,InfoAtNode);
%   fe_mknl('gstate -struct',model,Case)
%   nd=feval(fe_mknl('@getPosFromNd'),[],DOF); DofPos=feval(nd.getPosFcn,nd,DOF);

%	E. Balmes, J. Leclere, C. Delforge
%       Copyright (c) 2001-2021 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       Use fe_mknl('cvs') for revision information

%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*NOSEM>

if nargin==0; help fe_mknl
elseif isstruct(varargin{1}) % Assemble the mass and stiffness

 model=varargin{1}; carg=2;
 if isa(model,'v_handle'); model=model.GetData;end
 if carg<=nargin; Cam=varargin{carg};carg=carg+1;else; Cam='';end
 CAM=Cam; [CAM,Cam,RunOpt.Gstate]=comstr('-gstate',[-25 3],CAM,Cam);
 if ischar(Cam); st=horzcat('assemble',Cam);else; st='assemble';end
 if ~isempty(strfind(lower(st),'not')); RunOpt.NoT=1;else; RunOpt.NoT=0;end
 if ~isempty(strfind(Cam,'keep')); st1='initKeep';else; st1='init';end
 if RunOpt.NoT; st1=[st1 'NoT'];end
 if RunOpt.Gstate; st1=[st1 '-gstate'];end
 [Case,model.DOF]=fe_mknl(st1,model);
 m=fe_mknl(st,model,Case,2);
 k=fe_mknl(st,model,Case,1);

 model.K={m,k};  model.Klab={'m','k'};if ~isfield(model,'Opt');model.Opt=1;end
 model.Opt(2,1:2)=[2 1]; model.Opt(2,3:size(model.Opt,2))=0;
 if nargout<2; out=model;
 elseif nargout==2; out=k; 
   if RunOpt.NoT;out1=model.DOF; else;out1=Case.DOF;end
 elseif nargout==3; out=m; out1=k;  
   if RunOpt.NoT;out2=model.DOF; else;out2=Case.DOF;end
 end
 return; 

else
 [CAM,Cam]=comstr(varargin{1},1);carg=2;model=[];
 if carg<=nargin; model=varargin{carg};carg=carg+1;
 elseif comstr(Cam,'@'); out=eval(CAM);return;
 elseif comstr(Cam,'cvs')
  out='$Revision: 1.249 $  $Date: 2024/07/11 17:08:37 $';
  return;
 end
 if isa(model,'v_handle'); model=model.GetData;end
end

% -----------------------------------------------------------------------
%% #InitOptimProduct : init and optimize numbering for product evaluations
if comstr(Cam,'initoptimproduct'); [CAM,Cam]=comstr(CAM,9);

 [Case,model.DOF]=fe_mknl('init',model);
 k=fe_mknl('assemble',model,Case,1);
 [i1,i2,st]=fe_time('optimproduct',k);
 model=stack_set(model,'info','FastestProduct',{i2,st});
 Case.DOF=Case.DOF(i1);Case.T=Case.T(:,i1);
 r1=stack_get(model,'case');
 mode=stack_set(model,'case',r1{2},Case);
 %[Case,model.DOF]=fe_mknl('init',model,Case);
 out=model; out1=Case;

% -----------------------------------------------------------------------
%% #Init Initialize all the info needed for assembly 
elseif comstr(Cam,'init'); [CAM,Cam]=comstr(CAM,5);


RunOpt=struct('LastError','','Gstate',0,'InitFailed',0, ...
              'NeedSeCon',0,'NodePos',[]);

% Build cell array with 
% DofPos Pointers IntegOpt Constit State ElMap
[EGroup,nGroup]=getegroup(model.Elt);

[CAM,Cam,RunOpt.Gstate]=comstr('-gstate',[-25 3],CAM,Cam);
[CAM,Cam,RunOpt.NodePos]=comstr('-nodepos',[-25 3],CAM,Cam);
%% #Init.Phase_1  : generation of the full list of DOFS - - - - - - - - - -2

if ~isempty(strfind(Cam,'keep'))
 i1=strfind(Cam,'keep');CAM(i1+[0:3])='';[CAM,Cam]=comstr(CAM,1);
 if isempty(model.DOF); error('Must have DOFs for initkeep');end
 [mdof,r1,nw,model.Elt,eltid]=feutil('getdof',model); 
else
 [model.DOF,r1,nw,model.Elt,eltid]=feutil('getdof',model);
 mdof=model.DOF;  i1=find(mdof>0);
end
if of_mk('mwIndex')==4&&max(model.DOF)/2^31*100>1
   error('max NodeId too high for 32 bit Matlab');
end
% Some cleaning of where the basis gets resolved needs to be done
if carg>nargin; [Case,NNode]=fe_case(model,['gett' CAM]);
else; Case=varargin{carg};carg=carg+1;
  if size(Case.T,1)~=length(model.DOF)||size(Case.T,2)~=length(Case.DOF)
   [Case,NNode]=fe_case(model,['gett' CAM],Case);
  else
   [Case.Node,Case.bas,NNode,Case.cGL]=feutil(['getnodebas -force' CAM],model);
  end
end
node=Case.Node; mdof=model.DOF;  [nd,nde,N]=fe_mknl('nd',model,Case);

%% #Init.Phase_2 get PL,IL and global nodes
pl=fe_mat('getpl',model); il=fe_mat('getil',model);
if ~isempty(strfind(Cam,'nocon')); RunOpt.NoConn=1;
else;RunOpt.NoConn=0;conn=int32(zeros(size(model.Elt,2),size(model.Elt,1)));
end
RunOpt.EltOrient=stack_get(model,'info','EltOrient');
r1=stack_get(model,'MAP','#Group.*');
if ~isempty(r1); RunOpt.EltOrient(end+(1:size(r1,1)),1:3)=r1;end
r1=stack_get(model,'pro');
if ~isempty(r1); RunOpt.EltOrient(end+(1:size(r1,1)),1:3)=r1;end

RunOpt.Warn={}; EC=[];
if ~isempty(RunOpt.EltOrient); 
 [eltid,model.Elt]=feutil('eltidfix',model);
end

Case.GroupInfo=cell(nGroup,8);

for jGroup=1:nGroup  % Loop on groups

 Case.jGroup=jGroup;pointers=zeros(7,1);RunOpt.GroupInit='groupinit';

%% #Init.Phase_3 Build DofPos connectivity table 

  [ElemF,opt]= feutil('getelemf',model.Elt(EGroup(jGroup),:),jGroup);
  EGID=opt(1); cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;
  [eCall,opt]= fe_super('call',ElemF,model,0); 
  NodePos=[];
  i1=[];try; if ~strcmp(ElemF,'SE');i1=feval(ElemF,'dofcall');end;end
  if ischar(i1);  eval(i1); % Catch for in element building of DofPos
  else; % standard methodology with fixed number of DOF per elt
   i1=fe_super('dofp',ElemF,model); % EltConst may be assigned here
   if iscell(i1);i2=i1;
   elseif opt(1)==1 % unique superelement
    i4=zeros(length(i1),1);
    i3=find(i1>0); i4(i3)=feval(nd.getPosFcn,nd,i1(i3))-1; %full(nd(round(i1(i3)*100)-100))-1;
    i3=find(i1<0); 
    if ~isempty(i3)
      i4(i3)=full(nde(round(-i1(i3)*1000-1000)))-1; % zzz nde to be revised
    end
    i2=int32(i4); cEGI=[];
   elseif ~any(i1<0)
    i2=round(rem(i1(:),1)*100-100);
    i2=model.Elt(cEGI,fix(i1))'*100+i2(:,ones(length(cEGI),1));
    i2(i2<1)=99;  % this allows for variable node numbers (zeros)
    i2=int32(feval(nd.getPosFcn,nd,(i2+100)/100)-1);  %i2=int32(full(nd(i2))-1);
   elseif opt(1)==2 % generic superelement case with element DOFs
    i4=zeros(length(i1),length(cEGI));
    i3=find(i1(:)>0); i2=round(rem(i1(i3),1)*100-100);
    i2=model.Elt(cEGI,fix(i1(i3)))'*100+i2(:,ones(length(cEGI),1));
    i2(i2<1)=99;  % this allows for variable node numbers (zeros)
    i4(i3,:)=feval(nd.getPosFcn,nd,(i2+100)/100)-1; %i4(i3,:)=full(nd(i2))-1;
    i3=find(i1(:)<0);  i2=round(-i1(i3)*1000-2000);
    i2=ones(length(i3),1)*(eltid(cEGI)*1000)'+i2(:,ones(length(cEGI),1));
    i4(i3,:)=full(nde(i2))-1; % zzz nde to be revised
    i2=int32(i4);
   end
  end % local implementation or external approach
  if iscell(i2)&&length(i2)==1; i2=i2{1};end
  if iscell(i2);
   Case.GroupInfo{jGroup,1}=i2;
  else;Case.GroupInfo{jGroup,1}=int32(i2); % DofPos
   Case.DofPerElt(jGroup)=size(i2,1);% size of DofPos
   if RunOpt.NoConn % Deal with SE connectivity
   elseif strcmp(ElemF,'SE')&&(length(cEGI)<10||length(i2)>1000); 
       RunOpt.NeedSeCon=1; % Large enough for SE don't use MeshGraph
   elseif ~isempty(i2)&&~isempty(cEGI);conn(1:size(i2,1),cEGI)=i2;
   end
  end
%% #InitPhase_4 Get MatId and ProId and call elem('integinfo',[MatId ProId],pl,il)
% to generate integ and constit

  if EGID>0
  % Init IntegOpt and Constit, one COLUMM per property
  i1=fe_super('prop',ElemF);
  if any(i1(1:2)) 
   i2=find(i1(1:2)<=size(model.Elt,2)&i1(1:2));
   if isempty(i2); i1=[];
   elseif i1(1)==0 % Just Proid
       i1=model.Elt(cEGI,i1(i2)); [i1,i2,i3]=unique(i1,'rows');i1=[zeros(size(i1)) i1];
   else; i1=model.Elt(cEGI,i1(i2)); [i1,i2,i3]=unique(i1,'rows');
   end
  else; i3=1; % xxx what to do when some prop don't exist
  end

  integ=i1'; constit=zeros(1,size(integ,2));  gstate=[];elmap=[];r2=[];SE=[];
  for j1=1:size(integ,2); 
   if any(strcmp(ElemF,{'SE','rigid'}));
    r2(1:5,j1)=[0 0 j1 j1 0]';SE='done';
   else;
   SE=[];RunOpt.LastError='';
   try % [integ,constit,elmap]= ... call property functions 'buildconst'
    [r1,i1,elmap]=feval(ElemF,'integinfo',integ(1:2,j1),pl,il,model,Case);
   catch
    RunOpt.LastError=lasterr;
    try; 
     [EltConst,i1,model]=fe_super('InitSE',ElemF,model);
     if isfield(EltConst,'SE'); 
      r1=[];i1=[];elmap=[];Case.GroupInfo{jGroup,8}=EltConst;SE=EltConst.SE;
     elseif ~exist(ElemF,'var')
      RunOpt.LastError=sprintf( ...
       'Group %i (%s) not a valid superelement',jGroup,ElemF);
      error(RunOpt.LastError);
     else; error(RunOpt.LastError);
     end
    catch; 
     if sp_util('diag'); 
      if isempty(SE)
        [r1,i1,elmap]=feval(ElemF,'integinfo',integ(1:2,j1),pl,il,model,Case);
      else; 
       if any(strcmp(ElemF,{'celas','mass2'}))&&length(integ)==1&&integ==0
        warning('celas or mass2 no integ')
       else;       error(RunOpt.LastError);
       end
      end
     else
      if any(strcmp(ElemF,{'celas','mass2'}))%&&length(integ)==1&&integ==0
       % elements that do not forcely require integ entries
      else
       RunOpt.InitFailed=1;
       sdtw('_nb', ...
        '%s(''IntegInfo'') (build integ,constit) failed (group %i)',ElemF,jGroup);
       sdtw('_nb',RunOpt.LastError); r1=[];i1=[];
      end
      elmap=[];
     end
    end
   end
   if isstruct(elmap) % elmap contains things that need propagation
    st=fieldnames(elmap);RunOpt.GstateFcn='';
    for j2=1:length(st);
     if comstr(st{j2},'RunOpt_'); % Provide a RunOpt.(field)
      RunOpt.(st{j2}(8:end))=elmap.(st{j2});
     % example : eval('Case.GroupInfo{jGroup,8}.gstate=''keyboard''')
     elseif comstr(st{j2},'eval');eval(elmap.(st{j2}));
     elseif comstr(st{j2},'Case_');error('Please report bug'); 
     end
    end
    if isfield(elmap,'elmap');elmap=elmap.elmap;else;elmap=[];end
   elseif ~isempty(elmap)&&size(elmap,1)~=size(Case.GroupInfo{jGroup,1},1)
     error('ElMap and DofPos are inconsistent for group %i',jGroup);
   end % elmap handling
   if ischar(r1); RunOpt.Warn{end+1}=r1;r1=[];
   elseif any(any(r1))||any(strcmp(ElemF,{'celas','mass2'}))||~any(integ(1:2,j1))
   elseif isfield(SE,'Opt')&&any(SE.Opt(1)==[1 2])
   elseif integ(1)~=0
    sdtw('_nb','Null constit for [MatId %i ProId %i]',integ(1:2,j1));
   end
   integ(1:size(i1,1),j1)=i1;
   % allow for element dependent constit at this level
   if size(integ,2)>=1&&size(r1,2)==length(cEGI); 
     constit=r1;r2=zeros(5,max(length(cEGI),1));
     r2(3,:)=1; i3=1:length(i3);
     if ~isempty(cEGI);r2(4,:)=1:length(cEGI);end % Right columns
   % standard case with property dependent constit
   elseif size(r1,2)>1;% allow constit redefine (*distribution)
       constit=r1; r2(3:4,:)=j1;r2(5,:)=0;i3=1:length(i3);r2=repmat(r2,1,length(i3));
   else; constit(1:size(r1,1),j1)=r1; r2(1:5,j1)=[0 0 j1 j1 0]';
   end
   end
  end % loop on unique ID pairs
  if ~isempty(integ)
   r2(3,:)=(r2(3,:)-1)*size(integ,1);
   r2(4,:)=(r2(4,:)-1)*size(constit,1);
   pointers([1 2 6 7 10],1:length(i3))=r2(:,i3);
   if isfield(Case,'skipLin')&&~isempty(cEGI) 
    %% Acually set the pointers to activation
    in1=ismember(cEGI,Case.skipLin{1});
    if any(in1); pointers(3,in1)=-1; end
   end
  else;pointers=[];
  end
%% #InitPhase_5 let other initializations be performed for various element types
  RunOpt.NLdata=[];
  if ~isempty(RunOpt.EltOrient)
      [Case.GroupInfo{jGroup,7},gstate]=fe_mknl('orientmap',model,cEGI,RunOpt, ...
          eltid,Case.GroupInfo{jGroup,7});
  end
  eCall='';
  if isempty(SE);
  try; % nominal groupinit is a call to a 'constants' command
   eCall=feval(ElemF,RunOpt.GroupInit,integ,constit);
   eval(eCall); % Sets EltConst,pointers,InfoAtNode, sdtweb p_solid('const')
  catch Err; % in the case of superelements something else is needed
   if sp_util('diag')>11; eval(eCall);
   elseif sp_util('diag')
     error('Error %s\n during GroupInit use sdtdef(''diag'',12) for debug', ...
     Err.message);
   elseif ~isempty(eCall); 
    RunOpt.Warn{end+1}=sprintf( ...
    '%s(''GroupInit'') failed (group %i) with\n    %s',ElemF,jGroup,Err.message);
   end
  end
  end % standard groupinit/constants or superelement
  if size(pointers,2)==1&&length(cEGI)>1
      pointers=pointers(:,ones(1,length(cEGI)));
  end
  EC=Case.GroupInfo{jGroup,8};
  if ~RunOpt.NodePos
  elseif ~isempty(NodePos);EC.NodePos=NodePos;
  else;
    cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;
    EC.NodePos=fe_mknl('NodePos',NNode,model.Elt,cEGI,ElemF); 
  end
  if isfield(gstate,'dir') % Gstate initialized in pro.gstate
    NodePos=fe_mknl('NodePos',NNode,model.Elt,cEGI,ElemF); 
    Case.GroupInfo{jGroup,2}=pointers; Case.GroupInfo{jGroup,8}=EC;
    Case=elem0('VectFromDirGState',model,gstate,Case); %elem0('vectfromdir')
    gstate=Case.GroupInfo{jGroup,5};
  end
  Case.GroupInfo{jGroup,3}=int32(integ);
  if isfield(EC,'gstate'); % allow for gstate generation here
     gstate=EC.gstate;EC=rmfield(EC,'gstate');
     % this callback is sometimes used to generated InfoAtNode
     % if NoConn is not needed
     if ischar(gstate)&&~RunOpt.NoConn;eval(gstate);end 
  end
  if isfield(EC,'ElMap')
      if ~isequal(elmap,EC.ElMap)&&~isempty(elmap);
          error('Inconsistent ElMap %s, probable property problem', ...
              sprintf('jGroup=%i(%s)',jGroup,ElemF));
      else;Case.GroupInfo{jGroup,6}=EC.ElMap;EC=rmfield(EC,'ElMap');
          
      end
  else;Case.GroupInfo{jGroup,6}=elmap;
  end
  Case.GroupInfo{jGroup,2}=pointers;
  Case.GroupInfo{jGroup,3}=int32(integ);
  Case.GroupInfo{jGroup,4}=constit;
  Case.GroupInfo{jGroup,5}=gstate;
  if ~isempty(EC)&&isempty(fieldnames(EC));EC=[];end
  if ~isempty(RunOpt.NLdata);EC.NLdata=RunOpt.NLdata;end
  Case.GroupInfo{jGroup,8}=EC;
 end % EGID>0
  i1=fe_super('node',ElemF); 

end % jGroup
%% #Init_after_group_loop
if size(Case.GroupInfo,2)<8; Case.GroupInfo{1,8}=[];end
if ~isempty(RunOpt.Warn); 
    RunOpt.InitFailed=1;Case.MatGraph=RunOpt.Warn;
elseif RunOpt.NoConn;
 if RunOpt.InitFailed;Case.InitFailed=1;end
elseif RunOpt.InitFailed==0
 if ~any(any(conn))||length(mdof)<=6 % bug in meshgraph
    if RunOpt.NeedSeCon&&~any(any(conn));Case.MatGraph=spalloc(length(mdof),length(mdof),0);
    else;Case.MatGraph=sparse(ones(length(mdof)));
    end
 else
  if of_mk('mwIndex')==8; fun=@int64; else; fun=@int32;end
  if ~sp_util('issdt')||(max(Case.DofPerElt)<150&&~sdtkey('isdev')) % verify for 
    Case.MatGraph=of_mk('meshGraph', ...
	      size(mdof,1),nGroup,size(conn,1), ...
	      feval(fun,EGroup),feval(fun,diff(EGroup)-1), ...
              feval(fun,Case.DofPerElt),feval(fun,conn));
  else
   conn(:,EGroup(1:end-1))=-1; % In case one wants to change strategy one day
   for jGroup=1:nGroup; 
    cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
    conn(Case.DofPerElt(jGroup)+1:end,cEGI)=-1;
   end
   JJ=repmat(1:size(conn,2),size(conn,1),1);
   B=sparse(double(conn)+2,JJ,double(conn>=0),length(mdof)+1,size(conn,2)); 
   B(1,:)=[];
   Case.MatGraph=B*B'; %figure(1);spy(Case.MatGraph);
  end
 end
 if RunOpt.NeedSeCon;eval('Case=fesuper(''SeConn'',model,Case,RunOpt);');end

end

if nargout==1; Case.mDOF=mdof;end
if RunOpt.Gstate % initialize Group state information
 model.DOF=mdof;[Case,dc]=elem0('initgstate',model,Case);
 out=Case; out1=mdof;out2=dc;
else % Third argument is for use in fe_mk
 if isfield(RunOpt,'GstateFcn'); eval(RunOpt.GstateFcn);end
 out=Case; out1=mdof; out2={nd,nde,nw,eltid,node,Case.bas};
end % Group state information


%% #Assemble  ----------------------------------------------------------------
elseif comstr(Cam,'assemble'); [CAM,Cam]=comstr(CAM,9);

% model,MatDes  or model,Case,MatDes or model,Case,def,Matdes
% or model,Case,MatDes,def
if carg>nargin; Case=fe_case(model,'getcase');
else; Case=varargin{carg};carg=carg+1;
end
if ~isstruct(Case); MatDes=Case;Case=fe_case(model,'getcase'); 
elseif carg>nargin; error('You must specify a desired matrix type');
else; MatDes=varargin{carg};carg=carg+1;
end
if isstruct(MatDes); def=MatDes;MatDes=varargin{carg};carg=carg+1;
elseif ischar(MatDes);
    def=MatDes; eval(iigui(def,'MoveFromCaller'));
    MatDes=varargin{carg};carg=carg+1;
elseif carg<=nargin; def=varargin{carg};carg=carg+1;
elseif isfield(Case,'mDOF');def=struct('def',zeros(size(Case.mDOF)),'DOF',Case.mDOF); 
else; def=struct('def',zeros(size(model.DOF)),'DOF',model.DOF);
end

[EGroup,nGroup]=getegroup(model.Elt);
if ~isfield(Case,'GroupInfo'); Case=fe_mknl('init',model,Case); end
NNode=sparse(Case.Node(:,1),1,1:size(Case.Node,1));
if isfield(def,'def')&&issparse(def.def);def.def=full(def.def);end
jMat=1; % Provision for when multiple outputs will be supported
% energy computations - - - - - - - - - - - - - - - - - - - - - - -
if any(MatDes(jMat)==[250:299]) 
 
 k=zeros(size(model.Elt,1),size(def.def,2));
 OutType='ener'; MatDes(jMat)=MatDes(jMat)-250; 
 RunOpt.Ener=zeros(size(model.Elt,1),size(def.def,2));

% 300-400 modify internal state in Case
elseif MatDes(jMat)>=300&&MatDes(jMat)<400

 OutType='state'; 
 if isfield(def,'def'); FieldAtDof=def.def; else; FieldAtDof=def;end 
 
% RHS right hand side inits - - - - - - - - - - - - - - -
elseif MatDes(jMat)>=100&&MatDes(jMat)<200  

 OutType='rhs';if carg<=nargin; NNode=varargin{carg}; carg=carg+1;end
 o1=zeros(size(model.DOF,1),1);

% Standard matrix inits - - - - - - - - - - - - - - - - - - - - - - - -
else

 if ~isfield(Case,'MatGraph')||isempty(Case.MatGraph)
   error('Case.MatGraph is not filled, the fe_mknl(''init'') phase failed');
 elseif iscell(Case.MatGraph)
  error('Init failed with\n%s',sprintf('  %s\n',Case.MatGraph{:}));
 end
 k=Case.MatGraph; of_mk( 'fillvalue', k, 0 );
 OutType='mat';

end

% these are standard variable names as stated in elem0
if isfield(Case,'pl');pl=Case.pl;else;pl=fe_mat('getpl',model); end
if isfield(Case,'il');il=Case.il;else;il=fe_mat('getil',model); end
elt=model.Elt; node=Case.Node;eltid=[];
t=cputime; st1='';if isfield(Case,'bas'); bas=Case.bas;end
get_nodeE=elem0('@get_nodeE');
%% Loop on groups - - - - - - - - - - - - - - - -
for jGroup=1:nGroup  
  %fprintf('jgroup%i ',jGroup);evalin('base','eval(RO.st);');feutilb('tops');
  [ElemF,i1]= getegroup(model.Elt(EGroup(jGroup),:),jGroup);
  EGID=i1(1); Case.jGroup=jGroup;
  inode=fe_super('node',ElemF);

  % See #AssembleCall sdtweb elem0('call')
  DofPos=Case.GroupInfo{jGroup,1};
  if EGID<0||isempty(DofPos); continue;end 
  pointers=Case.GroupInfo{jGroup,2};
  integ=Case.GroupInfo{jGroup,3};
  constit=Case.GroupInfo{jGroup,4};
  gstate=Case.GroupInfo{jGroup,5};
  elmap=int32(Case.GroupInfo{jGroup,6});
  InfoAtNode=Case.GroupInfo{jGroup,7};
  EltConst=Case.GroupInfo{jGroup,8};
  if isfield(EltConst,'eltg');
      EltConst.eltg=model.Elt(EGroup(jGroup)+1:EGroup(jGroup+1)-1,:)';
  end
  switch OutType
  case {'mat','mat_of','mat_og'}
    if isfield(EltConst,'SE') % superelement init
     fHandle=fe_super('call',ElemF,model,MatDes);
     if isempty(eltid); eltid=feutil('eltid',model);end
    else;
     try
      [fHandle,SymFlag]=feval(ElemF,'matcall',integ,constit,model,Case);
     catch
      SymFlag=1; fHandle=feval(ElemF,'call',integ,constit,model,Case);
     end
    end 
    pointers(5,:)=MatDes(jMat);
    if isempty(fHandle)
      EltConst.mdl=stack_get(model,'SE',ElemF,'getdata');
      if isempty(EltConst.mdl); fHandle=feval(ElemF,'call');
      else; fHandle='mat_se';EltConst.material='SE';
      end 
    end
    if ischar(fHandle)&&strcmp(fHandle,'mat_of')  % standard of_mk call
       OutType='mat_of';SymFlag=1;
       if length(constit)<2; error('Initialization failed');end
    elseif ischar(fHandle)&&strcmp(fHandle,'mat_og')  % of_mk matrix integration
       OutType='mat_og';SymFlag=0;
       if length(constit)<2; error('Initialization failed');end
    elseif ischar(fHandle)&&strcmp(fHandle,'mat_se')  % of_mk matrix integration
       OutType='mat';SymFlag=0;
    else  % m file elements returning full element matrix
      fullmap=[];
      %fullmap=find(ones(Case.DofPerElt(jGroup),Case.DofPerElt(jGroup)));
      %fullmap=int32(reshape(fullmap, ...
      %   Case.DofPerElt(jGroup),Case.DofPerElt(jGroup)));
      OutType='mat';
    end
    % This should now be handled by pointer(3,:) set to -1;
    if ismember(MatDes,[1 5])&&isfield(EltConst,'NLdata')&& ...
          isfield(EltConst.NLdata,'keepLin')&&EltConst.NLdata.keepLin==0
         continue; % skip group (may also be done by nl_spring assembleInit)
    end
    opt=MatDes; % obsolete kept for compatibility with old elements
    if ~isequal(model.DOF,Case.DOF)&&~isempty(def.def)&&isequal(def.DOF,Case.DOF)
     def.def=Case.T*def.def; def.DOF=model.DOF;
     if MatDes(1)==5; 
      sdtw('_nb','def resizing is not acceptable for NL rhs computations')
     end
    end

  case 'ener' % #AssembleCallEner
    try
      [fHandle,SymFlag]=feval(ElemF,'matcall',integ,constit,model,Case);
      pointers(5,:)=MatDes(jMat);SEopt=2;
      %if ~strcmp(fHandle,'mat_og'); error('Use call');end
    catch; [fHandle,SEopt]=fe_super('call',ElemF,MatDes);
    end
    if isempty(fHandle) ;cEGI=[];
     if ~any(strcmpi(ElemF,{'rigid'}))
       disp([ElemF ' : element function/superelement not found and ignored']);
     end
    elseif comstr(fHandle,'mat_og')
    elseif comstr(fHandle,ElemF)
      disp(fHandle);cEGI=[];
    elseif SEopt(1,1)==1; cEGI=EGroup(jGroup); % single superelement
    end
    if i1(1)<=0; cEGI=[];end % do not compute energies
    opt=MatDes; % obsolete kept for compatibility with old elements
    if isfield(EltConst,'SE')
        fHandle=fe_super('call',ElemF,model,opt);
        if isempty(eltid); eltid=feutil('eltid',model);end
    end
 
  case 'state' % 300-399 modifications

    try;fHandle=feval(ElemF,'state');
    catch;fHandle=''; 
    end
    pointers(5,:)=int32(MatDes);gstate=Case.GroupInfo{jGroup,5};

  case {'rhs','rhs_of','rhs_og'}

    pointers(5,:)=int32(MatDes(jMat));% desired output
    if isfield(EltConst,'RhsDefinition')||~isempty(constit)&&constit(1)<0
      fHandle=elem0('rhsmat_og'); 
    else; [fHandle,SEopt]=fe_super('rhscall',ElemF,MatDes);
    end
    if strcmp(fHandle,'rhs_of')
       OutType='rhs_of'; pointers(30,:)=0; % space for face definitions 
       if length(constit)<2; 
           error('RHS initialization failed, Group%i',jGroup);
       end
    elseif strcmp(fHandle,'rhs_og'); OutType='rhs_og'; 
    else; OutType='rhs';
    end

  otherwise;  error('Not a supported output');
  end
  if isempty(fHandle); cEGI=[];  % in case the element is not assembled
  else 
    cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;  
    NodePos=fe_mknl('NodePos',NNode,elt,cEGI,ElemF);
  end
  switch OutType
  case 'mat' %% #AssembleM .m file element matrix assembly

     for jElt=1:length(cEGI)
      [nodeE,nodeEt]=get_nodeE(Case.Node,NodePos,jElt, ...
          InfoAtNode,DofPos,EltConst); 
      eval([fHandle ';']);
      if isequal(k1,Inf);break;end
      if isempty(k1);
      elseif isempty(fullmap)
       if size(k1,1)>100&&(~issparse(k1)||nnz(k1)/numel(k1)>.75)
        [II,JJ,KK]=find(k1); % xxx this should rather concern an asmsparse optim
        k=k+sparse(double(DofPos(II,jElt))+1,double(DofPos(JJ,jElt))+1,KK,size(k,1),size(k,2));
       else
        i1=[size(k1,1);0;0;DofPos(:,jElt)]; % NDDL, IsSymVal, IsSymK,DofPos(zero based)
        of_mk('asmsparse',k,int32(i1),k1,[]);
       end
      else
       i1=[Case.DofPerElt(jGroup);0;0;DofPos(:,jElt)];
       of_mk('asmsparse',k,int32(i1),k1,fullmap);
      end
     end
  case 'mat_og' %  #AssembleMat_OG of_mk MatrixIntegration assembly
    pointers=int32(pointers);
    if isfield(EltConst,'VectMap') 
      EltConst.VectMap=int32(EltConst.VectMap);
    end
    % Possibly propagate element information
    if sp_util('diag')>10||~isfield(EltConst,'MatrixIntegrationRule') 
      st1='Trying matrix by matrix mex';
    else;
     try;
      icase=int32([Case.DofPerElt(jGroup);SymFlag;0]);
      of_mk('matrixintegration',DofPos,NodePos,Case.Node, ...
        pointers,integ,constit,gstate, ...
        elmap,InfoAtNode,EltConst,def.def, ...
         k,icase);st1=''; if icase(1)==-999;st1='zero';end 
     catch; st1='Normal failed';
     end
    end
    if isempty(st1)||isequal(st1,'zero')
    elseif sp_util('diag')>11||~isfield(EltConst,'MatrixIntegrationRule') 
      st1='Trying matrix by matrix mex';
    else
     try;
      for jElt=1:length(cEGI)
       ke=of_mk('matrixintegration',jElt,NodePos,Case.Node, ...
       pointers,integ,constit,gstate, ...
       elmap,InfoAtNode,EltConst,def);
       i1=[Case.DofPerElt(jGroup);SymFlag;0;DofPos(:,jElt)];
       of_mk('asmsparse',k,int32(i1),ke,elmap);
      end; st1='';
     catch; st1='mat by mat failed';
     end
    end
    i3=MatDes(jMat);
    if sp_util('diag')>11|| (~isempty(st1)&&i3<=length(EltConst.MatrixIntegrationRule)&& ...
            ~isempty(EltConst.MatrixIntegrationRule{i3}))
      % loop on elements is done in ELEM0 for development purposes
      elem0('matrixintegration',DofPos,NodePos,Case.Node, ...
      pointers,integ,constit,gstate, ...
      elmap,InfoAtNode,EltConst,def.def, ...
      k,int32([Case.DofPerElt(jGroup);SymFlag;0]));
    end

  case 'ener' % #AssembleEnerDo
     pointers=int32(pointers);
     if isequal(fHandle,'mat_og'); 
      RunOpt.cEGI=int32(cEGI);
      EltConst.VectMap=int32(EltConst.VectMap);
      try;
       of_mk('matrixintegration',DofPos,NodePos,Case.Node, ...
       pointers,integ,constit,gstate, ...
       elmap,InfoAtNode,EltConst,def.def, ...
       RunOpt,int32([Case.DofPerElt(jGroup);-1;0;]));
      catch;
       elem0('matrixintegration',DofPos,NodePos,Case.Node, ...
       pointers,integ,constit,gstate, ...
       elmap,InfoAtNode,EltConst,def.def, ...
       RunOpt,int32([Case.DofPerElt(jGroup);-1;0;]));
      end
     else % standard element
     for jElt=1:length(cEGI)
      i1=NodePos(:,jElt); nodeE=Case.Node(i1(i1~=0),[5:7 1]);
      if isequal(fHandle,'mat_of')
       point=pointers(:,jElt);
       k1=of_mk(ElemF,int32(point),integ,constit,nodeE);
       k1=k1(elmap);
      else;eval(fHandle);
      end      
      in1=double(Case.GroupInfo{jGroup,1}(:,jElt))+1;
      if size(k1,1)<size(in1,1);in1(size(k1,1)+1:end)=[];end
      in2=find(in1); in1=in1(in2);r1=zeros(size(in1));
      k1=k1(in2,in2);r1=[];
      for j2=1:size(def.def,2)
         r1=def.def(in1,j2);
         RunOpt.Ener(cEGI(jElt),j2) = real(r1'*k1*r1);
      end
      end % mat_og or other element
      if cputime-t>15 % display every 15 seconds
         st1 = comstr(st1,-7,sprintf('Done %i elements',cEGI(jElt)));
         t=cputime;
      end
     end

  case 'rhs'
     if isempty(fHandle) % ignore these elements
     elseif isfield(def,'def')% strategy with field at DOFs
      for jElt=1:length(cEGI)
       i1=NodePos(:,jElt); nodeE=Case.Node(i1(i1~=0),[5:7 1]);
       defe=def.def(double(DofPos(:,jElt))+1,1); point=pointers(:,jElt);
       eval(fHandle)
       in1=double(Case.GroupInfo{jGroup,1}(:,jElt))+1;
       in2=find(in1);in1=in1(in2);
       o1(in1,1)=o1(in1,1)+be(in2);
      end
     else;
      for jElt=1:length(cEGI)
       i1=NodePos(:,jElt); nodeE=Case.Node(i1(i1~=0),[5:7 1]);
       defe=def(:,i1);  point=pointers(:,jElt);
       eval(fHandle)
       in1=double(Case.GroupInfo{jGroup,1}(:,jElt))+1;
       in2=find(in1);in1=in1(in2);
       o1(in1,1)=o1(in1,1)+be(in2);
      end
     end
  case 'rhs_og' %sdtweb elem0('rhs_og')
      o1=elem0('rhs_og',DofPos,NodePos,Case.Node,pointers,integ, ...
       constit,gstate,elmap,InfoAtNode,EltConst,def,o1);

  case 'state' % Modification of element internal state
     for jElt=1:length(cEGI)
      i1=NodePos(:,jElt);
      nodeE=Case.Node(i1,[5:7 1]);
      FieldAtEltDof=FieldAtDof(double(Case.GroupInfo{jGroup,1}(:,jElt))+1);
      point=pointers(:,jElt);
      eval(fHandle); 
      if isempty(jElt);break;end
     end

  case 'rhs_of'  % Legacy elements for RHS computation

     if isstruct(def)
      for jElt=1:length(cEGI)
       i1=NodePos(:,jElt);  nodeE=Case.Node(i1,[5:7 1]);
       point=pointers(:,jElt);
       be=of_mk(ElemF,int32(point),integ,constit,nodeE,gstate, ...
         full(def.def(DofPos(:,jElt)+1,:)),EltConst,InfoAtNode);

       in1=double(Case.GroupInfo{jGroup,1}(:,jElt))+1;
       in2=find(in1);in1=in1(in2);
       o1(in1,1)=o1(in1,1)+be(in2);
      end
     else
      pointers(1,:)=int32(size(DofPos,1));
      for jElt=1:length(cEGI)
       i1=NodePos(:,jElt);  nodeE=Case.Node(i1,[5:7 1]);
       point=pointers(:,jElt);
       be=of_mk(ElemF,int32(point),integ,constit,nodeE,gstate, ...
         def(:,NodePos(:,jElt)),EltConst,InfoAtNode);

       in1=double(Case.GroupInfo{jGroup,1}(:,jElt))+1;
       in2=find(in1);in1=in1(in2);
       o1(in1,1)=o1(in1,1)+be(in2);
      end
     end

  case 'mat_of' % Legacy elements with fixed integration rule
     if isempty(elmap); error('Legacy elements use elmap');end
     if isunix;warning('OpenFEM:obsolete','Legacy elements will be discontinued');end
     for jElt=1:length(cEGI)
      nodeE=Case.Node(NodePos(:,jElt),[5:7 1]);
      point=pointers(:,jElt);
      ke=of_mk(ElemF,int32(point),integ,constit,nodeE);
      i1=[Case.DofPerElt(jGroup);SymFlag;0;DofPos(:,jElt)];
      of_mk('asmsparse',k,int32(i1),ke,elmap);
     end

  otherwise;  error('Not a supported output');
  end

end % loop on groups
switch OutType
case {'mat','mat_of','mat_og'}
  if size(Case.T,1)==size(k,1)&&isempty(strfind(Cam,'not'));
     out=tkt_femknl(Case.T,k); 
  elseif isempty(strfind(Cam,'nop')); 
   out=k+sparse(0); %spalloc(size(k,1),size(k,2),0);
  else; out=k;
  end
  if nargout>=2; out1=Case;end
  if nargout>=3; out2=def;end
case 'ener';  out=RunOpt.Ener;
case {'rhs','rhs_of','rhs_og'}
   out=o1;
case 'state';
otherwise;error('Not a known OutType output');
end

%% #Nd Build Nd from mdof, nd=fe_mknl('nd',model,Case) ----------------------
% DofPos=feval(nd.getPosFcn,nd,DOF);
% (old: DofPos=full(nd(round(DOF*100)-100)))
elseif comstr(Cam,'nd'); [CAM,Cam]=comstr(CAM,5);

Case=varargin{carg};
if isfield(model,'DOF'); mdof=model.DOF;else; mdof=model;end
N=length(mdof);
if isempty(mdof)||all(mdof<0); nd=[];
else
 i2=fix(max(mdof)); 
 if isfield(Case,'Node');i2=max(i2,max(Case.Node(:,1)));end 
 if all(mdof>0);i1=mdof;else;i1=mdof(mdof>0);end
 nd=getPosFromNd([],i1,i2);
 if max(max(nd.nd))>N; error('Something wrong with MDOF'); end
 % older strategy that generated potential un-necessary memory peaks:
 % i3=rem(i1,1);i3=i3*100;i3=round(i3);i1=fix(i1); nd=sparse(i3,i1,1:length(i1),100,i2);clear i1 i2 i3;
end

if isempty(mdof)||min(mdof)>0; nde=[];
else % element DOFs
   i1=find(mdof<0);
   i2=-round(mdof(i1)*1000);
   i2=[rem(i2,1000) fix(i2/1e3)];
   i3=max(i2(:,2)+1);if i3>fix(2^31/1000); error('EltId too large');end
   nde=sparse(i2(:,1),i2(:,2),i1,1000,max(i2(:,2)+1));
end
out=nd; out1=nde; out2=N;
    
% -----------------------------------------------------------------------
%% #NodePos =fe_mknl('NodePos',NNode_or_Node,elt,cEGI,ElemF)
elseif comstr(Cam,'nodepos'); [CAM,Cam]=comstr(CAM,5);

    NNode=model; 
    elt=varargin{carg};carg=carg+1; 
    cEGI=varargin{carg};carg=carg+1; 
    if carg<=nargin;ElemF=varargin{carg};carg=carg+1; 
    else; ElemF=getegroup(elt(cEGI(1)-1,:),0);
    end
    if ~issparse(NNode)&&size(NNode,2)==7;
        NNode=sparse(NNode(:,1),1,1:size(NNode,1));
    end
    inode=fe_super('node',ElemF);
    NodePos=elt(cEGI,inode)'; 
    if min(min(NodePos))==0; 
     i1=find(NodePos);NodePos(i1)=NNode(NodePos(i1));
    else;NodePos=reshape(full(NNode(NodePos)),length(inode),length(cEGI));
    end
    out=int32(NodePos);

% -----------------------------------------------------------------------
%% #OrientMap : deals with material orientation
% [Case.GroupInfo{jGroup,7},gstate]= 
%    fe_mknl('orientmap',model,cEGI,RunOpt,eltid,InfoAtNode);
elseif comstr(Cam,'orientmap'); [CAM,Cam]=comstr(CAM,5);
     
      
    cEGI=varargin{carg};carg=carg+1;
    RunOpt=varargin{carg};carg=carg+1;
    eltid=varargin{carg};carg=carg+1;
    InfoAtNode=varargin{carg};carg=carg+1;
    %[i1,i2]=ismember(eltid(cEGI),RunOpt.EltOrient.EltId);
    
    r1=RunOpt.EltOrient;out=[];out1=[];
    if any(strcmpi(r1(:,1),'MAP')) % if a MAP is given
      jGroup=evalin('caller','jGroup');
      i1=strcmpi(RunOpt.EltOrient(:,2),sprintf('Group_%i',jGroup));
      if any(i1); out=RunOpt.EltOrient{i1,3};end
    end
    %if the info,EltOrient if given
    if strcmp(r1{1,1},'info')&&strcmp(r1{1,2},'EltOrient'); r1=r1{1,3};
    elseif isempty(out)&& any(strcmpi(r1(:,1),'pro'))
       integ=evalin('caller','integ'); 
       for j1=find(strcmp(r1(:,1),'pro'))'
         if isfield(r1{j1,3},'MAP')&&length(integ)>1&&r1{j1,3}.il(1)==integ(2)
           out=r1{j1,3}.MAP; % MAP given in pro.MAP with proper ID
           if isfield(out,'bas');r1=out;break; end
         end
         if 1==2 % isfield(r1{j1,3},'NLdata')&&length(integ)>1&& ...
            %     r1{j1,3}.il(1)==integ(2)
             %% Propagate NLdata to EC, actually handle by nl_spring
             assignin('caller','NLdata',r1{j1,3}.NLdata);
             evalin('caller','RunOpt.NLdata=NLdata;');
         end
         if isfield(r1{j1,3},'gstate')&&length(integ)>1&& ...
                 r1{j1,3}.il(1)==integ(2)
             out1=r1{j1,3}.gstate;break;
         end
         
       end
    end
    
    if isfield(r1,'EltId') 
      % early return if no EltId field
      [i1,i2]=ismember(eltid(cEGI),r1.EltId);if ~any(i1); return;end
      i3=(i2==0);
      if ~any(i3) % All given
      elseif ~any(r1.EltId==0) % Append global
        r1.bas(end+1,1:12)=[0 1 0   0 0 0  1 0 0  0 1 0];
        i2(i3)=size(r1.bas,1);
      else; % Use default associated with 0 EltId
        i2(i3)=r1.bas(find(r1.EltId==0,1,'first'),1);
      end
    elseif isfield(out,'NodePos') % If InfoAtNode pre-built, check
        pointers=evalin('caller','pointers'); 
        if size(out.NodePos,2)~=size(pointers,2)
         jGroup=evalin('caller','jGroup');
         error('Mismatch size(NodePos,2)/nElt in MAP, group %i',jGroup);
        end
    elseif isfield(out,'dir');r1=out;
    else; return;
    end
    ElemF=evalin('caller','ElemF'); in1=feval(ElemF,'node');
    if isfield(r1,'dir')% local orientation map
      Case=evalin('caller','Case'); NNode=evalin('caller','NNode'); % uses Case.Node
      elt=evalin('caller','model.Elt');
      jGroup=evalin('caller','jGroup'); pointers=evalin('caller','pointers'); 
      NodePos=evalin('caller','NodePos'); 
      if isempty(NodePos); 
         NodePos=fe_mknl('NodePos',NNode,elt,cEGI,ElemF);
         assignin('caller','NodePos',NodePos)
      end
      EC=struct('NodePos',NodePos); % force use of Case.Node with EC.NodePos
      out=elem0('VectFromDirAtNode',model,r1,EC,Case.Node); 
      
    elseif isfield(r1,'bas')&&size(r1.bas,1)==1; % constant orient 
     out=struct('data',r1.bas(1,7:12)', ...
      'NodePos',int32(repmat(ones(1,length(cEGI)),length(in1),1)), ...
      'lab',{{'v1x','v1y','v1z','v2x','v2y','v2z'}});
    elseif isfield(r1,'bas') 
     %% one orient per element, sdtweb elem0('EltId2Map')
     out=struct('data',r1.bas(i2,7:12)', ...
      'NodePos',int32(repmat(1:length(cEGI),length(in1),1)), ...
      'lab',{{'v1x','v1y','v1z','v2x','v2y','v2z'}});
    elseif isfield(out,'NodePos') % MAP has been defined in PRO entry
    else; error('Unexpected');
    end
    if isempty(InfoAtNode)
    elseif isequal(InfoAtNode.NodePos,out.NodePos)
     r1=InfoAtNode;clear InfoAtNode;
     r1.data(end+(1:size(out.data,1)),:)=out.data;
     r1.lab(end+(1:length(out.lab)))=out.lab;out=r1;
    elseif size(out.data,2)==1
     r1=InfoAtNode;clear InfoAtNode;
     r1.data(end+(1:size(out.data,1)),:)=repmat(out.data,1,size(r1.data,2));
     r1.lab(end+(1:length(out.lab)))=out.lab;out=r1;
    else
          error('InfoAtNode concatenation case not yet supported');
    end
% -----------------------------------------------------------------------
%% #MapMerge Combines multiple InfoAtNodes
elseif comstr(Cam,'mapmerge'); [CAM,Cam]=comstr(CAM,5);

r1=model; % reference
InfoAtNode=varargin{carg};carg=carg+1; % Things to add/modify
if isempty(r1); out=InfoAtNode;
elseif isequal(InfoAtNode.NodePos,r1.NodePos) % Same Columns
     out=InfoAtNode;clear InfoAtNode;
     r1.data(end+(1:size(out.data,1)),:)=out.data;
     r1.lab(end+(1:length(out.lab)))=out.lab;
     out=r1;
elseif size(r1.data,2)==1 % Single column in MAP1
     out=r1; r1=InfoAtNode;clear InfoAtNode;
     out.data=[repmat(out.data,1,size(r1.data,2));r1.data];
     out.lab=[out.lab(:)' r1.lab(:)'];
elseif size(InfoAtNode.data,2)==1 % Single column in MAP2
     out=r1;r1=InfoAtNode; clear InfoAtNode;
     out.data(end+(1:size(r1.data,1)),:)=repmat(r1.data,1,size(out.data,2));
     out.lab(end+(1:length(r1.lab)))=r1.lab;
    
else
          error('InfoAtNode concatenation case not yet supported');
end

% -----------------------------------------------------------------------
%% #Gstate fe_mknl('gstate -struct',model,Case)
elseif comstr(Cam,'gstate'); [CAM,Cam]=comstr(CAM,7);

 Case =varargin{carg};carg=carg+1;
 [EGroup,nGroup]=getegroup(model.Elt);
 eltid=feutil('eltidfix',model);
 ind=Case.jGroup;
 [CAM,Cam,RO.Struct]=comstr('-struct',[-25 3],CAM,Cam);
 matdes=comstr(CAM,[-1 1]);
 for jGroup=ind
 gstate=Case.GroupInfo{jGroup,5}; EC=Case.GroupInfo{jGroup,8};
 Nw=EC.Nw; if length(Nw)>1;Nw=Nw(2);end
 if Nw>size(EC.w,1);error('Mismatch');end
 cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;
 i1=size(gstate); 
 try
  
  if isempty(gstate)||(isfield(gstate,'Y')&& ...
    ~isequal(size(gstate.Y),[length(EC.StrainLabels{matdes}), ...
     Nw,length(cEGI)]))
     gstate=zeros(length(EC.StrainLabels{matdes}),Nw,length(cEGI));
  elseif isfield(gstate,'Y'); gstate=gstate.Y;
  end
  gstate=struct('X',{{EC.StrainLabels{matdes}(:),EC.w(1:Nw,:),eltid(cEGI)}}, ...
        'Xlab',{{sprintf('Stress %i',matdes),'Gauss','EltId'}}, ...
        'Y',gstate);
  if isfield(EC,'StressOut')&&matdes==1 % Allow stress component extraction
   gstate.X=[EC.StressOut.X gstate.X(2:end)];
   gstate.Xlab=[EC.StressOut.Xlab gstate.Xlab(2:end)];
   gstate.Y=zeros(cellfun(@(x)size(x,1),gstate.X));
  end
  Case.GroupInfo{jGroup,5}=gstate;
 end
 end
 out=Case;
%% #BlockStrain : precompute all jacobians
elseif comstr(Cam,'blockstrain'); [CAM,Cam]=comstr(CAM,7);
 
[Case,model.DOF]=fe_mknl('init NoT',fe_case(model,'reset'));
[EGroup,nGroup]=getegroup(model.Elt);
NNode=sparse(Case.Node(:,1),1,1:size(Case.Node,1));
out=cell(1,nGroup);
for jGroup=1:nGroup
 [ElemF,i1,ElemP]= getegroup(model.Elt(EGroup(jGroup),:),jGroup);
 cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
 NodePos=fe_mknl('NodePos',NNode,model.Elt,cEGI,ElemP);
 EC=Case.GroupInfo{jGroup,8};
 EC.J=zeros(9,length(EC.jdet));
 rule=Case.GroupInfo{jGroup,2}(4,1);
 BS=EC.N';BS(:,:,2)=EC.Nr';BS(:,:,3)=EC.Ns';BS(:,:,4)=EC.Nt';
 constit=Case.GroupInfo{jGroup,4};
 dd=double(full(EC.ConstitTopology{1}));dd(dd~=0)=constit(dd(dd~=0));
 BS=struct('c',permute(BS,[1 3 2]),'DofPos',Case.GroupInfo{jGroup,1}, ... 
     'wjdet',zeros(10,EC.Nw,length(cEGI)), ...
     'Lambda',dd,'rho',[]);
 for jElt=1:length(cEGI)
  EC.nodeE=model.Node(NodePos(:,jElt),[5:7 1]);
  of_mk('buildndn',3,EC);
  if jElt==1&&nnz(EC.J(:))==0;error('Zero jacobian is wrong');end
  BS.wjdet(:,:,jElt)=[(EC.jdet.*EC.w(:,4))';EC.J];
  %r2=EC.NDN(:,EC.Nw+1:EC.Nw:end);
  %r3=[EC.Nr(1,:);EC.Ns(1,:);EC.Nt(1,:)]'*reshape(EC.J(:,1),3,3)  
 end
 out{jGroup}=BS;
end
out1=model;out2=Case;


%% #End ----------------------------------------------------------------------
else; error('''%s'' is unknown',CAM);
end

%% #SubFunc ------------------------------------------------------------------
%% #tkt_femknl - - -----------------------------------------------------------
function K=tkt_femknl(T,K)

if sp_util('issdt'); K=feutilb('tkt',T,K);
else
 Tt=T';
 K=Tt*K*T;
end

%% #getPosFromNd - - ---------------------------------------------------------
function out=getPosFromNd(nd,DOF,varargin)
% strategy to find position of DOF in a global mdof vector, using a pseudo object
% DofPos=feval(nd.getPosFcn,nd,DOF);
% DofPos=feval(nd.getPosFcn,NodeId,[ind|dofext]);
% DofPos=feval(nd.getPosFcn,nd,DOF,'offset',81); % eg: fsc dof .19
% slave=feval(nd.setDofListNd,nd); % to generate a partial dofpos list
if isempty(nd); % constructor
 nd=struct('NNode',[],'nd',[],'getPosFcn',@getPosFromNd,...
  'setDofListNd',@setDofListNd);%,'size',[]); 
end
if isempty(nd.NNode) % init: build and output nd
 % feval(fe_mknl('@getPosFromNd'),[],DOF,max_node)
 node=unique(fix(DOF)); out=nd;
 out.NNode=sparse(node,1,1:length(node));
 if ~isempty(varargin)&&varargin{1}>length(out.NNode);
     out.NNode(varargin{1},1)=0;
 end
 out.nd=sparse(round(rem(DOF,1)*100),out.NNode(fix(DOF))+1,1:length(DOF),...
  100,max(out.NNode)+1);
 
else % output index of DOF in mdof using nd
 if ~isempty(varargin); i1=varargin{1}; carg=2; else; i1=''; carg=1; end
 try; in1=nd.NNode(fix(DOF))+1; % columns in nd.nd based on the nodeIndex from NodeId
 catch err 
  % robustness to cases (mostly fesuper) when nd was built with an
  % incomplete DOF vector. Strategy is to rethrow dofpos=0 for these DOF
  % this is realized by locally increasing the size of nd.nd and nd.NNode
  i2=fix(DOF)==0; 
  if any(i2); DOF(i2)=max(fix(DOF))+1+DOF(i2); end
  if length(nd.NNode)<max(max(fix(DOF))) % some DOF were not present at build
   if sp_util('diag')>0; warning('getPosFromNd, wrong DOF input'); end
   [node,unu]=find(nd.NNode); node=unique([node;fix(DOF(:))]);
   nd.NNode=sparse(node,1,1:length(node)); nd.nd(1,length(node)+1)=0;
   in1=nd.NNode(fix(DOF))+1;
  else; % should not happen but verbose fail if coding problem
   error('getPosFromNd failed with:\n%s',err.getReport);
  end
 end
 if isnumeric(i1) % get dofPos based on fieldDof and nodeId
  if any(i1<1); i1=round(i1*100); end % accept dofExt instead of fieldDof
  out=full(nd.nd(i1,in1)); 
 else % get dofPos based on dof input (usual)
  if ~isempty(i1)
   % accept custom offset to get other fieldDof for same nodes
   % eg: get .19 DOF from mech_dof, offset = 81, see fsc and fe_cyclic
   if comstr(lower(i1),'offset'); i1=varargin{carg}; carg=carg+1;
   else; error('getPosFromNd: option %s not recognized',i1);
   end
  else; i1=100;
  end
  % robust output to the DOF vector shape
  if any(size(DOF)==1); r1=size(DOF); DOF=DOF(:); else; r1=[]; end
  out=full(nd.nd(round((in1+rem(DOF,1))*100)-i1)); % dofPos
  if ~isempty(r1); out=reshape(out,r1(1),r1(2)); end
 end
end
%% #SetDofListNd - - --------------------------------------------------------3
function nd=setDofListNd(nd)
% tranform nd object into a sub_dofpos list object % fe_mpc slave strategy
 nd.getListedDofPos=@getListedDofPos;
 nd.appendDof=@appendDof;
 nd.ndG=nd.nd;
 nd.nd=sparse(size(nd.ndG,1),size(nd.ndG,2));
%% #getListedDofPos - - -----------------------------------------------------3
function out=getListedDofPos(nd)
 [unu1,unu2,out]=find(nd.nd);
%% #appendDof - - -----------------------------------------------------------3
function nd=appendDof(nd,DOF) % sdtweb fe_mpc('addSlave')
 % get DofPos of dof, and add it to nd list
 i1=rem(DOF,1); i2=nd.NNode(fix(DOF))+1;
 dofPos=full(nd.ndG(round((i2+i1)*100)-100));
 nd.nd=nd.nd+sparse(round(i1*100),i2,dofPos,100,size(nd.nd,2));
