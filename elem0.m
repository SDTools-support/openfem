function [out,out1,out2]=elem0(CAM,varargin);

% Default element function used for elements that conform strictly to an
% OpenFEM standard
%
% Documented commands are
%  VectFromDir
%  elem0('GaussObserve',rule,integ,constit,model,Case,cEGI)
%  TensorT : tensor transformations
%
% See <a href="matlab: sdtweb _taglist elem0">TagList</a>
% See sdtweb elem0


%       Etienne Balmes
%       Copyright (c) 2001-2025 by INRIA and SDTools, All rights reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license


if ischar(CAM); [CAM,Cam]=comstr(CAM,1);
elseif isa(CAM,'int32'); Cam='callback';
end
%#ok<*ASGLU,*CTCH,*NASGU,*NOSEM>

if comstr(Cam,'callmat_og');
%% #Generic generic callbacks for element functions
%% #callmat_og --------------------------------------------------------------2

 out='[k1,m1]=elem0(''mat_og'',jElt,NodePos,Case.Node,pointers,integ,constit,gstate,elmap,InfoAtNode,EltConst,def);';

% Elements of the Mat_Og family use a [1.01 ... N.01 1.02 ... N.02 ...]
% internal DOF sorting because this is appropriate for BLAS calls
% The standard numbering returned by elements is [1.01 1.02 ... 2.01 2.02 ...]
% elem0('elmapmat_og',[Nnode NdofPerNode]);
elseif comstr(Cam,'elmapmat_og');

  ind=varargin{1};
  if length(ind)>2 % Provides the VectMap 
     out=reshape(1:length(ind)^2,length(ind),length(ind));
     out(ind,ind)=out;
  else
    i2=prod(ind); i1=reshape(1:i2,ind(2),ind(1))';
    out=reshape(1:i2^2,i2,i2); out(i1,i1)=out;
  end

elseif comstr(Cam,'integinfo');[CAM,Cam]=comstr(CAM,10);
%% #IntegInfo ---------------------------------------------------------------2
% standard integinfo call for og elements

   if ~isempty(Cam); opt=comstr(CAM,-1);end
   [out,out1,out2]= ...
     p_solid('buildconstit',[varargin{1};opt(:)],varargin{2:end});
 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'groupinit');[CAM,Cam]=comstr(CAM,10);
%% #GroupInit standard call for og elements ---------------------------------2

 if comstr(Cam,'og'); [CAM,Cam]=comstr(CAM,3);end
 if isempty(Cam)&&ischar(varargin{1})
  out=sprintf('[Case.GroupInfo{jGroup,8},pointers,Case.GroupInfo{jGroup,7}]=%s(''Constants'',pointers,integ,constit,model,Case,cEGI,RunOpt);',varargin{1});
 elseif comstr(Cam,'shell')
 out=sprintf('[Case.GroupInfo{jGroup,8},pointers,Case.GroupInfo{jGroup,7}]=elem0(''Constants%s p_shell'',pointers,integ,constit,model,Case,cEGI,RunOpt);',varargin{1});
 else
  try; ElemP=evalin('caller','ElemP'); 
    out=sprintf('[Case.GroupInfo{jGroup,8},pointers,Case.GroupInfo{jGroup,7}]=elem0(''Constants%s'',pointers,integ,constit,model,Case,cEGI,RunOpt);',ElemP);
  catch;error('GroupInit failed'); 
  end
 end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'constants');[CAM,Cam]=comstr(CAM,10);
%% #Constants standard Constants call for og elements - -2
% typically shell or property type switch on shell or 2D

integ=varargin{2};constit=varargin{3}; out='empty'; out1=[]; 
RunOpt.Type='';
if constit(1)<0|| ...% surface elements - - - - - - - - - - - - - - - - -
   ( length(integ)==5&&integ(4)&&integ(3)/integ(4)==3 )
 st={'quad4','q4p';'tria3','t3p';'tria6','t6p';'quad9','quad9'; ...
     'quadb','q8p'};
 Cam=st{strcmp(sscanf(Cam,'%s',1),st(:,1)),2};
 [out,out1,out2]=feval(Cam,'constants',varargin{:});

elseif comstr(Cam,'shell'); RunOpt.Type=@p_shell;[CAM,Cam]=comstr(CAM,6);
elseif length(constit)>2&&any(remi(constit(3),[],1)==[0 5])&& ...
        sp_util('issdt')&&integ(3)>=6*integ(4) % Allow more than 6 for piezo
% shell map bypass if type is 0 or 5
   eval('[out,out1,out2]=q4cs([''constants'' CAM],varargin{:});')    
   return
elseif comstr(Cam,'quad4') % constant inits for quad4 - - - - - - -
 RunOpt.Type=@p_shell;
 if length(constit)>2&&constit(3)==4; % MITC4 call
 elseif ~any(integ(3)==[5 6]*double(integ(4)));RunOpt.Type=@p_solid;
 elseif constit(1)==-1&&~strcmp(fe_mat('typep',constit(2)),'p_shell') % redirect
   RunOpt.Type=fe_mat('typep',constit(2));
 else
   [out,i2]=p_shell('const',integrules('q4q',109),varargin{2:end});
   out1=varargin{1};out1(4,:)=i2; % Tell MatrixIntegration this is 3dsurf 
   out2=[]; % normals
 end
elseif comstr(Cam,'tria3') % constant inits for tria3 - - - - - - -
 RunOpt.Type=@p_shell;% 2D shell
elseif comstr(Cam,'tria6') % constant inits for tria6 - - - - - - -
 RunOpt.Type=@p_shell;% 2D shell
elseif comstr(Cam,'quad9') % constant inits for quad9 shell/2d - - - - - - -
 RunOpt.Type=@p_shell;% 2D shell
elseif comstr(Cam,'quadb') % constant inits for quadb - - - - - - -
 RunOpt.Type=@p_shell;% 2D shell
else; i1=find(Cam==' '); 
    if ~isempty(i1);RunOpt.Type=eval(sprintf('@%s',Cam(i1+1:end)));
    else;RunOpt.Type=@p_solid;
    end
end
if ~isequal(out,'empty')
else
  if ~isequal(RunOpt.Type,@p_shell)
  elseif ~any(integ(3)==[5 6]*double(integ(4)));RunOpt.Type=@p_solid;
  elseif constit(1)==-1&&~strcmp(fe_mat('typep',constit(2)),'p_shell') % redirect
    RunOpt.Type=fe_mat('typep',constit(2));
  end
  [out,i2,out2]=feval(RunOpt.Type,'const',CAM,varargin{2:end});
  %[out,i2]=feval(RunOpt.Type,'const',CAM,varargin{2:end});out2=[];
  if i2==2 % possibly force 23D integ
      Case=varargin{5};
      if isfield(Case,'Node')&&any(Case.Node(:,7)); i2=23;end
  end
  out1=varargin{1};out1(4,:)=i2; % Tell MatrixIntegration this is 3dsurf 
end


elseif comstr(Cam,'dofcall');
%% #dofcall -----------------------------------------------------------------2
% standard DOF building for variable field elements
 out='[i2,i3]=p_solid(''BuildDof'',model,cEGI,nd,ElemF);';
elseif comstr(Cam,'matcall');
%% #Matcall that allows selection of mat_og and mat_of strategies - - -2

 if ischar(varargin{1})
  out=sprintf(['[k1,m1]=%s(nodeE,elt(cEGI(jElt),:),', ...
     'pointers(:,jElt),integ,constit,elmap);'],varargin{1});out1=0;
 elseif varargin{1}(end,1)==-9999 % Elements issued from Modulef in of_mk_sub.c
  out='mat_of'; out1=1; % call SymFlag
 else % Newer elements handled by of_mk.c MatrixIntegration
   out='mat_og'; out1=0; % mat_og and non symmetric
 end
 
elseif comstr(Cam,'callback');
%% #callback ----------------------------------------------------------------2
% This is used for testing purposes in of_mk('matrixint') development
   %of_mk('matrixintegration',DofPos,NodePos,Case.Node, ...
   %     pointers,integ,constit,gstate, ...
   %     elmap,InfoAtNode,EltConst,def.def, ...
 [integ,constit,gstate,elmap,InfoAtNode,EltConst,def,jElt,DofPos]=deal(varargin{:});
 point=CAM(:,jElt(1)+1);
 nodeE=EltConst.nodeE;Nw=EltConst.Nw;defe=EltConst.defe;
 Nnode=EltConst.Nnode;ke=EltConst.ke; Ndof=size(ke,1);
 % M-file implementation of general non linear problem
 if isfield(EltConst,'FieldDofs')
   EC=EltConst;
   dd=eval(EC.ConstitTopology{point(5)});
   sp_util('setinput',EC.constit,dd,zeros(1));   
   sp_util('setinput',jElt,-1,zeros(1));
   out=EC;
 else % M File implementation of elastic with temperature sensitivity
 ci_ts_egt=[1 5 9 8 7 4]';
 ind_ts_eg=[1 6 5 6 2 4 5 4 3]; % e_11 22 33 23 31 12
 ind_ts_egt=[1 6 5 6 2 4 5 4 3]'; % e_11 22 33 23 31 12
 jElt=jElt(1)+1;
 if size(nodeE,2)>4&&(size(constit,1)<10||constit(39+double(point(7))));
     at=zeros(3);
 else;at=[];
 end
 
 for jW=0:Nw-1; % HERE IS THE LOOP ON INTEGRATION POINTS 

   r1=EltConst.NDN(:,[Nw 2*Nw 3*Nw]+jW+1);
   U_ij=reshape(defe(:,1),3,Nnode)*r1; F_ij=U_ij+eye(3);

  i1=EltConst.ConstitTopology{1};%i1=i1(ind_ts_eg,ind_ts_eg);
  if isfield(EltConst,'CTable');
      z=ConstitInterp(nodeE,EltConst.NDN(:,jW+1),constit,EltConst.CTable);
      of_time('cinterp',nodeE,EltConst.NDN,constit,EltConst.CTable,double(jW));
      if any(abs(z-constit)>eps*constit);error('Mismatch');end
  end
  if size(constit,1)<10
     % table for rho eta E nu G alpha
     r3=constit((1:9)+double(point(6+1)));
     nu=r3(4); E=r3(3); 
     G=r3(6); if G==0; G=E/2/(1+nu);end
     r1=nu./(1-nu); % n/(1-n) 
     r2=E.*(1-nu)./(1+nu)./(1-2*nu); % E(1-n)/(1+n)(1-2*n)
     d2wde2=zeros(6,6);
     d2wde2([1 2 3 7 8 9 13 14 15])= ...
      r2*[1 r1 r1  r1 1 r1  r1 r1 1]';
    d2wde2([22 29 36])=G; % G
    if size(r3,1)>8;at([1 5 9])=r3(8);T0=r3(9);end
  else;d2wde2=constit(double(i1)+double(point(6+1)));
    if ~isempty(at);at=reshape(constit((39:47)+double(point(7))),3,3);
        T0=constit(48+double(point(7)));
    end
  end
  Sigma=d2wde2(ind_ts_eg,ind_ts_eg)*U_ij(:);
  if  ~isempty(at); % add thermal stress      
    r1=EltConst.NDN(:,jW+1)'*nodeE(:,5)-T0; % dT
    r1=at*r1;r1=d2wde2*r1(ci_ts_egt);   % Stress at integration point
    Sigma=Sigma-r1(ind_ts_eg);
  end
  if ~isempty(gstate)
   Sigma=Sigma+gstate(jW*6+ind_ts_eg,jElt);
  end
  % assignin('base','d2wde2_elas',d2wde2);
  %DD(ci_ts_eg,ci_ts_eg)
  Mecha3DInteg(ke,F_ij,d2wde2,reshape(Sigma,3,3), ...
   def,EltConst.w(:,4),EltConst.jdet,EltConst.NDN,Nnode,Ndof,Nw,jW)

 end
 end % Callback implementation
 
%% #call --------------------------------------------------------------------2
elseif comstr(Cam,'call');[CAM,Cam]=comstr(CAM,5);

 out1=0; out=''; % Call, SymFlag
 if nargin==1; ID=0;else; ID=varargin{1};end
 if nargin>2&&any(varargin{2}(1,1)==[-1 -3])|| ... % generic topology holder
         length(ID)>3&&~any(ID(3)==[5 6]*double(ID(4)))
   [out,out1]=elem0('matcall',varargin{:});
 else
 if nargin<3; i1=0;else;i1=remi(varargin{2}(3,1),[],1);end% Formulation
 switch Cam % topology specific redirection
 case 'quad4'  % #CallQuad4 -3
  if nargin>2&&size(varargin{2},1)==9        
   out=mitc4('call');
  elseif nargin>2&&any(i1==[0 5])&&varargin{2}(1)>=0&&sp_util('issdt');
   eval('out=q4cs(''call'');');
  end
 case 'quadb'  % #CallQuadb -3
  if nargin>2&&any(i1==[0 5])&&varargin{2}(1)>=0&&sp_util('issdt');
   eval('out=q4cs(''call'');');
  end
 case 'tria3' % #CallTria3 -3
  if nargin>2&&any(i1==5)&&varargin{2}(1)>=0&&sp_util('issdt');
   eval('out=q4cs(''call'');');
  end
 case 'tria6'
  if nargin>2&&any(i1==[0 5])&&varargin{2}(1)>=0&&sp_util('issdt');
   eval('out=q4cs(''call'');');
  end
 end % topology specific redirection
 end
 if isempty(out) % default for fe_mk (obsolete)
   out=sprintf(['[k1,m1]=%s(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),', ...
     'integ,constit,elmap,InfoAtNode,EltConst,def);'],Cam);
 end
%% #rhsmat_og ------------------------------------------------------------------
elseif comstr(Cam,'rhsmat_og'); out='rhs_og';

%% #rhs_og ---------------------------------------------------------------------
% Matlab implementation of single element matrix integration call. 
% This is only meant for development purposes
elseif comstr(Cam,'rhs_og');

 [DofPos,NodePos,node, ...
      pointers,integ,constit,gstate, ...
      elmap,InfoAtNode,EltConst,def,out]=deal(varargin{:});

Nw=EltConst.Nw;Nnode=EltConst.Nnode; Ndof=size(DofPos,1);
%if any(pointers(4,1)==[3 31 32 23]); Ndim=3; else; Ndim=2;end
Ndim=size(DofPos,1)/EltConst.Nnode;
if ~isfield(EltConst,'VectMap')||numel(EltConst.VectMap)~=size(DofPos,1)
  i1=reshape(1:Ndof,Ndof/Nnode,Nnode)'; EltConst.VectMap=i1(:);
end
i1=[]; 
jdef=1;RunOpt.ndef=1;
try
% Volume load with field at DOF
if isfield(def,'def')&&pointers(5,1)==100&&isempty(InfoAtNode)&& ...
        isfield(EltConst,'NodePos')
 i2=unique(round(rem(def.DOF(double(DofPos(:,1))+1),1)*100))/100;
 if length(i2)*EltConst.Nnode==size(DofPos,1)
  RunOpt.ndef=size(def.def,2);out(1,RunOpt.ndef)=0;
  InfoAtNode=struct('NodePos', ...
    reshape(int32(1:numel(EltConst.NodePos)),size(EltConst.NodePos)), ...
    'data',[],'lab',{{'x','y','z'}});
  Case=def.Case;jGroup=1;
  if size(Case.GroupInfo,1)>1&&Case.jGroup>1
    Case.GroupInfo=Case.GroupInfo(Case.jGroup,:);
    Case.DofPerElt=Case.DofPerElt(Case.jGroup);
  end
  EltConst.material='';
  EltConst.RhsDefinition=int32( ...
   [101 0                1 0     0 0 -1    0 Nw;
    101 EltConst.Nnode   1 0     1 0 -1    0 Nw;
    101 EltConst.Nnode*2 1 0     2 0 -1    0 Nw]);
  while jdef
   InfoAtNode.data=reshape(full(def.def(double(DofPos)+1,jdef)),3,[]);
   b1=zeros(size(def.def,1),2);
   of_mk('matrixintegration',DofPos,NodePos,Case.Node, ...
        int32(pointers),integ,constit,gstate, ...
        elmap,InfoAtNode,EltConst,b1, ...
         [],int32([Case.DofPerElt(jGroup);0;0]));st1='';
   out(:,jdef)=out(:,jdef)+b1(:,2);
   if jdef==RunOpt.ndef;jdef=0;else;jdef=jdef+1;end
  end
 end
 return;%z=out;out=zeros(size(out));jdef=1;'z'
end
end % try compiled call for volume at DOF

while jdef
for jElt=1:size(DofPos,2);

 point=pointers(:,jElt);
 EltConst.nodeE=node(NodePos(:,jElt),5:7);
 of_mk('buildndn',point(4),EltConst);
 inde=double(DofPos(:,jElt))+1; 
 inde=inde(EltConst.VectMap);in2=find(inde>0);
if isfield(def,'dir')&&size(def.dir,2)==size(DofPos,2)
 defe=def.dir(:,jElt*ones(1,Nnode));
elseif isfield(def,'AtNode') % field at nodes
  defe=def.AtNode(double(EltConst.NodePos(:,jElt)),:)';
elseif isfield(def,'def') % field at DOFs
  defe=double(DofPos(:,jElt))+1;defe(defe~=0)=def.def(defe(defe~=0),jdef);
  defe=reshape(defe,Ndim,EltConst.Nnode);
  RunOpt.ndef=size(def.def,2);
elseif ~isempty(def) % correct for missing dimensions in input
  defe=def(1:min(Ndim,size(def,1)),NodePos(:,jElt));
  if size(def,1)<Ndim; defe(Ndim,1)=0;end
end
if size(out,2)<RunOpt.ndef;out(1,RunOpt.ndef)=0;end

switch double(point(5))
case 100  % Volume load (not proportional to density)

 for jW=0:Nw-1
  r1=EltConst.NDN(:,ones(Ndim,1)*(jW+1))*diag((defe*EltConst.NDN(:,jW+1))* ...
   EltConst.jdet(jW+1)* EltConst.w(jW+1,4));
  out(inde(in2),jdef)=out(inde(in2),jdef)+r1(in2);
 end % loop on jW
case 101  % Pressure load

 % sign was fixed to usual convention on 15/11/06
 defe=-EltConst.bas(7:9,:)*diag(defe(1,:)*EltConst.NDN(:,1:Nw));
 for jW=0:Nw-1
  r1=EltConst.NDN(:,ones(Ndim,1)*(jW+1))*diag(defe(:,jW+1)* ...
   EltConst.jdet(jW+1)* EltConst.w(jW+1,4));
  out(inde(in2),jdef)=out(inde(in2),jdef)+r1(in2);
 end % loop on jW

case 102 % #rhs_ogDens inertia load (proportional to density) -2
  if ~isempty(i1)
  elseif isfield(EltConst,'DensPos'); i1=EltConst.DensPos;
  else; i1=1; 
  end
 for jW=0:Nw-1
  r1=EltConst.NDN(:,[jW jW jW]+1)*diag((defe*EltConst.NDN(:,jW+1))* ...
   EltConst.jdet(jW+1)* EltConst.w(jW+1,4)*constit(i1+point(7)));
  % it is assumed (NOT CHECKED) that rho(1) contains density
  out(inde(in2),jdef)=out(inde(in2),jdef)+r1(in2);
 end % loop on jW

case 103 % Initial stress load (F = \int \sigma_0 \epsilon dv)

 % row, NDN, DDL, NwStart NwRule
 r2=EltConst.StrainDefinition{1}; Nstrain=length(EltConst.StrainLabels{1});
 Nnode=EltConst.Nnode;
 if isstruct(gstate);gstate=reshape(gstate.Y,[],size(gstate.Y,3));end
 for jW=0:Nw-1
  for ji=0:size(r2,1)-1 % Strain definition rows
   jj=(r2(ji+1,2)-1)*Nw+jW; % NDN column
   jk=(r2(ji+1,3)-1)*Nnode; % DOF block
   ib=double(DofPos(EltConst.VectMap(jk+[1:Nnode]),jElt))+1;
   r1=EltConst.NDN(:,jj+1)*gstate(jW*Nstrain+r2(ji+1,1),jElt);
   out(ib,jdef)=out(ib,jdef)+r1*EltConst.jdet(jW+1)* EltConst.w(jW+1,4);
  end

 end % loop on jW
case 104  % 2D pressure load
 defe=EltConst.bas(3:4,:)*diag(defe(1,:)*EltConst.NDN(:,1:Nw));
 for jW=0:Nw-1
  r1=EltConst.NDN(:,ones(Ndim,1)*(jW+1))*diag(defe(:,jW+1)* ...
   EltConst.jdet(jW+1)* EltConst.w(jW+1,4));
  out(inde(in2),jdef)=out(inde(in2),jdef)+r1(in2);
 end % loop on jW

case 105  % 2D line (surface load)
 for jW=0:Nw-1
  r1=EltConst.NDN(:,ones(Ndim,1)*(jW+1))*diag((defe*EltConst.NDN(:,jW+1))* ...
   EltConst.jdet(jW+1)* EltConst.w(jW+1,4));
  out(inde(in2),jdef)=out(inde(in2),jdef)+r1(in2);
 end % loop on jW

otherwise; sdtw(' point(5)=%i not a supported load type',point(5));break;
end
end % loop on elements of the group
if jdef==RunOpt.ndef;jdef=0;else;jdef=jdef+1;end
end % loop on deformations

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'mat_og');

  varg=varargin; 
  if ~isa(varg{4},'int32'); 
   varg{4}=int32(varg{4});evalin('caller','pointers=int32(pointers);');
  end
  if ~isa(varg{2},'int32'); 
   varg{2}=int32(full(varg{2}));evalin('caller','pointers=int32(pointers);');
  end
  [out,out1]=of_mk('matrixintegration',varg{:});
  if ~isempty(varg{8});
   out=out(varg{8});if ~isempty(out1); out1=out1(varg{8});end
  end

% -----------------------------------------------------------------------
elseif comstr(Cam,'matrixintegration');
%% #MatrixIntegration --------------------------------------------------------
% Matlab implementation of single element matrix integration call. 
% THIS IS ONLY MEANT FOR DEVELOPMENT PURPOSES

[DofPos,NodePos,node, ...
      pointers,integ,constit,gstate, ...
      elmap,InfoAtNode,EltConst,def,k,AsmPoint]=deal(varargin{:});

if isempty(k); 
 jElt=DofPos; DofPos=evalin('caller','DofPos');
elseif isa(k,'struct');
 cEGI=k.cEGI; out1=k.Ener;
end

Nw=EltConst.Nw;Nnode=EltConst.Nnode; 
try;Ndof=Nnode*length(EltConst.DofLabels);catch;Ndof=0;end
NdofPerField=size(EltConst.NDN,1);
if ~isfield(EltConst,'material'); EltConst.material='';end
Ncondense=0;if isfield(EltConst,'Ncondense'); Ncondense=EltConst.Ncondense;  end
if ischar(gstate);gstate=[];end
if isfield(EltConst,'N')
    if ~isfield(EltConst,'jdet'); error('Missing jdet field');end
end

EC=struct('v1x',0); % reproduce some of the C EltConst structure
if isfield(InfoAtNode,'lab')&&any(strcmp(InfoAtNode.lab{1},'v1x'))
    EC.v1x=find(strcmp(InfoAtNode.lab{1},'v1x'))+4;% give column number of orientation map start
end

for jElt=1:size(DofPos,2) % loop on elements

point=pointers(:,jElt);
if jElt==1;[EltConst.nodeE,EltConst.nodeEt]=get_nodeE(node,NodePos,jElt,InfoAtNode,DofPos,EltConst);
else;EltConst.nodeE=get_nodeE(node,NodePos,jElt,InfoAtNode,DofPos,EltConst);
end
if isfield(EltConst,'N');of_mk('buildndn',point(4),EltConst);end

if any(point(5)==5)  % 3D Geometric nonlinear stiffness matrix

 Debug=0;
 try; 
  if evalin('base','exist(''Debug'',''var'')');Debug=evalin('base','Debug');end
 end

 defe=def(double(DofPos(:,jElt))+1,1); %1xyz 2xyz, ...
 ke=zeros(length(defe));

 % Step 1 : build the local gradient vector at integration points
 % lestate (9 * ninteg) : here u_{i,j}

% constants 
 ind_ts_eg=[1 6 5 6 2 4 5 4 3]; % e_11 22 33 23 31 12
 ind_ts_egt=[1 6 5 6 2 4 5 4 3]'; % e_11 22 33 23 31 12
 if isempty(EltConst.material);
   if isequal(EltConst.DofLabels(:)',{'u','v','w'});
    EltConst.material='Elastic3DNL';
   end
 end
switch EltConst.material
case 'hyperelastic'
 ci_ts_eg=[1 5 9 8 7 4];ci_ts_egt=ci_ts_eg';

 for jW=0:Nw-1; % HERE IS THE LOOP ON INTEGRATION POINTS 

   r1=EltConst.NDN(:,[Nw 2*Nw 3*Nw]+jW+1);
   U_ij=reshape(defe,3,Nnode)*r1; F_ij=U_ij+eye(3);

   [C,I,dIdc,d2I3dcdc,d2I2dcdc]=elemCalc(F_ij); 
   [dWdI,d2WdI2]=EnHeart(integ,constit,I);
   % Equation (15)
   Sigma=reshape(2*dI1dc*dWdI,3,3);

  if Debug==1 % Test of rivlin cube
   L=evalin('base','L');
   test(1)=norm(diag(L)+eye(3)-F_ij);
   test(2)=norm(C-(eye(3)+diag(L))^2);
   r1=0;for ji=1:3;for jj=1:ji-1; r1=r1+(1+L(ji))^2*(1+L(jj))^2;end;end
   test(3)=norm(I-[sum((1+L).^2) r1 prod((1+L).^2)]);
   test(4)=norm(dI2dc- (I(1)*eye(3)-diag((1+L).^2)));
   test(5)= norm(dI3dc- (I(3)*diag((1+L).^(-2))));
   for ji=1:3; 
    r1(ji)=2*(dWdI(1)+dWdI(2)*(I(1)-(1+L(ji))^2) + dWdI(3)*I(3)*(1+L(ji))^-2);
   end
   test(6)=norm(Sigma-diag(r1))/norm(r1);
   if norm(test,'inf')>1e-5; error('Error on state check');end
   assignin('base','Sigma',Sigma);
  end % Test of rivlin cube

   % #d2wde2 Equation (16)
   d2wde2=dWdI(2)*d2I2dcdc+dWdI(3)*d2I3dcdc + ...
    d2WdI2(1)*dI1dc(ci_ts_egt)*dI1dc(ci_ts_eg) + ...
    d2WdI2(2)*dI2dc(ci_ts_egt)*dI1dc(ci_ts_eg) + ...
    d2WdI2(3)*dI3dc(ci_ts_egt)*dI1dc(ci_ts_eg) + ...
    d2WdI2(4)*dI1dc(ci_ts_egt)*dI2dc(ci_ts_eg) + ...
    d2WdI2(5)*dI2dc(ci_ts_egt)*dI2dc(ci_ts_eg) + ...
    d2WdI2(6)*dI3dc(ci_ts_egt)*dI2dc(ci_ts_eg) + ...
    d2WdI2(7)*dI1dc(ci_ts_egt)*dI3dc(ci_ts_eg) + ...
    d2WdI2(8)*dI2dc(ci_ts_egt)*dI3dc(ci_ts_eg) + ...
    d2WdI2(9)*dI3dc(ci_ts_egt)*dI3dc(ci_ts_eg); 
   d2wde2=d2wde2*4;
    %d2wde2(ind_ts_eg,ind_ts_eg)
   if Debug>0 && min(eig(d2wde2))<0; error('Negative constitutive law');end
   % dW/de = dW/di di/de

  %DD(ci_ts_eg,ci_ts_eg)
   Mecha3DInteg(ke,F_ij,d2wde2,Sigma, ...
    def,EltConst.w(:,4),EltConst.jdet,EltConst.NDN,Nnode,Ndof,Nw,jW)

  end % END THE LOOP ON INTEGRATION POINTS 
  %1;assignin('base','d2wde2_hyper',d2wde2);
case 'Elastic3DNL'

 ci_ts_egt=[1 5 9 8 7 4]';
 if size(EltConst.nodeE,2)>4&&(size(constit,1)<10|| ...
         (length(constit)>38&&constit(39+double(point(7)))));
     at=zeros(3);
 else;at=[];
 end
 
 for jW=0:Nw-1; % HERE IS THE LOOP ON INTEGRATION POINTS 

   r1=EltConst.NDN(:,[Nw 2*Nw 3*Nw]+jW+1);
   U_ij=reshape(defe,3,Nnode)*r1; F_ij=U_ij+eye(3);

  i1=EltConst.ConstitTopology{1};%i1=i1(ind_ts_eg,ind_ts_eg);
  if isfield(EltConst,'CTable');
      constit=ConstitInterp(EltConst.nodeE,EltConst.NDN(:,jW+1),constit,EltConst.CTable);
  end
  if size(constit,1)<10
     % table for rho eta E nu G alpha
     r3=constit((1:9)+double(point(6+1)));
     nu=r3(4); E=r3(3); 
     G=r3(6); if G==0; G=E/2/(1+nu);end
     r1=nu./(1-nu); % n/(1-n) 
     r2=E.*(1-nu)./(1+nu)./(1-2*nu); % E(1-n)/(1+n)(1-2*n)
     d2wde2=zeros(6,6);
     d2wde2([1 2 3 7 8 9 13 14 15])= ...
      r2*[1 r1 r1  r1 1 r1  r1 r1 1]';
    d2wde2([22 29 36])=G; % G
    if size(r3,1)>8;at([1 5 9])=r3(8);T0=r3(9);end
  else;d2wde2=constit(double(i1)+double(point(6+1)));
    if ~isempty(at);at=reshape(constit((39:47)+double(point(7))),3,3);
        T0=constit(48+double(point(7)));
    end
  end
  if EC.v1x 
      bas=EltConst.NDN(:,jW+1)'*EltConst.nodeE(:,EC.v1x+(0:5));
      bas=sp_util('basis',bas(1:3),bas(4:6));
      if 1 % check of formulas
       cLG=bas'; F_iGjL=F_ij*cLG'; 
       F_iLjL=cLG*F_ij*cLG'  % F_iLjL-reshape(toGL2(cLG')'*F_ij(:),3,3)
       % F_ij-reshape(toGL2(cLG)'*F_iLjL(:),3,3)
       e_iLjL=(F_iGjL'*F_iGjL-eye(3))/2       
       %e_ij=cLG'*e_iLjL*cLG;
       d2wde2LL=constit(double(i1)+double(point(6+1)))%L
       SigmaLL=reshape(d2wde2LL(ind_ts_eg,ind_ts_eg)*e_iLjL(:),3,3)%L
       Sigma=cLG'*SigmaLL*cLG; %;Sigma=Sigma(:)
       d2wde2=of_mk('translam',d2wde2LL,bas);% go to G since G gradient 
       
       % grad_l = toGL2(bas)'*grad
       DDG=LdDD(F_ij,d2wde2,Sigma);DDG(abs(DDG)<1e-3)=0 %Sigma GG used here
       % This is computed in mkl_utils
       DDLa=LdDD(F_iLjL,d2wde2LL,SigmaLL);DDLa(abs(DDLa)<1e-3)=0 %Sigma LL used here
       %dGL=toGL4(toGL2(cLG')');reshape(dGL*DDLa(:),9,9)-DDG
       tGL=toGL2(cLG'); %tGL*DDLa*tGL'-DDG
       %tGL'*DDG*tGL-DDLa
       F_ij=tGL*DDLa*tGL'; Sigma=zeros(3);
       
       %DDL=toGL2(bas)*DDG*toGL2(bas)';DDL(abs(DDL)<1e-3)=0
       % eval(iigui({'DDLa','DDL','DDG','bas'},'SetInBaseC')) 
      %% [U1,s1]=svd(DDL); [U2,s2]=svd(DDLa); T=U2/U1;T(abs(T)<1e-5)=0;
      % T'*DDLa*T-DDL
      elseif 1==1 
       %% use local basis and post transform to DDG
       d2wde2LL=constit(double(i1)+double(point(6+1))); % Local
       cGL=bas;
       tGL=toGL2(cGL);F_iLjL=reshape(tGL'*F_ij(:),3,3);
       e_iLjL=(F_iLjL'*F_iLjL-eye(3))/2;
       SigmaLL=reshape(d2wde2LL(ind_ts_eg,ind_ts_eg)*e_iLjL(:),3,3);%L
       DDL=LdDD(F_iLjL,d2wde2LL,SigmaLL);
       DD=tGL*DDL*tGL'; F_ij=DD; Sigma=zeros(3);% replace F_ij for correct Mecha3D integ
       DDL(abs(DDL)<1e-3)=0 %Sigma LL used here
       
          
      else % Change d2wde2 to global
       d2wde2=of_mk('translam',d2wde2,bas);
       e_ij=(F_ij'*F_ij-eye(3))/2; % e_ijG=cLG'*e_ij*cLG
       Sigma=d2wde2(ind_ts_eg,ind_ts_eg)*e_ij(:);reshape(Sigma,3,3) % G
      end
  else
    e_ij=(F_ij'*F_ij-eye(3))/2; 
    Sigma=d2wde2(ind_ts_eg,ind_ts_eg)*e_ij(:); % G
       DDG=LdDD(F_ij,d2wde2,Sigma); DDG(abs(DDG)<1e-3)=0
  end
  if  ~isempty(at); % add thermal stress      
    r1=EltConst.NDN(:,jW+1)'*EltConst.nodeE(:,5)-T0; % dT
    r1=at*r1;r1=d2wde2*r1(ci_ts_egt);   % Stress at integration point
    Sigma=Sigma-r1(ind_ts_eg);
  end
  if isempty(gstate)
  elseif isstruct(gstate)
   Sigma=Sigma+reshape(gstate.Y(ind_ts_eg,jW+1,jElt),size(Sigma));
  else
   Sigma=Sigma+gstate(jW*6+ind_ts_eg,jElt);
  end
  % assignin('base','d2wde2_elas',d2wde2);
  %DD(ci_ts_eg,ci_ts_eg)
  Mecha3DInteg(ke,F_ij,d2wde2,reshape(Sigma,3,3), ...
   def,EltConst.w(:,4),EltConst.jdet,EltConst.NDN,Nnode,Ndof,Nw,jW)

 end % END THE LOOP ON INTEGRATION POINTS 
case 'fs_matrix'  % Follower pressure (fluid/structure coupling)

 i1=(-Nnode+1:0)';Axr=zeros(3);Axs=zeros(3);
 ke=zeros(Nnode*4);F=def(:,2);fe=zeros(Nnode*4,1);

 for jW=0:Nw-1; % HERE IS THE LOOP ON INTEGRATION POINTS 

   
   p=EltConst.N(jW+1,:)*defe(EltConst.VectMap(end+i1));
   eDofPos=EltConst.VectMap([ ...
     EltConst.Nnode*1+i1 EltConst.Nnode*2+i1 EltConst.Nnode*3+i1]);
%evalin('caller',['fe_c(model.DOF([' sprintf('%i ',double(DofPos(eDofPos(:),1))+1) ']))'])
   x=defe(eDofPos)+EltConst.nodeE(:,1:3);
   axr=EltConst.Nr(jW+1,:)*x;
   axs=EltConst.Ns(jW+1,:)*x;
   Axr([6 7 2])=axr;Axr([8 3 4])=-axr; 
   dq=zeros(4*Nnode,3);
   dq(Nnode+i1,1)=EltConst.N(jW+1,:)';
   dq(Nnode*2+i1,2)=EltConst.N(jW+1,:)';
   dq(Nnode*3+i1,3)=EltConst.N(jW+1,:)';
   dqr=zeros(3,4*Nnode);
   dqr(1,Nnode+i1)=EltConst.Ns(jW+1,:);
   dqr(2,Nnode*2+i1)=EltConst.Ns(jW+1,:);
   dqr(3,Nnode*3+i1)=EltConst.Ns(jW+1,:);
   ke=ke+dq*(EltConst.w(jW+1,4)*Axr)*dqr;
   Axr([6 7 2])=axs;Axr([8 3 4])=-axs; 
   dqr(1,Nnode+i1)=EltConst.Nr(jW+1,:);
   dqr(2,Nnode*2+i1)=EltConst.Nr(jW+1,:);
   dqr(3,Nnode*3+i1)=EltConst.Nr(jW+1,:);
   ke=ke+p*dq*(-EltConst.w(jW+1,4)*Axr)*dqr;
  % RHS computation : P x,r ^ x,s \delta v
  fe=p*(EltConst.w(jW+1,4)*cross(axr,axs)'*EltConst.N(jW+1,:))';fe=fe(:);
  i2=double(DofPos(eDofPos,jElt))+1; F(i2)=F(i2)+fe;
 end % END THE LOOP ON INTEGRATION POINTS 
 sp_util('setinput',def,F,int32(1:length(def)),int32(2*ones(1)));

otherwise;
  error('''%s'' not a valid EltConst.material',EltConst.material);
end 

% Linear fluid structure coupling (mass[2] and stiffness[1] matrices)
elseif isequal(EltConst.material,'fs_matrix') % - - - - - - - - - - - - -

 ke=zeros(Ndof);
 for jW=0:Nw-1; % HERE IS THE LOOP ON INTEGRATION POINTS 
   % ({u}.{n}) p ds
   for ji=0:2; 
     coef=EltConst.bas(7+ji,jW+1)*EltConst.jdet(jW+1)*EltConst.w(jW+1,4);
     if coef~=0; jj=Nnode*ji+[1:Nnode]; jk=Nnode*3+[1:Nnode];
      if point(5)==1;
       ke(jj,jk)=ke(jj,jk)-coef*EltConst.N(jW+1,:)'*EltConst.N(jW+1,:);
      elseif point(5)==2
       ke(jk,jj)=ke(jk,jj)+coef*EltConst.N(jW+1,:)'*EltConst.N(jW+1,:);
      end
     end
   end
 end
elseif isequal(EltConst.material,'KmitcGen') % - - - - - - - - - - - - -

ke=p_mitc(EltConst,EltConst.nodeE,constit,DofPos);elmap=[];
	
else 
%% #MatrixInteg.standard methodology for MatrixIntegration rule - - - - - - - -

rule=double(EltConst.MatrixIntegrationRule{point(5)});
if isempty(rule);break;end
% check constit topology
%DD=double(EltConst.ConstitTopology{point(5)});DD(find(DD))=constit(DD(find(DD)));disp([jElt DD(1)])
if length(integ)<3; integ(3)=EltConst.Nnode*length(EltConst.DofLabels);end

% Stress computation, possibly bypass k building
% Dof1 Dof2 NDN1 NDN2 Constit StepConstit StepNW NwIni
if isstruct(gstate)&&~isempty(def); 
    if isfield(EC,'v1x');EltConst.v1x=EC.v1x;end;EC=EltConst;
    if jElt==1;
      RO.isEltId=isequal(gstate.Xlab{4},'EltId') ;
      RO.haswjdet=isfield(gstate,'wjdet');
    end
    r1=stress_observe(EC,rule,constit,point,gstate,jElt,RO); EltConst=EC;
    if isempty(r1);rule=[];k=[];end% skip matrix building
end
if isempty(rule); ke=[];
else; ke=zeros(double(integ(3))+Ncondense);
for j2=1:size(rule,1) % terms of rule
  wstep=rule(j2,6);
  for jw=1:double(rule(j2,7));
    i1=rule(j2,5); if i1<0; coef=-1; else;coef=1;end
    i1=abs(i1);i2=jw+rule(j2,8);
    coef=coef*EltConst.w(i2,4)*EltConst.jdet(i2)* ...
     constit((jw-1)*wstep+i1+1+double(point(7)));
    if coef~=0
       for j4=1:NdofPerField
         iblock=(rule(j2,2)+j4-1)*size(ke,1)+rule(j2,1)+1;
            for j3=1:NdofPerField
             ke(iblock)=ke(iblock)+coef* ...
             EltConst.NDN(j3+(jw+rule(j2,3)-1)*NdofPerField)* ...
             EltConst.NDN(j4+(jw+rule(j2,4)-1)*NdofPerField);
             iblock=iblock+1;
            end
       end
    end
  end % jw
end % rule
if size(ke,2)>double(integ(3))+Ncondense;error('Size change');end
end% Volume integration if defined
if size(def,2)>1&&isfield(EltConst,'RhsDefinition')
  rrule=double(EltConst.RhsDefinition);Mrule=size(rrule,1);
  F=def(:,2); 
  Ndof=size(InfoAtNode,1)/EltConst.Nnode;
  jdet=EltConst.jdet;w=EltConst.w(:,4);Nnode=EltConst.Nnode;
  for jRule=1:Mrule
   r3=rrule(jRule,:);
   switch r3(1)
   case 101 % Volume load defined in InfoAtNode
   %101(1) InfoAtNode1(2) InStep(3) NDNOff1(4) 
   % FDof1(5) NDNCol(6) NormalComp(7)
   % w1(8) nwStep(9)
   if isempty(InfoAtNode)||isstruct(InfoAtNode); break;end
    for jw=r3(8)+(0:r3(9)-1)
%     r1=EltConst.NDN(:,jw+r3(4)+1)'* ...
%       InfoAtNode(r3(2)+[0:Nnode-1]*r3(3)+1,jElt);
     r1=EltConst.NDN(:,jw+r3(4)+1)'*EltConst.nodeE(:,5+r3(2));
     r1=r1*jdet(jw+1)*w(jw+1);
     if r3(7)>=0; r1=r1*EltConst.bas(r3(7)+1,jw+1);end
     in1=double(DofPos(EltConst.VectMap(Nnode*r3(5)+[1:Nnode]),jElt))+1;
     F(in1)=F(in1)+  EltConst.NDN(:,jw+r3(6)+1)*r1;
    end % jw
   case 999 % Volume load defined with fields in NodeE, XXXXX
   %101(1) nodeE_col(2) Step==1 NDNOff1(4) 
   % FDof1(5) NDNCol(6) NormalComp(7)
   % w1(8) nwStep(9)
    for jw=r3(8)+(0:r3(9)-1)
     r1=EltConst.NDN(:,jw+r3(4)+1)'*EltConst.nodeE(:,4+r3(2)+1);
     r1=r1*jdet(jw+1)*w(jw+1);
     if r3(7)>=0; r1=r1*EltConst.bas(r3(7)+1,jw+1);end
     in1=double(DofPos(EltConst.VectMap(Nnode*r3(5)+(1:Nnode)),jElt))+1;
     F(in1)=F(in1)+  EltConst.NDN(:,jw+r3(6)+1)*r1;
    end % jw
   otherwise; sdtw('_nb','% Rhs type not supported');
   end
  end % jRule

  sp_util('setinput',def,F,int32(1:length(def)),int32(2*ones(1)));
end % volume integration if needed

end % matrix type

if Ncondense; % possibly condense matrix if needed
 i1=double(EltConst.CondenseInd)+1; ke=ke(i1,i1);
 ind=1:size(ke,1)-Ncondense;cind=size(ke,1)-Ncondense+1:size(ke,1);
 ke=ke(ind,ind)-ke(ind,cind)* (ke(cind,cind)\ke(cind,ind));
end

%       ke=of_mk('matrixintegration',jElt,NodePos,Case.Node, ...
%       pointers,integ,constit,gstate, ...
%       elmap,InfoAtNode,EltConst,def);
%       i1=[Case.DofPerElt(jGroup);SymFlag;0;DofPos(:,jElt)];
%       of_mk('asmsparse',k,int32(i1),ke,elmap);
 if isempty(k)
 elseif isfield(k,'Ener')
  in1=double(DofPos(EltConst.VectMap,jElt))+1; in2=(find(in1)); in1=in1(in2);
  vect=zeros(size(DofPos,1),1); if all(EltConst.jdet<0); ke=-ke;end
  for j2=1:size(def,2)
   vect(in2)=def(in1,j2); out1(cEGI(jElt),j2)=real(vect'*ke*vect);
  end
 elseif ~isempty(ke); % assemble element matrix if needed
  i1=[AsmPoint;DofPos(:,jElt)]; of_mk('asmsparse',k,int32(i1),ke,elmap);
 end
 out=ke;
end % loop on elements
if isa(k,'struct')
  sp_util('setinput',k.Ener,out1, ...
   int32(1:size(k.Ener,1)),int32(1:size(k.Ener,2)));
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Matlab implementation of single element stress computation call. 
%% #StressObserve (stress_og) --------------------------------------------------
elseif comstr(Cam,'stress_og');

[NodePos,jGroup,Case,def,pointers,integ,constit,EC,InfoAtNode]= ...
  deal(varargin{:});
if isfield(EC,'Tying') 
 % Case that needs reinterpolation (call the element assembly)
 C1=Case;C1.GroupInfo=C1.GroupInfo(jGroup,:);
 C1.GroupInfo{end}.defe=zeros(24,1); mo1=evalin('caller','model');
 mo1=fe_case(mo1,'reset');mo1.Elt=feutil('selelt group',mo1,jGroup);
 [C1,m1.DOF]=fe_mknl('init',mo1);C1.GroupInfo{end}.defe=zeros(24,1);
 C1=fe_mknl('gstate-struct',mo1,C1); % fill gstate structure
 C1.GroupInfo{5}.Y(1,1,1,size(def.def,2))=0; % multi-def
 C1.GroupInfo{8}.defe=zeros(24,size(def.def,2));
 fe_mknl('assemble',mo1,C1,def,1);
 out=reshape(C1.GroupInfo{5}.Y,[],size(C1.GroupInfo{1},2),size(def.def,2)); 
 % stress (Nc*Nw,Nelt,Ndef)
 out1=EC.Nw;out2=[];
 return;
end
node=Case.Node;
if comstr(version,'7'); rule=EC.StressRule{1};
else; rule=double(EC.StressRule{1});
end
out=[];DofPos=[]; 
if ~isfield(EC,'nodeEt'); EC.nodeEt=zeros(20,1,'int32');  end
%EC.gstate=zeros(9,EC.Nw)
for jElt=1:size(NodePos,2)

 EC.nodeE=get_nodeE(node,NodePos,jElt,InfoAtNode,DofPos,EC);
 point=int32(pointers(:,jElt));Ns=length(EC.StrainLabels{1});
 if size(EC.N,2)>size(EC.NDN,1) % ignore additional shape functions
    i1=size(EC.NDN,1)+1;EC.N(:,i1:end)=[];EC.Nr(:,i1:end)=[];
    if isfield(EC,'Ns');EC.Ns(:,i1:end)=[];end
    if isfield(EC,'Nt');EC.Nt(:,i1:end)=[];end
 end
 
ke=[];nw=EC.Nw;
Nnode=min(EC.Nnode,size(EC.N,2));Ndof=numel(EC.VectMap)/Nnode;
i3=point(4)==23&&max(rule(:,4))==2;  % project on surface 
state=zeros(Ndof*Nnode,1); %defE
i5=double(Case.GroupInfo{jGroup,1}(:,jElt))+1;i4=find(i5);i5=i5(i4);
%model=evalin('caller','model');fe_c(model.DOF(i3))
state(i4,1:size(def.def,2))=def.def(i5,:);

if point(5)==5&&isfield(EC,'material')&&strcmpi(EC.material,'Elastic3DNL')
  EC.gstate=zeros(6,EC.Nw); 
  ke=of_mk('StressObserve',EC,int32(rule),constit,EC.nodeE,point,state);
  %if 1==2 % Verify formulas (missing transpose checks)
  %   q=reshape(state,3,[])';
  %   F=eye(3)+q'*EC.NDN(:,9:8:end); e=(F'*F-eye(3))/2;
  %   [reshape(q'*EC.NDN(:,9:8:end),[],1) F(:) e(:) ...
  %       reshape((F*F'-eye(3))/2,[],1) reshape((F*F')/2,[],1)]
  %end  
else
 try;
  if sp_util('diag')<=10&&(~i3||isfield(EC,'CTable')); % i3 non supported cases
   ke=of_mk('StressObserve',EC,int32(rule),constit,EC.nodeE,point);
  end
 end
end
if isempty(ke)
 point=double(point);
 ke=zeros(Ns*nw,numel(EC.VectMap));
 % #StressRule [row, NDNBloc, ConstitPos, DofStart, NwTot,NwOff,coef]
 EC.constit=constit(point(7)+(1:size(constit,1)));EC.constit0=EC.constit+0;
 of_mk('buildndn',point(4),EC);% supports constitInterp
 if sdtdef('verm')<805;rule=double(rule);end
 for jConst=1:size(rule,1) % loop on elements of rule
  coef=EC.constit(rule(jConst,3)+1)*rule(jConst,7);
  for jN=1:Nnode;
   for jW=(1:nw)+rule(jConst,6);
    i1=(jW-1-rule(jConst,6))*rule(jConst,5)+rule(jConst,1)+1;
    if i3 % project the 3D field on given vector
      r3=EC.bas(rule(jConst,4)*3+int32(1:3),jW)';
      i2=(1:3)+(jN-1)*Ndof;
      ke(i1,i2)=ke(i1,i2)+EC.NDN(jN,rule(jConst,2)*nw+jW)*coef*r3;
    else % standard computation
     i2=rule(jConst,4)+(jN-1)*Ndof+1;
     ke(i1,i2)= ke(i1,i2)+EC.NDN(jN,rule(jConst,2)*nw+jW)*coef;
    end
   end
  end
 end
end % m-file implementation
if isempty(out);out=zeros(Ns*EC.Nw,size(def.def,2)*size(NodePos,2));end

if size(EC.nodeE,2)>=10&&isfield(InfoAtNode,'lab')&&any(strcmp(InfoAtNode.lab,'v1x'))
   i1=find(strcmp(InfoAtNode.lab,'v1x'))+4;
   r1=EC.N*EC.nodeE(:,i1+(0:5));
   r2=ke*state;r2=reshape(r2,[],EC.Nw,size(r2,2));
   if size(r2,1)==6 % 3d reorientation
    for jW=1:EC.Nw;
         bas=sp_util('basis',r1(jW,1:3),r1(jW,4:6));
         r2(:,jW,:)=elem0('tensortEng3',bas,[])*squeeze(r2(:,jW,:));
    end
   elseif size(r2,1)==9 % Shell reorient
   end
   r2=reshape(r2,[],size(r2,3));out(:,jElt:size(NodePos,2):end)=r2;
elseif size(ke,2)==size(state,1);out(:,jElt:size(NodePos,2):end)=ke*state;
else; out(:,jElt:size(NodePos,2):end)=ke;
end
if nargout<2
else
  if jElt==1; out2=zeros(max(nw,size(EC.jdet,1)),size(NodePos,2));end
  try;out2(:,jElt)=EC.jdet.*EC.w(:,4);end
end
end % loop on elements
if size(def.def,2)>1
 out=reshape(out,size(out,1),size(NodePos,2),size(def.def,2));
end
out1=nw;

elseif comstr(Cam,'initgstate');
%% #InitGtate ------------------------------------------------------------------

 model=varargin{1};Case=varargin{2};
 [EGroup,nGroup]=getegroup(model.Elt);
 out1=struct('DOF',model.DOF,'def',zeros(length(model.DOF),2), ...
   'lab',{{'CurrentState';'RHS'}});
 for j1=1:size(Case.Stack,1)
  switch comstr(Case.Stack{j1,1},-27)
  case 'fvol'  % Place the volume loads at the beginning of InfoAtNode
   data=Case.Stack{j1,3};
   if ~ischar(data.sel); 
    sdtw('_nb','Ignoring FVol(%s) defined by elements',Case.Stack{j1,2});
    break;
   end
   [i1,elt]=feutil(horzcat('findelt',data.sel),model);
   % field specified at nodes
   if isfield(data,'def')&&isfield(data,'DOF')
     if ischar(data.def);r2=elem0('VectFromDir',model,data);
       if ~isequal(r2.DOF,model.DOF);error('Mismatch');end
       r1=r2.def;
     elseif max(fix(data.DOF))==0; % constant field
      r1=zeros(size(model.DOF));
      for j2=1:length(data.DOF); 
        r1(fe_c(model.DOF,data.DOF(j2),'ind'))=data.def(j2);
      end
     else; % field specified at DOFs
      r1=fe_c(model.DOF,data.DOF,data.def(:,1)')';
     end
     for jGroup=1:nGroup;
      if isa(Case.GroupInfo{jGroup,8},'struct')
       % fill in the volume forces
       cEGI = intersect([EGroup(jGroup)+1:EGroup(jGroup+1)-1]',i1);
       if ~isempty(cEGI)
        gstate=struct('data',full(r1.'),'NodePos',...
            int32(double(Case.GroupInfo{jGroup,1})+1), ...
            'lab',{Case.Stack(j1,2)});
        %cEGI=cEGI-EGroup(jGroup);
        if ~isempty(Case.GroupInfo{jGroup,7});
         error('gstate augment not implemented');
        end
        Case.GroupInfo{jGroup,7}=gstate;
       end
      end
     end  % fill in the volume forces
   elseif isfield(data,'dir');
    r2=elem0('VectFromDir',model,data);
    for jGroup=1:nGroup
     if isa(Case.GroupInfo{jGroup,8},'struct')
       % fill in the volume forces
       cEGI = intersect([EGroup(jGroup)+1:EGroup(jGroup+1)-1]',i1);
       if ~isempty(cEGI)
        [ElemF,EGID]= feutil('getelemf',model.Elt(EGroup(jGroup),:),jGroup);
        inode=fe_super('node',ElemF);
        NodePos=model.Elt(cEGI,inode)'; 
        if any(Case.GroupInfo{jGroup,2}(4,1)==[2 12]); ind=1:2;
        else; ind=1:3;
        end
        gstate=reshape(r2(ind,NodePos),length(ind)*size(NodePos,1), ...
          size(NodePos,2));
        cEGI=cEGI-EGroup(jGroup);
        if any(size(Case.GroupInfo{jGroup,7})<size(gstate))
         Case.GroupInfo{jGroup,7}(size(gstate,1),size(gstate,2))=0;
        end
        Case.GroupInfo{jGroup,7}(1:size(gstate,1),cEGI)= ...
        Case.GroupInfo{jGroup,7}(1:size(gstate,1),cEGI)+gstate(:,cEGI);
       end
     end
    end % jGroup

   end
  end
 end
 out=Case;

% -----------------------------------------------------------------------
elseif comstr(Cam,'nl_');[CAM,Cam]=comstr(CAM,4);
%% #Nl --------------------------------------------------------

carg=1;
model=varargin{carg};carg=carg+1;
Case=varargin{carg};carg=carg+1;
if comstr(Cam,'mat')
%% NL_mat 

jGroup=Case.jGroup;
EC=Case.GroupInfo{jGroup,8};
EC.constit=Case.GroupInfo{jGroup,4};
EC.unl=zeros(4*length(EC.DofLabels)*EC.Nw*size(Case.GroupInfo{jGroup,1},2),1);
EC.vnl=zeros(size(EC.unl));
IA=Case.GroupInfo{jGroup,7};
  DofPos=Case.GroupInfo{jGroup,1};

[EGroup,nGroup]=getegroup(model.Elt);
[ElemF,i1,ElemP]= getegroup(model.Elt(EGroup(jGroup),:),jGroup);
NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));

    cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;  
    NodePos=fe_mknl('NodePos',NNode,model.Elt,cEGI,ElemF);

u=zeros(size(Case.mDOF)); iunl=zeros(1,2); ivnl=zeros(1,2);
for jElt=1:length(cEGI)
 [EC.nodeE,EC.nodeEt]=get_nodeE(Case.Node,NodePos,jElt,IA,DofPos,EC); 
 of_mk('buildndn',31,EC);
 EC.nodeG=(EC.N*EC.nodeE)';
 ue=reshape(u(DofPos(EC.VectMap,jElt)+1),[],length(EC.DofLabels));
 %reshape(fe_c(Case.mDOF(DofPos(EC.VectMap,jElt)+1)),[],3)
 ep=ue'*EC.NDN;
 
end
dbstack; keyboard

else; error('NL_%s not implemented');
end
% -----------------------------------------------------------------------
%% #VectFromDir builds FieldAtNode or def from analytic expressions ------------
%  elem0('VectFromDir',model,data,EC);
%  out=elem0('VectFromDirAtNode',model,r1,EC,Case.Node); % sdtweb fe_mknl('OrientMap')
elseif comstr(Cam,'vect{');
  S=sdth.findobj('_sub,~',CAM);
  st=S(2).subs;data=struct('dir',{{}},'DOF',[]);
  carg=1;model=varargin{carg};carg=carg+1;
  if carg<=nargin&&isfield(varargin{carg},'T');Case=varargin{carg};carg=carg+1;
  else; Case=[];
  end
  RO=struct;
  for j1=1:length(st)
   if strncmpi(st{j1},'x',1); data.dir{end+1}=st{j1}(2:end);data.DOF(end+1)=.01;
   elseif strncmpi(st{j1},'y',1); data.dir{end+1}=st{j1}(2:end);data.DOF(end+1)=.02;
   elseif strncmpi(st{j1},'z',1); data.dir{end+1}=st{j1}(2:end);data.DOF(end+1)=.03;
   elseif strncmpi(st{j1},'sel',1); 
     mo1=model;mo1.Elt=feutil(['selelt' st{j1}(4:end)],mo1);mo1.DOF=[];
     if isempty(mo1.Elt); RO.DOF=[];
     else; 
      mo1=fe_case(mo1,'reset');RO.DOF=feutil('getdof',mo1);
     end
   end
  end
  out=elem0('VectFromDirAtDof',model,data,RO.DOF);
  if isfield(Case,'DOF')&&isfield(Case,'mDOF')
      out=feutilb('placeindof',Case.DOF,out);
      out.DOF=Case.mDOF;out.def=Case.T*out.def;
      'xxxeb problem t_contact scldrange'
  end

elseif comstr(Cam,'vectfromdir');[CAM,Cam]=comstr(CAM,12);

RunOpt=[]; 
carg=1;model=varargin{carg};carg=carg+1;
data=varargin{carg};carg=carg+1;
EC=[];if carg<nargin;EC=varargin{carg};carg=carg+1;end
if ~isstruct(EC)&&size(EC,2); RunOpt.DOF=EC;EC=[];end % Possibly DOF for AtDof
if carg<nargin; node=varargin{carg}; carg=carg+1; else; node=model.Node; end

r1=[];r2=[];
% Build the field at nodes (r2) based in .dir field
if ~isfield(data,'dir')&&isfield(data,'def');
 if comstr(Cam,'atnode') % build field at node output
    i1=unique(round(rem(data.DOF,1)*100));
    r1=zeros(length(i1),size(model.Node,1));
    in1=find(data.DOF<1);
    if ~isempty(in1); % if 0.id gives reference value
       for j1=1:length(in1);
           i2=(i1==round(data.DOF(in1(j1))*100));r1(i2,:)=data.def(in1(j1));
       end
    end
    in1=(data.DOF>1);
    for j1=1:length(i1) % affect values when declared
     [r4,ind,r2]=fe_c(model.Node(:,1)+i1(j1)/100,data.DOF(in1),data.def(in1,:)');
     r1(j1,ind)=full(r2(ind));
    end
    r2=struct('data',r1,'NodePos',EC.NodePos, ...
        'lab',{strrep(feutil('stringdof',i1/100),'.','d')});

 elseif comstr(Cam,'atnode') % build field at integration points
     error('Not implemented');
 elseif max(fix(data.DOF))==0; % constant field
      r1=zeros(size(model.DOF));
      for j2=1:length(data.DOF); 
        ind=fe_c(model.DOF,data.DOF(j2),'ind');
        if isempty(ind)&&length(data.DOF)==1
            r2=data;r2.AtNode=zeros(size(node,1),1);r2.AtNode(:,1)=data.def;
        elseif isempty(ind); error('Multiple external fields not handled');
        elseif ischar(data.def) % eval build load
          NNode=sparse(node(:,1),1,1:size(model.Node,1));
          in1=full(NNode(fix(model.DOF(ind))));
          x=node(in1,5);y=node(in1,6);z=node(in1,7);r=sqrt(x.^2+y.^2);
          r1(ind,:)=eval(data.def);
        else;r1(ind,:)=data.def(j2);
        end
      end
 else % field specified at DOFs
      r1=fe_c(model.DOF,data.DOF,data.def')';
 end
 if isempty(r2); r2=data;r2.def=r1; r2.DOF=model.DOF;end
% r2=data.dir(:)'; r2=r2(ones(size(model.Node,1),1),:)';
% .dir {} support analytic evaluation of fields
% .dir [] numeric values of x y z 
elseif (isnumeric(data.dir)&&numel(data.dir)==3)|| ...
    isa(data.dir,'cell')||ischar(data.dir) 
 if isnumeric(data.dir); 
   r2=data.dir(:)'; r2=r2(ones(size(model.Node,1),1),:)';
   data.dir=num2cell(data.dir);
 elseif comstr(Cam,'gstate')
  Case=EC; 
  if carg>nargin-1 % call from fe_mknl
     ind=Case.jGroup;  RO.Ty=1;
  else % From command line (see t_thermal)
   ind=varargin{carg};carg=carg+1; RO.Ty=0;
   [EGroup,nGroup]=getegroup(model.Elt);
   eltid=feutil('eltidfix',model);
   NNode=sparse(node(:,1),1,1:size(node,1));
  end
  for jGroup=ind(:)'    
   EC=Case.GroupInfo{jGroup,8};
   i1=0;
   if RO.Ty
     cEGI=evalin('caller','cEGI');
     eltid=evalin('caller','eltid');
     NodePos=evalin('caller','NodePos');i1=1;
   else
    cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
    NodePos=fe_mknl('NodePos',NNode,model.Elt,cEGI);
   end
   InfoAtNode=Case.GroupInfo{jGroup,7}; 
   DofPos=Case.GroupInfo{jGroup,1};
   jElt=1;
   RunOpt.Rule=Case.GroupInfo{jGroup,2}(1,4);
   [nodeE,EC.nodeEt]=get_nodeE(node,NodePos,jElt,InfoAtNode,DofPos,EC);
   matdes=1;
   %RunOpt.N=EC.N(:,1:EC.Nnode); 
   %if RunOpt.Rule==13&&strcmp(EC.type,'beam1')
   % RunOpt.N=[1-EC.w(:,1) EC.w(:,1)];
   %end
   gstate=struct('X',{{EC.StrainLabels{matdes}(:),EC.w,eltid(cEGI)}}, ...
        'Xlab',{{sprintf('Stress %i',matdes),'Gauss','EltId'}}, ...
        'Y',zeros(length(EC.StrainLabels{matdes}),size(EC.w,1), ...
            length(cEGI)));
   r2=zeros(size(gstate.Y,1),size(gstate.Y,2));
   [i1,i2]=ismember(data.lab,gstate.X{1});
   if ~all(i1); error('Labels %s mismatched',comstr(data.lab(~i1),-30));end
   % build weighted master displacement observation
   for jElt=1:length(cEGI) % loop on contact masters 
    EC.nodeE=get_nodeE(node,NodePos,jElt,InfoAtNode,DofPos,EC);
    of_mk('buildndn',RunOpt.Rule,EC);
    n1=EC.N*nodeE; x=n1(:,1); y=n1(:,2); z=n1(:,3);
    r2(i2,:)=field_eval(data,n1,EC.nodeEt)';
    sp_util('setinput',gstate.Y,r2,numel(r2)*(jElt-1));
   end
   Case.GroupInfo{jGroup,5}=gstate; 
  end % jGroup loop
  out=Case; return; % return properly with Case filled 
 else % AtNode
  r2=field_eval(data,node);
  if ischar(r2); % coming back to get normals
   elt=feutil('addelt',evalin('caller','ElemF'), ...
       evalin('caller','elt(cEGI,:)'));
   if strcmpi(r2,'normal')
    data.MAP=feutilb('shellmapnodepos',node,elt,[]);
    data.MAP.NodeId=reshape(node(evalin('caller','EC.NodePos'),1), ...
       [],size(elt,1)-1);  % size(elt,1)-1,[])'; /!\ not equivalent 
    data.NodeId=reshape(node(EC.NodePos),size(EC.NodePos,1),[]);
   else
      eval(r2); % Callback to generate MAP at Nodes
   end
   [r2,EC.NodePos]=field_eval(data,node);
  end  
  r2=permute(r2,[2 1 3:length(size(r2))]);
 end % AtNode/gstate

 if comstr(Cam,'atnode') % Build InfoAtNode
   r2=struct('data',r2,'NodePos',EC.NodePos,'lab',[]);
   [r3,i1,i2]=unique(r2.data','rows'); % Pack Columns
   r2.data=r3';r2.NodePos=reshape(i2(EC.NodePos),size(EC.NodePos));
   if isfield(data,'lab'); r2.lab=data.lab; else; r2.lab=fe_c(data.DOF);end
   % pack field : eliminate unused values
   [i1,i2,i3]=unique(r2.NodePos(:));
   %r3=r2;r2.data(:,r2.NodePos(:,:))-r3.data(:,r3.NodePos(:,:));norm(ans)
   r2.NodePos=reshape(int32(i3),size(r2.NodePos,1),[]);
   r2.data=r2.data(:,i1);
   
 elseif comstr(Cam,'atdof') % AtDof
   if ~isfield(RunOpt,'DOF');
     if ~isfield(model,'DOF')||isempty(model.DOF);
         model.DOF=feutil('getdof',model);
     end
     RunOpt.DOF=model.DOF;
   end
   r1=repmat(node(:,1)',size(r2,1),1)+repmat(data.DOF(:),1,size(node,1));
   r2=struct('def',reshape(r2,size(r2,1)*size(r2,2),[]),'DOF',r1(:));
   if isempty(strfind(Cam,'used'));r2=feutil('placeindof',RunOpt.DOF,r2);end
   if isfield(data,'sig') % 'sig','dt.5m:Table{0 .5m 1,1 0 0}'
     r3=sdtsys('urnsig',data.sig); r3.X{1}(:,2)=r3.Y;r3.Xlab{1}={'Time';'u'};
     r2.def=r2.def.*r3.Y';r2.data=r3.X{1}; data.data=r3.X{1};data.Xlab={'DOF',r3.Xlab{1}};
   end
   data.DOF=r2.DOF;data.def=r2.def; r2=data;
 end
% constant vector for element
elseif isfield(data,'EltId')&&size(data.dir,2)==length(data.EltId)
 r2=data;r2.dir=r2.dir;
else
  error('Not a supported volume load direction format');
end
if isstruct(r2)
 st1=setdiff(fieldnames(data),[fieldnames(r2);'MAP']);
 for j1=1:length(st1); r2.(st1{j1})=data.(st1{j1}); end
end
out=r2;

%% #EltId2Map - - -2
%  InfoAtNodeMap=elem0('eltid2Map',model,R1,group)
%  R1.X={EltId,lab} R1.Y=values
%  need cross reference to sdtweb fe_mknl('orientmap')
elseif comstr(Cam,'eltid2map');[CAM,Cam]=comstr(CAM,10);
    
carg=1;
model=varargin{carg};carg=carg+1;
MAP=varargin{carg};carg=carg+1;
eltid=feutil('eltidfix;',model);
if carg<=length(varargin);jGroup=varargin{carg};carg=carg+1;else; jGroup=1;end

if isfield(MAP,'il') || isfield(MAP,'bas')
 %% Get the MAP from a pro.MAP entry
 if isa(model,'v_handle');model=model.GetData;end
 if isfield(MAP,'il')
  model.Elt=feutil('selelt proid',model,MAP.il(1));
 end
 Case=fe_mknl('init NoT NoCon',model);
 r2=feutil('getcg',model);
 [EGroup,nGroup]=getegroup(model.Elt);
 for jGroup=1:nGroup
  cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
  if jGroup==1
   out=Case.GroupInfo{jGroup,7};  out.EltId=[];out.vertex=[];out.data=[];
  end
  r3=Case.GroupInfo{jGroup,7}; 
  if isempty(r3); continue;end  
   if ~isfield(out,'lab');out.lab=r3.lab;end
   if ~isequal(out.lab,r3.lab);
       warning('group%i ignored because of lab mismatch',jGroup);
   end
   r3=reshape(r3.data(:,r3.NodePos),size(r3.data,1),size(r3.NodePos,1),[]);
   r3=squeeze(mean(r3,2));
   out.data=[out.data r3];
   out.EltId=[out.EltId;eltid(cEGI)];
   out.vertex=[out.vertex;r2(cEGI,:)];
 end
 out.NodePos=int32(1:length(out.EltId));
 return
end
if ~isfield(MAP,'Y')&&isfield(MAP,'bas'); % Format with EltId,Bas
 r2=MAP;
 MAP=struct('X',{{r2.EltId,{'v1x','v1y','v1z','v2x','v2y','v2z'}'}}, ...
    'Xlab',{{'EltId','lab'}},'Y',r2.bas(:,7:12), ...
    'name','Ael_Frame');
end
[EGroup,nGroup]=getegroup(model.Elt);
[ElemF,i1,ElemP]=getegroup(model.Elt(EGroup(jGroup),:),jGroup);
cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
[i1,i2]=ismember(eltid(cEGI),MAP.X{1});
if any(i1==0);error('Some EltID are not present in the MAP');end
%NodePos=fe_mknl('NodePos',model.Node,model.Elt,cEGI,ElemF);
out=struct('data',MAP.Y(i2,:)', ...
  'NodePos',int32(repmat(1:length(i2),length(fe_super('node',ElemF)),1)),...
  'lab',MAP.X(2),'EltId',eltid(cEGI));
 
%% #Condense -----------------------------------------------------------------
elseif comstr(Cam,'condense')

carg=1;
i1=varargin{carg};carg=carg+1;
ke=varargin{carg};carg=carg+1;
if carg<nargin; me=varargin{carg};carg=carg+1;else; me=[];end

for j1=size(ke,1):-1:i1
 ke(j1,j1)=1/ke(j1,j1); r1=ke(j1,j1); % set pointer kcc_i
 ind=1:j1-1;
 ke(ind,j1)=-ke(ind,j1)*r1; vect=ke(ind,j1); % set pointer k_ic
 ke(ind,ind)=ke(ind,ind)-1/r1*vect*vect'; % do a dger call
 if ~isempty(me); 
  me(ind,ind)=me(ind,ind)-vect*me(j1,ind)-me(ind,j1)*vect'+ vect*me(j1,j1)*vect';
 end
end
out=ke(1:i1-1,1:i1-1); out1=me(1:i1-1,1:i1-1);

%% #TensortT : tensor transformations ----------------------------------------
elseif comstr(Cam,'tensort');[CAM,Cam]=comstr(CAM,8);

% tensor transforms
% general "matrix" tensor transform

% s_il= [TGL_ji * TGL_kl] * s_jk
carg=1;
TGL=varargin{carg};carg=carg+1;
sig=varargin{carg};carg=carg+1;

if carg<=length(varargin) 
% vector of components provided, return transform 
% elem0('tensort',TGL,2,[1 1 1;2 2 2;3 1 2]) % 2D strain
% elem0('tensort',TGL,1,[1 1 1;2 2 2;3 1 2]) % 2D stress
% elem0('tensort',zeros(2),1,[1 1 1;2 2 2;3 1 2]) % 2D stress
% elem0('tensort',zeros(3),2,[1 1 1;2 2 2;3 3 3;4 2 3;5 3 1;6 1 2]) % 3D stress

 list=varargin{carg};carg=carg+1;i3=size(TGL,1);
 out=zeros(max(list(:,1)));if nnz(TGL)==0;out=cell(size(out)); out1='';end 
 a=1; % C or Matlab index style
 for j1=1:size(list,1);
  ja=list(j1,1); ji=list(j1,2); jl=list(j1,3);
  for j2=1:size(list,1)
  jb=list(j2,1); jj=list(j2,2); jk=list(j2,3);
  if iscell(out); 
      out{ja,jb}=[out{ja,jb} sprintf('T[%i]*T[%i]',(jj-a)+(ji-1)*i3,(jk-a)+(jl-1)*i3)];
  else; out(ja,jb)=out(ja,jb)+TGL(jj,ji)*TGL(jk,jl);
  end
  % transform of sigma (where the vector off-diagonal is term)
  if jj~=jk&&sig(1)==1 % account for symmetry
   if iscell(out)
      out{ja,jb}=[out{ja,jb} sprintf('+T[%i]*T[%i]',(jk-a)+(ji-1)*i3,(jj-a)+(jl-1)*i3)];
   else; out(ja,jb)=out(ja,jb)+TGL(jk,ji)*TGL(jj,jl);
   end
  end
  % transform of epsilon (where off-diagonal is sum of terms)
  if sig(1)==2&&ji~=jl % account for symmetry
   if iscell(out)
      out{ja,jb}=[out{ja,jb} sprintf('+T[%i]*T[%i]',(jj-a)+(jl-1)*i3,(jk-a)+(ji-1)*i3)];
   else; out(ja,jb)=out(ja,jb)+TGL(jj,jl)*TGL(jk,ji);
   end
  end
  end
 end
 if iscell(out) 
    % write transform out'*sig*out A_ij B_jk C_kl
   out1=out'; st=' ';st1={};i3=6;
   for ji=1:size(out,1)
   for jl=1:size(out,1)
      st(end+[0 1])=');';st1{end+1}=st;
      if jl>3; st=sprintf('S[%i]=(',ji-a+(jl-1)*i3);
      else; st=sprintf('S[%i]=(',ji-a+(jl-1)*i3);
      end
   for jj=1:size(out,1)
   for jk=1:size(out,1)
       st=sprintf('%s(%s)*(%s)*%s+',st,out1{ji,jj},out{jk,jl}, ...
           sprintf('IN[%i]',jj-a+(jk-1)*i3));
   end
   end
   end
   end
   clear out
   st(end+[0 1])=');';st1{end+1}=st;st1(1)=[];%st1=sort(st1);
   if nargout==0;
    fprintf('%s\n',st1{:});
    st=sprintf('%s\n',st1{:});% Show in matlab format
    st=textscan(st,'%s','whitespace','[]');st=st{1};
    for j1=1:length(st);
      if all(ismember(st{j1},'0123456789'))
       st{j1}=sprintf('(%i)',str2double(st{j1})+1);
      elseif any(ismember(st{j1},'='));
          st{j1-2}=sprintf('\n%s',st{j1-2});
      end
    end
    fprintf('%s',st{:});fprintf('\n\n');
   else;out=sprintf('%s\n',st1{:});
   end
   
 end
 % transform 3D engineering strain
elseif comstr(Cam,'eng')

 % sig=[s11 s22 s33 s23 s13 s12]
 out=toGLeng3s(TGL);
 
 % epsi=[e11 e22 e33 2*e23 2*e13 2*e12]
 out1=toGLeng3e(TGL');
 if ~isempty(sig)
   out=out'*sig*out1; % xxx 
 elseif comstr(Cam,'eng3s') % transform 3D engineering stress
 elseif comstr(Cam,'eng3e') % transform 3D engineering strain
     out=out1;
 end


elseif isempty(sig)
 %% nominal transform of constitutive law
 out=toGL4(TGL);
else % actually transform the tensor
 N=size(sig,1);out=zeros(N);
 for ji=1:N
     for jj=1:N;
         for jk=1:N;
             for jl=1:N;
  out(ji,jl)=out(ji,jl)+TGL(jj,ji)*TGL(jk,jl)*sig(jj,jk);
             end
         end
     end
 end
end

%% #GaussObserve ---------------------------------------------------------------
% elem0('GaussObserve',rule,integ,constit,model,Case,cEGI)
elseif comstr(Cam,'gaussobserve');
 [CAM,Cam]=comstr(CAM,13);
 [CAM,Cam,RunOpt.Pos]=comstr('-pos',[-25 3],CAM,Cam); 
 RunOpt.wst=sdtw('_off','MATLAB:lang:badlyScopedReturnValue');
 carg=1;
 opt=varargin{carg};carg=carg+1;
 integ=varargin{carg};carg=carg+1;
 constit=varargin{carg};carg=carg+1; 
 model=varargin{carg};carg=carg+1;
 Case=varargin{carg};carg=carg+1;
 if carg<nargin; cEGI=varargin{carg};carg=carg+1;
 else;[EGroup,nGroup]=getegroup(model.Elt);
     jGroup=Case.jGroup; cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
 end
 if isempty(opt);opt=Case.GroupInfo{Case.jGroup,8};
     integ=Case.GroupInfo{Case.jGroup,3};
     constit=Case.GroupInfo{Case.jGroup,4};
 end
 
 if isequal(opt.N,eye(size(opt.N,1))); RunOpt.nI=1; else; RunOpt.nI=0; end
 
 %[EGroup,nGroup]=getegroup(model.Elt);
 %jGroup=min(find(EGroup>cEGI(1)))-1; if jGroup~=Case.jGroup;error('mismatch');end%#ok
 out2=Case.GroupInfo{Case.jGroup,7};InfoAtNode=out2;
 if ischar(opt) % Allow for integrule switching here
  if size(integ,1)>4; 
    opt=integrules(opt,double(integ(5,1)));
    if any(integ(5,:)~=integ(5,1));
     sdtw('Cannot deal with multiple integration strategies in the same group');
    end
  else; opt=integrules(opt);
  end
 end
 % constit contains the element property row (see help p_contact)
 
 %opt=integrules(opt.type,-1); 
 if any(opt.xi(:,3)); RunOpt.Rule=3*ones(1);
 elseif any(opt.xi(:,2)); RunOpt.Rule=23*ones(1);
 else;RunOpt.Rule=13*ones(1);
 end
  
 Nw=opt.Nw; if length(Nw)>1;Nw=Nw(2);end 
 if size(opt.N,1)<Nw;error('Mismatch');end
 %match=struct('Node',zeros(length(cEGI)*Nw,3));
 NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
 DofPos=Case.GroupInfo{Case.jGroup,1}; 
 if any(DofPos(:)==-1);
     i1=model.Elt(cEGI(any(DofPos==-1,1)),1:4);
     warning('Cannot have Fixed Dofs in the contact group');
 end
 II=[];JJ=[];KK=[];
 ji=[-2 zeros(1)];jj=[-2 zeros(1)]; jk=[-2 zeros(1)];jn=[-2 zeros(1)];
 opt.wjdet=zeros(length(cEGI)*Nw,1);
 opt.bas=zeros(9,size(opt.N,1));opt.J=zeros(4,size(opt.N,1));
 
 %Ntheta=zeros(3,Nw*length(cEGI)*size(opt.N,2)*3);r3.in=1;
 r3.eltid=feutil('eltid',model);opt.EltId=zeros(length(cEGI)*Nw,4+RunOpt.nI);
 
 NodePos=fe_mknl('NodePos',NNode,model.Elt,cEGI);node=model.Node;
 jElt=1;[nodeE,opt.nodeEt]=get_nodeE(node,NodePos,jElt,InfoAtNode,DofPos,opt);
 
 if RunOpt.Rule==13&&strcmp(opt.type,'beam1')
   RunOpt.N=[1-opt.w(1:Nw,1) opt.w(1:Nw,1)];
 else;RunOpt.N=opt.N(1:Nw,1:opt.Nnode);
 end
 RunOpt.AtNode=(size(RunOpt.N,1)==size(nodeE,1)&& ...
     isequal(RunOpt.N,eye(size(RunOpt.N)))); 
 match=struct('Node',[],'jdet',[]); if RunOpt.nI; match.ID=[]; end
 RunOpt.iv=[]; 
 if isfield(InfoAtNode,'lab')&&~strcmp(opt.type,'beam1')
     i1=strcmpi(InfoAtNode.lab,'v1x'); RunOpt.iv1=[]; RunOpt.ibas=[];
     if any(i1);RunOpt.v1=cell(size(cEGI));RunOpt.iv1=find(i1)+[4:6];RunOpt.ibas(1,end+1)=1;end
     i1=strcmpi(InfoAtNode.lab,'v2x');RunOpt.iv2=[];
     if any(i1);RunOpt.v2=cell(size(cEGI)); RunOpt.iv2=find(i1)+[4:6];RunOpt.ibas(1,end+1)=2; end
     i1=strcmpi(InfoAtNode.lab,'v3x');RunOpt.iv3=[];
     if any(i1);RunOpt.v3=cell(size(cEGI)); RunOpt.iv3=find(i1)+[4:6];RunOpt.ibas(1,end+1)=3; end
     if RunOpt.Rule==23 % specific, keep normal, v1 projected, v2 tangent
      if ~isempty(RunOpt.iv3);RunOpt.iv=[RunOpt.iv3 RunOpt.iv1]; RunOpt.ibas=[2 3 1];
      else % to ensure use of normal as v3, better skip than generating garbage data
       RunOpt.iv=[];
       sdtw('_nb','skipping MAP for group %i, no v3 for rule 23',Case.jGroup)
      end
     else;  RunOpt.iv=[RunOpt.iv1 RunOpt.iv2 RunOpt.iv3];
     end
     if length(RunOpt.ibas)<3; RunOpt.ibas=[RunOpt.ibas setdiff(1:3,RunOpt.ibas)]; end
     opt.ebas=zeros(9,length(cEGI));
 end
 RunOpt.NField=3;
 % build weighted master displacement observation
 for jElt=1:length(cEGI) % loop on contact masters 
  opt.nodeE=get_nodeE(node,NodePos,jElt,InfoAtNode,DofPos,opt);
  %[nodeE,nodeEt]= comstr(nodeEt,-32)
  of_mk('buildndn',RunOpt.Rule,opt);
  in2=jElt*Nw+[-Nw+1:0];
  n1=RunOpt.N*opt.nodeE; 
  if size(match.Node,1)<length(cEGI)*Nw;match.Node(length(cEGI)*Nw,3)=0;
   if isfield(match,'ID'); match.ID=0*match.Node(:,1); end
  end
  match.Node(in2,:)=n1(:,1:3); if RunOpt.nI; match.ID(in2)=opt.nodeE(:,4); end
  if RunOpt.AtNode; match.ID(in2,1)=opt.nodeE(:,4);end
  match.jdet(in2,:)=opt.jdet(1:Nw);
  i3=ones(size(opt.N,2),1)*((jElt-1)*RunOpt.NField*Nw);
  opt.wjdet(in2)=opt.jdet(1:Nw,1).*opt.w(1:Nw,4);% weight -> press
  
  try;
    if isfield(opt,'FormFcn'); feval(opt.FormFcn,'GaussObserve')
    elseif isfield(opt,'DofLabels')
      i1=1:length(opt.DofLabels)*opt.Nnode;
      i1=reshape(double(DofPos(i1,jElt))+1,length(opt.DofLabels),opt.Nnode)';
      RunOpt.NField=size(i1,2);
      i3=ones(size(opt.N,2),1)*((jElt-1)*RunOpt.NField*Nw);
    else;i1=reshape(double(DofPos(:,jElt))+1,[],opt.Nnode)';
    end
  catch;i1=fix(size(DofPos,1)/opt.Nnode)*opt.Nnode;
      i1=reshape(double(DofPos(1:i1,jElt))+1,[],opt.Nnode)';
  end
  
  for jw=1:Nw % loop on contact points
   jNode=(jElt-1)*Nw+jw; 
   if ~isempty(RunOpt.iv);  % orientation map defined
     bas=(opt.N(jw,:)*opt.nodeE(:,RunOpt.iv));
     try; bas=basis(bas(1:3),bas(4:6));
     catch err;  try; bas=basis(bas(1:3),[0 0 0],1); catch; err.rethrow; end
     end
     opt.ebas(:,jNode)=reshape(bas(:,RunOpt.ibas),9,1); %bas(:);
   elseif isfield(opt,'bas')&&size(opt.bas,2)==Nw
       opt.ebas(:,jNode)=opt.bas(:,jw);
   elseif size(n1,2)>=10  % orientation map defined for contact xxx iv2
    sdtw('_ewt','Is this still used ???')
      bas=sp_util('basis',n1(jw,8:10),n1(jw,5:7));
      opt.ebas(:,jNode)=bas(:);
   end    %if no normal ebas will be given by the match later on
   opt.EltId(jNode,1:4)=[r3.eltid(cEGI(jElt)) match.Node(jNode,:)];
   if RunOpt.nI; opt.EltId(jNode,5)=match.ID(jNode); end
   if isempty(II); % was 3, now size(i1,2)
       II=zeros(1,length(cEGI)*Nw*size(opt.N,2)*size(i1,2));
       JJ=zeros(1,length(cEGI)*Nw*size(opt.N,2)*size(i1,2));
       KK=zeros(1,length(cEGI)*Nw*size(opt.N,2)*size(i1,2));
   end
   if RunOpt.Pos
   elseif RunOpt.Rule==13&&strcmp(opt.type,'beam1')
     r=opt.w(jw);Nu=RunOpt.N(jw,:);
     r1=Nu*opt.nodeE; L=r1(14); bas=reshape(r1(5:13),3,3);
     N_w=opt.N(jw,:).*[1 L 1 L];
     if jw==1;
       i1=i1';opt.wjdet(in2)=L.*opt.w(:,4); match.jdet(in2,:)=L;% u Jacobian
     end
     %[1-3*r.^2+2*r.^3 L*(r-2*r.^2+r.^3)  3*r.^2-2*r.^3   L*(-r.^2+r.^3)];
     b=bas'* ...
      [[Nu(1) 0 0;0 N_w(1) 0;0 0 N_w(1)]*bas [0 0 0;0 0 N_w(2);0 -N_w(2) 0]*bas ...
      [Nu(2) 0 0;0 N_w(3) 0;0 0 N_w(3)]*bas [0 0 0;0 0 N_w(4);0 -N_w(4) 0]*bas];
    sp_util('setinput',KK,b(1,:),jk,'KK');sp_util('setinput',II,i1,ji,'II');
    if length(i3)~=12;i3(end+1:12)=i3(1);end
    i3=i3+1;sp_util('setinput',JJ,i3,jj,'JJ');
    sp_util('setinput',KK,b(2,:),jk,'KK');sp_util('setinput',II,i1,ji,'II');
    i3=i3+1;sp_util('setinput',JJ,i3,jj,'JJ');
    sp_util('setinput',KK,b(3,:),jk,'KK');sp_util('setinput',II,i1,ji,'II');
    i3=i3+1;sp_util('setinput',JJ,i3,jj,'JJ');
    
   else % standard interp displacement for volumes
    % reshape(KK,size(opt.N,2),size(i1,2),Nw,length(cEGI))
    if RunOpt.Rule==3&&jElt==1&&jw==1
     KK=reshape(repmat(reshape(repmat(opt.N(1:Nw,:),1,size(i1,2))',[],1), ...
          1,length(cEGI)),[],1);
     II=reshape(repmat(reshape( ...
        permute(reshape(double(DofPos+1),[],opt.Nnode,length(cEGI)),[2 1 3]), ...
        [],length(cEGI)),Nw,1),[],1);
     JJ=repmat(1:length(cEGI)*size(i1,2)*Nw,size(opt.N,2),1);
     %sparse(z-reshape(JJ(1:numel(z)),size(z)))
    else
    %sp_util('setinput',KK,repmat(opt.N(jw,:)',1,size(i1,2)),jk,'KK');
     for j1=1:size(i1,2)
      sp_util('setinput',KK,opt.N(jw,:),jk,'KK');
      sp_util('setinput',II,i1(:,j1),ji,'II');
      i3=i3+1;sp_util('setinput',JJ,i3,jj,'JJ');
     end
    end
   end
  end % jw
 end % jElt
 ind=ji(2)+1:length(II); II(ind)=[];JJ(ind)=[];KK(ind)=[];
 if ~RunOpt.Pos&&isfield(model,'DOF') % Nodexyz x Gauss xyz
     opt.CIN=sparse(II,JJ,KK,length(model.DOF),length(cEGI)*Nw*RunOpt.NField);
 end
 % fecom('showmap',struct('vertex',match.Node,'normal',opt.ebas(1:3,:)','opt',2))
 sdtw('_set',RunOpt.wst);

 opt.nodeE=zeros(size(nodeE));
 out=opt; out1=match;
%% #Obsolete ---------
%% #DiagConst : obsolete save of Constants (replaced by 4 argument call) ----2
elseif comstr(Cam,'diagconst');
 carg=1;
 if evalin('base','exist(''DiagConst'',''var'')')
   r1=evalin('base','DiagConst');
   if carg+1==length(varargin)
    r1.(varargin{carg})=varargin{carg+1};
   else;r1.(varargin{carg})=varargin(carg+1:end);
   end
   assignin('base','DiagConst',r1);
 else
   disp(varargin{carg});disp(varargin{carg+1});
 end
elseif comstr(Cam,'mooney');error('use elem0(''@EnHeart'')');
% [out,out1]=EnHeart(varargin{:});

%% #end ------------------------------------------------------------------------
elseif comstr(Cam,'cvs')
    out='$Revision: 1.279 $  $Date: 2025/11/18 08:32:49 $'; return;
elseif comstr(Cam,'@');out=eval(CAM);
else; error('''%s'' not supported',CAM);
end

%% #SubFunc
%% feval(elem0('@dispDD'),Case,1)
function dispDD(Case,jgroup)

EC=Case.GroupInfo{jgroup,end};constit=Case.GroupInfo{jgroup,4};
for j1=1:length(EC.ConstitTopology)
 DD=double(EC.ConstitTopology{j1});if isempty(DD); continue;end
 DD(find(DD))=constit(DD(find(DD)));
 disp([DD eig(DD)])
 
end

function out=toGLeng3s(T,s)
out=[T(1)^2 T(2)^2 T(3)^2 2*T(2)*T(3) 2*T(1)*T(3) 2*T(1)*T(2);
      T(4)^2 T(5)^2 T(6)^2 2*T(5)*T(6) 2*T(4)*T(6) 2*T(4)*T(5);
      T(7)^2 T(8)^2 T(9)^2 2*T(8)*T(9) 2*T(7)*T(9) 2*T(7)*T(8);
      T(4)*T(7) T(5)*T(8) T(6)*T(9) (T(5)*T(9)+T(8)*T(6)) (T(4)*T(9)+T(7)*T(6)) (T(4)*T(8)+T(7)*T(5));
      T(1)*T(7) T(2)*T(8) T(3)*T(9) (T(2)*T(9)+T(8)*T(3)) (T(1)*T(9)+T(7)*T(3)) (T(1)*T(8)+T(7)*T(2));
      T(1)*T(4) T(2)*T(5) T(3)*T(6) (T(2)*T(6)+T(5)*T(3)) (T(1)*T(6)+T(4)*T(3)) (T(1)*T(5)+T(2)*T(4));
     ];
if nargin>1; out=out*s;end

function out=toGLeng3e(T,e)
%% 
if isfield(T,'data')
 out=zeros(36*T.Nw,size(T.data,2));
 for j1=1:size(T.data,2)
  %      bas=EC.NDN(:,jW+1)'*EC.nodeE(:,EC.v1x+(0:5)); % if interp in vol
     bas=sp_util('basis',T.data(1:3,j1),T.data(4:6,j1)); % constant in elt
     TenLG=inv(toGLeng3e(bas'));
     out(:,j1)= repmat(reshape(TenLG,[],1),T.Nw,1);
 end
 %out=vhandle.matrix.stressCutDDG(reshape(out,36,[]));
 i1=(1:numel(out)/36)*6;
 [II,JJ,KK]=find(ones(6));
 i1=repmat(i1,6,1)+repmat((-5:0)',1,size(i1,2));
 II=i1(II,:);JJ=i1(JJ,:);
 out=sparse(II,JJ,out(:));
 return
else
 out= [T(1,1)^2 T(1,2)^2 T(1,3)^2 T(1,2)*T(1,3) T(1,1)*T(1,3) T(1,1)*T(1,2);
      T(2,1)^2 T(2,2)^2 T(2,3)^2 T(2,2)*T(2,3) T(2,1)*T(2,3) T(2,1)*T(2,2);
      T(3,1)^2 T(3,2)^2 T(3,3)^2 T(3,2)*T(3,3) T(3,1)*T(3,3) T(3,1)*T(3,2);
      2*T(2,1)*T(3,1) 2*T(2,2)*T(3,2) 2*T(2,3)*T(3,3) (T(2,2)*T(3,3)+T(3,2)*T(2,3)) (T(2,1)*T(3,3)+T(3,1)*T(2,3)) (T(2,1)*T(3,2)+T(3,1)*T(2,2));
      2*T(1,1)*T(3,1) 2*T(1,2)*T(3,2) 2*T(1,3)*T(3,3) (T(1,2)*T(3,3)+T(3,2)*T(1,3)) (T(1,1)*T(3,3)+T(3,1)*T(1,3)) (T(1,1)*T(3,2)+T(3,1)*T(1,2));
      2*T(1,1)*T(2,1) 2*T(1,2)*T(2,2) 2*T(1,3)*T(2,3) (T(1,2)*T(2,3)+T(2,2)*T(1,3)) (T(1,1)*T(2,3)+T(2,1)*T(1,3)) (T(1,1)*T(2,2)+T(2,1)*T(1,2));
     ];
end
if nargin>1; out=out*e;end

function out=toGL2(T,F_ij) 
 % #toGL2 Second order tensor transform
 % T_ik Sigma_kl T_jl
 % F_iLjL=cLG*F_ij*cLG' = reshape(toGL2(bas)'*F_ij(:),3,3)
out=[
T(1,1)*T(1,1) T(1,2)*T(1,1) T(1,3)*T(1,1) T(1,1)*T(1,2) T(1,2)*T(1,2) T(1,3)*T(1,2) T(1,1)*T(1,3) T(1,2)*T(1,3) T(1,3)*T(1,3) 
T(2,1)*T(1,1) T(2,2)*T(1,1) T(2,3)*T(1,1) T(2,1)*T(1,2) T(2,2)*T(1,2) T(2,3)*T(1,2) T(2,1)*T(1,3) T(2,2)*T(1,3) T(2,3)*T(1,3) 
T(3,1)*T(1,1) T(3,2)*T(1,1) T(3,3)*T(1,1) T(3,1)*T(1,2) T(3,2)*T(1,2) T(3,3)*T(1,2) T(3,1)*T(1,3) T(3,2)*T(1,3) T(3,3)*T(1,3) 
T(1,1)*T(2,1) T(1,2)*T(2,1) T(1,3)*T(2,1) T(1,1)*T(2,2) T(1,2)*T(2,2) T(1,3)*T(2,2) T(1,1)*T(2,3) T(1,2)*T(2,3) T(1,3)*T(2,3) 
T(2,1)*T(2,1) T(2,2)*T(2,1) T(2,3)*T(2,1) T(2,1)*T(2,2) T(2,2)*T(2,2) T(2,3)*T(2,2) T(2,1)*T(2,3) T(2,2)*T(2,3) T(2,3)*T(2,3) 
T(3,1)*T(2,1) T(3,2)*T(2,1) T(3,3)*T(2,1) T(3,1)*T(2,2) T(3,2)*T(2,2) T(3,3)*T(2,2) T(3,1)*T(2,3) T(3,2)*T(2,3) T(3,3)*T(2,3) 
T(1,1)*T(3,1) T(1,2)*T(3,1) T(1,3)*T(3,1) T(1,1)*T(3,2) T(1,2)*T(3,2) T(1,3)*T(3,2) T(1,1)*T(3,3) T(1,2)*T(3,3) T(1,3)*T(3,3) 
T(2,1)*T(3,1) T(2,2)*T(3,1) T(2,3)*T(3,1) T(2,1)*T(3,2) T(2,2)*T(3,2) T(2,3)*T(3,2) T(2,1)*T(3,3) T(2,2)*T(3,3) T(2,3)*T(3,3) 
T(3,1)*T(3,1) T(3,2)*T(3,1) T(3,3)*T(3,1) T(3,1)*T(3,2) T(3,2)*T(3,2) T(3,3)*T(3,2) T(3,1)*T(3,3) T(3,2)*T(3,3) T(3,3)*T(3,3) 
];
if nargin>1; out=reshape(T'*F_ij(:),3,3);end

function out=toGL4(TGL);
 %% #toGL4 4th order tensor transform
 N=size(TGL,1); out=zeros(N^2);
 for ji=0:N-1;
  for jj=0:N-1
   for jk=0:N-1
    for jl=0:N-1;
      out(ji+N*jl+1,jj+jk*N+1)=TGL(jj+1,ji+1)*TGL(jk+1,jl+1);
    end
   end
  end
 end



%% #EnHeart hyperelastic models : movde to  m_hyper EnHyper ------------------
function [dWdI,d2WdI2]=EnHeart(integ,constit,I) %#ok<INUSL>
 error('Moved to m_hyper EnHyper')

%% #ConstitInterp --------------------------------------------------------------
function     constit=ConstitInterp(nodeE,N,constit,r1);
%  of_time('cinterp',nodeE,EltConst.NDN,constit,EltConst.CTable,double(jW));
%  if any(abs(z-constit)>eps*constit);error('Mismatch');end

% table interpolation, sdtweb fe_mat('ctable')
for j1=0:r1(1)-1 % loop on tables
    r2=r1(j1*7+(2:8));
    if r2(5)==-3
     constit(r2(7))=N'*nodeE(:,r2(6));
    else
     X=reshape(r1(r2(4)+[1:2*r2(5)]),r2(5),[]); 
     constit(r2(7))=of_time('lininterp',X,N'*nodeE(:,r2(6)),r2);
    end
end

%% #stress_observe -------------------------------------------------------------
function rule=stress_observe(EC,rule,constit,point,gstate,jElt,RO);

if nargin==0
 eval(iigui({'EC','rule','constit','point','gstate','jElt'},'getincaller'));
end
if jElt==1
  EC=integrules('stressrule',EC);
  if isfield(EC,'v1x')
   InfoAtNode=evalin('caller','InfoAtNode');
   if isfield(InfoAtNode,'lab')
    EC.nodeEt(EC.v1x+(0:length(InfoAtNode.lab)-1))=comstr(InfoAtNode.lab,-32);
   end
  end
  if isfield(EC,'nodeEt');EC.nodeEt(end+1:size(EC.nodeE,2))=0;
  end
  assignin('caller','EC',EC);
end
rule=point(5); 
if rule>length(EC.StressRule);
 fprintf('No StressRule for MatDes%i\n',rule);return;
end
rule=EC.StressRule{rule}; 
if isfield(EC,'FormFcn'); 
  ke=feval(EC.FormFcn,'stress_observe',EC);
else % Expecting Ke (stress x ng ) x DofPos order
    % DofPos is usually (fields@n1 , fieldns@nNode) 
 ke=of_mk('StressObserve',EC,int32(rule),constit,EC.nodeE,point);
end

if RO.isEltId;%isequal(gstate.Xlab{4},'EltId') % observation matrix

 i2=size(gstate.Y);i2(end)=1;
 sp_util('setinput',gstate.Y,ke,(jElt-1)*prod(i2));
 if RO.haswjdet;%isfield(gstate,'wjdet');
     sp_util('setinput',gstate.wjdet,EC.jdet.*EC.w(:,4),(jElt-1)*size(gstate.wjdet,1));
 end
 rule=[];%assignin('caller','rule',[]) % Skipp matrix building in this case
end


%% #Mecha3DInteg -------------------------------------------------------------
%% #LdDD the transport of d2wde2 to initial configuration : sdtweb feform#fehyper -browser 
function DD=LdDD(F_ij,d2wde2,Sigma)
   ind_ts_eg=[1 6 5 6 2 4 5 4 3];
   % energy : (dW_n,j F_ni) (d2wde)_(ij kl) (dW_m,l F_mk)
   % strain indices n,j m,l, DD indices  i k
   % DD= [F_ni (d2wde)_(ij kl) F_mk]_(nj ml)
   DD=zeros(9); 
   for jn=0:2;
   for jj=0:2;
   for jm=0:2;
   for jl=0:2; 
     for ji=0:2;
     for jk=0:2
      DD(jn+3*jj+9*jm+27*jl+1)=DD(jn+3*jj+9*jm+27*jl+1)+ ...
       F_ij(jn+ji*3+1)* ...
       d2wde2(ind_ts_eg(ji+3*jj+1)+6*(ind_ts_eg(jk+3*jl+1)-1))* ...
       F_ij(jm+jk*3+1);
     end
     end
     if jn==jm; 
       DD(jn+3*jj+9*jm+27*jl+1)=DD(jn+3*jj+9*jm+27*jl+1)+ ...
        Sigma(jj+3*jl+1);
     end
   end
   end
   end
   end

function out=Mecha3DInteg(ke,F_ij,d2wde2,Sigma,def,w,jdet,NDN,Nnode,Ndof,Nw,jW)

try
   %ind_ts_eg=[1 6 5 6 2 4 5 4 3];
   %d2wde2=d2wde2*0;Sigma=evalin('base','Sigma0');Sigma=Sigma(ind_ts_eg);
   %ind_eg_ts=[1 5 9 6 3 2];fprintf('%.8f \n',Sigma(ind_eg_ts))
  Be=zeros(size(ke,1),1);
  if ~isempty(Be);error('No intention to reimplement');end
  %mEC=evalin('caller','EltConst');mEC.ke=ke; mEC.Be=Be;
  %of_mk('Mecha3DInteg',mEC,F_ij,d2wde2,Sigma,int32([Nnode Ndof Nw jW]));
catch
  if numel(F_ij)==81; DD=F_ij;F_ij=zeros(3); % Possibly provide DD
  else;
    DD=LdDD(F_ij,d2wde2,Sigma);
  end
%1;
if isempty(ke);out=DD; return;end
try;
   error('Not active');
    for ji=0:2
    for jj=0:2
    for jk=0:2
    for jl=0:2;
    coef=DD(ji+jj*3+jk*9+jl*27+1)*jdet(jW+1)*w(jW+1);
    of_mk('k<-k+a*x*y',ke,NDN,NDN,coef, ...
     int32([Nnode*ji+Ndof*Nnode*jk ... % block in stiffness matrix
            Nnode*[Nw*(1+jj)+jW Nw*(1+jl)+jW]   ... % columns in NDN
            Ndof Nnode Nnode 1 1]));
    end
    end
    end
    end
  if 1==2 
   %% Generate observation matrix
   c=zeros(size(DD,1),size(ke,2));
   for ji=0:2
   for jj=0:2
   for jk=0:2
   for jl=0:2;
    in1=Nnode*ji+[1:Nnode];in2=Nnode*jk+[1:Nnode];
    x=NDN(:,Nw*(jj+1)+jW+1);y=NDN(:,Nw*(jl+1)+jW+1);
    c(1+jk+3*jl,in2)=c(1+jk+3*jl,in2)+y';
   end
   end
   end
   end
   i1=evalin('caller','EltConst.VectMap');
   c=evalin('caller','toGL2(bas)')*c(:,i1);
   eval(iigui({'c'},'SetInBaseC'))
   evalin('base','full(MJ.data.c)')
  end
catch
   for ji=0:2
   for jj=0:2
   for jk=0:2
   for jl=0:2;
    coef=DD(ji+jj*3+jk*9+jl*27+1)*jdet(jW+1)*w(jW+1);
    in1=Nnode*ji+[1:Nnode];in2=Nnode*jk+[1:Nnode];
    x=NDN(:,Nw*(jj+1)+jW+1);y=NDN(:,Nw*(jl+1)+jW+1);
    ke(in1,in2)=ke(in1,in2)+coef*x*y';
   end
   end
   end
   end
   assignin('caller','ke',ke);
end

   %if norm(k1-k)/norm(k)>eps; error(1);end
end

% Computation of \int_{\Omega} \Sigma : \delta e = \int \Sigma : (F^T.Dpde)
% Sigma_ij F_kj D_ki = r1_ik : (D_ki = N_,i(jW) u_k)
if size(def,2)>1
  F=def(:,2); 
  DofPos=evalin('caller','DofPos(EltConst.VectMap,jElt)');
  DofPos=reshape(double(DofPos)+1,Nnode,Ndof/Nnode);
  r1=Sigma*F_ij';
  for ji=0:2
  for jk=0:2;
   coef=r1(ji+3*jk+1)*jdet(jW+1)*w(jW+1);
   fe=NDN(:,Nw*(ji+1)+jW+1);
   F(DofPos(:,jk+1))=F(DofPos(:,jk+1))+coef*fe; % In C direct add to def(:,2)
  end
  end
  sp_util('setinput',def,F,int32(1:length(def)),int32(2*ones(1)));
end

%% #get_nodeE ----------------------------------------------------------------
function [nodeE,nodeEt]=get_nodeE(node,NodePos,jElt,InfoAtNode,DofPos,EC);

persistent LastNodeEt
i1=NodePos(:,jElt);i2=(i1==0);

if any(i2);nodeE=zeros(size(NodePos,1),4);
    i2=~i2;nodeE(i2,:)=node(i1(i2),[5:7 1]);
else;nodeE=node(NodePos(:,jElt),[5:7 1]);
end
nodeEt=[];
if isfield(InfoAtNode,'data') % fields at nodes
  nodeE=[nodeE InfoAtNode.data(:,InfoAtNode.NodePos(:,jElt))']; 
  if nargout==1;return;
  elseif jElt>1;nodeEt=LastNodeEt;
  elseif isfield(EC,'nodeEt') 
    r1=int32(comstr(InfoAtNode.lab,-32));nodeEt(end+[1:length(r1)])=r1;
  end
elseif isempty(InfoAtNode)
elseif size(InfoAtNode,1)==size(DofPos,2)&&size(InfoAtNode,2)==Nnode
  nodeE=[nodeE InfoAtNode(jElt,:)']; % for each element all nodes given
elseif ~isempty(InfoAtNode); 
   warning('OpenFEM:ToDo','get_nodeE InfoAtNode implementation needed');
end
if nargout>1
 if isempty(nodeEt);nodeEt=int32([120 121 122 25673]);end %comstr({'x','y','z','Id'},-32);
 if jElt==1;LastNodeEt=nodeEt;end
 if isfield(EC,'nodeEt'); 
    sp_util('setinput',EC.nodeEt,int32(nodeEt),zeros(1));
 end
end
%% #field_eval Analytic evaluation of field based on nodeE or equivalent
function  [r2,NodePos]=field_eval(data,node,nodeEt)

 if nargin==2
  if isfield(data,'MAP') % A map is provided and used to define other fields
   node(~ismember(node(:,1),reshape(data.MAP.NodeId,[],1)),:)=[];
   NNode=sparse(node(:,1),1,1:size(node,1));
   i1=full(NNode(data.MAP.NodeId));
   NodePos=reshape(full(NNode(data.NodeId)),size(data.NodeId,1),[]);
   for j1=1:length(data.MAP.lab)
      r1=data.MAP.data(j1,data.MAP.NodePos);
      eval(sprintf('%s(i1,1)=r1;',data.MAP.lab{j1}));
   end
  elseif isfield(data,'FieldFcn') % FieldFcn callback is used to generate fields
     r2=data.FieldFcn;return; % Build and call field_eval again
  else
   if iscell(data.dir);  % Robustness to dir as char or value
   elseif ischar(data.dir);data.dir={data.dir};
   else; data.dir=num2cell(data.dir);
   end
   % possible need to define normals
   st=[data.dir{cellfun('isclass',data.dir,'char')}];
   if ~isempty(st)&&~isempty(strfind(st,'v3')); 
     r2='normal';return; % Build data.MAP and call field_eval again
   end 
  end
  x=node(:,5);y=node(:,6);z=node(:,7);r=sqrt(x.^2+y.^2);
  if isfield(data,'SetField'); % Dynamic definition of fields
      for j1=1:2:length(data.SetField);
          eval(sprintf('%s=data.SetField{j1+1};',data.SetField{j1+1}))
      end
  end
 else
    st=comstr(nodeEt,-32);
    for j1=1:length(st);
        if ~isempty(st); eval(sprintf('%s=node(:,j1);',st{j1}));end
    end
 end
 if ~iscell(data.dir); data.dir={data.dir};end
 % For compatibility reasons, vector .data is always a single field
 if min(size(data.dir))==1; data.dir=reshape(data.dir,1,[]);end
 r2=zeros(size(node,1),length(data.dir));
 for j2=1:size(data.dir,1) % Def Columns 
 for j1=1:size(data.dir,2) % Def rows/fields
  r3=data.dir{j2,j1};
  if isnumeric(r3); r2(:,j1,j2)=r3;
  elseif isa(r3,'function_handle');r2(:,j1,j2)=feval(r3,node);
  elseif any(r3=='x'|r3=='y'|r3=='z')
    r3=eval(r3); r3=r3(:);
    if size(r3,1)~=size(r2,1)&&size(r3,1)~=1 
     error('String giving an inconsistent evaluation of volume load at nodes');
    end
    r2(:,j1,j2)=r3(:);
  else
   try; r2(:,j1,j2)=str2num(r3); %#ok<ST2NM>
   catch
    try; r2(:,j1,j2)=eval(r3); 
    catch; r2(:,j1,j2)=feval(r3,node);
    end
   end
  end
 end %j1
 end %j2

function  [C,I,dIdc,d2I3dcdc,d2I2dcdc]=elemCalc(F_ij);
   %% #elemCalc large transform 
   C=F_ij'*F_ij;
   I=[trace(C) (trace(C)^2-sum(C(:).^2))/2 det(C)]; % invariants
   CI=inv(C);

   d2I3dcdc=[0 C(9) C(5) -C(8) 0    0;
             C(9) 0 C(1)    0 -C(7) 0;
             C(5) C(1) 0    0  0   -C(4);
             -C(8) 0 0     -C(1)/2  C(4)/2 C(7)/2
              0 -C(7) 0     C(2)/2 -C(5)/2 C(8)/2
              0  0   -C(4)  C(3)/2 C(6)/2 -C(9)/2];
   d2I2dcdc=[0 1 1 0 0 0;1 0 1  0 0 0;1 1 0 0 0 0;0 0 0 -.5 0 0;
           0 0 0 0 -.5 0;0 0 0 0 0 -.5];
   dI2dc=I(1)*eye(3)-C; dI3dc=I(3)*CI; %#ok<MINV>
   % d2I2dcdc*C(ci_ts_eg)'
   dI1dc=eye(3);

   dIdc=[dI1dc(:) dI2dc(:) dI3dc(:)];

 
 
%% #field_eval Analytic evaluation of field based on nodeE or equivalent
function  data=field_REST2V123(node,elt,data,EC) %#ok<DEFNU>

 NNode=sparse(node(:,1),1,1:size(node,1));
 keyboard
   data.SetField={'v1x',v1(:,1)}; 
   
 %% #get_lambda
 function dd=get_lambda(Case,jGroup,matdes) %#ok<DEFNU>
 
 EC=Case.GroupInfo{jGroup,8};
 if matdes>length(EC.ConstitTopology);dd=[];
 else
  dd=double(EC.ConstitTopology{matdes});
  constit=Case.GroupInfo{jGroup,4};constit(end+1:max(dd(:)))=NaN;
  dd(dd~=0)=constit(dd(dd~=0));
 end
