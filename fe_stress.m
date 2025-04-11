function [out,out1]=fe_stress(varargin)

%FE_STRESS compute stress or energies
%
%       Syntax: result = fe_stress('Command',MODEL,DEF)
%       
%	MODEL can be specified as
%         - four arguments NODE, ELT, PL, IL (see : doc node ... and FE_MK)
%         - a structure with fields .Node, .Elt, .pl, .il 
%         - a database wrapper with those fields
%       DEF (the deformations) can be specified using
%         - two arguments MODE (deformations) and MDOF (DOF definition vector)
%         - a structure with fields .def, .DOF
%
%	Accepted commands are linked to the following calls
%	
%	- RESULT = fe_stress('ener (element selection)',MODEL,DEF)
%	  where the strain and kinetic energies in the elements selected with
%         (element selection) (see the FEMESH FINDELT command for details)
%         are computed for the deformations in DEF and the result is returned
%         in the structure RESULT with fields .StrainE, .KinE and .IndInElt
%         which specifies which elements were selected.
%
%         For backward compatibility, if no element selection is given
%         FE_STRESS returns [StrainE,KinE] as two arguments.
%         To select all elements, use the 'ener groupall' command.
%
%	- RESULT=fe_stress('stress (CritFcn) (rest)',MODEL,DEF,EltSel) 
%         Supported criteria are (CritFcn)
%          sI, sII, sIII principal stresses from max to min
%         Supported restitution are
%          1 (sum at nodes wheighted by number of connected elements)
%          AtCenter, AtNode, ...
%        hexa8('test eig stress') gives an example
%
%       - def=fe_stress('expand',model,Case) extrapolates stresses from
%        integration point to nodes
%
% See <a href="matlab: sdtweb _taglist fe_stress">TagList</a>
%	See also help   feplot, fecom, 
%            demo   d_ubeam

%	Etienne Balmes
%       Copyright (c) 2001-2025 by SDTools and INRIA, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       Use fe_stress('cvs') for revision information


%% check of inputs

CAM=varargin{1}; [CAM,Cam]=comstr(CAM,1); carg=2;
%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*NOSEM>


% check the model input
Case=[];
if carg>nargin && strcmp(Cam,'cvs')
 out='$Revision: 1.71 $  $Date: 2025/04/07 17:07:31 $';return;
elseif carg>nargin;model=[];
else; model=varargin{carg}; carg=carg+1;
end
if ~isempty(model)&&isa(model,'double') % all arguments given
  if size(model,2)~=7; error('improper node definition'); end
  model=struct('Node',model,'Elt',varargin{carg}, ...
  'pl',varargin{carg+1},'il',varargin{carg+2});
  carg=carg+3;

elseif isfield(model,'Y')&&isfield(model,'X')
elseif isstruct(model)
  if ~isfield(model,'Node')||~isfield(model,'Elt')
    error('Not a proper model definition');
  end
  model.pl=fe_mat('getpl',model); model.il=fe_mat('getil',model);
elseif isa(model,'sdth');cf=model; model=cf.mdl;
end
if isa(model,'v_handle');model=model.GetData;end
if isfield(model,'Node');
    NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
end

% check the deformation input
if carg<=nargin 
 def=varargin{carg}; carg=carg+1;
 if isa(def,'double') % all arguments given
   def=struct('def',def,'DOF',varargin{carg}); carg=carg+1;
 elseif isstruct(def)
  % do some checks
 elseif isa(tr,'sdth'); error('Not a supported case');
 end
 if isfield(def,'GroupInfo') % really an expand call
  Case=def;
  if carg<=nargin&&isfield(varargin{carg},'DOF');
      def=varargin{carg}; carg=carg+1;
  else; def=[];
  end
 else
  if isempty(def.DOF); error('No DOF definition vector');end
  if isempty(def.def); error('No deformation given');end
  if size(def.def,1)~=size(def.DOF,1)
   error('DEF.DEF and DEF.DOF are not compatible');
  end
 end
else;def=[];
end

%% #ener ek=fe_stress('ener -matdes 5',mdl,def); ------------------------------
if comstr(Cam,'en') 

[CAM,Cam]=comstr(CAM,'energy','%c');

if isfield(model,'Opt')&&model.Opt(1)==3&&isfield(model,'mind')  % Catch SDT type 3 superelement
 if carg<=nargin&&ischar(varargin{carg});
   eval(sprintf( ...
   'out=upcom(model,''ener%s'',''%s'',def);',CAM,varargin{carg}));
 else;eval(sprintf('out=upcom(model,''ener%s'',def);',CAM));
 end
 return;
end
[CAM,Cam,RunOpt.dens]=comstr('dens',[-25 3],CAM,Cam);
[CAM,Cam,RunOpt.SubEner]=comstr('-sub',[-25 3],CAM,Cam);
[CAM,Cam,RunOpt.Curve]=comstr('-curve',[-25 3],CAM,Cam); % New format for energy
[CAM,Cam,RunOpt.MatDes]=comstr('-matdes',[-25 1],CAM,Cam);
[CAM,Cam,RunOpt.silent]=comstr(';',[-25 3],CAM,Cam);

if ~isempty(RunOpt.MatDes);
elseif strncmpi(Cam,'k',1); RunOpt.MatDes=1; [CAM,Cam]=comstr(CAM,2);
elseif strncmpi(Cam,'m',1); RunOpt.MatDes=2; [CAM,Cam]=comstr(CAM,2);
else; RunOpt.MatDes=1;
end
if isempty(Cam)&&carg<=nargin&&ischar(varargin{carg}); % selection in other arg
 CAM=varargin{carg};carg=carg+1;[CAM,Cam]=comstr(CAM,1);
end
if ~isempty(Cam)
 el0=model.Elt;[ind,model.Elt,ind2] = feutil(['findelt' CAM],model);
end
[eltid,model.Elt]=feutil('eltidfix;',model.Elt);
[EGroup,nGroup]=getegroup(model.Elt);
if nGroup == 0;		error('no element group specified in ELT'); end
opt=0;RO.vol=[];

% Verify that DOF numbering is consistent with mdof
if isfield(Case,'GroupInfo')
else
 [Case,model.DOF]=fe_mknl('init NoGetT nocon nocgl',model);
end
if isfield(Case,'InitFailed')&&Case.InitFailed
    error('fe_mknl init failed');
elseif isequal(def.DOF,model.DOF)
elseif isfield(def,'TR')&&isequal(def.TR.DOF,model.DOF)
 if isa(def,'v_handle'); def=def.GetData; end
 def=fe_def('exp',def); % xxx would need better low level handling for performance
elseif sp_util('issdt'); 
 i1=fe_c(model.DOF,def.DOF,'ind',2);% Check if some model DOFs are not present in def.DOF
 if ~isempty(i1)
  sdtw('_nb','Model elements linked to model.DOF that are not in def.DOF have been removed');
  n1=unique(fix(model.DOF(i1,1)));
  model.Elt=feutil('removeelt withnode',model,n1);
  % Rebuild Case and model.DOF
  [Case,model.DOF]=fe_mknl('init NoGetT nocon nocgl',model);
  [eltid,model.Elt]=feutil('eltidfix;',model.Elt);
 end
 i1=fe_c(def.DOF,model.DOF,'ind'); % Remove def.DOF than are not in model.DOF 
 def=fe_def('subdofind',def,i1);
 % Old line below can add 0 to some DOFs, now elts removal
 % def=feutil('placeindof',model.DOF,def); 
else; def.def=fe_c(model.DOF,def.DOF,def.def')';def.DOF=model.DOF;
end
if RunOpt.MatDes==5; % zero prestress and RHS
 q0=stack_get(model,'curve','StaticState','getdata');
 if isempty(q0);
     def.def=[zeros(size(def.def,1),2) def.def];% 
     def.data=[zeros(2,size(def.data,2));def.data];
 else; if size(q0.def,2)<2;q0.def(1,2)=0;end
     if ~isequal(q0.DOF,model.DOF);
      r1=fe_c(q0.DOF,model.DOF,'dof',1);
      if isequal(r1,model.DOF) % DOF were only permuted, reser ordering
       q0=feutilb('placeindof',model.DOF,q0);
      else; error('Expecting model.DOF in StaticState');
      end
     end
     def.def=[q0.def(:,1:2) def.def];% 
     def.data=[zeros(2,size(def.data,2));def.data];
 end
end

if RunOpt.SubEner % split energies sdtweb fe_caseg('enersub')
  out=fe_caseg('enersub',model,Case,def,RunOpt);return
else
  out=fe_mknl('assemble',model,Case,250+RunOpt.MatDes,def);
  try; % Try computing the volumes/masses
   d1=def;d1.def=zeros(size(def.def,1),1);
   d1.def(fe_c(d1.DOF,.01,'ind'),1)=1; % RB mode for density
   out1=fe_mknl('assemble',model,Case,252,d1);clear d1;
   if any(~isfinite(out1)) % Problem with 10002
    el1=feutil('selelt eltind',model,find(~isfinite(out1)));
    feutil('infoelt',el1)
    error('Volume problem');
   end
   r2=out1(:,end);RO.Mass=sum(r2);
   if RO.Mass==0;warning('No mass, skipping density computation'); 
       r2=ones(size(r2));
   else;r2=r2/RO.Mass; ind=find(r2);RO.vol=r2;
   end
  end
  if RunOpt.dens; % energy density
   out(ind,:)=out(ind,:)./r2(ind,ones(size(out,2),1));
  end
end
if isreal(def.def);% modified 12/05/06  the 1/2 factor was in doc but not software
     out=out/2;
else; out=out/4; % if complex integration over a cycle
end  
if any(any(out<0)); 
    for j1=1:size(out,2);
        r1=out(:,j1); ind=find(r1<0);
        ind(-r1(ind)>norm(r1,'inf')*sqrt(eps))=[];
        if ~isempty(ind); out(ind,j1)=0;end
    end
    if any(any(out<0));% clean up negative energies
        r1=-min(out,[],1)./mean(out,1);
        i1=find(r1);r1=[i1;r1(i1)];
        st=sprintf('  Shape %i : min/mean=%.3e \n',r1);
        if ~RunOpt.silent; sdtw('_nb','Some energies are negative\n%s',st); end
    end
end

ind2=find(eltid);if max(ind2)>size(out,1);ind2=ind2(ind2<size(out,1));end
out=struct('data',out(ind2,:),'EltId',eltid(ind2),'eltind',ind2);
if ~isempty(RO.vol);out.vol=RO.vol(ind2);
    if RO.Mass;out.Mass=RO.Mass;end
end
if RunOpt.MatDes==5
 out.data(:,1:2)=[];% static state energies
end
if RunOpt.Curve; % Update to SDT curve format
    mpid=feutil('mpid',model);
    if ~isfield(def,'data'); def.data=(1:size(def.def,2))'; end
    [r1,def]=fe_def('subref',def);if RunOpt.MatDes==5;def.data(1:2,:)=[];end
    if ~isfield(RO,'vol')||isempty(RO.vol);RO.vol=zeros(length(eltid),1);end
    out=struct('X',{{[eltid(ind2) RO.vol(ind2) mpid(ind2,:)],def.data}}, ...
        'Xlab',{{{'EltId';'vol';'MatId';'ProId';'GroupId'},def.Xlab{2}}}, ...
        'Y',out.data,'Mass',RO.Mass); 
end
switch RunOpt.MatDes
case 1; out.name='Ener k'; case 2;out.name='Ener m'; 
otherwise; out.name=sprintf('Ener %i',RunOpt.MatDes');
end
if RunOpt.dens; out.name=[out.name '_dens'];end

elseif comstr(Cam,'expand')
%% #expand extrapolation field at integration points to stress at nodes ------
% def=fe_stress('expand -FieldExp',model,Case)

[EGroup,nGroup]=getegroup(model.Elt);
mdof=feutil('findnode groupall',model);
NNode=sparse(mdof+1,1,1:length(mdof));
conn=int32(zeros(size(model.Elt,2),size(model.Elt,1)));
if nargin==2||isempty(Case)
    model=fe_case(model,'reset');model=stack_rm(model,'mat');
    [Case,model.DOF]=fe_mknl('initnocon NoT',model);
end

% Matrix Init
opt=[];
for jGroup=1:nGroup;
 [ElemF,opt,ElemP]= feutil('getelemf',model.Elt(EGroup(jGroup),:),jGroup);
 if any(strcmp(ElemF,{'celas','mass1','mass2'}))|| ...
         isempty(Case.GroupInfo{jGroup,8})
     continue;
 end
 cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;
 i1=fe_super('node',ElemF);
 Case.DofPerElt(jGroup)=length(i1);
 i2=model.Elt(cEGI,i1);
 i2=int32(reshape(full(NNode(i2+1)-1),size(i2,1),size(i2,2))');
 Case.GroupInfo{jGroup,1}=i2;  % DofPos
 pointers=Case.GroupInfo{jGroup,2};  % Pointers
 pointers(1,:)=length(i1)^2;  % OutSize 1
 pointers(5,:)=1;             % desired matrix (pseudo-mass)
 pointers(6:7,:)=0;           % integ constit offset
 Case.GroupInfo{jGroup,2}=int32(pointers);  % Pointers
  
 Case.GroupInfo{jGroup,3}=int32([1;1;length(i1)]);  % integ
 Case.GroupInfo{jGroup,4}=1;  % Constit
 Case.GroupInfo{jGroup,6}=reshape(1:length(i1)^2,length(i1),length(i1));
    % ElMap
 
 conn(1:size(i2,1),cEGI)=i2;

 opt=Case.GroupInfo{jGroup,8};  %EltConst
 opt.StrainDefinition={[1 1 1 1 length(i1)]};
 opt.StrainLabels={{'s'}}; opt.ConstitTopology={int32(1)};
 opt.NDNLabels={'',',x',',y',',z'};
 opt=integrules('matrixrule',opt);
 Case.GroupInfo{jGroup,8}=opt;
 
end
if of_mk('mwIndex')==8; fun=@int64; else; fun=@int32;end
Case.MatGraph=of_mk('meshGraph', ...
	      size(mdof,1),nGroup,size(conn,1), ...
	      feval(fun,EGroup),feval(fun,diff(EGroup)-1),feval(fun,Case.DofPerElt), ...
	      feval(fun,conn)); %#ok<FVAL>
k=Case.MatGraph; of_mk( 'fillvalue', k, 0 );

% Loop to assemble pseudo-mass and associated RHS
elt=model.Elt;NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
b=zeros(size(k,1),size(Case.GroupInfo{jGroup,5},3));
out=struct('def',[],'DOF',mdof+.99,'label','Strain extended at nodes', ...
  'lab',{{}});
for jGroup=1:nGroup;

  [ElemF,i1]= feutil('getelemf',model.Elt(EGroup(jGroup),:),jGroup);
  if any(strcmp(ElemF,{'celas','mass1','mass2'}))|| ...
         isempty(Case.GroupInfo{jGroup,8})
     continue;
  end
  EGID=i1(1); 
  inode=fe_super('node',ElemF);
  cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;

  % See documentation elem0#call
  DofPos=Case.GroupInfo{jGroup,1};
  pointers=Case.GroupInfo{jGroup,2};
  integ=Case.GroupInfo{jGroup,3};
  constit=Case.GroupInfo{jGroup,4};
  gstate=Case.GroupInfo{jGroup,5};RO.FieldExp=0;
  if ~isempty(strfind(Cam,'fieldexp'))
    i3=size(gstate.Y);RO.FieldExp=1;
    gstate.Y=reshape(eye(i3(2)*i3(3)),[],i3(2),i3(3));
  end
  elmap=int32(Case.GroupInfo{jGroup,6});
  InfoAtNode=Case.GroupInfo{jGroup,7};
  EltConst=Case.GroupInfo{jGroup,8};
  NodePos=fe_mknl('NodePos',NNode,elt,cEGI,ElemF);
  of_mk('matrixintegration',DofPos,NodePos,Case.Node, ...
       pointers,integ,constit,gstate, ...
       elmap,InfoAtNode,EltConst,[],k,int32([Case.DofPerElt(jGroup);0;0]));
  Nw=EltConst.Nw;
  for jElt=1:length(cEGI)
      % Two ways to store stresses : (jW,jElt,jType) or (jType*jW,jElt)
      EltConst.nodeE=model.Node(NodePos(:,jElt),5:7);
      of_mk('buildndn',pointers(4,jElt),EltConst);
      r1=diag(sparse(EltConst.jdet.*EltConst.w(:,4)))*EltConst.N;
      if isfield(gstate,'Y');
        be=reshape(permute(squeeze(gstate.Y(:,:,jElt,:)),[1 3 2]),[],size(r1,1));
        be=(be*r1)';
      elseif size(gstate,3)>1
       be=(squeeze(gstate(:,jElt,:))'*r1)';
      else
       be=(reshape(gstate(:,jElt),size(gstate,1)/Nw,Nw)*r1)';
      end
      if isempty(out.lab);out.lab={};
       for j2=1:size(be,2);out.lab{j2,1}=sprintf('Field %i',j2);end
      end
      i2=double(DofPos(:,jElt))+1;
      if size(b,2)<size(be,2);b(1,size(be,2))=0;end
      b=b+sparse(i2,1:size(be,1),1,size(b,1),size(be,1))*be;
  end

end
if RO.FieldExp
   clear gstate integ constit
   out=struct('cNodeGauss',sparse(k)\sparse(b),'DOF',mdof);
else
    kd=ofact(k); out.def=kd\b; ofact('clear',kd);
end

%% ---------------------------------------------------------------
elseif comstr(Cam,'thermal');[CAM,Cam]=comstr(CAM,8);
%% #Thermal Computation of thermal induced stresses in volume elements of the *b family

RunOpt.Model=0;
if comstr(Cam,'model')
    [Case,model.DOF]=fe_mknl('init',model);RunOpt.Model=1;
    [eltid,model.Elt]=feutil('eltidfix;',model);
elseif comstr(Cam,'gstate'); Case=def;def=[];
else; [Case,model.DOF]=fe_mknl('init',model);
end
if isempty(def);def=m_elastic('thermaldef',Case);end
[EGroup,nGroup]=getegroup(model.Elt);
if ~isfield(Case,'GroupInfo'); Case=fe_mknl('init',model,Case); end
NNode(Case.Node(:,1))=1:size(Case.Node,1);

% build a defT defined at nodes
defT=zeros(size(Case.Node,1),1);
ind=fe_c(def.DOF,.20,'ind');
defT(NNode(fix(def.DOF(ind))),1)=def.def(ind,1);
elt=model.Elt;node=Case.Node;
ci_ts_eg=[1 5 9 8 7 4]';

for jGroup=1:nGroup  % Loop on groups

  % See documentation elem0#call
  DofPos=Case.GroupInfo{jGroup,1};
  pointers=Case.GroupInfo{jGroup,2};
  integ=Case.GroupInfo{jGroup,3};
  constit=Case.GroupInfo{jGroup,4};
  elmap=int32(Case.GroupInfo{jGroup,6});
  InfoAtNode=Case.GroupInfo{jGroup,7};
  EltConst=Case.GroupInfo{jGroup,8};

  if ~isfield(EltConst,'Nw'); continue;end
  gstate=zeros(6*EltConst.Nw,size(Case.GroupInfo{1,1},2));
  cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;
  %NodePos=fe_mknl('NodePos',NNode,elt,cEGI,ElemF);
  NodePos=reshape(NNode(elt(cEGI,1:EltConst.Nnode)'), ...
   EltConst.Nnode,length(cEGI));

  if size(constit,1)<48; error('Thermal properties not defined');end
  idd=find(EltConst.ConstitTopology{1});

  for jElt=1:size(DofPos,2) %
   point=pointers(:,jElt);
   DD=double(EltConst.ConstitTopology{1});
   DD(idd)=constit(DD(idd)+point(7));
   at=reshape(constit((39:47)+point(7)),3,3);

   EltConst.nodeE=node(NodePos(:,jElt),5:7);

   if (point(4)~=2 && point(4)~=3 && point(4)~=23) 
    error('Not a supported NDN strategy'); 
   end
   of_mk('buildndn',point(4),EltConst);

   Nstrain=size(DD,1);

   for jW=0:EltConst.Nw-1
    % Temperature at integration point
    r1=(defT(NodePos(:,jElt))-constit(48+point(7)))'*EltConst.NDN(:,jW+1);
    % Stress at integration point
    r1=at*r1;r1=DD*r1(ci_ts_eg);
    gstate(jW*Nstrain+(1:Nstrain),jElt)=r1;
   end % loop on integration points
  end % loop on elements
  if RunOpt.Model
      gstate=struct('X',{{EltConst.StrainLabels{1},EltConst.w, ...
          eltid(cEGI)}},'Xlab',{{'Comp','Gauss','EltId'}}, ...
          'Y',reshape(gstate,6,size(EltConst.w,1),[]));
      if size(integ,2)>1;
          error(['Use disjoint groups model.Elt=feutilb(' ...
              '''SeparatebyMat'',model.Elt);'])
      else
       model=feutil(sprintf('setpro%i',integ(2)),model,'gstate',gstate);
      end
  end
  Case.GroupInfo{jGroup,5}=gstate;
  Case.GroupInfo{jGroup,8}.material='Elastic3DNL';
end
b=fe_mknl('assemble NoT',model,Case,103);
if RunOpt.Model
 model=fe_case(model,'DofLoad','ThermalLoad',struct('def',b,'DOF',model.DOF));
 out=model;
else
 Case=stack_set(Case,'DofLoad','ThermalLoad',struct('def',b,'DOF',model.DOF));
 out=Case; out1=model.DOF;
end


elseif comstr(Cam,'st') 
%% #stress : stress computation ----------------------------------------------
st='stress';st(length(CAM)+1:end)=[];CAM(Cam(1:length(st))==st)=[];
[CAM,Cam]=comstr(CAM,1);

CritFcn=stack_get(model,'info','StressCritFcn','getdata');
if isempty(CritFcn); 
    [CAM,Cam,CritFcn]=comstr('critfcn',[-25 4],CAM,Cam);
end

% Attempt to deal with where the stress is evaluated.
RunOpt=struct('typ',1,'header','At node'); DIRS=[];% sum on nodes
[CAM,Cam,RunOpt.Curve]=comstr('-curve',[-25 3],CAM,Cam);
i1=strfind(Cam,'atinteg');
if ~isempty(i1);
 RunOpt.typ=3;CAM(i1+(0:6))='';RunOpt.header='At integ';
end
[CAM,Cam,i1]=comstr('atcenter',[-25 3],CAM,Cam);
[CAM,Cam,RunOpt.MatDes]=comstr('matdes',[-25 2],CAM,Cam);
if isempty(RunOpt.MatDes);RunOpt.MatDes=1;end

if i1; RunOpt.typ=2;RunOpt.header='At center';end
[CAM,Cam,i1]=comstr('atnode',[-25 3],CAM,Cam);if i1;RunOpt.typ=1;end
[CAM,Cam,i1]=comstr('gstate',[-25 3],CAM,Cam);if i1; RunOpt.typ=3; end
[CAM,Cam,RunOpt.gpos]=comstr('gpos',[-25 3],CAM,Cam);
[CAM,Cam,RunOpt.Post]=comstr('-post',[-25 4],CAM,Cam);
[CAM,Cam]=comstr(CAM,1);

%if RunOpt.typ==3
if comstr(Cam,'siii');    CritFcn='r1=Principal(3,r1,[],TensorTopology);'; 
elseif comstr(Cam,'sii');    CritFcn='r1=Principal(2,r1,[],TensorTopology);'; 
elseif comstr(Cam,'sidir');
  CritFcn='[r1,dir]=Principal(21,r1,bas,TensorTopology);DIRS(cEGI(jElt),1:length(dir),jDef)=dir;'; bas=eye(3);
elseif comstr(Cam,'si');    CritFcn='r1=Principal(1,r1,[],TensorTopology);'; 
elseif comstr(Cam,'mises'); CritFcn='r1=Principal(4,r1,[],TensorTopology);'; 
elseif comstr(Cam,'sxx');   CritFcn='r1=Principal(11,r1,[],TensorTopology);'; 
elseif comstr(Cam,'syy');   CritFcn='r1=Principal(12,r1,[],TensorTopology);'; 
elseif comstr(Cam,'szz');   CritFcn='r1=Principal(13,r1,[],TensorTopology);'; 
elseif comstr(Cam,'-comp'); 
 [CAM,Cam,i1]=comstr('-comp',[-25 1],CAM,Cam);
 st=sprintf('% i',i1); CritFcn=sprintf('r1=r1([%s],:);',st);
elseif comstr(Cam,'none');   CritFcn='';
elseif RunOpt.typ==3; 
elseif isempty(CritFcn); sdtw('_nb','Using max principal stress');
end
if isempty(CritFcn)&&RunOpt.typ~=3; 
    CritFcn='r1=Principal(1,r1,[],TensorTopology);'; 
end
if ischar(CritFcn)&&~isempty(CritFcn);
  RunOpt.header(2,1:length(CritFcn))=CritFcn;
end

[eltid,model.Elt]=feutil('eltidfix;',model.Elt);
if carg<=nargin && ischar(varargin{carg})
 [ind,elt] = feutil('findelt',model,varargin{carg});carg=carg+1;
else;elt=model.Elt;
end

%mods=model;mods.Elt=elt;Case=fe_mknl('initnocon',mods);

[EGroup,nGroup]=getegroup(elt);
if nGroup == 0;	error('no element group specified in ELT'); end

NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
mode=zeros(size(model.Node,1),3);

model.Elt=elt;
% Verify that DOF numbering is consistent with mdof
Case=fe_case(model,'getcase');
if ~isfield(model,'DOF')||isempty(model.DOF);model.DOF=feutil('getdof',model);end
Case.DOF=model.DOF; Case.T=speye(length(model.DOF));
Case=stack_rm(Case,'fixdof');Case=stack_rm(Case,'mpc');
[Case,mdof]=fe_mknl('initnocon NoT',model,Case);node=Case.Node;
if isequal(def.DOF,mdof)
elseif sp_util('issdt'); def=feutil('placeindof',mdof,def);
else; def.def=fe_c(mdof,def.DOF,def.def')';def.DOF=mdof;
end

switch RunOpt.typ(1)
case 1 % nodal stress field
 out=zeros(size(model.Node,1),size(def.def,2)); out1=zeros(size(model.Node,1),1);
case 2 % element stress value
 out=zeros(size(elt,1),size(def.def,2)); out1=[];
case 3 % Stress at integration points
otherwise; error('Not a valid stress option');
end
RunOpt.wjdet=cell(nGroup,1);
% loop over element groups - - - - - - - - - - - - - - - - - - - - - -
for jGroup = 1:nGroup

[ElemF,gOpt,ElemP]= feutil('getelemf',elt(EGroup(jGroup),:),jGroup);

EltConst=Case.GroupInfo{jGroup,8};
if size(Case.GroupInfo,2)>7 && ...
  isfield(EltConst,'MatrixIntegrationRule')&& ...
  isfield(EltConst,'StrainDefinition')&&~isequal(ElemP,'beam1')
  try;
   if RunOpt.typ==1;   ;EltConst=integrules('stressrule',EltConst,'node');
   elseif RunOpt.typ==2;EltConst=integrules('stressrule',EltConst,'center');
   else; EltConst=integrules('stressrule',EltConst);
   end
   ElemP='mat_og';
  catch;  
   if sp_util('diag')>=10; rethrow(lasterror);
   else;disp(lasterr);sdtw('_nb','Stress rule problem');
   end
  end
end

cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;
if RunOpt.typ(1)==3; out=zeros(1,length(cEGI),size(def.def,2));out1=[];end
NodePos=fe_mknl('NodePos',node,elt,cEGI,ElemF);
pointers=Case.GroupInfo{jGroup,2};
integ=Case.GroupInfo{jGroup,3};
constit=Case.GroupInfo{jGroup,4};
gstate=Case.GroupInfo{jGroup,5};
elmap=int32(Case.GroupInfo{jGroup,6});
InfoAtNode=Case.GroupInfo{jGroup,7};

if ~isempty(integ)&&~isa(integ,'int32'); error('Improper class');end

prop=fe_super('prop',ElemF);i1=fe_super('node',ElemF);

try;
if gOpt(1)<0 % not a display elements
elseif strcmp(ElemP,'mat_og')
   pointers(5,:)=RunOpt.MatDes(1);
   [stress,Nw,wjdet]=elem0('stress_og',NodePos,jGroup,Case,def, ...
    pointers,integ,constit,EltConst,InfoAtNode); % stress (nc*nw,nelt,ndef)
   RunOpt.wjdet{jGroup}=wjdet;
   Nelt=size(NodePos,2); Ns=size(stress,1)/Nw;
   if isfield(EltConst,'TensorTopology');TensorTopology=EltConst.TensorTopology;
   elseif Ns==6; TensorTopology=[1 6 5;6 2 4;5 4 3];
   elseif Ns==3;  TensorTopology=[1 3;3 2];
   else; TensorTopology=[];
   end

   if RunOpt.typ==3&&isempty(CritFcn); % gstate
     if RunOpt.gpos
      r2=reshape(EltConst.N*reshape(model.Node(NodePos,5:7),size(NodePos,1),[]), ...
        size(EltConst.N,1),[],3);
       Case.GroupInfo{jGroup,5}=struct('stress',stress,'GaussCoor',r2);       
     else;
       Case.GroupInfo{jGroup,5}=stress;
     end
   elseif isempty(strfind(CritFcn,'r1'));eval(CritFcn);
   else% do some post-processing on elements
     for jElt=1:Nelt;
     for jDef=1:size(def.def,2);
      r1=reshape(stress(:,jElt+Nelt*(jDef-1)),Ns,Nw);
      eval(CritFcn); 
      if RunOpt.typ~=3;
      [out,out1]=post_expand(out,out1,NodePos(:,jElt),jDef,cEGI(jElt),r1,RunOpt.typ);
      else; Case.GroupInfo{jGroup,5}(:,jElt,jDef)=r1;
      end
     end
     end
   end % CritFcn or Not
else

switch ElemF

case {'hexa8','penta6','hexa20','tetra4','tetra10','penta15'} % 3-D elements
%% 3-D elements - - - - - - - - - - - - - - - - - - - - - - -

 % hexa8 : iopt=[ndof SizeOfOut 200 Unused telt tcar noref]
 %         for car format is fe_mat('of_mk')
 iopt=[length(fe_super('dof',ElemF)) length(i1)*6 200 0 0 2];
 TensorTopology=[1 4 6;4 2 5;6 5 3];
 for jElt=1:length(cEGI)
    nodeE=Case.Node(NodePos(:,jElt),[5:7 1]); i2=NodePos(:,jElt);
    %i2=NNode(elt(cEGI(jElt),i1));nodeE=node(i2,[5:7 1]);
    
    point=pointers(:,jElt); point(5)=200;
    point(1)=double(integ(double(point(6))+4))*6;
    state=zeros(3*size(NodePos,1),1);
    i3=double(Case.GroupInfo{jGroup,1}(:,jElt))+1;i4=find(i3);i3=i3(i4);
    for jDef=1:size(def.def,2)
     state(i4)=def.def(i3,jDef);
     bi=of_mk(ElemF,int32(point),integ,constit,nodeE,[],state);
     r1=reshape(bi,6,length(i1));
     eval(CritFcn);
     [out,out1]=post_expand(out,out1,i2,jDef,cEGI(jElt),r1,RunOpt.typ);
    end
 end

case {'t3p','q4p','t6p','q5p','q8p'}
%% 2-D elements - - - - - - - - - - - - - - - - - - - - - - -

 % t3p : iopt=[ndof SizeOfOut 200 Unused telt tcar noref]
 %         for car format is fe_mat('of_mk') 
 % 3 for sigma xx yy xy
 
 for jElt=1:length(cEGI)
    nodeE=Case.Node(NodePos(:,jElt),[5:7 1]); i2=NodePos(:,jElt);
    point=pointers(:,jElt); point(5)=200;
    % axi : S_rr, S_rz, S_zz, S_tt
    % 2d  : S_xx, S_yy, S_xy 
    % r1=reshape(bi,3,1);
    if   integ(double(point(6))+7)==3
     point(1)=4;TensorTopology=[1 2 5;2 3 5;5 5 4]; %axi S_rr, S_rz, S_zz, S_tt
    else;point(1)=3; TensorTopology=[1 3;3 2];% plane Sxx yy xy
    end 

    state=zeros(2*size(NodePos,1),1);
    i3=double(Case.GroupInfo{jGroup,1}(:,jElt))+1;i4=find(i3);i3=i3(i4); %#ok<FNDSB>
    for jDef=1:size(def.def,2)
     state=def.def(i3,jDef);
     bi=of_mk(ElemF,int32(point),integ,constit,nodeE, ...
        [0 0 0 ...           % alpha
        zeros(1,10)],state, ... % deformation at element
        zeros(10));  % Ti
     r1=bi(:); eval(CritFcn);
     [out,out1]=post_expand(out,out1,i2,jDef,cEGI(jElt),r1,RunOpt.typ);
    end
 end
case {'dktp'}
%% 2-D plate elements - - - - - - - - - - - - - - - - - - - - - - -

TensorTopology=[1 3;3 2];
for jDef=1:size(def.def,2) 
  r1=feutil('dof2mode',full(def.def(:,jDef)),def.DOF);
  mode(NNode(r1(:,1)),1:size(r1,2)-1)=r1(:,2:end);
  for jElt=1:length(cEGI)
    i2=NNode(elt(cEGI(jElt),i1));nodeE=node(i2,[5:7 1]);
    point=pointers(:,jElt); point(1)=3; point(5)=200;
    bi=of_mk(ElemF,int32(point),integ,constit,nodeE, ...
         0,mode(i2,[3 4 5])',0,0);
    %bi=of_mk(ElemF,int32(point),integ,constit,nodeE, ...
    %    [0],mode(i2,[3 4 5 1])',[0],[0]); % deformation at element
    r1=bi(:);
    eval(CritFcn);
    [out,out1]=post_expand(out,out1,i2,jDef,cEGI(jElt),r1,RunOpt.typ);
  end
end

case {'quad4'}  % this IS WRONG AND FOR TESTING ONLY

TensorTopology=[1 3;3 2];

for jDef=1:size(def.def,2) 
  r1=feutil('dof2mode',full(def.def(:,jDef)),def.DOF);
  mode(NNode(r1(:,1)),1:size(r1,2)-1)=r1(:,2:end);
  pointers(5,:)=210;
  for jElt=1:length(cEGI)
    nodeE=Case.Node(NodePos(:,jElt),[5:7 1]); i2=NodePos(:,jElt);
    [bi,bas]=quad4(nodeE,elt(jElt,:),pointers(:,jElt),integ,constit, ...
        [],mode(i2,1:3)');

    % 2d  : S_xx, S_yy, S_xy 
    % r1=reshape(bi,3,1);
    r1=bi(:);
    eval(CritFcn);
    [out,out1]=post_expand(out,out1,i2,jDef,cEGI(jElt),r1,RunOpt.typ);
  end
end
case {'beam1t'}  % handling using C1.GroupInfo{5}
 mo1=model;mo1.Elt=model.Elt(EGroup(jGroup):EGroup(jGroup+1)-1,:);
 [C1,mo1.DOF]=fe_mknl('initNoT',mo1); 
 d1=feutilb('placeindof',model.DOF,def);
 C1=fe_mknl('gstate-struct',model,C1,d1);
 for j1=1:size(d1.def,2)
  k=fe_mknl('assemble',model,C1,fe_def('subdef',d1,j1),1);
  if j1==1;r1=C1.GroupInfo{5};
  else; r1.Y(:,:,:,j1)=C1.GroupInfo{5}.Y;r1.Xlab{4}='def';
   if isfield(def,'data')&&size(def.data,1)==size(def.def,2)
    r1.X{4}=def.data;
   else;r1.X{4}=(1:size(def.def,2))';
   end
  end
 end
 if RunOpt.typ(1)~=3;
  sdtw('_nb','Using -gstate the only currently accepted option for %s',ElemF);
  RunOpt.typ(1)=3;
 end
 out=r1; % Place in Case.GroupInfo{jGroup,5}
case {'bar1'}  % handling using C1.GroupInfo{5}
 d1=feutilb('placeindof',model.DOF,def); r1=[];
 for jElt=1:length(cEGI)
  p1=pointers(:,jElt);p1(5)=-1;
  nodeE=Case.Node(NodePos(:,jElt),[5:7 1]); i2=NodePos(:,jElt);
  bi=bar1(nodeE,elt(cEGI(jElt),:),p1,integ,constit,elmap);
  idof=Case.GroupInfo{jGroup,1}(:,jElt)+1; bi(idof<1)=[];idof(idof<1)=[];
  if isempty(r1);
   r1=struct('X',{{'Sxxe',eltid(cEGI),1:size(def.def,2)}}, ...
    'Xlab',{{'Comp','EltId','def'}}, ...
    'Y',reshape(bi'*d1.def(idof,:),1,1,size(d1.def,2)));
  else
   r1.Y(:,jElt,:)=bi'*d1.def(idof,:);
  end
 end
 out=r1;
otherwise
 try;
   eCall=feval(ElemF,'call');
 catch
  sdtw('_nb','stress computation not supported %s',ElemF);break;
 end
end; 
 if RunOpt.typ(1)==3; Case.GroupInfo{jGroup,5}=out;end
end % display or not
catch; % elements
  sdtw('_nb','stress computation failed for group %i(%s)',jGroup,ElemF);
end

end % loop on element groups

switch RunOpt.typ(1)  % Post process
case 1         % at node
 i1=find(out1); 
 if ~isstruct(out)
  out(i1,:)=out(i1,:)./out1(i1,ones(size(out,2),1))*size(out,2);
  out=struct('DOF',model.Node(:,1)+.99,'data',out,'header',RunOpt.header);
 else; warning('Problem with stress at node computation');
 end
case 2         % at center
 i1=feutil('eltidfix;',elt);
 i2=find(i1);
 out=struct('EltId',i1(i2),'data',out(i2,:),'header',RunOpt.header);
 if ~isempty(DIRS); out.DIRS=DIRS(i2,:,:); end
case 3; out=Case;
  eltid=feutil('eltidfix;',model);
  if RunOpt.Curve;out=cell(size(Case.GroupInfo,1),3);end
  for jGroup=1:size(Case.GroupInfo,1)
    gstate=Case.GroupInfo{jGroup,5}; EC=Case.GroupInfo{jGroup,8};
    if isfield(gstate,'GaussCoor'); 
        GaussCoor=gstate.GaussCoor;gstate=gstate.stress;
    else; GaussCoor=[];
    end
    if isstruct(gstate)||isempty(EC)||isempty(gstate); continue;end;
    cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;
    i1=size(gstate);
    r1=struct('X',{{EC.StrainLabels{1}(:),EC.w,eltid(cEGI)}}, ...
        'Xlab',{{'Stress',{'r';'s';'t';'w'},'EltId'}}, ...
        'Y',reshape(gstate,[i1(1)/EC.Nw EC.Nw i1(2:end)]));
    if ~isempty(GaussCoor);r1.GaussCoor=GaussCoor;end
    try;r1.wjdet=RunOpt.wjdet{jGroup};end
    if length(i1)>3; r1.X{4}=(1:i1(3))'; r1.Xlab{4}='def';end
    if RunOpt.Curve;
        out(jGroup,1:3)={'curve',sprintf('Stress_g_%i',jGroup),r1};
    else; out.GroupInfo{jGroup,5}=r1;
    end
  end
end
if ~isempty(RunOpt.Post);eval(RunOpt.Post);end
%%
elseif comstr(Cam,'schmid')
%% #Schmid initialisation of schmid tensor for sliding planes of
%  anisotropic cristals

[CAM,Cam]=comstr(CAM,7);
r1={'name','m','n';
    'CubicCenteredFace', ...% Provided by Nicolas Ranc
    [-1, 0,-1,-1,0,1, 0,1,1,-1,1,0; % m
      0,-1, 1, 0,1,1,-1,1,0, 1,0,1;
      1, 1, 0, 1,1,0, 1,0,1, 0,1,1]/sqrt(2), ...
     [1,1,1, 1, 1, 1,-1,-1,-1, 1, 1, 1; % n
      1,1,1,-1,-1,-1, 1, 1, 1, 1, 1, 1;
      1,1,1, 1, 1, 1, 1, 1, 1,-1,-1,-1]/sqrt(3)
    };
i1=strcmpi(r1(:,1),CAM);
if ~any(i1);
   fprintf('Supported cristal types in Schmid command');
   disp(r1);
   return
end
m=r1{i1,2};n=r1{i1,3};
out=struct('X',{{{},[]}},'Xlab',{{'plane','Comp'}}, ...
    'Y',zeros(size(m,2),6),'name',r1{i1,1});
[TensorTopology,out.X{2}]=GetTopo('meca3d'); % Tensor SDT convention
[i1,i2]=ismember(1:6,TensorTopology(:));

for j1=1:size(m,2)
    ten=(m(:,j1)*n(:,j1)'+n(:,j1)*m(:,j1)')/2;% tensor notation
    nm=[n(:,j1);m(:,j1)];
    st={'x','y','z';'x','y','z'}';
    coef=[sqrt(nnz(n(:,j1)))*[1 1 1] sqrt(nnz(m(:,j1)))*[1 1 1]];
    for j2=1:6 % build labels
        if nm(j2)==0;st{j2}='';
        elseif abs(nm(j2)*coef(j2)-1)<1e-5;st{j2}=['-' st{j2}];
        elseif abs(nm(j2)*coef(j2)+1)<1e-5
        else;st{j2}=[sprintf('%.1f',coef(j2)*nm(j2)) st{j2}];
        end
    end
    out.X{1}{j1,1}=sprintf('n %s%s%s,m %s%s%s',st{:});
    ten=ten(i2);ten(4:6)=ten(4:6)*2; % eng. notation for double contraction
    out.Y(j1,:)=ten';
end
if isfield(model,'Y') % transform precomputed stresses
    
 r1=model;r2=out;
 out=r1;out.Xlab{1}='Schmid shear';out.X{1}=r2.X{1};
 i1=size(r1.Y);
 out.Y=reshape(r2.Y*reshape(r1.Y,6,[]),[size(r2.Y,1) i1(2:end)]);
 
elseif nargout==0 % display the result cleanly
 coef=sqrt(nnz(n(:,1))*nnz(m(:,1)));
 %r1=[num2cell((1:size(out.Y,1))') out.X{1} num2cell(out.Y*coef)];
 %comstr(r1,-17,'text',{'%2i','%20s','%3.0f','%3.0f','%3.0f','%3.0f','%3.0f','%3.0f'})
 comstr([[{''};out.X{1}] [out.X{2}';num2cell(out.Y*sqrt(6))]],-17,'tab')
 clear out
  
end

elseif comstr(Cam,'tensortopology'); 
%% #TensorTopology [a,b]=feval(fe_stress('@GetTopo'),'Meca3d')
    if nargout==2;[out,out1]=GetTopo(comstr(CAM,15));
    else;out=GetTopo(comstr(CAM,15));
    end
%% ------------------------------------------------------------------------
elseif comstr(CAM,'@');out=eval(CAM);
else;error([CAM ' : not a available command']);
end


function [TensorTopology,lab]=GetTopo(type);
%% #TensorTopology see sdtweb p_solid('conv'); sdtweb feform#feelas3d
% feval(fe_stress('@GetTopo'),'meca3d')

switch lower(type);
    case {'meca3d','mecha3d'}; % SDT convention
     TensorTopology=[1 6 5;6 2 4;5 4 3];
     lab={'Sxx';'Syy';'Szz';'Syz';'Szx';'Sxy'}; 
     % lab=cellfun(@(x)['S' x(2:end)],p_solid('ConvStrain3D'),'uniform',0);
     % lab(TensorTopology)
     lab=lab(:);
    case {'mec3dm'};  % Verified Apr 13 
     % NASTRAN, Modulef, SAMCEF, ANSYS
     % info changed here should be reflected in formulations.tex
     TensorTopology=[1 4 6;4 2 5;6 5 3];
     lab={'Sxx';'Syy';'Szz';'Sxy';'Syz';'Szx'};
    case {'mec3da'};  % Verified Apr 13 : ABAQUS
      TensorTopology=[1 4 5;4 2 6;5 6 3];
      lab={'Sxx','Syy','Szz','Sxy','Sxz','Syz'}; % lab(TensorTopology)
    case {'piezo3d'}; 
     TensorTopology=[];%1 6 5;6 2 4;5 4 3];
     lab={'Sxx';'Syy';'Szz';'Syz';'Szx';'Sxy';'Ex';'Ey';'Ez'}; 
    case 'meca2d'; TensorTopology=[1 3;3 2];
      lab={'Sxx';'Syy';'Sxy'};
    otherwise
    if isempty(type)
        TensorTopology={'Meca3D','Meca2D'};lab={};
    else
        error('TensorTopology %s unknown',Cam);
    end
end

function [out,dir,lab] = Principal(ind,r1,bas,TensorTopology); 
%% #Principal : usual implementation of stress criteria ----------------------
% C1.Y=feval(fe_stress('@Principal'),4,C1.Y,[],'mecha3d');

if ischar(TensorTopology);TensorTopology=GetTopo(TensorTopology);end
r2=zeros(size(TensorTopology));
if length(size(r1))>2  % Apply for multiple deformations
    RO.dim=size(r1);RO.dim(1)=1;
    r1=reshape(r1,size(r1,1),[]);
    [out,dir,lab] = Principal(ind,r1,bas,TensorTopology);
    out=reshape(out,RO.dim);
    dir=[]; % multi def dir output not tested
    return;
end
persistent mexV
if isempty(mexV)
 i1=of_mk('cvs'); i1=sscanf(i1{1,2},'$Revision: %g');
 if i1>1.212; mexV=2;else; mexV=1;end
end
dir=[];
if size(r1,1)==6 % 3-D stress

  % Sxx Syy Szz Sxy Syz Szx
%   if ind==4; lab='Mises';
%    parfor (j1=1:size(r1,2),4)
%     r2=r1(:,j1);r2=r2(TensorTopology);
%     if isfinite(r2(1));s=sort(eig(r2)); else; s=[]; end %#ok<PFTUS>
%     out(j1)=sqrt(((s(1)-s(2))^2+(s(2)-s(3))^2+(s(3)-s(1))^2)/2);
%    end
%   end
  
  if ind==4&&mexV>1
    out=of_mk('StressCrit',r1,'VonMises');lab='Mises';
  else
    out=zeros(1,size(r1,2));
  for j1=1:size(r1,2)
    r2(1:9)=r1(TensorTopology,j1);
    if ~isfinite(r2(1));out(j1)=NaN;continue;
    else;    s=sort(eig(r2)); 
    end
    switch ind
    case 1; out(j1)=s(1);lab='si';
    case 2; out(j1)=s(2);lab='sii';
    case 3; out(j1)=s(3);lab='siii';
    case {4,41}; out(j1)=sqrt(((s(1)-s(2))^2+(s(2)-s(3))^2+(s(3)-s(1))^2)/2);
        lab='Mises';
    case 11; out(j1)=r1(1,j1);lab='sxx';
    case 12; out(j1)=r1(2,j1);lab='syy';
    case 13; out(j1)=r1(3,j1);lab='szz';
    end
  end
  end
elseif size(r1,1)==4 % 2-D axi stress

  r1(5)=0;out=zeros(1,size(r1,2));% S_rr, S_rz, S_zz, S_tt
  for j1=1:size(r1,2)
    r2(1:9)=r1(TensorTopology,j1);s=sort(eig(r2)); 
    switch ind 
    case 4;
        out(j1)=sqrt(((s(1)-s(2))^2+(s(2)-s(3))^2+(s(3)-s(1))^2)/2);
        lab='Mises';
    case 11; out(j1)=r1(1,j1);lab='sxx';
    case 12; out(j1)=r1(2,j1);lab='syy';
    case 13; out(j1)=r1(3,j1);lab='sxy';
    case 21;  
     [dir,s]=eig(r2);[s,i2]=max(abs(diag(s)));
     dir=dir(:,i2);out(j1)=s(end);
     dir=bas*[dir;0];lab='si+d';
    otherwise; out(j1)=s(ind); lab=sprintf('s%i',ind);
    end
  end


elseif size(r1,1)==3 % 2-D stress

  if any(ind==3); error('Can''t compute sIII for 2-D stress');end
  % Sxx yy xy
  out=zeros(1,size(r1,2));
  for j1=1:size(r1,2)
    r2(1:4)=r1(TensorTopology,j1);s=sort(eig(r2)); 
    switch ind 
    case 4; 
     s(3)=0;out(j1)=sqrt(((s(1)-s(2))^2+(s(2)-s(3))^2+(s(3)-s(1))^2)/2);
     lab='Mises';
    case 11; out(j1)=r1(1,j1);lab='sxx';
    case 12; out(j1)=r1(2,j1);lab='syy';
    case 13; out(j1)=r1(3,j1);lab='sxy';
    case 21;  
     [dir,s]=eig(r2);[s,i2]=max(abs(diag(s)));
     dir=dir(:,i2(end));out(j1)=s(end);
     dir=bas*[dir;0];lab='si+d';
    otherwise; out(j1)=s(ind); lab=sprintf('s%i',ind);
    end
  end
  
else % general stress evalution

  out=zeros(1,size(r1,2));r2=zeros(size(TensorTopology));
  i2=find(TensorTopology(:));i3=TensorTopology(i2);
  for j1=1:size(r1,2)
    r2(i2)=r1(i3,j1);s=sort(eig(r2)); 
    switch ind 
    case 4
        if length(s)>3;error('Not an expected case');end
        out(j1)=sqrt(((s(1)-s(2))^2+(s(2)-s(3))^2+(s(3)-s(1))^2)/2);
        lab='Mises';
    case 11; out(j1)=r1(1,j1);lab='s1';
    case 12; out(j1)=r1(2,j1);lab='s2';
    case 13; out(j1)=r1(3,j1);lab='s3';
    case 21;  
     [dir,s]=eig(r2);[s,i2]=max(abs(diag(s)));
     dir=dir(:,i2);out(j1)=s(end);
     dir=bas*[dir;0];lab='si+d';
    otherwise; out(j1)=s(ind); lab=sprintf('s%i',ind);
    end
  end

end

function   [out,out1]=post_expand(out,out1,i2,jDef,jElt,r1,opt);
%% #PostExpand ------------------------------------------------------------

    switch opt(1)
    case 1     % reexpand to nodes
      if size(r1,1)==1&&length(i2)~=length(r1);r1=mean(r1);end
      out(i2,jDef)=out(i2,jDef)+r1(:);
      out1(i2, 1)=out1(i2, 1)+1;
    case 2     % Give at element center
      out(jElt,jDef)=out(jElt,jDef)+mean(r1);
    case 3     % Give at element int. points
      jElt=evalin('caller','jElt');
      out(1:numel(r1),jElt,jDef)=r1;
    end

% -------------------------------------------------------------------------
% Computation of stress state associated with a deformation mode
%function [out,nw]=rule_in_mfile(NodePos,jGroup,Case,def,pointers,integ,constit,EltConst)
%   ke=rule_in_mfile(NodePos,jGroup,Case,def,pointers,integ,constit,EltConst);

