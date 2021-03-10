function [out,out1,out2]=beam1t(CAM,varargin)

%beam1t 2 node beam with pretension. For non-linear cable statics and dynamics
%
%       This element has an internal state where each colum of 
%       Case.GroupInfo{5} gives the local basis and element tension
%        [bas(:) T]

%	Etienne Balmes, Jean-Michel Leclere
%       Copyright (c) 2001-2015 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

%#ok<*NOSEM,*ASGLU,*NASGU>

persistent getKFcn vof

if ischar(CAM)
 if ~isempty(getKFcn)
 elseif exist('getBeamK_mex','file');getKFcn=@getBeamK_mex;
 else;getKFcn=@getBeamK;
 end
 if isempty(vof); vof=0;end
 
 [CAM,Cam]=comstr(CAM,1);
%% #IntegInfo - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 % Build Constit, Integ, and Elmap for later integration
 if comstr(Cam,'integinfo')% sdtweb p_solid('constitbeam')
  [out,out1,out2]=beam1('integinfo',varargin{:});
%% #GroupInit - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'groupinit') 
           % State  ElementConstants        Reindex
   out='[gstate,Case.GroupInfo{jGroup,8},Case.GroupInfo{jGroup,7}]=beam1t(''BasisAndConstants'',model,node,model.Elt(EGroup(jGroup)+1:EGroup(jGroup+1)-1,:),integ,Case);';

%% #BasisAndConstants - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'basisandconstants') % 

node=varargin{2};elt=varargin{3}; 
integ=varargin{4};Case=varargin{5}; 
constit=evalin('caller','constit');
state=zeros(11,size(elt,1));
NNode=evalin('caller','NNode');%NNode=sparse(node(:,1),1,1:size(node,1));
out2=Case.GroupInfo{Case.jGroup,7};
lab={'v1x','v2x','v3x','v1y','v2y','v3y','v1z','v2z','v3z','L','ten'};
if isempty(out2); % allow propagation of data from gstate
    out2=evalin('caller','gstate');
    if isfield(out2,'data')&& size(out2.data,1)>=11&& ...
       size(out2.data,2)==size(elt,1)
     if ~all(strcmp(lab,out2.lab(1:11)));
         error('Mismatch on InfoAtNodes labels');
     end
     state=out2.data;lab=out2.lab;
    else;out2=[];
    end
end

for jElt=1:size(elt,1) 
 x=node(NNode(elt(jElt,[1 2])),5:7);
 if size(elt,2)<5; x0=[1.5 1.5 1.5];  % if nothing,careful elt(6:7) may exist
 elseif elt(jElt,5)<0||rem(elt(jElt,5),1)||(size(elt,2)>5&&elt(jElt,6))|| ...
        (size(elt,2)>6&&elt(jElt,7)) % normal given
   if size(elt,2)<7;elt(1,7)=0;end;x0 = x(1,:)+elt(jElt,5:7);     
 else;x0 = find(node(:,1)==elt(jElt,5));
  if isempty(x0); x0=[1.5 1.5 1.5];
  else;x0 = node(x0,5:7); end
 end
 L = norm(x(2,:)-x(1,:));
 x=sp_util('basis',x(2,:)-x(1,:),x0-x(1,:))';%cLG
 state(1:9,jElt)=x(:);
 state(10,jElt)=L;
end

out=[]; EC=integrules('beam1',3); % constants
%EC.StrainDefinition{1}=[1 2 1 3 1]; % Should have all components
EC.MatrixIntegrationRule={};
EC.StrainLabels={{'e_xx','s_y','s_z','phi_x','k_y','k_z'}, ...
    {'tx','ty','tz','rx','ry','rz'}};
EC.ConstitTopology={int32(1),[]};
EC.StrainDefinition{2}=[1 1 1 1 3;2 1 2 1 3;3 1 3 1 3];
EC.Be=zeros(12,1);EC.defe=zeros(12,1);EC.VectMap=int32(reshape(1:12,6,2));
EC.eltg=zeros(1,10);
if vof||isfield(out2,'Formulation')&&strcmpi(out2.Formulation,'of_mk')
 EC.material='callback';EC.fHandle=@beam1t;
 EC.VectMap=reshape(1:12,6,2); %EC.Formulation='of_mk'; 
 EC.nodeE=zeros(2,17);EC.nodeEt=zeros(1,size(EC.nodeE,2),'int32');
 EC.ke=zeros(12); EC.Be=zeros(12,1);
end
r2=out2;
out2=struct('data',state,'NodePos',int32([1;1]*[1:size(state,2)]), ...
 'lab',{lab});
EC.DofLabels={'u','v','w','rx','ry','rz'};
if ~isstruct(r2)
elseif size(state,2)~=size(r2.NodePos,2); 
    error('Inconsistent group size')
else; 
  for j1=1:length(r2.lab) % merge given fields, really at elt
      i1=find(strcmpi(r2.lab{j1},out2.lab));
      if isempty(i1);i1=size(out2.data,1)+1;out2.lab{i1}=r2.lab{j1};end
      r3=mean(reshape(r2.data(j1,r2.NodePos),size(r2.NodePos)),1);
      out2.data(i1,:)=r3;
  end
  if isfield(r2,'data0')&&ischar(r2.data0)&&strcmpi(r2.data0,'need')
    out2.data0=out2.data;out2.data0(1:10,:)=0;
  end
end
%sdtweb p_solid('constitbeam')
if constit(6)~=0; % Thermal expansion 
  r2=stack_get(Case,'DofSet','ThermalState','getdata');
  if isempty(r2);r2=stack_get(Case,'DofSet','InfoAtNode','getdata');end
  if ~isempty(r2); % Place thermal field as additional DOF
    out2.lab{12}='T';
    r1=sparse(fix(r2.DOF),1,r2.def);
    out2.data(12,out2.NodePos(1,:))=full(r1(elt(:,1)))';
  end      
end
EC.CTable=[];
[pro,il]=matgui('getstackils',evalin('caller','model'),integ(2));pro=pro{3};
if ~isempty(intersect(constit(8:11),[-2 -3])) % ConstitInterpolation in PRO.(tag) field
 % sdtweb m_elastic.m#541
 pro.EC.nodeEt=int32(comstr([{'x','y','z','ID'} out2.lab],-32));
 pro.EC.MatPropertyUnit={'E','nu','Rho','G','eta','alpha','T0','J','I1', ...%I1(9)
     'I2','A','k1','k2','lump','NSM'}';
 pro.EC.constit=constit;
 if length(pro.il)==1; 
   pro.il=evalin('caller',sprintf('il(il(:,1)==%i,:)',pro.il(1)));
 end
 r1=feval(fe_mat('@field_interp'),pro);
 EC.CTable=r1.CTable;
 EC.material='callback';EC.fHandle=@beam1t;
 EC.VectMap=(1:12)'; 
 EC.nodeEt=pro.EC.nodeEt;EC.nodeE=zeros(2,length(EC.nodeEt));
 EC.ke=zeros(12);
end
try;
 if isfield(pro,'StressOut')
  if size(integ,2)>1; fprintf('Single IL assumed');error('Not supported');end
  mat=feutil(sprintf('GetPl%i -struct1',integ(1)),evalin('caller','model'));
  if length(pro.il)==1; pro.il=il(il(:,1)==pro.il,:);end
  EC.StressOut=p_beam('StressObserve',mat,pro);%mdl,[100 100]);
 end
end
out1=EC;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% #Map returns a normal map for first bending plane
 elseif  comstr(Cam,'map');

  model=varargin{1}; if isa(model,'v_handle');model=model.GetData;end
  model.Elt=feutil('setgroup beam1 name beam1t',model);
  [eltid,model.Elt]=feutil('eltidfix',model);
  Case=fe_mknl('initnocon',model);
  [EGroup,nGroup]=getegroup(model.Elt);
  out=struct('ID',[],'normal',[],'vertex',[],'n2',[]);
  NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
  for jGroup=1:nGroup
     ElemF= feutil('getelemf',model.Elt(EGroup(jGroup),:),jGroup);
     cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;
     if strcmp(ElemF,'beam1t');
      out.ID=[out.ID;eltid(cEGI)]; ind=[1:length(cEGI)];
      out.vertex(end+[1:length(cEGI)],1:3)= ...
       (model.Node(NNode(model.Elt(cEGI,1)),5:7)+ ...
         model.Node(NNode(model.Elt(cEGI,2)),5:7))/2;
      if isempty(Case.GroupInfo{jGroup,5}) % New strategy in InfoAtNode
          r1=Case.GroupInfo{jGroup,7};
          out.normal(end+ind,1:3)=r1.data([2 5 8],r1.NodePos(1,:))';
          out.n2(end+ind,1:3)=r1.data([3 6 9],r1.NodePos(1,:))';
      else
       out.normal(end+ind,1:3)=Case.GroupInfo{jGroup,5}([2 5 8],:)';
       out.n2(end+ind,1:3)=Case.GroupInfo{jGroup,5}([3 6 9],:)';
      end
     end
  end
  if sp_util('issdt')&&nargout==0;cf=feplot;cf.model=model;fecom(cf,'showmap',out);end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% #MatCall - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'matcall');
     out=beam1t('call'); 
     try;
       if nargin>4
        C1=varargin{4};EC=C1.GroupInfo{C1.jGroup,8};
        if isfield(EC,'fHandle');out='mat_og';end
       end
     end
     out1=0; % CallSymFlag
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'call')  % call for matrix assembly
  out='[k1,m1]=beam1t(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,gstate,InfoAtNode.data(:,jElt),Case.GroupInfo{jGroup,8});';
 elseif comstr(Cam,'rhscall') % call for load assembly
  out=     'be=beam1t(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,gstate,InfoAtNode.data(:,jElt),Case.GroupInfo{jGroup,8},defe);';
 elseif comstr(Cam,'state');   % call for state update
   out= ['beam1t(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ, '...
       ' constit,gstate,InfoAtNode,FieldAtEltDof,0);'];

%% #SetProMap : uses particular solution to define tensions as a map - - - - -
elseif comstr(Cam,'setpromap');[CAM,Cam]=comstr(CAM,5); 
  model=varargin{1};def=varargin{2};
  [Case,DOF]=fe_mknl('init',model); 
  if ~isfield(model,'DOF');model.DOF=DOF;
  elseif ~isequal(model.DOF,DOF); error('DOF mismatch');
  end
  if ~isequal(def.DOF,DOF);def=feutilb('placeindof',DOF,def);end
  fe_mknl('assemble',model,Case,300,def.def); 
  [EGroup,nGroup,names]=getegroup(model.Elt);
  for jGroup=find(strcmpi(names,'beam1t'))'
   if size(Case.GroupInfo{jGroup,3},2)>1;
       error('Single ProId per group expected');
   end
   MAP=Case.GroupInfo{jGroup,7};
   model=feutil(sprintf('setpro %i',Case.GroupInfo{jGroup,3}(2)), ...
       model,'MAP',MAP);
  end
  out=model; return;

%% #ViewTen : display tensions as color - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'viewten');[CAM,Cam]=comstr(CAM,5); 

 carg=1;
 model=varargin{carg};carg=carg+1;
 if ~isfield(varargin{carg},'T')
  Case=fe_mknl('init',model);
 else;Case=varargin{carg};carg=carg+1;
 end
 def=varargin{carg};carg=carg+1; 
 if ~isequal(model.DOF,def.DOF);error('Mismatch');end 
 eltid=feutil('eltidfix;',model);[EGroup,nGroup]=getegroup(model.Elt);
 out=struct('X',{{eltid*0}},'Xlab',{{'EltId','def'}}, ...
     'Y',zeros(length(eltid),size(def.def,2)));

 for j1=1:size(def.def,2)
  fe_mknl('assemble',model,Case,300,def.def(:,j1)); 
  for jGroup=1:size(Case.GroupInfo,1)
   cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
   out.Y(cEGI,j1)=Case.GroupInfo{jGroup,7}.data(11,:);
   out.X{1}(cEGI)=eltid(cEGI);
  end
 end
 i1=out.X{1}~=0;out.X{1}=out.X{1}(i1);out.Y=out.Y(i1,:);
 if nargout==0; fecom('colordataelt',out);clear out;end

%% #Test - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -----
 elseif comstr(Cam,'test');[CAM,Cam]=comstr(CAM,5); 

 if comstr(Cam,'mat')
   model=femesh('testbeam1');model.Elt(1,1:7)=[Inf abs('beam1t')];
   model=fe_mknl(model,'NoT');k=model.K([2 1]);
   out=model.K([2 1]);  
   return;  
 elseif comstr(Cam,'load')
   [out,out1]=femesh('teststruct beam1 load');
   out.Elt=feutil('setgroupall name beam1t',out);
   out1=fe_load(out);
   return;
 end
 [CAM,Cam,i1]=comstr('divide',[-25 1],CAM,Cam);
 if isempty(i1); i1=50;end
 model=femesh(sprintf('testbeam1 divide %i',i1));
  model.Elt=feutil('set group 1 name  beam1t',model);

  model=fe_case(model, 'reset', ...
   'fixdof','pinned',[1+[.01;.02;.03];2+[.01;.02;.03]], ...
   'fixdof','2D',[.01; .02;.04;.06]);
  model.il=[112 fe_mat('p_beam','SI',1) [1e-12 1e-12 1e-12] .1];
  model.pl(1,5)=1000; % 1000 kg/m3

  if comstr(Cam,'eig')
    def=fe_eig(model,[5 10]);
    if ~isempty(strfind(Cam,'stress'))
        [C1,model.DOF]=fe_mknl('init',model); 
        def=feutilb('placeindof',model.DOF,def);
        C1=fe_mknl('gstate-struct',model,C1,def);
        k=fe_mknl('assemble',model,C1,fe_def('subdef',def,1),1);
        %out1=fe_stress('stressgstate',model,def)
        r1=C1.GroupInfo{5};
        figure(1);plot(squeeze(r1.Y(2,:,:))')
    else;out1=def; 
    end
    out=model;
    return;
  end
  model=feutil('setpro 112',model,'MAP', ...
    struct('dir',{{'1e8'}},'lab',{{'ten'}}));

  if nargout==1; out=model; return; end
  
  [Case,model.DOF]=fe_mknl('init',model);
  %T=1e8; Case.GroupInfo{1,5}(11,:)=T;

  k=fe_mknl('assemble',model,Case,1);
  m=fe_mknl('assemble',model,Case,2);

  g=(m*sum(fe_c(Case.DOF,.03))');
  kd=ofact(k); q=kd\g;ofact('clear',kd)

  i1=feutil('findnode groupall',model);
  x=sort(model.Node(i1,5));
  i2=fe_c(Case.DOF,.03,'ind');
  x=model.Node(fix(Case.DOF(i2)),5);
  q(i2,2)=model.il(1,6)*model.pl(1,5)*1/T*(x.*(1-x))/2;

  figure(1);plot(x,q(i2,1)./q(i2,2))
  ylabel('Computed/true deflection');xlabel('x(m)');

  model.il=p_beam('dbval 112 circle 1e-2');

  [Case,model.DOF]=fe_mknl('init',model);
  T=1e3; Case.GroupInfo{1,5}(11,:)=T;
  k=fe_mknl('assemble',model,Case,1);
  f=fe_c(Case.DOF,28.03)'*100;
  kd=ofact(k); q=kd\f;ofact('clear',kd)
  r1=q(i2,1);

  model=fe_case(model,'fixdof','2D',[.01;.03;.04;.05]);
  [Case,model.DOF]=fe_mknl('init',model);
  T=1e3; Case.GroupInfo{1,5}(11,:)=T;
  k=fe_mknl('assemble',model,Case,1);
  f=fe_c(Case.DOF,28.02)'*100;
  kd=ofact(k); q=kd\f;ofact('clear',kd)
  r1(:,2)=q(i2,1);
  figure(102);h=plot(x,r1);set(h(2),'linestyle',':','marker','x');

 if 1==2 % arch test to correct normal computations - - - - - - -

 FEnode=[1 0 0 0  -1 -1 0];
 femesh('object mass 1'); 
 femesh('rev 20 o 0 0 0 180 -1.0 1.0 0')
 femesh plotel0
 FEel0(2:end,3:4)=1; FEel0(2:end,7)=1;
 model=femesh('model0');
 MAP=beam1t('map',model);

 end% arch test to correct normal computations - - - - - - -
 
 elseif comstr(CAM,'cvs')
  out='$Revision: 1.77 $  $Date: 2017/05/02 14:52:16 $'; return;
 elseif comstr(Cam,'mexon');
   if exist('getBeamK_mex','file');getKFcn=@getBeamK_mex;
   else; getKFcn=@getBeamK;
   end
 elseif comstr(Cam,'mexoff');getKFcn=@getBeamK;
 elseif comstr(Cam,'vof');vof=comstr(CAM(4:end),-1);
 elseif comstr(Cam,'@');out=eval(CAM);
 else;out=beam1(CAM);
 end

return
end % of standard calls with one input argument

FromOf_Mk=0;
%% #Assemble element matrix assembly - - - - - - - - - - - - - - - - - - ----
if nargin==10
 point=CAM;pointers=point;FromOf_Mk=1;
 %[integ,constit,gstate,elmap,InfoAtNode,EC,def,jElt,DofPos]=deal(varargin{:});
 integ=varargin{1};constit=varargin{2};gstate=varargin{3};
 elmap=varargin{4};InfoAtNode=varargin{5};EC=varargin{6};
 def=varargin{7};jElt=varargin{8};DofPos=varargin{9};
 constit=constit+0;
 nElt=jElt;
 state=EC.nodeE(1,5:15)';
 elt=EC.eltg(:,jElt(1)+1)';
 %elt=evalin('caller',sprintf('model.Elt(cEGI(%i),:)',jElt(1)+1));
 if ~isempty(EC.CTable)
  of_time('cinterp',EC.nodeE,[.5;.5],constit,EC.CTable,zeros(1));
  if sp_util('diag')>10
      fprintf('J %.1e I1 %.1e I2 %.1e A %.1e\n',constit(8:11));
  end
 end
 
else
 nodeE=CAM; 
 elt=varargin{1}; 
 point=varargin{2};point=double(point);
 integ=varargin{3};
 constit=varargin{4}; % sdtweb p_solid('constitbeam')
 gstate=varargin{5}; 
 state=varargin{6};EC=varargin{7};

end
typ=int32(point(5));

% basis function - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if any([300 302]==typ) 
 %% #stateUpdate - - - - - - - - - - - - - - - - - -
 % Only performed outside of of_mk
 %  fe_mknl('assemble',model,Case,300,def.def); % updates InfoAtNode T
 EC=evalin('caller','EltConst'); if ~isfield(EC,'CTable');EC.CTable=[];end
 InfoAtNode=varargin{6};
 FieldAtEltDof=varargin{7}; 
 jElt=varargin{8};
 if jElt>0
  state=InfoAtNode.data(:,jElt);L=state(10); 
  if ~isempty(EC.CTable) % 'xxx need optimize of nodeE'
   constit1=constit+0;
   nodeE=[nodeE(:,1:4) ...
      InfoAtNode.data(:,InfoAtNode.NodePos(:,evalin('caller','jElt')))'];
   of_time('cinterp',nodeE,[.5;.5],constit1,EC.CTable,zeros(1));
   constit=constit1;
  end
  cLG=reshape(state(1:9),3,3);
  r1=constit(point(7)+1)*constit(point(7)+11)/L; % EA/L
  r1=(-cLG(1,:)*FieldAtEltDof(1:3)+cLG(1,:)*FieldAtEltDof(7:9))*r1;
  if typ==302 % Possibly update orientation
   x=diff(nodeE(:,1:3)+FieldAtEltDof([1 2 3;7 8 9]));L=norm(x);
   cLG=sp_util('basis',x,cLG(2,:))'; % Should be incremental rotation xxx
  end
  r1=[cLG(:);L;r1];
  if isfield(InfoAtNode,'data0');r1=r1+InfoAtNode.data0(:,jElt);end
  sp_util('setinput',InfoAtNode.data,r1,int32(1:11),int32(jElt));
 else % Actually do the loop here
   cEGI=evalin('caller','cEGI');NodePos=evalin('caller','NodePos');
   Case=evalin('caller','Case');jGroup=evalin('caller','jGroup');
   FieldAtDof=evalin('caller','FieldAtDof');
   for jElt=1:length(cEGI)
      i1=NodePos(:,jElt);
      nodeE=Case.Node(i1,[5:7 1]);
      FieldAtEltDof=FieldAtDof(double(Case.GroupInfo{jGroup,1}(:,jElt))+1);
      state=InfoAtNode.data(:,jElt);L=state(10); 
      cLG=reshape(state(1:9),3,3);
      if ~isempty(EC.CTable); % tested in OSCAR
       constit1=constit+0;
       nodeE=[nodeE(:,1:4) InfoAtNode.data(:,InfoAtNode.NodePos(:,jElt))'];
       of_time('cinterp',nodeE,[.5;.5],constit1,EC.CTable,zeros(1));
       r1=constit1(point(7)+1)*constit1(point(7)+11)/L;
      else;      r1=constit(point(7)+1)*constit(point(7)+11)/L; % EA/L
      end
      r1=(-cLG(1,:)*FieldAtEltDof(1:3)+cLG(1,:)*FieldAtEltDof(7:9))*r1;
      if typ==302 % Possibly update orientation
       x=diff(nodeE(:,1:3)+FieldAtEltDof([1 2 3;7 8 9]));L=norm(x);
       cLG=sp_util('basis',x,cLG(2,:))'; % Should be incremental rotation xxx
      end
      r1=[cLG(:);L;r1];
      if isfield(InfoAtNode,'data0');r1=r1+InfoAtNode.data0(:,jElt);end
      sp_util('setinput',InfoAtNode.data,r1,int32(1:11),int32(jElt));
   end
  assignin('caller','jElt',[]);  return
 end
 out=[];
 return; 
elseif isstruct(gstate) 
  %% #Stress computation
  defe=evalin('caller','full(def.def(DofPos(:,jElt)+1,:))');
  %fe_c(model.DOF(DofPos(:,jElt)+1))
  L=state(10); bas=reshape(state(1:9),3,3)/L^2; Tpin=[];
  if strcmp(gstate.Xlab{end-1},'DOF');defe=eye(12);end
  defe(1:3,:)=bas*defe(1:3,:); % local coordinates
  defe(7:9,:)=bas*defe(7:9,:);
  bas=L*bas;
  defe(4:6,:)=bas*defe(4:6,:);
  defe(10:12,:)=bas*defe(10:12,:);
  
  r=EC.w(:,1);
  %state(11)*L*wx*defe([2 6 8 12],:) % pre-stress
  if typ==2 % Mass strains (displacement)
      r2=[[1-r r]*defe([1 7],:) % disp z
      L^2*EC.N*defe([2 6 8 12],:) % disp y
      L^2*EC.N*diag([1 -1 1 -1])*defe([3 5 9 11],:) % disp z
      [1-r r]*defe([4 10],:) % disp z
      L*EC.Nr*defe([2 6 8 12],:) % bending y
      L*EC.Nr*diag([1 -1 1 -1])*defe([3 5 9 11],:);% bend z, theta_y = - dz/dx
      ];
     r2=reshape(r2,EC.Nw,6,size(defe,2));
  elseif typ==1 % Stiffness sdtweb p_solid constitbeam (Timoshenko)
      dd=buildDD(constit,point(7));
      r2=[ones(size(r))*[-L L]*defe([1 7],:) % u,x % missing L May 11
      ones(size(r))/L*[12  6 -12  6]*defe([2 6 8 12],:) % shear s_y = y'''
      ones(size(r))/L*[12 -6 -12 -6]*defe([3 5 9 11],:) % shear s_z
      ones(size(r))*[-1 1]*defe([4 10],:) % phi,x
      EC.Nrr*defe([2 6 8 12],:) % bending y
      EC.Nrr*diag([1 -1 1 -1])*defe([3 5 9 11],:);% bend z, theta_y = - dz/dx
      ];
     r2=reshape(r2,EC.Nw,6,size(defe,2));
     for j1=1:6;r2(:,j1,:)=r2(:,j1,:)*dd(j1);end
  end
  r2=permute(r2,[2 1 3]);%squeeze(r2(4,:,:))
  
  if isfield(EC,'StressOut')&&typ==1
      r2=EC.StressOut.cta*reshape(r2,6,[]);
  else; r2=reshape(r2,size(gstate.X{1},1),[]);
  end
  
% z=evalin('caller','full(def.def(DofPos(:,jElt)+1,:))');
%disp(reshape(r2(4,:),size(gstate.X{2},1),12)*z)
%r3=integrules('beam1',linspace(0,1,30)'*[1 0 0 0]);
%figure(1);subplot(211);
%plot(r3.w(:,1),r3.N*diag([1 -1 1 -1])*defe([3 5 9 11],:),[0 1],defe([3 9],:),'x')
%subplot(212);plot(r3.w(:,1),r3.N*defe([2 6 8 12],:),[0 1],defe([2 8],:),'x')
%http://en.wikiversity.org/wiki/Nonlinear_finite_elements/Euler_Bernoulli_beams
  if isfield(gstate,'Y') % Return axial strain
   jElt=evalin('caller','jElt');i2=size(gstate.Y);i2(end)=[];
   sp_util('setinput',gstate.Y,r2,(jElt-1)*prod(i2));
  end
elseif ~any([0 1 2 3 5 8 100 102 103]==typ)
 warning('Matrix type %i not supported by beam1t',typ);
end

% stiffness - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

L=state(10); cLG=reshape(state(1:9),3,3); 
Tpin=[-1e100 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0];

if size(elt,2)>9&&any(elt(9:10)) % Condense pin flag
    i1=abs(sprintf('%i',elt(1,9)))-48; i2=abs(sprintf('%i',elt(1,10)))-48;
    i1=[i1(i1~=0) i2(i2~=0)+6];
    if ~isempty(i1)
      s2=state;s2(1:9)=[1 0 0 0 1 0 0 0 1];
      k=getKFcn(constit(double(point(7))+(1:size(constit,1))), ...
        EC.defe,EC.w,EC.Nr,EC.Nrr,s2,typ,Tpin);%local basis k
      i2=1:12;i2(i1)=0;i2=find(i2);
      Tpin=zeros(12); Tpin([i2(:);i1(:)],i2)= ...
        [eye(length(i2),length(i2));-k(i1,i1)\k(i1,i2)];
    end   
end

if any([0 1 5]==typ)
   k=getKFcn(constit(double(point(7))+(1:size(constit,1))), ...
        EC.defe,EC.w,EC.Nr,EC.Nrr,state,typ,Tpin);
     % pin flag handling (condense DOF a put zero stiffness)
   if typ==5;
     %if any(elt(1:2)==355)
     % h=evalin('base','h');
     % h=[h;diff(EC.nodeE(:,1:3))/norm(diff(EC.nodeE(:,1:3)))  cLG(1,:) state(11)];
     % assignin('base','h',h);
     %end
     sp_util('setinput',EC.Be,k*EC.defe+ ...
       [cLG(1,:)*-state(11) zeros(1,3) cLG(1,:)*state(11) zeros(1,3)]',zeros(1)); 
   end
   out=k; out1=[];
end

if any([0 2]==typ) % mass - - - - - - - - - - - - - - - - - - - - -


   rho=constit(point(7)+3);
   ie=constit(double(point(7))+(8:11));
   %A=constit(point(7)+11);
   m = zeros(12,12);

    % consistent mass
    %ind=[1 7];m=zeros(12);m(ind,ind)=ones(2);sprintf('% i',find(m))
    r1=rho*ie(4)*L/3; % Rho * A * L
    m([1 7 73 79]) = [1 .5 .5 1]*r1;

    r1=rho*sum(ie(2:3))*L/3; % Rho * (i1+i2) * L
    m([40 46 112 118]) = [1 .5 .5 1]*r1;

    L2=L^2; r1=ie(2); %I1
    if r1; ind = [2 6 8 12];
      r1=rho*ie(4)*L/210;  % rho A L /210
       m(ind,ind) = r1 * ...
	[78    11*L   27   -6.5*L; 11*L    2*L2    6.5*L -1.5*L2
	 27    6.5*L  78   -11*L ; -6.5*L -1.5*L2 -11*L   2*L2];
    end
    r1=ie(3); %I2 constit(point(7)+10);
    if r1; ind = [3 5 9 11];
      r1=rho*ie(4)*L/210; % rho A L /210
      m(ind,ind) =  r1 * ...
	[78   -11*L   27    6.5*L;-11*L    2*L2   -6.5*L -1.5*L2
	 27   -6.5*L  78    11*L ; 6.5*L  -1.5*L2  11*L   2*L2];
    end
    if Tpin(1)~=-1e100; m=Tpin'*m*Tpin; end

   % Lumped mass (sum over elements)
   if length(constit)>=point(7)+14 &&constit(point(7)+14)==1
    m=of_mk('xkx_trans',cLG,m); if off; m = tr'*m*tr;end % coordinate transform
    %r2=[0 0 0 1 0 0   0 0 0 1 0 0;
    %    0 0 0 0 1 0  0 0 -L 0 1 0;0 0 0 0 0 1   0 L 0 0 0 1];
    r2=[zeros(3);cLG';cLG'*[0 0 0;0 0 l;0 -l 0];cLG'];
    r3=sum(r2'*m*r2)/2;if any(r3)<0; r3=diag(r2'*m*r2)'/2;end
    r1=[(m(1)+m(7))*[1 1 1] r3 (m(1)+m(7))*[1 1 1] r3]; 
    m=diag(r1);
   else
     m=of_mk('xkx_trans',cLG,m); % coordinate transformation
   end

   if typ==2; out=m; out1=[];  
   else; out1=m; end
end

if typ==3; out=zeros(12,12); out1=[]; end 

if any([100 102]==typ) 
%% #Volume_load - - - - - - - - - - - - - - - - - -

  if typ==102; % multiply by section inertia  
    %rho=constit(point(7)+3);
    %A=constit(point(7)+11);
    dire=reshape(varargin{8},[],2);dire=dire(1:3,:); % global coord
    dire=dire*constit(point(7)+3); % A is taken into account in the loop
  else;
    dire=varargin{8}; % element coordinates as colums with 3 components
  end 
  dxds=L;
  Nu=[1-EC.w(:,1) EC.w(:,1)];
  Nw=EC.N; Nw(:,[2 4])=Nw(:,[2 4])*L; % rotation DOFs

  r1=cLG(1,:)'*cLG(1,:); r2=cLG(2,:)'*cLG(2,:);  r3=cLG(3,:)'*cLG(3,:);

  b1=zeros(12,1);
  for j1=1:size(Nu,1) % loop on integration points

    d=zeros(12,3);
    d(1:3,:)=r1*Nu(j1,1);% node 1 translation
    d(7:9,:)=r1*Nu(j1,2);% node 2 translation
    % Work associated with transverse displacement
    d(1:3,:)=d(1:3,:)+r2*Nw(j1,1); d(4:6,:)  =r2*Nw(j1,3);
    d(7:9,:)=d(7:9,:)+r2*Nw(j1,2); d(10:12,:)=r2*Nw(j1,4);

    d(1:3,:)=d(1:3,:)+r3*Nw(j1,1); d(4:6,:)  =r3*Nw(j1,3);
    d(7:9,:)=d(7:9,:)+r3*Nw(j1,2); d(10:12,:)=r3*Nw(j1,4);

    F = dire*Nu(j1,1:2)';
    b1=b1+d*(F*dxds*constit(point(7)+11)*EC.w(j1,4));%A
  end % loop on integration points

  out=b1;
elseif any(103==typ) 
%% Initial thermal load - - - - - - - - - - - - - - - -
    
 InfoAtNode=varargin{6};
 %FieldAtEltDof=varargin{7}; 
 r1=constit(point(7)+[6 7 1 11]); % alpha T0
 % E A alpha Dt
 out=[-cLG(1,:),0,0,0,cLG(1,:),0,0,0]'*r1(1)*r1(3)*r1(4)* ...
     (InfoAtNode(12)-r1(2));
elseif any(8==typ) 
%% Centrifugal softening, sdtweb t_beam('ng')
  %rho=constit(point(7)+3);
  %A=constit(point(7)+11); 
  dd=EC.ConstitTopology{8};
  dd=constit(dd); % rho A Omega^2

  L=state(10); bas=reshape(state(1:9),3,3); 
  k=zeros(12); J=L;
  for jw=1:size(EC.w,1)
    r=EC.w(jw,1);N=EC.N(jw,:).*[1 L 1 L];
    Ng=bas'*[ ... % shape functions with global input
        [1-r 0 0;0  N(1) 0  ;   0 0       N(1)  ]*bas ... % tA(local)=bas*tA
        [0 0 0  ;0  0    N(2);0 -N(2) 0]*bas ...  % rA
        [r 0 0;  0  N(3) 0;     0 0       N(3)  ]*bas ...    % tB
        [0 0 0;  0  0    N(4);0 -N(4) 0]*bas];    % rB
   k=k+Ng'*(dd*J*EC.w(jw,4))*Ng;
  end
  out=k;out1=[];

elseif FromOf_Mk&&typ<100
  sp_util('setinput',EC.ke,out,zeros(1));out=EC;
end

function dd=buildDD(constit,p7)


dd=[constit(p7+1)*constit(p7+11); %EA
    constit(p7+4)*constit(p7+11) % GA
    constit(p7+4)*constit(p7+11) % GA
    constit(p7+4)*constit(p7+8) % GJ
    constit(p7+1)*constit(p7+9); % EI 
    constit(p7+1)*constit(p7+10); % EI 
    ];

%% k=getBeamK(constit(double(point(7))+(1:size(constit,1))),L,cLG,EC.defe,EC.w,EC.Nr,EC.Nrr,state)

function  k=getBeamK(constit,defe,ECw,ECNr,ECNrr,state,typ,Tpin)

%coder.varsize('constit',[15 1],[0 0]);
%coder.varsize('defe',[12 1],[0 0]);
%coder.varsize('ECw',[5 4],[1 0]);
%coder.varsize('state',[11 1]);
%coder.varsize('typ',[1 1]);

   L=state(10); cLG=reshape(state(1:9),3,3);
   k = zeros(12,12);
   r1=constit(1)*constit(11)/L; %EA/L
   ind = [1 7 73 79];     k(ind) = [1 -1 -1 1]*r1; % 1 7
   % 4 10, torsionnal stiffness GJ/l
   r1=constit(4)*constit(8)/L; %GJ/L
   ind = [40 46 112 118]; k(ind) = [1 -1 -1 1]*r1; % 4 10

   wx=ECNr;
   if size(wx,1)==3;
    wx([1 2 3 7 8 9])=wx([1 2 3 7 8 9])/L;
    wxx=ECNrr/L;wxx([1 2 3 7 8 9])=wxx([1 2 3 7 8 9])/L;
   else
     wx(:,[1 3])=wx(:,[1 3])/L;
     wxx=ECNrr/L;wxx(:,[1 3])=wxx(:,[1 3])/L;
   end
   if typ==5;r1=(state(11)+cLG(1,:)*(defe(1:3)-defe(7:9)))*L;
   else;r1=state(11)*L;                               % T  * L
   end
   r2=constit(1)*constit(9)*L; % EI * L
   dk=zeros(4);
   for j1=1:size(ECw,1)
    dk=dk+ECw(j1,4)*r1*wx(j1,:)'*wx(j1,:)+  ... % Pre-tension stiffness
       ECw(j1,4)*r2*wxx(j1,:)'*wxx(j1,:);       % bending stiffness
   end
   k([14 18 20 24 62 66 68 72 86 90 92 96 134 138 140 144])=dk;%[2 6 8 12];

   if size(wx,1)==3;
     wx([1 2 3 7 8 9])=-wx([1 2 3 7 8 9]);  % theta_y = - dz/dx
     wxx([1 2 3 7 8 9])=-wxx([1 2 3 7 8 9]);
   else;
      wx(:,[1 3])=-wx(:,[1 3]);  % theta_y = - dz/dx
      wxx(:,[1 3])=-wxx(:,[1 3]);
   end
   r2=constit(1)*constit(10)*L; % EI * L
   dk=zeros(4);
   for j1=1:size(ECw,1)
    dk=dk+ECw(j1,4)*r1*wx(j1,:)'*wx(j1,:)+  ... % Pre-tension stiffness
       ECw(j1,4)*r2*wxx(j1,:)'*wxx(j1,:);       % bending stiffness
   end
   k([27 29 33 35 51 53 57 59 99 101 105 107 123 125 129 131])=dk;%[3 5 9 11];

   if Tpin(1)~=-1e100;
       k=Tpin'*k*Tpin;
   end
   k=of_mk('xkx_trans',cLG,k); % coordinate transformation
   %for j1=1:4;
   %  k(:,j1*3+(-2:0))=k(:,j1*3+(-2:0))*cLG;
   %end
   %for j1=1:4;
   %  k(j1*3+(-2:0),:)=cLG'*k(j1*3+(-2:0),:);
   %end
