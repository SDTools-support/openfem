function [out,out1,out2,out3]=p_zt(varargin)

%P_ZT zero thickness property
%
%  Material property defined using 2D elastic material in local coordinates
%    directions are tangent1, tangent2, normal. Thus a typically isotropic
%    in plane behavior would be given by
%       [103 fe_mat('m_elastic','SI',4)  k_t1  0 k_t2 0 0 k_n  RhoS(9) eta(10) ...
%           a1 a2 a3   c_t1(15) 0 c_t2 0 0 c_n(19)]
%       m_elastic('propertyunittypecell',4)
%
%       See sdtweb      fem (handling materials section), pl, fe_mat, p_shell
%       See also help   fe_mat


%       Etienne Balmes, with discussion with Phuor Ty
%       Copyright (c) 2001-2022 by SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*NOSEM>

if nargin<1; help p_solid;return; end
if ischar(varargin{1}); [CAM,Cam]=comstr(varargin{1},1); il=[]; carg=2;
else;il=varargin{1};[CAM,Cam]=comstr(varargin{2},1); carg=3;
end

% -------------------------------------------------------------------------
if comstr(Cam,'default') % #Default -1
 
 model=[];if carg<=nargin;model=varargin{carg};carg=carg+1;end
 if isempty(model); out=p_solid('database');out=out(1);
 else;              out=fe_mat(['defaultil' CAM(8:end)],model); % sdtweb fe_mat('default')
 end
 
%% #Info ---------------------------------------------------------------------
elseif comstr(Cam,'info')

 r1=p_solid('database');fprintf('%s\n',r1.name);

%% #Dbval --------------------------------------------------------------------
elseif comstr(Cam,'dbval')

 while 1==1
  i1=strfind(comstr(Cam,-27),'-unit'); out1={};
  if ~isempty(i1)
   [Unit,i2,i3,i4]=sscanf(CAM(i1+5:end),'%s',1);
   i4=i1+(0:4+i4);CAM(i4)=''; [CAM,Cam]=comstr(CAM,1);
  else;Unit='';
  end
  i2=strfind(comstr(Cam,-27),'-punit');
  if ~isempty(i2)
   [PUnit,i3,i4,i5]=sscanf(CAM(i2+6:end),'%s',1);
   i5=i2+[0:5+i5];CAM(i5)=''; [CAM,Cam]=comstr(CAM,1);
  else;PUnit='';
  end
  
  if ischar(CAM); [i1,CAM,Cam]=comstr(CAM,'dbval','%i');else; i1=[];end
  if isempty(CAM)&&carg<=nargin; st=varargin{carg};carg=carg+1;
  else;st=CAM;end
  if isempty(st);
  elseif ischar(st); mat=p_solid('database',st);
  elseif isnumeric(st)
   [typ,st1,i4]=fe_mat('typep',st(2));
   mat=struct('il',st,'name',sprintf('%i',st(1)),'type',typ,'unit',st1);
  end
  if ~isempty(PUnit)
   r1=fe_mat(sprintf('convert %s %s',mat.unit,PUnit),mat.il(1:2));
   mat.il(2)=r1(2); mat.unit=PUnit;
  end
  if ~isempty(Unit)
   mat.il=fe_mat(sprintf('convert %s %s',mat.unit,Unit),mat.il);mat.unit=Unit;
  end
  r1=mat.il; if length(i1)==1; r1(1)=i1;end
  if ~isempty(il); i2=find(il(:,1)==r1(1)); else;i2=[];end
  if isempty(i2); i2=size(il,1)+1;end
  il(i2,1:length(r1))=r1; %#ok<AGROW>
  if carg>nargin; break;end
  CAM=varargin{carg};carg=carg+1;if ischar(CAM);[CAM,Cam]=comstr(CAM,1);end
 end
 out=il;


%% #DataBase -----------------------------------------------------------------
elseif comstr(Cam,'database') 

  st=comstr(CAM,9);
  if isempty(st)&&carg<=nargin; st=varargin{carg}; carg=carg+1;end
  [MatId,i2,i3,i4]=sscanf(st,'%i',1); if i2~=1; MatId=1;end
  st=comstr(st,i4);

  % MatId Typ                     CordM Integ Stres Isop Field
  out.il=[MatId fe_mat('p_solid','SI',1) 0 -3 ]; 
  out.name='Default for topology';
  out.type='p_solid';
  out.unit='SI';
  out1='0Thick';


%% #PropertyUnitType  --------------------------------------------------------
elseif comstr(Cam,'propertyunittype')

 if nargin==1;out=1; return; end % return subtypes ID
 i1=varargin{carg};
 out1={};
 switch i1 % PropertySubType
 case 1 % [ProId type Coordm In Stress Isop Fctn  ]
   st={ ...
   'ProId'   0  'sdtweb(''p_solid'')';
   'Type'    0  '';
   'COORDM'  0  'Coordinates system id';
   'IN'      0  'Integration rule'};
 otherwise; st={'ProId' 0 'sdtweb(''p_zt'')'; 'Type', 0, ''};
 end
 if ~isempty(strfind(Cam,'cell')); out=st; else; out=[st{:,2}]; end

%% #SubTypeString ------------------------------------------------------------
elseif comstr(Cam,'subtypestring')

 i1=varargin{carg}; carg=carg+1;
 switch i1
 case 1;  out='p_solid';
 otherwise; out='p_solid';
 end

%% #const : defines EltConst based on model content --------------------------
% [EltConst,NDNDim]=p_solid('Const',ElemF,integ,constit,model,Case,cEGI,RunOpt);
elseif comstr(Cam,'const')

 EC=varargin{carg};carg=carg+1;
 integ=varargin{carg};carg=carg+1;
 if carg<=nargin; constit=varargin{carg};carg=carg+1;else;constit=[];end
 if carg<=nargin; model=varargin{carg};carg=carg+1;else;model=[];end
 if carg<=nargin; Case=varargin{carg};carg=carg+1;else;Case=[];end
 if carg<=nargin; cEGI=varargin{carg};carg=carg+1;else;cEGI=[];end
 if carg<=nargin; RunOpt=varargin{carg};carg=carg+1;else;RunOpt=[];end
 
 r1={3,'t3p';4,'q4p';6,'t6p';8,'q8p'};
 r1=r1{vertcat(r1{:,1})==length(feval(EC,'node'))/2,2};
 EC=integrules(r1,integ(5,1)); % Flat shape functions
 % Recombine in two layers of flat functions
 z=zeros(size(EC.N)); EC.N=[EC.N z;z EC.N]; % 
 EC.Nr=[EC.Nr z;z EC.Nr];EC.Ns=[EC.Ns z;z EC.Ns];
 EC.xi=[EC.xi;EC.xi]; EC.w=[EC.w;EC.w]; % Replicate rule for each layer
 EC.NDN=zeros(size(EC.N,2),3*size(EC.w,1)); 
 EC.jdet=zeros(size(EC.w,1),1); 
 EC.type=sprintf('zt%i',EC.Nw);% EC.Nw=half to allow duplication
 EC.Nnode=size(EC.N,2);
 
 EC.bas=zeros(9,size(EC.w,1)); % Keep local coordinates
 EC.nodeE=zeros(size(EC.N,2),4); 
 EC.CTable=[1 [7 0 0 -2 0 0 0]]; % sdtweb m_elastic ctable
 out2=[]; if ~isempty(Case);out2=Case.GroupInfo{Case.jGroup,7};end

 if size(integ,1)<6||any(integ(5:6,1)==0)||integ(6)>size(EC.w,1);
       rule=[1 EC.Nw]; % By default start at 1 use all points
 else; rule=integ(5:6,1);rule=rule(:)';
 end
 if ~any(rule); 
   rule=[1 EC.Nw];
   sdtw('integ is assumed to give standard integration rule in integ(5:6)');
 elseif rule(1)==0||rule(end)>size(EC.w,1);error('Inconsistent rule'); 
 end
 
 if constit(1)==-1 % Default 
     EC.DofLabels={'u','v','w'};
     EC.material='linear'; EC.bas=zeros(9,size(EC.N,1));
     r2=rule+[rule(end) 0]; % Second layer shape functions
     EC.StrainDefinition{1}= ... % row,NDNBloc,DOF,NwStart,NwTot
      [1 -1 1 rule; 2 -1 2 rule;3 -1 3 rule;  % disp first layer
       1  1 1 r2;   2  1 2 r2;3  1 3 r2  ];
     EC.StrainLabels(1:2)={{'t1','t2','g'},{'um','vm','wm'}};
     %i1=reshape(1:9,3,3);i1=i1([3 1 2],[3 1 2]);comstr(i1(:)',-30)
     %EC.ConstitTopology{1}=int32(reshape(6+[9,7,8,3,1,2,6,4,5],3,3));
     EC.ConstitTopology{1}=int32(reshape(6+(1:9),3,3));
     if length(constit)==24
      EC.StrainDefinition{3}=EC.StrainDefinition{1};
      EC.StrainLabels{3}={'t1','t2','g'};
      EC.CTable=[2 [7 0 0 -2 0 0 0] [16 0 0 -2 0 0 0]];
      EC.ConstitTopology{3}=int32(reshape(15+(1:9),3,3));
     end
     
     EC.StrainDefinition{2}= ... % row,NDNBloc,DOF,NwStart,NwTot
      [1  1 1 rule; 2  1 2 rule; 3  1 3 rule; % disp first layer
       1  1 1 r2;   2  1 2 r2;   3  1 3 r2];
     EC.ConstitTopology{2}=int32(eye(3)*5);

     Ndof=3*EC.Nnode; % VectMap translations
     EC.VectMap=int32([1:3:Ndof 2:3:Ndof 3:3:Ndof])';
     out1=23*ones(1);% Use 2D surface rule
     
 else; error('Not implemented');
 end
 EC=integrules('matrixrule',EC);
 EC.ConstitTopology{5}=[]; EC.MatrixIntegrationRule{5}=[];
 EC=feutil('rmfield',EC,'material'); % Generic does not use material
 
 % integ(1)==-1 is done at the very beginning of the command
 
 if nargout==0
  integrules('texstrain',EC); % Does not work here due to repetition
  try; EC=integrules('stressrule',EC);integrules('texstress',EC);end
 else; out=EC;
 end


% -------------------------------------------------------------------------
% #BuildConstit Implementation of elastic constitutive law building 3D
%[constit,integ,Inits,Data]=p_solid('buildconstit Nfield Nnode',ID,pl,il, ...
%    model,Case); This is called by the element 'integinfo' command
elseif comstr(Cam,'buildconstit');

  RunOpt=struct('warn',{{}},'Dim',0,'ProStack',[],'DoElmap',1);
  try;RunOpt.Dim=comstr(Cam(13:end),-1); end
  if isempty(RunOpt.Dim); RunOpt.Dim=[3 8]; end
   
  ID=varargin{carg};carg=carg+1;out1=int32(ID);out3=struct;
  pl=varargin{carg};carg=carg+1;   
  if isempty(pl);mat=[];else;mat=pl(pl(:,1)==ID(1),:);end
  il=varargin{carg};carg=carg+1; 
  if isempty(il);pro=[];else;pro=il(il(:,1)==ID(2),:);end
  if isempty(mat); 
   RunOpt.warn{end+1}=sprintf('MatId %i is not matched',ID(1));
   st='m_null';unit=1;typ=1;
  else;  [st,unit,typ]=fe_mat('typem',mat(2));
  end
  RunOpt.MatT=fe_mat('typemstring',mat(2));
  RunOpt.ProT=fe_mat('typepstring',pro(2));
  
  %% #p_zt.2  - - - - - - - - - - - - - - - - - - - - - - - - - -
  if strcmp(RunOpt.ProT,'p_zt.1')&&strcmp(RunOpt.MatT,'m_elastic.4')
   out3.ConstitLab={'-1','il(2)','pl(1)','pl(2)','RhoS/4','eta', ...
        'D11','D21','D31','D12','D22','D32','D13','D23','D33'};
   ID(3:4)=RunOpt.Dim(2)*[3 1]; % NDOF,NNode
   ID(5)=pro(3);  % Integrule
   if length(mat)<10;mat(10)=0;end
   constit=[-1;pro(2);mat(1);mat(2);mat(9)/4;mat([10 3 4 6 4 5 7 6 7 8])']; % p_solid will call p_zt
   if length(mat)>15
    constit=[constit;mat([15;16;18;16;17;19;18;19;20])'];
    out3.cc=reshape(constit(16:24),3,3);
   end
   out3.dd=reshape(constit(7:15),3,3);
  else; error('%s Not implemented',RunOpt.ProT)
  end
  out=constit(:);
  out1=int32(ID(:));
  out2=elem0('elmapmat_og',[double(out1(4)) double(out1(3))/double(out1(4))]);
  if ~isempty(RunOpt.warn);sdtw('_nb','%s\n',RunOpt.warn{:});end

% -------------------------------------------------------------------------
%% #BuildDof comming from p_solid
elseif comstr(Cam,'builddof')

  RunOpt=varargin{carg};carg=carg+1;
  pl=varargin{carg};carg=carg+1;
  il=varargin{carg};carg=carg+1;
  model=varargin{carg};carg=carg+1;
  RunOpt.FieldDofs=[1 2 3]; % 3 DOF case all the time
  out=RunOpt;
elseif comstr(Cam,'cleangauss')
%% #cleanGauss
  match=varargin{carg};carg=carg+1;
  C1=varargin{carg};carg=carg+1;
  ind=reshape(1:size(match.Node,1),size(C1.X{2},1)/2,[]);ind(:,2:2:end)=[];
  ind=ind(:);in2=ind+size(C1.X{2},1)/2;
  st={'Node','jdet','wjdet','ebas'};
  for j1=1:length(st);
   r2=match.(st{j1});
   if size(r2,1)==length(ind)*2;r2=r2(ind,:);else;r2=r2(:,ind);end
   match.(st{j1})=r2;
  end
  st=setdiff(fieldnames(match),[st,{'Lambda','DOF'}]);
  for j1=1:length(st);
   r2=match.(st{j1});match.(st{j1})=r2(ind,:);
  end
  r2=match.Lambda{4};r2.N=r2.N(1:size(r2.N,1)/2,:);r2.w(size(r2.N,1)+1:end,:)=[];
  r2.mpid(in2,:)=[];
  match.Lambda{4}=r2;
  out=match;
  
%% #View
elseif comstr(Cam,'view');[CAM,Cam]=comstr(CAM,5);
    
  cf=varargin{carg};carg=carg+1;
  RO=varargin{carg};carg=carg+1;
  if ~isfield(RO,'MatId')
   il=cf.mdl.il;
   for j1=1:size(il,1); dbstack; keyboard; end
  end
  
  model=cf.mdl.GetData;
  RA=struct('type','Gauss','sel',sprintf('MatId %i',RO.MatId));
  sel=fe_caseg('StressCut-SelOut',RA,model);

  mo3=model;mo3.Elt=feutil(sprintf('selelt matid %i',RO.MatId),model);
  mo3.Elt=feutil('addelt','quad4',mo3.Elt(isfinite(mo3.Elt(:,1)),[1:4 9:10]));
  if isfield(RO,'il');mo3.il=RO.il;
  else;mo3.il=[RO.MatId(1) fe_mat('p_contact','SI',1) 0 -3];
  end
  [i1,mo3.Elt]=feutil('eltidfix;',mo3);
  mo3=feval(fe_fmesh('@GaussMesh'),'2d',mo3);
  mo3=feutil('getpatchnew',mo3);sel.fs=mo3.fs; sel.ifs=mo3.ifs;
  sel=sdsetprop(sel,'fsProp','FaceColor','interp','edgecolor','none');
  sel.StressObs.CritFcn='r1=r1(3,:,:).'';lab=''gap'';';
  sel.f1=[];sel.if1=[];
  fe_caseg('stresscut',sel,cf) % Overlay view and nominal mesh

  if isfield(RO,'sel1')
   cf.sel(1)={sprintf('matid~=%i',RO.MatId),'coloredgek -edgealpha.2'};
  end
% -------------------------------------------------------------------------
%% #Test basic test : p_zt('testElt')
elseif comstr(Cam,'test');[CAM,Cam]=comstr(CAM,5);

if comstr(Cam,'elt')    
%% #TestElt Do not integrate example into femesh -2
% If you define your new p_function no need for new geometry
model=struct('Node',[1 0 0 0  0 0 0;2 0 0 0  1 0 0;3 0 0 0  1 1 0;4 0 0 0  0 1 0;
 5 0 0 0  0 0 0; 6 0 0 0  1 0 0;7 0 0 0  1 1 0;8 0 0 0  0 1 0], ...
 'Elt',feutil('addelt','hexa8b',[1:8 103 114]), ...
 'pl',[103 fe_mat('m_elastic','SI',4) 1:6  -1 -2 -11 -12 -13  -100 (1:6)*10], ...
 'il',[114 fe_mat('p_zt','SI',1) 0 -3]);

[Case,model.DOF]=fe_mknl('init',model);
EC=Case.GroupInfo{end};
constit=Case.GroupInfo{4};
if ~isequal(constit(EC.ConstitTopology{1}),[1,2,4;2,3,5;4,5,6]);error('change');end
if ~isequal(constit(EC.ConstitTopology{3}),[1,2,4;2,3,5;4,5,6]*10);error('change');end

% Verify EC
EC.nodeE=model.Node(:,[5:7 1]);
z=integrules('buildndn',23,EC);
if norm(z.jdet-.25);error('Mismatch');end

sp_util('diag',12); % Matlab implement (but no orient)
[Case,model.DOF]=fe_mknl('init',model);k1=fe_mknl('assemble',model,Case,1);figure(1);spy(k1)
rb=feutilb('geomrb',model,[0 0 0],model.DOF);
if norm(k1*rb.def)>1e-10;error('Mismatch');end

sp_util('diag',0); % Mex implemented with orient
[Case,model.DOF]=fe_mknl('init',model);k2=fe_mknl('assemble',model,Case,1);
if norm(k2-k1,'inf')>1e-10;error('Mismatch');end


% Now check in different plane
mo2=model;mo2.Node(:,5:7)=model.Node(:,[6 7 5]);
[Case,mo2.DOF]=fe_mknl('init',mo2);k3=fe_mknl('assemble',mo2,Case,1);
i2=1:3;i3([2 3 1])=1:3;
if norm(full([k2(i2,i2)-k3(i3,i3)]));error('Reordering problem -1');end
rb=feutilb('geomrb',mo2,[0 0 0],mo2.DOF);
if norm(k3*rb.def)>1e-10;error('Mismatch -1');end
%save o:/balmes/xxx1.mat 

% GV: call fe_stress
try;
 R1=fe_stress('stress-gstate',model,rb);C1=R1.GroupInfo{1,5};
 r1=reshape(C1.Y(:,1,:,:),[],6);
end
% Lets call this again check assy is noy compromized (.CTable)
mo2=model;mo2.Node(:,5:7)=model.Node(:,[6 7 5]);
[Case,mo2.DOF]=fe_mknl('init',mo2);k3=fe_mknl('assemble',mo2,Case,1);
i2=1:3;i3([2 3 1])=1:3;
if norm(full([k2(i2,i2)-k3(i3,i3)]));error('Reordering problem -2');end
rb=feutilb('geomrb',mo2,[0 0 0],mo2.DOF);
if norm(k3*rb.def)>1e-10;error('Mismatch -2');end

  

else
 %% #Simple consistence check with stiffness in 3 directions -2
 mo0=femesh('teststruct hexa8;');
 if comstr(Cam,'tet'); mo0=feutil('hexa2tetra',mo0);
  if comstr(Cam,'tet10'); mo0=feutil('lin2quad',mo0); [CAM,Cam]=comstr(CAM,6);
  else; [CAM,Cam]=comstr(CAM,4);
  end
 end
 model=mo0;model.Elt=feutilb('SeparatebyMat',model.Elt);
 RA=struct('type','group','sel1',1,'sel2',2, ...
     'PostFcn',{{'p_zt','unjoin'}},'pzid',10);
 model=feutil('unjoin',model,RA);
 
 model=feutil('setpro',model,[10 fe_mat('p_zt','SI',1) 0 -3]);
 % Stiff 3rd direction enforces sliding. 1 & 2 give sliding stiffness
 pl=[10 fe_mat('m_elastic','SI',4) 2e8 0 3e8 0 0 1e14 1e-10];
 pl(15:20)=[2e8 0 3e8 0 0 1e14]; % viscous coupling
 model=feutil('setmat',model,pl);
 model=fe_case(model,'reset');model.DOF=[];
 if comstr(Cam,'model'); out=model; return; end
 
 [C2,model.DOF]=fe_mknl('init',model);
 k3=fe_mknl('assemble',model,C2,1);
 c3=fe_mknl('assemble',model,C2,3);
 if 1==1 % Check viscous surface damping and weight 
  mo2=model;mo2.Elt=feutil('selelt group3',model);mo2.Elt(3:end,:)=[];
  mo2.il(end,3)=5; % large rule
  mo2.K=[];[mo2,C2]=fe_case(mo2,'assemble -matdes 1 2 3 -SE');
  if normest(mo2.K{1}-mo2.K{3})>1e-5; error('Mismatch');end
  c=sum(fe_c(mo2.DOF,.01));S=1/3*1/9;
  if abs(c*mo2.K{2}*c'*1e10/S-1)>1e-4; error('Wrong weight');end

 elseif 1==2
  EC=C2.GroupInfo{end};
  dd=feutil('GetDD',[10 10],model)
  %  rb=feutilb('geomrb',model,[0 0 0],model.DOF);
  % k3*rb.def
 end
 
 def=fe_eig(fe_case(model,'reset'),[5 15 1e3]);
 if any(feutilb('dtkt',def.def,k3)<-1e-4)
  error('p_zt assy problem')
 end
 
 if sdtdef('isinteractive')
  %cf=feplot(model,def); fecom('ch7');fecom colordataevala
 end
 mo1=model;mo1.Elt=feutil('selelt group3',mo1);
 R1=fe_stress('stress-gstate',mo1,def);C1=R1.GroupInfo{1,5};
 
 % Check stress evaluation
 rb=feutilb('geomrb',model);
 rb.def(fe_c(rb.DOF,feutil('findnode matid 100',model),'ind'),:)=0;
 mo2=model;mo2.Elt=feutil('selelt proid10',model);
 R1=fe_stress('stress-gstate',mo2,rb);C1=R1.GroupInfo{1,5}; 
 r1=reshape(C1.Y(:,1,:,:),[],6);
 % slide 1 = y, slide 2=z, gap=x
 r0=[0 2e8 0;0 0 3e8;1e14 0 0];
 if norm(r1(1:3,1:3)-r0)>1e-1; 
   fprintf('\n\n');disp([r1(1:3,1:3)./r0]);
   error('Mismatch');
 end
 %r2=reshape(C1.Y(:,2,:,:),[],6);
 %r3=reshape(C1.Y(:,3,:,:),[],6);
 if 1==1 % Prepare for CLIMA
  model=feutil('setmat 10 rho=1',model);feutilb('_write',model);
  r2=fe_caseg('StressCut-SelOut',struct('type','Gauss','sel','MatId10'),model);
  sens=feutil('placeindof',rb.DOF,r2.StressObs);
  r3=reshape(sens.cta*rb.def,size(C1.Y));r3=reshape(r3(:,1,:,:),[],6);
  if norm(r3-r1,'inf')>.1;error('Mismatch'); end
  
  cf=feplot(model,rb);
  p_zt('view',cf,struct('MatId',10,'il',[10 fe_mat('p_contact','SI',1) 0 -3]));
  
  %r1=ctc_utils('getobs',mo2,rb)
  %[model,Case]=fe_case(mo3,nl_spring('AssembleCall'));

 end
 
end
elseif nargin>2&&ischar(varargin{3});
  obj=varargin{1};RO=varargin{2};[CAM,Cam]=comstr(varargin{3},1);carg=4;
 if comstr(Cam,'unjoin')
 %% #Unjoin : build p_zt faces during unjoin
 model=obj;
 i2=feutil('findnode groupall',model.Node,RO.el2); % Nodes of second group
 if isnumeric(RO.sel1); % common faces 
  elt=feutil('selelt group & selface & innode',model,RO.sel1,i2); 
 else
  elt=feutil(sprintf('selelt %s & selface & innode',RO.sel1),model,i2); 
 end
 
 RO.nind=sparse(RO.inode(:,1),1,RO.inode(:,2));
 [EGroup,nGroup]=getegroup(elt);
 if ~isfield(RO,'ztid');RO.ztid=10;end
 if length(RO.ztid)<2;RO.ztid=ones(1,2)*RO.ztid;end
 if isempty(elt);sdtw('_nb','Empty face of group %i for unjoin',RO.sel1);end
 
 mo2=model;mo2.Elt=elt;mo2=feutil('optimdegen',mo2);
%  r1=fe_quality('measAspectRatio -silent2',mo2);r1=cellfun(@(x,y)y(~isfinite(x)),r1.data,r1.EltId,'uni',0);
%  [mo2.Elt,elt]=feutil('removeelt eltid',mo2,vertcat(r1{:}));
 elt=mo2.Elt;
 model=addZt(model,elt,RO);

 out=model;
 elseif comstr(Cam,'surf2zt')
 %% #Surf2Zt : add zt element on surface
 model=obj;
 if ~isfield(RO,'ztid');RO.ztid=[10 10];end
 elt=feutil('selelt',model,RO.sel);
 n1=feutil('getnode groupall',model.Node,elt);i1=n1(:,1);
 if ~isfield(RO,'N0');RO.N0=max(model.Node(:,1))-min(n1(:,1));end
 n1(:,1)=RO.N0+n1(:,1);
 model.Node=[model.Node;n1];
 RO.nind=sparse(i1,1,n1(:,1));
 model=addZt(model,elt,RO);
 out=model;

 else;sdtw('''%s'' not known',CAM);
 end 
  
%% #End ----------------------------------------------------------------------
elseif comstr(Cam,'tablecall');out='';
elseif comstr(Cam,'cvs');
 out='$Revision: 1.25 $  $Date: 2022/10/26 12:49:16 $'; return;
else;sdtw('''%s'' not known',CAM);
end
end %fcn

function model=addZt(model,elt,RO);
 %% addZt Add layer from renumbering in RO.nind
 nind=RO.nind; 
 [EGroup,nGroup]=getegroup(elt);
 for jGroup=1:nGroup
   [ElemF,i1,ElemP]= getegroup(elt(EGroup(jGroup),:),jGroup);
   cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
   % use first order topo anyways
   if strncmpi(ElemP,'quad',4)% Add zt elements that have hexa8 topology
    model.Elt=feutil('addelt',model.Elt,'hexa8', ... 
     [ elt(cEGI,1:4) reshape(full(nind(elt(cEGI,1:4))),length(cEGI),[]) ...
           ones(length(cEGI),1)*RO.ztid(:)']);
   elseif strncmpi(ElemP,'tria',4)
    model.Elt=feutil('addelt',model.Elt,'penta6', ... 
     [ elt(cEGI,1:3) full(nind(elt(cEGI,1:3))) ...
           ones(length(cEGI),1)*RO.ztid(:)']);
   else; error('%s Not yet implemented',ElemP)
   end
 end
end
