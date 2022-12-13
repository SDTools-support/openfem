function [out,out1,out2,out3]=p_solid(varargin)

%P_SOLID volume element property
%
%       Syntax : pro= p_solid('default') 
%                pro= p_solid('database matid (command)') 
%                il = p_solid('dbval ProId (command)');
%                il = p_solid('dbval -unit MM ProId (command)');
%
%       Accepted commands for database and dbval calls are
%         d3 2 : 2x2x2 integration rule for linear volumes (hexa8 ... )
%         d3 3 : 3x3x3 integration rule for quadratic volumes (hexa20 ... )
%         d2 2 : 2x2x2 integration rule for linear volumes (q4p ... )
%         d2 3 : 3x3x3 integration rule for quadratic volumes (q8p ... )
%
%       Properties supported by p_solid are
%       Subtype 1 : 3D solid volume element: hexa8, hexa20, penta6, tetra4, ...
%            [ProId Type CordM Integ Stress Isop] 
%          CordM : material coordinate system (unused)
%          Integ : if non zero, selects the integration rule in 
%            integrules(ElemP,Integ) for example integrules(hexa8('parent'),3)
%            if Integ=0 or property undefined old of_mk_subs elements are used 
%          Stress : stress integration rule (unused)
%          ISOP : geometry interpolation flag (unused)
%
%       Subtype 2 : 2D volume element : q4p, q8p, t3p, t6p, ... 
%         [ProId Type Form N Integ]
%          with 
%           Type = fe_mat('p_solid','SI',2)
%           Form : formulation (0 plane strain, 1 plane stress, 2 axisymetric)
%           N    : Fourier harmonic for axisymetric elements that support it
%           Integ: see above
%
%       Subtype 3 : ND-1 volume, surface load or pressure
%         [ProId fe_mat('p_solid','SI',3) Integ Type Ndof1 ...]
%           Integ see above
%           Type: 1 volume force, 2 density of volume force
%                 3 pressure,     4 fluid/structure coupling
%                 5 2D force      6 2D pressure
%
%       See sdtweb      fem (handling materials section), pl, fe_mat, p_shell
%       See also help   fe_mat


%       Jean-Michel Leclere, Etienne Balmes
%       Copyright (c) 2001-2022 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*NOSEM>

if nargin<1; help p_solid;return; end
if ischar(varargin{1}); [CAM,Cam]=comstr(varargin{1},1); il=[]; carg=2;
else; il=varargin{1};[CAM,Cam]=comstr(varargin{2},1); carg=3;
end
persistent eM
if isempty(eM)
  eM.topType=fe_mat('@topType');
  eM.tomType=fe_mat('@tomType');
  eM.TensorT=m_elastic('@MechaTensorT');
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
out1={}; % See if unit specified
 while 1==1
  i1=strfind(Cam,'-unit'); out1={};
  if ~isempty(i1)
   [Unit,i2,i3,i4]=sscanf(CAM(i1+5:end),'%s',1);
   i4=i1+(0:4+i4-1);CAM(i4)=''; [CAM,Cam]=comstr(CAM,1);
  else; Unit='';
  end
  i2=strfind(Cam,'-punit');
  if ~isempty(i2)
   [RO.PUnit,i3,i4,i5]=sscanf(CAM(i2+6:end),'%s',1);
   i5=i2+[0:5+i5];CAM(i5)=''; [CAM,Cam]=comstr(CAM,1);
  else; RO.PUnit='';
  end
  if ischar(CAM); [CAM,Cam,RO.isop]=comstr('isop',[-25 1],CAM,Cam);
  else; RO.isop=[];
  end

  if ischar(CAM); [i1,CAM,Cam]=comstr(CAM,'dbval','%i');else; i1=[];end
  if isempty(CAM)&&carg<=nargin; st=varargin{carg};carg=carg+1;
  else; st=CAM;end
  if isempty(st);
  elseif ischar(st); mat=p_solid('database',st);
  elseif isnumeric(st); [typ,st1,i4]=eM.topType(st(2));
   mat=struct('il',st,'name',sprintf('%i',st(1)),'type',typ,'unit',st1);
  end
  if ~isempty(RO.PUnit)
   r1=fe_mat(sprintf('convert %s %s',mat.unit,RO.PUnit),mat.il(1:2));
   mat.il(2)=r1(2); mat.unit=RO.PUnit;
  end
  if ~isempty(Unit)
   mat.il=fe_mat(sprintf('convert %s %s',mat.unit,Unit),mat.il);mat.unit=Unit;
  end
  if ~isempty(RO.isop);mat.il(6)=RO.isop;end
  if length(i1)==1; mat.il(1)=i1;end
  r1=mat.il;
  if ~isempty(il); i2=find(il(:,1)==r1(1)); else; i2=[];end
  if isempty(i2); i2=size(il,1)+1;end
  il(i2,1:length(r1))=r1; %#ok<AGROW>
  out1(end+1,1:3)={'pro',sprintf('%i_%s',mat.il(1),mat.name),mat};%#ok<AGROW> 
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
  out.il=[MatId fe_mat('p_solid','SI',1) 0 -3     0     1    0]; 
  out.name='Default for topology';
  out.type='p_solid';
  out.unit='SI';
  
  % MatId Typ                     CordM Integ Stres Isop Field
  mat.il=[MatId fe_mat('p_solid','SI',1) 0 2     0     1    0]; 
  mat.name='Full 2x2x2';
  mat.type='p_solid';
  mat.unit='SI';out(end+1)=mat;

  mat.il=[MatId fe_mat('p_solid','SI',1) 0 2 0 1 0]; 
  mat.name='Reduced shear';
  mat.type='p_solid';
  mat.unit='SI';
  out(end+1)=mat;

  mat.il=[MatId fe_mat('p_solid','SI',2) 1 0 0]; 
  mat.name='Plane stress';
  mat.type='p_solid';
  mat.unit='SI';
  out(end+1)=mat;
  % Smarter parameter input
  if sp_util('issdt')
   [RO,Cam,st]=cingui('paramedit -DoClean',[ ...
   'COORDM(0#%g#"Material CID") ISOP(0#%g#"Integration scheme") '...
   ],st);
  else; % OpenFEM only values ignored 
      RO=struct('COORDM',0,'ISOP',0);
  end
  
  i1=find(strncmpi(st,{out.name},length(st)));
  if ~isempty(i1);out=out(i1);
  elseif comstr(comstr(st,-27),'d3')
   r1=comstr(st(3:end),[-1 3]);
   if length(r1)==1; r1=[RO.COORDM r1 0 RO.ISOP];end
   out=struct('il',[MatId fe_mat('p_solid','SI',1) r1(:)'], ...
         'name',horzcat('Solid',st),'type','p_solid','unit','SI');
  elseif comstr(comstr(st,-27),'d2')
   if strncmpi(st,'d23',3); r1=comstr(st(4:end),[-1 2]);
    if length(r1)==1; r1=[231 0 r1 0];end % Use surface gradient
   else; r1=comstr(st(3:end),[-1 3]);
    if length(r1)==1; r1=[0 0 r1 0];end
   end
   out=struct('il',[MatId fe_mat('p_solid','SI',2) r1(:)'], ...
         'name',horzcat('Solid',st),'type','p_solid','unit','SI');
  elseif comstr(comstr(st,-27),'fsc')
   r1=comstr(st(4:end),[-1 -3]);r1(2)=4;
   out=struct('il',[MatId fe_mat('p_solid','SI',3) r1(:)'], ...
         'name',st,'type','p_solid','unit','SI');
  elseif strncmpi(st,'wallimp',7)
   r1=[-3 8 comstr(st(8:end),-1)];
   out=struct('il',[MatId fe_mat('p_solid','SI',3) r1(:)'], ...
         'name',st,'type','p_solid','unit','SI');
  elseif ~isempty(st)
   r1=str2num(st);  %#ok<ST2NM>
   if length(r1)<1; error('Not a consistent database call');end
   out=[];
   out.il=[MatId fe_mat('p_solid','SI',1) r1(:)']; 
   out.name=['Solid' st];
   out.type='p_solid';
   out.unit='SI';
  end

  out1='Solid';

%% #Conv convention labels ---------------------------------------------------
% see also sdtweb fe_stress('TensorTopology')
% WARNING these conventions are SOFTWARE DEPENDENT, sdtweb feform#feelas3d
elseif comstr(Cam,'conv')

out=struct( ...
    'Stress3D',{{'Sxx','Syy','Szz','Syz','Szx','Sxy'}}, ...
    'Stress3DL',{{'\sigma_{xx}','\sigma_{yy}','\sigma_{zz}', ...
    '\sigma_{yz}','\sigma_{zx}','\sigma_{xy}'}}, ...
    'Strain3D',{{'Exx','Eyy','Ezz','Gyz','Gzx','Gxy'}}, ...
    'Strain3DL',{{'\epsilon_x','\epsilon_y','\epsilon_z', ...
     '\gamma_{yz}','\gamma_{zx}','\gamma_{xy}'}}, ...
    'Strain2D',{{'Exx','Eyy','Gxy'}}, ...
    'Strain2DL',{{'\epsilon_x','\epsilon_y','\gamma_{xy}'}});
 
[CAM,Cam]=comstr(CAM,5); if ~isempty(CAM); out=out.(CAM);end

%% #PropertyUnitType  --------------------------------------------------------
elseif comstr(Cam,'propertyunittype')

 if nargin==1;out=1:3; return; end % return subtypes ID
 i1=varargin{carg};
 out1={};
 switch i1 % PropertySubType
 case 1 % [ProId type Coordm In Stress Isop Fctn  ]
   st={ ...
   'ProId'   0  'sdtweb(''p_solid'')';
   'Type'    0  '';
   'COORDM'  0  'Coordinates system id';
   'IN'      0  'Integration rule';
   'STRESS'  0  'Location for stress output';
   'ISOP'    0  'Integration scheme';
   'FCTN'    0  ''};
 case 2 % 2D   [ProId Type Form N]
   st={ ...
   'ProId'   0  'sdtweb(''p_solid'')';
   'Type'    0  '';
   'Form'    0  'Formulation';
   'N'       0  'Fourier harmonic'
   'In'      0  'Integration rule'};
 case 3 % ND-1 coupling element
   st={ ...
   'ProId'   0  'sdtweb(''p_solid'')';
   'Type'    0  '';
   'Integ'   0  'Integration rule';
   'Form'    0  'Formulation'};
   out1={'Ndof%i'        0  'NodeId DOF coupling';};
 otherwise; st={'ProId' 0 'sdtweb(''p_solid'')'; 'Type', 0, ''};
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
 if carg<=nargin; constit=varargin{carg};carg=carg+1;else; constit=[];end
 if carg<=nargin; model=varargin{carg};carg=carg+1;else; model=[];end
 if carg<=nargin; Case=varargin{carg};carg=carg+1;else; Case=[];end
 if ~isempty(constit)&&constit(1)==-1; 
  st=eM.topType(constit(2,1));[out,out1,out2]=feval(st,'const',varargin{2:end});
  return;
 end
 if isstruct(varargin{end});RO=varargin{end};else;RO=struct;end;RO.mtt='';
 if ~ischar(EC) % Allow for integrule switching here
 else % Classical init
  if size(integ,1)>6&&integ(6,1)==0&&integ(5,1)==0&&integ(7,1)~=0
    EC=integrules(EC,double(integ(7,1)));RO.rule=integ(7,1);
    if any(integ(7,:)~=integ(7,1));
     warning('OpenFEM:MULTINTEG', ...
        'Cannot deal with multiple integration strategies in the same group');
    end
  elseif size(integ,1)>4;
    EC=integrules(EC,double(integ(5,1)));RO.rule=integ(5,1);
    if any(integ(5,:)~=integ(5,1));
     warning('OpenFEM:MULTINTEG', ...
        'Cannot deal with multiple integration strategies in the same group');
    end
  else; EC=integrules(EC);
  end
 end
 % already initialized, this is used for formulations contained in PRO
 try;
 if isfield(Case,'GroupInfo')&&isstruct(Case.GroupInfo{Case.jGroup,8}) 
    r1=Case.GroupInfo{Case.jGroup,8};
    st=setdiff(fieldnames(EC),fieldnames(r1));
    st=setdiff(fieldnames(EC),{'DofLabels'});
    for j1=1:length(st); r1.(st{j1})=EC.(st{j1});end
    EC=r1; 
    if ~isfield(EC,'nodeE'); % Missing InfoAtNode fields
        EC.nodeE=zeros(size(EC.N,2),4); 
    end
    if ~isfield(EC,'defe'); 
        EC.defe=[]; 
    end
    if isfield(EC,'ke')
    elseif isfield(EC,'FieldDofs');    EC.ke=zeros(length(EC.FieldDofs)*size(EC.N,2));
    else;EC.ke=zeros(length(EC.DofLabels)*size(EC.N,2));
    end
 end
 end
 
 out2=[]; if ~isempty(Case);out2=Case.GroupInfo{Case.jGroup,7};end
 Ndof=0;
 if size(integ,1)<3; Ndof=0;
 elseif integ(3); Ndof=integ(3,1);
 elseif ~any(integ(3:4,1))% back propagate number of DOFs if needed
     if isfield(EC,'DofLabels');Ndof=length(EC.DofLabels)*EC.Nnode; 
     elseif ~any(EC.xi(:,3)); Ndof=2*EC.Nnode; 
     elseif norm(EC.xi(:,2:3),'inf')==0; Ndof=EC.Nnode; 
     else;  Ndof=3*EC.Nnode;
     end
 end
 EC.Nw=size(EC.N,1);out.Nnode=size(EC.N,2);
 if size(integ,1)<6||any(integ(5:6,1)==0)||integ(6)>size(EC.w,1);
       rule=[1 EC.Nw]; % By default start at 1 use all points
 else; rule=integ(5:6,1);rule=rule(:)';
 end
 if ~any(rule); 
   rule=[1 EC.Nw];
   sdtw('integ is assumed to give standard integration rule in integ(5:6)');
 elseif rule(1)==0||rule(end)>size(EC.w,1);error('Inconsistent rule'); 
 end
 % integ(1)==-1 is done at the very beginning of the command
 RO.CaseUrn=sprintf('PerNode%g',Ndof/EC.Nnode);
 if ~isempty(constit)&&constit(1)==-3
   %% #interface_element - - - - - - - - - - - - - - - - - - - - - - - - -2
   if ~isfield(EC,'StrainString')      
   switch constit(4,1)
   case 4 % Fluid/structure coupling
     EC.DofLabels={'u','v','w','p'};
     EC.material='fs_matrix'; EC.bas=zeros(9,size(EC.N,1));
     Ndof=4*EC.Nnode;
     EC.VectMap=int32([1:4:Ndof 2:4:Ndof 3:4:Ndof 4:4:Ndof])';
     EC.ElMap=zeros(Ndof);%reshape(1:Ndof^2,Ndof,Ndof);
     EC.ElMap(EC.VectMap,EC.VectMap)=reshape(1:Ndof^2,Ndof,Ndof);
     out1=23;
   case {1,2,3} % volume or pressure load on surface
     EC.DofLabels={'u','v','w'};
     EC.material=''; EC.bas=zeros(9,size(EC.N,1));
     Ndof=3*EC.Nnode;
     EC.VectMap=int32([1:3:Ndof 2:3:Ndof 3:3:Ndof])';
     out1=23;
   case {5,6} % volume or pressure load on 2D line
     EC.DofLabels={'u','v'};
     EC.material=''; EC.bas=zeros(4,size(EC.N,1));
     Ndof=2*EC.Nnode;
     EC.VectMap=int32([1:2:Ndof 2:2:Ndof])';
     out1=12;
   case {8} % acoustic impedance on a wall
     EC.DofLabels={'p'};
     EC.bas=zeros(9,size(EC.N,1));
     EC.StrainDefinition{3}= [1 1 1 rule];
     EC.StrainLabels{3}={'p'}; EC.ConstitTopology{3}=int32(5);
     EC.VectMap=int32(reshape(1:EC.Nnode,1,EC.Nnode)'); 
     EC=integrules('matrixrule',EC);
     out1=23;
      
   otherwise
      sdtw('Surface element %i not implemented',constit(4,1));
      out1=23;
   end
   else % formulation within pro
      EC=integrules('strainstring',EC,rule);
      if any(EC.xi(:,3));out1=3;elseif any(EC.xi(:,2));out1=2;else;out1=13;end
   end
%% #2D_ELASTIC SOLID - - - - - - - - - - - - - - - - - - - - - - - - - - - -2
 elseif strcmpi(RO.CaseUrn,'PerNode2')
  RO.rule=rule; [EC,RO]=EC_Elas2D(EC,RO);      
 
 %% #3D_SOLID elastic rule - - - - - - - - - - - - - - - - - - - - - - - - - - 
 elseif strcmpi(RO.CaseUrn,'PerNode3')

 EC.DofLabels={'u','v','w'};
 if ~any(EC.xi(:,3)) 
  %% Surface strain
  if 1==1 % non flat requires 23 rule and orientation
   EC.StrainDefinition{1}= ...
    [1 2 1 rule; 2 2 2 rule;3 2 3 rule; % 
     4 3 1 rule; 5 3 2 rule;6 3 3 rule]; % shear
   EC.StrainLabels{1}={'u,x','v,x','w,x','u,y','v,y','w,y'};
   EC.ConstitTopology{1}=int32(3*eye(6));
   EC.StrainLabels{2}={'u','v','w'};
   EC.StrainDefinition{2}= [1 1 1 rule;2 1 2 rule;3 1 3 rule];
   EC.ConstitTopology{2}=int32(eye(3));
  else
   EC.StrainDefinition{1}= ...
    [1 2 1 rule; 2 3 2 rule; % e_xx = N,x u, e_yy, e_zz
     3 3 1 rule;3 2 2 rule]; % shear
   EC.StrainLabels{1}=p_solid('ConvStrain2D');
   EC.ConstitTopology{1}=int32(diag([3 10 24]));'zzz need fixing in buildconstit'
  end
  %dd=double(EC.ConstitTopology{1});dd(dd~=0)=constit(dd(dd~=0))
  EC=integrules('matrixrule',EC);
  EC.bas=zeros(9,size(EC.w,1));
  out1=23; % Tell that BuildNDN rule is 3D
  EC.VectMap=int32(reshape(1:3*EC.Nnode,3,EC.Nnode)'); 
 
 elseif size(integ,1)<9||remi(integ(9,1),[],3)==0 % Isop is not 1xx
  RO.rule=rule; [EC,RO]=EC_Elas3D(EC,RO,integ,constit);

 elseif size(integ,1)>8&&remi(integ(9,1),[],3)==1 % ID(9)=1xx isop==100 choice 
 %% #Strain=ui,j and user defined non-linearity -2
  RO.rule=rule; RO.RunOpt=varargin{end};
  [EC,RO]=EC_Elas3Dld(EC,RO,integ,constit); 

 else; error('Not a valid case');
 end
 %% #ConstAcoustic2D fluid acoustic : as many DOFs as nodes - - - - - - - - - - - - - - - -
 elseif integ(3)==double(integ(4))&&~any(EC.w(:,3)) ...
   &&(size(integ,1)<=4 ||integ(7)==2||(size(integ,1)>8&&integ(9,1)==2))
 EC.DofLabels={'p'}; out1=2;% Tell that BuildNDN rule is 2D

 % Strain energy
 % define the deformation vector: row, NDN, DDL, NwStart NwRule
 EC.StrainDefinition{1}=[1 2 1 rule;2 3 1 rule];
 EC.StrainLabels{1}={'p,x','p,y'};
 EC.ConstitTopology{1}=int32(diag([3 3]));
 EC.StrainDefinition{5}=EC.StrainDefinition{1};
 EC.StrainLabels{5}=EC.StrainLabels{1};
 EC.ConstitTopology{5}=EC.ConstitTopology{1};
 % Kinetic energy
 EC.StrainDefinition{2}= [1 1 1 rule];
 EC.StrainLabels{2}={'p'}; EC.ConstitTopology{2}=int32(1);
 EC.VectMap=int32(reshape(1:EC.Nnode,1,EC.Nnode)'); 

 EC=integrules('matrixrule',EC);
 out1=2; % Tell that BuildNDN rule is 2D
 %% #ConstSurfHeat  - - - - - - - - - - - - - - - -
 elseif integ(3)==double(integ(4))&&~any(EC.w(:,3)) ...
   &&(size(integ,1)<=4 ||integ(7)==231||(size(integ,1)>8&&integ(9,1)==2))
 EC.DofLabels={'p'}; out1=23;% Tell that BuildNDN rule is 2D

 % Strain energy
 % define the deformation vector: row, NDN, DDL, NwStart NwRule
 EC.StrainDefinition{1}=[1 2 1 rule;2 3 1 rule];
 EC.StrainLabels{1}={'p,x','p,y'};
 EC.ConstitTopology{1}=int32(diag([3 3]));
 EC.StrainDefinition{5}=EC.StrainDefinition{1};
 EC.StrainLabels{5}=EC.StrainLabels{1};
 EC.ConstitTopology{5}=EC.ConstitTopology{1};
 % Kinetic energy
 EC.StrainDefinition{2}= [1 1 1 rule];
 EC.StrainLabels{2}={'p'}; EC.ConstitTopology{2}=int32(1);
 EC.VectMap=int32(reshape(1:EC.Nnode,1,EC.Nnode)'); 

 EC=integrules('matrixrule',EC);
 out1=231; % Tell that BuildNDN rule is 2D

 %% #Const_3D_fluid : as many DOFs as nodes - - - - - - - - - - - - - - - - -
 % sdtweb p_solid BuildAcoustic3D
 elseif integ(3)==double(integ(4))&& ...
    (size(integ,1)<=4||integ(7)==3||(size(integ,1)>8&&integ(9,1)==3))

 EC.DofLabels={'p'};
 out1=size(EC.NDN,2)/EC.Nw-1;  % Tell that BuildNDN rule is 3D
 if out1==1;out1=13; % Thermal/fluid on a line
 % Strain energy
 % define the deformation vector: row, NDN, DDL, NwStart NwRule
 EC.StrainDefinition{1}=[1 2 1 rule];
 EC.StrainLabels{1}={'p,x'};
 EC.ConstitTopology{1}=int32(diag([3 3 3])); % 1/rho
 else
 % Strain energy
 % define the deformation vector: row, NDN, DDL, NwStart NwRule
 EC.StrainDefinition{1}=[1 2 1 rule;2 3 1 rule;3 4 1 rule];
 EC.StrainLabels{1}={'p,x','p,y','p,z'};
 EC.ConstitTopology{1}=int32(diag([3 3 3]));
 EC.StrainDefinition{5}=EC.StrainDefinition{1};
 EC.StrainLabels{5}=EC.StrainLabels{1};
 EC.ConstitTopology{5}=EC.ConstitTopology{1};

 end
 % Kinetic energy
 EC.StrainDefinition{2}= [1 1 1 rule];
 EC.StrainLabels{2}={'p'}; EC.ConstitTopology{2}=int32(1);% 1/rho/C2
 EC.VectMap=int32(reshape(1:EC.Nnode,1,EC.Nnode)'); 

 EC=integrules('matrixrule',EC);

 elseif strcmpi(RO.CaseUrn,'PerNode3.4')&&any(RO.rule==[30003 20002]); 
    [EC,RO]=feval(m_hyper('@elasUP'),'EC',EC,RO,integ,constit);
 else
  fprintf('Not a known p_solid(''const'') case'); error('Not a valid case');
 end % 2D or 3D - - - - - - - - - - - - - - - - - - - - - - - - - - -
  % Adjust Number of Dof and Node in integ if needed
 if size(integ,1)>3&&~any(integ(3:4,1))% back propagate number of DOFs if needed
     if size(integ,1)<4; error('Improper init');end
     sp_util('setinput',integ, ...
         [Ndof;EC.Nnode],int32(3:4),int32(1));
     EC.ElMap=elem0('elmapmat_og',[EC.Nnode Ndof/EC.Nnode]);
 end

 if size(EC.N,2)~=EC.Nnode;
  try;
      opt1=p_piezo('ReducedShear',EC,integ,constit);EC=opt1;
  end
 end
 % InfoAtNode=Case.GroupInfo{Case.jGroup,7}
 if ~isempty(constit)&&isfield(EC,'MatrixIntegrationRule')
  if ~isfield(Case,'GroupInfo')||isempty(Case.GroupInfo{Case.jGroup,7})
   % With no material basis keep terms needed for the current law.
   i1=find(any(constit,2))-1; 
   for j0=1:length(EC.MatrixIntegrationRule);
     r1=EC.MatrixIntegrationRule{j0};
     if ~isempty(r1)
     EC.MatrixIntegrationRule{j0}=r1(ismember(double(r1(:,5)),i1),:);
     end
   end
  else % Strategies for constitutive law interpolation
   EC=ctableGen(EC,Case.GroupInfo{Case.jGroup,7},ConstitLab(RO.mtt));
  end
 else;EC.MatrixIntegrationRule={};
 end
 if isfield(RO,'NdnDim');out1=RO.NdnDim;end
 if size(EC.NDN,2)==3*length(EC.jdet)&&out1==3;
  error('You have declared a 3D element while your EltConst.NDN is 2D');
 elseif isfield(RO,'NdnDim');out1=RO.NdnDim;
 elseif out1==2&&isfield(Case,'Node')&&~isempty(Case.Node)&& ... % force 23 rule
         any(Case.Node(:,7)); out1=23;
 end
 if nargout==0
  integrules('texstrain',EC);
  try; EC=integrules('stressrule',EC);integrules('texstress',EC);end
 else; out=EC;
 end


%% -------------------------------------------------------------------------
elseif comstr(Cam,'buildconstit');
% #BuildConstit  initializes constit,integ
% Implementation of elastic constitutive law building 3D
%[constit,integ,Inits,Data]=p_solid('buildconstit Nfield Nnode',ID,pl,il, ...
%    model,Case); This is called by the element 'integinfo' command

  RunOpt=struct('warn',{{}},'Dim',0,'ProStack',[],'DoElmap',1);
  ID=varargin{carg};carg=carg+1;out1=int32(ID);out3=struct;
  pl=varargin{carg};carg=carg+1; 
  if isempty(pl);mat=[];else;mat=pl(pl(:,1)==ID(1),:);end
  il=varargin{carg};carg=carg+1; 
  if isempty(il);pro=[];else;pro=il(il(:,1)==ID(2),:);end
  % Default dimension building
  if isempty(pro)||length(ID)==2
    [CAM,Cam]=comstr(CAM,13);
    if ~isempty(Cam);RunOpt.Dim=comstr(Cam,-1); end
  end
  if length(ID)>=4&&ID(4); RunOpt.Dim=[ID(3)/ID(4) ID(4)]; end
  if ~isempty(pro)
  elseif RunOpt.Dim(1)==3; % Default for legacy elements
      pro=[ID(2) fe_mat('p_solid','SI',1) 0 0 0 0 ]; 
  elseif RunOpt.Dim(1)==2;
      pro=[ID(2) fe_mat('p_solid','SI',2) 0 0 0 0 ]; 
  end
  
  
  if size(mat,1)>1; 
   error('Duplicate material with ID %i',ID(1));%mat=mat(1,:);
  end
  if all(size(pro)>[0 1])
    [RunOpt.ProF,RunOpt.ProU,RunOpt.pTyp]=eM.topType(pro(1,2));
  else; RunOpt.ProF='m_null';
  end
  model=[]; Case=[];out2=[];
  % allow propagation of more complex material information as stack entries
  if carg<=nargin; 
   model=varargin{carg};carg=carg+1;
   r1=stack_get(model,'mat'); 
   for j1=1:size(r1,1)
    try
     if isfield(r1{j1,3},'pl')&&r1{j1,3}.pl(1)==ID(1); RunOpt.MatStack=r1{j1,3};end
    end
   end
   r1=stack_get(model,'pro'); 
   for j1=1:size(r1,1)
    try
        if isfield(r1{j1,3},'il')&&r1{j1,3}.il(1)==ID(2); RunOpt.ProStack=r1{j1,3};end
    end
   end
   if carg<=nargin; Case=varargin{carg};carg=carg+1;end
  end

  if isempty(mat); 
   RunOpt.warn{end+1}=sprintf('MatId %i is not matched',ID(1));
   st='m_null';unit=1;typ=1;
  else;  [st,unit,typ]=eM.tomType(mat(2));
   %% #ElasIsoDefaults for elastic istropic material -2
   if strcmp(st,'m_elastic')&&typ==1 % Standard isotropic
    if length(mat)<6||mat(6)==0; mat(6)=mat(3)/2/(1+mat(4));end % G
    if mat(3)==0; mat(3)=mat(6)*2*(1+mat(4));end % E=0 use G info
    if size(mat,2)<7;mat(7)=0;end
   end
  end
  stt=sprintf('%s.%i',st,typ); 
  RunOpt.ProT=sprintf('%s.%i',RunOpt.ProF,RunOpt.pTyp);

  RunOpt.Failed=0;    out2=[]; 
  if isfield(RunOpt,'MatStack')&&isfield(RunOpt.MatStack,'type')
      st2=feval(RunOpt.MatStack.type,'propertyunittype -cell',typ);
      if isfield(RunOpt.MatStack,'constit') % mat.constit and mat.integ predefined
        out2=struct;
        out=RunOpt.MatStack.constit; 
        out1=RunOpt.MatStack.integ; % Predefined integ
        out2.elmap= ...
            elem0('elmapmat_og',[double(out1(4)) double(out1(3))/double(out1(4))]);
        if isfield(RunOpt.MatStack,'eval0');eval(RunOpt.MatStack.eval0);end
        % Used to forward eval assembly process
        if isfield(RunOpt.MatStack,'eval1');out2.eval=RunOpt.MatStack.eval1;end
        return
      elseif ~isempty(intersect(st2(3:end,1),fieldnames(RunOpt.MatStack))) 
       %% #interpfield if some fieldnames of mat correspond to pl entries -2
       constit=RunOpt.MatStack.pl(:); 
       if length(constit)==1; constit=pl(pl(:,1)==constit,:)';end
       r2=stack_get(Case,'DofSet','ThermalState','getdata');
       if isempty(r2);r2=stack_get(Case,'DofSet','InfoAtNode','getdata');end
       if ~isempty(r2); % initialize material interp
        out2=struct('RunOpt_NodePos',1, ...
         'RunOpt_GstateFcn', ...
         'm_elastic(''AtNodeGstate'');Case.GroupInfo{jGroup,7}=InfoAtNode;');
       else;out2=[];RunOpt=rmfield(RunOpt,'MatStack');
       end      
      end % One of the fields is interpolated
  end
  %% #BuildElastic3D solid  - - - - - - - - - - - - - - - - - - - - - - - - - -
  if strcmp(RunOpt.ProT,'p_solid.1')
   out3.ConstitLab=ConstitLab('m_elastic.1');
   if ~pro(3)
   elseif ~isfield(Case,'bas'); error('Missing Case.bas for material');
   else; out3.bas=Case.bas(Case.bas(:,1)==pro(3),:);
   end
   if strcmp(stt,'m_elastic.1') % Standard isotropic
   if isfield(RunOpt,'MatStack')&&isfield(RunOpt.MatStack,'E')&& ...
        isfield(RunOpt.MatStack,'old')
      % obsolete generation of constit columns for each element
      [constit,out2]=m_elastic('formulaENG2DD',mat([3 4 6]),RunOpt,model,Case);
   elseif RunOpt.Dim(1)==3&&(length(pro)<4||pro(4)==0) % legacy elastic solid
     [constit,out1]=fe_mat('plil of_mk 3d',ID,mat,il);
     warning('OpenFEM:legacy',sprintf( ...
      'OBSOLETE fixed rule 3D p_solid Mat %i /Pro %i\n%s',ID(1:2), ...
      'you should define a rule with p_solid(''dbval d3 -3'')')); %#ok<SPWRN>
     out1(3)=prod(RunOpt.Dim(1:2));out1(4)=RunOpt.Dim(2); 
     % elmap is now done in the 'constants' call
     out=constit(:); out2=[]; 
     return;  % return Constit, Integ, Elmap
   end
   if ~isempty(out2); constit(20:end)=[]; 
       % Test on constit length in sdtweb m_elastic IsotropicInterp
   else; % standard without thermal state usage
     % 22/08/06 formula corrected for consistency when G ne E/2(1+n)
    if ismember(pro(4),[20002 30003]);RunOpt.Cb=m_hyper('@elasUP');else;RunOpt.Cb='';end
    if ~isempty(RunOpt.Cb)
        [dd,out2]=feval(RunOpt.Cb,'elmapdd',RunOpt,model,Case);
    else
     dd=m_elastic('formulaENG2DD',mat([3 4 6]),RunOpt,model,Case);
     if pro(3);dd=m_elastic('formulaPlAniso -1',dd,out3.bas);end;out3.dd=dd;
    end
     
     %E=mat(3);n=mat(4);G=E/2/(1+n);
     %inv([1/E -n/E -n/E 0 0 0;-n/E 1/E -n/E 0 0 0;-n/E -n/E 1/E 0 0 0;
     %     0 0 0 1/G 0 0; 0 0 0 0 1/G 0;0 0 0 0 0 1/G])-dd
     if size(dd,1)>30 % Interpolation done in Eng2dd
      constit=dd;
     else
      if any(diag(out3.dd)<0);
       sdtw('Negative diagonal term in constitutive law MatId %i',mat(1));
      end
      constit=[mat(5);mat(7);out3.dd(:)]; % 1:38 fully anisotropic constitutive law
     end
     if size(mat,2)>=9&&any(mat(8:9))  % Thermal expansion and ref temp
      at=eye(3)*mat(8); constit=[constit;at(:);mat(9)];  % 39:48
      out3.at=at; out3.T0=mat(9);
     end
     if sp_util('diag')>1; elem0('DiagConst',sprintf('D%i',ID(1)),out3.dd);end
   end
   if length(ID)<4||RunOpt.Dim(1)==6; ID(3:4)=RunOpt.Dim(2)*[3;1];end% not 6 for surf strain

   %% #BuildElas #ElasAniso3D anisotropic mat. -3
     % m_elastic('propertyunittypecell',3)
   elseif strcmp(stt,'m_elastic.3') 

    if length(mat)<25; mat(25)=0;end
    dd=triu(ones(6));dd(dd~=0)=3:23;dd=dd+triu(dd,1)';
    dd=mat(dd);
    if pro(3);
     dd=reshape(eM.TensorT(reshape(out3.bas(7:15),3,3),dd),6,6);
     % dd=m_elastic('formulaPlAniso -1',dd,out3.bas);
    end
    out3.dd=dd;
    constit=[mat(24);mat(25);dd(:)];% rho,eta,dd
    %m_elastic('propertyunittypecell',3);ans([24 32])
    if sp_util('diag')>1; fprintf('\n D(MatId %i)\n',mat(1)); disp(dd); end
    if size(mat,2)>=32; % at(26:31) here thermal expansion tested by Boris
     at=mat([26 27 29;27 28 30;29 30 31]);constit=[constit;at(:);mat(32)];
     if sp_util('diag')>1; elem0('DiagConst',sprintf('At%i',ID(1)),at);end
     out3.at=at;
    end
    ID(3:4)=RunOpt.Dim(2)*[3;1];
   % #BuildElas #ElasOrtho3 3D orthotropic materials -3
   elseif strcmp(stt,'m_elastic.6')  % 
    if length(mat)<17; mat(17)=0; end
    RunOpt.MatId=mat(1);dd=m_elastic('formulaortho',mat(3:11),RunOpt,model,Case);
    if pro(3);dd=m_elastic('formulaPlAniso -1',dd,out3.bas);end;out3.dd=dd;
    if sp_util('diag')>1; elem0('DiagConst',sprintf('D%i',ID(1)),dd);end
    constit=[mat(12);mat(17);dd(:)]; % 1:38 fully anisotropic constitutive law
    ID(3:4)=RunOpt.Dim(2)*[3;1];
   %% #BuildAcoustic3D fluid Acoustic - - - - - - - - -
   elseif strcmp(stt,'m_elastic.2')
       
    if length(mat)<5; mat(5)=0.; end
    out3.ConstitLab={'1/rho/C2' 'eta' '1/rho'};
    constit = [1/mat(1,3)/mat(1,4)^2 mat(1,5) 1/mat(1,3)]';
    ID(3:4)=RunOpt.Dim(2);pro(1,6)=3;% Was ID(7)
    if sp_util('diag')>1; 
         elem0('DiagConst',sprintf('D%i',ID(1)),constit,'1/rho/C^2, eta, 1/rho');
    end
   % #BuildElasShell for 3D : -3
   elseif strcmp(stt,'m_elastic.5') % 
     error('MatId%i m_elastic subtype 5 is not supported for volumes',ID(1));
   %% #Build_m_elastic_Failed -3
   else; RunOpt.Failed=1;
   end
   r2=stack_get(Case,'DofSet','ThermalState','getdata');
   if ~isempty(r2);
       out2=struct('RunOpt_NodePos',1, ...
         'RunOpt_GstateFcn', ...
         'm_elastic(''AtNodeGstate'');Case.GroupInfo{jGroup,7}=InfoAtNode;');
   end

   % Integ Stress Isop saved in integ if non zero (5:6) is used for rule
   if size(pro,1)&&size(pro,2)>3;
     r2=pro(1,4:end);  ID(6+(1:size(r2,2)))=r2; 
   end
  % --------------------------------------------------------------------
  %% #Build2D solid ---------------------------------------------------------2
  elseif isequal(RunOpt.ProF,'p_solid')&&RunOpt.pTyp==2
    if strcmp(st,'m_elastic')&& ...
         RunOpt.Dim(1)==2&&(length(pro)<5||pro(5)==0) % legacy elastic solid
     %error('OpenFEM legacy elements have been discontinued');
     warning('OpenFEM:legacy',sprintf( ...
      'OBSOLETE fixed rule 2D p_solid Mat %i /Pro %i\n%s',ID(1:2), ...
      'you should define a rule with p_solid(''dbval d2 -3'')')); %#ok<SPWRN>
     [constit,out1]=fe_mat('plil of_mk 2d',ID,mat,il);
     out1(3)=prod(RunOpt.Dim(1:2));out1(4)=RunOpt.Dim(2); 
     % elmap is now done in the 'constants' call
     out=constit(:); out2=[]; 
     return;  % return Constit, Integ, Elmap
    % #BuildElastic2D Standard isotropic fe_mat('pu','m_elastic',1) -3
    elseif strcmp(stt,'m_elastic.1')
        
     %if o2(8)==0 o2(8)=o2(7);end 
     if size(mat,2)<7; mat(7)=0;end;constit = mat([5 7]); % [rho eta]
     E =mat(3); nu=mat(4); ID(3:4)=RunOpt.Dim(2)*[2;1];
     out3.ConstitLab={'rho','eta','D11','D21','D31','D12','D22','D32', ...
        'D13','D23','D33'};

     switch pro(3) % Formulations
     case 1 % - - - - - - - - - -  - - - - - - - plane stress
       C=E/(1.-nu^2);C=[C C*nu C 0. 0. C*(1-nu)/2];
       out3.dd=C([1 2 4;2 3 5;4 5 6]);
     case 0 % - - - - - - - - - -  - - - - - - - plane strain
       C=E*(1-nu)/(1+nu)/(1-2*nu);
       C=[C C*nu/(1-nu) C 0. 0. C*(1-2*nu)/2/(1-nu)];
       out3.dd=C([1 2 4;2 3 5;4 5 6]);
     case 2 % - - - - - - - - - -  - - - - - - - axisymetric
       error('Axisymmetry is not supported by the *b element family');
     end % Formulations 
     constit=[constit(:);out3.dd(:)];
     if sp_util('diag')>1; 
      elem0('DiagConst',sprintf('D%i',ID(1)),reshape(constit(3:11),3,3));
     end
   % #Build2DAnisotropic : car = [rho eta E11 ... E33 a1 a2 a3 T0] -3
   elseif strcmp(stt,'m_elastic.4') 
    if length(mat)<14; mat(14)=0.; end; r1=mat([9 14 3:8 10:13]);
    constit=r1([1 2 3 4 6 4 5 7 6 7 8]); ID(3:4)=RunOpt.Dim(2)*[2;1];
    out3.dd=reshape(constit(3:11),3,3);
     if sp_util('diag')>1; 
      elem0('DiagConst',sprintf('D%i',ID(1)),reshape(constit(3:11),3,3));
     end
   % 2D membrane from orthotropic shell properties
   elseif strcmp(stt,'m_elastic.5')
    if length(mat)<19; mat(19)=0.; end; 
    %S=[1/E1(i) -Nu12(i)/E1(i) 0;-Nu12(i)/E1(i) 1/E2(i) 0 ;
    %    0 0  1/G12(i)];
    r1=pinv([1/mat(3) -mat(5)/mat(3) 0;-mat(5)/mat(3) 1/mat(4) 0;
            0 0 1/mat(6)]);
    constit=[mat([7 18]) r1(:)' mat([19 19 19 12])]; % incorrect alpha
    if sp_util('diag')>1; 
      elem0('DiagConst',sprintf('D%i',ID(1)),reshape(constit(3:11),3,3));
    end
    ID(3:4)=RunOpt.Dim(2)*[2;1];
   elseif strcmp(stt,'m_elastic.6')
     % See formula 6.13
      [constit1,integ1,un1,dd]=p_solid('buildconstit',ID,pl, ...
       [1e6 fe_mat('p_solid','SI',1) 0 1 ]); 
      dd=dd.dd; 
      constit=dd([1:2 4:6],[1:2 4:6]);
      for j1=1:2;
       for j2=1:2
        constit(j1,j2)=constit(j1,j2)-dd(j1,3)*dd(j1,3)/dd(3,3);
       end
      end
      out3.ConstitLab=ConstitLab('m_elastic.6');
      out3.dd=constit; 
      constit=[mat([12 17])';constit(:)];
   % #BuildAcoustic2D -3
   elseif strcmp(stt,'m_elastic.2')
    if length(mat)<5; mat(5)=0.; end
    out3.ConstitLab=ConstitLab('m_elastic.2');
    constit = [1/mat(1,3)/mat(1,4)^2 mat(1,5) 1/mat(1,3)]';
    out3.constit=constit;
    if length(ID)<4;ID(4)=RunOpt.Dim(2); end
    ID(3)=ID(4);ID(7)=2; 
    try;if pro(3)==231;ID(7)=231;end;end % handle curved 2D
    if sp_util('diag')>1;
     fprintf('\n 1/rho/C^2 = %g, eta= %g, 1/rho = %g\n',constit); 
    end
   else; RunOpt.Failed=1;
   end
   % Integration rule saved in integ
   if size(pro,1)&&size(pro,2)>=5;  ID(5)=pro(1,5); end  
  %% #BuildSurface ND-1 element (p_solid.3) ---------------------------2
  elseif isequal(RunOpt.ProF,'p_solid')&&RunOpt.pTyp==3
   % Integration rule saved in integ
   RunOpt.warn={};if length(ID)<4;ID(4)=RunOpt.Dim(2); end;ID(5)=pro(1,3);
   constit=[-3 pro(2:end)];constit=constit(:);
   if pro(1,4)==4 % fluid/structure 
     ID(3)=4*ID(4); % 4 dof per node
   elseif any(pro(1,4)==[1 2 3]); ID(3)=3*ID(4); % pressure, vol (3dof per node)
   elseif any(pro(1,4)==[5 6 7]); ID(3)=2*ID(4); % pressure, vol (2dof per node)
       if pro(1,4)==7;sdtw('_ewt','Need document');end
   elseif any(pro(1,4)==8); ID(3)=ID(4); % pressure, vol (2dof per node)
      if size(mat,2)<6;error('Real part of impedance required');end
      constit(5)=1/prod(mat([3 4 6])); % 1/rho/C/(Real(Z)) wall impedance
      if ~isfinite(constit(5));error('Problem');end
   else; sdtw('_err','%i Not a valid edge element formulation',pro(1,4));
   end
  %% #BuildShell -2
  elseif isequal(RunOpt.ProF,'p_shell')
     try;
       varargin{2}(4,:)=RunOpt.Dim(2);
       if nargout>3;
         eval(sprintf('[out,out1,out2,out3]=%s(varargin{:});',RunOpt.ProF));
       else;eval(sprintf('[out,out1,out2]=%s(varargin{:});',RunOpt.ProF));
       end
       return;
     catch;RunOpt.Failed=1;
     end
  %% #BuildBeam -2
  elseif isequal(RunOpt.ProF,'p_beam')||isequal(RunOpt.ProF,'p_bar')
   out2=[];  % #beam element map - - - - - - - - - - - - - - - - -2
   
   if isfield(RunOpt,'ProStack')&&isfield(RunOpt.ProStack,'constit')
        out2=struct;
        out=RunOpt.ProStack.constit; out1=RunOpt.ProStack.integ;
        out2.elmap= ...
         elem0('elmapmat_og',[double(out1(4)) double(out1(3))/double(out1(4))]);
        return
   end
   if ~strcmp(st,'m_elastic')||typ~=1;
       error('m_elastic(1) is required for p_beam');
   end
   if mat(6)==0; mat(6)=mat(3)/2/(1+mat(4)); end
   if length(mat)<9; mat(9)=0;end     
    %if pro(2)~=1 error('Not a supported beam1 property'); end
   if RunOpt.pTyp~=1; pro=p_beam('ConvertTo1',pro);RunOpt.pTyp=1; end
   if length(pro)<10 ;	pro(10)=0; end

  ElemF='p_beam';
  if evalin('caller','exist(''ElemF'',''var'')');ElemF=evalin('caller','ElemF');end
  if isequal(ElemF,'bar1')
   %     [pe,ie]=fe_mat(1,varargin{1},varargin{2},varargin{3});
   %          E*A nu       eta  rho*A         A lump
   % old : constit = [pe(1)*ie(4)  0    pe(3)*ie(4) ie(4) ie(7)];
   constit = [mat(3)*pro(6)  0    mat(5)*pro(6) pro(6) pro(9:10)];
   ID(3:4)=[6;2];RunOpt.DoElmap=0;
  else
   % #ConstitBeam [E nu Rho G eta alpha(6) T0  J I1 I2 A(11) k1 k2 lump NSM] -2
   % NSM']
   constit = [mat(3:9) pro(3:10)];
   ID(3:4)=[12;2];RunOpt.DoElmap=0;
   out3=struct('ConstitLab',{{'E','nu','rho','G','eta','alpha','T0', ...
       'J','I1','I2','A','k1','k2','lump','NSM'}});
  end
  out=constit(:);
 %% #BuildOther : propagate to other subfunctions -2
  elseif isequal(RunOpt.ProF,'m_null')&&~isempty(RunOpt.warn);
      out=sprintf('%s\n',RunOpt.warn{:});return;
  elseif isequal(RunOpt.ProF,'m_null')&&isempty(pro);
      out=sprintf('Missing ProId %i',ID(2));return;
  elseif isfield(RunOpt.ProStack,'EC')  % Formulation contained in Pro
     constit=[-3 pro(2:end)];constit=constit(:);
     out2=struct('eval','Case.GroupInfo{jGroup,8}=elmap.EC;elmap=[];', ...
           'EC',RunOpt.ProStack.EC);
  else; RunOpt.Failed=1; % 2D, 3D surface
  end

  if RunOpt.Failed % Attempts to call other p_* functions
    try; % Should be in p_* but if user wants to limit foot print p_* OK too
     if ~any(strcmp(RunOpt.ProF,{'p_solid','p_null'}))
      eval(sprintf('[out,out1,out2]=%s(varargin{:});',RunOpt.ProF));
     elseif ~any(strcmp(st,{'m_elastic'}))
      if nargout<4;eval(sprintf('[out,out1,out2]=%s(varargin{:});',st));
      else;eval(sprintf('[out,out1,out2,out3]=%s(varargin{:});',st));
      end
     else; 
         dbstack
         RunOpt.warn{end+1}=sprintf('%s.%i failed',RunOpt.ProF,RunOpt.pTyp);
         error('p_solid error');
     end
    catch
     if nargout<4;st1='[out,out1,out2]';else; st1='[out,out1,out2,out3]';end
     try; 
         st(1)='p';
         eval(sprintf('%s=%s(varargin{:});',st1,st));% try p_ if m_ failed
     catch;
       if sp_util('diag')>10; 
        eval(sprintf('%s=%s(varargin{:});',st1,st));
       end
       error( ...
        horzcat('Call to %s(''BuildConstit'') failed', ...
         ' for subtype (%i) with \n%s\nUse sdtdef(''diag'',12) to debug '), ...
         st,typ,lasterr); %#ok<LERR>
     end
    end
    if length(out1)<4; out1(4,1)=0; elseif out1(4,1)==0; 
    else
      out2=elem0('elmapmat_og', ...
      [double(out1(4)) double(out1(3))/double(out1(4))]); 
    end
    return;
  end% Attempts to call other p_* functions
  if any(size(constit)==1); out=constit(:); else; out=constit;end
  out1=int32(ID(:)); % integ
  if length(out1)<4; out1(4,1)=0; 
  elseif out1(4,1)==0;
    if RunOpt.DoElmap;
        error('Expecting non-zero number of DOF/Node in integ(3:4)');
    end
  elseif isstruct(out2)
   out2.elmap= ...
    elem0('elmapmat_og',[double(out1(4)) double(out1(3))/double(out1(4))]);
  elseif RunOpt.DoElmap % Sort : node/field
   if out1(4)==0; error('Expecting non-zero number of nodes');end
   out2=elem0('elmapmat_og',[double(out1(4)) double(out1(3))/double(out1(4))]);
  end
  if ~isempty(RunOpt.warn);sdtw('_nb','%s\n',RunOpt.warn{:});end

% -------------------------------------------------------------------------
elseif comstr(Cam,'builddof')
%% #BuildDof Advanced multi-physic dof building

  model=varargin{carg};carg=carg+1;
  cEGI=varargin{carg};carg=carg+1;
  nd=varargin{carg};carg=carg+1;
  ElemF=varargin{carg};carg=carg+1;
  RunOpt=struct('FieldDofs',[],'pro',fe_super('prop',ElemF), ...
   'node',fe_super('node',ElemF),'PropertyDof',[],'warn',{{}}, ...
   'OldFieldDofs',[]);
  i4=[];
 %proid
  if ~isfield(model,'pl');   model.pl=zeros(0,1); end
  if ~isfield(model,'il');   model.il=zeros(0,1); end
  RunOpt.EGID=0;
  try; [st,RunOpt.EGID]=feutil('getelemf',model.Elt(cEGI(1)-1,:),1);end %#ok<TRYNC>
  i1=unique(model.Elt(cEGI,RunOpt.pro(RunOpt.pro(1:2)~=0)),'rows');
  for j1=1:size(i1,1) % loop on MatId ProId pairs
   if isempty(model.pl);pl=[];
   else; pl=model.pl(model.pl(:,1)==i1(j1,1),:);
   end
   if isempty(pl); [pl,mat]=fe_mat(sprintf('getpl %i',i1(j1,1)),model);end
   if isempty(pl);
      typ='m_elastic';st='SI';i2=1;
      RunOpt.missPl=sprintf( ...
           'model.pl should contain MatId %i for BuildDof',i1(j1,1));
   else; [typ,st,i2]=eM.tomType(pl(1,2));RunOpt.missPl='';
   end
   if isempty(typ); error('Missing property type in pl=%s', ...
           comstr(pl,-30));
   end
   % default if undefined is a p_solid for backward compat
   il=model.il; pro=[];
   if ~isempty(il); il=il(il(:,1)==i1(j1,2),:);end
   if isempty(il); [il,pro]=fe_mat(sprintf('getil %i',i1(j1,2)),model);end
   if isempty(il)&&any(strcmp(ElemF,{'quad4','tria3','quadb','tria6', ...
    'beam1','beam1t'}))
      RunOpt.warn{end+1}=sprintf( ...
       'model.il should contain ProId %i for BuildDof',i1(j1,2));
     RunOpt.pFcn='p_shell';RunOpt.pUnit='SI';RunOpt.pTyp=1;
   elseif isempty(il)
     RunOpt.pFcn='p_solid';RunOpt.pUnit='SI';RunOpt.pTyp=0;
     RunOpt.warn{end+1}=sprintf( ...
       'model.il should contain ProId %i for BuildDof',i1(j1,2));
   else; [RunOpt.pFcn,RunOpt.pUnit,RunOpt.pTyp]=eM.topType(il(1,2));
   end
   RunOpt.Cb='';
   switch RunOpt.pFcn % property type
   case 'p_solid'
     RunOpt.FieldDofs=1:3; % this is the default
     if RunOpt.pTyp==3 
       %% #BuildDof.Surface elements for coupling
       if il(1,4)==4;  RunOpt.FieldDofs=[1:3 19]; % fluid structure coupling
       elseif any(il(1,4)==[1 2 3]);RunOpt.FieldDofs=1:3; % pressure 3D
       elseif any(il(1,4)==[5 6 7]);RunOpt.FieldDofs=1:2; % pressure 2D
       elseif any(il(1,4)==8);RunOpt.FieldDofs=19; % wall impedance
       else; sdtw('_nb','p_solid(3) type %i not supported',il(1,4));
       end
       RunOpt.warn={};
     elseif isequal(typ,'m_elastic')
      %% #BuildDof.m_elastic

      if i2==2; RunOpt.FieldDofs=19;
      elseif i2==4||RunOpt.pTyp==2; RunOpt.FieldDofs=[1 2];  % Plane elements
      elseif ~isempty(il)&&any(il(1,4)==[20002 30003]) % u-p hyperelastic
        RunOpt.Cb=m_hyper('@elasUP');
      else% if any(i2)==[1 3]; 
        RunOpt.FieldDofs=1:3;
      end
     elseif isequal(typ,'m_piezo');
       RunOpt.FieldDofs=[1 2 3 21];
     else; 
      try;  RunOpt.FieldDofs=feval(typ,'FieldDofs',RunOpt,i2);
      catch;
        sdtw('(p=%s and m=%s) is unknown',RunOpt.pFcn,typ);
      end
     end 
   % This calls a user p_*.m function. The call is supposed to fill in the
   % RunOpt.FieldDofs based on the contents of pl and il
   case 'p_beam'
      if strncmpi(ElemF,'bar',3); RunOpt.FieldDofs=[1 2 3];
      else;RunOpt.FieldDofs=[1 2 3 4 5 6];
      end
   case 'p_bar'
       RunOpt.FieldDofs=[1 2 3];
   otherwise 
      if isfield(pro,'EC');RunOpt.FieldDofs=pro.EC.FieldDofs;
      else
       try; 
        RunOpt.warn={};
        RunOpt=feval(RunOpt.pFcn,'BuildDofOpt',RunOpt,model.pl,il,model);
       catch
        if ~isempty(RunOpt.missPl)
            RunOpt.warn{end+1}=RunOpt.missPl;
        end
       end
      end
   end
   if ~isempty(RunOpt.Cb); RunOpt=feval(RunOpt.Cb,'dof',RunOpt,model,cEGI);end
   if isempty(RunOpt.OldFieldDofs);RunOpt.OldFieldDofs=RunOpt.FieldDofs;
   elseif ~isequal(RunOpt.FieldDofs(:),RunOpt.OldFieldDofs(:))
      error('Field DOFs cannot be variable within the same group %s', ...
        sprintf('MatId%i ProId%i \n use feutilb(''SeparatebyProp'')', ...
        i1(j1,1:2)))
   end
  end % Loop on Pairs
  if ~isempty(RunOpt.warn)&&RunOpt.EGID(1)>0;
    sdtw('_nb','%s\n',RunOpt.warn{:});
  end
  
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  % Below here is a common to all element functions - - - - - - - -
  RunOpt.FieldDofs=RunOpt.FieldDofs(:);
  if isempty(nd); % build the DOF list for use in feutil('getdof')
   i2=unique(model.Elt(cEGI,RunOpt.node))*100;i2=i2(:);
   i3=RunOpt.FieldDofs(:)';
   i2=(i2(:,ones(size(i3)))+i3(ones(size(i2)),:))';i2=i2(:);
   out=[i2;round(RunOpt.PropertyDof*100)]; % Nodal DOF + Property(/variable) DOF
   out1=[];% return nodal and element DOFs
   evalin('caller',sprintf( ...
    'out2=out2+length(cEGI)*%i^2;', ...
    length(RunOpt.FieldDofs)*length(RunOpt.node)+length(RunOpt.PropertyDof)));
  else % build the DofPos for use in fe_mknl
   Case=evalin('caller','Case');
   %    DofPos=reshape(int32(full(nd(RunOpt.FieldDofs(:), model.Elt(cEGI,RunOpt.node)'))-1), ...
   %       length(RunOpt.FieldDofs)*length(RunOpt.node),length(cEGI));
   if isempty(RunOpt.FieldDofs)&&isfield(RunOpt,'VariableDof')
    DofPos=reshape(full(nd.NNode(model.Elt(cEGI,RunOpt.node)')),[],length(cEGI));
    DofPos=DofPos(fix(RunOpt.VariableDof),:)+repmat(rem(RunOpt.VariableDof,1),1,length(cEGI));
    DofPos=int32(full(nd.nd(round(DofPos*100)))-1);out=DofPos; out1=[]; return;
   else
    DofPos=reshape(int32(feval(nd.getPosFcn,nd,model.Elt(cEGI,RunOpt.node)',...
     RunOpt.FieldDofs(:))-1), ...
     length(RunOpt.FieldDofs)*length(RunOpt.node),length(cEGI));
   end
   if ~isempty(RunOpt.PropertyDof); % 
    i4=RunOpt.PropertyDof; % i4=int32(full(nd(21,fix(i4)))-1);
    i4=int32(feval(nd.getPosFcn,nd,fix(i4),21)-1); i4=i4(:);
    DofPos(end+(1:length(i4)),:)=i4(:,ones(size(DofPos,2),1));
   end
   out=DofPos;out1=[];
  end % Build DOF list or DofPos


%% #End ----------------------------------------------------------------------
elseif comstr(Cam,'test');out='';
elseif comstr(Cam,'tablecall');out='';
elseif comstr(Cam,'@');out=eval(CAM);
elseif comstr(Cam,'cvs');
 out='$Revision: 1.267 $  $Date: 2022/12/01 17:48:32 $'; return;
else; sdtw('''%s'' not known',CAM);
end

%% #ctableGen EC=ctableGen(EC,Case.GroupInfo{Case.jGroup,7},ConstitLab(RO.mtt))
% xxx D11(T) with T stored in InfoAtNode
function EC=ctableGen(EC,IA,clab,RA)
 if ischar(clab); clab=ConstitLab(clab);end
     
 if nargin==4
   %EC=feval(p_solid('@ctableGen'),EC,InfoAtNode,clab,mat)
   % define nodeE column and offset for each interpolated field
   % 0 in third column is ref (0)
     r3=find(strcmpi(IA.lab,'T'))+4;
     constit=RA.constit;mat=RA.mat;
     if strcmpi(RA.type,'nonlin_elas_iso')
      RA.list={'E',r3,0;'Nu',r3,0;'G',r3,0;
         'Alpha',r3,0;'Rho',r3,0}; 
      [st,i1,i2]=fe_mat('type',constit(2));
      st1=feval(mat.type,'propertyunittype cell',i2);st1=st1(:,1);
     elseif strcmpi(RA.type,'nonlin_elas_aniso')
      RA.list=clab(:);RA.list(:,2)={r3}; RA.list(:,3)={0};i2=3;
      st1=clab(:);
     else; error('Not implemented');
     end

     r3=RA.list;
     r1=cell(1); 
     for j1=1:size(st1,1);st1{j1}=sscanf(st1{j1},'%s',1);end     
     st=fieldnames(mat);
     for j1=1:length(st)
        r2=mat.(st{j1});
        if isfield(r2,'X')&&isfield(r2,'Y') % this is a table
           j3=find(strcmp(r3(:,1),st{j1}));
           j4=find(strcmp(st1,st{j1})); 
           if isempty(j4);sdtw('_nb','Field %s not found in constit',st{j1});end
           r1{j1,1}=[0 0 0 0 length(r2.X) r3{j3,2} j4];
           r1{j1,2}=[r2.X(:)-r3{j3,3};r2.Y(:)]; % [X-ref;Y]
        else; st{j1}='';
        end
     end
     st(cellfun('isempty',st))=[]; 
     fprintf('Interpolating MatId %i ( %s) as functions of T\n',mat.pl(1), ...
           sprintf('%s ',st{:}));
     i1=find(~cellfun('isempty',r1(:,1)));
     r2=vertcat(r1{i1,1});i2=numel(r2)+1;
     for j1=1:size(r2,1); r2(j1,4)=i2;i2=i2+numel(r1{i1(j1),2});end
     EC.CTable=[size(r2,1);reshape(r2',[],1);vertcat(r1{i1,2})];
     % #CTable constit tables are stored using 
     %[Ntables
     % Current_values (7 x Ntables) giving 
     %  [i1 xi si xstartpos Nx nodeEfield constit-pos(1based)]
     % tables ]
     %  xstartpos==-1 : use nodeE field value
     %  [constitoff 0 0  xstartpos==-2 0 0 0 0]: local orient 3x3 constit
     %  [constitoff 0 0  xstartpos==-3 0 0 0 0]: local orient 6x6 mecha

     if strcmpi(RA.type,'nonlin_elas_iso')
      EC.MatrixIntegrationRule{1}=[];
     end
 elseif isfield(IA,'lab')&&length(intersect(IA.lab,clab))>1
 % -1 Map of D11(constit label) stored in InfoAtNode (not a table based interp)
     [i1,i2]=ismember(clab,[{'x';'y';'z';'NodeId'};IA.lab(:)]);
     r2=zeros(nnz(i1),7);r2(:,4)=-1; % No table
     r2(:,6)=i2(i1); % nodeE field
     r2(:,7)=find(i1); % constit-pos(1based)
     EC.CTable=[size(r2,1);reshape(r2',[],1)]; % sdtweb m_elastic('ctable')
     
  if isfield(EC,'nodeE') % Extend nodeE size
     EC.nodeEt(size(EC.nodeE,2)+(1:length(IA.lab)))=comstr(IA.lab,-32);
     EC.nodeEt=int32(EC.nodeEt);
     EC.nodeE(1,length(EC.nodeEt))=0;
  end

 end
   


%% #ConstitLab : constit row labels
% lab=feval(p_solid('@ConstitLab'),'m_elastic.1')
function out=ConstitLab(CAM)

switch lower(CAM)
case {'m_elastic.1','m_elastic.3'} 
 % Constit vector for m_elastic.1 (anistropic solid)
 out={'rho','eta','D11','D21','D31','D41','D51','D61', ...
        'D12','D22','D32','D42','D52','D62', ...
        'D13','D23','D33','D43','D53','D63', ...
        'D14','D24','D34','D44','D54','D64', ...
        'D15','D25','D35','D45','D55','D65', ...
        'D16','D26','D36','D46','D56','D66', ... % 1:38
        'k11','k21','k31','k12','k22','k32','k13','k23','k33', ...
        'T0'};
case 'm_elastic.2';out={'1/rho/C2','eta','1/rho'};
case 'm_elastic.6'; out={'rho','eta','D11','D21','D31','D12','D22','D32', ...
        'D13','D23','D33'};
case {'p_shell.1'} % sdtweb p_shell shellconstit
    out={'rho_h','eta','f','d','h','k','i12t3','nsm', ...
        'D11','D21','D31','D41','D51','D61', ...
        'D12','D22','D32','D42','D52','D62', ...
        'D13','D23','D33','D43','D53','D63', ...
        'D14','D24','D34','D44','D54','D64', ...
        'D15','D25','D35','D45','D55','D65', ...
        'D16','D26','D36','D46','D56','D66', ... % 9:44
        'DS11','DS21','DS12','DS22', ... % 45:48
        'Ddril'};
otherwise;  out={};
end


function  [EC,RO]=EC_Elas2D(EC,RO);     
 %% #EC_Elas2D Strain energy
 % define the deformation vector: row, NDN, DDL, NwStart, NwRule
 rule=RO.rule;
 EC.StrainDefinition{1}= ...
   [1 2 1 rule; 2 3 2 rule; % e_xx = N,x u, e_yy
    3 3 1 rule;3 2 2 rule]; % shear
 EC.StrainLabels{1}=p_solid('ConvStrain2D');
 EC.ConstitTopology{1}=int32(reshape(3:11,3,3));
 % fprintf('D=\n');disp(constit(EC.ConstitTopology{1}));
 EC.DofLabels={'u','v'};     Ndof=2*EC.Nnode;
 RO.NdnDim=2; % Tell that BuildNDN rule is 2D
 % Kinetic energy
 EC.StrainDefinition{2}= [1 1 1 rule;2 1 2 rule];
 EC.StrainLabels{2}={'u','v'};
 EC.ConstitTopology{2}=int32(eye(2));
 EC=integrules('matrixrule',EC);
 % \int f v (see elem0 for format)
 EC.RhsDefinition=int32( ...
   [101 0 2 0     0 0 -1    rule+[-1 0];
    101 1 2 0     1 0 -1    rule+[-1 0]]);
 %display with :integrules('texstrain',EC);
 EC.VectMap=int32(reshape(1:2*EC.Nnode,2,EC.Nnode)'); 
 
function [EC,RO]=EC_Elas3D(EC,RO,integ,constit);      
 % #EC_Elas3D #Elast3dStrainDef Strain energy -2
 % define the deformation vector: row, NDN, DDL, NwStart, NwRule
 rule=RO.rule;
 EC.StrainDefinition{1}= ...
   [1 2 1 rule; 2 3 2 rule;3 4 3 rule; % e_xx = N,x u, e_yy, e_zz
   4 4 2 rule;4 3 3 rule;5 4 1 rule;5 2 3 rule;6 3 1 rule;6 2 2 rule]; % shear
 EC.StrainLabels{1}=p_solid('ConvStrain3D'); RO.mtt='m_elastic.1';
 EC.ConstitTopology{1}=int32(reshape(3:38,6,6));
 RO.NdnDim=3; % Tell that BuildNDN rule is 3D

 % Kinetic energy
 EC.StrainDefinition{2}= [1 1 1 rule;2 1 2 rule; 3 1 3 rule];
 EC.StrainLabels{2}={'u','v','w'};
 if size(constit,1)<20 % reinterpolated material
   try;i1=fe_mat('getpos',constit(2,1),'rho'); %sdtweb elem0('rhs_og')
   catch; i1=1;
   end
   if ~isempty(i1)
     EC.DensPos=i1;EC.ConstitTopology{2}=int32(eye(3)*EC.DensPos);
   else;EC.ConstitTopology{2}=int32(eye(3));
   end
 else; EC.ConstitTopology{2}=int32(eye(3));
 end
 EC=integrules('matrixrule',EC);
 % \int f v (see elem0 for format)
 EC.RhsDefinition=int32( ...
   [101 0 3 0     0 0 -1    rule+[-1 0];
    101 1 3 0     1 0 -1    rule+[-1 0];
    101 2 3 0     2 0 -1    rule+[-1 0]]);
 %display with :integrules('texstrain',EC);

 EC.VectMap=int32(reshape(1:3*EC.Nnode,3,EC.Nnode)'); 
 EC.material='Elastic3DNL';
  if isfield(RO,'GstateFcn')&&~isempty(RO.GstateFcn)&& ...
         ~comstr(RO.GstateFcn,'Case=')
    EC.gstate=RO.GstateFcn;
  end

 % Switch on NL constitutive materials - - - - - - - - - - - - - - -
 if size(integ,1)>6; 
 switch integ(7,1) 
 % this should match StrategyType in of_mk.c matrix assembly
 case 3; 
   %% Tested in openfem/demos/RivlinCube % probably completely obsolete in 2022
   if isempty(EC.ConstitTopology{1})||length(constit)<38
    EC.material='hyperelastic';
    EC.ConstitTopology{1}=[];EC.MatrixIntegrationRule{1}=[];
    EC.ConstitTopology{2}=int32(eye(3));EC.MatrixIntegrationRule{2}(:,5)=0;
    RO.NdnDim = 31;
   end
 case 105;
   EC.material='pre_stress';
   RO.NdnDim=31;
 case 901; 
   EC.material='heart';
   cEGI=evalin('caller','varargin{end-1}');
   gstate=zeros(7+EC.Nnode*length(EC.DofLabels)+1,EC.Nw*length(cEGI));
   EC.gstate = gstate;
   RO.NdnDim = 31;
 end;
 end
 

function [EC,RO]=EC_Elas3Dld(EC,RO,integ,constit);
 % #EC_Elas3Dld large deformation using ui,j strain -2
 rule=RO.rule;
 EC.StrainDefinition{1}= ...
   [1 2 1 rule; 2 3 1 rule;3 4 1 rule; % F11 = N,x u, F21=N,y u 
    4 2 2 rule; 5 3 2 rule;6 4 2 rule;
    7 2 3 rule; 8 3 3 rule;9 4 3 rule];
 EC.StrainLabels{1}={'u1,1','u1,2','u1,3','u2,1','u2,2','u2,3','u3,1','u3,2','u3,3'};
 i1=int32(reshape(3:38,6,6));
 ind_ts_eg=[1 6 5 6 2 4 5 4 3];
 EC.ConstitTopology{1}=int32(i1(ind_ts_eg,ind_ts_eg));

 % Kinetic energy
 EC.StrainDefinition{2}= [1 1 1 rule;2 1 2 rule; 3 1 3 rule];
 EC.StrainLabels{2}={'u','v','w'};EC.ConstitTopology{2}=int32(eye(3));
 EC=integrules('matrixrule',EC);
 % \int f v (see elem0 for format)
 EC.RhsDefinition=int32( ...
   [101 0 3 0     0 0 -1    rule+[-1 0];
    101 1 3 0     1 0 -1    rule+[-1 0];
    101 2 3 0     2 0 -1    rule+[-1 0]]);
 EC.VectMap=int32(reshape(1:3*EC.Nnode,3,EC.Nnode)'); 
 RunOpt=RO.RunOpt;
 if isfield(RunOpt,'GstateFcn')&&~isempty(RunOpt.GstateFcn)&& ...
         ~comstr(RunOpt.GstateFcn,'Case=')
    EC.gstate=RunOpt.GstateFcn;
 end
 if ~isfield(RunOpt,'EltOrient')
 elseif size(RunOpt.EltOrient,1)>0&&isfield(RunOpt.EltOrient{1,3},'NLdata');
     EC.NLdata=RunOpt.EltOrient{1,3}.NLdata;
     if ~isfield(EC.NLdata,'constit'); EC.NLdata.constit=constit;end
 end
 EC.material='gennl';
 RO.NdnDim=3; % Tell that BuildNDN rule is 3D


