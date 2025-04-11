function [out,out1,out2,out3]=p_heat(varargin)

%P_HEAT property function of heat and other diffusion equations
%
%       Syntax : pro= p_heat('default') 
%                pro= p_heat('database matid (command)') 
%                il = p_heat('dbval ProId (command)');
%                il = p_heat('dbval ProId -unit MM (command)');
%                il = p_heat('dbval ProId -punit MM (command)');
%
%       Supported ELEMENT property subtypes are
%       1 : volume element for heat diffusion (dimension DIM)
%           [ProId fe_mat('p_heat','SI',1) CordM Integ DIM] 
%       2 : surface element for heat exchange (dimension DIM-1)
%           [ProId fe_mat('p_heat','SI',2) CordM Integ DIM] 
%         Integ : is rule number in integrules
%         DIM is problem dimension 2 or 3 D
%       3 : thermoelastic coupling : UNDER DEVELOPMENT
%           [ProId fe_mat('p_heat','SI',3) ElasMatId ThermoMatId Integ T0] 
%
%
% model=p_heat('SetFace',model,Selelt,pl,il);
%
%       See <a href="matlab: sdtweb _taglist p_heat">TagList</a>
%       See sdtweb     fem (handling materials section), pl, fe_mat
%       See also help  fe_mat m_heat

%	Etienne Balmes, Frederic Bourquin, A. Nassiopoulos, Jean-Philippe Bianchi
%       Copyright (c) 2001-2025 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

if nargin<1; help p_heat;return; end
if ischar(varargin{1}); [CAM,Cam]=comstr(varargin{1},1); il=[]; carg=2;
else il=varargin{1};[CAM,Cam]=comstr(varargin{2},1); carg=3;
end
%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*NOSEM>

%% #Default ------------------------------------------------------------------
if comstr(Cam,'default')

  out=p_heat('database');
  out=out(1);

%% #dbval --------------------------------------------------------------------
elseif comstr(Cam,'dbval')

 while 1==1
  i1=strfind(comstr(Cam,-27),'-unit'); out1={};
  if ~isempty(i1)
   [Unit,i2,i3,i4]=sscanf(CAM(i1+5:end),'%s',1);
   i4=i1+[0:4+i4];CAM(i4)=''; [CAM,Cam]=comstr(CAM,1);
  else Unit='';
  end
  i2=strfind(comstr(Cam,-27),'-punit');
  if ~isempty(i2)
   [PUnit,i3,i4,i5]=sscanf(CAM(i2+6:end),'%s',1);
   i5=i2+[0:5+i5];CAM(i5)=''; [CAM,Cam]=comstr(CAM,1);
  else PUnit='';
  end
  [CAM,Cam]=comstr(CAM,6);
  [i1,CAM,Cam]=comstr(CAM,'','%i');
  if isempty(CAM)&&carg<=nargin; st=varargin{carg};carg=carg+1;
  else st=CAM;end
  [mat,st1,i2]=p_heat('database',st,varargin{carg:end});carg=carg+i2-3;
  if ~isempty(PUnit)
   r1=fe_mat(sprintf('convert %s %s',mat.unit,PUnit),mat.il(1:2));
   mat.il(2)=r1(2); mat.unit=PUnit;
  end
  if ~isempty(Unit)
   mat.il=fe_mat(sprintf('convert %s %s',mat.unit,Unit),mat.il);mat.unit=Unit;
  end
  r1=mat.il; if length(i1)==1; r1(1)=i1;end
  if ~isempty(il); i2=find(il(:,1)==r1(1)); else i2=[];end
  if isempty(i2); i2=size(il,1)+1;end
  il(i2,1:length(r1))=r1;
  if carg>nargin; break;end
  CAM=varargin{carg};carg=carg+1;if ischar(CAM);[CAM,Cam]=comstr(CAM,1);end
 end
 out=il;

%% #Database -----------------------------------------------------------------
elseif comstr(Cam,'database')

  st=comstr(CAM,9);
  if isempty(st)&&carg<=nargin; st=varargin{carg}; carg=carg+1;end
  [MatId,i2,i3,i4]=sscanf(st,'%i',1); if i2~=1; MatId=1;end
  st=comstr(st,i4);

  out.il=[MatId fe_mat('p_heat','SI',1) 0 2 3];
    out.name='default heat (volume diffusion)';
    out.type='p_heat';
    out.unit='SI';

  out(end+1).il=[MatId fe_mat('p_heat','SI',1) 0 2 3];
    out(end).name='d3 (diffusion)';
    out(end).type='p_heat'; out(end).unit='SI';

  out(end+1).il=[MatId fe_mat('p_heat','SI',2) 0 2 3];
    out(end).name='d3 (surface exchange)';
    out(end).type='p_heat'; out(end).unit='SI';
    
  out1='Heat properties';out2=carg;


  i1=find(strncmpi(st,{out.name},length(st)));
  st0={'diffusion';'surface exchange'};
  if ~isempty(i1);out=out(i1);
  elseif comstr(comstr(st,-27),'d3') % DIM=3
   %[ProId fe_mat('p_heat','SI',...) CoordM Integ DIM] 
   r1=comstr(st(3:end),[-1 3]);
   i0=1; % subtype
   if isempty(r1); r1=-3; end
   if length(r1)==2; i0=r1(2); r1=r1(1); end
   if length(r1)==1; r1=[0 r1 3]; end
   out=struct('il',[MatId fe_mat('p_heat','SI',i0) r1(:)'], ...
              'name',horzcat('heat',st,' ',st0{i0}),'type','p_heat','unit','SI');
  elseif comstr(comstr(st,-27),'d2') % DIM+2
   %[ProId fe_mat('p_heat','SI',...) CoordM Integ DIM] 
   r1=comstr(st(3:end),[-1 3]);
   i0=1; % subtype
   if isempty(r1); r1=-3; end
   if length(r1)==2; i0=r1(2); r1=r1(1); end
   if length(r1)==1; r1=[0 r1 2];end
   out=struct('il',[MatId fe_mat('p_heat','SI',i0) r1(:)'], ...
         'name',horzcat('heat',st,' ',st0{i0}),'type','p_heat','unit','SI');
%   elseif ~isempty(st)
%    r1=str2num(st);  %#ok<ST2NM>
%    if length(r1)<1; error('Not a consistent database call');end
  elseif comstr(lower(st),'thermoelas') % DIM+2
   r1=comstr(st(11:end),-1);
   out=struct('il',[MatId fe_mat('p_heat','SI',3) r1(:)'], ...
         'name',sprintf('ThermoElas%i',MatId),'type','p_heat','unit','SI');
  end
  
% -------------------------------------------------------------------------
%% #BuildConstit [constit,integ]=p_heat('buildconstit',ID,pl,il);
elseif comstr(Cam,'buildconstit')

ID=varargin{carg};carg=carg+1;
pl=varargin{carg};carg=carg+1;
il=varargin{carg};carg=carg+1;il0=il;
il=il(il(:,1)==ID(2),:);
RunOpt=evalin('caller','RunOpt');
[RunOpt.ProF,RunOpt.ProU,RunOpt.pTyp]=fe_mat('typep',il(1,2));
out3=struct;
switch RunOpt.pTyp 
case {1,2} % Volume/surface
 pl=pl(pl(:,1)==ID(1),:); if size(pl,2)<6;pl(1,6)=0;end
 % Here one defines how a CONSTIT vector is specified - - - - -
 if il(1,5)==3; K=diag(pl(1,3)*[1 1 1]);
 else; K=diag(pl(1,3)*[1 1]);
 end
 % -1 typ     rho*C           alpha    K  
 out3.ConstitLab={'-1','typ','rhoC','alpha', ...
     'k11','k21','k31','k12','k22','k32','k13','k23','k33'};
 constit=[-1;pl(1,2);pl(1,4)*pl(1,5);pl(1,6);K(:)];
 % xxx For composites, anisotropic alpha should be possible

 [st,i1,i2]=fe_mat('typep',il(1,2));
 ID=ID(:); ID(3:4,1)=0; % only one DOF per node
 ID(5:7)=[il(1,4);il(1,5);i2];  % Integrule, DIM, subtype
case 3 % Thermoelastic coupling
 model=varargin{5};Case=varargin{6};
 % thermo
 [r4,i2,i4,r2]=p_solid(Cam,[il(4);1],pl, ...
     [1 fe_mat('p_solid','SI',1) 0 -3 3],model,Case);
 % elas call to build constit
 [r4,i3,i4,r3]=p_solid(Cam,[il(3);1],pl, ...
     [1 fe_mat('p_solid','SI',1) 0 -3 0],model,Case);
 %Compute C:alpha 
 if ~isfield(r3,'at'); 
     error('Thermal expansion not defined for MatId %i',il(3))
 else; r3.at=r3.at([1 5 9 8 7 4]');% to Eng : xx yy zz yz zx xy
 end
 Dalpha=r3.dd*r3.at;
 out3.ConstitLab={'-1','typ','rhoC','T0','da1','da2','da3','da4','da5','da6', ...
     'ad1','ad2','ad3','ad4','ad5','ad6'};
 if length(il)<6; il(6)=0; warning('Set T0 at 0 C');end
     
 constit=[-1;il(1,2);0;il(6);Dalpha(:);Dalpha(:)*(il(6)+273.15)];
 ID=ID(:); ID(3:4)=[0]; 
 ID(3:4)=RunOpt.Dim(2)*[4;1]; % 4 DOF per node
otherwise; error('p_heat stype=%i not implemented',RunOpt.pTyp)
end
out=constit(:);     % Return CONSTIT
out1=int32(ID(:));  % Return integ
out2=[]; % elmap not used

%% #Const EltConst=p_solid('ConstSolid',ElemF,integ,constit); ----------------
elseif comstr(Cam,'const')

 EC=varargin{carg};carg=carg+1;
 integ=varargin{carg};carg=carg+1;
 if carg<=nargin; constit=varargin{carg};carg=carg+1;
 else constit=[];
 end

 if ischar(EC) % Allow for integrule switching here
  if size(integ,1)>4; 
    EC=integrules(sscanf(EC,'%s',1),double(integ(5,1)));
    if any(integ(5,:)~=integ(5,1));
     warning('OpenFEM:INTEG',...
       'Cannot deal with multiple integration strategies in the same group');
    end
  else; EC=integrules(EC);
  end
 end
 % Adjust Number of Dof and Node in integ if needed
 if size(integ,1)>3&&~any(integ(3:4,1))% back propagate number of DOFs if needed
     if ~any(EC.xi(:,3)); i1=2; elseif norm(EC.xi(:,2:3),'inf')==0; i1=1; 
     else; i1=3; end
     if size(integ,1)<4; error('Improper init');end
     sp_util('setinput',integ, ...
         [EC.Nnode;EC.Nnode],int32(3:4),int32(1));
     EC.ElMap=elem0('elmapmat_og',[EC.Nnode 1]);
 end
 out2=[]; % Possibly InfoAtNode init
 
 EC.Nw=size(EC.N,1);EC.Nnode=size(EC.N,2);
 if size(integ,1)<9;rule=[1 EC.Nw]; % By default start at 1 use all points
 else; rule=integ(8:9,1);rule=rule(:)';
 end
 if ~any(rule); 
   rule=[1 EC.Nw];
   sdtw('integ is assumed to give standard integration rule in integ(5:6)');
 end
 EC.DofLabels={'T'}; 
% Thermoelastic coupling
if size(integ,1)==4&&fe_mat('typesub',constit(2))==3
 EC.DofLabels={'u','v','w','T'};
 % Strain energy N,i (C:alpha) N
 % sdtweb p_solid('Elast3dStrainDef')
 % ep  N,i N : row, NDN, DDL, NwStart NwRule
 EC.StrainDefinition{1}= ...
   [1 2 1 rule; 2 3 2 rule;3 4 3 rule; % e_xx = N,x u, e_yy, e_zz
   4 4 2 rule;4 3 3 rule;5 4 1 rule;5 2 3 rule;6 3 1 rule;6 2 2 rule; % shear
   7 1 4 rule]; 
 EC.StrainLabels{1}=p_solid('ConvStrain3D');EC.StrainLabels{1}{end+1}='T';
 r1=zeros(7,7); r1(1:6,7)=4+(1:6); EC.ConstitTopology{1}=int32(r1);
 
 % Damping (K_theta u)
 EC.StrainDefinition{3}=EC.StrainDefinition{1};
 EC.StrainLabels{3}=EC.StrainLabels{1}; 
 r1=zeros(7,7); r1(7,1:6)=10+(1:6); EC.ConstitTopology{3}=int32(r1);
 
 EC.RhsDefinition=[]; % \int f v (see elem0 for format)
 if sp_util('diag')>10
   dd=full(EC.ConstitTopology{1});dd(dd~=0)=constit(dd(dd~=0));
   disp(dd);
 end
 
 out1=3;

% 3D \int KgradT gradv - - - - - - - - - - - - - - - - - -
elseif isequal(double(integ(6:7,1)),[3;1])
  % Strain energy
  % define the deformation vector: row, NDN, DDL, NwStart NwRule
  EC.StrainDefinition{1}=[1 2 1 rule;2 3 1 rule;3 4 1 rule];
  EC.StrainLabels{1}={'T,x','T,y','T,z'};
  EC.ConstitTopology{1}=int32(reshape(4+[1:9],3,3));
  % Kinetic energy
  EC.StrainDefinition{2}= [1 1 1 rule];
  EC.StrainLabels{2}={'T'}; EC.ConstitTopology{2}=int32(3);
  % \int f v (see elem0 for format)
  EC.RhsDefinition=int32([101 0 1 0     0 0 -1    rule+[-1 0]]);
  out1=3;
 % 2D \int KgradT gradv - - - - - - - - - - - - - - - - - -
 elseif isequal(double(integ(6:7,1)),[2;1]) 
  % Strain energy
  % define the deformation vector: row, NDN, DDL, NwStart NwRule
  EC.StrainDefinition{1}=[1 2 1 rule;2 3 1 rule];
  EC.StrainLabels{1}={'T,x','T,y'};
  EC.ConstitTopology{1}=int32(reshape(4+[1:4],2,2));
  % Kinetic energy
  EC.StrainDefinition{2}= [1 1 1 rule];
  EC.StrainLabels{2}={'T'}; EC.ConstitTopology{2}=int32(3);
  % \int f v
  EC.RhsDefinition=int32([101 0 1 0     0 0 -1    rule+[-1 0]]);
  %EC.Rhs
  out1=2;
 % 3D \int_{\partial \omega} \alpha T v
 elseif isequal(double(integ(6:7,1)),[3;2]) 
  EC.StrainDefinition{1}= [1 1 1 rule];EC.StrainDefinition{2}=[];
  EC.StrainLabels{1}={'T'};
  EC.ConstitTopology{1}=int32(eye(2)*4);
  % \int_{\partial \omega} (g+\alpha \theta_ext) v
  EC.RhsDefinition=int32([101 0 1 0     0 0 -1    rule+[-1 0]]);
  out1=23;
 % 2D \int_{\partial \omega} \alpha T v
 elseif isequal(double(integ(6:7,1)),[2;2]) 
  EC.StrainDefinition{1}= [1 1 1 rule];EC.StrainDefinition{2}=[];
  EC.StrainLabels{1}={'T'};
  EC.ConstitTopology{1}=int32(4);
  % \int_{\partial \omega} g v
  %101(1) InfoAtNode1(2) InStep(3) NDNOff1(4)   FDof1(5) NDNCol(6)
  %NormalComp(7) w1(8) nwStep(9)
  EC.RhsDefinition=int32([101 0 1 0     0 0 -1    rule+[-1 0]]);
  out1=13;
 else; error('Not a supported case');
 end
 EC.VectMap=int32(reshape(1:EC.Nnode*length(EC.DofLabels), ...
   length(EC.DofLabels),EC.Nnode)');
 EC=integrules('matrixrule',EC);
 %EC.material='multi';

 if ~isempty(constit)
  % Only keep terms needed for the current law.
  i1=find(any(constit,2))-1; 
  for j0=1:length(EC.MatrixIntegrationRule);
    r1=EC.MatrixIntegrationRule{j0}; 
    if ~isempty(r1)
    EC.MatrixIntegrationRule{j0}=r1(ismember(double(r1(:,5)),i1),:);
    end
  end
 end
 if nargout==0
  integrules('texstrain',EC);
  try; EC=integrules('stressrule',EC);integrules('texstress',EC);end
 else; out=EC;
 end
% -------------------------------------------------------------------------
%% #BuildDofOpt RunOpt=p_heat('BuildDofOpt',RunOpt,pl,il);
elseif comstr(Cam,'builddofopt')

RunOpt=varargin{carg};carg=carg+1;
pl=varargin{carg};carg=carg+1;
il=varargin{carg};carg=carg+1;

switch RunOpt.pFcn % property type
case 'p_heat'
 if any(RunOpt.pTyp==1:2) % standard heat
  RunOpt.FieldDofs=20; RunOpt.PropertyDof=[];
 elseif RunOpt.pTyp==3 % Thermoelastic
  RunOpt.FieldDofs=[1:3 20]; RunOpt.PropertyDof=[];
 end
end
out=RunOpt;

%% #PropertyUnitType ---------------------------------------------------------
elseif comstr(Cam,'propertyunittype')

if ~isempty(strfind(Cam,'mat')) % m_heat
 error('Should now be a direct m_heat call')
%  if nargin==1;out=1; return; end % return subtypes ID
%  i1=varargin{carg}; out1={};
%  switch i1
%  case 1
%    st=...
%   {'MatId'    0  'sdtweb(''p_heat'')';
%   'Type'     0  '';
%   'k'        0  'conductivity';
%   'rho'      0  'mass density';
%   'C'        0  'heat capacity'
%   'alpha'    0 'heat exchange coefficient'};
%  otherwise; error('Not implemented subtype %i',i1);
%  end
else
 if nargin==1;out=1:3; return; end % return subtypes ID
 i1=varargin{carg}; out1={};
 switch i1
 case {1,2} % [ProId fe_mat('p_heat','SI',1) CordM Integ DIM] 
 st={ ...
         'ProId'   0   'sdtweb(''p_heat'')';
         'Type'    0   '';
         'CoordM'  0   '';
         'Integ'   0   '';
         'DIM'     0   '';
          };
 case 3
 st={ ...
         'ProId'   0   'sdtweb(''p_heat'')';
         'Type'    0   '';
         'ElasMatId'  0   '';
         'TherMatId'   0   '';
         'Integ'     0   '';
         'T0'     0   '';
          };
         
 otherwise; st={'ProId' 0 'sdtweb(''p_heat'')'; 'Type', 0, ''};
 end
end
if ~isempty(strfind(Cam,'cell')); out=st; else; out=[st{:,2}]; end

% -------------------------------------------------------------------------
elseif comstr(Cam,'subtypestring')

 i1=varargin{carg}; carg=carg+1;
 switch i1
 case 1;  out='3D volume';
 case 2;  out='3D surface element';
 case 3;  out='Thermoelastic coupling';
 otherwise; out='?';
 end

%% #Test ---------------------------------------------------------------------
elseif comstr(Cam,'test')

   model=femesh('testq8p divide 5 5');
   model.pl=m_heat('dbval 100 Aluminum');
   model.il=[113 fe_mat('p_heat','SI',1) 0 -3 2];
   
   [m,k,mdof]=fe_mknl(model);
%% #SetFace ------------------------------------------------------------------
elseif comstr(Cam,'setface')
    
 %model=p_heat('SetFace',model,Selelt,pl,il);
 model=varargin{carg}; carg=carg+1;
 SelElt=varargin{carg}; carg=carg+1;
 pl=varargin{carg}; carg=carg+1;
 il=varargin{carg}; carg=carg+1;
 
 [CAM,Cam,RunOpt.Load]=comstr('-load',[-25 2],CAM,Cam);
 
 % Add exchange elements:
 elt=feutil(sprintf('SelElt %s',SelElt),model);
 mpid=feutil('mpid',elt); i0=find(mpid(:,1)>0); 
 mpid(i0,1)=pl(1); mpid(i0,2)=il(1);
 elt=feutil('mpid',elt,mpid);
 model.Elt(end+1:end+size(elt,1),1:size(elt,2))=elt;
 
 % Deal with properties:
 if ~isfield(model,'il'); model.il=[]; end
 if ~isfield(model,'pl'); model.pl=[]; end
 % MATID - - -
 MatId=fe_mat('Getpl',model); MatId=MatId(:,1);
 i0=ismember(pl(1),MatId); % is this matid already exist ?
 if length(pl)>1 % pl is given
  if i0 % MatId already exist
   model=feutil(sprintf('RemoveMat%i',pl(1)),model);
   sdtw('_nb','MatId %i already exist : replaced.',pl(1))
  end
  model.pl(end+1,1:length(pl))=pl;
 else % only MatId given
  if ~i0 % MatId doesn't exist
   pl0=m_heat('default'); pl0=pl0.pl; pl0(1)=pl(1); pl=pl0;
   sdtw('_nb','MatId %i doesn''t exist : default is used.',pl(1))
   model.pl(end+1,1:length(pl))=pl;
  end
 end
 % PROID - - -
 ProId=fe_mat('Getil',model); if ~isempty(ProId);ProId=ProId(:,1);end
 i0=find(ProId==il(1)); % is this proid already exist ?
 if length(il)>1 % il is given
  if ~isempty(i0) % ProId already exist
   model=feutil(sprintf('RemovePro%i',il(1)),model);
   sdtw('_nb','ProId %i already exist : replaced.',il(1))
  end
  model.il(end+1,1:length(il))=il;
 else % only ProId given
  if ~isempty(i0)&&length(il)>1&&~strcmp(fe_mat('typep',il(i0,2)),'p_heat');
   sdtw('_nb','ProId %i already exist : replaced.',il(1))
   model=feutil(sprintf('RemovePro%i',il(1)),model);i0=[];
  end
  if isempty(i0) % ProId doesn't exist
   il0=p_heat('dbval d3 -3 2'); il0(1)=il(1); il=il0;
   sdtw('_nb','Adding ProId %i default.',il(1))
   model.il(end+1,1:length(il))=il;
  end
 end
 
 % Load - - -
 if ~isempty(RunOpt.Load)
  pl0=fe_mat('Getpl',model); alpha=pl0(pl0(:,1)==pl(1),6);
  if isempty(alpha); 
   error('MatId %i must be defined to define heat exchange load'); 
  end
  data2=struct('sel',sprintf('ProId%i',il(1)),...
      'def',RunOpt.Load*alpha,'DOF',.20);
  Case=fe_case(model,'getcase');
  st0=fe_def('StackLab',Case.Stack,'','Surface_load');
  model=fe_case(model,'Fvol',st0,data2);
 else
  sdtw('_nb',horzcat('You should build Load corresponding to heat ',...
    'exchange surface using:\n',...
    'model=fe_case(model,''Fvol'',''name'',',...
    'struct(''sel'',''ProId%i'',''def'',alpha*T,''DOF'',.20);'),il(1));  
 end
 
 out=model;
 
%% #End -------------------------------------------------------------------------
elseif comstr(Cam,'tablecall');out='';
elseif comstr(Cam,'cvs')
 out='$Revision: 1.55 $  $Date: 2025/04/07 17:08:29 $'; return;
else sdtw('''%s'' not known',CAM);
end
%% 
