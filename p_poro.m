function [out,out1,out2,out3]=p_poro(varargin)

%P_PORO material/property function poroelastic element support
%
%       Supported ELEMENT property subtypes are
%       1 : poroelastic volume 
%           [ProId fe_mat('p_poro','SI',1) Integ DIM] 
%         Integ : is rule number in integrules
%         DIM is problem dimension 2 or 3 D
%
%       Supported MATERIAL property subtypes see m_poro 
%
%       See <a href="matlab: sdtweb _taglist p_poro">TagList</a>
%       See also help fe_mat
%                doc  p_shell, fem (handling materials section), pl, fe_mat

%       Etienne Balmes
%       Copyright (c) 1990-2025 by SDTools, All Rights Reserved.
%	Use p_poro('cvs') for revision information

if nargin<1; help p_poro;return; end
if ischar(varargin{1}); [CAM,Cam]=comstr(varargin{1},1); il=[]; carg=2;
else; il=varargin{1};[CAM,Cam]=comstr(varargin{2},1); carg=3;
end

%% #propertyunittype ---------------------------------------------------------
if comstr(Cam,'propertyunittype')

 if nargin==1;out=1; return; end % return subtypes ID
 i1=varargin{carg};
 out1={};
 switch i1
 case 1 % [ProId fe_mat('p_poro','SI',1) Integ DIM] 
     st={...
         'ProId'             0   'help p_poro';
         'Type'              0   '';
         'Integ'             0   'Integration rule number';
         'DIM'               0   'Problem dimension'};

 otherwise; st={'ProId' 0 'help p_poro'; 'Type', 0, ''};
 end
 if ~isempty(strfind(Cam,'cell')); out=st; else; out=[st{:,2}]; end

%% #Default ------------------------------------------------------------------
elseif comstr(Cam,'default')

  out=p_poro('database');
  out=out(1);

%% #Dbval --------------------------------------------------------------------
elseif comstr(Cam,'dbval')

 while 1==1
  i1=strfind(comstr(Cam,-27),'-unit'); out1={};
  if ~isempty(i1)
   [Unit,i2,i3,i4]=sscanf(CAM(i1+5:end),'%s',1);
   i4=i1+[0:4+i4];CAM(i4)=''; [CAM,Cam]=comstr(CAM,1);
  else; Unit='';
  end
  i2=strfind(comstr(Cam,-27),'-punit');
  if ~isempty(i2)
   [PUnit,i3,i4,i5]=sscanf(CAM(i2+6:end),'%s',1);
   i5=i2+[0:5+i5];CAM(i5)=''; [CAM,Cam]=comstr(CAM,1);
  else; PUnit='';
  end
  [CAM,Cam]=comstr(CAM,6);
  [i1,CAM,Cam]=comstr(CAM,'','%i');
  if isempty(CAM)&&carg<=nargin; st=varargin{carg};carg=carg+1;
  else; st=CAM;end
  [mat,st1,i2]=p_piezo('database',st,varargin{carg:end});carg=carg+i2-3;
  if ~isempty(PUnit)
   r1=fe_mat(sprintf('convert %s %s',mat.unit,PUnit),mat.il(1:2));
   mat.il(2)=r1(2); mat.unit=PUnit;
  end
  if ~isempty(Unit)
   mat.il=fe_mat(sprintf('convert %s %s',mat.unit,Unit),mat.il);mat.unit=Unit;
  end
  r1=mat.il; if length(i1)==1; r1(1)=i1;end
  if ~isempty(il); i2=find(il(:,1)==r1(1)); else; i2=[];end
  if isempty(i2); i2=size(il,1)+1;end
  il(i2,1:length(r1))=r1;
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
  out.il=[MatId fe_mat('p_poro','SI',1) fe_mat('m_poro','SI',1) 1 2];
  % k rho C alpha 
  out.name='default poro';
  out.type='p_poro';
  out.unit='SI';

  out1='Poro properties';out2=carg;


  i1=strmatch(comstr(st,-27),comstr({out.name},-27));
  if ~isempty(i1); out=out(i1);
  elseif comstr(st,'shell') % piezo shell
   r1=comstr(st(6:end),-1,[1 1000 .01 0]);
   out=out(1);out.il(2+[1:length(r1)])=r1;
  end
%% #constitTopology
elseif comstr(Cam,'constittopology2d')
   dd=reshape([3:11]+2,3,3); dd(4,4)=14;dd(5,5)=14; out=int32(dd);
elseif comstr(Cam,'constittopology3d')
   dd=reshape([3:11]+2,3,3); dd(4,4)=14;dd(5,5)=14; out=int32(dd);

% #buildconstit --------------------------------------------------------------
% [constit,integ,Inits,Data]=p_poro('buildconstit',ID,pl,il);
elseif comstr(Cam,'buildconstit')

ID=varargin{carg};carg=carg+1;
pl=varargin{carg};carg=carg+1;
il=varargin{carg};carg=carg+1;
mat=pl(pl(:,1)==ID(1),:);
pro=il(il(:,1)==ID(2),:);

[st,i1,i2]=fe_mat('typep',pro(1,2));
if ~strcmp(st,'p_poro'); 
 error(sprintf('p_poro does not support rules %s',st));
end
ID=ID(:); 
ID(5:7)=[pro(1,3);pro(1,4);i2];  % Integrule, DIM, subtype
if pro(4)==2 % 2D
  if ID(4)==0;ID(4)=evalin('caller','RunOpt.Dim(2)');end
  [constit,integ,elmap,r3]=p_solid('buildconstit', ...
   [mat(3);max(il(:,1)+1);ID(4)*2;ID(4)], ...
   pl,[max(il(:,1))+1 fe_mat('p_solid','SI',2) 0 0 pro(3)]);%PlaneStrain N=0 Integ 
   out3=r3;out3.ConstitLab=[{'-1','m_type'} out3.ConstitLab(1:11) {'1'}];
   constit=[-1;mat(1,2);constit(1:11);1];% should have a non unity value here
   
   dd=p_poro('ConstitTopology2d');dd(find(dd))=constit(dd(find(dd)));
   out3.dd=dd;
   if sp_util('diag');fprintf('Constitutive strain matrix\n');disp(dd);end
  ID(3)=3*ID(4); % there are 3 DOF per node

elseif pro(4)==3
 error('3D not implemented yet');
 ID(3)=4*ID(4); % there are 4 DOF per node
else; error('Not a valid dimension');
end

out=constit(:);     % Return CONSTIT
out1=int32(ID(:));  % Return integ
out2=[]; % Elmap not used here

% -------------------------------------------------------------------------
% #Const EltConst=p_solid('ConstSolid',ElemF,integ,constit);
elseif comstr(Cam,'const')

 opt=varargin{carg};carg=carg+1;
 integ=varargin{carg};carg=carg+1;
 if carg<=nargin; constit=varargin{carg};carg=carg+1;
 else; constit=[];
 end
 out2=[]; % allow for MAPs, ...
 
 if ischar(opt) % Allow for integrule switching here
  if size(integ,1)>4; 
    opt=integrules(opt,double(integ(5,1)));
    if any(integ(5,:)~=integ(5,1));
     warning('Cannot deal with multiple integration strategies in the same group');
    end
  else; opt=integrules(opt);
  end
 end

 opt.Nw=size(opt.N,1);opt.Nnode=size(opt.N,2);
 if size(integ,1)<9;rule=[1 opt.Nw]; % By default start at 1 use all points
 else; rule=integ(8:9,1);rule=rule(:)';
 end
 if ~any(rule); 
   rule=[1 opt.Nw];
   sdtw('integ is assumed to give standard integration rule in integ(5:6)');
 end
 if integ(6)==2;     opt.DofLabels={'u','v','p'}; 
 elseif integ(6)==3; opt.DofLabels={'u','v','w','p'}; 
 else; error('Not a known dimension');
 end

 % 3D \int KgradT gradv - - - - - - - - - - - - - - - - - -
 if isequal(double(integ(6:7,1)),[3;1])
  error('Not a supported case');
 % 2D \int  - - - - - - - - - - - - - - - - - -
 elseif isequal(double(integ(6:7,1)),[2;1]) 
  % Strain energy
if sp_util('diag')
    opt.StrainString={{'u,x';'v,y';'u,y+v,x';'p,x';'p,y';'u';'v';'p'}, ...
        {'u';'v';'p'}};
    opt.ConstitTopology={'p_poro(''ConstitTopology2d'');'
        'eye(3)'};
    opt.FieldDofs=[1 2 19];
    a=integrules('StrainString',opt,rule);disp(a.StrainString{1});
end
  % define the deformation vector: row, NDN, DDL, NwStart NwRule
  opt.StrainDefinition{1}= ...
   [1 2 1 rule;2 3 2 rule; % e_xx = N,x u, e_yy=N,y v
    3 3 1 rule;3 2 2 rule; % shear 
    4 2 3 rule; 5 3 3 rule % N,x p, N,y p  % also contains C1
    ];
  opt.StrainLabels{1}=p_poro('convPoro2D');
  opt.ConstitTopology{1}=p_poro('ConstitTopology2d');

  % Kinetic energy
  opt.StrainDefinition{2}= [1 1 1 rule;2 1 2 rule;3 1 3 rule];
  opt.StrainLabels{2}={'u','v','p'};
  opt.ConstitTopology{2}=int32(diag([3 3 14])); % rho fluid=1

  opt=integrules('matrixrule',opt); % sum with Solid and Fluid Mass/Stiff

  % append rule to integrate : du . p,i C1
  % DDL1 DDL2 NDNi NDNj              constit ConstStep Nw W0
  r1=[integ(4)*[0 2] size(opt.w,1)*[0 1]  2 0 4 0; %du*p,x
      integ(4)*[1 2] size(opt.w,1)*[0 2]  2 0 4 0; %dv*p,y
     ];
 opt.MatrixIntegrationRule{1}(end+[1:2],:)=int32(r1); 
 
  % Append rule for C2
  r1=[integ(4)*[2 0] size(opt.w,1)*[1 0]  2 0 4 0; %dp,x*u
      integ(4)*[2 1] size(opt.w,1)*[2 0]  2 0 4 0; %dp,y*v
     ];
 opt.MatrixIntegrationRule{2}(end+[1:2],:)=int32(r1);

  % \int f v
  %opt.Rhs
  out1=2;
 else; error('Not a supported case');
 end
 opt.VectMap=reshape(1:opt.Nnode*length(opt.DofLabels), ...
   length(opt.DofLabels),opt.Nnode)';
 opt.material='multi';

 if ~isempty(constit)
  % Only keep terms needed for the current law.
  i1=find(any(constit,2))-1; 
  for j0=1:length(opt.MatrixIntegrationRule);
    r1=opt.MatrixIntegrationRule{j0}; if ~isempty(r1)
    opt.MatrixIntegrationRule{j0}=r1(find(ismember(double(r1(:,5)),i1)),:);
    end
  end
 end
 if nargout==0
  integrules('texstrain',opt);
  try; opt=integrules('stressrule',opt);integrules('texstress',opt);end
 else; out=opt;
 end
 
% #BuildDofOpt RunOpt=p_poro('BuildDofOpt',RunOpt,pl,il); --------------------
elseif comstr(Cam,'builddofopt')

RunOpt=varargin{carg};carg=carg+1;
pl=varargin{carg};carg=carg+1;
il=varargin{carg};carg=carg+1;
% xxx missing 2D/3D switch
RunOpt.FieldDofs=[01;02;19]; RunOpt.PropertyDof=[];
out=RunOpt;
% -------------------------------------------------------------------------
% convention labels
elseif comstr(Cam,'conv')

out=struct( ...
    'Poro2DL',{{'\epsilon_x','\epsilon_y','\gamma_{xy}', ...
    'p,x','p,y'}}, ...
    'Poro2D',{{'Exx','Eyy','Gxy','px','py'}});
 
[CAM,Cam]=comstr(CAM,5); if ~isempty(CAM); out=getfield(out,CAM);end


%% #SubTypeString ------------------------------------------------------------
elseif comstr(Cam,'subtypestring')

 i1=varargin{carg}; carg=carg+1;
 switch i1
 case 1;  out='Poro';
 otherwise; error('Unknown');
 end

%% #Formula ------------------------------------------------------------------
elseif comstr(Cam,'formula')

pl=varargin{2};
w=varargin{3}; w=w(:)*2*pi;s=1i*w;% Frequency Hz

[st,i1,typ]=fe_mat('type',pl(2));

switch typ % Type of poroelastic description

case 1
r1=m_poro('propertyunittype cell',typ);
r1(:,2)=num2cell(pl(1:size(r1,1))');r1(:,3)=[];r1=r1';
r1=struct(r1{:}); % structure with fields for constants

ga = 1./s;w=imag(s);
rho1 = r1.rhoS*(1 - r1.phip) ; % skeleton density
% viscosite dynamique %%%a identifier eta = 1.84e-5

% Calcul de K(w)
aa = (1 + 1i.*r1.rho0.*(w.*r1.Npr.*(r1.lambda2).^2)/(16.*r1.eta)).^0.5 ;
bb = 8.*r1.eta ./(1i.*(r1.lambda2).^2.*r1.Npr.*w.*r1.rho0) ;
cc = (1 + bb.*aa).^-1 ;
kk  = r1.qq.*r1.P0./(r1.qq - (r1.qq - 1).*cc) ;
% Transmissibilite
Q = (1 - r1.phip).*kk;
R = r1.phip.*kk;

b  = r1.sigmaR.*(r1.phip.^2)*...
    (1 + (4i.*r1.alpha.^2.*r1.eta.*r1.rho0.*w)./ ...
    (r1.sigmaR.*r1.lambda.^2.*r1.phip.^2)).^0.5 ;
rhoa  = r1.phip.*r1.rho0.*(r1.alpha - 1) ;
rho11 = rho1 + rhoa + ga.*b ;
rho12 = -rhoa - ga.*b ;
rho22 = r1.phip.*r1.rho0 + rhoa + ga.*b ;

% la masse solide equivalente
rhotild = rho11 -  rho12.^2./rho22 ;
% nomber d'onde fluide equivalent
Kf = w.*(rho22./R).^0.5 ;
% Coefficient de couplage
gamatild = r1.phip.*(rho12./rho22 -Q./R) ;
betatild = -w.^2.*rho22.*gamatild./r1.phip.^2;
kapa = r1.phip.^2./(w.^2.*rho22);

%Sortie
% SO = [rhoeq' gama' kapa' kk.^2' beta'];
rhoeq = rhotild(:);
gama = gamatild(:) ;
kapa = kapa(:) ;
K = Kf(:);
beta = betatild(:);
coef_afaib =R(:) ;

out=struct('data',[-w.^2 rhoeq gama kapa K.^2 beta coef_afaib r1.rho0*ones(size(w))], ...
    'lab',{{'-w.^2' 'rhoeq' 'gama' 'kapa' 'K' 'beta' 'coef_afaib' 'rho0'}});

otherwise; error('SubType %i not implemented',typ);
end

%% #Test -------------------------------------------------------------------------
elseif comstr(Cam,'test')

    
    
RunOpt.Layers=[0 2.1 32.1 52.1];
RunOpt.Divisions=[3 5 10];

% create mesh - - - - - - - - - - - - - - - - - - - - - - -
femesh('reset'); 
FEnode=[1 0 0 0 0 0 0; 2 0 0 0  100 0 0];
femesh('object beamline 1 2');
femesh('divide 10');
r1=[];
for j1=1:length(RunOpt.Layers)-1;
 r1=[r1 linspace(RunOpt.Layers(j1),RunOpt.Layers(j1+1),RunOpt.Divisions(j1))];
end
r1=unique(round(r1*1e6))/1e6;
femesh('extrude 0  0 1 0',r1);

% one property by layer
model=femesh('model0');
for j1=1:length(RunOpt.Layers)-1;
 i1=feutil('findelt innode {y>= & y <=}',model, ...
   RunOpt.Layers(j1),RunOpt.Layers(j1+1));
 RunOpt.NodeSets{j1}=feutil('findnode y>= & y <=',model, ...
   RunOpt.Layers(j1),RunOpt.Layers(j1+1));
 model.Elt(i1,5:6)=j1;
end
model.Elt=feutilb('SeparatebyProp',model);
model.Elt=feutil('orient',model);
model.Elt=feutil('setgroup 1:3 name q4p',model);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Isotropic solid (bottom plate) :  Thesis of O. Tanneau page 48 
model.pl=m_elastic('dbval 1',  ...
 [1 fe_mat('m_elastic','SI',1) 2.62e9 .45 1100 0 .0624]); % 
data=struct('sel','y==0','type','edge', ... 
             'eltsel','matid 1','dir',{{0 1 0}},'DOF',[.01;.02;.03]);
model=fe_case(model,'FSurf','Pressure',data);

% Elastic properties for poroelastic
model.pl=m_elastic(model.pl,'dbval 21', ...
  [21 fe_mat('m_elastic','SI',1)  292e3 .2 11.2 0 .06]);
i1=find(model.pl(:,1)==2); if isempty(i1); i1=size(model.pl,1)+1;end
pl2=m_poro('dbval 2 Ref');pl2(3)=21; % elast matid
model.pl(i1,1:length(pl2))=pl2;

% Properties for Air
model.pl=m_elastic(model.pl,'dbval 3 air');
% xxx : redo once m_poro has unit conversion
% model.unit='MM';model.pl=fe_mat(['convertSI' model.unit],model.pl);
model.unit='SI'; model.Node(:,5:7)=model.Node(:,5:7)/1000;

model.il=p_solid('dbval 1 d2 2','dbval 2 d2 2','dbval 3 d2 2');
i1=find(model.il(:,1)==2); if isempty(i1); i1=size(model.pl,1)+1;end
model.il(i1,1:4)=[2 fe_mat('p_poro','SI',1) 2 2];

[constit,integ,r1,r3]=p_poro('buildconstit',[2;2;12;4],model.pl,model.il);
disp(r3.dd);disp(r3.ConstitLab);

% Boundary conditions
model=fe_case(model,'FixDof','SlidingEdge','x==0|x=100-dof 1');

if nargout>0; out=model;
else;
 cf=feplot;cf.model=model;fecom('colordatamat');
 % - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 [Case,model.DOF]=fe_mknl('init',model);
 i1=[fe_c(model.DOF,(1:3)'/100,'ind')
    fe_c(model.DOF,(19)'/100,'ind')];
 k=fe_mknl('assemble NoT',model,Case,1); % assemble stiffness
 m=fe_mknl('assemble NoT',model,Case,2); % assemble mass
 if sp_util('issdt');figure(1);ii_plp('spy',k(i1,i1))
 else;figure(1);spy(k(i1,i1))
 end
end

%% #split  ------------------------------------------------------
% SE = p_poro('split',SE,struct('Tag','p2'));
elseif comstr(Cam,'split');

SE=varargin{carg};carg=carg+1;
RO=varargin{carg};carg=carg+1;

SE2=struct;
RO.imat=find(~cellfun('isempty',regexp(SE.Klab,RO.Tag))); % Poro MK tag
lab=SE.Klab(RO.imat);
RO.lab= {[lab{1} 's'],[lab{2} 's'], ...
    [lab{1} 'f'],[lab{2} 'f'],[lab{2} 'c1'],[lab{2} 'c2']};
DOF=SE.DOF(find(full(diag(SE.K{RO.imat(1)})|diag(SE.K{RO.imat(2)}))));% PoroDof

i_s=fe_c(SE.DOF,fe_c(DOF,(1:6)'/100,'dof'),'ind');
Ts=fe_c(SE.DOF,SE.DOF(i_s));
i_f=fe_c(SE.DOF,fe_c(DOF,.19,'dof'),'ind');
Tf=fe_c(SE.DOF,SE.DOF(i_f));
% split into Ms Ks Mf Kf C1 C2
Ms=feutilb('tkt',Ts,SE.K{RO.imat(1)}(i_s,i_s));
Mf=feutilb('tkt',Tf,SE.K{RO.imat(1)}(i_f,i_f));
Ks=feutilb('tkt',Ts,SE.K{RO.imat(2)}(i_s,i_s));
Kf=feutilb('tkt',Tf,SE.K{RO.imat(2)}(i_f,i_f));
C1=SE.K{RO.imat(2)}(i_s,i_f); C1=Ts'*C1*Tf+Tf'*C1'*Ts;
C2=SE.K{RO.imat(1)}(i_s,i_f); C2=Ts'*C2*Tf+Tf'*C2'*Ts;

SE.K(RO.imat)=[];SE.K=[SE.K(:)' {Ms Ks Mf Kf C1 C2}]; 
SE.Klab(RO.imat)=[];SE.Klab=[SE.Klab(:)' RO.lab]; 
out=SE;

%% #end ------------------------------------------------------
elseif comstr(Cam,'tablecall');out='';
elseif comstr(Cam,'info'); matgui('info','p_poro');
elseif comstr(Cam,'cvs')
  out='$Revision: 1.20 $  $Date: 2025/04/07 17:08:33 $';
else; error('''%s'' not known',CAM);
end


