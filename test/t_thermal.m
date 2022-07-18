function out=t_thermal(varargin); 

% Tests/developement of thermal capabilities in the mat_og element family
% Field interpolation
% Prestress setting

if nargin==0
   
  cd(fileparts(fileparts(which('fe_mk'))));eval('sdtcheck path');cd ..
  t_thermal('base')  % Example with uniform temperature setting
  t_thermal('map')
  t_thermal('SetGState') % Test of analytic gstate
  t_thermal('TexRule')
  % SDT tests 
  t_thermal('FieldInterp')
  t_thermal('ConstitInterp') % Interpolate constit value in InfoAtNode
  t_thermal('rivlin') 
  t_thermal('map') % rerun as SDT to get final part
  return
  
end

%#ok<*NOSEM,*NASGU,*ASGLU>
[CAM,Cam]=comstr(varargin{1},1);carg=2;

%% #Base --------------------------------------------------------------------
if comstr(Cam,'base') 

feplot;sdtdef('diag',0);
model=femesh('teststruct hexa8b');

% Define thermal expansion coefficients and reference temperatures
[i1,i2]=intersect(model.pl(:,1),[100;112]);
model.pl(i2,8:9)=[1e-5 20;2e-5 20];


model.Elt(2:end,9)=100;
model.Elt(feutil('findelt withnode {z==0}',model),9)=112;
[model.Elt,elt]=feutil('removeElt matid 112',model);
model=feutil('addelt',model,elt);
model.Elt=feutil('setgroup2 proid 112',model);
model.il=[];model=p_solid('default',stack_rm(model,'pro'));
if exist('feutilb','file')
% model.Elt=feutilb('SeparatebyMat',model.Elt);
end

% Initialize case and compute temperature generated stresses

% Uniform temperature of 30 degrees (nominal 20)
defT=struct('def',ones(size(model.Node,1),1)*30,'DOF',model.Node(:,1)+.20);
model=fe_case(model,'DofSet','ThermalState',defT);
[Case,model.DOF]=fe_mknl('init',model);
C2=fe_stress('thermal',model,defT);C0=C2;
Case.GroupInfo(:,5)=C2.GroupInfo(:,5);
b=fe_mknl('assemble',model,Case,103); % reference thermal pre-stress
Case.GroupInfo(:,5)={[]};Case.GroupInfo{end}.material='Elastic3DNL';
dc=struct('DOF',model.DOF,'def',zeros(length(model.DOF),2));
sdtdef('diag',12)
'xxxNot reimplemented sdtdef 12'; k=fe_mknl('assemble',model,Case,dc,1);
dc.def(:,2)=0;
sdtdef('diag',0)
[k1,C2,dc]=fe_mknl('assemble',model,Case,dc,5);

k=k1;
RO.kend=full(k1(end)); % last value pre_stressed
RO.k0end=full(k(end)); % last value pre_stressed

if norm(b+dc.def(:,2))/norm(b)>1e-10; error('Mismatch'); end
if mean(diag(k))/mean(diag(k1))<1; 
    error('Thermal stress with no displacement lowers stiffness')
end

sdtdef('diag',0);
EC=rmfield(Case.GroupInfo{end},'RhsDefinition');
C2.GroupInfo{end}=EC;
de=dc;de.def(:,2)=0;[k2,C2,de]=fe_mknl('assemble',model,C2,de,5);
%sparse(de.def(:,2))


%cf=feplot;cf.model=model;cf.def={b,model.DOF};fecom colordataa

% Compute deformation under thermal load
% The dc.def(:,2) computes the internal load, balanced with external (so
% negative sign if used as RHS in thermal)

q=ofact(k1,-Case.T'*dc.def(:,2)); 
def=struct('def',Case.T*q,'DOF',model.DOF);

% For a clean fe_simul example 
if sp_util('issdt')
   sdtweb('_link','sdtweb ofdemos(''ThermalCube'')','Clean example')
end
if 1==2 % Non linear static computation
  opt=fe_simul('nlstaticopt');opt.RelTol=-1e-3;
  d2=fe_time(opt,model);
end
% formal setting would be useful too
% mo0=feutil('setpro 1',mo0,'gstate', ...
%       struct('dir',{{RunOpt.T/RunOpt.il600(6)*cos(r1/180*pi),...
%                      RunOpt.T/RunOpt.il600(6)*sin(r1/180*pi)}},...
%                     'lab',{{'Exx','Eyy'}}));

%% #Rivlin -------------------------------------------------------------------
elseif comstr(Cam,'rivlin') 
% RivlinCube Expansion - - - - - - - - - - - - - - - - - - -
[model,dc,Case]=ofdemos('ThermalCube');
dc1=struct('DOF',model.DOF,'def',zeros(length(model.DOF),2));
% model=stack_rm(model,'pro');model.il=[];model=p_solid('default;',model);
% k1 has thermal expand in Case.GroupInfo{1,7}
[k1,C2,dc]=fe_mknl('assemble NoT',model,Case,dc1,5); 
[S0,C0]=fe_case(model,'assemble -matdes 2 1 5 -NoT -SE');
%normest(k1-S0.K{2})
%normest(k1-S0.K{3})

% Propagate temperature varying properties with table handling
sdtdef('diag',12)
% table for rho eta E(T) nu(T) G(T) alpha(T)
C2=Case;C2.GroupInfo{4}=model.pl';% E nu alpha
C2.material='Elastic3DNL';
dc1=struct('DOF',model.DOF,'def',zeros(length(model.DOF),2));
[k2,C2,dc1]=fe_mknl('assemble NoT',model,C2,dc1,5); 
if normest(k2-k1)/normest(k2)>1e-10; error('Mismatch');end

%% Isotropic elastic properties with interpolation
% sdtweb m_elastic('IsotropicInterp')
% [model,dc,Case]=ofdemos('ThermalCube');
mat=struct('pl',S0.pl,'type','m_elastic','unit','SI', ...
    'Rho',struct('X',[20;40],'Xlab',{{'T'}},'Y',S0.pl([5;5])'.*[.9;1.1]), ...
    'E',struct('X',[20;40],'Xlab',{{'T'}},'Y',S0.pl([3;3])'.*[.9;1.1]), ...
    'Nu',struct('X',[20;40],'Xlab',{{'T'}},'Y',S0.pl([4 4])'.*[.9;1.1]));
mat.pl(8:9)=0;
model=stack_set(model,'mat','Steel',mat);model=stack_rm(model,'pro');
model.pl=[];
model=fe_case(model,'DofSet','ThermalState', ...
    struct('DOF',model.Node(:,1)+.2,'def',ones(size(model.Node,1),1)*30));
[SE,CE]=fe_case(model,'assemble -matdes 2 1 5 -NoT -SE');
if normest(SE.K{3}-S0.K{2})/normest(S0.K{2})>1e-10; error('Mismatch MatDes 5');end
if normest(SE.K{1}-S0.K{1})/normest(S0.K{1})>1e-10; error('Mismatch 2');end
if ~isfield(CE.GroupInfo{end},'CTable'); error('Missing interp');end

% Now redo with 
mat.pl=S0.pl; model=stack_set(model,'mat','Steel',mat);
% now check using the of_mk element wise callback
sdtdef('diag',8); % the init is done for element callback
[C3,model.DOF]=fe_mknl('init',model);
EC=C3.GroupInfo{1,8};
sdtdef('diag',0);[C4,model.DOF]=fe_mknl('init',model);

if 1==2 % strategy with callback should be reactivated
 dc1=struct('DOF',model.DOF,'def',zeros(length(model.DOF),2));
 [k2,C3,dc1]=fe_mknl('assemble NoT',model,C3,dc1,5);
 if normest(k2-k1)/normest(k2)>1e-10; error('Mismatch');end
end
[k2,C4,dc1]=fe_mknl('assemble NoT',model,C4,dc1,5);
if normest(k2-k1)/normest(k2)>1e-10; full(diag(k2)./diag(k1))
    error('Mismatch 2');
end

% should use Elastic3DNL with matrix type 5
EC.material='Elastic3DNL';EC.pEC=zeros(max(of_mk('pecgetcall')),2,'int32');
C3.GroupInfo{end}=EC;C3.GroupInfo{1,5}=[]; C3.GroupInfo{7}.lab={'T'};
dc1=struct('DOF',model.DOF,'def',zeros(length(model.DOF),2));
sdtdef('diag',0);[k2,C3,dc1]=fe_mknl('assemble NoT',model,C3,dc1,5);
C3.GroupInfo{end}.pEC
if normest(k2-k1)/normest(k2)>1e-10; error('Mismatch 3');end

if 1==2
   d_ubeam
   
   mo2=cf.mdl.GetData;
   mo2.Elt=feutil('selelt selface',mo2);
   mo2.Elt=feutil('set group 1 name q4p mat 1 pro 1',mo2);
   mo2.pl=m_elastic('dbval 1 Strain');
   mo2.il=p_solid('dbval 1 d3 -1');
   [C2,mo2.DOF]=fe_mknl('init',mo2);
   integrules('texstrain',C2.GroupInfo{end});
   C2.GroupInfo{2}(4,:)=int32(21);
   sdtdef('diag',12)
   C2=fe_stress('stress gstate',mo2,def);
   
end

%% #SetGState Test analytic set of gstate ------------------------------------
elseif comstr(Cam,'setgstate') 

 [mo1,dc,C1]=ofdemos('ThermalCube');
 d2=dc;d2.def=d2.def(:,1);%d2=fe_def('subdef',dc,1);
 C2=fe_stress('stress gstate',mo1,d2);data=C2.GroupInfo{5};
 data.X{1}
 
 % sdtweb('vectfromdir')
 data=struct('dir',{{'1','3'}},'lab',{{'Exx','Ezz'}});
 [EGroup,nGroup]=getegroup(mo1.Elt);jGroup=1;
 cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
 eltid=feutil('eltid',mo1);
 C3=elem0('VectFromDirGState',mo1,data,C2,C2.Node,1);
  
elseif comstr(Cam,'texrule') 
% Display the integration rule formula - - - - - - - - - - - - - 
 
 EC=integrules('hexa8',-1); % 3D solid
 EC.StressLabels={p_solid('convStress3DL'),{'u','v','w'}}
 p_solid('const',EC,[1 1 24 8]',[1 1],[50 1 1 1])
 
 
% #Map Test orientation MAPs -------------------------------------------------
elseif comstr(Cam,'map') 

model=femesh('testhexa8b'); [i1,model.Elt]=feutil('eltidfix',model);

nodeE=model.Node(:,5:7);opt=integrules('hexa8',-1);
integrules('buildndn',3,opt,nodeE);N1=opt.NDN+0;
bas=basis('rotate',[],'rz=80;',1);cGL=reshape(bas(7:15),3,3);

nodeE(:,5:7)=repmat(cGL(:,1)',8,1); % interp xe
nodeE(:,8:10)=repmat(cGL(:,2)',8,1); % interp ye
integrules('buildndn',32,opt,nodeE);N2=opt.NDN+0;
if norm(N1(:,2:4)*cGL-N2(:,2:4))>1e-4; error('Mismatch');end
opt.nodeE=nodeE;of_mk('buildndn',int32(32),opt);N2=opt.NDN+0;
if norm(N1(:,2:4)*cGL-N2(:,2:4))>1e-4; error('Mismatch');end

model.pl=[100 fe_mat('m_elastic','SI',6) [1 2 3 0 0 0 1 1 1]*1e9];
model.il=p_solid('dbval 111 d3 -1');
data=struct('EltId',i1(i1~=0),'bas',basis('rotate',[],'rz=0;',1));
model=stack_set(model,'info','EltOrient',data);
[C1,model.DOF]=fe_mknl('init',model);
sdtdef('diag',12);k=fe_mknl('assembleNoT',model,C1,1);
RO.K={of_mk('xkx_trans',reshape(data.bas(7:15),3,3),full(k))};
full(diag(RO.K{1}(1:3,1:3)))*16e-9
data=struct('EltId',i1(i1~=0),'bas',basis('rotate',[],'rz=90;',1));
cGL=reshape(data.bas(7:15),3,3);
model=stack_set(model,'info','EltOrient',data);
[C1,model.DOF]=fe_mknl('init',model);
C1.GroupInfo{4}(4:end)=0;
C1.GroupInfo{2}(4,:)=32;sdtdef('diag',12);k=fe_mknl('assembleNoT',model,C1,1);
full(diag(k(1:3,1:3)))*16e-9

%% #Fully_anisotropic_material - - - - - - - - - - - - - - - - - - - - -2
model=femesh('testhexa8b'); [i1,model.Elt]=feutil('eltidfix',model);
data=struct('EltId',i1(i1~=0),'bas',basis('rotate',[],'rz=10;rx=20;ry=30',1));
%data=struct('EltId',i1(i1~=0),'bas',basis('rotate',[],'rz=0;rx=0;ry=0',1));
cGL=reshape(data.bas(7:15),3,3);
model=stack_set(model,'info','EltOrient',data);
DD=rand(6);[DD,r1]=svd(DD);DD=DD*r1*DD';
model.il=p_solid('dbval 111 d3 -1'); % only one integ point for check
model.pl=[100 fe_mat('m_elastic','SI',3) ...
    DD([1 7:8 13:15 19:22 25:29 31:36])];
% use the orientation map
[C1,model.DOF]=fe_mknl('init',model);
sdtdef('diag',12);k=fe_mknl('assembleNoT',model,C1,5);
sdtdef('diag',0);k2=fe_mknl('assembleNoT',model,C1,5);
if norm(full(k-k2))>sqrt(eps)*norm(full(k)); error('Mismatch k1 k2');end
% move the mesh
mo1=model;mo1.Node(:,5:7)=mo1.Node(:,5:7)*cGL;
mo1.Stack={};[C1,model.DOF]=fe_mknl('init',mo1);
sdtdef('diag',0);k1=fe_mknl('assembleNoT',mo1,C1,1);
k1=of_mk('xkx_trans',cGL',full(k1)); % back transform to global coord
%[full(k(1:6,1:6));zeros(1,6);full(k1(1:6,1:6))]
if norm(full(k-k1))/normest(k)>sqrt(eps); error('Mismatch k k1');end

model=femesh('testhexa8b'); [i1,model.Elt]=feutil('eltidfix',model);
data=struct('EltId',i1(i1~=0),'bas',basis('rotate',[],'rz=90;rx=90',1));
%data=struct('EltId',i1(i1~=0),'bas',basis('rotate',[],'rz=0;rx=0;ry=0',1));
cGL=reshape(data.bas(7:15),3,3);
model=stack_set(model,'info','EltOrient',data);
DD=diag(1:6)';
model.il=p_solid('dbval 111 d3 -1'); % only one integ point for check
model.pl=[100 fe_mat('m_elastic','SI',3) ...
    DD([1 7:8 13:15 19:22 25:29 31:36])];
% use the orientation map
[C1,model.DOF]=fe_mknl('init',model);k2=fe_mknl('assembleNoT',model,C1,5);


%% #EltOrient+D11(T) Now add a thermal field and an interpolated value -2
% (requires DW2DE2 computation from PL in C)
% sdtweb m_elastic IsotropicInterp % Constit(T)
% sdtweb p_solid('constitutive law interpolation') % InfoAtNode

if ~sp_util('issdt'); return; end
sdtw('_ewt','xxxEB: need to revise !')
% this entry generates CTable entry in GroupInfo, but value itself is not exploited
defT=struct('def',ones(size(model.Node,1),1)*5,'DOF',model.Node(:,1)+.20);
model=fe_case(model,'DofSet','ThermalState',defT);

%model.pl=m_elastic('dbval 100 steel');
%feval(p_solid('@ConstitLab'),'m_elastic.1')
mat=struct('pl',model.pl,'type','m_elastic','unit','SI', ...
    'D11',struct('X',[-10;20],'Xlab',{{'T'}},'Y',model.pl([3;3])'.*[.5;1.5]));
model=stack_set(model,'mat','orthotropic',mat);
[C1,model.DOF]=fe_mknl('init',model); % m_elastic('AtNodeGstate')
% C1.GroupInfo{8}
if ~isfield(C1.GroupInfo{8},'CTable')
  fprintf('ErrorTODO : interpolate G11 for orthotropic material');
elseif abs(C1.GroupInfo{1,4}(3)-DD(1)/2)>1e-10
    error('Mismatch interp DD(-10)');
    %feval(elem0('@ConstitInterp'),C1.GroupInfo{8}.Ctable)
end

% Back at 5degrees =need to declare T or RefT in stack otherwise min value is taken
%sdtweb fe_mat('field_interp')
model=stack_set(model,'info','T',5);
k3=fe_mknl('assembleNoT',model,C1,5);
if norm(full(k3-k2))/normest(k2)>sqrt(eps); error('Mismatch k k2');end
% sd pathc; sdtweb of_mk_pre.c#nonlin_elas

if 1==2 % test orientation propagation for contact
t_contact('contactcubeholemaster-int-1');
cf=feplot;
MAP=feutilb(sprintf('shellmapnodepos -xe %.15g %.15g 0',[1 
0]),cf.mdl,'group3');
cf.Stack{'MAP','Group_3'}=MAP;

q0=fe_time(cf.Stack{'TimeOptStat'},cf.mdl.GetData);
end 

%% #ConstitInterp constitutive interpolation -----------------
% This uses SDT curve interpolation - - - - - - - - - - - -
elseif comstr(Cam,'constitinterp') 

if ~sp_util('issdt')
 clear functions
 cd(fileparts(fileparts(which('fe_mk'))));eval('sdtcheck path');cd ..
end
model=d_mesh('RveConstitInterp');
SE=fe_mknl(stack_rm(model,'pro'),'NoT');

[C1,model.DOF]=fe_mknl('init',model);
m=fe_mknl('assemble',model,C1,2);
if norm(m-SE.K{1},'inf')>1e-10; error('Mismatch on m interp');end
n1=sortrows(feutil('getnode x==0 & y==0',model),7);
i1=fe_c(SE.DOF,n1(:,1)+.03,'ind');
model.il(model.il(:,1)==111,4)=-1;[C1,model.DOF]=fe_mknl('init',model);
for j1=1:(20*sdtdef('isinteractive')+10) % Repeat for omp verification
 k=fe_mknl('assemble',model,C1,1);
 r1=full([diag(SE.K{2}(i1,i1)) diag(k(i1,i1))]);% SE uniform, k variable
 if j1==1;r0=r1;
 elseif norm(r1-r0)>1e-8;
     [r1./r0]
     error('OMP change %i',j1);
  end
 r2=r1(:,2)./r1(:,1)./n1(:,7);
 if r2(2)>1;disp(r2);error('Mismatch on k variation,step %i',j1);end
end
if sdtdef('isinteractive'); figure(1);plot(r2); end
%def=fe_eig(model);C1=fe_stress('stress -gstate',model,def);

% now verify stress value
data=struct('sel','groupall','dir',{{0,0,'z.^1'}},'DOF',[.01;.02;.03]);
r2=elem0('VectFromDirAtDof',model,data);
C1=fe_stress('stress -gstate',model,r2);C2=C1.GroupInfo{5};
C2.Y=squeeze(C2.Y)/210e9;
if norm(C2.Y(3,:)-(.05:.1:.95))>1e-10;error('Mismatch');end
%feplot(model,r2);fecom('scc1')


%% Now interpolate on a plate (not yet working)
model=femesh('testquad4 divide 10 10 1');
model=feutil('setmat 100 eta .1',model);% feutilb('_write',model)

data=struct('dir',{{'7800','1-x'}}, ...
    'lab',{{'rho','eta'}});
[dd,constit]=feutil('getdd',[100 110],model); 
ind=find(~cellfun(@isempty,regexp(dd.ConstitLab,'D\d')));
r2=cellfun(@(x)sprintf('%.15g*x',x),num2cell(dd.dd),'uni',0);
data.dir=[data.dir r2(:)'];
data.lab=[data.lab dd.ConstitLab(ind)];
pro=struct('il',110,'type','p_solid','MAP',data);
model=stack_set(model,'pro','WithMap',pro);
[SE,CE]=fe_case(model,'Assemble -matdes 2 1 4 -SE -NoT');
%cf=feplot(model);cf.def=struct('def',diag(SE.K{3}),'DOF',SE.DOF)

% sdtweb of_mk_pre.c#constitInterp

%% #FieldInterp Test generic field interpolation -----------------
% This uses SDT curve interpolation - - - - - - - - - - - -
elseif comstr(Cam,'fieldinterp') 

if ~sp_util('issdt')
 clear functions
 cd(fileparts(fileparts(which('fe_mk'))));eval('sdtcheck path');cd ..
end

model=femesh('testquad4');
mat=m_elastic('database steel');mat.pl(1)=100;
st1=feval(mat.type,'propertyunittype cell',1);
i1=model.pl(:,1)==100;pl=model.pl(i1,:);model.pl(i1,3:end)=0;
mat=struct('pl',pl, ...
    'type','m_elastic','unit','SI', ...
    'E',struct('X',{[-10;20]},'Xlab',{{'T'}},'Y',pl([3;3])'), ...
    'Nu',struct('X',[20;40],'Y',(pl([4 4]).*[.9 1.1])'));
model=stack_set(model,'mat','Steel1',mat,'info','T',30);
pro = p_shell('database 100 MindLin .1');
pro.h=struct('X',{[1;2]},'Xlab',{{'T'}},'Y',pl([3;3])');
model.pl(model.pl(:,1)==100,:)=[];
r1=fe_mat('getpl',model); % sdtweb fe_mat field_interp
if norm(r1(r1(:,1)==100,:)-pl)>1e-10; error('Mismatch');end

mat.Nu.Xlab={'p'};
%mat.E.X={mat.E.X(:)'};%mat.E.Y=mat.E.Y;
model=stack_set(model,'mat','Steel1',mat,'info','p',30);
r1=fe_mat('getpl',model);if norm(r1(r1(:,1)==100,:)-pl)>1e-10; error('Mismatch');end

% Interp 
model=stack_set(model,'pro','thick',pro,'info','T',1.5);
r1=fe_mat('getil',model);

% test pre-emptive stack entry
mdl=femesh('testhexa8');
feutil('info',mdl)
pl0=mdl.pl(1,:);
pl0(1,3)=7777777;
mdl=stack_set(mdl,'mat','matid100',struct('pl',pl0,'type','p_solid'));
pl2=mdl.pl; mdl.pl=[]; pl=fe_mat('getpl',mdl); % a Warning is expected here
if pl(pl(:,1)==100,3)==pl2(pl2(:,1)==100,3)
    error('stack is pre-emptive');
end


%% #End
else error('%s unknown',CAM);
end



