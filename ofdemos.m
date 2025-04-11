function [out,out1,out2]=ofdemos(varargin);

% OpenFEM meshing scripts used by various demos
% This is used to provide shorted demos for inclusion in the manual
%
% ofdemos AnnularHeat   : 2D heat equation example with analytic comparison
% ofdemos HeatStatic    : 3D static heat problem
% ofdemos HeatTransient : 3D transient heat problem
% ofdemos ThermalCube   : 3D thermally pre-stressed cube
%
% See <a href="matlab: sdtweb _taglist ofdemos">TagList</a>

%	E. Balmes, A. Nassiopoulos, ...
%       Copyright (c) 2001-2025 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       Use ofdemos('cvs') for revision information

if nargin==0; sdtweb('_taglist','ofdemos');return; end
[CAM,Cam]=comstr(varargin{1},1);carg=2;

%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*NOSEM>

%% #Heat Equation demos --------------------------------------------------------------------------
%% #HeatRed : reduced heat transient based on snapshot (REQUIRES SDT)
if comstr(Cam,'heatred')
    
model=ofdemos('heatstatic'); %See sdtweb ofdemos('heatstatic')
model=ofdemos('heatLoadPen',model); % Loading case with penalized heat source

% Transient using full model - - - - - - - - - - - - - - - - - - - - - - - - - 
opt=struct('Method','euler','Opt',[0.5 0 0 5e-3 100]);% [theta 0 t0 deltaT Nstep]
def=fe_time(opt,model);

cf=feplot(model,def); 
fecom ch2;
fecom('colordataA-colorbartitle"T"'); fecom('ColorScaleOneInstant');

%% Generate reduction basis through snapshot, here actually 4 vectors

[SE,CE,Load]=fe_case(model,'assemble -matdes 2 3 1 -SE -NoT -load');

dred=fe_def('subdef',def,1:10:size(def.def,2));
[T,fr]=fe_norm(dred.def,SE.K{1},SE.K{3});
TR=struct('def',T,'DOF',SE.DOF);

% Prepare for new solution based on reduced model (here 4 DOFs)
SE.K=feutilb('tkt',T,SE.K); SE.DOF=(1:size(T,2))'+.99;
Load.def=T'*Load.def; Load.DOF=SE.DOF;
Load=feutil('rmfield',Load,'curve','ID','bset');
CE.Stack={'DofLoad','RedLoad',Load}; CE.DOF=SE.DOF; SE.Case=CE; 
mo2=stack_set(model,'SE','MVR',SE);
d2=fe_time(opt,mo2); % Run using 4 DOFs

% expand and display
d3=d2;d3.def=TR.def*d2.def;d3.DOF=TR.DOF; cf.def=d3;

d3.def(:,1)=d3.def(:,2);def.def(:,1)=def.def(:,2); % first step not at 0

%% display in iiplot
iiplot('curveinit',{'curve','ReducedTransient',d3;
    'curve','FullTransient',def});
iicom('ch',{'DOF',3.20})
    
    
%% #HeatStatic Static heat equation demo by A. Nassiopoulos
elseif comstr(Cam,'heatstatic')
    
femesh('reset')
FEnode=[1 0 0 0  0 0 0; 2 0 0 0  0 5 0;
        3 0 0 0  0 5 4; 4 0 0 0  0 0 4];
FEel0=[Inf abs('quad4');1 2 3 4 1 1];
femesh(';divide 4 5;addsel;');
femesh(';selgroup1;extrude 10 1 0 0;orientel0'); % group 1: volume elements
FEelt=FEel0;
% Gamma1 : % group 2: surface elements (x=0)
femesh(';selelt selface&withoutnode{x~=0};addsel'); 
% Gamma2 : % group 3: surface elements (x=10)
femesh(';selelt selface&withoutnode{x~=10};addsel'); 
% Gamma 3:6, group 4, all other faces
femesh(';selelt selface&withnode{x~=0}&withnode{x~=10};addsel');

femesh('set group1 mat 1 pro 1')
femesh('set group2 mat 2 pro 2 name q4p')
femesh('set group3 mat 3 pro 3 name q4p')
femesh('set group4 mat 4 pro 4 name q4p')
model=femesh('model');
                                   
model.pl=[1 fe_mat('m_heat','SI',1) 400 1 1 1   % k rho C alpha  
          2 fe_mat('m_heat','SI',1) 400 1 1 1   % k rho C alpha  
          3 fe_mat('m_heat','SI',1) 400 1 1 1e10  % k rho C (alpha_0: not used)
          4 fe_mat('m_heat','SI',1) 400 1 1 1   % k rho C alpha  
          ];
model.il=[1 fe_mat('p_heat','SI',1) 0 -3 3;   % volume integration
          2 fe_mat('p_heat','SI',2) 0 -3 3;   % exchange surface integration
          3 fe_mat('p_heat','SI',2) 0 -3 3;   % exchange surface integration
          4 fe_mat('p_heat','SI',2) 0 -3 3];  % exchange surface integration
      
% sort nodes for cleaner outputs
[i1,i2]=sortrows(model.Node(:,5:7));model.Node=model.Node(i2,:);
out=model;

%% --------------------------------------------------------------------------
%% #HeatTransient Transient heat equation demo by A. Nassiopoulos
elseif comstr(Cam,'heattransient')


% Material constants
alpha = 10;
diffu = 240;
rho = 2700;
heat = 900;

% Number of computed modes
Nm = 50;

disp('Mesh generation')
femesh('reset')
femesh

FEnode=[]; FEelt=[]; FEel0=[];
%
h=1.5;  % hauteur tablier 
L=3; L2=2; % 1/2 largeur

e=.5; % epaisseurs de beton
%
node = [0   0.    -h    ; 
        0   L2-1.2*e  -h    ; 
        0   L2    -h-e  ; 
        0   0.    -h-e  ];
model=feutil('Objectquad 1 1',node,5,4); % matid=1, proid=1
node = [0   L2    -h-e ;
        0   L     0.     ;
        0   L-2.2*e   -e ;
        0   L2-1.2*e  -h     ];
model=feutil('Objectquad 1 1',model,node,5,4); % matid=1, proid=1
node = [0   L         0. ;
        0   L-2.2*e   -e ;
        0   0.         -e ;     
        0   0          0. ];
model=feutil('Objectquad 1 1',model,node,4,5); % matid=1, proid=1

FEnode=model.Node; FEel0=model.Elt;
femesh('addsel');
femesh(';symsel 0 0 1 0;addsel;join group1:6');
FEel0=FEelt;

femesh(';selgroup1;extrude 20 0.2 0 0;orientel0'); % group 1: volume elements

FEelt=FEel0;
femesh(';selelt selface&withoutnode{x~=0};addsel'); % group 2 : surface x=0
femesh(';selelt selface&withoutnode{x~=4};addsel'); % group 3 : surface x=4
FEel0=FEelt;
femesh(';selelt selface&withoutnode{z~=0};addsel'); % surface z=0 (temporarily group 4)
femesh(';selelt selface&withoutnode{z~=-2};addsel'); % surface z=-h (temporarily group 5)
femesh(';selelt selface&withoutnode{group2|group3|group4|group5}&withnode{y>=2 & y<=3};addsel'); % surface y=L (temporarily group 6)
femesh(';selelt selface&withoutnode{group2|group3|group4|group5}&withnode{y<=-2 & y>=-3};addsel'); % surface y=L (temporarily group 7)
femesh('join group5:7') % group 5 : outer surfaces
FEel0=FEelt;
femesh(';selelt selface&withoutnode{z~=-0.5};addsel'); % surface z=-e (temporarily group 6)
femesh(';selelt selface&withoutnode{z~=-1.5};addsel'); % surface z=-h+e (temporarily group 7)
femesh(';selelt selface&withoutnode{group2|group3|group5|group6}&withnode{y>=1.4 & y<=1.9};addsel'); % surface y=L (temporarily group 7)
femesh(';selelt selface&withoutnode{group2|group3|group5|group6}&withnode{y<=-1.4 & y>=-1.9};addsel'); % surface y=L (temporarily group 8)
femesh('join group6:9') % group 6 : inner surfaces
FEel0=FEelt;

femesh('set group1 mat 1 pro 1')
femesh('set group2 mat 2 pro 2 name q4p')
femesh('set group3 mat 3 pro 2 name q4p')
femesh('set group4 mat 4 pro 2 name q4p')
femesh('set group5 mat 5 pro 2 name q4p')
femesh('set group6 mat 6 pro 2 name q4p')

model=femesh('model');

% ------------------------------------------------------------------------
% properties
model.pl=[1 fe_mat('m_heat','SI',1) diffu rho heat alpha   % k rho C alpha  
          2 fe_mat('m_heat','SI',1) diffu rho heat alpha   % k rho C alpha
          3 fe_mat('m_heat','SI',1) diffu rho heat alpha   % k rho C alpha
          4 fe_mat('m_heat','SI',1) diffu rho heat alpha   % k rho C alpha
          5 fe_mat('m_heat','SI',1) diffu rho heat alpha   % k rho C alpha
          6 fe_mat('m_heat','SI',1) diffu rho heat alpha   % k rho C alpha
          ];
model.il=[1 fe_mat('p_heat','SI',1) 0 2 3;   % volume integration
          2 fe_mat('p_heat','SI',2) 0 2 3];  % exchange surface integration

% Constants
% Surface loading matrix : f(i,j) = \int_{\Gamma} N_i * N_j d\gamma

r1=struct('pl',[],'il',[]);
r1.pl=[1 fe_mat('m_heat','SI',1) 0 1 1 0;      % k rho C alpha  
          2 fe_mat('m_heat','SI',1) 0 1 1 0;    % k rho C alpha
          3 fe_mat('m_heat','SI',1) 0 1 1 0;    % k rho C alpha
          4 fe_mat('m_heat','SI',1) 0 1 1 1;    % k rho C alpha
          5 fe_mat('m_heat','SI',1) 0 1 1 0;    % k rho C alpha
          6 fe_mat('m_heat','SI',1) 0 1 1 0];   % k rho C alpha
r1.il=   [1 fe_mat('p_heat','SI',1) 0 2 3; ...   % volume integration
           2 fe_mat('p_heat','SI',2) 0 2 3];  % exchange surface integration
model=stack_set(model,'info','PLIL_ForEdgeLoad',r1);
   
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
out=model;if nargout==0;feplot(model);end

%% --------------------------------------------------------------------------
%% #HeatLoad load cases for heat demos
elseif comstr(Cam,'heatload');[CAM,Cam]=comstr(CAM,9);
    
if comstr(Cam,'stat')
 model=varargin{carg};
%% #HeatLoadStat : first example of static load -3
% Build the loading conditions on nodes of the interface
% This contains g+\alpha \theta_ext which will be integrated on the surface

% **** data ****
RO=struct('TOnFace1',25, ...
          'Alpha',1, ...
          'TExt',20, ...
          'F',40);

% **************

% boundary conditions on face x==0 (group 2)
data1 = struct('sel','group2',...
               'eltsel','group2',...
               'def',RO.TOnFace1,'DOF',.20);
model = fe_case(model,'Fvol','Surface_load1',data1);
% FVol on a surface does a a surface integral

% boundary conditions on face x==10
mdof0=feutil('findnode group3',model)+.20;
data2 = struct('DOF',mdof0,'def',ones(size(mdof0))*RO.TExt);
%data2=struct('DOF',.20,'def',RO.TExt); % xxx eb should test formal
model = fe_case(model,'DofSet','Surface_load2',data2);

% boundary conditions on other faces
node=feutil('getnode group4',model);
%data3 = struct('sel','group4','def',25 - node(:,5).^2/20,'DOF',node(:,1)+.20);
data3=struct('sel','group4','def','25-x.^2/20','DOF',.20);
model = fe_case(model,'Fvol','Surface_load3',data3);

% constant internal heat source
data4 = struct('sel','group1','def',[RO.F],'DOF',.20);
model = fe_case(model,'Fvol','Internal source',data4);
model=stack_set(model,'info','LoadParam',RO);

out=model;

elseif comstr(Cam,'pen')
%% #HeatLoadPen : same load as sdtweb('ofdemos','HeatLoadStat') -3
%  using "penalization" with high alpha rather than DOFset

mo3=varargin{carg};
mo3=ofdemos('heatLoadStat',mo3); % Original loading
RO=stack_get(mo3,'info','LoadParam','get');

%Surface_load2 : x=10 : turn dofset as high alpha surface
mo3=fe_case(mo3,'Remove','Surface_load2'); % remove DOFSet (imposed temp)
alpha=1e10; % High alpha as penalization
mo3.pl(mo3.pl(:,1)==3,6)=alpha; % change material property alpha
% build corresponding load:
data2 = struct('sel','group3',...
               'eltsel','group3',...
               'def',RO.TExt*alpha,'DOF',.20);
mo3 = fe_case(mo3,'Fvol','Surface_load2',data2);
out=mo3;

%% 
else;error('HeatLoad%s unknown',CAM); 
end


%% --------------------------------------------------------------------------
%% #AnnularHeat 2D annular heat demo by S. Pagano
elseif comstr(Cam,'annularheat')
    
    
femesh('reset');FEnode=[1 0 0 0 .5 0 0;2 0 0 0 1 0 0];
femesh('objectbeamline',[1 2]);
femesh(';divide20;rev 60 o 0 0 0 360 0 0 1');
if exist('RunOpt','var')&&isfield(RunOpt,'tria')&&RunOpt.tria
    femesh(';quad2tria;addsel;set group1 name t3p;');
else;femesh(';addsel;set group1 name q4p;'); %#ok<NOSEM>
end
femesh(';selelt seledge;addsel;set group 2 pro 2 name bar1');
model=femesh('model');
model.pl=[1 fe_mat('m_heat','SI',1) 1 0 0 1e-10   % k rho C alpha  
%          2 fe_mat('m_heat','SI',1) 400 1 1 1e10  % k rho C (alpha_0: not used)
          ];
model.il=[1 fe_mat('p_heat','SI',1) 0 -3 2;   % volume integration
          2 fe_mat('p_heat','SI',2) 0 -3 2];  % exchange surface integration


R1=struct('TOnFace1','-exp(x)/0.5*(cos(y)*x-sin(y)*y)', ...  % 1
          'Alpha',1, ...
          'TExt','exp(x)*cos(y)', ...
          'F',0);

% boundary conditions on face x==0  COND NEUMANN
node = feutil('getnode r==0.5',model);
dataD_1 = struct('sel','eltname bar1 & innode {r==0.5}',...
            'def',-exp(node(:,5)).*(cos(node(:,6)).*node(:,5) -...
             sin(node(:,6)).*node(:,6))/0.5,'DOF',node(:,1)+.20);
model = fe_case(model,'Fvol','Surface_load1',dataD_1);


% boundary conditions on other faces r ==1.0 COND DIRICHLET
node = feutil('getnode r==1.0',model);
dataD_1 = struct('sel','eltname bar1 & innode {r==1.0}','def',exp(node(:,5)).*cos(node(:,6)),'DOF',node(:,1)+.20);
model = fe_case(model,'DofSet','Surface_load3',dataD_1);

% constant internal heat source
dataV_1 = struct('sel','group1','def',[R1.F],'DOF',.20);
model = fe_case(model,'Fvol','Internal source',dataV_1);

[Case,model.DOF]=fe_mknl('init -gstate',model);
% sp_util('diag',12);  % very slow elem0.m file implementation, just check
% the right-hand side is in dc.def(:,2) which is modified even though it is
% not returned by fe_mknl
sp_util('diag',0)
dc = struct('def',zeros(length(model.DOF),2),'DOF',model.DOF);
Case.GroupInfo{2,7}.data=full(Case.GroupInfo{2,7}.data);
k  = fe_mknl('assemble not',model,Case,dc,1);
sp_util('diag',0); % modified here

% display GroupInfo on edge
cEGI=feutil('findelt eltname bar1',model);
r1=Case.GroupInfo{2,7};
% computation of the response to internal forces and enforced displacement
[adof,ind,c]=fe_c(model.DOF,Case.DofIn);
cind = setdiff(1:length(model.DOF),ind);
resp=struct('def',zeros(size(model.DOF)),'DOF',model.DOF);
resp.def(ind,1)=Case.TIn(ind,:);
F=dc.def(cind,2)-k(cind,ind)*Case.TIn(ind,:);
resp.def(cind,1)=ofact(k(cind,cind),F);

% Plot
[i1,n1]=feutil('findnode r==0.5',model);
r1=fe_c(resp.DOF,i1)*resp.def;
r1(:,2)=exp(n1(:,5)).*cos(n1(:,6)); % reference analytic solution
figure(1);h=plot(atan2(n1(:,5),n1(:,6))/pi*180,r1,'+');
legend('FEM','Analytic');set(h(2),'marker','o')
if norm(diff(r1,[],2))/norm(r1(:,1))>1e-3;error('Mismatch');end

feplot(model,resp);fecom('colordata20')

%% --------------------------------------------------------------------------
%% #ThermalCube pres-stress of a cube demo by E. Balmes, tested with SDT only
elseif comstr(Cam,'thermalcube')

    
% RivlinCube Expansion - - - - - - - - - - - - - - - - - - -
model=femesh('test hexa8b divide 3 3 3');
model=fe_case(model,'reset','FixDof','x','x==0 -dof 1', ...
   'FixDof','y','y==0 -dof 2','FixDof','z','z==0 -dof 3');
model.pl=m_elastic('dbval 100 steel');
model.pl(8:9)=[1.2e-2 20]; % very large expansion for test purposes

if sp_util('issdt')  % High level NL static implementation in SDT

 defT=struct('sel','groupall','dir',{{30}},'DOF',.20);
 model=feutil('setpro 111',model,'MAP',defT);
 opt=fe_simul('nlstaticopt');opt.RelTol=-1e-3;
 dc=fe_time(opt,model);
 feplot(model,dc); % dc.def(:,3)=d2.def; fecom('ch1 3')
 if nargout>0; [Case,model.DOF]=fe_mknl('init',model);end
else % Low level OpenFEM implementation

 % Uniform temperature of 30 degrees (nominal 20)
 defT=struct('sel','groupall','dir',{{30}},'DOF',.20);
 model=fe_case(model,'DofSet','ThermalState',defT);
 model.DOF=feutil('getdof',model);

 % Low level illustration with all steps done
 dc=struct('DOF',model.DOF,'def',zeros(length(model.DOF),2));
 [Case,model.DOF]=fe_mknl('init',model);
 Case.GroupInfo{7}.data(2,:)=0; Case.GroupInfo{7}.lab={'T','p'};% check multi-field
 [k1,C2,dc]=fe_mknl('assemble NoT',model,Case,dc,5);
 dc.def(:,1)=ofact({k1,C2.T},-dc.def(:,2));
 dc.lab={'Deformation','Internal Load'};
 feplot(model,dc);fecom('ch2')

end
if nargout>0; out=model;out1=dc;out2=Case;end
%% --------------------------------------------------------------------------
% #Composite Simple composite plate
elseif comstr(Cam,'composite')

% Mesh the geometry
model=femesh(strcat('testquad4 struct divide 2 10 back'));
% old init should not be an issue
%model=femesh('testquad4');
%model=feutil('divide-OLD 2 10',model);
model.Node(:,5:6)=model.Node(:,5:6)*diag([.15 .01]);

if ~isempty(strfind(Cam,'tria'));
 model.Elt=feutil('quad2tria',model.Elt);
 if ~isempty(strfind(Cam,'tria6'));model=feutil('lin2quad',model);end
 % for triangles, an element orientation field is needed
 data=struct('sel','groupall','dir',{{1,0,0}},'DOF',[.01;.02;.03]);
 model=p_shell('setTheta',model,data);model.name='OrientedTriangles';
elseif ~isempty(strfind(Cam,'quadb'));
 model=feutil('lin2quad',model);
 % for triangles, an element orientation field is needed
 data=struct('sel','groupall','dir',{{1,0,0}},'DOF',[.01;.02;.03]);
 model=p_shell('setTheta',model,data);model.name='quadb';
else % 
 % for quadrangles, in this mesh, the element orientation is 0 for all Elts
 model.name='quad4';
end

% single ply with MAT100, thickness 1 and orientation 30 deg
model.il = p_shell('dbval 110 laminate 100 1 30'); 

% sample orthotropic material
model.pl=[100 fe_mat('m_elastic','SI',5) ...
      38.095e9  9.4e9 0.3   7.5e9 0    0     1640]; 
model.pl(7:8)=model.pl(6);fprintf('Setting G1z=G2z=G12\n'); 

% Old for reference
%model.pl=m_elastic('dbval 100 steel','dbval 110 steel', ...
%  'dbval 200 lamina .27 3e9 .4 1200 0  790e9 .3 1780 0');
  
model=fe_case(model,'fixdof','base','x==0');

out=model; 
%% --------------------------------------------------------------------------
% #MatInterp : position dependent material properties
elseif comstr(Cam,'matinterp')

model=femesh('test hexa8b divide 1 1 8');
model=fe_case(model,'reset','FixDof','x','x==0 -dof 1', ...
   'FixDof','y','y==0 -dof 2','FixDof','z','z==0 -dof 3');
mat=struct('pl',feutil('getpl 100',model), ...
    'type','m_elastic','unit','SI', ...
    'E',struct('X',[0;max(model.Node(:,7))], ...
         'Xlab',{{'T'}},'Y',210e9*[1;20]));
model=stack_set(model,'mat','Steel',mat);
model=fe_case(model,'DofSet','ThermalState', ...
    struct('DOF',model.Node(:,1)+.20,'def',model.Node(:,7)));
[SE,CE]=fe_case(model,'assemble -MatDes 2 5 -SE -NoT');

% Verify z variation of color
i1=fe_c(SE.DOF,.03,'ind');
d1=struct('def',diag(SE.K{2}(i1,i1)),'DOF',SE.DOF(i1));
i1=fe_c(d1.DOF,feutil('findnode z==0|z==1',model),'ind');
d1.def(i1)=d1.def(i1)*2; % top bottom counted once, *2 for proper color
cf=feplot(SE,d1);
fecom colordataevalz

% --------------------------------------------------------------------------
%% #Patch9 : 9 element patch
elseif comstr(Cam,'patch9')
% Prepare a 4 by 2 quad plate
mo1=femesh('testquad4 divide 4 2');
% Objective: add a central element in the plate
% Get central node
n1=feutil('getnodegroupall',mo1);
r1=mean(n1(:,5)); r2=mean(n1(:,6));
n1=feutil('findnode',mo1,sprintf('x==%.15g & y==%.15g',r1,r2));
% select elements connected to central node
% one will treat separately elements connected to n1 to add a central element
mo2=mo1; [mo1.Elt,mo2.Elt]=feutil('RemoveElt withnode',mo2,n1);
mo1.Node=feutil('GetNodeGroupAll',mo1); % clean nodes
mo2.Node=feutil('GetNodeGroupAll',mo2); % clean nodes
[mo2,n2]=feutil('unjoin',mo2,sprintf('withnode{x<%.15g}',r1),...
 sprintf('withnode{x>%.15g}',r2));
% n2 should give three unjoined nodes in the center along y axis
% distinguish the central node and the side nodes
i1=ismember(n2(:,1),n1);  i2=find(~i1); i1=find(i1);
% Modify coordinates of duplicated central nodes to open the new element
NN=sparse(mo2.Node(:,1),1,1:size(mo2.Node,1)); % NodeId to NodeInd matrix
mo2.Node(NN(n2(i1,1)),5)=mo2.Node(NN(n2(i1,1)),5)-.15; % open
mo2.Node(NN(n2(i1,2)),5)=mo2.Node(NN(n2(i1,2)),5)+.25; % open
% Add the central element, quad4
mo2=feutil('addelt',mo2,'quad4',[n2(i2(1),1) n2(i1,1) n2(i2(2),1) n2(i1,2)]);
mo2.Elt=feutil('Orient;',mo2); % check clean orientation of new element
% Merge modified elements to old model
mo1=feutil('AddTestMerge-NoOri;',mo1,mo2);
% Cleanup model
mo1=feutil('JoinAll',mo1); % set one single element group
mo1.Elt=feutil('SetGroupAll mat 100 pro 110',mo1); % clean mat/pro assignment
mo1.Elt=feutil('orient 1 n 0 0 100;',mo1);
%mo1=feutil('renumber',mo1,1:size(mo1.Node,1)); % reorder nodes
out=mo1;

if 1==2 % this is the old version
model=femesh('testquad4');
model=feutil('divideElt-OLD 4 2',model);
model.Node(10,5)=0.25;
model.Node(16,:)=[16 0 0 0 0.75 0.5 0];
model.Elt(7+1,1)=16;
model.Elt(6+1,4)=16;
model.Elt(9+1,1:6)=[9 16 11 10 100 110];
model.Node(10,5:6)=model.Node(10,5:6)+[0.1 0];
out=model; 
end
%% --------------------------------------------------------------------------
%% #Patch4 : 4 element patch
elseif comstr(Cam,'patch4')

model=femesh('testquad4 divide 2 2');
% Get central node and move it
n1=feutil('getnodegroupall',model);
r1=mean(n1(:,5)); r2=mean(n1(:,6));
n1=feutil('findnode',model,sprintf('x==%.15g & y==%.15g',r1,r2));
NN=sparse(model.Node(:,1),1,1:size(model.Node,1)); % NodeId to NodeInd matrix
model.Node(NN(n1),5:6)=[1/3 2/3];
% Orient elements
MAP=feutil('getnormalElt MAP -dir1',model);
MAP.normal(:,1)=0;MAP.normal(:,2)=1;MAP.normal(:,3)=0;
model=p_shell('setTheta -strategy 2',model,MAP);
out=model; 

if 1==2 % this is the old version
 model=femesh('testquad4');
 model=feutil('divide-OLD 3 3 ',model);
 model.Node(13,5:6)=[0.5 0.75];
 model.Elt([2 3 4 5 8],:)=[];
 MAP=feutil('getnormalElt MAP -dir1',model);
 MAP.normal(:,1)=0;MAP.normal(:,2)=1;MAP.normal(:,3)=0;
 model=p_shell('setTheta -strategy 2',model,MAP);
 out=model;
end


%% --------------------------------------------------------------------------
%% #LS : level set tests
elseif comstr(Cam,'ls'); [CAM,Cam]=comstr(CAM,3);
    
if nargin>=carg;RO=varargin{carg};else;RO=[];end
if comstr(Cam,'2d');[CAM,Cam]=comstr(CAM,3);
    %% #Ls2D
    r1=d_shm('MeshCfg',struct,'PimmShellA');r1=r1.param.MeshCfg.data{1};
    r2=r1;r2.list(3:end,:)=[];
    if ~isempty(RO); %initial dim [400 300 8]
      r2.list{2,3}.lx=RO.dim(1);r2.list{2,3}.ly=RO.dim(2);r2.list{2,3}.lc=RO.dim(3);
    end
    mdl=d_piezo('MeshPlate',r2);li=r1.list(3:end,3);
    for i1=1:length(li);li{i1}.mpid=i1*[100 1000];end
elseif comstr(Cam,'3d');[CAM,Cam]=comstr(CAM,3);
    %% #Ls3D
    if ~isempty(RO); mdl=femesh(sprintf('testhexa8 divide %d %d %d',RO.dim));
    else;mdl=femesh('testhexa8 divide ');%initial divide 3 3 3
    end
    li={struct('shape','sphere','xc',.0,'yc',.0,'zc',0.3,'rc',.2,'mpid',[200 300])};
 %% add an OrientMap
    eltid=feutil('eltidfix;',mdl);
    r1=eltid;r1(:,2)=1; 
    cg=feutil('getcg',mdl);cg(:,3)=0;
    cg=diag(sparse(1./sqrt(cg(:,1).^2+cg(:,2).^2)))*cg;
    r1(:,[7:8 10 11])=[cg(:,2) -cg(:,1) cg(:,1) cg(:,2)];% Horizontal circle
    i1=eltid~=0;r1(1,12)=0;
    r1=struct('EltId',eltid(i1),'basid',eltid(i1),'bas',r1(i1,:));
    mdl=stack_set(mdl,'info','EltOrient',r1);
    %feplot(mdl);fecom('showmap','EltOrient');
    
else
      error('LS%s unknown.',CAM);
end
[eltid,mdl.Elt]=feutil('eltidfix;',mdl);
out=mdl; out1=li;

%% --------------------------------------------------------------------------
%% #Hyper : test for hyper-elasticity
elseif comstr(Cam,'hyper')

dbstack; keyboard;% xxx copy from HBV17

%% #end ----------------------------------------------------------------------
elseif comstr(Cam,'cvs');
    out='$Revision: 1.33 $  $Date: 2025/04/07 17:08:17 $';
else error('%s unknown',CAM);
end
