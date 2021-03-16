function out=t_fe_time(varargin);

if nargin==0
   
  % currently tested with SDT only
  if isempty(which('ii_fin'));cd(fileparts(fileparts(which('fe_mk'))));eval('sdtcheck path');cd ..;end
  t_fe_time('old')
  
  if isempty(which('ii_fin'));cd(fileparts(fileparts(which('fe_mk'))));eval('sdtcheck path');cd ..;end
  t_fe_time('sdt')
  t_fe_time('theta')
  t_fe_time('explicitnewmark')
  return
  
end

%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*DEFNU,*NOSEM>
[CAM,Cam]=comstr(varargin{1},1);carg=2;
%%  ------------------------------------------------------------
if comstr(Cam,'theta') % Test of theta method
    
model=fe_time('demo bar'); q0=[];
data=struct('DOF',2.01,'def',1e6,...
  'curve',fe_curve('test ricker 10e-4 101 1 100e-4'));
model = fe_case(model,'DOFLoad','Point load 1',data);

%model.alpha (ancienne valeur de theta)
opt=struct('Opt',[.5 0 0 1e-4 1000],'Method','theta');
def=fe_time(opt,model);

if 1==2
 figure(1);plot(sum(def.def))
 feplot(model,def);fecom('ColorDataEvalX')
end

%sdtweb o:/vermot/cur_doc/R10_AcouFren_Lot1-2.pdf


%% #ExplicitNewmark ----------------------------------------------------------
elseif comstr(Cam,'explicitnewmark')

 model=fe_time('demo bar10'); % build the model
 % set the time options in model.Stack
 model=fe_time('TimeOpt Newmark 0 .5 0 1e-4 100',model);
 def=fe_time(model); % compute the response
  

%% #Old ------------------------------------------------------------
elseif comstr(Cam,'old')
% - with/without output time vector (OutputFcn):
% - output DOFs (OutInd) :
% - u, v, or a output :
% - observation matrix associated to u, v or a :
% - output evaluation string (OutputFcn with OutputInit) :

model=fe_time('demo bar'); q0=[];
data=struct('DOF',2.01,'def',1e6,...
  'curve',fe_curve('test ricker 10e-4 101 1 100e-4'));
model = fe_case(model,'DOFLoad','Point load 1',data);

model=fe_time('timeopt newmark .25 .5 0 1e-4 10000',model);

TimeOpt=stack_get(model,'info','TimeOpt','getdata');


% - with output time vector (OutputFcn):
r1=TimeOpt;r1.OutputFcn=[1e-4 5*1e-4 1e-3]; 
model=stack_set(model,'info','TimeOpt',r1);
def=fe_time(model);
if ~isequal(r1.OutputFcn,def.data(:)');
 error('OutputFcn : steps do not correspond');
end

% - output DOFs (OutInd) :
r1=TimeOpt;r1.OutInd=1:10;
model=stack_set(model,'info','TimeOpt',r1);
def=fe_time(model);
% - u, v, or a output :
% XXX 

% - observation matrix associated to u, v or a :
% XXX

% - output evaluation string (OutputFcn) :
r1=TimeOpt;r1.OutputFcn='out.def(:,j1+1)=0*u;'; 
model=stack_set(model,'info','TimeOpt',r1);
def=fe_time(model);
if norm(def.def)~=0; error('OutputFcn not taken into account'); end


%    without OutputInit :
model=stack_set(model,'info','TimeOpt',r1);
def=fe_time(model);
%    with OutputInit :
r1.OutputInit='accel=zeros(length(mdof),nt);out=struct(''DOF'',mdof,''def'',accel,''data'',t(:),''fun'',[0 4]);';
model=stack_set(model,'info','TimeOpt',r1);
def=fe_time(model);

% test of_time newmarkinterp
op2=fe_time('TimeOpt Newmark .25 .5 0 1e-4 100');
dt=op2.Opt(4); t=op2.Opt(3)+dt*[0:op2.Opt(5)]; 
op2.OutputFcn=t(2:end)-dt/4;op2.OutInd=1; def=fe_time(op2,model);
 
op3=op2;op3.OutputFcn=t(1:end)+dt/4;d2=fe_time(op3,model);
% shifted -+ quarter step
r2=of_time('lininterp',[d2.data d2.def(1,:)'],def.data,zeros(1,3));
figure(1); plot(def.data,def.def(1,:),'r-+',d2.data,d2.def(1,:),'b--', ...
    def.data,r2,'bo');
if (norm(def.def)-norm(r2))/norm(d2.def)>.1;
    error('on Newmark interp');
end

% select U, V, A
TimeOpt.NeedUVA=[1 1 0];
def=fe_time(TimeOpt,model);

% c_u
Case=fe_case('gett',model);
i1=feutil('findnode x>30',model);
TimeOpt=fe_time('TimeOpt Newmark .25 .5 0 1e-4 100');
TimeOpt.c_u=fe_c(Case.DOF,i1+.01);
def=fe_time(TimeOpt,model);


% curves
model=femesh('test hexa8');
model=fe_case(model,'FixDOF','base','z==0');
TimeOpt=struct('Method','Newmark','Opt',[.25 .5 0 1e-4 200]);
model=stack_set(model,'info','TimeOpt',TimeOpt);

c1=fe_curve('test ricker 10e-4 101 1 100e-4');
c2=fe_curve('test ricker 20e-4 101 1 100e-4');
model=fe_case(model,'DofLoad','ecraz',...
  struct('DOF',[5 6 7 8]'+0.03,'def',[-1 -1 -1 -1]',...
  'curve',c1)); % Single load
def=fe_time(model);

% --------------------------------------------------------- HHT
model=fe_time('demo bar');
TimeOpt=fe_time('TimeOpt hht .05 Inf 0 1e-4 100');
TimeOpt.NeedUVA=[1 1 0];
def=fe_time(TimeOpt,model); 
% compare with newmark
TimeOpt=fe_time('TimeOpt newmark .25 .5 0 1e-4 100');
TimeOpt.NeedUVA=[1 1 0];
def2=fe_time(TimeOpt,model); 
figure;plot([def2.v(:,end) def.v(:,end)]);
legend('NM','HHT');
% ---------------------------------------------------------


%% 
%% #SDT Below are SDT based tests---------------------------------------------
elseif comstr(Cam,'sdt')
 if mkl_utils&& (isempty(gcbf) ||  ...
   (~isempty(gcbf)&&strcmp(get(gcbf,'tag'),'test_sdt')));
  if isempty(which('ii_fin'))
   wd=pwd;
   cd(fileparts(fileparts(which('fe_mk'))));eval('sdtcheck path');cd ..
   clear functions
  end

%% First attempt at energy computations
model=fe_time('demo bar'); q0=[];
data=struct('DOF',2.01,'def',1e6,...
  'curve',fe_curve('test ricker dt=10e-4 A=1e-4'));
model = fe_case(model,'DOFLoad','Point load 1',data);

TimeOpt=fe_time('timeopt nlnewmark .25 .5 0 1e-4 30');
TimeOpt.nf=1;
TimeOpt.Residual='r=-fc; [EM,EK]=mkl_utils(''residual'',r,model.K,u,v,a,opt.Rayleigh);assignin(''caller'',''EM'',EM);assignin(''caller'',''EK'',EK);';
TimeOpt.Rayleigh=[0 0];%[1e-6 0.1];
TimeOpt.OutputFcn= 'out.def(:,j1+1)=u;out.v(:,j1+1)=v;if exist(''EM'',''var'');out.data(j1+1,2:3)=[EM EK];end';
def=fe_time(TimeOpt,model);
figure(1);plot(def.data(:,2:3))

%---------------------------------------------------------------------

d_ubeam
q0=[];
data=struct('DOF',[349.01;349.02],'def',1e6*eye(2),...
  'curve',{{fe_curve('test ricker .005     101 1 .1'),...
            fe_curve('test ricker .005     101 1 .1')}});
%                              duration (200 Hz)    total time
% 1st mode : 75 Hz
model = fe_case(model,'DOFLoad','Point load 1',data);

mo1=fe_mknl(model);
Case=fe_mknl('init',mo1);
phi=fe_eig(model);
phi.def=Case.T'*phi.def(:,1:3);phi.data=phi.data(1:3);

mo1=stack_set(mo1,'case','Case 1',Case);
mo1.K{3}=mo1.K{2};
mo1.K{2}=struct('def',mo1.K{1}*phi.def,'data',2*2*pi*[0.01*75 0.02*127 0.03*155]);

mo1.Klab={'m','c','k'};
mo1.Opt=[1 0 0;2 3 1];

if 1==2
TimeOpt=fe_time('timeopt newmark .25 .5 0 1e-3 1e3'); % 10 sec
TimeOpt.Residual='r = (ft(j1,:)*fc''-v''*(c.def*diag(c.data)*(c.def)'')-u''*k-a''*model.K{1})'';';
TimeOpt.Residual='r0=(ft(j1,:)*fc'')''; mkl_utils(''residual'',r0,model.K,u,v,a);';
TimeOpt.Residual='r=-(ft(j1,:)*fc'')''; mkl_utils(''residual'',r,model.K,u,v,a);r=-r;';
else
TimeOpt=fe_time('timeopt nlnewmark .25 .5 0 5e-4 5e3'); % 10 sec
TimeOpt.RelTol=-1e-6;
TimeOpt.Residual='r = model.K{1}*a+model.K{3}*u-fc;';
TimeOpt.Residual='r0=(ft(j1,:)*fc'')''; mkl_utils(''residual'',r0,model.K,u,v,a);';
TimeOpt.Residual='r=-fc; mkl_utils(''residual'',r,model.K,u,v,a);';
%TimeOpt.Residual='r = model.K{1}*a+model.K{2}*v+model.K{3}*u-fc;';
end
% Formula to evaluate the periodicity error : w^2*dt^2/12 : (127*2*pi)^2*1e-6/12

def=fe_time(TimeOpt,mo1);
ind=fe_c(def.DOF,[349.01],'ind');figure(1);subplot(2,1,1);plot(def.def(ind,:))

% select ouput => iiplot
f=[0:length(def.data)-1]/length(def.data)/diff(def.data(1:2));
figure(1);Y=fft(def.def(ind,:));subplot(2,1,2);semilogy(f,abs(Y));ii_plp(phi.data)

% no damping
TimeOpt=fe_time('timeopt nlnewmark .25 .5 0 5e-4 5000'); % 10 sec
TimeOpt.RelTol=-1e-6;
def2=fe_time(TimeOpt,model);

% select ouput => iiplot
f2=[0:length(def.data)-1]/length(def2.data)/diff(def2.data(1:2));
Y2=fft(def2.def(ind,:));

figure(4);in1=1:1000;semilogy(f(in1),abs(Y(in1)),f2(in1),abs(Y2(in1)))
ii_plp(phi.data)

end %~isempty(gcbf)&&strcmp(get(gcbf,'tag'),'test_sdt');

% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------
elseif 1==2
% ----------------------------------------------------------------------------
% to test extrapolation of curves

   curve_1=fe_curve('TestRicker 10e-4 120 3. 1e-3'); % here curve to extrapolate
   model=fe_curve(model,'set','ricker_1',curve_1)

   model = fe_case(model,'AddToCase 1','KeepDof','kd',.01);
   model = fe_case(model,'AddToCase 1','FixDof','fd',1);
   model=fe_mk(model);
   Case=fe_case(model,'getcase');

   def=fe_time('newmark .25 .5 0 1e-4 100',model,Case,[]);
   def_v=def;def_v.def=def_v.v; def_v.DOF=def.DOF+.01;
   feplot(model,def_v); fecom(';view2;animtime;ch30;scd3');
%% ----------------------------------------------------------------------------
else; error('%s unknown',CAM)
end
