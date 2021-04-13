function [out,out1,out2,out3] =fe_time(varargin)

%FE_TIME computes time integration with Newmark, HHT, Euler, Theta,
%        and Discontinuous Galerkin schemes
%
% Syntax : def=fe_time(model)
%          def=fe_time(TimeOpt,model)
%          model=fe_time('TimeOpt...',model)
%          TimeOpt=fe_time('TimeOpt...')
%          [def,model,opt]=fe_time(model)
%
%       Accepted solvers : newmark (damping supported), 
%                          dg (no damping supported),
%                          nlnewmark (damping supported) : newmark + newton 
%                          staticnewton
%                          HHT-alpha (through newmark scheme)
%                          theta (linear implmentation)
%       TimeOpt is a data structure with fiels
%
%                  .Method : 'nlnewmark', 'dg' or 'newmark'
%                  .Opt :  [beta gamma t0 deltaT  Nstep Nf Thres]
%                  .MaxIter : maximum number of iterations
%                  .Jacobian (only for non-linear). The default jacobian is
%                   'ki=ofact(model.K{3}+2/dt*model.K{2}+4/dt^2*m);'
%                  .JacobianUpdate (only for nlnewmark) : default is 0
%                     0 modified newmark (no update in Newton iterations)
%                     1 update in Newton iteration
%                  .Residual (only for non linear). The default residual is
%                      'r = model.K{1}*a+model.K{2}*v+model.K{3}*u-fc;'
%                  .InitAcceleration : optional field to specify
%                           a particular method to init
%                           the initial guess of acceleration field
%                  .OutputFcn  String or time vector of output steps
%                  .SaveTimes  : optionnal time vector, saves time steps
%                  .TimeVector (if exists TimeVector is used 
%                                instead of time parameters)  
%                  .RelTol (default getpref('OpenFEM','THRESHOLD',1e-6))
%                  .NeedUVA [NeedU NeedUV NeedA] Need u,v,a? 0 if no, 1 if yes
%
%       Other input arguments are :
%
%       q0    : initial conditions
%               data structure containing 1 or 2 fields
%               fiels are : .DOF specifying dofs
%                           .def and .v specifying displacements and velocity 
%                                       at t=0
%               if q0 is empty, zeros initial conditions are taken
%
%       See sdtweb      fe_time, fe_mk, fe_load, fe_case
%       See also help   fe_mk, fe_case

%	Jean-Michel Leclere, Etienne Balmes, Jean-Philippe Bianchi, 
%	Guillaume Vermot des Roches
%       Copyright (c) 2001-2020 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       Use fe_time('cvs') for revision information

RunOpt=struct('Out','Standard');
%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*NOSEM>

% - - - - - - - - - - - - - - - - - - -  - get command
if ischar(varargin{1})

 [CAM,Cam]=comstr(varargin{1},1);carg=2;
 opt=struct('Method','','Opt',[]);
 if comstr(Cam,'timeopt'); 
   [CAM,Cam]=comstr(CAM,8);RunOpt.Out='TimeOpt';
   [CAM,Cam,i1]=comstr('reset',[-25 3],CAM,Cam);if i1; RunOpt.Out='TimeOptReset'; end
   [CAM,Cam,i1]=comstr('set',[-25 3],CAM,Cam);if i1; RunOpt.Out='TimeOptSet'; end
   [st,i1,un1,i2]=sscanf(CAM,'%s',1); 
   if exist(st,'file');% fe_time('timeopt nl_solve NlNewmark')
       out=feval(st,['TimeOpt' CAM(i2:end)],varargin{carg:end});return;
   end
 end
 if comstr(Cam,'newmark')    
   opt.Method='Newmark'; [CAM,Cam]=comstr(CAM,8);
 elseif comstr(Cam,'cvs');
  out='$Revision: 1.353 $  $Date: 2021/04/12 16:52:14 $';return;
 elseif comstr(Cam,'nlnewmark') 
   opt.Method='NLnewmark'; [CAM,Cam]=comstr(CAM,10);
 elseif comstr(Cam,'hht');
   opt.Method='NLnewmark'; opt.HHTalpha=5e-2; [CAM,Cam]=comstr(CAM,4);
 elseif comstr(Cam,'dg')     
   opt.Method='dg';     [CAM,Cam]=comstr(CAM,3);
 elseif comstr(Cam,'staticnewton')     
   opt.Method='staticnewton';     [CAM,Cam]=comstr(CAM,12);
 elseif comstr(Cam,'static');       
  opt.Method=Cam; [CAM,Cam]=comstr(CAM,7);
 elseif comstr(Cam,'theta');
  opt.Method='theta'; opt.theta=.5; [CAM,Cam]=comstr(CAM,6);
 %% #BackModel : modify init acceleration to allow return of assembled model
 % opt.Method='back' is the preferred methodology
 elseif comstr(Cam,'backmodel')
  model=[]; opt=varargin{carg}; carg=carg+1;
  if isfield(opt,'Elt');
      model=opt;opt=stack_get(model,'info','TimeOpt','get');
  end
  if strcmpi(opt.Method,'nlnewmark')
    opt.InitAcceleration= ...
        ['model.K=feutilb(''MakeSparse'',model.K);out=model;' ...
         'if nargout>1;out1=Case;end;return;'];
  else; error('Back model not yet implemented for %s',opt.Method)
  end
  if isempty(model); out=opt;
  else;[out,out1]=fe_time(opt,model);
  end
  return
 %------------------------------------------------------------------------
 %------------------------------------------------------------------------
 %% #Demo --------------------------------------------------------------------
 elseif comstr(Cam,'demo'); [CAM,Cam]=comstr(CAM,5);
  
  if comstr(Cam,'bar'); 
  %% #DemoBar : simple bar example ------------------------------------------
   % Run with fe_time('demobar -run');
   [CAM,Cam]=comstr(CAM,4);
   [CAM,Cam,RunOpt.run]=comstr('-run',[-25 3],CAM,Cam);
   opt=comstr(CAM,[-1 500]); % number of elements
   
   model=struct('Node',[1 0 0 0    0   0  0;2 0 0 0    50  0  0], ...
       'Elt',feutil('addelt','bar1',[1 2 1 1 1]), ...
       'pl',[1 fe_mat('m_elastic',1,1) 2.1e10 0 2500 0], ...
       'il',[1 fe_mat('p_beam','SI',1) 0 0 0 1]);
   model=feutil(sprintf('divide%i',opt(1)),model);   
   
   % #LoadCurve Define load and time dependence   -2
   data=struct('DOF',2.01, ...  % Application DOF
       'def',1e6, ...           % Amplitude
       'curve','input');
   model = fe_case(model,'DOFLoad','Point load 1',data);
   %model = fe_curve(model,'set','input','testTestRicker dt=1e-2 A=2');
   model = fe_curve(model,'set','input','teststep 10e-4');

   % Define boundary conditions
   model = fe_case(model,'FixDof','Edge_and_1D',[1;.02;.03;.04;.05;.06]);
   model=stack_set(model,'info','Rayleigh',[1 0  0 1e-5]); % group 1 beta=.01
   model.name='demobar'; model.unit='SI';
   out=model; out1=fe_case(model,'getcase');
   if RunOpt.run
    opt=fe_time('TimeOpt Newmark .25 .5 0 1e-4 400');
    def=fe_time(opt,model); % compute the response
    def.DOF=def.DOF+0.02;
    if nargout==0; 
     cf=feplot(model,fe_def('subdef',def,'1:10:end'));fecom(';showFiCEvalz;chend')
    else; out=model;out1=def;
    end
   end
   return; 
   
  elseif comstr(Cam,'2d') 
  %% #Demo2D : simple 2D square ----------------------------------------------
   FEnode=[1 0 0 0    0   0  0; 2 0 0 0    10  0  0];
   FEel0=[Inf abs('bar1')]; FEel0(2,:)=[1 2 1 1 1];
   pl=[1 fe_mat('m_elastic','SI',1) 1e10 0. 2500 0]; 
   il=[1 fe_mat('p_solid','SI',2) 1 0 -3];
   femesh('extrude 1 0 10 0'); femesh(';divide100 100;quad2tria');
   femesh('set groupa1 name t3p matid1 proid1')
   model=struct('Node',FEnode,'Elt',FEel0,'pl',pl,'il',il);
   [nodeID,node]=feutil('findnode x==0 & y==0',model);
   data=struct('DOF',[nodeID+.01 nodeID+.02]','def',[1e6 1e6]');
   data.curve=struct('X',linspace(0,10e-4,10)','Y',ones(10,1));
   model = fe_case(model,'DOFLoad','Point load 1',data);
   model = fe_case(model,'FixDof','fdX','x==0 & y>0 -DOF 2');
   model = fe_case(model,'FixDof','fdY','y==0 & x>0 -DOF 1');
   model.name='demo2D'; model.unit='SI';
   out=model; return;  
   
  else;error('fe_time demo : unknown');
  end
 %------------------------------------------------------------------------
 %------------------------------------------------------------------------
 %------------------------------------------------------------------------
 %% #OptimProduct ------------------------------------------------------------
 elseif comstr(Cam,'optimproduct'); [CAM,Cam]=comstr(CAM,13);

  % [ind,flag]=fe_time('optimproduct',k); % flag=1 => transpose

  k=varargin{carg}; carg=carg+1;
  
  if comstr(version,'7');
   st1={'','symamd','colamd','colperm','symrcm'};
  else;st1={'','symamd','colmmd','colperm','symmmd','symrcm'};
  end

  st=''; t1=Inf; i3=10; u=ones(length(k),1);ut=u';

  for j1=1:length(st1);
   if isempty(st1{1,j1}); i2=1:length(k);else; i2=feval(st1{1,j1},k); end
   k2=k(i2,i2);
   while 1
    t=cputime; for j2=1:i3; r1=ut*k2; end;st1{3,j1}=cputime-t;
    if st1{3,j1}>1; break; else; i3=i3*10;end
   end
   t=cputime; for j2=1:i3; r1=k2*u; end;st1{2,j1}=cputime-t;
  end
  r2=[[st1{2,:}];[st1{3,:}]]';r3=min(r2(:));
  r2=(r2-ones(size(r2,1),2)*min(r3))/r3*100;
  [r3,i2]=min(r2(:)); 
  out1=0; if i2>size(r2,1); out1=1; i2=i2-size(r2,1);end % transpose
  j1=i2;
  for j1=1:size(r2,1); 
   fprintf('%10s %6.1f %6.1f %% slower\n',st1{1,j1},r2(j1,:));
  end
  if isempty(st1{1,j1}); i2=1:length(k);else; i2=feval(st1{1,j1},k); end
  out=i2; out2=st1{1,j1};
  return;
 
 %% --------------------------------------------------------------------------
 %% #Set : allow property editing on existing options
 elseif comstr(Cam,'set'); sdtw('_ewt','Obsolete should use TimeOptMethod ...');
  opt=varargin{carg};carg=carg+1;
  RO=varargin{carg};carg=carg+1;
  if isfield(RO,'dt') % Adjust dt
    if ~isfield(RO,'tend');RO.tend=opt.Opt(4)*opt.Opt(5);end
    opt.Opt(4)=RO.dt; 
  end
  if isfield(RO,'tend'); opt.Opt(5)=fix(RO.tend/opt.Opt(4));end
  out=opt; return;
  
 elseif comstr(Cam,'@');out=eval(CAM);return;
 else; 
  error('not a known type of time integration computation');
 end
 %% #Initialization_procedure ------------------------------------------------
 %% #Handle_TimeOpt - - - - - - - - - - - - - - - - - - - - - - - - - - - - -2
 opt.Opt=comstr(Cam,-1);
 if strncmpi(RunOpt.Out,'TimeOpt',7)&&nargin>1
  % [opt|mdl]=fe_time('TimeOpt[,Set,Reset] NLNewmark',[struct('Opt',[.25 .5 0 1e-5,99]),mdl],[,mdl])
  RB=[]; out=[]; RO=varargin{carg};carg=carg+1;
  if isfield(RO,'Node')||isfield(RO,'Elt'); % allow one arg as model=implicit set
   out=RO; RO=struct; if isempty(strfind(lower(RunOpt.Out),'set')); RunOpt.Out='TimeOptSet'; end
  elseif isfield(RO,'Opt');opt.Opt=RO.Opt; % xxx maybe other fields to handle
  elseif ~isstruct(RO)
   if ~isempty(RO);sdtw('_nb','bad opt input for TimeOpt call');end
   RO=struct;
  end
  if carg<=nargin&&isempty(out) % model as second argument
   out=varargin{carg}; carg=carg+1; RB=stack_get(out,'info','TimeOpt','get');
  end
  if isstruct(RB)&&~strncmpi(RunOpt.Out,'TimeOptReset',12);
   RO=sdth.sfield('addmissing',RO,RB); % keep ungiven fields from input model
  end
  RO=sdth.sfield('addmissing',opt,RO);
  opt=d_fetime(['TimeOpt' opt.Method],RO);
  if ~isempty(strfind(lower(RunOpt.Out),'set'))
   %strncmpi(RunOpt.Out,'TimeOptSet',10)||strncmpi(RunOpt.Out,'TimeOptReset',12)
   out=stack_set(out,'info','TimeOpt',opt);
  else; out=opt;
  end
  return
  
 elseif isempty(opt.Opt) && nargin>=4; opt.Opt=varargin{5};
 end
 if isfield(opt,'HHTalpha') % nlnewmark scheme with modified residual
  opt.HHTalpha=opt.Opt(1); % and num. damping
  opt.Opt(1:2)=[.25*(1+opt.Opt(1))^2 .5+opt.Opt(1)];
 elseif isfield(opt,'theta') % theta method
  if isempty(opt.Opt); opt.Opt(1)=opt.theta; end
  opt=rmfield(opt,'theta');
  opt.Opt=[opt.Opt(1) opt.Opt];
 end
 %if length(opt.Opt)<5 error('bad format for TimeOpt.Opt'); end

elseif isstruct(varargin{1})

 if isfield(varargin{1},'Method'); opt=varargin{1};carg=2;
 elseif sp_util('issdt'); 
   model=varargin{1};opt=fe_def('defTimeOpt',model);carg=1;
 else
   model=varargin{1};
   try; opt=stack_get(model,'info','TimeOpt','getdata');carg=1;
   catch; error('No method found');
   end
 end

else; error('First input argument must be a string or a data structure');
end

if carg>nargin; out=opt; return; end % returns TimeOpt

%% #get_model,_case_and_q0 - - - - - - - - - - - - - - - - - - - - - - - - -2
model=varargin{carg}; carg=carg+1;
q0=[];CaseName='Case 1';
if nargin>=carg
 Case = varargin{carg}; carg=carg+1;
 if isfield(Case,'def')&&~isfield(Case,'T')&&~isfield(Case,'Stack')
   q0=Case; [Case,CaseName]=fe_case(model,'getcase');
 elseif isempty(Case);      [Case,CaseName]=fe_case(model,'getcase');
 elseif ischar(Case);    [Case,CaseName]=fe_case(model,['getcase' Case]);
 else  
     try;CaseName=model.Stack{strcmpi(model.Stack(:,1),'case'),2};end
 end
else; [Case,CaseName]=fe_case(model,'getcase');
end
if isempty(CaseName);CaseName='Case 1';end
% q0=stack_get(model,'','q0','getdata'); % done later
% if isempty(q0) && carg<=nargin; 
%  q0=varargin{carg}; carg=carg+1; 
%  %warning('OpenFEM:q0','Taking initial condition given in input argument');
% end

if isequal(RunOpt.Out,'TimeOpt')
 model=stack_set(model,'info','TimeOpt',opt);
 out=model;return;
end

if ~isfield(opt,'Silent'); opt.Silent=0; end
% #build_time_vector - - - - - - - - - - - - - - - - - - - - - - - - - - - -2
t=[];
if ~isempty(strfind(lower(opt.Method),'static'))
 opt.nt=1; nt=1;
else
 if ~isfield(opt,'TimeVector')
  if isstruct(opt.Opt); t=0; nt=1;
  elseif ~isempty(opt.Opt);dt=opt.Opt(4); t=opt.Opt(3)+dt*[0:opt.Opt(5)]; nt=length(t);
  else;nt=1;
  end
 else; t=opt.TimeVector;  nt=length(t);
 end
 opt.nt=nt;
end

if ~isfield(opt,'MaxIter');opt.MaxIter=100;end

% #assemble_the_model_if_necessary - - - - - - - - - - - - - - - - - - - - -2
Load=[];
if sp_util('issdt')
 if isfield(opt,'AssembleCall');st=opt.AssembleCall;
 else;st='assemble -fetime';if isfield(opt,'NoT')&&opt.NoT;st=[st 'NoT'];end
 end
 if opt.Silent; st=[st ';']; end
 if ~isempty(strfind(lower(st),'load'))
   [model,Case,Load]=fe_case(stack_set(model,'case',CaseName,Case),st);
 elseif ~isempty(st);[model,Case,Load]=fe_case(stack_set(model,'case',CaseName,Case),st);
 end
else
 if isfield(model,'DOF') && isfield(model,'K')
 elseif ~isfield(opt,'NoT')||~opt.NoT;  model=fe_mknl(model); 
 else;model=fe_mknl(model,'NoT'); 
 end
 if length(model.K)<3&&isequal(model.Opt(2,:),[2 1])
   model.K={model.K{1}  ...
    spalloc(size(model.K{1},1),size(model.K{1},1),0) model.K{2}};
   model.Opt=[1 0 0 ;2 3 1];
 end
end

if isempty(Case.DOF); Case=fe_case(model,'gett',Case);end
mdof=Case.DOF;

% #assemble f(t)=fc.*ft - - - - - - - - - - - - - - - - - - - - - - - - - - -2
if isempty(Load);[fc,ft]=fe_load('buildu time',model,Case,t,opt);
else;[fc,ft]=fe_load('builduNoBuild time',model,Case,Load,t,opt);
end

if isempty(fc)
 if isfield(opt,'Residual'); disp('fc=residual')
 else; error('You must provide fc for linear solvers');
 end
elseif isfield(fc,'DOF');
 if isempty(fc.def); fc=spalloc(length(mdof),1,0);
    if ~isequal(size(ft),[opt.nt 1]);ft=zeros(opt.nt,1);end
 else; fc=fe_c(mdof,fc.DOF,fc.def')'; 
 end
end


u=zeros(length(mdof),1);v=zeros(length(mdof),1);
opt1=[];if isfield(opt,'Opt');opt1=opt.Opt;end;q1=[];

%% #check initial conditions - - - - - - - - - - - - - - - - - - - - - - - -2
q0=stack_get(model,'','q0','getdata');
if isempty(q0) && carg<=nargin; 
 q0=varargin{carg}; carg=carg+1; 
 %warning('OpenFEM:q0','Taking initial condition given in input argument');
end

if isempty(q0) % zeros initial conditions 
   u=zeros(length(mdof),1);v=u;
elseif isstruct(q0) && isfield(q0,'def') && isfield(q0,'DOF') 
   if size(q0.def,1)==length(q0.DOF); q1=q0.def;
   else
    q1(1:length(q0.DOF),1)=q0.def(:,1);
    if  isfield(q0,'v');  q1(1:length(q0.DOF),2)=q0.v; 
    else;q1(1:length(q0.DOF),2)=0.;
    end
   end
   % reorder CI dofs if necessary
   if (length(mdof)~=length(q0.DOF)  || norm(q0.DOF-mdof))
    q1=full(fe_c(mdof,q0.DOF,q0.def')');
  end
   u=q1(:,1); if size(q1,2)>=2; v=q1(:,2); else; v=zeros(size(u));end
elseif size(q0,1)==length(mdof)
   u=q0(:,1); if size(q0,2)>=2; v=q0(:,2); else; v=zeros(size(u));end
else;error('problem on initial conditions')
end
if size(q1,2)==3;a=q1(:,3);else; a=zeros(length(mdof),1);end

% #Give_standard_names_to_matrices - - - - - - - - - - - - - - - - - - - - -2
try;
 m=model.K{ismember(model.Opt(2,:),[2 20])};
 k=model.K{(model.Opt(2,:)==1|model.Opt(2,:)==5)};
 if any(model.Opt(2,:)==3); c=model.K{(model.Opt(2,:)==3)};
 else; c=spalloc(size(k,1),size(k,2),1);
 end
catch
 [m,c,k]=deal(model.K{1:3});
end
if isa(model,'v_handle'); model=model.GetData;end

% #initialize_the_output - - - - - - - - - - - - - - - - - - - - - - - - - -2
if ~isfield(opt,'NeedUVA'); opt.NeedUVA=[1 0 0]; end
% need refinenement : u only, it is imposed by a static case
if ~isempty(strfind(opt.Method,'static')); opt.NeedUVA=[1 0 0]; end

if (isfield(opt,'c_u') || isfield(opt,'c_v') || isfield(opt,'c_a'))  ... 
   && isfield(opt,'OutputFcn') && ~isa(opt.OutputFcn,'double')
 warning('Not compatible output options : observation matrix and OutputFcn')
end
if (isfield(opt,'c_u') || isfield(opt,'c_v') || isfield(opt,'c_a'))  ... 
   && isfield(opt,'OutInd')
 error('Not compatible output options : observation matrix and OutInd')
end

if ~exist('tout','var'); tout=[]; end
if isfield(opt,'OutputInit') && ~isempty(opt.OutputInit)
 %% This is the current default where things really occur
 eval(opt.OutputInit);
end%else
if (isfield(opt,'c_u') || isfield(opt,'c_v') || isfield(opt,'c_a')) && ...
    ~isfield(opt,'OutputFcn')
 %% A somewhat older variant that should dissapear
 out=struct('data',t(:),'fun',[0 4]);
 opt.OutputFcn='';
 if isfield(opt,'OutputFcn') && isa(opt.OutputFcn,'double')
  tout=opt.OutputFcn(:);tout=tout(tout>=min(t)&tout<=max(t));
 end
 if isfield(opt,'c_u')
  opt.OutputFcn=[opt.OutputFcn 'out.def(:,j1+1)=opt.c_u*u;']; out.def=[];
  if isfield(opt,'lab_u'); out.lab_u=opt.lab_u; end
 elseif opt.NeedUVA(1)==1; 
  opt.OutputFcn=[opt.OutputFcn 'out.def(:,j1+1)=u;'];
  out.def=[]; out.DOF=mdof;
 end
 if isfield(opt,'c_v')
  opt.OutputFcn=[opt.OutputFcn 'out.v(:,j1+1)=opt.c_v*v;']; out.v=[];
  if isfield(opt,'lab_v'); out.lab_v=opt.lab_v; end
 elseif opt.NeedUVA(2)==1; 
  opt.OutputFcn=[opt.OutputFcn 'out.v(:,j1+1)=v;'];
  out.v=[]; out.DOF=mdof;
 end
 if isfield(opt,'c_a')
  opt.OutputFcn=[opt.OutputFcn 'out.a(:,j1+1)=opt.c_a*a;']; out.a=[];
  if isfield(opt,'lab_a'); out.lab_a=opt.lab_a; end
 elseif opt.NeedUVA(3)==1; 
  opt.OutputFcn=[opt.OutputFcn 'out.a(:,j1+1)=a;'];
  out.a=[]; out.DOF=mdof;
 end

elseif ~isfield(opt,'OutputFcn') || isempty(opt.OutputFcn) || ...
    (ischar(opt.OutputFcn)&&~isfield(opt,'OutputInit')) || ...
    (ischar(opt.OutputFcn)&&isempty(opt.OutputInit)) 
   %% Obsolete calls 
   if isfield(opt,'OutInd'); opt.PostDim=[length(opt.OutInd) nt];
   else; opt.PostDim=[length(mdof) nt];
   end
   out=struct('DOF',mdof,'fun',[0 4]);
   if exist('t','var'); out.data=t(:);end
   if opt.NeedUVA(1)==1; out.def=[]; end
   if opt.NeedUVA(2)==1; out.v=[]; end
   if opt.NeedUVA(3)==1; out.a=[]; end

   if ~isfield(opt,'OutputFcn') || isempty(opt.OutputFcn)
    opt.OutputFcn='';
    if opt.NeedUVA(1)==1; opt.OutputFcn=[opt.OutputFcn 'out.def(:,j1+1)=u;']; end
    if opt.NeedUVA(2)==1; opt.OutputFcn=[opt.OutputFcn 'out.v(:,j1+1)=v;']; end
    if opt.NeedUVA(3)==1; opt.OutputFcn=[opt.OutputFcn 'out.a(:,j1+1)=a;']; end
   end
   if isfield(opt,'OutInd') 
    out.OutInd=int32(opt.OutInd);out.DOF=out.DOF(opt.OutInd); 
    out.cur=[0 t(1)];
    opt.OutputFcn='';
    if opt.NeedUVA(1)==1; 
     opt.OutputFcn=[opt.OutputFcn 'out.def(:,j1+1)=u(opt.OutInd);']; 
    end
    if opt.NeedUVA(2)==1; 
     opt.OutputFcn=[opt.OutputFcn 'out.v(:,j1+1)=v(opt.OutInd);']; 
    end
    if opt.NeedUVA(3)==1; 
     opt.OutputFcn=[opt.OutputFcn 'out.a(:,j1+1)=a(opt.OutInd);']; 
    end
   end
   
elseif isa(opt.OutputFcn,'double') % Vector of times is given
  tout=opt.OutputFcn(:);tout=tout(tout>=min(t)&tout<=max(t));
  if isfield(opt,'OutInd'); opt.PostDim=[length(opt.OutInd) length(tout)];
  else; opt.PostDim=[length(Case.DOF) length(tout)];
  end
  out=struct('DOF',Case.DOF,  ...
             'data',tout,    'cur',[0 t(1)],   ...
             'OutInd',int32(1:length(Case.DOF)), ...
             'fun',[0 4]);
  if opt.NeedUVA(1)==1; out.def=[]; end
  if opt.NeedUVA(2)==1; out.v=[]; end
  if opt.NeedUVA(3)==1; out.a=[]; end

  if isfield(opt,'OutInd') 
   out.OutInd=int32(opt.OutInd); out.DOF=out.DOF(opt.OutInd);
  end
end
if isfield(opt,'PostDim'); % allocate output
 uva=zeros(opt.PostDim);
 if opt.NeedUVA(1)==1; out.def=uva; end
 if opt.NeedUVA(2)==1; out.v=uva+0; end % +0 to duplicate data
 if opt.NeedUVA(3)==1; out.a=uva+0; end % +0 to duplicate data
 clear uva
end

cput0=cputime; % ki=[];
if sp_util('diag');fprintf('fe_time(''%s'') \n',opt.Method);end
if ~exist('ki','var'); ki=[]; end
% ----------------------------------------------------------------------------
switch lower(opt.Method)
%% #Newmark ------------------------------------------------------------------
case 'newmark' 

  opt=checkOutputFcn(opt,out,model); % out may be reassigned
  [opt,SaveTimes]=init_opt_linear(opt); 
 
  % beta = 0, gamma=.5 % Explicit
  % beta =.25, gamma=.5 % Classical values for implicit 

  dt0=0.; beta=opt1(1); gamma=opt1(2);
  % initial conditions
  j1=0;tc=t(1);dt=0.;
  Case.uva=initUVA(opt,out,model);
  eval(opt.OutputFcn);
  if ~isstruct(fc);fc=full(fc);end; 
  if isempty(opt.Residual); 
   opt.Residual='r=(ft(j1,:)*fc''-v''*c-u''*k)'';';
  end
  if ~isstruct(fc);fc=full(fc);end 
  
  if beta==0 %isequal([beta gamma],[0 0.5]) 
  %% #Newmark_Explicit -2
   % explicit (centered difference method shown in Geradin/Rixen
   % figure 7.4.1)
   % Possibility of using numerical damping using gamma>1/2
   % beta=(1+alpha^2)/4 and gamma=1/2+alpha
   %   amp error -alpha\omega^2h^2/2, freq error \omega^2 h^2/12
   % ISMA14_JointsDamping_paper.pdf z=-1/(2dt omega) ln(1-(g-1/2)omega^2dt^2) 

   v12=v+(t(2)-t(1))/2*a; 
   
   dt=t(2)-t(1); dt0=dt; ki=newmark_update(model,[],dt,dt0,opt1,opt);
   if max(abs(diff(diff(t))))<eps*100&&ki.ty(1)==6.1&& ...
           ~isempty(strfind(opt.Residual,'mkl'))
   %% exploit lumped mass gvdr test
    ki=ki.dinv(:); gh=gamma*dt; g1h=dt*(1-gamma); j1=0; ft=ft'; uva=Case.uva.uva;
    while (j1<length(t)-1); j1=j1+1; % exploit lumped mass gvdr test
     tc=t(j1+1);opt.tc=tc;dt=t(j1+1)-t(j1);
     ct=fc*ft(:,j1);
     of_time('storelaststep',model,Case,u,v,a,ct);
     u = u + dt*v12;
     eval(opt.Residual);
     a=ki.*r;
     v = v + gh*a + g1h*uva(:,3); %Case.uva.uva(:,3);
     v12 = v + (dt/2)*a;
     eval(opt.OutputFcn);dt0=dt;
     if ~isempty(SaveTimes) && (tc>=SaveTimes(1)||j1==length(t)-1);
      SaveTimes(1)=[];eval(opt.SaveFcn);
     end
     if ~isfinite(norm(u,'inf'))
      sdtw('_nb','The simulation has diverged (t(%i)=%g)',j1,tc); break
     end
    end % loop
   
   else 
   %% Non lumped mass loop on times - - - - - - - - - - - - - - - - -

   for j1=1:1:length(t)-1 
    tc=t(j1+1);dt=t(j1+1)-t(j1); opt.tc=tc;
    ct=ft(j1,:);if isstruct(fc); ct=fc;else;ct=(ct*fc')'; end
    of_time('storelaststep',model,Case,u,v,a,ct);
 
    % update iteration matrix if necessary
    if abs(dt0/dt-1)>1e-6 || isempty(ki);
      ki=newmark_update(model,ki,dt,dt0,opt1,opt);display(ki);
    end

    % increment displacement
    u = u + dt*v12;
    % increment velocity : assume v ~ v12
    % update acceleration
    eval(opt.Residual); % Residual must use v12
    a = ki\r; % mechanical resolution of new acceleration
    % correct v interpolation with new acc.
    v = v + dt*gamma*a + dt*(1-gamma)*Case.uva.uva(:,3); %v = v + dt/2*(a + Case.uva.uva(:,3)); 
    % increment half step velocity
    v12 = v + dt/2*a; % v12 = v12 + dt*a;
    
    % % was (jpb pass)
    % u=u+dt*v12; tc12=tc-dt/2; 
    % eval(opt.Residual); 
    % aa=a; v12 = v12 + dt*a; a = ki\r; v=v+dt*(a+aa)/2;
    
    eval(opt.OutputFcn);  dt0=dt;
    if ~isempty(SaveTimes) && (tc>=SaveTimes(1)||j1==length(t)-1);
     if ~isfinite(norm(u,inf)); 
      ofact('clear'); error('The simulation has diverged');
     end
     SaveTimes(1)=[];eval(opt.SaveFcn);  
    end

   end % loop on times - - - - - - - - - - - - - - - - -
   
   end
   
  else
  %% #implicit_case -2

    j1=0;
    while(j1<length(t)-1) % loop on times - - - - - - - - - - - - - - - - -
    j1=j1+1;

    tc=t(j1+1);dt=t(j1+1)-t(j1);opt.tc=tc;
    ct=ft(j1,:);if isstruct(fc); ct=fc;else;ct=(ct*fc')'; end
    % update iteration matrix if necessary
    if abs(dt0/dt-1)>1e-6 || isempty(ki);
      ki=newmark_update(model,ki,dt,dt0,opt1,opt);
    end

    of_time('storelaststep',Case.uva,u,v,a);
    % predictions
    u = u + dt*v + (.5-beta)*dt^2*a; v = v + (1-gamma)*dt*a;

    % update acceleration
    eval(opt.Residual);
    a = ki\r;
    % corrections
    u = u + beta*dt^2*a;  v = v + gamma*dt*a;
    eval(opt.OutputFcn);  dt0=dt;
    if ~isempty(SaveTimes) && (tc>=SaveTimes(1)||j1==length(t)-1);
     if ~isfinite(norm(u,inf)); 
       ofact('clear'); error('The simulation has diverged');
     end
     SaveTimes(1)=[];eval(opt.SaveFcn);  
    end

   end % loop on times - - - - - - - - - - - - - - - - -

  end % of explicit/implicit choice

  out=feutil('rmfield',out,{'cur','OutInd'});
  if cputime-cput0>30;
    fprintf('Elapsed CPU time : %i s\n',fix(cputime-cput0));
  end


% ----------------------------------------------------------------------------
%% #DG #Discrete #Galerkin ---------------------------------------------------
case 'dg' % 

 if ~isfield(opt,'RelTol'); opt.RelTol=sdtdef('OpenFEM.THRESHOLD-safe',1e-6);end
 
 opt=init_opt_linear(opt);
 if ~isstruct(fc);fc=full(fc);end ; u1=u;v1=v; u2=u;v2=v; a=v*0; 
 % update force
 f1=(ft(1,:)*fc')';f2=(ft(1,:)*fc')';

 % initial conditions
 j1=0; eval(opt.OutputFcn);
 
 ki=[];dt0=0.; if ~any(fc);error('No load');end

 for j1=1:1:length(t)-1
  tc=t(j1+1);dt=t(j1+1)-t(j1); 
  ki=dg_update(model,ki,dt,dt0);   
  % update forces
  f1=(ft(j1,:)*fc')';f2=(ft(j1,:)*fc')';
  %  update RHS
  f1s=dt/2*f1-dt/6*f2+5/3*m*v2-2/3*dt*k*u2;
  f2s=dt/2*f1+dt/2*f2+m*v2-dt*k*u2;
  % multicorrection
  j2=0;
  while (1)
   % update residual and correction
   v1s=ki\(f1s-2/3*m*v2);v2s=ki\(f2s-dt^2/3*model.K{3}*v1);
   r1=norm(f1s-ki*v1s-2/3*m*v2s)/norm(f1s);
   r2=norm(f2s-dt^2/3*model.K{3}*v1s-ki*v2s)/norm(f2s);
   v1=v1s; v2=v2s;
   if (r1<opt.RelTol && r2<opt.RelTol);  break; end
   if (r1>1e6 || r2>1e6);  error(1); ;end
   j2=j2+1;
  end %while
  fprintf('Step %i number of iterations =  %i\n',j1,j2);
  % update displacement
  u1=u1+dt/6*(v1-v2); u2=u2+dt/2*(v1+v2);
  u=u2; v=v2; 
  eval(opt.OutputFcn);
  dt0=dt;
 end %for j1

  
% ----------------------------------------------------------------------------
% #StaticNewton non linear statics -------------------------------------------
case 'staticnewton' % 

 dt=1;dt0=1;if isempty(ft); ft=1; end  % gradual load
 [opt,out]=init_static_newton_calls(opt,out,model);
 if ~isstruct(fc);fc=full(fc);end ;if ~exist('ki','var'); ki=[]; end
 opt.tc=0;
 opt.nt=size(ft,1);
 
 RunOpt=Follow_Update(RunOpt,opt); % First Init of follow if opt.Follow filled
 
 for j1=0:size(ft,1)-1
  opt.nf=0; 
  ct=ft(j1+1,:);if isstruct(fc); ct=fc;else;ct=(ct*fc')'; end % f_ext(n+1)
  [u,v,a,ki,opt] = feval(opt.IterFcn,ki,ct,u,0*u,[],dt,dt0,j1,model,opt,Case,j1);
  eval(opt.OutputFcn);
  if RunOpt.Follow
    RunOpt=Follow_Update(RunOpt,opt); % update follow and init at 1st call
  end
 end % j1   / ft
 if ischar(opt.SaveFcn);eval(opt.SaveFcn);
 else;feval(opt.SaveFcn);
 end

% ----------------------------------------------------------------------------
%% #NLNewmark : implicit + Newton --------------------------------------------
case 'nlnewmark' 

  opt=checkOutputFcn(opt,out,model); % out may be reassigned
  dt0=0.; beta=opt1(1); gamma=opt1(2); %   ki=[];
  % initial conditions
  j1=0; tc=t(1);dt=0.;opt.tc=tc;%IterFcn=@iterNewton;
  RunOpt=Follow_Update(RunOpt,opt); % First Init of follow if opt.Follow filled

  [opt,SaveTimes]=init_newton_calls(opt);
  % initial acceleration
  if ~isfield(opt,'InitAcceleration') || isempty(opt.InitAcceleration)
    sdtw('_nb','Using default acceleration init k\(f-k*u-c*v)');
    a=(ft(1,:)*fc')' - model.K{3}*u;
    if issparse(model.K{2});a=a-model.K{2}*v;end
    if ~isnumeric(m);r2=diag(m.GetData);else;r2=diag(m);end
    if size(m,1)<10;if ~isnumeric(m);a=full(m.GetData\a);else;a=full(m\a);end
    elseif nnz(r2)<size(m,1);a=zeros(size(a));
    else;a=ofact(m,full(a));
    end 
  else;                  eval(opt.InitAcceleration);
      if ~isempty(strfind(opt.InitAcceleration,';return;'));return;end
  end

  % Lookup figure 7.5.1 in Geradin/Rixen
  u=full(u);v=full(v);a=full(a);if ~isstruct(fc);fc=full(fc);end ; j1=0;
  Case.uva=initUVA(opt,out,model);
  eval(opt.OutputFcn); 
  if ~isfield(opt,'RelTol');opt.RelTol=sdtdef('OpenFEM.THRESHOLD-safe',1e-6);end
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  while j1<length(t)-1


   j1=j1+1;tc=t(j1+1);dt=tc-t(j1); 
   ct=ft(j1,:);if isstruct(fc); ct=fc;else;ct=(ct*fc')'; end
   of_time('storelaststep',model,Case,u,v,a,ct,tc);
   
   % predictions
   u = u + dt*v + (.5-opt1(1))*dt^2*a;
   v = v + (1-opt1(2))*dt*a;
   a = zeros(length(mdof),1);
   % - - - - - - - -  iterNewton Newton iterations => correction on u, v and a 
   
   [u,v,a,ki,opt] = feval(opt.IterFcn,ki,ct,u,v,a,dt,dt0,tc,model,opt,Case,j1);
   % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   eval(opt.OutputFcn);  dt0=dt; % typically fills out.def(:,j1+1), .v, .a

   if RunOpt.Follow
    RunOpt=Follow_Update(RunOpt,opt); % update follow and init at 1st call
   end

   if ~isempty(SaveTimes) && (tc>=SaveTimes(1)||j1==length(t)-1);
    if ~isfinite(norm(u,inf)); 
     ofact('clear'); error('The simulation has diverged'); 
    end
    SaveTimes(1)=[];
    if ischar(opt.SaveFcn);eval(opt.SaveFcn); 
    % fprintf('Time = %f          Newton iterations : %i \n',tc,opt.ite(1));
    else;feval(opt.SaveFcn);
    end
   end
  end % for j1

%% #ThetaMethod ---------------------------------------------------------------
case 'theta' % 
    
  opt=checkOutputFcn(opt,out,model); % out may be reassigned
  j1=0; tc=t(1);dt=0.;dt0=dt;
  if ~isfield(opt,'IterFcn')||isempty(opt.IterFcn); opt.IterFcn=@iterTheta; 
  elseif ~isa(opt.IterFcn,'function_handle'); opt.IterFcn=eval(sprintf('@%s',opt.IterFcn));
  end
  RunOpt=Follow_Update(RunOpt,opt); % First Init of follow if opt.Follow filled
  if ~isfield(opt,'IterInit'); opt.IterInit=''; end
  if ~isfield(opt,'JacobianUpdate'); opt.JacobianUpdate=0; end
  if ~isfield(opt,'PreCond'); opt.PreCond=0; end
  if ~isfield(opt,'diag'); opt.diag=0; end
  if ~isfield(opt,'Jacobian')%||isempty(strfind(opt.Jacobian,'theta'))
    opt.Jacobian='[ki,opt.kr]=theta_jacobian(model,ki,dt,dt0,opt.Opt);'; 
  end
  if ~isfield(opt,'SaveTimes'); SaveTimes=[]; else; SaveTimes=opt.SaveTimes; end
  if ~isfield(opt,'RHS')||isempty(opt.RHS)
   opt.RHS='Fn1=opt.kr*v+model.K{3}*(-dt*u)+(dt*opt.Opt(1))*fn1+(dt*(1-opt.Opt(1)))*fn;';
  end
  if ~isfield(opt,'nf'); opt.nf=0; end

  st='out.def(1:size(u,1),j1+1)=u;';
  try if ~isempty(strfind(opt.OutputFcn,'interp'));st='';end;end
  if ~isempty(st);opt.OutputFcn=st;end
  
  % Look PhD Thesis Xavier Lorang SNCF
  u=full(u);v=full(v);if ~isstruct(fc);fc=full(fc);end ;
  Case.uva=initUVA(opt,out,model);
  eval(opt.OutputFcn);
  
  % F_n+1=(M-dt^2 theta(1-theta)K-dt(1-theta)C) vn - dt K qn
  %       +dt * theta f_ext_n+1 + dt (1-theta) f_ext_n
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  while j1<length(t)-1

   j1=j1+1;tc=t(j1+1); dt=tc-t(j1);   
   % Compute force t_(n+1):fn1, and keep force t_(n): fn
   ct=ft(j1,:);if isstruct(fc); ct=fc;else;ct=(ct*fc')'; end % f_ext(n+1)
   if j1==1;  fn=ct; else;fn=fn1; end % f_ext(n)
   fn1=ct;  
   
   of_time('storelaststep',model,Case,u,v,a,ct,tc);
   if isempty(ki)||opt.JacobianUpdate; eval(opt.Jacobian); end
   % Constant RHS - - - - - - - - - - - - -
   eval(opt.RHS); %Fn1=opt.kr*v+model.K{3}*(-dt*u)+(dt*opt.Opt(1))*fn1+(dt*(1-opt.Opt(1)))*fn;
   v0=v; % keep previous time step vel for displacement updating
   % Iter - - - - - - - - - - - - - - - - -
   [u,v,a,ki,opt] = opt.IterFcn(ki,Fn1,u,v,a,dt,dt0,tc,model,opt,Case,j1);
   % Update final - - - - - - - - - - - - - 
   u=u+dt*(opt.Opt(1)*v+(1-opt.Opt(1))*v0); % Update displacement
   % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   eval(opt.OutputFcn);  dt0=dt; % typically fills out.def(:,j1+1), .v, .a
   
   if RunOpt.Follow
    RunOpt=Follow_Update(RunOpt,opt); % update follow and init at 1st call
   end
   
   if ~isempty(SaveTimes) && (tc>=SaveTimes(1)||j1==length(t)-1);
    if ~isfinite(norm(u,inf)); 
     ofact('clear'); error('The simulation has diverged'); 
    end
    SaveTimes(1)=[];eval(opt.SaveFcn); 
    fprintf('Time = %f            Theta-method iterations : %i \n',tc,opt.ite(1));
   end
  end % for j1
  
%% #euler : euleur solver for first order differential equations -----------
case 'euler' 

 opt=checkOutputFcn(opt,out,model); % out may be reassigned
 [opt,SaveTimes]=init_opt_linear(opt);
 
 dt0=0.; theta=opt1(1);
 % initial conditions
 j1=0;tc=t(1);dt=0.;
 Case.uva=initUVA(opt,out,model);
 eval(opt.OutputFcn);
 if ~isstruct(fc);fc=full(fc);end;
 if isempty(opt.Residual);
  if theta>0
   opt.Residual='r=(ft(j1,:)*fc''+(u+(1-opt.Opt(1))*dt*v)''*m/(opt.Opt(1)*dt))'';';
   %    elseif theta==1 % Backward Euler
   %     opt.Residual='r=(ft(j1,:)*fc''+u''*m/dt)'';';
  elseif theta==0 % Forward euler
   error('Forward Euler (explicit) is not available')
  end
 end
 if ~isfield(opt,'Jacobian')||isempty(opt.Jacobian)
  opt.Jacobian=...
   'ki=ofact(1/(opt.Opt(1)*dt)*model.K{1}.GetData+model.K{3}.GetData);';
  if ~isa(model.K{1},'v_handle');
   opt.Jacobian=strrep(opt.Jacobian,'.GetData','');
  end
 end
 
 if theta>0 % 1=backward euler
  j1=0;
  while(j1<length(t)-1)
   j1=j1+1;
   tc=t(j1+1);dt=t(j1+1)-t(j1);
   % update iteration matrix if necessary
   if abs(dt0/dt-1)>1e-6 || isempty(ki); eval(opt.Jacobian); end
   
   of_time('storelaststep',Case.uva,u,v,a);
   eval(opt.Residual);
   un=u; % last u
   u = ki\r;
   v=(u-un-(1-theta)*dt*v)/(theta*dt);
   eval(opt.OutputFcn);  dt0=dt;
   if ~isempty(SaveTimes) && (tc>=SaveTimes(1)||j1==length(t)-1);
    if ~isfinite(norm(u,inf));
     ofact('clear'); error('The simulation has diverged');
    end
    SaveTimes(1)=[];eval(opt.SaveFcn);
   end
   
  end % loop on times - - - - - - - - - - - - - - - - -
  
 else
  error('Not implemented');
 end
 
  out=feutil('rmfield',out,{'cur','OutInd'});
  if cputime-cput0>30;
    fprintf('Elapsed CPU time : %i s\n',fix(cputime-cput0));
  end
  
%% #Back returns model
case 'back' 
  out=model;out1=Case; out2=opt; 
  try; out3={Load,fc,ft,q0};end
  return;
  
case 'stabcheck';
%% #StabCheck : use of transition matrix to verify stability
% Newmark see PhD Vermot equation (4.52), Geradin Rixen (7.3.15)
1;
K=cellfun(@(x)x.GetData,model.K,'uni',0);
M=K{model.Opt(2,:)==2};C=K{model.Opt(2,:)==3};K=K{model.Opt(2,:)==1};
b=opt.Opt(1); g=opt.Opt(2); h=opt.Opt(4);
H1=[M+g*h*C g*h*K; b*h^2*C M+b*h^2*K];
H0=[(1-g)*h*C-M (1-g)*h*K; (.5-b)*h^2*C-h*M (.5-b)*h^2*K-M];
if size(H1,1)>500
 [U,d]=eigs(H1,H0,40,'smallestabs');d=diag(d);
else
 [U,d]=eig(full(H1),full(H0),'vector');
end
 
[d2,i1]=sort(abs(d)); 
figure(11);semilogy(1./d2);axis tight;
line(get(gca,'xlim'),[1 1],'linestyle',':');
i2=find(d2<1-sqrt(eps)); % unstable elements
def=struct('def',U(1:size(K,1),i2),'DOF',model.DOF,'data',d2(i2));
if isfield(model,'TR')
 if size(def.def,1)~=size(model.TR.def,2);return;end
 def.def=model.TR.def*def.def; def.DOF=model.TR.DOF;
 def.def(abs(def.def)<1e-8*norm(def.def,'Inf'))=0;
end
if isempty(def.def);fprintf('\nAll poles are stable\n')
else
 cf=feplot; if isempty(cf.mdl.Node); cf.model=model;end
 if 1==2
   cf=feplot;feplot(cf,'initmodel-back',model);fesuper('SEBuildSel',cf,{'s1','groupall'},def);   fecom('colordataEvalZ');
 else
  cf.def=def; %fecom(';colordataa;coloralpha');
 end
end
return;% set Breakpoint here for further analysis

%% #External method -----------------------------------------------------------
otherwise;
 %try; 
 if ~opt.Silent;  fprintf('Using ''%s'' method\n',opt.Method); end
%   if comstr(opt.Method,'static'); [CAM,Cam]=comstr(opt.Method,7); st='';
%   else;[CAM,i1,i2,i3]=sscanf(opt.Method,'%s',1);st=opt.Method(i3:end);
%   end
  if comstr(opt.Method,'static'); i1=7; else; i1=1; end
  r1=comstr(opt.Method,i1);
  [CAM,i1,i2,i3]=sscanf(r1,'%s',1);st=r1(i3:end);
  out=feval(CAM,['fe_time' st],model,Case,opt); 
 %catch;
 % error(sprintf('''%s'' not know',opt.Method)); 
 %end
end
if isa(ki,'ofact');ofact('clear',ki);end

out1=[];out2=[];
if exist('RunOpt','var')&&isfield(RunOpt,'Follow')
 cingui('TimerStop');
end
if isfield(opt,'FinalCleanupFcn')&&~isempty(opt.FinalCleanupFcn)
 try;eval(opt.FinalCleanupFcn);
 catch;fprintf('FinalCleanupFcn "%s" failed',opt.FinalCleanupFcn);
 end
else;% FinalCleanup stored as opt.Stack{'FinalCleanupFcn','name',r1}
  st=stack_get(opt,'FinalCleanupFcn');
  for j1=1:size(st,1);evt=st{j1,3}; feval(evt.Cb{:});end
end
if nargout>1;out1=model; out1.Case=Case; out2=opt; out3=Load; end

%% #SubFunc --------------------------------------------------------------------
% ------------------------------------------------------------------------------
%% #iterNewton : standard iterations for Newmark and Static newton -------------
function [u,v,a,ki,opt] = ...
    iterNewton(ki,fc,u,v,a,dt,dt0,tc,model,opt,Case,j1) %#ok<INUSL>

 % residual
 if ~isfield(opt,'nf'); opt.nf=0; end
 r=[];opt.tc=tc; eval(opt.IterInit);opt.ite=0;
 if opt.nf==0 % default init of residual
   eval(opt.Residual); r0=full(r);
   if opt.nf==0&&~isstruct(fc)
    if norm(full(fc)); opt.nf=norm(fc); 
    elseif any(u); opt.nf=max(1,norm(r0)); % with q0, r0 may be small
    else; opt.nf=norm(r0);
    end
   end
   if opt.nf<10*eps; opt.nf=1;end
 elseif isempty(r); eval(opt.Residual); % if IterInit has residual do it once
 end % default init of residual
 ite=0;opt1=opt.Opt;r2=0;r3=0;
 while(1)

  % update iteration matrix
  if opt.JacobianUpdate==0; if isempty(ki);eval(opt.Jacobian);end 
  elseif ite==0||opt.JacobianUpdate>0;  eval(opt.Jacobian); 
  end
  % compute corrections
  dq = ki\r;
  % corrections
  u = u - dq;
  if ~isempty(a) % newmark update (not needed for static)
   v = v - (opt1(2)/opt1(1))/dt*dq;
   a = a - 1/opt1(1)/dt^2*dq;
  end
  if opt.RelTol<0 % test convergence on displacement
    r2=norm(u);r3=norm(dq);
    if r2==0&&r3==0;break;
    elseif r3/r2<-opt.RelTol; break;
    elseif ~isfinite(r2);
        warning('Error : diverged');evalin('caller','j1=length(t)+1;'); 
        ofact('clear'); break;
    else;eval(opt.Residual);ite=ite+1;
    end
  else % test convergence on residual
   eval(opt.Residual);
   if opt.RelTol&&norm(r)/opt.nf<abs(opt.RelTol) ;break; end
   ite=ite+1;
  end 
  
  if ite>abs(opt.MaxIter(1));
   if opt.MaxIter(2) % StopTol at maxiter
    if (opt.MaxIter(2)<0&&(norm(dq)/norm(u))>-opt.MaxIter(2))||...
      (opt.MaxIter(2)>0&&norm(r)/opt.nf>opt.MaxIter(2))
     %ofact('clear',ki);
     sdtw('_nb','Max number of iterations reached outside stop tolerance tc %g tol %g',...
      tc,norm(r)/opt.nf); r2=norm(u);r3=norm(dq); break;
     %error('Max number of iterations reached outside stop tolerance');
    else; r2=norm(u);r3=norm(dq); break;
    end
   else
    fprintf('dq(%g)/q(%g)=%g\n',norm(dq),norm(u),norm(dq)/norm(u));
    fprintf('r(%g)/opt.nf(%g)=%g\n',norm(r),opt.nf,norm(r)/opt.nf);
    sdtw('_nb','Max number of iterations (%i) reached tc(%i)=%g',ite-1,j1,tc);
    if sp_util('diag')>1||any(~isfinite(r))||any(r>1e200)
     ofact('clear',ki);ki=[];evalin('caller','j1=length(t)+1;');
     warning('Stop');
    end;
    break;
   end
  elseif ite>abs(opt.MaxIter(1)) && opt.MaxIter(1)<0
   ofact('clear',ki);error('Max number of iterations reached');
  end
 end % while
 opt.ite=ite;eval(opt.IterEnd);

% ------------------------------------------------------------------------------
%% iterNewton_Sec --------------------------------------------------------------
function [u,v,a,ki,opt] = ...
    iterNewton_Sec(ki,fc,u,v,a,dt,dt0,tc,model,opt,Case,j1) %#ok<INUSL>

if ~isfield(opt,'RelTol');opt.RelTol=sdtdef('OpenFEM.THRESHOLD-safe',1e-6);end
 % residual
 if ~isfield(opt,'nf'); opt.nf=0; end
 r=[];eval(opt.IterInit);
 if opt.nf==0 % default init of residual
   eval(opt.Residual); r0=full(r);
   if opt.nf==0;
    if norm(full(fc)); opt.nf=norm(fc); else; opt.nf=norm(r0); end
   end
   if opt.nf==0; opt.nf=1;end
 elseif isempty(r); eval(opt.Residual); % if IterInit has residual do it once
 end % default init of residual
 ite=0;opt1=opt.Opt; r1=[0*u 0*u];
 while(1)

  % test convergence 
  if opt.RelTol<=0  % on displacement
    r2=norm(u);r3=norm(dq);
    if r2==0&&r3==0;break;
    elseif r3/r2<-opt.RelTol; break;
    elseif ~isfinite(r2);
        warning('Diverged');evalin('caller','j1=length(t)+1;');
        ofact('clear'); break;
    else;eval(opt.Residual);ite=ite+1;
    end
  else % on residual
      if norm(r)/opt.nf<opt.RelTol ;break; end
      eval(opt.Residual); ite=ite+1;
  end 
  % update iteration matrix
  if opt.JacobianUpdate==0; 
      if isempty(ki);eval(opt.Jacobian);end 
  elseif ite==0||opt.JacobianUpdate>0;  eval(opt.Jacobian); 
  end
  % compute corrections
  dq = ki\r;
  % Line Search for dq optimization
  rho=[0 1]; ite_sec=1;
  while(1)
   eval(opt.ResSec);
   r2=dq'*r1;
   r3=(diff(rho)/diff(r2))*r2(2);
   rho=[rho(2) rho(2)-r3];
   if abs(diff(rho))<eps; break
   elseif ite_sec>opt.MaxIterSec; break
   elseif isnan(r3)||~isfinite(r3)||r3>1e10
    if (diff(r2)<eps&&ite_sec>1)
     sdtw('_nb','Too much secant iterations ite %3i: cut-back @ ite_sec=%i',ite,ite_sec);
     rho(2)=rho(1);  opt.MaxIterSec=opt.MaxIterSec-1; break
    else
     rho(2)=1; break;
    end
   elseif abs(rho(2))<abs(opt.RelTol);
       sdtw('_nb','Control tends toward Zero, rho = %.2g @ ite %3i',[rho(2) ite]);
   end
   ite_sec=ite_sec+1;
  end
  dq=rho(2)*dq;
  % corrections
  u = u - dq;
  if ~isempty(a) % newmark update (not needed for static)
   v = v - (opt1(2)/opt1(1))/dt*dq;
   a = a - 1/opt1(1)/dt^2*dq;
  end
  % iteration control
  if ite>abs(opt.MaxIter(1)); 
   fprintf('dq(%g)/q(%g)=%g\n',norm(dq),norm(u),norm(dq)/norm(u));
   fprintf('r(%g)/opt.nf(%g)=%g\n',norm(r),opt.nf,norm(r)/opt.nf);
   sdtw('_nb','Max number of iterations (%i) reached tc(%i)=%g',ite-1,j1,tc);
   if sp_util('diag')>1||any(~isfinite(r));
       ofact('clear',ki);evalin('caller','j1=length(t)+1;');error('Stop');
   end;
   break; 
  elseif ite>abs(opt.MaxIter(1)) && opt.MaxIter(1)<0
   ofact('clear',ki);error('Max number of iterations reached');
  end 
 end % while
 opt.ite=ite-double(opt.RelTol>0);eval(opt.IterEnd);

% ------------------------------------------------------------------------------
%% #iterTheta, linear version --------------------------------------------------
 function [u,v,a,ki,opt] = iterTheta(ki,Fn1,u,v,a,dt,dt0,tc,model,opt,Case,j1); %#ok<INUSL>
  if isempty(ki);eval(opt.Jacobian);end % Defines ki and kr
  % of_time('storelaststep',model,Case,u,v,a,ct);
  % ct contains [F(t_n+1) F(t_n)]
  if opt.PreCond; v = 1/dt * (ki\(dt*Fn1));
  else; v = ki\Fn1;
  end
  a=[]; opt.ite=1; %a unused and should be left free

% ----------------------------------------------------------------------------
%% #init_newton_calls --------------------------------------------------------
function  [opt,SaveTimes]=init_newton_calls(opt);

  if ~isfield(opt,'Jacobian') ||  isempty(opt.Jacobian)
   opt.Jacobian='ki=basic_jacobian(model,ki,dt,dt0,opt.Opt);'; 
  end
  if ~isfield(opt,'Residual') || isempty(opt.Residual)
   if ~isfield(opt,'HHTalpha')
    opt.Residual='r = model.K{1}*a+model.K{2}*v+model.K{3}*u-fc;'; 
   else
    opt.Residual=horzcat('r=model.K{1}*a+',...
      'model.K{2}*(opt.HHTalpha*Case.uva.uva(:,2)+(1-opt.HHTalpha)*v)+',...
      'model.K{3}*(opt.HHTalpha*Case.uva.uva(:,1)+(1-opt.HHTalpha)*u)+',...
      'opt.HHTalpha*Case.uva.uva(:,4)-(1-opt.HHTalpha)*fc;');
   end
  end
  if ~isfield(opt,'IterFcn')
   opt.IterFcn=@iterNewton;
  else
   if comstr(opt.IterFcn,'iterNewton_Sec') % secant control method param.
      if ~isfield(opt,'ResSec')
       r1=textscan(opt.Residual,'%s',2,'delimiter','=;');r1=r1{1}{2};
       st={r1 r1};
       RunOpt.fields={'u' 'v' 'a';
    '(u-rho(%i)*dq)' '(v - (opt1(2)/opt1(1))/dt*rho(%i)*dq)' '(a - 1/opt1(1)/dt^2*rho(%i)*dq)'};
       RunOpt.fielddisp={',%s,','*%s','%s''*'};
       for j1=1:length(RunOpt.fields)
        for j2=1:length(RunOpt.fielddisp)
         for j3=1:2
           st{j3}=strrep(st{j3},sprintf(RunOpt.fielddisp{j2},RunOpt.fields{1,j1}),...
                 sprintf(RunOpt.fielddisp{j2},sprintf(RunOpt.fields{2,j1},j3)));
         end
        end % Attempt to guess the ResSec from the Residual
       end
       opt.ResSec=sprintf('r1(:,1)=%s;r1(:,2)=%s;',st{:});
      end
      if ~isfield(opt,'MaxIterSec'); opt.MaxIterSec=3; end
   end
   if ~isa(opt.IterFcn,'function_handle')
       opt.IterFcn=eval(sprintf('@%s',opt.IterFcn));
   end
  end
% r = (ft(j1,:)*fc'-v'*c-u'*k)';
  if ~isfield(opt,'JacobianUpdate'); opt.JacobianUpdate=0; end
  if isfield(opt,'SaveTimes') 
   pw0=evalin('base','pwd');  SaveTimes=opt.SaveTimes;
   if ~isfield(opt,'SaveFcn');
    opt.SaveFcn=sprintf('save %s out', ...
       fullfile(pw0,'fe_time_save.mat'));
   end
  else; opt.SaveTimes=[]; SaveTimes=[];
  end
  if ~isfield(opt,'IterInit'); opt.IterInit='';end
  if ~isfield(opt,'IterEnd');  opt.IterEnd='';end
  if length(opt.MaxIter)<2; opt.MaxIter(2)=0; end

% ------------------------------------------------------------------------------
%% #init_static_newton_calls ---------------------------------------------------
function  [opt,out]=init_static_newton_calls(opt,out,model);

  eval(iigui({'fc','ft','u','v','a','Case'},'GetInCaller')); t=0;
  opt=checkOutputFcn(opt,out,model); 
  if ~isfield(opt,'Jacobian') ||  isempty(opt.Jacobian)
    opt.Jacobian='ki=basic_jacobian(model,ki,0.,0.,opt.Opt);'; 
  end
  if ~isfield(opt,'Residual') || isempty(opt.Residual)
    opt.Residual='r = model.K{3}*u-fc;'; 
  end
  if ~isfield(opt,'IterFcn') || isempty(opt.IterFcn)
    opt.IterFcn=@iterNewton;
  else
   if comstr(opt.IterFcn,'iterNewton_Sec') % secant control method param.
    if ~isfield(opt,'ResSec')
        r1=textscan(opt.Residual,'%s',2,'delimiter','=;');r1=r1{1}{2};
        st={r1 r1};
        RunOpt.fields={'u' 'v' 'a';
     '(u-rho(%i)*dq)' '(v - (opt1(2)/opt1(1))/dt*rho(%i)*dq)' '(a - 1/opt1(1)/dt^2*rho(%i)*dq)'};
        RunOpt.fielddisp={',%s,','*%s','%s''*'};
        for j1=1:length(RunOpt.fields)
         for j2=1:length(RunOpt.fielddisp)
          for j3=1:2
            st{j3}=strrep(st{j3},sprintf(RunOpt.fielddisp{j2},RunOpt.fields{1,j1}),...
                  sprintf(RunOpt.fielddisp{j2},sprintf(RunOpt.fields{2,j1},j3)));
          end
         end % Attempt to guess the ResSec from the Residual
        end
        opt.ResSec=sprintf('r1(:,1)=%s;r1(:,2)=%s;',st{:});
%      opt.ResSec=...
%   'r1(:,1)=model.K{3}*(u-rho(1)*dq)-fc; r1(:,2)=model.K{3}*(u-rho(2)*dq)-fc;';
    end
    if ~isfield(opt,'MaxIterSec'); opt.MaxIterSec=3; end    
    opt.IterFcn=eval(sprintf('@%s',opt.IterFcn));
   end
  end
  if ~isfield(opt,'SaveFcn'); opt.SaveTimes=[]; opt.SaveFcn=''; end
  % r = fc'--u'*k)';
  if ~isfield(opt,'JacobianUpdate'); opt.JacobianUpdate=0; end
  if ~isfield(opt,'IterInit'); opt.IterInit='';end
  if ~isfield(opt,'IterEnd');  opt.IterEnd='';end
  if ~isfield(opt,'RelTol'); opt.RelTol=sdtdef('OpenFEM.THRESHOLD-safe',1e-6);end
  if length(opt.MaxIter)<2; opt.MaxIter(2)=0; end
  opt.Opt=[1. 0.];% to use basic_jacobian

  % Clean up out
  if size(out.data,1)~=size(ft,1); 
      out.data=ft; 
      if isfield(out,'def');out.def(1,end+1:size(ft,1))=0;end % Prealloc
      out.fun(2)=1; % Static
  end

% ------------------------------------------------------------------------------
%% #basic_jacobian -------------------------------------------------------------
function  ki=basic_jacobian(model,ki,dt,dt0,opt); %#ok<DEFNU>

if isfield(opt,'oProp');oProp=opt.oProp;
else; oProp=stack_get(model,'','oProp','g');
end
if isempty(oProp);oProp={};end

if dt==0 || abs(dt0/dt-1)>1e-10 || (isnumeric(ki) &&isempty(ki))
  if (isa(ki,'ofact') && ~isempty(ki)); ofact('clear',ki); end
  if dt==0; dt=1; end
  ki=model.K{3}; if ~isnumeric(ki);ki=ki.GetData;end
  if length(model.K{1})<10
      ki=full( ...
        ki+(opt(2)/opt(1))/dt*model.K{2} + 1/opt(1)/dt^2*model.K{1});
  elseif issparse(model.K{2})||isa(model.K{2},'v_handle')||isa(model.K{2},'vhandle.matrix')
      ki=ki+(opt(2)/opt(1))/dt*model.K{2} + 1/opt(1)/dt^2*model.K{1};
      if nnz(ki)==size(ki,1);ki=ofact(ki,'diag');else;ki=ofact(ki,oProp{:});end
  elseif isfield(model.K{2},'K')
      ki=ki+(opt(2)/opt(1))/dt*model.K{2}.K + 1/opt(1)/dt^2*model.K{1};
      if nnz(ki)==size(ki,1);ki=ofact(ki,'diag');else;ki=ofact(ki,oProp{:});end
  else
      ki=ki + 1/opt(1)/dt^2*model.K{1}; % auto detect diagonal
      if nnz(ki)==size(ki,1);ki=ofact(ki,'diag');else;ki=ofact(ki,oProp{:});end
  end
end

% ------------------------------------------------------------------------------
% #theta_jacobian - - ----------------------------------------------------------
function  [ki,kr]=theta_jacobian(model,ki,dt,dt0,opt); %#ok<DEFNU>
% J=M+dt*theta*C+dt^2*theta^2*K
% Kr=(M-dt^2 theta(1-theta)K-dt(1-theta)C)
if isfield(opt,'oProp');oProp=opt.oProp;
else; oProp=stack_get(model,'','oProp','g');
end
if isempty(oProp);oProp={};end

if dt==0 || abs(dt0/dt-1)>1e-10 || (isnumeric(ki) &&isempty(ki))
  if (isnumeric(ki) && ~isempty(ki)); ofact('clear',ki); end
  if dt==0; dt=1; end
  ki=model.K{1}; if ~isnumeric(ki);ki=ki.GetData;end
  if issparse(model.K{2})||isa(model.K{2},'v_handle')||isa(model.K{2},'vhandle.matrix')
     kr=ki;
     kr=kr-dt*(1-opt(1))*model.K{2} - (dt^2*opt(1)*(1-opt(1)))*model.K{3};
     ki=ki+dt*opt(1)    *model.K{2} + (dt*opt(1))^2*model.K{3};
     if nnz(ki)==size(ki,1);ki=ofact(ki,'diag');else;ki=ofact(ki,oProp{:});end
  else;error('Not implemented');
  end
end

% ------------------------------------------------------------------------------
%% #init_opt_linear - - --------------------------------------------------------
function [opt,SaveTimes]=init_opt_linear(opt)

 if ~isfield(opt,'Residual') || isempty(opt.Residual)
    opt.Residual='';% 'r=(fc*ft(j1)-model.K{2}*v-model.K{3}*u);';
 end
 if isfield(opt,'SaveTimes') 
  pw0=evalin('base','pwd');  SaveTimes=opt.SaveTimes;
  if ~isfield(opt,'SaveFcn');
   opt.SaveFcn=sprintf('save %s out', ...
      fullfile(pw0,'fe_time_save.mat'));
  end
 else; opt.SaveTimes=[]; SaveTimes=[];
 end

% ------------------------------------------------------------------------------
%% #newmark_update - - ---------------------------------------------------------
function  ki=newmark_update(model,ki,dt,dt0,Opt,opt);   
% Jacobian formulation of Newmark for implicit linear, or explicit (acc formulated)
if isfield(opt,'oProp');oProp=opt.oProp;
else; oProp=stack_get(model,'','oProp','g');
end
if isempty(oProp);oProp={};end
if isfield(opt,'Jacobian')&&~isempty(opt.Jacobian)
 Case=evalin('caller','Case');
 u=evalin('caller','u');v=evalin('caller','v');% for nl_spring jacobian call
 eval(opt.Jacobian);
elseif norm((dt-dt0)/dt)>1e-10 || (isnumeric(ki) &&isempty(ki))
 ofact('clear',ki);
 kiter=model.K{1}; if ~isnumeric(kiter);kiter=kiter.GetData;end
 if Opt(1)==0 % beta=0, explicit : rethrow Mass as integrator : nothing to do
 else
  % rethrow usual implicit Newmark Jacobian formulated in acceleration :
  % ki = M + gamma h C + beta h^2 K
  c=model.K{2};if ~isnumeric(c); c=c.GetData; end
  if isa(model.K{2},'struct') % Damping is specific (modal, nl_modaldmp)
   kiter=kiter+Opt(1)*dt^2*model.K{3}; % beta = Opt(1) (stiffness term)
   if isfield(model.K{2},'K'); kiter=kiter+Opt(2)*dt*model.K{2}.K; end % gamma =Opt(2)
  else; kiter=kiter+Opt(2)*dt*c+Opt(1)*dt^2*model.K{3}; % usual Opt(1)=beta, Opt(2)=gamma
  end
 end
 % factorization
 if nnz(kiter)==size(kiter,1)&&isequal(kiter,diag(diag(kiter)))||...
   max(max(abs(kiter-diag(diag(kiter)))))/max(abs(diag(kiter)))<100*eps;
  ki=ofact(kiter,'diag');
 elseif size(kiter,1)<10;ki=full(kiter);
  if any(diag(ki)==0);error('Zero terms on the diagonal of ki');end
 else; ki=ofact(kiter,oProp{:});
 end
end

% ------------------------------------------------------------------------------
%% #dg_update - - --------------------------------------------------------------
function  kiter=dg_update(model,kiter,dt,dt0);   

if norm((dt-dt0)/dt)>1e-10 || (isnumeric(kiter) &&isempty(kiter))
 kiter=model.K{1}; if ~isnumeric(kiter);kiter=kiter.GetData;end
 kiter=(kiter+(dt^2/6)*model.K{3});
end

% ------------------------------------------------------------------------------
%% #checkOutputFcn - - ---------------------------------------------------------
function  opt=checkOutputFcn(opt,out,model);

if isfield(model,'FNL')&&~isempty(model.FNL); st=',model.FNL';else;st='';end
%toutc=0;
if ~isempty(which('nl_spring')) % Bypass for NL vibration toolbox
   u=evalin('caller','u');t=evalin('caller','t');Case=evalin('caller','Case');
   ft=[]; 
   if evalin('caller','exist(''ft'',''var'')');ft=evalin('caller','ft');end
   evalin('caller','clear out');
   feval('nl_spring','OutputInitCheck');
   if isfield(opt,'cv')&&isa(opt.cv,'vhandle.chandle')%Do not reinit
   elseif isfield(opt,'chandle')&&isequal(opt.chandle,'model.cv')
    % Initial attempts at using chandle in fe_time (lem20 tests)
    opt.cv=nl_spring('nlchandlemodel',model);opt.FNL=model.FNL;%copy to model.FNL
    opt.Residual='r=-full(fc);mkl_utils(''residual'',r,opt.cv,u,v,a,opt,Case);';
   end

   assignin('caller','out',out);assignin('caller','opt',opt);
   %assignin('caller','toutc',toutc); % xxx when is this needed
   %return not needed
elseif ~isfield(opt,'OutputFcn') || isempty(opt.OutputFcn)
 if isfield(opt,'OutInd')
  if sp_util('issdt')
   opt.OutputFcn= [...
    'of_time(''interp'',out,beta,gamma,Case.uva,a,tc-dt,tc ' st ')'] ;
  else % openFEM only supports Newmark for this functionality
   opt.OutputFcn= [...
    'of_time(''newmarkinterp'',out,beta,gamma,Case.uva.uva,a,tc-dt,tc ' st ')'] ;
  end
 end
elseif isa(opt.OutputFcn,'double') % Vector of times is given
 % xxx suppress t=evalin('caller','t'); toutc=t(1);assignin('caller','toutc',toutc);
 if sp_util('issdt')
  opt.OutputFcn= [...
   'of_time(''interp'',out,beta,gamma,Case.uva,a,tc-dt,tc ' st ')'] ;
 else % openFEM only supports Newmark for this functionality
  opt.OutputFcn= [...
   'of_time(''interp'',out,beta,gamma,Case.uva.uva,a,tc-dt,tc ' st ')'] ;
 end
end
i3=0;
if ~isempty(strfind(lower(opt.OutputFcn),'interp'))&& ...
        (~isfield(out,'cur')||~isfield(out,'OutInd'))
    sdtw('_nb','Missing out.cur and .OutInd, adding standard fields');
    out.OutInd=int32(1:length(out.DOF))'; out.cur=[0 out.data(1)];
    if ~isempty(st) % model.FNL is actually used
     out.DOF=[out.DOF; model.FNLDOF];
     out.def=zeros(length(out.DOF),length(out.data));
    end
end
if (isfield(opt,'c_u')||isfield(opt,'c_v')) ...
        &&~isempty(strfind(opt.OutputFcn,'interp'))
 warning('SDT:CuCv_Implement','c_u c_v support not implemented in of_time interp')
end
if isfield(out,'def')&&~isfield(out,'Xlab');out.Xlab={'DOF','step'};end

% ------------------------------------------------------------------------------
%% #initUVA - - ----------------------------------------------------------------
function uva=initUVA(opt,out,model); %#ok<INUSL,*INUSD>

  if evalin('caller','exist(''u'',''var'')'); u=evalin('caller','u');end
  if evalin('caller','exist(''v'',''var'')'); v=evalin('caller','v');end
  if evalin('caller','exist(''a'',''var'')'); a=evalin('caller','a');end
  try; uva=[u v a 0*a];end
  if isfield(opt,'HHTalpha')||any(strfind(lower(opt.Method),'explicit'));   
   fc=[]; ft=[];
   if evalin('caller','exist(''fc'',''var'')'); fc=evalin('caller','fc');end
   if evalin('caller','exist(''ft'',''var'')'); ft=evalin('caller','ft');end
   if isstruct(fc); fc=[]; % xxx does this ever happen ?
   elseif isempty(ft); fc=sum(fc,2);
   else; fc=(ft(1,:)*fc')';
   end
   if ~isempty(fc); uva(:,4)=-full(fc); end %[uva -full(fc)]; end
  end
  if isa(uva,'double');uva=full(uva);end
  uva=struct('uva',uva,'FNL',[],'tc',zeros(1,3));
  % in case of non linearities, allow saving for interpolation
  if isfield(model,'FNL'); uva.FNL=model.FNL*0; end 
  %if isfield(opt,'ite');uva.ite=opt.ite*0; end % save iteration info
  
% ------------------------------------------------------------------------------
%% #Follow_Update - - ----------------------------------------------------------------
function RunOpt=Follow_Update(RunOpt,opt); %#ok<*INUSD>
  
 % Init - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 if ~isfield(RunOpt,'Follow') % first call= init Runopt
  if isfield(opt,'Follow'); 
   RunOpt.Follow=1;
   RunOpt.iFollow=1; % indice of the row to be filled in opt.Follow (needed by blocksave strategy)
  else;RunOpt.Follow=0;
  end
  return % Needed because start must occure after 1st iteration (opt.ite needed)
 end
 
 % Start - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 if RunOpt.iFollow==1 % first init of follow tracker
  RunOpt.jFollow=int32(1:length(opt.ite));
  t=evalin('caller','t'); ft=evalin('caller','ft'); 
  if isnumeric(opt.Follow); opt.Follow='cingui(''TimerStartPlot'')'; end
  if ischar(opt.Follow) % e.g. 'cingui(''TimerStartPlot'')'
   if evalin('caller','isfield(Case,''uva'')')
    uva=evalin('caller','Case.uva'); % needed for uva.tc init
   else;uva=[]; % static new ton case for example...
   end
   eval(opt.Follow);     
   assignin('caller','opt',opt);
  end
 end
 
 % Update - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 try
  sp_util('setinput',opt.Follow,opt.ite,int32(RunOpt.iFollow),RunOpt.jFollow);
  RunOpt.iFollow=RunOpt.iFollow+1; % indice of the row to fill in opt.Follow
 catch err
  sdtw('_nb','Follow tracker crashed with following error:')
  err.getReport
  RunOpt.Follow=0; % disable follow
 end
 
