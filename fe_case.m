function [out,out1,out2,out3]=fe_case(varargin) %#ok<STOUT>

%FE_CASE FEM computation case construction.
%
%	Syntax : Case = fe_case(Case,'PropertyName','Description',PropData)
%                fe_case(model,'command' ...)
%
%       FEM computation cases describe boundary conditions, constraints,
%       loads, and parametric design points. FE_CASE is used to initialize
%       (when Case is not provided as first argument) or modify cases 
%       (Case is provided).
%
%       This function scans trough CASES (see sdtweb('case')) defined 
%       in model.Stack and builds left hand side arguments needed for 
%       response evaluation.
%       
%	In each CASE, information needed to build the load is given in the
%        .Stack, the result is defined using DOFs declared in CASE.DOF
%
%       Stack entries for you give 'PropertyName','Label',PropData are
%        for boundary conditions
%         'KeepDof', 'FixDof', 'MPC'
%        for loads
%         'DOFLoad', 'DOFSet', 'Fvol', 'FSurf'
%        for physical parameters
%         'par'
%
%	Accepted commands are
%         Assemble[...]  calls used to assemble the matrices of a model
%            see on-line documentation with sdtweb('fe_case')
%         GetCase        returns the current case fe_case(model,'getcase')
%         GetDataName    returns data for case entry Name
%         GetT           returns a congruent transformation matrix 
%                        which verifies constraints
%         GetSensDof     retruns unique SensDof from model case
%         Reset          removes all entries in the current case
%         Remove         model=fe_case(model,'remove','name') removes entry
%         SetCurve       fe_case(model,'SetCurve','EntryName',curve) sets
%                        the .curve field of a particular entry
%
%   See sdtweb   sdt/case, fe_case
%	See also     fe_mk

%       Etienne Balmes
%       Copyright (c) 1996-2024 by SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*NOSEM>
if nargin==1 && comstr(varargin{1},'cvs')
 out='$Revision: 1.162 $  $Date: 2024/11/21 18:14:12 $'; return;
end

if nargin==0&&nargout==1
 out=struct('Stack',[],'T',[],'DOF',[]); out.Stack={};
 return;
                             % ('Command',model, ... )
elseif nargin>1 && ischar(varargin{1}) && isfield(varargin{2},'Elt')
 [CAM,Cam]=comstr(varargin{1},1);
 Case=varargin{2}; carg=3;
elseif nargin>1 && ischar(varargin{2}) && isfield(varargin{1},'Stack')
 [CAM,Cam]=comstr(varargin{2},1);
 Case=varargin{1}; carg=3;
elseif isa(varargin{1},'sdth');Case=varargin{1}.mdl;carg=2;
  if carg<=nargin; [CAM,Cam]=comstr(varargin{carg},1);carg=carg+1;
  else;CAM=''; Cam='';end
elseif ~isstruct(varargin{1})&&~isjava(varargin{1}) % ('Command' ...
 Case=fe_case;Case.Stack={}; carg=1;CAM=''; Cam='';
else                          % (model,'Command')
  Case=varargin{1};carg=2; 
  if carg<=nargin; [CAM,Cam]=comstr(varargin{carg},1);carg=carg+1;
   if isfield(Case,'parent') %isempty(CAM)&&carg<=nargin&&ischar(varargin{carg}) % cbk
    model=clean_get_uf('feplotcf',Case.parent); 
    if ~isempty(model); Case=model.mdl; else; Case=[]; end
    [CAM,Cam]=comstr(varargin{carg},1);carg=carg+1;
   end
  else; CAM=''; Cam='';
  end
end
if isfield(Case,'Node')||isfield(Case,'Elt')||~isfield(Case,'T')
     model=Case; if carg<=nargin; r1=varargin{carg}; else;r1=[];end
     if comstr(Cam,'getcase')
      [out,out1,CAM,Cam,model]=get_case(CAM,Cam,model);
      return;
     elseif comstr(Cam,'getdata') % single entry stack_get
      [Case,CaseName,st,st,model]=get_case('getcase','getcase',model);
      [CAM,Cam]=comstr(CAM,8);
      st='multi';if comstr(Cam,'-hdf');st='hdf';[CAM,Cam]=comstr(CAM,5);end
      if isempty(CAM); CAM=varargin{carg};carg=carg+1;end
      if iscell(CAM);if length(CAM)<2; error('Must provided typ,name');end
           out=stack_get(Case,CAM{:},st);
      else;out=stack_get(Case,'',CAM,st);
      end
      return;
 %% #Stack
     elseif comstr(Cam,'stack_get') % #stack_get multi entry
      [Case,CaseName,st,st,model]=get_case('getcase','getcase',model);
      out=stack_get(Case,varargin{carg:end});return;
     elseif comstr(Cam,'stack_set') % #stack_set multi entry
      if isequal(Cam(end),';'); sil=1; else; sil=0; end
      [Case,CaseName,st,st,model]=get_case('getcase','getcase',model);
      li=varargin(carg:end);
      if iscell(li)&&isequal(size(li),[1 1]); li=li{1}; end
      if size(li,1)>1; li=reshape(li',1,[]);  end
      li(1:3:end)=cleanUpperCType(li(1:3:end),sil);
      if ~isempty(Case.Stack);Case.Stack(:,1)=cleanUpperCType(Case.Stack(:,1),sil);end
      if strcmp(CaseName,'model.Case')
       if isempty(strfind(Cam,'new')); Case=stack_set(Case,li{:});
       else; Case=fe_def('stacknew',Case,li{:});
       end
       model.Case=Case; out=model;
      else
       if isempty(strfind(Cam,'new'))
        out=stack_set(model,'case',CaseName,stack_set(Case,li{:}));
       else; out=stack_set(model,'case',CaseName,fe_def('stacknew',Case,li{:}));
       end
      end
      return;
     elseif comstr(Cam,'stack_rm') % #stack_rm
      [Case,CaseName,st,st,model]=get_case('getcase','getcase',model);
      [Case,CaseRm]=stack_rm(Case,varargin{carg:end});
      out=stack_set(model,'case',CaseName,Case);
      out1=CaseRm;return;
     elseif comstr(Cam,'setdata')
      [Case,CaseName,st,st,model]=get_case('getcase','getcase',model);
      [CAM,Cam]=comstr(CAM,8);
      Case=stack_set(Case,'',CAM,varargin{carg});carg=carg+1;
      out=stack_set(model,'case',CaseName,Case);
      return;
     elseif isfield(r1,'Stack')&&isfield(r1,'T')&&isfield(r1,'DOF')
      Case=r1; carg=carg+1; CaseName='';
     elseif isfield(model,'Case')&&isfield(model.Case,'Stack')
       Case=model.Case;CaseName='model.Case';
     else
      [Case,CaseName,st,st,model]=get_case('','',model);
     end

else 
  model=[]; CaseName=''; 
  if comstr(Cam,'getcase');   out=Case; out1='Case 1'; return; end
end
%if carg==0 carg=1;end
if isempty(Case); Case=fe_case; end % default empty case
if isempty(Cam)&&carg==1&&ischar(varargin{1}) 
 [CAM,Cam]=comstr(varargin{1},1);carg=2;
end

epsl=sdtdef('epsl');RO=struct;

while ~isempty(Cam) % loop on arguments
 

%% #Add - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if comstr(Cam,'add'); [CAM,Cam]=comstr(CAM,4);

  [Case,CaseName,CAM,Cam,model]=get_case(CAM,Cam,model);

%% #Assemble [m,k,load] - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'assemble'); [CAM,Cam]=comstr(CAM,9);

if isempty(model)&&carg<=nargin
 model=varargin{carg}; carg=carg+1; 
 [Case,CaseName,CAM,Cam,model]=get_case(CAM,Cam,model);
end

if sp_util('issdt')
  eval(sprintf('[out,out1,out2,out3]=fe_caseg(''assemble%s'',model,Case,CaseName,nargout);',CAM));
else
  error('Not implemented yet');
end % OpenFEM

return;

% Case (return a particular case)  - - - - - - - - - - - - - - - -
elseif comstr(Cam,'case')

if length(Cam)>4
     [Case,CaseName,CAM,Cam,model]=get_case(CAM,Cam,model);
end
out=Case; return

elseif comstr(Cam,'merge')
%% #Merge : Merge common type of data (mpc, fixdof,...)- - - - - - - - - - - -
 st=varargin{carg}; carg=carg+1; % List of names or regexp (begin with #)
 st2=varargin{carg}; carg=carg+1; % Name of the returning mpc    
 [Case,r1]=stack_rm(Case,'',st,'get');
 typ=unique(r1(:,1));
 if length(typ)>1; error('Can only merge constraints of same type'); end
 typ=typ{1};
 if strcmpi(typ,'mpc')
  %% Merge.mpc
  r2=r1{1,3};
  for j1=2:length(r1) % Loop to merge each mpc
   r2=fe_mpc('mpcmerge',r2,r1{j1,3});
  end
 elseif strcmpi(typ,'fixdof')
  %% Merge.fixdof
  r2=[];
  for j1=1:length(r1) % Loop to merge each mpc
   r2=[r2;getFixDof(model,Case,r1(j1,:))];
  end
 else; error('Merging "%s" constraints not handled yet',typ);
 end

 Case=stack_set(Case,typ,st2,r2);

   %% Direct call with two mpcs%% #Get ----------------------------------------------------------------------
elseif comstr(Cam,'get'); [CAM,Cam]=comstr(CAM,4);


if comstr(Cam,'c')
%% #GetCase (get case given input) - - - - - - - - - - - - - - - - - -

 error('You should never get here');

elseif comstr(Cam,'t'); [CAM,Cam]=comstr(CAM,2);
%% #GetT (get congruent transformation matrix) - - - - - - - - - - - -

  if carg<=nargin; Case=varargin{carg};carg=carg+1; end

  if ~isfield(model,'DOF')||isempty(model.DOF)
    if isempty(model)&&isfield(Case,'DOF')&&~isempty(Case.DOF)
      Case.T=speye(length(Case.DOF)); out=Case; return;
    elseif ~isfield(model,'Elt')
     error('Element definitions are needed for case  T matrix building');
    end
    model.DOF=feutil('getdof',model);
  end
  if 1==2%nargout==1; % should really never happen but used in fe_mpc to force getnode 
      if isfield(Case,'cGL');Case=rmfield(Case,'cGL');end
  else % normally the GetNode is done here and not in fe_mpc
    [Case.Node,Case.bas,NNode,Case.cGL]=feutil(['getnodebasUsed-force' CAM],model);
    out1=NNode;
  end
  out2=model.DOF;
  % early return if transformation is defined
  if isfield(Case,'T')&&length(model.DOF)==size(Case.T,1)&& ...
     isfield(Case,'DOF')&&length(Case.DOF)==size(Case.T,2)&& ...
     ~comstr(Cam,'new')
       out=Case;
       if comstr(Cam,'dof'); out=Case.DOF;end
       return;
  end

  %% deal with possible MPC/rigid connections

  DOF=model.DOF; T=speye(length(DOF));
  if carg<=nargin % dof to keep given as argument
   i1=varargin{carg};carg=carg+1;
   if ~isempty(i1); DOF=DOF(i1); T=T(:,i1); end
  end

  RunOpt.fixdof=[];
  for j1=1:size(Case.Stack,1) % loop on things to account for
    switch lower(Case.Stack{j1,1})
    case 'fixdof'
         if comstr(Cam,'nogett'); continue;end
         r1=Case.Stack(j1,:); 
         r1=getFixDof(model,Case,r1);
         RunOpt.fixdof=[RunOpt.fixdof;r1(:)];
     end
  end

  if ~isfield(model,'Node')
  elseif comstr(Cam,'nogett')
  else
   [Case,pdof]=fe_mpc(model,Case,DOF,RunOpt);
   T=Case.T;DOF=Case.DOF;
   if ~isempty(strfind(Cam,'-mdof')); Case.mDOF=pdof; end
  end

  % Deal with DofSet entries (This has not been checked for multiple)
  if ~isempty(Case.Stack); 
   for j1=find(strcmpi('dofset',Case.Stack(:,1)))'
         r1=Case.Stack{j1,3};
         r1=safeDofSet(r1,model,Case);
         i1=fe_c(DOF,r1.DOF,'ind'); 
         if isscalar(r1.DOF)&&fix(r1.DOF)==0;r2=[];
         else;r2=fe_c(r1.DOF,DOF,'dof',2);
         end
         
         if ~isempty(r2)
          % DofSet may have been fixed so put them back in
           i2=fe_c(Case.DOF,r2,'ind');
           if ~isempty(i2)
             DOF=fe_c(Case.DOF,[DOF;r1.DOF],'dof');
             i1=fe_c(DOF,r1.DOF,'ind');
           end
          r2=fe_c(r1.DOF,DOF,'dof',2);
          if ~isempty(r2)&&~strcmpi(Case.Stack{j1,2},'thermalstate')
           fprintf( ...
           'The following DOFSet are not in the model DOF list\n%s\n', ...
           sdtw('_clip 50 1','%s ',fe_c(r2)));
          end
         end % check elimination of fixed DOFs
         if isempty(i1) && ...
           isempty(intersect(unique(round(rem(model.DOF,1)*100)), ...
             unique(round(rem(r1.DOF,1)*100))))
          if ~strcmpi(Case.Stack{j1,2},'thermalstate')
           sdtw('_nb','DofSet %s entry affects no free DOF',Case.Stack{j1,2})
          end
         else
          if isfield(Case,'TIn')&&~isempty(Case.TIn);
           error('Multiple DofSet entries are not supported');
          end
          if ~isempty(i1);Case.TIn=T(:,i1); Case.DofIn=DOF(i1); end
          if ~isfield(r1,'def')||~isfield(Case,'TIn');
          elseif isfield(r1,'data')&&isfield(r1,'def')&& ...
                  size(r1.def,2)==size(r1.data,1) % Enforce history provided
          elseif size(r1.def,1)==size(Case.TIn,2)
              Case.TIn=Case.TIn*sparse(r1.def);
          elseif size(r1.def,1)>size(Case.TIn,2)
              error('DofSet affects non present DOFs');
          else
              r2=feutilb('placeindof',Case.DofIn,r1);
              Case.TIn=Case.TIn*sparse(r2.def);
          end
          if isfield(r1,'adof');Case.AdofIn=r1.adof;end
          if isfield(r1,'KeepDof')&&r1.KeepDof==1 
            % Do not eliminate DofSet (used for time integration)
            % xxxT(:,fe_c(DOF,r1.DOF,'ind'))=0; %#ok<SPRIX> % No Shape in T but keep
          else
           i1=fe_c(DOF,r1.DOF,'ind');if ~isempty(i1);DOF(i1)=[];T(:,i1)=[];end
          end
         end  
   end
  end

  for j1=1:size(Case.Stack,1) % loop on things to account for
    switch lower(Case.Stack{j1,1})
    case 'fixdof'
         i1=fe_c(DOF,RunOpt.fixdof,'ind',2); DOF=DOF(i1); T=T(:,i1);
    case {'dofset','keepdof','par','rigid','dofload','fsurf', ...
          'fvol','sensdof', ...
          'sensstrain','cyclic','rbe3','mpc','nastran','info'}
    case 'pred'
    %% #pred partial reduction (direct implementation rather than trough MPC) -3
    d1=Case.Stack{j1,3}; 
    i1=fe_c(DOF,d1.DOF,'ind');i2=fe_c(d1.DOF,DOF(i1),'ind');
    i3=setdiff(1:size(DOF,1),i1);
    T=[T(:,i3) T(:,i1)*d1.def(i2,:)];
    DOF=[DOF(i3);d1.adof];

    case 'pcond'
      %% Pcond.do -3
      RunOpt.pcond=Case.Stack{j1,3}; % Save Preconditioner callback
    otherwise; sdtw('_nb','%s not supported by fe_case GetT',Case.Stack{j1,1});
    end % of supported case
  end
  if isfield(RunOpt,'pcond');eval(RunOpt.pcond);end
  if ~isempty(strfind(Cam,'mdof'));Case.mDOF=model.DOF;end
  Case.T=T; Case.DOF=DOF; out=Case;
  if comstr(Cam,'dof'); out=Case.DOF;end
  return;

elseif comstr(Cam,'sensdof'); [CAM,Cam]=comstr(CAM,8);
 %% #GetSensDof : if several SensDof, look for Test else take first one---2
 % [wire,name]=fe_case(cf.mdl,'getSensDof'); 
 r2=fe_case(model,'stack_get','SensDof');
 if isempty(r2); r2=[]; st='';
 elseif size(r2,1)==1; st=r2{1,2};r2=r2{1,3};
 else % Several SensDof, take 'Test' or the first one
  i1=strcmpi(r2(:,2),'Test');
  if any(i1); st=r2{i1,2};r2=r2{i1,3};
  else % Several Sens
   st=r2{1,2};r2=r2{1,3};
   sdtw('_nb','Several SensDof found, first one taken : %s',st)
  end
 end
 out=r2; out1=st; return;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
else;error('''get%s'' not a known command',CAM);
end

%% #Set
%% #Cyclic ------------------------------------------------------------------2
elseif comstr(Cam,'cyclic'); [CAM,Cam]=comstr(CAM,7);

   if comstr(Cam,'default')
    name='new cyclic';adof=[];
   else
    if carg>nargin; error('You must specify some data'); end
    name=varargin{carg}; carg=carg+1;
    if ~ischar(name); error('You must specify a name of  FixDof case entry');end
    data=varargin{carg}; carg=carg+1;
   end
   Case=stack_set(Case,'cyclic',name,data); Case.T=[];

%% #UN=0 --------------------------------------------------------------------2
elseif comstr(Cam,'un=0'); [CAM,Cam]=comstr(CAM,7);

    name=varargin{carg}; carg=carg+1;
    data=varargin{carg}; carg=carg+1;
    
    r1=data.normal; r1=r1./(max(abs(r1),[],2)*[1 1 1]);
    ind=find(any(abs(r1+1)<sqrt(eps),2));r1(ind,:)=-r1(ind,:);
   DOF=[data.ID(:)+.01;data.ID(:)+.02;data.ID(:)+.03];
   data=struct('DOF',DOF,'c',[diag(r1(:,1)) diag(r1(:,2)) diag(r1(:,3))]);
   Case=stack_set(Case,'mpc',name,data); Case.T=[];


%% #FixDOF ------------------------------------------------------------------2
elseif comstr(Cam,'fixdof'); [CAM,Cam]=comstr(CAM,7);

   if comstr(Cam,'default')
    name='new FixDof';adof=[];
   else
    if carg>nargin; error('You must specify some data'); end
    name=varargin{carg}; carg=carg+1;
    if ~ischar(name); error('You must specify a name of  FixDof case entry');end
   adof=varargin{carg}; carg=carg+1;
   end
   if isnumeric(adof); adof=adof(:);end
   if ~isempty(Case.Stack) 
    i1=find(strcmpi('FixDof',Case.Stack(:,1)));
    if ~isempty(i1); Case.Stack(i1,1)={'FixDof'};end
   end
   if ischar(adof);
         i1=strfind(comstr(adof,-27),'-id'); 
         if ~isempty(i1);
             [ID,i2,i3,i4]=sscanf(adof(i1+3:end),'%i');
             adof(i1:i1+1+i4)=''; adof=comstr(adof,1);
             adof=struct('data',adof,'ID',ID);
         end
   end
   Case=stack_set(Case,'FixDof',name,adof); Case.T=[];

%% #AutoSPC -----------------------------------------------------------------2
elseif comstr(Cam,'autospc') 

 eval('Case=fe_caseg(CAM,model,Case);'); % SDT extension

%% #KeepDof -----------------------------------------------------------------2
elseif comstr(Cam,'keepdof'); [CAM,Cam]=comstr(CAM,8);

   if comstr(Cam,'default')
    name='new KeepDof';adof=[];
   else
    if carg>nargin; error('You must specify some data'); end
    name=varargin{carg}; carg=carg+1;
    if ~ischar(name); error('You must specify a name of a KeepDof case entry');end
    adof=varargin{carg}; carg=carg+1;
   end
   if ~isempty(Case.Stack) 
    i1=find(strcmpi('keepdof',Case.Stack(:,1)));
    if ~isempty(i1); Case.Stack(i1,1)={'KeepDof'};end
   end

   Case=stack_set(Case,'KeepDof',name,adof(:));Case.T=[];

%% #sensDOF SDT extension ->sdtweb fe_caseg('sensdof') ----------------------------2
%% #sens SDT extension ->sdtweb fe_caseg('sens') ----------------------------2
elseif comstr(Cam,'sens')

 eval('[Case,carg]=fe_caseg(CAM,model,Case,varargin(carg:end),carg);')
 if ischar(carg);eval(carg);return; end
 

%% #MPC ---------------------------------------------------------------------2
elseif comstr(Cam,'mpc'); [CAM,Cam]=comstr(CAM,4);
 
   if comstr(Cam,'default')
    name='new MPC';adof=[];
   else
    if carg>nargin; error('You must specify some data'); end
    name=varargin{carg}; carg=carg+1;
   end
    
   if ~ischar(name); error('You must specify a name of MPC case entry');end
   adof=varargin{carg}; carg=carg+1;
   if isa(adof,'cell')&&size(adof,2)==3
    adof=struct('SID',adof{1},'DOF',adof{2},'c',adof{3});
   elseif isfield(adof,'c')&&isfield(adof,'DOF')
   elseif isfield(adof,'type')&&isfield(adof,'sel')
   else;error('Not a proper MPC format');
   end
   Case=stack_set(Case,'mpc',name,adof);Case.T=[];
 %% #PRED ---------------------------------------------------------------------2
elseif comstr(Cam,'pred'); [CAM,Cam]=comstr(CAM,5);
 
   if carg>nargin; error('You must specify some data'); end
   name=varargin{carg}; carg=carg+1; 
   if ~ischar(name); error('You must specify a name of MPC case entry');end
   adof=varargin{carg}; carg=carg+1;
   if ~isfield(adof,'def')||~isfield(adof,'DOF')||~isfield(adof,'adof')
     error('Expecting .def, .DOF, .adof fields')
   end
   Case=stack_set(Case,'pred',name,adof);Case.T=[];

%% #Rigid -------------------------------------------------------------------2
elseif comstr(Cam,'rigid'); [CAM,Cam]=comstr(CAM,6);

 if comstr(Cam,'default'); name='Rigid'; data=[]; % #RigidDefault, generate empty stack entry -3
 else % additional data provided
  if carg>nargin; error('You must specify a name of a Rigid case entry'); end
  name=varargin{carg}; carg=carg+1;
  if ~ischar(name); error('You must specify a name of a Rigid case entry');end
  if carg>nargin; error('You must specify some data'); end
  r1=varargin{carg}; carg=carg+1;
  if comstr(Cam,'append'); % #RigidAppend, use preexisting data -3
   data=stack_get(Case,'rigid',name,'getdata');
  else; data=[];
  end
  if isfield(r1,'Elt'); data=r1; r1=[];% allow large def support
  elseif r1(1)==Inf % list of elts is provided with header
  else % list of the form : [MasterNode slaveDOF slaveNode_1 slaveNode_1...]
   r1=r1(:)';
   r1=[ones(size(r1,2)-2,1)*r1(1) r1(3:end)' ones(size(r1,2)-2,1)*r1(2)];
   r1=feutil('addelt',[],'rigid',r1); % add header in this case only
  end
  % append to data and stack
  if ~isempty(r1)
   if isfield(data,'Elt'); data=data.Elt;end % robust to .Sel (to be reset anyway)
   data=feutil('addelt',data,r1);
  end
 end
 Case=stack_set(Case,'rigid',name,data);Case.T=[];

%% #DOFLoad -----------------------------------------------------------------2
elseif comstr(Cam,'dofload'); [CAM,Cam]=comstr(CAM,8);
 
   if comstr(Cam,'default')
    name='new dofload';adof=struct('name',name,'DOF',[],'def',[]);
   else
    if carg>nargin; 
        error('fe_case(model,''DofLoad'',''name'',data), missing name'); 
    end
    name=varargin{carg}; carg=carg+1;
    if ~ischar(name); error('You must specify a name of DOFLoad case entry');end
    if carg<=nargin; adof=varargin{carg}; carg=carg+1;
    else; adof=[];carg=carg+1; % Needed to allow empty adof
    end
   end

   if isfield(adof,'DOF') && isfield(adof,'def')
    adof.name=name;
    if size(adof.DOF,1)~=size(adof.def,1)
       error('.DOF and .def fields not consistent with DOF load case');
    end
    Case=stack_set(Case,'DOFLoad',name,adof);
   elseif isfield(adof,'DOF') && isfield(adof,'MasterNodes')
    % Rbe3 based generation of load, allows moments on volumes
    data=adof;
    data.Node=unique(fix(data.DOF));
    mo1=feutil('rmfield',model,'Stack');
    r2=data.MasterNodes(:,1);r2=r2(:)';r2=r2([1 1 1],:);
    r2(1,:)=1; % Weight
    r2(2,:)=123; % MasterDOF
    r2=[1 data.Node 123456 r2(:)'];
    mo1.Elt=feutil('addelt','mass1',[data.Node;data.MasterNodes(:)]*ones(1,7));
    mo1=fe_case(mo1,'rbe3','slave',r2);
    [C1,NNode,mo1.DOF]=fe_case(mo1,'gett');
    adof=struct('def',C1.T(fe_c(mo1.DOF,data.DOF,'ind'),:)','DOF',C1.DOF, ...
      'vertex',mo1.Node(NNode(data.Node),5:7));
    if isfield(data,'name');adof.name=data.name;end
    Case=stack_set(Case,'DOFLoad',name,adof);
    % Verify consistency of computed work
    % rb=feutilb('geomrb',mo1,adof.vertex,adof.DOF);adof.def'*rb.def

   elseif isfield(adof,'DOFi') && isfield(adof,'cbush')
    % resolve cbush observation matrix 
    % given DOF are internal directions
    % see ABAQUS connector load
    elt=feutil(sprintf('selelt %s',adof.cbush),model);
    if isa(model,'v_handle'); mo1=model.GetData; else; mo1=model; end
    mo1.Elt=elt; mo1.Node=feutil('getnodegroupall',mo1);
    mo1=fe_case(mo1,'reset');
    C1=fe_mknl('initnocon',mo1);
    cta = cbush(mo1.Node,mo1.Elt(2,:),[],mo1.il,[-2 1 1],C1);
    cta.def=cta.cta'; cta.def=cta.def(:,adof.DOFi)*diag(adof.def);
    cta.info=sdth.sfield('addmissing',adof.info,rmfield(adof,'info'));
    cta.info.bas=cta.bas;
    cta=feutil('rmfield',cta,'cta','bas');
    Case=stack_set(Case,'DOFLoad',name,cta);

   elseif ischar(adof)&&strncmpi(adof,'rb',2) % rb{selection,origin}
     %  'Case{DofLoad,Base,"inelt{proid111&selface&facing .9 0 0 -1000}"}'
    r1=sdth.findobj('_sub:',adof);
    st2=r1(2).subs; if ~iscell(st2);st2={st2};end
    adof=struct('type','rigid','sel',st2{1},'ori',[0 0 0]);
    i1=sdtm.regContains(st2,'^dir');
    if any(i1); adof.dir=comstr(comstr(st2{i1},4),-1);end
    i1=sdtm.regContains(st2,'^curve');
    if any(i1); adof.curve={comstr(st2{i1},6)};end
    if length(r1)>2;
        adof.ori=r1(3).subs{1};
    end
    i1=sdtm.regContains(st2,'^KeepDof','i');if any(i1);adof.KeepDof=1;end
    Case=stack_set(Case,'DOFLoad',name,adof);
   elseif isempty(Cam)&&isnumeric(adof)&&size(adof,2)==1
    sdtw('_nb','Assuming unit DofLoad entries');
    r1=struct('DOF',adof,'def',eye(size(adof,1)), ...
     'name',name,'lab',{fe_curve('doflab load',adof)});
    Case=stack_set(Case,'DOFLoad',name,r1);
   else
     try
      st=['[Case,carg]=fe_caseg([''DofLoad''' ...
            ' CAM],model,Case,varargin(carg-2:end),carg-2);'];
      eval(st); if ischar(carg);eval(carg);return; end
     catch;
      if sp_util('issdt'); eval(st);end
      error('For DOFLoad data.DOF data.def fields must be defined')
     end
   end

elseif comstr(Cam,'dofset'); [CAM,Cam]=comstr(CAM,7);
%% #DOFSet ------------------------------------------------------------------2
   if ~isempty(Cam)&&Cam(end)==';'; sil=1; Cam(end)=''; else; sil=0; end
   if comstr(Cam,'default')
    name='new DOFSet';adof=struct('name',name,'DOF',[],'def',[]);
   else
    if carg>nargin; error('You must specify some data'); end
    name=varargin{carg}; carg=carg+1;
    if ~ischar(name); error('You must specify a name of DOFSet case entry');end
    adof=varargin{carg}; carg=carg+1;
   end

   if isfield(adof,'DOF') && isfield(adof,'def')
    adof.name=name;
    if size(adof.DOF,1)~=size(adof.def,1)&&~ischar(adof.def)
       error('.DOF and .def fields not consistent with DOF load case');
    end
    Case=stack_set(Case,'DOFSet',name,adof);
   elseif isfield(adof,'DOF') && isfield(adof,'dir')
    adof.name=name;
    Case=stack_set(Case,'DOFSet',name,adof);
   elseif isnumeric(adof)&&size(adof,2)==1
    if ~sil; sdtw('_nb','Assuming unit DofSet entries'); end
    if any(rem(adof,1)==0); % eliminate wild cards
      if isfield(model,'DOF')&&~isempty(model.DOF);DOF=model.DOF;
      else; DOF=feutil('getdof',model);
      end
      adof=fe_c(DOF,adof,'dof');
    end
    r1=struct('DOF',adof,'def',eye(size(adof,1)),'name',name);
    Case=stack_set(Case,'DOFSet',name,r1);
   elseif ischar(adof)&&strncmpi(adof,'rb',2) % rb{selection,origin}
     %  'Case{FixDof,Base,"inelt{proid111&selface&facing .9 0 0 -1000}"}'
    r1=sdth.findobj('_sub:',adof);
    st2=r1(2).subs; if ~iscell(st2);st2={st2};end
    adof=struct('type','rigid','sel',st2{1},'ori',[0 0 0]);
    i1=sdtm.regContains(st2,'^dir');
    if any(i1); adof.dir=comstr(comstr(st2{i1},4),-1);end
    i1=sdtm.regContains(st2,'^curve');
    if any(i1); adof.curve={comstr(st2{i1},6)};end
    if length(r1)>2;
        adof.ori=r1(3).subs{1};
    end
    i1=sdtm.regContains(st2,'^KeepDof','i');if any(i1);adof.KeepDof=1;end
    Case=stack_set(Case,'DOFSet',name,adof);
   else
     error('For DOFSet data.DOF data.def fields must be defined')
   end

%% #Fvol --------------------------------------------------------------------2
elseif comstr(Cam,'fvol'); [CAM,Cam]=comstr(CAM,5);
 
   if comstr(Cam,'default')
    name='new FVol';
    r1=struct('name',name,'sel','GroupAll','dir',{{'0','0','0'}});
   else
    if carg>nargin; error('You must specify some data'); end
    name=varargin{carg}; carg=carg+1;
    if ~ischar(name); error('You must specify a name of Fvol case entry');end
    r1=varargin{carg}; carg=carg+1;
    r1.name=name;
   end

   if ~isfield(r1,'sel')
    error('You must specify .sel (element selection) for volume force');
   end
   % either a def/DOF or dir array
   if isfield(r1,'def')&&isfield(r1,'DOF');
   elseif ~isfield(r1,'dir'); r1.dir={{'1','0','0'}};
   end

   Case=stack_set(Case,'FVol',name,r1);

%% #FSurf -------------------------------------------------------------------2
elseif comstr(Cam,'fsurf'); [CAM,Cam]=comstr(CAM,6);

   if comstr(Cam,'default')
    name='new FSurf';
    r1=struct('name',name,'presel','','sel','','def',[1;1;1],'DOF',(1:3)/100);
   else
    if carg>nargin; error('You must specify some data'); end
    name=varargin{carg}; carg=carg+1;
    if ~ischar(name); error('You must specify a name of FSurf case entry');end
    r1=varargin{carg}; carg=carg+1;
    r1.name=name;
   end

   if ~isfield(r1,'sel')&&~isfield(r1,'set')&&~isfield(r1,'Elt')
    error('You must specify a face selection for surface force definition');
   end
   if ~isfield(r1,'DOF'); r1.DOF=[];end
   Case=stack_set(Case,'FSurf',name,r1);

%% #Par (SDT parameters) ----------------------------------------------------2
elseif comstr(Cam,'par'); [CAM,Cam]=comstr(CAM,4);

   if comstr(Cam,'stack'); [CAM,Cam]=comstr(CAM,6);end%ParStack compat
   if comstr(Cam,'add'); [CAM,Cam]=comstr(CAM,4);end

   if sdtm.Contains(Cam,'default')
    if comstr(Cam,'m'); name='new mass parameter';
         r1=struct('sel','groupall','coef',[2  1 .5 2 1]);
    else; name='new stiffness parameter';
         r1=struct('sel','groupall','coef',[1  1 .5 2 1]);
    end
   elseif comstr(Cam,'coef');
      des=stack_get(Case,'par');
      coef=varargin{carg};carg=carg+1;
      if size(coef,1)==1&&size(coef,2)==size(des,1);coef=coef(:);end
      if ~any(size(coef,2)==[1 5]);error('Expecting coef 1 or 5 column');end
      for j1=1:size(coef,1);
         r2=des{j1,3}.coef; 
         if size(coef,2)==1; r2(1,2)=coef(j1);else;r2=coef(j1,:); end
         des{j1,3}.coef=r2;
      end
      Case=stack_set(Case,des);r1=[];
   elseif comstr(Cam,'reset');Case=stack_rm(Case,'par');r1=[];
   else
     name=varargin{carg};carg=carg+1;
     if isstruct(varargin{carg})
      r1=varargin{carg};carg=carg+1;
      if ~isfield(r1,'coef'); r1.coef=[1 1 .5 2 1];end
      r1.name=name;
     else;r1=struct('sel','groupall','coef',[1 1 .5 2 1],'name',name);
     end
     
     if comstr(Cam,'m');r1.coef(1)=2;[CAM,Cam]=comstr(CAM,2);
     elseif comstr(Cam,'kg');r1.coef(1)=5;[CAM,Cam]=comstr(CAM,3);
     elseif comstr(Cam,'k');r1.coef(1)=1;[CAM,Cam]=comstr(CAM,2);
     elseif comstr(Cam,'t');r1.coef(1)=3;[CAM,Cam]=comstr(CAM,2);
     elseif comstr(Cam,'cut');r1.coef(1)=100;[CAM,Cam]=comstr(CAM,2);
     elseif comstr(Cam,'c');r1.coef(1)=3.1;[CAM,Cam]=comstr(CAM,2);
     elseif comstr(Cam,'ik'); r1.coef(1)=4; [CAM,Cam]=comstr(CAM,3);
     elseif comstr(Cam,'0'); r1.coef(1)=0; r1.sel='EltInd0'; [CAM,Cam]=comstr(CAM,2);
     end
     [RO,st,CAM]=cingui('paramedit -DoClean',[ ...
      'nom(#%g#"Nominal parameter") '...
      'matid(#%g#"parameter type as matid") '...
      'min(#%g#"Minimum") max(#%g#"Maximum")  scale(#%s#"Type of scale") '...
      'zCoef(#%s#"value for zCoefFcn") '...
      ],{struct,CAM});Cam=lower(CAM);
     if ~isempty(RO.matid);r1.coef(1)=RO.matid;end
     if ~isfield(r1,'sel'); r1.sel='groupall'; end
     
     if ~isempty(strfind(RO.scale,'lo'));
         r1.coef(5)=2; 
         if ~isempty(RO.min);RO.min=log10(RO.min);end
         if ~isempty(RO.max);RO.max=log10(RO.max);end
         if ~isempty(RO.nom);RO.nom=log10(RO.nom);end
     end    
     
     if ~isempty(RO.min);r1.coef(3)=RO.min;end
     if ~isempty(RO.max);r1.coef(4)=RO.max;end
     if ~isempty(RO.nom);r1.coef(2)=RO.nom;end
     % Assuming order nom, min, max, scale
     r1.coef(2:end)=comstr(CAM,[-1 r1.coef(2:end)]);
     if ~isempty(RO.zCoef);r1.zCoef=RO.zCoef; end % allow zcoef init
     if carg<=nargin
       [CAM,Cam]=comstr(varargin{carg},1);carg=carg+1;
       if ~isempty(CAM); r1.sel=CAM;end
     end
   end
   if ~isempty(r1);Case=stack_set(Case,'par',name,r1);end
   
elseif comstr(Cam,'setcurve') % #SetCurve -2
 % unified setcurve strategy, couples the old setcurve, fe_curve set LoadCurve
 % and new features (04/04/2014)
 if nargin<4;
  error('fe_case curve needs at least the case stack entry names and the curve names')
 end
 name=varargin{carg}; carg=carg+1;
 st=varargin{carg}; carg=carg+1;
 if strcmpi(name,'?')
  i1=ismember(lower(Case.Stack(:,1)),{'dofload','dofset','fvol','fsurf'});
  if nnz(i1)>1; error('? only works for a single load');
  else; i1=find(i1); 
    if nnz(i1)==0&&isfield(model,'Load');r1=model.Load;RO.out='model.Load';
    else; r1=Case.Stack{i1,3};end
  end
 else
  [r1,i1]=stack_get(Case,'',name,'Get'); % fetch concerned case entry
 end
 if length(i1)>1
  error('Several case types exist with the name %s, cannot SetCurve',name)
 end
 if isempty(r1);  % check EntryName
  error('Provided case entry %s does not exist in model',name)
 elseif ~isstruct(r1) % check that entry name supports curve
  error('Case entry %s, of type %s does not support curve entries',...
   name,Case.Stack{i1,1});
 end
 % assign curves
 if ischar(st)&&comstr(lower(st),'remove')
  r1=feutil('rmfield',r1,'curve');
 else
  if ~iscell(st); st={st}; end
  if carg<=nargin&&~isnumeric(varargin{carg}) % curve is directly provided
   R1=varargin{carg}; carg=carg+1;
   if ~iscell(R1); R1={R1}; end  % add curves to model stack
   if length(R1)~=length(st);
    error('number of curve names (%i) does not match number of input curves (%i)',...
     length(st),length(R1));
   end
   for j1=1:length(R1); 
      if ischar(R1{j1});% Force transfor to struct
        model=fe_curve(model,'set',st{j1},R1{j1});
      else;
       model=stack_set(model,'curve',st{j1},R1{j1}); 
      end
   end
  elseif all(cellfun(@isstruct,st)) % curves are directly provided
  else  % some curves are already stacked in model, check that they exist
   if ~isfield(model,'Stack'); model.Stack={}; end
   i2=find(strncmpi(model.Stack(:,1),'curve',5));
   if isempty(i2); sdtw('_nb','no curve found in model')
   else % check unstacked curve string input wheter direct fe_curve command or not
    st1=st(~cellfun(@isstruct,st));
    i3=ismember(st1(:)',model.Stack(i2,2));
    if ~all(i3) % unmatched strings
     st2=lower(cellfun(@(x) sprintf('test%s',x),fe_curve('testlist'),'uni',0));
     st1=st1(~i3);
     for j1=1:length(st1)
      if ~any(cellfun(@(x) strncmpi(st1{j1},x,length(x)),st2)) % check with fe_curve list
       sdtw('_nb','Some input curve names could not be found in model: %s',st1{j1})
      end
     end
    end
   end % curve check
  end
  % additional argument provides a position in curve field for input curves
  if carg<=nargin; 
   ch=varargin{carg}; carg=carg+1;
   if length(st)>1&&length(st)~=length(ch)
    error('Number or curves (%i) and assignment index (%i) mismatch',...
     length(st),length(ch))
   end
  else; ch=[]; 
  end
  if isfield(r1,'curve')&&~isempty(r1.curve)  % warn if replacement
   if strcmp(name,'?')&&length(st)==length(r1.curve);
   else; sdtw('_nb','Modifying existing case entry %s curve',name);
   end
  elseif ~isempty(ch); r1.curve=cell(1,length(ch));
  end
  if isempty(ch); r2=st;
  else; r2=r1.curve; r2(ch)=st;
  end
  if isfield(r1,'def')&&size(r1.def,2)>length(r2) % warn if possible incoherence
   r3=[];
   try;
    if isscalar(r2);r3=stack_get(model,'',r2{1},'g');end
    if ~isempty(r3)&&size(r1.def,2)==size(r3.Y,2)
    end
   catch; r3=[];
   end
   if isempty(r3)
    sdtw('_nb','%s defines %i loads and %i curves are defined',name, ...
     size(r1.def,2),length(r2)); %#ok<WNTAG>
   end
  end
  if any(cellfun(@isempty,r2))&&length(r2)>1
   sdtw('_nb','some input curves are not defined')
  end
  r1.curve=r2;
 end
 if isfield(RO,'out')&&strcmpi(RO.out,'model.Load');model.Load=r1;
 else;Case.Stack{i1,3}=r1;
 end
   
%% #Grav : assembles a gravity load --------------------------------------2
% fe_case(model,'grav','name',struct('dir',[0 0 1]))
elseif comstr(Cam,'grav');
  name=varargin{carg};carg=carg+1;r1=varargin{carg};carg=carg+1;
  try 
   if isfield(r1,'SE'); SE=r1.SE; r1=rmfield(r1,'SE');
   else; SE=fe_case('assemble -matdes 2 -NoT-SE',model);
   end
  rb=feutilb('geomrb',SE,[0 0 0],SE.DOF);
  r1.def=SE.K{1}*rb.def(:,1:3)*r1.dir(:);% Distributed inertia loads
  r1.DOF=SE.DOF;
  Case=stack_set(Case,'DofLoad',name,r1);
 catch; sdtw('_nb','Gravity computation failed');
 end
%% #SetCase --------------------------------------------------------------------
elseif comstr(Cam,'set');[CAM,Cam]=comstr(CAM,4);
% 'SetCase' seems to be obsolete
  [Case,CaseName,CAM,Cam,model]=get_case(CAM,Cam,model);
  Case=fe_case;

%% #Remove -------------------------------------------------------------------
elseif comstr(Cam,'remove')

  st=varargin{carg};carg=carg+1;
  Case=stack_rm(Case,'',st);

%% #Info ---------------------------------------------------------------------
elseif comstr(Cam,'info')

  fprintf('\nCase : %s\n\n',CaseName)
  disp(Case)
  disp(Case.Stack)
  if nargout==0; clear out; return; 
  else;out=Case; 
  end

%% #Pcond : handling of pre-continionner information  ------------------------
elseif comstr(Cam,'pcond')
  Case=stack_set(Case,'pcond',varargin{carg:carg+1});carg=carg+2;
%% #Reset --------------------------------------------------------------------
elseif comstr(Cam,'reset')

  Case=fe_case;

%% #Connection ---------------------------------------------------------------
elseif comstr(Cam,'connection')||comstr(Cam,'rbe3') % advanced SDT connections
  name=varargin{carg}; carg=carg+1;
  eval('model=fe_caseg(sprintf(''%s -name "%s"'',CAM,name),model,varargin{carg:end});');
  Case=stack_get(model,'case',CaseName,'getdata');
  while carg<=nargin&&~ischar(varargin{carg});carg=carg+1;end
else
 if ~ischar(Cam); error('A string command was expected');
 elseif comstr(Cam,'@'); out=eval(CAM);return;
 else; error('%s not a supported CASE entry',CAM);
 end 
end

if carg<=nargin; [CAM,Cam]=comstr(varargin{carg},1);carg=carg+1;
else;Cam='';
end

end % loop on input arguments - - - - - - - - - - - - - - - - - - -
%% #buildu ->sdtweb fe_load buildu

if ~isempty(model)&&~isempty(CaseName)
  if strcmp(CaseName,'model.Case');model.Case=Case;
  else; model=stack_set(model,'case',CaseName,Case);
  end
  if ~isfield(RO,'doSensMatch')
  elseif ischar(RO.doSensMatch)&&~isempty(RO.doSensMatch)
   model=fe_case(model,'SensMatch',RO.doSensMatch);
  elseif RO.doSensMatch
   model=fe_case(model,'SensMatch');
  end
  out=model;
else;out=Case;
end


% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%% #SubFunc ------------------------------------------------------------------
%% #get_case - - -------------------------------------------------------------
function  [Case,CaseName,CAM,Cam,model]=get_case(CAM,Cam,model);

if length(Cam)<4; i1=[];else;i1=strfind(Cam,'case');end
i2=1; Case=stack_get(model,'case'); CaseName='';

if ~isempty(i1) % the keyword case is present  
    [i2,st1,st2]=comstr(Cam(i1(1):end),'case','%i');
    if isempty(i2) && comstr(st2,'case') % Deal with GetCaseCase(i)
     [i2,st1,st2]=comstr(st2,'case','%i');
    end
    if isempty(i2)&&length(Cam)>i1+3
      [CaseName,CAM,Cam]=comstr(CAM(i1:end),'case','%s');
      if ~isempty(Case)&&any(strcmp(CaseName,Case(:,2)))
        i2 =find(strcmp(CaseName,Case(:,2)));
      else;i2=size(Case,1)+1;
      end
    else
      [CAM,Cam]=comstr([CAM(1:i1-1) st1],1);
      if isempty(i2); i2=1;end
    end
else;i2=1;
end
if isempty(CaseName)&&isfield(model,'Case')&&isfield(model.Case,'Stack')
       Case=model.Case;CaseName='model.Case';
elseif i2>size(Case,1)
  if i2>1&&isempty(CaseName)
    sdtw('Case %i undefined, using Case 1',i2);
  end
  Case=fe_case;if isempty(CaseName); CaseName='Case 1';end
else;CaseName=Case{i2,2}; Case=Case{i2,3};
end

if isempty(CaseName); CaseName=sprintf('Case %i',i2); end
if isempty(Case); Case=fe_case; end

%% #getFixDof
function r1=getFixDof(model,Case,r1);
 i2=0; r2=r1; r1=r1{3};
 if isnumeric(r1) 
 elseif isfield(r1,'data') 
  if isfield(r1,'DOF')&&max(r1.DOF)<1; i2=r1.DOF;end
  r1=r1.data;
 elseif isfield(r1,'DOF');r1=r1.DOF; 
 elseif ischar(r1); %allow data to contain the selection string
 else;error('Not supported FixDof format');   
 end
 if ischar(r1) 
  i1=strfind(comstr(r1,-27),'-dof');
  if ~isempty(i1) % deal with -dof if exist
   i2=comstr(r1(i1+4:length(r1)),[-1 0])/100;r1=r1(1:i1-1);
  end
  try; r1=feutil(['findnode' r1(:)'],model); catch; r1=[]; end
  if isempty(r1)
   sdtw('_nb','fixdof %s, could not find nodes',r2{2})
  else
  i2=i2(:)';r1=r1(:,ones(size(i2)))+i2(ones(size(r1,1),1),:);
  r1=r1(:);
  end
 end
%% #safeDofSet (robust implicit definition s)
function r1=safeDofSet(r1,model,Case);
 if ~isfield(r1,'DOF')&&isfield(r1,'sel')&&strcmpi(r1.type,'rigid')
   d1=feutilb('geomrb',model,r1.ori,model.DOF);
   n1=feutil(['findnode ' r1.sel],model);
   d1=fe_def('subDOF',d1,n1);
   if isfield(r1,'dir');d1=fe_def('subdef',d1,r1.dir);end
   r1=sdth.sfield('addmissing',d1,r1);
 elseif ~isempty(r1.DOF)&&min(r1.DOF)<1&&isfield(Case,'Node')
     r1=elem0('VectFromDirAtDofUsed',Case,r1);
 end

%% #cleanUpperCType: clean string names of case types - - --------------------
function out=cleanUpperCType(st,sil);
persistent names
if isempty(names)
 names={'FixDof','DofLoad','DofSet','FSurf','FVol','mpc','rbe3','par','rigid',...
  'SensDof','cyclic','info','map','pcond','pred'};
end
[i1,i2]=ismember(lower(st),lower(names));
if ~all(i1)&&~sil; 
  sdtw('_nb','type %s not recognized as a standard case type', ...
    comstr(st(~i1),-30)); 
end
out=st; out(i1)=names(i2(i1));
