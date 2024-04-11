function [o1,o2]=fe_load(varargin)

%FE_LOAD Construction of loads (input shape matrices)
%
%	Syntax : load = fe_load(model,Case)
%
%       This function evaluates the CASE (see sdtweb('case')) using the
%       defined MODEL (see sdtweb('fem'))
%       
%       CASE is a structure with fields CASE.Stack and Case.DOF
%
%	Supported loads are DofLoad, DofSet, FVol and FSurf see more details
%         in sdtweb('fe_load')
%
%       See sdtweb     fe_load, sdt/case, femk, fe_case
%       See also help  fe_case, fe_c

%	Etienne Balmes
%       Copyright (c) 2001-2024 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*NOSEM>
epsl=sdtdef('epsl');
if nargin==0; help('fe_load'); return; end
% string commands - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if ischar(varargin{1}) 
 [CAM,Cam]=comstr(varargin{1},1);carg=2;
%% #buildu building the time dependence - - - - - - - - - - - - - - - - - - -
if comstr(Cam,'buildu'); [CAM,Cam]=comstr(CAM,7);

 if carg<=nargin; model=varargin{carg}; carg=carg+1; end
 if carg<=nargin; Case=varargin{carg}; carg=carg+1; end

 if comstr(Cam,'nobuild');o1=varargin{carg}; carg=carg+1;o2=[];
     [CAM,Cam]=comstr(CAM,8);
 elseif isfield(model,'Load')&&isfield(model.Load,'def') % Preassembled load
     o1=model.Load; 
 else; o1=fe_load(model,Case);o2=[];
 end
 if comstr(Cam,'time');

  if carg<=nargin; t=varargin{carg}; carg=carg+1; end
  if carg<=nargin; opt=varargin{carg}; carg=carg+1; end
  if ~isstruct(opt); error('data must be a structure'); end
  RunOpt=struct('bset',0); 
  % if static case with no curve stepping defined
  if isempty(t) && isempty(opt.Opt)&& ...
    (~isfield(o1,'curve')||isempty(o1.curve)|| ...
     (iscell(o1.curve)&&all(cellfun('isempty',o1.curve))))
   if ~isfield(o1,'def')
   elseif size(o1.def,1)==size(Case.T,2); o1=o1.def; 
   elseif isempty(o1.def);o1=zeros(size(Case.DOF)); % No load given set 0
   else;o1=Case.T'*o1.def;
   end
   if isfield(model,'NL')&&~isempty(model.NL)&&strcmp(model.NL{1,2},'bset')
    o2=struct('t',model.NL{1,3}.table(:,1),'Cb',@(ft,fc,j1)zeros(size(o1,1),1));
   else;   o2=ones(1,size(o1,2));
   end
   return;
  elseif isempty(t);t=0; % if static with steps use t
  end
  try;
   if isstruct(o1) && isfield(o1,'def'); opt.nc=size(o1.def,2);
   else;                                opt.nc=size(o1,2); 
   end
   ft=[]; opt.nt=length(t);% each column for each fc.def
   % DofSet application in time domain
   if ~isstruct(o1) 
   elseif isfield(o1,'BIN') &&  isempty(o1.def)&&isfield(o1.bset,'curve')
    o1.curve=o1.bset.curve; 
    if norm(o1.BIN{1},'inf')||norm(o1.BIN{2},'inf');
        warning('Ignoring Mci,Cci in DofSet and assuming displacement input');
    end
    o1.def=-o1.BIN{3}; % BIn = T'*(K)*TIn hence as load need negative sign
    if size(o1.def,1)~=size(o1.DOF,1)&&size(o1.def,1)==size(Case.DOF,1)
      o1.DOF=Case.DOF;
    end
   elseif isempty(o1.curve)&&isfield(o1,'bset')&&isfield(o1.bset,'curve')
     o1.curve=o1.bset.curve;RunOpt.bset=1;
   end
   if isstruct(o1) && isfield(o1,'curve') &&  ~isempty(o1.curve) && ...
     (~isfield(o1,'adof')||~iscell(o1.adof)) % restrain to supported cases
    ft=o1.curve;
    if iscell(ft);
     if ~RunOpt.bset&&length(ft)~=size(o1.def,2)
       disp(o1); 
       error('Inconsistent : %i loads and %i time signal', ...
         size(o1.def,2),length(ft));
     end
     for j1=1:length(ft);
      if ~isempty(t)&&~isequal(t,0) % If curve possibly get t vector there
      else
       if isempty(ft{j1})
        warning('Empty load.curve{%i} seems inconsistent in this context',j1);
       end
       if isfield(ft{j1},'X');t=ft{j1}.X;if iscell(t);t=t{1};end;end
      end
      if ischar(ft{j1})&&~isempty(ft{j1}) % either struct or stack curve name or interpreted by fe_curve
       r1=stack_get(model,'curve',ft{j1},'get');
       if isempty(r1)&&isfield(model,'nmap')&&model.nmap.isKey(ft{j1})
        r1=model.nmap(ft{j1}); if ischar(r1); r1=sdtsys('urnsig',r1,t);end
       end
       if isempty(r1);
         st1=sdth.findobj('_sub:',ft{j1});
         r1=stack_get(model,'',st1(1).subs,'get');
         if ~isempty(r1)&&length(st1)>1% Input:Lab
           i3=strcmpi(r1.X{2},st1(2).subs);r1.X{2}=r1.X{2}(i3,:);r1.Y=r1.Y(:,i3);
           ft{j1}=r1;
         else; % Test ricker dt=1e-3
         end
       else
        if ischar(r1); r1=fe_curve(r1);end%         'TestRicker dt=.01 A=1'
         if ~isfield(r1,'name');r1.name=ft{j1};end;ft{j1}=r1; 
         if isequal(t,0)&&isfield(r1,'X');try;t=r1.X{1};end;end% Use curve for stat
         if isfield(ft{j1},'Y')&&isa(ft{j1}.Y,'function_handle')
           ft={struct('Cb',ft{j1}.Y,'t',t)}; break; % allow ft callback
         end
       end
      elseif isempty(ft{j1}); ft{j1}=[]; % empty, set to numeric
      end
      ft{j1}=fe_curve('returny -extrapby0;',ft{j1},t,model);
      % ft{j1}=ft{j1}(:); removed mar 17 (can be inconsistent)
      % -extrapby0  removed Apr 14 : but default in returny if no extrap given
     end
     ft=[ft{:}];
     if length(opt.Opt)>=6&&opt.Opt(6) ;
      % opt(6) was Nf the optional number of time step of the input force
       sdtw('_nb',['opt.Opt(6) support will be discontinued ' ...
           'in a future release of OpenFEM'])
         ft(opt.Opt(6)+1:end)=0;
     end
    else; ft=fe_curve('returny',ft,t,model);% -extrapby0  removed 
    end
   elseif isfield(opt,'TimeVector');  ft=zeros(opt.nt,opt.nc);ft(1:opt.nt,:)=1;
   elseif length(opt.Opt)>=6&&opt.Opt(6) ;
          ft=zeros(opt.nt,opt.nc);ft(opt.Opt(6)+1:end)=0;  % Step of fixed length
   else;  ft=ones(opt.nt,opt.nc);
   end
   o2=ft;
  catch;
   o1=o1.def; 
   if ~isequal(opt.nt,ft{j1}.X);   
     o2=ones(opt.nt,1);
     fprintf('%s\n',lasterr); warning('Building load case : Failed');
   else ;  o2=ft{j1}.Y;
   end
  end
 end % comstr(Cam,'time')
%% #Init define element of a surface load as additional surface elements
elseif comstr(varargin{1},'init')

 model=varargin{carg};carg=carg+1;RunOpt=[];
 [Case,CaseName]=fe_case(model,'getcase');
 if ~isfield(Case,'Stack')||isempty(Case.Stack);o1=model;return;end
 % initialize surface load entries - - - - - - - - - - - - - - - -
 for j1=find(strcmpi('fsurf',Case.Stack(:,1)))'
     r1=Case.Stack{j1,3};
     [mo1,RunOpt]=GetSurfSel(model,r1,RunOpt);
     if isempty(mo1.Elt); error('No elements were selected in face');end
     [C1,DOF]=fe_mknl('init keep nogett nocon -nodepos',mo1); % xxx with empty case should be optimized
     r1.Case=C1; Case.Stack{j1,3}=r1;
 end % stack entries
 % need to initialize volume load entries - - - - - - - - - - - - - - - -
 o1=stack_set(model,'case',CaseName,Case);

elseif comstr(varargin{1},'cvs')
 o1='$Revision: 1.176 $  $Date: 2024/04/03 15:21:08 $'; return;
elseif comstr(varargin{1},'@');o1=eval(varargin{1});
else;error('%s unknown',CAM);
end

return;
else % string commands or traditionnal model,Case
 carg=1;Cam=''; CAM='';
 model=varargin{carg};carg=carg+1;
end

%% #init ---------------------------------------------------------------------
%-----------------------------------------------------------------------------
if ~isstruct(model) % - - - - - - - - - - - - - - - - - - - - - - - - -
  error('you must define a MODEL and a CASE to call FE_LOAD');
else
 RunOpt=struct('NoT',0,'OnSe',0); Case=[]; CaseName='';% Give load on active DOFs
 while carg<=nargin
  r1=varargin{carg}; carg=carg+1;
  if ischar(r1)&&  strcmpi(r1,'not'); RunOpt.NoT=1;
  elseif ischar(r1); [Case,CaseName]=fe_case(model,['getcase' r1]);
  elseif isfield(r1,'Stack'); Case=r1;
  end
 end % Loop on arguments
 if isempty(Case); [Case,CaseName]=fe_case(model,'getcase');end
 % really only case was given
 if isempty(Case.Stack)&&~isfield(model,'Elt') 
         Case=model;model=struct('DOF',Case.DOF);
 elseif isa(Case,'cell')&&size(Case,2)==3; Case=struct('Stack',{Case});
 elseif isstruct(Case)&&isfield(Case,'Stack')&&size(Case.Stack,2)==3
 elseif isfield(model,'NL')&&size(model.NL,2)>1&&any(strcmpi(model.NL(:,2),'bset'))
 else
  if ~isfield(model,'nmap')||~isKey(model.nmap,'NoLoadWarn')
     sdtw('_nb','No load information seems to be specified in the case');
  end
 end
end % ~isstruct

% building the load - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if isfield(model,'Node')
 NNode='NNode=sparse(model.Node(:,1),1,1:length(model.Node(:,1)));';
 node=model.Node;
else; node=1;
end
if ~isfield(model,'DOF')||isempty(model.DOF); 
  model.DOF=feutil('getdof',model);
end

mdof=model.DOF; i1=find(mdof>0);
%nd=sparse(round(rem(mdof(i1),1)*100),fix(mdof(i1)),1:length(i1),100, ...
%    max([fix(max(mdof(i1)))+1 max(node(:,1))]));
N=length(mdof);
%if max(max(nd))>N error('Something wrong with MDOF'); end
b=[]; bset=struct('DOF',[],'def',[]);bset.lab={};
lab={}; curve={}; bi=[];dofi=[]; r1=[];  pl=[]; ID=[];
if isempty(Case.Stack);RunOpt.list={};
else;RunOpt.list=comstr(Case.Stack(:,1),-27);
end
RunOpt.warn={};
% Starting with 2006b, one should initialize, unless LoadInit used to use old
if ~isempty(intersect(RunOpt.list,{'fvol','fsurf'}))&& ...
  sdtdef('OpenFEM.LoadInit',1)
  if isempty(CaseName)&&isfield(model,'Stack')&&~isempty(model.Stack)
    i1=find(strcmpi(model.Stack(:,1),'case'));
    if ~isempty(i1); CaseName=model.Stack{i1,2};end
  end
  if isempty(CaseName);CaseName='Case1';end
  model=fe_load('init',model,RunOpt);[Case,CaseName]=fe_case(model,'getcase');
end
RunOpt.root={};
%% #loop_on_loads -1
for j1=1:size(Case.Stack,1) 

r1=Case.Stack{j1,3};

switch comstr(Case.Stack{j1,1},-27)
case 'fvol' % #fvol -2

 if ischar(r1.sel); [i1,elt]=feutil(['findelt' r1.sel],model);
 elseif ~isfinite(r1.sel(1)); elt=r1.sel; 
 else; error('Data.Sel must contain elements or a selection');
 end
 if isempty(elt); error('Element selection for volume load is empty');end
 if ischar(NNode);eval(NNode);end
 b1=volume_load(model,Case,elt,r1,NNode);

case 'fsurf' % #fsurf -2

 if isfield(r1,'Case') % case entry that was initialized properly
  b1=zeros(size(model.DOF));
  r2=elem0('VectFromDir',model,r1,model.DOF);
  for jGroup=1:size(r1.Case.GroupInfo,1);
   DofPos=r1.Case.GroupInfo{jGroup,1};
   pointers=r1.Case.GroupInfo{jGroup,2};
   constit=r1.Case.GroupInfo{jGroup,4};
   EltConst=r1.Case.GroupInfo{jGroup,8};
   if constit(1)~=-3; pointers(5,:)=int32(100);% surface element
   else;
    switch constit(4,1)
     case 1; pointers(5,:)=int32(100); % volume force
     case 2; pointers(5,:)=int32(102); % volume force, prop to density
     case 3; pointers(5,:)=int32(101); % pressure
      if ~isfield(EltConst,'bas');EltConst.bas=zeros(9,size(EltConst.w,1));end
      i1=double(DofPos(:,1))+1;
      if ~isfield(r2,'AtNode')&& ...
        ~any(unique(round(rem(r2.DOF(i1(i1>0)),1)*100))==19)
       i2=fe_c(r2.DOF,.19,'ind');r2.AtNode=zeros(size(r1.Case.Node,1),1);
       if ischar(NNode);eval(NNode);end
       if nnz(r2.def); % The vectFromDir does not place something on xyz
        r2.AtNode(NNode(fix(r2.DOF(i2))),:)=r2.def(i2,:);
       else;  i2=fe_c(r1.DOF,.19,'ind');
        i3=fix(r1.DOF(i2)); if any(i3==0);error('Not a valid case');end
        r2.AtNode(NNode(i3),:)=r1.def(i2,:);
       end
      end
     case 5; pointers(5,:)=int32(104); % pressure
      if ~isfield(EltConst,'bas');EltConst.bas=zeros(4,size(EltConst.w,1));end
     case 6; pointers(5,:)=int32(105); % pressure
      if ~isfield(EltConst,'bas');EltConst.bas=zeros(4,size(EltConst.w,1));end
     otherwise; sdtw('constit(4,1)=%i not know',constit(4,1));
    end
   end
   if isstruct(r2);r2.Case.jGroup=jGroup;end % r2 = def
   b1=elem0('rhs_og',DofPos,EltConst.NodePos,r1.Case.Node,pointers, ...
    r1.Case.GroupInfo{jGroup,3:7},EltConst,r2,b1);
  end % jGroup
 else % obsolete calls
  if isempty(pl); pl=fe_mat('getpl',model); il=fe_mat('getil',model);end
  if ischar(NNode);eval(NNode);end
  b1=surface_load(model,pl,il,r1,NNode);
 end
 
case 'dofload' % #dofload -2
 if isfield(r1,'lab')&&ischar(r1.lab); r1.lab={r1.lab}; end
 [r2,r3,b1]=fe_c(model.DOF,r1.DOF);
 if size(r1.def,1)==length(r3)&&(~isfield(r1,'OnSe')||r1.OnSe==1); 
     b1=b1'*r1.def;
 elseif isempty(r2)||isfield(r1,'SE')||(isfield(r1,'OnSe')&&r1.OnSe==0)
  r2=stack_get(model,'SE'); 
  if isempty(r2);
   sdtw('_nb','Load ''%s'' affects no DOF',Case.Stack{j1,2});
   b1=spalloc(length(model.DOF),1,0);
  else;
      eval('r3=fesuper(''sesens -fgetdof'',model,r1,Case);RunOpt.OnSe=1;');
      b1=r3.def;
      r1=r3; 
  end
 else % actually some DOFs are missing
  b1=b1'*r1.def(fe_c(r1.DOF,r2,'ind'),:);
  RunOpt.warn{end+1}=sprintf('%i DOFs were deleted in %s:%s', ...
      length(r1.DOF)-length(r2),Case.Stack{j1,1:2});
 end

case 'dofset' % #dofset -2
  r1=feval(fe_case('@safeDofSet'),r1,model,Case);
  ind=find(rem(round(r1.DOF*100),100)==0); % deal with wild card spec
  if ~isempty(ind);
   i1=fe_c(model.DOF,r1.DOF(ind),'ind');
   r1.DOF(ind)=[];r1.def(ind,:)=[];r1.DOF=[bset.DOF;model.DOF(i1)];
   r1.def(end+(1:length(i1)),1:length(i1))=eye(length(i1));
  end
 if nargout==2 % return a DOFsets as a second ouput

  if isempty(bset.DOF); bset.DOF=r1.DOF;
  else; 
   %bset.DOF=unique(round([bset.DOF;r1.DOF]*100))/100; % wrong for recent matlab!!!
   dof=fe_c(r1.DOF,bset.DOF,'ind',2); % new dof
   if ~isempty(dof); bset.DOF=[bset.DOF;dof]; end
  end
  if max(r1.DOF)>1; 
    st=intersect({'def','v','a'},fieldnames(r1));
    ind=fe_c(bset.DOF,r1.DOF,'ind');
    for j2=1:length(st);
     if ~isfield(bset,st{j2});%bset.(st{j2})=[];end
      bset.(st{j2})=double.empty(length(bset.DOF),0);
     end
     if issparse(bset.(st{j2}))||issparse(r1.(st{j2})) % just in case
      bset.(st{j2})=[bset.(st{j2}) fe_c(bset.DOF,r1.DOF)'*r1.(st{j2})];
     else;  bset.(st{j2})(ind,end+(1:size(r1.def,2)))=r1.(st{j2}); % xxx not sure it performs better
     end
    end
  end
  st=intersect(fieldnames(r1),{'data','KeepDof','sdof','type'});
  for j2=1:length(st); bset.(st{j2})=r1.(st{j2});end%Case with DOF kept
  if ~isfield(r1,'def');
  elseif size(r1.def,2)==1;
   st1=fe_curve('datatypecell','displacement');
   bset.lab(end+1,1:3)={Case.Stack{j1,2} st1{2:3}};
  elseif isfield(bset,'data');bset=feutil('rmfield',bset,'lab');
  elseif size(r1.def,2)>1
   st1=fe_curve('datatypecell','displacement');
   for j2=1:size(r1.def,2) 
    bset.lab(end+1,1:3)={sprintf('%s-u%i',Case.Stack{j1,2},j2) st1{2:3}}; 
   end
  end
  if isfield(r1,'curve');bset.curve=r1.curve;end
  if isfield(r1,'ID');bset.ID=r1.ID;end
  b1=[];

 elseif ~isfield(model,'K')||isempty(model.K)|| ...
   size(model.K{1},1)~=size(model.DOF,1)
    sdtw('_nb',['%s ''%s'' entry ignored because model.K was' ...
      ' not defined or wrong size'], Case.Stack{j1,1:2}); b1=[];
 elseif isfield(Case,'TIn') && size(model.DOF,1)==size(Case.TIn,1)
  b1=spalloc(length(model.DOF),0,1); 
  for j2=1:length(model.K)
   k1=model.K{j2}; if isa(k1,'v_handle'); k1=k1.GetData; end
   try; r2=k1*Case.TIn; 
   catch;eval('r2=feutilb(''a*b'',k1,Case.TIn);')
   end
   b1(:,end+[1:size(r2,2)])=r2;
  end
 else
  ind=fe_c(model.DOF,r1.DOF,'ind');cind=1:length(model.DOF);cind(ind)=0;
  if ~isfield(r1,'def')
  elseif length(ind)==size(r1.def,1)
   cind=find(cind);
   b1=zeros(length(model.DOF),0); 
   for j2=1:length(model.K)
    k1=model.K{j2}; k1=k1(:,ind)'; k1=k1(:,cind)'; % model.K{j2}(cind,ind)
    r2=k1*r1.def; b1(cind,end+[1:size(r2,2)])=r2;
   end
  else; fprintf(' %s %s ignored',Case.Stack{j1,1:2});b1=[];
  end
 end

case {'fixdof','keepdof','par','rigid','sensdof','mpc','nastran', ...
   'rbe3','sensstrain','cyclic','pcond','pred'}
 b1=[];
 otherwise % #other -2
  b1=[]; fprintf('''%s'' entry not supported by fe_load\n',Case.Stack{j1,1})
end

 if     isempty(b1);continue;
 elseif isfield(r1,'Curve'); st1=r1.Curve;
 elseif isfield(r1,'curve'); st1=r1.curve;
 else; st1='';
 end

 % labels are ported from load definition
 if isfield(r1,'lab')&&size(r1.lab,1)==size(b1,2); 
   lab(size(b,2)+(1:size(r1.lab,1)),1:size(r1.lab,2))=r1.lab;ind=size(b,2); 
   if iscell(st1)&&length(st1)==1&&size(b1,2)>1
    sdtw('Using the same curve for all inputs');st1=st1(ones(1,size(b1,2)));
   end
   for j2=1:size(b1,2) 
    if iscell(st1); curve{ind+j2}=st1{j2};else; curve{ind+j2}=st1;end
   end
 elseif size(b1,2)==1; 
  ind=size(b,2); lab{ind+1,1}=Case.Stack{j1,2};  
  if iscell(st1); curve(ind+[1:length(st1)])=st1;
   else; curve{ind+1}=st1;end
 elseif size(b1,2)>1
  ind=size(b,2); RunOpt.root{end+1}=Case.Stack{j1,2};
  for j2=1:size(b1,2) 
   lab{ind+j2,1}=sprintf('%s:u%i',Case.Stack{j1,2},j2); 
   if isfield(r1,'lab');
    try;lab{ind+j2,1}=sprintf('%s:%s',Case.Stack{j1,2},r1.lab{j2});end
   end
   if isempty(st1)&&isempty(curve) % no load defined
   elseif iscell(st1) curve{ind+j2}=st1{min(j2,length(st1))};
   else; curve{ind+j2}=st1;
   end
  end
 end
 %if isfield(Case,'T')&~isempty(Case.T)&size(b1,1)==size(Case.T,1)&~RunOpt.NoT 
 %  b1=Case.T'*b1;
 %end % done below and should thus not be done twice (see clean up below)
 if isfield(r1,'ID'); % Fix the ID if needed
  if length(r1.ID)<size(b1,2); 
    if isempty(r1.ID); r1.ID=1;end
    r1.ID(end+1:size(b1,2))=max([r1.ID(:);max(ID)])+ ...
     [1:size(b1,2)-length(r1.ID)];
  end
  r1.ID=r1.ID(:); ID(end+[1:size(b1,2)],1)=r1.ID(1:size(b1,2));
 end

 if isempty(b); b=b1; else; b=[b b1];end
 b1=[];
end
%% #loop_on_stack -1
if length(RunOpt.root)==1&&all(strncmp(lab(:,1),RunOpt.root{1},length(RunOpt.root{1})))
 lab=regexprep(lab,['^' RunOpt.root{1},':'],'');
end

if ~isempty(RunOpt.warn);fprintf('%s\n',RunOpt.warn{:});end

o1=struct('DOF',model.DOF,'def',b,'lab',{lab},'fun',[0 13], ...
  'ID',ID,'curve',{curve}); if RunOpt.OnSe;o1.OnSe=1;end
% fun(2)=13 is for UFF DataType excitation force
if ~RunOpt.NoT; Case=fe_case(model,'gett',Case);end

o2=bset;
if ~RunOpt.NoT&&isfield(Case,'T') % Project force
 if size(Case.T,1)==size(o1.def,1)&&~isequal(round(o1.DOF*100),round(Case.DOF*100))
  o1.def=Case.T'*o1.def; o1.DOF=Case.DOF;
 elseif size(Case.T,2)==size(o1.def,1)
  o1.DOF=Case.DOF;
 end
end


%% #volume_load --------------------------------------------------------------
function b=volume_load(model,Case,elt,r1,NNode)

b=zeros(size(model.DOF,1),1);
[EGroup,nGroup]=getegroup(elt);
node=model.Node;

if ~isequal(model.Elt,elt)
 model.Elt=elt;model=fe_case(model,'reset');
 [Case,mdof]=fe_mknl('initnoconkeep',model);
elseif ~isfield(Case,'GroupInfo')||isempty(Case.GroupInfo)
 [Case,mdof]=fe_mknl('initnoconkeep',model);
else; mdof=model.DOF;
end
if isfield(r1,'dir')
 r2=elem0('VectFromDir',model,r1);
 %r2=elem0('VectFromDirAtDof',model,r1,model.DOF);
 jGroup=1; 
 % Some adjustments if there is a problem with the field at node
 [ElemF,i1,ElemP]= feutil('getelemf',elt(EGroup(jGroup),:),jGroup);
 switch ElemF
 case {'q4p','t3p','q5p','q8p','q8a','q9a','t3a','t6a','t6p'}
   r2=r2(1:2,:);
 case 'mitc4'
   r2(5,:)=0; 
 end
else; r2=elem0('VectFromDir',model,r1,model.DOF);
end

for j1=1:size(Case.GroupInfo,1);
  if isempty(Case.GroupInfo{j1,5}); Case.GroupInfo{j1,5}=zeros(1,100);end
end
if isfield(r1,'MatDes'); b=fe_mknl('assemble',model,Case,r1.MatDes,r2,NNode);
else; b=fe_mknl('assemble',model,Case,100,r2,NNode);
end
b=fe_c(mdof,model.DOF,b')';

return 

% This is meant to disappear after 2006b
%% #surface_load -------------------------------------------------------------
% see sdtweb elem0('rhs_og')
function b=surface_load(model,pl,il,r1,NNode) %#ok<INUSL>

if isa(model,'v_handle');model=model.GetData;end
elt=model.Elt; node=model.Node;
RunOpt=struct('FaceInd',[]);

if isfield(r1,'eltsel') 
 [i1,elt]=feutil('findelt',node,model.Elt,[],r1.eltsel);
 if isempty(elt)
  error('Element selection for surface load selects no elements');
 end
% dir defined on elements translate to volume load format
elseif isfield(r1,'set')&&~isfield(r1,'def') 

 [elt,i1,i2]=feutil('getedgepatch',model,r1.set);
 [eltid,elt]=feutil('eltidfix',elt);r1.EltId=eltid(i2);
 if size(r1.dir,2)~=length(r1.EltId); error('Inconsistent face and dir');end
 mpid=feutil('mpid',elt);i1=unique(mpid(:,2));i1=i1(i1~=0);model.il=[];
 for j1=1:length(i1); 
  model.il(j1,1:4)=[i1(j1) fe_mat('p_solid','US',3) -3 1];
 end
 model.Elt=elt;  [Case,CaseName]=fe_case(model,'getcase'); 
 r1.sel='groupall';
 Case.Stack={'FVol','vol',r1};
 Case.T=[];Case.GroupInfo={};Case.MatGraph=[];Case.DofPerElt=[];Case.DOF=[];
 model=stack_set(model,'case',CaseName,Case); b=fe_load(model);
 mdof=evalin('caller','mdof'); 
 if ~isequal(b.DOF,mdof); b=fe_c(mdof,b.DOF,b.def')';else; b=b.def; end
 return;
end
model.Elt=elt; 

% Generate surface properties - - - - - - - - - - - - - - - - -
[Case,CaseName]=fe_case(model,'getcase');
i1=(strcmpi('fsurf',Case.Stack(:,1)));
Case.Stack=Case.Stack(i1,:);Case.T=[];
[Case,mdof]=fe_mknl('initnocon',model,Case);
b=zeros(size(mdof,1),1);

[EGroup,nGroup]=getegroup(elt);

% Obsolete selfacei call
if isfield(r1,'set')&&~isempty(r1.set)
 ind=[];nind=[];
 if ischar(r1.set); FaceSet=stack_get(model,'set',r1.set,'GiveData');
 else; FaceSet=r1.set;
 end
 if ~isfield(FaceSet,'data')
  error('Set ''%s'' not found',r1.set);
 end
 
elseif isa(r1.sel,'cell')&&size(r1.sel,1)==1&& ...
                           strcmp(comstr(r1.sel{1,2},-27),'selfacei')
 ind=r1.sel{1,4};
% Node numbers are given
elseif isnumeric(r1.sel); ind=r1.sel;
% Node selection command is given
elseif ischar(r1.sel) && ~isempty(r1.sel) 
  ind=feutil(['findnode' r1.sel],node,elt);
  if isempty(ind); error('No node selected for surface load'); end
% Node selection given as a Stack
elseif isa(r1.sel,'cell')&&size(r1.sel,2)==4
  ind=feutil('findnode',node,elt,[],r1.sel);
else
 error('Currently surface load must be defined with selfacei selection');
end

eltid=feutil('eltidfix;',elt); 
if isempty(ind)
else %if ~isempty(ind) % Edge/face set defined by nodes
 nind(ind)=ind; if max(node(:,1))>length(nind);nind(max(node(:,1)))=0;end
 if isempty(ind);disp(r1); error('No node selected for surface load'); end
 for jGroup=1:nGroup
  EltConst=Case.GroupInfo{jGroup,8};
  % find the associated parent element
  [ElemF,i1,ElemP]= feutil('getelemf',elt(EGroup(jGroup),:),jGroup);
  if isfield(EltConst,'MatrixIntegrationRule') ...
       && ~isempty(EltConst.MatrixIntegrationRule)
   % *b family - - - - - - - - - -
   r2=feutil('selelt selface & innode',node,elt,[],ind);
   i2=find(r2(:,3)==r2(:,4)); 
   if ~isempty(i2)
    i2=intersect(i2,feutil('findelt eltname quad4',node,r2));
    r2(end+1,1:6)=[Inf abs('tria3')];
    r2(end+[1:length(i2)],1:5)=r2(i2,[1:3 5 6]);r2(i2,:)=[];
   end
   if ~isfield(r1,'Elt'); r1.Elt=r2;
   else;r1.Elt(end+[1:size(r2,1)],1:size(r2,2))=r2;
   end
  else % MODULEF strategy with face numbers - - - - - - - - - - - - - - - -
   if i1<0; cEGI=[];  % Ignore groups for display only 
   else; cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;
   end

   if any(strcmp(ElemP, ...  % 2-D elements
    {'q4p','t3p','q5p','q8p','q8a','q9a','t3a','t3p','t6a','t6p'})) || ...
    any(strcmp(ElemF, ...  % 2-D elements
    {'q4p','t3p','q5p','q8p','q8a','q9a','t3a','t3p','t6a','t6p'}))
   i3=feval(ElemP,'edge')';dims=2;
   else; % 3D elements
     i3=feval(ElemP,'face')'; dims=3;
   end
   if ~isempty(cEGI);
    for j1=1:size(i3,2)
     i1=reshape(nind(elt(cEGI,i3(:,j1))),length(cEGI),size(i3,1));
     ind=find(all(i1,2)); % indices of included faces
     if ~isempty(ind); eltid(cEGI(ind),j1+1)=1; end
    end
   end % Loop on faces
   [i3,i4]=find(eltid(:,2:end));
   FaceSet=struct('ID',0,'data',[eltid(i3,1) i4]);
   eltid=eltid(:,1);
  end % *b or MODULEF
 end % Loop on groups
end % Define face/edge set from nodes

% *b family strategy - - - - - - - - - - - - - - - - - - - - - - - - - - -
% WARNING THIS IS REALLY EXPERIMENTAL
if isfield(r1,'Elt'); 
 error('This is now implemented in the init phase');
else % MODULEF strategy - - - - - - - - - - - - - - - - - - - - - -

% build sparse matrix with loading
if isfield(r1,'def')
  spb=sparse(round(rem(r1.DOF,1)*100),fix(r1.DOF)+1,r1.def, ...
   100,max(ceil(max(r1.DOF)),size(node,1))); % sparse reindexing of load
  i1=find(spb(:,1)); for j1=i1(:)'; spb(j1,:)=spb(j1,1); end % xxx should do better
  % then use NNode to point in spb instead of direct
end
faces=sparse(FaceSet.data(:,1),FaceSet.data(:,2),1);

% loop on elements - - - - - - - - - - - - - - - - - - - - - - - - - 
for jGroup=1:nGroup

 % find the associated parent element
 [ElemF,i1,ElemP]= feutil('getelemf',elt(EGroup(jGroup),:),jGroup);
 if i1<0; cEGI=[];  % Ignore groups for display only 
 else; cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;
  [i2,i3]=intersect(eltid(cEGI),FaceSet.data(:,1));
  cEGI=cEGI(i3);
 end
 pointers=Case.GroupInfo{jGroup,2};
 integ=Case.GroupInfo{jGroup,3};
 constit=Case.GroupInfo{jGroup,4};

if ~isempty(cEGI) % Legacy mode that calls modulef elements

  if any(strcmp(ElemF, ...  % 2-D elements
    {'q4p','t3p','q5p','q8p','q8a','q9a','t3a','t3p','t6a','t6p'}))
   FaceIndex=feval(ElemP,'edge')';
  else; FaceIndex=feval(ElemP,'face')'; 
  end
 i4=feval(ElemP,'node');
 i2=max(max(elt(cEGI,i4))); if i2>length(nind); nind(i2)=0;end
 iopt=[0 0 length(feval(ElemF,'dof')) length(i4) 100 0 2 zeros(1,2*size(i3,2))];
 iopt=int32(iopt);

vol=zeros(3,double(iopt(5)));

% Switch on element type - - - - - - - - - - - - - - - - - - - - - - - -

switch ElemF
case {'tetra4','hexa8','penta6','tetra10','hexa20','penta15'}
  % idof : help of_mk 
  % [coor nodeid], volf(3*noe), surf(3 ddl x nnof=4 x nface), 
  %                           pres(nnof=4 x nface)
if length(constit)<4; error('constit initialization failed');end

for jElt = 1:length(cEGI)

  nodeE=node(NNode(elt(cEGI(jElt),i4)),[5:7 1]);
  % define loaded faces
  point=pointers(:,jElt);point(5)=100; point(9:20)=0;

  % flags for loaded surface and pressure faces
  ind=full(faces(eltid(cEGI(jElt)),:));
  if length(ind)<size(FaceIndex,2); ind(size(FaceIndex,2))=0;end
  if length(ind)>size(FaceIndex,2); ind=ind(1:size(FaceIndex,2));end
  point(9:8+size(FaceIndex,2))=ind;
  point(9+size(FaceIndex,2):8+2*size(FaceIndex,2))=ind;

  % build the sur and pre values
  sur = zeros(3*size(FaceIndex,1),size(FaceIndex,2)); i5=size(sur,1);
  pre = zeros(size(FaceIndex,1),size(FaceIndex,2)); 
  i1=reshape(elt(cEGI(jElt),FaceIndex),size(FaceIndex,1),size(FaceIndex,2));
  ind=find(ind);

  sur(1:3:i5,ind)=reshape(spb(1,NNode(i1(:,ind))),size(i1,1),length(ind));
  sur(2:3:i5,ind)=reshape(spb(2,NNode(i1(:,ind))),size(i1,1),length(ind));
  sur(3:3:i5,ind)=reshape(spb(3,NNode(i1(:,ind))),size(i1,1),length(ind));

  pre(:,ind)=reshape(spb(19,NNode(i1(:,ind))),size(i1,1),length(ind));

  bi=of_mk(ElemF,int32(point),integ,constit,nodeE,[sur(:);pre(:)],vol);
  if ~isempty(bi)
    in1=double(Case.GroupInfo{jGroup,1}(:,cEGI(jElt)-EGroup(jGroup)))+1;
    in2=find(in1);in1=in1(in2);    b(in1,1)=b(in1,1)+bi(in2);
  end
end % loop on elements

case {'q4p','t3p','q5p','q8p','q8a','q9a','t3a','t6a','t6p'}
% 2-D elements - - - - - - - - - - - - - - - - - - - - - - - - - -
% (coor[noe,2], fomega[2,noe], fgamma[2,2*nb_edge] (fx,fy), ...
%  press[2*nb_edg], ...
%  no_ref[n_edge,2](int) (i,1) = 0 fgamma is 0, (i,2)=0 press is 0
%  alpha[3] (xx,xy,yy), theta[noe], car, iopt[1](int), be
%
% xxx [idof,bi]=of_mk(ElemF,[ndof noe 100 0 2 iopt no_ref], ...
%            node(i1,[5:7 1]),[fgamma(:),press(:);T],vol);
% The material properties are used for thermal problems.
if length(constit)<4; error('constit initialization failed');end

for jElt = 1:length(cEGI)

  nodeE=node(NNode(elt(cEGI(jElt),i4)),[5:7 1]);

  % define loaded faces
  point=pointers(:,jElt);point(5)=100; point(9:20)=0;

  % flags for loaded surface and pressure faces
  ind=full(faces(eltid(cEGI(jElt)),:));
  if length(ind)<size(FaceIndex,2); ind(size(FaceIndex,2))=0;end
  if length(ind)>size(FaceIndex,2); ind=ind(1:size(FaceIndex,2));end
  point(9:8+size(FaceIndex,2))=ind;
  point(9+size(FaceIndex,2):8+2*size(FaceIndex,2))=ind;

  % build the sur and pre values
  sur = zeros(2*size(FaceIndex,1),size(FaceIndex,2)); i5=size(sur,1);
  pre = zeros(size(FaceIndex,1),size(FaceIndex,2)); temp=pre;
  i1=reshape(elt(cEGI(jElt),FaceIndex),size(FaceIndex,1),size(FaceIndex,2));
  ind=find(ind);
  
  sur(1:2:i5,ind)=reshape(spb(1,NNode(i1(:,ind))),size(i1,1),length(ind));
  sur(2:2:i5,ind)=reshape(spb(2,NNode(i1(:,ind))),size(i1,1),length(ind));
  pre(:,ind)=reshape(spb(19,NNode(i1(:,ind))),size(i1,1),length(ind));
  temp(:,ind)=reshape(spb(20,NNode(i1(:,ind))),size(i1,1),length(ind)); 

  bi=of_mk(ElemF,int32(point),integ,constit,nodeE, ...
    [sur(:);pre(:);temp(:)],vol);
 
  if any(~isfinite(bi)); error(1);end

  if ~isempty(bi)
    in1=double(Case.GroupInfo{jGroup,1}(:,cEGI(jElt)-EGroup(jGroup)))+1;
    in2=find(in1);in1=in1(in2);    b(in1,1)=b(in1,1)+bi(in2);
  end

end % loop on elements

case {'beam1','bar1','beam3','mass1','celas','rigid','mass2','cbush'}
% 1-D elements - - - - - - - - - - - - - - - - - - - - - - - - - -
 sdtw('Surface Load not defined for : %s (group %i) \n',ElemF,jGroup);
 idof=[];bi=[];

otherwise 

 sdtw('Surface Load not yet supported for : %s \n',ElemF);
 idof=[];bi=[];

end
end % of element type and cEGI is not empty

end % of loop on jGroup
end % strategy

b=fe_c(model.DOF,mdof,b')';

%% #GetSurfSel ---------------------------------------------------------------
function     [mo2,RunOpt]=GetSurfSel(model,r1,RunOpt);
 mo2=[];
 if isfield(r1,'set') % select surface elements & dir
    [elt,i1,i2]=feutil('getedgepatch',model,r1.set);
    [eltid,elt]=feutil('eltidfix',elt);r1.EltId=eltid(i2(i2~=0));
    if isfield(r1,'dir')&&size(r1.dir,2)~=length(r1.EltId); 
         error('Inconsistent face and dir');
    end
    RunOpt.LoadType=3;
 elseif isfield(r1,'Elt'); elt=r1.Elt;RunOpt.LoadType=1;%volume force
 else
  mo1=model; if isa(mo1,'v_handle');mo1=mo1.GetData;end
  if isfield(r1,'eltsel'); mo1.Elt=feutil('selelt',mo1,r1.eltsel);end
  if ~isfield(r1,'EltSel'); 
  elseif isa(r1.EltSel,'double');mo1.Elt=r1.EltSel;
  else;mo1.Elt=feutil('selelt',mo1,r1.EltSel);
  end
  if ~isfield(r1,'type');r1.type='';end
  if strcmpi(r1.type,'edge')
    mo1.Elt=feutil('selelt seledge',mo1); 
    if isequal(unique(round(rem(r1.DOF,1)*100)),19); RunOpt.LoadType=5;
    else;  RunOpt.LoadType=6;
    end
  elseif strcmpi(r1.type,'resultant')
    mo2=mo1;mo2.DOF=[];mo2.name='';mo2=fe_case(mo2,'reset');
    if ~isfield(r1,'matdes');r1.matdes=1;end
    mo2=fe_case(sprintf('assemble -matdes %i -se -NoT',r1.matdes),mo2);
    mo1.K=mo2.K;mo1.Klab=mo2.Klab;mo1.DOF=mo2.DOF;
    mo2=mo1; 
    mo1.Elt=feutil('selelt selface',mo1);
    if isempty(mo1.Elt); mo1.Elt=mo2.Elt; % Allow beams
       mo2.name='beam'; % beam load
       mo2.n1=feutil(['getnode groupall &' r1.sel],mo2);
       mo2.adof=fe_c(mo2.DOF,mo2.n1(:,1),'dof');
       return
    else;  mo2.name='surf';  RunOpt.LoadType=1;%volume force
    end 
  else;  %need to fix problem with coincident volume face & surface elt
    mo1.Elt=feutil('selelt selface',mo1);
    if isequal(unique(round(rem(r1.DOF,1)*100)),19); RunOpt.LoadType=3;
    else;  RunOpt.LoadType=1; % sdtweb p_solid : type 1 volume force
    end
  end
  elt=[]; ind=[];
  if isfield(r1,'SurfSel') % SurfSel is either elements or node sel
    if isa(r1.SurfSel,'double'); elt=r1.SurfSel;else; ind=r1.SurfSel;end
  else;ind=r1.sel;
  end
  if ischar(ind)||iscell(ind);
   sti=ind; 
   try ind=feutil('findnode',mo1,ind); % test adv findelt if missed
   catch; 
    if ischar(ind); ind=feutil(sprintf('findnode inelt{%s}',sti),mo1); end
   end
  end
  if isempty(elt);elt=feutil('seleltinnode',mo1,ind);end
 end
 % remove empty groups
 i1=[find(~isfinite(elt));size(elt,1)+1];elt(i1(diff(i1)==1),:)=[];
 if isempty(mo2)
  mo2=r1;mo2.Node=model.Node; mo2.Elt=elt;
  if isfield(RunOpt,'ClipNode')&&RunOpt.ClipNode
     mo2.Node=feutil('getnodegroupall',mo2);
  end
  if ~isfield(model,'DOF'); model.DOF=feutil('getdof',model);end
  mo2.DOF=model.DOF;
 else;
  if ~strcmpi(r1.type,'resultant');error('Not a valid case');
  else
    mo2.n1=feutil('getnodegroupall',mo1.Node,elt);mo2.El1=elt;
    r2=fe_c(mo1.DOF,mo2.n1(:,1),'dof');
    mo2=fe_case(mo2,'reset','FSurf','Groupall', ...
        struct('Elt',elt,'DOF',r2,'def',eye(length(r2))));
    mo2.adof=r2;mo2.L2=fe_load(mo2);
  end
 end
 % fix surface load info
 mpid=feutil('mpid',mo2.Elt);i1=unique(mpid(:,2));i1=i1(i1~=0);mo2.il=[];
 for j2=1:length(i1); 
  mo2.il(j2,1:4)=[i1(j2) fe_mat('p_solid','US',3) -3 RunOpt.LoadType];
 end
