function ks = ofact(k,ind,varargin); 

%OFACT Ofact matrix class constructor
%
%	ks = ofact(k)  builds a ofact object from sparse matrix k
%           ks\b can then be used to solve problem k q = b in an optimized
%           fashion, including repeated solves.
%
%       ofact(k,b) replaces k\b for single solves with proper factor cleanup
%       ofact({{k,T},b} replaces T* ( (T'*k*T) \ (T'*b) ) which is a solve
%         with constrain elimination
%
%       ofact                      lists available methods for this system
%       ofact('method')            returns the current method
%       ofact('MethodName') sets current method
%       ofact('silent') attempts to configure current method to be silent
%
%       stack_set(model,'info','oProp',mklserv_utils('oprop','RealSym'))
%       sdtcheck('PatchMkl') % donwloads MKL_PARDISO solver
% 
%       OBSOLETE methods still supported by ofact are
%
%       ks = ofact(k,ind) builds a ofact matrix from a block of k. The
%	   indices are supposed to be ordered and are a subset of all indices
%       ks = ofact(k.data,k.ind) builds a ofact matrix from its components
%       ks = ofact(k,'lu') stores sparse 'lu' factors for static computations
%       ks = ofact(k,'ch') stores sparse Cholesky factors for static
%         computations
%       ks = ofact(k,'de') uses the 'lu' factorization for small and complex
%          matrices and the ofact form otherwise
%       ks= ofact(k,RuOpt.oProp{:}) : transmits properties
%
%	See also help fe_reduc, fe_eig
%		 sdtweb  ofact

%       Etienne Balmes
%       Copyright (c) 2001-2022 by INRIA and SDTools, All Rights Reserved.
%       $Revision: 1.88 $  $Date: 2025/11/05 15:47:07 $
%       Use under OpenFEM trademark.html license and LGPL.txt library license

%#ok<*NOSEM>
persistent methods silent

if isempty(methods) % initialize method selection

	methods=struct('name','Undefined ofact method');
	methods=ofact('methodde');
 silent=0;

end


if nargin==0

	if nargout==0; ofact('methodlist');
	else;ks=methods.name;
	end
	return
	%ks.ty = 0; ks.ind =[]; ks.data=[]; ks.dinv=[]; ks.l=[]; ks.u=[];

elseif ischar(k)
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%% #TextCommands - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -1
[CAM,Cam]=comstr(k,1);
if comstr(Cam,'close')||comstr(Cam,'clear')
 %% #clear, close  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -2
 met=[];ks=[];
 if nargin<2; met=methods;ind=[];
 else
  try; met=[]; %#ok<TRYNC>
   ks=ind;
   if isa(ks,'ofact');met=ks.method;
   elseif isstruct(ks)&&isfield(ks,'method');met=ks.method;
   end
  end
  if isempty(met); met=methods; end
 end
 if isempty(ks)||isnumeric(ks)
 elseif isa(ks.data,'ofact');ofact('clear',ks.data);% TkT clear
 end 
 if comstr(methods.name,'psldlt'); sky_fig('close');
 elseif ~isstruct(met)
 elseif ~isempty(met.Clear)
  if isa(ind,'ofact'); ks=ind;
  else
   if isempty(ind); ind=-1; end
   ks=empty_fact; ks.ty(2)=ind(1);
  end
  if ischar(met.Clear);eval(met.Clear);clear ks
  else; ks.ty=met.Clear(ks.ty,ks.method);
  end
 end

elseif comstr(Cam,'silent') % - - - - - - - - - - - - - - - - - - -
 %% #silent: Attempt to make current method silent ------------------------ -2
 if nargin>1; st=ind; else; st='on'; end
 if isequal(lower(st),'on'); silent=0; else; silent=1; end
 try;
  st=methods.name;
  if ~strncmp(fliplr(st),'slitu_',6); st=sprintf('%s_utils',methods.name); end
  if ~isempty(which(st));methods=feval(st,'silent',empty_fact,methods,silent);end
 catch; sdtw('Silent failed');
 end

elseif comstr(Cam,'numfact')
%% #NumFact --------------------------------------------------------------- -2
		ks=varargin{1};
		ks=feval([methods.name '_utils'],k,ind,ks);

elseif ~isempty(strfind(Cam,'fact')) % - - - - - - - - - - - - - - - -
%% #Fact ------------------------------------------------------------------ -2
		ks=empty_fact; ks.method=methods;
		ks=feval([methods.name '_utils'],k,ind,ks,varargin{:});

elseif comstr(Cam,'symrenumber') % - - - - - - - - - - - - - - - -

		if isfield(methods,'SymRenumber'); ks=methods.SymRenumber;
		else;ks='symmmd'; end

elseif comstr(Cam,'_sel') % - - - - - - - - - - - - - - - - 
%% #_Sel: keyword based solver selection -2
% ofact('_sel','spfmex 32 .1')
% ofact('_sel','mklserv_utils -silent')
[RO,unu1,st]=cingui('paramedit -DoClean',[ ...
  'opt(#%s#"further options to be declared to ofact, in the form par1 val1... ")' ...
  ],{struct,ind});
 i1=regexp(st,'[0-9]');
 if any(i1); st1=st(i1(1):end); st=st(1:i1(1)-1); else; st1=''; end
 if isempty(RO.opt); ofact(st);
 else % this has not been used yet
  opt=textscan(RO.opt,'%s','delimiter',' ','multipledelimsasone',1); opt=opt{1};
  i1=cellfun(@any,regexp(opt,'[0-9].*')); % make sure numeric values are correct
  opt(i1)=cellfun(@str2double,opt(i1),'uniformoutput',0);
  ofact(st,opt{:});
 end
 if ~isempty(st1) % direct numeric input for low level handling
  if comstr(st,'spfmex') % could be in param
   i1=comstr(st1,-1);
   ofact('spfmex','FactBuild',...
    sprintf('ks=spfmex_utils(''fact'',k,ks,[1 0 size(k,1)/%i %g*size(k,1)]);',i1));
  else; sdtw('_nb','options %s not interpreted for solver %s',st1,st);
  end
 end
 if nargout>0; ks=ofact('method'); ks=ks.name; end
 
elseif comstr(Cam,'stats')
 %% #Stats: provide info on operations -2

 k1=ind;
 % memory, peak memory, solve time
 if k1.ty(1)==5&&isfield(k1.method,'param')&&isfield(k1.method.param,'stats')
  ks={[k1.ty(1:3) k1.method.param.stats]};
  ks{1,2}=sprintf('%s Fact(%i) - Neq: %i | Memory[MB]: %g | Peak Memory[MB]: %g | Time[s]: %g',...
   k1.method.header,k1.ty(2:3),k1.method.param.stats);
 else; ks={[],''};
 end

elseif comstr(Cam,'_iter') % - - - - - - - - - - - - - - - -
    error('place holder for iterative methods');
elseif comstr(Cam,'@') % - - - - - - - - - - - - - - - -
    ks=eval(CAM);return;
elseif comstr(Cam,'cvs') ;
    ks='$Revision: 1.88 $  $Date: 2025/11/05 15:47:07 $';return;
elseif comstr(Cam,'oprop');
%% #oProp : deal with automated oProp building -2
    if length(Cam)>5; fname=comstr(CAM,6);CAM='oprop';Cam='oprop';
      if nargin>1; r1=[{ind},varargin];else;r1={};end
    elseif nargin==1;fname=methods.name;r1={};
    else;
     fname=ind;r1=varargin; 
    end
    f2=which([fname '_utils']);
    if ~isempty(f2)
      [wd,f2]=fileparts(f2);
      [un1,wd]=fileparts(wd);
      if strcmpi(f2,'ofact'); % subfunction located here
       out=feval([fname '_utils'],CAM,r1{:});
      elseif strcmp(wd,'@ofact');
       ks=struct('ty',-1,'ind',[],'data',[],'dinv',[],'l',[],'u',[],'method',[]);
       ks=class(ks,'ofact');
       out=feval(f2,CAM,r1{:},ks);
      else;% a=ofact('opropmklserv','Sym')
       out=feval(f2,CAM,r1{:});
      end
    end
    ks=out;
    
else %if strncmpi(k,'method',6) % - - - - - - - - - - - - - - - - - - -
 %% #Defaults -2
		k0=empty_fact;
		i1=strfind(k,'complex');
		if ~isempty(i1); NeedComplex=1; k(i1+[0:6])='';else;NeedComplex=0; end


		ks=[];
		ks=spfmex_utils(['method' CAM],k0);ks.param=[];
  ks(end+1).name='mklserv_utils';
  if sp_util('issdt')&&exist('mklserv_utils','file')==2&&exist('mklserv_client','file')
   r2=mklserv_utils('method');r2.param.param(67)=0;% silent
   st=fieldnames(r2);
   for j1=1:length(st);ks(end).(st{j1})=r2.(st{j1});end
   ks=ks([2 1]); % make mklserv default
  else;
   ks(end).header='MKL/PARDISO : Download with sdtcheck(''patchMkl'')';
   ks(end).SymRenumber='';
   ks(end).FactBuild='';
   ks(end).Available=exist('psldlt','file')==3;
   ks(end).HandlesComplex=0;
  end
		ks(end+1).name='lu';
		ks(end).SymRenumber='symmmd';
		ks(end).header='MATLAB sparse LU solver';
		ks(end).FactBuild='ks=lu_fact(k,ks);';
		ks(end).Available=1;
		ks(end).HandlesComplex=1;
        
		ks(end+1).name='ldl';
		ks(end).header='Matlab LDL';
		ks(end).SymRenumber='';
		ks(end).FactBuild='ks=ldl_fact(k,ks);';
		ks(end).Solve ='q=k.data.S*k.data.P*(k.u\( k.data.D \ (k.l\ (k.data.P''*k.data.S*b))));';
		ks(end).Available=~isempty(which('ldl'));
		ks(end).HandlesComplex=0;

		ks(end+1).name='chol';
		ks(end).SymRenumber='symmmd';
		ks(end).header='MATLAB sparse Cholesky solver';
		ks(end).FactBuild='ks=chol_fact(k,ks);';
		ks(end).Available=1;
		ks(end).HandlesComplex=1;

		ks(end+1).name = 'umfpack';
		ks(end).SymRenumber = '';
		ks(end).header = 'UMFPACK solver';
		ks(end).FactBuild = 'ks = umf_fact( k, ks );';
		ks(end).Available = exist('umfpack','file')==3||sdtdef('verm')>=900;
		ks(end).HandlesComplex = 1;
		ks(end).Solve ='cInt=k.data.P*(k.data.R\b);q=k.data.Q*(k.u\(k.l\cInt));clear cInt;';
		ks(end).Clear ='';

		ks(end+1).name='sp_util';
		ks(end).SymRenumber='symrcm';
		ks(end).header='SDT skyline solver';
		ks(end).FactBuild='ks=spu_fact(k,ks);';
		ks(end).Available=(sp_util>5);
		ks(end).HandlesComplex=0;

		ks=AppendToList(ks,'mumps_utils',k0);
		ks=AppendToList(ks,'mtaucs_utils',k0);
		ks=AppendToList(ks,'pardiso_utils',k0);
  if strncmpi(k,'method',6); [k,Cam]=comstr(k,7);
  else; [k,Cam]=comstr(k,1); 
  end
  RO.silent=0;if ~isempty(k)&&k(end)==';';RO.silent=1;Cam(end)='';k(end)='';end
  % Eliminate real solvers if complex matrix
  if nargin==2&&isnumeric(ind) && ~isreal(ind)
   ks=ks([ks.HandlesComplex]~=0);
  end
  
  if isempty(Cam); st1=Cam; % current method
  else; % allow tokens in string
   st1=textscan(Cam,'%s','delimiter',' \t\b\n'); st1=st1{1};
   if isempty(st1); st1=Cam; else; st1=st1{1}; end
  end
  i1=find(strcmpi(st1,{ks.name}));
   
   if nargin>=2&&ischar(ind); args=[{ind} varargin{:}];else;args={};end
   carg=1;
   % Init for file located elsewhere
   st1=sscanf(k,'%s',1);if strcmpi(st1,'ch');st1='chol';end
   if strcmpi(st1,'de')
   elseif isempty(i1)&&(exist(st1,'file')||exist(st1,'builtin'))  
     % File located elsewhere
     ks=AppendToList(ks,k,k0,args{:});carg=length(args)+1;
     i1=find(strcmpi(st1,{ks.name}));
   end
   if comstr(Cam,'list') % MethodList - - - - - - - - - - - - - - - -
			if nargout==1; out=ks([ks.Available]~=0);ks={out.name};return;
			end

			fprintf('\nAvailable factorization methods for OFACT object\n\n');
			ind=1:length(ks);
			for j1=1:length(ind)
				if ks(j1).Available
					if strcmp(methods.name,ks(j1).name)
						fprintf('->%12s : %s\n',ks(j1).name,ks(j1).header);
					else
						fprintf('%14s : %s\n',ks(j1).name,ks(j1).header);
					end
				else
					fprintf('%14s : %s (NOT AVAILABLE ON THIS MACHINE)\n', ...
						['*' ks(j1).name],ks(j1).header);
				end
			end

			fprintf('\nuse ofact(''method MethodName'') to select.\n');

			clear ks; return;
    % - - - - - - - - - - - - - - - - - - - - - - - -
    % ofact('method target','command')
    elseif length(i1)==1&&ks(i1).Available

     ks=ks(i1(1));
     while carg<=length(args)
      if isfield(ks,args{carg}) % set a field
       ks=setfield(ks,args{carg},args{carg+1});carg=carg+2; %#ok<SFLD>
      else;
       st=ks.name;if isempty(strfind(st,'utils'));st=[st '_utils'];end %#ok<AGROW>
       eval(sprintf('ks=%s(args{carg},k0);',st));carg=carg+1;
      end
     end

  % - - - - - - - - - - - - - - - - - - - - - - - -
   elseif ~isempty(Cam)&&~any(strcmp(Cam,{'de','sde','inv','diag'}))
       error('''%s'' unknown',Cam);
   else % method command not matched
    r1=methods;
    if ~isfield(r1,'Available')||~r1.Available %% #MethodsDe : default is first
        for j1=1:length(ks)
            if ks(j1).Available; methods=ks(j1);methods=ks(j1);break;end
        end
    end
    ks=methods;
   end
  
  if ~isempty(strfind(Cam,'-silent')); RO.cbksil=1; else; RO.cbksil=0; end
  if ~isempty(strfind(Cam,'-v')); RO.cbksil=-1; end
  if nargout==0
    methods=ks;
    if ~RO.silent&&RO.cbksil<=0;disp(ks);end
    if iscell(ks); return; % ofact('mklserv_utils','oProp')
    elseif strcmp(ks.name,'pardiso')
        warning('Pardiso is set for symmetric matrices by default')
    end
    clear ks;
  end
  if RO.cbksil>0; ofact('silent'); elseif nargout==0; ofact('silent','off'); end
%else; error('%s not implemented',CAM);      
end
return; % end string commands
end
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% #actual_object_building -1

if nargin==2&&ischar(ind)&&~ischar(k)&&~any(strcmpi(ind,{'diag','lu','inv'})); 
    un1=ofact(['method' ind]);
elseif ~any(exist(methods.name,'file')==[2 3 5 6]);un1=ofact('methodsde');
end


%% #ofact(k),res
if isa(k,'ofact'); ks = k;
 if nargin==2; try; o1=ks\ind; ks=o1; return; catch; error('k is already a ofact object'); end; end
 i1=find(strcmpi(varargin,'movefromcaller'));
 % Handle resolution with RHS memory duplication control
 if ~isempty(i1)&&length(varargin)>i1&&varargin{i1+1} % this is a MoveFromCaller call
  if isfield(k.method.param,'MoveFromCaller')&&(k.method.param.MoveFromCaller&&...
    k.ty(1)==5&&~isempty(k.method.Solve)&&~ischar(k.method.Solve))
   % available for external solver which support a MoveFromCaller call
   b=inputname(2); if ~isempty(b); evalin('caller',sprintf('clear %s',b)); end
   ks=[]; ks=k.method.Solve('ind',k.ty,k.method);
  else % then do a standard solve
   ks=k\ind;
  end
 else; error('This is a wrong Ofact call.');
 end

%% #ofact(ksky,ind)  give the sp_util fields in the input
elseif nargin==2 && size(k,2)==1 && ind(end)==size(k,1)

	if ~isreal(k)
		error('Complex ofact matrices are not handled by sp_util');
	end
	ks.ty = 0; ks.ind=ind; ks.data=k;
	if ks.ind(end)~=length(ks.data)||ks.data(end)~=0
		warning('.ind and .data fields don''t seem to be consistent');
	end
	ks.dinv=[]; ks.l=[]; ks.u=[]; ks.method=ofact('method sp_util');
	ks=class(ks,'ofact');

elseif ~iscell(k)&& size(k,1)~=size(k,2) && nargin==1
	error('Ofact matrices are always square');
%% #ofact({k,T} ...) % Eliminate constraints and solve/factor
elseif ~isa(k,'double')&&~isa(k,'sparse')&&~isa(k,'v_handle')

  if iscell(k) % {k,T},b  % eliminate constraints and solve
   if nargin>2&&isfield(varargin{1},'iter')
     T=k{2};Tt=T'; k=Tt*k{1}*T;  ind=Tt*ind;	
     ks=ofact(k,ind,varargin{:});
     ks=T*ks;
   elseif nargin==1||(~isnumeric(ind)||isempty(ind))
     % #ofact({k,T},[],oProp{:})factor that integrates contraint elimination
     T=k{2};
     if isa(k{1},'ofact'); 
       kd=k{1};if size(kd,1)~=size(T,2);error('Mismatch');end
     else
      varg=varargin; if ischar(ind); varg=[ind,varg]; end
      Tt=T';k=Tt*k{1}*T; kd=ofact(k,[],varg{:},'T',T);
      % Tt=T';k=Tt*k{1}*T; kd=ofact(k,ind,varargin{:},'T',T);
     end
     if isfield(kd.method,'TktSolve')&&~isempty(kd.method.TktSolve)
      ks=kd; ks.ind=T;
      i2=find(cellfun(@(x)ischar(x)&&strcmpi(x,'indDofSet'),varg));
      if ~isempty(i2);
          ks.method.indDofSet=varg{i2+1};
      end
     else
   	  ks=struct('ty',[11 size(T,1),size(T,1)], ...
 		'ind',T,'data',kd,'dinv',[],'l',[],'u',[]);
    
      ks.method=kd.method; ks.method.name='tkt';
      ks.method.header=sprintf('TKT %s',kd.method.header);
      ks.method.Solve=@tkt_solve;ks.method.Clear='';
	  ks=class(ks,'ofact');
     end
   elseif length(varargin)>1&&ischar(varargin{1})&&strcmpi(varargin{1},'method')&&isnumeric(ind)
    %% basic implementation with contraint elimination and clear, provided method  
    T=k{2};Tt=T'; k=Tt*k{1}*T;  ind=Tt*ind;	
    kd=ofact(k,varargin{:}); ks=kd\ind; ofact('clear',kd);
    ks=T*ks;
   elseif isfield(methods,'TktSolve')&&~isempty(methods.TktSolve) 
       eval(methods.TktSolve);
   else % basic implementation of solver with constrain elimination
    T=k{2};Tt=T'; k=Tt*k{1}*T;  ind=Tt*ind;	
    kd=ofact(k); ks=kd\ind; ofact('clear',kd);
    ks=T*ks;
   end
  else;	error('Ofact objects are created from sparse or full matrices')
  end 

%% #ofact(k)
elseif nargin==1 % standard call

	ks=struct('ty',[5 0 0], ...
		'ind',[],'data',[],'dinv',[],'l',[],'u',[]);
	ks.method=methods;
	if ~isreal(k)&&~ks.method.HandlesComplex
		ks.method=ofact(['method de'],k);
	end
	ks = class(ks,'ofact');
    if isa(k,'v_handle');k=k.GetData;end
    if ischar(methods.FactBuild);eval(methods.FactBuild);
    else;ks=methods.FactBuild(k,ks,methods);
        ks = class(ks,'ofact');
    end
    
elseif nargin==2 
%% #ofact(k,'method') standard call with method selection

	if ischar(ind)&& strcmp(ind,'diag') % #diag -2

		r1=full(diag(k)); if nnz(r1)==length(r1); ind=[];else;ind=find(r1);end
		ks=struct('ty',[6 size(k)], ...
			'ind',ind,'data',[],'dinv',1./r1(r1~=0),'l',[],'u',[], ...
			'method','diag');
		if isempty(ind); ks.ty(1)=6.1;end
		ks = class(ks,'ofact');
    elseif ischar(ind)&& strcmp(ind,'inv') % full inverse

		ks=struct('ty',[6.2 size(k)], ...
			'ind',[],'data',[],'dinv',inv(full(k)),'l',[],'u',[], ...
			'method','diag');
		ks = class(ks,'ofact');

	elseif ischar(ind) % #other -2

		ks=struct('ty',[5 0 0], ...
			'ind',[],'data',[],'dinv',[],'l',[],'u',[]);
		ks.method=ofact(['method ' ind]);
		if ~isreal(k)&&~ks.method.HandlesComplex
			ks.method=ofact(['method de'],k);
		end
		ks = class(ks,'ofact');
		r2=ks.method.FactBuild; 
        if ischar(r2);eval(r2);
        else;  ks=feval(r2,k,ks,ks.method); 
           if isstruct(ks);ks = class(ks,'ofact');end
        end

	elseif isnumeric(ind) && min(size(k))==1  % OBSOLETE ofact call
		ks=struct('ty',[0 0 0], ...
			'ind',ind,'data',k,'dinv',[],'l',[],'u',[]);
		ks.method=ofact('method sp_util');
		if ~isfield(ks,'dinv'); ks.dinv=[]; end
		ks = class(ks,'ofact');
		% Single solve with clear object
	elseif isnumeric(ind)&&size(ind,1)==size(k,2) % #directSolve -2

		kd=ofact(k); ks=kd\ind; ofact('clear',kd);

	elseif isnumeric(ind) % OBSOLETE ofact call
		k=k(ind,ind); sp_util('spind',k,ind);
		ks = sp_util('sp2sky',k,rind,length(ind));
		if ~isfield(ks,'dinv'); ks.dinv=[]; end
		ks = class(ks,'ofact');
	end

else
 %% #ofact(k,...)_method with argument passing -------------------------------1
 if isstruct(varargin{1})&&isfield(varargin{1},'iter')
  %% #ofact.iterative -2
  RO=varargin{1};
  if isfield(RO,'pcond') % 'pcond',ofact('@cgs_pcond')
   RO=feval(RO.pcond{1},k,RO,RO.pcond{2:end});
  end
  b=ind; clear ind %ks=ind;
  if isfield(RO,'multirhs')&&RO.multirhs
   ks=feval(RO.iter,k,b,RO.iterOpt{:});
  else
   ks=b;
   for j1=1:size(b,2)
    ks(:,j1)=feval(RO.iter,k,b(:,j1),RO.iterOpt{:});
   end
  end
  return;
 end
 ks=struct('ty',[5 0 0], ...
  'ind',[],'data',[],'dinv',[],'l',[],'u',[]);
 if ischar(ind)&&strcmpi(ind,'method') %ofact(k,'method','value') ofact(k,oProp{:})
  ks.method=ofact(['method' varargin{1}]);carg=2;
 elseif nargin>3&&strcmpi(varargin{1},'method') %ofact(k,b,'method','value')
  ks.method=ofact(['method' varargin{2}]);carg=3;
 else;ks.method=methods;carg=1;% use current default method
 end
 k1=0;
 % ofact(k,'movefromcaller',1)
 if strcmpi(ind,'movefromcaller')&&~isempty(varargin)&&varargin{carg}; k1=1;carg=carg+1;
 else % ... 'movefromcaller',1)
  i1=find(strcmpi(varargin,'movefromcaller'));
  if ~isempty(i1)&&length(varargin)>i1&&varargin{i1+1};  k1=1;  end
 end
 if k1; % MoveFromCaller Mode
  k1=inputname(1);
  if ~isempty(k1); evalin('caller',sprintf('clear %s',k1)); end
 end
 
 if ~isreal(k)&&~ks.method.HandlesComplex
  error('%s does not handle complex',ks.method.name)
 end
 method=ks.method; ks = class(ks,'ofact');
 % %/!\ FactBuild may clear k !!
 s_k=size(k);
 if isa(k,'v_handle');k=k.GetData;end
 if ischar(method.FactBuild);eval(method.FactBuild);
  %     elseif carg==1
  %       ks=method.FactBuild(k,ks,method,ind,varargin{:});
  %       ks = class(ks,'ofact');
 else;
  ks=method.FactBuild(k,ks,method,varargin{carg:end});
  ks = class(ks,'ofact');
 end
 i1=find(cellfun(@(x)ischar(x)&&strcmp(x,'T'),varargin));
 if ~isempty(i1)&&(~isfield(ks.method,'TktSolve')||isempty(ks.method.TktSolve)||ischar(ks.method.TktSolve))
  % TkT solver not yet present
  T=varargin{i1+1};ks.method.TktSolve=T; ks.ty(7)=size(T,1);
 end
 
 if ~isnumeric(ind)||isempty(ind)
 elseif size(ind,1)==s_k(2) %q4=ofact(k+1e3*b,b,'method','lu')
  kd=ks; ks=kd\ind; ofact('clear',kd);
 else; error(sprintf('Inconsistent size : k(%i*%i)\\b(%i,%i)',ks.ty([3 3]),size(ind))) %#ok<SPERR>
 end
end
% #.ty contains
%  [Type FactorNumber size(k,1)]
%% #SubFunc ------------------------------------------------------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function out=sky_fig(Cam)

gf=findall(0,'tag','ofact');
if isempty(gf)
	gf=figure('name','ofact','visible','off','tag','ofact', ...
		'integerhandle','off','HandleVisibility','off', ...
		'CloseRequestFcn','ofact(''close'')');
end

uf=get(gf,'userdata');
if ~isfield(uf,'LDLT'); uf.LDLT=[]; end

if comstr(Cam,'new')
	uf.LDLT(end+1)=length(uf.LDLT)+1;
	out=uf.LDLT(end);
elseif comstr(Cam,'close')
	for j1=1:length(uf.LDLT)
		psldlt('close',j1);
	end
end

%% #lu_utils  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function out=lu_utils(varargin) %#ok<DEFNU>

[CAM,Cam]=comstr(varargin{1},1);carg=2;
if comstr(Cam,'oprop');[CAM,Cam]=comstr(Cam,6);
 out={ ...
 'Gen',{'method','lu'}
     };
 if isempty(Cam)&&carg<=nargin&&ischar(varargin{carg});% a=ofact('oprop','spfmex','Sym')
     [CAM,Cam]=comstr(varargin{carg},1);carg=carg+1;
 end
 if ~isempty(Cam)
     out=out{strcmpi(out(:,1),CAM),2};
 end
else; error('%s',CAM);
end

%% #lu_fact  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -3
function ks=lu_fact(k,ki) %#ok<DEFNU>

ks.ty = 2; ks.ind=[];  ks.data=[]; ks.dinv=[];
if issparse(k); [ks.l,ks.u] = lu(k,sdtdef('luThres')); else; [ks.l,ks.u] = lu(k); end
ks.method=ki.method;
ks = class(ks,'ofact');


%% #ldl_fact - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function ks=ldl_fact(k,ki) %#ok<DEFNU>

[L,D,P,S] = ldl(k); 
ks=struct('ty',[10, 0, size( k, 1 )], 'ind',[],'data', ...
    struct('P',P,'D',D,'S',S),'dinv',[],'l',L,'u',L');
ks.method=ki.method;
ks = class(ks,'ofact');

% A=inv(k.data.S)*k.data.P*k.l*k.data.D*k.u*k.data.P'*inv(k.data.S)
%P'*S*A*S*P = L*D*L'. 

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function ks=chol_fact(k,ki) %#ok<DEFNU>

ks.ty = 3; ks.ind=[]; ks.data=[]; ks.dinv=[]; ks.l=[];
ks.u = chol(k); ks.l = ks.u';
ks.method=ki.method;
ks = class(ks,'ofact');

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function ks=spu_fact(k,ki) %#ok<INUSD,DEFNU>


ks = sp_util('sp2sky',k);
ks.method=ofact('sp_util');
ks = class(ks,'ofact');

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function ks=empty_fact

ks=struct('ty',[0 0 0], ...
	'ind',0,'data',[],'dinv',[],'l',[],'u',[]);
ks.method='';
ks = class(ks,'ofact');

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function ks=psldlt_fact(k,ki) %#ok<DEFNU>

i1=sky_fig('new');
ks=struct('ty',[5 i1 size(k,1)], ...
	'ind','psldlt','data',[],'dinv',[],'l',[],'u',[]);
psldlt('factor',k,i1);
ks.method=ki.method;
ks = class(ks,'ofact');

%% #umf_utils  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function out=umf_utils(varargin) %#ok<DEFNU>

[CAM,Cam]=comstr(varargin{1},1);carg=2;
if comstr(Cam,'oprop');[CAM,Cam]=comstr(Cam,6);
 out={ ...
 'Gen',{'method','umfpack'}
     };
 if isempty(Cam)&&carg<=nargin&&ischar(varargin{carg});% a=ofact('oprop','spfmex','Sym')
     [CAM,Cam]=comstr(varargin{carg},1);carg=carg+1;
 end
 if ~isempty(Cam)
     out=out{strcmpi(out(:,1),CAM),2};
 end
else; error('%s',CAM);
end

%% #umf_fact - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -3
function ks=umf_fact(k,ki) %#ok<DEFNU>

if exist('umfpack','file');[L,U,P,Q,R]=umfpack(k);
else;if ~issparse(k);k=sparse(k);end
    [L,U,P,Q,R]=lu(k);
end
ks=struct('ty',[10, 0, size( k, 1 )], 'ind',[],'data',struct('P',P,'Q',Q,'R',R),'dinv',[],'l',L,'u',U);
ks.method=ki.method;
ks = class(ks,'ofact');

%% #ichol_pcond - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function RO=ichol_pcond(k,RO,varargin) %#ok<DEFNU>
if ~isreal(k);k=real(k);end
if nargin>3
 M=ichol(k,varargin{:}); %Incomplete Cholesky factorization as preconditioner
else
 %% attempt at find
 alpha=max(sum(abs(k),2)./diag(k))-2;
 try
    opts=struct('type','nofill','michol','on');
    tic;M=ichol(k,opts);fprintf('ichol(k,%s) \nDone in %.1f s\n',comstr(opts,-30),toc);
 catch 
 try
     tic;M=ichol(k);%opts('type','nofill')
     fprintf('ichol(k) \nDone in %.1f s\n',toc)
 catch
 try
    opts=struct('type','ict','droptol',1e-3);
    tic;M=ichol(k,opts);fprintf('ichol(k,%s) \nDone in %.1f s\n',comstr(opts,-30),toc);
 catch
 try
    opts=struct('type','ict','droptol',1e-3,'diagcomp',0.3);
    tic;M=ichol(k,opts);fprintf('ichol(k,%s) \nDone in %.1f s\n',comstr(opts,-30),toc);
 catch
 try
   opts=struct('type','ict','droptol',1e-3,'diagcomp',alpha);
   tic;M=ichol(k,opts);fprintf('ichol(k,%s) \nDone in %.1f s\n',comstr(opts,-30),toc);
 catch; warning('Did not find a precondition choice');
                assignin('base','matrix_to_investigate',k)
 end
 end    
 end 
 end    
 end
end



if strcmpi(char(RO.iter),'gmres');RO.iterOpt{4}=M; 
else; RO.iterOpt{3}=M;  RO.iterOpt{4}=M'; % Left and right precond for sym
end

%% #tkt_solve  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function q=tkt_solve(kd,b) %#ok<DEFNU>

T=kd.ind; kr=kd.data;  
q = T * (kr\ (T'*b));

%% #AppendToList : adds external methods to list
function ks=AppendToList(ks,st,k0,varargin);
    try; 
        [st1,i2,i3,i4]=sscanf(st,'%s',1); % allow argument passing for default
        if exist(st1,'file'); 
        elseif exist(['ofact/' st1 '.m'],'file')
        else;return;
        end
        r1=feval(st1,['method' st(i4:end)],k0,varargin{:});
        st=fieldnames(r1);ks(end+1).name=st;
        for j1=1:length(st)
            ks(end).(st{j1})=r1.(st{j1});
        end
    catch;
        fprintf('%s not found , %s',st,lasterr);
    end

%% #check_k - - --------------------------------------------------------------
% Verify the presence of zero diagonal terms
function k=check_k(k); %#ok<DEFNU>
if ischar(k); eval(iigui(k,'MoveFromCaller')); end
if ~issparse(k); k = sparse(k); end
i1=find(diag(k)==0);
if ~isempty(i1)
    %i1=~any(full(k(:,i1)),1));
    error(sprintf('%i zero terms on diagonal :%s\n', ...
       length(i1),sdtw('_clip 50 1','%i ',i1))); 
end
