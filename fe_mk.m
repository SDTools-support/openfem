function [m,k,mdof]=fe_mk(varargin)

%FE_MK	assembly of full and reduced FE model matrices
%
%       Syntax: MODEL      = fe_mk(MODEL,'PropertyName',PropertyValue, ...)
%               ...        = fe_mk(node,elt,pl,il,[],adof,opt)
%               [m,k,mdof] = fe_mk( ... ,[0 OtherOptions])
%               [mat,mdof] = fe_mk( ... [Type OtherOptions])
%
%	A finite element model is characterized a data structure MODEL 
%       see sdtweb('fem') or the corresponding arguments node, elt, pl, il.
%
%       Optional/output fields are
%
%       .bas   coordinate system definitions
%       .DOF   DOF definition vector for output(see sdtweb('mdof')). 
%              If defined in the input and opt(2)==1, model.DOF is 
%              left unchanged              
%       .Stack data structure stack where fe_mk looks in particular for
%              CASE definitions which give boundary condition, constraints
%              and loading (see sdtweb('case'))
%
%       If MODEL is used as an output argument, the assembled matrices are
%       stored in the MODEL.K field with the types given in MODEL.Opt(2,:)
%       You can then assemble with a different CASE using 
%           [m,k]=fe_case(MODEL,'assemble')
%
%       Accepted property/value pairs are
%
%	'options',opt   gives a vector with options for model assembly
%	  opt(1) matrix type. 0 mass and stiffness (default), 1 stiffness,
%                2 mass, ...
%	  opt(2) in a standard assembly 
%                0 eliminates DOFs in adof that are not used or not in an
%                  existing model.DOF
%                1 keeps all DOFs in model.DOF/adof with the same order,
%                2 assembles on all possible DOF with CASE based DOF 
%                  elimination (see sdtweb('mpc'))
%	  opt(3) Assembly method (0: standard, 2 disk, 4: type 3 superelement)
%	  opt(4) Numbering optimization 0: (default) nothing done 
%         SYMRCM otherwise. , 1: symmmd, 2: ofact method
%
%       Note that in many applications, fe_mknl is faster than fe_mk
%        it however is less general in terms of superelement support.
%
%	See sdtweb     fe_case, eltfun, elem0
%	See also help  fe_stress, upcom, fe_mat, femesh, feplot, fe_c

%	Etienne Balmes
%       Copyright (c) 2001-2012 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

%#ok<*NOSEM>

if nargin==1&&ischar(varargin{1})
 [CAM,Cam]=comstr(varargin{1},1);
 if comstr(Cam,'cvs');m='$Revision: 1.73 $  $Date: 2024/10/18 09:32:21 $'; 
 elseif comstr(Cam,'mattype');
  i1=comstr(CAM(8:end),[-1]);
  st={1 'K','stiffness';
      2 'M','mass';
      3 'C','viscous damping';
      4 'B','hysteretic damping';
      5 'KG','tangent stiffness'};
   if isempty(i1)&&nargout==0;disp(st);
   else; i2=find(i1==[st{:,1}]);
    if isempty(i2); m=num2str('MatId %i',i1);else;m=st{i2,2};end
   end
 end
 return;
end

%check of inputs

fp=0;
if nargin==0; help fe_mk; return; else; model=varargin{1}; carg=2; end

% deal with backward compatibility issues - - - - - - - - - - - - - - - - -
if ~isstruct(model) 
 % (node,elt,pl,il,tr,adof,opt)
 model=struct('Node',model,'Elt',varargin{2},'pl',[],'il',[],'Stack',[], ...
 'DOF',[]);
 if nargin>2; model.pl=varargin{3};end;if nargin>3; model.il=varargin{4}; end

 % Define the case0
 Case=fe_case; CaseName='Case 1';
 if nargin>4; tr=varargin{5}; else tr=[]; end
 if nargin>5; adof=varargin{6}; else adof=[]; end
 if nargin>6; opt=varargin{7}; else opt=[0 0];end

 if length(opt)<3 
  if size(model.Elt,1)<5e3; opt(3)=0; else opt(3)=2; end
 elseif ~any([0:4]==opt(3));   error('not a valid value for opt(3)');
 end


 if ~isempty(tr) % projection basis
   Case.T=tr; model.DOF=adof; DOF0=adof;
   Case.DOF=-1-[1:size(tr,2)]'/1000;
 elseif ~isempty(adof)&&opt(2)==1 % impose DOFs
   model.DOF=adof; DOF0=adof;
 elseif ~isempty(adof)&&opt(2)~=1
   Case=fe_case(Case,'KeepDof','Retained DOFs',adof);DOF0=[];
 else; DOF0=[];
 end
 if ~isempty(Case.Stack)||~isempty(Case.T)
   CaseName='default fe_mk'; model=stack_set(model,'case',CaseName,Case);
 end
 % NOTE : one does not WANT to put rigid elements in the stack

else % this is the standard call starting with SDT 5.0 - - - - - - - - - - -

 opt=[0 0];
 [Case,CaseName]=fe_case(model,'getcase'); % deal with case selection
 ind=1:nargin;ind(1:carg-1)=0; 
 for j1=carg:nargin;
   if ischar(varargin{j1})&& comstr(comstr(comstr(varargin{j1},-27),1),'opt')
     opt=varargin{j1+1};ind(j1+[0 1])=0;
   end
 end
 if ~isempty(ind(ind~=0))
  Case=fe_case(Case,varargin{ind(ind~=0)});
 end
 if ~isfield(model,'DOF'); model.DOF=[]; end
 DOF0=model.DOF;

end % end of SDT 4.1/5.0 version calls
if isa(model,'v_handle'); model=model.GetData;end

% check of options

if length(opt)<3 % assembly method memory or disk
 if size(model.Elt,1)<20000; opt(3)=0; else opt(3)=2; end
elseif ~any([0:4]==opt(3));   error('not a valid value for opt(3)');
end
if length(opt)<4; opt(4)=0; 
elseif ~any([0 1 2 3]==opt(4)); error('not a valid value for opt(4)');
end
opt(1)=fix(opt(1)); opt=opt(:)'; st1='';

if ~isempty(model.Node); NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
else NNode=[];end

if opt(3)==4
  if ~isfield(model,'file')
   error('.file must be a field for type 3 superelement assembly');
  end
end
if nargout==3 && opt(1)~=0 
   disp('Warning argument count is not consistent with opt(1)');
elseif nargout==2 && opt(1)==0
  opt(1)=1;
end

% ----------------------------------------------------------------------
% definition of the groups

[EGroup,nGroup]=getegroup(model.Elt);
if nGroup == 0;		error('no element group specified in ELT'); end

% determination of DOFs
% nind: reindexing for nodal DOF, eind : reindexing for element DOFs


[Case,model.DOF,r1]=fe_mknl('initnocon',model,Case);
if isfield(Case,'InitFailed'); error('Init failed cannot assemble');end
nd=r1{1};nde=r1{2}; nw=r1{3}; eltid=r1{4}; node=r1{5};bas=r1{6};

if size(Case.T,1)~=length(model.DOF)||size(Case.T,2)~=length(Case.DOF)
   [Case,NNode]=fe_case(model,'gett',Case);
end
if ~isequal(node(:,1),model.Node(:,1)); % renumber with implicit nodes
    NNode=sparse(node(:,1),1,1:size(node,1));
end
% Memory allocation

st='';im=[];
% sparse full order matrix
 if opt(1,3)==1 
   sdtw('_nb','fe_mk : Symmetric real M&K  method is obsolete, using standard');
   opt(1,3)=0;
 elseif opt(1,3)==3  % old method by addition
   sdtw('_nb','fe_mk : Additive method is obsolete, using standard');opt(1,3)=0;
 end

 if opt(1,3)==0 % standard method in memory
    if nw(1)==0; nw(1)=2;end
    kie=zeros(nw(1),2);ik=int32([0 0]); 
    if opt(1,1)==0; mie=kie+0;im=int32([0 0]);else mie=[]; end; 
    
 elseif any(opt(1,3)==[2 4]) % standard method on disk
                           % standard method for type 3 superelement
    st = char(tempname,tempname);
    fidk = fopen(st(1,:),'w+');
    if opt(1,1)==0; fidm = fopen(st(2,:),'w+');
    else st=st(1,:); end
    nw(2:3)=[0 0];
    ik=0;im=0; kie=[]; mie=[];
    %mind=ones(size(model.Elt,1),1)*[1 0 1 0]; mind(EGroup(1:end-1),:)=0;
 end

n0=0; % shift for type 3 superelement with sub-matrix

% these are standard variable names as stated in elem0
pl=fe_mat('getpl',model); il=fe_mat('getil',model);
elt=model.Elt;N=length(model.DOF);
mdof=model.DOF; mind=zeros(size(elt,1),4);
def=[];

% ----------------------------------------------------------------------------
% loop over element groups
for jGroup = 1:nGroup

[ElemF,i1]= feutil('getelemf',model.Elt(EGroup(jGroup),:),jGroup);
EGID=i1(1); cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;

if EGID <0
 fprintf('\nGroup %i (EDID %i) is for display only and was ignored', ...
  jGroup,EGID);fp=1;
 mind(n0+cEGI,[1 3])=1;
else  

  DofPos=Case.GroupInfo{jGroup,1};
  pointers=Case.GroupInfo{jGroup,2};
  if opt(1)==0&&size(pointers,2)>1; pointers(2,:)=pointers(1,:);end
  pointers(5,:)=opt(1);
  integ=Case.GroupInfo{jGroup,3};
  constit=Case.GroupInfo{jGroup,4};
  gstate=Case.GroupInfo{jGroup,5};
  elmap=Case.GroupInfo{jGroup,6};
  InfoAtNode=Case.GroupInfo{jGroup,7};
  EltConst=Case.GroupInfo{jGroup,8};
  inode=fe_super('node',ElemF);

  step='call';
  try
      NodePos=int32(reshape(full(NNode(elt(cEGI,inode)')),length(inode),length(cEGI)));
      [eCall,SymFlag]=feval(ElemF,'matcall',integ,constit);SEopt=2;
      if strncmp(eCall,'mat_',4)||~isempty(findstr(eCall,'k1'))
      else;step='need';
      end
  catch;step='failed'; 
  end
  if ~strcmp(step,'call');
      %[eCall,SEopt]=fe_super('call',ElemF,opt); % V2006
      [eCall,SEopt]=fe_super('call',ElemF,model,opt);
      %SymFlag=1; fHandle=feval(ElemF,'call');
  end

  if isempty(eCall) 
     mind(cEGI,[1 3])=1;
     cEGI=[];
     if ~comstr(ElemF,'rigid')
     disp([ElemF ' : element function/superelement not found and ignored' 7]);
     end
  elseif comstr(eCall,ElemF)
    disp(eCall);cEGI=[];error('This call is not valid');
   elseif SEopt(1,1)==1; cEGI=EGroup(jGroup);n0=n0+1; % single superelement
  end
  t=cputime;

  %[fHandle,SymFlag]=feval(ElemF,'matcall');
  m1=[];k1=[];  % inits needed by some elements
  if strcmp(eCall,'mat_og'); 
   pointers=int32(pointers);NodePos=int32(NodePos);
  end

for jElt=1:length(cEGI)

 switch eCall
 case 'mat_of'
    nodeE=node(NodePos(:,jElt),[5:7 1]);
    point=pointers(:,jElt);
    if point(5)==0
      [k1,m1]=of_mk(ElemF,int32(point),integ,constit,nodeE);
      k1=reshape(k1(elmap),size(elmap,1),size(elmap,2));
      m1=reshape(m1(elmap),size(elmap,1),size(elmap,2));
    else 
      k1=of_mk(ElemF,int32(point),integ,constit,nodeE);
      k1=reshape(k1(elmap),size(elmap,1),size(elmap,2));
      m1=[];
    end
    if k1(1)<0      
      sdtw('_nb',sprintf(['element with index %i has negative orientation\n' ...
      'Use model.Elt=feutil(''orient'',model) before calling fe_mk\n' ...
      'fe_mk pursuing with ke=0 for now'], ...
      cEGI(jElt)));k1=k1*0; 
    end
 case 'mat_og'
    [k1,m1]=of_mk('matrixintegration',jElt,NodePos,Case.Node, ...
       pointers,integ,constit,gstate, ...
       elmap,InfoAtNode,EltConst); 
    if isfield(EltConst,'me')&&nnz(m1)==0; m1=EltConst.me;end
    if ~isempty(elmap);k1=k1(elmap);if ~isempty(m1);  m1=m1(elmap);end;end
    if k1(1)<0      
      sdtw('_nb',sprintf(['element with index %i has negative orientation\n' ...
      'Use model.Elt=feutil(''orient'',model) before calling fe_mk\n' ...
      'fe_mk pursuing with ke=0 for now'], ...
      cEGI(jElt)));k1=k1*0; 
    end
 otherwise
  if SEopt(1)==2  
      i2=elt(cEGI(jElt),inode);
      nodeE=node(NNode(i2(i2~=0)),[5:7 1]);
  elseif ~isempty(NodePos);nodeE=node(NodePos(:,jElt),[5:7 1]);
  end
  try; eval([eCall ';']);
  catch
    if ~isempty(strfind(lasterr,'orient'))
      error(['Some elements have negative orientation\n' ...
      'Use model.Elt=feutil(''orient'',model) before calling fe_mk']);
    elseif sp_util('diag')>=10;eval(eCall);
    else; fprintf('While evaluating\n%s\n',eCall);rethrow(lasterror);
    end
  end
 end %Mat_of or not

 if ~isreal(m1)||~isreal(k1)
        error('Assembly does not support complex matrices');
 end

  if jElt>length(cEGI)&&isempty(i1); break;
  elseif ~isempty(k1)

   if jElt>length(cEGI);   in1=int32(full(nd(round(i1(:)*100)-100))-1);
   else;                   in1=DofPos(:,jElt);
   end

   if opt(1,3)==0 % standard method in memory
     sp_util('insertinkie',k1,in1,kie,ik,N);
     if opt(1,1)==0; sp_util('insertinkie',m1,in1,mie,im,N); end
   elseif opt(1,3)==2||opt(1,3)==4 % standard on disk

    if opt(1,3)==4; 
      if length(unique(NodePos(:,jElt)))<size(NodePos,1) % element is degenerate
        [in1,in2,in3]=unique(in1);
        [II,JJ,KK]=find(m1);m1=sparse(in3(II),in3(JJ),KK);
        [II,JJ,KK]=find(k1);k1=sparse(in3(II),in3(JJ),KK);
      end
      k1=triu(k1); if opt(1,1)==0;  m1=triu(m1);end
    end
    if issparse(k1) && size(kie,1)<nnz(m1)
     if opt(1)==0; kie=zeros(max(nnz(m1),nnz(k1)),2);
     else kie=zeros(nnz(k1),2);end
    else kie=zeros(numel(k1),2);
    end

    ic=int32([0 0]);sp_util('insertinkie',k1,in1,kie,ic,N);
    fwrite(fidk,kie(1:double(ic(1)),:)','float64');

    mind(n0+cEGI(jElt),3:4)=ftell(fidk)/16+[-double(ic(1))+1 0];
    if opt(1)==0
     ic=int32([0 0]);sp_util('insertinkie',m1,in1,kie,ic,N);
     fwrite(fidm,kie(1:double(ic(1)),:)','float64');
     mind(n0+cEGI(jElt),1:2)=ftell(fidm)/16+[-double(ic(1))+1 0];
    end

   end % type of assembly
   if jElt>length(cEGI); break; end
  end % of there are DOFs for the element
  if jElt>length(cEGI); break; 
  elseif cputime-t>5 && rem(cEGI(jElt),100)==0 
   st1 = comstr(st1,-7,sprintf('Done %i elements',cEGI(jElt)));fp=1;t=cputime;
   drawnow;
  end

end % loop on elements of group
end % of real or display group 
end % loop on element groups


% ---------------------------------------------------------------------------
% Do some cleaning up at the end

%assembly of symmetric mass and stiffness matrices
clear NNode eind nind 
if length(ik)>1;ik=ik(1);end;if length(im)>1;im=im(1);end
if ~isreal(kie)||~isreal(mie)
         error('Cannot handle complex matrices are no longer supported');
elseif opt(1,3)==0 % mass and stiffness in memory using a single index

 if ~sp_util('issdt')
  k=sp_util('mind',kie(1:double(ik),1),kie(1:double(ik),2),N);
  if opt(1,1)==0  
   m=sp_util('mind',mie(1:double(im),1),mie(1:double(im),2),N);
  end
 else
  if opt(1,1)==0; m=sp_util('mind',mie,im,N); end
  k=sp_util('mind',kie,ik,N);
 end
elseif opt(1,3)==2
  
 fseek(fidk,0,-1);kie=fread(fidk,'float64');
 ik=length(kie)/2;kie=reshape(kie,2,ik)';
 k = sp_util('mind',kie,int32(ik),N);
 fclose(fidk);

 if opt(1,1)==0
   fseek(fidm,0,-1);kie=fread(fidm,'float64');
   ik=length(kie)/2;kie=reshape(kie,2,ik)';
   m = sp_util('mind',kie,int32(ik),N); fclose(fidm);
 end
 for j1=1:size(st,1); delete(st(j1,:));end

elseif opt(1,3)==4 % Saving the type 3 superelement

  model.Opt=[3;0;0];model.mind = mind;
  upcom(model,'save clear');
  ik=ftell(fidk)/16;
  fseek(fidk,0,-1); ki=fread(fidk,ik,'float64',8);
  fseek(fidk,8,-1); ke=fread(fidk,ik,'float64',8);
  fclose(fidk);

  st1={'ki','ke','DOF'};
  if opt(1,1)==0
   im=ftell(fidm)/16;
   fseek(fidm,0,-1); mi=fread(fidm,im,'float64',8);
   fseek(fidm,8,-1); me=fread(fidm,im,'float64',8);
   st1(end+[1:2])={'mi','me'};
   fclose(fidm);
  end

  if any(opt(2)==[1 2])  % assigned DOF positions

    [i1,i2]=sort(round([model.DOF;DOF0]*1000)); iout=zeros(size(model.DOF));
    i3=find(~diff(i1)); iout(i2(i3))=i2(i3+1)-length(model.DOF);

   N=length(model.DOF);ki=sort([iout(remi(ki,N))  iout(fix((ki-1)/N)+1)],2);
   ki=ki(:,1)+(ki(:,2)-1)*length(DOF0);
   if any(iout==0); i3=find(ki); ki=ki(i3);ke=ke(i3); 
    error('Problem with mind xxx');
   end
   if any(strcmp('mi',st1))
    N=length(model.DOF);mi=sort([iout(remi(mi,N))  iout(fix((mi-1)/N)+1)],2);
    mi=mi(:,1)+(mi(:,2)-1)*length(DOF0);
    if any(iout==0); i3=find(mi); mi=mi(i3);me=me(i3); end
   end

   DOF=DOF0;model.DOF=DOF0;

  else 
   try;
    r1=sort([remi(ki,N) fix((ki-1)/N)+1],2); ki=r1(:,1)+(r1(:,2)-1)*N;
    if any(strcmp('mi',st1))
      r1=sort([remi(mi,N) fix((mi-1)/N)+1],2); mi=r1(:,1)+(r1(:,2)-1)*N;
    end
   catch; % if failure due to memory use blocks
    in1=1:1e4;
    while in1(length(in1))<length(ki)
      r1=sort([remi(ki(in1),N) fix((ki(in1)-1)/N)+1],2);
      ki(in1)=r1(:,1)+(r1(:,2)-1)*N; in1=in1+in1(length(in1));
    end
    if any(strcmp('mi',st1))
     in1=1:1e4;
     while in1(length(in1))<length(mi)
      r1=sort([remi(mi(in1),N) fix((mi(in1)-1)/N)+1],2);
      mi(in1)=r1(:,1)+(r1(:,2)-1)*N; in1=in1+in1(length(in1));
     end
    end
   end % catch
    DOF=model.DOF;
  end

  save(fullfile(model.wd,model.file),st1{:},'-append');
  clear ke ki

  for j1=1:size(st,1); delete(st(j1,:));end
  if fp; fprintf(1,'\n'); end
  m=model;

  % With the new sparse solvers we no longer need profiling
  %disp('UPCOM profiling (this may take a while)');
  %m=upcom(m,'profile fix0');
  return;

end


% need or not to eliminate 0 DOFs and profile  - - - - - - - - - - - -


if ~isfield(model,'Opt'); model.Opt=[1]; end
if ~isfield(model,'K');   model.K={};    end
if opt(1)==0; model.K={m,k};model.Klab={'M','K'};model.Opt(2,1:2)=[2 1]; 
else model.K={k};model.Klab={fe_mk(sprintf('mattype %i',opt(1)))};
 model.Opt(2,:)=0;model.Opt(2,1)=opt(1);
end

if     opt(2)==2 % place DOF but ignore case information

 if ~isequal(model.DOF,DOF0)&&~isempty(DOF0)
  
   i1=fe_c(DOF0,model.DOF,'ind');T=speye(length(model.DOF));
   if length(i1)~=size(model.DOF)  % existing case DOFs not existing in DOF0
    i2=fe_c(model.DOF,DOF0,'ind');
    i1=fe_c(DOF0,model.DOF(i2),'ind');
    [II,JJ,TT]=find(T(:,i2));
   else [II,JJ,TT]=find(T);
   end
   T=sparse(II,i1(JJ),TT,size(T,1),length(DOF0));
   model.DOF=DOF0;  model.K=tkt_femk(T,model.K);
   i1=[];
 end
 mdof=model.DOF;

elseif opt(2)==1 % place DOF

 if ~isequal(model.DOF,DOF0)&&~isempty(DOF0)

   i1 = find(abs(rem(DOF0,1))<1e-6 | abs(fix(DOF0)-DOF0)<1e-6);
   if ~isempty(i1); error('you cannot use wild cards and place DOFs'); end
   i1=fe_c(DOF0,Case.DOF,'ind');
   if length(i1)~=size(Case.DOF)  % existing case DOFs not existing in DOF0
    i2=fe_c(Case.DOF,DOF0,'ind');
    i1=fe_c(DOF0,Case.DOF(i2),'ind');
    [II,JJ,TT]=find(Case.T(:,i2));
   else [II,JJ,TT]=find(Case.T);
   end
   T=sparse(II,i1(JJ),TT,size(Case.T,1),length(DOF0));
   model.K=tkt_femk(T,model.K);
   Case.DOF=DOF0; model.DOF=Case.DOF; Case.T=speye(length(Case.DOF));
   i1=[];

 end 
 mdof=model.DOF;

elseif ~isempty(Case.Stack)
    if opt(2)==1
     sdtw('_nb', ...
      'fe_mk: you cannot use opt(2)==1 and fix/keep dof at the same time')
    end
    model.K=tkt_femk(Case.T,model.K); mdof=Case.DOF; i1=[];
end

N=size(model.K{1},1);

if     opt(1,4)==0; i1=1:N; % nothing is done
elseif opt(1,4)==1&&comstr(version,'7'); i1=symamd(model.K{1}); 
elseif opt(1,4)&&sdtdef('verm')<7.7; i1=symmmd(model.K{1});   % symmmd sorting
elseif opt(1,4)==3; i1=1:N; % nothing is done
else % if opt(1,4)==2
   st=ofact('SymRenumber');
   switch  st
   case 'symrcm'; 
     if  exist('symrcm','builtin');i1 = symrcm(model.K{1});
     else;i1 = sparsfun('symrcm',model.K{1}); 
     end
   case 'symmmd'; 
    if comstr(version,'7');i1=symamd(model.K{1});else;i1=symmmd(model.K{1});end
   case '';       i1=1:N;
   otherwise;     eval('i1=%s(model.K{1});',st);
   end
end

% eliminate DOFs with zero stiffness contributions
if opt(1,2)==0&&any([0 1]==opt(1,1))
  i2=diag(model.K{1}~=0);if length(model.K)>1; i2=i2+diag(model.K{2}~=0);end
  i2=i2(i1);i1=i1(i2~=0);
elseif 1==2 % do the auto-spc here xxx
end

if ~isempty(i1)&&opt(1,2)~=2
  for j1=1:length(model.K);model.K{j1}=model.K{j1}(i1,i1);end
  mdof=mdof(i1); 
  if opt(1,2)==0;Case.T=Case.T(:,i1); Case.DOF=Case.DOF(i1);end
end

if nargout<2 % standard output starting with SDT 4.2 - - - - - - - - - - - -

  model=stack_set(model,'case',CaseName,Case);
  m=model;

% old format Standard output 
elseif nargout==2;  m=model.K{1}; k=mdof;
else

 m=model.K{1};k=model.K{2};

end

if fp; fprintf(1,'\n'); end

% ---------------------------------------------------------------------------
function K=tkt_femk(T,K)

if sp_util('issdt'); K=feutilb('tkt',T,K);
else
 Tt=T';
 for j1=1:length(K)
   %k=K{j1}+spalloc(size(K{j1},1),size(K{j1},1),0);
   [II,JJ,KK]=find(K{j1});k=sparse(II,JJ,KK,size(K{j1},1),size(K{j1},1));
   k=Tt*k; 
   k=k*T;
   K{j1}=k;
 end
end
