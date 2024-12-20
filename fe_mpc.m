function [out,out1]=fe_mpc(varargin)

%FE_MPC	generic handling of multiple point constraints and rigid links
%
% Syntax : [Case,mdof] = fe_mpc(model,Case,mdof);
%
%       Multiple point constraints are linear constraints between
%        DOFs of a model.
% 
% Utilities
%
% model=fe_mpc('fixrbe3Alt',model) % used to transform from alternate to
%  nominal RBE3 format
%  RBE3 nominal format 
%   Rbe3ID NodeIdSlave DofSlave Weight1 DofMaster1 NodeId1 Weight2 ...  
%  RBE3 alternate format
%   Rbe3ID NodeIdSlave DofSlave Weight DofMaster NodeId1 NodeId2 ...
%
%  model=fe_mpc('Rbe3Id',model) generates unique identifiers for mutiple RBE3
%  model=fe_mpc('DofSetMerge',mode,'name1','nam2') : merge DofSet entries
%
%   See sdtweb    eltfun, elem0
%	See also help rigid, celas, bar1, beam1, ...


%       Etienne Balmes, Guillaume Vermot des Roches
%       Copyright (c) 1990-2023 by SDTools, All Rights Reserved.

% Get the model, figure out true node locations, ...

%#ok<*NOSEM,*NASGU>

model=varargin{1};carg=2;
if ~ischar(model)
elseif comstr(varargin{1},'cvs')
 out='$Revision: 1.133 $  $Date: 2024/12/06 18:06:43 $'; return;
elseif comstr(lower(varargin{1}),'fixrbe3')
  %% #fixRBE3 ----------------------------------------------------------------
 r1=varargin{2};
 % Transform to nominal format 
 %  Rbe3ID NodeIdSlave DofSlave (Weight1 DofMaster1 NodeId1) ...
 % From alternate
 %  Rbe3ID NodeIdSlave DofSlave Weight NodeId1 NodeId2 ...
 if isfield(r1,'Elt')&&isfield(r1,'Stack')
  model=r1;[Case,CaseName]=fe_case(model,'getcase');
  for j1=find(strcmpi(Case.Stack(:,1),'rbe3'))'
    Case.Stack{j1,3}=fe_mpc(varargin{1},Case.Stack{j1,3},model);
  end
  out=stack_set(model,'case',CaseName,Case);
  return;
 end
 out=r1;
 RunOpt.Alternate=~isempty(strfind(lower(varargin{1}),'alt'));
 if isfield(r1,'MasterSel') % implicit defintion resolution
    model=varargin{3};
    if ischar(r1.MasterSel);i1=feutil(['findnode' r1.MasterSel],model);
    else; i1=r1.MasterSel; 
    end
    if ischar(r1.SlaveSel);i2=feutil(['findnode' r1.SlaveSel],model);
    else;i2=r1.SlaveSel;
    end
    if length(i2)>1; error('RBE3, SlaveSel %s should contain only one node',r1.SlaveSel); end
    if isempty(i1)||isempty(i2); out=[]; out1=RunOpt.Alternate; return; end
    r1=[[1;i2;r1.DOF] [[1;r1.MasterDOF]*ones(1,length(i1));i1(:)']];
    r1=r1(:)';
 elseif isfield(r1,'data');
   if isfield(r1,'Alternate');RunOpt.Alternate=r1.Alternate;end
   r1=r1.data;
 end
 r2=r1(:,8:3:end); r2=r2(:);r2(r2==0)=[]; % DofMaster
 r2=setdiff(unique(sprintf('%i',r2(:))),'123456');
 if rem(max(find(any(r1,1))),3)~=0; %#ok<MXFND> % Number of values not multiple of 3
    RunOpt.Alternate=1;
 elseif ~isempty(r2)&&~RunOpt.Alternate % xxx allow integer weights
  %any(rem(r1(:),1)) ... % Non integer weight
   %   ||~isempty(r2) % Cannot be a rbe3
         RunOpt.Alternate=1;      
 end
 if RunOpt.Alternate==-1 
   %% accept groups by weight/f(distance)
   data=r1;
   for j1=1:size(data,1)
    r1=data(j1,:);i1=find(r1,1,'last');r1=r1(1:i1);
    if all(r1(4:3:end)==1); continue;end % this is nominal format
    i2=find(r1<1); 
    if ~isempty(i2)
    elseif r1(4)==1;i2=4;% just one weight at 1
    else; error('Not expected');
    end 
    i3=cell(length(i2),1); i2(end+1)=length(r1)+1; 
    for j2=1:length(i2)-1;
      i4=r1(i2(j2)+2:i2(j2+1)-1)';
      i3{j2}=reshape([repmat(r1(i2(j2)+[0 1]),length(i4),1) i4]',1,[]);
    end
    r1=[r1(1:3) horzcat(i3{:})];data(j1,1:length(r1))=r1;
   end
   r1=data;
 elseif RunOpt.Alternate==1
   if varargin{1}(end)~=';' % Silent transform
     fprintf('Transforming %i to nominal RBE3 format (fe_mpc fixRbe3Alt)\n',...
         r1(:,1));
   end
   data=r1;
   for j1=1:size(data,1) % transform to nominal format
    r1=data(j1,:); i1=find(r1,1,'last');
    if any(r1(1:i1)==0); 
     sdtw('_nb','zero entries inside a rbe3 definition is not possible, they will be suppressed in sequence: \n%s',...
     sdtw('_clip 60 1','%g ',r1(1:i1)));
     r1(r1==0)=[]; i1=length(r1);
    end
    if RunOpt.Alternate;
     r2=r1(1:3)';r2(3,2:i1-4)=r1(6:i1);r2(2,2:end)=r1(5);r2(1,2:end)=r1(4);
     r1=r2(:)';data(j1,1:length(r1))=r1;
    end
   end
   r1=data;
 end
 if isfield(out,'data'); out.data=r1; 
 else; out=r1; 
 end
 out1=RunOpt.Alternate; return;

elseif comstr(lower(varargin{1}),'rbe3c')
%% #RBE3C data=fe_mpc('rbe3c',model,CaseEntryName) ---------------------------
 model=varargin{carg};carg=carg+1;
 name=varargin{carg};carg=carg+1;
 data=fe_case(model,'getdata',name);
 data=fe_mpc('fixrbe3',data,model);
 if isstruct(data); data=data.data; end
 i1=data(:,[2 6:3:end]);i1=unique(i1(:));i1(~i1)=[];
 mo1=struct('Node',feutil('getnode',model,i1),'Elt', ...
     feutil('objectmass',i1));
 mo1.Elt(1,1:6)=[Inf abs('celas')];mo1.Elt(2:end,3)=123456;
 if isfield(model,'bas')&&~isempty(model.bas); mo1.bas=model.bas; end
 mo1=fe_case(mo1,'rbe3',name,data);mo1.DOF=feutil('getdof',mo1);
 Case=fe_case(mo1,'gett');
 i1=fe_c(mo1.DOF,data(:,2),'ind');      % Slave DOFs
 i2=fe_c(Case.DOF,mo1.DOF(i1),'ind',2); % Master DOFs
 r1=Case.T(i1,i2); 
 i3=any(r1,2);
 out=struct('c',[speye(nnz(i3)) -r1(i3,:)], ...
     'DOF',[mo1.DOF(i1(i3));Case.DOF(i2)],'slave',1:nnz(i3));
 i1=~any(out.c,1); % remove spurious DOF that can be linked to 123 links instead of 123456
 if any(i1); out.c(:,i1)=[]; out.DOF(i1)=[]; end
 out.Stack={'info','origtyp','rbe3'};
 return
 
elseif comstr(lower(varargin{1}),'rigidc')
%% #RigidC data=fe_mpc('rigidc',model,CaseEntryName) xxx should call rigid
 % recover constraint from Case
 model=varargin{carg};carg=carg+1;
 name=varargin{carg};carg=carg+1;
 data=fe_case(model,'getdata',name);
 if isstruct(data); elt=data.Elt; else; elt=data; end
 if ~any(~isfinite(elt(:,1))); elt=feutil('addelt',[],'rigid',elt); end
 i1=isfinite(elt(:,1)); i2=rem(elt(i1,3),1)~=0; i1=find(i1);
 if any(i2); elt(i1(i2),3)=round(elt(i1(i2),3)*100); end % xxx old openfem format?
 n1=feutil('getnode',model,unique(reshape(elt(isfinite(elt(:,1)),1:2),[],1)));
 mo1=struct('Node',n1,'Elt',feutil('objectmass',n1(:,1)));
 mo1=fe_case(mo1,'stack_set',{'rigid','rigid',elt});
 if isfield(model,'bas')&&~isempty(model.bas); mo1.bas=model.bas; end
 mo1.DOF=feutil('getdof',mo1); Case=fe_case(mo1,'gett');
 i1=fe_c(mo1.DOF,Case.DOF,'ind',2); 
 c=Case.T'; c=c(:,i1)'; c=[-c speye(length(i1))];
 out=struct('c',c,'DOF',[Case.DOF;mo1.DOF(i1)],...
  'slave',length(Case.DOF)+[1:length(i1)]','Stack',{{'info','origtyp','rigid'}});
 return
 
elseif comstr(lower(varargin{1}),'mpcc')
 %% #MPCC data=fe_mpc('mpcc',model,CaseEntryName)
 % data is already here, only pack it up and clean
 model=varargin{carg};carg=carg+1;
 name=varargin{carg};carg=carg+1;
 data=fe_case(model,'getdata',name);
 out=data;
 try
  if ~isfield(data,'c')||isempty(data.c) % resolve if implicit
   [un1,data]=fe_caseg('connection',model,name,data);
  end
  if ~isfield(data,'slave');
   [data.c,data.slave]=feutil('fixmpcmaster',data.c);
 end
 catch err
  sdtw('_nb','failed to cleanup MPC %s, with:',name);
  err.getReport;
 end
 out=feutil('rmfield',data,setdiff(fieldnames(data),...
  {'c','DOF','slave','wjdet','Stack'})); 
 out=stack_set(out,'info','origtyp','mpc');
return
 
elseif comstr(lower(varargin{1}),'fixdofc')
 %% #FixDofC data=fe_mpc('fixdofc',model,CaseEntryName)
 % recover constraint from Case
 model=varargin{carg};carg=carg+1;
 name=varargin{carg};carg=carg+1;
 data=fe_case(model,'getdata',name);
 if 1==1 % resolve DOF to be constrained
  dof=feval(fe_case('@getFixDof'),model,[],fe_case(model,'stack_get','fixdof',name));
  if ~all(rem(dof,1))||any(fix(dof)==0) % refine if DOF are not explicitly defined (nodes)
   adof=feutil('getdof',model);
   dof=fe_c(adof,dof,'dof',1);
  end
 else % old not robust nor efficient
  dof=feutil('getdof',model);
  mo1=struct('Node',model.Node,'Elt',feutil('objectmass',model.Node(:,1)));
  if isfield(model,'Stack'); mo1.Stack=model.Stack; end
  mo1=fe_case(mo1,'reset');
  mo1=fe_case(mo1,'stack_set',{'FixDof','FixDof',data});
  if isfield(model,'bas')&&~isempty(model.bas); mo1.bas=model.bas; end
  mo1.DOF=dof; Case=fe_case(mo1,'gett');
  i1=fe_c(dof,Case.DOF,'ind',2);
  dof=dof(i1);
 end
 out=struct('c',speye(length(dof)),'DOF',dof,'Stack',...
  {{'info','sel',data;'info','origtyp','FixDof'}});
 return

else % other char commands
 [CAM,Cam]=comstr(varargin{1},1);
 %% #MPCMerge ----------------------------------------------------------------
 if comstr(Cam,'mpcmerge')     % combine two MPC
   data1=varargin{carg}; carg=carg+1;
   if iscell(data1) 
    %% Merge CE cards from Ansys at the end (fast merge many)
    DOF=cellfun(@(x)round(x.DOF(1)*100),data1);[~,i3]=unique(DOF);
    if length(i3)~=length(data1); 
        sdtw('_nb','MpcMerge: removing repeated contraints');
    end
    data1=data1(i3); % Skip repeated

    DOF=cellfun(@(x)round(x.DOF*100),data1,'uni',0);
    [DOF,i1,i2]=unique(vertcat(DOF{:}));DOF=DOF(:)/100;
    r1=cellfun(@(x)x.c,data1,'uni',0);j0=1; 
    for j1=1:size(r1,2);
        r1{2,j1}=ones(size(r1{1,j1}))*j1;
        r1{3,j1}=data1{j1}.DOF(1);j0=j0+length(r1{2,j1});
    end
    out=struct('c',sparse(horzcat(r1{2,:}),i2,horzcat(r1{1,:})), ...
        'DOF',DOF,'slave',fe_c(DOF,horzcat(r1{3,:})','ind'));
    return
   end
   data2=varargin{carg}; carg=carg+1;
   if isempty(data1); out=data2; return;
   elseif isempty(data2); out=data1; return; 
   end
   out=data1;
   out.DOF=[out.DOF;data2.DOF(~ismember(round(data2.DOF*1000),round(data1.DOF*1000)))];
   r1=feutil('placeindof',out.DOF,data1); out.c=r1.c;
   r1=feutil('placeindof',out.DOF,data2); out.c=[out.c;r1.c];
   if isfield(data1,'slave')&&~isfield(data2,'slave')
    [data2.c,data2.slave]=feutil('fixMpcMaster',data2.c);
   elseif ~isfield(data1,'slave')&&isfield(data2,'slave')
    [data1.c,data1.slave]=feutil('fixMpcMaster',data1.c);
   end
   if isfield(data1,'slave'); 
       out.slave=[data1.slave(:);fe_c(out.DOF,data2.DOF(data2.slave),'ind')]; 
   end
   if isfield(data2,'ID')
    if isfield(out,'ID'); out.ID=cat(1,out.ID,data2.ID);
    else;out.ID=data2.ID;
    end
   end
   if isfield(out,'Stack'); out.Stack=[]; end % Remove info,sel
 %% #Rbe3 utilities ----------------------------------------------------------
 elseif comstr(Cam,'rbe3');[CAM,Cam]=comstr(CAM,5);  
  if comstr(Cam,'id')
   %% #RBE3Id: fix rbe3 to unique IDs
   model=varargin{carg};carg=carg+1; 
   r1=fe_case(model,'stack_get','rbe3');
   i1=0; 
   for j1=1:size(r1,1)
     r2=r1{j1,3}; if isstruct(r2);r2=r2.data;end
     r2(:,1)=i1+(1:size(r2,1))';i1=r2(end,1);
     if isstruct(r1{j1,3});r1{j1,3}.data=r2;else;r1{j1,3}=r2;end
   end
   out=fe_case(model,'stack_set',r1);
   
  elseif comstr(Cam,'split')
   %% #RBE3Split: generate different case entries for each lines of RBE3
   model=varargin{carg};carg=carg+1; 
   if carg<=nargin; li=varargin{carg}; carg=carg+1; else; li=''; end
   r1=fe_case(model,'stack_get','rbe3',li);
   CaseStack=fe_case(model,'stack_get');
   for j1=1:size(r1,1)
    model=fe_case(model,'stack_rm',r1{j1,1},r1{j1,2});
    r2=r1{j1,3}; if isstruct(r2); r2=r2.data; end
    for j2=1:size(r2,1)
     r3=r2(j2,:); i2=find(r3==0,1,'first');
     if ~isempty(i2); r3=r3(1:i2-1); end
     model=fe_case(model,'stack_set',r1{j1,1},sprintf('%s_%i',r1{j1,2},j2),r3);
    end
   end
   
   out=model;   
   
  elseif comstr(Cam,'masterdofclean')
   %% #RBE3MasterDOFClean: edit rotation DOF if master node is only on volume elements
   model=varargin{carg};carg=carg+1;
   if carg<=nargin; li=varargin{carg}; carg=carg+1; else; li=''; end
   if carg<=nargin; RE=varargin{carg}; carg=carg+1; else; RE=[]; end % EltNodeCon (issdt)
   if ~isempty(RE)&&sp_util('issdt'); % recover groups and mpid to get a direct ElemF check
    NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
    RE=feval(feutilb('@EltNodeCon'),model.Elt,NNode,struct('useAllNodes',1));
    RE.cNodeElt=RE.cNodeElt';
    RE.VolTypes=cellfun(@(x)fe_mat('type',['m_' x],'SI',1),...
     {'hexa8','hexa20','hexa27','tetra4','tetra10','penta6','penta15','pyra5','pyra13'}','uni',1);
   end
   r1=fe_case(model,'stack_get','rbe3',li);
   %CaseStack=fe_case(model,'stack_get');
   for j1=1:size(r1,1)
    r2=r1{j1,3};
    r2=fe_mpc('fixrbe3',r2,model);
    n1=r2(6:3:end);
    %% check whether masterdof 456 exist and remove if not
    if isempty(RE) % high level calls (can be sub optimal with large numbers of ops)
     el1=feutil(sprintf('selelt withnode{nodeid %s}',num2str(n1(:)')),model);
     if isempty(el1);  sdtw('_nb','%s is a RBE3 involving only spurious nodes', r1{j1,2});
     else
      [EGroup,nGroup]=getegroup(el1); i0=1;
      for jGroup=1:nGroup
       ElemF=getegroup(el1(EGroup(jGroup),:),jGroup); %model.Elt(EGroup(jGroup),:),jGroup);
       if strncmpi(ElemF,'hexa',4)||strncmpi(ElemF,'tetra',5)||...
         strncmpi(ElemF,'penta',5)||strncmpi(ElemF,'pyra',4)
       else; i0=0; break;
       end
      end
     end
    else % direct check using cNodeElt and ElemCodes
     i0=all(ismember(RE.ElemCode(any(RE.cNodeElt(:,n1),2)),RE.VolTypes));
    end
    if i0;r2(5:3:end)=123; end
    r1{j1,3}=r2;
   end
   model=fe_case(model,'stack_set',r1);
   out=model;
   
  elseif comstr(Cam,'torbe2')
   %% #Rbe3ToRbe2: transform an RBE3 connection to an RBE2 version
   % question to multiple slave handling, less critical in this way
   model=varargin{carg};carg=carg+1; 
   if carg<=nargin; li=varargin{carg}; carg=carg+1; else; li=''; end
   r1=fe_case(model,'stack_get','rbe3',li);
   
   for j1=1:size(r1,1)
    model=fe_case(model,'stack_rm',r1{j1,1},r1{j1,2});
    r2=r1{j1,3}; 
    if isstruct(r2); 
     if isfield(r2,'data');r2=r2.data; 
     else;r2=fe_mpc('fixRbe3;',r2,model); elt=[];
     end
    end
    r2=fe_mpc('fixRbe3;',r2); elt=[];
    for j2=1:size(r2,1)
     r3=r2(j2,:); i3=find(r3==0,1,'first'); if ~isempty(i3); r3=r3(1:i3-1); end
     r4=[r3(6:3:end)' r3(5:3:end)'];
     r4=[r3(2)+zeros(size(r4,1),1) r4];
     elt=[elt;r4];
    end
    elt=feutil('addelt',[],'rigid',elt);
    model=fe_case(model,'rigid',sprintf('%s_rbe2',r1{j1,2}),elt);
   end
   out=model;
   
   
  end % END RBE3
   
%% #CleanUsed : remove constraints that use nodes no longer present ----------
elseif comstr(Cam,'cleanused')  
     
mo1=varargin{carg};carg=carg+1;
if carg<=nargin&&size(varargin{carg},2)==7;n1=varargin{carg};carg=carg+1;
else; n1=feutil('optimmodel',mo1);
end
[Case,CaseName]=fe_case(mo1,'getcase');
out1=struct('RemovedNodes',[],'Master',[],'slave',[]);
for j1=1:size(Case.Stack,1)
 r1=Case.Stack{j1,3};
 switch lower(Case.Stack{j1,1})
 case 'rigid'
  if ~isempty(r1)
   i2=[false;~ismember(r1(2:end,1),n1(:,1))|~ismember(r1(2:end,2),n1(:,1))];
   if nargout>1; 
    r2=r1(~i2,:);r2(1,:)=[];
    out1.RemovedNodes=[out1.RemovedNodes;unique(r2(:,1:2))];
   end
   r1(i2,:)=[];
  end
  if size(r1,1)==1;Case.Stack{j1,3}=[];else;Case.Stack{j1,3}=r1;end
 case 'mpc'
  i2=fe_c(r1.DOF,n1(:,1),'ind',2);i3=any(r1.c(:,i2),2);
  if nargout>1; 
   r2=r1.c(~i3,:);r2(:,i2)=[];r3=r1.DOF;r3(i2)=[];
   out1.RemovedNodes=[out1.RemovedNodes;unique(fix(r3(any(r2,1))))];
  end
  r1.c(i3,:)=[]; r1.c(:,i2)=[];r1.DOF(i2)=[]; r1=feutil('rmfield',r1,'slave');
  if isempty(r1.c);Case.Stack{j1,3}=[];else;Case.Stack{j1,3}=r1;end
 case 'rbe3'
  r1=fe_mpc('fixrbe3',r1,model);
  if isfield(r1,'data')
   i2=any(reshape(ismember(r1.data(:,6:3:end),n1(:,1)),size(r1.data,1),[]),2);
   if nargout>1
    out1.RemovedNodes=[out1.RemovedNodes;setdiff(r1.data(~i2,6:3:end),0)];
   end
   r1.data=r1.data(i2,:);
  else
   i2=any(reshape(ismember(r1(:,6:3:end),n1(:,1)),size(r1,1),[]),2);
   if nargout>1
    out1.RemovedNodes=[out1.RemovedNodes;setdiff(r1(~i2,6:3:end),0)];
   end
   r1=r1(i2,:);
  end
  Case.Stack{j1,3}=r1;
 case 'fixdof' % Actually no need to remove
 case 'dofload' 
   i1=fe_c(r1.DOF,n1(:,1),'ind');
   if isempty(i1); Case.Stack{j1,3}=[];end
 otherwise
    fprintf('CleanUsed %s,%s not implemented\n',Case.Stack{j1,1:2});
 end
end
Case.Stack(cellfun('isempty',Case.Stack(:,3)),:)=[];
mo1=stack_set(mo1,'case',CaseName,Case);
out=mo1;

%% #FixDofBas2mpc : transform local basis fixdof to mpc ----------------------
elseif comstr(Cam,'fixdofbas2mpc')  
 mo1=varargin{carg};carg=carg+1;
 name=varargin{carg};carg=carg+1;
 r1=fe_case(mo1,'getdata',name);
 if isfield(r1,'data');r1=r1.data; end; i1=fix(r1);
 DOF=feutil('getdof',fix(r1),unique(round(rem(r1,1)*100))/100);
 cGL = basis('trans e',mo1.bas,mo1.Node,DOF);
 data=struct('c',cGL(:,fe_c(DOF,r1,'ind'))','DOF',DOF);
 [data.c,data.slave]=feutil('fixMpcMaster',data.c);
 mo1=fe_case(mo1,'remove',name,'mpc',name,data);
 out=mo1; 
 
 %% #DofSetMerge : combine multiple DOFsets into 1 ---------------------------
 elseif comstr(Cam,'dofsetmerge')  
     
  model=varargin{carg}; carg=carg+1;
  r1={};
  [Case,CaseName]=fe_case(model,'getcase');
  while carg<=nargin
   [r2,i1]=stack_get(Case,'DofSet',varargin{carg});carg=carg+1;
   Case.Stack(i1,:)=[];
   r1=[r1;r2]; 
  end
  out=r1{1,3}; 
  if isfield(out,'name')&&~isfield(out,'lab');out.lab={out.name};end
  for j1=2:size(r1,1)
    d2=r1{j1,3};      
    DOF=[out.DOF;d2.DOF(~ismember(round(d2.DOF*1000),round(out.DOF*1000)))];
    out=feutil('placeindof',DOF,out);
    d2=feutil('placeindof',DOF,d2);
    if isfield(d2,'name')&&~isfield(d2,'lab');d2.lab={d2.name};end
    out=fe_def('appenddef',out,d2);
  end
  if ~isfield(out,'name');out.name='Merge';end
  Case=stack_set(Case,'DofSet',out.name,out);
  model=stack_set(model,'case',CaseName,Case);
  out=model;
  
 %% #MassAtMaster ------------------------------------------------------------
 elseif comstr(Cam,'massatmaster')  
  model=varargin{carg}; carg=carg+1;
  Case=fe_case(model,'getcase');
  r1=stack_get(Case,'mpc');
  DOF=feutil('getdof',model);r3=cell(size(r1,1),1);
  for j1=1:size(r1,1)
   r2=r1{j1,3};r4=r2.DOF;
   if ~isfield(r2,'slave'); 
     fprintf('Ignoring MPC %s with no slave\n',r1{j1,2});
   else
    r4(r2.slave)=[];
    r3{j1}=fe_c(r4,DOF,'dof',2);
   end
  end
  r1=stack_get(Case,'dofset');
  for j1=1:size(r1,1)
   r2=r1{j1,3};r4=r2.DOF;r3{end+1}=fe_c(r4,DOF,'dof',2); %#ok<AGROW>
  end

  r3=vertcat(r3{:});
  r3=unique([fix(r3) round(rem(r3,1)*100)],'rows');
  if ~isempty(r3)
   [r2,i1]=unique(r3(:,1));for j1=1:length(i1);r2(j1,r3(i1(j1),2)+1)=1e-10;end 
   model=feutil('addelt',model,'mass1',r2);
  end
  % model=stack_set(model,'case',CaseName,Case);
  out=model;
  
 %% #End other commands
 else;error('%s unknown',CAM);
 end
 return
end
%% #Resolutions --------------------------------------------------------------
% Get the Case
Case=[]; if carg<=nargin;Case=varargin{carg};carg=carg+1; end
if isempty(Case)||~isfield(Case,'Stack')
 eval('Case=fe_case(model,''getcase'');')
end
if isempty(Case)||~isfield(Case,'Stack')
 Case=struct('Stack',[]); Case.Stack=cell(0,3);
end


if ~isfield(model,'Node'); error('model must contain nodes'); end  
if isfield(Case,'Node')&&isfield(Case,'cGL') % coming from gett NodeBas done
  node=Case.Node;
  if ~isequal(model.DOF,Case.DOF) % avoid rotating unused DOFs
      node(~ismember(node(:,1),fix(model.DOF))&node(:,3),:)=[];
  end
  [cGL,mech_dof]=basis('trans l',Case.bas,node);  
elseif isfield(model,'bas')&&~isempty(model.bas)
  [node,bas]=basis('nodebas',model);
  [cGL,mech_dof]=basis('trans l',bas,node);
else 
  i1=(1:6)/100; 
  if isfield(Case,'Node');node=Case.Node;else; node=model.Node;end
  mech_dof=i1(ones(size(node,1),1),:)+node(:,ones(6,1));
  mech_dof=mech_dof';mech_dof=mech_dof(:);cGL=[];
end


% Get the active DOFs if any and append any non mech_dof to mech_dof

if carg<=nargin; mdof=varargin{carg};carg=carg+1; 
elseif isfield(model,'DOF')&&~isempty(model.DOF)
 mdof=model.DOF;
else
 mdof=feutil('getdof',model);
end
if carg<=nargin; RunOpt=varargin{carg};carg=carg+1; 
else; RunOpt=struct('fixdof',[]);
end

RunOpt.NonMechInd=fe_c(mdof,[13:99]'/100,'ind');
if ~isempty(RunOpt.NonMechInd) 
 if ~isempty(cGL); % non mech DOF are not transformed yet
  i1=1:length(RunOpt.NonMechInd);
  cGL=[cGL spalloc(size(cGL,1),length(RunOpt.NonMechInd),0); 
       sparse(i1,length(mech_dof)+i1,1)];
 end
 mech_dof=[mech_dof;mdof(RunOpt.NonMechInd)]; 
end

% Allow rigid in Elt rather than case
if isfield(model,'Elt') 
   [i2,el1]=feutil('findelt egid>=0 & eltname rigid ',model);
   if ~isempty(el1)
    Case=stack_rm(Case,'rigid','rigid in .Elt');
    Case.Stack(end+1,1:3)={'rigid','rigid in .Elt',el1};
   end
else;el1=[];
end

% join all keepdof entries
if ~isempty(Case.Stack); ind=find(strcmpi('keepdof',Case.Stack(:,1)));
else;ind=[];end
if ~isempty(ind)
   i1=[];
   for j1=ind(:)'
         r1=Case.Stack{j1,3};  
         if isnumeric(r1) 
         elseif isfield(r1,'data'); r1=r1.data;
         elseif ischar(r1);r1=feutil(['findnode' r1(:)'],model);
         else;error('Not supported KeepDof format');   
         end
         i1=[i1;r1(:)];
   end
   ikeep=fe_c(mdof,i1,'ind');
else; ikeep=[]; 
end

NNode=sparse(node(:,1),1,1:size(node,1));

if ~isempty(mech_dof)
 N=length(mech_dof);
 if ~isempty(Case.Stack) % avoid if unused
  nd=feval(fe_mknl('@getPosFromNd'),[],mech_dof, ...
   max([fix(max(mech_dof))+1 max(node(:,1))]));
  slave=feval(nd.setDofListNd,nd);
 else; slave=[];
 end
else;N=0;slave=[];
end
II=[]; JJ=[]; TT=[];  % contains rows columns and values of constraint matrix


for j0=1:size(Case.Stack,1)

switch lower(Case.Stack{j0,1})

case 'rigid' % #deal_with_rigid_element - - - - - - - - - - - - - - - - - - - 

% Increment slave DOF list

elt=Case.Stack{j0,3};if isfield(elt,'Elt');elt=elt.Elt;end
cEGI=find(isfinite(elt(:,1)));

i1=abs(num2str(elt(cEGI,3)))-48;
i1=(i1/100+elt(cEGI,2*ones(size(i1,2),1))).*sign(i1);


slave=AddSlave('slave',i1(:));

i0=find(elt(cEGI,3)>0);
if ~isempty(i0) % standard rigid links
  % tx_slave = tx_master + x(3) * ry_master - x(2)* rz_master
  % ty_slave = ty_master - x(3) * rx_master + x(1)* rz_master
  % tz_slave = tz_master + x(2) * rx_master - x(1)* ry_master

  i1=[elt(cEGI,1)];
  i1=[i1+.01 i1+.02 i1+.03 i1+.04 i1+.05 i1+.06 ... % diagonal terms
    i1+.06 i1+.04 i1+.05 ... % terms in +x
    i1+.05 i1+.06 i1+.04 ... % terms in -x
  ];
  i2=[elt(cEGI,2)];i2=[i2+.01 i2+.02 i2+.03 i2+.04 i2+.05 i2+.06 ... % diagonal
   i2+.02 i2+.03 i2+.01 ... % terms in +x
   i2+.03 i2+.01 i2+.02 ... % terms in -x
  ];
  x = node(NNode(elt(cEGI,2)),5:7)-node(NNode(elt(cEGI,1)),5:7);
  i2=feval(slave.getPosFcn,slave,i2(:));
  r1=[ones(length(cEGI)*6,1);x(:);-x(:)];
  i3=find(i2&r1);

  II=[II;i2(i3)];TT=[TT;r1(i3)];
  i4=feval(nd.getPosFcn,nd,i1(i3));
  JJ=[JJ;i4(:)];

end
i0=find(elt(cEGI,3)<0);
if ~isempty(i0) % DOF equality
  error('DOF equality not reimplemented'); 
%  i1=[elt(cEGI,1)];
%  i1=[i1+.01 i1+.02 i1+.03 i1+.04 i1+.05 i1+.06 ... % diagonal terms
%    i1+.06 i1+.04 i1+.05 ... % terms in +x
%    i1+.05 i1+.06 i1+.04 ... % terms in -x
%  ];
%  i2=[elt(cEGI,2)];i2=[i2+.01 i2+.02 i2+.03 i2+.04 i2+.05 i2+.06 ... % diagonal
%   i2+.02 i2+.03 i2+.01 ... % terms in +x
%   i2+.03 i2+.01 i2+.02 ... % terms in -x
%  ];
%  x = node(NNode(elt(cEGI,2)),5:7)-node(NNode(elt(cEGI,1)),5:7);
%  i2=full(slave(round(i2(:)*100)-100));i3=find(i2);
%  r1=[ones(length(cEGI)*6,1);x(:);-x(:)];

%  II=[II;i2(i3)];JJ=[JJ;full(nd(round(i1(i3)*100)-100))];TT=[TT;r1(i3)];
end

case 'mpc';    
%% #MPC Multiple point constraints - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     r1=Case.Stack{j0,3}; 
     if ~isfield(r1,'c')&& isfield(r1,'DOF')
       error('Not a proper MPC format');
     elseif (~isfield(r1,'c')||isempty(r1.c))&&isfield(r1,'type')
       eval('[model,r1]=fe_caseg(''connection'',model,Case.Stack{j0,3},r1);');
     end
     % implementation of NASTRAN rbe3
     %elseif ~isempty(findstr(lower(Case.Stack{j0,2}),'rbe3'))
     %  warning('RBE3 currently ignored');
     % other case of MPC
       i1=[]; 
       if ~isfield(r1,'slave')||isempty(r1.slave);
          [r1.c,r1.slave]=feutil('fixMpcMaster',r1.c);
       else;[r1.c,r1.slave]=feutil('fixMpcMaster',r1.c,r1.slave);
       end
       try; % eliminate fixed dof
        dofs=[]; %mech_dof(slave.getListedDofPos(slave));
        if ~isempty(RunOpt.fixdof); dofs=[dofs;RunOpt.fixdof]; end
        if ~isempty(dofs);    i2=1; ic=0; icm=size(r1,2);
         while ~isempty(i2)&&ic<=icm % fix mpc may not find a free master with its technique
          i2=fe_c(r1.DOF(r1.slave),dofs,'ind');  ic=ic+1;
          if 1==2%~isempty(i2)
           % Eliminate fixed rows
           i4=fe_c(r1.DOF,RunOpt.fixdof,'ind',2);% non fixed cols
           i2=~any(r1.c(i2,i4),2);
           if nnz(i2); r1.c(i2,:)=[];r1.slave(i2)=[];end
           i2=fe_c(r1.DOF(r1.slave),RunOpt.fixdof,'ind');
          end
          if ~isempty(i2)
           i3=r1.slave(i2); r2=fe_c(r1.DOF,r1.DOF(i3));
           r1.c=r1.c-r1.c(:,i3)*r2; % remove component on fixed slave
           r1.slave=[]; fprintf('Fixed mpc.slave of ''%s''',Case.Stack{j0,2});
           [r1.c,r1.slave]=feutil('fixMpcMaster',r1.c); break
          end
         end
        end
       end

       i1=1:size(r1.c,1);i2=r1.slave;
       if isempty(i2); error('No slave in %s',Case.Stack{j0,2})
       elseif norm(r1.c(:,i2)-speye(size(r1.c,1)),'inf')>1e-10
          error('mpc,%s Slave DOFs are not independent',Case.Stack{j0,2});
       end
       i3=[];i3(flipud(i1(:)))=flipud(i2(:));
       r2=sparse(i3,1:length(i3),1,size(r1.c,2),size(r1.c,1))-r1.c';

       slave=AddSlave('slave',r1.DOF(i3));
       [i2,i1,r2]=find(r2);
       i4=i2; i2=feval(nd.getPosFcn,nd,r1.DOF(i2)); % all dofs
       i5=feval(nd.getPosFcn,nd,r1.DOF(i3)); % slave
       if any(i5==0);
         error('Problem with MPC DOFs not in model.DOF %s\n entry %s -> ignored\n', ...
             comstr(unique(fe_c(r1.DOF(i3(i5==0))),'stable'), ...
             -30,struct('NoClip',2)),Case.Stack{j0,2})
       end
       i1=i5(i1); 
       if nnz(i2)~=length(i2); 
         warning('Problem with MPC DOFs not in model.DOF %s\n entry %s -> ignored', ...
             comstr(unique(fe_c(r1.DOF(i4(i2==0))),'stable'), ...
             -30,struct('NoClip',2)),Case.Stack{j0,2})
       else; II=[II;i1];JJ=[JJ;i2];TT=[TT;r2]; %#ok<AGROW>
       end
case 'rbe3' % #RBE3 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

     r1=fe_mpc('fixrbe3',Case.Stack{j0,3},model);
     if isfield(r1,'data');r1=r1.data;end
     for j1=1:size(r1,1)
try;
       % THIS NEEDS DOCUMENTATION/FORMALIZATION XXX
       % i2(:,1) slave [Id;NodeId;Dof], i2(:,2:end) =[weight;DOF;Node] 
       % Constraint on loads: trans: Fs=1/w Sum(w Fm), Moments Ms=1/(su wj rj^2) Sum(wj rj Mm)
       % Fslave=Sc' Fmaster, [Sc' -I] {Fs;Fm} = 0  ~ T'*F = 0
       %  Find kinematic basis for qs for load to verify constraint
       % as reduction T'T=[Sc*Sc' -Sc; -Sc' I];
       % guyan like kinematic basis [ -(Sc*Sc')^-1 Sc ; I]
       % beware of unconstrained master rotations
       try;
        i2=r1(j1,:); i2(~i2)=[];        
        i2=reshape(i2,3,length(i2)/3);
       catch
        error('RBE3 format is unknown');
       end
       if size(i2,2)==1; 
        sdtw('_nb','%i %i %i : single node RBE3',i2);continue
       end
       i1=NNode(i2(3,2:end)); i3=NNode(i2(2));
       if any(i1==0)||i3==0; 
        if i3==0; i3=i2(2); else; i3=[]; end;
        i3=[i3 i2(3,find(i1==0)+1)]; %#ok<AGROW>
        fprintf('Nodes %s not defined\n',sprintf(' %i',i3));
        sdtw('_err','Error in mpc %i',j1);
       end
       Lx=node(i1,5:7)-node(i3*ones(size(i2,2)-1,1),5:7); % ri
       Wj=i2(1,2:end); % 05/07/2017 now account for weights (gv)
       
       % Lc2=mean(sqrt(sum(Lx.^2,2)))^2; W=[1 1 1 Lc2 Lc2 Lc2];
       Kci=zeros(6,6*size(Lx,1)); Kcc=zeros(6);ind=zeros(1,size(Kci,2));
       CDOF=ones(6,1)*[i2(2,1) i2(3,2:end)]+ ...
         [.01;.02;.03;.04;.05;.06]*ones(1,size(i2,2));
       CDOF=CDOF(:); Sci=[eye(6) zeros(6,size(Kci,2))];
       % k=k=[Sc' -I];k=k'*k;k-[Sc*Sc' -Sc;-Sc' eye(6)]
       T=zeros(size(Lx,1)*6,6); Tw=T;
       for j2=1:size(Lx,1)
        i3=(j2-1)*6; 
        Sc=eye(6); Sc([17 6 10])=-Lx(j2,:); Sc([12 16 5])=Lx(j2,:); % 1,rj
        Wi=Wj(j2)*eye(6); %Wi=diag(W*i2(1,j2+1)); 
        Kcc=Kcc+Sc*Wi*Sc';
        Kci(:,i3+[1:6])=-Wi*Sc; %Sci(:,i3+[7:12])=Sc;
        i4=abs(sprintf('%i',i2(2,j2+1)))-48;ind(i3+i4)=1;
        i4=setdiff(1:6,i4);T(i3+i4,i4)=eye(length(i4));
        Tw(i3+i4,i4)=sqrt(Wj(j2))*eye(length(i4));
       end
       % condense on master DOFs, in1 contains slave DOFs (center + rot)
       % 
       if 1==1
        K1=Kci*T;K1=Kcc+K1+K1'+Tw'*Tw;ind=logical(ind);
        r2=null(full(K1),'r'); in1=1:6;
        if ~isempty(r2); in1=setdiff(in1,find(any(abs(r2(4:6,:)),2))+3);end
        Q=-pinv(full(K1(in1,in1)))*Kci(in1,ind);
        [in1,in3]=intersect(abs(sprintf('%i',i2(3,1)))'-48,in1);
        slave=AddSlave('slave',CDOF(in1));
        [i1,i2,r2]=find(Q(in1,:));ind=find(ind)+6; %zzzgvdr
        i3=feval(nd.getPosFcn,nd,CDOF(in1));i1=i3(i1);    % slave
        i3=feval(nd.getPosFcn,nd,CDOF(ind));i2=i3(i2);    % master

       else % Old /!\ slower and very intense on null and pinv if many master DOF
        Kci=sparse(Kci);Kcc=sparse(Kcc); K=[Kcc Kci;Kci' speye(length(Kci))];
        in1=[1:6 find(~ind)+6]; K1=K(in1,in1);
        r2=null(full(K1),'r');
        if ~isempty(r2); in1=setdiff(in1,find(any(abs(r2(4:6)),2))+3);end
        in2=1:size(K,1); in2(in1)=0; in2=find(in2);
        Q=-pinv(full(K(in1,in1)))*K(in1,in2);
        [in1,in3]=intersect(abs(sprintf('%i',i2(3,1)))'-48,in1);
        slave=AddSlave('slave',CDOF(in1));
        [i1,i2,r2]=find(Q(1:length(in1),:));
        i3=feval(nd.getPosFcn,nd,CDOF(in1));i1=i3(i1);    % slave
        i3=feval(nd.getPosFcn,nd,CDOF(in2));i2=i3(i2);    % master
       end

       if isempty(II); II=i1;JJ=i2;TT=r2;
       else; II=[II;i1];JJ=[JJ;i2];TT=[TT;r2];
       end

catch;
  if sp_util('diag')>10;rethrow(lasterror);end
  fprintf('RBE3 (%s,%i) computation failed\n',Case.Stack{j0,2},r1(j1,1));
end 
     end
case {'keepdof','par','dofload','fsurf','fvol','sensdof','dofset', ...
      'cyclic','fixdof','nastran','sensstrain','info','pcond','pred'}

otherwise
 sdtw('_nb','''%s'' not yet supported by fe_mpc',Case.Stack{j0,1});
end



end % loop on stack

% #add_masters_to_the_II,JJ,TT_list - - - - - - - - - - - - - - - - - - - - - 
% one adds masters used to define constraints and masters in mdof

i0=(1:N)'; 
i0(fe_c(mech_dof,mdof,'indskipcheck',2))=0; % eliminate masters not in mdof
i1=find(sparse(JJ,1,1));i0(i1)=i1; % reintroduce masters used for constraints
% [i1,i2,i3]=find(slave);i0(i3)=0; % get rid of slaves
if ~isempty(slave); i0(feval(slave.getListedDofPos,slave))=0; end
i1=find(i0);II=[II;i1]; JJ=[JJ;i1];% keep the remaining DOFs
TT=[TT;ones(size(i1))];

% close things up ------------------------------------------------------------

T=sparse(II,JJ,TT,N,N);

% recurse while some masters are slaves
i1=mech_dof(any(T)); 
if isempty(slave); i2=[]; i3=0;
else; i2=feval(slave.getPosFcn,slave,i1);i2=i2(i2~=0); i3=length(i2);
end

while ~isempty(i2)
 % slave masters
 T1t=T(i2,:);
 T=T+T(:,i2)*(T1t-sparse(1:length(i2),i2,1,length(i2),N));

 i1=mech_dof(any(T)); 
 i2=feval(slave.getPosFcn,slave,i1);i2=i2(i2~=0); 
 if length(i2)==i3
  T1=T(:,i2);
  i2(max(max(T1,[],1),0)-min(min(T1,[],1),0)<eps)=[];clear T1
  if length(i2)==i3;
      disp(fe_c(mech_dof(i2)))
      error('The above slave masters could not be eliminated');
  else; i3=length(i2);
  end
 else;i3=length(i2);
 end
end

if ~isempty(ikeep) % some keepdof
  if ~isempty(cGL)  % allow for non mechanical DOFs in T
   if size(cGL,1)==size(T,2);T=cGL'*T*cGL; 
    else;error('This is not yet supported');
   end
  end
  i1=sum(fe_c(mech_dof,mdof(ikeep)));
  i1=find(any(T).*i1);T=T(:,i1);
else
  if size(cGL,1)>1  % allow for non mechanical DOFs in T
   i1=find(any(T)); 
   % Built q_local= T q_{master,local}, desired q_global= T q_{master,local}
   if size(cGL,1)==size(T,2);T=cGL*T(:,i1);
   else;error('This is not yet supported');
   end
  else
   i1=find(any(T));T=T(:,i1);
  end

end

i2=fe_c(mech_dof,mdof,'indskipcheck');
i3=fe_c(mdof,mech_dof,'indskipcheck',2);

if isempty(i3) % no mdof is not in mech_dof
 T=T';T=T(:,i2);Case.T=T';
 %Case.T=T(i2,i1); 
 Case.DOF=mech_dof(i1);
else
 Case.DOF=[mech_dof(i1);mdof(i3)];
 T=[T(i2,:) spalloc(length(i2),length(i3),0);
         spalloc(length(i3),length(i1),0) speye(length(i3))];
 i4=size(T);
 [II,JJ,TT]=find(T);i1=[fe_c(mdof,mech_dof(i2),'ind');i3];
 II=i1(II); Case.T=sparse(II,JJ,TT,i4(1),i4(2));
 %Case.T([fe_c(mdof,mech_dof(i2),'ind');i3],:)=T;
end
out=Case; out1=mdof;

%% #AddSlave  ----------------------------------------------------------------
function slave=AddSlave(slave,i1);
if ischar(slave); eval(iigui(slave,'MoveFromCaller')); end
i1=i1(i1>0);
if length(unique(round(i1*100)))~=length(i1) 
  disp(sdtw('_clip 60 1','%s ',fe_c(find(sparse(round(i1*100),1,1)>1)/100)));
  st1=evalin('caller','Case.Stack{j0,2}');
  error('Same DOF is slave more than once, ''%s''',st1);
end

if any(feval(slave.getPosFcn,slave,i1)) %any(slave(round(i1*100)-100))
 r1=i1(feval(slave.getPosFcn,slave,i1)~=0);
  fprintf('\n%s\n\n', ...
   sdtw('_clip 80 1','%s ',fe_c(fix(r1)+round(rem(r1,1)*100)/100)));
 try;st1=evalin('caller','Case.Stack{j0,2}');catch;st1='';end
  error('''%s'' some slave DOFs are already slaves ',st1); 
else
 slave=feval(slave.appendDof,slave,i1);
end

