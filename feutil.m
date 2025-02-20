function [out,out1,out2,out3,out4]=feutil(varargin);

%FEUTIL Finite element model handling utilities
%
% <a href="matlab: sdtweb feutil">Documentation for feutil</a>
%
%	This function groups a number of mesh handling utilities with all
%	arguments given (unlike feplot, femesh, fecom which use context 
%	information).
%	sdtdef('epsl') gives the mesh precision (smallest distance considered)
%	whose default is 1e-6.
%
%	Accepted calls are (see the manual for details)
%
%       [AllNode,ind]=femesh('AddNode',OldNode,NewNode);
%       Pairs = feutil('closestnode',node1,node2) -> closest node matching
%       Mode      = feutil('dof2mode',def,dof) -> dof to [node disp] form
%       [def,dof] = feutil('mode2dof',mode)    -> [node disp] to dof form
%       [model,def]=feutil('DefLocalToGlocal',model,def) deformation in global
%       elt=feutil('divide group 1 withnode{x>10}',model) ->
%                                                    divide group 1 in 2 groups
%       [node,elt]=feutil('divide HexaTransition',node,elt,nodeId);
%       [EGID,eltnames]=feutil('egid',elt) or 
%            elt=feutil('egid',elt,EGID,eltnames)
%       [EltId]  = feutil('eltid [,fix]',elt)
%       [EltInd] = feutil('eltind',elt)
%       [ind,snode] = feutil('findnode (selection)',node,elt,elt0, ...)
%	  See the femesh findnode command for more details
%       [ind,els,indwithead] = feutil('findelt (selection)',node,elt,el0, ...)
%	  See the femesh findelt command for more details
%
%	[Line,nInf,node_ind]  = feutil('get Line',node,elt,nInf)
%	  node empty => use indices
%	Patch = feutil('get Patch',node,elt) : node empty => use node numbers
%	[ElemF,opt,ElemP]= feutil('getelemf',elt(EGroup(jGroup),:),[jGroup])
%	[mdof,tr,nword]  = feutil('get DOF',elt,adof,tr)
%	[nor,cg]  = feutil('get normal[elt,node]',node,elt)
%	Map  = feutil('get normalMap[elt,node]',node,elt)
%
%	feutil('info elt Name',elt) information on model description matrix
%       elt=feutil('joingroup1:2',model.Elt)
%	[el_match,point,el0_normal]=feutil('match',node,elt,el0) find faces 
%         in elt matching faces in el0.
%	[out]=feutil('[MatId ProId MPId]',elt) mat/pro ID for each element
%	elt=feutil('mpid',elt,mpid) set mat/pro ID for each element
%       Object : object constructors ( see more details in sdtweb('femesh') )
%          HoleInPlate, BeamLine, Mass, Beam, Quad, Hexa, Circle, Cylinder
%       [node,elt]=feutil('optimmodel',model) removes unused nodes from FEnode
%       model=feutil('optimeltcheck',model) : ill conditionning fixes
%       [node,elt]=feutil('optimnodenum',model) optimizes numbering
%       elt=feutil('orient',model) reorient elements if necessary
%       elt=feutil('removeelt (selection)',model) 
%	st=feutil('stringdof',mdof) string labels for DOFs
%       elt=feutil('trace2elt',ldraw) trace line for elt format
% 
%   See sdtweb     femesh
%	See also help  feplot
% 
% <a href="matlab: sdtweb feutil">Documentation for feutil</a>


%       Etienne Balmes, Guillaume Vermot des Roches, Jean-Philippe Bianchi
%       Copyright (c) 2001-2024 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       For revision information use feutil('cvs')

persistent OldRefQuad InfoMode Silent 
if isempty(OldRefQuad)
   OldRefQuad=sdtdef('OpenFEM.OldRefQuad-safe',0);
   Silent=0; %Silent=pmat(zeros(1));
end
%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*NOSEM,*NBRAK,*STREMP>
try;epsl=sp_util('epsl');catch;epsl=sdtdef('epsl'); end
ModelStack={};
[CAM,Cam]=comstr(varargin{1},1);carg=2;
if Silent==0; RSil=1; else; RSil=0; end
if ~sp_util('diag')&&Cam(end)==';'; Silent=1; end %of_time(-1,Silent,1);

if comstr(Cam,'silent'); out=Silent;
%% #Add ----------------------------------------------------------------------
elseif comstr(Cam,'add'); [CAM,Cam] = comstr(CAM,4);

%% #AddNode - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% ('addnode',OldNode,NewNode) - - - - - - - - - - - - - - - - - - - - - - - -
if comstr(Cam,'node'); [CAM,Cam] = comstr(CAM,5);

 r1=varargin{carg};carg=carg+1;
 r2=varargin{carg};carg=carg+1; 

 if ~isempty(r1) && any(~isfinite(r1(end,:))); r1 = r1(1:end-1,:); end
 if ~isempty(r2) && any(~isfinite(r2(end,:))); r2 = r2(1:end-1,:); end
 if size(r2,2)==3; r2=[[1:size(r2,1)]' zeros(size(r2,1),3) r2];end

 out=r1;  out1 = zeros(size(r2,1),1);

 [epsl,CAM,Cam]=test_epsl(epsl,CAM,Cam);
 if comstr(Cam,'from');   [opt,CAM,Cam]=comstr(Cam,'from','%i'); opt=opt-1;
 elseif ~isempty(out); opt=max(out(:,1)); else;opt=0; end

 if ~isempty(strfind(Cam,'new'))
  opt=comstr(CAM(strfind(Cam,'new')+4:end),[-1 opt]);
  out1=size(out,1)+[1:size(r2,1)]';
  out(end+[1:size(r2,1)],1:7)=[opt+[1:size(r2,1)]' r2(:,2:7)];

 elseif size(r2,1)>10||sp_util>5.0007 % try a new strategy based on sortrows

  % catch for SDT strategy to do the same thing
  if sp_util>5.0007;
   try;eval('[out,out1]=feutilb([''addnode_'' CAM],r1,r2,epsl,opt);');return;end
  end
  % i3 contains a list corresponding to sort of [r1;r2]
  if ~isempty(r1); r3=[r1(:,5:7);r2(:,5:7)];
  else; r3=[r2(:,5:7)];
  end
  [u,s,v]=svd(r3,0);r3=r3*v; % orient along main mesh directions
  [r4,i3]=sort(r3(:,1));r3=r3(i3,:);ind=i3; 

  inx=find(abs(diff(r3(:,1)))<epsl); % x is matching
  while ~isempty(inx)
   % one segment of matching x
   if isscalar(inx); i2=1; else;i2=min(find(diff(inx)>1));end %#ok<MXFND>
   if isempty(i2); i2=length(inx);end
   iny=inx(1):inx(i2)+1;
   [r4,i4]=sort(r3(iny,2));r3(iny,:)=r3(iny(i4),:); 
   i3(iny)=i3(iny(i4)); ind(iny)=ind(iny(i4));
   iny=iny(abs(diff(r4))<epsl);

   while ~isempty(iny) % y is matching
    if isscalar(iny); i5=1; else;i5=min(find(diff(iny)>1));end %#ok<MXFND>
    if isempty(i5); i5=length(iny);end
    inz=iny(1):iny(i5)+1;
    [r4,i4]=sort(r3(inz,3)); r3(inz,:)=r3(inz(i4),:); 
    i3(inz)=i3(inz(i4));ind(inz)=ind(inz(i4));
    inz=inz(abs(diff(r4))<epsl);

    while ~isempty(inz) % z is matching
     if length(inz)>1; i6=min(find(diff(inz)>1));else;i6=1;end %#ok<MXFND>
     if isempty(i6); i6=length(inz);end
     i7=inz(1):inz(i6)+1;
     [r4,i4]=min(i3(i7)); % verify coincidence
     r4=sqrt(sum((r3(i7,:)-ones(size(i7(:)))*r3(i7(i4),:)).^2,2))<epsl;
     for j1=find(r4(:)'==0&i3(i7)'>size(r1,1)) % when failed do a direct search
       r5=sqrt(sum((r3(i7,:)-ones(size(i7(:)))*r3(i7(j1),:)).^2,2));
       i5=find(r5<epsl);i3(i7(j1))=min(i3(i7(i5)));
     end % when failed do a direct search
     i7=i7(r4~=0);  i3(i7)=min(i3(i7));
     inz=inz(i6+1:end);
    end

    iny=iny(i5+1:end);
   end % y is matching
   inx=inx(i2+1:end);
  end % x is matching

  % check coincident nodes
  i4=find(i3<size(r1,1));  

  i4=find(i3>size(r1,1));
  if ~isempty(i4) % 
    i4=unique(i3(i4));
    out=[r1;opt+[1:length(i4)]' r2(i4-size(r1,1),2:end)];
    % new indices of r2 in out
    nind=[[1:size(r1,1)]';i4];nind(nind)=1:length(nind);
    out1(ind)=nind(i3); out1=out1(size(r1,1)+1:end);
  else
    out=r1;
    nind=[1:size(r1,1)];nind(nind)=1:length(nind);
    out1(ind)=nind(i3);
    out1=out1(size(r1,1)+1:end);
  end
  
 else % old strategy that scans through each node
  opt=opt(1)+1;
  for j1=1:size(r2,1)
    if isempty(out); out1(j1)=1;out=[opt r2(j1,2:7)];opt=opt+1;
    else
     if sdtdef('verm')>=704
      [r3,i1]=min(sqrt(sum(bsxfun(@minus,out(:,5:7),r2(j1,5:7)).^2,2)));
      if r3<epsl;out1(j1) = i1(1);
      else;out = [out;opt(1) r2(j1,2:7)];opt=opt+1; out1(j1,1)=size(out,1);
      end
     else % this old old code does not provide the closest node in general !!
      i1 = find(abs(out(:,5)-r2(j1,5))<epsl);
      if ~isempty(i1);   i1 = i1(abs(out(i1,6)-r2(j1,6))<epsl); end
      if ~isempty(i1);   i1 = i1(abs(out(i1,7)-r2(j1,7))<epsl); end
      if ~isempty(i1);   out1(j1) = i1(1);
      else
        out = [out;opt(1) r2(j1,2:7)];opt=opt+1; out1(j1,1)=size(out,1);
      end
     end % old
    end % of ~isempty(out)
  end
 end % 'new or not'

%% #AddTest - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% ('addtest',model,model2) - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'test'); 
    
RunOpt.Silent=strcmp(Cam(end),';');
if RunOpt.Silent; CAM(end)='';Cam(end)=''; RunOpt.Silent=RunOpt.Silent&&~sp_util('diag'); end
[CAM,Cam] = comstr(CAM,5);
[CAM,Cam,RunOpt.noorig]=comstr('noori',[-25 3],CAM,Cam); 

model=varargin{carg};carg=carg+1;
el0=varargin{carg};carg=carg+1; 
i1=strfind(Cam,'keeptest'); 
if ~isempty(i1); RunOpt.KeepTest=1;else;RunOpt.KeepTest=0;end
RunOpt.EGID=[];i1=strfind(Cam,'-egid');
if ~isempty(i1)
 [RunOpt.EGID,i2,i3,i4]=sscanf(CAM(i1+5:end),'%i',1);
 CAM(i1:i1+3+i4)='';[CAM,Cam]=comstr(CAM,1);
end
NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
el0.OrigNodeId=el0.Node(:,1);
if isfield(el0,'bas'); [el0.Node,el0.bas]=basis('nodebas',el0);end
FEnode=model.Node; FEelt=model.Elt;
[CAM,Cam,i1]=comstr('epsl',[-25 2],CAM,Cam); if ~isempty(i1);epsl=i1;end

if isempty(el0); opt(2)=Inf;  st1='same';% - - - - - - - - - - - -  same nodes 
elseif comstr(Cam,'combine') % combine node sets. Nodes with same number must
% be at same location or error

i2=find(~ismember(el0.Node,model.Node,'rows'));
[r3,i4,i5]=intersect(el0.Node(i2,1),model.Node(:,1));
if ~isempty(r3); 
 for j1=1:length(i4)
  r3(j1)=norm(el0.Node(i2(i4(j1)),:)-model.Node(i5(j1),:));
  if r3(j1)<epsl; i4(j1)=0;i5(j1)=0;r3(j1)=0; end
 end
 i4=i4(i4~=0);i5=i5(i5~=0);
 if ~isempty(i4)
  sdtw('_clip 60 1','%i ',[model.Node(i5,1)]')
  error('Some nodes with identical numbers are different');
 end
 i2=find(~ismember(el0.Node(:,1),model.Node(:,1),'rows'));
 FEnode(end+[1:length(i2)],:)=el0.Node(i2,:);
else % append new nodes
 FEnode(end+[1:length(i2)],:)=el0.Node(i2,:);
end
nind=sparse(FEnode(:,1),1,1:size(FEnode,1)); i1=nind(el0.Node(:,1));
opt=[-1 0];st1='combine';

elseif ~isempty(strfind(Cam,'merge')) %  #AddTestMerge : merge nodes - - - 
 
 [FEnode,i1]=feutil(sprintf('AddNode epsl %.15g',epsl),model.Node,el0.Node);
 [CAM,Cam,i2]=comstr('merge',[-25 3],CAM,Cam);
 [CAM,Cam,i2]=comstr('-edge',[-25 3],CAM,Cam);
 if i2 % mid edge nodes are forced to be common
  i2=feutil('selelt seledge all & eltname beam3',model);i2(~isfinite(i2),:)=[];
  i3=feutil('selelt seledge all & eltname beam3',el0.Node,el0.Elt);i3(~isfinite(i3),:)=[];
  if ~isempty(i2)&&~isempty(i3)
   nind=sparse(el0.Node(:,1),1,1:size(el0.Node,1));
   [i4,i5]=ismember(sort(FEnode(i1(nind(i3(:,1:2)))),2),sort(i2(:,1:2),2),'rows'); % used edges
   el0.Node(nind(i3(i4,3)),2:7)=model.Node(NNode(i2(i5(i4),3)),2:7);
   [FEnode,i1]=feutil(sprintf('AddNode epsl %.15g',epsl),model.Node,el0.Node);
  end
 end

 r2=FEnode(i1,5:7)-el0.Node(:,5:7);
 if max(sqrt(sum(r2.^2,2)))>epsl*sqrt(3)
  FEnode=model.Node;
  [FEnode,i1]=feutil(sprintf('AddNode epsl %.15g',epsl),FEnode,FEnode);
  if length(unique(i1))<size(FEnode,1);
   sdtw('_nb','Renumbering non-merged nodes in initial model');
   model=feutil('renumber',model,FEnode(i1,1)); FEnode=model.Node;
   [FEnode,i1]=feutil(sprintf('AddNode epsl %.15g',epsl),FEnode,el0.Node);
  else;
   error('This case is seen as a bug');
  end
 end
 % attempt to use el0 node numbers if NodeId is free
 [i2,i3]=setdiff(FEnode(i1,1),model.Node(:,1)); % new node numbers
 %[i2,i4]=setdiff(el0.Node(i3,1),model.Node(:,1)); % was
 [i2,i4]=setdiff(el0.Node(i3,1),FEnode(:,1)); % NodeId unused in FEnode
 FEnode(i1(i3(i4)),1)=el0.Node(i3(i4),1);
 opt=[-1 0]; st1='merge';

else                        %  - - - - - - - - - - - - - - - -  node shift 
 r2 = max(el0.Node(:,1)); if length(NNode)<r2; NNode(r2)=0;end
 opt = comstr(CAM,[-1 0]); opt(2)=0; i1=NNode(el0.Node(:,1)); 
 if opt(1)==0
   if all(i1) && norm(FEnode(i1,2:7)-el0.Node(:,2:7))<epsl 
     st1='sets identical -> unchanged';opt(2)=Inf;
   elseif any(i1); opt(1)= max(FEnode(:,1));
       st1=sprintf('shift by %i',opt(1));
   else;st1='sets disjoint and combined'; 
   end
 else;st1=sprintf('shift by %i',opt(1));
 end
end
if ~RunOpt.Silent&&~Silent;
    fprintf('Addtest node %s\n',st1);
end

if isinf(opt(2))  % same nodes
  i1=NNode(el0.Node(:,1));
elseif opt(1)==-1
elseif opt(1)~=0  % node shift
  el0.Node(:,1)=el0.Node(:,1)+opt(1);
  if length(NNode)<max(el0.Node(:,1)); NNode(max(el0.Node(:,1)),1)=0;end
  i1=NNode(el0.Node(:,1));
  if all(i1) && norm(FEnode(i1,5:7)-el0.Node(:,5:7))<epsl 
    warning('test nodes kept unchanged');opt(2)=Inf;
  end
end


if ~isempty(FEnode)&&~isfinite(FEnode(end,5)); FEnode=FEnode(1:end-1,:);end
% NNode giving old node positions in i1 (new node numbers)
% i1 gives the associated node numbers
if opt(1)==-1 % merged nodes
  NNode=[]; NNode=sparse(el0.Node(:,1)+1,1,1:size(i1,1));
elseif isfinite(opt(2))&&isfinite(opt(1))  % new test nodes added
  i1 = size(FEnode,1)+[1:size(el0.Node,1)]';FEnode=[FEnode;el0.Node];
  NNode=[];NNode(el0.Node(:,1)-opt(1)+1)=1:size(el0.Node,1);
elseif ~isempty(el0.Node)&&~isfinite(opt(2))
  NNode=[];NNode(el0.Node(:,1)-opt(1)+1)=1:size(el0.Node,1);
else;error('Bug NNode should be redefined');
end

elt=el0.Elt;
i4=size(FEelt,1);
if i4&&isempty(strfind(Cam,'keeptest'));
 FEelt=feutil('removeelt EGID -1',FEnode,FEelt);
 if size(FEelt,1)~=i4
    sdtw('_nb','Removed test wire-frame element groups');
 end
end


if isempty(elt)
 FEelt(end+1,1:8)=[Inf abs('mass2') 0 -1];
 FEelt(end+[1:size(el0.Node,1)],1)=el0.Node(:,1);
elseif isfinite(elt(1,1)) % this is a trace line
 i3=find(elt(:,83:end));
 i2=elt(:,83:end);i3=find(i2); i2(i3)=FEnode(i1(NNode(i2(i3)+1)),1);
 elt(:,83:end)=i2; 
 FEel0=feutil('trace2elt',elt); femesh('addsel',';');

else
    [EGroup,nGroup]=getegroup(elt);
    for jGroup=1:nGroup
     ElemF= feutil('getelemf',elt(EGroup(jGroup),:),jGroup);
     cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;
     if ~isempty(RunOpt.EGID); 
      elt(EGroup(jGroup),length(ElemF)+3)=RunOpt.EGID;
     end
     i2=fe_super('node',ElemF); % node indices
     i3=full(NNode(elt(cEGI,i2)+1)); i4=i3~=0; i3(i4)=FEnode(i1(i3(i4)),1);
     %i3(i3~=0)=FEnode(i1(i3(i3~=0)),1); % horrible sparse line indexing !
     elt(cEGI,i2)=reshape(i3,length(cEGI),length(i2));
    end
    FEelt(end+[1:size(elt,1)],1:size(elt,2))=elt;
end

out1=i1;  % Return new nodes
% Attempt to merge the model information

MergePlIl=fe_mat('@MergePlIl'); % add non intersecting and different pl/il entries
model=feval(MergePlIl,'pl',model,el0); % no renum. performed
model=feval(MergePlIl,'il',model,el0); % no renum. performed

% Refine the OrigNumbering Output
if ~isequal(out1,el0.OrigNodeId) 
 i1=stack_get(model,'info','OrigNumbering','getdata');
 i1=[i1;int32([el0.OrigNodeId FEnode(out1,1)])];
 model=stack_set(model,'info','OrigNumbering',i1);
end
% attempt to combine stack and nmap information
try; if RunOpt.Silent; sts=';'; else;sts=''; end
 model=feutilb(['AddTestStack' sts],model,el0);
catch;
 if isfield(el0,'Stack') && ~isempty(el0.Stack)
  if ~isfield(model,'Stack'); model.Stack=el0.Stack; 
  else; model=stack_set(model,el0.Stack); 
  end
 end
end

if RunOpt.noorig;
 model=stack_rm(model,'info','OrigNumbering');
end

out=model; out.Elt=FEelt; out.Node=FEnode;

%% #AddSet - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% ('AddSet',model,name,selection) - - - - - - - - - - - - - - - - - - - - - -
% ('AddSet-Append',model,name,{selection},{subnames}) - - - - - - - - - - - - 
elseif comstr(Cam,'set'); [CAM,Cam] = comstr(CAM,4);
 % Global strategy old/new version
 % Append triggers new set format
 % If set already exists with data field, convertion of exiting
 % No append: old behavior: old format and erase existing
 RO=struct('Append',0,'get',0,'ID',[],'typ',{{}},'Types',{{}},'cell',0,...
  'SColor',[],'cur',1,'New',0,'doNodes',1); % default
 [CAM,Cam,RO.FromInd]=comstr('fromind',[-25 3],CAM,Cam); % assume eltind provided in num eltsel
 [CAM,Cam,RO.New]=comstr('new',[-25 3],CAM,Cam);
 [CAM,Cam,i1]=comstr('-facelin',[-25 3],CAM,Cam);if i1;RO.FaceCmd='facelin';end
 model=varargin{carg};carg=carg+1;
 RO.name=varargin{carg};carg=carg+1;
 if ~isempty(RO.name); r2=stack_get(model,'set',RO.name,'get'); else; r2=[]; end
 if carg<=nargin; sel=varargin{carg};carg=carg+1;else; sel=''; end
 % Find ID if given
 i1=strfind(Cam,'-id');
 if ~isempty(i1);
  [r1,i2,i3,i4]=sscanf(CAM(i1+3:end),'%i'); if ~isempty(r1); RO.ID=r1; end
  CAM(i1:i1+1+i4)=''; [CAM,Cam]=comstr(CAM,1);
 end
 i1=strfind(Cam,'-cell');
 if ~isempty(i1); RO.cell=1; CAM(i1:i1+3)=''; [CAM,Cam]=comstr(CAM,1); end
 i1=strfind(Cam,'-get');
 if ~isempty(i1); RO.get=1; CAM(i1:i1+3)=''; [CAM,Cam]=comstr(CAM,1); end
 %% AddSetInit
 i1=strfind(Cam,'-append');
 if ~isempty(i1);
  RO.Append=1; CAM(i1:i1+6)=''; [CAM,Cam]=comstr(CAM,1);
  [CAM,Cam,RO.NoNodes]=comstr('-nonodes',[-25 3],CAM,Cam); RO.doNodes=~RO.NoNodes;
  RO.GName=RO.name;
  [eltid,model.Elt]=feutil('eltidfix',model);
  EEId=sparse(eltid+1,1,1:length(eltid));
  model.Node=fesuper('FNode',model);
  NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
  if RO.doNodes
   [u1,NEConn]=feval(feutilb('@EltNodeCon'),model.Elt,NNode,struct('useAllNodes',1));
  end
  if isempty(sel);[sel,RO.name,RO.typ,ConvFcn]=allEltSel(model);
  elseif carg<=nargin; %RO.name=varargin{carg}; carg=carg+1;
   li={'RO.name','RO.typ','ConvFcn'};
   while ~isempty(li)&&carg<=nargin
    eval(sprintf('%s=varargin{carg}; carg=carg+1; li{1}=[];',li{1}));
   end
  else; RO.name={}; % xxx strategy to robust name incrementation
  end
  if ~isempty(r2)&&isfield(r2,'data'); % convert set to new format into r1
   RO.SetNames={r2.name}; RO.info={r2.type};
   if isfield(r2,'ConvFcn'); ConvFcn={r2.ConvFcn}; end
  elseif ~isempty(r2); % exploit exiting set
   [ie,je,ke]=find(r2.SConn);
   if ~isequal(eltid,r2.EltId) % recast EltId
    if max(r2.EltId)+1>length(EEId);EE=EEId;EE(max(r2.EltId)+1)=0; else;EE=EEId; end
    ie=full(EE(r2.EltId(ie)+1)); iie=ie~=0;
    ie=ie(iie); je=je(iie); ke=ke(iie); clear EE iie
   end
   if ~isfield(r2,'NConn'); RO.doNodes=0; end
   if RO.doNodes
    [in,jn]=find(r2.NConn);
    if ~isequal(model.Node(:,1),r2.NodeId) % recast NodeId
     if length(r2.NodeId)>length(NNode);NN=NNode;NN(length(r2.NodeId))=0;else;NN=NNode; end
     in=full(NN(r2.NodeId(in))); iin=in~=0;
     in=in(iin); jn=jn(iin); clear NN iin
    end
   else; in=zeros(3,1); jn=zeros(3,1);
   end
   if isscalar(ie) %nnz(r2.SConn)==1; 
    ie(end+1,1)=zeros(1); je(end+1,1)=zeros(1); ke(end+1,1)=zeros(1); 
   end
   if RO.doNodes&&isscalar(in); in(end+1,1)=zeros(1); jn(end+1,1)=zeros(1); end
   if isfield(r2,'SColor'); RO.SColor=r2.SColor; end
    RO.SetNames=r2.SetNames; RO.Types=r2.type;
  else; % initialize from scratch
   ie=zeros(3,1); je=zeros(3,1); ke=zeros(3,1); 
   in=zeros(3,1); jn=zeros(3,1);  RO.SetNames={};
  end
  if any(ie); i1=find(ie~=0,1,'last'); RO.cur=size(r2.SConn,2)+1;
  else; i1=zeros(1); 
   if isfield(r2,'SConn')&&~isempty(r2.SConn); RO.cur=size(r2.SConn,2)+1; end
  end
  ije=[-1 i1+0]; jje=[-1 i1+0]; kje=[-1 i1+0];
  if any(in); i1=find(in~=0,1,'last'); else; i1=zeros(1); end
  ijn=[-1 i1+0]; jjn=[-1 i1+0]; 
  %if i1; RO.cur=max(jn)+1; end
 else % old behavior
  r2=[]; eltid=[]; RO.Out={};
 end
 
 if ~iscell(sel); sel={sel}; end
 if ~iscell(RO.name); RO.name={RO.name}; end

 if isempty(RO.typ)
  r2={'NodeId','EltId','FaceId','EdgeId'}; if ~RO.doNodes; r2(1)=[]; end
  i1=ismember(lower(r2),lower(strtrim(Cam)));
  if any(i1); RO.typ=r2{i1}; else; error('Set type %s unknown.',CAM); end
 end
 if ~iscell(RO.typ); RO.typ=repmat({RO.typ},1,length(sel)); end

 %% AddSetLoop
 if ~exist('ConvFcn','var'); ConvFcn=cell(length(sel),1); end
 for j1=1:length(sel)
  r3=sel{j1}; r4=[]; id=1; RO.Types{end+1}=RO.typ{j1}; stCv='';
  switch RO.typ{j1}
   case {'EltId','EltSet'}
    if isempty(eltid); [eltid,model.Elt]=feutil('eltidfix;',model); end
    if ischar(r3);r3=eltid(feutil(sprintf('FindElt %s',r3),model)); 
    elseif RO.FromInd; r3=eltid(r3); % eltind provided
    end
    
   case {'NodeId','NodeSet'}
    if ischar(r3);r3=feutil(sprintf('FindNode %s',r3),model); end
    if RO.Append; r4=r3; r3=[]; end % revert node/elt sel    
    
   case {'FaceId','FaceSet'} % do not do ConvFcn for Append
    if isempty(r3); r3='selface';end;  
    if ischar(r3)||(isnumeric(r3)&&any(~isfinite(r3(:,1))))
     if ischar(r3)
      if any(r3=='@'); stCv=sscanf(r3(find(r3=='@')+1:end),'%s',1); end
      if RO.Append % do not apply ConvFcn
       i1=find(r3=='@'); r3(i1+[0:length(stCv)])=[];
       sdtw('_nb','ConvFcn (set:%s) cannot be used for Append',RO.name{j1});
      end
      elt=feutil(sprintf('SelElt%s',r3),model);
      if ~isempty(strfind(lower(r3),'facelin'));RO.FaceCmd='facelin';end
     else; elt=r3;
     end
     r3=FaceIdFromSelface(elt); r3(r3(:,1)==0,:)=[];
    elseif RO.Append&&~isempty(ConvFcn{j1}) % r3 was previously defined but with ConvFcn
     % old implement cvsnum<=1.611
     try; r3=FaceIdFromSelface(r3,'resolve',model,ConvFcn{j1}); catch; r3=[]; end
     if isempty(r3)||~any(r3(:,2));
      sdtw('_nb','failed handling of face set %s,%s',RO.typ{j1},RO.name{j1});
      continue
     end
    end
    
   case {'EdgeId','EdgeSet'}
    id=2; % binary identifier
    if isempty(r3); r3='seledge'; end
    if ischar(r3)||(isnumeric(r3)&&any(~isfinite(r3(:,1))))
     if ischar(r3)
      if any(r3=='@'); stCv=sscanf(r3(find(r3=='@')+1:end),'%s',1); end
      if RO.Append % do not apply ConvFcn
       i1=find(r3=='@'); r3(i1+[0:length(stCv)])=[];
       sdtw('_nb','ConvFcn (set:%s) cannot be used for Append',RO.name{j1});
      end
      elt=feutil(sprintf('SelElt%s',r3),model);
     else; elt=r3;
     end
     [EGroup,nGroup]=getegroup(elt);r3=cell(nGroup,1);
     for jGroup=1:nGroup
      [ElemF,i1,ElemP]= getegroup(elt(EGroup(jGroup),:),jGroup);
      cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
      i1=feval(ElemP,'prop'); r4=elt(cEGI,:); % assumed beam1 or beam3
      r3{jGroup}=elt(cEGI,i1(2)+[1 2]); % output is not conform to beam elts
     end
     r3=vertcat(r3{:});
    end
  end % typ
  %% AddSetGen
  if RO.Append
   r1=RO.name(j1,:); if length(r1)<2; r1{1,2}=RO.typ{j1}; end
   try
   % Increase Global set, elements and nodes for all
   if ~isempty(r3) % increase SConn
    if (max(r3)+1)>length(EEId) % reset (if set gives rmed EltId)
     EEId=sparse(eltid+1,1,1:length(eltid),max(r3)+1,1);
    end
    if size(r3,2)>1 % Face/Edge identifiers
     r3=catSubEid(r3,id);
     kej=r3(:,2); r3=r3(:,1);
    else;kej=ones(length(r3),1);
    end
    sp_util('setinput',ie,full(EEId(r3(:,1)+1)),ije,'ie');
    sp_util('setinput',je,RO.cur+0*r3(:,1),jje,'je');
    sp_util('setinput',ke,kej,kje,'ke');
    if isempty(r4)&&RO.doNodes
     r5=EEId(r3(:,1)+1);
     if any(r5)
      [un1,un2,r5]=find(r5); % r5 is a vector, remove null entries
      [ii,jj,kk]=find(NEConn(:,r5)); r4=unique(ii); clear ii jj kk
     end
    end
     %r4=find(logical(sum(NEConn(:,EEId(r3(:,1)+1)),2))); end
   end
   if ~isempty(r4)&&RO.doNodes
    if max(r4)>length(NNode) % reset (if set gives rmed NodeId)
     NNode=sparse(model.Node(:,1),1,1:size(model.Node,1),max(r4),1);
    end
    sp_util('setinput',in,full(NNode(r4)),ijn,'in'); % increase NConn
    sp_util('setinput',jn,RO.cur+0*r4(:,1),jjn,'jn');
   end
   catch; sdtw('_nb','Failed to handle  set %s,%s',r1{:});
   end
   RO.SetNames=[RO.SetNames;r1];
   RO.cur=RO.cur+1;
  else % old behavior, generate data field
   r1=struct('type',RO.typ{j1},'name',RO.name{j1},'data',r3,'ID',RO.ID);
   if isfield(RO,'FaceCmd');r1.FaceCmd=RO.FaceCmd;end
   if ~isempty(stCv); r1.ConvFcn=stCv; end
   if ~isempty(r1.name); 
    if RO.New&&sp_util('issdt');[model,r1.name]=fe_def('stacknew',model,'set',r1.name,r1);
    else;model=stack_set(model,'set',r1.name,r1); 
    end
   elseif ~RO.get; sdtw('_nb','empty set name, not supported');
   end
   RO.Out{end+1}=r1;
  end
 end % loop
 %% AddSetOutput
 if RO.Append
  % Reorder to get contiguous groups (SetNames(:,2))
  [r3,i3]=unique(RO.SetNames(:,2)); i3=sort(i3); r3=RO.SetNames(i3,2);
  [i4,i5]=ismember(RO.SetNames(:,2),r3); i4=find(i4); i6=sort(i5,'ascend');
  st=cell(size(RO.SetNames,1),2);
  for j1=1:length(r3); st(i6==j1,:)=RO.SetNames(i5==j1,:); end
  [u1,i6]=sort(i5,'ascend'); i4=i4(i6);
  i4=sparse(i4,1,1:length(i4));
  % now build set 
  iie=ie~=0; 
  r1=struct('type',{reshape(RO.Types(i6),[],1)},'name',RO.GName,'SetNames',{st},...
   'EltId',eltid,'NodeId',model.Node(:,1),'SConn',...
   sparse(ie(iie),full(i4(je(iie))),ke(iie),length(eltid),size(st,1)));
  if RO.doNodes
   iin=in~=0; kn=ones(length(find(iin)),1);
   r1.NConn=logical(sparse(in(iin),full(i4(jn(iin))),kn,size(model.Node,1),size(st,1)));
  end
  if ~isempty(RO.SColor)
   RO.SColor=[RO.SColor;.8*ones(length(r1.type)-size(RO.SColor,1),3)];
   r1.SColor=RO.SColor(i6,:);
  end
  if RO.get; out=r1;
  else % retrhow model and possibly generated set
   out=stack_set(model,'set',RO.GName,r1); if nargout>1; out1=r1; end
  end
 else
  if ~RO.cell&&isscalar(RO.Out); RO.Out=RO.Out{1}; end
  if RO.get; out=RO.Out; if nargout>1&&exist('elt','var'); out1=elt; end
  elseif nargout==0; clear out;
  else; out=model; if nargout>1; out1=RO.Out; end
  end
 end
  
%% #AddElt - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
elseif comstr(Cam,'elt'); [CAM,Cam] = comstr(CAM,4);
 out=varargin{carg};carg=carg+1;
 if ischar(out); out=[]; carg=carg-1;end
 elt=varargin{carg};carg=carg+1;
 if ischar(elt); elt=[Inf abs(elt)];
  r2=varargin{carg};carg=carg+1;
  if ~isempty(r2)&&any(~isfinite(r2(:,1)));
    error(' ''name'',val cannot contain Inf in val');
  end
  elt(end+(1:size(r2,1)),1:size(r2,2))=r2;
 end
 
 % apply new eltids to added elements and output convertion tables
 if ~isempty(strfind(Cam,'-newid'))
  if isempty(elt); out1=[];
  else
   if isempty(out); eltid=0; 
   elseif isfield(out,'Elt')&&isempty(out.Elt); eltid=0; 
   else; eltid=feutil('eltidskipcheck;',out); 
   end
   elid1=feutil('eltidskipcheck;',elt); elid1(:,2)=0;
   in1=isfinite(elt(:,1)); elid1(in1,2)=max(eltid)+(1:sum(in1));
   elt=feutil('eltid-elt',elt,elid1(:,2));
   out1=elid1;
  end
 end
 
 if isfield(out,'Elt'); out.Elt(end+(1:size(elt,1)),1:size(elt,2))=elt;
 else;out(end+(1:size(elt,1)),1:size(elt,2))=elt;
 end
 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
else;error('Not a known add command');
end

%% #ClosestNode --------------------------------------------------------------
elseif comstr(Cam,'closestnode');  [CAM,Cam]=comstr(CAM,6);

r1=varargin{carg};carg=carg+1;
r2=varargin{carg};carg=carg+1;

if size(r1,2)==7; out1=r1(:,1); r1=r1(:,5:7); else; out1=[1:size(r1,1)]';end
if size(r2,2)==7; i2=r2(:,1); r2=r2(:,5:7); else; i2=[1:size(r2,1)]';end

i1=zeros(size(r1,1),1);
for j1=1:size(r1,1)
 r3=(r2-r1(ones(size(r2,1),1)*j1,:)).^2*[1;1;1];
 [r3,i1(j1)]=min(r3);
 out1(j1,2)=i2(i1(j1));
end
out=i1;
%% #Elt,Mat.ID related commands 
%% #EltInd ------------------------------------------------------------------2
elseif comstr(Cam,'eltind');  [CAM,Cam]=comstr(CAM,7);

 elt = varargin{carg};carg=carg+1;
 if isfield(elt,'Elt'); elt=elt.Elt; 
 elseif isfinite(elt(1))&&size(elt,2)==7; elt = varargin{carg};carg=carg+1;end
 out=[1:size(elt,1)]';

 [EGroup,nGroup]=getegroup(elt);
 for jGroup=1:nGroup
   [ElemF,opt]= feutil('getelemf',elt(EGroup(jGroup),:),jGroup);
   [i2,i3]=fe_super('prop',ElemF);
   if i3==2; out(EGroup(jGroup))=0; end
 end

%% #EltId -------------------------------------------------------------------2
elseif comstr(Cam,'eltid');  [CAM,Cam]=comstr(CAM,6);

model = varargin{carg};carg=carg+1; %if isa(model,'v_handle'); model=model.GetData; end
if isfield(model,'Elt'); elt=model.Elt; else; elt=model; end
if isfinite(elt(1))&&size(elt,2)==7; elt = varargin{carg};carg=carg+1;end

i1=strfind(Cam,'-elt'); i2=strfind(Cam,'-model');
if ~isempty(i1); 
 RunOpt.Back='elt';CAM(i1+[0:3])='';[CAM,Cam]=comstr(CAM,1);
elseif ~isempty(i2)
 if ~sp_util('issdt'); error('Command EltId-model requires SDT');
 elseif ~isstruct(model); error('Command EltId-model requires a model as input')
 else; RunOpt.Back='model'; RunOpt.EltId0=zeros(size(elt,1),1);
 end
else; RunOpt.Back='eltid';
end
i1=strfind(Cam,'fix');
if ~isempty(i1); 
 RunOpt.Fix=1;CAM(i1+[0:2])='';[CAM,Cam]=comstr(CAM,1);
else; RunOpt.Fix=0;
end
i1=strfind(Cam,'-skipcheck');
if ~isempty(i1); 
 RunOpt.skipCheck=1;CAM(i1+[0:9])='';[CAM,Cam]=comstr(CAM,1);
else; RunOpt.skipCheck=0;
end

if comstr(Cam,'from')
 i1=comstr(Cam(5:end),[-1]); if isempty(i1); error('Not a valid EltIdFrom');end
 RunOpt.Mode='Renumber'; RunOpt.From=i1;  RunOpt.EltId=[]; RunOpt.Fix=1;

elseif carg<=nargin; % return a string list by EltId
  RunOpt.EltId= varargin{carg};carg=carg+1;
  i3=find(RunOpt.EltId);cind(RunOpt.EltId(i3))=i3; st='';
  RunOpt.Mode='list';
%elseif ~isempty(strfind(Cam,'face')); RunOpt.EltId=[];st='';RunOpt.Mode='face';
else; RunOpt.EltId=[];st='';RunOpt.Mode='nominal';
end
i1=strfind(Cam,'silent');
if ~isempty(i1)&&~isempty(Cam)
 Cam(i1+[0:5])='';CAM(i1+[0:5])='';if ~sp_util('diag');RunOpt.Silent=1;end
elseif any(Cam==';')&&~sp_util('diag');
 Cam(Cam==';')='';CAM(CAM==';')='';if ~sp_util('diag');RunOpt.Silent=1;end
else; RunOpt.Silent=0;
end

[EGroup,nGroup]=getegroup(elt); eltrows=isfinite(elt(:,1));
out=zeros(size(elt,1),1);out1=out+0;
for jGroup=1:nGroup 
   [ElemF,opt]= feutil('getelemf',elt(EGroup(jGroup),:),jGroup); 
   cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1; if isempty(cEGI); continue; end
   [i2,i3]=fe_super('prop',ElemF,model);
   if i3==1 % Unique superelement
    i2=min(find(elt(EGroup(jGroup),:)==0));  %#ok<MXFND>
    if ~isempty(i2)&&size(elt,2)>=i2+2 
      out(EGroup(jGroup))=elt(EGroup(jGroup),i2+2);
    elseif isempty(i2); i2=size(elt,2);
    end
    out1(EGroup(jGroup))=i2+2;eltrows(EGroup(jGroup))=1;
   else % Generic element/superelement
    if length(i2)<3||~i2(3);   out1(cEGI)=0;out(cEGI)=0; 
    elseif i2(3)>size(elt,2); out1(cEGI)=i2(3);out(cEGI)=0; 
    else;out1(cEGI)=i2(3);  %out(cEGI)=elt(cEGI,i2(3));
     %sp_util('setinput',out1,i2(3)+0*elt(cEGI,i2(3)),EGroup(jGroup))
     sp_util('setinput',out,elt(cEGI,i2(3)),EGroup(jGroup))
     %sp_util('setinput',elt,out,[-3 i2(2)*size(elt,1)+EGroup(jGroup)-1 1]); % need stop
     %if strcmp(RunOpt.Mode,'face');out(cEGI)=elt(cEGI,i2(3)+[0 1])*[1;.1];
     %elseif strcmp(RunOpt.Mode,'edge');out(cEGI)=elt(cEGI,i2(3)+[0 1])*[1;.01];
     %end
    end 
    if isempty(RunOpt.EltId) 
    elseif strcmp(RunOpt.Back,'elt')||strcmp(RunOpt.Back,'model') % setting the modes
     if strcmp(RunOpt.Back,'model'); RunOpt.EltId0(cEGI)=elt(cEGI,i2(3)); end
     elt(cEGI,i2(3))=RunOpt.EltId(cEGI);
    elseif strcmp(RunOpt.Mode,'list')
      i2=out(cEGI); i3=find(i2); 
      if max(i2(i3))>length(cind); cind(max(i2(i3)))=0;end 
      if any(cind(i2(i3))) 
        st=sprintf('%s\n%s (%i)',st,ElemF,opt(1)); 
        for j1=find(cind(i2(i3)))' 
          st=sprintf('%s\n%s',st,sprintf('%i ',elt(cEGI(i3(j1)),:))); 
        end 
      end 
    end 
   end % Superelement type 
end 
if ~isempty(RunOpt.EltId)
 if strcmp(RunOpt.Back,'elt'); out=elt; % default for SDT and OpenFEM
 elseif strcmp(RunOpt.Back,'model')
  model.Elt=elt;r1=[];
  if isfield(model,'Stack')&&~isempty(model.Stack) % sdt specific call
   try
    RA=struct('Silent',0);
    i1=isfinite(elt(:,1)); r1=RunOpt.EltId0; i2=r1==0; r1(i1&i2)=max(r1)+1;
    eltind=sparse(r1(i1),1,RunOpt.EltId(i1));
    nind=sparse(model.Node(:,1),1,model.Node(:,1));
    model=feval(feutilb('@renumber_stack'),model,nind,model,eltind,RA);
   catch; sdtw('_nb','EltId-Elt failed stack renumber attempt');
   end
  end
  if isfield(model,'nmap')&&~isempty(r1);
   i1=isfinite(elt(:,1));
   model.nmap=vhandle.nmap.renumber(model.nmap,'Elts',[r1(i1) RunOpt.EltId(i1)]);
  end
  out=model;
 end
 return;
end
%i4=find(sparse(out+1,1,(out~=0)*1.1)>1.1)-1;
out2=[];
if ~RunOpt.Fix&&strcmp(RunOpt.Mode,'nominal')&&~RunOpt.skipCheck
 
 if RunOpt.Silent||Silent % In fact we dont care if in silent mode...
  %i4=find(sparse(out+1,1,(out~=0)*1.1)>1.1)-1;
  if length(unique(out(eltrows)))<sum(eltrows) %~isempty(i4)
   %ind=1:length(i4);
   %sdtw('_nb','Repeated EltId');
  else
   r5=out; r5(EGroup(1:nGroup))=max(r5)+(1:nGroup);
   [i4,i5]=sort(r5); ind=i5(find(~diff(i4))+1);
   in1=out(ind)~=0; % dissmiss zero values, as maybe duplicated but with no use
   if any(in1)
    ind=ind(in1);
    if length(ind)>100; ind=ind([1:10 end+[-9:0]]);
     fprintf('Displaying first 10 and last 10 only\n');
    end
    for j1=ind(:)'
     fprintf('EltInd %i ( %s)\n',j1,sprintf('%i ',out(j1))); %find(out==i4(j1))))
    end
    warning(['Repeated EltId ' sprintf(' %i',out(ind))]); %i4(ind))]);
   end
  end
 end
 
elseif RunOpt.Fix % fix

% i2=find(RunOpt.EltId); RunOpt.EltId(i2)=i1+[0:length(i2)-1]';
  if strcmp(RunOpt.Mode,'Renumber');
   i4=find(out); out(i4)=RunOpt.From+[0:length(i4)-1]'; out2=i4;
  %elseif ~isempty(find(sparse(out+1,1,(out~=0)*1.1)>1.1)-1) % was i4
  elseif length(unique(out(eltrows)))<sum(eltrows) % 4x faster (264509 elts)
   if ~RunOpt.Silent&&~Silent; 
       fprintf('Fixed repeated eltid\n');
   end
   [i5,i4]=sort(out); i1=find(diff(i5)==0&i5(1:length(i5)-1)~=0)+1;
   out(i4(i1))=max(out)+[1:length(i1)]'; out2=[out2;i4(i1)];
  end
  i1=find(out==0&eltrows); % fix zero values
  if ~isempty(i1)
    if ~RunOpt.Silent&&~Silent;
        fprintf('Fixed zero eltid\n');
    end
    out(i1)=max(out)+[1:length(i1)]'; out2=[out2;i1];
  end

  if ~isempty(out2) % some eltid have been set
   i2=find(out1); i3=i2+size(elt,1)*(out1(i2)-1);
   if numel(elt)<max(i3); elt(1,ceil(max(i3)/size(elt,1)))=0; end
   elt(i3)=out(i2);
  end
  out1=elt;
  
else; out1=elt; 
end
if nargout==0 && ~isempty(st); disp(st);
elseif strcmp(RunOpt.Back,'elt'); out=elt;
end

elseif comstr(Cam,'eltsetreplace')
 %% #EltSetReplace
 % renew sets in a meta-set if replaced elements
 model=varargin{carg}; carg=carg+1; % model
 r1=varargin{carg}; carg=carg+1; % meta-set from old model
 RO=varargin{carg}; carg=carg+1; % EEIdx, keepSets
 if ~isstruct(RO); RO=struct('EEIdx',RO); end
 
 % Extend EltId, NodeId
 eid=feutil('EltId',model);
 %nid=model.Node(:,1);
 
 r2=setdiff(eid,RO.EEIdx(:,2)); r2(r2==0)=[];
 if ~isempty(r2); EEIdx=[r2(:,[1 1]);RO.EEIdx]; else; EEIdx=RO.EEIdx;end
 EEIdx(any(EEIdx==0,2),:)=[];
 EEId=sparse(r1.EltId+1,1,1:length(r1.EltId));
 r1.SConn=r1.SConn(EEId(EEIdx(:,1)+1),:);
 r1.EltId=EEIdx(:,2);
 
 % rebuild NConn for EltSets
 % For NodeId sets, do nothing as link to element replacement is ambiguous
 
 out1=feutil('rmfield',r1,'NConn','NodeId');
 out=stack_set(model,'set','_gsel',out1);
 i1=~ismember(r1.SetNames(:,2),'NodeId');
 for j1=find(i1(:)')
  [out,r2]=stack_rm(out,'set',r1.SetNames{j1,1},'get');
  if ~isempty(r2)
   try
    if isfield(r2{1,3},'ID')&&~isempty(r2{1,3}.ID); stid=sprintf('-ID %i',r2{1,3}.ID); else; stid=''; end
    out=feutil(sprintf('AddSet%s %s',r1.SetNames{j1,2},stid),out,r1.SetNames{j1,1},...
     sprintf('setname _gsel:"%s"',r1.SetNames{j1,1}));
   catch; sdtw('_nb','ElSetReplace problem with %s: %s',r1.SetNames{j1,2},r1.SetNames{j1,1});
   end
  end
 end
 

%% #EGID --------------------------------------------------------------------2
elseif comstr(Cam,'egid')

 elt = varargin{carg};carg=carg+1;
 if isfield(elt,'Elt'); elt=elt.Elt; 
 elseif isfinite(elt(1))&&size(elt,2)==7; elt = varargin{carg};carg=carg+1;end

 if carg<=nargin; opt=varargin{carg}; carg=carg+1; else;opt=[];end
 % Specify name transcription table
 if ~isempty(strfind(Cam,'-vol'));
   RunOpt.NameTable={'quad4','q4p';'tria3','t3p';'quadb','q8p';'tria6','t6p'};
 else;RunOpt.NameTable={};
 end
 if carg<=nargin; out1=varargin{carg}; carg=carg+1; % eltname
 else; out1={};
 end 
 [EGroup,nGroup]=getegroup(elt);
 if ~isempty(opt)&&length(opt)~=nGroup
  error('you must provide nGroup EGID values to set');
 end
 if ~isempty(out1)&&length(out1)~=nGroup
  error('you must provide nGroup EltName values to set');
 end

 out=zeros(nGroup,1);
 for jGroup=1:nGroup
   [ElemF,i1]= feutil('getelemf',elt(EGroup(jGroup),:),jGroup);
   if ~isempty(RunOpt.NameTable)
    i2=strmatch(ElemF,RunOpt.NameTable(:,1),'exact');
    if ~isempty(i2);
     out1{jGroup}=RunOpt.NameTable{i2,2};
     elt(EGroup(jGroup),1:length(out1{jGroup})+2+length(i1))= ...
       [Inf abs(out1{jGroup}) 0 i1(:)'];
    end
   end
   if isempty(i1); i1=jGroup;end; out(jGroup)=i1(1); 
   if length(out1)<jGroup||isempty(out1{jGroup})
    out1{jGroup}=ElemF;
   end
   if ~isempty(opt);
     elt(EGroup(jGroup),1:length(out1{jGroup})+3)= ...
      [Inf abs(out1{jGroup}) 0 opt(jGroup)];
   end
 end
 if ~isempty(opt)||~isempty(RunOpt.NameTable); out=elt;end


%% #MatId  #ProId #MPID ------------------------------------------------------2
elseif comstr(Cam,'matid')||comstr(Cam,'proid')||comstr(Cam,'mpid')

 i1=strfind(Cam,'force');
 if ~isempty(i1); RO.force=1; CAM(i1+0:5)=''; Cam=lower(CAM);
 else; RO.force=0;
 end
 model=[];elt = varargin{carg};carg=carg+1;
 if isempty(elt)
 elseif isfield(elt,'Elt'); model=elt; elt=elt.Elt;  
 elseif isfinite(elt(1))&&size(elt,2)==7; elt = varargin{carg};carg=carg+1;
 end

 if ~isempty(strfind(Cam,'new'))&&~isempty(model)&&... % with new: assess content of changes
   isfield(model,'nmap')&&isKey(model.nmap,'Map:MChange')
  try  % add elements present in change cards to avoid missing ids
   mchM=model.nmap('Map:MChange');
   keys=mchM.keys; elc=[];
   for j1=1:length(keys)
    r1=mchM(keys{j1});
    if isequal(r1.type,'elt');elc=feutil('AddElt',elc,r1.Elt); end
   end
   if ~isempty(elc); elt=feutil('addelt',elt,elc); end
  catch; warning('MChange handling failed');
  end
 end
 
 [EGroup,nGroup]=getegroup(elt);
 if carg<=nargin; r1=varargin{carg};carg=carg+1;out=elt;
 else;r1=[];out=zeros(size(elt,1),3); end

 if ~isempty(model)&&~isempty(r1)&&(comstr(Cam,'matid')||comstr(Cam,'proid'))
   % renumber all matid proid, including pl il stack  mat pro
   try; model=feutilb(Cam,model,r1);
   catch; sdtw('_nb','%s renumbering failed',Cam); return;
   end
   out=model;return
 end
 
 for jGroup=1:nGroup
   cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
   [ElemF,opt]= feutil('getelemf',elt(EGroup(jGroup),:),jGroup);
   [i2,i3]=fe_super('prop',ElemF); i2(4)=0;
   if i3~=2 % single superelement
    out(EGroup(jGroup),3)=jGroup;
   elseif size(r1,2)>=2 %MPID [Matid Proid GroupId] ~isempty(r1) % set value
    if i2(1)&&(RO.force||any(r1(cEGI,1))); out(cEGI,i2(1))=r1(cEGI,1); end
    if i2(2)&&(RO.force||any(r1(cEGI,2))); out(cEGI,i2(2))=r1(cEGI,2); end    
   else
    if i2(1)&&size(elt,2)>=i2(1); out(cEGI,1)=elt(cEGI,i2(1)); end
    if i2(2)&&size(elt,2)>=i2(2); out(cEGI,2)=elt(cEGI,i2(2)); end
    out(cEGI,3)=jGroup;
   end
 end
 st='';
 if ~isempty(r1)
 elseif comstr(Cam,'matid'); [CAM,Cam]=comstr(CAM,6);out=out(:,1); st='pl'; 
 elseif comstr(Cam,'proid'); [CAM,Cam]=comstr(CAM,6);out=out(:,2); st='il';
 end
 if comstr(Cam,'new')&&~isempty(st);
  if isempty(out); out=1; else; out=max(out)+1; end
  if isstruct(model)
   r1=fe_mat(sprintf('get%s',st),model); 
   if isempty(r1); r1=1; else; r1=max(r1(:,1))+1; end
  else; r1=1;
  end
  out=max(out,r1);
 end

%% #Dof2Mode -----------------------------------------------------------------
elseif comstr(Cam,'dof2mode');  [CAM,Cam]=comstr(CAM,5);

 def= varargin{carg};carg=carg+1;
 mdof= varargin{carg};carg=carg+1;
 % NodeId ux uz uz rx ry rz
 i1=find(sparse(round(mdof),1,1));
 r2=fe_c(mdof,[i1+.01;i1+.02;i1+.03],'place');
 r3=fe_c(mdof,[i1+.04;i1+.05;i1+.06],'place');

 out=zeros(length(i1),4,size(def,2));
 if nnz(r3); r2=[r2;r3];i2=6; else;i2=3; end

 for j2=1:size(def,2)
  out(:,1,j2)=i1; 
  out(:,2:i2+1,j2)=reshape([r2]*full(def(:,j2)),length(i1),i2);
 end

%% #Mode2DOF -----------------------------------------------------------------
elseif comstr(Cam,'mode2dof');  [CAM,Cam]=comstr(CAM,9);

mode=varargin{carg}; carg=carg+1;
if ~isa(mode,'cell');mode={mode};end 
if isa(mode,'cell')
 def=[]; dof=[]; 
 r1=mode{1};
 RunOpt.Ndef=length(mode);if length(size(r1))>2;RunOpt.Ndef=size(r1,3);end
 i1=r1(:,1)';
 if size(r1,2)==4; dof=[i1+.01;i1+.02;i1+.03];
 else;dof=[1:size(r1,2)-1]'/100;
      dof=dof(:,ones(size(i1)))+i1(ones(size(dof)),:);
 end
 dof=dof(:); i1=i1(:); def=zeros(length(dof),length(mode));
 def(:,1)=reshape(r1(:,2:end,1)',length(dof),1);
 for j1=2:RunOpt.Ndef
    if j1>length(mode);r1=mode{1}(:,:,j1);else;r1=mode{j1};end
    if norm(r1(:,1)-i1); error('Variable node not supported'); end
    def(:,j1)=reshape(r1(:,2:end)',length(dof),1);
 end
end
if ~isempty(strfind(Cam,'struct'));out=struct('def',def,'DOF',dof);
else; out=def; out1=dof; 
end

%% #Divide -------------------------------------------------------------------
elseif comstr(Cam,'divide');  [CAM,Cam]=comstr(CAM,7);

 model=[];
 [carg,FEnode,FEelt,FEel0,ModelStack]=get_nodeelt(varargin,carg,ModelStack);
%% #DivideInParts - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if comstr(Cam,'inparts');
  mo1=struct('Node',FEnode,'Elt',FEelt);
  mo1.Elt=feutil('selelt seledgeall',feutil('quad2lin',mo1));
  mo1.Elt=feutil('divideingroups',mo1);
  [EGroup,nGroup]=getegroup(mo1.Elt);
  out=cell(1,nGroup); 
  for jGroup=1:nGroup; 
      i1=mo1.Elt(EGroup(jGroup)+1:EGroup(jGroup+1)-1,1:2);
      out{jGroup}=feutil('selelt withnode',model,i1(:));
  end
 elseif comstr(Cam,'in');[CAM,Cam]=comstr(CAM,9);
%% #DivideInGroups - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  [EGroup,nGroup]=getegroup(FEelt);
  NNode=sparse(FEnode(:,1)+1,1,1:size(FEnode,1));
  elt=[];RunOpt.Connected=[];
  if ~isempty(CAM); % connected to
   RunOpt.Connected=comstr(CAM,-1);RunOpt.nGroup=0; RunOpt.indOut=[];
  end 

  for jGroup = 1:nGroup %loop on element groups
   ElemF= feutil('getelemf',FEelt(EGroup(jGroup),:),jGroup);
   cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;
   if strncmpi(ElemF,'SE',2)
    r1=ones(length(cEGI),1)*FEelt(EGroup(jGroup),:);
    r1=[r1;FEelt(cEGI,:)]; r2=length(cEGI); %#ok<AGROW>
    r1=r1( reshape([1:r2;r2+1:2*r2],[],1),:);
    elt=[elt;r1]; %#ok<AGROW>
    
   else
    i1=fe_super('node',ElemF); % node indices
    i2 = [length(i1) length(cEGI) length(i1)*length(cEGI)];
    i3 = 1:i2(2);
    % build matrix of node use for elements
    
    i2=[2:length(i1) 1];
    in1=reshape(full(NNode(FEelt(cEGI,i1)+1)),[],1);
    in2=reshape(full(NNode(FEelt(cEGI,i2)+1)),[],1);
    i2=(in1&in2);
    i2 = sparse(in1(i2),in2(i2),1,max(in1),max(in2));
    i2=triu(i2+i2'); i3=etree(i2);
    i4=find(~i3); i3(i4)=-[1:length(i4)]; % no leaf
    i4=find(i3>0);while ~isempty(i4); i3(i4)=i3(i3(i4));i4=find(i3>0);end;
    i3=-i3;i4=i3(NNode(FEelt(cEGI,i1(1))+1));
    
    for j1=1:max(i4)
     i5=find(i4==j1);
     if isempty(i5)
     elseif ~isempty(RunOpt.Connected)
      if ~isempty(intersect(reshape(FEelt(cEGI(i5),i1),[],1),RunOpt.Connected))
       elt=[elt;FEelt(EGroup(jGroup),:);FEelt(cEGI(i5),:)]; %#ok<AGROW>
       RunOpt.indOut=[RunOpt.indOut;reshape(cEGI(i5),[],1)];
      end
     else
      elt=[elt;FEelt(EGroup(jGroup),:);FEelt(cEGI(i5),:)]; %#ok<AGROW>
     end
    end
   end % ElemF
  end % of loop on groups
  if ~isempty(strfind(Cam,'sort')) % Sort by group size
    [EGroup,nGroup]=getegroup(elt);  
    [i1,i2]=sort(diff(EGroup),'descend');ind=cell(nGroup,1);
    for jGroup=1:nGroup;
        ind{jGroup}=EGroup(i2(jGroup)):EGroup(i2(jGroup)+1)-1;
    end
    elt=elt(horzcat(ind{:}),:);
  end
  out=elt; if ~isempty(RunOpt.Connected); out1=RunOpt.indOut; end
  
%% #DivideGroup  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'group')  % 'DivideGroups'  - - - - - - - - - - - 

  [opt,CAM,Cam]=comstr(CAM,'group','%i');
  if length(opt)~=1
   error('DivideGroup can only be applied to a single group');
  end

  [EGroup,nGroup]=getegroup(FEelt);
  if opt>nGroup || opt<1; error('Not a valid group number');end
  if carg<=nargin; r1 = varargin{carg};carg=carg+1; else;r1=[];end
  [i1,elt]=feutil(['findelt' CAM],FEnode, ...
     FEelt(EGroup(opt):EGroup(opt+1)-1,:),FEel0,varargin{carg:end});
  if ~isempty(i1) && ~(length(i1)==EGroup(opt+1)-EGroup(opt))
   i2 = [EGroup(opt)+1:EGroup(opt+1)-1]';i2(i1-1)=0;i2=i2(i2~=0);
   if isempty(i2); sdtw('%s selects all elements, no division',CAM)
   else;i1 = [[1:EGroup(opt)]';i1+EGroup(opt)-1;EGroup(opt);i2; ...
        [EGroup(opt+1):size(FEelt,1)]'];
      FEelt = FEelt(i1,:); 
   end
  elseif ~Silent; warning('selection lead to no division'); 
  end
  out=FEelt;

% #HexaSurf ->sdtweb fe_shapeoptim('refinehexasurf'); -2
%% #DivideHexaTransition - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'hexatransition')  % 'DivideHexaTransition'

  [EGroup,nGroup]=getegroup(FEelt);
  NNode=sparse(FEnode(:,1),1,1:size(FEnode,1));
  elt=[];

  i1=varargin{carg};carg=carg+1; 
  cEGI=feutil('findelt withnode',FEnode,FEelt,[],i1);
  % Reference division
  xi = [-1 -1 -1;1 -1 -1;1 1 -1;-1 1 -1;-1 -1 1;1 -1 1;1 1 1;-1 1 1;
        -1/3 -1/3 0; 1/3 -1/3 0; 1/3  1/3 0; -1/3  1/3 0; -1/3 -1/3 1;
        1/3 -1/3 1;   1/3  1/3 1; -1/3  1/3 1  ];
  i2=[1 5 13 9 2 6 14 10; 2 6 14 10 3 7 15 11; 3 7 15 11 4 8 16 12;
        4 8 16 12 1 5 13 9; 9:16; 1:4 9:12 ];

  r1=ones(size(xi,1)*length(cEGI),3); % new nodes
  i3=ones(size(i2,1)*length(cEGI),10); % new elts

  for jElt=1:length(cEGI)
    nodeE=FEnode(NNode(FEelt(cEGI(jElt),1:8)),5:7);
    r2=integrules('hexa8',[xi zeros(size(xi,1),1)],nodeE);
    r1(16*jElt+[-15:0],1:3)=r2.N*nodeE;
    i3(6*jElt+[-5:0],1:10)=[16*(jElt-1)+i2 ones(6,1)*FEelt(cEGI(jElt),9:10)];
  end % jElt
  [FEnode,ind]=feutil('addnode',FEnode,r1);
  elt(end+[1:size(i3,1)],1:size(i3,2))=[ind(i3(:,1:8)) i3(:,9:10)];


  FEelt(cEGI(2:end,:),:)=[];
  FEelt=FEelt([1:cEGI(1)-1 cEGI(1)*ones(1,size(elt,1)) cEGI(1)+1:size(FEelt,1)],:);
  FEelt(cEGI(1)+[0:size(elt,1)-1],1:10)=elt;


  out=FEnode;out1=FEelt;
 
%% #DivideElt - model=feutil('divideelt (div)',model,div1,div2,div3)
% Moving of existing divide command in femesh
 else
    
  if comstr(Cam,'elt');[CAM,Cam]=comstr(CAM,4);end
  carg=2;model=varargin{carg};carg=carg+1;
  FEnode=model.Node; FEel0=model.Elt;

  [CAM,Cam,RunOpt.NewQuadDivide]=comstr('-new',[-25 3],CAM,Cam);
  [CAM,Cam,RunOpt.disp]=comstr('-disp',[-25 3],CAM,Cam); 
  [CAM,Cam,RunOpt.getInput]=comstr('-getinput',[-25 3],CAM,Cam); % output RO instead of performing cell call
  % get inline dv options
  if comstr(Cam,'{') % Accept {1:3,1} division
   opt=eval(Cam);
   if ~isempty(opt);dv1=opt{1}(:);end
   if length(opt)>1;dv2=opt{2}(:);end
   if length(opt)>2;dv3=opt{3};else;dv3=1;end
   opt=[];
  else;
   opt =comstr(CAM,[-1 1 1 1]);     opt=opt+1;
   % build dvi from further arguments
   dv1=linspace(0,1,opt(1))';dv2=linspace(0,1,opt(2))';dv3=linspace(0,1,opt(3))';
   i1=0;
   while carg<=nargin&&isa(varargin{carg},'double')
   r1=varargin{carg};carg=carg+1; i1=i1+1; if isempty(r1); r1=[0 1]; end
   if size(r1,1)>1
    r1=r1';
    if size(r1,2)>0; dv1=r1(1:max(find(r1(:,1))),1); else;dv1=[0;1];end
    if size(r1,2)>1; dv2=r1(1:max(find(r1(:,2))),2); else;dv2=[0;1];end
    if size(r1,2)>2; dv3=r1(1:max(find(r1(:,3))),3); else;dv3=[0;1];end
    warning('Use multiple arguments for each direction of irregular divide');
    break
   else
    if i1==1; dv1=r1(:); elseif i1==2; dv2=r1(:); elseif i1==3; dv3=r1(:); end
   end
   end
  end
  % loop on groups to get refinement if needed
  opt(1:3) = [length(dv1) length(dv2) length(dv3)];
  [EGroup,nGroup]=getegroup(FEel0);
  NNode=sparse(FEnode(:,1),1,1:size(FEnode,1));
  elt=[];node=FEnode;
  if carg>nargin||~isstruct(varargin{carg}); RO=struct;
   if carg<=nargin&&isequal(varargin{carg},';');  carg=carg+1;
    if ~sp_util('diag'); RunOpt.Silent=1; RunOpt.SiC=';'; end
   end
  else;RO=varargin{carg};carg=carg+1;end % (gvdr 15/05/2015)
  % implementation using sdtweb feutil('RefineCell')
  for jGroup = 1:nGroup %loop on element groups
    [ElemF,i1,ElemP]= feutil('getelemf',FEel0(EGroup(jGroup),:),jGroup);
    if ~isfield(RO,ElemP) % build refine cell input if necessary
     r2=isoCellRef(ElemP,{opt,dv1,dv2,dv3},RunOpt);
     if ~isempty(r2);RO.(ElemP)=r2;end
    end
  end
   
  if RunOpt.getInput; out=RO; return; % early output of RefineCell input (see ObjectDivide)
  else
    mo1=feutil('RefineCell',model,RO);
    if isfield(RO,'replace')&&RO.replace==2;out=mo1;
    else;     model.Node=mo1.Node; model.Elt=mo1.Elt;out=model;
    end
  end
   % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 end % 'divide... ' commands

%% #Find ---------------------------------------------------------------------
elseif comstr(Cam,'find');  [CAM,Cam]=comstr(CAM,5);

RunOpt.Silent=~isempty(Cam)&&strcmp(Cam(end),';');RunOpt.SiC='';
if RunOpt.Silent; CAM(end)='';Cam(end)='';
 if ~sp_util('diag'); RunOpt.SiC=';'; else; RunOpt.Silent=0; end
end

%% #FindNode - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% feutil('findnode',FEnode,FEelt,FEel0,varargin) % Cam,'findnode'
if comstr(Cam,'node'); [CAM,Cam]=comstr(CAM,5);
%i4 current node set, i5 operator positions
ModelStack={}; SetStack={};
[carg,FEnode,FEelt,FEel0,ModelStack,model]=get_nodeelt(varargin,carg,ModelStack);
[epsl,CAM,Cam]=test_epsl(epsl,CAM,Cam);
NNode=[];if ~isempty(FEnode);NNode=sparse(FEnode(:,1),1,1:size(FEnode,1));end
if comstr(Cam,'stack'); [CAM,Cam]=comstr(CAM,6);st2='return';else;st2=''; end
if carg<=nargin&&isa(varargin{carg},'cell'); Stack=varargin{carg};carg=carg+1;
elseif isempty(CAM)&&carg<=nargin&&isnumeric(varargin{carg})
 ind=varargin{carg};ind(ind>length(NNode))=[];
 out=full(NNode(ind)); out1=FEnode(out,:); return;
else;
 varg=varargin(carg:end); if isfield(model,'nmap'); varg{1,end+1}=model.nmap; end
 Stack=BuildSelStack(CAM,FEnode,FEelt,FEel0,1,varg{:});
end % of Stack is given as argument
if comstr(st2,'return'); out=Stack;return; end

if ~strcmp(Stack{end,1},'}'); Stack{end+1,1}='}';end
out=[]; j1=0;

%% loop on stack - - - - - - - - - - - - - - - - - - - - - - - - - - -
while j1 <size(Stack,1)-1

 j1=j1+1; Cam=comstr(Stack{j1,2},-27); Bole=Stack{j1,1};

 % finds nodes betwen the planes defined by two nodes
 if comstr(Cam,'between'); [st,Cam]=comstr(Cam,8);

  if length(Stack{j1,4})~=2; error('improper between entry'); end
  r1=FEnode(NNode(Stack{j1,4}),5:7);r1(2,:)=r1(2,:)-r1(1,:);
  r2=norm(r1(2,:)); if r2; r1(2,:)=r1(2,:)/r2;end
  r1=(FEnode(:,5:7)-r1(ones(size(FEnode,1),1),:))*r1(2,:)';
  i4=find(r1+epsl>=0&r1-epsl<=r2);

 % finds nodes in a group   - - - - - - - - - - - -
 elseif comstr(Cam,'group'); [st,Cam]=comstr(Cam,6);

 % model description matrix can be given as second argument to femesh
 if ~isempty(Cam)&&Cam(1)=='a'; elt=FEel0; else;elt=FEelt; end
 [EGroup,nGroup]=getegroup(elt);
 if strcmpi(Stack{j1,4},'all'); opt=1:nGroup;
 else;opt=IEqTest(1:nGroup,Stack{j1,3},Stack{j1,4}); end

 i1=[];
 for jGroup=opt; 
  if jGroup>0&&jGroup<=nGroup
   ElemF= feutil('getelemf',elt(EGroup(jGroup),:),jGroup);
   if strcmp(ElemF,'SE'); i2=[];SEopt=2;
    if ~RunOpt.Silent&&~Silent
     sdtw('_nb',['SE node selection not implemented, '...
      'nodes in SE will not be considered in the findnode result']);
    end
   else; [i2,SEopt]=fe_super('node',ElemF);
   end
   if SEopt(1,1)==2; i2 = elt(EGroup(jGroup)+1:EGroup(jGroup+1)-1,i2);end
   i1 =[i1;i2(:)]; %#ok<AGROW>
  end
 end
 i4 = find(sparse(i1+1,1,i1))-1;

 % GID Node group identification  - - - - - -
 elseif comstr(Cam,'gid'); [st,Cam]=comstr(Cam,4);

   opt=Stack{j1,4}; if ischar(opt); opt=str2num(opt);end %#ok<ST2NM>
   st=Stack{j1,3}; i1=FEnode(:,4);
   if comstr(st,'~='); i4=find(~ismember(i1,opt));
   elseif comstr(st,'>='); i4=find(i1>=opt(1));
   elseif comstr(st,'<='); i4=find(i1<=opt(1));
   elseif comstr(st,'>');  i4=find(i1>opt(1));
   elseif strcmp(st,'<');  i4=find(i1<opt(1));
   else 
    i4=find(ismember(i1,opt));
   end
   i4=FEnode(i4,1);

 % NodeId  - - - - - - - - - - - -
 elseif comstr(Cam,'nodeid'); [st,Cam]=comstr(Cam,7);

   opt=Stack{j1,4}; if ischar(opt); opt=str2num(opt);end %#ok<ST2NM>
   st=Stack{j1,3}; i1=FEnode(:,1);
   if comstr(st,'~=') 
     if max(opt)>length(NNode); NNode(max(opt),1)=0; end
     i4=full(NNode(opt));i4=i4(i4~=0); i1=1:size(FEnode,1);
     if ~isempty(i4) ; i1(i4)=0; end
     i4=find(i1);
   elseif comstr(st,'>='); i4=find(i1>=opt(1));
   elseif comstr(st,'<='); i4=find(i1<=opt(1));
   elseif comstr(st,'>');  i4=find(i1>opt(1));
   elseif strcmp(st,'<');  i4=find(i1<opt(1));
   else
     if max(opt)>length(NNode); NNode(max(opt),1)=0; end
     i4=full(NNode(opt));i4=i4(i4~=0);
   end
   i4=FEnode(i4,1);

 % Egid  - - - - - - - - - - - -
 elseif comstr(Cam,'egid')||comstr(Cam,'eltid')||comstr(Cam,'eltname') ...
        || comstr(Cam,'matid') || comstr(Cam,'proid')
   [r1,elt]=feutil('findelt',FEnode,FEelt,FEel0,Stack(j1,:));
   if isempty(elt);i4=[];
     if ~RunOpt.Silent&&~Silent;sdtw('No element selected for %s',Cam); end     
   else 
     i4=feutil('findnode',FEnode,elt,[],{'','Group','==','All'});
   end
 
 % inelt  - - - - - - - - - - - -
 elseif comstr(Cam,'inel') 

  ind=SubStackIndex(Stack,j1+1);

  if  Cam(5)=='t'
   mo1=struct('Node',FEnode,'Elt',FEelt,'El0',FEel0,'Stack',[]);
   mo1.Stack=ModelStack; if isfield(model,'nmap'); mo1.nmap=model.nmap; end
   [r1,elt]=feutil(['findelt' RunOpt.SiC],mo1,Stack(ind,:));
  elseif  Cam(5)=='0'
   mo1=struct('Node',FEnode,'Elt',FEel0,'El0',[],'Stack',[]);
   mo1.Stack=ModelStack;if isfield(model,'nmap'); mo1.nmap=model.nmap; end
   [r1,elt]=feutil('findelt',mo1,Stack(ind,:));
  else;error('Not a valid InElt command')
  end
  if isempty(elt); i4=[];     
      if ~RunOpt.Silent&&~Silent;sdtw('No element selected for %s',Cam); end     
  else;  i4=feutil('findnode',FEnode,elt,[],{'','Group','==','All'});
  end
  j1=max(ind);

 % notin - - - - - -
 elseif comstr(Cam,'notin') 

  ind=SubStackIndex(Stack,j1+1);
  mo1=struct('Node',FEnode,'Elt',FEelt,'El0',FEel0,'Stack',[]);
  mo1.Stack=ModelStack;if isfield(model,'nmap'); mo1.nmap=model.nmap; end
  [r1,elt]=feutil('findelt',mo1,Stack(ind,:));
  %[r1,elt]=feutil('findelt',FEnode,FEelt,FEel0,Stack(ind,:));
  if isempty(elt); i4=[];     
      if ~RunOpt.Silent&&~Silent;sdtw('No element selected for %s',Cam); end     
  else;i4=feutil('findnode',FEnode,elt,[],{'','Group','==','All'});
  end
  in1=FEnode(:,1);
  %i4=feutil('findnode',FEnode,elt,[],{'','Group','==','All'});
  in1(full(NNode(i4)))=0; i4=in1(in1~=0);clear in1;
  j1=max(ind);

 % finds nodes in a cylinder r x_0 y_0 z_0 nx ny nz - - - - - -
 elseif comstr(Cam,'cyl');
   st=Stack{j1,3};opt=Stack{j1,4};
   if ischar(opt); opt=str2num(opt);end %#ok<ST2NM>
   opt(5:7)=opt(5:7)/norm(opt(5:7)); % normalize vector
   r1=(FEnode(:,5:7)-opt(ones(size(FEnode,1),1),2:4))* ...
       sp_util('basis',opt(5:7),[0 0 0]);
   i3=sum((r1(:,2:3)).^2,2)-opt(1)^2;

   if     comstr(st,'>='); i4=find(i3>-epsl);
   elseif comstr(st,'<='); i4=find(i3<epsl);
   elseif comstr(st,'>');  i4=find(i3>0);
   elseif strcmp(st,'<');  i4=find(i3<0);
   else;                  i4=find(abs(i3)<epsl);
   end
   if length(opt)>7; i4(r1(i4,1)<opt(8))=[];end
   if length(opt)>8; i4(r1(i4,1)>opt(9))=[];end
   
   if ~isempty(i4); i4=FEnode(i4,1); end
     
 % finds nodes on a plane 
 elseif comstr(Cam,'plane'); [st,Cam]=comstr(Cam,6);% - - - - - -

   st=Stack{j1,3};opt=Stack{j1,4};
   if ischar(opt); opt=str2num(opt);end %#ok<ST2NM>
   i3=(FEnode(:,5:7)-opt(ones(size(FEnode,1),1),1:3))*opt(4:6)';

   if     comstr(st,'>='); i4=find(i3>-epsl);
   elseif comstr(st,'<='); i4=find(i3<epsl);
   elseif comstr(st,'>');  i4=find(i3>0);
   elseif strcmp(st,'<');  i4=find(i3<0);
   else;                  i4=find(abs(i3)<epsl);
   end
   if ~isempty(i4); i4=FEnode(i4,1); end

 % finds nodes in a set defined in the model stack
 elseif comstr(Cam,'set'); [st,Cam]=comstr(Cam,4); st='';
  try
   i4=FeutilMatchSet(ModelStack,Stack(j1,:));
   if isstruct(i4);sdtw('_nb','Assuming ''%s'' NodeId',Stack{j1,4});i4=i4.data;end
  catch err
   if comstr(Cam,'f'); i4=[]; % SetF: setnamesafe 
   else; err.rethrow;
   end
  end

 % finds nodes within a radius
 elseif comstr(Cam,'rad'); [st,Cam]=comstr(Cam,4); st='';

   st=Stack{j1,3};opt=Stack{j1,4};
   if ischar(opt); opt=str2num(opt);end %#ok<ST2NM>
   i2=length(opt);

   if i2<2; error('feutil findnode rad: must specify radius and node');
   elseif i2==2
    i2=find(FEnode(:,1)==opt(2));
    if isempty(i2); opt(2:4)=[0 0 0];else;opt(2:4)=FEnode(i2(1),5:7);end
   elseif i2<4; opt(4)=0;
   end
   i3=FEnode(:,5:7)-opt(ones(size(FEnode,1),1),2:4);
   i3=sum((i3').^2)-opt(1)^2;

   if     comstr(st,'=='); i4=find(abs(i3)<epsl);
   elseif comstr(st,'>='); i4=find(i3>=-epsl);
   elseif comstr(st,'>');  i4=find(i3>0);
   elseif strcmp(st,'<');  i4=find(i3<0);
   else;                  i4=find(i3<=epsl);
   end
   if ~isempty(i4); i4=FEnode(i4,1); end

 % finds nodes by coordinates
 elseif ~isempty(Cam) && any(Cam(1)=='xyzr')
   st=Stack{j1,3};opt=Stack{j1,4};
   if ischar(opt)&&isempty(intersect(opt,'xyzr')); opt=str2num(opt);end %#ok<ST2NM>
   i2=find('xyzr'==Cam(1));
   if i2<4; r1=FEnode(:,4+i2);else; r1=sqrt(sum(FEnode(:,5:6).^2,2));end
   if ischar(opt);x=FEnode(:,5);y=FEnode(:,6);z=FEnode(:,7);
     r=sqrt(sum(FEnode(:,5:6).^2,2));
     if isequal(st,'==')
         i4=find(abs(eval(horzcat(Cam,'-(',opt,')')))<epsl);
     elseif isequal(st,'~=')
         i4=find(abs(eval(horzcat(Cam,'-(',opt,')')))>epsl);
     elseif isequal(st,'>=')
         i4=find(eval(horzcat(Cam,'-(',opt,')'))>epsl);
     elseif isequal(st,'<=')
         i4=find(eval(horzcat(Cam,'-(',opt,')'))<-epsl);
     else;i4=find(eval(horzcat(Cam,st,opt)));
     end
   elseif comstr(st,'=='); i4=find(abs(r1-opt)<epsl);
   elseif comstr(st,'~='); i4=find(abs(r1-opt)>epsl);
   elseif comstr(st,'>='); i4=find(r1-opt>=-epsl);
   elseif comstr(st,'>');  i4=find(r1-opt>0);
   elseif strcmp(st,'<');  i4=find(r1-opt<0);
   else;                  i4=find(r1-opt<=epsl);
   end
   if ~isempty(i4); i4=FEnode(i4,1); end

 elseif strcmp(Stack{j1,1},'}')
 elseif comstr(Cam,'name');
   [n1,st3]=sdth.urn('nmap.Node',model,Stack{j1,4});
   i4=n1(:,1);
 elseif comstr(Cam,'distfcn');
   % fecom('shownodemark','distfcn{sphere{.1125 .0375 .075, .02}}')
   r3=comstr(Stack{j1,4},1);if r3(8)=='"';r3([8 end])='';end % distfcn"{}"
   [~,r3]=sdtm.urnPar(r3,'{}{}');
   r3.distFcn=lsutil('gen',[],r3.Other(1));
   i4=FEnode(r3.distFcn(FEnode(:,5:7))>0,1);
 elseif comstr(Cam,'dist');
   % Use distance function Dist(SurfDist,4003874)<100
   n1=regexp(Stack{j1,3},'([^,]*),(.*)','tokens');n1=n1{1};
   r1=ModelStack{strcmpi(ModelStack(:,2),n1{1}),3};
   r1=r1.dist.idToDist(r1.dist.C1,str2double(n1{2}));
   i4=fix(r1.DOF(eval(sprintf('r1.def%s',Stack{j1,4}))));

 else % give coordinates  with possible infinity

   opt=Stack{j1,4}; if ischar(opt); opt=str2num(opt);end %#ok<ST2NM>
   i4=1:size(FEnode,1);
   if isfinite(opt(1,1)); i4=i4(abs(FEnode(i4,5)-opt(1))<epsl); end
   if ~isempty(i4)&& isfinite(opt(1,2))
         i4 = i4(abs(FEnode(i4,6)-opt(2))<epsl); end
   if ~isempty(i4)&& isfinite(opt(1,3))
         i4 = i4(abs(FEnode(i4,7)-opt(3))<epsl); end
   i4=FEnode(i4,1);

 end

 % boolean operations on node sets
 if isa(i4,'v_handle'); i4=i4.GetData; end
 if strcmp(Bole,'&') % keep those present in all
   if isempty(i4); out=i4;
   elseif isempty(out); 
   else
    i4=find(sparse(i4,ones(size(i4)),ones(size(i4))));
    i4=sparse(i4,1,i4);if max(out)>length(i4); i4(max(out))=0;end
    out = nonzeros(i4(out));
   end
 elseif strcmp(Bole,'&~'); out=setdiff(out,i4);
 elseif isempty(out); out=double(i4);
 elseif strcmp(Bole,'|') % keep those present in any
   out=[double(i4);out]; %#ok<AGROW>
   out=find(sparse(out,ones(size(out)),ones(size(out))));
 end
end % of j1 loop on arguments

if ~isempty(out)&&out(end)==size(FEnode,1)&&any(~isfinite(FEnode(end,5:7)))
  out=out(1:end-1);
end
if nargout==2 % return the nodes too
  NNode=sparse(FEnode(:,1),1,1:size(FEnode,1));
  out=out(out<=max(FEnode(:,1))); out1=[]; % avoid error for max dim
  if ~isempty(out) % be robust to setname calls with non existing nodes
   i1=NNode(out); i1(i1==0)=[]; 
   if ~isempty(i1); out1=FEnode(i1,:); end
  end
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'el'); [CAM,Cam]=comstr(CAM,4);
 % THIS IS AN AUTODOC TEST DO NOT ERASE, ASK GM
 %{
  #FUNREF
  {'cmd','FindElt'
   'purpose','Find list of elements from element selection string'
   'syntax','@STX.FindElt'
   'in1','@STX.cmd'
   'in2','@DATA.model'
   'out1','@DATA.EltInd'
   'out2','@DATA.Elt'
   'opt',{}
  }
 %}

%{
```FUNREF
feutil(findelt)                       % Find list of elements from element selection string
stx: @STX.FindElt
opt: epsl(%g)                         % Evaluation tolerance for equality logical operators
out: eltind@DATA.EltInd               % Indices of selected elements in the element description matrix
out: elt@DATA.Elt                     % Description matrix of selected element
inp: cmd@STX.CMD                      % SDT command+
inp: model@DATA.model                 % SDT model in which elements are seeked
```
%}

%{
```FUNREF -2
feutil(FindElt):  Find list of elements from element selection string
```MD
Any other description text xxx
```INP
cmd (STX.CMD+FindElt): string with 'findelt' followed by element selectors using @STX.FindElt format
model (DATA.model): SDT model in which elements are seeked
RO (DATA.RO?): Optional argument providing command options as name/value pairs using a MATLAB structure 
```OPT
epsl (1e-6%g): Evaluation tolerance for equality logical operators
```OUT
eltind (DATA.EltInd): Indices of selected elements in the element description matrix
elt (DATA.Elt): Description matrix of selected element
%}
% #FindElt  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [ind,elt,IndWithHeaders]=feutil('findelt',FEnode,FEelt,FEel0)


%i4 current element set, out final element set, i5 operator positions
%i7 used groups

ModelStack={};
% variable model typically created here
RunOpt.ImplicitNode=1;RunOpt.CheckedNodes=[0 0];
[carg,node,elt,el0,ModelStack,model]=get_nodeelt(varargin,carg,ModelStack,RunOpt);
[r1,CAM,Cam]=test_epsl(epsl,CAM,Cam); 
if r1~=epsl; RunOpt.epsl=sprintf('epsl %.15g',r1);epsl=r1;
else;RunOpt.epsl='';
end

if comstr(Cam,'stack');[CAM,Cam]=comstr(CAM,6);st2='return';else;st2=''; end
if isempty(elt)&&isempty(el0)
 Stack={'','Empty model'}; RunOpt.LastOp=''; nGroup=0;
elseif carg<=nargin&&isa(varargin{carg},'cell')
  Stack=varargin{carg};carg=carg+1;
  if size(Stack,1)==1&&any(strcmp(Stack{1},{'|','&','{'})); Stack{1}='';end
  RunOpt.LastOp='';  %[EGroup,nGroup]=getegroup(elt);
else % stack is not given
 varg=varargin(carg:end);
 if isfield(model,'nmap')&&(isempty(varg)||~isa(varg{end},'vhandle.nmap'))
  varg{1,end+1}=model.nmap;
 end
 if carg>nargin&&isempty(Cam); error('You must give a selection argument');
 elseif ~isempty(CAM)
   RunOpt.LastOp='';Stack=BuildSelStack(CAM,node,elt,el0,1,varg{:});%varargin{carg:end});
 elseif ischar(varargin{carg}); 
   [CAM,Cam]=comstr(varargin{carg},1); carg=carg+1;
   RunOpt.LastOp='';Stack=BuildSelStack(CAM,node,elt,el0,1,varg{:});%varargin{carg:end});
 else
   r1=varargin{carg};carg=carg+1;RunOpt.LastOp=''; 
   if isfield(r1,'type'); Stack={'{' 'Set',[],r1; '}' '' '' ''};
   else; error('Not a known entry selection');
   end
 end
end % is stack given
if strcmp(st2,'return'); out=Stack;return; end
if ~strcmp(Stack{end,1},'}') && ~strcmp(Stack{1,2},'Empty model')
  Stack{end+1,1}='}';
end
j1=0; out=[];RunOpt.Transformed=0; RunOpt.EltIdFixed=0;
%[EGroup,nGroup]=upEG(EGroup,nGroup,elt); % to update if needed
% RunOpt.Transform set to j1 if current operation transforms type
% no combination to previous EltInd is possible if elt has been transformed
% current elements are assumed to be the new selection result
EGroup=[]; nGroup=[];

while j1<size(Stack,1)-1 % loop on elt sel stack- -  - - - - - - - - - - - - - -

 j1=j1+1; Bole=Stack{j1,1}; 
 st=Stack{j1,2};opts=Stack{j1,3};opt=Stack{j1,4};

 % Node Based selection
 if ischar(st)&&~isempty(strfind(st,'Node'))
  if strcmp(Stack{j1,3},'{}')
    ind=SubStackIndex(Stack,j1+1);
    if ischar(node);eval(node);end
    mo1=struct('Node',node,'Elt',elt,'El0',el0,'Stack',[], ...
        'CheckedNodes',RunOpt.CheckedNodes);
    mo1.Stack=ModelStack; if isfield(model,'nmap'); mo1.nmap=model.nmap; end
    if strcmp(opts,'{}'); 
        opt=feutil(['findnode' RunOpt.epsl RunOpt.SiC],mo1,Stack(ind,:)); 
    end
    j1=max(ind);
  else; opt=Stack{j1,4};
  end
  if any(opt<0); NNode=(sparse(-opt+1,1,1:length(opt)));
  else; NNode=(sparse(opt+1,1,1:length(opt)));end
 end
 % reindexing for options

 i4=[];
 if comstr(lower(st),'seledge')
    [EGroup,nGroup,elt]=upEG(EGroup,nGroup,elt,j1,out);
    if ischar(node);eval(node);end
    [CAM,Cam,i1]=comstr('-lin',[-25 3],Stack{j1,4},lower(Stack{j1,4}));
    [elt,CAM]=feutil(['getedgeline' CAM ';'],node,elt);
    if i1;elt=feutil('quad2lin',node,elt);elt=elt.Elt;end
    if isempty(elt); out=[]; %EGroup=[];nGroup=[];
    else; out=find(isfinite(elt(:,1))); %[EGroup,nGroup]=getegroup(elt);
    end
    EGroup=[];nGroup=[];
    RunOpt.LastOp='edge';RunOpt.Transformed=j1;

 elseif comstr(lower(st),'$');eval(Stack{j1,4}(2:end));
 elseif comstr(lower(st),'set')

   %i4=(0:min(find(~strcmp(Stack(j1:end,2),'Set')))-2)+j1;j1=i4(end); %#ok<MXFND>
   try;  [i4,elt]=FeutilMatchSet(ModelStack,Stack(j1:end,:),elt,RunOpt);
   catch err % allow safe mode
    if comstr(lower(st),'setf'); i4=[]; % allow empty, % do not alter elt if we are combining
    if ismember(Bole,{'|','&','&~'}); else; elt=[];  end
    else; err.rethrow; 
    end
   end
   EGroup=[];nGroup=[]; %[EGroup,nGroup]=getegroup(elt);
   if isequal(i4,'faces')||isequal(i4,'edges'); % faces were selected    
    if isempty(elt);
     ste=Stack{j1,4}; if isnumeric(ste); ste=num2str(ste(:)'); end
     if isstruct(ste)&&isfield(ste,'name');ste=ste.name; end
     if ~RunOpt.Silent&&~Silent;sdtw('_nb','No face/edge was selected for ''set'',''%s''',ste);end
     out=[]; i4=[];
    else
     out=find(isfinite(elt(:,1))); i4=[];
    end
    RunOpt.LastOp='setface';RunOpt.Transformed=j1;
   elseif isequal(i4,0) % no elements ?
    sdtw('_nb','No element was selected for ''set'',''%s''',Stack{j1,4});
    i4=[];
   % else % elements were selected
   end

 elseif comstr(lower(st),'selface') 
    % #lowlevel_selface -3
    [EGroup,nGroup,elt]=upEG(EGroup,nGroup,elt,j1,out);
    if isempty(elt); %error('No element to select'); %end
     if comstr(lower(st),'selfacef'); out=[];
     else;error('No element to select');
     end
    else
    if ischar(node);eval(node);end
    [elt,CAM]=feutil('getedgepatch',node,elt,[],Stack{j1,4});
    EGroup=[];nGroup=[]; i4=[]; %[EGroup,nGroup]=getegroup(elt);
    if isempty(elt); out=[];else;out=find(isfinite(elt(:,1)));end
    end
    RunOpt.LastOp='edge';RunOpt.Transformed=j1;
    
 elseif comstr(lower(st),'connectedto') % #lowlevel_ConnectedTo NodeId -3
    %[EGroup,nGroup,elt]=upEG(EGroup,nGroup,elt,j1,out);
    out=sort(out); % out can be empty if transformation prior to connectedTo
    [u1,u2,el1]=upEG(EGroup,nGroup,elt,j1,out); i5=find(isfinite(el1(:,1)));
    %size(out), size(i5) i5 consistent with out
    if isempty(el1); error('No element to select'); end
    if ischar(node);eval(node);end
    [u1,i4]=feutil(['divideingroups' sprintf('%i ',Stack{j1,4})],node,el1);
    %EGroup=[];nGroup=[]; %[EGroup,nGroup]=getegroup(elt);
    if isempty(i4); %out=[]; 
       if ~RunOpt.Silent&&~Silent;sdtw('_nb','No connected element');end
    %else; %out=find(isfinite(elt(:,1))); No due to transformation !
    elseif ~isempty(out); i4=out(ismember(i5,i4)); % we have elt(out(ismember(i5,i4)),:) == el1(i4,:)
    end
    RunOpt.LastOp='connectedto';%RunOpt.Transformed=0;
 elseif isempty(st)&&any(strcmp(Stack{j1,1},{'&','|'})) % #lowlevel_ConnectedTo NodeId -3
    error('&& and || are not accepted in FindElt commands')
 else % standard element selection (loop on groups)

 %[EGroup,nGroup]=getegroup(elt);RunOpt.LastOp='';
 [EGroup,nGroup]=upEG(EGroup,nGroup,elt); RunOpt.LastOp='';
 eltid=[];
 for jGroup = 1:nGroup %loop on element groups
   [ElemF,EGID]= feutil('getelemf',elt(EGroup(jGroup),:),jGroup);EGID=EGID(1);
   [i2,i6]=fe_super('prop',ElemF);
   cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1; i3 = [];

   switch lower(st)
   case 'group' % group
      i3=IEqTest(jGroup,opts,opt);
      if ~isempty(i3)
          i3 = [1:EGroup(jGroup+1)-EGroup(jGroup)-1];
          if i6(1)==1; i4(end+1,1)=EGroup(jGroup); end %#ok<AGROW>
      end
   case 'egid' % egid
      i3=IEqTest(EGID(1),opts,opt);
      if ~isempty(i3)
          i3 = [1:EGroup(jGroup+1)-EGroup(jGroup)-1];
          if i6(1)==1; i4(end+1,1)=EGroup(jGroup); end %#ok<AGROW>
      end
   case 'eltind' % eltind 
       if i6~=2; i3=IEqTest(EGroup(jGroup),opts,opt); 
       else;i3=IEqTest(cEGI,opts,opt);end 
   case {'eltid','name'} % eltid
    if comstr(lower(st),'name')
     [i3,st3]=sdth.urn('nmap.Elts',model,Stack{j1,4});
     opt=i3{:};
    end
    if isempty(eltid); eltid=feutil('eltid;',elt); end
    if i6~=2; i3=IEqTest(eltid(EGroup(jGroup)),opts,opt);
    else;i3=IEqTest(eltid(cEGI),opts,opt); end
   case 'eltname' % eltname
      if comstr(lower(opt),'se:');
          RunOpt.SubSe=opt(4:end);opt='SE';st1=opts;opts='';
      elseif jGroup==1;RunOpt.SubSe='';
      end
      if strcmp(opts,'~=')
       if any(opt==' '); i7=isempty(strfind(opt,ElemF));else; i7=~comstr(ElemF,opt);end
       if i7
        i3 = [1:EGroup(jGroup+1)-EGroup(jGroup)-1];
        if i6(1)==1; i4(end+1,1)=EGroup(jGroup); end %#ok<AGROW>
       end
      elseif isempty(opts)
       if comstr(ElemF,opt)
        i3 = [1:EGroup(jGroup+1)-EGroup(jGroup)-1];
        if i6(1)==1; i4(end+1,1)=EGroup(jGroup); end %#ok<AGROW>
        if ~isempty(RunOpt.SubSe);
         if isequal(RunOpt.SubSe(1),'#') % allow regexp on name
          i7=~cellfun(@isempty,regexp(cellfun(@(x)fesuper('s_',x),...
           num2cell(elt(cEGI(i3),1)),'uni',0),RunOpt.SubSe(2:end),'once'));
         else; i7=elt(cEGI(i3),1)==fesuper(['s_',RunOpt.SubSe]);
         end
         switch st1
         case '~=';  i3(i7)=[];  %i3(elt(cEGI(i3),1)==fesuper(['s_',RunOpt.SubSe]))=[];
         otherwise;  i3(~i7)=[]; %i3(elt(cEGI(i3),1)~=fesuper(['s_',RunOpt.SubSe]))=[];
         end
        end
       elseif ~isempty(RunOpt.SubSe)&&isequal(st1,'~=')
        % accept all other elts if eltname ~=SE:sename
        i3 = [1:EGroup(jGroup+1)-EGroup(jGroup)-1];
        if i6(1)==1; i4(end+1,1)=EGroup(jGroup); end %#ok<AGROW>
       end
      end
      
   case 'facing' % facing Node angle
    i3=fe_super('face',ElemF);
    if length(opt)<4; opt(4)=0; end
    if size(i3,1)<=1 % Faces and beams
     if ischar(node);eval(node);end
     [nor,cg]=feutil('getnormal;',node, ...
      elt(EGroup(jGroup):EGroup(jGroup+1)-1,:));
     cg=opt(ones(size(cg,1),1),2:4)-cg;r3=sqrt(sum(cg.*cg,2));
     i3=find(any(nor,2));r3=sum(nor(i3,:).*cg(i3,:),2)./r3(i3);
     switch opts
     case '>';  i3=i3(r3>opt(1));
     case '>=';  i3=i3(r3>=opt(1));
     case '<';  i3=i3(r3<opt(1));
     case '<=';  i3=i3(r3<=opt(1));
     case '<>';  i3=i3(abs(r3)>opt(1));
     case '><';  i3=i3(abs(r3)<opt(1));
     otherwise; error('''%s'' not a supported opts value',opts);
     end
     i3=i3-1;
    else
     warning('Facing selection can only be applied to flat elements')
     i3=[];
    end

   case 'matid' % mat
     if i2(1)~=0
      if i2(1)>size(elt,2); i3=[]; % No matid is given
      else
        i2=elt(cEGI,i2(1)); i2(rem(i2,1)~=0)=0; % integer only
        i3=IEqTest(i2,opts,opt);
      end
     end
   case 'proid' % pro
     if i2(2)==0||i2(2)>size(elt,2);
     else;  
       i2=elt(cEGI,i2(2)); i2(rem(i2,1)~=0)=0; % integer only
       i3=IEqTest(i2,opts,opt);
     end
   % with/without/in node single superelement
   case {'withnode','withoutnode','innode'} 

    if i6(1)==1 % unique superelement
     i2=fe_super('node',ElemF);
     if max(max(i2))+1>length(NNode); NNode(max(max(i2))+1,1)=0;end
     i2 = reshape(full(NNode(i2+1)),size(i2,1),size(i2,2))';
     if     comstr(st,'WithNode')&&any(i2);  i4(end+1,1)=EGroup(jGroup);%#ok<AGROW>
     elseif comstr(st,'InNode')  &&all(i2);  i4(end+1,1)=EGroup(jGroup);%#ok<AGROW>
     elseif comstr(st,'Without') &&~any(i2); i4(end+1,1)=EGroup(jGroup);%#ok<AGROW>
     end
   % with/in node standard element
    elseif strcmp(ElemF,'SE');
     cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;i3=zeros(size(cEGI));
     for jElt=1:length(cEGI)
      eval('i2=fesuper(''serenumber -nodeid'',model,elt(cEGI(jElt),:));');
      if max(max(i2))+1>length(NNode); NNode(max(max(i2))+1,1)=0;end
      i2=NNode(i2+1);
      if     comstr(st,'WithNode')&& any(i2);  i3(jElt)=1;
      elseif comstr(st,'InNode')  && all(i2); i3(jElt)=1;
      elseif comstr(st,'Without') &&~any(i2); i3(jElt)=1;
      end
     end
     i3=find(i3);
    else
     i2 = elt(EGroup(jGroup)+1:EGroup(jGroup+1)-1,fe_super('node',ElemF));
     i3=i2==0; % Detect elements with node entries (ni) set to 0
     if max(max(i2))+1>length(NNode); NNode(max(max(i2))+1,1)=0;end
     i2 = reshape(full(NNode(i2+1)),size(i2,1),size(i2,2));
     i2(i3)=NaN; % make sure ni=0 does not interfere with detection
     if     comstr(st,'WithNode'); i3 = find(any(i2,2));  % with
     elseif comstr(st,'InNode');   i3 = find(all(i2,2));  % in
     elseif comstr(st,'Without');  i3 = find(~any(i2,2)); % without
     end
    end
   otherwise; error('''%s'' not a supported FindElt command',st);
   end
   if ~isempty(i3);  i4=[i4;cEGI(i3(:))']; end %#ok<AGROW>
 end% of jGroup loop
 end % or seledge

 % bolean operations on element sets
 if RunOpt.Transformed==j1  % specific handling:
  % EltInd (i4) is irrelevant, we assume elt and out have been handled
  %(strcmp(RunOpt.LastOp,'edge')||strcmp(RunOpt.LastOp,'setface'))
  if strcmp(Bole,'&~')
   error(['Boolean &~ cannot be used for element selection with '...
    'type transformation, use token :exclude instead']);
  end
 elseif strcmp(Bole,'&')   % keep those present in all
   if isempty(i4); out=i4;
   elseif isempty(out);
   else
    i4=find(sparse(i4,ones(size(i4)),ones(size(i4))));
    i4=sparse(i4,1,i4);if max(out)>length(i4); i4(max(out))=0;end
    out = nonzeros(i4(out));
   end
 elseif strcmp(Bole,'&~'); out=out(~ismember(out,i4)); %setdiff(out,i4);
 elseif isempty(out); out=i4;
 elseif strcmp(Bole,'|')   % keep those present in any
   %out=[i4;out];out=find(sparse(out,ones(size(out)),ones(size(out)))); %#ok<AGROW>
   out=[i4;out(~ismember(out,i4))];
 end

end % of j1 loop on Stack

if nargout>1 % return the elements too
 [EGroup,nGroup]=upEG(EGroup,nGroup,elt);
  i7=zeros(nGroup,1); % specific treatment at output
  for jGroup = 1:nGroup
   if any(out>EGroup(jGroup)&out<EGroup(jGroup+1)); i7(jGroup)=1;end
  end
  if RunOpt.Transformed
   out2=sort([out;EGroup(i7~=0)']); out1=elt(out2,:);
   out=[]; out2=[1:size(elt,1)]';
  elseif nGroup==0; out2=[]; out1=[];
  else; out2=sort([out;EGroup(i7~=0)']); out1=elt(out2,:);
  end
elseif RunOpt.Transformed&&~RunOpt.Silent&&~Silent
  sdtw('output element selection is not a sub-part of the model');
end

else;error('Find%s unknown',CAM);
end % subcommand selection - - - - - - - - - - - - - -

%% #FixMPCMaster -------------------------------------------------------------
% recombines a constraint equation so that MASTER DOFs can be defined
elseif comstr(Cam,'fixmpcmaster')

 c=varargin{carg};carg=carg+1; RO=struct('typ','mat','ok',0);
 if isstruct(c); RO.typ='struct';
 else
  c=struct('c',c); RO.typ='mat';
  if carg<=nargin;  % an initial list of masters is provided
   c.slave=varargin{carg}; carg=carg+1;
  end
 end
  
 if isfield(c,'slave') % slave is known, just check c
  out1=c.slave;out1=out1(:);
  [i1,i2]=find(c.c(:,out1)-speye(size(c.c,1)));
  if isempty(i1); out=c.c; RO.ok=1;% the slave set is consistent
  else
   [II,JJ,KK]=find(c.c(:,out1)); % gvdr: further checks to avoid slave ordering problems
   if length(KK)==size(c.c,1)&&all(KK==1) % maybe the identity is not well ordered
    [i2,i3]=sort(II); i3=out1(i3); % reorder slave set to get identity
    [i2,i4]=find(c.c(:,i3)-speye(size(c.c,1))); % re check new slave consistency
    if isempty(i2); out=c.c; out1=i3; RO.ok=1; end % only modify output if identity was found
   end % i1 is kept unchanged from original check
   clear II JJ KK
  end
  if ~RO.ok % need to clean up slave set
   i1=unique(i1);
   % masters that are not correct and dependent of incorrect ones
   ind=[];
   while ~isequal(i1,ind);
    ind=i1; i1=find(any(c.c(:,any(c.c(ind,:),1)),2));
   end
   % consider what is left separately
   [i2,i3]=setdiff(1:size(c.c,1),i1);
   out=struct('i0',i3,'slave0',out1(i3),'c0',c.c(i3,:));c.c=c.c(i1,:);
   try;
    out1(i3)=[];
    [l,u]=lu(c.c(:,out1),sdtdef('luThres'));
    c1=u\ ( l\c.c);
    if norm(c1(:,out1)-speye(length(out1)),'inf')>1e-6
     fprintf('Independence problem on MPC');error('Go to catch');
    end
    c1(:,out1)=speye(length(out1));
    out1=[out.slave0;out1(:)]; out=[out.c0;c1]; RO.ok=1; %return;
   catch; out1=[];
   end
  end
 else; out=[]; out1=[];
 end

 if ~RO.ok % no slave provided or could not exploit it
  ind=find(any(c.c,1)); r1=c.c(:,ind)';i2=1:size(r1,2);
  if nnz(r1(:,i2))==i2(end)&&norm(r1(i2,:)-eye(length(i2)),'inf')<1e-16
  else % mixed slaves, but need to be robust to master with coef>1
   [i3,i4]=find(r1==1);
   if length(unique(i4))==size(r1,2)&&length(i3)==size(r1,2)&&...
     norm(r1(i3,:)-eye(length(i4)),'inf')<1e-16;i2=i3;
   else;[r2,i2]=max(abs(r1),[],1); % check that there is no obvious solution
   end
  end
  if length(unique(i2))==size(r1,2) % all max independent, does not check whether only one max
   [l,u]=lu(r1(i2,:)); r1=(r1/u)/l;
   out1=ind(i2);
  else
   out1=zeros(size(c.c,1),1);
   for j1=1:size(c.c,1); % select master and eliminate from constraints
    [r2,i2]=max(abs(r1(:,j1))); out1(j1)=ind(i2); % was out1(end+1)
    r3=r1(i2,j1); if r3~=0; r1(:,j1)=r1(:,j1)/r3;end
    r2=r1(i2,:);r2(j1)=0;r1=r1-r1(:,j1)*r2;
    if ~ismember(out1(j1),out1(1:j1-1))
     r2=zeros(size(r2));r2(j1)=1;r1(i2,:)=r2;
    end
   end
  end
  c.c(:,ind)=r1';
  if ~isempty(out);
   out1=[out.slave0;out1(:)]; out=[out.c0;c.c];
  else; out=c.c;if ~isempty(out1); out1=out1(:);end
  end
  % remove unused constraints
  c.c=out';ind=any(c.c,1);if nnz(ind); out=c.c(:,ind)';out1=out1(ind);end
 end
 
 % output
 if isequal(RO.typ,'struct'); %out.slave=out1;
  out=struct('DOF',c.DOF,'c',out,'slave',out1);
 elseif isequal(RO.typ,'mat') % legacy behavior
 end
 
%% #Geom
elseif comstr(Cam,'geom');[CAM,Cam]=comstr(CAM,5);
    
%% #GeomDilat ----------------------------------------------------------------
if comstr(Cam,'dilat')

 node=varargin{carg};carg=carg+1;
 if isfield(node,'Node'); node=node.Node;end

 if nargin<carg; nu=0.3; else;nu=varargin{carg};carg=carg+1; end

 out=struct('DOF',[node(:,1)+0.01;node(:,1)+0.02;node(:,1)+0.03], ...
            'def',[node(:,5);-nu*node(:,6);-nu*node(:,7)],'lab',{{'Dilatation'}}  ); 
 
%% #GeomRB -------------------------------------------------------------------
elseif comstr(Cam,'rb')

 if sp_util('issdt')
  warning('You should rather use feutilb(''geomrb'')')
 end
 model=varargin{carg};carg=carg+1;
 if ~isfield(model,'Node'); error('You must provide a model for GEOMRB');end

 [node,bas,NNode]=feutil('getnodebas',model);
 i1=unique(node(:,1));i1=i1(:)'; 
 rdof=[i1+.01;i1+.02;i1+.03];rdof=rdof(:);

 r2=zeros(length(rdof),6); r2(1:3:end,1)=1; r2(2:3:end,2)=1;r2(3:3:end,3)=1;
 i1=fe_c(rdof,node(:,1)+.01,'ind');
 r2(i1+1,4)=-node(:,7);r2(i1+2,4)=node(:,6);
 r2(i1+2,5)=-node(:,5);r2(i1+0,5)=node(:,7);
 r2(i1+0,6)=-node(:,6);r2(i1+1,6)=node(:,5);
 out=struct('def',r2,'DOF',rdof,'name','RB modes');
 
%% #GeomSmoothEdge : place edge mid-node based on interpolated normals -------
elseif comstr(Cam,'smoothedge')
    
 model=varargin{carg};carg=carg+1;
 sel=varargin{carg};carg=carg+1;
 
 mo1=model;mo1.Elt=feutil(['selelt' sel],mo1);
 MAP=fe_sens('normalvel -map',mo1);
 elt=feutil('getedgelineall',mo1);
 cEGI=find(isfinite(elt(:,1)));
 NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
 n2=sparse(MAP.ID,1,1:size(MAP.ID,1));
 
 for jElt=1:length(cEGI)
  i2=full(NNode(elt(cEGI(jElt),1:3)));
  i3=n2(elt(cEGI(jElt),1:2));
  n1=midPointNormal( ...
      model.Node(i2(1),5:7),model.Node(i2(2),5:7), ... %p_a,p_b
      MAP.normal(i3(1),:),MAP.normal(i3(2),:)); %n_a,n_b
  %[model.Node(i2(1:2),5:7);model.Node(i2(3),5:7);n1]
  model.Node(i2(3),5:7)=n1;
  %feplot(model);fecom('shownodemark',model.Node(i2,1))
 end 
 out=model;
 
else; error('Geom%s unknown',CAM);
end
%% #Get ----------------------------------------------------------------------
elseif comstr(Cam,'get');  [CAM,Cam]=comstr(CAM,4);

%% #GetDOF - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% feutil('GetDof',model) 
if comstr(Cam,'dof'); [CAM,Cam]=comstr(CAM,4);

model=varargin{carg};carg=carg+1;
ModelStack={};RunOpt=struct('ImplicitNode',1);
if isnumeric(model)
 if isfinite(model(1))&&size(model,2)==7 % Node elt
  carg=carg-1;
  [carg,FEnode,FEelt,FEel0,ModelStack]=get_nodeelt(varargin,carg,ModelStack,RunOpt);
   model=struct('Node',FEnode,'Elt',FEelt);carg=carg+1;
 else % backward compatible elt only
   model=struct('Node',[],'Elt',model,'Stack',{{}});
 end
elseif isfield(model,'Elt')
 carg=carg-1;
 [carg,FEnode,FEelt,FEel0,ModelStack]=get_nodeelt(varargin,carg,ModelStack,RunOpt);
end

if any(size(model.Elt)==1)&&isfinite(model.Elt(1))
 % getdof: from node and dof list
 r1=model.Elt; r2=varargin{carg}; carg=carg+1;
 out=reshape(r2*ones(1,length(r1)),[],1)+...
  reshape(ones(length(r2),1)*r1',[],1); % sorted by the first arg.
 
else % getdof from model

if isa(model,'v_handle');model=model.GetData;end
if ~isempty(Cam) % element selection
 [out,model.Elt]=feutil(['findelt' CAM],FEnode,FEelt,FEel0,varargin{carg:end});
end
if ~isfield(model,'Stack');model.Stack=ModelStack;end
[eltid,model.Elt]=feutil('eltidfix;',model);
elt=model.Elt;r1=fe_case(model,'stack_get','rigid');
if ~isempty(r1)
 r3=cell(size(r1,1),1);
 for j1=1:size(r1,1) % robust to .Sel
  if isstruct(r1{j1,3}); r2=r1{j1,3}.Elt; else; r2=r1{j1,3}; end
  if ~isempty(r2) % authorize exitence of empty rigid entries
   if ~any(~isfinite(r2(:,1))); r3{j1}=feutil('addelt','rigid',r2);
   else; r3{j1}=r2; % robust if no header
   end
  end
 end
 r3(cellfun(@(x)size(x,1),r3)==0,:)=[]; 
 i1=max(cellfun(@(x)size(x,2),r3)); for j1=1:size(r3,1);r3{j1}(1,end+1:i1)=0;end
 elt=feutil('addelt',elt,vertcat(r3{:}));
end
[EGroup,nGroup]=getegroup(elt);

% determination of DOFs in the unrestrained model
% i3 generic node set, i4 generic dof set, out2 word count

ndof=[];edof=[]; out2=0; nd=[];

for jGroup=1:nGroup
  [ElemF,i1]=feutil('getelemf',elt(EGroup(jGroup),:),jGroup);EGID=i1(1);
  cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;
  i4=[];try;if ~strcmp(ElemF,'SE'); i4=feval(ElemF,'dofcall');end;end
  % see sdtweb p_solid('builddof')
  if EGID<0  % display group
  elseif ischar(i4); 
  % there is a dof call strategy wich is assumed to return 
  % i2 (nodal DOF*100), i3 element DOF*1000
  % beware this must increment out2 (word count) too
   eval(i4); ndof=[ndof;i2(:)];edof=[edof;abs(i3(:))]; %#ok<AGROW>
  elseif EGID>=0 % not a display group
   [i3,SEopt]=fe_super('node',ElemF,model);
   i4=fe_super('dof',ElemF,model);
   i4=i4(:)*100;i3=i3(:);
   if size(SEopt,2)<3||SEopt(1,3)==0; i6=length(i4)^2;else;i6=SEopt(1,3);end
   if     SEopt(1,1)==2; out2 = out2 + length(cEGI)*i6;  % generic
   elseif SEopt(1,1)==1; out2 = out2 + i6;  % unique
   end
   i1=find(i4<0);
   if ~isempty(i1) % element DOFs
      if SEopt(1,1)==2 && ~isempty(cEGI) % generic
        i2=eltid(cEGI);i2=i2(:)'; 
	    i2=i2(ones(size(i1)),:)*1e3-i4(i1,ones(size(cEGI)))*10-1000;
      else;i2=-i4(i1)*10+(999000+EGID*1000);
      end;% unique
      edof=[edof;i2(:)]; %#ok<AGROW>
   end
   i1=find(i4>0);
   if (~isempty(i1)&&strcmp(ElemF,'rigid'))||strcmp(ElemF,'celas')
   % master and slave DOFs in rigid & celas
     i2=[];
     [i3,i8,i7]=unique(elt(cEGI,3:4),'rows');
     for j1=1:size(i3,1)
      i8=cEGI(i7==j1);
      if i3(j1,1)>0
       i4=[1:6];i5=[abs(sprintf('%i',i3(j1,1)))-48];
      elseif i3(j1,1)<0 && i3(j1,2)==0
       i4=[abs(sprintf('%i',-i3(j1,1)))-48];i5=i4;
      elseif i3(j1,1)<0 
       i4=[abs(sprintf('%i',abs(i3(j1,1))))-48];
       i5=[abs(sprintf('%i',abs(i3(j1,2))))-48];
      elseif i3(j1)==0; 
       sdtw('_nb','Missing DOFs for %s (elements %i:%i)',ElemF,cEGI([1 end]));
       i4=[];i5=[];
      end
      if ~isempty(i4);
       i4=elt(i8,ones(1,length(i4)))*100+i4(ones(length(i8),1),:);
       i5=elt(i8,2*ones(1,length(i5)))*100+i5(ones(length(i8),1),:);
       i2=[i2;i4(:);i5(:)]; %#ok<AGROW>
      end
     end
     i2=unique(i2(i2>100)); ndof=[ndof;i2]; %#ok<AGROW>
   elseif ~isempty(i1) % nodal DOFs
      if SEopt(1,1)==2 && ~isempty(cEGI) % generic
	    i3=i3(:);nind=sparse(i3,1,1:length(i3));%(i3,1)=[1:length(i3)]';
        i4=[nind(fix(i4(i1)/100)) round(rem(i4(i1),100))]';
        if max(i3)>size(elt,2);elt(1,max(i3))=0;end
        i2 = elt(cEGI,i3)*100;
        if strcmp(ElemF,'mass1') 
         i2=i2(:,i4(1,:)).*(elt(cEGI,2:7)~=0)+i4(2*ones(size(i2,1),1),:);
        elseif strcmp(ElemF,'mass2') 
         i5=elt(cEGI,:); if size(i5,2)<13;i5(1,13)=0;end
         i6=any(i5(:,[3:8 10:13]),2); % with inertia or offset
         i2=[i5(i6,i4(1,:))*100+i4(2*ones(nnz(i6),1),:)];
         i4(:,4:6)=[];i5(i6,:)=[];
         i7=[i5(:,i4(1,:))*100+i4(2*ones(size(i5,1),1),:)];
         i2=[i2(:);i7(:)];
        else
           i2=i2(:,i4(1,:))+i4(2*ones(size(i2,1),1),:);
        end
      else;i2=round(i4(i1));
      end
      i2=unique(i2(:)); ndof=[ndof;i2]; %#ok<AGROW>
   end
  end % of EGID>=0
end
Case=fe_case(model,'getcase');
r1=stack_get(Case,'rbe3');
for j1=1:size(r1,1) % DOF slave RBE3 DOFs
  r3=r1{j1,3}; if isfield(r3,'data');r3=r3.data;end
  if isstruct(r3);continue;end % implicit rbe3
  r2=num2cell(r3(:,2:3));
  for j2=1:size(r2,1);
      r2{j2,1}=r2{j2,1}*100+abs(sprintf('%i',abs(r2{j2,2})))'-48;
  end
  ndof=[ndof;vertcat(r2{:,1})]; 
end
ndof=unique(round(ndof))/100; ndof=ndof(ndof>1);
out = [ndof; -unique(edof)/1000];
out1=[];

if carg<=nargin&&~isempty(varargin{carg})
    if carg+1<=nargin&&~isempty(varargin{carg+1})
      r1=varargin{carg};i1=fe_c(r1,out,'ind');
      out1=varargin{carg+1};out1=out1(i1,:);out=r1(i1);
    else;if ~isempty(out); out=fe_c(out,varargin{carg},'dof'); end;out1=[];
    end    
end

if nargout>2; out3=model.Elt;out4=eltid;end
end % from model, or from node/dof list

%% #GetEdge  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% elt=feutil('getedge[Patch,Line]')
elseif comstr(Cam,'edge'); [CAM,Cam]=comstr(CAM,5);

model=[];out1=[];
[carg,node,elt,el0,ModelStack]=get_nodeelt(varargin,carg,ModelStack);
if comstr(Cam,'patch') 
 %% #getEdgePatch - - - - - - - - - - - - - - - - - - - -
[CAM,Cam]=comstr(CAM,6);
% SelEdgePatch replaces the selected group(s0 by a quad4 group containing
% edge patches
 i4=[];
 [EGroup,nGroup]=getegroup(elt); i6=[];
 % note : this will not work if two face sets have the same number of nodes

 %[eltid,elt]=feutil('eltidfix;',elt); % done later now
 RunOpt.ConvFcn=''; RunOpt.iNum=6; % Face Id numbering column in i5 to output
 FaceStack={}; RunOpt.FaceId=[];out2=[];r1=[];RunOpt.All=0; RunOpt.SNeg=0;
 if carg>nargin
 elseif isfield(varargin{carg},'type');r1=varargin{carg};carg=carg+1;
 else; st=varargin{carg};
  if ischar(st);
   [st,un1] = comstr('pos',[-25 3],st,lower(st)); % selfacepos is default
   [st,un1,RunOpt.SNeg] = comstr('neg',[-25 3],st,lower(st)); % selfaceneg: flip elts
   [st,un1,RunOpt.IdCols] = comstr('idcols',[-25 2],st,lower(st));
   if ~isempty(RunOpt.IdCols)
   elseif ~isempty(strfind(st,'pro'));RunOpt.IdCols=4;
   elseif strncmpi(st,'b',1); RunOpt.IdCols=2:4;st(1)='';%selfaceb keeps separates by mpg
   end
   [st,un1,RunOpt.trim] = comstr('-trim',[-25 32],st,lower(st));
   if RunOpt.trim==1; RunOpt.trim=60; end
   [st,un1,i1] = comstr('lin',[-25 3],st,lower(st));if i1;RunOpt.FaceCmd='facelin';end
  end
  if ~ischar(st)
  elseif strncmp(st,'@',1);RunOpt.ConvFcn=st(2:end);
  elseif strcmpi(st,'-all');RunOpt.All=1;
  else;r1=stack_get(model,'set',st,'getdata');
  end
 end
 if ~isfield(RunOpt,'trim')||isempty(RunOpt.trim); RunOpt.trim=0; end
 
 if sp_util('issdt')&&RunOpt.trim 
  %% surface trimming
  % assume already based on a selface result, post-treat
  % split surfaces as groups, but keep face identifiers: orig EltId FaceNum
  mo1=struct('Node',node,'Elt',elt);
  elid0=feutil('EltId;',elt); [elid1,mo1.Elt]=feutil('EltIdFix;',mo1.Elt);
  mo1=feutilb(sprintf('SurfaceAsQuadGroup-isFaceSel-storeInd %g',RunOpt.trim),mo1,'');
  [mo1,ind]=stack_rm(mo1,'info','NewEltInd','get'); ind=ind{3};
  if sum(~isfinite(mo1.Elt(:,1)))>1;
   el1=feutil('SelElt SelEdge',mo1); % find nodes on the edge
   if ~isempty(el1) % no edge, nothing to trim
   n1=unique(reshape(el1(isfinite(el1(:,1)),1:2),[],1));
   mpid=feutil('mpid',mo1.Elt);
   i1=feutil('FindElt WithoutNode{NodeId}',mo1,n1);
   r1=unique(mpid(i1,3)); % groups not touching the surface edges
   r2=setdiff(unique(mpid(:,3)),[0;r1(:)]);
   i1=feutil(sprintf('FindElt group %s & WithNode{NodeId}',sprintf('%i ',r2)),mo1,n1);
   ind(i1)=[]; elt=elt(unique([find(~isfinite(elt(:,1)));ind(ind>1)]),:);
   elt=feutil('OptimEmptyGroups',elt);
   end
  end
  
 else
  %% transform into face selection
  [eltid,elt]=feutil('eltidfix;',elt);
 RunOpt.SubcEGI=[];
 if isfield(r1,'type')&&strcmpi(r1.type,'FaceId');
   RunOpt.FaceId=r1.data;%[i1,i2]=ismember(r1.data(:,1),eltid);
   if isfield(r1,'ConvFcn'); RunOpt.ConvFcn=r1.ConvFcn;end
   if isfield(r1,'FaceCmd'); RunOpt.FaceCmd=r1.FaceCmd;end
   %RunOpt.FaceId(:,1)=i2;RunOpt.FaceId(~i1,:)=[];
   out2=zeros(size(RunOpt.FaceId,1),1);RunOpt.iNum=6;
   % set with ConvFcn, recover SDT numbering
   if ~RunOpt.All % use SubcEGI to avoid generating FaceStackAll for nothing
    RunOpt.SubcEGI=find(ismember(eltid,r1.data(:,1)));
   end
   
 elseif ~isempty(RunOpt.ConvFcn);
  RunOpt.iNum=1; % command with @, need to convert output to external convention
 end
 for jGroup = 1:nGroup
  
   cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1; 
   if ~isempty(RunOpt.SubcEGI); cEGI=intersect(cEGI,RunOpt.SubcEGI); end
   if isempty(cEGI); continue; end
   
   [ElemF,i1,ElemP]= feutil('getelemf',elt(EGroup(jGroup),:),jGroup);
   [r1,i1]=feutil_get_face(ElemF,RunOpt); 
   if RunOpt.iNum==1; i1=i1(i1); end % conversion direction @:SDT2ext, r1L.ConvFcn;ext2SDT
   i1=i1(i1);
   i2=fe_super('prop',ElemF);
   if isempty(r1)
   else
    % decompose in face subsets for variable node numbers
    r2={};r3=(1:size(r1,1));
    while ~isempty(r1)&&size(r1,2)>1;
     i3=find(r1(:,end-1)~=r1(:,end));
     if ~isempty(i3); r2(end+1,1:3)={r1(i3,:) i1(i3) r3(i3(:))};%#ok<AGROW>
         r1(i3,:)=[];i1(i3)=[];r3(i3)=[];
     end
     r1(:,end)=[];
    end

    if ~isempty(r1);r2(end+1,1:2)={r1,i1,[1:length(i1)]'};end %#ok<AGROW>
    
    for j1=1:size(r2,1);
     r1=r2{j1,1}; % loop on number of node per face
     if length(FaceStack)<size(r1,2); FaceStack{2,size(r1,2)}=[];end %#ok<AGROW>
     % i5 [facenumber jGroup matid proid eltid origfacenum]
     i5=repmat(r2{j1,3},1,length(cEGI));i5=i5(:);i5(:,2)=jGroup; 
     if i2(1); i3=elt(cEGI,i2(1)*ones(1,size(r1,1)))';i5(:,3)=i3(:); end
     if i2(2); i3=elt(cEGI,i2(2)*ones(1,size(r1,1)))';i5(:,4)=i3(:); end
     if i2(3)&&size(elt,2)>=i2(3)
            i3=elt(cEGI,i2(3)*ones(1,size(r1,1)))';i5(:,5)=i3(:);
     end
     if ~isempty(cEGI); i5(:,6)=repmat(r2{j1,2}(:)',1,length(cEGI)); 
     else; i5=double.empty(0,6);
     end
     i3=reshape(elt(cEGI,r1')',size(r1,2),size(r1,1)*length(cEGI));
     FaceStack{1,size(r1,2)}=[FaceStack{1,size(r1,2)};i3']; %#ok<AGROW>
     FaceStack{2,size(r1,2)}=[FaceStack{2,size(r1,2)};i5];  %#ok<AGROW>
    end
   end
 end
 elt=[];
 for j1=1:size(FaceStack,2)
  if ~isempty(FaceStack{1,j1})
   i4=FaceStack{1,j1};
   if RunOpt.All; i1=1:size(i4,1);i3= true(size(i1));i2=[];
   elseif isempty(RunOpt.FaceId)&&~isempty(RunOpt.IdCols)
    % eliminate duplicated faces with given Ids  SelFacePro
    [i2,i1]=sortrows(sort([i4 FaceStack{2,j1}(:,RunOpt.IdCols)],2)); 
    i3=find(~any(diff(i2),2));i1(i3)=0;i1(i3+1)=0;i1=i1(i1~=0);i2=[];
   elseif isempty(RunOpt.FaceId) % Just duplicate cols      
    [i2,i1]=sortrows(sort(i4,2)); % eliminate duplicated faces
    i3=find(~any(diff(i2),2));i1(i3)=0;i1(i3+1)=0;i1=i1(i1~=0);i2=[];
   else;% RunOpt.FaceId, comes from a FaceId set, identify in elements
    % handle here SNeg from FaceId sets   
    [i1,i2]=ismember(FaceStack{2,j1}(:,[5 1]),abs(RunOpt.FaceId(:,1:2)),'rows');
    if size(RunOpt.FaceId,2)>2 % Propagate subgroup info to MatId
       FaceStack{2,j1}(i1,3)=RunOpt.FaceId(i2(i1),3);
    end
    i1=find(i1); i2=i2(i1);i3=size(elt,1)+1+[1:length(i1)];
    if any(RunOpt.FaceId(:,2)<0) 
     % we do not do SNeg,  but fliplr on the fly those with negative FaceId
     in1=RunOpt.FaceId(i2,2)<0; in1=i1(in1);
     i4(in1,:)=fliplr(i4(in1,:));
     FaceStack{2,j1}(in1,RunOpt.iNum)=-FaceStack{2,j1}(in1,RunOpt.iNum);
    end
   end
   i4=i4(i1,:); r4=FaceStack{2,j1}(i1,[3:5 RunOpt.iNum]);
   if RunOpt.SNeg % flip topo and refer to -1*FaceId
    i4=fliplr(i4); r4(:,end)=-r4(:,end);
   end
   switch size(i4,2)
   case 3
     elt(end+[1:length(i1)+1],1:7) = [Inf abs('tria3') 0; i4 r4];  %#ok<AGROW>
           %i4(i1,:)  FaceStack{2,j1}(i1,[3:5 RunOpt.iNum])]; %#ok<AGROW>
   case 4
     elt(end+[1:length(i1)+1],1:8) = [Inf abs('quad4') 0 0; i4 r4]; %#ok<AGROW>
           %i4(i1,:)  FaceStack{2,j1}(i1,[3:5 RunOpt.iNum])]; %#ok<AGROW>
   case 6
     elt(end+[1:length(i1)+1],1:10) = [Inf abs('tria6') 0 0 0 0; i4 r4];  %#ok<AGROW>
           %i4(i1,:)  FaceStack{2,j1}(i1,[3:5 RunOpt.iNum])]; %#ok<AGROW>
   case 8
     %i4=i4(i1,:); % Here need to revert due to post-pro
     if RunOpt.SNeg; i4=fliplr(i4); end
     ind=find(i4(:,6)==i4(:,7)&i4(:,6)==i4(:,8)&i4(:,6)&i4(:,7)&i4(:,8));
     in1=setdiff(1:size(i4,1),ind);
     if ~isempty(in1)
      i3(in1)=size(elt,1)+1+[1:length(in1)];
      i5=i4(in1,1:8); if RunOpt.SNeg; i5=fliplr(i5); end
      elt(end+[1:length(in1)+1],1:12)=[Inf abs('quadb') zeros(1,6) ; i5  r4(in1,:)]; %#ok<AGROW>
           %i4(in1,1:8)  FaceStack{2,j1}(i1(in1),[3:5 RunOpt.iNum])]; %#ok<AGROW>
     end
     if ~isempty(ind) % degenerate tria6
      i3(ind)=size(elt,1)+1+[1:length(ind)];
      i5=i4(ind,1:6); if RunOpt.SNeg; i5=fliplr(i5); end
      elt(end+[1:length(ind)+1],1:10)=[Inf abs('tria6') zeros(1,4) ; i5  r4(ind,:)];  %#ok<AGROW>
           %i4(ind,1:6)  FaceStack{2,j1}(i1(ind),[3:5 RunOpt.iNum])]; %#ok<AGROW>
     end
    case 9
     elt(end+[1:length(i1)+1],1:13) = [Inf abs('quad9') zeros(1,7); i4 r4];  %#ok<AGROW>
           %i4(i1,:)  FaceStack{2,j1}(i1,[3:5 RunOpt.iNum])]; %#ok<AGROW>
    case 2 % faces for 2D/shell elements
     elt(end+[1:length(i1)+1],1:6) = [Inf abs('beam1') ; i4 r4];  %#ok<AGROW>
      %i4(i1,:)  FaceStack{2,j1}(i1,[3:5 RunOpt.iNum])]; %#ok<AGROW>
   otherwise; error('Not supported get face with %i nodes',size(i4,2));
   end
   if ~isempty(i2); out2(i2)=i3;end
  end
 end
 end 

% 'getEdgeLine[ProId,MatId,All,In]'  - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'line') 
 %% #GetEdgeLine - - - - - - - - - - - - - - - - - - -

 isSilent=CAM(end)==';';if isSilent;CAM(end)='';end
 [CAM,Cam] = comstr(CAM,5);
 
 % GetEdgeLine replaces the selected group by a beam1 group containing 
 % segments at the edge of the of the current group

 if comstr(Cam,'g'); [CAM,Cam]=comstr(CAM,'group','%c');    opt=1;
 elseif comstr(Cam,'p'); [CAM,Cam]=comstr(CAM,'proid','%c');opt=2;
 elseif comstr(Cam,'m'); [CAM,Cam]=comstr(CAM,'matid','%c');opt=3;
 elseif comstr(Cam,'a'); [CAM,Cam]=comstr(CAM,'all','%c');opt=5;
  if any(strfind(Cam,'-nouni')); opt=6; end
 elseif comstr(Cam,'i'); opt=4; ind=varargin{carg};carg=carg+1;
     [a,elt]=feutil('eltidfix;',elt);
 else;opt=0;
 end

 i4=[];[EGroup,nGroup]=getegroup(elt); i6=[];LD2=zeros(nGroup,1);
 for jGroup = 1:nGroup
   ElemF= feutil('getelemf',elt(EGroup(jGroup),:),jGroup);
   cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;LD2(jGroup,1)=0;
   [NodeLine,SEopt]=fe_super('edge',ElemF);NodeLine=NodeLine';
     
   if ~isempty(cEGI); i3=elt(cEGI,NodeLine(:))';
      % case of a single superelement
   else; [i3]=fe_super('node',ElemF);i3=i3(NodeLine);     
   end
   i5=ones(size(NodeLine,2)*size(i3,2),1)*jGroup;
   i5(:,5)=reshape([1:size(NodeLine,2)]'*ones(1,size(i3,2)),[],1);
   i1=fe_super('prop',ElemF);
   if i1(1)&&size(elt,2)>=i1(1)
    i2=elt(cEGI,i1(1)*ones(1,size(NodeLine,2)))';i5(:,2)=i2(:); 
   end
   if i1(2)&&size(elt,2)>=i1(2) 
    i2=elt(cEGI,i1(2)*ones(1,size(NodeLine,2)))';i5(:,3)=i2(:);
   end
   if i1(3)&&size(elt,2)>=i1(3)
    i2=elt(cEGI,i1(3)*ones(1,size(NodeLine,2)))';i5(:,4)=i2(:);
   end

   i3=reshape(i3,size(NodeLine,1),size(NodeLine,2)*length(cEGI))';
   i7=~all(i3,2);
   if any(i7);fprintf('Removing %i edges with zero nodes\n',nnz(i7));
      i3(i7,:)=[]; i5(i7,:)=[]; 
   end
   %i4(end+1:end+size(i3,1),1:size(i3,2))=i3;
   %i6(end+1:end+size(i5,1),1:size(i5,2))=i5;
   i4(1:size(i3,2),end+1:end+size(i3,1))=i3';
   i6(1:size(i5,2),end+1:end+size(i5,1))=i5';
 end % loop on groups
 i4=i4'; i6=i6';

 if ismember(opt,[5 6]) % all Edges
  [i3,i1,i2]=unique(sort(i4,2),'rows');
  if opt==5; i4=i4(i1,:);i6=i6(i1,:); % single instance
  else; i6(:,end+1)=i2; % keep duplicates with markers "edgeid"
  end
  
 elseif opt==4 % EdgesInNodes

   nind=zeros(max(node(:,1)),1);nind(ind)=1:length(ind);
   i3=find(all(nind(i4),2)); 
   i4=i4(i3,:); i6=i6(i3,:); 
   [i4,i1]=sortrows(sort(i4,2));  i6=i6(i1,:);
   i2=find(~any(diff(i4),2)); i6(i2,1)=i6(i2,4); % EltNumber in col 1
   nind=[1:size(i4,1) 0]; nind(i2)=0; j1=1;
   while any(nind)
    i3=i2+j1; if ~isempty(i3); i5=find(nind(i3)); else;i5=[];end
    if ~isempty(i5);   j1=j1+1;i6(i2(i5),j1)=i6(i3(i5),4); nind(i3(i5))=0;end
   end
   i6=[i2(:) i6(i2,1:j1)]; i4=i4(i2,:);

 elseif opt % some gpm selection

  ind=find(sparse(i6(:,opt)+1,1,1))-1; % selection values
  i3=[]; i5=[];
  for j1=1:length(ind) % edge of each group/mat/prop
    in2=find(i6(:,opt)==ind(j1));
    [i1,i2]=single_line(i4(in2,:),i6(in2,:));
    i3=[i3;i1];i5=[i5;i2]; %#ok<AGROW>
  end
  % now eliminate duplications & create new number of interface boundaries
  [i3,i1]=sortrows(i3);   i2=find(~any(diff(i3),2));
  i5(i1(i2),opt)=sort([i5(i1(i2),opt) i5(i1(i2+1),opt)],2)* ...
          [1;10^ceil(log10(length(ind)))];
  in2=1:size(i3,1);in2(i2+1)=0;i2=find(in2);
  i3=i3(i2,:);i5=i5(i1(i2),:);
  i5(:,2)=i5(:,opt); % new numbers as materials
  i4=i3;i6=i5;
 else % standard selection of edges
  [i4,i6]=single_line(i4,i6);
 end
 
 if isempty(i4)
   if ~isSilent&&~Silent; sdtw('GetEdge: selected group has no edges'); end
   elt=[];
 else 
    if size(i4,2)<3; i4(:,3)=0;end
    i1=find(i4(:,3)==0&i4(:,1));%elt=[];
    elt=zeros(size(i4,1),max(6,3+size(i6,2)-1)); ie=0;
    if ~isempty(i1)
     %elt = [Inf abs('beam1')];
     %elt(2:length(i1)+1,1:1+size(i6,2))=[i4(i1,1:2) i6(i1,2:end)];
     elt(1,1:6)=[Inf abs('beam1')];
     elt(2:1+length(i1),1:2)=i4(i1,1:2); elt(2:1+length(i1),3:(1+size(i6,2)))=i6(i1,2:end);
     i4(i1,:)=0; ie=1+length(i1);
    end
    i1=all(i4,2); %find(all(i4,2));
    if any(i1) %~isempty(i1)
      %elt(end+1,1:6) = [Inf abs('beam3')];
      %elt(end+1:end+length(i1),1:2+size(i6,2))=[i4(i1,[1:3]) i6(i1,2:end)];
      elt(ie+1,1:6)=[Inf abs('beam3')];
      elt(ie+2:ie+1+nnz(i1),1:3)=i4(i1,1:3); 
      elt(ie+2:ie+1+nnz(i1),4:(2+size(i6,2)))=i6(i1,2:end);
      i4(i1,:)=0; ie=ie+1+nnz(i1);
    end
    elt(ie+1:end,:)=[];

    if any(any(i4))
     sdtw('Edges with 2 or 3 nodes are the only supported');
    end

 end

 if ~isempty(elt);
     elt(isfinite(elt(:,1))&elt(:,1)==elt(:,2),:)=[];%eliminate nodes
 end
 
else;error('Not a known SelEdge command');
end % of SelEdgeLine

[CAM,Cam]=comstr(CAM,1); 
if ~isempty(CAM)&&CAM(1)=='{'&&CAM(end)=='}';[CAM,Cam]=comstr(CAM(2:end-1),1);end
out1=CAM;out=elt; 

%% #GetElemF - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%	[ElemF,opt]= feutil('getelemf',elt(EGroup(jGroup),:),jGroup)
elseif comstr(Cam,'elemf'); [CAM,Cam]=comstr(CAM,6);

 if comstr(Cam,'list')
  %% #GetElemFList -2
  elt=varargin{carg}; carg=carg+1;
  elt=elt(~isfinite(elt(:,1)),:); 
  out=cell(size(elt,1),1); out1=out; out2=out;
  for j1=1:length(out); [out{j1},out1{j1},out2{j1}]=feutil('getelemf',elt(j1,:)); end
  
 else
  %% #GetElemFBase -2
   st=varargin{carg};carg=carg+1; 
   i1 = [min([find(~st) find(st==32)]) length(st)+1];
   out = char(st(2:i1(1)-1));
   if nargout>1 
    out1 = st(i1(1)+1:length(st)); 
    if carg<=nargin; r1=varargin{carg};carg=carg+1; else;r1=0; end
    if isempty(out1); out1=r1; elseif out1(1)==0; out1=r1(1); end
   end
   if nargout>2    ;out2=fe_super('parent',out);  end
 end

%% #GetLine  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% this should be maintained for edge selection
% feutil('GetLine',node,elt)
elseif comstr(Cam,'line'); [CAM,Cam]=comstr(CAM,5);

[carg,node,elt,ModelStack]=get_nodeelt(varargin,carg,ModelStack);

if carg>nargin; nInf=0;else;nInf=varargin{carg};carg=carg+1;end
if ~isempty(elt)
  EGroup = [find(elt(:,1)==Inf)' size(elt,1)+1];nGroup = length(EGroup)-1;
else	nGroup=0;EGroup=[];elt([1 83])=[1 0]; 
end

if nGroup==0 % trace line plot - - - - - - - - - - - - - - - - - - - -

   LD2=zeros(size(elt,1),1);
   for j1 = 1:size(elt,1);
      i1=elt(j1,[83:82+elt(j1,1)]);LD2(j1,1:length(i1)+1) = [length(i1)+1 i1];
   end

elseif nGroup>0 % element plot - - - - - - - - - - - - - - - - - - - -

LD2(nGroup,1)=0; out=struct('Group',[]);out.Group={};out.Line={};

for jGroup = 1:nGroup
   ElemF= feutil('getelemf',elt(EGroup(jGroup),:),jGroup);
   cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;
   [NodeLine,SEopt]=fe_super('line',ElemF);NodeLine=NodeLine(:)';

   if (~isempty(NodeLine)&&~isempty(cEGI)) || SEopt(1,1)==1

     if length(NodeLine)>1
      NodeLine=[NodeLine(1:length(NodeLine)-1);NodeLine(2:length(NodeLine))];
     else;NodeLine=[NodeLine;NodeLine]; %#ok<AGROW>
     end
     NodeLine=NodeLine(:,NodeLine(1,:)&NodeLine(2,:));
     NodeLine=NodeLine(:);
     
     if ~isempty(cEGI); i3=elt(cEGI,NodeLine)'; % case of a unique superelement
     else; [i3]=fe_super('node',ElemF);i3=i3(NodeLine,1); end

     if ~isempty(i3)

        i3=reshape(i3,2,numel(i3)/2);
        i3=sort(i3(:,all(i3)));% ignoring nodes set to zero
        [i1,i2]=find(sparse(i3(1,:),i3(2,:),i3(1,:)));
        i1=[i1 i2 zeros(size(i2))]'; i1=i1(:);
        LD2(jGroup,1:length(i1)) = [length(i1) i1(1:length(i1)-1)'];
        out.Group{end+1}=sprintf('Group %i',jGroup);
        out.Line{end+1}=i1(1:length(i1)-1)';
     elseif isempty(i3);        LD2(jGroup,1) = [0];
     end % all(i3) the nodes are in FEnode
   else;disp(sprintf('No ''line''  plot for group %i (%s)',jGroup,ElemF));
   end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end
end % jGroup and type of nGroup

if ~isempty(strfind(Cam,'struct'))
 return;
elseif nargout>2
 i1 =LD2(:,2:end);i1=i1(:);out2=zeros(max(i1)+1,1);out2(i1+1)=1;
 out2=find(out2)-1;out2=out2(2:end);
end
if ~isempty(node) % if node is defined give in indices rather than numbers
 nInf=size(node,1);
 if isfinite(node(nInf,5))
   node=[node;max(node(:,1))+1 0 0 0 Inf Inf Inf];nInf=nInf+1;
 else;node(nInf,1)=max(node(1:nInf-1,1))+1;
 end
 NNode=sparse(node(:,1),1,1:size(node,1));
 for j1=1:size(LD2,1); 
  if LD2(j1,1)>0
   i1 =  LD2(j1,2:LD2(j1,1)); i2=find(i1==0); i1(i2)=i1(i2)+node(nInf,1);
   LD2(j1,1:length(i1)+1)=[length(i1)+1 reshape(full(NNode(i1)),1,length(i1))];
  end
 end
 if nargout>2; out2=[full(NNode(out2))';nInf];end
end
out = LD2; out1=nInf;

elseif comstr(Cam,'patch'); [CAM,Cam]=comstr(CAM,6);
%% #GetPatch feutil('GetPatch',model) - - - - - - - - - - - - - - - - - - - - 

model=[];
[carg,node,elt,el0,ModelStack]=get_nodeelt(varargin,carg,ModelStack);
if isempty(model); model=struct('Stack',{ModelStack});end

nInf=size(node,1);
if carg>nargin; 
elseif isstruct(varargin{carg});RunOpt=varargin{carg};
else;nInf=varargin{carg};carg=carg+1;
end
[nElt,i1]=size(elt);nInf=size(node,1);
if nElt>0; [EGroup,nGroup]=getegroup(elt);% elements have been defined
else;	nGroup=0;EGroup=[];elt([1 83])=[1 0];  
end
[CAM,Cam,RunOpt.Sci]=comstr('sci',[-25 3],CAM,Cam);
if ~isfield(RunOpt,'LinFace')
 [CAM,Cam,RunOpt.LinFace]=comstr('-linface',[-25 3],CAM,Cam);
end
RunOpt.LinSt='';if RunOpt.LinFace;RunOpt.LinSt='-linface';end
[CAM,Cam,RunOpt.KeepAll]=comstr('-keepall',[-25 3],CAM,Cam);
RunOpt.Post={};
if comstr(Cam,'new') % 'GetPatchNew returns a feplot selection

 out=struct('opt',[0 0 0 0 0],'selelt','', ...
     'vert0',[],'cna',[],'fs',[],'ifs',[],'f1',[],'if1',[],'f2',[],'if2',[],...
     'fsProp',[],'f2Prop',[],'f1Prop',[]);
 out.fsProp={'FaceColor',[1 1 1],'edgecolor',[0 0 0],'facelighting','phong'};
 out.f2Prop={'FaceColor','none','edgecolor',[0 0 0],'linewidth',2};
 out.f1Prop={'FaceColor','none','edgecolor',[0 0 0],'marker','o'};
 out.Node=[];out.MPID=[]; out.ScaleColorMode='one';
end
if nGroup==0 % trace line plot - - - - - -
 i2 =elt; elt=[];
 for j1=1:size(i2,1);
  i1=[i2(j1,82+[1:i2(j1,1)-1]);i2(j1,82+[2:i2(j1,1)])];
  i1=i1(:,~any(i1==nInf)&~any(i1==0))';
  elt=[elt;Inf abs('beam1'); i1 zeros(size(i1,1),4)]; %#ok<AGROW>
 end
 EGroup = [find(elt(:,1)==Inf)' size(elt,1)+1];nGroup = length(EGroup)-1;
end

LD2=[];i7=[];

for jGroup = 1:nGroup
   [ElemF,i1,ElemP]= getegroup(elt(EGroup(jGroup),:),jGroup);
   cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;
   EGID = i1(1); i7(end+1)=EGID; %#ok<AGROW>
   if RunOpt.Sci
     [NodeSurf,SEopt]=fe_super('sci_face',ElemF,model);
     if size(NodeSurf,2)==4;
         NodeSurf=[NodeSurf(:,[1 2 3]);NodeSurf(:,[3 4 1])];
     end
   elseif RunOpt.LinFace
     st={'beam3','beam1';'tetra10','tetra4';'hexa20','hexa8';'hexa27','hexa8';
         'tria6','tria3';'quadb','quad4';'penta15','penta6';'pyra13','pyra5'};
     i2=strcmp(st(:,1),ElemP);
     if ~any(i2);[NodeSurf,SEopt]=fe_super(['patch' RunOpt.LinSt],ElemF,model);
     else;[NodeSurf,SEopt]=fe_super(['patch' RunOpt.LinSt],st{i2,2},model);
     end
   else;[NodeSurf,SEopt]=fe_super(['patch' RunOpt.LinSt],ElemF,model);
   end
% new version of PATCH matrix corresponds to patch faces property
% xxx  GV tris to ensure out.ifj as a col. vector, seemed to produce uncoherent
% dimensions in general, at least recent MATLAB
if SEopt(1,1)==3  % locally built patch
  i2=NodeSurf;i7 = ones(size(i2,1),1); % needs revising
elseif ~isempty(NodeSurf)&&(SEopt(1,1)==1||~isempty(cEGI));  i2=[];
  if SEopt(1,1)==2 % generic superelement
   i2=[]; i7=[];
   for j1=1:size(NodeSurf,1) 
       i1 = NodeSurf(j1,:);
       i1 = reshape(elt(cEGI,i1),length(cEGI),length(i1));
       i2(end+[1:size(i1,1)],1:size(i1,2))=i1; i7 = [i7;cEGI(:)]; %#ok<AGROW>
   end % of loop on Patches linked to the element
  else % single superelement
    if ~isempty(setdiff(NodeSurf(:),node(:,1)));  % XXX FESUPER CHECK
     i2=reshape(node(NodeSurf,1),size(NodeSurf,1),size(NodeSurf,2));
    else; i2=NodeSurf;
    end
    i7 = ones(size(i2,1),1);
  end
  % for volumes elimination of duplicated patches within same group
  if size(i2,2)>2 && size(NodeSurf,1)>1
   i1 = sort(i2,2); [i1,i3]=sortrows(i1);
   i4=find(~any(diff(i1),2));
   i5=ones(size(i2,1),1);
   if ~RunOpt.KeepAll;
       i5(i3(i4))=0;i5(i3(i4+1))=0; % Eliminate duplicated
   end
   i5=find(i5);
   i2=i2(i5,:);i7=i7(i5);
  end
elseif isfield(NodeSurf,'fs')
else 
   fprintf('No ''patch'' plot for group %i (%s)\n',jGroup,ElemF);
   i2=[]; i7=[];
end

 if comstr(Cam,'new') % returns a feplot selection
   if     size(i2,2)==1  % Points
    out.f1(end+[1:size(i2,1)],1)=i2;
    out.if1(end+[1:length(i7)],1)=i7;
   elseif size(i2,2)==2 % Lines
    i3=find(any(i2==0,2)|i2(:,2)==i2(:,1));
    if ~isempty(i3) % Show degenerate as points
      out.f1(end+[1:length(i3)],1)=max(i2(i3,:),[],2);
      out.if1(end+[1:length(i3)],1)=i7(i3);
      i2(i3,:)=[]; i7(i3)=[];
    end
    if ~isempty(i2) % if some non degenerates
     out.f2(end+[1:size(i2,1)],1:2)=i2; out.if2(end+[1:length(i7)],1)=i7;
    end
   elseif  EGID==-4 % level set edge
     RunOpt.Post(end+1,1:3)= ...
         {'lsutil','edgeSelf2',struct('elt',elt(i7,:),'model',model)};
   elseif size(i2,2)>2 % standard faces
     if ~isempty(out.fs) && size(out.fs,2)<size(i2,2) 
       out.fs(:,end+1:size(i2,2))= ...
          out.fs(:,end)*ones(1,size(i2,2)-size(out.fs,2));
     end
     % filling of zeros xxx needs a systematic retest
     i1=find(~i2);
     while ~isempty(i1)
       if any(i1<size(i2,1)); 
        error('First node in face cannot be 0, jGroup=%i, ElemF=%s', ...
          jGroup,ElemF);
       end
       i2(i1)=i2(i1-size(i2,1)); i1=find(~i2); 
     end
     if size(i2,2)<size(out.fs,2)
       i2(:,end+1:size(out.fs,2))=i2(:,end)*ones(1,size(out.fs,2)-size(i2,2));
     end
     out.fs(end+[1:size(i2,1)],1:size(i2,2))=i2;
     out.ifs(end+[1:length(i7)],1)=i7;
   end  


 else % need something for old/feplot compatibility
   i2(:,5)=i7;i2(:,6)=jGroup;
   LD2(end+[1:size(i2,1)],1:size(i2,2))=[i2]; %#ok<AGROW>
 end
end % jGroup  - - - - -

if comstr(Cam,'new') % returns a feplot selection
  %NNode=sparse(node(:,1)+1,1,1:size(node,1));

  if isfield(RunOpt,'AllNodes') % Keep all nodes
   if isfield(RunOpt,'NNode'); NNode=RunOpt.NNode;
   else;NNode=sparse(out.Node,1,1:size(out.Node,1));   end
   if isfield(RunOpt,'Node'); out.Node=RunOpt.Node; 
   else; out.Node=node(:,1); end
   if isfield(RunOpt,'vert0'); out.vert0=RunOpt.vert0; 
   else; out.vert0=node(:,5:7);end
   %i2=[1:length(out.Node)]';
  else  % get nodes used by the faces
   if isfield(RunOpt,'NNp'); NNp=RunOpt.NNp;
   else;NNp=sparse(node(:,1)+1,1,1:size(node,1));
   end
   i2=[out.fs(:);out.f1(:);out.f2(:)];
   out.Node=find(sparse(i2+1,1,i2))-1;
   if max(out.Node)>length(NNp); NNp(max(out.Node)+1)=0; end
   i2 =NNp(out.Node+1);  i1=i2<1;
   if any(i1)
    % Display 10 last nodes in .Elt but not in .Node
    i1=find(i1);if length(i1)>10; i1=i1(end-9:end); end
    fprintf('NodeId [%s] in .Elt but not in .Node\n',sprintf('%i ',out.Node(i1)));
    error('NodeId in .Elt but not in .Node');
   end
   out.vert0=node(i2,5:7);
   NNode=sparse(out.Node,1,1:size(out.Node,1));
  end
%   %   if max(out.Node)>length(NNode); NNode(max(out.Node)+1)=0; end
%   %   i2 =NNode(out.Node+1);  
%   i1=i2<1;
%   if any(i1)
%     fprintf('%i %i %i %i %i %i\n',out.Node(i1));
%     error('NodeId in .Elt but not in .Node');
%   end
%   out.vert0=node(i2,5:7);
%   % renumber faces
%   NNode=sparse(out.Node,1,1:size(out.Node,1));
%   %NNode(out.Node)=[1:length(out.Node)]';
  if ~isempty(out.fs) 
   out.opt(1,3)=1;i1=find(out.fs);out.fs(i1)=NNode(out.fs(i1));
   if nnz(out.fs)~=numel(out.fs);error('Some faces have zero nodes');end
  end
  if size(out.fs,1)==size(out.vert0,1)&&~isempty(out.fs)% duplicate face for color control
    out.fs=out.fs([1:end end],:); out.ifs=out.ifs([1:end end]);
  end

  if ~isempty(out.f1)
   out.opt(1,5)=1;i1=find(out.f1);  out.f1(i1)=NNode(out.f1(i1));
   if nnz(out.f1)~=numel(out.f1);
     error('Some point faces have zero nodes');
   end
  end
  if ~isempty(out.f2) 
   out.opt(1,4)=1;i1=find(out.f2);  out.f2(i1)=NNode(out.f2(i1));
   if nnz(out.f2)~=numel(out.f2);
     error('Some 2 node faces have zero nodes');
   end
   i1=find(~out.f2(:,2));if ~isempty(i1); out.f2(i1,2)=out.f2(i1,1);end
   i1=find(~out.f2(:,1));if ~isempty(i1); out.f2(i1,1)=out.f2(i1,2);end
  end
  for j1=1:size(RunOpt.Post,1)
   out=feval(RunOpt.Post{j1,1}, ...
       out,struct('CAM','GetPatch','jPost',j1),RunOpt.Post{j1,2:end});
  end
  try;
   if isfield(varargin{2},'IndWithHeaders')
    i2=varargin{2}.IndWithHeaders;
    if isfield(out,'ifs'); out.ifs=i2(out.ifs); end 
    if isfield(out,'if1'); out.if1=i2(out.if1); end 
    if isfield(out,'if2'); out.if2=i2(out.if2); end
    if isfield(varargin{2},'MPID'); out.MPID=varargin{2}.MPID;
    else; i3=feutil('mpid',elt);out.MPID(1:max(i2),:)=0;out.MPID(i2,:)=i3;
    end
   else;  out.MPID=feutil('mpid',elt);
   end
   if sdtdef('verm')>=7.0;
    st=intersect({'fs';'ifs';'f2';'if2';'f1';'if1';'MPID'},fieldnames(out));
    for j1=1:length(st);out.(st{j1})=int32(out.(st{j1}));end
   end
  end
else  % not new
  if ~isempty(node) % if node is defined in indices rather than numbers
   nInf=size(node,1); if isfinite(node(nInf,5)); nInf=nInf+1; end
   NNode=sparse(node(:,1)+1,1,1:length(node(:,1)));
   i2=find(LD2(:,1:4));LD2(i2)=full(NNode(LD2(i2)+1));
  end
  i2=cumsum(~isfinite(elt(:,1)));
  out=[1 1 size(elt,1) zeros(1,size(LD2,2)-3);LD2];
end

%% #GetEG - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% feutil('GetEG',node,NNode,elt_head,elt_elt)
elseif comstr(Cam,'eg'); [CAM,Cam]=comstr(CAM,3);

node  = varargin{carg};carg=carg+1;
NNode = varargin{carg};carg=carg+1;
elt   = varargin{carg};carg=carg+1;
[ElemF,EGID]= feutil('getelemf',elt(1,:),1);cEGI=2:size(elt,1);

% nodes linked to group

[i2,SEopt]=fe_super('node',ElemF); if SEopt(1,1)==2; i2 = elt(cEGI,i2);end
i2=i2(:); out.Node = node(NNode(sparse(i2+1,1,i2)~=0),:);

NNode=[]; NNode=sparse(out.Node(:,1)+1,1,1:size(out.Node,1));
nInf=size(out.Node,1)+1; 
out.Elt=elt;

% Basic Patch representation of Group - - - - - - - - - - - - - - - - - -

LD2=[1 1 size(elt,1)];
[Faces,SEopt]=fe_super('patch',ElemF);
if isempty(Faces) % in case there is no patch definition 
  [NodeLine,SEopt]=fe_super('line',ElemF);NodeLine=NodeLine(:)';
  if length(NodeLine)>1
  NodeLine=[NodeLine(1:length(NodeLine)-1);NodeLine(2:length(NodeLine))];
  else;NodeLine=[NodeLine;NodeLine];
  end
  Faces=NodeLine(:,NodeLine(1,:)&NodeLine(2,:))';
end

% Patch matrix : Object vertex matrix 
if ~isempty(Faces)&&(SEopt(1,1)==1||~isempty(cEGI))

  if SEopt(1,1)==2 % generic element or superelement
   i2=[]; i7=[];
   for j1=1:size(Faces,1) 
       i1 = Faces(j1,:);
       i1 = reshape(elt(cEGI,i1),length(cEGI),length(i1));
       i2(end+[1:size(i1,1)],1:size(i1,2))=i1; i7 = [i7 cEGI]; %#ok<AGROW>
   end % of loop on Patches linked to the element
  else % single superelement
    i2=fe_super('node',ElemF); i2=i2(:,1); i1=find(Faces);
    Faces(i1)=i2(Faces(i1));  i2=Faces; i7 = ones(1,size(i2,1));
  end
  % elimination of duplicated patches within same group
  if  size(i2,2)>2  && ~any(strcmp(fe_super('parent',ElemF), ...
    {'beam3','tria3','tria6','quad4','quadb'}))
   i1 = sort(i2,2); [i1,i3]=sortrows(i1);
   i4=find(~any(diff(i1),2));i5=ones(size(i2,1),1);i5(i3(i4))=0;i5(i3(i4+1))=0;
   i5=find(i5);i2=i2(i5,:);i7=i7(i5);
  elseif size(i2,2)==1; i2(:,2)=0;
  elseif size(i2,2)==2; i2(:,3)=0;
  end
end
i1=find(i2);i2(i1)=NNode(i2(i1)+1);  % switch to node indices

% treat faces with nodes at zero here

% output the various face matrices
i1 = sum((isfinite(i2)&i2),2)';

i3 = find(i1>2);  if ~isempty(i3); out.ps=i2(i3,:); out.es=i7(i3);
out.opt(1,1)=1; end
i3 = find(i1==2); if ~isempty(i3); out.p2=i2(i3,1:2); out.e2 = i7(i3); 
out.opt(1,2)=1;end
i3 = find(i1==1); 
if ~isempty(i3); out.p1=i2(i3,1); out.e1 = i7(i3); 
 if any(out.p1==0); error(1); end
 out.opt(1,3)=1;
else;out.opt(1,3)=0; 
end
if ~any(i1)
  sdtw('Group %i (%s) is not displayed',varargin{5},ElemF);
end
out.opt(1,5)=0;

%% #GetNormal [,Map][Elt,Node]' - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'normal'); [CAM,Cam]=comstr(Cam,7);

  model=[];
  % allow for selection
  if isfield(varargin{carg},'Elt')&&carg<nargin&&ischar(varargin{carg+1}) 
    model=varargin{carg};carg=carg+1; node=model.Node;
    elt=feutil(['selelt' varargin{carg}],model);
  else
    [carg,node,elt,ModelStack]=get_nodeelt(varargin,carg,ModelStack);
  end
  % Possible parameters
  if carg<=nargin&&isstruct(varargin{carg}); RunOpt=varargin{carg};carg=carg+1;
  else; RunOpt=struct;
  end
  [CAM,Cam,RunOpt.Map]=comstr('map',[-25 3],CAM,Cam);
  [CAM,Cam,RunOpt.Dir]=comstr('-dir',[-25 1],CAM,Cam);
  RunOpt.Silent=~isempty(Cam)&&Cam(end)==';'&&~sp_util('diag');
  if ~isfield(RunOpt,'method'); RunOpt.method='base';end
  [EGroup,nGroup]=getegroup(elt);

  if strncmpi(RunOpt.method,'infoatnode',10)
    % Element orientation MAP based on 'info','EltOrient'
    %feutil('getnormal',mo2,struct('method','InfoAtNode','keep',{{'v1x','v1y','v1z'}}))
    model.Elt=elt; C2=fe_mknl('init',model);
    eltid=feutil('eltidfix;',model);
    for jGroup=1:nGroup
     InfoAtNode=C2.GroupInfo{jGroup,7};
     [ElemF,i1,ElemP]=getegroup(model.Elt(EGroup(jGroup),:),jGroup);
     cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
     if ~isfield(RunOpt,'keep'); RunOpt.keep=InfoAtNode.lab;end
     if jGroup==1; out=InfoAtNode; out.EltId=eltid(cEGI); 
     elseif ~isempty(InfoAtNode); error('Not yet implemented');
     end
    end
    if nargout==0; fecom('showmap',out); clear out;end
    return
  end
  
  if ~isempty(node);NNode=sparse(node(:,1)+1,1,1:size(node,1));end
  %NNode=[];NNode(node(:,1)+1)=1:size(node,1);
  out=zeros(size(elt,1),3); out1=zeros(size(elt,1),3);
  out2=zeros(size(elt,1),1); out3=zeros(size(elt,1),3);
  if comstr(Cam,'node'); 
    if isempty(node);r3=spalloc(size(elt,1),1,1);RunOpt.AtNode=1;     
    else;r3=spalloc(size(elt,1),max(node(:,1)),1);RunOpt.AtNode=1;
    end
  else;r3=[]; RunOpt.AtNode=0;
  end
  for jGroup = 1:nGroup 
     %% loop on groups
     [ElemF,i1,ElemP]=feutil('getelemf',elt(EGroup(jGroup),:),jGroup);
     i1=fe_super('node',ElemP);
     cEGI = [EGroup(jGroup)+1:EGroup(jGroup+1)-1]';
     NodePos=reshape(full(NNode(elt(cEGI,i1)+1))',length(i1),length(cEGI));

     if comstr(ElemP,'quad4')||comstr(ElemP,'quadb') || ...
        comstr(ElemP,'tria3')||comstr(ElemP,'tria6')||comstr(ElemP,'quad9')
      if strcmp(ElemP,'quad4');RunOpt.Rot=8;else;RunOpt.Rot=0;end
      if size(elt,2)<RunOpt.Rot; elt(1,RunOpt.Rot)=0;end
      try  % Updated with SDT Call
       mo1=struct('Node',node,'Elt', ...
           elt(EGroup(jGroup):EGroup(jGroup+1)-1,:)); 
       if RunOpt.Silent||Silent;[un1,mo1.Elt]=feutil('eltidfix;',mo1);end
       if isfield(model,'Stack');mo1.Stack=model.Stack;end
       r1=RunOpt; r1.AtNode=0;r1.ejdet=ones(size(cEGI));
       mo1=feutilb('shellmap',mo1,r1);
       out(cEGI,:)=mo1.normal;out2(cEGI)=mo1.ejdet;
       %if RunOpt.AtNode;r3(cEGI)=1;
       out1(cEGI,:)=mo1.vertex;
      catch
       % Old OpenFEM MAP
      for jElt=1:length(cEGI) %loop on elements of group
        i2=NodePos(:,jElt); r1 = node(i2(i2~=0),5:7);
        r2 = mean(r1,1); out1(cEGI(jElt),1:3)=r2;
        r1=r1-r2(ones(size(r1,1),1),:);
        [bas,r1,out2(cEGI(jElt),1)]= basis(r1);
        if RunOpt.Rot&&abs(elt(cEGI(jElt),8))>1e-8
            r1=elt(cEGI(jElt),RunOpt.Rot)*pi/180;
            bas=bas*[cos(r1) -sin(r1) 0;sin(r1) cos(r1) 0;0 0 1];
        end
        out (cEGI(jElt),1:3)=bas(:,3)'; % normal
        if RunOpt.Dir;out(cEGI(jElt),1:3)=bas(:,RunOpt.Dir)';
        else;    out3(cEGI(jElt),1:3)=bas(:,1)'; % x direction
        end
      end % loop on elements of group
      end
     elseif any(strcmp(ElemP,{'beam1','bar1','beam3','line3'}))
     %% orient 2-d elements by pointing along axis
     %elseif comstr(ElemP,'beam3') 
     % [i2,i3]=min(std(node(:,5:7),[],1));r1=zeros(1,3);r1(i2)=1;% plane normal

      z0=std(node(:,5:7),[],1);[z0,ind]=min(z0);
      z0=full(sparse(1,ind,1,1,3)); %#ok<ACCUM>
      if isempty(RunOpt.Dir);RunOpt.Dir=3;end
      for jElt=1:length(cEGI) %loop on elements of group
        i2=NodePos(:,jElt);r1 = node(i2(i2~=0),5:7);
        r2 = mean(r1,1); out1(cEGI(jElt),1:3)=r2; % segment center
        if size(r1,1)>1;
          r1=diff(r1(1:2,:),[],1); out2(cEGI(jElt))=norm(r1);
          p=sp_util('basis',r1,z0); %basis(-diff(r1,[],1),[ 0 0 0],1);
          out(cEGI(jElt),1:3)=p(:,RunOpt.Dir)';
        % xxx for cbush should recover ori from cbush -2 call
        end
      end % loop on elements of group
      if ~isempty(r3);
         if size(NodePos,1)==3;r4=[1 1 1];else; r4=[1 1];end 
         r3=r3+sparse(cEGI(:,ones(size(i1))),elt(cEGI,i1), ...
           out2(cEGI)*r4,size(elt,1),max(node(:,1)));
      end
      continue; % Different strategy for r3
     elseif ~isempty(i1)
      %%
      for jElt=1:length(cEGI) %loop on elements of group
        i2=NodePos(:,jElt);r1 = node(i2(i2~=0),5:7);
        r2 = sum(r1,1)/size(r1,1); out1(cEGI(jElt),1:3)=r2;
      end % loop on elements of group
     end
     % centroid
     if ~isempty(r3); 
         r3=r3+sparse(cEGI(:,ones(size(i1))),elt(cEGI,i1), ...
           out2(cEGI(:,ones(size(i1)))),size(elt,1),max(node(:,1)));
     end
  end % loop on groups
  if RunOpt.Map  % New Normal map format

    if ~isempty(r3) % MAP at nodes

     i1=find(any(r3)); r3=r3(:,i1);  % r3=r3*diag(sum(r3).^(-1));
     r3=out'*r3;r3=r3*diag(sparse(sum(r3.^2).^(-.5)));
     out=struct('normal',r3','ID',i1(:),'opt',2);

    else % MAP at elements
     eltid=feutil('eltidfix;',elt);
     i1=find(eltid & any(out,2));
     out=struct('normal',out(i1,:),'ID',eltid(i1),'vertex',out1(i1,:),'opt',1);
    end
   
  else  % old nor/cg outputs

    if ~isempty(r3) % determine normal on each node
     i1=find(any(r3)); r3=r3(:,i1);   r3=r3*diag(sum(r3).^(-1));
     r3=out'*r3;r3=r3*diag(sparse(sum(r3.^2).^(-.5)));
     out=[i1(:) i1(:) r3'];
    end

  end

%% #GetCG #GetWeight - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% cg=feutil('GetCG',model) get element center
% weight=feutil('getweight',model) get integration at not weighting
elseif comstr(Cam,'cg')||comstr(Cam,'weight')

  % allow for selection
  if isfield(varargin{carg},'Elt')&&carg<nargin&&ischar(varargin{carg+1}) 
   model=varargin{carg};carg=carg+1; node=model.Node;
   elt=feutil(['selelt' varargin{carg}],model);
  else
    [carg,node,elt,ModelStack]=get_nodeelt(varargin,carg,ModelStack);
  end
  [EGroup,nGroup]=getegroup(elt);
  NNode=sparse(node(:,1)+1,1,1:size(node,1));
  out=zeros(size(elt,1),1); i3=0;
  if comstr(Cam,'weight'); i3=1; out=zeros(size(node,1),1); end
      
  for jGroup = 1:nGroup %loop on groups
     [ElemF,i1,ElemP]=feutil('getelemf',elt(EGroup(jGroup),:),jGroup);
     i1=fe_super('node',ElemP);
     cEGI = [EGroup(jGroup)+1:EGroup(jGroup+1)-1]'; 
     NodePos=reshape(full(NNode(elt(cEGI,i1)+1))',length(i1),length(cEGI));
     if i3
      r1=integrules(ElemP,-2); 
      if all(any(r1.xi,1)); r2=3;elseif any(r1.xi(:,2));r2=23;else;r2=13;end
      for jElt=1:length(cEGI) %loop on elements of group
        i2=NodePos(:,jElt); r1.nodeE = node(i2,5:7);
        of_mk('buildndn',r2,r1)
        out(i2,1)=out(i2,1)+r1.jdet.*r1.w(:,4);
      end
     else % take the mean
     for jElt=1:length(cEGI) %loop on elements of group
        i2=NodePos(:,jElt); r1 = node(i2(i2~=0),5:7);
        out(cEGI(jElt),1:3)=mean(r1,1);
     end % loop on elements of group
     end
  end % loop on groups

%% #GetNodeBas - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [node,bas,NNode]=feutil('GetNodeBas',model) in global coordinate system
elseif comstr(Cam,'nodebas')

model=varargin{carg};carg=carg+1;
if isfield(model,'bas') % empty .bas is acceptable for use in SE
  if isfield(model,'CheckedNodes')&&model.CheckedNodes(1)==1
   out=model.Node; out1=model.bas;
  else;[out,out1]=basis(CAM,model); % allow argument passing
  end
  if nargout>3 % Return local to global transformation
    if ~isempty(strfind(Cam,'nocgl')); out3=[];
    elseif ~isfield(model,'DOF')||isempty(model.DOF);out3=basis('trans l',out1,out);
    else;out3=basis('trans e',out1,out,model.DOF);
    end
  end
else; out=model.Node; out1=[]; out3=[];
end
if nargout>2; 
    if isempty(out);out2=[];else; out2=sparse(out(:,1),1,1:size(out,1));end
end

%% #getNode - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% node=feutil('GetNode (sel)',model) in global coordinate system
elseif comstr(Cam,'node')

[out1,out]=feutil(['findnode' CAM(5:end)],varargin{2:end});

%% #GetPl - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% pl=feutil('Getpl (MatId) -struct1',model)
elseif comstr(Cam,'pl');[CAM,Cam]=comstr(CAM,3);

[CAM,Cam,i2]=comstr('struct',[-25 31],CAM,Cam);
model=varargin{carg};carg=carg+1;

if comstr(Cam,'m_') % ('GetPl m_elastic')
 pl=model.pl; %feutil('getpl',model); 
 for j1=1:size(pl,1)
  [st,i1,i2]=fe_mat('type',pl(j1,2));
  if ~comstr(st,Cam);pl(j1,1)=0;end
 end
 out=pl(pl(:,1)~=0,:); return
end

i1=comstr(CAM,-1); % Possibly given matid
if isempty(i1);i1=feutil('matid',model);i1=unique(i1);i1=i1(i1~=0);end
pl=fe_mat('GetPl',model);
if carg<=nargin; mo1=varargin{carg};carg=carg+1; pl=fe_mat('GetPl',mo1); end
if ~isempty(pl)&&~isempty(i1);out=pl(ismember(pl(:,1),i1),:);
else; out=pl;end
if i2==-1 % display
  if isempty(pl);i2=i1;else;i2=setdiff(i1,pl(:,1));  end
  if ~isempty(i2); fprintf('Missing MatId %s\n',sprintf(' %i',i2));end
  feutilb('_writepl',pl); clear out
elseif i2
  out1=cell(size(out,1),1);
  for j1=1:size(out,1)
    [st,i3,i4]=fe_mat('type',out(j1,2));
    st=feval(st,'PropertyUnitType -cell',i4);st(:,2)={[]};
    i3=min(size(st,1),size(out,2));
    st(1:i3,2)=num2cell(out(j1,1:i3)');
    if i2==1; st=st(:,1:2)';out1{j1}=struct(st{:});% Struct
    elseif i2==2 ;out1{j1}=st; % Cell
    end
  end
  out=out1;if isscalar(out);out=out{1};end
end

%% #GetIl - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% il=feutil('Getil (ProId)',model)
elseif comstr(Cam,'il');[CAM,Cam]=comstr(CAM,3);

model=varargin{carg};carg=carg+1;
[CAM,Cam,i2]=comstr('struct',[-25 31],CAM,Cam);

if comstr(Cam,'p_') % ('Getil p_shell')
 if ~isfield(model,'il');il=[];out=[];return;else;il=model.il; end%feutil('getil',model); 
 for j1=1:size(il,1)
  [st,i1,i3]=fe_mat('typep',il(j1,2));
  if ~comstr(st,Cam);il(j1,1)=0;end
 end
 out=il(il(:,1)~=0,:); return
end

i1=comstr(CAM,-1); 
if isempty(i1);i1=feutil('proid',model);i1=unique(i1);i1=i1(i1~=0);end
il=fe_mat('GetIl',model);
if carg<=nargin; mo1=varargin{carg};carg=carg+1; il=fe_mat('GetIl',mo1); end
if ~isempty(il);out=il(ismember(il(:,1),i1),:);else; out=il;end

if i2==-1 % display
  if isempty(il);i2=i1;else;i2=setdiff(i1,il(:,1));  end
  if ~isempty(i2); fprintf('Missing ProId %s\n',sprintf(' %i',i2));end
  feutilb('_writeil',il); clear out
elseif i2
  out1=cell(size(out,1),1);
  for j1=1:size(out,1)
    [st,i3,i4]=fe_mat('typep',out(j1,2));
    st=feval(st,'PropertyUnitType -cell',i4);st(:,2)={[]};
    i3=min(size(st,1),size(out,2));
    st(1:i3,2)=num2cell(out(j1,1:i3)');
    if i2==1; st=st(:,1:2)';out1{j1}=struct(st{:});% Struct
    elseif i2==2 ;out1{j1}=st; % Cell
    end
  end
  out=out1;if isscalar(out);out=out{1};end
end


%% #GetDD : return detailed constitutive law structure - - - - - - - - - - - -
% dd=feutil('GetDD',[MatId ProId 3 1],model)
elseif comstr(Cam,'dd');[CAM,Cam]=comstr(CAM,3);
i1=varargin{carg};carg=carg+1;
model=varargin{carg};carg=carg+1; % Possibly case provided

il=feutil(sprintf('getil %i',i1(2)),model);
if isempty(il);
    il=p_solid(sprintf('dbval %i D3 3',i1(2)));model.il=il;
    i1(3:4)=[3 1]; % Default dims
end
[st,st1,i2]=fe_mat('typep',il(2));
if length(i1)>2; st2=sprintf('BuildConstit %i %i',i1(3:4));
else;st2='BuildConstit';
end
[constit,integ,i1,dd]= ...
    feval(st,st2,i1,fe_mat('getpl',model),model.il,model,varargin{carg:end});
out=dd; out1=constit; out2=integ;

%% #GetPro used to define element properties ->sdtweb fe_mat('getpro') - - - - - - -
elseif comstr(Cam,'pro');
  eval(iigui({'fe_mat',nargout},'OutReDir')); % now in fe_mat

%% #GetMat used to define material properties - - - - - - - - - - - - - - - - 
elseif comstr(Cam,'mat');
  eval(iigui({'fe_mat',nargout},'OutReDir')); % now in fe_mat

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
else;error('''Get%s'' unknown',CAM);
end % subcommand selection - - - - - - - 

%% #DefLocalToGlobal ---------------------------------------------------------
% [model,def]=feutil('DefLocalToGlobal',model,def) deformation in global
elseif comstr(Cam,'deflocaltoglobal')

model=varargin{carg};carg=carg+1;
def=varargin{carg};carg=carg+1;
if ~isfield(model,'DOF'); 
   model.DOF=feutil('getdof',model); 
end
[node,bas,NNode,r1]=feutil('getnodebas',model);
[r1,mdof]=basis('trans',bas,node,def.DOF);
i1=find(any(r1,2));r1=r1(i1,:);mdof=mdof(i1);
def.def=r1*def.def; def.DOF=mdof; 
model.Node(:,3)=0;

out=model; out1=def;


%% #Info ---------------------------------------------------------------------
elseif comstr(Cam,'info');  [CAM,Cam]=comstr(CAM,5);

if isempty(Cam) % guess info about what
 r1=varargin{carg};
 if isstruct(r1)||(~isempty(r1)&&~isfinite(r1(1))) 
   Cam=['elt' Cam];CAM=['elt' CAM];
 end
end

%% #InfoElt - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if comstr(Cam,'elt'); [CAM,Cam]=comstr(CAM,4); 
    
  elt = varargin{carg};carg=carg+1;
  if isstruct(elt)&&isfield(elt,'Elt'); elt=elt.Elt; end
  if isempty(elt)||ischar(elt); st=sprintf('\nModel %s is empty\n',CAM);
  else
   [EGroup,nGroup]=getegroup(elt);
   st=sprintf('\nModel %s with %i element(s) in %i group(s)\n', ...
        CAM,size(elt,1)-nGroup,nGroup);

   table=cell(nGroup,6);
   for jGroup = 1:nGroup %loop on element groups
    [ElemF,i1]= feutil('getelemf',elt(EGroup(jGroup),:),jGroup);
    [i3,SEopt]=fe_super('node',ElemF); i4 =fe_super('prop',ElemF);
    i3 = EGroup(jGroup+1)-EGroup(jGroup)-1;
    table{jGroup,1}=sprintf('%i',jGroup);table{jGroup,6}='';
    table{jGroup,4}=ElemF;
    if SEopt(1,1)==2
     table{jGroup,3}=sprintf('%i',i3);
     if i1(1)<0
      if i1(1)==-1 ; table{jGroup,2}='test wire frame';
      elseif i1(1)==-2; table{jGroup,2}='test/FEM node links';
      elseif i1(1)==-3; table{jGroup,2}='rotation interpolation links';
      else;            table{jGroup,2}=sprintf('EGID %i',i1(1));
      end
     elseif i1(1)~=jGroup; table{jGroup,2}=sprintf('%i',i1(1));
     else;table{jGroup,2}='';
     end
     if i4(1)~=0 && i4(1)<=size(elt,2) % MatId
      i5 = elt(EGroup(jGroup)+1:EGroup(jGroup+1)-1,i4(1));
      i5=unique(fix(i5));
      table{jGroup,5}=sdtw('_clip 20 1','%i ',i5);
     end
     if i4(2)~=0 && i4(2)<=size(elt,2) % proid
      i5=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
      if ~isempty(i5); i5 = elt(i5,i4(2));
       i5=unique(fix(i5));%i5 = find(sparse(fix(i5+1),1,i5+1))-1; 
       table{jGroup,6}=sdtw('_clip 20 1','%i ',i5);
      end
     end
    elseif SEopt(1,1)==1 % single superelement
     if i1(1)~=jGroup; table{jGroup,2}=sprintf('%i',i1(1)); end
     if SEopt(1,2)==0; table{jGroup,5}='single SE';
     else
      table{jGroup,5}='single SE';
      table{jGroup,6}=sprintf('%i',SEopt(1,2));
     end
    else
     table{jGroup,5}='Element Not Found';
    end
   end % of jGroup loop
   if nargout<2 % text display

     st2='Group ';i1=ones(size(table,1),1);
     st1=[st2(i1,:) char(table{:,1}) st2(i1,end)];
     i2=size(st1,2);
     for j1=1:size(table,1)
      if ~isempty(table{j1,2})
       st2=sprintf('(%s) ',table{j1,2});
       st1(j1,i2+[1:length(st2)])=st2;
      end
     end
     st2=' : '; st1=[st1 st2(i1,:) char(table{:,3})  ...
         st2(i1,1) char(table{:,4})];
     i2=size(st1,2);
     for j1=1:size(table,1)
      st2='';
      if ~isempty(table{j1,5})
       st2=sprintf('  MatId %s ',table{j1,5});
      end
      if ~isempty(table{j1,6})
       st2=sprintf('%s  ProId %s ',st2,table{j1,6});
      end
      if ~isempty(st2); st1(j1,i2+[1:length(st2)])=st2;end
     end
     st1(st1==0)=' ';
     st1=cellstr(st1);st1=sprintf('\n%s',st1{:});
     st=[st st1];

   end
  end % model is empty

 if nargout==0; disp(st);
 elseif nargout==1; out=st; 
 elseif nargout==2; out=st; out1=table;
 end

%% #InfoNode - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% feutil('infomode',struct('NodePolar',1))
% feutil('infomode',struct('iicom','chlink'))
elseif comstr(Cam,'node');[CAM,Cam]=comstr(CAM,5);

model=[];RunOpt.ImplicitNode=1;
[carg,FEnode,FEelt,ModelStack]=get_nodeelt(varargin,carg,ModelStack,RunOpt);

[CAM,Cam,RunOpt.Case]=comstr('case',[-25 3],CAM,Cam);
  opt = comstr(CAM,[-1 1]);
  [EGroup,nGroup,RunOpt.ElemP]=getegroup(FEelt);
  [i1,i2]=find(any(FEelt==opt(1),2));
  i1=sort([i1;EGroup(diff(EGroup)==1)']);
  if ischar(FEnode) % implicit nodes
   st1=', may be implicit';
   if isfield(model,'Node');
       FEnode=model.Node;i3=find(FEnode(:,1)==opt(1));
   end
  else;i3=find(FEnode(:,1)==abs(opt(1)));st1=' !!!!';
  end
  fprintf(1,'\n');
  st3='';
  try;
   if isfield(model,'nmap')&&isKey(model.nmap,'Map:Nodes')
    r2=model.nmap('Map:Nodes'); st3=r2(abs(opt(1))).safe;
    if ~ischar(st3);sdtw('_ewt','clean get name');st3='';end
   end
  end
  if isempty(i3)
   fprintf(1,'node %i not found%s\n\n',opt(1),st1);
   out=[];
  elseif isfield(InfoMode,'NodePolar')
   out=sprintf('node(p) %s [%i,%i %i %i,%.12g %.12g %.12g]\n\n', ...
       st3,basis('rect2cyl',FEnode(i3,:)));
   fprintf(1,out);
  elseif isfield(InfoMode,'iicom')
   st1='xyz';
   for j1=1:3
    sdtweb('_links',sprintf('iicom(''ch'',{''DOF'',%i.0%i})',FEnode(i3),j1), ...
      sprintf('%i%s',FEnode(i3),st1(j1)));
   end
   out='';
  else
   out=sprintf('node %s [%i,%i %i %i,%.12g %.12g %.12g]\n\n',st3,FEnode(i3,:));
   fprintf(1,out);
  end
  jGroup=0;
  for j1=1:length(i1)
    if max(find(EGroup<=i1(j1)))~=jGroup
      jGroup = max(find(EGroup<=i1(j1)));
      [ElemF,i3,ElemP]= getegroup(FEelt(EGroup(jGroup),:),jGroup);
      i3=fe_super('node',ElemF); 
    end
    i2=FEelt(i1(j1),:);
    if EGroup(jGroup)==i1(j1) % single superelement
      [i2]=fe_super('node',ElemF);
      if ~isempty(i2) && any(i2==opt(1))
       fprintf(1,'e %i (g %i : %s) single superelement\n',i1(j1),jGroup,ElemF);
      end
    elseif ~any(i2(i3)==opt(1))||~any(i2==opt(1))
    else
     fprintf(1,'e %i (g %i : %s), ',i1(j1),jGroup,ElemF);
     i2=FEelt(i1(j1),:);fprintf(1,' %i',i2(1:length(i3)));
     try;try;i4=feval(ElemF,'prop');catch;i4=feval(ElemF,'prop');end
         i4(i4==0|i4>size(i2,2))=[];
         if ~isempty(i4)&&i4(1)==i3(end)+1;
             fprintf(1,',%s ,',sprintf(' %i',i2(i4)));
         else;i4=i3;
         end
     catch;i4=i3;
     end
     fprintf(1,' %g',i2(i4(end)+1:end));
     fprintf(1,'\n')
    end    
  end
  if RunOpt.Case
   % fecom infonodeCase11003
   Case=fe_case(model,'getcase');
   for j1=1:size(Case.Stack,1)
    r1=Case.Stack{j1,3};
    if isstruct(r1)&&isfield(r1,'MasterSel')
        if ischar(r1.MasterSel); r2=feutil(['findnode' r1.MasterSel],model);
        else; warning('%s not handled',comstr(r1.MasterSel,-30));r2=[];
        end
        if isnumeric(r1.SlaveSel);r3=r1.SlaveSel;
        else; r3=feutil(['findnode' r1.SlaveSel],model);
        end
        if ismember(opt,r2)||ismember(opt,r3)
          disp(comstr(Case.Stack(j1,:),-30))
        end
        continue;
    end
    switch lower(Case.Stack{j1,1})
    case 'rbe3'
        otherwise 
            1;
    end
   end
  end
  if any(strcmp(RunOpt.ElemP,'SE')); eval('fesuper(''SeInfoNode'',model,opt);');end
%% #InfoMode - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'mode')
 InfoMode=varargin{2};
 
else;error('%s unknown',CAM);
end % 'info' subcommand selection - - - - - - - - - - - -

%% #Join ---------------------------------------------------------------------
elseif comstr(Cam,'join');   [CAM,Cam]=comstr(CAM,5);

model=[];elt=varargin{carg};if isfield(elt,'Elt'); model=elt;elt=model.Elt;end
[EGroup,nGroup]=getegroup(elt);
if comstr(Cam,'group')
   opt=comstr(CAM(6:length(CAM)),-1);st=[];
   if ~isempty(opt); opt=round(opt(:));opt=opt(opt>0&opt<=nGroup);end
   if isscalar(opt); disp(['Cannot join a single group' 7]); end
   st = feutil('getelemf',elt(EGroup(opt(1)),:),opt(1));
   i4=EGroup(opt(1));
elseif comstr(Cam,'all')
   [i1,st]=feutil('egid',elt);st=unique(st);
   for j2=1:length(st);elt=feutil(['join' st{j2}],elt);end
   if ~isempty(model);model.Elt=elt;out=model;else;out=elt;end
   return;
else; st = char(CAM); opt=[0 0];i4=[];
  if isempty(st)
      st = feutil('getelemf',elt(EGroup(1),:),1);
  end
end
if length(opt)>1
  i3=[]; i2=[];
  for jGroup = 1:nGroup %loop on element groups
   if strcmp(feutil('getelemf',elt(EGroup(jGroup),:),jGroup),st)
    if opt(1)==0 || any(opt==jGroup)
     i2=[i2 EGroup(jGroup)+1:EGroup(jGroup+1)-1]; %#ok<AGROW>
    else;i3=[i3 EGroup(jGroup):EGroup(jGroup+1)-1];%#ok<AGROW>
    end
    if isempty(i4); i4=EGroup(jGroup); end
   elseif any(opt==jGroup)
      fprintf('Cannot join group %i : different element type\n',jGroup)
      i3=[i3 EGroup(jGroup):EGroup(jGroup+1)-1];%#ok<AGROW>
   else % not a group to join
     i3=[i3 EGroup(jGroup):EGroup(jGroup+1)-1];%#ok<AGROW>
   end
  end % of loop on jGroup
end

if nGroup==1
elseif isempty(i2); fprintf('Found no groups of ''%s'' elements\n',st);
else
  i2=[i3(i3<min(i2)) i4 i2 i3(i3>min(i2))];
  elt=elt(i2,:);
end
if ~isempty(model);model.Elt=elt;out=model;else;out=elt;end


%% #MatComp ------------------------------------------------------------------
elseif comstr(Cam,'matcomp')  

r1 = varargin{carg};carg=carg+1; r2 = varargin{carg};carg=carg+1;
i1 = find(r1); i2=find(r2); 
if all(size(i1)==size(i2)) && all(i1==i2)
else; st='different sparsity'; end

d=r1-r2; 
[i3,i4] = find(r1-r2);i5=zeros(1,length(i3));
for j1=1:length(i3)
 i5(j1)= (d(i3(j1),i4(j1))^2/r1(i3(j1),i3(j1))*r1(i4(j1),i4(j1)));
end
if ~isempty(i3); st=sprintf('%g ',max(i5));
else;st='matrices are identical'; end
if nargout==1; out=st; else;fprintf('%s\n',st); clear out; end


elseif comstr(Cam,'nnode')  

r1=varargin{carg};carg=carg+1;
if isfield(r1,'Node'); r1=r1.Node;end
out=sparse(r1(:,1),1,1:size(r1,1));

%% #Node : manipulations of nodes positions in a model -----------------------
% feutil('Node [trans,rot,rb]',model,RO)
% trans x y z
% rot x1 x2 x3 n1 n2 n3 theta
% mir o x1 x2 x3 n1 n2 n3
% mir x ; mir y ; mir z
% mir eq a b c d  => aX+bY+cZ+d=0
% mir node  [x y z;x y z; ...]
% mir nodeid [nodeid list]
% feutil('nodeDefShift',model,def);
% rb 4*4 full rigidbody matrix

elseif comstr(Cam,'node'); [CAM,Cam]=comstr(CAM,5);
 
 model=[]; % Model can be assigned below
 if isfield(varargin{2},'Node') % model provided
  [carg,FEnode,FEel0,r1,ModelStack]=get_nodeelt(varargin,carg,ModelStack);
  if isempty(model);model=struct('Node',FEnode,'Elt',FEel0);end
 else % Node provided
  model=varargin{2}; carg=carg+1;
 end
 % Input parsing conv
 if carg>nargin||~isstruct(varargin{carg});RO=struct;
 else;RO=varargin{carg};carg=carg+1;
 end
 if ~isempty(CAM)&&strcmpi(CAM(end),';'); RO.Silent=1; CAM(end)=[]; end
 [RO,st,CAM]=paramEdit([ ...
  'trans(#%g#"Translation")' ...
  'rot(#%g#"Rotation")' ...
  'mir(#3#"Mirror")' ...
  'bas(#3#"Basis")' ...
  '-sel("all"#%s#"Selection on which displacement is applied")' ...
  'DefShift(#3#"shift def")' ...
  'Silent(#3#"Nor warning")' ...
  ],{RO,CAM}); Cam=lower(CAM);
 % end input
 
 %% #NodeTranslation
 [out,out1,CAM]= doNode(model,RO,CAM);



%% #Object -------------------------------------------------------------------
elseif comstr(Cam,'object'); [CAM,Cam] = comstr(CAM,7);

%% #ObjectBeamLine - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if comstr(Cam,'beamline'); [CAM,Cam] = comstr(CAM,9);

if carg<=nargin&&isstruct(varargin{carg})&&isfield(varargin{carg},'Node'); 
      model=varargin{carg}; carg=carg+1;
else; model=[];
end

if strfind(Cam,'-egid') % Get The EGID if possible
  i1=strfind(Cam,'-egid'); 
  [EGID,i2,i5,i3]=sscanf(CAM(i1+5:end),'%i',1);
  EGID=[0 EGID];
  CAM(i1:i1+i3+3)='';[CAM,Cam] = comstr(CAM,1);
else; EGID=[];
end
if isempty(Cam)||Cam(1)~='-'; st='beam1';
else;[st,i1,i2,i3]=sscanf(CAM(2:end),'%s',1);CAM=CAM(i3+1:end);
end

opt = comstr(CAM,[-1]); 

if isempty(opt)&&carg<=nargin;   opt=varargin{carg};carg=carg+1; end

opt=opt(:);opt(end+1)=0;
if isempty(opt); error('You need at least two nodes for a beam line');end
if carg<=nargin; i2=varargin{carg};carg=carg+1;else;i2=[1 1];end
if ~isnumeric(i2)||size(i2,1)~=1; i2=[1 1];end

elt=[Inf abs(st) EGID];
elt(2:length(opt),1:2+length(i2))=[opt(1:length(opt)-1) opt(2:length(opt)) ...
       ones(length(opt)-1,1)*i2];
i1=find(elt(2:end-1,1)==0&elt(3:end,2)==0 ...
        &elt(2:end-1,2)==elt(3:end,1));elt(i1+1,1)=elt(i1+1,2);
elt=elt(elt(:,1)&elt(:,2),:);

out=elt;
if ~isempty(model); 
  model.Elt(end+1:end+size(out,1),1:size(out,2))=out;
  out=model.Elt;
end

%% #ObjectCircle - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% model=feutil('ObjectCircle x y z r nx ny nz Nseg',model)
elseif comstr(Cam,'circle');[CAM,Cam]=comstr(CAM,7);
 
 r1=comstr(Cam,[-1]);
 if carg<=nargin&&isstruct(varargin{carg}); model=varargin{carg};carg=carg+1;
 else; model=struct('Node',[],'Elt',[]);
 end

 if length(r1)~=8; error(' x y z r nx ny nz Nseg must be provided');end
 r1(8)=ceil(r1(8));
 r2=linspace(0,2*pi,r1(8)+1);  p=basis(r1(5:7),[0 0 0],0);p=p(:,[2 3 1]);
 r2=r1(4)*p*[cos(r2);sin(r2);zeros(size(r2))];
 r2=r2'+ones(size(r2,2),1)*r1(1:3);

 [model.Node,i2]=feutil('addnode',model.Node,r2);
 i2=model.Node(i2,1);
 i1=[i2(1:end-1) i2(2:end) ones(size(i2,1)-1,1)*[1 1 0 0]];
 model.Elt(end+[1:size(i1,1)+1],1:6)=[Inf abs('beam1');i1];
 out=model;

%% #ObjectHoleInPlate - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%-- MC mod - begin
% CenterNode# Edge1Node# Edge2Node# r1 r2 nDiv1 nDiv2 nQuadrant  
elseif comstr(Cam,'holeinplate'); [CAM,Cam]=comstr(CAM,12);

 if carg<=nargin&&isstruct(varargin{carg}); model=varargin{carg};carg=carg+1;
 else; model=struct('Node',[],'Elt',[]);
 end
 
 NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
 if comstr(Cam,'o')
   opt = comstr(CAM(2:end),-1);
   node =[0 0 0 0 opt(1:3); model.Node(NNode(opt(4:5)),:)];
   opt=[0 opt(4:end)]; 
 else; opt = comstr(CAM,-1); node = model.Node(NNode(opt(1:3)),:);
 end
 if length(opt)~=8; 
   error('CenterNode# Edge1Node# Edge2Node# r1 r2 nDiv1 nDiv2 nQuadrant must be provided');
 end

% i1 basis for the hole, i3 combinations of basis vectors
i1 = [node(2,5:7)-node(1,5:7);node(3,5:7)-node(1,5:7)];
i1 = [i1(1,:)/norm(i1(1,:))*opt(4);i1(2,:)/norm(i1(2,:))*opt(5);i1];

i3 = [ones(opt(6)+1,1)   linspace(0,1,opt(6)+1)';
      linspace(1-1/opt(7),-1,2*opt(7))'  ones(2*opt(7),1);
      -ones(2*opt(6),1) linspace(1-1/opt(6),-1,2*opt(6))';
      linspace(-1+1/opt(7),1,2*opt(7))' -ones(2*opt(7),1);
      ones(opt(6),1)   linspace(-1+1/opt(6),0,opt(6))'];

%i3 = [atan2(i3(:,2),i3(:,1)) i3];
%i3 = [cos(i3(:,1)) sin(i3(:,1)) i3(:,[2 3])];
% Modified for regular angular spacing by M. Corus

ang_i3_1=linspace(0,pi/4,opt(6)+1);
ang_i3_2=linspace(0,pi/4,opt(7)+1);

ang_i3=[ang_i3_1 ...
        ang_i3_2(2:end)+pi/4 ang_i3_2(2:end)+pi/2 ...
        ang_i3_1(2:end)+3*pi/4 ang_i3_1(2:end)+pi ...
        ang_i3_2(2:end)+5*pi/4 ang_i3_2(2:end)+3*pi/2 ...
        ang_i3_1(2:end)+7*pi/4]';
i3 = [cos(ang_i3) sin(ang_i3) i3];



if opt(8)>=4
elseif opt(8)>=3.5; i3 = i3(1:size(i3,1)-  opt(6),:);
elseif opt(8)>=3  ; i3 = i3(1:size(i3,1)-  opt(6)-  opt(7),:);
elseif opt(8)>=2.5; i3 = i3(1:size(i3,1)-  opt(6)-2*opt(7),:);
elseif opt(8)>=2  ; i3 = i3(1:size(i3,1)-2*opt(6)-2*opt(7),:);
elseif opt(8)>=1.5; i3 = i3(1:size(i3,1)-3*opt(6)-2*opt(7),:);
elseif opt(8)>=1  ; i3 = i3(1:size(i3,1)-3*opt(6)-3*opt(7),:);
elseif opt(8)>=0.5; i3 = i3(1:size(i3,1)-3*opt(6)-4*opt(7),:);
end

i2 = [i3(:,1)*i1(1,:)+i3(:,2)*i1(2,:);i3(:,3)*i1(3,:)+i3(:,4)*i1(4,:)];
i2 = i2([1:size(i2,1)/2;size(i2,1)/2+1:size(i2,1)],:);
i2 = [[1:size(i2,1)]' zeros(size(i2,1),3) i2+node(ones(size(i2,1),1),5:7)];


[model.Node,i3]=feutil('AddNode',model.Node,i2);
i1 = size(i2,1);
if max(opt(4:5))<sqrt(sum((node(2,5:7)-node(1,5:7)).^2 ...
       +(node(3,5:7)-node(1,5:7)).^2))
    el0 = [2:2:i1-2; 1:2:i1-2; 3:2:i1 ; 4:2:i1]';
else
    el0 = [1:2:i1-2; 2:2:i1-2; 4:2:i1; 3:2:i1]';
end;
el0 = reshape(model.Node(i3(el0),1),size(el0,1),size(el0,2));
el0 = [Inf abs('quad4'); el0 ones(size(el0,1),1)*[1 1]];
model.Elt(end+[1:size(el0,1)],1:6)=el0;
out=model;

%-- MC mod - end

%% #ObjectArc - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% contributed by Mathieu Corus
% model=feutil('ObjectArc x_c y_c z_c x_1 y_1 z_1 x_2 y_2 z_2 Nseg obt',model)
% x_c y_c z_c : coordinate of the center of the arc
% x_1 y_1 z_1 : coordinate of the first end
% x_2 y_2 z_2 : coordinate of the second end
% Nseg        : number of division
% obt         : 1 if arc is to be obtuse, -1 if not.
elseif comstr(Cam,'arc');[CAM,Cam]=comstr(CAM,4);
 
  r1=comstr(Cam,[-1]);
  if carg<=nargin&&isstruct(varargin{carg}); model=varargin{carg};carg=carg+1;
  else; model=struct('Node',[],'Elt',[]);
  end
  if length(r1)~=11; error(' x_c y_c z_c x_1 y_1 z_1 x_2 y_2 z_2 Nseg obt must be provided');end
  
  cent=r1(1:3);
  no1=r1(4:6);
  no2=r1(7:9);
  N_seg=r1(10)+1;
  sh_obt=r1(11);
  v1=no1-cent;v2=no2-cent;
  r1=norm(v1);r2=norm(v2);
  v1=v1/r1;
  v2=v2/r2;
  
  c_th=v1*v2';
  th=acos(c_th);
  if sh_obt==-1 
    th=th-2*pi; 
  end;
  ori=cross(v1,v2);ori=ori/norm(ori);
  v2=cross(ori,v1);
  rr=linspace(r1,r2,N_seg)';
  th=linspace(0,th,N_seg)';
  Node=[(1:N_seg)' zeros(N_seg,3) cent(ones(N_seg,1),:)+...
      rr(:,[1 1 1]).*v1(ones(N_seg,1),:).*cos(th(:,[1 1 1]))+...
      rr(:,[1 1 1]).*v2(ones(N_seg,1),:).*sin(th(:,[1 1 1]))];
  Elt=[(1:N_seg-1)'  (1:N_seg-1)'+1 ones(N_seg-1,2) zeros(N_seg-1,2)];
  [model.Node,i2]=feutil('addnode',model.Node,Node);
  i2=model.Node(i2,1);
 i1=[i2(1:end-1) i2(2:end) ones(size(i2,1)-1,1)*[1 1 0 0]];
 model.Elt(end+[1:size(i1,1)+1],1:6)=[Inf abs('beam1');i1];
 out=model;
 
%% #ObjectDivide - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% model=feutil('ObjectDivide ...',model,selection,...)
elseif comstr(Cam,'divide');[CAM,Cam]=comstr(CAM,7);
% New call to refine cell (gvdr 15/05/2015), allows MPC generation and KN
model=varargin{carg};carg=carg+1;
sel=varargin{carg}; carg=carg+1;
[CAM,Cam,RO.mpc]=comstr('-mpc',[-25 3],CAM,Cam);
[CAM,Cam,RO.KN]=comstr('knownnew',[-25 3],CAM,Cam);
[CAM,Cam,RO.noSData]=comstr('-nosdata',[-25 3],CAM,Cam);
mo1=model;if isa(mo1,'v_handle');mo1=mo1.GetData;end
[i1,elt]=feutil(sprintf('FindElt %s',sel),mo1); mo1.Elt=elt;
R1=feutil(['DivideElt-getInput' CAM],mo1,varargin{carg:end});
R1.sel=i1; % Provide a subEGI
st=''; if RO.KN; st='KnownNew'; end
if RO.mpc; st=[st '-MPC']; end
[mo1,out1,out2]=feutil(sprintf('RefineCell-Replace%s;',st),model,R1); % out1 oldelt, mo1 new 
out=model; out.Node=mo1.Node; out.Elt=mo1.Elt; % keep other fields
if RO.mpc; out=fe_case(out,'stack_set',fe_case(mo1,'stack_get')); end
if ~RO.noSData; out=stack_set(out,'info','newcEGI',out2); end
% [model.Elt,FEel0]=feutil(['removeelt' varargin{carg}],model);carg=carg+1;
% mo1=model;if isa(mo1,'v_handle');mo1=mo1.GetData;end
% mo1.Elt=FEel0;mo1=feutil(['DivideElt' CAM],mo1,varargin{carg:end});
% model.Node=mo1.Node;
% model.Elt(end+[1:size(mo1.Elt,1)],1:size(mo1.Elt,2))=mo1.Elt;
% out=model;

%% #ObjectDisk - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% model=feutil('ObjectDisk x y z r nx ny nz NsegT NsegR',model)
elseif comstr(Cam,'disk');[CAM,Cam]=comstr(CAM,5);
 [CAM,Cam,nodeg]=comstr('-nodeg',[-25 3],CAM,Cam);
 r1=comstr(Cam,[-1]);
 if carg<=nargin&&isstruct(varargin{carg}); model=varargin{carg};carg=carg+1;
 else; model=struct('Node',[],'Elt',[]);
 end
 if length(r1)~=9; error(' x y z r nx ny nz NsegT NsegR must be provided');end
 r1(8:9)=ceil(r1(8:9)); r3=linspace(0,r1(4),r1(9)+1)';
 r2=linspace(0,2*pi,r1(8)+1);  p=basis(r1(5:7),[0 0 0],0);p=p(:,[2 3 1]);
 r2=[cos(r2);sin(r2);zeros(size(r2))];
 i1=length(r3)*length(r2);
 r2=[reshape(r3*p(1,:)*r2,i1,1)+r1(1) reshape(r3*p(2,:)*r2,i1,1)+r1(2) ...
  reshape(r3*p(3,:)*r2,i1,1)+r1(3)];
 
 [model.Node,i2]=feutil('addnode',model.Node,r2);
 i2=model.Node(i2,1);
 i1=length(r3);i1=[1:i1-1;2:i1;i1+2:2*i1;i1+1:2*i1-1];i1=i1(:);
 i3=[0:r1(8)-1]*length(r3);
 i1=i1(:,ones(r1(8),1))+i3(ones(size(i1,1),1),:);
 i1=reshape(i2(i1),4,numel(i1)/4)';i1(:,5:6)=1;
 if nodeg % transform degenerate quad4 into tria3
  i2=i1(:,1)==i2(1);
  i3=i1(~i2,:); i2=i1(i2,:); i2(:,4:5)=i2(:,5:6); i2(:,6)=0;
  model.Elt(end+[1:size(i1,1)+2],1:6)=[Inf abs('quad4');i3;Inf abs('tria3');i2];
 else; model.Elt(end+[1:size(i1,1)+1],1:6)=[Inf abs('quad4');i1];
 end
 out=model;

elseif comstr(Cam,'qdisk');[CAM,Cam]=comstr(CAM,6);
%% #ObjectQDisk : disk meshed as quad - -
% model=feutil('ObjectQDisk x y z crad hrad ND ran(i)',model)
% crad (center corner), hrad (outer corner) Nd (divide center edge), ran (diamters)
% http://www.truegrid.com/QualityMesh.html : idea of CLIMA4 corner mesh
%  feutil('ObjectQdisk 0 0 0  3 4 4  8 13  14');

 if carg<=nargin&&isstruct(varargin{carg}); model=varargin{carg};carg=carg+1;
 else; model=struct('Node',[],'Elt',[]);
 end
if any(Cam=='{') 
 %% #qdisk.2022 using fe_shapeoptim -22
   [st,RO]=sdtm.urnPar(CAM,'{rAnulus%g,Lc%g}{box%ug,h%ug,mat%sow}');

 r1=[0;abs(RO.rAnulus(:))];[r2,RO.ipos]=feutil(sprintf('refineline %g',RO.Lc),r1);
 RO.ipos(1)=[];
 RO.Nc=round(r1(2)*pi/2/RO.Lc/2);
 mo1=struct('Node',[(1:RO.Nc+1)'*[1 0 0 0] r2(1:RO.Nc+1)/sqrt(2)*[1 0 0]], ...
     'Elt',feutil(['objectBeamLine ' ,comstr(1:RO.Nc+1,-30)],[1 1]));
 mo1=feutil('extrude 0 0 1 0',mo1,r2(1:RO.Nc+1)/sqrt(2));
 [mo1.Node,i2]=feutil('addnode',mo1.Node,r2(RO.Nc+1:end)*[1 1 0]/sqrt(2));
 [mo1.Node,i2(:,2)]=feutil('addnode',mo1.Node,r2(RO.Nc+1:end)*[1 0 0]/sqrt(2));
 [mo1.Node,i2(:,3)]=feutil('addnode',mo1.Node,r2(RO.Nc+1:end)*[0 1 0]/sqrt(2));
 mo2=mo1; 
 mo2.Elt=feutil('addelt','quad4', ...
     [i2(1:end-1,1) i2(2:end,1) i2(2:end,2) i2(1:end-1,2) 2*ones(size(i2,1)-1,2);
      i2(1:end-1,3) i2(2:end,3) i2(2:end,1) i2(1:end-1,1) 3*ones(size(i2,1)-1,2) ...
     ]);
mo2=feutil(sprintf('divide %i 1',RO.Nc),mo2); mo2=feutil('addelt',mo2,mo1.Elt);
r2=feutilb('geolinetopo',mo2,struct('starts',i2(RO.ipos-RO.Nc,1:2)));
r3=feutilb('geolinetopo',mo2,struct('starts',i2(RO.ipos-RO.Nc,[3 1])));
r2=horzcat(r2{:}); r2(1,:)=[];r2=flipud([horzcat(r3{:});r2]);
RO.idiag=feutilb('geolinetopo',mo2,struct('starts',i2(RO.ipos([1 end])-RO.Nc,1)'));
RO.idiag=RO.idiag{1};

MAP=struct('ID',[],'normal',[]); % Build maps 
for j1=1:length(RO.rAnulus)
 n2=(0:size(r2,1)-1)';n2=n2/n2(end)*90;
 if RO.rAnulus(j1)>0; n2=RO.rAnulus(j1)*[cosd(n2) sind(n2) n2*0];
 else; n2=mo2.Node(r2(:,j1),5:7);% Do not move if negative
 end
 MAP.normal=[MAP.normal;n2-mo2.Node(r2(:,j1),5:7)];
 MAP.ID=[MAP.ID;r2(:,j1)];
end
RO.contour=r2; 
MAP.ID(end+1)=1;MAP.normal(end+1,:)=0;

RB.SelStep={'model.Elt=feutil(''selelt seledgeall'',model);'
 'model.Elt=feutil(''selelt selface'',model);'}; 
mo3=fe_shapeoptim('deform',mo2,MAP,RB); % Force radii

r3=feutilb('geolinetopo',mo2,struct('starts', ...
    [r2(1,3) 1;1 r2(end,3);r2(1,3) RO.idiag(end);RO.idiag(end) r2(end,3)]));
if 1==2
 MAP2=MAP; MAP2.normal=MAP.normal*0;
 i3=setdiff(r3{1},MAP.ID);n2=mo3.Node(i3,5:7);MAP2.normal=[MAP2.normal;n2*diag([1 0 0])-n2];MAP2.ID=[MAP2.ID;i3];
 i3=setdiff(r3{2},MAP.ID);n2=mo3.Node(i3,5:7);MAP2.normal=[MAP2.normal;n2*diag([0 1 0])-n2];MAP2.ID=[MAP2.ID;i3];
 mo3=fe_shapeoptim('deform',mo3,MAP2); % force straight edges
else
    i3=setdiff(r3{1},MAP.ID);mo3.Node(i3,6)=0; 
    i3=setdiff(r3{2},MAP.ID);mo3.Node(i3,5)=0; 
    RO.iedge=[r3{3} r3{4}];
end
if isfield(RO,'box')&&~isempty(RO.box) % Possibly add outer box
 mo1=mo3;
 mo2=mo1;mo2.Elt=feutil(['objectbeamline' sprintf(' %i',RO.iedge(:,1))],[4 4]);
 mo2=feutil('extrude 0  1 0 0',mo2, ...
     feutil(sprintf('Refineline %g',RO.Lc),[0 RO.box(1)-max(mo2.Node(:,5))]));
 mo1.Elt=feutil('addelt',mo1.Elt,mo2.Elt);mo2.Elt=mo1.Elt;
 mo2.Elt=feutil('SelElt seledge & innode {y>=}',mo2,max(mo2.Node(:,6)));
 mo2.Elt(2:end,3:4)=4;
 mo2=feutil('extrude 0 0 1 0',mo2,feutil(sprintf('Refineline %g',RO.Lc),[0 RO.box(2)-max(mo2.Node(:,6))]));
 mo2.Elt=feutil('addelt',mo1.Elt,mo2.Elt);
 mo3=mo2;
end
if isfield(RO,'h')&&~isempty(RO.h)
 mo3=feutil('extrude 0 0 0 1',mo3,feutil(sprintf('Refineline %g',RO.Lc),RO.h));
end
model=mo3;

else
 %% initial version from CLIMA
 r1=comstr(Cam,-1);RO.cRad=r1(4); RO.hRad=r1(5);RO.ND=r1(6);
 RO.rAnulus=r1(7:end);RO.NR=length(RO.rAnulus);
    % build the bolt head structure in ref plane, 1sq head, 1sq core
  mo2=struct('Node',[1 0 0 0 RO.cRad RO.cRad 0;2 0 0 0 RO.cRad -RO.cRad 0;
   3 0 0 0 -RO.cRad -RO.cRad 0; 4 0 0 0 -RO.cRad RO.cRad 0;
   5 0 0 0 RO.hRad RO.hRad 0;6 0 0 0 RO.hRad -RO.hRad 0;
   7 0 0 0 -RO.hRad -RO.hRad 0; 8 0 0 0 -RO.hRad RO.hRad 0],...
   'Elt',feutil('addelt','quad4', [1 5 6 2;4 8 5 1;3 7 8 4;2 6 7 3]));
  mo2=feutil(sprintf('divide%i %i',RO.ND*2,RO.NR),mo2); % refine mesh
  elt=mo2.Elt; mo2.Elt=feutil('addelt','quad4',[1 2 3 4]); % mesh interior
  mo2=feutil('addelt',feutil(sprintf('divide%i %i',RO.ND*[2 2]),mo2),elt);
  i1=reshape(elt(2:end,3),2*RO.ND,RO.NR,[]);
  
  % Transform external mesh into cylindrical shape: use polar coordinates
  mo2.Node=basis('rect2cyl',mo2.Node);
  mo2.Node(:,5)=mo2.Node(:,5)*.9/sqrt(2); %*min(RO.rAnulus); % main square at .8
  NNode=sparse(mo2.Node(:,1),1,1:size(mo2.Node,1));
  for j1=1:RO.NR % loop on external radii
   n1=mo2.Node(NNode(reshape(i1(:,j1,:),[],1)),5:7);
   n1(:,1)=RO.rAnulus(j1); r2=linspace(-180,180,size(n1,1)+1);r2(1)=[];
   [r3,i3]=sort(n1(:,2));r2=linspace(-180,180,length(r3)+1);
   if r3(1)<-179.9;r2(end)=[];else;r2(1)=[];end% Distribute angles evenly
   n1(i3,2)=r2(:);
   mo2.Node(NNode(reshape(i1(:,j1,:),[],1)),5:7)=n1;
  end
  % Generate screw by extruding in z dimension
  mo2.Elt=feutil('joinall',mo2.Elt);
  mo2.Node=basis('cyl2rect',mo2.Node);
  mo2.Node(:,5:7)=mo2.Node(:,5:7)+ones(size(mo2.Node,1),1)*reshape(r1(1:3),1,3);
  if isempty(model.Elt)&&isempty(model.Node);model=mo2;
  else
      model=feutil('addtest merge',model,mo2);
  end
end
if nargout==0;disp(comstr(RO,-30)); feplot(model);
else;  out=model;
end
if nargout>1;out1=RO;end

%% #ObjectAnnulus  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% model=feutil('ObjectAnnulus x y z r1 r2 nx ny nz NsegT NsegR',model)
elseif comstr(Cam,'annulus');[CAM,Cam]=comstr(CAM,8);
 
 r1=comstr(Cam,[-1]);
 if carg<=nargin&&isstruct(varargin{carg}); model=varargin{carg};carg=carg+1;
 else; model=struct('Node',[],'Elt',[]);
 end
 if length(r1)~=10; error(' x y z r1 r2 nx ny nz NsegT NsegR must be provided');end
 r1(9:10)=ceil(r1(9:10));r3=linspace(r1(4),r1(5),r1(10)+1)';
 r2=linspace(0,2*pi,r1(9)+1);  p=basis(r1(6:8),[0 0 0],0);p=p(:,[2 3 1]);
 r2=[cos(r2);sin(r2);zeros(size(r2))];
 i1=length(r3)*length(r2);
 r2=[reshape(r3*p(1,:)*r2,i1,1)+r1(1) reshape(r3*p(2,:)*r2,i1,1)+r1(2) ...
     reshape(r3*p(3,:)*r2,i1,1)+r1(3)];

 [model.Node,i2]=feutil('addnode',model.Node,r2);
 i2=model.Node(i2,1);
 i1=length(r3);i1=[1:i1-1;2:i1;i1+2:2*i1;i1+1:2*i1-1];i1=i1(:);
 i3=[0:r1(9)-1]*length(r3);
 i1=i1(:,ones(r1(9),1))+i3(ones(size(i1,1),1),:);
 i1=reshape(i2(i1),4,numel(i1)/4)';i1(:,5:6)=1;
 model.Elt(end+[1:size(i1,1)+1],1:6)=[Inf abs('quad4');i1];
 out=model;

%% #ObjectCylinder - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% model=feutil('ObjectCylinder x1 y1 z1 x2 y2 z2 R NsegR NsegT',model)
elseif comstr(Cam,'cylinder');[CAM,Cam]=comstr(CAM,9);
 
 r1=comstr(Cam,[-1]);
 if carg<=nargin&&isstruct(varargin{carg}); model=varargin{carg};carg=carg+1;
 else; model=struct('Node',[],'Elt',[]);
 end

 if length(r1)~=9; 
  error('x1 y1 z1 x2 y2 z2 R NsegR NsegT must be provided');
 end

 r2=linspace(0,2*pi,r1(8)+1);
 p=basis(r1(4:6)-r1(1:3),[0 0 0],0);p=p(:,[2 3 1]);
 r2=r1(7)*p*[cos(r2);sin(r2);zeros(size(r2))];
 r2=r2'+ones(size(r2,2),1)*r1(1:3);

 % build nodes and topology
 r3=zeros(size(r2,1)*r1(9),3);r4=ones(size(r2,1),1)*(r1(4:6)-r1(1:3))/r1(9);
 i2=size(r2,1);ind=[1:i2]';r3(ind,:)=r2;i1=ones((i2-1)*(r1(9)),6);
 for j1=1:r1(9); 
  i1((1:i2-1)+(i2-1)*(j1-1),[1 4 2 3])=[ ...
    ind(1:end-1)+i2*(j1-1) ind(1:end-1)+i2*j1 ...
    ind(2:end)+i2*(j1-1) ind(2:end)+i2*j1];
  r3(ind+j1*size(r2,1),:)=r2+r4*j1; 
 end

 [model.Node,i2]=feutil('addnode',model.Node,r3);
 i2=model.Node(i2,1); i1(:,1:4)=reshape(i2(i1(:,1:4)),size(i1,1),4);
 model.Elt(end+[1:size(i1,1)+1],1:6)=[Inf abs('quad4');i1];
 out=model;


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% #ObjectWireFrame - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
elseif comstr(Cam,'wireframe');
 
elt=[];  i2=[1 1];
while carg<=nargin

 opt=varargin{carg};carg=carg+1;
 opt=opt(:);opt(end+1)=0; %#ok<AGROW>
 if isempty(elt); elt=[Inf abs('beam1') 0 -1];
 else; elt(end+1,:)=elt(1,:); %#ok<AGROW>
 end

 opt=[opt(1:length(opt)-1) opt(2:length(opt)) ...
       ones(length(opt)-1,1)*i2];
 i1=find(opt(1:end-1,1)==0&opt(2:end,2)==0 ...
        &opt(1:end-1,2)==opt(2:end,1));opt(i1,1)=opt(i1,2);
 opt=opt(opt(:,1)&opt(:,2),:);
 elt(end+[1:size(opt,1)],1:size(opt,2))=opt;

end
out=elt;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%% #ObjectMass - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'mass'); 

if carg<=nargin&&isstruct(varargin{carg})&&isfield(varargin{carg},'Node'); 
      model=varargin{carg}; carg=carg+1;
else; model=[];
end

if length(Cam)>4 && ~isempty(strfind(Cam,'egid'))
 i2=strfind(Cam,'egid');
 [i1,i3,i4,i5]=sscanf(CAM(i2+4:end),'%d',1);
 [CAM,Cam]=comstr(CAM([1:i2-1 i2+3+i5:end]),1);
else;i1=0;
end

opt = comstr(comstr(CAM,5),[-1]);
if isempty(opt)&&carg<=nargin; opt=varargin{carg};carg=carg+1; end
opt=opt(:);opt=opt(opt~=0);

out=[Inf abs('mass1') 0]; if i1; out=[out i1];end % give EGID
opt=opt(:); 
if carg<=nargin
 r1=varargin{carg};carg=carg+1; r1=r1(:)';
 opt(:,1+[1:length(r1)])=r1(ones(size(opt,1),1),:);
end
out(1+[1:size(opt,1)],1:size(opt,2)) = opt;
if ~isempty(model); 
  model.Elt(end+1:end+size(out,1),1:size(out,2))=out;
  out=model.Elt;
end

%% #ObjectGet - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [model,node,carg]=feutil('objectget',carg,varargin{:});
elseif comstr(Cam,'get'); [CAM,Cam]=comstr(CAM,4);

 if isstruct(varargin{carg+2}); model=varargin{carg+2};carg=carg+1;
 else; model=struct('Node',[],'Elt',[]);
 end

 node=varargin{carg+2};carg=carg+1; 
 if size(node,2)==7; node=node(:,5:7);
 elseif size(node,2)==1 % nodeId given
  node=feutil('getnode',model,node); node=node(:,5:7);
 elseif size(node,2)==3
 else; error('Bad format for nodes');
 end
 out=model; out1=node; out2=carg;

%% #ObjectBeam - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% model=feutil('objectbeam',model,[Corners],dvx,dvy);
elseif comstr(Cam,'beam'); [CAM,Cam]=comstr(CAM,5);

 i3=comstr(Cam,[-1 1 1]);
 [model,node,carg]=feutil('objectget',carg,varargin{:});

 dv1=varargin{carg};carg=carg+1; 
 if isscalar(dv1); opt(1)=dv1+1;dv1=linspace(0,1,opt(1));
 else; opt(1)=length(dv1);
 end
 [r1,i1]=feutil('objectrefcube',dv1);r1(:,4)=0;
 rule=integrules('bar1',r1); 
 [model.Node,i2]=feutil('addnode',model.Node,rule.N*node);
 i2=model.Node(i2,1); i1=reshape(i2(i1),size(i1,1),size(i1,2));
 i1(:,3)=i3(1);i1(:,4)=i3(2);i1(1,6)=0;

 model.Elt(end+[1:size(i1,1)+1],1:6)=[Inf abs('beam1');i1];
 out=model;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%% #ObjectQuad - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% model=feutil('objectquad',model,[Corners],dvx,dvy);
elseif comstr(Cam,'quad'); [CAM,Cam]=comstr(CAM,5);

 i3=comstr(Cam,[-1 1 1]);
 [model,node,carg]=feutil('objectget',carg,varargin{:});

 dv1=varargin{carg};carg=carg+1; 
 if isscalar(dv1); opt(1)=dv1+1;dv1=linspace(0,1,opt(1));
 else; opt(1)=length(dv1);
 end
 dv2=varargin{carg};carg=carg+1; 
 if isscalar(dv2); opt(2)=dv2+1;dv2=linspace(0,1,opt(2));
 else; opt(2)=length(dv2);
 end

 [r1,i1]=feutil('objectrefcube',dv1,dv2); % Reference nodes and connectivity
 r1(:,3:4)=.5;rule=integrules('quad4',r1*2-1); 
 if size(node,1)==3 % orig dx,dy
     node=node(ones(4,1),:)+[0 0 0;node(2,:);sum(node(2:3,:));node(3,:)];
 end
 if isempty(model.Node) % keep order
   model.Node=[(1:size(rule.N,1))'*[1 0 0 0] rule.N*node];
   i2=model.Node(:,1);
 else;[model.Node,i2]=feutil('addnode',model.Node,rule.N*node);
     i2=model.Node(i2,1);
 end
  i1=reshape(i2(i1),size(i1,1),size(i1,2));
 i1(:,5)=i3(1); i1(:,6)=i3(2);
 model.Elt(end+[1:size(i1,1)+1],1:6)=[Inf abs('quad4');i1];
 out=model;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%% #ObjectHexa - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% model=feutil('objecthexa',[Oxyz;0Axyz;OBxyz;OCxyz],dvx,dvy,dvz);
elseif comstr(Cam,'hexa'); [CAM,Cam]=comstr(CAM,5);

 i3=comstr(Cam,[-1 1 1]);
 [model,node,carg]=feutil('objectget',carg,varargin{:});

 if size(node,1)==4 %[Oxyz;0Axyz;OBxyz;OCxyz]
  node=node(ones(8,1),:)+ ...
  [0 0 0;node(2,:);node(2,:);0 0 0;0 0 0;node(2,:);node(2,:);0 0 0]+ ...
  [0 0 0;0 0 0;node([3 3],:);zeros(2,3);node([3 3],:)]+...
  [zeros(4,3);node([4 4 4 4],:)];
 end

 dv1=varargin{carg};carg=carg+1; 
 if isscalar(dv1); opt(1)=dv1+1;dv1=linspace(0,1,opt(1));
 else; opt(1)=length(dv1);
 end
 dv2=varargin{carg};carg=carg+1; 
 if isscalar(dv2); opt(2)=dv2+1;dv2=linspace(0,1,opt(2));
 else; opt(2)=length(dv2);
 end
 dv3=varargin{carg};carg=carg+1; 
 if isscalar(dv3); opt(3)=dv3+1;dv3=linspace(0,1,opt(3));
 else; opt(3)=length(dv3);
 end

 [r1,i1]=feutil('objectrefcube',dv1,dv2,dv3); % Reference nodes and connectivity
 r1(:,4)=.5; rule=integrules('hexa8 node',r1*2-1); 
 [model.Node,i2]=feutil('addnode',model.Node,rule.N*node);
 i2=model.Node(i2,1); i1=reshape(i2(i1),size(i1,1),size(i1,2));
 i1(:,9)=i3(1); i1(:,10)=i3(2);
 model.Elt(end+[1:size(i1,1)+1],1:10)=[Inf abs('hexa8') 0 0 0 0;i1];
 out=model;


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%% #ObjectRefCube : generates nodes and connectivity of a reference cube - - -
elseif comstr(Cam,'refcube'); 

 dv1=varargin{carg};carg=carg+1; 
 if isscalar(dv1); opt(1)=dv1;dv1=linspace(0,1,opt(1));
 else; opt(1)=length(dv1);
 end
 if carg>nargin; dv2=[]; else; dv2=varargin{carg};carg=carg+1; end
 if isscalar(dv2); opt(2)=dv2;dv2=linspace(0,1,opt(2));
 else; opt(2)=length(dv2);
 end
 if carg>nargin; dv3=[]; else; dv3=varargin{carg};carg=carg+1; end
 if isscalar(dv3); opt(3)=dv3;dv3=linspace(0,1,opt(3));
 else; opt(3)=length(dv3);
 end
 if length(unique(dv1))~=length(dv1);sdtw('degenerate r division');end
 if length(unique(dv2))~=length(dv2);sdtw('degenerate s division');end
 if length(unique(dv3))~=length(dv3);sdtw('degenerate t division');end
 i1=sign(diff(dv1));
 if ~all(i1(1)==i1);error('Sign change in r division');end
 if ~isempty(dv2)
  i1=sign(diff(dv2));
  if ~all(i1(1)==i1);error('Sign change in s division');end
 end
 if ~isempty(dv3)
  i1=sign(diff(dv3));
  if ~all(i1(1)==i1);error('Sign change in t division');end
 end
if isempty(dv2) % 1D generation of a beam

  out=[dv1(:)*[1 0 0]];
  out1=[1:size(out,1)-1;2:size(out,1)]';

elseif isempty(dv3) % 2D generation of a quad4  - - - - - - - - - - -

 [r1,r2]=ndgrid(dv1,dv2);
 out=[r1(:) r2(:)*[1 0]];
 i4 = [1:opt(1)-1];i4=i4(ones(opt(2)-1,1),:)';
 i5=[0:opt(2)-2]*(opt(1));i5=i5(ones(opt(1)-1,1),:);
 i4=i4(:)+i5(:); 
 if OldRefQuad;i4 =[i4 i4+opt(1) i4+opt(1)+1 i4+1];
 else; i4 =[i4 i4+1 i4+opt(1)+1 i4+opt(1)]; % preserve x direction
 end
 out1=i4;

else  % 3D generation of an hexa8  - - - - - - - - - -
 [r1,r2,r3]=ndgrid(dv1,dv2,dv3);
 out=[r1(:) r2(:) r3(:)];

 % row of elements
 i4 = [1:opt(1)-1;2:opt(1)]';
 i4=[i4 i4(:,[2 1])+opt(1)];i4=[i4 i4+opt(1)*opt(2)];
 % face of elements
 i5=[0:opt(1):opt(1)*(opt(2)-1)-1];i5=i5(ones(size(i4,1),1),:);i5=i5(:);
 i6=[1:size(i4,1)]';i6=i6(:,ones(opt(2)-1,1)); i6=i6(:);
 i4=i4(i6,:)+i5(:,ones(size(i4,2),1));
 % volume of elements
 i5=[0:opt(1)*opt(2):opt(1)*opt(2)*(opt(3)-1)-1];
 i5=i5(ones(size(i4,1),1),:);i5=i5(:);
 i6=[1:size(i4,1)]';i6=i6(:,ones(opt(3)-1,1)); i6=i6(:);
 i4=i4(i6,:)+i5(:,ones(size(i4,2),1));

 out1=i4;

end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% #ObjectHoleInBlock - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% model=feutil('ObjectHoleInBlock 
%    x0 y0 z0  nx1 ny1 nz1  nx3 ny3 nz3 dim1 dim2 dim3 r nd1 nd2 nd3 ndr',
%    model)
elseif comstr(Cam,'holeinblock');[CAM,Cam]=comstr(CAM,12);
 opt=comstr(Cam,[-1]);
 %opt=str2num(Cam); %#ok<ST2NM>
 % param :
 x0=opt(1); y0=opt(2); z0=opt(3); n1=opt(4:6); n3=opt(7:9);
 dim1=opt(10); dim2=opt(11); dim3=opt(12); r=opt(13); 
 nd1=ceil(opt(14)/2);  nd2=ceil(opt(15)/2); nd3=opt(16);
 if length(opt)>16; ndr=opt(17); 
 else; ndr=fix((dim1/2-r)/(2*pi*r/(4*nd1+nd2)))+1; 
 end
 % cyl in cube geom
 model=struct('Node',...
  [1 0 0 0  0 0 -dim3/2;  2 0 0 0  dim1/2 0  -dim3/2;
   3 0 0 0  0 dim2/2  -dim3/2],'Elt',[]);
 model=feutil(sprintf('ObjectHoleInPlate 1 2 3 %g %g %g %g 4',r,r,nd1,nd2),model);
 model=feutil(sprintf('divide 1 %i',ndr),model);
 model=feutil(sprintf('Extrude %i 0 0 %g',nd3,dim3/nd3),model);
 % rotate / translate:
 n1=n1/norm(n1); n3=n3/norm(n3); n2=cross(n1,-n3); n2=n2/norm(n2);
 p=[n1(:) n2(:) n3(:)]';
 model.Node(:,5:7)=model.Node(:,5:7)*p;
 model.Node(:,5)=model.Node(:,5)+x0;
 model.Node(:,6)=model.Node(:,6)+y0;
 model.Node(:,7)=model.Node(:,7)+z0;
 
 out=model;
 
else;error('%s unknown',CAM);
end % subcommand selection - - - - - - - - - - - - - - -

% ----------------------------------------------------------------------------
%% #Optim : optimisation functions -------------------------------------------
elseif comstr(Cam,'optim'); [CAM,Cam] = comstr(CAM,6);

model=[]; % possible assign in get_nodeelt
[carg,FEnode,FEelt,FEel0,ModelStack]=get_nodeelt(varargin,carg,ModelStack);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% #OptimModel : removes unused nodes from FEnode - - - - - - - - - - - - - - -
if comstr(Cam,'model')

  [EGroup,nGroup]=getegroup(FEelt);
  if ~isempty(model);FEnode=model.Node;end
  NNode=sparse(FEnode(:,1),1,1:size(FEnode,1));

  i2 =[]; % nodes used
  for jGroup = 1:nGroup %loop on element groups
   ElemF= feutil('getelemf',FEelt(EGroup(jGroup),:),jGroup);
   i3=fe_super('node',ElemF);
   i3 = FEelt(EGroup(jGroup)+1:EGroup(jGroup+1)-1,i3);
   i2 = [i2;unique(i3(:))]; %#ok<AGROW>
  end % of jGroup loop
  % indices in FEnode of the retained nodes
  out=FEnode(ismember(FEnode(:,1),i2),:); out1=FEelt;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% #OptimNodeNum : optimizes model node numbering - - - - - - - - - - - - - - - 
elseif comstr(Cam,'nodenum')

  [EGroup,nGroup]=getegroup(FEelt);
  NNode=sparse(FEnode(:,1),1,1:length(FEnode(:,1)));
  i2 =[]; i4=size(FEnode,1); 
  for jGroup = 1:nGroup %loop on element groups
   ElemF= feutil('getelemf',FEelt(EGroup(jGroup),:),jGroup);
   i3=fe_super('node',ElemF);
   i3 = FEelt(EGroup(jGroup)+1:EGroup(jGroup+1)-1,i3);
   i5 = find(i3); i3(i5)=NNode(i3(i5));
   i5=[];
   for j1=1:size(i3,2)
    i5=[i5 (i3(:,ones(size(i3,2)-j1+1,1)*j1)-1)*i4+i3(:,j1:end)]; %#ok<AGROW>
   end; i2=[i2;i5(:)]; %#ok<AGROW>
  end % of jGroup loop

  k=reshape(sparse(i2,1,1,i4^2,1),i4,i4); k = k+k';
  i3=symrcm(k); 
  %subplot(211);spy(k); subplot(212);spy(k(i3,i3))
  FEnode=FEnode(i3,:);
  out=FEnode; out1=FEelt;

elseif comstr(Cam,'emptygroups')
%% #OptimEmptyGroups: removes empty group entries - - - - - - - - - - - - - -
if ~isempty(FEelt) % clean up empty groups
 i1=find(~isfinite(FEelt(:,1))); i2=find(diff(i1)==1);
 if ~isempty(i2); FEelt(i1(i2),:)=[]; end
 if ~isfinite(FEelt(end,1)); FEelt(end,:)=[]; end
end
out=FEelt;
  
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% #OptimDegen: transforms degenerated elements to their low level counterpart
% Developped with Eric.Monteiro@ensam.eu
elseif comstr(Cam,'degen')
    
RO.conv={'hexa','penta6',[1 5 6 4 8 7;1 6 2 4 7 3;1 4 8 2 3 7
    5 1 4 6 2 3;2 1 5 3 4 8;3 2 6 4 1 5], ...
    'penta15',[1 5 6 4 8 7 13 17 14 12 20 18 16 19 15
    1 6 2 4 7 3 17 14 9 12 18 10 19 15 11
    1 4 8 2 3 7 12 16 20 9 11 19 10 15 18
    5 1 4 6 2 3 13 12 16 17 9 11 14 10 15
    2 1 5 3 4 8 9 13 17 10 12 20 11 16 19
    3 2 6 4 1 5 10 14 18 11 9 17 12 13 20];
    'quad','tria3',[],'tria6',[]}; % Degenerate face
if isempty(model);out=RO;return;end

[EGroup,nGroup]=getegroup(model.Elt);
RO.new={};
for jGroup=1:nGroup
  EGroup=getegroup(model.Elt); % Update if removed element somwhere
  [ElemF,i1,ElemP]= getegroup(model.Elt(EGroup(jGroup),:),jGroup);
  i1=fe_super('node',ElemF);
  cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
  % 0/ find elements with repeated nodes
  i2=sum(diff(sort(model.Elt(cEGI,i1),2)',[],1)~=0,1)+1; %#ok<UDIM>
  cEGI(i2==length(i1))=[]; % Eliminate non degenerate
  if isempty(cEGI); continue;
  elseif strncmpi(ElemP,'hexa',4)
    % 1/ find degenerate face index  
    if strcmp(ElemP,'hexa8')
      RO.renum=RO.conv{1,3}; RO.in=1:length(RO.renum);RO.new{end+1}=RO.conv{1,2};
    else
      RO.renum=RO.conv{1,5}; RO.in=1:length(RO.renum);RO.new{end+1}=RO.conv{1,4};
    end
    elt=model.Elt(cEGI,1:8)';
    i3=fe_super('faces','hexa8')';
    i2=sum(double(diff(sort(reshape(elt(i3,:),size(i3,1),[])))~=0))+1;
    i2=reshape(i2,size(i3,2),[])'; % Number of distinct nodes in each face
    elt=model.Elt(cEGI,:); 
    for jElt=1:size(elt,1)
      elt(jElt,RO.in)=elt(jElt,RO.renum(find(i2(jElt,:)==2),:));
    end
    if strcmp(ElemP,'hexa8');elt=elt(:,[1:6 9:end]);
    else; elt=elt(:,[1:15 21:end]);
    end
    model.Elt(cEGI,1)=-1e100; RO.new{end+1}=elt;
  elseif strncmpi(ElemP,'quad',4) % xxx Missing tria6
    RO.renum=RO.conv{2,3}; RO.in=1:3;RO.new{end+1}='tria3';
    elt=model.Elt(cEGI,:); 
    for jElt=1:size(elt,1)
      i2=elt(jElt,1:4); if length(unique(i2))<3;elt(jElt,:)=Inf;continue;end
      i2=find(diff(i2)==0); 
      if isempty(i2); elt(jElt,:)=[elt(jElt,[1:3 5:end]) 0];
      else; elt(jElt,:)=[elt(jElt,[1:i2-1 i2+1:end]) 0];
      end
    end
    elt(~isfinite(elt(:,1)),:)=[];
    model.Elt(cEGI,1)=-1e100; RO.new{end+1}=elt;
    
  elseif strncmpi(ElemP,'tetra',4) % first attempt at flat tetra
    elt=model.Elt(cEGI,:); 
    elt(elt(:,1)==elt(:,2)|elt(:,2)==elt(:,3)|elt(:,3)==elt(:,4)|elt(:,1)==elt(:,4),:)=[];
    model.Elt(cEGI,1)=-1e100; 
    if ~isempty(elt);RO.new{end+1}=ElemP;RO.new{end+1}=elt;end

  else; warning('%i degenerate %s ignored',length(cEGI),ElemP);
  end
end
model.Elt(model.Elt(:,1)==-1e100,:)=[];
for j1=1:2:length(RO.new); 
    model.Elt=feutil('addelt',model.Elt,RO.new{j1+(0:1)});
end
out=model;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% #OptimEltCheck: checks the conformity of elements (twisted things, etc) - - - 
% check/fix ill-conditionned hexa8 / penta6 / quad4
% Contributed by Mathieu Corus
elseif comstr(Cam,'eltcheck')
  
  if nargin==1;test_fe_check;return;
  elseif nargin==2
    model=varargin{2}; if ~isstruct(model); error('Please input a structure with .Node and .Elt fields'); end;
    bypass=0;
  elseif nargin==3
    if isnumeric(varargin{2}) && isnumeric(varargin{3})
      model.Node=varargin{2};
      model.Elt=varargin{3};
      bypass=0;
    else
      model=varargin{2}; if ~isstruct(model); error('Please input a structure with .Node and .Elt fields'); end;
      bypass=varargin{3};
    end;
  elseif nargin==4
    if isnumeric(varargin{2}) && isnumeric(varargin{3}) && isnumeric(varargin{4})
      model.Node=varargin{2};
      model.Elt=varargin{3};
      bypass=varargin{4};
    else
      error('Bad argument types');
    end;
  else
    error('Too many input arguments');
  end;

NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));

%-- loop on hexa8 / penta6 / quad4 groups

[EGroup,nGroup]=getegroup(model.Elt);
ind_new_deg_quad=[];
new_deg_quad=[];
for jGroup=1:nGroup
  [ElemF,i1,ElemP]= getegroup(model.Elt(EGroup(jGroup),:),jGroup);
  cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
  switch ElemP
      case {'hexa8','penta6'}
       i1=feval(ElemP,'node');
       for jElt=cEGI
         model.Elt(jElt,i1)=elt_check(model.Elt(jElt,i1),'hexa8',model,NNode,bypass);
       end
      case 'quad4'
       i1=feval(ElemP,'node');
       for jElt=cEGI
        nodes=elt_check(model.Elt(jElt,i1),ElemP,model,NNode,bypass);
        if size(nodes,1)==2 
         ind_new_deg_quad=[ind_new_deg_quad;jElt];
         new_deg_quad=[new_deg_quad ; nodes];
        else;model.Elt(jElt,i1)=nodes;
       end
       end
  end
end

%-- add new lines in model.Elt if necessary --%

if ~isempty(ind_new_deg_quad)
  ind_new_deg_quad=ind_new_deg_quad+(0:(length(ind_new_deg_quad)-1))';
  for i1=1:length(ind_new_deg_quad);
    model.Elt=model.Elt([1:ind_new_deg_quad(i1) ind_new_deg_quad(i1):end],:);
    model.Elt(ind_new_deg_quad(i1)+(0:1),1:4)=new_deg_quad(2*i1-1:2*i1,:);
  end;
end;
out=model;


else;error('Optim%s unknown',CAM);
end % subcommand selection - - - - - - - - - - - - - -

% ----------------------------------------------------------------------------
%% #NewMid : structured generation of new nodes ------------------------------
% node,elt,face
% .Cyl=1 use cylindrical coordinates
% .NodeNormals =struct('ID','normal')
elseif comstr(Cam,'newmid')
    
 node=varargin{carg};carg=carg+1;
 elt=varargin{carg};carg=carg+1;
 face=varargin{carg};carg=carg+1;
 if isstruct(elt)
   NNode=[];RunOpt=elt;
   elt=vertcat(RunOpt.I{:,2})';if isempty(elt);out=node;out1=RunOpt;return;end
 else;elt=elt';  % new nodes based on elements
  NNode=sparse(node(:,1),1,1:size(node,1));
  if carg<=nargin; RunOpt=varargin{carg}; carg=carg+1; 
  else; RunOpt=struct; 
  end
  r1={'KnownNew','Cyl'}; i1=find(~ismember(r1,fieldnames(RunOpt)));
  for j1=1:length(i1); RunOpt.(r1{i1(j1)})=0; end
  %RunOpt=struct('KnownNew',0,'Cyl',0);
 end
 out1=zeros(size(elt,2),0);
 for j1=1:length(face) % face or edge
  i1=sort(reshape(elt(face{j1},:),size(face{j1},1),[]),1)';
  [i1,i2,i3]=unique(i1,'rows');
  if ~isempty(NNode);r1=NNode(i1');else;r1=i1';end
  if ~RunOpt.Cyl;r1=reshape(node(r1,5:7),size(i1,2),[],3);
  else;r1=reshape(basis('rect2cyl',node(r1,5:7)),size(i1,2),[],3);
     i2=find(diff(r1(:,:,2),[],1)<-200);r1(1,i2,2)=r1(1,i2,2)-360;
     i2=find(diff(r1(:,:,2),[],1)>200);r1(1,i2,2)=r1(1,i2,2)+360;
     i2=find(r1(1,:,1)<epsl);r1(1,i2,2)=r1(2,i2,2); % on axis
     i2=find(r1(2,:,1)<epsl);r1(2,i2,2)=r1(1,i2,2);
  end
  r1=squeeze(sum(r1,1))/size(i1,2);if size(r1,2)==1;r1=r1';end
  if RunOpt.Cyl;r1=basis('cyl2rect',r1);end
  if isfield(RunOpt,'KnownEdges')
    [i4,i5]=ismember(sort(node(i1),2),RunOpt.KnownEdges(:,1:2),'rows');
    if any(i4); error('Need implement reuse of known edge');
    end
  end
  if RunOpt.KnownNew;
      i2=size(node,1)+[1:size(r1,1)]';
      node(i2,1)=max(node(:,1))-i2(1)+1+i2;
      node(i2,5:7)=r1;
  else;[node,i2]=feutil('AddNode',node,r1);
  end
  i2=node(i2,1);
  i2=reshape(i2(i3),size(face{j1},2),[])';
  if ~isempty(i2);out1(:,end+[1:size(i2,2)])=i2; end% node number for each face
 end
 if isscalar(face)&&length(face{1})==2; 
     RunOpt.EdgeI=struct('Edge',elt','NodeId',out1);
 end
 if isfield(RunOpt,'I')
   for jGroup=1:size(RunOpt.I,1)
       if ~isempty(RunOpt.I{jGroup,1})
           i2=1:size(RunOpt.I{jGroup,2},1);i3=out1(i2);out1(i2)=[];
           RunOpt.I{jGroup,3}=reshape(i3,RunOpt.I{jGroup,3},[])';
       end
   end
   out1=RunOpt;
 end
 out=node;

    
% ----------------------------------------------------------------------------
%% #Orient -------------------------------------------------------------------
elseif comstr(Cam,'orient'); [CAM,Cam]=comstr(CAM,7);

  [carg,node,elt,el0,ModelStack]=get_nodeelt(varargin,carg,ModelStack);
if comstr(Cam,'permute')
 %% #OrientPermute xyz  
 i1=comstr(Cam,8);i1=abs(i1)-abs('x')+1;i1=i1(:);
 node(:,i1+4)=node(:,5:7);model.Node=node;
 [Case,ind]=stack_get(model,'case');
 for j1=1:size(Case,1)
  r1=Case{j1,3};
  for j2=1:size(r1.Stack,1)
    r2=r1.Stack{j2,3};
    switch lower(r1.Stack{j2,1})
        case 'fixdof';
            if isstruct(r2);error('Not implemented');
            else; r2=fix(r2)+i1(round(rem(r2,1)*100))/100;
            end
        otherwise
            if isfield(r2,'DOF')
                r2.DOF=fix(r2.DOF)+i1(round(rem(r2.DOF,1)*100))/100;
            end
    end
    r1.Stack{j2,3}=r2;
  end
  Case{j1,3}=r1;
 end
 model.Stack(ind,:)=Case;
 out=model;
 return;
end
%% #Orient standard
  [EGroup,nGroup]=getegroup(elt);
  NNode=sparse(node(:,1)+1,1,1:size(node,1));

  % give element groups and/or nodes
  RunOpt=struct('neg',1,'Group',[]);
  [CAM,Cam,i1]=comstr('-neg',[-25 3],CAM,Cam); if i1;RunOpt.neg=-1;end
  [CAM,Cam,i1]=comstr('nlist',[-25 3],CAM,Cam);
  [CAM,Cam,RunOpt.RmWarp]=comstr('-rmwarp',[-25 31],CAM,Cam);
  [CAM,Cam,RunOpt.RmFlat]=comstr('-rmflat',[-25 31],CAM,Cam);
  if i1; RunOpt.GiveList=2; RunOpt.List=[];
  else;  RunOpt.GiveList=0;
  end
  [CAM,Cam,i1]=comstr('gauss',[-25 3],CAM,Cam);
  if i1;RunOpt.GiveList=1; end
  [CAM,Cam,RunOpt.Back]=comstr('-back',[-25 3],CAM,Cam);
  i1=strfind(Cam,'n'); if isempty(i1); i1=length(CAM)+1; end
  RunOpt.Group=comstr(Cam(1:i1-1),[-1]); i4=0;
  RunOpt.Warped=[]; RunOpt.Flat=[];
  RunOpt.jdet=zeros(size(elt,1),3);
  RunOpt.silent=any(Cam==';');
  if isempty(RunOpt.Group) % orientation based on k-matrix sign - - - - - - 

  i5=zeros(size(elt,1),1);
  RunOpt.Flat=zeros(size(elt,1),1);

  for jGroup = 1:nGroup %loop on groups
     [ElemF,i1,ElemP]=feutil('getelemf',elt(EGroup(jGroup),:),jGroup);
     i1=fe_super('node',ElemP);
     cEGI = [EGroup(jGroup)+1:EGroup(jGroup+1)-1]';
     %i2(DOF) i3(node)
     switch ElemP
      case {'tria3','t3p'}; ElemP=ElemF; [in2,in1]=tria3('flip');  i2=6; i3=3;
      case {'tria6','t6p'}; ElemP=ElemF; [in2,in1]=tria6('flip');  i2=12;i3=6;
      case {'quad4','q4p'}; ElemP=ElemF; [in2,in1]=quad4('flip');  i2=8; i3=4;
      case {'quadb','q8p'}; ElemP=ElemF; [in2,in1]=quadb('flip');  i2=16; i3=8;
      case {'quad9','q9a'}; ElemP=ElemF; [in2,in1]=quad9('flip');  i2=18;i3=9;
      case 'q5p';                        [in2,in1]=q5p('flip');    i2=10;i3=5;
      case {'tetra4','flui4'} ;          [in2,in1]=tetra4('flip'); i2=12;i3=4;
      case 'tetra10';                    [in2,in1]=tetra10('flip');i2=30;i3=10;
      case {'penta6','flui6'};           [in2,in1]=penta6('flip'); i2=18;i3=6;
      case 'penta15' ;                   [in2,in1]=penta15('flip');i2=45;i3=15;
      case {'hexa8','flui8'} ;           [in2,in1]=hexa8('flip');  i2=24;i3=8;
      case {'hexa20'} ;               [in2,in1]=hexa20('flip');  i2=60;i3=20;
      case 'pyra5';[in2,in1]=pyra5('flip');  i2=15;i3=5;
      case 'pyra13';[in2,in1]=pyra13('flip');  i2=39;i3=13;
      case {'beam1','bar1','mass1','mass2'}
      otherwise 
       if ~Silent;sdtw('Orientation test not implemented for %s',ElemP);end
     end
     if strcmp(ElemP(end),'b'); ElemP=ElemP(1:end-1); end
     % this should be rewritten using calls of BuildNDN 
     if any(strcmp(ElemP,...
             {'pyra5','penta6','penta15','tetra4','tetra10','hexa8','hexa20'}))
      if RunOpt.GiveList==1; EltConst=integrules(ElemP,-3); 
      else; EltConst=integrules(ElemP,'node'); 
      end; % if change keep same for GiveList
      NodePos=reshape(full(NNode(elt(cEGI,i1)+1)),length(cEGI),[])';
      for jElt=1:length(cEGI) %loop on elements of group
        EltConst.nodeE = node(NodePos(:,jElt),[5:7 4]);%xxx 4? also penta?
        of_mk('buildndn',3,EltConst);
        RunOpt.jdet(cEGI(jElt),1:3)= ...
            [min(EltConst.jdet) max(EltConst.jdet) mean(EltConst.jdet)];
        if all(EltConst.jdet>=0); %k(1)<0  
          r1=(max(EltConst.nodeE)-min(EltConst.nodeE));r1(4)=[];
          r1=[(EltConst.jdet'*EltConst.w(:,4)) max(r1)^3];
          if ~r1(1)||r1(1)<1e-8*r1(2); 
            RunOpt.Flat(cEGI(jElt))=1; 
          end
        elseif RunOpt.GiveList==2
         RunOpt.List(end+[1 2],[1:length(i1)+1])= ...
          [cEGI(jElt) elt(cEGI(jElt),i1);0 EltConst.jdet(:)'];
        elseif all(EltConst.jdet<=0)
         if RunOpt.RmWarp==2 % Do not remove negative volume
         else 
          elt(cEGI(jElt),in1)=elt(cEGI(jElt),in2);i5(cEGI(jElt))=1; 
         end
        else; 
         RunOpt.Warped(end+1,1)=cEGI(jElt);
        end
      end % loop on elements of group
     elseif any(strcmp(ElemP,{'t3p','t6p','q4p','q8p','q5p'}))
      for jElt=1:length(cEGI) %loop on elements of group
        r1 = node(NNode(elt(cEGI(jElt),i1)+1),[5:7 4]);
        if ~any(r1(:,2))&&any(r1(:,3)); r1=r1(:,[1 3 2 4]); end
        k=of_mk(ElemP,int32([i2*(i2+1)/2 0 0 0 0 0 0 0 0]), ...
           int32([0 0 i2 i3 1 0 1 2]),[7800 2e11 .3],r1);
        if k(1)<0;
          elt(cEGI(jElt),in1)=elt(cEGI(jElt),in2);i5(cEGI(jElt))=1;
        end
      end % loop on elements of group
     end % elements to run
  end % loop on groups
  RunOpt.FlatEltInd=find(RunOpt.Flat)~=0;
  RunOpt.Flat=elt(RunOpt.FlatEltInd,:);
  if ~isempty(RunOpt.Warped)&&~RunOpt.Back;
   if ~RunOpt.silent;
    disp(elt(RunOpt.Warped,:));
    try;
     fprintf('%s %s '';\n','cf.sel=''eltind',sprintf('%i ',RunOpt.Warped))
    end
   end
   if RunOpt.RmWarp
     if ~RunOpt.silent; fprintf('Removing warped elements');end
     elt(RunOpt.Warped,:)=[];
   elseif ~RunOpt.Back; fprintf('Found warped elements');
   else; error('There were some warped volumes, list above');
   end
  end
  if ~isempty(RunOpt.Flat)&&~RunOpt.Back;
   if RunOpt.RmFlat
     if ~RunOpt.silent; fprintf('Removing flat elements');end
     elt(RunOpt.FlatEltInd,:)=[];
   elseif ~RunOpt.Back; fprintf('Found flat elements');
   else; 
    RunOpt.Flat(:,min(find(~any(RunOpt.Flat,1))):end)=[]; %#ok<MXFND>
    disp(RunOpt.Flat);error('There were some zero volume elements, list above');
   end
  end
  
  if RunOpt.GiveList==2; out=RunOpt.List;
  else; 
   out=elt; out1=find(i5); i5=nnz(i5);
   if nargout>2; out2=RunOpt;end
   if RunOpt.Back;out2=RunOpt;
   elseif i5&&~any(Cam==';'); fprintf('%i elements reoriented\n',i5);
   end
  end
  % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  else % surface orientation based on nodes given ...

  out1=[];
  onode = comstr(Cam(i1+1:end),[-1]);
  if isempty(onode)&&carg<=nargin && isa(varargin{carg},'double')
    onode=varargin{carg};carg=carg+1;
    if isempty(onode)&&carg<=nargin
      onode=varargin{carg};carg=carg+1;
    end
    if size(onode,2)==7; onode=onode(:,5:7);end
  elseif rem(numel(onode),3)==0 
    onode=reshape(onode,3,numel(onode)/3)';
  end
  if ~any(Cam==';');fprintf('Orientation nodes \n');disp(onode);end
  [EGroup,nGroup]=getegroup(elt);
  NNode=sparse(node(:,1)+1,1,1:length(node(:,1)));
  for jGroup = RunOpt.Group(RunOpt.Group<=nGroup) %loop on groups
   [ElemF,i1,ElemP]=feutil('getelemf',elt(EGroup(jGroup),:),jGroup);
   i1=fe_super('node',ElemP);
   cEGI = [EGroup(jGroup)+1:EGroup(jGroup+1)-1]';
   NodePos=reshape(NNode(elt(cEGI,i1)+1,:),length(cEGI),[])';
   if comstr(ElemP,'quad4')||comstr(ElemP,'quadb')|| ...
      comstr(ElemP,'tria')
      for jElt=1:length(cEGI) %loop on elements of group
        i2=NodePos(:,jElt);i2(i2==0)=[];
        r1 = node(i2,5:7); r2 = mean(r1);r1=r1-r2(ones(size(r1,1),1),:);
        [bas,r1] = basis(r1);
        out1(end+1,1:6)=[r2 bas(:,3)']; %#ok<AGROW>
        r1 = onode-r2(ones(size(onode,1),1),:);
        [r2,i3]=sort(sum(r1.^2,2));

        % change orientation of surfaces based on normal
        if r1(i3(1),:)*bas(:,3)*RunOpt.neg<0  
          if comstr(ElemP,'quadb')
            elt(cEGI(jElt),[1:8])=elt(cEGI(jElt),[3 2 1 4 6 5 8 7]);
          elseif comstr(ElemP,'quad4')
            if ~diff(elt(cEGI(jElt),[3:4])) % degenerate
             elt(cEGI(jElt),[1 3 4])=elt(cEGI(jElt),[3 1 1]);
            else;elt(cEGI(jElt),[1 3])=elt(cEGI(jElt),[3 1]);
            end
          elseif comstr(ElemP,'tria6')
            elt(cEGI(jElt),1:6)=elt(cEGI(jElt),[1 3 2 6 5 4]);
          else % tria3
            elt(cEGI(jElt),[1 3])=elt(cEGI(jElt),[3 1]);
          end
          i4 = i4+1; out1(end,4:6)=[-bas(:,3)']; 
        end
      end
   elseif comstr(ElemP,'beam3')
      for jElt=1:length(cEGI) %loop on elements of group
        i2=NodePos(:,jElt);
        r1 = node(i2(1:2),5:7);
        r2 = mean(r1);r1=r1-r2(ones(size(r1,1),1),:);
        bas = basis(r1(2,:),[0 0 0],1);
        out1(end+1,1:6)=[r2 bas(:,1)']; %#ok<AGROW>
        r1 = onode-r2(ones(size(onode,1),1),:);
        [r2,i3]=sort(sum(r1'.^2));
        % change orientation of lines based on axis
        if r1(i3(1),:)*bas(:,1)*RunOpt.neg<0
            elt(cEGI(jElt),1:2)=elt(cEGI(jElt),[2 1]);
            out1(end,4:6)=bas(:,1)';  i4 = i4+1; 
        end
      end
   elseif any(strcmp(ElemP,{'beam1','bar1'}))
   else
       error('I don''t know how to orient %s',ElemP);
   end % supported
  end % loop on groups
  if ~any(Cam==';');fprintf('%i  element(s) reoriented\n',i4);end
  out=elt;

  end %of orientation type

%% #Refine : REFINE functions ------------------------------------------------
elseif comstr(Cam,'refine'); [CAM,Cam] = comstr(CAM,7);

%% #RefineLine MaxLength  - - - - - - - - - - - - - -
if comstr(Cam,'line'); [CAM,Cam]=comstr(CAM,5);
 [CAM,Cam,tol]=comstr('-tolmerge',[-25 2],CAM,Cam);
 opt=comstr(CAM,[-1 1]);
 r1=varargin{carg};carg=carg+1;
 r1=r1(:)';out=num2cell(r1);
 for j1=1:length(r1)-1
  r2=linspace(r1(j1),r1(j1+1),abs(ceil(diff(r1(j1:j1+1))/opt))+1);
  out{j1}=r2(1:end-1);
 end
 out{end}=r1(end);
 out=horzcat(out{:})';
 if ~isempty(tol) % tolMerge nodes that are too close
   out(1+find(diff(out(1:end-1))<tol))=[];
 end
 if nargout>1;[i1,out1]=ismember(r1,out); end% Return initial positions
  
%% #RefineBeam - - - - - - - - - - - - - -
elseif comstr(Cam,'beam'); [CAM,Cam]=comstr(CAM,5); RO=struct;
 [CAM,Cam,RO.MergeNew]=comstr('-mergenew',[-25 3],CAM,Cam);
 [CAM,Cam,RO.Pin]=comstr('-pin',[-25 3],CAM,Cam);
 if RO.MergeNew; RO.KN=''; else; RO.KN='KnownNew'; end
 if comstr(Cam,'uni'); [CAM,Cam]=comstr(CAM,4);
  %% #RefineBeamUni -3
  % should be a call to object divide, when new implementation performed
  model=varargin{carg}; carg=carg+1;
  if carg>nargin; varg={'groupall'}; else; varg=varargin(carg:end); end
  [model,el1,i1]=feutil(['ObjectDivide-noSData' RO.KN CAM],model,varg{:});
  if RO.Pin % handle pin flags
   mo1=struct('Node',model.Node,'Elt',el1); % removed elements
   elt=feutil('selelt eltname beam',mo1);
   if size(elt,2)<10; elt(end,10)=0; end
   if ~isempty(elt)&&any(any(elt(isfinite(elt(:,1)),9:10))) % pin flags to handle
    el2=elt(isfinite(elt(:,1)),[1 2 9 10]);
    for j1=1:2 % handle pins
     i2=el2(:,j1+2)~=0;r2=unique(el2(:,j1+2)); r2(r2==0)=[];
     if ~isempty(r2)
      for j2=1:length(r2) % loop on pin values
       i3=el2(ismember(el2(:,j1+2),r2(j2)),j1);
       model.Elt(i1(ismember(model.Elt(i1,j1),i3)),j1+8)=r2(j2);
      end
     end
    end % pin 1/2
   end
  end % pin handling
  out=model;
  
 else
  %% #RefineBeamMaxLength -3
  opt=comstr(CAM,[-1 1]);
  model=[];
  [carg,FEnode,FEel0,elt,ModelStack]=get_nodeelt(varargin,carg,ModelStack);
  NNode=sparse(FEnode(:,1),1,[1:length(FEnode(:,1))]); elt=[];
  if carg<=nargin&&ischar(varargin{carg});
   [RO.Keep,FEel0]=feutil(['RemoveElt' varargin{carg}],model);carg=carg+1;
   if isempty(FEel0);error('Empty selection for refine');end
  end
  [EGroup,nGroup]=getegroup(FEel0);
  for jGroup = 1:nGroup %loop on element groups
   [ElemF,i1,ElemP]= feutil('getelemf',FEel0(EGroup(jGroup),:),jGroup);
   if comstr(ElemP,'beam1')&&~strcmp(ElemF,'celas')
    elt=[elt;FEel0(EGroup(jGroup),:)]; %#ok<AGROW>
    cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;
    for jElt=cEGI
     l = norm([1 -1]*FEnode(NNode(FEel0(jElt,1:2)),5:7));
     i1=fix((l-epsl)/opt)+2;
     if i1>1 % clean sdtweb('_tracker',700)
      r1=linspace(0,1,i1)'; i2=full(NNode(FEel0(jElt,1:2)));
      r1([1 end])=[]; i1=[];
      if ~isempty(r1) % add new intermediate nodes
       node=FEnode(i2(1)*ones(length(r1),1),5:7) ...
        +r1*diff(FEnode(i2,5:7));
       [FEnode,i1]=feutil(sprintf('addnode%s',RO.KN),FEnode,node);i1=FEnode(i1,1);
      end
      i1=[FEnode(i2(1));i1;FEnode(i2(2))]; %#ok<AGROW>
      el1=[i1(1:end-1) i1(2:end) ones(size(i1,1)-1,1)*FEel0(jElt,3:end)];
      if size(el1,2)<10; el1(end,10)=0; end
      if comstr(ElemF,'beam') % handle pin flags
       for j1=1:2 % loop p1/p2
        if any(el1(:,8+j1))
         r2=FEel0(jElt,j1+[0 8]);el1(:,j1+8)=0;
         el1(ismember(el1(:,j1),r2(1)),j1+8)=r2(2);
        end
       end
      end
      if size(elt,2)<size(el1,2)
       if any(any(el1(:,9:10))); elt(end+1:size(el1,1),1:size(el1,2))=el1;
       else; elt=[elt;el1(:,1:size(elt,2))]; %#ok<AGROW>
       end
      else; elt=[elt;el1]; %#ok<AGROW>
      end
     else; elt=[elt;FEel0(jElt,:)]; %#ok<AGROW>
     end
    end
   else;elt=[elt;FEel0(EGroup(jGroup):EGroup(jGroup+1)-1,:)]; %#ok<AGROW>
   end % type of element
  end % of loop on groups
  model.Node=FEnode; model.Elt=elt;
  if isfield(RO,'Keep');model.Elt=feutil('addelt',RO.Keep,model.Elt);end
  out=model;
 end

elseif comstr(Cam,'toquad');[CAM,Cam]=comstr(CAM,7);
%% #RefineToQuad : triangles to 3 quad, penta to 3 hexa, hexa to 8 hexa
% see sdtweb t_femesh('RefDivide2')
R1.tria3=struct('edge',{{[4 1 2;5 2 3;6 3 1],[.5]}},'face',{{[7 1:3],[]}}, ...
     'Elt',{{'quad4',[1,4,7,6;2,5,7,4;3,6,7,5]}});
R1.quad4=struct('edge',{{[5 1 2;6 2 3;7 3 4;8 4 1],[.5]}},'face',{{[9 1:4],[]}}, ...
     'Elt',{{'quad4',[1 5 9 8;5 2 6 9;9 6 3 7;8 9 7 4]}});
R1.tetra4=struct('edge',{{[5 1 2;6 2 3;7 3 1;8 1 4;9 2 4;10 3 4],[.5]}}, ...
    'face',{{[11 1 2 3;12 2 3 4;13 1 3 4;14 1 2 4],[]}}, ...
    'volume',{{[15 1 2 3 4],[]}},'Elt',{{'hexa8', ...
    [1,5,11,7,8,14,15,13;2,6,11,5,9,12,15,14;3,7,11,6,10,13,15,12;
     8,4,9,14,13,10,12,15]}});
R1.penta6=struct('edge',{{[7 1 2;8 2 3;9 3 1;10 4 5;11 5 6;12 6 4;
    13 1 4; 14 2 5; 15 3 6],[.5]}}, ...
    'face',{{[16 1 2 3 3;17 4 5 6 6; 18 1 4 6 3; 19 1 2 5 4; 20 2 3 6 5],...
    [1/3 1/3 1/3 0;1/3 1/3 1/3 0; 1/4 1/4 1/4 1/4; 1/4 1/4 1/4 1/4;
    1/4 1/4 1/4 1/4]}},'volume',{{[21 1:6],[]}},'Elt',{{'hexa8', ...
    [1,7,16,9,13 19 21 18 ; 13 19 21 18 4,10,17,12;
    2,8,16,7,14 20 21 19;14 20 21 19 5,11,17,10;
    3,9,16,8,15 18 21 20;15 18 21 20 6,12,17,11]}});
R1.hexa8=struct('edge',{{[9 1 2;10 2 3;11 3 4;12 4 1;13 1 5;14 2 6;
      15 3 7;16 4 8;17 5 6; 18 6 7; 19 7 8; 20 8 5],[.5]}}, ...
    'face',{{[21,1,4,3,2;22,1,5,8,4;23,1,2,6,5;24,5,6,7,8;25,2,3,7,6;
       26,3,4,8,7],[]}},'volume',{{[27 1:8],[]}},'Elt',{{'hexa8', ...
    [1 9 21 12 13 23 27 22;9 2 10 21 23 14 25 27;
     12 21 11 4 22 27 26 16;21 10 3 11 27 25 15 26
     13 23 27 22 5 17 24 20;23 14 25 27 17 6 18 24
     22 27 26 16 20 24 19 8;27 25 15 26 24 18 7 19]}});
R1.pyra5=struct('edge',{{[6 1 2; 7 2 3; 8 3 4; 9 4 1; 10 1 5; 11 2 5; 12 3 5;13 4 5],.5}},...
 'face',{{[14 1 4 3 2; 15 4 1 5 5; 16 1 2 5 5; 17 2 3 5 5; 18 3 4 5 5],...
 [1/4 1/4 1/4 1/4;1/3 1/3 1/3 0;1/3 1/3 1/3 0; 1/3 1/3 1/3 0; 1/3 1/3 1/3 0]}},...
 'volume',{{[19 1:5],[1/6 1/6 1/6 1/6 1/3]}},'Elt',{{'hexa8',...
 [1 6 14 9 10 16 19 15; 6 2 7 14 16 11 17 19; 7 3 8 14 17 12 18 19;8 4 9 14 18 13 15 19];
 'pyra5',[10 16 11 19 5; 11 17 12 19 5; 12 18 13 19 5;13 15 10 19 5 ]}});
if ~isempty(strfind(Cam,'-data')); out=R1; % output structure
else
 mo1=varargin{carg}; carg=carg+1;
 out=feutil(sprintf('RefineCell %s',CAM),mo1,R1);
 if ~isempty(strfind(Cam,'-keep'))
  st=setdiff(fieldnames(mo1),fieldnames(out));
  for j1=1:length(st); out.(st{j1})=mo1.(st{j1}); end
 end
end

elseif comstr(Cam,'cell'); [CAM,Cam]=comstr(CAM,5);
%% #RefineCell : Generic element wise refinement procedure
% model=feutil('RefineCell',model,RO);
ModelStack={}; SetStack={};
[carg,FEnode,FEelt,FEel0,ModelStack]=get_nodeelt(varargin,carg,ModelStack);
if carg<=nargin; RO=varargin{carg}; carg=carg+1; else; RO=struct; end
if ~isempty(strfind(Cam,';'))&&~sp_util('diag'); RO.Silent=';'; else; RO.Silent=''; end
[RO,st,CAM]=paramEdit([ ...
  'KnownNew(#3#"dont check new nodes")' ...
  'allElt(#3#"use all elts")' ...
  'replace(#31#"0 : refined area, 1 append, 2 also preserve Stack,prop ")' ...
  'mpc(#3#"") keepEP(#3#"") keepSets(#3#"") given(#3#"") mpcALL(#3#"")' ...
  'noEOri(#3#"force base behavior, ignore EltOrient")' ...
  ],{RO,CAM}); Cam=lower(CAM);
% Check ElemP input in RO, use ToQuad if nothing given or allElt required
R1={'tria3','quad4','penta6','hexa8'};
if ~RO.given&&((RO.allElt&&~all(ismember(R1,fieldnames(RO))))||...
  (isempty(intersect(fieldnames(RO),R1))))
 R1=feutil('refineToQuad-data');
 r1=fieldnames(R1); r1=setdiff(fieldnames(R1),fieldnames(RO));
 for j1=1:length(r1);RO.(r1{j1})=R1.(r1{j1}); end;
end
% Check if an element set is given for subselection and face orientation
if ~isfield(RO,'set'); RO.set=[]; 
elseif isstruct(RO.set); RO.set=RO.set.data;
elseif ischar(RO.set); % get set in stack
 RO.set=stack_get(struct('Stack',{ModelStack}),'set',RO.set,'get');
 RO.set=RO.set.data;
end
% allow on-the-fly element selection here
if isempty(RO.set)&&(isfield(RO,'sel')||carg<=nargin)
 if isfield(RO,'sel'); sel=RO.sel; else; sel=varargin{carg}; carg=carg+1; end
 if isnumeric(RO.sel); RO.subEGI=RO.sel(:);
 else; RO.subEGI=feutil('findelt', ...
         struct('Node',FEnode,'Elt',FEelt,'Stack',{ModelStack}),sel);
 end
elseif ~isempty(RO.set)
 eltid=feutil('eltid;',FEelt); EEid=sparse(eltid+1,1,1:length(eltid));
 RO.subEGI=full(EEid(RO.set(:,1)+1));
end
if isfield(RO,'subEGI')&&(RO.mpc||RO.mpcALL)
 if RO.mpcALL; RO.MPCedge=[]; RO.MPCface=[]; RO.MPCvolume=[];
 else
 % setdiff edge/faces of submodel to orig model
 mo1=struct('Node',FEnode,'Elt',FEelt); [u1,mo1.Elt]=feutil('EltIdFixx;',mo1);
 el1=feutil('addsetFaceId-get',mo1,'f1','selface');
 %mo2=mo1; [mo2.Elt,elt]=feutil('removeelt eltid',mo2,unique(RO.set(:,1)));
 mo2=mo1; [mo2.Elt,elt]=feutil('removeelt eltind',mo2,RO.subEGI);
 el2=feutil('addsetFaceId-get',mo2,'f2','selface');
 el2.data=setdiff(el2.data,el1.data,'rows');
 if ~isempty(el2.data) % found faces: volume model
  mo1=stack_set(mo1,'set','f1',el2);
  elt=feutil('selelt setname f1',mo1);
  RO.MPCface=elt(isfinite(elt(:,1)),1:4);
  elt=feutil('selelt setname f1 & seledgeall',...
   stack_set(mo1,'set','f1',el2));
 else % shell model, do setdiff on edges instead of faces
  el1=feutil('selelt selface & seledge',mo1);
  el2=feutil('selelt selface & seledge',mo2);
  if size(el1,2)<size(el2,2); el1(end,size(el2,2))=0;
  elseif size(el2,2)<size(el1,2); el2(end,size(el1,2))=0;
  end
  elt=setdiff(el2(isfinite(el2(:,1)),:),el1(isfinite(el1(:,1)),:),'rows');
 end
 RO.MPCedge=elt(isfinite(elt(:,1)),1:2);
 end
elseif ~isfield(RO,'subEGI'); RO.subEGI=[];
elseif isempty(RO.subEGI); sdtw('_nb','no element to refine');
 out=struct('Node',FEnode,'Elt',FEelt); return
end
% Sequential addition of nodes per edge, then face, then volume.
% This allows avoiding generating duplicate nodes
% Cell refinement is defined for an element type ElemP
% giving additional nodes by edge, then face, then volume in increasing
% index. New nodes are computed using an operator performing
% weighted sums of initial cell coordinates. In no weights are given,
% arithmetic average is used.
% edge,face,volume share the same format as 
% {[newId cellNodeId_forAv] [weights]}
% new elements are added with the base topology given
% if an element selection or set is given, refine is limited to selection
% for non symmetric refinement an orientation has to be given through a subset
% of reference faces, with field shift and faces, that provides a way to
% match the reference cell nodes and the current cell nodes
% for several reference faces, it is needed to give field orient
% that also provides the face ordering sequence in the reference cell
[EGroup,nGroup,names]=getegroup(FEelt);
RO.Points=cell(nGroup,1); RO.pro=cell(nGroup,1); 
RO.Edges=cell(nGroup,2); RO.Faces=cell(nGroup,2); RO.Volumes=cell(nGroup,2);
RO.eCoef=cell(nGroup,1); RO.fCoef=cell(nGroup,1); RO.vCoef=cell(nGroup,1);
r2={'Edges','edge';'Faces','face';'Volumes','volume'};
%% Initialize group wise data
for jGroup=1:nGroup
 [ElemF,i1,ElemP]= getegroup(FEelt(EGroup(jGroup),:),jGroup);
 if strcmpi(ElemP,'SE'); continue; end
 try; epn=feval(ElemP,'node');
 catch; epn=1:str2num(regexprep(ElemP,'\D',''));
 end
 try; epp=feval(ElemP,'prop'); 
 catch; epp=epn(end)+1:3; 
 end
 cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
 if ~isempty(RO.subEGI) % restrain cEGI
  [cEGI,isub]=intersect(cEGI,RO.subEGI);
 end
 if ~isfield(RO,ElemP)||isempty(RO.(ElemP)) 
  if isempty(RO.Silent)&&~Silent; fprintf('Refine %s ignored\n',ElemP); end
   continue;
 end
 if isempty(cEGI); continue; end
 r1=RO.(ElemP);i2=epn; %feval(ElemP,'node');
 i1=epp; %feval(ElemP,'prop'); 
 if max(i1)>size(FEelt,2);FEelt(1,max(i1))=0;end
 RO.pro{jGroup}=FEelt(cEGI,i1(i1~=0)); % store mat/pro data
 RO.Points{jGroup}=FEelt(cEGI,epn); %feval(ElemP,'node')); % store initial nodes
 if isfield(r1,'shift')&&size(RO.set,2)>1 % unsymm refine
  for j1=1:size(r1.faces,1) % robust no doublons on face list
   [u1,i1]=unique(r1.faces(j1,:)); r1.faces(j1,length(i1)+1:end)=0;
  end
  % in this case need to handle elements by selected face
  r3=unique(abs(RO.set(ismember(RO.set(:,1),eltid(cEGI)),2:end)),'rows'); % faces id to treat
  % EltNum will store the renumbered base cell topology 
  if ~isfield(r1,'EltNum') % initialize once (multi-groups)
   r4=num2cell(repmat(size(r1.faces,1),size(r3,2),1)); r1.EltNum=cell(r4{:});
  end
 else; r3=0;
 end
 % loop on face Id (one pass if distinction not necessary)
 for j3=1:size(r3,1)
  ic=0;
  if r3(j3) % renum r1 to corresponding face
   %% Init renumbering for cell orientation
   % renumbering is based on an identification between the reference cell nodes
   % and a rotated reference cell to match the good FaceId
   n1=sort(unique(r1.faces(:))); n1(n1==0)=[]; n1(:,2)=0; r4=r3(j3,:);
   in5=1:size(r1.faces,1);
   % place ref nodes in cur node index
   % shift: reference FaceId to identify, 
   % r5 is a reordered r1.faces with coherence regarding current set face selection scheme 
   if length(r1.shift)>1
    % find face order that preserves common features, not necesserily a "shift"
    r5=r1.faces;
    r5(r1.shift,:)=r1.faces(r4,:); 
    i5=find(all(r1.faces,2)); i5=intersect(1:size(r5,1),i5);
    r5(setdiff(i5,r1.shift),:)=r1.faces(setdiff(i5,r4),:);
    [u1,in5]=ismember(r1.faces,r5,'rows');
   else % with one faceid to work, we can shift an preserve all properties
    r6=r1.shift(1)-r4(1); % shift offset
    r5=circshift(r1.faces,r6); % shift to get in the reference cell
    if ~isequal(r5==0,r1.faces==0) % issue with unsym face pattern
     i5=all(r1.faces,2); % stra limited to full faces references
     r5=r1.faces; r6=r1.shift(1)-r4(1)+sum(i5(r1.shift(1):r4(1))==0);
     r5(i5,:)=circshift(r1.faces(i5,:),r6);
    end
   end
   for j4=1:size(r3,2) % loop on reference faces to identify nodes for renum.
    % if r6<0; r6=size(r1.faces,1)+r6; end
    %i5=r5(remi(r3(j3,j4)+r6,size(r1.faces,1)),:); % slave face
    i5=r5(r1.shift(j4),:); % slave face
    n1(i5,j4+1)=r1.faces(r1.shift(j4),:); % ref nodes in cur node index
   end
   % find a coherent renumbering if several faces are used
   if size(n1,2)>2
    i5=false(size(n1,1),1); 
    for j5=1:length(i5); 
     r6=diff(n1(j5,n1(j5,:)~=0),[],2);
     if ~isempty(r6)&&any(r6); i5(j5)=true; end
    end
    if any(i5); %diff(n1(:,2:end),[],2)) % coherence to resolve
     % rotate first face until common feature is in agreement with all faces
     n2=n1(all(n1(:,1:end),2),1); % common nodes on elt cfg
     n3=n1(ismember(n1(:,2),n1(:,3)),2);
     for j5=4:size(n1,2); n3=n3(ismember(n3,n1(:,j5))); end
     n3(n3==0)=[]; % common nodes on ref
     if ~isempty(n3) % what if two disjoint faces are used ?
      % rotate face1 until common nodes are the same at the same topo place
      r7=r1.faces(r1.shift(1),:); in1=ismember(r7,n3); ic=0;
      r8=r5(r1.shift(1),:); r80=r8;
      while ~all(ismember(r8(in1),n2))&&ic<length(r7)
       ic=ic+1; r8=circshift(r80,[0 ic]); %7.14 compat xxx
      end
     n1(r8,2)=r7;
     end
     % does not seem necessary to do further checks
     
    else % complete on col 2 to help
     for j5=find(n1(:,2)==0&any(n1(:,2:end),2))';
      n1(j5,2)=n1(j5,1+find(n1(j5,2:end),1,'first')); 
     end
    end
   end
   if any(n1(:,2)==0)
    % now find other nodes by getting another face with the less common nodes
    % with the reference face(s). ok for all elts where another face always
    % contains the remaining points
    r6=unique(r5(r1.shift(1),:));
    i3=sum(~ismember(r5,r6)&r5~=0,2); 
    i3=sort(find(i3==max(i3)),'ascend'); %,1,'first'); 
    r6b=unique(r1.faces(r1.shift(1),:));
    i3b=sum(~ismember(r1.faces,r6b)&r1.faces~=0,2); 
    i3b=sort(find(i3b==max(i3b)),'ascend'); %,1,'first'); 
    while any(n1(:,2)==0)&&~isempty(i3)
     r7=r5(i3(1),:);
     [i4,i5]=ismember(r7,r6); i4=~i4&r7~=0; %&~ismember(r7,find(n1(:,2))); %i4=intersect(find(~i4&r7~=0),find(n1(:,2)'==0));
     %%n1(r5(i3,~i4),2)=circshift(r1.faces(i3(1),i4),ic,2);
     %n1(r5(i3(1),i4),2)=circshift(r1.faces(i3(1),i4),[0 ic]);%7.14 compat
     n2=circshift(r1.faces(i3b(1),:),[0 ic]);
     n1(r5(i3(1),i4),2)=n2(1,i4); %circshift(r1.faces(i3b(1),i4),[0 ic]);%7.14 compat
     %%n1(r5(i3(1),i4),2)=circshift(r1.faces(iface(i3(1)),i4),[0 ic]);%7.14 compat
     i3(1)=[]; i3b(1)=[];
    end
   end
   n1=n1(:,1:2);
   % add refined cell new nodes that will not be actually renumbered in the ref. topo
   r8=[max(cellfun(@(x)max(max(x)),r1.Elt(:,2))) max(max(r1.faces))];
   n1=[n1;[[r8(2)+1:r8(1)]' [r8(2)+1:r8(1)]']]; %#ok<AGROW>
   % prepare renumbering: applied to edge/face/volume node selection and 
   % on new elements connectivity
   NN=sparse(n1(:,2),1,n1(:,1)); % cleaned strategy with multi faces
   st=sprintf('%i,',r3(j3,:)); %st(end)='';
   for j2=1:size(r1.Elt,1)
    eval(sprintf('r1.EltNum{%sj2}=full(reshape(NN(r1.Elt{j2,2}),size(r1.Elt{j2,2},1),[]));',st));
   end
   % how eltnum was handled before support for multiple face reference selection
   %    NN=sparse([r1.faces(r1.shift,:) r1.faces(i3,~i4) r8(2)+1:r8(1)],1,...
   %     [r5(r1.shift,:) r5(i3,~i4) r8(2)+1:r8(1)]);
   %r1.EltNum{r3(j3),2}=full(NN(r1.Elt{2})); % new elts for rotated face
  end
  for j2=1:size(r2,1) % loop on edge/face,volume, 
   %% fill eCoef,fCoef,vCoef
   if isfield(r1,r2{j2,2}) % edge/face,volume given
    stc=sprintf('%sCoef',r2{j2,2}(1));
    if isempty(RO.(stc){jGroup}) % not filled yet (if r3loop)
     r4=r1.(r2{j2,2}){2}; % base coef
     if length(r4)>1 % non uniform coef need replication for all
      RO.(stc){jGroup}=repmat(r4,length(cEGI),1);
     else; RO.(stc){jGroup}=r4; % will be uniform
     end
    end
    %% fill Edges,Faces,Volumes, will provide base elt nodes for newCoor
    r4=r1.(r2{j2,2}){1}(:,2:end)'; % nodes used for newCoor interp
    if r3(j3) % renum due to face orientation
     r4=full(NN(r4));
     if isempty( RO.(r2{j2,1}){jGroup,1}) % preallocate here (r3 loop:several passes)
      RO.(r2{j2,1}){jGroup,1}=zeros(size(r4,2)*length(cEGI),size(r4,1));
      RO.(r2{j2,1}){jGroup,2}=r1.(r2{j2,2}){1}(:,1);
     end
     %% fill for elements of correct FaceId selection
     i3=cEGI(ismember(cEGI,EEid(RO.set(ismember(abs(RO.set(:,2:end)),r3(j3,:),'rows'),1)+1))); % good cEGI
     i4=reshape(1:(length(cEGI)*size(r4,2)),[],length(cEGI))';
     i4=reshape(i4(ismember(cEGI,i3),:)',[],1); % clean index in Edges/Faces/Volumes
     RO.(r2{j2,1}){jGroup,1}(i4,:)=reshape(FEelt(i3,r4)',...
      size(r1.(r2{j2,2}){1},2)-1,[])';
     
    else % no renum needed direct fill
     RO.(r2{j2,1})(jGroup,1:2)={reshape(FEelt(cEGI,r4)',...
      size(r1.(r2{j2,2}){1},2)-1,[])',r1.(r2{j2,2}){1}(:,1)};
    end % if r3
   end % if edge/face/volume
  end % loop on edge/face,volume, 
 end % loop on FaceId
 if isfield(r1,'EltNum'); RO.(ElemP)=r1; end % new renum elts to store
end % group wise initialization
NNode=sparse(FEnode(:,1),1,1:size(FEnode,1));
%% Add Edges/Faces/Volumes
for j2=1:size(r2,1)
 if any(~cellfun(@isempty,RO.(r2{j2,1}))) % something to be done
  stm=sprintf('MPC%s',r2{j2,2}); stmc=sprintf('%scell',stm);
  stc=sprintf('%sCoef',r2{j2,2}(1));
  if isfield(RO,stm) % also generate MPC
   RO.(stm)=sort(RO.(stm),2);
   [n1,i2,i3,RO.(stmc)]=newCoor(FEnode,RO.(r2{j2,1}),RO.(stc),NNode);
  else;[n1,i2,i3]=newCoor(FEnode,RO.(r2{j2,1}),RO.(stc),NNode); % no mpc need
  end
  i2=0; % add nodes
  if RO.KnownNew; FEnode=[FEnode;n1(~ismember(n1(:,1),FEnode(:,1)),:)]; %#ok<AGROW>
  else
   [FEnode,i4]=feutil('addnode',FEnode,n1); n1=FEnode(i4,:); % keep initial n1 order
   if isfield(RO,stmc); RO.(stmc){1}=n1(:,1); end % also renum data for MPC
  end
  for jGroup=1:nGroup % prepare node ordering for final topology in Points
   if isempty(RO.(r2{j2,1}){jGroup,2}); continue; end
   i2=i2(end)+(1:size(RO.Points{jGroup},1)*length(RO.(r2{j2,1}){jGroup,2}));
   RO.Points{jGroup}(:,RO.(r2{j2,1}){jGroup,2})=reshape(n1(i3(i2),1), ...
    length(RO.(r2{j2,1}){jGroup,2}),[])';
  end
 end
end
if issparse(FEnode); FEnode=full(FEnode); end
%% Create elements
out=struct('Node',FEnode,'Elt',[]); EEIdx=[];
for jGroup=1:nGroup
 [ElemF,i1,ElemP]= getegroup(FEelt(EGroup(jGroup),:),jGroup);
 cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
 if ~isempty(RO.subEGI); [cEGI,isub]=intersect(cEGI,RO.subEGI); end % restrain cEGI
 if isempty(cEGI); continue;
 elseif ~isfield(RO,ElemP); % Keep non refined elements
  out=feutil('addelt',out,FEelt(EGroup(jGroup):EGroup(jGroup+1)-1,:));
  continue;
 end
 r1=RO.(ElemP); 
 
 for j2=1:size(r1.Elt,1) % loop on new elemf
  i1=RO.Points{jGroup}'; % get 
  if isfield(r1,'shift')&&size(RO.set,2)>1
   % unsymm refine: loop on FaceId to add renumbered connectivity
   i2=zeros(length(cEGI)*size(r1.Elt{j2,2},1),size(r1.Elt{j2,2},2)); % pre allocate output
   i3=ismember(RO.set(:,1),eltid(cEGI)); % find corresponding EltId
   r3=unique(abs(RO.set(i3,2:end)),'rows'); % faces to treat
   for j3=1:size(r3,1)
    i4=find(abs(RO.set(i3,2))==r3(j3)); % current eltid
    st=sprintf('%i,',r3(j3,:)); %st(end)='';
    eval(sprintf('r4=r1.EltNum{%sj2};',st)); % renumbered topology
    i4=ismember(cEGI,EEid(RO.set(ismember(abs(RO.set(:,2:end)),r3(j3,:),'rows'),1)+1)); % corresponding cEGI
    % index in i2 (newElt x nodes) for current cEGI only
    in1=reshape(1:(length(cEGI)*size(r4,1)),[],length(cEGI))';
    in1=reshape(in1(i4,:)',[],1);
    i2(in1,:)=reshape(i1(r4',i4),size(r4,2),[])';
   end
   i1=i2; clear i2
   
  else % direct generation if no shift
   i1=reshape(i1(r1.Elt{j2,2}',:),size(r1.Elt{j2,2},2),[])'; % i1 NewElts x nodes
  end
  
  % assign ElemF, initial mat/pro and output
  i2=reshape(repmat(RO.pro{jGroup}',size(r1.Elt{j2,2},1),1), ...
   size(RO.pro{jGroup},2),[])';
  if ~RO.keepEP&&size(r1.Elt,1)==1 % maintain ElemF over ElemP if transformation is consistent
   stE=r1.Elt{1};
   if isequal(ElemP,stE); stE=ElemF; % consistent transformation
   elseif isempty(RO.Silent)&&~Silent
    sdtw('_nb',['Group %i elements of type %s are transformed into %s\n'...
     'Use option -keepEP to avoid this warning'],jGroup,ElemF,stE);
   end
  else; stE=r1.Elt{j2,1};
  end
  if ~isempty(i2); i3=feval(stE,'prop'); i1(:,i3)=i2; end % EltId is kept
  [out]=feutil('addelt',out,stE,i1); %[i1 i2]);
 end % loop on new elemf
 
end
%% Replace elements in base model if not all elements were transformed
out1=[]; out2=[];
if ~isempty(RO.subEGI)&&RO.replace&&~isempty(out.Elt)
 eltid=feutil('eltid;',FEelt);
 [out.Elt,EEIdx]=feutil('addelt-newId;',FEelt,out.Elt);% FEelt initial, Out.Elt new
 [i2,out.Elt]=feutil('eltidfix;',out); % Refined use different eltid, fix not needed)
 if ~any(eltid); eltid=i2(1:size(FEelt,1)); end
 eltid=setdiff(i2,eltid); r5=[]; % new eltId
 if RO.noEOri; data=[]; else; data=stack_get(model,'info','EltOrient','get'); end %#ok<NODEF>
 if RO.keepSets % store meta-set for recast
  r5=feutil('addseteltid-append-NoNodes-get',stack_rm(model,'set','_gsel'),'_gsel');
 end
 if ~isempty(data); 
  model.Node=out.Node;model.Elt=out.Elt;
  RB=struct('Type','eltorient','Out','EltOrient', ...
     'EltSel',['EltId~=' sprintf(' %i',RO.set(:,1))], ...%sprintf('EltInd>%i',size(FEelt,1)), ...
     'KeepOld',1);%);
  [model.Elt,RB.RemElt]=feutil('removeelt EltId',out,RO.set(:,1)); 
  out=fe_shapeoptim('interp',model,RB); out1=RB.RemElt;
 elseif ~isempty(RO.set)
  [out.Elt,out1]=feutil('removeelt EltId',out,RO.set(:,1));
 elseif ~isempty(RO.subEGI)
  [out.Elt,out1]=feutil('RemoveElt EltInd',out,RO.subEGI); % initial FEelt order
 end
 if nargout>2; % also provide cEGI of new elts
  out2=feutil('findelt eltid',out,eltid);
 end
 if RO.replace==2; % Actually keep initial model properties
  model.Node=out.Node;model.Elt=out.Elt;out=model; 
  if isfield(out,'DOF');out=rmfield(out,'DOF');end
 end
 % replace eltsets
 if ~isempty(r5); 
  out=feutil('EltSetReplace',stack_set(out,stack_get(ModelStack,'set')),r5,EEIdx); 
 end
 %% Handle MPC if needed, consider DOF as 1 2 3
 for j2=1:size(r2,1)
  stm=sprintf('MPC%s',r2{j2,2}); stmc=sprintf('MPC%scell',r2{j2,2});
  if isfield(RO,stm)&&isfield(RO,stmc)
   r4=RO.(stmc); % {new,master,[coef]}
   if size(r4{2},2)==2 % Possibly degenerate edge
    i5=r4{2}(:,1)==r4{2}(:,2);r4{1}(i5,:)=[];r4{2}(i5,:)=[];r4{3}(i5,:)=[];
    RO.(stmc)=r4;
   end
   if length(RO.(stm))<1; 
    i1=[1:length(RO.(stmc){1})]'; % all wanted xxx should be free faces/edges
   else; i1=find(ismember(RO.(stmc){2},RO.(stm),'rows')); % sub set wanted
   end
   dof=feutil('getdof',[RO.(stmc){1}(i1);unique(RO.(stmc){2}(i1,:))],[1:3]'/100);
   nd=feval(fe_mknl('@getPosFromNd'),[],dof);
   r3=sum(cellfun(@(x)size(x,2),RO.(stmc)(1:2)));
   II=zeros(length(i1)*3*r3,1);  JJ=0*II;  KK=0*JJ; 
   ii=zeros(1);ij=zeros(1);ik=zeros(1);
   for j3=1:3 % loop on DOF ext
    ii=of_time([-1 ii],II,...
     reshape(length(i1)*(j3-1)+[1:length(i1)]'*ones(1,r3),[],1));
    ij=of_time([-1 ij],JJ, feval(nd.getPosFcn,nd,...
     [RO.(stmc){1}(i1);reshape(RO.(stmc){2}(i1,:),[],1)]+j3/100) );
    ik=of_time([-1 ik],KK,[ones(length(i1),1);-reshape(RO.(stmc){3}(i1,:),[],1)]);
   end
   out=fe_case(out,'stack_setnew','mpc',stm,struct('DOF',dof,...
    'c',sparse(II,JJ,KK,3*length(i1),length(dof)),...
    'slave',fe_c(dof,feutil('getdof',[1:3]'/100,RO.(stmc){1}(i1)),'ind',1)));
  end
 end
elseif RO.replace&&isempty(out.Elt)
 out=struct('Node',FEnode,'Elt',FEelt);
 if RO.keepSets; out=stack_set(out,stack_get(ModelStack,'set')); end
end

%% #RefineEnd -2
else;error('Refine%s unknown',CAM);
end

%% #Lin2 : linear to/from transformations

%% #quad2tria -----------------------------------------------------------2
elseif comstr(Cam,'quad2tria');  [CAM,Cam] = comstr(CAM,7);

 model=[];FEel0=varargin{carg};carg=carg+1;
 if isstruct(FEel0);model=FEel0;FEel0=FEel0.Elt;end
 [EGroup,nGroup]=getegroup(FEel0);

  i2 =[]; % i2 first group nodes, i3 second group nodes
  for jGroup = 1:nGroup %loop on element groups
   ElemF= feutil('getelemf',FEel0(EGroup(jGroup),:),jGroup);
   if strcmp(ElemF,'quad4')
     i3 = FEel0(EGroup(jGroup)+1:EGroup(jGroup+1)-1,1:6);
     i3 = reshape([i3(:,[1 2 3 5 6]) i3(:,[1 3 4 5 6])]',5,2*size(i3,1))';
     elt=[Inf abs('tria3') zeros(1,size(FEel0,2)-6)
      i3 zeros(size(i3,1),size(FEel0,2)-5)];
     elt(elt(:,1)==elt(:,2)|elt(:,2)==elt(:,3)|elt(:,1)==elt(:,3),:)=[]; % degen
     FEel0 = [FEel0(1:EGroup(jGroup)-1,:);elt
      FEel0(EGroup(jGroup+1):size(FEel0,1),:)];
      EGroup=getegroup(FEel0);
   end
  end % of jGroup loop
  if isstruct(model);model.Elt=FEel0;out=model; else; out=FEel0; end


%% #hexa2tetra -----------------------------------------------------------2
elseif comstr(Cam,'hexa2tetra')
  [CAM,Cam]=comstr(CAM,11);
  [CAM,Cam,RunOpt.KnownNew]=comstr('knownnew',[-25 3],CAM,Cam);
  
  [carg,FEnode,FEel0,elt,ModelStack]=get_nodeelt(varargin,carg,ModelStack);
  [EGroup,nGroup]=getegroup(FEel0);
  i2 =[]; % i2 first group nodes, i3 second group nodes

  for jGroup = 1:nGroup %loop on element groups
   ElemF= feutil('getelemf',FEel0(EGroup(jGroup),:),jGroup);
   if strcmp(ElemF,'hexa8')
     i4=[]; cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
     [FEnode,i2]=feutil('NewMid',FEnode,FEel0(cEGI,:), ...
         {[1 2 3 4;5 6 7 8;1 2 6 5;2 3 7 6;3 4 8 7;1 4 8 5]',(1:8)'},RunOpt);

     
     i3=[FEel0(cEGI,1:8) i2 FEel0(cEGI,9:10)]';
     i3=i3([1 2 9 15 16 17;2 3 9 15 16 17;3 4 9 15 16 17;4 1 9 15 16 17;
         6 5 10 15 16 17;7 6 10 15 16 17;8 7 10 15 16 17;5 8 10 15 16 17;
         2 1 11 15 16 17;6 2 11 15 16 17;5 6 11 15 16 17;1 5 11 15 16 17;
         3 2 12 15 16 17;7 3 12 15 16 17;6 7 12 15 16 17;2 6 12 15 16 17;
         4 3 13 15 16 17;8 4 13 15 16 17;7 8 13 15 16 17;3 7 13 15 16 17;
         1 4 14 15 16 17;4 8 14 15 16 17;8 5 14 15 16 17;5 1 14 15 16 17]',:);
     i3=reshape(i3,6,[])';
     FEel0 = [FEel0(1:EGroup(jGroup)-1,:);
      Inf abs('tetra4') zeros(1,size(FEel0,2)-7)
      i3 zeros(size(i3,1),size(FEel0,2)-6)
      FEel0(EGroup(jGroup+1):size(FEel0,1),:)];
      EGroup=getegroup(FEel0);
   end
  end % of jGroup loop

model.Node=FEnode; model.Elt=FEel0; out=model;

%% #hexa2pyra -----------------------------------------------------------2
elseif comstr(Cam,'hexa2pyra')
  [CAM,Cam]=comstr(CAM,11);
  [CAM,Cam,RunOpt.KnownNew]=comstr('knownnew',[-25 3],CAM,Cam);
  
  [carg,FEnode,FEel0,elt,ModelStack]=get_nodeelt(varargin,carg,ModelStack);
  [EGroup,nGroup]=getegroup(FEel0);
  i2 =[]; % i2 first group nodes, i3 second group nodes

  for jGroup = 1:nGroup %loop on element groups
   ElemF= feutil('getelemf',FEel0(EGroup(jGroup),:),jGroup);
   if strcmp(ElemF,'hexa8')
     i4=[]; cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
     [FEnode,i2]=feutil('NewMid',FEnode,FEel0(cEGI,:), ...
         {[],(1:8)'},RunOpt);
     i3=hexa8('faces'); i3=i3(:,[1 4 3 2]);
     i3=[FEel0(cEGI,i3(1,:)) i2 FEel0(cEGI,9:10)
         FEel0(cEGI,i3(2,:)) i2 FEel0(cEGI,9:10)
         FEel0(cEGI,i3(3,:)) i2 FEel0(cEGI,9:10)
         FEel0(cEGI,i3(4,:)) i2 FEel0(cEGI,9:10)
         FEel0(cEGI,i3(5,:)) i2 FEel0(cEGI,9:10)
         FEel0(cEGI,i3(6,:)) i2 FEel0(cEGI,9:10)];
     FEel0=[FEel0(1:EGroup(jGroup)-1,:);
      Inf abs('pyra5') zeros(1,size(FEel0,2)-6)
      i3 zeros(size(i3,1),size(FEel0,2)-size(i3,2))
      FEel0(EGroup(jGroup+1):size(FEel0,1),:)];
      EGroup=getegroup(FEel0);
   end
  end % of jGroup loop
model.Node=FEnode; model.Elt=FEel0; out=model;

elseif comstr(Cam,'hexa2penta');  [CAM,Cam] = comstr(CAM,7);
%% #hexa2penta -----------------------------------------------------------2
  [CAM,Cam,RunOpt.KnownNew]=comstr('knownnew',[-25 3],CAM,Cam);
  [CAM,Cam,RunOpt.face]=comstr('face',[-25 31],CAM,Cam);
  [carg,FEnode,FEel0,elt,ModelStack]=get_nodeelt(varargin,carg,ModelStack);
  [EGroup,nGroup]=getegroup(FEel0);

  i2 =[]; % i2 first group nodes, i3 second group nodes
  i4=[];

  for jGroup = 1:nGroup %loop on element groups
   ElemF= feutil('getelemf',FEel0(EGroup(jGroup),:),jGroup);
   if strcmp(ElemF,'hexa8')
     i4=[]; cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
     if RunOpt.face==2; FEel0(cEGI,1:8)=FEel0(cEGI,[1 2 6 5 4 3 7 8]);end
     [FEnode,i2]=feutil('NewMid',FEnode,FEel0(cEGI,:), ...
      {[1 2 3 4;5 6 7 8]'},RunOpt);
     for jElt=1:length(cEGI)
       i3=[FEel0(cEGI(jElt),1:8) i2(jElt,:)];
       i5=[i3(:,[1 2 end-1 5 6 end]);
           i3(:,[2 3 end-1 6 7 end]);
           i3(:,[3 4 end-1 7 8 end]);
           i3(:,[4 1 end-1 8 5 end])];
       i5(:,7:8)=ones(4,1)*FEel0(cEGI(jElt),9:10);
       i4=[i4;i5];
     end

     FEel0 = [FEel0(1:EGroup(jGroup)-1,:);
      Inf abs('penta6') zeros(1,size(FEel0,2)-7)
      i4 zeros(size(i4,1),size(FEel0,2)-8)
      FEel0(EGroup(jGroup+1):size(FEel0,1),:)];
      EGroup=getegroup(FEel0);
   end
  end % of jGroup loop
  if ~any(Cam==';');
   sdtw('_nb','Be carefull using hexa2penta, always check mesh');
  end
  model.Node=FEnode; model.Elt=FEel0; out=model;

%% #penta2tetra : (contributed by A. Sternschuss) ------------------------2
elseif comstr(Cam,'penta2tetra'); [CAM,Cam]=comstr(CAM,12);
  [CAM,Cam,RunOpt.KnownNew]=comstr('knownnew',[-25 3],CAM,Cam);
  [carg,FEnode,FEel0,elt,ModelStack]=get_nodeelt(varargin,carg,ModelStack);
  [EGroup,nGroup]=getegroup(FEel0);

  i2 =[]; % i2 first group nodes, i3 second group nodes
  

  for jGroup = 1:nGroup %loop on element groups
   ElemF= feutil('getelemf',FEel0(EGroup(jGroup),:),jGroup);
   if strcmp(ElemF,'penta6');
     i4=[]; cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
     [FEnode,i2]=feutil('NewMid',FEnode,FEel0(cEGI,:), ...
        {[1 2 5 4;2 3 6 5;3 1 4 6]'},RunOpt);
    
     i3=[FEel0(cEGI,1:6) i2 FEel0(cEGI,7:8)]';
     i3=i3([1 2 3 7 10 11;2 8 3 7 10 11;3 1 7 9 10 11;3 7 8 9 10 11;
         4 6 5 7 10 11;5 6 8 7 10 11;6 7 4 9 10 11;6 8 7 9 10 11;
         1 4 7 9 10 11;2 7 5 8 10 11;3 8 6 9 10 11]',:);
     i3=reshape(i3,6,[])';
     FEel0 = [FEel0(1:EGroup(jGroup)-1,:);
      Inf abs('tetra4') zeros(1,size(FEel0,2)-7)
      i3 zeros(size(i3,1),size(FEel0,2)-6)
      FEel0(EGroup(jGroup+1):size(FEel0,1),:)];
      EGroup=getegroup(FEel0);
   end
  end % of jGroup loop
  model.Node=FEnode; model.Elt=FEel0; out=model;

%% #Lin2Quad #Quad2Lin -------------------------------------------------------
% order transformation (works to go from linear to quadratic elements)
% lin2quadCyl for cylindrical coordinates
elseif comstr(Cam,'lin2quad')||comstr(Cam,'quad2lin')

 [CAM,Cam,RunOpt.sel]=comstr('-sel',[-25 4],CAM,Cam);
 if ~isempty(RunOpt.sel) % done on selection
  [CAM,Cam,RunOpt.mpc]=comstr('-mpc',[-25 3],CAM,Cam);
  model=varargin{carg}; carg=carg+1;
  mo1=model;
  [model.Elt,mo1.Elt]=feutil('removeelt',model,RunOpt.sel);

  n1=feutil('getnodegroupall',mo1);
  mo1=feutil(CAM,mo1);
  n2=feutil('getnodegroupall',model);
  
  model=feutil('addtestcombine;',model,...
   struct('Node',feutil('getnodegroupall',mo1),'Elt',mo1.Elt));
  
  if RunOpt.mpc
   mo2=mo1;
   mo2.Elt=feutil(sprintf('selelt selface & withnode{nodeid %s}',num2str(n2(:,1)')),mo2);
   [mo2,r1]=StraightenEdge(mo2,mo2.Elt); % r1: edge
   
   % Now generate the MPC: 1*q_mid - 0.5*q_1 - 0.5*q_2 = 0
   dof=feutil('getdof',[.01;.02;.03],unique(r1(:)));
   [i1,i2]=ismember(feutil('getdof',[.01;.02;.03],...
    reshape(circshift(r1',1),1,[])'),dof);
   II=reshape(ones(3,1)*[1:3*size(r1,1)],[],1); JJ=i2;
   KK=reshape([1;-.5;-.5]*ones(1,3*size(r1,1)),[],1);
   c=sparse(II,JJ,KK,3*size(r1,1),length(dof));
   [c,r1]=feutil('fixMpcMaster',c);
   c=struct('DOF',dof,'c',c,'slave',r1);
   
   model=fe_case(model,'mpc','QuadTrans',c);
  end
  
  out=model;
  return
 end
 
 RunOpt.Stack={'tetra4','tetra10',[1 2;2 3;3 1;1 4;2 4;3 4];
  'hexa8','hexa20',[1 2;2 3;3 4;4 1;1 5;2 6;3 7;4 8;5 6;6 7;7 8;8 5];
  'tria3','tria6',[1 2;2 3;3 1];
  'beam1','beam3',[1 2];
  'quad4','quadb',[1 2;2 3;3 4;4 1];
  'penta6','penta15',[1 2;2 3;3 1;1 4;2 5;3 6;4 5;5 6;6 4];
  'pyra5','pyra13',[1 2;2 3;3 4;1 4;1 5;2 5;3 5;4 5];
  };

 %if 1==2
 % i3=[1 2;2 3;3 4;4 1;1 5;2 6;3 7;4 8;5 6;6 7;7 8;8 5];
 % i4=max(i3(:));
 % r1=sparse(i3,[1:size(i3,1)]'*[1 1],.5);st={'%g '};st=st(ones(i4,1));
 % fprintf(['      ' st{:} ';\n'],full(r1)); %#ok<ACCUM>
 %end
 [epsl,CAM,Cam]=test_epsl(epsl,CAM,Cam);
 [CAM,Cam,RunOpt.KnownNew]=comstr('knownnew',[-25 3],CAM,Cam);
 [CAM,Cam,RunOpt.Cyl]=comstr('cyl',[-25 3],CAM,Cam);
 if comstr(Cam,'quad'); RunOpt.ToLin=1;else;RunOpt.ToLin=0;end
 
 [carg,FEnode,FEelt,FEel0,ModelStack]=get_nodeelt(varargin,carg,ModelStack);
 if carg<=nargin&&ischar(varargin{carg}) % transform only one type of element
  ind=find(strcmpi(varargin{carg},RunOpt.Stack(:,1)));carg=carg+1;
  if isempty(ind);ind=find(strcmpi(varargin{carg},RunOpt.Stack(:,2)));end
  RunOpt.Stack=RunOpt.Stack(ind,:);
 end
 [EGroup,nGroup]=getegroup(FEelt);
 if carg<=nargin&&isstruct(varargin{carg});
     RunOpt=sdth.sfield('addmissing',RunOpt,varargin{carg});carg=carg+1;
 elseif carg>nargin||any(varargin{carg}>nGroup);RunOpt.nGroup=1:nGroup;
 else; RunOpt.nGroup=varargin{carg};carg=carg+1; % selected group numbers
 end
 if ~isfield(RunOpt,'nGroup');RunOpt.nGroup=1:nGroup;end
 NNode=sparse(FEnode(:,1),1,1:size(FEnode,1));
 RunOpt.I=cell(max(RunOpt.nGroup),5);
 for jGroup = RunOpt.nGroup %loop on element groups
   [ElemF,i1,ElemP]=feutil('getelemf',FEelt(EGroup(jGroup),:),jGroup);
   if RunOpt.ToLin
    ind=find(strcmpi(ElemP,RunOpt.Stack(:,2)));
    if ~isempty(ind) % there is a match
      cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;
      r1=RunOpt.Stack{ind,3}; inode=feval(ElemP,'node');  
      RunOpt.I(jGroup,1:2)={RunOpt.Stack{ind,1} ...
       [feval(RunOpt.Stack{ind,1},'node') feval(ElemP,'prop')]};
      i2=max(feval(ElemP,'prop')); % add zeros in FEelt if no prop
      if size(FEelt,2)>i2&&...
        (strncmpi(ElemP,'tria',4)||strncmpi(ElemP,'quad',4))&&any(FEelt(cEGI,i2+1))
       RunOpt.I{jGroup,2}(1,end+1)=i2+1; i2=i2+1; % keep shell theta
      end
      if size(FEelt,2)<i2; FEelt(end,i2)=0; end
    end
   else
    ind=find(strcmpi(ElemP,RunOpt.Stack(:,1)));
    if ismember(ElemF,{'celas','rigid','cbush'}); continue;end
    if ~isempty(ind) % there is a match
      cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;
      r1=RunOpt.Stack{ind,3}; inode=feval(ElemP,'node');  r2=3;
      if strncmpi(ElemP,'tria',4)||strncmpi(ElemP,'quad',4); r2=4; end % shell theta
      RunOpt.I(jGroup,1:5)={RunOpt.Stack{ind,2} [] [] length(inode) ...
         FEelt(cEGI,length(inode)+1:min(size(FEelt,2),length(inode)+r2))};
      i1=RunOpt.Stack{ind,3}';RunOpt.I{jGroup,3}=size(i1,2);% generate edges
      i2=FEelt(cEGI,inode)'; if any(i2==0); continue;end
      i2=reshape(full(NNode(i2)),length(inode),length(cEGI)); 
      i2=sort(reshape(i2(i1,:),2,length(cEGI)*size(i1,2)))';
      RunOpt.I{jGroup,2}=i2;
    end
   end
 end
 if ~RunOpt.ToLin;[FEnode,RunOpt]=feutil('NewMid',FEnode,RunOpt,{[1;2]});end
  
 for jGroup = RunOpt.nGroup %loop on element groups
  if ~isempty(RunOpt.I{jGroup,1});
   cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;
   ElemF=RunOpt.I{jGroup,1};
   FEelt(EGroup(jGroup),1:length(ElemF)+2)=[Inf abs(ElemF) 0];
   if RunOpt.ToLin
    r1=zeros(length(cEGI),size(FEelt,2));
    i1=RunOpt.I{jGroup,2};r1(:,1:length(i1))=FEelt(cEGI,i1);
    FEelt(cEGI,:)=r1;
   else;
    i4=[RunOpt.I{jGroup,3} RunOpt.I{jGroup,5}];
    FEelt(cEGI,RunOpt.I{jGroup,4}+[1:size(i4,2)])=i4;
   end
  end % no match
 end % jGroup
 if nargout==1; model.Node=FEnode; model.Elt=FEelt;out=model;
  if ~isempty(strfind(Cam,'-optim'))
   if ~sp_util('issdt'); sdtw('_nb','-optim option requires feutilb');return;end
   if RunOpt.ToLin
    n1=model.Node; model.Node=feutilb('GetUsedNodes',model);
    r1=n1(:,1); r1(:,2)=0; NNode=sparse(n1(:,1),1,1:size(n1,1));
    r1(NNode(model.Node(:,1)),2)=model.Node(:,1); clear n1 NNode
    nind=sparse(r1(:,1),1,r1(:,2)); RA=struct('Silent',';','clean',1,'keepCase',1);
    model=feval(feutilb('@renumber_stack'),model,nind,model,[],RA);
   end
  end
  out=model;
 else
  out=FEnode; out1=FEelt;out2=[];
  if isfield(RunOpt,'EdgeI'); out2=RunOpt.EdgeI;end
 end
  
%% #quad42quadb ----------------------------------------------------------2
elseif comstr(Cam,'quad42quadb');  [CAM,Cam] = comstr(CAM,7);

  [carg,FEnode,FEelt,FEel0,ModelStack]=get_nodeelt(varargin,carg,ModelStack);
  [out,out1]=feutil('lin2quad',FEnode,FEelt,[],'quad4');

%% #quad42quad9 : q42q9 --------------------------------------------------2
elseif comstr(Cam,'quad42quad9');  [CAM,Cam] = comstr(CAM,7);

  [carg,FEnode,FEelt,FEel0,ModelStack]=get_nodeelt(varargin,carg,ModelStack);
  [EGroup,nGroup]=getegroup(FEelt);

  r1=[.5 .5 0 0;0 .5 .5 0;0 0 .5 .5;.5 0 0 .5;.25 .25 .25 .25];
  i2 =[]; % i2 first group nodes, i3 second group nodes
  for jGroup = 1:nGroup %loop on element groups
   ElemF= feutil('getelemf',FEelt(EGroup(jGroup),:),jGroup);
   if strcmp(ElemF,'quad4')
     FEelt(EGroup(jGroup),1:7)=[Inf abs('quad9') 0];
     for j1=EGroup(jGroup)+1:EGroup(jGroup+1)-1
       NNode=sparse(FEnode(:,1),1,1:length(FEnode(:,1)));
       r2 = r1*FEnode(NNode(FEelt(j1,1:4)),5:7);
       [FEnode,i2]=feutil('AddNode',FEnode,r2);
       r2=[FEelt(j1,1:4) FEnode(i2(:),1)' FEelt(j1,5:max(find(FEelt(j1,:))))];
       FEelt(j1,1:length(r2))=r2;
     end
   end
  end % of jGroup loop
  out=FEnode;out1=FEelt;

%% #quad42q5p ----------------------------------------------------2
elseif comstr(Cam,'quad42q5p');  [CAM,Cam] = comstr(CAM,7);

  [carg,FEnode,FEelt,FEel0,ModelStack]=get_nodeelt(varargin,carg,ModelStack);
  [EGroup,nGroup]=getegroup(FEelt);

  r1=[.25 .25 .25 .25];
  i2 =[]; % i2 first group nodes, i3 second group nodes
  for jGroup = 1:nGroup %loop on element groups
   ElemF= feutil('getelemf',FEelt(EGroup(jGroup),:),jGroup);
   if strcmp(ElemF,'quad4')
     FEelt(EGroup(jGroup),1:7)=[Inf abs('q5p') 0 0 0];
     for j1=EGroup(jGroup)+1:EGroup(jGroup+1)-1
       NNode=sparse(FEnode(:,1),1,1:length(FEnode(:,1)));
       r2 = r1*FEnode(NNode(FEelt(j1,1:4)),5:7);
       [FEnode,i2]=feutil('AddNode',FEnode,r2);
       r2=[FEelt(j1,1:4) FEnode(i2(:),1)' FEelt(j1,5:max(find(FEelt(j1,:))))];
       FEelt(j1,1:length(r2))=r2;
     end
   end
  end % of jGroup loop
  out=FEnode;out1=FEelt;

%% #Remove -------------------------------------------------------------------
elseif comstr(Cam,'remove'); [CAM,Cam] = comstr(CAM,7);

 if comstr(Cam,'elt') % #RemoveElt - - - - - - - - - - - - - - - - - - - - - -

  [i4,out1]=feutil(['Find elt' CAM(4:end)],varargin{carg:end});
  RunOpt=struct('ImplicitNode',1);
  [carg,FEnode,FEelt,FEel0,ModelStack]=get_nodeelt(varargin,carg,ModelStack,RunOpt);

  out=FEelt(setdiff(1:size(FEelt,1),i4),:);
  % remove empty groups
  i1=find(~isfinite(out(:,1)));i1=intersect(i1,i1-1);if isempty(i1);i1=[];end
  if ~isfinite(out(end,1));i1(end+1,1)=size(out,1);end
  if ~isempty(i1);out(i1,:)=[];end

 elseif comstr(Cam,'pro')||comstr(Cam,'mat') % #RemoveMat #RemovePro - - -
  % removing pl/il given by ID
  model=varargin{carg}; carg=carg+1;
  [CAM,Cam,RO.all]=comstr('-all',[-25 3],CAM,Cam);
  [i0,i0,RO.pro]=comstr('pro',[-25 3],CAM,Cam);
  [i0,i0,RO.mat]=comstr('mat',[-25 3],CAM,Cam);
  [i0,i0,RO.unused]=comstr('-unused',[-25 3],CAM,Cam);
  if carg<=nargin; i1=varargin{carg}; carg=carg+1; else; i1=[]; end
  r1={'mat','pro';'pl','il';2,5;1,2}; % pos in Rayleigh mat, pos in mpid
  [data,RO.iR]=stack_get(model,'info','Rayleigh','Get');
  mpid=feutil('mpid',model);
  for j1=1:size(r1,2) % loop mat/pro
   if RO.(r1{1,j1})
    if RO.all % get all pl or il IDs
     r2=fe_mat(sprintf('get%s',r1{2,j1}),model);
     if nnz(r2); RO.(r1{1,j1})=r2(:,1); else; RO.(r1{1,j1})=r2; end
    elseif RO.unused % remove unused
     r2=fe_mat(sprintf('get%s',r1{2,j1}),model);
     r3=unique(mpid(:,r1{4,j1})); r3=r3(r3~=0);
     r3=addImplicitMat(model,mpid,r1{4,j1},r3);
     if isempty(r2); r2=r3; else; r2=setdiff(r2(:,1),r3); end
     RO.(r1{1,j1})=r2;
    elseif isempty(i1) % ID specified in string
     [CAM,Cam,RO.(r1{1,j1})]=comstr(r1{1,j1},[-25 1],CAM,Cam);
    else; RO.(r1{1,j1})=i1; % ID specified as additional arg
    end
    if ~isempty(RO.(r1{1,j1}))
     if isfield(model,r1{2,j1})&&~isempty(model.(r1{2,j1})) %pl/il field
      model.(r1{2,j1})(ismember(model.(r1{2,j1})(:,1),RO.(r1{1,j1})),:)=[];
      if isempty(model.(r1{2,j1})); model.(r1{2,j1})=[]; end
     end
     [val,ival]=stack_get(model,r1{1,j1});
     if ~isempty(ival) % deal with stack entries
      % keep NLpro as some are not related to elements if unused
      if RO.unused; i2=find(cellfun(@(x) ~isfield(x,'NLdata'),val(:,3)));
      else; i2=[1:size(val,1)]'; 
      end
      model.Stack(ival(i2(cellfun(@(x) ...
       ismember(x.(r1{2,j1})(1),RO.(r1{1,j1})),val(i2,3)))),:)=[];
     end
     if ~isempty(data)&&size(data,2)>=r1{3,j1} % deal with info,Rayleigh entries
      i2=ismember(data(:,r1{3,j1}),RO.(r1{1,j1}));
      if any(i2);data(i2,:)=[]; end % mat/pro is referenced in Rayleigh -> remove
     end
     if isfield(model,'nmap')&&sp_util('issdt')
      r5=RO.(r1{1,j1}); r5=r5(:); r5(:,2)=-1;
      model.nmap=vhandle.nmap.renumber(model.nmap,r1{1,j1},r5);
     end
    end
   end%matpro loop
  end %for
  if ~isfield(model,'Stack')||isempty(model.Stack); model.Stack={}; end 
  if ~isempty(RO.iR); model=stack_set(model,'info','Rayleigh',data); end
  out=model;

 else;error('Remove%s unknown',CAM);
 end % subcommand selection - - - - - - - - - - - - - -

%% #Renumber -NoOri ----------------------------------------------------------
elseif comstr(Cam,'renumber'); [CAM,Cam] = comstr(CAM,9);
if ~isempty(strfind(Cam,';'))||Silent; RO.silent=';'; else; RO.silent=''; end
model=varargin{carg};carg=carg+1;
if ~isfield(model,'Elt') % detect if def
 if isfield(model,'DOF')
  model.Node=unique(fix(model.DOF),'stable');
 elseif isfield(model,'dof')
  model.Node=unique(fix(reshape(model.dof(:,1:min(size(model.dof,2),2)),...
   [],1)),'stable');
 else; error('input %s is not a model or a def curve',inputname(carg-1));
 end
 if carg<=nargin&&size(varargin{carg},2)==2||isscalar(varargin{carg})
 elseif isempty(RO.silent)
  sdtw('_nb',['Renumbering of a def structure may generate unwanted results\n'...
   'regarding coherence with its orginal model. Robust results can be obtained\n' ...
   'by providing an OrigNumbering input, or by performing a global shift'])
 end
end

if carg>nargin % renumber nodes from 1 to N
  i1=1:size(model.Node,1); NNode=sparse(model.Node(:,1)+1,1,i1);
elseif carg+1<=nargin% .Node already numbered only fix NNode
   i1=varargin{carg};carg=carg+1;
   NNode=varargin{carg};carg=carg+1;
else % provide new node numbers in i1
  i1=varargin{carg};carg=carg+1; 
  if size(i1,1)==1&&size(i1,2)==size(model.Node,1);i1=i1(:);end
  if size(i1,1)==size(model.Node,1)&&size(i1,2)~=2 % possible renumber
   NNode=sparse(model.Node(:,1)+1,1,i1); 
  elseif isfield(i1,'NNode'); 
   % struct('Renum',i1,'NNode',sparse(double(i1(:,1)),1,double(i1(:,2))))
   NNode=i1.NNode;i1=i1.Renum;
   i2=model.Node(:,1);[i3,i4]=ismember(i2,i1(:,1));
   i2(i3)=i1(i4(i3),2);
   i1=i2;
  elseif size(i1,2)==2
   %NNode=sparse(double(i1(:,1)),1,double(i1(:,2)));
   i2=model.Node(:,1);[i3,i4]=ismember(i2,i1(:,1));
   i2(i3)=i1(i4(i3),2);i1=i2;
   NNode=sparse(model.Node(:,1)+1,1,i1); 
  elseif isscalar(i1);
   i1=model.Node(:,1)+i1; NNode=sparse(model.Node(:,1)+1,1,i1);
  else; error('Not a valid case');
  end
end

% deal with local bases
% renumbering of reference to nodes would be good
if isfield(model,'bas')  
 if isa(model.bas,'v_handle'); model.bas=model.bas.GetData; end
        [node,model.bas]=basis('nodebas',model);
end

% renumber nodes in elements
if isfield(model,'Elt')
 [EGroup,nGroup]=getegroup(model.Elt);
 for jGroup = 1:nGroup %loop on element groups
  ElemF= feutil('getelemf',model.Elt(EGroup(jGroup),:),jGroup);
  i3=fe_super('node',ElemF); ind=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
  i2 = model.Elt(ind,i3); % nodes
  model.Elt(ind,i3) =reshape(full(NNode(i2+1)),size(i2,1),size(i2,2));
 end % of jGroup loop
end

if ~isempty(i1)
 model=stack_set(model,'info','OrigNumbering',int32([model.Node(:,1) i1(:)]));
 model.Node(:,1)=i1;
 if length(unique(i1))~=length(i1) % If some nodes are duplicated eliminate
   i2=[]; i2(flipud(i1))=length(i1):-1:1; 
   model.Node=model.Node(i2(i2~=0),:);
 end
end

try; 
    if isfield(model,'Stack')&&(size(model.Stack,1)>1||~isequal(model.Stack{1,2},'OrigNumbering'))     
     model=feutilb(['AddTestStack' RO.silent],model);
    end
end

% deal with deformations
if isfield(model,'DOF')
  i3=full(NNode(fix(model.DOF)+1));
  if any(i3(:)==0)&&isempty(RO.silent);warning('zero DOF in model.DOF renumbering');end
  model.DOF=rem(model.DOF,1)+i3;
end
if isfield(model,'dof')
 i3=min(size(model.dof,2),2);
 model.dof(:,1:i3)=rem(model.dof(:,1:i3),1)+full(NNode(fix(model.dof(:,1:i3))+1));
end
if isfield(model,'tdof')
 model.tdof(:,1)=rem(model.tdof(:,1),1)+full(NNode(fix(model.tdof(:,1))+1));
 if size(model.tdof,2)>1
  model.tdof(:,2)=sign(model.tdof(:,2)).*full(NNode(abs(model.tdof(:,2))+1));
 end
end
if isfield(model,'def')&&isfield(model.def,'DOF')
  model.def.DOF=rem(model.def.DOF,1)+full(NNode(fix(model.def.DOF)+1));
end
if isfield(model,'TR')&&isfield(model.TR,'DOF')
  model.TR.DOF=rem(model.TR.DOF,1)+full(NNode(fix(model.TR.DOF)+1));
  if isfield(model.TR,'adof')
   model.TR.adof=rem(model.TR.adof,1)+full(NNode(fix(model.TR.adof)+1));
  end
end
if isfield(model,'nmap')&&sp_util('issdt')
 model.nmap=vhandle.nmap.renumber(model.nmap,'Nodes',NNode);
end
if ~isempty(strfind(Cam,'-noori'));
    model=stack_rm(model,'info','OrigNumbering');
end
if ~isfield(model,'Elt')
 model=feutil('rmfield',model,{'Node'});
end
out=model;

%% #Rotate #RotateSel #RotateNode --------------------------------------------
% Rotate<Sel,Node> o x y z Ang nx ny nz  - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'rotate'); [CAM,Cam] = comstr(CAM,7);
 
 if comstr(Cam,'node'); [CAM,Cam] = comstr(CAM,5);
  RunOpt.Cam='node'; % sel or node
 elseif comstr(Cam,'sel'); [CAM,Cam] = comstr(CAM,4);
  RunOpt.Cam='sel'; % sel or node
 end
 
 model=[]; % Model can be assigned below
 [carg,FEnode,FEel0,r1,ModelStack]=get_nodeelt(varargin,carg,ModelStack);

if ~isempty(strfind(Cam,'o'))
    i1=strfind(Cam,'o');
    opt=comstr(CAM([1:i1-1 i1+1:length(CAM)]),[-1 0 0 0  0   1 0 0]);
    i5=opt(1:3); opt=[0 opt(4:7)];
else
    opt=comstr(CAM,[-1 0  0  1  0  0]);
    i1=find(FEnode(:,1)==opt(1));
    if isempty(i1) 
     fprintf('feutil Rotate : origin set to 0,0,0\n');
     i5=[0 0 0];
    else;i5=FEnode(i1(1),5:7);
    end
end
if opt(2)==0; error('Need a non-zero angle for Rotate');end
opt(2)=opt(2)/180*pi;
if norm(opt(3:5))==0; error('Need a non-zero direction for Rotate');end

if comstr(RunOpt.Cam,'node') % RotateNode
 
 if carg<=nargin;
  stnode=varargin{carg}; carg=carg+1;
  if ischar(stnode)
   NodeId=feutil(sprintf('FindNode %s',stnode),model);
  else;NodeId=stnode;
  end
 else; NodeId=FEnode(:,1);
 end
  
 bas=basis('rotate',[],...
  sprintf('r=%.15g;n=[%.15g %.15g %.15g]',opt(2)/pi*180,opt(3:5)),1);
 bas(4:6)=i5(:)'; % translation
 i1=ismember(FEnode(:,1),NodeId);
 node=FEnode(i1,:); 
 node(:,5)=node(:,5)-i5(1); node(:,6)=node(:,6)-i5(2); node(:,7)=node(:,7)-i5(3);
 node=basis('gnode',bas,node); % pour calculer les noeuds tournes
 out.Node=FEnode; out.Node(i1,:)=node;
 out.Elt=FEel0; out.Stack={};
else % RotateSel
 i5(3,:)=opt(3:5)/norm(opt(3:5));
 i5(4:6,1:3)=[0 -i5(3,3) i5(3,2);i5(3,3) 0 -i5(3,1);-i5(3,2) i5(3,1) 0];

 % i5 [origin;tx ty tz;rx ry rz];
 i6=[eye(3) i5(1,:)';0 0 0 1]* ...
    [cos(opt(2))*eye(3)+i5(3,:)'*i5(3,:)*(1-cos(opt(2)))+ ...
	     sin(opt(2))*i5(4:6,1:3) zeros(3,1);zeros(1,3) 1]* ...
    [eye(3) -i5(1,:)';0 0 0 1];
 %i6 = i6*[eye(3) i5(2,:)'*j1;0 0 0 1]; %translations

  out=feutil('repeatSel',FEnode,FEel0,[],i6);
end
  
  
 if ~isempty(model);
  model.Node=out.Node;model.Elt=out.Elt;
  model=stack_set(model,out.Stack);
  out=model;
 end

%% #Trans --------------------------------------------------------------------
elseif comstr(Cam,'trans');[CAM,Cam] = comstr(CAM,6);
%%  #TransSel dx dy dz - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if comstr(Cam,'sel')||comstr(Cam,'el')
 model=[]; % Model can be assigned below
 [carg,FEnode,FEel0,r1,ModelStack]=get_nodeelt(varargin,carg,ModelStack);
 if Cam(1)=='s'; [CAM,Cam] = comstr(CAM,4);else;[CAM,Cam] = comstr(CAM,3);end
 opt = comstr(CAM,-1);
 if isempty(opt)
   opt=varargin{carg}; carg=carg+1;
 end
 if length(opt)<3; opt(3)=0; end
 opt=opt(1:3); i6 = [eye(3) opt(:);0 0 0 1]; %translations

 out=feutil('repeatSel',FEnode,FEel0,[],i6);
 if ~isempty(model);
  model.Node=out.Node;model.Elt=out.Elt;
  model=stack_set(model,out.Stack);
  out=model;
 end
 
else;error('Trans%s unknown',CAM);
end

% ----------------------------------------------------------------------------
%% #SymSel Orig nx ny nz -----------------------------------------------------
elseif comstr(Cam,'symsel'); [CAM,Cam] = comstr(CAM,7);

 model=[]; % Model can be assigned below 
 [carg,FEnode,FEel0,r1,ModelStack]=get_nodeelt(varargin,carg,ModelStack);
 if comstr(Cam,'o')
    opt=comstr(CAM(2:length(CAM)),[-1 0 0 0 0 0 0]);
    i5=opt(1:3); opt=opt(3:6);
 else
    opt =  comstr(CAM,[-1 0 0 0 0]);i1=find(FEnode(:,1)==opt(1));
    if isempty(i1); i5=[0 0 0];else;i5=FEnode(i1(1),[5:7]); end
 end

 if norm(opt(2:4))==0; error('Need a non-zero direction for Sym');end
 i5(3,:)=opt(2:4)/norm(opt(2:4));
 i5(4:6,1:3)=[0 -i5(3,3) i5(3,2);i5(3,3) 0 -i5(3,1);-i5(3,2) i5(3,1) 0];

 % i5 [origin;0 0 0;tx ty tz;rx ry rz];
 i6=[eye(3) i5(1,:)';0 0 0 1]* ...
    [eye(3)-2*i5(3,:)'*i5(3,:) zeros(3,1);zeros(1,3) 1]* ...
    [eye(3) -i5(1,:)';0 0 0 1];
 %i6 = i6*[eye(3) i5(2,:)'*j1;0 0 0 1]; %translations
 out=feutil('repeatSel',FEnode,FEel0,[],i6);
 if ~isempty(model);
  model.Node=out.Node;model.Elt=out.Elt;
  model=stack_set(model,out.Stack);
  out=model;
 end
 
%% #Repeat -------------------------------------------------------------------
elseif comstr(Cam,'repeat'); [CAM,Cam] = comstr(CAM,7);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%%  #RepeatSel :  nITE tx ty tz Orig# AnleDeg nx ny nz ite0
if comstr(Cam,'sel');  [CAM,Cam]=comstr(CAM,4);
 model=[]; % Model can be assigned below
 [carg,FEnode,FEel0,r1,ModelStack]=get_nodeelt(varargin,carg,ModelStack);

 if carg<=nargin && all(size(varargin{carg})==4)
   i6=varargin{carg};carg=carg+1; opt(1)=-1;
 else
  if ~isempty(strfind(Cam,'o'))
    i1=strfind(Cam,'o');
    opt=comstr(CAM([1:i1-1 i1+1:length(CAM)]),[-1 1 0 0 0  0 0 0   0   1 0 0]);
    i5=opt(5:7); opt=[opt(1:4) 0 opt(8:10)];
  else
    opt=comstr(CAM,[-1  1   0  0  0   0      0     1  0  0]);
    i3=find(FEnode(:,1)==opt(5));
    if isempty(i3); i5(1,1:3)=[0 0 0];
    else;i5(1,1:3)=[FEnode(i3(1),5:7)];end
  end

  % Repeats element groups in FEel0
  elt=[];
  % i5 [origin;tx ty tz;rx ry rz];
  i5(2,1:3)=opt(2:4);
  i3=opt(7:9); if norm(i3)~=0; i3=i3/norm(i3); else;opt(6)=0; end;i5(3,:)=i3;
  i5(4:6,1:3)=[0 -i5(3,3) i5(3,2);i5(3,3) 0 -i5(3,1);-i5(3,2) i5(3,1) 0];
  opt(6)=opt(6)/180*pi;
    if opt(6)~=0 %rotations included
        i6=[eye(3) i5(1,:)';0 0 0 1]* ...
           [cos(opt(6))*eye(3)+i5(3,:)'*i5(3,:)*(1-cos(opt(6)))+ ...
	    sin(opt(6))*i5(4:6,1:3) zeros(3,1);zeros(1,3) 1]* ...
           [eye(3) -i5(1,:)';0 0 0 1];
    else;i6=eye(4);
    end
    i6 = i6*[eye(3) i5(2,:)';0 0 0 1]; %translations
 end

 [EGroup,nGroup]=getegroup(FEel0);
 elt=[]; r3=zeros(1,nGroup); i0=0;

 for jGroup=1:nGroup;

  ElemF= feutil('getelemf',FEel0(EGroup(jGroup),:),jGroup);
  elt = [elt;FEel0(EGroup(jGroup),:)];
  celle = FEel0(EGroup(jGroup):EGroup(jGroup+1)-1,:);

  % define celle node numbers in terms of indices in FEnode
  iNode=fe_super('node',ElemF);iNode = iNode(:)';
  NNode=sparse(FEnode(:,1),1,1:length(FEnode(:,1)));

  % i1 indices in FEnode of all the needed nodes, i2 nodes without repeat
  i1 = celle(2:size(celle,1),iNode);
  i1 = reshape(NNode(i1(:)),size(i1,1),size(i1,2));
  i2 = find(sparse(i1(:),1,i1(:))); ci2(i2) = 1:length(i2);
  % loop on nITE
  j2 = 1:opt(1); if opt(1)==-1; j2 = 2; end % this is for translations
  % r1 contains current cell positions as a quadrivector
  r1=[FEnode(i2,5:7) ones(length(i2),1)]; if opt(1)==-1; r1 = (i6*r1')'; end
  if length(opt)>=10&&opt(10)<0 % % move back if needed
    for j1=1:abs(opt(10));r1=(i6\r1')';end
  end  % move back if needed
  i0=i0(end)+[1:size(r1,1)];
  for j1 = j2-1
    [FEnode,i4]=feutil('AddNode',FEnode,r1(:,1:3));
    r3(i0,j1+1)=FEnode(i4,1); % set node group value to cell number
    i4=reshape(FEnode(i4(ci2(i1)),1),size(i1,1),size(i1,2));
    i5 = celle(2:size(celle,1),:); i5(:,iNode) = i4; 
    elt = [elt; i5];
    r1 = (i6*r1')'; % modified cell nodes
  end

 end % jGroup of loop on groups
 out=struct('Node',FEnode,'Elt',elt,'Stack',[]);
 out.Stack={'info','CellNodes',unique(r3,'rows')};
 if ~isempty(model);
  model.Node=out.Node;model.Elt=out.Elt;
  model=stack_set(model,out.Stack);
  out=model;
 end
else;sdtw('Repeat%s unknown',CAM);
end % subcommand selection - - - - - -

% ----------------------------------------------------------------------------
%% #Rev : REVOLUTION functions REV -------------------------------------------
elseif comstr(Cam,'rev');  [CAM,Cam] = comstr(CAM,4);
 
 model=[]; % Model can be assigned below
 [carg,FEnode,FEel0,r1,ModelStack]=get_nodeelt(varargin,carg,ModelStack);
 if ~isempty(strfind(Cam,'-trans'))
  i1=strfind(Cam,'-trans');
  RunOpt=struct('type','trans');CAM(i1+[0:5])='';[CAM,Cam] = comstr(CAM,1);
 elseif ~isempty(strfind(Cam,'-tria'))
  i1=strfind(Cam,'-tria');
  RunOpt=struct('type','tria');CAM(i1+[0:4])='';[CAM,Cam] = comstr(CAM,1);
 else;RunOpt=struct('type','');
 end
 [CAM,Cam,RunOpt.OptimDegen]=comstr('optimdegen',[-25 3],CAM,Cam);
 [CAM,Cam,RunOpt.KnownNew]=comstr('knownnew',[-25 3],CAM,Cam);
 [CAM,Cam,RunOpt.NoOrient]=comstr('noorient',[-25 3],CAM,Cam);
 % Rev               nRep Orig# AngleDeg nx ny nz tx ty tz
 if ~isempty(strfind(Cam,'o'))
   i1 = strfind(Cam,'o');
   opt = comstr(CAM([1:i1-1 i1+1:length(CAM)]), ...
             [-1 10 0 0 0  360    0  0  1  0  0  0]);
   i4 = [1 0 0 0 opt(2:4)];opt=[opt(1) 0 opt(5:length(opt))];
 else
  opt = comstr(CAM,[-1 10    0     360    0  0  1  0  0  0]);
  if opt(1)==1&&opt(3)==360; error('Rev: nRep=1 and Angle=360 is not valid');end
  i4=find(FEnode(:,1)==opt(2));
  if isempty(i4); i4 = [1  0 0 0  0 0 0]; else;i4 = FEnode(i4,:); end
 end


 if carg<=nargin&&isa(varargin{carg},'double')
  ind=varargin{carg};carg=carg+1;ind=ind(:)';
  i7 =ind;  opt(1)=length(ind)-1;
 elseif carg<=nargin&&isfield(varargin{carg},'map')
  RunOpt=sdth.sfield('addmissing',varargin{carg},RunOpt);carg=carg+1;
  ind=[0 1];i7=[0 1];opt(3)=0;
 else;ind=linspace(0,1,opt(1)+1); i7=0:opt(1);
 end

 % revolution/extrusion of the groups in FEel0
 [EGroup,nGroup]=getegroup(FEel0);
 elt=[];RunOpt.ind0=ind;RunOpt.i7=i7;
 RunOpt.opt0=opt;
 for jGroup=1:nGroup
  [ElemF,opt1,ElemP]= feutil('getelemf',FEel0(EGroup(jGroup),:),jGroup);
  cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;
  ind=RunOpt.ind0;i7=RunOpt.i7;opt=RunOpt.opt0;
  
  if     strcmp(ElemP,'beam1') 
     iNode=[1 2];elt(size(elt,1)+1,1:6)=[Inf abs('quad4')]; %#ok<AGROW>
     cEGI(~all(FEel0(cEGI,iNode),2))=[];if isempty(cEGI); continue;end
     r1=FEel0(cEGI,[3:4]);
  % transition ring using tria6
  elseif     strcmp(ElemP,'beam3') &&strcmp(RunOpt.type,'trans');
     iNode=[1 2 3];elt(size(elt,1)+1,1:6)=[Inf abs('tria6')]; %#ok<AGROW>
     ind=sort([ind ind(1:end-1)+diff(ind)/2]);
     i7=[i7(:) i7(:)+[diff(i7(:))/2;0]]';i7=i7(:);i7=i7(1:end-1)';
     opt(6:9)=opt(6:9)/2;opt(1)=opt(1)*2;
     if size(FEel0,2)<10; FEel0(1,10)=0;end;r1=FEel0(cEGI,[4:5]);
  elseif     strcmp(ElemP,'beam3') &&strcmp(RunOpt.type,'tria');
     iNode=[1 2 3];elt(size(elt,1)+1,1:6)=[Inf abs('tria6')]; %#ok<AGROW>
     ind=sort([ind ind(1:end-1)+diff(ind)/2]);
     i7=[i7(:) i7(:)+[diff(i7(:))/2;0]]';i7=i7(:);i7=i7(1:end-1)';
     opt(1)=opt(1)*2;%opt(6:9)=opt(6:9)/2;
     if size(FEel0,2)<10; FEel0(1,10)=0;end;r1=FEel0(cEGI,[4:5]);
  elseif     strcmp(ElemP,'beam3') 
     iNode=[1 2 3];elt(size(elt,1)+1,1:6)=[Inf abs('quadb')]; %#ok<AGROW>
     ind=sort([ind ind(1:end-1)+diff(ind)/2]);
     i7=[i7(:) i7(:)+[diff(i7(:))/2;0]]';i7=i7(:);i7=i7(1:end-1)';
     opt(1)=opt(1)*2;%opt(6:9)=opt(6:9)/2;
     if size(FEel0,2)<10; FEel0(1,10)=0;end;r1=FEel0(cEGI,[4:5]);
  elseif     strcmp(ElemP,'mass1') 
     iNode=[1];elt(size(elt,1)+1,1:6)=[Inf abs('beam1')]; %#ok<AGROW>
     r1=ones(length(cEGI),1)*[1 1];
  elseif strcmp(ElemP,'tria3') 
     iNode=[1:3];elt(size(elt,1)+1,1:7)=[Inf abs('penta6')]; %#ok<AGROW>
     r1=FEel0(cEGI,4:5);
  elseif strcmp(ElemP,'tria6') 

     iNode=[1:6];elt(size(elt,1)+1,1:8)=[Inf abs('penta15')]; %#ok<AGROW>
     ind=sort([ind ind(1:end-1)+diff(ind)/2]);
     i7=[i7(:) i7(:)+[diff(i7(:))/2;0]]';i7=i7(:);i7=i7(1:end-1)';
     opt(1)=opt(1)*2;%opt(6:9)=opt(6:9)/2; corrected 2009
     if size(FEel0,2)<17; FEel0(1,17)=0;end
     r1=FEel0(cEGI,7:8);

  elseif strcmp(ElemP,'quad4')
     iNode=[1:4];elt(size(elt,1)+1,1:6)=[Inf abs('hexa8')]; %#ok<AGROW>
     r1=FEel0(cEGI,5:6);
  elseif strcmp(ElemP,'quadb')
     iNode=[1:8];elt(size(elt,1)+1,1:7)=[Inf abs('hexa20')]; %#ok<AGROW>
     ind=sort([ind ind(1:end-1)+diff(ind)/2]);
     i7=[i7(:) i7(:)+[diff(i7(:))/2;0]]';i7=i7(:);i7=i7(1:end-1)';
     opt(1)=opt(1)*2;%opt(6:9)=opt(6:9)/2; corrected 2009
     if size(FEel0,2)<10; FEel0(1,10)=0;end
     r1=FEel0(cEGI,9:10);
  elseif strcmp(ElemP,'quad9'); iNode=[1:4];
     sdtw('_nb','quad9 extruded into hexa8');
     iNode=[1:4];elt(size(elt,1)+1,1:6)=[Inf abs('hexa8')]; %#ok<AGROW>
     r1=FEel0(cEGI,10:11);
  else;                    error([ElemP ' not supported by Rev']);
  end

  NNode=sparse(FEnode(:,1),1,1:length(FEnode(:,1)));
  % i1 angles for revolution meridians i3 revolution basis i4 origin
  i1 = -ind*opt(3)*pi/180;

  % i2 initial nodes
  i2=NNode(FEel0(cEGI,iNode)');i2 = FEnode(i2,5:7);
  i6 = [1:size(i2,1)]'*ones(1,length(ind));
  i7 = i7(ones(1,size(i2,1)),:);  i1=i1(ones(1,size(i2,1)),:);

  if opt(3)~=0 %rotate the nodes
    i2 = i2 - i4(ones(length(iNode)*length(cEGI),1),5:7);
    
    i3 = basis(opt(4:6),[1 0 0],1);    i3 = i3(:,[2 3 1]);
    i2 = i2*i3;
    i2(:,1:2) = [sqrt(i2(:,1).^2+i2(:,2).^2) atan2(i2(:,1),i2(:,2))];

    % i2 rotated nodes
    i2 = i4(ones(size(i2,1)*(opt(1)+1),1),5:7) +  ...  % origin
         i2(i6,3)*i3(:,3)' + ... %vertical pos
         (i2(i6,1)*ones(1,3)) .* ...
                           (sin(i1(:)+i2(i6,2)) *i3(:,1)' + ...
                            cos(i1(:)+i2(i6,2)) *i3(:,2)');
  else 
    i2=i2(i6,:);
  end

  if size(r1,1)==length(cEGI) % propagate properties
   r1=[r1(:,ones(1,length(ind))) r1(:,2*ones(1,length(ind)))];
   r1=reshape(r1,size(r1,1)*length(ind),2);
  end

  if isfield(RunOpt,'map') % Generate the extrusion layers
    i3=repmat(RunOpt.map.normal(full(RunOpt.map.nind(FEel0(cEGI,iNode)')),:),length(ind),1);
    i2=i2+diag(sparse(i7(:)))*i3;     
  else
   i2=i2+i7(:)*opt(1,7:9); %translate nodes
  end
  if 1==1
   if ~isempty(RunOpt.KnownNew)&&RunOpt.KnownNew % no control allow zero extrudes and unclean degen
    % first zero occurence in i7 should be kept
    in1=find(~any(i7),1,'first');
    if ~isempty(in1)
     in2=(size(i2,1)*(in1-1))+(1:size(i7,1)); in3=setdiff(1:size(i2,1),in2);
     i3=i2(:,1);
     [FEnode,i3(in2)]=feutil('AddNode',FEnode,i2(in2,:));
     [FEnode,i3(in3)]=feutil('AddNodeNew',FEnode,i2(in3,:));
     i2=i3;
    else; [FEnode,i2]=feutil('AddNodeNew',FEnode,i2);
    end
   else % clean
      [r2,RO.i1,RO.i2]=unique(i2,'rows');
      [FEnode,i2]=feutil('AddNode',FEnode,r2);
      i2=i2(RO.i2);
   end
  else
   [FEnode,i2]=feutil('AddNode',FEnode,i2);
  end
  i2=reshape(FEnode(i2,1),length(iNode),length(i2)/length(iNode))';
  i3=1:size(i2,1)-length(cEGI);  i6=size(elt,1);

    if strcmp(ElemP,'beam1')
      i5=find(i2(i3,1)-i2(i3,2));
      if length(i5)~=length(i3)
       warning('Removing degenerate beam from extrusion');
       i3=i3(i5); % xxx deal with MatId propagation
      end
      elt(i6+[1:length(i3)],1:6)= ...
              [i2(i3,1:2) i2(i3+length(cEGI),[2 1]) r1(i3,:)]; %#ok<AGROW>
    elseif strcmp(ElemP,'beam3') &&strcmp(RunOpt.type,'trans');
     i5=1:4:size(i2,1);
     for j1=i5(1:end-1)-1
      elt(end+[1:3],1:6)=[
      i2(j1+1,[1 2]) i2(j1+3,1) i2(j1+1,3) i2(j1+2,[3 1]) 
      i2(j1+1,2) i2(j1+5,2) i2(j1+3,1) i2(j1+3,2)  i2(j1+4,3)  i2(j1+2,3)   
     i2(j1+3,1) i2(j1+5,[2 1]) i2(j1+4,3) i2(j1+5,3) i2(j1+4,1) ]; %#ok<AGROW>
    end

    elseif strcmp(ElemP,'beam3') &&strcmp(RunOpt.type,'tria');

      i3=[length(cEGI):2*length(cEGI):size(i2,1)-1]';i5=[-length(cEGI)+1:0];
      i3=i3(:,ones(size(i5)))+i5(ones(size(i3)),:);i3=i3(:);i5=length(cEGI);

      elt(i6+[1:2*length(i3)],1:8)= ...
       [i2(i3,1:2) i2(i3+2*i5,1)    i2(i3,3) i2(i3+i5,3)  ...
        i2(i3+i5,1) r1(i3,:) ;
        i2(i3,2) i2(i3+2*i5,[2 1])  i2(i3+i5,2) i2(i3+2*i5,3) ...
        i2(i3+i5,3) r1(i3,:)]; %#ok<AGROW>
    
    elseif strcmp(ElemP,'beam3')
      i3=[length(cEGI):2*length(cEGI):size(i2,1)-1]';i5=[-length(cEGI)+1:0];
      i3=i3(:,ones(size(i5)))+i5(ones(size(i3)),:);i3=i3(:);i5=length(cEGI);

      elt(i6+[1:length(i3)],1:10)= [i2(i3,1:2) i2(i3+2*i5,[2 1]) ...
                              i2(i3,3) i2(i3+i5,2) i2(i3+2*i5,3) i2(i3+i5,1) ...
                               r1(i3,:)]; %#ok<AGROW>
    elseif strcmp(ElemP,'tria6')

      i3=[length(cEGI):2*length(cEGI):size(i2,1)-1]';i5=[-length(cEGI)+1:0];
      i3=i3(:,ones(size(i5)))+i5(ones(size(i3)),:);i3=i3(:);i5=length(cEGI);
      elt(i6+[1:length(i3)],1:17)= [i2(i3,1:3) i2(i3+2*i5,1:3) ...
                                    i2(i3,4:6) i2(i3+i5,1:3) ...
                                    i2(i3+2*i5,4:6)  r1(i3,:)]; %#ok<AGROW>

    elseif strcmp(ElemP,'mass1')
      elt(i6+[1:length(i3)],1:4)= ...
              [i2(i3,1) i2(i3+length(cEGI),1) r1(i3,:)]; %#ok<AGROW>
    elseif strcmp(ElemP,'tria3')
      elt(i6+[1:length(i3)],1:8)= ...
              [i2(i3,:) i2(i3+length(cEGI),:) r1(i3,:)]; %#ok<AGROW>

    elseif strcmp(ElemP,'quad4')||strcmp(ElemP,'quad9')
      elt(i6+[1:length(i3)],1:10)= ...
              [i2(i3,:) i2(i3+length(cEGI),:) r1(i3,:)]; %#ok<AGROW>

    elseif strcmp(ElemP,'quadb')
      i3=[length(cEGI):2*length(cEGI):size(i2,1)-1]';i5=[-length(cEGI)+1:0];
      i3=i3(:,ones(size(i5)))+i5(ones(size(i3)),:);i3=i3(:);i5=length(cEGI);
      elt(i6+[1:length(i3)],1:22)= ...
       [i2(i3,1:4) i2(i3+2*i5,1:4) ... % corner nodes
              i2(i3,5:8) i2(i3+i5,1:4) i2(i3+2*i5,5:8) r1(i3,:)]; %#ok<AGROW>
      %              r1(1:length(i3),:)]; was modified everywhere
    end
 end % of loop on jGroup

 if ~RunOpt.NoOrient % Used in sdtweb fe_shapeoptim SDMThickSensMap for negative volume
  [elt,i1,r1]=feutil('orient-back;',FEnode,elt);
 end
 
 out=struct('Node',FEnode,'Elt',elt,'Stack',[]);
 if RunOpt.OptimDegen % Check and remove (or transform) some degenerated elts
   out=feutil('OptimDegen',out);
 end
 if ~isempty(model);
  model.Node=out.Node;model.Elt=out.Elt;
  model=stack_set(model,out.Stack);
  out=model;
 end
 
%% #Extrude ------------------------------------------------------------------
elseif comstr(Cam,'extrude'); [CAM,Cam]=comstr(CAM,8);

 % extrude nRep tx ty tz
 [CAM,Cam,knew]=comstr('knownnew',[-25 3],CAM,Cam);
 if knew; knew='KnownNew'; else; knew=''; end
 opt=comstr(CAM,[-1 1 0 0 1]);
 out=feutil(sprintf('Rev %s %i 0 0 0 0 0 %.20g %.20g %.20g',knew,opt(1:4)), ...
            varargin{2:end});

%% #Sel(node,elt) #SelElt, #Selnode-----------------------------------------
elseif comstr(Cam,'sel');  [CAM,Cam] = comstr(CAM,4);


  if comstr(Cam,'elt{')
    % Sequence of selections
    % cf.sel={'{eltname~=flui & seledge &innode{x<y-515},eltnameflu&inode{x<y-515}}','showfipro'}
    [~,RB]=sdtm.urnPar(CAM,'{}{}');out=[];
    for j1=1:length(RB.Other)
     [r1,elt]=feutil(['findelt' RB.Other{j1}],varargin{2:end});
     if j1==1;out=elt;
     else; out(end+(1:size(elt,1)),1:size(elt,2))=elt;
     end
    end
  else
   [r1,out]=feutil(['find' CAM],varargin{2:end});  
  end


%% #Set ----------------------------------------------------------------------
elseif comstr(Cam,'set');  [CAM,Cam] = comstr(CAM,4);

 % #SetGroup - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if comstr(Cam,'group'); [CAM,Cam] = comstr(CAM,6);

  [carg,FEnode,elt,el0,ModelStack]=get_nodeelt(varargin,carg,ModelStack);
  if carg<=nargin; st1=varargin{carg};carg=carg+1; else;st1='elt';end
  [EGroup,nGroup]=getegroup(elt);

  RunOpt=struct('MatId',[],'ProId',[],'EGID',[]);
  i2=[strfind(Cam,'mat') length(Cam)+1];i1(1)=i2(1);
  i2=[strfind(Cam,'pro') strfind(Cam,'sec') length(Cam)+1];i1(2)=i2(1);
  i2=[strfind(Cam,'eg') length(Cam)+1];i1(3)=i2(1);
  i2=[strfind(Cam,'name') length(Cam)+1];i1(4)=i2(1);
  i1(5)=length(Cam)+1;

  % select group
  if ~isletter(Cam(1)); opt=comstr(Cam(1:min(i1)-1),[-1]);else;opt=[];end
   if isempty(opt & min(i1)~=1) % group selected by name
    [st1,st] = comstr(CAM(1:min(i1)-1),1);
    if strcmp(st,'all'); opt=1:nGroup;
    else
     for jGroup=1:nGroup
      [ElemF]= feutil('getelemf',elt(EGroup(jGroup),:));
      if comstr(ElemF,st); opt(end+1)=jGroup; end
     end
    end
   end
  opt=fix(opt);opt=opt(opt<=nGroup);%contains selected groups 
  if isempty(opt); error(['SetGroup: not a group in ' st1]); end

  i2=i1(4)+4:min([i1(i1>i1(1)) i1(end)])-1; % mat
  if i1(4)<i1(5) % name
   [i5,i2,i3,i4]=sscanf(CAM(i1(4)+4:length(CAM)),'%s',1);
   CAM(i2)=' ';
  else;i5=[];
  end

  i2=i1(1)+3:min([i1(i1>i1(1)) i1(end)])-1; % mat
  if ~isempty(i2) 
   i3=comstr(comstr(Cam(i2),'id','%c'),[-1 0]);RunOpt.MatId=fix(i3(1));
  else;i3=0;
  end
  i2=i1(2)+3:min([i1(i1>i1(2)) i1(end)])-1;% pro
  if ~isempty(i2) 
   i4 = comstr(comstr(Cam(i2),'id','%c'),[-1 0]);RunOpt.ProId=fix(i4(1));
  else;i4=0;
  end
  i2=i1(3)+2:min([i1(i1>i1(3)) i1(end)])-1;% egid
  if ~isempty(i2) 
   i8=comstr(comstr(Cam(i2),'id','%c'),[-1 0]);RunOpt.EGID=fix(i8(1));
  end

  for jGroup = opt %loop on element groups
   [ElemF,i7]= feutil('getelemf',elt(EGroup(jGroup),:),jGroup);
   [i6,SEopt]=fe_super('prop',ElemF);i6=[i6(:);0;0;0];
   cEGI = [EGroup(jGroup)+1:EGroup(jGroup+1)-1]';
   if ~isempty(RunOpt.MatId);
    if i6(1); elt(cEGI,i6(1))=ones(size(cEGI))*RunOpt.MatId;
    else; warning(sprintf('MatId not defined for %s',ElemF));end
   end
   if ~isempty(RunOpt.ProId);
    if i6(2); elt(cEGI,i6(2))=ones(size(cEGI))*RunOpt.ProId; 
    elseif SEopt(1)==1; fesuper(sprintf('set %s ProId %i',ElemF,RunOpt.ProId))
    else; warning(sprintf('ProId not defined for %s',ElemF));
    end
   end
   if ~isempty(RunOpt.EGID) % egid is given
     i7(1)=RunOpt.EGID;elt(EGroup(jGroup),1:2+length(ElemF)+length(i7))=...
      [Inf abs(ElemF) 0 i7];
   end
   if ~isempty(i5)
    i7 = elt(EGroup(jGroup),2+length(ElemF):end);
    if isempty(i7);  elt(EGroup(jGroup),2:end)=0;
                    elt(EGroup(jGroup),2:length(i5)+1)=abs(i5);
    else;elt(EGroup(jGroup),2:length(i5)+2+length(i7))=[abs(i5) 0 i7];
    end
   end
  end

  st=['For group(s) ' sprintf('%i ',opt)];
  if ~isempty(RunOpt.MatId); st=[st sprintf('  MatId set to %i  ',i3)]; end
  if ~isempty(RunOpt.ProId); st=[st sprintf('  ProId set to %i  ',i4)]; end
  if ~isempty(RunOpt.EGID); st=[st sprintf('  EGID  set to %i  ',i8)]; end
  if ~isempty(i5); st=[st sprintf('  Name set to %s  ',i5)]; end
  out1=st;
  out=elt;

 %% #SetPro used to define element properties ->sdtweb fe_mat('setpro') - - - - - - -
 elseif comstr(Cam,'pro');
  eval(iigui({'fe_mat',nargout},'OutReDir')); % now in fe_mat

 %% #SetMat used to define material properties - - - - - - - - - - - - - - - - 
 elseif comstr(Cam,'mat');
  eval(iigui({'fe_mat',nargout},'OutReDir')); % now in fe_mat
 %% #SetSel - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'sel');
  i1=regexpi(Cam,'(matid|mat)\s*(\d*)','tokens');
  if ~isempty(i1);RO.MatId=str2double(i1{1}{2});else; RO.MatId=[];end
  i1=regexpi(Cam,'(proid|pro)\s*(\d*)','tokens');
  if ~isempty(i1);RO.ProId=str2double(i1{1}{2});else;RO.ProId=[];end
  model=[];[carg,FEnode,elt,el0,ModelStack]=get_nodeelt(varargin,carg,ModelStack);
  if ~isstruct(model);
    error('Expecting a model structure');
  end
  if carg<=nargin; st1=varargin{carg};carg=carg+1; 
  else;error('Expecting a selection string');
  end
  mpid=feutil('mpid',model);
  i1=feutil(['findelt' st1],model);
  if ~isempty(RO.MatId);mpid(i1,1)=RO.MatId;end
  if ~isempty(RO.ProId);mpid(i1,2)=RO.ProId;end
  st='mpid'; if ~isempty(strfind(Cam,'force')); st=[st 'force']; end
  out=feutil(st,model,mpid);
 %% SetEnd
 else;error('Set%s not a valid command',CAM);
 end

%% #String -------------------------------------------------------------------
elseif comstr(Cam,'s');  CAM=comstr(CAM,'string','%c');Cam=lower(CAM);

%% #StringDof - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if comstr(Cam,'d'); CAM=comstr(CAM,'dof','%s');Cam=lower(CAM);

    if 1==2
     for j1=99:-1:13; st{j1}=sprintf('.%i',j1);end
     st([1:12 19 20 21])={'x','y','z','\theta x','\theta y','\theta z', ...
        '-x','-y','-z','-\theta x','-\theta y','-\theta z','p','T','V'}';
     disp(comstr(st,-30,struct('NoCLip',2)))
    else;      
     if ~isempty(Cam)&&Cam(end)==':';
      % From fe_sens : lab
      st={'+X','+Y','+Z','RX','RY','RZ', ...
       '-X','-Y','-Z','-RX','-RY','-RZ','Fx','Fy','Fz', ...
       'Mx','My','Mz','p','T','V','.22','.23','.24','.25','.26','.27', ...
       '.28','.29','.30','.31','.32','.33','.34','.35','.36','.37', ...
       '.38','.39','.40','.41','.42','.43','.44','.45','.46','.47','.48', ...
       '.49','.50','.51','.52','.53','.54','.55','.56','.57','.58','.59', ...
       '.60','.61','.62','.63','.64','.65','.66','.67','.68','.69','.70', ...
       '.71','.72','.73','.74','.75','.76','.77','.78','.79','.80','.81', ...
       '.82','.83','.84','.85','.86','.87','.88','.89','.90','.91','.92', ...
       '.93','.94','.95','.96','.97','.98','.99'};
      st1=st;
      st2='%i:%s';
     elseif comstr(Cam,'_d') % stringdof_f scalar disp labels
      % From fe_sens : lab
      st={'x','+y','z','rx','ry','rz', ...
       '-x','-y','-z','-rx','-ry','-rz','Fx','Fy','Fz', ...
       'Mx','My','Mz','p','HFLU','V','.22','.23','.24','.25','.26','.27', ...
       '.28','.29','.30','.31','.32','.33','.34','.35','.36','.37', ...
       '.38','.39','.40','.41','.42','.43','.44','.45','.46','.47','.48', ...
       '.49','.50','.51','.52','.53','.54','.55','.56','.57','.58','.59', ...
       '.60','.61','.62','.63','.64','.65','.66','.67','.68','.69','.70', ...
       '.71','.72','.73','.74','.75','.76','.77','.78','.79','.80','.81', ...
       '.82','.83','.84','.85','.86','.87','.88','.89','.90','.91','.92', ...
       '.93','.94','.95','.96','.97','.98','.99'};
      st1=st;
      st2='%i:%s';
         
     elseif comstr(Cam,'_f') % stringdof_f scalar force labels
     st={'Fx','Fy','Fz','Mx','My','Mz','-Fx','-Fy','-Fz','-Mx','-My','-Mz', ...
      'Fx','Fy','Fz', ...
      'Mx','My','Mz','av','T','V','.22','.23','.24','.25','.26','.27', ...
      '.28','.29','.30','.31','.32','.33','.34','.35','.36','.37', ...
      '.38','.39','.40','.41','.42','.43','.44','.45','.46','.47','.48', ...
      '.49','.50','.51','.52','.53','.54','.55','.56','.57','.58','.59', ...
      '.60','.61','.62','.63','.64','.65','.66','.67','.68','.69','.70', ...
      '.71','.72','.73','.74','.75','.76','.77','.78','.79','.80','.81', ...
      '.82','.83','.84','.85','.86','.87','.88','.89','.90','.91','.92', ...
      '.93','.94','.95','.96','.97','.98','.99'};
      st2='%i:%s';
     elseif comstr(Cam,'_v') % stringdof_v scalar velocity labels
     st={'vx','vy','vz','vrx','vry','vrz','-vx','-vy','-vz','-vrx','-vry','-vrz', ...
      'vFx','vFy','vFz', ...
      'vMx','vMy','vMz','vp','vT','vV','.22','.23','.24','.25','.26','.27', ...
      '.28','.29','.30','.31','.32','.33','.34','.35','.36','.37', ...
      '.38','.39','.40','.41','.42','.43','.44','.45','.46','.47','.48', ...
      '.49','.50','.51','.52','.53','.54','.55','.56','.57','.58','.59', ...
      '.60','.61','.62','.63','.64','.65','.66','.67','.68','.69','.70', ...
      '.71','.72','.73','.74','.75','.76','.77','.78','.79','.80','.81', ...
      '.82','.83','.84','.85','.86','.87','.88','.89','.90','.91','.92', ...
      '.93','.94','.95','.96','.97','.98','.99'};
      st2='%i:%s';
         
     else
      % Model
     st={'x','y','z','\theta x','\theta y','\theta z', ...
      '-x','-y','-z','-\theta x','-\theta y','-\theta z','Fx','Fy','Fz', ...
      'Mx','My','Mz','p','T','V','.22','.23','.24','.25','.26','.27', ...
      '.28','.29','.30','.31','.32','.33','.34','.35','.36','.37', ...
      '.38','.39','.40','.41','.42','.43','.44','.45','.46','.47','.48', ...
      '.49','.50','.51','.52','.53','.54','.55','.56','.57','.58','.59', ...
      '.60','.61','.62','.63','.64','.65','.66','.67','.68','.69','.70', ...
      '.71','.72','.73','.74','.75','.76','.77','.78','.79','.80','.81', ...
      '.82','.83','.84','.85','.86','.87','.88','.89','.90','.91','.92', ...
      '.93','.94','.95','.96','.97','.98','.99'};
      st1={'x','y','z','rx','ry','rz', ...
       '-x','-y','-z','-rx','-ry','-rz','Fx','Fy','Fz', ...
       'Mx','My','Mz','p','T','V','.22','.23','.24','.25','.26','.27', ...
       '.28','.29','.30','.31','.32','.33','.34','.35','.36','.37', ...
       '.38','.39','.40','.41','.42','.43','.44','.45','.46','.47','.48', ...
       '.49','.50','.51','.52','.53','.54','.55','.56','.57','.58','.59', ...
       '.60','.61','.62','.63','.64','.65','.66','.67','.68','.69','.70', ...
       '.71','.72','.73','.74','.75','.76','.77','.78','.79','.80','.81', ...
       '.82','.83','.84','.85','.86','.87','.88','.89','.90','.91','.92', ...
       '.93','.94','.95','.96','.97','.98','.99'};
      st2='%i%s';
     end
    end
    
    r1 = varargin{carg};carg=carg+1; r1=r1(:);map='';
%  [n2,st2]=sdth.urn('nmap.Node',mo1,fix(sens.tdof));
    if isa(r1,'cell')
      out=zeros(length(r1),1); 
      for j1=1:length(r1)
        [i2,st1] = comstr(r1{j1},'','%i');
        i3=strmatch(st1,st); if isempty(i3)||i3(1)>100; i3=99;end
        out(j1,1)=i2+i3(1)/100;
      end
    else
      out=cell(length(r1),1); 
      i2=abs(r1);i2 = [fix(i2) round(rem(i2,1)*100)];
      ind=[r1<0 i2(:,1)<1&i2(:,2)>0 i2(:,2)==0 i2(:,1)<0]; % Each case
%       for j1=1:3; ind(:,j1+1)=ind(:,j1+1)&~any(ind(:,1:j1); end
      % Case1
      %test{j1}=sprintf('e%i.%3i\n',i2(ind(:,1),1),round(rem(-r1(ind(:,1)),1)*1000));
      %Case2
      if carg<=nargin&&isa(varargin{carg},'containers.Map')
       map=varargin{carg};carg=carg+1; 
       if strcmpi(map.KeyType,'char');
        map=containers.Map(map.values,map.keys);
       end
      elseif carg<=nargin&&isfield(varargin{carg},'nmap')
        map=varargin{carg}.nmap;carg=carg+1;
        if map.isKey('Map:Nodes');
            map=map('Map:Nodes'); 
        else;map=[];
        end
      end
      if ~isempty(map)
        if ~strcmp(map.KeyType,'int32');map=containers.Map(map.values,map.keys);end
        st=st1;
      end

      % To speed-up, use strsplit(sprintf(xxx),'\n') rather than for loop
      i2=abs(r1);i2 = [fix(i2) round(rem(i2,1)*100)];
      
      ind2=false(size(r1,1),1); % Indices of already treated dofs
      i3=~ind2 & r1<0; % Not treated+negative dof
      if any(i3)
       r3=strsplit(sprintf('e%i.%3i\n',i2(i3,1),round(rem(-r1(i3),1)*1000)),'\n');
       out(i3)=r3(1:end-1); ind2=ind2|i3;
      end
      i3=~ind2 & i2(:,1)<1&i2(:,2)>0; % Not treated+negative id+positive decimal
      if any(i3)
       out(i3)=st(i2(i3,2));ind2=ind2|i3;
      end
      i3=~ind2 & i2(:,2)==0; % Not treated+integer
      if any(i3)
       r3=strsplit(sprintf('%i*\n',i2(i3,1)),'\n');
       out(i3)=r3(1:end-1); ind2=ind2|i3;
      end
      if ~isempty(map); % use map to generate 
       i3=~ind2;
       r2=[map(int32(i2(i3,1))) st(i2(i3,2))]';
       r3=strsplit(sprintf('%s:%s\n',r2{:}),'\n');
       out(i3)=r3(1:end-1); ind2=ind2|i3;
      end
      i3=~ind2; % Not treated (all dofs whith standard 1.01 and no map)
      r2=[num2cell(i2(i3,1)) st(i2(i3,2))']';
      r3=strsplit(sprintf([st2 '\n'],r2{:}),'\n');
      out(i3)=r3(1:end-1); ind2=ind2|i3;

% Old for loop below
%       for j1=1:length(r1)
%         i2=abs(r1(j1));i2 = [fix(i2) round(rem(i2,1)*100)];
%         if r1(j1)<0
%            out{j1}=sprintf('e%i.%3i',i2(1),round(rem(-r1(j1),1)*1000));
%         elseif i2(1)<1&&i2(2)>0; out{j1}=sprintf('%s',st{i2(2)});
%         elseif i2(2)==0; out{j1}=sprintf('%i*',i2(1));
%         elseif i2(1)<0;out{j1}=sprintf('e%i.%3i',abs(i2));
%         elseif ~isempty(map);
%           out{j1}=sprintf('%s:%s',map(int32(i2(1))),st{i2(2)});
%         else 
%           out{j1}=sprintf(st2,i2(1),st{i2(2)});
%         end
%       end
    end
    if ~isempty(strfind(Cam,'-str'))&&isscalar(out); out=out{1}; end

% StringIODOF FLOR:180:+Z / FRNT:15:+Z
elseif comstr(Cam,'iodof')
 
else;error('String%s unknown',CAM);
end % subcommand selection - - - - - - - - - - - - - -

%% #Trace2Elt ----------------------------------------------------------------
elseif comstr(Cam,'trace2elt'); CAM=comstr(CAM,'trace2elt','%s');Cam=lower(CAM);

   elt=varargin{carg};carg=carg+1;
   r1=[Inf abs('beam1') 0 -1];r2=[Inf abs('mass2') 0 -1]; 
   out1={};
   for j1 = 1:size(elt,1)
        i2=[elt(j1,83:end-1);elt(j1,84:end)];
        i2 = i2(:,all(i2))';
        i4=[0 elt(j1,83:end) 0];
        i3=find(i4(1:end-2)==0&i4(2:end-1)~=0&i4(3:end)==0);
        if ~isempty(i3); i2=[i2;i4(i3+1)'*[1 1]]; end
        i3=find(i2(:,1)~=i2(:,2)); 
        if ~isempty(i3)
          r1(end+[1:size(i3,1)],1:4)=[i2(i3,:) ones(size(i3))*[j1 elt(j1,2)]];
        end
        i3=find(i2(:,1)==i2(:,2)); 
        if ~isempty(i3)
          r2(end+[1:size(i3,1)],[1 14 15])=[i2(i3,1) ...
                                        ones(size(i3))*[j1 elt(j1,2)]];
        end
        out1{j1,1}=deblank(char(elt(j1,3:82)));
   end; 

   if size(r1,1)>1; out = [r1]; else;out=[]; end
   if size(r2,1)>1; if size(out,2)<15; out(1,15)=0; end;out = [out;r2]; end
   
%% #Unjoin -------------------------------------------------------------------   
elseif comstr(Cam,'unjoin'); [CAM,Cam] = comstr(CAM,7);
    
 model=varargin{carg};carg=carg+1; %xxx check struct
 
 RunOpt=struct('type','group','NodeSel','groupall');
 sel1 = fix(comstr(CAM,-1)); % Unjoin 1 2 (group=
 if isempty(sel1)  
  sel1=varargin{carg};carg=carg+1;
  if isstruct(sel1) % Unjoin,model,struct()
   r1=setdiff(fieldnames(RunOpt),fieldnames(sel1));
   for j1=1:length(r1); sel1.(r1{j1})=RunOpt.(r1{j1}); end
   RunOpt=sel1;
  else % Unjoin,model,sel1,sel2[,NodeSel]
   sel2=varargin{carg};carg=carg+1;
   RunOpt.type='eltid'; RunOpt.sel1=sel1; RunOpt.sel2=sel2;
  end
 else; sel2=sel1(2); sel1=sel1(1);RunOpt.sel1=sel1; RunOpt.sel2=sel2;
 end
 if carg<=nargin&&ischar(varargin{carg});RunOpt.NodeSel=varargin{carg};carg=carg+1; end

 NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
 [EGroup,nGroup]=getegroup(model.Elt);

 % i2 first group nodes, i3 second group nodes
 switch lower(RunOpt.type)
  case 'group'
   i2=feutil(sprintf('findnode group %i ',RunOpt.sel1),model);
   [ind,elt]=feutil(sprintf('findelt group %i',RunOpt.sel2),model);
  case {'eltid','eltind'}
   if ischar(RunOpt.sel1)
    i2=feutil(sprintf('findnode inelt{%s}',RunOpt.sel1),model);
   else;i2=feutil(sprintf('findnode inelt{%s}',RunOpt.type),model,RunOpt.sel1);
   end
   if ischar(RunOpt.sel2)
    [ind,elt]=feutil(sprintf('findelt %s',RunOpt.sel2),model);
   else;[ind,elt]=feutil(sprintf('findelt %s',RunOpt.type),model,RunOpt.sel2);
   end
  otherwise; error('Unjoin entry for input selection type %s unknown.',RunOpt.type);
 end
 i3=feutil(sprintf('findnode %s',RunOpt.NodeSel),model.Node,elt); % Nodes in second select
 
%  if ischar(sel1)
%    i2=feutil(sprintf('findnode inelt {%s}',sel1),model);
%  else
%    i2=feutil(['findnode group' sprintf('%i ',sel1)],model);
%  end
%  if ~ischar(sel2); sel2=sprintf('group %i',sel2);end
%  [ind,elt]=feutil(sprintf('findelt %s',sel2),model);
%  if carg<=nargin&&ischar(varargin{carg}); st=varargin{carg};carg=carg+1;
%  else; st='groupall';
%  end
%  i3=feutil(['findnode' st],model.Node,elt); % Nodes in second select
  
 [i4,i5]=intersect(i3,i2); % Nodes to duplicate
 [model.Node,i8] = feutil('AddNodeNew',model.Node,model.Node(NNode(i4),:));
 out1=[i4 model.Node(i8,1)];
 NNode=sparse(model.Node(:,1),1,model.Node(:,1));
 NNode(i4)=model.Node(i8,1); % New nodes
 for jGroup=1:nGroup % place new nodes in elements
  ElemF= feutil('getelemf',model.Elt(EGroup(jGroup),:),jGroup);
  cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
  i6=intersect(cEGI,ind);
  if ~isempty(i6)
   i7=fe_super('node',ElemF);
   model.Elt(i6,i7)=reshape(full(NNode(model.Elt(i6,i7))),length(i6),length(i7));
  end
 end
 if isfield(RunOpt,'PostFcn')
  RunOpt.inode=out1;
  RunOpt.el2=elt; 
  model=feval(RunOpt.PostFcn{1},model,RunOpt,RunOpt.PostFcn{2:end});
 end
 out=model;
 
%% #CVS ----------------------------------------------------------------------
elseif comstr(Cam,'cvs')

 out='$Revision: 1.766 $  $Date: 2025/02/19 10:29:01 $';

elseif comstr(Cam,'@'); out=eval(CAM);
 
%% #RMField ------------------------------------------------------------------
elseif comstr(Cam,'rmfield');[CAM,Cam]=comstr(CAM,8); % Remove fields from a structure if they exist

 if sp_util('issdt')
  out=sdtm.rmfield(varargin{2:end});

 else % simple OpenFEM version
  r1=varargin{carg}; carg=carg+1;
  if carg<=nargin; st=varargin{carg}; carg=carg+1; else; out=r1; return; end
  if ~isa(st,'cell')&&~isempty(st); st={st varargin{carg:end}};carg=nargin+1;end
  if ~isstruct(r1);out=r1;return;end;st=intersect(fieldnames(r1),st);
  if ~isempty(st); out=rmfield(r1,st);else; out=r1;end
 end

%% ------------------------------------------------------------------------
elseif exist('feutilb','file'); eval(iigui({'feutilb',nargout},'OutReDir'));

%% ------------------------------------------------------------------------
else;sdtw('%s unknown',CAM);
end
if RSil; Silent=0; end %of_time(-1,Silent,0);

%% #SubFunc

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
function [st,i1,carg]=ineqarg(Cam,carg,varargin);
   if isnumeric(Cam)
    st='=='; i1=double(Cam);
   else
    st=''; Cam=comstr(Cam,1);
    if ~isempty(Cam)&&any(Cam(1)=='~<>='); st=Cam(1);Cam=Cam(2:end);end;
    if ~isempty(Cam)&&any(Cam(1)=='='); st=[st Cam(1)];Cam=Cam(2:end);end;
    i1=comstr(Cam,-1);
   end

   if isempty(i1)&&carg<=nargin-2
    i1=varargin{carg};carg=carg+1;i1=i1(:)';
   elseif isempty(Cam); i1=[];
   elseif ischar(Cam)&&(Cam(1)=='>' || Cam(1)=='<')
    st(2)=Cam(1);i1=comstr(Cam(2:end),-1);
   end

% ------------------------------------------------------------------------
%% #IEqTest integer inequality tests
function [i4]=IEqTest(val,opts,opt);

   if ischar(opt); opt=comstr(opt,-1);end

   if comstr(opts,'~=') 
     %i0=-min(min(opt),min(val))+1; NNode(i0+opt)=1:length(opt);
     %val=val+i0; if max(val)>length(NNode) NNode(max(val),1)=0; end
     %i4=NNode(val);i4=find(i4==0);
     i4=find(~ismember(val,opt));
   elseif comstr(opts,'>='); i4=find(val>=opt(1));
   elseif comstr(opts,'<='); i4=find(val<=opt(1));
   elseif comstr(opts,'>');  i4=find(val>opt(1));
   elseif comstr(opts,'<');  i4=find(val<opt(1));
   else % ==
     %i0=-min(min(opt),min(val))+1;NNode(i0+opt)=1:length(opt);
     %val=val+i0;if max(val)>length(NNode) NNode(max(val),1)=0; end
     %i4=NNode(val); i4=find(i4);
     i4=find(ismember(val,opt));
   end

function [out,out1,CAM]=doNode(model,RO,CAM)
%% #doNode implement node commands with minimal checks
 if nargin<3;CAM='';end
 if isfield(RO,'trans')&&~isempty(RO.trans);
  RB=eye(4);
  RB(:,4)=[RO.trans(:);1];
  RO.trans=[]; RO.rb=RB;
  [out,out1]=doNode(model,RO,CAM); % Call with proper Rigid Body matrix
  %% #NodeRotation
 elseif isfield(RO,'rot')&&~isempty(RO.rot)
  % http://fr.wikipedia.org/wiki/Matrice_de_rotation
  X=reshape(RO.rot(1:3),[],1); % point on the axis of rotation
  N=RO.rot(4:6)/norm(RO.rot(4:6)); % Orientation of the axis of rotation
  theta=RO.rot(7);
  Q=[0 -N(3) N(2);N(3) 0 -N(1);-N(2) N(1) 0];
  R=eye(3)+sin(pi/180*theta)*Q+(1-cos(pi/180*theta))*Q*Q;
  RB=[R -R*X+X;0 0 0 1]; RO.rot=[]; RO.rb=RB;
  [out,out1]=doNode(model,RO,CAM); % Call with proper Rigid Body matrix
  %% #NodeBas
 elseif isfield(RO,'bas')&&RO.bas
  % Parse input for bas commands
  [RO,st,CAM]=paramEdit([ ...
   'oap(#%g#"Origin Axis Plane as node id")' ...
   'lsplane(#3#"Best plane approximation of wireframe")' ...
   ],{RO,CAM}); Cam=lower(CAM);
  if ~isempty(RO.oap) % Origin axis plane definition
   n1=feutil('getnode',model,sprintf('nodeid %i %i %i',RO.oap));
   n1=n1(:,5:7);
   if RO.lsplane % Project the nodes on the best plane approximation
    plane=LSPlane(model.Node(:,5:7)); 
   else % Plane defined by the three points
    plane=LSPlane(n1); 
   end
   % Project the three points on the best plane
   n1=proj2plane(n1,plane);
   trans=n1(1,:);
   n2=n1-trans;
   x=n2(2,:)/norm(n2(2,:)); % xdir in current frame 
   z=plane(1:3)/norm(plane(1:3)); % zdir in current frame
   y=cross(z,x); % ydir in current frame
   % Matrix from expected frame to current frame
   T=[x' y' z' trans';0 0 0 1];
   % Inverse to get matrix from current frame to expected one
   RB=inv(T);

   RO=rmfield(RO,'bas');RO.rb=RB;
   [out,out1]=doNode(model,RO,CAM); % Call with proper Rigid Body matrix
  else
   error('unknown command %s',RO.bas)
  end
  %% #NodeMirror
 elseif isfield(RO,'mir')&&RO.mir
  % Parse input for mir commands
  [RO,st,CAM]=paramEdit([ ...
   'x(#3#"Symmetry plane x=0")' ...
   'y(#3#"Symmetry plane y=0")' ...
   'z(#3#"Symmetry plane z=0")' ...
   'o(#%g#"Point + normal")' ...
   'eq(#%g#"four coeff of plane equation")' ...
   'node(#%g#"Best plane for points list")' ...
   'nodeid(#%g#"Best plane for points list found in model.Node")' ...
   ],{RO,CAM}); Cam=lower(CAM);
  if RO.x % Symmetry plane x=0
   plane=[1 0 0 0];
  elseif RO.y % Symmetry plane y=0
   plane=[0 1 0 0];
  elseif RO.z % Symmetry plane z=0
   plane=[0 0 1 0];
  elseif ~isempty(RO.o) % Symmetry plane passing by (x1,x2,x3) and with normal (n1,n2,n3)
   r1=RO.o;
   plane=[r1(4) r1(5) r1(6) -r1(1)*r1(4)-r1(2)*r1(5)-r1(3)*r1(6)];
  elseif ~isempty(RO.eq); % Symmetry plane with eq aX+bY+cZ+d=0;
   plane=RO.eq(:)';
  elseif ~isempty(RO.node);% Symmetry plane passing by given node list [x y z]
   r1=RO.node;
   plane=LSPlane(r1);
  elseif ~isempty(RO.nodeid); % Symmetry plane passing by given nodeid
   r1=RO.nodeid;
   [un1,r1]=feutil('findnode',model,r1); % points on the mirror
   plane=LSPlane(r1(:,5:7));
  else
   error('unknown command %s',RO.mir)
  end
  % Translate in 4*4 matrix transformation
  % http://fr.wikipedia.org/wiki/Distance_d%27un_point_%C3%A0_un_plan
  % AH = lambda*n
  % mirror = init + 2*AH
  RB=eye(4)-2/norm(plane(1:3))^2*[plane(1:3)'*plane;0 0 0 0];
  RO=rmfield(RO,'mir'); RO.rb=RB;
  [out,out1]=doNode(model,RO,CAM); % Call with proper Rigid Body matrix
  %% #NodeRb (4x4 transformation given) Compute rb if not already provided
 elseif isfield(RO,'rb')
  RB=RO.rb;
  % Test if rigid body
  if round(1000*abs(det(RB)))~=1000
   warning('Transformation does not conserve length')
  end
  if isfield(model,'Node') % model provided
   if strcmpi(RO.sel,'all');node=model.Node(:,1);
   else; node=feutil('findnode',model,RO.sel);
   end
  else % node field or coor provided
   node=model;
  end
  % Do the RB transformation
  if isfield(model,'Node') % model provided, find the selection and apply RB transformation
   ind=ismember(model.Node(:,1),node);
   temp=model.Node(ind,5:7);
   temp=[temp ones(size(temp,1),1)];
   temp=(RB*temp')';
   model.Node(ind,5:7)=temp(:,1:3);
   if isfield(model,'tdof')&&size(model.tdof,2)==5 % Only apply the rotations
    %[un1,ind]=fe_c(model.tdof,node);
    ind=ismember(fix(model.tdof(:,1)),node);
    model.tdof(ind,3:5)=(RB(1:3,1:3)*model.tdof(ind,3:5)')';
   end
  elseif size(model,2)>3 % Node field directly provided
   temp=model(:,5:7);
   temp=[temp ones(size(temp,1),1)];
   temp=(RB*temp')';
   model(:,5:7)=temp(:,1:3);
  elseif size(model,2)==3 % Coor of nodes directly provided
   temp=model;
   temp=[temp ones(size(temp,1),1)];
   temp=(RB*temp')';
   model=temp(:,1:3);
  end
  % Warning to do
  if det(RB)<0&&~RO.Silent
   sdtw('_nb','Use feutil(''Orient'') to flip nodes in elt + scan FaceId set in the stack');
  end 
  out=model;out1=RB;
  %% #NodeDefShift Apply def field to Nodes
 elseif isfield(RO,'DefShift')
  if isnumeric(RO.DefShift)&&isfield(RO,'def');def=RO;
  elseif isnumeric(RO.DefShift)&&isfield(RO,'T');
    % called from IterInitLagUp
     r2=fe_c(RO.DOF(:,1), ...
      [model.Node(:,1)+.01;model.Node(:,1)+.02;model.Node(:,1)+.03],'place');  
     out=r2*RO.T; return;
  else; def=RO.DefShift;
  end
  i1=fe_c(def.DOF,0.01,'ind'); %dof ind in direction x
  i2=fe_c(def.DOF,0.02,'ind'); %dof ind in direction y
  i3=fe_c(def.DOF,0.03,'ind'); %dof ind in direction z
  [un1,ind]=ismember(fix(def.DOF),model.Node(:,1)); % Row indice of each dof in model.Node
  ind=sub2ind(size(model.Node),...
   [ind(i1);ind(i2);ind(i3)],... % rowSub : node number
   [repmat(5,length(i1),1);repmat(6,length(i2),1);repmat(7,length(i3),1)]); % colSub : x y or z
  if ~isfield(RO,'coeff'); RO.coeff=1; end
  model.Node(ind)=model.Node(ind)+RO.coeff*def.def([i1;i2;i3]);
  out=model; out1=[];
 end
 
%% #get_nodeelt  clean import of model arguments --------------------------
% WHEN CHANGING CHANGE IN FEUTILB
function [carg,node,elt,el0,ModelStack,omodel]=get_nodeelt(varg,carg,ModelStack,RunOpt);

  node=varg{carg};carg=carg+1;ModelStack={};model=[];
  if nargin<4;RunOpt=[];end
  if nargout==6;omodel=[];end
  if isstruct(node)&&isfield(node,'Node') % This is a model

   model=node; 
   if isfield(model,'mdl'); 
    model=model.mdl; 
    if nargout==6; omodel=model;else;assignin('caller','model',model);end
   elseif isfield(model,'Elt'); 
    if nargout==6; omodel=model;else;assignin('caller','model',model);end
   end
   if isfield(model,'Stack'); ModelStack=model.Stack;end
   if ~isfield(node,'Elt'); error('Fields .Node & .Elt must be defined');end
   elt=model.Elt; 
   if isfield(model,'El0'); el0=model.El0; else;el0=[];end
   
   if ~isempty(ModelStack)&&any(strcmp(ModelStack(:,1),'SE'))
      if isfield(model,'CheckedNodes')&&model.CheckedNodes(1)==1
         node=model.Node;
         evalin('caller','RunOpt.CheckedNodes=ones(2,1);');
      elseif isfield(RunOpt,'ImplicitNode'); 
          node='[node,bas]=basis(''nodebas'',model);';
      else; 
          [node,bas]=basis('nodebas',model);
      end
   elseif isfield(model,'bas')&&~isempty(model.Node)&&~isempty(model.bas) ...
       &&any(model.Node(:,2))
      if isfield(RunOpt,'ImplicitNode'); 
          node='[node,bas]=basis(''nodebas'',model);';
      else; [node,bas]=basis('nodebas',model);
      end
   else; node=model.Node;
   end
  elseif ~isempty(node)&&isnumeric(node)&&~isfinite(node(1)) % elt only
   elt=node; node=[]; el0=[];
  elseif ~isempty(node)&&size(node,2)~=7
    error('Node matrices must have 7 columns');
  else
   elt = varg{carg};carg=carg+1;
   if nargout<4
   elseif carg<=length(varg)&&(isempty(varg{carg})||~isfinite(varg{carg}(1)))
      el0 = varg{carg};carg=carg+1;
   else; el0=[]; 
   end
  end

% ------------------------------------------------------------------------
% return non repeated lines in a list
function [list,prop]=single_line(list,prop);

% eliminate degenerate lines
if isempty(list);return;end
i1=find(list(:,1)==list(:,2));list(i1,:)=[];prop(i1,:)=[];

% flip lines to have lowest number node first

% fill in zeros
for j1=size(list,2)-1:-1:2 
 i1=find(~any(list(:,j1+1:end),2)&list(:,1)>list(:,j1));
 if ~isempty(i1); list(i1,1:j1)=list(i1,j1:-1:1);end
end

i3=1;
while ~isempty(i3)
   [r2,i1]=sortrows(sort(list,2));
   i2=1:size(list,1);i3=find(~any(diff(r2),2));
   i2(i3)=0;i2(i3+1)=0;i2=find(i2);
   list=list(i1(i2),:);prop=prop(i1(i2),:);
end


% ------------------------------------------------------------------------
% #BuildSelStack ---------------------------------------------------------
function [Stack,carg]=BuildSelStack(CAM,node,elt,el0,carg,varargin)
persistent nmap inCall
if isempty(inCall); inCall=0; end % recursive call counter to keep nmap variable through the call stack
inCall=inCall+1;
Args=varargin; if inCall==0; nmap=[]; end
if ~isempty(Args)&&isa(Args{end},'vhandle.nmap'); nmap=Args{end}; Args(end)=[];end

% end of edge selection - - - - - - - - - - - - - - - - - - - - - - - - - -

if isempty(CAM)&&~isempty(Args)&&ischar(Args{1});
 CAM=Args{1};carg=carg+1;
 [epsl,CAM,Cam]=test_epsl(-1,CAM,lower(CAM));
 if epsl~=-1;assignin('caller','epsl',epsl);end
end
i7=[]; CAM=[CAM ' '];out=[]; opt=[];
%[EGroup,nGroup]=getegroup(elt);

in1=strfind(CAM,'"'); % skip specical characters in double quoted strings
if rem(length(in1),2); error('double quotes " must come in pair'); end
ind=[find(CAM=='&'|CAM=='|'|CAM=='{') length(CAM)+1]; 
for j1=1:2:length(in1)-1; ind(ind>in1(j1)&ind<in1(j1+1))=[];end
j1=0;Stack=cell(10,4);jend=0;

% ind keeps track of what has been read
while j1<=length(CAM)  % loop on the selected arguments to build stack
try;
 [st,Cam]=comstr(CAM(j1+1:ind(1)-1),1); 
 if ind(1)<length(CAM)&&CAM(ind(1))=='{'
  [st1,CAM(ind(1):end)]=comstr(CAM(ind(1):end),-8,'{}');
  ind = ind(ind>=ind(1)+length(st1));
 else;st1='';
 end
 if j1; jend=jend+1;Stack{jend,1}=CAM(j1); else;Stack{1}='{'; jend=1; end
 if comstr(Cam,'~'); 
  if strcmp(Stack{jend,1},'&'); Stack{jend,1}=[Stack{jend,1} '~']; 
  else; error('Invalid boolean %s~',Stack{jend,1});
  end
  [st,Cam]=comstr(CAM(j1+1:ind(1)-1),2); 
 end
 % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if comstr(Cam,'between')  
   [opt,st,Cam]=comstr(Cam,'between','%i');
   if isempty(opt)&&carg<=length(Args) && ~ischar(Args{carg})
    opt=Args{carg}; carg=carg+1;
   end
   Stack(jend,2:4)={'between',[],opt};
 % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'egid');  [st,Cam]=comstr(st,5);
   [opts,opt,carg]=ineqarg(Cam,carg,Args{:});i1 = 5;
   if isempty(opt)&&carg<=length(Args) && ~ischar(Args{carg})
    opt=Args{carg}; carg=carg+1;
   end
   Stack(jend,2:4)={'EGID',opts,opt};
 % EltName - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'eltname')
   [opt,Cam]=comstr(st,8); i1 = 7; opts='';
   if comstr(opt,'~='); opts='~='; opt=comstr(opt,3);
   elseif comstr(opt,'==');opts=''; opt=comstr(opt,3);
   end
   Stack(jend,2:4)={'EltName',opts,opt};

 % EltId - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'eltid'); [st,Cam]=comstr(Cam,'eltid','%c');
   [opts,opt,carg]=ineqarg(Cam,carg,Args{:});i1 = 9;
   if isempty(opt)&&carg<=length(Args) && ~ischar(Args{carg})
    opt=Args{carg}; carg=carg+1;
   end
   Stack(jend,2:4)={'EltId',opts,opt};
 % EltInd - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'eltind'); [st,Cam]=comstr(Cam,'eltind','%c');
   [opts,opt,carg]=ineqarg(Cam,carg,Args{:});i1 = 8;
   if isempty(opt)&&carg<=length(Args) && ~ischar(Args{carg})
    opt=Args{carg}; carg=carg+1;
   end
   Stack(jend,2:4)={'EltInd',opts,opt};

 % Facing- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'facing'); [st,Cam]=comstr(Cam,'facing','%c');
   [opts,opt,carg]=ineqarg(Cam,carg,Args{:});
   if isempty(opt)&&carg<=length(Args) && ~ischar(Args{carg})
    opt=Args{carg}; carg=carg+1;
   end
   Stack(jend,2:4)={'Facing',opts,opt};

 % GID - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'gid'); [st,Cam]=comstr(st,4);

   st1='GID'; 
   [opts,opt,carg]=ineqarg(Cam,carg,Args{:});
   if isempty(opt)&&carg<=length(Args) && ~ischar(Args{carg})
    opt=Args{carg}; carg=carg+1;
   end
   Stack(jend,2:4)={st1,opts,opt};

 % Group - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'group'); [st,Cam]=comstr(st,6);
  [EGroup,nGroup]=getegroup(elt);
   st1='Group'; 
   if comstr(Cam,'al'); opt=1:nGroup; opts='==';
   elseif comstr(Cam,'a') 
      st1='GroupA'; [opts,opt,carg]=ineqarg(Cam(2:end),carg,Args{:});
   else;[opts,opt,carg]=ineqarg(Cam,carg,Args{:});
   end
   if isempty(opt)&&carg<=length(Args) && ~ischar(Args{carg})
    opt=Args{carg}; carg=carg+1;
   end
   Stack(jend,2:4)={st1,opts,opt};

 elseif comstr(Cam,'in');  [st,Cam]=comstr(st,3);% - - - - - - - - - - - -
 % InNode - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if comstr(Cam,'node')

   [st,Cam]=comstr(st,5); opt = round(comstr(st,[-1]));
   if     isempty(opt)&&isempty(st1)&&carg<=length(Args) && ~ischar(Args{carg})
     opt=Args{carg}; carg=carg+1;   Stack(jend,2:4)={'InNode','',opt};
   elseif isempty(opt)&&(~isempty(st1)||~isempty(st))
    if isempty(st1); st1=st; else;st1=st1(2:end-1); end
    %opt=feutil(['findnode' st1(2:end-1)],node,elt,el0);
    Stack(jend,2:4)={'InNode','{}',''};
    [Stack2,carg]=BuildSelStack(st1,node,elt,el0,carg,Args{:});
    Stack(jend+[1:size(Stack2,1)],1:4)=Stack2;jend=jend+size(Stack2,1);
   elseif ~isempty(opt)&&isempty(st1)
    Stack(jend,2:4)={'InNode','',opt};
   elseif ~isempty(opt)&&~isempty(st1)
    error('InNode : Recursion is only valid with string argument');
   else;warning('Incomplete InNode Selection');
   end

 % 'InElt' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  elseif comstr(Cam,'elt')

   if isempty(st1); error('InElt requires specification of a selection');end
   if  Cam(end)=='t'; Stack(jend,2:3)={'inelt','{}'};
   elseif  Cam(end)=='0'; Stack(jend,2:3)={'inel0','{}'};
   else;error('Not a valid FindNode InElt Command');
   end
   st2=st1(2:end-1); 
   if isempty(strtrim(st2)); st2='groupall'; % inf loop otherwise
    if ~feutil('Silent'); sdtw('_nb','Empty selection to inElt will use groupall'); end
   end
   [Stack2,carg]=BuildSelStack(st2,node,elt,el0,carg,Args{:});
   Stack(jend+[1:size(Stack2,1)],1:4)=Stack2;jend=jend+size(Stack2,1);
  else
    error('in %s : improper feutil selection command',Cam);
  end

 % MatId  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'mat');  [st1,Cam]=comstr(Cam,'matid','%c');
  if comstr(Cam,'name')
   [st1,Cam]=comstr(st,'matname','%c');
   try; Cam=nmap.('Map:MatName').(st1); catch; error('%s is not a MatName',st1); end
  end
  [opts,opt,carg]=ineqarg(Cam,carg,Args{:});i1 = 4;
  if isempty(opt)&&carg<=length(Args) && ~ischar(Args{carg})
   opt=Args{carg}; carg=carg+1;
  end
  Stack(jend,2:4)={'MatId',opts,opt};

 % NodeId - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'nodeid'); [st,Cam]=comstr(st,7);

   [st,opt,carg]=ineqarg(Cam,carg,Args{:});
   Stack(jend,2:4)={'nodeid',st,opt};

 % 'NotIn' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  elseif comstr(Cam,'notin')

   if isempty(st1); error('NotIn requires specification of a selection');end
   Stack(jend,2:3)={'notin','{}'};
   [Stack2,carg]=BuildSelStack(st1(2:end-1),node,elt,el0,carg,Args{:});
   Stack(jend+[1:size(Stack2,1)],1:4)=Stack2;jend=jend+size(Stack2,1);

 % Plane - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'plane'); [st,Cam]=comstr(st,6);

  st1='==';if any(Cam(1)=='~<>='); st1=Cam(1);Cam=Cam(2:end);end;
  if ~isempty(Cam)&&any(Cam(1)=='='); st1=[st1 Cam(1)];Cam=Cam(2:end);end;
  [st,Cam]=comstr(Cam,1);
  if comstr(Cam,'o');    opt=comstr(Cam(2:length(Cam)),[-1 0 0 0 0 0 0]);
  elseif isempty(Cam)&&carg<=length(Args);  opt=Args{carg};carg=carg+1;
  else;   opt =  comstr(Cam,[-1 0 0 0 0]);i1=find(node(:,1)==opt(1));
          if isempty(i1); i1=[0 0 0];else;i1=node(i1(1),[5:7]); end
          opt=[i1 opt(2:4)];
  end

  Stack(jend,2:4)={'plane',st1,opt};
 % cyl==r o X Y Z nx ny nz - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'cyl'); [st,Cam]=comstr(st,4);

  st1='==';if any(Cam(1)=='~<>='); st1=Cam(1);Cam=Cam(2:end);end;
  if ~isempty(Cam)&&any(Cam(1)=='='); st1=[st1 Cam(1)];Cam=Cam(2:end);end;
  [st,Cam]=comstr(Cam,1);
  [st,Cam,i1]=comstr('o',[-25 2 3],st,Cam);
  if ~isempty(i1);opt=comstr(st,[-1 0 0 0 0]);opt=[opt(1) i1(:)' opt(2:end)];
  elseif isempty(Cam)&&carg<=length(Args);  opt=Args{carg};carg=carg+1;
  else;   opt =  comstr(Cam,[-1 0 0 0 0 0]);i1=find(node(:,1)==opt(2));
          if isempty(i1); i1=[0 0 0];else;i1=node(i1(1),[5:7]); end
          opt=[opt(1) i1 opt(3:5)];
  end
  Stack(jend,2:4)={'cyl',st1,opt};

 % ProId  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'pro'); [st1,Cam]=comstr(Cam,'proid','%c');
  if comstr(Cam,'name')
   [st1,Cam]=comstr(st,'proname','%c');
   try; Cam=nmap('Map:ProName').(st1); catch; error('%s is not a ProName',st1); end
  end
  [opts,opt,carg]=ineqarg(Cam,carg,Args{:});i1 = 3;
  if isempty(opt)&&carg<=length(Args) && ~ischar(Args{carg})
   opt=Args{carg}; carg=carg+1;
  end
  Stack(jend,2:4)={'ProId',opts,opt};
 % Rad - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'rad'); [st,Cam]=comstr(st,4); st='';

   if any(Cam(1)=='~<>='); st=Cam(1);Cam=Cam(2:end);end;
   if ~isempty(Cam)&&any(Cam(1)=='='); st=[st Cam(1)];Cam=Cam(2:end);end;

   if isempty(Cam)&&carg<=length(Args)
    i1=Args{carg};carg=carg+1;i1=i1(:)'; i2 = length(i1);
   else;[i1,i2,i3,i6]=sscanf(Cam,'%g %g %g %g');i1=i1'; 
   end
   Stack(jend,2:4)={'rad',st,i1};

 % SelEdge - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'seledge')
  [st,Cam]=comstr(Cam,'seledge','%c');
  Stack(jend,2:4)={'SelEdge',[],st};
 % SelFace  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'selface')||comstr(Cam,'safeselface') 
  if comstr(Cam,'safe'); sts='SelFaceF'; stc='safeselface';
  else; sts='SelFace'; stc='selface';
  end
  [st,Cam]=comstr(Cam,stc,'%c');
  Stack(jend,2:4)={sts,[],st};
 % SelName  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'selname'); 

   Cam=st;Cam(1:7)='';[st,i1,i2,i3]=sscanf(Cam,'%s',1);Cam=Cam(i3:end);
   if st(1)=='"'; st(1)=[];% with blanks
    i2=min(find(st=='"'));
    if ~isempty(i2)&&i2>1;  % cannot allow empty name with quotes
     st=st(1:i2-1);Cam=comstr(horzcat(st(i2+1:end),Cam),1);
     if ~isempty(Cam)&&Cam(1)=='"';Cam=comstr(Cam(2:end),1);end
    else
     i2=strfind(Cam,'"'); %i2=min(find(Cam=='"'));
     if length(i2)>1; % allow double quotes in setnames extremities
      i3=find(diff(i2)>1,1,'first'); if isempty(i3); i3=length(i2); end
      i2=i2(i3);
     end
     st=horzcat(st(1:end),Cam(1:i2-1));
     Cam=comstr(Cam(i2+1:end),1);
    end
   end
   Stack(jend,2:4)={'Sel',[],st};
 % SetName  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'setname')||comstr(Cam,'safesetname') 

  if comstr(Cam,'safe'); sts='SetF'; st=comstr(st,5); %[CAM,Cam]=comstr(CAM,5);
  else; sts='Set'; 
  end
   Cam=st;Cam(1:7)='';[st,i1,i2,i3]=sscanf(Cam,'%s',1);Cam=Cam(i3:end);
   if st(1)=='"'; st(1)=[];% with blanks
    i2=find(st=='"',1,'first'); st0=comstr(st(i2+1:end),1);
    if ~isempty(i2)&&i2>1;  % cannot allow empty name with quotes
     st=st(1:i2-1);Cam=comstr(horzcat(st(i2+1:end),Cam),1);
     if ~isempty(Cam)&&Cam(1)=='"';Cam=comstr(Cam(2:end),1);end
    else
     i2=strfind(Cam,'"'); %i2=min(find(Cam=='"')); 
     if length(i2)>1; % allow double quotes in setnames extremities
      i3=find(diff(i2)>1,1,'first'); if isempty(i3); i3=length(i2); end
      i2=i2(i3);
     end
     st=horzcat(st(1:end),Cam(1:i2-1));
     Cam=comstr(Cam(i2+1:end),1);
    end
    if ~isempty(st0)&&comstr(st0,':')
     [u1,u2,st0]=comstr(':',[-25 4 1],st0,lower(st0));
     st=[st ':' st0];
    end
   elseif any(strfind(st,':con'))||any(strfind(st,':f'))||any(strfind(st,':exclude'))
    st=[st Cam]; % allow  blanks if ":" present as can provide options
   end
   Stack(jend,2:4)={sts,[],st}; %'Set'

 % Set  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'set'); [st,Cam]=comstr(Cam,'set','%c'); sts='';
  i1=strfind(Cam,':'); if ~isempty(i1); sts=Cam(i1(1):end); Cam=Cam(1:i1(1)-1); end
  [opts,opt,carg]=ineqarg(Cam,carg,Args{:});i1 = 3;
  if isempty(opt)&&carg<=length(Args) && ~ischar(Args{carg})
   opt=Args{carg}; carg=carg+1;
  end
  if isempty(opts)&&~isempty(sts); opts=sts; end
  Stack(jend,2:4)={'Set',opts,opt};

 % WithoutNode - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'without');[st,Cam]=comstr(st,8); % - - - - - - - - - -
  if ~comstr(Cam,'node')
    error('without %s : improper feutil(''findelt'') command',Cam);
  end
  [st,Cam]=comstr(st,5); opt = round(comstr(st,[-1]));   
  if     isempty(opt)&&isempty(st1)&&carg<=length(Args) && ~ischar(Args{carg})
    opt=Args{carg}; carg=carg+1;Stack(jend,2:4)={'WithoutNode','',opt};
  elseif isempty(opt)&&~isempty(st1)
   Stack(jend,2:4)={'WithoutNode','{}',''};
   [Stack2,carg]=BuildSelStack(st1(2:end-1),node,elt,el0,carg,Args{:});
   Stack(jend+[1:size(Stack2,1)],1:4)=Stack2;jend=jend+size(Stack2,1);
  elseif ~isempty(opt)&&isempty(st1)
   Stack(jend,2:4)={'WithoutNode','',opt};
  elseif ~isempty(opt)&&~isempty(st1)
   error('WithoutNode : Recursion is only valid with string argument');
  end

 elseif comstr(Cam,'$');  [st,Cam]=comstr(st,1); % - - - - - - - - - -
   Stack(jend,2:4)={'$','',st};
 elseif comstr(Cam,'with');  [st,Cam]=comstr(st,5); % - - - - - - - - - -
  if ~comstr(Cam,'node')
    error('with %s : improper feutil(''findelt'') command',Cam);
  end
  [st,Cam]=comstr(st,5);
  if isscalar(st)&&isequal(st,'(');opt=[];
  else;  opt = round(comstr(st,-1)); 
  end
  if     isempty(opt)&&isempty(st1)&&carg<=length(Args) && ~ischar(Args{carg})
    opt=Args{carg}; carg=carg+1;    Stack(jend,2:4)={'WithNode','',opt};
  elseif isempty(opt)&&~isempty(st1)
   %opt=feutil(['findnode' st1(2:end-1)],node,elt,el0);
   Stack(jend,2:4)={'WithNode','{}',''};
   [Stack2,carg]=BuildSelStack(st1(2:end-1),node,elt,el0,carg,Args{:});
   Stack(jend+[1:size(Stack2,1)],1:4)=Stack2;jend=jend+size(Stack2,1);
  elseif ~isempty(opt)&&isempty(st1)
   Stack(jend,2:4)={'WithNode','',opt};
  elseif ~isempty(opt)&&~isempty(st1)
   error('WithNode : Recursion is only valid with string argument');
  end

 elseif comstr(Cam,'connectedto');  [st,Cam]=comstr(st,12); % - - - - - - - - - -
  opt = round(comstr(st,[-1])); 
  Stack(jend,2:4)={'ConnectedTo','',opt};

 % xyzr - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif ~isempty(Cam) && any(Cam(1)=='xyzr') % FindNode xyz

   i2=find('xyzr'==Cam(1)); st1='xyzr';   [st,Cam]=comstr(st,2); st='';
   if any(Cam(1)=='~<>='); st=Cam(1);Cam=Cam(2:end);end;
   if ~isempty(Cam)&&any(Cam(1)=='='); st=[st Cam(1)];Cam=Cam(2:end);end;
   i1=comstr(Cam,[-1]); 

   if isempty(i1)&&~isempty(Cam); i1=Cam;
   elseif isempty(i1)&&carg<=length(Args);  i1=Args{carg};carg=carg+1;
   end
   if isempty(i1); error('xyzr ><= ... must specify value'); end
   if ~ischar(i1); i1=i1(1);end 
   Stack(jend,2:4)={st1(i2),st,i1}; %#ok<FNDSB>
 % name - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif strncmpi(Cam,'name',4) % name using nmap
    Stack(jend,2:4)={'name','',comstr(st(5:end),1)};
 % Dist{ using lsutil distFcn  - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'distfcn'); 
   if length(Cam)==7;st=Args{carg-1};
   else; st=CAM;
   end
   Stack(jend,2:4)={'distfcn',[],st};
   
 % Dist - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'dist'); 
   r2=regexpi(st,'[Dd]ist\(([^)]*))([\w><=].*)','tokens');r2=r2{1};
   Stack(jend,2:4)={'dist',r2{1},r2{2}};
 % give coordinates  with possible infinity - - - - - - - - - - - - - - - -
 elseif ~isempty(Cam)&&~all(Cam==' ') 

   if isempty(st)&&carg<=length(Args); opt=Args{carg};carg=carg+1;
     Stack(jend,[2 4])={'Position',opt};

   elseif ~isempty(st) &&any(comstr(st(1),-27)=='0123456789.I+-') 
     opt=comstr(st,[-1]);
     if length(opt)==3||any(rem(opt,1)~=0);
       opt(end+1:3)=Inf;Stack(jend,[2 4])={'Position',opt};
     else
       Stack(jend,[2 3 4])={'NodeId','==',opt(:)};
     end
   else; error('''%s'' is not a valid FindNode/FindElt command',st);
   end
   

 end

catch
 inCall=inCall-1;
  if j1<2
   error('''%s ...'' is not a valid selection format', ...
     CAM(1:min(j1+20,size(CAM,2))));
  else
   error(''' ... %s'' is not a valid selection format', ...
     CAM(j1:min(j1+20,size(CAM,2))));
  end
 end

 j1=ind(1);ind=ind(2:end); 

end % of j1 loop on arguments
jend=jend+1;Stack{jend,1}='}';jend=jend+1;
if jend<=size(Stack,1);Stack(jend:end,:)=[];end
inCall=inCall-1;
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
function ind=SubStackIndex(Stack,i0);

   i1=find(strcmp(Stack(:,1),'{'));i2=find(strcmp(Stack(:,1),'}'));

   i4=[1 i0]; 
   while i4(1)~=0 
       i3=min(i1(i1>i4(2)));
       i5=min(i2(i2>i4(2)));
       if i3<i5 % opening first
         i4(1)=i4(1)+1; i4(2)=i3;
       else % i5 is closing
         i4(1)=i4(1)-1; i4(2)=i5;
       end
   end
   ind=i0:i4(2);

% ------------------------------------------------------------------------
%% #test_epsl ----------------------------------------------------------------
function [epsl,CAM,Cam]=test_epsl(epsl,CAM,Cam);

if isempty(epsl);try;epsl=sp_util('epsl');catch;epsl=sdtdef('epsl'); end;end 
i1=strfind(Cam,'epsl'); 
if ~isempty(i1);
  [epsl,i2,i3,i4]=sscanf(Cam(i1+4:end),'%g');
  CAM(i1+[0:i4+2])=''; [CAM,Cam] = comstr(CAM,1);
end

%% #upEG: update EGroup, nGroup if reset
function [EGroup,nGroup,elt]=upEG(EGroup,nGroup,elt,j1,out);
if isempty(EGroup)||isempty(nGroup); [EGroup,nGroup]=getegroup(elt); end
if nargin>3
 if j1>1  % keep current elements
  i7=zeros(nGroup,1);
  for jGroup = 1:nGroup
   if any(out>EGroup(jGroup)&out<EGroup(jGroup+1));i7(jGroup)=1;end
  end
  i7=sort([out;EGroup(i7~=0)']); elt=elt(i7,:);
 end
end
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%% #FeutilMatchSet. For sets see sdtweb('model'), one supports SelN and Set
function [i4,elt]=FeutilMatchSet(ModelStack,Stack,elt,RO);

j1=1;i4=[];st1='';
if isnumeric(Stack{j1,4}) % Selection by ID number in SelV

 % Build the list of ID values
 i1=zeros(size(ModelStack,1),2);
 for j2=1:size(ModelStack,1)
   r1=ModelStack{j2,3};
   try; if isfield(r1,'ID')&&~isempty(r1.ID);i1(j2,1)=r1.ID(1);end;end
   i2=find(strcmpi(ModelStack{j2,1},{'set','seln','sele'}));
   if ~isempty(i2); i1(j2,2)=i2(1);end
 end
 i2=find(i1(:,2)==1&i1(:,1));
 if ~isempty(Stack{j1,3})&&ischar(Stack{j1,3})&&Stack{j1,3}(1)==':';
  st1=Stack{j1,3}(2:end); Stack{j1,3}=[]; % xxx combination of options ?
 end
 opt=IEqTest(i1(i2,1),Stack{j1,3},Stack{j1,4});opt=i2(opt);

elseif ischar(Stack{j1,4});
 st=Stack{j1,4};
 if any(st==':');
     opt=find(st==':');st1=comstr(st(opt(1)+1:end),1);st=st(1:opt(1)-1);
 end
 if ~isempty(st)&&st(1)=='#';
    opt=find(~cellfun('isempty',regexp(ModelStack(:,2),st(2:end))));
 else;
    opt=find(strcmpi(st,ModelStack(:,2)));
 end
 if isempty(opt)%&&any(strcmpi(Stack{j1,4},ModelStack(:,2))); % misleading ":" token
  r5=cellfun(@(x)strncmpi(Stack{j1,4},x,length(x)),ModelStack(:,2));
  if any(r5) % check if we have external : legit tokens
   st5=ModelStack(r5,2); i5=cellfun(@length,st5)<length(Stack{j1,4});
   st6=st5; st6(:)={''};
   st6(i5)=cellfun(@(x)Stack{j1,4}(x+1),cellfun(@length,st5(i5),'uni',0),'uni',0);
   opt=find(r5); 
   i6=cellfun(@isempty,st6); 
   if any(i6); opt=opt(i6); st1='';
   else
    i6=ismember(st6,':');opt=opt(i6);
    if sum(i6)==1;st=ModelStack{opt(1),2}; st1=comstr(Stack{j1,4},length(st)+2);end
   end
   %opt=find(strcmpi(Stack{j1,4},ModelStack(:,2)));
  end
 end
 if isempty(opt); error('Name %s is not a set',Stack{j1,4})
 elseif length(opt)>1;
  opt=opt(ismember(ModelStack(opt,1),{'set','seln','sele'})); % robust to types
  if length(opt)>1; fprintf('Multiple match to %s\n',st); end
 end
 if isempty(opt)||~any(ismember(ModelStack(opt,1),{'set','seln','sele'})) %~isempty(setdiff(ModelStack(opt,1),{'set','seln','sele'}))
  error('Name is not a set, seln or sele');
 end
elseif isstruct(Stack{j1,4})
  opt=Stack{j1,4};
end
if isempty(opt); error('Set selection failed');
elseif nargin==2 % Match a set pointing to nodes

 i4=[];model=[];
 for j2=1:length(opt)
  r1=ModelStack{opt(j2),3};
  if isfield(r1,'SConn') % resolve mate-set based selection
   if ~isempty(st1) % NodeId sets with sub names % 'setname _PSE:PSE_1', setname _PSE:PSE_1:PSE_2
    st2=textscan(st1,'%s','delimiter',':'); st2=st2{1};
    i1=ismember(lower(r1.SetNames(:,1)),lower(st2));
    if ~any(i1);i1=str2double(st1);% Allow numeric index _PSE:2
     if ~isfinite(i1);error('%s not valid',st1);end
    end
   else; i1=ismember(r1.type,'NodeId'); % all
   end
   if isempty(i1)||~any(i1); error('No elements in set %s',Stack{j1,4}); end
   r1.type=r1.type(i1);  [r3,i3]=unique(r1.type);
   if length(r3)>1
    sdtw('_nb','selection %s in set %s yields mutliple set types, set to NodeId',st1,Stack{j1,4})
    r1.type='NodeId';
   else; r1.type=r3{1};
    if ~strcmpi(r1.type,'nodeid'); error('set %s does not refer to NodeId type',st1); end
   end
   r1.data=r1.NodeId(logical(sum(r1.NConn(:,i1),2)));
  end
  
  st2='';if isfield(r1,'type');st2=r1.type; end
  if iscell(st2)&&isscalar(unique(st2));st2=st2{1};
  elseif ischar(st2); % ok (old format)
  else;error('Not implemented');
  end
  if ~isfield(r1,'type')
   if isfield(r1,'data');r1=r1.data(:,1);
   elseif ~isfinite(r1(1))&&size(r1,2)>3% Elements (possible for master of contact)
     r1=feutil('findNodegroupall',r1);
   end
  elseif strcmpi(st2,'FaceId')||strcmpi(st2,'EltId')
   if isempty(model);model=evalin('caller','model');end
   st2=ModelStack{opt(j2),2};
   if ~isempty(st1);st2=[st2 ':' st1];end %#ok<AGROW>
   r1=feutil(sprintf('findnode inelt{ setname "%s"}',st2),model);
  elseif strcmpi(r1.type,'NodeId'); r1=r1.data;
  end
  if ischar(r1); r1=feutil('findnode',evalin('caller','model'),r1);end
  if ~isempty(st1)
   if comstr(st1,'exclude')
    model=evalin('caller','model'); r1=setdiff(model.Node(:,1),r1);
   end
  end
  i4=[i4;r1];
 end

elseif ~isstruct(opt)&&isscalar(opt)&&ischar(ModelStack{opt,3})% char selection

  node=evalin('caller','node');
  if comstr(ModelStack(opt,1),'seln'); % THIS DOES NOT WORK ?
   i4=feutil(horzcat('findnode',ModelStack{opt,3}),node,elt);
  else % 
   i4=feutil(horzcat('selelt',ModelStack{opt,3}),node,elt);
  end

else % standard cases with matching - - - - - - - - - - -

 if nargin>3&&isfield(RO,'EltIdFixed'); else; RO=struct('EltIdFixed',-1); end
 RO.MultiCall=0;
 if ~isfield(RO,'Transformed'); RO.Transformed=0; end
 %RO.EltId=feutil('eltid;',elt); RO.EltIdFixed=1; RO.Transformed=0; xxxforce
 if RO.EltIdFixed>0&&~RO.Transformed; %EltId=feutil('EltIdSkipCheck;',elt);
  EltId=RO.EltId; RO.MultiCall=1;
 else % base case do fix, but try to communicate to do it once if possible
  EltId=feutil('eltidfix;',elt);assignin('caller','EltId',EltId);
  if ~RO.EltIdFixed; RO.EltIdFixed=1; RO.EltId=EltId; assignin('caller','RunOpt',RO); end
 end
 %[EGroup,nGroup]=getegroup(elt);
 i4=[];
 if isstruct(opt);r2=opt;
 elseif length(opt)>1;
     r2=ModelStack{opt(1),3};
     r2.data=cellfun(@(x)x.data,ModelStack(opt,3),'UniformOutput',false);
     if ~all(cellfun(@isnumeric,r2.data))
      r2.data=sprintf('setname"%s" |',ModelStack{opt,2}); r2.data(end)=''; 
     elseif length(unique(cellfun(@(x)size(x,2),r2.data)))>1
      sdtw('_nb','multiple data types in matched sets, keeping type %s of first matched one %s',...
       r2.type,ModelStack{opt(1),2})
      r2.data=r2.data(ismember(cellfun(@(x)x.type,ModelStack(opt,3),'uni',0),r2.type));
      r2.data=vertcat(r2.data{:});
     else;     r2.data=vertcat(r2.data{:});
     end
 else;r2=ModelStack{opt,3};
 end
 if isfield(r2,'data')&&isa(r2.data,'v_handle'); r2=sdthdf('hdfreadref',r2); end
 st1=lower(st1);
 if isfield(r2,'SConn') % New set format SConn, NConn, SetNames
  if ~isempty(st1) % EltId sets with sub names
   if comstr(st1,'con') % use node connectivity to get corresponding sets
    st1=comstr(st1,4); % setname face:con 65 // setname face:con~=65
    if comstr(st1,'~='); i3=1; st1=comstr(st1,3); else; i3=0; end
    i2=comstr(st1,-1);if isempty(i2);error('No connected node');end
    if isfield(r2,'NNode'); NNode=r2.NNode; 
    else; NNode=sparse(r2.NodeId,1,1:length(r2.NodeId)); % xxx perf.
    end
    Nc=r2.NConn'; i1=full(logical(sum(Nc(:,NNode(i2)),2)));
    if i3; i1=~i1; end % revert if ~=
   else; % 'setname _PSE:PSE_1', setname _PSE:PSE_1:PSE_2
    st2=textscan(st1,'%s','delimiter',':'); st2=st2{1};
    for j5=1:length(st2) % allow double quotes to escape tokens
     if isequal(st2{j5}(1),'"')&&isequal(st2{j5}(end),'"'); st2{j5}([1 end])=''; end; 
    end
    i1=ismember(lower(r2.SetNames(:,1)),st2);
    if ~any(i1);i1=str2double(st1);% Allow numeric index _PSE:2
     if ~isfinite(i1);error('%s not valid',st1);end
    end
   end
  else; i1=~ismember(r2.type,'NodeId'); % all
  end
  if isempty(i1)||~any(i1); error('No elements in set %s',Stack{j1,4}); end
  r2.type=r2.type(i1);  [r3,i3]=unique(r2.type);
  if length(r3)>1
   sdtw('_nb','selection %s in set %s yields mutliple set types, set to EltId',st1,Stack{j1,4})
   r2.type='EltId';
  else; r2.type=r3{1};
  end
  
  if strncmpi(r2.type,'FaceId',3);
   r3=r2.SConn(:,i1);[i3,j3,k3]=find(r3); % here will loop on face types
   for j5=[1 -1] % loop if signed FaceId, treat as separated and add sign if needed
    i5=sign(k3)==sign(j5);
    if any(i5)
     i5=find(i5);
     st3=sprintfc('%i',abs(k3(i5))); r4=cellfun(@abs,unique(st3),'uni',0);
     r4=cellfun(@char,num2cell(unique(cat(2,r4{:}))),'uni',0);
     r2.data=zeros(length(i3(i5))*length(r4),2); i0=1; % Preallocate worst case
     for j4=1:length(r4)
      r5=r2.EltId(i3(i5(~cellfun(@isempty,strfind(st3,r4{j4})))));
      r2.data(i0+(0:length(r5)-1),1)=r5(:);
      r2.data(i0+(0:length(r5)-1),2)=sign(j5).*str2double(r4{j4});
      i0=i0+length(r5);
     end
    end
   end
   r2.data=r2.data(r2.data(:,1)>0,:);
  else % direct EltId
   r2.data=r2.EltId(logical(sum(r2.SConn(:,i1),2)));
  end
  
 else % Standard set format (and old NodeCon format)
  if comstr(st1,'con'); % name:con allow specification of connected to node
   if ~isfield(r2,'NodeCon');error('Missing set.NodeCon');end
   st1=comstr(st1,4);
   if comstr(st1,'~='); i1=1; st1=comstr(st1,3); else; i1=0; end
   i2=comstr(st1,-1);if isempty(i2);error('No connected node');end
   i2=find(any(r2.NodeCon(i2,:),1));
   i2=ismember(r2.data(:,3),i2);
   if i1; r2.data=r2.data(~i2,:);else; r2.data=r2.data(i2,:);end
  elseif comstr(st1,'f') % allow face id instead of nodeid
   st1=comstr(st1,2);
   if comstr(st1,'~='); i1=1; st1=comstr(st1,3); else; i1=0; end
   i2=comstr(st1,-1);
   if i1; i2=setdiff(1:size(r2.NodeCon,2),i2); end
   r2.data=r2.data(ismember(r2.data(:,3),i2),:);
  elseif comstr(st1,'exclude') % setname toto:exclude
   if size(r2.data,2)>1 % exclusion on face selection
    if isfield(r2,'ConvFcn')
     r2=FaceIdFromSelface(r2,'resolve',evalin('caller','model.Elt'));
    end
    elidf=FaceIdFromSelface(elt);
    if ~any(elidf(:,2)) % we still have volume elts
     el1=feutil('selelt selface-All',elt); elidf=FaceIdFromSelface(el1);
    end
    elidf(elidf(:,1)==0,:)=[]; % gets [oldEltId fid] with 0 at Inf lines
    [r2.data,i5]=setdiff(elidf,abs(r2.data),'rows');
   else % exclusion with elements
    r2.data=setdiff(EltId,r2.data); r2.data(r2.data==0)=[];
   end
  elseif comstr(st1,'underlying') % come back to underlying element if existing
   if size(r2.data,2)>1
    r2.data=r2.data(:,1); r2.type='EltId';
   end
  end % specific processing of sets
 end % set format
 
if isfield(r2,'type')&&strncmpi(r2.type,'FaceId',3);
  j1=j1+1;
  while j1<=size(Stack,1)
   Bole=Stack{j1,1};
   if ~isequal(Stack{j1,2},'Set')||~ischar(Stack{j1,4}); break;end
   try; r3=ModelStack{strcmpi(ModelStack(:,2),Stack{j1,4}),3};
   catch % not trivial, hard resolve and see if compatible
    model=evalin('caller','model');
    [i5,el1]=FeutilMatchSet(ModelStack,Stack(j1:end,:),elt); clear model
    if strcmp(i5,'faces'); % gen r3
     r3=FaceIdFromSelface(el1); r3(r3(:,1)==0,:)=[];
     r3=struct('type','FaceId','data',r3);
    else; sdtw('_nb','Resolution problem for %s',Stack{j1,4});
    end
    j1=size(Stack,1);
   end
   if strcmp(r3.type,'FaceId')
    if isfield(r2,'ConvFcn')
     r2=FaceIdFromSelface(r2,'resolve',evalin('caller','model.Elt'));
    end
    if isfield(r3,'ConvFcn')
     r3=FaceIdFromSelface(r3,'resolve',evalin('caller','model.Elt'));
    end
    if strcmp(Bole,'&')
     r2.data=intersect(r2.data,r3.data,'rows');
    elseif strcmp(Bole,'|')
     r2.data(end+[1:size(r3.data,1)],:)=r3.data;
    elseif strcmp(Bole,'&~')
     r2.data=setdiff(r2.data,r3.data,'rows');
    else; error('Not a valid case');
    end
   end
   j1=j1+1;
  end
  if comstr(st1,'sub') %name:sub (i)
    if size(r2.data)<3;error('expecting face subset identifier');end
    i4=comstr(st1(4:end),-1);r2.data=r2.data(ismember(r2.data(:,3),i4),:);
  end
  if j1>1; evalin('caller',sprintf('j1=j1+%i;',j1-2)); end
  % now gen elements, trying to be robust to current Transform or not
  elidf=FaceIdFromSelface(elt);
  if RO.Transformed&&any(elidf(:,2)) % asssumed transformed we have an ID
   if isfield(r2,'ConvFcn')
    r2=FaceIdFromSelface(r2,'resolve',evalin('caller','model.Elt'));
   end
   elt(isfinite(elt(:,1))&~ismember(elidf,r2.data,'rows'),:)=[];
   elt=feutil('OptimEmptyGroups',elt); i4='faces';
  else % underlying volume elements
   elt(~ismember(EltId,[r2.data(:,1);0]),:)=[];
   elt=feutil('getedgepatch',struct('Node',[],'Elt',elt),r2); i4='faces';
  end
  
elseif isfield(r2,'type')&&strncmpi(r2.type,'NodeId',3)&&isfield(r2,'weight')
 wrn=''; i4='faces';
 if iscell(r2.data);
  stn=r2.data; i5=cellfun(@isnumeric,stn);
  if any(i5); stn(i5)=cellfun(@(x)sprintf('nodeid %s',num2str(x(:)')),stn(i5),'uni',0); end
  stn=sprintf('%s |',stn{:}); stn(end)=''; 
 elseif isnumeric(r2.data); stn=sprintf('nodeid %s',num2str(r2.data(:)')); 
 else; stn=r2.data;
 end
 el1=feutil(sprintf('selelt withnode{%s} & selface & innode{%s}',...
  stn,stn),evalin('caller','model'));
 if isempty(el1)
  wrn=sprintf('node based face set "%s" does not not show interior nodes',ModelStack{opt,2});
  el1=feutil(sprintf('selelt withnode{%s} & selface & withnode{%s}',...
   stn,stn),evalin('caller','model'));
  if ~isempty(el1); wrn=[wrn ', surface extended to contain declared nodes']; end
 end
 if ~isempty(el1); elt=el1;
 else; elt=[];
  wrn=sprintf('could not find a surface using node base definition of set "%s"',ModelStack{opt,2});
  [i4,elt]=feutil(sprintf('findelt withnode{%s}',stn),evalin('caller','model'));
  elt=evalin('caller','model.Elt');
  % xxx not sure this mock-up should be done
  %if isempty(elt); elt=feutil('addelt',[],'mass1',double(r2.data(:))); i4=[]; RO.Transformed=j1; end
 end
 if ~isempty(wrn);sdtw('_nb',wrn); end
  
elseif ischar(r2.data) % selection is provided
 if comstr(Stack{j1,2},'SetF')
  r2.data=strrep(r2.data,'setname','safesetname');
  r2.data=strrep(r2.data,'safesafesetname','safesetname'); % remove duplicates
 end
 i4=feutil('findelt',evalin('caller','model'),r2.data);
elseif size(r2.data,2)==2 % Matching faces based on set information
 [EGroup,nGroup]=getegroup(elt);
  for jGroup=1:nGroup
   [ElemF,i1,ElemP]= feutil('getelemf',elt(EGroup(jGroup),:),jGroup);
   if any(strcmp(ElemF, ...  % 2-D elements
     {'q4p','t3p','q5p','q8p','q8a','q9a','t3a','t3p','t6a','t6p'})) || ...
     strcmpi(r2.type,'edgeid')
    FaceIndex=feval(ElemP,'edge')';
    if size(FaceIndex,1)==2; i4(end+1,1:6)=[Inf abs('beam1')]; %#ok<AGROW>
    else; i4(end+1,1:6)=[Inf abs('beam3')];%#ok<AGROW>
    end
   else
    FaceIndex=feutil_get_face(ElemP,r2)';
    if size(FaceIndex,1)==4; i4(end+1,1:6)=[Inf abs('quad4')];%#ok<AGROW>
    else; error('Not implemented');
    end
   end
   if length(opt)>1; error('Multiple set selection not implemented');end
   if size(FaceIndex,2)<max(r2.data(:,2)); FaceIndex(end,max(r2.data(:,2)))=0; end
   nind=sparse(EltId+1,1,1:size(elt,1)); ind=full(nind(r2.data(:,1)+1));
   r1=ind(:,ones(size(FaceIndex,1),1))+ ...
    (FaceIndex(:,abs(r2.data(:,2)))-1)'*size(elt,1);
   r1(any(r1<0,2),:)=[]; % mixed elt in sel: tria group cannot match quad ref
   if ~isempty(r1)
    i4(end+[1:size(r1,1)],1:size(r1,2))=elt(r1);
   end
  end % loop on groups
  elt=i4; i4='faces'; if strcmpi(r2.type,'edgeid'); i4='edges'; end
else % set matching
 if isfield(r2,'type')&&~strncmpi(r2.type,'eltid',5)
  sdtw('_nb','selection on elements with a set of type %s',r2.type);
 end
 if RO.EltIdFixed>0&&~RO.Transformed&&RO.MultiCall
  if ~isfield(RO,'EEid') % gen once if possible
   RO.EEid=sparse(EltId+1,1,1:length(EltId));
   assignin('caller','RunOpt',RO);
  end
  i5=r2.data; i5(i5==0)=[]; 
  if isempty(i5); i4=[];
   if ~feutil('Silent'); sdtw('_nb','No valid EltId in set %s',ModelStack{opt,2}); end
  else; i4=full(RO.EEid(i5+1)); i4(i4==0)=[];
  end
 else % no EEid
  i4=find(ismember(EltId,r2.data)&isfinite(elt(:,1)));
 end
end
end % type of output




%% #elt_check ---------------------------------------------------------
%-- Check Penta6 / Hexa8 element --%
%-- Additionnal functions to check meshes --%
%-- Mathieu Corus -- %
function t_node=elt_check(t_node,elt_type,model,NNode,bypass);

switch comstr(elt_type,-27)
  case 'hexa8'
    
    t_node3=unique(t_node);
    %-- check if hexa8 is correct --%
    rls=[0 0 0 1
      1 0 0 1
      0 1 0 1
      1 1 0 1
      0 0 1 1
      1 0 1 1
      0 1 1 1
      1 1 1 1
      0.5 0.5 0 1
      0.5 0.5 1 1
      0.5 0 0.5 1
      0.5 1 0.5 1
      0 .5 0.5 1
      1 0.5 0.5 1];
    rls(:,1:3)=rls(:,1:3)*2;
    rls(:,1:3)=rls(:,1:3)-ones(14,3);
    rules=integrules('hexa8',rls); % might want to add points
    rules.nodeE=model.Node(NNode(t_node),5:7); of_mk('buildndn',3,rules);
    if bypass==0    
      if all(rules.jdet>0); return; end; %-- if hexa8 is correct, do nothing --%
    end;
    disp(['Processing Hexa8 with nodes ' num2str(t_node(:)') ' ...']);
    pause(.01)
    %-- if not, find the best hexa8 possible --%
    
    % get convex hull for point cloud --%
    [hull,v]=convhulln(model.Node(NNode(t_node),5:7));
    vec1=model.Node(NNode(t_node(hull(:,2))),5:7)-model.Node(NNode(t_node(hull(:,1))),5:7);
    vec2=model.Node(NNode(t_node(hull(:,3))),5:7)-model.Node(NNode(t_node(hull(:,2))),5:7);
    nor=cross(vec1,vec2); for i1=1:length(nor); nor(i1,:)=nor(i1,:)./norm(nor(i1,:)); end;
    
    %-- True Hexa8 element --%
    %-- not degenerated element --%
    if length(t_node3)==8
      % get the planest "face" 
      tt=unique(hull(:)); if length(tt)<8; disp(['Warning : Hexa8 element with nodes ' num2str(t_node(:)') ...
            ' is not convex']); end; 
      ori=[];
      for i1=1:length(nor)-1; ori(i1,1:length(nor)-i1)=nor(i1,:)*nor(i1+1:end,:)'; end;
      [val,x1]=max(ori); 
      max_val=find(val==max(val));
      t_node2=[]; deter=[];
      for k1=1:length(max_val) %-- if maximum is not unique, get best hexa8 --%                   
        
        [val,x1]=max(ori);
        y1=max_val(k1);x1=x1(y1);y1=y1+x1;
        face1=hull(x1,:)';face2=hull(y1,:)';
        n1=fe_c(face1,face2,'dof',2);n2=fe_c(face2,face1,'dof',2);
        if isscalar(n1) && isscalar(n2)  
          face=[fe_c(face1,face2,'dof',1)' n1 n2];face=face([1 3 2 4]);  
          
          %-- check face orientation --%                   
          oth_no=fe_c([1:8]',face(:),'dof',2);
          mid_node_f=sum(model.Node(NNode(t_node(face)),5:7))/size(face(:),1);
          mid_node_o=sum(model.Node(NNode(t_node(oth_no)),5:7))/size(oth_no(:),1);
          ve=model.Node(NNode(t_node(face([2 4]))),5:7)-...
            model.Node(NNode(t_node(face([1 1]))),5:7);
          if cross(ve(1,:),ve(2,:))*(mid_node_o-mid_node_f)' <= 0; face=face(end:-1:1); end;
          
          %-- find opposite face --%
          vec1=model.Node(NNode(t_node(oth_no)),5:7)-...
            model.Node(NNode(t_node(face(ones(length(oth_no),1)))),5:7);
          for i1=1:length(vec1); vec1(i1,:)=vec1(i1,:)./norm(vec1(i1,:)); end;
          ind=[];                
          for i1=1:length(hull)
            if isempty(fe_c(face',hull(i1,:)','ind')); ind=[ind i1]; end;
          end;
          
          if length(ind)==2
            face1=hull(ind(1),:)';face2=hull(ind(2),:)';
            n1=fe_c(face1,face2,'dof',2);n2=fe_c(face2,face1,'dof',2);
            face2=[fe_c(face1,face2,'dof',1)' n1 n2];face2=face2([4 2 3 1]);
            
            %-- check face orientation --%                   
            oth_no=fe_c([1:8]',face2(:),'dof',2);
            mid_node_f=sum(model.Node(NNode(t_node(face2)),5:7))/size(face2(:),1);
            mid_node_o=sum(model.Node(NNode(t_node(oth_no)),5:7))/size(oth_no(:),1);
            ve=model.Node(NNode(t_node(face2([2 4]))),5:7)-...
              model.Node(NNode(t_node(face2([1 1]))),5:7);
            if cross(ve(1,:),ve(2,:))*(mid_node_o-mid_node_f)' >= 0; face2=face2(end:-1:1); end;
            
            deter=[];
            rules=integrules('hexa8',rls);
            for i1=1:4
              rules.nodeE=model.Node(NNode(t_node([face face2([i1:end 1:i1-1])])),5:7);
              of_mk('buildndn',3,rules);
              if ~all(rules.jdet>=0);  deter=[deter -1]; 
              else;deter=[deter sum(rules.jdet)]; 
              end;
            end;
            
            [val,ind]=max(deter);
            t_node2(k1,:)=t_node([face face2([ind:end 1:ind-1])]);
          end; %-- if length(ind)==2
        end; %-- if length(n1) == 1 && length(n2) == 1
      end; %--for k1=1:length(max_val)
      
      t_node2=t_node2(sum(abs(t_node2'))~=0,:);
      deter=[];
      rules=integrules('hexa8',rls);
      for i1=1:size(t_node2,1)
        rules.nodeE=model.Node(NNode(t_node2(i1,:)),5:7);
        of_mk('buildndn',3,rules);
        if ~all(rules.jdet>=0) ; deter=[deter -1]; else;deter=[deter sum(rules.jdet)]; end;
      end;
      [val,ind]=max(deter); 
      if val > 0; t_node=t_node2(ind,:); else;disp(['Warning : failed to find adequate Hexa8 with nodes '  ...
            num2str(t_node(:)')]); 
        disp(['     -> Try to build Penta6 elements instead']);
        disp(' ');
      end;
      
      
      
    elseif length(t_node3)==7; disp(['Warning : 7 nodes degenerated Hexa8 is not supported - element with nodes '...
          num2str(t_node(:)')]);
    elseif length(t_node3)==6 %-- assemble Penta6 as degenerated Hexa8 --%
      t_node=elt_check('penta6',t_node3,model,NNode);
      t_node=t_node([1 2 3 3 4 5 6 6]);
    else;error('Not supported yet')
    end; %--length(t_node3)==8
    
  case 'penta6'
    
    t_node3=unique(t_node);    
    %-- check if penta6 is correct --%
    alt_node=[];
    if length(t_node)==8
      alt_node=t_node;
    elseif length(t_node)==6
      alt_node=t_node([1 2 3 3 4 5 6 6]);
    end;
    if ~isempty(alt_node)
      rules=penta6b('constants',[],[],[]);
      rules.nodeE=model.Node(NNode(alt_node),5:7);
      of_mk('buildndn',3,rules);
      if all(rules.jdet>=0); return; end; %-- if penta6 is correct, do nothing --%
    end;
    %-- if not, find the best hexa8 possible --%
    % get convex hull for point cloud --%
    hull=convhulln(model.Node(NNode(t_node3),5:7));
    vec1=model.Node(NNode(t_node3(hull(:,2))),5:7)-model.Node(NNode(t_node3(hull(:,1))),5:7);
    vec2=model.Node(NNode(t_node3(hull(:,3))),5:7)-model.Node(NNode(t_node3(hull(:,2))),5:7);
    nor=cross(vec1,vec2); for i1=1:length(nor); nor(i1,:)=nor(i1,:)./norm(nor(i1,:)); end;
    
    % get the best opposite faces --% 
    tt=unique(hull(:)); if length(tt)~=6; error(['Penta6 element with nodes ' num2str(t_node(:)') ...
          'is not valid']); end; 
    ori=[];
    for i1=1:length(nor)-1;ori(i1,1:length(nor)-i1)=nor(i1,:)*nor(i1+1:end,:)'; end;
    [val,x1]=min(ori); 
    %-- check for all "opposite" faces --%
    min_val=find(val<0);
    t_node2=[];
    
    for k1=1:length(min_val) %-- if maximum is not unique, get best hexa8 --%                   
      
      [val,x1]=min(ori);
      y1=min_val(k1);x1=x1(y1);y1=y1+x1;
      face=hull(x1,:)';face2=hull(y1,:)';
      
      %-- check faces orientation --%    
      
      oth_no=fe_c([1:6]',face(:),'dof',2);
      mid_node_f=sum(model.Node(NNode(t_node3(face)),5:7))/size(face(:),1);
      mid_node_o=sum(model.Node(NNode(t_node3(oth_no)),5:7))/size(oth_no(:),1);
      ve=model.Node(NNode(t_node3(face([2 3]))),5:7)-...
        model.Node(NNode(t_node3(face([1 1]))),5:7);
      if cross(ve(1,:),ve(2,:))*(mid_node_o-mid_node_f)' <= 0; face=face(end:-1:1); end;
      face=face([1 2 3 3])';
      
      oth_no=fe_c([1:6]',face2(:),'dof',2);
      mid_node_f=sum(model.Node(NNode(t_node3(face2)),5:7))/size(face2(:),1);
      mid_node_o=sum(model.Node(NNode(t_node3(oth_no)),5:7))/size(oth_no(:),1);
      ve=model.Node(NNode(t_node3(face2([2 3]))),5:7)-...
        model.Node(NNode(t_node3(face2([1 1]))),5:7);
      if cross(ve(1,:),ve(2,:))*(mid_node_o-mid_node_f)' >= 0; face2=face2(end:-1:1); end;
      
      deter=[];faceb=[];
      for i1=1:3
        if i1==1; faceb(i1,:)=face2([i1:end end])'; 
        else;faceb(i1,:)=face2([i1:end 1:i1-1 i1-1])';
        end;
        
        rules=penta6b('constants',[],[],[]);
        rules.nodeE=model.Node(NNode(t_node3([face faceb(i1,:)])),5:7);
        of_mk('buildndn',3,rules);
        if ~all(rules.jdet>=0);  deter=[deter -1]; 
        else; deter=[deter sum(rules.jdet)]; 
        end;
      end;
      [val,ind]=max(deter);
      t_node2(k1,:)=t_node3([face faceb(ind,:)]);
    end; %--k1=1:length(min_val)
    
    t_node2=t_node2(find(sum(abs(t_node2'))),:);
    deter=[];
    rules=penta6b('constants',[],[],[]);
    for i1=1:size(t_node2,1)
      rules.nodeE=model.Node(NNode(t_node2(i1,:)),5:7);
      of_mk('buildndn',3,rules);
      if ~all(rules.jdet>=0);  deter=[deter -1]; else;deter=[deter sum(rules.jdet)]; end;
    end;
    [val,ind]=max(deter); 
    if val > 0; t_node=t_node2(ind,:); else;warning(['Failed to find adequate Penta6 with nodes '  ...
          num2str(t_node(:)')]);
    end;
    t_node=t_node([1 2 3 5 6 7]);
  case 'quad4'
    t_node4=unique(t_node);
    if length(t_node4)<4; return; end;
    t_node4=unique(t_node);
    sens_rot=0;  %-- check convexity --%
    compt=0;
    vec=model.Node([1 1 1]*NNode(t_node(1)),5:7)-model.Node(NNode(t_node(2:4)),5:7);
    si=[cross(vec(1,:),vec(2,:))*cross(vec(2,:),vec(3,:))'];
    ve=diff(model.Node(NNode(t_node([1:end 1])),5:7),1,1);
    ve=[ve;ve(1,:)];for i1=1:size(ve,1); ve(i1,:)=ve(i1,:)./norm(ve(i1,:)); end;
    for k1=1:4 
      n(k1,:)=cross(ve(k1,:),ve(k1+1,:)); 
    end;
    nn=n*n';
    ii=mod(find(nn(1,:)<0)+1,4); %-- find inflexion point

    if isscalar(ii)  %-- non convex element
        if ii==0; ii=ii+4; end;
        disp(['Warning : Quad4 element with nodes ' num2str(t_node(:)') ' is not convex']);
        disp('      -> It will be divided into two degenerated quad4 elements')
        t_node=t_node([ii:4 1:(ii-1)]);
        t_node=[[t_node(1:3) t_node(3)];[t_node([1 3 4]) t_node(4)]];
    else
      while sens_rot==0
        vec=model.Node([1 1 1]*NNode(t_node(1)),5:7)-model.Node(NNode(t_node(2:4)),5:7);
        si=cross(vec(1,:),vec(2,:))*cross(vec(2,:),vec(3,:))';
        if si>0; sens_rot=1; 
        elseif si<0; t_node=t_node(randperm(4)); 
        elseif si==0; t_node=t_node([2 3 4 1]);
        end;
        compt=compt+1;
        if compt>=1000 
          disp(['Warning : failed to generate a proper Quad4 element with nodes ' ...
              num2str(t_node(:)')]);break; 
        end;
      end;
    end;
  
end; %-- switch lower(elt_type)
  
%% #feutil_get_face 
function [r1,i1]=feutil_get_face(ElemF,RunOpt)
i1=[];
if isfield(RunOpt,'ConvFcn')&&~isempty(RunOpt.ConvFcn);
 try;r1=feval(RunOpt.ConvFcn,['conv face.' ElemF]);
  if isempty(r1);r1=fe_super('face',ElemF);end
 catch;sdtw('_nb','%s ConvFace.%s failed',char(RunOpt.ConvFcn),ElemF);
  r1=fe_super('face',ElemF);
 end
 try; i1=feval(RunOpt.ConvFcn,['conv faceNum.' ElemF]);
 catch;sdtw('_nb','%s ConvFaceNum.%s failed',char(RunOpt.ConvFcn),ElemF);
 end
elseif isfield(RunOpt,'FaceCmd');%Allow FaceLin
    r1=fe_super(RunOpt.FaceCmd,ElemF);
else; r1=fe_super('face',ElemF);
end
if isempty(i1); i1=[1:size(r1,1)]'; end

%% #FaceIdFromSelface - - ----------------------------------------------------
function r3=FaceIdFromSelface(elt,model,varargin);
if nargin==1 % recover EltId FaceId in an elt from selface selection
 % elidf=FaceIdFromSelface(elt);
 [EGroup,nGroup]=getegroup(elt);r3=zeros(size(elt,1),2);
 for jGroup=1:nGroup
  [ElemF,i1,ElemP]= getegroup(elt(EGroup(jGroup),:),jGroup);
  cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
  if strcmp(ElemP,'SE')||strncmpi(ElemP,'beam',4)||strncmpi(ElemP,'mass',4);
   continue; % skip SE and 1D/0D elts
  end
  try; i1=feval(ElemP,'prop'); catch; i1=[]; end
  if ~isempty(i1)
   if size(elt,2)<i1(3)+1; elt(end,i1(3)+1)=0; end
   r3(cEGI,:)=elt(cEGI,i1(3)+[0 1]); % r3{jGroup}
  end
 end
 
else
 if ischar(model) % do convertion
  % data=FaceIdFromSelface(data,'f2e')
  % data=FaceIdFromSelface(data,'e2f')
  % data=FaceIdFromSelface(data,'resolve',model);
  % data.data=FaceIdFromSelface(data.data,'resolve',model,data.ConvFcn);
  if strncmpi(model,'f2e',3) % [EltId FaceId; ....] to [EltId FId1 ...FIdn; ...
   r1=elt; if ~isfinite(r1(1)); r1=FaceIdFromSelface(r1); end
   [u1,i1]=sort(r1(:,1),'ascend');  r1=r1(i1,:);
   r5=max(abs(r1(:,2)))+1; r1(r1(:,2)==0,2)=r5;
   r2=sparse(r1(:,1)+1,abs(r1(:,2)),ones(size(r1,1),1));
   r3=sum(r2(:,1:r5-1),2); % number of faces per underlying element
   r4=unique(r3); r4(r4==r5|r4==0)=[];
   out=double.empty(0,max(r4));
   for j3=1:length(r4) % loop on found faces number found
    r6=r1(ismember(r1(:,1),find(r3==r4(j3))-1),:); % sub list for current case
    r6=reshape(r6',[],1);
    r6(setdiff(3:2:numel(r6),1:2*r4(j3):numel(r6)))=[];
    r6=reshape(r6,r4(j3)+1,[])';
    out(end+1:end+size(r6,1),1:size(r6,2))=r6;
   end
   r3=out;

  elseif strncmpi(model,'e2f',3) % [EltId FId1 ...FIdn; ...] to [EltId FaceId; ....]
   [II,JJ,KK]=find(elt); in1=find(JJ==1); r3=double.empty(0,2); 
   for j1=2:max(JJ)
    in2=JJ==j1; r3(end+1:end+sum(in2),:)=[KK(in1(II(in2))) KK(in2)];
   end

  elseif strncmpi(model,'resolve',7)
   r3=elt;  ConvFcn=''; data=elt; if isfield(data,'data'); data=data.data; end
   if isnumeric(data); data=double(data); end
   if isfield(elt,'ConvFcn'); ConvFcn=r3.ConvFcn;
   elseif length(varargin)>1; ConvFcn=varargin{2};
   end
   if ~isempty(ConvFcn)
    if isfield(elt,'ElemP'); r5=elt.ElemP;
    elseif ~isempty(varargin); r5=FaceIdFromSelface(data,varargin{1});
    end
    li=fieldnames(r5); % list of element topo to convert
    for j1=1:length(li)
     r6=feval(ConvFcn,['conv faceNum.' li{j1}]); % Abq->SDT
     r6=sparse(r6,1,1:length(r6));
     data(r5.(li{j1}),2)=full(r6(abs(data(r5.(li{j1}),2)))).*sign(data(r5.(li{j1}),2));
    end
    if isstruct(elt); elt.data=data; elt=feutil('rmfield',elt,'ConvFcn'); r3=elt;
    else; r3=data;
    end
   else; error('resolve without ConvFcn');
   end
  else; error('FaceIdFromSel %s unknown.',model);
  end
  
 else % recover element type info to easen convertion operations
  % data.ElemP=FaceIdFromSelface(data.data,model);
  % data=FaceIdFromSelface(data,model);
  % data=FaceIdFromSelface(data,elt);
  if isstruct(elt); data=elt.data; % FaceId set
  elseif isa(elt,'int32')||size(elt,2)==2; data=elt;
  else; data=FaceIdFromSelface(elt); % elt from selface
  end
  if isstruct(model); el0=model.Elt; else; el0=model; end
  r4=struct;
  [EG,nG]=getegroup(el0); eltid=feutil('eltidskipcheck;',el0);
  for jG=1:nG
   [ElemF,i1,ElemP]= getegroup(el0(EG(jG),:),jG);
   cEGI=EG(jG)+1:EG(jG+1)-1;
   i1=ismember(data(:,1),eltid(cEGI));
   if any(i1)
    if isfield(r4,'ElemP'); r5=r4.ElemP; else; r5=[]; end
    r4.(ElemP)=[r5;find(i1)];
   end
  end
  if isstruct(elt); r3=elt; r3.ElemP=r4;  else; r3=r4;  end
  
 end
end

%% #EdgeQuad : quadratic edge utilities
function out=EdgeQuad(varargin) %#ok<DEFNU>

[CAM,Cam]=comstr(varargin{1},1);carg=2;
if comstr(Cam,'init')
  %% #EdgeQuad.init
  % quadratic edge problem
  model=varargin{carg};carg=carg+1;
  if carg<=nargin; RO=varargin{carg};carg=carg+1;else;RO=struct;end
  RO.Elt0=model.Elt;
  if isfield(RO,'sel');
      model.Elt=feutil(['selelt ' RO.sel],model);
  end
  RO=sdth.sfield('addmissing',RO,feutil('getnormalmap Node',model));
  RO=sdth.sfield('addmissing',RO, ...
      feval(feutilb('@OEdgeElt'),struct,model));  % ,'adjacentsurf'
  EC=integrules('beam3','inline');
  
  [II,JJ]=find(RO.OEdgeElt);RO.Edge=unique(sort(int32([II JJ]),2),'rows')';
  Nr=EC.Nr([0;1]);
  RO.nEdge=zeros(size(RO.Edge,2),3);
  %if isempty(RO.normal); RO=EdgeQuad('LineNormal',model,RO); end
  RO.nind=sparse(RO.ID,1,1:length(RO.ID),size(model.Node,1),1); 
  for jElt=1:size(RO.Edge,2)
   n1=model.Node(RO.Edge(:,jElt),5:7); n3=mean(n1);
   normal=RO.normal(RO.nind(RO.Edge(:,jElt)),:);
   doObs=@(n3)sum(normal.*(Nr*[n1;n3]),2);
   p=sp_util('basis',diff(n1),mean(normal));
   p(:,1)=p(:,1)*norm(diff(n1));
   tgt=Nr*[n1;n3]; 
   % doObs(n3)
   % make tangents |_ to normal and place node in segment/element normal
   % cos(A) = nAx * (Nr(1,1) n1(1,1) +Nr(1,2) n1(2,1) + Nr(1,3)*a 
   %          nAy * (Nr(1,1) n1(1,2) +Nr(1,2) n1(2,2) + Nr(1,3)*b
   % 
   Dir=diff(n1);Dir=Dir/norm(Dir);
   %b=sum(normal.*(Nr*[n1;n3]),2);
   %A=normal.*Nr(:,[3 3 3])*p(:,2); r3=(pinv(A)*-b); 
   %[norm(doObs(n3)) norm(doObs(n3+p(:,2)'*r3))]
   
   b=sum(normal.*(Nr*[n1;n3]),2);
   A=normal.*Nr(:,[3 3 3])*p(:,1:2);   r3=(pinv(A)*-b); n3b=n3+(p(:,1:2)*r3)';
   %[norm(doObs(n3)) norm(doObs(n3b))]   
   %fecom('shownodemark',[n1;n3b])
   if abs(r3(1))>.16; % allow slight edge movement
    n3=n3+sign(r3(1))*.16*p(:,1)'; % xxx do y 
    b=sum(normal.*(Nr*[n1;n3]),2);
    A=normal.*Nr(:,[3 3 3])*p(:,2); r3=(pinv(A)*-b); n3b=n3+p(:,2)'*r3;
    [norm(doObs(mean(n1))) norm(doObs(n3)) norm(doObs(n3b))]
    fecom('shownodemark',[n1;mean(n1);n3b])   
   else;
     n3b=n3+(p(:,1:2)*r3)';%doObs(n3b)
   end   
   %fecom('shownodemark',[n1;n3])
   %fecom('shownodemark',[n1;n3;n1+1e-2*Nr*[n1;n3]])
   %fecom('shownodemark',[n1;n1+1e-3*Nr*[n1;mean(n1)]])
   RO.nEdge(jElt,:)=n3b;
  end
  if nargout==0;
    fecom('shownodemark', RO.nEdge)
  elseif isfield(RO,'EdgeMid')
    model.Node(RO.EdgeMid(:,3),5:7)=RO.nEdge;model.Elt=RO.Elt0;
    out=model;
  end
elseif comstr(Cam,'linenormal')
 %% #EdgeQuad.LineNormal : at nodes : mean of tangents
 model=varargin{carg};carg=carg+1;
 if carg<=nargin; RO=varargin{carg};carg=carg+1;else;RO=struct;end
 
 r1=model.Node(RO.Edge(2,:),5:7)-model.Node(RO.Edge(1,:),5:7);
 r1=diag(sparse(1./sqrt(sum(r1.^2,2))))*r1;
 r2=sparse(double(RO.Edge(1,:)),1:size(RO.Edge,2),1,size(model.Node,1),size(RO.Edge,2))- ...
    sparse(double(RO.Edge(2,:)),1:size(RO.Edge,2),1,size(model.Node,1),size(RO.Edge,2))
else;error('%s',CAM);
end

%% #LSPlane : least square plane search - - ----------------------------------
function out=LSPlane(node)

plane=cross(node(2,:)-node(1,:),node(3,:)-node(1,:)); %normal
plane(4)=-plane*node(1,:)';

if size(node,1)>3 %Optim
 plane=fminsearch(@(x) sum(abs(Point2Plane(node,x))),plane);
end

 out=plane/norm(plane(1:3));
 err=abs(Point2Plane(node,plane)); % absolute distance of all points to plane
 if length(err)>3 % If more than 3 points, mean square estimation of the best plane => error display
  display(sprintf('Estimation of plane : mean error = %f ; std error = %f',mean(err),std(err)));
 end


%% #LSPlane : least square plane search - - ----------------------------------
function out=proj2plane(node,plane)
 d=Point2Plane(node,plane);
 out=node-d*plane(1:3)/norm(plane(1:3));



%% #midPointNormal(p_a,p_b,nor_a,nor_b)
% contributed by Illoul Amran, ENSAM PIMM
function out=midPointNormal(p_a,p_b,nor_a,nor_b)

    v_ab=p_b-p_a;
    
    v_a=cross(v_ab,nor_a);
    v_b=cross(v_ab,nor_b);
    v_c=v_a+v_b;
    med=cross(v_c,v_ab);
    
    Y=med/norm(med);
    x_b=norm(v_ab);
    X=v_ab/x_b;
    
    x_nor_a=nor_a*X';
    y_nor_a=nor_a*Y';
    
    x_nor_b=nor_b*X';
    y_nor_b=nor_b*Y';
    
    T1=-x_nor_a/y_nor_a;
    T3=-x_nor_b/y_nor_b;
    L=x_b;
    
    A=(3*T3*L-T1*L)/(4*T3-4*T1);
    B_a=-T1*(2*(L-2*A)-L)/4;
    %B_b=-T3*(2*(L-2*A)+L)/4;
    if A<0||A>x_b||abs(B_a)>x_b; % Doo nothing if too far
        out=mean([p_a;p_b]);%p_a+.5*X;
    else;    out=p_a+A*X+B_a*Y;
    end
     
 
 %% #Point2Plane : signed distance from node [x y z] to plane [a b c d]- - 
function out=Point2Plane(node,plane)

out=1/norm(plane(1:3))*plane*[node';ones(1,size(node,1))];
out=out(:);


%% #newCoor: new node coordinates operator from existing ones - - ------------
% [n1,i2,i3]=newCoor(FEnode,RO.Faces,RO.fCoef,NNode);
function [n1,i2,i3,r5]=newCoor(FEnode,r1,coef,NNode)
% Initialize coefficients to iso if not given
i1=max(cellfun(@(x)size(x,2),r1(:,1)));
for jGroup=1:size(r1,1) % complete empty data with zeros for each group
 r2=r1{jGroup,1}; if ~isempty(r2); r2(1,end+1:i1)=0; r1{jGroup,1}=r2; end
 if isempty(coef{jGroup}); % guess coef to apply to points
  if ~isempty(r2); coef{jGroup}=diag(sparse((1./sum(r2>0,2))))*logical(r2);end
 elseif isscalar(coef{jGroup}) % coef given apply to points
  if ~isempty(r2) % robust to large data sets
   coef{jGroup}=diag(sparse(coef{jGroup}*ones(size(r2,1),1)))*logical(r2);
  end
 else; coef{jGroup}(1,end+1:i1)=0; % only clean sizes
 end
end
% Avoid duplicate nodes (sort [node coef] matrix and do unique)
[i1,i2]=sort(vertcat(r1{:,1}),2); r2=vertcat(coef{:,1});
for j1=1:size(i2,1); r2(j1,:)=r2(j1,i2(j1,:)); end % clean coef ordering xxx
if size(i1,2)==2; % deal with degenerate edges
    i4=i1(:,1)==i1(:,2);
    r2(i4,:)=ones(nnz(i4),1)*[1 0];
end
epsl=sp_util('epsl');i1=[i1 round(r2/epsl)]; 
[i4,i2,i3]=unique(i1,'rows'); i1(:,end-size(r2,2)+1:end)=r2;
% would have used 'stable' option but too recent, 
% want to keep initial node order:
[i2,i4]=sort(i2); i4=sparse(i4,1,1:length(i4)); i3=full(i4(i3));
r3=i1(i2,size(r2,2)+1:end); i1=i1(i2,1:end-size(r2,2)); 
i5=(r3(:,1)==1); 
% r3 is coef matrix, i1 is corresponding node index
n1=zeros(size(i1,1),7); % preallocate
n1(:,1:4)=(max(FEnode(:,1))+(1:size(n1,1))')*[1 0 0 0]; % fill in Ids
for j1=1:size(i1,2); % compute coordinates
 i4=r3(:,j1)~=0; % nodes concerned
 if any(i4)
  n1(i4,5:7)=n1(i4,5:7)+diag(sparse(r3(i4,j1)))*FEnode(NNode(i1(i4,j1)),5:7);
 end
end
% data for MPC generation : new nodes, origins, weights
n1(i5,1)=i1(i5,1);
if nargout>3; r5={n1(:,1) i1 r3}; end 

%% #isoCellRef: RefineCell input for specific elements coming from DivideElt 
function out=isoCellRef(ElemP,opt,RO)
% build input to RefineCell based on desired refinement applied to reference cell
if nargin<3; RO=struct('Silent',1,'NewQuadDivide',1); 
else % handle defaults in any case
 if ~isfield(RO,'Silent'); RO.Silent=1; end; 
 if ~isfield(RO,'NewQuadDivide'); RO.NewQuadDivide=1; end
end
if iscell(opt); [opt,dv1,dv2,dv3]=deal(opt{:});end
switch ElemP
case 'beam1' 
 out=struct('edge',{{[[2+(1:opt(1)-2)]' ones(opt(1)-2,1)*[1 2]],...
  reshape([ 1-dv1(2:end-1)   dv1(2:end-1)  ],[],2)}},...
  'Elt',{{'beam1',...
  reshape([1 reshape((ones(2,1)*(3:opt(1))),1,[]) 2],[],opt(1)-1)'}});
 
case 'quad4'
 opt=opt-1;
 if RO.NewQuadDivide;
  i3 = dv1(:,ones(size(dv2)))*2-1; i2 = dv2(:,ones(size(dv1)))'*2-1;
 else
  i2 = dv1(:,ones(size(dv2)))*2-1; i3 = dv2(:,ones(size(dv1)))'*2-1;
 end
 xi = [-1 -1 0;-1 1 0;1 1 0;1 -1 0];
 n = (1+i2(:)*xi(:,1)').*(1+i3(:)*xi(:,2)')/4; % isoparametric shape fcn
 % elts
 i4 = [1:opt(1)];i4=i4(ones(opt(2),1),:)';
 i5=[0:opt(2)-1]*(opt(1)+1);i5=i5(ones(opt(1),1),:);
 i4=i4(:)+i5(:); i4 =[i4 i4+opt(1)+1 i4+opt(1)+2 i4+1];
 out=isoCellStruct(n,i4,'quad4',4,2,4); % na,el1,elemp,inode,nedge,nface
  
case 'quadb'
 opt=opt-1; % xxx dv1, dv2, dv3 is not used !!!
 % basic cell
 xi = [-1 -1 0;1 -1 0;1 1 0;-1 1 0;0 -1 0;1 0 0;0 1 0;-1 0 0];
 coef=xi(:,1:2)+1;i1 = 2./opt;
 % w points in isoparametric coordinates
 % xi isoparametric coordinates of nodes
 % na shape functions are w points
 w=[coef(:,1)/opt(1) coef(:,2)/opt(2)]-1; % basic cell
 i2 = [1:size(w,1)]'; i2=i2(:,ones(opt(1),1)); i2=i2(:);
 i3 = [0:opt(1)-1]*i1(1);i3=i3(ones(size(w,1),1),:);i3=i3(:);
 w=w(i2,:);w(:,1)=w(:,1)+i3;				% repeat in x
 i2 = [1:size(w,1)]'; i2=i2(:,ones(opt(2),1)); i2=i2(:);
 i3 = [0:opt(2)-1]*i1(2);i3=i3(ones(size(w,1),1),:);i3=i3(:);
 w=w(i2,:);w(:,2)=w(:,2)+i3;			% repeat in y
 w(1,3)=0;
 na = [(1+w(:,1)*xi(1:4,1)').*(1+w(:,2)*xi(1:4,2)')/4 ...
  (1-w(:,1).^2).*(1-w(:,2))/2      (1-w(:,2).^2).*(1+w(:,1))/2 ...
  (1-w(:,1).^2).*(1-w(:,2))/2      (1-w(:,2).^2).*(1+w(:,1))/2];
 na(:,1:4)=na(:,1:4)-na(:,5:8)*[1 0 0 1;0 0 1 1;0 1 1 0;1 1 0 0]'/2;
 na(abs(na)<10*eps)=0;
 % build elements coherent with generated na
 i4=[1:8]; i4=i4(ones(opt(1)*opt(2),1),:);
 i5=[0:opt(1)*opt(2)-1]'*8;i5=i5(:,ones(1,8));
 el1 = i4+i5; % each cell is independently numbered
 % Base topology has been generated, now remove doublon nodes, and clean renum
 n1=na*xi; [n2,i2]=feutil('AddNode;',n1,n1);
 i3=i2; i4=i2-[1:length(i2)]'==0; % doublon nodes
 i5=i2(i4); i6=[1:length(i5)]'; % renum kept nodes only in logical order
 i6=sparse(i5,1,i6); i3= full(i6(i3));
 na=na(i4,:); el1=i3(el1);
 % now organize for RefineCell call
 out=isoCellStruct(na,el1,'quadb',8,3,8);

case 'quad9'
 opt=opt-1; % xxx dv1, dv2, dv3 is not used !!!
 % basic cell
 xi = [-1 -1 0;1 -1 0;1 1 0;-1 1 0;0 -1 0;1 0 0;0 1 0;-1 0 0;0 0 0];
 coef=xi(:,1:2)+1;i1 = 2./opt;
 % w points in isoparametric coordinates
 % xi isoparametric coordinates of nodes
 % na shape functions are w points
 w=[coef(:,1)/opt(1) coef(:,2)/opt(2)]-1; % basic cell
 i2 = [1:size(w,1)]'; i2=i2(:,ones(opt(1),1)); i2=i2(:);
 i3 = [0:opt(1)-1]*i1(1);i3=i3(ones(size(w,1),1),:);i3=i3(:);
 w=w(i2,:);w(:,1)=w(:,1)+i3;				% repeat in x
 i2 = [1:size(w,1)]'; i2=i2(:,ones(opt(2),1)); i2=i2(:);
 i3 = [0:opt(2)-1]*i1(2);i3=i3(ones(size(w,1),1),:);i3=i3(:);
 w=w(i2,:);w(:,2)=w(:,2)+i3;			% repeat in y
 w(1,3)=0;
 na = [(1+w(:,1)*xi(1:4,1)').*(1+w(:,2)*xi(1:4,2)')/4 ...
  (1-w(:,1).^2).*(1-w(:,2))/2      (1-w(:,2).^2).*(1+w(:,1))/2 ...
  (1-w(:,1).^2).*(1-w(:,2))/2      (1-w(:,2).^2).*(1+w(:,1))/2];
 na(:,1:4)=na(:,1:4)-na(:,5:8)*[1 0 0 1;0 0 1 1;0 1 1 0;1 1 0 0]'/2;
 na(abs(na)<10*eps)=0;
 % build elements coherent with generated na
 i4=[1:9]; i4=i4(ones(opt(1)*opt(2),1),:);
 i5=[0:opt(1)*opt(2)-1]'*9;i5=i5(:,ones(1,9));
 el1 = i4+i5;
 % Base topology has been generated, now remove doublon nodes, and clean renum
 n1=na*xi(1:8,:); [n2,i2]=feutil('AddNode;',n1,n1);
 i3=i2; i4=i2-[1:length(i2)]'==0; % doublon nodes
 i5=i2(i4); i6=[1:length(i5)]'; % renum kept nodes only in logical order
 i6=sparse(i5,1,i6); i3= full(i6(i3));
 na=na(i4,:); el1=i3(el1);
 % now organize for RefineCell call
 out=isoCellStruct(na,el1,'quad9',9,3,8);
 
case 'hexa8'
if isstruct(opt)
 %% #isoCellRefAutoHexa8 -3
 %  mo1=femesh('teststruct hexa27 divide 1 1 1')
 %% feval(feutil('@isoCellRef'),'hexa8',mo1);
 if ~isfield(opt,'Node')
  opt.Node=[(1:size(opt.EC.xi,1))'*[1 0 0 0] opt.EC.xi];
  opt.Elt=feutil('addelt',opt.EC.type,[1:size(opt.EC.xi,1)]);
 end
 r1=opt.Node(:,5:7);r1=r1-min(r1); r1=r1*diag(1./max(r1))*2-1;%unit cube
 EC=integrules(ElemP);
 r1=fix(r1/sp_util('epsl'));xi=fix(EC.xi/sp_util('epsl'));
 out=struct('edge',{{[],[]}},'face',{{[],[]}},'volume',{{[],[]}},'Elt',{{}}); 
 [i1,i2]=ismember(r1,xi,'rows');
 i2(~i1)=size(xi,1)+1:size(r1,1);
 opt.Stack={};opt=feutil('renumber',opt,[opt.Node(:,1) i2]); 
 [EGroup,nGroup]=getegroup(opt.Elt);
 for jGroup=1:nGroup
   ElemF= getegroup(opt.Elt(EGroup(jGroup),:),jGroup);
   cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
   out.Elt(jGroup,1:2)={ElemF,opt.Elt(cEGI,feval(ElemF,'node'))};
 end

 i3=feval(ElemP,'edges');opt.icorner=i2(i1);
 r1(i1,:)=[]; i2(i1)=[];% Remove corners
 getR=@(x,y)[y(2)-x x-y(1)]/diff(y);
 for j1=1:size(i3,1)
  r2=xi(i3(j1,:),:);i4=find(any(diff(r2),1));
  r3=getR(r1(:,i4),r2(:,i4));  
  i4=1:size(r1,1);%i4=find(abs(sum(r3,2)-1)<1e-6);r3=r3(i4,:);
  i5=sum((r1(i4,:)- r3*r2).^2,2)>1;i4(i5)=[];r3(i5,:)=[]; % not on edge
  i5=[i2(i4) repmat(i3(j1,:),length(i4),1)];
  r1(i4,:)=[]; i2(i4)=[];% Removed used nodes
  out.edge{1}(end+(1:size(i5,1)),1:size(i5,2))=i5; 
  out.edge{2}(end+(1:size(i5,1)),1:size(r3,2))=r3; 
 end
 %i3=feval(ElemP,'faces');
 i3=[1,2,3,4;1,4,8,5;1,2,6,5;5,6,7,8;2,3,7,6;4,3,7,8];%order for rst
 N=feval(feutilb('@NFromRST'),'quad4',[]);
 for j1=1:size(i3,1) % Find within unit square
  r2=xi(unique(i3(j1,:),'stable'),:);
  r3=N(r1(:,(any(diff(r2),1)))*sp_util('epsl'));
  i4=1:size(r1,1);i5=sum((r1(i4,:)- r3*r2).^2,2)>1;
  i4(i5)=[];r3(i5,:)=[]; % not on edge
  i5=[i2(i4) repmat(i3(j1,:),length(i4),1)];
  r1(i4,:)=[]; i2(i4)=[];% Removed used nodes
  out.face{1}(end+(1:size(i5,1)),1:size(i5,2))=i5; 
  out.face{2}(end+(1:size(i5,1)),1:size(r3,2))=r3; 
 end
 % Now interior nodes
 r3=feval(feutilb('@NFromRST'),'hexa8',r1*sp_util('epsl'));
 i4=[i2 repmat(1:size(xi,1),length(i2),1)];
 out.volume={i4,r3};
 
else
 i1=dv1(:)*2-1;i1= i1(:,ones(1,opt(3)*opt(2))); i3=i1(:);
 i1= ones(opt(1),1)*(dv2(:)*2-1)';i1=i1(:)*ones(1,opt(3)); i3(:,2)=i1(:);
 i1= ones(opt(1)*opt(2),1)*(dv3(:)*2-1)'; i3(:,3)=i1(:);
 % this is an isoparametric mapping with linear elements
 i2 = [-1 -1 -1;1 -1 -1;1 1 -1;-1 1 -1;-1 -1 1;1 -1 1;1 1 1;-1 1 1];
 coef =[0 0 0;-1 -1 -1;1 -1 -1;1 1 -1;-1 1 -1;-1 -1 1;1 -1 1;1 1 1;-1 1 1];
 na=zeros(size(i3,1),8);
 for j1=1:8
  n=i3.*coef(ones(size(i3,1),1)*(j1+1),1:3)+1;
  na(:,j1)=(n(:,1).*n(:,2).*n(:,3)/8);
 end
 % row of elements
 i4 = [1:opt(1)-1;2:opt(1)]';
 i4=[i4 i4(:,[2 1])+opt(1)];i4=[i4 i4+opt(1)*opt(2)];
 % face of elements
 i5=[0:opt(1):opt(1)*(opt(2)-1)-1];i5=i5(ones(size(i4,1),1),:);i5=i5(:);
 i6=[1:size(i4,1)]';i6=i6(:,ones(opt(2)-1,1)); i6=i6(:);
 i4=i4(i6,:)+i5(:,ones(size(i4,2),1));
 % volume of elements
 i5=[0:opt(1)*opt(2):opt(1)*opt(2)*(opt(3)-1)-1];
 i5=i5(ones(size(i4,1),1),:);i5=i5(:);
 i6=[1:size(i4,1)]';i6=i6(:,ones(opt(3)-1,1)); i6=i6(:);
 i4=i4(i6,:)+i5(:,ones(size(i4,2),1));
 el1=i4;
 % now organize for RefineCell call
 out=isoCellStruct(na,i4,'hexa8',8,2,4);
end

case 'hexa20'
 opt=opt-1;
 % basic cell on -1 1 interval
 coef=[-1 -1 -1;1 -1 -1;1 1 -1;-1 1 -1;-1 -1 1;1 -1 1;1 1 1;-1 1 1
  0 -1 -1;1 0 -1;0 1 -1;-1 0 -1;-1 -1 0;1 -1 0;1 1 0;-1 1 0; %9:16
  0 -1 1;1 0 1;0 1 1;-1 0 1];
 i4=(coef+1)/2;
 i2=i4(:,1)*dv1(2:end)'+(1-i4(:,1))*dv1(1:end-1)';
 i4=[i2(:) reshape(i4(:,2*ones(1,opt(1))),20*opt(1),1) ...
  reshape(i4(:,3*ones(1,opt(1))),20*opt(1),1)] ; % x repeat
 i2=i4(:,2)*dv2(2:end)'+(1-i4(:,2))*dv2(1:end-1)';
 i4=[reshape(i4(:,ones(1,opt(2))),size(i4,1)*opt(2),1) i2(:) ...
  reshape(i4(:,3*ones(1,opt(2))),size(i4,1)*opt(2),1)] ; %y repeat
 i2=i4(:,3)*dv3(2:end)'+(1-i4(:,3))*dv3(1:end-1)';
 i4=[reshape(i4(:,ones(1,opt(3))),size(i4,1)*opt(3),1) ...
  reshape(i4(:,2*ones(1,opt(3))),size(i4,1)*opt(3),1) i2(:)] ; % z repeat
 i4=i4*2-1;
 %initial i4; mapped i3; isoparametric mapping with linear elements
 coef=[0 0 0;coef]; na=zeros(size(i4,1),8);
 for j1=1:8
  n=i4.*coef(ones(1,size(i4,1))*(j1+1),1:3)+1;
  na(:,j1)=n(:,1).*n(:,2).*n(:,3)/8;
 end
 i4=[1:20]; i4=i4(ones(opt(1)*opt(2)*opt(3),1),:);
 i5=[0:opt(1)*opt(2)*opt(3)-1]'*20;i5=i5(:,ones(1,20));
 el1 = i4+i5;
 % Base topology has been generated, now remove doublon nodes, and clean renum
 n1=na*coef(2:9,:); [n2,i2]=feutil('AddNode;',n1,n1);
 i3=i2; i4=i2-[1:length(i2)]'==0; % doublon nodes
 i5=i2(i4); i6=[1:length(i5)]'; % renum kept nodes only in logical order
 i6=sparse(i5,1,i6); i3= full(i6(i3));
 na=na(i4,:); el1=i3(el1);
 % now organize for RefineCell call
 out=isoCellStruct(na,el1,'hexa20',20,2,4);

case 'tria3'
 % isoparametric refine on reference cell
 na = ones(opt(1),1)*linspace(1,2,opt(1));
 na=[nonzeros(triu(na)) nonzeros(triu(na'))];
 na=[2-na(:,1) na(:,2)-1];na(:,3)=1-sum(na,2);
 % build new elements
 i3=cumsum(1:opt(1));i3=i3(:); el1=[];
 while length(i3)>1
  el1=[el1;i3(1:end-1) i3(2:end) i3(2:end)-1 ];
  if length(i3)>2; el1=[el1;i3(2:end-1)-1 i3(2:end-1) i3(3:end)-1 ]; end
  i3=i3(2:end)-1;
 end
 % now organize for RefineCell call
 out=isoCellStruct(na,el1,'tria3',3,2,3);
 
 % for tetra only do central refine (works once if reasonable initial shapes)
case 'tetra4'
 out=struct('volume',{{[5 1 2 3 4],.25}},...
 'Elt',{{'tetra4',[1 2 3 5;1 3 4 5;1 4 2 5;2 4 3 5]}});
case 'tetra10'
 out=struct('volume',...
 {{[11 1 2 3 4;12 1 2 3 4;13 1 2 3 4;14 1 2 3 4;15 1 2 3 4];
 [0 1/4 1/4 1/4;0 1/8 1/8 1/8;0 5/8 1/8 1/8;0 1/8 5/8 1/8;0 1/8 1/8 5/8]}},...
 'Elt',{{'tetra10',[1 2 3 11 5 6 7 12 13 14; 1 3 4 11 7 10 8 12 14 15;
1 4 2 11 8 9  5 12 15 13; 2 4 3 11 9 10 6 13 15 14 ]}}); 

otherwise; out=[];
end

%% #isoCellStruct: build RefineCell input from a provided reference cell topology
function out=isoCellStruct(na,el1,ElemP,inode,nedge,nface)
% na: coef matrix to generate new nodes from initial nodes
% el1: refined cell topology assuming ordered nodeid in na
% ElemP: element type in the new topology
% inode: number of nodes of the initial cell (not usually corellated to na)
% nedge: nodes on an edge (2 or 3 usually)
% nface: nodes on a face (3 or 4 usually)

% see also sdtweb feutil buildRCIQuad

% detect initial cell nodes, edge dependent, face dependent nodes and other
r1=logical(na); out=struct;
i1=sum(r1,2)==1; i2=sum(r1,2)==nedge; i5=sum(r1,2)==nface; i3=~i1&~i2&~i5;
i1=find(i1); i2=find(i2); i3=find(i3); i5=find(i5);
r1=1:size(r1,1); % base renumbering
% reorder initial nodes to keep reference cell as conform
r2=na(i1,:); i6=0*i1;
for j1=1:length(i1); i6(j1)=find(r2(:,j1)); end
i1=i1(i6);
if length(i1)<inode % not all reference cell nodes are used (2nd order elements)
 r1(length(i1)+1:end)=inode+(1:(length(r1)-length(i1)));
end
NN=sparse([i1;i2;i5;i3],1,r1); % renum between na and edge/face ordering
% renum final elts
out.Elt={ElemP,full(NN(el1))}; 
if size(el1,1)==1; out.Elt{2}=reshape(out.Elt{2}(:)',1,[]); end
% edge, based on nedge nodes
if ~isempty(i2)
 ed=na(i2,:)'; ed=ed(:); i4=(ones(size(i2,1),1)*(1:size(na,2)))';
 i6=ed~=0; ed=reshape(ed(i6),nedge,[])'; i4=reshape(i4(i6),nedge,[])';
 out.edge={[inode+[1:size(i4,1)]' i4],ed};
end
% face, based on nface nodes
if ~isempty(i5)
 ed=na(i5,:)'; i4=(ones(size(i5,1),1)*(1:size(na,2)))';
 i6=ed~=0; ed=reshape(ed(i6),nface,[])'; i4=reshape(i4(i6),nface,[])';
 out.face={[inode+length(i2)+[1:size(i4,1)]' i4],ed};
end
% volume, no assumption on the number of nodes to be used
if ~isempty(i3)
 ed=na(i3,:); i4=(ones(size(i3,1),1)*(1:size(na,2)));
 out.volume={[inode+length(i2)+length(i5)+[1:size(i4,1)]' i4],ed};
end

%% #buildRCIQuad: build RefineCell input for quad version of lin elts - - ----
function out=buildRCIQuad(typ,RO)
switch typ
 case 'hexa8' % build a base model, apply transform to lin, then identify input
  model=struct('Node',...
   [1 0 0 0 -1 -1 -1;   2 0 0 0 1 -1 -1;   3 0 0 0 1 1 -1;   4 0 0 0 -1 1 -1;
   5 0 0 0 -1 -1 1;     6 0 0 0 1 -1 1;    7 0 0 0 1 1 1;    8 0 0 0 -1 1 1
   9 0 0 0 0 -1 -1;    10 0 0 0 1 0 -1;   11 0 0 0 0 1 -1;  12 0 0 0 -1 0 -1;
   13 0 0 0 -1 -1 0;   14 0 0 0 1 -1 0;   15 0 0 0 1 1 0;   16 0 0 0 -1 1 0;
   17 0 0 0 0 -1 1;    18 0 0 0 1 0 1;    19 0 0 0 0 1 1;   20 0 0 0 -1 0 1],...
   'Elt',[Inf abs('hexa20') zeros(1,20-7);
   1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20]);
  [eltid,model.Elt]=feutil('EltIdFix;',model);
  model=feutil('quad2lin',model);
  if isfield(RO,'shift') % clean input
   mo1=feutil('RefineCell',model,struct('hexa8',RO,'set',...
    [1 RO.shift])); %[1+0*RO.shift(:) RO.shift(:)]));
  else; mo1=feutil('RefineCell',model,struct('hexa8',RO));
  end
  mo1=feutil('lin2quad',mo1);
  i4=mo1.Node(:,5:7); i4(abs(i4)<1e-15)=0; mo1.Node(:,5:7)=i4; % clean num err
  coef=[0 0 0;i4]; na=zeros(size(i4,1),8);
  for j1=1:8 % now get coef from mapping with linear elements
   n=i4.*coef(ones(1,size(i4,1))*(j1+1),1:3)+1;
   na(:,j1)=n(:,1).*n(:,2).*n(:,3)/8;
  end
  na(abs(na)<1e-15)=0; % clean num err
  out=feval(feutil('@isoCellStruct'),na,mo1.Elt(2:end,1:20),'hexa20',8,2,4);  
  
 otherwise; error('buildRCIQuad no implemented for %s',typ)
end
r1=intersect(fieldnames(RO),{'shift','faces','orient'});
if ~isempty(r1); for j1=1:length(r1); out.(r1{j1})=RO.(r1{j1}); end; end


%% #StraightenEdge - - -------------------------------------------------------
function [model,edge]=StraightenEdge(model,elt);
[EG,nG]=getegroup(elt);NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
for jG=1:nG
 ElemF=getegroup(elt(EG(jG),:),jG); cEGI=EG(jG)+1:EG(jG+1)-1;
 if ismember(ElemF,{'tetra10','hexa20','penta15','quadb','tria6','beam3'})
  %edge=reshape(elt(isfinite(elt(:,1)),feval(ElemF,'edge')),[],3);
  edge=reshape(elt(cEGI,feval(ElemF,'edge')),[],3);
  edge(:,1:2)=sort(edge(:,1:2),2); edge=unique(edge,'rows');
  model.Node(NNode(edge(:,3)),5:7)=...
   .5*(model.Node(NNode(edge(:,1)),5:7)+model.Node(NNode(edge(:,2)),5:7));
 end
end

%% #allEltSel: aggregate mpid and set information for addset
function [sel,name,typ,ConvFcn]=allEltSel(model)
mpid=feutil('mpid',model); eltid=feutil('eltid',model);
i1=isfinite(model.Elt(:,1)); mpid=mpid(i1,:); eltid=eltid(i1,:);
group={'MatId','ProId','GroupId'}; r1=cellfun(@(x) sprintf('%s %%i',x),group,'uni',0);
i1={unique(mpid(:,1)) unique(mpid(:,2)) unique(mpid(:,3))};
r1=cellfun(@(x,y) sprintfc(x,y),r1,i1,'uni',0);
typ={'EltId','EltId','EltId'};
for j2=1:3;i1{j2}=cellfun(@(x) eltid(mpid(:,j2)==x),num2cell(i1{j2}),'uni',0);end
r2=stack_get(model,'set'); r2=r2(cellfun(@isstruct,r2(:,3)),:);
r3={'eltid','EltId','EltId';'faceid','FaceId','FaceId';
 'edgeid','EdgeId','EdgeId';'nodeid','NodeId','NodeId'};
ConvFcn=cellfun(@(x,y) repmat({x},length(y),1),{'','',''},i1,'uni',0);
for j1=1:size(r3,1)
 r4=r2(cellfun(@(x) isequal(lower(x.type),r3{j1,1}),r2(:,3)),:);
 i2=cellfun(@(x)~isnumeric(x.data),r4(:,3));
 if any(i2)
  r5=r4(i2,2); r4=r4(~i2,:);
  sdtw('_nb','%s set(s) %s is(are) not resolved, thus not used for EltSet',r3{j1,3},...
   sprintf('%s ',r5{:}))
 end
 r1{end+1}=r4(:,2); group{end+1}=r3{j1,2}; typ{end+1}=r3{j1,3};
 i1{end+1}=cellfun(@(x)double(x.data),r4(:,3),'uni',0);
 i2=cellfun(@(x) isfield(x,'ConvFcn'),r4(:,3));
 st=cell(length(i2),1); st(:)={''};
 if any(i2); st(i2)=cellfun(@(x) x.ConvFcn,r4(i2,3),'uni',0); end
 ConvFcn{end+1}=st;
end
group=cellfun(@(x,y)repmat({x},length(y),1),group,i1,'uni',0);
typ=cellfun(@(x,y) repmat({x},length(y),1),typ,i1,'uni',0);
sel=cat(1,i1{:}); name=[cat(1,r1{:}) cat(1,group{:})]; typ=cat(1,typ{:});
ConvFcn=cat(1,ConvFcn{:});

%% #catSubEid: transform Face/Edge selection to agregated format
function r3=catSubEid(r3,id)
if id==2; r3(:,2)=10.^(r3(:,2)-1); coef=1; else; coef=10; end % bin.num identifiers
r3=sortrows(r3,2); r5=r3; [i3,i4]=unique(r3(:,1));
i4=sort(i4); r3=r3(i4,:); i5=setdiff(1:size(r5,1),i4(:)');
i6=sparse(r3(:,1),1,1:size(r3,1));
while ~isempty(i5) % remaining elts identifiers to add
 r5=r5(i5,:); [i3,i4]=unique(r5(:,1)); i4=sort(i4); r6=r5(i4,:);
 r10=min(r3(i6(r6(:,1)),2),r6(:,2)); r11=max(r3(i6(r6(:,1)),2),r6(:,2));
 r3(i6(r6(:,1)),2)=coef*r10+r11;
 i5=setdiff(1:size(r5,1),i4(:)');
end

%% #addImplicitMat: recover MatId in advanced il entries (composites...)
function r3=addImplicitMat(model,mpid,i1,r3);
if i1==1 % mat/pro xxx
 il=fe_mat('getil',model);
 if ~isempty(il)
  r4=unique(mpid(ismember(mpid(:,i1),r3),2));
  il=il(ismember(il(:,1),r4),:);
 end
 for j1=1:size(il,1)
  [st,unit,i2]=fe_mat('type',il(j1,2));
  if ~exist(st,'file'); st(1)='p'; end
  try; [r5,r6]=feval(st,'propertyunittypecell',i2);
  catch;
   st1=setdiff({'p','m'},st(1));
   [r5,r6]=feval(sprintf('%s%s',st1{1},st(2:end)),'propertyunittypecell',i2);
  end
  if size(il,2)<length(r5)+length(r6); il(end,length(r5)+length(r6))=0; end
  i2=~cellfun(@isempty,strfind(r5(:,1),'MatId'));
  if any(i2); r3=[r3;reshape(setdiff(unique(il(j1,i2)),0),[],1)]; end
  if ~isempty(r6)
   i2=find(~cellfun(@isempty,strfind(r6(:,1),'MatId')));
   il1=il(j1,:); il1(1:size(r5,1))=[];
   for j2=1:length(i2)
    r3=[r3;reshape(setdiff(unique(il1(i2(j2):size(r6,1):end)),0),[],1)];
   end
  end
 end
end

function out=toTRg(TR,cf)
%% #toTRg transform to global
if ischar(TR);  eval(iigui(TR,'MoveFromCaller')); end
if isfield(TR,'isGlobal')&&TR.isGlobal; out=TR; return;end
  r1=stack_get(cf,'info','GNodeBas','g');
  if ~isempty(r1) % Not present in wire-frame
   cGL = basis('trans le',r1.bas,r1.Node,TR.DOF);
   TR.def=cGL*TR.def;TR.isGlobal=1; TR.cGL=sparse(cGL);
  end
  out=TR;

function  [RO,st,CAM]=paramEdit(fmt,ROCAM)

if exist('cingui','file')
[RO,st,CAM]=cingui('paramedit -DoClean',fmt,ROCAM);
else; 
 % simplified paramedit for openfem
 RO=ROCAM{1};CAM=ROCAM{2}; st=CAM;
 st=textscan(regexprep(fmt,'(\([^\)]*)\)',' '),'%s');st=st{1};
 Cam=lower(CAM);
 for j1=1:length(st)
  [CAM,Cam,RO.(st{j1})]=comstr(st{j1},[-25 3],CAM,Cam);
 end
end
Cam=lower(CAM);
