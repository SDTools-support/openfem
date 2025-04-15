function [EGroup,nGroup,ElemP]=getegroup(elt,jGroup)

%GETEGROUP get element group positions and number of groups
%
%   Syntax : [EGroup,nGroup]=getegroup(elt);
%            [EGroup,nGroup,names]=getegroup(elt);
%            [ElemF,i1,ElemP]=getegroup(model.Elt(EGroup(jGroup),:),jGroup);
%
%   Typical format for a group loop
%
%            for jGroup=1:nGroup
%             [ElemF,i1,ElemP]= getegroup(model.Elt(EGroup(jGroup),:),jGroup);
%              cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
%
%              for jElt=1:length(cEGI)
%              end
%            end
%  For node renumbering use
%   NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));


%	Etienne Balmes
%       Copyright (c) 2001-2025 by INRIA and SDTools
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       All Rights Reserved.

%#ok<*NOSEM>

if ischar(elt)
  if strcmp(elt,'cvs');
    EGroup='$Revision: 1.19 $  $Date: 2025/04/14 12:57:01 $';
  elseif strcmpi(elt,'names')  % #Names getegroup('names',model)
    [st,nGroup,EGroup]=getegroup(jGroup);
  elseif strncmpi(elt,'eltconnect',10); 
  %% build component connectivity matrix
  model=jGroup;[EGroup,nGroup]=getegroup(model.Elt);
  NNode=sparse(model.Node(:,1)+1,1,1:size(model.Node,1)); % allow 0 node
  RO.Use=comstr(elt(11:end),-1);
  if isempty(RO.Use);error('Not specified connected node');end
  RO.N=length(NNode);
  RO.ind=zeros(size(model.Elt,1),1);
  for jGroup = 1:nGroup %loop on element groups
   ElemF= feutil('getelemf',model.Elt(EGroup(jGroup),:),jGroup);
   cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;
   i1=fe_super('node',ElemF); % node indices
   i2 = [length(i1) length(cEGI) length(i1)*length(cEGI)];
   i3 = 1:i2(2);
   % build matrix of node use for elements

   i2=[2:length(i1) 1];
   in1=reshape(full(NNode(model.Elt(cEGI,i1)+1)),[],1);
   in2=reshape(full(NNode(model.Elt(cEGI,i2)+1)),[],1);
   i2=(in1&in2);
   i2 = sparse(in1(i2),in2(i2),1,RO.N,RO.N);
   if jGroup==1;RO.Topo=i2;else;RO.Topo=RO.Topo+i2;end
   RO.ind(cEGI,1)=model.Elt(cEGI,i1(1))+1; % nodeid+1 for first node
  end % of loop on groups
   i2=triu(RO.Topo+RO.Topo'); i3=etree(i2);
   i4=find(~i3); i3(i4)=-[1:length(i4)]; % no leaf
   i4=find(i3>0);while ~isempty(i4); i3(i4)=i3(i3(i4));i4=find(i3>0);end;
   i3=-i3; % leaf numbers for each node 
   i4=find(RO.ind);
   RO.ind(i4)=i3(full(NNode(RO.ind(i4))));      
   EGroup=find(ismember(RO.ind,full(i3(NNode(RO.Use+1))))); 
   % Return eltind connected to comp
  end
elseif nargin>1&&ischar(jGroup)&&strcmpi(jGroup,'ElemF')
 % getegroup(mo1.Elt,'elemf')
 ElemF={}; [EGroup,nGroup]=getegroup(elt);
 for jGroup=1:nGroup
  ElemF=unique(cat(1,ElemF,feutil('getelemf',elt(EGroup(jGroup),:))));
 end
 EGroup=ElemF;
 
elseif nargin>1
   st=elt; i1 = [min([find(~st) find(st==32)]) length(st)+1];
   out = char(st(2:i1(1)-1));out1=[];
   if nargout>1 
    out1 = st(i1(1)+1:length(st)); 
    if isempty(out1); out1=jGroup; elseif out1(1)==0; out1=jGroup(1); end
   end
   if nargout>2    ;ElemP=fe_super('parent',out);  end
   EGroup=out; nGroup=out1;
    
elseif isempty(elt)
 EGroup=[1 1]; nGroup=0;
else
  tmp = reshape(find(~isfinite(elt(:,1))),1,[]); %tmp=tmp(:)';
  nGroup = length(tmp); EGroup=[tmp size(elt,1)+1]; 
  if nGroup == 0; warning('no element group specified in ELT'); end %#ok<WNTAG>
  if nargout>2
    ElemP=cell(nGroup,1);
    for jGroup=1:nGroup
           st=elt(EGroup(jGroup),:); 
           i1 = [min([find(~st) find(st==32)]) length(st)+1];
           ElemP{jGroup} = char(st(2:i1(1)-1));
    end
  end
end
