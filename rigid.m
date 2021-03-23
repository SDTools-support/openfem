function [out,out1,out2]=rigid(CAM,varargin);

%RIGID	non-standard element function for the hanling of linearized
%	rigid links
%
%   Syntax : [T,cdof] = rigid(node,elt,mdof)
%            [T,cdof] = rigid(Up)
%
%   RIGID elements are used for a direct elimination of linearized rigid
%       links. A penalization approach is implemented with CELAS elements.
%
%   When using the input structure Up (see UPCOM), Up.Elt is used unless
%       ELT is given, Up.DOF is used unless mdof is given. If coordinate 
%       systems are defined in field Up.bas, PID (position coordinate system) 
%       and DID (displacement coordinate system) are properly handled.
%
%   RIGID elements are ignored during assembly by FE_MK. For the resulting
%       model [m,k,mdof], a basis of the subspace verifying 
%       the constraints is generated using [T,cdof] = rigid(node,elt,mdof) 
%       where you can leave all elements.
%       The constrained model is given by mc=T'*m*T; kc=T'*k*T which are
%       associated to the DOFs cdof. To go back to initial coordinates use
%       the relation  q = T*qc.
%
%	Rigid element property rows follow the format
%         [n1 n2 DofSel 0 0 EltID]
%       where
%	 a DofSel of 123 links only translations, 123456 all DOFs, etc.
%	 A negative DofSel indicates a DOF to DOF connection where the
%	 appropriate DOFs of n2 are eliminated. 
%
%	Warning : when seeking rigid rotations, be sure that the master node
%	  n1 has rotational DOFs. 
%
%   See sdtweb     eltfun, elem0
%	See also help  celas, bar1, beam1, ...


%       Etienne Balmes
%       Copyright (c) 2001-2017 by INRIA and SDTools
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       All Rights Reserved.

%#ok<*ASGLU,*NOSEM>
if comstr(CAM,'cvs')
 out='$Revision: 1.17 $  $Date: 2017/11/24 15:31:56 $'; return;
end
% standard calls with one input argument
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1); 
 if comstr(Cam,'integinfo')

   out=[];out1=[];out2=[];

 elseif  comstr(Cam,'groupinit');   out = '';
 elseif comstr(Cam,'matcall');   out='';out1='';
 elseif comstr(Cam,'call');   out='';
 elseif comstr(Cam,'rhscall');   out='';
 elseif comstr(Cam,'state');   out='';
 elseif comstr(Cam,'node'); out = 1:2;
 elseif comstr(Cam,'prop'); out = [0 0 6];
 elseif comstr(Cam,'dof')
   out =[1+(1:6)/100 2+(1:6)/100]';
 elseif comstr(Cam,'line');  out = [1 2];
 elseif comstr(Cam,'patch'); out = [1 2];
 elseif  comstr(Cam,'face'); out=[];
 elseif comstr(Cam,'sci_face'); out = [1 2 2];
 elseif  comstr(Cam,'edge'); out=[1 2];
 elseif  comstr(Cam,'parent')
   out = 'beam1';
 elseif  comstr(Cam,'num2dof')
   r1=varargin{1};
   out=abs(sprintf('%i',abs(r1)))-48;
 end
return
end % of standard calls with one input argument

% building the rigid matrix
node=CAM;elt=varargin{1};

cGL=[]; 
if ~isstruct(node) 
 Up=[]; mdof=varargin{2};
else
 Up=node; 
 if ~isfield(Up,'Node'); error('First input must give nodes'); end  
 if isfield(Up,'bas')&&~isempty(Up.bas)
  [node,bas]=basis(Up.Node,Up.bas);
  cGL=basis('trans l',bas,Up.Node);
 else; node=Up.Node; 
 end
 if (nargin<2 || isempty(elt))&&isfield(Up,'Elt'); elt=Up.Elt; end
 if nargin<3 && isfield(Up,'DOF')&&~isempty(Up.DOF); mdof=Up.DOF; 
 else; mdof=feutil('getdof',Up); 
 end
end

[EGroup,nGroup]=getegroup(elt);
NNode=sparse(node(:,1),1,1:size(node,1));
RunOpt.ReturnT=1;
if isequal(varargin{end},'c');RunOpt.ReturnT=0;
  RunOpt.II=[];RunOpt.JJ=[];RunOpt.KK=[];RunOpt.slave=[];j0=0;
else;T = speye(length(mdof),length(mdof));
end
ind = round(rem(mdof,1)*100);ind=find(ind>0&ind<7);
nim = 6;
if isempty(ind) % no mechanical DOF

else % some translation DOFs

nind=sparse(fix(mdof(ind))*nim+round(rem(mdof(ind),1)*100)-1,1, ...
              ind,ceil(max(mdof(ind))*nim+nim),1);

cEGI=[];
for jGroup = 1:nGroup %loop on element groups
   [ElemF,i1]= feutil('getelemf',elt(EGroup(jGroup),:),jGroup);
   if comstr(ElemF,'rigid') && (isempty(i1) || i1(1)>=0 )
    cEGI=[cEGI EGroup(jGroup)+1:EGroup(jGroup+1)-1]; %#ok<AGROW>
   end
end

% find repeated slave nodes
i2=elt(cEGI,2);i2=i2(:);
i2=find(sparse(i2+1,1,ones(size(i2,1),size(i2,2)))>1)-1;i2=i2(i2~=0);
if ~isempty(i2)
  st=sprintf('%i ',i2);
  if length(i2)==1; warning(['Slave node' st ' is repeated']);
  else;  warning(['Slave nodes' st ' are repeated']); end
end

% build the constraint matrix 

for j1=cEGI
      i2 = elt(j1,3);i2=abs(sprintf('%i',abs(elt(j1,3))))-48;
      %i2=[rem(i2,10) rem(fix(i2/10),10) rem(fix(i2/100),10) ...
      %   rem(fix(i2/1e3),10) rem(fix(i2/1e4),10) rem(fix(i2/1e5),10)];
      %i2 = find(sparse(abs(i2)+1,1,i2))-1; 

      if elt(j1,3)>0
       % i3 master, i4 slave
       i3=elt(j1,1)*nim+[0:5];   if max(i3)>length(nind); nind(max(i3))=0;end
       i3=full(nind(i3));i3i=find(i3);
       x = node(NNode(elt(j1,[1 2])),5:7);x=x(2,:)-x(1,:);
       x = [0 x(3) -x(2) ;-x(3) 0 x(1);x(2) -x(1) 0];

       i4=elt(j1,2)*nim+i2-1;   if max(i4)>length(nind); nind(max(i4))=0;end
       i4=full(nind(i4));i4i=find(i4);
       if RunOpt.ReturnT;
        r = sparse([eye(3,3) x;zeros(3,3) eye(3,3)]);
        if length(i2)>length(i3i)
           warning('Missing DOFs in Master node %i',elt(j1,1));
        elseif ~any(any(T(:,i3(i3i))))
           warning('Master node is fixed %i %i %i',elt(j1,[1 2 3]));
        end
        T(:,i4(i4i))=T(:,i3(i3i))*r(i2(i4i),i3i)';
        if any(any(T(i4(i4i,:))))
           warning('Slave node problem %i %i %i',elt(j1,[1 2 3]));
        end
       else
        tr=[[-eye(3,3) -x;zeros(3,3) -eye(3,3)] eye(6)];
        tr=tr(i2,:);
        [II,JJ,KK]=find(tr);
        i5=[i3;i4];JJ=i5(JJ);
        RunOpt.II=[RunOpt.II;II+j0];RunOpt.JJ=[RunOpt.JJ;JJ];
        RunOpt.KK=[RunOpt.KK;KK];j0=j0+max(II);
        RunOpt.slave=[RunOpt.slave;i4];
       end
      else %negative value for direct DOF match
       i3=elt(j1,1)*nim+i2-1;   if max(i3)>length(nind); nind(max(i3))=0;end
       i3=full(nind(i3));i3i=find(i3);
       i4=elt(j1,2)*nim+i2-1;   if max(i4)>length(nind); nind(max(i4))=0;end
       i4=full(nind(i4));i4i=find(i4);
       if RunOpt.ReturnT;
          if ~any(any(T(:,i3(i3i))))
           warning('Master node is fixed %i %i %i',elt(j1,[1 2 3]));
          end
          T(:,i4(i4i))=T(:,i3(i3i));
          if any(any(T(i4(i4i),:)))
           warning('Slave node problem %i %i %i',elt(j1,[1 2 3]));
          end
       else;
        i4=i4(i4i); i3=i3(i3i);
        [II,JJ,KK]=find([-eye(length(i4)) eye(length(i4))]);
        i5=[i3 i4];JJ=i5(JJ);
        RunOpt.II=[RunOpt.II;II+j0];RunOpt.JJ=[RunOpt.JJ;JJ];
        RunOpt.KK=[RunOpt.KK;KK];j0=j0+max(II);
        RunOpt.slave=[RunOpt.slave;i4];
       end
      end
end % of loop on rigid elements
end % there are mechanical DOFs

if RunOpt.ReturnT;
 T=T';i1=find(any(T));T=T(:,i1);cdof=mdof(i1);
 if ~isempty(cGL);  T=cGL(i1,i1)'*T*cGL(i1,i1); end
 if isstruct(Up)&&nargout<2
  T=struct('Stack',[],'T',T,'DOF',cdof);
 end
 out=T;out1=cdof;
else % return equivalent constraint
 out=struct('c',sparse(RunOpt.II,RunOpt.JJ,RunOpt.KK,j0,length(mdof)), ...
     'DOF',mdof,'slave',RunOpt.slave);
end
