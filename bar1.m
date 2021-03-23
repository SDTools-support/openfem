function [out,out1,out2]=bar1(CAM,varargin)

%BAR1 element function for the 2-node 6-DOF axial bar element
%
%	As all element functions (see ELEM0), BAR1 is called by FE_MK for
%	model assembly, FEPLOT for structural deformation visualization, ...
%
%	In an model description matrix a group of BAR1 elements starts with a
%	header row [Inf  abs('bar1') 0 ...] followed by element property rows
%       ELT following the format
%	    [n1 n2 MatId ProId EltId]
%         with
%	   n1,n2  node numbers of the tips of the bar element
%	   MatId  material property identification number
%	   ProId  element property identification number
%	   EltId  element identifier (optional, default 0)
%
%       PL material property matrix. Isotropic materials [MatId 1  E  nu rho G]
%	   (See FE_MAT) are the only supported.
%
%       IL element property matrix. Rows associated to BEAM1 elements follow
%	   the format
%	    [ProId Type 0 0 0  A]
%          with
%	     Type : 1 for standard beam definition (no other type supported)
%            A    : section area
%
%	See sdtweb      bar1, eltfun, elem0
%	See also help   beam1, elem0

%       Etienne Balmes
%       Copyright (c) 2001-2010 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

%#ok<*ASGLU,*NASGU>
% standard calls with one input argument
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1);
 % Build Constit, Integ, and Elmap for later integration
 if comstr(Cam,'integinfo')
  ElemF='bar1'; % Bypass because bars can have p_beam properties
  [out,out1,out2]= ...
     p_solid('buildconstit 1 2',[varargin{1}],varargin{2:end});
 elseif comstr(Cam,'matcall') 
     constit=varargin{2};
     if ~isempty(constit)&&constit(1)<0;[out,out1]=elem0(CAM,varargin{:});
     else; out=bar1('call'); out1=0; % CallSymFlag
     end
 elseif comstr(Cam,'call')  % call for matrix assembly
   out='[k1,m1]=bar1(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,elmap);';
 elseif comstr(Cam,'rhscall') % call for load assembly
   out='be=bar1(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,elmap,defe);';
   %out='be=bar1(DofPos,NodePos,node,pointers,integ,constit,elmap,InfoAtNode,EltConst,defe,out);'

 elseif comstr(Cam,'groupinit');out=elem0(CAM,'bar1');
 elseif comstr(Cam,'constants')
  if nargin<4; p_solid('constsolid','bar1',[],[]);return;
  elseif any(varargin{3}(1)==[-1 -3]) % [-1 line elt, -3 surface elt]
   [out,i2,out2]=p_solid('constsolid','bar1',varargin{2:end});
   out1=varargin{1};out1(4,:)=i2; % Tell MatrixIntegration this is 3d 
  else; out2=[]; out=[];out1=varargin{1};out1(4,:)=13;
  end
 elseif  comstr(Cam,'node');  out = [1 2];
 elseif  comstr(Cam,'prop');  out = [3 4 5];
 elseif comstr(Cam,'dofcall');out=elem0('dofcall'); % variable field elements
 elseif  comstr(Cam,'dof');   out=[1.01 1.02 1.03 2.01 2.02 2.03]';
 elseif  comstr(Cam,'line');   out = [1 2];
 elseif  comstr(Cam,'face');   out =[];
 elseif  comstr(Cam,'sci_face'); out = [1 2 2];
 elseif  comstr(Cam,'edge');   out =[1 2];
 elseif  comstr(Cam,'patch');  out = [1 2];
 elseif  comstr(Cam,'parent');   out = 'beam1';
 elseif comstr(Cam,'state');  out='';% call for state update

 % Basic matrix test of the element - - - - - - - - - - - - - - - - - -
 elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,5);

   model=femesh('testbar1 divide 1');
   [constit,integ,elmap]=bar1('integinfo',[100;112],model.pl,model.il);
   [k,m]=bar1(model.Node,model.Elt(2,:),[0 0 0 0 0 0 0 0 0],int32(integ),constit,elmap);
   integ(5)=4; constit(2)=.5;
   k1=bar1(model.Node,model.Elt(2,:),[0 0 0 0 0 0 0 0 0],int32(integ),constit,elmap);

   out=stack_cell(k,m);
   if nargout==0;disp('TestMat passed');end

 elseif  comstr(Cam,'test')

   [out,out1]=femesh(strcat(['teststruct bar1' Cam]));
   C1=fe_stress('stress-gstate',out,out1);data=C1.GroupInfo{1,5};
   
 elseif comstr(CAM,'cvs')
  out='$Revision: 1.36 $  $Date: 2016/06/24 06:49:33 $'; return;

 end
return
end % of standard calls with one input argument

% element matrix assembly - - - - - - - - - - - - - - - - - - - - - - - - -

node=CAM; 
elt=varargin{1}; 
point=varargin{2};if strncmp(version,'6',1);point=double(point);end
integ=varargin{3}; 
constit=varargin{4};
elmap=varargin{5};
typ=point(5);

if size(node,2)~=4
 NNode(node(:,1))=1:size(node,1);
 node=node(NNode(elt(1,1:2)),[5:7 1]);
end

x = [find(node(:,4)==elt(1,1));find(node(:,4)==elt(1,2))];
x = node(x,1:3); x0 = find(node(:,4)==elt(1,5));

if isempty(x0); x0=[1.5 1.5 1.5]; else; x0 = node(x0,1:3); end
l = norm(x(2,:)-x(1,:)); l2 = l^2;
x(1,:) = x(2,:)-x(1,:);x(2,:) = x0-x(2,:);
x=basis(x(1,:),x(2,:),1);

% standard bar element
if any([0 1 5]==typ)	% stiffness
  k = zeros(6,6);ind = [1 4];k(ind,ind) = constit(1+point(7))/l*[1 -1;-1 1];
  % coordinate transformation - - - - - - - - - - - - - - - - - - -
  x=x';  Coo = zeros(6,6); for j1 = 1:3:6; Coo(j1+(0:2),j1+(0:2)) = x; end
  out = Coo'*k*Coo; out1=[]; 
end
if any([0 2]==typ)  % mass
   m = zeros(6,6);
   m([1 8 15 22 29 36])  = constit(3+point(7))*l/3;
   m([4 11 18 19 26 33]) = constit(3+point(7))*l/6;

   % Lumped mass (sum over elements)
   if length(constit)>=point(7)+5 &&constit(point(7)+5)==1
    m=diag(sum(m)); 
   end
   if typ==2; out=m; out1=[];  
   else; out1=m; end
end
if any(3==typ)	% vscous damping
  out = zeros(6,6); out1=[];

end

if any(100==typ) % volumic load - - - - - - - - - - - - - - - - - -

  EltConst=evalin('caller','EltConst');
  if isempty(EltConst)
   EltConst=integrules('bar1',3); assignin('caller','EltConst',EltConst);
  end
  dire=varargin{6}';
  EltConst.nodeE=node;of_mk('buildndn',13,EltConst); b1=zeros(2,3);
  for jw=1:size(EltConst.w,1) % loop on integration points
    F=EltConst.N(jw,:)*dire;
    b1=b1+EltConst.jdet(jw)*EltConst.w(jw,4)* ...
      constit(point(7)+4)*EltConst.N(jw,:)'*F;
  end % loop on integration points
  b1=b1';out=b1(:);
elseif typ==-1
  EL=constit(1+point(7))/constit(4+point(7))/l;
  out=[x(:,1);-x(:,1)]*EL;return;
end


if ~any([0 1 2 3 5 100]==typ)
 error('Matrix type %i not supported by bar1',typ)
end

%------------------------------------------------------------------
