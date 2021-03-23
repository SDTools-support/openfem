function [out,out1,out2]=quadb(CAM,varargin);

%QUADB	element function of a 8-node 40/48-DOF plate/shell element
%
%	As all element functions (see ELEM0), QUADB is called by FE_MK for
%	model assembly, FEPLOT for structural deformation visualization, ...
%
%	QUADB header rows follow the format [Inf  abs('quadb') 0 EGID ...]
%	QUADB element property rows follow the format
%	    [n1 n2 n3 n4 n5 n6 n7 n8 MatId ProId EltId (Theta Zoff)]
%         with
%          n1 ... n8  identification numbers for the element nodes
%          Theta  angle of material x axis with n1-n2 line 
%          Zoff   optional off-set of the element reference plane with 
%                 respect to the element nodes in the positive z-elt direction 
%	   MatId  material property identification number
%	   ProId  element property identification number
%
%      For material see
%         M_ELASTIC elastic materials
%      For formulations see
%         P_SHELL  shell elements
%         P_SOLID  2D or 3D surface elements
%
%	   Shell formulations
%        0 (default) : preferred formulation (currently 2)
%        1 : 8 tria3 thin plate elements with condensation of central node
%        2 : isoparametric thick plate with reduced integration.
%			 For non-flat elements, formulation 1 is used.
%     Standard tests available with quadb('testeig') (mat,eig,load) 
%
%	See sdtweb     quadb, eltfun, elem0
%	See also help  quad4, tria3

%	Etienne Balmes
%       Copyright (c) 2001-2015 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

persistent TriaEltConst
%#ok<*ASGLU,*NASGU,*NOSEM>
% standard calls with one input argument
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1); ElemF='quadb'; ElemP='q8p'; 
 % Build Constit, Integ, and Elmap for later integration
 if comstr(Cam,'integinfo');
  [out,out1,out2]= ...
   p_shell('buildconstit',[varargin{1};48;8],varargin{2},varargin{3});
 elseif comstr(Cam,'cvs');
  out='$Revision: 1.37 $  $Date: 2015/03/06 11:57:14 $'; return;
 elseif comstr(Cam,'matcall'); %[out,out1]=elem0('matcall','quadb');
     [out,out1]=elem0('callquadb',varargin{:});
 elseif comstr(Cam,'call'); [out,out1]=elem0('matcall','quadb');
 elseif comstr(Cam,'rhscall') % call for load assembly
   out=quad4('rhscall');
 elseif comstr(Cam,'groupinit');out=elem0('groupinitogShell','quadb');
 elseif comstr(Cam,'constants');
   [out,out1,out2]=elem0('ConstantsShellq8p',varargin{:});
 elseif comstr(Cam,'node');     out = [1:8];
 elseif comstr(Cam,'prop');     out = [9 10 11];
 elseif comstr(Cam,'dofcall');  out=elem0('dofcall'); % variable field elements
 elseif comstr(Cam,'dof');
   m = [1:8];k=[1:6]'/100;m=m(ones(6,1),:)+k(:,ones(8,1));out=m(:); return;
 elseif comstr(Cam,'line');     out=[1 5 2 6 3 7 4 8 1];
 elseif comstr(Cam,'patch');    out = [1 5 2 6 3 7 4 8]; return;
 elseif  comstr(Cam,'edge');   out = q8p('edge');
 elseif  comstr(Cam,'face');   out = [1:8];
 elseif  comstr(Cam,'flip');   out=[1 4 3 2 8 7 6 5]; out1=1:8; 
 elseif  comstr(Cam,'sci_face'); out = [1 5 8 8;5 2 6 6;6 3 7 7;7 4 8 8;5 6 7 8];
 elseif comstr(Cam,'parent');   out = 'quadb';

 % Basic matrix test of the element - - - - - - - - - - - - - - - - - -
 elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,5); 

   model=femesh('testquadb');
   NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
   nodeE=model.Node(NNode(model.Elt(2,1:8)),[5:7 1]);
   [constit,integ,elmap]=quadb('integinfo',[100;110],model.pl,model.il);
   constit(3)=1; % use 8 triangles for reference (eig is with reduced shear)
   [k,m]=quadb(nodeE,model.Elt(2,:),[0 0 0 0 0 0 0 0 0], ...
       int32(integ),constit,elmap);
   integ(5)=4; constit(2)=.5;
   k1=quadb(nodeE,model.Elt(2,:),[0 0 0 0 0 0 0 0 0],int32(integ),constit,elmap);
   out=stack_cell(k,m);
   disp('TestMat passed');

 elseif  comstr(Cam,'test');   [CAM,Cam] = comstr(Cam,5);

  [out,out1]=femesh(strcat(['teststruct quadb' Cam]));

 % redirect the rest to elem0 - - - - - - - - - - - - - - - - - - - -
 else;
   if nargin==1; varg={ElemP}; else; varg=varargin;end
   if nargout<=1;    out=elem0(CAM,varg{:});
   else;      [out,out1]=elem0(CAM,varg{:});
   end
 end

return
end % of standard calls with one input argument

% -----------------------------------------------------------------------------
% element matrix assembly - - - - - - - - - - - - - - - - - - - - - - - - - -
node=CAM; 
elt=varargin{1}; 
point=varargin{2};
integ=varargin{3};
constit=varargin{4};constit=constit(double(point(7))+[1:size(constit,1)])';
elmap=varargin{5};
typ=point(5);
if isa(elmap,'int32'); elmap = double(elmap); end


% find properties

if size(node,1)~=8
 NNode=sparse(node(:,1),1,1:size(node,1));
 node=node(NNode(elt(1,1:8)),[5:7 1]);
end

if constit(3)~=0
elseif ~any(constit([45 46])); 
  % if no shear sdtweb p_shell('ShellConstit')
  constit(3)=1; 
else; constit(3)=2; %current best formulation
end 
% define nodes x and local basis bas

jElt=1;

if size(node,2)==10 % xyz id v1xyz v3xyz
 p=sp_util('basis',mean(node(:,5:7)),mean(node(:,8:10)));
 p=p*[1 0 0;0 0 1;0 -1 0];bas=p;
 x =node(:,1:3)*p;
else; % Orient along line
 i1 = 1:8; x=node(:,1:3);x=x-ones(8,1)*x(1,1:3);
 bas = sp_util('basis',x(2,:),x(3,:));  x=x*bas;
end


% ---------------------------------------------------------------------------
if constit(5)==0; disp('element ignored : thickness is equal to zero');

   out=[0];out11=[0];return;

% ---------------------------------------------------------------------------
elseif any(constit(3)==[0 1])  || ...  % formulation based on 8 tria3 elements
 (constit(3)==2 && max(abs(x(:,3)))>max(max(abs(x)))*sqrt(eps))  % not flat

  node(size(node,1)+1,:)=mean(node,1);
  i2=[elt(1,[9 10]) 0];
  if isempty(TriaEltConst); TriaEltConst=integrules('tria3',-1);end
  m=zeros(54,54);k=m; 
  
% [m2,k2,idof2] = fe_mk(cnode,celt,pl,il,[],idof,[0 1]);

  i3=[0 0 0;5 8 1;6 5 2;7 6 3;8 7 4;8 5 9;5 6 9;6 7 9;7 8 9];
  constit=varargin{4};

  for j1=2:9
    %[i1,k1,m1] = tria3(cnode,celt(j1,:),pl,il,[0 1 1]);
    [k1,m1] = tria3(node(i3(j1,:),:),i3(j1,:),[0 0 0 0 0 0 0 0 0],integ,constit,elmap,[],TriaEltConst,[]);
    i4=[-5:0]';i4=i4(:,[1 1 1]);i5=i3(j1,:);i5=i5(ones(6,1),:)*6+i4;i5=i5(:);
    k(i5,i5)=k(i5,i5)+k1;
    m(i5,i5)=m(i5,i5)+m1;
  end
 
  % condensation of centernode

  i1=diag(k);i1=find(i1>eps*mean(i1));
  if length(i1)==size(k,1)
   i1=1:48;i2=49:54;
   tr=[eye(length(i1),length(i1));-pinv(real(k(i2,i2)))*real(k(i2,i1))];
  else
   i2=zeros(54,1);i2(i1)=1;i1=find(i2(1:48));i2=48+find(i2(49:54));
   tr=zeros(54,48);
   tr([i1;i2],i1)=[eye(length(i1),length(i1));-real(k(i2,i2))\real(k(i2,i1))];
  end

  k=tr'*k*tr; m=tr'*m*tr;
  if     typ==0;  out=k; out1=m;
  elseif typ==1;  out=k; out1=[];
  elseif typ==2;  out=m; out1=[];
  elseif typ==100; error('Load computation implemented in quad4');
  else error('Not a supported element matrix type');
  end
 
% ---------------------------------------------------------------------------
elseif constit(3)==2

% w   : [r s t weight] quadrature rule
% xi  : isoparametric coordinates of nodes
% na  : (integration points) x (shape functions)
% nar,s : (integration points) x (partial derivatives with respect to r,s)
%[m,k,mdof] = fe_mk(FEnode,FEelt,pl,il,[],[.01 .02 .03 .04 .05]',0);


% generated with : integrules('quadb')
FEw = [ -5.7735026918962584e-001 -5.7735026918962584e-001 0.0000000000000000e+000 , ...
1.0000000000000000e+000 ;
5.7735026918962584e-001 -5.7735026918962584e-001 0.0000000000000000e+000 , ...
1.0000000000000000e+000 ;
-5.7735026918962584e-001 5.7735026918962584e-001 0.0000000000000000e+000 , ...
1.0000000000000000e+000 ;
5.7735026918962584e-001 5.7735026918962584e-001 0.0000000000000000e+000 , ...
1.0000000000000000e+000 ;
-7.7459666924148340e-001 -7.7459666924148340e-001 0.0000000000000000e+000 , ...
3.0864197530864201e-001 ;
0.0000000000000000e+000 -7.7459666924148340e-001 0.0000000000000000e+000 , ...
4.9382716049382713e-001 ;
7.7459666924148340e-001 -7.7459666924148340e-001 0.0000000000000000e+000 , ...
3.0864197530864201e-001 ;
-7.7459666924148340e-001 0.0000000000000000e+000 0.0000000000000000e+000 , ...
4.9382716049382713e-001 ;
0.0000000000000000e+000 0.0000000000000000e+000 0.0000000000000000e+000 , ...
7.9012345679012341e-001 ;
7.7459666924148340e-001 0.0000000000000000e+000 0.0000000000000000e+000 , ...
4.9382716049382713e-001 ;
-7.7459666924148340e-001 7.7459666924148340e-001 0.0000000000000000e+000 , ...
3.0864197530864201e-001 ;
0.0000000000000000e+000 7.7459666924148340e-001 0.0000000000000000e+000 , ...
4.9382716049382713e-001 ;
7.7459666924148340e-001 7.7459666924148340e-001 0.0000000000000000e+000 , ...
3.0864197530864201e-001 ;
];
FEna = [ 9.6225044864937770e-002 -1.6666666666666663e-001 -9.6225044864937589e-002 , ...
-1.6666666666666663e-001 5.2578342306320847e-001 1.4088324360345802e-001 , ...
1.4088324360345802e-001 5.2578342306320847e-001 ;
-1.6666666666666663e-001 9.6225044864937770e-002 -1.6666666666666663e-001 , ...
-9.6225044864937589e-002 5.2578342306320847e-001 5.2578342306320847e-001 , ...
1.4088324360345802e-001 1.4088324360345802e-001 ;
-1.6666666666666663e-001 -9.6225044864937589e-002 -1.6666666666666663e-001 , ...
9.6225044864937770e-002 1.4088324360345802e-001 1.4088324360345802e-001 , ...
5.2578342306320847e-001 5.2578342306320847e-001 ;
-9.6225044864937589e-002 -1.6666666666666663e-001 9.6225044864937770e-002 , ...
-1.6666666666666663e-001 1.4088324360345802e-001 5.2578342306320847e-001 , ...
5.2578342306320847e-001 1.4088324360345802e-001 ;
4.3237900077244512e-001 -9.9999999999999964e-002 -3.2379000772445002e-002 , ...
-9.9999999999999964e-002 3.5491933384829660e-001 4.5080666151703308e-002 , ...
4.5080666151703308e-002 3.5491933384829660e-001 ;
-9.9999999999999978e-002 -9.9999999999999978e-002 -9.9999999999999978e-002 , ...
-9.9999999999999978e-002 8.8729833462074170e-001 1.9999999999999996e-001 , ...
1.1270166537925830e-001 1.9999999999999996e-001 ;
-9.9999999999999964e-002 4.3237900077244512e-001 -9.9999999999999964e-002 , ...
-3.2379000772445002e-002 3.5491933384829660e-001 3.5491933384829660e-001 , ...
4.5080666151703308e-002 4.5080666151703308e-002 ;
-9.9999999999999978e-002 -9.9999999999999978e-002 -9.9999999999999978e-002 , ...
-9.9999999999999978e-002 1.9999999999999996e-001 1.1270166537925830e-001 , ...
1.9999999999999996e-001 8.8729833462074170e-001 ;
-2.5000000000000000e-001 -2.5000000000000000e-001 -2.5000000000000000e-001 , ...
-2.5000000000000000e-001 5.0000000000000000e-001 5.0000000000000000e-001 , ...
5.0000000000000000e-001 5.0000000000000000e-001 ;
-9.9999999999999978e-002 -9.9999999999999978e-002 -9.9999999999999978e-002 , ...
-9.9999999999999978e-002 1.9999999999999996e-001 8.8729833462074170e-001 , ...
1.9999999999999996e-001 1.1270166537925830e-001 ;
-9.9999999999999964e-002 -3.2379000772445002e-002 -9.9999999999999964e-002 , ...
4.3237900077244512e-001 4.5080666151703308e-002 4.5080666151703308e-002 , ...
3.5491933384829660e-001 3.5491933384829660e-001 ;
-9.9999999999999978e-002 -9.9999999999999978e-002 -9.9999999999999978e-002 , ...
-9.9999999999999978e-002 1.1270166537925830e-001 1.9999999999999996e-001 , ...
8.8729833462074170e-001 1.9999999999999996e-001 ;
-3.2379000772445002e-002 -9.9999999999999964e-002 4.3237900077244512e-001 , ...
-9.9999999999999964e-002 4.5080666151703308e-002 3.5491933384829660e-001 , ...
3.5491933384829660e-001 4.5080666151703308e-002 ;
];FEnar = [ -6.8301270189221941e-001 -2.2767090063073980e-001 -1.8301270189221933e-001 , ...
-6.1004233964073118e-002 9.1068360252295921e-001 3.3333333333333326e-001 , ...
2.4401693585629242e-001 -3.3333333333333326e-001 ;
2.2767090063073980e-001 6.8301270189221941e-001 6.1004233964073118e-002 , ...
1.8301270189221933e-001 -9.1068360252295921e-001 3.3333333333333326e-001 , ...
-2.4401693585629242e-001 -3.3333333333333326e-001 ;
-6.1004233964073118e-002 -1.8301270189221933e-001 -2.2767090063073980e-001 , ...
-6.8301270189221941e-001 2.4401693585629242e-001 3.3333333333333326e-001 , ...
9.1068360252295921e-001 -3.3333333333333326e-001 ;
1.8301270189221933e-001 6.1004233964073118e-002 6.8301270189221941e-001 , ...
2.2767090063073980e-001 -2.4401693585629242e-001 3.3333333333333326e-001 , ...
-9.1068360252295921e-001 -3.3333333333333326e-001 ;
-1.0309475019311125e+000 -3.4364916731037087e-001 -1.3094750193111251e-001 , ...
-4.3649167310370857e-002 1.3745966692414835e+000 1.9999999999999996e-001 , ...
1.7459666924148337e-001 -1.9999999999999996e-001 ;
-3.4364916731037087e-001 3.4364916731037087e-001 -4.3649167310370829e-002 , ...
4.3649167310370829e-002 0.0000000000000000e+000 1.9999999999999996e-001 , ...
0.0000000000000000e+000 -1.9999999999999996e-001 ;
3.4364916731037087e-001 1.0309475019311125e+000 4.3649167310370857e-002 , ...
1.3094750193111251e-001 -1.3745966692414835e+000 1.9999999999999996e-001 , ...
-1.7459666924148337e-001 -1.9999999999999996e-001 ;
-3.8729833462074170e-001 -3.8729833462074170e-001 -3.8729833462074170e-001 , ...
-3.8729833462074170e-001 7.7459666924148340e-001 5.0000000000000000e-001 , ...
7.7459666924148340e-001 -5.0000000000000000e-001 ;
0.0000000000000000e+000 0.0000000000000000e+000 0.0000000000000000e+000 , ...
0.0000000000000000e+000 0.0000000000000000e+000 5.0000000000000000e-001 , ...
0.0000000000000000e+000 -5.0000000000000000e-001 ;
3.8729833462074170e-001 3.8729833462074170e-001 3.8729833462074170e-001 , ...
3.8729833462074170e-001 -7.7459666924148340e-001 5.0000000000000000e-001 , ...
-7.7459666924148340e-001 -5.0000000000000000e-001 ;
-4.3649167310370857e-002 -1.3094750193111251e-001 -3.4364916731037087e-001 , ...
-1.0309475019311125e+000 1.7459666924148337e-001 1.9999999999999996e-001 , ...
1.3745966692414835e+000 -1.9999999999999996e-001 ;
4.3649167310370829e-002 -4.3649167310370829e-002 3.4364916731037087e-001 , ...
-3.4364916731037087e-001 0.0000000000000000e+000 1.9999999999999996e-001 , ...
0.0000000000000000e+000 -1.9999999999999996e-001 ;
1.3094750193111251e-001 4.3649167310370857e-002 1.0309475019311125e+000 , ...
3.4364916731037087e-001 -1.7459666924148337e-001 1.9999999999999996e-001 , ...
-1.3745966692414835e+000 -1.9999999999999996e-001 ;
];FEnas = [ -6.8301270189221941e-001 -6.1004233964073118e-002 -1.8301270189221933e-001 , ...
-2.2767090063073980e-001 -3.3333333333333326e-001 2.4401693585629242e-001 , ...
3.3333333333333326e-001 9.1068360252295921e-001 ;
-6.1004233964073118e-002 -6.8301270189221941e-001 -2.2767090063073980e-001 , ...
-1.8301270189221933e-001 -3.3333333333333326e-001 9.1068360252295921e-001 , ...
3.3333333333333326e-001 2.4401693585629242e-001 ;
2.2767090063073980e-001 1.8301270189221933e-001 6.1004233964073118e-002 , ...
6.8301270189221941e-001 -3.3333333333333326e-001 -2.4401693585629242e-001 , ...
3.3333333333333326e-001 -9.1068360252295921e-001 ;
1.8301270189221933e-001 2.2767090063073980e-001 6.8301270189221941e-001 , ...
6.1004233964073118e-002 -3.3333333333333326e-001 -9.1068360252295921e-001 , ...
3.3333333333333326e-001 -2.4401693585629242e-001 ;
-1.0309475019311125e+000 -4.3649167310370857e-002 -1.3094750193111251e-001 , ...
-3.4364916731037087e-001 -1.9999999999999996e-001 1.7459666924148337e-001 , ...
1.9999999999999996e-001 1.3745966692414835e+000 ;
-3.8729833462074170e-001 -3.8729833462074170e-001 -3.8729833462074170e-001 , ...
-3.8729833462074170e-001 -5.0000000000000000e-001 7.7459666924148340e-001 , ...
5.0000000000000000e-001 7.7459666924148340e-001 ;
-4.3649167310370857e-002 -1.0309475019311125e+000 -3.4364916731037087e-001 , ...
-1.3094750193111251e-001 -1.9999999999999996e-001 1.3745966692414835e+000 , ...
1.9999999999999996e-001 1.7459666924148337e-001 ;
-3.4364916731037087e-001 4.3649167310370829e-002 -4.3649167310370829e-002 , ...
3.4364916731037087e-001 -1.9999999999999996e-001 0.0000000000000000e+000 , ...
1.9999999999999996e-001 0.0000000000000000e+000 ;
0.0000000000000000e+000 0.0000000000000000e+000 0.0000000000000000e+000 , ...
0.0000000000000000e+000 -5.0000000000000000e-001 0.0000000000000000e+000 , ...
5.0000000000000000e-001 0.0000000000000000e+000 ;
4.3649167310370829e-002 -3.4364916731037087e-001 3.4364916731037087e-001 , ...
-4.3649167310370829e-002 -1.9999999999999996e-001 0.0000000000000000e+000 , ...
1.9999999999999996e-001 0.0000000000000000e+000 ;
3.4364916731037087e-001 1.3094750193111251e-001 4.3649167310370857e-002 , ...
1.0309475019311125e+000 -1.9999999999999996e-001 -1.7459666924148337e-001 , ...
1.9999999999999996e-001 -1.3745966692414835e+000 ;
3.8729833462074170e-001 3.8729833462074170e-001 3.8729833462074170e-001 , ...
3.8729833462074170e-001 -5.0000000000000000e-001 -7.7459666924148340e-001 , ...
5.0000000000000000e-001 -7.7459666924148340e-001 ;
1.3094750193111251e-001 3.4364916731037087e-001 1.0309475019311125e+000 , ...
4.3649167310370857e-002 -1.9999999999999996e-001 -1.3745966692414835e+000 , ...
1.9999999999999996e-001 -1.7459666924148337e-001 ;
];


 % going back to local gradient information
  
  xr = FEnar*x(:,1); xs = FEnas*x(:,1);
  yr = FEnar*x(:,2); ys = FEnas*x(:,2);

  jdet = xr.*ys-xs.*yr;  jdet=jdet*sign(jdet(1));
  if any(jdet<0); disp('quad4 negative Jacobian');end

  i2=ones(size(FEna,2),1); % number of shape functions
  nax = [ ys(:,i2).*FEnar-yr(:,i2).*FEnas]; %true divided by jdet
  nay = [-xs(:,i2).*FEnar+xr(:,i2).*FEnas];

 % stiffness matrix assembly - - - - - - - - - - - - - - - - - - - - - - - -
 if any([1 0]==typ);

  k = spalloc(48,48,48);
  BM=zeros(3,48);%spalloc(3,48,16);
  BB=zeros(3,48);%spalloc(3,48,16);
  BS=zeros(2,48);%spalloc(2,48,16);

  db=constit([30 36 42;31 37 43;32 38 44]);
  dm=constit([9 15 21;10 16 22;11 17 23]);
  ds=constit([45 47;46 48]);

  % see Hugues p. 150 and 322
  for j1=5:13
     BM(1,1:6:48)=          nax(j1,:)            ;
     BM(2,2:6:48) =                    nay(j1,:) ;
     BM(3,[1:6:48 2:6:48])=[nay(j1,:)  nax(j1,:)];

     BB(1, 5:6:48) =                  -nax(j1,:) ;
     BB(2, 4:6:48) =        nay(j1,:)            ;
     BB(3,[4:6:48 5:6:48])=[nax(j1,:) -nay(j1,:)];
     k = k + BM'*(dm*FEw(j1,4)/jdet(j1))*BM ...
           + BB'*(db*FEw(j1,4)/jdet(j1))*BB;
  end % loop on integration points

  for j1=1:4 %5:13 % for reduced integration use j1=1:4;
     BS(1,[3:6:48 5:6:48]) = [nax(j1,:)/jdet(j1)             FEna(j1,:)];
     BS(2,[3:6:48 4:6:48]) = [nay(j1,:)/jdet(j1) -FEna(j1,:)           ];
     k = k + BS'*(ds*FEw(j1,4)*jdet(j1))*BS;
  end % loop on integration points

  if constit(4)~=-1  % drilling dof xxx some fixing needed
   ind  = [6 12 18 24]; if constit(4)==0; constit(4)=1; end
   i1=diag(k);k(ind,ind)=constit(4)*1e-6*mean(i1([4 5 10 11 16 17])) ...
                  *([1 -.5 0 -.5;-.5 1 -.5 0;0 -.5 1 -.5;-.5 0 -.5 1]);
   ind  = [30 36 42 48]; 
   k(ind,ind)=constit(4)*1e-6*mean(i1([4 5 10 11 16 17])) ...
                  *([1 -.5 0 -.5;-.5 1 -.5 0;0 -.5 1 -.5;-.5 0 -.5 1]);
  end

  out=of_mk('xkx_trans',bas,full(k)); out1=[];

 end;

if any([0 2]==typ); % mass matrix assembly

  m = spalloc(48,48,48);  B=zeros(3,48);%spalloc(3,48,12);
  r1 = constit(1); 
  %xxx if length(il)>11 & il(11) r1=r1+il(11);end 
  for j1=5:13
     B(1,1:6:48)= FEna(j1,:);
     B(2,2:6:48)= FEna(j1,:);
     B(3,3:6:48)= FEna(j1,:);
     m = m + (r1*FEw(j1,4)*jdet(j1)*B')*B;
  end % loop on integration points
  m=of_mk('xkx_trans',bas,full(m));

  if typ==2; out=m; out1=[];  
  else out1=m; end

 end % mass matrix


% ---------------------------------------------------------------------------
else error('unknown formulation (%i)',constit(3)); 
end




