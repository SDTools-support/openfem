function [out,out1,out2]=tetra10(CAM,varargin);

%TETRA10 element function for a 10-node 12-DOF isoparametric solid element(3p2d)
%
%	As all element functions (see ELEM0), TETRA10 is called by FE_MK for
%	model assembly, FEPLOT for structural deformation visualization, ...
%
%	In an model description matrix a group of TETRA10 elements starts with 
%	a header row [Inf  abs('tetra10') 0 ...] followed by element property 
%	rows ELT following the format
%	    [n1 ... n10 MatId ProId EltId]
%         with
%          n1 ... n10  identification numbers for the element nodes
%	   MatId  material property identification number
%	   ProId  element property identification number
%
%       PL material property matrix. 
%       3-D isotropic and orthotropic materials are supported.  
%       See M_ELASTIC and FE_MAT.
%
%       IL TETRA10 elements do not use element properties.
%
%       Standard tests available with tetra10('testeig') (mat,eig,load) 
%
% 	See sdtweb     eltfun, elem0
%	See also help  hexa20, hexa8, penta6


%       Etienne Balmes, Jean-Michel Leclere, Marina Vidrascu  
%       Copyright (c) 2001-2011 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

if comstr(CAM,'cvs')
 out='$Revision: 1.33 $  $Date: 2011/05/09 15:16:32 $'; return;
end
% standard calls with one input argument
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1);
 % Build Constit, Integ, and Elmap for later integration
 if comstr(Cam,'integinfo')

  %constit integ,elmap                 ID,pl,il
  [out,out1,out2]=p_solid('buildconstit 3 10',varargin{:});

 elseif comstr(Cam,'matcall'); [out,out1]=elem0(CAM,varargin{:});
 elseif comstr(Cam,'dofcall');out=elem0('dofcall'); % variable field elements

 elseif comstr(Cam,'call')
   out='be=tetra10(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,gstate,defe,EltConst,InfoAtNode);';

 elseif comstr(Cam,'rhscall') % call for load assembly
   out='rhs_of';

 elseif comstr(Cam,'groupinit');out=elem0(CAM,'tetra10');
 % Here one defines the constants needed for element integration
 elseif comstr(Cam,'constants');

  if nargin<3; p_solid('constsolid','tetra10',[],[]);return;
  elseif varargin{2}(size(varargin{2},1),1)==-9999; % old of_mk_sub.c elements
    i2=30;i1=find(triu(ones(i2,i2)));r1=zeros(i2,i2);r1(i1)=1:length(i1);
    out=struct('ElMap',r1+tril(r1',-1));   out2=[];
    out1=varargin{1};out1(4,:)=3;% Tell of_mk('MatrixIntegration') this is 3d
  else;
   [out,i2,out2]=p_solid('constsolid','tetra10',varargin{2:end});
   out1=varargin{1};out1(4,:)=i2; % Tell MatrixIntegration this is 3d 
  end

 elseif comstr(Cam,'node');   out = [1:10];
 elseif comstr(Cam,'prop');   out = [11:13];
 elseif comstr(Cam,'dof')
   out=[1:10];i1=[.01;.02;.03];
   out=out([1 1 1],:)+i1(:,ones(10,1));
   out=out(:);
 elseif comstr(Cam,'line');    out = [1 5 2 6 3 7 1 8 4 9 2 0 4 10 3];
 elseif comstr(Cam,'patch')
  out=[1 7 3 6 2 5;1 8 4 10 3 7;1 5 2 9 4 8;2 6 3 10 4 9];
 elseif  comstr(Cam,'parent'); out = 'tetra10';
 elseif  comstr(Cam,'edge')  
	 out = [1 2 5; 2 3 6; 3 1 7; 1 4 8; 2 4 9; 3 4 10];	 
 elseif  comstr(Cam,'face') 
	 out = [1 3 2 7 6 5; ...
	 				 1 4 3 8 10 7; ...
					 1 2 4 5 9 8; ...
					 2 3 4 6 10 9];  
 elseif  comstr(Cam,'flip'); out=[1 3 2 4 7 6 5 8 10 9];out1=[1:10]; % volume flip
 elseif  comstr(Cam,'sci_face');
          out = [1 5 8;5 2 9;9 4 8;5 9 8;
			2 6 9;6 3 10;10 4 9;6 10 9;
			1 8 7;8 4 10;10 3 7;8 10 7;
			2 5 6;5 1 7;7 3 6;5 7 6];

 elseif  comstr(Cam,'state');out='';% State should be handled by of_mk
 elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,5);

   model=femesh('testtetra10');
   [m,k,mdof]=fe_mknl(model,'NoT');out=stack_cell(k,m);
   disp('TestMat passed');

 elseif  comstr(Cam,'test'); [CAM,Cam] = comstr(CAM,5);

  [out,out1]=femesh(strcat(['teststruct tetra10' Cam]));

 % Basic tests of the element - - - - - - - - - - - - - - - - - - - - - -
 else sdtw('''%s'' unknown',CAM);  
 end

return
end % of standard calls with one input argument

% -----------------------------------------------------------------------------
% element matrix assembly - - - - - - - - - - - - - - - - - - - - - - - - - -

node=CAM; 
elt=varargin{1}; 
point=varargin{2};
integ=varargin{3};
constit=varargin{4};
elmap=varargin{5};
if isa(elmap,'int32'); elmap = double(elmap); end


if size(node,2)~=4
 NNode(node(:,1))=1:size(node,1);
 node=node(NNode(elt(1,1:10)),[5:7 1]);
end

if integ(point(6)+7)~=3
 typ=point(5);
 if (typ==0) % mass and stiffness
  [k1,m1]=of_mk('tetra10',int32(point),integ,constit,node);
  k=reshape(k1(elmap),size(elmap,1),size(elmap,2)); m=reshape(m1(elmap),size(elmap,1),size(elmap,2));
  if k(1)<0; error('Reorient elements using feutil(''orient'')'); end
  out=k; out1=m;
 elseif typ>99
  warning('OpenFEM:DIAG','This standard call should be done in fe_mknl');
  out=of_mk('tetra10',int32(point),varargin{3:4},node,varargin{[5 6 7 8]});
 else
   k1=of_mk('tetra10',int32(point),integ,constit,node);
   k=reshape(k1(elmap),size(elmap,1),size(elmap,2));
   out=k; out1=[];
 end


 
% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------
end


