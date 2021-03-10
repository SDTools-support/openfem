function [out,out1,out2]=penta6(CAM,varargin);

%PENTA6	element function of a 6-node 18-DOF isoparametric solid element (3r1d)
%
%	As all element functions (see ELEM0), PENTA6 is called by FE_MK for
%	model assembly, FEPLOT for structural deformation visualization, ...
%
%	In an model description matrix a group of PENTA6 elements starts with a
%	header row [Inf  abs('penta6') 0 ] followed by element property rows
%       ELT following the format
%	    [n1 ... n6 MatId ProId]
%         with
%	   n1 ... n6  identification numbers for the element nodes
%	   MatId  material property identification number
%	   ProId  section property identification number (only used by upcom)
%
%	Note that a degenerate HEXA8 with nodes [n1 n2 n3 n4 n5 n5 n6 n6] can
%       be used instead of a PENTA6
%
%       PL material property matrix. 
%       3-D isotropic and orthotropic materials are supported.  
%       See M_ELASTIC and FE_MAT.
%
%       IL PENTA6 elements do not use section properties.
%
%       Standard tests available with penta6('testeig') (mat,eig,load) 
%
% 	See sdtweb      eltfun, elem0
%	See also help   hexa20, hexa8, tetra4


%       Etienne Balmes, Jean-Michel Leclere, Marina Vidrascu  
%       Copyright (c) 2001-2011 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

% standard calls with one input argument
if ischar(CAM)
 [CAM,Cam]=comstr(CAM,1);
 if comstr(Cam,'integinfo')

  %constit integ,elmap                 ID,pl,il
  [out,out1,out2]=p_solid('buildconstit 3 6',varargin{:});

 elseif comstr(Cam,'matcall'); [out,out1]=elem0(CAM,varargin{:});
 elseif comstr(Cam,'dofcall');out=elem0('dofcall'); % variable field elements
 elseif comstr(Cam,'call')
   out='[k1,m1]=penta6(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,elmap);';

 elseif comstr(Cam,'rhscall') % call for load assembly
   out='rhs_of';

 elseif comstr(Cam,'groupinit');out=elem0(CAM,'penta6');
 % Here one defines the constants needed for element integration
 elseif comstr(Cam,'constants');

  if nargin<3; p_solid('constsolid','penta6',[],[]);return;
  elseif varargin{2}(size(varargin{2},1),1)==-9999; % old of_mk_sub.c elements
    i2=18;i1=find(triu(ones(i2,i2)));r1=zeros(i2,i2);r1(i1)=1:length(i1);
    out=struct('ElMap',r1+tril(r1',-1));   out2=[];
    out1=varargin{1};out1(4,:)=3;% Tell of_mk('MatrixIntegration') this is 3d
  else;
   [out,i2,out2]=p_solid('constsolid','penta6',varargin{2:end});
   out1=varargin{1};out1(4,:)=i2; % Tell MatrixIntegration this is 3d 
  end
 elseif comstr(Cam,'node');   out = [1:6];
 elseif  comstr(Cam,'prop');  out = [7 8 9];
 elseif comstr(Cam,'dof')
   out =[1:6];k=[1:3]'/100;out=out(ones(3,1),:)+k(:,ones(6,1));out=out(:);
 elseif comstr(Cam,'line')
   out = [1:3 1 4:6 4 0 3 6 0 2 5];
 elseif comstr(Cam,'patch')
   out = [1 3 2 2;4:6 6;2 3 6 5;1 2 5 4;1 4 6 3];
 elseif  comstr(Cam,'parent')
   out = 'penta6';
 elseif  comstr(Cam,'edge')   
   out = [1 2; 2 3; 3 1; 1 4; 2 5; 3 6; 4 5; 5 6; 6 4]; 
 elseif  comstr(Cam,'face')   
   out = [1 3 2 2; 1 4 6 3; 1 2 5 4; 4 5 6 6; 2 3 6 5]; 
 elseif  comstr(Cam,'flip');   out=[1 3 2 4 6 5]; out1=1:6; 
 elseif comstr(Cam,'sci_face')
   out = [1 2 5 4;2 3 6 5;1 4 6 3;4 5 6 6;1 3 2 2];
 elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,5);

   model=femesh('testpenta6');
   [m,k,mdof]=fe_mknl(model,'NoT');out=stack_cell(k,m);
   disp('TestMat passed');


 elseif  comstr(Cam,'test'); [CAM,Cam] = comstr(CAM,5);

   [out,out1]=femesh(strcat(['teststruct penta6' Cam]));
 elseif  comstr(Cam,'state');out='';% State should be handled by of_mk
 elseif comstr(CAM,'cvs')
  out='$Revision: 1.32 $  $Date: 2011/05/09 13:41:47 $'; return;

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
 node=node(NNode(elt(1,1:6)),[5:7 1]);
end

if integ(point(6)+7)~=3
 typ=point(5);
 if (typ==0) % mass and stiffness
  [k1,m1]=of_mk('penta6',int32(point),integ,constit,node);
  k=reshape(k1(elmap),size(elmap,1),size(elmap,2)); m=reshape(m1(elmap),size(elmap,1),size(elmap,2));
  if k(1)<0; error('Reorient elements using feutil(''orient'')'); end
  out=k; out1=m;
 elseif typ>99
  warning('OpenFEM:DIAG','This standard call should be done in fe_mknl');
  out=of_mk('penta6',int32(point),varargin{3:4},node,varargin{[5 6 7 8]});
 else
   k1=of_mk('penta6',int32(point),integ,constit,node);
   k=reshape(k1(elmap),size(elmap,1),size(elmap,2));
   out=k; out1=[];
 end

end
% -----------------------------------------------------------------------------
% -----------------------------------------------------------------------------
% -----------------------------------------------------------------------------


