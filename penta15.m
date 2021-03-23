function [out,out1,out2]=penta15(CAM,varargin);

%PENTA15 element function for the 15-node 45-DOF pentaedral 3D element (3r2d)
%
%
%	As all element functions (see ELEM0), PENTA15 is called by FE_MK for
%	model assembly, FEPLOT for structural deformation visualization, ...
%
%	In an model description matrix a group of PENTA15 elements starts with a
%	header row [Inf  abs('penta15') 0 ] followed by element property rows
%       ELT following the format
%	    [n1 ... n15 MatId ProId]
%         with
%	   n1 ... n15  identification numbers for the element nodes
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
%       IL PENTA15 elements do not use section properties.
%
%       Standard tests available with penta15('testeig') (mat,eig,load) 
%
% 	See sdtweb      eltfun, elem0
%	See also help   hexa20, hexa8, tetra4, penta6


%       Etienne Balmes, Jean-Michel Leclere, Marina Vidrascu  
%       Copyright (c) 2001-2011 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

% standard calls with one input argument
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1);
 % Build Constit, Integ, and Elmap for later integration
 if comstr(Cam,'integinfo')

  %constit integ,elmap                 ID,pl,il
  [out,out1,out2]= ...
   p_solid('buildconstit 3 15',varargin{:});

 elseif comstr(Cam,'matcall'); [out,out1]=elem0(CAM,varargin{:});
 elseif comstr(Cam,'dofcall');out=elem0('dofcall'); % variable field elements
 elseif comstr(Cam,'call')
   out='[k1,m1]=penta15(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,elmap);';

 elseif comstr(Cam,'rhscall') % call for load assembly
   out='rhs_of';

 elseif comstr(Cam,'groupinit');out=elem0(CAM,'penta15');
 % Here one defines the constants needed for element integration
 elseif comstr(Cam,'constants');

  if nargin<3; p_solid('constsolid','penta15',[],[]);return;
  elseif varargin{2}(size(varargin{2},1),1)==-9999; % old of_mk_sub.c elements
    i2=45;i1=find(triu(ones(i2,i2)));r1=zeros(i2,i2);r1(i1)=1:length(i1);
    out=struct('ElMap',r1+tril(r1',-1));   out2=[];
    out1=varargin{1};out1(4,:)=3;% Tell of_mk('MatrixIntegration') this is 3d
  else;
   [out,i2,out2]=p_solid('constsolid','penta15',varargin{2:end});
   out1=varargin{1};out1(4,:)=i2; % Tell MatrixIntegration this is 3d 
  end

 elseif  comstr(Cam,'node');  out = [1:15];
 elseif  comstr(Cam,'prop');  out = [16 17 18];
 elseif  comstr(Cam,'dof')
    out=[1:15];k=[1;2;3]/100;out=out(ones(3,1),:)+k(:,ones(15,1));
    out=out(:);
 elseif  comstr(Cam,'line');   out = [1 7 2 8 3 9 1 0 1 10 4 0 ...
                                       2 11 5 0 3 12 6 0 4 13 5 14 6 15 4];
 elseif  comstr(Cam,'patch');  
             out = [1 9 3 8 2 7 7 7; 1 10 4 15  6 12 3 9; 
                    1 7 2 11 5 13 4 10; 4 13 5 14 6 15 15 15; 
                    2 8 3 12 6 14 5 11];  
elseif  comstr(Cam,'edge');   
     out = [1 2 7; 2 3 8; 3 1 9; ...
             1 4 10; 2 5 11; 3 6 12; ...
             4 5 13; 5 6 14; 6 4 15];
elseif  comstr(Cam,'face');   
          out = [1 3 2 9 8 7 7 7; ...
                  1 4 6 3 10 15 12 9; ...
                  1 2 5 4 7 11 13 10; ...
                  4 5 6 13 14 15 15 15; ...
                  2 3 6 5 8 12 14 11];
elseif  comstr(Cam,'flip');   out=[1 3 2 4 6 5 9 8 7 10 12 11 15 14 13]; 
                              out1=1:15; 
elseif comstr(Cam,'sci_face');
          out = [1 7 10 10; 7 2 11 11; 11 5 13 13;13 4 10 10;7 11 13 10;
			2 8 11 11;8 3 12 12;12 6 14 14;14 5 11 11;8 12 14 11;
			1 10 9 9;10 4 15 15;15 6 12 12;12 3 9 9;9 10 15 12;
			5 14 13 13;14 6 15 15;15 4 13 13;13 14 15 15;
			2 7 8 8;7 1 9 9;9 3 8 8;8 7 9 9];
 elseif  comstr(Cam,'parent');   out = 'penta15'; 
 elseif  comstr(Cam,'state');out='';% State should be handled by of_mk

 % Basic matrix test of the element - - - - - - - - - - - - - - - - - -
 elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,5);

   model=femesh('testpenta15');
   [m,k,mdof]=fe_mknl(model,'NoT');out=stack_cell(k,m);
   disp('TestMat passed');

 elseif  comstr(Cam,'test'); [CAM,Cam] = comstr(CAM,5);

   [out,out1]=femesh(strcat(['teststruct penta15' Cam]));
 elseif comstr(CAM,'cvs')
  out='$Revision: 1.35 $  $Date: 2011/05/09 13:41:47 $'; return;

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
 node=node(NNode(elt(1,1:15)),[5:7 1]);
end

if integ(point(6)+7)~=3
 typ=point(5);
 if (typ==0) % mass and stiffness
  [k1,m1]=of_mk('penta15',int32(point),integ,constit,node);
  k=reshape(k1(elmap),size(elmap,1),size(elmap,2)); m=reshape(m1(elmap),size(elmap,1),size(elmap,2));
  if k(1)<0; error('Reorient elements using feutil(''orient'')'); end
  out=k; out1=m;
 elseif typ>99
  warning('OpenFEM:DIAG','This standard call should be done in fe_mknl');
  out=of_mk('penta15',int32(point),varargin{3:4},node,varargin{[5 6 7 8]});
 else
   k1=of_mk('penta15',int32(point),integ,constit,node);
   k=reshape(k1(elmap),size(elmap,1),size(elmap,2));
   out=k; out1=[];
 end


% -----------------------------------------------------------------------------
% -----------------------------------------------------------------------------
% -----------------------------------------------------------------------------
end


