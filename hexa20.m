function [out,out1,out2]=hexa20(CAM,varargin);

%HEXA20 20-node 60-DOF isoparametric solid element (3q2r)
%
%	As all element functions (see ELEM0), HEXA20 is called by FE_MK for
%	model assembly, FEPLOT for structural deformation visualization, ...
%
%	In an model description matrix a group of HEXA20 elements starts with a
%	header row [Inf  abs('hexa20') 0 ...] followed by element property rows
%       ELT following the format
%	    [n1 ... n20 MatId ProId EltId]
%         with
%	   n1 ... n20  identification numbers for the element nodes
%	   MatId  material property identification number
%	   ProId  element property identification number (only used by upcom)
%	   EltId  optional element identifier
%
%       PL material property matrix. 
%       3-D isotropic and orthotropic materials are supported.  
%       See M_ELASTIC and FE_MAT.
% 
%       Standard tests available with hexa20('testeig') (mat,eig,load) 
%
%       IL HEXA20 elements do not use element properties.
%
%  See <a href="matlab: sdtweb _taglist hexa20">TagList</a>
% 	See sdtweb     eltfun, elem0
%	 See also help  hexa8, penta6, tetra4

%       Etienne Balmes, Jean-Michel Leclere, Marina Vidrascu  
%       Copyright (c) 2001-2025 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

if comstr(CAM,'cvs')
 out='$Revision: 1.36 $  $Date: 2025/04/07 17:07:43 $'; return;
end

% standard calls with one input argument
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1);
 % Build Constit, Integ, and Elmap for later integration
 if comstr(Cam,'integinfo')
  %% IntegInfo
  %constit integ,elmap                 ID,pl,il
  [out,out1,out2]=p_solid('buildconstit 3 20',varargin{:});

 elseif comstr(Cam,'matcall'); [out,out1]=elem0(CAM,varargin{:});
 elseif comstr(Cam,'dofcall');out=elem0('dofcall'); % variable field elements
 elseif comstr(Cam,'call')
   out='[k1,m1]=hexa20(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,elmap);';

 elseif comstr(Cam,'rhscall') % call for load assembly
   out='rhs_of';

 elseif comstr(Cam,'groupinit');out=elem0(CAM,'hexa20');
 % Here one defines the constants needed for element integration
 elseif comstr(Cam,'constants');

  if nargin<3; p_solid('constsolid','hexa20',[],[]);return;
  elseif varargin{2}(size(varargin{2},1),1)==-9999; % old of_mk_sub.c elements
    i2=60;i1=find(triu(ones(i2,i2)));r1=zeros(i2,i2);r1(i1)=1:length(i1);
    out=struct('ElMap',r1+tril(r1',-1));   out2=[];
    out1=varargin{1};out1(4,:)=3;% Tell of_mk('MatrixIntegration') this is 3d
  else;
   pointers=varargin{1};
   [out,i2,out2]=p_solid('constsolid','hexa20',varargin{2:end});
   if length(i2)==1;pointers(4,:)=i2; end; out1=pointers;% Tell MatrixIntegration this is 3d 
  end

 elseif comstr(Cam,'node');  out = [1:20];
 elseif  comstr(Cam,'prop');  out = [21 22 23];
 elseif comstr(Cam,'dof')
   out =[1:20];k=[1:3]'/100;
   out=out(ones(3,1),:)+k(:,ones(20,1));out=out(:);
 elseif comstr(Cam,'line')
   out = ...
    [1 9 2 10 3 11 4 12 1 13 5 17 6 18 7 19 8 20 5 0 2 14 6 0 3 15 7 0 4 16 8];
 elseif comstr(Cam,'patch')
   out = [1 12 4 11 3 10 2 9;1 13 5 20 8 16 4 12;1 9 2 14 6 17 5 13;
          5 17 6 18 7 19 8 20;2 10 3 15 7 18 6 14;3 11 4 16 8 19 7 15];
 elseif  comstr(Cam,'edge');   out = [1 2 9; 2 3 10; 3 4 11; ...
				4 1 12; 1 5 13; 2 6 14; ...
				3 7 15; 4 8 16; 5 6 17; ...
				6 7 18; 7 8 19; 8 5 20];
 elseif  comstr(Cam,'face');   out = [1 4 3 2 12 11 10 9; ...
				1 5 8 4 13 20 16 12; ...
				1 2 6 5  9 14 17 13; ...
				5 6 7 8 17 18 19 20; ...
				2 3 7 6 10 15 18 14; ...
				4 8 7 3 16 19 15 11]; %3 4 8 7 11 16 19 15];
elseif  comstr(Cam,'flip');   out=[5:8 1:4 17:20 13:16 9:12]; 
                              out1=1:20; 
 elseif comstr(Cam,'sci_face');
   out = [1 9 13 13;9 2 14 14;14 6 17 17;17 5 13 13;9 14 17 13;
		2 10 14 14;10 3 15 15;15 7 18 18;18 6 14 14;10 15 18 14;
		6 18 17 17;18 7 19 19;19 8 20 20;20 5 17 17;18 19 20 17;
		1 12 9 9;12 4 11 11;11 3 10 10;10 2 9 9;12 11 10 9;
		4 16 11 11;16 8 19 19;19 7 15 15;15 3 11 11;16 19 15 11;
		1 13 12 12;13 5 20 20;20 8 16 16;16 4 12 12;13 20 16 12];
 elseif comstr(Cam,'parent');out = 'hexa20';
 elseif  comstr(Cam,'state');out='';% State should be handled by of_mk

 elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,5);

   model=femesh('testhexa20');
   [m,k,mdof]=fe_mknl(model,'NoT');out=stack_cell(k,m);
   disp('TestMat passed');

 % Basic tests of the element - - - - - - - - - - - - - - - - - - - - - -
 elseif  comstr(Cam,'test'); [CAM,Cam] = comstr(CAM,5);

   [out,out1]=femesh(strcat(['teststruct hexa20' Cam]));

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
iopt=int32(varargin{3});
constit=varargin{4};
elmap=varargin{5};
if isa(elmap,'int32'); elmap = double(elmap); end

if size(node,2)~=4
 NNode(node(:,1))=1:size(node,1);
 node=node(NNode(elt(1,1:20)),[5:7 1]);
end

if iopt(point(6)+7)~=3
 typ=point(5);

 if (typ==0) % mass and stiffness
  [k1,m1]=of_mk('hexa20',int32(point),iopt,constit,node);
  k=reshape(k1(elmap),size(elmap,1),size(elmap,2)); m=reshape(m1(elmap),size(elmap,1),size(elmap,2));
  if k(1)<0; error('Reorient elements using feutil(''orient'')'); end
  out=k; out1=m;
 elseif typ>99
  warning('OpenFEM:DIAG','This standard call should be done in fe_mknl');
  out=of_mk('hexa20',int32(point),varargin{3:4},node,varargin{[5 6 7 8]});
 else
   k1=of_mk('hexa20',int32(point),iopt,constit,node);
   k=reshape(k1(elmap),size(elmap,1),size(elmap,2));m=[];
   out=k; out1=[];
 end

end
% -----------------------------------------------------------------------------
% -----------------------------------------------------------------------------
% -----------------------------------------------------------------------------

