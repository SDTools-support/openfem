function [out,out1,out2]=hexa8(CAM,varargin)

%HEXA8	8-node 24-DOF isoparametric volume
%
%	In an model description matrix a group of HEXA8 elements starts with a
%	header row [Inf  abs('hexa8') 0 ...] followed by element property rows
%       ELT following the format
%	    [n1 ... n8 MatId ProId EltId]
%         with
%	   n1 ... n8  identification numbers for the element nodes
%		Repeated nodes forming a degenerate elements are accepted
%	        PENTA6 is then called  
%	   MatId  material property identification number
%	   ProId  element property identification number (only used by upcom)
%	   EltId  optional element identifier
%
%       PL material property matrix. 
%       3-D isotropic and orthotropic materials are supported.  
%       See M_ELASTIC and FE_MAT.
%
%       Standard tests available with hexa8('testeig') (mat,eig,load) 
%
% 	See sdtweb     eltfun, elem0
%	See also help  hexa20, penta6, tetra4

%       Etienne Balmes, Jean-Michel Leclere, Marina Vidrascu  
%       Copyright (c) 2001-2015 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

%#ok<*ASGLU>

if comstr(CAM,'cvs')
 out='$Revision: 1.44 $  $Date: 2019/05/14 17:26:29 $'; return;
end
% standard calls with one input argument
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1);
 % Build Constit, Integ, and Elmap for later integration
 if comstr(Cam,'integinfo')

  %constit integ,elmap                 ID,pl,il
  [out,out1,out2]=p_solid('buildconstit 3 8',varargin{:});

 elseif comstr(Cam,'matcall'); [out,out1]=elem0(CAM,varargin{:});
 elseif comstr(Cam,'dofcall');out=elem0('dofcall'); % variable field elements
 elseif comstr(Cam,'call')  % call for matrix assembly
   out='[k1,m1]=hexa8(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,elmap);';
 elseif comstr(Cam,'rhscall') % call for load assembly
   out='rhs_of';
 elseif comstr(Cam,'groupinit');out=elem0(CAM,'hexa8');
 % Here one defines the constants needed for element integration
 elseif comstr(Cam,'constants')
  if nargin<3;error('TexStress now handled in integrules');
  elseif varargin{2}(size(varargin{2},1),1)==-9999 % old of_mk_sub.c elements
    i2=24;i1=find(triu(ones(i2,i2)));r1=zeros(i2,i2);r1(i1)=1:length(i1);
    out=struct('ElMap',r1+tril(r1',-1));   out2=[];
    out1=varargin{1};out1(4,:)=3;% Tell of_mk('MatrixIntegration') this is 3d
  else
   [out,i2,out2]=p_solid('const','hexa8',varargin{2:end});
   out1=varargin{1};out1(4,:)=i2; % Tell MatrixIntegration this is 3d 
  end

 elseif comstr(Cam,'node');   out = 1:8;
 elseif  comstr(Cam,'prop');  out = [9 10 11];
 elseif comstr(Cam,'dof')
   out =1:8;k=(1:3)'/100;out=out(ones(3,1),:)+k(:,ones(8,1));out=out(:);
 elseif comstr(Cam,'line')
   out = [1 2 3 4 1 5 6 7 8 5 0 2 6 0 3 7 0 4 8];
 elseif comstr(Cam,'patch')
   out = [1 4 3 2;5 1 2 6;5 6 7 8;4 8 7 3;5 8 4 1;2 3 7 6];
 elseif  comstr(Cam,'edge');   out = [1 2; 2 3; 3 4; 4 1; 1 5; 2 ...
                    6; 3 7; 4 8; 5 6; 6 7; 7 8; 8 5];    
 elseif  comstr(Cam,'face');   out = [1 4 3 2; 1 5 8 4; 1 2 6 5; ...
                    5 6 7 8; 2 3 7 6; 4 8 7 3]; %3 4 8 7
 elseif  comstr(Cam,'flip'); out=[5:8 1:4];out1=1:8; % negative volume flip
 elseif comstr(Cam,'sci_face'); out =  [1 2 3 4;5 8 7 6;3 7 8 4;2 6 7 3;1 5 6 2;1 4 8 5];
 elseif  comstr(Cam,'parent')
   out = 'hexa8';
 elseif  comstr(Cam,'state');out='';% State should be handled by of_mk

 % Basic matrix test of the element - - - - - - - - - - - - - - - - - -
 elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,5); 

   model=femesh('testhexa8');
   if isunix;model.il=p_solid('dbval 111 d3 -3');end% legacy fail
   [m,k,mdof]=fe_mknl(model,'NoT');out=stack_cell(k,m);
   disp('TestMat passed');

 elseif  comstr(Cam,'test'); [CAM,Cam] = comstr(CAM,5);

   [out,out1]=femesh(strcat(['teststruct hexa8' Cam]));
 
 % Basic tests of the element - - - - - - - - - - - - - - - - - - - - - -
 else;sdtw('''%s'' unknown',CAM);  
 end

return
end % of standard calls with one input argument

% ----------------------------------------------------------------------------
% element matrix assembly - - - - - - - - - - - - - - - - - - - - - - - - - -
% Legacy calls to of_mk_subs elements

node=CAM; 
elt=varargin{1}; 
point=double(varargin{2});
integ=varargin{3};
constit=varargin{4};


if size(node,2)~=4
 NNode(node(:,1))=1:size(node,1);
 node=node(NNode(elt(1,1:8)),[5:7 1]);
end

if integ(point(6)+7)~=3
 typ=point(5);
 elmap=varargin{5};if isa(elmap,'int32'); elmap = double(elmap); end

 if (typ==0) % mass and stiffness
  [k1,m1]=of_mk('hexa8',int32(point),integ,constit,node);
  k=reshape(k1(elmap),size(elmap,1),size(elmap,2));
  m=reshape(m1(elmap),size(elmap,1),size(elmap,2));
  if k(1)<0; error('Reorient elements using feutil(''orient'')'); end
  out=k; out1=m;
 elseif typ>99
  warning('OpenFEM:Elem','This standard call should be done in fe_mknl');
  out=of_mk('hexa8',int32(point),varargin{3:4},node,varargin{[5 6 7 8]});
 else
   k1=of_mk('hexa8',int32(point),integ,constit,node);
   k=reshape(k1(elmap),size(elmap,1),size(elmap,2));
   out=k; out1=[];
 end

end
% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------
