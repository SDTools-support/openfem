function [out,out1,out2]=t3p(CAM,varargin);

%T3P element function for the 3-node 6-DOF 
%    triangular 2D element (2P1D_DP and 2P1D_CP)
%
%	In an model description matrix a group of T3P elements starts with a
%	header row [Inf  abs('t3p') 0 ...] followed by element property rows
%       ELT following the format
%	    [n1 ... n3 MatId ProId EltId]
%         with
%	   n1 ... n3  identification numbers for the element nodes
%	   MatId  material property identification number
%	   ProId  element property identification number (only used by upcom)
%	   EltId  optional element identifier
%
%     PL material property rows are either isotropic materials (Subtype 1)
%           [MatId type E nu rho G eta alpha T0]
%           type=fe_mat('m_elastic','SI',1)
%      or 2-D anisotropic materials (Subtype 4)
%           [MatId type E11 E12 E22 E13 E23 E33 rho eta a1 a2 a3]
%           type = fe_mat('m_elastic','SI',4)
%     See m_elastic for details on PL 
%
%     IL Element property rows for 2-D elements follow the format
%      [ProId Type Form N]
%     with 
%       Type = fe_mat('p_solid','SI',2)
%       Form : formulation (0 plane strain, 1 plane stress, 2 axisymetric)
%       N    : Fourier harmonic for axisymetric elements that support it
%     See p_solid for details on IL 
%
%     Standard tests available with q4p('testeig') (mat,eig,load) 
%
% 	See sdtweb     eltfun, elem0
%	See also help  t3p, ...

%       Jean-Michel Leclere, Amine Hassim  
%       Copyright (c) 2001-2008 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

if comstr(CAM,'cvs')
 out='$Revision: 1.31 $  $Date: 2009/05/28 16:42:00 $'; return;
end
% standard calls with one input argument
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1);
 % Build Constit, Integ, and Elmap for later integration
 if comstr(Cam,'integinfo')

  %constit integ,elmap                 ID,pl,il
  [out,out1,out2]=p_solid('buildconstit 2 3',varargin{:});

 elseif comstr(Cam,'matcall'); [out,out1]=elem0(CAM,varargin{:});
 elseif comstr(Cam,'dofcall');out=elem0('dofcall'); % variable field elements
 elseif comstr(Cam,'call')
   out='[k1,m1]=t3p(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,elmap);';

 elseif comstr(Cam,'rhscall') % call for load assembly
   out='rhs_of';
 elseif comstr(Cam,'groupinit');out=elem0(CAM,'t3p');
 % Here one defines the constants needed for element integration
 elseif comstr(Cam,'constants');

  if nargin<3; p_solid('constsolid','t3p',[],[]);return;
  elseif varargin{2}(size(varargin{2},1),1)==-9999; % old of_mk_sub.c elements
    i2=6;i1=find(triu(ones(i2,i2)));r1=zeros(i2,i2);r1(i1)=1:length(i1);
    out=struct('ElMap',r1+tril(r1',-1));   out2=[];
    out1=varargin{1};out1(4,:)=2;% Tell of_mk('MatrixIntegration') this is 2d
  else;
   [out,i2,out2]=p_solid('constsolid','t3p',varargin{2:end});
   out1=varargin{1};out1(4,:)=i2; % Tell MatrixIntegration this is 3d 
  end

 elseif  comstr(Cam,'node');   out = [1 2 3];
 elseif  comstr(Cam,'prop');   out = [4 5 6];
 elseif  comstr(Cam,'dof');    out=  [1.01 1.02 2.01 2.02 3.01 3.02]';
 elseif  comstr(Cam,'line');   out = [1 2 3 1];
 elseif  comstr(Cam,'patch');  out = [1 2 3];
 elseif  comstr(Cam,'edge');   out = [1 2; 2 3; 3 1];
 elseif  comstr(Cam,'face');   out = [1 2 3];
 elseif  comstr(Cam,'flip');   out=[2 3];out1=[3 2]; 
 elseif  comstr(Cam,'sci_face'); out = [1 2 3];
 elseif  comstr(Cam,'parent'); out = 'tria3'; 
 elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,5);

   model=femesh('testt3p');
   if nargin==2; model.pl=varargin{1}; end
   if ~isempty(strfind(Cam,'_0'));     model.il(:,3)=0;
   elseif ~isempty(strfind(Cam,'_1')); model.il(:,3)=1;
   elseif ~isempty(strfind(Cam,'_2')); model.il(:,3)=2;
   end
   [m,k,mdof]=fe_mknl(model,'NoT');
   out=stack_cell(k,m);
   disp('TestMat passed');

 elseif  comstr(Cam,'test');   [CAM,Cam] = comstr(Cam,5);

   if nargin==2 % specified pl
     [out,out1]=femesh(strcat(['teststruct t3p' Cam]),varargin{1});
   else
     [out,out1]=femesh(strcat(['teststruct t3p' Cam]));
   end

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
 node=node(NNode(elt(1,1:3)),[5:7 1]);
end

typ=point(5); 

if (typ==0) % mass and stiffness
  [k1,m1]=of_mk('t3p',int32(point),integ,constit,node);
  k=reshape(k1(elmap),size(elmap,1),size(elmap,2));
  m=reshape(m1(elmap),size(elmap,1),size(elmap,2));
  if k(1)<0; error('Reorient elements using feutil(''orient'')'); end
  out=k; out1=m;
elseif typ>99
  warning('OpenFEM:DIAG','This standard call should be done in fe_mknl');
  out=of_mk('t3p',int32(point),varargin{3:4},node,varargin{[5 6 7 8]});
else
   k1=of_mk('t3p',int32(point),integ,constit,node);
   k=reshape(k1(elmap),size(elmap,1),size(elmap,2));
   out=k; out1=[];
end

return

% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------




