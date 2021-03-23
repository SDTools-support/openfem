function [out,out1,out2]=q4p(CAM,varargin);

%Q4P element function for the 4-node 8-DOF quadrangular 2D element (2Q1D)
%
%	In an model description matrix a group of Q4P elements starts with a
%	header row [Inf  abs('q4p') 0 ...] followed by element property rows
%       ELT following the format
%	    [n1 ... n4 MatId ProId EltId]
%         with
%	   n1 ... n4  identification numbers for the element nodes
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

%       Etienne Balmes, Jean-Michel Leclere
%       Copyright (c) 2001-2009 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

% standard calls with one input argument
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1);
 % Build Constit (if isotropic Constit=[rho eta E nu])
 %       Integ   
 %   and Elmap 
 % for later integration
 if comstr(Cam,'integinfo')

  %constit integ,elmap                 ID,pl,il
  [out,out1,out2]=p_solid('buildconstit 2 4',varargin{:});

 elseif comstr(Cam,'matcall'); [out,out1]=elem0(CAM,varargin{:});
 elseif comstr(Cam,'dofcall');out=elem0('dofcall'); % variable field elements
 elseif comstr(Cam,'call')
   out='[k1,m1]=q4p(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,elmap);';

 elseif comstr(Cam,'rhscall') % call for load assembly
   out='rhs_of';

 elseif comstr(Cam,'groupinit');out=elem0(CAM,'q4p');
 % Here one defines the constants needed for element integration
 elseif comstr(Cam,'constants');

  if nargin<3; p_solid('const','q4p',[],[]);return;
  elseif varargin{2}(size(varargin{2},1),1)==-9999; % old of_mk_sub.c elements
    i2=8;i1=find(triu(ones(i2,i2)));r1=zeros(i2,i2);r1(i1)=1:length(i1);
    out=struct('ElMap',r1+tril(r1',-1));   out2=[];
    out1=varargin{1};out1(4,:)=2;% Tell of_mk('MatrixIntegration') this is 2d
  else;
   [out,i2,out2]=p_solid('const','q4p',varargin{2:end});
   out1=varargin{1};out1(4,:)=i2; % Tell MatrixIntegration this is 3d 
  end
 elseif  comstr(Cam,'node');  out = [1 2 3 4];
 elseif  comstr(Cam,'prop');  out = [5 6 7]; 
 elseif  comstr(Cam,'dof');   out=[1.01 1.02 2.01 2.02 3.01 3.02 4.01 4.02]';
 elseif  comstr(Cam,'line');   out = [1 2 3 4 1];
 elseif  comstr(Cam,'patch');  out = [1 2 3 4];
 elseif  comstr(Cam,'edge');   out = [1 2; 2 3; 3 4; 4 1];
 elseif  comstr(Cam,'face');   out = [1 2 3 4];
 elseif  comstr(Cam,'flip');  out = [2 3 4]; out1 = [4 3 2]; 
 elseif  comstr(Cam,'sci_face'); out = [1 2 3 4];
 elseif  comstr(Cam,'parent'); out = 'quad4'; 

 % Basic matrix test of the element - - - - - - - - - - - - - - - - - -
 elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,5);

   model=femesh('testq4p');
   if nargin==2; model.pl=varargin{1}; end
   if ~isempty(strfind(Cam,'_0'));     model.il(:,3)=0;
   elseif ~isempty(strfind(Cam,'_1')); model.il(:,3)=1;
   elseif ~isempty(strfind(Cam,'_2')); model.il(:,3)=2;
   end
   % model.il(:,5)=-3; no legacy
%   [constit,iopt,elmap]=q4p('integinfo',[model.Elt(2,5:6)';8;4],model.pl,model.il);
%   [k,m]=q4p(model.Node,model.Elt(2,:),[36 36 0 0 0 0 0 0 0],int32(iopt),constit,elmap);
%   iopt(5)=4; constit(2)=.5;
%   k1=q4p(model.Node,model.Elt(2,:),[36 36 0 0 4 0 0 0 0],int32(iopt),constit,elmap);
   [m,k,mdof]=fe_mknl(model,'NoT');
   out=stack_cell(k,m);
   if nargout==0;disp('TestMat passed');end

  elseif  comstr(Cam,'test');   [CAM,Cam] = comstr(Cam,5);

   [out,out1]=femesh(strcat('teststruct q4p',Cam),varargin{:});
 elseif comstr(CAM,'cvs')
  out='$Revision: 1.36 $  $Date: 2017/05/09 10:55:14 $'; return;

 % Basic tests of the element - - - - - - - - - - - - - - - - - - - - - -
 else sdtw('''%s'' unknown',CAM);  
 end

  return
end % of standard calls with one input argument

% -----------------------------------------------------------------------------
% element matrix assembly - - - - - - - - - - - - - - - - - - - - - - - - - -
% OBSOLETE of_mk_subs.c elements - - - - - - - - - - - - - - - - - - - - - - 
node=CAM; 
elt=varargin{1}; 
point=varargin{2};
integ=varargin{3};
constit=varargin{4};
elmap=varargin{5};
if isa(elmap,'int32'); elmap = double(elmap); end

if size(node,2)~=4
 NNode(node(:,1))=1:size(node,1);
 node=node(NNode(elt(1,1:4)),[5:7 1]);
end

typ=point(5); 

if (typ==0) % mass and stiffness

  [k1,m1]=of_mk('q4p',int32(point),integ,constit,node);
  k=reshape(k1(elmap),size(elmap,1),size(elmap,2)); 
  if m1(end);m=reshape(m1(elmap),size(elmap,1),size(elmap,2)); 
  else; m=zeros(8,8);m(1:9:64)=m1(1:8); % m is diagonal
  end
  if k(1)<0; error('Reorient elements using feutil(''orient'')'); end
  out=k; out1=m;

elseif typ>99
  warning('OpenFEM:DIAG','This standard call should be done in fe_mknl');
  out=of_mk('q4p',int32(point),varargin{3:4},node,varargin{[5 6 7 8]});
else

   k1=of_mk('q4p',int32(point),integ,constit,node);
   if typ==2&&k1(end)==0;k=zeros(8,8);k(1:9:64)=k1(1:8);% m is diagonal
   else; k=reshape(k1(elmap),size(elmap,1),size(elmap,2));
   end
   out=k; out1=[];

end

return

% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------

