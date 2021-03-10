function [out,out1,out2]=q5p(CAM,varargin);

%Q5P element function for the 5-node 10-DOF 
%    quadrangular incompressible 2D element (5noe)
%
%	In an model description matrix a group of Q5P elements starts with a
%	header row [Inf  abs('q5p') 0 ...] followed by element property rows
%       ELT following the format
%	    [n1 ... n5 MatId ProId EltId]
%         with
%	   n1 ... n5  identification numbers for the element nodes
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
%     Standard tests available with q5p('testeig') (mat,eig,load) 
%
% 	See sdtweb     eltfun, elem0
%	See also help  t3p, ...


%	Jean-Michel Leclere, Amine Hassim  
%       Copyright (c) 2001-2009 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

if comstr(CAM,'cvs')
 out='$Revision: 1.21 $  $Date: 2009/05/28 16:42:00 $'; return;
end
% standard calls with one input argument
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1);
 % Build Constit, Integ, and Elmap for later integration
 if comstr(Cam,'integinfo')

  [constit,iopt] = ...
     fe_mat('plil of_mk 2d',varargin{1},varargin{2},varargin{3}); %solid
  iopt([3:4 8])=[10 5 varargin{3}(3)+1]; % [nddl nnode flag_stres/strain]
  tmp_varar = varargin{2}(3:4);
  if size(tmp_varar,1)~=1; tmp_varar = tmp_varar'; end
  constit=[constit(1:2) tmp_varar];

  % element map
  i1=find(triu(ones(10,10))); r1=zeros(10,10);r1(i1)=1:length(i1);
  out2=r1+tril(r1',-1);
  out1=iopt(:); out=constit(:); 

 elseif comstr(Cam,'matcall')

  out=@of_mk;  k=1; % SymFlag

 elseif comstr(Cam,'call')
   out='[k1,m1]=q5p(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,elmap);';

 elseif comstr(Cam,'rhscall') % call for load assembly
   out='rhs_of';

  elseif comstr(Cam,'groupinit');   out = '';
  elseif  comstr(Cam,'node')  
   out = [1 2 3 4 5];
  elseif  comstr(Cam,'prop')  
   out = [6 7 8]; 
  elseif  comstr(Cam,'dof')   
   out=[1.01 1.02 2.01 2.02 3.01 3.02 4.01 4.02 5.01 5.02]';
  elseif  comstr(Cam,'line')   
   out = [1 2 3 4 1];
  elseif  comstr(Cam,'patch')  
   out = [1 2 3 4 1];
  elseif  comstr(Cam,'edge')   
   out = [1 2; 2 3; 3 4; 4 1];
  elseif  comstr(Cam,'face')   
   out = [1 2 3 4];
  elseif  comstr(Cam,'flip');   out=[2 3 4]; out1=[4 3 2]; 
  elseif  comstr(Cam,'sci_face')
   out = [1 2 5 5;2 3 5 5;3 4 5 5;4 1 5 5];
  elseif  comstr(Cam,'parent')   
   out = 'q5p'; 
 elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,5);

   model=femesh('testq5p');
   if nargin==2; model.pl=varargin{1}; end
   if ~isempty(strfind(Cam,'_0'));     model.il(:,3)=0;
   elseif ~isempty(strfind(Cam,'_1')); model.il(:,3)=1;
   elseif ~isempty(strfind(Cam,'_2')); model.il(:,3)=2;
   end
   [m,k,mdof]=fe_mknl(model,'NoT');out=stack_cell(k,m);
   disp('TestMat passed');

  elseif  comstr(Cam,'test');   [CAM,Cam] = comstr(Cam,5);

   if nargin==2 % specified pl
     [out,out1]=femesh(strcat(['teststruct q5p' Cam],varargin{1}));
   else
     [out,out1]=femesh(strcat(['teststruct q5p' Cam]));
   end

 % Basic tests of the element - - - - - - - - - - - - - - - - - - - - - -
 else sdtw('''%s'' unknown',CAM);  end

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
 node=node(NNode(elt(1,1:5)),[5:7 1]);
end

typ=point(5);
if (typ==0) % mass and stiffness
  [k1,m1]=of_mk('q5p',int32(point),integ,constit,node);
  k=reshape(k1(elmap),size(elmap,1),size(elmap,2)); 
  m=reshape(m1(elmap),size(elmap,1),size(elmap,2));
  if k(1)<0; error('Reorient elements using feutil(''orient'')'); end
  out=k; out1=m;
elseif typ>99
  warning('This standard call should be done in fe_mknl');
  out=of_mk('q5p',int32(point),varargin{3:4},node,varargin{[5 6 7 8]});
else
   k1=of_mk('q5p',int32(point),integ,constit,node);
   k=reshape(k1(elmap),size(elmap,1),size(elmap,2));
   out=k; out1=[];
end

return

% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------

