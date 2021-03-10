function [out,out1,out2]=q9a(CAM,varargin);

%Q9A element function for the 9-node 27-DOF 
%    quadrangular AXI element (aq2c)
%	In an model description matrix a group of Q4P elements starts with a
%	header row [Inf  abs('q4a') 0 ...] followed by element property rows
%       ELT following the format
%	    [n1 ... n9 MatId ProId EltId]
%         with
%	   n1 ... n9  identification numbers for the element nodes
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
%     Standard tests available with q9a('testeig') (mat,eig,load) 
%
% 	See sdtweb     eltfun, elem0
%	See also help  q8p, ...

%	Jean-Michel Leclere, Amine Hassim  
%       Copyright (c) 2001-2009 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

if comstr(CAM,'cvs')
 out='$Revision: 1.22 $  $Date: 2009/05/28 16:42:00 $'; return;
end
% standard calls with one input argument
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1);
 % Build Constit, Integ, and Elmap for later integration
 if comstr(Cam,'integinfo')

  [constit,iopt] = ...
     fe_mat('plil of_mk 2d',varargin{1},varargin{2},varargin{3}); %solid
  iopt([3:4])=[27 9]; if length(iopt)<8; iopt(8)=0;end % [nddl NFourier]
  pl=varargin{2}; pl=pl(pl(:,1)==varargin{1}(1),:);
  constit=[constit(1:2) pl(3:4)];

  % element map
  i1=find(triu(ones(27,27))); r1=zeros(27,27);r1(i1)=1:length(i1);
  out2=r1+tril(r1',-1);
  out1=iopt(:); out=constit(:); 

 elseif comstr(Cam,'matcall')

  out=@of_mk;  k=1; % SymFlag

 elseif comstr(Cam,'call')
   out='[k1,m1]=q9a(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,elmap);';

 elseif comstr(Cam,'rhscall') % call for load assembly
   out='be=q9a(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,elmap,gstate);';

  elseif  comstr(Cam,'groupinit');   out=elem0(CAM,'q9a');
  elseif comstr(Cam,'constants');

  if nargin<3; p_solid('const','q9p',[],[]);return;
  elseif varargin{2}(size(varargin{2},1),1)==-9999; % old of_mk_sub.c elements
    i2=27;i1=find(triu(ones(i2,i2)));r1=zeros(i2,i2);r1(i1)=1:length(i1);
    out=struct('ElMap',r1+tril(r1',-1));   out2=[];
    out1=varargin{1};out1(4,:)=2;% Tell of_mk('MatrixIntegration') this is 2d
  else;
   [out,i2,out2]=p_solid('const','q9p',varargin{2:end});
   out1=varargin{1};out1(4,:)=i2; % Tell MatrixIntegration this is 3d 
  end
 elseif  comstr(Cam,'node');  out = [1:9];
  elseif  comstr(Cam,'prop');  out = [10 11 12]; 
  elseif  comstr(Cam,'dof')   
    out = [1:9];k=[1;2;3]/100;out=out(ones(3,1),:)+k(:,ones(9,1));
    out=out(:);
  elseif  comstr(Cam,'line');   out = [1 5 2 6 3 7 4 8 1];
	elseif  comstr(Cam,'patch');  out = [5 2 6 3 7 4 8 1];
  elseif  comstr(Cam,'edge');   out = [1 2 5; 2 3 6; 3 4 7; 4 1 8];
  elseif  comstr(Cam,'face');   out = [1 2 3 4 5 6 7 8 9];
  elseif  comstr(Cam,'sci_face'); out = [1 5 9 8;5 2 6 9;6 3 7 9;7 4 8 9];
  elseif  comstr(Cam,'parent'); out = 'quad9'; 
 elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,5);

   model=femesh('testq9a');
   [m,k,mdof]=fe_mknl(model,'NoT');
   out=stack_cell(k,m);
   disp('TestMat passed');

  elseif  comstr(Cam,'test');  [CAM,Cam] = comstr(Cam,5);

   pl=[100 1 210e9 .285 7800 210e9/2/(1.285)];
   il=[110 2 2 10];
   if ~exist('CAM','var'); Cam='eig'; end
 
   femesh testquad4;femesh quad42quad9;
   femesh('divide10 10');
   femesh('set groupa1 name q9a matid100 proid110')
   mdof=feutil('getdof',FEel0);%mdof=fe_c(mdof,femesh('findnode x==0'),'dof',2);
   [m,k,mdof]=fe_mk(FEnode,FEel0,pl,il,[],mdof,0);
   model=struct('Node',FEnode,'Elt',FEel0,'pl',pl,'il',il,'DOF',mdof);

   if comstr(Cam,'eig') % matrix/eigenvalue test
    [md1,f1]=fe_eig(m,k,[0.5 8 1e3 11]);
    def=struct('def',md1,'DOF',model.DOF,'data',f1/2/pi);

   elseif comstr(Cam,'vol')||comstr(Cam,'sur')||comstr(Cam,'pre')
    error('RHS for q9a not implemented');
   else error('Not a valid test');
   end
 
  end
  return
end % of standard calls with one input argument



% element matrix assembly - - - - - - - - - - - - - - - - - - - - - - - - - -
node=CAM; 
elt=varargin{1}; 
point=varargin{2};
iopt=varargin{3};
constit=varargin{4};
elmap=varargin{5};
if isa(elmap,'int32'); elmap = double(elmap); end


if size(node,2)~=4
 NNode(node(:,1))=1:size(node,1);
 node=node(NNode(elt(1,1:9)),[5:7 1]);
end

typ=point(5);
if (typ==0) % mass and stiffness

  [k1,m1]=of_mk('q9a',int32(point),iopt,constit,node);
  k=reshape(k1(elmap),size(elmap,1),size(elmap,2)); 
  m=reshape(m1(elmap),size(elmap,1),size(elmap,2));
else
   k1=of_mk('q9a',int32(point),iopt,constit,node);
   k=reshape(k1(elmap),size(elmap,1),size(elmap,2));m=[];
end
if k(1)<0; error('Reorient elements using feutil(''orient'')'); end
out=k;
out1=m;

% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------


