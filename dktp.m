function [out,out1,out2]=dktp(CAM,varargin);

%DKTP element function for the 3-node 9-DOF 
%    triangular 3D element 
%
%	In an model description matrix a group of DKTP elements starts with a
%	header row [Inf  abs('dktp') 0 ...] followed by element property rows
%       ELT following the format
%	    [n1 n2 n3 MatId ProId EltId]
%         with
%	   n1 ... n3  identification numbers for the element nodes
%	   MatId  material property identification number
%	   ProId  element property identification number (only used by upcom)
%	   EltId  optional element identifier
%
%     PL material property rows are either isotropic materials
%           [MatId type E nu rho G eta alpha T0]
%           type=fe_mat('m_elastic',1,1)
%      or 2-D anisotropic materials
%           [MatId type E11 E12 E22 E13 E23 E33 rho eta a1 a2 a3]
%           type = fe_mat('m_elastic',1,4)
%
%     IL Element property rows for 2-D elements follow the format
%      [ProId Type  f   d 0 h]
%     with 
%       Type = fe_mat('p_solid',1,2)
%       f    : not used        	    
%       d    : not used
%       h    : thickness
%
%     Standard tests available with dktp('testeig') (mat,eig,load) http://www.mathworks.fr/store/default.do
%
% 	See sdtweb     eltfun, elem0
%	See also help  t3p, ...

%	E. Balmes, J. Leclere, H. Assime
%       Copyright (c) 2001-2009 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

if comstr(CAM,'cvs')
 out='$Revision: 1.21 $  $Date: 2009/07/29 14:37:19 $'; return;
end

% standard calls with one input argument
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1);
 % Build Constit, Integ, and Elmap for later integration
 if comstr(Cam,'integinfo')

if 1==2
  [constit,iopt] = ...
     fe_mat('plil of_mk 2d',varargin{1},varargin{2},varargin{3}); %solid
end
 
  pl=varargin{2}; pl=pl((pl(:,1)==varargin{1}(1)),:);
  % if length(pl)<7; pl(7)=0.; end
  il=varargin{3}; il=il((pl(:,1)==varargin{1}(1)),:);
  % constit = [pl(5) il(6) pl([3 4]) il(6)]; % constit =[rho h E nu h]
  iopt(3:6)=[9 3 0 0];
  [st,i1,i2]=fe_mat('type',pl(2));
  if strcmp(st,'m_elastic')&&i2==1; % isotrope pl = [matid x E nu rho]
    if length(pl)<5; pl(5)=0.; end;
    constit = [pl(5) il(6) pl([3 4]) il(6)]; 
    iopt(8)=1; % flag isotrope
  else % anisotrope pl =[ matid x E11 E12 E22 E13 E23 E33 rho]
    if length(pl)<9;pl(9)=0.;end;
    constit = [pl(9) il(6) pl([3:8]) il(6)] ;
    iopt(8)=3; % flag anisotrope
  end
  out1=iopt(:); out=constit(:); 

  if isempty(constit); error('Material properties are not defined');end
  % element map
  i1=find(triu(ones(9,9))); r1=zeros(9,9);r1(i1)=1:length(i1);
  out2=r1+tril(r1',-1);

 elseif comstr(Cam,'matcall');  out='mat_of';out1=1; % SymFlag
 elseif comstr(Cam,'call')
   out='[k1,m1]=dktp(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,elmap);';

 elseif comstr(Cam,'rhscall') % call for load assembly
   out='be=dktp(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,gstate,defe,EltConst,InfoAtNode);';

  elseif  comstr(Cam,'groupinit');   out = '';
  elseif  comstr(Cam,'node');   out = [1 2 3];
  elseif  comstr(Cam,'prop');   out = [4 5 6];
  elseif  comstr(Cam,'dof')    
   out =[1:3];k=[3 4 5]'/100;
   out=out(ones(3,1),:)+k(:,ones(3,1));out=out(:);
  elseif  comstr(Cam,'line');   out = [1 2 3 1];
  elseif  comstr(Cam,'patch');  out = [1 2 3];
  elseif  comstr(Cam,'edge');   out = [1 2; 2 3; 3 1];
  elseif  comstr(Cam,'face');   out = [1 2 3];
  elseif comstr(Cam,'sci_face'); out = [1 2 3];
  elseif  comstr(Cam,'parent'); out = 'tria3'; 

 % Basic matrix test of the element - - - - - - - - - - - - - - - - - -
  elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,5);
   model=femesh('testdktp');
   [constit,iopt,elmap]=dktp('integinfo',[100;110],model.pl,model.il);
   [k,m]=dktp(model.Node,model.Elt(2,:),[45 45 0 0 0 0 0 0 0],int32(iopt),constit,elmap);

 iopt(5)=4; constit(2)=.5;
   k1=dktp(model.Node,model.Elt(2,:),[45 45 0 0 0 0 0 0 0],int32(iopt),constit,elmap);

   out=stack_cell(k,m);
   disp('TestMat passed');
  elseif  comstr(Cam,'test');   [CAM,Cam] = comstr(Cam,5);

% Basic tests of the element - - - - - - - - - - - - - - - - - - - - - -

%----------------------------
% materiau isotrope
%----------------------------
femesh('testdktp');femesh('divide2');model=femesh('model0');
model.Elt(1,1:6)=[Inf abs('t3p') 0 0];model.Elt=feutil('orient',model);
model.Elt(1,1:6)=[Inf abs('dktp') 0];
model.Elt(2:end,1:3)=model.Elt(2:end,[1 3 2]);
model = fe_case(model,'FixDof','CL encastrement 1','x==0');
data =  struct('sel','groupall','dir',[0 0  -1e-3]);
model = fe_case(model,'AddToCase 1','Fvol','Charge surf',data);
if ~sp_util('issdt') % openfem
    Load=fe_load(model,'Case1');
    [Case,model.DOF]=fe_mknl('init',model);
    k=fe_mknl('assemble',model,Case,1);
    q=ofact(k,Load.def);
    def = struct('def',Case.T*q,'DOF',model.DOF);
else % sdt
    def=fe_simul('static',model);
end
sig=fe_stress('stress gstate',model,def);

%------------------------------
% materiau anisotrope
%------------------------------
femesh('testdktp');femesh('divide2');model=femesh('model0');
E=210e+9; nu=0.285; E1=E; E2=E; nu12=nu;G12=E/(2*(1+nu));
%E1  = 180e+9; E2  = 10.3e+9; G12 = 7.17e+9; nu12 = 0.28;
nu21=nu12*E2/E1;E11=E1/(1-nu12*nu21); E22=E2/(1-nu12*nu21); E12=nu12*E22; 
E13=0; E23=0;E33=G12;
model.pl=[100 fe_mat('m_elastic','SI',4) E11 E12 E22 E13 E23 E33];
model = fe_case(model,'FixDof','CL encastrement 1','x==0');
data =  struct('sel','groupall','dir',[0 0  -1e-3]);
model = fe_case(model,'AddToCase 1','Fvol','Charge surf',data);
if ~sp_util('issdt') % openfem
    Load=fe_load(model,'Case1');
    [Case,model.DOF]=fe_mknl('init',model);
    k=fe_mknl('assemble',model,Case,1);
    q=ofact(k,Load.def);
    def = struct('def',Case.T*q,'DOF',model.DOF);
else % sdt
    def=fe_simul('static',model);
end
sig=fe_stress('stress gstate',model,def);

%---------------------------------------
% script avec maillage de type nopo
%---------------------------------------
fname=fullfile(fileparts(which('ofutil')),'test','f.nopo');
if exist(fname,'file')
  model = nopo(['read -p 2D ' fname]);
model.Elt(1,:)=[Inf abs('dktp') 0];
model.Elt(2:end,4)=100;model.Elt(2:end,5)=110;
%model.pl = [100 fe_mat('m_elastic','SI',1) 200e+9 0.3];
E=200e+9; nu=0.3; E1=E; E2=E; nu12=nu;G12=E/(2*(1+nu));
nu21=nu12*E2/E1; E11=E1/(1-nu12*nu21); E22=E2/(1-nu12*nu21); E12=nu12*E22; 
E13=0; E23=0;E33=G12;
model.pl=[100 fe_mat('m_elastic','SI',4) E11 E12 E22 E13 E23 E33];
model.il = [110 fe_mat('p_solid','SI',2)  0  0  0  0.01 0 0 ];
FEnode=model.Node;i1=femesh('FindNode x==0 | x==1 | y==0 | y==1')+.03;
model = fe_case(model,'FixDof','CL encastrement 1',i1);
data =  struct('sel','groupall','dir',[0 0  -1e+5]);
model = fe_case(model,'AddToCase 1','Fvol','Charge surf',data);
if ~sp_util('issdt') % openfem
    Load=fe_load(model,'Case1');
    [Case,model.DOF]=fe_mknl('init',model);
    k=fe_mknl('assemble',model,Case,1);
    q=ofact(k,Load.def);
    def = struct('def',Case.T*q,'DOF',model.DOF);
else % sdt
    def=fe_simul('static',model);
end
sig=fe_stress('stress gstate',model,def);
end
out=model;out1=def;


end
return
end % of standard calls with one input argument

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
if (typ==0) 
  [k1,m1]=of_mk('dktp',int32(point),integ,constit,node);
  k=reshape(k1(elmap),size(elmap,1),size(elmap,2)); m=diag(m1(1:9));
  out=k; out1=m;
elseif typ==100
  point(5)=100; point(25)=0;varargin{6}=varargin{6}(3,:);
  out=of_mk('dktp',int32(point),varargin{3:4},node,varargin{[5 6 7 8]});

else

   k1=of_mk('dktp',int32(point),integ,constit,node);
   k=reshape(k1(elmap),size(elmap,1),size(elmap,2));m=[];
   out=k; out1=[];

end
% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------




