function [out,out1,out2]=pyra5(CAM,varargin)%(node,elt,pl,il,opt,NNode)

%PYRA5	5-node 15-DOF solid element
%
% 	See sdtweb     eltfun, elem0
%	See also help  hexa20, hexa8, penta6


%       Etienne Balmes
%       Copyright (c) 2001-2018 by INRIA and SDTools
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       All Rights Reserved.

% standard calls with one input argument
if ischar(CAM)
 [CAM,Cam]=comstr(CAM,1);
 if comstr(Cam,'integinfo')
  %constit integ,elmap                 ID,pl,il
  [out,out1,out2]=p_solid('buildconstit 3 5',varargin{:});
 elseif comstr(Cam,'matcall'); [out,out1]=elem0(CAM,varargin{:});
 elseif comstr(Cam,'dofcall');out=elem0('dofcall'); % variable field elements
 elseif comstr(Cam,'call'); out=elem0('callmat_og');
 elseif comstr(Cam,'rhscall'); out=elem0('rhsmat_og'); % call for load
 elseif comstr(Cam,'groupinit');out=elem0(CAM,'pyra5');
 elseif comstr(Cam,'constants')
  if nargin<3;error('TexStress now handled in integrules');
  elseif varargin{2}(size(varargin{2},1),1)==-9999 % old of_mk_sub.c elements
    i2=24;i1=find(triu(ones(i2,i2)));r1=zeros(i2,i2);r1(i1)=1:length(i1);
    out=struct('ElMap',r1+tril(r1',-1));   out2=[];
    out1=varargin{1};out1(4,:)=3;% Tell of_mk('MatrixIntegration') this is 3d
  else
   [out,i2,out2]=p_solid('const','pyra5',varargin{2:end});
   out1=varargin{1};out1(4,:)=i2; % Tell MatrixIntegration this is 3d 
  end
 elseif comstr(Cam,'node');   out = 1:5;
 elseif comstr(Cam,'dof')
   out=[1+[1:3]'/100; 2+[1:3]'/100; 3+[1:3]'/100; 4+[1:3]'/100;5+[1:3]'/100];
 elseif comstr(Cam,'line');    out = [1 2 3 4 1 5 4 0 3 5 2];
 elseif  comstr(Cam,'edge');   out = [1 2; 2 3; 3 4; 4 1; 1 5; 2 5;3 5;4 5];
 elseif  comstr(Cam,'face');   out = [1 4 3 2;4 1 5 5; 1 2 5 5; ...
                                      2 3 5 5; 3 4 5 5];
 elseif  comstr(Cam,'patch');   out = [1 4 3 2;4 1 5 5; 1 2 5 5; ...
                                      2 3 5 5; 3 4 5 5];
 elseif comstr(Cam,'sci_face');  out = [1 4 3 2;4 1 5 5; 1 2 5 5; ...
                                      2 3 5 5; 3 4 5 5];                             
 elseif  comstr(Cam,'parent'); out = 'pyra5';
 elseif  comstr(Cam,'prop');  out = [6 7 8];
 elseif  comstr(Cam,'flip'); out=[1 4 3 2 5];out1=1:5; % negative volume flip

 elseif  comstr(Cam,'testmat')
   model=femesh('testhexa8b');
   model.Elt=feutil('addelt','pyra5',[1:5 100 111]);
   [m,k,mdof]=fe_mknl(model,'NoT');out=stack_cell(k,m); %#ok<*ASGLU>
   
 elseif  comstr(Cam,'test'); [CAM,Cam] = comstr(CAM,5);

   [out,out1]=femesh(strcat(['teststruct pyra5' Cam]));

 elseif comstr(Cam,'cvs')
  out='$Revision: 1.14 $  $Date: 2018/11/12 12:02:51 $'; return;
 else;error('Not an implemented method %s',CAM);
 end
return
end % of standard calls with one input argument

