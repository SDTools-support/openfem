function [out,out1,out2]=flui4(CAM,varargin);

%flui4	3D linear acoustics. Should really be a tetra4 with 
%       appropriate material properties
%
%       PL material property matrix for fluids see m_elastic subtype 2
%
%	See also m_elastic, p_solid

%       Adrien Bobillot, Etienne Balmes
%       Copyright (c) 2001-2009 by INRIA and SDTools, All rights reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

if comstr(CAM,'cvs')
 out='$Revision: 1.33 $  $Date: 2016/10/18 13:14:39 $'; return;
end

% standard calls with one input argument
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1);
 % Build Constit, Integ, and Elmap for later integration
 if comstr(Cam,'integinfo')

  il=varargin{3};
  if isempty(il)||~any(il(:,1)==varargin{1}(2))
   varargin{3}=[varargin{1}(2) fe_mat('p_solid','SI',1) 0 0 0 0];
  end
  %constit integ                          ID,pl,il
  [out ,   out1] = p_solid('buildconstit',[varargin{1};4;4],varargin{2:end});
  %elmap
  i1=[1:4];out2=reshape(1:4^2,4,4);out2(i1,i1)=out2;

 elseif comstr(Cam,'matcall'); out='mat_og'; out1=0; % mat_og and non symmetric
 elseif comstr(Cam,'call'); out=elem0('callmat_og');
 elseif comstr(Cam,'rhscall'); out=elem0('rhsmat_og'); % call for load
 elseif comstr(Cam,'dof');     out=[1 2 3 4]+.19;

 % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 % Callback to define element constants during the fe_mknl init phase
 elseif comstr(Cam,'groupinit');out=elem0(CAM,'tetra4');

 % Here one defines the constants needed for element integration
 elseif comstr(Cam,'constants');

  if nargin>3; out=p_solid('constfluid','tetra4',varargin{2:3});
  else; p_solid('constfluid','tetra4',[102 111 4 4],[7.4e-6;0;.81]);return;
  end
  out1=varargin{1};out1(4,:)=3; % Tell of_mk('MatrixIntegration') this is 3d 

 % Basic matrix test of the element - - - - - - - - - - - - - - - - - -
 elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,5);

   model=femesh('testflui4');model.Elt=feutil('orient',model);
  [m,k,mdof]=fe_mknl(model);
  Case=fe_mknl('init',model);
  if size(Case.GroupInfo{end}.StrainDefinition{1},1)~=3
   error('Mismatch');
  end
  out=stack_cell(full(k),full(m));
  if nargout==0;disp('TestMat passed');end


 elseif  comstr(Cam,'test'); [CAM,Cam] = comstr(Cam,5);
   % Basic tests of the element - - - - - - - - - - - - - - - - - - - - - -
   [out,out1]=femesh(strcat(['teststruct flui4' Cam]));
   %fecom colordataa;
 else
  if nargout==0; tetra4(CAM,varargin{:});
  elseif nargout==1; out=tetra4(CAM,varargin{:});
  elseif nargout==2; [out,out1]=tetra4(CAM,varargin{:});
  elseif nargout==3; [out,out1,out2]=tetra4(CAM,varargin{:});
  end
 end

return
end % of standard calls with one input argument
