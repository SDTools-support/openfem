function [out,out1,out2]=flui6(CAM,varargin);

%flui4	3D linear acoustics. Should really be a penta6 with 
%       appropriate material properties
%
%       PL material property matrix for fluids see m_elastic subtype 2
%
%	See also m_elastic, p_solid

%       Adrien Bobillot, Etienne Balmes
%       Copyright (c) 2001-2009 by INRIA and SDTools, All rights reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

if comstr(CAM,'cvs')
 out='$Revision: 1.26 $  $Date: 2009/05/25 17:19:22 $'; return;
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
  [out ,   out1] = p_solid('buildconstit',[varargin{1};6;6],varargin{2:end});
  %elmap
  i1=[1:6];out2=reshape(1:6^2,6,6);out2(i1,i1)=out2;

 elseif comstr(Cam,'matcall'); out='mat_og'; out1=0; % mat_og and non symmetric
 elseif comstr(Cam,'call'); out=elem0('callmat_og');
 elseif comstr(Cam,'rhscall'); out=elem0('rhsmat_og'); % call for load assembly

 % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 % Callback to define element constants during the fe_mknl init phase
 elseif comstr(Cam,'groupinit');out=elem0(CAM,'penta6');

 % Here one defines the constants needed for element integration
 elseif comstr(Cam,'constants');

  constit=varargin{3}; 
  out=p_solid('constfluid','penta6',varargin{2},varargin{3});
  out1=varargin{1};out1(4,:)=3; % Tell of_mk('MatrixIntegration') this is 3d 

 elseif comstr(Cam,'dof');     out =[1:6]+.19;
 % Basic matrix test of the element - - - - - - - - - - - - - - - - - -
 elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,5);

   model=femesh('testflui6');
   [m,k,mdof]=fe_mknl(model);
   out=stack_cell(full(k),full(m));
   disp('TestMat passed');

 elseif  comstr(Cam,'test'); [CAM,Cam] = comstr(Cam,5);
   % Basic tests of the element - - - - - - - - - - - - - - - - - - - - - -
   [out,out1]=femesh(strcat(['teststruct flui6' Cam]));
   %fecom colordataa;

 else
  if nargout==0; penta6(CAM,varargin{:});
  elseif nargout==1; out=penta6(CAM,varargin{:});
  elseif nargout==2; [out,out1]=penta6(CAM,varargin{:});
  elseif nargout==3; [out,out1,out2]=penta6(CAM,varargin{:});
  end
 end

 return
end % of standard calls with one input argument

