function [out,out1,out2]=flui8(CAM,varargin);

%flui4	3D linear acoustics. Should really be a hexa8 with 
%       appropriate material properties
%
%       PL material property matrix for fluids see m_elastic subtype 2
%
%	See also m_elastic, p_solid

%       Adrien Bobillot, Etienne Balmes
%       Copyright (c) 2001-2009 by INRIA and SDTools, All rights reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

if comstr(CAM,'cvs')
 out='$Revision: 1.24 $  $Date: 2009/05/25 17:19:22 $'; return;
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
  [out ,out1]=p_solid('buildconstit',[varargin{1};8;8],varargin{2:end});
  %elmap
  i1=[1:8];out2=reshape(1:8^2,8,8);out2(i1,i1)=out2;

 elseif comstr(Cam,'matcall'); out='mat_og'; out1=0; % mat_og and non symmetric
 elseif comstr(Cam,'call'); out=elem0('callmat_og');
 elseif comstr(Cam,'rhscall'); out=elem0('rhsmat_og'); % call for load
 % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 % Callback to define element constants during the fe_mknl init phase
 elseif comstr(Cam,'groupinit');out=elem0(CAM,'hexa8');

 % Here one defines the constants needed for element integration
 elseif comstr(Cam,'constants');

  constit=varargin{3}; 
  out=p_solid('constfluid','hexa8',varargin{2},varargin{3});
  out1=varargin{1};out1(4,:)=3; % Tell of_mk('MatrixIntegration') this is 3d 

 % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 elseif comstr(Cam,'groupinit');   out = '';
 elseif comstr(Cam,'dof');     out =[1:8]'+.19;
 elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,5);

   model=femesh('testflui8 divide 1 1 1');
   [m,k,mdof]=fe_mknl(model);
   out=stack_cell(full(k),full(m));
   disp('TestMat passed');
 elseif  comstr(Cam,'test'); [CAM,Cam] = comstr(Cam,5);
   % Basic tests of the element - - - - - - - - - - - - - - - - - - - - - -
   [out,out1]=femesh(strcat(['teststruct flui8' Cam]));
   %fecom colordataa;

 else
  if nargout==0; hexa8(CAM,varargin{:});
  elseif nargout==1; out=hexa8(CAM,varargin{:});
  elseif nargout==2; [out,out1]=hexa8(CAM,varargin{:});
  elseif nargout==3; [out,out1,out2]=hexa8(CAM,varargin{:});
  end
 end

return
end % of standard calls with one input argument

