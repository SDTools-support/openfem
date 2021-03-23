function [out,out1,out2]=tetra10(CAM,varargin);

%TETRA10B 10 node, 4 sided volume element supporting 
%
%     Should be replaced by TETRA10 elements with appropriate 
%       P_SOLID, P_HEAT or P_PIEZO element property entries
%
%       See sdtweb     tetra10, integrules, p_solid

%	Etienne Balmes, Jean-Michel Leclere
%       Copyright (c) 2001-2009 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

if comstr(CAM,'cvs')
 out='$Revision: 1.17 $  $Date: 2009/05/28 16:42:00 $'; return;
end
if ischar(CAM);

 [CAM,Cam]=comstr(CAM,1);
 if comstr(Cam,'call'); out=elem0('callmat_og');
 elseif comstr(Cam,'rhscall'); out=elem0('rhsmat_og'); % call for load
 elseif  comstr(Cam,'testmat')
  model=femesh('testtetra10b');model.il(:,4)=-3;
  m1=fe_mk(model);
  out=m1.K([2 1]);  
  if nargout==0;disp('TestMat passed');end
 elseif  comstr(Cam,'test'); [CAM,Cam] = comstr(CAM,5);
   [out,out1]=femesh(strcat('teststruct tetra10b',Cam));
 else; 
  if nargout==0; tetra10(CAM,varargin{:});
  elseif nargout==1; out=tetra10(CAM,varargin{:});
  elseif nargout==2; [out,out1]=tetra10(CAM,varargin{:});
  elseif nargout==3; [out,out1,out2]=tetra10(CAM,varargin{:});
  end
 end % commands

 return;

 end % of standard calls with one input argument

