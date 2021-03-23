function [out,out1,out2]=penta15b(CAM,varargin);

%PENTA15B 15 node, 5 sided volume element supporting 
%
%     Should be replaced by PENTA15 elements with appropriate 
%       P_SOLID, P_HEAT or P_PIEZO element property entries
%
%       See sdtweb  penta15, integrules, p_solid

%	Etienne Balmes, Jean-Michel Leclere
%       Copyright (c) 2001-2009 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

if ischar(CAM);

 [CAM,Cam]=comstr(CAM,1);
 if comstr(Cam,'call'); out=elem0('callmat_og');
 elseif comstr(Cam,'rhscall'); out=elem0('rhsmat_og'); % call for load
 elseif  comstr(Cam,'testmat')
  model=femesh('testpenta15bb');model.il(:,4)=-3;
  m1=fe_mk(model);
  out=m1.K([2 1]);  
  if nargout==0;disp('TestMat passed');end
 elseif  comstr(Cam,'test'); [CAM,Cam] = comstr(CAM,5);
   [out,out1]=femesh(strcat('teststruct penta15b',Cam));
 elseif comstr(CAM,'cvs')
  out='$Revision: 1.17 $  $Date: 2009/05/28 16:42:00 $'; return;
 else; 
  if nargout==0; penta15(CAM,varargin{:});
  elseif nargout==1; out=penta15(CAM,varargin{:});
  elseif nargout==2; [out,out1]=penta15(CAM,varargin{:});
  elseif nargout==3; [out,out1,out2]=penta15(CAM,varargin{:});
  end
 end % commands

 return;

 end % of standard calls with one input argument
