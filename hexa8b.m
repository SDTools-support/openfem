function [out,out1,out2,out3]=hexa8b(CAM,varargin);

%HEXA8B 8 node, 6 sided volume element supporting 
%
%     Should be replaced by HEXA8 elements with appropriate 
%       P_SOLID, P_HEAT or P_PIEZO element property entries
%
%       See sdtweb   hexa8, integrules, p_solid

%	Etienne Balmes, Jean-Michel Leclere
%       Copyright (c) 2001-2009 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

if comstr(CAM,'cvs')
 out='$Revision: 1.32 $  $Date: 2009/05/28 16:42:00 $'; return;
end
if ischar(CAM);

 [CAM,Cam]=comstr(CAM,1);
 if comstr(Cam,'call'); out=elem0('callmat_og');
 elseif comstr(Cam,'rhscall'); out=elem0('rhsmat_og'); % call for load
 elseif  comstr(Cam,'testmat')
  model=femesh('testhexa8b'); model.il(:,4)=-3;
  m1=fe_mknl(model,'NoT');out=m1.K([2 1]);  
  if nargout==0;disp('TestMat passed');end
 elseif  comstr(Cam,'test'); [CAM,Cam] = comstr(CAM,5);
   [out,out1]=femesh(strcat('teststruct hexa8b',Cam));
 else; 
  if nargout==0; hexa8(CAM,varargin{:});
  elseif nargout==1; out=hexa8(CAM,varargin{:});
  elseif nargout==2; [out,out1]=hexa8(CAM,varargin{:});
  elseif nargout==3; [out,out1,out2]=hexa8(CAM,varargin{:});
  end
 end % commands

 return;

 end % of standard calls with one input argument
