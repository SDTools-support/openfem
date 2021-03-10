function [out,out1,out2]=t6pb(CAM,varargin);

%T6PB element function for the 6-node 12-DOF triangle
%
%     Should be replaced by T6P elements with appropriate 
%       P_SOLID, P_HEAT or P_PIEZO element property entries
%
%       See sdtweb   t6p, integrules, p_solid


%	Jean-Michel Leclere, Amine Hassim  
%       Copyright (c) 2001-2009 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

if comstr(CAM,'cvs')
 out='$Revision: 1.16 $  $Date: 2009/05/28 16:42:00 $'; return;
end
if ischar(CAM);

 [CAM,Cam]=comstr(CAM,1);
 if comstr(Cam,'call'); out=elem0('callmat_og');
 elseif comstr(Cam,'rhscall'); out=elem0('rhsmat_og'); % call for load
 elseif  comstr(Cam,'testmat')
  model=femesh('testt6pb');model.il(:,5)=-3;
  if ~isempty(strfind(Cam,'_0'));     model.il(:,3)=0;
  elseif ~isempty(strfind(Cam,'_1')); model.il(:,3)=1;
  end
  m1=fe_mk(model);
  out=m1.K([2 1]);  
  if nargout==0;disp('TestMat passed');end
 elseif  comstr(Cam,'test'); [CAM,Cam] = comstr(CAM,5);
   [out,out1]=femesh(strcat('teststruct t6pb',Cam));
 else; 
  if nargout==0; t6p(CAM,varargin{:});
  elseif nargout==1; out=t6p(CAM,varargin{:});
  elseif nargout==2; [out,out1]=t6p(CAM,varargin{:});
  elseif nargout==3; [out,out1,out2]=t6p(CAM,varargin{:});
  end
 end % commands

 return;

 end % of standard calls with one input argument
