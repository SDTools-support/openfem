function [out,out1,out2]=line2(CAM,varargin)

%lin2 topology holder for two node line elements
%
%
%	See sdtweb     bar1, eltfun, elem0
%	See also help  bar1, p_beam, elem0

%	Etienne Balmes
%       Copyright (c) 2001-2015 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

% standard calls with one input argument
%#ok<*NOSEM,*ASGLU>
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1);
 % Build Constit, Integ, and Elmap for later integration
 if comstr(Cam,'integinfo')
  %constit integ,elmap                 ID,pl,il
  [out,out1,out2]=p_solid('buildconstit 1 2',varargin{:});
 elseif comstr(Cam,'dofcall');out=elem0('dofcall'); % variable field elements
 elseif comstr(Cam,'call')||comstr(Cam,'matcall');  % call for matrix assembly
     [out,out1]=elem0(CAM,varargin{:});
 elseif comstr(Cam,'rhscall') % call for load assembly
  
   out='be=node1(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,elmap,node,defe);';

 elseif  comstr(Cam,'state');   out = ''; % no state evaluation yet
 elseif comstr(Cam,'groupinit');out=elem0(CAM,'line2');
 elseif comstr(Cam,'constants');     
  if nargin<3; p_solid('constsolid','line2',[],[]);return;
  else; %sdtweb p_solid('constsolid')
   [out,i2,out2]=p_solid('constsolid','line2',varargin{2:end});
   out1=varargin{1};out1(4,:)=i2; % Tell MatrixIntegration this is 3d 
  end

     
 elseif  comstr(Cam,'node');   out = [1 2];
 elseif  comstr(Cam,'prop');   out = [3 4 5];
 elseif  comstr(Cam,'dof');    out=[];
 elseif  comstr(Cam,'edge');   out =1:2;
 elseif  comstr(Cam,'face');   out =[];
 elseif  comstr(Cam,'sci_face'); out = [1 2 2];  
 elseif  comstr(Cam,'line');   out =[1 2];
 elseif  comstr(Cam,'patch');  out= [1 2];
 elseif  comstr(Cam,'parent'); out = 'bar1';

 % Basic matrix test of the element - - - - - - - - - - - - - - - - - -
 elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,5); 
 elseif  comstr(Cam,'test');   [CAM,Cam] = comstr(Cam,5);
   out={};
 % Basic tests of the element - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(CAM,'cvs')
  out='$Revision: 1.3 $  $Date: 2024/10/01 10:04:58 $'; return;
 else sdtw('''%s'' unknown',CAM);
 end

return
end % of standard calls with one input argument

if varargin{4}(1)==-1; 
    constit=varargin{4};
    [out,out1]=feval(fe_mat('typep',constit(3)),CAM,varargin{:});
end
