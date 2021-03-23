function [out,out1,out2]=line3(CAM,varargin);

%line3 topology holder for 3 node line elements
%
%
%	See sdtweb     bar1, eltfun, elem0
%	See also help  bar1, p_beam, elem0

%	Etienne Balmes
%       Copyright (c) 2001-2012 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

% standard calls with one input argument
%#ok<*ASGLU,*NOSEM,*NASGU>

if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1);
 % Build Constit, Integ, and Elmap for later integration
 if comstr(Cam,'integinfo')
  %constit integ,elmap                 ID,pl,il
  [out,out1,out2]=p_solid('buildconstit 0 3',varargin{:});
 elseif comstr(Cam,'dofcall');out=elem0('dofcall'); % variable field elements
 elseif comstr(Cam,'call')||comstr(Cam,'matcall');  % call for matrix assembly
     [out,out1]=elem0(CAM,varargin{:});
 elseif comstr(Cam,'rhscall') % call for load assembly
   out='be=node1(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,elmap,node,defe);';
 elseif  comstr(Cam,'state');   out = ''; % no state evaluation yet
 elseif comstr(Cam,'groupinit');out=elem0(CAM,'line3');
 elseif comstr(Cam,'constants');     
  if nargin<3; p_solid('constsolid','line3',[],[]);return;
  else; %sdtweb p_solid('constsolid')
   [out,i2,out2]=p_solid('constsolid','line3',varargin{2:end});
   out1=varargin{1};out1(4,:)=i2; % Tell MatrixIntegration this is 3d 
  end

     
 elseif  comstr(Cam,'node');   out = [1 2 3];
 elseif  comstr(Cam,'prop');   out = [4 5 6];
 elseif  comstr(Cam,'dof');    out=[];
 elseif  comstr(Cam,'edge');   out =1:3;
 elseif  comstr(Cam,'face');   out =[];
 elseif  comstr(Cam,'sci_face'); out = [1 2 2];  
 elseif  comstr(Cam,'line');   out =[1 3 2];
 elseif  comstr(Cam,'patch');  out= beam3('patch');
 elseif  comstr(Cam,'parent'); out = 'beam3';

 % Basic matrix test of the element - - - - - - - - - - - - - - - - - -
 elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,5);  
 elseif  comstr(Cam,'test');   [CAM,Cam] = comstr(Cam,5);
 % Basic tests of the element - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(CAM,'cvs')
  out='$Revision: 1.2 $  $Date: 2012/10/25 17:10:27 $'; return;
 else sdtw('''%s'' unknown',CAM);
 end

return
end % of standard calls with one input argument

if varargin{4}(1)==-1; 
    constit=varargin{4};
    [out,out1]=feval(fe_mat('typep',constit(3)),CAM,varargin{:});
end
