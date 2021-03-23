function [out,out1,out2]=node1(CAM,varargin)

%node1 topology holder for one node elements
%
%
%	See sdtweb     bar1, eltfun, elem0
%	See also help  bar1, p_beam, elem0

%	Etienne Balmes
%       Copyright (c) 2001-2012 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

% standard calls with one input argument
%#ok<*ASGLU,*NASGU>
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1);
 % Build Constit, Integ, and Elmap for later integration
 if comstr(Cam,'integinfo')

  %constit integ,elmap                 ID,pl,il
  [out,out1,out2]=p_solid('buildconstit 0 1',varargin{:});
 elseif comstr(Cam,'dofcall');out=elem0('dofcall'); % variable field elements
 elseif comstr(Cam,'call')||comstr(Cam,'matcall');  % call for matrix assembly
   out='[k1,m1]=node1(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,elmap);';
   out1=[];
 elseif comstr(Cam,'rhscall') % call for load assembly
  
   out='be=node1(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,elmap,node,defe);';

 elseif  comstr(Cam,'state');   out = ''; % no state evaluation yet
 elseif comstr(Cam,'groupinit');out=elem0(CAM,'beam1');
 elseif comstr(Cam,'constants');out1=[]; out2=[]; out=[];
 elseif  comstr(Cam,'node');   out = 1;
 elseif  comstr(Cam,'prop');   out = [2 3 4];
 elseif  comstr(Cam,'dof');    out=1+(1:6)/100;
 elseif  comstr(Cam,'edge');   out =1;
 elseif  comstr(Cam,'face');   out =[];
 elseif  comstr(Cam,'sci_face'); out = [1 1 1];  
 elseif  comstr(Cam,'line');   out =1;
 elseif  comstr(Cam,'patch');  out= 1;
 elseif  comstr(Cam,'parent'); out = 'mass1';

 % Basic matrix test of the element - - - - - - - - - - - - - - - - - -
 elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,5);  
 elseif  comstr(Cam,'test');   [CAM,Cam] = comstr(Cam,5);
 % Basic tests of the element - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(CAM,'cvs')
  out='$Revision: 1.4 $  $Date: 2012/10/25 08:55:14 $'; return;
 else sdtw('''%s'' unknown',CAM);
 end

return
end % of standard calls with one input argument

if varargin{4}(1)==-1; 
    constit=varargin{4};
    [out,out1]=feval(fe_mat('typep',constit(3)),CAM,varargin{:});
end
