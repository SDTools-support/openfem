function [out,out1,out2]=tria6(CAM,varargin);

%TRIA6	element function of a 6-node plate element
%
%	As all element functions (see ELEM0), TRIA6 is called by FE_MK for
%	model assembly, FEPLOT for structural deformation visualization, ...
%
%	TRIA6 header rows follow the format [Inf  abs('tria6') 0 EGID ...]
%	TRIA6 element property rows follow the format
%	    [n1 ... n6 MatId ProId (EltId Theta Zoff T1 ... T6)]
%         with
%          n1 ... n6  identification numbers for the element nodes
%	   MatId  material property identification number
%	   ProId  element property identification number
%	   EltId  optional element identification number
%	   Theta  unused
%          Zoff   optional offset from the surface of the nodes to the 
%		  reference plane
%	   T1 ... T6  optional thickness at nodes (used instead of IL entry)
%		  Currently the mean of Ti is used.


% Etienne Balmes
% Copyright (c) 2001-2012 by INRIA and SDTools, All Rights Reserved
% Use under OpenFEM trademark.html license and LGPL.txt library license

% standard calls with one input argument
%#ok<*ASGLU,*NASGU>
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1); 
 % Build Constit, Integ, and Elmap for later integration
 if comstr(Cam,'integinfo')

    %constit integ,elmap                 ID,pl,il
  [out,out1,out2]= ...
   p_shell('buildconstit',[varargin{1};36;6],varargin{2:end});

 elseif comstr(Cam,'matcall')||comstr(Cam,'call');  %Matrix assembly in fe_mknl
  [out,out1]=elem0('calltria6',varargin{:}); 
 elseif comstr(Cam,'rhscall');out=quad4('rhscall'); % call for load assembly
   
 elseif comstr(Cam,'groupinit');out=elem0('groupinitogShell','tria6');
 elseif comstr(Cam,'constants'); error('was moved to elem0');
 elseif  comstr(Cam,'node');    out = [1:6];
 elseif comstr(Cam,'dofcall');out=elem0('dofcall'); % variable field elements
 elseif  comstr(Cam,'dof');
   out = [1+[1:6]'/100;2+[1:6]'/100;3+[1:6]'/100; 
          4+[1:6]'/100;5+[1:6]'/100;6+[1:6]'/100];
 elseif  comstr(Cam,'line'); out = [1 4 2 5 3 6 1]; return; 
 elseif  comstr(Cam,'patch');   out = [1 4 2 5 3 6];return; 
 elseif  comstr(Cam,'prop');   out = [7 8 9];
 elseif  comstr(Cam,'parent');  out = 'tria6';
 elseif  comstr(Cam,'edge');   out = [1 2 4; 2 3 5; 3 1 6];
 elseif  comstr(Cam,'face');   out = [1 2 3 4 5 6];
 elseif  comstr(Cam,'flip');   out=[1 3 2 6 5 4];out1=1:6; 
 elseif  comstr(Cam,'sci_face');  out = [1 4 6;4 2 5;5 3 6;4 5 6];
 elseif  comstr(Cam,'cvs')
  out='$Revision: 1.22 $  $Date: 2012/03/13 18:15:42 $';
 elseif  comstr(Cam,'testmat'); 
   model=femesh('testt6p');
   model.Elt(1,1:7)=[Inf abs('tria6') 0];model.Elt(2,8)=110;
   model.il(model.il(:,1)==110,3)=5;
   [m,k,mdof]=fe_mknl(model); 
   out=stack_cell(full(k),full(m));
   disp('TestMat passed');
 elseif  comstr(Cam,'test'); [CAM,Cam] = comstr(Cam,5);
   [out,out1]=femesh(strcat(['teststruct tria6' Cam]));
 end
 return
end % of standard calls with one input argument

warning('OpenFEM:FORMULATION', ...
    'No TRIA6 matrix in OpenFEM yet, but implemented in SDT/Q4CS');
out='error';
