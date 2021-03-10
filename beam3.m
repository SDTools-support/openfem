function [out,out1,out2]=beam3(node,varargin);

%BEAM3 element function for the 3-node 18DOF beam element
%
%	In an model description matrix a group of BEAM3 elements starts with a
%	header row [Inf  abs('beam3') 0 ...] followed by element property rows
%       ELT following the format
%	    [n1 n2 n3 MatId ProId EltId]
%         with
%	   n1,n2,n3  node numbers of the tips of the bar element
%	   MatId  material property identification number
%	   ProId  element property identification number
%	   EltId  element identifier (optional, default 0)
%
%       CURRENTLY FOR PLOT ONLY
%
%	See sdtweb     bar1, eltfun, elem0
%	See also help  beam1, elem0

%       Etienne Balmes
%       Copyright (c) 2001-2011 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

% standard calls with one input argument
if ischar(node)
 [CAM,Cam]=comstr(node,1);ElemP='beam3';
 if comstr(Cam,'integinfo');[out,out1,out2]=elem0('integinfo 1 3',varargin{:});
 elseif      comstr(Cam,'cvs')
 out='$Revision: 1.15 $  $Date: 2011/03/07 17:34:15 $'; return;
 elseif      comstr(Cam,'call')
   out = ['[i1,k1,m1] = bar1(node(NNode(elt(cEGI(jElt),[1 2])),:),elt(cEGI(jElt),:),pl,il,[opt(1) jGroup jElt]);'];
 elseif comstr(Cam,'groupinit');out=elem0(CAM,'beam3');
 elseif comstr(Cam,'constants');
  if nargin<4; p_solid('constsolid','beam3',[],[]);return;
  elseif varargin{3}(1)==-3
   [out,i2,out2]=p_solid('constsolid','beam3',varargin{2:end});
   out1=varargin{1};out1(4,:)=i2; % Tell MatrixIntegration this is 3d 
  else; out2=[]; out=[];out1=varargin{1};out1(4,:)=13;
  end
 elseif  comstr(Cam,'node');    out = [1 2 3];
 elseif  comstr(Cam,'prop');    out = [4 5 6];
 elseif  strcmp(Cam,'dof');     out=[1+[1:6]/100 2+[1:6]/100 3+[1:6]/100];
 elseif  comstr(Cam,'edge');    out = [1 2 3];
 elseif  comstr(Cam,'line');    out = [1 3 2];
 elseif  comstr(Cam,'patch');   out = [1 3;3 2];
 elseif comstr(Cam,'sci_face'); out = [1 2 3];
 elseif  comstr(Cam,'parent');  out = ElemP;
 elseif  comstr(Cam,'test'); disp('beam3 test : just for display')
 % redirect the rest to elem0 - - - - - - - - - - - - - - - - - - - -
 else
   if nargin==1; varg={ElemP}; else; varg=varargin;end
   if nargout<=1;    out=elem0(CAM,varg{:});
   else;      [out,out1]=elem0(CAM,varg{:});
   end
 end
return
end % of standard calls with one input argument

error('Element matrix computations available through of_mk only');



