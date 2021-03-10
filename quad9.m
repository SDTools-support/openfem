function [out,out1,out2]=quad9(CAM,varargin);
%QUAD9 support for a 9 node surface topology (Q2)
%


%       Etienne Balmes
%       Copyright (c) 2001-2009 by INRIA and SDTools
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       All Rights Reserved.

if comstr(CAM,'cvs')
 out='$Revision: 1.21 $  $Date: 2020/06/19 12:41:39 $'; return;
end
% standard calls with one input argument
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1);
 % Build Constit, Integ, and Elmap for later integration
 if comstr(Cam,'integinfo')

  %constit integ,elmap                 ID,pl,il
  [out,out1,out2]=p_solid('buildconstit 2 9',varargin{1:end});

 elseif comstr(Cam,'matcall'); [out,out1]=elem0(CAM,varargin{:});
 elseif comstr(Cam,'dofcall');out=elem0('dofcall'); % variable field elements
 elseif comstr(Cam,'call'); out=elem0('callmat_og');
 elseif comstr(Cam,'rhscall'); out=elem0('rhsmat_og'); % call for load
 elseif comstr(Cam,'groupinit');out=elem0(CAM,'quad9');
 % Here one defines the constants needed for element integration
 elseif comstr(Cam,'constants');

  if nargin<3; p_solid('const','quad9',[],[]);return;
  elseif varargin{2}(size(varargin{2},1),1)==-9999; % old of_mk_sub.c elements
    i2=18;i1=find(triu(ones(i2,i2)));r1=zeros(i2,i2);r1(i1)=1:length(i1);
    out=struct('ElMap',r1+tril(r1',-1));   out2=[];
    out1=varargin{1};out1(4,:)=2;% Tell of_mk('MatrixIntegration') this is 2d
  else;
   [out,i2,out2]=p_solid('const','quad9',varargin{2:end});
   out1=varargin{1};out1(4,:)=i2; % Tell MatrixIntegration this is 3d 
  end
 elseif comstr(Cam,'node'); out = [1:9];
 elseif comstr(Cam,'prop'); out = [10 11 12];
 elseif comstr(Cam,'dofcall');out=elem0('dofcall'); % variable field elements
 elseif  comstr(Cam,'dof')
   out = [1:9];k=[1;2;3]/100;out=out(ones(3,1),:)+k(:,ones(9,1));
   out=out(:);
 elseif comstr(Cam,'line');   out = [1 5 2 6 3 7 4 8 1];
 elseif comstr(Cam,'face');   out = [1 5 9 8;5 2 6 9;8 9 7 4;9 6 3 7];
 %elseif comstr(Cam,'face');   out = 1:9;
 elseif comstr(Cam,'edge');   
     out = [1 2 5; 2 3 6; 3 4 7; 4 1 8];
 elseif comstr(Cam,'patch');  out = [1 5 2 6 3 7 4 8];
 elseif comstr(Cam,'sci_face'); out = [1 5 9 8;5 2 6 9;6 3 7 9;7 4 8 9];
 elseif comstr(Cam,'parent'); out = 'quad9';
 elseif  comstr(Cam,'flip');   out=[1 4 3 2 8 7 6 5 9]; out1=1:9; 
  % Basic matrix test of the element - - - - - - - - - - - - - - - - - -
 elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,5);

   model=femesh('testquad9');
   if nargin==2; model.pl=varargin{1}; end
   model.il=[110 fe_mat('p_shell','SI',1) 0 0 0  1e-1];
   if ~isempty(strfind(Cam,'2d'));model.il=[110 fe_mat('p_solid','SI',2) 0];end
   [m,k,mdof]=fe_mknl(model,'NoT');out=stack_cell(k,m);
   if nargout==0;disp('TestMat passed');end
 elseif  comstr(Cam,'test'); 
   out='No Test';
 end
 return

end % of standard calls with one input argument

warning('OpenFEM:NotImplemented','quad9 not implemented')
% OBSOLETE of_mk_subs.c elements - - - - - - - - - - - - - - - - - - - - - - 

