function [out,out1,out2]=pyra13(CAM,varargin)%(node,elt,pl,il,opt,NNode)

%PYRA13	13-node pyramid volume element
%
% 	See sdtweb     eltfun, elem0
%	See also help  hexa20, hexa8, penta6


%       Etienne Balmes
%       Copyright (c) 2001-2019 by INRIA and SDTools
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       All Rights Reserved.

% standard calls with one input argument
if ischar(CAM)
 [CAM,Cam]=comstr(CAM,1);
 if comstr(Cam,'integinfo')
  %constit integ,elmap                 ID,pl,il
  [out,out1,out2]=p_solid('buildconstit 3 13',varargin{:});
 elseif comstr(Cam,'matcall'); [out,out1]=elem0(CAM,varargin{:});
 elseif comstr(Cam,'dofcall');out=elem0('dofcall'); % variable field elements
 elseif comstr(Cam,'call'); out=elem0('callmat_og');
 elseif comstr(Cam,'rhscall'); out=elem0('rhsmat_og'); % call for load
 elseif comstr(Cam,'groupinit');out=elem0(CAM,'pyra13');
 elseif comstr(Cam,'constants')
  if nargin<3;error('TexStress now handled in integrules');
  elseif varargin{2}(size(varargin{2},1),1)==-9999 % old of_mk_sub.c elements
    i2=24;i1=find(triu(ones(i2,i2)));r1=zeros(i2,i2);r1(i1)=1:length(i1);
    out=struct('ElMap',r1+tril(r1',-1));   out2=[];
    out1=varargin{1};out1(4,:)=3;% Tell of_mk('MatrixIntegration') this is 3d
  else
   [out,i2,out2]=p_solid('const','pyra13',varargin{2:end});
   out1=varargin{1};out1(4,:)=i2; % Tell MatrixIntegration this is 3d 
  end
 elseif comstr(Cam,'node');   out = 1:13;
 elseif comstr(Cam,'dof');out=feutil('getdof',(1:13)',(1:3)'/100);
 elseif comstr(Cam,'line');    out = [1 6 2 7 3 8 4 9 1 10 5 13 4 0 3 12 5 11 2];
 elseif  comstr(Cam,'edge');   out = [1 2 6; 2 3 7; 3 4 8; 4 1 9; 1 5 10; 2 5 11;3 5 12;4 5 13];
 elseif  comstr(Cam,'face')
     out = [1 4 3 2 6 7 8 9; 4 1 5 9 10 13 13 13; 
            1 2 5 6 11 10 10 10; 2 3 5 7 12 11 11 11;3 4 5 8 13 12 12 12];
 elseif  comstr(Cam,'patch')
     out = [1 9 4 8 3 7 2 6;4 9 1 10 5 5 5 13; 
            1 6 2 11 5 5 5 10;2 7 3 12 5 5 5 11;3 8 4 13 5 5 5 12];
 elseif comstr(Cam,'sci_face');  out = [1 4 3 2;4 1 5 5; 1 2 5 5; ...
                                      2 3 5 5; 3 4 5 5];                             
 elseif  comstr(Cam,'parent'); out = 'pyra13';
 elseif  comstr(Cam,'prop');  out = 14:16;
 elseif  comstr(Cam,'flip'); out=[1 4 3 2 5 9 8 7 6 10 13 12 11];out1=1:13; % negative volume flip

 elseif  comstr(Cam,'testmat')
   EC=integrules('pyra13');model=struct('Node',[(1:13)'*[1 0 0 0] EC.xi],'Elt',[]);
   model.Elt=feutil('addelt','pyra13',[1:13 100 111]);
   model.pl=m_elastic('dbval 100 Steel');
   model.il=p_solid('dbval 111 d3 -3');
   [m,k,mdof]=fe_mknl(model,'NoT');out=stack_cell(k,m); %#ok<*ASGLU>
   
 elseif  comstr(Cam,'test'); [CAM,Cam] = comstr(CAM,5);

   [out,out1]=femesh(strcat(['teststruct pyra13' Cam]));

 elseif comstr(Cam,'cvs')
  out='$Revision: 1.2 $  $Date: 2019/01/29 17:22:04 $'; return;
 else;error('Not an implemented method %s',CAM);
 end
return
end % of standard calls with one input argument

