function [out,out1,out2]=hexa27(CAM,varargin);

%HEXA27 27-node,, 6 sided volume element supporting 
%       - fully anisotropic elasticity (hyperelasticity beta version)
%       - integration rule selection
%       - pressure, volume, inertial and thermal loads
%       - geometric non linearity
%
%       Elements of the same family
%         hexa8, hexa20, hexa27, 
%         tetra4, tetra10, 
%         penta6, penta15,
%
%       See sdtweb hexa27

%       Etienne Balmes, Jean-Michel Leclere
%       Copyright (c) 2001-2019 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

% standard calls with one input argument
%#ok<*NOSEL,*ASGLU>

if ischar(CAM)
 [CAM,Cam]=comstr(CAM,1);

 if comstr(Cam,'integinfo')
 %constit integ                          ID,pl,il
  [out ,out1,out2] = p_solid('buildconstit 3 27',varargin{:});

 % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 % Matrix assembly in fe_mknl

 elseif comstr(Cam,'matcall'); out='mat_og'; out1=0; % mat_og and non symmetric
 elseif comstr(Cam,'call'); out=elem0('callmat_og');
 elseif comstr(Cam,'rhscall'); out=elem0('rhsmat_og'); % call for load
 elseif comstr(Cam,'groupinit');out=elem0(CAM,'hexa27');
 elseif comstr(Cam,'dofcall');out=elem0('dofcall'); % variable field elements

 % Here one defines the constants needed for element integration
 elseif comstr(Cam,'constants'); 

  if nargin<3; p_solid('constsolid','hexa27',[],[]);return;
  elseif varargin{2}(size(varargin{2},1),1)==-9999; % old of_mk_sub.c elements
    i2=60;i1=find(triu(ones(i2,i2)));r1=zeros(i2,i2);r1(i1)=1:length(i1);
    out=struct('ElMap',r1+tril(r1',-1));   out2=[];
    out1=varargin{1};out1(4,:)=3;% Tell of_mk('MatrixIntegration') this is 3d
  else;
   [out,i2,out2]=p_solid('constsolid','hexa27',varargin{2:end});
   out1=varargin{1};out1(4,:)=i2; % Tell MatrixIntegration this is 3d 
  end

 % end of GroupInit - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 % Basic tests of the element - - - - - - - - - - - - - - - - - - - - - -
 elseif  comstr(Cam,'test'); [CAM,Cam] = comstr(CAM,5);
   [out,out1]=femesh(strcat(['teststruct hexa27' Cam]));
 % end Basic tests of the element - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'node');  out = 1:27;
 elseif  comstr(Cam,'prop');  out = [28 29 30];
 elseif comstr(Cam,'dof')
   out =(1:27);k=(1:3)'/100;
   out=out(ones(3,1),:)+k(:,ones(27,1));out=out(:);
 elseif comstr(Cam,'line');
   out = ...
    [1 9 2 10 3 11 4 12 1 13 5 17 6 18 7 19 8 20 5 0 2 14 6 0 3 15 7 0 4 16 8];
 elseif comstr(Cam,'patch');
   out = [  4 11 3 10 2 9 1 12;
            5 13 1 9  2 14   6 17;
            5 17 6 18 7 19 8 20;
            4 16 8 19 7 15 3 11;
            5 20 8 16 4 12 1 13;
            2 10 3 15 7 18 6 14];

 elseif  comstr(Cam,'edge');  out = [1 2 9; 2 3 10; 3 4 11; ...
				     4 1 12; 1 5 13; 2 6 14; ...
				     3 7 15; 4 8 16; 5 6 17; ...
				     6 7 18; 7 8 19; 8 5 20];
 elseif  comstr(Cam,'facelin'); % first order face split
   out = [1 12 21 9; 12 4 11 21;9 21 10 2;21 11 3 10
          1 13 22 12;13 5 20 22;12 22 16 4;22 20 8 16;
          1 9 23 13;9 2 14 23;13 23 17 5;23 14 6 17
          5 17 24 20;17 6 18 24;20 24 19 8;24 18 7 19
          2 10 25 14;10 3 15 25;14 25 18 6;25 15 7 18
          4 16 26 11;16 8 19 26;11 26 15 3;26 19 7 15
          ];
 elseif  comstr(Cam,'vollin'); % first order volume split
    out=[1 9 21 12 13 23 27 22;9 2 10 21 23 14 25 27
        21 10 3 11 27 25 15 26;12 21 11 4 22 27 26 16
        13 23 27 22 5 17 24 20;23 14 25 27 17 6 18 24
        27 25 15 26 24 18 7 19;22 27 26 16 20 24 19 8];
    

 elseif  comstr(Cam,'face'); 
     out = [1 4 3 2 12 11 10 9  21 ; ...
				      1 5 8 4 13 20 16 12 22 ; ...
				      1 2 6 5  9 14 17 13 23 ; ...
			              5 6 7 8 17 18 19 20 24 ; ...
				      2 3 7 6 10 15 18 14 25 ; ...
				      3 4 8 7 11 16 19 15 26 ];
 elseif comstr(Cam,'sci_face');
   out = [1 9 13 13;9 2 14 14;14 6 17 17;17 5 13 13;9 14 17 13;
          2 10 14 14;10 3 15 15;15 7 18 18;18 6 14 14;10 15 18 14;
	  6 18 17 17;18 7 19 19;19 8 20 20;20 5 17 17;18 19 20 17;
	  1 12 9 9;12 4 11 11;11 3 10 10;10 2 9 9;12 11 10 9;
	  4 16 11 11;16 8 19 19;19 7 15 15;15 3 11 11;16 19 15 11;
	  1 13 12 12;13 5 20 20;20 8 16 16;16 4 12 12;13 20 16 12];

 elseif  comstr(Cam,'parent');out = 'hexa27';
 elseif  comstr(Cam,'state');out='';% State should be handled by of_mk
 elseif comstr(CAM,'cvs')
 out='$Revision: 1.11 $  $Date: 2020/06/19 12:41:39 $'; return;
 else; error('%s unknown',CAM);
 end
return
end % of standard calls with one input argument

