function [out,out1,out2]=penta18(CAM,varargin);

%PENTA18 element function for the 18-node 36-DOF pentaedral 3D element (C3D15V)
%
% Abaqus element : display only (OpenFEM numbering)
%
% Abaqus description :
% Stress/displacement variable node elements
% Face node number	Corner nodes on face
% 16	1 - 4 - 5 - 2
% 17	2 - 5 - 6 - 3
% 18	3 - 6 - 4 - 1
%
% 	See sdtweb      eltfun, elem0
%	See also help   hexa20, hexa8, tetra4, penta6

%	Etienne Balmes, Jean-Michel Leclere, Marina Vidrascu  
%       Copyright (c) 2001-2020 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

if comstr(CAM,'cvs')
 out='$Revision: 1.3 $  $Date: 2020/02/26 08:47:25 $'; return;
end
% standard calls with one input argument
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1);
 
 if  comstr(Cam,'node');  out = [1:18];
 elseif  comstr(Cam,'prop');  out = [19 20 21];
 elseif  comstr(Cam,'dof')
    out=[1:18];k=[1;2;3]/100;out=out(ones(3,1),:)+k(:,ones(18,1));
    out=out(:);
 elseif  comstr(Cam,'line');   [1 7 2 8 3 9 1 0 1 10 4 0 ...
                                       2 11 5 0 3 12 6 0 4 13 5 14 6 15 4];
 elseif  comstr(Cam,'patch');  
             out = [1 9 3 8 2 7 7 7 7; 
                    1 10 4 15 6 12 3 9 18; 
                    1 7 2 11 5 13 4 10 16; 
                    4 13 5 14 6 15 15 15 15; 
                    2 8 3 12 6 14 5 11 17];  
 elseif  comstr(Cam,'edge');   
      out = [1 2 7; 2 3 8; 3 1 9; ...
             1 4 10; 2 5 11; 3 6 12; ...
             4 5 13; 5 6 14; 6 4 15];
 elseif  comstr(Cam,'face');   
          out = [1 3 2 9 8 7 7 7 7; ...
                 1 4 6 3 10 15 12 9 18; ...
                 1 2 5 4 7 11 13 10 16; ...
                 4 5 6 13 14 15 15 15 15; ...
                 2 3 6 5 8 12 14 11 17];
 elseif  comstr(Cam,'flip');   out=[1 3 2 4 6 5 9 8 7 10 12 11 15 14 13 18 17 16]; 
                              out1=1:18; 

 elseif  comstr(Cam,'parent');   out = 'penta15'; 

 % Basic tests of the element - - - - - - - - - - - - - - - - - - - - - -
 else sdtw('''%s'' unknown',CAM);  end

return
end % of standard calls with one input argument



