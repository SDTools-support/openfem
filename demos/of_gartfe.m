
% Demonstration of the finite element mesh building capabilities of the SDT
% for the case of the GARTEUR SM-AG-19 Testbed.
%
% See also demos d_truss, demo_fe, beambar, d_plate, d_ubeam
%          doc   fem, dfeplot

femesh('reset');

FEelt=[];
FEnode = [1 0 0 0  0 0 0;2 0 0 0  0 0 .15;
          3 0 0 0 0.4 1.0 .176;4 0 0 0 0.4 0.9 0.176];
% fuselage 
femesh('objectbeamline 1 2');
femesh('extrude 0  1.0 0.0 0.0',[linspace(0,.55,5) linspace(.65,1.4,6) 1.5]);
femesh('addsel;');

% vertical tail
femesh('objectbeamline',femesh('findnode z==.15 & x>=1.4'));
femesh(';extrude 3 0 0 .1;addsel;');
% vertical horizontal tail
femesh('objectbeamline',femesh('findnode z==.45'));
femesh('extrude 0  0.0 0.2 0.0',[-1 -.5 0 .5 1]);
femesh('addsel;');

% right drum
femesh(';objectbeamline 3 4;extrude 1 .4 0 0');
femesh('divide',[0 2/40 15/40 25/40 1],[0 .7 1]);
femesh('addsel;');

% left drum
femesh(';symsel 1 0 1 0;addsel;');

% wing

femesh('objectbeamline',femesh('findnode y==1 & x>=.55 & x<=.65'));
femesh('divide',[0 1-.762 1]);
femesh('extrude 0  0.0 -1.0 0.0',[0 0.1 linspace(.15,.965,9) ...
     linspace(1.035,1.85,9) 1.9 2.0]);
femesh('addsel;');

% connection plate

femesh('objectbeamline',femesh('findnode y==0.035 | y==-0.035 & x==.55'));
femesh(';divide 2;transsel -.02 0 0;');
femesh('extrude 0 1 0 0',[0 .02 .12 .14]);
i1 = femesh('findnode group6 & groupa1');
femesh(';transsel 0.0 0.0 -0.026;addsel;');

% stiff links for the connection
femesh('object mass',i1);
femesh('extrude 1 0 0 -.026');
femesh('set groupa 1 name celas');
% set connected DOFs and spring value
FEel0(2:end,3)=12345;FEel0(2:end,4)=0;FEel0(2:end,7)=1e12;
femesh('addsel;');

% make a group of the part covered by the constraining layer
femesh('divide group 6 InNode {x>.55 & y<=.85 & y>=-.85}');

% tip masses
femesh('object mass',femesh('findnode y==0.93 | y==-0.93& x==0.42'));
femesh('addsel;');FEelt(femesh('findeltgroup10'),2:4)= 0.2; 
femesh('object mass',femesh('findnode z==.45 & y==0'));
femesh('addsel;');FEelt(femesh('findeltgroup11'),2:4)= 0.5; 
femesh(';join mass1;');

% orient plates that will need an off-set
femesh('orient 4:8 n 0 0 3');

FEelt(femesh('findeltgroup4:5'),9)= 0.005; % drums (positive off-set)
FEelt(femesh('findeltgroup6:7'),9)=-0.005; % wing 
FEelt(femesh('findeltgroup  8'),9)= 0.008; % wing 

femesh('plotelt'); fecom(';sub 1 1;view3');


pause;%------------------------------------------------

femesh(';set group1 mat1 pro3;set group2:7 mat1 pro1;');
femesh(';set group8 mat2 pro2;setgroup6 pro4;');
femesh('selgroup1:10');

model=femesh('model');
model.pl=[m_elastic('dbval 1 aluminum');
    m_elastic('dbval 2 steel')];

model.il = [1 fe_mat('p_shell','SI',1)     2 1 0     .01
      2 fe_mat('p_shell','SI',1)     2 1 0     .016
      3 fe_mat('p_shell','SI',1)     2 1 0     .05
      4 fe_mat('p_shell','SI',1)     2 1 0     .011];


[m,k,mdof] = fe_mk(model);
[md0,f0] = fe_eig(m,k,[5 20 1e2 11]);f0=f0/2/pi;
feplot(FEnode,FEelt,md0,mdof,2);
fecom(';ch 7;colordatagroup;showpatch');
feplot

%-----------------------------------------------------------------

%	Etienne Balmes  07/14/93, 07/13/00
%       Copyright (c) 1990-2004 by SDTools,
%       All Rights Reserved.

