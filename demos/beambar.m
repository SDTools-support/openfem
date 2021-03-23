
%
% This example illustrates the use of mixed element types with FEMESH
%

femesh('reset'); % initialize FEMESH


FEnode = [1  0 0 0  0 0 0;2  0 0 0  1 0 0;
          3  0 0 0  0 1 0;4  0 0 0  1 1 0];
FEelt=[];

% make an outer frame with beam elements

femesh(';object beamline 1 2 0 4 3 1;addsel;info FEelt');

% make diagonals with bar elements

femesh(';object beamline 1 4 0 2 3;addsel;set group2 name bar1');

% repeat the cell 10 times and keep the result

femesh(';selgroup1:2;repeatsel 10 1 0 0;');
FEelt=FEel0;

% add a beam at the end of the last bay

femesh('objectbeamline',femesh('findnode x==10'));
femesh addsel;


% Define the properties of both groups - - - - - - - - - - - - - - -


% Set the MatID and SecID values for elements

femesh(';setgroup1 3 mat1sec1;setgroup2 mat2 sec2');

% compute modes and display the result
model=femesh('model');
%           MatId   MatType                 E        nu    rho
model.pl=[  1   fe_mat('m_elastic','SI',1)  72e+9    0.3   2700 ;
            2   fe_mat('m_elastic','SI',1)  210e+9   0.3   7800 ];

typ=fe_mat('p_beam','SI',1);
%        SecId SecType  J       I1      I2      A (beam)
model.il=[1     typ        5e-09   5e-09   5e-09   2e-05   % longerons/battens
%        SecId SecType  0       0       0       A (bar)
         2     typ        0       0       0       1e-5];
model=fe_case(model,'fixdof','2D',[.03;.04;.04]);
[m,k,mdof]  = fe_mk(model);
[mode,freq] = fe_eig(m,k,[4 10 1e3]);
femesh plotelt;
feplot('initdef',mode,mdof,freq/2/pi);
fecom(';scaledef1;view2;ch4');

%-----------------------------------------------------------------

%	Etienne Balmes  07/14/93, 07/13/98
%       Copyright (c) 1990-2004 by SDTools,
%       All Rights Reserved.

