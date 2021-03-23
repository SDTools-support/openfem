%=========================================================%
%                       D_TRUSS                           %
%=========================================================%


%---------------------------------------------------------%
% 1. Geometry declaration with femesh                     %
% See section 2.1.2 of the tutorial                       %
%---------------------------------------------------------%
femesh('test2bay')
femesh('removeelt group2');
femesh('divide group 1 InNode 1 4')
femesh('set group1 name bar1');
femesh(';selgroup2 1;repeatsel 10 1 0 0;addsel');
femesh(';rotatesel 1 60 1 0 0;addsel;')
femesh(';selgroup3:4;rotatesel 2 -60 1 0 0;addsel;')
femesh(';selgroup3:8');
% export FEnode and FEel0 geometry in model
model=femesh('model0');  
mflag=~system('medit');
if mflag medit('write truss',model); 
else     femesh('plotel0');
end


%---------------------------------------------------------%
% 2. Handling material and element properties             %
% See section 2.1.4 of the tutorial                       %
%---------------------------------------------------------%
model.pl=[ 1   fe_mat('m_elastic','SI',1)  7.2e+10  0.3   2700 ];
model.il = [ ...
% ProId SecType                 J      I1     I2       A
       1 fe_mat('p_beam','SI',1) 5e-9   5e-9   5e-9   2e-5  0 0 % longerons
       p_beam('dbval 2','circle 4e-3') % circular section 4 mm
       p_beam('dbval 3','rectangle 4e-3 3e-3') %rectangular section 4 x 3 mm
];


%---------------------------------------------------------%
% 3. Boundary conditions and constraints                  %
% See section 2.2.1 of the tutorial                       %
%---------------------------------------------------------%
i1 = femesh('findnode x==10');
model = fe_case(model,'FixDof','base',i1);


%---------------------------------------------------------%
% 4. Assembly                                             %
% See section 2.3.1 of the tutorial                       %
%---------------------------------------------------------%
model = fe_mk(model);


%---------------------------------------------------------%
% 5. Normal modes                                         %
% See section 2.3.3 of the tutorial                       %
%---------------------------------------------------------%
[md1,f1] = fe_eig(model.K{1},model.K{2},[1 4 0 11]);
def = struct('def',md1,'data',f1,'DOF',model.Stack{3}.DOF);
if mflag medit('write visu/truss',model,def,3);
else     feplot(model,def); fecom('ch4');
end