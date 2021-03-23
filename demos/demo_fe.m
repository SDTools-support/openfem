%=========================================================%
%                       DEMO_FE                           %
%=========================================================%


%---------------------------------------------------------%
% 1. Geometry declaration with femesh                     %
% See section 3.1.2 of the tutorial                       %
%---------------------------------------------------------%
clear model
femesh('reset');
FEelt=[];
FEnode=[1 0 0 0  0 0 0;2 0 0 0    0 1 0;
        3 0 0 0  1 0 0;4 0 0 0    1 1 0];
femesh('objectbeamline 1 3 0 2 4 0 3 4 0 1 4')

femesh(';addsel;transsel 1 0 0;addsel;info');
% export FEnode and FEelt geometry in model
model=femesh('model');  
feplot(model);
fecom('view2');
pause;

%---------------------------------------------------------%
% 2. Handling material and element properties             %
% See section 3.1.4 of the tutorial                       %
%---------------------------------------------------------%
model.pl=[ 1   fe_mat('m_elastic','SI',1)  7.2e+10  0.3   2700 ];
model.il = [ ...
% ProId SecType                 J      I1     I2       A
       1 fe_mat('p_beam','SI',1) 5e-9   5e-9   5e-9   2e-5  0 0 % longerons
       p_beam('dbval 2','circle 4e-3') % circular section 4 mm
       p_beam('dbval 3','rectangle 4e-3 3e-3') %rectangular section 4 x 3 mm
];

mpid = feutil('mpid',model.Elt);
mpid(:,2) = [0 1 1 2 3 0 1 1 2 3]';
model.Elt = feutil('mpid',model.Elt,mpid);


%---------------------------------------------------------%
% 3. Boundary conditions and constraints                  %
% See section 3.2.1 of the tutorial                       %
%---------------------------------------------------------%
model=fe_case(model,'SetCase1', ...         % defines a new case
  'FixDof','2-D motion',[.03 .04 .05]', ... 
  'FixDof','Clamp edge',[1 2]');  


%---------------------------------------------------------%
% 4. Assembly                                             %
% See section 3.3.1 of the tutorial                       %
%---------------------------------------------------------%
model = fe_mk(model);
model.DOF = model.Stack{3}.DOF;


%---------------------------------------------------------%
% 5. Normal modes                                         %
% See section 3.3.3 of the tutorial                       %
%---------------------------------------------------------%
[md1,f1] = fe_eig(model.K{1},model.K{2},[1 4 0 11]);
def = struct('def',md1,'data',f1,'DOF',model.DOF);


%---------------------------------------------------------%
% 6. Visualization of deformed structure                  %
% See section 3.4 of the tutorial                         %
%---------------------------------------------------------%

if ~system('medit'); medit('write visu/fe',model,def);
else; feplot(model,def); fecom('ch3');
end
