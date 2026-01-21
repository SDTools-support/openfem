
close('all');clear('all');
if exist('sdth','file');ofutil('path');end

% this example illustrates the use of the FEMESH preprocessor to build a
% solid model of a U-beam, compute the associated modes, and display strain
% energy levels

femesh('reset');
FEnode=[1 0 0 0  -.5 -.5 0;2  0 0 0  -.5+1/6 -.5 0;3 0 0 0  -.5 .5-1/6 0
        4 0 0 0  -.5+1/6 .5-1/6 0;5 0 0 0  -.5 .5 0;6 0 0 0 -.5+1/6 .5 0
        7 0 0 0 .5-1/6 .5 0;8 0 0 0 .5 .5 0;9 0 0 0 .5-1/6 .5-1/6 0
        10 0 0 0 .5 .5-1/6 0;11 0 0 0 .5-1/6 -.5 0;12 0 0 0 .5 -.5 0];

FEelt=[Inf abs('quad4');4 6 5 3 1 1;9 10 8 7 1 1];
FEel0=[Inf abs('quad4');1 2 4 3 1 1];
femesh(';divide 5 1;addsel;');
FEel0=[Inf abs('quad4');11 12 10 9 1 1];
femesh(';divide 5 1;addsel;');
FEel0=[Inf abs('quad4');4 6 7 9 1 1];
femesh(';divide 4 1;addsel;');
femesh('join group 1:4');

femesh(';selgroup1;extrude 10 0 0 .25;orientel0');

% This section is to impose a cantilevered boundary condition

model=femesh('model0');
model.pl =  m_elastic('dbval 1 steel');
model=p_solid('default',model); % Define integration rule
model=fe_case(model,'fixdof','base','z==0');


% This section is to assemble, compute modes and show the result

[m,k,mdof]=fe_mk(model);
[md1,f1] = fe_eig(m,k,[4 10 0 11]);
def=struct('def',md1,'DOF',mdof,'data',f1/2/pi)
feplot(model.Node,model.Elt,md1,mdof,1); axis auto

StrainEnergy = fe_stress('ener',model,def);
fecom(';color face flat;color edge w');
feplot('init cdef',StrainEnergy);


disp('pause');pause

% Now we will apply a load on the edge

data=struct('sel','x==-.5', ... 
    'eltsel','withnode {z>1.25}','def',1,'DOF',.19);
model=fe_case(model,'Fsurf','Surface load',data);


%view load
Load = fe_load(model); 
feplot(FEnode,FEel0,Load.def,Load.DOF,2);
fecom('view3');axis auto


disp('pause');pause


%view response
Case=fe_case(model,'gett');

def.def=ofact(k,Load.def); feplot(FEnode,FEel0,def,mdof,1);
Stress = fe_stress('stress mises',model,def);
feplot(FEnode,FEel0,def,mdof,1,Stress.data); axis auto

feplot('init cdef',Stress)
fecom(';color face interp;color edge k');
axis auto

disp('pause');pause

% Loads
% See section 3.2.2 of the tutorial   

% Volume forces
data  = struct('sel','GroupAll','dir',[1 0 0]);
model = fe_case(model,'FVol','Volume load',data);
Load  = fe_load(model);
feplot(model.Node,model.Elt,Load.def,Load.DOF,2);

disp('pause');pause

% Surfacic forces
data=struct('sel','x==-.5', ... 
   'eltsel','withnode {z>1.25}','def',1,'DOF',.19);
model=fe_case(model,'Fsurf','Surface load',data);
Load = fe_load(model); 
feplot(model,Load);

disp('pause');pause
 

% 2 loads
data  = struct('DOF',[207.01;241.01;207.03],'def',[1 0;-1 0;0 1]);
model = fe_case(model,'reset','DOFLoad','Point load 1',data);
data  = struct('DOF',365.03,'def',1);
model = fe_case(model,'DOFLoad','Point load 2',data);
Load  = fe_load(model);
feplot(model,Load);


%       Etienne Balmes
%       Copyright (c) 1996-2007 by SDTools
%       All Rights Reserved.

