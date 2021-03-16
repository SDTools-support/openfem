% test_medit.sci
% use :
%  test_medit() : runs tests for medit
%  test_medit('clean') : run tests and clean directory


function test_medit(cam)

if ~isempty(gcbf)&&strcmp(get(gcbf,'tag'),'test_sdt');return;end

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
pl = [1 1 2e11 .30 7800 (190e9/2/(1+.29))];
mdof = femesh('finddof group1',FEel0);
i1 = femesh('findnode z==0');
[mdof] = fe_c(mdof,i1,'dof',2);

% This section is to assemble, compute modes and show the result
[m,k,mdof] = fe_mk(FEnode,FEel0,pl,[],[],mdof,[0 1]);

% Now we will apply a load on the edge

model=struct('Node',FEnode,'Elt',FEel0,'pl',pl,'il',[],'DOF',mdof);

data=struct('sel','x==-.5', ... 
    'eltsel','withnode {z>1.25}','def',1,'DOF',.19);
model=fe_case(model,'Fsurf','Surface load',data);


% view load
Load = fe_load(model); 

%--------------------------------%
%        LOAD ANIMATION          %
%--------------------------------%
%disp('Load animation');
medit('write visu/open',model,Load,'a',[1 10 3]);
% animation running in medit : press mouse right button, menu Animation,
% submenu Play sequence

% view response
def = k\Load.def;
Stress = fe_stress('stress mises',FEnode,FEel0,pl,[],def,mdof);

%------------------------------------%
%       CONSTRAINT ANIMATION         %
%------------------------------------%
%disp('Constraint animation');
defs = struct('def',def,'DOF',mdof);
medit('write visu/opencol',model,defs,Stress,'a',[1 10 1e8]);
% coloring display in medit : press mouse right button, menu Data, Toggle metric
% animation running in medit : press mouse right button, menu Animation,
% submenu Play sequence

%---------------------------------%
%           MESH DISPLAY          %
%---------------------------------%
disp('Mesh display');
medit('write visu/open_mail',model);

%---------------------------------------%
%       MESH DISPLAY AND COLORING       %
%---------------------------------------%
disp('Colored mesh display');
medit('write visu/opencol_mail',model,[],Stress);
% coloring display in medit : press mouse right button, menu Data, Toggle metric

if nargin>0
	if comstr(cam,'clean')
		system('rm visu/open*');
	end
end

