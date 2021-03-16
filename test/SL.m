%
% Scordelis-Lo roof test-problem with MITC4 elements
%

%ofact('method umfpack');

disp(' ');
disp('Scordelis-Lo roof test-problem');
disp('Reference value (defl. point B comp. with Modulef): -3.58893');
disp(' ');

a=which('fe_mk');a=a(1:end-7);fname=fullfile(a,'test','sl.nopo');   
model=nopo(sprintf('read -p COQUE %s',fname));
model.pl=[100 fe_mat('m_elastic',1,1) 3.e6 0 0 0];

il5=[110 fe_mat('p_shell',1,1) 1 -1 1 3.]; %5 dofs
il6=[110 fe_mat('p_shell',1,1) 1 1 1 3.];  %6 dofs

femesh(model);
femesh('set group1 matid100 proid110');
femesh('sel group1');

model5=femesh('model0');
model5.il=il5;

data=struct('sel','groupall','dir',[0 0 -.625]);
model5=fe_case(model5,'FVol','Gravity',data);

model6=model5; model6.il=il6;

model5 = fe_case(model5,'fixdof','Ref5','y==0 -DOF 2 4');
model5 = fe_case(model5,'fixdof','Ref3','x==300 -DOF 2 3');
model5 = fe_case(model5,'fixdof','Ref4','x==0 -DOF 1 5');
model6 = fe_case(model6,'fixdof','Ref5','y==0 -DOF 2 5');
model6 = fe_case(model6,'fixdof','Ref3','x==300 -DOF 2 3');
model6 = fe_case(model6,'fixdof','Ref4','x==0 -DOF 1 4');

model5=fe_case(model5,'fixdof','DDL6',.06);

[Case,model5.DOF]=fe_mknl('init',model5);
% Define normals manually
mitc4normals=[zeros(size(FEnode,1),1) FEnode(:,6)./300 FEnode(:,7)./300]';
Case.GroupInfo{1,7}.normal=mitc4normals;
NNode=sparse(model5.Node(:,1),1,1:size(model5.Node,1));
Case.GroupInfo{1,7}.EltInd=full(NNode(model.Elt(2:end,1:4))');
k=fe_mknl('assemble',model5,Case,1);
def = fe_load(model5,Case);
defor5=ofact(k,def.def);
disp(' ');
fprintf('\n Result with 5 dofs: %15.5f\n',full(fe_c(Case.DOF,17.03)*defor5));

[Case,model6.DOF]=fe_mknl('init',model6);
Case.GroupInfo{1,7}.normal=mitc4normals;
Case.GroupInfo{1,7}.EltInd=full(NNode(model.Elt(2:end,1:4))');
k=fe_mknl('assemble',model6,Case,1);
def = fe_load(model6,Case);
kd=ofact(k); defor6=kd\def.def;ofact('clear',kd);
fprintf('\n Result with 6 dofs: %15.5f\n',full(fe_c(Case.DOF,17.03)*defor6));

