function out=t_quad4(varargin)

% test QUAD4 element :
%
% test matrices, different formulations, degenerated elts, RHS 
%
%
%

if nargin==0
 % This test really runs in SDT although it tests OpenFEM
 if ~sp_util('issdt');
  cd(fileparts(fileparts(which('fe_mk'))));eval('sdtcheck path');
  clear variables global;comgui('close all');cinguj('initSwing');
 end
 t_quad4('eig')
 t_quad4('load')
 t_quad4('degen')
 t_quad4('tr')
 t_quad4('integ')
 
 return
end

CAM=varargin{1};
if strncmpi(CAM,'eig',3)
 %% #EIG --------------------------------------------------- eigenvalue tests
femesh('reset');
model=femesh('testquad4 divide6 6');
FE.pl=m_elastic(model.pl,'dbval 99 Air');FE.il=p_solid(model.il,'dbval 99 Full');
FE.pl=model.pl(end:-1:1,:);FE.il=model.il(end:-1:1,:);
model=struct('Node',FEnode,'Elt',FEel0,'pl',FE.pl,'il',FE.il);
[i1,i2]=find(model.il==110);

freq=[];
for j1=0:2 % formulation
 for j2=-1:0 % drilling

  model.il(i1,3)=j1;model.il(i1,4)=j2;
  [m,k,mdof]=fe_mk(model);
  [Case,model.DOF]=fe_mknl('init',model);
  k1=fe_mknl('assemble',model,Case,1);  
  if any(size(k)-size(k1)) 
    ind=fe_c(Case.DOF,.06,'ind',2); 
    k2=k1(ind,ind); 
  else k2=k1;
  end
  if norm(diag(k)-diag(k2))/norm(diag(k))>sqrt(eps) 
   error('Problem in NL assembly');
  end
  if j2~=-1;def=fe_eig(model,[5 30 1e3 11]);
  else; def=fe_eig(fe_case(model,'fixdof','dril',.06),[5 30 1e3 11]);
  end
  r1=sort(def.data);%r1(1:9)
  freq(end+1,1:7)=r1(1:7)';

 end % j2
end %j1
format short e;
freq % displays 7 first frequencies
if norm(freq(:,1:6),'inf')>1;error('Rb problem');end
if norm(freq(:,7)/freq(1,7)-1,'inf')>.02;error('First mode problem');end

elseif strncmpi(CAM,'load',4)
 %% #LOAD --------------------------------------------------- load tests
% load
femesh('reset');
model=femesh('testquad4 divide6 6');
FE.pl=m_elastic(model.pl,'dbval 99 Air');FE.il=p_solid(model.il,'dbval 99 Full');
FE.pl=model.pl(end:-1:1,:);FE.il=model.il(end:-1:1,:);
model=struct('Node',FEnode,'Elt',FEel0,'pl',FE.pl,'il',FE.il);

model=fe_case(model,'fixdof','Edge','x==0');
data=struct('sel','groupall','dir',[0 0 9.81]);
model=fe_case(model,'FVol','Gravity',data);
[i1,i2]=find(model.il==110);

for j1=0:2 % formulation
 for j2=-1:0 % drilling

  model.il(i1,3)=j1;model.il(i1,4)=j2;
  [m,k,mdof]=fe_mk(model);
  [Case,model.DOF]=fe_mknl('init',model);
  k1=fe_mknl('assemble',model,Case,1);  
  m1=fe_mknl('assemble',model,Case,2);  
  if any(size(k)-size(k1)) 
    ind=fe_c(Case.DOF,.06,'ind',2); 
    k2=k1(ind,ind); 
  else k2=k1;
  end
  if norm(diag(k)-diag(k2))/norm(diag(k))>sqrt(eps) 
   error('Problem in NL assembly');
  end
  def = fe_load(model,Case);i1=find(diag(k1)+0);
  def.def(i1,:)=ofact(k1(i1,i1),def.def(i1,:)); 
  feplot(model,def);

 end % j2
end %j1

elseif strncmpi(CAM,'degen',5)
 %% #DEGEN --------------------------------------------------- degenerate elements
femesh('reset');
model=femesh('testquad4');
FE.pl=m_elastic(model.pl,'dbval 99 Air');FE.il=p_solid(model.il,'dbval 99 Full');
FE.pl=model.pl(end:-1:1,:);FE.il=model.il(end:-1:1,:);
femesh('testquad4 divide6 6'); femesh quad2tria;
FEel0(2:end,1:7)=FEel0(2:end,[1 2 3 3 4 5 6]);
femesh('set groupa1 name quad4')
model=struct('Node',FEnode,'Elt',FEel0,'pl',FE.pl,'il',FE.il);

model=fe_case(model,'fixdof','Edge','x==0');
data=struct('sel','groupall','dir',[0 0 9.81]);
model=fe_case(model,'FVol','Gravity',data);

[m,k,mdof]=fe_mk(model);
def = fe_load(model);
kd=ofact(k);def.def=kd\def.def; ofact('clear',kd);
feplot(model,def);


elseif strncmpi(CAM,'tr',2)
 %% #TR --------------------------------------------------- 
femesh('reset');
model=femesh('testquad4');
model.Node(:,5:7)=[0 .1 0; 1 -.1 0;1 0 1;0 0 1]; offset=.0;
[tr,x]=sp_util('quad4tr',model.Node(:,5:7),offset);
cf=feplot(2);cf.model=model;

d1=struct('def',zeros(24,1),'DOF',quad4('dof'));
%d1.def(4:6:24,1)=1;d1.def(3:6:24,1)=node(:,2);d1.def(2:6:end)=-node(:,3);
d1.def([1 19],1)=1;
cf.def=d1;fecom(cf,'undef');

smap=fevisco('strainmap',model,d1);

m2=struct('Elt',model.Elt,'Node',model.Node);m2.Node(:,5:7)=x;m2.Node(:,7)=0;
cg=feplot(3);cg.model=m2;cg.def={tr*d1.def,d1.DOF};fecom(cg,'undef');


elseif strncmpi(CAM,'integ',5) % --------------------------------------------------- 
%% #Integ: test the 'integinfo' calls
femesh('reset');model=femesh('testquad4');
model.Elt(1,1:6)=[Inf abs('q4p') 0 0];
fe_mknl('init',model);

model=femesh('testquad9');fe_mknl('init',model);

%% ---------------------------------------------------------------------------
elseif strncmpi(CAM,'cvs',3)
 out='$Revision: 1.17 $  $Date: 2021/12/13 18:20:08 $';
else; error('%s unknown',CAM);
end


