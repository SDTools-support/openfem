% test TRIA3 element :
%
% test matrices, different formulations, degenerated elts, RHS 
%
%
%

 if ~isempty(gcbf)&&strcmp(get(gcbf,'tag'),'test_sdt');
  wd=pwd;
  cd(fileparts(fileparts(which('fe_mk'))));eval('sdtcheck path');
  clear functions variables global;comgui('close all');cinguj('initSwing');
 end



% --------------------------------------------------- eigenvalue tests
femesh('reset');
model=femesh('testquad4');
femesh(';testquad4 divide6 6;quad2tria;set groupa1 name tria3 matid100 proid110');
FE.pl=m_elastic(model.pl,'dbval 99 Air');FE.il=p_solid(model.il,'dbval 99 Full');
FE.pl=model.pl(end:-1:1,:);FE.il=model.il(end:-1:1,:);
model=struct('Node',FEnode,'Elt',FEel0,'pl',FE.pl,'il',FE.il);
[i1,i2]=find(model.il==110);

freq=[];
for j2=-1:0 % drilling

  model.il(i1,4)=j2;
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
  if j2==-1
   def=fe_eig(fe_case(model,'fixdof','drill',.06),[105 10 1e3 11]);
  else
   def=fe_eig(model,[105 10 1e3 11]);
  end
  freq(end+1,1:7)=def.data(1:7)';

end % j2
freq % displays 7 first frequencies


% --------------------------------------------------- load tests
% load
femesh('reset');
model=femesh('testquad4');
femesh(';testquad4 divide6 6;quad2tria;set groupa1 name tria3 matid100 proid110');
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
  [Case,model.DOF]=fe_mknl('init',model);
  [m,k,mdof]=fe_mk(model);
  k1=fe_mknl('assemble',model,Case,1);  
  if any(size(k)-size(k1)) 
    ind=fe_c(Case.DOF,.06,'ind',2); 
    k2=k1(ind,ind); 
  else k2=k1;ind=fe_c(Case.DOF,mdof,'ind');
  end
  if norm(diag(k)-diag(k2))/norm(diag(k))>sqrt(eps) 
   error('Problem in NL assembly');
  end
  def = fe_load(model);
  def.def(ind,:)=ofact(k,def.def(ind,:));
  feplot(model,def);

 end % j2
end %j1

% --------------------------------------------------- 
