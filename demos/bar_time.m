
% simple test : bar in traction-compression

 FEnode=[1 0 0 0    0   0  0
         2 0 0 0    50  0  0];
 FEel0=[Inf abs('bar1')];
 FEel0(2,:)=[1 2 1 1 1];
 pl=m_elastic('dbval 1 steel');
 il=[1 1 0 0 0 1];
 femesh(';divide500');

 model=struct('Node',FEnode,'Elt',FEel0,'pl',pl,'il',il);
 data=struct('DOF',2.01,'def',1e6);

 model = fe_case(model,'AddToCase 1','DOFLoad','Point load 1',data);
 model = fe_case(model,'AddToCase 1','KeepDof','kd',.01);
 model = fe_case(model,'AddToCase 1','FixDof','fd',1);

 model=fe_mk(model);
 mdof=model.DOF;

 q0=[];
 com.Method='newmark';
 com.Opt=[.25 .5 3e-4 1e-4 100 10];
 com.Residual='';
 com.OutputFcn='';
 com.NeedUVA=[1 1 0];

 def1=fe_time(com,model,'Case 1',q0);

 disp('With Newmark')
 if sp_util('issdt')
  cf=feplot;cf.model=model;cf.def={def1.v,def1.DOF+.01};
  fecom(';view2;scd1e1;animtime;ch99');
 else
  feplot(model.Node,model.Elt,def1.v,def1.DOF+.01,2)
  fecom(';view2;ch99;scaledef1e1');
 end

 disp('--------------------------------------------------------pause')
 pause;

 com.Method='dg'; model=stack_set(model,'info','q0',q0);
 def2=fe_time(com,model);

 disp('With DG')
 if sp_util('issdt')
  cf=feplot;cf.model=model;cf.def={def2.v,def2.DOF+.01};fecom(';view2;scd1e1;animtime');
 else
  feplot('initmodel',model);feplot('initdef',def2.v,def2.DOF+.01);
  fecom(';view2;ch99;scaledef1e1');
 end

%------------------------------------------------------------------------

