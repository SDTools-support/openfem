

% Generate the model

pl=m_elastic('dbval 100 steel');
il=p_solid('dbval 110 d3');
femesh('reset');
femesh('testhexa8 divide10 10 10');
%femesh('divide2 2 2');
femesh hexa2tetra;
femesh('set groupa1 name tetra4 matid100 proid110');


%mdof=feutil('getdof',FEel0);
%FEel0=feutil('orient',FEnode,FEel0);
%adof=fe_c(mdof,femesh('findnode x==0'),'dof',2);
%[m,k,mdof]=fe_mk(FEnode,FEel0,pl,il,[],adof,[0 1]);


model=struct('Node',FEnode,'Elt',FEel0,'pl',pl,'il',il);
data=struct('sel','groupall','dir',[0 0 9.81]);
model=fe_case(model,'fixdof','base','x==0', ...
    'FVol','Gravity',data);

if ~isempty(gcbf)&strcmp(get(gcbf,'tag'),'test_sdt')
 sdtw('No need to do performance test every day');
 return
else
 profile on
 timing(0); [Case,model.DOF]=fe_mknl('init',model); 
 timing(1);k=fe_mknl('assemble',model,Case,1);timing(2);
 %profile on;k=fe_mknl('assemble',model,Case,1);profile report
 [k1,mdof]=fe_mk(model,'options',1);timing(3);
 timing % Preprocess, NL reassembly, FE_MK assembly
 try; profile report;
 catch;disp('profile failed');
 end
end

Case=fe_case(model,'gett');
def = fe_load(model,Case);
%profile on;def = fe_load(model,Case1);profile report


ofact('method umfpack')
kd=ofact(k);
def.def=kd\def.def;
feplot(model,def);
