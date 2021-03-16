

%-- Test example --%

function test_fe_check;

FEnode=[];FEelt=[];
FEnode=[1 0 0 0 -1 0 0;2 0 0 0 1 0 0];
femesh('ObjectBeamLine 1 2');femesh('divide 10');
FEelt=FEel0;
femesh('extrude 10 0 -.2 0');
ind=unique(round(rand(15,1)*(length(FEel0(:,1))-1)))+1;
for i1=ind 
FEel0(i1,1:4)=FEel0(i1,randperm(4)); 
end;

ind=unique(round(rand(30,1)*(length(FEnode)-12)))+12;
dx=(rand(length(ind),1)*2-1)*.1;
dy=(rand(length(ind),1)*2-1)*.1;
FEnode(ind,5:6)=FEnode(ind,5:6)+[dx dy];


femesh('add FEelt FEel0');
femesh('selgroup 1');
femesh('extrude 10 0 .2 0');
femesh('extrude 10 0 0 .2');
FEel0=feutil('orient',FEnode,FEel0);
ind=unique(round(rand(40,1)*(length(FEel0(:,1))-1)))+1;
for i1=ind 
tt=[1 randperm(6)+1 8];
FEel0(i1,1:8)=FEel0(i1,tt); 
end;
femesh('add FEelt FEel0');

ind=unique(round(rand(300,1)*(length(FEnode)-1)))+1;
dx=(rand(length(ind),1)*2-1)*.05;
dy=(rand(length(ind),1)*2-1)*.05;
dz=(rand(length(ind),1)*2-1)*.05;
FEnode(ind,5:7)=FEnode(ind,5:7)+[dx dy dz];

model.Node=FEnode;model.Elt=FEelt;
model_clean=feutil('optim eltcheck',model);
if sp_util('issdt')
 cf=feplot(2);cf.model=model;
 disp('Initial model. Contains twisted / non convex elements.');
 cg=feplot(3);cg.model=model_clean;
end
