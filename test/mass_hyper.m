format short e;

ElemF='hexa20b';
%ElemF='tetra4b';

div=3;
nrand=0; L=[.1 .2 .3];% Lambda values (see theoretical description)
%nrand=.2;L=L/100; % randomize node positions to check element consistence

r1=0;for ji=1:3;for jj=1:ji-1; r1=r1+(1+L(ji))^2*(1+L(jj))^2;end;end
I=[sum((1+L).^2) r1 prod((1+L).^2)];

[dWdI,d2WdI2]=feval(elem0('@EnHeart'),[],[],I);
for ji=1:3; 
    r1(ji)=2*(dWdI(1)+dWdI(2)*(I(1)-(1+L(ji))^2) + dWdI(3)*I(3)*(1+L(ji))^-2);
end
Sigma=diag(r1);

femesh('reset');
model=femesh(horzcat('teststruct',ElemF, ...
  sprintf(' divide %i %i %i rand %.15g back',div,div,div,nrand)));
model.Node(:,7)=model.Node(:,7)/max(abs(model.Node(:,7)));

% Boundary conditions
model=fe_case(model,'reset','FixDof','x','x==0 -dof 1', ...
   'FixDof','y','y==0 -dof 2','FixDof','z','z==0 -dof 3');
model.Elt(2:end,length(feval(ElemF,'node'))+1)=100; % all elements with MatId 100
model.pl=[100 fe_mat('type','m_hyper','SI',2) 3 .2 2. 1e-6];
%model.pl=m_hyper('dbval 100 Ref2'); % this is where the material is defined
%model.pl=m_hyper('dbval 100 mooney'); % this is where the material is defined


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Now compute loads and start a newton

% Verify the load based on a surface computation
% (tracking load not implemented yet)
data=struct('sel','x==1','def',(1+L(1))*Sigma(1,1),'DOF',.19);
model=fe_case(model,'FSurf','xs',data);
data=struct('sel','y==1','def',(1+L(2))*Sigma(2,2),'DOF',.19);
model=fe_case(model,'FSurf','ys',data);
%data=struct('sel','z==1','def',0.05,'DOF',.19);
data=struct('sel','z==1','def',(1+L(3))*Sigma(3,3),'DOF',.19);
model=fe_case(model,'FSurf','zs',data);
Load=fe_load(model); F1=sum(Load.def,2);
% save('fp.mat','F1'); 
if sdtdef('isinteractive');feplot(model);end

[Case,model.DOF]=fe_mknl('init',model);

%F1 = zeros(size(Case.DOF));
%F1(9:10) = 0.001;
%old = sdtdef('diag');

% Start the newton
if sp_util('issdt'); ofact('spfmex');ofact('silent');
elseif exist('pardiso','file');ofact pardiso; profile on
end
%sdtdef('diag',12);

dc=struct('def',zeros(length(model.DOF),2),'DOF',model.DOF);
mass=fe_mknl('assemble',model,Case,dc,2);

