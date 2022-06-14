%-------------------------------------------------------------
if ~exist('RO','var'); RO=struct;end
% Element ElemF you can select : tetra4b, tetra10b, hexa8b, hexa20b....
if ~isfield(RO,'ElemF');RO.ElemF='hexa20';end; 
% div fix number of subdivisions in the mesh. Increase to refine mesh
if ~isfield(RO,'div'); RO.div=3; end
if ~isfield(RO,'nrand');RO.nrand=0; end % Use uniform mesh
RO.C1=.3; RO.C2=.2; RO.K=.3; %C1=0.3MPa, C2=0.2MPa, K=0.3MPa

 
% Simple cube with displacement field u1 = l1x1, u2 = l2x2, u3 = l3x3.
% you can change the L values
if ~isfield(RO,'L');RO.L=[.1 .2 .3];end % Lambda values (see theory)


% ------------------------------------------------------------
% Select another solver : UMFPACK,PARDISO  activate only if provided
% ------------------------------------------------------------
%ofact pardiso


%--------------------------------------------------------------------------
%---------------------------------------------------------------------------
% FIRST PASS: COMPUTATION WITH EXTERNAL PRESURE LOAD
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

disp('COMPUTATION WITH EXTERNAL PRESSURE LOAD')


%---------------------------------------------------------------------------
% COMPUTATION OF SIGMA 2ND PIOLA-KIRCHHOFF TENSOR
%--------------------------------------------------------------------------

r1=0;for ji=1:3;for jj=1:ji-1; r1=r1+(1+RO.L(ji))^2*(1+RO.L(jj))^2;end;end

I=[sum((1+RO.L).^2) r1 prod((1+RO.L).^2)]; %INVARIANT COMPUTATION
[dWdI,d2WdI2]=feval(m_hyper('@EnHyper'),[],[0 RO.C1 RO.C2 RO.K 0]',I); %ENERGY DERIVATIVES

for ji=1:3 % compute Signam
 r1(ji)=2*(dWdI(1)+dWdI(2)*(I(1)-(1+RO.L(ji))^2) + dWdI(3)*I(3)*(1+RO.L(ji))^-2);
end
Sigma=diag(r1);
feval(m_hyper('@checkNL'),struct('opt',[0 0 0 0 RO.C1 RO.C2 0 RO.K]','unl',reshape(diag(RO.L),[],1)))

%-------------------------------------------------------------------------
% DEFINE THE STRUCTURE MODEL: .Node, .Elt, .pl, .il
%------------------------------------------------------------------------

% Mesh
model=femesh(sprintf('testStruct %s %s',RO.ElemF, ...
  sprintf('divide %i %i %i rand %.15g back',RO.div*[1 1 1],RO.nrand)));
model.Node(:,7)=model.Node(:,7)/max(abs(model.Node(:,7)));

% Boundary conditions
model=fe_case(model,'reset','FixDof','x','x==0 -dof 1', ...
   'FixDof','y','y==0 -dof 2','FixDof','z','z==0 -dof 3');

%Material property set all elements to matid 100 and proid 111
model.Elt=feutil('set group1 matid 100 proid 111',model);
model.pl=m_hyper('dbval 100 Ref'); % this is where the material is defined
model.il=p_solid('dbval 111 d3 3');

%------------------------------------------------------------------------
% COMPUTE EXTERNAL PRESSURE LOAD (see first pass in doc)
%-----------------------------------------------------------------------

data=struct('sel','x==1','def',-(1+RO.L(1))*Sigma(1,1),'DOF',.19);
model=fe_case(model,'FSurf','xs',data);
data=struct('sel','y==1','def',-(1+RO.L(2))*Sigma(2,2),'DOF',.19);
model=fe_case(model,'FSurf','ys',data);
data=struct('sel','z==1','def',-(1+RO.L(3))*Sigma(3,3),'DOF',.19);
model=fe_case(model,'FSurf','zs',data);
Load=fe_load(model); F1=sum(Load.def,2);



%DISPLAY THE MODEL
feplot(model);


%---------------------------------------------------------------------------
% NEWTON RESOLUTION
%--------------------------------------------------------------------------

%INITALISATION
[Case,model.DOF]=fe_mknl('init',model);
dc=struct('def',zeros(length(model.DOF),3),'DOF',model.DOF);
ind=fe_c(model.DOF,Case.DOF,'ind');check=[0 Inf];
maxnewton = 10;

% Start the newton
profile on

q=[];dq=[];
%DETERMINE NUMBER OF INCREMENT
% disp(comstr(feutil('rmfield',curve,'Edit','PlotInfo','Xlab'),-30))
curve=struct('X',{{[0;0.25;0.5;0.75;1]}}, ...
'Y',[0;0.25;0.5;0.75;1],'name','ramp','Extrap','flat');
if RO.nrand; curve=fe_curve('testramp 0 1 1',linspace(0,1,50));end

disp('Nbitertotal  rezidual  norm(dq)   norm(q) incrementstep')

iter=0;
for jStep=1:length(curve.Y)  % Outer loop on force increments
 F=F1*curve.Y(jStep);check(2)=Inf;
 for j1=1:maxnewton;  % Inner loop on convergence at current step
  dc.def(:,2)=0; Debug=0;
  % assemble tamgent matrix and rhs
  k=fe_mknl('assemble not',model,Case,dc,5);
  if ~any(dc.def(:,2))&&any(dc.def(:,1));error('Problem with RHS computation');end
  R=F-Case.T'*(dc.def(:,2)); %+k*dc.def(:,1)); 

  check=[iter norm(R) norm(dq) norm(q) curve.Y(jStep) ]; 
  disp(check);%drawnow;
  iter=iter+1; 
  if check(2)<1e-6;break; %convergence test
  elseif j1 == maxnewton 
      sdtw('newton did not converge step:%d lastrezidual:%e',jStep,check(2));
      check(1)=Inf;break;
  end
  ocheck=check;

  %resolution by lu factorization
  dq=ofact(k(ind,ind),R);
  q=dc.def(ind,1);

  %solution update
  dc.def(ind,1)=q+dq; %/max([norm(dq),norm(q),1])*.5;
 end % inner loop
 if ~isfinite(check(1)); break;end
end  % outer loop
try;profile report;end

%BUILD THEORETICAL DISPLACEMENT
def=struct('def',zeros(size(model.DOF)),'DOF',model.DOF);
NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
i1=fe_c(model.DOF,.01,'ind');
def.def(i1,1)=model.Node(NNode(fix(model.DOF(i1))),5)*RO.L(1);
i1=fe_c(model.DOF,.02,'ind');
def.def(i1,1)=model.Node(NNode(fix(model.DOF(i1))),6)*RO.L(2);
i1=fe_c(model.DOF,.03,'ind');
def.def(i1,1)=model.Node(NNode(fix(model.DOF(i1))),7)*RO.L(3);

%COMPARE THEORETICAL AND COMPUTED DISPLACEMENTS FOR FIRST PASS
err1=norm(dc.def(:,1)-def.def(:,1));
fprintf('\nDifference with ref. solution for external pressure load: %g\n',err1)


%//////////////////////////////////////////////////////////////////////////////////
disp('---------------------------------------------------------')
disp('---------------------------------------------------------')
disp(' PRESS A KEY TO RUN 2ND PASS ')
pause
%--------------------------------------------------------------------------
%---------------------------------------------------------------------------
% 2ND PASS: COMPUTATION WITH FOLLOWER PRESSURE 
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

disp('COMPUTATION WITH FOLLOWER PRESSURE')

%--------------------------------------------------------------------
% Simple cube with displacement field u1 = l1x1, u2 = l2x2, u3 = l3x3,
% l1 = l2 = l3 (see theory)
L=RO.L(1)*ones(1,3);

%---------------------------------------------------------------------------
% COMPUTATION OF SIGMA 2ND PIOLA-KIRCHHOFF TENSOR
%--------------------------------------------------------------------------

%
r1=0;for ji=1:3;for jj=1:ji-1; r1=r1+(1+L(ji))^2*(1+L(jj))^2;end;end

%INVARIANTS COMPUTATION
I=[sum((1+L).^2) r1 prod((1+L).^2)];

%ENERGY DERIVATIVES
[dWdI,d2WdI2]=feval(m_hyper('@EnHyper'),[],[0 RO.C1 RO.C2 RO.K 0]',I); 
%SIGMA
for ji=1:3;     
 r1(ji)=2*(dWdI(1)+dWdI(2)*(I(1)-(1+L(ji))^2) + dWdI(3)*I(3)*(1+L(ji))^-2);
end
Sigma=diag(r1);


%-------------------------------------------------------------------------
% DEFINE THE STRUCTURE MODEL: .Node, .Elt, .pl, .il
%------------------------------------------------------------------------

% Mesh
femesh('reset');
model=femesh(horzcat('teststruct',RO.ElemF, ...
  sprintf('divide %i %i %i rand %.15g back',RO.div*[1 1 1],RO.nrand)));
model.Node(:,7)=model.Node(:,7)/max(abs(model.Node(:,7)));

% Boundary conditions
model=fe_case(model,'reset','FixDof','x','x==0 -dof 1', ...
   'FixDof','y','y==0 -dof 2','FixDof','z','z==0 -dof 3');

%Material property
model.Elt(2:end,length(feval(RO.ElemF,'node'))+1)=100; % all elements with MatId 100
model.pl=m_hyper('dbval 100 Ref'); % this is where the material is defined
model.il=p_solid('dbval 111 d3 3');

%------------------------------------------------------------------------
% FOLOWER PRESSURE (see 2nd pass in doc RivlinCube)
%-----------------------------------------------------------------------

%ACTIVATE FOLLOWER PRESSURE PROPERTY
 elt=feutil('selelt selface & innode {x==1 | y==1 | z==1}',model);
 elt=feutil('egid-vol',elt);
 elt=feutil('set group1 proid 200',model.Node,elt);
 model.Elt(end+(1:size(elt,1)),1:size(elt,2))=elt;
 model.il=p_solid(model.il,'dbval 200',[200 fe_mat('p_solid','SI',3) 0 4 0]);

%COMPUTE PRESSURE see doc...
 p1=-(1+L(1))/(1+L(2))/(1+L(3))*Sigma(1,1);
 p2=-(1+L(2))/(1+L(1))/(1+L(3))*Sigma(2,2);
 p3=-(1+L(3))/(1+L(1))/(1+L(2))*Sigma(3,3);


%DISPLAY THE MODEL
feplot(model);

%---------------------------------------------------------------------------
% NEWTON RESOLUTION
%--------------------------------------------------------------------------

%INITALISATION
[Case,model.DOF]=fe_mknl('init',model);
dc=struct('def',zeros(length(model.DOF),3),'DOF',model.DOF);
 indn=femesh('findnode x==1');
 indp1=fe_c(model.DOF,indn+.19,'ind');
 indn=femesh('findnode y==1');
 indp2=fe_c(model.DOF,indn+.19,'ind');
 indn=femesh('findnode z==1');
 indp3=fe_c(model.DOF,indn+.19,'ind');
 io=fe_c(Case.DOF,[.01;.02;.03],'ind');
 ind=fe_c(model.DOF,Case.DOF,'ind');check=[0 Inf];
 maxnewton = 10;

% Start the newton
profile on

q=[];dq=[];

%DETERMINE NUMBER OF INCREMENT
i1=5;if RO.nrand; i1=50;end
curve=struct('X',{{(0:i1-1)'}},'Y',linspace(0,1,i1)');


disp('Nbitertotal  rezidual  norm(dq)   norm(q) incrementstep')

iter=0;
for jStep=1:length(curve.Y)  % Outer loop on force steps

 dc.def(indp1,1)=p1*curve.Y(jStep);
 dc.def(indp2,1)=p2*curve.Y(jStep);
 dc.def(indp3,1)=p3*curve.Y(jStep);
 check(2)=Inf;
 for j1=1:maxnewton;  % Inner loop on convergence at current step
   dc.def(:,2)=0; Debug=0;
   k=fe_mknl('assemble not',model,Case,dc,5);
   R=-Case.T(:,fe_c(Case.DOF,[.01;.02;.03],'ind'))'*(dc.def(:,2));

   check=[iter norm(R) norm(dq) norm(q) curve.Y(jStep) ]; 
   disp(check);%drawnow;
   iter=iter+1; 
   if check(2)<1e-6;break; %convergence test
   elseif j1 == maxnewton 
       sdtw('newton did not converge step:%d lastrezidual:%e',jStep,check(2));
       check(1)=Inf;break;
   end
   ocheck=check;

   %resolution by lu factorization
   dq=ofact(k(ind(io),ind(io)),R);
   q=dc.def(ind(io),1);

   %solution update
   dc.def(ind(io),1)=q+dq; %/max([norm(dq),norm(q),1])*.5;
 end % inner loop
 if ~isfinite(check(1)); break;end
end  % outer loop
try;profile report;end

%BUILD THEORETICAL DISPLACEMENT
def=struct('def',zeros(size(model.DOF)),'DOF',model.DOF);
NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
i1=fe_c(model.DOF,.01,'ind');
def.def(i1,1)=model.Node(NNode(fix(model.DOF(i1))),5)*L(1);
i1=fe_c(model.DOF,.02,'ind');
def.def(i1,1)=model.Node(NNode(fix(model.DOF(i1))),6)*L(2);
i1=fe_c(model.DOF,.03,'ind');
def.def(i1,1)=model.Node(NNode(fix(model.DOF(i1))),7)*L(3);

%COMPARE THEORETICAL AND COMPUTED DISPLACEMENTS FOR 2ND PASS
icomp=fe_c(model.DOF,[.01;.02;.03],'ind');
err2=norm(dc.def(icomp,1)-def.def(icomp,1));
fprintf('\nDifference with ref. solution for follower pressure load: %g\n',err2)


disp('---------------------------------------------------------')

%disp(sprintf('EXTERNAL PRESSURE LOAD: error = %e',err1))
%disp(sprintf('FOLLOWER PRESSURE:      error = %e',err2))
