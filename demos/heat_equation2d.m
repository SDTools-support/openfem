
% The theoretical part of this example can be found in the OpenFEM manual
% "Application examples : Heat equation"

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Mesh the model - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

femesh reset
femesh('testquad4 divide 4 10');
femesh(';addsel;selelt seledge;addsel')
femesh(';selgroup1;quad2tria;addsel;selgroup 3 2');
FEelt=FEel0;
femesh('dividegroup2 innode {x==0}')

femesh('set group1 mat 1 pro 1')
femesh('set group2 mat 1 pro 2 name bar1')
femesh('set group3 mat 1 pro 2 name bar1')
model=femesh('model');model.Node=model.Node*diag([1 1 1 1  10 1 1]);
                                   
model.pl=[1 fe_mat('m_heat','SI',1) 400 1 1 1   % k rho C alpha  
          2 fe_mat('m_heat','SI',1) 400 1 1 1e10  % k rho C (alpha_0: not used)
          ];
model.il=[1 fe_mat('p_heat','SI',1) 0 -3 2;   % volume integration
          2 fe_mat('p_heat','SI',2) 0 -3 2];  % exchange surface integration


RO=struct('TOnFace1',25, ...
          'Alpha',1, ...
          'TExt',20, ...
          'F',40);

% Face x==0 : enforced flux. 
data1 = struct('sel','eltname bar1 & innode {x==0}',...
            'def',RO.TOnFace1,'DOF',.20);
model = fe_case(model,'Fvol','Surface_load1',data1);
        % the surface integrals correspond to volume integrals on
        % the surface elements representing the boundary

% Face x==10 : Dirichlet (enforced temperature), no need for alpha
mdof0=feutil('findnode x==10',model)+.20;
model.Elt=feutil('removeelt eltname bar1 & innode x==10',model); 
data2 = struct('DOF',mdof0,'def',ones(size(mdof0))*RO.TExt);
 % data2=struct('DOF',.20,'def',RO.TExt); % xxx eb should work
model = fe_case(model,'DofSet','Surface_load2',data2);

% Top and bottom face : given flux
node=feutil('getnode group4',model);
%data3 = struct('sel','group4','def',25 - node(:,5).^2/20,'DOF',node(:,1)+.20);
data3=struct('sel','eltname bar1 & innode {y==0 | y==1}','def','25-x.^2/20','DOF',.20);% xxx eb should
% workr
model = fe_case(model,'Fvol','Surface_load3',data3);

% constant internal heat source
data4 = struct('sel','group1','def',[RO.F],'DOF',.20);
model = fe_case(model,'Fvol','Internal source',data4);

% Assemble matrices and RHS - - - - - - - - - - - - - - - - - - - - - - - - -

DOF=feutil('getdof',model);[Case,model.DOF]=fe_mknl('init -gstate',model);

 if exist('RunOpt','var')&&isfield(RunOpt,'debug')&& RunOpt.debug
  sdtdef('diag',12);  % very slow elem0.m file implementation, just check
  % the right-hand side is in dc.def(:,2) which is modified even though it is
  % now returned by fe_mknl
  dc2=struct('def',zeros(length(model.DOF),2),'DOF',model.DOF);
  [k,C2,dc2]=fe_mknl('assemble not',model,Case,dc2,1);
  sdtdef('diag',0); % modified here
 end

dc=struct('def',zeros(length(model.DOF),2),'DOF',model.DOF);
[k,Case,dc]=fe_mknl('assemble not',model,Case,dc,1);
if ~any(dc.def(:,2)); error('Problem with assembly');end

% computation of the response to internal forces and enforced displacement
ind=fe_c(model.DOF,Case.DofIn,'ind');cind=setdiff(1:length(model.DOF),ind);
resp=struct('def',zeros(size(model.DOF)),'DOF',model.DOF);
resp.def(ind,1)=Case.TIn(ind,:);
F=dc.def(cind,2)-k(cind,ind)*Case.TIn(ind,:);
resp.def(cind,1)=ofact(k(cind,cind),F);

% View results  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

feplot(model,resp); try;fecom('colordata20');end

n1=feutil('getnode y==1',model);[r1,i2]=sortrows(n1(:,5:7));n1=n1(i2,:);
figure(10); % (figure 4.1 of the manual)
plot(n1(:,5),fe_c(resp.DOF,n1(:,1))*resp.def);
%title('Temperature on the y=z=0 edge');
xlabel('x coordinate'); ylabel('temperature')


%	Etienne Balmes, Frederic Bourquin, A. Nassiopoulos
%       Copyright (c) 2001-2008 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       $Revision: 1.8 $  $Date: 2020/02/26 08:47:43 $
