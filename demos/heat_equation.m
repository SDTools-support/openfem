
% The theoretical part of this example can be found in the OpenFEM manual
% "Application examples : Heat equation"

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Mesh the model - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

model=ofdemos('heatstatic');

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Build the loading conditions on nodes of the interface
% This contains g+\alpha \theta_ext which will be integrated on the surface

% **** data ****
RO=struct('TOnFace1',25, ...
          'Alpha',1, ...
          'TExt',20, ...
          'F',40);

% **************

% boundary conditions on face x==0
data1 = struct('sel','group2',...
            'eltsel','group2',...
            'def',RO.TOnFace1,'DOF',.20);
model = fe_case(model,'Fvol','Surface_load1',data1);
        % the surface integrals correspond to volume integrals on
        % the surface elements representing the boundary

% boundary conditions on face x==10
mdof0=feutil('findnode group3',model)+.20;
data2 = struct('DOF',mdof0,'def',ones(size(mdof0))*RO.TExt);
%data2=struct('DOF',.20,'def',RO.TExt); % xxx eb should test formal
model = fe_case(model,'DofSet','Surface_load2',data2);

% boundary conditions on other faces
node=feutil('getnode group4',model);
data3 = struct('sel','group4','def',25 - node(:,5).^2/20,'DOF',node(:,1)+.20);
data3=struct('sel','group4','def','25-x.^2/20','DOF',.20);
% workr
model = fe_case(model,'Fvol','Surface_load3',data3);

% constant internal heat source
data4 = struct('sel','group1','def',[RO.F],'DOF',.20);
model = fe_case(model,'Fvol','Internal source',data4);

% Assemble matrices and RHS - - - - - - - - - - - - - - - - - - - - - - - - -

[Case,model.DOF]=fe_mknl('init -gstate',model);

% try; if RunOpt.debug
%  sdtdef('diag',12);  % very slow elem0.m file implementation, just check
%  % the right-hand side is in dc.def(:,2) which is modified even though it is
%  % now returned by fe_mknl
%  dc2=struct('def',zeros(length(model.DOF),2),'DOF',model.DOF);
%  k=fe_mknl('assemble not',model,Case,dc2,1);
%  sdtdef('diag',0); % modified here
% end;end

dc=struct('def',zeros(length(model.DOF),2),'DOF',model.DOF);
k=fe_mknl('assemble not',model,Case,dc,1);
if norm(dc.def(:,2))==0; 
    warning('RHS should be computed with k'); 
    Load=fe_load(model,'NoT');
    dc.def(:,2)=fe_c(dc.DOF,Load.DOF,sum(Load.def,2)')';
    %dc.def(:,3)=fe_c(dc.DOF,Load.DOF,sum(Load.def,2)')';
end

% computation of the response to internal forces and enforced displacement
ind=fe_c(model.DOF,Case.DofIn,'ind');cind=setdiff(1:length(model.DOF),ind);
resp=struct('def',zeros(size(model.DOF)),'DOF',model.DOF);
resp.def(ind,1)=Case.TIn(ind,:);
F=dc.def(cind,2)-k(cind,ind)*Case.TIn(ind,:);
resp.def(cind,1)=ofact(k(cind,cind),F);

% View results  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

feplot(model,resp); try;fecom('colordata20');end

[i1,n1]=feutil('findnode z==2 & y==2',model);
figure(10);subplot(311); % (figure 4.1 of the manual)
plot(n1(:,5),fe_c(resp.DOF,i1)*resp.def);
%title('Temperature on the y=z=0 edge');
xlabel('x coordinate'); ylabel('temperature')

% This was used to generate the manual figure
% set(10,'position',[520 911 560 199]);
% convert('-pautops2pdf','../tex/plots/heat_eq.eps')

[i1,n1]=feutil('findnode x==5 & y==2',model);
subplot(312);plot(n1(:,7),fe_c(resp.DOF,i1)*resp.def)
title('Temperature at points x=5, y=2');

[i1,n1]=feutil('findnode x==5 & z==2',model);
subplot(313);plot(n1(:,6),fe_c(resp.DOF,i1)*resp.def)
title('Temperature at points x=5, z=2');

feutil('infoelt',model)


%	Etienne Balmes, Frederic Bourquin, A. Nassiopoulos
%       Copyright (c) 2001-2008 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       $Revision: 1.17 $  $Date: 2020/02/26 08:47:39 $
