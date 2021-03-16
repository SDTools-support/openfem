

% This TEST actually runs with SDT although the capabilities being tested
% are strictly OpenFEM
%       Copyright (c) 2002-2005 by SDTools, All Rights Reserved.
%       $Revision: 1.17 $  $Date: 2015/03/03 16:21:49 $


if ~sp_util('issdt');
 cd(fileparts(fileparts(which('fe_mk'))));eval('sdtcheck path');
 ofact('spfmex');
end


wd='c:/tmp/sdt/bucher/gyrotests';
if exist(wd,'dir') % Test given by Olivier Nicolas (EDF-R&D)

cd(wd);
model=ufread('blade_3d_tot.unv');
model.pl=m_elastic('dbval 1',[1 fe_mat('m_elastic','SI',1) 7.1e10 .3 7820]);
model.Elt=feutil('orient',model);
model.Elt(1,1:8)=[Inf abs('hexa8b') 0];
r1=stack_get(model,'','BAS_NO','getdata');
model=fe_case(model,'fixdof','BAS_NO',r1{4});

r1=[      1      8.85636E+00       2.30056E-09 
         2      8.85637E+00       2.19752E-09 
         3      8.88810E+00       2.27680E-09 
         4      9.14240E+00       2.21560E-09 
         5      9.14241E+00       1.97715E-09];
r2=[    1      3.13901E+01       2.02885E-10 
         2      3.19127E+01       1.88720E-10 
         3      3.19127E+01       1.98735E-10 
         4      3.34147E+01       1.95090E-10 
         5      3.34147E+01       1.63615E-10 ];

[m,k,model.DOF]=fe_mknl(model);
def=fe_eig({m,k,model.DOF},[6 5 0]);  norm(def.data(r1(:,1))-r1(:,2))

[c_g,k_g,k_e]=fe_cyclic('assemble',model,[0 0 157]);

def2=fe_eig({m,k_e,model.DOF},[5 5 0]);  
norm(def2.data(r2(:,1))-r2(:,2))
[def2.data(r2(:,1)) r2(:,2)]

end % Gyro tests
wd='c:/tmp/sdt/bucher/gyrotests';
if exist(wd,'dir') % Test given by Hadar Raz at technion
 cd(wd);
 Up=samcef('read cyl2.bdf');

 [m,k]=upcom(Up,'assemble');
 def=fe_eig({m,k,Up.DOF},[5 20 1e3]);

 model=feutil('rmfield',Up,{'file','wd','mind','Stack','copt','eltid','Opt'});
 model.Elt(1,1:9)=[Inf abs('hexa20b') 0];
 def1=fe_eig(model,[5 20 1e3]);
 model.Elt(1,1:9)=[Inf abs('hexa20') 0 0];
 def3=fe_eig(model,[5 20 1e3]);

 [def.data def1.data def3.data]
 copyfile(fullfile(fileparts(which('naswrite')),'dmap','mode.dat'),'.')


end % Gyro tests

% Simple test of a rotating beam

feplot;
model=femesh('teststruct hexa8b divide 30 1 1');
model.Node(:,5:7)=model.Node(:,5:7)*diag([10 .1 .1/.4]);
feplot(model)

model.Elt(1,1:8)=[Inf abs('hexa8b') 0];
model=fe_case(model,'reset');
[Case,model.DOF]=fe_mknl('init',model);

k=fe_mknl('assemble',model,Case,1);k0=k;
d1=struct('def',model.DOF*0,'DOF',model.DOF);
k=fe_mknl('assemble',model,Case,d1,5);
[norm(k-k0,'inf') k0(1) k(1)]


m=fe_mknl('assemble',model,Case,2);
range=linspace(0,30,20); % rotation speed in rad/s
fr=[];
for j1=1:length(range)
 [c_g,k_g,k_e]=fe_cyclic('assemble NoT',model,[0 0 range(j1)]);
 k_g=feutil('tkt',Case.T,k_g);
 k_e=feutil('tkt',Case.T,k_e);
 if j1==1; fr(:,j1)=fe_eig(m,k_g+k_e,[5 10 1e3]);
 else; fr(:,j1)=fe_eig(m,k_g+k_e,[5 10 0]);
 end
end
% Plot the campbell diagram
figure(1);plot(range,fr);set(gca,'ylim',[0 max(fr(:,2))])


% computation of gyroscopic matrices needs to update
% the case information and associated integration
% rule

% r1 <- Omega
omega=[1 2 3]; 
r1=[0 3 5;6 0 1;2 4 0];r2=[omega(:);-omega(:)];
r1(find(r1))=r2(r1(find(r1)))

EltConst=Case.GroupInfo{8};
EltConst.ConstitTopology{6}=reshape(1:9,3,3);
EltConst.StrainDefinition{6}=EltConst.StrainDefinition{2};
EltConst.StrainLabels{6}=EltConst.StrainLabels{2};
EltConst=integrules('matrixrule',EltConst);
C2=Case;C2.GroupInfo{8}=EltConst;
C2.GroupInfo{2}(7,:)=0; % no offset in Constit
C2.GroupInfo{4}=r1(:); c_g=fe_mknl('assemble',model,C2,6); %
C2.GroupInfo{4}=r1^2 ; k_g=fe_mknl('assemble',model,C2,6); %

% In SDT use
%[c_g,k_g]=fe_cyclic('assemble',model,[1 2 3]);

% ----------------------------------------------------------------------------------
if 1==2% NASTRAN COMPARISONS

cd c:/tmp/sdt/bucher/GyroTests
edits={'set','PARAM','POST','-2'};
nas2up('editbulk','mode_rotating.dat',edits,'nas_rot.dat');
nas2up('joball','nas_rot.dat')

end
