function [out,out1]=t_xfem(varargin);

% Test Level set and XFEM developments
%
%
% Contributed by Eric Monteiro and Etienne Balmes : ENSAM / PIMM

%       Copyright (c) 2001-2023 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       For revision information use feutil('cvs')

%#ok<*ASGLU,*NASGU,*STOUT,*TRYNC,*NOSEM>

if nargin==0
  %% Tests with no output artuments
  if ~isempty(strfind(which('fe_mk'),'sdt.cur'))
    cd(fileparts(fileparts(which('fe_mk'))));eval('sdtcheck path');%Test with SDT
    clear variables global;comgui('close all');cinguj('initSwing');
  end
  t_xfem('piezomesh')
  t_xfem('volumetest')
  t_xfem('piezo2')
  %t_xfem('bug') % This should not be executed in daily tests
  %t_xfem('cutelem'); % This should not be executed in daily tests
  t_xfem('remap') % Actually another test
  t_xfem('triafc'); % face2tria with underlying volume transform
  t_xfem('volcut;')
  return
  
end
[CAM,Cam]=comstr(varargin{1},1);carg=2;

%% #Piezomesh ---------------------------
if comstr(Cam,'piezomesh')
  
  model=femesh('testquad4 divide 50 50');
  def=[];li={};
  for j1=1:11
    li{j1}=struct('shape','rect','xc',.5,'yc',.5,'lx',.1,'ly',.2,'alpha',30*(j1-1));
  end
  def=lsutil('gen',model,li);
  lsutil('viewLS',model,def)  
  
  %idef=6;model.Node(:,7)=def.def(:,idef);feplot(model,def);fecom(';colordata98;colorbar');set(gcf,'renderer', 'zbuffer');
  
  %case SHM
  % sdtweb _tracker pje 1337
  r1=d_shm('MeshCfg',struct,'PimmShellA');r1=r1.param.MeshCfg.data{1};
  r2=r1;r2.list(3:end,:)=[];mo2=d_piezo('MeshPlate',r2);
  
  li=r1.list(3:end,3);def=lsutil('gen',mo2,li);
  cf=feplot(mo2);cf.def=struct('def',max(def.def,[],2),'DOF',def.DOF);
  fecom(';colordata98;colorbar;colorscale 2sided');
  if sdtdef('verm')<804;set(gcf,'renderer', 'zbuffer');end
  fecom('colormap','jet(5)')
  
  % Mesh connectivity line
  
  lsutil('EdgeCut',mo2,def,struct('method',@(x)max(x,[],2),'morph',.5))
  conn=feval(feutilb('@levNodeCon'),[],mo2,'econ');
  
  
elseif comstr(Cam,'piezo2')
  %% #Piezo2 -------------------------------------------------------------------
  RO=struct('dim',[400 300 8],'tolE',.3);
  [model,li]=ofdemos('LS2d',RO);lsutil('ViewLs',model,li);
  mo3=model; for j1=1:length(li); mo3=lsutil('cut',mo3,li(j1),RO);end;
  lsutil('ViewLs',mo3,li);
  %def=fe_eig(mo3);
  
elseif comstr(Cam,'pad')
 %% #Pad: pad on disc case (functional stability)
 model=d_contact('tutoSimpleBrakeLS -s3 -reset mo3')

elseif comstr(Cam,'d2Tri')
  %% #d2Tri : distance to triangulated surface

 model=femesh('testhexa8 -divide 10 10 10');
 model=feutil('hexa2tetra',model);
 cf=feplot(model); 
 
 mo2=model;mo2.Elt=feutil('selelt innode{x>=.5} & selface & innode {x==.5}',mo2);
 mo2.shape='Tri';li={mo2};
 lsutil('viewls',cf,li);fecom('colorscaleinstant')

 mo2=model;mo2.Elt=feutil('selelt innode{x>=.5&y>=.5&z>=.1 & z<=.9} & selface',mo2);
 mo2.shape='Tri';li={mo2};
 cf.sel={'innode {z<=.501}'};
 lsutil('viewls',cf,li);fecom('colorscaleinstant')
 
  
elseif comstr(Cam,'volcut')
  %% #VolCut -------------------------------------------------------------------
  RO=struct('dim',[10 10 40],'tolE',.1);
  if ~exist('Cam','var');Cam=' ';end
  [model,li]=ofdemos('LS3d',RO);if Cam(end)~=';';lsutil('ViewLs',model,li);end
  
  %mo2=lsutil('split',model,li,RO);if Cam(end)~=';';lsutil('ViewLs',mo2,li);end
  %RO.mpid=[2 2];
  mo3=lsutil('cut',model,li,RO);if Cam(end)~=';';lsutil('ViewLs',mo3,li);end
  [EGroup,nGroup]=getegroup(mo3.Elt);%comstr(EGroup,-30)
  if ~isequal(EGroup,[1,3956,4003,4160,4213,4262,4287]);
    lsutil('ViewLs',mo3,li)
    error('Change of behaviour');
  end
  feplot('ShowFiMat')
  
  %% #VolCut_cylinder
  li={struct('shape','cyl','xc',.5,'yc',.5,'zc',1,'nx',0,'ny',0,'nz',-1, ...
    'rc',.2,'z0',-.4,'z1',.4,'mpid',[200 300])};
  mo3=lsutil('cut',model,li,RO);if Cam(end)~=';';lsutil('ViewLs',mo3,li);end
  cf=feplot;cf.sel={'innode {x>=.5}','colordatamat -edgealpha1'}
  
  %fecom showline;fecom('showmap','EltOrient');
  
  %% #VolCut_isosurface in femesh
 lsutil('ViewLs',model,li); 
 feval(lsutil('@iso_sel'),'init',struct('cf',cf)); 
 feval(lsutil('@iso_sel'),'def',struct('cf',cf,'def',cf.def)); 
 if 1==2

    sel=feval(lsutil('@iso_sel'),'init',cf);sel.off=0;sel.step=.005; 
    feval(lsutil('@isoSurf'),sel,cf.def,cf.ga);
    %   
    st=feval(lsutil('@charLs'),cf.def,elt);
    st1=unique(st(~map.isKey(st))); for j1=1:length(st1);map(st1{j1})=[];end
    cellfun(@(x)map(x),st,'uni',0)  
    % compute dist and charLS of those elements
    
 end

 model=femesh('testhexa8');model=feutil('hexa2tetra',model);model=feutil('lin2quad',model);
 model=feutil('renumber',model,model.Node(:,1)*2);
 
 li={struct('shape','toplane','xc',0,'yc',0,'zc',.3,'nx',0,'ny',0,'nz',1,'mpid',[2 2])};
 mo4=lsutil('cut',model,li,struct('Etol',.1));
 feplot(mo4); fecom colordatapro-alpha.4;
 fecom('shownodemark',feutil('findnode groupall',mo4))

 
elseif comstr(Cam,'cut');[CAM,Cam]=comstr(CAM,4);
  %% #Cut -------------------------------------------------------------------
  RO=struct('dim',[10 10 40],'tolE',.5);
  if ~exist('Cam','var');Cam=' ';end
  [model,li]=ofdemos('LS3d',RO);if Cam(end)~=';';lsutil('ViewLs',model,li);end
  mo2=lsutil('cut',model,li,RO);if Cam(end)~=';';lsutil('ViewLs',mo2,li);end
  cf=feplot;cf.sel={'innode {y<=.2}','colordatamat -edgealpha1'}
  keyboard
  
  
elseif comstr(Cam,'bug');[CAM,Cam]=comstr(CAM,4);
  %% #BugA : problems with DampBench ------------------------------------------
  %   
  if comstr(Cam,'corner')
    %t_xfem('bugcorner');
    RO=struct('dim',[10,10,10],'tolE',0.2,'ls',[],'vf',0.7);
    RO.ls={struct('shape','box','xc',.5,'yc',.5,'zc',0.,'nx',[1. 0 0],'ny',[0 1. 0],'nz',[0 0 1.], ...
      'lx',sqrt(RO.vf)/2,'ly',sqrt(RO.vf)/2,'lz',1.4,'mpid',[200 300])};
    model=ofdemos('LS3d',RO);model=lsutil('cut',model,RO.ls,RO);lsutil('ViewLs',model,RO.ls);
  elseif comstr(Cam,'mat')
    %t_xfem('bugmat');
    RO=struct('dim',[20,20,20],'tolE',0.1,'ls',[],'vf',0.5,...
        'onSurf','findnode z==1 | z==0|x==0 | y==0 | x==1 | y==1');
    RO.ls={struct('shape','cyl','xc',.5,'yc',.5,'zc',0.,'nx',0,'ny',0,'nz',1, ...
        'rc',sqrt(RO.vf/pi),'z0',-.4,'z1',1.4,'mpid',[200 300])};
    model=ofdemos('LS3d',RO);model=lsutil('cut',model,RO.ls,RO);  
    feplot(model);fecom('colordata mat')
  elseif comstr(Cam,'penta')
    %t_xfem('bugpenta');
    RO=struct('dim',7*[1,1,1],'tolE',.1); [model,li]=ofdemos('LS3d',RO);
    model.pl=m_elastic('dbval 100 steel','dbval 200 aluminum');
    model.il=p_solid('dbval 111 d3 -3','dbval 300 d3 -3');
    model=feutil('hexa2penta',model);% Does not propagate eltorient as expected
    lsutil('ViewLs',model,li); 
    %model=stack_rm(model,'info','EltOrient');
    mo3=lsutil('cut',model,li,RO);lsutil('ViewLs',mo3,li);
    if sdtdef('isinteractive');fecom('showline');fecom('showmap','EltOrient');end
  elseif comstr(Cam,'quality')
    %t_xfem('bugquality');
    RO=struct('dim',[3,3,3],'tolE',0.2,'ls',[],'vf',0.5,'Fixed','z==0 | z==1');%'onSurf','findnode z==0 | z==1'
    RO.ls={struct('shape','cyl','xc',0.5,'yc',0.5,'zc',0.,'nx',0,'ny',0,'nz',1, ...
      'rc',sqrt(RO.vf/pi),'z0',-1.,'z1',2.,'mpid',[200 300])};
    
    model=ofdemos('LS3d',RO);
    model.pl=m_elastic('dbval 100 steel','dbval 200 aluminum');
    model.il=p_solid('dbval 111 d3 -3','dbval 300 d3 -3');
    model=lsutil('cut',model,RO.ls,RO);
    cf=feplot(model);cf.model=model;fe_quality(cf.mdl);fe_quality('view',cf);
    %fe_quality('view',cf,'EltName tetra4')
    %Elt=feutil('Orient',model);ind=find(sum(model.Elt~=Elt,2));if~isempty(ind);keyboard;end
    keyboard    
    
  elseif comstr(Cam,'failedcoal')
   f1=sdtcheck('PatchFile',struct('fname','pad.igs','in','demo_squeal_abq.zip'));
   %fname=comgui('cdO:\distrib\contrib\demo_squeal_abq\pad.igs');
   RO=struct( ... % Predefine materials
    'pl',m_elastic('dbval -unit TM 1 steel'), ...
    'sel','selelt seledge', ... % Elements to retain at end
    'Run','-1 -order 1 -clmax 3 -clmin 2 -v 0');  %RunCommand
   mo1=fe_gmsh('write',f1,RO);
   % Volume meshing of disk
   model=d_contact('discbrake -nopad -noplot');
   model.Node(:,5:7)=model.Node(:,5:7)*1000;% .igs MM/model is SI
   % model=feutil('objectdivide 5 1 2',model,'eltname hexa');
   model=feutil('objectdivide 2 2 1',model,'eltname hexa');
   
   RC=struct('tolE',.3,'Fixed','inelt{seledge}',...
    'onSurf',sprintf('findnode z==%.15g',max(model.Node(:,7))),'normal',[0 0 1],...
    'keepOrigMPID',1,'NoDegen',1);
   [mo3,li]=lsutil('SurfRemesh',model,mo1,RC);lsutil('ViewLs',mo3,li);fecom('coloredgealpha1');
   
   % Problem: degen elements with duplicated nodes
   cf=feplot(mo3)
   cf.sel='withnode 18 & withnode 4423'
   elt=feutil('selelt withnode 18 & withnode 4423',mo3)
   n3=feutil('getnode 18 4423',mo3)

    
  else
    
    sdtw('_ewt1_23-Oct-2018','Need to cleanup')
    if ~sd('_user','balmes'); return; end
    
    wd='C:\Users\monteiroe\Dropbox\Collab\SDT_trackers\';
    if ~exist(wd,'dir'); wd=sdeb('wd/sdtdata/collab/SDT_trackers');end
    RO=struct('tolE',.3);load(fullfile(wd,'LsBug.mat'),'mdl');
    mdl=feutil('rmfield',mdl,'K','Klab','TR');
    RO.Fixed=feutil('findnode inelt{seledge}',mdl); % Nodes that cannot be moved
    %RO.Usable=feutil('findnode inelt{selface}',mdl); % Superset of cutable edge nodes
    
    li={struct('shape','cyl','xc',-12.5,'yc',0,'zc',-30,'nx',0,'ny',-1,'nz',0, ...
      'rc',12.5,'z0',-30,'z1',30,'mpid',[200 300])};
    
    mo3=lsutil('cut',feutil('quad2lin epsl .01',mdl),li,RO);%lsutil('ViewLs',mo3,li);
    mo4=feutil('Lin2Quad epsl .01',mo3);lsutil('ViewLs',mo4,li);
    fecom('colordatamat -edgealpha.1');fecom('view1')
    
  end
  
elseif comstr(Cam,'remap')
  %% #ScriptRemap-------------------------------------------------------------------
  sdtkey('cvsnum > 1.059','m_piezo')
  
  RO=struct('ndiv',[50 10 10],'tolE',.1);
  [model,li]=t_xfem('Meshremap',RO); 
  mo3=model;for j1=1:length(li); mo3=lsutil('cut',mo3,li(j1),RO);end
  cf=feplot(mo3);fecom('colordatamat-edgealpha.1');
  p_piezo('ElectrodeDOF',mo3)
  
elseif comstr(Cam,'volumetest')
  %% #VolumeTest : LS display tests ---------------------------
  
  RO=struct('dim',[10 10 40],'tolE',.1);
  [model,li]=ofdemos('LS3d',RO); % li is sphere

  cf=feplot(model);
  lsutil('viewls',cf,li);
  lsutil('viewls',cf,{struct('shape','seg','orig',[.2 .2 0],'normal',[0 0 1],'z0',.1,'z1',.5')});
  lsutil('viewls',cf,{struct('shape','plane','orig',[.2 .2 .2],'normal',[0 1 1])});
  
  %def=lsutil('gen',model,li);
  %fecom('shownodemark',def.DOF(max(def.def,[],2)>0))
  if sdtdef('isinteractive')
    model=lsutil('eltset',model,struct('list',{li},'CritFcn',@(x)x>=0, ...
      'sel','innode'));
    cf=feplot(model);
    sdtroot([],struct('cf',cf.opt(1),'setname','_lsset'),'initEltSet-reset')
  end
  
  %% LevelLines selection
  cf=feplot(model);
  data=struct('dir',{{'x'}},'DOF',[.98]);
  def=elem0('VectFromDirAtDof',model,data,model.Node(:,1)+.98);
  def.LevelList=[.11 .34 .53];
  cf.def=def;cf.sel='reset';fecom('colordata98');
  cf.ua.clim=[-.01 .01];fecom('colormap','jet(5)');feplot
  
  % Angular position with respect to cylinder
  d2=lsutil('gen',model,{struct('shape','cyla','xc',.5,'yc',.5,'zc',0, ...
    'nx',0,'ny',0,'nz',1)});d2.LevelList=-170:45:170;
  cf.def=d2;cf.sel='reset';fecom('colordata98');cf.ua.clim=[];feplot
  
  RB=struct('Sel','selface');
  RB.def=feutilb('geomrb',model);
  %RB.def=struct('DOF',feutil('getdof',model.Node(:,1),(1:3)'/100));
  %RB.def.def=RB.def.DOF*0;
  cf.def=RB.def;cf.sel='reset';fecom('colordataevalA -alpha.1')
  
  RB.def=d2;sel=lsutil('edgeSelLevelLines',model,RB);cf.SelF{2}=sel;
  %sel=lsutil('edgeSelLevelLines',model,d2,RB);cf.SelF{2}=sel;
%   RB=struct('Elt','selface','gen',{{struct('shape','cyla','xc',.5,'yc',.5,'zc',0, ...
%     'nx',0,'ny',0,'nz',1)}},'LevelList',-170:45:170);
%   sel=lsutil('edgeSelLevelLines',model,RB);cf.SelF{2}=sel;
  cf.o(1)='sel 1 def 1 ch 1 ty1 scc 0.07'; % mesh
  cf.o(2)='sel 2 def 1 ch 1 ty1 scc 0.07'; % mesh
  cla(cf.ga);feplot
  
  
  %feplot(model,def);
  %fecom(';colordata98;colorbar');%;colorscale 2sided;fecom coloralpha');
  %if sdtdef('verm')<804;set(gcf,'renderer', 'zbuffer');end
  %fecom('shownodemark',def.DOF(def.def<0))
  
elseif comstr(Cam,'volsurf')
  %% #VolSurf : generate a line on the surface of a volume
  
  
  model=femesh('testhexa8 struct');
  RB=lsutil('SurfFromPoly',feutilb('surfaceasquad',model,'selface & innode{z>.3}'),struct('start',[5 48],'zlim',[.32,.34]));
  li={struct('shape','DSurf','distFcn',RB.distFcn)};
  lsutil('viewls',model,li)
  RO=struct('tolE',.3,'newTol',.1); mo3=lsutil('cut',model,li,RO);
  lsutil('viewls',mo3,li)
  
elseif comstr(Cam,'mesh');[CAM,Cam]=comstr(CAM,5);
  %% #Mesh -------------------------------------------------------------------
  if nargin>=carg;RO=varargin{carg};else;RO=[];end
  
  if comstr(Cam,'remap');[CAM,Cam]=comstr(CAM,6);
    %% #MeshRemap : geometry for Marc
    
    %geometry %xx to be checked
    RA=struct('dim',[300 165 46.5]*1e-3,'ply1',[1e-4*ones(12,1),[45;0;45;90;-45;0;45;0;45;90;-45;0]],...
      'ply2',[1e-4*ones(8,1),[45;0;45;-45;45;0;45;-45]],'piezo',[[1e-4;1e-3;1e-4],[0;0;0]]);
    %materials %xx to be checked
    
    %% Step1 mesh
    
    mo1=femesh(sprintf('testquad4 divide %d %d 1',RO.ndiv([2 1])  ));mo1.Elt=feutil('setgroup all pro 200',mo1);
    mo1.Node(:,5)=mo1.Node(:,5).*RA.dim(1);mo1.Node(:,6)=mo1.Node(:,6).*RA.dim(2)/2;
    mo2=mo1;mo2.Node(:,6)=mo2.Node(:,6)-RA.dim(2)/2;mo1=feutil('addtest Merge-Edge-NoOri',mo1,mo2);
    %mesh stiffner matid=110 proid=210
    mo2=femesh(sprintf('testquad4 divide %d %d 1',RO.ndiv([3 1])  ));
    mo2.Elt=feutil('setgroup all mat 110',mo2);mo2.Elt=feutil('setgroup all pro 210',mo2);
    mo2.Node(:,5)=mo2.Node(:,5).*RA.dim(1);mo2.Node(:,6)=mo2.Node(:,6).*RA.dim(3);
    mo2.Node(:,6:7)=mo2.Node(:,[7 6]);model=feutil('addtest Merge-Edge-NoOri',mo1,mo2);
    [EltId,model.Elt]=feutil('EltIdFix',model.Elt);
    
    %% Step2 do materials
    %sections  %xx to be checked
    sec1=sprintf('100 %5.2e %6.2f ', RA.ply1');sec2=sprintf('110 %5.2e %4.2f ', RA.ply2');
    secp=sprintf('%d %5.2e %6.2f ', [[1;3;1],RA.piezo]');
    model.il=p_shell(['dbval 200 laminate ',sec1],['dbval 210 laminate ',sec2]);
    
    %pl=m_piezo(pl,'dbval 3 -elas 12 PIC_255');%3-PZT
    E1=[1.58 1.59]'*1e9;G1=[160 67.6]'*1e9; nu1=E1./G1./2-1;
    model.pl=[[100;110] fe_mat('m_elastic','SI',1)*[1;1] E1 [0.3;0.3] G1];
    model.unit='SI';
    %% Step3 add piezo properties
    
    [model,r1]=m_piezo('patchaddpro',model,'BaseId 200 +Noliac.NCE51.OD0TH1.in');% OD0 because rc defined
    [model,r2]=m_piezo('patchaddpro',model,'BaseId 200 +Noliac.NCE51.OD0TH1.in');% OD0 because rc defined
    [model,r3]=m_piezo('patchaddpro',model,'BaseId 200 +Noliac.NCE51.OD0TH1.in');% OD0 because rc defined
    
    feutilb('_write',model)
    
    %il=p_shell(il,);
    %    il=p_shell(il,['dbval 301 laminate ',[sec1,secp]]);il=p_shell(il,['dbval 311 laminate ',[sec2,secp]]);
    %    il=p_piezo(il,sprintf('dbval 300 shell 301 1 %d 0 2 %d 0',size(RA.ply1,1)+[1 3]));
    %    il=p_piezo(il,sprintf('dbval 310 shell 311 1 %d 0 2 %d 0',size(RA.ply2,1)+[1 3]));
    %mesh panel matid=100 proid=200
    
    %add level-set
    li={struct('shape','sphere','xc',RA.dim(1)/3,'yc',RA.dim(2)/4,'zc',0.,'rc',.02,'mpid',[r1.MatId r1.ProId])
      struct('shape','sphere','xc',2*RA.dim(1)/3,'yc',-RA.dim(2)/4,'zc',0.,'rc',.02,'mpid',[r2.MatId r2.ProId])
      struct('shape','sphere','xc',3*RA.dim(1)/4,'yc',0.,'zc',RA.dim(3)/2,'rc',.01,'mpid',[r3.MatId r3.ProId])
      struct('shape','sphere','xc',RA.dim(1)/3,'yc',0.,'zc',RA.dim(3)/2,'rc',.01,'mpid',[r3.MatId r3.ProId])
      struct('shape','sphere','xc',RA.dim(1)/5,'yc',RA.dim(2)/5,'zc',0.,'rc',.02,'mpid',[r2.MatId r2.ProId])};
    %out=lsutil('cut',model,li,RO);
    
    
    
  else;
    [out,li]=ofdemos(['LS' CAM],varargin{2:end});
  end
  out=model; if nargout>1;out1=li;end
  
  
elseif comstr(Cam,'phyline'); [CAM,Cam]=comstr(CAM,8);
 %% #PhyLine: cut physical line with provided seed in model
 if comstr(Cam,'sq')
  mo1=feutil('ObjectAnnulus.5 .5 0 .5 .1 0 0 1 20 20');
  % square to be inserted
  % case 4: much less  mix coarse/fine
  mo2=femesh('testquad4 divide 7 2'); 
  mo2.Node(:,5:7)=mo2.Node(:,5:7)/10+ones(size(mo2.Node,1),1)*[.23 .76 0];
  
  li=lsutil('surfFromRectMesh',mo2);
  fecom('shownodemark',li.contour(:,5:7),'marker','o')
  RC.tolE=.1; mo3=lsutil('cut',mo1,{li},RC);lsutil('ViewLs',mo3,li);fecom('shownodemark',li.contour(:,5:7),'marker','o')
  % need to cleanup corners
  
  %% t_xfem phylinelshape attempts with lines  
  % examples : d_contact simplebrake, aubefan08 xxx
  
  mo2=femesh('testquad4 divide 7 2'); 
  mo2.Node(:,5:7)=mo2.Node(:,5:7)/10+ones(size(mo2.Node,1),1)*[.23 .26 0];
  li=lsutil('SurfFromPoly',mo3,mo2);
  
  %lsutil('ViewLs',mo1,li); % display the level set
  
  mo6=lsutil('cut',mo1,{li})
  lsutil('ViewLs',mo6,li); % display the level set
  
 end
 
  
elseif comstr(Cam,'triafc')
 
 %% basic topo tests
 cf=feplot
 model=femesh('testhexa8')
 model=feutil('hexa2penta',model);
 [u1,model.Elt]=feutil('eltidfix;',model);
 sel=feutil('addsetfaceid -get',model,'skin','selface & eltname quad');
 mo1=lsutil('cutface2tria',model,sel)
 
 model=femesh('testpenta6');
 [u1,model.Elt]=feutil('eltidfix;',model);
 mo1=lsutil('cutface2tria',model,[1 2]); feplot(cf,'initmodel',mo1); fecom(cf,'showpa')
 mo2=lsutil('cutface2tria',model,[1 3]);feplot(cf,'initmodel',mo2); fecom(cf,'showpa')
 mo3=lsutil('cutface2tria',model,[1 5]);feplot(cf,'initmodel',mo3); fecom(cf,'showpa')
 mo3=lsutil('cutface2tria',model,[1 2;1 3]);feplot(cf,'initmodel',mo3); fecom(cf,'showpa')
 mo3=lsutil('cutface2tria',model,[1 5;1 3]);feplot(cf,'initmodel',mo3); fecom(cf,'showpa')
 mo3=lsutil('cutface2tria',model,[1 2;1 3;1 5]);feplot(cf,'initmodel',mo3); fecom(cf,'showpa')



 model=femesh('testhexa8');
 [u1,model.Elt]=feutil('eltidfix;',model);
 mo1=lsutil('cutface2tria',model,[1 1]); feplot(cf,'initmodel',mo1); fecom(cf,'showpa')
 mo1=lsutil('cutface2tria',model,[1 2]); feplot(cf,'initmodel',mo1); fecom(cf,'showpa')
 mo1=lsutil('cutface2tria',model,[1 3]); feplot(cf,'initmodel',mo1); fecom(cf,'showpa')
 mo1=lsutil('cutface2tria',model,[1 4]); feplot(cf,'initmodel',mo1); fecom(cf,'showpa')
 mo1=lsutil('cutface2tria',model,[1 5]); feplot(cf,'initmodel',mo1); fecom(cf,'showpa')
 mo1=lsutil('cutface2tria',model,[1 6]); feplot(cf,'initmodel',mo1); fecom(cf,'showpa')
 mo1=lsutil('cutface2tria',model,[1 1;1 2]); feplot(cf,'initmodel',mo1); fecom(cf,'showpa')
 mo1=lsutil('cutface2tria',model,[1 1;1 4]); feplot(cf,'initmodel',mo1); fecom(cf,'showpa')
 mo1=lsutil('cutface2tria',model,[1 1;1 2;1 4]); feplot(cf,'initmodel',mo1); fecom(cf,'showpa')
 mo1=lsutil('cutface2tria',model,[1 1;1 2;1 3]); feplot(cf,'initmodel',mo1); fecom(cf,'showpa')
 mo1=lsutil('cutface2tria',model,[1 1;1 2;1 3;1 4]); feplot(cf,'initmodel',mo1); fecom(cf,'showpa')
 mo1=lsutil('cutface2tria',model,[1 1;1 3;1 4;1 6]); feplot(cf,'initmodel',mo1); fecom(cf,'showpa')
 mo1=lsutil('cutface2tria',model,[1 1;1 2;1 3;1 4;1 6]); feplot(cf,'initmodel',mo1); fecom(cf,'showpa')
 mo1=lsutil('cutface2tria',model,[1 1;1 2;1 3;1 4;1 5;1 6]); feplot(cf,'initmodel',mo1); fecom(cf,'showpa')
 
 %% tria face cut
 
 RO=struct('dim',[10 10 40],'tolE',.1);
 [model,li]=ofdemos('LS3d',RO);li{1} % Spherical cut
 RO.keepOrigMPID=1;
 [mo3]=lsutil('cut',model,li,RO);
 
 mo4=lsutil('cutface2tria',stack_rm(mo3,'info','EltOrient'),li);
 
 [u1,mpid]=lsutil('mpid',mo4,li);
 
 mo4=feutil('addseteltid',mo4,'cuts',sprintf('eltind %s',num2str(find(mpid(:,1)>0)')));
 
 mo5=mo4;
 [mo4.Elt,mo5.Elt]=feutil('removeelt setname cuts',mo4);
 mo5.Elt=feutil('selelt selface',mo5);
 mo5.Node=feutil('getnodegroupall',mo5);
 
 [i1,st]=sdtdef('in','OpenFEM.tetgen');
 if i1
  mo6=fe_tetgen('mesh -arg pqJYY -keepnum',mo5);
  mo7=feutil('addtestMerge;',mo4,mo6);
 end
 
  %% #CVS ----------------------------------------------------------------------
elseif comstr(Cam,'cvs')
  out='$Revision: 1.68 $  $Date: 2023/01/03 08:27:06 $';
elseif comstr(Cam,'@'); out=eval(CAM);
  %% ------------------------------------------------------------------------
else; error('%s unknown',CAM);
end
