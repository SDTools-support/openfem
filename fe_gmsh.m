function [out,out1]=fe_gmsh(varargin)

% FE_GMSH: gmsh interface, utilities
%
% INSTALL GSMSH : 
%    1. download from GMSH from  http://www.geuz.org/gmsh/
%    2. Tell OpenFEM where to find GMSH using 
%        setpref('OpenFEM','gmsh','/path_to_binary/gmsh.exe')
%
% fe_gmsh('addline',model,sel);
%  Nodes are stored in model.Node
%  Other geometric info GM=stack_get(model,'geom','GMSH','getdata');
% r1=fe_gmsh('lineloops',elt); % finds line loops
%
% fe_gmsh('write FileName.geo [-lc (length)] [-run options]',model)
% Typical options would be -run -2 -v 0 (for 2D meshing non verbode mode)
%
% mode=fe_gmsh('read FileName.msh')   
%

% This function use modifications of the stlread binary STL file reader written 
% by Francis Esmonde-White :
% 
% Copyright (c) 2011, Francis Esmonde-White
% Copyright (c) 2011, Eric Johnson
% All rights reserved.
% 
% This function use stlread code written by Francis Esmonde-White
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the The MathWorks, Inc. nor the names
%       of its contributors may be used to endorse or promote products derived
%       from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%

%       Etienne Balmes, Jean-Philippe Bianchi, Guillaume Vermot des Roches,
%       Guillaume Martin
%       Copyright (c) 2001-2019 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       For revision information use fe_gmsh('cvs')

eM=eMethods;
obj=[]; evt=[];
if nargin==0; 
  if ispref('OpenFEM','gmsh')
   sdtweb('_link',sprintf('system(''%s'')',getpref('OpenFEM','gmsh')),'open')
  end
  sdtweb('fe_gmsh'); 
  return;
elseif ~ischar(varargin{1})
 obj=varargin{1}; evt=varargin{2}; [CAM,Cam]=comstr(varargin{3},1); carg=4;
else; [CAM,Cam]=comstr(varargin{1},1); carg=2;
end
if sp_util('diag');fprintf('fe_gmsh %s \n',CAM);end
%#ok<*NOSEM,*NASGU,*ASGLU,*CTCH>

% ------------------------------------------------------------------------
if comstr(Cam,'set'); [CAM,Cam] = comstr(CAM,4);
 %% #Set------------------------------------------------------------------1 
 if isfield(evt,'RO');RO=evt.RO;
  if isa(obj,'sdth'); UI=feval(obj.data.MainFcn,'PARAMUI',obj); cf=obj;
  elseif isstruct(obj); UI=obj; % Assume UI direclty
  else; errror('Miss UI');
  end
  if isfield(evt,'CAM');CAM=evt.CAM;Cam=lower(CAM);end
 elseif isjava(obj)||isfield(obj,'ob')
   [UI,ua]=eM.initFig(obj); % CinCell callback
   [RO,uo,CAM,Cam]=clean_get_uf('getuo',['SetStruct' CAM],obj,evt);
 else;error('Report EB');
 end
 if ~isfield(UI,'MainFcn'); UI.MainFcn=@sdtroot; end
 out={UI.MainFcn,'InitGMSH'}; % Refresh command
 PARAM=feval(UI.MainFcn,'paramvh',UI);
 if comstr(Cam,'gmsh'); [CAM,Cam]=comstr(CAM,5);
  %% #SetGMSH-------------------------------------------------------------2
  [r1j,r1,st,PARAM,u0]=sdth.eMethods.defaultSet(UI,PARAM,RO,'GMSH',...
   {'Parent','MeshDim','MeshOrder','clmin','clmax','clscale','netgen'...
   'highorder','PostCb'});
  R1=fe_def('cleanentry',r1);
  j1=0;
  while j1<length(st); j1=j1+1; val=RO.(st{j1});
   if strcmpi(st{j1},'FilePut')||strcmpi(st{j1},'FileName'); 
    %% #FilePut : Button callback to get file-----------------------------3
    RO2=struct('title','Choose file to mesh','osdic','FFImGMSHIn','UI',UI);
    if strcmpi(st{j1},'FileName'); RO2.FromScript=1; 
     % Only check that file exists + fill LastWd
     RO2.defname=val;
    end 
    [fname,wd]=sdt_dialogs('getfile',RO2); %GUI
    FileName=fullfile(wd,fname);
    sdcedit(r1j,'FileName',FileName);
   elseif strcmpi(st{j1},'mesh'); 
    %% #Mesh : Compute mesh with GMSH + Execute PostCb--------------------3
    if isempty(R1.FileName); sdtw('_IsGui','Select file first !'); return; end
    if ~strcmpi(val,'FromBusy')
     r3=struct('cbk',...
      {{@fe_gmsh,UI,struct('RO',struct('Mesh','FromBusy'),'CAM','GMSH'),'Set'}}, ...
            'Message',sprintf('Meshing %s with GMSH',R1.FileName), ...
            'FailMessage','Meshing failed !', ...
            'Busy','Blocking');
     feval(cinguj('@hgFe'),[],r3); return;
    else
     RO2=struct('ext','.m'); if fe_gmsh('ver')<400; RO2.ext='.msh'; end
     R1.MeshDim=eval(R1.MeshDim(1)); % 1, 2 or 3
     if strcmpi(R1.MeshOrder,'linear'); R1.MeshOrder=1;
     elseif strcmpi(R1.MeshOrder,'quadratic'); R1.MeshOrder=2;
     end
     sel={'tria3' 'tria6';'tetra4' 'tetra10'};
     if R1.MeshDim>1; 
      RO2.sel=sprintf('selelt eltname %s',sel{R1.MeshDim-1,R1.MeshOrder});
     end
     RO2.Run=sprintf('-%i -order %i ',R1.MeshDim,R1.MeshOrder);
     if ~isempty(R1.clmin); RO2.Run=[RO2.Run sprintf('-clmin %g ',R1.clmin)]; end
     if ~isempty(R1.clmax); RO2.Run=[RO2.Run sprintf('-clmax %g ',R1.clmax)]; end
     if ~isempty(R1.clscale); RO2.Run=[RO2.Run sprintf('-clscale %g ',R1.clscale)]; end
     if ~isempty(R1.netgen); RO2.Run=[RO2.Run '-optimize_netgen ']; end
     if ~isempty(R1.highorder); RO2.Run=[RO2.Run '-optimize_ho ']; end
     % Display command in console
     fprintf('fname=''%s'';\nRO=%s;\nmo1=fe_gmsh(''write'',fname,RO);\n',R1.FileName,comstr(RO2,-30));
     % Execute callback
     mo1=fe_gmsh('write',R1.FileName,RO2);
    end
    if isempty(R1.PostCb)
     cf=clean_get_uf('feplotcf',UI.gf);
     cf.mdl=mo1; fecom(cf,'ShowFiMDef');
    else;
     R1.model=mo1;hgfeval(R1.PostCb,UI,R1)
    end
   end
  end
  cingui('resize',UI.gf);
 else; error('Command "Set%s" unknown',CAM);
 end


elseif comstr(Cam,'add'); [CAM,Cam] = comstr(CAM,4);
%% #Add_ -----------------------------------------------------------------------
model=varargin{carg};carg=carg+1;
r2=varargin{carg};carg=carg+1;
RO=struct('seed',0);

[CAM,Cam,RO.Loop]=comstr('-loop',[-25 1],CAM,Cam);
if carg<=nargin
 data=varargin(carg:end); carg=nargin+1; 
 try data=struct(data{:});
  r3=fieldnames(data); for j1=1:length(r3); RO.(r3{j1})=data.(r3{j1}); end
 catch; error('Wrong additional parameter definition') 
 end
end

GM=stack_get(model,'geom','GMSH','getdata');
if isempty(GM)
 GM=struct('Line',[],'Circle',[],'LineLoop',[], ...
   'PlaneSurface',[],...  % col1 >0 edge, col ... negative = holes, col1<0 ruled
   'SurfaceLoop',[],'Volume',[],'TransfiniteLine',[],'EmbeddedLine',[]);
end

if comstr(Cam,'fullcircle')
%% #AddFullCircle center,edge,normal - - - - - - - - - - - - - - - - - - - - - - - - 
% provides center, one edge node, normal - - - - - - - - - - - - - - - - - -
 %r2=varargin{carg};carg=carg+1;
 if ~isequal(size(r2),[3 3]);error('You must provide center,edgen,normal');end
 r2(3,:)=r2(3,:)/norm(r2(3,:));
 r3=r2(1,:)-r2(2,:);r4=cross(r3,r2(3,:));
 r2=r2([1 1 1 1 1],:)+[0 0 0;r3;r4;-r3;-r4];

 [model.Node,i1]=feutil('addnode',model.Node,r2); i1=model.Node(i1,1); % NodeId
 GM.Circle(end+(1:4),1:3)=i1([2 1 3;3 1 4;4 1 5;5 1 2]);
 RO.LoopInfo={'Circle',size(GM.Circle,1)+(-3:0)};
 
elseif comstr(Cam,'hole')
%% #AddHole center,edge,axis (with proper length) - - -
 r4=r2; r2=r4.node;
 if ~isequal(size(r2),[3 3]);error('You must provide center,edge,length');end
 i1=size(GM.LineLoop,1);
 [model,r1]=fe_gmsh(sprintf('AddFullCircle -loop%i',i1+1),model,r2);
 [model,r3]=fe_gmsh('AddExtrude',model,struct('dir',r2(3,:),'lineloop',i1+1));
 % Remove r3.iLoop(1) from initial node surface 
 % and add r3.iLoop(2) as flat at end, if not open hole
 GM=stack_get(model,'geom','GMSH','getdata');
 for j1=1:length(r4.iSurf)
  if isempty(GM.PlaneSurface)||~any(GM.PlaneSurface(:,1)==r4.iSurf(j1))
   i2=size(GM.PlaneSurface,1)+1;i3=r4.iSurf(j1);
  else;i2=find(GM.PlaneSurface(:,1)==r4.iSurf(j1));
     i3=GM.PlaneSurface(i2,:);i3(i3==0)=[];
  end
  GM.PlaneSurface(i2,1:length(i3)+1)=[i3 -r3.iLoop(j1)]; % add hole in original
  GM.PlaneSurface(GM.PlaneSurface(:,1)<0& ...
      abs(GM.PlaneSurface(:,1)==-r4.iSurf(j1)),:)=[];
 end
 if length(r4.iSurf)==1
  GM.PlaneSurface(end+1,1)=r3.iLoop(2); % Add ruled surface edge of hole
 end
 
 %i3=-setdiff(-GM.PlaneSurface(GM.PlaneSurface(:,1)<0,1),r4.iSurf);
 %GM.PlaneSurface(end+1,1)=r3.iLoop(2);
 model=stack_set(model,'geom','GMSH',GM);
 
elseif comstr(Cam,'circlearc') % 
%% #AddCircleArc - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 %r2=varargin{carg};carg=carg+1;
 [CAM,Cam,i1]=comstr('tangent',[-25 31],CAM,Cam); 
 if ~isempty(i1)&&i1>0
  T=r2(1,:)'; % tangent at first point
  if i1==1; A=r2(2,:)'; B=r2(3,:)'; % tangent T to first node
  elseif i1==2; A=r2(3,:)'; B=r2(2,:)'; % tangent T to 2nd node
  end
  n=cross((B-A),T);n=n/norm(n);
  G=[T';2*(B-A)';n']; % [G]{x}={S}   % {AC}{T}=0, norm(AC)=norm(BC), (AC)^n=0
  S=[dot(A,T);sum(B.^2)-sum(A.^2);dot(A,n)];
  X=pinv(G)*S; % center
  r2(1,:)=X(:)';
 end
 if ~isequal(size(r2,1),3)
  error('You must provide three nodes');
 end
 if size(r2,2)==3 % x y z provided
  [model.Node,i1]=feutil('addnode',model.Node,r2); i1=model.Node(i1,1); % NodeId
 elseif size(r2,2)==1 % NodeId provided
  i1=r2;
 else
  error('You must provide 3 nodes as X Y Z or as NodeId')   
 end
 GM.Circle(end+1,1:3)=i1([2 1 3]);RO.LoopInfo={'Circle',size(GM.Circle,1)};

% Define a surface between two circles - - - - - - - - - - - - - - - - - -
% provide : center, node_radius_1, node_radius_2, normal
elseif comstr(Cam,'disk') % 
%% #adddisk - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

 circle=stack_get(model,'geom','circle','getdata');
 line=stack_get(model,'geom','line','getdata');
 LineLoop=stack_get(model,'geom','LineLoop','getdata');

 %r2=varargin{carg};carg=carg+1;
 if ~isequal(size(r2),[4 3]); error('You must provide 3 nodes and a normal');end

 r3=r2(2:3,:)-r2([1 1],:); % two radius
 for j1=1:3 
  r3(end+1,:)=cross(r2(4,:),r3(end-1,:)); %#ok<AGROW>
  r3(end+1,:)=cross(r2(4,:),r3(end-1,:)); %#ok<AGROW>
 end
 r3=r2(ones(size(r3,1),1),:)+r3;
 [model.Node,i1]=feutil('addnode',model.Node,[r2(1,:);r3]);
 i1=model.Node(i1,1); % NodeId
 circle(end+[1:8],1:3)=i1([2 1 4;4 1 6;6 1 8;8 1 2;3 1 5;5 1 7;7 1 9; ...
    9 1 3]);
 %line(end+[1:4],1:2)=i1([2 3;4 5;6 7;8 9]);
 LineLoop(end+[1:2],:)=size(circle,1)-8+[1 2 3 4;5 6 7 8];

 model=stack_set(model,'geom','circle',circle);
 model=stack_set(model,'geom','line',line);
 model=stack_set(model,'geom','LineLoop',LineLoop);
 out=model;


% Add a straightline between two nodes - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'line') % 
%% #AddLine_ - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if ischar(r2)||iscell(r2)||~isfinite(r2(1)) % Line is specified as a selection
 if iscell(r2);r1=r2;r3=[];
 elseif ischar(r2)
  elt=feutil(horzcat('selelt',r2),model);
  if comstr(Cam,'line3')
   r3=[]; r2=fe_gmsh('lineloops-elt',elt); r4=1;
   i1=find(~isfinite(r2(:,1))); r1=cell(1,length(i1)); i1(end+1)=size(r2,1)+1;
   for j1=1:length(i1)-1
    r1{j1}=r2(i1(j1)+1:i1(j1+1)-1,1:3);
   end
  else; r3=[];r1=fe_gmsh('lineloops',elt);r4=0;
  end
 else;r3=[];r1=fe_gmsh('lineloops',r2); r4=0;
 end
 RO.LoopInfo={};
 for j1=1:length(r1) % list of continuous lines
  if r4 % #AddLine3 - - - - - -
   r2=r1{j1};
   ind=size(GM.Line,1)+[1:size(r2,1)];
   GM.Line(ind,1:size(r2,2))=r2;
  else % #AddLine - - - - - - -
   r2=r1{j1}(:);
   ind=size(GM.Line,1)+[1:length(r2)-1];
   GM.Line(ind,1:2)=[r2(1:end-1) r2([2:end])];
  end
  if isfield(RO,'seed')&&RO.seed % definie specific seed subdivision on lines
   GM.TransfiniteLine=vertcat(GM.TransfiniteLine,[ind' RO.seed(ones(length(ind),1))]);
  end
  if isfield(RO,'embed')&&RO.embed
   GM.EmbeddedLine=vertcat(GM.EmbeddedLine,ind'); r3=[];
  else
   if (size(r2,2)>1&&r2(1,1)==r2(end,2))||r2(1)==r2(end) % LoopLine
    if isfield(RO,'Loop')&&~isempty(RO.Loop)
      RO.LoopInfo(end+1,1:2)={'line',ind};
    else; % Old strategy always addline loop
      GM.LineLoop(end+1,1:length(ind))=ind(:)'; 
    end
    r3(end+1)=size(GM.LineLoop,1); %#ok<AGROW>
   end
  end
 end % list of continuous lines
 if ~isempty(r3); GM.PlaneSurface(end+1,1:length(r3))=r3(:)'; end
elseif isequal(size(r2),[2 3]); % Line specified using two nodes
 [model.Node,i1]=feutil('addnode',model.Node,r2);
 GM.Line(end+1,1:2)=model.Node(i1,1)'; %i1(:)'; % store nodeid
elseif size(r2,2)>=2 % by nodeid [n1 n2;n2 n3 ...]
 %    NNode=[];NNode(model.Node(:,1))=1:size(model.Node,1);
 %    i1=NNode(r2);
 GM.Line(end+1:end+size(r2,1),1:size(r2,2))=r2;
 RO.LoopInfo={'line',size(GM.Line,1)+[-size(r2,1)+1:0]};
 if RO.seed % definie specific seed subdivision on lines : obsolete use col4
  ind=size(GM.Line,1)+1+[-size(r2,1):-1]';
  GM.TransfiniteLine=vertcat(GM.TransfiniteLine,[ind RO.seed(ones(length(ind),1))]);
 end
 %  end
else; error('Not a valid AddLine input format');
end

%% #AddExtrude : extrude a line loop as a series of ruled surfaces - - - - - -
elseif comstr(Cam,'extrude') 
 if isfield(r2,'iPlane')
  %% #AddExtrudeSurf : extrude a surface loop as a series of ruled surfaces - - - - - -
  r1=GM.PlaneSurface(r2.iPlane,:);r1(r1==0)=[];
  r2.SurfaceLoop={r2.iPlane};
  for j1=1:length(r1)
    [model,r3]= ...
        fe_gmsh('AddExtrude',model,struct('lineloop',abs(r1(j1)), ...
        'dir',r2.dir,'AddSurf',1));
    r2.Close(1,j1)=r3.iLoop(2);
    r2.SurfaceLoop{1,end+1}=r3.iRuled;
    %GM=stack_get(model,'geom','GMSH','getdata')
  end
  GM=stack_get(model,'geom','GMSH','getdata');
  GM.PlaneSurface(end+1,1:length(r2.Close))=r2.Close;
  r2.SurfaceLoop=horzcat(r2.SurfaceLoop{:},size(GM.PlaneSurface,1));
  GM.SurfaceLoop(end+1,1:length(r2.SurfaceLoop))=r2.SurfaceLoop;
  GM.Volume(end+1,1)=size(GM.SurfaceLoop,1);

 else
 
 r1=GM.LineLoop(r2.lineloop,:);
 RO.ruled=[];  
 for j1=1:2:length(r1)
  if isempty(r1{j1}); break;end
  r3=r1{j1+1};
  if strcmpi(r1{j1},'line')
   n1=feutil('getnode',model,GM.Line(r3,1:2)'); 
   [model.Node,i1]=feutil('addnode',model.Node, ...
      n1(:,5:7)+ones(size(n1,1),1)*r2.dir(:)');
   i2=model.Node(i1);i1=n1(:,1);
   i3=[i1(1:2:end) i1(2:2:end) i2(2:2:end) i2(1:2:end)];
   for j2=1:size(i3,1); 
    r4={'line' i3(j2,1:2) 'line' i3(j2,2:3) 'line' i3(j2,3:4) ...
        'line' i3(j2,[4 1])};
    GM=AddRuledSurface(GM,r4);RO.ruled(end+1)=size(GM.LineLoop,1);
    r3(j2)=-GM.LineLoop{end,6};
   end
   r1{j1+1}=r3;
  elseif strcmpi(r1{j1},'circle')% start,center,end
   n1=feutil('getnode',model,GM.Circle(r3,1:3)'); 
   [model.Node,i1]=feutil('addnode',model.Node, ...
      n1(:,5:7)+ones(size(n1,1),1)*r2.dir(:)');
   i2=model.Node(i1);i1=n1(:,1); i1=reshape(i1,3,[])';i2=reshape(i2,3,[])';
   r3=zeros(size(i1,1),1);
   for j2=1:size(i2,1)
    r4={'circle' i1(j2,:) 'line' [i1(j2,3) i2(j2,3)] 'circle' i2(j2,[3 2 1]) ...
        'line' [i2(j2,1) i1(j2,1)]};
    GM=AddRuledSurface(GM,r4);RO.ruled(end+1)=size(GM.LineLoop,1);
    r3(j2)=-GM.LineLoop{end,6};
   end
   r1{j1+1}=r3;
  else; error('Not a known extrude type');
  end
 end
 GM.LineLoop(end+1,1:length(r1))=r1;
 RO.iLoop=[r2.lineloop;size(GM.LineLoop,1)];
 if isfield(r2,'AddSurf')&&r2.AddSurf
    RO.iRuled=size(GM.PlaneSurface,1)+(-length(RO.ruled)+1:0);
 end
 end % extrude type 
 
elseif comstr(Cam,'surf') % 
%% #AddSurf Add a surface - - - - - - - - - - - - - - - - - - - - - - - - - - 
if comstr(Cam,'surfaceloop') % 
%% #AddSurfaceLoop - - - - - - - - - - - - - - - - - - - - - - - - - - 
 GM.SurfaceLoop(end+1,1:length(r2))=r2;
 RO.iLoop=size(GM.SurfaceLoop,1);
else   
 %r2=varargin{carg};carg=carg+1;
 if ischar(r2) % Line is specified as a selection
  elt=feutil(horzcat('selelt',r2),model);
 else; elt=r2; 
 end
 [cEGI,un1,name]=getegroup(elt); 
 if length(cEGI)>2; error('Volume to mesh has several groups'); end
 % dirty for triangles
 fun=eval(['@' name{1}]);
 cEGI=2:size(elt,1);r1=fun('edges')';
 r2=elt(cEGI,fun('nodes'))';r2=r2(r1,:);
 r2=reshape(r2,2,numel(r2)/2)';i1=size(r2,1);
 [r3,i2,i3]=unique(sort(r2,2),'rows'); i2=any(r3(i3,:)~=r2,2);i3(i2)=-i3(i2);
 GM.Line=r3;
 GM.LineLoop=reshape(i3,[],length(cEGI))';
 GM.PlaneSurface=(1:length(cEGI))';
 GM.SurfaceLoop=1:size(GM.PlaneSurface,1);
 GM.Volume=1;
end

elseif comstr(Cam,'plane') % 
%% #AddPlaneSurface - - - - - - - - - - - - - - - - - - - - - - - - - - 
 GM.PlaneSurface(end+1,1:length(r2))=r2;
 RO.iPlaneSurface=size(GM.PlaneSurface);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
else;sdtw('''Add%s'' unknow',CAM); % subcommand selection - - - - - - - 
end

% Append to LineLoop accepts multiple and possibly negative to reverse line
if ~isempty(RO.Loop)&&max(RO.Loop)>0&&isfield(RO,'LoopInfo')
 for j1=1:length(RO.Loop)
  i2=RO.Loop(j1);if ~isfield(GM,'LineLoop');GM.LineLoop={};end
  if ~iscell(GM.LineLoop)
    if isempty(GM.LineLoop);GM.LineLoop={};
    else;GM.LineLoop={'line',GM.LineLoop};
    end
    i1=1;  
  elseif abs(i2)>size(GM.LineLoop,1)
   i1=1;  
  else; % Seek if new line loop
     i1=find(cellfun('isempty',GM.LineLoop(abs(i2),:)),1,'first');
     if isempty(i1);i1=size(GM.LineLoop,2)+1;end
  end
  r2=RO.LoopInfo; if i2<0;for j2=2:2:length(r2);r2{j2}=-r2{j2};end;end
  GM.LineLoop(abs(i2),i1-1+(1:length(r2)))=r2;
 end
end
out=stack_set(model,'geom','GMSH',GM);
if nargout>1;out1=RO; end

elseif comstr(Cam,'lineloops')
%% #LineLoops ------------------------------------------------------------------

elt=varargin{carg};carg=carg+1;

[ElemF,st2,ElemP]= feutil('getelemf',elt(1,:),1);
elt(~isfinite(elt(:,1)),:)=[];
i1=strfind(Cam,'-elt');
if ~isempty(i1); RunOpt.out='elt';
else;RunOpt.out='loops';
end
if strcmp(ElemP,'beam3') % allow for mid-node export
      RunOpt.edge=sort(elt(:,1:2),2); RunOpt.edge3=elt(:,1:3);
end

  j2=[];r1=cell(1,0);
  conn=sparse(elt(:,1:2),[1:size(elt,1)]'*[1 1],ones(size(elt,1),1)*[1 2]);
  while any(elt(:));
   if isempty(j2);  % new LineLoop
     j2=min(find(any(conn)));
     r1{end+1}=elt(j2,1:2);elt(j2,:)=0;conn(:,j2)=0; %#ok<AGROW>
   else  % append at the end
    r3=conn(r1{end}(1),:);j2=[];
    if any(r3); 
     [j1,j2,r3]=find(r3); j2=j2(1);r1{end}=[elt(j2,3-r3) r1{end}];
     conn(:,j2)=0;elt(j2,:)=0;
    end
    r3=conn(r1{end}(end),:);
    if any(r3); 
     [j1,j2,r3]=find(r3); j2=j2(1);r1{end}=[r1{end} elt(j2,3-r3)];
     conn(:,j2)=0;elt(j2,:)=0;
    end
   end
  end % while any elt

if strcmp(RunOpt.out,'elt') % return oriented elements
 out=[];
 for j1=1:length(r1)
  out(end+1,1:length(ElemP)+1)=[Inf abs(ElemP)]; %#ok<AGROW>
  r2=[r1{j1}(1:end-1)' r1{j1}(2:end)'];
  if strcmp(ElemP,'beam3')
    [r3,i3,i4]=intersect(sort(r2,2),RunOpt.edge,'rows');
    i2=[];i2(i3)=1:length(i3);%isequal(sort(r3(i2,:),2),sort(r2,2))
    r3=RunOpt.edge3(i4,:);r3=r3(i2,:);r3(:,1:2)=r2; r2=r3;
  end
  out(end+[1:size(r2,1)],1:size(r2,2))=r2;
 end
else
 out=r1;
end

elseif comstr(Cam,'write'); [CAM,Cam] = comstr(CAM,6);
%% #Write ----------------------------------------------------------------------
model=varargin{carg};carg=carg+1;
if carg<=nargin;RunOpt=varargin{carg};carg=carg+1; else; RunOpt=struct;end
if sp_util('issdt')
  [RunOpt,st,CAM]=cingui('paramedit -DoClean',[ ...
         'lc(15#%g# characteristic mesh size")' ...
         'multiple(#3# ")' ...
         'keepContour(#3#"keep contour")' ...
         'import(#3#"keep contour")' ...
         'splineb(#3#"keep contour")' ...
         'spline(#3#"keep contour")' ...
      ],{RunOpt,CAM}); Cam=lower(CAM);
elseif isempty(fieldnames(RunOpt)); % OpenFEM no param in string
 [CAM,Cam,RunOpt.lc]=comstr('-lc',[-25 2],CAM,Cam);
 [CAM,Cam,RunOpt.multiple]=comstr('-multiple',[-25 3],CAM,Cam);
 [CAM,Cam,RunOpt.keepContour]=comstr('-keepcontour',[-25 3],CAM,Cam);
 [CAM,Cam,RunOpt.import]=comstr('-import',[-25 3],CAM,Cam);
 [CAM,Cam,RunOpt.splineb]=comstr('-splineb',[-25 3],CAM,Cam);
 [CAM,Cam,RunOpt.spline]=comstr('-spline',[-25 3],CAM,Cam);
end
if RunOpt.splineb; RunOpt.spline=1; end
RunOpt.ver=fe_gmsh('ver');
if isfield(RunOpt,'ext');
elseif RunOpt.ver<400;RunOpt.ext='.msh';
else;RunOpt.ext='.m'; % Now use Matlab format by default
end
i1=strfind(Cam,'-run');
if ~isempty(i1);
 RunOpt.Run=CAM(i1+4:end);[CAM,Cam]=comstr([CAM(1:i1-1)],1);
end

if ischar(model); fname=model; % Conversion
elseif isempty(Cam) % Write to text
  fid=1; fname=''; 
  if ~isempty(RunOpt.Run);fname=fullfile(sdtdef('tempdir'),'mesh.geo');end
else;fname=CAM;
end
[wd,fname,ext]=fileparts(fname); if isempty(ext);ext='.geo';end
if ~isempty(wd);elseif strcmpi(fname,'del');wd=sdtdef('tempdir');
else;wd=pwd;
end
RunOpt.fgeo=fullfile(wd,horzcat(fname,ext));
RunOpt.oName=fullfile(wd,horzcat(fname,RunOpt.ext));
if ~strcmpi(RunOpt.ext,'.msh'); % Specify output format
 RunOpt.Run=sprintf('%s -o "%s"',RunOpt.Run,RunOpt.oName); 
end
if ischar(model) 
%% this will be a convertion
 fid=-1;
 if isempty([strfind(RunOpt.Run,'-1') strfind(RunOpt.Run,'-2') ...
   strfind(RunOpt.Run,'-3')])
  sdtw('_nb',['When converting it is preferable to specify mesh dimension'...
   '(-1,-2 or -3) in the -run options']);
  RunOpt.Run=[RunOpt.Run ' -3'];
 end
 RunOpt.Run=sprintf('%s -o "%s"',RunOpt.Run,RunOpt.oName); 
else
%% classical write
fid=fopen(RunOpt.fgeo,'w');
GM=stack_get(model,'geom','GMSH','getdata');
if isempty(RunOpt.lc); RunOpt.lc=sdtdef('OpenFEM.DefaultElementLength-safe',.1);end
if strcmp(ext,'.msh') % write a geometry file
    return;
elseif strcmp(ext,'.stl') % #writeSTL: write a stl file - - - - - - -
    
[EGroup,nGroup]=getegroup(model.Elt);
MAP=feutil('getnormal map',model);
[eltid,model.Elt]=feutil('eltidfix',model);
NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
for jGroup=1:nGroup
  cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
  fprintf(fid,'solid Created by fe_gmsh\n');
  for jElt=1:length(cEGI)
    fprintf(fid,'facet normal %.5g %.5g %.5g\n  outer loop\n', ...
        MAP.normal(MAP.ID==eltid(cEGI(jElt)),:));
    fprintf(fid,'  vertex %.5g %.5g %.5g\n', ...
        model.Node(NNode(model.Elt(cEGI(jElt),1:3)),5:7)');
    fprintf(fid,'endloop\nendfacet\n')
  end
end
fprintf(fid,'solid Created by fe_gmsh\n');
fclose(fid);
return;
end %
ind=find(model.Node(:,4)==0);if ~isempty(ind); model.Node(ind,4)=RunOpt.lc;end
%% #Write_Nodes - - - - - - - - 

fprintf(fid,'Point(%i)={%.10g,%.10g,%.10g,%.10g};\n', ...
 model.Node(:,[1 5 6 7 4])');
% #Write_Circles / Id from 1 to size(GM.Circle,1) - - - - -
i0=0;
r1=[];if isfield(GM,'Circle');r1=GM.Circle;end
if ~isempty(r1)
 r1=[i0+[1:size(r1,1)]' r1]; i0=r1(end,1); RunOpt.iCircle=r1(:,1);
 fprintf(fid,'Circle(%i)={%i,%i,%i};\n',r1');
end
%% #Write_Lines / Id from size(GM.Circle,1)+1 to size(GM.Circle,1)+size(GM.Line,1) - - -
% http://geuz.org/gmsh/doc/texinfo/#Command_002dline-options
r1=[];RunOpt.NoTransLine=[];if isfield(GM,'Line');r1=GM.Line;end
if ~isempty(r1)
 GM.iLine=[];
 r2=i0+[1:size(r1,1)]'; %i0=r2(end);
 if isfield(GM,'TransfiniteLine')&&~isempty(GM.TransfiniteLine)
  GM.iTransfiniteLine=r2(GM.TransfiniteLine(:,1));
 end
 if size(r1,2)>2&&RunOpt.spline; i2=find(r1(:,3)); else; i2=find(r1(:,1)); end
 if RunOpt.spline
  r2=[i0+[1:size(r1(i2),1)]' r1(i2,[1 3 2])]; i0=r2(end,1);
  if RunOpt.splineb; stSpline='BSpline(%i)={%i,%i,%i};\n'; 
  else;stSpline='Spline(%i)={%i,%i,%i};\n'; 
  end
  fprintf(fid,stSpline,r2'); GM.iLine=r2(:,1);
  i2=find(r1(:,3)==0); RunOpt.NoTransLine=[]; %r2(:,1);
 end
 if ~isempty(i2)
  r1=[i0+[1:size(r1(i2),1)]' r1(i2,:)]; i0=r1(end,1);
  fprintf(fid,'Line(%i)={%i,%i};\n',r1(:,1:3)'); GM.iLine=[GM.iLine;r1(:,1)];
 end
end
% #Write_Transfinite_Lines - - - - - - r1=[Index n1 n2 n3 ...]
if RunOpt.ver>250; st1='Using Progression 1'; else;st1=''; end
if size(r1,2)>4 % Reimplement with n1,n2,n3/O,ndes
  i3=r1(r1(:,5)>0,[1 5])';
  fprintf(fid,sprintf('Transfinite Line{%%i} = %%i %s;\n',st1),i3);
else
 r1=[];if isfield(GM,'TransfiniteLine');r1=GM.TransfiniteLine;end
 if ~isempty(r1)
  r1=[GM.iTransfiniteLine GM.TransfiniteLine(:,2)];
  if RunOpt.keepContour
   r2=setdiff([size(GM.Line,1)+1:size(GM.Circle,1)+size(GM.Line)]',...
    [GM.iTransfiniteLine;RunOpt.NoTransLine]);
   r1=[r1;r2 2*ones(length(r2),1)];
  end
  fprintf(fid,sprintf('Transfinite Line{%%i} = %%i %s;\n',st1),r1');
 elseif RunOpt.keepContour;
  r1=setdiff(GM.iLine,RunOpt.NoTransLine);
  if ~isempty(r1)
   r1=[r1 repmat(2,length(r1),1)];
   fprintf(fid,sprintf('Transfinite Line{%%i} = %%i %s;\n',st1),r1');
  end
 end
end

%% #Write_LineLoops - - - - - - - - -
r1=[];if isfield(GM,'LineLoop');r1=GM.LineLoop;end;
GM.iLineLoop=i0;
if ~isempty(r1)&&iscell(r1) % LineLoop=cell array {'type',ind,'type',ind,...}
 if iscell(r1{1})
  r2=cell(size(r1,1),max(cellfun(@length,r1)));
  for j1=1:size(r2,1); r2(j1,1:size(r1{j1},2))=r1{j1}; end
  r1=r2;
 end
 for j1=1:size(r1,1)
  i2=i0+j1; % id of lineloop
  for j2=1:2:nnz(double(~cellfun('isempty',r1(j1,:)))) % loop on non empty
   if isempty(r1{j1,j2}); break; end
   i3=r1{j1,j2+1}; % Id in gmsh
   if comstr(lower(r1{j1,j2}),'line'); 
       i3=sign(i3(:)).*GM.iLine(abs(i3)); 
   end% i3=i3+size(GM.Circle,1); end
   % shift i3 indices for other line type if add any xxx...
   i2=[i2 i3(:)'];
  end
  fprintf(fid,'Line Loop(%i)={%i%s};\n',i2(1:2),sprintf(',%i',i2(3:end)));
 end
 GM.iLineLoop=i0+1:i0+size(r1,1);
 i0=i0+size(r1,1); % inc id 
elseif ~isempty(r1) % LineLoop=matrix
 GM.iLine=GM.iLine(:);
 i2=r1~=0; 
 if all(all(i2))
  r1(i2)=reshape(reshape(sign(r1(i2)),[],1).*...
  GM.iLine(abs(r1(i2))),size(r1,1),size(r1,2));
 else
  r1(i2)=sign(r1(i2)).*GM.iLine(abs(r1(i2)));
 end
 r1=[i0+[1:size(r1,1)]' r1];i0=r1(end,1);
 for j1=1:size(r1,1)
  i2=r1(j1,r1(j1,:)~=0);
  fprintf(fid,'Line Loop(%i)={%i%s};\n',i2(1:2),sprintf(',%i',i2(3:end)));
 end
 GM.iLineLoop=r1(:,1)'; i0=r1(end,1);
end
% #Write_PlaneSurfaces and RuleSurface one per closed line loop - - - - -
r1=[];if isfield(GM,'PlaneSurface');r1=GM.PlaneSurface; end
GM.iPlaneSurface=[];
if ~isempty(r1)
 if RunOpt.multiple
  r1=[i0+[1:size(r1,1)]' GM.iLineLoop'];i0=abs(r1(end,1));
  fprintf(fid,'Plane Surface(%i)={%i};\n',r1');
  GM.iPlaneSurface=r1(:,1);
 else
  for j1=1:size(r1,1)
   i0=i0+1; % increment gmsh ID
   i2=r1(j1,:);i2=i2(i2~=0);
   if isempty(i2)
   elseif i2(1)<0 % Write ruled surface
    GM.iPlaneSurface=[GM.iPlaneSurface GM.iLineLoop(-i2)];%i0=i0-1;
    fprintf(fid,'Ruled Surface(%i) = {%i};\n',GM.iLineLoop(-i2)*[1 1]);
   else % write plane surface
    i2=(abs(i2)-1+min(GM.iLineLoop)).*sign(i2);
    st=repmat('%i,',1,length(i2)); st(end)='';
    st=sprintf(st,i2);
    fprintf(fid,'Plane Surface(%i)={%s};\n',i0,st);
    GM.iPlaneSurface=[GM.iPlaneSurface i0];
   end
  end
 end
end

% #Write_Embedded_Lines in surfaces (limited to surface 1) - - - - -
r1=[];if isfield(GM,'EmbeddedLine');r1=GM.EmbeddedLine;end
if ~isempty(r1)&&isfield(GM,'iPlaneSurface')&&~isempty(GM.iPlaneSurface)
 r1=[r1 GM.iPlaneSurface(1)*ones(size(r1))]';
 fprintf(fid,'Line {%i} In Surface {%i};\n',r1);
end
% #Write_SurfaceLoops - - - - - - - - 
r1=[];if isfield(GM,'SurfaceLoop');r1=GM.SurfaceLoop;end
GM.iSurfaceLoop=[];
if ~isempty(r1)
 r2=r1;r2(r2~=0)=GM.iPlaneSurface(r2(r2~=0));
 r1=[i0+[1:size(r1,1)]' r2.*sign(r1)];i0=abs(r1(end,1));
 for j1=1:size(r1,1)
  i2=r1(j1,:);i2=i2(i2~=0);
  if length(i2)==2;st='';else;st=sprintf(',%i',i2(3:end));end
  fprintf(fid,'Surface Loop(%i)={%i%s};\n',i2(1:2),st);
 end
 GM.iSurfaceLoop=r1(:,1);
end
% #Write_Volumes - - - - - - - - - - 
r1=[];if isfield(GM,'Volume');r1=GM.Volume;end
if ~isempty(r1)
 r1(r1~=0)=GM.iSurfaceLoop(r1(r1~=0));
 r1=[i0+[1:size(r1,1)]' r1];i0=abs(r1(end,1));
 for j1=1:size(r1,1)
  i2=r1(j1,:);i2=i2(i2~=0);
  if length(i2)==2;st='';else;st=sprintf(',%i',i2(3:end));end
  fprintf(fid,'Volume(%i)={%i%s};\n',i2(1:2),st);
 end
end
if isfield(GM,'Post'); %% #WritePost commands at the end of file
    fprintf(fid,'%s\n',GM.Post{:});
end
end

if fid~=-1; fclose(fid);end
if isfield(RunOpt,'Run'); 
%% #Run : actually run results del.geo to delete on exit - - - - - - - - - - - - -
 RunOpt.Run=sprintf('"%s" "%s" %s', ...
  sdtdef('OpenFEM.gmsh-safe','gmsh.exe'),RunOpt.fgeo,RunOpt.Run);
 fprintf('\nStarting GMSH ...'); 
 pw0=pwd;cd(wd);[i1,RunOpt.log]=system(RunOpt.Run);cd(pw0);
 if ~isempty(RunOpt.oName);
   fprintf('done, reading %s \n',RunOpt.oName);
   out=fe_gmsh(sprintf('read %s',RunOpt.oName));
 else;  fprintf('done\n');
 end
 % .sel select some elements (volumes)
 if isfield(RunOpt,'sel')&&~isempty(RunOpt.sel)
  if ~comstr(lower(RunOpt.sel),'selelt')
      RunOpt.sel=['selelt',RunOpt.sel];
  end
  elt=feutil(RunOpt.sel,out);
  if ~isempty(elt);out.Elt=elt; 
  else;error('Empty selection %s',RunOpt.sel);
  end
 end
 % .pl set properties
 if isfield(RunOpt,'pl'); 
   out.pl=RunOpt.pl; out.il=[];
   out=stack_set(out,'info','GmshId',feutil('proid',out));
   out.Elt=feutil(sprintf('set groupall Mat%i Pro%i',out.pl([1 1])),out);
   out=p_solid('default;',out);
 end
 if ~isempty(RunOpt.log);out=stack_set(out,'info','MeshLog',RunOpt.log);end
 
 if strcmp(fname,'del.') % Delete temporary files if del
   delete(fullfile(wd,'del.geo'));delete(fullfile(wd,'del.msh'))
 end
end

elseif comstr(Cam,'read'); [CAM,Cam] = comstr(CAM,5);
%% #Read -----------------------------------------------------------------------
[CAM,Cam,RunOpt.Vol]=comstr('-vol',[-25 3],CAM,Cam);
if ~isempty(CAM)
elseif nargin>1&&ischar(varargin{2});CAM=varargin{2};
elseif nargin>1&&isstruct(varargin{2})&&isfield(varargin{2},'FileName');
 CAM=varargin{2}.FileName;
else;
    [fname,wd]=uigetfile({'*.msh','GMSH mesh :';'*.stl;*.STL','STL file';...
     '*.*','All Files'});
    if isnumeric(fname); return; end
    CAM=fullfile(wd,fname);
end
if ~exist(CAM,'file')
 error('Can''t find file %s',CAM)
end
fid=fopen(CAM,'r');[wd,fname,ext]=fileparts(CAM);

%% #Read.MSH files - - - - - - - - - - - - - - - -
if comstr(ext,'.msh') 

st=fgetl(fid);
model=struct('Node',[],'Elt',[]);

while ischar(st);

st1=sscanf(st,'%s',1);
switch upper(st1)
case {'$MESHFORMAT'} % Reading nodes
  RunOpt.Format=fscanf(fid,'%f');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
case {'$NOD','$NODES'} % Reading nodes

if RunOpt.Format(1)>=4
 if RunOpt.Format(1)>=4.1; RunOpt.nI=fscanf(fid,'%i',4);
 else; RunOpt.nI=fscanf(fid,'%i',2); % entity numnodes
 end
 while size(model.Node,1)<RunOpt.nI(2)
 i1=fscanf(fid,'%i',4);%tagEntity(int) dimEntity(int) parametric(int; see below) numNodes
 if i1(3)~=0; error('Not implemented');
 elseif i1(4)==0; continue;
 else % single node 
  r1=textscan(fid,'%n%n%n%n',i1(4));r1=horzcat(r1{:});
  % else
  %r1(:,2:4)=fscanf(fid,'%g %g %g',[3 i1(4)])';r1(:,5)=0;
 end
 r1(:,5)=i1(1)+i1(2)/10; % id of entity / dim 
 model.Node=[model.Node;r1(:,[1 5 5 5 2 3 4])]; 
 end
else
 i1=fscanf(fid,'%i',1); % number of nodes
 r1=fscanf(fid,'%i %g %g %g',[4 i1])';r1(:,5)=0;
 model.Node=r1(:,[1 5 5 5 2 3 4]); 
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
case {'$ENTITIES'} % Reading entities in 4.0 format 
  r1=fscanf(fid,'%g');  
case {'$ELM'} % Reading elements format 1.
 %elm-number elm-type reg-phys reg-elem number-of-nodes node-number-list
 %elm-number elm-type number-of-tags node-number-list
 ne=fscanf(fid,'%i',1); % number of elements
 i1=fscanf(fid,'%i');
 while ~isempty(i1)
  switch i1(2)
  case 1;  ElemF='beam1';   i3=7;   ind=1:2;
  case 2;  ElemF='tria3';   i3=8;   ind=1:3;
  case 3;  ElemF='quad4';   i3=9;   ind=1:4;
  case 4;  ElemF='tetra4';  i3=9;   ind=1:4;
  case 5;  ElemF='hexa8';   i3=13;  ind=1:8;
  case 6;  ElemF='penta6';  i3=11;  ind=1:6;
  case 7;  ElemF='pyra5';   i3=10;  ind=1:5;
  case 8;  ElemF='beam3';   i3=8;   ind=1:3;
  case 9;  ElemF='tria6';   i3=11;  ind=1:6;
  case 10; ElemF='quad9';   i3=14;  ind=1:9;
  case 11; ElemF='tetra10'; i3=15;  ind=[1:8 10 9];
  case 12; ElemF='hexa20';  i3=25;  ind=1:20;
  case 13; ElemF='penta15'; i3=20;  ind=1:15;
  case 14; ElemF='pyra13';  i3=18;  ind=1:13;
  case 15; ElemF='mass2';   i3=6;   ind=1;
      otherwise; error(1)
  end

  i2=floor(length(i1)/i3);i2=reshape(i1(1:i2*i3),i3,i2)';
  i3=find(diff(i2(:,2))~=0); if ~isempty(i3);i2=i2(1:min(i3),:); end
  i1(1:numel(i2))=[];
  model.Elt(end+1,1:length(ElemF)+1)=[Inf abs(ElemF)];
  model.Elt(end+[1:size(i2,1)],1:size(i2,2)-3)=i2(:,[6:end 3 4]);
 end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
case {'$ELEMENTS'} % Reading elements GMSH

if RunOpt.Format(1)>=4; 
 st=fgetl(fid);
 RunOpt.iElt=sscanf(st,'%i'); % first line
% numEntityBlocks(unsigned long) numElements(unsigned long)
% tagEntity(int) dimEntity(int) typeEle(int; see below) numElements(unsigned long)
% tag(int) numVert(int) ...
     
 i1=fscanf(fid,'%i');RunOpt.Last='';jElt=1;
 while ~isempty(i1)
  switch i1(3) % i3 gives number of values
  case 1;  ElemF='beam1';   ind=[1:2];
  case 2;  ElemF='tria3';   ind=[1:3];
  case 3;  ElemF='quad4';   ind=[1:4];
  case 4;  ElemF='tetra4';  ind=[1:4];
  case 5;  ElemF='hexa8';   ind=[1:8];
  case 6;  ElemF='penta6';  ind=[1:6];
  case 7;  ElemF='pyra5';   ind=[1:5];
  case 8;  ElemF='beam3';   ind=[1:3];
  case 9;  ElemF='tria6';   ind=[1:6];
  case 10; ElemF='quad9';   ind=[1:9];
  case 11; ElemF='tetra10'; ind=[1:8 10 9];
  case 12; ElemF='hexa20';  ind=[1:20];
  case 13; ElemF='penta15'; ind=[1:15];
  case 14; ElemF='pyra13';  ind=[1:13];
  case 15; ElemF='mass2';   ind=[1];
      otherwise; error('Problem Element %s',comstr(i1,-30))
  end
  if ~isequal(ElemF,RunOpt.Last);
    model.Elt(jElt,1:length(ElemF)+1)=[Inf abs(ElemF)];jElt=jElt(end)+1;
    RunOpt.Last=ElemF;
  end
  i2=reshape(i1(4+(1:(1+length(ind))*i1(4))),[],i1(4))';
  i2=i2(:,[1 1+ind]); % reorder nodes for SDT convention
  i2=i2(:,[2:end 1 1 1]); 
  i2(:,end-2)=i1(1);% tagEntity-MatId
  i2(:,end-1)=i1(2);% DimEntity->ProId
  jElt=jElt(end)+(0:size(i2)-1);model.Elt(jElt,1:size(i2,2))=i2; 
  jElt=jElt(end)+1;
  i1(1:4+size(i2,1)*(size(i2,2)-2))=[];
 end
    
else; 
 %number-of-elements
 %elm-number elm-type number-of-tags < tag > ... node-number-list
    
 ne=fscanf(fid,'%i',1); % number of elements
 i1=fscanf(fid,'%i');
 while ~isempty(i1)
  switch i1(2)
  case 1;  ElemF='beam1';   i3=7;   ind=[1:2];
  case 2;  ElemF='tria3';   i3=8;   ind=[1:3];
  case 3;  ElemF='quad4';   i3=9;   ind=[1:4];
  case 4;  ElemF='tetra4';  i3=9;   ind=[1:4];
  case 5;  ElemF='hexa8';   i3=13;  ind=[1:8];
  case 6;  ElemF='penta6';  i3=11;  ind=[1:6];
  case 7;  ElemF='pyra5';   i3=10;  ind=[1:5];
  case 8;  ElemF='beam3';   i3=9;   ind=[1:3];
   %ElemF='beam1';   i3=8;   ind=[1:2];  % corrected Mar, 2009,was probably wrong 
  case 9;  ElemF='tria6';   i3=11;  ind=[1:6];
  case 10; ElemF='quad9';   i3=14;  ind=[1:9];
  case 11; ElemF='tetra10'; i3=15;  ind=[1:8 10 9];
  case 12; ElemF='hexa20';  i3=25;  ind=[1:20];
  case 13; ElemF='penta15'; i3=20;  ind=[1:15];
  case 14; ElemF='pyra13';  i3=18;  ind=[1:13];
  case 15; ElemF='mass2';   i3=6;   ind=[1];
      otherwise; error(1)
  end

  i3=3+i1(3)+length(ind);
  %in1=[i1(3)+3+1:i3 4:i1(3)+3];
  in1=[i1(3)+3+ind(:)' 4:i1(3)+3];
  i2=floor(length(i1)/i3);i2=reshape(i1(1:i2*i3),i3,i2)';
  i3=find(diff(i2(:,2))~=0); if ~isempty(i3);i2=i2(1:min(i3),:); end
  i1(1:numel(i2))=[];
  model.Elt=feutil('addelt',model.Elt,ElemF,i2(:,in1));
 end
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
otherwise; 
  if comstr(st1,'$END')||comstr(st1,'$End')
  elseif  ~isempty(st1); fprintf('''%s'' ignored\n',st1);
  end
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
st=fgetl(fid);
end % loop on file
if RunOpt.Vol % Keep only volumes
 model.Elt=feutil('removeelt eltname mass',model);
 model.Elt=feutil('removeelt eltname beam',model);
 model.Elt=feutil('removeelt eltname quad',model);
 model.Elt=feutil('removeelt eltname tria',model);
end

elseif strcmpi(ext,'.geo') 
%% #ReadGeo
model=struct('Node',[],'Elt',[]);
St=' ';GM=struct('Node',{cell(1000,5)},'Line',{cell(1000,1)}, ...
    'LineLoop',{cell(1000,1)});
RO.iNode=1;RO.iLine=1;
RO.beam=[];
while 1
 St=fgetl(fid);if ~ischar(St);break;end;[St,st]=comstr(St,1);
 if comstr(st,'point')
   st=textscan(St,'Point%[^)]%*c=%*c%n,%n,%n,%[^}]',St);
   st{1}=st{1}{1};st{1}(1)='';
   if length(st)<5||isempty(st{5});st{5}='';else;st{end}=st{end}{end};end
   GM.Node(RO.iNode,:)=st;RO.iNode=RO.iNode+1;
 elseif comstr(st,'spline')
   st=textscan(St,'Spline%[^)]%*c=%[^}]',St);
   st{1}=st{1}{1};st{1}(1)='';st{end}=st{end}{end};
   GM.Line(RO.iLine,1:3)={'spline',st{1:2}};RO.iLine=RO.iLine+1;  
 elseif comstr(st,'line')
  if comstr(st,'line loop')
    disp(St);
  else
   st=textscan(St,'Line%[^)]%*c=%*c%n,%n,%[^}]',St);
   st{1}=st{1}{1};st{1}(1)='';
   GM.Line(RO.iLine,1:3)={'line',st{1},[st{2:3}]};RO.iLine=RO.iLine+1;  
   RO.beam(end+1,:)=[st{2:3}];
  end
 elseif comstr(st,'circle') % {n1,cent,n2}
   st=textscan(St,'Circle%[^)]%*c=%*c%[^,],%[^,],%[^}]',St);
   for j1=1:4;st{j1}=st{j1}{1};end;st{1}(1)='';
   GM.Line(RO.iLine,1:5)={'circle',st{1:4}};RO.iLine=RO.iLine+1;
 else; 
    r2=[];
    try; 
      eval(sprintf('r2.%s',St));
    end
    if isstruct(r2);st1=fieldnames(r2);RO.(st1{1})=r2.(st1{1});
    else
        disp(St)
    end
 end
end
r1=reshape(vertcat(GM.Node{:,2:4}),[],3);
model.Node=[(1:size(r1,1))'*[1 0 0 0] r1];
[i2,i3]=ismember(cellfun(@(x)sprintf('%i',x),num2cell(RO.beam(:)), ...
    'UniformOutput',false),GM.Node(:,1));RO.beam=reshape(i3,[],2); %#ok<ASGLU>
model=feutil('addelt',model,'beam1',RO.beam);
model=stack_set(model,'geom','GMSH',GM);

%% #ReadStl
elseif comstr(lower(ext),'.stl')
 if fid == -1 
     error('File could not be opened, check name or path.')
 end
 ftitle=fread(fid,80,'uchar=>schar'); % Read file title
 if strcmpi(char(ftitle'),'solid') % ASCII file
  error('ASCII STL file not handled yet');
 end
 [vert,face]=stlread(fid);
 model.Node=feutil('AddNode -nocoal',[],vert); %Keep all node even within epsl
 model.Elt=feutil('AddElt',[],'tria3',face);
elseif comstr(lower(ext),'.bdf')
 model=nasread(fopen(fid));
elseif comstr(lower(ext),'.m')
 st=fscanf(fid,'%c');
 eval(st);
 model=struct('Node',[(1:size(msh.POS,1))'*[1 0 0 0] msh.POS], ...
     'Elt',[]);
 li={'LINES3','beam3';'TRIANGLES6','tria6';'TETS10','tetra10';...
  'LINES','beam1';'TRIANGLES','tria3';'QUADS','quad4';...
  'HEXAHEDRA','hexa8'};
 i1=ismember(li(:,1),fieldnames(msh));
 st=li(i1,:);
 for j1=1:size(st,1)
  if strcmpi(st{j1,1},'TETS10');
   msh.(st{j1,1})(:,[9 10])=msh.(st{j1,1})(:,[10 9]); % Permute nodes for SDT convention
  end
  model.Elt=feutil('addelt',model.Elt,st{j1,2},msh.(st{j1,1}));
 end

else; sdtw('''%s'' unknow format',ext); % subcommand selection - - - - - - - 
end

if fid~=1; fclose(fid);end
if nargout==0; feplot(model);else; out=model;end


elseif comstr(Cam,'ver')
%% #VER : check gmsh version -------------------------------------------------
fname=sdtdef('OpenFEM.gmsh');
if isdir(fname);error('''%s'' should be a path to the gmsh executable',fname);
end
try [i1,out]=system(sprintf('%s -info',fname));  
catch; out=[]; i1=0;
end
if ~isempty(strfind(out,'error using'))%||i1; 
 disp(out)
 warning(['The call to gmsh using ''%s'' generates an error.' ...
        ' You should check your gmsh setup ' ...
        ' (see preference in getpref(''OpenFEM'',''gmsh''))'],... 
        sdtdef('OpenFEM.gmsh'));
end
if ~isempty(strfind(out,'command not found'))||...
   ~isempty(strfind(out,'pas reconnu')); 
 out=[]; 
end
out=out(strfind(out,'Version'):end);
if isempty(out)
 error('\nfe_gmsh could not get gmsh version from ''%s'' \n  %s', ...
  sdtdef('OpenFEM.gmsh'),'you should check your gmsh setup')
end

out=textscan(out,'%s %s','delimiter','\t\b:-');
out=textscan(strrep(out{2}{...
 strncmp('Version',cellfun(@deblank,out{1},'uniformoutput',0),1)},'-svn',''),...
 '%s','delimiter','\t\b.'); out=cellfun(@str2double,(out{1}));
out=sum(out.*flipud(logspace(0,length(out)-1,length(out))'));

%% #end ----------------------------------------------------------------------
elseif comstr(Cam,'cvs')
 out='$Revision: 1.99 $  $Date: 2022/06/03 16:49:29 $';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
else ; sdtw('''%s'' unknow',CAM); % subcommand selection - - - - - - - 
end % function

%% #eMethods - - -------------------------------------------------------------
function out=eMethods
persistent eM
if isempty(eM)
 r1=dbstack;eM=[];
 if sum(ismember({r1.name},'eMethods'))>1||sum(ismember({r1.name},mfilename))>3;
 else
  eM=struct('initFig',sdtroot('@initFig'),...
   'defaultSetCf',sdtroot('@defaultSetCf'), ...
   'defaultSet',sdtroot('@defaultSet'), ...
   'GetQual',id_rc('@GetQual'),...
   'rms',id_rc('@rms'),...
   'tabChange',cinguj('@tabChange'),...
   'isInDock',iimouse('@isInDock'),...
   'dockGroup',iimouse('@dockGroup'),...
   'iimouse_scroll',iimouse('@iimouse_scroll'),...
   'test_exist',idcom('@test_exist'),...
   'idcom_ident',idcom('@idcom_ident'),...
   'check_cf',iicom('@check_cf'));
 end
end
out=eM;

%% #AddRuledSurface a closed LineLoop and associated ruled surface
function  GM=AddRuledSurface(GM,i3); 

for j1=1:2:length(i3)
 if strcmpi(i3{j1},'line') % Add line with existing check
  i1=i3{j1+1};i2=find(ismember(GM.Line,i1,'rows'));
  if isempty(i2); i2=-find(ismember(GM.Line,i1([2 1]),'rows'));end
  if isempty(i2);
      i2=size(GM.Line,1)+1;GM.Line(i2,:)=i3{j1+1};
  end
  i3{j1+1}=i2;
 elseif strcmpi(i3{j1},'circle')
  i1=i3{j1+1};i3{j1+1}=zeros(size(i1,1),1);
  for j2=1:size(i1,1)
   i2=find(ismember(GM.Circle,i1(j2,:),'rows'));
   if isempty(i2); i2=-find(ismember(GM.Circle,i1(j2,[3 2 1]),'rows'));end
   if isempty(i2);
      i2=size(GM.Circle,1)+1;GM.Circle(i2,:)=i1(j2,:);
   end
   if length(i2)~=1;error('Incorrect length');end
   i3{j1+1}(j2)=i2;
  end
 else; error('Not valid');
 end
end
GM.LineLoop(end+1,1:length(i3))=i3;
GM.PlaneSurface(end+1,1)=-size(GM.LineLoop,1);

%% #stlread
function [v, f, n, c] = stlread(fid, verbose)
%[v, f, n, c, stltitle] = stlread(filename, verbose);
% v contains the vertices for all triangles [3*n x 3].
% f contains the vertex lists defining each triangle face [n x 3].
% n contains the normals for each triangle face [n x 3].
% c is optional and contains color rgb data in 5 bits [n x 3].
% stltitle contains the title of the specified stl file [1 x 80].
%  http://www.esmonde-white.com/home/diversions/matlab-program-for-loading-stl-files
%
% Based on code originally written by:
%    Doron Harlev
% and combined with some code by:
%    Eric C. Johnson, 11-Dec-2008
%    Copyright 1999-2008 The MathWorks, Inc.
%
% Re-written and optimized by Francis Esmonde-White, May 2010.
use_color=(nargout>=4);
if ~exist('verbose','var')
 verbose = false;
end
numFaces=fread(fid,1,'int32'); % Read number of Faces
T = fread(fid,inf,'uint8=>uint8'); % read the remaining values
if verbose
 fprintf('\nTitle: %s\n', stltitle);
 fprintf('Number of Faces: %d\n', numFaces);
 disp('Please wait...');
end
% Each facet is 50 bytes
%  - Three single precision values specifying the face normal vector
%  - Three single precision values specifying the first vertex (XYZ)
%  - Three single precision values specifying the second vertex (XYZ)
%  - Three single precision values specifying the third vertex (XYZ)
%  - Two color bytes (possibly zeroed)

% 3 dimensions x 4 bytes x 4 vertices = 48 bytes for triangle vertices
% 2 bytes = color (if color is specified)
trilist = 1:48;
ind = reshape(repmat(50*(0:(numFaces-1)),[48,1]),[1,48*numFaces])+repmat(trilist,[1,numFaces]);
Tri = reshape(typecast(T(ind),'single'),[3,4,numFaces]);
n=squeeze(Tri(:,1,:))';
n=double(n);
v=Tri(:,2:4,:);
v = reshape(v,[3,3*numFaces]);
v = double(v)';
f = reshape(1:3*numFaces,[3,numFaces])';
if use_color
 c0 = typecast(T(49:50),'uint16');
 if (bitget(c0(1),16)==1)
  trilist = 49:50;
  ind = reshape(repmat(50*(0:(numFaces-1)),[2,1]),[1,2*numFaces])+repmat(trilist,[1,numFaces]);
  c0 = reshape(typecast(T(ind),'uint16'),[1,numFaces]);
  
  r=bitshift(bitand(2^16-1, c0),-10);
  g=bitshift(bitand(2^11-1, c0),-5);
  b=bitand(2^6-1, c0);
  c=[r; g; b]';
 else
  c = zeros(numFaces,3);
 end
end
[v, f]=patchslim(v, f);
if verbose
 disp('Done!');
end

%% #patchslim
function [vnew, fnew]=patchslim(v, f)
% This function finds and removes duplicate vertices.
% [v, f]=patchslim(v, f)
% Francis Esmonde-White, May 2010

if ~exist('v','var')
 error('The vertex list (v) must be specified.');
end
if ~exist('f','var')
 error('The vertex connectivity of the triangle faces (f) must be specified.');
end

[vnew, indexm, indexn] =  unique(v, 'rows');
fnew = indexn(f);
