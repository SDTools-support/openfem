function [out,out1,out2,out3,out4]=lsutil(varargin)

%LSUTIL Level set utilities
%
% Supported analytic level sets :
%    rect(lx,ly,xc,yc,alpha) circ(rc,xc,yc)
%    box(lx,ly,lz,xc,yc,zc,nx,ny,nz); sphere(rc,xc,yc,zc))
%    cyl(xc,yc,zc,nx,ny,nz,rc,z0,z1)
%    toPlane(xc,yc,zc,nx,ny,nz)
%
%  Note for remeshing, the convention is : - outside, + inside
%
% Contributed by Eric Monteiro, Etienne Balmes : ENSAM / PIMM
%                Guillaume Vermot des Roches : SDTools

%       Copyright (c) 2001-2020 by SDTools & ENSAM, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       For revision information use feutil('cvs')

%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*NOSEM,*NBRAK,*SEPEX>

if ischar(varargin{1});obj=[];evt=[];[CAM,Cam]=comstr(varargin{1},1);carg=2;
else;obj=varargin{1};evt=varargin{2};[CAM,Cam]=comstr(varargin{3},1);carg=4;
end

%% #Gen : analytic generation of typical level set functions -----------------
if comstr(Cam,'gen'); [CAM,Cam] = comstr(CAM,4);
 %if (+ outside, - inside): max for intersection, min for union
 %if (- outside, + inside): min for intersection, max for union
 %SDT convention: - outside/exterior, + inside/interior
 model=varargin{carg};carg=carg+1;
 R1=varargin{carg};carg=carg+1;
 
 if iscell(R1)&&ischar(R1{1})
  R1=cell2struct(R1(2:end,:),R1(1,:),2);
  phi=arrayfun(@(x)lsutil('gen',model,x),R1,'uni',0);
  out=struct('def',horzcat(phi{:}),'DOF',model.Node(:,1)+.98);
 elseif iscell(R1)
  %% expected {R1}
  phi=cell(length(R1),1);
  for j1=1:numel(R1)
   if numel(R1)==1&&isfield(R1{1},'distFcn');
       phi{j1}=R1{1}.distFcn(model.Node(:,5:7));continue;
   end
   if isfield(R1{j1},'idc') % idc : allow definition of center as NodId
    r2=model.Node(model.Node(:,1)==R1{j1}.idc,5:7);
    R1{j1}.xc=r2(1);R1{j1}.yc=r2(2);R1{j1}.zc=r2(3);
   end
   try;
    st1=lower(R1{j1}.shape);if strncmpi(st1,'to',2);st1(1:2)='';end
    st1=sprintf('@dTo%s',st1);st1(5)=upper(st1(5));
    %% #generic handling of dTo(Shape) calls
    fun=eval(st1);r2=functions(fun);
    if ~isempty(r2.file);
     try; R1{j1}=feval(fun,'init',R1{j1},model);
      if isfield(R1{j1},'pl')&&size(R1{j1}.pl,1)==1;
        R1{j1}.mpid=R1{j1}.pl(1)*[1 1];
      end
      if isfield(R1{j1},'il')&&size(R1{j1}.il,1)==1;
       if ~isfield(R1{j1},'mpid');R1{j1}.mpid=R1{j1}.il(1)*[1 1];
       else; R1{j1}.mpid(2)=R1{j1}.il(1);
       end
      end
     catch;
       sdtw('_ewt','%s init failed',st1)
     end
     fun=eval(sprintf('@(xyz)%s(xyz,R1{j1})',st1(2:end))); doPost=1;
    else; fun=@(x)lsutil('gen',x,R1{j1}); doPost=0;
    end
    if ~isfield(R1{j1},'distFcn');R1{j1}.distFcn=fun;else;fun=R1{j1}.distFcn;end
    phi{j1}=fun(model.Node(:,5:7));
    if doPost
     if isfield(R1{j1},'LevelList')
      phi{j1}=defLevelList(phi{j1},R1{j1}.LevelList);
     end
    end
   catch;keyboard;
   end
  end
  out=struct('def',horzcat(phi{:}),'DOF',model.Node(:,1)+.98);
  if length(R1)==1&&isfield(R1{1},'distFcn');
   if isfield(R1{1},'mpid');out.mpid=R1{1}.mpid;end
   out.distFcn=R1{1}.distFcn;
   if isfield(R1{1},'contour');out.contour=R1{1}.contour;end
   if isfield(R1{1},'sel');out.sel=R1{1}.sel;end
  else;out.distFcn=@(xyz)dToCell(xyz,R1,Cam);
  end
    
 else
  %% Base functions
  if isstruct(model);x=model.Node(:,5);y=model.Node(:,6);z=model.Node(:,7);
  else; xyz=model;x=xyz(:,1);y=xyz(:,2);z=xyz(:,3);
  end
  if isfield(R1,'idc') % idc : allow definition of center as NodId
   r2=model.Node(model.Node(:,1)==R1.idc,5:7);
   R1.xc=r2(1);R1.yc=r2(2);R1.zc=r2(3);
  end
  if ischar(R1)
   if strcmpi(R1,'x'); phi={x};
   elseif strcmpi(R1,'y'); phi={y};
   elseif strcmpi(R1,'z'); phi={z};
   else; error('%s not yet implemented');
   end
   out=struct('def',horzcat(phi{:}),'DOF',model.Node(:,1)+.98);
  elseif comstr(R1.shape,'rect');
   %% #Rect -2
   %X-dir : -s(x-xd) + c(y-yd) = 0
   %Y-dir :  c(x-xd) + s(y-yd) = 0
   error('Should used DistFcn');%% #Sphere -2
   out=dToRect([x y z],R1);
   
  elseif comstr(R1.shape,'box')
   %% #Box -2
   % Intersection of distances to the six planes of the boxes
   phi1=dToPlane([x y z],[R1.xc R1.yc R1.zc]+R1.lx*R1.nx,-R1.nx);
   phi2=dToPlane([x y z],[R1.xc R1.yc R1.zc]-R1.lx*R1.nx,R1.nx);
   phi3=dToPlane([x y z],[R1.xc R1.yc R1.zc]+R1.ly*R1.ny,-R1.ny);
   phi4=dToPlane([x y z],[R1.xc R1.yc R1.zc]-R1.ly*R1.ny,R1.ny);
   phi5=dToPlane([x y z],[R1.xc R1.yc R1.zc]+R1.lz*R1.nz,-R1.nz);
   phi6=dToPlane([x y z],[R1.xc R1.yc R1.zc]-R1.lz*R1.nz,R1.nz);
   out=min([phi1 phi2 phi3 phi4 phi5 phi6],[],2);
   
  elseif comstr(R1.shape,'circ'); out=dToCirc([x y z],R1);
   %% #Circ -2
   
  elseif comstr(R1.shape,'sphere');
   error('Should used DistFcn');%% #Sphere -2
   
  elseif strcmpi(R1.shape,'cyl');
   %% #Cyl -2
   % plane: nx*(x-xp) + ny*(y-yp) + nz*(z-zp)
   % axis: |CM x n| - R
   xx=x-R1.xc;yy=y-R1.yc;zz=z-R1.zc;
   xc=[R1.xc;R1.yc;R1.zc]*[1 1]+[R1.nx;R1.ny;R1.nz]*[R1.z0 R1.z1];
   
   phi1= R1.rc - sqrt( (yy.*R1.nz-zz.*R1.ny).^2 + (zz.*R1.nx-xx.*R1.nz).^2 ...
    + (xx.*R1.ny-yy.*R1.nx).^2 ) ; %lateral surface
   phi2=R1.nx*(x-xc(1,1))+R1.ny*(y-xc(2,1))+R1.nz*(z-xc(3,1));%bottom (z0)
   phi3=-(R1.nx*(x-xc(1,2))+R1.ny*(y-xc(2,2))+R1.nz*(z-xc(3,2)));%top (z1)
   out=min([phi1,phi2,phi3],[],2);
   
   if isfield(R1,'toAxis')&&R1.toAxis % keep score as distance to axis signed by being inside cyl
    out=abs(phi1).*sign(out);
   end
   
  elseif strcmpi(R1.shape,'cyla');
   %% #CylA: angular position (degrees) with respect to a cylinder -2
   x=x-R1.xc;y=y-R1.yc;z=z-R1.zc;
   p=sp_util('basis',[R1.nx R1.ny R1.nz],[0 0 0]);
   r2=[x y z]*p(:,2:3);
   out=atan2(r2(:,1),r2(:,2))*180/pi;
   
  elseif strncmpi(R1.shape,'toseg',5);
   %% #ToSeg: distance to segment -2
   error('Obsolete replaced by cell array call');
   
  elseif strncmpi(R1.shape,'toplane',5);
   %% #ToPlane: distance to plane -2
   out=dToPlane([x y z],[R1.xc R1.yc R1.zc],[R1.nx,R1.ny,R1.nz]);
   if isfield(R1,'lc');out=R1.lc-out;end
   
  elseif strncmpi(R1.shape,'tes',3);
   %% #Tes: distance to 2D or 3D tessellation -2
   out=R1.distInt(x,y,z);
   
  elseif strncmpi(R1.shape,'cnem',4);
   %% #cnem: PIMM CNEM -2
   out=NaN*ones(size(x,1),1);
   if isfield(R1,'box')
    nidin=x>R1.box(1,1) & x<R1.box(2,1) & ...
     y>R1.box(1,2) & y<R1.box(2,2) & z>R1.box(1,3) & z<R1.box(2,3);
   else
    nidin=(1:size(x,1))';
   end
   XYZ=[x(nidin,:),y(nidin,:),z(nidin,:)];
   Interpol=m_cnem3d_interpol(R1.xyz,[],XYZ,0);
   d1=Interpol.interpolate(R1.val);d1(Interpol.In_Out==0)=NaN;
   out(nidin)=d1;
   
   
  elseif isfield(R1,'distFcn')
   out=R1.distFcn([x y z]);
  else;error('gen%s unknown',R1.shape);
  end
  if isfield(R1,'rsc'); out=out*R1.rsc; end % rescaling
  if isfield(R1,'LevelList');out=defLevelList(out,R1.LevelList); end
  
  %if any(~isfinite(out));error('Problem with finite value');end
 end
 
 %% Use the max a single level set
 if ~isempty(strfind(Cam,'-max')) %#ok<STREMP>
  if ~isstruct(out);out=struct('def',out,'DOF',model.Node(:,1)+.98);end
  [out.def,i2]=max(out.def,[],2);
  if nargout>1;out1=i2;end
 end
 if isnumeric(out);out(abs(out)<sp_util('epsl'))=0;end
 
 % ----------------------------------------------------------------------------
 %% #Edge : #EdgeUtilities  edge utilitilies ----------------------------
elseif comstr(Cam,'edge');[CAM,Cam]=comstr(CAM,5);
 
 %% #EdgeCut : genarate information about how edges are cut
 if comstr(Cam,'cut')
  
  model=varargin{carg};carg=carg+1;
  def=varargin{carg};carg=carg+1;
  RO=varargin{carg};carg=carg+1;
  if ~isfield(RO,'conn')||isempty(RO.conn)
   mo1=model; 
   if isfield(RO,'Elt');
    if isnumeric(RO.Elt);mo1.Elt=RO.Elt;
    elseif ischar(RO.Elt);mo1.Elt=feutil(['selelt ' RO.Elt],mo1);
    end
   else;[eltid,model.Elt]=feutil('eltidfix;',model);
   end
   RO.conn=feval(feutilb('@levNodeCon'),[],mo1,'econ');
   RO.conn.edges(any(RO.conn.edges==0,2),:)=[];
  end
  if size(def.def,2)>1
   if ~isfield(RO,'method');def.def=max(def.def,[],2);
   else;  def.def=feval(RO.method,def.def);
   end
  end
  def=feutilb('placeindof',model.Node(:,1)+.98,def);
  
  edges=RO.conn.edges;% With 3 node last is middle, uses nodeNumber
  %if isfield(RO,'Usable');edges=edges(all(ismember(edges,RO.Usable),2),:);end
  
  NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
  c1=sparse(1:size(edges,1),NNode(edges(:,1)),1,size(edges,1),size(def.DOF,1));
  c2=sparse(1:size(edges,1),NNode(edges(:,2)),1,size(edges,1),size(def.DOF,1));
  
  r1=c1*def.def; r2=c2*def.def;

  if ~isempty(strfind(Cam,'all'));i1=1:size(r1,1);else; i1=r1.*r2<=0;end %#ok<STREMP> 
  r1=r1(i1);r2=r2(i1);i1=full(NNode(edges(i1,:)));
  
  % cf=feplot(model,def);fecom('colordata98');set(cf.ga,'clim',[-.1 .1])
  if nargout>1;out1=model;end
  RO.cEGI=find(any(RO.conn.ElNoCon(model.Node(i1(:,1)),:)& ...
   RO.conn.ElNoCon(model.Node(i1(:,2)),:),1))';
  RO.edges=i1;RO.values=[r1 r2];% RO.edges now nodeindex based
  RO.NNode=NNode;
  
  if isfield(RO,'morph') % Move cut nodes below threshold
   r=r2./(r2-r1);
   [r3,i3]=sort(min(r,1-r));
   for j1=1:length(i3)
    i2=i3(j1); r3=r2(i2)./(r2(i2)-r1(i2));
    if r1(i2)>0%&&r3>.5%#ok<BDSCI> % if first node > move to edge
     model.Node(i1(i2,1),5:7)= r3*model.Node(i1(i2,1),5:7)+ ...
      (1-r3)*model.Node(i1(i2,2),5:7);
     r1(i2,1)=0;
    elseif r2(i2)>0%&&r3<.5 %#ok<BDSCI> % move second node to edge
     model.Node(i1(i2,2),5:7)= r3*model.Node(i1(i2,1),5:7)+ ...
      (1-r3)*model.Node(i1(i2,2),5:7);r2(i2,1)=0;
    else;
     disp([i1(i2,:) r1(i2) r2(i2)])
    end
   end
   cf=feplot(4);feplot(model);
  elseif isfield(RO,'tolE')
   RO.NNode=NNode;
   [model,RO,def]=tolMoveNode(model,RO,def);% def may be wrong
   
   if nargout>1;out1=model;end
  elseif isfield(RO,'tolX');error('Obsolete');
  end
  if ~isfield(RO,'newTol');RO.newTol=.01;end
  
  %no need to be cut
  if ~(any(def.def>0) && any(def.def<0))
   RO.cEGI=[];out=RO;out1=model;out2=def;return;
  end
  
  % List of elements with cuts. ElNoCon ususe node numbers and eltindex
  if isfield(RO,'ed2elt')
   %% #ed2elt define EdgeIndToEltInd and EltIndToEdgeInd (in RO.edges, and model.Elt) -3
   cElEd=(RO.conn.ElNoCon(model.Node(RO.edges(:,1),1),RO.cEGI)+ ...
    RO.conn.ElNoCon(model.Node(RO.edges(:,2)),RO.cEGI))';
   [I,J]=ind2sub(size(cElEd),find(cElEd==2));
   s1=accumarray(J,1);RO.ed2elt=mat2cell(I,s1,1);
   [a,ia]=sort(I);s1=accumarray(a,1);RO.elt2ed=mat2cell(J(ia),s1,1);
  end
  
  if isfield(RO,'cf')
   if isa(RO.cf,'sdth');cf=RO.cf; else;cf=feplot(RO.cf);end
   cf.sel={['eltind' sprintf('%i ',RO.cEGI)],'colorfacew'};
  end
  r1=stack_get(model,'','EltOrient','get');
  if ~isempty(r1);RO.EltOrient=r1;end
  RO.nextId=max(RO.conn.EltId)+1;
  out=RO;      if nargout>2;out2=def;end
  
  
 elseif comstr(Cam,'sellevellines')
  %% #EdgeSelLevelLines
  if ~isempty(evt)
   [CAM,Cam]=comstr(evt.CAM,1);
   sel=obj;RO=varargin{carg};carg=carg+1;
   if comstr(Cam,'usednode')% From fesuper
    out=unique([sel.Node(sel.Node>0);RO.Edges(:)]);
   elseif comstr(Cam,'serestitfix')
    %% #SelEdgeSeRestitFix From fesuper -2
    data=sel.cna{1};
    % c for dependent DOFs
    % clean each restit
    RO.DOF=reshape(data.DOF,[],3);
    if ~isfield(RO,'Nend');
     RO.Nend=max(1e5,10^ceil(log10(max(ceil(data.DOF)))));
     r2=sprintf('%i',max(ceil(data.DOF)));
     RO.Nend=max(1e5,(str2double(r2(1))+1)*10^(length(r2)-1));
    end
    i1=sort(fe_c(data.DOF,unique(RO.Edges),'ind'));
    ind=sort(fe_c(RO.DOF(:,1),unique(RO.Edges),'ind',2));%KeptNodes
    cind=setdiff(1:size(RO.DOF,1),ind); % Edge Nodes
    out=data;
    out.DOF=[RO.DOF(ind,:);
     reshape(feutil('getdof',(1:3)'/100, ...
     RO.Nend+(1:size(RO.cna{1},1)/3)'),[],3)];
    i4=(1:length(ind))';
    RO.nind=sparse([ind;ind+size(RO.DOF,1);ind+size(RO.DOF,1)*2],1, ...
     [i4;i4+size(out.DOF,1);i4+size(out.DOF,1)*2],length(data.DOF),1);
    i4=(length(ind)+1:size(out.DOF,1))';
    RO.ncind=[i4;i4+size(out.DOF,1);i4+size(out.DOF,1)*2];
    % Check KeptInd
    if isfield(data,'KeptInd')
     i3=int32([cind;cind+size(RO.DOF,1);cind+2*size(RO.DOF,1)]);
     for j1=1:size(data.KeptInd,1)
      if ~isempty(intersect(out.KeptInd{j1,1},i3(:)));
       warning('KeptInd mix not implemented');
       %fecom('shownodemark',data.DOF(data.KeptInd{j1}))
      end
     end
    end
    % Adjust Restit
    for j1=1:size(data.Restit,1)
     i4=int32(full(RO.nind(data.Restit{j1,1})));
     if all(i4);out.Restit{j1,1}=i4;
     elseif ~any(i4)
      out.Restit{j1,1}=int32(RO.ncind);
      out.Restit{j1,3}=(RO.cna{1}+RO.cna{2})*out.Restit{j1,3};
     else; error('Not implemented')
      %fecom('shownodemark',data.DOF(data.Restit{j1,1}))
     end
     if isfield(out,'cGL')&&isfield(out.cGL{j1},'vert0')
      out.cGL{j1}.vert0=reshape(sel.vert0(out.Restit{j1,1}),[],3);
     end
    end
    out.DOF=out.DOF(:); out=feutil('rmfield',out,'KeptInd');
    sel.cna{1}=out; out=sel;
   elseif comstr(Cam,'cna')
    %% Cna build from get_def_vertices
    % feval(sel.SelFcn{1},sel,struct('CAM','cna'),sel.SelFcn{2:end},def);
    
    def=varargin{carg}; carg=carg+1;
    %[sel.Node,un1,i2]=unique(sel.Edges(:,1));
    r2=sel;r2.Node=fix(RO.DOF(1:size(RO.DOF,1)/3));
    r2=fe_def('iselcna',[],def,r2,'new');
    sel.cna{1}=(RO.cna{1}+RO.cna{2})*r2.cna{1};
    %  r=sel.EdgeR; r=diag(sparse([r;r;r]));
    %  sel=fe_def('iselcna',[],def,sel,'new');
    %  sel.cna{1}=(r* ...
    %      sparse(1:length(i2)*3, ...
    %      [i2;i2+length(sel.Node);i2+2*length(sel.Node)],1))* ...
    %        sel.cna{1};
    %  r2=sel;[r2.Node,un1,i3]=unique(sel.Edges(:,2));
    %  r2=fe_def('iselcna',[],def,r2,'new');
    %  sel.cna{1}=sel.cna{1}+((speye(size(r))-r)* ...
    %    sparse(1:length(i3)*3, ...
    %      [i3;i3+length(r2.Node);i3+2*length(r2.Node)],1))*r2.cna{1};
    %  sel.Node=sel.Edges(:,1);
    out=sel;
   else;error('SelLevelLines%s',CAM);
   end
  elseif strcmpi(Cam(end+[-2:0]),'cna')
   dbstack;keyboard
   
  else
  %% basic generation of level lines
   model=varargin{carg}; carg=carg+1;cf=[];
   if isa(model,'sdth');cf=model;model=cf.mdl.GetData;
   elseif isa(model,'v_handle');model=model.GetData;
   end
   
   RO=varargin{carg}; carg=carg+1;

   Range=struct('val',[],'lab',{{'gen','LevelList','ProId'}}, ...
          'param',struct('gen',struct('type','pop','AutoAdd',1), ...
          'sel',struct('type','pop','AutoAdd',1,'choices',{{''}})));
   if isfield(RO,'subs'); evt=RO;
   elseif isfield(RO,'gen')% sdtweb t_feplot LevelSet
    if ~isfield(RO,'LevelList');RO.LevelList=0;end
    evt=struct('type','{}','subs',{{{'g',RO.gen},RO.LevelList}});
   end
   %% Interpret URN
   r1=stack_get(model,'info','SelLevelLines','get');
   if isfield(r1,'Elt');mpid=feutil('mpid',r1.Elt);else; mpid=[];end
   for j1=1:length(evt)
       st1=evt(j1).subs{1};
       if ischar(st1)&&any(st1=='=')% {y=val,x>v2,ByProId} 
         [Range,val]=sdtm.range('popMerge',Range,'gen',regexprep(st1,'[\W]*=.*','')); 
         val(:,2)=str2double(regexprep(st1,'.*=',''));
         [Range,val(:,4)]=sdtm.range('popMerge',Range,'sel',evt(j1).subs{2}); 
         r2=val;Range.lab{4}='sel';
       else % {y,val,ByProId}
        [Range,val]=sdtm.range('popMerge',Range,'gen',st1);
        r2=evt(j1).subs{2}; if ischar(r2);r2=comstr(r2,-1);end;r2=r2(:);
        r2(:,2)=val;r2=r2(:,[2 1]); 
       end
       st2=''; if length(evt(j1).subs)>2;st2=evt(j1).subs{3};end
       if isempty(st2);r2(:,3)=0; 
       elseif strcmpi(st2,'byproid')
         %% {x,level,ByProId}
         i2=unique(mpid(:,2));i2(1,:)=[];i3=reshape(repmat(i2(:)',size(r2,1),1),[],1);
         r2=repmat(r2,length(i2),1);r2(:,3)=i3;
       else
         [Range,val]=sdtm.range('popMerge',Range,'ScanMode',evt(j1).subs{3});
         r2(:,3)=val;
       end
       Range.val(end+(1:size(r2,1)),1:size(r2,2))=r2;
   end
   if length(Range.lab)==4;Range.val(Range.val(:,4)==0,4)=1;end% allow mix sel non sel
   if isfield(RO,'subs');RO=struct('Range',Range);
   else; RO.Range=Range;end
   %%
   r1=stack_get(model,'info','SelLevelLines','get');
   if ~isempty(r1);RO=sdth.sfield('addmissing',RO,r1);end
   if isfield(RO,'Elt');% Possibly provide sub model
    if ischar(RO.Elt);
     RO.Elt=feutil(['selelt',RO.Elt],model);
     if isfield(RO,'UsableNodes');
        [n2,RO.Elt]=feutil('quad2lin',model.Node,RO.Elt);
        RO.Elt=feutil('selelt innode',model.Node,RO.Elt,RO.UsableNodes);
     end
    end
   elseif isfield(RO,'Sel')&&ischar(RO.Sel)
    RO.Elt=feutil(['selelt',RO.Sel],model);
    [eltid,RO.Elt]=feutil('eltidfix;',model.Node,RO.Elt);
   end
   if isfield(RO,'gen')
    d1=lsutil('gen',model,RO.gen);
    r2=d1.def;%r2=r2*diag(1./max(abs(r2)));
    def.def=prod(r2,2);def.DOF=d1.DOF;
   elseif isfield(RO,'def'); def=RO.def;
   elseif carg<=nargin; def=RO;RO=varargin{carg};carg=carg+1;
   else; def=struct('def',model.Node(:,5),'DOF',model.Node(:,1)+.01);
   end
   if ~isfield(RO,'oedges')||(isfield(RO,'reset')&&RO.reset) % First time init
    RO=sdth.sfield('rmfield',RO,'reset');
    RO=lsutil('edgecutall',model,def,RO);
    if isempty(RO.edges); error('No edge cut');end
    [r3,i2]=unique(sort(RO.edges,2),'rows');
    RO.edges=RO.edges(i2,:);RO.values=RO.values(i2,:);%Index based
    RO.oedges=RO.edges; RO.ovalues=RO.values;
   elseif ~isempty(def) % Second time update values
    r1=[];r1(RO.NNode(fix(def.DOF)),1)=def.def(:,1); 
    RO.ovalues=r1(RO.oedges(:,1:2));
   end   
  if ~isempty(def);RO.def=def; end
  if ~isempty(cf);
   if ~isfield(RO,'SelFcn');RO.SelFcn={@lsutil,'edgeSelLevelLines'};end
   stack_set(cf,'info','SelLevelLines',RO);
   if ~isfield(RO,'LevelList');return;end
  end
  out=[]; RO.Node=model.Node; 
  if isfield(RO,'ScanMode')
   if strncmpi(RO.ScanMode,'bypro',5)
     i1=setdiff(unique(mpid),0);  li=cell(length(i1),1);
     for j1=1:length(li); li{j1}=find(mpid==i1(j1));end
   else
     mpid=feutil('proid',RO.Elt);
     i1=comstr(RO.ScanMode(6:end),-1);li=cell(length(i1),1);
     for j1=1:length(li); li{j1}=find(mpid==i1(j1));end
   end
  elseif isfield(RO,'cEGI')&&~isempty(RO.cEGI)
     li={RO.cEGI};
  end
  if ~isfield(RO,'Nend');
     RO.Nend=max(1e5,10^ceil(log10(max(model.Node(:,1)))));
  end

  for jPar=1:size(RO.Range.val,1)
    evt=fe_range('valCell',RO.Range,jPar,struct('Table',2));
    if isfield(def,'gen')&&isequal(def.gen,evt.gen)
    else; 
     if iscell(evt.gen);d1=lsutil('gen',model,evt.gen);
     elseif isstruct(evt.gen);d1=lsutil('gen',model,{evt.gen});
     else; d1=lsutil('gen',model,evt.gen);
     end
     r2=d1.def;%r2=r2*diag(1./max(abs(r2)));
     def.def=prod(r2,2);def.DOF=d1.DOF;def.gen=evt.gen;RO.def=def;
    end
    if isempty(evt.ProId)||evt.ProId==0; RO.cEGI=find(isfinite(RO.Elt(:,1)));
    else; RO.cEGI=find(mpid(:,2)==evt.ProId);
    end
    if isfield(evt,'sel')&&~isempty(evt.sel)
      % Possibly restrict cut {z==.95,ProId2,ByProId}
      i1=feutil(['findelt' evt.sel],RO);
      RO.cEGI=intersect(RO.cEGI,i1); if isempty(RO.cEGI);continue;end
    end
   sel=isoContour(RO,evt);
   out=sdtm.feutil.MergeSel('merge',{out,sel});
  end% Range
  if isempty(out); error('Nothing selected');end
  out.Node=out.StressObs.EdgeN(:,1);
  i1=out.StressObs.r>.5; out.Node(i1)=out.StressObs.EdgeN(i1,2);
  if isfield(out,'mdl') % Renumber based on closest edge node
   out.mdl=feutil('renumber -noOri;',out.mdl,out.Node);out.mdl.Node(:,2)=999;
  end
  
  end 
  out.StressObs.InitFcn={@lsutil,'EdgeCna'};out.cna=[];
 elseif comstr(Cam,'self2')
  %% #EdgeSelF2 edit selection to include level lines
  if ~isempty(obj);% Called from feutil getpatch (allow renumber)
   sel=obj;RO=varargin{carg};carg=carg+1;
  else; RO=varargin{carg};carg=carg+1;
   sel=varargin{carg};carg=carg+1;
  end
  model=RO.model;
  if isempty(RO.elt);
   warning('Expecting non empty selection');out=sel;return;
  end
  if isfield(RO,'NNode');NNode=RO.NNode;
  else; NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
  end
  RO.ielt=reshape(full(NNode(RO.elt(:,1:4))),[],4);
  
  [RO.iedges,i2,i3]=unique(RO.elt(:,5:6));%localize unique edges
  RO.elt(:,5:6)=reshape(i3,[],2); % Allow renumber
  r2=reshape(RO.elt(:,8:9),[],1);r=r2(i2,:);
  r2=[RO.ielt(:,1:2);RO.ielt(:,3:4)];RO.iedges=r2(i2,:);
  if ~isfield(sel,'Node');sel.Node=[];end
  r2=struct('iEdge0',size(sel.Node,1),'Edges',model.Node(RO.iedges), ...
   'r',r,'Elt',RO.elt);
  r4=diag(sparse([r;r;r])); r=diag(sparse(r));
  RO.ioff=size(sel.vert0,1);
  sel.vert0=[sel.vert0;r*model.Node(RO.iedges(:,1),5:7)+ ...
   (speye(size(r))-r)*model.Node(RO.iedges(:,2),5:7)];
  %line(sel.vert0(:,1),sel.vert0(:,2),sel.vert0(:,3),'linestyle','none','marker','o');
  sel.f2=[sel.f2;RO.ioff+RO.elt(:,5:6)];
  sel.Node=[sel.Node;-r2.Edges(:,1)];
  
  
  i1=unique(RO.elt(:,1:4));
  def=struct('DOF',feutil('getdof',(1:3)'/100,i1),'def',zeros(length(i1)*3,1));
  r3=sel;r3.Node=model.Node(RO.iedges(:,1));r3=fe_def('iselcna',[],def,r3,'new');
  r2.DOF=def.DOF;r2.cna={r4*r3.cna{1}};
  r3=sel;r3.Node=model.Node(RO.iedges(:,2));r3=fe_def('iselcna',[],def,r3,'new');
  r2.cna{2}=(speye(size(r4))-r4)*r3.cna{1};
  r2.ivert=RO.ioff+(1:size(RO.iedges,1));
  %r=sel.EdgeR; r=diag(sparse([r;r;r]));
  %sel=fe_def('iselcna',[],def,sel,'new');
  %sel.cna{1}=(r* ...
  %    sparse(1:length(i2)*3, ...
  %    [i2;i2+length(sel.Node);i2+2*length(sel.Node)],1))* ...
  %      sel.cna{1};
  %r2=sel;[r2.Node,un1,i3]=unique(sel.Edges(:,2));
  %r2=fe_def('iselcna',[],def,r2,'new');
  %sel.cna{1}=sel.cna{1}+((speye(size(r))-r)* ...
  %  sparse(1:length(i3)*3, ...
  %    [i3;i3+length(r2.Node);i3+2*length(r2.Node)],1))*r2.cna{1};
  if isfield(evt,'jPost');
   sel.SelFcn(evt.jPost,1:3)={@lsutil,'edgeSelLevelLines',r2};
  else;sel.SelFcn={@lsutil,'edgeSelLevelLines',r2};
  end
  sel.f2Prop={'FaceColor','none','EdgeColor','k','marker','none','linewidth',2};
  sel.if2=zeros(size(sel.f2,1),1,'int32');
  sel.opt(1,2)=sdtdef(1); % stamp
  out=sel;
  
  %%
 elseif comstr(Cam,'gensel')
  %% #EdgeGenSel save colormap as level set
  % d1=lsutil('edgeGenSel',cf,RO)
  cf=varargin{carg};carg=carg+1;
  RO=varargin{carg};carg=carg+1;
  if ~isfield(RO,'ob');RO.ob=cf.o(1,1);end
  
  sel=cf.sel;
  r1=get(RO.ob,'facevertexcdata');
  def=struct('DOF',sel.Node(:,1)+.98,'def',r1);
  if isfield(RO,'LevelList')&&ischar(RO.LevelList)&&strcmpi(RO.LevelList,'auto')
   r1=get(cf.ga,'clim');r1=linspace(r1(1),r1(2),10);r1([1 end])=[];
   RO.LevelList=r1;
  end
  def.def=defLevelList(def.def,RO.LevelList(:));
  
  def.DOF=(1:size(sel.vert0,1))'+.98;
  mo2=struct('Node',[(1:size(sel.vert0,1))'*[1 0 0 0] sel.vert0],'Elt', ...
   feutil('addelt','quad4',sel.fs));
  sel2=lsutil('edgeSelLevelLines',mo2,struct('def',def));
  r2=def;r2.DOF=fix(def.DOF)+.01;
  sel2=feval(sel2.SelFcn{1},sel2,struct('CAM','cna'),sel2.SelFcn{2:end},r2);
  r3=reshape(sel2.cna{1}*get(RO.go,'facevertexcdata'),[],3)*[1;0;0];
  sel2=sdsetprop(sel2,'f2Prop','facevertexcdata',r3, ...
   'linewidth',2,'edgecolor','interp');
  
  if isfield(RO,'cf')
   %% possibly generate new figure with appropriate level set lines
   
   figure(RO.cf);RO.ga=get(RO.cf,'currentaxes');
   if isempty(RO.ga); RO.ga=axes('parent',RO.cf); else;cla(RO.ga);end
   set(RO.ga,'dataaspectratio',[1 1 1]);sel2.cna={};
   set(RO.ga,'clim',get(cf.ga,'clim'))
   set(ancestor(RO.ga,'figure'),'colormap',get(cf.opt(1),'colormap'));
   
   [sel2,uaob]=feutilg('initPatch_object',sel2,RO.ga,zeros(1,3));
   s2=sel; s2.f2=[];
   s2=sdsetprop(s2,'fsProp','facecolor','w','facealpha',.1,'edgecolor','k','edgealpha',.05);
   [s2,uaob]=feutilg('initPatch_object',s2,RO.ga,zeros(1,3));
   colorbar
   
  end
  
  out=sel2;
 elseif comstr(Cam,'cna') 
  %% #EdgeCna /Cta : intialize observation
   def=[];sel=[];eval(iigui({'def','sel'},'GetInCaller'))
   r1=sel.StressObs; i1=reshape(r1.EdgeN',[],1);
   if isequal(varargin{end},'cta')
     % display of sensors assuming dir predefined and sensors at end
     % sdtweb fegui SaveSel.AddTestNode
     if ~isequal(def.tdof,r1.tdof); error('Not expected');end
     out=sparse(size(sel.dir,1)+(-size(def.tdof,1)+1:0),1:size(def.tdof,1),1);

   else
    % cna : display of FEM shape
    if isfield(def,'TR');DOF=def.TR.DOF;
    elseif ~isfield(def,'DOF');out=[];return;
    else; DOF=def.DOF;
    end
    r2=fe_c(DOF(:,1),[i1+.01;i1+.02;i1+.03],'place');
    i2=(1:length(r1.r))'; i3=i2(end);
    r2=sparse([i2 i2 i2+i3 i2+i3 i2+2*i3 i2+2*i3], ...
        [i2*2-1,i2*2, i2*2-1+2*i3,i2*2+2*i3  i2*2-1+4*i3,i2*2+4*i3], ...
        [1-r1.r r1.r 1-r1.r r1.r 1-r1.r r1.r],i3*3,size(r2,1))*r2;
    if isfield(def,'TR');r2=r2*def.TR.def;end
    out=r2;
   end

 else; error('Edge%s',CAM);
 end
 
 %% #eltset: build element sets based on levelsets results
elseif comstr(Cam,'eltset'); [CAM,Cam]=comstr(CAM,7);
 %  model=lsutil('eltset',model,struct('list',{li},'CritFcn',@(x)x>=0, ...
 %     'sel','innode'));
 model=varargin{carg}; carg=carg+1;
 RO=varargin{carg}; carg=carg+1;
 if iscell(RO); RO=struct('list',{RO}); end
 [RO,st,CAM]=cingui('paramedit -DoClean',[ ...
  'CritFcn("@(x)x>=0"#%s#"criterion function on level set")' ...
  'sel("innode"#%s#"element selection for node list")' ...
  'postSel(#%s#"additional selection after lsset")' ...
  'name("_lsset"#%s#"meta-set name")' ...
  ],{RO,CAM}); Cam=lower(CAM);
 if ischar(RO.CritFcn)&&isequal(RO.CritFcn(1),'@');
  RO.CritFcn=eval(RO.CritFcn);
 end
 d1=lsutil('gen',model,RO.list);
 if size(RO.list,2)>1 % table cell input
  RO.lab=RO.list(1,:); RO.list=RO.list(2:end,:);
 end
 RO.set=cell(size(d1.def,2),2);
 for j1=1:size(d1.def,2)
  R1=RO.list(j1,:);
  if isstruct(R1{1}); R1=R1{1};
  else; R1=cell2struct(R1(:),RO.lab(:));
  end
  if isfield(R1,'label'); lab=R1.label;
  else; lab=sprintf('%s_%i',R1.shape,j1);
  end
  stn=sprintf('_nset_%s_%i',RO.name,j1);
  model=feutil('AddSetNodeId',model,stn,...
   unique(fix(d1.DOF(feval(RO.CritFcn,d1.def(:,j1)))))  );
  if isempty(RO.postSel); stp='';
  else; stp=RO.postSel; if ~ismember(stp(1),{'&','|'}); stp=['&' stp]; end
  end
  RO.set(j1,:)={ sprintf('%s{setname %s}%s',RO.sel,stn,stp),lab};
 end
 model=feutil('addseteltid-append',model,RO.name,RO.set(:,1),RO.set(:,2));
 out=model;
 
 % ----------------------------------------------------------------------------
 %% #surf : #SurfaceUtilities  surface utilitilies ----------------------------
elseif comstr(Cam,'surf');[CAM,Cam]=comstr(CAM,5);
 
 %% #SurfStreamLine (obsolete ?)
 if comstr(Cam,'stream')
  sdtw('_ewt','report eb');
  model=varargin{carg};carg=carg+1;
  MAP=varargin{carg};carg=carg+1; % MapAtNode
  if carg<=nargin; RO=varargin{carg};carg=carg+1;else; RO=struct;end
  if ~isfield(RO,'start'); RO.start=MAP.ID(1);% Start at nodes
  end
  if ~isfield(RO,'Out');RO.Out='line';end
  if ~isfield(RO,'dl');RO.dl=.01;end
  
  RO=OEdgeElt(RO,model);
  if isempty(MAP)&&isfield(RO,'dir') % Fixed orientation
   MAP=struct('ID',model.Node(:,1), ...
    'normal',repmat(RO.dir(:)',size(model.Node,1),1));
  end
  if isfield(MAP,'ID')
   RO.MapInd=sparse(MAP.ID,1,1:size(MAP.ID,1));
  else; % InfoAtNode MapjElt
   i1=feutil('eltid',model);
   [i1,i2]=ismember(i1,MAP.EltId);
   RO.MapjElt=sparse(find(i1),1,i2(i1));
  end
  
  RO.Lines=cell(size(RO.start,1),1);  RO.cut=[];
  for j3=1:length(RO.Lines)
   %% Loop on lines
   Ls=zeros(7,1000); RO.iLine=[-2 zeros(1)];
   %sp_util('setinput',Ls,1:3,RO.iLine,'Ls');
   
   % SurfStreamLine
   %  1. Cur -> Ref : get origin, theta_i at corners
   %  2. RefLine : pn+1 = p_n + dl * v1(pn) until exit
   %  3. x(physical), next-element
   
   if size(RO.start,2)==7 % Previous seed
    n1=RO.start(j3,1:3);jElt=RO.start(j3,7);n1p=n1;
   else
    in1=RO.NNode(RO.start(j3));
    n1=model.Node(in1,5:7); n1p=n1;
    [un1,un1,jElt]=find(RO.OEdgeElt(:,in1),1);
   end
   ElemP='quad4'; inode=1:4;
   EC=integrules('quad4');
   EC.jdet=ones(1); EC.bas=zeros(9,1);EC.Nw=1;EC.NDN=zeros(4,3);EC.J=zeros(4,1);
   EC.nodeEt=int32(comstr({'x','y','z','Id','v1x','v1y','v1z'},-32));
   sdtw('_ewt','integrules quad4 inline');
   RO.N=@(w)(1+w(:,1)*EC.xi(:,1)').*(1+w(:,2)*EC.xi(:,2)')/4;
   RO.Nr=@(w)[EC.xi(1:4,ones(1,size(w,1))*1)'.*(1+w(:,2)*EC.xi(:,2)')/4 ];
   RO.Ns=@(w)[EC.xi(1:4,ones(1,size(w,1))*2)'.*(1+w(:,1)*EC.xi(:,1)')/4 ];
   
   while 1
    %% Seek until end
    i2=RO.NNode(model.Elt(jElt,inode));
    n2=model.Node(i2,:);
    if isfield(MAP,'ID') % Using fieldatnode
     EC.nodeE=[n2(:,[5:7 1]) MAP.normal(RO.MapInd(n2(:,1)),:)];
    else; % Using InfoAtNode
     r3=MAP.data(:,MAP.NodePos(:,RO.MapjElt(jElt)))';
     if size(r3,1)==1;r3=repmat(r3,size(n2,1),1);end
     EC.nodeE=[n2(:,[5:7 1]) r3];
    end
    p = FindRST(ElemP,n2(:,5:7)-repmat(n1,length(inode),1),1e-4)'; % first point within element
    while 1==1
     EC.N=RO.N(p); EC.Nr=RO.Nr(p); EC.Ns=RO.Ns(p); EC.w=[p 0 0];
     of_mk('buildndn',23,EC);n1=EC.N*EC.nodeE(:,1:3);
     sp_util('setinput',Ls,[n1 norm(n1-n1p) p jElt],RO.iLine,'Ls');n1p=n1;
     if any(abs(p)>1); break;end
     d=reshape(EC.J,2,2)';  p=p+d(1,:)*RO.dl;
     % r1=diff(RO.N([p;p2])*n2(:,5:7))/RO.dl
    end
    [un1,i1]=max(abs(p));
    if i1==1; if p(i1)>=1; i2=2;else; i2=4; end % edge 2
    else; if p(2)>=1; i2=3; else;i2=1;end
    end
    if strcmpi(RO.Out,'morph') % Place exit node on intersection
     if abs(p(i1))>1% should actually go back to edge
      d=d(1,:);d=d/d(i1);p=p-d*(p(i1)-sign(p(i1)));
     end
     %p(p>1)=1;p(p<-1)=-1; n1=RO.N(p)*EC.nodeE(:,1:3);
     RO.iLine(2)=RO.iLine(2)-7;
     sp_util('setinput',Ls,[n1 norm(n1-n1p) p jElt],RO.iLine,'Ls');n1p=n1;
     i1=[1 2 4 3];i1=i1(round((p+1)/2)*[1;2]+1);i1=model.Elt(jElt,i1);
     RO.cut(end+1)=i1;
     model.Node(RO.NNode(i1),5:7)=n1;
    elseif strcmpi(RO.Out,'seed') % seeds along a line
     j2=fix(sum(Ls(4,:))/RO.lc)+1;
     if ~isfield(RO,'OutInd'); RO.OutInd=0;
     elseif j2>length(RO.OutInd); RO.OutInd(j2)=RO.iLine(2);
     end
    end
    i1=quad4('edge');i1=i1(i2,[2 1]);i1=RO.NNode(model.Elt(jElt,i1));
    jElt=RO.OEdgeElt(i1(1),i1(2));
    if jElt==0; break;end
   end
   
   % Now do a line within current element
   Ls=Ls(:,1:RO.iLine(2)/size(Ls,1))';Ls(:,4)=cumsum(Ls(:,4));
   RO.Lines{j3}=Ls;
  end
  if strcmpi(RO.Out,'morph') % Place exit node on intersection
   out=model; out1=RO;
   if isfield(RO,'ga');cf=feplot(ancestor(RO.ga,'figure'));feplot(cf,model);end
  elseif strcmpi(RO.Out,'seed') % equally spaced seeds
   i3=(RO.OutInd)/size(Ls,2)+1; if i3(end)>size(Ls,1);i3(end)=i3(end-1);end
   Ls=Ls(i3,:); out=Ls;
  elseif strcmpi(RO.Out,'line2sided') % equally spaced seeds
   if isfield(MAP,'ID');MAP.normal=-MAP.normal;
   else;MAP.data=-MAP.data;
   end
   RB=RO;RB.Out='line'; RB=feutil('rmfield',RB,'ga');
   RB=feutilb('surfstream',model,MAP,RB);
   for j1=1:length(RB);
    RB{j1}=[flipud(RB{j1});RO.Lines{j1}];
    l1=sqrt(sum((RB{j1}(2:end,1:3)-RB{j1}(1:end-1,1:3)).^2,2));
    RB{j1}(:,4)=cumsum([0;l1]);
   end
   RO.Lines=RB;out=RB;
  else;out=RO.Lines;
  end
  
  if isfield(RO,'ga')
   %figure(2);delete(findall(ga,'type','line'))
   if ~isfield(RO,'LineProp');
    RO.LineProp={'color','k','tag','stream','linewidth',2};
   end
   for j1=1:length(out)
    Ls=RO.Lines{j1}; if isempty(Ls);continue;end
    h=line(Ls(:,1),Ls(:,2),Ls(:,3),RO.LineProp{:});
   end
  end
  
  
  %% #SurfFromPoly : use geodesic approximation to draw on surface
 elseif comstr(Cam,'frompoly')
  
  sdtw('_ewt','report EB for transition to cutLineSurf'); 
  model=varargin{carg};carg=carg+1;
  RO=varargin{carg};carg=carg+1;
  if ~isfield(RO,'tol');RO.tol=1e-4;end
  if ~isfield(RO,'dl');RO.dl=.01;end
  if ~isfield(RO,'start')
   i1=fe_gmsh('LineLoops',feutil('selelt seledge',RO));i1=i1{1}(:);
   if isfield(RO,'CornerMove') % Place corner nodes at right position
    error('fe_shapeoptim');
    NNode=sparse(model.Node(:,1),1,1:size(model.Node));
    n1=feutil('getnode',RO,i1); if ~isfield(RO,'epsl');RO.epsl=1;end
    i1=feutilb(sprintf('addnode -nearest epsl%g',RO.epsl),model.Node,n1(:,5:7));
    i2=feutilb('geolinetopo',model,struct('starts',[i1(1:end-1,1) i1(2:end,1)]));
    for j1=1:4 % Interpolate along line
     n2=model.Node(i2{j1},5:7); n2=n2-repmat(n2(1,:),size(n2,1),1);
     d=diff(n2([1 end],:));d=d(:)/norm(d)^2;r=n2*d;
     model.Node(NNode(i2{j1}),5:7)=(1-r)*n1(j1,5:7)+r*n1(j1+1,5:7);
    end
   end
   RO.start=[i1(1:end-1) i1(2:end)];
  end
  %   RA=struct('start',[5 95;95 36;36 23;23 5],'ga',cf.ga,'LineProp',{{'linestyle',':','marker','o'}},'dl',.03);
  
  i1=feutilb('surfStream',model,[],RO);RO.Sl=vertcat(i1{:});
  r1=RO.Sl(:,1:3);
  RO.ze=RO.Sl(:,8:10);RO.xe=diff([RO.Sl(:,1:3);RO.Sl(1,1:3)]);
  i1=all(abs(diff(r1))<RO.tol,2);
  if all(abs(diff(r1([1 end],:)))<RO.tol);i1(end+1)=1;end
  r1(i1,:)=[];
  % sdtweb lsutil dToPoly.init
  RO.ze(i1,:)=[];RO.xe(i1,:)=[];RO.L=sqrt(sum(RO.xe.^2,2));
  RO.xe=diag(sparse(1./RO.L))*RO.xe; RO.dir=cross(RO.xe,RO.ze);%seg dir
  RO.dirn=RO.dir(1:end-1,:)+RO.dir([end 1:end-2],:); % Node dir
  RO.dirn=diag(sparse(1./sqrt(sum(RO.dirn.^2,2))))*RO.dirn;

  RO.DT = delaunayTriangulation(r1(:,1:3));
  RO.distFcn=@(xyz)dToPoly(xyz,RO); RO.shape='Poly';
  out=RO;  out1=model;
  
  %% #SurfRemesh : create a surface line on initial mesh
 elseif comstr(Cam,'remesh')
  
  model=varargin{carg};carg=carg+1;
  mo1=varargin{carg};carg=carg+1;
  RO=varargin{carg};carg=carg+1;

  if isfield(RO,'epsl'); RO.epsl0=sp_util('epsl');sp_util('epsl',RO.epsl);end

  if iscell(mo1);mo1=mo1{1};end;[EGroup,nGroup]=getegroup(mo1.Elt);
  if nGroup>1
   %% With multiple groups do each one sequentially
   for jGroup=1:nGroup
    mo2=mo1;mo2.Elt=mo1.Elt(EGroup(jGroup):EGroup(jGroup+1)-1,:);
    [model,li]=lsutil('surfremesh',model,mo2,RO);
    if jGroup==1;out1=li;end
   end
   out=model;
   return;
  end
  mo1=dToPoly('init',mo1);
  li=struct('shape','Poly','distFcn',@(xyz)dToPoly(xyz,mo1));
  % cf=feplot;lsutil('ViewLs',model,li)
  mo3=lsutil('cut',model,{li},RO);
  RB=li.distFcn('ro');
  if isfield(RB,'MatId');RB.mpid=[RB.MatId RB.ProId];
  elseif isfield(RO,'mpid');RB.mpid=RO.mpid;
  end
  if isfield(RB,'mpid')
      %% Define surface associated with interior
      i1=li.distFcn(mo3.Node(:,5:7)); i1=(isfinite(i1)&i1>=0);i1=mo3.Node(i1,1);
      try;  r2=li.distFcn('RO'); % Allow onSurfCheck for each cut
       if isfield(r2,'onSurf')&&ischar(r2.onSurf) 
        i1=intersect(i1,feutil(r2.onSurf,mo3));
       elseif isfield(RO,'onSurf')&&ischar(RO.onSurf) 
        i1=intersect(i1,feutil(RO.onSurf,mo3));
       end
      end
      elt=feutil('selelt selface & innode',mo3,i1);
      if ~isempty(elt)
       elt=feutil(sprintf('set groupall matid %i proid%i',RB.mpid),elt);
       mo3=feutil('addelt',mo3,elt);
      end
  end
  if isfield(RO,'NoDegen')&&RO.NoDegen
   mo3=feutil('OptimDegen',mo3);
  end
  out=mo3; %lsutil('ViewLs',mo3,li);
  if isfield(RO,'epsl0'); sp_util('epsl',RO.epsl0);end

  out1=li;
 elseif comstr(Cam,'fromrectmesh')
  %% #SurfFromRectMesh
  mo2=varargin{carg};carg=carg+1;
  if carg<=nargin; RO=varargin{carg};carg=carg+1; else;RO=struct;end
  if ~isfield(RO,'ContourSel'); RO.ContourSel='seledge';end
  mo2.Elt=feutil(['selelt', RO.ContourSel],mo2);
  mo2.Node=feutil('getnodegroupall',mo2);
  n1=fe_gmsh('lineloops',mo2.Elt);n1=feutil('getnode',mo2,unique(n1{1}));
  %fecom('shownodemark',mo2.Node(:,5:7),'marker','o')
  
  li=struct('shape','rect','xc',(max(mo2.Node(:,5))+min(mo2.Node(:,5)))/2,...
   'yc',(max(mo2.Node(:,6))+min(mo2.Node(:,6)))/2,...
   'lx',max(mo2.Node(:,5))-min(mo2.Node(:,5)),'ly',max(mo2.Node(:,6))-min(mo2.Node(:,6)),'alpha',0);
  dfun=lsutil('@dToRect');
  li=struct('shape','Fcn','distFcn',@(x)dfun(li,x),'contour',n1);
  out=li;
 %% #SurfStick : project nodes on surface (used by rail19)
 elseif comstr(Cam,'stick')
  mo1=varargin{carg};carg=carg+1;
  if carg<=nargin; RO=varargin{carg};carg=carg+1; else;RO=struct;end
  if ~isfield(RO,'distFcn'); error('Missing distFcn');end
  
  n1=feutil(['getnode' RO.sel],mo1);
  n2=RO.distFcn(struct('stick',n1(:,5:7)));%dToPoly
  NNode=sparse(mo1.Node(:,1),1,1:size(mo1.Node,1));
  mo1.Node(NNode(n1(:,1)),5:7)=n2; 
  out=mo1; 
  
  %% #SurfEnd
 else;error('Surf%s unknown',CAM);
 end
 
elseif comstr(Cam,'3dintersect'); [CAM,Cam]=comstr(CAM,12);
 %% #3DIntersect: recover local surface selections from macro 3D selections
 
 mo1=varargin{carg}; carg=carg+1;
 RO=varargin{carg}; carg=carg+1;
 [CAM,Cam,RO.forceGFrame]=comstr('-forcegframe',[-25 3],CAM,Cam);
 if ~isfield(RO,'tol')||isempty(RO.tol)||~RO.tol
  [CAM,Cam,RO.tol]=comstr('-tol',[-25 2],CAM,Cam);
 end
 
 % initialize level sets
 li={'shape','xc','yc','zc','lx','ly','lz','nx','ny','nz'};
 % Box containing components
 for j1=1:length(RO.sel)
  n1=feutil(sprintf('getnode inelt{%s}',RO.sel{j1}),mo1); r1=mean(n1(:,5:7));
  if RO.forceGFrame; v=eye(3);
  else;  [u,s,v]=svd(bsxfun(@minus,n1(:,5:7),r1),0);
  end
  r2=max(abs([n1(:,5)-r1(1) n1(:,6)-r1(2) n1(:,7)-r1(3)]*(v')));
  li(end+1,:)=...
   {'box',r1(1),r1(2),r1(3),r2(1),r2(2),r2(3),v(1,:),v(2,:),v(3,:)};
 end
 
 % Generate boxes and recover the intersection
 RA=struct('list',{li},'CritFcn',@(x)x>=-RO.tol,'sel','withnode','postSel','eltname ~=SE');
 mo1=stack_rm(mo1,'set','_lsset');
 mo1=lsutil('eltset',mo1,RA);
 r1=stack_get(mo1,'set','_lsset','get');
 % Now define the area for connection as intersection of boxes
 stSel=sprintf('%s_EltSel',RO.name);
 mo1=feutil('addseteltid',mo1,stSel,r1.EltId(all(r1.SConn,2)));
 
 for j1=1:length(RO.sel) % generate restrained selections
  RO.sel{j1}=sprintf('setname %s & %s',stSel,RO.sel{j1});
 end
 
 out=mo1; out1=RO;
 
 
 
 
 
elseif comstr(Cam,'split');[CAM,Cam]=comstr(CAM,6);
 %% #split : refine elements cut by a level set -----------------------------
 model=varargin{carg};carg=carg+1;li=varargin{carg};carg=carg+1;
 if nargin>=carg;RO=varargin{carg};else RO=[];end
 
 while 1
  jGroup=0;
  if isfield(li,'def')&&isfield(li,'distFcn');def=li;
  else; def=lsutil('gen-max',model,li);% lsutil('ViewLs',model,li);
  end
  [RO,model,def]=lsutil('EdgeCut',model,def,RO);
  [EGroup,nGroup]=getegroup(model.Elt);mo2=struct('Node',model.Node,'Elt',[]);
  %i1=RO.edges((RO.r>=RO.newTol|RO.r<=1-RO.newTol)&all(RO.values~=0,2),:);
  %RO.cEGI=find(any(RO.conn.ElNoCon(model.Node(i1(:,1)),:)& ...
  %  RO.conn.ElNoCon(model.Node(i1(:,2)),:),1))';
  
  proid=feutil('mpid',model); proid(:,3:end)=[];
  proid(:,3)=feutil('eltid',model); % allow eltorient propag
  x=model.Node(:,5);y=model.Node(:,6);z=model.Node(:,7);
  icut=ones(nGroup,1);icut(1:jGroup)=0;listE=[];
  while jGroup<nGroup; %
   jGroup=jGroup+1;
   [ElemF,opt,ElemP] = feutil('GetElemF',model.Elt(EGroup(jGroup),:),jGroup);
   cEGI=(EGroup(jGroup)+1:EGroup(jGroup+1)-1)';
   idelt=RO.cEGI(ismember(RO.cEGI,cEGI));
   if isempty(idelt);
    icut(jGroup)=0;continue;
   else
    % XXX cases where added nodes should be coalesced ???
    elt=model.Elt(idelt,:);
    elti=full(RO.NNode(elt(:,1:length(feval(ElemP,'node')))));
    gid=proid(idelt,:);
    if size(elt,1)==1;ndir=1;elti=elti';else ndir=2;end
    lsC=charLs(def,elti);
    %i3=cellfun(@(x)length(setdiff(unique(x),'0'))<=1,lsC);% No need to cut
    i3=cellfun(@(x)~(any(x=='+')&any(x=='-')),lsC);% No need to cut
    idelt(i3)=[];elt(i3,:)=[];elti(i3,:)=[];lsC(i3,:)=[];gid(i3,:)=[];
    if isempty(idelt);icut(jGroup)=0;continue;end
    if strcmp(ElemP,'hexa8')
     if size(elt,1)==1;ndir=1;else ndir=2;end
     [mo2,ind]=newOnLS(x,y,z,elti,ndir,mo2,def);ind=[ind gid];
     newElt=[[elt(:,1:4) ind];[elt(:,[2 6 7 3]) ind];[elt(:,[6 5 8 7]) ind]
      [elt(:,[5 1 4 8]) ind];[elt(:,[5 6 2 1]) ind]; [elt(:,[4 3 7 8]) ind];];
     icut(jGroup)=1;ElemP='pyra5';
    elseif strcmp(ElemP,'pyra5')
     if size(elt,1)==1;ndir=1;else ndir=2;end
     xc=[sum(x(elti(:,1:4)),ndir) sum(y(elti(:,1:4)),ndir) sum(z(elti(:,1:4)),ndir)]/4;
     [mo2.Node,ind]=feutil('AddNode',mo2.Node,xc);ind=mo2.Node(ind,1); %KnownNew
     newElt=[[elt(:,1:2) ind elt(:,5) gid];[elt(:,2:3) ind elt(:,5) gid]
      [elt(:,3:4) ind elt(:,5) gid];[elt(:,[4 1]) ind elt(:,5)] gid];
     ElemP='tetra4';icut(jGroup)=0;
    elseif any(strcmp(ElemP,{'tetra4','tria3','quad4','quadb'}))
     icut(jGroup)=0;continue;
    elseif strcmp(ElemP,'penta6')
     if size(elt,1)==1;ndir=1;else ndir=2;end
     xc=[sum(x(elti(:,1:6)),ndir) sum(y(elti(:,1:6)),ndir) sum(z(elti(:,1:6)),ndir)]/6;
     [mo2.Node,ind]=feutil('AddNode',mo2.Node,xc);ind=mo2.Node(ind,1); %KnownNew
     newElt=[ind elt(:,4:6) gid;elt(:,1:3) ind gid];ElemP='tetra4';
     newElt1=[elt(:,[2 5 6 3]) ind gid;elt(:,[1 4 5 2]) ind gid;elt(:,[4 1 3 6]) ind gid];ElemP1='pyra5';
     icut(jGroup)=1;
    else
     warning('element %s not supported',ElemF);ElemP='';
    end
    [mo2,RO]=safeAddElt(mo2,RO,ElemP,newElt);
    if exist('newElt1','var');[mo2,RO]=safeAddElt(mo2,RO,ElemP1,newElt1);clear newElt1 ElemP1;end
   end
   listE=[listE;idelt];
  end
  model.Node=mo2.Node;
  model.Elt=feutil('removeelt eltind',model,listE);% Remove refined elts
  if isfield(RO,'EltOrient');
   model=stack_set(model,'info','EltOrient',RO.EltOrient);
  end
  jGroup=nnz(~isfinite(model.Elt(:,1)));
  if ~isempty(mo2.Elt);model=feutil('addelt',model,mo2.Elt);end % add new elts
  if ~sum(icut)||isfield(RO,'iter'); break;end
  RO.conn=[];
 end
 [eltid,model.Elt]=feutil('eltidfix',model); % This should not be necessary with safeAddElt
 out=model;
 if nargout<2
 elseif isfield(def,'def')&&isfield(def,'distFcn');out1=def;
 else; out1=lsutil('gen-max',model,li);
 end
 if nargout>2;out2=RO;end
 
 %% #Cut_: surface build ------------------------------------------------------
elseif comstr(Cam,'cut');[CAM,Cam]=comstr(CAM,4);
 if comstr(Cam,'face2tria')
  %% #cutface2tria: force cut surface to be tria based (for tet remesh) -----2
  % discuss with GVDR 
  mo3=varargin{carg}; carg=carg+1;
  sel=varargin{carg}; carg=carg+1;
  
  if iscell(sel) % level sets
   [u1,mpid]=lsutil('mpid',mo3,sel);
   mo3=feutil('addseteltid',mo3,'cuts',sprintf('eltind %s',num2str(find(mpid(:,1)>0)')));
   mo3=feutil('addsetnodeid',mo3,'cutsn','inelt{setname cuts}');
   mo3=feutil('addsetfaceid',mo3,'cutqi','setname cuts & selface & eltname quad');
   mo3=feutil('addsetfaceid',mo3,'cutqo',...
    'setname cuts:exclude & selface & innode{setname cutsn} & eltname quad');
   elt=feutil('selelt setname cutqi | setname cutqo',mo3);
  elseif isstruct(sel) % a set structure is given
   elt=feutil('selelt',mo3,sel);
  elseif ischar(sel) % FindElt string
   elt=feutil(sprintf('selelt %s',sel),mo3);
  elseif isnumeric(sel)
   if isfinite(sel(1))
    if size(sel,2)>1; typ='FaceId'; else; typ='EltId'; end
    elt=feutil('selelt',mo3,struct('type',typ,'data',sel));
   else; elt=sel; % elts directly given
   end
  else; error('incorrect selection format');
  end
  
  %  % had to flip evertything... used
  %  for jj=fieldnames(R1)'
  %   for j1=1:size(R1.(jj{1}).Elt,1)
  %    R1.(jj{1}).Elt{j1,2}=R1.(jj{1}).Elt{j1,2}(:,feval(R1.(jj{1}).Elt{j1,1},'flip'));
  %   end
  %  end
  %  comstr(R1,-30,struct('NoClip',1))
  
  eltid=feutil('eltid;',mo3); EEid=sparse(eltid+1,1,1:length(eltid));
  r1=feval(feutil('@FaceIdFromSelface'),elt,'f2e'); % [EltId, FId1,... ;] list
  % now rule on faces found to apply correct transform
  r3=sum(logical(r1),2)-1; r4=unique(r3); % number of faces per underlying element
  RO={};
  for j3=1:length(r4) % loop on found faces number found
   r6=r1(ismember(r3,r4(j3)),:); r6(:,r6(1,:)==0)=[];
   switch r4(j3)
    case 1 % one quad face, support pyra, penta, hexa
     R1=struct('pyra5',struct('shift',1,'faces',pyra5('faces'),...
      'face',{{[6 1:4],[1/4]}},...
      'Elt',{{'tetra4',[1 6 2 5;  2 6 3 5;  3 6 4 5;  4 6 1 5]}}), ...
      ...{{'tetra4',[1 2 6 5; 2 3 6 5; 3 4 6 5; 4 1 6 5]}}),...
      'penta6',struct('shift',[2],'faces',....penta6('faces'),...
      [1 3 2 2;1 4 6 3;2 5 4 1;4 5 6  6;3 6 5 2],...
      'face',{{[7 1 4 6 3],[1/4]}},'volume',{{[8 1:6],[1/6]}},...
      'Elt',{{'tetra4',[1 7 4 8;  4 7 6 8;  6 7 3 8;  3 7 1 8;  1 3 2 8;  4 6 5 8];
      'pyra5',[1 4 5 2 8;  2 5 6 3 8]}}), ...
      ...{{'tetra4',[1 4 7 8; 4 6 7 8; 6 3 7 8; 3 1 7 8; 1 2 3 8; 4 5 6 8];
      ... 'pyra5',[1 2 5 4 8; 2 3 6 5 8]}}),...
      'hexa8',struct('shift',1,'faces',hexa8('faces'),...
      'face',{{[9 1:4],[1/4]}},'volume',{{[10 1:8],[1/8]}},...
      'Elt',{{'tetra4',[1 9 4 10;  4 9 3 10;  3 9 2 10;  2 9 1 10];
      'pyra5',[1 4 8 5 10;  1 5 6 2 10;  5 8 7 6 10;  2 6 7 3 10;  4 3 7 8 10]}}));
     ...{{'tetra4',[1 4 9 10; 4 3 9 10; 3 2 9 10; 2 1 9 10];
      ... 'pyra5',[1 5 8 4 10; 1 2 6 5 10; 5 6 7 8 10; 2 3 7 6 10; 4 8 7 3 10]}}));
      R1.set=r6; RO{end+1}=R1;
     
    case 2 % two quad faces, support penta, hexa (2cases)
     r7=r6; % second case
     el0=feutil('selelt eltind',mo3,full(EEid(r6(:,1)+1))); el0id=feutil('eltid;',el0);
     [EG,nG]=getegroup(el0);
     for jG=1:nG
      [ElemF,u1,ElemP]= getegroup(el0(EG(jG),:),jG);
      cEGI=EG(jG)+1:EG(jG+1)-1;
      switch ElemP
       case 'penta6'; r7(ismember(r7(:,1),el0id(cEGI)),1)=0; % always first case
       case 'hexa8'
        i2=feval(ElemP,'faces');
        [i3,i4]=ismember(r6(:,1),el0id(cEGI)); i4=i4(i3); i3=find(i3);
        for jElt=1:length(i3)
         n1=i2(r6(i3(jElt),2:end),:);
         if any(sparse(n1(:),1,ones(numel(n1),1))==2); r7(i3(jElt),1)=0;% first case adjacent faces (one shared node)
         else; r6(i3(jElt),1)=0; % second case two opposite faces (no common nodes)
         end
        end
      end
     end % loop nG
     r6(r6(:,1)==0,:)=[]; r7(r7(:,1)==0,:)=[];
     if ~isempty(r6)  % Two faces /adjacent
      R1=struct('penta6',struct('shift',[2 3],'faces',...penta6('faces'),...
       [1 3 2 2;1 4 6 3;2 5 4 1;4 5 6  6;3 6 5 2],...
       'face',{{[7 1 4 6 3; 8 1 2 5 4],1/4}},'volume',{{[9 1:6],1/6}},...
       'Elt',{{'tetra4',[1 7 4 9;  4 7 6 9;  6 7 3 9;  3 7 1 9;  1 8 2 9;
       2 8 5 9;  5 8 4 9;  4 8 1 9;  4 6 5 9;  1 3 2 9];
       'pyra5',[2 5 6 3 9]}}), ...
       ...,{{'tetra4',[1 4 7 9; 4 6 7 9; 6 3 7 9; 3 1 7 9; 1 2 8 9; 2 5 8 9;
       ...5 4 8 9; 4 1 8 9; 4 5 6 9; 1 2 3 9];
       ...'pyra5',[2 3 6 5 9]}}),...
       'hexa8',struct('shift',[5 6],'faces',hexa8('faces'),...
       'face',{{[9 2 3 7 6; 10 3 4 8 7],1/4}},'volume',{{[11 1:8],1/8}},...
       'Elt',{{'tetra4',[2 9 3 11;  3 9 7 11;  7 9 6 11;  6 9 2 11;  4 10 8 11;
       8 10 7 11;  7 10 3 11;  3 10 4 11];
       'pyra5',[1 2 3 4 11;  1 4 8 5 11;  1 5 6 2 11;  5 8 7 6 11]}}));
      %      {{'tetra4',[2 3 9 11; 3 7 9 11; 7 6 9 11; 6 2 9 11;
      %       4 8 10 11; 8 7 10 11; 7 3 10 11; 3 4 10 11];
      %       'pyra5',[1 4 3 2 11; 1 5 8 4 11; 1 2 6 5 11; 5 6 7 8 11]}}));
      R1.set=r6; RO{end+1}=R1;
     end
     if ~isempty(r7)  % two faces no common
      R1=struct('hexa8',struct('shift',[1 4],'faces',hexa8('faces'),...
       'face',{{[9 1 4 3 2; 10 5 6 7 8],1/4}},'volume',{{[11 1:8],1/8}},....
       'Elt',{{'tetra4',[1 9 4 11;  4 9 3 11;  3 9 2 11;  2 9 1 11;  5 10 6 11;
       6 10 7 11;  7 10 8 11;  8 10 5 11];
       'pyra5',[1 4 8 5 11;  1 5 6 2 11;  2 6 7 3 11;  4 3 7 8 11]}}));
      %      {{'tetra4',[1 4 9 11; 4 3 9 11; 3 2 9 11; 2 1 9 11;
      %       5 6 10 11; 6 7 10 11; 7 8 10 11; 8 5 10 11];
      %       'pyra5',[1 5 8 4 11; 1 2 6 5 11; 2 3 7 6 11; 4 8 7 3 11]}}));
      R1.set=r7; RO{end+1}=R1;
     end
     
    case 3 % three faces, support penta, hexa(2cases)
     r7=r6; % second case
     el0=feutil('selelt eltind',mo3,full(EEid(r6(:,1)+1))); el0id=feutil('eltid;',el0);
     [EG,nG]=getegroup(el0);
     for jG=1:nG
      [ElemF,u1,ElemP]= getegroup(el0(EG(jG),:),jG);
      cEGI=EG(jG)+1:EG(jG+1)-1;
      switch ElemP
       case 'penta6'; r7(ismember(r7(:,1),el0id(cEGI)),1)=0; % always first case
       case 'hexa8'
        i2=feval(ElemP,'faces');
        [i3,i4]=ismember(r6(:,1),el0id(cEGI)); i4=i4(i3); i3=find(i3);
        for jElt=1:length(i3)
         n1=i2(r6(i3(jElt),2:end),:);
         if any(sparse(n1(:),1,ones(numel(n1),1))==3); r6(i3(jElt),1)=0;% second case one common node
         else; r7(i3(jElt),1)=0; % first case no common node
         end
        end
      end
     end % loop nG
     r6(r6(:,1)==0,:)=[]; r7(r7(:,1)==0,:)=[];
     if ~isempty(r6) % Three faces /adjacent (no common node)
      R1=struct('penta6',struct('shift',[2 3 5],'faces',...penta6('faces'),...
       [1 3 2 2;1 4 6 3;2 5 4 1;4 5 6  6;3 6 5 2],...
       'face',{{[7 1 4 6 3;8 1 2 5 4;9 2 3 6 5],[1/4]}},'volume',{{[10 1:6],[1/6]}},...
       'Elt',{{'tetra4',[1 7 4 10;  4 7 6 10;  6 7 3 10;  3 7 1 10;  1 8 2 10;
       2 8 5 10;  5 8 4 10;  4 8 1 10;  2 9 3 10;  3 9 6 10;  6 9 5 10;  5 9 2 10;
       1 3 2 10;  4 6 5 10]}}), ...
       ...{{'tetra4',[1 4 7 10; 4 6 7 10; 6 3 7 10; 3 1 7 10; 1 2 8 10; 2 5 8 10;
       ...5 4 8 10; 4 1 8 10; 2 3 9 10; 3 6 9 10; 6 5 9 10; 5 2 9 10; 1 2 3 10; 4 5 6 10]}}),...
       'hexa8',struct('shift',[1 2 4],'faces',hexa8('faces'),...
       'face',{{[9 1 4 3 2; 10 1 5 8 4; 11 5 6 7 8],[1/4]}},'volume',{{[12 1:8],[1/8]}},...
       'Elt',{{'tetra4',[1 9 4 12;  4 9 3 12;  3 9 2 12;  2 9 1 12;  1 10 5 12;
       5 10 8 12;  8 10 4 12;  4 10 1 12;  5 11 6 12;  6 11 7 12;  7 11 8 12;
       8 11 5 12];
       'pyra5',[1 5 6 2 12;  2 6 7 3 12;  4 3 7 8 12]}}));
      ...{{'tetra4',[1 4 9 12; 4 3 9 12; 3 2 9 12; 2 1 9 12;
       ... 1 5 10 12; 5 8 10 12; 8 4 10 12; 4 1 10 12;
       ... 5 6 11 12; 6 7 11 12; 7 8 11 12; 8 5 11 12];
       ... 'pyra5',[1 2 6 5 12;2 3 7 6 12;4 8 7 3 12]}}));
       R1.set=r6; RO{end+1}=R1;
     end
     if ~isempty(r7) % Three faces one common node
      R1=struct('hexa8',struct('shift',[4 5 6],'faces',hexa8('faces'),...
       'face',{{[9 5 6 7 8; 10 2 3 7 6;11 4 8 7 3],[1/4]}},'volume',{{[12 1:8 ],[1/8]}},...
       'Elt',{{'tetra4',[5 9 6 12;  6 9 7 12;  7 9 8 12;  8 9 5 12;  2 10 3 12;
       3 10 7 12;  7 10 6 12;  6 10 2 12;  4 11 8 12;  8 11 7 12;  7 11 3 12;
       3 11 4 12];
       'pyra5',[1 2 3 4 12;  1 4 8 5 12;  1 5 6 2 12]}}));
      %      ,{{'tetra4',[ 5 6 9 12; 6 7 9 12; 7 8 9 12; 8 5 9 12;
      %       2 3 10 12; 3 7 10 12; 7 6 10 12; 6 2 10 12;
      %       4 8 11 12; 8 7 11 12; 7 3 11 12;3 4 11 12];
      %       'pyra5',[1 4 3 2 12; 1 5 8 4 12; 1 2 6 5 12]}}));
      R1.set=r7; RO{end+1}=R1;
     end
     
    case 4 % support hexa(2cases)
     r7=r6; % second case
     el0=feutil('selelt eltind',mo3,full(EEid(r6(:,1)+1))); el0id=feutil('eltid;',el0);
     [EG,nG]=getegroup(el0);
     for jG=1:nG
      [ElemF,u1,ElemP]= getegroup(el0(EG(jG),:),jG);
      cEGI=EG(jG)+1:EG(jG+1)-1;
      switch ElemP
       case 'hexa8'
        i2=feval(ElemP,'faces');
        [i3,i4]=ismember(r6(:,1),el0id(cEGI)); i4=i4(i3); i3=find(i3);
        for jElt=1:length(i3)
         n1=i2(r6(i3(jElt),2:end),:);
         if any(sparse(n1(:),1,ones(numel(n1),1))==3); r6(i3(jElt),1)=0;% second case one shared edge
         else; r7(i3(jElt),1)=0; % first case no shared edge
         end
        end
      end
     end % loop nG
     r6(r6(:,1)==0,:)=[]; r7(r7(:,1)==0,:)=[];
     if ~isempty(r6) % four faces no common edge/node
      R1=struct('hexa8',struct('shift',[2 3 5 6],'faces',hexa8('faces'),...
       'face',{{[9 1 5 8 4; 10 1 2 6 5; 11 2 3 7 6; 12 4 8 7 3],[1/4]}},...
       'volume',{{[13 1:8],[1/8]}},...
       'Elt',{{'tetra4',[1 9 5 13;  5 9 8 13;  8 9 4 13;  4 9 1 13;
       1 10 2 13;  2 10 6 13;  6 10 5 13;  5 10 1 13;  2 11 3 13;
       3 11 7 13;  7 11 6 13;  6 11 2 13;  4 12 8 13;  8 12 7 13;
       7 12 3 13;  3 12 4 13];
       'pyra5',[1 2 3 4 13;  5 8 7 6 13]}}));
      %      {{'tetra4',[1 5 9 13; 5 8 9 13; 8 4 9 13; 4 1 9 13;
      %       1 2 10 13; 2 6 10 13; 6 5 10 13; 5 1 10 13;
      %       2 3 11 13; 3 7 11 13; 7 6 11 13; 6 2 11 13;
      %       4 8 12 13; 8 7 12 13 ; 7 3 12 13; 3 4 12 13];
      %       'pyra5',[1 4 3 2 13; 5 6 7 8 13]}}));
      R1.set=r6; RO{end+1}=R1;
     end
     if ~isempty(r7) % four faces one commone edge/node
      R1=struct('hexa8',struct('shift',[1 2 3 4],'faces',hexa8('faces'),...
       'face',{{[9 1 4 3 2; 10 1 5 8 4; 11 1 2 6 5; 12 5 6 7 8],[1/4]}},...
       'volume',{{[13 1:8],[1/8]}},...
       'Elt',{{'tetra4',[1 9 4 13;  4 9 3 13;  3 9 2 13;  2 9 1 13;
       1 10 5 13;  5 10 8 13;  8 10 4 13;  4 10 1 13;  1 11 2 13;
       2 11 6 13;  6 11 5 13;  5 11 1 13;  5 12 6 13;  6 12 7 13;
       7 12 8 13;  8 12 5 13];
       'pyra5',[2 6 7 3 13;  4 3 7 8 13]}}));
      %      {{'tetra4',[1 4 9 13; 4 3 9 13; 3 2 9 13; 2 1 9 13;
      %       1 5 10 13; 5 8 10 13; 8 4 10 13; 4 1 10 13;
      %       1 2 11 13; 2 6 11 13; 6 5 11 13; 5 1 11 13;
      %       5 6 12 13; 6 7 12 13; 7 8 12 13; 8 5 12 13];
      %       'pyra5',[2 3 7 6 13; 4 8 7 3 13]}}));
      R1.set=r7; RO{end+1}=R1;
     end
     
    case 5 % support hexa(id is unchanged face)
     R1=struct('hexa8',struct('shift',[6],'faces',hexa8('faces'),...
      'face',{{[9 1 4 3 2; 10 1 5 8 4; 11 1 2 6 5; 12 5 6 7 8; 13 2 3 7 6;],[1/4]}},...
      'volume',{{[14 1:8],[1/8]}},...
      'Elt',{{'tetra4',[1 9 4 14;  4 9 3 14;  3 9 2 14;  2 9 1 14;
      1 10 5 14;  5 10 8 14;  8 10 4 14;  4 10 1 14;  1 11 2 14;
      2 11 6 14;  6 11 5 14;  5 11 1 14;  5 12 6 14;  6 12 7 14;
      7 12 8 14;  8 12 5 14;  2 13 3 14;  3 13 7 14;  7 13 6 14;
      6 13 2 14];
      'pyra5',[4 3 7 8 14]}}));
     %     {{'tetra4',[1 4 9 14; 4 3 9 14; 3 2 9 14; 2 1 9 14;
     %      1 5 10 14; 5 8 10 14;8 4 10 14;4 1 10 14;
     %      1 2 11 14; 2 6 11 14; 6 5 11 14;5 1 11 14;
     %      5 6 12 14; 6 7 12 14; 7 8 12 14; 8 5 12 14;
     %      2 3 13 14; 3 7 13 14; 7 6 13 14;6 2 13 14];
     %      'pyra5',[4 8 7 3 14]}}));
     r7=true(size(r6,2),size(r6,1));  % reverse faceid
     r7(reshape(bsxfun(@plus,(0:6:(6*size(r6,1)-1))',r6(:,2:end)),[],1))=false;
     r6(:,2)=mod(find(r7),6); r6(:,3:end)=[];
     R1.set=r6;  RO{end+1}=R1;
     
    case 6 % support hexa
     R1=struct('hexa8',...'shift',[],'faces',hexa8('faces'),...
      struct('face',{{[9 1 4 3 2; 10 1 5 8 4; 11 1 2 6 5; 12 5 6 7 8;
      13 2 3 7 6; 14 4 8 7 3],[1/4]}},...
      'volume',{{[15 1:8],[1/8]}},...
      'Elt',{{'tetra4',[1 9 4 15;  4 9 3 15;  3 9 2 15;  2 9 1 15;
      1 10 5 15;  5 10 8 15;  8 10 4 15;  4 10 1 15;  1 11 2 15;
      2 11 6 15;  6 11 5 15;  5 11 1 15;  5 12 6 15;  6 12 7 15;
      7 12 8 15;  8 12 5 15;  2 13 3 15;  3 13 7 15;  7 13 6 15;
      6 13 2 15;  4 14 8 15;  8 14 7 15;  7 14 3 15;  3 14 4 15]}}));
     R1.set=r6(:,1); RO{end+1}=R1; % no shift needed
     
   end % matched faces supported cases
  end % loop on number of matched faces
  
  out=mo3; out1=[]; % old elts
  for j1=1:length(RO)
   [mo3,o1,o2]=feutil('RefineCell-replace-keepSets;',mo3,RO{j1}); %KnownNew
   out1=feutil('addelt',out1,o1);
  end
  
  % also transform quad4 to tria in a similar fashion
  R1=struct('quad4',struct('volume',{{[5 1:4],1/4}},...
   'Elt',{{'tria3',[1 2 5;2 3 5; 3 4 5; 4 1 5]}}));
  R1.set=feutil('eltid;',elt); R1.set(R1.set==0)=[];
  [mo3,o1]=feutil('RefineCell-replace-keepSets;',mo3,R1); %KnownNew
  if ~isempty(o1); out1=feutil('addelt',out1,o1); end
  
  out.Node=mo3.Node; out.Elt=mo3.Elt;
  out=stack_set(out,stack_get(mo3,'set'));
 elseif comstr(Cam,'linesurf')
 %% #cutLineSurf : form a polygonal line on a surface -----2
  
  model=varargin{carg};carg=carg+1;
  RO=varargin{carg};carg=carg+1;if iscell(RO);RO=RO{1};end
  RO.NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
  if size(RO.line,1)>1&&size(RO.line,2)==7; RO.line=RO.line(:,5:7);
  elseif size(RO.line,2)~=3; RO.line=model.Node(RO.NNode(RO.line),5:7);
  end
  %if isfield(RO,'StickCorner')&&strcmpi(RO.StickCorner,'nearest')
  %  i2=feutilb('addnode -nearest epsl 2',model.Node(:,5:7),RO.line);
  %  model.Node(i2,5:7)=RO.line;
  % end
  if isfield(RO,'CutSel')
   [RO.NonCutElt,model.Elt]=feutil(['RemoveElt' RO.CutSel],model); 
   model=feutil('joinall',model);
  end
  sel=feutil('getpatchnew',model);
  RO=sdth.sfield('addmissing',RO, ...
      feval(feutilb('@OEdgeElt'),struct,model,'adjacentsurf'));  
  RO.vert=sel.vert0; RO.Node=sel.Node; % prepare for dynamic
  RO.elt=sel.fs.'; RO.FindRST=feutilb('@FindRST'); %RO.Edge nodeindex
  RO.q4=integrules('quad4','inline');RO.t3=integrules('tria3','inline');
  if ~isfield(RO,'diag');RO.diag=0;end
  % start 1 nearest node 
  % 2. NotANode -> findRS (no possible for now)
  %    IsANode -> ClosestEdge (smallest angle), then EdgeCut
  % 3. tolMoveNode
  RO.nmap=containers.Map; RO.incuts=[];
  for j1=[2:size(RO.line,1)] % deal with corners xxx closed line
   [un1,RO.inearest]=min(sum((RO.vert-repmat(RO.line(j1,:),size(RO.vert,1),1)).^2,2));
   [un1,i2]=find(RO.Edge==RO.inearest);% i2 edgenumbers
   i3=ismember(RO.EltEdge,i2);i3=find(any(i3,1)); % i3 associated RO.elt cols
   RO.curMatchE=i3; %iz=RO.Node(RO.elt(:,1:4));fecom('shownodemark',iz(:))
   n2=RO.line(j1,:);RO=findRS(RO,n2,struct('do','moveinsert'));
  end
  RO.start=RO.line(1,:); RO.end=RO.line(2,:); RO.dir=RO.end-RO.start;RO.dir=RO.dir/norm(RO.dir);
  RO.cuts=RO.start; RO.isFirstNode=1;
  RO.cutelt=[]; RO.ls=nan(size(RO.vert,1),1);RO.incuts(1)=RO.inearest(1);
  RO.iecuts=[]; cEGI=find(isfinite(model.Elt(:,1)));
 while 1
  % i2 nodes connected by edges, i3 elements with these edges
  [un1,i2]=find(RO.Edge==RO.inearest);i3=ismember(RO.EltEdge,i2);i3=find(any(i3,1));
  if isempty(i2); dbstack; keyboard;end
  RO.curMatchE=i3;
  i2=RO.Edge(:,i2);i2(i2==RO.inearest)=0;i2=sum(i2,1);
  r1=RO.vert(i2,:)-RO.start;r2=sqrt(sum(r1.^2,2));r1=diag(sparse(1./r2))*r1;
  % Point within next elt (step less than shortest edge)
  n2=RO.start+RO.dir*min(r2)*RO.tolE;RO=findRS(RO,n2,struct('do','match'));
  % fecom('shownodemark',{RO.vert(RO.NNode(reshape(RO.elt(:,RO.curMatchE),[],1)),:);n2},'marker','o');
  % approximate levelset for cut element
  ns=RO.vert(RO.elt(:,RO.curMatchE),:);
  p1=sp_util('basis',ns'); p1=cross(RO.dir,p1(:,3));p1=p1/norm(p1);% Local basis   
  r3=(ns-repmat(RO.start,size(ns,1),1))*p1';r3(abs(r3)<RO.epsl)=0;
  st=sprintf('f%i.%i.%i.%i',sort(RO.elt(:,RO.curMatchE)));
  if RO.nmap.isKey(st)&&~RO.isFirstNode; % element contains intra element corner
       RO.cuts(end+1,:)=RO.end; RO.start=RO.end;RO.incuts(size(RO.cuts,1))=RO.nmap(st);
       RO.line(1,:)=[];if size(RO.line,1)==1; break;end;RO.end=RO.line(2,:); 
       RO.dir=RO.end-RO.start;RO.dir=RO.dir/norm(RO.dir); RO.isFirstNode=1;
       continue;
  end
  RO.ls(RO.elt(:,RO.curMatchE))=r3; 
  % EdgeExitPoint
  i3=RO.EltEdge(:,RO.curMatchE);   %RO.curMatchEdgefecom('shownodemark',reshape(RO.Edge(:,i3),[],1))
  for j3=1:length(i3)
   i4=RO.Edge(:,i3(j3));st=sprintf('e%i.%i',sort(i4));
   if RO.nmap.isKey(st);continue;end
   n2=RO.vert(i4,:);
   r=findR(n2,[RO.start;RO.end]); %fecom('shownodemark',[n2;RO.start;RO.end])
   if norm(n2(1,:)-RO.start)<RO.epsl||norm(n2(2,:)-RO.start)<RO.epsl
     % too close to start node (not a proper cut)
      continue;
   elseif r(4)>=-RO.epsl; % this is a cut edge
     n3=[1-r(1) r(1)]*n2; if (n3-RO.start)*RO.dir(:)<0; continue;end%wrong dir
     if r(1)<RO.tolE; %sufficiently close for movenode1
      i2=RO.Edge(1,i3(j3)); RO.vert(i2,:)=n3; RO.inearest=i2;
     elseif abs(1-r(1))<RO.tolE; %sufficiently close for movenode
      i2=RO.Edge(2,i3(j3)); RO.vert(i2,:)=n3; RO.inearest=i2;
     else; % Now need to add node with edge label 
      RO.vert(end+1,:)=n3; i2=size(RO.vert,1);RO.nmap(st)=i2;
      if r(1)<.5; RO.inearest=RO.Edge(1,i3(j3));
      else;RO.inearest=RO.Edge(2,i3(j3));
      end      
     end
     RO.start=RO.vert(i2,:);  RO.isFirstNode=0;
     RO.cuts(end+1,:)=RO.start; RO.incuts(size(RO.cuts,1))=i2;
     RO.iecuts(end+1)=cEGI(RO.curMatchE);
     if RO.diag;fecom('shownodemark',{RO.cuts,[RO.start;RO.end]},'marker','o'); end
     if norm(RO.start-RO.end)<RO.epsl; % Done segment
       RO.line(1,:)=[];if size(RO.line,1)==1; break;end;RO.end=RO.line(2,:); 
       RO.dir=RO.end-RO.start;RO.dir=RO.dir/norm(RO.dir);
     end
     break;
   end
  end
  if size(RO.line,1)==1; break;end
 end
 out=model; i1=length(RO.nmap);
 if i1~=0 % Some nodes were added
  i2=max(out.Node(:,1))+(1:i1)'; 
  out.Node=[out.Node;[i2*[1 0 0 0],RO.vert(end+[-i1+1:0],:)]];
 end
 out.Node(RO.NNode(RO.Node),5:7)=RO.vert(1:length(RO.Node),:);%Moved nodes
 def=struct('def',RO.ls,'DOF',out.Node(:,1)+.98);def.def(end+1:length(def.DOF),1)=0;
 def.def(RO.incuts)=0;
 if isfield(RO,'NonCutElt')&&~isempty(RO.NonCutElt);
     out=feutil('addelt',out,RO.NonCutElt);
 end
 for j1=1:size(model.Node,1);RO.nmap(sprintf('n%i',j1))=j1; end
 if isfield(RO,'ProId')
   % fill sign field by connectivity for selection
   while 1 
     r1=def.def(RO.Edge); 
     i2=~isfinite(r1(1,:))&isfinite(r1(2,:))&r1(2,:)~=0;
     if any(i2); def.def(RO.Edge(1,i2))=def.def(RO.Edge(2,i2));continue;end
     i2=~isfinite(r1(2,:))&isfinite(r1(1,:))&r1(1,:)~=0;
     if any(i2); def.def(RO.Edge(2,i2))=def.def(RO.Edge(1,i2));continue;end
     break;
   end
 end
 RO.def=def.def; RO.inEltInd=RO.iecuts;
 %cf=feplot(out); cf.sel=['eltind', sprintf('%i ',RO.inEltInd)];fecom('shownodemark',RO.cuts,'marker','o')
 %lsutil('viewls',out,def); fecom('shownodemark',RO.cuts,'marker','o')
 out=genericNameMap('docut',out,RO);
 return
  
 else%if comstr(Cam,'cut');[CAM,Cam]=comstr(CAM,4);
  %% #cutRemesh : remesh with an interface in levelset -----------------------2
  model=varargin{carg};carg=carg+1;
  li=varargin{carg};carg=carg+1;
  if nargin>=carg;RO=varargin{carg};else;RO=struct;end
  if isfield(RO,'sel')
    model.Elt=feutil(['selelt' RO.sel],model);
  end
  if ~isfield(RO,'tolE'); RO.tolE=.1; end
  if ~isempty(strfind(Cam,'-keeporigmpid')); RO.keepOrigMPID=1; end
  if ~isfield(RO,'keepSets'); RO.keepSets=[0 0]; end
    
  if isfield(model,'vert0') % Actually a selection start by going to model
   RO.InitMdl=model;sel=model;
   mo1=struct('Node',[sel.Node*[1 0 0 0] sel.vert0],'Elt',[]);
   if size(model.fs,2)<4;model.fs(:,4)=model.fs(:,3);end
   i1=model.fs(:,3)==model.fs(:,4);sel.fs=sel.Node(sel.fs);
   mo1=feutil('addelt',mo1,'tria3',[sel.fs(i1,1:3) sel.MPID(sel.ifs(i1),:)]);
   i1=~i1;
   mo1=feutil('addelt',mo1,'quad4',[sel.fs(i1,1:4) sel.MPID(sel.ifs(i1),:)]);
   [eltid,mo1.Elt]=feutil('eltidfix;',mo1);
   model=mo1;
  end
  if isfield(RO,'onSurf')
   st=RO.onSurf;RO.onSurf=feutil(RO.onSurf,model);RO.onSurfSel=st;
  end
  if ~isfield(RO,'Fixed');RO.Fixed=[];
  elseif ischar(RO.Fixed);RO.Fixed=feutil(['findnode',RO.Fixed],model);
  end
  
  if isfield(RO,'doSplit')&&RO.doSplit;% Start by splitting intersected elts
   error('.doSplit obsolete : now on the fly splitting');
   if any(RO.keepSets)&&~isfield(RO,'gset')
    RO.gset=feutil('addseteltid-append-get -NoNodes',model,'_gset');
   end
   [model,un1,RO]=lsutil('split',model,li,RO);
   RO=feutil('rmfield',RO,'conn','edges','values','NNode','ed2elt','elt2ed');
   if isfield(RO,'onSurf');RO.onSurf=feutil(st,model);end
  end
  if isfield(RO,'sel') % Possibly perform on sub model
   [RO.Elt,model.Elt]=feutil(['removeelt' RO.sel],model);
  end
  
  cpt=0;def=lsutil('gen-max',model,li);
  while 1
   
   cpt=cpt+1;   
   %search elements separated by cut edges (once only !)
   % define EdgeIndToEltInd and EltIndToEdgeInd (in RO.edges, and model.Elt) -xxxeb
   RO=feutil('rmfield',RO,'elt2ed','ed2elt','conn','edges');RO.ed2elt=[];
   %% if tolE decide on moving nodes with acceptable tolerances
   [RO,model,def]=lsutil('EdgeCut',model,def,RO);
   
   %def.sel=li{1}.sel;lsutil('ViewLs',model,def);
   [model,RO,def]=tolMoveNode(model,RO,def); %xxx already in edgecut
   % fecom('shownodemark',def.DOF(def.def==0))
   [model,isplit,def]=check_split(model,def);
   if isplit;break;end
  end
  
  %% #create_new_nodes of those outside tol -3
  i1=RO.r>RO.newTol & RO.r<=1-RO.newTol;
  if ~isfield(def,'distFcn');error('Obsolete');%xnew=cleanNew(model,RO,li,i1);
  else; [xnew,def]=cleanNew(model,RO,def,i1);
  end
  
  
  mo3=struct('Node',[],'Elt',[]);
  if isfield(RO,'disJnodes')
   % need to addnodes by node subentries
   id2=zeros(size(xnew,1),1); mo3.Node=model.Node; in0=zeros(size(RO.edges,1),1);
   for j1=1:length(RO.disJnodes)
    in1=find(any(ismember(RO.edges,RO.disJnodes{j1}),2)); in0(in1)=1;
    in2=find(ismember(mo3.Node(:,1),RO.disJnodes{j1}));
    [n3,id3]=feutil('addnode',mo3.Node(in2,:),xnew(in1,:));
    % new nodes to be added with knownnew and renewed index
    in3=id3>length(in2);
    n4=n3(id3(in3),:); % new nodes to add
    % new nodes to be added with knownnew and renewed index
    [mo3.Node,id4]=feutil('addnodeknownnew',mo3.Node,xnew(in1(in3),:));
    % coalesced nodes to be matched to initial index
    id3(in3)=id4; id3(~in3)=in2(id3(~in3));
    id2(in1)=id3;
   end
   if any(id2==0)
    [mo3.Node,id3]=feutil('addnode',mo3.Node,xnew(id2==0,:));
    id2(id2==0)=id3;
   end
  else
   [mo3.Node,id2]=feutil('addnode',model.Node,xnew);
  end
  def.DOF=[def.DOF;mo3.Node(id2,1)+.98];def.def=[def.def;zeros(size(id2))];%known new
  % Row due to later formatting
  RO.idnew=zeros(1,size(RO.edges,1));RO.idnew(i1)=mo3.Node(id2); idnew=RO.idnew;
  RO.NNode=sparse(mo3.Node(:,1),1,1:size(mo3.Node,1));
  
  % asmany rows as edges, new id number for edges that have
  %feplot(mo3);hold on;plot3(xnew(:,1),xnew(:,2),xnew(:,3),'r*')
  %'xxx'
  %def.def(abs(def.def)<norm(RO.values,'inf')*1e-5)=0;
  
  [EGroup,nGroup]=getegroup(model.Elt);
  listE=spalloc(size(model.Elt,1),1,size(RO.cEGI,1));
  tria=[];quad=[];tetra=[];hexa=[];pyra=[];penta=[];Beam1=[];
  idedges=reshape(model.Node(RO.edges,1),size(RO.edges));
  for jGroup=1:nGroup
   [ElemF,opt,ElemP] = feutil('GetElemF',model.Elt(EGroup(jGroup),:),jGroup);
   if opt(1)<0;continue;end
   eltind=(EGroup(jGroup)+1:EGroup(jGroup+1)-1)';prop=feval(ElemP,'prop');
   
   in0=ismember(RO.cEGI,eltind);cEGI=RO.cEGI(in0);
   if isempty(cEGI);continue;end
   elt=model.Elt(cEGI,:);%elt(1,end+1:prop(3))=0;
   if size(elt,2)<max(prop); elt(end,max(prop))=0; end
   elt2ed=RO.elt2ed(in0,:);
   
   i1=reshape(full(RO.NNode(elt(:,feval(ElemP,'node')))),size(elt,1),[]);
   sevLS=reshape((sign(def.def(i1))),size(i1));% sign of level set
   
   lsC=charLs(sevLS);% No need to cut
   %i3=cellfun(@(x)length(setdiff(unique(x),'0'))<=1,lsC);
   i3=cellfun(@(x)~(any(x=='+')&any(x=='-')),lsC);
   RO.cEGI=setdiff(RO.cEGI,cEGI(i3));
   elt(i3,:)=[];elt2ed(i3)=[];lsC(i3)=[];sevLS(i3,:)=[];cEGI(i3)=[];
   if isempty(cEGI);continue;end
   
   % RO.r(vertcat(elt2ed{:}))
   % cf=feplot(model);cf.sel=['eltind ' sprintf('%i ',cEGI)];
   if sum(strcmp(ElemP,{'quad4','quadb'}))
    %% #quad4 -3
    genericCut(reQuad,elt,ElemP,sevLS,idedges,RO); listE(cEGI)=1;
    
   elseif strcmp(ElemP,'tria3')
    %% #tria3 -3   
    genericCut(reTria,elt,ElemP,sevLS,idedges,RO); listE(cEGI)=1;
    
   elseif sum(strcmp(ElemP,{'tetra4'}))
    %% #tetra4 -3
    genericCut(reTetra,elt,ElemP,sevLS,idedges,RO); listE(cEGI)=1;
    
   elseif sum(strcmp(ElemF,{'pyra5'})) %% #pyra5 -3
    i1=genericCut(rePyra,elt,ElemP,sevLS,idedges,RO); listE(cEGI(i1))=1;
    
   elseif any(strcmp(ElemP,'beam1'));
    %% #beam1 -3
    
    %one cut
    allcases={[1 -1;-1 1],[1 2],'[curelt(:,[1 3]) curprop;curelt(:,[3 2]) curprop]'};
    
    RE=struct('cases',vertcat(allcases{:,1}), ...
     'cindex',reshape([1;1]*(1:size(allcases,1)),[],1));
    RE.newelt=[{'Beam1'};allcases(:,3)];RE.CutEdges=allcases(:,2);
    genericCut(RE,elt,ElemP,sevLS,idedges,RO);
    
   elseif any(strcmp(ElemP,'mass1'));
   elseif any(strcmp(ElemP,'hexa8'));
    %% #hexa8 -3
    i1=genericCut(reHexa,elt,ElemP,sevLS,idedges,RO);
    listE(cEGI(i1))=1;
    
   else;
    %% #penta6 -3
    i1=genericCut(rePenta,elt,ElemP,sevLS,idedges,RO);
    listE(cEGI(i1))=1;
    %     lsC=charLs(sevLS);
    %     i3=cellfun(@(x)length(setdiff(unique(x),'0'))<=1,lsC);% No need to cut
    %     if any(~i3)
    %      i1=genericCut(rePenta,elt,ElemP,sevLS,idedges,RO);listE(cEGI(i1))=1;
    %      sdtw([ElemP,' elements are not expected use lsutil(''split'') first']);
    %      disp(elt(~i3,:));%cf=feplot;cf.sel=['eltid ', sprintf('%i ',elt(~i3,11))]
    %      %warning('element not implemented');keyboard
    %      %mo4=model;mo4.Elt=[inf,abs(ElemF),zeros(1,size(elt,2)-length(ElemF)-1);elt(~i3,:)];
    %      %RA=rmfield(RO,{'conn','ed2elt','elt2ed','NNode','idnew','values','edges','EEIdx','cEGI'});
    %      %[mo5,toto1,toto2]=lsutil('split',mo4,li,RA); mo6=lsutil('cut',mo4,li,RA);
    %     end
   end
  end
  %add new elements to model
  [mo3,RO]=safeAddElt(mo3,RO,'tria3',tria);
  [mo3,RO]=safeAddElt(mo3,RO,'quad4',quad);
  [mo3,RO]=safeAddElt(mo3,RO,'tetra4',tetra);
  [mo3,RO]=safeAddElt(mo3,RO,'pyra5',pyra);
  [mo3,RO]=safeAddElt(mo3,RO,'penta6',penta);
  [mo3,RO]=safeAddElt(mo3,RO,'hexa8',hexa);
  [mo3,RO]=safeAddElt(mo3,RO,'beam1',Beam1);
  % feplot(mo3);fecom('shownodemark',xnew)
  if isfield(RO,'PostCutFcn')
   if ischar(RO.PostCutFcn);eval(RO.PostCutFcn);
   else; dbstack; keyboard;
   end
  end
  if ~isempty(mo3.Elt)
   %remove cut elements from model
   if any(RO.keepSets)&&~isfield(RO,'gset')
    RO.gset=feutil('addseteltid-append-get -NoNodes',model,'_gset');
   end
   [model.Elt,out2]=feutil('removeelt eltind',model,find(listE));
   model.Node=mo3.Node;
   model=feutil('addelt',model,mo3.Elt);%merge models
   if any(RO.keepSets) % RO.EEIdx has list of [oldEltId newEltId]
    model=feutil('EltSetReplace',model,RO.gset,RO);
   end
   out=model;
  else
   out=model; out2=[];
   if ~isempty(RO.cEGI) % Some elements need cutting and empty model
    warning('empty model');dbstack;%keyboard;
   end
  end
  if isfield(RO,'onSurf') % CheckEdges
   RO.mTol=-0e-3;RO.onSurf=feutil(RO.onSurfSel,model);
   if 1==2 % gv unsure edges check is needed
   i1=RO.edges(all(RO.values>=RO.mTol,2),:);
   i1=i1(ismember(i1(:,1),RO.onSurf)&ismember(i1(:,2),RO.onSurf),:);
   i2=fix(def.DOF(def.def>=RO.mTol));i2=intersect(i2,RO.onSurf);
   i2=[i1(:);i2];
   else
    i2=fix(def.DOF(def.def>=RO.mTol));i2=intersect(i2,RO.onSurf);
   end
   [i3,out.Elt]=feutil('eltidfix;',out);
   if isfield(RO,'keepOrigMPID')&&RO.keepOrigMPID
   else; out=lsutil('mpid',out,li,RO);
   end
   out.Elt=feutil('orient;',out);
   if ~isempty(i2)
    out=feutil('AddSetFaceId',out,'onSurf',['selface & innode' sprintf('%i ',i2)]);
   end
  elseif 1==2
   %lsutil('ViewLs',model,li);fecom('shownodemark',[i1(:);i2],'marker','o','color','r')
   %fecom('coloredgealpha1'); 'xxx'
   cf=feplot(out);
   cf.sel=['selface & innode' sprintf('%i ',i2)];
   cf.sel='setname onSurf';
  elseif isfield(RO,'InitMdl');
   %   lsutil('viewls',model,li)
   elt=out.Elt;mpid=feutil('mpid',elt);
   r1=zeros(size(model.Elt,1),3);r1(isfinite(out.Elt(:,1)),:)=1e6;
   out.Elt=feutil('mpid',out.Elt,r1);
   out=lsutil('mpid',out,li,RO);r1=feutil('mpid',out);
   %mpid(r1(:,1)==1e6,:)=1e6;
   %out.Elt=feutil('mpid',out.Elt,mpid);
   %out.Elt=feutil('removeelt matid 1e6',out);
   out.Elt=feutil('removeelt eltind',out,find(r1(:,1)==1e6));
  elseif isfield(RO,'keepOrigMPID')&&RO.keepOrigMPID
  else;
   try;out=lsutil('mpid',out,{def},RO);
   catch;warning('MPID failed');
   end
  end
  if isfield(RO,'sel') % Possibly perform on sub model
   out.Elt=feutil('addelt',RO.Elt,out.Elt);
  end
  
  if nargout>1;out1=lsutil('gen-max',out,li);end
  
 end % cut
 
elseif comstr(Cam,'mpid');[CAM,Cam]=comstr(CAM,4);
 %% #mpid : assign mpid based on level set ----------------------------------
 model=varargin{carg};carg=carg+1;li=varargin{carg};carg=carg+1;
 if nargin>=carg;RO=varargin{carg};else RO=struct;end
 
 [def,mid]=lsutil('gen-max',model,li);%warning('obsolete');
 RO.mpid=zeros(size(li,1),2);
 [EGroup,nGroup]=getegroup(model.Elt);mpid=feutil('mpid',model);
 mpid(:,1:2)=-mpid(:,1:2);
 for j1=1:size(li,1);
  if isfield(li{j1},'MatId');   RO.mpid(j1,1:2)=[li{j1}.MatId li{j1}.ProId];
  elseif isfield(li{j1},'mpid');RO.mpid(j1,1:2)=li{j1}.mpid;
  else; 
    try; r1=def.distFcn('ro');catch;r1.mpid=[j1 j1];end
    if isfield(r1,'MatId'); RO.mpid(j1,1:2)=[r1.MatId r1.ProId];%p_piezo 
    else; RO.mpid(j1,1:2)=r1.mpid(1:2);
    end
    
  end
 end
 if ~isfield(RO,'NNode'); RO.NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));end
 
 for jGroup=1:nGroup
  [ElemF,opt,ElemP] = feutil('GetElemF',model.Elt(EGroup(jGroup),:),jGroup);
  cEGI=(EGroup(jGroup)+1:EGroup(jGroup+1)-1)'; i1=feval(ElemP,'node');
  if strcmpi(ElemF,'celas');continue;end %
  i2=reshape(full(RO.NNode(model.Elt(cEGI,i1))),size(cEGI,1),[]);
  if 1==1
   dmoy=~all(reshape(def.def(i2),size(i2))>=0,2);% elements inside (d>0)
  elseif 1==1
   dmoy=mean(reshape(def.def(i2),size(i2,1),[]),2);
   %if sum(dmoy==0);dbstack; keyboard;end
   dmoy=(dmoy<=0 | isnan(dmoy)); % remove elements outside
  end
   cEGI(dmoy)=[];i2(dmoy,:)=[];
   [dmax,idmax]=max(abs(def.def(i2)),[],2);
   if isempty(i2); continue;
   elseif size(i2,1)==1;IND = sub2ind(size(model.Elt),cEGI,idmax(1));
   else; IND = sub2ind(size(model.Elt),cEGI,idmax);
   end
   IND=full(RO.NNode(model.Elt(IND)));
   try
    mpid(cEGI,1:2)=RO.mpid(mid(IND),:);
   catch; fprintf('Problem with MPID\n');
   end
 end
 out1=mpid; mpid=abs(mpid);
 model.Elt=feutil('mpid',model,mpid);
 try; 
   r1=def.distFcn('ro');
   if isfield(r1,'bas');
       model=r1.distFcn('ebas',r1,model);
   end
 end
 for j1=1:size(li,1);
  if isfield(li{j1},'PostMpidFcn')
   feval(li{j1}.PostMpidFcn{:});
  end
 end
 out=model;
 
%% #View 
elseif comstr(Cam,'view');[CAM,Cam]=comstr(CAM,5);
 
 if comstr(Cam,'ls')
  %% #ViewLs -------------------------------------------------------------------
  model=varargin{carg};carg=carg+1;def=varargin{carg};carg=carg+1;
  if isa(model,'sdth');cf=model;
   model=cf.mdl.GetData;
  else; cf=feplot;
   cf.SelF{1}=[];feplot(cf,'initmodel',model);
  end
  
  if ~isfield(def,'def')
   def=lsutil('gen-max',model,def);
  end
  cf.SelF{1}.cna{1}=spalloc(numel(cf.SelF{1}.vert0),size(def.def,1),0);
  feplot(cf,'initdef',def);
  fecom(cf,'colordata98 -edgealpha.1;');
  fecom(';colorbar;colorscale 2sided instant');
  if sdtdef('verm')<804;set(gcf,'renderer', 'zbuffer');end
  if sdtdef('verm')>804;fecom('colormap',parula(4))
  else; fecom('colormap',jet(4))
  end
  if isfield(def,'sel');
    cf.sel(1)={['withnode ' sprintf('%i ',fix(def.DOF(def.def>0)))], ...
        'colordata 98 -edgealpha.5 -alpha1'};
    cf.o(1)='sel 1 def 1 ch 1 ty1 scc 1'; % mesh
    cf.sel(2)={'groupall', ...
        'colordata 98 -edgealpha.5 -alpha.2'};
    cf.o(2)='sel 2 def 1 ch 1 ty1 scc 1'; % mesh
    sel=def.sel;
    sel=sdsetprop(sel,'fsProp','EdgeAlpha',.0,'FaceAlpha',.2,'facecolor','r');
    cf.SelF{3}=sel;
    cf.o(3)='sel 3 def 0 ch 0 ty1 scc 1'; % mesh
  end
  cf.ua.axProp={'@EndFcn','set(gf,''clim'',[-1 1])'};feplot
  
 elseif comstr(Cam,'scat')
  %% #ViewScat : surface mesh of yarn -------------------------------------
   model=varargin{carg};carg=carg+1;
   if carg<=nargin;RO=varargin{carg};carg=carg+1;else;RO=struct;end
   if isfield(RO,'box');
    r1=min(RO.box(2,:));
    i1=feutil('findnode x>=&x<=',model,RO.box(1)-r1,sum(RO.box(1:2,1))+r1);
    model.Elt=feutil('selelt withnode',model,i1);
   end
   cf=feplot(model); 
   fecom('colordatamat -alpha.9 -edgealpha.1');
   cf.SelF{1}=sdsetprop(cf.SelF{1},'f2Prop','edgealpha',.1);
   feplot
 else; error('View%s',CAM);
 end
 
 
 %% #CVS ----------------------------------------------------------------------
elseif comstr(Cam,'cvs')
 out='$Revision: 1.139 $  $Date: 2021/12/23 16:10:03 $';
elseif comstr(Cam,'@'); out=eval(CAM);
 %% ------------------------------------------------------------------------
else;error('%s unknown',CAM);
end

%% #SubFunc ------------------------------------------------------------------

%% #d(Function) level set (distance functions) -2
%% #dToCell : combination of level set functions -3
function out=dToCell(xyz,R1,Cam);

if ischar(xyz)
 if strcmpi(xyz,'init');out=r1;
 elseif strcmpi(xyz,'ro');out=R1;
 else; error('%s unknown',xyz);
 end
 return
else
 if isstruct(xyz);a=R1;xyz=R1;R1=a;end
end
 phi=cellfun(@(x)feval(x.distFcn,xyz),R1,'uni',0);phi=horzcat(phi{:});
 if ~isempty(strfind(Cam,'max'));out=max(phi,[],2);end
out(abs(out)<sp_util('epsl'))=0;

%% #dToRect : 2d rectangular distance -3
function out=dToRect(xyz,R1,model);%#ok<DEFNU>


if ischar(xyz)
 if strcmpi(xyz,'init')
  if ~isfield(R1,'contour')
     r2=[R1.xc R1.yc 0; R1.lx 0 0;0 R1.ly 0];
     r2(2:3,:)=[cosd(R1.alpha) sind(R1.alpha); ...
      -sind(R1.alpha) cosd(R1.alpha)]*r2(2:3,:)/2;
     R1.contour=r2(ones(4,1),:)+[1 1;1 -1;-1 -1;-1 1]*r2(2:3,:);
  end    
  out=R1; 
 elseif strcmpi(xyz,'ro');out=R1;
 else; error('%s unknown',xyz);
 end
 return
end
if isstruct(xyz);a=xyz;xyz=R1;R1=a;end
%if ~isfield(R1,'alpha');R1.alpha=0;end
R=[cosd(R1.alpha) -sind(R1.alpha);sind(R1.alpha) cosd(R1.alpha)];
xd=R*[R1.lx -R1.lx;R1.ly -R1.ly]/2+[R1.xc;R1.yc]*[1 1];
phi1=-[(xyz(:,1)-xd(1,1)) (xyz(:,2)-xd(2,1))]*R; %right top
phi2=[(xyz(:,1)-xd(1,2)) (xyz(:,2)-xd(2,2))]*R; %left bottom
out=min([phi1 phi2],[],2);
out(abs(out)<sp_util('epsl'))=0;


%% #dToCirc : 2dcircle (x-xc)^2 + (y-yc)^2 - R^2 =0  -3
function out=dToCirc(xyz,R1,model);
if isstruct(xyz);a=xyz;xyz=R1;R1=a;end
if ischar(xyz)
 if strcmpi(xyz,'init')
  out=R1; 
 elseif strcmpi(xyz,'ro');out=R1;
 else; error('%s unknown',xyz);
 end
 return
end

out= R1.rc - sqrt( (xyz(:,1)-R1.xc).^2 + (xyz(:,2)-R1.yc).^2 ) ;
out(abs(out)<sp_util('epsl'))=0;
%out=struct('DOF',model.Node(:,1)+.98,'def',out);

%% #dToSphere: distance to segment --------------------- - -3
function out=dToSphere(xyz,R1,model); %#ok<DEFNU>

   % sphere: (x-xc)^2 + (y-yc)^2 + (z-zc)^2 - R^2 =0
if ischar(xyz)
 if strcmpi(xyz,'init')
  st=setdiff({'rc','xc','yc','zc'},fieldnames(R1));
  if ~isempty(st);error('Expecting %s',comstr(st,-30));end
  out=R1; 
 elseif strcmpi(xyz,'ro');out=R1;
 else; error('%s unknown',xyz);
 end
 return
end

out=R1.rc - sqrt( (xyz(:,1)-R1.xc).^2 + (xyz(:,2)-R1.yc).^2 + (xyz(:,3)-R1.zc).^2 ) ;
out(abs(out)<sp_util('epsl'))=0;

%% #dToSeg: distance to segment ---------------------- -3
function out=dToSeg(xyz,R1,model); %#ok<DEFNU>

if ischar(xyz)
 if strcmpi(xyz,'init')
     if isfield(R1,'elt') 
      % Give first and second node
      n1=feutil('getnode',model,R1.elt(:));
      R1.orig=n1(1:size(R1.elt,1),5:7);
      R1.normal=n1(size(R1.elt,1)+1:end,5:7)-R1.orig;
      R1.z0=zeros(size(R1.orig,1),1);
      R1.z1=sqrt(sum(R1.normal.^2,2));
     elseif ~isfield(R1,'orig')
      R1.orig=[R1.xc R1.yc R1.zc];
     end
     if ~isfield(R1,'normal')
      R1.normal=[R1.nx R1.ny R1.nz];
     end
     if ~isfield(R1,'z0');R1.z0=0;end
     if ~isfield(R1,'z1');R1.z1=0;end
     if ~isfield(R1,'tip');R1.tip='round';end
   out=R1; 
   return;
 elseif strcmpi(xyz,'ro');out=R1;
 else; error('%s unknown',xyz);
 end
 return
end

r1=xyz-repmat(R1.orig,size(xyz,1),1);
%p=sp_util('basis',[R1.nx,R1.ny,R1.nz],[0 0 0]);
r1=r1*sp_util('basis',R1.normal,[0 0 0]);
out=sqrt(sum(r1(:,2:3).^2,2)); % Dist to seg
switch lower(R1.tip)
case 'flat' % Flat top 
 if isfield(R1,'rc');out=out-R1.rc;end % negative inside
 out(abs(out)<sp_util('epsl'))=0; 
 i2=r1(:,1)<R1.z0&out<0;  i3=r1(:,1)>R1.z1&out<0; % fecom('shownodemark',xyz(i1,:))
 r2=-r1(i2,1)+R1.z0; r3=r1(i3,1)-R1.z1; out(i2)=r2; out(i3)=r3; 
 %i3=abs(r2)<abs(out(i1));i1=find(i1);out(i1(i3))=r2(i3);
 %r2=r1(i1,1)-R1.z1; i3=abs(r2)<abs(out(i1));i1=find(i1);out(i1(i3))=r2(i3);
 %fecom('shownodemark',[R1.orig;R1.orig+R1.normal])
 
otherwise % Rounded
 i1=r1(:,1)<R1.z0; out(i1)=sqrt(out(i1).^2+(r1(i1,1)-R1.z0).^2);% rounded top
 i1=r1(:,1)>R1.z1; out(i1)=sqrt(out(i1).^2+(r1(i1,1)-R1.z1).^2);% rounded bot
 if isfield(R1,'rc');out=R1.rc-out;end
end


%% #dToPlane: distance to plane - - --------------------------------------- -3
function out=dToPlane(xyz,orig,normal);

if ischar(xyz)
 if strcmpi(xyz,'init')
  R1=orig;
     if ~isfield(R1,'orig')
      R1.orig=[R1.xc R1.yc R1.zc];
     end
     if ~isfield(R1,'normal')
      R1.normal=[R1.nx R1.ny R1.nz];
     end
  out=R1; 
 elseif strcmpi(xyz,'ro');out=R1;
 else; error('%s unknown',xyz);
 end
 return
end

if isstruct(orig);normal=orig.normal; orig=orig.orig;end
r1=xyz-ones(size(xyz,1),1)*orig;
p=sp_util('basis',normal,[0 0 0]);p=p(:,1);
out=r1*p;


%% #dToTri : distance to oriented triangulated surface -3
function [out]=dToTri(xyz,RO,model)

if ischar(xyz) 
 switch lower(xyz)
 case 'init'
    model=RO; % Surface mesh
    model.Elt=feutil('selelt selface',model); model=feutil('quad2tria',model);
    n1=feutil('getnode inelt {eltname tria}',model);
    model.Elt=feutil('selelt eltname tria',model);
    model.Node=feutil('getnode groupall',model);model.Stack={};
    model=feutil('renumber-noori',model,(1:size(model.Node,1))');
    cEGI=find(isfinite(model.Elt(:,1)));
    model.en=sparse(reshape(model.Elt(cEGI,1:3),[],1),[cEGI;cEGI;cEGI],1)';
    model.normal=feutil('getnormal',model);% face normals
    model.faces=model.Elt(cEGI,1:3)'; %[n1,n2,n3,e1,e2,e3,f1] % index for normals
    
    r2=feval(feutilb('@OEdgeElt'),struct,model,'adjacentsurf');% For surface adjacent element
    model=sdth.sfield('addmissing',model,r2);model.normal=model.normal(cEGI,:);
    %model.NodeNormal=r2.NodeNormal;
    %model.normal=[r2.NodeNormal r2.EdgeNormal model.normal(cEGI,:)'];%node normal
    %model.faces(end+(1:3),:)=reshape(i4+size(r1,1),3,[]);% index of edgenormal
    %model.faces(end+1,:)=size(r1,1)+size(r2,1)+(1:size(model.faces,2));% index of facenormal
    r1=model.Node(:,5:7);r2=std(r1); 
    if any(r2==0); r1(end+1,:)=max(r1)+[1 2 3];end
    model.DT=delaunayTriangulation(r1(:,1),r1(:,2),r1(:,3));
    model.sel=feutil('getpatchnew',model); 
    model.shape='Tri';
    out=model;
 case 'ro'
    out=RO; out.distFcn=@dToTri;
 case 'ebas' 
  % set element wise coordinates based on nearest node
  def=evalin('caller','def');
  n1=model.Node(ismember(model.Node(:,1),fix(def.DOF(def.def>=0))),:);
  cEGI=feutil('findelt innode',model,fix(def.DOF(def.def>=0)));
  r2=struct('Node',model.Node,'Elt', ...
      feutil('selelt eltind',model,cEGI));
  if isempty(r2.Elt);out=model;return;end
  r1=feutil('getcg',r2);r1(~isfinite(r2.Elt(:,1)),:)=[];
  i2=RO.Node(RO.DT.nearestNeighbor(r1),4);
  bas=RO.bas(ismember(RO.bas(:,1),i2),:);
  [i3,i4]=ismember(i2,bas(:,1));
  bas(:,1)=(1:size(bas,1))'+RO.mpid(2)*1000;
  mpid=feutil('mpid',model); mpid(cEGI,2)=bas(i4,1);
  model.Elt=feutil('mpid',model,mpid);
  if ~isfield(model,'bas');model.bas=[];end
  model.bas=[model.bas;bas];
  out=model;
 otherwise 
         dbstack; keyboard;
 end
else
 i1=RO.DT.nearestNeighbor(xyz);tol=1e-4;
 out=NaN(size(xyz,1),1);  i3=1; 
 for j1=1:length(i1)
  cur=struct; 
  i2=find(RO.en(:,i1(j1))~=0)-1;
  cur.Adjacent=RO.Adjacent(:,i2);% adjacent
  cur.faces=RO.faces(:,i2);% node
  cur.EltEdge=RO.EltEdge(:,i2);
  for j2=1:size(cur.faces,2)% In any face
   r2=triRS(RO.Node(cur.faces(:,j2),5:7)',xyz(j1,:)');
   r=r2(1);s=r2(2);t=1-r-s;z=r2(3);
   %r=r2(2);s=r2(3);t=1-r-s;z=r2(1);
   if r>-tol&&s>-tol&&1-r-s>-tol % Within triangle
    if ~isfinite(out(j1));out(j1)=z;
    elseif abs(z)<abs(out(j1)); out(j1)=z;
    end
   end      
  end  
  if ~isfinite(out(j1)); % Not in front of any face, seek edges
   i5=unique(cur.EltEdge);
   cur.edge=RO.Edge(:,i5);cur.normal=RO.EdgeNormal(:,i5);
   i3=cur.edge(1,:)~=i1(j1)&cur.edge(2,:)~=i1(j1);
   cur.edge(:,i3)=[];cur.normal(:,i3)=[];
   for j2=1:size(cur.edge,2)% In any edge
     n1=RO.Node(cur.edge(:,j2),5:7);
     x=xyz(j1,:)-n1(1,:); n1=n1(2,:)-n1(1,:); 
     r=x*n1'/(n1(1)*n1(1)+n1(2)*n1(2)+n1(3)*n1(3));
     if r>-tol&&r<1+tol;  % facing the segment
      r2=(x-r*n1)*cur.normal(:,j2);
      if ~isfinite(out(j1));out(j1)=r2;
      elseif abs(r2)<abs(out(j1)); out(j1)=r2;
      end
     end
   end
  end %Edges
  if ~isfinite(out(j1)); % Not in front of any face,edge, seek Nodes   
   r2=(xyz(j1,:)-RO.Node(i1(j1),5:7))*RO.NodeNormal(:,i1(j1));
   if ~isfinite(out(j1));out(j1)=r2;
   elseif abs(r2)<abs(out(j1)); out(j1)=r2;
   end
  end %Nodes
 end% loop on nodes
 out(abs(out)<sp_util('epsl'))=0;
end
%% #findR : generic segment search -4
function out = findR(ns,np);

if size(np,1)==1; % Single point : position with respect to two planes
    % positive interior min abs or dist to other
  dbstack; keyboard;
elseif size(np,1)==2 % ls for two segments
  r2=np(1,:)-ns(1,:); r4=diff(ns); 
  % 1. give plane between ns and np(1,:)
  p=sp_util('basis',r4,r2);r4=norm(r4);
  % 2. find projection of np(2,:) on this plane
  r3=(np(2,:)-ns(1,:))*p(:,1:2);r2=r2*p(:,1:2);
  % 3. find intersection of two segments
  r=r3(2)/(r3(2)-r2(2));%r*r2+(1-r)*r3
  r=(r*r2(1)+(1-r)*r3(1))/r4;
  out=zeros(4,1); out(1)=r; 
  if r>1; out(4)=-r+1;elseif r<0; out(4)=r; 
  else; out(4)=min(r,1-r);
  end
    
else; error('Not implemented');
end
%% #findRS : generic multielement surface search -4
function [out,out1] = findRS(RO,np,RA);

out=zeros(4,length(RO.curMatchE));
for jElt=1:length(RO.curMatchE);
  i1=RO.elt(:,RO.curMatchE(jElt));
  if i1(3)==i1(4); 
      triRS(RO.vert(i1,:))
  else
    ns=RO.vert(i1,:)-np([1 1 1 1],:);
    p1=sp_util('basis',ns'); % Local basis   
    rst=RO.FindRST('quad4',ns*p1); 
    if abs(rst(1))<=1&&abs(rst(2))<=1; % inside
      r2=[1-rst;rst-1];  rst(4)=min(abs(r2)); 
    elseif abs(rst(1))>abs(rst(2))
     if rst(1)<=-1; rst(4)=rst(1)+1; % Facing an r edge
     elseif rst(1)>=1; rst(4)=-rst(1)+1; % Facing an r edge
     end
    else
     if rst(2)<=-1; rst(4)=rst(2)+1; % Facing an s edge
     elseif rst(2)>=1; rst(4)=-rst(2)+1; % Facing an s edge
     end
    end
    rst(3)=RO.q4.N(rst(1:2)')*ns*p1(:,3); 
    out(:,jElt)=rst;    % columns give [r,s,z,dist]
  end
end
if nargin<3
elseif strcmpi(RA.do,'moveinsert')
 %% select the matching element and possibly insert node
 i3=find(out(4,:)>=0);
 if isempty(i3);[r4,i3]=max(out(4,:));end
 if length(i3)>1; [r4,i4]=min(r2(4,i3));i3=i3(i4);end
 if isempty(i3); out=RO;return;end % If not connected to an elt then new node
 RO.curMatchE=RO.curMatchE(i3);
 out=out(:,i3);i1=RO.elt(:,RO.curMatchE);
 if out(1)>RO.tolE*2 % x+ Edge
    if out(2)>RO.tolE*2; dbstack; keyboard;
    elseif out(2)<-RO.tolE*2;dbstack; keyboard;
    else; dbstack; keyboard;
    end
 elseif out(1)<-RO.tolE*2 %x- edge
   dbstack; keyboard;
 else; % x center band
    if out(2)>RO.tolE*2; dbstack; keyboard;
    elseif out(2)<-RO.tolE*2;dbstack; keyboard;
    else;st=sprintf('f%i.%i.%i.%i',sort(i1));i2=size(RO.vert,1)+1;
    end
 end
 RO.vert(i2,:)=np; RO.nmap(st)=i2; % Name the node
 RO.incuts(end+1)=i2;
 out=RO; 
elseif strcmpi(RA.do,'match')
 % select the matching element
 i3=find(out(4,:)>=0);
 if isempty(i3);[r4,i3]=max(out(4,:));end
 if length(i3)>1; [r4,i4]=min(r2(4,i3));i3=i3(i4);end
 out=out(:,i3);
 RO.curMatchE=RO.curMatchE(i3); out=RO; 
end

%% #triRS FindRST_tria Find the r,s -4
function out = triRS(ns,np);

%% safe coordinate search in triangle
%figure(1);plot(ns(1,[1:3 1]),ns(2,[1:3 1]),'-',np(1),np(2),'o');text(ns(1),ns(2),'1');text(ns(1,2),ns(2,2),'2');
np=np-ns(:,1); % vector from node 1
xe=ns(:,2)-ns(:,1); %lx2=(xe(1)*xe(1)+xe(2)*xe(2)+xe(3)*xe(3));

ye=ns(:,3)-ns(:,1); %ly2=(ye(1)*ye(1)+ye(2)*ye(2)+ye(3)*ye(3));
ze=-[ye(2)*xe(3)-ye(3)*xe(2);ye(3)*xe(1)-ye(1)*xe(3);ye(1)*xe(2)-ye(2)*xe(1)];
lz=sqrt(ze(1)*ze(1)+ze(2)*ze(2)+ze(3)*ze(3));
ze=ze/lz; 

out=[xe ye ze]\np;
%pos=[xe ye ze]\np;out=[pos([3 1 2]);0;0;ze;0];%[xp,yp,zp,g,r s ielt wjdet n3x n3y n3z]

%% #quadRS FindRST_tria Find the r,s -4
function out = quadRS(ns,np);

dbstack; keyboard;
%% safe coordinate search in triangle
%figure(1);plot(ns(1,[1:3 1]),ns(2,[1:3 1]),'-',np(1),np(2),'o');text(ns(1),ns(2),'1');text(ns(1,2),ns(2,2),'2');
np=np-ns(:,1); % vector from node 1
xe=ns(:,2)-ns(:,1); %lx2=(xe(1)*xe(1)+xe(2)*xe(2)+xe(3)*xe(3));

ye=ns(:,3)-ns(:,1); %ly2=(ye(1)*ye(1)+ye(2)*ye(2)+ye(3)*ye(3));
ze=-[ye(2)*xe(3)-ye(3)*xe(2);ye(3)*xe(1)-ye(1)*xe(3);ye(1)*xe(2)-ye(2)*xe(1)];
lz=sqrt(ze(1)*ze(1)+ze(2)*ze(2)+ze(3)*ze(3));
ze=ze/lz; 

out=[xe ye ze]\np;
%pos=[xe ye ze]\np;out=[pos([3 1 2]);0;0;ze;0];%[xp,yp,zp,g,r s ielt wjdet n3x n3y n3z]



%% #dToInterp : distance to interpolant -3
function out=dToInterp(xyz,RO,model) %#ok<DEFNU>

if ischar(xyz)&&strcmpi(xyz,'init')
  def=stack_get(model,'curve','DInfo',3);
  R1=RO;
  lim1=[min(model.Node(:,5:7));max(model.Node(:,5:7))];
  out=struct('shape','interp','distFcn',[], ...
      'interp',[],'box',lim1,'mpid',[model.pl(1) model.il(1)]);
  if any(strcmp(R1.dtype,{'nearest';'linear';'natural'}))
   out.interp=scatteredInterpolant(model.Node(:,5:7),def.def(:,1),R1.dtype,'none'); 
   out.distFcn=@(xyz)ddToInterp(xyz,out);
  elseif strcmp(R1.dtype,'cnem')   
   Interpol=m_cnem3d_interpol(model.Node(:,5:7),[]);
   out.interp=@(x,y,z)feval(lsutil('@cnem'),Interpol,def.def(:,1),x,y,z);     
   out.distFcn=@(xyz)ddToCnem(xyz,out);
  else;warning('%unknown method %s',R1.dtype);dbstack;keyboard
  end  
  out.sel=feutil('getpatchnew',model); 
  out.sel=sdsetprop(out.sel,'f1Prop','visible','off');
  out.sel=sdsetprop(out.sel,'f2Prop','visible','off');
  return
end

function out=ddToInterp(xyz,R1)
   %% #ddToInterp: interpolation from scattered -4
   
   if isfield(R1,'box')
    nidin=xyz(:,1)>R1.box(1,1) & xyz(:,1)<R1.box(2,1) & ...
     xyz(:,2)>R1.box(1,2) & xyz(:,2)<R1.box(2,2) & xyz(:,3)>R1.box(1,3) & xyz(:,3)<R1.box(2,3);
    d=NaN*ones(size(xyz,1),1);
    d1=R1.interp(xyz(nidin,1),xyz(nidin,2),xyz(nidin,3));d(nidin)=d1;
   else
    d=R1.interp(xyz(:,1),xyz(:,2),xyz(:,3));
   end
   out=d;

function out=dToCnem(xyz,RO)
   %% #dToCnem Interp: interpolation from scattered -3
 dbstack; keyboard;

%% #CNEM -4
function out=cnem(Interpol,val,x,y,z)
Interpol=Interpol.set_point([x,y,z],0);
out=Interpol.interpolate(val);
 
%% #dToPoly : delaunay defined distance -3
function [out,out1]=dToPoly(xyz,RO,model)
% RO.distFcn=@(xyz)distFromPoly(DT,xyz)
Cam='';
if isstruct(xyz)
 r2=xyz; Cam=fieldnames(xyz);xyz=r2.(Cam{1});Cam=Cam{1};
end
if ischar(xyz)
 if comstr(xyz,'init')   
  %% Create the levelset
  mo1=RO; 
  % default single line
  if ~isfield(mo1,'Elt');mo1.Elt=feutil('objectbeamline',mo1.Node(:,1)');end
  r1=fe_gmsh('lineloops',feutil('selelt seledge',mo1));r1=r1{1};
  NNode=sparse(mo1.Node(:,1),1,1:size(mo1.Node,1));r1=full(NNode(r1));
  mo1.isClosed=r1(1)==r1(end);
  n1=mo1.Node(r1,5:7);
  if ~isfield(mo1,'Type');RO.Type='3D';mo1.Type='3D';
  elseif strcmpi(RO.Type,'RevLine')
   p=sp_util('basis',RO.axis,[0 0 0]);mo1.orig=mo1.orig(:)';
   n1=n1-repmat(RO.orig,size(n1,1),1);
   n1=[n1*p(:,1) sqrt(sum((n1*p(:,2:3)).^2,2))];n1(:,3)=0;%z,r coords
  elseif strcmpi(RO.Type,'ExtLine')
   p=sp_util('basis',mo1.axis,[0 0 0]);p=p(:,[2:3 1]); p(:,3)=0; 
   n1=n1*p; 
   mo1.p=p(:,1:2); % Keep for later ls evaluation
  end
  r2=diff(n1);
  mo1.L=sqrt(sum(r2.^2,2));r2=diag(sparse(1./mo1.L))*r2;
  mo1.xe=r2;
  if mo1.isClosed;mo1.xe(end+1,:)=mo1.xe(1,:);mo1.L(end+1,:)=mo1.L(1,:);
  else;mo1.xe(end+1,:)=mo1.xe(end,:);mo1.L(end+1,:)=mo1.L(end,:);
  end
  
  if isfield(RO,'normal');% normal to surface
   mo1.ze=repmat(mo1.normal(:)',size(mo1.xe,1),1);
  else % Attempt to build normal field
   mo1.ze=ones(size(mo1.xe,1),3); 
   for j1=1:size(r2,1);
    j2=j1;p=0;
    while all(p==0);
     j2=j2+1; if j2>size(r2,1);j2=1;end
     p=r2(j2,:)-r2(j1,:);
     if all(abs(p)<1e-8);p=0;continue;end
     p=p/norm(p);
    end
    r1=sp_util('basis',r2(j1,:),p);
    if j1>1&&mo1.ze(j1-1,:)*r1(:,3)<0; mo1.ze(j1,:)=-r1(:,3)';
    else; %if min(cross(r2(j1,:)',p))<0;dbstack; keyboard;end
     mo1.ze(j1,:)=r1(:,3)';
    end
    %if r1(end)<0;keyboard;end
   end
   if mo1.isClosed;mo1.ze(end,:)=mo1.ze(1,:);% last point is same as first
   else; mo1.ze(end,:)=mo1.ze(end-1,:);end
  end
  %fecom('showmap',struct('vertex',mo1.Node(:,5:7),'normal',mo1.ze))
  mo1.dir=cross(mo1.xe,mo1.ze);% normal to local line plane
  if mo1.isClosed
   mo1.dirn=mo1.dir(1:end-1,:)+mo1.dir([end 1:end-2],:);
   mo1.dirn=diag(sparse(1./sqrt(sum(mo1.dirn.^2,2))))*mo1.dirn;
   n1(end,:)=[];
  else
   r1=sparse([1;1]*(1:size(n1,1)),[1:length(mo1.L) length(mo1.L);1 1:length(mo1.L)], ...
       [mo1.L' mo1.L(end);mo1.L(1) mo1.L']);
   mo1.dirn=diag(1./sum(r1',1))*r1*mo1.dir(1:end-1,:); %#ok<UDIM>
   mo1.L(end+1)=mo1.L(end);
  end
  if any(strcmpi(mo1.Type,{'RevLine','ExtLine'}));
      n1(:,3)=[];mo1.xe(:,3)=[];mo1.dirn(:,3)=[]; mo1.dir(:,3)=[];
  end
  mo1.DT = delaunayTriangulation([n1;1e4*ones(1,size(n1,2))]);
  mo1.distFcn=@(xyz)dToPoly(xyz,mo1);
  out=mo1;
 elseif strncmpi(xyz,'ro',2);out=RO;
 else; error('%s Not valid',xyz);
 end
 if 1==2
  figure(10);clf;
  n1=out.DT.Points(1:end-1,:);d1=out.dir(1:end-1,:)*5;%d1=out.xe(1:end-1,:)*5;
  n1=out.DT.Points(1:end-1,:);n1=(n1(1:end-1,:)+n1(2:end,:))/2;d1=out.dir(1:end-2,:)*5;%segorient
  mo1.dirn=mo1.dir(1:end-1,:)+mo1.dir([end 1:end-2],:);
  mo1.dirn=diag(sparse(1./sqrt(sum(mo1.dirn.^2,2))))*mo1.dirn;
  n1=out.DT.Points(1:end-1,:);d1=mo1.dirn*5;%node
  figure(10);clf;
  plot3(n1(:,1),n1(:,2),n1(:,3),'+', ...
      n1(:,1)+d1(:,1),n1(:,2)+d1(:,2),n1(:,3)+d1(:,3),'o')
  set(gca,'dataaspectratio',[1 1 1]);view(2);
 end
 return
    
elseif strcmpi(RO.Type,'RevLine')
 %% #dToPoly.Revline : distance to flat line with revolution -3
 p=sp_util('basis',RO.axis,[0 0 0]);
 xyz=xyz-repmat(RO.orig,size(xyz,1),1);
 xyz=[xyz*p(:,1) sqrt(sum((xyz*p(:,2:3)).^2,2))];

elseif strcmpi(RO.Type,'ExtLine')
 %% #dToPoly.Extline : distance to extruded flat line -3
 a=xyz;
 xyz=xyz*RO.p; % Project to plane
 
else;xyz=xyz(:,1:size(RO.DT.Points,2)); 
end
[i1,out1]=nearestNeighbor(RO.DT,xyz);
%r2=RO.DT.Points; line((r2(:,[1 1 1])+RO.orient(:,1)*[0 1 Inf])',(r2(:,[1 1 1]*2)+RO.orient(:,2)*[0 1 Inf])',(r2(:,[1 1 1]*3)+RO.orient(:,3)*[0 1 Inf])')
%r2=RO.DT.Points; line((r2(:,[1 1 1])+RO.normal(:,1)*[0 1 Inf])',(r2(:,[1 1 1]*2)+RO.normal(:,2)*[0 1 Inf])',(r2(:,[1 1 1]*3)+RO.normal(:,3)*[0 1 Inf])')
%  < (point - nearest) ^ tangent , normal >
n1=RO.DT.Points(i1,:)-xyz;

out=NaN(size(n1,1),1); i0=i1; % normal at nodes
%z1=sum(n1.*RO.ze(i1,:),2);% vertical distance to plane at nearest node
for j2=1:2
 if j2==1; o1=-RO.xe(i1,:);else; o1=RO.xe(i1,:);end
 % fecom('shownodemark',[RO.DT.Points([i1;i1+1],:);xyz],'marker','o','color','r')
 r=sum(n1.*o1,2)./RO.L(i1);%local coord first seg
 i3=find(r>=0&r<=1);% edge dist
 d3=sum(RO.dir(i1(i3),:).*n1(i3,:),2); 
 i4=abs(d3)<abs(out(i3))|~isfinite(out(i3));
 out(i3(i4))=d3(i4);out1(i3(i4),2)=i1(i3(i4));
 i1=i1-1; i1(i1<=0)=length(RO.L);% Do second segment 
end
i4=~isfinite(out); out(i4)=sum(RO.dirn(i0(i4),:).*n1(i4,:),2);

if isfield(RO,'zlim')&&~RO.planar % Introduce
 dbstack; keyboard;
 out(z1<RO.zlim(1)|z1>RO.zlim(2))=NaN; 
elseif ~isempty(Cam)&&strcmpi(Cam,'stick')
 out=(diag(sparse(out))*RO.dirn(out1(:,2),:)+xyz)*RO.p'+(a*RO.axis(:))*RO.axis(:)';
 return;
end
out(abs(out)<sp_util('epsl'))=0;
% z=RO.DT.Points; line(z(:,1),z(:,2),z(:,3),'marker','o')
%% #defLevelList :  -3 
function  out=defLevelList(out,LevelList);
r2=[min(out) max(out)]*[1.01 -.01;-.01 1.01];% range
r3=(LevelList-r2(1))/diff(r2); r3(r3<0|r3>1)=[];
if isempty(r3);out=ones(size(out));
else
 r4=(out-r2(1))/diff(r2); out=(r4-r3(1));
 for j1=2:length(r3);out=out.*(r4-r3(j1));end
end

%% #perm_nid
function [elt,sevLS,listPerm]=perm_nid(ind,elt,sevLS,allcase);
listPerm=[];
for i1=1:length(allcase)
 if isempty(allcase{i1});continue;end
 in1=find(ismember(sevLS,[allcase{i1};-1*allcase{i1}],'rows'));
 if ~isempty(in1)
  listPerm=[listPerm;in1 i1*ones(size(in1))];
  elt(in1,ind)=circshift(elt(in1,ind),[0 i1]); %i1,2);
  sevLS(in1,:)=circshift(sevLS(in1,:),[0 i1]); %i1,2);
 end
end

function [va,vb]=checkTol(va,vb);

if nargin==1
 aa=abs(va);ma=max(aa)*1e-10;
 va(aa<ma)=0;
else
 ma=max(abs(va),abs(vb))*1e-10;
 if abs(va)<ma;va=0;end
 if abs(vb)<ma;vb=0;end
end
function  [nc,vc,contour]=onContour(nc,vc,contour,dis);
if ~isempty(contour);
 [r1,i1]=min(sum((contour-repmat(nc,size(contour,1),1)).^2,2));
 %[r1 norm(na-nb)*RO.tolE]
 if sqrt(r1)<dis; % Use a countour point if close enough
  nc=contour(i1,:);vc=0; contour(i1,:)=[];
  %li.distFcn(nc);
 else;
  1;
 end
end

%% #newOnLs : add internal node but force on LS 
function [mo2,ind]=newOnLS(x,y,z,elti,ndir,mo2,def);
 r1=reshape(def.def(elti),size(elti));
 xc=[sum(x(elti),ndir) sum(y(elti),ndir) sum(z(elti),ndir)]/8;
 for j1=1:size(xc,1) % Seek to place xc on LS
  [va,i2]=sort(abs(r1(j1,:)),'descend');
  na=[x(elti(j1,i2)) y(elti(j1,i2)) z(elti(j1,i2))];vd=NaN; 
  % xxx should actually then select element quality
  for j2=1:size(na,1)
   [nc,vc]=segSearch(na(j2,:),xc(j1,:),def,0,1.5);
   if vc==0; nd=nc; vd=0; break;
   elseif abs(vc)<abs(vd); nd=vc; vd=vc; 
   end
  end
  if vd~=0&&~isempty(nd);dbstack; keyboard;end
  xc(j1,:)=nd;
 end
 [mo2.Node,ind]=feutil('AddNode',mo2.Node,xc); %KnownNew
 ind=mo2.Node(ind,1);% gid];

%% #segSearch robust placement of a new node on an edge
function [nc,vc]=segSearch(na,nb,def,r0,r1);

va=def.distFcn(na); vb=def.distFcn(nb); nb0=nb;vb0=vb;
d=(nb-na)';d=d/(d(1)*d(1)+d(2)*d(2)+d(3)*d(3));
nc=na; vc=va; ite=4;% start point arbitrary
while 1
  if vc==vb; 
    figure(10);clf; plot(r,vd); 
    error('Inconsistent variation of distance'); 
  end
  r=vc/(vc-vb); nd=nb*r+nc*(1-r);vd=def.distFcn(nd);
  r=(nd-na)*d;
  if r>r0&&r<r1; 
       if vd==0; nc=nd;vc=0;break;
       elseif abs(vd)<abs(vc); nb=nc;vb=vc; nc=nd;vc=vd; 
       elseif sign(vc)~=sign(vd); % went over LS
         r=linspace(0,1,50)';
         nd=r*nc+(1-r)*nd;vd=def.distFcn(nd);
         [un1,i1]=sort(abs(vd));nc=nd(i1(1),:);vc=vd(i1(1));
         nb=nd(i1(2),:);vb=vd(i1(2));
       elseif sign(va)~=sign(vb0);% Something went wrong restart
         r=linspace(0,1,100)';
         nd=r*na+(1-r)*nb0;vd=def.distFcn(nd);
         [un1,i1]=sort(abs(vd));nc=nd(i1(1),:);vc=vd(i1(1));
         nb=nd(i1(2),:);vb=vd(i1(2));
       else; 
           dbstack; keyboard;
       end
       if vd==0; break;end % on ls
  elseif ite&&sign(va)~=sign(vb0) % divide segment
         r=linspace(0,1,100)';
         nd=r*na+(1-r)*nb0;vd=def.distFcn(nd);
         [un1,i1]=sort(abs(vd));nc=nd(i1(1),:);vc=vd(i1(1));
         nb=nd(i1(2),:);vb=vd(i1(2));
         ite=0;   
  else;%if r>=r1;
     nc=[];vc=NaN;   break; % out of segment
  end
  if vb==0; nc=nb; vc=vb;break;end
end
    
 
%% #cleanNew robust placement of a new node on an edge
function [xnew,li]=cleanNew(model,RO,li,i1);

% call 1 for new nodes
if ~isfield(RO,'ind2')
 n1=model.Node(RO.edges(i1,1),5:7);
 n2=model.Node(RO.edges(i1,2),5:7);r=RO.r(i1);
else % Call 2 for tol move nodes
 n1=model.Node(RO.ind1,5:7);n2=model.Node(RO.ind2,5:7);
 r=RO.r;
end
if ~isfield(RO,'tolVc');RO.tolVc=1e-5;end
xnew= (r*[1 1 1]).*n1+((1-r)*[1 1 1]).*n2;
if iscell(li)&&length(li)==1&&isfield(li{1},'distFcn')
 li=li{1};
elseif ~isfield(li,'distFcn');return
end
val=li.distFcn(xnew);
if ~isfield(li,'contour');contour=[];
elseif size(li.contour,2)==7;contour=li.contour(:,5:7);
else;contour=li.contour;
end
for j1=find(val==0)' % nearest nodes for accurate contour
 nc=xnew(j1,:);na=n1(j1,:);nb=n2(j1,:);
 if ~isempty(contour);
  [r1,i1]=min(sum((contour-repmat(nc,size(contour,1),1)).^2,2));
  %[r1 norm(na-nb)*RO.tolE]
  if sqrt(r1)<norm(na-nb)*RO.tolE; % Use a countour point if close enough
   xnew(j1,:)=contour(i1,:);vc=0; val(j1)=0; contour(i1,:)=[];
   %li.distFcn(nc);
  end
 end
end
for j1=find(val~=0)' % Move nodes were edge lininterp incorrect
  na=n1(j1,:);nb=n2(j1,:); 
 [nc,vc]=segSearch(na,nb,li,0,1);
 [nc,vc,contour]=onContour(nc,vc,contour,norm(n1(j1,:)-n2(j1,:))*RO.tolE);
 if vc~=0;dbstack; keyboard;end
 xnew(j1,:)=nc;val(j1)=vc;
end
if nargout>0;li.contour=contour;end
% fecom('shownodemark',xnew,'marker','o')
if isfield(li,'def')&&isfield(RO,'i0');% Set level of moved node to 0
 li.def(RO.i0)=0;
end
%
1;
%cf=feplot;axes(cf.ga);lin=line(xnew(:,1),xnew(:,2),xnew(:,3),'marker','o','linestyle','none')
%lsutil('ViewLs',model,li);p=patch('vertices',model.Node(:,5:7),'faces',RO.edges(i1,:),'edgecolor','r','linewidth',2)
% Should not show edges with


%% #tolMoveNode : tolerance/set node moving
function  [model,RO,def]=tolMoveNode(model,RO,def)
% ls should never be reevaluated at exit of this function

if ~isfield(RO,'set0');RO.set0=[];
else; def.def(ismember(fix(def.DOF),RO.set0))=0;
end
[i1,i2]=ismember(fix(def.DOF),model.Node(:,1));%sparse([model.Node(:,1)-fix(def.DOF)])
def.def=def.def(i2);def.DOF=def.DOF(i2);
i1=RO.edges; NNode=RO.NNode; % .edge NodeIndex based here
% Possibly move some nodes to the edges
[i1,i2]=unique(sort(i1,2),'rows'); %remove duplicate edges & order them
c1=sparse(1:size(i1,1),i1(:,1),1,size(i1,1),size(def.DOF,1));
c2=sparse(1:size(i1,1),i1(:,2),1,size(i1,1),size(def.DOF,1));
r1=c1*def.def;r2=c2*def.def;
i3=r2==r1&r1;
if any(i3); % Fully flat edge
 [r1(i3) r2(i3)]
 if sdtdef('isinteractive')
  %dbstack; keyboard;
  sdtw('_ewt','lsutil tolMoveNode specific case to handle')
 else; sdtw('_nb','lsutil tolMoveNode specific case to handle')
 end
end
r=r2./(r2-r1); % r=1 level set on node 1
% Sort by increasing proximity to edge
%[un1,i3]=sort(min([r 1-r],[],2));r1=r1(i3);r2=r2(i3);r=r(i3);i1=i1(i3,:);
% Long term the exact position of moved node on edge may be needed
in2=find(abs(r)<RO.tolE); % r close to 0 LS node 2
if isfield(RO,'onSurf'); % Nodes can only be moved on the surface
 i3=any(reshape(~ismember(model.Node(i1(in2,1:2)),RO.onSurf),[],2),2);
 in2(i3)=[];
end
if isfield(RO,'Fixed');
 % Some nodes may be forced to be fixed
 if ~isfield(RO,'tolF');RO.tolF=1e-3;end % tolF set LS to 0
 i3=ismember(model.Node(i1(in2,2)),RO.Fixed);
 i4=i1(in2(r(in2(i3))<RO.tolF),2);def.def(i4)=0;% On edge set=0;
 RO.set0=unique([RO.set0;fix(def.DOF(i4))]);
 in2(i3)=[];
end
[r3,i3]=sort(r(in2));in2=in2(i3);% Use closest to end
ind2=i1(in2,2);ind1=i1(in2,1);% Left right node indices
[ind2,i3]=unique(ind2,'first');ind1=ind1(i3);r3=r3(i3);
if ~isfield(RO,'tolVc');RO.tolVc=1e-4;end
RC=struct('ind1',ind1,'ind2',ind2,'r',r3,'i0',ind2,'tolE',RO.tolE,'tolVc',RO.tolVc);
[model.Node(ind2,5:7),def]=cleanNew(model,RC,def);

r1=c1*def.def; r2=c2*def.def;r=r2./(r2-r1);
% r on right ege between [1-tol2 - 1-tol1 ] % move to right
in1=abs(r)>1-RO.tolE;in1=find(in1);
if isfield(RO,'onSurf'); % Nodes can only be moved on the surface
 i3=any(reshape(~ismember(model.Node(i1(in1,1:2)),RO.onSurf),[],2),2);
 in1(i3)=[];
end
if isfield(RO,'Fixed')&&~isempty(RO.Fixed); % Some nodes may be forced to be fixed
 i3=ismember(model.Node(i1(in1,1)),RO.Fixed);
 i4=i1(in1(1-r(in1(i3))<RO.tolF),1);def.def(i4)=0;% On edge set=0;
 RO.set0=unique([RO.set0;fix(def.DOF(i4))]);
 in1(i3)=[];
end
[r3,i3]=sort(r(in1),'descend');in1=in1(i3);% Use closest to end
ind2=i1(in1,2);ind1=i1(in1,1);% Left right node indices
[ind1,i3]=unique(ind1,'first');ind2=ind2(i3);r3=r3(i3);

RC=struct('ind1',ind1,'ind2',ind2,'r',r3,'i0',ind1,'tolE',RO.tolE,'tolVc',RO.tolVc);
[model.Node(ind1,5:7),def]=cleanNew(model,RC,def);
r1=c1*def.def; r2=c2*def.def;r=r2./(r2-r1);
if isfield(RO,'onSurf') % Do not move surface nodes
 RO.edges=i1;
 i2=~ismember(model.Node(RO.edges(:,1),1),RO.onSurf);r1(i2)=0;
 i2=~ismember(model.Node(RO.edges(:,2),1),RO.onSurf);r2(i2)=0;
 r=r2./(r2-r1); r(~isfinite(r))=0;
 def.def(~ismember(fix(def.DOF),RO.onSurf),:)=0;
end
RO.edges=i1; RO.values=[r1 r2]; RO.r=r;
if isfield(def,'contour')&&~isempty(def.contour)
 % Fprce move to nearest node on unused contour point
 i1=norm(RO.values,'inf');
 if i1==0; i1=norm(max(def.contour)-min(def.contour),'inf');end
 i1=feutilb(sprintf('addnode-nearest epsl %g',i1), ...
  model.Node(:,5:7),def.contour);
if any(i1) % xxEM
 model.Node(i1,5:7)=def.contour;def.def(i1,:)=0;def.contour=[];
 [i2,i3]=ismember(model.Node(i1,1),RO.edges);RO.values(i3(i2))=0;
end
end
1;%[nnz(def.def==0) nnz(def.distFcn(model.Node(:,5:7))==0)]

%% #charLs : strings representing level set signs
function lsC=charLs(def,elt);

if nargin==2
 lsC=repmat('0',size(elt,1),size(elt,2));
 r2=def.def(elt);lsC(r2>0)='+';lsC(r2<0)='-';
else
 lsC=repmat('0',size(def,1),size(def,2));
 lsC(def>0)='+';lsC(def<0)='-';
end
lsC=cellstr(lsC);
%% #safeAddElt
function [mo2,RO]=safeAddElt(mo2,RO,ElemP,elt);

if isempty(ElemP)||isempty(elt); return;end
i1=feval(ElemP,'prop');
i2=RO.nextId+(0:size(elt,1)-1)';RO.nextId=RO.nextId+size(elt,1);
if isfield(RO,'EltOrient')&&size(elt,2)>=i1(3)
 eltid=elt(:,i1(3)); % Original eltid
 if all(eltid)
  nind=sparse(RO.EltOrient.EltId,1,1:size(RO.EltOrient.EltId,1));
  i3=max(eltid);if i3>length(nind);nind(i3)=0;end
  i3=full(nind(eltid));
  if ~isfield(RO.EltOrient,'BasId');
   RO.EltOrient.BasId=RO.EltOrient.basid;
   RO.EltOrient=rmfield(RO.EltOrient,'basid');
  end
  i4=(i3==0); % No orient add new identity
  if any(i4);
   [i3,i4]=unique(eltid(i4));
   RO.EltOrient.BasId=[RO.EltOrient.BasId;i3];
   RO.EltOrient.EltId=[RO.EltOrient.EltId;i3];
   r2=i3;r2(:,[2 7 11 15])=1;
   RO.EltOrient.bas=[RO.EltOrient.bas;r2(:,1:size(RO.EltOrient.bas,2))];
   nind=sparse(RO.EltOrient.EltId,1,1:size(RO.EltOrient.EltId,1));
   i3=full(nind(eltid));
  end
  RO.EltOrient.EltId=[RO.EltOrient.EltId;i2];
  RO.EltOrient.BasId=[RO.EltOrient.BasId;i2];
  RO.EltOrient.bas=[RO.EltOrient.bas;RO.EltOrient.bas(i3,:)];
 end
end

if size(elt,2)<i1(3); elt(end,i1(3))=0; end
if ~isfield(RO,'EEIdx'); RO.EEIdx=[elt(:,i1(3)) i2];
else;
 [i3,i4]=ismember(elt(:,i1(3)),RO.EEIdx(:,2));
 if any(i3) % multiple elt transform requires update from orginial EltId
  r5=[elt(:,i1(3)) i2];
  r5(i3,1)=RO.EEIdx(i4(i3),1);
  RO.EEIdx=[RO.EEIdx;r5];
 else; RO.EEIdx=[RO.EEIdx;elt(:,i1(3)) i2];
 end
end
elt(:,i1(3))=i2;
mo2=feutil('AddElt',mo2,ElemP,elt);
% r1=[Inf abs(ElemP)];
% mo2.Elt(end+1,1:length(r1))=r1;
% mo2.Elt(end+(1:size(elt,1)),1:size(elt,2))=elt;

%% #ListGenericCuts re(Shape) -2
%% #reTria : generic cut of triangular elements -3
function RE=reTria
allcases={[0 -1 1;0 1 -1],[2 3],'[curelt(:,[1 2 4]) curprop;curelt(:,[1 4 3]) curprop]',[]
 [-1 0 1;1 0 -1],[1 3],'[curelt(:,[1 2 4]) curprop;curelt(:,[2 3 4]) curprop]',[]
 [-1 1 0;1 -1 0],[1 2],'[curelt(:,[1 4 3]) curprop;curelt(:,[2 3 4]) curprop]',[]
 [-1 -1 1;1 1 -1],[1 3;2 3],'[curelt(:,[3 4 5]) curprop]','[curelt(:,[1 2 5 4]) curprop]'
 [-1 1 -1;1 -1 1],[1 2;2 3],'[curelt(:,[2 5 4]) curprop]','[curelt(:,[1 4 5 3]) curprop]'
 [1 -1 -1;-1 1 1],[1 2;1 3],'[curelt(:,[1 4 5]) curprop]','[curelt(:,[2 3 5 4]) curprop]'
 };
RE=struct('cases',vertcat(allcases{:,1}), ...
 'cindex',reshape([1;1]*(1:size(allcases,1)),[],1));
RE.newelt=[{'tria','quad'};allcases(:,3:4)];
RE.CutEdges=allcases(:,2);RE.ElemP='tria3';
if nargout==0; genericDisplay(RE);end % feval(lsutil('@reTria'))

%% #reQuad : generic cut of quadrangle elements -3
function RE=reQuad
allcases=cell(27,4);
%LS on 2 nodes + 1 cut
allcases(1,:)={[1 -1 0 0;-1 1 0 0],[1 2],'[curelt(:,[1 5 4]) curprop;curelt(:,[5 2 3]) curprop;curelt(:,[5 3 4]) curprop]',[]};
allcases(2,:)={[1 -1 0 1;-1 1 0 -1],[1 2],'[curelt(:,[5 2 3]) curprop]','[curelt(:,[1 5 3 4]) curprop]'};
allcases(3,:)={[1 -1 -1 0;-1 1 1 0],[1 2],'[curelt(:,[1 5 4]) curprop]','[curelt(:,[5 2 3 4]) curprop]'};
%
allcases(4,:)={[0 1 -1 0;0 -1 1 0],[2 3],'[curelt(:,[1 2 5]) curprop;curelt(:,[5 3 4]) curprop;curelt(:,[1 5 4]) curprop]',[]};
allcases(5,:)={[0 1 -1 -1;0 -1 1 1],[2 3],'[curelt(:,[1 2 5]) curprop]','[curelt(:,[1 5 3 4]) curprop]'};
allcases(6,:)={[1 1 -1 0;-1 -1 1 0 ],[2 3],'[curelt(:,[5 3 4]) curprop]','[curelt(:,[1 2 5 4]) curprop]'};
%
allcases(7,:)={[0 0 1 -1;0 0 -1 1],[3 4],'[curelt(:,[1 2 5]) curprop;curelt(:,[2 3 5]) curprop;curelt(:,[1 5 4]) curprop]',[]};
allcases(8,:)={[1 0 -1 1;-1 0 1 -1 ],[3 4],'[curelt(:,[2 3 5]) curprop]','[curelt(:,[1 2 5 4]) curprop]'};
allcases(9,:)={[0 1 1 -1;0 -1 -1 1],[3 4],'[curelt(:,[1 5 4]) curprop]','[curelt(:,[1 2 3 5]) curprop]'};
%
allcases(10,:)={ [1 0 0 -1;-1 0 0 1],[4 1],'[curelt(:,[1 2 5]) curprop;curelt(:,[2 3 5]) curprop;curelt(:,[5 3 4]) curprop]',[]};
allcases(11,:)={ [1 0 -1 -1;-1 0 1 1],[4 1],'[curelt(:,[1 2 5]) curprop]','[curelt(:,[2 3 4 5]) curprop]'};
allcases(12,:)={ [1 1 0 -1;-1 -1 0 1 ],[4 1],'[curelt(:,[5 3 4]) curprop]','[curelt(:,[1 2 3 5]) curprop]'};
%
allcases(13,:)={[ 1 -1 1 0;-1 1 -1 0],[1 2;2 3],'[curelt(:,[1 5 4]) curprop;curelt(:,[5 2 6]) curprop;curelt(:,[6 3 4]) curprop;curelt(:,[5 6 4]) curprop]',[]};
allcases(14,:)={[ 1 -1 1 1;-1 1 -1 -1],[1 2;2 3],'[curelt(:,[1 5 4]) curprop;curelt(:,[5 2 6]) curprop;curelt(:,[6 3 4]) curprop;curelt(:,[5 6 4]) curprop]',[]};
%
allcases(15,:)={[0 1 -1 1;0 -1 1 -1],[2 3;3 4],'[curelt(:,[1 2 5]) curprop;curelt(:,[5 3 6]) curprop;curelt(:,[6 4 1]) curprop;curelt(:,[1 5 6]) curprop]',[]};
allcases(16,:)={ [1 1 -1 1;-1 -1 1 -1 ],[2 3;3 4],'[curelt(:,[1 2 5]) curprop;curelt(:,[5 3 6]) curprop;curelt(:,[6 4 1]) curprop;curelt(:,[1 5 6]) curprop]',[]};
%
allcases(17,:)={[1 0 1 -1;-1 0 -1 1],[3 4;1 4],'[curelt(:,[1 2 6]) curprop;curelt(:,[2 3 5]) curprop;curelt(:,[5 4 6]) curprop;curelt(:,[2 5 6]) curprop]',[]};
allcases(18,:)={[1 1 1 -1;-1 -1 -1 1],[3 4;1 4],'[curelt(:,[1 2 6]) curprop;curelt(:,[2 3 5]) curprop;curelt(:,[5 4 6]) curprop;curelt(:,[2 5 6]) curprop]',[]};
%
allcases(19,:)={[1 -1  0 -1;-1 1 0 1],[1 4;1 2],'[curelt(:,[1 6 5]) curprop;curelt(:,[6 2 3]) curprop;curelt(:,[3 4 5]) curprop;curelt(:,[5 6 3]) curprop]',[]};
allcases(20,:)={[1 -1 -1 -1;-1 1 1 1],[1 4;1 2],'[curelt(:,[1 6 5]) curprop;curelt(:,[6 2 3]) curprop;curelt(:,[3 4 5]) curprop;curelt(:,[5 6 3]) curprop]',[]};
%
allcases(21,:)={[1 -1 -1 1;-1 1 1 -1],[1 2;3 4],[],'[curelt(:,[1 5 6 4]) curprop;curelt(:,[5 2 3 6]) curprop]'};
allcases(22,:)={[1 1 -1 -1;-1 -1 1 1],[2 3;1 4],[],'[curelt(:,[1 2 5 6]) curprop;curelt(:,[6 5 3 4]) curprop]'};
%
allcases(23,:)={[1 -1 1 -1;-1 1 -1 1],[1 2;2 3;3 4;1 4],'[curelt(:,[1 5 8]) curprop;curelt(:,[5 2 6]) curprop;curelt(:,[6 3 7]) curprop;curelt(:,[7 4 8]) curprop]','[curelt(:,[5 6 7 8]) curprop]'};
% 0 is on edge and -1 outside, a cut is needed (corner case)
% decision to cut should actually be external
%allcases(24,:)={[1 0 -1 0;-1 0 1 0;-1 0 0 0;0 0 -1 0],[],'[curelt(:,[1 2 4]) curprop;curelt(:,[2 3 4]) curprop]',[]};
%allcases(25,:)={[0 1 0 -1;0 -1 0 1;0 0 0 -1;0 0 0 -1],[],'[curelt(:,[1 2 3]) curprop;curelt(:,[1 3 4]) curprop]',[]};
allcases(24,:)={[1 0 -1 0;-1 0 1 0],[],'[curelt(:,[1 2 4]) curprop;curelt(:,[2 3 4]) curprop]',[]};
allcases(25,:)={[0 1 0 -1;0 -1 0 1],[],'[curelt(:,[1 2 3]) curprop;curelt(:,[1 3 4]) curprop]',[]};
allcases(26,:)={[-1 0 0 0;0 0 -1 0],[],'[curelt(:,[1 2 4]) curprop;curelt(:,[2 3 4]) curprop]',[]};
allcases(27,:)={[0 0 0 -1;0 0 0 -1],[],'[curelt(:,[1 2 3]) curprop;curelt(:,[1 3 4]) curprop]',[]};


RE=struct('cases',vertcat(allcases{:,1}),  'cindex',reshape([1;1]*(1:size(allcases,1)),[],1)); %xxEM
% 'cindex',[reshape([1;1]*(1:23),[],1);24*ones(4,1);25*ones(4,1)]);
RE.newelt=[{'tria','quad'};allcases(:,3:4)];
RE.CutEdges=allcases(:,2);RE.ElemP='quad4';
if nargout==0; genericDisplay(RE);end % feval(lsutil('@reQuad'))

%% #reTetra : generic cut of tetra elements -3
function RE=reTetra
allcases=cell(25,5);
%one cut (6)
allcases(1,:)={[1 -1 0 0;-1 1 0 0],[1 2],'[curelt(:,[1 5 3 4]) curprop;curelt(:,[5 2 3 4]) curprop]',[],[]};
allcases(2,:)={[1 0 -1 0;-1 0 1 0],[1 3],'[curelt(:,[1 2 5 4]) curprop;curelt(:,[5 2 3 4]) curprop]',[],[]};
allcases(3,:)={[1 0 0 -1;-1 0 0 1],[1 4],'[curelt(:,[1 2 3 5]) curprop;curelt(:,[5 2 3 4]) curprop]',[],[]};
allcases(4,:)={[0 1 -1 0;0 -1 1 0],[2 3],'[curelt(:,[1 2 5 4]) curprop;curelt(:,[1 5 3 4]) curprop]',[],[]};
allcases(5,:)={[0 1 0 -1;0 -1 0 1],[2 4],'[curelt(:,[1 2 3 5]) curprop;curelt(:,[1 5 3 4]) curprop]',[],[]};
allcases(6,:)={[0 0 1 -1;0 0 -1 1],[3 4],'[curelt(:,[1 2 5 4]) curprop;curelt(:,[1 2 3 5]) curprop]',[],[]};
%two cuts (12)
allcases(7,:)={[1 -1 -1 0;-1 1 1 0],[1 2;1 3],'[curelt(:,[1 5 6 4]) curprop]','[curelt(:,[6 5 2 3 4]) curprop]',[]};
allcases(8,:)={[1 -1 0 -1;-1 1 0 1],[1 2;1 4],'[curelt(:,[1 5 3 6]) curprop]','[curelt(:,[4 2 5 6 3]) curprop]',[]};%
allcases(9,:)={[1 -1 0 1;-1 1 0 -1],[1 2;2 4],'[curelt(:,[5 2 3 6]) curprop]','[curelt(:,[4 6 5 1 3]) curprop]',[]};%
allcases(10,:)={[1 -1 1 0;-1 1 -1 0],[1 2;2 3],'[curelt(:,[5 2 6 4]) curprop]','[curelt(:,[1 5 6 3 4]) curprop]',[]};
allcases(11,:)={[1 0 -1 -1;-1 0 1 1],[1 3;1 4],'[curelt(:,[1 2 5 6]) curprop]','[curelt(:,[3 4 6 5 2]) curprop]',[]};
allcases(12,:)={[1 0 -1 1;-1 0 1 -1],[1 3;3 4],'[curelt(:,[5 2 3 6]) curprop]','[curelt(:,[1 5 6 4 2]) curprop]',[]};%
allcases(13,:)={[0 1 -1 -1;0 -1 1 1],[2 3;2 4],'[curelt(:,[1 2 5 6]) curprop]','[curelt(:,[4 3 5 6 1]) curprop]',[]};%
allcases(14,:)={[0 1 -1 1;0 -1 1 -1],[2 3;3 4],'[curelt(:,[1 5 3 6]) curprop]','[curelt(:,[4 6 5 2 1]) curprop]',[]};%
allcases(15,:)={[0 1 1 -1;0 -1 -1 1],[2 4;3 4],'[curelt(:,[1 5 6 4]) curprop]','[curelt(:,[3 2 5 6 1]) curprop]',[]};%
allcases(16,:)={[1 0 1 -1;-1 0 -1 1],[1 4;3 4],'[curelt(:,[5 2 6 4]) curprop]','[curelt(:,[6 5 1 3 2]) curprop]',[]};
allcases(17,:)={[1 1 -1 0;-1 -1 1 0],[1 3;2 3],'[curelt(:,[5 6 3 4]) curprop]','[curelt(:,[1 2 6 5 4]) curprop]',[]};
allcases(18,:)={[1 1 0 -1;-1 -1 0 1],[1 4;2 4],'[curelt(:,[5 6 3 4]) curprop]','[curelt(:,[5 6 2 1 3]) curprop]',[]};%
%three cuts (4)
allcases(19,:)={[1 -1 -1 -1;-1 1 1 1],[1 2;1 3;1 4],'[curelt(:,[1 5 6 7]) curprop]',[],'[curelt(:,[5 6 7 2 3 4]) curprop]'};
allcases(20,:)={[1 -1 1 1;-1 1 -1 -1],[1 2;2 3;2 4],'[curelt(:,[5 2 6 7]) curprop]',[],'[curelt(:,[1 3 4 5 6 7]) curprop]'};
allcases(21,:)={[1 1 -1 1;-1 -1 1 -1],[1 3;2 3;3 4],'[curelt(:,[5 6 3 7]) curprop]',[],'[curelt(:,[5 6 7 1 2 4]) curprop]'};
allcases(22,:)={[1 1 1 -1;-1 -1 -1 1],[1 4;2 4;3 4],'[curelt(:,[5 6 7 4]) curprop]',[],'[curelt(:,[1 2 3 5 6 7]) curprop]'};
%four cuts (3)
allcases(23,:)={[1 -1 -1 1;-1 1 1 -1],[1 2;1 3;2 4;3 4],[],[],'[curelt(:,[6 3 8 5 2 7]) curprop;curelt(:,[1 5 6 4 7 8]) curprop]'};
allcases(24,:)={[1 -1 1 -1;-1 1 -1 1],[1 2;1 4;2 3;3 4],[],[],'[curelt(:,[3 7 8 1 5 6]) curprop;curelt(:,[5 2 7 6 4 8]) curprop]'};
allcases(25,:)={[1 1 -1 -1;-1 -1 1 1],[1 3;1 4;2 3;2 4],[],[],'[curelt(:,[5 7 3 6 8 4]) curprop;curelt(:,[1 5 6 2 7 8]) curprop]'};

RE=struct('cases',vertcat(allcases{:,1}), ...
 'cindex',reshape([1;1]*(1:size(allcases,1)),[],1));
RE.newelt=[{'tetra','pyra','penta'};allcases(:,3:5)];
RE.CutEdges=allcases(:,2);  RE.ElemP='tetra4';
if nargout==0; genericDisplay(RE);end

%% #rePyra : generic cut of pyra elements -3
function RE=rePyra

allcases=cell(90,6);
%zero cut: LS on opposite nodes of base
allcases(1,:)={['-0+00';'+0-00'],[], ...
 [],'[curelt(:,[1 2 4 5]) curprop;curelt(:,[2 3 4 5]) curprop] ',[],[]};
allcases(2,:)={['0-0+0';'0+0-0'],[], ...
 [],'[curelt(:,[1 2 3 5]) curprop;curelt(:,[1 3 4 5]) curprop]',[],[]};
%two cuts: LS on peak + cut base
allcases(3,:)={['--++0';'++--0'],[1 4;2 3], ...
 [],[],'[curelt(:,[1 2 7 6 5]) curprop;curelt(:,[6 7 3 4 5]) curprop]',[]};
allcases(4,:)={['-++-0';'+--+0'],[1 2;3 4], ...
 [],[],'[curelt(:,[1 6 7 4 5]) curprop;curelt(:,[6 2 3 7 5]) curprop]',[]};
%four cuts
allcases(5,:)={['++++-';'----+'],[5 1; 5 2; 5 3; 5 4], ...
 '[curelt(:,[1:4 6:9]) curprop]',[],'[curelt(:,[6:9 5]) curprop]',[]};
%two cuts: 2 nodes on LS + cut of opposite face
allcases(6,:)={['00--+';'00++-'],[5 3; 5 4], [], [] ...
 '[curelt(:,[1 2 6 7 5]) curprop]','[curelt(:,[1 4 7 2 3 6]) curprop]'};
allcases(7,:)={['-00-+';'+00+-'],[1 5; 4 5], [], [] ...
 '[curelt(:,[2 3 7 6 5]) curprop]','[curelt(:,[1 6 2 4 7 3]) curprop]'};
allcases(8,:)={['--00+';'++00-'],[5 1; 5 2], [], [] ...
 '[curelt(:,[6 7 3 4 5]) curprop]','[curelt(:,[2 7 3 1 6 4]) curprop]'};
allcases(9,:)={['0--0+';'0++0-'],[2 5; 3 5], [], [] ...
 '[curelt(:,[1 6 7 4 5]) curprop]','[curelt(:,[1 6 2 4 7 3]) curprop]'};
%1 cut: 2 nodes on LS +  cut of opposite face
allcases(10,:)={['0-++0';'0+--0'],[2 3], [], ...
 '[curelt(:,[1 2 6 5]) curprop]', '[curelt(:,[1 6 3 4 5]) curprop]',[]};
allcases(11,:)={['0++-0';'0--+0'],[4 3], [], ...
 '[curelt(:,[1 6 4 5]) curprop]', '[curelt(:,[ 1 2 3 6 5]) curprop]',[]};
allcases(12,:)={['+0-+0';'-0+-0'],[4 3], [], ...
 '[curelt(:,[2 3 6 5]) curprop]', '[curelt(:,[ 1 2 6 4 5]) curprop]',[]};
allcases(13,:)={['-0++0';'+0--0'],[4 1], [], ...
 '[curelt(:,[1 2 6 5]) curprop]', '[curelt(:,[6 2 3 4 5]) curprop]',[]};
allcases(14,:)={['+-0+0';'-+0-0'],[1 2], [], ...
 '[curelt(:,[6 2 3 5]) curprop]', '[curelt(:,[1 6 3 4 5]) curprop]',[]};
allcases(15,:)={['--0+0';'++0-0'],[1 4], [], ...
 '[curelt(:,[6 3 4 5]) curprop]', '[curelt(:,[1 2 3 6 5]) curprop]',[]};
allcases(16,:)={['-++00';'+--00'],[1 2], [], ...
 '[curelt(:,[1 6 4 5]) curprop]', '[curelt(:,[6 2 3 4 5]) curprop]',[]};
allcases(17,:)={['--+00';'++-00'],[2 3], [], ...
 '[curelt(:,[6 3 4 5]) curprop]', '[curelt(:,[1 2 6 4 5]) curprop]',[]};
%three cuts (corners)
allcases(18,:)={['-++++';'+----'],[1 2; 1 4; 1 5], ...
 [],['[curelt(:,[1 6 7 8]) curprop;curelt(:,[8 3 5 6]) curprop; curelt(:,[7 3 5 8]) curprop;'...
 'curelt(:,[4 3 5 7]) curprop;curelt(:,[7 3 8 6]) curprop;curelt(:,[6 2 3 5]) curprop]'],[],[]};
allcases(19,:)={['+-+++';'-+---'],[2 1; 2 3; 2 5], ...
 [],['[curelt(:,[2 7 6 8]) curprop;curelt(:,[7 3 4 5]) curprop; curelt(:,[8 7 4 5]) curprop;' ...
 'curelt(:,[1 6 4 5]) curprop;curelt(:,[6 7 4 8]) curprop; curelt(:,[6 8 4 5]) curprop]'],[],[]};
allcases(20,:)={['++-++';'--+--'],[3 2; 3 4; 3 5], ...
 [],['[curelt(:,[3 7 6 8]) curprop;curelt(:,[2 5 6 1]) curprop;curelt(:,[5 8 6 1]) curprop;'...
 'curelt(:,[8 7 6 1]) curprop;curelt(:,[5 7 8 1]) curprop;curelt(:,[4 7 5 1]) curprop]'],[],[]};
allcases(21,:)={['+++-+';'---+-'],[4 1; 4 3; 4 5], ...
 [],['[curelt(:,[4 6:8]) curprop;curelt(:,[5 8 7 2]) curprop;curelt(:,[7 3 5 2]) curprop;'...
 'curelt(:,[6 7 8 2]) curprop;curelt(:,[8 5 6 2]) curprop;curelt(:,[6 5 1 2]) curprop]'],[],[]};
%1 cut: LS on 3 nodes of base
allcases(22,:)={['-000+';'+000-'],[1 5], [], ...
 '[curelt(:,[2 3 6 5]) curprop;curelt(:,[3 4 6 5]) curprop]','[curelt(:,[1:4 6]) curprop]',[]};
allcases(23,:)={['0-00+';'0+00-'],[2 5], [], ...
 '[curelt(:,[3 4 6 5]) curprop;curelt(:,[1 6 4 5]) curprop]', '[curelt(:,[1:4 6]) curprop]', []};
allcases(24,:)={['00-0+';'00+0-'],[3 5], [], ...
 '[curelt(:,[1 6 4 5]) curprop;curelt(:,[1 2 6 5]) curprop]','[curelt(:,[1:4 6]) curprop]',[]};
allcases(25,:)={['000-+';'000+-'],[4 5], [], ...
 '[curelt(:,[1 2 6 5]) curprop;curelt(:,[2 3 6 5]) curprop]','[curelt(:,[1:4 6]) curprop]',[]};
%1cut: LS on 2 nodes of base
allcases(26,:)={['-0+0+';'+0-0-'],[1 5], [], ...
 '[curelt(:,[1 2 4 6]) curprop;curelt(:,[2 3 4 5]) curprop;curelt(:,[2 4 6 5]) curprop]',[],[]};
allcases(27,:)={['0-0++';'0+0--'],[2 5], [], ...
 '[curelt(:,[1 2 3 6]) curprop;curelt(:,[1 3 4 5]) curprop;curelt(:,[1 6 3 5]) curprop]',[],[]};
allcases(28,:)={['+0-0+';'-0+0-'],[3 5], [], ...
 '[curelt(:,[2 3 4 6]) curprop;curelt(:,[4 2 6 5]) curprop;curelt(:,[1 2 4 5]) curprop]',[],[]};
allcases(29,:)={['0+0-+';'0-0+-'],[4 5], [], ...
 '[curelt(:,[1 3 4 6]) curprop;curelt(:,[1 3 6 5]) curprop;curelt(:,[1 2 3 5]) curprop]',[],[]};
%2cut: LS on node 1 + 2 cuts
allcases(30,:)={['0+---';'0-+++'],[2 3;2 5], [], ...
 '[curelt(:,[1 2 6 7]) curprop;curelt(:,[1 3 4 5]) curprop]', '[curelt(:,[3 6 7 5 1]) curprop]', []};
allcases(31,:)={['0--+-';'0++-+'],[3 4;4 5], [], ...
 '[curelt(:,[1 6 4 7]) curprop;curelt(:,[1 2 3 5]) curprop]', '[curelt(:,[3 5 7 6 1]) curprop]', []};
%2cut: LS on node 2 + 2 cuts
allcases(32,:)={['-0+++';'+0---'],[1 4;1 5], [], ...
 '[curelt(:,[1 2 6 7]) curprop;curelt(:,[2 3 4 5]) curprop]', '[curelt(:,[4 5 7 6 2]) curprop]', []};   %penta
allcases(33,:)={['-0+--';'+0-++'],[3 4;3 5], [], ...
 '[curelt(:,[2 3 6 7]) curprop;curelt(:,[1 2 4 5]) curprop]', '[curelt(:,[6 7 5 4 2]) curprop]', []};
%2cut: LS on node 3 + 2 cuts
allcases(34,:)={['++0-+';'--0+-'],[1 4;4 5], [], ...
 '[curelt(:,[3 4 6 7]) curprop;curelt(:,[1 2 3 5]) curprop]', '[curelt(:,[6 7 5 1 3]) curprop]', []};
allcases(35,:)={['+-0++';'-+0--'],[1 2;2 5], [], ...
 '[curelt(:,[2 3 6 7]) curprop;curelt(:,[1 3 4 5]) curprop]', '[curelt(:,[1 5 7 6 3]) curprop]', []};
%2cut: LS on node 4 + 2 cuts
allcases(36,:)={['+--0-';'-++0+'],[1 2;1 5], [], ...
 '[curelt(:,[1 6 4 7]) curprop;curelt(:,[3 4 2 5]) curprop]','[curelt(:,[2 6 7 5 4]) curprop]',[]};
allcases(37,:)={['++-0+';'--+0-'],[2 3;3 5], [], ...
 '[curelt(:,[6 3 4 7]) curprop;curelt(:,[1 2 4 5]) curprop]', '[curelt(:,[2 5 7 6 4]) curprop]',[]};
%2cut: LS on node 5 + 2 cuts
allcases(38,:)={['-+++0';'+---0'],[1 2;1 4],[], ...
 '[curelt(:,[1 6 7 5]) curprop;curelt(:,[2 3 6 5]) curprop;curelt(:,[6 3 7 5]) curprop;curelt(:,[7 3 4 5]) curprop]',[],[]};   %penta
allcases(39,:)={['+-++0';'-+--0'],[1 2;2 3],[], ...
 '[curelt(:,[1 6 4 5]) curprop;curelt(:,[2 7 6 5]) curprop;curelt(:,[7 3 4 5]) curprop;curelt(:,[6 7 4 5]) curprop]',[],[]};   %penta
allcases(40,:)={['++-+0';'--+-0'],[2 3;3 4], [], ...
 '[curelt(:,[1 7 4 5]) curprop;curelt(:,[1 6 7 5]) curprop;curelt(:,[1 2 6 5]) curprop;curelt(:,[6 3 7 5]) curprop]',[],[]};
allcases(41,:)={['+++-0';'---+0'],[3 4;1 4],[], ...
 '[curelt(:,[1 2 7 5]) curprop;curelt(:,[2 6 7 5]) curprop;curelt(:,[2 3 6 5]) curprop;curelt(:,[7 6 4 5]) curprop]',[],[]};   %penta
%2cut: LS on 2 nodes (5+i) + 2 cuts
allcases(42,:)={['+-0-0';'-+0+0'],[1 2;1 4],[], ...
 '[curelt(:,[1 6 7 5]) curprop;curelt(:,[2 3 6 5]) curprop;curelt(:,[6 3 7 5]) curprop;curelt(:,[7 3 4 5]) curprop]',[],[]};   %penta
allcases(43,:)={['+-+00';'-+-00'],[1 2;2 3],[], ...
 '[curelt(:,[1 6 4 5]) curprop;curelt(:,[2 7 6 5]) curprop;curelt(:,[7 3 4 5]) curprop;curelt(:,[6 7 4 5]) curprop]',[],[]};   %penta
allcases(44,:)={['0+-+0';'0-+-0'],[2 3;3 4], [], ...
 '[curelt(:,[1 7 4 5]) curprop;curelt(:,[1 6 7 5]) curprop;curelt(:,[1 2 6 5]) curprop;curelt(:,[6 3 7 5]) curprop]',[],[]};
allcases(45,:)={['+0+-0';'-0-+0'],[3 4;1 4],[], ...
 '[curelt(:,[1 2 7 5]) curprop;curelt(:,[2 6 7 5]) curprop;curelt(:,[2 3 6 5]) curprop;curelt(:,[7 6 4 5]) curprop]',[],[]};   %penta
%3cut: LS on 3 nodes (5+i) + 1 cut
allcases(46,:)={['+-000';'-+000'],[1 2],[], ... %hexa
 '[curelt(:,[1 6 4 5]) curprop;curelt(:,[6 2 3 5]) curprop;curelt(:,[6 3 4 5]) curprop]',... %tetra
 [],[]};   %penta
allcases(47,:)={['0+-00';'0-+00'],[2 3],[], ... %hexa
 '[curelt(:,[1 2 6 5]) curprop;curelt(:,[1 6 4 5]) curprop;curelt(:,[6 3 4 5]) curprop]',... %tetra
 [],[]};   %penta
allcases(48,:)={['00+-0';'00-+0'],[3 4],[], ... %hexa
 '[curelt(:,[1 6 4 5]) curprop;curelt(:,[1 2 6 5]) curprop;curelt(:,[2 3 6 5]) curprop]',... %tetra
 [],[]};   %penta
allcases(49,:)={['+00-0';'-00+0'],[1 4],[], ... %hexa
 '[curelt(:,[1 2 6 5]) curprop;curelt(:,[2 3 6 5]) curprop;curelt(:,[3 4 6 5]) curprop]',... %tetra
 [],[]};   %penta
%6cut:
i1=49;
i1=i1+1;allcases(i1,:)={['+-+--';'-+-++'],[1 2;2 3;3 4;1 4;1 5;3 5], [], ... %hexa
 '[curelt(:,[8 4 9 5]) curprop;curelt(:,[6 2 7 5]) curprop;curelt(:,[1 6 9 10]) curprop;curelt(:,[6 5 9 10]) curprop;curelt(:,[7 3 8 11]) curprop;curelt(:,[7 8 5 11]) curprop]',... %tetra
 '[curelt(:,[6:9 5]) curprop]', []};   %penta
i1=i1+1;allcases(i1,:)={['-+-+-';'+-+-+'],[1 2;2 3;3 4;1 4;2 5;4 5], [], ... %hexa
 '[curelt(:,[1 6 9 5]) curprop;curelt(:,[7 3 8 5]) curprop;curelt(:,[6 2 7 10]) curprop;curelt(:,[8 4 9 11]) curprop;;curelt(:,[8 9 5 11]) curprop;curelt(:,[6 7 5 10]) curprop]',... %tetra
 '[curelt(:,[6:9 5]) curprop]', []};   %penta
%5cuts
i1=i1+1;allcases(i1,:)={['-+++-';'+---+'],[1 2;1 4;2 5;3 5;4 5], [], ... %hexa
 '[curelt(:,[1 6 7 5]) curprop;curelt(:,[2 3 6 8]) curprop;curelt(:,[3 4 7 10]) curprop;curelt(:,[3 10 8 9]) curprop;curelt(:,[8 9 10 5]) curprop]',... %tetra
 '[curelt(:,[6 7 10 8 3]) curprop;curelt(:,[6 8 10 7 5]) curprop]', []};   %penta
i1=i1+1;allcases(i1,:)={['+-++-';'-+--+'],[3 2;2 1;3 5;4 5;1 5], [], ... %hexa
 '[curelt(:,[2 6 7 5]) curprop;curelt(:,[3 4 6 8]) curprop;curelt(:,[4 1 7 10]) curprop;curelt(:,[4 10 8 9]) curprop;curelt(:,[8 9 10 5]) curprop]',... %tetra
 '[curelt(:,[6 7 10 8 4]) curprop;curelt(:,[6 8 10 7 5]) curprop]', []};   %penta
i1=i1+1;allcases(i1,:)={['++-+-';'--+-+'],[3 4;2 3;4 5;1 5;2 5], [], ... %hexa
 '[curelt(:,[3 6 7 5]) curprop;curelt(:,[4 1 6 8]) curprop;curelt(:,[1 2 7 10]) curprop;curelt(:,[1 10 8 9]) curprop;curelt(:,[8 9 10 5]) curprop]',... %tetra
 '[curelt(:,[6 7 10 8 1]) curprop;curelt(:,[6 8 10 7 5]) curprop]', []};   %penta
i1=i1+1;allcases(i1,:)={['+++--';'---++'],[1 4;3 4;1 5;2 5;3 5], [], ... %hexa
 '[curelt(:,[6 7 4 5]) curprop;curelt(:,[1 2 6 8]) curprop;curelt(:,[2 3 7 10]) curprop;curelt(:,[2 10 8 9]) curprop;curelt(:,[8 9 10 5]) curprop]',... %tetra
 '[curelt(:,[6 7 10 8 2]) curprop;curelt(:,[6 8 10 7 5]) curprop]',[]};   %penta
%4cuts
i1=i1+1;allcases(i1,:)={['+-+-0';'-+-+0'],[1 2;2 3;3 4;1 4], [], ... %hexa
 '[curelt(:,[1 6 9 5]) curprop;curelt(:,[6 2 7 5]) curprop;curelt(:,[7 3 8 5]) curprop;curelt(:,[8 4 9 5]) curprop]',... %tetra
 '[curelt(:,[6:9 5]) curprop]', []};   %penta
%
i1=i1+1;allcases(i1,:)={['0+-+-';'0-+-+'],[2 3;3 4;2 5;4 5], [], ... %hexa
 '[curelt(:,[1 2 6 8]) curprop;curelt(:,[1 7 4 9]) curprop;curelt(:,[1 5 8 9]) curprop]',... %tetra
 '[curelt(:,[6 8 9 7 1]) curprop]', '[curelt(:,[6 3 7 8 5 9]) curprop]'};   %penta
i1=i1+1;allcases(i1,:)={['+0+--';'-0-++'],[3 4;1 4;3 5;1 5], [], ... %hexa
 '[curelt(:,[1 2 7 9]) curprop;curelt(:,[2 3 6 8]) curprop;curelt(:,[2 5 8 9]) curprop]',... %tetra
 '[curelt(:,[6 8 9 7 2]) curprop]', '[curelt(:,[6 4 7 8 5 9]) curprop]'};   %penta
i1=i1+1;allcases(i1,:)={['-+0+-';'+-0-+'],[1 4;1 2;4 5;2 5], [], ... %hexa
 '[curelt(:,[2 3 7 9]) curprop;curelt(:,[3 4 6 8]) curprop;curelt(:,[3 5 8 9]) curprop]',... %tetra
 '[curelt(:,[6 8 9 7 3]) curprop]', '[curelt(:,[6 1 7 8 5 9]) curprop]'};   %penta
i1=i1+1;allcases(i1,:)={['+-+0-';'-+-0+'],[1 2;2 3;1 5;3 5], [], ... %hexa
 '[curelt(:,[3 4 7 9]) curprop;curelt(:,[4 1 6 8]) curprop;curelt(:,[4 5 8 9]) curprop]',... %tetra
 '[curelt(:,[6 8 9 7 4]) curprop]', '[curelt(:,[6 2 7 8 5 9]) curprop]'};   %penta
%3cuts
i1=i1+1;allcases(i1,:)={['0-+--';'0+-++'],[2 3;3 4;3 5], [], ... %hexa
 '[curelt(:,[1 6 7 5]) curprop;curelt(:,[1 7 4 5]) curprop;curelt(:,[1 2 6 5]) curprop;curelt(:,[6 3 7 8]) curprop;curelt(:,[6 7 5 8]) curprop]',... %tetra
 [], []};   %penta
i1=i1+1;allcases(i1,:)={['-0-+-';'+0+-+'],[3 4;1 4;4 5], [], ... %hexa
 '[curelt(:,[2 6 7 5]) curprop;curelt(:,[2 7 1 5]) curprop;curelt(:,[2 3 6 5]) curprop;curelt(:,[6 4 7 8]) curprop;curelt(:,[6 7 5 8]) curprop]',... %tetra
 [], []};   %penta
i1=i1+1;allcases(i1,:)={['+-0--';'-+0++'],[1 4;1 2;1 5], [], ... %hexa
 '[curelt(:,[3 6 7 5]) curprop;curelt(:,[3 7 2 5]) curprop;curelt(:,[3 4 6 5]) curprop;curelt(:,[6 1 7 8]) curprop;curelt(:,[6 7 5 8]) curprop]',... %tetra
 [], []};   %penta
i1=i1+1;allcases(i1,:)={['-+-0-';'+-+0+'],[1 2;2 3;2 5], [], ... %hexa
 '[curelt(:,[4 6 7 5]) curprop;curelt(:,[4 7 3 5]) curprop;curelt(:,[4 1 6 5]) curprop;curelt(:,[6 2 7 8]) curprop;curelt(:,[6 7 5 8]) curprop]',... %tetra
 [], []};   %penta
%2 cuts: LS on 1 edge + 2 cuts
i1=i1+1;allcases(i1,:)={['00-+-';'00+-+'],[3 4;4 5], [], ... %hexa
 '[curelt(:,[1 6 4 7]) curprop;curelt(:,[1 2 6 7]) curprop;curelt(:,[1 5 2 7]) curprop]',... %tetra
 '[curelt(:,[3 5 7 6 2]) curprop]', []};   %penta
i1=i1+1;allcases(i1,:)={['00+--';'00-++'],[3 4;3 5], [], ... %hexa
 '[curelt(:,[2 3 6 7]) curprop;curelt(:,[1 2 6 7]) curprop;curelt(:,[1 2 7 5]) curprop]',... %tetra
 '[curelt(:,[4 6 7 5 1]) curprop]', []};   %penta
i1=i1+1;allcases(i1,:)={['-00+-';'+00-+'],[1 4;4 5], [], ... %hexa
 '[curelt(:,[3 4 6 7]) curprop;curelt(:,[2 3 6 7]) curprop;curelt(:,[2 3 7 5]) curprop]',... %tetra
 '[curelt(:,[1 6 7 5 2]) curprop]', []};   %penta
i1=i1+1;allcases(i1,:)={['+00--';'-00++'],[1 4;1 5], [], ... %hexa
 '[curelt(:,[1 2 6 7]) curprop;curelt(:,[2 3 6 7]) curprop;curelt(:,[2 3 7 5]) curprop]',... %tetra
 '[curelt(:,[4 5 7 6 3]) curprop]', []};   %penta
i1=i1+1;allcases(i1,:)={['+-00-';'-+00+'],[1 2;1 5], [], ... %hexa
 '[curelt(:,[4 1 6 7]) curprop;curelt(:,[3 4 6 7]) curprop;curelt(:,[3 4 7 5]) curprop]',... %tetra
 '[curelt(:,[6 7 5 2 3]) curprop]', []};   %penta
i1=i1+1;allcases(i1,:)={['-+00-';'+-00+'],[1 2;2 5], [], ... %hexa
 '[curelt(:,[2 3 6 7]) curprop;curelt(:,[3 4 6 7]) curprop;curelt(:,[3 4 7 5]) curprop]',... %tetra
 '[curelt(:,[1 5 7 6 4]) curprop]', []};   %penta
i1=i1+1;allcases(i1,:)={['0+-0-';'0-+0+'],[2 3;2 5], [], ... %hexa
 '[curelt(:,[1 2 6 7]) curprop;curelt(:,[4 1 6 7]) curprop;curelt(:,[4 1 7 5]) curprop]',... %tetra
 '[curelt(:,[3 6 7 5 4]) curprop]', []};   %penta
i1=i1+1;allcases(i1,:)={['0-+0-';'0+-0+'],[2 3;3 5], [], ... %hexa
 '[curelt(:,[3 4 6 7]) curprop;curelt(:,[4 1 6 7]) curprop;curelt(:,[4 1 7 5]) curprop]',... %tetra
 '[curelt(:,[2 5 7 6 1]) curprop]', []};   %penta
%2 cuts (with possible incompatible faces)
i1=i1+1;allcases(i1,:)={['0+0+-';'0-0-+'],[2 5;4 5], [], ... %hexa
 '[curelt(:,[1 6 3 7]) curprop;curelt(:,[1 2 3 6]) curprop;curelt(:,[1 3 4 7]) curprop;curelt(:,[1 5 6 7]) curprop;curelt(:,[3 5 7 6]) curprop]',... %tetra
 [], []};   %penta
i1=i1+1;allcases(i1,:)={['+0+0-';'-0-0+'],[1 5;3 5], [], ... %hexa
 '[curelt(:,[2 4 6 7]) curprop;curelt(:,[1 2 4 6]) curprop;curelt(:,[2 3 4 7]) curprop;curelt(:,[2 7 6 5]) curprop;curelt(:,[4 6 7 5]) curprop]',... %tetra
 [], []};   %penta
%3 cuts (incompatible faces)
i1=i1+1;allcases(i1,:)={['+++0-';'---0+'],[1 5;2 5;3 5], [], ... %hexa
 '[curelt(:,[6:8 5]) curprop;curelt(:,[6 8 4 5]) curprop]',... %tetra
 '[curelt(:,[1 6 8 3 4]) curprop]', '[curelt(:,[1:3 6:8]) curprop]'};   %penta
i1=i1+1;allcases(i1,:)={['0+++-';'0---+'],[2 5;3 5;4 5], [], ... %hexa
 '[curelt(:,[6:8 5]) curprop;curelt(:,[6 8 1 5]) curprop]',... %tetra
 '[curelt(:,[2 6 8 4 1]) curprop]', '[curelt(:,[2:4 6:8]) curprop]'};   %penta
i1=i1+1;allcases(i1,:)={['+0++-';'-0--+'],[3 5;4 5;1 5], [], ... %hexa
 '[curelt(:,[6:8 5]) curprop;curelt(:,[6 8 2 5]) curprop]',... %tetra
 '[curelt(:,[3 6 8 1 2]) curprop]', '[curelt(:,[3:4 1 6:8]) curprop]'};   %penta
i1=i1+1;allcases(i1,:)={['++0+-';'--0-+'],[4 5;1 5;2 5], [], ... %hexa
 '[curelt(:,[6:8 5]) curprop;curelt(:,[6 8 3 5]) curprop]',... %tetra
 '[curelt(:,[4 6 8 2 3]) curprop]', '[curelt(:,[4 1 2 6:8]) curprop]'};   %penta
i1=i1+1;allcases(i1,:)={['++---';'--+++'],[1 5;1 4;2 5;2 3], [], [],... %tetra
 '[curelt(:,[3 4 7 9 5]) curprop;curelt(:,[6 8 9 7 5]) curprop]', '[curelt(:,[1 7 6 2 9 8]) curprop]'};   %penta
i1=i1+1;allcases(i1,:)={['++--+';'--++-'],[4 5;1 4;3 5;2 3], [], [],... %tetra
 '[curelt(:,[1 2 9 7 5]) curprop;curelt(:,[6 7 9 8 5]) curprop]', '[curelt(:,[3 9 8 4 7 6]) curprop]'};   %penta
i1=i1+1;allcases(i1,:)={['+--+-';'-++-+'],[1 5;1 2;4 5;3 4], [], [],... %tetra
 '[curelt(:,[2 3 9 7 5]) curprop;curelt(:,[6 7 9 8 5]) curprop]', '[curelt(:,[4 9 8 1 7 6]) curprop]'};   %penta
i1=i1+1;allcases(i1,:)={['+--++';'-++--'],[2 5;1 2;3 5;3 4], [],  [],... %tetra
 '[curelt(:,[1 7 9 4 5]) curprop;curelt(:,[6 8 9 7 5]) curprop]', '[curelt(:,[2 7 6 3 9 8]) curprop]'};   %penta
%3 cuts (incompatible faces)
i1=i1+1;allcases(i1,:)={['++-0-';'--+0+'],[2 3;1 5;2 5], [], ... %hexa
 '[curelt(:,[4 7 8 5]) curprop;curelt(:,[1 8 4 7]) curprop]',... %tetra
 '[curelt(:,[1 2 6 4 8]) curprop;curelt(:,[3 6 8 5 4]) curprop]', []};   %penta
i1=i1+1;allcases(i1,:)={['++0--';'--0++'],[1 4;1 5;2 5], [], ... %hexa
 '[curelt(:,[3 5 7 8]) curprop;curelt(:,[2 3 7 8]) curprop]',... %tetra
 '[curelt(:,[1 2 3 6 7]) curprop;curelt(:,[4 5 7 6 3]) curprop]', []}; %penta
i1=i1+1;allcases(i1,:)={['+0-+-';'-0+-+'],[3 4;1 5;4 5], [], ... %hexa
 '[curelt(:,[2 8 7 5]) curprop;curelt(:,[2 7 8 1]) curprop]',... %tetra
 '[curelt(:,[1 2 6 4 8]) curprop;curelt(:,[5 8 6 3 2]) curprop]', []};   %penta
i1=i1+1;allcases(i1,:)={['+-0+-';'-+0-+'],[1 2;1 5;4 5], [], ... %hexa
 '[curelt(:,[3 8 7 5]) curprop;curelt(:,[1 8 7 6]) curprop]',... %tetra
 '[curelt(:,[1 6 3 4 8]) curprop;curelt(:,[6 7 5 2 3]) curprop]', []};   %penta
i1=i1+1;allcases(i1,:)={['0++--';'0--++'],[3 4;2 5;3 5], [], ... %hexa
 '[curelt(:,[1 7 8 5]) curprop;curelt(:,[1 8 7 2]) curprop]',... %tetra
 '[curelt(:,[1 2 3 6 8]) curprop;curelt(:,[4 6 8 5 1]) curprop]', []};   %penta
i1=i1+1;allcases(i1,:)={['-++0-';'+--0+'],[1 2;2 5;3 5], [], ... %hexa
 '[curelt(:,[4 7 8 5]) curprop;curelt(:,[6 8 7 2]) curprop]',... %tetra
 '[curelt(:,[2 3 4 6 8]) curprop;curelt(:,[1 5 7 6 4]) curprop]', []};   %penta
i1=i1+1;allcases(i1,:)={['0-++-';'0+--+'],[2 3;3 5;4 5], [], ... %hexa
 '[curelt(:,[1 7 8 5]) curprop;curelt(:,[4 7 8 1]) curprop]',... %tetra
 '[curelt(:,[3 4 1 6 7]) curprop;curelt(:,[2 5 7 6 1]) curprop]', []};   %penta
i1=i1+1;allcases(i1,:)={['-0++-';'+0--+'],[1 4;3 5;4 5], [], ... %hexa
 '[curelt(:,[2 7 8 5]) curprop;curelt(:,[3 7 8 2]) curprop]',... %tetra
 '[curelt(:,[2 3 4 6 8]) curprop;curelt(:,[1 6 8 5 2]) curprop]', []};   %penta

RE=struct('cases',{cellstr(vertcat(allcases{:,1}))}, ...
 'cindex',reshape([1;1]*(1:size(allcases,1)),[],1));
RE.newelt=[{'hexa','tetra','pyra','penta'};allcases(:,3:6)];
RE.CutEdges=allcases(:,2); RE.ElemP='pyra5';
if nargout==0; genericDisplay(RE);end % feval(lsutil('@rePyra')) %findedge(RE);


%% #rePenta : generic cut of penta elements -3
function RE=rePenta 

allcases=cell(13,6);
%LS on 2 nodes + 2 cuts on opposite side
allcases(1,:)={['+-0+-0';'-+0-+0'],[1 2;4 5],[],[],[],...
 '[curelt(:,[1 7 3 4 8 6]) curprop;curelt(:,[7 2 3 8 5 6]) curprop]'};   %penta
allcases(2,:)={['0+-0+-';'0-+0-+'],[2 3;5 6],[],[],[],...
 '[curelt(:,[1 2 7 4 5 8]) curprop;curelt(:,[1 7 3 4 8 6]) curprop]'};   %penta
allcases(3,:)={['-0+-0+';'+0-+0-'],[1 3;4 6],[],[],[],...
 '[curelt(:,[1 2 7 4 5 8]) curprop;curelt(:,[2 3 7 5 6 8]) curprop]'};   %penta
%3 cuts
allcases(4,:)={['+++---';'---+++'],[1 4;2 5;3 6],[],[],[],...
 '[curelt(:,[1:3 7:9]) curprop;curelt(:,[7:9 4:6]) curprop]'};   %penta
%4 cuts
allcases(5,:)={['+--+--';'-++-++'],[1 2;1 3;4 5;4 6], ...
 '[curelt(:,[7 2 3 8 9 5 6 10]) curprop]', [],[],... %pyra
 '[curelt(:,[1 7 8 4 9 10]) curprop]'};   %penta
allcases(6,:)={['-+--+-';'+-++-+'],[1 2;2 3;4 5;5 6], ...
 '[curelt(:,[1 7 8 3 4 9 10 6]) curprop]', [],[],... %pyra
 '[curelt(:,[7 2 8 9 5 10]) curprop]'};   %penta
allcases(7,:)={['--+--+';'++-++-'],[1 3;2 3;4 6;5 6], ...
 '[curelt(:,[1 2 8 7 4 5 10 9]) curprop]', [],[],... %pyra
 '[curelt(:,[7 8 3 9 10 6]) curprop]'};   %penta
%3 cuts (alternatives are possble)
allcases(8,:)={['+-----';'-+++++'],[1 2;1 3;1 4], [], ... %hexa
 '[curelt(:,[1 7 8 9]) curprop]','[curelt(:,[2 3 8 7 9]) curprop]',... %pyra
 '[curelt(:,[9 2:6]) curprop]'};   %penta
allcases(9,:)={['-+----';'+-++++'],[1 2;2 3;2 5], [], ... %hexa
 '[curelt(:,[2 8 7 9]) curprop]','[curelt(:,[1 7 8 3 9]) curprop]',... %pyra
 '[curelt(:,[1 9 3:6]) curprop]'};   %penta
allcases(10,:)={['--+---';'++-+++'],[1 3;2 3;3 6], [], ... %hexa
 '[curelt(:,[3 7 8 9]) curprop]','[curelt(:,[1 2 8 7 9]) curprop]',... %pyra
 '[curelt(:,[1 2 9 4:6]) curprop]'};   %penta
allcases(11,:)={['---+--';'+++-++'],[4 5;4 6;1 4], [], ... %hexa
 '[curelt(:,[4 8 7 9]) curprop]','[curelt(:,[5 7 8 6 9]) curprop]',... %pyra
 '[curelt(:,[1:3 9 5:6]) curprop]'};   %penta
allcases(12,:)={['----+-';'++++-+'],[4 5;5 6;2 5], [], ... %hexa
 '[curelt(:,[5 7 8 9]) curprop]','[curelt(:,[4 6 8 7 9]) curprop]',... %pyra
 '[curelt(:,[1:3 4 9 6]) curprop]'};   %penta
allcases(13,:)={['-----+';'+++++-'],[4 6;5 6;3 6], [], ... %hexa
 '[curelt(:,[6 8 7 9]) curprop]','[curelt(:,[4 7 8 5 9]) curprop]',... %pyra
 '[curelt(:,[1:3 4:5 9]) curprop]'};   %penta

% allcases(1,:)={['--+---';'++-+++'],[1 3;2 3;3 6], ...
%  [], ... %hexa
%  '[curelt(:,[3 7 8 9]) curprop]',... %tetra
%  '[curelt(:,[1 2 8 7 9]) curprop]',... %pyra
%  '[curelt(:,[1 2 9 4:6]) curprop]'};   %penta

RE=struct('cases',{cellstr(vertcat(allcases{:,1}))}, ...
 'cindex',reshape([1;1]*(1:size(allcases,1)),[],1));
RE.newelt=[{'hexa','tetra','pyra','penta'};allcases(:,3:6)];
RE.CutEdges=allcases(:,2); RE.ElemP='penta6';
if nargout==0; genericDisplay(RE);end % feval(lsutil('@rePenta'))


%% #reHexa : generic cut of hexa elements -3
function RE=reHexa

allcases=cell(54,6);
%four cuts
allcases(1,:)={['++++----';'----++++'],[1 5;2 6;3 7;4 8], ...
 '[curelt(:,[1:4 9:12]) curprop;curelt(:,[9:12 5:8]) curprop]',[],[],[]};
allcases(2,:)={['+--++--+';'-++--++-'],[1 2;3 4;5 6;7 8], ...
 '[curelt(:,[1 9 10 4 5 11 12 8]) curprop;curelt(:,[9 2 3 10 11 6 7 12]) curprop]',[],[],[]};
allcases(3,:)={['++--++--';'--++--++'],[1 4;2 3;5 8;6 7], ...
 '[curelt(:,[1 2 10 9 5 6 12 11]) curprop;curelt(:,[9 10 3 4 11 12 7 8]) curprop]',[],[],[]};
%cut faces
allcases(4,:)={['0+0-0+0-';'0-0+0-0+'],[],[],[],[], ...
 ['[curelt(:,[1:3 5:7]) curprop;curelt(:,[1 3 4 5 7 8]) curprop]']};
allcases(5,:)={['+0-0+0-0';'-0+0-0+0'],[],[],[],[], ...
 ['[curelt(:,[1 2 4 5 6 8]) curprop;curelt(:,[2:4 6:8]) curprop]']};
allcases(6,:)={['00--++00';'00++--00'],[],[],[],[], ...
 ['[curelt(:,[1 8 5 2 7 6]) curprop;curelt(:,[1 4 8 2 3 7]) curprop]']};
allcases(7,:)={['--0000++';'++0000--'],[],[],[],[], ...
 ['[curelt(:,[1 4 5 2 3 6]) curprop;curelt(:,[5 4 8 6 3 7]) curprop]']};
allcases(8,:)={['-00-0++0';'+00+0--0'],[],[],[],[], ...
 ['[curelt(:,[1 5 2 4 8 3]) curprop;curelt(:,[2 5 6 3 8 7]) curprop]']};
allcases(9,:)={['0--0+00+';'0++0-00-'],[],[],[],[], ...
 ['[curelt(:,[1 6 2 4 7 3]) curprop;curelt(:,[1 5 6 4 8 7]) curprop]']};
% two cuts: edges on LS + cut of opposite face
allcases(10,:)={['00--++++';'00++----'],[3 7; 4 8],...
 '[curelt(:,[1 2 9 10 5:8]) curprop]',[],[], '[curelt(:,[1 4 10 2 3 9]) curprop]'};
allcases(11,:)={['00--++--';'00++--++'],[6 7; 5 8],...
 '[curelt(:,[1:4 10 9 7 8]) curprop]',[],[], '[curelt(:,[1 10 5 2 9 6]) curprop]'};
allcases(12,:)={['-00-++++';'+00+----'],[1 5; 4 8],...
 '[curelt(:,[9 2 3 10 5:8]) curprop]',[],[], '[curelt(:,[1 9 2 4 10 3]) curprop]'};
allcases(13,:)={['-00--++-';'+00++--+'],[5 6; 8 7],...
 '[curelt(:,[1:4 5 9 10 8]) curprop]',[],[], '[curelt(:,[2 9 6 3 10 7]) curprop]'};
allcases(14,:)={['--00++++';'++00----'],[1 5; 2 6],...
 '[curelt(:,[9 10 3 4 5:8]) curprop]',[],[], '[curelt(:,[1 4 9 2 3 10]) curprop]'};
allcases(15,:)={['--00--++';'++00++--'],[5 8; 6 7],...
 '[curelt(:,[1:4 5 6 10 9]) curprop]',[],[], '[curelt(:,[3 10 7 4 9 8]) curprop]'};
allcases(16,:)={['0--0++++';'0++0----'],[2 6; 3 7],...
 '[curelt(:,[1 9 10 4 5:8]) curprop]',[],[], '[curelt(:,[1 9 2 4 10 3]) curprop]'};
allcases(17,:)={['0--0+--+';'0++0-++-'],[5 6; 8 7],...
 '[curelt(:,[1:4 9 6 7 10]) curprop]',[],[], '[curelt(:,[1 5 9 4 8 10]) curprop]'};
allcases(18,:)={['0++-0++-';'0--+0--+'],[4 3 ; 8 7],...
 '[curelt(:,[1 2 3 9 5 6 7 10]) curprop]',[],[], '[curelt(:,[1 9 4 5 10 8]) curprop]'};
allcases(19,:)={['0-++0-++';'0+--0+--'],[2 3; 6 7],...
 '[curelt(:,[1 9 3 4 5 10 7 8]) curprop]',[],[], '[curelt(:,[1 2 9 5 6 10]) curprop]'};
allcases(20,:)={['-0++-0++';'+0--+0--'],[1 4; 5 8],...
 '[curelt(:,[9 2 3 4 10 6 7 8]) curprop]',[],[], '[curelt(:,[5 10 6 1 9 2 ]) curprop]'};
allcases(21,:)={['-0+--0+-';'+0-++0-+'],[4 3; 8 7],...
 '[curelt(:,[1 2 9 4 5 6 10 8]) curprop]',[],[], '[curelt(:,[2 3 9 6 7 10]) curprop]'};
allcases(22,:)={['++0-++0-';'--0+--0+'],[1 4 ; 5 8],...
 '[curelt(:,[1 2 3 9 5 6 7 10]) curprop]',[],[], '[curelt(:,[4 9 3 8 10 7]) curprop]'};
allcases(23,:)={['-+0--+0-';'+-0++-0+'],[1 2; 5 6],...
 '[curelt(:,[1 9 3 4 5 10 7 8]) curprop]',[],[], '[curelt(:,[9 2 3 10 6 7]) curprop]'};
allcases(24,:)={['-++0-++0';'+--0+--0'],[1 2 ; 5 6],...
 '[curelt(:,[9 2 3 4 10 6 7 8]) curprop]',[],[], '[curelt(:,[1 9 4 5 10 8]) curprop]'};
allcases(25,:)={['++-0++-0';'--+0--+0'],[2 3; 6 7],...
 '[curelt(:,[1 2 9 4 5 6 10 8]) curprop]',[],[], '[curelt(:,[4 9 3 8 10 7]) curprop]'};
allcases(26,:)={['--++00++';'++--00--'],[1 4 ; 2 3],...
 '[curelt(:,[9 10 3 4 5:8]) curprop]',[],[], '[curelt(:,[1 9 5 2 10 6]) curprop]'};
allcases(27,:)={['----00++';'++++00--'],[3 7 ; 4 8],...
 '[curelt(:,[1:4 5 6 9 10]) curprop]',[],[], '[curelt(:,[5 10 8 6 9 7]) curprop]'};
allcases(28,:)={['----+00+';'++++-00-'],[1 5 ; 4 8],...
 '[curelt(:,[1:4 9 6 7 10]) curprop]',[],[], '[curelt(:,[9 5 6 10 8 7 ]) curprop]'};
allcases(29,:)={['+--++00+';'-++--00-'],[1 2 ; 4 3],...
 '[curelt(:,[1 9 10 4 5:8]) curprop]',[],[], '[curelt(:,[9 6 2 10 7 3]) curprop]'};
allcases(30,:)={['++--++00';'--++--00'],[1 4 ; 2 3],...
 '[curelt(:,[1 2 10 9 5:8]) curprop]',[],[], '[curelt(:,[9 4 8 10 3 7 ]) curprop]'};
allcases(31,:)={['----++00';'++++--00'],[1 5 ; 2 6],...
 '[curelt(:,[1 :4 9 10 7 8]) curprop]',[],[], '[curelt(:,[9 8 5 10 7 6]) curprop]'};
allcases(32,:)={['----0++0';'++++0--0'],[2 6 ; 3 7],...
 '[curelt(:,[1:4 5 9 10 8]) curprop]',[],[], '[curelt(:,[5 6 9 8 7 10 ]) curprop]'};
allcases(33,:)={['-++-0++0';'+--+0--0'],[1 2 ; 4 3],...
 '[curelt(:,[9 2 3 10 5:8]) curprop]',[],[], '[curelt(:,[1 5 9 4 8 10]) curprop]'};
%4 cuts by a plane
allcases(34,:)={['-+++-+++';'+---+---'],[1 2 ; 1 4; 5 6; 5 8],[],[],[], ...
 '[curelt(:,[1 9 10 5 11 12]) curprop;curelt(:,[9 2 3 11 6 7]) curprop;curelt(:,[9 3 10 11 7 12]) curprop;curelt(:,[10 3 4 12 7 8]) curprop]'};
allcases(35,:)={['-+---+--';'+-+++-++'],[2 1 ; 2 3; 5 6; 6 7],[],[],[], ...
 '[curelt(:,[11 12 6 9 10 2]) curprop;curelt(:,[1 9 4 5 11 8]) curprop;curelt(:,[11 8 12 9 4 10]) curprop;curelt(:,[10 3 4 12 7 8]) curprop]'};
allcases(36,:)={['--+---+-';'++-+++-+'],[3 2 ; 3 4; 7 6; 7 8],[],[],[], ...
 '[curelt(:,[9 3 10 11 7 12]) curprop;curelt(:,[1 2 9 5 6 11]) curprop;curelt(:,[5 8 12 1 4 10]) curprop;curelt(:,[1 9 10 5 11 12]) curprop]'};
allcases(37,:)={['---+---+';'+++-+++-'],[4 1 ; 4 3; 8 5; 8 7],[],[],[], ...
 '[curelt(:,[1 2 9 5 6 11]) curprop;curelt(:,[9 2 10 11 6 12]) curprop;curelt(:,[10 2 3 12 6 7]) curprop;curelt(:,[11 8 12 9 4 10]) curprop]'};
allcases(38,:)={['--++++++';'++------'],[1 4 ; 1 5; 2 3; 2 6],[],[],[], ...
 '[curelt(:,[1 9 10 2 11 12]) curprop;curelt(:,[11 7 3 9 8 4]) curprop;curelt(:,[12 7 11 10 8 9]) curprop;curelt(:,[6 7 12 5 8 10]) curprop]'};
allcases(39,:)={['+--+++++';'-++-----'],[2 1 ; 2 6; 3 4; 3 7],[],[],[], ...
 '[curelt(:,[9 10 2 11 12 3]) curprop;curelt(:,[1 5 9 4 8 11]) curprop;curelt(:,[9 5 10 11 8 12]) curprop;curelt(:,[5 6 10 8 7 12]) curprop]'};
allcases(40,:)={['++--++++';'--++----'],[3 2; 3 7; 4 1; 4 8],[],[],[], ...
 '[curelt(:,[9 10 3 11 12 4]) curprop;curelt(:,[2 6 9 1 5 11]) curprop;curelt(:,[9 6 10 11 5 12]) curprop;curelt(:,[5 12 8 6 10 7]) curprop]'};
allcases(41,:)={['-++-++++';'+--+----'],[1 2; 1 5; 4 3; 4 8],[],[],[], ...
 '[curelt(:,[1 10 9 4 12 11]) curprop;curelt(:,[9 6 2 11 7 3]) curprop;curelt(:,[10 6 9 12 7 11]) curprop;curelt(:,[5 6 10 8 7 12]) curprop]'};
allcases(42,:)={['++++--++';'----++--'],[5 1; 5 8; 6 2; 6 7],[],[],[], ...
 '[curelt(:,[5 9 10 6 11 12]) curprop;curelt(:,[1 4 9 2 3 11]) curprop;curelt(:,[9 4 10 11 3 12]) curprop;curelt(:,[10 4 8 12 3 7]) curprop]'};
allcases(43,:)={['+++++--+';'-----++-'],[6 2; 6 5; 7 3; 7 8],[],[],[], ...
 '[curelt(:,[1 9 2 4 11 3]) curprop;curelt(:,[1 10 9 4 12 11 ]) curprop;curelt(:,[1 5 10 4 8 12 ]) curprop;curelt(:,[9 10 6 11 12 7]) curprop]'};
allcases(44,:)={['++++++--';'------++'],[7 6; 7 3; 8 5; 8 4],[],[],[], ...
 '[curelt(:,[9 7 10 11 8 12]) curprop;curelt(:,[2 10 3 1 12 4]) curprop;curelt(:,[2 9 10 1 11 12]) curprop;curelt(:,[2 6 9 1 5 11]) curprop]'};
allcases(45,:)={['++++-++-';'----+--+'],[5 1; 5 6; 8 4; 8 7],[],[],[], ...
 '[curelt(:,[5 10 9 8 12 11]) curprop;curelt(:,[1 9 2 4 11 3]) curprop;curelt(:,[9 10 2 11 12 3]) curprop;curelt(:,[2 10 6 3 12 7]) curprop]'};
%2 cuts
allcases(46,:)={['+++00--0';'---00++0'],[2 6;3 7],'[curelt(:,[1:4 5 9 10 8]) curprop]',[],[],'[curelt(:,[8 10 7 5 9 6]) curprop]'};




if 1==2
 %3 cuts by a plane (corners) %cut diag
 allcases(46,:)={['-+++++++';'+-------'],[1 2 ; 1 4; 1 5 ; 1 7],[],...
  ['[curelt(:,[1 9:11]) curprop;curelt(:,[9 12 6 2]) curprop;curelt(:,[4 8 10 12]) curprop;'...
  'curelt(:,[10 8 11 12]) curprop;curelt(:,[5 11 8 12]) curprop;curelt(:,[9 12 11 6]) curprop;'...
  'curelt(:,[10 12 11 9]) curprop;curelt(:,[4 12 10 3]) curprop;curelt(:,[10 12 9 3]) curprop;'...
  'curelt(:,[11 5 6 12]) curprop;curelt(:,[9 12 2 3]) curprop]'],...
  '[curelt(:,[2 6 7 3 12]) curprop;curelt(:,[3 7 8 4 12]) curprop;curelt(:,[8:-1:5 12]) curprop]',[]};
 %cut diag
 allcases(47,:)={['+-++++++';'-+------'],[2 1; 2 3; 2 6; 2 8],[],...
  ['[curelt(:,[2 10 9 11]) curprop;curelt(:,[12 5 1 9]) curprop;curelt(:,[9 1 12 4]) curprop;'...
  'curelt(:,[11 9 12 10]) curprop;curelt(:,[3 10 12 4]) curprop;curelt(:,[10 9 12 4]) curprop;'...
  'curelt(:,[7 6 12 11]) curprop;curelt(:,[12 5 11 6]) curprop;curelt(:,[12 5 9 11]) curprop;'...
  'curelt(:,[7 11 12 10]) curprop;curelt(:,[7 10 12 3]) curprop]'],...
  '[curelt(:,[3 7 8 4 12]) curprop;curelt(:,[1 4 8 5 12]) curprop;curelt(:,[8:-1:5 12]) curprop]',[]};
 %cut diag
 allcases(48,:)={['++-+++++';'--+-----'],[3 2; 3 4; 3 7; 3 5],[],...
  ['[curelt(:,[9 3 10 11]) curprop;curelt(:,[1 12 2 9]) curprop;curelt(:,[7 12 8 11]) curprop;'...
  'curelt(:,[10 12 8 4]) curprop;curelt(:,[1 12 10 4]) curprop;curelt(:,[11 12 8 10]) curprop;'...
  'curelt(:,[6 12 7 11]) curprop;curelt(:,[1 12 9 10]) curprop;curelt(:,[9 12 11 10]) curprop;'...
  'curelt(:,[6 12 11 9]) curprop;curelt(:,[2 6 9 12]) curprop]'],...
  '[curelt(:,[1 4 8 5 12]) curprop;curelt(:,[1 5 6 2 12]) curprop;curelt(:,[8:-1:5 12]) curprop]',[]};
 %cut diag
 allcases(49,:)={['+++-++++';'---+----'],[4 3 ; 4 1; 4 8; 4 6],[],...
  ['[curelt(:,[9 4 10 11]) curprop;curelt(:,[10 12 2 9]) curprop;curelt(:,[10 5 1 12]) curprop;'...
  'curelt(:,[11 12 10 9]) curprop;curelt(:,[9 7 12 3]) curprop;curelt(:,[11 8 12 7]) curprop;'...
  'curelt(:,[9 12 2 3]) curprop;curelt(:,[11 5 10 12]) curprop;curelt(:,[11 8 5 12]) curprop;'...
  'curelt(:,[11 7 12 9]) curprop;curelt(:,[10 12 1 2]) curprop]'],...
  '[curelt(:,[1 5 6 2 12]) curprop;curelt(:,[2 6 7 3 12]) curprop;curelt(:,[8:-1:5 12]) curprop]',[]};
 %cut diag
 allcases(50,:)={['++++-+++';'----+---'],[5 1 ; 5 6; 5 8; 5 3],[],...
  ['[curelt(:,[9:11 5]) curprop;curelt(:,[4 11 9 12]) curprop;curelt(:,[4 9 1 12]) curprop;'...
  'curelt(:,[11 10 9 12]) curprop;curelt(:,[4 8 11 12]) curprop;curelt(:,[11 8 7 12]) curprop;'...
  'curelt(:,[2 10 6 12]) curprop;curelt(:,[10 7 6 12]) curprop;curelt(:,[11 7 10 12]) curprop;'...
  'curelt(:,[1 9 2 12]) curprop;curelt(:,[2 9 10 12]) curprop]'],...
  '[curelt(:,[1:4 12]) curprop;curelt(:,[2 6 7 3 12]) curprop;curelt(:,[4 3 7 8 12]) curprop]',[]};
 %cut diag
 allcases(51,:)={['+++++-++';'-----+--'],[6 5 ; 6 7; 6 2; 6 4],[],...
  ['[curelt(:,[9 11 10 6]) curprop;curelt(:,[3 12 10 7]) curprop;curelt(:,[7 12 10 8]) curprop;'...
  'curelt(:,[1 5 9 12]) curprop;curelt(:,[12 5 9 8]) curprop;curelt(:,[10 12 9 8]) curprop;'...
  'curelt(:,[3 2 11 12]) curprop;curelt(:,[2 1 11 12]) curprop;curelt(:,[10 11 9 12]) curprop;'...
  'curelt(:,[11 1 9 12]) curprop;curelt(:,[3 11 10 12]) curprop]'],...
  '[curelt(:,[1:4 12]) curprop;curelt(:,[3 7 8 4 12]) curprop;curelt(:,[5 1 4 8 12]) curprop]',[]};
 %cut diag
 allcases(52,:)={['++++++-+';'------+-'],[7 6 ; 7 8; 7 3; 7 1],[],...
  ['[curelt(:,[9 11 10 7]) curprop;curelt(:,[3 12 11 4]) curprop;curelt(:,[6 12 9 2]) curprop;'...
  'curelt(:,[12 5 10 8]) curprop;curelt(:,[12 8 10 4]) curprop;curelt(:,[9 12 11 2]) curprop;'...
  'curelt(:,[11 12 10 4]) curprop;curelt(:,[2 12 11 3]) curprop;curelt(:,[9 5 12 6]) curprop;'...
  'curelt(:,[10 5 12 9]) curprop;curelt(:,[10 9 12 11]) curprop]'],...
  '[curelt(:,[1:4 12]) curprop;curelt(:,[1 4 8 5 12]) curprop;curelt(:,[1 5 6 2 12]) curprop]',[]};
 %cut diag
 allcases(53,:)={['+++++++-';'-------+'],[8 5 ; 8 7; 8 4; 8 2],[],...
  ['[curelt(:,[9:11 8]) curprop;curelt(:,[7 6 10 12]) curprop;curelt(:,[11 1 4 12]) curprop;'...
  'curelt(:,[11 12 4 3]) curprop;curelt(:,[7 12 10 3]) curprop;curelt(:,[11 10 12 3]) curprop;'...
  'curelt(:,[6 12 9 10]) curprop;curelt(:,[11 9 12 10]) curprop;curelt(:,[11 9 1 12]) curprop;'...
  'curelt(:,[5 1 9 12]) curprop;curelt(:,[5 12 9 6]) curprop]'],...
  '[curelt(:,[1:4 12]) curprop;curelt(:,[2 6 7 3 12]) curprop;curelt(:,[1 5 6 2 12]) curprop]',[]};
 %cut diag do other
 % 'xxx'
 % allcases(54,:)={['++00+0-0';'--00-0+0'],[1 7],[],...
 %    '[curelt(:,[3 7 8 6]) curprop;curelt(:,[2 3 9 6]) curprop;curelt(:,[3 8 4 9]) curprop;curelt(:,[5 8 6 9]) curprop]',...
 %     '[curelt(:,[1:4 9]) curprop;curelt(:,[1 5 6 2 9]) curprop;curelt(:,[1 4 8 5 9]) curprop]',[]};
end


% allcases(1,:)={['+++00--0';'---00++0'],[2 6;3 7], ...
%  '[curelt(:,[1:4 5 9 10 8]) curprop]',... %hexa
%  [],... %tetra
%  [],... %pyra
%  '[curelt(:,[8 10 7 5 9 6]) curprop]'}; %penta





% one hexa and one extruded edge xxx
% allcases(10,:)={['++++-000';'----+000'],[1 5], '[curelt(:,[1:4 9 6 7 8]) curprop]',[],'[curelt(:,[5:8 9]) curprop]',[]}; %em

RE=struct('cases',{cellstr(vertcat(allcases{:,1}))}, ...
 'cindex',reshape([1;1]*(1:size(allcases,1)),[],1));
RE.newelt=[{'hexa','tetra','pyra','penta'};allcases(:,3:6)];
RE.CutEdges=allcases(:,2); RE.ElemP='hexa8';
if nargout==0; genericDisplay(RE);end % feval(lsutil('@reHexa'))


%% #genericNameMap : name based cutting
% feval(lsutil('@genericNameMap'));
function out=genericNameMap(varargin);

persistent cutMap
if isempty(cutMap)
 %% standard orient cuts
 cutMap=containers.Map; 
 RE=reQuad; 
 for j1=1:size(RE.CutEdges,1)
   r2=RE.CutEdges{j1};
   n2=cell(4+size(r2,1),1);
   r1=struct;
   for j2=1:4; n2{j2}=@(x)sprintf('n%i',x(j2));end 
   for j2=1:size(r2,1)
    n2{j2+4}=@(x)sprintf('e%i.%i',sort(x(r2(j2,:))));
   end
   r1.Node=n2;r1.elt={};
   r3=RE.newelt(j1+1,:);
   for j2=1:length(r3); 
     if ~isempty(r3{j2})
       curelt=1:length(r1.Node); curprop=[];
       r1.elt(end+[1:2])={RE.newelt{1,j2}; eval(r3{j2})};
     end
   end
   st=charLs(RE.cases(RE.cindex==j1,:));
   for j2=1:length(st); cutMap(st{j2})=r1;end
 end
 %% refine cuts : needs implementation xxxEB 
 doAdd=@(x)feval(lsutil('@genericNameMap'),'combinations',x);
 r1=struct('name','+0--c','node',{{1,2,3,4,[3 4],[1:4]}}, ...
      'elt',{{'tria',[2 6 5],'quad4',[1 5 6 4;2 3 4 6]}});
 doAdd(r1);
 r1=struct('name','0--+c','node',{{1,2,3,4,[1 2],[1:4]}}, ...
      'elt',{{'tria',[1 5 6],'quad4',[5 2 3 6;1 6 3 4]}});
 doAdd(r1);
 
 if 1==2
  cutMap('q1121')=struct('node',{{1,2,3,4,[3 4]}},'elt',{{'tria',[1 2 5;2 3 5;1 5 4]}});
  
  % NodeName e(n1).(n2) with (n1<2) name of midside node
  % e(n1).(n2).(ind=2/3)
  % f(n1).(n2).(n3) for face nodes
  % v(n1).(n2).(n3) for face nodes
  % sdtweb feutil#refineCell
  % sdtweb lsutil face2tria
  % quad4('orders')
 end
end
[CAM,Cam]=comstr(varargin{1},1);
if comstr(Cam,'docut')
 %% #docut
 model=varargin{2};carg=3; RO=varargin{carg};carg=carg+1;
 [EGroup,nGroup]=getegroup(model.Elt);
 for jGroup=1:nGroup
  [ElemF,i1,ElemP]=getegroup(model.Elt(EGroup(jGroup),:),jGroup);
  i2=feval(ElemF,'prop');i2(3:end)=[];
  cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
  elt=full(RO.NNode(model.Elt(cEGI,feval(ElemP,'node'))));
  if isfield(RO,'withNodeInd')
   i1=~any(ismember(elt,RO.withNodeInd),2);cEGI(i1)=[];elt(i1,:)=[];
  end
  if isfield(RO,'inEltInd')
   i1=~ismember(cEGI,RO.inEltInd);cEGI(i1)=[];elt(i1,:)=[];
  end
      
  st=charLs(RO,elt);
  for j1=1:size(elt);
   st2=sprintf('f%i.%i.%i.%i',sort(elt(j1,:)));
   % fecom('shownodemark',{model.Node(elt(j1,:),1);RO.cuts})
   if RO.nmap.isKey(st2);
     st{j1}(end+1)='c';st2=st{j1};
     if ~any(st2=='+')||~any(st2=='-'); % an edge is cut and corner added defg
      i3=find(cellfun(@(x)RO.nmap.isKey(x),edgeName(elt(j1,:),feval(ElemP,'Edge'))));
      st{j1}=[st{j1}(1:4) sprintf('c.e%i',i3)];
      %st2(quad4('edges'))
     end
     if ~cutMap.isKey(st{j1});error('Unknown cut %s',st{j1});end
   end
  end
  i1=cellfun(@(x)~cutMap.isKey(x),st);cEGI(i1)=[];st(i1)=[];elt(i1,:)=[];% no map
  r3=struct('tria',{cell(length(cEGI),1)},'quad',{cell(length(cEGI),1)});
  for jElt=1:length(cEGI)
   r1=cutMap(st{jElt});
   st1=cellfun(@(x)x(elt(jElt,:)),r1.Node,'uni',0);
   %cellfun(@(x)RO.nmap.isKey(x),st1)
   curelt=model.Node(cellfun(@(x)RO.nmap(x),st1),1)';
   prop=model.Elt(cEGI(jElt),i2);
   for j2=1:2:length(r1.elt)
     i3=curelt(r1.elt{j2+1});
     r3.(r1.elt{j2}){jElt}=[i3 repmat(prop,size(i3,1),1)];
   end
  end
  model.Elt(cEGI,:)=[]; 
  i1=vertcat(r3.tria{:}); 
  if ~isempty(i1);model=feutil('addelt',model,'tria3',i1);end
  i1=vertcat(r3.quad{:}); 
  if ~isempty(i1);model=feutil('addelt',model,'quad4',i1);end
 end
 if isfield(RO,'ProId')
  mpid=feutil('mpid',model);
  mpid(feutil('findelt withnode',model,model.Node(RO.def>0)),2)=RO.ProId;
  model.Elt=feutil('mpid',model,mpid);%feplot(model);fecom showfipro
 end
 out=model; 
 
elseif cutMap.isKey(CAM);
    disp(cutMap(CAM));
elseif comstr(Cam,'combinations')
 %% #combinations
 r1=varargin{2};
 st1='';st1(abs('+-0c'))='-+0c';
 st=strrep(r1.name,'c','');
 if length(st)==4 % quad xxx need to use q+0-- for quad/tetra
  i1=quad4('orders');
  cutMap(r1.name)=r1; cutMap(st1(abs(r1.name)))=r1;
  for j1=2:size(i1,1)
   r2=r1; nind=sparse(i1(j1,:),1,i1(1,:));
   r2.name(1:4)=r2.name(i1(j1,:));
   r2.node=cellfun(@(x)full(nind(x)),r2.node,'uni',0);
   cutMap(r1.name)=r2; cutMap(st1(abs(r2.name)))=r2;
  end
 else;
  dbstack; keyboard; 
 end
 
else; error('%s unknown / not known key',CAM);
end

%#edgeName(elt(j1,:),ElemP(Edge))
function out=edgeName(elt,i1);

out=cell(1,size(i1,1));
for j1=1:size(i1,1); out{j1}=sprintf('e%i.%i',sort(elt(i1(j1,:))));end

%% #genericDisplay : see all cases
function genericDisplay(RE);

i1=unique(RE.cindex);
EC=integrules(RE.ElemP);
if ~iscell(RE.cases);RE.cases=charLs(RE.cases);end
map=containers.Map({'tetra','penta','tria','hexa','pyra','quad'}, ...
 {'tetra4','penta6','tria3','hexa8','pyra5','quadb'});
cf=sdth.urn('feplot(2)');
cingui('plotwd',cf,'@OsDic(SDT Root)',{'FnTM','ImToFigN','ImSw80','WrW33c'});

for j1=1:length(i1)
 mo1=struct('Node',[(1:size(EC.xi,1))'*[1 0 0 0] EC.xi],'Elt',[]);
 j2=find(RE.cindex==i1(j1),1,'first');
 i2=RE.CutEdges{i1(j1)};
 if ~isempty(i2)
  n2=(EC.xi(i2(:,1),:)+EC.xi(i2(:,2),:))/2;
  mo1.Node(end+(1:size(n2,1)),[1 5:7])=[(1:size(n2,1))'+size(mo1.Node,1) n2];
 end
 curprop=[1 1]; curelt=mo1.Node(:,1)';
 
 for j3=1:size(RE.newelt,2) % Add new elements
  if isempty(RE.newelt{i1(j1)+1,j3});continue;end
  mo1=feutil('addelt',mo1,map(RE.newelt{1,j3}),eval(RE.newelt{i1(j1)+1,j3}));
 end
 mo1.name=RE.cases{j2};
 mo1.pl=m_elastic('dbval 1 steel');mo1.il=p_solid('dbval 1 d3 -3');
 cf.model=mo1;fecom('colorfacew-edgealpha.2 -alpha.6');cf.os_('LgMdlName')
 if max(mo1.Node(:,7));fecom('view3');end;feplot;fecom('textnode');
 fprintf('Orient for view %s\n',mo1.name); keyboard;
 %sdth.urn('feplot(2).Report');
 %cf=feplot(mo1);cf.model=mo1;fe_quality(cf.mdl);fe_quality('view',cf);keyboard
end


%% #genericCut : implements actual building of new elements
function i1=genericCut(RE,elt,ElemP,sevLS,idedges,RO)
%RE=iso_zero(ElemP);
[RC,RE]=combiLS(ElemP);
RE.iso0=[]; % default not needed
if iscell(RE.cases)
 sevLS=charLs(sevLS);
 [in1,RE.cloc]=ismember(sevLS,RE.cases);
 i1=RE.cloc~=0; RE.cloc(i1)=RE.cindex(RE.cloc(i1));
 if ~all(i1);
  fprintf('Not implemented cases\n');
  disp(unique(sevLS(~i1)));
  if 1==2
   % DT = delaunayTriangulation(P)
   model=evalin('caller','model');def=evalin('caller','def');
   cf=feplot(model);lsutil('viewls',cf,def);set(cf.ga,'clim',[-1 1]);
  end
  if 1==2
   i3=find(~i1)';Node=evalin('caller','mo3.Node');
   for j3=1:length(i3)
    curelt=elt(i3(j3),feval(ElemP,'node'));
    curprop=elt(i3(j3),feval(ElemP,'prop'));
    i4=RO.idnew(ismember(RO.edges,sort(curelt(feval(ElemP,'edges')),2),'rows'));
    i4(i4==0)=[];i4=[curelt i4];
    n1=Node(RO.NNode(i4),5:7);DT=delaunayTriangulation(n1);
    figure(1);clf;tetramesh(DT,'FaceAlpha',.3);drawnow;pause
    r1=reshape(i4(DT.ConnectivityList),[],4);
    r1=[r1 repmat(curprop,size(r1,1),1)];assignin('caller','r1',r1);
    evalin('caller', sprintf('%s=[%s;r1];','tetra','tetra'));
    %cf.sel={['innode',sprintf('%i ',curelt(:))],'colorfacew -alpha.5'};
    %lsutil('viewls',cf,def);set(cf.ga,'clim',[-1 1]);
    %sevLS(i3(1))
    % Need to cut edges, then do AutoDelaunay for each part
   end
  end
 end
else
 [in1,RE.cloc]=ismember(sevLS,RE.cases,'rows');
 i1=RE.cloc~=0; RE.cloc(i1)=RE.cindex(RE.cloc(i1));
end

RE.icase=setdiff(unique(RE.cloc),0);
inode=feval(ElemP,'node');
for j3=RE.icase(:)' % Existing cases
 in1=RE.cloc==j3; % Elements with current cut type
 
 curelt=elt(in1,inode);curprop=elt(in1,feval(ElemP,'prop'));
 i4=size(RE.CutEdges{j3},1);if i4>0;curelt(end,end+i4)=0;end%add zeros for each cut edge to put new nodes
 for j4=1:i4 % i4=number of cut edges
  [t1,i6]=ismember(sort(curelt(:,RE.CutEdges{j3}(j4,:)),2),idedges,'rows');
  if ~any(i6)||any(RO.idnew(i6)==0);
   disp(RE.cases(RE.cindex==j3))
   cf=feplot(evalin('caller','model'));
   cf.sel={['innode',sprintf('%i ',curelt(:))],'colorfacew -alpha.5'};
   fecom('shownodemark',sort(curelt(:,RE.CutEdges{j3}(j4,:)),2))
   error('Problem');
  end
  curelt(:,length(inode)+j4)=RO.idnew(i6)';
 end
 for j4=1:size(RE.newelt,2)
  st=RE.newelt{j3+1,j4};
  if ~isempty(st);RE.newelt{j3+1,j4}=eval(st);
  else;RE.newelt{j3+1,j4}=[];
  end
 end
 for j4=1:size(RE.iso0,2)
  st=RE.iso0{j3+1,j4};
  if ~isempty(st);RE.iso0{j3+1,j4}=st(curelt);
  else;RE.iso0{j3+1,j4}=[];
  end
 end
end
for j4=1:size(RE.newelt,2)
 % Append to various elements : tetra/pyra/penta
 r1=RE.newelt(2:end,j4);r1(cellfun(@ischar,r1))=[];
 st=RE.newelt{1,j4}; if strcmpi(st,'beam1');st='Beam1';end
 assignin('caller','r1',r1);
 evalin('caller', ...
  sprintf('%s=[%s;vertcat(r1{:})];',st,st));
end

%% #iso_sel : feplot isosurface selection -2
function [out,out1]=iso_sel(varargin) %#ok<DEFNU>

if nargin>1&&ischar(varargin{1})
 [CAM,Cam]=comstr(varargin{1},1);carg=2;
 
 if comstr(Cam,'init')
 %% #iso_sel.Init
 model=varargin{carg};carg=carg+1;
 if isfield(model,'cf');
   cf=comgui('guifeplot',model.cf); RO=model;model=cf.mdl.GetData;
 end
 if isa(model,'sdth'); cf=model;model=cf.mdl.GetData;end
 model.Node=feutil('getnodeGroupall',model);
 [EGroup,nGroup]=getegroup(model.Elt);
 sel=struct('cut',struct('elt',[],'map',[],'mode','clim')); st1='';
 NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
 for jGroup=1:nGroup
  [ElemF,i1,ElemP]= getegroup(model.Elt(EGroup(jGroup),:),jGroup);
  cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
  if isfield(cf.data,'isoMap')&&~isempty(cf.data.isoMap)
   sel.cut.map=cf.data.isoMap;
  elseif ~strcmpi(st1,ElemP);
   RE=feval(lsutil('@iso_zero'),ElemP,struct('silent',1));
   if jGroup>1;dbstack; keyboard;end
   sel.cut.map=RE.isoMap;st1=ElemP; 
   cf.data.isoMap=sel.cut.map;
  end
  elt=full(NNode(model.Elt(cEGI,feval(ElemP,'nodes'))'));
  sel.cut.elt(1:size(elt,1),end+(1:size(elt,2)))=elt;
 end
 sel.Node=model.Node(:,1);sel.vert0=model.Node(:,5:7);
 sel.cna={fe_c(cf.def.DOF,feutil('getdof',(1:3)'/100,model.Node(:,1)))};
 sel.fvcs='r1=r1(:,3);';
 sel.opt=[3 0 1 0 0 0];% NodeBasedColor,tTag, just surface
 sel.off=-.5; sel.step=.01; % default middle of clim
 sel.ScaleColorMode='one'; sel.DispFcn=@isoSurf;
 if ~isempty(cf)
  cf.vfields.SelF{2}=sel; 
  %isoSurf(sel,cf);
 end
 if nargout>0; out=sel;end
 elseif comstr(Cam,'vh8')
  %% iso_sel.vh8 view to help generating cuts
  % feval(lsutil('@iso_sel'),'vh8','++++---+')
  st1=varargin{carg};carg=carg+1;
  def=struct('def',abs(st1)-abs('+')-1);
  def.def=def.def(:);def.DOF=(1:length(def.def))'+.01; 
  def.lab={st1};
  if length(def.def)==8; ElemP='hexa8';
  else; dbstack; keyboard;
  end
  cf=feplot(2);cf.model=femesh(['test' ElemP]);
  cf.def=def; 
  fecom(cf,';ShowFiCevalX;showpatch; ColorAlpha.5');fecom('scc 1e-10');
  fecom('textnode','groupall','fontsize',10)
  EC=integrules(ElemP);i1=feval(ElemP,'edges');st2=sort(st1(i1),2); %#ok<FVAL>
  i2=i1(st2(:,1)=='+'&st2(:,2)=='-',:);r2=[];
  n1=(EC.xi(i2(:,1),:)+EC.xi(i2(:,2),:))/4+.5;
  st2=cellfun(@(x)sprintf('%i',x+length(st1)),num2cell(1:size(n1,1)),'uni',0);
  text(n1(:,1),n1(:,2),n1(:,3),st2,'fontsize',15);
  
 elseif comstr(Cam,'def')
  %% #iso_sel.def initialize view
  RO=varargin{carg};carg=carg+1;
  cf=comgui('guifeplot',RO.cf);
  sel=cf.SelF{2};
  def=feutilb('placeindof',sel.Node+rem(RO.def.DOF(1),1),RO.def);
  sel.fvcs=def.def;sel.ScaleColorMode='one'; 
  if ~isfield(RO,'off')
    r1=[min(sel.fvcs(:)) max(sel.fvcs(:))]; 
    RO.off=mean(r1); if ~isfield(RO,'step');RO.step=diff(r1)/50;end
  end
  sel.off=-RO.off; sel.step=RO.step; % value of iso surface and scroll step %sel.off=-1e8
  cf.SelF{2}=sel; 
  cf.SelF{1}=sdsetprop(cf.SelF{1},'fsProp','facealpha',.1);
  cf.o(1)='sel 1 def 1 ty1'; % mesh
  cf.o(2)='sel 2 def 1 ty10'; % mesh
  fecom(cf,'scaleequal');
  return
 else; error('%s',CAM);
 end
elseif nargin==3&&isequal(varargin{3},'Scroll'); 
 obj=varargin{1};evt=varargin{2};cf=get(obj,'userdata');
 go=findobj(obj.CurrentAxes,'type','patch','tag','iso');
 if isempty(go); return; end
 sel=go.UserData; 
 sel.off=sel.off+evt.VerticalScrollCount*evt.VerticalScrollAmount*sel.step;
 cf=get(ancestor(go,'figure'),'userdata');feplot(cf);
  %isoSurf(sel,cf);
end

%% #isoSurf : generation of a surface patch for cut -2
function ob=isoSurf(sel,cf,j1)

elt=sel.cut.elt;  map=sel.cut.map; 
    fs=ones(size(elt,1),4); js=0;  % Prealloc
    vert=sel.vert0; jn=0; 
ua=cf.ua;ga=cf.ga; 
% [CrntObj(=1 to use wire) CrntSel ObjStamp ObjTy]
sel.opt(1,3)=1; % Force surface view
if ~isfield(sel,'step');sel.step=sel.off/50;end 
[r1,ua,lab,sel,col]= feval(feplot('@get_def_vertices'),cf,ua,sel,[1 2 0 1]);
% def=cf.DefF{crd};def=sel.cna{1}*real(def.def(:,ch)*uaob(9)); % observe color
% r1=get(ga,'clim'); def=(def-r1(1))/diff(r1);
if isempty(col.fsProp); dbstack; keyboard;end

def=col.fsProp{1,2}; 
if strcmpi(sel.ScaleColorMode,'one') % Use a physical offset
 def=def+sel.off; cf.data.IsoSurf=-sel.off; 
else
 if ~isfinite(col.clim(1));col.clim=[min(def) max(def)];end
 def=(def-col.clim(1))/diff(col.clim)+sel.off;
 RO.clim=[min(def) max(def)];
 if RO.clim(1)<0||RO.clim(2)>1; 
  sdtw('_nb','color problem %s',comstr([col.clim;RO.clim],-30));
  col.clim=RO.clim; 
  def=(def-col.clim(1))/diff(col.clim)+sel.off;
 end
end
    st0=repmat('0',1,size(elt,1));
    for jelt=1:size(elt,2)
     curelt=elt(:,jelt);
     st=st0; r3=def(curelt); st(r3>0)='+';st(r3<0)='-';
     if ~map.isKey(st); 
       i2=nnz(st=='-'); i3=nnz(st=='+');st1=st;
       if i2>i3; st(st=='0')='+'; elseif i3>i2; st(st=='0')='-'; end
       if map.isKey(st)
       elseif i2 && i3; 
         st2=st;if st2(1)=='-'; i4=st2=='+';i3=st2=='-';st2(i4)='-';st2(i3)='+';end
         fprintf('%i feval(lsutil(''@iso_sel''),''vh8'',''%s'')\n',jelt,st2);
         continue,
       else
        continue;
       end
     end
     r2=map(st); 
     curelt=curelt';
     if ~isempty(r2.CutEdges);
      n1=r1(curelt(r2.CutEdges(:,1)),:); % r1 current node pos
      n2=r1(curelt(r2.CutEdges(:,2)),:);
      r=r3(r2.CutEdges); r=r(:,1)./(r(:,1)-r(:,2));
      nc=[r1(curelt,:);diag(1-r)*n1+diag(r)*n2];
      nc=nc(r2.elt,:);
     else;
      nc=r1(curelt(r2.elt),:);
     end
     i1=1:size(nc,1);
     vert(jn+i1,:)=nc;
     fs(js+(1:size(r2.elt,2)),:)=reshape(jn+i1,size(fs,2),[])';
     js=js+size(r2.elt,2);jn=jn+size(nc,1);
    end

ob=ua.ob(j1,1);
if ~ob||~ishandle(ob)||~strcmp(get(ob,'type'),'patch')
  ob=patch('tag','iso','parent',ga);
end
setappdata(ga,'SdtScroll',{@iso_sel,'Scroll'});
if all(fs(1,:)==1)&&all(fs(:)==1)
 sdtw('_nb',sprintf('Empty isoSurf not visible, .off=%.1f, dclim(%.1f:%.1f)', ...
     sel.off,min(def),max(def)));
end
%ob=ua.ob(2,1); 
%if ~ishandle(ua.ob(2,1))||~isa(ua.
set(ob,'vertices',vert,'faces',fs,'facecolor','w', ...
    'edgecolor','k','edgealpha',.01,'tag','iso', ...
    'userdata',sel);
out=ob; 
%set(ga,'axis','auto')

%% #iso_zero -2
function [RE]=iso_zero(ElemF,RO) %#ok<DEFNU>
% RE=feval(lsutil('@iso_zero'),'quad4')
%re=feval(lsutil('@iso_zero'),'quad4')
if ~exist('ElemF','var');ElemF='tria3';end
if nargin<2&&nargout==0; RO=struct('cf',sdth.urn('feplot(2)'));cf=RO.cf;
elseif nargin==1||~isfield(RO,'cf');RO.cf=[];cf=[];
end
[RC,RE]=combiLS(ElemF); % All cases of current 
EC=integrules(RE.ElemP);
map=containers.Map({'tetra','penta','tria','hexa','pyra','quad'}, ...
 {'tetra4','penta6','tria3','hexa8','pyra5','quad4'});
if ~isempty(cf)
 cingui('plotwd',cf,'@OsDic(SDT Root)',{'FnTM','ImToFigN','ImSw80','WrW33c'});
end
i1=unique(RE.cindex);
allcases=cell(1+size(i1,1),3);allcases(1,:)={'Beam1','tria','quad'};
for j1=1:length(i1)
 mo1=struct('Node',[(1:size(EC.xi,1))'*[1 0 0 0] EC.xi],'Elt',[]);
 j2=find(RE.cindex==i1(j1),1,'first');
 i2=RE.CutEdges{i1(j1)};
 if ~isempty(i2)
  n2=(EC.xi(i2(:,1),:)+EC.xi(i2(:,2),:))/2;
  mo1.Node(end+(1:size(n2,1)),[1 5:7])=[(1:size(n2,1))'+size(mo1.Node,1) n2];
 end
 curprop=[1 1]; curelt=mo1.Node(:,1)';
 
 for j3=1:size(RE.newelt,2) % Add new elements
  if isempty(RE.newelt{i1(j1)+1,j3});continue;end
  mo1=feutil('addelt',mo1,map(RE.newelt{1,j3}),eval(RE.newelt{i1(j1)+1,j3}));
 end
 if isempty(mo1.Elt);mo1=feutil('addelt',mo1,RE.ElemP,feval(RE.ElemP,'node'));  end
 mo1.name=RE.cases{j2};lsn=abs(mo1.name);lsn(lsn==48)=0;lsn(lsn==45)=-1;lsn(lsn==43)=1;
 mo1.pl=[m_elastic('dbval 1 steel');m_elastic('dbval 2 aluminum')];mo1.il=p_solid('dbval 1 d3 -3');
 def=struct('DOF',(1:size(mo1.Node,1))'+.98,'def',zeros(size(mo1.Node,1),1));
 def.def(1:size(EC.xi,1))=lsn;%lsutil('viewls',mo1,def);

 %fix bug with quadb
 [EGroup,nGroup]=getegroup(mo1.Elt);curcase=cell(nGroup,3);
 for jGroup=1:nGroup
  [ElemF,~,ElemP]= getegroup(mo1.Elt(EGroup(jGroup),:),jGroup);
  elt=mo1.Elt(EGroup(jGroup)+1:EGroup(jGroup+1)-1,feval(ElemP,'node'));
  switch ElemP
  case {'tria3','quad4'}; type1='edges';
  case {'tetra4','pyra5','penta6','hexa8'}; type1='faces';
  otherwise
   dbstack;keyboard
  end
  in1=feval(ElemP,type1); 
  edges=reshape(elt(:,in1')',size(in1,2),[])';
  [un1,un2,i3]=unique(sort(edges,2),'rows');edges=edges(i3,:);
  %feplot(mo1);fecom('textnode')
  lsE=def.def(edges);
  if strcmpi(type1,'faces')&&size(lsE,2)==4
   in2=find(sum(double(~lsE),2)==3);% faces with a triangle remove on nodes
   for j2=in2(:)'; 
     in3=find(lsE(j2,:)~=0); in4=in3+1; in4(in4>4)=1;
     if length(unique(edges(j2,:)))==4; % Only for non triangular
      edges(j2,in3)=edges(j2,in4);
     end
   end
   lsE=def.def(edges); in1=all(~lsE,2);
  else
   in1=all(~lsE,2);
  end
  if all(~in1);% cases that are not handled
   if ~isfield(RO,'silent');  fprintf([mo1.name,'\n']);continue;end
  end
  curcase(jGroup,size(edges,2)-1)={edges(in1,:)};
 end
 
 for type1=1:3;
  st=unique(vertcat(curcase{:,type1}),'rows')';
  if length(st)==1; dbstack; keyboard;
  elseif ~isempty(st)
   allcases{j1+1,type1}=struct('CutEdges',RE.CutEdges{j1},'elt',st);
   %str1=sprintf(sprintf('curelt(:,[%s]);',repmat('%d ',1,type1+1)),st);
   %allcases(j1+1,type1)={sprintf('[%s]',str1(1:end-1))};
  end
  %allcases(j1+1,type1)={vertcat(curcase{:,type1})};
 end

if nargout==0
 for type1=1:3
 if ~isempty(allcases{j1+1,type1});mo1=feutil('addelt',mo1,allcases{1,type1},allcases{j1+1,type1});end
 end
 cf.model=mo1;fecom('colordatamat -alpha.3');cf.os_('LgMdlName')
 if max(mo1.Node(:,7));fecom('view3');end;feplot;fecom('textnode');
 fprintf('Orient for view %s\n',mo1.name); pause(0.1);
 %sdth.urn('feplot(2).Report');
 %cf=feplot(mo1);cf.model=mo1;fe_quality(cf.mdl);fe_quality('view',cf);keyboard
end
end % Loop on cases 
RE.iso0=allcases;
r1=RE.iso0(RE.cindex+1,3);i1=cellfun(@(x)~isempty(x),r1);
r1=r1(i1); st3=RE.cases(i1);
for j1=1:length(r1); % Keep one instance of each face
  [un1,i2]=unique(sort(r1{j1}.elt,1)','rows'); 
  r1{j1}.elt=r1{j1}.elt(:,i2);
end
RE.isoMap=containers.Map(RE.cases(i1),r1);
%% attempt at auto fill missing
st1=setdiff(RC.lsC,RE.isoMap.keys);

for j1=1:length(st1)
 if j1==1; 
  i1=feval(RE.ElemP,'edges');EC=integrules(RE.ElemP);
  EC.faces=feval(RE.ElemP,'faces');
 end
 st2=sort(st1{j1}(i1),2);%st3=st1{j1}(EC.faces)
 i2=i1(st2(:,1)=='+'&st2(:,2)=='-',:);r2=[];
 st2=st1{j1}; if st2(1)=='-'; i4=st2=='+';i3=st2=='-';st2(i4)='-';st2(i3)='+';end
 if 1==2% size(i2,1)==3 % triangle
   r2=struct('CutEdges',i2,'elt',size(EC.xi,1)+[1:3 3]');
 elseif any(strcmp(st2,{'++++++-+','+++++++-','+++++-++','++++-+++', ...
         '+-------','+-++++++','++-+++++','+++-++++'}))
   r2=struct('CutEdges',i2,'elt',size(EC.xi,1)+[1:3 3]');
 %elseif size(i2,1)>4; disp(st2)
 elseif any(strcmp(st2,{'++-+----','--+-++++'}))
   r2=struct('CutEdges',i2,'elt',[9 10 13 12;13 12 11 11]');
 elseif strcmp(st2,'+++--+--'); r2=struct('elt',[9 10 11 12;11 12 14 13]');
 elseif any(strcmp(st2,{'+-++---+'}))
   r2=struct('CutEdges',i2,'elt',[9 10 12 11;11 12 13 14]');
 elseif strcmp(st2,'+---++++');r2=struct('elt',[9 10 13 11;11 13 12 12]');
 elseif strcmp(st2,'++++-+--');r2=struct('elt',[9 10 11 11;9 10 12 11]');
 elseif strcmp(st2,'++++---+');r2=struct('elt',[9 10 12 13;12 10 11 11]');
 elseif strcmp(st2,'++-+++--');r2=struct('elt',[9 10 11 12;12 11 13 13]');
 elseif strcmp(st2,'++-++--+');r2=struct('elt',[9 10 12 11;10 12 13 13]');
 elseif strcmp(st2,'+++++---');r2=struct('elt',[9 10 13 12;10 11 13 13]');
 elseif strcmp(st2,'++++--+-');r2=struct('elt',[9 10 12 13;9 11 13 13]');
 elseif strcmp(st2,'++-++---');r2=struct('elt',[9 11 13 10;10 13 14 12]');
 elseif strcmp(st2,'+++--++-');r2=struct('elt',[9 10 11 12;9 12 13 13]');
 elseif strcmp(st2,'+++-++--');r2=struct('elt',[9 10 13 11;11 13 12 12]');
 elseif strcmp(st2,'+-++--++');r2=struct('elt',[9 10 13 11;10 13 12 12]');
 elseif strcmp(st2,'-+++--+-');r2=struct('elt',[9 10 12 14;9 14 13 11]');
 elseif strcmp(st2,'+-+++--+');r2=struct('elt',[9 10 12 12;10 12 13 11]');
 elseif strcmp(st2,'-+++--++');r2=struct('elt',[9 10 13 13;9 13 12 11]');
 elseif strcmp(st2,'-++--+--');r2=struct('elt',[9 12 13 11;9 11 10 10]');
 elseif strcmp(st2,'++---+--');r2=struct('elt',[9 10 11 12;9 12 13 13]');
 elseif strcmp(st2,'--++---+');r2=struct('elt',[9 10 12 11;10 12 13 13]');
 elseif strcmp(st2,'+--+---+');r2=struct('elt',[10 11 13 12;10 11 9 9]');
 elseif strcmp(st2,'-+++-++-');r2=struct('elt',[9 11 13 12;9 11 10 10]');
 elseif strcmp(st2,'+++-----');r2=struct('elt',[9 10 12 13;10 12 11 11]');
 elseif strcmp(st2,'+-++----');r2=struct('elt',[9 10 12 13;10 12 11 11]');
 elseif strcmp(st2,'+++++-+-');r2=struct('elt',[9 11 12 12;10 14 13 13]');
 elseif strcmp(st2,'++++-+-+');r2=struct('elt',[9 10 12 11;9 10 13 14]');
 elseif strcmp(st2,'+++-++-+');r2=struct('elt',[9 10 13 11;10 13 14 12]');
 elseif strcmp(st2,'+++-++-+');r2=struct('elt',[9 10 13 11;10 13 14 12]');
 elseif strcmp(st2,'+---++-+');r2=struct('elt',[9 10 12 14;9 14 13 11]');
 elseif strcmp(st2,'+--++-++');r2=struct('elt',[9 10 13 12;10 13 11 11]');
 elseif strcmp(st2,'+---++--');r2=struct('elt',[9 10 12 11;10 12 13 13]');
 elseif strcmp(st2,'+---+--+');r2=struct('elt',[9 10 13 12;10 11 13 13]');
 elseif strcmp(st2,'++--+++-');r2=struct('elt',[9 10 13 12;9 12 11 11]');
 elseif strcmp(st2,'++--+---');r2=struct('elt',[9 10 13 12;9 12 11 11]');
 elseif strcmp(st2,'+--++---');r2=struct('elt',[9 10 13 12;10 13 11 11]');
 elseif strcmp(st2,'+--+++-+');r2=struct('elt',[9 10 12 11;10 12 13 13]');
 elseif strcmp(st2,'++--++-+');r2=struct('elt',[9 10 11 13;9 13 12 12]');
 elseif strcmp(st2,'+--+++--');r2=struct('elt',[9 11 13 10;10 13 14 12]');
 elseif strcmp(st2,'++-++++-');r2=struct('elt',[9 10 11 11;12 13 14 14]');
 elseif strcmp(st2,'+----+--');r2=struct('elt',[9 10 11 11;12 13 14 14]');
 elseif strcmp(st2,'++--+--+');r2=struct('elt',[10 12 14 13;14 13 11 11]');
 elseif strcmp(st2,'+--+-+++');r2=struct('elt',[9 10 11 15;10 15 13 13;13 15 14 12]');
 elseif strcmp(st2,'+++--+++');r2=struct('elt',[9 10 11 13;9 13 14 12]');
 elseif strcmp(st2,'+++---+-');r2=struct('elt',[9 10 11 12;9 12 13 14]');
 elseif strcmp(st2,'+------+');r2=struct('elt',[9 10 11 11;12 13 14 14]');
 elseif strcmp(st2,'++-+--+-');r2=struct('elt',[9 10 13 13;11 12 15 16;11 16 14 14]');
 elseif strcmp(st2,'++---++-');r2=struct('elt',[10 11 13 14;9 10 11 11;13 14 12 12]'); % ?
 elseif strcmp(st2,'+-++++-+');r2=struct('elt',[9 10 12 14;9 11 13 14]');
 elseif strcmp(st2,'++---+-+');r2=struct('elt',[12 15 16 16;9 10 11 13;9 13 14 14]');
 elseif strcmp(st2,'++---+++');r2=struct('elt',[9 10 11 14;9 14 15 12;15 12 13 13]');
 elseif strcmp(st2,'++-++-++');r2=struct('elt',[9 10 12 12;11 13 14 14]');
 elseif strcmp(st2,'+---+++-');r2=struct('elt',[9 11 12 10;10 12 13 14]');
 elseif strcmp(st2,'+-++-+++');r2=struct('elt',[9 11 14 10;10 14 13 12]');
 elseif strcmp(st2,'++--+-++');r2=struct('elt',[9 10 13 12;11 14 15 15]');
 elseif strcmp(st2,'+-++-++-');r2=struct('elt',[10 12 14 15;9 10 11 13]');
 %elseif any(strcmp(st1{j1},{'-+++--++'}))
 %  r2=struct('CutEdges',i2,'elt',[9 10 13 13;13 9 12 ]');
 elseif ~any(st1{j1}=='0') % Things to cut again
  fprintf('feval(lsutil(''@iso_sel''),''vh8'',''%s'')\n',st2);
 end
 if ~isempty(r2);
  if ~isfield(r2,'CutEdges'); r2.CutEdges=i2;end
  RE.isoMap(st1{j1})=r2; 
  st2=st1{j1};i2=st2=='+';i3=st2=='-';st2(i2)='-';st2(i3)='+';RE.isoMap(st2)=r2;
 end
end

out=RE;

function sel=isoContour(RO,evt); 
%% #isoContour line on a surface 
   val=evt.LevelList; 
   if 1==2 % Just nodes on edges
    sel.f1=int32((1:size(sel.vert0,1))');
    sel.f1Prop={'FaceColor','none','EdgeColor','k','marker','o','linewidth',[2]};
    sel.fvcs=[];sel.opt=[0 29 0 0 1];
   elseif 1==2 
    %% Now edges in elements with two edges allowing straight line
     RO.edges=RO.oedges; RO.values=RO.ovalues-val;
     RO.values(abs(RO.values)<RO.newTol)=0;
     r3=RO.values(:,1).*RO.values(:,2); i3=(r3>=0)|(RO.edges(:,1)==RO.edges(:,2));
     RO.edges(i3,:)=[];RO.values(i3,:)=[];
     r3(i3)=[]; r4=RO.values(:,2)./r3;
     r=RO.values; r=r(:,1)./(r(:,1)-r(:,2)); r(~isfinite(r))=0;
     sel.vert0=diag(sparse(1-r))*RO.Node(RO.edges(:,1),5:7)+diag(sparse(r))*RO.Node(RO.edges(:,2),5:7);
     r2=max(RO.Node(:,1));
     r2=sparse(RO.Node(RO.edges(:,1)),RO.Node(RO.edges(:,2)), ...
      1:size(RO.edges,1),r2,r2);
     RO.iedge=r2+r2';
     if isempty(RO.cEGI); RO.cEGI=find(isfinite(RO.Elt(:,1)));end
     [EGroup,nGroup]=getegroup(RO.Elt);
     sel.f2=zeros(length(RO.cEGI),2);j1=1;
     for jGroup=1:nGroup
      [ElemF,i1,ElemP]= getegroup(RO.Elt(EGroup(jGroup),:),jGroup);
      if strcmpi(ElemP,'SE');continue;end;i3=feval(ElemP,'edges');
      cEGI=intersect(EGroup(jGroup)+1:EGroup(jGroup+1)-1,RO.cEGI);
      n1=sort(reshape(RO.Elt(cEGI,i3(:,1:2)')',2,[]),1)';
      n2=sort(reshape(RO.Node(RO.edges(:,1:2),1),[],2),2);
      i1=reshape(ismember(n1,n2,'rows'),size(i3,1),length(cEGI));
      i2=sum(double(i1),1)<=1; cEGI(i2)=[];i1(:,i2)=[]; 
      for jElt=1:length(cEGI)
       %% should avoid loop xxx
       i6=reshape(RO.Elt(cEGI(jElt),i3(:,1:2)),[],2);  i5=1;
       for j3=1:size(i6,1)
        i4=RO.iedge(i6(j3,1),i6(j3,2));
        if i4~=0;sel.f2(j1,i5)=i4;i5=i5+1;end
        if i5>2;break;end
       end
       j1=j1+1;
      end
     end
     sel.f2=unique(sel.f2,'rows');
    
    sel.f2(~all(sel.f2,2),:)=[];sel.if2=zeros(size(sel.f2,1),1,'int32');
    sel.f2Prop={'FaceColor','none','EdgeColor','k','marker','none','linewidth',2};
    sel.f1Prop={};sel.fsProp={};sel.fvcs=[];sel.opt=[0 29 0 1 0];
   elseif 1==1 % 
    %% Refine based on nodal values
     [EGroup,nGroup]=getegroup(RO.Elt);
     sel.f2=[];j1=1; sel.vert0=[];RO.edges=[];sel.StressObs.r=[];
     for jGroup=1:nGroup
      [ElemF,i1,ElemP]= getegroup(RO.Elt(EGroup(jGroup),:),jGroup);
      if strcmpi(ElemP,'SE');continue;end;i3=feval(ElemP,'nodes');
      cEGI=intersect(EGroup(jGroup)+1:EGroup(jGroup+1)-1,RO.cEGI);
      if isempty(cEGI);continue;end
      elt=reshape(full(RO.NNode(RO.Elt(cEGI,i3))),length(cEGI),length(i3));
      if strcmpi(ElemP,'tria6')||strcmpi(ElemP,'tria3')
          elt(:,4:end)=[];
          st1={'++-',[1 3 2 3];'--+',[1 3 2 3]
               '+-+',[1 2 3 2];'-+-',[1 2 3 2]
               '+--',[1 3 1 2];'-++',[1 3 1 2]
               '0+-',[1 1 2 3];'0-+',[1 1 2 3]
               '+0-',[2 2 1 3];'-0+',[2 2 1 3]
               '+-0',[3 3 1 2];'-+0',[3 3 1 2]
               '00+',[1 1 2 2];'00-',[1 1 2 2]
               '0+0',[1 1 3 3];'0-0',[1 1 3 3]
               '+00',[2 2 3 3];'-00',[2 2 3 3]
               };
          %map=containers.Map(st1(:,1),st1(:,2));
      elseif strcmpi(ElemP,'quad4')||strcmpi(ElemP,'quadb')
          elt(:,5:end)=[];
          st1={'++--',[1 4 2 3];'--++',[1 4 2 3];
              '+--+',[1 2 4 3];'-++-',[1 2 4 3];
              '+++-',[4 1 3 4];'---+',[4 1 3 4];
              '-+++',[1 2 4 1];'+---',[1 2 4 1];
              '+-++',[1 2 2 3];'-+--',[1 2 2 3];
              '++-+',[2 3 3 4];'--+-',[2 3 3 4];
              '++00',[3 3 4 4];'--00',[3 3 4 4]
              '+00+',[2 2 3 3];'-00-',[2 2 3 3]
              '00++',[1 1 2 2];'00--',[1 1 2 2]
              '0++0',[1 1 4 4];'0--0',[1 1 4 4]
              '0+0-',[1 1 3 3];'0-0+',[1 1 3 3]
              '+0-0',[2 2 4 4];'-0+0',[2 2 4 4]
              '0--+',[1 1 3 4];'0++-',[1 1 3 4];
              '+0--',[2 2 4 1];'-0++',[2 2 4 1];
              '-+0-',[3 3 1 2];'+-0+',[3 3 1 2];
              '--+0',[4 4 2 3];'++-0',[4 4 2 3];
              '+--0',[4 4 1 2];'-++0',[4 4 1 2];
              '0+--',[1 1 2 3];'0-++',[1 1 2 3];
              '-0+-',[2 2 3 4];'+0-+',[2 2 3 4];
              '--0+',[3 3 4 1];'++0-',[3 3 4 1];
              };
          %setdiff(unique(st),st1(:,1))
      else;error('Not implemented yet')
      end
      r2=RO.def.def(elt)-val; r2(abs(r2)<RO.newTol)=0;st=charLs(r2);
      [i1,i2]=ismember(st,st1(:,1));
      i2(~i1)=[];st(~i1)=[];cEGI(~i1)=[];r2(~i1,:)=[];elt(~i1,:)=[];
      i3=(vertcat(st1{i2,2})-1)*size(RO.Elt,1);
      i5=RO.Elt(repmat(cEGI,1,size(i3,2))+i3);% first segment, second segment
      i5=full(RO.NNode(i5));
      [i7,un1,i6]=unique(sort(reshape(i5',2,[]))','rows');i6=reshape(i6,2,[])'+size(sel.vert0,1);
      RO.edges=[RO.edges;i7];
      r=RO.def.def(i7)-val; r=r(:,1)./(r(:,1)-r(:,2)); r(~isfinite(r))=0;
      sel.vert0=[sel.vert0;
        diag(sparse(1-r))*RO.Node(i7(:,1),5:7)+diag(sparse(r))*RO.Node(i7(:,2),5:7)];
      sel.StressObs.r=[sel.StressObs.r;r];
      sel.f2=[sel.f2;i6];
      % fecom('shownodemark',sel.vert0,'marker','o','color','r')
     end
     sel.f2=unique(sort(sel.f2,2),'rows');
    
    sel.f2(~all(sel.f2,2),:)=[];sel.if2=zeros(size(sel.f2,1),1,'int32');
    sel.f2Prop={'FaceColor','none','EdgeColor','k','marker','none','linewidth',2};
    sel.f1Prop={};sel.fsProp={};sel.fvcs=[];sel.opt=[0 29 0 1 0];
    
   end
   sel.StressObs.EdgeN=RO.Node(RO.edges);
   if isfield(RO,'ToFace')&&~isempty(sel.f2)
    %% #isContour.Coarse : coarsen  -3
    [un1,i1]=min(std(sel.vert0));idir=setdiff(1:3,i1);
    if length(RO.ToFace)<2
      i2=sparse(sel.f2(:,1),sel.f2(:,2),1);%% Fastversion of LineLoops
      i2(end+1:size(i2,2),1)=0;i2(1,end+1:size(i2,1))=0;conn=i2+i2'; 
      i2=sdtm.feutil.k2PartsVec(conn);
    else
      %% coarsen lines provide LC, cosMin 
      for j3=1:4
      i2=sparse(sel.f2(:,1),sel.f2(:,2),1);%% Fastversion of LineLoops
      i2(end+1:size(i2,2),1)=0;i2(1,end+1:size(i2,1))=0;conn=i2+i2'; 
      i2=sdtm.feutil.k2PartsVec(conn);
      if length(RO.ToFace)<3; RO.ToFace(2:3)=[30 .9]; end
      i3=conn; i3(:,sum(conn)~=2)=0;
      [II,JJ]=find(i3);JJ(1:2:end)=[];% 
      r1=sel.vert0(II(1:2:end),:)-sel.vert0(JJ,:); l1=sqrt(sum(r1.^2,2));
      r2=sel.vert0(JJ,:)-sel.vert0(II(2:2:end),:); l2=sqrt(sum(r2.^2,2));
      [r3,i3]=sortrows([-round(sum(r1.*r2,2)./l1./l2*1000) sort([l1 l2],2)]);
      i4=[reshape(II,2,[]);JJ']';i4=i4(i3,:);
      r3(r3(:,1)>-1000*RO.ToFace(3)|r3(:,3)>RO.ToFace(2),:)=[];i4(size(r3,1)+1:end,:)=[];
      % Sorted possibly edge sequence, remove center points
      for j1=1:size(i4,1)
        if ~all(i4(j1,:)); continue;end % Already merged
        l1=norm(sel.vert0(i4(j1,1),:)-sel.vert0(i4(j1,3),:));
        l2=norm(sel.vert0(i4(j1,2),:)-sel.vert0(i4(j1,3),:));
        if l1==0||l2==0||l1+l2>RO.ToFace(2); continue
          %if l1+l2<RO.ToFace(2)&&r3(j1)<-990;% further merge, is wrong
          %  i6=setdiff(i4(j1,:),i4(j1,1)); if isempty(i6);continue;end
          %  sel.f2(sel.f2==i6)=i4(j1,1); i4(i4==i6)=i4(j1,1);
          %else; continue;
          %end
        elseif l1<l2; sel.f2(sel.f2==i4(j1,3))=i4(j1,1); i4(i4==i4(j1,3))=i4(j1,1);
        else;  sel.f2(sel.f2==i4(j1,3))=i4(j1,2); i4(i4==i4(j1,3))=i4(j1,2);
        end
        sel.f2(sel.f2(:,1)==sel.f2(:,2),:)=[];
      end
      end
      i3=unique(sel.f2(:)); nind=sparse(i3,1,1:length(i3));
      sel.vert0=sel.vert0(i3,:);sel.f2=full(nind(sel.f2));
      sel.if2=zeros(size(sel.f2,1),1,'int32');
      sel.StressObs=struct('EdgeN',sel.StressObs.EdgeN(i3,:),'r',sel.StressObs.r(i3));

      %disp(length(unique(sel.f2))/size(sel.vert0,1))
    end
    %% #isContour.2DFill : generate a face  -3
    try
     i2=sdtm.feutil.k2PartsVec(conn);% find connected nodes
     r3=sparse(sel.f2,i2(sel.f2),1);
     [i3,i4]=find(r3==1); 
     r3(:,unique(i4))=[]; i3=find(any(r3',1));% Remove open loops
    
     nind=sparse(i3,1,1:length(i3),size(sel.vert0,1),1);
     i4=full(nind(sel.f2)); i4(any(i4==0,2),:)=[];
     warning('off','MATLAB:triangulation:EmptyTri2DWarnId');
     if ~isempty(i3)
      DT=delaunayTriangulation(sel.vert0(i3,idir(1)),sel.vert0(i3,idir(2)),i4);
     else; DT=struct('Points',[1]);
     end
     if isempty(DT.ConnectivityList); sel.fs=[];
     elseif size(DT.Points,1)==length(i3)&&~isempty(sel.f2)
      i4=DT.ConnectivityList(DT.isInterior,:);
      sel.fs=reshape(i3(i4),size(i4));
      sel=sdsetprop(sel,'fsProp','edgecolor','none','facecolor','white');
     else; 
      sel.fs=[];
     end
    catch; 
      sel.fs=[];
    end
    warning('on','MATLAB:triangulation:EmptyTri2DWarnId');
    
    % cf=feplot;sdtm.feutil.SelPatch(sel,cf.ga,'reset')
   else; sel.fs=[];
   end
    % remove unused nodes
    i3=find(ismember(1:size(sel.vert0),sel.f2(:)));
    nind=sparse(i3,1,1:length(i3),size(sel.vert0,1),1);
    sel.f2=full(nind(sel.f2)); 
    if ~isempty(sel.fs);sel.fs=reshape(full(nind(sel.fs)),size(sel.fs));end
    sel.vert0=sel.vert0(i3,:);
    if ~isfield(evt,'ProId')
    elseif ~isempty(sel.fs)
      sel.mdl=struct('Node', ...
        [(RO.Nend+(1:size(sel.vert0))')*[1 0 0 0] sel.vert0], ...
        'Elt',feutil('addelt','tria3',[sel.fs+RO.Nend ones(size(sel.fs,1),2)*evt.ProId]));
      sel.ifs=(2:size(sel.fs,1)+1)';
    else
      sel.mdl=struct('Node', ...
        [(RO.Nend+(1:size(sel.vert0))')*[1 0 0 0] sel.vert0], ...
        'Elt',feutil('addelt','beam1',[sel.f2+RO.Nend ones(size(sel.fs,1),2)*evt.ProId]));
      sel.if2=(2:size(sel.fs,1)+1)';
    end
    %cf=feplot;sdtm.feutil.SelPatch(sel,cf.ga,'reset')
1;

%% #cutCurElt 
function out=cutCurElt(curelt,i1,pro);

if isequal(i1,0); out=[];
elseif size(i1,2)==1; out=curelt(:,i1);
else
    dbstack; keyboard;
end

%% #addContourCutLine : adds a line in a single cut contour
function addContourCutLine %#ok<DEFNU>

eval(iigui({'RO','Beam1','idedges'},'GetInCaller'))

if size(RO.values,1)~=2;
 error('Not an expected case');
end
n1=[0 0];
for j1=1:2
 i1=RO.values(j1,:)==0;
 if nnz(i1);n1(j1)=idedges(j1,i1);
 else; n1(j1)=RO.idnew(j1);end
end
assignin('caller','Beam1',n1);
if isempty(Beam1)
 evalin('caller','[mo3,RO]=safeAddElt(mo3,RO,''beam1'',Beam1);');
else
 evalin('caller','mo3.Elt(end+(1:size(Beam1,1)),1:2)=Beam1;');
end


%% #combiLS: compute all cases of cut
function [out,out1]=combiLS(ElemF)
%RE=feval(lsutil('@combiLS'),'tetra4')
if nargin==0;ElemF='pyra5';end
%create grid
v=[-1,0,1];nd=3;
switch lower(ElemF)
 case {'tria3','tria'}
  [out1,out2,out3]=ndgrid(v);nn=3;RE=reTria;nd=2;
 case {'quad4','tetra4','quad'}
  [out1,out2,out3,out4]=ndgrid(v);nn=4;
  if strfind(ElemF,'quad');RE=reQuad;nd=2;else;RE=reTetra;end
 case {'pyra5','pyra'}
  [out1,out2,out3,out4,out5]=ndgrid(v);nn=5;RE=rePyra;
 case {'penta6','penta'}
  [out1,out2,out3,out4,out5,out6]=ndgrid(v);nn=6;RE=rePenta;
 case {'hexa8','hexa'}
  [out1,out2,out3,out4,out5,out6,out7,out8]=ndgrid(v);nn=8;RE=reHexa;
 otherwise
  warning('element %s not supported',ElemF);dbstack;keyboard;
end
if ~iscell(RE.cases);RE.cases=charLs(RE.cases);end
%compute combi
val=zeros(3^nn,nn);
for j1=1:nn;val(:,j1)=eval(sprintf('reshape(out%i,[],1,1);',j1));end
%find combi that need to be cut
icut=(any(val'==-1)&any(val'==1))';
%
allcases=val(icut,:);lsc=charLs(allcases);
[a,ia,ib]=intersect(allcases,-allcases,'rows','stable');
ipair=unique(sort([ia ib],2),'rows');
%
dval=val(~icut,:);
[a,ia,ib]=intersect(dval,-dval,'rows','stable');
in0=unique(sort([ia ib],2),'rows');
in1=in0(:,1)~=in0(:,2);in0=in0(in1,:);% remove only 0
in1=sum(dval(in0(:,1),:)==0,2)>=nd;in0=in0(in1,:);
%
out1=RE;out1.cases(end+(1:numel(in0)),1)=charLs(dval(in0',:));
out1.cindex=reshape([1;1]*(1:size(out1.cases,1)/2),[],1);
out1.newelt{1+size(out1.cases,1)/2,end}=[];
out1.CutEdges{size(out1.cases,1)/2,end}=[];
%output : .todo list remaining cases 
out=struct('ElemP',RE.ElemP,'all',val,...
 'lsN',allcases,'lsC',{lsc},'only',ipair,...
 'done',{RE.cases},'todo',{setdiff(lsc(ipair(:,1),:),RE.cases)});

%% #check_split
function [model,out1,out2]=check_split(model,def)

% this is only implemented with distFcn present to avoid problems
%check split need
[eltid,model.Elt]=feutil('eltidfix;',model);
mo3=struct('Node',[],'Elt',[]);mo3.Node=model.Node;
[EGroup,nGroup]=getegroup(model.Elt);listE=[];
NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
for jGroup=1:nGroup
 [ElemF,opt,ElemP] = feutil('GetElemF',model.Elt(EGroup(jGroup),:),jGroup);
 if opt(1)<0;continue;end
 if strcmp(ElemP,'hexa8');RE=reHexa;
 elseif strcmp(ElemP,'penta6');RE=rePenta;
 else;continue
 end
 cEGI=(EGroup(jGroup)+1:EGroup(jGroup+1)-1)';prop=feval(ElemP,'prop');
 elt=full(NNode(model.Elt(cEGI,feval(ElemP,'node'))));sevLS=sign(def.def(elt));
 in1=any(sevLS==-1,2)&any(sevLS==1,2)&~any(isnan(sevLS),2);
 if any(in1)
  celt=cEGI(in1);lsC=charLs(sevLS(in1,:));[in2,RE.cloc]=ismember(lsC,RE.cases);
  if any(~in2)
   [mo3,RA]=safeAddElt(mo3,struct('nextId',size(mo3.Elt,1)+1),ElemF,model.Elt(celt(~in2),:));
   listE=[listE;celt(~in2)];
  end
 end
end

out1=isempty(listE);
if ~isempty(mo3.Elt);
 n1=mo3.Node(mo3.Node(:,1)==max(mo3.Node(:,1)),:);
 mo3.Node=feutil('getnode groupall',mo3);
 if max(mo3.Node(:,1))~=n1(1);mo3.Node(end+1,:)=n1;end
 [mo4,un1,RO]=lsutil('split',mo3,def,struct('iter',1));
 n2=mo4.Node(~ismember(mo4.Node(:,1),mo3.Node(:,1)),:);
 r1=def.distFcn(n2(:,5:7));
 model.Elt=feutil('removeelt eltind',model,listE);
 if size(def.def,1)~=size(model.Node,1);error('Not expected');end
 def.DOF=[def.DOF;n2(:,1)+.98];def.def=[def.def;r1];
 model.Node=[model.Node;n2];  
 model=feutil('addelt',model,mo4.Elt);
 [eltid,model.Elt]=feutil('eltidfix;',model);
 %if ~any(r1);out1=def; return;end% did all cuts
end
out2=def;

%[model,un1,RO]=lsutil('split',model,li,RO);
%RO=feutil('rmfield',RO,'conn','edges','values','NNode','ed2elt','elt2ed');
%RO.conn=feval(feutilb('@levNodeCon'),[],model,'econ');

