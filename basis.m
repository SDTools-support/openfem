function [out,out1,out2]=basis(varargin)

%BASIS  coordinate transformation utilities
%
%	Syntax: see below
%
%       p = basis(x,y) P is the orthonormal 3 by 3 matrix whose first column
%          is in the direction of X, second column is along Y (but orthogonal
%          to X), and third column forms a direct orthonormal basis
%
%	   If x and y are collinear y is selected along the smallest component
%	   of x. A message is passed unless a third argument exists (call of 
%          the form basis(x,y,1))
%
%      [nodeGlob,bas]=basis(nodeLocal,bas) ...=basis('nodebas',model)
%          transforms nodeLocal to a standard nodeGlobal node definition matrix. 
%          BAS is transformed to eliminate all recursive definitions.
%
%       nodeGlobal = basis('gnode',bas,nodeLocal) given an single coordinate
%          system definition BAS, associated nodes nodeLocal are transform to
%          the global coordinate system.
%
%       cGL = basis('trans [ ,l][ ,t][,e]',bas,node,DOF) given a set of 
%          coordinate system definitions BAS, the displacement coordinate 
%          defintions in NODE(:,3) are used to create a local to global 
%          coordinate DOF transformation. The modifier L is used to specify 
%          that nodes are given in local rather than global coordinates.
%
%       [p,nodeL]  = basis(node) With 2 outputs, NODE gives coordinates of a 
%          QUAD4 element, P is the local element basis, and NODEL the node 
%          positions in the element coordinate system
%
%       bas=basis('rotate',bas,rot,id)
%       node=basis('cyl2rect',node); cylindrical to rect transform
%       node=basis('rect2cyl',node); rectangular to cylindrical transform
%       [xl yl zl]=basis('bunge',[phi_1 Phi phi_2]) % Crystallography angles
%
%	See sdtweb   basis, node, elt


%	Etienne Balmes
%       Copyright (c) 2001-2021 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license.

if ischar(varargin{1});
 [CAM,Cam]=comstr(varargin{1},1);x=[];carg=2;
else; x=varargin{1}; carg=2;Cam='';
end
%#ok<*NASGU,*ASGLU,*NOSEM>

%% #basis_construction (classical function) ------------------------------------
if numel(x)==3  
 %% should really call sp_util('basis',a,b)
 y = varargin{2}; 

 x=x(:);y=y(:);
 x=x/sqrt(x'*x);
 y=y-(x'*y)*x;
 if norm(y)<eps*norm(x)
   if nargin==2; warning('OpenFEM:BASCOL','collinear vectors'); end
   [y,i1]=min(abs(x)); y=zeros(3,1);y(i1(1))=1;y=y-(x'*y)*x; 
 end
 y=y/sqrt(y'*y);
 z=-[y(2)*x(3)-y(3)*x(2);y(3)*x(1)-y(1)*x(3);y(1)*x(2)-y(2)*x(1)];
 out=[x y z];

%% #element_basis_given_faces --------------------------------------------------
elseif isempty(Cam)&&nargin==1 %& size(x,1)==4  

 if size(x,2)==7;     x=x(:,5:7);
 elseif size(x,2)~=3; error('Not a valid node matrix');end

 if size(x,1)>3&&size(x,1)~=6 % quad (4 or b)
   if size(x,1)==8;x=x(1:4,:);end
   [out,out1,out2]=sp_util('basis',x');

  % out2 = mean(x,1);node = x(1:4,:)-out2(ones(4,1),:);
  % x = node(3,:)-node(1,:); x=x/sqrt(x*x');
  % y = node(4,:)-node(2,:); y=y/sqrt(y*y');
  % x = x-y; x=x'/sqrt(x*x');
  % y=y';y=y-(x'*y)*x; y=y/sqrt(y'*y);
  % z=-[y(2)*x(3)-y(3)*x(2);y(3)*x(1)-y(1)*x(3);y(1)*x(2)-y(2)*x(1)];
  % p=[x y z]; out=p; out1=node*p; x=node;

   if nargout>2 %return surface
      FEnar=[-.25 .25 .25 -.25]; FEnas=[-.25 -.25 .25 .25];
      xr = FEnar*x(:,1); xs = FEnas*x(:,1);
      yr = FEnar*x(:,2); ys = FEnas*x(:,2);
      zr = FEnar*x(:,3); zs = FEnas*x(:,3);
      jdet = sqrt((xr.*ys-xs.*yr)^2+(yr.*zs-ys.*zr)^2+(zr.*xs-zs.*xr)^2);
      out2=jdet;
   end

 elseif size(x,1)==3||size(x,1)==6 % triangle
   out2 = mean(x,1); 
   node = x-out2(ones(size(x,1),1),:);
   x = node(2,:)-node(1,:); x=x/sqrt(x*x');
   y = node(3,:)-node(2,:); r1=norm(y);y=y/r1;
   x = (x-y)'; r1=norm(x); 
   if r1==0; [r1,i1]=min(abs(y)); x(i1)=1; else; x=x/r1;end
   y=y';y=y-(x'*y)*x; 
   r1=norm(y); if r1==0; [r1,i1]=min(abs(x)); y(i1)=1; else; y=y/r1;end
   z=-[y(2)*x(3)-y(3)*x(2);y(3)*x(1)-y(1)*x(3);y(1)*x(2)-y(2)*x(1)];
   p=[x y z]; 
   out=p; out1=node*p;
   if nargout>2 %return surface
      FEw=[1/3 1/3 0 1/2];  % ASSUMES 3 NODE  or straight 6 node TRIANGLE
      xi = [0 0 0;1 0 0;0 1 0];
      FEnar=ones(size(FEw,1),1)*[-1 1 0];
      FEnas=ones(size(FEw,1),1)*[-1 0 1];
      xr = FEnar*node(1:3,1); xs = FEnas*node(1:3,1);
      yr = FEnar*node(1:3,2); ys = FEnas*node(1:3,2);
      zr = FEnar*node(1:3,3); zs = FEnas*node(1:3,3);
      jdet = sqrt((xr.*ys-xs.*yr)^2+(yr.*zs-ys.*zr)^2+(zr.*xs-zs.*xr)^2);
      out2=jdet;
   end
 end

 %% #Bas -----------------------------------------------------------
elseif comstr(Cam,'bas')
 bas=varargin{carg};
 model=struct('Node',[],'bas',bas);
 [un1,bas]=basis('nodebas-force',model);
 out=bas;

%% #NodeBas -Force : nodes from local to global ------
elseif comstr(Cam,'nodebas')||(nargin==2&&size(x,2)==7)

if ~isempty(Cam);
 model=varargin{carg};carg=carg+1; 
 if isnumeric(model); FEnode=model;bas=varargin{carg};carg=carg+1;
 elseif ~isempty(stack_get(model,'SE')); 
     [FEnode,bas]=fesuper(['fnode' Cam(8:end)],model);
 elseif ~isfield(model,'Node');error('Missing node');
 else
   FEnode=model.Node; if isfield(model,'bas');bas=model.bas;else;bas=[];end
 end
else;
 FEnode=x;  bas=varargin{carg};carg=carg+1;
end
% types : CORD2R 1, CORD2C 2, CORD2S 3, CORD1R 4, CORD1C 5 , CORD1S 6 
% CORD Final is [CID RID(0) Type Origin ex ey ez]

if isa(bas,'cell') 
  r1=zeros(length(bas),1);
  for j2=1:length(bas);  r1(j2,1:length(bas{j2}))=bas{j2}; end
elseif isa(bas,'double'); r1=bas;
else; error('Not a proper basis definition');
end
if ~isempty(strfind(Cam,'force'))
    if isempty(bas);out1=bas; out=FEnode;return;end
elseif isempty(bas)||isempty(FEnode)||~any(FEnode(:,2)); 
    out1=bas; out=FEnode;return;
end

NNode=[];if size(r1,2)<15; r1(1,15)=0; end

% work the coordinate changes to global coordinates

if ~isempty(r1) 

 i4=(1:size(r1,1))'; i5=find(i4&r1(:,3)==0); % systems that can be transformed
 i1=length(i5); 
 while(any(i4)) % not all have been transformed
  r0=r1;
  if isempty(FEnode); % Calls from basis('bas') to resolve recursive bas def
   nind=sparse(max(max(r1(:,1)),max(i5))+2,size(FEnode,1));
  else
   nind=sparse(FEnode(:,2)+2,1:size(FEnode,1),FEnode(:,2), ...
     max(max(r1(:,1)),max(i5))+2,size(FEnode,1));
  end
  if nnz(nind)==0;nind=[];end
 for j1 = i5(:)' % loop on systems defined from global

   i2=r1(j1,:);
   % go to everything given form if needed
   if i2(2)==10
    error('Use [model.Node,model.bas]=feutilb(''NodeBas'',model) for ANSYS format (not basis)');
    % [node,bas]=feutilb('nodebas',model.Node,model.bas);
   elseif any(i2(7:12))&&~any(i2(13:15)) % Ai Bi Ci given
     r1(j1,:)=abc2bas(i2);
     i4(j1)=0;
   elseif ~any(i2(7:15)) % node numbers are given
     if isempty(NNode);NNode=sparse(FEnode(:,1),1,1:size(FEnode,1));end
     i2=FEnode(NNode(r1(j1,4:6)),2);
     i2=[find(r1(:,1)==i2(1)) find(r1(:,1)==i2(2)) find(r1(:,1)==i2(3))];
     if isempty(i2)|| ...
        (~any(r1(i2,3))&&all(i4(i2)==0)) % refering to transformed basis
      x = FEnode(NNode(r1(j1,4:6)),5:7);
      r3 = basis(x(2,:)-x(1,:),x(3,:)-x(1,:)); r3 = r3(:,[2 3 1]);
      r1(j1,1:15) = [r1(j1,1:2)  0 x(1,:) r3(:)'];
      i4(j1)=0;
     end
   else; i4(j1)=0; % everything is given
   end

   % transform nodes that depend on this system
   if i4(j1)==0&&~isempty(nind)
     ind=find(nind(r1(j1,1)+2,:)); %ind=find(FEnode(:,2)==r1(j1,1));
     %ind=[];if r1(j1,1)<=size(indb,1) ind=find(indb(r1(j1,1),:));end
     if ~isempty(ind)
        FEnode(ind,:)=basis('gnode',r1(j1,:),FEnode(ind,:));
        FEnode(ind,2)=0;
     end
   end

   % transform systems that depend on this one
   for j3= find(r1(:,3)==r1(j1,1))'
    i2=r1(j3,:);
    if any(i2(7:12))&&~any(i2(13:15)) % Ai Bi Ci given
     r2=basis('gnode',r1(j1,:),reshape(r1(j3,4:12),3,3)');%A,B,C nodes in glob
     i2=[0 0 0 r2(1,5:7)  r2(2,5:7)  r2(3,5:7)];
     i3 = basis(i2(7:9)-i2(4:6),i2(10:12)-i2(4:6)); i3 = i3(:,[2 3 1]);
     r1(j3,3:15)=[0 r2(1,5:7) i3(:)'];
    elseif  any(i2(7:9))&&any(i2(13:15))&&any(i2(10:12))% p given
     i3=bas2cGL(r1(j1,:)); %i3=bas2cGL(i2); % Reference basis 
     r1(j3,3:15)=[0 r1(j1,4:6)+i2(1,4:6)*i3' reshape(i3*bas2cGL(i2),1,[])];
    else
     error('Not a supported dependent basis transformation');
    end
   end

 end % j1 : systems defined from global
 i5=find(i4&r1(:,3)==0);
 if i1==length(i5)&&isequal(r0,r1);
  sdtw('Reference Coord Systems %s are not defined', ...
      sprintf(' %i',unique(r1(r1(:,3)~=0,3))));
  break
 else; i1=length(i5);
 end
 end
end

out1=r1; out=FEnode;

%% #GNode transform nodes given coordinate definition ------------------
elseif comstr(Cam,'gnode')

% base definition
% [Id Type RefId OriginXYZ Xvector Yvector Zvector]
r1  = varargin{carg};carg=carg+1; node=varargin{carg};carg=carg+1;
%try;out=sp_util('gnode',r1,node);return;end

if ischar(r1);r1=basis('rotate',[],r1,1);end
if size(node,2)==3; node=[zeros(size(node,1),4) node];end

if size(r1,2)<16||r1(16)==0; scale=1; else;scale=r1(16);end
if size(r1,2)<15||~any(r1(13:15)); error('Basis format.'); end

if r1(2)==1 % rectangular
  node(:,5:7)=r1(ones(size(node,1),1),4:6)+ ...
           node(:,5:7)*reshape(r1(7:15)*scale,3,3)';
elseif r1(2)==2 % CORD2C
  node(:,5:7)=r1(ones(size(node,1),1),4:6)+ ...
           [node(:,5).*cos(node(:,6)/180*pi) ...
            node(:,5).*sin(node(:,6)/180*pi) node(:,7)] ...
            *reshape(r1(7:15)*scale,3,3)';
elseif r1(2)==3 % CORD2S
  node(:,6:7)=node(:,6:7)*pi/180;
  node(:,5:7)=r1(ones(size(node,1),1),4:6)+ node(:,[5 5 5]).* ...
           [sin(node(:,6)).*cos(node(:,7)) sin(node(:,6)).*sin(node(:,7)) ...
            cos(node(:,6))] ...
            *reshape(r1(7:15)*scale,3,3)';
else;error('Not a valid type');
end
out=node;

%% #Trans ----------------------------------------------------------------------
% trans build the local to global coordinate transformation matrix
elseif comstr(Cam,'trans')

[CAM,Cam]=comstr(CAM,6); if isempty(Cam); Cam=' ';end
[CAM,Cam,i1]=comstr('-force',[-25 3],CAM,Cam); % Avoid -force in string

bas=varargin{carg};carg=carg+1;
FEnode=varargin{carg};carg=carg+1;
if carg<=nargin; odof=varargin{carg};carg=carg+1; else;odof=[]; end

if any(Cam=='l') % nodes are given in local coordinates
  [FEnode,bas]=basis(FEnode,bas);
end
if ~isempty(odof)&&any(Cam=='e')
 i2=unique(round(rem(odof,1)*100));
 if isempty(intersect(i2,[4 5 6])); Cam(end+1)='t';end
 FEnode(~ismember(FEnode(:,1),fix(odof(:,1))),:)=[]; % keep used Nodes
end
mdof = FEnode(:,1);
if any(Cam=='t'); mdof=[mdof+.01 mdof+.02 mdof+.03]'; ND=3;
else
   mdof=[mdof+.01 mdof+.02 mdof+.03 mdof+.04 mdof+.05 mdof+.06]';ND=6;
end
mdof=mdof(:);

i3=0;
if ~isempty(bas); if size(bas,2)<15; bas(1,15)=0;end;bas=bas(:,1:15);end
bas=[0 1 0  0 0 0  1 0 0  0 1 0  0 0 1;bas];

i1=find(sparse(FEnode(:,3)+1,1,FEnode(:,3)+1))-1;i1(i1==0)=[]; 

i2=nnz(FEnode(:,3))*3*ND;II=zeros(i2,1); JJ=zeros(i2,1);KK=zeros(i2,1);
ji=[-1 0]; ji(2)=0; jj=[-1 0]; jj(2)=0; jk=[-1 0]; jk(2)=0;
for j1=i1(:)'  % loop on coordinates
  r1=find(bas(:,1)==j1); 
  if isempty(r1); sdtw('Coordinate sytem %i is undefined',j1);
  elseif length(r1)>1; error('Repeated BasID');
  else
   r1 = bas(r1,:);
   if ~any(r1(13:15))&&r1(2)~=10;  
    r1=abc2bas(r1);
    %error('Basis definitions must be transformed use ''TransL''');
   end
   tr= reshape(r1(7:15),3,3);i2 = find(FEnode(:,3)==j1);
   switch r1(2)
   case 1 % rectangular
     ind=ND*(i2-1)+1;tr=tr(:)';tr=tr(ones(size(ind)),:);
     sp_util('setinput',II,ind  ,ji,'II');sp_util('setinput',JJ,ind,jj,'JJ');
     sp_util('setinput',II,ind+1,ji,'II');sp_util('setinput',JJ,ind,jj,'JJ');
     sp_util('setinput',II,ind+2,ji,'II');sp_util('setinput',JJ,ind,jj,'JJ');
     sp_util('setinput',II,ind  ,ji,'II');sp_util('setinput',JJ,ind+1,jj,'JJ');
     sp_util('setinput',II,ind+1,ji,'II');sp_util('setinput',JJ,ind+1,jj,'JJ');
     sp_util('setinput',II,ind+2,ji,'II');sp_util('setinput',JJ,ind+1,jj,'JJ');
     sp_util('setinput',II,ind  ,ji,'II');sp_util('setinput',JJ,ind+2,jj,'JJ');
     sp_util('setinput',II,ind+1,ji,'II');sp_util('setinput',JJ,ind+2,jj,'JJ');
     sp_util('setinput',II,ind+2,ji,'II');sp_util('setinput',JJ,ind+2,jj,'JJ');
     sp_util('setinput',KK,tr,jk,'KK');
     if ND==6; ind=ind+3;
     sp_util('setinput',II,ind  ,ji,'II');sp_util('setinput',JJ,ind,jj,'JJ');
     sp_util('setinput',II,ind+1,ji,'II');sp_util('setinput',JJ,ind,jj,'JJ');
     sp_util('setinput',II,ind+2,ji,'II');sp_util('setinput',JJ,ind,jj,'JJ');
     sp_util('setinput',II,ind  ,ji,'II');sp_util('setinput',JJ,ind+1,jj,'JJ');
     sp_util('setinput',II,ind+1,ji,'II');sp_util('setinput',JJ,ind+1,jj,'JJ');
     sp_util('setinput',II,ind+2,ji,'II');sp_util('setinput',JJ,ind+1,jj,'JJ');
     sp_util('setinput',II,ind  ,ji,'II');sp_util('setinput',JJ,ind+2,jj,'JJ');
     sp_util('setinput',II,ind+1,ji,'II');sp_util('setinput',JJ,ind+2,jj,'JJ');
     sp_util('setinput',II,ind+2,ji,'II');sp_util('setinput',JJ,ind+2,jj,'JJ');
     sp_util('setinput',KK,tr,jk,'KK');
     elseif 1==2
     %ind=ND*(i2-1)+1;in1=i3+[1:length(ind)];
     %T1(in1,1:2)=[ind   ind  ];T1(in1,3)=tr(1);in1=in1+length(ind);
     %T1(in1,1:2)=[ind+1 ind  ];T1(in1,3)=tr(2);in1=in1+length(ind);
     %T1(in1,1:2)=[ind+2 ind  ];T1(in1,3)=tr(3);in1=in1+length(ind);
     %T1(in1,1:2)=[ind   ind+1];T1(in1,3)=tr(4);in1=in1+length(ind);
     %T1(in1,1:2)=[ind+1 ind+1];T1(in1,3)=tr(5);in1=in1+length(ind);
     %T1(in1,1:2)=[ind+2 ind+1];T1(in1,3)=tr(6);in1=in1+length(ind);
     %T1(in1,1:2)=[ind   ind+2];T1(in1,3)=tr(7);in1=in1+length(ind);
     %T1(in1,1:2)=[ind+1 ind+2];T1(in1,3)=tr(8);in1=in1+length(ind);
     %T1(in1,1:2)=[ind+2 ind+2];T1(in1,3)=tr(9);
     elseif ND==6
      %in1=in1+length(ind);ind=ind+3;
      %T1(in1,1:2)=[ind   ind  ];T1(in1,3)=tr(1);in1=in1+length(ind);
      %T1(in1,1:2)=[ind+1 ind  ];T1(in1,3)=tr(2);in1=in1+length(ind);
      %T1(in1,1:2)=[ind+2 ind  ];T1(in1,3)=tr(3);in1=in1+length(ind);
      %T1(in1,1:2)=[ind   ind+1];T1(in1,3)=tr(4);in1=in1+length(ind);
      %T1(in1,1:2)=[ind+1 ind+1];T1(in1,3)=tr(5);in1=in1+length(ind);
      %T1(in1,1:2)=[ind+2 ind+1];T1(in1,3)=tr(6);in1=in1+length(ind);
      %T1(in1,1:2)=[ind   ind+2];T1(in1,3)=tr(7);in1=in1+length(ind);
      %T1(in1,1:2)=[ind+1 ind+2];T1(in1,3)=tr(8);in1=in1+length(ind);
      %T1(in1,1:2)=[ind+2 ind+2];T1(in1,3)=tr(9);
     end
     %i3=in1(end);
   case 2 % cylindrical
    for j2=i2(:)'
     r2=(FEnode(j2,5:7)-r1(4:6))*tr;
     r3=atan2(r2(2),r2(1)); 
     r3=tr*[cos(r3) -sin(r3) 0;sin(r3) cos(r3) 0;0 0 1];
     ind=ND*(j2-1)+1;
     sp_util('setinput',II,ind  ,ji,'II');sp_util('setinput',JJ,ind,jj,'JJ');
     sp_util('setinput',II,ind+1,ji,'II');sp_util('setinput',JJ,ind,jj,'JJ');
     sp_util('setinput',II,ind+2,ji,'II');sp_util('setinput',JJ,ind,jj,'JJ');
     sp_util('setinput',II,ind  ,ji,'II');sp_util('setinput',JJ,ind+1,jj,'JJ');
     sp_util('setinput',II,ind+1,ji,'II');sp_util('setinput',JJ,ind+1,jj,'JJ');
     sp_util('setinput',II,ind+2,ji,'II');sp_util('setinput',JJ,ind+1,jj,'JJ');
     sp_util('setinput',II,ind  ,ji,'II');sp_util('setinput',JJ,ind+2,jj,'JJ');
     sp_util('setinput',II,ind+1,ji,'II');sp_util('setinput',JJ,ind+2,jj,'JJ');
     sp_util('setinput',II,ind+2,ji,'II');sp_util('setinput',JJ,ind+2,jj,'JJ');
     sp_util('setinput',KK,r3,jk,'KK');
     if ND==6; ind=ind+3;
     sp_util('setinput',II,ind  ,ji,'II');sp_util('setinput',JJ,ind,jj,'JJ');
     sp_util('setinput',II,ind+1,ji,'II');sp_util('setinput',JJ,ind,jj,'JJ');
     sp_util('setinput',II,ind+2,ji,'II');sp_util('setinput',JJ,ind,jj,'JJ');
     sp_util('setinput',II,ind  ,ji,'II');sp_util('setinput',JJ,ind+1,jj,'JJ');
     sp_util('setinput',II,ind+1,ji,'II');sp_util('setinput',JJ,ind+1,jj,'JJ');
     sp_util('setinput',II,ind+2,ji,'II');sp_util('setinput',JJ,ind+1,jj,'JJ');
     sp_util('setinput',II,ind  ,ji,'II');sp_util('setinput',JJ,ind+2,jj,'JJ');
     sp_util('setinput',II,ind+1,ji,'II');sp_util('setinput',JJ,ind+2,jj,'JJ');
     sp_util('setinput',II,ind+2,ji,'II');sp_util('setinput',JJ,ind+2,jj,'JJ');
     sp_util('setinput',KK,tr,jk,'KK');
     end
     % T1(i3+[1:9],1:3) = [ND*(j2-1)+1+[0;1;2;0;1;2;0;1;2] ...
     %                     ND*(j2-1)+1+[0;0;0;1;1;1;2;2;2] r3(:)];
     % i3=i3+9;
     % if ND==6 
     %   T1(i3+[1:9],1:3) = [ND*(j2-1)+4+[0;1;2;0;1;2;0;1;2] ...
     %                       ND*(j2-1)+4+[0;0;0;1;1;1;2;2;2] r3(:)];
     % i3=i3+9;
     % end
    end
   case 3 % spherical
     error('Spherical coordinates where never tested')
%      for j2=i2(:)'
%       r2=(FEnode(j2,5:7)-r1(4:6))*tr;
%       r3=atan2(r2(2),r2(1)); 
%       r3=tr*[cos(r3) -sin(r3) 0;sin(r3) cos(r3) 0;0 0 1];
%       ind = (j2-1)*6+[1:3]; T1(ind,ind)=r3; ind=ind+3;T1(ind,ind)=r3;
%      end
   otherwise; error('Unsupported coordinate transformation');
  end
 end % the basis was found
end
%T1=sparse(II,JJ,KK,length(mdof),length(mdof));
IIn=find(FEnode(:,3)==0)*ND;
if ND==3; IIn=[IIn;IIn-1;IIn-2];else; IIn=[IIn;IIn-1;IIn-2;IIn-3;IIn-4;IIn-5];end
%T1=T1+sparse(IIn,IIn,1,length(mdof),length(mdof));
if isempty(II); T1=sparse(IIn,IIn,1,length(mdof),length(mdof));
else
 sp_util('setinput',II,IIn,ji,'II');
 sp_util('setinput',JJ,IIn,jj,'JJ');
 sp_util('setinput',KK,ones(length(IIn),1),jk,'KK');
 clear IIn
 T1=sparse(II(1:ji(2)),JJ(1:jj(2)),KK(1:jk(2)),length(mdof),length(mdof)); clear II JJ KK
end


if      isempty(odof); out=T1;
elseif length(odof)==length(mdof) && norm(odof-mdof)<.01; out=T1;
elseif any(Cam=='e'); %with a return a square transformation
    ind=fe_c(mdof,odof,'ind'); out=T1(ind,ind);
    if length(ind)<length(odof) % some DOFs were not rotated
        in1=fe_c(odof,mdof(ind),'ind');
        [II,JJ,KK]=find(out); II=in1(II); JJ=in1(JJ);
        out=sparse(II,JJ,KK,length(odof),length(odof));
        in1=setdiff(1:length(odof),in1);
        out=out+sparse(in1,in1,1,length(odof),length(odof));
    end
else;out=fe_c(odof,mdof,T1);
end
out1=mdof;

%% #Rect2Cyl -------------------------------------------------------------------
% node=basis('rect2cyl',node); rectangular to cylindrical transform
elseif comstr(Cam,'rect2cyl')

 node=varargin{carg};r1=[];
 if ~isempty(strfind(Cam,'fix'))
    if size(node,2)==3;i2=1:2;else;i2=5:6;end     
    r1=mean(node(:,i2));r1=atan2(r1(2),r1(1));  
    r2=[cos(r1) sin(r1);-sin(r1) cos(r1)]';
    node(:,i2)=node(:,i2)*r2;r1=r1/pi*180;
 end
 if size(node,2)==3; 
  node(:,1:2)=[sqrt(sum(node(:,1:2).^2,2)) atan2(node(:,2),node(:,1))/pi*180];
 elseif size(node,2)==7
  node(:,5:6)=[sqrt(sum(node(:,5:6).^2,2)) atan2(node(:,6),node(:,5))/pi*180];
 else;error('Not a valid node size');
 end
 if isempty(r1);elseif size(node,2)==3;node(:,2)=node(:,2)+r1;
 else;node(:,6)=node(:,6)+r1;
 end
 out=node;

%% #Cyl2Rect -----------------------------------------------------------------
% node=basis('cyl2rect',node); cylindrical to rect transform
elseif comstr(Cam,'cyl2rect')

 node=varargin{carg};
 if size(node,2)==3; 
  r1=node(:,2)*pi/180;node(:,1:2)=node(:,[1 1]).*[cos(r1) sin(r1)];
 elseif size(node,2)==7
  r1=node(:,6)*pi/180;node(:,5:6)=node(:,[5 5]).*[cos(r1) sin(r1)];
 else;error('Not a valid node size');
 end
 out=node;
 
%% #Rotate : rot uses 'tx,rx ' -----------------------------------------------
% node=basis('rotate',bas,rot,id); rotation of a given basis
elseif comstr(Cam,'rotate')
    
   bas=varargin{carg};carg=carg+1;
   r1=varargin{carg};carg=carg+1;
   i2=varargin{carg};carg=carg+1;
   if isstruct(bas); model=bas;bas=model.bas;  % accept model 
     if size(bas,2)<15;[node,bas]=basis('nodebas-force',model);end
   end
   if isempty(bas); i1=[];else; i1=find(ismember(bas(:,1),i2));end
   if isempty(i1); 
       for j1=1:length(i2); 
           bas(end+1,1:15)=[i2(j1) 1 0  0 0 0 1 0 0 0 1 0 0 0 1]; %#ok<AGROW>
       end
       i1=find(ismember(bas(:,1),i2));
   end

   for j1=1:length(i1)
    if bas(i1(j1),2)~=1; error('Rectangular is the only supported here');end
    if bas(i1(j1),3)~=0; error('basis must be fully resolved');end
    r2=reshape(bas(i1(j1),7:15),3,3);
    if ~all(any(r2,1)); error('All three basis vectors expected');end
    tx=0;ty=0;tz=0; % possibliy define translation
    if ischar(r1)
      rx=0;ry=0;rz=0; eval([r1 ';']);
      if exist('r','var') % r n : random axis
         r=r*pi/180; % angle
         n=n/norm(n); % axis
         r3=cos(r)*eye(3)+...
           (1-cos(r))*[n(1)^2,n(1)*n(2),n(1)*n(3); n(1)*n(2),n(2)^2,n(2)*n(3); n(1)*n(3),n(2)*n(3),n(3)^2]+...
            sin(r)*[0 -n(3) n(2); n(3) 0 -n(1); -n(2) n(1) 0]; % rotation matrix/n
      else %rx ry rz
         rx=rx*pi/180; ry=ry*pi/180; rz=rz*pi/180;
         r3=rotZYX(rx,ry,rz); 
      end
      r2=r2*r3;
    elseif all(size(r1)==3); r2=r2*r1; 
    end % do rotation
    bas(i1(j1),7:15)=r2(:)';
    bas(i1(j1),4:6)=bas(i1(j1),4:6)+[tx ty tz]; % origin at opposite
   end
   out=bas;
%% #Bunge ----------------------------------------------------------------------
% Material angles used in cristalography
% cGL = [xl yl zl]=basis('bunge',[phi_1 Phi phi_2])  % radians
elseif comstr(Cam,'bunge')

 r1=varargin{carg};
 out=[cos(r1(1)) -sin(r1(1)) 0;sin(r1(1)) cos(r1(1)) 0;0 0 1]*...%phi_1
     [1 0 0; 0 cos(r1(2)) -sin(r1(2)); 0 sin(r1(2)) cos(r1(2))]* ... %Phi
     [cos(r1(3)) -sin(r1(3)) 0;sin(r1(3)) cos(r1(3)) 0;0 0 1]; % phi_2
 out(abs(out)<1e-13)=0;
 
%% #PlaneFromN: generate full basis from a vector, thus adding a basis of the plane to a given normal
elseif comstr(Cam,'planefromn'); [CAM,Cam]=comstr(CAM,11);
 r1=varargin{carg}; carg=carg+1;
 
 out=cross(r1,ones(size(r1,1),1)*[1 1 1]);
 i1=~any(out,2); % is colinear to 1 1 1, redo other value
 if any(i1);  out(i1,:)=cross(r1(i1,:),ones(sum(i1),1)*[1 0 0]); end
 if ~isempty(CAM) % only one component required
  out=out(:,abs(Cam(1))-119);
 end
 
%% #struct2row : robust structure input 
elseif comstr(Cam,'struct2row')
    
 RO=varargin{carg};carg=carg+1;
 model=varargin{carg};carg=carg+1;
 
 if isfield(RO,'pxb');  r3=RO.pxb-RO.ori;RO.xb=r3/norm(r3);end
 if isfield(RO,'pyb');  r3=RO.pyb-RO.ori;RO.yb=r3/norm(r3);end
 if isfield(RO,'pzb');  r3=RO.pzb-RO.ori;RO.zb=r3/norm(r3);end
 st=sort(intersect(fieldnames(RO),{'xb','yb','zb'}));
 if isempty(st); RO.xb=[1 0 0];RO.yb=[0 1 0];st={'xb','yb'};
     warning('Assuming global basis');
 end
st=sprintf('%s:%s',st{1:2});out=[];
 
 switch st
 case 'xb:yb'
  RO.zb=cross(RO.xb,RO.yb);
 case 'xb:zb'; 
  RO.yb=cross(RO.zb,RO.xb);
 otherwise; error('%s not yet implemented',st)
 end
 if ~isfield(RO,'RelTo')||strcmpi(RO.RelTo,'B0')
 else
   bas=sdth.urn(sprintf('nmap.bas.%s',RO.RelTo),model);
   %bas=sdth.urn(sprintf('nmap.bas.%s',RO.RelTo),model);
   bas(2,1:12)=[bas(1)+1 1 bas(1) RO.ori(:)' RO.ori(:)'+RO.zb(:)' RO.ori(:)'+RO.xb(:)'];
   [node,bas]=basis('nodebas-force',struct('Node',[],'Elt',[],'bas',bas));
   out=bas(2,:);   
 end
 if isempty(out);out=[0 1 0 RO.ori(:)' RO.xb(:)' RO.yb(:)' RO.zb(:)'];end
 out1=RO;
%% #LRtraj ----------------------------------------------------------------------
% large angle trajectory as linear combination of constant shapes
elseif comstr(Cam,'lrtraj')
    
model=varargin{carg};carg=carg+1;
RO=varargin{carg};carg=carg+1;
n1=model.Node;
if isfield(RO,'UsedNodes')
  n1(~ismember(n1(:,1),RO.UsedNodes),:)=[];
end
if size(RO.def,1)~=6; error('Expecting master node translation and rotation');end
if ischar(RO.rotFcn);RO.rotFcn=eval(RO.rotFcn);end% '@rotZYX'
qLR=RO.def;qLR(12,1)=0;
for j1=1:size(RO.def,2)
  rot=RO.rotFcn(RO.def(4:6,j1))-eye(3);
  qLR(4:end,j1)=rot(:);
end
RO.def=qLR; 

RO.TR=struct('def',[],'DOF',feutil('getdof',n1(:,1),(1:3)'/100));
i1=size(n1,1);in1=1:3:3*i1;
RO.TR.def=sparse(repmat([in1;in1+1;in1+2],4,1), ...
  [(1:12)'*ones(1,i1)], ...
  [ones(3,i1);n1(:,[5 5 5 6 6 6 7 7 7])']);
RO.DOF=fix(RO.DOF(1))+[1:3 13:21]'/100; 
out=RO;

%% #CVS ------------------------------------------------------------------------
elseif comstr(Cam,'@');out=eval(CAM);
elseif comstr(Cam,'cvs')
   out='$Revision: 1.82 $  $Date: 2022/02/01 20:06:57 $'; return;
else
 error('Not a valid call')
end

%% #SubFunc --------------------------------------------------------------------

%% #bas2cGL(bas,Node) : [xb yb zb] columns  -3
function cGL=bas2cGL(bas,Node)

if size(bas,1)>1;error('Not an expected case');end
i1=[any(bas(7:9)) any(bas(13:15)) any(bas(10:12))];
if all(i1)
  cGL=reshape(bas(7:15),3,3);
else;
 error('recursive bases requires a definition using U,V,W');
 dbstack; keyboard
end

%% #abc2bas - - ----------------------------------------------------------------
function i2=abc2bas(i2)
 i3 = basis(i2(7:9)-i2(4:6),i2(10:12)-i2(4:6)); 
 i3 = i3(:,[2 3 1]);
 i2(7:15)=i3(:);

%% #rotZYX  using Bryan / aerospace convention angles feval(basis('@rotZYX'),rand(1,3))
function r3=rotZYX(rx,ry,rz); 
  if nargin==3
  elseif length(rx)==6
   ry=rx(5);rz=rx(6);rx=rx(4);
  elseif length(rx)==3
   ry=rx(2);rz=rx(3);rx=rx(1);
  else; error('Invalid call');
  end
  
  cx=cos(rx);sx=sin(rx); cy=cos(ry);sy=sin(ry); cz=cos(rz); sz=sin(rz);
  r3=[cy*cz -sz*cx+cz*sy*sx sz*sx+cz*sy*cx
       sz*cy cz*cx+sz*sy*sx -cz*sx+sz*sy*cx
       -sy  cy*sx cy*cx];
  % r3=[cos(rz) -sin(rz) 0;sin(rz) cos(rz) 0;0 0 1]*...
  %    [cos(ry) 0 sin(ry);0 1 0; -sin(ry) 0 cos(ry)]*...
  %    [1 0 0; 0 cos(rx) -sin(rx); 0 sin(rx) cos(rx)];
  %if norm(r3-r3b)>1e-15;error(1);end

% #irotZYX Inverse of Bryant ang : from rotation matrix to angles. GerCard (4.128)
function ang=irotZYX(r)
 %a=@(x)feval(basis('@irotZYX'),feval(basis('@rotZYX'),x))-x;a(rand(1,3))
 % atan2(sin,cos) 
 
 theta=atan2(r(2,1),r(1,1));cz=cos(theta);sz=sin(theta);
 psi=atan2(-r(3,1),r(2,1)*sz+r(1,1)*cz);
 phi=atan2(r(1,3)*sz-r(2,3)*cz,-r(1,2)*sz+r(2,2)*cz);
 ang=[phi psi theta];
  
%% #rotXYZ  cardan angles feval(basis('@rotZYX'),rand(1,3))
function r3=rotXYZ(rx,ry,rz); 
  if nargin==3
  elseif length(rx)==6
   ry=rx(5);rz=rx(6);rx=rx(4);
  elseif length(rx)==3
   ry=rx(2);rz=rx(3);rx=rx(1);
  else; error('Invalid call');
  end
  
  cx=cos(rx);sx=sin(rx); cy=cos(ry);sy=sin(ry); cz=cos(rz); sz=sin(rz);
  r3=[   cy*cz          -cy*sz            sy;
   cx*sz+sx*sy*cz    cx*cz-sx*sy*sz  -sx*cy;
   sx*sz-cx*sy*cz  sx*cz+cx*sy*sz     cx*cy ];
 
%% #rotZYXt
function r3=rotZYXt(rx,ry,rz);  %#ok<DEFNU>
 if nargin<3||~isstruct(rz)
  if nargin==3
  elseif length(rx)==6
   ry=rx(5);rz=rx(6);rx=rx(4);
  elseif length(rx)==3
   ry=rx(2);rz=rx(3);rx=rx(1);
  else; error('Invalid call');
  end
   r3=([cos(rz) -sin(rz) 0;sin(rz) cos(rz) 0;0 0 1]*...
      [cos(ry) 0 sin(ry);0 1 0; -sin(ry) 0 cos(ry)]*...
      [1 0 0; 0 cos(rx) -sin(rx); 0 sin(rx) cos(rx)])';
 else
   r4=ry;r2=rx;r3=rz;
   if nnz(r4(1:3)) % Used in Oscar translations
       r2=r2+ones(size(r2,1),1)*r4(1:3)';
       r4(1:3)=0; 
   end
   if isfield(r3,'q0');error('Obsolete');end
     %fprintf('Dof %.2f, alpha %.2f\n',def.DOF(r3.ind(2)),r4*180/pi)
   if nnz(r4(4:6)); % large rotation
       r4=feval(r3.RotFcn,r4);
       r2=(r2*r4+r3.vert0*(r4-eye(3)));
   end
   r3=r2;
 end
%% #rotZYXTt Rotation then translation
function r3=rotZYXTt(rx,ry,rz);  %#ok<DEFNU>
 if nargin<3||~isstruct(rz)
  if nargin==3
  elseif length(rx)==6
   ry=rx(5);rz=rx(6);rx=rx(4);
  elseif length(rx)==3
   ry=rx(2);rz=rx(3);rx=rx(1);
  else; error('Invalid call');
  end
   r3=([cos(rz) -sin(rz) 0;sin(rz) cos(rz) 0;0 0 1]*...
      [cos(ry) 0 sin(ry);0 1 0; -sin(ry) 0 cos(ry)]*...
      [1 0 0; 0 cos(rx) -sin(rx); 0 sin(rx) cos(rx)])';
 else
   r4=ry;r2=rx;r3=rz;
   if isfield(r3,'q0');error('Obsolete');end
     %fprintf('Dof %.2f, alpha %.2f\n',def.DOF(r3.ind(2)),r4*180/pi)
   trans=r4(1:3)';
   if nnz(r4(4:6)); % large rotation
       r4=feval(r3.RotFcn,r4);
       r2=(r2*r4+r3.vert0*(r4-eye(3)));
   end
   if nnz(r4(1:3)) % Used in Oscar translations
       r2=r2+ones(size(r2,1),1)*trans;
   end
   r3=r2;
 end

  
%% #spin
function x1=spin(x) 

x1=[0 -x(3) x(2);x(3) 0 -x(1);-x(2) x(1) 0]; % (3.18)
 
%% #largeRot [R,T]=largeRot(psi) see Geradin-Cardona
% some documentation in project/engine.tex 
function [R,T]=largeRot(psi); %#ok<DEFNU>
      npsi=norm(psi);
      if npsi<eps;R=eye(3);T=eye(3);
      else
       a=psi/norm(psi);a=spin(a);b=a*a;
       R=eye(3)+sin(npsi)*a+(1-cos(npsi))*b; % (4.5) or geradin-cardona
       T=eye(3)+(cos(npsi)-1)/npsi*a+(1-sin(npsi)/npsi)*b;
      end

 
% largeRotToEuler       
