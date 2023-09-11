function [out,out1,out2]=integrules(varargin); 

% This function is used to generate integration constants used by elements
%
% integrules(TopologyName,Type,position)  generates integration rules. 
% The second optional argument type lets you select between various 
% integration rules
%
% Supported topologies are
%  bar1 (1D linear)
%  beam1 (1D cubic), beam3 (1D 3 node quadratic)
%  quad4 (2D bi-linear), quadb (2d quadratic)
%  tria3 (2D affine), tria6 (2D quadratic)
%  tetra4, tetra10
%  penta6, penta15
%  hexa8, hexa20, hexa27
%
% Supported Gauss point families can be listed using 
%  integrules('gauss')
% For a given rule family (1d, q2d, t2d, t3d, p3d, h3d), available rules 
%  and ID are listed with integrules('gauss h3d')
%

%	Etienne Balmes, Jean Michel Leclere
%       Copyright (c) 2001-2020 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       integrules('cvs') for revision information

%#ok<*NOSEM,*NASGU>
if nargin==0;
 integrules('gauss');return
end
out = [];
[CAM,Cam]=comstr(varargin{1},1); carg=2;

if comstr(Cam,'cvs') 

 out='$Revision: 1.152 $  $Date: 2023/09/04 13:56:10 $';
 
 
%% #BuildNDN Jacobian computations for Surfaces in 3D ------------------------
elseif comstr(Cam,'buildndn') 

typ=varargin{carg};carg=carg+1;
if any(typ==[2 3 31 23 13 32])&&nargout==1 % 2D/3D NDN building
 opt=varargin{carg};carg=carg+1;
 if carg<=nargin; opt.nodeE=varargin{carg};carg=carg+1;end
 if carg<=nargin;error('Obsolete');end;%of_mk('buildndn',typ,opt,nodeE,varargin{carg:end});
 of_mk('buildndn',typ,opt);
else

opt=varargin{carg};carg=carg+1;
if carg<=nargin;nodeE=varargin{carg};carg=carg+1;else;nodeE=opt.nodeE;end
st1='';
try;
 if typ==23
   if ~isfield(opt,'bas');st1='NeedM';else;st1='Nz';end
 elseif (typ~=32 && (isfield(opt,'bas')&&norm(opt.bas)==0) && typ~=13 );
  st1='NeedM';
 end
 if carg<=nargin;opt.nodeEt=varargin{carg};carg=carg+1;end
 if carg<=nargin;warning('Error TODO please report');end
 opt.nodeE=nodeE;of_mk('buildndn',typ,opt);
 if strcmp(st1,'Nz')&&~any(opt.bas(:));st1='NeedM';end
 if any(~isfinite(opt.NDN(:))); st1='Not finite';end

catch; st1='NeedM';
end
%w=opt;st1='NeedM'
if strcmp(st1,'NeedM');
typ=abs(typ);
Nw=size(opt.N,1); % number fo integration points
Nnode=size(opt.N,2);
if ~isfield(opt,'bas'); opt.bas=zeros(9,size(opt.N,1));end

switch typ
case 23; % surface rule
 opt.NDN(:,1:Nw)=opt.N'; R=zeros(2);
 for jw=1:Nw

  % local basis on the surface, sdtweb t_rigid('offset')
  % [dx_e;dy_e]=J[dr;ds]
  x=opt.Nr(jw,:)*nodeE(:,1:3);xr=norm(x); x=x/xr; xs=0;
  y=opt.Ns(jw,:)*nodeE(:,1:3);  a=-x*y'; y=y+a*x;
  ys=norm(y);yr=-a; y=y/ys;

  if size(nodeE,2)>4 % local orientation v1x, v1y, v1z in 5:7
   x=x(:);y=y(:);
   xe=opt.NDN(:,jw)'*nodeE(:,5:7); R(1)=xe*x; R(2)=xe*y; 
   r1=sqrt(R(1)^2+R(2)^2);R(1)=R(1)/r1; R(2)=R(2)/r1; R(3)=-R(2);R(4)=R(1);
   r1=[x y]*R;x=r1(:,1); y=r1(:,2);
   z=-[y(2)*x(3)-y(3)*x(2);y(3)*x(1)-y(1)*x(3);y(1)*x(2)-y(2)*x(1)];
   % N,xe ye = N,rs * xy,rs * R(xe,x)
   opt.NDN(:,jw+[Nw 2*Nw])= ...
    [opt.Nr(jw,:)' opt.Ns(jw,:)']*[1/xr  -(yr/xr/ys); 0 1/ys]*R;
   opt.bas(:,jw)=[x;y;z];
  else
   z=-[y(2)*x(3)-y(3)*x(2);y(3)*x(1)-y(1)*x(3);y(1)*x(2)-y(2)*x(1)];
   opt.jdet(jw)=xr*ys;
   %J=[1/xr -yr/xr/ys;0 1/ys];Ji=[xr yr;xs ys];
   %opt.NDN(:,jw+(1:2)*Nw)=[opt.Nr(jw,:)' opt.Ns(jw,:)']*J;
   opt.NDN(:,Nw+jw)=opt.Nr(jw,:)'/xr; % N,x(jw)
   opt.NDN(:,2*Nw+jw)=opt.Ns(jw,:)'/ys-opt.Nr(jw,:)'*(yr/xr/ys); % N,x(jw)
   opt.bas(:,jw)=[x';y';z];
  end

 end
case 32 % oriented local basis

 opt.NDN(:,1:Nw)=opt.N';cof=zeros(3,3);

 for jw=1:Nw

  % local basis on the surface
  xr=opt.Nr(jw,:)*nodeE(:,1);
  xs=opt.Ns(jw,:)*nodeE(:,1);
  xt=opt.Nt(jw,:)*nodeE(:,1);
  yr=opt.Nr(jw,:)*nodeE(:,2);
  ys=opt.Ns(jw,:)*nodeE(:,2);
  yt=opt.Nt(jw,:)*nodeE(:,2);
  zr=opt.Nr(jw,:)*nodeE(:,3);
  zs=opt.Ns(jw,:)*nodeE(:,3);
  zt=opt.Nt(jw,:)*nodeE(:,3);
  xe=opt.NDN(:,jw)'*nodeE(:,5:7);
  ye=opt.NDN(:,jw)'*nodeE(:,8:10);
  
      cof(1,1) = ys.*zt - yt.*zs; % xxx transpose
      cof(1,2) = yt.*zr - yr.*zt;
      cof(1,3) = yr.*zs - ys.*zr;
      cof(2,1) = zs.*xt - zt.*xs;
      cof(2,2) = zt.*xr - zr.*xt;
      cof(2,3) = zr.*xs - zs.*xr;
      cof(3,1) = xs.*yt - xt.*ys;
      cof(3,2) = xt.*yr - xr.*yt;
      cof(3,3) = xr.*ys - xs.*yr;

      jdet = xr.*cof(1)+xs.*cof(2)+xt.*cof(3);
      Nrst=[opt.Nr(jw,:)' opt.Ns(jw,:)' opt.Nt(jw,:)'];
      Nxyze=Nrst*(cof*basis(xe,ye)/jdet); % dxe/dx^T=basis(xe,ye);
      opt.NDN(:,jw+[Nw 2*Nw 3*Nw])=Nxyze;
 end

end
%norm(opt.NDN-opt2.NDN)
end % catch
end % Typ 23
out=opt;

if nargout==2
 Nw=size(opt.N,1); % number for integration points
 out1=zeros(Nw*3,opt.Nnode*3);
 % 3 dof per node
 %                  j           k
 for jw=1:Nw; 
   x=opt.bas(1:3,jw);y=opt.bas(4:6,jw);z=opt.bas(7:9,jw);
   J=reshape(opt.J(:,jw),2,2);%J=(J'/J)';
   Nxyt=(opt.NDN(:,(1:2)*Nw+jw)*(J'/J)')';% from of_mk buildndn
   %b=(opt.NDN(:,(1:2)*Nw+jw))'; % integrule value (sign change on jacobian)
   for j1=1:3;
    for j3=1:3
    %                  u_ik                      Ny_i
    out1(j1+3*(jw-1),j3:3:end)= ...
       [-y(j1)*z(j3)+y(j3)*z(j1)/2   x(j1)*z(j3)-x(j3)*z(j1)/2]*Nxyt;
    end
   end
 end
 out1(abs(out1)<2.2e-14)=0;% remove round off errors at 100 eps
end

%% #MITC4 --------------------------------------------------------------------
elseif comstr(Cam,'mitc4')

  OfMkPol=of_mk('mitcinit');
  fprintf('\n%%Written by integrules(''mitc4'')\n\n');
  writetoscreen(OfMkPol)

%% #MIT4 Integration rules for the MITC4 element -----------------------------
elseif comstr(Cam,'mit4') 

 if nargin<2; w='def'; else; w=varargin{carg};carg=carg+1;end
 if isequal(w,'center'); out=integrules('quad4',w); return;end
 % 4 nodes and 4 tying points
 xi = [-1 -1 0;1 -1 0;1 1 0;-1 1 0; ...
       0 1 0;-1 0 0;0 -1 0;1 0 0];

 % points of interest : integration and tying points
   w=1/sqrt(3);
   w=[-w -w 0 1;w -w 0 1;w w 0 1;-w w 0 1;
         0 1 0 0;-1 0 0 0;0 -1 0 0;1 0 0 0];

   r=w(:,1); s=w(:,2); z=zeros(size(r)); o=ones(size(r));

   N=[(1-r).*(1-s)/4 (1+r).*(1-s)/4  (1+r).*(1+s)/4  (1-r).*(1+s)/4  ...
      (1+s)/2 (1-r)/2 (1-s)/2 (1+r)/2];
 
   Nr=[(s-1)/4 (1-s)/4 (1+s)/4 (-1-s)/4 z -o/2 z o/2];
   Ns=[(r-1)/4 (-1-r)/4 (1+r)/4 (1-r)/4 o/2 z -o/2 z];

  out=struct('N',N,'Nr',Nr,'Ns',Ns,'Nw',size(w,1),'NDN',zeros(4,size(N,1)*3), ...
    'jdet',zeros(size(N,1),1),'w',w,'NDNLabels',{{'',',x',',y'}}, ...
   'DofLabels',{{'u','v','w','\theta_u','\theta_v','\theta_w'}}, ...
   'Nnode',4,'xi',xi,'type','mit4','PerField',4,...
   'Gmesh',[Inf 116 114 105 97 51;1 2 3 0 0 0;3 4 1 0 0 0;
    1 4 6 0 0 0;1 2 7 0 0 0;1 8 3 0 0 0;4 3 5 0 0 0],'Gedge',[]);

%% #Strain_ ----------------------------------------------------------------
%% #StrainString -2
% Analyses the strain string to build a strain definition
elseif comstr(Cam,'strainstring')

  opt=varargin{carg};carg=carg+1;
  r3=varargin{carg};carg=carg+1; %  NwStart NwRule
  if ~isfield(opt,'FieldDofs')
      fprintf('Trying to guess FieldDofs');
      opt.FieldDofs=1:length(opt.DofLabels);
      opt.FieldDofs(strcmpi(opt.DofLabels,'p'))=19;
      opt.FieldDofs(strcmpi(opt.DofLabels,'t'))=20;
  end
  opt.StrainDefinition=cell(opt.StrainString);
  for j1=1:length(opt.StrainString)
    st=opt.StrainString{j1}; rule=[];
    for j2=1:length(st)
       st1=st{j2}; coef=1;
       while ~isempty(st1);
        if comstr(st1,'+');coef=1; st1(1)='';continue;end
        if comstr(st1,'+');coef=-1; st1(1)='';continue;end
        i1=regexp(st1,'[+-,]','once'); 
        if isempty(i1); % interpolated field
          rule(end+1,:)=[j2 1*coef strmatch(st1,opt.DofLabels) r3];
          st1='';
        elseif st1(i1)==',';
          rule(end+1,:)=[j2 1 strmatch(st1(1:i1-1),opt.DofLabels) r3];
          st1(1:i1)='';i1=regexp(st1,'[+-,]','once');
          if isempty(i1);i1=length(st1)+1;end
          rule(end,2)=coef*find(strcmp([',' st1(1:i1-1)],opt.NDNLabels));
          st1(1:i1-1)='';
        else error('Not expected');
        end
       end
    end
    opt.StrainDefinition{j1}=rule;
  end
  opt=integrules('matrixrule',opt);
  opt.VectMap=int32(reshape(1:length(opt.FieldDofs)*opt.Nnode,[],opt.Nnode)'); 
  out=opt;

%% #StrainNdn - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'strainndn')

   opt=varargin{carg};carg=carg+1;
   j0=varargin{carg};carg=carg+1;
   r1=opt.StrainDefinition{j0}; 

   DofPerNode=length(opt.DofLabels); Nnode=opt.Nnode;
   Ndof=DofPerNode*Nnode; Nstrain=max(r1(:,1)); Nw=size(opt.N,1);
   Nshape=size(opt.N,2);

   i1=[1 cumprod([Nstrain Nw Ndof])]; i2=[1 cumprod([Nnode Nw 3])];
   c=spalloc(i1(end),i2(end),1);

   for jw=1:Nw; 
    for jStrain=1:size(r1,1); 
    for jNode=1:Nnode
     %             j(epsilon)  jw     jDOF
     c(i1(1:3)*([r1(jStrain,1);jw;r1(jStrain,3)+(jNode-1)*DofPerNode]-1)+1, ...
       i2(1:3)*([jNode;jw;abs(r1(jStrain,2))]-1)+1)=sign(r1(jStrain,2));
    end
    end
   end
   
   in1=varargin{carg};carg=carg+1; % integration points to eliminate
   % Now eliminate strain assembly at the tying points
   ind=[];
   for jw=in1(:)';
       for jStrain=1:Nstrain
           for jNode=1:Nnode; 
    ind=[ind i1(1:3)*([jStrain;jw;1+(jNode-1)*DofPerNode]-1)+1];
    ind=[ind i1(1:3)*([jStrain;jw;2+(jNode-1)*DofPerNode]-1)+1];
    ind=[ind i1(1:3)*([jStrain;jw;3+(jNode-1)*DofPerNode]-1)+1];
    ind=[ind i1(1:3)*([jStrain;jw;4+(jNode-1)*DofPerNode]-1)+1];
    ind=[ind i1(1:3)*([jStrain;jw;5+(jNode-1)*DofPerNode]-1)+1];
    ind=[ind i1(1:3)*([jStrain;jw;6+(jNode-1)*DofPerNode]-1)+1];
           end
       end
   end

   out=c;out1=Nstrain;out2=ind;

%% #MitShellStrain ----------------------------------------------------------
elseif comstr(Cam,'mitshellstrain') 


  opt=varargin{carg};carg=carg+1;

  % 1. need to write the e at all integration points as c*opt.NDN(:)
   DofPerNode=length(opt.DofLabels); Nnode=opt.Nnode;
   Ndof=DofPerNode*Nnode; Nw=size(opt.N,1);
   Nshape=size(opt.N,2);

   [c,Nstrain,ind]=integrules('StrainNDN',opt,1,5:8);

  % 2. then create the linear combinations to interpolate the shear strain
  %    \gamma_r (r,s) = N_A(r,s) \gamma_r(A) + N_C  \gamma_r(C)
  %    \gamma_s (r,s) = N_D(r,s) \gamma_r(D) + N_B  \gamma_r(B)

  % A,B,C,D -> jw=5,6,7,8
  i1=[1 cumprod([Nstrain Nw Ndof])]; i2=[1 cumprod([Nnode Nw 3])];
  for j1=1:2
   if j1==1; jStrain=7; DOF=[3 5]; TyingW=[5 7];% gamma_r
   else      jStrain=8; DOF=[3 4]; TyingW=[6 8];% gamma_s
%   if j1==1; jStrain=7; DOF=[3 5]; TyingW=[6 8];% gamma_r
%   else      jStrain=8; DOF=[3 4]; TyingW=[5 7];% gamma_s
   end
  for jNode=1:Nnode;
      for jw=1:4;
   % Do not use the local value of gamma_r but combination of tying points
   c(i1(1:3)*([jStrain;jw;DOF(1)+(jNode-1)*DofPerNode]-1)+1,:)= ...
    opt.N(jw,TyingW(1))* ...
     c(i1(1:3)*([jStrain;TyingW(1);DOF(1)+(jNode-1)*DofPerNode]-1)+1,:)+ ...
    opt.N(jw,TyingW(2))* ...
     c(i1(1:3)*([jStrain;TyingW(2);DOF(1)+(jNode-1)*DofPerNode]-1)+1,:);
   c(i1(1:3)*([jStrain;jw;DOF(2)+(jNode-1)*DofPerNode]-1)+1,:)= ...
    opt.N(jw,TyingW(1))* ...
     c(i1(1:3)*([jStrain;TyingW(1);DOF(2)+(jNode-1)*DofPerNode]-1)+1,:)+ ...
    opt.N(jw,TyingW(2))* ...
     c(i1(1:3)*([jStrain;TyingW(2);DOF(2)+(jNode-1)*DofPerNode]-1)+1,:);
      end
  end 
  end % case gamma_r or gamma_s
  c(ind,:)=[];

  opt.NDN=[opt.N(:,1:Nnode)' opt.Nr(:,1:Nnode)' opt.Ns(:,1:Nnode)'];
  c=c*opt.NDN(:);
  c=reshape(c,Nw/2*Nstrain,Ndof);

  opt.MatrixIntegrationRule{1}=struct('c',c,'Nstrain',Nstrain, ...
    'Ndof',Ndof,'typ','rs');
  % eliminate tying point shape functions that are no longer needed
  opt.N=opt.N(:,1:Nnode); opt.Nr=opt.Nr(:,1:Nnode);opt.Ns=opt.Ns(:,1:Nnode);
  opt.Nnode=Nnode;opt.NDN=zeros(Nnode,3*Nw);

  % mass rule
  [c,Nstrain,ind]=integrules('StrainNDN',opt,2,5:8);c(ind,:)=[];
  opt.MatrixIntegrationRule{2}=struct('c',c,'Nstrain',Nstrain,'Ndof',Ndof);

  out=opt;

%% #Bar1, Line2 Integration rules for 1D [0,1] linear segment ----------------
elseif comstr(Cam,'bar1')||comstr(Cam,'line2') 

 xi = [0 0 0;1 0 0];
 if nargin<2; w='def'; else; w=varargin{carg};carg=carg+1;end
 [w,iN,Gdata]=QuadraturePoints('1d',w,xi,varargin{carg:end});

 r=w(:,1);N=[1-r r]; Nr=ones(size(r))*[-1 1];
 out=struct('w',w,'xi',xi(iN,:), ...
  'N',N(:,iN),'Nr',Nr(:,iN),'Nrr',zeros(size(r)), ...
  'NDN',zeros(2,2*size(w,1)),'jdet',zeros(size(w,1),1), ...
  'NDNLabels',{{'',',x'}},'type','bar1','Nnode',2,...
  'Gmesh',Gdata(1),'Gedge',Gdata(2));

%% #beam3 Integration rules for 1D [0,1] quadratic segment ------------------
elseif comstr(Cam,'beam3') 

 xi = [0 0 0;1 0 0;.5 0 0];
 iN=@(w)[1-3*w(:,1)+2*w(:,1).^2 2*w(:,1).^2-w(:,1) (w(:,1)-w(:,1).^2)*4 ];
 iNr=@(w)[4*w(:,1)-3 4*w(:,1)-1 4-8*w(:,1)];
 iNrr=@(w)ones(size(w,1),1)*[4 4 -8];
 
 if nargin<2; w='def'; else; w=varargin{carg};carg=carg+1;end
 if ischar(w)&&strcmpi(w,'inline')
  out=struct('N',iN,'Nr',iNr,'Nrr',iNrr,'xi',xi,'type','beam3');     
 else
  [w,jN,Gdata]=QuadraturePoints('1d',w,xi,varargin{carg:end});
  r=w(:,1);
  % r=[0;.5;1];inv([ones(size(r)) r r.^2])'
  N=iN(w);Nr=iNr(w);Nrr=iNrr(w);
  out=struct('w',w,'N',N(:,jN),'Nr',Nr(:,jN),'Nrr',Nrr(:,jN),'xi',xi(jN,:), ...
  'type','beam3','NDN',zeros(size(N,2),3*size(N,1)), ...
  'jdet',zeros(size(N,1),1),'Nnode',3,'NDNLabels',{{'',',x'}},...
  'Gmesh',Gdata(1),'Gedge',Gdata(2));
 end

%% #beam1 Integration rules for 1D [0,1] cubic segment -----------------------
elseif comstr(Cam,'beam1') 

 xi = [0 0 0;1 0 0];
 if nargin<2; w='def'; else; w=varargin{carg};carg=carg+1;end
 [w,iN,Gdata]=QuadraturePoints('1d',w,xi,varargin{carg:end});

 N=[1-3*w(:,1).^2+2*w(:,1).^3 w(:,1)-2*w(:,1).^2+w(:,1).^3 ...
      3*w(:,1).^2-2*w(:,1).^3   -w(:,1).^2+w(:,1).^3];

 Nr=[-6*w(:,1)+6*w(:,1).^2 1-4*w(:,1)+3*w(:,1).^2 ...
       6*w(:,1)-6*w(:,1).^2   -2*w(:,1)+3*w(:,1).^2];

 Nrr=[-6+12*w(:,1) -4+6*w(:,1) 6-12*w(:,1)   -2+6*w(:,1)];
 out=struct('w',w,'N',N,'Nr',Nr,'Nrr',Nrr,'xi',xi,'type','beam1', ...
   'Nw',size(w,1),'NDN',zeros(size(N,2),3*size(N,1)),'jdet',zeros(size(N,1),1), ...
   'Nnode',2,'NDNLabels',{{'',',x'}},...
  'Gmesh',Gdata(1),'Gedge',Gdata(2));
 %k=zeros(4);for j1=1:3 k=k+out.Nrr(j1,:)'*out.Nrr(j1,:)*w(j1,4);end

%% #quad4, #q4p Integration rules for 2D [-1,1 x -1,1] square ----------------
elseif comstr(Cam,'quad4') ||comstr(Cam,'q4p')

 xi = [-1 -1 0;1 -1 0;1 1 0;-1 1 0];

 iN =@(w)(1+w(:,1)*xi(:,1)').*(1+w(:,2)*xi(:,2)')/4;% ...
         %  1-w(:,1).^2 1-w(:,2).^2]
 iNr=@(w)xi(1:4,ones(1,size(w,1))*1)'.*(1+w(:,2)*xi(:,2)')/4;%...
        % -2*w(:,1) w(:,1)*0];
 iNs=@(w)xi(1:4,ones(1,size(w,1))*2)'.*(1+w(:,1)*xi(:,1)')/4;%...
        % w(:,1)*0 -2*w(:,2)];
 if nargin<2; w='def'; else;  w=varargin{carg};carg=carg+1; end
 if ischar(w)&&strcmpi(w,'inline')
  out=struct('xi',xi,'N',iN,'Nr',iNr,'Ns',iNs,'type','quad4');
 else
  [w,jN,Gdata]=QuadraturePoints('q2d',w,xi,varargin{carg:end});
  N=iN(w);Nr=iNr(w); Ns=iNs(w); 
  out=struct('N',N(:,jN),'Nr',Nr(:,jN),'Ns',Ns(:,jN), ...
      'Nw',size(w,1),'NDN',zeros(size(N,2),3*size(N,1)), ...
      'NDNLabels',{{'',',x',',y'}}, ...
      'jdet',zeros(size(N,1),1),'w',w,'Nnode',4,'xi',xi(jN,:),'type','quad4',...
      'Gmesh',Gdata(1),'Gedge',Gdata(2));
   if carg<=nargin&&isempty(strfind(Cam,'new')) % OBSOLETE
       out=integrules('2d',out,varargin{carg});carg=carg+1;
   end % OBSOLETE
 end
%% #q4q Integration rules for 2D [-1,1 x -1,1] square ------------------------
elseif comstr(Cam,'q4q')

 xi = [-1 -1 0;1 -1 0;1 1 0;-1 1 0];
 if nargin<2; w='def'; else;  w=varargin{carg};carg=carg+1; end
 [w,iN,Gdata]=QuadraturePoints('q2d',w,xi,varargin{carg:end});

 N = [(1+w(:,1)*xi(:,1)').*(1+w(:,2)*xi(:,2)')/4 ...
           1-w(:,1).^2 1-w(:,2).^2];
 Nr=[xi(1:4,ones(1,size(w,1))*1)'.*(1+w(:,2)*xi(:,2)')/4 ...
         -2*w(:,1) w(:,1)*0];
 Ns=[xi(1:4,ones(1,size(w,1))*2)'.*(1+w(:,1)*xi(:,1)')/4 ...
         w(:,1)*0 -2*w(:,2)];
 out=struct('N',N,'Nr',Nr,'Ns',Ns, ...
      'Nw',size(w,1),'NDN',zeros(size(N,2),3*size(N,1)), ...
      'NDNLabels',{{'',',x',',y'}}, ...
      'jdet',zeros(size(N,1),1),'w',w,'Nnode',4,'xi',xi,'type','q4q',...
      'Gmesh',Gdata(1),'Gedge',Gdata(2));

%% #quad9 --------------------------------------------------------------------
elseif comstr(Cam,'quad9')||comstr(Cam,'q9p')||comstr(Cam,'q9a')
% comstr(Cam,'q2')||
  xi = [-1 -1 0;1 -1 0;1 1 0;-1 1 0;0 -1 0;1 0 0;0 1 0;-1 0 0;0 0 0];

 if nargin<2; w='def'; else; w=varargin{carg};carg=carg+1;end 
 [w,iN,Gdata]=QuadraturePoints('q2d',w,xi);


 r=w(:,1); s=w(:,2);
 N=[1/4*r.*s-1/4*r.^2.*s+1/4*r.^2.*s.^2-1/4*r.*s.^2 ...
-1/4*r.*s-1/4*r.^2.*s+1/4*r.^2.*s.^2+1/4*r.*s.^2 ...
1/4*r.*s+1/4*r.*s.^2+1/4*r.^2.*s.^2+1/4*r.^2.*s ...
-1/4*r.*s+1/4*r.^2.*s+1/4*r.^2.*s.^2-1/4*r.*s.^2 ...
-1/2*s+1/2*s.^2+1/2*r.^2.*s-1/2*r.^2.*s.^2 ...
1/2*r-1/2*r.*s.^2+1/2*r.^2-1/2*r.^2.*s.^2 ...
1/2*s+1/2*s.^2-1/2*r.^2.*s-1/2*r.^2.*s.^2 ...
-1/2*r+1/2*r.*s.^2+1/2*r.^2-1/2*r.^2.*s.^2 ...
(-1+r.^2).*(-1+s.^2)];
%
Nr=[1/4*s-1/2*r.*s+1/2*r.*s.^2-1/4*s.^2 ...
-1/4*s-1/2*r.*s+1/2*r.*s.^2+1/4*s.^2 ...
 1/4*s+1/4*s.^2+1/2*r.*s.^2+1/2*r.*s ...
-1/4*s+1/2*r.*s+1/2*r.*s.^2-1/4*s.^2 ...
r.*s-r.*s.^2 ...
1/2-1/2*s.^2+r-r.*s.^2 ...
-r.*s-r.*s.^2 ... 
-1/2+1/2*s.^2+r-r.*s.^2 ...
2*r.*(-1+s.^2)];
%
Ns=[1/4*r-1/4*r.^2+1/2*r.^2.*s-1/2*r.*s ...
    -1/4*r-1/4*r.^2+1/2*r.^2.*s+1/2*r.*s ...
 1/4*r+1/2*r.*s+1/2*r.^2.*s+1/4*r.^2 ...
 -1/4*r+1/4*r.^2+1/2*r.^2.*s-1/2*r.*s ...
 -1/2+s+1/2*r.^2-r.^2.*s ...
 -r.*s-r.^2.*s ...
 1/2+s-1/2*r.^2-r.^2.*s ...
 r.*s-r.^2.*s ...
 2*(-1+r.^2).*s];

  out=struct('N',N,'Nr',Nr,'Ns',Ns,'Nw',size(w,1), ...
      'Nnode',9,'NDN',zeros(size(N,2),size(N,1)*3), ...
      'NDNLabels',{{'',',x',',y'}}, ...
      'jdet',zeros(size(w,1),1),'w',w, ...
      'DofLabels',{{'u','v','w','ru','rw'}},'type','quad9','xi',xi,...
      'Gmesh',Gdata(1),'Gedge',Gdata(2));

%% #quadb, #q8p --------------------------------------------------------------
elseif comstr(Cam,'quadb')||comstr(Cam,'q8p')

  xi = [-1 -1 0;1 -1 0;1 1 0;-1 1 0;0 -1 0;1 0 0;0 1 0;-1 0 0];

 if nargin<2; w='def'; else; w=varargin{carg};carg=carg+1;end
 [w,iN,Gdata]=QuadraturePoints('q2d',w,xi);

% r=w(:,1); s=w(:,2); xxx
% See Hugues 132-135                          na = 1/4 (1+xa x) (1 + ya y)
  N = [(1+w(:,1)*xi(1:4,1)').*(1+w(:,2)*xi(1:4,2)')/4 ...
      (1-w(:,1).^2).*(1-w(:,2))/2 ...   % n5 = 1/2 (1-x2  ) (1-y)
      (1-w(:,2).^2).*(1+w(:,1))/2 ...   % n6 = 1/2 (1+x)    (1-y2)
      (1-w(:,1).^2).*(1+w(:,2))/2 ...   % n7 = 1/2 (1-x2)   (1+y)
      (1-w(:,2).^2).*(1-w(:,1))/2 ];   % n8 = 1/2 (1-x)    (1-y2)

  Nr=[xi(1:4,ones(1,size(w,1))*1)'.*(1+w(:,2)*xi(1:4,2)')/4 ...
       ( -w(:,1)).*(1-w(:,2))       (1-w(:,2).^2)/2 ...
       ( -w(:,1)).*(1+w(:,2))      -(1-w(:,2).^2)/2];

  Ns=[xi(1:4,ones(1,size(w,1))*2)'.*(1+w(:,1)*xi(1:4,1)')/4 ...
       -(1-w(:,1).^2)/2     ( -w(:,2)).*(1+w(:,1)) ...
        (1-w(:,1).^2)/2     ( -w(:,2)).*(1-w(:,1))];

  % if n5 n8 present
  % n1 = n1 -1/2 (n5+n8)
  % n2 = n2 -1/2 (n5+n6)
  % n3 = n3 -1/2 (n6+n7)
  % n4 = n4 -1/2 (n7+n8)
  % [eye(4) -[1 0 0 1;1 1 0 0;0 1 1 0;0 0 1 1]/2;zeros(4) eye(4)];comstr(ans',-30)

  N (:,1:4)=N (:,1:4)-N (:,5:8)*[1 0 0 1;1 1 0 0;0 1 1 0;0 0 1 1]'/2;
  Nr(:,1:4)=Nr(:,1:4)-Nr(:,5:8)*[1 0 0 1;1 1 0 0;0 1 1 0;0 0 1 1]'/2;
  Ns(:,1:4)=Ns(:,1:4)-Ns(:,5:8)*[1 0 0 1;1 1 0 0;0 1 1 0;0 0 1 1]'/2;

  out=struct('N',N,'Nr',Nr,'Ns',Ns,'Nw',size(w,1), ...
      'Nnode',8,'NDN',zeros(size(N,2),size(N,1)*3), ...
      'NDNLabels',{{'',',x',',y'}}, ...
      'jdet',zeros(size(w,1),1),'w',w, ...
      'DofLabels',{{'u','v','w','ru','rw'}},'type','quadb','xi',xi,...
      'Gmesh',Gdata(1),'Gedge',Gdata(2));

%% #Tria3 Integration rules for 2D [0,0 1,0 0,1] triangle ---------------------------
elseif comstr(Cam,'tria3')||comstr(Cam,'t3p')

 % Shape functions are 1-r-s, r and s
 xi = [0 0 0;1 0 0;0 1 0];
 iN = @(w)[1-w(:,1)-w(:,2) w(:,1) w(:,2)];
 iNr=@(w)ones(size(w,1),1)*[-1 1 0];
 iNs=@(w)ones(size(w,1),1)*[-1 0 1];
 if nargin<2; w='def'; else; w=varargin{carg};carg=carg+1;end
 if ischar(w)&&strcmpi(w,'inline')
   out=struct('xi',xi,'N',iN,'Nr',iNr,'Ns',iNs,'type','tria3');
 else
  [w,jN,Gdata]=QuadraturePoints('t2d',w,xi);
  N=iN(w); Nr=iNr(w); Ns=iNs(w); nw=size(N,1);
  out=struct('N',N,'Nr',Nr,'Ns',Ns,'w',w,'Nw',size(w,1), ...
      'jdet',zeros(nw,1),'NDN',zeros(3,nw*3),'NDNLabels',{{'',',x',',y'}}, ...
      'Nnode',3,'xi',xi,'type','tria3',...
      'Gmesh',Gdata(1),'Gedge',Gdata(2));
  if carg<=nargin % OBSOLETE
      out=integrules('2d',out,varargin{carg});carg=carg+1;
  end % OBSOLETE
 end
 
%% #tria6 6 node triangle - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'tria6') ||comstr(Cam,'t6p')

 xi=[0 0 0;1 0 0;0 1 0;.5 0 0;.5 .5 0;0 .5 0];
 if nargin<2; w='def'; else; w=varargin{carg};carg=carg+1;end
 [w,jN,Gdata]=QuadraturePoints('t2d',w,xi);

 r=w(:,1)*2-1; s=w(:,2)*2-1;  z=zeros(size(r));

 % shape between [-1 1]
 N=[   .5*(r+s).*(r+s+1)  .5*(1+r).*r .5*(1+s).*s   ...
          -(1+r).*(r+s)      (1+r).*(1+s) -(1+s).*(r+s)   ];  

 Nr=[   .5*(2*r+2*s+1) .5*(1+2*r) z...
         -(2*r+s+1)     (1+s) -(1+s)  ]*2;

 Ns=[.5*(2*r+2*s+1) z .5*(1+2*s)  ...
          -(1+r)   (1+r) -(2*s+r+1) ]*2;

 out=struct('N',N,'Nr',Nr,'Ns',Ns, ...
      'Nw',size(w,1),'Nnode',6,'NDN',zeros(size(N,2),size(N,1)*3), ...
      'NDNLabels',{{'',',x',',y'}}, ...
      'jdet',zeros(size(w,1),1),'w',w,'DofLabels',{{'u','v','w'}}, ...
      'xi',xi,'type','tria6',...
      'Gmesh',Gdata(1),'Gedge',Gdata(2));

% ---------------------------------------------------------------------
elseif comstr(Cam,'tetra4')

 xi = [0 0 0;1 0 0;0 1 0;0 0 1];
 if nargin<2; w='def'; else; w=varargin{carg};carg=carg+1;end
 [w,iN,Gdata]=QuadraturePoints('t3d',w,xi);

  % Shape functions [1-r-s-t r s t] 
  r=w(:,1);  s=w(:,2);  t=w(:,3);  v=1-r-s-t;z=v*0;o=ones(size(z));

  N = [v r s t];
  Nr= [-o o z z];
  Ns= [-o z o z];
  Nt= [-o z z o];

 out=struct('N',N,'Nr',Nr,'Ns',Ns,'Nt',Nt, ...
      'Nw',size(w,1),'Nnode',4,'NDN',zeros(size(N,2),size(N,1)*4), ...
      'NDNLabels',{{'',',x',',y',',z'}}, ...
      'jdet',zeros(size(w,1),1),'w',w,'DofLabels',{{'u','v','w'}},...
      'xi',xi,'type','tetra4',...
      'Gmesh',Gdata(1),'Gedge',Gdata(2));

% ---------------------------------------------------------------------
elseif comstr(Cam,'tetra10')

  xi = [0 0 0; 1 0 0; 0 1 0; 0 0 1;
       .5 0 0;.5 .5 0; 0 .5 0; 0 0 .5; .5 0 .5;  0 .5 .5   ];
 if nargin<2; w='def'; else; w=varargin{carg};carg=carg+1;end
 [w,iN,Gdata]=QuadraturePoints('t3d',w,xi);

  % See also DOC code aster R30101b
  % v=1-r-s-t
  % [ v(2v-1) r(2r-1) s(2s-1) t(2t-1)  
  %   4rv  4rs  4sv  4tv  4rs  4st ]
  r=w(:,1);  s=w(:,2);  t=w(:,3);  v=1-r-s-t; z=zeros(size(s));

  N = [ v.*(2*v-1)  r.*(2*r-1)  s.*(2*s-1)  t.*(2*t-1)   ...
           4*r.*v      4*s.*r      4*s.*v     4*t.*v    4*r.*t    4*s.*t ];
  Nr = [ -4*v+1  4*r-1  z z   4*(v-r)  4*s -4*s     -4*t   4*t   z];
  Ns = [-4*v+1  z  4*s-1 z    -4*r     4*r  4*(v-s) -4*t    z   4*t];
  Nt = [-4*v+1  z z   4*t-1   -4*r     z   -4*s    4*(v-t) 4*r  4*s];
   


  out=struct('N',N,'Nr',Nr,'Ns',Ns,'Nt',Nt, ...
      'Nw',size(w,1),'Nnode',10,'NDN',zeros(size(N,2),size(N,1)*4), ...
      'NDNLabels',{{'',',x',',y',',z'}}, ...
      'jdet',zeros(size(w,1),1),'w',w,'DofLabels',{{'u','v','w'}}, ...
      'xi',xi,'type','tetra10',...
      'Gmesh',Gdata(1),'Gedge',Gdata(2));

% ---------------------------------------------------------------------
elseif comstr(Cam,'penta6')

  xi=[0 0 0; 1 0 0; 0 1 0; 
      0 0 1; 1 0 1; 0 1 1 ];
 if nargin<2; w='def'; else; w=varargin{carg};carg=carg+1;end
 [w,iN,Gdata]=QuadraturePoints('p3d',w,xi);

  r=w(:,1); s=w(:,2); t=w(:,3); z=zeros(size(t));
  N = [ (1-r-s).*(1-t)   r.*(1-t)    s.*(1-t)  ...
           (1-r-s).*t       r.*t        s.*t   ]  ;
  Nr = [ -(1-t)   (1-t)   z  -t         t      z   ]  ;
  Ns = [ -(1-t)   z    (1-t) -t       z    t   ]  ;
  Nt = [ -(1-r-s)   -r    -s  ...
           (1-r-s)     r     s   ]  ;

 out=struct('N',N,'Nr',Nr,'Ns',Ns,'Nt',Nt, ...
      'Nw',size(w,1),'Nnode',6,'NDN',zeros(size(N,2),size(N,1)*4), ...
      'NDNLabels',{{'',',x',',y',',z'}}, ...
      'jdet',zeros(size(w,1),1),'w',w,'DofLabels',{{'u','v','w'}}, ...
      'xi',xi,'type','penta6',...
      'Gmesh',Gdata(1),'Gedge',Gdata(2));

% ---------------------------------------------------------------------
elseif comstr(Cam,'penta15')

  xi = [0 0 0;1 0 0;0 1 0;
        0 0 1;1 0 1;0 1 1;
        .5 0 0;.5 .5 0; 0 .5 0;
        0 0 .5;1 0 .5; 0 1 .5;
        .5 0 1;.5 .5 1; 0 .5 1];

 if nargin<2; w='def'; else; w=varargin{carg};carg=carg+1;end
 [w,iN,Gdata]=QuadraturePoints('p3d',w,xi);

  % shape
  r=w(:,1); s=w(:,2); t=w(:,3); u=1-r-s; z=zeros(size(u));

  N = [ u.*(t-1).*(2*t-2*u+1) ...  % Bottom triangle
           r.*(t-1).*(2*t-2*r+1) ...
           s.*(t-1).*(2*t-2*s+1) ...
           u.*t.*(2*t+2*u-3) ... % Top triangle
           r.*t.*(2*t+2*r-3) ...
           s.*t.*(2*t+2*s-3) ... 
           4*u.*r.*(1-t) ... % bottom triangle mid nodes
           4*r.*s.*(1-t) ...
           4*s.*u.*(1-t) ... 
           4*u.*t.*(1-t) ... % mid-sides
           4*r.*t.*(1-t) ...
           4*s.*t.*(1-t) ...
           4*u.*r.*t ... % top triangle mid nodes
           4*r.*s.*t ...
           4*s.*u.*t ];  
  %
  Nr = [ (t-1).*(-2*t+4*u-1) ...  % Bottom triangle
           (t-1).*(2*t-4*r+1) ...
           z ...
           t.*(-2*t-4*u+3) ... % Top triangle
           t.*(2*t+4*r-3) ...
           z ... 
           4*(1-t).*(u-r) ... % bottom triangle mid nodes
            4*s.*(1-t) ...
           -4*s.*(1-t) ... 
           -4*t.*(1-t) ... % mid-sides
           4*t.*(1-t) ...
           z ...
            4*(-r+u).*t ... % top triangle mid nodes
            4*s.*t ...
           -4*s.*t ];  
  %
  Ns = [ (t-1).*(-2*t+4*u-1) ...  % Bottom triangle
           z ...
           (t-1).*(2*t-4*s+1) ...
           t.*(-2*t-4*u+3) ... % Top triangle
           z ...
           t.*(2*t+4*s-3) ... 
           4*-r.*(1-t) ... % bottom triangle mid nodes
           4* r.*(1-t) ...
           4*(u-s).*(1-t) ... 
           -4*t.*(1-t) ... % mid-sides
           z ...
            4*t.*(1-t) ...
           -4*r.*t ... % top triangle mid nodes
            4*r.*t ...
            4*(u-s).*t ];

  Nt = [u.*(4*t-2*u-1) ...  % Bottom triangle
           r.*(4*t-2*r-1) ...
           s.*(4*t-2*s-1) ...
           u.*(4*t+2*u-3) ... % Top triangle
           r.*(4*t+2*r-3) ...
           s.*(4*t+2*s-3) ... 
           -4*u.*r ... % bottom triangle mid nodes
           -4*r.*s ...
           -4*s.*u ... 
           4*u.*(1-2*t) ... % mid-sides
           4*r.*(1-2*t) ...
           4*s.*(1-2*t) ...
           4*u.*r ... % top triangle mid nodes
           4*r.*s ...
           4*s.*u ];  

  out=struct('N',N,'Nr',Nr,'Ns',Ns,'Nt',Nt, ...
      'Nw',size(w,1),'Nnode',15,'NDN',zeros(size(N,2),size(N,1)*4), ...
      'NDNLabels',{{'',',x',',y',',z'}}, ...
      'jdet',zeros(size(w,1),1),'w',w,'DofLabels',{{'u','v','w'}}, ...
      'xi',xi,'type','penta15',...
      'Gmesh',Gdata(1),'Gedge',Gdata(2));

%% #hexa8 ---------------------------------------------------------------------
elseif comstr(Cam,'flui8')||comstr(Cam,'hexa8')

 xi = [-1 -1 -1;1 -1 -1;1 1 -1;-1 1 -1;-1 -1 1;1 -1 1;1 1 1;-1 1 1];
 % xi = [ 1 1 1; -1 1 1; -1 -1 1; 1 -1 1; 1 1 -1; -1 1 -1; -1 -1 -1; 1 -1 -1];
 typ=0;
 if nargin<2; w='def'; else; 
   w=varargin{carg};carg=carg+1;
   if length(w)==1&&w(end)>1e4; typ=fix(w/1e4); w=rem(w,1e4);end
 end 
 [w,iN,Gdata]=QuadraturePoints('h3d',w,xi);
 if isstruct(w);out=w;return;end

 if typ==1 % alternate shape functions 10002 using bubble
  r=w(:,1);s=w(:,2);t=w(:,3);z=zeros(size(t));
  N = [(1+r*xi(:,1)').*(1+s*xi(:,2)').*(1+t*xi(:,3)')/8 ... % classical
      (1-r.^2) (1-s.^2) (1-t.^2) % quadratic shapes for shear protection
      ];
  Nr=[xi(1:8,ones(1,size(w,1))*1)'.*(1+s*xi(:,2)').* ...
        (1+t*xi(:,3)')/8 -2*r z z];
  Ns=[xi(1:8,ones(1,size(w,1))*2)'.*(1+r*xi(:,1)').* ...
        (1+t*xi(:,3)')/8 z -2*s z];
  Nt=[xi(1:8,ones(1,size(w,1))*3)'.*(1+r*xi(:,1)').* ...
        (1+s*xi(:,2)')/8 z z -2*t];
 else
  N = (1+w(:,1)*xi(:,1)').*(1+w(:,2)*xi(:,2)').*(1+w(:,3)*xi(:,3)')/8;
  Nr=xi(1:8,ones(1,size(w,1))*1)'.*(1+w(:,2)*xi(:,2)').* ...
        (1+w(:,3)*xi(:,3)')/8;
  Ns=xi(1:8,ones(1,size(w,1))*2)'.*(1+w(:,1)*xi(:,1)').* ...
        (1+w(:,3)*xi(:,3)')/8;
  Nt=xi(1:8,ones(1,size(w,1))*3)'.*(1+w(:,1)*xi(:,1)').* ...
        (1+w(:,2)*xi(:,2)')/8;
 end
 out=struct('N',N,'Nr',Nr,'Ns',Ns,'Nt',Nt,'Nw',size(w,1), ...
      'NDN',zeros(size(N,2),4*size(N,1)),'NDNLabels',{{'',',x',',y',',z'}}, ...
      'jdet',zeros(size(N,1),1),'w',w,'Nnode',8, ...
      'DofLabels',{{'u','v','w'}},'xi',xi,'type','hexa8',...
      'Gmesh',Gdata(1),'Gedge',Gdata(2));

 if carg<=nargin % Possibly compute the derivatives OBSOLETE
      out=integrules('3d',out,varargin{carg});carg=carg+1;
 end

%% #hexa20 ------------------------------------------------------------------
elseif comstr(Cam,'hexa20') %

  xi = [-1 -1 -1;1 -1 -1;1 1 -1;-1 1 -1;-1 -1 1;1 -1 1;1 1 1;-1 1 1;
        0 -1 -1;1 0 -1;0 1 -1;-1 0 -1;-1 -1 0;1 -1 0;1 1 0;-1 1 0;
        0 -1  1;1 0  1;0 1  1;-1 0  1];

 typ=0;
 if nargin<2; w='def'; else; 
   w=varargin{carg};carg=carg+1;
   if length(w)==1&&w(end)>1e4; typ=fix(w/1e4); w=rem(w,1e4);end
 end 
 [w,iN,Gdata]=QuadraturePoints('h3d',w,xi);
  

if 1==2
  % on [-1 1] ref. in Zienkiewicz vol.1 p133 and sdt4.1
  % corner nodes : (1+r)(1+s)(1+t)(r+s+t-2)/8
  % transformation : r=(r+1)/2; s=(s+1)/2; t=(t+1)/2; => on [0,1]
  r=w(:,1); s=w(:,2); t=w(:,3);N=zeros(1,8);
  for j1=1:8 % corner nodes
   ri=xi(j1,1);si=xi(j1,2); ti=xi(j1,3);
   N(:,j1)=(1+r*ri).*(1+s*si).*(1+t*ti).*(r*ri+s*si+t*ti-2)/8;
   %Nr(:,j1)=(1+s*si).*(1+t*ti)*(ri+2*r+s*si+t*ti-2)/8;
  end
  r2=r.^2; s2=s.^2; t2=t.^2;
  for j1=11:20 % midside nodes
   ri=xi(j1,1);si=xi(j1,2); ti=xi(j1,3);
   ri2=abs(ri)==0;si2=abs(si)==0;ti2=abs(ti)==0;
   N(:,j1)=(1+r*ri-ri2*r2).*(1+s*si-si2*s2).*(1+t*ti-ti2*t2)/4;
  end
end

  N =[(1+w(:,1)*xi(1:8,1)').*(1+w(:,2)*xi(1:8,2)').* ...
         (1+w(:,3)*xi(1:8,3)')/8 ... %na = 1/8 (1+xa x) (1 + ya y) (1 + za z)
   ... % n9  = 1/4 (1-x2 )(1-y)(1-z)
      (1-w(:,1).^2).*(1-w(:,2)).*(1-w(:,3))/4 ... 
   ... % n10 = 1/4 (1-y2 )(1+x)(1-z)
      (1-w(:,2).^2).*(1+w(:,1)).*(1-w(:,3))/4 ... 
   ... % n11 = 1/4 (1-x2 )(1+y)(1-z)
      (1-w(:,1).^2).*(1+w(:,2)).*(1-w(:,3))/4 ... 
   ... % n12 = 1/4 (1-y2 )(1-x)(1-z)
      (1-w(:,2).^2).*(1-w(:,1)).*(1-w(:,3))/4 ... 
   ... % n13 = 1/4 (1-z2 )(1-y)(1-x)
      (1-w(:,3).^2).*(1-w(:,2)).*(1-w(:,1))/4 ... 
   ... % n14 = 1/4 (1-z2 )(1+x)(1-y)
      (1-w(:,3).^2).*(1+w(:,1)).*(1-w(:,2))/4 ... 
   ... % n15 = 1/4 (1-z2 )(1+x)(1+u)
      (1-w(:,3).^2).*(1+w(:,1)).*(1+w(:,2))/4 ... 
   ... % n16 = 1/4 (1-z2 )(1-x)(1+y)
      (1-w(:,3).^2).*(1+w(:,2)).*(1-w(:,1))/4 ... 
   ... % n17 = 1/4 (1-x2 )(1-y)(1+z)
      (1-w(:,1).^2).*(1-w(:,2)).*(1+w(:,3))/4 ... 
   ... % n18 = 1/4 (1-y2 )(1+x)(1+z)
      (1-w(:,2).^2).*(1+w(:,1)).*(1+w(:,3))/4 ... 
   ... % n19 = 1/4 (1-x2 )(1+y)(1+z)
      (1-w(:,1).^2).*(1+w(:,2)).*(1+w(:,3))/4 ... 
   ... % n20 = 1/4 (1-y2 )(1-x)(1+z)
      (1-w(:,2).^2).*(1-w(:,1)).*(1+w(:,3))/4];

  Nr=[xi(1:8,ones(1,size(w,1))*1)'.*(1+w(:,2)*xi(1:8,2)').* ...
        (1+w(:,3)*xi(1:8,3)')/8 ...
      -w(:,1).*(1-w(:,2)).*(1-w(:,3))/2 ... 
      (1-w(:,2).^2).*(1-w(:,3))/4 ... 
      -w(:,1).*(1+w(:,2)).*(1-w(:,3))/2 ... 
      -(1-w(:,2).^2).*(1-w(:,3))/4 ... 
   ...
      -(1-w(:,3).^2).*(1-w(:,2))/4 ... 
       (1-w(:,3).^2).*(1-w(:,2))/4 ... 
       (1-w(:,3).^2).*(1+w(:,2))/4 ... 
      -(1-w(:,3).^2).*(1+w(:,2))/4 ... 
   ...
      -w(:,1).*(1-w(:,2)).*(1+w(:,3))/2 ... 
       (1-w(:,2).^2).*(1+w(:,3))/4 ... 
      -w(:,1).*(1+w(:,2)).*(1+w(:,3))/2 ... 
      -(1-w(:,2).^2).*(1+w(:,3))/4];

  Ns=[xi(1:8,ones(1,size(w,1))*2)'.*(1+w(:,1)*xi(1:8,1)').* ...
        (1+w(:,3)*xi(1:8,3)')/8 ...
      -(1-w(:,1).^2).*(1-w(:,3))/4 ... 
      -w(:,2).*(1+w(:,1)).*(1-w(:,3))/2 ... 
       (1-w(:,1).^2).*(1-w(:,3))/4 ... 
      -w(:,2).*(1-w(:,1)).*(1-w(:,3))/2 ... 
   ...
      -(1-w(:,3).^2).*(1-w(:,1))/4 ... 
      -(1-w(:,3).^2).*(1+w(:,1))/4 ... 
       (1-w(:,3).^2).*(1+w(:,1))/4 ... 
       (1-w(:,3).^2).*(1-w(:,1))/4 ... 
   ...
      -(1-w(:,1).^2).*(1+w(:,3))/4 ... 
      -w(:,2).*(1+w(:,1)).*(1+w(:,3))/2 ... 
       (1-w(:,1).^2).*(1+w(:,3))/4 ... 
      -w(:,2).*(1-w(:,1)).*(1+w(:,3))/2];
  Nt=[xi(1:8,ones(1,size(w,1))*3)'.*(1+w(:,1)*xi(1:8,1)').* ...
        (1+w(:,2)*xi(1:8,2)')/8 ...
      -(1-w(:,1).^2).*(1-w(:,2))/4 ... 
      -(1-w(:,2).^2).*(1+w(:,1))/4 ... 
      -(1-w(:,1).^2).*(1+w(:,2))/4 ... 
      -(1-w(:,2).^2).*(1-w(:,1))/4 ... 
   ...
      -w(:,3).*(1-w(:,2)).*(1-w(:,1))/2 ... 
      -w(:,3).*(1+w(:,1)).*(1-w(:,2))/2 ... 
      -w(:,3).*(1+w(:,1)).*(1+w(:,2))/2 ... 
      -w(:,3).*(1+w(:,2)).*(1-w(:,1))/2 ... 
   ...
      (1-w(:,1).^2).*(1-w(:,2))/4 ... 
      (1-w(:,2).^2).*(1+w(:,1))/4 ... 
      (1-w(:,1).^2).*(1+w(:,2))/4 ... 
      (1-w(:,2).^2).*(1-w(:,1))/4];
   i1 = [     9     9    10    11    13    14    15    16
             12    10    11    12    17    17    18    19
             13    14    15    16    20    18    19    20];
  for j1=1:8
    N (:,j1)=N (:,j1)-N (:,i1(:,j1))*[.5;.5;.5];  
    Nr(:,j1)=Nr(:,j1)-Nr(:,i1(:,j1))*[.5;.5;.5];  
    Ns(:,j1)=Ns(:,j1)-Ns(:,i1(:,j1))*[.5;.5;.5];  
    Nt(:,j1)=Nt(:,j1)-Nt(:,i1(:,j1))*[.5;.5;.5];  
  end

 out=struct('N',N,'Nr',Nr,'Ns',Ns,'Nt',Nt,'Nw',size(w,1), ...
      'NDN',zeros(size(N,2),size(N,1)*4), ...
      'NDNLabels',{{'',',x',',y',',z'}}, ...
      'jdet',zeros(size(w,1),1),'w',w,'DofLabels',{{'u','v','w'}}, ...
      'Nnode',20,'xi',xi,'type','hexa20',...
      'Gmesh',Gdata(1),'Gedge',Gdata(2));
 if typ>0 % Also add the hexa8  (u-p formulation of incompressible material)
    r2=integrules('hexa8',typ); % integrules('hexa20',30003)
    i2=size(out.N,1)+(1:size(r2.N,1));i3=1:size(r2.N,2);
    st={'N','Nr','Ns','Nt'};
    for j2=1:length(st);out.(st{j2})(i2,i3)=r2.(st{j2});end
    out.w(i2,:)=r2.w;
    out.NDN=zeros(20,4*size(out.N,1));out.jdet=zeros(size(out.w,1),1);
    out.Nw(2)=r2.Nw;
 end
 if carg<=nargin % Possibly compute the derivatives OBSOLETE
      out=integrules('3d',out,varargin{carg});carg=carg+1;
 end

%% #hexa21 ------------------------------------------------------------------
elseif comstr(Cam,'hexa21') %

 % xi gefdyn
   xi = [1 1 1;-1 1 1;-1 -1 1;1 -1 1;1 1 -1;-1 1 -1;-1 -1 -1;1 -1 -1;
        0 1 1;-1 0 1;0 -1 1;1 0 1;0 1 -1;-1 0 -1;0 -1 -1;1 0 -1;
        1 1 0; -1 1 0; -1 -1 0;1 -1 0; 0 0 0];
 if nargin<2; w='def'; else; w=varargin{carg};carg=carg+1;end
 [w,iN,Gdata]=QuadraturePoints('h3d',w,xi);
 if ischar(w); w=lower(w); end

  N =[(1+w(:,1)*xi(1:8,1)').*(1+w(:,2)*xi(1:8,2)').* ...
         (1+w(:,3)*xi(1:8,3)')/8 ... %na = 1/8 (1+xa x) (1 + ya y) (1 + za z)
   ... % n9  = 1/4 (1-x2 )(1+y)(1+z)
      (1-w(:,1).^2).*(1+w(:,2)).*(1+w(:,3))/4 ... 
   ... % n10 = 1/4 (1-y2 )(1-x)(1+z)
      (1-w(:,2).^2).*(1-w(:,1)).*(1+w(:,3))/4 ... 
   ... % n11 = 1/4 (1-x2 )(1-y)(1+z)
      (1-w(:,1).^2).*(1-w(:,2)).*(1+w(:,3))/4 ... 
   ... % n12 = 1/4 (1-y2 )(1+x)(1+z)
      (1-w(:,2).^2).*(1+w(:,1)).*(1+w(:,3))/4 ... 
   ... % n13 = 1/4 (1-x2 )(1+y)(1-z)
      (1-w(:,1).^2).*(1+w(:,2)).*(1-w(:,3))/4 ... 
   ... % n14 = 1/4 (1-y2 )(1-x)(1-z)
      (1-w(:,2).^2).*(1-w(:,1)).*(1-w(:,3))/4 ... 
   ... % n15 = 1/4 (1-x2 )(1-y)(1-z)
      (1-w(:,1).^2).*(1-w(:,2)).*(1-w(:,3))/4 ...  % n16 = 1/4 (1-y2 )(1+x)(1-z)
      (1-w(:,2).^2).*(1+w(:,1)).*(1-w(:,3))/4 ...  % n17 = 1/4 (1-z2 )(1+x)(1+y)
      (1-w(:,3).^2).*(1+w(:,1)).*(1+w(:,2))/4 ... % n18 = 1/4 (1-z2 )(1-x)(1+y)
      (1-w(:,3).^2).*(1-w(:,1)).*(1+w(:,2))/4 ... % n19 = 1/4 (1-z2 )(1-x)(1-y)
      (1-w(:,3).^2).*(1-w(:,1)).*(1-w(:,2))/4 ... % n20 = 1/4 (1-z2 )(1+x)(1-y)
      (1-w(:,3).^2).*(1+w(:,1)).*(1-w(:,2))/4 ... % n21 = 1/4 (1-z2 )(1-x2)(1-y2)
      (1-w(:,3).^2).*(1-w(:,1).^2).*(1-w(:,2).^2)];

Nr=[xi(1:8,ones(1,size(w,1))*1)'.*(1+w(:,2)*xi(1:8,2)').* ...
        (1+w(:,3)*xi(1:8,3)')/8 ...
   ... % n9  = 1/4 (1-x2 )(1+y)(1+z)
      -w(:,1).*(1+w(:,2)).*(1+w(:,3))/2 ... 
   ... % n10 = 1/4 (1-y2 )(1-x)(1+z)
      -(1-w(:,2).^2).*(1+w(:,3))/4 ... 
   ... % n11 = 1/4 (1-x2 )(1-y)(1+z)
      -w(:,1).*(1-w(:,2)).*(1+w(:,3))/2 ... 
   ... % n12 = 1/4 (1-y2 )(1+x)(1+z)
      (1-w(:,2).^2).*(1+w(:,3))/4 ... 
   ... % n13 = 1/4 (1-x2 )(1+y)(1-z)
      -w(:,1).*(1+w(:,2)).*(1-w(:,3))/2 ... 
   ... % n14 = 1/4 (1-y2 )(1-x)(1-z)
      -(1-w(:,2).^2).*(1-w(:,3))/4 ... 
   ... % n15 = 1/4 (1-x2 )(1-y)(1-z)
      -w(:,1).*(1-w(:,2)).*(1-w(:,3))/2 ... 
   ... % n16 = 1/4 (1-y2 )(1+x)(1-z)
      (1-w(:,2).^2).*(1-w(:,3))/4 ... 
   ... % n17 = 1/4 (1-z2 )(1+x)(1+y)
      (1-w(:,3).^2).*(1+w(:,2))/4 ... 
   ... % n18 = 1/4 (1-z2 )(1-x)(1+y)
      -(1-w(:,3).^2).*(1+w(:,2))/4 ... 
   ... % n19 = 1/4 (1-z2 )(1-x)(1-y)
      -(1-w(:,3).^2).*(1-w(:,2))/4 ... 
   ... % n20 = 1/4 (1-z2 )(1+x)(1-y)
      (1-w(:,3).^2).*(1-w(:,2))/4 ...
... % n21 = 1/4 (1-z2 )(1-x2)(1-y2)
      -2*(1-w(:,3).^2).*w(:,1).*(1-w(:,2).^2)];
    


  Ns=[xi(1:8,ones(1,size(w,1))*2)'.*(1+w(:,1)*xi(1:8,1)').* ...
        (1+w(:,3)*xi(1:8,3)')/8 ...
   ... % n9  = 1/4 (1-x2 )(1+y)(1+z)
      (1-w(:,1).^2).*(1+w(:,3))/4 ... 
   ... % n10 = 1/4 (1-y2 )(1-x)(1+z)
      -w(:,2).*(1-w(:,1)).*(1+w(:,3))/2 ... 
   ... % n11 = 1/4 (1-x2 )(1-y)(1+z)
      -(1-w(:,1).^2).*(1+w(:,3))/4 ... 
   ... % n12 = 1/4 (1-y2 )(1+x)(1+z)
      -w(:,2).*(1+w(:,1)).*(1+w(:,3))/2 ... 
   ... % n13 = 1/4 (1-x2 )(1+y)(1-z)
      (1-w(:,1).^2).*(1-w(:,3))/4 ... 
   ... % n14 = 1/4 (1-y2 )(1-x)(1-z)
      -w(:,2).*(1-w(:,1)).*(1-w(:,3))/2 ... 
   ... % n15 = 1/4 (1-x2 )(1-y)(1-z)
      -(1-w(:,1).^2).*(1-w(:,3))/4 ... 
   ... % n16 = 1/4 (1-y2 )(1+x)(1-z)
      -w(:,2).*(1+w(:,1)).*(1-w(:,3))/2 ... 
   ... % n17 = 1/4 (1-z2 )(1+x)(1+y)
      (1-w(:,3).^2).*(1+w(:,1))/4 ... 
   ... % n18 = 1/4 (1-z2 )(1-x)(1+y)
      (1-w(:,3).^2).*(1-w(:,1))/4 ... 
   ... % n19 = 1/4 (1-z2 )(1-x)(1-y)
      -(1-w(:,3).^2).*(1-w(:,1))/4 ... 
   ... % n20 = 1/4 (1-z2 )(1+x)(1-y)
      -(1-w(:,3).^2).*(1+w(:,1))/4 ...
... % n21 = 1/4 (1-z2 )(1-x2)(1-y2)
      -2*(1-w(:,3).^2).*(1-w(:,1).^2).*w(:,2)];



  Nt=[xi(1:8,ones(1,size(w,1))*3)'.*(1+w(:,1)*xi(1:8,1)').* ...
        (1+w(:,2)*xi(1:8,2)')/8 ...

    ... % n9  = 1/4 (1-x2 )(1+y)(1+z)
      (1-w(:,1).^2).*(1+w(:,2))/4 ... 
   ... % n10 = 1/4 (1-y2 )(1-x)(1+z)
      (1-w(:,2).^2).*(1-w(:,1))/4 ... 
   ... % n11 = 1/4 (1-x2 )(1-y)(1+z)
      (1-w(:,1).^2).*(1-w(:,2))/4 ... 
   ... % n12 = 1/4 (1-y2 )(1+x)(1+z)
      (1-w(:,2).^2).*(1+w(:,1))/4 ... 
   ... % n13 = 1/4 (1-x2 )(1+y)(1-z)
      -(1-w(:,1).^2).*(1+w(:,2))/4 ... 
   ... % n14 = 1/4 (1-y2 )(1-x)(1-z)
      -(1-w(:,2).^2).*(1-w(:,1))/4 ... 
   ... % n15 = 1/4 (1-x2 )(1-y)(1-z)
      -(1-w(:,1).^2).*(1-w(:,2))/4 ... 
   ... % n16 = 1/4 (1-y2 )(1+x)(1-z)
      -(1-w(:,2).^2).*(1+w(:,1))/4 ... 
   ... % n17 = 1/4 (1-z2 )(1+x)(1+y)
      -w(:,3).*(1+w(:,1)).*(1+w(:,2))/2 ... 
   ... % n18 = 1/4 (1-z2 )(1-x)(1+y)
      -w(:,3).*(1-w(:,1)).*(1+w(:,2))/2 ... 
   ... % n19 = 1/4 (1-z2 )(1-x)(1-y)
      -w(:,3).*(1-w(:,1)).*(1-w(:,2))/2 ... 
   ... % n20 = 1/4 (1-z2 )(1+x)(1-y)
      -w(:,3).*(1+w(:,1)).*(1-w(:,2))/2 ...
... % n21 = 1/4 (1-z2 )(1-x2)(1-y2)
      -2*w(:,3).*(1-w(:,1).^2).*(1-w(:,2).^2)];


  
  i1 = [     9     9    11    11    13   13  15    15    
             12    10   10    12    16   14  14    16    
             17    18   19    20    17   18  19    20    ];
  for j1=1:8
    N (:,j1)=N (:,j1)-N (:,i1(:,j1))*[.5;.5;.5];  
    Nr(:,j1)=Nr(:,j1)-Nr(:,i1(:,j1))*[.5;.5;.5];  
    Ns(:,j1)=Ns(:,j1)-Ns(:,i1(:,j1))*[.5;.5;.5];  
    Nt(:,j1)=Nt(:,j1)-Nt(:,i1(:,j1))*[.5;.5;.5];  
  end


%
% si le noeuds 21 existe
%

 for j1=1:8
    N (:,j1)=N (:,j1)-N (:,21)*.125;  
    Nr (:,j1)=Nr (:,j1)-Nr (:,21)*.125;  
    Ns (:,j1)=Ns (:,j1)-Ns (:,21)*.125;  
    Nt (:,j1)=Nt (:,j1)-Nt (:,21)*.125;   
end
 for j1=9:20
    N (:,j1)=N (:,j1)-N (:,21)*.25;  
    Nr (:,j1)=Nr (:,j1)-Nr (:,21)*.25;  
    Ns (:,j1)=Ns (:,j1)-Ns (:,21)*.25;  
    Nt (:,j1)=Nt (:,j1)-Nt (:,21)*.25;   
end

 out=struct('N',N,'Nr',Nr,'Ns',Ns,'Nt',Nt, ...
      'NDN',zeros(size(N,2),size(N,1)*4),'NDNLabels',{{'',',x',',y',',z'}}, ...
      'jdet',zeros(size(FEw,1),1),'w',FEw,'DofLabels',{{'u','v','w'}}, ...
      'Nnode',20,'Nw',size(FEw,1),'type','hexa20',...
      'Gmesh',Gdata(1),'Gedge',Gdata(2));

%% #hexa27 -------------------------------------------------------------------
elseif comstr(Cam,'hexa27')

 xi=[-1,1,1,-1, -1,1,1,-1, 0,1,0,-1, -1,1,1,-1,0,1,0,-1,0.,-1.,0.,0.,1.,0.,0.;
    -1, -1,1,1,  -1,-1,1,1, -1,0,1,0,-1,-1,1,1,-1,0,1,0,0.,0.,-1.,0.,0.,1.,0.;
   -1,-1,-1,-1, 1,1,1,1,  -1,-1,-1,-1, 0,0,0,0, 1,1,1,1,-1.,0.,0.,1.,0.,0.,0]';
 if nargin<2; w='def'; else; w=varargin{carg};carg=carg+1;end
 [w,iN,Gdata]=QuadraturePoints('h3d',w,xi);
 if isfield(w,'N');out=w;return;end

l=zeros(27,size(w,1));lr=zeros(27,size(w,1));ls=zeros(27,size(w,1));lz=ls;
for jw=1:size(w,1); r=w(jw,1); s=w(jw,2);z=w(jw,3);

% Fonctions de base 
l(1,jw)=-1/8*(1-r)*r*(1-s)*s*(1-z)*z;
l(2,jw)=1/8*(1+r)*r*(1-s)*s*(1-z)*z;
l(3,jw)=-1/8*(1+r)*r*(1+s)*s*(1-z)*z;
l(4,jw)=1/8*(1-r)*r*(1+s)*s*(1-z)*z;
l(5,jw)=1/8*(1-r)*r*(1-s)*s*(1+z)*z;
l(6,jw)=-1/8*(1+r)*r*(1-s)*s*(1+z)*z;
l(7,jw)=1/8*(1+r)*r*(1+s)*s*(1+z)*z;
l(8,jw)=-1/8*(1-r)*r*(1+s)*s*(1+z)*z;
l(9,jw)=.25*(1-r^2)*(1-s)*s*(1-z)*z;
l(10,jw)=-.25*(1+r)*r*(1-s^2)*(1-z)*z;
l(11,jw)=-.25*(1-r^2)*(1+s)*s*(1-z)*z;
l(12,jw)=.25*(1-r)*r*(1-s^2)*(1-z)*z;
l(13,jw)=.25*(1-r)*r*(1-s)*s*(1-z^2);
l(14,jw)=-.25*(1+r)*r*(1-s)*s*(1-z^2);
l(15,jw)=.25*(1+r)*r*(1+s)*s*(1-z^2);
l(16,jw)=-.25*(1-r)*r*(1+s)*s*(1-z^2);
l(17,jw)=-.25*(1-r^2)*(1-s)*s*(1+z)*z;
l(18,jw)=.25*(1+r)*r*(1-s^2)*(1+z)*z;
l(19,jw)=.25*(1-r^2)*(1+s)*s*(1+z)*z;
l(20,jw)=-.25*(1-r)*r*(1-s^2)*(1+z)*z;
l(21,jw)=-.5*(1-r^2)*(1-s^2)*z*(1-1.*z);
l(22,jw)=-.5*(1-1.*r)*r*(1-s^2)*(1-z^2);
l(23,jw)=-.5*(1-r^2)*s*(1-1.*s)*(1-z^2);
l(24,jw)=.5*(1-r^2)*(1-s^2)*z*(1+1.*z);
l(25,jw)=.5*(1+1.*r)*r*(1-s^2)*(1-z^2);
l(26,jw)=.5*(1-r^2)*s*(1+1.*s)*(1-z^2);
l(27,jw)=(1-r^2)*(1-s^2)*(1-z^2);

% Derivees en r
lr(1,jw)=1/8*r*(1-s)*s*(1-z)*z-1/8*(1-r)*(1-s)*s*(1-z)*z;
lr(2,jw)=1/8*r*(1-s)*s*(1-z)*z+1/8*(1+r)*(1-s)*s*(1-z)*z;
lr(3,jw)=-1/8*r*(1+s)*s*(1-z)*z-1/8*(1+r)*(1+s)*s*(1-z)*z;
lr(4,jw)=-1/8*r*(1+s)*s*(1-z)*z+1/8*(1-r)*(1+s)*s*(1-z)*z;
lr(5,jw)=-1/8*r*(1-s)*s*(1+z)*z+1/8*(1-r)*(1-s)*s*(1+z)*z;
lr(6,jw)=-1/8*r*(1-s)*s*(1+z)*z-1/8*(1+r)*(1-s)*s*(1+z)*z;
lr(7,jw)=1/8*r*(1+s)*s*(1+z)*z+1/8*(1+r)*(1+s)*s*(1+z)*z;
lr(8,jw)=1/8*r*(1+s)*s*(1+z)*z-1/8*(1-r)*(1+s)*s*(1+z)*z;
lr(9,jw)=-.50*r*(1-s)*s*(1-z)*z;
lr(10,jw)=-.25*r*(1-s^2)*(1-z)*z-.25*(1+r)*(1-s^2)*(1-z)*z;
lr(11,jw)=.50*r*(1+s)*s*(1-z)*z;
lr(12,jw)=-.25*r*(1-s^2)*(1-z)*z+.25*(1-r)*(1-s^2)*(1-z)*z;
lr(13,jw)=-.25*r*(1-s)*s*(1-z^2)+.25*(1-r)*(1-s)*s*(1-z^2);
lr(14,jw)=-.25*r*(1-s)*s*(1-z^2)-.25*(1+r)*(1-s)*s*(1-z^2);
lr(15,jw)=.25*r*(1+s)*s*(1-z^2)+.25*(1+r)*(1+s)*s*(1-z^2);
lr(16,jw)=.25*r*(1+s)*s*(1-z^2)-.25*(1-r)*(1+s)*s*(1-z^2);
lr(17,jw)=.50*r*(1-s)*s*(1+z)*z;
lr(18,jw)=.25*r*(1-s^2)*(1+z)*z+.25*(1+r)*(1-s^2)*(1+z)*z;
lr(19,jw)=-.50*r*(1+s)*s*(1+z)*z;
lr(20,jw)=.25*r*(1-s^2)*(1+z)*z-.25*(1-r)*(1-s^2)*(1+z)*z;
lr(21,jw)=1.0*r*(1-s^2)*z*(1-1.*z);
lr(22,jw)=.5*r*(1-s^2)*(1-z^2)-.5*(1-1.*r)*(1-s^2)*(1-z^2);
lr(23,jw)=1.0*r*s*(1-1.*s)*(1-z^2);
lr(24,jw)=-1.0*r*(1-s^2)*z*(1+1.*z);
lr(25,jw)=.5*r*(1-s^2)*(1-z^2)+.5*(1+1.*r)*(1-s^2)*(1-z^2);
lr(26,jw)=-1.0*r*s*(1+1.*s)*(1-z^2);
lr(27,jw)=-2*r*(1-s^2)*(1-z^2);

% Derivees en s
ls(1,jw)=1/8*(1-r)*r*s*(1-z)*z-1/8*(1-r)*r*(1-s)*(1-z)*z;
ls(2,jw)=-1/8*(1+r)*r*s*(1-z)*z+1/8*(1+r)*r*(1-s)*(1-z)*z;
ls(3,jw)=-1/8*(1+r)*r*s*(1-z)*z-1/8*(1+r)*r*(1+s)*(1-z)*z;
ls(4,jw)=1/8*(1-r)*r*s*(1-z)*z+1/8*(1-r)*r*(1+s)*(1-z)*z;
ls(5,jw)=-1/8*(1-r)*r*s*(1+z)*z+1/8*(1-r)*r*(1-s)*(1+z)*z;
ls(6,jw)=1/8*(1+r)*r*s*(1+z)*z-1/8*(1+r)*r*(1-s)*(1+z)*z;
ls(7,jw)=1/8*(1+r)*r*s*(1+z)*z+1/8*(1+r)*r*(1+s)*(1+z)*z;
ls(8,jw)=-1/8*(1-r)*r*s*(1+z)*z-1/8*(1-r)*r*(1+s)*(1+z)*z;
ls(9,jw)=-.25*(1-r^2)*s*(1-z)*z+.25*(1-r^2)*(1-s)*(1-z)*z;
ls(10,jw)=.50*(1+r)*r*s*(1-z)*z;
ls(11,jw)=-.25*(1-r^2)*s*(1-z)*z-.25*(1-r^2)*(1+s)*(1-z)*z;
ls(12,jw)=-.50*(1-r)*r*s*(1-z)*z;
ls(13,jw)=-.25*(1-r)*r*s*(1-z^2)+.25*(1-r)*r*(1-s)*(1-z^2);
ls(14,jw)=.25*(1+r)*r*s*(1-z^2)-.25*(1+r)*r*(1-s)*(1-z^2);
ls(15,jw)=.25*(1+r)*r*s*(1-z^2)+.25*(1+r)*r*(1+s)*(1-z^2);
ls(16,jw)=-.25*(1-r)*r*s*(1-z^2)-.25*(1-r)*r*(1+s)*(1-z^2);
ls(17,jw)=.25*(1-r^2)*s*(1+z)*z-.25*(1-r^2)*(1-s)*(1+z)*z;
ls(18,jw)=-.50*(1+r)*r*s*(1+z)*z;
ls(19,jw)=.25*(1-r^2)*s*(1+z)*z+.25*(1-r^2)*(1+s)*(1+z)*z;
ls(20,jw)=.50*(1-r)*r*s*(1+z)*z;
ls(21,jw)=1.0*(1-r^2)*s*z*(1-1.*z);
ls(22,jw)=1.0*(1-1.*r)*r*s*(1-z^2);
ls(23,jw)=-.5*(1-r^2)*(1-1.*s)*(1-z^2)+.5*(1-r^2)*s*(1-z^2);
ls(24,jw)=-1.0*(1-r^2)*s*z*(1+1.*z);
ls(25,jw)=-1.0*(1+1.*r)*r*s*(1-z^2);
ls(26,jw)=.5*(1-r^2)*(1+1.*s)*(1-z^2)+.5*(1-r^2)*s*(1-z^2);
ls(27,jw)=-2*(1-r^2)*s*(1-z^2);

% Derivees en z
lz(1,jw)=1/8*(1-r)*r*(1-s)*s*z-1/8*(1-r)*r*(1-s)*s*(1-z);
lz(2,jw)=-1/8*(1+r)*r*(1-s)*s*z+1/8*(1+r)*r*(1-s)*s*(1-z);
lz(3,jw)=1/8*(1+r)*r*(1+s)*s*z-1/8*(1+r)*r*(1+s)*s*(1-z);
lz(4,jw)=-1/8*(1-r)*r*(1+s)*s*z+1/8*(1-r)*r*(1+s)*s*(1-z);
lz(5,jw)=1/8*(1-r)*r*(1-s)*s*z+1/8*(1-r)*r*(1-s)*s*(1+z);
lz(6,jw)=-1/8*(1+r)*r*(1-s)*s*z-1/8*(1+r)*r*(1-s)*s*(1+z);
lz(7,jw)=1/8*(1+r)*r*(1+s)*s*z+1/8*(1+r)*r*(1+s)*s*(1+z);
lz(8,jw)=-1/8*(1-r)*r*(1+s)*s*z-1/8*(1-r)*r*(1+s)*s*(1+z);
lz(9,jw)=-.25*(1-r^2)*(1-s)*s*z+.25*(1-r^2)*(1-s)*s*(1-z);
lz(10,jw)=.25*(1+r)*r*(1-s^2)*z-.25*(1+r)*r*(1-s^2)*(1-z);
lz(11,jw)=.25*(1-r^2)*(1+s)*s*z-.25*(1-r^2)*(1+s)*s*(1-z);
lz(12,jw)=-.25*(1-r)*r*(1-s^2)*z+.25*(1-r)*r*(1-s^2)*(1-z);
lz(13,jw)=-.50*(1-r)*r*(1-s)*s*z;
lz(14,jw)=.50*(1+r)*r*(1-s)*s*z;
lz(15,jw)=-.50*(1+r)*r*(1+s)*s*z;
lz(16,jw)=.50*(1-r)*r*(1+s)*s*z;
lz(17,jw)=-.25*(1-r^2)*(1-s)*s*z-.25*(1-r^2)*(1-s)*s*(1+z);
lz(18,jw)=.25*(1+r)*r*(1-s^2)*z+.25*(1+r)*r*(1-s^2)*(1+z);
lz(19,jw)=.25*(1-r^2)*(1+s)*s*z+.25*(1-r^2)*(1+s)*s*(1+z);
lz(20,jw)=-.25*(1-r)*r*(1-s^2)*z-.25*(1-r)*r*(1-s^2)*(1+z);
lz(21,jw)=-.5*(1-r^2)*(1-s^2)*(1-1.*z)+.5*(1-r^2)*(1-s^2)*z;
lz(22,jw)=1.0*(1-1.*r)*r*(1-s^2)*z;
lz(23,jw)=1.0*(1-r^2)*s*(1-1.*s)*z;
lz(24,jw)=.5*(1-r^2)*(1-s^2)*(1+1.*z)+.5*(1-r^2)*(1-s^2)*z;
lz(25,jw)=-1.0*(1+1.*r)*r*(1-s^2)*z;
lz(26,jw)=-1.0*(1-r^2)*s*(1+1.*s)*z;
lz(27,jw)=-2*(1-r^2)*(1-s^2)*z;

end %jw

out=struct('N',l','Nr',lr','Ns',ls','Nt',lz', ...
  'NDN',zeros(27,4*size(w,1)), ...
      'NDNLabels',{{'',',x',',y',',z'}}, ...
      'jdet',zeros(size(w,1),1),'w',w,'DofLabels',{{'u','v','w'}}, ...
      'Nnode',27,'Nw',size(w,1),'xi',xi,'type','hexa27',...
      'Gmesh',[],'Gedge',[]);

%% #pyra5 ---------------------------------------------------------------------
elseif comstr(Cam,'pyra5')

 xi = [1 0 0;0 1 0;-1 0 0;0 -1 0;0 0 1];
 typ=0;
 if nargin<2; w='def'; 
 else; 
   w=varargin{carg};carg=carg+1;
   if length(w)==1&&w(end)>1e4; typ=fix(w/1e4); w=rem(w,1e4);end
 end 
 [w,iN,Gdata]=QuadraturePoints('y3d',w,xi);
 r=w(:,1); s=w(:,2);t=w(:,3);
 N = [(-r+s+t-1).*(-r-s+t-1)./(1-t)/4 ...
      (-r-s+t-1).*(r-s+t-1)./(1-t)/4 ...
      ( r+s+t-1).*( r-s+t-1)./(1-t)/4 ...
      ( r+s+t-1).*(-r+s+t-1)./(1-t)/4 ...
      t];
 Nr=[(1+r-t)./(1-t)./2 -r./(1-t)./2 (r+t-1)./(1-t)./2 -r./(1-t)./2 zeros(size(w,1),1)];
 Ns=[-s./(1-t)./2 (1+s-t)./(1-t)./2 -s./(1-t)./2 (s+t-1)./(1-t)./2 zeros(size(w,1),1)];
 Nt=[(r.^2-s.^2-t.^2+2.*t-1)./(1-t).^2./4  -(r.^2-s.^2+t.^2-2.*t+1)./(1-t).^2./4 ...
   (r.^2-s.^2-t.^2+2.*t-1)./(1-t).^2./4  -(r.^2-s.^2+t.^2-2.*t+1)./(1-t).^2./4  ones(size(w,1),1)];
 out=struct('N',N,'Nr',Nr,'Ns',Ns,'Nt',Nt,'Nw',size(w,1), ...
      'NDN',zeros(size(N,2),4*size(N,1)),'NDNLabels',{{'',',x',',y',',z'}}, ...
      'jdet',zeros(size(N,1),1),'w',w,'Nnode',5, ...
      'DofLabels',{{'u','v','w'}},'xi',xi,'type','pyra5',...
      'Gmesh',Gdata(1),'Gedge',Gdata(2));

 if carg<=nargin % Possibly compute the derivatives OBSOLETE
      out=integrules('3d',out,varargin{carg});carg=carg+1;
 end
%% #pyra13 ---------------------------------------------------------------------
elseif comstr(Cam,'pyra13')
%EC=integrules('pyra5');i1=pyra5('edges');[EC.xi;(EC.xi(i1(:,1),:)+EC.xi(i1(:,2),:))/2];comstr(ans,-30)

 xi = [1,0,0;0,1,0;-1,0,0;0,-1,0;0,0,1;0.5,0.5,0;-0.5,0.5,0;-0.5,-0.5,0;
     0.5,-0.5,0;0.5,0,0.5;0,0.5,0.5;-0.5,0,0.5;0,-0.5,0.5];
 typ=0;
 if nargin<2; w='def'; 
 else; 
   w=varargin{carg};carg=carg+1;
   if length(w)==1&&w(end)>1e4; typ=fix(w/1e4); w=rem(w,1e4);end
 end 
 [w,iN,Gdata]=QuadraturePoints('y3d',w,xi); 
 r=w(:,1); s=w(:,2);t=w(:,3);
 
 N = [ (-r+s+t-1.).*(-r-s+t-1.).*(r-0.5) ./ (2.*(1.-t)) ...
   (-r-s+t-1.).*(r-s+t-1.).*(s-0.5) ./ (2.*(1.-t)) ...
   (r-s+t-1.).*(r+s+t-1.).*(-r-0.5) ./ (2.*(1.-t)) ...
   (r+s+t-1.).*(-r+s+t-1.).*(-s-0.5) ./ (2.*(1.-t)) ...
   2.*t.*(t-0.5) ...
   -(-r+s+t-1.).*(-r-s+t-1.).*(r-s+t-1.) ./ (2.*(1.-t)) ...
   -(-r-s+t-1.).*(r-s+t-1.).*(r+s+t-1.) ./ (2.*(1.-t)) ...
   -(r-s+t-1.).*(r+s+t-1.).*(-r+s+t-1.) ./ (2.*(1.-t)) ...
   -(r+s+t-1.).*(-r+s+t-1.).*(-r-s+t-1.) ./ (2.*(1.-t)) ...
   t.*(-r+s+t-1.).*(-r-s+t-1.) ./ (1.-t) ...
   t.*(-r-s+t-1.).*(r-s+t-1.) ./ (1.-t) ...
   t.*(r-s+t-1.).*(r+s+t-1.) ./ (1.-t) ...
   t.*(r+s+t-1.).*(-r+s+t-1.) ./ (1.-t) ];
 
 Nr=[0.5*(-3*r.^2 + 4*r.*t - 3.0*r + s.^2 - t.^2 + 1.0*t)./(t - 1) ...
 0.5*r.*(2*s - 1.0)./(t - 1) ...
 0.5*(3*r.^2 + 4*r.*t - 3.0*r - s.^2 + t.^2 - 1.0*t)./(t - 1) ...
 -0.5*r.*(2*s + 1.0)./(t - 1), zeros(size(s)) ...
 0.5*(3*r.^2 - 2*r.*s - 2*r.*t + 2.0*r - s.^2 + 2*s.*t - 2.0*s - t.^2 + 2.0*t - 1.0)./(t - 1) ...
 0.5*(-3*r.^2 - 2*r.*s - 2*r.*t + 2.0*r + s.^2 - 2*s.*t + 2.0*s + t.^2 - 2.0*t + 1.0)./(t - 1) ...
 0.5*(-3*r.^2 + 2*r.*s - 2*r.*t + 2.0*r + s.^2 + 2*s.*t - 2.0*s + t.^2 - 2.0*t + 1.0)./(t - 1) ...
 0.5*(3*r.^2 + 2*r.*s - 2*r.*t + 2.0*r - s.^2 - 2*s.*t + 2.0*s - t.^2 + 2.0*t - 1.0)./(t - 1) ...
 t.*(-2*r + 2*t - 2.0)./(t - 1.0),  2*r.*t./(t - 1.0), t.*(-2*r - 2*t + 2.0)./(t - 1.0), 2*r.*t./(t - 1.0)];
  
 Ns=[0.5*s.*(2*r - 1.0)./(t - 1) ...
 0.5*(r.^2 - 3*s.^2 + 4*s.*t - 3.0*s - t.^2 + 1.0*t)./(t - 1) ...
 -0.5*s.*(2*r + 1.0)./(t - 1) ...
 0.5*(-r.^2 + 3*s.^2 + 4*s.*t - 3.0*s + t.^2 - 1.0*t)./(t - 1) ...
 zeros(size(s)) ...
 0.5*(-r.^2 - 2*r.*s + 2*r.*t - 2.0*r + 3*s.^2 - 2*s.*t + 2.0*s - t.^2 + 2.0*t - 1.0)./(t - 1) ...
 0.5*(-r.^2 + 2*r.*s - 2*r.*t + 2.0*r + 3*s.^2 - 2*s.*t + 2.0*s - t.^2 + 2.0*t - 1.0)./(t - 1) ...
 0.5*(r.^2 + 2*r.*s + 2*r.*t - 2.0*r - 3*s.^2 - 2*s.*t + 2.0*s + t.^2 - 2.0*t + 1.0)./(t - 1) ...
 0.5*(r.^2 - 2*r.*s - 2*r.*t + 2.0*r - 3*s.^2 - 2*s.*t + 2.0*s + t.^2 - 2.0*t + 1.0)./(t - 1) ...
 2*s.*t./(t - 1.0), t.*(-2*s + 2*t - 2.0)./(t - 1.0), 2*s.*t./(t - 1.0), t.*(-2*s - 2*t + 2.0)./(t - 1.0)];
 
 Nt=[1.0*(2.0*r.^3 - 1.0*r.^2 - 2.0*r.*s.^2 - 2.0*r.*t.^2 + 4.0*r.*t - 2.0*r + 1.0*s.^2 + 1.0*t.^2 - 2.0*t + 1.0)./(4.0*t.^2 - 8.0*t + 4.0) ...
 (-2.0*r.^2.*s + 1.0*r.^2 + 2.0*s.^3 - 1.0*s.^2 - 2.0*s.*t.^2 + 4.0*s.*t - 2.0*s + 1.0*t.^2 - 2.0*t + 1.0)./(4.0*t.^2 - 8.0*t + 4.0) ...
 (-2.0*r.^3 - 1.0*r.^2 + 2.0*r.*s.^2 + 2.0*r.*t.^2 - 4.0*r.*t + 2.0*r + 1.0*s.^2 + 1.0*t.^2 - 2.0*t + 1.0)./(4.0*t.^2 - 8.0*t + 4.0) ...
 1.0*(2.0*r.^2.*s + 1.0*r.^2 - 2.0*s.^3 - 1.0*s.^2 + 2.0*s.*t.^2 - 4.0*s.*t + 2.0*s + 1.0*t.^2 - 2.0*t + 1.0)./(4.0*t.^2 - 8.0*t + 4.0) ...
 4.0*t - 1.0 ...
 (-2.0*r.^3 + 2.0*r.^2.*s + 2.0*r.*s.^2 - 2.0*r.*t.^2 + 4.0*r.*t - 2.0*r - 2.0*s.^3 - 2.0*s.*t.^2 + 4.0*s.*t - 2.0*s + 4.0*t.^3 - 12.0*t.^2 + 12.0*t - 4.0)./(4.0*t.^2 - 8.0*t + 4.0) ...
 1.0*(2.0*r.^3 + 2.0*r.^2.*s - 2.0*r.*s.^2 + 2.0*r.*t.^2 - 4.0*r.*t + 2.0*r - 2.0*s.^3 - 2.0*s.*t.^2 + 4.0*s.*t - 2.0*s + 4.0*t.^3 - 12.0*t.^2 + 12.0*t - 4.0)./(4.0*t.^2 - 8.0*t + 4.0) ...
 1.0*(2.0*r.^3 - 2.0*r.^2.*s - 2.0*r.*s.^2 + 2.0*r.*t.^2 - 4.0*r.*t + 2.0*r + 2.0*s.^3 + 2.0*s.*t.^2 - 4.0*s.*t + 2.0*s + 4.0*t.^3 - 12.0*t.^2 + 12.0*t - 4.0)./(4.0*t.^2 - 8.0*t + 4.0) ...
 (-2.0*r.^3 - 2.0*r.^2.*s + 2.0*r.*s.^2 - 2.0*r.*t.^2 + 4.0*r.*t - 2.0*r + 2.0*s.^3 + 2.0*s.*t.^2 - 4.0*s.*t + 2.0*s + 4.0*t.^3 - 12.0*t.^2 + 12.0*t - 4.0)./(4.0*t.^2 - 8.0*t + 4.0) ...
 1.0*(1.0*r.^2 + 2.0*r.*t.^2 - 4.0*r.*t + 2.0*r - 1.0*s.^2 - 2.0*t.^3 + 5.0*t.^2 - 4.0*t + 1.0)./(1.0*t.^2 - 2.0*t + 1.0) ...
 (-1.0*r.^2 + 1.0*s.^2 + 2.0*s.*t.^2 - 4.0*s.*t + 2.0*s - 2.0*t.^3 + 5.0*t.^2 - 4.0*t + 1.0)./(1.0*t.^2 - 2.0*t + 1.0) ...
 1.0*(1.0*r.^2 - 2.0*r.*t.^2 + 4.0*r.*t - 2.0*r - 1.0*s.^2 - 2.0*t.^3 + 5.0*t.^2 - 4.0*t + 1.0)./(1.0*t.^2 - 2.0*t + 1.0) ...
 (-1.0*r.^2 + 1.0*s.^2 - 2.0*s.*t.^2 + 4.0*s.*t - 2.0*s - 2.0*t.^3 + 5.0*t.^2 - 4.0*t + 1.0)./(1.0*t.^2 - 2.0*t + 1.0)];
 
 out=struct('N',N,'Nr',Nr,'Ns',Ns,'Nt',Nt,'Nw',size(w,1), ...
      'NDN',zeros(size(N,2),4*size(N,1)),'NDNLabels',{{'',',x',',y',',z'}}, ...
      'jdet',zeros(size(N,1),1),'w',w,'Nnode',13, ...
      'DofLabels',{{'u','v','w'}},'xi',xi,'type','pyra13',...
      'Gmesh',Gdata(1),'Gedge',Gdata(2));

 if carg<=nargin % Possibly compute the derivatives OBSOLETE
      out=integrules('3d',out,varargin{carg});carg=carg+1;
 end



%% #MatrixRule ---------------------------------------------------------------
% Display and generation of matrix integration rules based on StrainDefinition
% and ConstitTopology
elseif comstr(Cam,'matrixrule')  

 opt=varargin{carg};carg=carg+1;
 sdx=opt.NDNLabels;
 su=opt.DofLabels;
 if isfield(opt,'PerField'); PerField=opt.PerField;
 else; PerField=size(opt.N,2);
 end

for j0=1:length(opt.StrainDefinition)
 if isempty(opt.StrainDefinition{j0}); opt.MatrixIntegrationRule{j0}=[];
 else
 if isfield(opt,'StrainLabels');ss=opt.StrainLabels{j0};
   dd=double(opt.ConstitTopology{j0});
 else; ss=opt.StrainString{j0};
   if length(opt.ConstitTopology)<j0; break;end
   constit=1:1024;dd=double(eval(opt.ConstitTopology{j0}));
 end
 ep=double(opt.StrainDefinition{j0});
 out=[];
 % Dof1 Dof2 NDN1 NDN2 Constit StepConstit StepNW NwIni
 for j1=1:size(dd,1);
  for j1b=find(ep(:,1)==j1)';
   for j2=1:size(dd,2); 
    for j2b=find(ep(:,1)==j2)';

   i2=sign(ep(j1b,2)*double(ep(j2b,2)))*dd(j1,j2);
   if i2~=0
    out(end+1,1:8)=[ ...
     (ep(j1b,3)-1)*PerField,(ep(j2b,3)-1)*PerField, ... % DDL1, DDL2
     (abs(ep(j1b,2))-1)*opt.Nw+abs(ep(j1b,4))-1, ... % NDNi
     (abs(ep(j2b,2))-1)*opt.Nw+abs(ep(j2b,4))-1, ... % NDNj
     sign(ep(j1b,2)*ep(j2b,2))*(dd(j1,j2)-1), ... % constit
     0,min(ep(j2b,5),ep(j2b,5)),ep(j2b,4)-1]; % ConsStep Nw W0
    if nargout==0
     fprintf( ...
      '   %2i %2i   %2i %2i  %3i    %i %i %i %% G(%s,%s) * %s%s * %s%s  \n',...
      out(end,:), ...
      ss{j1},ss{j2},su{abs(ep(j1b,3))},sdx{abs(ep(j1b,2))}, ...
      su{abs(ep(j2b,3))},sdx{abs(ep(j2b,2))});
    end
    end
    end
   end
  end
 end
 opt.MatrixIntegrationRule{j0}=int32(out);
 end % empty matrix or not
end
opt.MatrixIntegrationRule(end+1:3)={[]}; % Always viscous damping
if nargout>0;out=opt;end

%% #texstrain ----------------------------------------------------------------
% Display and generation of stress integration rules based on StrainDefinition
% and ConstitTopology
elseif comstr(Cam,'texstrain')

opt=varargin{carg};carg=carg+1;
fid=1; fprintf(fid,'\n');

for j0=1:length(opt.StrainDefinition)

 fprintf(fid,'\\begin{equation}\n');
 st=opt.StrainLabels{j0}; if isempty(st); continue;end
 st=sprintf('%s \\\\',st{:});
 fprintf('\\ve{\\ba{c}%s\\ea}\n',st(1:end-2))

 r1=opt.StrainDefinition{j0}; st=opt.StrainLabels{j0};
 r2=cell(length(unique(r1(:,1))),length(unique(r1(:,3))));
 for j1=1:size(r1,1)
  if r1(j1,2)>0;r2{r1(j1,1),r1(j1,3)}=sprintf('N%s',opt.NDNLabels{r1(j1,2)});
  else;r2{r1(j1,1),r1(j1,3)}=sprintf('-N%s',opt.NDNLabels{-r1(j1,2)});
  end
 end
 st='c'; fprintf(fid,'=\\ma{\\ba{%s}\n',st(ones(1,size(r2,1))));
 for j1=1:numel(r2); if isempty(r2{j1}); r2{j1}='0';end;end
 st=[' %s &']';st=st(:,ones(size(r2,2),1));st=st(:)';st(end+[0:5])='\\\\\n';
 r2=r2';st=sprintf(st,r2{:});st(end+[-2:1])='\ea}';
 fprintf(fid,'%s',st);
 st=opt.DofLabels;st=sprintf(' %s \\\\',st{:});
 fprintf('\n\\ve{\\ba{c}%s\\ea}\n\\end{equation}\n\n',st(1:end-2))

end
clear out;

  
%% #StressRule ---------------------------------------------------------------
% Display and generation of stress integration rules based on StrainDefinition
% and ConstitTopology
elseif comstr(Cam,'stressrule')  

opt=varargin{carg};carg=carg+1;
if carg<=nargin % Possibly evaluate stress elsewhere than integration point
 st=varargin{carg};carg=carg+1;
 if strcmp(opt.type,'mit4'); opt.type='quad4';end
 r1=integrules(opt.type,st);
 st={'N','Nr','Ns','Nt','NDN','jdet','w','Nnode','Nw'};
 for j1=1:length(st)
   try; 
    if isfield(r1,st{j1});eval(sprintf('opt.%s=r1.%s;',st{j1},st{j1}));end
   end
 end
end

if ~isfield(opt,'StressRule'); opt.StressRule={};end
nNode=size(opt.N,2);
for j0=1:length(opt.StrainDefinition)
 ss=opt.StrainLabels{j0};
 dd=double(opt.ConstitTopology{j0});
 ep=double(opt.StrainDefinition{j0});%[StrainRow,+-NDNBloc,DOF,NwStart,NwTot]
 out=[]; if isempty(ep);continue;end
 rule=[];
 % sigma = dd * ep
 for j1=1:size(dd,1);
 for j2=1:size(dd,2);
 for j2b=find(ep(:,1)==j2)';
  % StressRowValue NDNBlock ConstitCoefPos DOFOffSet nStres NwOffset sign]
  if dd(j1,j2)~=0
   rule(end+1,:)=[j1-1 (abs(ep(j2b,2))-1) dd(j1,j2)-1  ...
       (ep(j2b,3)-1) size(dd,1) ep(j2b,4)-1 sign(ep(j2b,2))];
  end % dd~=0
 end    
 end 
 end
 opt.StressRule{j0}=int32(rule);
end
if isfield(opt,'bas');opt.bas(9,opt.Nw)=0;end
if isfield(opt,'J');opt.J(4,opt.Nw)=0;end

if nargout>0;out=opt;else; integrules('texstress',opt); end

%% #texmacros ---------------------------------------------------------------
% Display and generation of stress integration rules based on StrainDefinition
% and ConstitTopology 
elseif comstr(Cam,'texmacros')


st={'You may need to define the following macros'
 '\newcommand{\ba}[1]{\begin{array}{#1}}'
'\newcommand{\ea}{\end{array}}'
'\newcommand{\ve}[1]{\left\{#1\right\}}'
'\newcommand{\ma}[1]{\left[#1\right]}';''};
fprintf('%s\n',st{:})
clear out; 

%% #texstress ----------------------------------------------------------------
% Display and generation of stress integration rules based on StrainDefinition
% and ConstitTopology, call is tested in t_thermal
elseif comstr(Cam,'texstress')

EC=varargin{carg};carg=carg+1;
fid=1; integrules('texmacros');fprintf(fid,'\n');

for j0=1:length(EC.StressRule)

 fprintf(fid,'\\begin{equation}\n');
 if isfield(EC,'StressLabels')&&j0<=length(EC.StressLabels);
     st=EC.StressLabels{j0};st=sprintf('%s \\\\',st{:});
 else;
   st=EC.StrainLabels{j0};st=sprintf('\\sigma %s \\\\',st{:});
 end
 fprintf('\\ve{\\ba{c}%s\\ea}\n',st(1:end-2));

 r1=double(EC.StressRule{j0});
 r2=cell(length(unique(r1(:,1))),length(unique(r1(:,4))));
 for j1=1:size(r1,1)
  st2=r2{r1(j1,1)+1,r1(j1,4)+1}; 
  [i1,i2]=find(EC.ConstitTopology{j0}==(r1(j1,3)+1));i1=i1(1);i2=i2(1);
  st3=sprintf('d_{%i,%i} N%s',i1,i2,EC.NDNLabels{r1(j1,2)/EC.Nw+1});
  if isempty(st2); r2{r1(j1,1)+1,r1(j1,4)+1}=st3;
  else;  r2{r1(j1,1)+1,r1(j1,4)+1}=sprintf('%s+%s',st2,st3);
  end
 end

 st='c'; fprintf(fid,'=\\ma{\\ba{%s}\n',st(ones(1,size(r2,1))));
 for j1=1:numel(r2); if isempty(r2{j1}); r2{j1}='0';end;end
 st=[' %s &']';st=st(:,ones(size(r2,2),1));st=st(:)';st(end+[0:5])='\\\\\n';
 r2=r2';st=sprintf(st,r2{:});st(end+[-2:1])='\ea}';
 fprintf(fid,'%s',st);
 st=EC.DofLabels;st=sprintf(' %s \\\\',st{:});
 fprintf('\n\\ve{\\ba{c}%s\\ea}\n\\end{equation}\n\n',st(1:end-2))

end
clear out;

% 2D Tensor transformation code generation
if 1==2

 st=cell(4,4); N=2;
 for j1=0:N-1;
     for j2=0:N-1
         for j3=0:N-1
             for j4=0:N-1;
  st1=st{j1*N+j4+1,j2*N+j3+1}; if isempty(st1); st1='';else; st1(end+1)='+';end
  st{j1*N+j4+1,j2*N+j3+1}=sprintf('%sJ[%i*%i+%i]*J[%i*%i+%i]',st1, ...
   j2,N,j1,j3,N,j4);
             end
         end
     end
 end

end

%% #PlotNodePos ----------------------------------------------------------
% Display node positions of a given software -----------------------------
elseif comstr(Cam,'plotnodepos') 

   ElemF=varargin{carg};carg=carg+1;
   model=femesh(['test' ElemF]);
   % force to order acording to .elt
   i1=feval(ElemF,'node');i1=model.Elt(2,i1);
   i2=[];i2(i1)=1:length(i1);model=feutil('renumber',model,i2);
   cg=feplot(100);feplot(model);fecom('showpatch');
   
   fecom(cg,'SetProp sel(1).fsProp','FaceAlpha',0.300);
   fecom(cg,'textnode','groupall','fontsize',20);fecom('triax');
    
% Node positions for gefdyn-------------------------------------------------
% OBSOLETE
elseif comstr(Cam,'gef'); [CAM,Cam]=comstr(CAM,4);

if nargin<2; w='def'; else; w=varargin{carg};carg=carg+1;end

if comstr(Cam,'bar1');      xi=[0 0 0;1 0 0];
elseif comstr(Cam,'beam1'); xi=[0 0 0;1 0 0];
elseif comstr(Cam,'beam3'); xi=[0 0 0;1 0 0;.5 0 0];
elseif comstr(Cam,'quad4'); xi=[1 1 0;-1 1 0;-1 -1 0;1 -1 0]; 
elseif comstr(Cam,'quadb');
 xi= [1 1 0;-1 1 0;-1 -1 0;1 -1 0;0 1 0;-1 0 0;0 -1 0; 1 0 0];
elseif comstr(Cam,'tria3');
 xi=[1 0 0;0 1 0;0 0 0]; %GEFDYN
elseif comstr(Cam,'tria6');
 xi=[1 0 0;0 1 0;0 0 0;.5 0.5 0;0 .5 0;0.5 0 0]; % GEFDYN 
else; error('unknown');
end
out=integrules(horzcat(Cam,'new'),w,xi);

% Integration rules for 2D jacobian -----------------------------
% OBSOLETE use BuildNdN
elseif comstr(Cam,'2d') 

  out=varargin{carg};carg=carg+1;
  x=varargin{carg};carg=carg+1;
      xr = out.Nr*x(:,1); xs = out.Ns*x(:,1);
      yr = out.Nr*x(:,2); ys = out.Ns*x(:,2);
      jdet = xr.*ys-xs.*yr;  jdet=jdet*sign(jdet(1));
      if any(jdet<0); disp('quad4 negative Jacobian');end
      i3=ones(size(out.Nr,2),1); % number of shape functions
      out.Nx = [ ys(:,i3).*out.Nr-yr(:,i3).*out.Ns]./jdet(:,i3);
      out.Ny = [-xs(:,i3).*out.Nr+xr(:,i3).*out.Ns]./jdet(:,i3);;
      out.jdet=jdet;

% Integration rules for 3D jacobian -----------------------------
% OBSOLETE use BuildNdN
elseif comstr(Cam,'3d') 

  out=varargin{carg};carg=carg+1;
  x=varargin{carg};carg=carg+1;
      Nr=out.Nr; Ns=out.Ns; Nt=out.Nt; 
      xr = Nr*x(:,1); xs = Ns*x(:,1); xt = Nt*x(:,1);
      yr = Nr*x(:,2); ys = Ns*x(:,2); yt = Nt*x(:,2);
      zr = Nr*x(:,3); zs = Ns*x(:,3); zt = Nt*x(:,3);

      cof11 = ys.*zt - yt.*zs;
      cof12 = yt.*zr - yr.*zt;
      cof13 = yr.*zs - ys.*zr;
      cof21 = zs.*xt - zt.*xs;
      cof22 = zt.*xr - zr.*xt;
      cof23 = zr.*xs - zs.*xr;
      cof31 = xs.*yt - xt.*ys;
      cof32 = xt.*yr - xr.*yt;
      cof33 = xr.*ys - xs.*yr;

      jdet = xr.*cof11+xs.*cof12+xt.*cof13;
      jdet=jdet*sign(jdet(1));
      if any(jdet<0); warning('negative Jacobian');end
      i3=ones(size(Nr,2),1); % number of shape functions
      out.jdet=jdet;jdet=jdet(:,i3);
      out.Nx=[Nr.*cof11(:,i3)+Ns.*cof12(:,i3)+Nt.*cof13(:,i3)]./jdet;
      out.Ny=[Nr.*cof21(:,i3)+Ns.*cof22(:,i3)+Nt.*cof23(:,i3)]./jdet; 
      out.Nz=[Nr.*cof31(:,i3)+Ns.*cof32(:,i3)+Nt.*cof33(:,i3)]./jdet; 

% ---------------------------------------------------------------------
elseif comstr(Cam,'gauss');[CAM,Cam]=comstr(CAM,6);

 if nargin==1&&isempty(Cam)
  out={'1d' 'lines';'q2d' '2D quadrangles';
       't2d','2D triangles';
       't3d','3D tetrahedrons '
       'p3d','3D pentahedrons (prisms)';
       'h3d','3D hexahedrons (brick)';
       'y3d','3D pyramids';
       };
  if nargout==0; 
   out=out';
   if exist('sdtweb','file')
     for j1=1:size(out,2);
         sdtweb('_link',sprintf('integrules(''Gauss %s'')',out{1,j1}),out{2,j1});
     end
   else;fprintf('\nAvailable rule families\n\n')
    fprintf('integrules(''Gauss %-4s'') %% %s\n',out{:});
   end
   clear out
  end
 elseif nargin==1
  out=QuadraturePoints(CAM,'list');
 else
  [out,out1]=QuadraturePoints(CAM,varargin{carg:end});
 end

%% #SeekInFile h125, h64 ----------------------------------------------------
else
 fname=sscanf(CAM,'%s',1);
 if exist(fname,'file');
  out=feval(fname,'rule',varargin{2:end});
  return
 elseif comstr(CAM,'@');out=eval(CAM);return;
 end
 error('''%s'' unknown',CAM);
end

%% #QuadraturePoints --------------------------------------------------------
function [w,iN,Gdata]=QuadraturePoints(topo,w,xi,varargin)

if nargin<2; w=-3;
elseif ischar(w); 
 w=lower(w); 
 if comstr(w,'node'); w=-2;
 elseif comstr(w,'center'); w=-1;
 elseif comstr(w,'def'); w=-3;
 end 
elseif length(w)==1; w=double(w); if isequal(w,0);w=-3;end
end
if nargin<3; xi=[];iN=[];
elseif nargin>3  % Allow for reordering of reference nodes
 xi1=varargin{1};carg=2;
 [xi2,i1,i2]=intersect(xi,xi1,'rows');
 iN=[];iN(i1)=i2;
else; iN=1:size(xi,1);
end % reference node reordering
Gdata={[],[]};
if ~ischar(w)&&size(w,2)==4; 

else;
switch topo

case '1d' 
%% #l1d - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 rules={ ...
  -3,[],'order default',[],[];
  -2,[xi ones(size(xi,1),1)/max(1,size(xi,1))],'node',[],[]
  -1,[.5 0 0 1],'center',[],[]
  2, [.5-sqrt(3)/6 0 0 .5; .5+sqrt(3)/6 0 0 .5],'2',[],[]; % order 2 integration GEFDYN
    % order 3 integration GEFDYN
  3, [.5-sqrt(3/20) 0 0 5/18
       .5 0 0 8/18
       .5+sqrt(3/20) 0 0 5/18],'3',[],[]};
  % order 4 integration GEFDYN
   r1=[ -.8611363115941 0 0 .3478548451375
        -.3399810435849 0 0 .6521451548625
         .3399810435849 0 0 .6521451548625
         .8611363115941 0 0 .3478548451375];r1(:,1)=r1(:,1)+1;r1=r1/2;
 rules(end+1,1:5)={4,r1,'4',[],[]};
 % a=;b=1/3*sqrt(5-2*sqrt(10/7))
 r1=[-1/3*sqrt(5+2*sqrt(10/7)) 0 0 (322-13*sqrt(70))/900
     -1/3*sqrt(5-2*sqrt(10/7)) 0 0 (322+13*sqrt(70))/900
    0 0 0 128/225
     1/3*sqrt(5-2*sqrt(10/7)) 0 0 (322+13*sqrt(70))/900
     1/3*sqrt(5+2*sqrt(10/7)) 0 0 (322-13*sqrt(70))/900
     ];r1(:,1)=r1(:,1)+1;r1=r1/2;
 rules(end+1,1:5)={5,r1,'5',[],[]};
 r1=[-.932469514203152 0 0 .171324492379170
     -.661209386466265 0 0 .360761573048139
     -.238619186083197 0 0 .467913934572691
      .238619186083197 0 0 .467913934572691
      .661209386466265 0 0 .360761573048139
      .932469514203152 0 0 .171324492379170
     ];r1(:,1)=r1(:,1)+1;r1=r1/2;
 rules(end+1,1:5)={6,r1,'6',[],[]};

 if isequal(w,-3) % defaults
   if size(xi,1)==2; w=2;else; w=3; end
 elseif length(w)==2&&w(1)==-1 % equal rule at a given number of points
   w=linspace(xi(1),xi(2),w(2))'*[1 0 0];
   w([1 end],4)=1/(2*size(w,1)-2);w([2:end-1],4)=2/(2*size(w,1)-2);
   return
 elseif ischar(w)&&strcmp(w,'list'); w=rules; return;
 elseif w>6 % equal rule at a given number of points
   rules(end+1,1:3)={w,gauss(w,0,1),sprintf('Gauss%i',w)};
   
 end

case 'q2d' 
%% #q2D quadrangle - - - - - - - - - - - - - - - - - - - - - -

if 1==2
 EC=integrules('beam1',3);w=EC.w(:,1)*2-1;
 [c,d]=meshgrid(EC.w(:,4)*2,EC.w(:,4)*2);[a,b]=meshgrid(w,w); 
 a=a';b=b';c=c';d=d';w=[a(:) b(:) zeros(numel(b),1) c(:).*d(:)];
 comstr(w,-30,struct('NoCplip',2))
end

FEw=1/sqrt(3);    xg13 = -sqrt(3/5); xg23 = 0. ; xg33 = sqrt(3/5) ;
    wg1 = 5/9; wg2 = 8/9; wg3 = 5/9;
rules ={ ...
  -3,[],'order default',[],[];
  -2,[xi ones(size(xi,1),1)*(4/max(1,size(xi,1)))],'node',...
     [Inf 113 117 97 100 52;1 2 3 4 0 0],[]; % nodes
  -1, [0 0 0 4],'center',[],[1;1;1;1];
  % order 2 integration GEFDYN
  102,  [-FEw -FEw 0 1;-FEw FEw 0 1;FEw -FEw 0 1;FEw FEw 0 1],'gefdyn 2x2',...
   [Inf 113 117 97 100 52;1 3 4 2 0 0],[1 3;3 4;4 2;2 1];
  % order 2 integration standard
  2,  [-FEw -FEw 0 1;FEw -FEw 0 1;FEw FEw 0 1;-FEw FEw 0 1],'standard 2x2',...
   [Inf 113 117 97 100 52;1 2 3 4 0 0],[1 2;2 3;3 4;4 1];
  % Q4WT standard 2x2 , center, midsides
  109, [-FEw -FEw 0 1;FEw -FEw 0 1;FEw FEw 0 1;-FEw FEw 0 1;0 0 0 4;
       -1 -1 0 1;1 -1 0 0;1 1 0 0;-1 1 0 0],'Q4WT',...
       [Inf 116 114 105 97 51;
       1 2 5 0 0 0;
       5 2 3 0 0 0;
       5 3 4 0 0 0;
       5 4 1 0 0 0;
       Inf 113 117 97 100 52
       1 6 7 2 0 0;
       2 7 8 3 0 0;
       3 8 9 4 0 0;
       4 9 6 1 0 0],[];
  % order 3 integration GEFDYN
  103,  [ xg13 xg13 0 wg1*wg1
          xg13 xg23 0 wg1*wg2
          xg13 xg33 0 wg1*wg3
          xg23 xg13 0 wg2*wg1
          xg23 xg23 0 wg2*wg2
          xg23 xg33 0 wg2*wg3
          xg33 xg13 0 wg3*wg1
          xg33 xg23 0 wg3*wg2
          xg33 xg33 0 wg3*wg3],'gefdyn 3x3',...
          [Inf 116 114 105 97 51
          1 2 5 4 0 0;
          4 5 8 7 0 0;
          5 6 9 8 0 0;
          6 5 2 3 0 0],[1 4 7;7 8 9;6 9 3;3 2 1];
      % order 4 integration GEFDYN
  104, [  -.8611363115941 -.8611363115941 0 .3478548451375*.3478548451375
         -.8611363115941 -.3399810435849  0 .6521451548625*.3478548451375
         -.8611363115941 .3399810435849  0 .6521451548625*.3478548451375
         -.8611363115941 .8611363115941  0 .3478548451375*.3478548451375

         -.3399810435849 -.8611363115941  0 .3478548451375*.6521451548625
         -.3399810435849 -.3399810435849  0 .6521451548625*.6521451548625
         -.3399810435849 .3399810435849  0 .6521451548625*.6521451548625
         -.3399810435849 .8611363115941  0 .3478548451375*.6521451548625
         
         .3399810435849 -.8611363115941  0 .3478548451375*.6521451548625
         .3399810435849 -.3399810435849  0 .6521451548625*.6521451548625
         .3399810435849 .3399810435849  0 .6521451548625*.6521451548625
         .3399810435849 .8611363115941  0 .3478548451375*.6521451548625
         
         .8611363115941 -.8611363115941 0  .3478548451375*.3478548451375
         .8611363115941 -.3399810435849 0  .6521451548625*.3478548451375
         .8611363115941 .3399810435849 0  .6521451548625*.3478548451375
         .8611363115941 .8611363115941 0  .3478548451375*.3478548451375] ...
       ,'gefdyn 4x4',...
        [Inf 116 114 105 97 51;
        1 2 6 5 0 0;2 3 7 6 0 0;3 4 8 7 0 0;
        5 6 10 9 0 0;6 7 11 10 0 0; 7 8 12 11 0 0;
        10 14 13 9 0 0;10 11 15 14 0 0;11 12 16 15 0 0],...
        [1 5 9 13;13 14 15 16;16 12 8 4;4 3 2 1]
        };  % --- ordre gefdyn

% 9 point rule plate rule
r9=[-sqrt(3/5) sqrt(3/5)  0 25/81 40/81 64/81];
rules(end+1,1:5)={9, r9([1 1 3 4;2 1 3 4;2 2 3 4;1 2 3 4;
  1 3 3 5;3 1 3 5;2 3 3 5;3 2 3 5;3 3 3 6]),'9 point',...
  [Inf 113 117 97 100 52;1 5 9 6 0 0;9 5 4 8 0 0;9 8 3 7 0 0;9 7 2 6 0 0],...
  [1 6 2;2 7 3;3 8 4;4 5 1]};
% 3x3 standard rule
  w2=[sqrt(3/5) 5/9 0 8/9]; % 3x3 integ.
  rules(end+1,1:5)={3, ...
   [-w2(1) -w2(1) 0 w2(2)^2;0 -w2(1) 0 w2(2)*w2(4);w2(1) -w2(1) 0 w2(2)^2;
    -w2(1) 0 0 w2(2)*w2(4);0 0 0 w2(4)^2;w2(1) 0 0 w2(2)*w2(4);
    -w2(1)  w2(1) 0 w2(2)^2;0  w2(1) 0 w2(2)*w2(4);w2(1)  w2(1) 0 w2(2)^2],...
  'standard 3x3',...
  [Inf 113 117 97 100 52;1 4 5 2 0 0;5 4 7 8 0 0;5 8 9 6 0 0;5 6 3 2 0 0],...
  [1 2 3;3 6 9;9 8 7;7 4 1]};
% 2x2 standard rule
   r1=1/sqrt(3);
  rules(end+1,1:5)={2, ...
     [ -r1  -r1 0  1; r1  -r1 0  1; r1   r1 0  1; -r1   r1 0  1 ], ...
    'standard 2x2',[Inf 113 117 97 100 52;1 2 3 4 0 0],[1 2;2 3;3 4;4 1]};
% combine 4 and 9 point rules

  w1=1/sqrt(3);w1=[-w1 -w1 0 1;w1 -w1 0 1;-w1 w1 0 1;w1 w1 0 1]; %2x2 integ.
  w2=[sqrt(3/5) 5/9 0 8/9]; % 3x3 integ.
  rules(end+1,1:5)={13, ...
   [w1;-w2(1) -w2(1) 0 w2(2)^2;0 -w2(1) 0 w2(2)*w2(4);w2(1) -w2(1) 0 w2(2)^2;
     -w2(1) 0 0 w2(2)*w2(4);0 0 0 w2(4)^2;w2(1) 0 0 w2(2)*w2(4);
   -w2(1)  w2(1) 0 w2(2)^2;0  w2(1) 0 w2(2)*w2(4);w2(1)  w2(1) 0 w2(2)^2],...
   '2x2 and 3x3',...
   [Inf 116 114 105 97 51;5 8 1 0 0 0;8 11 3 0 0 0;5 1 6 0 0 0;1 8 9 0 0 0;
   8 3 9 0 0 0;3 11 12 0 0 0;1 2 6 0 0 0;1 9 2 0 0 0;9 3 12 0 0 0;
   9 12 14 0 0 0;6 2 7 0 0 0;2 9 10 0 0 0;10 9 4 0 0 0;4 12 13 0 0 0;
   7 2 10 0 0 0;10 4 13 0 0 0],...
   [5 6 7;7 10 13;13 12 11;11 8 1]};
 rules(end+1,1:3)={19,[],'Gauss 9x9'};
 if isequal(w,-3) % defaults
   if size(xi,1)==4; w=2;else; w=9; end
 elseif length(w)==1&&any(w==[4 5 6 7 8 10 11 12 13 19]); % defaults
   if w==19;i1=find(vertcat(rules{:,1})==19);r1=gauss(9,0,1);st1='Gauss 9x9';
   else; i1=size(rules,1)+1; r1=gauss(w,0,1);st1=sprintf('Gauss%ix%i',w,w);end 
   %QuadraturePoints('1d',w);
   [c,d]=meshgrid(r1(:,4)*2,r1(:,4)*2);r1=r1*2-1;[a,b]=meshgrid(r1(:,1),r1(:,1)); 
   a=a';b=b';c=c';d=d';r1=[a(:) b(:) zeros(numel(b),1) c(:).*d(:)];
   rules(i1,1:3)={w,r1,st1};
 elseif strcmp(w,'list'); w=rules; return;
 end

case 't2d' 
%% #t2D triangle - - - - - - - - - - - - - - - - - - - - - -

  r3=[1/6 2/3 0 1/6];r4=[1/2 0 1/6];
rules={
  -3,[],'order default',[],[];
  -2,[xi ones(size(xi,1),1)*(.5/max(1,size(xi,1)))],'node',...
  [Inf 116 114 105 97 51;1 2 3 0 0 0],[]; % nodes
  -1,[1/3 1/3 0 1/2],'center',[],[1;1;1];
   3,r3([1 2 3 4;1 1 3 4;2 1 3 4]),'standard 3 point',...
   [Inf 116 114 105 97 51;1 2 3 0 0 0],[2 3;3 1;1 2];
 103,r4([1 2 1 3;1 1 2 3;2 1 1 3]),'GEFDyn 3 point',[],[]; 
 203,[.5 0 0 1/6;.5 .5 0 1/6;0 .5 0 1/6],'midside',...
 [Inf 116 114 105 97 51;1 2 3 0 0 0],[];
};
  % 6 points between [-1 1] documentation Code_Aster R30101b
  p1=0.11169079483905; p2=0.0549758718227661;
  a=0.445948490915965; b=0.091576213509771;
  w1=[2*b-1  1-4*b  0  p2;
     2*b-1  2*b-1  0  p2;
     1-4*b  2*b-1  0  p2;
     1-4*a  2*a-1  0  p1;
     2*a-1  1-4*a  0  p1; 
     2*a-1  2*a-1  0  p1 ];
   w1(:,1:2)=(w1(:,1:2)+1)/2;
 rules(end+1,1:5)={6,w1,'6 point aster',...
  [Inf 116 114 105 97 51;1 4 6 0 0 0;6 4 2 0 0 0;6 2 5 0 0 0;6 5 3 0 0 0],...
  [2 5 3;3 6 1;1 4 2]};


 if isequal(w,-3) % defaults
   if size(xi,1)==3; w=3;else; w=6;end
 elseif strcmp(w,'list'); w=rules; return;
 elseif length(w)==1&&~any(vertcat(rules{:,1})==w)&&w>2
   % regular grid
   r=repmat(linspace(0,1,w),w,1)+1;s=repmat(linspace(0,1,w)',1,w)+1;
   i1=find(triu(r));
   w=[2-r(i1) s(i1)-1 ones(size(i1))*[0 0 1/length(i1)]];
   %figure(1);plot(w1(:,1),w1(:,2),'o')
   return
 end

case 't3d' 
%% #t3D tetrahedron - - - - - - - - - - - - - - - - - - - - - -

  xg1 = 1/3 ; xg2 = 0.6 ; xg3 = 0.2  ; 
  wg1 = - 0.28125; wg2 = 0.2604166667;
  xG1 = 1/3 ; xG2 = 0.0597158717D0; xG3 = 0.4701420641D0; xG4 = xG3 ;
  xG5 = 0.7974269853D0; xG6 = 0.1012865073D0; xG7 = 0.1012865073D0;
  wG1 = 0.1125D0; wG2 = 0.0661970764D0 ; wG3 = wG2 ; wG4 = wG2 ; wG5 = 0.0629695903D0; wG6 = wG5 ; wG7 = wG5 ;

rules={
  -3,[],'order default',[],[];
  -2,[xi ones(size(xi,1),1)*(1/6/max(1,size(xi,1)))],'node',[],[]; % nodes
  -1,[.25 .25 .25 1/6],'center',[],[];
 % ordre d'integration 2
 102,[ 1/2 0 1/2 1/6; 0 1/2 1/2 1/6; 1/2 0 1/2 1/6],'gefdyn 2nd order',[],[];
  %
  % ordre d'integration 3
  %
 103,[ xg1 xg1 xg1 wg1
         xg2 xg3 xg3 wg2
         xg3 xg2 xg3 wg2
         xg3 xg3 xg2 wg2],'gefdyn order 3',[],[];
  %
  % ordre d'integration 4
  %
 104,[xG1 xG1 xG1 wG1
        xG2 xG3 xG3 wG2
        xG3 xG2 xG3 wG2
        xG3 xG3 xG2 wG2
        xG5 xG6 xG6 wG5
        xG6 xG5 xG6 wG5
        xG6 xG6 xG5 wG5],'gefdyn order 4',[],[]
 % 5 Point rule. Hughes page 174
  5,[1/4 1/4 1/4   -4/30;
     1/2 1/6 1/6   3/40;
     1/6 1/2 1/6   3/40;
     1/6 1/6 1/2   3/40;
     1/6 1/6 1/6   3/40],'5 point hughes',[],[];
};

% Default 4 point rules Code_Aster R30101b
w1=[(5-sqrt(5))/20 (5+3*sqrt(5))/20 1/24];
rules(end+1,1:3)={4,w1([1 1 1 3;1 1 2 3;1 2 1 3;2 1 1 3]),'4 point'};
% Default 15 point rules Code_Aster R30101b
      % This is ordered as in MODULEF

  w1=[(7-sqrt(15))/34 (13+3*sqrt(15))/34 (7+sqrt(15))/34 (13-3*sqrt(15))/34 ...
     (5-sqrt(15))/20 (5+sqrt(15))/20];
  a=(2665+14*sqrt(15))/226800; b=(2665-14*sqrt(15))/226800; c=5/567;

  w1=[  1/4 1/4 1/4         8/405;
       w1([ 1 1 1 ])        a ; %
       w1([ 2 1 1 ])        a ;
       w1([ 1 2 1 ])        a ;
       w1([ 1 1 2 ])        a ;
       w1([ 3 3 3 ])        b ; %
       w1([ 4 3 3 ])        b ;
       w1([ 3 4 3 ])        b ;
       w1([ 3 3 4 ])        b ;
       w1([ 5 6 6 ])        c ; %
       w1([ 5 5 6 ])        c ;
       w1([ 6 5 5 ])        c ;
       w1([ 6 6 5 ])        c ;
       w1([ 6 5 6 ])        c ;
       w1([ 5 6 5 ])        c ];
rules(end+1,1:3)={15,w1,'15 point'};

 if isequal(w,-3) % defaults
   if size(xi,1)==4; w=4;else; w=15;end
 elseif strcmp(w,'list'); w=rules; return;
 end

case 'y3d' 
%% #y3d 3D pyramid - - - - - - - - - - - - - - - - - - - - - -
 h1=0.1531754163448146;h2=0.6372983346207416;r2=2/15/.999;
 
rules={-3,[.5 0 h1 2/15;0 .5 h1 2/15;-.5 0 h1 2/15;0 -.5 h1 2/15; 0 0 h2 2/15],'order default',[],[]% FPG5
       -2,[1 0 0 r2 ;0 1 0 r2;-1 0 0 r2 ;0 -1 0 r2;0 0 1 r2]*.999,'AtNode',[],[]};

case 'p3d' 
%% #p3D pentahedron - - - - - - - - - - - - - - - - - - - - - -

rules={
 -3,[],'order default',[],[];
 -2,[xi ones(size(xi,1),1)*(.5/max(1,size(xi,1)))],'node',[],[];
 -1,[1/3 1/3 .5 1/2],'center',[],[];
};
   w1=1/sqrt(3);w1=[-w1 -w1 -w1 1;w1 -w1 -w1 1;-w1 w1 -w1 1;w1 w1 -w1 1;
               -w1 -w1 w1 1;w1 -w1 w1 1;-w1 w1 w1 1;w1 w1 w1 1;];
   w1(:,1:3)=(w1(:,1:3)+1)/2;w1(:,4)=w1(:,4)/8;
rules(end+1,1:3)={15,w1,'2x2x2'};
% Modulef rule (documented in OpenFEM)
  w1=[1/6 4/6 .5-.5/sqrt(3) .5+.5/sqrt(3) 1/12];  
  w1=w1([1 1 3 5;
       2 1 3 5;
       1 2 3 5;
       1 1 4 5;
       2 1 4 5;
       1 2 4 5]);
rules(end+1,1:3)={6,w1,'modulef 6'};

  b=[(9-2*sqrt(15))/21 (9+2*sqrt(15))/21 (6+sqrt(15))/21 (6-sqrt(15))/21  ...
     .5*(1-sqrt(3/5)) .5 .5*(1+sqrt(3/5))     1/3];
  w1=[(155-sqrt(15))/2400 5/18 (155+sqrt(15))/2400 9/80 8/18 ];
 
  w1=[ b([4 4 5]) w1(1)*w1(2) ;
        b([2 4 5]) w1(1)*w1(2) ;
        b([4 2 5]) w1(1)*w1(2) ;
        b([3 1 5]) w1(3)*w1(2) ;
        b([3 3 5]) w1(3)*w1(2) ;
        b([1 3 5]) w1(3)*w1(2) ;
        b([8 8 5]) w1(4)*w1(2) ;
        b([4 4 6]) w1(1)*w1(5) ;
        b([2 4 6]) w1(1)*w1(5) ;
        b([4 2 6]) w1(1)*w1(5) ;
        b([3 1 6]) w1(3)*w1(5) ;
        b([3 3 6]) w1(3)*w1(5) ;
        b([1 3 6]) w1(3)*w1(5) ;
        b([8 8 6]) w1(4)*w1(5) ;
        b([4 4 7]) w1(1)*w1(2) ;
        b([2 4 7]) w1(1)*w1(2) ;
        b([4 2 7]) w1(1)*w1(2) ;
        b([3 1 7]) w1(3)*w1(2) ;
        b([3 3 7]) w1(3)*w1(2) ;
        b([1 3 7]) w1(3)*w1(2) ;
        b([8 8 7]) w1(4)*w1(2) ];

rules(end+1,1:3)={21,w1,'modulef 21'};
% Gefdyn penta6, penta15
 w1=[ .3333333333  .3333333333  -.774596669241483  .06250000000000
      .0597158717  .4701420641  -.774596669241483  .03677615355236
      .4701420641  .0597158717  -.774596669241483  .03677615355236
      .4701420641  .4701420641  -.774596669241483  .03677615355236
      .7974269853  .1012865073  -.774596669241483  .03498310570690
      .1012865073  .7974269853  -.774596669241483  .03498310570690
      .1012865073  .1012865073  -.774596669241483  .03498310570690
      .3333333333  .3333333333   .0                .1
      .0597158717  .4701420641   .0                .05884184568378
      .4701420641  .0597158717   .0                .05884184568378
      .4701420641  .4701420641   .0                .05884184568378
      .7974269853  .1012865073   .0                .05597296913103
      .1012865073  .7974269853   .0                .05597296913103
      .1012865073  .1012865073   .0                .05597296913103
      .3333333333  .3333333333   .774596669241483  .06250000000000
      .0597158717  .4701420641   .774596669241483  .03677615355236
      .4701420641  .0597158717   .774596669241483  .03677615355236
      .4701420641  .4701420641   .774596669241483  .03677615355236
      .7974269853  .1012865073   .774596669241483  .03498310570690
      .1012865073  .7974269853   .774596669241483  .03498310570690
      .1012865073  .1012865073   .774596669241483  .03498310570690 ];
  w1=[w1(:,1:2) .5*(w1(:,3)+1) 2*w1(:,4)];
  rules(end+1,1:3)={121,w1,'gefdyn penta'};

 if isequal(w,-3) % defaults
   if size(xi,1)==6; w=6;else; w=21;end
 elseif strcmp(w,'list'); w=rules; return;
 end

case 'h3d' 
%% #h3D hexahedron - - - - - - - - - - - - - - - - - - - - - -

w1=1/sqrt(3);
    % order 3 integration
    xg13 = -.7745966692415;     xg23 = 0.;     xg33 = .7745966692415;
    wg1= 5/9; wg2=8/9; wg3=5/9;

    xG14 = -.8611363115941; xG24 = -.3399810435849; xG34 = .3399810435849; xG44 = .8611363115941;
    wG1 = .3478548451375 ; wG2 = .6521451548625; wG3 = .6521451548625; wG4 = .3478548451375 ;

rules={
 -3,[],'order default',[],[];
 -2,[xi ones(size(xi,1),1)*(8/max(1,size(xi,1)))],'node',[],[];
 -1,[0 0 0 8],'center',[],[];
% GEFDYN
%     %  des tests a ajouter !!!
    % order 2 integration
 102,[-w1 -w1 -w1 1;-w1 -w1 w1 1;-w1 w1 -w1 1;-w1 w1 w1 1;
               w1 -w1 -w1 1;w1 -w1 w1 1;w1 w1 -w1 1;w1 w1 w1 1;],'gefdyn 2',[],[];
 103, [xg13 xg13 xg13 wg1*wg1*wg1
         xg13 xg13 xg23 wg1*wg1*wg2
         xg13 xg13 xg33 wg1*wg1*wg3
         xg13 xg23 xg13 wg1*wg2*wg1
         xg13 xg23 xg23 wg1*wg2*wg2
         xg13 xg23 xg33 wg1*wg2*wg3
         xg13 xg33 xg13 wg1*wg3*wg1
         xg13 xg33 xg23 wg1*wg3*wg2
         xg13 xg33 xg33 wg1*wg3*wg3
         
         xg23 xg13 xg13 wg2*wg1*wg1
         xg23 xg13 xg23 wg2*wg1*wg2
         xg23 xg13 xg33 wg2*wg1*wg3
         xg23 xg23 xg13 wg2*wg2*wg1
         xg23 xg23 xg23 wg2*wg2*wg2
         xg23 xg23 xg33 wg2*wg2*wg3
         xg23 xg33 xg13 wg2*wg3*wg1
         xg23 xg33 xg23 wg2*wg3*wg2
         xg23 xg33 xg33 wg2*wg3*wg3
         
         xg33 xg13 xg13 wg3*wg1*wg1
         xg33 xg13 xg23 wg3*wg1*wg2
         xg33 xg13 xg33 wg3*wg1*wg3
         xg33 xg23 xg13 wg3*wg2*wg1
         xg33 xg23 xg23 wg3*wg2*wg2
         xg33 xg23 xg33 wg3*wg2*wg3
         xg33 xg33 xg13 wg3*wg3*wg1
         xg33 xg33 xg23 wg3*wg3*wg2
         xg33 xg33 xg33 wg3*wg3*wg3],'gefdyn 3',[],[];
    % order 4 integration
  104, [   xG14 xG14 xG14 wG1*wG1*wG1;
           xG14 xG14 xG24 wG1*wG1*wG2
           xG14 xG14 xG34 wG1*wG1*wG3
           xG14 xG14 xG44 wG1*wG1*wG4
           xG14 xG24 xG14 wG1*wG2*wG1
           xG14 xG24 xG24 wG1*wG2*wG2
           xG14 xG24 xG34 wG1*wG2*wG3
           xG14 xG24 xG44 wG1*wG2*wG4
           xG14 xG34 xG14 wG1*wG3*wG1
           xG14 xG34 xG24 wG1*wG3*wG2
           xG14 xG34 xG34 wG1*wG3*wG3
           xG14 xG34 xG44 wG1*wG3*wG4
           xG14 xG44 xG14 wG1*wG4*wG1
           xG14 xG44 xG24 wG1*wG4*wG2
           xG14 xG44 xG34 wG1*wG4*wG3
           xG14 xG44 xG44 wG1*wG4*wG4
           xG24 xG14 xG14 wG2*wG1*wG1
           xG24 xG14 xG24 wG2*wG1*wG2
           xG24 xG14 xG34 wG2*wG1*wG3
           xG24 xG14 xG44 wG2*wG1*wG4
           xG24 xG24 xG14 wG2*wG2*wG1
           xG24 xG24 xG24 wG2*wG2*wG2
           xG24 xG24 xG34 wG2*wG2*wG3
           xG24 xG24 xG44 wG2*wG2*wG4
           xG24 xG34 xG14 wG2*wG3*wG1
           xG24 xG34 xG24 wG2*wG3*wG2
           xG24 xG34 xG34 wG2*wG3*wG3
           xG24 xG34 xG44 wG2*wG3*wG4
           xG24 xG44 xG14 wG2*wG4*wG1
           xG24 xG44 xG24 wG2*wG4*wG2
           xG24 xG44 xG34 wG2*wG4*wG3
           xG24 xG44 xG44 wG2*wG4*wG4
           xG34 xG14 xG14 wG3*wG1*wG1
           xG34 xG14 xG24 wG3*wG1*wG2
           xG34 xG14 xG34 wG3*wG1*wG3
           xG34 xG14 xG44 wG3*wG1*wG4
           xG34 xG24 xG14 wG3*wG2*wG1
           xG34 xG24 xG24 wG3*wG2*wG2
           xG34 xG24 xG34 wG3*wG2*wG3
           xG34 xG24 xG44 wG3*wG2*wG4
           xG34 xG34 xG14 wG3*wG3*wG1
           xG34 xG34 xG24 wG3*wG3*wG2
           xG34 xG34 xG34 wG3*wG3*wG3
           xG34 xG34 xG44 wG3*wG3*wG4
           xG34 xG44 xG14 wG3*wG4*wG1
           xG34 xG44 xG24 wG3*wG4*wG2
           xG34 xG44 xG34 wG3*wG4*wG3
           xG34 xG44 xG44 wG3*wG4*wG4
           xG44 xG14 xG14 wG4*wG1*wG1
           xG44 xG14 xG24 wG4*wG1*wG2
           xG44 xG14 xG34 wG4*wG1*wG3
           xG44 xG14 xG44 wG4*wG1*wG4
           xG44 xG24 xG14 wG4*wG2*wG1
           xG44 xG24 xG24 wG4*wG2*wG2
           xG44 xG24 xG34 wG4*wG2*wG3
           xG44 xG24 xG44 wG4*wG2*wG4
           xG44 xG34 xG14 wG4*wG3*wG1
           xG44 xG34 xG24 wG4*wG3*wG2
           xG44 xG34 xG34 wG4*wG3*wG3
           xG44 xG34 xG44 wG4*wG3*wG4
           xG44 xG44 xG14 wG4*wG4*wG1
           xG44 xG44 xG24 wG4*wG4*wG2
           xG44 xG44 xG34 wG4*wG4*wG3
           xG44 xG44 xG44 wG4*wG4*wG4],'gefdyn 4',[],[];
};

   w1=[5/18  8/18  5/18];
   b=[(1-sqrt(3/5))/2 .5 (1+sqrt(3/5))/2 (1-sqrt(3/5))/2];
   w1=[ b([1 1 1]) w1(1)*w1(1)*w1(1) ;
         b([1 1 2]) w1(1)*w1(1)*w1(2) ;
         b([1 1 3]) w1(1)*w1(1)*w1(3) ;
         b([1 2 1]) w1(1)*w1(2)*w1(1) ;
         b([1 2 2]) w1(1)*w1(2)*w1(2) ;
         b([1 2 3]) w1(1)*w1(2)*w1(3) ;
         b([1 3 1]) w1(1)*w1(3)*w1(1) ;
         b([1 3 2]) w1(1)*w1(3)*w1(2) ;
         b([1 3 3]) w1(1)*w1(3)*w1(3) ; 
         b([2 1 1]) w1(2)*w1(1)*w1(1) ;
         b([2 1 2]) w1(2)*w1(1)*w1(2) ;
         b([2 1 3]) w1(2)*w1(1)*w1(3) ;
         b([2 2 1]) w1(2)*w1(2)*w1(1) ;
         b([2 2 2]) w1(2)*w1(2)*w1(2) ;
         b([2 2 3]) w1(2)*w1(2)*w1(3) ;
         b([2 3 1]) w1(2)*w1(3)*w1(1) ;
         b([2 3 2]) w1(2)*w1(3)*w1(2) ;
         b([2 3 3]) w1(2)*w1(3)*w1(3) ;
         b([3 1 1]) w1(3)*w1(1)*w1(1) ;
         b([3 1 2]) w1(3)*w1(1)*w1(2) ;
         b([3 1 3]) w1(3)*w1(1)*w1(3) ;
         b([3 2 1]) w1(3)*w1(2)*w1(1) ;
         b([3 2 2]) w1(3)*w1(2)*w1(2) ;
         b([3 2 3]) w1(3)*w1(2)*w1(3) ;
         b([3 3 1]) w1(3)*w1(3)*w1(1) ;
         b([3 3 2]) w1(3)*w1(3)*w1(2) ;
         b([3 3 3]) w1(3)*w1(3)*w1(3) ;
       ];
   w1(:,1:3)=2*w1(:,1:3)-1;w1(:,4)=w1(:,4)*8;
 rules(end+1,1:3)={3,w1,'3x3x3'};
 %    Z=[-.774596669241483 0  .774596669241483]';
 %   W=[5/9,8/9,5/9]';
 %   i1=[1 1 1;1 1 2;1 1 3;1 2 1;1 2 2;1 2 3;1 3 1;1 3 2;1 3 3;
 %    2 1 1;2 1 2;2 1 3;2 2 1;2 2 2;2 2 3;2 3 1;2 3 2;2 3 3;
 %    3 1 1;3 1 2;3 1 3;3 2 1;3 2 2;3 2 3;3 3 1;3 3 2;3 3 3];
 %   w=[Z(i1) prod(W(i1),2)];
    
   w1=1/sqrt(3);w1=[-w1 -w1 -w1 1;w1 -w1 -w1 1;-w1 w1 -w1 1;w1 w1 -w1 1;
               -w1 -w1 w1 1;w1 -w1 w1 1;-w1 w1 w1 1;w1 w1 w1 1;];
 rules(end+1,1:3)={2,w1,'2x2x2'};
 if isequal(w,-4) % spectral element
   rules(end+1,1:3)={-4,feval(q16p('@N_Nr_Ns_Nt'),round(size(xi,1).^(1/3))),'spectral'};
 elseif isequal(w,-3) % defaults
   if size(xi,1)==8; w=2; else; w=3;end
 elseif strcmp(w,'list'); w=rules; return;
 end

 
otherwise; error('Not supported');
end

 if size(w,2)==4  % standard rule selection
 elseif length(w)==1; 
   for j1=1:size(rules,1)
    if isequal(w,rules{j1,1}); w=rules{j1,2};Gdata=rules(j1,4:5); break;end
   end
   if length(w)==1;
     fprintf('No match for %g using %g,''%s''\n',w,rules{end,1},rules{end,3})
       w=rules{end,2};Gdata=rules(end,4:5);
   end
 else
   error('not a valid case');
 end % standard rule selection


end; % w provided or not


% ---------------------------------------------------------------------
function writetoscreen(in)


fprintf('%s = [ ',inputname(1));
for j1=1:size(in,1)
 i0=0;
 for j2=1:size(in,2);
  st1=sprintf('%20.16e ',in(j1,j2)); i0=i0+1;
  st1=strrep(st1,'e-00','e-'); st1=strrep(st1,'e+00','e+');
  st1=strrep(st1,'00000000000000e','e');
  fprintf(st1);
  if i0==3&&j2<size(in,2); fprintf(', ...\n'); 
      i0=0; end
 end
 fprintf(';\n');
end
 fprintf('];%%%s\n',inputname(1));


function w=gauss(N,a,b)
%
% Derived from function posted by Greg von Winckel - 02/25/2004
% EC=integrules('beam1',3);EC.w-feval(integrules('@gauss'),3,0,1)
%
% This script is for computing definite integrals using Legendre-Gauss 
% Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
% [a,b] with truncation order N

N=N-1;
N1=N+1; N2=N+2;
xu=linspace(-1,1,N1)';
% Initial guess
y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);
% Legendre-Gauss Vandermonde Matrix
L=zeros(N1,N2);
% Derivative of LGVM
Lp=zeros(N1,N2);
% Compute the zeros of the N+1 Legendre Polynomial
% using the recursion relation and the Newton-Raphson method
y0=2;
% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0))>eps
    L(:,1)=1;    Lp(:,1)=0;    
    L(:,2)=y; Lp(:,2)=1;   
    for k=2:N1
        L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);       
    y0=y; y=y0-L(:,N2)./Lp;
end
if nargin==1; a=-1; b=1; end% Linear map from[-1,1] to [a,b]
x=(a*(1-y)+b*(1+y))/2;      
% Compute the weights
w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
w=sortrows([x*[1 0 0] w],1);
