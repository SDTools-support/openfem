function [out,out1,out2]=beam1(CAM,varargin);

%BEAM1   element function for the 2-node 12 DOF beam element
%
%	As all element functions (see ELEM0), BEAM1 is called by FE_MK for
%	model assembly, FEPLOT for structural deformation visualization, ...
%
%	In model description matrices, a group of BEAM1 elements starts with a
%	header row [Inf  abs('beam1') 0 ...] followed by element property rows
%       ELT following the format
%	   [n1 n2 MatId ProId  (vx vy vz EltId p1 p2 x1 y1 z1 x2 y2 z2)]
%	   [n1 n2 MatId ProId  (nr 0  0  EltId p1 p2 x1 y1 z1 x2 y2 z2)]
%         with required values
%	   n1,n2  node numbers of the tips of the beam element
%	   MatId  material property identification number
%	   ProId  section property identification number
%	   NR     number of node not in the beam direction defining the first
%		  bending plane (default node is 1.5 1.5 1.5) or the alternate
%		  format
%	   vx,vy,vz components of a vector in plane 1 (not collinear to the
%		  beam axis)
%	  and optional values
%	   EltId  element identifier (optional, default 0)
%	   p1,p2  pin flags (list of released DOFs, e.g. 456 is pinned)
%	   x1, ... components of offset vector at node 1
%	   x2, ... components of offset vector at node 2
%
%      PL material property matrix. Isotropic materials, see HELP M_ELASTIC
%         are the only supported.
%      IL element property matrix, see P_BEAM
%
%	See sdtweb      bar1, eltfun, elem0
%	See also help   bar1, p_beam, elem0


%	Etienne Balmes
%       Copyright (c) 2001-2024 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license


%#ok<*NOSEM,*NASGU,*ASGLU>
% standard calls with one input argument
persistent Beam1In
if isempty(Beam1In)
   % Default is now to ignore inertia correction. 
   % To revert to previous status : setpref('OpenFEM','Beam1In',1)
   Beam1In= sdtdef('OpenFEM.Beam1In-safe',0); 
end
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1);
 %%  #BuildConstit, Integ, and Elmap for later integration
 if comstr(Cam,'integinfo')
  ElemF='beam1';  % Bypass because bars can have p_beam properties
  [out,out1,out2]= ... % sdtweb p_solid('constitbeam')
     p_solid('buildconstit 1 2',[varargin{1}],varargin{2:end});
 elseif comstr(Cam,'matcall');  
     constit=varargin{2};
     if ~isempty(constit)&&constit(1)<0;[out,out1]=elem0(CAM,varargin{:});
     else; out=beam1('call'); out1=0; % CallSymFlag
     end
 elseif comstr(Cam,'dofcall');out=elem0('dofcall'); % variable field elements
 elseif comstr(Cam,'call')  % call for matrix assembly
   out='[k1,m1]=beam1(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,elmap,node);';
 elseif comstr(Cam,'rhscall') % call for load assembly
   out='be=beam1(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,elmap,node,defe);';

 elseif  comstr(Cam,'state');   out = ''; % no state evaluation yet
 elseif comstr(Cam,'groupinit');out=elem0(CAM,'beam1');
 elseif comstr(Cam,'constants');
  if nargin<4; p_solid('constsolid','bar1',[],[]);return;end
  integ=varargin{2};
  if varargin{3}(1)==-3||(length(integ)>3&&integ(3)==integ(4)) % scalar
   [out,i2,out2]=p_solid('constsolid','bar1',varargin{2:end});
   out1=varargin{1};out1(4,:)=i2; % Tell MatrixIntegration this is 3d 
  else; out2=[]; out=[];out1=varargin{1};out1(4,:)=13; 
  end
 elseif  comstr(Cam,'node');   out = [1 2];
 elseif  comstr(Cam,'prop');   out = [3 4 8];
 elseif  comstr(Cam,'dof');    out=[1+(1:6)/100 2+(1:6)/100];
 elseif  comstr(Cam,'edge');   out =[1 2];
 elseif  comstr(Cam,'face');   out =[];
 elseif  comstr(Cam,'sci_face'); out = [1 2 2];  
 elseif  comstr(Cam,'line');   out =[1 2];
 elseif  comstr(Cam,'patch');  out= [1 2];
 elseif  comstr(Cam,'beam1in'); 
   if nargin>1; Beam1In=varargin{1};end;out=Beam1In;
 elseif  comstr(Cam,'parent'); out = 'beam1';

 % Basic matrix test of the element - - - - - - - - - - - - - - - - - -
 elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,5); 

   model=femesh('testbeam1');
   model=fe_mk(model);
   out=model.K([2 1]);  
   
   % some orient check
   wd=getpref('SDTools','sdtdata','~');
   model.il(4,4:5)=0.1333*1e-7*[0.5 2];
   r1=[0 0 1; 3 0 0; 0 1 0];
   for j1=1:size(r1,1)
    model.Elt(2,5:7)=r1(j1,1:3);
    model=fe_mk(model);
    %K0=model.K{2}; save(sprintf('O:\\sdtdata\\ref_test\\beam1\\test%i_K.mat',j1),'K0')
    try
     if ~isempty(strfind(wd,'ref_test'));wd=fileparts(wd);end
     fname=fullfile(wd,sprintf('/ref_test/beam1/test%i_K.mat',j1));
     load(fname,'K0')
     if norm(full(model.K{2}-K0))/norm(full(K0))>5e-6; 
      error('beam1 stiffness mismatch, test %i',j1); 
     end
    catch
     sdtw('_nb','XXX can''t find %s. Test ignored.',fname); 
    end
   end
      
   if nargout==0;disp('TestMat passed');end

 elseif  comstr(Cam,'test');   [CAM,Cam] = comstr(Cam,5);

   [out,out1]=femesh(strcat(['teststruct beam1' Cam]));

 % Basic tests of the element - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(CAM,'cvs')
  out='$Revision: 1.83 $  $Date: 2024/10/16 08:16:35 $'; return;
 else;sdtw('''%s'' unknown',CAM);
 end

return
end % of standard calls with one input argument

% element matrix assembly - - - - - - - - - - - - - - - - - - - - - - - - - -
nodeE=CAM; 
elt=varargin{1}; 
point=double(varargin{2});
integ=varargin{3};
%constit=varargin{4};
constit=varargin{4}; %sdtweb p_solid('constitbeam')
elmap=varargin{5};
node=varargin{6};
typ=point(5);

if isa(elmap,'int32'); elmap = double(elmap); end

%if nargin<6 NNode(node(:,1))=1:size(node,1); end
%idof = [(1:6)/100+elt(1,1) (1:6)/100+elt(1,2)];
if size(nodeE,2)>4 
 node=nodeE; NNode(node(:,1))=1:size(node,1);
 nodeE = nodeE(NNode(elt(1,1:2)),[5:7 1]);
end
x=nodeE(:,1:3);

if size(elt,2)<5; x0=[1.5 1.5 1.5];  % if nothing,careful elt(6:7) may exist
elseif elt(5)<0||rem(elt(5),1)||(size(elt,2)>5&&elt(6))|| ...
        (size(elt,2)>6&&elt(7)) % normal given
   if size(elt,2)<7;elt(1,7)=0;end;x0 = x(1,:)+elt(5:7);     
else;x0 = find(node(:,1)==elt(5));
  if isempty(x0); x0=[1.5 1.5 1.5];
  else;x0 = node(x0,5:7); end
end

% off-set if any - - - - - - - - - - - - - - - - - - - - - - - - -
tr = eye(12,12);off=0;

r3=elt(11:min(13,length(elt)));
if ~isempty(r3)&&any(r3)
  r3(end+1:3)=0;
  tr([51 61 38])= elt(11:13);%[3 5;1 6;2 4];ans(:,1)+12*(ans(:,2)-1)
  tr([62 39 49])=-elt(11:13);%[2 6;3 4;1 5];ans(:,1)+12*(ans(:,2)-1);ans'
  off=1;
end
r3=elt(14:min(16,length(elt)));
if ~isempty(r3)&&any(r3)
  r3(end+1:3)=0;
  tr([129 139 116])= r3;%[3 5;1 6;2 4]+6;ans(:,1)+12*(ans(:,2)-1)
  tr([140 117 127])=-r3;%[2 6;3 4;1 5]+6;ans(:,1)+12*(ans(:,2)-1);ans'
  off=1;
end

% basis function - - - - - - - - - - - - - - - - - - - - - - - - - - - -

l = norm(x(2,:)-x(1,:)); l2 = l^2;
if l==0; error('beam with zero length, NodeId %i-%i',elt(1),elt(2)); end

x=sp_util('basis',x(2,:)-x(1,:),x0-x(2,:))';


%x(1,:) = x(2,:)-x(1,:);x(2,:) = x0-x(2,:);
%x(1,:) = x(1,:)/norm(x(1,:));
%x(2,:) = x(2,:) - (x(1,:)*x(2,:)')*x(1,:);
%x(2,:) = x(2,:)/norm(x(2,:));
%x(3,:) = [ -x(2,2)*x(1,3)+x(1,2)*x(2,3) ...
%            -x(2,3)*x(1,1)+x(1,3)*x(2,1) ...
%            -x(2,1)*x(1,2)+x(1,1)*x(2,2)];

%[pe,ie]=fe_mat(1,elt(3:4),pl,il);

if ~any([0 1 2 3 5 7 8 9 70 100]==typ)
 warning('OpenFEM:UNSUP','Matrix type %i not supported by beam1',typ)
else
 pe=constit(point(7)+(1:7));pe=pe(:)';
 ie=constit(point(7)+(8:15));ie=ie(:)';
end
Tpin=[];

%% #stiffness - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if any([0 1 5]==typ) || (size(elt,2)>9&&any(elt(9:10)))

   k = zeros(12,12);
   ind = [1 7 73 79];     k(ind) = [1 -1 -1 1]*pe(1,1)*ie(1,4)/l; % 1 7
   % 4 10, torsionnal stiffness GJ/l
   ind = [40 46 112 118]; k(ind) = [1 -1 -1 1]*pe(1,4)*ie(1,1)/l; % 4 10
   ind = [14 18 20 24 62 66 68 72 86 90 92 96 134 138 140 144]; %[2 6 8 12];
   if ie(5)==0 % bernoulli
      k(ind) = (pe(1,1)*ie(1,2)/l)* ...
	[12/l2    6/l   -12/l2   6/l     6/l   4  -6/l  2 ...
	 -12/l2  -6/l    12/l2  -6/l     6/l   2  -6/l  4];
   else % timoshenko with reduced integration
      k(ind) = ...
        (pe(1,1)*ie(1,2)/l)*[0 0 0 0 0 1 0 -1 0 0 0 0 0 -1 0 1] + ... %flexion
        (pe(4)*ie(4)*ie(5)/l)*[1 l/2 -1 l/2   l/2 l2/4 -l/2 l2/4 ...
                           -1 -l/2 1 -l/2   l/2 l2/4 -l/2 l2/4];
   end      

   ind = [3 5 9 11];
   if ie(6)==0 % bernoulli
     k(ind,ind) = (pe(1,1)*ie(1,3)/l)* ...
	[12/l2   -6/l   -12/l2  -6/l;-6/l   4   6/l  2;
	 -12/l2   6/l    12/l2   6/l;-6/l   2   6/l  4];
   else % timoshenko with reduced integration
      k(ind,ind) = ...
        pe(1,1)*ie(1,3)/l*[0 0 0 0;0 1 0 -1;0 0 0 0;0 -1 0 1] + ... % flexion
        pe(4)*ie(4)*ie(6)/l*[1 -l/2 -1 -l/2;-l/2 l2/4 l/2 l2/4;         % shear
                           -1 l/2 1 l/2;-l/2 l2/4 l/2 l2/4];
   end

   % pin flag handling (condense DOF a put zero stiffness)
   if size(elt,2)>9
    i1=abs(sprintf('%i',elt(1,9)))-48; i2=abs(sprintf('%i',elt(1,10)))-48;
    i2=[i1(i1~=0) i2(i2~=0)+6];
    if ~isempty(i2)
      i1=i2(i2~=0); i2=1:12;i2(i1)=0;i2=find(i2);
      Tpin([i2(:);i1(:)],1:length(i2))=[eye(length(i2),length(i2));-k(i1,i1)\k(i1,i2)];
      k1=zeros(12,12); k1(i2,i2)=Tpin'*k*Tpin; k=k1;
    end
   end

   k=of_mk('xkx_trans',x,k);  if off; k = tr'*k*tr;end % coordinate transformation
   out=k; out1=[];
end

if any([0 2]==typ) 
%% #mass - - - - - - - - - - - - - - - - - - - - -

   m = zeros(12,12); r2=(pe(3) * l); % RhoL
   ind = [1 7];    m(ind,ind) = [1 .5;.5 1]*(ie(4)*r2/3);
   ind = [4 10];   m(ind,ind) = [1 .5;.5 1]*(ie(2)+ie(3))*r2/3; % correct 23/4/98
   if ie(2); 
    c1= r2*ie(4)/210;  % rho A L/210
    c2= pe(3)*ie(2)/30/l; if Beam1In==0;c2=0;end
    ind = [2 6 8 12]; m(ind,ind) =  c1 * ... % rho A L/210
	[78    11*l   27   -6.5*l; 11*l    2*l2    6.5*l -1.5*l2
	 27    6.5*l  78   -11*l ; -6.5*l -1.5*l2 -11*l   2*l2]+ ...
     c2 * ... % secondary inertia effect, Lalanne p. 92
     [36 3*l -36 3*l;3*l 4*l2 -3*l -l2;-36 -3*l 36 -3*l;3*l -l^2 -3*l 4*l^2];
   end
   if ie(3); 
    c2= pe(3)*ie(3)/30/l; if Beam1In==0;c2=0;end       
    ind = [3 5 9 11]; m(ind,ind) =  r2*ie(4)/210 * ...
	[78   -11*l   27    6.5*l;-11*l    2*l2   -6.5*l -1.5*l2
	 27   -6.5*l  78    11*l ; 6.5*l  -1.5*l2  11*l   2*l2]+ ...
     c2* ... % secondary inertia effect, Lalanne p. 92
     [36 -3*l -36 -3*l;-3*l 4*l2 3*l -l2;-36 3*l 36 3*l;-3*l -l^2 3*l 4*l^2];
   end
   if ie(8) % Handle non structural mass
       ind = [1 7];    m(ind,ind) = m(ind,ind)+[1 .5;.5 1]/3*ie(8)*l;
       ind = [2 8];    m(ind,ind) = m(ind,ind)+[1 .5;.5 1]/3*ie(8)*l;
       ind = [3 9];    m(ind,ind) = m(ind,ind)+[1 .5;.5 1]/3*ie(8)*l;
   end
   if ~isempty(Tpin);  m1=zeros(12,12); m1(i2,i2)=Tpin'*m*Tpin; m=m1; end

   if ie(7)==0 % consistent mass
     m=of_mk('xkx_trans',x,m); % coordinate transformation
     if off; m = tr'*m*tr;end % coordinate transform   
   elseif ie(7)==2 % simple lump mass : half mass at node
    r1=[.5 .5 .5 0 0 0 .5 .5 .5 0 0 0]*(ie(4)*r2+ie(8)*l);    m=diag(r1);
   % lumped mass (sum over elements)
   else %if ie(7)==1%length(constit)>=point(7)+13 &&constit(point(7)+13)==1
    m=of_mk('xkx_trans',x,m); if off; m = tr'*m*tr;end % coordinate transform
    %r2=[0 0 0 1 0 0   0 0 0 1 0 0;
    %    0 0 0 0 1 0  0 0 -L 0 1 0;0 0 0 0 0 1   0 L 0 0 0 1];
    r2=[zeros(3);x';x'*[0 0 0;0 0 l;0 -l 0];x'];
    r3=sum(r2'*m*r2)/2;if any(r3<0); r3=diag(r2'*m*r2)'/2;end
    r1=[(m(1)+m(7))*[1 1 1] r3 (m(1)+m(7))*[1 1 1] r3]; 
    m=diag(r1);
   end
   if typ==2; out=m; out1=[];  
   else;out1=m; end
elseif typ==70
   out=zeros(12);out1=[]; 
   if 1==2  % this is used for verification purposes
    rule=integrules('beam1',3); % dof 2 6 8 12
    %r1=diag([1 -1 1 -1]);rule.N=rule.N*r1;rule.Nr=rule.Nr*r1;rule.Nrr=rule.Nrr*r1;
    k=zeros(4);%integrules('buildndn',13,rule,[0 0 0 1;1 0 0 1]);
    rule.jdet(:)=1/l;
    for j1=1:size(rule.w,1)
     % 2 rho I_polar (sum Iy Iz)  Ny,x  Nz,x  expression for non rotating frame
     k=k+rule.Nr(j1,:)'*rule.Nr(j1,:)*rule.jdet(j1)*rule.w(j1,4);
    end
   end
   % See Lalanne & Ferraris second edition page 92 (they rotate around y)
   C1e = [   0      -36     3*l       0      0     36    3*l       0;
        36        0       0     3*l    -36      0      0     3*l;
      -3*l        0       0  -4*l*l    3*l      0      0     l*l;
         0     -3*l   4*l*l       0      0    3*l   -l*l       0;
         0       36    -3*l       0      0    -36   -3*l       0;
       -36        0       0    -3*l     36      0      0    -3*l;
      -3*l        0       0     l*l    3*l      0      0  -4*l*l;
         0     -3*l    -l*l       0      0    3*l  4*l*l       0]';
  out([2 3 5 6 8 9 11 12],[2 3 5 6 8 9 11 12])=(pe(3)*mean(ie(2:3))/15/l)*C1e;
  out=of_mk('xkx_trans',x,out); % coordinate transformation
end

if ismember(typ,[3 6 8 9 7]); out=zeros(12,12); out1=[]; end 

if any(100==typ) % volumic load - - - - - - - - - - - - - - - - - -

  dire=varargin{7};x=x'; % element coordinates as colums

  if length(ie)<8; ie(8)=0; end

  FEw = [-sqrt(3/5) 0 0 .555555555555556;
                   0 0 0 .888888888888889;
           sqrt(3/5) 0 0 .555555555555556];
  %FEw=linspace(-1,1,50)';
  xi = [-1 0 0;1 0 0];
  FEna = [(1+FEw(:,1)*xi(:,1)')/2 ...         % linear shape functions
         .25*[2-3*FEw(:,1)+FEw(:,1).^3 ...  % cubic shape functions
               (1-FEw(:,1)-FEw(:,1).^2+FEw(:,1).^3) ...
               2+3*FEw(:,1)-FEw(:,1).^3 ...
               (-1-FEw(:,1)+FEw(:,1).^2+FEw(:,1).^3)]];
  dxds = ones(size(FEw,1),1)*xi(:,1)'/2*[0;l];
  FEna(:,[4 6])=FEna(:,[4 6])./dxds(:,[1 1]); % rotation dofs

  %plot(FEw,FEna)
  r1=x(:,1)*x(:,1)';  r2=x(:,2)*x(:,2)';  r3=x(:,3)*x(:,3)';

  b1=zeros(12,1);
  for j1=1:size(FEw,1) % loop on integration points
    d=zeros(12,3);
    d(1:3,:)=r1*FEna(j1,1);% node 1 translation
    d(7:9,:)=r1*FEna(j1,2);% node 2 translation
    if ie(5)==0 
       d(1:3,:)=d(1:3,:)+r2*FEna(j1,3); d(4:6,:)  =r2*FEna(j1,4);
       d(7:9,:)=d(7:9,:)+r2*FEna(j1,5); d(10:12,:)=r2*FEna(j1,6);
    else 
       d(1:3,:)=d(1:3,:)+r2*FEna(j1,1); 
       d(7:9,:)=d(7:9,:)+r2*FEna(j1,2);
    end
    if ie(6)==0 
       d(1:3,:)=d(1:3,:)+r3*FEna(j1,3); d(4:6,:)  =r3*FEna(j1,4);
       d(7:9,:)=d(7:9,:)+r3*FEna(j1,5); d(10:12,:)=r3*FEna(j1,6);
    else 
       d(1:3,:)=d(1:3,:)+r3*FEna(j1,1); 
       d(7:9,:)=d(7:9,:)+r3*FEna(j1,2);
    end
    F = dire*FEna(j1,1:2)';
    b1=b1+d*(F*dxds(j1)*ie(4)*FEw(j1,4));
  end % loop on integration points

  out=b1;

end
