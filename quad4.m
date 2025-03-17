function [out,out1,out2]=quad4(CAM,varargin);

%QUAD4	element function of a 4-node 20/24-DOF plate/shell element
%
%	As all element functions (see ELEM0), QUAD4 is called by FE_MK for
%	model assembly, FEPLOT for structural deformation visualization, ...
%
%	QUAD4 header rows follow the format [Inf  abs('quad4') 0 EGID ...]
%       QUAD4 element property rows follow the format
%           [n1 n2 n3 n4 MatId ProId EltId (Theta Zoff T1 T2 T3 T4)]
%         with
%          n1 ... n4  identification numbers for the element nodes
%                 if n4 == 0 or n4==n3 a TRIA3 is used instead of a quad4
%          MatId  material property identification number
%          ProId  element property identification number
%          Theta  angle of material x axis with n1-n2 line (unused)
%          Zoff   optional off-set of the element reference plane with 
%                 respect to the element nodes in the positive z-elt direction 
%          T1 ... T4 optional thickness at nodes. The mean of the 4 values is
%                 currently used
%
%       PL material property matrix. Isotropic material (see m_elastic)
%
%       IL element property matrix, see p_shell. 
%       Supported formulations (il(3)) for quad4 :
%          0 (default) : preferred formulation (currently 2)
%          1 : 4 tria3 thin plate elements with condensation
%	       of central node.
%		- If the last node is repeated a single tria3 is used
%		  you can thus mix 3 and 4 node plate/shells.
%		- for non flat shells, you should not set d to -1 or be
%		  prepared to face numerical conditioning problems
%	   2 : Q4WT for membrane and Q4gamma (MITC4) for bending
%              For non-flat quads formulation 1 is always used
%		- If the last node is repeated a single tria3 is used
%		  you can thus mix 3 and 4 node plate/shells.
%		- for non flat shells, you should not set d to -1 or be
%		  prepared to face numerical conditioning problems
%          4 : calls MITC4
%          5 : calls Q4CS (SDT element for composite shell)
%
%       Standard tests available with quad4('testeig') (mat,eig,load) 
%
%	See sdtweb     quad4, eltfun, elem0
%	See also help  quadb, tria3

%	Etienne Balmes
%       Copyright (c) 2001-2025 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       For revision information use quad4('cvs')

persistent TriaEltConst
%#ok<*NOSEM,*NASGU,*ASGLU>

% standard calls with one input argument
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1);
 % Build Constit, Integ, and Elmap for later integration
 if comstr(Cam,'integinfo')

  %constit integ,elmap                 ID,pl,il
  %See sdtweb p_shell('ShellConstit')
  [out,out1,out2]= ...
   p_shell('buildconstit',[varargin{1};24;4],varargin{2},varargin{3});
  %[pe,ie,dm,db,ds]=fe_mat(2,varargin{1},varargin{2},varargin{3});
  %if ie(3)==4; [out,out1,out2]=mitc4(CAM,varargin{:});% MITC4
  %elseif ie(3)==5; eval('[out,out1,out2]=q4cs(CAM,varargin{:});');
  %end
 % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 % Callback to define element constants during the fe_mknl init phase
 elseif comstr(Cam,'groupinit');out=elem0('groupinitogShell','quad4');
 elseif comstr(Cam,'constants'); error('was moved to elem0');
 elseif comstr(Cam,'matcall')||comstr(Cam,'call');  %Matrix assembly in fe_mknl
  [out,out1]=elem0('callquad4',varargin{:}); 
 % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'rhscall') % call for load assembly
   out='be=quad4(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,gstate,Case.GroupInfo{jGroup,8},defe);';

 % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'node'); out = 1:4;
 elseif comstr(Cam,'prop'); out = [5 6 7];
 elseif comstr(Cam,'dofcall');out=elem0('dofcall'); % variable field elements
 elseif comstr(Cam,'dof')
   out = 1:4;k=(1:6)'/100;out=out(ones(6,1),:)+k(:,ones(4,1));
   out=out(:);
 elseif comstr(Cam,'edge');   out = [1 2;2 3;3 4;4 1];
 elseif comstr(Cam,'face');   out = [1 2 3 4];
 elseif  comstr(Cam,'flip');  out = [1 2 3 4]; out1 = [1 4 3 2]; 
 elseif  comstr(Cam,'orders');  out = [1 2 3 4;2 3 4 1;3 4 1 2;4 1 2 3]; 
 elseif comstr(Cam,'line');   out = [1 2 3 4 1];
 elseif comstr(Cam,'patch');  out = [1 2 3 4];
 elseif comstr(Cam,'parent'); out = 'quad4';
 elseif comstr(Cam,'sci_face'); out = [1 2 3 4];

 % Basic matrix test of the element - - - - - - - - - - - - - - - - - -
 elseif  comstr(Cam,'testmat')

  model=femesh('testquad4');
  if ~isempty(strfind(Cam,'_1'));     model.il(:,3)=1;
  elseif ~isempty(strfind(Cam,'_2')); model.il(:,3)=2;
  elseif ~isempty(strfind(Cam,'_5')); model.il(:,3)=5;
  end
  if nargin==2; model.pl=varargin{1}; end
  [m,k,mdof]=fe_mknl(model); 
  out=stack_cell(k,m); if nargout==0;disp('TestMat passed');end

 elseif  comstr(Cam,'test'); [CAM,Cam] = comstr(CAM,5); 

   [out,out1]=femesh(strcat(['teststruct quad4' Cam]));


 % Basic tests of the element - - - - - - - - - - - - - - - - - - - - - -

 elseif comstr(Cam,'cvs');out='$Revision: 1.76 $  $Date: 2025/03/12 09:02:29 $';
 else; sdtw('''%s'' unknown',CAM); 
 end

return
end % of standard calls with one input argument

%% element matrix assembly - - - - - - - - - - - - - - - - - - - - - - - - - -
node=CAM; 
elt=varargin{1}; 
point=double(varargin{2});
integ=varargin{3};
% sdtweb p_shell('ShellConstit')
constit=varargin{4};constit=constit(point(7)+(1:size(constit,1))');
elmap=varargin{6};
typ=point(5);
if isa(elmap,'int32'); elmap = double(elmap); end

if size(integ,1)>6 %MITC4
 %constit [rho eta E nu T(1:4)]
  [k1,m1] = mitc4(CAM,varargin{:});
end

% find properties
if size(node,2)==7&&size(node,1)~=4
 NNode(node(:,1))=1:size(node,1); node=node(NNode(elt(1,1:4)),[5:7 1]);
end

if length(elt)>12 
    if any(elt(10:13)); t = elt(10:13);t=mean(t(t~=0)); 
    else;t = []; end
else;t=[];
end


% ---------------------------------------------------------------------------
if constit(5)==0; 
   sdtw('_nb','quad4 %i %i %i %i %i %i: thickness is equal to zero',elt(1:6));
   out=0;out1=0;return;

% ---------------------------------------------------------------------------
elseif typ>99 % loads, stresses, etc.

if any(typ==[100 102]) % 100 volume load, 102 inertial load - - - - - - - - -

i1=1:size(node,1);
if size(node,2)~=4; error('4 column nodE must be provided');end

 dirn=varargin{7}'; out=zeros(6*size(node,1),1);
 EltConst=evalin('caller','EltConst');
 if ~isempty(EltConst)&&size(EltConst.N,2)==size(node,1);
 elseif size(node,1)==4
   EltConst=integrules('quad4',-3);EltConst.VectMap=int32(reshape(1:24,6,4)');
   assignin('caller','EltConst',EltConst);
 else; error('Not an expected eltconst')
 end
 if size(dirn,1)~=size(EltConst.N,2);dirn=dirn(EltConst.VectMap);end
 EltConst.nodeE=node;of_mk('buildndn',23,EltConst); 
 b1=zeros(size(node,1),size(dirn,2));
 if typ==100; coef=constit(5);% h 
 elseif typ==102; coef=constit(1);% rho*h
 end
 for jw=1:size(EltConst.w,1) % loop on integration points
    F=EltConst.N(jw,:)*dirn; % dirn (Nnode,nfield)
    b1=b1+EltConst.jdet(jw)*EltConst.w(jw,4)* ...
      coef*EltConst.N(jw,:)'*F;
 end % loop on integration points
 out(EltConst.VectMap(:,1:size(b1,2)))=b1(:);
elseif typ==210 % membrane stress

 dirn=varargin{6};[bas,x]=basis(node(:,1:3)); 
 constit=[78000 0 2.2857e11  6.5141e10  2.2857e11 0 0  8.1712e10 0 0 0 0];
 integ=int32([0 0 8 4 0 0 1 3]);
 %constit=[78000 0 1.0 0.0 1.0 0.0 0.0 1.0];iopt=int32([0 0 8 4 200 0 1 1]);
 modeN=bas(:,1:2)'*dirn;
 out=of_mk('q4p',int32([3 0 0 0 200  0 0 0 0 0]),integ, ...
  constit,x,zeros(1,13),bas(:,1:2)'*dirn, ... % deformation at element
  zeros(10));

 %elmap=[1 3;3 2];idof(elmap)
 out1=bas;

else;error('Typ%i not implemented yet',typ);
end

% ---------------------------------------------------------------------------
elseif typ==3||typ==9; out=zeros(24); out1=[]; % Viscous matrix / unsym stiffness
% ---------------------------------------------------------------------------
elseif any(constit(3)==[0 1 2 3 4])
% See Hugues 132-135 : na = 1/4 (1+xa x) (1 + ya y)
% Q4WT for membrane see Batoz & Dhatt volume 1 pp 258-261
% Q4gamma (also known as MITC4) see Batoz & Dhatt volume 2 pp 326-330

% w   : [r s t weight] quadrature rule
% xi  : isoparametric coordinates of nodes
% na  : (integration points) x (shape functions)
% nar,s : (integration points) x (partial derivatives with respect to r,s)
if isempty(TriaEltConst); TriaEltConst=integrules('tria3',-1);end

if ~all(elt(1:4))||~all(diff(sort(elt(1,1:4)))) % actually a triangle

  i1=find(sparse(elt(1,1:4),1,1)); elt(1,1:3)=i1(:)';
  x = node(1:3,:); point(7)=0;
  [k,m] = tria3(x,elt([1:3 5:length(elt)]),point,integ,constit,...
    [],[],TriaEltConst,[]);
  x=x(:,1:3)-ones(size(x,1),1)*mean(x(:,1:3));
  bas = basis(x(2,:)/norm(x(2,:))+x(3,:)/norm(x(3,:)),x(3,:)); x=x*bas;
  k(24,24)=0;m(24,24)=0;
  if     typ==0;  out=k; out1=m;
  elseif typ==1;  out=k; out1=[];
  elseif typ==2;  out=m; out1=[];
  else;error('Not a supported element matrix type');
  end

else

[bas,x]=basis(node(:,1:3));

if size(integ,1)>5 % MITC4 call

elseif any(constit(3)==[1 3]) || abs(x(size(x,1),3))>max(max(abs(x)))*sqrt(eps)
%---------------------------------------------------------------------
% not a flat element use 4 tria3 elements
    % disp(sprintf('Warning not a flat element (%g)',x(end,3)))

  node(5,:)=[mean(node(:,1:3),1) max(node(:,4))+1];

  i3=[0 0 0;1 2 5;2 3 5;3 4 5;4 1 5]; m=zeros(30,30);k=m;
  for j1=2:5
    [k1,m1] = tria3(node(i3(j1,:),:),i3(j1,:),[0 0 0 0 0 0 0 0 0],integ,constit,elmap,[],TriaEltConst,[]);

    i4=(-5:0)';i4=i4(:,[1 1 1]);i5=i3(j1,:);i5=i5(ones(6,1),:)*6+i4;i5=i5(:);
    k(i5,i5)=k(i5,i5)+k1;
    m(i5,i5)=m(i5,i5)+m1;
  end 
  % condensation of centernode
  %def=feutilb('geomrb',node(:,[4 1 1 1 1 2 3]));svd(def.def'*k*def.def)

  if size(elt,2)>8 
      if elt(9)~=0 % z-offset
            r = bas(:,3)*elt(9);
            r = [eye(3,3) [0 r(3) -r(2);-r(3) 0 r(1);r(2) -r(1) 0];zeros(3,3) eye(3,3)];
            tr = zeros(size(k,1),size(k,2));
            for j1=6:6:size(k,1); tr(j1+(-5:0),j1+(-5:0))=r;end
            k=tr'*k*tr; m=tr'*m*tr;
      end
  end


  i1=diag(k);i1=find(i1>eps*mean(i1));
  if length(i1)==size(k,1)
   i1=1:24;i2=25:30;
   tr=[eye(length(i1),length(i1));-pinv(real(k(i2,i2)))*real(k(i2,i1))];
  else
   i2=zeros(30,1);i2(i1)=1;i1=find(i2(1:24));i2=24+find(i2(25:30));
   tr=zeros(30,24);
   tr([i1;i2],i1)=[eye(length(i1),length(i1));-real(k(i2,i2))\real(k(i2,i1))];
  end
  k=tr'*k*tr; m=tr'*m*tr;

  if     typ==0;  out=k; out1=m;
  elseif typ==1;  out=k; out1=[];
  elseif typ==2;  out=m; out1=[];
  else;error('Not a supported element matrix type');
  end



else % this is a flat element use Q4WT ----------------------------------------

EltConst=varargin{7};

if 1==2
 EltConst.nodeE=x;of_mk('buildndn',23,EltConst);

elseif 1==1 % Obsolete code kept for reference for a while
 %x=x.*(1+rand(size(x))/100);of_mk('buildndn',2,EltConst,x);
 % going back to local gradient information
  xi = [-1 -1 0;1 -1 0;1 1 0;-1 1 0];
  FEnar=EltConst.Nr(:,1:4); FEnas=EltConst.Ns(:,1:4); FEw=EltConst.w;
  FEna=EltConst.N(:,1:4);
  xr = FEnar*x(:,1); xs = FEnas*x(:,1);
  yr = FEnar*x(:,2); ys = FEnas*x(:,2);

  jdet = xr.*ys-xs.*yr;  jdet=jdet*sign(jdet(1));
  if any(jdet<0); 
   fprintf('\nquad4 (%i %i %i %i %i %i) negative Jacobian',elt(1:6));
  end
  %isequal(jdet,EltConst.jdet)
  d1=diag(1./jdet);
  i2=ones(size(FEnar,2),1); % number of shape functions
  nax = d1*[ ys(:,i2).*FEnar-yr(:,i2).*FEnas];
  nay = d1*[-xs(:,i2).*FEnar+xr(:,i2).*FEnas];

  nbx = -2*d1*[ ys.*FEw(:,1) -yr.*FEw(:,2)]; %additional quad shape for Q4WT
  nby = -2*d1*[-xs.*FEw(:,1)  xr.*FEw(:,2)];
  %NDN=[EltConst.N' [nax nbx]' [nay nby]'];norm(NDN-EltConst.NDN)
end % end of obsolete code

  k = spalloc(24,24,48);
 % stiffness matrix assembly - - - - - - - - - - - - - - - - - - - - - - - -

 if any([0 1]==typ);

  kbb = zeros(4,4);kab = zeros(24,4);
  BMb=zeros(3,4);idof=[elt(1,1)+.01];
  %BB=zeros(3,24); BM=zeros(3,24); BS=zeros(2,24);
  BM=spalloc(3,24,16);BB=spalloc(3,24,16); BS=spalloc(2,24,24);

  %          rho*h       eta       f d h k 12I/t3 nsm   dd(9:44) ds(
  db=constit([30 36 42;31 37 43;32 38 44]);
  dm=constit([9 15 21;10 16 22;11 17 23]);
  ds=constit([45 47;46 48]);

  %ija1 = reshape([-xi(1:4,1).*yr(6:9) -xi(1:4,2).*ys(6:9) ...
  %                xi(1:4,1).*xr(6:9) xi(1:4,2).*xs(6:9)]',2,2,4);

  % see Hugues p. 150 and 322
  for j1=1:4
     %BM(1,1:6:24)=          nax(j1,:)            ;
     %BM(2,2:6:24) =                    nay(j1,:) ;
     %BM(3,[1:6:24 2:6:24])=[nay(j1,:)  nax(j1,:)];
     %BB(1, 5:6:24) =                  -nax(j1,:) ;
     %BB(2, 4:6:24) =        nay(j1,:)            ;
     %BB(3,[4:6:24 5:6:24])=[nax(j1,:) -nay(j1,:)];
     %BS(1,3:6:24)=nax(j1,:); BS(2,3:6:24)=nay(j1,:);
     %BMb(1,[1 3])=          nbx(j1,:)            ;
     %BMb(2,[2 4]) =                    nby(j1,:) ;
     %BMb(3,[1 3 2 4])=[nby(j1,:)  nbx(j1,:)];
     r1=nax(j1,:);
       BM([1 19 37 55])=r1; BM([6 24 42 60])=r1;
       BB([12 30 48 66])=r1;BB([13 31 49 67]) = -r1;
       BS([5 17 29 41])=r1; 
     r1=nay(j1,:);
       BM([5 23 41 59])  =r1; BM([3 21 39 57])=r1;
       BB([11 29 47 65]) =r1; BB([15 33 51 69])=-r1;
       BS([6 18 30 42])=r1;
     r1=nbx(j1,:); BMb([1 7])=r1; BMb([6 12])=r1; 
     r1=nby(j1,:); BMb([5 11])=r1; BMb([3 9])=r1;
     % this is the MITC4 or (also called Q4gamma) model
     r1=FEnar(j1,1:4); r2=FEnas(j1,1:4);
     %ja=reshape([ys(j1)*r1;-xs(j1)*r1;-yr(j1)*r2;xr(j1)*r2],2,2,4);
     %ija1 = reshape([-xi(1:4,1).*yr(6:9) -xi(1:4,2).*ys(6:9) ...
     %             xi(1:4,1).*xr(6:9) xi(1:4,2).*xs(6:9)]',2,2,4);
     for j2 = 1:4
      r3 = [ys(j1)*r1(j2) -yr(j1)*r2(j2);-xs(j1)*r1(j2) xr(j1)*r2(j2)] * ...
         [-xi(j2,1)*yr(5+j2) xi(j2,1)*xr(5+j2)
          -xi(j2,2)*ys(5+j2) xi(j2,2)*xs(5+j2)];
      BS(j2*12-[5 4 3 2]) = r3(:)/jdet(j1);
     end
     % FEw 1 so it is omitted
     dm1 = dm*jdet(j1);
     k = k + BM'*(dm1)*BM + ...
           BB'*full((db*jdet(j1))*BB) + ...
           BS'*full((ds*jdet(j1))*BS);
     kbb = kbb + BMb'*(dm1)*BMb;
     kab = kab + BM' *(dm1)*BMb;

  end % loop on integration points

  k = k - kab*inv(kbb)*kab';% condensation of additional nodes linked to Q4WT

%  for j1=5%1:4 % for reduced integration use j1=5;
%     BS(1,[3:6:24 5:6:24]) = [nax(j1,:)/jdet(j1)             FEna(j1,:)];
%     BS(2,[3:6:24 4:6:24]) = [nay(j1,:)/jdet(j1) -FEna(j1,:)           ];
%     k = k + BS'*(ds*FEw(j1,4)*jdet(j1))*BS;
%  end % loop on integration points

  % Deal with drilling DOF
  % drilling DOF \theta_z-(v,x-u,y)/2
  if constit(4)~=-1 
   %i1=diag(k);if il(4)==0 il(4)=1; end
   i1=diag(k);if constit(4)==0; constit(4)=1; end
   kc=constit(49);
  %rb=feutilb('geomrb',[(1:4)'*[1 0 0 0] x],[0 0 0],feutil('getdof',(1:4)',(1:6)'/100));
   for j1=1:4
    DD=zeros(1,24);DD(6:6:end)=FEna(j1,:);
    DD(1:6:end)=+nay(j1,:)/2;DD(2:6:end)=-nax(j1,:)/2;
    k=k+(jdet(j1)*kc)*(DD'*DD);
   end

   % tests for drilling DOFs. The objective is to have zero energy
   % for the rigid body rotation around drilling DOF and non-zero but small
   % otherwise. One defines r1 with 4 motion of corner nodes associated with
   % unit drilling rotations.
   %
   % rb=zeros(24,1);rb([1 7 13 19 2 8 14 20 6 12 18 24])= ...
   %   [-x(:,2);x(:,1);[1;1;1;1]];
   
   end

 end;
 if any([0 2]==typ); % mass matrix assembly

  m = spalloc(24,24,48);  B=spalloc(3,24,12); BR = spalloc(2,24,12);
  r1=constit(1); % rho * h + NSM added in sdtweb p_shell('buildconstit')
  %xxx if length(il)>11 & il(11) r1=r1+il(11);end  
  for j1=1:4
     B(1,1:6:24)= FEna(j1,:);
     B(2,2:6:24)= FEna(j1,:);
     B(3,3:6:24)= FEna(j1,:);
     m = m + (r1*FEw(j1,4)*jdet(j1)*B')*B;
     BR(1,4:6:24)= FEna(j1,:);
     BR(2,5:6:24)= FEna(j1,:);
     %m = m + (pl(5)*il(6)^3/12*FEw(j1,4)*jdet(j1)*BR')*BR;
     m = m + (r1*constit(5)^2/12*FEw(j1,4)*jdet(j1)*BR')*BR;
  end % loop on integration points
  if typ==2; k=m; m=[];end

 else;m=[];
 end % mass matrix
 if ~isempty(k); k=of_mk('xkx_trans',bas',full(k)); end
 if ~isempty(m); m=of_mk('xkx_trans',bas',full(m)); end

 if size(elt,2)>8&&elt(9)~=0 % z-offset
        r = bas(:,3)*elt(9);
        r = [eye(3,3) [0 r(3) -r(2);-r(3) 0 r(1);r(2) -r(1) 0];zeros(3,3) eye(3,3)];
        tr = zeros(size(k,1),size(k,2));
        for j1=6:6:size(k,1); tr(j1+[-5:0],j1+[-5:0])=r;end
        k=tr'*k*tr; if ~isempty(m); m=tr'*m*tr;end
 end
 out=k;out1=m;

end % off choice between 1 quad and 4 triangles
end % off choice between quad and triangle

if typ==100
 error(1);
end

% ---------------------------------------------------------------------------
else;error('unknown formulation');
end

