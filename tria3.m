function [out,out1,out2]=tria3(CAM,varargin);

%TRIA3	element function for 3-node 15/18-DOF plate/shell element (dktp)
%
%	As all element functions (see ELEM0), TRIA3 is called by FE_MK for
%	model assembly, FEPLOT for structural deformation visualization, ...
%
%	TRIA3 header rows follow the format [Inf  abs('tria3') 0 EGID ...]
%	TRIA3 element property rows follow the format
%	    [n1 n2 n3 MatId ProId (EltId Theta Zoff T1 T2 T3)]
%         with
%          n1 ... n3  identification numbers for the element nodes
%	   MatId  material property identification number
%	   ProId  element property identification number
%	   EltId  optional element identification number
%	   Theta  unused
%          Zoff   optional offset from the surface of the nodes to the 
%		  reference plane
%	   T1 ... T3  optional thickness at nodes (used instead of IL entry)
%		  Currently the mean of Ti is used.
%
%       PL material property matrix. Isotropic materials [MatId 1  E  nu rho G]
%	   (See FE_MAT) are the only supported.
%
%       IL element property matrix. Rows associated to TRIA3 elements follow
%	   the format
%	    [Id# Type f   d 0 h   k 0 12I/T^3 0 NSM]
%          with
%	     Type : 1 for standard plate definition (no other type supported)
%	     f    : formulation (0 : thin plate (default))
%		    T3 for membrane, DKT for flexion
%	     d    : -1 no drilling stiffness, If d==0 d is set to 1.
%		    d>0 drilling DOF stiffness coefficient
%	     h    : thickness
%            k    : (unused by tria3) shear factor for thick plate model
%            12I/T^3 : Ratio of bending moment of inertia to nominal T^3/12
%                      (default 1)
%            NSM  : non structural mass per unit area
%
%	See sdtweb    quad4, eltfun, elem0
%	See also help quad4, quadb


%	Etienne Balmes
%       Copyright (c) 2001-2019 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

%#ok<*NASGU,*ASGLU,*NOSEM>

if comstr(CAM,'cvs')
 out='$Revision: 1.54 $  $Date: 2019/10/25 16:05:23 $'; return;
end
% standard calls with one input argument
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1);
 % Build Constit, Integ, and Elmap for later integration
 if comstr(Cam,'integinfo')

    %constit integ,elmap                 ID,pl,il
  [out,out1,out2]= ...
   p_shell('buildconstit',[varargin{1};18;3],varargin{2:end});

 elseif comstr(Cam,'matcall')||comstr(Cam,'call');  %Matrix assembly in fe_mknl
  [out,out1]=elem0('calltria3',varargin{:}); 
 elseif comstr(Cam,'rhscall') % call for load assembly
   out='be=tria3(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,elmap,EltConst,defe);';

 elseif comstr(Cam,'groupinit');out=elem0('groupinitogShell','tria3');
 elseif comstr(Cam,'constants'); error('was moved to elem0');
 elseif  comstr(Cam,'node');    out = 1:3;
 elseif  comstr(Cam,'prop');    out = [4 5 6];
 elseif comstr(Cam,'dofcall');
     out=elem0('dofcall'); % variable field elements
 elseif  comstr(Cam,'dof');
   out = [1+(1:6)'/100;2+(1:6)'/100;3+(1:6)'/100];
 elseif  comstr(Cam,'edge');    out = [1 2;2 3;3 1]; return; 
 elseif  comstr(Cam,'face');    out = [1 2 3];
 elseif  comstr(Cam,'flip');   out=[1 2 3];out1=[1 3 2]; 
 elseif  comstr(Cam,'sci_face');  out = [1 2 3];
 elseif  comstr(Cam,'line');    out = [1 2 3 1]; return; 
 elseif  comstr(Cam,'patch');   out = [1 2 3];return; 
 elseif  comstr(Cam,'parent');  out = 'tria3';

 % Basic matrix test of the element - - - - - - - - - - - - - - - - - -
 elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,5); 

  model=femesh('testtria3');
  [m,k,mdof]=fe_mknl(model);
  out=stack_cell(k,m);
  disp('TestMat passed');

 elseif  comstr(Cam,'test');   [CAM,Cam] = comstr(Cam,5);

   [out,out1]=femesh(strcat(['teststruct tria3' Cam]));
 
 % Basic tests of the element - - - - - - - - - - - - - - - - - - - - - -
 else;sdtw('''%s'' unknown',CAM);
 end

return
end % of standard calls with one input argument

% element matrix assembly - - - - - - - - - - - - - - - - - - - - - - - - - -
node=CAM; 
elt=varargin{1}; 
point=varargin{2};
integ=varargin{3}; 
constit=varargin{4};constit=constit(double(point(7))+(1:size(constit,1)))';
elmap=varargin{5};
EltConst=varargin{7};
typ=point(5);
if isa(elmap,'int32'); elmap = double(elmap); end

k = zeros(18,18);m=k; kb=zeros(9,9);

if size(node,1)==3; z = node(:,1:3);
else % Obsolete node
 NNode=sparse(node(:,1),1,1:size(node,1));
 i1 = NNode(elt(1:3)); z = node(i1,5:7);
end


if length(elt)>10  
 if any(elt(9:11)); t = mean(elt(9:11)); else;t=[]; end
else;t=[];
end

  db=constit([30 36 42;31 37 43;32 38 44]);
  dm=constit([9 15 21;10 16 22;11 17 23]);
  ds=constit([45 47;46 48]);

if size(node,2)==10 % xyz id v1xyz v3xyz (true normal must be kept)
 z = node(:,1:3);p=sp_util('basis',z(2,:)-z(1,:),z(3,:)-z(1,:));
 p=sp_util('basis',mean(node(:,5:7)),p(:,3));%mean(node(:,8:10)));
 p=p*[1 0 0;0 0 1;0 -1 0];

else;% Orient along line
 p=sp_util('basis',z(2,:)-z(1,:),z(3,:)-z(1,:));
end
x =z*p(:,1); y = z*p(:,2);


if any(typ==[0 1 2])
% ---------------------------------------------------------------------------- 
% computation of the membrane stiffness properties

EltConst.nodeE=[x y y*[0 0]];of_mk('buildndn',2,EltConst)

Bm=zeros(3,6); 
Bm([1 3 5 6 7 9 11 12 13 15 17 18])= ...
  EltConst.NDN([4 7 7 4 5 8 8 5 6 9 9 6]);
k([1 2 7 8 13 14],[1 2 7 8 13 14])=Bm'*(dm*(EltConst.jdet/2)*Bm);

if 1==2  % strain matrix bb in u/v coordinates. OBSOLETE
 J = det([ones(3,1) x y]); area=J/2;	% triangle area
 b = [y(2)-y(3) y(3)-y(1) y(1)-y(2)]; %Nx*jDet
 c = [x(3)-x(2) x(1)-x(3) x(2)-x(1)]; %Ny*jDet
 Bm1 = [b(1)   0    b(2)   0   b(3)   0
        0    c(1)     0  c(2)    0  c(3)
        c(1) b(1)   c(2) b(2)  c(3) b(3)];
 % membrane stiffness matrix
 k([1 2 7 8 13 14],[1 2 7 8 13 14]) = Bm1'*dm/(4*area)*Bm1; 
end

r1=constit(1); %r1 = pl(5)*il(6); % rho * h
%xxx if length(il)>11 & il(11); r1=r1+il(11);end
mm = (r1*EltConst.jdet/18)*ones(3,3);
m(1:6:18,1:6:18) = mm;    m(2:6:18,2:6:18) = mm;

% ----------------------------------------------------------------------------
% computation of the bending stiffness properties

if 1==1
 % triangle area
 J=EltConst.jdet(1);
 b = [y(2)-y(3) y(3)-y(1) y(1)-y(2)]; %Nx*jDet
 c = [x(3)-x(2) x(1)-x(3) x(2)-x(1)]; %Ny*jDet
 w=[1/6 1/6 1/6]';
 r=[1/6 2/3 1/6]';
 s=[1/6 1/6 2/3]';
 N = [2/3 1/6 1/6;1/6 2/3 1/6;1/6 1/6 2/3];
 Pk = [4 1 4;4 4 1;1 4 4]/9;Pm=Pk(:,[3 1 2]);
 Pkr = 4*[N(:,1)-r s -s]; Pks = 4*[-r r N(:,1)-s];
 Pmr = Pkr(:,[3 1 2]); Pms=Pks(:,[3 1 2]);
 Nr=[-1 1 0];Nr=Nr(ones(3,1),:);
 Ns=[-1 0 1];Ns=Ns(ones(3,1),:);

 Lk=[sqrt(c(3)^2+b(3)^2) sqrt(c(1)^2+b(1)^2) sqrt(c(2)^2+b(2)^2)];
 Lm=[Lk(3) Lk(1) Lk(2)];
 Ck=[c(3)/Lk(1) c(1)/Lk(2) c(2)/Lk(3)];
 Cm=[c(2)/Lk(3) c(3)/Lk(1) c(1)/Lk(2)];
 Sk=-[b(3)/Lk(1) b(1)/Lk(2) b(2)/Lk(3)];
 Sm=-[b(2)/Lk(3) b(3)/Lk(1) b(1)/Lk(2)];


 Lk=Lk(ones(3,1),:);Lm=Lm(ones(3,1),:);
 Ck=Ck(ones(3,1),:);Cm=Cm(ones(3,1),:);
 Sk=Sk(ones(3,1),:);Sm=Sm(ones(3,1),:);

  
  Nxr = [ 1.5*(Ck.*Pkr./Lk-Cm.*Pmr./Lm) ...
           0.75*(Ck.*Sk.*Pkr+Cm.*Sm.*Pmr) ...
           Nr-0.75*((Ck.^2).*Pkr+(Cm.^2).*Pmr) ];

  Nyr = [ 1.5*(Sk.*Pkr./Lk-Sm.*Pmr./Lm) ...
           -Nr+0.75*((Sk.^2).*Pkr+(Sm.^2).*Pmr) ...
           -Nxr(:,4:6)];

  Nxs = [ 1.5*(Ck.*Pks./Lk-Cm.*Pms./Lm) ...
           0.75*(Ck.*Sk.*Pks+Cm.*Sm.*Pms) ...
           Ns-0.75*((Ck.^2).*Pks+(Cm.^2).*Pms) ];

  Nys = [ 1.5*(Sk.*Pks./Lk-Sm.*Pms./Lm) ...
           -Ns+0.75*((Sk.^2).*Pks+(Sm.^2).*Pms) ...
           -Nxs(:,4:6)  ];
 for j1=1:3

  %BB=zeros(3,9); % DOFs 1z,2z,3z,1tx, ... 1ty
  %BB(1,3:3:end)=-EltConst.NDN(:,4+j1)'; % -N,x
  %BB(2,2:3:end)=EltConst.NDN(:,4+j1)';  %  N,y
  %BB(3,2:3:end)=EltConst.NDN(:,4+j1)';  %  N,y
  
  Bb = ... 
     [b(2)*Nxr(j1,:)+b(3)*Nxs(j1,:); c(2)*Nyr(j1,:)+c(3)*Nys(j1,:);
     c(2)*Nxr(j1,:)+c(3)*Nxs(j1,:)+b(2)*Nyr(j1,:)+b(3)*Nys(j1,:)]/J;

  kb=kb+Bb'*(w(j1)*J*db)*Bb;
  
 end
  ind = [3 9 15 4 10 16 5 11 17];
  Nwxy = [N (Pk.*Lk.*Sk-Pm.*Lm.*Sm)/8 -(Pk.*Lk.*Ck-Pm.*Lm.*Cm)/8];
  r1=constit(1);% rho * h + NSM added in sdtweb p_shell('buildconstit') 
  %xxx if length(il)>11 & il(11) r1=r1+il(11);end
  m(ind,ind)=(w(1)*r1*J)*(Nwxy'*Nwxy);
  k(ind,ind) = kb;
else
 
 % constit for DKTP anisotropi [rho*h eta tri(db)]
 ind = [3 9 15 4 10 16 5 11 17];
 [k1,m1]=of_mk('dktp',int32([45 45 0 0 0 0 0 0 0]), ...
   int32([0 0 9 3 0 0 3 3 3]),[constit(1);0;db([1;4;5;7;8;9])],[x y y*0]);
 elmapb=[ ...
     1     2     4     7    11    16    22    29    37
     2     3     5     8    12    17    23    30    38
     4     5     6     9    13    18    24    31    39
     7     8     9    10    14    19    25    32    40
    11    12    13    14    15    20    26    33    41
    16    17    18    19    20    21    27    34    42
    22    23    24    25    26    27    28    35    43
    29    30    31    32    33    34    35    36    44
    37    38    39    40    41    42    43    44    45
];

  k(ind,ind)=k1(elmapb);m(ind,ind)=m1(elmapb);

end % 1==2

% ----------------------------------------------------------------------------
% assembly and coordinate transformation of stiffness matrix

  if constit(4)~=-1
      
   %r1=integrules('tria3',[0 0 0],[x y y*0]);
   r2=zeros(3,18);
   r2([16 35 54])=-1;
   r2([4 5 6])=EltConst.NDN(4);
   r2([22 23 24])=EltConst.NDN(5);
   r2([40 41 42])=EltConst.NDN(6);

   %r2(1,[6 2 8 14])=[-1 r1.Nx];
   %r2(2,[12 2 8 14])=[-1 r1.Nx];
   % r2(3,[18 2 8 14])=[-1 r1.Nx];

   if constit(4)==0; constit(4)=1; end
   %i1=diag(k);kc=constit(4)*1e-6*mean(i1([4 5 10 11 16 17]));
   kc=constit(4)*1e-6*mean(k([58 77 172 191 286 305]));
   k=k+((r2'*kc)*r2);
   %def=feutilb('geomrb',[node(:,[4 4 4 4]) x y y*0]);svd(def.def'*k*def.def)

  end

  k=of_mk('xkx_trans',p',k); m=of_mk('xkx_trans',p',m);
  %d2=feutilb('geomrb',[node(:,[4 4 4 4]) node(:,1:3)]);feutilb('tkt',d2.def,k)
  if size(elt,2)>7
   if elt(8)~=0
    r = p(:,3)*elt(8);
    r = [eye(3,3) [0 r(3) -r(2);-r(3) 0 r(1);r(2) -r(1) 0];zeros(3,3) eye(3,3)];
    tr = zeros(size(k,1),size(k,2));for j1=6:6:size(k,1); tr(j1+(-5:0),j1+(-5:0))=r;end
    k=tr'*k*tr; m=tr'*m*tr;
   end
  end
end % of matrix assembly

out1=[];
if     typ==0;   out=k;out1=m;
elseif any(typ==[1 5]);  out=k;
elseif typ==2;  out=m; 
elseif typ==3;  out=zeros(18); % no viscous damping
elseif typ==100 % volume load

 EltConst=varargin{6};
 dirn=varargin{7}'; out=zeros(18,1);
 EltConst=evalin('caller','EltConst');
 if isempty(EltConst)%||size(EltConst.N,1)<3
   EltConst=integrules('tria3',-3); EltConst.VectMap=int32(reshape(1:18,6,3)');
   assignin('caller','EltConst',EltConst);
 end
 EltConst.nodeE=node;of_mk('buildndn',23,EltConst); b1=zeros(3,3);
 for jw=1:size(EltConst.w,1) % loop on integration points
    F=EltConst.N(jw,:)*dirn;
    b1=b1+EltConst.jdet(jw)*EltConst.w(jw,4)* ...
      constit(5)*EltConst.N(jw,:)'*F;
 end % loop on integration points
 out(EltConst.VectMap(:,1:3))=b1(:);

else;error('Not a supported element matrix type');
end



