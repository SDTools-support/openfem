function [out,out1,out2]=p_mitc(varargin)

%p_mitc testmitc4
%p_mitc testmitc9
% This file supports
%   MitcGen : generic assembly
%   p_mitc('command') ...
%   p_mitc(NumberOfNodes) builds constants
%   Standard commands of p_functions
%    'BuildDofOpt' : select 5 or 6 DOF depending on wether
%      drilling==-1 il(4) in p_mitc entry
%    'BuildConstit'
%      return integ=[MatId;ProId;0;0];
%             constit=[-1;type;E;nu]
%    'Const'
%      return the EltConst  the normal map building is done in p_shell
%      since this is a classical problem that should be true for all
%      shells.

%       Dominique Chapelle, Etienne Balmes
%       Copyright (c) 2001-2009 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license


if isstruct(varargin{1}) % Standard call for matrix computation
	%K=KmitcGen(node,const,EltConst);
 if isstruct(varargin{2});out=KmitcGen(varargin{:});
 else % prepare transformation from what is available in elem0/of_mk
   %ke=p_mitc(EltConst,nodeE,constit,DofPos);
   nodeE=varargin{2};constit=varargin{3};DofPos=varargin{4};
   r1=struct('pos',nodeE(:,1:3)','a3',nodeE(:,5:7)','t',nodeE(:,8)');
   r2=struct('E',constit(3),'nu',constit(4));
   ke=KmitcGen(r1,r2,varargin{1});
   if size(ke,1)<size(DofPos,1); %size(r1.pos,2)*6;
	 ke(end+1,end+1)=0;
	 ind=reshape(1:size(r1.pos,2)*5,5,[]);ind(6,:)=length(ke);
	 ke=ke(ind,ind);%figure(1);spy(ke);
   end
   out=ke;
 end

%ke=p_mitc(EltConst,nodeE,constit,DofPos);elmap=[];
% test of callback from of_mk for each element
elseif isa(varargin{1},'int32')
    
[pointers,integ,constit,gstate,elmap,InfoAtNode,EltConst,def,jElt,DofPos]=deal(varargin{:});
 jElt=jElt(1)+1;point=pointers(:,jElt);
 nodeE=EltConst.nodeE;Nw=EltConst.Nw;%defe=EltConst.defe;
 if isfield(EltConst,'pEC')&&size(nodeE,2)>8
   r1=struct('pos',nodeE(:,1:3)','a3',nodeE(:,EltConst.pEC(3)+int32(1:3))', ...
      't',nodeE(:,EltConst.pEC(3)+1)');
 else
   r1=struct('pos',nodeE(:,1:3)','a3',nodeE(:,5:7)','t',nodeE(:,8)');
 end
   r2=struct('E',constit(3),'nu',constit(4));
   ke=KmitcGen(r1,r2,EltConst);
   if size(ke,1)<size(DofPos,1); %size(r1.pos,2)*6;
	 ke(end+1,end+1)=0;
	 ind=reshape(1:size(r1.pos,2)*5,5,[]);ind(6,:)=length(ke);
	 ke=ke(ind,ind);%figure(1);spy(ke);
   end
   i1=int32(1:size(ke,1));
   sp_util('setinput',EltConst.ke,ke,i1,i1);
   return;
	
elseif ischar(varargin{1})

    
[CAM,Cam]=comstr(varargin{1},1);carg=2;

% --------------------------------------------------------------
% 'test' ...
if strncmpi(Cam,'test',4);[CAM,Cam]=comstr(Cam,5);

out=[];
% 'TestModel ElemF' - - - - - - - - - - - - - - - - - - - - -
if comstr(Cam,'model')||isempty(Cam);
    
 try;[CAM,Cam]=comstr(Cam,6); if isempty(Cam);Cam='quad4'; end
 catch;Cam='quad4'; end
 model=femesh(['test' Cam 'divide 3 1']);
 model.il=[110 fe_mat('p_mitc','SI',1) 4 -1 0 .01];
 model.il=p_solid(model.il,'dbval 99 d3 3'); % test multi prop
 [Case,model.DOF]=fe_mknl('init',model);
    
 Case.GroupInfo{2}(4,:)=0; % ndn rule, need to deal with Mass

 sp_util('diag',0) % utiliser of_mk
 k1=fe_mknl('assemble',model,Case,1);
 ind=fe_c(model.DOF,.06,'ind',2);k1=k1(ind,ind);

 
 
% TestMitc4 - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif any(Cam=='4')   
X=[0,.1,.1,0];[a,b]=constMitc4;
Y=[0,0,.1,.1];
node=struct('pos',[100.*X ; 100.*cos(pi/2*Y) ; 100.*sin(pi/2*Y)]);
node.a3=[zeros(1,4) ; cos(pi/2*Y) ; sin(pi/2*Y)];

node.t=1.*ones(1,4);

const.E=3.e7;const.nu=.3;

EltConst=p_mitc(4);
K=p_mitc(node,const,EltConst);

fname='K_mitc4';
if exist(fname,'file');load(fname)
kdiag=diag(K_mitc4);kref=sqrt(kdiag*kdiag');
disp(['Rel. Diff. w.r.t. ref: '  num2str(max(max(abs(K_mitc4-K)./kref))) ' (should be <1.e-14)']);
clear K
end

% TestMitc9 - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif any(Cam=='9')

X=[0,.1,.1,0,.05,.1,.05,0,.05];%[a,b]=constMitc4;
Y=[0,0,.1,.1,0,.05,.1,.05,.05];
node=struct('pos',[100.*X ; 100.*cos(pi/2*Y) ; 100.*sin(pi/2*Y)]);
node.a3=[zeros(1,9) ; cos(pi/2*Y) ; sin(pi/2*Y)];
node.t=1.*ones(1,9);

const.E=3.e7; const.nu=.3;

EltConst=p_mitc(9);K=p_mitc(node,const,EltConst);

fname='K_mitc9';
if exist(fname,'file');load(fname)
kdiag=diag(K_mitc9);
kref=sqrt(kdiag*kdiag');
disp(['Rel. Diff. w.r.t. ref: '  num2str(max(max(abs(K_mitc9-K)./kref))) ' (should be <1.e-14)']);
clear K
end

if 1==1
	model=femesh('testquad9');
	model.Node=sortrows(model.Node);model.Node(:,5:7)=node.pos';
	model.pl=[100 fe_mat('m_elastic','SI',1) const.E const.nu];
	model.il=[110 fe_mat('p_mitc','SI',1) 4   -1 0 node.t];

	r1=struct('data',node.a3,'NodePos',int32([1:9])','lab',{{'nx','ny','nz','t'}});
	r1.data(4,:)=1; % thickness at node
    model=stack_set(model,'MAP','Shell_1',r1);
    
    [Case,model.DOF]=fe_mknl('init',model);

	sp_util('diag',12) % utiliser elem0
	k1=fe_mknl('assemble',model,Case,1);
	ind=fe_c(model.DOF,.06,'ind',2);k1=k1(ind,ind);
 if exist(fname,'file')
  disp(['Diff. w.r.t. ref: '  num2str(max(max(abs(K_mitc9-k1)))) ...
     ' (should be <1.e-7)']);
 end
end

    
% End Test Mitc9 - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'of_mk_lib')
    
hfile = [matlabroot '\extern\include\matrix.h'];
loadlibrary('libmx', hfile)
y = rand(4, 7, 2);
calllib('libmx', 'mxGetNumberOfElements', y)
calllib('libmx', 'mxGetClassID', y)
cdh
hfile = [matlabroot '\extern\include\matrix.h'];
loadlibrary('libmx', hfile, 'mfilename', 'mxproto')



unloadlibrary libmx



hfile = fullfile(fileparts(which('ofutil')),'mex','of_mk.c');

wd=fullfile(fileparts(which('ofutil')),'mex');cd(wd)

st='of_time';unloadlibrary(st);
delete([st '.' mexext]);mex -v -DOSTYPEmexw32 of_time.c
[a,b]=loadlibrary('of_time.mexw32','of_time.c','includepath',wd,'addheader','of_time_interp.c');
libisloaded(st)
m=libfunctions('of_time.mexw32','-full')
st='of_time';unloadlibrary(st);delete([st '.' mexext]) 

    
else; error('Not a known test');
end

% Build the Field per DOF - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'builddofopt');
    
RunOpt=varargin{carg};carg=carg+1;
pl=varargin{carg};carg=carg+1;
il=varargin{carg};carg=carg+1;
if il(4)==-1;RunOpt.FieldDofs=[1:5]; % no drilling  
else; RunOpt.FieldDofs=[1:6]; 
end
out=RunOpt;

% Build the constit vector - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'buildconstit');
  %[constit,ID]=feval(st,'buildconstit',ID,pl,il0);
  ID=varargin{carg};carg=carg+1;
  pl=varargin{carg};carg=carg+1;
  il=varargin{carg};carg=carg+1;
  pl=pl(pl(:,1)==ID(1),:);
  il=il(il(:,1)==ID(2),:);
  constit=[-1;il(2);pl(3);pl(4);il(6)];
  out=constit; out1=ID; out2=[];
  
% Build the normal map ... - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'const');
 ElemF=sscanf(varargin{2},'%s',1);
 [st1,st,integ,constit,model,Case,cEGI,RunOpt]=deal(varargin{:});
 r2=constit;r2(1)=0;

 [EC,out1,InfoAtNode]=p_shell('Const',st,integ,r2,model,Case,cEGI,RunOpt);
 if 1 %size(integ,1)>3&&~any(integ(3:4,1))% back propagate number of DOFs if needed
     i1=5; if size(integ,1)<4; error('Improper init');end% dof per elt
     sp_util('setinput',integ, ...
         [EC.Nnode*i1;EC.Nnode],int32(3:4),int32(1));
     EC.ElMap=elem0('elmapmat_og',[EC.Nnode i1]);
 end

 
 if size(InfoAtNode.data,1)<4; % make sure thickness is in the last dim.
     % only works with single dim
     InfoAtNode.data(4,:)=constit(5);InfoAtNode.lab{4}='t';
 end

 r1=p_mitc(length(feval(ElemF,'node'))); % reads EC and returns in r1
 r1.jdet=zeros(r1.Nw,1);
 if 1==1 % adjustments for callback testing in elem0
  r1.material='callback';r1.fHandle=@p_mitc;
  r1.nodeE=zeros(EC.Nnode,8);EltConst.defe=zeros(EC.Nnode*6);
  r1.MatrixIntegrationRule=cell(1,2); 
  r1.VectMap=EC.VectMap;r1.ElMap=EC.ElMap;r1.ke=zeros(size(EC.ElMap));
 end
 
 r1.pEC=zeros(40,1,'int32'); % copy of EC C structure for callbacks
 out=r1;  
 
 out1=23; % NDN rule : xxx need to use a 23 for mass
 out2=InfoAtNode;

elseif comstr(Cam,'database');% #DataBase
    out=[];out1=[];
elseif comstr(Cam,'cvs');
  out='$Revision: 1.9 $  $Date: 2012/03/16 16:02:01 $'; return;
else; error('%s unknow',CAM);
end
    
    
else% Call for constant array computation (only at init) fills EltConst

    
    
    EltConst.Nnode=varargin{1};
	EltConst.DofLabels={};
	EltConst.material='KmitcGen';
    clear interp;
    fun=eval(sprintf('@constMitc%i;',EltConst.Nnode));
	[integr,interp]=feval(fun); 
	% This is what is really needed by the Matlab code
	
	EltConst.interp=interp;
	EltConst.integr=integr;
	
	% Preallocate, build tables for C use
	EltConst.Nw=integr.n(1);
	EltConst.NodeE=zeros(EltConst.Nnode,3+1+3*3);
	EltConst.NDN=zeros(EltConst.Nnode,EltConst.Nw*3); % xxx may need adjust
	% place the interp struct in something compatible with C
	st=fieldnames(interp);
	EltConst.TyingInt=[];
	EltConst.NTying=[];
	EltConst.TyingLam=[];
	for j1=1:length(st);
		r1=getfield(interp,st{j1});
		EltConst.NTying=[EltConst.NTying r1.n];
		if r1.n>0
			EltConst.TyingInt=[EltConst.TyingInt r1.int];
			EltConst.TyingLam=[r1.lam(:).val r1.lam(:).r r1.lam(:).s];
		end
	end
	% Place the integr in something compatible with z
	% EltConst=
	% .w [gauss rule first [r s 0 wi], z rule last [0 0 zi wzi]]
	%	NOTE: r,s not given for now (0 instead)
	% EltConst.Nw -> end of wi (i.e. last row)
	% .N (Ns * Nw)
	% .Nr (Ns * Nw)
	% .Ns (Ns * Nw)
	EltConst.w=zeros(sum(integr.n),4);
	EltConst.w(1:integr.n(1),4)=integr.wrs;
	EltConst.w(integr.n(1)+1:end,3)=integr.z;
	EltConst.w(integr.n(1)+1:end,4)=integr.wz;
	for j1=1:EltConst.Nw
		EltConst.N(:,j1)=integr.lam(j1).val;
		EltConst.Nr(:,j1)=integr.lam(j1).r;
		EltConst.Ns(:,j1)=integr.lam(j1).s;
	end
	
	out=EltConst;

end

end



function K=KmitcGen(node,const,EltConst)
%
% Generic MITC element stiffness computation
%
% Input:
%	node: nodal data structure
%		node.pos: node positions (3,Nnodes)
%		node.a3: normals (3,Nnodes)
%		node.t: thickness (1,Nnodes)
%	const: mechanical constants
%		const.E: Young
%		const.nu: Poisson
%	integr: integration point data (Q2)
%		integr.n: number of in-plane and transverse Gauss points (1,2)
%		integr.lam(1:Ngauss): shape functions at Gauss points (.val,.r,.s fields)
%		integr.wrs(1:Ngauss): (2D) Gauss weights
%		integr.z: through-thickness Gauss points coordinates
%		integr.wz: corresponding weights
%	interp: MITC tying points data
%		interp.xx.n: number of tying points for xx-strain
%		interp.xx.lam(1:Ntying): shape functions at tying points (.val,.r,.s fields)
%		interp.xx.int: tying shape functions at Gauss points (Ngauss,Ntying)
%
% Output: K element stiffness matrix
%	


dimK=size(node.pos,2)*5;
K=zeros(dimK,dimK); % preallocated in of_mk


integr=EltConst.integr; interp=EltConst.interp;
Nnode=size(node.pos,2); % name in of_mk

% Precompute some nodal arrays
node.t2a3=zeros(3,Nnode);
node.t2V1=zeros(3,Nnode);
node.t2V2=zeros(3,Nnode);

for ji=1:Nnode % Loop on nodes  % -> use basis
	node.t2a3(:,ji)=node.t(ji)/2*node.a3(:,ji);
	% On-the-fly computation of (V1,V2) vectors
	V1=cross(node.a3(:,ji),[1;0;0]);
	normV1=norm(V1);
	if normV1<1.e-3
		V1=cross(node.a3(:,ji),[0;1;0]);
		normV1=norm(V1);
	end
	V1=V1/normV1;
	V2=cross(node.a3(:,ji),V1);
	node.t2V1(:,ji)=node.t(ji)/2*V1;
	node.t2V2(:,ji)=node.t(ji)/2*V2;
end


% allocate all things here (just to figure out actual size!)
Errc=zeros(1,dimK);
Essc=zeros(1,dimK);
Ersc=zeros(1,dimK);
Erzc=zeros(1,dimK);
Eszc=zeros(1,dimK);
Ersz=zeros(5,dimK);

Err=zeros(interp.rr.n,dimK);
Ess=zeros(interp.ss.n,dimK);
Ers=zeros(interp.rs.n,dimK);
Erz=zeros(interp.rz.n,dimK);
Esz=zeros(interp.sz.n,dimK);

% Stiffness computation loops:
%    * Across thickness
%       - over tying points
%       - over in-plane integration points

% Loop across thickness
for ji=1:integr.n(2)
	zc=integr.z(ji);
	% Loop over reinterpolated strains
	% 1- Err
    for jj=1:interp.rr.n
			Err(jj,:)=shell_strain('rr',interp.rr.lam(jj),zc,node);
    end
	% 2- Ess
	for jj=1:interp.ss.n
			Ess(jj,:)=shell_strain('ss',interp.ss.lam(jj),zc,node);
	end
	% 3- Ers
    for jj=1:interp.rs.n
			Ers(jj,:)=shell_strain('rs',interp.rs.lam(jj),zc,node);
    end
	% 4- Erz
    for jj=1:interp.rz.n
			Erz(jj,:)=shell_strain('rz',interp.rz.lam(jj),zc,node);
    end
	% 5- Esz
	for jj=1:interp.sz.n
			Esz(jj,:)=shell_strain('sz',interp.sz.lam(jj),zc,node);
	end
	% Loop over in-plane integration points
	%
	for jj=1:integr.n(1)
		% Strain interpolation (or not -> else)
		if ~isempty(Err)
			Errc=interp.rr.int(jj,:)*Err;
		else
			Errc=shell_strain('rr',integr.lam(jj),zc,node);
		end
		if ~isempty(Ess)
			Essc=interp.ss.int(jj,:)*Ess;
		else
			Essc=shell_strain('ss',integr.lam(jj),zc,node);
        end
        if ~isempty(Ers)
			Ersc=interp.rs.int(jj,:)*Ers;
		else
			Ersc=shell_strain('rs',integr.lam(jj),zc,node);
        end
        if ~isempty(Erz)
			Erzc=interp.rz.int(jj,:)*Erz;
	    else
			Erzc=shell_strain('rz',integr.lam(jj),zc,node);
        end
        if ~isempty(Esz)
			Eszc=interp.sz.int(jj,:)*Esz;
	    else
			Eszc=shell_strain('sz',integr.lam(jj),zc,node);
        end
		Ersz=[Errc ; Essc ; Ersc ; Erzc ; Eszc];

		% Constitutive matrix
		%  1- geometric coefficients
		Phi_r=(node.pos+zc*node.t2a3)*integr.lam(jj).r;
		Phi_s=(node.pos+zc*node.t2a3)*integr.lam(jj).s;
		tc=node.t*integr.lam(jj).val; % current thickness value
		grr=Phi_r'*Phi_r; % covariant metric tensor
		gss=Phi_s'*Phi_s;
		grs=Phi_r'*Phi_s;
		sg=sqrt(grr*gss-grs^2); % square root of g
		Gcon=inv([grr grs;grs gss]);% (2x2) contravariant metric tensor
		gcrr=Gcon(1,1);
		gcss=Gcon(2,2);
		gcrs=Gcon(1,2);
		%  2- mechanical coefs
		a1=const.E/(2*(1+const.nu));
		a2=2*const.nu/(1-const.nu);
		%
		C=zeros(5,5);
		C(1,1)=a1*(2+a2)*gcrr^2;
		C(1,2)=a1*(2*gcrs^2+a2*gcrr*gcss);
		C(2,1)=C(1,2);
		C(2,2)=a1*(2+a2)*gcss^2;
		C(3,3)=a1*(gcrr*gcss+(1+a2)*gcrs^2);
		C(1,3)=a1*(2+a2)*gcrr*gcrs;
		C(3,1)=C(1,3);
		C(2,3)=a1*(2+a2)*gcss*gcrs;
		C(3,2)=C(2,3);
		C(4:5,4:5)=a1*4/tc^2*Gcon;

		K=K+Ersz'*C*Ersz*sg*tc*integr.wz(ji)*integr.wrs(jj)/2;
	end
end


end



function e = shell_strain(type,lam,z,node)
%
% Computation of shell strains nodal vectors
%
% Input:
%	type: 'xx' string for strain component
%	lam: shape function values at current point (.val,.r,.s fields)
%	z: transverse coordinate at current point (in [-1,1])
%	node: nodal data
%
% Output:
%	e: nodal vector for strain xx (1,Ndof)
%

Nnode=size(node.pos,2);e=zeros(1,Nnode*5);

switch type
	case 'rr'
		Phi_r= (node.pos+z*node.t2a3)*lam.r;
		for ji=1:Nnode % Loop on nodes
			e(5*ji-4:5*ji)=[lam.r(ji)*Phi_r' lam.r(ji)*z*node.t2V1(:,ji)'*Phi_r lam.r(ji)*z*node.t2V2(:,ji)'*Phi_r];
		end
	case 'ss'
		Phi_s= (node.pos+z*node.t2a3)*lam.s;
		for ji=1:Nnode % Loop on nodes
			e(5*ji-4:5*ji)=[lam.s(ji)*Phi_s' lam.s(ji)*z*node.t2V1(:,ji)'*Phi_s lam.s(ji)*z*node.t2V2(:,ji)'*Phi_s];
		end
	case 'rs'
		Phi_r= (node.pos+z*node.t2a3)*lam.r;
		Phi_s= (node.pos+z*node.t2a3)*lam.s;
		for ji=1:Nnode % Loop on nodes
			e(5*ji-4:5*ji)=[lam.r(ji)*Phi_s'+lam.s(ji)*Phi_r' ...
				lam.r(ji)*z*node.t2V1(:,ji)'*Phi_s+lam.s(ji)*z*node.t2V1(:,ji)'*Phi_r  ...
				lam.r(ji)*z*node.t2V2(:,ji)'*Phi_s+lam.s(ji)*z*node.t2V2(:,ji)'*Phi_r];
		end
	case 'rz'
		% Commented: untruncated expressions (w.r.t. z)
		%Phi_r= (node.pos+z*node.t2a3)*lam.r;
		Phi_r0= node.pos*lam.r;
		Phi_z=node.t2a3*lam.val;
		for ji=1:Nnode % Loop on nodes
			%e=[e lam.r(i)*Phi_z'  lam.r(i)*z*node.t2V1(:,i)'*Phi_z+lam.val(i)*node.t2V1(:,i)'*Phi_r ...
			%	lam.r(i)*z*node.t2V2(:,i)'*Phi_z+lam.val(i)*node.t2V2(:,i)'*Phi_r];
			e(5*ji-4:5*ji)=[lam.r(ji)*Phi_z'  lam.val(ji)*node.t2V1(:,ji)'*Phi_r0 lam.val(ji)*node.t2V2(:,ji)'*Phi_r0];
		end
	case 'sz'
		% Commented: untruncated expressions (w.r.t. z)
		%Phi_s= (node.pos+z*node.t2a3)*lam.s;
		Phi_s0= node.pos*lam.s;
		Phi_z= node.t2a3*lam.val;
		for ji=1:Nnode % Loop on nodes
			%e=[e lam.s(i)*Phi_z'  lam.s(i)*z*node.t2V1(:,i)'*Phi_z+lam.val(i)*node.t2V1(:,i)'*Phi_s ...
			%	lam.s(i)*z*node.t2V2(:,i)'*Phi_z+lam.val(i)*node.t2V2(:,i)'*Phi_s];
			e(5*ji-4:5*ji)=[lam.s(ji)*Phi_z'  lam.val(ji)*node.t2V1(:,ji)'*Phi_s0 lam.val(ji)*node.t2V2(:,ji)'*Phi_s0];
		end
	case 'qq'
		error('Strain qq not implemented yet');
	otherwise
		error('Unknown strain type %s',type);

end

end




function [integr,interp] = constMitc4

% Interpolation (tying) and integration points data for MITC4
%

ga2=.5773502691896258;
rga2=[-ga2,ga2];
wga2=[1,1];

R12=[0,0];
S12=[-1,1];
R21=[-1,1];
S21=[0,0];

interp.rr.n=0;
interp.ss.n=0;
interp.rs.n=0;
interp.rz.n=2;
interp.sz.n=2;
interp.qq.n=0;

r1=integrules('quad4',[R12(:) S12(:)*[1 0 0]]);
interp.rz.lam=repmat(struct('val',zeros(4,1),'r',zeros(4,1),'s',zeros(4,1)),1,2);
for jk=1:2 % loop 2x1
	interp.rz.lam(jk).val=r1.N(jk,:)';
	interp.rz.lam(jk).r=r1.Nr(jk,:)';
	interp.rz.lam(jk).s=r1.Ns(jk,:)';
end

r1=integrules('quad4',[R21(:) S21(:)*[1 0 0]]);
interp.sz.lam=repmat(struct('val',zeros(4,1),'r',zeros(4,1),'s',zeros(4,1)),1,2);
for jk=1:2 % loop 2x1
	interp.sz.lam(jk).val=r1.N(jk,:)';
	interp.sz.lam(jk).r=r1.Nr(jk,:)';
	interp.sz.lam(jk).s=r1.Ns(jk,:)';
end

jk=0;
for ji=1:2
	s=rga2(ji);
	for jj=1:2
		jk=jk+1;
		r=rga2(jj);
		tying=tying_mitc4(r,s);
		interp.rz.int(jk,:)=tying.rz;
		interp.sz.int(jk,:)=tying.sz;
	end
end


% Integration points definitions
%

integr.n=[4 2]; % Number of in-plane and transverse Gauss points

% In-plane
r1=integrules('quad4');
for jk=1:4
	integr.lam(jk).val=r1.N(jk,:)';
	integr.lam(jk).r=r1.Nr(jk,:)';
	integr.lam(jk).s=r1.Ns(jk,:)';
end
integr.wrs=r1.w(:,4);

% Transverse
integr.z=rga2;
integr.wz=wga2;

end


function [integr,interp] = constMitc9

% Interpolation (tying) and integration points data for MITC9
%

r0=integrules('quad9');
ga=-r0.w(1,1);% Basic coordinate for 3 point Gauss quadrature

ga2=.5773502691896258;
rga2=[-ga2,ga2];
wga2=[1,1];

R23=[-ga2,ga2,-ga2,ga2,-ga2,ga2];
S23=[-ga,-ga,0,0,ga,ga];
R22=[-ga2,ga2,-ga2,ga2];
S22=[-ga2,-ga2,ga2,ga2];
R32=[-ga,0,ga,-ga,0,ga];
S32=[-ga2,-ga2,-ga2,ga2,ga2,ga2];

interp.rr.n=6;
interp.ss.n=6;
interp.rs.n=4;
interp.rz.n=6;
interp.sz.n=6;
interp.qq.n=0;

r1=integrules('quad9',[R23(:) S23(:)*[1 0 0]]);
interp.rr.lam=repmat(struct('val',zeros(9,1),'r',zeros(9,1),'s',zeros(9,1)),1,6);
for jk=1:6
	interp.rr.lam(jk).val=r1.N(jk,:)';
	interp.rr.lam(jk).r=r1.Nr(jk,:)';
	interp.rr.lam(jk).s=r1.Ns(jk,:)';
end
interp.rz.lam=interp.rr.lam;

r1=integrules('quad9',[R32(:) S32(:)*[1 0 0]]);
interp.ss.lam=repmat(struct('val',zeros(9,1),'r',zeros(9,1),'s',zeros(9,1)),1,6);
for jk=1:6
	interp.ss.lam(jk).val=r1.N(jk,:)';
	interp.ss.lam(jk).r=r1.Nr(jk,:)';
	interp.ss.lam(jk).s=r1.Ns(jk,:)';
end
interp.sz.lam=interp.ss.lam;

r1=integrules('quad9',[R22(:) S22(:)*[1 0 0]]);
interp.rs.lam=repmat(struct('val',zeros(9,1),'r',zeros(9,1),'s',zeros(9,1)),1,4);
for jk=1:4
	interp.rs.lam(jk).val=r1.N(jk,:)';
	interp.rs.lam(jk).r=r1.Nr(jk,:)';
	interp.rs.lam(jk).s=r1.Ns(jk,:)';
end


for jk=1:9
	tying=tying_mitc9(r0.w(jk,1),r0.w(jk,2));
	interp.rr.int(jk,:)=tying.rr;
	interp.ss.int(jk,:)=tying.ss;
	interp.rs.int(jk,:)=tying.rs;
	interp.rz.int(jk,:)=tying.rz;
	interp.sz.int(jk,:)=tying.sz;
end


% Integration points definitions
%

integr.n=[9 2]; % Number of in-plane and transverse Gauss points

% In-plane
for jk=1:9
	integr.lam(jk).val=r0.N(jk,:)';
	integr.lam(jk).r=r0.Nr(jk,:)';
	integr.lam(jk).s=r0.Ns(jk,:)';
end
integr.wrs=r0.w(:,4);


% Transverse
integr.z=rga2;
integr.wz=wga2;

end



function tying=tying_mitc4(r,s)
%
% Output: tying point shape functions at (r,s)
%	tying.xx: value of shape functions for strain xx (Ntying_xx,1)

tying=struct('rz',zeros(2,1),'sz',zeros(2,1));

tying.rz(1)=1/2*(1-s);
tying.rz(2)=1/2*(1+s);
%
tying.sz(1)=1/2*(1-r);
tying.sz(2)=1/2*(1+r);

end


function tying=tying_mitc9(r,s)
%
% Output: tying point shape functions at (r,s)
%	tying.xx: value of shape functions for strain xx (Ntying_xx,1)

ga=.774596669241483;
ga2=.5773502691896258;

tying=struct('rr',zeros(6,1),'ss',zeros(6,1),'rs',zeros(4,1),'rz',zeros(6,1),'sz',zeros(6,1));

tying.rr(1)=1/(4*ga2*ga^2)*(ga2-r)*s*(s-ga);
tying.rr(2)=1/(4*ga2*ga^2)*(ga2+r)*s*(s-ga);
tying.rr(3)=1/(2*ga2*ga^2)*(ga2-r)*(ga^2-s^2);
tying.rr(4)=1/(2*ga2*ga^2)*(ga2+r)*(ga^2-s^2);
tying.rr(5)=1/(4*ga2*ga^2)*(ga2-r)*s*(s+ga);
tying.rr(6)=1/(4*ga2*ga^2)*(ga2+r)*s*(s+ga);
%
tying.rs(1)=1/(4*ga2^2)*(ga2-r)*(ga2-s);
tying.rs(2)=1/(4*ga2^2)*(ga2+r)*(ga2-s);
tying.rs(3)=1/(4*ga2^2)*(ga2-r)*(ga2+s);
tying.rs(4)=1/(4*ga2^2)*(ga2+r)*(ga2+s);
%
tying.ss(1)=1/(4*ga2*ga^2)*r*(r-ga)*(ga2-s);
tying.ss(2)=1/(2*ga2*ga^2)*(ga^2-r^2)*(ga2-s);
tying.ss(3)=1/(4*ga2*ga^2)*r*(r+ga)*(ga2-s);
tying.ss(4)=1/(4*ga2*ga^2)*r*(r-ga)*(ga2+s);
tying.ss(5)=1/(2*ga2*ga^2)*(ga^2-r^2)*(ga2+s);
tying.ss(6)=1/(4*ga2*ga^2)*r*(r+ga)*(ga2+s);
%
tying.rz=tying.rr;
tying.sz=tying.ss;

end

