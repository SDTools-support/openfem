function [idof,k,m]=celas(node,varargin)

%CELAS	element function for scalar springs and penalized rigid links
%
%	As all element functions (see ELEM0), CELAS is called by FE_MK for
%	model assembly, FEPLOT for structural deformation visualization, ...
%
%	In an model description matrix a group of CELAS elements starts with a
%	header row [Inf  abs('celas') 0 ...] followed by element property rows
%       ELT following the format
%	    [n1 n2 DofId1 DofId2 ProId EltId  Kv Mv Cv Bv]
%         with
%	   n1,n2  node numbers of the connected nodes. Grounded springs are
%		  obtained by setting n1 or n2 to 0. 
%	   DofId  Identification of selected DOFs.
%		For RIGID LINKS, the first node defines the rigid body motion.
%		  DofId (>0) defines which DOFs of the slave are connected
%		  by the constraint. Thus [1 2 123 0 0 0 1e14] will only impose
%		  translations of node 2 are imposed by motion of node 1, while
%		  [1 2 123456 0 0 0 1e14] will also penalize the difference in
%		  rotations.
%		For SCALAR SPRINGS, DofId1 (negative) defines which DOFs of
%		  node 1 are connected to which of node 2. Use of negative DofId1 will
%    also only activate specified DOF at the node beyond the initial ones.
%    DofId2 can be used to specify different DOFs on the 2 nodes. For example:
%		  [1 2 -123 231 0 0 1e14] connects DOFs 1.01 to 2.02, etc.
%	   ProId  Property identification number (see format below)
%	   Kv     Stiffness value used as a weighting  associated with the 
%		  constraint. If Kv is zero (or not given), the default value
%		  in the element property declaration is used. If this
%		  is still 0, Kv is set to sdtdef('kcelas').
%
%	Element property rows for CELAS elements take the form
%	   [ProId fe_mat('p_spring','SI',1)  KvDefault]
%
%	Warning : when seeking rigid rotations, be sure that the master node
%	  n1 has rotational DOFs. 
%
%	See sdtweb     rigid, bar1, eltfun, elem0
%	See also help  rigid

%	Etienne Balmes, Guillaume Vermot des Roches
%       Copyright (c) 2001-2021 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license


% standard calls with one input argument
if ischar(node)
 [CAM,Cam]=comstr(node,1);
 if comstr(Cam,'integinfo'); idof=[];k=[];m=[];
 elseif comstr(node,'cvs')
  idof='$Revision: 1.42 $  $Date: 2021/12/10 10:11:02 $'; return;
 elseif comstr(node,'rhscall'); idof='';
 elseif comstr(node,'call')
   idof = ['[i1,k1,m1] = celas(nodeE, elt(cEGI(jElt),:),pl,' ...
       'il,[opt(1) jGroup jElt],Case);'];
 elseif comstr(node,'matcall')
   idof=celas('call');k=0;
 elseif  comstr(Cam,'groupinit');   idof = '';
 elseif  comstr(node,'node');    idof = [1 2];
 elseif  comstr(node,'prop');    idof = [0 5 6];
 elseif  comstr(node,'dof');    idof=[1+(1:6)/100 2+(1:6)/100];
 elseif  comstr(node,'line');    idof =[1 2];
 elseif  comstr(node,'patch');   idof= [1 2];
 elseif  comstr(node,'edge');   idof =[1 2];
 elseif  comstr(node,'face');    idof= [];
 elseif  comstr(node,'state');    idof= [];
 elseif  comstr(node,'sci_face'); idof = [1 2];
 elseif  comstr(node,'parent');  idof = 'beam1';
 elseif  comstr(node,'test')
   mdl.Node=[1  0 0 0  0 0 0; 2  0 0 0  0 0 1];
   mdl.Elt=[Inf abs('celas') 0; 1 2 -3 3 100 1 1e3]; %[n1 n2 DofId1 DofId2 ProId EltId  Kv Mv Cv Bv]
   mdl=fe_mknl(mdl);
   if mdl.K{2}(1,1)~=1e3; sdtw('_err','bad stiffness'); end
   mdl.Elt(2,7)=0; mdl=fe_mknl(mdl);
   if mdl.K{2}(1,1)~=sdtdef('kcelas'); sdtw('_err','bad stiffness'); end
   mdl.Elt(2,8)=100; mdl=fe_mknl(mdl);
   if mdl.K{2}(1,1)~=0; sdtw('_err','bad stiffness'); end
   mdl.Elt(2,7:9)=[10 10 30];mdl=fe_mknl(mdl);K=mdl.K;K{3}=fe_mknl('assemble',mdl,3);
   if ~isequal(mdl.K{1},mdl.K{2})||K{3}(1,1)~=30;sdtw('_err','celas matrix'); end
   mdl.il=[100 fe_mat('p_spring','SI',1) 10 10 30];
   mdl=fe_mknl(mdl);mdl.K{3}=fe_mknl('assemble',mdl,3);
   if ~isequal(mdl.K,K);sdtw('_err','celas matrix with pl'); end
   mdl.Elt(2,7:10)=[0 0 0 0]; mdl=fe_mknl(mdl);mdl.K{3}=fe_mknl('assemble',mdl,3);
   if ~isequal(mdl.K,K);sdtw('_err','celas matrix with pl'); end
 end
 return
elseif nargin==0; idof=[];return;
end %of standard calls with 1 input argument

[elt,pl,il,opt,Case]=deal(varargin{:});  %#ok<*ASGLU>

if opt(1)==-2&&size(elt,1)>1
 %% Multi celas observation
 out=struct('DOF',feutil('getdof',struct('Node',node,'Elt',elt)), ...
  'cta',[]); %,'bas',{{}});
 cEGI=find(isfinite(elt(:,1)));
 for jElt=1:length(cEGI)
  r1=celas(node,elt(cEGI(jElt),:),varargin{2:end});
  out.cta=[out.cta;fe_c(out.DOF,r1.DOF,r1.cta)];
  %out.bas{size(out.cta,1)/6}=bas;
 end
 out.Elt=elt(isfinite(elt(:,1)),:);
 idof=out;return;
end


i2=abs(sprintf('%i',abs(elt(1,3))))-48; % DOF selected for node 1
i3=abs(sprintf('%i',abs(elt(1,4))))-48; % DOF selected for node 2

if ~isempty(il); i1 = find(il(:,1)==elt(1,5));
 if ~isempty(i1)
      i1=il(i1,3:end);if length(i1)<4;i1(1,4)=0;end
      if i1(4);i1(4)=i1(1)*i1(4);end % In p_spring loss rather than Kv
      elt(1,7:10)=i1(1:4); 
 end
end
if size(elt,2)<10; elt(1,10)=0;end
if elt(1,7)==0&&all(elt(8:10)==0)&&opt(1)~=-2 % allow pure damper
   elt(1,7)=sdtdef('kcelas'); 
   sdtw('_nb','use default kcelas %g (%i-%i)',elt(1,[7 1 2]))
elseif elt(1,7)==-1;elt(1,7)=0;
end
elt(1,6+find(elt(1,7:9)<0))=0; % ignore negative values of k m c
if isempty(i2); warning('OpenFEM:NoDOF','No DOF selected'); end

idof=[];
if opt(1)==3; elt(1,7)=elt(1,9);     % viscous
elseif opt(1)==4;elt(1,7)=elt(1,10); % hysteretic
end

if elt(1,3)>0&&size(node,1)>1&&elt(1,2)  % penalized linearized rigid link
  
       x = node(:,1:3);x=x(2,:)-x(1,:);
       x = [0 x(3) -x(2) ;-x(3) 0 x(1);x(2) -x(1) 0];
       r = [eye(3,3) x;zeros(3,3) eye(3,3)];
       r = sparse([r(i2,:) -eye(length(i2),length(i2))]);
       idof =[elt(1,1)+(1:6)'/100; elt(1,2)+i2(:)/100];
       i3=[1:6 i2+6]; 

else % standard spring connection

       if isequal(i3,0); i3=i2; elseif isequal(i2,0);i2=i3;end
       if length(i3)~=length(i2)
         error('non-matching numbers of DOFs for standard spring connection');
       end
       r = sparse([eye(length(i2),length(i2)) -eye(length(i2),length(i2))]);
       if any(elt(1,1:2)==0)||opt(1)==-2
        idof =[elt(1,1)+i2(:)/100; elt(1,2)+i3(:)/100];
        r(:,idof<1)=0;%k(:,i1)=0;k(i1,:)=0;
       end
       i3=[i2 i3+6]; 

end

if opt(1)==-2
 %% #Cta build return observability of celas
 i1=any(fix(idof),2);
 idof=struct('DOF',idof(i1),'cta',r(:,i1));
 if ~isempty(Case.cGL) % in this case need to expand to all active DOF
  jElt=evalin('caller','jElt');
  DofPos=Case.GroupInfo{Case.jGroup,1}(:,jElt)+1;
  i1=DofPos(DofPos>0); dof=Case.mDOF(i1);
  idof.cta=fe_c(dof,idof.DOF,idof.cta); idof.DOF=dof;
  idof.cta=idof.cta*Case.cGL(i1,i1)'; 
  
 end
 return;
 
else
 %% Assemble
 k=zeros(12,12);m=k;
 k(i3,i3)=elt(1,7)*(r'*r);
 m(i3,i3)=elt(1,8)*(r'*r);
 
 if ~isempty(Case.cGL)
  in1=evalin('caller','DofPos(:,jElt)');in1=double(in1)+1;
  in2=in1~=0; T=zeros(12);in1=in1(in2);T(in2,in2)=Case.cGL(in1,in1)';Tt=T';
  %if ~all(in2); k=k(in2,in2); m=m(in2,in2);in1=in1(in2);end
  %  T=Case.cGL(in1,in1); Tt=T';
  k=T*k*Tt;m=T*m*Tt;
 end
 
 if opt(1)==2; k=m;
 elseif opt(1)==5 % NL tangent stiffness
 elseif opt(1)>4; idof=[]; k=[];
 end
end

