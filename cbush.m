function [idof,k,m]=cbush(node,varargin);

%CBUSH element function for the 2-node spring/damper element
%
%       The element property row of a cbush element takes the form
%
%       n1 n2 MatId ProId EltId x1 x2 x3 EDID S OCID S1 S2 S3
%
%           EDID : coordinate system for element orientation
%           OCID : coordinate system for offset definition
%           Si   : components of offset in OCID
%
%	See sdtweb     p_spring
%	See also help  beam1, elem0

%	Etienne Balmes, Guillaume Vermot des Roches
%       Copyright (c) 2001-2021 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

%#ok<*ASGLU,*NOSEM>
%% #ElemF standard calls with one input argument
if ischar(node)
 [CAM,Cam]=comstr(node,1); 
 if comstr(Cam,'integinfo');%constit integ,elmap 
     idof=[varargin{1}];k=[];m=[];
 elseif comstr(node,'cvs');
  idof='$Revision: 1.54 $  $Date: 2022/02/09 15:19:15 $'; return;
 elseif      comstr(Cam,'call')
   idof = ['[i1,k1,m1] = cbush(nodeE,elt(cEGI(jElt),:),pl,il,' ...
       '[opt(1) jGroup jElt],Case);'];
 elseif comstr(node,'matcall');idof=cbush('call');k=0;
 elseif comstr(node,'rhscall'); idof='';
 elseif  comstr(Cam,'groupinit');  idof = '';
 elseif  comstr(Cam,'node');  idof = [1 2];
 elseif  comstr(Cam,'prop');  idof = [3 4 5];
 elseif  comstr(Cam,'dof');   idof=[1+(1:6)/100 2+(1:6)/100];
 elseif  comstr(Cam,'line');  idof = [1 2];
 elseif  comstr(Cam,'face');   idof =[];
 elseif  comstr(Cam,'edge');   idof =[1 2];
 elseif  comstr(Cam,'patch');  idof = [1 2];
 elseif  comstr(Cam,'parent');   idof = 'beam1';
 elseif  comstr(Cam,'sci_face'); idof = [1 2 2];
 elseif  comstr(Cam,'test')
  warning('OpenFEM:TEST','cbush : tests not implemented ');return
 end
return
end % of standard calls with one input argument

% element matrix assembly - - - - - - - - - - - - - - - - - - - - - - - - -
[elt,pl,il,opt,Case]=deal(varargin{:});bas=Case.bas;

if size(node,1)>2||size(node,2)==7
 %% #Obs
 if size(elt,1)==1&&isfinite(elt(1))
  i1=node(:,1);node=node([find(i1==elt(1));find(i1==elt(2))],:);
  nodeE=node(:,5:7);
 elseif opt(1)==-2 
  %% Multi cbush observation
  out=struct('DOF',feutil('getdof',struct('Node',node,'Elt',elt)), ...
      'cta',[],'bas',{{}}); 
  for jElt=find(isfinite(elt(:,1)))'
   r1=cbush(node,elt(jElt,:),varargin{2:end});
   out.cta=[out.cta;fe_c(out.DOF,r1.DOF,r1.cta)];
   out.bas{size(out.cta,1)/6}=bas;
  end
  out.Elt=elt(isfinite(elt(:,1)),:);
  idof=out;return;
 else;error('Not an expected case');
 end
elseif size(node,1)==1||~any(node(2,:)); nodeE=node([1;1],1:3);
else; nodeE=node(:,1:3);
end
%% #Ass
il=il(il(:,1)==elt(4),:); k=[]; % (k holder for type 4)
if isempty(il);error('Property %i not found',elt(4));
elseif length(il)<15; il(15)=0; 
end
stp=fe_mat('typepstring',il(2));
switch stp
 case 'p_spring.2' % usual
 case 'p_spring.1' % allow homogenous
  sdtw('_nb',['ProId %i is of type %s but defines a cbush.\n'...
   'k,m,c Eta properties are resp. duplicated in k1:6, m, c1:6, Eta'],il(1),stp)
  il0=il;   il(3:8)=il0(3); il(9:14)=il0(5); il(15)=il0(6); il(20)=il0(4);
 case 'p_spring.4' % full matrices
  il0=il; if length(il0)<76; il0(1,76)=0; end
  if opt(1)==3; k=reshape(il0(39:74),6,6)';
  else; k=reshape(il0(3:38),6,6)'; if opt(1)==4; k=k*il0(75); end % fill stiffness
  end
 otherwise; error('ProId %i type %s cannot define a cbush',il(1),stp);
end
if isempty(k); il(2+find(il(3:14)<0))=0; end % ignore negative terms

idof=[elt(1)+(1:6)/100 elt(2)+(1:6)/100]';
if  size(elt,2)<10; elt(10:14)=[.5 0 0 0 0]; end
% element coordinate system and S used
if size(elt,2)<11 || (size(elt,2)>10&&elt(11)<1)
   if elt(10)==0&&(size(elt,2)<11||elt(11)==0); elt(10)=.5; end
   if size(nodeE,1)==1; error('Element coordinates not defined with a single node');end
   S=[elt(10) 1-elt(10)]*nodeE;
elseif elt(11)>0 % offset defined (OCID in elt(11))
 % (neg values used to fix S=0 in non default mode)
   i2=find(bas(:,1)==elt(11));
   r1=basis([1 bas(i2,1) 0 0 elt(12:14)],bas(i2,:));
   S=r1(5:7);
end
if elt(9)==0;cLG=eye(3);
elseif elt(9)<1  % EDID Coordinate system defined by nodes 
  b=diff(nodeE);
  if any(elt(7:8))||rem(elt(6),1)~=0 % (x2,x3 ~=0 || x1 non integer)
    if nnz(b)==0;% Allow definition of local x_e for zero length
     cLG=sp_util('basis',elt(6:8),b)';  %b=b(:,[3 1 2])';
    else; % Non zero length hence x_e given by b, y defined in elt(6:8)
     cLG=sp_util('basis',b,elt(6:8))';
    end
  elseif nnz(b)==0 % Use ref to define x
    error('External NodeId reference not implemented')
  else % Use ref to define ye
   if elt(6)
    n0=evalin('caller','model.Node');n0=n0(n0(:,1)==elt(6),:);
    if any(n0(2:3));error('Global node expected');end
   else; n0=[0 0 0 0 nodeE(1,1:3)+1];
   end
   cLG=sp_util('basis',b,n0(5:7)-nodeE(1,1:3))';
  end
else
  % EDID(elt9) Externally defined coordinate system (pre-emptive)
  % xe is first column hence ' needed to observe local x component
  cLG=reshape(bas(bas(:,1)==elt(9),7:15),3,3)';
  if any(elt(6:8));error('Definition of DID is incompatible with x1,x2,x3');end
end

x=S-nodeE(1,:);
tr1=[1 0 0   0 x(3) -x(2);0 1 0 -x(3) 0 x(1);0 0 1 x(2) -x(1) 0;
     0 0 0  1 0 0;0 0 0  0 1 0;0 0 0  0 0 1];
x=S-nodeE(2,:);
tr2=[1 0 0   0 x(3) -x(2);0 1 0 -x(3) 0 x(1);0 0 1 x(2) -x(1) 0;
     0 0 0  1 0 0;0 0 0  0 1 0;0 0 0  0 0 1];

tr=[tr1 -tr2];
% standard cbush element
if any([0 1 4 5]==opt(1))	% stiffness or hysteretic damping
 if isempty(k);  k = diag(il(3:8)); if opt(1)==4; k=k*il(15);end; end
 k=of_mk('xkx_trans',cLG,k);%x'*k*x
 k=tr'*k*tr;
elseif opt(1)==3  % viscous damping
 if isempty(k); k = diag(il(9:14)); end
 k=of_mk('xkx_trans',cLG,k); k=tr'*k*tr;
else; k = zeros(12,12);
end

m=zeros(12,12);

if any([0 1 2 3 4 5 7 9 100]==opt(1)) % xxx NL tgt may require nodeE update 
elseif opt(1)==-1 % return CBUSH without any local basis reference

  be=cLG(:,3);[r1,i1]=max(abs(be));if be(i1)<0; be=-be;end
  be=basis(diff(nodeE),be);

 % n1 n2  MatId ProId EltId x1 x2 x3   EDID S   OCID   S1 S2 S3
 idof=[elt(1:5)              be(:,2)'  -2 0.5  -1  (S-nodeE(1,:))*be'];
 return;
elseif opt(1)==-2; % Return observability of cbush 
  idof=struct('DOF',idof,'cta',sparse([cLG zeros(3);zeros(3) cLG]*-tr), ...
      'bas',cLG); 
  return;
else; idof=[]; k=[];
 %error('Matrix type %i not supported by cbush',opt(1))
end
% 
if ~isempty(Case.cGL)
  in1=evalin('caller','DofPos(:,jElt)');in1=double(in1)+1;
  in2=in1~=0; T=zeros(12);in1=in1(in2);T(in2,in2)=Case.cGL(in1,in1)';%Tt=T';
  if ~all(diag(T)==1) % keep a valid test for grounded cbush
   if (~any(node(2,:))&&~isequal(diag(T),[ones(6,1);zeros(6,1)]))
    error('No handling of CBUSH with nodal DID, use CBUSH CID');
   end
  end
  %if ~all(in2); k=k(in2,in2); m=m(in2,in2);in1=in1(in2);end
  %T=Case.cGL(in1,in1); Tt=T';
  %k=T*k*Tt;m=T*m*Tt; % k_gg = cGL * k_L * cLG
end

