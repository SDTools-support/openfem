function [idof,k,m]=mass2(node,varargin)

%MASS2  concentrated mass element with offset
%
%	As all element functions (see ELEM0), MASS2 is called by FE_MK for
%	model assembly, FEPLOT for structural deformation visualization, ...
%
%	In an model description matrix a group of MASS2 elements starts with a
%	header row [Inf  abs('mass2') 0 ...] followed by element property rows
%       ELT following the format
%	    [n1 M I11 I21 I22 I31 I32 I33 EltID CID X1 X2 X3 MatId ProId]
%
%       PL,IL no material or section properties used by this element
%
%   See sdtweb     eltfun
%	See also help  mass1, bar1, beam1, ...

%	Etienne Balmes
%       Copyright (c) 2001-2017 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

%#ok<*NOSEM,*ASGLU>

% standard calls with one input argument
if ischar(node)
 [CAM,Cam]=comstr(node,1);  
 if comstr(Cam,'integinfo');   idof=[];k=[];m=[];
 elseif comstr(node,'cvs');
  idof='$Revision: 1.28 $  $Date: 2020/04/15 18:54:50 $'; return;
 elseif      comstr(node,'call')||comstr(node,'matcall')
   idof=['[i1,k1,m1] = mass2(node(NNode(elt(cEGI(jElt),1)),:),elt(cEGI(jElt),:),pl,il,[opt(1) jGroup jElt]);'];k=0;
 elseif  comstr(Cam,'groupinit');   idof = '';
 elseif  comstr(node,'node');   idof = 1;
 elseif  comstr(node,'prop');   idof = [14 15 9];
 elseif  comstr(node,'dof');    idof=1+(1:6)'/100;
 elseif  comstr(node,'line');   idof  = 1;
 elseif  comstr(node,'patch');  idof = 1;
 elseif comstr(node,'sci_face'); idof = [1 1 1];
 elseif  comstr(node,'parent'); idof = 'mass1';
 elseif  comstr(node,'pu'); 
 %% #pu : property unit type
  idof={'NodeId'    0  'sdtweb(''mass2'')';
  'M'     0  '';
  'I11'   11  '';
  'I21'   11  '';
  'I22'   11  '';
  'I31'   11  '';
  'I32'   11  '';
  'I33'   11  '';
  'EltId' 0  'EltId'
  'CID'    0 'Local coordinate system'
  'X1',4,'offset along x';'X2',4,'offset along y';'X3',4,'offset along z'
  'MatId' 0  'MatId';'ProId' 0  'ProId'
  };
  if comstr(node,'pustruct')&&nargin>1
   elt=feutil('selelt eltname mass2',varargin{1});
   elt(~isfinite(elt(:,1)),:)=[];
   while ~any(elt(:,end));elt(:,end)=[];end
   idof=idof(1:min(size(elt,2),size(idof,1)),1)';
   for j1=1:size(idof,2);idof(2,j1)=num2cell(elt(:,j1));end
   idof=struct(idof{:});
  end

 elseif  comstr(node,'viewbox'); 
 %% #viewbox : generate a view element associated with the 
  mdl=varargin{1};
  elt=feutil('selelt eltname mass2',mdl);elt(1,:)=[];
  [node,bas]=basis('nodebas-force',mdl);
  cGL=full(basis('transte',bas,node,elt(1)+(1:3)'/100));
  I=diag(elt([3 5 8]))*2-elt([3 4 6;4 5 7;6 7 8]); % Generate inertia tensor
  %I=elt([3 4 6;4 5 7;6 7 8]); % Generate inertia tensor
  [dim,I]=eig(I); % Find principal axis
  I=diag(I);
  % Apply inertia scaling + local to global transformation
  % From I, we find a, b, c length of box along X Y and Z (principal
  % directions) with a~sqrt(Y+Z-X), b~sqrt(X+Z-Y) and c~sqrt(X+Y-Z)
  dim=cGL*dim*diag([sqrt(I(2)+I(3)-I(1)) sqrt(I(1)+I(3)-I(2)) sqrt(I(1)+I(2)-I(3))]); 
  dim=dim'; % Store main directions in row
  n1=feutil('getnode',mdl,elt(1));
  n2=node(node(:,1)==n1(1),5:7);
  dim=dim/norm(dim,'inf')*norm(max(node(:,5:7))-min(node(:,5:7)),'inf')/2; % xxx coef is arbitrary
  
  % [Oxyz;0Axyz;OBxyz;OCxyz] (origin + directions as rows)
  mdl=feutil('objecthexa 999 999',mdl,[n2-sum(dim)/2;dim],1,1,1);  
  mdl.Elt(find(~isfinite(mdl.Elt(:,1)),1,'last'),8)=-1;
  elt=setdiff(mdl.Node(:,1),node(:,1));elt=elt(:);elt(:,2)=n1(1);elt(:,3)=123456;
  mdl=feutil('addelt',mdl,'rigid',elt(:,[2 1 3]));
  idof=mdl;
 else; idof=mass1(node,varargin{:}); 
 end
return
end % of standard calls with one input argument
[elt,pe,ie,opt]=deal(varargin{:});
if opt(1)==-2 
  %% Multi mass2 observation
  error('Not yet implemented');
  out=struct('DOF',feutil('getdof',struct('Node',node,'Elt',elt)), ...
      'cta',[],'bas',{{}}); 
  for jElt=find(isfinite(elt(:,1)))'
   r1=cbush(node,elt(jElt,:),varargin{2:end});
   out.cta=[out.cta;fe_c(out.DOF,r1.DOF,r1.cta)];
   out.bas{size(out.cta,1)/6}=bas;
  end
  out.Elt=elt(isfinite(elt(:,1)),:);
  idof=out;return;
end

% element matrix assembly - - - - - - - - - - - - - - - - - - - - - - - - - -

k=zeros(6,6);idof=(1:6)'/100+elt(1,1);
if length(elt)>9&&elt(10)>0;
 bas=evalin('caller','bas');
 bas=bas(bas(:,1)==elt(10),:);
 if bas(3)~=0; 
  error('Bases should be defined in absolute reference, this is a bug');
 end
else;bas=[];
end

if     opt(1)==0; k=zeros(6,6);  m=offset(node,elt,bas);
elseif opt(1)==2;  k=offset(node,elt,bas); m=[];
else; idof=[];m=[];k=zeros(6,6);
end

% ----------------------------------------------------------------
function m = offset(node,elt,bas); 

  r1=elt(1:8); r1([4 6 7])=-r1([4 6 7]);
  m = sparse([1 2 3 4 5 6 4 5 6 4 5 6],[1 2 3 4 4 4 5 5 5 6 6 6], ...
             r1([2 2 2 3 4 6 4 5 7 6 7 8]));
  
  if length(elt)<10; return;elseif length(elt)<13;elt(13)=0;end
  if any(elt(11:13))||elt(10)~=0
      p=[];
      if elt(10)==-1; warning('Need check that node is in global coord');
          x=elt(11:13)-node(node(:,1)==elt(1),5:7);
      elseif isempty(bas); x=elt(11:13);
      else % use basis to define both offset and local orientation
       n1=basis('gnode',bas,[2 bas(1) 0 0 0 0 0;1 bas(1) 0 0 elt(11:13)]);
       x=diff(n1(:,5:7));p=feval(basis('@bas2cGL'),bas,node);
      end
      r = [0 x(3) -x(2) ;-x(3) 0 x(1);x(2) -x(1) 0];
      if isempty(p);r = [eye(3,3) r;zeros(3,3) eye(3,3)];
      else; r = [p r*p;zeros(3,3) p]; end
      m=r'*m*r;
   end
