function [out,out1,out2]=mass1(node,varargin)

%MASS1  concentrated diagonal mass element
%
%	As all element functions (see ELEM0), MASS1 is called by FE_MK for
%	model assembly, FEPLOT for structural deformation visualization, ...
%
%	In an model description matrix a group of MASS1 elements starts with a
%	header row [Inf  abs('mass1') 0 ...] followed by element property rows
%       ELT following the format
%	    [n mxx myy mzz ixx iyy izz EltId]
%         with
%	   n      node number for the concentrated mass
%	   mxx myy mzz ixx iyy izz, the concentrated mass is associated to the
%	          DOFs .01 to .06 (xyz translations/rotations) of the conside-
%		  node and given by diag([mxx myy mzz ixx iyy izz])
%
%       PL,IL no material or section properties used by this element
%
%   See sdtweb     eltfun
%	See also help  bar1, beam1, ...

%	Etienne Balmes
%       Copyright (c) 2001-2020 by INRIA and SDTools
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       All Rights Reserved.

% standard calls with one input argument
if nargin>1; elt=varargin{1};end %#ok<*NASGU,*ASGLU>
if ischar(node)
 [CAM,Cam]=comstr(node,1); 
 if comstr(Cam,'integinfo')
   out=[];out1=[];out2=[];
 elseif comstr(node,'cvs')
  out='$Revision: 1.29 $  $Date: 2020/04/15 18:54:50 $'; return;
 elseif comstr(node,'matcall'); out=mass1('call');out1=0;
 elseif comstr(node,'call')
   out = ['if opt(1)==0; [k1,m1] =mass1(nodeE,elt(cEGI(jElt),:),opt); ' ...
          'else;k1=mass1(nodeE,elt(cEGI(jElt),:),opt); end'];
 elseif  comstr(Cam,'groupinit');  out = '';
 elseif  comstr(node,'node');   out = 1;
 elseif  comstr(node,'prop');   out = [0 9 8];
 elseif  comstr(node,'dof');    out=1+(1:6)'/100;
 elseif  comstr(node,'line');   out  = 1;
 elseif  comstr(node,'patch');  out = 1;
 elseif  comstr(node,'edge');   out  = [];
 elseif  comstr(node,'face');   out = [];
 elseif comstr(node,'sci_face'); out = [1 1 1];
 elseif  comstr(node,'parent'); out = 'mass1';
 elseif comstr(node,'state')||comstr(node,'rhscall');  out='';% call for state update
 end
return
end % of standard calls with one input argument

% element matrix assembly - - - - - - - - - - - - - - - - - - - - - - - - - -

%k=zeros(6,6);idof=[1:6]'/100+elt(1,1);
elt=varargin{1};opt=varargin{2}; if size(elt,2)<7;elt(1,7)=0;end
if isempty(opt)&&nargin==6;opt=varargin{4};end% nlspring
if     opt(1)==0; out=zeros(6,6);out1 = diag(elt(2:7));
elseif opt(1)==2; out = diag(elt(2:7));
elseif any(opt(1)==[1 3 5 9]); out=zeros(6,6);
elseif opt(1)==7 % gyro coupling, local rotating frame
 if any(elt(2:4)-mean(elt(2:4))) % rot inertia<>0 or not equal translation mass
  out=zeros(6);
  sdtw('_nb','ignoring mass1 for type 7 gyroscopic assembly');
 else
  if any(elt(5:7)) 
   sdtw('_nb','ignoring mass1 inertia for type 7 gyroscopic assembly');
  end
  model=evalin('caller','model'); [r1,r3]=feval('fe_cyclic','omega',model); %#ok<FVAL>
  out=2*elt(2)*r3;out(6,6)=0; 
 end
elseif opt(1)==8 % centrifugal softening, local rotating frame
 if any(elt(5:7))||any(elt(2:4)-mean(elt(2:4))) % rot inertia<>0 or not equal translation mass
  out=zeros(6);
  sdtw('_nb','ignoring mass1 for type 8 centrifugal softening assembly');
 else
  model=evalin('caller','model'); [r1,r3]=feval('fe_cyclic','omega',model); %#ok<FVAL>
  out=elt(2)*r3^2;out(6,6)=0; 
 end
elseif opt(1)==70 % global fixed frame
   if abs(elt(5)/elt(6)-1)<1e-5 % z rotation
       out=zeros(6);out([23 28])=+elt(7)*[-1 1];%out(4,5)=1;out(5,4)=1;find(out)
   elseif abs(elt(6)/elt(7)-1)<1e-5 % x rotation
       out=zeros(6);out([30 35])=+elt(5)*[-1 1];%out(5,6)=1;out(6,5)=1;find(out)
   else % y rotation
       out=zeros(6);out([34 24])=+elt(6)*[-1 1];%out(4,6)=1;out(6,4)=1;find(out)
   end
elseif opt(1)==-2 % Position based NL using elements
  elt=feutil('addelt','mass1',elt);
  out=struct('DOF',feutil('getdof',struct('Node',node,'Elt',elt)), ...
      'cta',[]); 
  out.cta=speye(length(out.DOF));

else;error('not a valid type of call')
end

