function [out,out1,out2]=mck3(CAM,varargin);

%MCK3 element function for the 3-node 3-DOF axial bar element
%
%	In an model description matrix a group of MCK3 elements starts with a
%	header row [Inf  abs('mck3') 0 ...] followed by element property rows
%       ELT following the format
%	    [n1 n2 n3 MatId ProId EltId M1 M2 K1 K2 C1]
%         with
%	   n1,n2,n3  node numbers of the tips of the element
%	   MatId     (optional) material property identification number 
%	   ProId     (optional)  element property identification number
%	   EltId     element identifier (optional, default 0)
%      M1,M2     mass of nodes n1 and n2
%      K1,K2     stiffness
%      C1        viscous damping coeff
%          
%       12-12-01
%
%       Copyright (c) 1990-2009 by SDTools
%       All Rights Reserved.

if comstr(CAM,'cvs')
 out='$Revision: 1.8 $  $Date: 2009/05/28 16:42:00 $'; return;
end

% standard calls with one input argument
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1);
 % Build Constit, Integ, and Elmap for later integration
 if comstr(Cam,'integinfo')

  out=[]; out1=[]; out2=[]; 

 elseif comstr(Cam,'matcall'); out=mck3('call');  out1=0; % Call, SymFlag
 elseif comstr(Cam,'call')  % call for matrix assembly

  out='[k1,m1]=mck3(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,elmap);';

 elseif  comstr(Cam,'groupinit');out = '';
 elseif  comstr(Cam,'node');    out = [1:3];
 elseif  comstr(Cam,'prop');    out = [4 5 6];
 elseif  comstr(Cam,'dof');     out = [1.03 2.03 3.03]';
 elseif  comstr(Cam,'edge');    out = [1 2]; 
 elseif  comstr(Cam,'face');    out = [1 2];
 elseif  comstr(Cam,'sci_face');out = [1 2 3];
 elseif  comstr(Cam,'line');    out = [1 2]; 
 elseif  comstr(Cam,'patch');   out = [1 2];
 elseif  comstr(Cam,'parent');  out = 'beam1';
 % Basic matrix test of the element - - - - - - - - - - - - - - - - - -
 elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,5);
   model.Node=[1 0 0 0 0 0 0;2 0 0 0 0 0 1;3 0 0 0 0 0 .5];
   model.Elt=[Inf abs('mck3') 0 0 0 0 0 0 ;1 2 3 0 0 0 100 200 1000 2000 2];
   model.pl=[]; model.il=[];
   [constit,integ,elmap]=mck3('integinfo',[0;0],model.pl,model.il);
   [m,k]=mck3(model.Node,model.Elt(2,:),[0 0 0 0 0 0 0 0 0],int32(integ),constit,elmap);
 else sdtw('''%s'' unknown',CAM); 
 end
return
end % of standard calls with one input argument

% element matrix assembly - - - - - - - - - - - - - - - - - - - - - - - - -
node=CAM; 
elt=varargin{1}; 
point=varargin{2}; typ=point(5);
integ=varargin{3};
constit=varargin{4}; 
if nargin<6; elmap=[];
else elmap=varargin{5};
end

if any([0 1]==typ)	% - - - - - - - - stiffness
  k = zeros(3,3);
  k1= elt(9);
  k2= elt(10);
  k([1 2 4 5 6 8 9])=[k1 -k1 -k1 k1+k2 -k2 -k2 k2];
else   k = []; end

if any([0 2]==typ)  % - - - - - - - - mass
  m = zeros(3,3);
  m3 = min(elt(10)/(2*pi*1000)^2, .01*(elt(7)+elt(8))/2); % param 
  m([1 5 9])=[elt([7 8]) m3];
  if typ==2; k=m;m=[]; end
else   m = [];end

out=k; out1=m;

if any(3==typ)  % - - - - - - - - damping
  out = zeros(3,3);
  c = elt(11);
  out([1 3 7 9])=[c -c -c c];   
end

if ~any([0 1 2 3]==typ)
 error(sprintf('Matrix type %i not supported by mck3',typ))
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


