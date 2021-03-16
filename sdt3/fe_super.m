function [idof,o1,o2,o3,o4,o5,o6]=fe_super(node,elt,pl,il,opt,ElemF)

%FE_SUPER generic element function for super-element support
%
%	Syntax: out=fe_super('Command',ElemF)
%		[idof,k,m,c]=fe_super(node,elt,pl,il,opt,ElemF)
%
%	Usual elements are supported by element functions (beam1, bar1, ...)
%	FE_SUPER provides the alternative to store all the element properties
%	as global variables of the form SEName with fields
%
%	Opt, Node, DOF, Line, Patch, K, TR, Elt, Ref
%
%	For details on superelement fields look at  doc('fe_super')
%	For an example look at the D_CMS2 demonstration
%
%	See also help fesuper, fe_mk, d_cms2 (demo)
%                doc  fem, trd

%	E. Balmes
%       Copyright (c) 2001-2003 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       $Revision: 1.5 $  $Date: 2006/05/29 09:02:06 $


% calls with one input argument
if nargin==1 
  error('Not a valid call to fe_super');
% standard element function calls with one input argument
elseif  any(nargin==[2 3 4])

ElemF=elt; s=who('global',['SE' ElemF]);
if ~isempty(s) eval(['global SE' ElemF ';o1=isstruct(SE' ElemF ');']); 
else o1=0; end

if ~o1 % there is an element function associated - - - - - 

 o1=2; % all element functions are generic
 if node==1 | comstr(node,'call')
   eval(['idof=' ElemF '(''call'',0);'],'idof=[];');
 elseif node==1 | comstr(node,'rhscall')
   eval(['idof=' ElemF '(''rhscall'');'],'idof=[];');
 elseif comstr(node,'node')
   eval(['idof=' ElemF '(''node'');'],'idof=[];');
 elseif comstr(node,'prop')
   ans=[];   eval([ElemF '(''prop'');'],''); idof=ans;
   if isempty(idof) idof=length(fe_super('node',ElemF))+[1 2 3];
     if idof(1)==1 idof=[0 0 0];end
   end
 elseif  comstr(node,'dof')
   eval(['idof=' ElemF '(''dof'');'],'idof=[];');
 elseif  comstr(node,'face')
   eval(['idof=' ElemF '(''face'');'],'idof=[];');
 elseif  comstr(node,'sci_face')
   eval(['idof=' ElemF '(''sci_face'');'],'idof=[];');
   if isempty(idof)
     eval(['idof=' ElemF '(''face'');'],'idof=[];');
   end
 elseif comstr(node,'edge')
   eval(['idof=' ElemF '(''edge'');'],'idof=[];');
 elseif node==500 | comstr(node,'line')
   eval(['idof=' ElemF '(''line'');'],'idof=[];');
 elseif  node==501 | comstr(node,'patch')
   eval(['idof=' ElemF '(''patch'');'],'idof=[];');
 elseif  comstr(node,'parent')
   eval(['idof=' ElemF '(''parent'');'],'idof=ElemF;');
   if isempty(idof) idof=ElemF; end
 end

else % this is a true super-element - - - - - - - - - -

 ElemF=['SE' ElemF];
 eval(['global ' ElemF ';SEopt=' ElemF '.Opt;'],'');
 if node==1 | comstr(node,'call')
   if isempty(SEopt) idof=[ElemF ' no options defined'];
   else
    opt=pl; if opt(1)==0 idof='[i1,k1,m1]='; else idof='[i1,k1]=';end
    idof=[idof 'fe_super(node,elt(cEGI(jElt),:),pl,il,' ...
          '[opt(1) EGID jElt],ElemF);'];
   end

 elseif comstr(node,'node')    eval(['idof = ' ElemF '.Node;']);
   if SEopt(1,1)==1 idof=idof(:,1); else idof=1:size(idof,1);end
 elseif  comstr(node,'dof')    eval(['idof = ' ElemF '.DOF;']);
 elseif  comstr(node,'prop')
    if SEopt(1)==2 eval(['idof=' ElemF '.Node;']); idof=[0 size(idof,1)+[2 3]];
    else idof=[0 0 0]; end
 elseif  comstr(node,'line')   eval(['idof = ' ElemF '.Line;']);
 elseif  comstr(node,'patch')  eval(['idof = ' ElemF '.Patch;']);
 elseif  comstr(node,'parent') 
   idof=[];eval(sprintf('if isfield(%s,''Parent'') idof = %s.Parent;end',ElemF,ElemF));
   if isempty(idof) idof=ElemF; end
 end
 o1=SEopt;
end % of is a true superelement
return
end % of standard calls with one input argument

% element matrix assembly - - - - - - - - - - - - - - - - - - - - - - - - - -

error('Superlements are not supported by OpenFEM');

