function [Up,i1]=stack_set(Up,typ,name,val,varargin)

%STACK_SET stack handling routines for use by SDT
%
%       Syntax: Up=stack_set(Up,typ,name,val)
%               [Up,IndPos]=stack_set(Up,SubStack)
%
%       See also : stack_get, stack_rm

%       Etienne Balmes
%       Copyright (c) 2001-2017 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

%#ok<*NOSEM>
 if ~ischar(Up)
  if isa(Up,'sdth');
  elseif iscell(Up)&&(size(Up,2)==3||size(Up,1)==0)
   Up=struct('Stack',{Up});
   if nargin==2; varargin={Up,typ}; 
   elseif nargin==3; varargin={Up,typ,name};
   elseif nargin==4; varargin={Up,typ,name,val};
   elseif nargin==5; varargin={Up,typ,name,val,varargin};
   end
   [Up,i1]=stack_set(varargin{:}); Up=Up.Stack; return;
  elseif ~isfield(Up,'Stack')||~isa(Up.Stack,'cell'); Up.Stack={}; 
  end
 elseif strcmp(Up,'cvs')
  Up='$Revision: 1.16 $  $Date: 2021/09/14 12:56:21 $';return;
 else; error('You must provide a structure for stack_set'); 
 end
if nargin==2 % set substack
 % two step procedure to bypass pointer callback problems
 r1=struct('Stack',[]);try; r1.Stack=Up.Stack;catch;r1.Stack={};end
 i1=zeros(size(typ,1),1);
 if isempty(r1.Stack)
   r1.Stack=typ; i1=(1:size(typ,1))';  
 else % Overwrite stack
  for j1=1:size(typ,1)
   [r1,i1(j1)]=stack_set(r1,typ{j1,1:3});
  end
 end
 Up.Stack=r1.Stack;
elseif nargin>4
  r1=struct('Stack',[]);try; r1.Stack=Up.Stack;catch;r1.Stack={};end
  r1=stack_set(r1,typ,name,val);i1=zeros(length(varargin)/3,1);
  for j1=1:3:length(varargin);
   [r1,i1((j1+2)/3)]=stack_set(r1,varargin{j1+(0:2)});
  end
  Up.Stack=r1.Stack;

elseif isempty(Up.Stack) % set a new value
 Up.Stack(1,1:3)={typ,name,val};i1=1;
else % search if value exist
 r2=Up.Stack;
 if ~isempty(typ) 
  i1=find(strcmpi(typ,r2(:,1))&strcmpi(name,r2(:,2)));
 elseif ~isempty(name)
  i1=find(strcmp(name,r2(:,2)));
 else; error('You must provide at least a non empty typ or name');
 end
 if isempty(i1); i1=size(r2,1)+1;end
 if length(i1)==1
   if isempty(typ); Up.Stack(i1,2:3)={name,val};
   elseif isempty(name)
       if size(Up.Stack,1)>=i1&&ischar(Up.Stack{i1,2})
        Up.Stack(i1,[1 3])={typ,val};
       else; Up.Stack(i1,1:3)={typ,'',val};
       end
   else; Up.Stack(i1,1:3)={typ,name,val}; 
   end
 else
   sdtw('_nb','Multiple match (%s,%s) in stack_set, probably an error',typ,name);
   Up.Stack(i1(1),1:3)={typ,name,val};
 end

end
