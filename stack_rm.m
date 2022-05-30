function [Up,out1]=stack_rm(Up,st,st1,opt)

% stack handling routines for use by SDT
%
%       Syntax: Up=stack_rm(Up,typ,name);
%               Up=stack_rm(Up,typ);
%               Up=stack_rm(Up,'',name);
%               Up=stack_rm(Up,'',{namelist});
%               [Up,r1]=stack_rm(Up,typ,name,'Get'); %also outputs removed entries
%               [Up,r1]=stack_rm(Up,typ,name,'GetData');

%       Etienne Balmes, Guillaume Vermot des Roches
%       Copyright (c) 2001-2018 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       $Revision: 1.17 $  $Date: 2022/05/24 09:12:11 $

out1=[];
if isfield(Up,'Stack')||isfield(Up,'vfields');
elseif ischar(Up)&&strcmp(Up,'cvs');
    Up='$Revision: 1.17 $  $Date: 2022/05/24 09:12:11 $'; return;
else;Up.Stack={}; 
end
if ~isempty(Up.Stack)

 if ischar(st);[unu,st]=comstr(st,1);end %#ok<ASGLU>
 if isempty(st);i1=[1:size(Up.Stack,1)]; 
 elseif iscell(st)&&nargin==2
    i1=find(ismember(Up.Stack(:,1),st));
 elseif st(1)=='#'
  % stack_rm(Up,'#(?!word)^.*','') to only keep "word"
    i1=find(~cellfun('isempty',regexp(Up.Stack(:,1),st(2:end))));
 elseif nargin>3&&~isempty(strfind(lower(opt),'low'));
  i1=find(strcmp(lower(st),lower(Up.Stack(:,1)))); %#ok<STCI>
 else; i1=find(strcmpi(st,Up.Stack(:,1)));
 end

 if nargin==2; st1=''; 
 elseif iscell(st1) % Multi remove
  if nargin>3&&~isempty(strfind(opt,'get'))
   out1=cell(length(st1),3);
   for j1=1:length(st1);[Up,out1(j1,:)]=stack_rm(Up,st,st1{j1},'get');end
   if size(out1,1)==1&&~isempty(strfind(opt,'data')); out1=out1{1,3}; end
  else
   [u1,i2]=stack_get(Up,st,st1); Up.Stack(i2,:)=[]; % quicker stra
   %for j1=1:length(st1);Up=stack_rm(Up,st,st1{j1});end
  end
  return
 else; st1=comstr(st1,1); 
 end
 if isempty(i1) % no match for type
 elseif ~isempty(st1) % try matching 
  if st1(1)=='#';
   i2=find(~cellfun('isempty',regexp(Up.Stack(i1,2),st1(2:end))));
  else;i2=find(strcmp(st1,Up.Stack(i1,2)));
  end
  if ~isempty(i2); i1=i1(i2); else; i1=[];end
 %else;  i2=[1:length(i1)]; % remove all (implicit in i1)
 end

 if ~isempty(i1) 
  if nargin>3&&~isempty(strfind(lower(opt),'get'))
   out1=Up.Stack(i1,:);
   if ~isempty(strfind(lower(opt),'data')); out1=out1(:,3);
    if size(out1,1)==1; out1=out1{1}; end
   end
  end
  i2=true(size(Up.Stack,1),1); i2(i1)=false; %i2=[1:size(Up.Stack,1)];i2(i1)=0;i2=find(i2);
  Up.Stack=Up.Stack(i2,:);
 end

end 
