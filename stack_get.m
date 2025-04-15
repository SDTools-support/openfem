function [val,i1]=stack_get(Up,st,st1,GiveData)

%STACK_GET stack handling routines for use by SDT
%
%       Syntax: [StackRows,index]=stack_get(Up,typ);
%               [StackRows,index]=stack_get(Up,typ,name);
%               [data,index]     =stack_get(Up,typ,name,GiveData);
%        with GiveData='multi' substack return is allowed

%       Etienne Balmes, Guillaume Vermot des Roches
%       Copyright (c) 2001-2025 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license     

%#ok<*NOSEM>

name='';i1=[]; %#ok<NASGU>

if isa(Up,'v_handle');
 Up=Up.GetData;
 %h=Up.GetVHandleHandle; %if ~isempty(h)&&~iscell(h);Up=get(h,'userdata');end
 if ~isfield(Up,'Stack');Up.Stack={};end
elseif isfield(Up,'Stack')||isa(Up,'sdth');
elseif ischar(Up)&&strcmp(Up,'cvs');
    val='$Revision: 1.22 $  $Date: 2025/04/14 12:57:13 $'; return;
elseif iscell(Up)&&size(Up,2)==3
 Up=struct('Stack',{Up});
 if nargin==2; varargin={Up,st}; 
 elseif nargin==3; varargin={Up,st,st1};
 elseif nargin==4; varargin={Up,st,st1,GiveData};
 end
 [val,i1]=stack_get(varargin{:}); return;
elseif isstruct(Up); Up.Stack=cell(0,3);
else; val=cell(0,3); i1=[]; if nargin==4; val=[]; end
 if ~isempty(Up); sdtw('_nb','Up is not in a valid format'); end
 return;
end
r1=Up.Stack;
if isempty(r1); i1=[];
elseif nargin==1||isempty(st); i1=(1:size(r1,1))';
elseif st(1)=='#'
    i1=find(~cellfun('isempty',regexp(r1(:,1),st(2:end))));
else;i1=find(strcmpi(st,r1(:,1)));
end

if isempty(i1); val=cell(0,3);
elseif nargin>2
  if isempty(st1);i2=1:length(i1);
  elseif iscell(st1) % multiple match
    i2=false(size(i1));st3=Up.Stack(i1,2);
    i3=strncmp(st1,'#',1); % regexp entries
    i2(ismember(lower(st3),lower(st1(~i3))))=true;
    i2(any(cell2mat( cellfun(@(x) not(cellfun(@isempty,...
     regexp(st3,x(2:end),'once') )),st1(i3),'uni',0) ) ,2))=true;    
%     for j1=1:length(st1)
%       if st1{j1}(1)=='#';i2(~cellfun('isempty',regexp(st3,st1{j1}(2:end))))=1;
%       else; i2(strcmpi(st1{j1},st3))=1;
%       end
%     end 
    i2=find(i2);
  elseif st1(1)=='#';
   i2=find(~cellfun('isempty',regexp(r1(i1,2),st1(2:end))));
  else;i2=find(strcmpi(st1,r1(i1,2)));
  end
  if isempty(i2); val=cell(0,3); i1=[];
  else;i1=i1(i2);val=r1(i1,:);
  end
else
  val=r1(i1,:);
end

if nargin==4 
 if size(val,1)==1; val=val{1,3};
   if isequal(GiveData,'hdf');eval('val=sdthdf(''HdfGetData'',val);');end
 elseif isempty(val); val=[];
 elseif strcmpi(GiveData,'multi')
     warning('OpenFEM:Multicase','Returning sub-stack rather than value');
    if isequal(GiveData,'hdf');eval('val=sdthdf(''HdfGetData'',val);');end
 else;    error('Cannot return value for multiple matches');
 end
end
