function out=iigui(names,command)

%IIGUI	Advanced tricks linked to variable handling.
%
%	User calls to this function are obsolete, and simply replaced by
%       calls to iiplot (without arguments iiplot, initializes the proper
%	global variables)
%
%   eval(iigui({'v1','v2'},'SetInBaseC')) % set with comments
%    GetInBase, GetInCaller, SetInCaller are also available
%
%   eval(iigui({'fun(args)',nargout},'OutReDir')) % redirect to other fun
%   eval(iigui(model,'MoveFromCaller')) % move variable from caller workspace
%   eval(iigui({'n1','n2'},'MoveFromCaller')) % move variables from caller 


%       Etienne Balmes, Jean-Philippe Bianchi
%       Copyright (c) 1990-2014 SDTools, All Rights Reserved

%#ok<*NOSEM>
if nargin==0 
 names={'IDopt','IIpo','IIpo1','IIres','IIres1','XF', ...
        'XFdof','IIxf','IIxe','IIxh','IIxi','IIw'};
elseif nargin==2
  out=[];
  if strcmp(command,'emptyinbase') % #EmptyInBase - - - - - - - - - - - - - - 
    out=1;
    if evalin('base',sprintf('~exist(''%s'',''var'')',names))
    elseif evalin('base',sprintf('~isempty(%s)',names)); out=0;
    end
  elseif strncmpi(command,'setinbase',9) % #SetInBase - - - - - - - - - - - - 
     names=names(:)';
     if length(command)>9;
         fprintf('\nSetting base variables %s\n',strutil(names,-30))
     end
     names=names([1 1],:);
     out=sprintf('assignin(''base'',''%s'',%s);',names{:});
  elseif strncmpi(command,'setincaller',11) % #SetInCaller - - - - - - - - - - 
     names=names(:)';
     if length(command)>11
         fprintf('\nSetting caller variables %s\n',strutil(names,-30))
     end
     names=names([1 1],:);
     out=sprintf('assignin(''caller'',''%s'',%s);',names{:});
  elseif strcmpi(command,'getincaller') % #GetInCaller - - - - - - - - - - - -
     names=names(:)';names=names([1 1],:);
     out=sprintf('%s=evalin(''caller'',''%s'');',names{:});
  elseif strncmpi(command,'getinbase',9) % #GetInBase - - - - - - - - - - - - - -
     names=names(:)';names=names([1 1],:);
     out=sprintf('%s=evalin(''base'',''%s'');',names{:});
  elseif strcmp(command,'caller');  % #Caller - - - - - - - - - - - - - - - -
  elseif strcmpi(command,'movefromcaller') % #MoveFromCaller - - - - - - - - -
    st1=inputname(1);
    if isempty(st1) % Cell array of names for move
     names=names(:)';
     names=names([1  1 1],:);
     out=sprintf(['%s=evalin(''caller'',''%s'');' ...
        'evalin(''caller'',''clear %s'');'],names{:});
    else
     out=sprintf(['%s=evalin(''caller'',''%s'');' ...
        'evalin(''caller'',''clear %s'')'],inputname(1),names,names);
    end
  elseif strcmpi(command,'outredir')  % #OutReDir - - - - - - - - - - - - - -
    if length(names)>1; nout=names{2}; else nout=0; end
    if nout==0; st0='';
    else
     st0='out';
     if names{2}>=2 
      st0=sprintf('out,%s',sprintf('out%i,',1:names{2}-1)); st0(end)=''; 
     end
     st0=sprintf('[%s]=',st0);
    end
    if any(names{1}=='(');out=sprintf('%s%s;',st0,names{1});
    else;out=sprintf('%s%s(varargin{:});',st0,names{1}); 
    end
  else error('Not supported');
  end
 if ~isempty(out); return; end
%% ---------------------------------------------------------------------------  
elseif comstr(names,'cvs')
 out='$Revision: 1.6 $  $Date: 2015/02/04 14:36:26 $'; return;
elseif ischar(names); names = {names};
elseif isa(names,'cell') && ischar(names{1})
elseif isa(names,'sdth') && strcmp(names.type,'IDopt') && ~isempty(names.opt)
 out=names; cf=get(out.opt(1),'userdata');out.opt=[];
 out.data=cf.IDopt.data;return;
else out=names; return
end
if nargin<2; command='';end
if isequal(names,{'XF'})&&sdtkey('isdev'); error('Not valid');end

for j1=1:length(names)
 eval(['global ' names{j1}])

 if strcmp(command,'caller');  
 try; 
  i1 = evalin('caller',sprintf('exist(''%s'',''var'')',names{j1}));
  if i1; i1 = evalin('caller',names{j1});else i1=[];end

  if eval(sprintf('isequal(i1,%s)',names{j1}))
  % variable in caller workspace is correct, make sure it is global
   evalin('caller',sprintf('clear %s;global %s;',names{j1},names{j1}))
  elseif eval(sprintf('isempty(%s)',names{j1})) && ~isempty(i1)
  % global is empty and base is not
   evalin('caller',sprintf('clear %s;global %s;',names{j1},names{j1}))
   eval(sprintf('%s=i1;',names{j1}))
  else % Replace the global variable by the caller value
   eval(sprintf('%s=i1;',names{j1}))
   evalin('caller',sprintf('clear %s;global %s;',names{j1},names{j1}))
  end
 catch;  %#ok<CTCH>
   evalin('caller',sprintf('global %s;',names{j1}))
 end;
 else;  
  i1 = evalin('caller',sprintf('exist(''%s'',''var'')',names{j1}));
  if ~i1; evalin('caller',sprintf('global %s;',names{j1}));
  else; evalin('caller',sprintf('clear %s;global %s;',names{j1},names{j1})); end
 end

 try; 
  i1 = evalin('base',sprintf('exist(''%s'',''var'')',names{j1}));
  if i1; i1 = evalin('base',names{j1});else i1=[];end
  %  i1 = evalin('base',names{j1});

  if eval(sprintf('isequal(i1,%s)',names{j1}))
  % variable in base workspace is correct, make sure it is global
   evalin('base',sprintf('clear %s;global %s;',names{j1},names{j1}))
  elseif eval(sprintf('isempty(%s)',names{j1})) && ~isempty(i1)
  % global is empty and base is not
   evalin('base',sprintf('clear %s;global %s;',names{j1},names{j1}))
   eval(sprintf('%s=i1;',names{j1}))
  else
   evalin('base',sprintf('clear %s;global %s;',names{j1},names{j1}))
  end
 catch;  %#ok<CTCH>
   evalin('base',sprintf('global %s;',names{j1}))
 end

end

if nargout==1 && length(names)==1; out=evalin('base',names{1}); end

if nargin==0; comgui('gcf ii'); end
