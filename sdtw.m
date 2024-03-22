function out=sdtw(varargin)

% SDT extended warning handling
%    sdtw('_nb',message)      for no backtrace
%    sdtw('_IsGui',message)   selects warndlg if gcbf is defined
%    sdtw('_err',message)     error graphic/or in command window
%    sdtw('_clip length ReturnString','format',values)
%    sdtw('_cbkl','Message') % Open to line callback message
%    sdtw('_ew','Message') % error for devs, warning else
%    sdtw('_ewt','Message') % displays an ErrorToDo for devs only
%    sdtw('_ewtN_date','Message') % for devs displays ErrorToDo then error
%                            asking to cleanup when N months passed after date
%    warn=sdtw('_off',{'warnid1',... }) % deact warns at once and get their state
%    sdtw('_set',struct('state',{'on'|'off',...},'identifier',{'warnid1',...}); %set states at once
%
%
%    additional arguments handled using sprintf

%	Etienne Balmes, Guillaume Vermot des Roches
%       Copyright (c) 2001-2024 by SDTools and INRIA, All Rights Reserved
%       Use under OpenFEM trademark.html license and LGPL.txt library license

%#ok<*NASGU,*ASGLU,*NOSEM>

if nargin==1 && comstr(varargin{1},'cvs')
 out='$Revision: 1.49 $  $Date: 2024/02/06 11:40:53 $'; return;
end
if nargin==0; help sdtw; return; end
CAM=varargin{1}; carg=2; Cam=lower(CAM);

s=''; grw='';
st1=version;i1=str2double(st1(1:3));
[CAM,Cam,ref]=comstr('-ref',[-25 3],CAM,Cam);

if comstr(CAM,'_nb')
 %% #_nb -1
  if i1<=6.2; [s,f]=warning;
  else;s=warning('query','backtrace'); warning('off','backtrace');
  end
  %warning('on'); XXX why should this be necessary ?
  CAM=varargin{carg};carg=carg+1;

% return a potentially clipped list
elseif comstr(Cam,'_clip')
 %% #_clip -1
  opt=comstr(CAM(6:end),[-1 72 0]);
  st=varargin{carg};carg=carg+1;
  val=varargin{carg};carg=carg+1;
  out={}; i1=0;j1=1;
  if isa(val,'cell')
   while j1<=length(val)
    st1=sprintf(st,val{j1}); i1=i1+length(st1);
    if i1<opt(1)/2 || i1>1e3 || (j1>length(val)-j1+1); out{end+1}=st1; %#ok<AGROW>
    else;out{end+1}='... '; j1=length(val)-j1+1;i1=i1+1e3; %#ok<AGROW>
    end
    j1=j1+1;
   end
  else
   while j1<=length(val)
    st1=sprintf(st,val(j1)); i1=i1+length(st1);
    if i1<opt(1)/2 || i1>1e3 || (j1>length(val)-j1+1); out{end+1}=st1;%#ok<AGROW>
    else;out{end+1}='... '; j1=length(val)-j1+1;i1=i1+1e3;%#ok<AGROW>
    end
    j1=j1+1;
   end
  end
  if opt(2)==1; out=[out{:}];end % return a string
  return;
elseif comstr(Cam,'_isgui')
 %% #isGui -1
  if ~isempty(get(0,'CallbackObject')); grw='Warning'; else; grw=''; end
  if ~isempty(findobj(0,'type','figure','tag','test_sdt'));grw='';end
  CAM=varargin{carg};carg=carg+1;
elseif comstr(CAM,'_gr')
 %% #_gr -1
  grw=comstr(CAM,4); if isempty(grw); grw='warning'; end
  CAM=varargin{carg};carg=carg+1;
elseif comstr(Cam,'_err'); 
 %% #_err -1
 grw='error'; CAM=varargin{carg};carg=carg+1;
elseif comstr(Cam,'_ew');[CAM,Cam]=comstr(CAM,4); % 
 %% #_ewt  sdtweb _bp sdtw displayet %error for developpers -1
 if comstr(Cam,'t')
  if ~sp_util('issdt'); grw='warning'; CAM=varargin{carg};carg=carg+1;
  elseif sdtkey('isdev');
   if length(Cam)>1 % allow review timeout
    r1=regexp(Cam,'^t([0-9].*)_(.*)$','tokens'); r1=r1{1};
    r2=datenum(r1{2}); r1=str2double(r1{1});
    if (now-r2)/30.438>r1; dbstack, error('Need to cleanup %s',varargin{carg}); end
   end
   grw='displayet'; CAM=varargin{carg};carg=carg+1;
   dbstack;
  else; return;
  end
 else
  %% #_ew -1
  if sdtkey('isdev'); grw='error'; else;grw='warning';end
  CAM=varargin{carg};carg=carg+1;
 end
 
elseif comstr(CAM,'_cbk') % generate a clickable callback
 %% #_cbk -1
 st1='';
 if  comstr(CAM,'_cbkl') % callback line
  r1=dbstack;r1=r1(2);st=sprintf('sdtweb(''%s#%i'')',r1.file,r1.line);
  if carg<=nargin; st1=varargin{carg};end
 elseif carg<nargin;st=sprintf(varargin{carg:end});
 else;st=varargin{carg};
 end
 fprintf('\n<a href="matlab:%s">%s</a>  %% %s\n',st,st,st1);return;
 
elseif comstr(CAM,'_off')
 %% #_off: deactivate several warnings at once and provide their initial state
 st=varargin{carg}; carg=carg+1;
 if ~iscell(st); st={st}; end
 warn=cellfun(@(x)warning('query',x),st,'uni',0);
 warn=cat(1,warn{:});
 cellfun(@(x)warning('off',x),{warn.identifier});
 if nargout>0; out=warn; end
 return
 
elseif comstr(CAM,'_set')
 %% #_set: set statuses of several warnings at once using a multi dim struct
 warn=varargin{carg}; carg=carg+1; % warn format must be like the output of _off
 cellfun(@warning,{warn.state},{warn.identifier});
 return

elseif comstr(Cam,'_lerr')
 %% #_lerr: last error reprint
 sts='<a href="matlab:matlab.internal.language.introspective.errorDocCallback(''%s'',''%s'',%i)" style="font-weight:bold">%s</a>';
 stl='<a href="matlab: opentoline(''%s'',%i,0)">line %i</a>';

 if carg>nargin;err=lasterror; %#ok<LERR> % MException.last can only be called in a command
 else; err=varargin{carg}; carg=carg+1;  
 end
 st=err.getReport;
 stn=err.stack(1).name; sti=err.stack(1).line; stf=err.stack(1).file;
 i1=strfind(st,'</a>'); i1=i1(1)+4;
 st2=sprintf('  %s:\n\nError using %s (%s)\n',...
  err.message,sprintf(sts,stn,stf,sti,stn),sprintf(stl,stf,sti,sti));
 st=[st2 st(i1+1:end) sprintf('\n')];
 if nargout>0; out=st; 
 else; disp(st);%fprintf(2,st);
 end

 return

end 
%% #output -1
mid=[];
if carg<=nargin % sprintf input style or MSGID
 r1=varargin(carg:end); % For ref, cleanup \ not to be interpreted by warning
 if strncmpi(r1{end},'MSGID:',6); mid=r1{end}(7:end); r1(end)=[]; end
 if ref; i1=cellfun(@ischar,r1); r1(i1)=strrep(r1(i1),'\','\\'); end
 CAM=sprintf(CAM,r1{:}); 
end
if strcmp(grw,'error');
  if sp_util('issdt')&&sdtdef('isinteractive'); errordlg(CAM,'SDT:Error','modal'); drawnow
  else; error(CAM);
  end
elseif strcmp(grw,'displayet'); 
 fprintf(1,'ErrorToDo : %s\n',CAM);% sdtweb _bp sdtw ErrorTodo
elseif ~isempty(grw);
 if sp_util('issdt')&&sdtdef('isinteractive'); warndlg(CAM,grw);drawnow
 else; sdtw('_nb',CAM);
 end
elseif ref; warning('SDTWarning:nb',CAM); 
else; % Cannot be error (eb)
 if ~isempty(mid); warning(mid,CAM); % allow MSGID
 else; warning(CAM);
 end
end
if ~isempty(s); warning(s); end
