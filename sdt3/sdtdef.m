function out=sdtdef(in1,in2,in3); %#ok<INUSD>

% Defaults handling routine (internal)
%
%      Syntax : sdtdef info
%               sdtdef('ConstantName',Value)
%               ...
%
%      cf            : current GUI figure storage
%      DefBackCol    : default background color
%      epsl          : tolerance on node coincidence (1e-6)
%      info          : info on currently stored defaults
%      kcelas        : default penalization stiffness (1e14)
%      view          : list of standard views and associated tags
%                      set with sdtdef('view tag',value)

%	Etienne Balmes
%       Copyright (c) 1990-2023 by SDTools,
%       All Rights Reserved.

persistent II
%#ok<*NOSEM,*TRYNC>
%mlock
if nargin>=2&&isequal(in1,'diag')&&ischar(in2) % early return of diag warning
 try; i1=sp_util('diag');if i1==0;return;end;end  
end

if nargin==0   % return stamp
 if ~isfield(II,'stamp'); II.stamp=0; end
 out=II.stamp;
elseif ~ischar(in1)
 if in1>0  % return 
  if ~isfield(II,'stamp'); II.stamp=0; end
  II.stamp=II.stamp+1;
  out=II.stamp;
 else;error('Not a supported case');
 end
else
    
    CAM=in1; Cam=lower(CAM);
st1={'-setpref','-safe'};
for j1=1:length(st1)
 i1=strfind(Cam,st1{j1}); RO.(st1{j1}(2:end))=0;
 if ~isempty(i1); RO.(st1{j1}(2:end))=1; CAM(i1:i1+length(st1{j1})-1)=''; end
end
in1=CAM; Cam=lower(CAM);
i1=strfind(Cam,'.'); if ~isempty(i1); Cam=Cam(1:i1-1); end % group access

if strcmp(Cam,'diag') % - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if ~isfield(II,'diag'); II.diag=0; end
  if nargin==2
   if ~ischar(in2); II.diag=in2(1);   % Set value 
   elseif ~isempty(strfind(in2,'diagoff')); II.diag=0;
   elseif ~isempty(strfind(in2,'diagon'));  II.diag=1;
   elseif ~strcmp(in2,'tog');  
    if II.diag>1; evalin('caller','dbstack');end
    if II.diag;  disp(in2); end % Display
   else; % Tog
    II.diag=rem(II.diag(1)+1,2);   ssetpref('SDT','diag',II.diag);
   end
  elseif ~isempty(strfind(Cam,'diagoff')); II.diag=0;
  elseif ~isempty(strfind(Cam,'diagon')); II.diag=1;
  elseif nargin>2
   fprintf('%s ',in2);fprintf('\n');
  else  
    out=II.diag; 
  end
  try;sp_util('diag',II.diag);end

elseif strcmp(Cam,'isdeployed') % deployement test- - - - - - - - - - - - -
  if ~isfield(II,'IsDeployed');
   try; II.IsDeployed=isdeployed; catch; II.IsDeployed=0;end 
  end
 out=II.IsDeployed;
elseif strcmp(Cam,'info') % display current values - - - - - - - - - - - - -
  if isempty(II)
   disp('No default defined');
  else;disp(II); 
  end

elseif strcmp(Cam,'deflen') % default length for various figures
  if ~isfield(II,'DefLen'); II.DefLen=1; end
  if nargin==2&&nargout==0&&length(in2)==1
    II.DefLen(in2(1))=in2(2);
  else
    if nargin<2||in2(1)>length(II.DefLen); in2(1)=1;end
    out=II.DefLen(in2(1));
  end
 
elseif strcmp(Cam,'defbackcol') % default background color - - - - - - - - - -
  if ~isfield(II,'DefBackCol')
    II.DefBackCol(1,:)=[1 1 1]-round(mean(get(0,'defaultaxescolor')));
  end
  if nargin==2; II.DefBackCol(in2(1),:)=[1 1 1]-in2(2); 
  else;out=II.DefBackCol(1,:); end

elseif strcmp(Cam,'epsl') % length tolerance - - - - - - - - - - - - -
  if ~isfield(II,'epsl'); II.epsl=1e-6; end
  if nargin==2; II.epsl=in2(1); end
  try; sp_util('epsl',II.epsl);end 
  out=II.epsl;

elseif strcmp(Cam,'kcelas') % default penalization - - - - - - - - - - - - -
  if ~isfield(II,'kcelas'); II.kcelas=1e14; end
  if nargin==2; II.kcelas=in2(1); end
  out=II.kcelas;
  
elseif strncmp(Cam,'openfem',7);CAM=in1;CAM(1:8)=''; % allow reload of openfem - - - - - - - - - - -
  if ~isfield(II,'OpenFEM'); II.OpenFEM=struct;end
  if nargin==2 % If provided, set in II and 
   if ~isfield(II.OpenFEM,CAM)||~isequal(II.OpenFEM.(CAM),in2)
      setpref('OpenFEM',CAM,in2);II.OpenFEM.(CAM)=in2;
   end
  elseif ~isfield(II.OpenFEM,CAM)
   II.OpenFEM.(CAM)=getpref('OpenFEM',CAM,in2);
  end
  out=II.OpenFEM.(CAM);

elseif strncmp(Cam,'nostop',6) % #nostop: stop with demosdt| interactivity marker
 if ~isempty(strfind(Cam,'-setpref')); doPref=1; else; doPref=0; end
  if strncmp(builtin('version'),'5',1)
    if ~isfield(II,'nostop'); II.nostop=0; end
    if nargin==2; II.nostop=in2;
    else;out=II.nostop; end
  elseif nargin==1 
   if isfield(II,'nostop'); out=II.nostop;
   elseif ~doPref; out=0; II.nostop=0; %II.nostop;
   else; out=getpref('SDT','nostop',0);  II.nostop=out;
   end
  elseif nargin==2;
   if doPref; setpref('SDT','nostop',in2); end
   out=in2; II.nostop=out;
  end
elseif strncmp(Cam,'verm',4) % version - - - - - - - - - - - - - - - -

 st=builtin('version');i1=find(st=='.');
 i1=[0 i1 length(st)+1]; i1=i1(1:3); out=0;
 for j1=1:length(i1)-1
  out=out+10^(2*(length(i1)-j1)-2)*str2double(st(i1(j1)+1:i1(j1+1)-1));
 end

elseif strncmp(Cam,'luthres',7) % luThres
 % lu pivoting threshold strategy changed for some cases
 if sdtdef('verm')>903; out=.1; else; out=0; end

elseif strcmp(Cam,'rwarn') % - - - - - - - - - - - - - - - - - - - - - - - - -
 out=findall(0,'type','figure','tag','rwarn');
elseif strcmp(Cam,'isinteractive') % - - - - - - - - - - - - - - - - - - - - -
out=1;
if usejava('jvm') && ~feature('ShowFigureWindows');out=0;end

elseif strcmp(Cam,'cvs') % - - - - - - - - - - - - - - - - - - - - - - - - - -
 out='$Revision: 1.21 $  $Date: 2024/01/19 09:56:06 $';
else;error('%s unknown',CAM);
end
end

