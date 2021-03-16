function [CAM,Cam,st]=comstr(CAM,ind,opt,in4)

% Port of string manipulations routines common to OpenFEM and SDT

%	E. Balmes
%       Copyright (c) 2001-2008 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       $Revision: 1.13 $  $Date: 2017/02/24 08:13:04 $


if nargin==2&&ischar(ind) %initial string comparison

   if ~ischar(CAM); CAM=0;
   elseif ~isempty(ind) && length(CAM)>=length(ind)
    if all(CAM(1:length(ind))==ind); CAM=1;return; end
   end
   CAM=0; return;

elseif ischar(ind)&&nargin>2 % parameter extraction utility

   if length(CAM)<length(ind); ind=ind(1:length(CAM));end
   if ~isempty(ind)
    i1=min([find(lower(CAM(1:length(ind)))~=lower(ind)) ...
         length(ind)+1 length(CAM)+1]);
    CAM=CAM(i1:length(CAM));
   end
   if ~isempty(CAM)
    i1=find(CAM==','); if ~isempty(i1); CAM(i1)=char(32*ones(size(i1)));end
    i1=find(CAM==''''); if ~isempty(i1); CAM(i1)=char(32*ones(size(i1)));end
   end
   if ischar(opt) % format given
     if comstr(opt,'%s'); [i1,i2,i3,i4]=sscanf(CAM,opt,1);
     else [i1,i2,i3,i4]=sscanf(CAM,opt);end
     Cam=CAM(i4:length(CAM));CAM=i1(:)';
     if ~isempty(Cam)
      ind = find(Cam~=32);ind = [min(ind):max(ind)];
      if ~isempty(ind); Cam = Cam(ind); else Cam=[]; end
     end
     if nargout==3; st=lower(Cam);
     elseif nargout<=2 && ischar(CAM) 
      if strcmp(opt,'%c'); [CAM,Cam]=comstr(CAM,1);
      else Cam=lower(CAM); end
     end
   else Cam=lower(CAM); 
   end % return the string

elseif length(ind)==1 && ind(1)>0  % eliminate front and tail blanks

   CAM = CAM(ind:length(CAM));

   if ~isempty(CAM)
    ind = find(CAM~=32&CAM~=0);
    if ~isempty(ind)
       CAM = CAM(ind(1):ind(length(ind))); Cam=lower(CAM);
    else CAM='';Cam='';
    end
   else Cam=''; 
   end

elseif ind(1)==-1 % find numbers with error and default handling

   if ~isempty(CAM); i1 = find(CAM~=32&CAM~=10&CAM~=13);
     if ~isempty(i1); i1=min(i1):max(i1); end
   else i1=[]; 
   end
   if ~isempty(i1) && length(i1)<length(CAM); CAM=CAM(min(i1):max(i1)); end
   eval(['CAM=[' CAM '];'],'CAM=[];')

   if length(ind)>1
	i1 = length(CAM)+1:length(ind)-1;
	if ~isempty(i1); CAM(i1) = ind(i1+1); end	  
   end
elseif ind(1)==-7 % Platform dependent handling of running label

   if comstr(computer,'PC'); fprintf(1,'\n')
   else;fprintf(1,char(8*ones(1,length(CAM))));
   end
   fprintf(1,'%s',opt);
   CAM = opt;

elseif ind(1)==-8 % Extract string with matching {} 
  % [extract,CAM,Cam] 

  i1=find(CAM==opt(1)|CAM==opt(2));
  if isempty(i1); Cam=CAM;st=lower(Cam);CAM='';return; end
  i2=0; 
  for j1=i1
    if     CAM(j1)==opt(1); i2=i2+1; 
    elseif CAM(j1)==opt(2); i2=i2-1; end
    if i2<0; error('''%s'' contains non matching %s',CAM,opt);end
    if i2==0; break; end
  end
  Cam=CAM;Cam(i1(1):j1)=' '; st=lower(Cam); CAM=CAM(i1(1):j1);

elseif ind(1)==-27 % Extract string with matching {} 

  CAM=lower(CAM);

% THIS IS COPIED FROM SDT. ANY CHANGE SHOULD BE NOTIFIED TO SDTOOLS
% Remove Token - - - - - - -  - - - - - - - - - - - - - - - - - - - -
% out=comstr(token,-25,CAM,Cam)
elseif ind(1)==-25

  i1=strfind(in4,CAM);st=[];
  if isempty(i1)  % default returns
    if length(ind)<2
    elseif isempty(st)&&ind(2)==4; st=''; % empty text
    elseif ind(2)==2&&length(ind)>3; st=ind(4); 
    elseif ind(2)==3; st=0;
    end

  elseif length(ind)==1;  opt(i1+[0:length(CAM)-1])='';
  elseif ind(2)==1; % integer - - - - - - - - - - - - - -
   if length(ind)==2; 
    [st,i2,i3,i4]=sscanf(in4(i1+length(CAM):end),'%i');
    i2=i1+i4+length(CAM)-1;
    if length(opt)>i2&&opt(i2)==':'
     [st1,i3,i4,i5]=sscanf(in4(i2+1:end),'%i');
     if ~isempty(st1);st=[st:st1(1) st1(2:end)];i2=i2+i5;end
    end
   else;[st,i2,i3,i4]=sscanf(in4(i1+length(CAM):end),'%i',ind(3));
      i2=i1+i4+length(CAM)-1;
   end
   opt(i1:i2-1)='';

  elseif ind(2)==2; % double - - - - - - - - - - - - - -
   if length(ind)==2; 
    [st,i2,i3,i4]=sscanf(in4(i1+length(CAM):end),'%g');
   else;[st,i2,i3,i4]=sscanf(in4(i1+length(CAM):end),'%g',ind(3));
   end
   opt(i1+[0:i4+length(CAM)-2])='';
  elseif ind(2)==3; % present or no
   st=1;opt(i1+[0:length(CAM)-1])='';
  elseif ind(2)==4; % string
   if length(ind)==2; 
    [st,i2,i3,i4]=sscanf(opt(i1+length(CAM):end),'%s');
   else;[st,i2,i3,i4]=sscanf(opt(i1+length(CAM):end),'%s',ind(3));
   end
   opt(i1+[0:i4+length(CAM)-2])='';
  end
  [CAM,Cam]=comstr(opt,1);
% THIS IS COPIED FROM SDT. ANY CHANGE SHOULD BE NOTIFIED TO SDTOOLS


elseif ind(1)==-271 % Extract string with matching {} 

  CAM=upper(CAM);
% ------------------------------------------------------------------------
elseif ind(1)==-32 % char 2 string : limited to 6 char
    
i1=[1 256 65536 16777216 4294967296 1099511627776]';%sprintf('%.0f ',2.^([0:5]*8))
%if ~strcmp(comstr(abs('abcd')*i1(1:4),-100),'abcd');i1=flipud(i1);end
[st,i2,i3]=fopen(0);
if ~isempty(strfind(i3,'be'));i1=flipud(i1);end
if ischar(CAM);
    st=abs(CAM);
    if length(st)<length(i1); st(end+1:length(i1))=0;
    elseif length(st)>length(i1); error('String too long, not supported');end
    CAM=abs(st)*i1;
elseif iscell(CAM)
  out=zeros(size(CAM));for j1=1:length(out);out(j1)=comstr(CAM{j1},-32);end
  CAM=out;
else;
  out=cell(size(CAM));
  for j1=1:length(CAM)
   r1=rem(fix(double(CAM(j1))./[1 256 65536 16777216 4294967296]),256);
   out{j1}=comstr(char(r1),1);
  end
  if length(out)==1;CAM=out{1};else;CAM=out;end
end
% ------------------------------------------------------------------------
else ; error(sprintf(' %i unknown',ind(1))); 

end
