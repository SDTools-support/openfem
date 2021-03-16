function commode(FunName,Cam)

%COMMODE prompt and script modes for evaluation of text commands
%
%	Syntax: commode NameOfComFuntion
%		commode('NameOfComFuntion','CommandString');
%
%	COMMODE NameOfComFuntion (a single argument) enters a command prompt
%	  mode allowing keyboard entries of commands. Existing command
%	  functions are iicom, idcom, fecom, ... In the prompt mode, the MATLAB
%	  prompt '>>' is replaced by a 'NameOfComFcn>' prompt. Different
%	  command functions can be chained using a comma ',' as separator.
%	  IDCOM for example enters COMMODE using commode('idcom,iicom');
%
%	CommandStrings are text commands recognized by the command functions.
%	  COMMODE allows to chain commands by separating them with
%	  a semi-column ';'. For example commode('iicom','cax1;abs;cax2;pha')
%	  will display FRF magnitude and phase in the first two axes of the
%	  GUI drawing figure.
%
%	Most command functions will parse a command starting by ';' using
%	  COMMODE. Thus commode('iicom','cax1;abs') and iicom(';cax1;abs')
%	  are equivalent.
%
%	Reserved commands not passed to the command functions are
%
%	q or quit	exit COMMODE but not MATLAB
%	%comments	parts of commands after a % are eliminated as comments
%	script FileName reads FileName line by line and executes the lines as
%	                command chains.
%
%	See also COMSTR, IICOM, FECOM, FEMESH, ...

%       Etienne Balmes  02/02/94, 10/28/96
%       Copyright (c) 1990-2004 by SDTools,
%       All Rights Reserved.

if nargin==0 error('No function name specified');    end
iFun=[0 find(FunName==',') length(FunName)+1];
out = 0;

if nargin<2
  for j1=1:length(iFun)-1
    if exist(FunName(iFun(j1)+1:iFun(j1+1)-1))~=2
      error([FunName(iFun(j1)+1:iFun(j1+1)-1) 'is not on your path']);
    end
  end
  Cam='';FileIn=0;
elseif comstr(Cam,'script')
  [FileIn,message]=fopen(comstr(Cam,7),'r');
  if ~isempty(message)   disp(['commode command: ' Cam]);error(message)
  else Cam=[]; end
else   FileIn=-1;  end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
while out == 0

 if isempty(Cam)&FileIn==0 drawnow;Cam = input([FunName '> '],'s');
 elseif isempty(Cam)&FileIn>0  Cam=fgetl(FileIn); if Cam==-1 Cam=';q'; end;end

 if ~isempty(Cam)
   i1 = find(Cam=='%'); if ~isempty(i1) Cam=Cam(1:i1(1)-1); end
   ind = [0 find(Cam==';'|Cam==10|Cam==13)];
   in2 = find(Cam=='[');in3 = find(Cam==']');
 else ind=0; i1=[]; in2=[];in3=[];end

 % matching brackets - - - - - - -
 while ~isempty(in2)
  while ~isempty(in3)
     i1 = max(find(in2<in3(1)));
     if isempty(i1) disp('commode: Mismatched brackets');Cam=[];in3=[];
     else
       ind = ind(find(ind<in2(i1)|ind>in3(1)));
       in2=in2([1:i1-1 i1+1:length(in2)]);in3=in3(2:length(in3));
     end
   end %isempty(in3) (closing brackets)
   if ~isempty(in2)
      Cam1='';while isempty(Cam1)
        if FileIn==0 Cam1 = input(['> '],'s');
	elseif FileIn>0  Cam1=fgetl(FileIn);if Cam1==-1 Cam1=';q'; end; end
	% commented lines - - - - - - -
        i1 = find(Cam1=='%');
        if ~isempty(i1) Cam1=Cam1(1:i1(1)-1); end
        Cam1 = deblank(fliplr(deblank(fliplr(Cam1))));
      end % isempty(Cam1)
      if strfind(Cam,'...')==length(Cam)-2 % ... line continuation
           Cam  = [Cam(1:length(Cam)-3) ',' Cam1];
      else Cam  = [Cam ';' Cam1]; end

      ind = [0 find(Cam==';')];in2 = find(Cam=='[');in3 = find(Cam==']');
   end
 end

 % loop for command evaluation - - - - - - - - - - - - - - - - - - - - - - -
 if max(ind)~=length(Cam) ind = [ind length(Cam)]; end

 for j1=1:length(ind)-1
  Cam0 = comstr(Cam(ind(j1)+1:ind(j1+1)),1);
  % double occurences of '
  i1=find(Cam0=='''');while ~isempty(i1)
    Cam0=[Cam0(1:i1(1)) '''' Cam0(i1(1)+1:length(Cam0))];i1=i1(2:length(i1))+1;
  end

  if     isempty(Cam0)	Cam1=';';
  elseif Cam0(length(Cam0))==';'   Cam0=Cam0(1:length(Cam0)-1);Cam1=';';
  else                             Cam1=' ';end

  if isempty(Cam0)
  elseif strcmp(lower(Cam0),'q')|strcmp(lower(Cam0),'quit')  return
  elseif comstr(Cam0,'script')             commode(FunName,Cam0);
  else
    for j2=1:length(iFun)-1 %loop on UI command functions
      if iFun(j2)+1<iFun(j2+1)-1
        eval(['out1=' FunName(iFun(j2)+1:iFun(j2+1)-1) ...
            '(''' Cam0 ''',''' Cam1 ''');'],'disp(lasterr);out1=''done'';');
      else out1='done';eval([Cam0 Cam1],'disp(lasterr);out1=''done'';');
      end
      if ~strcmp(out1,'unknown') break; end
    end
    if strcmp(lower(out1),'unknown') disp(['unknown command: ' Cam0]);end
  end

 end %loop on j1

Cam=''; if FileIn==-1 out=1; end
end %loop on commands - - - - - - - - - - - - - - - - - - - - - - - - - - - -




