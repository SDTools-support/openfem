function ks=mtaucs_utils(varargin);


% Gateway function for inclusion of the TAUCS solver into the 
% OpenFEM ofact object
%
% http://www.tau.ac.il/~stoledo/taucs/
%
% TAUCS Version 1.0, November 29, 2001. 
% Copyright (c) 2001 by Sivan Toledo, Tel-Aviv Univesity, 
% stoledo@tau.ac.il. All Rights Reserved. 
%
%TAUCS License:
%
% Your use or distribution of TAUCS or any derivative code implies that you 
% agree to this License. 
% THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY EXPRESSED 
% OR IMPLIED. ANY USE IS AT YOUR OWN RISK. 
%
%Permission is hereby granted to use or copy this program, provided that the
% Copyright, this License, and the Availability of the original version is 
% retained on all copies. User documentation of any code that uses this code 
% or any derivative code must cite the Copyright, this License, the 
% Availability note, and "Used by permission." If this code or any derivative 
% code is accessible from within MATLAB, then typing "help taucs" must cite 
% the Copyright, and "type taucs" must also cite this License and the 
% Availability note. Permission to modify the code and to distribute 
% modified code is granted, provided the Copyright, this License, 
% and the Availability note are retained, and a notice that the code 
% was modified is included. This software is provided to you free of charge. 
%
%
% To modify factorization method use
% ofact('method mtaucs','FactBuild','ks=mtaucs_utils(''fact'',k,ks,[2 1 1]);')
% Methods are 1 LDLt 2 SNMF 3 SNLL 4 Symbolic 5 OOC

%       $Revision: 1.7 $  $Date: 2008/04/23 08:06:03 $


Cam=varargin{1};carg=2;
if comstr(Cam,'method')

   ks.name='mtaucs';
   ks.header='TAUCS sparse solver';
   ks.SymRenumber='';
   ks.FactBuild='ks=mtaucs_utils(''fact'',k,ks,[2 0 0]);';
   ks.Solve='q=mtaucs_utils(''solve'',k,full(b));';
   ks.Clear='mtaucs_utils(''clear'',ks);';
   ks.Available=exist('mtaucs','file')==3;
   ks.HandlesComplex=0;
   ks.param=[];

elseif comstr(Cam,'fact')

 k=varargin{carg};carg=carg+1;
 ki=varargin{carg};carg=carg+1;
 if carg>nargin opt=[1 1 size(k,1)/32  0.01*size(k,1)];
 else opt=varargin{carg};carg=carg+1;
 end
 st={'ldlt','snmf','snll','sym','ooc'};
 ks=ki;

 t=cputime;
 if ischar(k)  % Matrix stored in file
    i1=mtaucs('fact',k,opt);
    ks.method.name=sprintf('mtaucs',st{opt(1)});
 elseif opt(1)==15 % ooc
     k=tril(k);k=k+spalloc(size(k,1),size(k,2),0);
     ka=sp_util('sp2st',k);
     fid=fopen('kooc','w');
     fwrite(fid,[size(k) 2057],'int32');
     fwrite(fid,ka.jc,'int32');
     fwrite(fid,ka.ir,'int32');
     fwrite(fid,ka.pr,'float64');
     fclose(fid);
     fid=fopen('mat2taucs','w');fprintf(fid,'0');fclose(fid);
     st=which(sprintf('mtaucs%s.exe',mexext));
     if isempty(st) 
      st=fullfile(matlabroot,'bin','win32',sprintf('mtaucs%s.exe',mexext));
     end
     if ~exist(st,'file') error('ooc executable not found');end
     st=sprintf('!"%s" &',st); 
     i1=1; ks.method.name='taucs_ooc';
     ks.method.Solve=sprintf('q=mtaucs_utils(''solve'',k,full(b),%i);',opt(4));
     SendMessage('ReadAndFactor');
     delete('t2m');
     eval(st);
     t=clock; while ~exist('t2m','file') 
      pause(.1);if etime(clock,t)>30  error('Startup failure');end
     end
     WaitUntilReady;
 else
   if ~issparse(k) k = sparse(k); end
   i1=find(diag(k)==0);
   if ~isempty(i1) 
     if length(i1)<10
      error(sprintf('Zero terms on diagonal :%i\n',i1)); 
     else error(sprintf('%i zero terms on diagonal \n',length(i1))); end   
     return
   end
   if mtaucs>1.0001 
    i1=mtaucs('fact',k,opt);
   else
    k=tril(k);k=k+spalloc(size(k,1),size(k,2),0);i1=mtaucs('fact',k,opt);
   end
   ks.method.name=sprintf('mtaucs_%s',st{opt(1)});
 end

 if opt(2) fprintf('\nOptions %sm %.2f s',sprintf('%g ',opt),cputime-t); end

 ks.ty(2:3)=[i1 size(k,1)];
 ks.ind='mtaucs';

elseif comstr(Cam,'solve')

     k=varargin{2}; b=varargin{3};
     if issparse(b) b=full(b);end
     ks=zeros(size(b));
     for j1=1:size(b,2)
       if strcmp(k.method.name,'taucs_ooc')  % external program
        opt=varargin{4};
        fid=fopen('booc','w');
        fwrite(fid,size(b,1),'int32');
        fwrite(fid,b(:,j1),'float64');
        fclose(fid);
        SendMessage('Solve');  
        WaitUntilReady;
	fid=fopen('qooc','r'); N=fread(fid,1,'int32'); r1=fread(fid,N,'float64');
        fclose(fid); 
        if any(~isfinite(r1)) error(1);end
        ks(:,j1)=r1;
        SendMessage('Received');
       else % DLL
        ks(:,j1)=mtaucs('solve',k.ty(2),b(:,j1));
       end
     end

elseif comstr(Cam,'clear')

 ks=varargin{carg};carg=carg+1;

 if isempty(ks)  
   mtaucs('clear',-1);
 elseif isa(ks,'double') 
    try; mtaucs('clear',ks);catch; disp(lasterr);end
 elseif isfield(ks,'method')&isfield(ks.method,'name') & ...
   strcmp(ks.method.name,'taucs_ooc')
        fid=fopen('mat2taucs','w');fprintf(fid,'exit');fclose(fid);
 else
    try; mtaucs('clear',ks.ty(2)); catch;disp(lasterr);end
   if exist('mat2taucs')
        fid=fopen('mat2taucs','w');fprintf(fid,'exit');fclose(fid);
   end
 end

end

% -------------------------------------------------------------------------
function   WaitUntilReady(st);

persistent fid
if isempty(fid) fid=fopen('t2m','r'); end


fseek(fid,0,-1);st=fscanf(fid,'%s');
while ~comstr(st,'Waiting')
  pause(1);drawnow;
  fseek(fid,0,-1);st=fscanf(fid,'%s');
  if comstr(st,'failed') error(1);end
end
SendMessage('Received');
while ~strcmp(st,'ready')
  pause(1);drawnow;
  fseek(fid,0,-1);st=fscanf(fid,'%s');
  if comstr(st,'failed') error(1);end
end

% -------------------------------------------------------------------------
function  SendMessage(st);

 fid=fopen('mat2taucs','w');fprintf(fid,st);fclose(fid);
