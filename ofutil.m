function out=ofutil(varargin);

% OpenFEM developper utilities
%
% Accepted commands are
%
%  Path     : checks path consistency with possible removal of SDT
%  mexall   : compiles all needed DLL on UNIX and all C-mex files on Windows
%             On windows, run mex('-setup') and ofutil('of_mk') for Fortran mex
%  of_mk    : compiles of_mk.c (see openfem/mex directory)
%           ofutil('of_mk',otherfiles,mexoptions) is possible for developpers
%           to add other files in compilation and with specific mex options
%           giving an example: ofutil('of_mk','myfile.c','-DMYSYMBOL=MYVALUE')
%  nopo2sd  : compiles nopo2sd.c (located in openfem/mex directory)
%  sp_util  : compiles sp_ufil.c
%  zip      : creates a zip archive of the OpenFEM library 
%  hevea    : generates documentation with HEVEA
%  latex    : generates documentation with LaTeX
%   You may need to adjust the following preferences
%   setpref('OpenFEM','TexPath',fullfile(fileparts(which('ofutil')),'tex'))
%   setpref('OpenFEM','LatexFcn','!pdflatex')
%
%  validate : validate compatibility with SciLab
%
% HP-UX start MATAB R13 : setenv LD_PRELOAD /usr/lib/libF90.sl;matlab

%

%	E. Balmes, F. Genot, ...
%       Copyright (c) 2001-2024 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       Use ofutil('cvs') for revision information

[CAM,Cam]=comstR(varargin{1},1);


%Path Target verifications
%#ok<*ASGLU,*NOSEM>

Hom=fileparts(which('ofutil')); 

assignin('base','Hom',Hom);
if ~ispref('OpenFEM','TexPath') 
 setpref('OpenFEM','TexPath',fullfile(fileparts([Hom '.tex']),'tex'));
end
if ~ispref('OpenFEM','LatexFcn') 
 setpref('OpenFEM','LatexFcn','!pdflatex ')
end

if ~isempty(strfind(Cam,'-mexopt'));
  [CAM,Cam,RunOpt.mexopts]=comstR('-mexopt',[-25 4],CAM,Cam);
end

%% ---------------------------------------------------------------
if comstR(Cam,'cvs')

 out='$Revision: 1.174 $  $Date: 2024/01/17 09:01:04 $';

%% ---------------------------------------------------------------
elseif comstR(Cam,'cd'); [CAM,Cam]=comstR(CAM,3); 

if comstR(Cam,'tex'); cd(getpref('OpenFEM','TexPath'));
elseif comstR(Cam,'h'); cd(Hom);
end

% ---------------------------------------------------------------
elseif comstR(Cam,'cleana64')

pw0=pwd;
cd(fullfile(fileparts(which('ofutil')),'src'));
a=dir('*.c');
for j1=1:length(a)
 fid=fopen(a(j1).name);st=fscanf(fid,'%c');fclose(fid);
 st=strrep(st,'integer ','int32 ');
 st1=textscan(st,'%s','delimiter','\n'); st1=st1{1};
 in1=ismember(st1,'#include "../mex/of_def.h"');
 if sum(in1)>1||1
  in1=find(in1); st1(in1(2:end))=[]; st=sprintf('%s\n',st1{:});
 elseif ~any(in1)
  st=strrep(st,'#include "f2c.h"', ...
   sprintf('#include "f2c.h"\n#include "../mex/of_def.h"'));
 end
 fid=fopen(a(j1).name,'w');fprintf(fid,'%c',st);fclose(fid);
end
cd(fullfile(fileparts(which('ofutil')),'mex'));
a=dir('*.c');
for j1=1:length(a)
 fid=fopen(a(j1).name);st=fscanf(fid,'%c');fclose(fid);
 st=strrep(st,'integer ','int32 ');
 fid=fopen(a(j1).name,'w');fprintf(fid,'%c',st);fclose(fid);
end
cd(pw0);

%% #of_mk ---------------------------------------------------------------
elseif comstR(Cam,'of_mk') || comstR(Cam,'of_mk -move') 
% ofutil('of_mk -omp')
% Revised version where FORTRAN is converted with F2C

if comstR(Cam,'of_mk -move')
 Cam='of_mk';CAM='of_mk';
end
i1=strfind(Cam,'-newo');
if ~isempty(i1) % renew.o or .obbj files
 CAM(i1+(0:4))=''; Cam=lower(CAM); RunOpt.newo=1;
else; RunOpt.newo=0;
end

carg=2;
if carg<=nargin; otherfiles=varargin{carg};carg=carg+1;else; otherfiles='';end
of=[];

pw0=feval(ofutil('@MexTarget'),fileparts(which('ofutil')),'of_mk');

if carg<=nargin
 RunOpt.MexSwitch=varargin{carg};carg=carg+1;
elseif isfield(of,'MEXOPT'); RunOpt.MexSwitch=of.MEXOPT;
elseif isunix; RunOpt.MexSwitch='-v -Isrc';
else; RunOpt.MexSwitch='-v';
end

CAM(1:5)='';Cam(1:5)=''; 
[RunOpt,CAM,Cam]=VerSwitch(struct('MexSwitch',RunOpt.MexSwitch,'newo',RunOpt.newo), ...
    'blas',[],CAM,Cam);

 objectfiles = {'tab0d.o' ,'tab1d.o' ,'tab2d.o' ,'taba0d.o' ,'taba1d.o' ,'taba2d.o' ,'taba6d.o' ,'taba8d.o' ,'tabaxd.o' ,'tabszd.o' ,'ab0d.o' ,'ab1d.o' ,'ab4d.o' ,'ab5d.o' ,'hookax.o' ,'plmasd.o' ,'plonad.o' ,'etr3p1d.o' ,'etm3p1d.o' ,'ets3p1d.o' ,'etc3p1d.o' ,'etr3q1d.o' ,'etm3q1d.o' ,'ets3q1d.o' ,'etc3q1d.o' ,'etm3q2c.o' ,'etr3q2c.o' ,'ets3q2c.o' ,'ec3c2c.o' ,'etc3q2c.o' ,'etr3r1d.o' ,'etm3r1d.o' ,'ets3r1d.o' ,'etc3r1d.o' ,'er3c2c.o' ,'em3c2c.o' ,'es3d2c.o' ,'es3c2c.o' ,'emaq2c.o' ,'dpaq2c.o' ,'etr3r2c.o' ,'etm3r2c.o' ,'ets3r2c.o' ,'etc3r2c.o' ,'etr3p2c.o' ,'etm3p2c.o' ,'ets3p2c.o' ,'etc3p2c.o' ,'fobase.o' ,'dcopy.o' ,'etm2p1d.o' ,'etr2p1d.o' ,'ets2p1d.o' ,'etc2p1d.o' ,'etm2p2c.o' ,'etr2p2c.o' ,'ets2p2c.o' ,'etc2p2c.o' ,'etm2q1d.o' ,'etr2q1d.o' ,'ets2q1d.o' ,'etc2q1d.o' ,'etm2q2c.o' ,'etr2q2c.o' ,'ets2q2c.o' ,'etc2q2c.o' ,'e1ap1d.o' ,'etrap1d.o' ,'etmap1d.o' ,'etsap1d.o' ,'etcap1d.o' ,'e2ap2c.o' ,'etrap2c.o' ,'etmap2c.o' ,'etsap2c.o' ,'etcap2c.o' ,'e1aq1c.o' ,'etmaq1d.o' ,'etraq1d.o' ,'etsaq1d.o' ,'etcaq1d.o' ,'e2aq2c.o' ,'etmaq2c.o' ,'etraq2c.o' ,'eraq2c.o' ,'etsaq2c.o' ,'etcaq2c.o' ,'etmdktp.o' ,'etrdktp.o' ,'etsdktp.o' ,'etcdktp.o' ,'etm5noe.o' ,'etr5noe.o' ,'ets5noe.o' ,'etc5noe.o' ,'etrmit4.o' ,'etsmit4.o' ,'bremit.o' ,'etsmitx.o' ,'inmit4.o' ,'mitlin.o' ,'fonmit.o' ,'melmit.o' ,'replo1.o' ,'pcq12d.o' ,'canoq1.o' ,'chan57mod.o','chan56.o' ,'ddot.o'};

if strcmp(mexext,'mexmac')||strcmp(mexext,'mexmaci')
 objectfiles=setdiff(objectfiles,{'ddot.o','dcopy.o'});
end
if exist(fullfile(Hom,mexext),'dir');SrcWd=fullfile(Hom,mexext);
 st1={'of_mk_subs.c','f2c.h'};cd(SrcWd);
 for j1=1:length(st1);
  try;copyfile(fullfile('..','src',st1{j1}),st1{j1});
  catch; eval(sprintf('!cp -v %s %s',fullfile('..','src',st1{j1}),st1{j1}));
  end
 end
elseif exist(fullfile(Hom,'src'),'dir');SrcWd=fullfile(Hom,'src');
else;SrcWd=fullfile(Hom,'..','src');
end

cd(SrcWd); a=dir('*.f');

% If the .c file does not exist create it with f2c
for j1=1:length(objectfiles)
 if isempty(strfind(SrcWd,mexext)) && ...
  ~exist(fullfile(SrcWd,[objectfiles{j1}(1:end-2) '.c']),'file')
   if isunix
     system(sprintf('f2c %s',[objectfiles{j1}(1:end-2) '.f']));
   else
     system(sprintf('c:/scilab/bin/f2c %s',[objectfiles{j1}(1:end-2) '.f']));
   end
 end
end

% If the .object file does not exist create it with mex
if isunix;ext='.o'; else;ext='.obj';end

for j1=1:length(objectfiles)
 if ~exist(fullfile(SrcWd,[objectfiles{j1}(1:end-2) ext]),'file') || RunOpt.newo
  if ~exist(fullfile(SrcWd,[objectfiles{j1}(1:end-2) '.c']),'file') 
   sto=fullfile(fileparts(SrcWd),'src',[objectfiles{j1}(1:end-2) '.c']);
  else; sto=[objectfiles{j1}(1:end-2) '.c'];
  end
   eval(sprintf('mex -c -DOF_SMALL %s %s',RunOpt.MexSwitch,sto));
 end
end
cd(SrcWd);
st2={'of_mk.c','of_mk_pre.c','of_xkx.c','of_time_interp.c','hyper.c', ...
    'of_basis.c'};

if isunix % This is for compilation on unix systems

 if ~isempty(strfind(SrcWd,mexext))
   objectfiles2 = strcat(fullfile('..',mexext),filesep,objectfiles);
 elseif ~isempty(strfind(SrcWd,'../src'))
   objectfiles2 = strcat(fullfile('..','..','src'),filesep,objectfiles);
 else
   objectfiles2 = strcat(fullfile('..','src'),filesep,objectfiles);
 end

% On UNIX one is not supposed to need linking to BLASLIB
 st2=sprintf(' %s ',st2{:});
 st=sprintf('mex  %s %s %s %s %s\n!mv of_mk.%s %s',RunOpt.MexSwitch,st2, ...
  fullfile(SrcWd,'of_mk_subs.c'),otherfiles,sprintf('%s ',objectfiles2{:}),mexext,pw0);

    cd(fullfile(Hom,'mex'));
   
else % Now for windows assuming that f2c is located in Scilab

% Compile of_mk.c (now in openfem/mex) and of_mk_subs.c LOCATED IN SRC

  objectfiles=strrep(objectfiles,'.o ','.obj ');
  
  if ~isempty(strfind(SrcWd,'../src'))||~isempty(strfind(SrcWd,'..\src'))
  st2=sprintf('../openfem/mex/%s ',st2{:})
  st=sprintf('mex -v  %s %s  ../src/of_mk_subs.c %s %s %s\n!mv of_mk.dll %s', ...
   RunOpt.MexSwitch,st2,otherfiles, ...
   sprintf('%sbj ',objectfiles{:}),ofutil('blaslib'),pw0);

  else
   cd(SrcWd)
   st2=sprintf('../mex/%s ',st2{:});
   st=sprintf('mex -outdir ''%s'' %s %s %s of_mk_subs.c %s %s', ...
   pw0,RunOpt.MexSwitch,st2,otherfiles, ...
   sprintf('%sbj ',objectfiles{:}),ofutil('blaslib'));

  end

end
clear of_mk;assignin('base','st',st); 
try;evalin('base','eval(st)');
catch err;disp('Failed running the command below');disp(st);err.message,return;
end
MexTarget(Hom,'of_mk');
which('of_mk'); of_mk('cvs');disp(ans(:,2)) %#ok<NOANS>
if strcmpi(fileparts(which('of_mk')),pwd);cd ..;end % Do not finish in SDT3

%% ---------------------------------------------------------------
elseif comstR(Cam,'mexallopenfem')

ofutil('path')
ofutil('mexall')

%% #MexAll : generates all mex files into SDT3 -------------------------------
elseif comstR(Cam,'mexall')

st={'of_mk','nopo2sd','sp_util','of_time'};
tname=fullfile(tempdir,'mexall.log');diary(tname);
for j1=1:length(st);
 pw0=MexTarget(Hom);
 try; ofutil(st{j1});
 catch;
  st1='--------------------------------------------------------';
  fprintf('%s\n%s\n%s\n Error with %s\n %s\n%s\n%s\n',st1,st1,st1,st{j1},st1,st1,st1);
 end
end
if exist('sdeb','file'); feval('sdeb','ofPutZip');end
diary('off');fprintf('\n\nSaved log in %s\n\n',tname);

%% ---------------------------------------------------------------
elseif comstR(Cam,'of_time')

pw0=MexTarget(Hom,'of_time');
clear functions
[RunOpt,CAM,Cam]=VerSwitch([],'blas',[],CAM,Cam);

st=sprintf('mex %s "%s"',RunOpt.MexSwitch,fullfile(Hom,'mex','of_time.c'));
cd(pw0);disp(st);eval(st);
which('of_time');eval('of_time cvs');disp(ans(:,2)) %#ok<NOANS>
if strcmpi(fileparts(which('of_time')),pwd);cd ..;end % Do not finish in SDT3

%% ---------------------------------------------------------------
elseif comstR(Cam,'of_celt')

pw0=MexTarget(Hom,'of_celt');
clear functions
mex('-v',fullfile(Hom,'mex','of_celt.c'));
cd(pw0)

%% ---------------------------------------------------------------
elseif comstR(Cam,'nopo2sd')

pw0=feval(ofutil('@MexTarget'),Hom,'nopo2sd');
clear functions
mex('-v',fullfile(Hom,'mex','nopo2sd.c'));
cd(pw0)

if 1==2
 cd c:/balmes/sdt.old/openfem/test
 [nodes,elms]=readnopo('care','2D','v');
 [nodes,elms]=readnopo('q1','3D','v');
end

%% ---------------------------------------------------------------
elseif comstR(Cam,'medit')

pw0=MexTarget(Hom);
clear functions
mex('-v',fullfile(Hom,'mex','write_mesh.c'));
mex('-v',fullfile(Hom,'mex','write_bb.c'));
cd(pw0)

%% ---------------------------------------------------------------
elseif comstR(Cam,'sp_util')

[RunOpt,CAM,Cam]=VerSwitch([],[],[],CAM,Cam);

pw0=pwd;cd(MexTarget(Hom));
if exist('fe_caseg','file')~=2
 clear functions
 st=sprintf('mex %s "%s"',RunOpt.MexSwitch,fullfile(Hom,'mex','sp_util.c'));
 assignin('base','st',st);eval(st);
else
 warning('SDT provides it''s own version of sp_util, nothing DONE')
 return
end
cd(pw0)
which('sp_util');eval('st=sp_util(''cvs'')');disp(st(:,2))
if strcmpi(fileparts(which('sp_util')),pwd);cd ..;end % Do not finish in SDT3

%% ---------------------------------------------------------------
elseif comstR(Cam,'blaslib');[CAM,Cam]=comstR(CAM,8);

if comstR(Cam,'microsoft')
 if sdtdef('verm')>7.5
  setpref('OpenFEM','blaslib', ...
    fullfile(matlabroot,'extern','lib','win32','microsoft', ...
     'libmwblas.lib'));
 else;error('Not tested');
 end
 disp(getpref('OpenFEM','blaslib'));
 return;
elseif ispref('OpenFEM','blaslib')
  out=getpref('OpenFEM','blaslib');
elseif ~isunix&&strcmp(mexext,'mexw64')
 out=fullfile(matlabroot,'extern','lib','win64','microsoft','libmwlapack.lib');
elseif ~isunix
 if strcmpi(mexext,'mexw32')
     out=fullfile(matlabroot,'extern','lib','win32','microsoft', ...
     'libmwblas.lib');
    if exist(out,'file');return;end
 end 
 out=fullfile(matlabroot,'extern','lib','win32','microsoft', 'msvc60', ...
   'libmwlapack.lib');
else
 pw0=pwd;
 wd=fullfile(matlabroot,'bin');
 cd(wd); st=mexext;
 if strcmp(st,'mexglx')    
   wd=fullfile(matlabroot,'bin','glnx86','libmwlapack.so');
 elseif strcmp(st,'mexa64') 
   wd=fullfile(matlabroot,'bin','glnxa64','libmwlapack.so');
 elseif strcmp(st,'mexsg') 
   wd=fullfile(matlabroot,'bin','sgi','libmwlapack.so');
 elseif strcmp(st,'mexhp7') 
    wd=fullfile(matlabroot,'bin','hp700','libmwlapack.sl');
 elseif strcmp(st,'mexhpux') 
    wd=fullfile(matlabroot,'bin','hpux','libmwlapack.sl');
 elseif strcmp(st,'mexaxp') 
    wd=fullfile(matlabroot,'bin','axp','libmwlapack.so');
 elseif strcmp(st,'mexsol') 
    wd=fullfile(matlabroot,'bin','sol2','libmwlapack.so');
 elseif strcmp(st,'mexrs6') 
    wd=fullfile(matlabroot,'bin','ibm_rs','libmwlapack.a');
 elseif strcmp(st,'mexmac') 
    wd=fullfile(matlabroot,'bin','mac','libmwlapack.dylib');
 elseif strncmpi(st,'mexmac',6) 
    wd=fullfile(matlabroot,'bin',st(4:end),'libmwlapack.dylib');
 elseif strcmp(st,'mexw32') 
    wd=fullfile(matlabroot,'bin','win32','libmwlapack.dll');
 end
 out=wd;cd(pw0);
end
out=strcat('''',out,'''');

%% #zip : generate zip distribution ---------------------------------------
elseif comstR(Cam,'zip'); [CAM,Cam]=comstR(CAM,4);

 pw0=pwd;
 cd(tempdir);st=pwd;cd(fullfile(Hom,'..'));

 if isempty(CAM)
   st0=fullfile(st,'of2014.zip'); delete(st0)
   st1=sprintf(['!zip -rv %s openfem -x */*~ *.el *.dll */*.o*' ...
          ' mexa64 mexmaci mexmaci64 *.mtc* *.log *.out *CVS* */*/CVS .c* *#*'],st0);
   disp(st1);eval(st1)
 end
 try;eval(sprintf('!chmod 644 "%s"',st0));end
 cd(pw0);
 fprintf('\n\ncd "%s";scp of2014.zip sdtools.com:public_ftp/contrib/.\n\n',st) 

%% ---------------------------------------------------------------
elseif comstR(Cam,'lib');[CAM,Cam]=comstR(CAM,4);
    
 
if comstR(Cam,'of_mk')
  wd=pwd;
  cd(fullfile(fileparts(which('fe_mknl')),'mex'))
  %unloadlibrary of_mk
  [a,b]=loadlibrary('of_mk','of_EltConst.h');
  cd(wd);
  %libfunctions of_mk
  %libfunctionsview of_mk
elseif comstR(Cam,'ec')
  
 r3=varargin{2};
 if ~libisloaded('of_mk'); ofutil('libof_mk');end
 %z=libpointer('EltConst');
 r1=libstruct('EltConst'); val=int32(zeros(1,4));
 val=calllib('of_mk','GetpEC',r1,val);
 [r2,val]=calllib('of_mk','SetpEC',r1,r3);
 out=r1;

end
    
%% ---------------------------------------------------------------
elseif comstR(Cam,'validate'); [CAM,Cam]=comstR(CAM,4);

if isunix; Shell='unix'; else Shell='dos';end
pw0=pwd; 
cd(fileparts(which('ofutil')));

% Do not use persistent variables
[a,b]=feval(Shell,'grep persistent *.m */*.m');
if ~isempty(b)
 fprintf('----------------------------------------------------------\n');
 fprintf('persistent is used in :\n%s',b)
end

% Do not use evalin
[a,b]=feval(Shell,'grep evalin *.m */*.m');
if ~isempty(b)
 fprintf('----------------------------------------------------------\n');
  fprintf('evalin is used in :\n%s',b)
end

% Do not use the same function name more than once
[a,b]=feval(Shell,'grep function *.m sdt3/*.m @ofact/*.m');
ind=[0 find(b==10)];

st={}; 
for j1=1:length(ind)-1
 st1=b(ind(j1)+1:ind(j1+1)-1); if ~isempty(strfind(st1,'ofutil')); st1='';end
 st1=regexprep(st1,'(.*).m:','');
 i1=find(st1=='%'); if ~isempty(i1); st1(i1:end)='';end
 st1=strrep(st1,'_function','');
 st1=regexprep(st1,'''(.*)''','');
 if  isempty(strfind(st1,'function')); st1='';end
 st1=regexprep(st1,'%(.*)','');
 st1=regexprep(st1,'function \[(.*)\]=','');
 st1=regexprep(st1,'function (.*)=','');
 st1=regexprep(st1,'function','');
 st1=regexprep(st1,'\((.*)\)','');
 st1=regexprep(st1,';','');st1=comstR(st1,1);
 if ~isempty(st1); st{end+1}=st1;end 
end
fprintf('----------------------------------------------------------\n');
fprintf('Repeated functions called\n\n')
st=st(:); [st1,i1,i2]=unique(st);
i3=find(sparse(i2,1,1)>1);disp(char(st1(i3))); %#ok<FNDSB>

% Do use second arg for : ones, zeros, eyes
%                         sum, cumsum, prod, cumprod
st3={'ones(', 'zeros(', 'eyes(', 'sum(', 'prod('};

for j2=1:length(st3)
 [a,b]=feval(Shell,['grep ' st3{j2} ' *.m sdt3/*.m @ofact/*.m']);
 ind=[0 find(b==10)];
 st={}; 
 for j1=1:length(ind)-1
  st1=b(ind(j1)+1:ind(j1+1)-1);
  st2=st1;
  st1=st1(strfind(st1,st3{j2})+length(st3{j2}):end);
  % if min(findstr(st1,'(')) < min(findstr(st1,')'))    
  %   st1([min(findstr(st1,'(')) min(findstr(st1,')')) ])='';
  % end
  if min(strfind(st1,')')) < min(strfind(st1,','))    
   st{end+1}=st2;
  end
 end % j1
 fprintf('----------------------------------------------------------\n');
 fprintf('Bad number of input arguments for function %s : \n',st3{j2});
 st=st(:);disp(char(st(:)));
end %j2

% Do not use a(end,..) or (...,end)
%xxx

%




cd(pw0);

%% -----------------------------------------------------------------------
elseif comstR(Cam,'chmod'); [CAM,Cam]=comstR(CAM,6);
 
 a=fileparts(which(strcat('of_mk.',mexext)));
 if isunix
  st=sprintf('!chmod a+rx %s',fullfile(a,'*.mex*'));
 else
  st=sprintf('!cacls "%s" /E /G "%%username%%:R"',fullfile(a,'*.dll'));
 end
 disp(st);
 eval(st);

%% #Path ---------------------------------------------------------------
elseif comstR(Cam,'path'); [CAM,Cam]=comstR(CAM,5);

pw0=pwd;if isempty(CAM) && nargin==2; CAM=varargin{2};end
if ~isempty(CAM); cd(CAM);end
st=which('ofutil');[wd,st1]=fileparts(st);
try; eval('sdtcheck(''pathnone'');');end

if comstR(pwd, fileparts(which('sdtcheck')))
 error('You should not be in SDT home directory'); 
end

cd(wd);st = path; i1=strfind(st,[filesep 'sdt']);
i2 = [0 find(st==pathsep) length(st)+1];
for j1=1:length(i1)
  st1=st(i2(max(find(i2<i1(j1))))+1:i2(min(find(i2>i1(j1))))-1);
  if ~isempty(strfind(st1,path)); rmpath(st1); end
end

addpath(pwd);
st2={'sdt3','804','7.5','demos','test'};
for j1=1:length(st2)
 if j1==2 %804
  if sdtdef('verm')>=804
   st1=fullfile(pwd,st2{j1});addpath(st1);fprintf('%s\n',st1);
  end
 elseif j1==3%7.5
   if sdtdef('verm')<804
     st1=fullfile(pwd,st2{j1});addpath(st1);fprintf('%s\n',st1);   
   end
 elseif ~isempty(strfind(st2{j1},'6.0'))
  if exist(st2{j1},'dir')==7&&(comstR(version,'6')||comstR(version,'7')) % Matlab 6 7
    st1=fullfile(pwd,st2{j1});addpath(st1);fprintf('%s\n',st1);
  end
 elseif exist(st2{j1},'dir')==7
  st1=fullfile(pwd,st2{j1});addpath(st1);fprintf('%s\n',st1);
 end
end

i1=findall(0,'tag','rwarn');
if ~isempty(i1); delete(i1); end
% change ofact method if needed and possible
try;st=ofact('methodlist');ofact(st{1});end
cd(pw0);


%% ---------------------------------------------------------------
elseif comstR(Cam,'latex')

 [CAM,Cam,RunOpt.color]=comstr('-color',[-25 3],CAM,Cam);

 pw0=pwd;cd(getpref('OpenFEM','TexPath',fullfile(fileparts(which('ofutil')),'tex')));

 fprintf('Running latex manual in %s',pwd);
 st=fullfile(pwd,'..','..','tex');
 if exist(fullfile(st,'macros_tex.tex'),'file')
  !cp -v ../../tex/macros*.tex .
 end
 if ~RunOpt.color
  st=sprintf('%s manual',getpref('OpenFEM','LatexFcn','!pdflatex'));
  eval(st);
  try; eval('!makeindex manual');end
 else
  lat colortex
  lat _manual.tex
 end
 cd(pw0);

%% ---------------------------------------------------------------
elseif comstR(Cam,'hevea')

pw0=pwd;cd(getpref('OpenFEM','TexPath'));

if exist('manual.bbl','file')
 copyfile manual.bbl man_hevea.bbl
end
eval('!hevea mathaccents.hva man_hevea.tex');

pw2=fullfile(fileparts([getpref('OpenFEM','TexPath') '.tex']),'html');
if ~exist(pw2,'dir') 
 error(['You must create this directory : ' pw2]);
end


fid=fopen('man_hevea.html','r'); st=fscanf(fid,'%c'); fclose(fid);
delete man_hevea.html
st=strrep(st,'NAME="#','NAME="');
st=strrep(st,'<H1>Index</H1>','<A Name="Index"></a><H1>Index</H1>');
st=strrep(st,'<H1>References</H1>','<A Name="ref"></a><H1>References</H1>');
st2=fullfile(pw2,'man_hevea.html');if exist(st2,'file'); delete(st2);end
fid=fopen(st2,'w'); fprintf(fid,'%c',st); fclose(fid);

cd(pw2)
eval('!hacha -tocbis man_hevea.html');
delete man_hevea.html
eval('!cp ../tex/plots/*.gif .');
fprintf('Plots copied\n');
cd(pw0);

% ---------------------------------------------------------------
elseif comstR(Cam,'@');out=eval(CAM);
else;error('%s unknown',CAM);
end


%% #MexTarget ----------------------------------------------------------------
function pw0=MexTarget(Hom,Fcn) %#ok<INUSD>

pw0=pwd;
Vers=version; Vers=Vers(1:3); if Vers(1)=='6'||Vers(1)=='7'; Vers='6.0';end
st=fileparts(which('feplot'));

if exist('sd','file'); 
 eval('of=sd(''mexTarget'',Fcn);');
 if ischar(of);of=struct('MexTarget',of);end %#ok<NODEF>
 assignin('caller','of',of);
else
 if sdtdef('verm')>804&&exist(fullfile(Hom,'804'),'dir')
   of.MexTarget=fullfile(Hom,'804');
 elseif exist(fullfile(Hom,'sdt3'),'dir')
   of.MexTarget=fullfile(Hom,'sdt3');
 else;
  try('of.MexTarget=sd(''mextarget'',Fcn)');
  catch;  %#ok<CTCH>
   if exist(fullfile(st,Vers),'dir'); of.MexTarget=fullfile(st,Vers);
   else; error('of.MexTarget was not properly defined');
   end
  end
 end
end % if function sd exists

%if ispref('OpenFEM','MexTarget')
%  of.MexTarget=getpref('OpenFEM','MexTarget');
%end
fprintf('%%---------------\n%%---------------\n%%---------------\n');
fprintf('Target mex directory : %s\n',of.MexTarget)
fprintf('%%---------------\n%%---------------\n%%---------------\n');
pw0=of.MexTarget;

%% #OmpSwitch ----------------------------------------------------------------
function st=OmpSwitch(Cam,st);

i1=strfind(Cam,'-omp');
%[CAM,Cam,i1]=comstR('-omp',[-25 3],CAM,Cam);
if i1;
    st=[st ' -DOFOMP '];
    if ~isunix; st=[st ' COMPFLAGS="$COMPFLAGS /openmp "' ];
    else; % Verified that -fopenmp is also valid for MAC.
      st1=getenv('MATLAB_ARCH'); if isempty(st1);st1=getenv('ARCH');end
      wd=fullfile(matlabroot,'sys','os',st1);
      %wd='/opt/intel/clck/2017.1.016/provider/share/common/lib/intel64'
      if strncmpi(version,'7',1)
       st=[ sprintf(' -L"%s" -liomp5 ',wd) st ];
      else
       st=[' CFLAGS="\$CFLAGS -fopenmp " ' sprintf(' -L"%s" -liomp5 ',wd) st ];
      end
    end
end
i1=strfind(Cam,'-debug');
if i1;st=[' -g ' strrep([st ' '],' -O ',' ')];end


%% #VerSwitch ----------------------------------------------------------------
function [RunOpt,CAM,Cam]=VerSwitch(RunOpt,typ,of,CAM,Cam);
if nargin<4; CAM='-verbose'; Cam=lower(CAM);end
if exist('sd','file')&&isempty(RunOpt);eval('RunOpt=sd(''mexTarget'');');end

if isfield(RunOpt,'MEXOPT');RunOpt.MexSwitch=RunOpt.MEXOPT;end
if ~isfield(RunOpt,'MexSwitch');RunOpt.MexSwitch=' -v'; end
st=RunOpt.MexSwitch;
if isempty(strfind(st,'OSTYPE'));st=sprintf('%s -DOSTYPE%s ',st,mexext);end

if isempty(strfind(st,'-DFORMATLAB')); st=sprintf('%s -DFORMATLAB',st);end
if isempty(strfind(st,'MatlabVER'))
 i1=sdtdef('verm');if i1<10;i1=i1*10;end
 st=sprintf('%s -DMatlabVER=%.0f',st,i1);
end

RunOpt.ver=sdtdef('verm');
if isfield(RunOpt,'mexopts')&&~isempty(RunOpt.mexopts)
    st=sprintf('%s -f %s ',st,RunOpt.mexopts);
elseif ispref('OpenFEM','mexopts')&&~isempty(getpref('OpenFEM','mexopts'))
    st=sprintf('%s -f %s ',st,getpref('OpenFEM','mexopts'));
end
a=speye(1);a=whos('a');
if sdtdef('verm')>=904; 
 st=strrep(st,'-largeArrayDims',''); 
 if isempty(strfind(st,'R2018'));st=sprintf('-R2018a %s',st);end
elseif a.bytes==32; st=sprintf('-largeArrayDims %s',st);
end
i1=0; try; i1=exist('sp_util','file')&&sp_util('issdt'); end
if i1&&isempty(strfind(st,'SDT_ADD')); st=sprintf('-DSDT_ADD#yes %s',st); end

if nargin<2||isempty(strfind(typ,'blas'))||RunOpt.ver<7.5
    if RunOpt.ver>=7.1
     [st1,st2]=fileparts(ofutil('blaslib'));if st1(1)=='''';st1(1)='';end
     st=sprintf('-L"%s" -l%s %s',st1,st2(4:end),st);
    end
%elseif isfield(of,'BLASLIB')
%    st=sprintf('%s %s',of.BLASLIB,st);
elseif isunix
  % 7.5 -L/scratch_abelia/balmes/matlab75/bin/glnxa64 -lmwblas
  st=sprintf('-L"%s" -lmwblas %s', ...
       fullfile(matlabroot,'bin',getenv('ARCH')),st);
else % Windows
  st1=mexext;st1(1:4)='';
  st=sprintf('-L"%s" -lmwblas %s', ...
       fullfile(matlabroot,'extern','lib',['win' st1],'microsoft'),st);
end
i1=strfind(Cam,'-verbose');
if ~isempty(i1);st=[st ' -Dverbose ' ]; end

%[CAM,Cam,i1]=comstR('-debug',[-25 3],CAM,Cam);
i1=strfind(Cam,'-debug');
if ~i1;
 i1=strfind(Cam,'-g');
   %[CAM,Cam,i1]=comstR('-g',[-25 3],CAM,Cam);
end
if i1; st=[' -g ' strrep(st,' -O ',' ')];end
st=OmpSwitch(Cam,st);
RunOpt.MexSwitch=st;


%% #comstR ------------------------------------------------------------------
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
function [CAM,Cam,st]=comstR(CAM,ind,opt)

if nargin==2&&ischar(ind) %initial string comparison

   if ~isempty(ind) && length(CAM)>=length(ind)
    if all(CAM(1:length(ind))==ind); CAM=1;return; end
   end
   CAM=0; return;

elseif ischar(ind)&&nargin>2 % parameter extraction utility

   if length(CAM)<length(ind); ind=ind(1:length(CAM));end
   i1=min([find(comstR(CAM(1:length(ind)),-27)~=ind) length(ind)+1 length(CAM)+1]);
   CAM=CAM(i1:length(CAM));
   if ~isempty(CAM)
    i1=find(CAM==','); if ~isempty(i1); CAM(i1)=char(32*ones(size(i1)));end
    i1=find(CAM==''''); if ~isempty(i1); CAM(i1)=char(32*ones(size(i1)));end
   end
   if ischar(opt) % format given
     if comstR(opt,'%s'); [i1,i2,i3,i4]=sscanf(CAM,opt,1);
     else [i1,i2,i3,i4]=sscanf(CAM,opt);end
     Cam=CAM(i4:length(CAM));CAM=i1(:)';
     if ~isempty(Cam)
      ind = find(Cam~=32);ind = [min(ind):max(ind)];
      if ~isempty(ind); Cam = Cam(ind); else Cam=[]; end
     end
     if nargout==3; st=comstR(Cam,-27);
     elseif nargout==2 && ischar(CAM); Cam=comstR(CAM,-27); 
     end
   else Cam=comstR(CAM,-27); 
   end % return the string

elseif length(ind)==1 && ind(1)>0  % eliminate front and tail blanks

   CAM = CAM(ind:length(CAM));

   if ~isempty(CAM)
    ind = find(CAM~=32&CAM~=0);
    if ~isempty(ind)
       CAM = CAM(ind(1):ind(length(ind))); Cam=comstR(CAM,-27);
    else CAM='';Cam='';
    end
   else Cam=''; 
   end

end
