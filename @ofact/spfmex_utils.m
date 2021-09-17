function ks=spfmex_utils(varargin);

% Gateway function for inclusion of the SPFMEX solver into the 
% OpenFEM ofact object
%
% Standard Factor building method
%
%   ofact('spfmex','FactBuild', ...
%    'ks=spfmex_utils(''fact'',k,ks,[1 0 size(k,1)/32  0.01*size(k,1)]);');
%
% Call to the external SDT_Server. 
%  ofact('spfmex','FactBuild', ...
%   'ks=spfmex_utils(''servfact'',k,ks,[1 0 size(k,1)/32  0.01*size(k,1)]);');
%
%  To close the server down : pnet('closeall')
%
% You can set host and port using
%  setpref('SDT','ServerHost','localhost');
%  setpref('SDT','ServerPort',888);
%
%       $Revision: 1.47 $  $Date: 2021/09/15 15:14:53 $

Cam=varargin{1};carg=2;
%#ok<*NASGU,*ASGLU,*NOSEM>

%% #Solve --------------------------------------------------------------------
if comstr(Cam,'solve')

 k=varargin{carg};carg=carg+1;
 b=varargin{carg};carg=carg+1; 
 if size(b,1)~=size(k,2) ; error('RHS has a bad size');  
 elseif ~isreal(b) || ~spfmex('isreal',k.ty(2))
  if isempty(k.dinv) % complex numbering outside
    ks=spfmex('solve',k.ty(2),complex(full(b)));
  else % complex numbering here
    if ~isreal(b)
      ks(k.dinv,:)=spfmex('solve',k.ty(2),full(b(k.dinv,:)));
    else
      ks(k.dinv,:)=spfmex('solve',k.ty(2),complex(full(b(k.dinv,:))));
    end
  end
 elseif isempty(k.dinv)  % real numbered outside
  ks=spfmex('solve',k.ty(2),full(b));
 else % real numbering here
  if size(k.dinv,1)~=k.ty(3); error('Bad size of RHS'); end
  ks=zeros(size(b,1),size(b,2));
  %ks(k.dinv,:)=spfmex('solve',k.ty(2),full(b(k.dinv,:)));
  b=full(b(k.dinv,:)); spfmex('solve',k.ty(2),b,ks,int32(k.dinv-1));
 end

%% #Method -------------------------------------------------------------------
elseif comstr(Cam,'method')

   ks.name='spfmex';
   ks.header='SDT sparse LDLt solver';
   ks.SymRenumber='';
   if ~isempty(strfind(Cam,'silent')) %#ok<STREMP> 
    ks.FactBuild='ks=spfmex_utils(''fact'',k,ks,[0 0 size(k,1)/32  0.01*size(k,1)],1);';
   else
    ks.FactBuild='ks=spfmex_utils(''fact'',k,ks,[1 0 size(k,1)/32  0.01*size(k,1)]);';
   end
%   ks.Solve='q=spfmex(''solve'',k.ty(2),full(b));';
   ks.Solve={@doSolve};%'q=spfmex_utils(''solve'',k,full(b));';
   ks.Clear='spfmex_utils(''clear'',ks);';
   ks.Available=exist('spfmex','file')==3;
   ks.HandlesComplex=1;
   ks.TktSolve='ks=feutilb(''tktsolve'',k{:},ind);';

%% #Serv ---------------------------------------------------------------------
elseif comstr(Cam,'serv');[CAM,Cam]=comstr(Cam,5);

%% #ServFact : server factorization - - - - - - - - - - - - - - - - - - - - - 
if comstr(Cam,'fact')

 k=varargin{carg};carg=carg+1;
 ki=varargin{carg};carg=carg+1;

 method=ki.method;
 method.host=getpref('SDT','ServerHost','localhost');
 method.port=getpref('SDT','ServerPort',8888);
 
 [uo,go]=serv_get('uo');
 if ~isfield(uo,'con')||eval('pnet(uo.con,''status'')<=0','1')
  eval('con=pnet(''tcpconnect'',method.host,method.port);');
 else; con=uo.con;
 end

 if con==-1, % If server not running try starting it 
  if isunix 
   st=which(sprintf('sdt_serv_%s.exe',strrep(mexext,'mex','')));
   st=sprintf('unix(''xterm -T SDT_Server -e %s &'')',st);
  else;   st=which('sdt_serv.exe');
         st=sprintf('dos(''%s &'',''-echo'');',st);
  end
  st=getpref('SDT','ServerStart',st);
  eval(st);
  pause(1);  eval('con=pnet(''tcpconnect'',method.host,method.port);');
 end
 if con==-1; error('Bad url or server down.....'); end

 uo.con=con; set(go,'userdata',uo);
 eval('pnet(con,''setwritetimeout'',10);');
 eval('pnet(con,''setreadtimeout'',10);');

 eval('pnet(con,''printf'',''SPOOLESFACT\n'');')

 k=k+spalloc(size(k,1),size(k,1),0);
 eval('pnet(con,''write'',int32([size(k) nnz(k) 0]),''intel'');')
 eval('pnet(con,''sparsewrite'',k,''intel'');')
 fprintf('Matrix transmitted');

 try;
  st='';
  while isempty(st)  
   pause(1); % Try to do some SMART GUESS of how long to wait
   eval('st=pnet(con,''readline'',''noblock'');')
   eval('stat=pnet(con,''status'');');
   if stat<=0; error('The server has diconnected');end
  end
 end
  if ~strcmp(st,'DONE');error('Problem during factorization');end
  i1=double(eval('pnet(con,''Read'',1,''uint32'')')); % factor number

 eval('pnet(con,''printf'',''EMPTY\n'')');
 method.con=con;
 method.Solve='q=spfmex_utils(''servsolve'',k,full(b));';
 method.Clear='spfmex_utils(''servclear'',ks);';
 ki.method=method;
 ks=ki; ks.ty(2:3)=[i1 size(k,1)]; ks.ind='sdt_serv';

%% #ServSolve : server solve - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'solve') % -------------------------------------

 ki=varargin{carg};carg=carg+1;
 b=varargin{carg};carg=carg+1;

 con=ki.method.con;

  eval('pnet(con,''printf'',''SPOOLESSOLVE\n'')'); 
 
  cF=0;
  eval('pnet(con,''write'',int32([size(b) ki.ty(2) 0]),''intel'');')
  eval('pnet(con,''write'',b,''intel'');')

  dataclass='';
  while isempty(dataclass)  
   pause(1); % Try to do some SMART GUESS of how long to wait
   dataclass=eval('pnet(con,''readline'',1024,''noblock'');')
   stat=pnet(con,'status'); if stat<=0; error('The server has diconnected');end
  end

  datadims=double(eval('pnet(con,''Read'',1,''uint32'',''intel'')'));
  datasize=double(eval('pnet(con,''Read'',datadims,''uint32'',''intel'')'));
  fprintf('Receiving %i %i %s\n',datasize,dataclass);
  ks=pnet(con,'Read',datasize,'double','intel');

  if size(ks,1)~=size(b,1); error('Improper return');end
  pnet(con,'printf','EMPTY\n');


%% #ServClear : server Clear - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'clear') % -------------------------------------

  ki=varargin{carg};carg=carg+1;
  try;
   con=ki.method.con;
   pnet(con,'printf','SPOOLESCLEAR\n'); 
   pnet(con,'write',int32([ki.ty(2)]),'intel');
   pnet(con,'printf','EMPTY\n');
  catch
   pnet('closeall');
  end
  
%% #ServOn : server On - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'on')

ks=ofact('spfmex','FactBuild', ...
 'ks=spfmex_utils(''servfact'',k,ks,[1 1 size(k,1)/32  0.01*size(k,1)]);',...
 'Solve','q=spfmex_utils(''servsolve'',k,full(b));', ...
 'Clear','spfmex_utils(''servclear'',ks)');

%% #ServOff : server Off - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'off')

ks=ofact('spfmex','FactBuild', ...
 'ks=spfmex_utils(''fact'',k,ks,[1 0 size(k,1)/32  0.01*size(k,1)]);', ...
 'Solve','q=spfmex_utils(''solve'',k,full(b));', ...
 'Clear','spfmex_utils(''clear'',ks);');

%%
else; sdtw('''Serv%s'' unknown',Cam);
end

%% #Fact : -------------------------------------------------------------------
elseif comstr(Cam,'fact') % -----------------------------------------

 k=varargin{carg};carg=carg+1;
 ki=varargin{carg};carg=carg+1;
 if carg>nargin; opt=[1 0 size(k,1)/32  0.01*size(k,1)];
 else; opt=varargin{carg};carg=carg+1;
 end
 k=check_k(k);t=cputime;
 % opt=[matlabmsglvl msglvl maxdomainsize maxzeros]

 if spfmex==5.0001; i1=spfmex('fact',k); 
 elseif spfmex==5.0003; [i1,ind]=spfmex('fact',k,opt); 
 else; i1=spfmex('fact',k,opt); 
 end
 if opt(1); fprintf(' %.2f s',cputime-t); end
 ks=ki; ks.ty(2:3)=[i1 size(k,1)]; ks.ind='spfmex';
 if spfmex==5.0003; 
     ks.dinv=ind; 
     if isreal(k); ks.method.Solve={@doSolveRealDinv};end
 end

%% #Clear --------------------------------------------------------------------
elseif comstr(Cam,'clear')

 ks=varargin{carg};carg=carg+1;
 if size(ks,1)==0;  spfmex('clear',-1); % DO NOT USE ISEMPTY
 elseif isa(ks,'double') 
    try; spfmex('clear',ks);catch; disp(lasterr);end
 else
    try;  if length(ks.ty)>1&&ks.ty(1)<6;spfmex('clear',ks.ty(2)); end
    catch;disp(lasterr);
    end
 end

%% #Silent -------------------------------------------------------------------
elseif comstr(Cam,'silent') % - - - - - - - - - - - - - - - - - -

 k=varargin{carg};carg=carg+1;
 ki=varargin{carg};carg=carg+1;
 if carg<=nargin; silent=varargin{carg}; carg=carg+1; else; silent=1; end
 ki.FactBuild=sprintf('ks=spfmex_utils(''fact'',k,ks,[%i 0 size(k,1)/32  0.01*size(k,1)]);',silent);
 ks=ki;

%% #Job ----------------------------------------------------------------------
elseif comstr(Cam,'job');[CAM,Cam]=comstr(Cam,4);

%% #JobFact : server factorization - - - - - - - - - - - - - - - - - - - - - - 
if comstr(Cam,'fact')

 k=varargin{carg};carg=carg+1;
 ki=varargin{carg};carg=carg+1;
 [uo,go]=serv_get('job');

 k=check_k(k);t=cputime;
 uo=sdtjob('subStartSlave');  % starts the slave
 [wd,st]=fileparts(uo.job.FileName);
 tname=fullfile(wd,[st '.mat']);save(tname,'k','-append');
 fprintf('Slave Started');
 %sdtjob('subeval',uo,'z=fieldnames(uo);cingui([''rwarn'' sprintf(''\n%s'',z{:})])');
 %sdtjob('subeval',uo,'cingui([''rwarn'' sprintf(''\n%s'',uo.job.FileName)])');
 uo.job.Data(1)=comstr('fact',-32);uo.job.Data(2)=0;
 while uo.job.Data(2)~=2; pause(.1);end
 fprintf('Slave finished factorization');
 
 method.uo=uo;
 method.Solve='q=spfmex_utils(''jobsolve'',k,full(b));';
 method.Clear='spfmex_utils(''jobclear'',ks);';
 ki.method=method;
 ks=ki; ks.ty(2:3)=[uo.job.Data(101) size(k,1)]; 
 ks.ind='SpfmexJob';

%% #JobSolve : server solve - - - - - - - - - - - - - - - - - - - - - - - - - 
elseif comstr(Cam,'solve') % -------------------------------------

 ki=varargin{carg};carg=carg+1;
 b=varargin{carg};carg=carg+1;
 uo=ki.method.uo;
 b=full(b);if numel(b)>numel(uo.job.Data)-300;error('Buflength error');end
 if size(b,1)~=ki.ty(3);error('Row mistmatch');end
 for j1=1:size(b,2)
   if size(b,2)>1;uo.job.Data(300+[0:size(b,1)])=[size(b,1);b(:,j1)];
   else;uo.job.Data(300+[0:numel(b)])=[numel(b);b];
   end
  uo.job.Data(2)=1;%disp(mean(b))
  st=['i1=uo.job.Data(300);i1=300+[1:i1];b=reshape(uo.job.Data(i1),[],1);' ...
     'uo.job.Data(301)=NaN;uo.job.Data(i1)=uo.kd', ...
     num2str(ki.ty(2)) '\b;uo.job.Data(2)=0;'];
  if length(st)>200;error('Command too long');end
  sdtjob('subeval',uo,st);
  if j1==1;ks= uo.job.Data(300+[1:size(b,1)]);
  else; ks(:,j1)=uo.job.Data(300+[1:size(b,1)]);
  end
 end
 if ~isfinite(ks(1));error('Solve failed');end

%% #JobClear - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
elseif comstr(Cam,'clear') % -------------------------------------

  ki=varargin{carg};carg=carg+1;
  r1=ki.method; if ~isfield(r1,'uo'); return;end
  try;
   uo=r1.uo;
   st=['kd=uo.kd',num2str(ki.ty(2)) ';ofact(''clear'',kd)'];
   sdtjob('subeval',uo,st);
  end
  
%% #JobOn - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
elseif comstr(Cam,'on')

ks=ofact('spfmex','FactBuild', ...
 'ks=spfmex_utils(''jobfact'',k,ks,[1 1 size(k,1)/32  0.01*size(k,1)]);',...
 'Solve','q=spfmex_utils(''jobsolve'',k,full(b));', ...
 'Clear','spfmex_utils(''jobclear'',ks)');

%% #JobOff server Off - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
elseif comstr(Cam,'off');ks=spfmex_utils('servoff',varargin{carg});
%%
elseif comstr(Cam,'@');ks=eval(CAM);
else; sdtw('''Job%s'' unknown',Cam);
end

%% #Fact ---------------------------------------------------------------------
elseif comstr(Cam,'fact') % -----------------------------------------

 k=varargin{carg};carg=carg+1;
 ki=varargin{carg};carg=carg+1;
 if carg>nargin; opt=[1 0 size(k,1)/32  0.01*size(k,1)];
 else; opt=varargin{carg};carg=carg+1;
 end
 k=check_k(k);t=cputime;
 % opt=[matlabmsglvl msglvl maxdomainsize maxzeros]

 if spfmex==5.0001; i1=spfmex('fact',k); 
 elseif spfmex==5.0003; [i1,ind]=spfmex('fact',k,opt); 
 else; i1=spfmex('fact',k,opt); 
 end
 if opt(1); fprintf(' %.2f s',cputime-t); end
 ks=ki; ks.ty(2:3)=[i1 size(k,1)]; ks.ind='spfmex';
 if spfmex==5.0003; ks.dinv=ind; end

%% #OfactOptim -------------------------------------------------------------------------
elseif comstr(Cam,'ofactoptim'); [CAM,Cam]=comstr(Cam,11);
    [CAM,Cam,RunOpt.setOpt]=comstr('-setopt',[-25 3],CAM,Cam);out=[];
    if ~RunOpt.setOpt
     K=varargin{carg}; carg=carg+1;
     st1=ofact;
     if comstr(st1,'spfmex'); else; ofact('clear'); ofact('spfmex -silent'); end
     if carg<=nargin; RunOpt=varargin{carg}; carg=carg+1; end
     if carg<=nargin;kd=varargin{carg}; carg=carg+1;else;kd=ofact(1,'lu');end
     [i1,i2,RunOpt.fact]=comstr('fact',[-25 3],CAM,Cam);
     [i1,i2,RunOpt.solve]=comstr('solve',[-25 3],CAM,Cam);
     [CAM,Cam,RunOpt.refine]=comstr('-refine',[-25 3],CAM,Cam);
     [CAM,Cam,RunOpt.plot]=comstr('-plot',[-25 3],CAM,Cam);
     if ~RunOpt.solve&&~RunOpt.fact; RunOpt.solve=1; RunOpt.fact=1; end
     RunOpt=spfmex_utils('ofactoptim-setopt',RunOpt,kd);
     if RunOpt.solve; RunOpt.Nsolve=RunOpt.nCompt; else; RunOpt.Nsolve=1; end
     if RunOpt.fact;  RunOpt.Nfact=RunOpt.nCompt;  else; RunOpt.Nfact=1;end
     if ~isfield(RunOpt,'b')||length(RunOpt.b)~=length(K); RunOpt.b=rand(length(K),1); end
     r1=zeros(length(RunOpt.maxDomain),length(RunOpt.maxZeros),3);
     ofact('silent')
     for j1=1:length(RunOpt.maxDomain)
      for j2=1:length(RunOpt.maxZeros)
       r3=ofact('spfmex','FactBuild',...
        sprintf('ks=spfmex_utils(''fact'',k,ks,[1 0 size(k,1)/%i %g*size(k,1)]);',...
        RunOpt.maxDomain(j1),RunOpt.maxZeros(j2)));
       for j3=1:RunOpt.Nfact % !fact test
        t=cputime;
        k=ofact(K); 
        r1(j1,j2,1)=r1(j1,j2,1)+cputime-t;
       end
       if ~any(size(k)~=length(K))
        t=cputime;
        for j3=1:RunOpt.Nsolve; r2=k\RunOpt.b; end % !solve high level
        r1(j1,j2,2)=cputime-t;
        t=cputime;
        for j3=1:RunOpt.Nsolve; ks=spfmex('solve',k.ty(2),RunOpt.b); end
        r1(j1,j2,3)=cputime-t;  % !solve low level
       end
       ofact('clear',k);
      end
     end
     
    %if nargout>0;
     r5=r1; 
     r5(:,:,1)=r5(:,:,1)/RunOpt.Nfact;   
     r5(:,:,2:3)=r5(:,:,2:3)/RunOpt.Nsolve;   
     out=struct('Y',r1,...
         'X',{{RunOpt.maxDomain(:) RunOpt.maxZeros(:) {'fact';'solve';'solveLow'}}},...
         'Xlab',{{'maxDomains' 'maxZeros' 'computation'}});
    %end
 
    if RunOpt.plot  
     % plot all optim in the same iiplot
     % - Optim as table in iiplot :
     C1=out;
     C1.Y=1e3*r5; % ms %squeeze(r1(:,:,2:3))./RunOpt.Nsolve;
     C1.X{1}=num2cell(C1.X{1});C1.X{2}=num2cell(C1.X{2}); % needed for table display
     C1.PlotInfo={'sub','1 1';
                  'ua.axProp',{'ydir','normal'};'ua.XDimPos',[1,2];
                  'ua.PlotType','image';'ua.YFcn','r3=r3;';
                  'ua.Corner',[1,1.05];
                  'ua.TickInfo',{}};
     %C1=sdsetprop(C1,'PlotInfo','ua.TickInfo',{});
     C1=ii_plp('tickXCell',C1);C1=ii_plp('tickYCell',C1); % Display cell labels
     C1=sdsetprop(C1,'PlotInfo','ColorBar',{'units','normalized', ...
                     'position',[.95 .2 .02 .7],'YAxisLocation','left','FontSize',10, ...
                     '@xlabel',{'String','ms','FontSize',14}});
     iicom('curveinit','Optim',C1)
     % - Differences between high (k\b) and low level (spfmex solve) calls :
     if RunOpt.solve
      figure; 
      plot([reshape(squeeze(r1(:,:,2)),[],1) reshape(squeeze(r1(:,:,3)),[],1)]...
                            /RunOpt.nCompt,'+:');
      legend('High level call','Low level call')
     end
    end

    if any(~[RunOpt.solve RunOpt.fact])
     if RunOpt.fact; r2=squeeze(r1(:,:,1)); end
     if RunOpt.solve;r2=squeeze(r1(:,:,3)); end % Use Low level solve
    else
     r2=r1(:,:,1:2);
     r2(:,:,1)=r2(:,:,1)./max(max(r2(:,:,1)));
     r2(:,:,2)=r2(:,:,2)./max(max(r2(:,:,2)));
     r2=sum(r2,3);  
    end
    [i1,i2]=find(r2==min(min(r2)));
    i1=i1(1); i2=i2(1);
    st=sprintf(['unu1=ofact(''spfmex'',''FactBuild'',''ks=spfmex_utils' ...
     '(''''fact'''',k,ks,[1 0 size(k,1)/%g %.15g*size(k,1)]);'')'],...
        RunOpt.maxDomain(i1),RunOpt.maxZeros(i2));
    sdtweb('_link',st);eval(st);
    
    if RunOpt.refine
      opt=struct('maxDomain',unique(fix(linspace(RunOpt.maxDomain(max(i1-1,1)),...
             RunOpt.maxDomain(min(i1+1,length(RunOpt.maxDomain))),RunOpt.refineStep))),...
                 'maxZeros',linspace(RunOpt.maxZeros(max(i2-1,1)),...
             RunOpt.maxZeros(min(i2+1,length(RunOpt.maxZeros))),RunOpt.refineStep));
      out=spfmex_utils(sprintf('ofactoptim%s',CAM),K,opt,ofact(1,'lu'));
    end
    
     
    else % set default RunOpt
     out=struct('nCompt',100,...
                'maxDomain',2.^[3:7],...
                'maxZeros',logspace(-4,1,6),...
                'refineStep',5);
     if nargin>1&&~isempty(varargin{2})
       Opt=varargin{carg}; carg=carg+1;
       r1=fields(Opt);
       for j1=1:length(r1); out.(r1{j1})=Opt.(r1{j1}); end
     end
    end
    if ~isempty(out);ks=out;end
%% #Renumber
elseif comstr(Cam,'renumber');

 kd=varargin{carg};carg=carg+1;
 i1=kd.dinv; 
 if ~isempty(i1) 
   while carg<=nargin;
       st=varargin{carg};carg=carg+1;
       val=evalin('caller',st);
       if size(val,1)==length(i1)&&size(val,2)==length(i1)
         val=val(i1,i1);assignin('caller',st,val);
       elseif size(val,1)==length(i1)
         val=val(i1,:);assignin('caller',st,val);
       end
   end   
 end
 kd.dinv=1:length(i1);assignin('caller',inputname(2),kd);
  
%% #oProp
elseif strncmpi(Cam,'oprop',5);[CAM,Cam]=comstr(Cam,6);
 out={ ...
 'Sym',{'method','spfmex -silent'}
     };
 if isempty(Cam)&&ischar(varargin{carg});% a=ofact('oprop','spfmex','Sym')
     [CAM,Cam]=comstr(varargin{carg},1);carg=carg+1;
 end
 if ~isempty(Cam)
     out=out{strcmpi(out(:,1),CAM),2};
 end
 ks=out;   
else; error('%s unknown',Cam);
end

% #Serv_Get opens a server container figure
function [out,out1]=serv_get(Cam);

if comstr(Cam,'uo')
  gf=sdtdef('rwarn');
  go=findall(gf,'tag','sdt_serv');
  if isempty(go)
   go=uicontrol('Style','list','String','','tag','sdt_serv', ...
        'value',1,'parent',gf,  ...
        'callback','','visible','off','HandleVisibility','off');
  end
  uo=get(go,'userdata');
  out=uo; out1=go; 
elseif comstr(Cam,'job')
  gf=sdtdef('rwarn');go=findall(gf,'tag','JobContainer');
  if isempty(go);uo=[];else; uo=get(go,'userdata');end
  out=uo; out1=go; 
else; error('Not valid');
end

function k=check_k(k);

if ~issparse(k); k = sparse(k); end
i1=find(diag(k)==0);
if ~isempty(i1) 
     if length(i1)<10
      error('Zero terms on diagonal :%s\n',sprintf('%i ',i1));  %#ok<SPERR>
     else; error('%i zero terms on diagonal \n',length(i1)); 
     end   
     return
end

%% #doSolve
function ks=doSolveRealDinv(k,b);
 i1=k.ty;
 if size(b,1)~=i1(3) ; 
   m=k.method;
   if isfield(m,'TktSolve')&&size(m.TktSolve,1)==size(b,1)
    b=m.TktSolve'*b;   
    ks=zeros(size(b,1),size(b,2));
    b=b(k.dinv,:); 
    spfmex('solve',k.ty(2),b,ks,int32(k.dinv-1));
    ks=m.TktSolve*ks;    
   else; error('RHS has a bad size'); 
   end
 elseif isreal(b) % real numbering here
  ks=zeros(size(b,1),size(b,2));
  b=b(k.dinv,:); 
  spfmex('solve',k.ty(2),b,ks,int32(k.dinv-1));
  
 else;ks=doSolve(k,b);
 end

function ks=doSolve(k,b);
 i1=k.ty;
 if size(b,1)~=i1(3) ; error('RHS has a bad size');  
 elseif ~isreal(b) || ~spfmex('isreal',k.ty(2))
  if isempty(k.dinv) % complex numbering outside
    ks=spfmex('solve',k.ty(2),complex(full(b)));
  else % complex numbering here
    if ~isreal(b)
      ks(k.dinv,:)=spfmex('solve',k.ty(2),full(b(k.dinv,:)));
    else
      ks(k.dinv,:)=spfmex('solve',k.ty(2),complex(full(b(k.dinv,:))));
    end
  end
 elseif isempty(k.dinv)  % real numbered outside
  ks=spfmex('solve',k.ty(2),b);
 else % real numbering here
  if size(k.dinv,1)~=i1(3); error('Bad size of RHS'); end
  ks=zeros(size(b,1),size(b,2));
  %ks(k.dinv,:)=spfmex('solve',k.ty(2),full(b(k.dinv,:)));
  b=b(k.dinv,:); 
  spfmex('solve',k.ty(2),b,ks,int32(k.dinv-1));
 end


