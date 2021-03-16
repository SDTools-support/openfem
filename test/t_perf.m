

% Standardized performance test of OpenFEM and ofact sparse solvers
%
% Test larger models with
% RO.Gen={'div','10 20 100'};t_perf
% RO.Gen={'div','10 40 100'};t_perf
% RO=struct('Gen',{{'div','10 40 100'}},'methods',{{'spfmex','mklserv_utils'}});t_perf
% RO=struct('Gen',{{'div','10 40 100'}},'methods',{{'mklserv_utils -silent'}});t_perf

if ~isempty(gcbf)&&strcmp(get(gcbf,'tag'),'test_sdt');return;end

if ~exist('RO','var');RO=struct;end
if ~isfield(RO,'Gen'); RO.Gen={'div','10 10 100'};end
if ~isfield(RO,'methods')
  % select methods to test
 RO.methods=setdiff(ofact('methodlist'),{'sp_util'});
end

out=RO.Gen(1,:);femesh('reset');
femesh(sprintf('testhexa8 divide %s',out{strcmp(RO.Gen(:,1),'div'),2}));
model=femesh('model0');
% [model.Node,model.Elt]=feutil('lin2quad',model);
if 1==2
 %% OpenMP assembly performance
 model=femesh('testhexa8 divide 10 10 100 -rand.01');
 model=feutil('lin2quad',model);
 
 [C1,model.DOF]=fe_mknl('init',model);
 tic;
  maxNumCompThreads(8);of_mk('setomppro',8);
  k=fe_mknl('assemble NoT',model,C1,1);
 toc
 tic;
  maxNumCompThreads(1);of_mk('setomppro',1);
  k1=fe_mknl('assemble NoT',model,C1,1);
 toc
 if norm(k1-k,'inf')/norm(k,'inf')>1e-10;error('mismatch'); end
 
end
t=cputime;[m,k,mdof]=fe_mknl(model);
 out(end+1,1:2)={'Assembly',sprintf('%.2f',cputime-t)};
if size(k,1)>50e3;RO.methods=setdiff(RO.methods,{'lu','chol'});end


if strcmp(getenv('USERDOMAIN'),'epicea')
 RO.methods=setdiff(RO.methods,{'pardiso'});
end
% 
k=k+1e3*m;out(end+1,1:4)={'Method','Factor','10*FBS 1','FBS10'};


for j1=1:length(RO.methods);
 try;  %#ok<*NOSEM>
  ofact(['method ' RO.methods{j1}]); b=ofact(speye(10),ones(10,1));
  b=ones(size(k,1),1); 
  t=cputime; tic;kd=ofact(k); 
  out(end+1,1:2)={RO.methods{j1},sprintf('%.2f/%.2f',cputime-t,toc)};
  t=cputime;tic; 
  for j2=1:10; q=kd\b;end;out(end,3)={sprintf('%.2f/%.2f',cputime-t,toc)};
  b=repmat(b,1,10);t=cputime;tic; 
  q=kd\b;out(end,4)={sprintf('%.2f/%.2f',cputime-t,toc)};
  ofact('clear',kd);
 end
 disp(out);
end

% generation some information about the machine

try; st=getenv('ARCH'); if isempty(st); st=getenv('MATLAB_ARCH');end;end
out(end+1,1:3)={'ARCH',st,version};
st='';try;st=system('uname -a');end
if isempty(st);st=getenv('HOSTNAME');end
if length(st)>1; out(end+1,1:2)={'uname',st};end
st=getenv('PROCESSOR_IDENTIFIER');if length(st)>1; out(end+1,1:2)={'Processor',st};end

% display the result in a suitable format

st=sprintf('Result={\n'); i1=1;
for j1=1:size(out,1)
 if strcmp(out{j1,1},'ARCH');i1=0;end
 if i1;st=[st sprintf('''%-8s'' ',out{j1,:}) sprintf('\n')];
 else;st=[st sprintf('''%s'' ',out{j1,:}) sprintf('\n')];
 end
end
st=[st sprintf('};\n')];
disp(st);

return
% Below are save results on various machines ---------------------------

Result={ % Platane (bi-opteron) 64 bit windows
'div     ' '10 10 100' '        ' '        ' 
'Assembly' '7.23    ' '        ' '        ' 
'Method  ' 'Factor  ' '10*FBS 1' 'FBS10   ' 
'chol    ' '6.81    ' '7.56    ' '5.56    ' 
'lu      ' '52.59   ' '7.44    ' '4.28    ' 
'spfmex  ' '28.53   ' '2.63    ' '1.31    ' 
'ARCH' 'win64' '7.4.0.287 (R2007a)' '' 
'Processor' 'AMD64 Family 15 Model 37 Stepping 1, AuthenticAMD' '' '' 
};

Result={% Platane (bi-opteron) but 32 bit
'div     ' '10 10 100' '        ' '        ' 
'Assembly' '7.84    ' '        ' '        ' 
'Method  ' 'Factor  ' '10*FBS 1' 'FBS10   ' 
'chol    ' '6.83    ' '5.16    ' '3.70    ' 
'lu      ' '45.16   ' '5.55    ' '3.75    ' 
'mtaucs  ' '7.67    ' '1.92    ' '1.91    ' 
'spfmex  ' '13.33   ' '1.38    ' '0.58    ' 
'ARCH' 'win32' '7.3.0.267 (R2006b)' '' 
'Processor' 'AMD64 Family 15 Model 37 Stepping 1, AuthenticAMD' '' '' 
};


