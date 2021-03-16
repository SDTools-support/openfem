function [out,out1,out2]=basic_elt_test(varargin);

% ELEMENT TESTS
%
% eig, load, mat for all elements => comparisons with SDT5.1
% basic_elt_test('compare')
% basic_elt_test mat+ pyra5 % create ref values for new test
%
% test of basic commands (integinfo, call, etc.) for all elements
% basic_elt_test('calls')
% basic_elt_test('all')
% basic_elt_test('eig')
% basic_elt_test('mat')
% basic_elt_test('load')  basic_elt_test('load',ElemF)
% basic_elt_test('dilat') % checks stresses for uniform dilatation
% out = basic_elt_test('load energy')
% 
% of_mk('setomppro',2);for j1=1:10;basic_elt_test('load;','hexa20b');end

%#ok<*NOSEM,*NASGU,*ASGLU>


%	Etienne Balmes, J.M. Leclere
%       Copyright (c) 2001-2017 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       $Revision: 1.71 $  $Date: 2019/05/14 17:26:51 $

femesh('TestUseLegacy',1); % switch femesh test to legacy integration rule
warning('off','OpenFEM:legacy');

carg=1;
if carg<=nargin; [CAM,Cam]=comstr(varargin{carg},1);carg=carg+1;
else
 if ~isempty(gcbf)&&strcmp(get(gcbf,'tag'),'test_sdt');
  wd=pwd;
  cd(fileparts(fileparts(which('fe_mk'))));eval('sdtcheck path');
   clear global;comgui('close all');cinguj('initSwing');
  cd(wd);evalin('base','basic_elt_test(''all'');'); return;
 else;Cam='all';
 end
end

%delete basic_elt_test.log

% List of available tests - - - - - - - - - - - - - - - - - - - - - - - 
 list={ ...
     'q4p','eig_0','load_0','mat_0','visc'
     'q8p','eig_0','load_0','mat_0','visc'
     't3p','eig_0','load_0','mat_0','visc'
     't6p','eig_0','load_0','mat_0','visc'
     'q4p','eig_1','load_1','mat_1','visc'
     'q8p','eig_1','load_1','mat_1','visc'
     't3p','eig_1','load_1','mat_1','visc' 
     't6p','eig_1','load_1','mat_1','visc'
     'q4p','eig_2','load_2','mat_2','visc'
     'q8p','eig_2','load_2','mat_2','visc'
     't3p','eig_2','load_2','mat_2','visc' 
     't6p','eig_2','load_2','mat_2','visc'
     'hexa8','eig','load','mat','visc'
     'hexa20','eig','load','mat','visc'
     'penta6','eig','load','mat','visc'
     'penta15','eig','load','mat','visc'
     'tetra4','eig','load','mat','visc'
     'tetra10','eig','load','mat','visc'
     'q4pb','eig_0','load_0','mat_0','visc'
     'q8pb','eig_0','load_0','mat_0','visc'
     't3pb','eig_0','load_0','mat_0','visc'
     't6pb','eig_0','load_0','mat_0','visc'
     'q4pb','eig_1','load_1','mat_1','visc'
     'q8pb','eig_1','load_1','mat_1','visc'
     't3pb','eig_1','load_1','mat_1','visc'
     't6pb','eig_1','load_1','mat_1','visc'
     'pyra5','eig','load','mat','visc'
     'pyra13','eig','load','mat','visc'
     'hexa8b','eig','load','mat','visc'
     'hexa20b','eig','load','mat','visc'
     'penta6b','eig','load','mat','visc'
     'penta15b','eig','load','mat','visc'
     'tetra4b','eig','load','mat','visc'
     'tetra10b','eig','load','mat','visc'
     'tria3','eig','load','mat','visc'
     'tria6','','','','visc'
     'quad4','eig','load','mat','visc'
%
     'quad4','eig_1','load_1','mat_1','visc_1'
     'quad4','eig_2','load_2','mat_2','visc_2'
     'quad4','eig_5','load_5','mat_5','visc_5'
%
     'quadb','eig','','mat','visc'
     'quad9','','','mat','visc'
     'bar1','eig','load','mat','visc'
     'fsc3','','load','','visc'
     'flui4','eig','','mat','visc'
     'flui6','eig','','mat','visc'
     'flui8','eig','','mat','visc'
     'beam1','eig','load','mat','visc'
     'beam1t','eig','load','mat','visc'
     };
% List of available tests - - - - - - - - - - - - - - - - - - - - - - - 
if comstr(Cam,'list'); out=list; return; end
%-----------------------------------------------------------------------
% testmat for all elements

if ~isempty(strfind(Cam,'enerk')); flag='enerk';
else; flag='';
end
if ~isempty(strfind(Cam,'stress')); flag='stress';
else; flag='';
end
if ~isempty(strfind(Cam,'visc')); flag='visc';
else; flag='';
end
RO.Silent=1; 
if Cam(end)=='+';RO.Silent=0;CAM(end)='';Cam(end)='';
    dbstop if caught error
end

if comstr(Cam,'mat');     jTest=4;
elseif comstr(Cam,'eig'); jTest=2;
elseif comstr(Cam,'visc'); jTest=5;
%xxxelseif comstr(Cam,'load energy'); jTest=5;
elseif comstr(Cam,'load'); jTest=3;
else; jTest=[];
end
if carg<=nargin; st=varargin{carg};carg=carg+1; % do a partial test if argument
 i1=strmatch(st,list(:,1)); if ~isempty(i1); list=list(i1,:);end
end

ofact silent
if ~isempty(jTest) % basic series of tests (mat eig load)

out={Cam,'data','CheckValue','Ref'};
try; ref_elt_test;end % Loads references 

for j1=1:size(list,1);
 if ~isempty(list{j1,jTest})
   if comstr(Cam,'mat') % - - - - - - - - - - - - - - - - - - - - - -
     k=[];
     st=sprintf('k=%s(''test%s'');',list{j1,[1 4]}); 
     eval_disp(st,Cam);
     if iscell(k)
      r2=(1:size(k{1},1))';
      r3=[max(svd(full(k{1}))) max(svd(full(k{2}))) diag(r2'*k{1}*r2) diag(r2'*k{2}*r2)];
      %comstr(r3,-30)
      out(end+1,1:3)={st,k,[max(svd(full(k{1}))) max(svd(full(k{2})))]};
     else;out(end+1,1:3)={st,[],'Failed'};
     end
     if ~RO.Silent; fprintf('\n[%s]\n',sprintf('%16.16g ',out{end,3}));end
     try; i4=strmatch(out{end,1},Mat_ref(:,1));
          r1=Mat_ref{i4,2}; 
          if length(r1)>2; r2=(1:size(k{1},1))';
              out{end,3}=[out{end,3} feutilb('dtkt',r2,k)];
          end
          if ischar(r1); comstr(out{end,3},-30)
              out{end,4}='no_ref';
          elseif norm((r1-out{end,3})/norm(r1))>sdtdef('epsl'); out{end,4}=r1;
          end
     end
     if ~isempty(stu); out{end,4}=stu;end
   elseif comstr(Cam,'eig') % - - - - - - - - - - - - - - - - - - - - - -
     k=[];
     st=sprintf('[model,def]=%s(''test%s%s'');',list{j1,[1 2]},flag); 
     eval_disp(st,Cam);[eltid,model.Elt]=feutil('eltidfix;',model.Elt);
     r1=fe_stress('enerk',model,def);
     r1=fe_stress('stress',model,def);
     out(end+1,1:3)={st,def,def.data(1:5)'};
     if ~RO.Silent; fprintf('\n[%s]\n',sprintf('%16.16g ',def.data(1:5)'));end
     try; i4=strmatch(out{end,1},Eig_ref(:,1));
          r1=Eig_ref{i4,2}(:)'; 
          if ischar(r1);sdtw('No ref %s',Eig_ref{i4,1});
          elseif norm(r1-out{end,3})/norm(r1)>1e-6; out{end,4}=r1;
          end
     end
     if ~isempty(stu); out{end,4}=stu;end
   elseif comstr(Cam,'load') % - - - - - - - - - - - - - - - - - - - - - -
     st=sprintf('[model,def]=%s(''test%s%s'');',list{j1,[1 3]},flag);
     eval_disp(st,Cam);
     out(end+1,1:3)={st,def,sqrt(sum(full(def.def).^2,1))}; %#ok<AGROW,NODEF>
     if ~RO.Silent; fprintf('\n[%s]\n',sprintf('%16.16g ',sqrt(sum(full(def.def).^2,1))));end
     try; i4=strmatch(out{end,1},Load_ref(:,1)); %#ok<USENS>
          if isempty(i4);
          end 
          r1=Load_ref{i4,2};
          if ~ischar(r1) && ~isempty(find(~isfinite(r1)))
           r2=sqrt(sum(full(def.def).^2,1));
           r2=r2(1);r1=r1(1);
           out{end,3}=r2;
           if norm(r1-r2)/norm(r1)>sdtdef('epsl');out{end,4}=r1;end
          elseif ~ischar(r1)
           r2=sqrt(sum(full(def.def).^2,1));%comstr(r2,-30)
           if length(r1)==1; r2=norm(r2); end
           if norm(r1-r2)/norm(r1)>sdtdef('epsl');out{end,4}=r1;end
          else; sdtw('No ref %s',Load_ref{i4,1});
          end
     end

     if ~isempty(stu); out{end,4}=stu;end
   elseif comstr(Cam,'visc') 
     %% - - - - - - - - - - - - - - - - - - - - - -
     st=sprintf('model=(''teststruct%s -back'');',list{j1,1});
     eval_disp(st,Cam);
     [mat,mdof] = fe_mk(model,'options',3);
     [Case,model.DOF]=fe_mknl('init',model);
     mat1=fe_mknl('assemble',model,Case,3);
     out(end+1,1:3)={st,{mat,mat},{size(find(mat)),size(find(mat1))}};
     if ~isempty(stu); out{end,4}=stu;end
   end

 end
end

summary={};
for j1=2:size(out,1)
  if isempty(out{j1,4})
  elseif ischar(out{j1,4});
    summary{end+1}=sprintf(' %-30s %s',out{j1,1},out{j1,4}); %#ok<AGROW>
  else; summary{end+1}=sprintf(' %-30s DIFFER : [%s]\n current [%s]',out{j1,1},...
                       sprintf('%16.16g ',[out{j1,3}-out{j1,4}]), ...
                       sprintf('%16.16g ',[out{j1,3}])); %#ok<AGROW>
  end
end % j1
out1=summary;
 
if nargout==0; 

 if isempty(out1); 
   if ~isempty(Cam)&&Cam(end)==';';%fprintf('EVERYTHING WENT WELL');
   else
    fprintf('\n\nEVERYTHING WENT WELL\n\n');
   end
 else; fprintf('\n\nTHE FOLLOWING FAILED\n\n');fprintf('%s\n',out1{:});
 end
 fprintf('\n');
 if length(out)>10
   st=fullfile(sdtdef('tempdir'),'elt_test.mat');save(st,'out')
   fprintf('Result saved in : %s\n',st);
 end
 clear out;
end

%-----------------------------------------------------------------------
% Verification of dilatation for volume elements
elseif comstr(Cam,'dilat')


 % 3D  : Uniform dilatation check

 list={'hexa8b rand.1','hexa8 rand.1','hexa20 rand.1', ...
       'tetra4b rand.1'}';list=list(:);
 for j1=1:length(list)
  femesh('reset');
  model=femesh(horzcat(';test',list{j1,1},'struct back'));
  model.pl=m_elastic('dbval 100 steel','dbval 112 aluminum');
  model.pl(:,4)=.285;model.pl(2,6)=model.pl(2,3)/2/(1.285);
  nind=sparse(model.pl(:,1),1,1:size(model.pl)); mpid=feutil('mpid',model);
  def=feutil('geomdilat',model,model.pl(nind(mpid(2,1)),4));
  s1=fe_stress('stress sxx atcenter',model,def);
  r1=norm([s1.data./model.pl(nind(mpid(s1.EltId+1,1)),3)]-1);
  s1=fe_stress('stress syy atcenter',model,def);  r1(2)=norm(s1.data);
  list(j1,2:3)={r1(1) r1(2)};
 end
 out=list;

 % 2D  : Uniform dilatation check

 list={'q4pb rand.3','q4p rand.3'};list=list(:);
 for j1=1:length(list)
  femesh('reset');
  model=femesh(horzcat(';test',list{j1,1},'divide 4 4 struct back'));
  model.pl=m_elastic('dbval 100 steel','dbval 112 aluminum');
  model.pl(:,4)=.285;model.pl(2,6)=model.pl(2,3)/2/(1.285);
  nind=sparse(model.pl(:,1),1,1:size(model.pl)); mpid=feutil('mpid',model);
  def=feutil('geomdilat',model,model.pl(nind(mpid(2,1)),4));
  s1=fe_stress('stress sxx atcenter',model,def);
  r1=norm([s1.data./model.pl(nind(mpid(s1.EltId+1,1)),3)]-1);
  s1=fe_stress('stress syy atcenter',model,def);  r1(2)=norm(s1.data);
  list(j1,2:3)={r1(1) r1(2)};
 end
 out(end+[1:size(list,1)],1:size(list,2))=list;


%-----------------------------------------------------------------------
% 
elseif comstr(Cam,'compare')

 st0=list; ref_elt_test;
 [st0,i1]=intersect(st,st0);
 if nargin==1; 
  sdtw('You should use basic_elt_test_all to arrive here');return;
 end  
 out=varargin{carg}; carg=carg+1;

 % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - Matrices
 summary={}; Mat=local{2}; result_Mat=Mat;
 for j1=2:size(Mat,1)
  if ~isempty(Mat{j1,4})
   summary{end+1}=sprintf(' %30s DIFFER\n',Mat{j1,1}); %#ok<AGROW>
  end
 end % j1
 % this is an other check for properties
 a=femesh('testbeam1');
 if norm(a.pl-pl)>sdtdef('epsl');  error('Pb on PL properties'); end
 if norm(a.il-il)>sdtdef('epsl');  error('Pb on IL properties'); end

 % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - EIG
 Eig=out{3};

 result_Eig=Eig;
 for j1=1:size(Eig,1)
  i1=find(strcmp(Eig_ref,Eig{j1,1}));
  if ~isempty(i1) 
   try
    if  ~ischar(Eig_ref{i1,2}) && ~ischar(Eig{j1,2}) 
     result_Eig{j1,2}=(Eig_ref{i1,2}(:)'-Eig{j1,3}(:)')./Eig_ref{i1,2}(:)';
    else
     result_Eig{j1,2}='Error';
    end
   catch
     result_Eig{j1,2}='Error';
   end
  else
    result_Eig{j1,2}='not compared';
  end
 end % j1
 % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - LOAD
 Load = basic_elt_test('load');result_Load=Load;
 for j1=1:size(Load,1)
  i1=find(strcmp(Load_ref,Load{j1,1}));
  if ~isempty(i1) 
   try
    if  ~ischar(Load_ref{i1,2}) && ~ischar(Load{j1,2})
     i2=1; % look if column are not zeros
     for j2=1:size(Load{j1,2}.def,2)
      if norm(full(Load{j1,2}.def(:,j2)))==0; i2=0; break; end
     end
     if i2==0
      result_Load{j1,2}=sprintf('One column is zero, diff=%16.16g ',(Load_ref{i1,2}-Load{j1,3})./Load_ref{i1,2});
     else
      result_Load{j1,2}=(Load_ref{i1,2}-Load{j1,3})./Load_ref{i1,2};
     end
    else
     result_Load{j1,2}='Error';
    end
   catch
    result_Load{j1,2}='Error';
   end
  else
    result_Load{j1,2}='not compared';
  end
 end % j1
 % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - STRESS
if 1==2 %xxx
 Stress = basic_elt_test('stress');result_Stress=Stress;
 for j1=1:size(Stress,1)
  i1=find(strcmp(Stress_ref,Stress{j1,1}));
  if ~isempty(i1) 
   if  ~ischar(Stress_ref{j1,2}) && ~ischar(Stress{j1,2})
    result_Stress{j1,2}=Stress_ref{i1,2}-Stress{i1,3};
   else
    result_Stress{j1,2}='Error';
   end
  else
    result_Stress{j1,2}='not compared';
  end
 end % j1
end
 
 % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

 out={result_Mat,result_Eig,result_Load};

 % - - - - - - - - - - - - - - - - - - - - - - - - - - -  Comparisons output
 for j1=1:size(out{1},1) % Mat
  if ~ischar(out{1}{j1,2}) && ~isempty(find(out{1}{j1,2}>sdtdef('epsl')))
    fprintf('%s : Error ',out{1}{j1,1});
    fprintf('%f ',out{1}{j1,2});
    fprintf('\n');
  elseif ~ischar(out{1}{j1,2}) && isempty(find(out{1}{j1,2}>sdtdef('epsl'))) %#ok<EFIND>
    %fprintf('%s : comparison OK\n',out{1}{j1,1});
  elseif ischar(out{1}{j1,2}) && comstr(out{1}{j1,2},'not compared')
    fprintf('%s : not compared\n',out{1}{j1,1});
  else fprintf('%s : comparison Error\n',out{1}{j1,1});
  end
 end

 for j1=1:size(out{2},1) % Eig
  if ~ischar(out{2}{j1,2}) && ~isempty(find(out{2}{j1,2}>sdtdef('epsl'))) %#ok<EFIND>
    fprintf('%s : Error ',out{2}{j1,1});
    fprintf('%f ',out{2}{j1,2});
    fprintf('\n');
  elseif ~ischar(out{2}{j1,2}) && isempty(find(out{2}{j1,2}>sdtdef('epsl'))) %#ok<EFIND>
    %fprintf('%s : comparison OK\n',out{2}{j1,1});
  elseif ischar(out{2}{j1,2}) && comstr(out{2}{j1,2},'not compared')
    fprintf('%s : not compared\n',out{2}{j1,1});
  else fprintf('%s : comparison Error\n',out{2}{j1,1});
  end
 end

 for j1=1:size(out{3},1) % Load
  if ~ischar(out{3}{j1,2}) && ~isempty(find(out{3}{j1,2}>sdtdef('epsl'))) %#ok<EFIND>
    fprintf('%s : Error ',out{3}{j1,1});
    fprintf('%f ',out{3}{j1,2});
    fprintf('\n');
  elseif ~ischar(out{3}{j1,2}) && isempty(find(out{3}{j1,2}>sdtdef('epsl'))) %#ok<EFIND>
    %fprintf('%s : comparison OK\n',out{3}{j1,1});
  elseif ischar(out{3}{j1,2}) && comstr(out{3}{j1,2},'not compared')
    fprintf('%s : not compared\n',out{3}{j1,1});
  elseif ischar(out{3}{j1,2}) && comstr(out{3}{j1,2},'One column')
    fprintf('%s : %s \n',out{3}{j1,1},out{3}{j1,2});
  else fprintf('%s : comparison Error\n',out{3}{j1,1});
  end
 end
 return;

%-----------------------------------------------------------------------
% test each basic command of elements
elseif comstr(Cam,'call')

 st1={'call','rhscall','node','patch',...
      'dof','line','face','parent'};

 for j1=1:length(list)
  for j2=1:length(st1)
    st3=['r1=' list{j1,1} '(''' st1{j2} ''');'];
    try ;
     eval_disp(st3)
     if isempty(r1); r1=sprintf('Empty %s',st1{j2});else; r1=''; end
    catch
     r1=lasterr;
    end
  list{j1,j2+1}=r1;
  end %j2
 end %j1

 
 j2=2;while j2<=size(list,2) % empty columns
  r1=char(list(:,j2)); 
  if isempty(r1) list(:,j2)=[];st1(j2-1)=[];
  else; j2=j2+1;
  end
 end
 j2=1;while j2<=size(list,1) % empty rows
  if ~isempty(list{j2,2:end})
   r1=strcat(list{j2,2:end}); 
  else r1='';
  end
  if isempty(r1) list(j2,:)=[];
  else; j2=j2+1;
  end
 end

 st1={'' st1{:}};st1(end+[1:size(list,1)],:)=list;
 out=st1;
 out1=list(:,1);
 for j1=1:1:size(list,1);out1{j1}=sprintf('%-30s:%s',list{j1,:});end
 

%-----------------------------------------------------------------------
% test each basic command of elements
elseif comstr(Cam,'all')

if strcmp(get(gcbf,'tag'),'test_sdt');set(0,'recursionlimit',50);end
st='_';st=st(ones(1,60));
out=cell(2,5);out(3,1:4)={'call','mat','eig','load'};

 fprintf('\n%s\n\n Checking calls \n%s\n',st,st);
 [out{1,1},out{2,1}]=basic_elt_test('call');

 fprintf('\n%s\n\n Checking matrix tests \n%s\n',st,st);
 [out{1,2},out{2,2}]=basic_elt_test('mat');


 fprintf('\n%s\n\n Checking eig tests \n%s\n',st,st);
 [out{1,3},out{2,3}]=basic_elt_test('eig');

 fprintf('\n%s\n\n Checking load tests \n%s\n',st,st);
 [out{1,4},out{2,4}]=basic_elt_test('load');

 %basic_elt_test('compare');
 for j1=1:4
  out1=out{2,j1};
  if isempty(out1); fprintf('\n\n%s EVERYTHING WENT WELL\n\n',out{3,j1});
  elseif iscell(out1) 
   fprintf('\n\n%s THE FOLLOWING FAILED\n\n',out{3,j1});
   fprintf('%s\n',out1{:});
  else
   fprintf('\n\n%s THE FOLLOWING FAILED\n\n',out{3,j1});
   fprintf('%s\n',out1);
  end
 end

%clear out

elseif comstr(Cam,'OLD')

%ref_elt_test
%cd d:\dis_sdt\sdt5.1;sdtcheck path;load i:/tmp/xxx_ref.mat

for j1=1:size(Load_ref(:,1),1)
 try ;
  eval(Load_ref{j1,1}); 
  fprintf('%16.16f       %16.16f       %16.16f  \n',sqrt(sum(full(def.def).^2,1)));
 catch
  fprintf('%s : Problem\n',Load_ref(:,1))
 end
end


else; disp(sprintf('Error : Not a valid test => %s ',varargin{1}))
end


return;


%-----------------------------------------------------------------------
% Actually do the test for one specified element

st2=sprintf('%s test%s;',st,st1);

try ; eval_disp(st2)
catch
    disp(lasterr);
    fprintf(' ------------------ Error in %s test\n%s\n', st2)
end

if any(strcmp(st,{'q4p','q5p','q8p','t3p','t6p','q9a'}))
 st2=sprintf('%s test%s_0;',st,st1);
 try ; eval_disp(st2)
 catch
    disp(lasterr);
    fprintf(' ------------------ Error in %s test\n%s\n',st2)
 end
 st2=sprintf('%s test%s_1;',st,st1);
 try ; eval_disp(st2)
 catch
    disp(lasterr);
    fprintf(' ------------------ Error in %s test\n%s\n',st2)
 end
 st2=sprintf('%s test%s_2;',st,st1);
 try ; eval_disp(st2)
 catch
    disp(lasterr);
    fprintf(' ------------------ Error in %s test\n%s\n',st2)
 end

 % test anisotropic material
 % plane stress
 E=2.1e11;nu=.285;C=E/(1.-nu*nu);
 e=[C C*nu C 0. 0. C*(1-nu)/2];
 pl=[100 fe_mat('m_elastic','SI',4) e 7800 0 0 0 0 .1];
 st2=sprintf('k=%s (''test%s_0'');',st,st1)
 eval_disp(st2);
 st2=sprintf('k1=%s (''test%s_0'',pl);',st,st1)
 eval_disp(st2);
 if isa(k1,'cell')
   if norm(k{1}-k1{1}) error(['anisotropic element ' st2]); end
 end

 % plane strain
 unmnu=1.-nu;
 C=E*unmnu/(1+nu)/(1-2*nu);
 e=[C nu*C/unmnu C 0. 0. C*(1-2*nu)/(2*unmnu)];
 pl=[100 fe_mat('m_elastic','SI',4) e 7800 0 0 0 0 .1];
 st2=sprintf('k=%s (''test%s_1'');',st,st1);
 eval_disp(st2);
 st2=sprintf('k1=%s (''test%s_1'',pl);',st,st1);
 eval_disp(st2);
 if isa(k1,'cell')
   if norm(k{1}-k1{1}) error(['anisotropic element ' st2]); end
 end


end

return;




%-----------------------------------------------------------------------
if 1==2
 % anisotropic
  femesh('reset');
  model=femesh(';testq4p');
  model.il(3)=0;

  E=2.1e11;nu=.285;C=E/(1.-nu*nu);
  e=[C C*nu C 0. 0. C*(1-nu)/2];
  model.pl=[100 fe_mat('m_elastic','SI',4) e 7800 0 0 0 0 .1];

  [constit,iopt,elmap]=q4p('integinfo',[100;110],model.pl,model.il);
  [k,m]=q4p(model.Node,model.Elt(2,:),[36 36 0 0 0 0 0 0 0],int32(iopt),constit,elmap);


     unmnu=1.-nu;
     C=E*unmnu/(1+nu)/(1-2*nu);
     e=[C nu*C/unmnu C 0. 0. C*(1-2*nu)/(2*unmnu)];
pl=[100 fe_mat('m_elastic','SI',4) e 7800 0 0 0 0 .1];
end





%--------------------------------------------------------- RHS
if 1==2
% st={'q4p','q5p','q8p','t3p','t6p',...
%      'hexa8','hexa20','penta6',...
%     'tetra4','tetra10'};

% st={'quad4','tria3','quadb'};

%st={'bar1','beam1'};
%st={'flui4','flui6','flui8'};

 st={'q4p','q5p','q8p','t3p','t6p'};
 
errors={}; 

for j2=1:4

 for j1=1:length(st)

    st1=['[model,def1]=femesh('';teststruct' st{j1} 'load'');' ];
    st2=['[model,def2]=femesh('';teststruct' st{j1} 'load'');' ];
    st3=['[model,def3]=femesh('';teststruct' st{j1} 'load'');' ];

    femesh('reset');eval_disp(st1);
    femesh('reset');eval_disp(st2);
    femesh('reset');eval_disp(st3);

    if ~isequal(def1.def,def2.def) | ~isequal(def1.def,def3.def)
     errors{end+1,1}=st1;
     if ~isequal(def1.def(:,1:2),def2.def(:,1:2)) |...
                            ~isequal(def1.def(:,1:2),def3.def(:,1:2))
       errors{end,2}='gravity and/or surf';
     else
       errors{end,2}='pressure only';
     end
    end 
 end %j1
end %j2



end

%------------------------------------------ compare basic matrices
if 1==2
out=basic_elt_test('mat')
load basic_elt_test_matrices
result=out;

for j1=1:size(out,1)
  i1=find(strcmp(k,out{j1,1}));

  if ~isempty(i1) 
   if  ~ischar(out{j1,2});
    result{j1,2}=norm(out{j1,2}{1}-k{i1,2}{1});
    result{j1,3}=norm(out{j1,2}{2}-k{i1,2}{2});
   else
    result{j1,2}='Error';
    result{j1,3}='Error';
   end
  else
    result{j1,2}='not compared';
  end
end

end
%--------------------------------------------------------- RHS
if 1==2
out=basic_elt_test('load')
load basic_elt_test_load
result=out;

for j1=1:size(out,1)
  i1=find(strcmp(k,out{j1,1}));

  if ~isempty(i1) 
   if  ~ischar(out{j1,2});
    result{j1,2}=norm(out{j1,3}-k{i1,3});
   else
    result{j1,2}='Error';
   end
  else
    result{j1,2}='not compared';
  end
end

end
%--------------------------------------------------------- EIG
if 1==2
out=basic_elt_test('eig')
load basic_elt_test_load
result=out;

for j1=1:size(out,1)
  i1=find(strcmp(k,out{j1,1}));

  if ~isempty(i1) 
   if  ~ischar(out{j1,2});
    result{j1,2}=norm(out{j1,3}-k{i1,3});
   else
    result{j1,2}='Error';
   end
  else
    result{j1,2}='not compared';
  end
end

end


%---------------------------------------------------------
function a=eval_disp(st,Cam);

if nargin==1;silent=0;else;silent=~isempty(Cam)&&Cam(end)==';';end
if ~silent;fprintf('Begin... : %s  ',st); end
try; 
  evalin('caller',st); evalin('caller','stu='''';');
  if ~silent;fprintf('%s Finished...',st); end
catch;evalin('caller','stu=''Failed'';');fprintf('Failed :%s %s',st,lasterr); 
end
if ~silent;fprintf(' End\n'); end

%---------------------------------------------------------
