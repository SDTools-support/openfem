if ~isempty(gcbf)&&strcmp(get(gcbf,'tag'),'test_sdt');return;end


if 1==2
%% #Generate data for testing    
cd G:/balmes/sdt.cur;sdtcheck path
cd G:\balmes\bug

%addpath ('G:\balmes\bug\of_mk\x64\Debug\')
addpath ('G:\balmes\bug\of_mk\Debug\')
addpath ('C:\Users\ial\Downloads\comp_struct')
wd1='G:/balmes/bug/';

% For lounes : test that fails
z=load(comgui('cd',fullfile(wd1,'bug_omp_xxx1.mat')));
mo2=z.mo2;r1.matdes=1;
of_mk('setomppro',10)
mo2.K={}; mo2=fe_case(sprintf('assemble -matdes %i -se -NoT',r1.matdes),mo2);
%sd('compare',mo2,z.mo2)
sdthdf('compare',mo2,z.mo2)
%[match, er1, er2, erc] = comp_struct(mo2,z.mo2,0,0,1e-6,'a','b');
% of_mk('setomppro',2)
% mo2.K={}; mo2=fe_case(sprintf('assemble -matdes %i -se -NoT',r1.matdes),mo2);
% sdthdf('compare',mo2,z.mo2)
%[match, er1, er2, erc] = comp_struct(mo2,z.mo2,0,0,1e-6,'a','b');


clc;
clear all;
cd G:/balmes/sdt.cur;sdtcheck path
cd G:/balmes/bug

%addpath ('f:\balmes\bug\of_mk\x64\Debug\')
addpath ('G:\balmes\bug\of_mk\Debug\')

%% Generate the data for call to of_mk_exe
wd1='G:/balmes/bug/';

sdtweb fe_mknl matrixinteg
basic_elt_test mat q4pb
EvalFcn=['of_mk(''matrixintegration'',DofPos,NodePos,Case.Node,' ...
         'pointers,integ,constit,gstate,' ...
        'elmap,InfoAtNode,EltConst,def.def,' ...
         'k,int32([Case.DofPerElt(jGroup);SymFlag;0]));'];
cd D:\share\sdtdata\Lounes; save('q4pb.mat'); 

% load('bug_omp_call.mat');
% cd D:\share\sdtdata\Lounes;load q4pb.mat

varg={'matrixintegration',DofPos,NodePos,Case.Node, ...
    pointers,integ,constit,gstate,elmap,InfoAtNode,EltConst,def.def, ...
    k,int32([Case.DofPerElt(jGroup);SymFlag;0])};
save data varg

of_mk('setomppro',5); of_mk(varg{:})

load('bug_omp_xxx1.mat');
sdthdf('compare',k,mo2.K{1,1})

%sd of_mk-omp-debug
end



if 1==2 % test d'Etienne

 try;sdtcheck('pathnone');end
 cd /scratch/openfem;ofutil path;cd test;
%ofutil('of_mk-DOFOMP -f ../test/mexoptions ')
of_mk setomppro
of_mk('setomppro',2)

model=femesh('testhexa8b divide 2 2 2');
dbstop if caught error

tic;fe_mknl(model);toc 


model=femesh('testhexa8b divide 100 10 10');
of_mk('setomppro',1);tic;fe_mknl(model);toc
of_mk('setomppro',2);tic;fe_mknl(model);toc


end

