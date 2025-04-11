

% tests of coherence with values in INRIA contributed elements

list={'tetra4','etr3p1d','node';
 'penta6','etr3r1d','';
 'tetra10','etr3p2c','';
 'penta15','etr3r2c','';
 'tria6','etr2p2c',''}

for j0=1:4;

fid=fopen(fullfile(fileparts(which('ofutil')),'src',[list{j0,2} '.c'])); 
st=fscanf(fid,'%c');fclose(fid);

if ~isempty(list{j0,3}); rule=integrules(list{j0,[1 3]});
else; rule=integrules(list{j0,1});
end
i1=strfind(st,'Local variables');if ~isempty(i1);st=st(1:i1);end


Nnode=rule.Nnode; Nw=rule.Nw; Ndim=3;

i2=strfind(st,'poids'); i2=i2+min(find(st(i2(1):end)=='{'));
w=sscanf(strrep(st(i2(1):end),',',' '),'%g');
Nw=length(w);
i2=strfind(st,'dpx');if ~isempty(i2);
 i2=i2+min(find(st(i2(1):end)=='{'));
 dpx=sscanf(strrep(st(i2(1):end),',',' '),'%g');
 i2=strfind(st,'dpy');i2=i2+min(find(st(i2(1):end)=='{'));
 dpy=sscanf(strrep(st(i2(1):end),',',' '),'%g');
 dp=reshape(dpx,3,6);dp(:,:,2)=reshape(dpy,3,6);
else
 i2=strfind(st,'dp');if isempty(i2); i2=strfind(st,'e_');end
 i2=i2+min(find(st(i2(1):end)=='{'));
  dp=sscanf(strrep(st(i2(1):end),',',' '),'%g');
  dp=reshape(dp,Ndim,Nnode,Nw);
end
i2=strfind(st,'vp');if isempty(i2);i2=strfind(st,'e_');i2=i2(2);end
i2=i2+min(find(st(i2(1):end)=='{'));
[vp,j1,j2,j3]=sscanf(strrep(st(i2(1):end),',',' '),'%g');
vp=reshape(vp,Nnode,Nw)';


if norm(rule.w(:,4)-w)>1e-6; error('weight');end
r1=norm([squeeze(dp(1,:,:))'-rule.Nr])+norm([squeeze(dp(2,:,:))'-rule.Ns])+ ...
 norm([squeeze(dp(3,:,:))'-rule.Nt]);
if norm(r1)>1e-5; error(sprintf('%g DP difference',r1));end

r1=norm(abs(rule.N-vp));
if r1>1e-5; error(sprintf('%g N difference',r1));end

end % Loop on j0


%% Verification of derivatives by finite differences
l2={ 'hexa20','hexa8','tetra4','tetra10','penta6','penta15','quad4','tria3','tria6','quadb'}';
for j0=1:size(l2,1)

 ElemF=l2{j0,1};
 r1=integrules(ElemF);
 r2=integrules(ElemF,r1.w+ones(size(r1.w,1),1)*[1e-5 0 0 0]);
 if norm(((r2.N-r1.N)/1e-5)-r1.Nr,'inf')>1e-4; error('wrong Nr');end
 r2=integrules(ElemF,r1.w+ones(size(r1.w,1),1)*[0 1e-5 0 0]);
 if norm(((r2.N-r1.N)/1e-5)-r1.Ns,'inf')>1e-4; error('wrong Ns');end
 if isfield(r2,'Nt')
  r2=integrules(ElemF,r1.w+ones(size(r1.w,1),1)*[0 0 1e-5 0]);
  if norm(((r2.N-r1.N)/1e-5)-r1.Nt,'inf')>1e-4; error('wrong Nt');end
 end
end

%% Consistence of compilation

 ElemF='tetra10';ElemF(end+1)='b';femesh('reset');
 model=femesh(strcat('teststruct',ElemF,'divide 1 1 1'));

 %model.Elt=model.Elt(1:2,:);
  sp_util('diag',0);[Case,model.DOF]=fe_mknl('init',model);
 def=struct('def',model.DOF*0,'DOF',model.DOF);
 k1=fe_mknl('assemble',model,Case,1);
 m1=fe_mknl('assemble',model,Case,2);
 Case=fe_mknl('init',model);Case.GroupInfo{end}.material='';
 k3a=fe_mknl('assemble',model,Case,def,1);
 Case=fe_mknl('init',model);Case.GroupInfo{end}.material='Elastic3DNL';
 k3=fe_mknl('assemble',model,Case,def,5);

 i1=find(abs(k3)>1e-10*norm(diag(k3),'inf'));[r1,i2]=max(abs(k3(i1)-k3a(i1)));r1=r1/abs(k3(i1(i2)));
 if r1>1e-10; error('Inconsistent compilation');end

 def.def=rand(size(def.def))/10; % consistence of elem0 and of_mk
 z=def.def+0;
 sp_util('diag',0);k3=fe_mknl('assemble',model,Case,def,5);
 sp_util('diag',12);k3a=fe_mknl('assemble',model,Case,def,5);sp_util('diag',0);
 i1=find(abs(k3)>1e-10*norm(diag(k3),'inf'));[r1,i2]=max(abs(k3(i1)-k3a(i1)));r1=r1/abs(k3(i1(i2)));
 if r1>1e-10; error('Inconsistent compilation / elem0');end
 

% Consistence of mat_of and mat_og matrices

st={ 'hexa20','hexa8','tetra4','tetra10','penta6','penta15' };
st1={};
for j1=1:length(st)  

 ElemF=[st{j1} 'b'];femesh('reset');
 model=femesh(strcat('test',ElemF,'struct divide 4 4 4 load'));

 [Case,model.DOF]=fe_mknl('init',model);
 def=struct('def',model.DOF*0,'DOF',model.DOF);
 k1=fe_mknl('assemble',model,Case,1);
 m1=fe_mknl('assemble',model,Case,2);
 k3=fe_mknl('assemble',model,Case,def,5);
 LB=fe_load(model);

 model.Elt(1,1:length(ElemF)+3)=[Inf abs(ElemF) 0 0];
 [C2,model.DOF]=fe_mknl('init',model);
 k2=fe_mknl('assemble',model,C2,1);
 m2=fe_mknl('assemble',model,C2,2);
 L1=fe_load(model);

 st1(end+1,1:2)={ElemF [[norm(k1-k2,'inf') norm(k1-k3,'inf')]/norm(k2,'inf') norm(m1-m2,'inf')/norm(m1,'inf') full(sum(abs(LB.def-L1.def))./sum(abs(LB.def)))]};
 %sum(abs(L1.def-LB.def))
 %cf.def(1)=L1;cf.def(2)=LB;fecom show2def
end

% err:  K(*b/*) K(*b lin/NL) M(*b/*) Gravi, Surf, Pres
st1=st1';fprintf('\n\n');fprintf('%10s %12.5g  %12.5g %12.5g %6.1g %6.1g %6.1g\n',st1{:})


% 2D elements of the b family

st={ 'q4p','q8p','t6p','t3p'};st1=cell(length(st),2);
for j1=1:length(st)  

 ElemF=st{j1}; femesh('reset');
 model=femesh(strcat('teststruct',ElemF,'divide 3 3'));

 model.Elt(1,1:length(ElemF)+3)=[Inf abs(ElemF) 0 0];
 [C2,model.DOF]=fe_mknl('init',model);
 k2=fe_mknl('assemble',model,C2,1);
 m2=fe_mknl('assemble',model,C2,2);
 model.il(:,5)=-3;
 
 if j1==1; model.il(:,5)=-2;end  % modulef uses integration at nodes
 model.Elt(1,1:length(ElemF)+3)=[Inf abs(horzcat(ElemF,'b')) 0];
 [Case,model.DOF]=fe_mknl('init',model);
 def=struct('def',model.DOF*0,'DOF',model.DOF);
 k1=fe_mknl('assemble',model,Case,1);
 m1=fe_mknl('assemble',model,Case,2);

 %[diag(k1)./diag(k2) diag(m1)./diag(m2)];full(ans)

 st1(j1,1:2)={ElemF [[norm(k1-k2,'inf')]/norm(k2,'inf') norm(m1-m2,'inf')/norm(m1,'inf')]};

 %if norm(r2)>1e-10; error('inconsistent formulation');end
end

st1=st1';fprintf('\n\n');fprintf('%10s %10.5g  %10.5g\n',st1{:})

if 1==2
 ElemF='hexa8';ElemF='q4p';
 femesh('reset');
 model=femesh(strcat('teststruct',ElemF,'divide 3 3'));
 %model.Elt(1,1:length(ElemF)+3)=[Inf abs(horzcat(ElemF,'b')) 0];
 [Case,model.DOF,dc]=fe_mknl('init -gstate',model);
 dc.def(:,2)=0;k=fe_mknl('assemble not',model,Case,dc,1);
 sdtdef('diag',12);dc.def(:,2)=0;k=fe_mknl('assemble not',model,Case,dc,1);

end
