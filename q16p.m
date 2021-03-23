function [out,out1,out2]=q16p(CAM,varargin)

%Q16P element function 16 node plane element (4x4 grid)
%   Implements methods for spectral elements
%
% 	See sdtweb     eltfun, elem0
%	See also help  t3p, ...

%       Etienne Balmes  
%       Copyright (c) 2001-2020 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

%#ok<*NOSEM,*ASGLU,*NASGU>

% standard calls with one input argument
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1);
 % Build Constit, Integ, and Elmap for later integration
 if comstr(Cam,'integinfo')
  %constit integ,elmap                 ID,pl,il
  [out,out1,out2]=p_solid('buildconstit2 16',varargin{:});

 elseif comstr(Cam,'matcall'); [out,out1]=elem0(CAM,varargin{:});
 elseif comstr(Cam,'dofcall');out=elem0('dofcall'); % variable field elements
 elseif comstr(Cam,'call')
   out='[k1,m1]=q16p(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,elmap);';
 elseif comstr(Cam,'rhscall') % call for load assembly
   out='rhs_of';

 elseif comstr(Cam,'groupinit');out=elem0(CAM,'q16p');
 % Here one defines the constants needed for element integration
 elseif comstr(Cam,'constants');

  if nargin<3; p_solid('constsolid','q16p',[],[]);return;
  elseif varargin{2}(size(varargin{2},1),1)==-9999; % old of_mk_sub.c elements
    i2=16;i1=find(triu(ones(i2,i2)));r1=zeros(i2,i2);r1(i1)=1:length(i1);
    out=struct('ElMap',r1+tril(r1',-1));   out2=[];
    out1=varargin{1};out1(4,:)=2;% Tell of_mk('MatrixIntegration') this is 2d
  else; 
   [out,i2,out2]=p_solid('constsolid','q16p',varargin{2:end});
   out1=varargin{1};out1(4,:)=i2; % Tell MatrixIntegration this is nd 
  end

elseif  comstr(Cam,'rule')  % #Rule integration rule
 out=feval(t_spectral('@N_Nr_Ns'),4);
 
  elseif  comstr(Cam,'node');  out = 1:16;
  elseif  comstr(Cam,'prop');  out = (17:19); 
  elseif  comstr(Cam,'dof')   
    out = (1:16);k=[1;2]/100;out=out(ones(2,1),:)+k(:,ones(8,1));
    out=out(:);
  elseif  comstr(Cam,'line');  out = [ 1 7 11 4 15 12 8 5 6 9 ...
     13 16 3 14 10 2 6 0 5 1 0 7 8 9 10 0 14 13 12 11 0 15 16 ];
  elseif  comstr(Cam,'patch'); 
    out=[1 2 6 5;2 3 7 6;3 4 8 7;
         5 6 10 9;6 7 11 10;7 8 12 11;
         9 10 14 13;10 11 15 14;11 12 16 15];
  elseif  comstr(Cam,'edge');   out = [ 1 2 3 4;4 8 12 16;16 15 14 13;13 9 5 1];
  elseif  comstr(Cam,'face');   out = [1 4 16 14];
  elseif  comstr(Cam,'sci_face'); out = [1 4 16 14];
  elseif  comstr(Cam,'parent'); out = 'q16p'; 
  elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,5); 
   %% #TestMat
   if 1==2
   %% #TestGen
     model=femesh('testquad4 divide 3 3'); 
     model.Node(:,5:6)=model.Node(:,5:6)*2-1;
     comstr(model.Node(:,5:7),-30,struct('NoClip',2))
     r1=feutil('getPatch',model);comstr(r1(2:end,1:4),-30,struct('NoClip',1))
   end
   
   model=femesh('testq16p');
   if nargin==2; model.pl=varargin{1}; end
   if ~isempty(strfind(Cam,'_0'));     model.il(:,3)=0;
   elseif ~isempty(strfind(Cam,'_1')); model.il(:,3)=1;
   elseif ~isempty(strfind(Cam,'_2')); model.il(:,3)=2;
   end
   [m,k,mdof]=fe_mknl(model,'NoT');
   out=stack_cell(k,m);
   disp('TestMat passed');

  elseif  comstr(Cam,'q4to');  [CAM,Cam] = comstr(CAM,5);
  %% #Q4To(Elt) mesh spectral element from from quad4
  carg=1;model=varargin{carg};carg=carg+1; 
  if carg<=nargin-1; RO=varargin{carg};carg=carg+1;
  else; RO=struct('EC',integrules(CAM));
  end
  
  RO.base='q4p'; RO.iProp=feval(RO.base,'prop');
  RO.iNode=feval(RO.base,'node'); RO.Nnode=length(RO.iNode);
  
  [EGroup,nGroup]=getegroup(model.Elt);
  jGroup=1;   cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;
  r2=RO.EC.xi; r2(:,4)=0; r2=integrules(RO.base,r2);% New points
  r3=[reshape(repmat(model.Elt(cEGI,RO.iNode)',size(r2.N,1),1),RO.Nnode,[]);
      repmat(r2.N,length(cEGI),1)'];
  
  % new node given by [[NodeId Coef]x corners]
  % sort by NodeId to avoid repetition
  r4=r3';r4(:,RO.Nnode+1:end)=round(r4(:,RO.Nnode+1:end)*1e5);
  for j1=1:size(r4,1); 
    r5=r4(j1,1:RO.Nnode); r6=r4(j1,RO.Nnode+1:end); r5(r6==0)=0; % Unused weight for edges
    [r5,i2]=sort(r5);r6=r6(i2);
     r4(j1,:)=[r5 r6];
  end
  [r4,i1,i3]=unique(r4,'rows');
  % xxx Missing preservation of existing nodes with one non-zero in r4 row
  NNode=sparse(model.Node(:,1),1,1:size(model.Node));
  c=sparse(repmat(1:length(i1),4,1),r3(1:4,i1),r3(5:8,i1));
  n1=[(max(model.Node(:,1)*0)+(1:size(c,1))')*[1 0 0 0] c*model.Node(:,5:7)];
  model.Node=n1; %[model.Node;n1];  
  
  % Now deal with new elements
  model.Elt=feutil('addelt',sprintf('q%ip',size(RO.EC.xi,1)), ...
      [reshape(n1(i3,1),size(r2.N,1),[])'  model.Elt(cEGI,RO.iProp)]);
  out=model;
  
  elseif  comstr(Cam,'h8to');  [CAM,Cam] = comstr(CAM,5);
  %% #h8To(Elt) mesh spectral element from from hexa8
  [CAM,Cam,v]=comstr('v',[-25 31],CAM,Cam); 
  carg=1;model=varargin{carg};carg=carg+1; 
  if carg<=nargin-1; RO=varargin{carg};carg=carg+1;
  else; 
   if strcmpi(CAM,'h27');
     CAM='hexa27';ElemF=CAM;RO=struct('EC',integrules(CAM,-2));
   else;RO=struct('EC',integrules(CAM));ElemF=CAM;
   end
  end
  
  RO.base='hexa8'; RO.iProp=feval(RO.base,'prop');
  RO.iNode=feval(RO.base,'node'); RO.Nnode=length(RO.iNode);
  
if v~=0
  r2=feval(feutil('@isoCellRef'),'hexa8',struct('EC',RO.EC));
  r2=feutil('RefineCell-given -keepEP',model,struct('hexa8',r2));
  model.Node=r2.Node;model.Elt=r2.Elt;
else
  [EGroup,nGroup]=getegroup(model.Elt);
  jGroup=1;   cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;
  r2=RO.EC.xi; r2(:,4)=0; r2=integrules(RO.base,r2);% New points
  r3=[reshape(repmat(model.Elt(cEGI,RO.iNode)',size(r2.N,1),1),RO.Nnode,[]);
      repmat(r2.N,length(cEGI),1)'];
  
  % new node given by [[NodeId Coef]x corners]
  % sort by NodeId to avoid repetition
  r4=r3';r4(:,RO.Nnode+1:end)=round(r4(:,RO.Nnode+1:end)*1e5);
  for j1=1:size(r4,1); 
    r5=r4(j1,1:RO.Nnode); r6=r4(j1,RO.Nnode+1:end); r5(r6==0)=0; % Unused weight for edges
[r5,i2]=sort(r5);r6=r6(i2);
     r4(j1,:)=[r5 r6];
  end
  [r4,i1,i3]=unique(r4,'rows');
  % xxx Missing preservation of existing nodes with one non-zero in r4 row
  NNode=sparse(model.Node(:,1),1,1:size(model.Node));
  c=sparse(repmat(1:length(i1),RO.Nnode,1),r3(1:RO.Nnode,i1), ...
      r3(RO.Nnode+1:2*RO.Nnode,i1));
  n1=[(max(model.Node(:,1)*0)+(1:size(c,1))')*[1 0 0 0] c*model.Node(:,5:7)];
  model.Node=n1; %[model.Node;n1];  
  model.Elt(1,end+1:max(RO.iProp))=0;
  % Now deal with new elements
  model.Elt=feutil('addelt',ElemF, ...
      [reshape(n1(i3,1),size(r2.N,1),[])'  model.Elt(cEGI,RO.iProp)]);
end
  out=model;
elseif  comstr(Cam,'toh8'); 
%% #h(d)ToH8 : obtain hexa8 mesh 
 model=varargin{1};
 [EGroup,nGroup,names]=getegroup(model.Elt);
 elt=[];
  if 1==2
   mo1=femesh('testhexa8 divide 4 4 4');
   %mo1=femesh('testhexa8 divide 3 3 3');
   mo1.Node(:,5:7)=round(mo1.Node(:,5:7)*1000);
   [un1,i1]=sortrows(mo1.Node,[7 6 5]);i1(:,2)=(1:length(i1))';
   mo1=feutil('renumber',mo1,i1);feplot(mo1);fecom('textnode');
   comstr(mo1.Elt(2:end,1:8),-30,struct('NoClip',2))
  end
 for jGroup=1:nGroup
  [ElemF,i1,ElemP]= getegroup(model.Elt(EGroup(jGroup),:),jGroup);
  cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;

  if strcmpi(ElemP,'h125')
   i1=[1 2 7 6 26 27 32 31;2 3 8 7 27 28 33 32;3 4 9 8 28 29 34 33;
       4 5 10 9 29 30 35 34;6 7 12 11 31 32 37 36;7 8 13 12 32 33 38 37;
       8 9 14 13 33 34 39 38;9 10 15 14 34 35 40 39;11 12 17 16 36 37 42 41;
       12 13 18 17 37 38 43 42;13 14 19 18 38 39 44 43;14 15 20 19 39 40 45 44;
       16 17 22 21 41 42 47 46;17 18 23 22 42 43 48 47;18 19 24 23 43 44 49 48;
       19 20 25 24 44 45 50 49;26 27 32 31 51 52 57 56;27 28 33 32 52 53 58 57;
       28 29 34 33 53 54 59 58;29 30 35 34 54 55 60 59;31 32 37 36 56 57 62 61;
       32 33 38 37 57 58 63 62;33 34 39 38 58 59 64 63;34 35 40 39 59 60 65 64;
       36 37 42 41 61 62 67 66;37 38 43 42 62 63 68 67;38 39 44 43 63 64 69 68;
       39 40 45 44 64 65 70 69;41 42 47 46 66 67 72 71;42 43 48 47 67 68 73 72;
       43 44 49 48 68 69 74 73;44 45 50 49 69 70 75 74;51 52 57 56 76 77 82 81;
       52 53 58 57 77 78 83 82;53 54 59 58 78 79 84 83;54 55 60 59 79 80 85 84;
       56 57 62 61 81 82 87 86;57 58 63 62 82 83 88 87;58 59 64 63 83 84 89 88;
       59 60 65 64 84 85 90 89;61 62 67 66 86 87 92 91;62 63 68 67 87 88 93 92;
       63 64 69 68 88 89 94 93;64 65 70 69 89 90 95 94;66 67 72 71 91 92 97 96;
       67 68 73 72 92 93 98 97;68 69 74 73 93 94 99 98;69 70 75 74 94 95 100 99;
       76 77 82 81 101 102 107 106;77 78 83 82 102 103 108 107;
       78 79 84 83 103 104 109 108;79 80 85 84 104 105 110 109;
       81 82 87 86 106 107 112 111;82 83 88 87 107 108 113 112;
       83 84 89 88 108 109 114 113;84 85 90 89 109 110 115 114;
       86 87 92 91 111 112 117 116;87 88 93 92 112 113 118 117;
       88 89 94 93 113 114 119 118;89 90 95 94 114 115 120 119;
       91 92 97 96 116 117 122 121;92 93 98 97 117 118 123 122;
       93 94 99 98 118 119 124 123;94 95 100 99 119 120 125 124];
    i2=126:127;
  elseif strcmpi(ElemP,'h64'); 
    i1=[1 2 6 5 17 18 22 21;2 3 7 6 18 19 23 22;3 4 8 7 19 20 24 23;
        5 6 10 9 21 22 26 25;6 7 11 10 22 23 27 26;7 8 12 11 23 24 28 27;
        9 10 14 13 25 26 30 29;10 11 15 14 26 27 31 30;11 12 16 15 27 28 32 31;
        17 18 22 21 33 34 38 37;18 19 23 22 34 35 39 38;19 20 24 23 35 36 40 39;
        21 22 26 25 37 38 42 41;22 23 27 26 38 39 43 42;23 24 28 27 39 40 44 43;
        25 26 30 29 41 42 46 45;26 27 31 30 42 43 47 46;27 28 32 31 43 44 48 47;
        33 34 38 37 49 50 54 53;34 35 39 38 50 51 55 54;35 36 40 39 51 52 56 55;
        37 38 42 41 53 54 58 57;38 39 43 42 54 55 59 58;39 40 44 43 55 56 60 59;
        41 42 46 45 57 58 62 61;42 43 47 46 58 59 63 62;43 44 48 47 59 60 64 63];
    i2=65:66;
  elseif strcmpi(ElemP,'hexa27'); i1=hexa27('vollin');i2=28:29;
  else;
      i1=model.Elt(EGroup(jGroup):EGroup(jGroup+1)-1,:);i2=[];
  end
  if ~isempty(i2)% Now split elements
    i3=model.Elt(cEGI,feval(ElemP,'nodes'))';
    i3=reshape(i3(i1',:),8,[])';
    i2=reshape(repmat(model.Elt(cEGI,i2)',size(i1,1),1),length(i2),[]);
    elt=feutil('addelt',elt,'hexa8',[i3 i2']);
  else;
   elt(end+(1:size(i1,1)),1:size(i1,2))=i1;
  end 
 end
 out=model;out.Elt=elt;
 
  elseif  comstr(Cam,'test');  [CAM,Cam] = comstr(CAM,5);

   if comstr(Cam,'mesh');[CAM,Cam] = comstr(CAM,5);
    model=femesh(['testq4p' CAM]); 
    out=q16p('q4Toq16p',model);
   elseif nargin==2 % specified pl
     [out,out1]=femesh(strcat(['teststruct q16p' Cam],varargin{1}));
   else
     [out,out1]=femesh(strcat(['teststruct q16p' Cam]));
   end

 % Basic tests of the element - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(CAM,'@');out=eval(CAM);
 elseif comstr(CAM,'cvs')
  out='$Revision: 1.13 $  $Date: 2020/06/19 12:41:39 $';
 else sdtw('''%s'' unknown',CAM);  
 end

  return
end % of standard calls with one input argument

% -----------------------------------------------------------------------------
% element matrix assembly - - - - - - - - - - - - - - - - - - - - - - - - - -

error('This is only implemented as a mat_og element');

end



%% Legendre polynomials
function y=polylegendre(n)
m=n+2;
P=zeros(m);
P(1,:)=[1 zeros(1,m-1)];
P(2,:)=[0 1 zeros(1,m-2)];
for i=3:m
    Pi=zeros(1,m); Pi(1)=-(i-2)*P(i-2,1)/(i-1);
    for j=2:m-1
        Pj=((2*i-3)*P(i-1,j-1)-(i-2)*P(i-2,j))/(i-1);
        Pi(j)=Pj;
    end
    P(i,:)=Pi;
end
y=fliplr(P(n+1,1:n+1));
end

function y=coord(n) 
%% #coord(n) % Roots of (1) in see Sylvere MILLAN 2016
n=n-1;
P=polyder(polylegendre(n));
A=[1 0 -1];
y=sort(roots(conv(A,P)));
end

function [pN,pdN]=forme(ind,n)
% #forme(N,n);
lambda=coord(n);
C=zeros(n,n);
for i=1:n
    B=[zeros(1,n) 1];
for j=1:n
    if i==j;A=1;
    else;A=[1 -lambda(j)]./(lambda(i)-lambda(j));
    end
    B=conv(B,A);
end
B=B(length(B)-(n-1):end);
C(i,:)=B;
end
pN=C(ind,:);
if nargout==2
 for j1=1:size(pN,1)
  pdN(j1,1:size(pN,1)-1)=polyder(pN(j1,:));
 end
end
end

function EC=N_Nr_Ns(n)
% #N_Nr_Ns : n number of nodes in each direction

EC=struct;
pos=meshgrid(coord(n),coord(n));
xi=[reshape(pos',[],1) reshape(pos,[],1)*[1 0]];
EC=struct('N',[],'Nr',[],'Ns',[], ...
      'Nw',n^2,'NDN',[], ...
      'NDNLabels',{{'',',x',',y'}}, ...
      'jdet',[],'w',xi,'Nnode',[],'xi',xi, ...
      'type',sprintf('q%ip',size(xi,1)));


%% N N_(ij)(r_(ij),s_(ij))=N_ 2D Shapes functions. 
B=zeros(n^2,length(forme(1,n)));
%NR 
A=zeros(n^2,length(forme(1,n)));
B=zeros(n^2,length(forme(1,n)));
C=forme(1:n,n);
D=zeros(n,length(forme(1,n)));
for i=1:n
    D(i,2:length(forme(1,n)))=polyder(forme(i,n));
end
for i=1:n
    A(n*(i-1)+1:n*i,:)=D;
    for j=1:n+1
        B((i-1)*n+j,:)=C(i,:);
    end
end
E=zeros(n^2);
for i=1:n^2
    for k=1:n^2
        E(i,k)=polyval(A(k,:),xi(i,1))*polyval(B(k,:),xi(i,2));
    end
end
EC.Nr=E;

% Ns
A=zeros(n^2,length(forme(1,n)));
B=zeros(n^2,length(forme(1,n)));
C=forme(1:n,n);
D=zeros(n,length(forme(1,n)));
for i=1:n
    D(i,2:length(forme(1,n)))=polyder(forme(i,n));
end
for i=1:n
    A(n*(i-1)+1:n*i,:)=C;
    for j=1:n
        B((i-1)*n+j,:)=D(i,:);
    end
end
E=zeros(n^2);
for i=1:n^2
    for k=1:n^2
        E(i,k)=polyval(A(k,:),xi(i,1))*polyval(B(k,:),xi(i,2));
    end
end
EC.Ns=E; 

% W
A=polyval(polylegendre(n-1),xi(:,1));
B=2./(n*(n-1)*A.*A);
C=polyval(polylegendre(n-1),xi(:,2));
D=2./(n*(n-1)*C.*C);
EC.w(:,4)=B.*D;

if any(~isfinite(EC.w(:)));
    fprintf('Problem with weights');
    dbstack;keyboard;
end
EC.Nr(abs(EC.Nr)<1e-15)=0;
EC.Ns(abs(EC.Ns)<1e-15)=0;

EC.Nnode=size(EC.xi,1); EC.jdet=zeros(EC.Nnode,1);
EC.NDN=zeros(size(EC.N,1),size(EC.N,2)*3);

end


function EC=N_Nr_Ns_Nt(n) %#ok<DEFNU>
% #N_Nr_Ns_Nt : feval(q16p('@N_Nr_Ns_Nt'),3)

EC=struct;
pos=meshgrid(coord(n),coord(n));
xi=[reshape(repmat(reshape(pos',[],1),1,n),[],1) ...
    reshape(repmat(reshape(pos,[],1),1,n),[],1) ...
    reshape(ones(n^2,1)*coord(n)',[],1)];
EC=struct('N',[],'Nr',[],'Ns',[],'Nt',[], ...
      'Nw',n^3,'NDN',[], ...
      'NDNLabels',{{'',',x',',y','z',}}, ...
      'jdet',[],'w',xi,'Nnode',[],'xi',xi, ...
      'type',sprintf('h%i',size(xi,1)));
if strcmpi(EC.type,'h27');EC.type='hexa27';end
  % W
PL=polylegendre(n-1);
EC.w(:,4)=2./(n*(n-1)*polyval(PL,xi(:,1)).^2) .* ...
 2./(n*(n-1)*polyval(PL,xi(:,2)).^2) .* ...
 2./(n*(n-1)*polyval(PL,xi(:,3)).^2) ;
[N1,N1r]=forme(1:n,n);

for j1=1:n
for j2=1:n
for j3=1:n
  ind=1+(j1-1)+(j2-1)*n+(j3-1)*n^2;
  EC.N(:,ind)=polyval(N1(j1,:),xi(:,1))  .*polyval(N1(j2,:),xi(:,2)) .*polyval(N1(j3,:),xi(:,3));
  EC.Nr(:,ind)=polyval(N1r(j1,:),xi(:,1)).*polyval(N1(j2,:),xi(:,2)) .*polyval(N1(j3,:),xi(:,3));
  EC.Ns(:,ind)=polyval(N1(j1,:),xi(:,1)) .*polyval(N1r(j2,:),xi(:,2)).*polyval(N1(j3,:),xi(:,3));
  EC.Nt(:,ind)=polyval(N1(j1,:),xi(:,1)) .*polyval(N1(j2,:),xi(:,2)) .*polyval(N1r(j3,:),xi(:,3));
end
end
end
if norm(EC.N-eye(size(EC.N)))>1e-10;error('mismatch');end
EC.N=eye(size(EC.N));

EC.Nr(abs(EC.Nr)<1e-15)=0;
EC.Ns(abs(EC.Ns)<1e-15)=0;
EC.Nt(abs(EC.Nt)<1e-15)=0;
%norm(EC.Nr*xi-ones(size(xi,1),1)*[1 0 0],'inf')
%norm(EC.Ns*xi-ones(size(xi,1),1)*[0 1 0],'inf')
%norm(EC.Nt*xi-ones(size(xi,1),1)*[0 0 1],'inf')

EC.Nnode=size(EC.xi,1); EC.jdet=zeros(EC.Nnode,1);
EC.NDN=zeros(size(EC.N,1),size(EC.N,2)*4);
EC.NDN(:,1:size(EC.N,2))=EC.N;

end


