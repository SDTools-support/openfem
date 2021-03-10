function [out,out1,out2]=q25p(CAM,varargin)

%H125 element function : 5x5x5 brick element
%
% 	See sdtweb     eltfun, elem0
%	See also help  t3p, ...

%       Etienne Balmes  
%       Copyright (c) 2001-2018 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

%#ok<*NOSEM,*ASGLU,*NASGU>

% standard calls with one input argument
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1);
 % Build Constit, Integ, and Elmap for later integration
 if comstr(Cam,'integinfo')
  %constit integ,elmap                 ID,pl,il
  [out,out1,out2]=p_solid('buildconstit3 125',varargin{:});

 elseif comstr(Cam,'matcall'); [out,out1]=elem0(CAM,varargin{:});
 elseif comstr(Cam,'dofcall');out=elem0('dofcall'); % variable field elements
 elseif comstr(Cam,'call')
   out='[k1,m1]=h125(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,elmap);';
 elseif comstr(Cam,'rhscall') % call for load assembly
   out='rhs_of';

 elseif comstr(Cam,'groupinit');out=elem0(CAM,'h125');
 % Here one defines the constants needed for element integration
 elseif comstr(Cam,'constants');

  if nargin<3; p_solid('constsolid','h125',[],[]);return;
  elseif varargin{2}(size(varargin{2},1),1)==-9999; % old of_mk_sub.c elements
    i2=125*3;i1=find(triu(ones(i2,i2)));r1=zeros(i2,i2);r1(i1)=1:length(i1);
    out=struct('ElMap',r1+tril(r1',-1));   out2=[];
    out1=varargin{1};out1(4,:)=2;% Tell of_mk('MatrixIntegration') this is 2d
  else; 
   [out,i2,out2]=p_solid('constsolid','h125',varargin{2:end});
   out1=varargin{1};out1(4,:)=i2; % Tell MatrixIntegration this is nd 
  end

elseif  comstr(Cam,'rule') %% #Rule integration rule
 out=feval(q16p('@N_Nr_Ns_Nt'),5);

  elseif  comstr(Cam,'node');  out = 1:125;
  elseif  comstr(Cam,'prop');  out = 125+(1:3); 
  elseif  comstr(Cam,'dof')   
    out = (1:125);k=[1;2;3]/100;out=out(ones(3,1),:)+k(:,ones(125,1));
    out=out(:);
   if 1==2 % Display
      model=feutil('object hexa',[0 0 0;eye(3)],4,4,4);model.Elt=feutil('orient',model);
      cf=feplot(model);
      cf.sel='seledge all';i2=unique(sort(cf.sel.f2,2),'rows');i2(:,3)=0;sprintf('%i ',i2')
      cf.sel='selface all';i2=unique(cf.sel.fs,'rows');
      sprintf('%i %i %i %i; %i %i %i %i;%i %i %i %i;%i %i %i %i;\n',cf.sel.Node(cf.sel.fs'))
      
   end
  elseif  comstr(Cam,'line');   out = [ 1 2 0 1 6 0 1 26 0 2 3 0 2 7 0 2 27 0 ...
    3 4 0 3 8 0 3 28 0 4 5 0 4 9 0 4 29 0 5 10 0 5 30 0 6 7 0 6 11 0 6 31 0 ...
    7 8 0 7 12 0 7 32 0 8 9 0 8 13 0 8 33 0 9 10 0 9 14 0 9 34 0 10 15 0 ...
    10 35 0 11 12 0 11 16 0 11 36 0 12 13 0 12 17 0 12 37 0 13 14 0 13 18 ...
    0 13 38 0 14 15 0 14 19 0 14 39 0 15 20 0 15 40 0 16 17 0 16 21 0 ...
    16 41 0 17 18 0 17 22 0 17 42 0 18 19 0 18 23 0 18 43 0 19 20 0 19 24 ...
    0 19 44 0 20 25 0 20 45 0 21 22 0 21 46 0 22 23 0 22 47 0 23 24 ...
    0 23 48 0 24 25 0 24 49 0 25 50 0 26 27 0 26 31 0 26 51 0 27 28 ...
    0 27 32 0 27 52 0 28 29 0 28 33 0 28 53 0 29 30 0 29 34 0 29 54 ...
    0 30 35 0 30 55 0 31 32 0 31 36 0 31 56 0 32 33 0 32 37 0 32 57 ...
    0 33 34 0 33 38 0 33 58 0 34 35 0 34 39 0 34 59 0 35 40 0 35 60 ...
    0 36 37 0 36 41 0 36 61 0 37 38 0 37 42 0 37 62 0 38 39 0 38 43 ...
    0 38 63 0 39 40 0 39 44 0 39 64 0 40 45 0 40 65 0 41 42 0 ...
    41 46 0 41 66 0 42 43 0 42 47 0 42 67 0 43 44 0 43 48 0 43 68 0 ...
    44 45 0 44 49 0 44 69 0 45 50 0 45 70 0 46 47 0 46 71 0 47 48 0 47 72 ...
    0 48 49 0 48 73 0 49 50 0 49 74 0 50 75 0 51 52 0 51 56 0 51 76 0 52 53 ...
    0 52 57 0 52 77 0 53 54 0 53 58 0 53 78 0 54 55 0 54 59 0 54 79 0 55 60 0 ...
    55 80 0 56 57 0 56 61 0 56 81 0 57 58 0 57 62 0 57 82 0 58 59 0 58 63 ...
    0 58 83 0 59 60 0 59 64 0 59 84 0 60 65 0 60 85 0 61 62 0 61 66 0 61 86 ...
    0 62 63 0 62 67 0 62 87 0 63 64 0 63 68 0 63 88 0 64 65 0 64 69 0 64 89 ...
    0 65 70 0 65 90 0 66 67 0 66 71 0 66 91 0 67 68 0 67 72 0 67 92 0 68 69 ...
    0 68 73 0 68 93 0 69 70 0 69 74 0 69 94 0 70 75 0 70 95 0 71 72 0 71 96 0 ...
    72 73 0 72 97 0 73 74 0 73 98 0 74 75 0 74 99 0 75 100 0 76 77 0 76 81 0 ...
    76 101 0 77 78 0 77 82 0 77 102 0 78 79 0 78 83 0 78 103 0 79 80 0 79 ...
    84 0 79 104 0 80 85 0 80 105 0 81 82 0 81 86 0 81 106 0 82 83 0 82 87 ...
    0 82 107 0 83 84 0 83 88 0 83 108 0 84 85 0 84 89 0 84 109 0 85 90 ...
    0 85 110 0 86 87 0 86 91 0 86 111 0 87 88 0 87 92 0 87 112 0 88 89 ...
    0 88 93 0 88 113 0 89 90 0 89 94 0 89 114 0 90 95 0 90 115 0 91 92 0 ...
    91 96 0 91 116 0 92 93 0 92 97 0 92 117 0 93 94 0 93 98 0 93 118 0 ...
    94 95 0 94 99 0 94 119 0 95 100 0 95 120 0 96 97 0 96 121 0 97 98 ...
    0 97 122 0 98 99 0 98 123 0 99 100 0 99 124 0 100 125 0 101 102 0 ...
    101 106 0 102 103 0 102 107 0 103 104 0 103 108 0 104 105 0 104 109 0 ...
    105 110 0 106 107 0 106 111 0 107 108 0 107 112 0 108 109 0 108 113 0 ...
    109 110 0 109 114 0 110 115 0 111 112 0 111 116 0 112 113 0 112 117 ...
    0 113 114 0 113 118 0 114 115 0 114 119 0 115 120 0 116 117 0 116 121 ...
    0 117 118 0 117 122 0 118 119 0 118 123 0 119 120 0 119 124 0 120 125 ...
    0 121 122 0 122 123 0 123 124 0 124 125 0 ];

  elseif  comstr(Cam,'patch')
    out=[1 2 7 6; 1 26 27 2;1 6 31 26;2 3 8 7;
2 27 28 3; 3 4 9 8;3 28 29 4;4 5 10 9;
4 29 30 5; 5 30 35 10;6 7 12 11;6 11 36 31;
7 8 13 12; 8 9 14 13;9 10 15 14;10 35 40 15;
11 12 17 16; 11 16 41 36;12 13 18 17;13 14 19 18;
14 15 20 19; 15 40 45 20;16 17 22 21;16 21 46 41;
17 18 23 22; 18 19 24 23;19 20 25 24;20 45 50 25;
21 22 47 46; 22 23 48 47;23 24 49 48;24 25 50 49;
26 51 52 27; 26 31 56 51;27 52 53 28;28 53 54 29;
29 54 55 30; 30 55 60 35;31 36 61 56;35 60 65 40;
36 41 66 61; 40 65 70 45;41 46 71 66;45 70 75 50;
46 47 72 71; 47 48 73 72;48 49 74 73;49 50 75 74;
51 76 77 52; 51 56 81 76;52 77 78 53;53 78 79 54;
54 79 80 55; 55 80 85 60;56 61 86 81;60 85 90 65;
61 66 91 86; 65 90 95 70;66 71 96 91;70 95 100 75;
71 72 97 96; 72 73 98 97;73 74 99 98;74 75 100 99;
76 101 102 77; 76 81 106 101;77 102 103 78;78 103 104 79;
79 104 105 80; 80 105 110 85;81 86 111 106;85 110 115 90;
86 91 116 111; 90 115 120 95;91 96 121 116;95 120 125 100;
96 97 122 121; 97 98 123 122;98 99 124 123;99 100 125 124;
101 106 107 102; 102 107 108 103;103 108 109 104;104 109 110 105;
106 111 112 107; 107 112 113 108;108 113 114 109;109 114 115 110;
111 116 117 112; 112 117 118 113;113 118 119 114;114 119 120 115;
116 121 122 117; 117 122 123 118;118 123 124 119;119 124 125 120;
];
  elseif  comstr(Cam,'edge');   out = [ 1 5 6 7 2;2 12 17 22 3;
      3 25 24 23 4 ; 4 18 13 8 1 ];
  elseif  comstr(Cam,'face');   
    out=h125('patch');out=out(:,[1 4 3 2]);
  elseif  comstr(Cam,'sci_face'); out = 1:4;
  elseif  comstr(Cam,'parent'); out = 'h125'; 
  elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,5); 
   %% #TestMat
   if 1==2
   %% #TestGen
     model=femesh('testquad4 divide 4 4'); 
     model.Node(:,5:6)=model.Node(:,5:6)*2-1;
     comstr(model.Node(:,5:7),-30,struct('NoClip',2))
     r1=feutil('getPatch',model);comstr(r1(2:end,1:4),-30,struct('NoClip',1))
   end
   
   model=femesh('testh125');
   if nargin==2; model.pl=varargin{1}; end
   [m,k,mdof]=fe_mknl(model,'NoT');
   out=stack_cell(k,m);
   disp('TestMat passed');
  
  elseif  comstr(Cam,'test');  [CAM,Cam] = comstr(CAM,5);

   if comstr(Cam,'mesh');[CAM,Cam] = comstr(CAM,5);
    model=femesh(['testhexa8' CAM]); 
    out=q16p('h8Toh125',model);
   elseif nargin==2 % specified pl
     [out,out1]=femesh(strcat(['teststruct h125' Cam],varargin{1}));
   else
     [out,out1]=femesh(strcat(['teststruct h125' Cam]));
   end

 % Basic tests of the element - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(CAM,'cvs')
  out='$Revision: 1.4 $  $Date: 2020/06/11 16:29:23 $';
 else;sdtw('''%s'' unknown',CAM);  
 end

  return
end % of standard calls with one input argument

% -----------------------------------------------------------------------------
% element matrix assembly - - - - - - - - - - - - - - - - - - - - - - - - - -

error('This is only implemented as a mat_og element');
