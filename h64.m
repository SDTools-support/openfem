function [out,out1,out2]=h64(CAM,varargin)

%h64 element function : 4*4*4 brick element
%
% 	See sdtweb     eltfun, elem0
%	See also help  t3p, ...

%       Etienne Balmes  
%       Copyright (c) 2001-2019 by SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

%#ok<*NOSEM,*ASGLU,*NASGU>

% standard calls with one input argument
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1);
 % Build Constit, Integ, and Elmap for later integration
 if comstr(Cam,'integinfo')
  %constit integ,elmap                 ID,pl,il
  [out,out1,out2]=p_solid('buildconstit3 64',varargin{:});

 elseif comstr(Cam,'matcall'); [out,out1]=elem0(CAM,varargin{:});
 elseif comstr(Cam,'dofcall');out=elem0('dofcall'); % variable field elements
 elseif comstr(Cam,'call')
   out='[k1,m1]=h64(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,elmap);';
 elseif comstr(Cam,'rhscall') % call for load assembly
   out='rhs_of';

 elseif comstr(Cam,'groupinit');out=elem0(CAM,'h64');
 % Here one defines the constants needed for element integration
 elseif comstr(Cam,'constants');

  if nargin<3; p_solid('constsolid','h64',[],[]);return;
  elseif varargin{2}(size(varargin{2},1),1)==-9999; % old of_mk_sub.c elements
    i2=64*3;i1=find(triu(ones(i2,i2)));r1=zeros(i2,i2);r1(i1)=1:length(i1);
    out=struct('ElMap',r1+tril(r1',-1));   out2=[];
    out1=varargin{1};out1(4,:)=2;% Tell of_mk('MatrixIntegration') this is 2d
  else; 
   [out,i2,out2]=p_solid('constsolid','h64',varargin{2:end});
   out1=varargin{1};out1(4,:)=i2; % Tell MatrixIntegration this is nd 
  end

elseif  comstr(Cam,'rule') %% #Rule integration rule
 out=feval(q16p('@N_Nr_Ns_Nt'),4);

  elseif  comstr(Cam,'node');  out = 1:64;
  elseif  comstr(Cam,'prop');  out = 64+(1:3); 
  elseif  comstr(Cam,'dof')   
    out = (1:64);k=[1;2;3]/100;out=out(ones(3,1),:)+k(:,ones(64,1));
    out=out(:);
   if 1==2 % Display
      model=feutil('object hexa',[0 0 0;eye(3)],3,3,3);
      cf=feplot(model);
      cf.sel='seledge all';i2=unique(sort(cf.sel.f2,2),'rows');i2(:,3)=0;sprintf('%i ',i2')
      cf.sel='selface all';i2=unique(cf.sel.fs,'rows');
      sprintf('%i %i %i %i; %i %i %i %i;%i %i %i %i;%i %i %i %i;\n',cf.sel.Node(cf.sel.fs'))
      
   end
  elseif  comstr(Cam,'line');   out = [ 1 2 0 1 5 0 1 17 0 2 3 0 2 6  ...
   0 2 18 0 3 4 0 3 7 0 3 19 0 4 8 0 4 20 0 5 6 0 5 9 0 5 21 0 6 7 0 ...
   6 10 0 6 22 0 7 8 0 7 11 0 7 23 0 8 12 0 8 24 0 9 10 0 9 13 0 9 25 0 ...
   10 11 0 10 14 0 10 26 0 11 12 0 11 15 0 11 27 0 12 16 0 12 28 0 13 14 0 ...
   13 29 0 14 15 0 14 30 0 15 16 0 15 31 0 16 32 0 17 18 0 17 21 0 17 33 0 ...
   18 19 0 18 22 0 18 34 0 19 20 0 19 23 0 19 35 0 20 24 0 20 36 0 21 22 0 ...
   21 25 0 21 37 0 22 23 0 22 26 0 22 38 0 23 24 0 23 27 0 23 39 0 24 28 0 ...
   24 40 0 25 26 0 25 29 0 25 41 0 26 27 0 26 30 0 26 42 0 27 28 0 27 31 0 ...
   27 43 0 28 32 0 28 44 0 29 30 0 29 45 0 30 31 0 30 46 0 31 32 0 31 47 0 ...
   32 48 0 33 34 0 33 37 0 33 49 0 34 35 0 34 38 0 34 50 0 35 36 0 35 39 0 ...
   35 51 0 36 40 0 36 52 0 37 38 0 37 41 0 37 53 0 38 39 0 38 42 0 38 54 0 ...
   39 40 0 39 43 0 39 55 0 40 44 0 40 56 0 41 42 0 41 45 0 41 57 0 42 43 0 ...
   42 46 0 42 58 0 43 44 0 43 47 0 43 59 0 44 48 0 44 60 0 45 46 0 45 61 0 ...
   46 47 0 46 62 0 47 48 0 47 63 0 48 64 0 49 50 0 49 53 0 50 51 0 50 54 0 ...
   51 52 0 51 55 0 52 56 0 53 54 0 53 57 0 54 55 0 54 58 0 55 56 0 55 59 0 ...
   56 60 0 57 58 0 57 61 0 58 59 0 58 62 0 59 60 0 59 63 0 60 64 0 61 62 0 ...
   62 63 0 63 64 0];


  elseif  comstr(Cam,'patch'); 
    out=[1 2 6 5; 1 17 18 2;1 5 21 17;2 3 7 6;
     2 18 19 3; 3 4 8 7;3 19 20 4;4 20 24 8;
     5 6 10 9; 5 9 25 21;6 7 11 10;7 8 12 11;
     8 24 28 12; 9 10 14 13;9 13 29 25;10 11 15 14;
     11 12 16 15; 12 28 32 16;13 14 30 29;14 15 31 30;
     15 16 32 31; 17 33 34 18;17 21 37 33;18 34 35 19;
     19 35 36 20; 20 36 40 24;21 25 41 37;24 40 44 28;
     25 29 45 41; 28 44 48 32;29 30 46 45;30 31 47 46;
     31 32 48 47; 33 49 50 34;33 37 53 49;34 50 51 35;
     35 51 52 36; 36 52 56 40;37 41 57 53;40 56 60 44;
     41 45 61 57; 44 60 64 48;45 46 62 61;46 47 63 62;
     47 48 64 63; 49 53 54 50;50 54 55 51;51 55 56 52;
     53 57 58 54; 54 58 59 55;55 59 60 56;57 61 62 58;
     58 62 63 59; 59 63 64 60];out(:,:)=out(:,[1 4 3 2]);
  elseif  comstr(Cam,'edge');   out = [ 1 5 9 13; 13 29 45 61;
          61 57 53 49;49 33 17 1;1 2 3 4;13:16;61:64;49:52
          4 8 12 16;16 32 48 64;64 60 56 52;52 36 20 4];
  elseif  comstr(Cam,'face');   out=h64('patch');
  elseif  comstr(Cam,'sci_face'); out = h64('patch');
  elseif  comstr(Cam,'parent'); out = 'h64'; 
  elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,5); 
   %% #TestMat
   if 1==2
   %% #TestGen
     model=femesh('testquad4 divide 4 4'); 
     model.Node(:,5:6)=model.Node(:,5:6)*2-1;
     comstr(model.Node(:,5:7),-30,struct('NoClip',2))
     r1=feutil('getPatch',model);comstr(r1(2:end,1:4),-30,struct('NoClip',1))
   end
   
   model=femesh('testh64');
   if nargin==2; model.pl=varargin{1}; end
   [m,k,mdof]=fe_mknl(model,'NoT');
   out=stack_cell(k,m);
   disp('TestMat passed');
  
  elseif  comstr(Cam,'test');  [CAM,Cam] = comstr(CAM,5);

   if comstr(Cam,'mesh');[CAM,Cam] = comstr(CAM,5);
    model=femesh(['testhexa8' CAM]); 
    out=q16p('h8Toh64',model);
   elseif nargin==2 % specified pl
     [out,out1]=femesh(strcat(['teststruct h64' Cam],varargin{1}));
   else
     [out,out1]=femesh(strcat(['teststruct h64' Cam]));
   end

 % Basic tests of the element - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(CAM,'cvs')
  out='$Revision: 1.2 $  $Date: 2020/06/11 16:29:23 $';
 else;sdtw('''%s'' unknown',CAM);  
 end

  return
end % of standard calls with one input argument

% -----------------------------------------------------------------------------
% element matrix assembly - - - - - - - - - - - - - - - - - - - - - - - - - -

error('This is only implemented as a mat_og element');
