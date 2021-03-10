function [out,out1,out2]=q36p(CAM,varargin)

%Q36P element function 36 node plane element (5x5 grid)
%
% 	See sdtweb     eltfun, elem0
%	See also help  t3p, ...

%       Etienne Balmes  
%       Copyright (c) 2001-2015 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

%#ok<*NOSEM,*ASGLU,*NASGU>

% standard calls with one input argument
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1);
 % Build Constit, Integ, and Elmap for later integration
 if comstr(Cam,'integinfo')
  %constit integ,elmap                 ID,pl,il
  [out,out1,out2]=p_solid('buildconstit2 36',varargin{:});

 elseif comstr(Cam,'matcall'); [out,out1]=elem0(CAM,varargin{:});
 elseif comstr(Cam,'dofcall');out=elem0('dofcall'); % variable field elements
 elseif comstr(Cam,'call')
   out='[k1,m1]=q36p(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,elmap);';
 elseif comstr(Cam,'rhscall') % call for load assembly
   out='rhs_of';

 elseif comstr(Cam,'groupinit');out=elem0(CAM,'q36p');
 % Here one defines the constants needed for element integration
 elseif comstr(Cam,'constants');

  if nargin<3; p_solid('constsolid','q36p',[],[]);return;
  elseif varargin{2}(size(varargin{2},1),1)==-9999; % old of_mk_sub.c elements
    i2=16;i1=find(triu(ones(i2,i2)));r1=zeros(i2,i2);r1(i1)=1:length(i1);
    out=struct('ElMap',r1+tril(r1',-1));   out2=[];
    out1=varargin{1};out1(4,:)=2;% Tell of_mk('MatrixIntegration') this is 2d
  else; 
   [out,i2,out2]=p_solid('constsolid','q36p',varargin{2:end});
   out1=varargin{1};out1(4,:)=i2; % Tell MatrixIntegration this is nd 
  end

elseif  comstr(Cam,'rule') %% #Rule integration rule
 out=feval(t_spectral('@N_Nr_Ns'),6);

  elseif  comstr(Cam,'node');  out = 1:36;
  elseif  comstr(Cam,'prop');  out = (37:38); 
  elseif  comstr(Cam,'dof')   
    out = (1:25);k=[1;2]/100;out=out(ones(2,1),:)+k(:,ones(8,1));
    out=out(:);
  elseif  comstr(Cam,'line');  
    out=q36p('patch'); out=[out(:,1:2);out(:,2:3);out(:,3:4);out(:,[4 1])];
    out=unique(sort(out,2),'rows');
  elseif  comstr(Cam,'patch'); 
    out=[1 7 8 2];out=[out;out+1;out+2;out+3;out+4];
    out=[out;out+6;out+12;out+18;out+24];
  elseif  comstr(Cam,'edge');   out = [ 1 5 6 7 2;2 12 17 22 3;
      3 25 24 23 4 ; 4 18 13 8 1 ];
  elseif  comstr(Cam,'face');   out = [1 2 3 4];
  elseif  comstr(Cam,'sci_face'); out = [1:4];
  elseif  comstr(Cam,'parent'); out = 'q36p'; 
  elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,5); 
   %% #TestMat
   if 1==2
   %% #TestGen
     model=femesh('testquad4 divide 4 4'); 
     model.Node(:,5:6)=model.Node(:,5:6)*2-1;
     comstr(model.Node(:,5:7),-30,struct('NoClip',2))
     r1=feutil('getPatch',model);comstr(r1(2:end,1:4),-30,struct('NoClip',1))
   end
   
   model=femesh('testq36p');
   if nargin==2; model.pl=varargin{1}; end
   if ~isempty(strfind(Cam,'_0'));     model.il(:,3)=0;
   elseif ~isempty(strfind(Cam,'_1')); model.il(:,3)=1;
   elseif ~isempty(strfind(Cam,'_2')); model.il(:,3)=2;
   end
   [m,k,mdof]=fe_mknl(model,'NoT');
   out=stack_cell(k,m);
   disp('TestMat passed');
  
  elseif  comstr(Cam,'test');  [CAM,Cam] = comstr(CAM,5);

   if comstr(Cam,'mesh');[CAM,Cam] = comstr(CAM,5);
    model=femesh(['testq4p' CAM]); 
    out=q16p('q4Toq36p',model);
   elseif nargin==2 % specified pl
     [out,out1]=femesh(strcat(['teststruct q36p' Cam],varargin{1}));
   else
     [out,out1]=femesh(strcat(['teststruct q36p' Cam]));
   end

 % Basic tests of the element - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(CAM,'cvs')
  out='$Revision: 1.2 $  $Date: 2015/03/02 18:49:27 $';
 else sdtw('''%s'' unknown',CAM);  
 end

  return
end % of standard calls with one input argument

% -----------------------------------------------------------------------------
% element matrix assembly - - - - - - - - - - - - - - - - - - - - - - - - - -

error('This is only implemented as a mat_og element');
