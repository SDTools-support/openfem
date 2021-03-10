function [out,out1,out2]=quadc(CAM,varargin);

%QUADC composite plate
%
% 	See sdtweb     eltfun, elem0
%	See also help  t3p, ...

%	Jean-Michel Leclere, Amine Hassim  
%       Copyright (c) 2001-2020 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

if comstr(CAM,'cvs')
 out='$Revision: 1.16 $  $Date: 2020/02/26 08:47:29 $'; return;
end
% standard calls with one input argument
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1);
 % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 % Matrix integration
 if comstr(Cam,'matrix')

% -----------------------------------------------------------------------------
% element matrix assembly - - - - - - - - - - - - - - - - - - - - - - - - - -

[jElt,DofPos,NodePos,node, ...
      pointers,integ,constit,gstate, ...
      elmap,InfoAtNode,EltConst,def]=deal(varargin{:});

point=pointers(:,jElt);
nodeE=node(NodePos(:,jElt),:);
 out=[]; out1=[];

 i1 = 1:8; x=nodeE(:,5:7);x=x-ones(8,1)*x(1,1:3);
 bas = basis(x(2,:),x(3,:));  x=x*bas;

% stiffness matrix assembly - - - - - - - - - - - - - - - - - - - - - - - -
 if any([0 1]==point(5));
  %rules=varargin{end}; of_mk('buildndn2d',rules,CAM)

if 1==2;
rule=[1 9]; rule1=rule;
   r2=[1 -2 1 rule1; 2 3 2 rule1;3 3 1 rule1;3 -2 2 rule1;% e_xx =N,x u,e_yy,e_xy
     4 -2 5 rule1;5 3 4 rule1;6 2 4 rule1;6 -3 5 rule1; % curvatures
     7 2 3 rule;7 1 5 rule;8 3 3 rule;8 -1 4 rule];  % shear constraint
 EltConst.StrainDefinition{1}=r2(2,:);
 %EltConst.StrainDefinition{1}=z.StrainDefinition{1}(5:8,:);
 EltConst=integrules('matrixrule',EltConst);
end
   point(5)=1;
   ke=of_mk('matrixintegration',1,int32([1:8]'),nodeE, ...
       int32(point),integ,constit,gstate, ...
       elmap,InfoAtNode,EltConst,def);
   k=ke(elmap); out=of_mk('xkx_trans',bas,k);
if 1==2
   mdof=evalin('caller','model.DOF(DofPos(:,jElt)+1)');
   [i1,i2]=unique(DofPos(:,jElt)+1); [i3,i4]=unique(NodePos(:,jElt));
   rb=feutilb('geomrb',nodeE,[0 0 0],mdof(i2));
   sum(abs(out*rb.def))
end

 end;
 % mass matrix assembly - - - - - - - - - - - - - - - - - - - - - - - -
 if any([0 2]==point(5));
   point(5)=2;
   ke=of_mk('matrixintegration',1,int32([1:8]'),nodeE, ...
       int32(point),integ,constit,gstate, ...
       elmap,InfoAtNode,EltConst,def);
   if min(diag(ke))<0;error('Not positive definite');end
   ke=ke(elmap); k=of_mk('xkx_trans',bas,ke);
   if pointers(5)==0; out1=k;else;  out=k;end

 end

 % Build Constit, Integ, ElMap
 elseif comstr(Cam,'integinfo')

  [out,out1,out2]=quad4('integinfo',varargin{:});
  out2=elem0('elmapmat_og',[8 6]);out1(3:4)=[48;8];

 % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 % Matrix assembly in fe_mknl
 elseif comstr(Cam,'matcall')
 out=quadc('call');  out1=0; % Call, SymFlag

 % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'call')

  out='[k1,m1]=quadc(''matrix'',jElt,DofPos,NodePos,Case.Node,pointers,integ,constit,gstate,elmap,InfoAtNode,EltConst,def);';

 % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'rhscall') % call for load assembly
   error(1);
 % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 % Callback to define element constants during the fe_mknl init phase
 elseif comstr(Cam,'groupinit');

          % ElementConstants   
   out='[Case.GroupInfo{jGroup,8},pointers]=quadc(''Constants'',pointers,integ,constit);';

 % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'straintopo');
  
  switch varargin{1}
  case 1;  % StiffNess

  case 2;  % mass
   out=[1 1 1 5 9;2 1 2 5 9; 3 1 3 5 9];
   out1=eye(3);

  otherwise; error('Not a known matrix');
  end
 

 % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 % Here one defines the constants needed for element integration
 elseif comstr(Cam,'constants');   out = '';

  opt=integrules('quadb',3);
  if nargin>3; out=p_shell('constshell',opt,varargin{2:3});
  else; p_shell('constshell',opt,[],[]);return;
  end
  out1=varargin{1};out1(4,:)=23; % Tell of_mk('MatrixIntegration') this is 2d 

 % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 % Basic matrix test of the element - - - - - - - - - - - - - - - - - -
 elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,5);

   model=femesh('testquadb');
   model.il(find(model.il(:,1)==110),3:4)=[2 -1];
   [mb,kb,mdofb]=fe_mknl(model,'not');

   model.Elt(1,1:8)=[Inf abs('quadb') 0 0];[mb,kb,mdofb]=fe_mknl(model,'not');
   model.Elt(1,1:8)=[Inf abs('quadc') 0 0];[mc,kc,mdofc]=fe_mknl(model,'not');
   %constit(9:end)=0;constit([9 16 23])=1;
   
   
   kc=full(kc);kc(find(abs(kc)<sqrt(eps)))=0; 
   i1=fe_c(mdofb,mdofc,'ind');kc=kc(i1,i1);mc=mc(i1,i1);

   if max([norm(kb-kc,'inf')/norm(kb,'inf') norm(mb-mc,'inf')/norm(mc,'inf')])>1e-15
    i1=find(kc);r1=kc;r1(i1)=r1(i1)./kb(i1);r1(find(abs(r1-1)<sqrt(eps)))=0;

    error('quadc and quadb inconsistent');
   end

   out=stack_cell(k,m);
   disp('TestMat passed');

  elseif  comstr(Cam,'test');   [CAM,Cam] = comstr(Cam,5);

   if nargin==2 % specified pl
     [out,out1]=femesh(strcat('teststruct quadc',Cam),varargin{1});
   else
     [out,out1]=femesh(strcat(['teststruct quadc' Cam]));
   end
 else; % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if nargout==0; quadb(CAM,varargin{:});
  elseif nargout==1; out=quadb(CAM,varargin{:});
  elseif nargout==2; [out,out1]=quadb(CAM,varargin{:});
  elseif nargout==3; [out,out1,out2]=quadb(CAM,varargin{:});
  end
 end % commands
end % of standard calls with one input argument

