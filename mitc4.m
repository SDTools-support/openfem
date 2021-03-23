function [out,out1,out2]=mitc4(CAM,varargin);

%MITC4 element function for the 4-node 20-DOF quadrangular shell element
%
%	In an model description matrix a group of MITC4 elements starts with a
%	header row [Inf  abs('mitc4') 0 ...] followed by element property rows
%       ELT following the format
%	    [n1 n2 n3 n4 MatId ProId EltId]
%         with
%	   n1 ... n4  identification numbers for the element nodes
%	   MatId  material property identification number
%	   ProId  element property identification number (only used by upcom)
%	   EltId  optional element identifier
%
%     PL material property rows are detailed under m_elastic
%
%     IL Element property rows for 2-D elements follow the format
%      [ProId Type  f   d Rot h]
%     with 
%       Type = fe_mat('p_shell','SI',1)
%       f    : 3 (notel=1 standard case) 4 (notel=10 no interpolation) 
%       d    : -1 use a 5 DOF formulation, otherwise 6 DOF used
%       Rot  : rotation DOF map generation method (ichoix)
%       h    : thickness
%
%     Standard tests available with mitc4('testmat') (point,vol) 
%     Note that the normal maps are built using the fe_mknl init phase. To
%     modify your normal map by hand use
%
%      model=mitc4('testreturn');
%      Case=fe_mknl('init',model);
%      % You can change the map here
%      Case.GroupInfo{1,7}
%      k=fe_mknl('assemble',model,Case,1);
%
%
% 	See sdtweb     eltfun, elem0
%	See also help  dktp,quad4, ...


%	Marina Vidrascu, Dominique Chapelle, Etienne Balmes
%       Copyright (c) 2001-2009 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

if comstr(CAM,'cvs')
 out='$Revision: 1.32 $  $Date: 2009/05/28 16:42:00 $'; return;
end
% standard calls with one input argument
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1);
 if comstr(Cam,'integinfo') % form constitutive information

  [pe,ie,dm,db,ds]=fe_mat(2,varargin{1},varargin{2},varargin{3});

  %constit [rho eta E nu T(1:4) d]
  constit = [pe(5) 0 pe(3:4) ie([6 6 6 6]) 0];out=constit(:); 
  if length(constit)<5; 
   sdtw('Properties (%i,%i) where not found',varargin{1});
  end

  out1=[pe(1) ie(1) 24 4  0 0 0 1 ie(5) 6]'; % integ
  % IOPT in call_matrix
  %                 1  2  3 4 5 6 7 8 9
  if ie(3)==4;  out1(6+2)=10; end % NOTEL (1, 10 for debugging)
  if ie(5)==0;  out1(7+2)=1; end % ichoix (default 1) xxx needs documentation
  if ie(4)==-1; out1(8+2)=5; end % 56 DDL (defaut 6)

  % element map
  out2=[];
 % group initialization command
 elseif comstr(Cam,'groupinit') 

           % State                ElementConstants        Reindex
   out='[gstate,Case.GroupInfo{jGroup,8},elmap]=mitc4(''constants'',model.Node,model.Elt(EGroup(jGroup):EGroup(jGroup+1)-1,:),integ);';

 % Check Map existence and return integration constants 
 elseif comstr(Cam,'constants') % group initialization command

 elt=varargin{2};
 MAP=feutil('getnormal map node',varargin{1},elt); 
 NNode=sparse(MAP.ID,1,1:length(MAP.ID)); MAP.normal=MAP.normal';
 out=reshape(MAP.normal(:,full(NNode(elt(2:size(elt,1),1:4)'))), ...
   12,size(elt,1)-1);

out1=struct('const',[ 6.2200846792814624e-1 ;1.6666666666666663e-1 ;4.4658198738520435e-2 ;1.6666666666666663e-1 ;1.6666666666666663e-1 ;6.2200846792814624e-1 ;1.6666666666666663e-1 ;4.4658198738520435e-2 ;4.4658198738520435e-2 ;1.6666666666666663e-1 ;6.2200846792814624e-1 ;1.6666666666666663e-1 ;1.6666666666666663e-1 ;4.4658198738520435e-2 ;1.6666666666666663e-1 ;6.2200846792814624e-1 ;5.00e-1 ;5.00e-1 ;0.00e+0 ;0.00e+0 ;0.00e+0 ;0.00e+0 ;5.00e-1 ;5.00e-1 ;5.00e-1 ;0.00e+0 ;0.00e+0 ;5.00e-1 ;0.00e+0 ;5.00e-1 ;5.00e-1 ;0.00e+0 ;-3.9433756729740643e-1 ;-3.9433756729740643e-1 ;3.9433756729740643e-1 ;-1.0566243270259354e-1 ;1.0566243270259354e-1 ;1.0566243270259354e-1 ;-1.0566243270259354e-1 ;3.9433756729740643e-1 ;-3.9433756729740643e-1 ;-1.0566243270259354e-1 ;3.9433756729740643e-1 ;-3.9433756729740643e-1 ;1.0566243270259354e-1 ;3.9433756729740643e-1 ;-1.0566243270259354e-1 ;1.0566243270259354e-1 ;-1.0566243270259354e-1 ;-1.0566243270259354e-1 ;1.0566243270259354e-1 ;-3.9433756729740643e-1 ;3.9433756729740643e-1 ;3.9433756729740643e-1 ;-3.9433756729740643e-1 ;1.0566243270259354e-1 ;-1.0566243270259354e-1 ;-3.9433756729740643e-1 ;1.0566243270259354e-1 ;-1.0566243270259354e-1 ;3.9433756729740643e-1 ;1.0566243270259354e-1 ;-3.9433756729740643e-1 ;3.9433756729740643e-1 ;-5.00e-1 ;-2.50e-1 ;5.00e-1 ;-2.50e-1 ;0.00e+0 ;2.50e-1 ;0.00e+0 ;2.50e-1 ;0.00e+0 ;-2.50e-1 ;0.00e+0 ;-2.50e-1 ;5.00e-1 ;2.50e-1 ;-5.00e-1 ;2.50e-1 ;-2.50e-1 ;-5.00e-1 ;2.50e-1 ;0.00e+0 ;2.50e-1 ;0.00e+0 ;-2.50e-1 ;5.00e-1 ;-2.50e-1 ;0.00e+0 ;2.50e-1 ;-5.00e-1 ;2.50e-1 ;5.00e-1 ;-2.50e-1 ;0.00e+0  ]);%OfMkPol

  integ=varargin{3};
 % Elmap for 6 DOF element
 if size(integ,1)<10||(all(integ(10,:)==6)||all(integ(10,:)==7)) 
  i1=find(triu(ones(24))); r1=zeros(24);r1(i1)=1:length(i1);
  out2=r1+tril(r1',-1);
 elseif all(integ(10,:)==5) % Elmap for 5 DOF element
  i1=find(triu(ones(20))); r1=zeros(20);r1(i1)=1:length(i1);
  r1=r1+tril(r1',-1);
  out2=ones(24,24)*211; i1=[1:5 7:11 13:17 19:23];  out2(i1,i1)=r1;
 else
  error('You cannot mix 5 and 6 DOF elements in the same group');
 end

elseif comstr(Cam,'call')||comstr(Cam,'matcall')  % call for matrix assembly
    out = ['[k1,m1] = mitc4(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,elmap,Case.GroupInfo{jGroup,8},gstate(:,jElt));'];
    out1=0;
 elseif comstr(Cam,'rhscall') % call for load assembly
   out='be=mitc4(nodeE,elt(cEGI(jElt),:),pointers(:,jElt),integ,constit,elmap,Case.GroupInfo{jGroup,8},gstate(:,jElt),defe);';

  elseif  comstr(Cam,'node');  out = [1 2 3 4];
  elseif  comstr(Cam,'prop');  out = [5 6 7]; 
  elseif  comstr(Cam,'dof');
   out = [1:4];k=[1:6]'/100;out=out(ones(6,1),:)+k(:,ones(4,1));
   out=out(:);
  elseif  comstr(Cam,'line');   out = [1 2 3 4 1];
  elseif  comstr(Cam,'patch');  out = [1 2 3 4];
  elseif  comstr(Cam,'edge');  out = [1 2; 2 3; 3 4; 4 1];
  elseif  comstr(Cam,'face');   out = [1 2 3 4];
  elseif  comstr(Cam,'parent');   out = 'quad4'; 

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 elseif  comstr(Cam,'testmat'); [CAM,Cam] = comstr(CAM,9);

   model=femesh('testquad4'); model.il(1,4)=-1;
   if comstr(Cam,'6');  model.il(1,4)=0; end % 6 DOF test

   [constit,integ,elmap]=mitc4('integinfo',[100;110],model.pl,model.il);
   node=model.Node; NNode(node(:,1))=1:size(node,1); elt=model.Elt;jElt=2;

   [map,EltConst,elmap]=mitc4('Constants',model,node,elt,integ);

   integ(10)=6;[k,m]=mitc4(node(:,[5:7 1]),elt(jElt,:),[0 0 0 0 0 0 0 0 0],int32(integ),constit,elmap,EltConst,map);

 rb=[eye(6);eye(6);eye(6);eye(6)];
 i1=1:6:size(rb,1);
 rb(i1+1,4)=-node(:,7);rb(i1+2,4)=node(:,6);
 rb(i1+2,5)=-node(:,5);rb(i1+0,5)=node(:,7);
 rb(i1+0,6)=-node(:,6);rb(i1+1,6)=node(:,5);
 
 %mdof=mitc4('DOF');cf.def={rb,mdof}

 integ(10)=7;[k,m]=mitc4(node(:,[5:7 1]),elt(jElt,:),[0 0 0 0 0 0 0 0 0],int32(integ),constit,elmap,EltConst,map);

 if max(svd(rb'*k*rb))>1e-4
  r2=null(k); r2=r2*inv(r2(1:6,1:6));r2(find(abs(r2)<1e-8))=0;
  i1=find(sum(abs(rb'*k*rb))>1e-4); for j1=1:length(i1)
     disp(reshape(r2(:,i1(j1)),6,4))  % Should be rotation around X
  end
  error('Inconsistent rigid body modes');
 end

 %for j1=1:6; r2=r2-rb(:,j1)*r2(j1,:);end;r2(find(abs(r2)<1e-8))=0;
   
  out={k,m};
  disp('TestMat passed');


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  elseif  comstr(Cam,'test');   [CAM,Cam] = comstr(Cam,5);

   model=femesh('testquad4 divide 10 10');
   model.Elt(1,1:7)=[Inf abs('mitc4') 0];
   model=fe_case(model,'FixDof','encas','x==0', ...
    'dofload','z',3.03);
   [m,k,model.DOF]=fe_mk(model);
   if comstr(Cam,'return'); out=model;return
   elseif comstr(Cam,'eig') % punctual loads

    [md1,f1]=fe_eig(m,k,[5 10 1e3 11]);
    def=struct('def',md1,'DOF',model.DOF,'data',f1/2/pi);
    
   elseif comstr(Cam,'point') % punctual loads
    b=zeros(length(mdof),1);
    b(fe_c(model.DOF,[14.03 ]','ind'))=10;
    b(fe_c(model.DOF,[74.03  57.03]','ind'))=-1;
    def.DOF=mdof; def.def=k\b;
   elseif  comstr(Cam,'vol')  % volume load test
    data=struct('sel','groupall','dir',[0 0 9.81 ]);
    Case1=struct('Stack',{{'FVol','Gravity',data}});
    def = fe_load(model,Case1);
   else error('Not a valid test');return
   end

   if sp_util<5;  feplot(model.Node,model.Elt,def.def,def.DOF,2);
   else cf=feplot;cf.model=model;cf.def=def; fecom view3; end
   out=def;
  end
  return
end % of standard calls with one input argument

% element matrix assembly - - - - - - - - - - - - - - - - - - - - - - - - - -

%Written by integrules('mitc4')

node=CAM; 
elt=varargin{1}; 
point=varargin{2};
out1=[];
if any(point(5)==[0 1 2])

 % etrmit4.f call coor(1),constit(3),integ(8),constit(9)=pol,normal=coor(17),out
 % x contains [3x4 node coordinates, 4 node numbers 3x4 normals]
%     of_mk(ElemF,point,integ,constit,nodeE,gstate,
%                         defe,EltConst,InfoAtNode)

 if any(point(5)==[0 2]); % Simple mass with one integration point
   % THIS REALLY NEEDS TO BE OPTIMIZED AND FURTHER CHECKED
   m=zeros(24); 
   [bas,x]=basis(node(:,1:3));
   r1=integrules('2d',integrules('quad4',2),x);
   B=spalloc(3,24,12); 
   for jw=1:size(r1.w,1)
     B(1,1:6:24)= r1.N(jw,:);
     B(2,2:6:24)= r1.N(jw,:);
     B(3,3:6:24)= r1.N(jw,:);
     m = m + (varargin{4}(1)*mean(varargin{4}(5:7))*r1.jdet(jw)*r1.w(jw,4)*B')*B;
   end % loop on integration points
   if point(5)==0; out1=m; else out=m;return; end
 end
 point(5)=1;point(1)=300;
%integ=varargin{3};integ([6:8]+2)';integ(8:10)=[1 1 6];integ=int32(integ);point'
 k1=of_mk('mitc4',int32(point),varargin{3:4},node,[],[],varargin{6}.const,varargin{7});
 out=k1(varargin{5});

 if size(out)==24 
  r1=[1 1 1 -1 -1 1 1 1 1 -1 -1 1 1 1 1 -1 -1 1 1 1 1 -1 -1 1];
  out=diag(r1)*out*diag(r1);
 end

elseif point(5)==100 % volume load

  defe=varargin{8};defe(5,:)=0;
  be=of_mk('mitc4',int32(point),varargin{3:4},node,[],defe,varargin{6}.const,varargin{7});

 if varargin{3}(10)==5;  i1=[1:5 21 6:10 22 11:15 23 16:20 24];out=be(i1);
 else out=be; 
 end

elseif point(5)==2; out=zeros(24);
elseif point(5)~=1; sdtw('%i not a supported MITC4 matrix',point(5));
end







