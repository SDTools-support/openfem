function [out,out1,out2]=fsc3(CAM,varargin);

%FSC3 Non linear follower force test function.
%
%	Implementation is done in elem0.m fs_matrix formulations
%       Properties are [ProId fe_mat('p_solid','SI',3) Integ 4]
%

%	Etienne Balmes, Dominique Chapelle, Mathieu Alba
%       Copyright (c) 2001-2011 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*NOSEM>

if comstr(CAM,'cvs')
 out='$Revision: 1.28 $  $Date: 2014/05/30 06:53:36 $'; return;
end

% standard calls with one input argument
if ischar(CAM)

 [CAM,Cam]=comstr(CAM,1);
 if comstr(Cam,'constants');  error('This had moved to p_solid');
 elseif  comstr(Cam,'test');   [CAM,Cam] = comstr(Cam,5);


 if comstr(Cam,'simple')||comstr(Cam,'load')
   [model,def]=femesh('test tetra4b divide 3 3 3 struct');

   elt=feutil('selelt selface & innode {x==1}',model);
   elt=feutil('egid-vol',elt);
   elt=feutil('set group1 proid 200',model.Node,elt);
   model.Elt(end+[1:size(elt,1)],1:size(elt,2))=elt;
   RO.INTEG=6;
   model.il=p_solid(model.il,'dbval 200',[200 fe_mat('p_solid','SI',3) RO.INTEG 4 0]);

   [Case,model.DOF]=fe_mknl('init',model);
   fprintf('DOFs= ');fprintf(' %i ',unique(round(rem(model.DOF,1)*100))); 
   fprintf('\n');
   if norm(double(Case.GroupInfo{2,3})-([112 200 12 3 RO.INTEG]'))
       error('Init problem');
   end
   
   Case.GroupInfo{1,end}.material='Elastic3DNL'; 
   dc=struct('def',zeros(length(model.DOF),2),'DOF',model.DOF, ...
    'lab',{{'State','RHS'}});
   dc.def(fe_c(dc.DOF,.19,'ind'),1)=1;
   [k,C1,dc]=fe_mknl('assemble',model,Case,dc,5); % non linear follower
   try; feplot(model,dc);fecom('colordataa');end

  % test linear fluid/structure coupling

   sp_util('diag',12)
   mfs1=fe_mknl('assemble not',model,Case,2);
   kfs1=fe_mknl('assemble not',model,Case,1);
   sp_util('diag',0);
   mfs=fe_mknl('assemble not',model,Case,2);
   kfs=fe_mknl('assemble not',model,Case,1);
   [mf2,kf2,mdof]=fe_mk(model,'options',[0 2]);
   
   r1=[norm(full(kfs1-kfs))/norm(full(kfs1)) ...
       norm(full(mfs1-mfs))/norm(full(mfs1))];
   fprintf('Relative error on K and M :       %g        %g \n',...
            r1);
   in1=fe_c(model.DOF,.01,'ind');in2=fe_c(model.DOF,.19,'ind');
   if normest(mfs(in2,in1)+kfs(in1,in2)')>eps; error('Obsolete version');end
   if abs(sum(sum(mfs(in2,in1),2))-1/3)>1e-10; error('Incorrect surface');end
   
   if norm(r1)>1e-13; error('inconsistent assembly');end
   if nargout>1
       out=model; dc.def(:,3)=sum(mfs,2); dc.def(:,4)=sum(kfs,2);out1=dc;
   else
    ind=[in1;in2];
    figure(1);subplot(121);eval('ii_plp(''spy'',mfs(ind,ind));')
    figure(1);subplot(122);eval('ii_plp(''spy'',kfs(ind,ind));')
   end
   
else
    disp('obsolete test, run Rivlincube (2nd pass)');out='';
    
%   % test using Rivlin reference - - - - - - - - - - - - - - - - -
%    if isempty(CAM); ElemF='tetra4b';else; ElemF=CAM;end
%    L=[.1 .1 .1];
%    addpath(fullfile(fileparts(which('fsc3.m')),'demos'));
%    eval('RivlinCube');
%    elt=feutil('selelt selface & innode {x==1 | y==1 | z==1}',model)
%    elt=feutil('egid-vol',elt);
%    elt=feutil('set group1 proid 200',model.Node,elt);
%    model.Elt(end+[1:size(elt,1)],1:size(elt,2))=elt;
%    model.il=p_solid(model.il,'dbval 200',[200 fe_mat('p_solid','SI',3) 0 4 0]);
% 
%    model=fe_case(model,'remove','xs','remove','ys','remove','zs');
%    fe_case(model,'info')
% 
%    [Case,model.DOF]=fe_mknl('init',model);
% 
%    % exact solution in input for resolution
%    def=struct('def',zeros(size(model.DOF)),'DOF',model.DOF);
%    NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
%    i1=fe_c(model.DOF,.01,'ind');
%    def.def(i1,1)=model.Node(NNode(fix(model.DOF(i1))),5)*L(1);
%    i1=fe_c(model.DOF,.02,'ind');
%    def.def(i1,1)=model.Node(NNode(fix(model.DOF(i1))),6)*L(2);
%    i1=fe_c(model.DOF,.03,'ind');
%    def.def(i1,1)=model.Node(NNode(fix(model.DOF(i1))),7)*L(3);
% 
%    % Allow for user restart during development phase : 
%    % results in base workspace
%    % pressure
%    p1=-(1+L(1))/(1+L(2))/(1+L(3))*Sigma(1,1);
%    p2=-(1+L(2))/(1+L(1))/(1+L(3))*Sigma(2,2);
%    p3=-(1+L(3))/(1+L(1))/(1+L(2))*Sigma(3,3);
%    FEnode=model.Node;
%    indn=femesh('findnode x==1');
%    indp1=fe_c(model.DOF,indn+.19,'ind');
%    indn=femesh('findnode y==1');
%    indp2=fe_c(model.DOF,indn+.19,'ind');
%    indn=femesh('findnode z==1');
%    indp3=fe_c(model.DOF,indn+.19,'ind');
%    def.def(indp1,1)=p1;
%    def.def(indp3,1)=p3;
%    def.def(indp2,1)=p2;
%    
%    % residual computation with follower pressure instead of F1  
%    def.def(:,2)=0; k=fe_mknl('assemble not',model,Case,def,5);
%    R=-Case.T(:,fe_c(Case.DOF,[.01;.02;.03],'ind'))'*(def.def(:,2));
%    fprintf('\nExact residual with follower pressure %g\n',norm(R))
% 
% %F2=Case.T(:,fe_c(Case.DOF,[.01;.02;.03],'ind'))'*(def.def(:,2));
% %iz=find(abs(F1)>1e-5);
% %[F1(iz) F2(iz) F1(iz)./F2(iz)];full(ans)
% 
%    st={'model','Case','def','L','F1','Sigma'};
%    for j1=1:length(st);assignin('base',st{j1},eval(st{j1}));end



end % tests

 % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 else; 
  if nargout==0; t3p(CAM,varargin{:});
  elseif nargout==1; out=t3p(CAM,varargin{:});
  elseif nargout==2; [out,out1]=t3p(CAM,varargin{:});
  elseif nargout==3; [out,out1,out2]=t3p(CAM,varargin{:});
  end
 end
 return
end % of standard calls with one input argument
