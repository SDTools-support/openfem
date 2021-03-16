function out=t_constit(varargin);

% Test non regregression on material/property functions

if nargin==0
   
  t_constit('elastic')
  t_constit('beam')
  %% #p_beam ->sdtweb t_femk('p_beam')
  %% #shell #p_shell ->sdtweb t_plate('laminate')
  % basic p_fun ->sdtweb t_femk('p_fun')
  return
  
end
%#ok<*ASGLU,*NASGU,*STOUT,*TRYNC,*NOSEM>
[CAM,Cam]=comstr(varargin{1},1);

%% #Elastic : check building of elastic properties ---------------------------
if comstr(Cam,'elastic')
 
 model=femesh('testhexa8b');
 assignin('base','DiagConst',[]);
 % isotropic
 r1=[210e9 .3 7800];
 model.pl=[100 fe_mat('m_elastic','SI',1) r1]; % check G 
 sp_util('diag',2);
 Case=fe_mknl('init',model);
 %comstr(DiagConst.D100,-30)
 DiagConst=evalin('base','DiagConst');
 r1=blkdiag([r1(1)*(1-r1(2))/(1+r1(2))/(1-2*r1(2))*...
  [       1         r1(2)/(1-r1(2))  r1(2)/(1-r1(2))
   r1(2)/(1-r1(2))        1           r1(2)/(1-r1(2))
   r1(2)/(1-r1(2)) r1(2)/(1-r1(2))          1
  ]],...
  [r1(1)/2/(1+r1(2))*eye(3)]);
 r1=norm(DiagConst.D100-r1);
 if r1~=0; error('Something has changed, m_elastic, subtype 1'); end
%  r1=norm(DiagConst.D100-[282692307692.308,121153846153.846,121153846153.846,0,0,0;
%      121153846153.846,282692307692.308,121153846153.846,0,0,0;
%      121153846153.846,121153846153.846,282692307692.308,0,0,0;
%      0,0,0,80769230769.2308,0,0;0,0,0,0,80769230769.2308,0;
%      0,0,0,0,0,80769230769.2308;
%      ]);

 % Orthotropic (6)
  model.Elt(2,9)=102; 
  % E1(1) E2(2) E3(3) Nu23(4) Nu31(5) Nu12(6) G23(7) G31(8) G12(9) rho(10)
  r1=[6.5e3 2.5e3 4.8e2 2.5e-1 8e-2 1.8e-1 2.75e3 3.8e3 1.8e3 2.7e-9];
  model.pl=[102 fe_mat('m_elastic','SI',6) r1];
  Case=fe_mknl('init',model);
  % nu_ji=Ej/Ei nu_ij
  r1=pinv([  1/r1(1)      -r1(6)/r1(1)  -r1(5)/r1(3)   0         0       0  
      -r1(6)/r1(1)     1/r1(2)     -r1(4)/r1(2)   0         0       0
      -r1(5)/r1(3)   -r1(4)/r1(2)    1/r1(3)      0         0       0
      0                 0              0          1/r1(7)   0       0
      0                 0              0          0         1/r1(8) 0
      0                 0              0          0         0       1/r1(9)
   ]); % ^(-1); % compliance from the doc, 
  % make sure inversion is the same than in fomulaothro, pinv/inv, ...
  % sdtweb m_elastic('formulaortho')
  DiagConst=evalin('base','DiagConst');r1=norm(DiagConst.D102-r1);
 if r1~=0; error('Something has changed, m_elastic, subtype 6'); end

 % COMP12 : SVD of constitutive law gives
 %  isotropic : uniform compression, shear, ...
 %  fully orthotropic : different components
 [u,s]=eig(DiagConst.D102);[u diag(s)]
 [u,s]=eig(DiagConst.D100);[u diag(s)]
 
 % Orthotropic
  model.Elt(2,9)=103; % 
  model.pl % set mat with 103
  
 % compare isotropic to orthotropic law 
 pw0=pwd;
 cd(fileparts(fileparts(which('fe_mk'))));eval('sdtcheck path');
 cd(pw0);
 model=femesh('testhexa8b');
 K0=fe_caseg('Assemble -matdes 2 3 1 -reset -cdof -cell',model); % Isotropic 
 E=model.pl(model.pl(:,1)==100,3);
 Nu=model.pl(model.pl(:,1)==100,4);
 rho=model.pl(model.pl(:,1)==100,5);
 % E1(1) E2(2) E3(3) Nu23(4) Nu31(5) Nu12(6) G23(7) G31(8) G12(9) rho(10)
 r1=[E E E Nu Nu Nu E/(2*(1+Nu)) E/(2*(1+Nu)) E/(2*(1+Nu)) rho];
 model.pl=[100 fe_mat('m_elastic','SI',6) r1];
 K=fe_caseg('Assemble -matdes 2 3 1 -reset -cdof -cell',model); % Isotropic 
 if norm(K{1}-K0{1},'inf')/norm(K0{1},'inf')>1e-14||norm(K{2}-K0{2},'inf')/norm(K0{2},'inf')>1e-14 
  error('Something has changed, m_elastic, subtype 6 / subtype 1');
 end
 
%% #Gibson
 model=femesh('testhexa8b');
 % xxx inconsistent units SI, mm
 % E(1) Nu(2) rho(3)    a(4) b(5) theta(6) t(7) t'(8)
 
 RO.dim=[72e3 .3 2700   5e-3, 5e-3, 30, .1e-3, .2e-3];

 model.pl=m_elastic('dbval 100 gibson-punit SI  ',RO.dim, ...
    'dbval 201 Steel');
 [constit,integ,i1,dd]=p_solid('buildconstit',[100;111;3;1],model.pl, ...
    model.il,model);dd=dd.dd
 %  sdtweb m_elastic('gibson')
st=sprintf(['dbval 100 gibson -punit SI ' ...
    'E=%.15g Nu=%.15g rho=%.15g a=%g b=%g theta=%.15g t=%g tp=%g'],RO.dim);
 pl=m_elastic(st);
 if norm(pl(1,:)-model.pl(1,1:12))>1e-4;error('Mismatch');end
 
 
%% #CompShell ----------------------------------------------------------------
model=femesh('testquad4');
model.pl=m_elastic('dbval 1 Steel'); mat=model.pl;
model=feutil('setmat',model,[2 fe_mat('m_elastic','SI',5) 210e9 210e9 .285 ...
    mat(6)*[1 1 1] mat(5)]);
il = p_shell('dbval 1 kirchhoff .1 -f5','dbval 2 laminate z0=-.05 2 .1 0',...
  'dbval 3 laminate z0=-.05 1 .1 0');
[constit,integ,i1,dd1]=p_shell('buildconstit',[1;1;6;1],model.pl,il,model);
[constit,integ,i1,dd2]=p_shell('buildconstit',[2;1;6;1],model.pl,il,model);
dd2.layer(1)=1; 
st=sd('compare',dd1,dd2,struct('RelTol',1e-12));
if ~isempty(st);error('Mismatch');end

% sdtweb m_elastic('buildply')
[constit,integ,i1,dd2]=p_shell('buildconstit',[2;2;6;1],model.pl,il,model);
dd2.layer(1)=1; 
st=sd('compare',dd1,dd2,struct('RelTol',1e-12,'FieldOpt',{{'eta','skip'}}))
%if ~isempty(st);error('Mismatch');end
[constit,integ,i1,dd2]=p_shell('buildconstit',[2;3;6;1],model.pl,il,model);
dd2.layer(1)=1; 
st=sd('compare',dd1,dd2,struct('RelTol',1e-12,'FieldOpt',{{'eta','skip'}}))


%% #Beam --------------------------------------------------------------------
elseif comstr(Cam,'beam')
    
model=femesh('testbeam1');
model.pl=m_elastic('dbval 100 Steel');
Case=fe_mknl('init',model);

%sdtweb p_solid('constitbeam')
r1=[210000000000;0.285;7800;81712062256.8093;
        0;0;0;2.253333333333334e-08;1.33333333333333e-8;1.33333333333333e-8;
        0.0004;0;0;0;0];i1=find(r1); 
r1(6:7)=Case.GroupInfo{4}(6:7); % thermal not checked
if good_compare(Case.GroupInfo{4},r1)>1e-5; error('mismatch');end


%% -------------------------------------------------------------------------
else; sdtw('''%s'' not known',CAM);
end % commands

function out=good_compare(r1,r2)

if ~isequal(size(r1),size(r2));out=-1; return;end

i1=find(r1==0); if any(r2(i1));out=-2;return;end
r1(i1)=[];r2(i1)=[];out=max(abs(r1./r2-1));

