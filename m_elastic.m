function [out,out1,out2]=m_elastic(varargin)

%M_ELASTIC Elastic material function
%
%       Syntax : mat= m_elastic('default') 
%                mat= m_elastic('database name') 
%                mat= m_elastic('database -therm name') 
%                pl = m_elastic('dbval MatId name');
%                pl = m_elastic('dbval -unit MM MatId name');
%                pl = m_elastic('dbval -punit MM MatId name');
%                pl = m_elastic('dbval -therm MatId name');
%
%  Material subtypes supported by m_elastic are
%
%	1 : standard isotropic materials ([MatID typ E nu rho G Eta Alpha	T0])
%	2 : acoustic fluid ([MatId typ rho C eta])
%	3 : 3-D anisotropic materials
%	4 : 2-D anisotropic materials
%	5 : Orthotropic material for shell
%	6 : Orthotropic material 
%
%  For more details use  m_elastic('propertyunittypecell',3)
%
%       See sdtweb      m_elastic, pl, fem
%       See also help   fe_mat, p_shell, p_beam


%	Etienne Balmes, Jean-Michel Leclere, Corine Florens
%       Copyright (c) 2001-2020 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*NOSEM>

if nargin<1; help m_elastic;return; end
if ischar(varargin{1}); [CAM,Cam]=comstr(varargin{1},1); pl=[]; carg=2;
else; pl=varargin{1};[CAM,Cam]=comstr(varargin{2},1); carg=3;
end

%% #Info -----------------------------------------------------------------------
if comstr(Cam,'info')

 matgui('info','m_elastic');%r1=m_elastic('database');disp(char({r1.name}))

%% #Default --------------------------------------------------------------------
elseif comstr(Cam,'default')
    
 model=[];if carg<=nargin;model=varargin{carg};carg=carg+1;end
 if isempty(model); out=m_elastic('database'); out=out(1);
 else;              out=fe_mat('defaultpl',model);
 end
 
%% #DBVal-----------------------------------------------------------------------
elseif comstr(Cam,'dbval') 
 RO=struct;
 while 1==1
  out1={}; % See if unit specified
  [CAM,Cam,RO.punit]=comstr('-punit',[-25 4 1],CAM,Cam);
  [CAM,Cam,RO.unit]=comstr('-unit',[-25 4 1],CAM,Cam);
  [CAM,Cam,RO.therm]=comstr('-therm',[-25 3],CAM,Cam);
  if ischar(CAM); [i1,CAM,Cam]=comstr(CAM,'dbval','%i');else; i1=[];end
  if isempty(CAM)&&carg<=nargin; st=varargin{carg};carg=carg+1; 
  else; st=CAM;
  end
  if isempty(st)
  elseif ischar(st)||iscell(st)||isstruct(st); % Call database to build constit
   st1='database'; if RO.therm; st1=[st1 '-therm']; end
   if isstruct(st);st1=[st1 'structguess'];
    if ~isfield(st,'unit');st.unit=RO.unit;end
    if ~isfield(st,'punit');st.punit=RO.punit;end
   end
   [mat,st1,i2]=m_elastic(st1,st,RO,varargin{carg:end});
   carg=carg+i2-4;
   % mat=m_elastic('database',st,RO); 
  elseif isnumeric(st)
    [typ,st1,i4]=fe_mat('typem',st(2));
    mat=struct('pl',st,'name',sprintf('%i',st(1)),'type',typ,'unit',st1);
    if ~isempty(RO.punit);
     error('For numeric pl, it should be in proper unit and punit is not valid');
    end
    if ~isempty(RO.unit)
     mat.pl=fe_mat(sprintf('convert %s%s',mat.unit,RO.unit),mat.pl);
     mat.unit=RO.unit;
    end
  end
  if (length(i1)==1); mat.pl(1)=i1;end
  r1=mat.pl; 
  if ~isempty(pl) ; i2=find(pl(:,1)==r1(1)); else; i2=[];end
  if isempty(i2)  ; i2=size(pl,1)+1;end
  pl(i2,1:length(r1))=r1;  %#ok<AGROW>
  out1(end+1,1:3)={'mat',sprintf('%i_%s',mat.pl(1),mat.name),mat};%#ok<AGROW> 
  if carg>nargin; break;end
  CAM=varargin{carg};carg=carg+1;if ischar(CAM);[CAM,Cam]=comstr(CAM,1);end
 end
 out=pl;

%% #DataBase -------------------------------------------------------------------
elseif comstr(Cam,'database'); [CAM,Cam]=comstr(CAM,9);
  RO=struct; % Later RO will be given later
  [CAM,Cam,RO.therm]=comstr('-therm',[-25 3],CAM,Cam);
  st=CAM;%[st,unu]=comstr('-therm',[-25 3],st,lower(st));
  if ~isempty(st)||carg>nargin
  elseif ischar(varargin{carg}); 
      st=varargin{carg}; carg=carg+1;
  end
  [MatId,i2,i3,i4]=sscanf(st,'%i',1); if i2~=1; MatId=1;end
  st=comstr(st,i4);

  out.pl=[MatId fe_mat('type','m_elastic','SI',1) ... 
          210e9 .285 7800 210e9/2/(1.285) 0 13e-6 20];
  out.name='Steel';
  out.type='m_elastic';
  out.unit='SI';

  out(2).pl=[MatId fe_mat('type','m_elastic','SI',1) ...
   72e9 .3 2700 72e9/2/(1.3) 0 23.1e-6 20];
  out(2).name='Aluminum';
  out(2).type='m_elastic';
  out(2).unit='SI';

  out(3).pl=[MatId fe_mat('type','m_elastic','SI',1) 2.4e9 .3  1.8e3 2.4e9/2/1.3];
  out(3).name='Nylon';
  out(3).type='m_elastic';
  out(3).unit='SI';

  out(end+1).pl=[MatId fe_mat('type','m_elastic','SI',1) ...
    120e9 .33 4700 120e9/2/(1.33) 0 7e-6 20];
  out(end).name='Titanium';
  out(end).type='m_elastic';
  out(end).unit='SI';

  out(end+1).pl=[MatId fe_mat('type','m_elastic','SI',2) 1.225 330 0];
  out(end).name='Air';
  out(end).type='m_elastic';
  out(end).unit='SI';

  % Engineering toolbox : bulk modulus of Sea water 2.34 GPa, water 2.19 GPa
  %  K=rho C^2
  out(end+1).pl=[MatId fe_mat('type','m_elastic','SI',2) 1000 1500 0];
  out(end).name='Water';
  out(end).type='m_elastic';
  out(end).unit='SI';

  out(end+1).pl=[MatId fe_mat('type','m_elastic','SI',1) ...
    4e9 .4 1200 4e9/2/(1.4) 0 1e-4 20];
  out(end).name='Epoxy_Resin';
  out(end).type='m_elastic';
  out(end).unit='SI';

  out(end+1).pl=[MatId fe_mat('type','m_elastic','SI',1) ...
    3e9 .4 1200 3e9/2/(1.4) 0 8.5e-5 20];
  out(end).name='Polyamide_Resin';
  out(end).type='m_elastic';
  out(end).unit='SI';

  out(end+1).pl=[MatId fe_mat('type','m_elastic','SI',1) ...
          4.e10 .17 2200 4.e10/2/(1.17) .01 1.2e-5 20];
  out(end).name='Concrete';
  out(end).type='m_elastic';
  out(end).unit='SI';

  out(end+1).pl=[MatId fe_mat('type','m_elastic','SI',1) ...
          15e9 .3 1380 15e9/2/(0.3)]; 
  % found in "techniques de l'ingenieur", grande incertitude pour Nu
  out(end).name='Nomex';
  out(end).type='m_elastic';
  out(end).unit='SI';
  
  out(end+1).pl=[MatId fe_mat('type','m_elastic','SI',1) ...
   1.5e9 .42 900 1.5e9/2/(0.42) 0 8e-5 20];
  out(end).name='Polypropylene';
  out(end).type='m_elastic';
  out(end).unit='SI';

  out(end+1).pl=[MatId fe_mat('type','m_elastic','SI',1) 1 0.0 1 0.5]; 
  out(end).name='Strain';
  out(end).type='m_elastic';
  out(end).unit='SI';
  
  out(end+1).pl=[MatId fe_mat('type','m_elastic','SI',6) ...
          9e9 10e9 3e9 0.125 0.25 0.2 4e9 2e9 1.5e9 3000];
  out(end).name='Orthotropic_Example';
  out(end).type='m_elastic';
  out(end).unit='SI';

  % Silicon wafer properties, Hopcroft, from 
  % JOURNAL OF MICROELECTROMECHANICAL SYSTEMS, VOL. 19, NO. 2, APRIL 2010
  out(end+1).pl=[MatId fe_mat('type','m_elastic','SI',1) 130e9 .28 2330 79.6e9];
  out(end).name='SiliconIso';
  out(end).type='m_elastic';
  out(end).unit='SI';

  out(end+1).pl=[MatId fe_mat('type','m_elastic','SI',6) 169e9 169e9 130e9 ...
      .36 .28 .064 79.6e9 79.6e9 50.9e9 2330];
  out(end).name='SiliconOrtho';
  out(end).type='m_elastic';
  out(end).unit='SI';
  
  out(end+1).pl=[MatId fe_mat('type','m_elastic','SI',5) ...
          1 1 0  1 1  1  1];
  out(end).name='StrainShell';
  out(end).type='m_elastic';
  out(end).unit='SI';
  
  if ~RO.therm % cleanup thermal terms and trim pl lines off zeros at end
   % detection is based on the convert identifier 7 or 8
   for j1=1:size(out,2)
    [unu,unu1,styp]=fe_mat('typem',out(j1).pl(2));
    r1=m_elastic('propertyUnitTypeCell',styp);
    out(j1).pl(ismember(cell2mat(r1(:,2)),[7 8]))=0; 
    out(j1).pl=out(j1).pl(1:find(out(j1).pl,1,'last')); 
   end
  end
  
  i1=find(strcmpi(st,{out.name}));
  out1='Elastic';
  
 if isempty(i1) && isempty(st); return; end
 if carg<=nargin&&isstruct(varargin{carg});RO=varargin{carg};carg=carg+1;
 else; RO=struct;
 end

  % match a name 
if ~isempty(i1); out=out(i1);
  if isfield(RO,'punit')&&~isempty(RO.punit);
   error('-punit is only meaningful for input parameters, not fixed entries');
  end
elseif comstr(st,'lamina'); st=comstr(st,7);
%% #DatabaseLamina From Jones p.91 Estimate composite lamina properties - - - -
 if ~isempty(st)&&isempty(intersect(st,'abcdfhi'))
   opt=comstr(st,-1);
   % Dbval 100 lamina VolFrac Ef(2) nu_f rho_f G_m E_m(6) nu_m Rho_m G_m
   r1={'vf','Ef','nuf','rhof','Gf','Em','num','rhom','Gm'};
   r1(2,:)=num2cell(opt); for j1=1:size(r1,2);RO.(r1{1,j1})=r1{2,j1};end
   RO.rhoc=[];
 else
  [RO,st,CAM]=cingui('paramedit -DoClean',[ ...
   'vf(#%g#"fiber volume fraction") mf(#%g#"mass fraction")'...
   'rhof(#%g#"Fiber density") Ef(#%g#"Fiber modulus longitudinal") '...
   'nuf(#%g#"Fiber nu") Eft(#%g#"Fiber modulus transverse") '...
   'Gf(#%g#"Fiber shear modulus") '...
   'rhom(#%g#"Matrix density") Em(#%g#"Matrix modulus") '...
   'num(#%g#"Matrix nu") Gm(#%g#"Matrix shear modulus") ' ...
   'rhoc(#%g#"Composite density")' ...
   ],{RO,st});
 end
 %RO.Ef=640e9; RO.Eft=12e9; RO.Gf=20e9; RO.nuf=.3; RO.rhof=2120;
 %RO.Em=3e9; RO.num=.3; RO.rhom=1200; RO.Gm=2e9;
 % 
 if isfield(RO,'mf')&&~isempty(RO.mf)&&isempty(RO.vf) % Mass fraction
  % From mass fraction build volume fraction
  RO.vf=RO.rhom*RO.mf/(RO.rhof*(1-RO.mf)+RO.rhom*RO.mf);
  % if ~isempty(RO.rhoc); RO.vf=RO.rhoc*RO.mf/RO.rhof; end
 end
 if ~isempty(RO.rhoc)
    %vf*rhof+(1-vf)*rhom=rhoc ; vf(rhof-rhom)+rhom=rhoc
    RO.rhoc=RO.vf * RO.rhof + (1-RO.vf) * RO.rhom;
 end
 RO=punitCheck(RO);
 v=[RO.vf 1-RO.vf];%v_f v_m (volume fraction)
 if isempty(RO.Gf)||RO.Gf==0; RO.Gf=RO.Ef/2/(1+RO.nuf);end
 if isempty(RO.Gm)||RO.Gm==0; RO.Gm=RO.Em/2/(1+RO.num);end
 % v_f E_f(2) nu_f(3) rho_f(4) G_f(5) E_m(6) nu_m(7) rho_m(8) G_m(9)
 out=struct('pl',[MatId fe_mat('type','m_elastic',RO.punit,5) ...
   RO.Ef*v(1)+RO.Em*v(2) ... % E_1=E_f V_f + E_m V_m (3.7)
   RO.Ef*RO.Em/(v(2)*RO.Ef+v(1)*RO.Em) ... % E_2=E_fE_m/V_mE_f+V_fE_m 3.13
   RO.nuf*v(1)+RO.num*v(2) ... % \nu_{12}=V_m nu_m+V_f nu_f
   RO.Gm*RO.Gf/(v(2)*RO.Gf+v(1)*RO.Gm) ...% G_12=G_fG_m/V_mG_f+V_fG_m 3.27
   0 0 RO.rhof*v(1)+RO.rhom*v(2)  ... % rho
   ],'unit',RO.punit,'type','m_elastic','name',sprintf('lamina_%i',MatId));

elseif comstr(st,'structguess')
%% #DatabaseStructGuess - - - - - - - - 
  RO=punitCheck(RO);
  if length(intersect(fieldnames(RO),{'E11','E22','E33'}))==3
    r2=[fieldnames(RO) struct2cell(RO)];
    r2(:,1)=strrep(r2(:,1),'1','x');r2(:,1)=strrep(r2(:,1),'2','y');
    r2(:,1)=strrep(r2(:,1),'3','z');
    r2(:,1)=strrep(r2(:,1),'nu','PR');
    r2=r2';RO=struct(r2{:});
  end
  if length(intersect(fieldnames(RO),{'Exx','Eyy','Ezz'}))==3
   % a=[fieldnames(RO) struct2cell(RO)]';
   % [i1,i2]=ismember({'Exx','Eyy','Ezz', ...
   %     'Pryz','PRzx','Prxy','Gyz','Gxz','Gxy','rho'},a(1,:));
   if ~isfield(RO,'PRzx'); RO.PRzx=RO.Ezz/RO.Exx*RO.PRxz;end%Nu31=E3/E1*Nu13
  out=struct('pl',[MatId  fe_mat('m_elastic',RO.punit,6)  ...
      RO.Exx RO.Eyy RO.Ezz RO.PRyz RO.PRzx RO.PRxy RO.Gyz RO.Gxz RO.Gxy RO.rho], ...
      'unit',RO.punit,'type','m_elastic', ...
      'name','guess');
  if isfield(RO,'name');out.name=RO.name;end
  carg=carg+1;
    
  else; error('Not implemented');
  end
  
elseif comstr(st,'iso')
 RO=punitCheck(RO);
 out=struct('name','iso','pl',[1 fe_mat('m_elastic',RO.punit,1) comstr(comstr(st,4),-1)]);
elseif comstr(st,'aniso')
 RO=punitCheck(RO);
 if comstr(st,'aniso3')
  out=struct('name','iso','pl',[1 fe_mat('m_elastic',RO.punit,3) comstr(comstr(st,7),-1)]);
 elseif comstr(st,'aniso2')
  out=struct('name','iso','pl',[1 fe_mat('m_elastic',RO.punit,4) comstr(comstr(st,7),-1)]);
 end
elseif comstr(st,'acous')
 RO=punitCheck(RO);
 out=struct('name','iso','pl',[1 fe_mat('m_elastic',RO.punit,2) comstr(comstr(st,6),-1)]);
elseif comstr(st,'shellortho')
 RO=punitCheck(RO);
 out=struct('name','iso','pl',[1 fe_mat('m_elastic',RO.punit,5) comstr(comstr(st,11),-1)]);
  
else % assume values given
  try;
    [out,carg]=formula_homo(st,out,varargin,carg); 
  catch; warning('Formula %s failed',st);out=[];
  end
  if ~isfield(out,'pl')
    error('''%s'' not a supported material',st);
  end
end
if isfield(RO,'unit')&&~isempty(RO.unit);
 if ~isempty(RO.punit); stc=['convert' RO.punit]; else; stc='convert'; end
 out.pl=fe_mat([stc RO.unit],out.pl);out.unit=RO.unit;
end 
out2=carg;

%% #BuildConstit : ->sdtweb p_solid('buildConstit') --------------------------
elseif comstr(Cam,'buildconstit');[out,out1,out2]=p_solid(varargin{:});
%% #BuildPly -------------------------------------------------------------------
%[S,RhoH]=p_shell('buildply',mat,pl,il,cz,rhoh);
elseif comstr(Cam,'buildply')

mat=varargin{carg};carg=carg+1;
pl=varargin{carg};carg=carg+1;
il=varargin{carg};carg=carg+1;
cz=varargin{carg};carg=carg+1;
rhoh=varargin{carg};carg=carg+1;

[st,unit,typ]=fe_mat('type',mat(2));
S=zeros(5); % xx,yy,zz,xy,yz,xz sdtweb p_solid('conv');
if diff(cz)<0; error('z_i+1-z_i must be positive for plies');end
   if ~strcmp(st,'m_elastic')
       if exist(st,'file')
       [S,rhoh]=feval(st,'buildply',mat,pl,il0);
       else; error('Ply not implemented for %s(SubType=%i)',st,typ);
       end       
   elseif typ==5  
     % Orthotropic material for shell. See    help m_elastic
     % ply compliance matrix
     %S=[1/E1(i) -Nu12(i)/E1(i) 0 0 0;-Nu12(i)/E1(i) 1/E2(i) 0 0 0;
     %    0 0  1/G12(i) 0 0;0 0 0 1/G23(i) 0;0 0 0 0 1/G13(i)];
     %S=eye(5);S(1:2,1:2)=1;sprintf('%i ',find(S))
     S([1 2 6 7 13])= [1/mat(3) -mat(5)/mat(3) -mat(5)/mat(3) 1/mat(4) 1/mat(6)];
     if mat(6)==0; S(13)=0; end
     % The shear coefficient should be applied at the shell level
     if mat(7); S(19)=1/mat(7); 
     elseif mat(6); S(19)=1/mat(6); fprintf('Assuming G1z=G12\n');
     end; 
     if mat(8); S(25)=1/mat(8); 
     elseif mat(6); S(25)=1/mat(6); fprintf('Assuming G2z=G12\n');
     end
     rhoh=rhoh+diff(cz)*mat(9);
   elseif typ==1 % isotropic material
     % defaults for elastic istropic material
     if strcmp(st,'m_elastic')&&typ==1 % Standard isotropic
      if length(mat)<6||mat(6)==0; mat(6)=mat(3)/2/(1+mat(4));end % G
      if size(mat,2)<7;mat(7)=0;end
     end
     %S=[1/E -Nu/E 0 0 0;-Nu/E1 1/E 0 0 0;
     %    0 0  1/G 0 0;0 0 0 1/G 0;0 0 0 0 1/G];
     % Si=E*[1 nu;nu 1]/(1-nu^2)
     if ~any(mat([3 6])); sdtw('_nb','no stiffness in mat %i',mat(1)); 
     else
     S([1 2 6 7 13])= [1/mat(3) -mat(4)/mat(3) -mat(4)/mat(3) 1/mat(3) ...
      1/mat(6)];S(19)=1/mat(6);S(25)=1/mat(6); 
     end
     rhoh=rhoh+diff(cz)*mat(5);
   elseif strcmp(st,'m_elastic')&&typ==6 
    % orthotropic material, see SDT manual classical lamination theory
    [r1,r1,r1,dd]=p_solid('buildconstit 3 1',[mat(1);1;3;1],mat,...
        p_solid('dbval 1 d3 -3'));  dd=dd.dd;h=diff(cz);
    for j1=1:2;for j2=1:2;S(j1,j2)=dd(j1,j2)-dd(j1,3)*dd(j2,3)/dd(3,3);end;end
    S(13)=dd(36);  S(19)=dd(22); S(25)=dd(29); % xy, yz, xz 
    S=inv(S); 
    rhoh=rhoh+h*mat(12);
   elseif strcmp(st,'m_elastic')&&typ==4
    fprintf('Bending stiffness ignored for m_elastic.4\n');
    S(1:3,1:3)=inv(mat([3 4 6;4 5 7;6 7 8])); 
    h=diff(cz);rhoh=rhoh+h*mat(9);
   else; error('Ply not implemented for %s(SubType=%i)',st,typ);
   end       
       
   out=S; out1=rhoh;out2=[];

%% #PropertyUnitType -----------------------------------------------------------
elseif comstr(Cam,'propertyunittype')

if nargin==1;out=[1:6]; return;end
 %out=PropertyUnitType_elastic(varargin{2});
MaterialSubType=varargin{carg};
% the indices match rows in
% fe_mat('convertSIGM');[num2cell(1:size(ans,1))' ans(:,1:2)]
out1={};

switch MaterialSubType % fe_mat('unitlabel','SI')
case 1 % isotropic   [MatId typ E nu rho G eta alpha T0]
 st=...
 {'MatId'    0  'sdtweb(''m_elastic'')';
  'Type'     0  '';
  'E'        1  'Youngs Modulus';
  'Nu'       0  'Poisson''s ratio';
  'Rho'      3  'Density';
  'G'        1  'Shear Modulus';
  'Eta'      0  'Loss factor';
  'Alpha'    8  'Thermal expansion coef';
  'T0'       7  'Reference temperature'};

case 2 % acoustic fluid [MatId typ rho C eta]
 st=...
 {'MatId'    0  'sdtweb(''m_elastic'')';
  'Type'     0  '';
  'Rho'      3  'Density';
  'C'        5  'Velocity';
  'Eta'      0  'Loss factor'};

case 3 % 3-D anisotropic solid [MatId typ Gij ... rho eta A1... A6 T0 eta]
  % sdtweb p_solid('elasaniso3')
 st=...
 {'MatId' 0  'sdtweb(''m_elastic'')';
  'Type'  0  '';
  'G11'   1  '';
  'G12'   1  '';
  'G22'   1  '';
  'G13'   1  '';
  'G23'   1  '';
  'G33'   1  '';
  'G14'   1  '';
  'G24'   1  '';
  'G34'   1  '';
  'G44'   1  '';
  'G15'   1  '';
  'G25'   1  '';
  'G35'   1  '';
  'G45'   1  '';
  'G55'   1  '';
  'G16'   1  '';
  'G26'   1  '';
  'G36'   1  '';
  'G46'   1  '';
  'G56'   1  '';
  'G66'   1  '';
  'Rho'   3  'Density';
  'Eta'   0  'Loss factor'
  'A1'    8  ''; % 26
  'A2'    8  '';
  'A3'    8  '';
  'A4'    8  '';
  'A5'    8  '';
  'A6'    8  '';
  'T0'    7  'Reference temperature';
  }; 

case 4 % 2-D anisotropic solid 
       % [MatId typ E11 E12 E22 E13 E23 E33 rho eta a1 a2 a3 T0]
 st=...
 {'MatId' 0 'sdtweb(''m_elastic'')';
  'Type'  0 '';
  'E11'   1 '';
  'E12'   1 '';
  'E22'   1 '';
  'E13'   1 '';
  'E23'   1 '';
  'E33'   1 '';
  'Rho' 3 'Density';
  'Eta' 0 'Loss factor';
  'A1'  8 'Thermal expansion coef';
  'A2'  8 'Thermal expansion coef';
  'A3'  8 'Thermal expansion coef'
  'T0'       7  'Reference temperature'};

case 5 % Orthotropic material for shell
 %[MatId type E1 E2 nu12 G12 G1z G2z Rho A1 A2 TREF Xt Xc Yt Yc S Ge F12 STRN]
 %                                                  xxx not valid xxx
 st=...
 {'MatId' 0  'sdtweb(''m_elastic'')';
  'Type'  0  '';
  'E11'   1  '';
  'E12'   1  ''
  'nu12'  0  ''
  'G12'   1  ''
  'G1z'   1  ''
  'G2z'   1  ''
  'Rho'   3  'Density';
  'A1'  8 'Thermal expansion coef';
  'A2'  8 'Thermal expansion coef';
  'T0'       7  'Reference temperature';
  'xt'       1  'Allowable stress tension';
  'xc'       1  'Allowable stress compression';
  'yt'       1  'Allowable stress tension';
  'yc'       1  'Allowable stress compression';
  'Eta' 0 'Loss factor';
  'F12' 0 'F12';
  'STRN' 0 'STRN';
  };

 case 6 % Orthotropic material
 st={'MatId' 0  'sdtweb(''m_elastic'')';
  'Type'  0  '';
  'E1'   1  '';
  'E2'   1  ''
  'E3'   1  ''
  'nu23'  0  ''
  'nu31'  0  ''
  'nu12'  0  ''
  'G23'   1  ''
  'G31'   1  ''
  'G12'   1  ''
  'Rho'   3  'Density';
  'A1'    8  'x Thermal expansion';
  'A2'    8  'y Thermal expansion';
  'A3'    8  'z Thermal expansion';
  'T0'    7  'Reference temperature';
  'Eta'   0  'Loss factor'}; 

otherwise; st={'MatId' 0 'sdtweb(''m_elastic'')'; 'Type', 0, ''};
end

if ~isempty(strfind(Cam,'cell')); out=st; else; out=[st{:,2}]; end

%% #SubType ------------------------------------------------------------------
elseif comstr(Cam,'subtype');[CAM,Cam]=comstr(CAM,8);

st={'1 (Isotropic)','2 (Acoustic fluid)',...
    '3 (3-D anisotropic solid)','4 (2-D anisotropic solid)', ...
    '5 (Orthotropic for shell)','6 (Orthotropic material)'}; 
if carg<=nargin;
 i1=varargin{carg}; carg=carg+1;
 if ischar(i1);out=find(strncmpi(i1,st,1));
 else;
  try; out=st{i1};
  catch; out=sprintf('m_elastic %i',i1);
  end
 end
end

%% #Formula ------------------------------------------------------------------
elseif comstr(Cam,'formula');[CAM,Cam]=comstr(CAM,8);

%% #FormulaEng2dd - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
if comstr(Cam,'eng2dd')

 % formula from m_elastic documentation
 mat=varargin{carg};carg=carg+1;alpha=[];
 if carg<=nargin; RunOpt=varargin{carg};carg=carg+1;else;RunOpt=[];end
 if isfield(RunOpt,'MatStack')&&isfield(RunOpt.MatStack,'E')
  % case with variable property, currently only temperature
  model=varargin{carg};carg=carg+1;
  Case=varargin{carg};carg=carg+1; jGroup=Case.jGroup;

  [EGroup,nGroup]=getegroup(model.Elt);
  cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;

  [ElemF,opt]= feutil('getelemf',model.Elt(EGroup(jGroup),:),jGroup);
  EltConst=integrules(ElemF,-1);

  EGID=opt(1); cEGI = EGroup(jGroup)+1:EGroup(jGroup+1)-1;
  inode=fe_super('node',ElemF);
  NNode=sparse(Case.Node(:,1),1,1:size(Case.Node,1));
  NodePos=model.Elt(cEGI,inode)'; 
  if min(min(NodePos))==0; 
     i1=find(NodePos);NodePos(i1)=NNode(NodePos(i1));
  else;NodePos=reshape(NNode(NodePos),length(inode),length(cEGI));
  end
  NodePos=int32(full(NodePos));  
  if isempty(Case.Stack);i1=[];
  else;i1=find(strcmpi('dofset',Case.Stack(:,1)));
  end;r1=[];
  for j1=1:length(i1)
   DOF=Case.Stack{i1(j1),3}.DOF;
   if length(fe_c(DOF,.20,'ind'))==length(DOF); r1=Case.Stack{i1(j1),3};end
  end
  if ~isempty(r1);r1=feutil('placeindof',model.Node(:,1)+.20,r1);end
  RunOpt.T=zeros(size(EltConst.w,1),length(cEGI));
  for jElt=1:length(cEGI)
    RunOpt.T(:,jElt)=EltConst.N*r1.def(NodePos(:,jElt));
  end  
  RunOpt.T=RunOpt.T(:);
  fprintf('Temperature range %.2f %.2f',[min(RunOpt.T) max(RunOpt.T)])
  % now build the temperature dependent constit for each element
  r1=RunOpt.MatStack;
  cur=[1 0 0];cur(1)=0;
  E=of_time('lininterp',[r1.E.X r1.E.Y],RunOpt.T,cur);
  if ~isfield(r1,'nu');nu=mat(2)*ones(size(RunOpt.T));
  else; nu=of_time('lininterp',[r1.nu.X r1.nu.Y],RunOpt.T,cur);
  end
  if any(~isfinite(E))||any(~isfinite(nu));error('NaN');end
  if isfield(r1,'G')
   G=of_time('lininterp',[r1.G.X r1.G.Y],RunOpt.T,cur);
  else;G=E/2./(1+nu);
  end
  if isfield(r1,'alpha');   
   alpha=of_time('lininterp',[r1.alpha.X r1.alpha.Y],RunOpt.T,cur);
  end
  out1=struct('RunOpt_NodePos',1, ...
      'RunOpt_GstateFcn','fe_stress(''thermalgstate'')');
  % ... 'Case=fe_stress(''thermalgstate'',model,Case);');
 else % standard case with one material
   E=mat(1); nu=mat(2); G=mat(3);
 end
 
 % Basic formula for 3d isotropic elasticity

 r1=nu./(1-nu); % n/(1-n) 
 r2=E.*(1-nu)./(1+nu)./(1-2*nu); % E(1-n)/(1+n)(1-2*n)
 if ~isfinite(r2);error('Full incompressibility is not accepted');end
 out=zeros(36,length(E));
 for j1=1:length(E)
    out([1 2 3 7 8 9 13 14 15],j1)= ...
      r2(j1)*[1 r1(j1) r1(j1)  r1(j1) 1 r1(j1)  r1(j1) r1(j1) 1]';
    out([22 29 36],j1)=G(j1); % G
 end
 if 1==2 % FROM (K,G) #FromKG -3
  K=E./(3*(1-2*nu)); G=E./(2*(1+nu));
  out=zeros(36,length(K));
  for j1=1:length(K)
    out([1 2 3 7 8 9 13 14 15],j1)=...[
     K(j1)*[1 1 1 1 1 1 1 1 1]'+ G(j1)*2/3*[2 -1 -1 -1 2 -1 -1 -1 2]';
    out([22 29 36],j1)=G(j1)*[1 1 1]'; % G
  end
 end
 
 if ~isempty(alpha) % place diagonal matrix of alpha afterwards
  for j1=1:length(alpha)
    out([37 41 45],j1)=alpha(j1); % G
  end
 end

 if isfield(RunOpt,'MatStack')&&isfield(RunOpt.MatStack,'E')
  mat=evalin('caller','mat');
  out=[repmat(mat(5),1,size(out,2)); repmat(mat(7),1,size(out,2)); 
       out];
  if ~isempty(alpha)  % Thermal expansion and ref temp
      out(end+1,:)=mat(9);
  end
  out=reshape(out,size(out,1)*size(EltConst.w,1),length(cEGI));
 else; out=reshape(out,6,6);
 end
    %E=mat(3);n=mat(4);G=E/2/(1+n);
    %inv([1/E -n/E -n/E 0 0 0;-n/E 1/E -n/E 0 0 0;-n/E -n/E 1/E 0 0 0;
    %     0 0 0 1/G 0 0; 0 0 0 0 1/G 0;0 0 0 0 0 1/G])-dd

%% #FormulaOrthoPl - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - - -
elseif comstr(Cam,'orthopl') ;[CAM,Cam]=comstr(CAM,8);%orthotropl. 
 r1=varargin{carg};carg=carg+1;
 %E1,E2,E3, n12,n13,n23, G12,G13,G23 % abaqus style
r1(5)=r1(5)*r1(3)/r1(1); 
r1=r1([1 2 3 6 5 4 9 8 7]);
 r1=m_elastic('formula ortho',r1); r2=[];
 for j1=1:6
  for j2=1:j1
    r2(end+1)=r1(j2,j1); %#ok<AGROW>
  end
 end
 out=[1 fe_mat('m_elastic','US',3) r2];

%% #FormulaOrtho : build D orthotropic material - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'ortho') %orthotro. 

 % From Jones (Mechanics of composite materials) page 38
 % n_ji=E_j/E_i n_ij
 % r1=[E1(1) E2(2) E3(3) n23(4) n31(5) n12(6) G23(7) G31(8) G12(9)]
 r1=varargin{carg};carg=carg+1;
 out=formulaOrtho(r1);
 dd=[1/r1(1)      -r1(6)/r1(1) -r1(5)/r1(3)   0 0 0
     -r1(6)/r1(1) 1/r1(2)      -r1(4)/r1(2)   0 0 0
     -r1(5)/r1(3) -r1(4)/r1(2) 1/r1(3)        0 0 0
     0 0 0                                    1/r1(7) 0 0;
     0 0 0                                    0 1/r1(8) 0;
     0 0 0                                    0 0 1/r1(9)];
 out=pinv(dd);
 if any(diag(out)<0)
  RunOpt=varargin{3};
  if isfield(RunOpt,'MatId'); 
      fprintf('MatId %i check nu for positive dd\n',RunOpt.MatId);
  else
      fprintf('check nu for positive dd\n');
  end
  disp(dd);
  d=diag(dd);
  fprintf('xy %.1f zx %.1f yz %.1f should be <1\n', ...
   [-dd(1,2)/sqrt(d(1)*d(2)) -dd(1,3)/sqrt(d(1)*d(3)) -dd(2,3)/sqrt(d(2)*d(3))])
  [u,v]=eig(out);v=diag(v);v=abs(v); out=u*diag(v)*u';
 end
 
%% #FormulaPlAniso :  Build anisotropic law from 6x6 constitutive matrix - - - -
elseif comstr(Cam,'planiso') 
 [CAM,Cam]=comstr(Cam,8);ID=comstr(CAM,[-1 1]);
 RO=struct('pl',[ID(1) fe_mat('m_elastic','SI',3)]);    
 DD=varargin{carg};carg=carg+1;
 if isfield(DD,'Elt')
  model=DD;
  mpid=feutil('mpid',model);mpid=mpid(ismember(mpid(:,1),ID(1)),:);
  mpid=mpid(1,:);mpid(3:4)=[3 1];
  DD=feutil('getdd',mpid,model);DD=DD.dd;
  RO.prop=feutil(sprintf('getpl%i -struct',ID(1)),model);
  [st,i2,r3]=fe_mat('type',RO.prop.Type);
  if ~strcmpi(st,'m_elastic');error('Mismatch');end
  RO.pl=[ID(1) fe_mat('type',st,i2,3)]; % Preserve unit for output
 end
 if carg>nargin||isempty(varargin{carg});bas=eye(3); % No rotation
 elseif isstruct(varargin{carg}); 
   RO=sdth.sfield('addmissing',varargin{carg},RO);carg=carg+1;
   if isfield(RO,'bas');bas=RO.bas;else; bas=eye(3);end
 else
  bas=varargin{carg};carg=carg+1; % provide rotation basis
  if size(bas,1)==1; bas=reshape(bas(7:15),3,3); 
  elseif size(bas,2)>3 % not a 3x3 basis
   sdtw('_nb','Basis specified does not have a proper format, check bas ID %s and its use in p_solid COORDM',num2str(bas(:,1)'));
   if all(diff(bas(:,1))==0)&&size(bas,2)>=15
    sdtw('_nb','bas ID %i is duplicated, using first occurence',bas(1))
    bas=reshape(bas(1,7:15),3,3); 
   end
  end
 end
 % transformation command
 %st=elem0('tensort',zeros(3),2,[1 1 1;2 2 2;3 3 3;4 2 3;5 3 1;6 1 2]);
 IN=DD;S=[];T=bas; S=MechaTensorT(T,IN);
 %eval(strrep(strrep(st,'[','('),']','+1)'));
  
 if ID==-1; out=reshape(S,6,6); % Return rotated tensor
 else; i1=triu(ones(6));i1=find(i1)';
  out=RO.pl; S=S(i1); out(2+(1:length(S)))=S;
  if isfield(RO,'Rho');      
      st=m_elastic('propertyunittypecell',3);
      out((strcmpi(st(:,1),'Rho')))=RO.Rho; 
  elseif isfield(RO,'prop');
      st=m_elastic('propertyunittypecell',3);
      out((strcmpi(st(:,1),'Rho')))=RO.prop.Rho; 
  end
  %[ID fe_mat('m_elastic','SI',3) S(i1)];
 end
 
%% #FormulaLabToOrtho : build ortho from labels - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'labtoortho')  

    % sdtweb ans2sdt PRXZ
  r1=varargin{carg};carg=carg+1;
  st=lower({'MatId','T2pe','E1','E2','E3','nu23','nu31','nu12','G23', ...
      'G31','G12','rho','eta'});
  r1(1,:)=lower(r1(1,:));
  r1(1,:)=cellfun(@(x)strrep(strrep(strrep(x,'x','1'),'y','2'),'z','3'), ...
      r1(1,:),'uni',0);
  out=(1:size(r1,1)-1)';out(:,2)=fe_mat('m_elastic','US',6);
  for j1=1:size(st,2)
   i2=find(strcmpi(r1(1,:),st{j1}));
   if ~isempty(i2);
       out(:,j1)=vertcat(r1{2:end,i2});
   elseif ismember(st{j1},{'matid','t2pe','eta'})
   %elseif strcmpi(st{j1},'nu12')
   %  dbstack; keyboard;
   elseif strcmpi(st{j1},'nu31')% nu31=nu13/E1*E3
     i2=strcmpi(r1(1,:),'nu13');
     out(:,j1)=vertcat(r1{2:end,i2})./out(:,3).*out(:,5);
   elseif strcmpi(st{j1},'g31')%
     i2=strcmpi(r1(1,:),'g13');
     out(:,j1)=vertcat(r1{2:end,i2});
   else;error('%s',st{j1});
   end
  end

%% #FormulaDdTopro : material entry from sig=dd ep - - - - - - - - - - - - - - - - -
% see formula sdtweb feform#feelas3d
elseif comstr(Cam,'ddtoortho')  

Y=inv(varargin{carg});carg=carg+1; % compliance
%M.Xlab{1,2}={'E1', 'E2', 'E3','G23','G13','G12','nu23','nu13','nu12'}';
out1=struct('E1',1./Y(1,1),'E2',1/Y(2,2),'E3',1/Y(3,3), ...
    'G23',1/Y(4,4),'G13',1/Y(5,5),'G12',1/Y(6,6),'nu23',[],'nu13',[],'nu12',[]);
out1.nu23=-Y(2,3)*out1.E2;
out1.nu13=-Y(1,3)*out1.E1;
out1.nu12=-Y(1,2)*out1.E1;

out=[1 fe_mat('m_elastic','US',6) out1.E1 out1.E2 out1.E3 out1.nu23  ...
    -Y(3,1)*out1.E3 out1.nu12 out1.G23 out1.G13 out1.G12];
else
    
 RO=struct('lambda_mu_e_nu',@(E,nu)[E*nu/(1+nu)/(1-2*nu) E/2/(1+ nu)]);
 if isfield(RO,CAM);out=RO.(CAM);
 else
    error('Formula%s unknown',CAM);
 end
   
end
%% #AtNode ---------------------------------------------------------------------
elseif comstr(Cam,'atnode');[CAM,Cam]=comstr(CAM,7);

if comstr(Cam,'gstate'); 
 %% #AtNodeGstate Build the InfoAtNode group thermal MAP
 % initialized in sdtweb p_solid interpfield   
 model=evalin('caller','model');
 Case=evalin('caller','Case');
 EC=evalin('caller','EC');
 constit=evalin('caller','constit'); integ=evalin('caller','integ');
 RunOpt.cStep=size(constit,1);
 [def,st]=m_elastic('AtNodedef',Case);
 EC.material='Elastic3DNL';
 if size(constit,2)>1; 
  error('You AtNodeGState assume on property per group\n %s',...
      'use model.Elt=feutil(''SeparateByProp'',model.Elt)')
 end
 % build a defT defined at nodes
 cEGI=evalin('caller','cEGI'); %Node=evalin('caller','cEGI');
 r2=stack_get(model,'info','RefTemp','getdata');
 if ~isempty(r2);warning('You should use info,RefT stack entries');
 else; r2=stack_get(model,'info','RefT','getdata');
 end
 if ~isempty(r2); r2=r2(1);def.DOF(end+1)=.20;def.def(end+1,:)=r2;end
 InfoAtNode=elem0('VectFromDirAtNode',Case,def,EC);
 IA=Case.GroupInfo{Case.jGroup,7};
 if ~isempty(IA); InfoAtNode=fe_mknl('mapmerge',InfoAtNode,IA);end
 mat=stack_get(model,'mat'); 
 if size(integ,2)~=1;error('Not implemented for multimaterial group');end
 for j1=1:size(mat,1); 
         if mat{j1,3}.pl(1)==integ(1); mat=mat{j1,3};break;end
 end

 if length(constit)<20 
    %% #IsotropicInterp allow temperature interpolation
     % sdtweb('of_mk_pre.c','thermoelastic isotropic material')
     constit(10:end)=[]; if length(constit)<9;constit(9)=0;end
     if ~isempty(r2);constit(9,:)=r2;
         fprintf('RefTemp set to %g in MatId %i\n',r2,constit(1));
     end % set ref temp
     RA=struct('mat',mat,'type','nonlin_elas_iso','constit',constit, ...
         'IA',InfoAtNode);
     % define nodeE column and offset for each interpolated field
     EC=feval(p_solid('@ctableGen'),EC,InfoAtNode,'un1',RA);
     if sp_util('diag')==8 % set tables for of_mk element wise callback
       EC.material='callback';EC.fHandle='elem0';
       EC.nodeE=zeros(EC.Nnode,6);EC.ke=zeros(EC.Nnode*3);
       EC.defe=zeros(EC.Nnode*3,2);
     elseif 1==2 % Verify values taken by interp
       T=EC.N*InfoAtNode.data(InfoAtNode.NodePos);
       E=T;E(:)=of_time('lininterp',[mat.E.X mat.E.Y],T(:),zeros(1,3));
       G=T;G(:)=of_time('lininterp',[mat.G.X mat.G.Y],T(:),zeros(1,3));
       [min(E(:)) max(E(:)) min(G(:)) max(G(:))]
     end
     evalin('caller','Case.pl=pl;');
     EC.DensPos=5;
 elseif ~isempty(mat) 
    RA=struct('mat',mat,'type','nonlin_elas_aniso','constit',constit, ...
         'IA',InfoAtNode);
    EC=feval(p_solid('@ctableGen'),EC,InfoAtNode,'m_elastic.3',RA);
 end
 assignin('caller','InfoAtNode',InfoAtNode);
 assignin('caller','gstate','');
 if size(constit,1)~=RunOpt.cStep&&size(constit,2)~=1
   pointers=evalin('caller','pointers');
   pointers(7,:)=round(pointers(7,:)*size(constit,1)/RunOpt.cStep);
   assignin('caller','pointers',pointers);
 end
 
 assignin('caller','constit',constit);
 assignin('caller','EC',EC);
 evalin('caller','RunOpt.GstateFcn='''';');
 
elseif comstr(Cam,'def');
 Case=varargin{carg};carg=carg+1;
 ind=find(ismember(comstr(Case.Stack(:,1),-27),'dofset'));
 for j1=1:length(ind)
  r2=Case.Stack{ind(j1),3}; 
  if isequal(unique(round(rem(r2.DOF,1)*100)),20); def=r2;break;end 
 end
 if isempty(def); error('You must specify a thermal DofSet entry');end
 out=def;out1=Case.Stack{ind(j1),2};
end

%% #Test -----------------------------------------------------------------------
elseif comstr(Cam,'test');[CAM,Cam]=comstr(CAM,7);
 eval('t_constit(''elastic'')');
%% #end ------------------------------------------------------------------------
elseif comstr(Cam,'@');out=eval(CAM);
elseif comstr(Cam,'tablecall');out='';
elseif comstr(Cam,'cvs')
    out='$Revision: 1.159 $  $Date: 2021/03/09 07:04:46 $';
else; sdtw('''%s'' not known',CAM);
end % commands

%% #SubFunctions
% ---------------------------------- internal functions
% --------------------------------------------------------------------
%% #formula_homo : implement classical homogeneization formulas
function [out,carg] = formula_homo(CAM,out,varg,carg)

[CAM,Cam]=comstr(CAM,1);

if comstr(Cam,'gibson');[CAM,Cam]=comstr(CAM,7);
%% #Gibson ---------------------------------------------------------------2
% cf Florens thesis : O:\sdtdata\publications\theses\These_Florens.pdf 
% formulas of table 2.4 , p25, geometry page 45

 r1=[];G=[];if carg<=length(varg); r1=varg{carg};carg=carg+1; end
 % E(1) Nu(2) rho(3) a(4) b(5) theta(6) t(7) t'(8)
 if isnumeric(r1)&&length(r1)==8
  E=r1(1);  Nu=r1(2); rho=r1(3);
  a=r1(4); b=r1(5);
  theta=r1(6);  t=r1(7); tp=r1(8);
 else; % default sample aluminum core
  st=[
   'E(72e9#%g#"Young modulus") ' ...
   'G(#%g#"Shear modulus") ' ...
   'Nu(0.3#%g#"Poisson coeff") ' ...
   'rho(2700#%g#"Density of cell wall") ' ...
   'a(.005#%g#"Length of cell side") ' ...
   'b(.005#%g#"Length of sloping cell side") ' ...
   'theta(30#%g#"Angle of cell [deg]") ' ...
   ' t(.0001#%g#"Thickness of simple cell wall") ' ...
   'tp(.0002#%g#"Thickness of double cell wall")' ...
    ];
  [r1,st,CAM]=cingui('paramedit -DoClean',st,{struct,CAM});Cam=lower(CAM);
  st1=fieldnames(r1); 
  for j1=1:length(st1);eval(sprintf('%s=r1.%s;',st1{j1},st1{j1}));end
 end 
 if isempty(G)||G==0;G=E/2/(1+Nu); end
 theta=theta*pi/180; % Rad
  
 %[MatId typ E1 E2 E3 Nu23 Nu31 Nu12 G23 G31 G12 rho a1 a2 a3 T0 eta
  RO=evalin('caller','RO');
  if isempty(RO.punit);RO.punit='US';end
  RO.S_metal=4*a*t+2*b*tp; 
  RO.S_cell=(2*a*cos(theta))*(2*b+2*a*sin(theta));
  out=[1 fe_mat('m_elastic',RO.punit,6) ... % Orthotropic material
       E*(t/b)^3*(a/b+sin(theta))/(cos(theta)^3) ... % Ex
       E*(t/b)^3*(cos(theta))/((a/b+sin(theta))*sin(theta)^2) ... % Ey
       E/b*(2*t+(a/b)*tp)/((a/b+sin(theta))*cos(theta)) ... % Ez
       0 ... % Nu23 defined below
       0 ... % Nu_zx
       (a/b+sin(theta))*sin(theta)/(cos(theta)^2)  ... % Nu_xy
       G*t/b*cos(theta)/(a/b+sin(theta)) ... % G_yz
       G/b*(a/b*tp+2*t*sin(theta)^2)/(2*cos(theta)*(a/b+sin(theta))) ... % G31 Glow/Gup ? G31=G13 ?...
       E/(b^3)*(a/b+sin(theta))/((a/b)^2*(1/tp^3+2*a/b/t^3)*cos(theta)) ... % G12
       rho*RO.S_metal/RO.S_cell; ...
       ];
  % xxx EB Gxz low or up
  if 1==2
   out(6)=Nu*out(3)/out(4); %Nu23=Nu32*E2/E3
   out(6)=min(Nu,.99*sqrt(out(4)/out(5))); % n23<=sqrt(E3/E2)
   out(7)=min(Nu*out(3)/out(4),sqrt(out(5)/out(3))); % n31<=sqrt(E3/E1)
   out(8)=min(out(8),.99*sqrt(out(3)/out(4))); % n31<=sqrt(E3/E1)
  end
  if 1==1 % Enforce positivity
   for j1=1:3
   [a,b,c,dd1]=p_solid('buildconstit 3 1',[1;1],out,p_solid('dbval 1 d3 -3'));
   dd=pinv(dd1.dd);[u,s]=eig(dd);disp(diag(s)')
   dd=u*abs(s)*u'; out(6:8)=[-dd(2,3)/dd(2,2) -dd(3,1)/dd(1) -dd(1,2)/dd(1)]/2;
   end
   %[-dd(3,2)/sqrt(dd(2,2)*dd(3,3)) -out(6)/out(4)
   % -dd(2,1)/sqrt(dd(2,2)*dd(1,1)) 0
   % -dd(3,1)/sqrt(dd(3,3)*dd(1,1)) 0]
  end
  
  out=struct('pl',out,'name','Gibson','type','m_elastic','unit',RO.punit);
%% #LongFiber -2
% sdtweb comp13('rveUD')
elseif comstr(Cam,'longfiber')
    
  error('xxxeb');
Em=model.pl(2,3);                %Matrix Young's modulus
Ef=model.pl(1,3);                %Fiber Young's modulus
Num=model.pl(2,4);               %Matrix Poison's ration
Nuf=model.pl(1,4);               %Fiber Poison's ration

%Introduction of other elasicity moduli to simplify the writing

Gm=Em/(2*(1+Num));                    %matrix shear modulus
Gf=Ef/(2*(1+Nuf));                    %fiber shear modulus
Km=.5*Em/((1+Num)*(1-2*Num));         %matrix lateral compression modulus          
Kf=.5*Ef/((1+Nuf)*(1-2*Nuf));         %fiberlateral compression modulus 
km=Em/(3*(1-2*Num));                  %matrix bulk modulus
kf=Ef/(3*(1-2*Nuf));                  %fiber bulk modulus

%homogenized moduli

%mixture law for the longitudinal Young's mod
El= model.pl(1,3)*RO.vf+model.pl(2,3)*(1-RO.vf);             
Nult=model.pl(1,4)*RO.vf+model.pl(2,4)*(1-RO.vf);
%mixture law for Poisson ratio LT
Glt=Gm*(Gf*(1+RO.vf)+Gm*(1-RO.vf))/(Gf*(1-RO.vf)+Gm*(1+RO.vf));  
Kl=Km+RO.vf/(1/(kf-km+1/3*(Gf-Gm))+(1-RO.vf)/(km+4/3*Gm));
Gtt=Gm*(1+RO.vf/(Gm/(Gf-Gm)+(1-RO.vf)*(km+7/3*Gm)/(2*km+8/3*Gm)));

%homogenized stiffness constants with x axis of symmetry
if 1==2
c(1,1)=El+4*(Nult^2)*Kl;
c(1,2)=2*Kl*Nult;
c(2,2)=Gtt+Kl;
c(2,3)=-Gtt+Kl;
c(6,6)=Glt;

c(3,3)=c(2,2);
c(4,4)=.5*(c(2,2)-c(2,3));
c(5,5)=c(6,6);
c(1,3)=c(1,2);
end
%homogenized stiffness constants with z axis of symmetry
c(3,3)=El+4*(Nult^2)*Kl;
c(1,3)=2*Kl*Nult;
c(1,1)=Gtt+Kl;
c(1,2)=-Gtt+Kl;
c(5,5)=Glt;

c(2,2)=c(1,1);
c(6,6)=.5*(c(1,1)-c(1,2));
c(4,4)=c(5,5);
c(2,3)=c(1,3);

transiso=c-tril(c)+c';
out=transiso;
 
else;    error('Formula%s unknown',CAM);
end
%% #punitCheck
function RO=punitCheck(RO);
 if isfield(RO,'punit')&&~isempty(RO.punit);
 elseif isfield(RO,'unit')&&~isempty(RO.unit)
   fprintf('Assuming properties given in %s units\n',RO.unit);
   RO.punit=RO.unit;
 else; RO.punit='US';
 end
 
%% #formulaOrthoFun 
function out=formulaOrtho(r1,ver);
 if nargin==2; % abaqus version sdt expects nu23,nu31,nu12, moldflow gives nu23,nu13,nu12
   r1(5)=r1(5)/r1(1)*r1(3);
 end
 dd=[1/r1(1)      -r1(6)/r1(1) -r1(5)/r1(3)   0 0 0
     -r1(6)/r1(1) 1/r1(2)      -r1(4)/r1(2)   0 0 0
     -r1(5)/r1(3) -r1(4)/r1(2) 1/r1(3)        0 0 0
     0 0 0                                    1/r1(7) 0 0;
     0 0 0                                    0 1/r1(8) 0;
     0 0 0                                    0 0 1/r1(9)];
 out=pinv(dd);


%% #MechaTensorT : formula for coordinate transformation
% TGL basis, IN is constitutive law, sdtweb elem0('tensort')
function S=MechaTensorT(T,IN);

if 1==2
 S.tGL=cell(6); S.tLG=cell(6);
 % sG_il = tGL_(il,jk) sL_(jk) = vGL_ij vGL_lk  sigmaL_jk 
 % eL_il = tLG_ij,kl eG_jk = vGL_ji  vGL_kl  eG_jk
 for ai=1:3;for aj=1:3;for ak=1:3; for al=1:3 %#ok<ALIGN>
  [un1,il]=ismember([ai al],[1 1;2 2;3 3;2 3;3 1;1 2],'rows');if il==0;continue;end
  [un1,jk]=ismember(sort([aj ak]),sort([1 1;2 2;3 3;2 3;3 1;1 2],2),'rows');
  st1=sprintf('T(%i,%i)*T(%i,%i)',ai,aj,al,ak);
  if isequal(st1,S.tGL{il,jk});S.tGL{il,jk}=['2*' st1];
  else;S.tGL{il,jk}=[S.tGL{il,jk} '+' st1];
  end
  st1=sprintf('T(%i,%i)*T(%i,%i)',aj,ai,ak,al);
  if isequal(st1,S.tLG{il,jk});S.tLG{il,jk}=['2*(' st1 ')'];
  else;S.tLG{il,jk}=[S.tLG{il,jk} '+' st1];
  end
 end;end;end;end
 S.tGL=cellfun(@(x)x(2:end),S.tGL,'uni',0);
 S.tLG=cellfun(@(x)x(2:end),S.tLG,'uni',0);
 for j1=1:3;for j2=4:6;S.tLG{j1,j2}=sprintf('.5*(%s)',S.tLG{j1,j2});end;end
 for j1=1:3;for j2=4:6;S.tLG{j2,j1}=sprintf('2*(%s)',S.tLG{j2,j1});end;end
 fprintf('\n\n tGL=[\n')
 z=S.tGL';fprintf('%s %s %s %s %s %s \n',z{:})
 fprintf('\n\n tLG=[\n')
 z=S.tLG';fprintf('%s %s %s %s %s %s \n',z{:})
 %% Now 9 x 9 transform sG_il=vGL_ij vGL_lk  sigmaL_jk 
 S.tGL=cell(9,9);
 for ai=1:3;for aj=1:3;for ak=1:3; for al=1:3 %#ok<ALIGN>
  S.tGL{ai+3*(al-1),aj+3*(ak-1)}=sprintf('T(%i,%i)*T(%i,%i)',ai,aj,al,ak);
 end;end;end;end
 fprintf('\n\n tGL=[\n')
 z=S.tGL';fprintf('%s %s %s %s %s %s %s %s %s \n',z{:});fprintf(']\n\n ')
 
 
elseif nargin<=1 % Return tensor transforms (see Zhilong Liu report) PJE/cristal
S=struct('vGL',@(T)T,'vLG',@(T)T', ...
 'tGL',@(T)[T(1,1)*T(1,1) T(1,2)*T(1,2) T(1,3)*T(1,3) T(1,2)*T(1,3)+T(1,3)*T(1,2) T(1,1)*T(1,3)+T(1,3)*T(1,1) T(1,1)*T(1,2)+T(1,2)*T(1,1) 
T(2,1)*T(2,1) T(2,2)*T(2,2) T(2,3)*T(2,3) T(2,2)*T(2,3)+T(2,3)*T(2,2) T(2,1)*T(2,3)+T(2,3)*T(2,1) T(2,1)*T(2,2)+T(2,2)*T(2,1) 
T(3,1)*T(3,1) T(3,2)*T(3,2) T(3,3)*T(3,3) T(3,2)*T(3,3)+T(3,3)*T(3,2) T(3,1)*T(3,3)+T(3,3)*T(3,1) T(3,1)*T(3,2)+T(3,2)*T(3,1) 
T(2,1)*T(3,1) T(2,2)*T(3,2) T(2,3)*T(3,3) T(2,2)*T(3,3)+T(2,3)*T(3,2) T(2,1)*T(3,3)+T(2,3)*T(3,1) T(2,1)*T(3,2)+T(2,2)*T(3,1) 
T(3,1)*T(1,1) T(3,2)*T(1,2) T(3,3)*T(1,3) T(3,2)*T(1,3)+T(3,3)*T(1,2) T(3,1)*T(1,3)+T(3,3)*T(1,1) T(3,1)*T(1,2)+T(3,2)*T(1,1) 
T(1,1)*T(2,1) T(1,2)*T(2,2) T(1,3)*T(2,3) T(1,2)*T(2,3)+T(1,3)*T(2,2) T(1,1)*T(2,3)+T(1,3)*T(2,1) T(1,1)*T(2,2)+T(1,2)*T(2,1)], ...
  'tLG',@(T)[
T(1,1)*T(1,1) T(2,1)*T(2,1) T(3,1)*T(3,1) .5*(T(2,1)*T(3,1)+T(3,1)*T(2,1)) .5*(T(1,1)*T(3,1)+T(3,1)*T(1,1)) .5*(T(1,1)*T(2,1)+T(2,1)*T(1,1)) 
T(1,2)*T(1,2) T(2,2)*T(2,2) T(3,2)*T(3,2) .5*(T(2,2)*T(3,2)+T(3,2)*T(2,2)) .5*(T(1,2)*T(3,2)+T(3,2)*T(1,2)) .5*(T(1,2)*T(2,2)+T(2,2)*T(1,2)) 
T(1,3)*T(1,3) T(2,3)*T(2,3) T(3,3)*T(3,3) .5*(T(2,3)*T(3,3)+T(3,3)*T(2,3)) .5*(T(1,3)*T(3,3)+T(3,3)*T(1,3)) .5*(T(1,3)*T(2,3)+T(2,3)*T(1,3)) 
2*(T(1,2)*T(1,3)) 2*(T(2,2)*T(2,3)) 2*(T(3,2)*T(3,3)) T(2,2)*T(3,3)+T(3,2)*T(2,3) T(1,2)*T(3,3)+T(3,2)*T(1,3) T(1,2)*T(2,3)+T(2,2)*T(1,3) 
2*(T(1,3)*T(1,1)) 2*(T(2,3)*T(2,1)) 2*(T(3,3)*T(3,1)) T(2,3)*T(3,1)+T(3,3)*T(2,1) T(1,3)*T(3,1)+T(3,3)*T(1,1) T(1,3)*T(2,1)+T(2,3)*T(1,1) 
2*(T(1,1)*T(1,2)) 2*(T(2,1)*T(2,2)) 2*(T(3,1)*T(3,2)) T(2,1)*T(3,2)+T(3,1)*T(2,2) T(1,1)*T(3,2)+T(3,1)*T(1,2) T(1,1)*T(2,2)+T(2,1)*T(1,2) 
     ]);
 % sig=[s11 s22 s33 s23 s31 s12]
 if nargin==1
  T(abs(T)<1e-15)=0;
  S.vGL=T; S.vLG=S.vLG(T); S.tGL=S.tGL(T); S.tLG=S.tLG(T);
 end
 % sig_G= tGL * sig_L
 % epsi=[e11 e22 e33 2*e23 2*e13 2*e12]
  % S.tLG*S.tGL
else
    % see of_mk_pre.c TransformLambda
S(1)=((T(1)*T(1))*(T(1)*T(1))*IN(1)+(T(1)*T(1))*(T(4)*T(4))*IN(7)+(T(1)*T(1))*(T(7)*T(7))*IN(13)+(T(1)*T(1))*(T(4)*T(7)+T(7)*T(4))*IN(19)+(T(1)*T(1))*(T(7)*T(1)+T(1)*T(7))*IN(25)+(T(1)*T(1))*(T(1)*T(4)+T(4)*T(1))*IN(31)+(T(4)*T(4))*(T(1)*T(1))*IN(2)+(T(4)*T(4))*(T(4)*T(4))*IN(8)+(T(4)*T(4))*(T(7)*T(7))*IN(14)+(T(4)*T(4))*(T(4)*T(7)+T(7)*T(4))*IN(20)+(T(4)*T(4))*(T(7)*T(1)+T(1)*T(7))*IN(26)+(T(4)*T(4))*(T(1)*T(4)+T(4)*T(1))*IN(32)+(T(7)*T(7))*(T(1)*T(1))*IN(3)+(T(7)*T(7))*(T(4)*T(4))*IN(9)+(T(7)*T(7))*(T(7)*T(7))*IN(15)+(T(7)*T(7))*(T(4)*T(7)+T(7)*T(4))*IN(21)+(T(7)*T(7))*(T(7)*T(1)+T(1)*T(7))*IN(27)+(T(7)*T(7))*(T(1)*T(4)+T(4)*T(1))*IN(33)+(T(4)*T(7)+T(7)*T(4))*(T(1)*T(1))*IN(4)+(T(4)*T(7)+T(7)*T(4))*(T(4)*T(4))*IN(10)+(T(4)*T(7)+T(7)*T(4))*(T(7)*T(7))*IN(16)+(T(4)*T(7)+T(7)*T(4))*(T(4)*T(7)+T(7)*T(4))*IN(22)+(T(4)*T(7)+T(7)*T(4))*(T(7)*T(1)+T(1)*T(7))*IN(28)+(T(4)*T(7)+T(7)*T(4))*(T(1)*T(4)+T(4)*T(1))*IN(34)+(T(7)*T(1)+T(1)*T(7))*(T(1)*T(1))*IN(5)+(T(7)*T(1)+T(1)*T(7))*(T(4)*T(4))*IN(11)+(T(7)*T(1)+T(1)*T(7))*(T(7)*T(7))*IN(17)+(T(7)*T(1)+T(1)*T(7))*(T(4)*T(7)+T(7)*T(4))*IN(23)+(T(7)*T(1)+T(1)*T(7))*(T(7)*T(1)+T(1)*T(7))*IN(29)+(T(7)*T(1)+T(1)*T(7))*(T(1)*T(4)+T(4)*T(1))*IN(35)+(T(1)*T(4)+T(4)*T(1))*(T(1)*T(1))*IN(6)+(T(1)*T(4)+T(4)*T(1))*(T(4)*T(4))*IN(12)+(T(1)*T(4)+T(4)*T(1))*(T(7)*T(7))*IN(18)+(T(1)*T(4)+T(4)*T(1))*(T(4)*T(7)+T(7)*T(4))*IN(24)+(T(1)*T(4)+T(4)*T(1))*(T(7)*T(1)+T(1)*T(7))*IN(30)+(T(1)*T(4)+T(4)*T(1))*(T(1)*T(4)+T(4)*T(1))*IN(36));
S(7)=((T(1)*T(1))*(T(2)*T(2))*IN(1)+(T(1)*T(1))*(T(5)*T(5))*IN(7)+(T(1)*T(1))*(T(8)*T(8))*IN(13)+(T(1)*T(1))*(T(5)*T(8)+T(8)*T(5))*IN(19)+(T(1)*T(1))*(T(8)*T(2)+T(2)*T(8))*IN(25)+(T(1)*T(1))*(T(2)*T(5)+T(5)*T(2))*IN(31)+(T(4)*T(4))*(T(2)*T(2))*IN(2)+(T(4)*T(4))*(T(5)*T(5))*IN(8)+(T(4)*T(4))*(T(8)*T(8))*IN(14)+(T(4)*T(4))*(T(5)*T(8)+T(8)*T(5))*IN(20)+(T(4)*T(4))*(T(8)*T(2)+T(2)*T(8))*IN(26)+(T(4)*T(4))*(T(2)*T(5)+T(5)*T(2))*IN(32)+(T(7)*T(7))*(T(2)*T(2))*IN(3)+(T(7)*T(7))*(T(5)*T(5))*IN(9)+(T(7)*T(7))*(T(8)*T(8))*IN(15)+(T(7)*T(7))*(T(5)*T(8)+T(8)*T(5))*IN(21)+(T(7)*T(7))*(T(8)*T(2)+T(2)*T(8))*IN(27)+(T(7)*T(7))*(T(2)*T(5)+T(5)*T(2))*IN(33)+(T(4)*T(7)+T(7)*T(4))*(T(2)*T(2))*IN(4)+(T(4)*T(7)+T(7)*T(4))*(T(5)*T(5))*IN(10)+(T(4)*T(7)+T(7)*T(4))*(T(8)*T(8))*IN(16)+(T(4)*T(7)+T(7)*T(4))*(T(5)*T(8)+T(8)*T(5))*IN(22)+(T(4)*T(7)+T(7)*T(4))*(T(8)*T(2)+T(2)*T(8))*IN(28)+(T(4)*T(7)+T(7)*T(4))*(T(2)*T(5)+T(5)*T(2))*IN(34)+(T(7)*T(1)+T(1)*T(7))*(T(2)*T(2))*IN(5)+(T(7)*T(1)+T(1)*T(7))*(T(5)*T(5))*IN(11)+(T(7)*T(1)+T(1)*T(7))*(T(8)*T(8))*IN(17)+(T(7)*T(1)+T(1)*T(7))*(T(5)*T(8)+T(8)*T(5))*IN(23)+(T(7)*T(1)+T(1)*T(7))*(T(8)*T(2)+T(2)*T(8))*IN(29)+(T(7)*T(1)+T(1)*T(7))*(T(2)*T(5)+T(5)*T(2))*IN(35)+(T(1)*T(4)+T(4)*T(1))*(T(2)*T(2))*IN(6)+(T(1)*T(4)+T(4)*T(1))*(T(5)*T(5))*IN(12)+(T(1)*T(4)+T(4)*T(1))*(T(8)*T(8))*IN(18)+(T(1)*T(4)+T(4)*T(1))*(T(5)*T(8)+T(8)*T(5))*IN(24)+(T(1)*T(4)+T(4)*T(1))*(T(8)*T(2)+T(2)*T(8))*IN(30)+(T(1)*T(4)+T(4)*T(1))*(T(2)*T(5)+T(5)*T(2))*IN(36));
S(13)=((T(1)*T(1))*(T(3)*T(3))*IN(1)+(T(1)*T(1))*(T(6)*T(6))*IN(7)+(T(1)*T(1))*(T(9)*T(9))*IN(13)+(T(1)*T(1))*(T(6)*T(9)+T(9)*T(6))*IN(19)+(T(1)*T(1))*(T(9)*T(3)+T(3)*T(9))*IN(25)+(T(1)*T(1))*(T(3)*T(6)+T(6)*T(3))*IN(31)+(T(4)*T(4))*(T(3)*T(3))*IN(2)+(T(4)*T(4))*(T(6)*T(6))*IN(8)+(T(4)*T(4))*(T(9)*T(9))*IN(14)+(T(4)*T(4))*(T(6)*T(9)+T(9)*T(6))*IN(20)+(T(4)*T(4))*(T(9)*T(3)+T(3)*T(9))*IN(26)+(T(4)*T(4))*(T(3)*T(6)+T(6)*T(3))*IN(32)+(T(7)*T(7))*(T(3)*T(3))*IN(3)+(T(7)*T(7))*(T(6)*T(6))*IN(9)+(T(7)*T(7))*(T(9)*T(9))*IN(15)+(T(7)*T(7))*(T(6)*T(9)+T(9)*T(6))*IN(21)+(T(7)*T(7))*(T(9)*T(3)+T(3)*T(9))*IN(27)+(T(7)*T(7))*(T(3)*T(6)+T(6)*T(3))*IN(33)+(T(4)*T(7)+T(7)*T(4))*(T(3)*T(3))*IN(4)+(T(4)*T(7)+T(7)*T(4))*(T(6)*T(6))*IN(10)+(T(4)*T(7)+T(7)*T(4))*(T(9)*T(9))*IN(16)+(T(4)*T(7)+T(7)*T(4))*(T(6)*T(9)+T(9)*T(6))*IN(22)+(T(4)*T(7)+T(7)*T(4))*(T(9)*T(3)+T(3)*T(9))*IN(28)+(T(4)*T(7)+T(7)*T(4))*(T(3)*T(6)+T(6)*T(3))*IN(34)+(T(7)*T(1)+T(1)*T(7))*(T(3)*T(3))*IN(5)+(T(7)*T(1)+T(1)*T(7))*(T(6)*T(6))*IN(11)+(T(7)*T(1)+T(1)*T(7))*(T(9)*T(9))*IN(17)+(T(7)*T(1)+T(1)*T(7))*(T(6)*T(9)+T(9)*T(6))*IN(23)+(T(7)*T(1)+T(1)*T(7))*(T(9)*T(3)+T(3)*T(9))*IN(29)+(T(7)*T(1)+T(1)*T(7))*(T(3)*T(6)+T(6)*T(3))*IN(35)+(T(1)*T(4)+T(4)*T(1))*(T(3)*T(3))*IN(6)+(T(1)*T(4)+T(4)*T(1))*(T(6)*T(6))*IN(12)+(T(1)*T(4)+T(4)*T(1))*(T(9)*T(9))*IN(18)+(T(1)*T(4)+T(4)*T(1))*(T(6)*T(9)+T(9)*T(6))*IN(24)+(T(1)*T(4)+T(4)*T(1))*(T(9)*T(3)+T(3)*T(9))*IN(30)+(T(1)*T(4)+T(4)*T(1))*(T(3)*T(6)+T(6)*T(3))*IN(36));
S(19)=((T(1)*T(1))*(T(2)*T(3))*IN(1)+(T(1)*T(1))*(T(5)*T(6))*IN(7)+(T(1)*T(1))*(T(8)*T(9))*IN(13)+(T(1)*T(1))*(T(5)*T(9)+T(8)*T(6))*IN(19)+(T(1)*T(1))*(T(8)*T(3)+T(2)*T(9))*IN(25)+(T(1)*T(1))*(T(2)*T(6)+T(5)*T(3))*IN(31)+(T(4)*T(4))*(T(2)*T(3))*IN(2)+(T(4)*T(4))*(T(5)*T(6))*IN(8)+(T(4)*T(4))*(T(8)*T(9))*IN(14)+(T(4)*T(4))*(T(5)*T(9)+T(8)*T(6))*IN(20)+(T(4)*T(4))*(T(8)*T(3)+T(2)*T(9))*IN(26)+(T(4)*T(4))*(T(2)*T(6)+T(5)*T(3))*IN(32)+(T(7)*T(7))*(T(2)*T(3))*IN(3)+(T(7)*T(7))*(T(5)*T(6))*IN(9)+(T(7)*T(7))*(T(8)*T(9))*IN(15)+(T(7)*T(7))*(T(5)*T(9)+T(8)*T(6))*IN(21)+(T(7)*T(7))*(T(8)*T(3)+T(2)*T(9))*IN(27)+(T(7)*T(7))*(T(2)*T(6)+T(5)*T(3))*IN(33)+(T(4)*T(7)+T(7)*T(4))*(T(2)*T(3))*IN(4)+(T(4)*T(7)+T(7)*T(4))*(T(5)*T(6))*IN(10)+(T(4)*T(7)+T(7)*T(4))*(T(8)*T(9))*IN(16)+(T(4)*T(7)+T(7)*T(4))*(T(5)*T(9)+T(8)*T(6))*IN(22)+(T(4)*T(7)+T(7)*T(4))*(T(8)*T(3)+T(2)*T(9))*IN(28)+(T(4)*T(7)+T(7)*T(4))*(T(2)*T(6)+T(5)*T(3))*IN(34)+(T(7)*T(1)+T(1)*T(7))*(T(2)*T(3))*IN(5)+(T(7)*T(1)+T(1)*T(7))*(T(5)*T(6))*IN(11)+(T(7)*T(1)+T(1)*T(7))*(T(8)*T(9))*IN(17)+(T(7)*T(1)+T(1)*T(7))*(T(5)*T(9)+T(8)*T(6))*IN(23)+(T(7)*T(1)+T(1)*T(7))*(T(8)*T(3)+T(2)*T(9))*IN(29)+(T(7)*T(1)+T(1)*T(7))*(T(2)*T(6)+T(5)*T(3))*IN(35)+(T(1)*T(4)+T(4)*T(1))*(T(2)*T(3))*IN(6)+(T(1)*T(4)+T(4)*T(1))*(T(5)*T(6))*IN(12)+(T(1)*T(4)+T(4)*T(1))*(T(8)*T(9))*IN(18)+(T(1)*T(4)+T(4)*T(1))*(T(5)*T(9)+T(8)*T(6))*IN(24)+(T(1)*T(4)+T(4)*T(1))*(T(8)*T(3)+T(2)*T(9))*IN(30)+(T(1)*T(4)+T(4)*T(1))*(T(2)*T(6)+T(5)*T(3))*IN(36));
S(25)=((T(1)*T(1))*(T(3)*T(1))*IN(1)+(T(1)*T(1))*(T(6)*T(4))*IN(7)+(T(1)*T(1))*(T(9)*T(7))*IN(13)+(T(1)*T(1))*(T(6)*T(7)+T(9)*T(4))*IN(19)+(T(1)*T(1))*(T(9)*T(1)+T(3)*T(7))*IN(25)+(T(1)*T(1))*(T(3)*T(4)+T(6)*T(1))*IN(31)+(T(4)*T(4))*(T(3)*T(1))*IN(2)+(T(4)*T(4))*(T(6)*T(4))*IN(8)+(T(4)*T(4))*(T(9)*T(7))*IN(14)+(T(4)*T(4))*(T(6)*T(7)+T(9)*T(4))*IN(20)+(T(4)*T(4))*(T(9)*T(1)+T(3)*T(7))*IN(26)+(T(4)*T(4))*(T(3)*T(4)+T(6)*T(1))*IN(32)+(T(7)*T(7))*(T(3)*T(1))*IN(3)+(T(7)*T(7))*(T(6)*T(4))*IN(9)+(T(7)*T(7))*(T(9)*T(7))*IN(15)+(T(7)*T(7))*(T(6)*T(7)+T(9)*T(4))*IN(21)+(T(7)*T(7))*(T(9)*T(1)+T(3)*T(7))*IN(27)+(T(7)*T(7))*(T(3)*T(4)+T(6)*T(1))*IN(33)+(T(4)*T(7)+T(7)*T(4))*(T(3)*T(1))*IN(4)+(T(4)*T(7)+T(7)*T(4))*(T(6)*T(4))*IN(10)+(T(4)*T(7)+T(7)*T(4))*(T(9)*T(7))*IN(16)+(T(4)*T(7)+T(7)*T(4))*(T(6)*T(7)+T(9)*T(4))*IN(22)+(T(4)*T(7)+T(7)*T(4))*(T(9)*T(1)+T(3)*T(7))*IN(28)+(T(4)*T(7)+T(7)*T(4))*(T(3)*T(4)+T(6)*T(1))*IN(34)+(T(7)*T(1)+T(1)*T(7))*(T(3)*T(1))*IN(5)+(T(7)*T(1)+T(1)*T(7))*(T(6)*T(4))*IN(11)+(T(7)*T(1)+T(1)*T(7))*(T(9)*T(7))*IN(17)+(T(7)*T(1)+T(1)*T(7))*(T(6)*T(7)+T(9)*T(4))*IN(23)+(T(7)*T(1)+T(1)*T(7))*(T(9)*T(1)+T(3)*T(7))*IN(29)+(T(7)*T(1)+T(1)*T(7))*(T(3)*T(4)+T(6)*T(1))*IN(35)+(T(1)*T(4)+T(4)*T(1))*(T(3)*T(1))*IN(6)+(T(1)*T(4)+T(4)*T(1))*(T(6)*T(4))*IN(12)+(T(1)*T(4)+T(4)*T(1))*(T(9)*T(7))*IN(18)+(T(1)*T(4)+T(4)*T(1))*(T(6)*T(7)+T(9)*T(4))*IN(24)+(T(1)*T(4)+T(4)*T(1))*(T(9)*T(1)+T(3)*T(7))*IN(30)+(T(1)*T(4)+T(4)*T(1))*(T(3)*T(4)+T(6)*T(1))*IN(36));
S(31)=((T(1)*T(1))*(T(1)*T(2))*IN(1)+(T(1)*T(1))*(T(4)*T(5))*IN(7)+(T(1)*T(1))*(T(7)*T(8))*IN(13)+(T(1)*T(1))*(T(4)*T(8)+T(7)*T(5))*IN(19)+(T(1)*T(1))*(T(7)*T(2)+T(1)*T(8))*IN(25)+(T(1)*T(1))*(T(1)*T(5)+T(4)*T(2))*IN(31)+(T(4)*T(4))*(T(1)*T(2))*IN(2)+(T(4)*T(4))*(T(4)*T(5))*IN(8)+(T(4)*T(4))*(T(7)*T(8))*IN(14)+(T(4)*T(4))*(T(4)*T(8)+T(7)*T(5))*IN(20)+(T(4)*T(4))*(T(7)*T(2)+T(1)*T(8))*IN(26)+(T(4)*T(4))*(T(1)*T(5)+T(4)*T(2))*IN(32)+(T(7)*T(7))*(T(1)*T(2))*IN(3)+(T(7)*T(7))*(T(4)*T(5))*IN(9)+(T(7)*T(7))*(T(7)*T(8))*IN(15)+(T(7)*T(7))*(T(4)*T(8)+T(7)*T(5))*IN(21)+(T(7)*T(7))*(T(7)*T(2)+T(1)*T(8))*IN(27)+(T(7)*T(7))*(T(1)*T(5)+T(4)*T(2))*IN(33)+(T(4)*T(7)+T(7)*T(4))*(T(1)*T(2))*IN(4)+(T(4)*T(7)+T(7)*T(4))*(T(4)*T(5))*IN(10)+(T(4)*T(7)+T(7)*T(4))*(T(7)*T(8))*IN(16)+(T(4)*T(7)+T(7)*T(4))*(T(4)*T(8)+T(7)*T(5))*IN(22)+(T(4)*T(7)+T(7)*T(4))*(T(7)*T(2)+T(1)*T(8))*IN(28)+(T(4)*T(7)+T(7)*T(4))*(T(1)*T(5)+T(4)*T(2))*IN(34)+(T(7)*T(1)+T(1)*T(7))*(T(1)*T(2))*IN(5)+(T(7)*T(1)+T(1)*T(7))*(T(4)*T(5))*IN(11)+(T(7)*T(1)+T(1)*T(7))*(T(7)*T(8))*IN(17)+(T(7)*T(1)+T(1)*T(7))*(T(4)*T(8)+T(7)*T(5))*IN(23)+(T(7)*T(1)+T(1)*T(7))*(T(7)*T(2)+T(1)*T(8))*IN(29)+(T(7)*T(1)+T(1)*T(7))*(T(1)*T(5)+T(4)*T(2))*IN(35)+(T(1)*T(4)+T(4)*T(1))*(T(1)*T(2))*IN(6)+(T(1)*T(4)+T(4)*T(1))*(T(4)*T(5))*IN(12)+(T(1)*T(4)+T(4)*T(1))*(T(7)*T(8))*IN(18)+(T(1)*T(4)+T(4)*T(1))*(T(4)*T(8)+T(7)*T(5))*IN(24)+(T(1)*T(4)+T(4)*T(1))*(T(7)*T(2)+T(1)*T(8))*IN(30)+(T(1)*T(4)+T(4)*T(1))*(T(1)*T(5)+T(4)*T(2))*IN(36));
S(2)=((T(2)*T(2))*(T(1)*T(1))*IN(1)+(T(2)*T(2))*(T(4)*T(4))*IN(7)+(T(2)*T(2))*(T(7)*T(7))*IN(13)+(T(2)*T(2))*(T(4)*T(7)+T(7)*T(4))*IN(19)+(T(2)*T(2))*(T(7)*T(1)+T(1)*T(7))*IN(25)+(T(2)*T(2))*(T(1)*T(4)+T(4)*T(1))*IN(31)+(T(5)*T(5))*(T(1)*T(1))*IN(2)+(T(5)*T(5))*(T(4)*T(4))*IN(8)+(T(5)*T(5))*(T(7)*T(7))*IN(14)+(T(5)*T(5))*(T(4)*T(7)+T(7)*T(4))*IN(20)+(T(5)*T(5))*(T(7)*T(1)+T(1)*T(7))*IN(26)+(T(5)*T(5))*(T(1)*T(4)+T(4)*T(1))*IN(32)+(T(8)*T(8))*(T(1)*T(1))*IN(3)+(T(8)*T(8))*(T(4)*T(4))*IN(9)+(T(8)*T(8))*(T(7)*T(7))*IN(15)+(T(8)*T(8))*(T(4)*T(7)+T(7)*T(4))*IN(21)+(T(8)*T(8))*(T(7)*T(1)+T(1)*T(7))*IN(27)+(T(8)*T(8))*(T(1)*T(4)+T(4)*T(1))*IN(33)+(T(5)*T(8)+T(8)*T(5))*(T(1)*T(1))*IN(4)+(T(5)*T(8)+T(8)*T(5))*(T(4)*T(4))*IN(10)+(T(5)*T(8)+T(8)*T(5))*(T(7)*T(7))*IN(16)+(T(5)*T(8)+T(8)*T(5))*(T(4)*T(7)+T(7)*T(4))*IN(22)+(T(5)*T(8)+T(8)*T(5))*(T(7)*T(1)+T(1)*T(7))*IN(28)+(T(5)*T(8)+T(8)*T(5))*(T(1)*T(4)+T(4)*T(1))*IN(34)+(T(8)*T(2)+T(2)*T(8))*(T(1)*T(1))*IN(5)+(T(8)*T(2)+T(2)*T(8))*(T(4)*T(4))*IN(11)+(T(8)*T(2)+T(2)*T(8))*(T(7)*T(7))*IN(17)+(T(8)*T(2)+T(2)*T(8))*(T(4)*T(7)+T(7)*T(4))*IN(23)+(T(8)*T(2)+T(2)*T(8))*(T(7)*T(1)+T(1)*T(7))*IN(29)+(T(8)*T(2)+T(2)*T(8))*(T(1)*T(4)+T(4)*T(1))*IN(35)+(T(2)*T(5)+T(5)*T(2))*(T(1)*T(1))*IN(6)+(T(2)*T(5)+T(5)*T(2))*(T(4)*T(4))*IN(12)+(T(2)*T(5)+T(5)*T(2))*(T(7)*T(7))*IN(18)+(T(2)*T(5)+T(5)*T(2))*(T(4)*T(7)+T(7)*T(4))*IN(24)+(T(2)*T(5)+T(5)*T(2))*(T(7)*T(1)+T(1)*T(7))*IN(30)+(T(2)*T(5)+T(5)*T(2))*(T(1)*T(4)+T(4)*T(1))*IN(36));
S(8)=((T(2)*T(2))*(T(2)*T(2))*IN(1)+(T(2)*T(2))*(T(5)*T(5))*IN(7)+(T(2)*T(2))*(T(8)*T(8))*IN(13)+(T(2)*T(2))*(T(5)*T(8)+T(8)*T(5))*IN(19)+(T(2)*T(2))*(T(8)*T(2)+T(2)*T(8))*IN(25)+(T(2)*T(2))*(T(2)*T(5)+T(5)*T(2))*IN(31)+(T(5)*T(5))*(T(2)*T(2))*IN(2)+(T(5)*T(5))*(T(5)*T(5))*IN(8)+(T(5)*T(5))*(T(8)*T(8))*IN(14)+(T(5)*T(5))*(T(5)*T(8)+T(8)*T(5))*IN(20)+(T(5)*T(5))*(T(8)*T(2)+T(2)*T(8))*IN(26)+(T(5)*T(5))*(T(2)*T(5)+T(5)*T(2))*IN(32)+(T(8)*T(8))*(T(2)*T(2))*IN(3)+(T(8)*T(8))*(T(5)*T(5))*IN(9)+(T(8)*T(8))*(T(8)*T(8))*IN(15)+(T(8)*T(8))*(T(5)*T(8)+T(8)*T(5))*IN(21)+(T(8)*T(8))*(T(8)*T(2)+T(2)*T(8))*IN(27)+(T(8)*T(8))*(T(2)*T(5)+T(5)*T(2))*IN(33)+(T(5)*T(8)+T(8)*T(5))*(T(2)*T(2))*IN(4)+(T(5)*T(8)+T(8)*T(5))*(T(5)*T(5))*IN(10)+(T(5)*T(8)+T(8)*T(5))*(T(8)*T(8))*IN(16)+(T(5)*T(8)+T(8)*T(5))*(T(5)*T(8)+T(8)*T(5))*IN(22)+(T(5)*T(8)+T(8)*T(5))*(T(8)*T(2)+T(2)*T(8))*IN(28)+(T(5)*T(8)+T(8)*T(5))*(T(2)*T(5)+T(5)*T(2))*IN(34)+(T(8)*T(2)+T(2)*T(8))*(T(2)*T(2))*IN(5)+(T(8)*T(2)+T(2)*T(8))*(T(5)*T(5))*IN(11)+(T(8)*T(2)+T(2)*T(8))*(T(8)*T(8))*IN(17)+(T(8)*T(2)+T(2)*T(8))*(T(5)*T(8)+T(8)*T(5))*IN(23)+(T(8)*T(2)+T(2)*T(8))*(T(8)*T(2)+T(2)*T(8))*IN(29)+(T(8)*T(2)+T(2)*T(8))*(T(2)*T(5)+T(5)*T(2))*IN(35)+(T(2)*T(5)+T(5)*T(2))*(T(2)*T(2))*IN(6)+(T(2)*T(5)+T(5)*T(2))*(T(5)*T(5))*IN(12)+(T(2)*T(5)+T(5)*T(2))*(T(8)*T(8))*IN(18)+(T(2)*T(5)+T(5)*T(2))*(T(5)*T(8)+T(8)*T(5))*IN(24)+(T(2)*T(5)+T(5)*T(2))*(T(8)*T(2)+T(2)*T(8))*IN(30)+(T(2)*T(5)+T(5)*T(2))*(T(2)*T(5)+T(5)*T(2))*IN(36));
S(14)=((T(2)*T(2))*(T(3)*T(3))*IN(1)+(T(2)*T(2))*(T(6)*T(6))*IN(7)+(T(2)*T(2))*(T(9)*T(9))*IN(13)+(T(2)*T(2))*(T(6)*T(9)+T(9)*T(6))*IN(19)+(T(2)*T(2))*(T(9)*T(3)+T(3)*T(9))*IN(25)+(T(2)*T(2))*(T(3)*T(6)+T(6)*T(3))*IN(31)+(T(5)*T(5))*(T(3)*T(3))*IN(2)+(T(5)*T(5))*(T(6)*T(6))*IN(8)+(T(5)*T(5))*(T(9)*T(9))*IN(14)+(T(5)*T(5))*(T(6)*T(9)+T(9)*T(6))*IN(20)+(T(5)*T(5))*(T(9)*T(3)+T(3)*T(9))*IN(26)+(T(5)*T(5))*(T(3)*T(6)+T(6)*T(3))*IN(32)+(T(8)*T(8))*(T(3)*T(3))*IN(3)+(T(8)*T(8))*(T(6)*T(6))*IN(9)+(T(8)*T(8))*(T(9)*T(9))*IN(15)+(T(8)*T(8))*(T(6)*T(9)+T(9)*T(6))*IN(21)+(T(8)*T(8))*(T(9)*T(3)+T(3)*T(9))*IN(27)+(T(8)*T(8))*(T(3)*T(6)+T(6)*T(3))*IN(33)+(T(5)*T(8)+T(8)*T(5))*(T(3)*T(3))*IN(4)+(T(5)*T(8)+T(8)*T(5))*(T(6)*T(6))*IN(10)+(T(5)*T(8)+T(8)*T(5))*(T(9)*T(9))*IN(16)+(T(5)*T(8)+T(8)*T(5))*(T(6)*T(9)+T(9)*T(6))*IN(22)+(T(5)*T(8)+T(8)*T(5))*(T(9)*T(3)+T(3)*T(9))*IN(28)+(T(5)*T(8)+T(8)*T(5))*(T(3)*T(6)+T(6)*T(3))*IN(34)+(T(8)*T(2)+T(2)*T(8))*(T(3)*T(3))*IN(5)+(T(8)*T(2)+T(2)*T(8))*(T(6)*T(6))*IN(11)+(T(8)*T(2)+T(2)*T(8))*(T(9)*T(9))*IN(17)+(T(8)*T(2)+T(2)*T(8))*(T(6)*T(9)+T(9)*T(6))*IN(23)+(T(8)*T(2)+T(2)*T(8))*(T(9)*T(3)+T(3)*T(9))*IN(29)+(T(8)*T(2)+T(2)*T(8))*(T(3)*T(6)+T(6)*T(3))*IN(35)+(T(2)*T(5)+T(5)*T(2))*(T(3)*T(3))*IN(6)+(T(2)*T(5)+T(5)*T(2))*(T(6)*T(6))*IN(12)+(T(2)*T(5)+T(5)*T(2))*(T(9)*T(9))*IN(18)+(T(2)*T(5)+T(5)*T(2))*(T(6)*T(9)+T(9)*T(6))*IN(24)+(T(2)*T(5)+T(5)*T(2))*(T(9)*T(3)+T(3)*T(9))*IN(30)+(T(2)*T(5)+T(5)*T(2))*(T(3)*T(6)+T(6)*T(3))*IN(36));
S(20)=((T(2)*T(2))*(T(2)*T(3))*IN(1)+(T(2)*T(2))*(T(5)*T(6))*IN(7)+(T(2)*T(2))*(T(8)*T(9))*IN(13)+(T(2)*T(2))*(T(5)*T(9)+T(8)*T(6))*IN(19)+(T(2)*T(2))*(T(8)*T(3)+T(2)*T(9))*IN(25)+(T(2)*T(2))*(T(2)*T(6)+T(5)*T(3))*IN(31)+(T(5)*T(5))*(T(2)*T(3))*IN(2)+(T(5)*T(5))*(T(5)*T(6))*IN(8)+(T(5)*T(5))*(T(8)*T(9))*IN(14)+(T(5)*T(5))*(T(5)*T(9)+T(8)*T(6))*IN(20)+(T(5)*T(5))*(T(8)*T(3)+T(2)*T(9))*IN(26)+(T(5)*T(5))*(T(2)*T(6)+T(5)*T(3))*IN(32)+(T(8)*T(8))*(T(2)*T(3))*IN(3)+(T(8)*T(8))*(T(5)*T(6))*IN(9)+(T(8)*T(8))*(T(8)*T(9))*IN(15)+(T(8)*T(8))*(T(5)*T(9)+T(8)*T(6))*IN(21)+(T(8)*T(8))*(T(8)*T(3)+T(2)*T(9))*IN(27)+(T(8)*T(8))*(T(2)*T(6)+T(5)*T(3))*IN(33)+(T(5)*T(8)+T(8)*T(5))*(T(2)*T(3))*IN(4)+(T(5)*T(8)+T(8)*T(5))*(T(5)*T(6))*IN(10)+(T(5)*T(8)+T(8)*T(5))*(T(8)*T(9))*IN(16)+(T(5)*T(8)+T(8)*T(5))*(T(5)*T(9)+T(8)*T(6))*IN(22)+(T(5)*T(8)+T(8)*T(5))*(T(8)*T(3)+T(2)*T(9))*IN(28)+(T(5)*T(8)+T(8)*T(5))*(T(2)*T(6)+T(5)*T(3))*IN(34)+(T(8)*T(2)+T(2)*T(8))*(T(2)*T(3))*IN(5)+(T(8)*T(2)+T(2)*T(8))*(T(5)*T(6))*IN(11)+(T(8)*T(2)+T(2)*T(8))*(T(8)*T(9))*IN(17)+(T(8)*T(2)+T(2)*T(8))*(T(5)*T(9)+T(8)*T(6))*IN(23)+(T(8)*T(2)+T(2)*T(8))*(T(8)*T(3)+T(2)*T(9))*IN(29)+(T(8)*T(2)+T(2)*T(8))*(T(2)*T(6)+T(5)*T(3))*IN(35)+(T(2)*T(5)+T(5)*T(2))*(T(2)*T(3))*IN(6)+(T(2)*T(5)+T(5)*T(2))*(T(5)*T(6))*IN(12)+(T(2)*T(5)+T(5)*T(2))*(T(8)*T(9))*IN(18)+(T(2)*T(5)+T(5)*T(2))*(T(5)*T(9)+T(8)*T(6))*IN(24)+(T(2)*T(5)+T(5)*T(2))*(T(8)*T(3)+T(2)*T(9))*IN(30)+(T(2)*T(5)+T(5)*T(2))*(T(2)*T(6)+T(5)*T(3))*IN(36));
S(26)=((T(2)*T(2))*(T(3)*T(1))*IN(1)+(T(2)*T(2))*(T(6)*T(4))*IN(7)+(T(2)*T(2))*(T(9)*T(7))*IN(13)+(T(2)*T(2))*(T(6)*T(7)+T(9)*T(4))*IN(19)+(T(2)*T(2))*(T(9)*T(1)+T(3)*T(7))*IN(25)+(T(2)*T(2))*(T(3)*T(4)+T(6)*T(1))*IN(31)+(T(5)*T(5))*(T(3)*T(1))*IN(2)+(T(5)*T(5))*(T(6)*T(4))*IN(8)+(T(5)*T(5))*(T(9)*T(7))*IN(14)+(T(5)*T(5))*(T(6)*T(7)+T(9)*T(4))*IN(20)+(T(5)*T(5))*(T(9)*T(1)+T(3)*T(7))*IN(26)+(T(5)*T(5))*(T(3)*T(4)+T(6)*T(1))*IN(32)+(T(8)*T(8))*(T(3)*T(1))*IN(3)+(T(8)*T(8))*(T(6)*T(4))*IN(9)+(T(8)*T(8))*(T(9)*T(7))*IN(15)+(T(8)*T(8))*(T(6)*T(7)+T(9)*T(4))*IN(21)+(T(8)*T(8))*(T(9)*T(1)+T(3)*T(7))*IN(27)+(T(8)*T(8))*(T(3)*T(4)+T(6)*T(1))*IN(33)+(T(5)*T(8)+T(8)*T(5))*(T(3)*T(1))*IN(4)+(T(5)*T(8)+T(8)*T(5))*(T(6)*T(4))*IN(10)+(T(5)*T(8)+T(8)*T(5))*(T(9)*T(7))*IN(16)+(T(5)*T(8)+T(8)*T(5))*(T(6)*T(7)+T(9)*T(4))*IN(22)+(T(5)*T(8)+T(8)*T(5))*(T(9)*T(1)+T(3)*T(7))*IN(28)+(T(5)*T(8)+T(8)*T(5))*(T(3)*T(4)+T(6)*T(1))*IN(34)+(T(8)*T(2)+T(2)*T(8))*(T(3)*T(1))*IN(5)+(T(8)*T(2)+T(2)*T(8))*(T(6)*T(4))*IN(11)+(T(8)*T(2)+T(2)*T(8))*(T(9)*T(7))*IN(17)+(T(8)*T(2)+T(2)*T(8))*(T(6)*T(7)+T(9)*T(4))*IN(23)+(T(8)*T(2)+T(2)*T(8))*(T(9)*T(1)+T(3)*T(7))*IN(29)+(T(8)*T(2)+T(2)*T(8))*(T(3)*T(4)+T(6)*T(1))*IN(35)+(T(2)*T(5)+T(5)*T(2))*(T(3)*T(1))*IN(6)+(T(2)*T(5)+T(5)*T(2))*(T(6)*T(4))*IN(12)+(T(2)*T(5)+T(5)*T(2))*(T(9)*T(7))*IN(18)+(T(2)*T(5)+T(5)*T(2))*(T(6)*T(7)+T(9)*T(4))*IN(24)+(T(2)*T(5)+T(5)*T(2))*(T(9)*T(1)+T(3)*T(7))*IN(30)+(T(2)*T(5)+T(5)*T(2))*(T(3)*T(4)+T(6)*T(1))*IN(36));
S(32)=((T(2)*T(2))*(T(1)*T(2))*IN(1)+(T(2)*T(2))*(T(4)*T(5))*IN(7)+(T(2)*T(2))*(T(7)*T(8))*IN(13)+(T(2)*T(2))*(T(4)*T(8)+T(7)*T(5))*IN(19)+(T(2)*T(2))*(T(7)*T(2)+T(1)*T(8))*IN(25)+(T(2)*T(2))*(T(1)*T(5)+T(4)*T(2))*IN(31)+(T(5)*T(5))*(T(1)*T(2))*IN(2)+(T(5)*T(5))*(T(4)*T(5))*IN(8)+(T(5)*T(5))*(T(7)*T(8))*IN(14)+(T(5)*T(5))*(T(4)*T(8)+T(7)*T(5))*IN(20)+(T(5)*T(5))*(T(7)*T(2)+T(1)*T(8))*IN(26)+(T(5)*T(5))*(T(1)*T(5)+T(4)*T(2))*IN(32)+(T(8)*T(8))*(T(1)*T(2))*IN(3)+(T(8)*T(8))*(T(4)*T(5))*IN(9)+(T(8)*T(8))*(T(7)*T(8))*IN(15)+(T(8)*T(8))*(T(4)*T(8)+T(7)*T(5))*IN(21)+(T(8)*T(8))*(T(7)*T(2)+T(1)*T(8))*IN(27)+(T(8)*T(8))*(T(1)*T(5)+T(4)*T(2))*IN(33)+(T(5)*T(8)+T(8)*T(5))*(T(1)*T(2))*IN(4)+(T(5)*T(8)+T(8)*T(5))*(T(4)*T(5))*IN(10)+(T(5)*T(8)+T(8)*T(5))*(T(7)*T(8))*IN(16)+(T(5)*T(8)+T(8)*T(5))*(T(4)*T(8)+T(7)*T(5))*IN(22)+(T(5)*T(8)+T(8)*T(5))*(T(7)*T(2)+T(1)*T(8))*IN(28)+(T(5)*T(8)+T(8)*T(5))*(T(1)*T(5)+T(4)*T(2))*IN(34)+(T(8)*T(2)+T(2)*T(8))*(T(1)*T(2))*IN(5)+(T(8)*T(2)+T(2)*T(8))*(T(4)*T(5))*IN(11)+(T(8)*T(2)+T(2)*T(8))*(T(7)*T(8))*IN(17)+(T(8)*T(2)+T(2)*T(8))*(T(4)*T(8)+T(7)*T(5))*IN(23)+(T(8)*T(2)+T(2)*T(8))*(T(7)*T(2)+T(1)*T(8))*IN(29)+(T(8)*T(2)+T(2)*T(8))*(T(1)*T(5)+T(4)*T(2))*IN(35)+(T(2)*T(5)+T(5)*T(2))*(T(1)*T(2))*IN(6)+(T(2)*T(5)+T(5)*T(2))*(T(4)*T(5))*IN(12)+(T(2)*T(5)+T(5)*T(2))*(T(7)*T(8))*IN(18)+(T(2)*T(5)+T(5)*T(2))*(T(4)*T(8)+T(7)*T(5))*IN(24)+(T(2)*T(5)+T(5)*T(2))*(T(7)*T(2)+T(1)*T(8))*IN(30)+(T(2)*T(5)+T(5)*T(2))*(T(1)*T(5)+T(4)*T(2))*IN(36));
S(3)=((T(3)*T(3))*(T(1)*T(1))*IN(1)+(T(3)*T(3))*(T(4)*T(4))*IN(7)+(T(3)*T(3))*(T(7)*T(7))*IN(13)+(T(3)*T(3))*(T(4)*T(7)+T(7)*T(4))*IN(19)+(T(3)*T(3))*(T(7)*T(1)+T(1)*T(7))*IN(25)+(T(3)*T(3))*(T(1)*T(4)+T(4)*T(1))*IN(31)+(T(6)*T(6))*(T(1)*T(1))*IN(2)+(T(6)*T(6))*(T(4)*T(4))*IN(8)+(T(6)*T(6))*(T(7)*T(7))*IN(14)+(T(6)*T(6))*(T(4)*T(7)+T(7)*T(4))*IN(20)+(T(6)*T(6))*(T(7)*T(1)+T(1)*T(7))*IN(26)+(T(6)*T(6))*(T(1)*T(4)+T(4)*T(1))*IN(32)+(T(9)*T(9))*(T(1)*T(1))*IN(3)+(T(9)*T(9))*(T(4)*T(4))*IN(9)+(T(9)*T(9))*(T(7)*T(7))*IN(15)+(T(9)*T(9))*(T(4)*T(7)+T(7)*T(4))*IN(21)+(T(9)*T(9))*(T(7)*T(1)+T(1)*T(7))*IN(27)+(T(9)*T(9))*(T(1)*T(4)+T(4)*T(1))*IN(33)+(T(6)*T(9)+T(9)*T(6))*(T(1)*T(1))*IN(4)+(T(6)*T(9)+T(9)*T(6))*(T(4)*T(4))*IN(10)+(T(6)*T(9)+T(9)*T(6))*(T(7)*T(7))*IN(16)+(T(6)*T(9)+T(9)*T(6))*(T(4)*T(7)+T(7)*T(4))*IN(22)+(T(6)*T(9)+T(9)*T(6))*(T(7)*T(1)+T(1)*T(7))*IN(28)+(T(6)*T(9)+T(9)*T(6))*(T(1)*T(4)+T(4)*T(1))*IN(34)+(T(9)*T(3)+T(3)*T(9))*(T(1)*T(1))*IN(5)+(T(9)*T(3)+T(3)*T(9))*(T(4)*T(4))*IN(11)+(T(9)*T(3)+T(3)*T(9))*(T(7)*T(7))*IN(17)+(T(9)*T(3)+T(3)*T(9))*(T(4)*T(7)+T(7)*T(4))*IN(23)+(T(9)*T(3)+T(3)*T(9))*(T(7)*T(1)+T(1)*T(7))*IN(29)+(T(9)*T(3)+T(3)*T(9))*(T(1)*T(4)+T(4)*T(1))*IN(35)+(T(3)*T(6)+T(6)*T(3))*(T(1)*T(1))*IN(6)+(T(3)*T(6)+T(6)*T(3))*(T(4)*T(4))*IN(12)+(T(3)*T(6)+T(6)*T(3))*(T(7)*T(7))*IN(18)+(T(3)*T(6)+T(6)*T(3))*(T(4)*T(7)+T(7)*T(4))*IN(24)+(T(3)*T(6)+T(6)*T(3))*(T(7)*T(1)+T(1)*T(7))*IN(30)+(T(3)*T(6)+T(6)*T(3))*(T(1)*T(4)+T(4)*T(1))*IN(36));
S(9)=((T(3)*T(3))*(T(2)*T(2))*IN(1)+(T(3)*T(3))*(T(5)*T(5))*IN(7)+(T(3)*T(3))*(T(8)*T(8))*IN(13)+(T(3)*T(3))*(T(5)*T(8)+T(8)*T(5))*IN(19)+(T(3)*T(3))*(T(8)*T(2)+T(2)*T(8))*IN(25)+(T(3)*T(3))*(T(2)*T(5)+T(5)*T(2))*IN(31)+(T(6)*T(6))*(T(2)*T(2))*IN(2)+(T(6)*T(6))*(T(5)*T(5))*IN(8)+(T(6)*T(6))*(T(8)*T(8))*IN(14)+(T(6)*T(6))*(T(5)*T(8)+T(8)*T(5))*IN(20)+(T(6)*T(6))*(T(8)*T(2)+T(2)*T(8))*IN(26)+(T(6)*T(6))*(T(2)*T(5)+T(5)*T(2))*IN(32)+(T(9)*T(9))*(T(2)*T(2))*IN(3)+(T(9)*T(9))*(T(5)*T(5))*IN(9)+(T(9)*T(9))*(T(8)*T(8))*IN(15)+(T(9)*T(9))*(T(5)*T(8)+T(8)*T(5))*IN(21)+(T(9)*T(9))*(T(8)*T(2)+T(2)*T(8))*IN(27)+(T(9)*T(9))*(T(2)*T(5)+T(5)*T(2))*IN(33)+(T(6)*T(9)+T(9)*T(6))*(T(2)*T(2))*IN(4)+(T(6)*T(9)+T(9)*T(6))*(T(5)*T(5))*IN(10)+(T(6)*T(9)+T(9)*T(6))*(T(8)*T(8))*IN(16)+(T(6)*T(9)+T(9)*T(6))*(T(5)*T(8)+T(8)*T(5))*IN(22)+(T(6)*T(9)+T(9)*T(6))*(T(8)*T(2)+T(2)*T(8))*IN(28)+(T(6)*T(9)+T(9)*T(6))*(T(2)*T(5)+T(5)*T(2))*IN(34)+(T(9)*T(3)+T(3)*T(9))*(T(2)*T(2))*IN(5)+(T(9)*T(3)+T(3)*T(9))*(T(5)*T(5))*IN(11)+(T(9)*T(3)+T(3)*T(9))*(T(8)*T(8))*IN(17)+(T(9)*T(3)+T(3)*T(9))*(T(5)*T(8)+T(8)*T(5))*IN(23)+(T(9)*T(3)+T(3)*T(9))*(T(8)*T(2)+T(2)*T(8))*IN(29)+(T(9)*T(3)+T(3)*T(9))*(T(2)*T(5)+T(5)*T(2))*IN(35)+(T(3)*T(6)+T(6)*T(3))*(T(2)*T(2))*IN(6)+(T(3)*T(6)+T(6)*T(3))*(T(5)*T(5))*IN(12)+(T(3)*T(6)+T(6)*T(3))*(T(8)*T(8))*IN(18)+(T(3)*T(6)+T(6)*T(3))*(T(5)*T(8)+T(8)*T(5))*IN(24)+(T(3)*T(6)+T(6)*T(3))*(T(8)*T(2)+T(2)*T(8))*IN(30)+(T(3)*T(6)+T(6)*T(3))*(T(2)*T(5)+T(5)*T(2))*IN(36));
S(15)=((T(3)*T(3))*(T(3)*T(3))*IN(1)+(T(3)*T(3))*(T(6)*T(6))*IN(7)+(T(3)*T(3))*(T(9)*T(9))*IN(13)+(T(3)*T(3))*(T(6)*T(9)+T(9)*T(6))*IN(19)+(T(3)*T(3))*(T(9)*T(3)+T(3)*T(9))*IN(25)+(T(3)*T(3))*(T(3)*T(6)+T(6)*T(3))*IN(31)+(T(6)*T(6))*(T(3)*T(3))*IN(2)+(T(6)*T(6))*(T(6)*T(6))*IN(8)+(T(6)*T(6))*(T(9)*T(9))*IN(14)+(T(6)*T(6))*(T(6)*T(9)+T(9)*T(6))*IN(20)+(T(6)*T(6))*(T(9)*T(3)+T(3)*T(9))*IN(26)+(T(6)*T(6))*(T(3)*T(6)+T(6)*T(3))*IN(32)+(T(9)*T(9))*(T(3)*T(3))*IN(3)+(T(9)*T(9))*(T(6)*T(6))*IN(9)+(T(9)*T(9))*(T(9)*T(9))*IN(15)+(T(9)*T(9))*(T(6)*T(9)+T(9)*T(6))*IN(21)+(T(9)*T(9))*(T(9)*T(3)+T(3)*T(9))*IN(27)+(T(9)*T(9))*(T(3)*T(6)+T(6)*T(3))*IN(33)+(T(6)*T(9)+T(9)*T(6))*(T(3)*T(3))*IN(4)+(T(6)*T(9)+T(9)*T(6))*(T(6)*T(6))*IN(10)+(T(6)*T(9)+T(9)*T(6))*(T(9)*T(9))*IN(16)+(T(6)*T(9)+T(9)*T(6))*(T(6)*T(9)+T(9)*T(6))*IN(22)+(T(6)*T(9)+T(9)*T(6))*(T(9)*T(3)+T(3)*T(9))*IN(28)+(T(6)*T(9)+T(9)*T(6))*(T(3)*T(6)+T(6)*T(3))*IN(34)+(T(9)*T(3)+T(3)*T(9))*(T(3)*T(3))*IN(5)+(T(9)*T(3)+T(3)*T(9))*(T(6)*T(6))*IN(11)+(T(9)*T(3)+T(3)*T(9))*(T(9)*T(9))*IN(17)+(T(9)*T(3)+T(3)*T(9))*(T(6)*T(9)+T(9)*T(6))*IN(23)+(T(9)*T(3)+T(3)*T(9))*(T(9)*T(3)+T(3)*T(9))*IN(29)+(T(9)*T(3)+T(3)*T(9))*(T(3)*T(6)+T(6)*T(3))*IN(35)+(T(3)*T(6)+T(6)*T(3))*(T(3)*T(3))*IN(6)+(T(3)*T(6)+T(6)*T(3))*(T(6)*T(6))*IN(12)+(T(3)*T(6)+T(6)*T(3))*(T(9)*T(9))*IN(18)+(T(3)*T(6)+T(6)*T(3))*(T(6)*T(9)+T(9)*T(6))*IN(24)+(T(3)*T(6)+T(6)*T(3))*(T(9)*T(3)+T(3)*T(9))*IN(30)+(T(3)*T(6)+T(6)*T(3))*(T(3)*T(6)+T(6)*T(3))*IN(36));
S(21)=((T(3)*T(3))*(T(2)*T(3))*IN(1)+(T(3)*T(3))*(T(5)*T(6))*IN(7)+(T(3)*T(3))*(T(8)*T(9))*IN(13)+(T(3)*T(3))*(T(5)*T(9)+T(8)*T(6))*IN(19)+(T(3)*T(3))*(T(8)*T(3)+T(2)*T(9))*IN(25)+(T(3)*T(3))*(T(2)*T(6)+T(5)*T(3))*IN(31)+(T(6)*T(6))*(T(2)*T(3))*IN(2)+(T(6)*T(6))*(T(5)*T(6))*IN(8)+(T(6)*T(6))*(T(8)*T(9))*IN(14)+(T(6)*T(6))*(T(5)*T(9)+T(8)*T(6))*IN(20)+(T(6)*T(6))*(T(8)*T(3)+T(2)*T(9))*IN(26)+(T(6)*T(6))*(T(2)*T(6)+T(5)*T(3))*IN(32)+(T(9)*T(9))*(T(2)*T(3))*IN(3)+(T(9)*T(9))*(T(5)*T(6))*IN(9)+(T(9)*T(9))*(T(8)*T(9))*IN(15)+(T(9)*T(9))*(T(5)*T(9)+T(8)*T(6))*IN(21)+(T(9)*T(9))*(T(8)*T(3)+T(2)*T(9))*IN(27)+(T(9)*T(9))*(T(2)*T(6)+T(5)*T(3))*IN(33)+(T(6)*T(9)+T(9)*T(6))*(T(2)*T(3))*IN(4)+(T(6)*T(9)+T(9)*T(6))*(T(5)*T(6))*IN(10)+(T(6)*T(9)+T(9)*T(6))*(T(8)*T(9))*IN(16)+(T(6)*T(9)+T(9)*T(6))*(T(5)*T(9)+T(8)*T(6))*IN(22)+(T(6)*T(9)+T(9)*T(6))*(T(8)*T(3)+T(2)*T(9))*IN(28)+(T(6)*T(9)+T(9)*T(6))*(T(2)*T(6)+T(5)*T(3))*IN(34)+(T(9)*T(3)+T(3)*T(9))*(T(2)*T(3))*IN(5)+(T(9)*T(3)+T(3)*T(9))*(T(5)*T(6))*IN(11)+(T(9)*T(3)+T(3)*T(9))*(T(8)*T(9))*IN(17)+(T(9)*T(3)+T(3)*T(9))*(T(5)*T(9)+T(8)*T(6))*IN(23)+(T(9)*T(3)+T(3)*T(9))*(T(8)*T(3)+T(2)*T(9))*IN(29)+(T(9)*T(3)+T(3)*T(9))*(T(2)*T(6)+T(5)*T(3))*IN(35)+(T(3)*T(6)+T(6)*T(3))*(T(2)*T(3))*IN(6)+(T(3)*T(6)+T(6)*T(3))*(T(5)*T(6))*IN(12)+(T(3)*T(6)+T(6)*T(3))*(T(8)*T(9))*IN(18)+(T(3)*T(6)+T(6)*T(3))*(T(5)*T(9)+T(8)*T(6))*IN(24)+(T(3)*T(6)+T(6)*T(3))*(T(8)*T(3)+T(2)*T(9))*IN(30)+(T(3)*T(6)+T(6)*T(3))*(T(2)*T(6)+T(5)*T(3))*IN(36));
S(27)=((T(3)*T(3))*(T(3)*T(1))*IN(1)+(T(3)*T(3))*(T(6)*T(4))*IN(7)+(T(3)*T(3))*(T(9)*T(7))*IN(13)+(T(3)*T(3))*(T(6)*T(7)+T(9)*T(4))*IN(19)+(T(3)*T(3))*(T(9)*T(1)+T(3)*T(7))*IN(25)+(T(3)*T(3))*(T(3)*T(4)+T(6)*T(1))*IN(31)+(T(6)*T(6))*(T(3)*T(1))*IN(2)+(T(6)*T(6))*(T(6)*T(4))*IN(8)+(T(6)*T(6))*(T(9)*T(7))*IN(14)+(T(6)*T(6))*(T(6)*T(7)+T(9)*T(4))*IN(20)+(T(6)*T(6))*(T(9)*T(1)+T(3)*T(7))*IN(26)+(T(6)*T(6))*(T(3)*T(4)+T(6)*T(1))*IN(32)+(T(9)*T(9))*(T(3)*T(1))*IN(3)+(T(9)*T(9))*(T(6)*T(4))*IN(9)+(T(9)*T(9))*(T(9)*T(7))*IN(15)+(T(9)*T(9))*(T(6)*T(7)+T(9)*T(4))*IN(21)+(T(9)*T(9))*(T(9)*T(1)+T(3)*T(7))*IN(27)+(T(9)*T(9))*(T(3)*T(4)+T(6)*T(1))*IN(33)+(T(6)*T(9)+T(9)*T(6))*(T(3)*T(1))*IN(4)+(T(6)*T(9)+T(9)*T(6))*(T(6)*T(4))*IN(10)+(T(6)*T(9)+T(9)*T(6))*(T(9)*T(7))*IN(16)+(T(6)*T(9)+T(9)*T(6))*(T(6)*T(7)+T(9)*T(4))*IN(22)+(T(6)*T(9)+T(9)*T(6))*(T(9)*T(1)+T(3)*T(7))*IN(28)+(T(6)*T(9)+T(9)*T(6))*(T(3)*T(4)+T(6)*T(1))*IN(34)+(T(9)*T(3)+T(3)*T(9))*(T(3)*T(1))*IN(5)+(T(9)*T(3)+T(3)*T(9))*(T(6)*T(4))*IN(11)+(T(9)*T(3)+T(3)*T(9))*(T(9)*T(7))*IN(17)+(T(9)*T(3)+T(3)*T(9))*(T(6)*T(7)+T(9)*T(4))*IN(23)+(T(9)*T(3)+T(3)*T(9))*(T(9)*T(1)+T(3)*T(7))*IN(29)+(T(9)*T(3)+T(3)*T(9))*(T(3)*T(4)+T(6)*T(1))*IN(35)+(T(3)*T(6)+T(6)*T(3))*(T(3)*T(1))*IN(6)+(T(3)*T(6)+T(6)*T(3))*(T(6)*T(4))*IN(12)+(T(3)*T(6)+T(6)*T(3))*(T(9)*T(7))*IN(18)+(T(3)*T(6)+T(6)*T(3))*(T(6)*T(7)+T(9)*T(4))*IN(24)+(T(3)*T(6)+T(6)*T(3))*(T(9)*T(1)+T(3)*T(7))*IN(30)+(T(3)*T(6)+T(6)*T(3))*(T(3)*T(4)+T(6)*T(1))*IN(36));
S(33)=((T(3)*T(3))*(T(1)*T(2))*IN(1)+(T(3)*T(3))*(T(4)*T(5))*IN(7)+(T(3)*T(3))*(T(7)*T(8))*IN(13)+(T(3)*T(3))*(T(4)*T(8)+T(7)*T(5))*IN(19)+(T(3)*T(3))*(T(7)*T(2)+T(1)*T(8))*IN(25)+(T(3)*T(3))*(T(1)*T(5)+T(4)*T(2))*IN(31)+(T(6)*T(6))*(T(1)*T(2))*IN(2)+(T(6)*T(6))*(T(4)*T(5))*IN(8)+(T(6)*T(6))*(T(7)*T(8))*IN(14)+(T(6)*T(6))*(T(4)*T(8)+T(7)*T(5))*IN(20)+(T(6)*T(6))*(T(7)*T(2)+T(1)*T(8))*IN(26)+(T(6)*T(6))*(T(1)*T(5)+T(4)*T(2))*IN(32)+(T(9)*T(9))*(T(1)*T(2))*IN(3)+(T(9)*T(9))*(T(4)*T(5))*IN(9)+(T(9)*T(9))*(T(7)*T(8))*IN(15)+(T(9)*T(9))*(T(4)*T(8)+T(7)*T(5))*IN(21)+(T(9)*T(9))*(T(7)*T(2)+T(1)*T(8))*IN(27)+(T(9)*T(9))*(T(1)*T(5)+T(4)*T(2))*IN(33)+(T(6)*T(9)+T(9)*T(6))*(T(1)*T(2))*IN(4)+(T(6)*T(9)+T(9)*T(6))*(T(4)*T(5))*IN(10)+(T(6)*T(9)+T(9)*T(6))*(T(7)*T(8))*IN(16)+(T(6)*T(9)+T(9)*T(6))*(T(4)*T(8)+T(7)*T(5))*IN(22)+(T(6)*T(9)+T(9)*T(6))*(T(7)*T(2)+T(1)*T(8))*IN(28)+(T(6)*T(9)+T(9)*T(6))*(T(1)*T(5)+T(4)*T(2))*IN(34)+(T(9)*T(3)+T(3)*T(9))*(T(1)*T(2))*IN(5)+(T(9)*T(3)+T(3)*T(9))*(T(4)*T(5))*IN(11)+(T(9)*T(3)+T(3)*T(9))*(T(7)*T(8))*IN(17)+(T(9)*T(3)+T(3)*T(9))*(T(4)*T(8)+T(7)*T(5))*IN(23)+(T(9)*T(3)+T(3)*T(9))*(T(7)*T(2)+T(1)*T(8))*IN(29)+(T(9)*T(3)+T(3)*T(9))*(T(1)*T(5)+T(4)*T(2))*IN(35)+(T(3)*T(6)+T(6)*T(3))*(T(1)*T(2))*IN(6)+(T(3)*T(6)+T(6)*T(3))*(T(4)*T(5))*IN(12)+(T(3)*T(6)+T(6)*T(3))*(T(7)*T(8))*IN(18)+(T(3)*T(6)+T(6)*T(3))*(T(4)*T(8)+T(7)*T(5))*IN(24)+(T(3)*T(6)+T(6)*T(3))*(T(7)*T(2)+T(1)*T(8))*IN(30)+(T(3)*T(6)+T(6)*T(3))*(T(1)*T(5)+T(4)*T(2))*IN(36));
S(4)=((T(2)*T(3))*(T(1)*T(1))*IN(1)+(T(2)*T(3))*(T(4)*T(4))*IN(7)+(T(2)*T(3))*(T(7)*T(7))*IN(13)+(T(2)*T(3))*(T(4)*T(7)+T(7)*T(4))*IN(19)+(T(2)*T(3))*(T(7)*T(1)+T(1)*T(7))*IN(25)+(T(2)*T(3))*(T(1)*T(4)+T(4)*T(1))*IN(31)+(T(5)*T(6))*(T(1)*T(1))*IN(2)+(T(5)*T(6))*(T(4)*T(4))*IN(8)+(T(5)*T(6))*(T(7)*T(7))*IN(14)+(T(5)*T(6))*(T(4)*T(7)+T(7)*T(4))*IN(20)+(T(5)*T(6))*(T(7)*T(1)+T(1)*T(7))*IN(26)+(T(5)*T(6))*(T(1)*T(4)+T(4)*T(1))*IN(32)+(T(8)*T(9))*(T(1)*T(1))*IN(3)+(T(8)*T(9))*(T(4)*T(4))*IN(9)+(T(8)*T(9))*(T(7)*T(7))*IN(15)+(T(8)*T(9))*(T(4)*T(7)+T(7)*T(4))*IN(21)+(T(8)*T(9))*(T(7)*T(1)+T(1)*T(7))*IN(27)+(T(8)*T(9))*(T(1)*T(4)+T(4)*T(1))*IN(33)+(T(5)*T(9)+T(8)*T(6))*(T(1)*T(1))*IN(4)+(T(5)*T(9)+T(8)*T(6))*(T(4)*T(4))*IN(10)+(T(5)*T(9)+T(8)*T(6))*(T(7)*T(7))*IN(16)+(T(5)*T(9)+T(8)*T(6))*(T(4)*T(7)+T(7)*T(4))*IN(22)+(T(5)*T(9)+T(8)*T(6))*(T(7)*T(1)+T(1)*T(7))*IN(28)+(T(5)*T(9)+T(8)*T(6))*(T(1)*T(4)+T(4)*T(1))*IN(34)+(T(8)*T(3)+T(2)*T(9))*(T(1)*T(1))*IN(5)+(T(8)*T(3)+T(2)*T(9))*(T(4)*T(4))*IN(11)+(T(8)*T(3)+T(2)*T(9))*(T(7)*T(7))*IN(17)+(T(8)*T(3)+T(2)*T(9))*(T(4)*T(7)+T(7)*T(4))*IN(23)+(T(8)*T(3)+T(2)*T(9))*(T(7)*T(1)+T(1)*T(7))*IN(29)+(T(8)*T(3)+T(2)*T(9))*(T(1)*T(4)+T(4)*T(1))*IN(35)+(T(2)*T(6)+T(5)*T(3))*(T(1)*T(1))*IN(6)+(T(2)*T(6)+T(5)*T(3))*(T(4)*T(4))*IN(12)+(T(2)*T(6)+T(5)*T(3))*(T(7)*T(7))*IN(18)+(T(2)*T(6)+T(5)*T(3))*(T(4)*T(7)+T(7)*T(4))*IN(24)+(T(2)*T(6)+T(5)*T(3))*(T(7)*T(1)+T(1)*T(7))*IN(30)+(T(2)*T(6)+T(5)*T(3))*(T(1)*T(4)+T(4)*T(1))*IN(36));
S(10)=((T(2)*T(3))*(T(2)*T(2))*IN(1)+(T(2)*T(3))*(T(5)*T(5))*IN(7)+(T(2)*T(3))*(T(8)*T(8))*IN(13)+(T(2)*T(3))*(T(5)*T(8)+T(8)*T(5))*IN(19)+(T(2)*T(3))*(T(8)*T(2)+T(2)*T(8))*IN(25)+(T(2)*T(3))*(T(2)*T(5)+T(5)*T(2))*IN(31)+(T(5)*T(6))*(T(2)*T(2))*IN(2)+(T(5)*T(6))*(T(5)*T(5))*IN(8)+(T(5)*T(6))*(T(8)*T(8))*IN(14)+(T(5)*T(6))*(T(5)*T(8)+T(8)*T(5))*IN(20)+(T(5)*T(6))*(T(8)*T(2)+T(2)*T(8))*IN(26)+(T(5)*T(6))*(T(2)*T(5)+T(5)*T(2))*IN(32)+(T(8)*T(9))*(T(2)*T(2))*IN(3)+(T(8)*T(9))*(T(5)*T(5))*IN(9)+(T(8)*T(9))*(T(8)*T(8))*IN(15)+(T(8)*T(9))*(T(5)*T(8)+T(8)*T(5))*IN(21)+(T(8)*T(9))*(T(8)*T(2)+T(2)*T(8))*IN(27)+(T(8)*T(9))*(T(2)*T(5)+T(5)*T(2))*IN(33)+(T(5)*T(9)+T(8)*T(6))*(T(2)*T(2))*IN(4)+(T(5)*T(9)+T(8)*T(6))*(T(5)*T(5))*IN(10)+(T(5)*T(9)+T(8)*T(6))*(T(8)*T(8))*IN(16)+(T(5)*T(9)+T(8)*T(6))*(T(5)*T(8)+T(8)*T(5))*IN(22)+(T(5)*T(9)+T(8)*T(6))*(T(8)*T(2)+T(2)*T(8))*IN(28)+(T(5)*T(9)+T(8)*T(6))*(T(2)*T(5)+T(5)*T(2))*IN(34)+(T(8)*T(3)+T(2)*T(9))*(T(2)*T(2))*IN(5)+(T(8)*T(3)+T(2)*T(9))*(T(5)*T(5))*IN(11)+(T(8)*T(3)+T(2)*T(9))*(T(8)*T(8))*IN(17)+(T(8)*T(3)+T(2)*T(9))*(T(5)*T(8)+T(8)*T(5))*IN(23)+(T(8)*T(3)+T(2)*T(9))*(T(8)*T(2)+T(2)*T(8))*IN(29)+(T(8)*T(3)+T(2)*T(9))*(T(2)*T(5)+T(5)*T(2))*IN(35)+(T(2)*T(6)+T(5)*T(3))*(T(2)*T(2))*IN(6)+(T(2)*T(6)+T(5)*T(3))*(T(5)*T(5))*IN(12)+(T(2)*T(6)+T(5)*T(3))*(T(8)*T(8))*IN(18)+(T(2)*T(6)+T(5)*T(3))*(T(5)*T(8)+T(8)*T(5))*IN(24)+(T(2)*T(6)+T(5)*T(3))*(T(8)*T(2)+T(2)*T(8))*IN(30)+(T(2)*T(6)+T(5)*T(3))*(T(2)*T(5)+T(5)*T(2))*IN(36));
S(16)=((T(2)*T(3))*(T(3)*T(3))*IN(1)+(T(2)*T(3))*(T(6)*T(6))*IN(7)+(T(2)*T(3))*(T(9)*T(9))*IN(13)+(T(2)*T(3))*(T(6)*T(9)+T(9)*T(6))*IN(19)+(T(2)*T(3))*(T(9)*T(3)+T(3)*T(9))*IN(25)+(T(2)*T(3))*(T(3)*T(6)+T(6)*T(3))*IN(31)+(T(5)*T(6))*(T(3)*T(3))*IN(2)+(T(5)*T(6))*(T(6)*T(6))*IN(8)+(T(5)*T(6))*(T(9)*T(9))*IN(14)+(T(5)*T(6))*(T(6)*T(9)+T(9)*T(6))*IN(20)+(T(5)*T(6))*(T(9)*T(3)+T(3)*T(9))*IN(26)+(T(5)*T(6))*(T(3)*T(6)+T(6)*T(3))*IN(32)+(T(8)*T(9))*(T(3)*T(3))*IN(3)+(T(8)*T(9))*(T(6)*T(6))*IN(9)+(T(8)*T(9))*(T(9)*T(9))*IN(15)+(T(8)*T(9))*(T(6)*T(9)+T(9)*T(6))*IN(21)+(T(8)*T(9))*(T(9)*T(3)+T(3)*T(9))*IN(27)+(T(8)*T(9))*(T(3)*T(6)+T(6)*T(3))*IN(33)+(T(5)*T(9)+T(8)*T(6))*(T(3)*T(3))*IN(4)+(T(5)*T(9)+T(8)*T(6))*(T(6)*T(6))*IN(10)+(T(5)*T(9)+T(8)*T(6))*(T(9)*T(9))*IN(16)+(T(5)*T(9)+T(8)*T(6))*(T(6)*T(9)+T(9)*T(6))*IN(22)+(T(5)*T(9)+T(8)*T(6))*(T(9)*T(3)+T(3)*T(9))*IN(28)+(T(5)*T(9)+T(8)*T(6))*(T(3)*T(6)+T(6)*T(3))*IN(34)+(T(8)*T(3)+T(2)*T(9))*(T(3)*T(3))*IN(5)+(T(8)*T(3)+T(2)*T(9))*(T(6)*T(6))*IN(11)+(T(8)*T(3)+T(2)*T(9))*(T(9)*T(9))*IN(17)+(T(8)*T(3)+T(2)*T(9))*(T(6)*T(9)+T(9)*T(6))*IN(23)+(T(8)*T(3)+T(2)*T(9))*(T(9)*T(3)+T(3)*T(9))*IN(29)+(T(8)*T(3)+T(2)*T(9))*(T(3)*T(6)+T(6)*T(3))*IN(35)+(T(2)*T(6)+T(5)*T(3))*(T(3)*T(3))*IN(6)+(T(2)*T(6)+T(5)*T(3))*(T(6)*T(6))*IN(12)+(T(2)*T(6)+T(5)*T(3))*(T(9)*T(9))*IN(18)+(T(2)*T(6)+T(5)*T(3))*(T(6)*T(9)+T(9)*T(6))*IN(24)+(T(2)*T(6)+T(5)*T(3))*(T(9)*T(3)+T(3)*T(9))*IN(30)+(T(2)*T(6)+T(5)*T(3))*(T(3)*T(6)+T(6)*T(3))*IN(36));
S(22)=((T(2)*T(3))*(T(2)*T(3))*IN(1)+(T(2)*T(3))*(T(5)*T(6))*IN(7)+(T(2)*T(3))*(T(8)*T(9))*IN(13)+(T(2)*T(3))*(T(5)*T(9)+T(8)*T(6))*IN(19)+(T(2)*T(3))*(T(8)*T(3)+T(2)*T(9))*IN(25)+(T(2)*T(3))*(T(2)*T(6)+T(5)*T(3))*IN(31)+(T(5)*T(6))*(T(2)*T(3))*IN(2)+(T(5)*T(6))*(T(5)*T(6))*IN(8)+(T(5)*T(6))*(T(8)*T(9))*IN(14)+(T(5)*T(6))*(T(5)*T(9)+T(8)*T(6))*IN(20)+(T(5)*T(6))*(T(8)*T(3)+T(2)*T(9))*IN(26)+(T(5)*T(6))*(T(2)*T(6)+T(5)*T(3))*IN(32)+(T(8)*T(9))*(T(2)*T(3))*IN(3)+(T(8)*T(9))*(T(5)*T(6))*IN(9)+(T(8)*T(9))*(T(8)*T(9))*IN(15)+(T(8)*T(9))*(T(5)*T(9)+T(8)*T(6))*IN(21)+(T(8)*T(9))*(T(8)*T(3)+T(2)*T(9))*IN(27)+(T(8)*T(9))*(T(2)*T(6)+T(5)*T(3))*IN(33)+(T(5)*T(9)+T(8)*T(6))*(T(2)*T(3))*IN(4)+(T(5)*T(9)+T(8)*T(6))*(T(5)*T(6))*IN(10)+(T(5)*T(9)+T(8)*T(6))*(T(8)*T(9))*IN(16)+(T(5)*T(9)+T(8)*T(6))*(T(5)*T(9)+T(8)*T(6))*IN(22)+(T(5)*T(9)+T(8)*T(6))*(T(8)*T(3)+T(2)*T(9))*IN(28)+(T(5)*T(9)+T(8)*T(6))*(T(2)*T(6)+T(5)*T(3))*IN(34)+(T(8)*T(3)+T(2)*T(9))*(T(2)*T(3))*IN(5)+(T(8)*T(3)+T(2)*T(9))*(T(5)*T(6))*IN(11)+(T(8)*T(3)+T(2)*T(9))*(T(8)*T(9))*IN(17)+(T(8)*T(3)+T(2)*T(9))*(T(5)*T(9)+T(8)*T(6))*IN(23)+(T(8)*T(3)+T(2)*T(9))*(T(8)*T(3)+T(2)*T(9))*IN(29)+(T(8)*T(3)+T(2)*T(9))*(T(2)*T(6)+T(5)*T(3))*IN(35)+(T(2)*T(6)+T(5)*T(3))*(T(2)*T(3))*IN(6)+(T(2)*T(6)+T(5)*T(3))*(T(5)*T(6))*IN(12)+(T(2)*T(6)+T(5)*T(3))*(T(8)*T(9))*IN(18)+(T(2)*T(6)+T(5)*T(3))*(T(5)*T(9)+T(8)*T(6))*IN(24)+(T(2)*T(6)+T(5)*T(3))*(T(8)*T(3)+T(2)*T(9))*IN(30)+(T(2)*T(6)+T(5)*T(3))*(T(2)*T(6)+T(5)*T(3))*IN(36));
S(28)=((T(2)*T(3))*(T(3)*T(1))*IN(1)+(T(2)*T(3))*(T(6)*T(4))*IN(7)+(T(2)*T(3))*(T(9)*T(7))*IN(13)+(T(2)*T(3))*(T(6)*T(7)+T(9)*T(4))*IN(19)+(T(2)*T(3))*(T(9)*T(1)+T(3)*T(7))*IN(25)+(T(2)*T(3))*(T(3)*T(4)+T(6)*T(1))*IN(31)+(T(5)*T(6))*(T(3)*T(1))*IN(2)+(T(5)*T(6))*(T(6)*T(4))*IN(8)+(T(5)*T(6))*(T(9)*T(7))*IN(14)+(T(5)*T(6))*(T(6)*T(7)+T(9)*T(4))*IN(20)+(T(5)*T(6))*(T(9)*T(1)+T(3)*T(7))*IN(26)+(T(5)*T(6))*(T(3)*T(4)+T(6)*T(1))*IN(32)+(T(8)*T(9))*(T(3)*T(1))*IN(3)+(T(8)*T(9))*(T(6)*T(4))*IN(9)+(T(8)*T(9))*(T(9)*T(7))*IN(15)+(T(8)*T(9))*(T(6)*T(7)+T(9)*T(4))*IN(21)+(T(8)*T(9))*(T(9)*T(1)+T(3)*T(7))*IN(27)+(T(8)*T(9))*(T(3)*T(4)+T(6)*T(1))*IN(33)+(T(5)*T(9)+T(8)*T(6))*(T(3)*T(1))*IN(4)+(T(5)*T(9)+T(8)*T(6))*(T(6)*T(4))*IN(10)+(T(5)*T(9)+T(8)*T(6))*(T(9)*T(7))*IN(16)+(T(5)*T(9)+T(8)*T(6))*(T(6)*T(7)+T(9)*T(4))*IN(22)+(T(5)*T(9)+T(8)*T(6))*(T(9)*T(1)+T(3)*T(7))*IN(28)+(T(5)*T(9)+T(8)*T(6))*(T(3)*T(4)+T(6)*T(1))*IN(34)+(T(8)*T(3)+T(2)*T(9))*(T(3)*T(1))*IN(5)+(T(8)*T(3)+T(2)*T(9))*(T(6)*T(4))*IN(11)+(T(8)*T(3)+T(2)*T(9))*(T(9)*T(7))*IN(17)+(T(8)*T(3)+T(2)*T(9))*(T(6)*T(7)+T(9)*T(4))*IN(23)+(T(8)*T(3)+T(2)*T(9))*(T(9)*T(1)+T(3)*T(7))*IN(29)+(T(8)*T(3)+T(2)*T(9))*(T(3)*T(4)+T(6)*T(1))*IN(35)+(T(2)*T(6)+T(5)*T(3))*(T(3)*T(1))*IN(6)+(T(2)*T(6)+T(5)*T(3))*(T(6)*T(4))*IN(12)+(T(2)*T(6)+T(5)*T(3))*(T(9)*T(7))*IN(18)+(T(2)*T(6)+T(5)*T(3))*(T(6)*T(7)+T(9)*T(4))*IN(24)+(T(2)*T(6)+T(5)*T(3))*(T(9)*T(1)+T(3)*T(7))*IN(30)+(T(2)*T(6)+T(5)*T(3))*(T(3)*T(4)+T(6)*T(1))*IN(36));
S(34)=((T(2)*T(3))*(T(1)*T(2))*IN(1)+(T(2)*T(3))*(T(4)*T(5))*IN(7)+(T(2)*T(3))*(T(7)*T(8))*IN(13)+(T(2)*T(3))*(T(4)*T(8)+T(7)*T(5))*IN(19)+(T(2)*T(3))*(T(7)*T(2)+T(1)*T(8))*IN(25)+(T(2)*T(3))*(T(1)*T(5)+T(4)*T(2))*IN(31)+(T(5)*T(6))*(T(1)*T(2))*IN(2)+(T(5)*T(6))*(T(4)*T(5))*IN(8)+(T(5)*T(6))*(T(7)*T(8))*IN(14)+(T(5)*T(6))*(T(4)*T(8)+T(7)*T(5))*IN(20)+(T(5)*T(6))*(T(7)*T(2)+T(1)*T(8))*IN(26)+(T(5)*T(6))*(T(1)*T(5)+T(4)*T(2))*IN(32)+(T(8)*T(9))*(T(1)*T(2))*IN(3)+(T(8)*T(9))*(T(4)*T(5))*IN(9)+(T(8)*T(9))*(T(7)*T(8))*IN(15)+(T(8)*T(9))*(T(4)*T(8)+T(7)*T(5))*IN(21)+(T(8)*T(9))*(T(7)*T(2)+T(1)*T(8))*IN(27)+(T(8)*T(9))*(T(1)*T(5)+T(4)*T(2))*IN(33)+(T(5)*T(9)+T(8)*T(6))*(T(1)*T(2))*IN(4)+(T(5)*T(9)+T(8)*T(6))*(T(4)*T(5))*IN(10)+(T(5)*T(9)+T(8)*T(6))*(T(7)*T(8))*IN(16)+(T(5)*T(9)+T(8)*T(6))*(T(4)*T(8)+T(7)*T(5))*IN(22)+(T(5)*T(9)+T(8)*T(6))*(T(7)*T(2)+T(1)*T(8))*IN(28)+(T(5)*T(9)+T(8)*T(6))*(T(1)*T(5)+T(4)*T(2))*IN(34)+(T(8)*T(3)+T(2)*T(9))*(T(1)*T(2))*IN(5)+(T(8)*T(3)+T(2)*T(9))*(T(4)*T(5))*IN(11)+(T(8)*T(3)+T(2)*T(9))*(T(7)*T(8))*IN(17)+(T(8)*T(3)+T(2)*T(9))*(T(4)*T(8)+T(7)*T(5))*IN(23)+(T(8)*T(3)+T(2)*T(9))*(T(7)*T(2)+T(1)*T(8))*IN(29)+(T(8)*T(3)+T(2)*T(9))*(T(1)*T(5)+T(4)*T(2))*IN(35)+(T(2)*T(6)+T(5)*T(3))*(T(1)*T(2))*IN(6)+(T(2)*T(6)+T(5)*T(3))*(T(4)*T(5))*IN(12)+(T(2)*T(6)+T(5)*T(3))*(T(7)*T(8))*IN(18)+(T(2)*T(6)+T(5)*T(3))*(T(4)*T(8)+T(7)*T(5))*IN(24)+(T(2)*T(6)+T(5)*T(3))*(T(7)*T(2)+T(1)*T(8))*IN(30)+(T(2)*T(6)+T(5)*T(3))*(T(1)*T(5)+T(4)*T(2))*IN(36));
S(5)=((T(3)*T(1))*(T(1)*T(1))*IN(1)+(T(3)*T(1))*(T(4)*T(4))*IN(7)+(T(3)*T(1))*(T(7)*T(7))*IN(13)+(T(3)*T(1))*(T(4)*T(7)+T(7)*T(4))*IN(19)+(T(3)*T(1))*(T(7)*T(1)+T(1)*T(7))*IN(25)+(T(3)*T(1))*(T(1)*T(4)+T(4)*T(1))*IN(31)+(T(6)*T(4))*(T(1)*T(1))*IN(2)+(T(6)*T(4))*(T(4)*T(4))*IN(8)+(T(6)*T(4))*(T(7)*T(7))*IN(14)+(T(6)*T(4))*(T(4)*T(7)+T(7)*T(4))*IN(20)+(T(6)*T(4))*(T(7)*T(1)+T(1)*T(7))*IN(26)+(T(6)*T(4))*(T(1)*T(4)+T(4)*T(1))*IN(32)+(T(9)*T(7))*(T(1)*T(1))*IN(3)+(T(9)*T(7))*(T(4)*T(4))*IN(9)+(T(9)*T(7))*(T(7)*T(7))*IN(15)+(T(9)*T(7))*(T(4)*T(7)+T(7)*T(4))*IN(21)+(T(9)*T(7))*(T(7)*T(1)+T(1)*T(7))*IN(27)+(T(9)*T(7))*(T(1)*T(4)+T(4)*T(1))*IN(33)+(T(6)*T(7)+T(9)*T(4))*(T(1)*T(1))*IN(4)+(T(6)*T(7)+T(9)*T(4))*(T(4)*T(4))*IN(10)+(T(6)*T(7)+T(9)*T(4))*(T(7)*T(7))*IN(16)+(T(6)*T(7)+T(9)*T(4))*(T(4)*T(7)+T(7)*T(4))*IN(22)+(T(6)*T(7)+T(9)*T(4))*(T(7)*T(1)+T(1)*T(7))*IN(28)+(T(6)*T(7)+T(9)*T(4))*(T(1)*T(4)+T(4)*T(1))*IN(34)+(T(9)*T(1)+T(3)*T(7))*(T(1)*T(1))*IN(5)+(T(9)*T(1)+T(3)*T(7))*(T(4)*T(4))*IN(11)+(T(9)*T(1)+T(3)*T(7))*(T(7)*T(7))*IN(17)+(T(9)*T(1)+T(3)*T(7))*(T(4)*T(7)+T(7)*T(4))*IN(23)+(T(9)*T(1)+T(3)*T(7))*(T(7)*T(1)+T(1)*T(7))*IN(29)+(T(9)*T(1)+T(3)*T(7))*(T(1)*T(4)+T(4)*T(1))*IN(35)+(T(3)*T(4)+T(6)*T(1))*(T(1)*T(1))*IN(6)+(T(3)*T(4)+T(6)*T(1))*(T(4)*T(4))*IN(12)+(T(3)*T(4)+T(6)*T(1))*(T(7)*T(7))*IN(18)+(T(3)*T(4)+T(6)*T(1))*(T(4)*T(7)+T(7)*T(4))*IN(24)+(T(3)*T(4)+T(6)*T(1))*(T(7)*T(1)+T(1)*T(7))*IN(30)+(T(3)*T(4)+T(6)*T(1))*(T(1)*T(4)+T(4)*T(1))*IN(36));
S(11)=((T(3)*T(1))*(T(2)*T(2))*IN(1)+(T(3)*T(1))*(T(5)*T(5))*IN(7)+(T(3)*T(1))*(T(8)*T(8))*IN(13)+(T(3)*T(1))*(T(5)*T(8)+T(8)*T(5))*IN(19)+(T(3)*T(1))*(T(8)*T(2)+T(2)*T(8))*IN(25)+(T(3)*T(1))*(T(2)*T(5)+T(5)*T(2))*IN(31)+(T(6)*T(4))*(T(2)*T(2))*IN(2)+(T(6)*T(4))*(T(5)*T(5))*IN(8)+(T(6)*T(4))*(T(8)*T(8))*IN(14)+(T(6)*T(4))*(T(5)*T(8)+T(8)*T(5))*IN(20)+(T(6)*T(4))*(T(8)*T(2)+T(2)*T(8))*IN(26)+(T(6)*T(4))*(T(2)*T(5)+T(5)*T(2))*IN(32)+(T(9)*T(7))*(T(2)*T(2))*IN(3)+(T(9)*T(7))*(T(5)*T(5))*IN(9)+(T(9)*T(7))*(T(8)*T(8))*IN(15)+(T(9)*T(7))*(T(5)*T(8)+T(8)*T(5))*IN(21)+(T(9)*T(7))*(T(8)*T(2)+T(2)*T(8))*IN(27)+(T(9)*T(7))*(T(2)*T(5)+T(5)*T(2))*IN(33)+(T(6)*T(7)+T(9)*T(4))*(T(2)*T(2))*IN(4)+(T(6)*T(7)+T(9)*T(4))*(T(5)*T(5))*IN(10)+(T(6)*T(7)+T(9)*T(4))*(T(8)*T(8))*IN(16)+(T(6)*T(7)+T(9)*T(4))*(T(5)*T(8)+T(8)*T(5))*IN(22)+(T(6)*T(7)+T(9)*T(4))*(T(8)*T(2)+T(2)*T(8))*IN(28)+(T(6)*T(7)+T(9)*T(4))*(T(2)*T(5)+T(5)*T(2))*IN(34)+(T(9)*T(1)+T(3)*T(7))*(T(2)*T(2))*IN(5)+(T(9)*T(1)+T(3)*T(7))*(T(5)*T(5))*IN(11)+(T(9)*T(1)+T(3)*T(7))*(T(8)*T(8))*IN(17)+(T(9)*T(1)+T(3)*T(7))*(T(5)*T(8)+T(8)*T(5))*IN(23)+(T(9)*T(1)+T(3)*T(7))*(T(8)*T(2)+T(2)*T(8))*IN(29)+(T(9)*T(1)+T(3)*T(7))*(T(2)*T(5)+T(5)*T(2))*IN(35)+(T(3)*T(4)+T(6)*T(1))*(T(2)*T(2))*IN(6)+(T(3)*T(4)+T(6)*T(1))*(T(5)*T(5))*IN(12)+(T(3)*T(4)+T(6)*T(1))*(T(8)*T(8))*IN(18)+(T(3)*T(4)+T(6)*T(1))*(T(5)*T(8)+T(8)*T(5))*IN(24)+(T(3)*T(4)+T(6)*T(1))*(T(8)*T(2)+T(2)*T(8))*IN(30)+(T(3)*T(4)+T(6)*T(1))*(T(2)*T(5)+T(5)*T(2))*IN(36));
S(17)=((T(3)*T(1))*(T(3)*T(3))*IN(1)+(T(3)*T(1))*(T(6)*T(6))*IN(7)+(T(3)*T(1))*(T(9)*T(9))*IN(13)+(T(3)*T(1))*(T(6)*T(9)+T(9)*T(6))*IN(19)+(T(3)*T(1))*(T(9)*T(3)+T(3)*T(9))*IN(25)+(T(3)*T(1))*(T(3)*T(6)+T(6)*T(3))*IN(31)+(T(6)*T(4))*(T(3)*T(3))*IN(2)+(T(6)*T(4))*(T(6)*T(6))*IN(8)+(T(6)*T(4))*(T(9)*T(9))*IN(14)+(T(6)*T(4))*(T(6)*T(9)+T(9)*T(6))*IN(20)+(T(6)*T(4))*(T(9)*T(3)+T(3)*T(9))*IN(26)+(T(6)*T(4))*(T(3)*T(6)+T(6)*T(3))*IN(32)+(T(9)*T(7))*(T(3)*T(3))*IN(3)+(T(9)*T(7))*(T(6)*T(6))*IN(9)+(T(9)*T(7))*(T(9)*T(9))*IN(15)+(T(9)*T(7))*(T(6)*T(9)+T(9)*T(6))*IN(21)+(T(9)*T(7))*(T(9)*T(3)+T(3)*T(9))*IN(27)+(T(9)*T(7))*(T(3)*T(6)+T(6)*T(3))*IN(33)+(T(6)*T(7)+T(9)*T(4))*(T(3)*T(3))*IN(4)+(T(6)*T(7)+T(9)*T(4))*(T(6)*T(6))*IN(10)+(T(6)*T(7)+T(9)*T(4))*(T(9)*T(9))*IN(16)+(T(6)*T(7)+T(9)*T(4))*(T(6)*T(9)+T(9)*T(6))*IN(22)+(T(6)*T(7)+T(9)*T(4))*(T(9)*T(3)+T(3)*T(9))*IN(28)+(T(6)*T(7)+T(9)*T(4))*(T(3)*T(6)+T(6)*T(3))*IN(34)+(T(9)*T(1)+T(3)*T(7))*(T(3)*T(3))*IN(5)+(T(9)*T(1)+T(3)*T(7))*(T(6)*T(6))*IN(11)+(T(9)*T(1)+T(3)*T(7))*(T(9)*T(9))*IN(17)+(T(9)*T(1)+T(3)*T(7))*(T(6)*T(9)+T(9)*T(6))*IN(23)+(T(9)*T(1)+T(3)*T(7))*(T(9)*T(3)+T(3)*T(9))*IN(29)+(T(9)*T(1)+T(3)*T(7))*(T(3)*T(6)+T(6)*T(3))*IN(35)+(T(3)*T(4)+T(6)*T(1))*(T(3)*T(3))*IN(6)+(T(3)*T(4)+T(6)*T(1))*(T(6)*T(6))*IN(12)+(T(3)*T(4)+T(6)*T(1))*(T(9)*T(9))*IN(18)+(T(3)*T(4)+T(6)*T(1))*(T(6)*T(9)+T(9)*T(6))*IN(24)+(T(3)*T(4)+T(6)*T(1))*(T(9)*T(3)+T(3)*T(9))*IN(30)+(T(3)*T(4)+T(6)*T(1))*(T(3)*T(6)+T(6)*T(3))*IN(36));
S(23)=((T(3)*T(1))*(T(2)*T(3))*IN(1)+(T(3)*T(1))*(T(5)*T(6))*IN(7)+(T(3)*T(1))*(T(8)*T(9))*IN(13)+(T(3)*T(1))*(T(5)*T(9)+T(8)*T(6))*IN(19)+(T(3)*T(1))*(T(8)*T(3)+T(2)*T(9))*IN(25)+(T(3)*T(1))*(T(2)*T(6)+T(5)*T(3))*IN(31)+(T(6)*T(4))*(T(2)*T(3))*IN(2)+(T(6)*T(4))*(T(5)*T(6))*IN(8)+(T(6)*T(4))*(T(8)*T(9))*IN(14)+(T(6)*T(4))*(T(5)*T(9)+T(8)*T(6))*IN(20)+(T(6)*T(4))*(T(8)*T(3)+T(2)*T(9))*IN(26)+(T(6)*T(4))*(T(2)*T(6)+T(5)*T(3))*IN(32)+(T(9)*T(7))*(T(2)*T(3))*IN(3)+(T(9)*T(7))*(T(5)*T(6))*IN(9)+(T(9)*T(7))*(T(8)*T(9))*IN(15)+(T(9)*T(7))*(T(5)*T(9)+T(8)*T(6))*IN(21)+(T(9)*T(7))*(T(8)*T(3)+T(2)*T(9))*IN(27)+(T(9)*T(7))*(T(2)*T(6)+T(5)*T(3))*IN(33)+(T(6)*T(7)+T(9)*T(4))*(T(2)*T(3))*IN(4)+(T(6)*T(7)+T(9)*T(4))*(T(5)*T(6))*IN(10)+(T(6)*T(7)+T(9)*T(4))*(T(8)*T(9))*IN(16)+(T(6)*T(7)+T(9)*T(4))*(T(5)*T(9)+T(8)*T(6))*IN(22)+(T(6)*T(7)+T(9)*T(4))*(T(8)*T(3)+T(2)*T(9))*IN(28)+(T(6)*T(7)+T(9)*T(4))*(T(2)*T(6)+T(5)*T(3))*IN(34)+(T(9)*T(1)+T(3)*T(7))*(T(2)*T(3))*IN(5)+(T(9)*T(1)+T(3)*T(7))*(T(5)*T(6))*IN(11)+(T(9)*T(1)+T(3)*T(7))*(T(8)*T(9))*IN(17)+(T(9)*T(1)+T(3)*T(7))*(T(5)*T(9)+T(8)*T(6))*IN(23)+(T(9)*T(1)+T(3)*T(7))*(T(8)*T(3)+T(2)*T(9))*IN(29)+(T(9)*T(1)+T(3)*T(7))*(T(2)*T(6)+T(5)*T(3))*IN(35)+(T(3)*T(4)+T(6)*T(1))*(T(2)*T(3))*IN(6)+(T(3)*T(4)+T(6)*T(1))*(T(5)*T(6))*IN(12)+(T(3)*T(4)+T(6)*T(1))*(T(8)*T(9))*IN(18)+(T(3)*T(4)+T(6)*T(1))*(T(5)*T(9)+T(8)*T(6))*IN(24)+(T(3)*T(4)+T(6)*T(1))*(T(8)*T(3)+T(2)*T(9))*IN(30)+(T(3)*T(4)+T(6)*T(1))*(T(2)*T(6)+T(5)*T(3))*IN(36));
S(29)=((T(3)*T(1))*(T(3)*T(1))*IN(1)+(T(3)*T(1))*(T(6)*T(4))*IN(7)+(T(3)*T(1))*(T(9)*T(7))*IN(13)+(T(3)*T(1))*(T(6)*T(7)+T(9)*T(4))*IN(19)+(T(3)*T(1))*(T(9)*T(1)+T(3)*T(7))*IN(25)+(T(3)*T(1))*(T(3)*T(4)+T(6)*T(1))*IN(31)+(T(6)*T(4))*(T(3)*T(1))*IN(2)+(T(6)*T(4))*(T(6)*T(4))*IN(8)+(T(6)*T(4))*(T(9)*T(7))*IN(14)+(T(6)*T(4))*(T(6)*T(7)+T(9)*T(4))*IN(20)+(T(6)*T(4))*(T(9)*T(1)+T(3)*T(7))*IN(26)+(T(6)*T(4))*(T(3)*T(4)+T(6)*T(1))*IN(32)+(T(9)*T(7))*(T(3)*T(1))*IN(3)+(T(9)*T(7))*(T(6)*T(4))*IN(9)+(T(9)*T(7))*(T(9)*T(7))*IN(15)+(T(9)*T(7))*(T(6)*T(7)+T(9)*T(4))*IN(21)+(T(9)*T(7))*(T(9)*T(1)+T(3)*T(7))*IN(27)+(T(9)*T(7))*(T(3)*T(4)+T(6)*T(1))*IN(33)+(T(6)*T(7)+T(9)*T(4))*(T(3)*T(1))*IN(4)+(T(6)*T(7)+T(9)*T(4))*(T(6)*T(4))*IN(10)+(T(6)*T(7)+T(9)*T(4))*(T(9)*T(7))*IN(16)+(T(6)*T(7)+T(9)*T(4))*(T(6)*T(7)+T(9)*T(4))*IN(22)+(T(6)*T(7)+T(9)*T(4))*(T(9)*T(1)+T(3)*T(7))*IN(28)+(T(6)*T(7)+T(9)*T(4))*(T(3)*T(4)+T(6)*T(1))*IN(34)+(T(9)*T(1)+T(3)*T(7))*(T(3)*T(1))*IN(5)+(T(9)*T(1)+T(3)*T(7))*(T(6)*T(4))*IN(11)+(T(9)*T(1)+T(3)*T(7))*(T(9)*T(7))*IN(17)+(T(9)*T(1)+T(3)*T(7))*(T(6)*T(7)+T(9)*T(4))*IN(23)+(T(9)*T(1)+T(3)*T(7))*(T(9)*T(1)+T(3)*T(7))*IN(29)+(T(9)*T(1)+T(3)*T(7))*(T(3)*T(4)+T(6)*T(1))*IN(35)+(T(3)*T(4)+T(6)*T(1))*(T(3)*T(1))*IN(6)+(T(3)*T(4)+T(6)*T(1))*(T(6)*T(4))*IN(12)+(T(3)*T(4)+T(6)*T(1))*(T(9)*T(7))*IN(18)+(T(3)*T(4)+T(6)*T(1))*(T(6)*T(7)+T(9)*T(4))*IN(24)+(T(3)*T(4)+T(6)*T(1))*(T(9)*T(1)+T(3)*T(7))*IN(30)+(T(3)*T(4)+T(6)*T(1))*(T(3)*T(4)+T(6)*T(1))*IN(36));
S(35)=((T(3)*T(1))*(T(1)*T(2))*IN(1)+(T(3)*T(1))*(T(4)*T(5))*IN(7)+(T(3)*T(1))*(T(7)*T(8))*IN(13)+(T(3)*T(1))*(T(4)*T(8)+T(7)*T(5))*IN(19)+(T(3)*T(1))*(T(7)*T(2)+T(1)*T(8))*IN(25)+(T(3)*T(1))*(T(1)*T(5)+T(4)*T(2))*IN(31)+(T(6)*T(4))*(T(1)*T(2))*IN(2)+(T(6)*T(4))*(T(4)*T(5))*IN(8)+(T(6)*T(4))*(T(7)*T(8))*IN(14)+(T(6)*T(4))*(T(4)*T(8)+T(7)*T(5))*IN(20)+(T(6)*T(4))*(T(7)*T(2)+T(1)*T(8))*IN(26)+(T(6)*T(4))*(T(1)*T(5)+T(4)*T(2))*IN(32)+(T(9)*T(7))*(T(1)*T(2))*IN(3)+(T(9)*T(7))*(T(4)*T(5))*IN(9)+(T(9)*T(7))*(T(7)*T(8))*IN(15)+(T(9)*T(7))*(T(4)*T(8)+T(7)*T(5))*IN(21)+(T(9)*T(7))*(T(7)*T(2)+T(1)*T(8))*IN(27)+(T(9)*T(7))*(T(1)*T(5)+T(4)*T(2))*IN(33)+(T(6)*T(7)+T(9)*T(4))*(T(1)*T(2))*IN(4)+(T(6)*T(7)+T(9)*T(4))*(T(4)*T(5))*IN(10)+(T(6)*T(7)+T(9)*T(4))*(T(7)*T(8))*IN(16)+(T(6)*T(7)+T(9)*T(4))*(T(4)*T(8)+T(7)*T(5))*IN(22)+(T(6)*T(7)+T(9)*T(4))*(T(7)*T(2)+T(1)*T(8))*IN(28)+(T(6)*T(7)+T(9)*T(4))*(T(1)*T(5)+T(4)*T(2))*IN(34)+(T(9)*T(1)+T(3)*T(7))*(T(1)*T(2))*IN(5)+(T(9)*T(1)+T(3)*T(7))*(T(4)*T(5))*IN(11)+(T(9)*T(1)+T(3)*T(7))*(T(7)*T(8))*IN(17)+(T(9)*T(1)+T(3)*T(7))*(T(4)*T(8)+T(7)*T(5))*IN(23)+(T(9)*T(1)+T(3)*T(7))*(T(7)*T(2)+T(1)*T(8))*IN(29)+(T(9)*T(1)+T(3)*T(7))*(T(1)*T(5)+T(4)*T(2))*IN(35)+(T(3)*T(4)+T(6)*T(1))*(T(1)*T(2))*IN(6)+(T(3)*T(4)+T(6)*T(1))*(T(4)*T(5))*IN(12)+(T(3)*T(4)+T(6)*T(1))*(T(7)*T(8))*IN(18)+(T(3)*T(4)+T(6)*T(1))*(T(4)*T(8)+T(7)*T(5))*IN(24)+(T(3)*T(4)+T(6)*T(1))*(T(7)*T(2)+T(1)*T(8))*IN(30)+(T(3)*T(4)+T(6)*T(1))*(T(1)*T(5)+T(4)*T(2))*IN(36));
S(6)=((T(1)*T(2))*(T(1)*T(1))*IN(1)+(T(1)*T(2))*(T(4)*T(4))*IN(7)+(T(1)*T(2))*(T(7)*T(7))*IN(13)+(T(1)*T(2))*(T(4)*T(7)+T(7)*T(4))*IN(19)+(T(1)*T(2))*(T(7)*T(1)+T(1)*T(7))*IN(25)+(T(1)*T(2))*(T(1)*T(4)+T(4)*T(1))*IN(31)+(T(4)*T(5))*(T(1)*T(1))*IN(2)+(T(4)*T(5))*(T(4)*T(4))*IN(8)+(T(4)*T(5))*(T(7)*T(7))*IN(14)+(T(4)*T(5))*(T(4)*T(7)+T(7)*T(4))*IN(20)+(T(4)*T(5))*(T(7)*T(1)+T(1)*T(7))*IN(26)+(T(4)*T(5))*(T(1)*T(4)+T(4)*T(1))*IN(32)+(T(7)*T(8))*(T(1)*T(1))*IN(3)+(T(7)*T(8))*(T(4)*T(4))*IN(9)+(T(7)*T(8))*(T(7)*T(7))*IN(15)+(T(7)*T(8))*(T(4)*T(7)+T(7)*T(4))*IN(21)+(T(7)*T(8))*(T(7)*T(1)+T(1)*T(7))*IN(27)+(T(7)*T(8))*(T(1)*T(4)+T(4)*T(1))*IN(33)+(T(4)*T(8)+T(7)*T(5))*(T(1)*T(1))*IN(4)+(T(4)*T(8)+T(7)*T(5))*(T(4)*T(4))*IN(10)+(T(4)*T(8)+T(7)*T(5))*(T(7)*T(7))*IN(16)+(T(4)*T(8)+T(7)*T(5))*(T(4)*T(7)+T(7)*T(4))*IN(22)+(T(4)*T(8)+T(7)*T(5))*(T(7)*T(1)+T(1)*T(7))*IN(28)+(T(4)*T(8)+T(7)*T(5))*(T(1)*T(4)+T(4)*T(1))*IN(34)+(T(7)*T(2)+T(1)*T(8))*(T(1)*T(1))*IN(5)+(T(7)*T(2)+T(1)*T(8))*(T(4)*T(4))*IN(11)+(T(7)*T(2)+T(1)*T(8))*(T(7)*T(7))*IN(17)+(T(7)*T(2)+T(1)*T(8))*(T(4)*T(7)+T(7)*T(4))*IN(23)+(T(7)*T(2)+T(1)*T(8))*(T(7)*T(1)+T(1)*T(7))*IN(29)+(T(7)*T(2)+T(1)*T(8))*(T(1)*T(4)+T(4)*T(1))*IN(35)+(T(1)*T(5)+T(4)*T(2))*(T(1)*T(1))*IN(6)+(T(1)*T(5)+T(4)*T(2))*(T(4)*T(4))*IN(12)+(T(1)*T(5)+T(4)*T(2))*(T(7)*T(7))*IN(18)+(T(1)*T(5)+T(4)*T(2))*(T(4)*T(7)+T(7)*T(4))*IN(24)+(T(1)*T(5)+T(4)*T(2))*(T(7)*T(1)+T(1)*T(7))*IN(30)+(T(1)*T(5)+T(4)*T(2))*(T(1)*T(4)+T(4)*T(1))*IN(36));
S(12)=((T(1)*T(2))*(T(2)*T(2))*IN(1)+(T(1)*T(2))*(T(5)*T(5))*IN(7)+(T(1)*T(2))*(T(8)*T(8))*IN(13)+(T(1)*T(2))*(T(5)*T(8)+T(8)*T(5))*IN(19)+(T(1)*T(2))*(T(8)*T(2)+T(2)*T(8))*IN(25)+(T(1)*T(2))*(T(2)*T(5)+T(5)*T(2))*IN(31)+(T(4)*T(5))*(T(2)*T(2))*IN(2)+(T(4)*T(5))*(T(5)*T(5))*IN(8)+(T(4)*T(5))*(T(8)*T(8))*IN(14)+(T(4)*T(5))*(T(5)*T(8)+T(8)*T(5))*IN(20)+(T(4)*T(5))*(T(8)*T(2)+T(2)*T(8))*IN(26)+(T(4)*T(5))*(T(2)*T(5)+T(5)*T(2))*IN(32)+(T(7)*T(8))*(T(2)*T(2))*IN(3)+(T(7)*T(8))*(T(5)*T(5))*IN(9)+(T(7)*T(8))*(T(8)*T(8))*IN(15)+(T(7)*T(8))*(T(5)*T(8)+T(8)*T(5))*IN(21)+(T(7)*T(8))*(T(8)*T(2)+T(2)*T(8))*IN(27)+(T(7)*T(8))*(T(2)*T(5)+T(5)*T(2))*IN(33)+(T(4)*T(8)+T(7)*T(5))*(T(2)*T(2))*IN(4)+(T(4)*T(8)+T(7)*T(5))*(T(5)*T(5))*IN(10)+(T(4)*T(8)+T(7)*T(5))*(T(8)*T(8))*IN(16)+(T(4)*T(8)+T(7)*T(5))*(T(5)*T(8)+T(8)*T(5))*IN(22)+(T(4)*T(8)+T(7)*T(5))*(T(8)*T(2)+T(2)*T(8))*IN(28)+(T(4)*T(8)+T(7)*T(5))*(T(2)*T(5)+T(5)*T(2))*IN(34)+(T(7)*T(2)+T(1)*T(8))*(T(2)*T(2))*IN(5)+(T(7)*T(2)+T(1)*T(8))*(T(5)*T(5))*IN(11)+(T(7)*T(2)+T(1)*T(8))*(T(8)*T(8))*IN(17)+(T(7)*T(2)+T(1)*T(8))*(T(5)*T(8)+T(8)*T(5))*IN(23)+(T(7)*T(2)+T(1)*T(8))*(T(8)*T(2)+T(2)*T(8))*IN(29)+(T(7)*T(2)+T(1)*T(8))*(T(2)*T(5)+T(5)*T(2))*IN(35)+(T(1)*T(5)+T(4)*T(2))*(T(2)*T(2))*IN(6)+(T(1)*T(5)+T(4)*T(2))*(T(5)*T(5))*IN(12)+(T(1)*T(5)+T(4)*T(2))*(T(8)*T(8))*IN(18)+(T(1)*T(5)+T(4)*T(2))*(T(5)*T(8)+T(8)*T(5))*IN(24)+(T(1)*T(5)+T(4)*T(2))*(T(8)*T(2)+T(2)*T(8))*IN(30)+(T(1)*T(5)+T(4)*T(2))*(T(2)*T(5)+T(5)*T(2))*IN(36));
S(18)=((T(1)*T(2))*(T(3)*T(3))*IN(1)+(T(1)*T(2))*(T(6)*T(6))*IN(7)+(T(1)*T(2))*(T(9)*T(9))*IN(13)+(T(1)*T(2))*(T(6)*T(9)+T(9)*T(6))*IN(19)+(T(1)*T(2))*(T(9)*T(3)+T(3)*T(9))*IN(25)+(T(1)*T(2))*(T(3)*T(6)+T(6)*T(3))*IN(31)+(T(4)*T(5))*(T(3)*T(3))*IN(2)+(T(4)*T(5))*(T(6)*T(6))*IN(8)+(T(4)*T(5))*(T(9)*T(9))*IN(14)+(T(4)*T(5))*(T(6)*T(9)+T(9)*T(6))*IN(20)+(T(4)*T(5))*(T(9)*T(3)+T(3)*T(9))*IN(26)+(T(4)*T(5))*(T(3)*T(6)+T(6)*T(3))*IN(32)+(T(7)*T(8))*(T(3)*T(3))*IN(3)+(T(7)*T(8))*(T(6)*T(6))*IN(9)+(T(7)*T(8))*(T(9)*T(9))*IN(15)+(T(7)*T(8))*(T(6)*T(9)+T(9)*T(6))*IN(21)+(T(7)*T(8))*(T(9)*T(3)+T(3)*T(9))*IN(27)+(T(7)*T(8))*(T(3)*T(6)+T(6)*T(3))*IN(33)+(T(4)*T(8)+T(7)*T(5))*(T(3)*T(3))*IN(4)+(T(4)*T(8)+T(7)*T(5))*(T(6)*T(6))*IN(10)+(T(4)*T(8)+T(7)*T(5))*(T(9)*T(9))*IN(16)+(T(4)*T(8)+T(7)*T(5))*(T(6)*T(9)+T(9)*T(6))*IN(22)+(T(4)*T(8)+T(7)*T(5))*(T(9)*T(3)+T(3)*T(9))*IN(28)+(T(4)*T(8)+T(7)*T(5))*(T(3)*T(6)+T(6)*T(3))*IN(34)+(T(7)*T(2)+T(1)*T(8))*(T(3)*T(3))*IN(5)+(T(7)*T(2)+T(1)*T(8))*(T(6)*T(6))*IN(11)+(T(7)*T(2)+T(1)*T(8))*(T(9)*T(9))*IN(17)+(T(7)*T(2)+T(1)*T(8))*(T(6)*T(9)+T(9)*T(6))*IN(23)+(T(7)*T(2)+T(1)*T(8))*(T(9)*T(3)+T(3)*T(9))*IN(29)+(T(7)*T(2)+T(1)*T(8))*(T(3)*T(6)+T(6)*T(3))*IN(35)+(T(1)*T(5)+T(4)*T(2))*(T(3)*T(3))*IN(6)+(T(1)*T(5)+T(4)*T(2))*(T(6)*T(6))*IN(12)+(T(1)*T(5)+T(4)*T(2))*(T(9)*T(9))*IN(18)+(T(1)*T(5)+T(4)*T(2))*(T(6)*T(9)+T(9)*T(6))*IN(24)+(T(1)*T(5)+T(4)*T(2))*(T(9)*T(3)+T(3)*T(9))*IN(30)+(T(1)*T(5)+T(4)*T(2))*(T(3)*T(6)+T(6)*T(3))*IN(36));
S(24)=((T(1)*T(2))*(T(2)*T(3))*IN(1)+(T(1)*T(2))*(T(5)*T(6))*IN(7)+(T(1)*T(2))*(T(8)*T(9))*IN(13)+(T(1)*T(2))*(T(5)*T(9)+T(8)*T(6))*IN(19)+(T(1)*T(2))*(T(8)*T(3)+T(2)*T(9))*IN(25)+(T(1)*T(2))*(T(2)*T(6)+T(5)*T(3))*IN(31)+(T(4)*T(5))*(T(2)*T(3))*IN(2)+(T(4)*T(5))*(T(5)*T(6))*IN(8)+(T(4)*T(5))*(T(8)*T(9))*IN(14)+(T(4)*T(5))*(T(5)*T(9)+T(8)*T(6))*IN(20)+(T(4)*T(5))*(T(8)*T(3)+T(2)*T(9))*IN(26)+(T(4)*T(5))*(T(2)*T(6)+T(5)*T(3))*IN(32)+(T(7)*T(8))*(T(2)*T(3))*IN(3)+(T(7)*T(8))*(T(5)*T(6))*IN(9)+(T(7)*T(8))*(T(8)*T(9))*IN(15)+(T(7)*T(8))*(T(5)*T(9)+T(8)*T(6))*IN(21)+(T(7)*T(8))*(T(8)*T(3)+T(2)*T(9))*IN(27)+(T(7)*T(8))*(T(2)*T(6)+T(5)*T(3))*IN(33)+(T(4)*T(8)+T(7)*T(5))*(T(2)*T(3))*IN(4)+(T(4)*T(8)+T(7)*T(5))*(T(5)*T(6))*IN(10)+(T(4)*T(8)+T(7)*T(5))*(T(8)*T(9))*IN(16)+(T(4)*T(8)+T(7)*T(5))*(T(5)*T(9)+T(8)*T(6))*IN(22)+(T(4)*T(8)+T(7)*T(5))*(T(8)*T(3)+T(2)*T(9))*IN(28)+(T(4)*T(8)+T(7)*T(5))*(T(2)*T(6)+T(5)*T(3))*IN(34)+(T(7)*T(2)+T(1)*T(8))*(T(2)*T(3))*IN(5)+(T(7)*T(2)+T(1)*T(8))*(T(5)*T(6))*IN(11)+(T(7)*T(2)+T(1)*T(8))*(T(8)*T(9))*IN(17)+(T(7)*T(2)+T(1)*T(8))*(T(5)*T(9)+T(8)*T(6))*IN(23)+(T(7)*T(2)+T(1)*T(8))*(T(8)*T(3)+T(2)*T(9))*IN(29)+(T(7)*T(2)+T(1)*T(8))*(T(2)*T(6)+T(5)*T(3))*IN(35)+(T(1)*T(5)+T(4)*T(2))*(T(2)*T(3))*IN(6)+(T(1)*T(5)+T(4)*T(2))*(T(5)*T(6))*IN(12)+(T(1)*T(5)+T(4)*T(2))*(T(8)*T(9))*IN(18)+(T(1)*T(5)+T(4)*T(2))*(T(5)*T(9)+T(8)*T(6))*IN(24)+(T(1)*T(5)+T(4)*T(2))*(T(8)*T(3)+T(2)*T(9))*IN(30)+(T(1)*T(5)+T(4)*T(2))*(T(2)*T(6)+T(5)*T(3))*IN(36));
S(30)=((T(1)*T(2))*(T(3)*T(1))*IN(1)+(T(1)*T(2))*(T(6)*T(4))*IN(7)+(T(1)*T(2))*(T(9)*T(7))*IN(13)+(T(1)*T(2))*(T(6)*T(7)+T(9)*T(4))*IN(19)+(T(1)*T(2))*(T(9)*T(1)+T(3)*T(7))*IN(25)+(T(1)*T(2))*(T(3)*T(4)+T(6)*T(1))*IN(31)+(T(4)*T(5))*(T(3)*T(1))*IN(2)+(T(4)*T(5))*(T(6)*T(4))*IN(8)+(T(4)*T(5))*(T(9)*T(7))*IN(14)+(T(4)*T(5))*(T(6)*T(7)+T(9)*T(4))*IN(20)+(T(4)*T(5))*(T(9)*T(1)+T(3)*T(7))*IN(26)+(T(4)*T(5))*(T(3)*T(4)+T(6)*T(1))*IN(32)+(T(7)*T(8))*(T(3)*T(1))*IN(3)+(T(7)*T(8))*(T(6)*T(4))*IN(9)+(T(7)*T(8))*(T(9)*T(7))*IN(15)+(T(7)*T(8))*(T(6)*T(7)+T(9)*T(4))*IN(21)+(T(7)*T(8))*(T(9)*T(1)+T(3)*T(7))*IN(27)+(T(7)*T(8))*(T(3)*T(4)+T(6)*T(1))*IN(33)+(T(4)*T(8)+T(7)*T(5))*(T(3)*T(1))*IN(4)+(T(4)*T(8)+T(7)*T(5))*(T(6)*T(4))*IN(10)+(T(4)*T(8)+T(7)*T(5))*(T(9)*T(7))*IN(16)+(T(4)*T(8)+T(7)*T(5))*(T(6)*T(7)+T(9)*T(4))*IN(22)+(T(4)*T(8)+T(7)*T(5))*(T(9)*T(1)+T(3)*T(7))*IN(28)+(T(4)*T(8)+T(7)*T(5))*(T(3)*T(4)+T(6)*T(1))*IN(34)+(T(7)*T(2)+T(1)*T(8))*(T(3)*T(1))*IN(5)+(T(7)*T(2)+T(1)*T(8))*(T(6)*T(4))*IN(11)+(T(7)*T(2)+T(1)*T(8))*(T(9)*T(7))*IN(17)+(T(7)*T(2)+T(1)*T(8))*(T(6)*T(7)+T(9)*T(4))*IN(23)+(T(7)*T(2)+T(1)*T(8))*(T(9)*T(1)+T(3)*T(7))*IN(29)+(T(7)*T(2)+T(1)*T(8))*(T(3)*T(4)+T(6)*T(1))*IN(35)+(T(1)*T(5)+T(4)*T(2))*(T(3)*T(1))*IN(6)+(T(1)*T(5)+T(4)*T(2))*(T(6)*T(4))*IN(12)+(T(1)*T(5)+T(4)*T(2))*(T(9)*T(7))*IN(18)+(T(1)*T(5)+T(4)*T(2))*(T(6)*T(7)+T(9)*T(4))*IN(24)+(T(1)*T(5)+T(4)*T(2))*(T(9)*T(1)+T(3)*T(7))*IN(30)+(T(1)*T(5)+T(4)*T(2))*(T(3)*T(4)+T(6)*T(1))*IN(36));
S(36)=((T(1)*T(2))*(T(1)*T(2))*IN(1)+(T(1)*T(2))*(T(4)*T(5))*IN(7)+(T(1)*T(2))*(T(7)*T(8))*IN(13)+(T(1)*T(2))*(T(4)*T(8)+T(7)*T(5))*IN(19)+(T(1)*T(2))*(T(7)*T(2)+T(1)*T(8))*IN(25)+(T(1)*T(2))*(T(1)*T(5)+T(4)*T(2))*IN(31)+(T(4)*T(5))*(T(1)*T(2))*IN(2)+(T(4)*T(5))*(T(4)*T(5))*IN(8)+(T(4)*T(5))*(T(7)*T(8))*IN(14)+(T(4)*T(5))*(T(4)*T(8)+T(7)*T(5))*IN(20)+(T(4)*T(5))*(T(7)*T(2)+T(1)*T(8))*IN(26)+(T(4)*T(5))*(T(1)*T(5)+T(4)*T(2))*IN(32)+(T(7)*T(8))*(T(1)*T(2))*IN(3)+(T(7)*T(8))*(T(4)*T(5))*IN(9)+(T(7)*T(8))*(T(7)*T(8))*IN(15)+(T(7)*T(8))*(T(4)*T(8)+T(7)*T(5))*IN(21)+(T(7)*T(8))*(T(7)*T(2)+T(1)*T(8))*IN(27)+(T(7)*T(8))*(T(1)*T(5)+T(4)*T(2))*IN(33)+(T(4)*T(8)+T(7)*T(5))*(T(1)*T(2))*IN(4)+(T(4)*T(8)+T(7)*T(5))*(T(4)*T(5))*IN(10)+(T(4)*T(8)+T(7)*T(5))*(T(7)*T(8))*IN(16)+(T(4)*T(8)+T(7)*T(5))*(T(4)*T(8)+T(7)*T(5))*IN(22)+(T(4)*T(8)+T(7)*T(5))*(T(7)*T(2)+T(1)*T(8))*IN(28)+(T(4)*T(8)+T(7)*T(5))*(T(1)*T(5)+T(4)*T(2))*IN(34)+(T(7)*T(2)+T(1)*T(8))*(T(1)*T(2))*IN(5)+(T(7)*T(2)+T(1)*T(8))*(T(4)*T(5))*IN(11)+(T(7)*T(2)+T(1)*T(8))*(T(7)*T(8))*IN(17)+(T(7)*T(2)+T(1)*T(8))*(T(4)*T(8)+T(7)*T(5))*IN(23)+(T(7)*T(2)+T(1)*T(8))*(T(7)*T(2)+T(1)*T(8))*IN(29)+(T(7)*T(2)+T(1)*T(8))*(T(1)*T(5)+T(4)*T(2))*IN(35)+(T(1)*T(5)+T(4)*T(2))*(T(1)*T(2))*IN(6)+(T(1)*T(5)+T(4)*T(2))*(T(4)*T(5))*IN(12)+(T(1)*T(5)+T(4)*T(2))*(T(7)*T(8))*IN(18)+(T(1)*T(5)+T(4)*T(2))*(T(4)*T(8)+T(7)*T(5))*IN(24)+(T(1)*T(5)+T(4)*T(2))*(T(7)*T(2)+T(1)*T(8))*IN(30)+(T(1)*T(5)+T(4)*T(2))*(T(1)*T(5)+T(4)*T(2))*IN(36));
end
 
