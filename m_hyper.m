function [out,out1,out2]=m_hyper(varargin)

%M_HYPER hyperelastic material function
%
%       Syntax : mat= m_hyper('default') 
%                mat= m_hyper('database name') 
%                pl = m_hyper('dbval MatId name');
%                pl = m_hyper('dbval -unit MM MatId name');
%                pl = m_hyper('dbval -punit MM MatId name');
%
%       Material subtypes supported by m_hyper are
%
%	Subtype 1 : [MatId type rho wtype C_1 C_2 K Kv]
%     with
%      type  fe_mat('m_hyper','SI',1)
%      rho(density)
%      wtype(choice energy)
%      0 : Money rivlin + linear compression [MatId type rho wtype C_1 C_2 K Kv]
%        W = C_1(J_1-3) + C_2(J_2-3) + K(J_3-1)^2
%      1 : Money Rivlin + Ciarlet Gemonat  [MatId type rho wtype C_1 C_2 K Kv]
%        W = C_1(J_1-3) + C_2(J_2-3) + K(J_3-1) - (C_1 + 2C_2 + K)\ln(J_3) 
%           C_1  C_2  K (energy coefficients)
%      2 : Yeoh + Ciarlet Gemonat [MatId type rho wtype C10 C20 K Kv C_i0 -1]
%
%       See sdtweb     m_elastic, pl, fem
%       See also help  fe_mat, p_shell, p_beam

%	Subtype 2 : Nominal Hyperelastic law with viscosity for dynamic case
%	    [MatId type  rho wtype C_1 C_2 K eta dt]
%            eta = viscous coefficient
%            dt = time step

%	Etienne Balmes, Jean-Michel Leclere, Corine Florens
%       Copyright (c) 2001-2022 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*STREMP>

if nargin<1; help m_hyper;return; end
if ischar(varargin{1}); [CAM,Cam]=comstr(varargin{1},1); pl=[]; carg=2;
else; pl=varargin{1};[CAM,Cam]=comstr(varargin{2},1); carg=3;
end


%% #PropertyUnitType -- m_hyper('propertyunittype cell',1) ---------
if comstr(Cam,'propertyunittype')

 if nargin==1; out=1;return;end
 i1=varargin{carg};
 out1={};
 switch i1
 case 1
 st=...
 {'MatId'  0  'sdtweb(''m_hyper'')';
  'Type'   0  '';
  'Rho'    3  'Density';
  'Wtype'  0  'Energy type';
  'C1'     1  'Energy coefficient';
  'C2'     1  'Energy coefficient';
  'K'      1  'Bulk stiffness';
  'Kv'     1  'Bulk viscosity';
  'C3'     1  'Energy coefficient';
 };
 out1={         'Ca%i'            1 'Energy coefficient';
 };
 case 2

 otherwise; st={'MatId' 0; 'Type', 0};
 end
 if ~isempty(strfind(Cam,'cell')); out=st; else; out=[st{:,2}]; end 

% -------------------------------------------------------------------------

%% #Info ---------------------------------------------------------------------
elseif comstr(Cam,'info')

 r1=m_hyper('database');
 disp(char({r1.name}))

%% #Default ------------------------------------------------------------------
elseif comstr(Cam,'default')

  out=m_hyper('database'); out=out(1);

%% #DBval --------------------------------------------------------------------
elseif comstr(Cam,'dbval') 

 while 1==1
  i1=strfind(comstr(Cam,-27),'-unit'); out1={}; % See if unit specified
  if ~isempty(i1)
   [Unit,i2,i3,i4]=sscanf(CAM(i1+5:end),'%s',1);
   i4=i1+[0:4+i4];CAM(i4)=''; [CAM,Cam]=comstr(CAM,1);
  else; Unit='';
  end
  i2=strfind(comstr(Cam,-27),'-punit');
  if ~isempty(i2)
   [PUnit,i3,i4,i5]=sscanf(CAM(i2+6:end),'%s',1);
   i5=i2+[0:5+i5];CAM(i5)=''; [CAM,Cam]=comstr(CAM,1);
  else; PUnit='';
  end
  if ischar(CAM); [i1,CAM,Cam]=comstr(CAM,'dbval','%i');else; i1=[];end
  if isempty(CAM)&&carg<=nargin; st=varargin{carg};carg=carg+1;
  else st=CAM;end
  if isempty(st)
  elseif ischar(st); mat=m_hyper('database',st);
  elseif isnumeric(st)
   [typ,st1,i4]=fe_mat('typem',st(2));
   mat=struct('pl',st,'name',sprintf('%i',st(1)),'type',typ,'unit',st1);
  end
  if ~isempty(PUnit)
   r1=fe_mat(sprintf('convert %s %s',mat.unit,PUnit),mat.pl(1:2));
   mat.pl(2)=r1(2); mat.unit=PUnit;
  end
  if ~isempty(Unit)
   mat.pl=fe_mat(sprintf('convert %s %s',mat.unit,Unit),mat.pl);mat.unit=Unit;
  end
  if length(i1)==1; mat.pl(1)=i1;end
  r1=mat.pl;
  if ~isempty(pl); i2=find(pl(:,1)==r1(1)); else i2=[];end
  if isempty(i2); i2=size(pl,1)+1;end
  pl(i2,1:length(r1))=r1;  %#ok<AGROW>
  out1(end+1,1:3)={'mat',sprintf('%i_%s',mat.pl(1),mat.name),mat};%#ok<AGROW>
  if carg>nargin; break;end
  CAM=varargin{carg};carg=carg+1;if ischar(CAM);[CAM,Cam]=comstr(CAM,1);end
 end
 out=pl;

elseif comstr(Cam,'urn')
%% #URN based building of constants --------------------------------
% r1=m_hyper('urn','SimoA{1,1,0,30,3,f5 20,g .33 .33,rho2.33n}');
% r1=m_hyper('urn','PadA{2.5264,-0.9177,0.4711,1200,3,f .35,g .5688,rho1n,tyYeoh,unTM}')
% Youssera El Archi (these LMA 2022) : confined pressure test k 1280 MPa 
% Zhuravlev table 2.3 page 48 : {Yeoh,C10 2.5264, C20 -0.9177, C30 0.4711, g1 0.5688, tau 2.635}

RO=varargin{carg};carg=carg+1; 
if ischar(RO);RO=struct('urn',RO);end
[CAM,r1]=sdtm.urnPar(RO.urn,'{c1%ug,c2%ug,c3%ug,kappa%ug,kappav%ug}{f%ug,g%ug,rho%ug,ty%s,un%s}');

if ~isfield(r1,'g'); r1.g=[];r1.f=[];end % Possible relaxation cells
if ~isfield(r1,'un');r1.un='US';end

NL=struct('type','nl_inout','opt',[0 0 0],'MexCb',{{m_hyper('@hyper_g'),struct}}, ...
              'adofi',zeros(9*(length(r1.g)+1)+2,1)-.99);
NL=sdth.sfield('addselected',NL,r1,setdiff(fieldnames(r1), ...
    {'c1','c2','c3','kappa','kappav','rho','ty','un'}));
if isfield(r1,'ty')
 switch lower(r1.ty)
 case 'yeoh'
  % G=2*C10 Zhuravlev table 2.2
  NL.mtype=r1.ty;
  out=struct('pl',[],'NLdata',NL,'name',CAM,'unit',r1.un);
  out.pl=[1 fe_mat('m_hyper',r1.un,1) r1.rho 2 r1.c1 r1.c2 r1.kappa r1.kappav r1.c3 -1];         
 case 'moon'
  NL.mtype='Moon';
  out=struct('pl',[],'NLdata',NL,'name',CAM,'unit',r1.un);
  out.pl=[1 fe_mat('m_hyper',r1.un,1) r1.rho 1 r1.c1 r1.c2 r1.kappa r1.kappav r1.c3 -1];         
  
 otherwise; error('Not yet implemented')
 end
else % Mooney Rivlin
 RO.G=2*(r1.c1+r1.c2);  % G=mu=2/(1+nu)
 RO.E=9*r1.kappa*RO.G/(3*r1.kappa+RO.G); %3G, if kappa>>G
 RO.nu=3*r1.kappa-RO.E/6/r1.kappa;
 out=struct('pl',[],'NLdata',NL,'name',CAM,'unit',r1.un);
 out.pl=[1 fe_mat('m_hyper',r1.un,1) r1.rho 1 r1.c1 r1.c2 r1.kappa r1.kappav r1.c3];
end
if nargout==0; % Verifications
 out.check2=1; checkNL('ViewA',out);clear out; 
end

elseif comstr(Cam,'database'); [st,Cam]=comstr(CAM,9);
%% #Database -----------------------------------------------------------------

  if isempty(st)&&carg<=nargin; st=varargin{carg}; carg=carg+1;end
  [MatId,i2,i3,i4]=sscanf(st,'%i',1); if i2~=1; MatId=1;end
  st=comstr(st,i4);

  out.pl=[MatId fe_mat('type','m_hyper','SI',1) 1e-06 0 .3 .2 .3 0.1 1e-2];
  out.name='Ref';
  out.type='m_hyper';
  out.unit='SI';

  out(2).pl=[MatId fe_mat('type','m_hyper','SI',1) 2e-2 1 0.3 0.2 0.3 0.1 1e-2];
  out(2).name='Rivlin';
  out(2).type='m_hyper';
  out(2).unit='SI';
  
  out(3).pl=[MatId fe_mat('type','m_hyper','SI',1) 1e3 1 4.2e3 3.3e1 6.3e5 6.8e2 1e-3];
  out(3).name='phi';
  out(3).type='m_hyper';
  out(3).unit='SI';


  out(4).pl=[MatId fe_mat('type','m_hyper','SI',1) 1e-6 1 25. 0.25 1274. 1.2e1 1e-3];
  out(4).name='coeur';
  out(4).type='m_hyper';
  out(4).unit='SI';

  i1=find(strcmpi(st,{out.name}));
  out1='Hyper_Elastic';

  if isempty(i1) && isempty(st); return; end

  % match a name 
  if ~isempty(i1); out=out(i1(1));

  else % assume values given
    error('Not a supported material');
  end


%% #BuildConstit (initialize integ/constit) ---------------------------------
% Implementation of elastic constitutive law building 3D
elseif comstr(Cam,'buildconstit')

  ID=varargin{carg};carg=carg+1;
  pl=varargin{carg};carg=carg+1; mat=pl(pl(:,1)==ID(1),:);
  il=varargin{carg};carg=carg+1; pro=[];try;pro=il(il(:,1)==ID(2),:);end
  if isempty(pl); error('MatId is not matched');end

   [st,unit,typ]=fe_mat('type',mat(2));
   
   if strcmp(st,'m_hyper')&&typ==1 % Nominal hyperelastic behaviour

    % set viscous coef to zero for static case
    out=[mat(3:7) 0 1];out=out(:); 
    
   elseif strcmp(st,'m_hyper')&&(typ==2||typ==3)
   
   %propagate constits, viscosity enabled for dynamic behavior
    out=mat(3:end); out=out(:);
   else
    error('Volume law not implemented for this material %s(%i)',st,typ);
   end 
   % Integration rule saved in integ % beware that ISOP=100 needed for large transform
   if length(ID)<4; ID(4)=0;end; ID(3)=3*ID(4); % must be number of DOFs 
   if size(pro,1)&&size(pro,2)>3;r3=pro(1,4:end);ID(6+(1:length(r3)))=r3; end  
   % used for later NL constit type selection
   if length(ID)<7||ID(7)==0; ID(7)=3; end
   % this should match StrategyType in of_mk.c matrix assembly
   out1=int32(ID(:)); out2=[]; 

%% #FieldDOFs ----------------------------------------------------------------
elseif comstr(Cam,'fielddofs')
 out=[1 2 3];
%% #SubTypeString ------------------------------------------------------------
elseif comstr(Cam,'subtypestring')

 i1=varargin{carg}; carg=carg+1;
 switch i1
 case 1;  out='xxx';
 otherwise; out='m_hyper';
 end
 
%% #Test ---------------------------------------------------------------------
elseif comstr(Cam,'test')
    
k1=4200;
k2=33;
kb=630000;
mu=2*(k1+k2);
lamb=kb-4/3*(k1+k2);
E=mu*(3*lamb+2*mu)/(lamb+mu);
nu=lamb/2/(lamb+mu);

model_lin=femesh('testhexa8');
model_lin.pl=[100,fe_mat('m_elastic','SI',1),E,nu,1000];
model_lin.il=p_solid('dbval 111 d3 4');
[Case_lin,model_lin.DOF]=fe_mknl('init',model_lin);
Klin = fe_mknl('assemble',model_lin,Case_lin,1);
Mlin = fe_mknl('assemble',model_lin,Case_lin,2);

model_hyp=femesh('testhexa8');
model_hyp.pl=[100,fe_mat('m_hyper','SI',2),1000,1,k1,k2,kb,0,0.001];
model_hyp.il=p_solid('dbval 111 d3 4');
[Case_hyp,model_hyp.DOF]=fe_mknl('init',model_hyp);
ofVar = struct('def',zeros(size(model_hyp.DOF,1),3),'DOF',model_hyp.DOF);
Mhyp = fe_mknl('assemble',model_hyp,Case_hyp,2);
Khyp=fe_mknl('assemble',model_hyp,Case_hyp,ofVar,5);


if (normest(Mlin-Mhyp)/normest(Mlin) > 1e-14) || ...
        (normest(Klin-Khyp)/normest(Klin) > 1e-14)
    error('Mismatch in hyperelastic matrices');
end
%% #TableCall -------------------------------------------------------------------------
elseif comstr(Cam,'tablecall');out='';
elseif comstr(Cam,'@');out=eval(CAM);
elseif comstr(Cam,'cvs')
 out='$Revision: 1.37 $  $Date: 2022/06/09 19:08:08 $'; return;
else; sdtw('''%s'' not known',CAM);
end
% -------------------------------------------------------------------------
end


%% #hyper : implementation of methods associated with hyperelastic constutive law
% May require SDT-nlsim license to work
%% #hypertoOpt [NL,r1,checkFu]=
function [out,out1,out2]=hypertoOpt(r1,mo1,C1) 
  NL=r1; 
 if nargin>1&&isfield(C1,'GroupInfo')
  %% #Large_deformation_multi-material case hyper_form=3  -3
  sdtw('_ewt','report calling case expected obsolete')
  dt=0;
  if size(C1.GroupInfo,1)>1; error('Not implemented');end
  EC=C1.GroupInfo{1,8};
  NL.opt=[double(comstr('hyper',-32));0;dt; %1:3
    C1.GroupInfo{1,4}(:)]; % constit is stored from 4 to end
  pointers=C1.GroupInfo{1,2}; integ=C1.GroupInfo{1,3};
  if length(unique(integ(4,:)))~=1; error('Not handled');end
  NL.iopt=int32([-1 -1 size(r1.X{1},1) size(r1.X{2},1)]);
  NL.iopt(7)=3; % hyperform=3; 
  % r1.Elt(2:end,5) is eltid;  
  eltid=feutil('eltid;',mo1); cEGI=find(isfinite(mo1.Elt(:,1)));
  nind=sparse(eltid(cEGI),1,1:length(cEGI));
  i1=pointers(7,nind(r1.Elt(2:end,5)));% MatOffSet for each gauss
  
  NL.iopt(9+(1:numel(i1)))=int32(i1(:));
  if NL.iopt(3)==9&&size(EC.ConstitTopology{1},1)==6
   ind_ts_eg=[1 6 5 6 2 4 5 4 3];
   i2=EC.ConstitTopology{1};i2=i2(ind_ts_eg,ind_ts_eg);
   NL.ddg=vhandle.matrix.stressCutDDG(struct('alloc',[9 NL.iopt(4)]));
  elseif NL.iopt(3)==size(EC.ConstitTopology{1},1)
    i2=EC.ConstitTopology{1};
  end
  NL.ConstitTopology{1}=int32(i2); 
  %% dd=cblas_dgthr(numel(i2),constit,dd,topo)
  NL.c=NL.cta; NL.b=[]; NL=feutil('rmfield',NL,'cta');
  NL=sdth.sfield('orderfirst',NL,{'unl','opt','pjg','vnl'});
  NL.unl=zeros(9,size(NL.c,1)/9);
 elseif ~isfield(NL,'old')
  %% 2021 revised version with uMax calls
  pl=[]; 
  try;
   if isfield(r1,'pl');pl=r1.pl;
   else
    if nargin==1; mo1=evalin('caller','mo1');end
    pl=mo1.pl(mo1.pl(:,1)==r1.Sens.Lambda{1,4}.ID(1),:);
   end
  end

   if isfield(NL,'kappa_v');pl(8)=r1.kappa_v; end % No relaxation cells
   if isfield(NL,'c3');pl(9)=r1.c3; elseif length(pl)<9; pl(9)=0; end %not Carrol
   if isfield(NL,'c1');pl(5)=r1.c1; elseif isfield(NL,'C1');pl(5)=r1.C1; end
   if isfield(NL,'c2');pl(6)=r1.c2; elseif isfield(NL,'C2');pl(6)=r1.C2;end
   if isfield(NL,'kappa');pl(7)=r1.kappa; elseif isfield(NL,'K');pl(7)=r1.K; end
   if ~isfield(NL,'g');r1.g=[];r1.f=[];  end % No relaxation cells
   ncell=size(r1.g,1); dt=0;
   NL.opt=[double(comstr('uMaxw',-32));0;dt]; %1:3
   if any(pl(4)==[1 0])
    % c1 c2 c3 kappa, kappa_v (dIdc)
    tcell=[300 0 0];
    cval=[pl([5:6 9 7 8])';ones(3,1);zeros(15,1)];
   elseif any(pl(4)==[2 3]) % 2Yeoh 3 params no c1h, 3Carrol
    tcell=[300 pl(4) 0];
    cval=[pl([5:6 9 7 8])';ones(18,1)];%dIdc at end
   else
     error('Not yet implemented')
   end
   if ~isempty(r1.g); tcell(3)=length(r1.g);
     cval=[cval;reshape([r1.g(:)';r1.f(:)'],[],1)];
   end
   cval=[cval(:);ones(18,1)]; %dIdc at end
   NL.opt(9+(1:length(cval)))=cval;
   NL.iopt(9+(1:length(tcell)+1))=int32([tcell -1]);
  if ~isfield(NL,'Sens');NL.iopt(5)=size(NL.unl,2); % Ngauss
  else; NL.iopt(5)=size(NL.Sens.Node,1);
  end
  NL.MexCb{2}=struct('unl',zeros(9*(ncell+2)+2), ... % single gauss here
      'opt',[], ... % Keep empty for pointer handling
      'jg',int32([0 1]));
  % iopt=[Find,iu,Nstrain,Nistate,Ngauss,StoreType(6),hyper_form(7)
  %       un(8), un(9) ostart(i=0;i<Ng)]
  NL.snllab={'PK1xx';'PK1xy';'PK1xz';'PK1yx';'PK1yy';'PK1yz';'PK1zx';'PK1zy';'PK1zz'};
  NL.snl=[]; 
  NL.adof=(0.5:.01:.58)'; 
  NL.unllab={'guxx';'guxy' ;'guxz' ;'guyx' ;'guyy' ;'guyz' ;'guzx' ;'guzy' ;'guzz';
                'un'; 'I_3' ; 
                'H1xx';'H1xy' ;'H1xz' ;'H1yx' ;'H1yy' ;'H1yz' ;'H1zx' ;'H1zy' ;'H1zz';
                'H2xx';'H2xy' ;'H2xz' ;'H2yx' ;'H2yy' ;'H2yz' ;'H2zx' ;'H2zy' ;'H2zz';
                'H3xx';'H3xy' ;'H3xz' ;'H3yx' ;'H3yy' ;'H3yz' ;'H3zx' ;'H3zy' ;'H3zz'};
  NL.udof=(.59:.01:.96)';  
  if NL.iopt(3)==9
   NL.ddg=vhandle.matrix.stressCutDDG(struct('alloc',[9 NL.iopt(4)]));
  end

 elseif 1==2
  %% Hyperelastic case   
  mu=r1.c1+r1.c2; lambda=r1.kappa-2*mu/3;
  if isfield(r1,'dt_fact') % if dt_fact given need estimate
   dt=r1.lmin/sqrt((lambda+2*mu)/r1.rho)/r1.dt_fact;
  elseif isfield(r1,'opt')&&length(r1.opt)>2&&r1.opt(3) 
    % If dt given nonzero don't estimate (forced displacement)
    dt=r1.opt(3);
  else
  end
  if ~isfield(NL,'kappa_v');r1.kappa_v=0; end % No relaxation cells
  if ~isfield(NL,'g');r1.g=[];r1.f=[]; dt=0; end % No relaxation cells
  ncell=size(r1.g,1); 
          
  NL.opt=[double(comstr('hyper',-32));0;dt; %1:3
    r1.c1;r1.c2;r1.kappa;r1.kappa_v;0; ... %c1(4),c2,kappa,kappa_v,un1(8)
    r1.g; % 8+(1:Ncell) maxwell fractions
    exp(-dt.*r1.f);exp(-dt./2*r1.f)]; %pre-computed exponentials
  NL.MexCb{2}=struct('unl',zeros(9*(ncell+2)+2,1), ...
      'opt',[], ... % Keep empty for pointer handling
      'jg',int32([0 1]));
  % iopt=[Find,iu,Nstrain,Nistate,Ngauss,StoreType(6),hyper_form(7)
  %       un(8), un(9) ostart(i=0;i<Ng)]
  NL.FNLlab={'PK1xx';'PK1xy';'PK1xz';'PK1yx';'PK1yy';'PK1yz';'PK1zx';'PK1zy';'PK1zz'};
  NL.adof=(0.5:.01:.58)'; 
  NL.unllab={'guxx';'guxy' ;'guxz' ;'guyx' ;'guyy' ;'guyz' ;'guzx' ;'guzy' ;'guzz';
                'un'; 'I_3' ; 
                'H1xx';'H1xy' ;'H1xz' ;'H1yx' ;'H1yy' ;'H1yz' ;'H1zx' ;'H1zy' ;'H1zz';
                'H2xx';'H2xy' ;'H2xz' ;'H2yx' ;'H2yy' ;'H2yz' ;'H2zx' ;'H2zy' ;'H2zz';
                'H3xx';'H3xy' ;'H3xz' ;'H3yx' ;'H3yy' ;'H3yz' ;'H3zx' ;'H3zy' ;'H3zz'};
  NL.udof=(.59:.01:.96)';  
 end
  out=NL; out1=r1; out2=@hypercheckFu;
end

function [out,out1,out2]=hypercheckFu(varargin)
 %% #hypercheckFu possibly change relaxation times 
 NL=varargin{1};mo1=varargin{2};
 opt=stack_get(mo1,'info','TimeOpt','g');
 if length(NL.iopt)>=10&&NL.iopt(10)==300 % Nov 21 update don't use DT
 elseif opt.Opt(4)~=NL.opt(3) 
  dt=opt.Opt(4);r1=NL; error('Obsolete needs recheck')
  NL.opt=[double(comstr('hyper',-32));0;dt;r1.c1;r1.c2;r1.kappa;r1.kappa_v;0; ...
    r1.g; %maxwell fractions
    exp(-dt.*r1.f);exp(-dt./2*r1.f)]; %pre-computed exponentials
 end
 NL.ddg=vhandle.matrix.stressCutDDG(struct('alloc',[9 NL.iopt(5)]));
 NL.JacFcn=@hyperJacobian;
 %NL.snl=zeros(NL.iopt(3)*NL.iopt(4),1)
 if size(NL.Sens.X{1},1)~=9; error('Observation problem use isop 100');end
 out=NL; 
end

%% #hypertoTable : generic SDT table to allow automated formatting
% feval(dfr_3Dnl('@hyper_toTable'),mo1b.NL{end},mo1b)
% sdtweb dfr_ident LoadOneTest.simo
function [out,out1]=hypertoTable(NL,mo1)
  if isfield(NL,'plot')||isfield(NL,'table');
   %% Display Maxwell cell impact
   RO=NL; opt=stack_get(mo1,'','TimeOpt','g');
   r1=stack_get(mo1,'pro');
   for j1=1:size(r1,1)
    NL=r1{j1,3}; if ~isfield(NL,'NLdata')||~isfield(NL.NLdata,'MexCb');continue;end
    NL=NL.NLdata;
    ncell=(length(NL.MexCb{2}.unl)-2)/9-2;
    if isempty(opt);dt=NL.opt(3); else; dt=opt.Opt(4); end
    r2=struct('c1',NL.opt(4),'c2',NL.opt(5), ...
        'kappa',NL.opt(6),'kappa_v',NL.opt(7),'g',NL.opt(8+(1:ncell)), ...
        'f',-log(NL.opt(8+ncell+(1:ncell)))/dt);   
    if isfield(RO,'table')
     disp(comstr(r2,-30,struct('NoClip',1)))
    else; r2.cf=51;feval(nlutil('@Maxwell'),'view',r2);
    end
   end
   
  else
   labels=fe_mat('unitlabel',mo1.unit);
   out={'opt'   ,'  '  , ' ';
       'c1'   ,'1st Mooney Rivlin parameter' , labels{1} ;
       'c2'   ,'2nd Mooney Rivlin parameter' , labels{1};
       'kappa'   ,'Bulk modulus',  labels{1};
       'kappav'   ,'Bulk viscosity (numeric)',  labels{1};
       'un1','unused 1',labels{1};
       'g_i'   ,'Viscous fraction for each cell', '[]';
       'et_i'   ,'exp(-dt/tau_i)', '[]';
       'e2t_i'   ,'exp(-dt/2/tau_i)', '[]';
       'u_int'   ,' ',   '   ';
       'Snl1'   ,'Stress at last time-step' ,labels{1};
       'Hi'   ,'Viscous stress for each cell',labels{1} };
  end
  
end
%% #hyperJacobian implement jacobian -2
% see openfem/mex/hyper.c 
function [kj,cj]=hyperJacobian(varargin) 
 % [kj2,cj2]=feval(r2.JacFcn,r2,[],model,u,v,[],opt,Case,RunOpt);
 kj=[];cj=[]; NL=varargin{1}; 
 r2=struct('unl',NL.unl,'opt',NL.opt,'iopt',NL.iopt, ...
     'wjdet',NL.Sens.wjdet,'ddg',NL.ddg);
 ja=mkl_utils(struct('jac',1,'simo',r2));
 c=NL.c.GetData; 
 i2=reshape(1:size(c,1),[],NL.iopt(5));i2=i2(1:NL.iopt(3),:);
 if 1==2
  ci_ts_eg=[1 5 9 8 7 4];
  for j1=1:size(NL.ddg,1)/9
   dd=full(NL.ddg((j1-1)*9+ci_ts_eg,(j1-1)*9+ci_ts_eg));d=eig(dd);
   if any(d<0);dbstack; keyboard;end
  end
 end
 kj=feutilb('tkt',c(i2,:),NL.ddg);
end

%% #hyper_g  Single gauss point evaluation of hyper
function [out,out1]=hyper_g(RO)

  unl=RO.unl(:);
  F=unl(1:9)+[1;0;0;0;1;0;0;0;1]; % u_i,j+eye(3)
  Sn=reshape(unl(12:20),3,3);
  Hn=reshape(unl(21:end),3,3,[]);  

  ncell=(length(RO.opt)-8)/3;  
  g=RO.opt(9:9+ncell-1); et=RO.opt(9+ncell:9+2*ncell-1); 
  e2t=RO.opt(9+2*ncell:end); 
  F=(reshape(F,3,3)); 

  % Cauchy Green tensor
  C = F'*F; I3=det(F); %change for direct expression
  
  %if I3<0; out=F; out(1)=NaN; out1=zeros(size(unl,1)-9,1); return; end
  if I3<0; error('Negative I3 problem');end
  J23=I3^(-2/3);  iF=inv(F); J23Ft=J23*F'; iFt=iF';
  J23C=J23*C; I1b=J23C(1)+J23C(5)+J23C(9); 
  
  %% deviatoric hyperelastic stress
  c1=RO.opt(4); c2=RO.opt(5); 
  S_0 = 2*((c1+c2*I1b)*eye(3)-c2*(J23C)); %  dW=2*(c1*eye(3)+c2*(I1b*eye(3)-J23C));

  tau_0 = F*S_0*J23Ft; % Hyperelastic Kirchoff
  tau0b = tau_0-(tau_0(1)+tau_0(5)+tau_0(9))/3*eye(3); % Deviatoric part of tau_0

  out1=zeros(9,ncell+1);
  Snp1=iF*(tau0b/J23)*iFt;  %#ok<MINV> Piola Kirchoff 2 
  out1(:,1)=Snp1(:);
  hb=zeros(3,3); %ht=0;
  
  for jcell=1:ncell % Loop on viscous cells
    Htn = et(jcell)*Hn(:,:,jcell) - e2t(jcell)*Sn;
    Hnp1 = Htn + e2t(jcell) *(Snp1); 
    hp = F*Htn*J23Ft; 
    %   ht= ht+g(jcell)*trace(hp);
    hb = hb+ g(jcell)*hp-(g(jcell)*(hp(1)+hp(5)+hp(9))/3)*eye(3);
    out1(:,jcell+1)=Hnp1(:);
  end
  % Sum of cell fractions
  g_full = 1 + sum(g.*(e2t-1));
 
  % Compressibility p = Uo'(I3); % pression
  p0=RO.opt(6)/2*((I3)-1/I3);  
  I3n=unl(10); % (I3-1)_n invariant at last time step
  pv=RO.opt(7)*(I3-1-I3n); % kappa_v*(I3-1-I3n); damped pressures
  
  %Output management
  tau = I3*(pv+p0)*eye(3) + g_full*tau0b + hb; 
  FS = tau*(iFt); % Piola Kirchhoff 1 = FS
  
  out=FS(:); % snl_t+1^k+1
  un=unl(11); % currently unused state
  out1=[I3-1; un; out1(:)]; %out1(end+1)=p0; % unl_t+1^k+1
 %if RO.jg(1)==0; 
 %    disp([out([1])' out1([1 3+[0 9 18]])']);
 %end
 
end
%% #bv_hyper_g  Single gauss point evaluation of hyper -2
function [out,out1]=bv_hyper_g(RO)

  unl=RO.unl;
  F=reshape(unl(1:9)+[1;0;0;0;1;0;0;0;1],3,3); 
  Sn=reshape(unl(10:18),3,3);
  Hn=reshape(unl(19:end-2),3,3,[]);

  ncell=(length(RO.opt)-8)/3;  
  c1=RO.opt(4); c2=RO.opt(5); kappa=RO.opt(6);
  g=RO.opt(7:7+ncell-1); et=RO.opt(7+ncell:7+2*ncell-1); 
  e2t=RO.opt(7+2*ncell:end-2);
  etp=RO.opt(end-1); e2tp=RO.opt(end-2);
  g_inf=1-sum(g);

  C = F'*F; I3=sqrt(det(C)); 
  J23=I3^(-2/3);  iF=inv(F); J23Ft=J23*F';iFt=iF';
  J23C=J23*C; I1b=trace(J23C); 
  %deviatoric hyperelastic stress
  S_0 = 2*((c1+c2*I1b)*eye(3)-c2*(J23C)); %  dW=2*(c1*eye(3)+c2*(I1b*eye(3)-J23C));
  tau_0 = F*S_0*J23Ft; 
  tau0b = tau_0-trace(tau_0)/3*eye(3); 
  %viscous stress
  out1=zeros(9,ncell+1);
  Snp1=iF*(tau0b/J23)*iFt;  %#ok<MINV>
  out1(:,1)=Snp1(:);
  hb=zeros(3,3); %ht=0;
  
  for jcell=1:ncell 
    Htn = et(jcell)*Hn(:,:,jcell) - e2t(jcell)*Sn;
    Hnp1 = Htn + e2t(jcell) *(Snp1); 
    hp = F*Htn*J23Ft; 
    hb = hb+ g(jcell)*hp-(trace(g(jcell)*hp)/3)*eye(3);
    out1(:,jcell+1)=Hnp1(:);
  end
  g_full = g_inf + sum(g.*e2t);
 
  % Compressibility 
  % p = Uo'(I3); % pression
  % bulk visco
   p0=kappa/2*((I3)-(1));  pvn=unl(end);  
   pv= e2tp*p0-e2tp*pvn;%
  
  out1=out1(:); out1(end+1)=00; out1(end+1)=p0;
  
  %Output management
  tau = I3*pv*eye(3) + g_full*tau0b + hb; 
  FS = tau*(iFt); % Piola Kirchhoff 1 = FS
   % out = inv(F)*out; %PK2
  out=FS(:); % snl_t+1^k+1
  out1=out1(:); % unl_t+1^k+1

end

function [out,out1]=checkNL(RO,mo1)
%% #checkNL non regression tests. t_nlspring volmathyper

 if ischar(RO)
  if     strcmpi(RO,'funiaxial'); out=@(l)diag([l 1/sqrt(l) 1/sqrt(l)]);
  elseif strcmpi(RO,'fequibiaxial'); out=@(l)diag([l l l^-2]);
  elseif strcmpi(RO,'fshear'); out=@(l)[1 l-1 0;0 1 0;0 0 1];
  elseif strcmpi(RO,'fiso'); out=@(l)[l 0 0;0 l 0;0 0 l];
  elseif strcmpi(RO,'viewa')
  %% #checkNL.viewA : display stress for standard range of values 
   NL=mo1; 
   if isfield(NL,'li');
     RO=NL; NL=struct;
     for j1=1:length(RO.li)
       r2=m_hyper('urn',RO.li{j1}); RO.pl(j1,1:length(r2.pl))=r2.pl;
     end
   else; RO=NL; 
   end
   if isfield(RO,'ld');ld=RO.ld; 
   else; ld=unique([linspace(.2,2,100)';1+logspace(-3,-1,10)';1-logspace(-3,-1,10)']);
   end
   RO.to={checkNL('funiaxial');checkNL('fshear');checkNL('fequibiaxial');checkNL('fiso')};
   %RO.icomp={@(x)x(1,1)-x(2,2),@(x)x(1,2), @(x)x(1,1)-x(3,3),@(x)trace(x)};%See Rafael
   RO.icomp={@(x)x(1,1),@(x)x(1,2), @(x)x(1,1),@(x)trace(x)};%See Rafael
   out=struct('X',{{ld-1,{'uniaxial','shear','equibiaxial','iso'}}},'Y',zeros(size(ld)));
   RO.plot={};
   for j1=1:size(RO.pl,1) % Loop on materials to compare
    NL.pl=RO.pl(j1,:);
    for jpar=1:length(ld)
     for j2=1:length(RO.to)
      NL.unl=reshape(RO.to{j2}(ld(jpar))-eye(3),[],1);
      if jpar==1&&j2==1
          NL=feval(m_hyper('@hypertoOpt'),NL); 
          if isfield(RO,'check2');NL.check2=1;end
      end
      NL.silent=1; r2=checkNL(NL); if ~r2.stable;r2.Sigma(:)=NaN;end
      NL=feutil('rmfield',NL,'check2');
      out.Y(jpar,j2)=RO.icomp{j2}(r2.Sigma);
     end
    end
    RO.plot(1,end+(1:2))={out.X{1},out.Y};
    % display assymptotic values
    %% #check.NL.Consistence around 0 
    NL.unl(:)=0;r2=checkNL(NL); 
    if NL.iopt(11)==2
     G=2*(NL.opt(10));RO.mtype='Yeoh'; % mu=G=C10
    else % Money Rivlin
     G=2*(NL.opt(10)+NL.opt(11)); RO.mtype='Moon'; % mu=G=C1+C2
    end
    kappa=NL.opt(13); nu=(1-G/kappa)/2; E=2*G*(1+nu);
    % block kappa*eye(3) = G*3/(1-2*nu)
    dd=m_elastic('formulaENG2DD',[E nu G]);
    r3=reshape(dd(1:3,1:3),[],1)\reshape(r2.dd(1:3,1:3),[],1);
    nu=.5-((.5-nu)/r3); % Correct nu to match better
    fprintf('%s E=%.5g, nu=%.6f, G=%.5g, kappa=%.5g (top m_elastic, bot: m_hyper)\n',RO.mtype,E,nu,G,kappa)
    dd=m_elastic('formulaENG2DD',[E nu G]);
    
    disp([dd;r2.dd])
    1;
   end
   gf=10; figure(gf)
   out1=plot(RO.plot{:});out1=reshape(out1,[],size(RO.pl,1));
   if isfield(RO,'os'); sdth.os(gf,RO.os{:});end
   legend(out.X{2},'location','best'); xlabel('Strain');ylabel('Stress')
   grid on
   
  end
  return

 elseif ~isfield(RO,'unl') 
  %% Possibly init from RO
  F=RO.F;
  RA=RO;
  if isfield(RO,'NLdata')
   RO=RO.NLdata; RO.unl=zeros(9,1,3); if isfield(RA,'pl');RO.pl=RA.pl;end
  elseif ~isfield(RO,'unl'); RO.unl=reshape(RO.F-eye(3),[],1); end
  RO=feval(m_hyper('@hypertoOpt'),RO); %NL.opt
 else
  F=reshape(RO.unl(1:9),3,3)+eye(3); 
 end
 %% #checkNL.compute 
 out=struct('F',F);
 [C,I,dIdc,d2I3dcdc,d2I2dcdc]=feval(elem0('@elemCalc'),F); 
 ci_ts_eg=[1 5 9 8 7 4];ci_ts_egt=ci_ts_eg';dIdc(ci_ts_eg,:);
 if isfield(RO,'iopt')&&RO.iopt(10)==300
 %% see rafael (2.19) to (2.21) [dWdI,d2WdI2]
  if any(RO.iopt(11)==[0 1]);    constit=[double(RO.iopt(11)) RO.opt([10 11 13])'];
  elseif any(RO.iopt(11)==[2 3]);constit=[double(RO.iopt(11)) RO.opt([10 11 12 13])'];
  else; error('M check not implemented')
  end
  [dWdI,d2WdI2]=feval(m_hyper('@EnHyper'),[],constit,I);% WithLog
  if isfield(RO,'check2') % check tangent stiff
   dx=1e-7;RO.I=I;
   r2=zeros(3);I=RO.I; [dWdI,d2WdI2]=feval(m_hyper('@EnHyper'),[],constit,I);RO.d2WdI2=d2WdI2;RO.dWdI=dWdI;
   for j1=1:3
     I=RO.I;I(j1)=I(j1)+dx; 
     [r3,d2WdI2]=feval(m_hyper('@EnHyper'),[],constit,I);
     r2(:,j1)=(r3(:)-RO.dWdI(:))/dx;
   end
   fprintf('d2WdI2 derivative check\n');
   disp([r2 d2WdI2]);
  end
 else % Used for OpenFEM tests
  [dWdI,d2WdI2]=feval(m_hyper('@EnHyper'),[],[1 RO.opt(4:end)'],I);% WithLog
 end
 out.Sigma=reshape(2*dIdc*dWdI(:),3,3);
 %dIdc(ci_ts_eg,:)*d2WdI2, hyperelasticity Chapelle (40) there is a factor 4
 %4*[reshape(dWdI(2)*d2I2dcdc+dWdI(3)*d2I3dcdc,[],1) ...
 %   reshape(dIdc(ci_ts_eg,:)*d2WdI2*dIdc(ci_ts_eg,:)',[],1)];
 %ans(1:6,:)
 % Rafael (2.16)
 out.d2wde2=4*(dWdI(2)*d2I2dcdc+dWdI(3)*d2I3dcdc) + 4*dIdc(ci_ts_eg,:)*d2WdI2*dIdc(ci_ts_eg,:)'; 
 out.stable=1; 
 out.d2WdI2=d2WdI2;out.dWdI=dWdI(:)';out.I=I; out.dIdc=dIdc; 
 % Rafael (4.10), SDT sdtweb feform  % (6.54)
 out.dd9=feval(elem0('@LdDD'),out.F,out.d2wde2,out.Sigma);
 out.dd=out.dd9(ci_ts_eg,ci_ts_eg); out.ci_ts_eg=ci_ts_eg;

 [x,d]=eig(out.dd);d=diag(d);
 if any(d<0) 
     out.stable=0;
     if ~isfield(RO,'silent')
       i1=find(d<0); disp(d(i1)');disp(x(:,i1))
       warning('dd <=0 stability problem');
     end
     [x,d]=eig(d2WdI2);d=diag(d); 
     if any(d<0)&&~isfield(RO,'silent') 
      warning('d2WdI2 <=0 stability problem');
      disp(d');disp(x)
     end
 end

end


%% #EnHyper hyperelastic models m file checks --------------------------------
function [dWdI,d2WdI2]=EnHyper(integ,constit,I) %#ok<INUSL>
%[dWdI,d2WdI2]=feval(elem0('@EnHeart'),[],[],I);
%C1=0.3MPa, C2=0.2MPa, K=0.3MPa

% Cenerg : C1 C2 K
%if length(constit)<3;constit=[25. 0.25 1274.];end
%if length(constit)<3;constit=[0 1 3 2 3];end
if length(constit)<3;constit=[0 .3 .2 .3];end
if ~any(constit(2:end)); error('No constit');end

if constit(1)==0
% ----------------------------------
%% hperelastic type 0: C1*(J1-3)+C2*(J2-3)+K*(J3-1)^2%
kappa=constit(4); ktyp=0;
dWdI(1) = constit(2)*I(3)^(-1./3.);
dWdI(2) = constit(3)*I(3)^(-2./3.);
dWdI(3) = - 1./3.* constit(2)*I(1)*I(3)^(-4./3.) ...
         - 2./3.* constit(3)*I(2)*I(3)^(-5./3.) ;
d2WdI2=[0 0 -1./3.*constit(2)*I(3)^(-4./3.) ;
       0 0  -2./3.*constit(3)*I(3)^(-5./3.);
      -1./3.*constit(2)*I(3)^(-4./3.)  -2./3.*constit(3)*I(3)^(-5./3.) ...
       4/9*constit(2)*I(1)*I(3)^(-7./3.)...
             + 10./9.* constit(3)*I(2)*I(3)^(-8./3.)];
elseif constit(1)==1.1
%% #he_MooneyRivlin_CiarletGeymonat derived from Rafael (kappa wrong)
  c1=constit(2); c2=constit(3); kappa=constit(4); ktyp=0;
dWdI=[c1*I(3)^(-1/3);
    c2*I(3)^(-2/3);
    -c1*I(1)*I(3)^(-4/3)/3-2*c2*I(2)*I(3)^(-5/3)/3];
 %[-c1*I(1)*I(3)^(-4/3)/3 -2*c2*I(2)*I(3)^(-5/3)/3 kappa/2*(1-I(3)^(-1))]
  d2WdI2=[0,0,-c1*I(3)^(-4/3)/3 ; 
   0,0,-2*c2*(I(3)^-(5/3))/3;
   -c1*I(3)^(-4/3)/3, -2*c2*(I(3)^-(5/3))/3, 4*c1*I(1)*(I(3)^(-7/3))/9+...
   10*c2*I(2)*(I(3)^(-8/3))/9];


elseif constit(1)==1
% ----------------------------------------------------------
%% hyperelastic type 1: C1*(J1-3)+C2*(J2-3)+K*(J3-1-ln(J3)) 
% J3=I3^(-1/2) p=k (1-J)/J ref mex/hyper.c line 138
kappa=constit(4); ktyp=1;
dWdI(1) = constit(2)*I(3)^(-1./3.);
dWdI(2) = constit(3)*I(3)^(-2./3.);
dWdI(3) = - 1./3.* constit(2)*I(1)*I(3)^(-4./3.) ...
          - 2./3.* constit(3)*I(2)*I(3)^(-5./3.);  

d2WdI2=[0 0 -1./3.*constit(2)*I(3)^(-4./3.) ;
        0 0  -2./3.*constit(3)*I(3)^(-5./3.);
       -1./3.*constit(2)*I(3)^(-4./3.)  -2./3.*constit(3)*I(3)^(-5./3.) ...
        (4./9.*constit(2)*I(1)*I(3)^(-7./3.)...
              + 10./9.* constit(3)*I(2)*I(3)^(-8./3.))];
elseif constit(1)==2
 %% Yeoh C10 C20 C30 kappa kappav 
  c1=constit(2); c2=constit(3); c3=constit(4); kappa=constit(5); 
 va=I(1)*I(3)^(-1/3)-3; vb= (c1+2*c2*va+3*c3*va*va);va=(2*c2+6*c3*va);
 dWdI=vb*[I(3)^(-1/3) 0 -1/3*I(1)*I(3)^(-4/3)];
 d2WdI2=va*[I(3)^(-2/3) 0 -1/3*I(1)*I(3)^(-5/3);0 0 0;
     -1/3*I(1)*I(3)^(-5/3) 0 1/9*I(1)^2*I(3)^(-8/3)] + ...
     vb*[0 0 -1/3*I(3)^(-4/3); 0 0 0; -1/3*I(3)^(-4/3) 0 4/9*I(1)*I(3)^(-7/3)];
 %[va vb  I(1)^2*I(3)^(-8/3) I(1)*I(3)^(-7/3)  d2WdI2(end)]
 %[va*[I(3)^(-2/3) -1/3*I(1)*I(3)^(-5/3)]  vb*-1/3*I(3)^(-4/3) d2WdI2(1,[1 3])]

 ktyp=1;
elseif constit(1)==3
  c1=constit(2); c2=constit(3); c3=constit(4); kappa=constit(5); 
  dWdI=[c1*I(3)^(-1/3)+4*c2*I(1)^3*I(3)^(-4/3);
    c3*I(2)^(-1/2)*I(3)^(-1/3)/2;
    -c1*I(1)*I(3)^(-4/3)/3-4*c2*(I(1)^4)*I(3)^(-7/3)/3-c3*I(2)^(1/2)*I(3)^(-4/3)/3];
 
  d2WdI2=[12*c2*I(1)^2*I(3)^(-4/3),  0 ,-c1*I(3)^(-4/3)/3-16*c2*(I(1)^3)*I(3)^(-7/3)/3 ; 
   0,-c3*I(2)^(-3/2)*I(3)^(-1/3)/4,-c3*I(2)^(-1/2)*(I(3)^(-4/3))/6;
   -c1*I(3)^(-4/3)/3-16*c2*(I(1)^3)*I(3)^(-7/3)/3, -c3*I(2)^(-1/2)*(I(3)^(-4/3))/6, ...
   4*c1*I(1)*(I(3)^(-7/3))/9+ 28*c2*(I(1)^4)*I(3)^(-10/3)/9+4*c3*I(2)^(1/2)*(I(3)^(-7/3))/9];

end
%% #EnHyper.Compression (verified June 22)
if ktyp==1 % Ciarlet Gemonat
 dWdI(3)=dWdI(3)+0.5*kappa*(I(3)^(-1./2.)-I(3)^(-1));
 d2WdI2(9)=d2WdI2(9)-1./4.*kappa*I(3)^(-3./2.)+1/2*kappa*I(3)^(-2);
else % k/2 (I3^(1/2)-1)^2 
 dWdI(3)=dWdI(3)+kappa*(1-I(3)^(-1/2))/2;
 d2WdI2(9)=d2WdI2(9)+ kappa/4*I(3)^(-3/2);
end
1;
end
