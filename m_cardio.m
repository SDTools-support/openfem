function [out,out1,out2]=m_cardio(varargin)

%M_CARDIO myocardium contraction material function
%
%       Syntax : mat= m_elastic('default') 
%                mat= m_elastic('database name') 
%                pl = m_elastic('dbval MatId name');
%                pl = m_elastic('dbval -unit MM MatId name');
%                pl = m_elastic('dbval -punit MM MatId name');
%
%       Material subtypes supported by m_hyper are
%
%	Subtype 1 : xxx
%	    [MatId type rho Es mu k_0 s_0 dt alpha]
%         with
%           type  fe_mat('m_cardio','SI',1)
%       based on Hill-Maxwell model
%       Es: Young modulus associated to the series element 
%       mu: viscous coef related to the strain in contractile element
%       k_0: maximum stiffness for the contractile element
%       s_0: maximum stress for the contractile element
%       dt: time step
%       alpha: action potential factor
%	    rho xxx Es xxx mu xxx k_0 xxx s_0 xxx dt xxxx alpha xxxxx
%
%       See sdtweb      m_elastic, pl, fem
%       See also help   fe_mat, p_shell, p_beam


%	Etienne Balmes, Jean-Michel Leclere, Corine Florens, Mathieu Alba
%       Copyright (c) 2001-2012 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

if nargin<1; help m_cardio;return; end
if ischar(varargin{1}); [CAM,Cam]=comstr(varargin{1},1); pl=[]; carg=2;
else pl=varargin{1};[CAM,Cam]=comstr(varargin{2},1); carg=3;
end
if sp_util('diag');fprintf('m_cardio %s \n',CAM);end

%% #PropertyUnitType ---------------------------------------------------------
if comstr(Cam,'propertyunittype')

if nargin<2;out=1;return;end
i1=varargin{carg};
 out1={};
 switch i1
 %case 1 % isotropic   [MatId typ E nu rho G eta alpha T0]
 %               ind = [0     0   1 0  3   1 0   8     7]; % xxx 
 case 1
 st=...
 {'MatId'  0  'help m_cardio';
  'Type'   0  '';
  'C1'     0  '';
  'C2'     0  '';
  'K'      0  '';
 };
 otherwise; st={'MatId' 0 'help m_cardio'; 'Type', 0, ''};
 end
 if ~isempty(strfind(Cam,'cell')); out=st; else; out=[st{:,2}]; end


%% #Info ---------------------------------------------------------------------
elseif comstr(Cam,'info')

 r1=m_cardio('database');
 disp(char({r1.name}))

%% #Default ------------------------------------------------------------------
elseif comstr(Cam,'default')

  out=m_cardio('database'); out=out(1);

%% #DBval --------------------------------------------------------------------
elseif comstr(Cam,'dbval') 

 while 1==1
  i1=strfind(comstr(Cam,-27),'-unit'); out1={}; % See if unit specified
  if ~isempty(i1)
   [Unit,i2,i3,i4]=sscanf(CAM(i1+5:end),'%s',1);
   i4=i1+[0:4+i4];CAM(i4)=''; [CAM,Cam]=comstr(CAM,1);
  else Unit='';
  end
  i2=strfind(comstr(Cam,-27),'-punit');
  if ~isempty(i2)
   [PUnit,i3,i4,i5]=sscanf(CAM(i2+6:end),'%s',1);
   i5=i2+[0:5+i5];CAM(i5)=''; [CAM,Cam]=comstr(CAM,1);
  else PUnit='';
  end
  if ischar(CAM); [i1,CAM,Cam]=comstr(CAM,'dbval','%i');else; i1=[];end
  if isempty(CAM)&&carg<=nargin; st=varargin{carg};carg=carg+1;
  else st=CAM;end
  if isempty(st)
  elseif ischar(st); mat=m_cardio('database',st);
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
  pl(i2,1:length(r1))=r1;
  out1(end+1,1:3)={'mat',sprintf('%i_%s',mat.pl(1),mat.name),mat};
  if carg>nargin; break;end
  CAM=varargin{carg};carg=carg+1;if ischar(CAM);[CAM,Cam]=comstr(CAM,1);end
 end
 out=pl;
 
%% #Database -----------------------------------------------------------------
elseif comstr(Cam,'database'); [st,Cam]=comstr(CAM,9);

  if isempty(st)&&carg<=nargin; st=varargin{carg}; carg=carg+1;end
  [MatId,i2,i3,i4]=sscanf(st,'%i',1); if i2~=1; MatId=1;end
  st=comstr(st,i4);

  out.pl=[MatId fe_mat('type','m_cardio','SI',1) .3 1e-3 0.1 0.2 1e-3 1];
  out.name='Ref';
  out.type='m_cardio';
  out.unit='SI';

  out(2).pl=[MatId fe_mat('type','m_cardio','SI',1) 2e1 1e-2 0.4 0.6 5e-3 1];
  out(2).name='Rivlin';
  out(2).type='m_cardio';
  out(2).unit='SI';
  
  out(3).pl=[MatId fe_mat('type','m_cardio','SI',1) 5e6 3e3 1.5e5 2.7e5 1e-3 0.8];
  out(3).name='phi';
  out(3).type='m_cardio';
  out(3).unit='SI';
  

  out(4).pl=[MatId fe_mat('type','m_cardio','SI',1) 5e3 2e-1 150 300 1e-2 1];
  out(4).name='coeur';
  out(4).type='m_cardio';
  out(4).unit='SI';


  i1=strmatch(comstr(st,-27),comstr({out.name},-27),'exact');
  out1='Cardio_contract';

  if isempty(i1) && isempty(st); return; end

  % match a name 
  if ~isempty(i1); out=out(i1(1));

  else % assume values given
    error('Not a supported material');
  end


%% #buildConstit -------------------------------------------------------------
% Implementation of elastic constitutive law building 3D
elseif comstr(Cam,'buildconstit')

  ID=varargin{carg};carg=carg+1;
  pl=varargin{carg};carg=carg+1; mat=pl((pl(:,1)==ID(1)),:);
  il=varargin{carg};carg=carg+1; pro=[];try;pro=il((il(:,1)==ID(2)),:);end
  if isempty(pl); error('MatId is not matched');end

   [st,unit,typ]=fe_mat('type',mat(2));
   
   if strcmp(st,'m_cardio')&&typ==1 % Nominal hyperelastic behaviour

    % for now just propagate constants to a column of constit
    out=[mat(3:end)];out=out(:); 
    
   elseif strcmp(st,'m_cardio')&&typ==2 %not used
       out=[mat(3:end)];out=out(:); 
   else
    error('Volume law not implemented for this material %s(%i)',st,typ);
   end 
   % Integration rule saved in integ
   if length(ID)<4; ID(4)=0;end; ID(3)=3*ID(4); % must be number of DOFs 
   if size(pro,1)&&size(pro,2)>3;  ID(5)=pro(1,4); end  
   ID(7)=901; % used for later NL constit type selection
   % this should match StrategyType in of_mk.c matrix assembly
   out1=int32(ID(:));

%% #FieldDOFs ----------------------------------------------------------------
elseif comstr(Cam,'fielddofs')
 out=[1 2 3];
%% #SubTypeString ------------------------------------------------------------
elseif comstr(Cam,'subtypestring')

 i1=varargin{carg}; carg=carg+1;
 switch i1
 case 1;  out='xxx';
 otherwise; out='m_cardio';
 end
%% #TableCall ----------------------------------------------------------------
elseif comstr(Cam,'tablecall'); out='';
elseif comstr(Cam,'cvs')
 out='$Revision: 1.13 $  $Date: 2012/07/09 07:03:29 $';
else; sdtw('''%s'' not known',CAM)
 
end
% -------------------------------------------------------------------------
