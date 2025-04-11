function [out,out1,out2]=m_poro(varargin)

%M_PORO Poroelastic material function
%
%       Syntax : mat= m_poro('default') 
%                mat= m_poro('database name') 
%                pl = m_poro('dbval MatId name');
%                pl = m_poro('dbval -unit MM MatId name');
%                pl = m_poro('dbval -punit MM MatId name');
%
%  Material subtypes supported by m_poro are
%
%	1 : basic
%     [MatId fe_mat('m_poro','SI',1) ElasMatId   ...
%         phip alpha sigmaR lambda lambda2 rho0 rhoS eta P0 Npr qq ]
%    For more details use  m_poro('propertyunittypecell',1)
%
%       See <a href="matlab: sdtweb _taglist m_poro">TagList</a>
%       See sdtweb      m_poro, pl, fem
%       See also help   fe_mat, p_shell, p_beam


%	Etienne Balmes, Jean-Michel Leclere, Corine Florens
%       Copyright (c) 2001-2025 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC>

if nargin<1; help m_poro;return; end
if ischar(varargin{1}); [CAM,Cam]=comstr(varargin{1},1); pl=[]; carg=2;
else; pl=varargin{1};[CAM,Cam]=comstr(varargin{2},1); carg=3;
end

%% #Info ---------------------------------------------------------------------
if comstr(Cam,'info')

 r1=m_poro('database');
 disp(char({r1.name}))

%% #Default ------------------------------------------------------------------
elseif comstr(Cam,'default')
    
 model=[];if carg<=nargin;model=varargin{carg};carg=carg+1;end
 if isempty(model); out=m_poro('database'); out=out(1);
 else;              out=fe_mat('defaultpl',model);
 end
 
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
  else; st=CAM;end
  if isempty(st)
  elseif ischar(st); mat=m_poro('database',st); 
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

%% #Database -----------------------------------------------------------------
elseif comstr(Cam,'database'); [st,Cam]=comstr(CAM,9);

  if isempty(st)&&carg<=nargin; st=varargin{carg}; carg=carg+1;end
  [MatId,i2,i3,i4]=sscanf(st,'%i',1); if i2~=1; MatId=1;end
  st=comstr(st,i4);
 %  Thesis of O. Tanneau page 48 
  out.pl=[MatId fe_mat('type','m_poro','SI',1) ... 
           1  .98 1.03 6600  200e-6 380e-6  1.255 11.2   .06 101320 .71 1.4];
  out.name='Ref';
  out.type='m_poro';
  out.unit='SI';

  i1=find(strcmpi(st,{out.name}));
  out1='Poro';

  if isempty(i1) && isempty(st); return; end

  % match a name 
  if ~isempty(i1); out=out(i1);
  else % assume values given
    error('%s Not a supported material',CAM);
  end


%% #BuildConstit -------------------------------------------------------------
elseif comstr(Cam,'buildconstit');[out,out1,out2]=p_solid(varargin{:});

%% #PropertyUnitType ---------------------------------------------------------
elseif comstr(Cam,'propertyunittype')

if nargin==1;out=1; return;end
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
  'Elas'     0  'Elastic MatId';
  'phip'       0  'Porosity (no unit)';
  'alpha'       0  'Tortuosity (no unit)';
  'sigmaR'       0  'Resistivity [Ns/m4]';
  'lambda'       0  'Characteristic length [m]';
  'lambda2'       0  'Characteristic length [m]';
  'rho0'       0  'Air density [kg/m3]';
  'rhoS'       0  'Skeletton density [kg/m3]';
  'eta'       0  'Dynamic viscosity [kg/s/m]';
  'P0'       0  'Atmospheric pressure [N/m2]';
  'Npr'       0  'Prandt number [none]';
  'qq'       0  'Ratio of specific heat [none]';
  };

otherwise; st={'MatId' 0 'sdtweb(''m_elastic'')'; 'Type', 0, ''};
end

if ~isempty(strfind(Cam,'cell')); out=st; else; out=[st{:,2}]; end

%% #SubType ------------------------------------------------------------------
elseif comstr(Cam,'subtype');[CAM,Cam]=comstr(CAM,8);

st={'1 Basic'}; 
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

%% #FormulaEng2dd
if comstr(Cam,'eng2dd')

    keyboard

% -------------------------------------------------------------------------
else; sdtw('''Formula%s'' not known',CAM);
end
%% #end -------------------------------------------------------------------------
elseif comstr(Cam,'tablecall');out='';
elseif comstr(Cam,'cvs')
    out='$Revision: 1.6 $  $Date: 2025/04/07 17:08:05 $';
else; sdtw('''%s'' not known',CAM);
end % commands

