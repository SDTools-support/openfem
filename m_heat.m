function [out,out1,out2,out3]=m_heat(varargin)

% M_HEAT  thermic material function
%
%
%       Syntax : mat= m_heat('default') 
%                mat= m_heat('database name') 
%                pl = m_heat('dbval MatId name');
%                pl = m_heat('dbval -unit MM MatId name');
%                pl = m_heat('dbval -punit MM MatId name');
%                 
%       Material subtypes supported by m_heat are
%       Subtype 1 : 3D solid volume element: hexa8, hexa20, penta6, tetra4, ...
%            [MatId fe_mat('m_heat','SI',1) k rho C Hf]
%        k: conductivity
%        rho: mass density
%        C: ehat capacity
%        Hf: heat exchange coeffcient
%       Subtype 2 : anisotropic 
%          [MatId fe_mat('m_heat','SI',2) k1 k2 k3 rho C]


%       See sdtweb      fem (handling materials section), pl, fe_mat
%       See also help   p_heat


%       Etienne Balmes, Jean-Philippe Bianchi
%       Copyright (c) 2001-2013 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*NOSEM>

if nargin<1; help p_heat;return; end
if ischar(varargin{1}); [CAM,Cam]=comstr(varargin{1},1); pl=[]; carg=2;
else pl=varargin{1};[CAM,Cam]=comstr(varargin{2},1); carg=3;
end

%% #Default ------------------------------------------------------------------
if comstr(Cam,'default')
 
 model=[];if carg<=nargin;model=varargin{carg};carg=carg+1;end
 if isempty(model); out=m_heat('database');out=out(1);
 else;              out=fe_mat(['defaultil' CAM(8:end)],model);
 end
 
%% #Info ---------------------------------------------------------------------
elseif comstr(Cam,'info')

 r1=m_heat('database');fprintf('%s\n',r1.name);

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
  elseif ischar(st); mat=m_heat('database',st); 
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
  if ~isempty(pl) ; i2=find(pl(:,1)==r1(1)); else i2=[];end
  if isempty(i2)  ; i2=size(pl,1)+1;end
  pl(i2,1:length(r1))=r1;  %#ok<AGROW>
  out1(end+1,1:3)={'mat',sprintf('%i_%s',mat.pl(1),mat.name),mat};%#ok<AGROW> 
  if carg>nargin; break;end
  CAM=varargin{carg};carg=carg+1;if ischar(CAM);[CAM,Cam]=comstr(CAM,1);end
 end
 out=pl;
 
 
%% #Database -----------------------------------------------------------------
elseif comstr(Cam,'database')

  st=comstr(CAM,9);
  if isempty(st)&&carg<=nargin; st=varargin{carg}; carg=carg+1;end
  [MatId,i2,i3,i4]=sscanf(st,'%i',1); if i2~=1; MatId=1;end
  st=comstr(st,i4);

  % From Chaigne : steel k [11-46], rho 7800 C [460-625], alpha 14
  out.pl=[MatId fe_mat('m_heat','SI',1) 35 7850 476 15]; % inox steel
  % k rho C alpha 
  out.name='default heat';
  out.type='p_heat';
  out.unit='SI';

  % Material database :
  % Aluminium: Groult thesis, p78
  out(end+1).pl=[MatId fe_mat('type','m_heat','SI',1) ...
                 237 2700 900 15]; %k rho C alpha
  out(end).name='Aluminum'; out(end).type='m_heat';  out(end).unit='SI';
  
  % Steel (inox): Merlette Thesis , p127
  out(end+1).pl=[MatId fe_mat('type','m_heat','SI',1) ...
                 35 7850 476 15]; %k rho C alpha
  out(end).name='Steel'; out(end).type='m_heat';  out(end).unit='SI';
  
  
  
  out1='Heat properties';out2=carg;

  i1=strmatch(comstr(st,-27),comstr({out.name},-27));
  if ~isempty(i1); out=out(i1);
  elseif comstr(st,'shell') % piezo shell
   r1=comstr(st(6:end),-1,[1 1000 .01 0]);
   out=out(1);out.il(2+[1:length(r1)])=r1;
  end
  

%% #PropertyUnitType ---------------------------------------------------------
elseif comstr(Cam,'propertyunittype')

 if nargin==1;out=1:2; return; end % return subtypes ID
 i1=varargin{carg}; out1={};
 switch i1
 case 1
   st=...
  {'MatId'    0       'sdtweb(''p_heat'')';
   'Type'     0       '';
   'k'        16      'conductivity';
   'rho'      3       'mass density';
   'C'        17      'heat capacity' %J/kg-K in SI units or BTU/lbm-^\circR in British units. 
   'Hf'    16.004     'heat exchange coefficient'}; % W/m^2/K
 case 2
   st=...
  {'MatId'    0       'sdtweb(''p_heat'')';
   'Type'     0       '';
   'k1'       16      'conductivity d1';
   'k2'       16      'conductivity d2';
   'k3'       16      'conductivity d3';
   'rho'      3       'mass density';
   'C'        17      'heat capacity'
   };
 otherwise; 
  error('Not implemented subtype %i',i1);
 end
 
 if ~isempty(strfind(Cam,'cell')); out=st; else; out=[st{:,2}]; end
%% #SubType-------------------------------------------------------------------------
elseif comstr(Cam,'subtype');[CAM,Cam]=comstr(CAM,8);

st={'1 (Isotropic)'}; 
if carg<=nargin;
 i1=varargin{carg}; carg=carg+1;
 if ischar(i1);out=find(strncmpi(i1,st,1));
 else; 
  try; out=st{i1};
  catch; out=sprintf('m_heat %i',i1);
  end
 end
end

% --------------------------------------------------------------------
% -------------------------------------------------------------------------
%% #BuildConstit -------------------------------------------------------------------------
elseif comstr(Cam,'buildconstit')
 RunOpt=evalin('caller','RunOpt');
 if nargout<4;[out,out1,out2]=p_heat(varargin{:});
 else;[out,out1,out2,out3]=p_heat(varargin{:});
 end
elseif comstr(Cam,'tablecall');out='';
elseif comstr(Cam,'cvs');
 out='$Revision: 1.9 $  $Date: 2013/07/11 16:38:53 $'; return;
else; sdtw('''%s'' not known',CAM);
end
