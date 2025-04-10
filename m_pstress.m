function [out,out1,out2]=m_pstress(varargin)

% M_PSTRESS pre stress function
%
%       Syntax : mat= m_pstress('default') 
%                mat= m_pstress('database name') 
%                pl = m_pstress('dbval MatId name');
%                pl = m_pstress('dbval -unit MM MatId name');
%                pl = m_pstress('dbval -punit MM MatId name');
%
%
%       Material subtypes supported by m_pstress are
%
%	Subtype 1 : linear case
%	    [MatId type thetai sigma0 t]
%         with
%           type  fe_mat('m_pstress','SI',1)
%           thetai sigma0 (pre stress coefficients)
%           t (time)
%
%	Subtype 2 : non linear case
%	    [MatId type thetai sigma0 t]
%         with
%           type  fe_mat('m_pstress','SI',2)
%           thetai sigma0 (pre stress coefficients)
%           t (time)
%
%       See <a href="matlab: sdtweb _taglist m_pstress">TagList</a>
%       See sdtweb     m_elastic, pl, fem
%       See also help  fe_mat, p_shell, p_beam


%	Etienne Balmes, Jean-Michel Leclere, Corine Florens
%       Copyright (c) 2001-2025 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

if nargin<1; help m_pstress;return; end
if ischar(varargin{1}); [CAM,Cam]=comstr(varargin{1},1); pl=[]; carg=2;
else pl=varargin{1};[CAM,Cam]=comstr(varargin{2},1); carg=3;
end
if sp_util('diag');fprintf('m_pstress %s \n',CAM);end


%% #PropertyUnitType ---------------------------------------------------------
if comstr(Cam,'propertyunittype')

 if nargin==1; out=1;return;end
 i1=varargin{carg};
 out1={};
 switch i1
 case 1
 st=...
 {'MatId'   0  'help m_pstress';
  'Type'    0  '';
  'thetai'  0  'Pre stress coefficient';
  'sigma0'  0  'Pre stress coefficient';
  't'       0  'Time';
 };
 otherwise; st={'MatId' 0 'help m_pstress'; 'Type', 0, ''};
 end
 if ~isempty(strfind(Cam,'cell')); out=st; else; out=[st{:,2}]; end


%% #Info ---------------------------------------------------------------------
elseif comstr(Cam,'info')

 r1=m_pstress('database');
 disp(char({r1.name}))

%% #Default ------------------------------------------------------------------
elseif comstr(Cam,'default')

  out=m_pstress('database'); out=out(1);

%% #DBval --------------------------------------------------------------------
elseif comstr(Cam,'dbval') 

 while 1==1
  i1=strfind(comstr(Cam,-27),'-unit'); out1={}; % See if unit specified
  if ~isempty(i1)
   [Unit,i2,i3,i4]=sscanf(CAM(i1+5:end),'%s',1);
   i4=i1+[0:4+i4];CAM(i4)=''; [CAM,Cam]=comstr(CAM,1);
  else Unit=''; end
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
  elseif ischar(st); mat=m_pstress('database',st);
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
  if isempty(i2);  i2=size(pl,1)+1;end
  pl(i2,1:length(r1))=r1;
  out1(end+1,1:3)={'mat',sprintf('%i_%s',mat.pl(1),mat.name),mat};
  if carg>nargin; break;end
  CAM=varargin{carg};carg=carg+1;if ischar(CAM);[CAM,Cam]=comstr(CAM,1);end
 end
 out=pl;

%% #d-------------------------------------------------------------------------
elseif comstr(Cam,'database'); [st,Cam]=comstr(CAM,9);

  if isempty(st)&&carg<=nargin; st=varargin{carg}; carg=carg+1;end
  [MatId,i2,i3,i4]=sscanf(st,'%i',1); if i2~=1; MatId=1;end
  st=comstr(st,i4);

  out.pl=[MatId fe_mat('type','m_pstress','SI',1) 0.1 0.5 0.0];
  out.name='Ref';
  out.type='m_pstress';
  out.unit='SI';

%   out(2).pl=[MatId fe_mat('type','m_pstress','SI',1) 0.1 0.5 0.0];
%   out(2).name='Rivlin';
%   out(2).type='m_pstress';
%   out(2).unit='SI';
%   
%   out(3).pl=[MatId fe_mat('type','m_pstress','SI',1) 0.1 0.5 0.0];
%   out(3).name='phi';
%   out(3).type='m_pstress';
%   out(3).unit='SI';
% 
% 
%   out(4).pl=[MatId fe_mat('type','m_pstress','SI',1) 0.1 0.5 0.0];
%   out(4).name='coeur';
%   out(4).type='m_pstress';
%   out(4).unit='SI';

  i1=strmatch(comstr(st,-27),comstr({out.name},-27),'exact');
  out1='pre_stress';

  if isempty(i1) && isempty(st); return; end

  % match a name 
  if ~isempty(i1); out=out(i1(1));

  else % assume values given
    error('Not a supported material');
  end


%% #BuildConstit -------------------------------------------------------------
% Implementation of elastic constitutive law building 3D
elseif comstr(Cam,'buildconstit')

  ID=varargin{carg};carg=carg+1;
  pl=varargin{carg};carg=carg+1; mat=pl(find(pl(:,1)==ID(1)),:);
  il=varargin{carg};carg=carg+1; pro=[];try;pro=il(find(il(:,1)==ID(2)),:);end
  if isempty(pl); error('MatId is not matched');end

   [st,unit,typ]=fe_mat('type',mat(2));
   
   if strcmp(st,'m_pstress')&&typ==1 % ...

    % set viscous coef to zero for static case
    out=[mat(3:5) 0];out=out(:); 
    
   elseif strcmp(st,'m_pstress')&&typ==2 
   
   %propagate constits, viscosity enabled for dynamic behavior
    out=[mat(3:5) 1]; out=out(:);
   else
    error(sprintf('Volume law not implemented for this material %s(%i)',st,typ));
   end 
   % Integration rule saved in integ
   if length(ID)<4; ID(4)=0;end; ID(3)=3*ID(4); % must be number of DOFs 
   if size(pro,1)&&size(pro,2)>3;  ID(5)=pro(1,4); end  
   ID(7)=105; % used for later NL constit type selection
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
 otherwise; out='m_pstress';
 end
%% #T-------------------------------------------------------------------------
elseif comstr(Cam,'tablecall'); out='';
elseif comstr(Cam,'cvs')
 out='$Revision: 1.11 $  $Date: 2025/04/07 17:08:09 $';
else; sdtw('''%s'' not known',CAM);
end
 
end
% -------------------------------------------------------------------------

