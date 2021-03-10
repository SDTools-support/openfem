function [out,out1,out2]=p_spring(varargin)

%P_SPRING Spring element property function
%
%       Syntax : il = p_spring('default') 
%                il = p_spring('database ProId value') 
%                il = p_spring('dbval ProId value');
%                il = p_spring('dbval -unit MM ProId value');
%
% 	    Subtype 1 : point to point connection
%            [ProId type k m c eta S]
%            type  fe_mat('p_spring','SI',1)
%            stiffness, mass, damping, loss factor, stress coefficient
%       Subtype 2 : BUSH
%            [ProId Type k11 k22 k33 ... k66 c11 ... c66(12)
%               eta SA ST EA ET Mass Volume]
%            type  fe_mat('p_spring','SI',2)
%            ki : stiffness for each direction
%            ci : viscous damping for each direction
%            SA ST EA ET : stress/strain coef
%       Subtype 3 : PGAP
%            [ProId Type U0 F0 KA KB KT MU1 MU2 TMax MAR TRMin
%            type  fe_mat('p_spring','SI',3)
%            U0   : Offset
%            F0   : Pre-load
%            KA   : Axial stiff closed
%            KB   : Axial stiff opened
%            KT   : Transvers stiff
%            MU1  : Static friction
%            MU2  : Dynamic friction
%            TMax : Max penetration
%            MAR  : Max adjustment ratio
%       Subtype 4 : BUSH_FULLMAT
%            [ProId Type k11 k12 ... k21 k22 ...k66 c11 c12 ... c66 eta Mass]
%            type  fe_mat('p_spring','SI',4)
%            kij : stiffness matrix coefficients
%            cij : viscous damping matrix coefficients
%            eta : loss factor
%
%       Examples
%                p_spring('database 100 1e12')
%                p_spring('dbval 100 1e12')
%
%       See sdtweb     fem (handling materials section), pl, fe_mat, p_shell
%       See also help  fe_mat


%	Jean-Michel Leclere
%       Copyright (c) 2001-2020 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC>
if nargin<1; help p_spring;return; end
if ischar(varargin{1}); [CAM,Cam]=comstr(varargin{1},1); il=[]; carg=2;
else; il=varargin{1};[CAM,Cam]=comstr(varargin{2},1); carg=3;
end

%% #PropertyUnitType ---------------------------------------------------------
if comstr(Cam,'propertyunittype')

 if nargin==1;out=1:2; return; end % return subtypes ID
 i1=varargin{carg};
 out1={};
 switch i1 % PropertySubType
 case 1 % [matid type K M C eta S %the unit of S is unclear ]
  st={...
      'ProId',              0     'sdtweb(''p_spring'')';
      'Type',               0     '';
      'k'                   13    'Stiffness';
      'm'                   9     'Mass';
      'c'                   2.005 'Viscous damping';
      'Eta'                 0     'Loss factor';
      'S'                   0     'Stress coef'};
 %      'Hysteretic damping'  2.005 'Hysteretic damping';
 case 2 % PBUSH (Based on Nastran)
 % ProId Type k1(3) k2 k3 k4 k5 k6 b1(9) b2 b3 b4 b5 b6 eta 
 %  SA(16) ST EA ET m v
  st={...
      'ProId',              0        'sdtweb(''p_spring'')';
      'Type',               0        '';
      'k1'                  13       'Stiffness';
      'k2'                  13       'Stiffness';
      'k3'                  13       'Stiffness';
      'k4'                  13       'Stiffness';
      'k5'                  13       'Stiffness';
      'k6'                  13       'Stiffness';
      'c1'                  2.005    'Viscous damping';
      'c2'                  2.005    'Viscous damping';
      'c3'                  2.005    'Viscous damping';
      'c4'                  2.005    'Viscous damping';
      'c5'                  2.005    'Viscous damping';
      'c6'                  2.005    'Viscous damping';
      'Eta'                 0        'Loss factor';
      'SA'                  0        'Stress recovery coef for translations';
      'ST'                  0        'Stress recovery coef for rotations';
      'EA'                  0        'Strain recovery coef for translations';
      'ET'                  0        'Strain recovery coef for rotations';
      'm'                   9        'Mass';
      'v'                   11.004   'Volume'};
 case 3 % PGAP (Based on Nastran)
  st={...
      'ProId',              0        'sdtweb(''p_spring'')';
      'Type',               0        '';
      'U0'                  4       'Offset';
      'F0'                  2       'Pre-load';
      'KA'                  13       'Axial stiff closed';
      'KB'                  13       'Axial stiff opened';
      'KT'                  13       'Transvers stiff';
      'MU1'                 0        'Static friction';
      'MU2'                 0        'Dynamic friction';
      'TMax'                4    'Max penetration';
      'MAR'                 0    'Max adjustment ratio';
      'TRMin'               0    'TRMin';
      };
  case 4 % BUSHING_FULLMAT
   st={...
      'ProId',              0        'sdtweb(''p_spring'')';
      'Type',               0        '';
      'k11'                 13       'Stiffness';
      'k12'                 13       'Stiffness';
      'k13'                 13       'Stiffness';
      'k14'                 13       'Stiffness';
      'k15'                 13       'Stiffness';
      'k16'                 13       'Stiffness';
      'k21'                 13       'Stiffness';
      'k22'                 13       'Stiffness';
      'k23'                 13       'Stiffness';
      'k24'                 13       'Stiffness';
      'k25'                 13       'Stiffness';
      'k26'                 13       'Stiffness';
      'k31'                 13       'Stiffness';
      'k32'                 13       'Stiffness';
      'k33'                 13       'Stiffness';
      'k34'                 13       'Stiffness';
      'k35'                 13       'Stiffness';
      'k36'                 13       'Stiffness';
      'k41'                 13       'Stiffness';
      'k42'                 13       'Stiffness';
      'k43'                 13       'Stiffness';
      'k44'                 13       'Stiffness';
      'k45'                 13       'Stiffness';
      'k46'                 13       'Stiffness';
      'k51'                 13       'Stiffness';
      'k52'                 13       'Stiffness';
      'k53'                 13       'Stiffness';
      'k54'                 13       'Stiffness';
      'k55'                 13       'Stiffness';
      'k56'                 13       'Stiffness';
      'k61'                 13       'Stiffness';
      'k62'                 13       'Stiffness';
      'k63'                 13       'Stiffness';
      'k64'                 13       'Stiffness';
      'k65'                 13       'Stiffness';
      'k66'                 13       'Stiffness';
      'c11'                 2.005    'Viscous damping';
      'c12'                 2.005    'Viscous damping';
      'c13'                 2.005    'Viscous damping';
      'c14'                 2.005    'Viscous damping';
      'c15'                 2.005    'Viscous damping';
      'c16'                 2.005    'Viscous damping';
      'c21'                 2.005    'Viscous damping';
      'c22'                 2.005    'Viscous damping';
      'c23'                 2.005    'Viscous damping';
      'c24'                 2.005    'Viscous damping';
      'c25'                 2.005    'Viscous damping';
      'c26'                 2.005    'Viscous damping';
      'c31'                 2.005    'Viscous damping';
      'c32'                 2.005    'Viscous damping';
      'c33'                 2.005    'Viscous damping';
      'c34'                 2.005    'Viscous damping';
      'c35'                 2.005    'Viscous damping';
      'c36'                 2.005    'Viscous damping';
      'c41'                 2.005    'Viscous damping';
      'c42'                 2.005    'Viscous damping';
      'c43'                 2.005    'Viscous damping';
      'c44'                 2.005    'Viscous damping';
      'c45'                 2.005    'Viscous damping';
      'c46'                 2.005    'Viscous damping';
      'c51'                 2.005    'Viscous damping';
      'c52'                 2.005    'Viscous damping';
      'c53'                 2.005    'Viscous damping';
      'c54'                 2.005    'Viscous damping';
      'c55'                 2.005    'Viscous damping';
      'c56'                 2.005    'Viscous damping';
      'c61'                 2.005    'Viscous damping';
      'c62'                 2.005    'Viscous damping';
      'c63'                 2.005    'Viscous damping';
      'c64'                 2.005    'Viscous damping';
      'c65'                 2.005    'Viscous damping';
      'c66'                 2.005    'Viscous damping';
      'Eta'                 0        'Loss factor';
      'm'                   9        'Mass';
      };
     
 otherwise; st={'ProId' 0 'sdtweb(''p_spring'')'; 'Type', 0, ''};
 end
 if ~isempty(strfind(Cam,'cell')); out=st; else; out=[st{:,2}]; end

%% #Default ------------------------------------------------------------------
elseif comstr(Cam,'default')

  out=p_spring('database');
  out=out(1);

%% #DBval --------------------------------------------------------------------
elseif comstr(Cam,'dbval')

 while 1==1
  i1=strfind(comstr(Cam,-27),'-unit'); out1={};
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
  %if ischar(CAM); [i1,CAM,Cam]=comstr(CAM,'dbval','%i');else; i1=[];end
  %if isempty(CAM)&&carg<=nargin; st=varargin{carg};carg=carg+1;
  %else st=CAM;end
  if ischar(CAM); [CAM,Cam,i1]=comstr('dbval',[-25 1 1],CAM,Cam);
   if ~isempty(Cam);CAM=sprintf('%i %s',i1(1),CAM);end
  else; i1=[];
  end
  if isempty(CAM)&&carg<=nargin; st=varargin{carg};carg=carg+1;
  else; st=CAM;end
  if isempty(st)
  elseif ischar(st);
   [mat,st1,i2]=p_spring('database',st,varargin{carg:end});carg=carg+i2-3;
  elseif isnumeric(st)
   [typ,st1,i4]=fe_mat('typem',st(2));
   mat=struct('il',st,'name',sprintf('%i',st(1)),'type',typ,'unit',st1);
  end
  if ~isempty(PUnit)
   r1=fe_mat(sprintf('convert %s %s',mat.unit,PUnit),mat.il(1:2));
   mat.il(2)=r1(2); mat.unit=PUnit;
  end
  if ~isempty(Unit)
   mat.il=fe_mat(sprintf('convert %s %s',mat.unit,Unit),mat.il);mat.unit=Unit;
  end
  r1=mat.il; if length(i1)==1; r1(1)=i1;end
  if ~isempty(il); i2=find(il(:,1)==r1(1)); else i2=[];end
  if isempty(i2); i2=size(il,1)+1;end
  il(i2,1:length(r1))=r1;
  if carg>nargin; break;end
  CAM=varargin{carg};carg=carg+1;if ischar(CAM);[CAM,Cam]=comstr(CAM,1);end
 end
 out=il;

%% #database -----------------------------------------------------------------
elseif comstr(Cam,'database') 

  st=comstr(CAM,9);
  if isempty(st)&&carg<=nargin; st=varargin{carg}; carg=carg+1;end
  [MatId,i2,i3,i4]=sscanf(st,'%i',1); if i2~=1; MatId=1;end
  st=comstr(st,i4);

  out.il=[MatId fe_mat('p_spring',1,1) 1e15]; 
  out.name='stiff spring 1e15';
  out.type='p_spring';
  out.unit='SI';

  out(2).il=[MatId fe_mat('p_spring',1,1) 1e6]; 
  out(2).name='stiff spring 1e6';
  out(2).type='p_spring';
  out(2).unit='SI';

  i1=strmatch(st,{out.name},'exact');

  if ~isempty(i1); out=out(i1);
  elseif ~isempty(st)
   r1=str2num(st); 
   if length(r1)<1; error('Not a consistent database call');end
   il=[MatId fe_mat('p_spring',1,1) r1(:)']; 
   name=['Spring' st];
   type='p_spring';
   unit='SI';
   out=struct('name',name,'il',il,'type',type,'unit',unit);
  end

  out1='Spring/rigid connection';out2=carg;

%% #SubTYpeString -------------------------------------------------------------------------
elseif comstr(Cam,'subtypestring')

 i1=varargin{carg}; carg=carg+1;
 switch i1
 case 1;  out='p_spring';
 otherwise; out='p_spring';
 end

%% #TableCall -------------------------------------------------------------------------
elseif comstr(Cam,'tablecall');out='';
elseif comstr(Cam,'test');out='';
elseif comstr(Cam,'cvs');
 out='$Revision: 1.38 $  $Date: 2020/08/26 15:12:11 $'; return;
else; sdtw('''%s'' not known',CAM);
end
