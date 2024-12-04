function [out,out1,out2,out3]=p_shell(varargin)

%P_SHELL Shell element property function
%
%       Syntax : pro= p_shell('default') 
%                pro= p_shell('database ProId name') 
%                il = p_shell('dbval ProId name');
%                il = p_shell('dbval -unit MM ProId name');
%                il = p_shell('dbval -punit MM ProId name');
%                il = p_shell('SetDrill d');
%
%       Supported subtypes are
%       1 : standard isotropic
%	    [Id# Type f   d 0 h   k MID2 12I/T^3 MID3 NSM Z1 Z2 MID4]
%          with
%	     Type : fe_mat('p_shell','SI',1)
%            f    : formulation (see the element : quad4, ...)
%            d    : -1 no drilling stiffness, If d==0 d is set to 1.
%		    d>0 drilling DOF stiffness coefficient
%            h    : thickness
%            k    : shear factor (default value is 5/6)
%            12I/T^3 : Ratio of bending moment of inertia to nominal T^3/12
%                      (default 1)
%            NSM  : non structural mass per unit area
%            MID2, MID3, MID4 : materials for classical lamination thery
%            Z1, Z2 currently unused
%
%       2 : composite (follows NASTRAN's PCOMP format)
%           [ProId type   Z0 NSM SB d TREF GE f MatId1 T1 Theta1 SOUT1 ...]
%           Type fe_mat('p_shell','SI',2)
%           z0   distance from the reference plane to the bottom surface 
%                default -(total thickness)/2
%           SB   allowable shear stress (unused)
%           d    drilling coefficient (see above)
%           TREF reference temperature (unused)
%           GE   loss factor
%           f    formulation
%           MatId1,T1,Theta1,SOUT1 : material idenfier for the layer,
%                layer thickness, orientation, layer options (unused)
%
%       For automated generation of properties see, the P_SHELL
%       dbval and database commands.
%
%       See sdtweb      p_shell, fem (handling materials section), pl, fe_mat
%       See also help   fe_mat

%       Etienne Balmes, Corine Florens
%       Copyright (c) 2001-2020 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

if nargin<1; help p_solid;return; end
if ischar(varargin{1}); [CAM,Cam]=comstr(varargin{1},1); il=[]; carg=2;
else;il=varargin{1};[CAM,Cam]=comstr(varargin{2},1); carg=3;
end
%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*NOSEM>

%% #PropertyUnitType ----------------------------------------------------------
if comstr(Cam,'propertyunittype')

 if nargin==1;out=[1 2]; return;end % return subtypes ID
 i1=varargin{carg}; out1={};
 switch i1 % PropertySubType
 case 1 % standard fe_mat('unitlabel') 
     st={...
         'ProId'             0  'sdtweb(''p_shell'')';
         'Type'              0  '';
         'f'                 0  'Formulation';
         'd'                 0  'Drilling DOF'; 
         'O'                 0  'Unused';
         'h'                 4  'Thickness';
         'k'                 0  'Shear factor';
         'MID2'              0  'Bending material';
         'RatI12_T3'         0  'Bending Ratio T3/I12';
         'MID3'              0  'Transverse shear material';
         'NSM'               9  'Non structural mass per unit area';
         'Z1'                4  'Unused. Offset';
         'Z2'                4  'Unused. Offset';
         'MID4'              0  'Membrane/bending material';
          };
%   r2.value=[0 1 2 3 4 5];
%   r2.label={'Unknown','Kirchhoff','Mindlin','Unused', ...
%     'OpenFEM MITC4','SDTools MITC4'};

 case 2 % composite
     st={...
         'ProId'             0  'sdtweb(''p_shell'')';
         'Type'              0  '';
         'Z0'                4  'Position of bottom in mesh coordinates'
         'NSM'               9  'Non structural mass per unit area';
         'SB'                0  'Allowable shear stress of the bonding material';
         'FT'                0  'Failure Theory';
         'TREF'              7  'Reference temperature';
         'Eta'               0  'Hysteretic loss factor';
         'LAM'               0  'Laminate type'
         };
     out1={         'MatId%i'            0 'Material';
         'T%i'                4 'Thickness';
         'Theta%i'            0 'Orientation (deg)';
         'SOUT%i'             0 'Stress flag'};
 case 3 % General CLPT for split
     st={...
         'ProId'             0  'sdtweb(''p_shell'')';
         'Type'              0  '';
         'rhoh',0,'';'eta',0,'';'f',0,'';'d',0,''; 'h',0,''; 
         'A11' 0 '';'A12' 0 '';'A13' 0 '';
         'A22' 0 '';'A23' 0 '';'A33' 0 '';
         'B11' 0 '';'B12' 0 '';'B13' 0 '';
         'B21' 0 '';'B22' 0 '';'B23' 0 '';
         'B31' 0 '';'B32' 0 '';'B33' 0 '';
         'D11' 0 '';'D12' 0 '';'D13' 0 '';
         'D22' 0 '';'D23' 0 '';'D33' 0 '';
         'F11' 0 '';'F12' 0 '';'F22' 0 '';'Dr',0,''
         };
  % st=p_shell('propertyunittype -cell',3);st=st(:,1);st(1:2)='';
  % dd=cell(9,9); dd([1 10 19 11 20 21 27+[1 10 19 2 11 20 3 12 21] 30+[1 10 19 11 20 21] 61 70 71 81])=st
  % topo(i3);
 otherwise; st={'ProId' 0; 'Type', 0};
 end
 if ~isempty(strfind(Cam,'cell')); out=st; else; out=[st{:,2}]; end

%% #StressCrit ---------------------------------------------------------------
elseif comstr(Cam,'stresscrit');

    r1=varargin{carg};carg=carg+1;
    opt=varargin{carg};carg=carg+1;
    if length(opt)<2; opt(2)=0; end 
    TensorTopology=[1 3;3 2];
    r1=r1(1:3,:)+opt(1)*r1(4:6,:);  % (e11,e22,e12)+h*(k11,k22,k12)
    out=zeros(1,size(r1,2)); out1=zeros(2,size(r1,2));
    for j1=1:size(r1,2)
        r2=r1(:,j1); [r3,r2]=eig(r2(TensorTopology));
        switch opt(2) 
            case 0; [out(j1),i2]=max(diag(r2)); % max
            case -1; [out(j1),i2]=min(diag(r2)); % min
            otherwise; r2=diag(r2); i2=opt(2);out(j1)=r2(i2);
        end
        out1(:,j1)=r3(:,i2);
    end
    
%% #SetDrill -----------------------------------------------------------------
elseif comstr(Cam,'setdrill');

  opt=comstr(CAM(9:end),[-1 0]);
  if isempty(opt); opt=0; end
  il=varargin{carg}; carg=carg+1;
  for j1=1:size(il,1)
    [m_function,UnitCode,SubType]=fe_mat('typem',il(j1,2));
    if ~comstr(m_function,'m_shell')
    elseif SubType==1; il(j1,4)=opt; 
    end
  end
  out=il;

%% Dbval ---------------------------------------------------------------------
elseif comstr(Cam,'dbval')

 while 1==1
  i1=strfind(comstr(Cam,-27),'-unit'); out1={};
  if ~isempty(i1)
   [Unit,i2,i3,i4]=sscanf(CAM(i1+5:end),'%s',1);
   i4=i1+[0:4+i4];CAM(i4)=''; [CAM,Cam]=comstr(CAM,1);
  else;Unit='';
  end
  i2=strfind(comstr(Cam,-27),'-punit');
  if ~isempty(i2)
   [PUnit,i3,i4,i5]=sscanf(CAM(i2+6:end),'%s',1);
   i5=i2+[0:5+i5-1];CAM(i5)=''; [CAM,Cam]=comstr(CAM,1);
  else;PUnit='';
  end
  [CAM,Cam]=comstr(CAM,6);
  [i1,CAM,Cam]=comstr(CAM,'','%i');
  if isempty(CAM)&&carg<=nargin; st=varargin{carg};carg=carg+1;
  else;st=CAM;end
  [mat,st1,i2]=p_shell('database',st,varargin{carg:end});carg=carg+i2-3;
  if ~isempty(PUnit)
   r1=fe_mat(sprintf('convert %s %s',mat.unit,PUnit),mat.il(1:2));
   mat.il(2)=r1(2); mat.unit=PUnit;
  end
  if ~isempty(Unit)
   mat.il=fe_mat(sprintf('convert %s %s',mat.unit,Unit),mat.il);mat.unit=Unit;
  end
  r1=mat.il; if length(i1)==1; r1(1)=i1;end
  if ~isempty(il); i2=find(il(:,1)==r1(1)); else;i2=[];end
  if isempty(i2); i2=size(il,1)+1;end
  il(i2,1:length(r1))=r1; %#ok<AGROW>
  if carg>nargin; break;end
  CAM=varargin{carg};carg=carg+1;if ischar(CAM);[CAM,Cam]=comstr(CAM,1);end
 end
 out=il;

%% #DataBase -----------------------------------------------------------------
elseif comstr(Cam,'database')


  st=comstr(CAM,9);
  if isempty(st)&&carg<=nargin; st=varargin{carg}; carg=carg+1;end
  [MatId,i2,i3,i4]=sscanf(st,'%i',1); if i2~=1; MatId=1;end
  st=comstr(st,i4);

  out.il=[MatId fe_mat('p_shell',1,1) 1 1 0 .01]; 
  out.name='Kirch';
  out.type='p_shell';
  out.unit='SI';

  out(2).il=[MatId fe_mat('p_shell',1,1) 0 1 0 .01];
  out(2).name='Mind';
  out(2).type='p_shell';
  out(2).unit='SI';

  i1=find(strcmpi(st,{out.name}));
  if ~isempty(i1); out=out(i1);
  elseif ~isempty(st)
   [out,carg]=create_section_shell(st,out,varargin,carg);  out.il(1)=MatId;
  end
  out1='Shell properties';out2=carg;
%% #Default ------------------------------------------------------------------
elseif comstr(Cam,'default'); out=p_shell('database');out=out(1);

%% #Const : definition of EltConst, for SDT see sdtweb q4cs('constants') -----
elseif comstr(Cam,'const')

 if carg>nargin; EC=integrules('quadb');
 else; EC=varargin{carg};carg=carg+1;
 end
 if ~ischar(EC);
 elseif comstr(EC,'tria3'); EC=integrules(EC,-1); % tria3 at center
 else; EC=integrules(EC);
 end
 if carg>nargin; integ=[];else; integ=varargin{carg};carg=carg+1;end
 if carg<=nargin; constit=varargin{carg};carg=carg+1;else;constit=[];end
 if ~isempty(constit)&&constit(1)==-1; 
  st=fe_mat('typep',constit(2,1));[out,out1]=feval(st,'const',varargin{2:end});
  return;
 end
 if size(integ,1)>3&&~any(integ(3:4,1))% back propagate number of DOFs if needed
     i1=6; if size(integ,1)<4; error('Improper init');end% dof per elt
     sp_util('setinput',integ, ...
         [EC.Nnode*i1;EC.Nnode],int32(3:4),int32(1));
     EC.ElMap=elem0('elmapmat_og',[EC.Nnode i1]);
 end
   
 EC.Nw=size(EC.N,1);
 if ~isfield(EC,'Nnode'); EC.Nnode=size(EC.N,2);end
 if size(integ,1)<6;rule=[1 EC.Nw];
 else; rule=integ(5:6);rule=rule(:)';
 end
 EC.DofLabels={'u','v','w','\theta_u','\theta_v','\theta_w'};
 % Strain energy
 % define the deformation vector: row, NDN, DDL, NwStart, NwRule
   % define the deformation vector: row, NDN, DDL
   if size(EC.w,1)==13; rule1=[5 9];rule=[1 4];else; rule1=rule;end
   r2=[1 2 1 rule1; 2 3 2 rule1;3 3 1 rule1;3 2 2 rule1;% e_xx =N,x u,e_yy,e_xy
     4 2 5 rule1;5 -3 4 rule1;6 -2 4 rule1;6 3 5 rule1; % curvatures
     7 2 3 rule;7 -1 5 rule;8 3 3 rule;8 1 4 rule];  % shear constraint
   if size(EC.w,1)==13;
    r2=[r2; ...
    9 1 6 rule;9 1 6 rule;9 -2 2 rule;9 3 1 rule]; % dril 2\theta_z-(v,x-u,y)
   end
   % Topology
   % [rho eta Gii] in constit -> shift values of col 5 by two
    dd=[ones(6) zeros(6,2);zeros(2,6) ones(2)];dd(9,9)=1;
    dd(dd~=0)=[1:length(find(dd))]+8; % see integinfo call
 EC.StrainDefinition{1}=r2;
 EC.StrainLabels{1}={'e_x','e_y','e_xy','k_x','k_y','k_xy','s_x','s_y','d_z'};
 EC.ConstitTopology{1}=int32(dd);
 EC.TensorTopology=int32([1 3 0 0 0 0;3 2 0 0 0 0;0 0 4 6 0 0;0 0 6 5 0 0; ...
     0 0 0 0 7 0;0 0 0 0 0 8]);
 % Kinetic energy
 EC.StrainDefinition{2}= [1 1 1 rule;2 1 2 rule; 3 1 3 rule];
 EC.StrainLabels{2}={'u','v','w'};
 EC.ConstitTopology{2}=int32(eye(3));
  if size(constit,1)==50  % Attempt to have mass on rotations
   % constit(3)=constit(3)-remi(constit(3),[],2)*10;% need correction
   EC.StrainDefinition{2}= [1 1 1 rule;2 1 2 rule; 3 1 3 rule; ...
       ;4 1 4 rule;; 5 1 5 rule];    
   EC.StrainLabels{2}={'u','v','w','rx','ry'};
   EC.ConstitTopology{2}=int32(diag([1 1 1 50 50]));   
  end
 EC=integrules('matrixrule',EC);
 EC.VectMap=int32(reshape(1:6*EC.Nnode,6,EC.Nnode)'); 
 EC.pEC=zeros(40,1,'int32'); % copy of EC C structure for callbacks
 %display with :integrules('texstrain',EC);

 if ~isempty(constit)
  % Only keep terms needed for the current law.
  i1=find(any(constit,2))-1; 
  for j0=1:length(EC.MatrixIntegrationRule);
    r1=EC.MatrixIntegrationRule{j0};
    if ~isempty(r1)
     EC.MatrixIntegrationRule{j0}=r1(ismember(double(r1(:,5)),i1),:);
    end
  end
 end
 if nargout==0
  integrules('texstrain',EC);
  try; EC=integrules('stressrule',EC);integrules('texstress',EC);end
 else; out=EC;
 end
 out1=23; % 3d surface
 if nargout==3 % Build normal map - - - - - - - - - - - - - - - - - - - -
    model=varargin{carg};carg=carg+1;
    Case=varargin{carg};carg=carg+1;
    cEGI=varargin{carg};carg=carg+1;
    RunOpt=varargin{carg};carg=carg+1;
    elt=model.Elt([cEGI(1)-1;cEGI(:)],:);

    % obsolete Shell_%i
    out2=stack_get(model,'MAP',sprintf('Shell_%i',Case.jGroup),'getdata');
    if ~isempty(out2); error('Obsolete use MAP in Pro entry');end
    % 'MAP','Group' entries are obsolete and should diss
    r1=stack_get(model,'MAP',sprintf('Group_%i',Case.jGroup),'getdata');
    if ~isempty(r1); 
     sdtw('''MAP'',''GroupI'' entries are obsolete');
    else;r1=Case.GroupInfo{Case.jGroup,7};
    end
    [ElemF,i1,ElemP]=feutil('getelemf',elt(1,:),Case.jGroup);
    try; % improved map at nodes
     out2=feutilb('shellmapnodepos',Case.Node,elt,r1);  
    catch; % if not available use element normals at center
     sdtw('ShellMap failed : no composite support');
     MAP=feutil('getnormal map',Case.Node,elt); 
     out2=[MAP.normal MAP.normal MAP.normal MAP.normal]';
     % Currently field at DOFs that is replicated everywhere
     % will need to become a continuous normal when possible
     out2=struct('data',reshape(out2,3,numel(out2)/3), ...
      'NodePos',reshape(int32(1:4*length(cEGI)),4,length(cEGI)), ...
      'lab',{{'v3x','v3y','v3z'}});
    end    
 end

%% #BuildConstit [constit,integ]=p_shell('buildconstit',ID,pl,il);
elseif comstr(Cam,'buildconstit')

ID=varargin{carg};carg=carg+1;out1=int32(ID);
pl=varargin{carg};carg=carg+1;
il=varargin{carg};carg=carg+1;il0=il;
il=il(il(:,1)==ID(2),:);

if size(il,1)>1; error('BuildConstit assumes one p_shell');end

dd=zeros(6); ds=zeros(2); % dd: [A B;B D]  DS shear 
rhoh=0; etah=0; out1=int32(ID(:));out2=[];out3=struct;
RunOpt=struct;% Needed to store .Dim

if isempty(il); out=sprintf('Error ProId %i not found',ID(2));return;
elseif il(2)==2; st='p_shell';typ=1; unit='US'; % OBSOLETE
else; [st,unit,typ]=fe_mat('typep',il(2));
end
if strcmp(st,'p_shell')
if length(ID)<4;ID(3:4)=[6 1];else;ID(3)=6*ID(4); end

switch typ
case 1
%%  #BuildConstit.Pshell_1 OneLayer plate - - - - - - - - - - - - - - - - - -

 if isempty(pl); i1=[]; else; i1=find(pl(:,1)==ID(1)); end
 if isempty(i1); 
  if length(il)>10&&il(11)&&~il(6) % ok to integrate NSM only
   mat=zeros(1,6); mat(2)=fe_mat('m_elastic','US',1);
  else;out=sprintf('Error MatId %i not found',ID(1));return;
  end
 else; mat=pl(i1,:);mat=mat(1,:);
 end
 %mat=pl(i1,:);mat=mat(1,:);
 [st,unit,typ]=fe_mat('typem',mat(2));
 % Shear correction factor
 if length(il)<7; il(7)=5/6; elseif il(7)==0; il(7)=5/6;end
 il(15)=0;  
 
 if ~strcmp(st,'m_elastic');
   error('Shell only implemented for m_elastic materials');
 elseif typ==1 % Isotropic material
   if length(mat)<6||mat(6)==0; mat(6)=mat(3)/2/(1+mat(4));end % G
   if size(mat,2)<7;mat(7)=0;end
   o1 = mat; o2 = il;
   o3 = mat(3)/(1-mat(4)^2); % membrane
   % was modified for SDT 6, G is used directly for membrane
   %o3 = [o3 mat(4)*o3 0;      mat(4)*o3 o3 0;      0 0 (1-mat(4))/2*o3]*il(6);
   o3 = [o3 mat(4)*o3 0;      mat(4)*o3 o3 0;      0 0 mat(6)]*il(6);
   o4 = o3*il(6)^2/12;     % bending
   ds = eye(2)*mat(6)*il(7)*il(6); % shear
   out3.rhoh=mat(5)*il(6);out3.eta=mat(7);
 elseif typ==6 % See CLPT section in SDT manual
   [r1,r1,r1,dd]=p_solid('buildconstit 3 1',[mat(1);1;3;1],pl,...
        p_solid('dbval 1 d3 -3'));
   dd=dd.dd;
   o3=zeros(3); % membrane
   for j1=1:2;for j2=1:2;o3(j1,j2)=dd(j1,j2)-dd(j1,3)*dd(3,j2)/dd(3,3);end;end
   o3(3,3)=dd(6,6); o3=o3*il(6); % xy
   o4=o3*il(6)^2/12;             % bending
   ds = dd(4:5,4:5)*il(7)*il(6); % shear
   out3.eta=0;if length(pl)>16; out3.eta=pl(17);end
   out3.rhoh=mat(12)*il(6);
 elseif typ==5
  % Orthotropic material for shell (E1, E2, nu12, G12,G1z,G2z)
  % n12 E2 = n21 E1, 
  % m_elastic('propertyunittype-cell',5);find(strcmpi(ans(:,1),'eta'))
  % E=[1/E1 -N12/E1 0;-N12/E1 1/E2 0;0 0 1/G12] S, shear=[G1z 0;0 G2z]
  % See sdtweb m_elastic('buildply')
  %S=[1/E -Nu/E 0 0 0;-Nu/E1 1/E 0 0 0;
     %    0 0  1/G 0 0;0 0 0 1/G 0;0 0 0 0 1/G];
     % Si=E*[1 nu;nu 1]/(1-nu^2)
     dd=zeros(3);
     dd([1 2 4 5 9])=[1/mat(3) -mat(5)/mat(3) -mat(5)/mat(3) 1/mat(4) 1/mat(6)];
     o3=inv(dd)*il(6); 
     o4=o3*il(6)^2/12;             % bending
     if mat(7)==0; 
       mat(7)=mat(6);
       if mat(6);warning('No Shear using, non validated default %g, please report',mat(7));end
     end
     if mat(8)==0; 
       mat(8)=mat(7);
     end
     ds = [mat(7) 0;0 mat(8)]*il(7)*il(6); % shear
     out3.rhoh=mat(9)*il(6); 
     out3.eta=0;if length(pl)>16; out3.eta=pl(17);end
 else
  error(sprintf('%s SubType%i not supported for shell',st,typ)) %#ok<SPERR>
 end
 
 dd=[o3 o3*0;o3*0 o4];

 if il(9)~=0; o4=o4*il(9); end  % ratio for bending moment of inertia
 RO=struct('warn',{cell(4,1)});
 if ~any(il(9)==[0 1]); RO.warn{1}='12I/T^3'; end
 if il(8)&&il(8)~=ID(1) % Propagate different bending prop
      [un1,un1,un1,ddb]=p_shell('buildconstit',[il(8);ID(2:end)], ...
       varargin{3:end});
      dd(4:6,4:6)=ddb.dd(4:6,4:6); 
 end
 if il(10)&&il(10)~=ID(1) % Propagate different shear stiffness
      [un1,un1,un1,ddb]=p_shell('buildconstit',[il(10);ID(2:end)], ...
       varargin{3:end});
      ds=ddb.ds; 
 end
 if any(il(12:14))
     RO.warn{3}='Z1, Z2, MID4 not yet supported';
 end
 if any(~cellfun('isempty',RO.warn))
   fprintf('ProId %i : %s\n',ID(2),sprintf('%s, ',RO.warn{:}))
 end
 ddril=mean(dd([1 8 15]))*il(6)^2/12;% drilling DOF stiffness
 if il(4)==-1; ddril=0; elseif il(4)==0;ddril=ddril*1e-4;
 else; ddril=ddril*il(4);
 end
 % #ShellConstit : rho*h eta f(3) d h(5) k 12I/t3 nsm(8) dd(9:44) ds(45:48) ddril(49)
 %  rhoh^3/12(50)
 constit = [out3.rhoh out3.eta il([3 4 6 7 9  11 ]) dd(:)' ds(:)' ddril];
 constit(1)=constit(1)+constit(8); % Add NSM to rhoh
 out3.ConstitLab=feval(p_solid('@ConstitLab'),'p_shell.1'); 
 out3.dd=dd;out3.ds=ds;out3.ddril=ddril;
 out3.layer=[mat(1) -il(6)/2 il(6)/2 0]; % one layer
 if remi(constit(3),[],2)==1&&length(constit)==49 
   %% #Add_Rotation_Inertia
   constit(50)=constit(1)/constit(5)*constit(5).^3/12;
 end

case 2 
%% #BuildConstit.Pshell_2 Multi layer Composite specified as PCOMP entry - - - - - - - - - - - - - - - - -

r1=il(10:end);i1=ceil(length(r1)/4)*4;if length(r1)<i1;r1(i1)=0;end
r1=reshape(r1,4,length(r1)/4)';  % each row defines a ply
r1=r1(any(r1,2),:);
z=cumsum([0;r1(:,2)])+il(3);
other=[];

for j1=1:size(r1,1) % number of plies

   mat=pl(pl(:,1)==r1(j1,1),:);
   if isempty(mat);error('MatId=%i not found in pl',r1(j1,1));end
   [st,unit,typ]=fe_mat('type',mat(2));
   cz=z(j1:j1+1); 
   if any(exist(st,'file')==[2 3 6])
     [S,rhoh,r2]=feval(st,'buildply',mat,pl,il0,cz,rhoh);
     if ~isempty(r2); other=[other;r2(:)];end %#ok<AGROW>
   else; error('Ply not implemented for this type of material');
   end 
  out3.layer(j1,1:4)=[mat(1) cz(:)' r1(j1,3)]; % one layer
    %Stiffness in the principal coordinates
    C=pinv(S);
    % Stiffness matrix obtained in a arbitrary coordinate by rotation of
    % angle theta around the z axis of the main coordinates 
    % Q=(Tsigma)^-1*C*Tepsilon

    %coordinate transformation matrix
    teta=-r1(j1,3)*pi/180; ct=cos(teta); st=sin(teta);
    % ITsigma=inv(Tsigma) 
    iTsigma=[ct^2 st^2 2*ct*st 0 0;st^2 ct^2 -2*ct*st 0 0; ...
      -ct*st ct*st ct^2-st^2 0 0; 0 0 0 ct -st;0 0 0 st ct];
    Tepsilon=[ct^2 st^2 -ct*st 0 0;st^2 ct^2 ct*st 0 0; ...
      2*ct*st -2*ct*st ct^2-st^2 0 0;0 0 0 ct st;0 0 0 -st ct];
    Q=iTsigma*C*Tepsilon; % \sigma_global = Q epsilon_global
    Q1=Q([1 6 11;2 7 12;3 8 13]);
    Q2=Q([19 24;20 25]); 
    %Q1=[ct^2 st^2 -2*ct*st;st^2 ct^2 2*st*ct;st*ct -st*ct ct^2-st*2]* ...
    %    C([1 6 11;2 7 12;3 8 13])* ...
    %    [ct^2 st^2 st*ct;st^2 ct^2 -st*ct;-2*st*ct 2*st*ct ct^2-st^2];
    % 
    dd=dd+[diff(cz)*Q1 diff(cz.^2)/2*Q1;diff(cz.^2)/2*Q1 diff(cz.^3)/3*Q1];
    if typ==1; Q2=Q2*5/6;end % coincide with iso material in single layer
    ds=ds+diff(cz)*Q2;

end % loop on plies
% sdtweb t_fmesh('layers') for validation of consitency
%  if typ==1&&size(r1,1)==1;ds=ds*5/6;end % to coincide with iso material in single layer


% detect zero terms
%[i1,i2]=find(dd); r2=diag(dd); 
%dd(find(abs(dd(:)./sqrt(r2(i1).*r2(i2)))<eps));
% pay attention to limit ddril value to near coincidence to multi layer
 ddril=mean(dd([1 8 15]))*1e-4*diff(z([1 end]))^2/12;
 typ=2;
%          rho*h eta  f    d     h       k  12I/t3 nsm  dd(9) ds    dril(49)
constit = [rhoh etah il(8) 0 z(end)-z(1) 5/6 0      0  dd(:)' ds(:)' ddril];
constit(1)=constit(1)+constit(8); % Add NSM to rhoh
out3.dd=dd;out3.ds=ds;out3.ddril=ddril;
if remi(il(9),[],2)==1&&length(constit)==49 
   %% #Add_Rotation_Inertia
   constit(50)=constit(1)/constit(5)*constit(5).^3/12;
end

if ~isempty(other);constit=[constit(:);other(:)];end
case 3 
%% #BuildConstit.Pshell_3 stored configuration - - - -

%          rho*h eta  f    d     h       k  12I/t3 nsm  dd(9) ds    dril(49)
% #ShellConstit : rho*h eta f(3) d h(5) k 12I/t3 nsm(8) dd(9:44) ds(45:48) ddril(49)
i3=[8 9 10  14:16;9 11 12  17:19;10 12 13 20:22;
    14 17 20 23 24 25;15 18 21 24 26 27;16 19 22 25 27 28];
dd=il(i3); ds=il([29 30;30 31]);
out3.dd=il(i3); out3.ds=ds;out3.ddril=il(32);
constit=[il(3:7) 0 0 0 dd(:)' ds(:)' il(32)];

otherwise; error('Not a supported p_shell subtype');
end
if sp_util('diag')>1
 fprintf('\nMatId %i ProId %i, formulation %i, shell typ %i\n', ...
     ID(1:2),constit(3),typ)
 fprintf('Shell constitutive matrix membrane and bending\n')
 disp(dd)
 fprintf('Shear constitutive matrix\n')
 disp(ds)
 fprintf('Rho*h = %.3g, h= %.3g, dril = %.3g\n',constit([1 5 49]))
end

else % not a p_shell property
 [st,unit,typ]=fe_mat('typep',il(2));
 [constit,ID,out2]=feval(st,'buildconstit',ID,pl,il0);
end

if ischar(constit);out=constit;
else
  out=constit(:);     % Return CONSTIT
  ID=ID(:); if length(ID)<4;ID(4)=0;end
  out1=int32(ID(:));  % Return integ
  if isempty(out2);out2=struct('RunOpt_GroupInit','groupinitShell');end % elmap
end

% -------------------------------------------------------------------------
%% #BuildDofOpt RunOpt=p_shell('BuildDofOpt',RunOpt,pl,il);
elseif comstr(Cam,'builddofopt')

RunOpt=varargin{carg};carg=carg+1;
pl=varargin{carg};carg=carg+1;
il=varargin{carg};carg=carg+1;
RunOpt.FieldDofs=[1:6]; 
out=RunOpt;

%% SubType -------------------------------------------------------------------
elseif comstr(Cam,'subtype')

 i1=varargin{carg}; carg=carg+1;
 switch i1
 case 1;  out='Std_shell';
 case 2;  out='Comp_shell';
 otherwise; out='p_shell';
 end

%% #SetThick p_shell('setthick',model,RO,def) ------------------------------------
% RO=struct('sel','val'); def thick at nodes using DOF .19
elseif comstr(Cam,'setthick')
 model=varargin{carg};carg=carg+1;
 RO=varargin{carg};carg=carg+1;  
 RO.level=length(RO.ProId); % List of target ProId
 d1=varargin{carg};carg=carg+1;  
 [RO.iElt,elt]=feutil(['findelt' RO.sel],model);
 [EGroup,nGroup]=getegroup(elt);
 NNode=sparse(model.Node(:,1)+1,1,1:size(model.Node,1));
 d1=feutilb('placeindof',model.Node(:,1)+.19,d1);
 
 for jGroup = 1:nGroup %loop on groups to determine thickness
  [ElemF,i1,ElemP]=feutil('getelemf',elt(EGroup(jGroup),:),jGroup);
  cEGI = (EGroup(jGroup)+1:EGroup(jGroup+1)-1)';
  i1=fe_super('node',ElemP);
  i2=elt(cEGI,i1)+1;i2=full(reshape(NNode(i2),size(i2)));
  RO.thick(cEGI,1)=mean(reshape(d1.def(i2,1),size(i2)),2);
 end
 cEGI=isfinite(elt(:,1));r1=RO.thick(cEGI);
 cur=zeros(1,3);
 mpid=feutil('mpid',elt);
 il=feutil(sprintf('getil %i',mpid(2,2)),model);
 if ~isfield(RO,'pthick')
  RO.pthick=linspace(min(r1),max(r1),RO.level)';
 end
 il=repmat(il,length(RO.ProId),1);
 il(:,1)=RO.ProId(:);il(:,6)=RO.pthick(:);

 mpid(cEGI,2)=il(round((r1-RO.pthick(1))/diff(RO.pthick([1 end]))*(RO.level-1))+1,1);
 elt=feutil('mpid',elt,mpid);
 model.Elt(RO.iElt,1:size(elt,2))=elt(cEGI,:);
 out=model;
%% #SetShiftOff : transform shells with constant offset to composite ---------
elseif comstr(Cam,'setshiftoff')
 model=varargin{carg};carg=carg+1;
 [EGroup,nGroup]=getegroup(model.Elt);
 for jGroup=1:nGroup
  [ElemF,i1,ElemP]= getegroup(model.Elt(EGroup(jGroup),:),jGroup);
  if ~strncmpi(ElemP,'quad',4)&&~strncmpi(ElemP,'tria',4);continue;end
  cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
  i1=feval(ElemP,'prop');i2=feval(ElemP,'node');i1=[i1(1:2) max(i1)+i2];
  [i3,un1,i4]=unique(model.Elt(cEGI,i1),'rows');if ~any(i3(3:end)); continue;end
  for j2=1:size(i3,1)
   i2=model.il(:,1)==i3(j2,2);
   if strcmp(fe_mat('typepstring',model.il(i2,2)),'p_shell.1')
    il=p_shell(sprintf('dbval%i  laminate z0=%.15g %i %.15g 0', ...
        i3(j2,2),-i3(j2,4)+model.il(i2,6)/2,i3(j2,1),model.il(i2,6)));
    model.il(i2,:)=0;model.il(i2,1:length(il))=il;
   elseif abs((model.il(i2,3)-model.il(i2,11)/2)/i3(j2,4)+1)<1e-5
   else; error('Not yet implemented');
   end
   model.Elt(cEGI(i4==j2),i1(3:end))=0; % Remove offset
  end
 end
 out=model;
    
%% #SetTheta #Map ------------------------------------------------------------
elseif comstr(Cam,'settheta')

  model=varargin{carg};carg=carg+1;
  r1=varargin{carg};carg=carg+1;
  [eltid,model.Elt]=feutil('eltidfix;',model);
  [CAM,Cam,RunOpt.Strategy]=comstr('-strategy',[-25 1],CAM,Cam);

  if isfield(r1,'dir')
    st1='Dir must give 3 translations';
    try; if norm(r1.DOF-[.01;.02;.03])<.001; st1='';end;end
    if ~isempty(st1);error(st1);end
    r2=elem0('VectFromDirAtDofUsed',model,r1);
    r1.opt=2; r1.normal=reshape(r2.def,3,[])';
    r1.ID=fix(r2.DOF(1:3:end));RunOpt.NodeMap=1; 
  elseif isequal(r1,0); RunOpt.NodeMap=-1; % reset; 
  elseif isequal(r1.opt,2); RunOpt.NodeMap=1; % Map at Nodes
  elseif isequal(r1.opt,1); RunOpt.NodeMap=0;
  else; error('You must provide Node or Element MAP');
  end
  elt=model.Elt;[EGroup,nGroup]=getegroup(elt);
  NNode=sparse(model.Node(:,1)+1,1,1:size(model.Node,1));
  for jGroup = 1:nGroup %loop on groups
     [ElemF,i1,ElemP]=feutil('getelemf',elt(EGroup(jGroup),:),jGroup);
     i1=fe_super('node',ElemP);
     cEGI = (EGroup(jGroup)+1:EGroup(jGroup+1)-1)';
     NodePos=reshape(full(NNode(elt(cEGI,i1)+1))',length(i1),length(cEGI));

     if comstr(ElemP,'quad4')||comstr(ElemP,'quadb') || ...
        comstr(ElemP,'tria3')||comstr(ElemP,'tria6')||comstr(ElemP,'quad9')
      RunOpt.thetaPos= max(feval(ElemP,'prop'))+1;
      if size(elt,2)<RunOpt.thetaPos; elt(1,RunOpt.thetaPos)=0;end % Theta 
      if isequal(RunOpt.Strategy,2); EC=-2; else; EC=-1;end
      EC=p_shell('const',integrules(ElemP,EC));
      EC.bas=zeros(9,EC.Nw); % atcenter 
      if RunOpt.NodeMap==-1
       elt(cEGI,RunOpt.thetaPos)=0; continue;
      elseif RunOpt.NodeMap
        nind=sparse(r1.ID+1,1,1:size(r1.ID,1));
        RunOpt.Pos=...
         reshape(full(nind(elt(cEGI,i1)+1))',length(i1),length(cEGI));
          
      else; % Map at cG, r1.ID contains eltid 
          [i2,i3]=ismember(eltid(cEGI),r1.ID);cEGI(~i2)=[];
          RunOpt.Pos=i3(i2);
      end
      get_nodeE=elem0('@get_nodeE'); InfoAtNode=[];DofPos=[];
      for jElt=1:length(cEGI) %loop on elements of group
        [EC.nodeE,EC.nodeEt]=get_nodeE( ...
            model.Node,NodePos,jElt,InfoAtNode,DofPos,EC);
        of_mk('buildndn',23,EC);
        if RunOpt.NodeMap % use node map
         r3=reshape(mean(EC.bas,2),3,3)'* ...
             mean(r1.normal(RunOpt.Pos(:,jElt),:),1)'; % plane projection
         elt(cEGI(jElt),RunOpt.thetaPos)=atan2(r3(2),r3(1))/pi*180; % compute angle
        elseif ~isempty(i3)&&RunOpt.thetaPos
         r3=reshape(mean(EC.bas,2),3,3)'*r1.normal(RunOpt.Pos(jElt),:)'; % plane projection
         elt(cEGI(jElt),RunOpt.thetaPos)=atan2(r3(2),r3(1))/pi*180; % compute angle
        end
      end % loop on elements of group
     end
  end % loop on groups

out=model; out.Elt=elt;
%%
elseif comstr(Cam,'test');
    fprintf('Empty test see t_plate\n');
%% #end -----------------------------------------------------------------------
elseif comstr(Cam,'coefparam');out=[];
elseif comstr(Cam,'tablecall');out='';
elseif comstr(Cam,'cvs');
 out='$Revision: 1.139 $  $Date: 2024/11/29 08:36:51 $'; return;
else; sdtw('''%s'' not known',CAM);
end

% ---------------------------------- internal functions
% --------------------------------------------------------------------
%% #create_section_shel of standard section properties
function [out,carg] = create_section_shell(CAM,out,varg,carg)

[CAM,Cam]=comstr(CAM,1);
RunOpt=struct('Formulation',0,'Drilling',0);
i1=strfind(Cam,'-f'); 
if ~isempty(i1) % Formulation
 [RunOpt.Formulation,i2,i3,i4]=sscanf(Cam(i1+2:end),'%i',1);
 CAM=CAM([1:i1-1 i1+2+i4:end]);[CAM,Cam]=comstr(CAM,1);
end

RunOpt.unit='SI';typ=fe_mat('p_shell',RunOpt.unit,1);

if comstr(Cam,'kirchoff'); 
  if RunOpt.Formulation==0&&isempty(i1); RunOpt.Formulation=1;end
  r1=[RunOpt.Formulation RunOpt.Drilling 0 comstr(Cam,'kirchoff','%g')];
elseif comstr(Cam,'kirchhoff');  
  if RunOpt.Formulation==0&&isempty(i1); RunOpt.Formulation=1;end
  r1=[RunOpt.Formulation RunOpt.Drilling 0 comstr(Cam,'kirchhoff','%g')];
elseif  comstr(Cam,'mindlin')
  r1=[RunOpt.Formulation RunOpt.Drilling 0 comstr(Cam,'mindlin','%g')]; 
elseif  comstr(Cam,'laminate')
    
 [CAM,Cam,r2]=comstr('z0',[-25 2 1],CAM,Cam);
 r1=comstr(Cam,'laminate','%g');
 r1=reshape(r1,3,length(r1)/3);r1(4,:)=0;
 z=cumsum([0 r1(2,:)]); if isempty(r2);r2=-z(end)/2;end
 r1=[r2 0 0 0 0 0 0 r1(:)'];
 typ=fe_mat('p_shell','SI',2);
elseif  comstr(Cam,'honeycomb') % - - - - - - - - - - - - - - - - - -

[CAM,Cam]=comstr(CAM,10); opt=comstr(CAM,-1);
mat=varg{carg}; carg=carg+1;
if length(opt)<9; opt(9)=0; end % use G1z_min;

%opt=[hpeau1,hpeau2,hcellule,a,b,teta,t,t1,alpha_minmax]
%	Subtype 5 : Orthotropic material for shell
%           [MatId type E1 E2 nu12 G12 G1z G2z Rho A1 A2 TREF Xt Xc Yt Yc S ...
%             Ge F12 STRN]
%           type  fe_mat('m_elastic','SI',4)

[st,RunOpt.unit,typ]=fe_mat('typem',mat.pl(1,2));
if ~strcmp(st,'m_elastic')||(typ~=1)
 error('Honeycomb supposes an elastic material');
end

o6=opt(6)/180*pi;

r1=[ ...
 mat.pl(3)*(opt(7)/opt(5))^3*((opt(4)/opt(5)+sin(o6))/(cos(o6))^3); %E1
 mat.pl(3)*(opt(7)/opt(5))^3*cos(o6)/((opt(4)/opt(5)+sin(o6))*sin(o6)^2); %E2
 ((opt(4)/opt(5)+sin(o6)))*sin(o6)/cos(o6)^2; % nu12
 mat.pl(3)*(opt(4)/opt(5)+sin(o6))/(2*(opt(4))^2*cos(o6))*1/(opt(5)/(2*(opt(8))^3)+opt(4)/opt(7)^3); % G12
mat.pl(6)*(opt(4)*opt(8)+opt(5)*opt(7)*sin(o6))^2/(opt(5)*cos(o6)*(opt(4)+opt(5)*sin(o6))*(2*opt(4)*opt(8)+opt(5)*opt(7))) ; %G1z
 mat.pl(6)*cos(opt(6))/(opt(4)/opt(5)+sin(opt(6)))*opt(7)/opt(5)]'; % G2z

r1(end-1)=r1(end-1)*(1-opt(9))+opt(9)* ...
 mat.pl(6)*(opt(4)*opt(8)+2*opt(5)*opt(7)*(sin(o6))^2)/ ...
  (2*opt(5)*cos(o6)*(opt(4)+opt(5)*sin(o6)));
% voir comment choisir G1z en fonction de G1z_inf et G1z_sup.

 typ=fe_mat('p_shell',mat.unit,5);
%% #Merge : combine multiple laminates ProIdBase +/- ProidNext1 +- ...
elseif  comstr(Cam,'merge')
 model=varg{carg};carg=carg+1;
 i1=comstr(CAM(6:end),-1);
 RO=feutil('getdd',[0 i1(1)],model);
 if ~isfield(RO,'layer'); % Not a composite def 
  mpid=feutil('mpid',model);mpid=unique(mpid(:,1:2),'rows');
  mpid(mpid(:,2)~=i1(1),:)=[]; 
  if size(mpid,1)>1; error('Not a valid case');end
  RO=feutil('getdd',mpid,model);
 end
 for j1=2:length(i1)
  dd=feutil('getdd',[0 abs(i1(j1))],model);
  if i1(j1)>0; % Add to top
    dd.layer(:,2:3)=dd.layer(:,2:3)+RO.layer(end,3)-dd.layer(1,2);
    RO.layer=[RO.layer;dd.layer];
  else
    dd.layer(:,2:3)=dd.layer(:,2:3)+RO.layer(1,2)-dd.layer(end,3);
    RO.layer=[dd.layer;RO.layer];
  end
 end
 il=feutil(sprintf('getil %i',i1(1)),model);
 [m_fun,unit,substr,subi]=fe_mat('typep',il(1,2));
 if subi==1 %transform to composite
  il=[il(1) fe_mat('type',m_fun,unit,2) RO.layer(1,2) zeros(1,6)];
 end
 r1=[RO.layer(:,1) RO.layer(:,3)-RO.layer(:,2) RO.layer(:,4)*[1 0]]';
 il(9+(1:numel(r1)))=reshape(r1,[],1); % layers
 il(3)=RO.layer(1,2); % z0
 r1=il(3:end); typ=il(2);
 
else 
  warning('Using default shell property'); out=out(1); %#ok<WNTAG>
  return
end

 out=struct('il',[1 typ r1],'name',Cam,'type','p_shell', ...
   'unit',RunOpt.unit);

% EOF
