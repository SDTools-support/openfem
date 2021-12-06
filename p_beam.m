function [out,out1,out2,out3]=p_beam(varargin)

%P_BEAM Beam element property function
%
%       Syntax : pro= p_beam('default') 
%                pro= p_beam('database','name') 
%                il = p_beam('dbval ProId name');
%                il = p_beam('dbval -unit MM ProId name');
%                il = p_beam('dbval -punit MM ProId name val');
%
%       Subtype 1 : standard section
%               [ProId   type   J I1 I2 A   k1 k2 Lump NSM]
%          with
%	     type : fe_mat('p_beam','SI',1)
%	     J    : torsional stiffness parameter (often different from polar
%		    moment of inertia I1+I2)
%	     I1   : moment of inertia for bending plane 1 
%           see vx vy vz in beam1 documentation of details 
%	     I2   : moment of inertia for bending in the second plane which
%		    contains the beam and is orthogonal to plane 1
%        A    : section area
%        k1   : (optional) shear factor for direction 1 (when not 0 a
%		     Timoshenko beam element is used)
%        k2   : (optional) shear factor for direction 2
%        lump : 1 for lumped mass, 0 otherwise
%        NSM  : non structural mass (density of mass per unit length).
%
%        pro=p_beam('database','name') : returns struct array
%         for reference beam sections. If 'name' is provided a match
%         in the section database is made. The database supports
%         standard sections
%                    circle (r)      : full circular section
%                    rectangle (b h) : full rectangle
%         For example p_beam('database circle .1')
%
%       Subtype 3 : predefined section with geometric parameters
%
%       See sdtweb     p_beam, p_shell, fem (handling materials), pl, fe_mat
%       See also help  fe_mat

%	Etienne Balmes, Jean-Philippe Bianchi
%       Copyright (c) 2001-2020 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*NOSEM>

if nargin<1; help p_beam;return; end
obj=[];evt=[];
if ischar(varargin{1}); [CAM,Cam]=comstr(varargin{1},1); il=[]; carg=2;
elseif ~ischar(varargin{1})&&nargin==3&&ischar(varargin{3})
  obj=varargin{1};evt=varargin{2};[CAM,Cam]=comstr(varargin{3},1); carg=4;
else; il=varargin{1};[CAM,Cam]=comstr(varargin{2},1); carg=3;
end

if sp_util('diag');fprintf('p_beam %s \n',CAM);end

%% #propertyunittype ---------------------------------------------------------
if comstr(Cam,'propertyunittype')||comstr(Cam,'pu')

 if nargin==1; out=[1 3 4];return;end % return subtypes ID
 i1=varargin{carg}; 
 out1={};
 switch i1
 case 1 % 
     %     [ProId   type   J I1 I2 A   k1 k2]
     st={... 
         'ProId'  0      'sdtweb(''p_beam'')';
         'Type'   0      '';
         'J'      11     'Torsion';
         'I1'     11     'Bending'; 
         'I2'     11     'Bending';
         'A'      12     'Area';
         'k1'     0      'Shear factor';
         'k2'     0      'Shear factor';
         'lump'   0      'Lumped mass model'
         'NSM'    9.004  'Non structural mass density';
         'C1'     4      'Stress recovery y1'
         'C2'     4      'Stress recovery z1'
         'D1'     4      'Stress recovery y1'
         'D2'     4      'Stress recovery z1'
         'E1'     4      'Stress recovery y1'
         'E2'     4      'Stress recovery z1'
         'F1'     4      'Stress recovery y1'
         'F2'     4      'Stress recovery z1'
         'I12'    11      'Area product of inertia (unused)'
         'Jw'     11      'wrapping coefficient (unused)'
         'CGy',4,'unused'
         'CGz',4,'unused'
         'Shy',4,'unused'
         'Shz',4,'unused'
         % CGy,CGz,Shy,Shz
         };
 %case 2 % NASTRAN PBEAM 
 case 3 % NASTRAN PBEAML : GROUP TYPE VALUES
 % Lumped and NSM are not supported here due to incompatibility with
 % convert, eventually new subtype should be implemented
       st={... 
         'ProId'   0   'sdtweb(''p_beam'')';
         'Type'    0   '';
         'Group'  -1  '';
         'Section' -1  'comstr(''SectionName'',-32)';
         %'' Inf ''; % empty line as Dim%i has to be inserted here
         %'NSM'    9.004  'Non structural mass density';
         %'lump'   0      'Lumped mass model'
         }; 
       out1={'Dim%i'    4  'Dimension'};

 case 4 % SAP2000 tapered beam
  %- Local Coordinate is from 166 to 171 ( PDF page )
  %- Non-prismatic is from 177~180
    
       st={... 
         'ProId'        0 'sdtweb(''p_beam'')';
         'Type'         0 'fe_mat(''p_beam'',''SI'',4)';
         'StartSec'     0 'Index of start section'
         'EndSec'       0 'Index of end section'
         'SecVarType1'  0 '1 linear, 2 parabolic, 3 cubic'
         'SecVarType2'  0 '1 linear, 2 parabolic, 3 cubic'
         }; 
 case 5 % Ansys pipe properties
  st={         'ProId'        0 'sdtweb(''p_beam'')';
         'Type'         0 'fe_mat(''p_beam'',''SI'',5)';
     'Do',0,'';'Tw',0,'';'Nc',0,'';'Ss',0,'';'Nt',0,'';'Mint',0,'';'Mins',0,'';'Tins',0,''};
 otherwise; st={'ProId' 0 'sdtweb(''p_beam'')'; 'Type', 0, ''};
 end
 if ~isempty(strfind(Cam,'cell')); 
    if nargout==0;disp([num2cell((1:size(st,1))') st]);else;out=st;end
 else; out=[st{:,2}]; 
 end

%% #Default ------------------------------------------------------------------
elseif comstr(Cam,'default');

 if ~isempty(strfind(Cam,'smart'))
  model=evalin('caller','model');
  node=feutil('getnode groupall',model);
  r1=mean(mean(diff(sortrows(node(:,5:7)))))/100;
  out=p_beam(sprintf('dbval 1 circle %g',r1));
 else
  out=p_beam('database');out=out(1);
 end

%% #Dbval --------------------------------------------------------------------
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
   i5=i2+[0:5+i5];CAM(i5)=''; [CAM,Cam]=comstr(CAM,1);
  else;PUnit='';
  end
  [CAM,Cam]=comstr(CAM,6);
  [i1,CAM,Cam]=comstr(CAM,'','%i');
  if isempty(CAM)&&carg<=nargin; st=varargin{carg};carg=carg+1;
  else;st=CAM;end
  if ischar(st); mat=p_beam('database',st);
  elseif isnumeric(st)
   [typ,st1,i4]=fe_mat('typep',st(2));
   mat=struct('il',st,'name',sprintf('%i',st(1)),'type',typ,'unit',st1);
  end
  if ~isempty(PUnit)
   if ~strncmpi(PUnit,'US',2)
    r1=fe_mat(sprintf('convert %s %s',mat.unit,PUnit),mat.il(1:2));
    mat.il(2)=r1(2);
   end
   mat.unit=PUnit;
  end
  if ~isempty(Unit)
   mat.il=fe_mat(sprintf('convert %s %s',mat.unit,Unit),mat.il);mat.unit=Unit;
  end
  r1=mat.il; if length(i1)==1; r1(1)=i1;end
  if ~isempty(il); i2=find(il(:,1)==r1(1)); else;i2=[];end
  if isempty(i2); i2=size(il,1)+1; else; il(i2,:)=0; end %#ok<AGROW> % do not forget to clean up !
  il(i2,1:length(r1))=r1; %#ok<AGROW>
  if carg>nargin; break;end
  CAM=varargin{carg};carg=carg+1;if ischar(CAM);[CAM,Cam]=comstr(CAM,1);end
 end
 out=il;

%% #Database -----------------------------------------------------------------
elseif comstr(Cam,'database') 
  
  st=comstr(CAM,9);
  if ~isempty(st)
  elseif carg<=nargin; st=varargin{carg}; carg=carg+1;
     if ~ischar(st); 
         st=sprintf('%i %s %s',st(1),comstr(st(4),-32), ...
             sprintf('%.15g',st(5:end)));
     end
  end
  [MatId,i2,i3,i4]=sscanf(st,'%i',1); if i2~=1; MatId=1;end
  st=lower(comstr(st,i4));

  st1='circle .01 m';if comstr(st,'circle');st1=st;end
  out =feutil('rmfield',create_section_beam(st1,[]),'sectionmdl'); 
  out.il(1)=MatId; 
  st1='rectangle .1 .2 m'; if comstr(st,'rect'); st1=st;end
  out(2)=feutil('rmfield',create_section_beam(st1,[]),'sectionmdl'); 
  out(2).il(1)=MatId;
  %out(3) = struct('sectionmdl',[],'il',[MatId fe_mat('p_beam','SI',3) 0 comstr('ROD',-32) 1.1],'name','reftube','type','p_beam','unit','SI');
  out1='Beam section';
  i1=strmatch(sscanf(regexprep(st,'[.0123456789]',''),'%s',1),{out.name}); % find in dat
  if length(i1)==1; out=out(i1); return; end
  if ~sp_util('issdt') % OpenFem
   return % circle and rectangle
  end
  % Build all subtype 3
  list=create_section_beam('',[]); % list of all subtype3 section (nastran type)
  %i1=strmatch(sscanf(regexprep(st,'[.0123456789]',''),'%s',1),list(:,1),'exact'); % find in NASTRAN type database
  i1=[];
  if ~isempty(st); 
   for j1=1:size(list,1); if comstr(st,list{j1,1}); i1=j1; break; end; end
  end
  if isempty(i1); % return the whole database
   if ~isempty(st); sdtw('_nb','Can''t find %s in database',st); end
   i1=1:size(list,1); 
  end
  for j1=i1 % handle p_beam('info') predefined properties
    r2=list{j1,2}; 
    if ~isempty(r2)
     if any(st=='=');
      [r2,st1,st]=cingui('paramedit -DoClean',r2,{struct,st});st1=lower(st);
      st1=fieldnames(r2); 
     else;r2=fe_def('cleanentry',cingui('paramedit',r2));
      [st,st1,r3]=comstr(list{j1,1},[-25 2],st,st); % get values
      st1=fieldnames(r2); 
      for j2=1:min(length(st1),length(r3)); r2.(st1{j2})=r3(j2);end % Set values provided
     end
     [st2,i2]=setdiff(st1,{'NSM','lump'}); st1=st1(sort(i2)); % critical to keep fields ordered
     r3=struct(...
         'il',[MatId fe_mat('p_beam','SI',3) 0  ...
         comstr(list{j1,1}(1:min(5,length(list{j1,1}))),-32) ...
         zeros(1,length(st1))], ...
         'name',list{j1,1},'type','p_beam','unit','SI');
     for j2=1:length(st1); r3.il(4+j2)=r2.(st1{j2});end
     if r2.NSM||r2.lump % p_beam('pu-cell',1)
      r3.il=p_beam('convertTo1',r3.il);
      r3.il(9)=r2.lump; r3.il(10)=r2.NSM;
     end
    else;   r3=feutil('rmfield',create_section_beam(list{j1,1},[]),'sectionmdl');
    end
    if length(i1)==1; out=r3; else; out(j1+2)=r3; end
  end 
%% #constit ->sdtweb p_solid('BuildBeam')
   % ConstitLab={'E','nu','Rho','G','eta','alpha','T0','J','I1','I2','A','k1','k2','lump'}
%% #EltConst ->sdtweb beam1t('BasisAndConstants')
elseif comstr(Cam,'buildconstit') % 'buildconstit',ID,pl,il,model
   ElemF='beam1'; % not bar constit
   if nargout==4;
    [out,out1,out2,out3]=p_solid('buildconstit',varargin{2:end});
   else;[out,out1,out2]=p_solid('buildconstit',varargin{2:end});
   end
elseif comstr(Cam,'xxx') % 'buildconstit',ID,pl,il,model
   if nargout==4;
    [out,out1,out2,out3]=beam1t('basisandconstants',varargin{2:end});
   else;[out,out1,out2]=beam1t('basisandconstants',varargin{2:end});
   end
   dbstack;keyboard;
%% #BoxFromIl ----------------------------------------------------------------
elseif comstr(Cam,'boxfromil') 
 % find equivalent box section from I1 I2 and A
 il=varargin{carg}; carg=carg+1;
 [st,i1,i2]=fe_mat('type',il(2));
 if ~comstr(st(3:end),'beam'); sdtw('_err','Not a beam il'); end
 if i2==3 % nastran section : get il
  il=p_beam('ConvertTo1',il);
 elseif i2==2||i2>3
  sdtw('_err','Not available type')  
 end
 if il(4)~=il(5)
  sdtw('_nb','Something to be done to take in account I1 and I2 : XXX')
 end
 %out=[il(1) fe_mat('p_beam','SI',3) sqrt(A) sqrt(A)];
 out=p_beam(sprintf('database 1 BAR %.15g %.15g',sqrt(il(6)),sqrt(il(6)))); 

 %keyboard
 % voir nastran p90, donnant
 
%% #ConvertTo1 p_beam('ConvertTo1;',mo2.il) or p_beam('ConvertTo1;',sec) -----
elseif comstr(Cam,'convertto'); [CAM,Cam]=comstr(CAM,10);
 
 RunOpt.typ=comstr(CAM,[-1 1]);
 RunOpt.Silent=~isempty(Cam)&&CAM(end)==';';
 if ~isempty(strfind(Cam,'silent'));RunOpt.Silent=1;end
 %[CAM,Cam,RunOpt.struct]=comstr('-struct',[-25 3],CAM,Cam);
 il=varargin{carg}; carg=carg+1;
 if isfield(il,'Elt') 
   % A simple procedure to guess properties from arbitrary sections. 
   % Can be grossly overestimated for torsion.
   mo1=il; mo1.Elt=feutil('set groupall matid 1 proid 1',mo1);
   mo1.il=p_shell('dbval 1 mindlin 1e-6');
   mo1.pl=[1 fe_mat('m_elastic','US',1) 0 0 1];
   r1=feutilb('geomrb -cg',mo1,[0 0 0]);
   st='US';if isfield(mo1,'unit');st=mo1.unit;end
   out=[1 fe_mat('p_beam',st,1) ...
          flipud(sort(diag(r1.inertia)))'*1e6 r1.mass*1e6]; %#ok<FLPST>
   return
 elseif isfield(il,'il'); il=il.il; 
 end

 out1=[];out=il;
 for j1=1:size(il,1)
  [st,i1,i2]=fe_mat('type',il(j1,2));
  if ~comstr(st,'m_beam')
   if ~RunOpt.Silent;
       sdtw('_nb','ProId %i is not a p_beam property',il(j1,1)); 
   end
   continue;
  end
  if i2==2||i2>3
   sdtw('_err','ProId %i SubType %i is unknown',il(j1,1),i2)
  end
  if i2==RunOpt.typ
   if ~RunOpt.Silent
    sdtw('_nb','ProId %i is already SubType %i',il(j1,1),RunOpt.typ)
   end
   continue
  end
  if RunOpt.typ==1     % n to 1 - - - - - - - - - - - - - - - - - - - -
   r1=create_section_beam(sprintf('%s %s',comstr(il(j1,4),-32),...
                                 sprintf('%.15g ',il(j1,5:end))),[]);
   % /!\ make sure that values over il are cleared not to interfere with type 1 entries
   r1.il(1)=out(j1,1);out(j1,:)=0; 
   out(j1,1:size(r1.il,2))=r1.il; % subtype 1 il
   if isfield(r1,'sectionmdl'); out1=r1.sectionmdl; end
  elseif RunOpt.typ==3 % n to 3 - - - - - - - - - - - - - - - - - - - - 
  'boxfromil call XXXJP'
  else
   sdtw('_err','SubType %i conversion not available',RunOpt.typ)
  end
 end% 
   
%% #Info ---------------------------------------------------------------------
elseif comstr(Cam,'info') 
 
 out=create_section_beam('',[]);
if nargout==0; out=out';fprintf('%10s   %s\n',out{:});clear out;end
% -------------------------------------------------------------------------
elseif comstr(Cam,'subtype')

 i1=varargin{carg}; carg=carg+1;
 if ischar(i1); out=1;
 else
  switch i1
  case 1;  out='p_beam';
      otherwise; out='p_beam';
  end
 end
 
%% #Show3D -------------------------------------------------------------------
elseif comstr(Cam,'show3d') % callback for gui, full mdl 3D display

 Sel=[];
 if ~isempty(obj)
    [uf,gf]=clean_get_uf('feplotpro',obj);
    cf=clean_get_uf('feplotcf',gf);
    Sel=cf.sel;
    if ~isfield(Sel,'mdl')||~isfield(Sel.mdl,'name')|| ...
            ~strcmp(Sel.mdl.name,'Show3D');Sel=[];
    end
 elseif carg<=nargin&&isa(varargin{carg},'sdth') %p_beam('show3d',cf);
    cf=varargin{carg};carg=carg+1;[uf,gf]=clean_get_uf('feplotpro',cf,'open');
    Sel=cf.sel;
    if ~isfield(Sel,'mdl')||~isfield(Sel.mdl,'name')|| ...
            ~strcmp(Sel.mdl.name,'Show3D');Sel=[];
    end
 else
  try; 
   ga=evalin('caller','ga');ua=evalin('caller','ua');
   uf=evalin('caller','uf');
  catch
   if isempty(gcbo)||strcmp(get(gcbo,'tag'),'test_sdt')
    cf=feplot;[uf,gf]=clean_get_uf('feplotpro',cf,'open');
   else
    uo=get(gcbo,'userdata'); ga=uo.CallingAxes;ua=get(ga,'userdata');
    uf=clean_get_uf(ga);
   end
  end
 end
 if carg<=nargin&&~isfield(varargin{carg},'def');
     RO=varargin{carg};carg=carg+1;
 else; RO=struct;
 end
 if ~isempty(strfind(Cam,'reset'))||(isfield(RO,'reset')&&RO.reset);Sel=[];end
 if ~isfield(RO,'sel');RO.sel=1;end
 if isempty(Sel) % Init model 
  mdl=fe_fmesh(['volbeam' CAM],uf.mdl,RO); % get visualisation mdl
  mdl.name='Show3D';
  if isempty(mdl.Node); return; end
  Sel=feutil('GetPatchNew',mdl); Sel.mdl=mdl; Sel.opt(1,2)=sdtdef(1);
 end
 if carg<=nargin&&isfield(varargin{carg},'def')
  cf.vfields.SelF{RO.sel}=Sel;
  def=varargin{carg};carg=carg+1; 
  def=feutilb('placeindof',Sel.mdl.TR.adof,def);
  def.TR=Sel.mdl.TR;cf.def=[];cf.def=def;
 elseif ~isempty(RO.sel);
  cf.vfields.SelF{RO.sel}=Sel;
 else
  h=findall(uf.ParentFigure,'tag','fepro_stacksel');
  st=get(h,'state');st1=get(h,'OffCallback');st2=get(h,'OnCallback');
  eval('DisplaySel=feutilg(''@DisplaySel'');')
  set(h,'OffCallback','','OnCallback','','state','on');
  obj=DisplaySel(uf,struct('Sel',Sel));
  set(h,'state',st,'OffCallback',st1,'OnCallback',st2); 
  %set(h,'state','Off');
  %cf=uf.FeplotHandle; cf.vfields.SelF{1}=Sel; feplot
 end
 if ~isempty(obj);
     delete(findobj(cf.ga,'tag','StackSel'));feplot(cf);
     figure(cf.opt(1));iimouse('resetview')
 end
 if nargout>0; out=Sel; end

%% #showmap --------------------------------------------------------------
elseif comstr(Cam,'showmap') % callback for gui

   cf=feplot;
   
%% #showsection --------------------------------------------------------------
elseif comstr(Cam,'showsection') % callback for gui
 
 if ~isempty(obj)
  [uf,gf]=clean_get_uf('feplotpro',obj);cf=clean_get_uf('feplotcf',gf);
  [pro,il]=matgui('getstackil',uf.mdl);pro=pro(1,:);
 elseif carg<=nargin&&isa(varargin{carg},'sdth')
  cf=varargin{carg};carg=carg+1;uf=clean_get_uf('feplotpro',cf);
  [pro,il]=matgui('getstackil',uf.mdl);pro=pro(1,:);
 else
  try; 
   ga=evalin('caller','ga');ua=evalin('caller','ua');uf=evalin('caller','uf');
  catch
   uo=get(gcbo,'userdata'); ga=uo.CallingAxes;ua=get(ga,'userdata');
   uf=clean_get_uf(ga);
  end
  [pro,il]=matgui('getstackil',uf.mdl,ua.table{2,2}.value);
 end
 
 name=pro{2};pro=pro{3};pro.il=il;
 [il,sectionmdl]=p_beam('ConvertTo1',pro.il); % get visualisation mdl
 if isempty(sectionmdl);warning('SDT:NoSection','No section model for\n %s', ...
         feutilb('_writeil',pro.il))
     return;
 end
 pro.Sel=feutil('GetPatchNew',sectionmdl);pro.Sel.mdl.name='BeamSection';
 cf=uf.FeplotHandle;cf.mdl=stack_set(cf.mdl,'pro',name,pro);
 ga=cf.ga;cla(ga);
 [sel,uo]=feutilg('initpatch_object',pro.Sel,ga);
 set(uo(uo~=0),'tag','StackSel');
 axes(ga);axis auto; fecom view2 

%% #StressObserve ------------------------------------------------------------
% sdtweb p_beam EltConst % for EltConst init
elseif comstr(Cam,'stressobserve') 
 
 model=varargin{carg}; carg=carg+1;
 if isfield(model,'MatId')
  mat=model; pro=varargin{carg};carg=carg+1;
 else
  ID=varargin{carg}; carg=carg+1; % ID
  [pro,il]=matgui('getstackilS',model,ID(2));pro=pro{3};
  mat=feutil(sprintf('GetPl%i -struct1',ID(1)),model);
  %pr=feutil(sprintf('GetIl%i -struct1',ID(2)),model);
 end
 
 % get default Stress points if needed
 if ~isfield(pro,'StressOut');error('This should not occur');end
 if ~isfield(pro.StressOut,'Node')
  % Pre-defined A,B,C,D Nastran stress recovery points
  [il,sectionmdl]=p_beam('ConvertTo1',pro.il); % get visualisation mdl
  if isempty(sectionmdl)
   sdtw('_nb','No pre-defined stress points for ProId %i beams.',pro.il(1))
   out=[]; return
  else
   n1=feutil('GetNode setname StressRP',sectionmdl);
   pro.StressOut=struct('Node',n1,...
                           'lab',{{'C';'D';'E';'F'}});
  end
 end
 % Build observation matrix
 n1=pro.StressOut.Node;
 st=cellfun(@(x)sprintf('Sxx_%s',x),pro.StressOut.lab,'uni',0);
 out=struct('cta',[],'X',{{st}},'Xlab',{{'Component'}},'Node',n1);
 % C1D.X{1} : 'e_xx'    'k_y'    'k_z'    's_y'    's_z'    'phi_x'
 % beam1 : beam axis = x
 % p_beam : beam axis = z, 1st bending plane = y-z
 
 out.cta=[ones(size(n1,1),1) -mat.E*n1(:,6) mat.E*n1(:,7) zeros(size(n1,1),3)];


%% #End ----------------------------------------------------------------------
elseif comstr(Cam,'tablecall');out='';
elseif comstr(Cam,'cvs')
    out='$Revision: 1.115 $  $Date: 2021/11/29 11:03:45 $';
elseif comstr(Cam,'@');out=eval(CAM);
else; sdtw('''%s'' not known',CAM); %error('''%s'' not known',CAM);
end
%% #p_get_prop --------------------------------------------------------------
function [r1,out1]=p_get_prop(r1,st)

r2=fe_def('cleanentrycell',cingui('paramedit',st{2}));
st1=r2(:,1);
if length(r1)<length(st1)
   r3=horzcat(r2{:,2});r4=r3(length(r1)+1:length(r3));
   if any(r4)
    sdtw('_nb','%s section missing %s, using default',upper(st{1}), ...
       comstr(st1(length(r1)+find(r4)),-30))
   end
   r1(length(r1)+1:length(r3))=r4;
end
if nargout==2
 out1=r2(:,1:2)';out1(2,:)=num2cell(r1(1:size(out1,2)));out1=struct(out1{:});
end

%% #create_section_beam(CAM,) ------------------------------------------------------
% Creation of standard section properties
function out = create_section_beam(CAM,out)

 [CAM,Cam]=comstr(CAM,1);
 %   'circle' 'r(.1#%g#"Radius")'
 %   'rectangle','b(.1#%g#"width") h(.2#%g#"height")'
    
 st={'rod','r(1#%g#"Radius") NSM(0#%g#"NSM") lump(0#%g#"Lump")'
     'tube' 'r2(1#%g#"Outer radius") r1(.5#%g#"Inner radius") NSM(0#%g#"NSM") lump(0#%g#"Lump")'
     't2','d1(1#%g#"DIM1") d2(2#%g#"DIM2") d3(.2#%g#"DIM3") d4(.4#%g#"DIM4") NSM(0#%g#"NSM") lump(0#%g#"Lump")'
     't','d1(1#%g#"DIM1") d2(2#%g#"DIM2") d3(.2#%g#"DIM3") d4(.4#%g#"DIM4") NSM(0#%g#"NSM") lump(0#%g#"Lump")'
     'i','d1(1#%g#"DIM1") d2(1#%g#"DIM2") d3(1#%g#"DIM3") d4(.2#%g#"DIM4") d5(.2#%g#"DIM5") d6(.2#%g#"DIM6") NSM(0#%g#"NSM") lump(0#%g#"Lump")'
     'l','d1(.6#%g#"DIM1") d2(1#%g#"DIM2") d3(.2#%g#"DIM3") d4(.1#%g#"DIM4") NSM(0#%g#"NSM") lump(0#%g#"Lump")'
     'bar','b(1#%g#"width") h(2#%g#"height") NSM(0#%g#"NSM") lump(0#%g#"Lump")'
     'box','w(1#%g#"width") h(1.2#%g#"height") tv(.2#%g#"vertical thickness") th(.1#%g#"h thickness") NSM(0#%g#"NSM") lump(0#%g#"Lump")'
     'hat','d1(.5#%g#"height") d2(.5#%g#"thick") d3(.8#%g#"top") d4(1#%g#"edge width") NSM(0#%g#"NSM") lump(0#%g#"Lump")'
     'chan1','d1(.5#%g#"DIM1") d2(.5#%g#"DIM2") d3(.8#%g#"DIM3") d4(1#%g#"DIM4") NSM(0#%g#"NSM") lump(0#%g#"Lump")'
     'chan2','d1(.1#%g#"DIM1") d2(.1#%g#"DIM2") d3(1#%g#"DIM3") d4(1#%g#"DIM4") NSM(0#%g#"NSM") lump(0#%g#"Lump")'
     'chan','d1(.5#%g#"DIM1") d2(1#%g#"DIM2") d3(.1#%g#"DIM3") d4(.2#%g#"DIM4") NSM(0#%g#"NSM") lump(0#%g#"Lump")'};

if isempty(Cam); out=st; return; end
RS=[];

if comstr(Cam,'circ') % - - - - - - - - - - - - - - - - - - - - - - - CIRCLE

  r1=comstr(Cam,'circle','%g');
  if isempty(r1); r1=comstr(Cam,'circl','%g'); end
  if isempty(r1); r1=0.1; end
  out=[];
  out.sectionmdl=[];
  out.il=[1 fe_mat('p_beam',1,1) pi*[1/2 1/4 1/4]*r1(1)^4 pi*r1(1)^2 0 0];
  out.name=Cam;
  out.type='p_beam';
  out.unit='SI';

elseif comstr(Cam,'rect') %  - - - - - - - - - - - - - - - - - - - RECTANGLE 

  r1=comstr(Cam,'rectangle','%g');
  if isempty(r1); r1=comstr(Cam,'recta','%g'); end
  if isempty(r1); r1=[0.1 0.2]; end
  if length(r1)<2
    error('You must specify 2 dimensions for rectangle beam section');
  end
  out=[];
  out.sectionmdl=[];
  b=r1(1); h=r1(2); 
  out.il=[1 fe_mat('p_beam',1,1) 0 b*h^3/12 b^3*h/12 b*h 0 0];
  a=max(b,h)/2;b=min(b,h)/2;
  out.il(1,3)=a*b^3*(16/3-3.36*b/a*(1-b^4/(12*a^4))); % from wikipedia
  out.name=Cam;
  out.type='p_beam';
  out.unit='SI';

elseif comstr(Cam,'rod') % - - - - - - - - - - - - - - - - - - - - - - - - ROD
 % special visu parameter for sectionmdl mesh
 [CAM,Cam,nrev]=comstr('nrev',[-25 1],CAM,Cam); if isempty(nrev); nrev=8; end
 % 2004 nastran reference manual p90 nastran.pdf
 [r1,RS]=p_get_prop(comstr(Cam,'rod','%g'),st(strcmpi('rod',st(:,1)),:));
 if length(r1)<1; sdtw('_nb','ROD section needs 1 dimension'); r1=[1]; end
 mdl=[]; %no shape prepared
 mdl.Node=[1 0 0 0 0 0 0; 2 0 0 0 r1(1) 0 0];
 mdl.Elt=[Inf abs('beam1'); 1 2 1 1 0 0];
 mdl=feutil(sprintf('rev %i o 0 0 0 360 0 0 1',nrev),mdl);
 out=struct('sectionmdl',mdl,'il',[1 fe_mat('p_beam',1,1)  ...
     pi*[1/2 1/4 1/4]*r1(1)^4 pi*r1(1)^2 .9 .9],...
     'name',Cam,'type','p_beam','unit','SI');

elseif comstr(Cam,'tube') % - - - - - - - - - - - - - - - - - - - - - - - TUBE

 [r1,RS]=p_get_prop(comstr(Cam,'tube','%g'),st(strcmpi('tube',st(:,1)),:));
 
 if length(r1)<1; sdtw('_nb','TUBE section needs 2 dimensions'); r1=[1 .5]; end
 if r1(2)>r1(1); fprintf('%s does not verify ro<ri, switching\n',CAM);
     r1(1:2)=r1([2 1]);
 end
 mdl.Node=[1 0 0 0 r1(2) 0 0; 2 0 0 0 r1(1) 0 0];
 mdl.Elt=[Inf abs('beam1'); 1 2 1 1 0 0];
 mdl=feutil('rev 8 o 0 0 0 360 0 0 1',mdl);
 out=struct('sectionmdl',mdl,...
     'il',[1 fe_mat('p_beam',1,1) pi*[1/2 1/4 1/4]*(r1(1)^4-r1(2)^4) ...
     pi*(r1(1)^2-r1(2)^2) .5 .5],...
     'name',Cam,'type','p_beam','unit','SI');
 
elseif comstr(Cam,'t2') % - - - - - - - - - - - - - - - - - - - - - - - - - T2
 % 2004 nastran reference manual p90 nastran.pdf
 %r1=comstr(Cam,'t2','%g');
 [r1,RS]=p_get_prop(comstr(Cam,'t2','%g'),st(strcmpi('t2',st(:,1)),:));
 
 if length(r1)<4; sdtw('_nb','T2 section needs 4 dimensions'); r1=[1 2 .2 .4]; end
 out=create_section_beam([CAM([1 3:end]) sprintf('%.15g ',r1)],out);
 out.sectionmdl=feutil('RotateNode o 0 0 0 180 0 0 1',out.sectionmdl);
 out.name='t2';
 
elseif comstr(Cam,'t') % - - - - - - - - - - - - - - - - - - - - - - - - - - T
 % 2004 nastran reference manual p90 nastran.pdf
 [r1,RS]=p_get_prop(comstr(Cam,'t','%g'),st(strcmpi('t',st(:,1)),:));
 
 if length(r1)<4; sdtw('_nb','T section needs 4 dimensions'); r1=[1 2 .2 .4]; end
 node=[-r1(1)/2 r1(3)/2 0; r1(1)/2 r1(3)/2 0;
     r1(1)/2 -r1(3)/2 0; -r1(1)/2 -r1(3)/2 0];
 mdl=feutil('Objectquad 1 1',node,1,1);
 node=[-r1(4)/2 -r1(3)/2 0; r1(4)/2 -r1(3)/2 0;
     r1(4)/2 -r1(2)+r1(3)/2 0; -r1(4)/2 -r1(2)+r1(3)/2 0];
 mdl=feutil('Objectquad 1 1',mdl,node,1,1);
 il=feutilb('geomrbbeam1',mdl,[]); %[0,yc,0]);
 il(3)=(r1(3)^3*r1(1)+r1(4)^3*(r1(2)-r1(3)/2))/3; % J according to nastran reference
 il(7)=(r1(2)-r1(3))*r1(4)/il(6); % k1
 il(8)=r1(3)*r1(1)/il(6); % k2
 out=struct('sectionmdl',mdl,'il',il,...
     'name',Cam,'type','p_beam','unit','SI');

elseif comstr(Cam,'i') % - - - - - - - - - - - - - - - - - - - - - - - - - - I
 % 2004 nastran reference manual p90 nastran.pdf
 % Ansys : W1,W2 width of top and bottom flanges
 [r1,RS]=p_get_prop(comstr(Cam,'i','%g'),st(strcmpi('i',st(:,1)),:));
 
 if length(r1)<6; sdtw('_nb','I section needs 6 dimensions'); r1=[1 1 1 .2 .2 .2]; end
 node=[-r1(3)/2 r1(1)/2 0; r1(3)/2 r1(1)/2 0;
     r1(3)/2 r1(1)/2-r1(6) 0; -r1(3)/2 r1(1)/2-r1(6) 0];
 mdl=feutil('Objectquad 1 1',node,1,1);
 node=[-r1(2)/2 -r1(1)/2+r1(5) 0; r1(2)/2 -r1(1)/2+r1(5) 0;
     r1(2)/2 -r1(1)/2 0; -r1(2)/2 -r1(1)/2 0]; 
 mdl=feutil('Objectquad 1 1',mdl,node,1,1);
 node=[-r1(4)/2 r1(1)/2-r1(6) 0; r1(4)/2 r1(1)/2-r1(6) 0;
     r1(4)/2 -r1(1)/2+r1(6) 0; -r1(4)/2 -r1(1)/2+r1(6) 0];
 mdl=feutil('Objectquad 1 1',mdl,node,1,1);
 il=feutilb('geomrbbeam1',mdl,[]); % at centroid
 mdl=feutil('AddSetNodeID',mdl,'StressRP',[2 8 7 1]'); % Stress Recovery Points C D E F
 mdl2=mdl; % contourmdl;
 mdl2.Elt=feutil('Objectbeamline 1 2 4 10 12 6 8 7 5 11 9 3 1',mdl2.Node);
 mdl=stack_set(mdl,'info','contour',mdl2);
 il(3)=(r1(6)^3*r1(3)+r1(5)^3*r1(2)+r1(4)^3*(r1(1)-0.5*(r1(5)+r1(6))))/3; % J
 il(7)=r1(4)*(r1(1)-r1(5)-r1(6))/il(6); % k1
 il(8)=5*(r1(2)*r1(5)+r1(3)*r1(6))/(6*il(6)); % k2
 out=struct('sectionmdl',mdl,'il',il,...
     'name',Cam,'type','p_beam','unit','SI');

elseif comstr(Cam,'l') % - - - - - - - - - - - - - - - - - - - - - - - - - - L
 % 2004 nastran reference manual p90 nastran.pdf
 %r1=comstr(Cam,'l','%g');
 [r1,RS]=p_get_prop(comstr(Cam,'l','%g'),st(strcmpi('l',st(:,1)),:));
 if length(r1)<4; sdtw('_nb','L section needs 4 dimensions'); r1=[.6 1 .2 .1]; end
 node=[-r1(4)/2 r1(2)-r1(3)/2 0; r1(4)/2 r1(2)-r1(3)/2 0;
     r1(4)/2 -r1(3)/2 0; -r1(4)/2 -r1(3)/2 0];
 mdl=feutil('Objectquad 1 1',node,1,1);
 node=[r1(4)/2 r1(3)/2 0; r1(1)-r1(4)/2 r1(3)/2 0;
     r1(1)-r1(4)/2 -r1(3)/2 0; r1(4)/2 -r1(3)/2 0];
 mdl=feutil('Objectquad 1 1',mdl,node,1,1);
 il=feutilb('geomrbbeam1',mdl,[]); % at centroid
 il(3)=(r1(3)^3*(r1(1)-0.5*r1(4))+r1(4)^3*(r1(2)-0.5*r1(3)))/3; % J
 il(7)=(r1(2)-r1(3))*r1(4)/il(6); % k1
 il(8)=(r1(1)-r1(4))*r1(3)/il(6); % k2
 out=struct('sectionmdl',mdl,'il',il,...
     'name',Cam,'type','p_beam','unit','SI');

elseif comstr(Cam,'bar') % - - - - - - - - - - - - - - - - - - - - - - - - BAR
 % 2004 nastran reference manual p90 nastran.pdf
 %r1=comstr(Cam,'bar','%g');
 [r1,RS]=p_get_prop(comstr(Cam,'bar','%g'),st(strcmpi('bar',st(:,1)),:));
 if length(r1)<2; sdtw('_nb','BAR section needs 2 dimensions'); r1=[1 1]; end
 node=[-r1(1)/2  r1(2)/2 0;  r1(1)/2  r1(2)/2 0;
        r1(1)/2 -r1(2)/2 0; -r1(1)/2 -r1(2)/2 0];
 mdl=feutil('Objectquad 1 1',node,1,1);
 il=feutilb('geomrbbeam1',mdl,[]); % at centroid
 il(3)=r1(1)*r1(2)^3*(1/3-0.21*r1(2)/r1(1)*(1-r1(2)^4/(12*r1(1)^4))); %J, nastran ref
 il(7)=5/6; % k1
 il(8)=5/6; % k2
 out=struct('sectionmdl',mdl,'il',il,...
     'name',Cam,'type','p_beam','unit','SI');

elseif comstr(Cam,'box') % - - - - - - - - - - - - - - - - - - - - - - - - BOX
 % 2004 nastran reference manual p90 nastran.pdf
 %r1=comstr(Cam,'box','%g'); 
 [r1,RS]=p_get_prop(comstr(Cam,'box','%g'),st(strcmpi('box',st(:,1)),:));
 if length(r1)<4; sdtw('_nb','BOX section needs 4 dimensions'); r1=[1 1.2 .2 .1]; end
 node=[-r1(1)/2 r1(2)/2 0; r1(1)/2 r1(2)/2 0;
     r1(1)/2 r1(2)/2-r1(3) 0; -r1(1)/2 r1(2)/2-r1(3) 0];
 mdl=feutil('Objectquad 1 1',node,1,1);
 node=[r1(1)/2-r1(4) r1(2)/2-r1(3) 0; r1(1)/2 r1(2)/2-r1(3) 0;
     r1(1)/2  -r1(2)/2+r1(3) 0; r1(1)/2-r1(4) -r1(2)/2+r1(3) 0];
 mdl=feutil('Objectquad 1 1',mdl,node,1,1);
 node=[-r1(1)/2 r1(2)/2-r1(3) 0; -r1(1)/2+r1(4) r1(2)/2-r1(3) 0;
     -r1(1)/2+r1(4)  -r1(2)/2+r1(3) 0; -r1(1)/2 -r1(2)/2+r1(3) 0];
 mdl=feutil('Objectquad 1 1',mdl,node,1,1);
 node=[-r1(1)/2 -r1(2)/2+r1(3) 0; r1(1)/2 -r1(2)/2+r1(3) 0;
     r1(1)/2 -r1(2)/2 0; -r1(1)/2 -r1(2)/2 0];
 mdl=feutil('Objectquad 1 1',mdl,node,1,1);
 il=feutilb('geomrbbeam1',mdl,[]); % at centroid
 il(3)=(2*r1(4)*r1(3)*(r1(1)-r1(4))^2*(r1(2)-r1(3))^2)/(r1(1)*r1(4)+r1(2)*r1(3)-r1(4)^2-r1(3)^2); % J
 il(7)=(r1(2)-2*r1(3))*r1(4)*2/il(6); % k1
 il(8)=(r1(1)-2*r1(4))*r1(3)*2/il(6); % k2
 out=struct('sectionmdl',mdl,'il',il,...
     'name',Cam,'type','p_beam','unit','SI');
 
elseif comstr(Cam,'chan1'); % - - - - - - - - - - - - - - - - - - - - -  CHAN1
 % 2004 nastran reference manual p90 nastran.pdf
 %r1=comstr(Cam,'chan1','%g');
 [r1,RS]=p_get_prop(comstr(Cam,'chan1','%g'),st(strcmpi('chan1',st(:,1)),:));
 if length(r1)<4; sdtw('_nb','CHAN1 section needs 4 dimensions'); r1=[.5 .5 .8 1]; end
 node=[-r1(2)/2 r1(4)/2 0; r1(2)/2 r1(4)/2 0;
       r1(2)/2 -r1(4)/2 0; -r1(2)/2 -r1(4)/2 0];
 mdl=feutil('Objectquad 1 1',node,1,1);
 node=[r1(2)/2 r1(4)/2 0; r1(2)/2+r1(1) r1(4)/2 0;
       r1(2)/2+r1(1) r1(3)/2 0; r1(2)/2 r1(3)/2 0];
 mdl=feutil('Objectquad 1 1',mdl,node,1,1);
 node=[r1(2)/2 -r1(4)/2 0; r1(2)/2+r1(1) -r1(4)/2 0;
       r1(2)/2+r1(1) -r1(3)/2 0; r1(2)/2 -r1(3)/2 0];
 mdl=feutil('Objectquad 1 1',mdl,node,1,1);
 il=feutilb('geomrbbeam1',mdl,[]); %  at centroid
 il(3)=(2*(r1(1)+r1(2)/2)*(r1(4)/2-r1(3)/2)^3+(r1(4)/2+r1(3)/2)*r1(2)^3)/3;
 il(7)=r1(2)*r1(3)/il(6); % k1
 il(8)=2*0.5*(r1(4)-r1(3))*r1(1)/il(6); % k2
 out=struct('sectionmdl',mdl,'il',il,...
     'name',Cam,'type','p_beam','unit','SI');

elseif comstr(Cam,'chan2'); % - - - - - - - - - - - - - - - - - - - - -  CHAN2
 % 2004 nastran reference manual p90 nastran.pdf
 %r1=comstr(Cam,'chan2','%g');
 [r1,RS]=p_get_prop(comstr(Cam,'chan2','%g'),st(strcmpi('chan2',st(:,1)),:));
 if length(r1)<4; sdtw('_nb','CHAN2 section needs 4 dimensions'); r1=[.1 .1 1 1]; end
 node=[-r1(4)/2 r1(3)-r1(2)/2 0; -r1(4)/2+r1(1) r1(3)-r1(2)/2 0;
       -r1(4)/2+r1(1) -r1(2)/2 0; -r1(4)/2 -r1(2)/2 0];
 mdl=feutil('Objectquad 1 1',node,1,1);
 node=[-r1(4)/2+r1(1) r1(2)/2 0; r1(4)/2-r1(2) r1(2)/2 0;
       r1(4)/2-r1(2) -r1(2)/2 0; -r1(4)/2+r1(2) -r1(2)/2 0];
 mdl=feutil('Objectquad 1 1',mdl,node,1,1);
 node=[r1(4)/2-r1(1) r1(3)-r1(2)/2 0; r1(4)/2 r1(3)-r1(2)/2 0;
       r1(4)/2 -r1(2)/2 0; r1(4)/2-r1(1) -r1(2)/2 0];
 mdl=feutil('Objectquad 1 1',mdl,node,1,1);
 il=feutilb('geomrbbeam1',mdl,[]); % at centroid
 il(3)=(2*(r1(3)-r1(2)/2)*r1(1)^3+(r1(4)-r1(1))*r1(2)^3)/3; % J at centroid cf nastran reference
 il(7)=2*r1(1)*(r1(3)-r1(2))/il(6); % k1
 il(8)=r1(2)*(r1(4)-2*r1(1))/il(6); % k2
 out=struct('sectionmdl',mdl,'il',il,...
     'name',Cam,'type','p_beam','unit','SI');

elseif comstr(Cam,'chan'); % - - - - - - - - - - - - - - - - - - - - - - -CHAN
 % 2004 nastran reference manual p90 nastran.pdf
 %r1=comstr(Cam,'chan','%g');
 [r1,RS]=p_get_prop(comstr(Cam,'chan','%g'),st(strcmpi('chan',st(:,1)),:));
 if length(r1)<4; sdtw('_nb','CHAN section needs 4 dimensions'); r1=[.5 1 .1 .2]; end
 node=[-r1(3)/2 r1(2)/2-r1(4) 0; r1(3)/2 r1(2)/2-r1(4) 0;
       r1(3)/2 -r1(2)/2+r1(4) 0; -r1(3)/2 -r1(2)/2+r1(4) 0];
 mdl=feutil('Objectquad 1 1',node,1,1);
 node=[-r1(3)/2 r1(2)/2 0; r1(1)-r1(3)/2 r1(2)/2 0;
       r1(1)-r1(3)/2 r1(2)/2-r1(4) 0; -r1(3)/2 r1(2)/2-r1(4) 0];
 mdl=feutil('Objectquad 1 1',mdl,node,1,1);
 node=[-r1(3)/2 -r1(2)/2 0; r1(1)-r1(3)/2 -r1(2)/2 0;
       r1(1)-r1(3)/2 -r1(2)/2+r1(4) 0; -r1(3)/2 -r1(2)/2+r1(4) 0];
 mdl=feutil('Objectquad 1 1',mdl,node,1,1);
 il=feutilb('geomrbbeam1',mdl,[]); % at centroid
 il(3)=(2*(r1(1)-r1(3)/2)*r1(4)^3+(r1(2)-r1(4))*r1(3)^3)/3; % J at centroid cf nastran reference
 il(7)=r1(3)*(r1(2)-2*r1(4))/il(6); % k1
 il(8)=2*r1(4)*(r1(1)-r1(3))/il(6); % k2
 out=struct('sectionmdl',mdl,'il',il,...
     'name',Cam,'type','p_beam','unit','SI');
elseif comstr(Cam,'hat'); % - - - - - - - - - - - - - - - - - - - - - - -HAT

else;sdtw('%s unknown, using default beam property',CAM); out=out(1);
end
if isfield(RS,'NSM')&&RS.NSM; out.il(10)=RS.NSM;end
if isfield(RS,'lump')&&RS.lump; out.il(9)=RS.lump; end

%% --------------------------------------------------------------------
 
