function [o1,o2,o3,o4,o5]=fe_mat(varargin)

%FE_MAT general material/section property handling routine
%
%       Syntax:   out = fe_mat('unit')
%                 out = fe_mat('unitlabel',UnitSystemCode)
%                 typ=fe_mat('type','m_function',UnitCode,SubType)
%                 [m_function,UnitCode,SubType]=fe_mat('type',typ)
%                 [o1,o2,o3]=fe_mat(ElemP,ID,pl,il)
%                 out = fe_mat('convert si ba',pl);
%                 model=fe_mat('defaultil');
%                 model=fe_mat('defaultpl');
%                 rho=fe_mat('getMat MatId Rho',model);
%
%	This function is used by a number of element functions to support 
%	material property declarations. The input arguments are the
%	element function name or parent name (ElemP), ID =  [MatId ProId]
%	contains identification numbers for the current material/property,
%	PL and IL are the material/element property matrices.
%
%	Outputs depend on the element type and are detailed in the manual.
%
%	The general form for material property rows is
%	   [MatId Type Prop]
%	  with
%	   MatId a positive integer identifying a particular material property
%	   Type  identifying the type of material. Currently supported types
%	         1 standard isotropic
%	           no other type currently supported
%	   Prop  as many properties as needed
%
% See <a href="matlab: sdtweb _taglist fe_mat">TagList</a>
%	See sdtweb     pl, m_elastic, elem0, eltfun
%	See also help  m_elastic, fe_mk, beam1, bar1, tria3, ...

%       Etienne Balmes
%       Copyright (c) 2001-2025 by INRIA and SDTools
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       All Rights Reserved.

if comstr(varargin{1},'cvs')
 o1='$Revision: 1.235 $  $Date: 2025/10/28 18:11:22 $'; return;
end
%#ok<*NASGU,*ASGLU,*NOSEM>
if nargin==0; help fe_mat;return; end
% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
if ischar(varargin{1})

 [CAM,Cam]=comstr(varargin{1},1); carg=2;

if comstr(Cam,'p_'); Cam='type'; carg=1; end
if comstr(Cam,'m_'); Cam='type'; carg=1; end

if comstr(Cam,'steel');    o1 =m_elastic('database','steel');
elseif comstr(Cam,'alum'); o1 =m_elastic('database','aluminum');

%% #Get ----------------------------------------------------------------------
elseif comstr(Cam,'get');  [CAM,Cam]=comstr(CAM,4);
    
 %% #GetPl #GetIl - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 % Create the pl field no matter initial format
 % pl=fe_mat('getpl',model)
 if comstr(Cam,'pl')||comstr(Cam,'il'); 
  if comstr(Cam,'pl'); st1='pl'; st2='MatId'; st3='mat';
  else; st1='il'; st2='ProId'; st3='pro';
  end
  [CAM,Cam]=comstr(CAM,3); 
  [CAM,Cam,RO.lin]=comstr('-lin',[-25 3],CAM,Cam);
  [CAM,Cam,RO.used]=comstr('-used',[-25 3],CAM,Cam);
  opt=comstr(CAM,-1);
  model=varargin{carg};carg=carg+1;
  val=stack_get(model,st3);
  pl=[];  if isfield(model,st1); pl=model.(st1); end
  if isobject(pl); pl=sdth.GetData(pl); end
  if RO.used&&isfield(model,'Elt')&&~isempty(model.Elt)
   mpid=feutil(st2,model.Elt); opt=unique(mpid(isfinite(model.Elt(:,1)))); opt(opt==0)=[];
  end
  if ~isempty(pl)&&(RO.used||~isempty(opt));pl(~ismember(pl(:,1),opt),:)=[];end
  ind=1:size(val,1);
  if strcmpi(CAM,';'); field_interp(CAM);end
  for j1=1:size(val,1)
   r1=val{j1,3};
   if isfield(r1,st1)&&(~RO.lin||(RO.lin&&~isfield(r1,'NLdata')))
    r0=r1.(st1); if isa(r0,'v_handle'); r0=r0.GetData; end
    if length(r0)<3;ind(j1)=0;continue
    else
     [r1,model,i3]=field_interp(r1,model);
     r2=r1.(st1)(:)'; if ~isempty(opt)&&~any(r2(1)==opt);ind(j1)=0;continue;end
     r0(i3==1)=r2(i3==1);% Was interpolated
    end
    if isempty(pl)||isempty(r2); i2=[];else;i2=find(pl(:,1)==r2(1));end 
    if isscalar(r0)||isempty(i2)
    elseif length(r2)<=size(pl,2)
      if ~isequal(r0(:)',pl(i2,1:length(r0)))
       sdtw('_nb','Using ''%s'',''%s'' .%s field (differ from model.%s)', ...
           st3,val{j1,2},st1,st1)
      end
    elseif any(r2(size(pl,2)+1:end))|| ...
            ~isequal(r2(1:min(size(pl,2),length(r2))), ...
            pl(i2,1:min(size(pl,2),length(r2))))
       sdtw('_nb','Using ''%s'',''%s'' .%s field (differ from model.%s)', ...
          st3,val{j1,2},st1,st1)
    end
    if isempty(i2); pl(end+1,1:length(r2))=r2; %#ok<AGROW>
    elseif isscalar(r2) % .pl entry is pre-emptive
    elseif length(i2)>1
     r3=unique(pl(i2,:),'rows');
     if size(r3,1)==1; sdtw('Coalesced duplicated %s %s %i rows',st1,st2,r2(1))
      pl(i2(1),1:length(r2))=r2;pl(i2(2:end,:),:)=[];
     else; error('Duplicated %s %i rows in model.%s',st2,r2(1),st1);
     end
    else % replace the pl value by the stack value
     if size(pl,2)>length(r2); pl(i2,:)=0; end % /!\ reset the full .pl value
     pl(i2,1:length(r2))=r2;
    end
   end
  end
  o1=pl;o2=[];
  if any(ind);o2=val(ind~=0,:);if isscalar(opt);o2=o2{1,3};end;end
  if comstr(Cam,'of_mk'); error('You should use GetPlIlOfMk');end
  % for j1=1:size(o1,1)  r1=fe_mat('of_mk',o1(j1,:)); 
  %  o1(j1,:)=0;o1(j1,1:length(r1))=r1;
  % end
  %end
  
 %% #GetPro #GetMat : Get the given prop parameter - - - - - - - - - - - - - - 
 % rho=fe_mat('GetMat 100 rho',model)
 elseif comstr(Cam,'mat') || comstr(Cam,'pro');
  model=varargin{carg}; carg=carg+1;
  if comstr(Cam,'mat'); 
   RunOpt.st='pl'; RunOpt.st0='mat';
   [CAM,Cam,RunOpt.i1]=comstr('mat',[-25 31],CAM,Cam);   
  else
   RunOpt.st='il'; RunOpt.st0='pro';
   [CAM,Cam,RunOpt.i1]=comstr('pro',[-25 31],CAM,Cam);   
  end
  if isempty(RunOpt.i1); error('You must give %id',RunOpt.st0); end
  
  RunOpt.param=strtrim(Cam);
  if isempty(RunOpt.param);
  else
   % check if multiple parameters are asked
   st=textscan(RunOpt.param,'%s','delimiter',' ,','multipledelimsasone',true);
   st=st{1};
   if length(st)>1 % if then recursive call to the single param implementation
    o1=zeros(size(st));
    for j1=1:length(o1);
     o1(j1)=fe_mat(sprintf('Get%s %i %s',RunOpt.st0,RunOpt.i1,st{j1}),model);
    end
    return
   end
  end
  
  if nargout>1;[plil,o2]=fe_mat(sprintf('Get%s',RunOpt.st),model);
  else; plil=fe_mat(sprintf('Get%s',RunOpt.st),model);o2=[];
  end
  i1=plil(:,1)~=RunOpt.i1; plil(i1,:)=[];
  if ~isempty(o2);o2(~cellfun(@(x)ismember(x.(RunOpt.st)(1),RunOpt.i1),o2(:,3)),:)=[];end
  if isempty(plil); error('%sId %i not found',RunOpt.st0,RunOpt.i1); end
  o1=plil;
  if ~isempty(RunOpt.param)
   [st0,UnitCode,SubType]=fe_mat('type',plil(2));
   if comstr(RunOpt.st,'il'); st0(1)='p'; end
   [r1,st1]=feval(st0,'propertyunittypecell',SubType);
   if size(r1,1)<size(o1,2)&&~isempty(st1) % generate extensions
    r2=double(strcmpi(r1{end,1},sprintf(st1{end,1},1)));
    st2=cellfun(@(x)sprintfc(x,r2+1:(size(o1,2)-size(r1,1))),st1(:,1),'uni',0);
    st2=reshape(cat(1,st2{:}),[],1);
    r1(end+1:size(o1,2),1)=st2(1:size(o1,2)-size(r1,1));
   end
   %i1=find(strncmpi(RunOpt.param,r1(:,1),length(RunOpt.param)));
   i1=find(ismember(lower(r1(:,1)),lower(RunOpt.param)));
   if ~isempty(i1)
   elseif strcmp(RunOpt.param,'struct')||strcmp(RunOpt.param,'-struct')
    if size(plil,2)>size(r1,1);plil(size(r1,1)+1:end)=[];end
    r1=[r1(1:length(plil),1)';num2cell(plil)];r1=struct(r1{:});
    o1=r1;return;
   else
    error('Parameter %s is not a %s %i parameter',RunOpt.param,st0,SubType)
   end
   if i1>size(plil,2);o1=zeros(size(i1));
   else; o1=plil(1,i1);
   end
  end
  
 elseif comstr(Cam,'pos');  [CAM,Cam]=comstr(CAM,4); 
 %% #GetPos fe_mat('getpos(p)',type,name)
  r1=varargin{carg};carg=carg+1; 
  st=varargin{carg};carg=carg+1; 
  if size(r1,2)>1;r1=r1(1,2);end
  [st1,i1,i2]=fe_mat(['type' Cam],r1);
  if isempty(st1); o1=[];
  else
   if ~exist(st1,'file');o1=[];return;
   else
    st2=feval(st1,'propertyunittype cell',i2);
   end
   o1=find(strcmpi(st2(:,1),strrep(st,';','')));
   if length(o1)~=1&&(isempty(CAM)||CAM(end)~=';');
       error('Did not find ''%s'' in %s subtype %i',st,st1,i2);
   end
  end
 end
 
 %% #Set ---------------------------------------------------------------------
 elseif comstr(Cam,'set');  [CAM,Cam]=comstr(CAM,4);
 
  %% #SetPro used to define element properties - - - - - - - - - - - - - - - - - - - - - 
  if comstr(Cam,'pro');

   model=varargin{carg};carg=carg+1;
   [CAM,Cam,doCat]=comstr('cat',[-25 3],CAM,Cam); % cat fields in struct entries
   [CAM,Cam,i1]=comstr('pro',[-25 1],CAM,Cam);
   if isempty(i1)&&carg<=nargin&&isnumeric(varargin{carg}); 
     %model=feutil('setpro',model,il);
     il=varargin{carg};carg=carg+1;
     if ~isfield(model,'il')||isempty(model.il);model.il=[];     end
     for j1=1:size(il,1)
       if isempty(model.il);i2=1;else;i2=find(model.il(:,1)==il(j1,1)); end
       if isempty(i2);i2=size(model.il,1)+1;
       elseif length(i2)>1;warning('Replicated ProId %i',model.il(i2(1)));
       end
       model.il(i2(1),1:size(il,2))=il(j1,:);
     end
     if carg>nargin;o1=model; return;else; i1=il(:,1);end
   end
   if isempty(i1)&&carg<=nargin&& ...
           (ischar(varargin{carg})||isfield(varargin{carg},'il'))
     % new calls with full assign in arguments
      %fe_mat('setpro',cf.mdl,'m_elastic(dbval1 aluminum)')
      mlist=varargin(carg:end); if isscalar(mlist)&&iscell(mlist{1}); mlist=mlist{1}; end
      if ~isfield(model,'nmap'); model.nmap=vhandle.nmap([]); end
      if isKey(model.nmap,'Map:SetName')
       matM=[];nameM=model.nmap('Map:SetName');matcM=model.nmap('Map:SetColor');
      else
       nameM=[];matcM=model.nmap('Map:MatColor'); matM=model.nmap('Map:MatName');
      end
      for j1=1:length(mlist)
       matj=mlist{j1};
       if ischar(matj)
        st=matj; matj=sdtm.urnCb(matj);
        if iscell(matj); 
          matj=feval(matj{:}); 
        else; matj=eval(matj);
        end
       end
       if isnumeric(matj);
           matj={'pro','?',struct('il',matj)};
           st1=regexprep(st,'.*[Pp]ro[Dd]b\.([^)]*)\)','$1');
           if ~isempty(st)&&~isequal(st1,st);matj{3}.name=st1;end
       elseif isfield(matj,'il');
           matj={'pro','?',matj};
       end
       if ~isfield(matj{3},'name')
       elseif isequal(matM,[])
        nameM(sprintf('Mat:%i',matj{3}.il(1)))=matj{3}.name;
       else
        matM(matj{3}.name)=matj{3}.il(1); %#ok<AGROW>
       end
       if ~isfield(matj{3},'color')
       elseif isequal(matM,[]);
           matcM(sprintf('Mat:%i',matj{3}.il(1)))=matj{3}.color;
       else;  matcM(matj{3}.il(1))=matj{3}.color;
       end
       if isequal(matj{2},'?')||isempty(setdiff(fieldnames(matj{3}),{'il','name','color'}))
        if ~isfield(model,'il')||isempty(model.il);model.il=[];i2=1;
        else; i2=find(model.il(:,1)==matj{3}.il(1));
        end
        if isempty(i2);i2=size(model.il,1)+1;end
        model.il(i2,1:size(matj{3}.il,2))=matj{3}.il;
       else
        model=stack_set(model,matj);
       end
      end % mlist
      o1=model; return;
   end

   [pro,il]=matgui('getstackil',model,i1);% Get property but no fill of pro.il
   if ~isempty(pro)
   elseif any(strncmpi(varargin(cellfun(@ischar,varargin)),'nldata',6))
     %% allow NLdata for null properties
     pro={'pro',sprintf('nl%i',i1),struct('il',i1,'type','p_nl')};
   else
       error('ProId %i is not defined',i1);
   end
   if ~isempty(CAM) % SDT Set some values in il
     try;eval('[model,pro]=matgui([''setpro'',CAM],model,pro);');
     catch; sdtw('(%i).il setting failed',pro{3}.il(1));
     end
   end
   if size(pro,1)>1;error('Problem with multiple match');end
   r1=pro{3};
   if carg<=nargin&&isnumeric(varargin{carg})&&varargin{carg}(1)==i1
     il=varargin{carg};carg=carg+1;
     if length(r1.il)>1; r1.il=il; %feutil('setpro 1000',mdl,[1000 ...])
     else; 
       [r1.type,un1,un2]=fe_mat('typep',il(2));
       if ~isfield(model,'il')||isempty(model.il);model.il=[];i2=0;
       else; i2=find(model.il(:)==il(1));
       end
       if isempty(i2)||i2==0;i2=size(model.il,1)+1;end
       model.il(i2,1:length(il))=il;model.il(i2,length(il)+1:end)=0;
     end
   end

   while carg<nargin
       st=varargin{carg};r2=varargin{carg+1};carg=carg+2;
       if strcmpi(st,'map')||strcmpi(st,'InfoAtNode'); st1='MAP';
       elseif strcmpi(st,'gstate');st1='gstate'; 
       elseif strcmpi(st,'name');% Possibly change stack name
           [r3,i1]=stack_get(model,pro{1:2});r1.name=r2;
           if ~isempty(i1);model.Stack{i1,2}=r2;end
           pro{2}=r2; continue; 
       elseif strcmpi(st,'nldataedit');% Replace fields in nldata
           if ~isfield(r1,'NLdata');r1.NLdata=r2;
           else; r2=sdth.sfield('MergeI -empty',r1.NLdata,r2);end
           st1='NLdata';
       elseif strcmpi(st,'nldata');st1='NLdata';
       elseif strcmpi(st,'param');st1='param';
        if ~isfield(r1,'param'); r1.param=struct;
        elseif doCat %&&isfield(r1,'param')
         r2=sdth.sfield('AddMissing',r1.param,r2);
        end
       elseif strcmpi(st,'StressOut');
         st1='StressOut';
         if ischar(r2)&&strcmpi(r2,'default')
          r2=struct('type',r1.type);
         end
       else; error('SetOut %s not implemented',st);
       end
       if isempty(r2);r1=feutil('rmfield',r1,st1);else;r1.(st1)=r2;end
   end
   pro{3}=r1; model=stack_set(model,pro);
   o1=model;
   
  %% #SetMat used to define material properties - - - - - - - - - - - - - - - - 
  % ToDo accept struct input :  st=[fieldnames(RO) struct2cell(RO)]';
  % st(:,strcmpi(st(1,:),'matid'))=[];st=sprintf('%s=%.15g ',st{:});
  % model=feutil(sprintf('setmat %i %s',RO.MatId,st),model);    
  elseif comstr(Cam,'mat');
   [CAM,Cam,doCat]=comstr('cat',[-25 3],CAM,Cam);
   [CAM,Cam,i1]=comstr('mat',[-25 1],CAM,Cam);
   model=varargin{carg};carg=carg+1;
   if isempty(i1)&&carg<=nargin&&isnumeric(varargin{carg});
    %model=feutil('setmat',model,pl);
    pl=varargin{carg};carg=carg+1;
    if ~isfield(model,'pl')||isempty(model.pl);model.pl=[];     end
    for j1=1:size(pl,1)
     if isempty(model.pl);i2=1;else;i2=find(model.pl(:,1)==pl(j1,1)); end
     if isempty(i2);i2=size(model.pl,1)+1;end
     model.pl(i2,1:size(pl,2))=pl(j1,:);
    end
    o1=model; return;

   else
    if isempty(i1) % accept input with mat Stack Name
     % model=feutil('setmat"Acier" Rho=1.5e-09',model);
     [CAM,Cam]=comstr(varargin{1},4); % restart parsing, remove 'Set'
     [CAM,Cam,i1]=comstr('mat',[-25 4],CAM,Cam);

     if isempty(i1) % new calls with full assign in arguments
      %fe_mat('setmat',cf.mdl,'m_elastic(dbval1 aluminum)')
      mlist=varargin(carg:end); if isscalar(mlist)&&iscell(mlist{1}); mlist=mlist{1}; end
      if ~isfield(model,'nmap'); model.nmap=vhandle.nmap([]); end
      if isKey(model.nmap,'Map:SetName')||~isKey(model.nmap,'Map:MatName')
       matM=[];nameM=model.nmap('Map:SetName');
      else
       nameM=[];matcM=model.nmap('Map:MatColor'); matM=model.nmap('Map:MatName');
      end
      for j1=1:length(mlist)
       matj=mlist{j1};
       if ischar(matj)
        matj=sdtm.urnCb(matj);
        if iscell(matj); 
          try;
             r1=feval(matj{:}); 
             if iscell(r1)&&isequal(size(r1),[1 3])&&isfield(r1{3},'pl');matj=r1;
             elseif isfield(r1,'pl')||isnumeric(r1);matj=r1;
             else
                     error('Old call');
             end
          catch; 
            [un1,matj]=feval(matj{:});% xxx out2 is not ideal
          end
        else; matj=eval(matj);
        end
       end
       if isequal(matM,[])
        nameM(sprintf('Mat:%i',matj{3}.pl(1)))=matj{3}.name;
       else
        matM(matj{3}.name)=matj{3}.pl(1); %#ok<AGROW>
       end
       if isfield(matj{3},'color')
        matcM(matj{3}.pl(1))=matj{3}.color;
       end
       if isequal(matj{2},'?')||isempty(setdiff(fieldnames(matj{3}),{'pl','name','color'}))
        if ~isfield(model,'pl')||isempty(model.pl);model.pl=[];i2=1;
        else; i2=find(model.pl(:,1)==matj{3}.pl(1));
        end
        if isempty(i2);i2=size(model.pl,1)+1;end
        model.pl(i2,1:size(matj{3}.pl,2))=matj{3}.pl;
       else
        model=stack_set(model,matj);
       end
      end % mlist

      o1=model; return;
     end
    end% i1/pro recovery

    % we got an identifier (stack name/matid)
    % base case, other calls already returned
    [pro,il]=matgui('getstackpl',model,i1);
    % xxx setpro forwarding maybe not the best thing to authorize
    if isempty(pro);error('MatId %s is not defined',num2str(i1));end
    if ~isempty(CAM) % SDT Set some values in il
     try;eval('[model,pro]=matgui([''setpro'',CAM],model,pro);');
     catch; sdtw('(%i).il setting failed',pro{3}.il(1));
     end
    end
    r1=pro{3};
    % additionnal fields
    while carg<nargin
     st=varargin{carg};r2=varargin{carg+1};carg=carg+2;
     if strcmpi(st,'param');st1='param';
      if ~isfield(r1,'param'); r1.param=struct;
      elseif doCat;r2=sdth.sfield('AddMissing',r1.param,r2);
      end
     end
     if isempty(r2);r1=feutil('rmfield',r1,st1);else;r1.(st1)=r2;end
    end
    pro{3}=r1;

    model=stack_set(model,pro);
    o1=model;

   end % shortcut for pl vector assignment or not

   %%
  else;error('Set%s not a valid command',CAM);
  end

%% #convert : unit conversion support ----------------------------------------
% convert from a specified unit system to another
% pl1=fe_mat('convert -il',pl,PropertyUnitTypeVector);
% coef=fe_mat('convert SIMM',ucode)
% coef=fe_mat('convert SIMM',uname); r1=fe_mat('convertsimm','struct')
elseif comstr(Cam,'convert');  [CAM,Cam]=comstr(CAM,8);

  RO.silent=~isempty(CAM)&&CAM(end)==';';
  if isempty(Cam)&&nargin==1
   r1=fe_mat('unitsystem');
   st=cellstr(num2str((1:length(r1))'));
   st(:,2)={r1.name};st=st';
   fprintf('%s  %s\n',st{:})
   return;
  end
  [r1,lab]=fe_mat('unitsystems');pl=[]; RO.Des=[];
  if carg>nargin
  elseif ischar(varargin{carg}); % fe_mat('convertMMSI','charge')
    RO.Des=varargin{carg};carg=carg+1;
  elseif isfield(varargin{carg},'Des');RO=varargin{carg};carg=carg+1;
      if ~isfield(RO,'coef');RO.coef=[];end
      if ~isfield(RO,'lab');RO.lab={};end
  else;pl=varargin{carg};carg=carg+1; 
  end
  if isfield(pl,'Elt') 
   %% Clean up all units of a model
   model=pl;
   if comstr(Cam,'us') % From us to given unit
     model.unit=upper(CAM(3:4));
   end
   if isfield(model,'unit')
     if ~isfield(model,'pl')||isempty(model.pl)
     else
      for j1=1:size(model.pl,1)
       [m_function,UnitCode,SubType]=fe_mat('type',model.pl(j1,2)); 
       if UnitCode==9 % US do not convert
          model.pl(j1,2)=fe_mat(m_function,model.unit,SubType);
       else; pl=fe_mat(['convert',model.unit],model.pl(j1,:));
           model.pl(j1,1:length(pl))=pl;
       end
      end
     end
     if ~isfield(model,'il')||isempty(model.il)
     else
      for j1=1:size(model.il,1)
       [m_function,UnitCode,SubType]=fe_mat('typep',model.il(j1,2)); 
       if UnitCode==9 % US do not convert
          model.il(j1,2)=fe_mat(m_function,model.unit,SubType);
       else; il=fe_mat(['convert',model.unit],model.il(j1,:));
           model.il(j1,:)=0; model.il(j1,1:length(il))=il;
       end
      end
     end
   else; error('Not yet implemented')
   end
   o1=model;return;
  
  elseif carg<=nargin % PropertyUnitType given
    ind=varargin{carg};carg=carg+1; 
    [m_function,UnitCode,SubType]=fe_mat('type',pl(2)); 
  elseif isempty(pl)  % no argument lists possibilities
   ind=1:size(lab,1);%ind(10)=0;ind=find(ind);
   if any(strcmpi(RO.Des,{'struct','ulab'}));
   elseif ~isempty(RO.Des);
       % fe_mat('convertMMSI','charge')
       if ~iscell(RO.Des);RO.Des={RO.Des};end
       ind=zeros(1,length(RO.Des));
       for j1=1:length(RO.Des)
        ind(j1)=find(strncmpi(lab(:,9),RO.Des{j1},length(RO.Des{j1})));
       end
   end
   UnitCode=0;
  elseif size(pl,1)>1  % Multiple conversion : loop and return
    for j1=1:size(pl,1)
      r1=fe_mat(['convert' CAM],pl(j1,:));
      pl(j1,1:length(r1))=r1;
    end
    o1=pl;
    return; 
  elseif isscalar(pl) % Convert a given unit
   ind=[fix(pl) round(rem(pl,1)*1000)];
   ind(ind==0)=[]; RO.Need='coef';
  else % standard call
   [m_function,UnitCode,SubType]=fe_mat('type',pl(2)); 
   [CAM,Cam,i1]=comstr('il',[-25 3],CAM,Cam);
   if ~exist(m_function,'file')||i1; m_function(1)='p';end
   eval(sprintf('[ind,in2]= %s(''PropertyUnitType'',%i);',m_function,SubType));
   if ~isempty(in2);%#ok<NODEF>
    in2=horzcat(in2{:,2}); i2=find(~isfinite(ind)); %#ok<NODEF>
    if ~isempty(i2) % (case of p_beam subtype 3) sdtweb p_beam('propteryunittype')
     ind=[ind(1:i2-1) repmat(in2,1,ceil((size(pl,2)-length(ind)+1)/length(in2))) ind(i2+1:end)];
    else % standard: additional parameters are in the end
     ind=[ind repmat(in2,1,ceil((size(pl,2)-length(ind))/length(in2)))];
    end
   elseif any(~isfinite(ind));
    error('Please report to SDTools. [Propertyunittype %s, %i].',m_function,SubType);
   end
  
  end

  if length(pl)>1 % get unit types
    [st,i1,isub]=fe_mat('type',pl(2)); 
    if length(pl)<length(ind); pl(length(ind))=0; end
  end

  % get initial and final unit systems
  if strncmpi(Cam,'US',2)
   st=CAM; [CAM,Cam]=comstr(CAM,3);
   if strncmpi(Cam,'US',2)
    o1=repmat({1},size(lab,1),1);o1(:,end+1)=lab(:,9);
    if strcmpi(RO.Des,'struct');
        o1=o1(:,[2 1])';o1(1,:)=cellfun(@(x)comstr(x,-36),o1(1,:),'uni',0)';
        o1=struct(o1{:});
    end
    return;
   end
   error('%s bad unit system. Can''t convert',st);
  end
  if length(CAM)==2||(length(CAM)==3&&CAM(3)==';'); i1=UnitCode; % input code
  else;i1=find(strncmpi(Cam(1:2),{r1.name},2));[CAM,Cam]=comstr(CAM,3);
  end
  RO.inUnit=i1; 
  i3=find(strncmpi(Cam(1:2),{r1.name},2)); % output code
  if isempty(i3) || isempty(i1); 
      error('bad unit system. Can''t convert');
  end
  if exist('UnitCode','var') && UnitCode && i1~=UnitCode 
   %warning('There is a mismatch in the initial unit system');
  end

  if i1==0;i1=i3;RO.inUnit=i3;
  elseif i1==9; 
      if CAM(end)~=';';
       sdtw('_nb','Id=%i uses US, assuming unit=%i %s',pl(1),i3,r1(i3).name(1:2));
      end
      i1=i3; RO.inUnit=i1;% input is US 
  end
  RO.outUnit=i3;
  % length, force, temp, temp-offset, time
  r2=r1(RO.outUnit).data(3:7)./r1(RO.inUnit).data(3:7);           % basic unit conversion
  r3=reshape([lab{:,end}],length(r2),size(lab,1))';% rows in basic units

  if isempty(ind)
  elseif any(~isfinite(r2))&&any(ind>0)
    error('%s%s Conversion not defined',r1(RO.inUnit).name(1:2),r1(RO.outUnit).name(1:2));
  else
    ind(ind<=0)=10; % code 10 for no unit 
    in3=rem(ind,1)*1000; % denominator
    in2=round(in3); in3=round(rem(in3,1)*1000);
    
    r4=zeros(length(ind),size(r3,2));
    if any(in2) % compute denominator coeff for each unit
     r4(in2~=0,:)=r3(in2(in2~=0),:);% denominator units
     r4=prod((r2(ones(length(in2),1),:)).^r4,2)'; 
     r4(in2==0)=1; % make sure no unit change for no denominator
    else; r4=ones(1,length(ind));
    end
    if 1==2 % any(in3) % Second denominator
     r5=zeros(length(ind),size(r3,2));
     r5(in3~=0,:)=r3(in3(in3~=0),:);% denominator units
     r5=prod((r2(ones(length(in3),1),:)).^r5,2)'; 
     r5(in3==0)=1; % make sure no unit change for no denominator
     r4=r4.*r5;
    end
    r3=r3(fix(ind),:); % numerator units 
    r3=prod((r2(ones(length(ind),1),:)).^r3,2)';
    r3=r3./r4;

    if length(pl)>1; % update values & type Identifier
     o1(1:length(ind))=pl(1:length(ind)).*r3;
     o1(2)=fe_mat('type',m_function,min(RO.outUnit,9),isub);% if >9 use US
    elseif strcmpi(RO.Des,'ulab');o1=lab(:,[9 RO.inUnit])';
        o1(1,:)=strrep(o1(1,:),' ','_');o1=struct(o1{:});
    elseif strcmpi(RO.Des,'struct') % fe_mat('convertSIMM','struct')
     o1=[cellfun(@(x)comstr(x,-36),lab(:,9),'uni',0)';num2cell(r3)];
     o1=struct(o1{:});
     if nargout>1; o2=lab(:,1);end
     if nargout==2
      o2=lab(:,[9 RO.inUnit RO.outUnit]);
     end
    else% display the conversions
     pl=ones(size(ind));pl(2,:)=r3;pl(3,:)=ind(:)';
     o1=num2cell(pl');o1(:,4:5)=lab(ind,[RO.inUnit RO.outUnit]);
     if isfield(RO,'coef')
       RO.coef=vertcat(o1{:,2})./vertcat(o1{:,2});
       RO.lab=o1(:,5);
       o1=RO;
     elseif nargout==0; o1=o1(:,[3 1 4 2 5]);o1(:,6)=lab(ind,9);
     elseif isfield(RO,'Need')&&strcmpi(RO.Need,'coef')||(~isempty(strfind(Cam,'coef')))
       r1=o1{1,2}/o1{1,1}; 
       if size(o1,1)>1; r1=r1/(o1{2,2}/o1{2,1});end
       o1=r1; % coefficient for unit change 
     else;o1=o1(:,[1 4 2 5]);
     end
   end
  end

%% #Unit (s) handling --------------------------------------------------------
elseif comstr(Cam,'unit');  [CAM,Cam]=comstr(CAM,5);

if comstr(Cam,'label') 
%% #UnitLabel fe_mat('unitLabel','SI') - - - - - - - - - - - - - - - - - -

o1={'Pa','lbf/ft^2','kgf/m^2','pdl/ft^2','milli-N/mm^2',...
                         'centi-N/cm^2','lbf/in^2(psi)','kgf/m^2','pressure','N/mm^2','GPa'};
o1(2,:)={'N','lbf','kgf','pdl','milli-N','centi-N','lbf','kgf','force','N','milli-N'};
o1(3,:)={'kg/m^3','lbf-s^2/ft^4','kgf-s^2/m^4','lbm/ft^3',...
          'kg/mm^3','kg/cm^3','lbf-s^2/in^4','kgf-s^2/m^4','density','t/mm^3','t/mum^3'};
o1(4,:)={'m','ft','m','ft','mm','cm','in','mm','length','mm','mum'};
o1(5,:)={'m/s','ft/s','m/s','ft/s','mm/s','cm/s','in/s','mm/s','speed','mm/s', ...
    'mum/s'};
o1(6,:)={'m/s^2','ft/s^2','m/s^2','ft/s^2','mm/s^2','cm/s^2','in/s^2','mm/s^2', ...
    'acceleration','mm/s^2','mum/s^2'};
o1(7,:)={'^oC','^oF','^oC','^oF','^oC','^oC','^oF','^oC','temperature','^oC','^oC'};
o1(8,:)={'/^oC','/oF','/oC','/oF','/oC','/oC','/oF','/oC','thermal coef','/oC','/oC'};
o1(9,:)={'kg','lbf-sec^2/ft','kgf-sec^2/m','lbm','kg','kg',...
                 'lbf-sec^2/in','kgf-sec^2/mm','mass','t','t'};
o1(10,:)={' ',' ',' ',' ',' ',' ',' ',' ','no unit',' ',' '};
o1(11,:)={'m^4','ft^4','m^4','ft^4','mm^4','cm^4','in^4','mm^4', ...
     'inertia','mm^4','mum^4'};
o1(12,:)={'m^2','ft^2','m^2','ft^2','mm^2','cm^2','in^2','mm^2', ...
    'area','mm^2','mum^2'};
o1(13,:)={'N/m','lbf/ft','kgf/m','pdl/ft','milli-N/mm',...
          'centi-N/cm','lbf/in','kgf/mm','stiffness','N/mm','kN/mum'};
o1(14,:)={'C','C','C','C','C','C','C','C','charge','C','C'};
o1(15,:)={'m^3','ft^3','m^3','ft^3','mm^3','cm^3','in^3','mm^3', ...
    'volume','mm^3','mum^3'};
o1(16,:)={'W/m/K','lbf/s/^oF','kgf/s/^oC','pdl/s/^oF','milli-N/s/^oC','centi-N/s/^oC', ...
    'lbf/s/^oF','kgf/s/^oC','heat conductivity','W/m/K',' '};
o1(17,:)={'J/kg/K','ft^2/s^2/^oF','J/kg/K','ft^2/s^2/^oF','mm^2/s^2/^oC', ...
    'cm^2/s^2/^oC','in^2/s^2/^oF','mm^2/s^2/^oC','heat capacity','mm^2/s^2/^oC',' '};
o1(18,:)={'F/m','permitivity','permitivity','permitivity','permitivity', ...
    'permitivity','permitivity','permitivity','permitivity','permitivity','permitivity'};
o1(19,:)={'V','tension','tension','tension','mu-V', ...
    'tension','tension','tension','tension','tension','tension'};
o1(20,:)={'s','s','s','s','s','s','s','s','time','s','s'};
o1(21,:)={'W','?','?','?','?','?','?','?','power','milli-W','milli-W'};
o1(22,:)={'J','?','?','?','?','?','?','?','work','milli-J','milli-W'};
 
%o1(:,10)={'pressure','force','density','length','speed','acceleration',...
%         'temperature','thermal coeff','mass','no unit','inertia','area',...
%         'stiffness','charge','volume','heat conductivity','heat capacity', ...
%         'permittivity','tension'};

% length, force, temp, temp-offset, time
o1(:,end+1)={[-2 1 0 0 0];[0 1 0 0 0];[-4 1 0 0 2];[1 0 0 0 0];
          [1 0 0 0 -1];[1 0 0 0 -2]; % m/s ms-2
          [0 0 1 0 0];[0 0 -1 0 0];[-1 1 0 0 2];[0 0 0 0 0];... % ... none
          [4 0 0 0  0];[2 0 0 0 0];[-1 1 0 0 0];[0 0 0 0 0]; ... % m4,m2,N/m,C
          [3 0 0 0 0];...
          [0 1 -1 0 -1];[2 0 -1 0 -2];[-2 -1 0 0 0] % perm
          [1 1 0 0 0];[0 0 0 0 1]% volt
          [1 1 0 0 -1];[1 1 0 0 0]; % power,work
          };

% System of Units        Length  Time      Mass          Force       Temp(R:A)
%----------------------------------------------------------------------------- 
% SI-Meter (newton)      meter   sec   kilogram(kg)       newton(N)  deg C:K
% BG-Foot (pound f)      foot    sec   lbf-sec**2/ft     pound(lbf)  deg F:R
% MG-Meter (kilogram f)  meter   sec   kgf-sec**2/m   kilogram(kgf)  deg C:K
% BA-Foot (poundal)      foot    sec   pound(lbm)      poundal(pdl)  deg F:R
% MM-mm (milli-newton)   mm      sec   kilogram(kg)    milli-newton  deg C:K
% CM-cm (centi-newton)   cm      sec   kilogram(kg)    centi-newton  deg C:K
% IN-Inch (pound f)      inch    sec   lbf-sec**2/in     pound(lbf)  deg F:R
% GM-mm (kilogram f)     mm      sec   kgf-sec**2/mm  kilogram(kgf)  deg C:K
% US-User defined        ____    sec   _____________  _____________  _____
% TM-mm (ton- mm)        mm      sec   ton               N           deg C:K
% MU-mum (kN-mum)        mum     sec   ton              kN           deg C:K

  carg=2;
  if carg<=nargin
   st=varargin{carg};carg=carg+1;
   if isfield(st,'unit'); st=st.unit;end
   if isnumeric(st)&&st>0&&st<10;    i1=st;
   elseif ~ischar(st);               error('unit must be a string');
   else;i1=find(ismember({'SI','BG','MG','BA','MM','CM','IN','GM','US','TM','MU'},st));
   end
   if ~isempty(i1);  o1=o1(:,i1); else;o1=o1(:,end-1); end
   if carg<=nargin % Allow InitForming with Num.Den(3 digit den)
     i1=varargin{carg};carg=carg+1;if iscell(i1);i1=vertcat(i1{:});end
     o2=cell(size(i1));
     for j1=1:length(o2)
      i2=i1(j1);
      if i2==-1; o2{j1}=comstr(ua.table{i1(j1),2}.value,-32);
      elseif i2==0;o2{j1}=''; % No unit
      elseif rem(i2,1)==0; o2{j1}=sprintf('%s',o1{i2}); 
      elseif fix(i2)==0
       st2=o1{round(rem(i2,1)*1000)};o2{j1}=sprintf('1/%s',st2); 
       try
        st3=textscan(st2,'%s','whitespace','/'); st3=st3{1};
        if length(st3)==2;o2{j1}=sprintf('%s/%s',st3{[2 1]});end
       end
      else
       o2{j1}=sprintf('%s/%s',o1{fix(i2)},o1{round(rem(i2,1)*1000)}); 
      end
     end
     o1=o2; % return strings
   end
  elseif nargout==0; %unitlable call
   o1(:,end+1)=num2cell([1:size(o1,1)]');disp(o1(:,[end end-2 1 2]));
   clear o1
  end


elseif comstr(Cam,'system')  % #UnitSystem' - - - - - - - - - - - - - - - - - -
  [CAM,Cam]=comstr(CAM,7);
  % Code TempCode (1 absolu, 2 relative) length, force, temp, temp-offset, time
  % to convert to SI units (divide by factor)   
  o1=struct('name','SI Meter (newton)', ...
           'data',[1 2 1 1 1 1 1]);

  o1(2)=struct('name','BG Foot (pound f)', ...
           'data',[2 2 3.28083989501312350e+00  2.24808943099710480e-01 ...
                   1.80000000000000000          4.59670000000000002e+02 1]);

  o1(3)=struct('name','MG Meter (kilogram f)', ...
           'data',[3 2 1 1.01971621297792830e-01 1 2.73149999999999980e+02 1]);

  o1(4)=struct('name','BA Foot (poundal)', ...
           'data',[4 2 3.28083989501312350  7.23301385120989430  ...
                   1.80000000000000000   4.59670000000000020e+02 1]);

  o1(5)=struct('name','MM mm (milli-newton)', ...
           'data',[5 2 1.00000000000000000e+03  1.00000000000000000e+03  ...
                   1.00000000000000000  2.73149999999999980e+02 1]);

  o1(6)=struct('name','CM cm (centi-newton)', ...
           'data',[6 2 1.00000000000000000e+02  1.00000000000000000e+02  ...
                   1.00000000000000000    2.73149999999999980e+02 1]);

  o1(7)=struct('name','IN Inch (pound f)', ...
           'data',[7 2 3.93700787401574810e+01  2.24808943099710470e-01  ...
                   1.80000000000000000     4.59670000000000020e+02 1]);

  o1(8)=struct('name','GM mm (kilogram f)', ...
           'data',[8 2 1.00000000000000000e+03  1.01971621297792830e-01  ...
                   1.00000000000000000       2.73149999999999980e+02 1]);
  o1(9)=struct('name','US user','data',[9 2 0 0 0 0 0]);
  o1(10)=struct('name','TM mm (Newton)', ...
           'data',[9 2 1.00000000000000000e+03  1.00000000000000000e+00  ...
                   1.00000000000000000  2.73149999999999980e+02 1]);
  o1(11)=struct('name','MU microns GPa', ...
           'data',[10 2 1.00000000000000000e+06  1.00000000000000000e+3  ...
                   1.00000000000000000  2.73149999999999980e+02 1]); % micro-m, GPascal, C, s ==> nano Joule
            
  % http://www.sdrl.uc.edu/sdrl/referenceinfo/universalfileformats/file-format-storehouse/universal-dataset-number-164
  if carg<=nargin; % Deal with proper code selection
      i1=varargin{carg};carg=carg+1;
      if isfield(i1,'code'); o1=o1(i1.code(1));end
  end
  if nargout==2; o2=fe_mat('unitlabel'); end
  if comstr(Cam,'name'); o1={o1.name}; end

else;sdtw('''Unit%s'' not know',Cam);
end

%% #type (Sub) ---------------------------------------------------------------
% [m_fun,unit,substr,subi]=fe_mat('typep',pl(1,2))
% fe_mat('typemstring',pl(1,2))
% fe_mat('typemstringl',pl(1,2)) % Long string
elseif comstr(Cam,'type') 

 if strcmpi(Cam,'types') % search available material types
 o1=matgui('search m_*'); return
 end
 
 [CAM,Cam]=comstr(CAM,5); if isempty(Cam); Cam='m'; end
  st1={'SI','BG','MG','BA','MM','CM','IN','GM','US'};
  r1=varargin{carg};carg=carg+1;
  if isnumeric(r1)
   %% type to result
   if size(r1,2)>1; r1=r1(1,2);
   elseif size(r1,1)>1 % fe_mat('typemstring',model.pl(:,2))
    o1=cell(size(r1)); 
    for j1=1:size(r1,1); o1{j1}=fe_mat(varargin{1},r1(j1));end
    return;
   end 
   if Cam(1)=='p'; [o1,o2,o3]=topType(r1);
   else [o1,o2,o3]=tomType(r1);
    if r1==1; o1='m_elastic';o2=1;o3=1;
    elseif r1<0||r1> 4.503599627370496e+001
     o1='m_null';o2=0;o3=0;
    else
     try
      st=dec2base(round(1e14*r1),36);
      if length(st)>2
       o2=abs(st(end-1))-48; o3=abs(st(end))-48;;
       o1=['m_' lower(st(1:end-2))];
      else;o2=0; o3=0; o1='';
      end
     catch
      o1='m_null';o2=0;o3=0;
     end
    end  
   end % m or p
   if nargout==4  % return [m_fun,unit,substr,subi]
    if any(o2==(1:9)); o2=st1{o2}; else;o2='US';end;
    o4=o3;
    try; 
     if ~strcmp(o1,'m_null')&&~isempty(o1)&&exist(o1,'file')
      o3=feval(o1,'subtypestring',o3); 
     elseif exist(['p' o1(2:end)],'file');% Mat implemented in p file
      o3=feval(['p' o1(2:end)],'subtypestringmat',o3); 
     else;o3=o1;
     end
    catch; o3=o1; 
    end
   elseif comstr(Cam,'mstringl')||comstr(Cam,'pstringl')
    if any(o2==(1:9)); o2=st1{o2}; else;o2='US';end;
       o1=sprintf('%s.%i.%s',o1,o3,o2);
   elseif comstr(Cam,'mstring')||comstr(Cam,'pstring')
       o1=sprintf('%s.%i',o1,o3);
   elseif comstr(Cam,'sub'); o1=o3;
   end
  else
   %% code properties into type number
   if ~comstr(r1,'m_') && ~comstr(r1,'p_')
    error('Material functions should start with ''m_''');
   end
   r1=r1(3:end);
   if ~isempty(r1)&&r1(end)==';'; sil=1; r1(end)=[]; else; sil=0; end
   if carg<=nargin; st=varargin{carg};carg=carg+1;else;st='1'; end
   if isnumeric(st); st=num2str(st);
   elseif any(strncmpi(st,st1,2)); st=num2str(find(strncmpi(st,st1,2)));
   else
    st='9'; if ~sil;  fprintf('fe_mat coding unit system %s as US\n',st); end
   end
   if length(st)>1; error('Unit code is between 0 and 9');end
   r1(end+1)=st;
   if carg<=nargin; st=varargin{carg};carg=carg+1;else;st='1'; end
   if isnumeric(st); st=num2str(st);end
   if length(st)>1; error('Sub type code is between 0 and 9');end
   r1(end+1)=st;
   o1=base2dec(comstr(r1,-271),36)*1e-14;
  end

%% #GetPlIl - - - - - - - - - - - - -
elseif comstr(Cam,'plil'); [CAM,Cam]=comstr(CAM,5);% - - - - - - - - - - - - -

ID=varargin{2};pl=varargin{3};il=varargin{4};
st='';
if ~isempty(pl)
 i1=find(pl(:,1)==ID(1));
 if isempty(i1)&&ID(1)~=0; st=sprintf('MatId %i not found PL',ID(1));pl=[];
 elseif ID(1)==0; pl=[]; else;pl=pl(i1,:); end
end

if ~isempty(il)
 i2=find(il(:,1)==ID(2));
 if isempty(i2)&&ID(2)~=0 
   st=sprintf('%s, ProId %i not found IL',st,ID(2));il=[ID(2) 0 1];
 elseif ID(2)==0; il=[0]; else;il=il(i2,:);
 end
else
 if ID(2)==0; ID(2)=ID(1); end
 il=[ID(2) 0 1];
end

o3=[];
if ~isempty(pl)&&comstr(Cam,'of_mk')
  [o1,o2]=fe_mat(Cam,pl,il); % [Constit,iopt]
  if size(o2,1)==1&&length(o2)>1; o2=o2(:);end
  o2(end+1,1)=-9999; % This is a flag to use Modulef Elements
else
  if ~isempty(st); fprintf('%s\n',st);end
  o1=pl; o2=il;  if ~isempty(pl); pl=pl(3:end); end  
end

%% #defaultil #defaultPl
%% -------------------------------------------------------------------------
elseif comstr(Cam,'default'); [CAM,Cam]=comstr(CAM,8);
 [CAM,Cam,RunOpt.list]=comstr('-list',[-25 3],CAM,Cam);
 [CAM,Cam,RunOpt.Silent]=comstr(';',[-25 3],CAM,Cam);
 [CAM,Cam,RunOpt.Unit]=comstr('-unit',[-25 4],CAM,Cam);
  list={
    'bar1'    ,'m_elastic(''dbval steel -unitSI'')','[1 fe_mat(''p_beam'',''SI'',1) 0 0 0 1]' 
    'beam1'   ,'m_elastic(''dbval steel -unitSI'')','p_beam(''defaultSmart'')'
    'beam3'   ,'m_elastic(''dbval steel -unitSI'')','p_beam(''defaultSmart'')'
    'beam1t'  ,'m_elastic(''dbval steel -unitSI'')','p_beam(''defaultSmart'')'
    'celas'   ,'[]'                                ,'p_spring(''defaultSI'')'
    'cbush'   ,'[]'                                ,'p_spring(''defaultSI'')'
    'line2'   ,'m_elastic(''dbval air -unitSI'')'  ,'p_solid(''dbval d3 13 -unitSI'')'
    'flui4'   ,'m_elastic(''dbval air -unitSI'')'  ,'p_solid(''dbval d3 -3 -unitSI'')'
    'flui6'   ,'m_elastic(''dbval air -unitSI'')'  ,'p_solid(''dbval d3 -3 -unitSI'')'
    'flui8'   ,'m_elastic(''dbval air -unitSI'')'  ,'p_solid(''dbval d3 -3 -unitSI'')'
    'fsc3'    ,'[]'                        ,'p_solid(''dbval fsc1 -unitSI'')'
    'hexa20'  ,'m_elastic(''dbval steel -unitSI'')','p_solid(''dbval d3 -3 -unitSI'')'
    'hexa20b' ,'m_elastic(''dbval steel -unitSI'')','p_solid(''dbval d3 -3 -unitSI'')'
    'hexa8'   ,'m_elastic(''dbval steel -unitSI'')','p_solid(''dbval d3 -3 -unitSI'')'
    'hexa8b'  ,'m_elastic(''dbval steel -unitSI'')','p_solid(''dbval d3 -3 -unitSI'')'
    'penta15' ,'m_elastic(''dbval steel -unitSI'')','p_solid(''dbval d3 -3 -unitSI'')'
    'penta15b','m_elastic(''dbval steel -unitSI'')','p_solid(''dbval d3 -3 -unitSI'')'
    'penta6'  ,'m_elastic(''dbval steel -unitSI'')','p_solid(''dbval d3 -3 -unitSI'')'
    'penta6b' ,'m_elastic(''dbval steel -unitSI'')','p_solid(''dbval d3 -3 -unitSI'')'
    'pyra13'  ,'m_elastic(''dbval steel -unitSI'')','p_solid(''dbval d3 -3 -unitSI'')'
    'pyra5'   ,'m_elastic(''dbval steel -unitSI'')','p_solid(''dbval d3 -3 -unitSI'')'
    'q4p'     ,'m_elastic(''dbval steel -unitSI'')','p_solid(''dbval d2 -3 -unitSI'')'
    'q4pb'    ,'m_elastic(''dbval steel -unitSI'')','p_solid(''dbval d2 -3 -unitSI'')'
    'q8p'     ,'m_elastic(''dbval steel -unitSI'')','p_solid(''dbval d2 -3 -unitSI'')'
    'q8pb'    ,'m_elastic(''dbval steel -unitSI'')','p_solid(''dbval d2 -3 -unitSI'')'
    'quad4'   ,'m_elastic(''dbval steel -unitSI'')','p_shell(''dbval kirchhoff .1 -punitSI'')' %'p_solid(''dbval d2 -3'')'
    'quad9'   ,'m_elastic(''dbval steel -unitSI'')','p_shell(''dbval kirchhoff .1 -punitSI'')' %'p_solid(''dbval d2 -3'')'
    'quadb'   ,'m_elastic(''dbval steel -unitSI'')','p_shell(''dbval kirchhoff .1 -punitSI'')' %'p_solid(''dbval d2 -3'')'
    't3p'     ,'m_elastic(''dbval steel -unitSI'')','p_solid(''dbval d2 -3 -unitSI'')'
    't3pb'    ,'m_elastic(''dbval steel -unitSI'')','p_solid(''dbval d2 -3 -unitSI'')'
    't6p'     ,'m_elastic(''dbval steel -unitSI'')','p_solid(''dbval d2 -3 -unitSI'')'
    't6pb'    ,'m_elastic(''dbval steel -unitSI'')','p_solid(''dbval d2 -3 -unitSI'')'
    'tetra10' ,'m_elastic(''dbval steel -unitSI'')','p_solid(''dbval d3 -3 -unitSI'')'
    'tetra10b','m_elastic(''dbval steel -unitSI'')','p_solid(''dbval d3 -3 -unitSI'')'
    'tetra4'  ,'m_elastic(''dbval steel -unitSI'')','p_solid(''dbval d3 -3 -unitSI'')'
    'tetra4b' ,'m_elastic(''dbval steel -unitSI'')','p_solid(''dbval d3 -3 -unitSI'')'
    'tria3'   ,'m_elastic(''dbval steel -unitSI'')','p_shell(''dbval kirchhoff .1 -punitSI'')'
    'tria6'   ,'m_elastic(''dbval steel -unitSI'')','p_shell(''dbval kirchhoff .1 -punitSI'')'};
  [CAM,Cam,i1]=comstr('-solid',[-25 3],CAM,Cam);matM=[];
  if i1
   list(strncmpi(list(:,3),'p_shell',6),3)={'p_solid(''dbval d2 -3 -unitSI'')'};
  end
  model=[];if carg<=nargin;model=varargin{carg};carg=carg+1;end
  if carg<=nargin;RunOpt=sdth.sfield('addmissing',RunOpt,varargin{carg});carg=carg+1;end
  warn={};
  if isempty(RunOpt.Unit)&&isfield(model,'unit')&&~isempty(model.unit)
   li=fe_mat('unitsystem'); li={li.name}; li=cellfun(@(x)x(1:2),li,'uni',0);
   li(ismember(li,'US'))=[];
   if ismember(upper(model.unit),li); RunOpt.Unit=upper(model.unit);
   else; warn{1}=sprintf('model unit %s is unknown or user defined(US)',model.unit);
   end
  end
  if ~isempty(RunOpt.Unit)
   list(:,2)=strrep(list(:,2),'SI',upper(RunOpt.Unit));
   list(:,3)=strrep(list(:,3),'SI',upper(RunOpt.Unit));
  end
  if sp_util('issdt'); list(:,3)=strrep(list(:,3),'p_shell(''dbval','p_shell(''dbval -f5'); end
  if RunOpt.list; o1=list; return; end
 
 RunOpt.typ=Cam;
 if comstr(RunOpt.typ,'il');     % pro default
  RunOpt.mpid=2; RunOpt.st='ProId';
 elseif comstr(RunOpt.typ,'pl'); % mat default
  RunOpt.mpid=1; RunOpt.st='MatId';
  if isfield(RunOpt,'nmap')
     matM=useOrDefault(RunOpt.nmap,'Map:MatDB');  
     if ~isempty(matM);matM=cell(matM);
      for j1=1:size(matM,1)
       if ischar(matM{j1,2})
        st=sdtm.urnCb(matM{j1,2});matM{j1,2}=feval(st{:});
       end
      end
     end
  end
 else
  if isempty(model); sdtw('_err','you must give a model'); end
  model=fe_mat('defaultil',model); model=fe_mat('defaultpl',model);
  o1=model;return;
 end
 
 if isempty(model)
  if comstr(RunOpt.typ,'il'); o1=p_solid('database');o1=o1(1);
  else;          o1=m_elastic('database'); o1=o1(1);
  end
 else
  r1=feutil(RunOpt.st,model); i1=r1==0&isfinite(model.Elt(:,1));
  if any(i1) % some elts are not assigned IDs, asgn new per group
   r2=feutil([RunOpt.st 'new'],model);
   [EG,nG]=getegroup(model.Elt); i1=find(i1);
   mpid=feutil('mpid',model); 
   for jG=1:nG
    % skip matid/Proid for mass
    % skip proid celas if pro exist
    [ElemF,i3,ElemP]=feutil('getelemf',model.Elt(EG(jG),:));
    switch ElemF
     case {'mass1','mass2','rigid'}; continue % skip
     case 'celas'; % keep only if no pro is present
      if comstr(lower(RunOpt.st),'matid'); continue; end
      i2=intersect(i1,EG(jG)+1:EG(jG+1)-1);
      if ~isempty(i2)&&size(model.Elt,2)>6
       i2=intersect(i2,find(~any(model.Elt(i2,7:end))));
      end
     case 'cbush'
      if comstr(lower(RunOpt.st),'matid'); continue; end
      i2=intersect(i1,EG(jG)+1:EG(jG+1)-1);
     case {'quad4','quadb','tria3','tria6'}
      % skip if used with p_contact for SDT
      i2=intersect(i1,EG(jG)+1:EG(jG+1)-1); % unassgn in group
      if comstr(lower(RunOpt.st),'matid')&&sp_util('issdt')
       p1=feval(ElemF,'prop'); p1=model.Elt(i2,p1(2)); p2=unique(p1); p2(p2==0)=[];
       il1=matgui('getstackil',model,p2); % contact has stack entry anyways
       il1=il1(cellfun(@(x)isfield(x,'il')&&length(x.il)>1,il1(:,3)),3);
       il1=[cellfun(@(x)x.il(1),il1,'uni',0) cellfun(@(x)fe_mat('typepstring',x.il(2)),il1,'uni',0)];
       i4=strncmp(il1(:,2),'p_contact',9);
       if any(i4);        i2(ismember(p1,cell2mat(il1(i4,1))))=[];       end
      end
     otherwise; i2=intersect(i1,EG(jG)+1:EG(jG+1)-1); % unassgn in group
    end
    if ~isempty(i2);  % set id and increment
     warn{end+1}=sprintf('Some elements with no %s, setting new %s %i in group %i',...
      RunOpt.st,RunOpt.st,r2,jG);
%      sdtw('_nb','Some elements with no %s, setting new %s %i in group %i',...
%       RunOpt.st,RunOpt.st,r2,jG);
     r1(i2)=r2; r2=r2+1;
    end
   end % jG
   mpid(:,RunOpt.mpid)=r1;
   model.Elt=feutil('mpid',model,mpid);
  end
  i1=unique(r1); i1(~i1)=[];
  plil=fe_mat(sprintf('get%s',RunOpt.typ),model);  
  if isempty(plil); plil=double.empty(0,2); end
  if size(plil,2)<2;plil=[];else;plil(plil(:,2)==0,:)=[];end
  
  if ~isempty(plil); i1=setdiff(i1,plil(:,1));end
  mpid=feutil('mpid',model);[EGroup,nGroup,RunOpt.GroupElemF]=getegroup(model.Elt);
  for j1=1:length(i1); 
   ind1=find(mpid(:,RunOpt.mpid)==i1(j1));%feutil(sprintf('findelt %s %i',RunOpt.st,i1(j1)),model);
   [ElemF,i3]=feutil('getelemf',model.Elt(EGroup(mpid(ind1(1),3)),:),mpid(ind1(1),3));
   %ind2=feutil(sprintf('findelt eltname %s &%s~=%i',ElemF,RunOpt.st,i1(j1)),model);
   ind2=find(mpid(:,RunOpt.mpid)~=i1(j1)&ismember(mpid(:,3),find(ismember(RunOpt.GroupElemF,ElemF))));
   % check if there is the same elt with another property
   i4=[];
   if ~isempty(ind2)&&any(~ismember(mpid(ind2,RunOpt.mpid),i1)) 
     ind2=ind2(~ismember(mpid(ind2,RunOpt.mpid),i1)); 
     i4=find(plil(:,1)==mpid(ind2(1),RunOpt.mpid)); 
   end
   if ~isempty(i4) % same elt with valid prop found
     plilj1=plil(i4(1),:); plilj1(1)=i1(j1);
     plil(end+1,1:length(plilj1))=plilj1;  %#ok<AGROW>
   else % default
    if strcmpi(ElemF,'SE');
     warn{end+1}=sprintf('SE proid %i not defined',i1(j1));
     continue
    end
     i3=find(strcmpi(ElemF,list(:,1)));plilj1=[];
     if ~isempty(matM)
      %% automate assignation of material properties from matM
      for j2=1:size(matM,1)
        if isnumeric(matM{j2,2})&&matM{j2,2}(1)==i1(j1)
          plilj1=fe_mat(['convert' model.unit],matM{j2,2});
        elseif isfield(matM{j2,2},'pl')&&matM{j2,2}.pl(1)==i1(j1)
          plilj1=matM{j2,2};
          if ~isempty(setdiff(fieldnames(plilj1),{'pl','type','name','unit'}))
           plilj1=feval(plilj1.type,['convert' model.unit],plilj1);
           r2=plilj1;r2.pl=r2.pl(1);plilj1=plilj1.pl;
           model=stack_set(model,'mat',matM{j2,1},r2);
          else
           plilj1=fe_mat(['convert' model.unit],matM{j2,2}.pl);
          end
        end
      end
     end
     if ~isempty(plilj1)
     elseif isempty(i3) 
       sdtw('_nb','unknown element %s',ElemF)
       if isfield(model,'unit');  st=sprintf('-unit%s',model.unit);
       else; st=[];
       end
       if comstr(RunOpt.typ,'il')
        plilj1=p_solid(sprintf('dbval %s %i d3 -3',st,i1(j1)));
       else
        plilj1=m_elastic(sprintf('dbval %s %i steel',st,i1(j1)));
       end
     else % default of list
       eval(sprintf('plilj1=%s;',list{i3(1),RunOpt.mpid+1}))
       % check if somethingelse is more revelant than default:
    
       if comstr(ElemF,'quad') % 2D ?
        if comstr(RunOpt.typ,'il') 
           n1=feutil(sprintf('GetNode proid%i',i1(j1)),model);
           if ~isempty(n1)&&~any(n1(:,7)); % orientation q4p-like should be performed xxx
               plilj1=p_solid('dbval d2 -3'); 
           end 
        end
       else% keep default   
       end
       
       if isempty(plil); i4=[]; 
       elseif isempty(plilj1); continue;
       elseif isfield(plilj1,'il'); plilj1=plilj1.il;
           i4=find(plil(:,2)==plilj1(2)); 
       else; i4=find(plil(:,2)==plilj1(2)); 
       end
       if ~isempty(i4) % there is the same mat subtype existing in pl
         plilj1=plil(i4(1),:);
       end
       plilj1(1)=i1(j1);
     end
     if isstruct(plilj1); plilj1=plilj1.(RunOpt.typ); end

     plil(end+1,1:length(plilj1))=plilj1; 
     warn{end+1}=sprintf('Defining default %s %i',RunOpt.st,i1(j1));
   end
  end
  if ~RunOpt.Silent; 
   if sp_util('issdt'); cellfun(@(x) sdtw('_nb',x),warn,'uni',0);
   else; disp(warn);
   end
  end
  o1=model;
  if ~isempty(plil);
    if isfield(model,'unit')&&~isempty(model.unit)&&~strncmpi(model.unit,'us',2)
     plil=fe_mat(sprintf('Convert %s;',model.unit),plil);% convert to current unit
     [r1,i1]=stack_get(model,lower(RunOpt.st(1:3))); 
     if ~isempty(r1) % edit convertion in stack as well!
      ils=cellfun(@(x)x.(RunOpt.typ)(1),r1(:,3));
      for j1=1:length(ils) % il in stack
       i2=ismember(plil(:,1),ils(j1));
       if any(i2)
        if length(r1{j1,3}.(RunOpt.typ))>1
            r1{j1,3}.(RunOpt.typ)=plil(i2,:); 
        end
        r1{j1,3}.unit=model.unit;
       end
      end
      o1.Stack(i1,:)=r1;
     end
    end
    o1.(RunOpt.typ)=plil;
  end
 end
%% #pu : Property unit type formatting fe_mat('put','m_elastic',3)
elseif comstr(Cam,'pu') % - - - - - - - - - - - - - - - - - - - - - - - -
    
 r1=feval(varargin{2},'propertyunittypecell',varargin{3:end});
 if ~iscell(r1)
  for j1=1:length(r1); fe_mat('pu',varargin{2},j1);  end
  return;
 end
 r1{2,3}=sprintf('fe_mat(''%s'',''SI'',%i)',varargin{2},varargin{3});
 if nargout==0; 
    r1(:,4)=num2cell(1:size(r1,1));r1=r1(:,[4 1 2 3])';
    fprintf('%-3i%-5s %8.3f %s\n',r1{:});fprintf('\n');
 end

%% #of_mk
elseif comstr(Cam,'of_mk') % - - - - - - - - - - - - - - - - - - - - - - - -

% this command reformats any SDT Material entry to the row input needed by
% of_mk calls to element functions
%
% type 1 : plane stress isotropic
%      2 : plane strain isotropic or 3-D isotropic
%      3 : 3-D anisotropic
%      4 : 2-D anisotropic
%      12 : reserved for acoustic fluid
%
% These are not coherent with SDT reference defined in m_elastic

pl=varargin{carg};carg=carg+1;
il=varargin{carg};carg=carg+1;
[CAM,Cam]=comstr(Cam,6);

o2=[pl(1) il(1) zeros(1,6)];
if comstr(Cam,'2d') % - - - - - - - - - - - - - - - - - - - - - - - - - -

  % INTEG Formatting for 2D elements by Amine
  %[MatId ProId NDOF NNode Type 0 Formulation IoptStiff IoptMass N]
  % of_mk_subs.f compute_matrix IOPT indices
  %              1     2    3   4  5           6        7        8 
  %
  % Values of Formulation  :
  %           IoptStiff    : always 3 (anisotropic call)
  %           IoptMass     : currently unused
  %                  N     : only used by q9a eraq2c.f element   

  if ~any(il(3)==[0 1 2]) 
    error('%i: unknown 2-D element formulation',il(3));
  end
  o2(7)=il(3)+1; % Formulation (plane or axisymetric)
  if length(il)>3; o2(8)=il(4);end % NFourier

elseif comstr(Cam,'3d') % - - - - - - - - - - - - - - - - - - - - - - - - - -

  % INTEG formatting for 3D elements by Marina Vidrascu
  % [MatId ProId NDOF NNode Type 0 Formulation Iopt]
  % of_mk_subs.f compute_matrix IOPT indices
  %              1     2    3   4  5           6        
  % IoptStiff  2: isotropic, 1: orthotropic

else % other % - - - - - - - - - - - - - - - - - - - - - - - - - -

end

[st,UnitCode,SubType]=fe_mat('type',pl(2));

switch st
case 'm_elastic'
 
 switch SubType
 case 1 % standard isotropic of_mk 

  if size(pl,2)<7; pl(7)=0;end
  if comstr(Cam,'2d') % For 2D only : constit=[rho eta E11 ... E33 a1 a2 a3 T0]

     %if o2(8)==0 o2(8)=o2(7);end 
     o1 = pl([5 7]); % [rho eta]
     E =pl(3); nu=pl(4); 

     switch il(3) % Formulations
     case 1 % - - - - - - - - - -  - - - - - - - plane stress
       o2(8)=3; % anisotropic element is always called
       C=E/(1.-nu^2);r1=[C C*nu C 0. 0. C*(1-nu)/2];
       o1(3:8)=r1;
       o1(9:12)=0; % a1 a2 a3 T0
     case 0 % - - - - - - - - - -  - - - - - - - plane strain
       o2(8)=3; % anisotropic element is always called
       C=E*(1-nu)/(1+nu)/(1-2*nu);
       r1=[C C*nu/(1-nu) C 0. 0. C*(1-2*nu)/2/(1-nu)];
       o1(3:8)=r1;  o1(9:12)=0; % a1 a2 a3 T0
     case 2 % - - - - - - - - - -  - - - - - - - axisymetric
 
       o1 = pl([5 7 3 4]);% xxx only isotropic case => o2(8)=0;
       o1(9:12)=0; % a1 a2 a3 T0

     end % Formulations 


  elseif  pl(6)==0 || abs(pl(6)-pl(3)/2/(1+pl(4)))<sqrt(eps)*pl(3) 

     o1 = pl([5 7 3 4]); o2(8)=2; % 3D isotropic : constit=[rho eta E nu] 

  else % G is not E/2/(1+nu) use anisotropic
     % 22/08/06 formula corrected for consistency when G ne E/2(1+n)
     % was r1=pl(6);r2=pl(3)*pl(4)/(1+pl(4))/(1-2*pl(4));
     % was o1=[pl([5 7]) r2+[2*r1 0 2*r1 0 0 2*r1] r1 r1 r1];
     r1=pl(3)/2/(1+pl(4));r2=pl(3)*pl(4)/(1+pl(4))/(1-2*pl(4));
     o1=[pl([5 7]) r2+[2*r1 0 2*r1 0 0 2*r1] pl([6 6 6])];

     o2(8)=1; % 3D anisotropic
  end

 case 4 % 2-D Anisotropic : car = [rho eta E11 ... E33 a1 a2 a3 T0]

    if length(pl)<9; pl(14)=0.; end 
    o1=pl([9 14 3:8 10:13]);
    o2(8)=3; % anisotropic

 case 3 
 % 3-D Anisotropic : car = [rho eta G11 G12 G22 G13 G23 G33 G14 ...  G66]
 % this has not been rechecked, hence the warning
 warning('You should use variable integration rule elements non isotropic problems')
   if any(pl(2+[7 8 9 11:14 16:20]))
    % il(3) calls the old .m file assembly with fully anisotropic material
    % Constit has [ rho eta dd(:)'] 
    dd=zeros(6);dd(triu(ones(6))~=0)=pl(3:23);
    dd=dd+tril(dd',-1); if length(pl)<25; pl(25)=0; end
    o1=[pl(24:25) dd(:)'];
    o2(8)=4; % anisotropic
   else % actually orthotropic [... eta]
    if length(pl)<40; pl(40)=0;end
    o1=[pl(1) 3 pl([24 3:8 12 17 23 40])];
    o2(8)=3;        % orthotropic
   end

 otherwise
   warning('Not an expected situation')        
 end

otherwise; error('Material function not supported by of_mk');
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'info'); % lists all material functions and subtypes
    
  for j0=1:2;
   if j0==1;st=matgui('search m_*');
   else;st=matgui('search p_*');
   end
  for j1=1:length(st)
    for j2=1:9
     [r1,r2]=feval(st{j1},'PropertyUnitType cell',j2); 
     if length(r1)>2; 
         fprintf('%s Subtype %i\n',st{j1},j2);
         for j3=3:size(r1,1); fprintf('   (%i) %s\n',j3,r1{j3,1});end
     end
     if ~isempty(r2);
         for j3=1:size(r2,1)
          fprintf('   (%i+i*%i) %s\n',size(r1,1)+j3,size(r2,1),r2{j3,1}); 
         end
     end
    end
  end
  end
%% #EndCommands
elseif comstr(Cam,'@');o1=eval(CAM);
    if nargin>1;feval(o1,varargin{2:end});end
else;error('''%s'' Not a valid fe_mat command',CAM);
end


% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
% THIS SECTION IS OBSOLETE AND WILL EVENTUALLY DISAPPEAR IN P_* functions
else % ElemF is not a string

if nargin<4; error('Not the right number of inputs');end
ElemF=varargin{1};ID=varargin{2};pl=varargin{3};il=varargin{4};
if nargin>4; t=varargin{5}; end

% find the good properties

if ~isempty(pl)
 i1=find(pl(:,1)==ID(1));
 if isempty(i1)&&ID(1)~=0; warning('MatId %i not found PL',ID(1));pl=[];
 elseif ID(1)==0; pl=[]; else;pl=pl(i1,:); end
end

if ~isempty(il)
 i2=find(il(:,1)==ID(2));
 if isempty(i2)&&ID(2)~=0 ;warning('ProId %i not found IL',ID(1));il=[];
 elseif ID(2)==0; il=[]; else;il=il(i2,:); end
end

if size(pl,1)>1 
  pl=pl(1,:);fprintf(1,'taking first property with MatId %i\n',pl(1));
end
if size(il,1)>1; il=il(1,:);
  fprintf(1,'taking first property with ProId %i\n',il(1));
end

if ~isempty(pl); [st,i1,i2]=fe_mat('typem',pl(2));end
if isempty(pl) 
elseif pl(2)==1||i2==1 % standard isotropic material
   if length(pl)==5; 	pl=[pl 0]; 		 end
   if pl(6)==0; 		pl(6)=pl(3)/2/(1+pl(4)); end
else 
 error('MatId %i : %s,%i,%i not supported for beams',pl(1),st,i1,i2);
end

% element specific treatment of element properties

if ElemF==1; error('was moved to p_solid buildconstit')
elseif ElemF==2 % comstr(ElemF,'shell1')

   %  1    2   3      4   5   6      7     8     9      10   11 12 13  14
   % ProId 2 Formul Drill 0 Thick k(TS/T) MID2 12l/T^3 MID3 NSM Z1 Z2 MID4

   if isempty(pl); error('No declared material property');end
   il(15)=0;
   if il(2)==2
   else;[st,i1,i2]=fe_mat('typep',il(2)); 
    if ~strcmp(st,'p_shell')
     error('ProId %i is not a shell property (%s)',il(1),st);
    end
    if i2~=1; sdtw('p_shell subtype %i is supported by q4cs elements',i2);end 
   end
   if nargin>4 && ~isempty(t); il(6)=t; end
   if il(4)==0; il(4)=1; end % drilling DOF stiffness
   o1 = pl; o2 = il;
   o3 = pl(3)/(1-pl(4)^2); % membrane
   o3 = [o3 pl(4)*o3 0;      pl(4)*o3 o3 0;      0 0 (1-pl(4))/2*o3]*il(6);
   o4 = o3*il(6)^2/12;     % bending
   if il(9)~=0; o4=o4*il(9); end  % ratio for bending moment of inertia
   if any(il([8 10 12 13 14]))
     %sdtw('_nb','Unsupported properties in shell declaration');
   end
   if length(il)<7; il(7)=5/6; elseif il(7)==0; il(7)=5/6;end
   o5 = eye(2)*pl(6)*il(7)*il(6); % shear

%a=0; % a=0 --> Deformations planes
%     % a=1 --> Contraintes planes

%H1 = pl(1,1)*(1-a*pl(1,2))/((1+pl(1,2))*(1-pl(1,2)-a*pl(1,2)));
%Hm = H1*[1 pl(1,2)/(1-a*pl(1,2)) 0;
%         pl(1,2)/(1-a*pl(1,2)) 1 0;
%         0 0 (1-pl(1,2)-a*pl(1,2))/(2*(1-a*pl(1,2)))];

elseif ElemF==3 %comstr(ElemF,'solid1')

    o1 = pl; o2 = il;

    if pl(6)==0; pl(6)=(pl(3)/2/(1+pl(4))); end
    r1=pl(6);r2=pl(3)*pl(4)/(1+pl(4))/(1-2*pl(4));
    dd=zeros(6);
    dd([1 2 3 7 8 9 13 14 15])=r2+[2*r1 0 0  0 2*r1 0  0 0 2*r1];
    dd([22 29 36])=pl(6);
    o3=dd;

else;error('Not a supported material/property call');
end

% THIS SECTION IS OBSOLETE AND WILL EVENTUALLY DISAPPEAR IN P_* functions
end % ElemF is not a string

%% #SubFunc ------------------------------------------------------------------

function [o1,o2,o3]=topType(r1);
%% #topType
    if     r1==1; o1='p_beam';o2=1;o3=1;
    elseif r1==2; o1='p_shell';o2=1;o3=1;
    elseif r1<0||r1> 4.503599627370496e+001
     o1='p_null';o2=0;o3=0;
    else
     try
      st=dec2base(round(1e14*r1),36);
      if length(st)>2
       o2=abs(st(end-1))-48; o3=abs(st(end))-48;
       o1=['p_' lower(st(1:end-2))];
      else;o2=0; o3=0; o1='';
      end
     catch
      o1='p_null';o2=0;o3=0;
     end
    end

function [o1,o2,o3]=tomType(r1);
    if r1==1; o1='m_elastic';o2=1;o3=1;
    elseif r1<0||r1> 4.503599627370496e+001
     o1='m_null';o2=0;o3=0;
    else
     try
      st=dec2base(round(1e14*r1),36);
      if length(st)>2
       o2=abs(st(end-1))-48; o3=abs(st(end))-48;;
       o1=['m_' lower(st(1:end-2))];
      else;o2=0; o3=0; o1='';
      end
     catch
      o1='m_null';o2=0;o3=0;
     end
    end 
    
%% #MergePlIl Merge different pl
function  model=MergePlIl(type,model,el0);
% only non intersecting ids of el0 are added, no renumbering here
pl=[];
if isstruct(el0)
 if isfield(el0,type) && ~isempty(el0.(type));  pl=el0.(type); end
else; pl=el0;
end
if ~isfield(model,type); model.(type)=[]; end
if ~isempty(model.(type))&&~isempty(pl) % both have entries
 [i1,i2]=ismember(pl(:,1),model.(type)(:,1));i1=find(i1); % intersecting ids
 pl1=pl(i1,:);i3=false(size(i1));
 for j1=1:size(pl1,1) % check if same
  pl2=pl1(j1,:);pl2=pl2(1:find(any(pl2,1),1,'last'));
  if size(pl2,2)<=size(model.(type),2) % Same material
   i3(j1)=any(model.(type)(i2(i1(j1)),1:size(pl2,2))-pl2,2);
  end
 end
 if any(i1) % intersecting entries
  if any(i3) % material is not the same, warn and ignore
   if strcmpi(type,'pl');st='MatId';else;st='ProId';end
   sdtw('_nb','Ignoring %s %s in second model', ...
    st,comstr(pl(i1(i3),1),-30));
  end
  pl(i1,:)=[];
 end
end

if ~isempty(pl); model.(type)(end+(1:size(pl,1)),1:size(pl,2))=pl; end
%model.pl=unique(model.pl,'rows');

% -------------------------------------------------------------------------
%% #field_interp Interpolate fields if needed
% tested in t_thermal('fieldinterp')
% fe_mat('@field_interp',';');
function  [mat,model,i3]=field_interp(mat,model);
 persistent silent;
 if nargin==1&&ischar(mat);
     silent=mat;return;
 end
 % below should be moved to elem0
 if isfield(mat,'pl');RO.ty='pl';else;RO.ty='il';end
 r1=mat.(RO.ty); nodeEt=[];i3=[];
 if isscalar(r1);return;end
 if isfield(mat,'EC')&&isfield(mat.EC,'MatPropertyUnit'); % Self contained structure
     st1=mat.EC.MatPropertyUnit;st1=st1(:,1);
     for j1=1:size(st1,1);st1{j1}=sscanf(st1{j1},'%s',1);end
     if isfield(mat.EC,'nodeEt');nodeEt=mat.EC.nodeEt;end
 else; st1={}; % delay st1 generation to ensure it is needed
 end
 for j1=1:size(st1,1);st1{j1}=sscanf(st1{j1},'%s',1);end
 st=fieldnames(mat);r1={};
 for j1=1:length(st) % st {PropName FieldName Table}
  r2=mat.(st{j1});
  if isfield(r2,'X')&&isfield(r2,'Y') % field is a table
   if isempty(st1) % get fields of current mat
    [st1,unu1,i1]=fe_mat('type',mat.(RO.ty)(2));
    RO.constitLab=feval(p_solid('@ConstitLab'),sprintf('%s.%i',st1,i1));
    st1=feval(mat.type,'propertyunittype cell',i1);st1=st1(:,1);
    for j2=1:size(st1,1);st1{j2}=sscanf(st1{j2},'%s',1);end
   end
   if isfield(r2,'Xlab'); st{j1,2}=r2.Xlab{1}; else; st{j1,2}='T';end
   st{j1,3}=find(strcmp(st1,st{j1,1})); % Interpolate from plLab
   if isempty(st{j1,3}) % interpolate from constitLab
     st{j1,3}=find(strcmp(RO.constitLab,st{j1,1}));
   end
   if isempty(st{j1,3});
    if strcmpi(st{j1,1},'data')||strcmpi(st{j1,2},'eltid');% No warning if field is data
    elseif isempty(silent)||~strcmp(silent,';')
     sdtw('_nb','%s(%s) not in ConstitLab, skipped',st{j1,1:2});
    end
    continue;
   end
   if nargin==2 % pl interpolation
    r3=stack_get(model,'info',st{j1,2},'get');
    if isempty(r3);
     r3=stack_get(model,'info',sprintf('Ref%s',st{j1,2}),'get');
    end
    if (~isfield(RO,st{j1,2})||isempty(RO.(st{j1,2})))&&~isempty(r3)
     st{j1,4}=sprintf('@ %s=%g',st{j1,2},r3);
    end
    RO.(st{j1,2})=r3;
    X=r2.X;if iscell(X);X=X{1};end
    if isempty(RO.(st{j1,2})); RO.(st{j1,2})=min(X);% Define Field
     st{j1,4}=sprintf('@ %s=%g',st{j1,2},RO.(st{j1,2}));
    end
    if  st{j1,3}>size(mat.(RO.ty),2)
    elseif any(mat.(RO.ty)(st{j1,3})==[-2 -3]) % CTable table build
     st{j1,1}=''; st{j1,2}='';% Don't interp yet
    else
     r2.Extrap=' '; % Do not extrapolate by 0
     i3(st{j1,3})=1; % this was interpolated
     mat.(RO.ty)(st{j1,3})=fe_curve('returny',r2,RO.(st{j1,2}));
    end
   else; 
   % #CTable building [i1 xi si xstartpos Nx nodeEfield constit(pos_Matlab)]
    j3=find(strcmp(st1,st{j1}));
    r1{j1,1}=[0 0 0 0 length(r2.X) ...
     find(mat.EC.nodeEt==comstr(st{j1,2},-32)) ... % nodeEfield
     find(strcmp(st1,st{j1}))];
    r1{j1,2}=[r2.X(:);r2.Y(:)];
   end
  else; st{j1}='';
  end
 end
 st(cellfun('isempty',st(:,1)),:)=[]; 
 if nargin==1 &&~isempty(nodeEt) % allow fields obtained from field at node
   for j1=reshape(find(mat.EC.constit==-3),1,[])
     j2= find(nodeEt==comstr(st1{j1},-32)); % st1 predefined
     if isempty(j2); 
         error('''%s'' not found in InfoAtNode/nodeEt field',st1{j1});
     end
     r1{j1,1}=[0 0 0 0 -3 ...
          j2 ... % nodeEfield
           j1];r1{j1,2}=[];
   end
 end
 if ~isempty(r1) % Combine CTable for interp
    r1(cellfun('isempty',r1(:,1)),:)=[];
    i1=[0;cumsum(cellfun('length',r1(:,2)))];i1(end)=[];
    r2=vertcat(r1{:,1});r2(:,4)=i1+1+numel(r2);
    mat.CTable=[size(r1,1);reshape(r2',[],1);
        vertcat(r1{:,2})];
 end
 if ~isempty(st);
  if strcmp(RO.ty,'pl');st1='MatId';else;st1='ProId';end
  st2=st(:,1:2)'; if size(st,2)<4;st(:,4)={''};end
  st3=''; if isfield(mat,'name');st3=sprintf('(%s)',mat.name);end
  if ~isempty(silent)&&strcmp(silent,';')
      % feval(fe_mat('@field_interp'),';')
  else
   fprintf('%s=%i %s interpolated %s%s\n', ...
      st1,mat.(RO.ty)(1),st3,sprintf('%s(%s) ',st2{:}),sprintf(' %s',st{:,4}));
  end
 end
 

%============================================================
%Material               : 1 - GENERIC_ISOTROPIC_STEEL
%============================================================
%------------------------------------------------------------------------------
%Material Property      : 1 - MODULUS OF ELASTICITY
%  Dim : PRESSURE    Unit : PASCAL
%  Scalar Value :  2.068000E+11 
%Material Property      : 3 - MASS DENSITY
%  Dim : MASS DENSITY    Unit : KILOGRAM/METER^3
%  Scalar Value :  7.820000E+03 
%Material Property      : 7 - THERMAL EXPANSION REFERENCE TEMPERATURE
%  Dim : TEMPERATURE    Unit : KELVIN
%  Scalar Value :  2.950000E+02 
%------------------------------------------------------------------------------
%Material Property      : 1 - MODULUS OF ELASTICITY
%  Dim : PRESSURE    Unit : POUND FORCE/FOOT^2
%  Scalar Value :  4.319109E+09 
%Material Property      : 3 - MASS DENSITY
%  Dim : MASS DENSITY    Unit : LBF SEC^2/FT^4
%  Scalar Value :  1.517332E+01 
%Material Property      : 7 - THERMAL EXPANSION REFERENCE TEMPERATURE
%  Dim : TEMPERATURE    Unit : RANKINE
%  Scalar Value :  5.310000E+02 
%------------------------------------------------------------------------------
%Material Property      : 1 - MODULUS OF ELASTICITY
%  Dim : PRESSURE    Unit : KILOGRAM FORCE/M^2
%  Scalar Value :  2.108773E+10 
%Material Property      : 3 - MASS DENSITY
%  Dim : MASS DENSITY    Unit : KGF SEC^2/M^4
%  Scalar Value :  7.974181E+02 
%Material Property      : 7 - THERMAL EXPANSION REFERENCE TEMPERATURE
%  Dim : TEMPERATURE    Unit : KELVIN
%  Scalar Value :  2.950000E+02 
%------------------------------------------------------------------------------
%Material Property      : 1 - MODULUS OF ELASTICITY
%  Dim : PRESSURE    Unit : POUNDAL/FOOT^2
%  Scalar Value :  1.389632E+11 
%Material Property      : 3 - MASS DENSITY
%  Dim : MASS DENSITY    Unit : POUND MASS/FOOT^3
%  Scalar Value :  4.881870E+02 
%Material Property      : 7 - THERMAL EXPANSION REFERENCE TEMPERATURE
%  Dim : TEMPERATURE    Unit : RANKINE
%  Scalar Value :  5.310000E+02 
%------------------------------------------------------------------------------
%Material Property      : 1 - MODULUS OF ELASTICITY
%  Dim : PRESSURE    Unit : MILLINEWTON/MILLIMETER^2
%  Scalar Value :  2.068000E+08 
%Material Property      : 3 - MASS DENSITY
%  Dim : MASS DENSITY    Unit : KILOGRAM/MILLIMETER^3
%  Scalar Value :  7.820000E-06 
%Material Property      : 7 - THERMAL EXPANSION REFERENCE TEMPERATURE
%  Dim : TEMPERATURE    Unit : KELVIN
%  Scalar Value :  2.950000E+02 
%------------------------------------------------------------------------------
%Material Property      : 1 - MODULUS OF ELASTICITY
%  Dim : PRESSURE    Unit : CENTINEWTON/CENTIMETER^2
%  Scalar Value :  2.068000E+09 
%Material Property      : 3 - MASS DENSITY
%  Dim : MASS DENSITY    Unit : KILOGRAM/CENTIMETER^3
%  Scalar Value :  7.820000E-03 
%Material Property      : 7 - THERMAL EXPANSION REFERENCE TEMPERATURE
%  Dim : TEMPERATURE    Unit : KELVIN
%  Scalar Value :  2.950000E+02 
%------------------------------------------------------------------------------
%Material Property      : 1 - MODULUS OF ELASTICITY
%  Dim : PRESSURE    Unit : POUND FORCE/INCH^2
%  Scalar Value :  2.999380E+07 
%Material Property      : 3 - MASS DENSITY
%  Dim : MASS DENSITY    Unit : LBF SEC^2/IN^4
%  Scalar Value :  7.317372E-04 
%Material Property      : 7 - THERMAL EXPANSION REFERENCE TEMPERATURE
%  Dim : TEMPERATURE    Unit : RANKINE
%  Scalar Value :  5.310000E+02 
%------------------------------------------------------------------------------
%Material Property      : 1 - MODULUS OF ELASTICITY
%  Dim : PRESSURE    Unit : KILOGRAM FORCE/MM^2
%  Scalar Value :  2.108773E+04 
%Material Property      : 3 - MASS DENSITY
%  Dim : MASS DENSITY    Unit : KGF SEC^2/MM^4
%  Scalar Value :  7.974181E-10 
%Material Property      : 7 - THERMAL EXPANSION REFERENCE TEMPERATURE
%  Dim : TEMPERATURE    Unit : KELVIN
%  Scalar Value :  2.950000E+02 
%------------------------------------------------------------------------------

