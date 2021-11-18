function [out,out1,out2,out3]=fe_curve(varargin)

%FE_CURVE Curve handling and signal processing utilities 
%
%Supported curve handling  commands are
%
%    DimLab=fe_curve('DataType[Cell]',DesiredType,unit);
%    curve =fe_curve('GetCurve',model,curve_name);
%    out=fe_curve('GetId Id',model);
%    model=fe_curve(model,'Set',name,curve_def);
%    out=fe_curve('read /path/filename.cyt');
%    Y = fe_curve('returny',curve,X)
%    Y = fe_curve('returny -lin',curve,X) % nota : default is quadratic
%
%    curve=fe_curve('getcurve',model,curve_name)
%    list=fe_curve('getcurve',model)
%    model=fe_curve(model,'set',curve_name,data_structure)
%    model=fe_curve(model,'set',curve_name,string)
%
%    fe_curve('plot',curve);
%
%Signal processing  commands are
%    out=fe_curve('FFT',frame); % discrete fourier transform
%    out=fe_curve('TimeFreq',Input,Transfert);
%
%    out=fe_curve('BandPass Hz f_min f_max',frames);
%    out=fe_curve('BandPass Rd w_min w_max',frames);
%    out=fe_curve('window window_name',Nb_pts);
%    out=fe_curve('H1H2 input_channels',frame,window,weighing);
%                  window and weighing can be omitted
%    out=fe_curve('ResSpectrum Type',frames,freq,damp);
%    out=fe_curve('Noise',Nw_pt,sf,spec);
%    out=fe_curve('TestFrame');
%    out=fe_curve('TestFunc');% Fun=cos, sin, etc.
%    out=fe_curve('TestFunc',t);%t:time vector, ex: t=linspace(0,1,199)
%
% See also sdtweb('fe_curve'), sdtweb('curve')

%	Etienne Balmes, Mathieu Corus, J.-P. Bianchi, G. Vermot des Roches
%       Copyright (c) 2001-2019 by SDTools and INRIA, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

if nargin==0; CAM=''; 
else;CAM=varargin{1};  carg=2; 
end
model=[]; 
if ischar(CAM)
elseif (isfield(CAM,'Stack')||isfield(CAM,'Elt'));
   model=CAM;  CAM=varargin{2}; carg=3; 
elseif isfield(CAM,'Y');
    if ischar(CAM.Y)&&strncmpi(CAM.Y,'win',3);CAM='Window';carg=1;
    else;CAM='Test';carg=1;
    end
end

[CAM,Cam]=comstr(CAM,1);
%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*NOSEM>
 
%---------------------------------------------------------------
%% #Get : get properties in a curve
%% #GetX Fail safe procedure to get the x axis values
if comstr(Cam,'getx')

 val=varargin{carg};carg=carg+1;
 if carg<=nargin; RunOpt=varargin{carg};carg=carg+1;
 else; RunOpt=[];
 end
 if ~isfield(val,'X')||~isfield(val,'Y')
    error('Not a valid curve'); 
 end

 t=val.X;if iscell(t);t=t{1};end; if isstruct(t);t=t.data;end
 if size(t,1)==1&&size(t,2)>1; t=t(:);end
 N=length(t); 
 if N~=size(val.Y,1) 
    error('x-axis (%i) is not consistent with .Y rows (%i)', ...
    N,size(val.Y,1)); 
 end
 out=t;

 if nargout>1&&comstr(Cam,'getxtime') 
 %% [t,f]=fe_curve('GetXTime',noise); Build the frequency vector
  if isfield(RunOpt,'BufTime')
   ind=((t-t(1))<RunOpt.BufTime);out=t(ind);N=length(out);
   RunOpt.SkipInd=0;
   if ~isempty(RunOpt.SkipTime)
       RunOpt.SkipInd=max(find(t<RunOpt.SkipTime));
   end
  end
  r1=diff(out(:,1));
  if std(r1)/mean(r1)>1e-8; 
   error('Can''t compute FFT for irregular time steps')
  end
  f=1/mean(r1)*(0:length(out)-1)'/length(out);
  out1=f; 
  if isstruct(RunOpt);RunOpt.N=length(out);out2=RunOpt;
  else;out2=length(out);end
  
 elseif nargout>1&&comstr(Cam,'getxfreq') % Build the time vector
  f=t; r1=diff(f); 
  if std(r1)/mean(r1)>1e-10; 
   error('Can''t compute FFT for uneven frequency steps')
  end
  if f(1)~=0; 
      i1=f(1)/r1(1);
      ind=round(f([1 end])/r1(1));N=ind(2)+5;
      ind=ind(1):ind(2);
      sdtw('_nb','Padding with zeroes as f(1)~=0');f=(0:N-1)*r1(1);
      t=(0:N-1)'/(N*diff(f(1:2))); out=t;out1=f;
      out2=struct('N',length(t),'Padding',ind);return;
  end
  N=length(f); 
  if abs(val.Y(N,1)./conj(val.Y(2,1))-1)>1e-10; % is conjugate spectrum
      N=2*N-2;
  end
  t=(0:N-1)'/(N*diff(f(1:2))); out=t;out1=f;out2=length(t);

 end

%% #getia indices for spectrum halves and real values
elseif comstr(Cam,'getia')

 N=varargin{carg};carg=carg+1;
 if isstruct(N); out=N;N=N.N;
 else;out=struct('ia',[],'ib',[],'ic',[],'ind',[]);
 end
 if round(rem(N,2))==1; 
    out.ia=N:-1:(N+3)/2;out.ib=2:(N+1)/2;out.ind=1:(N+1)/2;out.ic=[1];
 else; 
   out.ia=N:-1:N/2+2;out.ib=2:N/2;out.ind=1:N/2+1;out.ic=[1 N/2+1];
 end

%% #Get curve from a stack ---------------------------------------------------
elseif comstr(Cam,'get')

  st=Cam; [CAM,Cam]=comstr(CAM,4);
  if isempty(model)&&nargin>=carg ; model=varargin{carg};carg=carg+1; end

  % get all curves (by name)
  if nargin<carg && comstr(Cam,'curve') 
   out={};
   out=stack_get(model,'curve'); out=out(:,3)';

   Case=fe_case(model,'getcase');
   st={'DOFLoad','FVol','FSurf'};
   for j3=1:length(st)
    C1=stack_get(Case,st{j3});
    for j1=1:size(C1,1)
     if isfield(C1{j1,3},'curve')
        i2=[];
        for j2=1:length(C1{j1,3}.curve); 
         if ~ischar(C1{j1,3}.curve{j2});i2(end+1)=j2;end %#ok<AGROW>
        end
        out(end+[1:length(i2)])=C1{j1,3}.curve(i2);
     end
    end % j1
   end % j3
   return;
  end % return all curves


  if comstr(Cam,'curve') %
   %% #GetCurve
    
    curve_name=varargin{carg};carg=carg+1;
    if ischar(curve_name)
      if isempty(curve_name); out={}; return; end
      [out,i1]=stack_get(model,'curve',curve_name);
      if isempty(i1)&&comstr(lower(curve_name),'test'); out=[]; return; end % fe_curve test curves
      if isempty(i1) ;
       try
        Case=fe_case(model,'getcase'); st={'DOFLoad','FVol','FSurf'};
        for j3=1:length(st)
          C1=stack_get(Case,st{j3});
          for j1=1:size(C1,1)
           if isfield(C1{j1,3},'curve')
            i2=[];
            for j2=1:length(C1{j1,3}.curve); 
             if ~ischar(C1{j1,3}.curve{j2});i2(end+1)=j2;end %#ok<AGROW>
            end
            out=C1{j1,3}.curve(i2);return;
           end
          end % j1
        end % j3
       catch
        error('Unknown curve name'); 
       end
      end 
      out1=out{1,2}; out=out{1,3};
    elseif isa(curve_name,'cell')
      out={}; out1={}; 
      for j1=1:length(curve_name)
        if isempty(curve_name{j1}); r1=[];
        else;[r1,st1]=fe_curve('getcurve',model,curve_name{j1});
        end
        if ~isempty(r1); out{end+1}=r1; out1{end+1}=st1;return;end %#ok<AGROW>
      end
    elseif isstruct(curve_name) && isfield(curve_name,'name')
      out=fe_curve('getcurve',model,curve_name.name);
    else
      out=fe_curve(model,sprintf('GetId %i',curve_name)); 
    end
    
  elseif comstr(Cam,'id'); [CAM,Cam]=comstr(CAM,3); % 'GetId'
   if comstr(Cam,'new')
    %% #GetIdNew: get an ID not used by curves in model
    r1=stack_get(model,'curve');
    i1=cell2mat(cellfun(@(x)x.ID,r1(cellfun(@(x)isfield(x,'ID'),r1(:,3)),3),'uni',0));
    if isempty(i1); out=1; else; out=max(i1+1); end
    
   else
    %% #GetId: recover curve by ID
    i1=comstr(CAM,-1);
    if nargin>=carg
    elseif isempty(model); error('You must give a model'); end %xxx
    if isempty(i1); error('ID not provided'); end
    [out,i2]=stack_get(model,'curve');
    out=out(cellfun(@(x)isfield(x,'ID')&&any(x.ID==i1),out(:,3)),:);
    out1=out(:,2); out=out(:,3); 
    if length(out)==1; out=out{1}; out1=out1{1}; end
   
   end
  else;sdtw('''Get%s'' unknown',CAM);
  end
  
%---------------------------------------------------------------
%% #Check format consistency and reformat if needed
elseif comstr(Cam,'check');
 sdtw('_nb','Check command is obsolete, use FixCurve instead')
 out=fe_curve('FixCurve',varargin{:});
 return
 
curve=varargin{carg};carg=carg+1;

% SDT xf format 
if isfield(curve,'w')&&isfield(curve,'xf');
 out=feutil('rmfield',curve,{'w','xf','dof'});
 out.X=iigui(curve.w); out.Y=iigui(curve.xf); 
 if isfield(curve,'dof'); out.data=iigui(curve.dof);
  if size(out.data,2)>=8%Universal file 58 : Z axis 
    out.Z=out.data(:,8);
  end
 end
else;out=curve;
end


%% #Set ----------------------------------------------------------------------
elseif comstr(Cam,'set'); [CAM,Cam]=comstr(CAM,4);

  curve_name=varargin{carg};carg=carg+1;
  curve=varargin{carg};carg=carg+1;
  
  %% #SetLoadCurve : add curve in a specified load case
  if comstr(Cam,'loadcurve')
   sdtw('_nb','fe_curve setloadCurve is obsolete, you should now use fe_case SetCurve')
   % model=fe_curve(model,'set LoadCurve','Point load 1',2,'step 1e-4*10');
   if isnumeric(curve);ch=curve; data=varargin{carg};carg=carg+1;
   else; data=curve;ch=1;
   end
   out=fe_case(model,'SetCurve',curve_name,data,ch);
   %    Case=fe_case(model,'getcase');
   %    [r1,i1]=stack_get(Case,'',curve_name);
   %    r1{1,3}.curve{ch}=data;
   %    fprintf('Set %s:%s.curve{%i}\n',r1{1,1:2},ch);
   %    out = fe_case(model,r1{1,1},r1{1,2},r1{1,3});

  else % in add curve in model.Stack

   if ischar(curve)
     if comstr(lower(curve),'test')
       curve=fe_curve(curve);
     elseif comstr(lower(curve),'step') % Step Tmax
      warning('Use TestStep instead of step'); % XXX
      curve=fe_curve(['Test' curve]);
     else;curve=struct('ID',1,'X',[],'Y',curve,'Interp','');
     end
     if length(curve)==1; curve.name=curve_name; end 
   elseif isstruct(curve) && isfield(curve,'name') ...
                          && ~comstr(curve.name,curve_name)
     sdtw('_nb',sprintf('Curve name changed from %s to %s',...
                         curve.name,curve_name));
     curve.name=curve_name;
   elseif isstruct(curve) && ~isfield(curve,'name')
     sdtw('_nb',sprintf('Curve name set to %s',curve_name));
     curve.name=curve_name;
   end
   if length(curve)==1&&(~isfield(curve,'type')||isempty(curve.type))
     curve.type='fe_curve'; % curve will be editable in feplot
   end
   out=stack_set(model,'curve',curve_name,curve);

  end 
%% DataType These are datatypes defined in universal file format ----------
elseif comstr(Cam,'typedata')||comstr(Cam,'datatype');[CAM,Cam]=comstr(CAM,9);
  
 if comstr(Cam,'fix') % go from old field format to cell array
  r1=varargin{carg};carg=carg+1;
  if isfield(r1,'label');out={r1.label r1.unit r1.type};
  elseif iscell(r1);out=r1;
  elseif ~isstruct(r1); out=DataType(r1);
  else;out=r1;
  end
 else
  out=DataType(varargin{carg:end});
  if nargout==0
     fprintf('Available data types\n\n');
     out=[num2cell(out.value(:,1)) out.label(:)]';fprintf('%4i : %s\n',out{:});
     clear out;return;
  end
 end
 if ~isempty(strfind(Cam,'cell'))&&~iscell(out);
   out={out.label out.unit out.type};
 end
 if ~isempty(Cam)
  [CAM,Cam,st]=comstr('-label',[-25 4],CAM,Cam);
  if ~isempty(st); out{1}=st; end
 end
%---------------------------------------------------------------
%% #fun #.fun  Analysis types are also defined in UFF 55
elseif comstr(Cam,'typeana')||comstr(Cam,'analysistype')
  
  out={0','Unknown',{}
    1,'Static',{'AnimStatic','ScaleColorOne'} 
    2,'Normal Mode',{'AnimFreq','AnimColor'}
    3,'Cplx eig. order 1',{'AnimFreq','AnimColor'}
    4,'Transient',{'AnimTime','ScaleColorOne'} 
    5,'Frequency Response',{'AnimFreq','ScaleColorOne'}
    6,'Buckling',{'AnimFreq','AnimColor'}  
    7,'Cplx eig. order 2',{'AnimFreq','AnimColor'}};
 if carg<=nargin % Default handling for feplot
   def=varargin{carg};carg=carg+1;
   if isfield(def,'fun')&&isnumeric(def.fun)&&length(def.fun)>1;
       r1=def.fun(2);
   else; r1=0;
   end
   i1=find(vertcat(out{:,1})==r1);
   if length(i1)~=1;i1=1;end
   out=out(i1,:);
   if isfield(def,'def')&&~isreal(def.def);
           out{1,3}{end+1}='AnimFreq';
   end
 end
 
%% #TypeField ---------------------------------------------------------------
elseif comstr(Cam,'typefield')
    
%function st=FieldType(i1)
st1={'Unknown','General','Stress','Strain','Element Force','Temperature', ...
    'Heat Flux','Strain Energy','Displacement','Reaction Force', ...
    'Kinetic Energy','Velocity','Acceleration','Strain Energy Density', ...
    'Kinetic Energy Density','Hydro-Static Pressure','Heat Gradient', ...
    'Code Checking Value','Coefficient Of Pressure'};
%0 Unknown, 1 Scalar, 2: Tx Ty Tz, 3: Tx Ty Tz Rx Ry Rz,
%4: Sxx Sxy Syy Sxz Syz Szz, 5: Sxx Syx Szx Sxy Syy Szy Sxz Syz Szz

if nargin==1; out=struct('label',{st1},'value',0:length(st1)-1);
else;i1=varargin{2};end
if nargin==1
    if nargout==0;out=num2cell(out.value(:));out(:,2)=st1(:);
        out=out';fprintf('%4i  %s\n',out{:});clear out;
    end
elseif isnumeric(i1)
  if length(i1)>1; i1=0;  elseif i1+1>length(st1); i1=0;end
  out=st1{i1+1};
elseif ischar(i1)
  out=strmatch(comstr(i1,-27),comstr(st1,-27));
  if isempty(out); out=0; else; out=out(1);end
end

   
%% #DofLab - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'doflab');[CAM,Cam]=comstr(CAM,7);

if comstr(Cam,'load');
 st1={'Force','N',[13 0 1 0 0];
      'Moment','N/m',[100 -1 1 0 0]
      'd*S','m^3',[100 3 0 0 0]};
else
 st1={DataType('Displacement -cell');
      DataType('Rotation -cell')
      DataType('Pressure -cell')};
 for j1=1:size(st1,1);st1(j1,1:3)=st1{j1,1};end
end

out=num2cell([1:12 19]');
i1=[1 2 3 7 8 9];out(i1,2)=st1(1,1);out(i1,3)=st1(1,2);out(i1,4)=st1(1,3);
i1=[4 5 6 10 11 12];out(i1,2)=st1(2,1);out(i1,3)=st1(2,2);out(i1,4)=st1(2,3);
i1=13;out(i1,2)=st1(3,1);out(i1,3)=st1(3,2);out(i1,4)=st1(3,3);

if nargin>1; % Return dof labels if needed
 DOF=varargin{carg};carg=carg+1;
 if size(DOF,2)==1; 
  out1=fe_c(DOF);out1=out1(:);
  i2=[out{:,1}];nind=sparse(i2,1,1:length(i2),99,1);
  i2=nind(round(rem(DOF,1)*100));i3=find(i2);
  out1(i3,2:3)=out(i2(i3),3:4);
 elseif size(DOF,2)==5; % 5 column sensor format [Tag NodeId Dir[xyz]]
  try
   out1=sprintf('%20.2f_(%20i)',DOF(:,1:2)');
   out1=strrep(out1,'.00','   ');
   out1=reshape(out1,43,length(out1)/43)';
   out1(:,all(out1==' '))=''; out1(out1=='_')=' ';
  catch; out1=num2str(DOF(:,2));
  end
  out1=cellstr(out1);  out1(:,2)={out{1,3}}; out1(:,3)={out{1,4}};
 end
 out=out1;
end
 
%---------------------------------------------------------------
% #Plot (obsolete) Basic plotting capabilities
elseif comstr(Cam,'plot')

 if sp_util('issdt');
  if ~isempty(strfind(Cam,'-noii'))||~isempty(strfind(Cam,'-gca'))
   eval('out=ii_plp(varargin{:});');
  else % now display in iiplot
   ci=sdth.urn('iiplot');iicom(ci,'sub'); 
   iicom(ci,'curveinit',varargin{carg});
  end
   return;
 end  
 if ishandle(varargin{carg})
   ca=varargin{carg};carg=carg+1;
   if ~isempty(get(ca,'children'))
     delete(get(ca,'children'));
     axes(ca);set(ca,'visible','on','fontsize',7);
   end;
 else 
   figure;
   ca=gca;
 end
  
  curve=varargin{carg};carg=carg+1;
  ind_curve=[];
  
  if isfield(curve,'Stack')
    curve_name=varargin{carg};carg=carg+1;
    for i1=1:size(curve.Stack)
      if ischar(curve.Stack{i1,1}) && ischar(curve.Stack{i1,2})
        if comstr(curve.Stack{i1,1},'curve') && ...
            comstr(curve.Stack{i1,2},curve_name)
          ind_curve=[ind_curve i1]; %#ok<AGROW>
        end;
      end;
    end;
    
    if isempty(ind_curve); error('Unknown curve name'); end
    
    if length(ind_curve)>1
      warning('OpenFEM:curve',['More than one curve named ' curve_name]);
      disp('First curve retained');
      ind_curve=ind_curve(1);
    end;
    curve=curve.Stack{ind_curve,3};
  end
  
  if isfield(curve,'PlotFcn')
    if ischar(curve.PlotFcn); eval(curve.PlotFcn);
    else;error('Not expected');end
  else;line(curve.X,curve.Y);
  end
  
 curve=fe_curve('datatype fix',curve);
 try
     if isfield(curve,'xunit');
    st=curve.xunit{1,1}; 
    if ~strcmpi(curve.xunit{1,2},'none');st=[st ' ' curve.xunit{1,2}];end
    xlabel(st);
     end
 end
  
 try
 if isfield(curve,'yunit');
    st=curve.yunit{1,1}; 
    if ~strcmpi(curve.yunit{1,2},'none');st=[st ' ' curve.yunit{1,2}];end
    ylabel(st);
 end
 end
  
 if isfield(curve,'name'); title(curve.name); end
  
%% #TimeFreq Compute Time response to a signal given the transfer    
elseif comstr(Cam,'timefreq'); [CAM,Cam]=comstr(CAM,5);
  
  val=varargin{carg};carg=carg+1;
  TF=varargin{carg};carg=carg+1;
  
  % Fourier transform of input signal
  if isempty(val); 
   w=fe_curve('getx',TF); RunOpt=fe_curve('getia',length(w));  
   if abs(TF.Y(RunOpt.ia(1),1)./TF.Y(RunOpt.ib(1),1)-1)>1e-10
    N=2*length(w); w=[0:N-1]*mean(diff(w));
   end
   val=struct('X',w,'Y',ones(size(w)));
  else;w=fe_curve('getx',val);   N=length(w);   
  end
  RunOpt=fe_curve('getia',N);  
  if size(val.Y,1)==N; ffta=fft(val.Y,[],1);
  else;  ffta=fft(val.Y.',[],1);
  end;
  
  if isfield(TF,'xf');r1=TF.w; xf=TF.xf;
  elseif isfield(TF,'Y');r1=fe_curve('getx',TF);xf=TF.Y;
  else;r1=w; xf=TF;
  end

  if size(xf,1)==N/2; xf(N,1)=0;end % pad symmetry with zeros
  if size(xf,1)~=N;xf=xf.';end
  if any(size(xf)==N);
  elseif size(xf,1)==N/2;
  else
   error('Inconsistent Time Signal and Transfer');
  end
 
  step=mean(diff(w));  w=2*pi/step/N*[0:N-1];time=[0:(N-1)]'*step;
  
  
  Xt=zeros(size(xf,1),size(xf,2));
  
  for j1=1:size(xf,2)
    X=sum(ffta.*squeeze(xf(:,j1,:)),2);
    %-- Before IFFT of X padd with conjugate f(i+N/2+1)=conj(f(N/2-i+1))    
    X(RunOpt.ia)=conj(X(RunOpt.ib));  X(RunOpt.ic)=real(X(RunOpt.ic));
    Xt(:,j1)=real(ifft(X));
  end
  
  out=struct('X',{{time}},'Xlab',{{'Time'}},'Y',Xt, ...
    'name','Time-Freq response');
  try
   if isfield(TF,'X');
       out.X(2:length(TF.X))=TF.X(2:end);
       out.Xlab(2:length(TF.Xlab))=TF.Xlab(2:end);
   end
   if isfield(TF,'dof');out.dof=TF.dof;end
  end
  if nargout==0; fe_curve('plot',out);end
  
  
%---------------------------------------------------------------
%% #ResSpectrum : see "Engineering Vibration", Daniel J. Inman
%--                         Prentice Hall, 1994           
elseif comstr(Cam,'resspectrum'); [CAM,Cam]=comstr(CAM,12);
  
  Stack={'Flexibility','H= 1./(s2+2*zeta*omj*s+omj^2)'
         'Admittance','H= s./(s2+2*zeta*omj*s+omj^2)'
         'Accelerance','H=s2./(s2+2*zeta*omj*s+omj^2)'
         'Disp transmis.','H=[2*zeta*omj*s+omj^2]./(s2+2*zeta*omj*s+omj^2)';
         'RDisp transmis.','H=s2./(s2+2*zeta*omj*s+omj^2)'};
  % xxx stiffness on disp input
  if nargin<4
    Stack=Stack';fprintf('%20s : %s\n',Stack{:}); return;
  end;
  
  val=varargin{carg};carg=carg+1;
  freq=varargin{carg};carg=carg+1;
  if carg<=nargin; damp=varargin{carg};carg=carg+1;
  else; damp=[1e-3 1-2 1e-1];
  end
  damp=damp(:); freq=freq(:);
 
  
  % use closest actual frequency to avoid leakage
  freq=sort(freq); 
  [t,F,N]=fe_curve('getxtime',val);
  for j1=1:length(freq); 
   [r1,i1]=min(abs(F-freq(j1)));freq(j1)=F(i1);
  end
  
  RunOpt=fe_curve('getia',N);
  %F(find(abs(F)<eps))=eps;
  s=1i*2*pi*F(RunOpt.ind);s=s(:);s2=s.^2; U=fft(val.Y);U=U(RunOpt.ind);

  % derivatives of output signal
  RunOpt.Deriv=[]; RunOpt.tplot=0;
  i1=strfind(Cam,'-v'); 
  if ~isempty(i1);CAM(i1+[0:1])=''; RunOpt.Deriv=s; end
  Cam=comstr(CAM,-27);i1=strfind(Cam,'-a'); 
  if ~isempty(i1);CAM(i1+[0:1])=''; RunOpt.Deriv=s2; end
  Cam=comstr(CAM,-27);i1=strfind(Cam,'+v'); 
  if ~isempty(i1);CAM(i1+[0:1])=''; RunOpt.Deriv=1./s; end
  Cam=comstr(CAM,-27);i1=strfind(Cam,'+a'); 
  if ~isempty(i1);CAM(i1+[0:1])=''; RunOpt.Deriv=1./s2; end
  Cam=comstr(CAM,-27);i1=strfind(Cam,'-tplot');
  if ~isempty(i1);CAM(i1+[0:5])=''; RunOpt.tplot=1;end
  [CAM,Cam]=comstr(CAM,1);


  j1=strmatch(Cam(1:min(3,length(Cam))),comstr(Stack(:,1),-27));
  if isempty(j1);j1=1;end
  RunOpt.HFcn=horzcat(Stack{j1,2},';');

  out=struct('X',freq(:)*ones(1,length(damp)), ...
             'Z',ones(size(freq(:)))*damp(:)', ...
             'name','Stack{j1,1}', ...
             'Xlab',{fe_curve('datatypecell','Frequency')});
  z=[]; 
  for j1=1:length(freq)
  for j2=1:length(damp)
    omj=2*pi*freq(j1); zeta=damp(j2);
    eval(RunOpt.HFcn); Y=H.*U; 
    if ~isempty(RunOpt.Deriv);Y=Y.*RunOpt.Deriv;end 
    Y(RunOpt.ia)=conj(Y(RunOpt.ib)); 
    Y(RunOpt.ic)=real(Y(RunOpt.ic)); y=real(ifft(Y));
    if RunOpt.tplot; z(:,end+1)=y;end %#ok<AGROW>
    out.Y(j1,j2)=max(abs(y));
  end
  end
  % Debug for tplot
  if RunOpt.tplot; figure(100);clf; plot(val.X,z);iimouse; end

  if ~isempty(strfind(Cam,'pseudov'))  % Pseudo-velocity
    out.Y=out.Y.*freq(:,ones(1,size(out.Y,2)));
  elseif ~isempty(strfind(Cam,'pseudoa')) % Pseudo-acceleration
    f2=freq(:).^2; out.Y=out.Y.*f2(:,ones(1,size(out.Y,2)));
  end
   
%---------------------------------------------------------------
%% #Ifft out =fe_curve('IFFT',frames); continuous inverse transform
elseif comstr(Cam,'ifft');[CAM,Cam]=comstr(CAM,5);
  
  frames=varargin{carg};carg=carg+1;
  [t,f,N]=fe_curve('getxfreq',frames);
  RunOpt=fe_curve('getia',N); if isstruct(N);N=N.N;end 
  if comstr(Cam,'shock'); RunOpt.Coef=diff(t(1:2));[CAM,Cam]=comstr(CAM,6);
  else; RunOpt.Coef=1/N;
  end
  out=frames;
  if ~iscell(out.X);
   %out.X={t(:) 1:size(frames.Y,2)}; out.Xlab={'Time','channel'};
   out.X={t}; out.Xlab=fe_curve('datatype','time'); %{'Frequency'};
   for j1=1:length(size(frames.Y))-1
     out.X{end+1}=1:size(frames.Y,j1+1); 
     out.Xlab{end+1}=sprintf('channel%i',j1);
   end
  else; out.X{1}=t;
  end
  if ~isfield(out,'Xlab')||~isa(out.Xlab,'cell')
%   out.Xlab={'Time','channel'};
   out.Xlab={'Time'};
   for j1=1:length(frames.X)-1
     out.Xlab{end+1}=sprintf('channel%i',j1);
   end
  end
%   r2=frames.Y;out.Y(length(t),1)=0;
%   for j1=1:size(r2,3); 
%    X=squeeze(frames.Y(:,:,j1));
%    X(RunOpt.ia,:)=conj(X(RunOpt.ib,:));  X(RunOpt.ic,:)=real(X(RunOpt.ic,:));
%    out.Y(:,:,j1)=real(ifft(X))/RunOpt.Coef;% compute continuous IFT
%   end

%  out.Y(length(t),1)=0; ??? XXX juste la 1ere colonne ?
   if isfield(RunOpt,'Padding');
      X=zeros(N,size(frames.Y,2));X(RunOpt.Padding,:)=frames.Y;
   else;X=frames.Y;
   end
   X(RunOpt.ia,:,:,:)=conj(X(RunOpt.ib,:,:,:));
   X(RunOpt.ic,:,:,:)=real(X(RunOpt.ic,:,:,:));
   out.Y=real(ifft(X))/RunOpt.Coef;% compute continuous IFT

% #FFT out =fe_curve('FFT',frames); -----
elseif comstr(Cam,'fft');[CAM,Cam]=comstr(CAM,4);
  
  frames=varargin{carg};carg=carg+1;
  if carg<=nargin&&isstruct(varargin{carg});RunOpt=varargin{carg};carg=carg+1;
  else; [CAM,Cam,RunOpt]=StandardOptions(CAM,Cam);
  end
  [t,f,RunOpt]=fe_curve('getxtime',frames,RunOpt);  
  RunOpt=fe_curve('getia',RunOpt);  
  if comstr(Cam,'shock'); RunOpt.Coef=diff(t(1:2));[CAM,Cam]=comstr(CAM,6);
  else;RunOpt.Coef=1/RunOpt.N; end
  if comstr(Cam,'full'); RunOpt.full=1;[CAM,Cam]=comstr(CAM,5);
  else;RunOpt.full=0; end
  if ~isempty(Cam); fc=comstr(Cam,[-1 1]); 
  else;fc=f(end);
  end

  ind=ceil(length(f)*fc/f(end));  if ind>=length(f); ind=length(f) ;end
  ind=1:ind;

  out=frames;i1=size(frames.Y);i1(1)=length(ind);
  if ~iscell(out.X);
%    out.X={f(ind) 1:size(frames.Y,2)};
%    out.Xlab={'Frequency','channel'};
   out.X={f(ind)}; out.Xlab=fe_curve('datatypecell','freq');
   for j1=1:length(size(frames.Y))-1
     out.X{end+1}=1:size(frames.Y,j1+1); 
     out.Xlab{end+1}=sprintf('channel%i',j1);
   end
  else; out.X{1}=f(ind);
  end
  if ~isfield(out,'Xlab')||~isa(out.Xlab,'cell')
   %
   out.Xlab=fe_curve('datatypecell','freq');
   for j1=1:length(frames.X)-1
     out.Xlab{end+1}=sprintf('channel%i',j1);
   end
  else;out.Xlab{1}='Frequency';
  end
  
%   r2=frames.Y;
%   for j1=1:size(r2,3); 
%    SPEC=fft(squeeze(frames.Y(:,j1)))*RunOpt.Coef; % DFT or ShockDFT
%    out.Y(:,:,j1)=SPEC(ind);
%   end
  
  if size(frames.Y,1)>RunOpt.N % Multiple buffers
     in1=[1:RunOpt.N]+RunOpt.SkipInd;
     in1=RunOpt.SkipInd+1:fix(RunOpt.N*(1-RunOpt.Overlay)):size(frames.Y,1);
     in2=cell(length(in1));
     for j1=1:length(in1);in2{j1}=in1(j1)+[0:RunOpt.N-1];end
     in2=reshape(horzcat(in2{:}),RunOpt.N,[]);
     in1=size(frames.Y);in1(1)=[]; out.Y=zeros([length(ind) in1]);   
    for j1=1:prod(in1)
     SPEC=fft(reshape(squeeze(frames.Y(in2,j1)),size(in2)));
     SPEC=mean(SPEC,2)*RunOpt.Coef; % DFT or ShockDFT
     out.Y(:,j1)=SPEC(ind);
    end
  else % classical direct DFT or ShockDFT
    SPEC=fft(frames.Y,[],1)*RunOpt.Coef; 
    out.Y=SPEC(ind,:,:,:);
  end
  
  if ~RunOpt.full; % if not full return the half spectrum
   i1=unique([RunOpt.ib RunOpt.ic]);
   out.X{1}=out.X{1}(i1);  out.Y=out.Y(i1,:,:,:);
  end
  
%---------------------------------------------------------------
%% #h1 [H1,H2,Coh]=fe_curve('h1h2 ref',frames,window);
elseif comstr(Cam,'h1h2'); [CAM,Cam]=comstr(CAM,5);
  
  frames=varargin{carg};carg=carg+1;
  if carg>nargin; RO.Window=[];
  elseif isstruct(varargin{carg});
      RO=varargin{carg};carg=carg+1;
  else;RO.Window=varargin{carg};carg=carg+1;
  end
  if carg<=nargin; pond=varargin{carg};carg=carg+1;
  else;pond=[];
  end
  
  if ~isfield(RO,'Stack');
   [CAM,Cam,RO.Stack]=comstr('-stack',[-25 3],CAM,Cam);
  end
  if isfield(RO,'in');i1=RO.in;
  elseif isfield(RO,'lab_in');
      RO.in=find(ismember(frames{1}.X{2}(:,1),RO.lab_in));
      i1=RO.in;
  else;i1=comstr(Cam,[-1 1]); 
  end
  if ~iscell(frames); frames={frames};end
  
  i2=1:size(frames{1}.Y,2);i2(i1)=0;i2=find(i2);    
  
  t=frames{1}.X; if iscell(t); t=t{1};end
  RO.N=length(t); 
  if ~isfield(RO,'Window')||isempty(RO.Window); RO.Window='None'; end;
  if ischar(RO.Window); RO.Window=fe_curve(['window' RO.Window],RO.N);end
  win=RO.Window(:);
  
  dt=diff(t([1 end]))/(length(t)-1);
  f=1/dt*(0:length(t)-1)'/length(t);
  iw=1:round(450/1024*length(f));RO.tolf=f(2)*1e-3;
  if isfield(RO,'fmax');iw(f(iw)>RO.fmax+RO.tolf)=[];end
  if isfield(RO,'fmin');iw(f(iw)<RO.fmin-RO.tolf)=[];end
  
  %-- check validity for multiple frames case --%
  for j1=2:length(frames) 
    if length(frames{1}.X)~=length(frames{j1}.X);
      error('Inconsistent frames length');
    end;
  end;
  RO.NFrame=length(frames); RO.Dim3=0;
  if RO.NFrame==1; RO.NFrame=size(frames{1}.Y,3);RO.Dim3=1;end
  
%%  SIMO case
if length(i1)==1  
    k1=size(frames{1}.Y);
    iout=setdiff(1:k1(2),i1); % input channel
       
    Gyu=zeros(length(iw),k1(2)); Gyy=zeros(length(iw),length(iout));Guu=zeros(length(iw),1);
    
    %-- Using PSD averaging
    %-- instead of linear spectra
        
    for j2=1:RO.NFrame % fill in 
      if RO.Dim3; r2=frames{1}.Y(:,:,j2);else;r2=frames{j2}.Y;end
      r2 = r2.*win(:,ones(1,size(r2,2)));     % apply time window here
      SPEC=fft(r2);                      % compute FFT
      Gyu=Gyu+SPEC(iw,:).*conj(SPEC(iw,i1*ones(1,k1(2))));
      Gyy=Gyy+SPEC(iw,iout).*conj(SPEC(iw,iout));
      Guu=Guu+SPEC(iw,1).*conj(SPEC(iw,1));
    end
    r3=angle(Gyu(:,iout));
    Hlog=zeros(length(iw),length(iout));
    for j2=1:RO.NFrame % fill in 
      if RO.Dim3; r2=frames{1}.Y(:,:,j2);else;r2=frames{j2}.Y;end
      r2 = r2.*win(:,ones(1,size(r2,2)));     % apply time window here
      SPEC=fft(r2);                      % compute FFT
      % log
      Hlog=Hlog+log(SPEC(iw,iout)./SPEC(iw,i1*ones(1,length(iout))).*exp(-1i*r3));
    end
    Hlog=exp(Hlog/length(frames)).*exp(1i*r3);
    
    H1=Gyu(:,iout)./(Gyu(:,i1*ones(size(iout))));
    H2=conj(Gyy./Gyu(:,iout));
    COH=real(H1./H2);
    out=struct('X',f(iw),'H1',H1,'H2',H2,'COH',COH,'Gyy',Gyy/length(frames)/length(iw),...
      'Guu',Guu/length(frames)/length(iw),'Gyu',Gyu/length(frames)/length(iw),...
      'Hlog',Hlog);
    
%% MIMO Case    
else
    %-- set indices to extract matrix element from vector
    %-- compute only the lower triangular part of cross correlation matrix
    i1=i1(:)';i2=i2(:)';
    ind_comp=[i1 i2];
    ind_1=ind_comp(ones(length(ind_comp),1),:);
    ind_1=ind_1(:);
    ind_2=ind_comp(ones(length(ind_comp),1),:)';
    ind_2=ind_2(:);
    
    if isempty(pond); % default
        pond=norm(frames{1}.Y,Inf);
        for j1 = 2:length(frames); pond = pond + norm(frames{j1}.Y,inf);end
        pond=diag(pond/length(frames));
    elseif pond==0; pond=ones(size(frames{1}.Y,2),1); 
    else; pond=1./(pond+eps);  %-- to avoid singularities 
    end
    
    if size(pond,2)<2 || size(pond,1)<2; pond=diag(pond(:)); end;
    
    if size(pond,1)<size(frames{1}.Y,2) 
      error('Ponderation error - check for weighing vector');
    end
    
    
    %-- Using PSD averaging
    %-- instead of linear spectra
    
    Sy=zeros(size(frames{1}.Y));Gaa=zeros(length(iw),length(ind_comp).^2);
    for j1=1:length(frames)
      Y=(frames{j1}.Y*diag(1./diag(pond))).*...
        win(:,ones(size(frames{j1}.Y,2),1));
      Sy=fft(Y,[],1);
      Gaa=Gaa+Sy(iw,ind_1).*conj(Sy(iw,ind_2)); % SM 18/01/2005
    end;
    
    %-- compute estimators
    li1li2=length(i1)*length(i2);
    H1=zeros(length(iw),li1li2);
    Hv=zeros(length(iw),length(i1)*length(i2));
    co=[];
    for j1=1:length(iw)
      gaa=zeros(length(ind_comp));
      gaa=reshape(Gaa(j1,:),length(ind_comp),length(ind_comp)); 
      gyu=gaa(i2,i1);
      guu=gaa(i1,i1);
      gyy=gaa(i2,i2);
      %-- H1 estimator --%
      H1_=gyu*pinv(guu);%(guu\eye(length(i1)));
      H1(j1,:)=H1_(:)';
      %-- Hv, Hnu estimator 
      %-- Using pseudo inverse solution for least square problem
      %gi=(([guu gyu']*[guu gyu']')\eye(length(i1)));
      gi=pinv(([guu gyu']*[guu gyu']'));
      Hv_=([gyu gyy]*[guu gyu']')*gi;
      Hv(j1,:)=Hv_(:)';
      [u,s,v]=svd([guu gyu']*[guu gyu']');
      co(j1,:)=diag(s)';
    end;
    
    if size(pond,2)>=2 && size(pond,1)>=2; pond=diag(pond); end;
    
    Pond(1,:,:)=(pond(:,ones(size(pond,1),1))./pond(:,ones(size(pond,1),1))').^2;
    pond_vec=sqrt(Pond(1,i2,i1));
    pond_vec=pond_vec(:)';
    H1=H1.*pond_vec(ones(length(iw),1),:);
    Hv=Hv.*pond_vec(ones(length(iw),1),:);
    
    out=struct('X',f(iw),'H1',H1,'Hv',Hv,'Gyy',gyy/length(frames)/length(iw),...
      'Guu',guu/length(frames)/length(iw), ...
      'Gyu',gyu/length(frames)/length(iw),'Pond',pond_vec,'SingVal',co);
end
%iw=true(size(out.X));r1=min(diff(out.X))*1e-5;%Round off tolerance
%if isfield(RO,'fmax');iw(out.X>RO.fmax+r1)=false;end
%if isfield(RO,'fmin');iw(out.X<RO.fmin-r1)=false;end
%for st=intersect(fieldnames(out),{'X','H1','H2','Hv','Hlog','Gyy','Guu','Guy','Gyu','COH'})'
% r2=out.(st{1});
% if size(r2,1)==size(out.X,1);out.(st{1})=r2(iw,:); end
%end
if isfield(RO,'Out')
  RO.Stack=1; 
end
if RO.Stack % Format as SDT stack
 if ~isfield(RO,'Out');RO.Out={'H1','H2','Hv','Hlog','COH','Gyy','Guu','Gyu'};
 elseif ischar(RO.Out);RO.Out={RO.Out};
 end
 C1=out; st=fieldnames(C1);
 st(~ismember(lower(st),lower(RO.Out)))=[];
 out=repmat(st(:),1,3);out(:,1)={'curve'};
 % Build dof
 if ~isfield(RO,'dof')&&isfield(frames{1},'dof')
  try;
   dof=frames{1}.dof;
   iout=setdiff(1:size(frames{1}.dof,1),i1)'; % input channel
   odof=dof(iout,1); % Output dofs
   idof=dof(i1,1); % Input dofs
   odof=odof(:)*ones(1,size(idof,1));
   idof=ones(size(odof,1),1)*idof(:)';
   RO.dof=[odof(:) idof(:)];
  catch; sdtw('_nb','dof building failed');
  end
 end
 for j1=1:size(out,1)
  r2=C1.(st{j1});
  C2=struct('X',{{C1.X,(1:size(r2,2))'}},...
      'Xlab',{{'Frequency','DOF'}},'Y',r2,'Ylab',out{j1,2});
  if isfield(RO,'dof')&&any(strcmpi(st{j1},{'h1','h2','hv','hlog','coh'}))
   C2.dof=RO.dof;
  end
  if size(C2.Y,1)~=size(C1.X,1)
  elseif isfield(RO,'in');
    RO.lab_out=frames{1}.X{2};RO.lab_out(RO.in,:)=[];
    C2.X{3}=frames{1}.X{2}(RO.in,:); C2.Y=reshape(C2.Y,size(C2.Y,1),[],size(C2.X{3},1));
    C2.X{2}=RO.lab_out;
    C2.Xlab(2:3)={'Out','In'};
  end
  out{j1,3}=C2;
  if size(out,1)==1; out=out{3}; end
 end
 
end
  
%---------------------------------------------------------------
%% #spectra [Spec]=fe_curve('Spectra ref_channels',frames,nb_pts,window,overlap);
elseif comstr(Cam,'spectra'); [CAM,Cam]=comstr(CAM,8);
  
  clear frames
  frames=varargin{carg};carg=carg+1;
  nb_pts=varargin{carg};carg=carg+1;
  
  if carg<=nargin; win=varargin{carg};carg=carg+1; else;win=[]; end;
  if carg<=nargin; ovlp=varargin{carg};carg=carg+1; else;ovlp=1; end;
  if ovlp>=1 && ovlp<=100; ovlp=ovlp/100;
  elseif ovlp > 100; error('Overlap can not exceed 100%%');
  end; 
  if isstruct(frames); frames{1}=frames; end;
  
  i1=comstr(Cam,[-1 1]);
  i1=i1(1); 
  i2=1:size(frames{1}.Y,2);i2(i1)=0;i2=find(i2);
  if isempty(i1); error('No reference channel(s) specified'); end;
  
  
  li1=length(i1);li2=length(i2);li12=size(frames{1}.Y,2);

  Time=frames{1}.X; N=length(Time);
  if isempty(win); win='None'; end;
  if ischar(win); win=fe_curve(['window' win],nb_pts);end;
  win=win(:);
  f=1/diff(frames{1}.X(1:2))*[0:nb_pts-1]/nb_pts;
  iw=1:round(450/1024*length(f));
  win=win(:,ones(li12,1));
  Sy=zeros(li12,nb_pts);
  
  for j1=1:length(frames)
    Gyy{j1}=zeros(length(iw),li12^2);
    
    for k1=1:floor((length(Time)-nb_pts)/((1-ovlp)*nb_pts))+1   
      Y=frames{j1}.Y((1:nb_pts)+round((1-ovlp)*nb_pts*(k1-1)),:).*...
        win;
      Sy=Sy+fft(Y',[],2);
    end;
    
    for l1=1:length(iw)
      Gt=Sy(:,l1)*Sy(:,l1)';
      Gyy{j1}(l1,:)=Gt(:).'./nb_pts./k1;
    end;
  end;
  
  out=struct('X',f(iw),'Y',Gyy,'xlabel',...
    fe_curve('datatype',18),'ylabel',fe_curve('datatype',0),...
    'name','Power spectra');
%% #ReturnY ---------------------------------------------------------------
elseif comstr(Cam,'returny')

  [CAM,Cam,i1]=comstr('-extrapby0',[-25 3 1],CAM,Cam);
  RunOpt.Extrap='';if i1; RunOpt.Extrap='By0';end
%   [CAM,Cam,RunOpt.Typ]=comstr('-lin',[-25 3 1],CAM,Cam);
%   if RunOpt.Typ==0; RunOpt.Typ=2; end
  RunOpt.Interp='';
  if comstr(Cam,'-linear');RunOpt.Interp='linear'; [CAM,Cam]=comstr(CAM,8);end
  if comstr(Cam,'-lin'); RunOpt.Interp='linear'; [CAM,Cam]=comstr(CAM,4); end
  if comstr(Cam,'-log'); RunOpt.Interp='log'; [CAM,Cam]=comstr(CAM,4); end
  if comstr(Cam,'-stair'); RunOpt.Interp='stair'; [CAM,Cam]=comstr(CAM,6); end

  if nargin>=carg;  curve=varargin{carg};carg=carg+1; 
   if isfield(curve,'Stack'); model=curve;curve=varargin{carg};carg=carg+1;end
  end
  if nargin<carg; error('Bad number of input argument'); end
  x=varargin{carg};carg=carg+1;
  if nargin>=carg; model=varargin{carg};carg=carg+1; end
  
  if size(x,1)==1;x=x(:);end

  if isa(curve,'cell') % array of curves 
    out=cell(length(curve,1));
    for j1=1:length(curve)
      out{j1}=fe_curve(Cam,curve{j1},x,model);
    end
    out=reshape([out{:}],length(x),length(out));
    return; 
  end

  if isempty(curve); out=ones(size(x)); return; end

  if ischar(curve); % get curve by name or build it
      r2=fe_curve('getcurve',model,curve);
      if ~isempty(r2);curve=r2;
      else; fprintf('Rebuilding curve : %s\n',curve);
          curve=fe_curve(curve,x);
      end
  end
  curve_name=''; 
  if ~isfield(curve,'X')||~isfield(curve,'Y')
   if isstruct(curve)&&isfield(curve,'type')
     if isfield(curve,'name');curve_name=curve.name;end 
     curve=fe_curve('test',curve,x);
   else
    curve_name=varargin{carg};carg=carg+1;
    [row,i1]=stack_get(model,'curve',curve_name);
    if isempty(i1); 
     sdtw('_nb','Unknown curve name'); out=ones(size(x));return;
    end 
    curve=row{3};
   end
  elseif isfield(curve,'name');curve_name=curve.name;
  elseif isfield(curve,'Ylab')&&ischar(curve.Ylab);curve_name=curve.Ylab;
  end
  %% Deal with inter-/extrapolation
  if isempty(curve.X); curve.X=x; end
  if ischar(curve.Y)
    if comstr(comstr(curve.Y,-27),'step') % Step Tmax
      r1=comstr(curve.Y(5:length(curve.Y)),-1);
      out=zeros(size(x)); out(x<r1)=1; 
    elseif isempty(strfind(curve.Y,'out'));
        %out=eval(sprintf('fe_curve(''%s'',x)',curve.Y));
        if isfield(curve,'Edit')
         out=fe_curve(curve,x);%out=fe_curve(sprintf('%s',curve.Y),curve.Edit,x);
        else
         out=fe_curve(sprintf('%s',curve.Y),x);
        end
        if ischar(out.Y);error('Not expected');end
        out=out.Y;
    elseif strncmpi(curve.Y,'testeval',8);eval(curve.Y(9:end));
    else;eval(curve.Y);
    end
  elseif iscell(curve.Y); % Anonymous functions (see SDT-HBM)
    out=zeros(size(x)); 
  else
   if ~isempty(RunOpt.Extrap) % If extrapolation defined
   elseif isfield(curve,'Extrap')&&~isempty(curve.Extrap)
     RunOpt.Extrap=curve.Extrap;
   else; RunOpt.Extrap='By0'; % default if nothing said
   end
   if strcmpi(RunOpt.Extrap,'By0')||strcmpi(RunOpt.Extrap,'zero')
    % extrapolate if necessary 
    r1=curve.X; if iscell(r1); r1=r1{1};end; i1=find(x>=max(r1));
    if ~isempty(i1)
     if length(i1)==1&&abs(max(r1)/x(i1)-1)<1e-10 % just roundoff
         x(i1)=max(r1);
     else
      if size(curve.Y,1)==1&&size(curve.Y,2)==length(r1);curve.Y=curve.Y(:);
      end
      r1=[r1(:);x(i1(:))];if iscell(curve.X);curve.X{1}=r1;else;curve.X=r1;end
      curve.Y(length(r1),1)=0;
      sdtw('_nb','curve extrapolation %s',curve_name)
     end
    end
   end
   if sp_util('issdt')
    if isempty(RunOpt.Interp)&&(~isfield(curve,'Interp')||isempty(curve.Interp))
      curve.Interp='linear';
    elseif ~isempty(RunOpt.Interp) % interp as asked in returny command
      curve.Interp=RunOpt.Interp;
    end
    eval('out=feutilb(''computeR'',curve,x,0);');
   elseif any(ismember(x,curve.X))
    [i1,i2]=ismember(x,curve.X);
    out=zeros(size(x));out(i1)=curve.Y(i2(i1));
    if nnz(i1)~=length(x); warning('Non full coincidence');end

   else
    error('Curve interpolation supported with SDT');
   end
  end
  
  
%---------------------------------------------------------------
% #Read val=fe_curve('read h:/sdtdata/gefdyn/cas_test/input.cyt');
elseif comstr(Cam,'read'); [CAM,Cam]=comstr(CAM,5);
  
  fin=0;
  if carg<=nargin; fin=varargin{carg};CAM=fopen(fin);Cam=comstr(CAM,-27);end
  if isempty(Cam)
    [fname,wd]=uigetfile(CAM,'Pick a curve containing a file');
    if ~ischar(name); return; end
    [wd,fname,ext]=fileparts(fullfile(wd,fname));
  else
    [wd,fname,ext]=fileparts(CAM);
  end
  if isempty(wd); wd=pwd;end
  FileName=fullfile(wd,[fname ext]);
  
  if comstr(ext,'.cyt') % Cyberquake text format
    
    [fin,mes]=fopen(FileName); if ~isempty(mes); error(mes);end
    
    st=fgetl(fin); 
    if comstr(comstr(st,-27),'#cyber')
      st=fgetl(fin);
      %[i1,st]=comstr(st,'#','%i');
      %st=st(find(st=='=')+1:end);i1(2)=comstr(st,-1);
      r1=fscanf(fin,'%g');
      val=struct('name',fname,'type','fe_curve', ...
        'X',r1(1)+[0:r1(3)-1]'*r1(2),'Y',r1(5:end));
      %[1 2 r1(3) 0 0 0],'data',r1([1:2 5:end]));
      if length(val.Y)==2*r1(3)
        val.X=val.Y(1:2:end); val.Y=val.Y(2:2:end);
      end
    end
    
  elseif comstr(ext,'.cyb') % Cyberquake binary format
    
    [fin,mes]=fopen(FileName,'r','b'); if ~isempty(mes); error(mes);end
    fseek(fin,0,-1);  i1=fread(fin,1,'int8');
    val.name=char(fread(fin,i1,'char')');
    
    i1=fread(fin,2,'int32');
    i1=fread(fin,6,'float32');
    i1=fread(fin,2,'int32');
    
    fclose(fin);
 % Dactron sig format very initial reading - - - - - - - - - - - - -
  elseif comstr(ext,'.sig') 

   if fin<2; [fin,mes]=fopen(FileName,'r'); end
   fseek(fin,0,1);filelen=ftell(fin); fseek(fin,0,-1);
   val=struct('version',{{fread(fin,1,'int32') char(fread(fin,8,'char')')}});
   % ? ? Blocksize ? ?
   val.i1=fread(fin,8,'int32')';
   val.header=reshape(fread(fin,30,'char'),10,3)';
   st1=fread(fin,70,'*char');
   ind=unique([1 find(st1<26)]);
   for j1=1:length(ind)-1;
     val.st1{1,j1}=st1(ind(j1)+1:ind(j1+1)-1);
   end
   i2=filelen-val.i1(3)*4-4-ftell(fin);i2=i2/4;
   r1=fread(fin,i2,'int32')';fseek(fin,-i2*4,0);
   r1(2,:)=fread(fin,i2,'float')';
   val.r1=r1;
   
   fseek(fin,filelen-val.i1(3)*4-4,-1);
   val.Y=fread(fin,val.i1(3),'float');
   i1=fread(fin,'int32'); if i1~=0; warning('unexpected end');end
   fclose(fin);
   
  else;sdtw('''Read%s'' is not a valid call',ext);val=[];
  end
  
  if isempty(val)
  elseif carg<=nargin; mdl=varargin{carg};
    out=stack_set(mdl,'curve',fname,val);
  else 
    out=val;
  end
  
%% #Window -----------------------------------------------------------------
elseif comstr(Cam,'window'); [N,CAM,Cam]=comstr(CAM,'window','%i');
  
  % Symmetric windows
  if isempty(N)&&carg<=nargin; N=varargin{carg};carg=carg+1;end
  R1=[];
  if isempty(N); 
      [st,i1,i2,i3]=sscanf(CAM,'%s',1);N=sscanf(CAM(i3:end),'%i',1);
  elseif isfield(N,'Edit') % A window structure is passed
      R1=N.Edit;[CAM,Cam]=comstr(N.Y,7);
      if carg<=nargin; N=varargin{carg};carg=carg+1;
      else; N=1024;end
  end
  if length(N)>1; N=length(N);end % a time vector was passed
  [CAM,Cam,RunOpt.Period]=comstr('-per',[-25 3],CAM,Cam);
  [CAM,Cam,RunOpt.norm]=comstr('-norm',[-25 3],CAM,Cam);
  out={'None','';'BoxCar','pt0(1#%i#"first point") nSample(N#%i#"length")';...
        'Exponential','nS0(0#%i#"number of 0 points") nS1(10#%i#"number of 1 flat top points") alpha(10#%g#"last point at exp(-alpha)") nS0f(0#%i#"number of final 0 points")';...
        'FlatTop','';'KaiserBessel',''; 'Triangular',''; ...
        'Hanning','';'Hamming','';'Blackman',''};
   
  if isempty(Cam)||isequal(Cam,'s') % List implemented windows
    if nargout==0
     out=out'; fprintf('%20s : %s\n',out{:});clear out
    end
  elseif isempty(N); % Return editable curve
    i1=strcmpi(out(:,1),Cam);
    out=struct('X',[],'Y',['Window' CAM], ...
      'Edit',cingui('paramEdit',out{i1,2}));
  elseif comstr(Cam,'none');out=ones(N,1);
  elseif comstr(Cam,'boxcar'); [CAM,Cam]=comstr(CAM,7);   
    if isempty(R1);
        R1=cingui('paramEdit',out{strcmpi(out(:,1),'boxcar'),2},CAM);
    end
    R1=fe_def('cleanentry',R1);
    opt=comstr(CAM,[-1 R1.pt0 R1.nSample]);pt0=opt(1); nSample=opt(2);
    out=[zeros(opt(1)-1,1);ones(opt(2),1);zeros(N-nSample-pt0+1,1)]; 
  elseif comstr(Cam,'exponential');[CAM,Cam]=comstr(CAM,12);
    if isempty(R1);
       R1=out{strcmpi(out(:,1),'exponential'),2};
       R1=cingui('paramEdit',R1,CAM);
    end
    R1=fe_def('cleanentry',R1);
    opt=comstr(CAM,[-1 R1.nS0 R1.nS1 R1.alpha R1.nS0f]);
    opt(1:2)=round(opt(1:2));  if length(opt)<4; opt(4)=0; end
    nS0=opt(1); nS1=opt(2); alpha=opt(3); nS0f=opt(4);
    Ne=N-nS0-nS1-nS0f-1;  
    if Ne<0; error('Improper number of zero and flat top points');end
    out=[zeros(nS0,1);ones(nS1,1); exp(-(0:Ne)'/Ne*alpha); zeros(nS0f,1)];
  elseif comstr(Cam,'flattop'); [CAM,Cam]=comstr(CAM,8); 
    out=1-1.93*cos(2*pi*((1:N)'-0.5)/N)+1.29*cos(4*pi*((1:N)'-0.5)/N)-...
        0.388*cos(6*pi*((1:N)'-0.5)/N)+0.0322*cos(8*pi*((1:N)'-0.5)/N);
  elseif comstr(Cam,'kaiserbessel');[CAM,Cam]=comstr(CAM,13);
    k=(1:N)';
    out=1-1.24*cos(2*pi*(k-0.5)/N)+0.244*cos(4*pi*(k-0.5)/N)-...
        0.00305*cos(6*pi*(k-0.5)/N);
  elseif comstr(Cam,'triangular');[CAM,Cam]=comstr(CAM,11);
    k=(1:N)';
    out=(N-1)/2-abs(k-0.5-N/2);
  elseif comstr(Cam,'han') % Hanning window
    out = (1 - cos(2*pi/(N-RunOpt.Period)*(0:N-1)'))/2; % don't finish window
  elseif comstr(Cam,'ham') % Hamming window
    out = .54 - .46*cos(2*pi/(N-RunOpt.Period)*(0:N-1)');
  elseif comstr(Cam,'bl') % Blackmann
    t=(0:N-1)'/(N-RunOpt.Period);
    out=0.42-0.5*cos(2*pi*t)+0.08*cos(4*pi*t); 
  else
    sdtw('_err','Unknown fe_curve windows: %s',Cam)  
  end
  
  if RunOpt.norm % Normalize as fct_101
    out=out*N/sum(out);  
  end

%% #Butter Order FreqHz -----------------------------------------------------
elseif comstr(Cam,'butter'); [opt,CAM,Cam]=comstr(CAM,'butter','%g');
  
if rem(opt(1),2)==1; opt(1)=opt(1)+1; end % only even orders are supported

% opt(2) = freq_filt/(F_ech/2);
%continuous poles
r1 = exp(1i*((1:2:opt(1)-1)/opt(1) + 1)*pi/2);r1=[r1;conj(r1)];
% low pass state-space
a=[];b=zeros(0,1);c=zeros(1,0); d=1; b1 = [1; 0]; c1 = [0 1]; 
for j1=1:size(r1,2)
    den = real(poly(r1(:,j1))); a1 = [-den(2) -den(3); 1 0];d1=0;
    a = [a zeros(size(a,1),2); b1*c a1];
    b = [b; b1*d]; c = [d1*c c1];d=d1*d;
end
a=a*opt(2)*2*pi;b=b*opt(2)*2*pi;
if nargout<2;out=struct('a',a,'b',b,'c',c,'d',d);
else; % continuous polynomial
  out1=poly(a); % den
  out=[zeros(1,opt(1)) polyval(out1,0)];
end

  
%% BandPass -- Perfect band pass filter by FFT/IFFT --%  
elseif comstr(Cam,'bandpass'); [CAM,Cam]=comstr(CAM,9);
  
  
  if ~isa(varargin{carg},'cell'); frames={varargin{carg}};carg=carg+1;
  else;frames=varargin{carg};
  end
  if carg<=nargin&&isfield(varargin{carg},'fmin')
    RunOpt=varargin{carg};
  else
   RunOpt.unit=1;
   if comstr(Cam,'hz'); [CAM,Cam]=comstr(CAM,3);
   elseif comstr(Cam,'rd'); [CAM,Cam]=comstr(CAM,3); RunOpt.unit=1/2/pi;
   end
   r1=str2num(Cam)*RunOpt.unit;  %#ok<ST2NM>
   RunOpt.fmin=min(r1);RunOpt.fmax=max(r1);
  end
  out=cell(1,length(frames));
  for j1=1:length(frames)
    r1=frames{j1}; 
    [t,f,RunOpt.N]=fe_curve('getxtime',r1);
    RunOpt=fe_curve('getia',RunOpt);  
    if isfield(r1,'X');r2=r1;
    else;r2=struct('X',{{t}},'Xlab',{{'Time'}},'Y',[]);
    end
    if iscell(r2.X);r2.X{1}=t;else;r2.X=t;end  
    
    if find(size(r1.Y)==RunOpt.N)==1; Y=fft(r1.Y);
    else;Y=fft((r1.Y).');
    end;
    Y(f < RunOpt.fmin | f > RunOpt.fmax,:)=0;
    Y(RunOpt.ia,:)=conj(Y(RunOpt.ib,:)); 
    Y(RunOpt.ic,:)=real(Y(RunOpt.ic,:));
    r2.Y=real(ifft(Y));
    out{j1}=r2;
  end
  
  if length(out)==1; out=out{1}; end;

%---------------------------------------------------------------
%-- fe_curve('zoomfft (CenterFreq) (ZoomFactor)',frames)
elseif comstr(Cam,'zoomfft'); [CAM,Cam]=comstr(CAM,8);

  if ~isa(varargin{carg},'cell'); frames={varargin{carg}};carg=carg+1;
  else;frames=varargin{carg};
  end;    

  r1=comstr(CAM,-1);
  if length(r1)~=2; 
   error('You must provide center frequency and zoom factor');
  end
  RunOpt.Zoom=round(ftarg(2));
  RunOpt.f0=ftarg(1);

  for i1=1:length(frames)
    
    t=frames{i1}.X(:);
    r1=frames{i1}.Y; if find(size(r1)==length(t))==2;r1=r1';end
    r1=fft(r1.*exp(-sqrt(-1)*2*pi*RunOpt.f0*t)); % shift and transform

     N=length(t)/RunOpt.Zoom; r1(N:size(r1,1)-N)=0;  % filter  
     r1(end:-1:end-N+1)=conj(r1(2:N+1));          % make symmetric
     r1(1)=real(r1(1));
 
     x=real(ifft(r1));
     x=x(1:RunOpt.Zoom:end);tx=t(1:RunOpt.Zoom:end);
     fx=[1/diff(tx(1:2))*[0:length(tx)-1]'/length(tx)]+RunOpt.f0;
     r1=fft(x)*RunOpt.Zoom;
 
     ind=round(length(fx)*[.05 .45]);ind=ind(1):ind(2);  
     out=struct('X',fx(ind),'Y',r1(ind,:), ...
        'Xlab',{fe_curve('datatypecell','Freq')});
    
  end
 
  
%% # ---Noise
elseif comstr(Cam,'noise'); [CAM,Cam]=comstr(CAM,6);
 if comstr(CAM,'cn'); [CAM,Cam]=comstr(CAM,3);
  %% Noise Complex Normal law
  if carg<=nargin; XF=varargin{carg}; carg=carg+1;
  else; error('Curve missing as second argument'); 
  end
  if carg<=nargin; RO=varargin{carg}; carg=carg+1;
  else; error('RO containing var or rel or extra options missing as third argument'); 
  end
  if isfield(XF,'xf'); xf=XF.xf;
  elseif  isfield(XF,'Y'); xf=XF.Y;   
  else; error('Unknow curve type');
  end
  r1=(randn(size(xf))+1i*randn(size(xf))); % Standard complex normal law
  if isfield(RO,'var')&&~isempty(RO.var) % Add noise with constant variance
   xf=xf+sqrt(RO.var/2)*r1;
  elseif isfield(RO,'rel')&&~isempty(RO.rel) % Add relative noise
   xf=xf+RO.rel*xf.*r1;
  end
  
  if isfield(XF,'xf'); XF.xf=xf;
  elseif  isfield(XF,'Y'); XF.Y=xf;   
  end
  
  out=XF;  
 else
  %% Noise Unit & given PSD noise generation from randomly phased cosines
  M=varargin{carg};carg=carg+1;
  if length(M)<2;M(2)=1;end
  RunOpt.Even=0; RunOpt.N=M(1);
  if mod(M(1),2)==0; M(1)=M(1)/2; RunOpt.Even=1; else;M(1)=(M(1)+1)/2; end
  
  fe=varargin{carg};carg=carg+1;%Sampling
  RunOpt.f_max=[]; RunOpt.Spec=[];
  
  if nargin==4
    if length(varargin{4})==1
      RunOpt.f_max=varargin{carg};carg=carg+1;
      if RunOpt.f_max>fe/2;error('Select a max frequency <= sampling rate / 2');end
    elseif isa(varargin{4},'double')
      RunOpt.Spec=varargin{4};
      if length(RunOpt.Spec)~=M(1)
        warning('resampling PSD spectrum...');
        if length(RunOpt.Spec)<M(1)
          RunOpt.Spec=interp1((0:length(RunOpt.Spec)-1)/(length(RunOpt.Spec)-1),...
            RunOpt.Spec,(0:M(1)-1)/(M(1)-1));
        else
          M_ech=lcm(M(1),length(RunOpt.Spec));
          RunOpt.Spec=interp1(round(linspace(1,M_ech,length(RunOpt.Spec))),...
            RunOpt.Spec,1:M_ech);
          RunOpt.Spec=RunOpt.Spec(round(linspace(1,M_ech,M(1))));
        end;
      end;
    else
      error('Bad argument type');
    end;
  end;
  
  if isempty(RunOpt.f_max) && isempty(RunOpt.Spec)
    X=ones(M).*exp(1i*rand(M)*2*pi);
  elseif ~isempty(RunOpt.f_max) && isempty(RunOpt.Spec)
    X=[ones(ceil(M(1)*(2*RunOpt.f_max/fe)),M(2)) ; ...
        zeros(M(1)-length(ones(ceil(M(1)*(2*RunOpt.f_max/fe)),1)),M(2))].*...
      exp(1i*rand(M)*2*pi);
  elseif isempty(RunOpt.f_max) && ~isempty(RunOpt.Spec)
    X=sqrt(RunOpt.Spec(:)).*exp(1i*rand(M)*2*pi);
  end;
  
  t=(0:RunOpt.N-1)'/fe;r2=fe_curve('getia',length(t)); 
  X(r2.ia,:,:,:)=conj(X(r2.ib,:,:,:));
  X(r2.ic,:,:,:)=real(X(r2.ic,:,:,:));
  out=struct('X',{{t (1:M(2))'}}, ...
      'Xlab',{{'Time','Channel'}},'Y',2*(fe)^(1/2)*real(ifft(X)));
  
  for j1=1:size(out.Y,2);r1=out.Y(:,j1);r1=r1-mean(r1);out.Y(:,1)=r1;end
 end
 out=fe_curve('fixcurve',out);
%% #Test ---------------------------------------------------------------
elseif comstr(Cam,'test'); [CAM,Cam]=comstr(CAM,5);

% Find ID if given
%% #TestDealWithInput
[CAM,Cam,RunOpt.ID]=comstr('-id',[-25 1],CAM,Cam);
R1=[];opt=[];if isempty(RunOpt.ID); RunOpt.ID=1;end
RunOpt.Interp=0; RunOpt.t=[];
if isempty(CAM)&&carg<nargin; 
    opt=varargin{carg};carg=carg+1;
    if isfield(opt,'type')&&~strcmpi(opt.type,'fe_curve');
       %fe_curve('Test',struct('type',R2.window),f)
     [CAM,Cam]=comstr(opt.type,1);
     if isfield(opt,'Edit');R1=opt.Edit;
     elseif isequal(fieldnames(opt),{'type'})
       opt=[];
     else; R1=opt;
     end
    elseif isfield(opt,'Y')&&ischar(opt.Y);[CAM,Cam]=comstr(opt.Y,5);
     if isfield(opt,'Edit');R1=opt.Edit;end
    end
    if carg<=nargin&&isnumeric(varargin{carg});
        RunOpt.t=varargin{carg};carg=carg+1;
    end
else;
 if nargin>=carg&&isstruct(varargin{carg}) % Edit field given
  RunOpt.RO=varargin{carg}; carg=carg+1;
  % fe_curve('test',struct('type','sine','f',10))
  if isempty(Cam)&&isfield(R1,'type');
      [CAM,Cam]=comstr(R1.type,1);
  end
  if carg<=nargin&&isnumeric(varargin{carg});
        RunOpt.t=varargin{carg};carg=carg+1;
  end
  [opt,R1,RunOpt]=TestList(Cam,RunOpt); 
 else
  RunOpt.NeedOut=1; if carg<=nargin; RunOpt.t=varargin{carg};carg=carg+1;end
  if strcmpi(Cam,'list'); out=TestList;return;end
  [opt,R1,RunOpt]=TestList(Cam,RunOpt); 
  if ~isempty(RunOpt)&&(strcmpi(RunOpt.CAM,'frame')||strcmpi(RunOpt.CAM,'acq'))
  elseif ~isfield(RunOpt,'t')||isempty(RunOpt.t);out=opt;return;% Return editable structure
  end 
  CAM=RunOpt.CAM;Cam=lower(CAM);
 end
end
    
if isempty(RunOpt.t)&&~isfield(opt,'Y')
   if nargout==0; fe_curve('list');
   else; out=fe_curve(list);end
   return
   
elseif comstr(Cam,'frame');   
%% #TestFrame frame - - - - - - - - - - - - - - - - - - - - - - - - - -
% Example where to test frames are built
  fs=512;  %-- sampling frequency
  ech_length=4;  %-- sample length (s)
  noise=fe_curve('Noise',fs*ech_length,fs);  %-- computes noise
  
  %-- build the curve associated to the time signal of noise
  if ~iscell(noise.X); noise.X={noise.X};end;X=noise.X{1};
  out{1}=struct('X',{noise.X},'Y',noise.Y,'Xlab',...
    {fe_curve('DataTypeCell','Time')},'Ylab',...
    {fe_curve('DataTypeCell','Excit. force')},'name','Input at DOF 2');  
  
  %-- Setup a 3DOF oscillator 
  
  Puls = [30 80 150]'*2*pi;  %-- natural frequencies
  Damp = [.02 .015 .01]';  %-- damping
  Amp = [1 2 -1;2 -1 1;-1 1 2];  %-- pseudo "mode shapes"
  Amp=Amp./det(Amp);
  
  C=[1 0 0];  %-- Observation matrix
  B=[0 1 0]';  %-- Command matrix
    
  freq=([0:length(X)-1]/length(X))*fs*2*pi;  % Compute frequency vector %
  %-- Eliminate frequencies corresponding to the aliased part of the noise
  %-- spectrum
  
  freq=freq(1:length(X)/2)';
  FRF=nor2xf(Puls,Damp,Amp*B,C*Amp,freq);%-- Compute transfer function 
  
  %-- Compute the time response to input noise
  Resp=fe_curve('TimeFreq',noise,struct('w',freq,'xf',FRF));
   
  %-- build the curve associated to the time signal of response
  
  out{2}=struct('X',Resp.X,'Y',Resp.Y,'Xlab',{{fe_curve('DataTypeCell','Time')}},...
    'Ylab',{fe_curve('DataTypeCell','Displacement')},'name','Output at DOF 1');
  
  if nargout==0
    h1=figure(1);fe_curve('Plot',h1,out{1});    
    figure(2);
    semilogy(freq/2/pi,abs(FRF));
    title('FRF DOF 2 / DOF 1')
    xlabel('Frequency Hz');
    ylabel('Transfert m/N');
    h2=figure(3);
    fe_curve('Plot',h2,out{2});
  end;
  return

%% #TestAcq : simulated acquisition
elseif comstr(Cam,'acq'); 
    
 out=fe_curve('testframe');    % 3 DOF system response
 r1=struct('X',{{out{1}.X{1} [out{1}.Ylab;out{2}.Ylab]}},...
  'Y',[out{1}.Y out{2}.Y],...
  'Xlab',{[out{2}.Xlab 'Channel']},...
  'Ylab',2,...
  'name','Time frames');
 out=r1; return
end

%% #GenericEDIT now the generic handling of options - - - - - - - - - - - - - -
if carg<=nargin;error('Expecting RunOpt.t');t=varargin{carg};t=t(:);   else; t=[];  end
[CAM,Cam,Tf]=comstr('-stoptime',[-25 2],CAM,Cam);
st=Cam; st(ismember(st,'0123456789.,'))=' ';
[st,i1,i2,i3]=sscanf(st,'%s',1);RunOpt.EmptyOpt=0;
if isstruct(opt);
elseif 1==2 % obsolete handled in TestList
 %sp_util('issdt'); % SDT gui handling
 r2=TestList(st);
 if isstruct(opt) % curve is provided
     out=opt;opt=[];R1=fe_def('cleanentry',R1);
 elseif ~isempty(r2); % known in list 
     [opt,Cam]=comstr(CAM,i3);
     if ~isempty(opt)&&all(ismember(lower(opt), ... % Just numbers : compat
         ['0123456789.ed+-: *' char(1:32)]));
         opt=comstr(opt,-1);
     else;opt=[];
     end
     if isempty(opt) % Not just numbers : current SDT
       if ~isstruct(R1);st={};
        if comstr(Cam,'eval')
           R2=cingui('paramedit',r2,['st"' comstr(CAM,i3) '"']);
        else; R2=cingui('paramedit',r2,comstr(CAM,i3)); % SDT style
         if length(CAM)>i3; RunOpt.EmptyOpt=0;end
        end
       else % Possibly reuse edit too
        if comstr(Cam,'eval')
           R2=cingui('paramedit-given',r2,['st"' comstr(CAM,i3) '"']);
        else; R2=cingui('paramedit-given',r2,comstr(CAM,i3)); % SDT style
         if length(CAM)>i3; RunOpt.EmptyOpt=0;end
        end
        st=setdiff(fieldnames(R1),[{'st'};fieldnames(R2)]);
        for j1=1:length(st); % propag values in R1 fe_curve('TestNoise',Edit,t)
         r3=R1.(st{j1});
         if isjava(r3)||isfield(r3,'value') % From .Edit
          R2.(st{j1})=r3;
         elseif isa(r3,'function_handle');
          R2.(st{j1}).type='function_handle';R2.(st{j1}).value=r3;
         else; % from data structure
          R2.(st{j1}).value=r3;
         end
        end 
       end
       R1=R2;
     else; 
       RunOpt.Edit=cingui('paramedit',r2); R1=[];% OpenFEM style with numeric
     end
 else % OpenFEM style with only numbers given
  opt=comstr(CAM,-1); R1=[];% Compat layer expected in each subfun
 end
 Cam=lower(CAM);
 if isstruct(R1);st=fieldnames(R1);st(strcmp(st,'st'))=[];else;st={};end
 
else;% OpenFEM only handling  
    opt=comstr(comstr(CAM,length(st)+1),-1);
end
RunOpt.SetVal=0;t=RunOpt.t(:);
if isstruct(R1);RunOpt.Edit=R1;
elseif ~isfield(RunOpt,'Edit');RunOpt.Edit=[];
end

%Begin switch on specific options for various curves - - - - - - - - - - -
%% TestList - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if comstr(Cam,'list'); [CAM,Cam]=comstr(CAM,5);
 % just provide a list of test functions
 if nargout==0; TestList;
 else; out=TestList; out=out(:,1);
 end
 
elseif comstr(Cam,'ramp');[CAM,Cam]=comstr(CAM,5);
 %% #TestRamp - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if isstruct(R1);R1=fe_def('cleanentry',R1);
 elseif length(opt)==2; 
    R1=struct('NStep',opt(1),'Yf',opt(2));RunOpt.SetVal=1;
 elseif length(opt)==3; 
    R1=struct('t0',opt(1),'t1',opt(2),'Yf',opt(3));RunOpt.SetVal=1;
 end
 if isfield(R1,'NStep')&&~isempty(R1.NStep)&&isfinite(R1.NStep) % NStep FinalValue :  old call
%  sdtw('_nb','TestRamp curve should be defined giving gradient')
  if isempty(t); t=(1:R1.NStep)'; end
  r1=1:R1.NStep;r1=r1(:);
  out=struct('X',t,'Y',r1/r1(end)*R1.Yf,'name','ramp',...
            'PlotInfo',{{'sub','1 1'; 'show','real'; 'scale','xlin;ylin'}});
 elseif ~isempty(R1.Yf) % t0 t1 X1   tabulated :
   if R1.t0>R1.t1;error('TestRamp : t0 should <t1');end
   out=struct('X',[0 R1.t0 R1.t1 2*R1.t1]','Y',[0 0 R1.Yf R1.Yf]', ...
      'name','ramp','Interp','linear','Extrap','flat','Edit',[]);
  if R1.t0==0; out.X=out.X(2:end);out.Y=out.Y(2:end); end
 else
   sdtw('_err','obsolete call. testramp should be called with t0 t1 x1')
 end
 RunOpt.Interp=1; out.name='ramp';

%% #TestStep Duration - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'step');[CAM,Cam]=comstr(CAM,5);

 if isstruct(R1);R1=fe_def('cleanentry',R1);
 else; R1=struct('t1',opt); RunOpt.SetVal=1;
 end
 if isempty(opt); opt=1; end % default
 out=struct('X',[0 R1.t1]','Y',[1 0]','name','step','Interp','stair',...
            'Extrap','flat','Edit',[]);
 RunOpt.Interp=1;

%% #TestSweep - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% TestSweep freq_min freq_max T_start T_stop % OpenFEM call
% TestSweep fmin=freq_min fmax=freq_max t0=T_start t1=T_stop % SDT call
% OBSOLETE NStep Tf freq_min, freq_max T_start T_stop 
% sdtweb fe_curve testlist % for parameters
elseif comstr(Cam,'sweep');[CAM,Cam]=comstr(CAM,6);

 if isstruct(R1);R1=fe_def('cleanentry',R1);
 else %if length(opt)<4;
  if length(opt)>4; t=linspace(0,opt(2),opt(1))'; %OBSOLETE t in Cam
  elseif length(opt)<2;error('You must provide f_min, f_max, t0,t1');
  else; opt=[0 0 opt]; % NStep Tf not given
  end   
  R1=struct('fmin',opt(3),'fmax',opt(4),'t0',opt(5),'t1',opt(6));
  RunOpt.SetVal=1;
 end
 
 % if length(opt)<5; opt(5)=0;end;if length(opt)<6; opt(6)=opt(2);end
 if ~isfield(R1,'fcn');R1.fcn='cos';end
 RunOpt.Edit.fcn=struct('type','string','value',char(R1.fcn)); 
 if ~isempty(t)
  if any(strcmpi(char(R1.fcn),{'cos','sin'}));r3=linsweep(t,R1);
  else; r3=feval(R1.fcn,t,R1);
  end
  if isfield(R1,'A');r3=R1.A*r3;end
  % Future need : exponential sweep : f(t)=f0 * k^t
  % x(t)=sin(2pi*f0(k^t-1/ln(k)+f0)
  out=struct('X',t,'Y',r3);
  out.name=sprintf('sweep %g-%g Hz',[R1.fmin R1.fmax]);
  out.type='tabulated';
 else
  out=struct('X',[],'Y',['TestSweep ' CAM],'Edit',[]);
 end
 
%% #TestRicker : SignalDuration Amplitude - - - - - - - - - - - - - - - - - - 
% OBSOLETE TestRicker SignalDuration Nstep Amplitude TotalTime
% TestRicker SignalDuration Amplitude % OpenFEM call (xxx is this still working ?)
% TestRicker dt=SignalDuration A=Amplitude t0=Beginning_time % SDT call 
elseif comstr(Cam,'ricker');
 
 if isstruct(R1);R1=fe_def('cleanentry',R1);
 else %if length(opt)<4;
     if length(opt)<3;opt(3)=0;end
     if length(opt)==4 % this is rather obsolete, ask for time vect in command
      % dt Nsteps A timstep
      sdtw('_nb','You should call TestRicker with a time vector')
      t=(0:opt(4):opt(2)*opt(4))';
      opt=opt([1 3]); opt(3)=0;
     end
     R1=struct('dt',opt(1),'A',opt(2),'t0',opt(3));
     RunOpt.SetVal=1;
 end

 if ~isempty(t)
  ts=R1.dt/2; tp=R1.dt/2;
  a=t-R1.t0; % shift time
  a = (a-ts).*(a-ts)/(tp*tp);
  y = R1.A*(1 - exp(1)*(a .* exp(-a))); y(t<=R1.t0|t>2*ts+R1.t0)=0.;
  out=struct('ID',RunOpt.ID,'X',t(:),'Xlab',{{'Time'}},'Y',y(:), ...
      'name','ricker');
  out.type='tabulated';
 else
  out=struct('ID',RunOpt.ID,'X',t(:),'Y','TestRicker','name','ricker','Edit',[]);
 end
 
%% #TestCosHan F0,N0,A - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% TestCosHan F0 N0 A % OpenFEM call
% TestCosHan f0=F0 n0=N0 A=A % SDT call
elseif comstr(Cam,'coshan');[CAM,Cam]=comstr(CAM,7);
 
 % Get R1:
 if isstruct(R1);R1=fe_def('cleanentry',R1);
 else % old OpenFEM call
  if length(opt)~=3; 
   error('TestCosHan must have 3 input parameters : F0 N0 A'); 
  else
   R1=struct('f0',opt(1),'n0',opt(2),'A',opt(3)); %  F0 N0 A
  end
  RunOpt.SetVal=1;
 end

 % Build curve:
 if isempty(t) % implicit definition
  out=struct('ID',RunOpt.ID,'X',[],'Y','TestCosHan','name','CosHan','Edit',[]);
 else % from time vector
  t=t(:);
  out=struct('ID',RunOpt.ID,'X',t,'Y', ...
    double(t<=R1.n0/R1.f0).*R1.A/2.* ... % n periods
      (1-cos(t*R1.f0*2*pi/R1.n0)).* ...   % hanning
      cos(t*R1.f0*2*pi), ... % cos
     'name','CosHan','Edit',[]);
 end
 
%% #TestEval  eval_str(t) - - - - - - - - - - - - - - - - - - - - - - - - - - - 
elseif comstr(Cam,'eval');[CAM,Cam]=comstr(CAM,5);  
  if ~isempty(Cam) 
  elseif isfield(opt,'EvalString');CAM=opt.EvalString;Cam=CAM;
  elseif isfield(opt,'Edit')&&isfield(opt.Edit,'st');CAM=opt.Edit.st;Cam=CAM;
  else; CAM='ones(size(t))';Cam='unit response';
  end
  if isfield(CAM,'value');CAM=CAM.value;end
  out=struct('X',[],'Y',sprintf('TestEval %s',CAM),'ID',RunOpt.ID, ...
      'name',Cam);
  if ~isempty(t);
    eval(['r1=' CAM ';']); out.X=t;out.Y=r1(:); 
  end
  
%% #TestSin #TestCos #TestTan #TestExp ... - - - - - - - - - - - - - - - - - -
% Test[Sin, Cos, ...] Period NStep Amplitude TotalTime -stoptime Tf % OpenFEM call
% Test[Sin, Cos, ...] Period Amplitude -stoptime Tf % OpenFEM call
% Test[Sin, Cos, ...] T=Period A=Amplitude -stoptime Tf % SDT call
elseif any(strncmp(Cam,{'sin','cos','tan','exp','tri','squ'},3)); 

  f1=Cam(1:3); % function name (sin, ...)
  if isstruct(R1);R1=fe_def('cleanentry',R1);
  else % old OpenFEM call
   if isempty(opt); R1=struct;
   elseif length(opt)==2; R1=struct('T',opt(1),'A',opt(2)); % T A
   elseif length(opt)==4; %Period NStep Amplitude TotalTime
    if isempty(t);  
     t=linspace(0.,opt(4),opt(2))';
     R1=struct('T',opt(1),'A',opt(3));
    else
     error(['Can''t deal with t vector and implicit t definition ',...
            '(4 parameter call) at the same time']);
    end
   end
   RunOpt.SetVal=1;
  end
  
  if isempty(Tf);Tf=0; end
  if isempty(t) % only def of the curve
   out=struct('X',[],'Y',sprintf('Test%s',CAM), ...
               'name',sprintf('%s(X)',f1));
  else%if any(exist(f1,'file')==[2 3 5 6])||exist(f1,'builtin') % return X Y
   if isempty(fieldnames(R1)) %isempty(opt)
    out=struct('X',t,'Y',feval(f1,t),...
      'name',sprintf('%s(X)',Cam),'type','tabulated');
   elseif isfield(R1,'T')&&~isempty(R1.T)&&isfield(R1,'A')
     out=struct('X',t,'Y',R1.A*feval(f1,2*pi*t/R1.T),...
               'name',sprintf('%s(X)',f1),'type','tabulated');
   elseif isfield(R1,'f'); 
     out=struct('X',t,'Y',R1.A*feval(f1,2*pi*t*R1.f),...
               'name',sprintf('%s(X)',f1),'type','tabulated');
   elseif isfield(R1,'f0')&&isfield(R1,'A')
     out=struct('X',t,'Y',R1.A*feval(f1,2*pi*t*R1.f0),...
               'name',sprintf('%s(X)',f1),'type','tabulated');
   else; error('incorrect parameter input')
   end
  end
  if ~isempty(Tf)&&Tf>0&&~isempty(out.X); out.Y(out.X>=Tf)=0; end

  out.ID=RunOpt.ID;
%% #TestBurstRandom - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
elseif comstr(Cam,'burstrandom'); [CAM,Cam]=comstr(CAM,12); % TestBurstRandom
 R2=fe_def('cleanentry',R1);
 
 i1=floor(length(t)*R2.Ratio);
 y=zeros(length(t),1);
 y(1:i1)=R2.A*randn(i1,1);
 out=struct('X',t,'Xlab',{{'Time' 's' []}},...
  'Y',y,'Ylab','Amplitude','name','BurstRandom','Edit',R1);

%% #TestNoise - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
elseif comstr(Cam,'noise'); [CAM,Cam]=comstr(CAM,6); % TestNoise
 R2=fe_def('cleanentry',R1);
 % there is no openFEM call ?
 if ~isempty(t)
  fs=1/diff(t(1:2)); % XXX t is assumed not to have uneven time step
  M=length(t);
  if mod(M,2)==0; M=M/2; RunOpt.Even=1; else;M=(M+1)/2; end
  f=linspace(0,fs/2,M+1);%linspace(0,fs/2,M)
  f(end)=[];
  Shape=fe_curve(['Test' R2.window],R2,f);
  %Shape=fe_curve(sprintf('Test%s',R2.window),f);
  %Shape=fe_curve('test',R1,f);
  out=fe_curve('noise',length(t),fs,Shape.Y);
  out.name='noise'; out.Edit=R1;
 else
  out=struct('ID',RunOpt.ID,'X',[],'Y','TestNoise','name','noise','Edit',[]);
 end

%% #TestBox - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
elseif comstr(Cam,'box'); [CAM,Cam]=comstr(CAM,4);

 if isstruct(R1);R1=fe_def('cleanentry',R1);
 else % OpenFEM calls
  if length(opt)==2; % min max
    R1=struct('min',opt(1),'max',opt(2),'A',1);
  elseif length(opt)==3; % A min max
    R1=struct('min',opt(2),'max',opt(3),'A',opt(1));
  else;error('You must provide A min max');
  end   
  RunOpt.SetVal=1;
 end
 
 if ~isempty(t)
  y=R1.A*ones(length(t),1);
  y(t<R1.min|t>R1.max)=0;
  out=struct('ID',RunOpt.ID,'X',t(:),'Y',y(:),'name','box');
  out.type='tabulated';
 else
  out=struct('ID',RunOpt.ID,'X',[],'Y','TestBox','name','box','Edit',[]);
 end

%% 
else
 if ~exist('out','var');out=opt;end
 if isfield(out,'Y');
   if iscell(out.X)&&isequal(out.X{1},t)
   elseif isequal(out.X,t)
   else; 
     if ~ischar(out.Y);sdtw('_nb','Reinterpolated Y curve with ''ReturnY''');end
     out.Y=fe_curve('returny',out,t);out.X=t;
   end
 else;sdtw('_nb','''Test%s'' is not a known command',CAM);
 end
end %  'TestFcn' 

if ~isempty(RunOpt.Edit);
    out.Edit=RunOpt.Edit;
 if RunOpt.SetVal
    st=fieldnames(R1); 
    for j1=1:length(st);
      r2=R1.(st{j1});out.Edit.(st{j1}).value=r2;
    end
 end
end
out.PlotInfo={'sub','1 1'; 'show','real'; 'scale','xlin;ylin'};
if ~isempty(t)&&RunOpt.Interp % interp
   out.Y=fe_curve('returny',out,t);  out.X=t;
   out=feutil('rmfield',out,'Interp');
end
if sp_util('issdt')
 if (~isempty(t)||~isempty(out.X))&&~isfield(out,'Xlab'); out.Xlab={'Time'}; end
 out=fe_curve('fixcurve',out);
end
if nargout==0&&~isempty(out); fe_curve('plot',out);end
  
%---------------------------------------------------------------
% This validates the structure of a given curve structure and makes any
% appropriate modification
elseif comstr(Cam,'fixcurve')
 %% #FixCurve: cleanup curve structure
  
  curve=varargin{carg}; carg=carg+1;
  st=fieldnames(curve);
  
  % This curve is a SDT <5.1 frequency/time response curve
  if length(intersect(st,{'w','xf'}))==2
    out = struct('X',iigui(curve.w),'Y',iigui(curve.xf), ...
      'zunit',[],'name',[],'type','fe_curve','unit','USER', ...
      'Interp','','PlotFcn',[]);
    if isfield(curve,'dof')
     out.data=iigui(curve.dof);
     if size(out.data,2)>=8; out.Z=out.data(:,8);end %Universal file 58 : Z axis
    end
    st=setdiff(st,{'w','xf','x','yd','yn','dof'});
    for j1=1:length(st); out=setfield(out,st{j1},getfield(curve,st{j1}));end
  else; out=curve;
  end
  if ~iscell(out.X)&&size(out.X,1)==1&&size(out.X,2)>1;out.X=out.X(:);end
  if ~ischar(out.Y)&&size(out.Y,1)==1&&size(out.Y,2)>1;out.Y=out.Y(:);end

  if ~isfield(out,'name'); out.name=''; end
  if isfield(curve,'name'); out.name=curve.name;end
  if isfield(curve,'x'); out.Xlab=curve.x; end
  try;
   if isfield(curve,'yn')&&isfield(curve,'yd');
   for j1=1:size(curve.yn,1)
    r2=curve.yn{j1,2}-curve.yd{j1,2};r2(1)=0;
    out.yunit(j1,1:3)={sprintf('%s/%s',curve.yn{j1,1},curve.yd{j1,1}), ...
     sprintf('%s/%s',curve.yn{j1,2},curve.yd{j1,2}),r2};
   end
   end
  end

  out=fe_curve('datatype fix',out);
  
  if isfield(out,'Xlab') % clean SDT curve format output
   if ~iscell(out.Xlab); out.Xlab={out.Xlab}; end
  end
  if isfield(out,'Y')&&~iscell(out.X); out.X={out.X(:)}; end
  
%% #Type ---------------------------------------------------------------
elseif comstr(Cam,'types')
   
out={'Off',''};
r1={'Eval'};
r1(:,2)={struct( ...
 'EvalString',cingui('ParamEdit','EvalString','1','%s','fcn(t)'))};
out(end+(1:size(r1,1)),1:2)=r1;
r1={'sine';'square';'triangle'};
r1(:,2)={struct( ...
 'Frequency',cingui('ParamEdit','Frequency','1','%g','Hz'),...
 'Amplitude',cingui('ParamEdit','Amplitude','1','%g','V'))};
out(end+(1:size(r1,1)),1:2)=r1;
r1={'random';'pseudo random';'DC'};
r1(:,2)={struct( ...
 'Amplitude',cingui('ParamEdit','Amplitude','1','%g','V'))};
out(end+(1:size(r1,1)),1:2)=r1;
r1={'burst random'};
r1(:,2)={struct( ...
 'Length',cingui('ParamEdit','Length','1','%g','n non zero samples'), ...
 'Amplitude',cingui('ParamEdit','Amplitude','1','%g','V'))};
out(end+(1:size(r1,1)),1:2)=r1;

%% #convert -unit "str" ---------------------------------------------------
elseif comstr(Cam,'convert'); % 'list'  - - - - - - - - - - - - - - -

[CAM,Cam,RO.unit]=comstr('unit',[-25 4],CAM,Cam);
def=varargin{carg};carg=carg+1;

if isfield(def,'yn')
 %st=fe_curve('')    
end
out=def;

%% #Resample : lininterp of samples at mean time step  -----------------------
elseif comstr(Cam,'resample'); % 'list'  - - - - - - - - - - - - - - -

 [CAM,Cam,RO.unit]=comstr('unit',[-25 4],CAM,Cam);
 R1=varargin{carg};carg=carg+1;
 
 tn=min(R1.X{1}):mean(diff(R1.X{1})):max(R1.X{1}); tn=tn(:);
 R1.Y=of_time('lininterp',[R1.X{1} R1.Y],tn,zeros(1,3)); % resample with regular time step
 R1.X{1}=tn;
 out=R1;

%---------------------------------------------------------------
elseif comstr(Cam,'list'); % 'list'  - - - - - - - - - - - - - - -
    
 out=TestList;
 if ~isempty(findstr(Cam,'-pop'));testlist=cat(1,{'Select',''},out); end
 if nargout==0
    out=out'; fprintf('%15s : %s\n',out{:});clear out
 end
%% #End -----------------------------------------------------------------
elseif comstr(Cam,'cvs')  
  out='$Revision: 1.235 $  $Date: 2021/10/29 11:01:30 $';
%---------------------------------------------------------------
else;error('''%s'' is not a known command',CAM);    
end;

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% #SubFunc -----------------------------------------------------------------
%% #TestList : analyze options of commond curves -----------------------------
function [out,out1,out2]=TestList(CAM,RunOpt);

 list={...
   'Step','t1(1#%g#"Final time")'
   'Ramp',['t0(0#%g#"Initial time") t1(1#%g#"Final time") ' ...
           'Yf(1#%g#"Final value") NStep(#%g#"number of steps")'];... 
   'Ricker',['dt("t(end)/10"#%g#"Duration") '...
             'A(.1#%g#"Amplitude") t0(0#%g#"Start time")']; ...
   'Sweep',['fmin(100#%g#"Initial frequency") ',...
            'fmax(110#%g#"Final frequency") ',...
            't0(0#%g#"Start time") t1(1#%g#"Stop time")',...
            'fcn("cos"#%s#"sin, cos, esweepcos, esweepsin") A(1#%g#"Amplitude") ']; ...
   'Noise','window(Box#%s#"Window") A(1#%g#"Amplitude") min(#%g#"min abscissa") max(#%g#"max abscissa")'; ...
   'Box',['A(1#%g#"Amplitude") min(0#%g#"min abscissa") ',...
          'max(1#%g#"max abscissa")'];... 
   'sine','f(1#%g#"Frequency (Hz)") A(1#%g#"Amplitude")'; ...
   'triangle','f(1#%g#"Frequency (Hz)") A(1#%g#"Amplitude")'; ...
   'square','f(1#%g#"Frequency (Hz)") A(1#%g#"Amplitude")'; ...
   'random','A(1#%g#"Amplitude") End0(0#%g#"end samples at 0")'; ... 
   'BurstRandom','A(1#%g#"Amplitude") Ratio(.6#%g#"ratio with noise")'; ... 
   'Eval','st("ones(size(t))"#%s#"fct(t) expression")'
   'sin','T(1#%g#"Period") A(1#%g#"Amplitude")'; ... 
   'cos','T(0#%g#"Period") A(1#%g#"Amplitude")'; ...
   'tan','T(0#%g#"Period") A(1#%g#"Amplitude")'; ...
   'exp','T(0#%g#"Period : A.exp(2*pi*t/T)") A(1#%g#"Amplitude")';
   'CosHan',['f0("1/t(2)/10"#%g#"Frequency (Hz)") n0(10#%g#"N period") ',...
             'A(1#%g#"Amplitude")']; ...
   };

if nargin==0&&nargout>0; out=list;return;
elseif isempty(CAM)&&nargin==2&&isfield(RunOpt,'RO')&&isfield(RunOpt.RO,'type');
  CAM=RunOpt.RO.type;
elseif nargin==0||isempty(CAM)
  list=list';
  fprintf(1,'\nAvailable Test curves:\n%s',sprintf('Test%-15s %s\n',list{:}));
  if nargout>0; out=[];out1=[];out2=[];end;return
end
i1=strcmpi(list(:,1),CAM);
if any(i1); tag=CAM;
else; % compatibility with tag(number)
   tag=sscanf(regexprep(CAM,'[()\s"\.\d]',' '),'%s',1);
   i1=strcmpi(list(:,1),tag);
end
if ~any(i1);list=[];else;list=list{i1,2};end
if nargin>1
  %% Possibly check Edit
  if ~isfield(RunOpt,'RO');RunOpt.RO=struct;end
  out=list;
  if isempty(list);Edit=struct;st2=CAM;
  else; 
      [Edit,st,st2]=cingui('paramedit',list,{RunOpt.RO,CAM(length(tag)+1:end)});
  end
  if comstr(CAM,'eval')&&~strcmpi(tag,CAM);
    % fe_curve('testEval -id1 sin(t)',linspace(0,pi,10))
      Edit.st=comstr(CAM,5);RunOpt.CAM='eval';
  elseif ~strcmpi(CAM,tag)&&any(ismember(st2,'0123456789'))&&~isempty(strfind(CAM,st2)); %#ok<STREMP>
     % No edit but some info, backward compat 
     % c1=fe_curve('test -ID 100 ricker 10e-4 100 1 100e-4');
     opt=comstr(CAM(length(tag)+1:end),-1);RunOpt.CAM=tag;
     st=fieldnames(Edit);
     for j1=1:min(length(st),length(opt)); Edit.(st{j1}).value=opt(j1);end
     %disp(comstr(out1,-30))
  else;RunOpt.CAM=st2;
  end
  out=struct('ID',RunOpt.ID,'X',[],'Y',['Test' tag],'Edit',Edit);
  out1=Edit;
  if isfield(Edit,'NStep')&&~isempty(Edit.NStep.value); 
    r1=fe_def('cleanentry',Edit);  RunOpt.t=linspace(r1.t0,r1.t1,r1.NStep)';
  end
  out2=RunOpt;
else; out=list;
end

%% #linsweep
function out=linsweep(t,R1); 

  ind=t>=R1.t0&t<=R1.t1;r2=t(ind);
  r3=zeros(size(t));
  if any(ind) % Linear sweep see wikipedia, f(t)=f0+kt
   ti=t(ind)-R1.t0;
   if ischar(R1.fcn)&&any(R1.fcn=='/');
       R1.fcn=eval([strrep(R1.fcn,'/','(''@') ''')']);
   end
   r4=2*pi*(R1.fmin*ti+(R1.fmax-R1.fmin)/2/(ti(end)-ti(1))*ti.^2);
   r3(ind)=feval(R1.fcn,r4);
   % see also : sdtweb fe_curve ecos
  end
out=r3;

%% #ecos : exponential sweep
function out=esweepcos(t,R1); %#ok<DEFNU>


    % Choice of m (in order that the phase property is valid)
    R1.m = ceil((2*pi*(R1.t1-R1.t0)*R1.fmin/log(R1.fmax/R1.fmin)+pi/2)/(2*pi));

    % Actual length in samples
    R1.T_act = (2*R1.m*pi - pi/2)*log(R1.fmax/R1.fmin)/(2*pi*R1.fmin) ;

    % Computation of the phase
    phi = 2*pi*R1.fmin*R1.T_act/log(R1.fmax/R1.fmin)* ...
        (exp(t/R1.T_act*log(R1.fmax/R1.fmin))-1)-pi/2;

    % Computation of the sweep
    out = cos(phi).*double(t<=R1.T_act);
    
%% #esin : exponential sweep
function out=esweepsin(t,R1);%#ok<DEFNU>

    % m parameter
    R1.m = round( ((R1.t1-R1.t0)*R1.fmin)/(log(R1.fmax/R1.fmin)) ); 
    R1.T_act = R1.m/R1.fmin*log(R1.fmax/R1.fmin);

    % Computation of the phase
    phi = 2*pi*R1.m*(exp(R1.fmin*t/R1.m)-1);
    out = sin(phi).*double(t<=R1.T_act);

    % computation of the sweep
    out = zeros(size(t)) ;
    out(t<=R1.T_act) = sin(phi(t<=R1.T_act));
    

%% #DataType -----------------------------------------------------------------
function r2 = DataType(r1,def,unit) %#ok<INUSD>

% def is a char, unit label provided
if nargin<2||isempty(def); def={0,'None'}; end
if nargin>0&&ischar(r1)&&~isempty(strfind(r1,'-cell'));
 OutTyp='cell';r1(strfind(r1,'-cell')+[0:4])='';r1=comstr(r1,1);
else;OutTyp='struct';
end

% Typ [Length Force Temp Time exponents]
DataTypes={[0 0 0 0 0],'Unknown','NONE';
  [1 0 0 0 0],'General','NONE';
  [2 -2 1 0 0],'Stress','Pa';
  [3  0 0 0 0],'Strain','';
  [5 0 0 1 0],'Temperature','^o';
  [6 1 1 0 0],'Heat flux','NONE';
  [8  1 0 0 0],'Displacement','m';
  [9  0 1 0 0],'React. force','N';
...  [9  0 1 0 0],'Force','N'; Not UFF label, use fe_curve('datatypecell-label"Force"',9)
  [11 1 0 0 -1],'Velocity','m/s';
  [12 1 0 0 -2],'Acceleration','ms^{-2}';
  [13 0 1 0 0],'Excit. force','N';
  [15 -2 1 0 0],'Pressure','Pa';
  [16 -1 1 0 0],'Mass','Kg';
  [17 0 0 0 1],'Time','s';
  [18 0 0 0 -1],'Frequency','Hz';
  [19 0 0 0 -1],'RPM','RPM'; 
  [20 0 0 0 0],'Order','?';    % ABOVE HERE UFF SPECIFICATION
  [100 -1 1 0 0],'Moment','N'; % BELOW HERE OPENFEM Specific
  [101 0 0 0 0],'Rotation','rad';
  };

i1=reshape([DataTypes{:,1}],5,size(DataTypes,1))';
if nargin==0
  r2=struct('label',[],'value',i1);r2.label=DataTypes(:,2);
  out=r2;   return;
elseif isnumeric(r1) % match against number
  if isempty(r1);i2=[];else;i2=find(i1(:,1)==r1(1)); end
  if ~isempty(i2); r1=DataTypes(i2,[2 3 1]); else; r1=def;end
elseif isstruct(r1)
elseif iscell(r1)&&size(r1,2)==3  % already cell format do nothing
elseif comstr(comstr(r1,-27),'_datatype')
  r2=struct('label',[],'value',i1);r2.label=DataTypes(:,2); return; 
elseif ischar(r1)
  i2=strmatch(comstr(r1,-27),comstr(DataTypes(:,2),-27));
  if isempty(i2); sdtw('_nb',sprintf('''%s'' not matched',r1));r1=def;
  else
   r1=DataTypes(i2,[2 3 1]);
   r1=struct('label',r1{1,1},'unit',r1{1,2},'type',r1{1,3});
  end
  try
  if ischar(def); r1.unit=def; 
  elseif ~isfield(def,'xuCell')||~strncmp(r1.label,'Fre',3)||~iscell(def.xuCell)
  elseif (any(def.xuCell{2}=='m')||(length(def.xuCell)>2&&def.xuCell{3}(2)==1))
     r1.unit=['1/' def.xuCell{2}];r1.type=[18 -1 0 0 0];
  elseif strcmpi(def.xuCell{2},'deg')||strcmpi(def.xuCell{2},'rev')
     r1.unit=['1/' def.xuCell{2}];r1.type=[18 0 0 0 0];
  end
  end 
else; r1=def;
end
if isa(r1,'cell') % Newer cell array format
  if size(r1,2)<3; r1(:,3)={[0 0 0 0 0]};end
  for j1=1:size(r1,1);
   if isempty(r1{j1,3})||(r1{j1,3}(1)==0 &&  def{1}); r1{j1,3}(1)=def{1}; end
   if  r1{j1,3}(1)~=0; r1{j1,3}(2:5)=i1(i1(:,1)==r1{j1,3}(1),2:5); end
   % MISSING label checks
  end
elseif isstruct(r1)&&isfield(r1,'Xlab')&&isfield(r1,'Y')
 %% Possibly fill x,yn,yd based on Xlab information
  try; r2=r1.Xlab{1}; 
      if iscell(r2);r3=DataType(r2{1});r3.unit=r2{2};
      else; r3=DataType(r2);
      end
      r1.x=r3;
  end
  try; % Assume of the forme 'displacement / force'
      r2=r1.Xlab{2}; 
      if iscell(r2);r3=DataType(regexprep(r2{1},'[ ]*(/.*)',''));
          r3.unit=regexprep(r2{2},'[ ]*(/.*)','');
      else; r3=DataType(r2);
      end
      r1.yn=r3;
      if iscell(r2);r3=DataType(regexprep(r2{1},'(.*/)[ ]*',''));
          r3.unit=regexprep(r2{2},'(.*/)[ ]*','');
      else; r3=DataType(r2);
      end
      r1.yd=r3;
  end
  
elseif isstruct(r1)
 % checks on type
 if r1.type(1)==0 &&  def{1}; r1.type(1)=def{1}; end
 if ~any(r1.type(1)==i1(:,1)); r1.type(1)=0; end
 if r1.type(1)~=0;  r1.type(2:5)=i1(find(i1(:,1)==r1.type(1),1),2:5); end

 if strcmpi(r1.label,'none')||strcmpi(r1.label,'unknown')||isempty(r1.label)
   r1.label=DataTypes{(i1(:,1)==def{1}),2};
 end
 if strcmpi(r1.label,'none')||isempty(r1.label) 
   r1.label=DataTypes{(i1(:,1)==r1.type(1)),2};
 end
 if strcmpi(r1.unit,'none')||isempty(r1.unit); r1.unit=def{2}; end
end

if strcmp(OutTyp,'cell')&&isstruct(r1);
 r2={r1.label,r1.unit,r1.type};
else; r2=r1;
end


%---------------------------------------------------------------
function  [CAM,Cam,RunOpt]=StandardOptions(CAM,Cam);

RunOpt=struct('type','');
if any(Cam=='-')
  [CAM,Cam,RunOpt.Window]=comstr('-window',[-25 4],CAM,Cam);
  [CAM,Cam,RunOpt.SkipTime]=comstr('-skiptime',[-25 2],CAM,Cam);
  [CAM,Cam,RunOpt.Skip]=comstr('-skip',[-25 1],CAM,Cam);
  [CAM,Cam,RunOpt.Overlay]=comstr('-overlay',[-25 2],CAM,Cam);  
  [CAM,Cam,RunOpt.BufTime]=comstr('-buftime',[-25 2],CAM,Cam);
end

%% #Tri : triangular input
function Y=tri(t);

Y=t/2/pi; Y=rem(Y,1); Y=min(Y,1-Y);

%---------------------------------------------------------------
%---------------------------------------------------------------
function [r1,list,cur]=get_curve(GEF,cur)

[r1,ind]=stack_get(GEF,'curve');
list=r1(:,2);r1=r1(:,3);

if isempty(r1); i1=0; list={'New ...'};r1=cell(1,1);
else;list{end+1,1}='New ...'; r1{end+1,1}=[];end

if nargin==1; cur=[];elseif cur<=length(ind); cur=ind(cur); else;cur=[];end
