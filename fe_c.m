function [o1,o2,o3]=fe_c(mdof,adof,c,opt)

%FE_C DOF selection and I/O shape matrix construction
%
%	Syntax : c            = fe_c(mdof,adof)
%	         c            = fe_c(mdof,adof,cin,opt)
%		 [adof,ind,c] = fe_c(mdof,adof,cin,opt)
%		 ind          = fe_c(mdof,adof,'ind')
%		 adof         = fe_c(mdof,adof,'dof')
%		 labels       = fe_c(mdof,adof,'dofs')
%		 b            = fe_c(mdof,adof)'
%
%	FE_C is used both
%        - to characterize sensors/actuators by the output C or input B matrix
%        - to select DOFs (in particular to impose boundary conditions)
%
%	The input arguments are
%
%	MDOF main DOF definition vector (see sdtweb('mdof')).
%	Each element of this column vector defines a DOF following the format 
%	   NodeId.DofId for nodal DOFs (with DOFs 01 to 99 accepted). For
%	      example 100.01 indicates the x translation at node 100)
%	      By default DOFs .01 to .06 are xyz translations/rotations.
%	      DOFs .07 to .12 are reserved for -x-y-z translations/rotations
%	      but should not be used in finite element analysis.
%	  -EltId.DofId  for element DOFs (DOFs 1 to 999 are accepted). The 
%	      EltId value is discussed sdtweb('mdof').
%	ADOF gives the DOF that are used for the initial definition of c
%            Accepted simplifications (wild cards) are
%	       10.0 all DOFs (01 to 06 in general) of node 10
%		0.1 x-translations at all nodes
%	      -10.0 all internal DOFs of element 10
%	CIN  describes different outputs (one per row) in the DOFs ADOF
%	     If not specified or empty CIN is taken to be the identity matrix
%	     If the strings 'dof' or 'ind' are used instead of CIN the first
%	     and only output argument is ADOF or IND respectively.
%	OPT   (optional)
%	     OPT = 1  keeps the DOFs in ADOF (this is the default)
%	     OPT = 2  keeps the DOFs that are not in ADOF
%	
%	The output arguments are
%
%	C    the output shape matrix in model DOFs described by MDOF
%	     (if only one output argument is asked, this argument is C and
%	      not the expanded ADOF)
%	ADOF expanded (without wild cards) version of the input
%	     argument ADOF (or its complementary if OPT==2)
%	IND  indices of the ADOF DOFs in MDOF (i.e. ADOF = MDOF(IND))
%
%	NOTE: FE_C supports the convention that nodal DOFs .07 to .12 are
%	     the opposite of nodal DOFs .01 to .06
%
%	See sdtweb      mdof, adof
%	See also help   fe_eig, fe_mk, feplot, fe_coor, fe_load, nor2ss

%	Etienne Balmes
%       Copyright (c) 2001-2020 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       $Revision: 1.29 $  $Date: 2020/02/06 16:28:53 $

%#ok<*NOSEM>
if comstr(mdof,'cvs')
 o1='$Revision: 1.29 $  $Date: 2020/02/06 16:28:53 $'; return;
end

if nargin==1
 o1 = feutil('stringdof',mdof); return;
end
if size(adof,1)==1 && size(adof,2)~=1; adof=adof(:);
   warning('active DOFs should be specified as a column');
end
if size(mdof,1)==1 && size(mdof,2)~=1; mdof=mdof(:);
   warning('mdof should be specified as a column');
elseif size(mdof,2)>1; mdof=mdof(:,1); 
end
if nargin<3; c = []; end
if nargin<4; opt=1; elseif isempty(find([1 2]==opt)); opt = 1; end

if ~isempty(adof)
  an=fix(adof(:,1));  % NodeId or EltId
  ad=round((adof(:,1)-an)*1000); % DofId
else;ad=[];an=[]; 
end



in2=find(ad>69&ad<121);
if ~isempty(in2); adof(in2)=adof(in2)-60;ad(in2)=ad(in2)-60; end
if ~isempty(mdof);
 mdofn=fix(mdof); mdofd=round((mdof-mdofn)*1000);
 if nnz(mdofd)==0; mdofd(:,1)=1;end
 in1=find(mdof>0&mdofd<121&mdofd>69);
 if ~isempty(in1); mdof(in1)=mdof(in1)-.06;mdofd(in1)=mdofd(in1)-60; end
 ind=find(mdof>0);cind=find(mdof<0);
 if isempty(cind)&&isequal(c,'indskipcheck');c='ind';
 elseif isempty(cind)
  try;   %i1=find(sparse(mdofn*1000+mdofd,1,1)>1); % 4-5 times quicker than catch
   i1=sort(mdofn*1000+mdofd); % 2x quicker than line above
   i2=i1; i1=find(diff(i2)==0); i3=1;
  catch % xxx32Bit ?
   [i1,i2]=find(sparse(mdofn,mdofd,ones(size(mdof,1),size(mdof,2)))>1);
   i1=i1*1000+i2; i3=0;
  end
  if ~isempty(i1) 
   if i3; i1=i2(i1+1); end
   fprintf('%.2f ',fix(i1/1000)+rem(i1,1000)/1000);fprintf('\n');
    if max(mdof)>2^31; error('Current max node number is 2^31/100')
    else; error('Repeated DOFs in mdof'); 
    end
  end
 else;  
  i1=diff(sort(mdof));if ~all(abs(i1)>1e-4); error('Repeated DOFs in mdof');end
 end
else; in1=[];ind=[]; mdofn=[];mdofd=[]; 
end


rdof = []; opt(2)=0;

% mdof=1+[1:6]'/100;mdof=[mdof;mdof+2];fe_c2(mdof,.01,'ind')

% no wild cards - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

i5=find(ad&an);
if ~isempty(i5)
 if length(mdofn)+length(an)<1e4
  %%[i6,i7]=sort([mdof;adof(i5)]);
  %%[i6,i7]=sortrows([mdofn mdofd;an(i5) ad(i5)]);
  [i6,i7]=sort([mdofn*1000+mdofd;an(i5)*1000+ad(i5)]);
  %position in i7 of retained mdof dofs
  i8=find(~diff(i6));
  i8=i8(i7(i8)<=length(mdof));
  [i9,i10]=sort(i7(i8+1)); % sort using first appearance in adof
  rdof=i7(i8(i10));        % indices of retained DOFs in mdof
  i9=i9-length(mdof);      % sorted existing positions in adof
 else % 2019/09/18 gv optim (4x) for large vectors due to progress of ismember in recent MATLAB
  [i6,i7]=ismember([an(i5)*1000+ad(i5)],[mdofn*1000+mdofd]);
  rdof=i7(i6); i9=find(i6);
  if ~isempty(rdof)
   [i6,i7]=sort(rdof); %rdof=rdof(sort(i7([1;1+find(diff(i6))])));
   i7(1+find(diff(i6)==0))=[]; rdof=rdof(sort(i7));
  end
 end
 opt(2)=1;
 if ~ischar(c)&&size(c,2)==size(adof,1)
  c=c(:,i5(i9));
 elseif ~ischar(c)
 elseif comstr(c,'place') % gvdr review 10/12/2015 optim with repeated
  i2=any(diff(sort(an*1000+ad))==0); %[i2,i4]=find(sparse(an,ad,1)>1);
  if i2 %~isempty(i2) % deal with repeated DOFs in adof.
   ii=sparse(mdofn*1000+mdofd,1,1:length(mdof));
   if comstr(c,'placei'); o1=full(ii(an*1000+ad));
   else;o1=sparse(i5,full(ii(an*1000+ad)),1,length(adof),length(mdof));
   end
  elseif comstr(c,'placei'); o1=rdof;
  else; o1=sparse(i5(i9),rdof,1,length(adof),length(mdof));
  end
  if ~isempty(in2); o1(in2,:)=-o1(in2,:); end
  if ~isempty(in1); o1(:,in1)=-o1(:,in1); end
  %     [i2,i4]=find(sparse(an,ad,1)>1);
  %     if ~isempty(i2) % OLD AND SLOOOOW
  %      % deal with repeated DOFs in adof.
  %        for j1=1:length(i2)
  %         i3=find(ismember([an ad],[i2(j1) i4(j1)],'rows'));
  %         %i3=find(adof==i2(j1));
  %         o1(i3(2:length(i3)),:)=o1(i3(ones(length(i3)-1,1)),:);
  %        end
  %     end
  return
 elseif comstr(c,'ind')
  o2=i5(i9);
 end %if some DOF eliminated
end

% NodeDOF ID wild card  - - - - - - - - - - - - - - - - - - - - - - - - - -

i5=find(ad>0&~an);
if ~isempty(i5)
  i6=[];i6(ad(i5),1)=[1:length(i5)]';
  i7=mdofd(ind);if length(i6)<max(i7); i6(max(i7))=0;end
  rdof=[rdof;ind(i6(i7)~=0)];
  opt(2)=0;
end

% EltDOF ID wild card  - - - - - - - - - - - - - - - - - - - - - - - - - -

i5=find(ad<0&~an);
if ~isempty(i5)
  i6=[];i6(-ad(i5),1)=[1:length(i5)]';
  i7=-mdofd(cind);%i7=-rem(mdof(cind),1000);
  if length(i6)<max(i7); i6(max(i7))=0;end
  rdof=[rdof;cind(i6(i7)~=0)];
end

% NodeId wild card - - - - - - - - - - - - - - - - - - - - - - - - - - - -

i5=find(~ad&an>0);
if ~isempty(i5)
  %fprintf('%.0f %.0f\n',2^31-2,2^48-1)
  i6=2147483646; if of_mk('mwIndex')==8; i6=281474976710655;end
  i6=sparse(an(i5),1,1:length(i5),i6,1);%i6=[];i6(an(i5),1)=[1:length(i5)]';
  i7=mdofn(ind); % i7=fix(mdof(ind)/1000);
  if ~isempty(i7)
    i8=find(i6(i7));i7=i6(i7(i8));  [i7,i9]=sort(i7); % declared node order
    rdof=[rdof;ind(i8(i9))];
  end
  opt(2)=0; if length(i5)==length(adof); opt(2)=1;end
end

% EltId wild card - - - - - - - - - - - - - - - - - - - - - - - - - - - -

i5=find(~ad&an<0);
if ~isempty(i5)
    i10=[round(rem(-mdof(cind),1e3)) round(-mdof(cind)/1e3)];
    i8=[round(rem(-adof(i5),1e3)) round(-adof(i5)/1e3)];
    i11=max([i10;i8])+1;i11=[i11(1) i11(2)*i11(1)];
    i12=sparse(i8(:,1)+i8(:,2)*i11(1),1,1:size(i8,1),i11(2),1);
    i7=i12(i10(:,1)+i10(:,2)*i11(1));
    rdof=[rdof;cind(i7~=0)];  opt(2)=0;
end

% Cleaning up  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% eliminates duplicates in rdof but sorts DOFs
if opt(2)~=1; rdof=find(sparse(rdof,1,rdof)); end

if opt(1)==2 % seek complementary 
  i1=1:size(mdof,1);i1(rdof)=zeros(1,length(rdof));rdof=find(i1);rdof=rdof(:);
end

if ischar(c) 
  o1=[];
  if     strcmp(c,'ind'); o1 = full(rdof); 
  elseif strcmp(c,'dof'); o1 = mdof(rdof); o2=[];
  elseif strcmp(c,'dofs') 
   o1 = feutil('stringdof',mdof(rdof)); 
   o2=[];
  end
elseif nargout==2;		o1=mdof(rdof);o2=full(rdof);
else

  if size(c,1)==0; c=speye(length(rdof),length(rdof)); end
  if size(c,2)~=length(adof) && ~isempty(in2)
    warning('Sign change for DOFs .01-.06 to .07-.12 not performed');
    in1=[];in2=[];
  end
  if length(rdof)~=size(c,2)
      error('adof does not correspond to the number of columns in c');
  end
  if ~isempty(in2);     c(:,in2)=-c(:,in2);      end
  if ~isempty(c)
    if nnz(c-speye(size(c,1),size(c,2)))||size(c,1)~=length(rdof)
     % if heavily populated use full matrix rather than big-sparse
     if nnz(c)/size(c,1)/length(mdof)>.5
       co=zeros(size(c,1),length(mdof));co(:,rdof) = c;
     else  
       % co = spalloc(size(c,1),length(mdof),1);co(:,rdof) = c;
       [i1,i2,co]=find(c);
       co=sparse(i1,rdof(i2),co,size(c,1),max(max(rdof),length(mdof)));
     end
     
    else
     co=sparse(1:size(c,1),rdof,1,size(c,1),length(mdof));
    end
    if ~isempty(in1);     co(:,in1)=-co(:,in1);      end
  else
    if size(c,1)==0; co = spalloc(length(adof),size(mdof,1),0);
    else;co=spalloc(size(c,1),size(mdof,1),0);
    end
  end
  if nargout<2; o1=co; else;o1=mdof(rdof);o2=full(rdof);o3=co; end
end

