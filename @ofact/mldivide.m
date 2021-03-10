function  q = mldivide(k,b)

% ofact/mldivide

%	Etienne Balmes
%       Copyright (c) 2001-2021 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       $Revision: 1.17 $  $Date: 2021/02/26 09:37:57 $

%#ok<*NOSEM> 

st=k.method; if isstruct(st);st=st.Solve;end
if ischar(b)
  if strcmpi(b,'fact');q='done';return;
  else;eval(iigui(b,'MoveFromCaller'));
  end
end
if iscell(st);q=feval(st{:},k,full(b));return; end

switch k.ty(1) % sdtweb ofact('.ty')
case 6.1
   if size(b,2)==1; q=b.*k.dinv; return;
   else
     for j1=1:size(b,2);  b(:,j1)=b(:,j1).*k.dinv; end
   end
   q=b;
case 6.2; q=k.dinv*b;% the inverse is actually stored
case 5
 if ischar(b); q='done';
 elseif ~isempty(k.method.Solve); 
    if ischar(st);eval(k.method.Solve)
    elseif ischar(b);eval(iigui(b,'MoveFromCaller'));q=st('b',k.ty,k.method);
    else; 
      if ~isempty(k.ind)&&size(b,1)==size(k.ind,1)
       b=k.ind'*b; % T'*b
       q=k.ind*st(b,k.ty,k.method); 
      else;q=st(b,k.ty,k.method); % allow external solver
      end
    end
 else
  if strcmp(k.ind,'spfmex')  
    q=spfmex('solve',k.ty(2),full(b));
  elseif strcmp(k.ind,'psldlt')
    q=zeros(size(b));
    for j1=1:size(b,2); q(:,j1)=psldlt('solve',full(b(:,j1)),k.ty(2)); end
  end
 end
case 6
   for j1=1:size(b,2);  b(k.ind,j1)=b(k.ind,j1).*k.dinv; end
   q=b;
case 0
 if ischar(b); sp_util('sky_dec',k);end;q=[];
case 2
 if ischar(b); q=[]; 
 elseif ~isempty(k.ind)&&size(b,1)==size(k.ind,1)
  b=k.ind'*b; % T'*b
  q=k.ind* (k.u \ (k.l \ b) );
 else; q = k.u \ (k.l \ b);
 end
case 3
 if ischar(b); q=[];
 elseif ~isempty(k.ind)&&size(b,1)==size(k.ind,1)
  b=k.ind'*b; % T'*b
  q=k.ind* (k.u \ (k.l \ b) );
 else; q = k.u \ (k.l \ b);
 end
otherwise
 if isempty(b)||ischar(b); q=b;
 elseif ~isempty(k.method.Solve); 
   st=k.method.Solve; if ischar(st);eval(st);else;q=feval(st,k,b);end
 else
  if issparse(b); b=full(b);end
  if k.ty(1)==0;  sp_util('sky_dec',k); end
  if isreal(b); q = sp_util('sky_inv',k,b);
  else q = sp_util('sky_inv',k,real(b))+1i*sp_util('sky_inv',k,imag(b));end
 end
end 
