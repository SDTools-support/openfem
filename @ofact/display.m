function display(k,val) %#ok<INUSD>

%       Display information about factored matrices

%	Etienne Balmes
%       Copyright (c) 2001-2020 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       $Revision: 1.11 $  $Date: 2020/12/03 08:31:49 $


if nargin==1; disp(' ');disp([inputname(1),' = ']); end

if k.ty(1)==2||k.ty(1)==3 % LU factors
   i1 = nnz(k.l)+nnz(k.u);st='LU'; if k.ty(1)==3; st='Cholesky';end
   st = sprintf('\n[ %i DOFs ] (%i elements %.1f MB) Sparse %s factors', ...
     size(k.l,1),i1,(i1*16+8*size(k.l,1))/1024^2,st);

elseif any(k.ty(1)==[6 6.1])
   st = sprintf('\n[ %i DOFs ] diag matrix factor',k.ty(3));
elseif any(k.ty(1)==6.2)
   st = sprintf('\n[ %i DOFs ] full inverse',k.ty(3));
elseif k.ty(1)==5
   r2=k.method; ty=k.ty;
   if isfield(r2,'TktSolve')&&~isempty(r2.TktSolve)&&length(ty)>6
    st = sprintf('\n[ %i/%i DOFs ] TkT Matrix %i %s', ...
      k.ty(3),k.ty(7),k.ty(2),r2.header);
   elseif isfield(r2,'header')
    st = sprintf('\n[ %i DOFs ] Matrix %i %s', ...
      k.ty(3),k.ty(2),r2.header);
   else
    st = sprintf('\n[ %i DOFs ] External %s solver matrix %i', ...
      k.ty(3),k.ind,k.ty(2));
   end 
elseif k.ty(1)==11 % Tkt
  st = sprintf( '\n[ %i DOFs ] %s matrix',k.ty(3), k.method.header );
elseif k.ty(1)==10
  st = sprintf( '\n[ %i DOFs ] %s matrix',k.ty(3), k.method.header );
else
  st = sprintf('\n[ %i DOFs ] (%i elements %.1f MB)', ...
     length(k.ind)-1,k.ind(end)-1, ...
     ((length(k.data)+length(k.ind)+1)*8/1024^2));
  if k.ty(1)==0; st=[st ' Nominal']; 
  elseif k.ty(1)==1; st=[st ' Factored']; 
  end

end
st = [st st(1)];disp(st)
