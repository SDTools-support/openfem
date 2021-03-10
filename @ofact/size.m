function a = size(k,b);

%       Etienne Balmes  07/28/92, 07/22/01
%       Copyright (c) 2001-2009 by INRIA and SDTools
%       $Revision: 1.3 $  $Date: 2009/10/09 15:49:17 $

%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       All Rights Reserved.

if k.ty(1)==2 || k.ty(1)==3
   if nargin==1; a = size(k.l);
   elseif nargin==2; a=size(k.l,b); end
elseif any(k.ty(1)==[5 6 6.1 6.2])
   if nargin==1; a = [1 1]*k.ty(3);
   elseif nargin==2; a = k.ty(3); end
else
   if nargin==1; a = [1 1]*length(k.ind)-1;
   elseif nargin==2; a = length(k.ind)-1; end
end