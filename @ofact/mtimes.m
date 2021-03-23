function [u]=mtimes(a,b);

%MTIMES matrix multiplication
%
%       Synopsis: u = sky_mul(k,ind,v)
%	U = SKY_MUL(M,IND,V) is the equivalent for a skyline matrix of
%	V.'*M for a full matrix.
%
%	See also: sky_inv, sky_dec, sky2sp

%	Etienne Balmes
%       Copyright (c) 2001-2008 by INRIA and SDTools, All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       $Revision: 1.4 $  $Date: 2008/10/28 08:56:32 $


if isa(a,'ofact'); nk=length(a.ind)-1;
    if nk==-1; error('ofact multiplication is not defined for %s solver', ...
            a.method.name); end
   if ~isa(b,'double'); error('Not a skyline/vector product');end
   nv = size(b,1); if nv~=nk; error('k and b must have as many rows');end
   if a.ty(1)~=0; error('Factored matrix product not implemented'); end
   u = skymul(a,b'); u = u';

elseif ~isa(b,'ofact'); error('Not a skyline/vector product');
else 
   if ~isa(a,'double'); error('Not a vector/skyline product');end
   nk = length(b.ind)-1;
   if nk==-1; error('ofact multiplication is not defined for %s solver', ...
           a.method.name); end
   nv = size(a,2); if nv~=nk; error('b and k must have as many columns');end
   u = skymul(b,a);

end

% do the multiplication - - - - - - - - - - - - - - - - - - - - - - - - 

function u = skymul(k,v);

 u = zeros(size(v));

for j1 = 1:length(k.ind)-1

	kl=j1+[k.ind(j1)-k.ind(j1+1)+1:0];
	u(:,j1)=v(:,kl)*k.data(k.ind(j1+1)-1:-1:k.ind(j1));
   if k.ind(j1)+1<k.ind(j1+1)
	kl=j1+[k.ind(j1)-k.ind(j1+1)+1:-1];
	ku=k.ind(j1+1)-1:-1:k.ind(j1)+1;
	u(:,kl)=u(:,kl)+v(:,j1)*k.data(ku)';
   end

end
