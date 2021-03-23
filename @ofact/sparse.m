function k = sparse(ks);

%SPARSE	Skyline to sparse conversion
%
%

%       $Revision: 1.2 $  $Date: 2003/09/03 06:55:49 $
 
 n=length(ks.ind)-1;nk=ks.ind(n+1)-1;

 k=ones(size(ks.data));
 for j1=1:n
  ind = ks.ind(j1):ks.ind(j1+1)-1;
  k(ind) = n*(j1-1)+j1+[0:-1:-length(ind)+1]';
 end

 k = sparse(k,1,ks.data,n^2,1);
 k = reshape(k,n,n);
 k = k + triu(k,1)';
