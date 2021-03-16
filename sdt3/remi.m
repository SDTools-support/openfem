function x=remi(x,y,pos)

%REMI	Remainder after division for use with integer indices
%	remi(x,y) is x - n.*y where n=fix(x./y) and the elements of
%	x-n.*y equal to 0 are set to y
%   remi(x,[],pos) returns the (pos) digit in x

%       Etienne Balmes  08/25/94, 04/02/97
%       Copyright (c) 1990-2009 by SDTools,
%       All Rights Reserved.
%       $Revision: 1.3 $  $Date: 2015/03/04 11:20:19 $

if nargin==2
    y=y(1);x = rem(x,y);i1 = find(x==0); x(i1)=x(i1)+y;
else
    pos=10.^(pos-1);
    for j1=1:size(x,1);x(j1,1:length(pos))=rem(fix(x(j1,1)./pos),10);end
end
