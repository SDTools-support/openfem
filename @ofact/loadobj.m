function out = load(obj)

%LOADOBJ load filter for OFACT  objects


%	Etienne Balmes
%       Copyright (c) 1990-2020 by SDTools, All Rights Reserved.
%       $Revision: 1.3 $  $Date: 2020/02/26 08:47:33 $

error(nargchk(1,1,nargin));

sdtw('ofact objects may not load properly');
out=obj;
