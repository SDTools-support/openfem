function out = saveobj(obj)

%SAVEOBJ save filter for OFACT objects


%	Etienne Balmes
%       Copyright (c) 1990-2005 by SDTools, All Rights Reserved.
%       $Revision: 1.3 $  $Date: 2007/05/21 16:25:16 $

error(nargchk(1,1,nargin));

sdtw(' %s not a supported ofact saveobj',obj.method.name);
r1='';

out=[];
