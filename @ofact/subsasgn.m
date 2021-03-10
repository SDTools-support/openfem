function obj = subsasgn(obj,S,B)

% ofact/subsasgn
%

%       Etienne Balmes  28/07/92, 14/03/01
%       Copyright (c) 2001-2003 by INRIA and SDTools,All Rights Reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       $Revision: 1.2 $  $Date: 2004/08/17 09:15:39 $

 switch S.type
 case '.'
  eval(sprintf('obj.%s=B;',S.subs));
 end