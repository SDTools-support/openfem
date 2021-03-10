function B = subsref(h,S)

% ofact/subsref

%       Etienne Balmes  07/28/92, 08/22/99
%       Copyright (c) 2001-2002 by INRIA and SDTools
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       All Rights Reserved.
%       $Revision: 1.3 $  $Date: 2004/08/17 09:15:39 $

switch S(1).type
case '.'
   B = eval(sprintf('h.%s',S(1).subs));
   if length(S)>1 B = subsref(B,S(2:end)); end
otherwise
   error('Not a valid ofact object reference');
end
