function out=skyline(varargin);

%OBSOLETE NAME FOR OFACT OBJECT
%       Copyright (c) 1990-2009 by SDTools, All Rights Reserved.
%       
if nargout==1;
    if isequal(varargin{1},'cvs');
        out='$Revision: 1.6 $  $Date: 2009/05/28 16:42:00 $';
    else;
      warning('SKYLINE is obsolete, use OFACT instead');
      out=ofact(varargin{:});
    end
elseif nargout==0; 
    warning('SKYLINE is obsolete, use OFACT instead');
    ofact(varargin{:});
end