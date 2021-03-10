function [out,out1]=nopo(varargin);

%NOPO imports nopo file (cf. Modulef)
%
%       Syntax: model      = nopo('read -v -p type FileName')
%               [Node,Elt] = nopo('read -v -p type FileName')
%       
%       '-v' is used for verbose output
%       '-p type' gives the type of problem described in the nopo
%       file. To allow proper translation to OpenFEM element names
%       TYPE can be:
%
%        '2D'               2D     
%        '3D'               3D  
%        'AXI'              Axisymmetric
%        'FOURIER'          Axisymmetric Fourier Expansion
%        'INCOMPRESSIBLE'   Incompressible
%        'PLAQUE'           Plate
%        'COQUE'            Shell
%                             
       
%       Author Frank GENOT 13/09/01, Etienne Balmes
%       Copyright (c) 2001-2009 by INRIA and SDTools
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       All Rights Reserved


if comstr(varargin{1},'cvs')
 out='$Revision: 1.13 $  $Date: 2009/05/28 16:42:00 $'; return;
end
[CAM,Cam]=comstr(varargin{1},1);carg=2;

if comstr(Cam,'read')

  [CAM,Cam]=comstr(CAM,5);
  if comstr(CAM,'-v'); [CAM,Cam]=comstr(CAM,3); verbose='v';
  else verbose='n';end
  if comstr(CAM,'-p'); [CAM,Cam]=comstr(CAM,3); 
   [problem,i2,i3,i4]=sscanf(CAM,'%s',1);[CAM,Cam]=comstr(CAM,i4);
   problem=comstr(problem,-271);
  else problem='3D'; end
  [CAM,Cam]=comstr(CAM,1);

  [wd,fname,ext]=fileparts(CAM);if isempty(ext);ext='.nopo';end
  fid=fopen(fullfile(wd,[fname ext]));
  if fid<0; error('File not found');
  else; fclose(fid); end

  [Node, Elt, redge, rface] = nopo2sd(fullfile(wd,fname), problem, verbose);
  rows = find(Elt(:,1) == -1); Elt(rows, 1) = Inf*ones(size(rows,1), 1);

  eltid=feutil('eltidfix',Elt);

  if nargout==2 out=Node; out1=Elt;
  else
   out=struct('Node',Node,'Elt',Elt);
   i1=unique(redge);i1=i1(find(i1));i2=unique(rface);i2=i2(find(i2));
   if ~isempty(i1)|~isempty(i2) out.Stack=cell(length(i1)+length(i2),3);end
   j0=1;
   for j1=1:length(i1)
    [i3,i4]=find(redge==i1(j1));
    out.Stack(j0,1:3)={'set',sprintf('Edge %i',i1(j1)), ...
     struct('ID',i1(j1),'data',[eltid(i3) i4])};
    j0=j0+1;
   end
   for j1=1:length(i2)
    [i3,i4]=find(rface==i2(j1));
    out.Stack(j0,1:3)={'set',sprintf('Face %i',i2(j1)), ...
     struct('ID',i2(j1),'data',[eltid(i3) i4])};
    j0=j0+1;
   end

  end


end
