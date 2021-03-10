function out = gid(varargin)


% model=gid('Read FileName');
%

%	Etienne Balmes
%       Copyright (c) 2001-2009 by INRIA and SDTools
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       All Rights Reserved.


[CAM,Cam]=comstr(varargin{1},1);

if comstr(Cam,'read')

FileName=comstr(CAM,5);

out.Node=[]; out.Elt=[];
fid = fopen(FileName); fseek(fid,0,-1);
flag=ffindstr(fid,'MESH','mesh');
if flag==0 || flag==1;  st=fgetl(fid); end % header

while(1)
 % - - - - - - - - - - - reading nodes
 if ~isempty(strfind(fgetl(fid),'oordinates'))
  node=fscanf(fid,'%i %f %f %f',[4 inf]);
  if ~isempty(strfind(st,'dimension 3'))||~isempty(strfind(st,'dimension = 3'))
   out.Node=[out.Node; node(1,:)' zeros(size(node,2),3) node(2:4,:)'];
  elseif ~isempty(strfind(st,'dimension 2')) || ...
         ~isempty(strfind(st,'dimension = 2'))
   out.Node=[out.Node;                            
             node(1,:)' zeros(size(node,2),3) node([3 4 2],:)']; 
  end
 else error('error on nodes')
 end
 if isempty(strfind(fgetl(fid),'end coordinates')) 
  error('error on end of nodes')
 end

 % - - - - - - - - - - - reading elements
 if ffindstr(fid,'lements')==0
  if ~isempty(strfind(st,'riangle'))
   % tria3 - - - - - - - - - - - - - - - - - - - - - - - - - 
   elt=fscanf(fid,'%i',[5 inf]);
   out.Elt(size(out.Elt,1)+1,1:5)=[Inf abs('tria3')];
   out.Elt(size(out.Elt,1)+[1:size(elt,2)],1:6)=elt([2:4 size(elt,1) size(elt,1) 1],:)';
  elseif ~isempty(strfind(st,'uadrilateral'))  
   % quad - - - - - - - - - - - - - - - - - - - - - - - - - 
   elt=fscanf(fid,'%i',[6 inf]);
   out.Elt(size(out.Elt,1)+1,1:6)=[Inf abs('quad4')];
   out.Elt(size(out.Elt,1)+[1:size(elt,2)],1:7)=elt([2:5 size(elt,1) size(elt,1) 1],:)';
  elseif ~isempty(strfind(st,'etrahedra'))  
   % tetra4  - - - - - - - - - - - - - - - - - - - - - - - - - 
   elt=fscanf(fid,'%i',[6 inf]);
   out.Elt(size(out.Elt,1)+1,1:7)=[Inf abs('tetra4')];
   out.Elt(size(out.Elt,1)+[1:size(elt,2)],1:7)=elt([2:5 size(elt,1) size(elt,1) 1],:)';
  elseif ~isempty(strfind(st,'exahedra'))  
   % hexa8  - - - - - - - - - - - - - - - - - - - - - - - - - 
   a1=ftell(fid); 
   fgetl(fid); st=fgetl(fid);
   i2=length(str2num(st));
   a2=ftell(fid);fseek(fid,-a2+a1,0); 

   elt=fscanf(fid,'%i',[i2 inf]);
   if i2==10
    out.Elt(size(out.Elt,1)+1,1:6)=[Inf abs('hexa8')];
    out.Elt(size(out.Elt,1)+[1:size(elt,2)],1:11)=elt([2:9 size(elt,1) size(elt,1) 1],:)';
   elseif i2==22
    out.Elt(size(out.Elt,1)+1,1:7)=[Inf abs('hexa20')];
    out.Elt(size(out.Elt,1)+[1:size(elt,2)],1:23)=...
                   elt([2:9 14:17 18:21 10:13 size(elt,1) size(elt,1) 1],:)';
   end
  elseif ~isempty(strfind(st,'risma'))  
   % penta6  - - - - - - - - - - - - - - - - - - - - - - - - - 
   elt=fscanf(fid,'%i',[8 inf]);
   out.Elt(size(out.Elt,1)+1,1:7)=[Inf abs('penta6')];
   out.Elt(size(out.Elt,1)+[1:size(elt,2)],1:9)=elt([2:7 size(elt,1) size(elt,1) 1],:)';
  else error('not supported element type')
  end
 else error('error on elements')
 end
  
 if isempty(strfind(fgetl(fid),'end elements'))
  error('error on end of elements')
 end

 flag= ffindstr(fid,'MESH','mesh');
 if flag==0 || flag==1; st=fgetl(fid);
 elseif flag==-1
  break
 end

end % end loop on MESH

fclose(fid);
if nargout==0; feplot(out); end

%  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'cvs'); 
    out='$Revision: 1.8 $  $Date: 2009/05/28 16:42:00 $';
else sdtw('''%s'' unknown',Cam);
end 