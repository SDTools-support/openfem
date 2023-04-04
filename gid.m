function out = gid(varargin)


% model=gid('Read FileName');
%

%	Etienne Balmes
%       Copyright (c) 2001-2023 by INRIA and SDTools
%       Use under OpenFEM trademark.html license and LGPL.txt library license
%       All Rights Reserved.


[CAM,Cam]=comstr(varargin{1},1);

if comstr(Cam,'read');[CAM,Cam]=comstr(CAM,5);
if comstr(Cam,'mts')
%% #ReadMts 
fin=fopen(varargin{2},'r');st=fscanf(fin,'%c');fclose(fin);
st=feval(sdtweb('@stringAsUx'),st);

%assignin('base','st',st)
%return
[st1,i2]=textscan(st,'%s',20,'delimiter','\n');%read lines
% st(1:i2+1)=''
st1=st1{1};ind=find(cellfun(@isempty,st1)==1,1,'first');

[st2,i2]=regexp(st1(1:ind-1),'^([^:]*):([^\n]*)','tokens');
st2=[vertcat(cellfun(@(x)x{1}{1},st2,'uni',0)) vertcat(cellfun(@(x)x{1}{2},st2,'uni',0))];
out=struct('header',{st2}); 
while isempty(st1{ind});ind=ind+1;end
st2=textscan(st1{ind},'%s','delimiter','\t'); out.ColumnName=st2{1}';ind=ind+1;
st2=textscan(st1{ind},'%s','delimiter','\t'); out.unit=st2{1}';ind=ind+1;
st(st==',')='.';
i2=cellfun(@length,st1(1:ind-1));i2=sum(i2)+length(i2);st(1:i2);
val=textscan(st(i2+1:end),repmat('%n',1,length(out.ColumnName)),'delimiter','\t');
out.X{1}=val{1};out.Xlab{1}={out.ColumnName{1},out.unit{1},[]};
out.X{2}=[out.ColumnName(2:end)' out.unit(2:end)'];out.Xlab{2}='comp';
out.Y=horzcat(val{2:end});
out=rmfield(out,'ColumnName');out=rmfield(out,'unit');
out.name=out.header{2,2};
if nargout==0; iicom('curveinit',out);end

else
%% #ReadGid

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
end
elseif comstr(Cam,'writemts');FileName=comstr(CAM,9); 
%% #WriteMts : gid('WriteMts FileName',C1,RunOpt) - - - - - - - - - - - - - -
fid = fopen(FileName,'w');
C1=varargin{2};
if ~isfield(C1,'header')
  if ~isfield(C1,'name');C1.name='curve';end
  C1.header={'FileType','Block-Arbitrary'
      'channels',num2str(size(C1.Y,2));'Description',C1.name;
      'Date',datestr(now,'yyyy-mm-dd HH:MM')
      'ActionList','Counter 1, Counter2'}; %#ok<TNOW1,DATST> 
end
if length(C1.X)<2; C1.X{2}={'Level 1','mm'};end

st=C1.header'; fprintf(fid,'%s:%s\n',st{:});
st=repmat('%s\t',1,size(C1.Y,2)+1);st(end)='n';
fprintf(fid,'\n');
fprintf(fid,st,C1.Xlab{1}{1},C1.X{2}{:,1});
fprintf(fid,st,C1.Xlab{1}{2},C1.X{2}{:,2});
st=repmat('%.8g\t',1,size(C1.Y,2)+1);st(end)='n';
fprintf(fid,st,[diff(C1.X{1}) C1.Y(1:end-1)]');
fprintf(fid,'\n');
fclose(fid);


%%  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'cvs') 
    out='$Revision: 1.10 $  $Date: 2023/03/20 17:42:03 $';
else; sdtw('''%s'' unknown',Cam);
end 
