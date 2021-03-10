function out=of2vtk(varargin)

%  Convert openFEM format to VTK
%
%  To create a new file:
%     of2vtk('write',FileName,model,NodeData,EltData);
%
%  To append data to an existing file
%     of2vtk('append',FileName,model,NodeData,EltData);
%  The previous command does not copy model data, but it is mandatory to
%  provide it since a data check is performed:
%     - length of NodeData must be equal to the number of nodes
%     - length of EltData  must be equal to the number of elements
%
%  NodeData is the form of a deformation (with fields def and DOF is
%           accepted)

%       Arnaud Sternchuss, Etienne Balmes
%       Copyright (c) 2001-2020 by INRIA and SDTools, All rights reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license


CAM=varargin{1};carg=2;
[CAM,Cam]=comstr(CAM,1);
if comstr(Cam,'write'); RunOpt.Cmd='write';[CAM,Cam]=comstr(CAM,6);
elseif comstr(Cam,'append')
  RunOpt.DataType=comstr(CAM,7); RunOpt.Cmd='append';
elseif comstr(Cam,'cvs')
 out='$Revision: 1.15 $  $Date: 2020/02/26 08:47:15 $'; return;
elseif comstr(Cam,'const') % constants for paraview
  % CellType,SdtName,ParaviewNumbers in SDT
  % See vtkCellType.h for more details
  out={25,'hexa20',1:8;
       26,'penta15',[1:6 7:9 13:15 10:12]};
    
  return
else;RunOpt.Cmd='write';RunOpt.Name=CAM;
end

if ischar(varargin{carg});RunOpt.Name=varargin{carg};carg=carg+1;
elseif ~isempty(CAM); RunOpt.Name=CAM;
else
 st1={'*.vtk', 'VTK legacy file'};
 [RunOpt.Name,wd]=uiputfile(st1,'Write to VTK legacy format');
 RunOpt.Name=fullfile(wd,RunOpt.Name);carg=carg-1;
end
[wd,obj,ext]=fileparts(RunOpt.Name); if isempty(ext);ext='.vtk';end
if isempty(wd); wd=pwd;end

model=varargin{carg};carg=carg+1; % always load model
model=feutil('renumber',model); % renumber for compatibility

switch RunOpt.Cmd
case 'write'
 % NEW FILE
 fidend = fopen(fullfile(wd,[obj ext]),'w');

Nodes = model.Node(:,5:7);
elt = model.Elt;
[EGroup,nGroup]=getegroup(elt);
ig = find(elt(:,1)==Inf);
nGroup=length(ig);
ig=[ig;size(elt,1)+1];

 
  % HEADER  
  fprintf(fidend,['# vtk DataFile Version 2.0','\n']);
  fprintf(fidend,'%s\n',RunOpt.Name);
  fprintf(fidend,['ASCII','\n']);
  fprintf(fidend,['DATASET UNSTRUCTURED_GRID','\n']);
  
  % POINTS
  %disp('Write Points')
  fprintf(fidend,'POINTS ');
  fprintf(fidend,['%i',' '],size(Nodes,1));
  fprintf(fidend,['float','\n']);
  fprintf(fidend,'%13.6f %13.6f %13.6f\n',Nodes');
  
  
  cell_type=zeros(nGroup,1);
  cell_type=zeros(nGroup,1);
  ne=zeros(nGroup,1);
  elemF=cell(nGroup,1);
  T=cell(nGroup,1);
  n_cell=0;
  size_cell = 0;
  % ELEMENTS IDENTIFICATION
  for jGroup=1:nGroup
    [ElemF,opt,ElemP]= feutil('getelemf',elt(EGroup(jGroup),:),jGroup);
    k=jGroup; elemF{jGroup} = ElemP;
    ne(jGroup)=size(elt(ig(jGroup)+1:ig(jGroup+1)-1,:),1);
    n_cell = n_cell+ne(jGroup);
    if (comstr(ElemP,'tria3') || comstr(ElemP,'t3p') || comstr(ElemP,'dktp'))
        size_cell = size_cell + ne(jGroup)*4;
        T{jGroup} = [3*ones(ne(jGroup),1) elt(ig(jGroup)+1:ig(jGroup+1)-1,1:3)-1];
        cell_type(jGroup)=5;
    elseif (comstr(ElemP,'beam1')  || comstr(ElemP,'bar1'))
        size_cell = size_cell + ne(jGroup)*3;
        T{jGroup} = [2*ones(ne(jGroup),1) elt(ig(jGroup)+1:ig(jGroup+1)-1,1:2)-1];
        cell_type(jGroup)=3;        
    elseif (comstr(ElemP,'quad4') || comstr(ElemP,'q4p') || comstr(ElemP,'mitc4'))
        size_cell = size_cell + ne(jGroup)*5;
        T{jGroup} = [4*ones(ne(jGroup),1) elt(ig(jGroup)+1:ig(jGroup+1)-1,1:4)-1];
        cell_type(jGroup)=9;
    elseif comstr(ElemP,'tetra4')
        size_cell = size_cell + ne(jGroup)*5;
        T{jGroup} = [4*ones(ne(jGroup),1) elt(ig(jGroup)+1:ig(jGroup+1)-1,1:4)-1];
        cell_type(jGroup)=10;
    elseif comstr(ElemP,'hexa8')
        size_cell = size_cell + ne(jGroup)*9;
        T{jGroup} = [8*ones(ne(jGroup),1) elt(ig(jGroup)+1:ig(jGroup+1)-1,1:8)-1];
        cell_type(jGroup)=12;
    elseif comstr(ElemP,'penta6')
        size_cell = size_cell + ne(jGroup)*7;
        T{jGroup} = [6*ones(ne(jGroup),1) elt(ig(jGroup)+1:ig(jGroup+1)-1,1:6)-1];
        cell_type(jGroup)=13;
    elseif comstr(ElemP,'mitc6')||comstr(ElemP,'tria6')
        size_cell = size_cell + ne(jGroup)*7;
        T{jGroup} = [6*ones(ne(jGroup),1) elt(ig(jGroup)+1:ig(jGroup+1)-1,1:6)-1];
        cell_type(jGroup)=22;
    elseif (comstr(ElemP,'q8p') || comstr(ElemP,'quadb'))
        size_cell = size_cell + ne(jGroup)*9;
        T{jGroup} = [8*ones(ne(jGroup),1) elt(ig(jGroup)+1:ig(jGroup+1)-1,1:8)-1];
        cell_type(jGroup)=23;
    elseif comstr(ElemP,'tetra10')
        size_cell = size_cell + ne(jGroup)*11;
        T{jGroup} = [10*ones(ne(jGroup),1) elt(ig(jGroup)+1:ig(jGroup+1)-1,1:10)-1];
        cell_type(jGroup)=24;
     elseif comstr(ElemP,'penta15')
         size_cell = size_cell + ne(jGroup)*16;
         T{jGroup} = [15*ones(ne(jGroup),1) ...
         elt(ig(jGroup)+1:ig(jGroup+1)-1,1:9)-1 ...
         elt(ig(jGroup)+1:ig(jGroup+1)-1,13:15)-1 ...
         elt(ig(jGroup)+1:ig(jGroup+1)-1,10:12)-1];
         cell_type(jGroup)=26;
    elseif comstr(ElemP,'hexa20')
        size_cell = size_cell + ne(jGroup)*21;
        T{jGroup} = [20*ones(ne(jGroup),1) elt(ig(jGroup)+1:ig(jGroup+1)-1,1:12)-1 elt(ig(jGroup)+1:ig(jGroup+1)-1,17:20)-1 elt(ig(jGroup)+1:ig(jGroup+1)-1,13:16)-1];
        cell_type(jGroup)=25;
    else
        error('Sorry, %s not yet implemented :-) !',ElemP)
    end
  end
  
  ref = [];
  fprintf(fidend,'CELLS ');
  fprintf(fidend,['%i',' '],n_cell);
  fprintf(fidend,['%i','\n'],size_cell);
  
  % ELEMENTS WRITING
  for jGroup=1:nGroup
	if cell_type(jGroup)==5
		fprintf(fidend,'%i %6i %6i %6i\n',T{jGroup}');
	elseif ( cell_type(jGroup)==3 )
		fprintf(fidend,'%i %6i %6i\n',T{jGroup}');        
	elseif (cell_type(jGroup)==9 ||cell_type(jGroup)==10)
		fprintf(fidend,'%i %6i %6i %6i %6i\n',T{jGroup}');
  	elseif (cell_type(jGroup)==12 || cell_type(jGroup)==23)
		fprintf(fidend,'%i %6i %6i %6i %6i %6i %6i %6i %6i\n',T{jGroup}');
        elseif (cell_type(jGroup)==13 || cell_type(jGroup)==22)
                fprintf(fidend,'%i %6i %6i %6i %6i %6i %6i\n',T{jGroup}');
  	elseif cell_type(jGroup)==24
		fprintf(fidend,'%i %6i %6i %6i %6i %6i %6i %6i %6i %6i %6i\n',T{jGroup}');
   	elseif cell_type(jGroup)==26
		fprintf(fidend,'%i %6i %6i %6i %6i %6i %6i %6i %6i %6i %6i %6i %6i %6i %6i %6i\n',T{jGroup}');
  	elseif cell_type(jGroup)==25
		fprintf(fidend,'%i %6i %6i %6i %6i %6i %6i %6i %6i %6i %6i %6i %6i %6i %6i %6i %6i %6i %6i %6i %6i\n',T{jGroup}');
	end
  end	

  fprintf(fidend,'CELL_TYPES ');
  fprintf(fidend,'%i\n',n_cell);
  
  for jGroup=1:nGroup
    L = cell_type(jGroup)*ones(1,ne(jGroup));
    %disp('Write Cells types')
    fprintf(fidend,'%i\n', L');
  end

case 'append'
 % APPEND TO AN EXISTING FILE 
 % look for existing stacks
 fidend = fopen(fullfile(wd,[obj ext]),'a+');
  i1=[];i2=[];lnum=0;
   while feof(fidend)==0
    lnum=lnum+1;
    curline=fgetl(fidend);
    switch RunOpt.DataType
    case 'pointdata'
    if ~isempty(strfind(curline,'POINT_DATA')); i1=line; break; end
    case 'celldata'
    if ~isempty(strfind(curline,'CELL_DATA'));  i2=line; break; end
    end
   end 
otherwise; sdtw('%s is not a supported command',RunOpt.Cmd);
end

if carg<=nargin % DEAL WITH OTHER ENTRIES
    nodevar =[]; eltvar=[];
    while carg<=nargin
        r1=varargin{carg};carg=carg+1;
        if isfield(r1,'def')&&isfield(r1,'DOF')
         r2=double(stack_get(model,'info','OrigNumbering','getdata'));
         % as many rows as nodes in the model
         if sp_util('issdt');
          r2=[r2(:,1)'+.01;r2(:,1)'+.02;r2(:,1)'+.03];
          eval('r1=feutilb(''placeindof'',r2(:),r1);');
          r1.def=reshape(r1.def,3,[],size(r1.def,2));
         else
          r1.def = feutil('dof2mode',r1.def,r1.DOF);
          [i1,i2]=ismember(r2(:,1),r1.def(:,1,1));
          r2=zeros(3,size(i2,1),size(r1.def,3));
          r2(:,i1,:)=permute(r1.def(i2(i1),2:4,:),[2 1 3]);
         end

         if ~isempty(strfind(RunOpt.Cmd,'write'))||isempty(i1)
         % write file, one stack per data
         fprintf(fidend,'\nPOINT_DATA ');
         fprintf(fidend,'%i\n',size(r1.def,2));
         end
         for j1=1:size(r1.def,3)
           if isfield(r1,'lab')&&~isempty(r1.lab{j1}) % name already exists
            st=r1.lab{j1};
            st(isspace(st))='_'; % vtk does not like blank spaces
           elseif isfield(r1,'data')&&j1<size(r1.data,1);
            st=sprintf('Mode%i_@_%.2f_Hz',j1,r1.data(j1,1));
           else;st=sprintf('%i ',j1);
           end
           fprintf(fidend,['\nVECTORS ' st ' float\n']);
           fprintf(fidend,'%e %e %e\n',r1.def(:,:,j1));
         end 

        elseif isfield(r1,'data')&&isfield(r1,'EltId')
         if ~isempty(strfind(RunOpt.Cmd,'write'))||isempty(i2)
         % write file, one stack per data
         fprintf(fidend,'\nCELL_DATA ');
         fprintf(fidend,'%i\n',size(r1.data,1));
         end
         for j1=1:size(r1.data,2)
           if isfield(r1,'lab')&&~isempty(r1.lab{j1}) % name already exists
            st=r1.lab{j1};
            st(isspace(st))='_'; % vtk does not like blank spaces
           else
            st=sprintf('Data%i',j1);
           end
           fprintf(fidend,['\nSCALARS ' st ' float\n']);
           fprintf(fidend,'LOOKUP_TABLE default\n');
           fprintf(fidend,'%e\n',r1.data(:,j1));     
         end 
         
       elseif size(r1.Data,1) == size(Nodes,1)
         nodevar = [nodevar carg-1];
       elseif size(r1.Data,1) == n_cell
         eltvar = [eltvar carg-1];
       end
    end
    if ~isempty(nodevar);writeVTKnodedata(fidend,varargin{nodevar});end
    if ~isempty(eltvar);writeVTKeltdata(fidend,varargin{eltvar});end 

end

fclose('all');

% ----------------------------------------------------------------

function writeVTKnodedata(varargin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    append data at node to an existant vtk file
%
%    io_writeVTKnodedata(vtkrootfile,var1,var2,...varn)
%
%    with var_i as follow:
%
%            var_i.Name = 'namegiventodata'
%            var_i.Data = datavalues   
%                    nnode rows & one column for scalars
%                    nnode rows & three columns for vectors 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%rootMeshFile = varargin{1};

  % APPEND FILE
  %obj = [rootMeshFile,'.vtk'];
  %fidend = fopen(obj,'a');
  fidend = varargin{1};

 fprintf(fidend,'\nPOINT_DATA ');
 fprintf(fidend,'%i\n',size(varargin{2}.Data,1));


for k=2:nargin
  if size(varargin{k}.Data,2)>1
     fprintf(fidend,'\nVECTORS %s float\n',varargin{k}.Name);
     fprintf(fidend,'%e %e %e\n',varargin{k}.Data');
  else
     fprintf(fidend,'SCALARS %s float\n',varargin{k}.Name);
     fprintf(fidend,'LOOKUP_TABLE default\n');
     fprintf(fidend,'%e\n',varargin{k}.Data');
  end
end

%fclose('all');


% -------------------------------------------------------------------

function writeVTKeltdata(varargin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    append data at node to an existant vtk file
%
%    io_writeVTKeltdata(vtkrootfile,var1,var2,...varn)
%
%    with var_i as follow:
%
%            var_i.Name = 'namegiventodata'
%            var_i.Data = datavalues   
%                    nelt rows & one column for scalars
%                    nelt rows & three columns for vectors 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%rootMeshFile = varargin{1};

  % APPEND FILE
  %obj = [rootMeshFile,'.vtk'];
  %fidend = fopen(obj,'a');
  fidend = varargin{1};


 fprintf(fidend,'\nCELL_DATA ');
 fprintf(fidend,'%i\n',size(varargin{2}.Data,1));


for k=2:nargin
  if size(varargin{k}.Data,2)>1
     fprintf(fidend,'\nVECTORS %s float\n',varargin{k}.Name);
     fprintf(fidend,'%e %e %e\n',varargin{k}.Data');
  else
     fprintf(fidend,'SCALARS %s float\n',varargin{k}.Name);
     fprintf(fidend,'LOOKUP_TABLE default\n');
     fprintf(fidend,'%e\n',varargin{k}.Data');
  end
end

%fclose('all');
