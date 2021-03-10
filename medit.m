function [indnum,scale] = medit(CAM,model,varargin)

% THIS FUNCTION DOES NOT PASS ITS TEST_MEDIT. ANYBODY IS WELCOME TO SEND
% FIXES SINCE NOBODY IN THE DEVELOPMENT TEAM CURRENTLY USES IT

%============================================================================%
%
% medit('write CAM',model)                        :     displays the mesh defined by model
%
% medit('write CAM',model,def,[opt])              :     displays the mesh defined by model in a window and the deformation definded by def in an other window
%
% medit('write CAM',model,def,'a',[opt])          :     animates the deformations defined by def on model
%
% medit('write CAM',model,[],strain)              :     displays the mesh defined by model and colors it with the help of strain
%
% medit('write CAM',model,def,strain,[opt])       :     displays the mesh defined model in a window and the deformation defined by def with colors due to strain in an other window
%
% medit('write CAM',model,def,strain,'a',[opt])   :     animates the deformations defined by def on model and colors them with the help of strain
%
%===========================================================================%
%
% CAM : filename where information for Medit will be written, no extension has to be given
% model : a structure defining a model. It must contain at least fields .Node and .Elt
% def : a structure defining deformations that users want to visualize. It must contain at least fields .def and .DOF
% strain : a structure defining deformations coloring. It must contain at least fields .data and .DOF (case .EltInd is not treated, use feplot)
% opt : an option vector. opt = [numdef nbimag scale] where :
%		- numdef is the number of the mode to display
%		- nbimag is the number of files to create the animation of deformations
%		- scale is the display scale (a parameter for increasing the deformations
% 
% The default system call 'medit' can be modified. For example on windows
% setpref('OpenFEM','MeditFcn','medit-2.3-win')

%============================================================================%
%       Author Claire Delforge, Etienne Balmes
%       Copyright (c) 2003-2005 by INRIA and SDTools, All Rights Reserved
%       Use under OpenFEM trademark.html license and LGPL.txt library license

if comstr(CAM,'cvs')
 indnum='$Revision: 1.14 $  $Date: 2019/01/15 18:11:29 $'; return;
end

indnum = []; scale = [];
MeditFcn=getpref('OpenFEM','MeditFcn','medit');MeditFcn(end+1)=' ';

if ~comstr(CAM,'write'); error('medit : first argument must be a character array beginning by ''write'''); end
if nargin<2; error('medit : at less 2 inputs are required'); end

[CAM,Cam] = comstr(CAM,6);

opt = [];
carg = 1;

if nargin==2;
	indic_case = 1;
else
	def = varargin{carg}; carg = carg+1;
	if nargin==3;
		indic_case = 2;
	else
		if ischar(varargin{carg})
			indic_case = 3;
			indic_anim = varargin{carg}; carg = carg+1;
			if nargin==5; opt = varargin{carg}; carg = carg+1; end
		elseif isempty(def)
			indic_case = 4;
			strain = varargin{carg}; carg = carg+1;
		elseif isstruct(varargin{carg}) & nargin==4
			indic_case = 5;
			strain = varargin{carg}; carg = carg+1;
		elseif nargin==4
			indic_case = 2;
			opt = varargin{carg}; carg = carg+1;
		elseif nargin==5 & ~ischar(varargin{3})
			indic_case = 5;
			strain = varargin{carg}; carg = carg+1;
			opt = varargin{carg}; carg = carg+1; 
		else
			indic_case = 6;
			strain = varargin{carg}; carg = carg+1;
			indic_anim = varargin{carg}; carg = carg+1;
			if ~comstr(indic_anim,'a'); error('medit : bad inputs'); end
			if nargin==6; opt = varargin{carg}; end
		end
	end
end


%===========================%
% get the options parameter %
%===========================%
if ~isempty(opt); numdef = opt(1);
else numdef = 1; end
if length(opt)>=2; nbimag = opt(2);
else nbimag = 20; end 
if length(opt)>=3; scale = opt(3); 
else scale = 10; end

%================%
% DISPLAY A MESH %
%================%
if indic_case==1  
	filename = strcat([CAM '.mesh']);
	indnum = writemesh(filename,model);
	system(strcat([MeditFcn filename]));
	scale = 0;

%==================================%
% DISPLAY A MESH AND A DEFORMATION %
%==================================%
elseif indic_case==2
	if size(def.def,2)>1; def.def = def.def(:,numdef); end 
	[indnum,scale] = writemesh_and_deformed(CAM,model,def,scale);
	disp('you can link the 2 medit windows : type <ALT> c in a window, <ALT> l in the other and move the structure of the first window.');
	system(strcat([MeditFcn CAM 'deformed.mesh ' CAM 'original.mesh']));

%=======================%
% ANIMATE A DEFORMATION %
%=======================%
elseif indic_case==3 
	if size(def.def,2)>1; def.def = def.def(:,numdef); end 
	[indnum,scale] = animate_deformed(CAM,model,def,scale,nbimag);
	system(strcat([MeditFcn CAM  ' -a 1 ' num2str(nbimag)]));

%=================================%
% DISPLAY A MESH WITH CONSTRAINTS %
%=================================%
elseif indic_case==4
	writemesh_and_constraint(CAM,model,strain);
	system(strcat([MeditFcn CAM]));

%=========================================%
% DISPLAY A MESH + DEFORMED + CONSTRAINTS %
%=========================================%
elseif indic_case==5
	if size(def.def,2)>1; def.def = def.def(:,numdef); end 
	[indnum,scale] = writemesh_and_defconst(CAM,model,def,scale,strain);
	system(strcat([MeditFcn CAM 'deformed.mesh ' CAM 'original.mesh']));

%===================================%
% ANIMATE DEFORMED WITH CONSTRAINTS %
%===================================%
elseif indic_case==6
	if size(def.def,2)>1; def.def = def.def(:,numdef); end 
	if size(strain.data,2)>1; strain.data = strain.data(:,numdef); end 
	[indnum,scale] = animate_deformed_constraint(CAM,model,def,scale,nbimag,strain);
	system(strcat([MeditFcn CAM  ' -a 1 ' num2str(nbimag)]));

end




%---------------------------------------------------------------------------------------------------------------------------------------------------%
%							writemesh										     %
%---------------------------------------------------------------------------------------------------------------------------------------------------%
function out = writemesh(filename,model)
%==============================================================%
% 			displays a mesh 			%
%==============================================================%
% filename : filename where information for Medit will be written (with extension)
% model : a structure defining a model. It must contain at least fields .Node and .Elt
% out : vector that gives the corresponding between nodes numerotation in model.Node and in medit
%       out(model.Node(:,1)) = numerotation in medit
% the vector that gives the corresponding between medit and model.Node is not given because it is obvious (model.Node(:,1))

out = [];

[EGroup,nGroup] = getegroup(model.Elt);

matquad = [];
mattri = [];

%========================%
% builds rface and redge %  for the moment, only rface
%========================%
% determination of rface and dimensions
dimrface = 0;
for i=1:nGroup
	ElemF{i} = feutil('getelemf',model.Elt(EGroup(i),:));
	tmpdimrface = size(fe_super('face',ElemF{i}),1);
	dimrface = max(dimrface,tmpdimrface);
end
rface = zeros(size(model.Elt,1),dimrface);
% building the out vector
eltid = feutil('eltidfix',model.Elt);
coreltid = zeros(1,max(eltid));
tmpeltid = find(eltid);
coreltid(eltid(tmpeltid)) = tmpeltid;
% building rface
if isfield(model,'Stack')
	for i=1:size(model.Stack,1)
		if comstr(model.Stack{i,1},'set') & comstr(model.Stack{i,2},'Face')
			indexrface = (model.Stack{i,3}.data(:,2)-1)*size(rface,1)+coreltid(model.Stack{i,3}.data(:,1))';
			rface(indexrface) = model.Stack{i,3}.ID;
		end
	end
end

%==============================%
% renumerotation of model.Node %
%==============================%
tmpindnode = model.Node(:,1);
model.Node(:,1) = [1:size(model.Node,1)]';
out = zeros(max(tmpindnode),1);
out(tmpindnode) = model.Node(:,1);


%% node writing %
[fd,err] = fopen(filename,'wb');
fprintf(fd,'MeshVersionFormatted 1\n');
fprintf(fd,'Dimension 3\n');

fprintf(fd,'Vertices %d\n',size(model.Node,1));
for i=1:size(model.Node,1)
	fprintf(fd,'%f %f %f %d\n',model.Node(i,5),model.Node(i,6),model.Node(i,7),model.Node(i,4));
end

addmatquad = [];
addmattri = [];
%================================%
% building of matquad and mattri %
%================================%
for jGroup=1:nGroup % loop on groups
	elt = model.Elt(EGroup(jGroup):EGroup(jGroup+1)-1,:);
	
	curface = fe_super('sci_face',ElemF{jGroup});
	curnode = fe_super('node',ElemF{jGroup});
	curfaceref = fe_super('face',ElemF{jGroup});
	if size(curfaceref,2)<=4; indic_pb1 = 1;
	else indic_pb1 = 0;
	end
	indquad = [];
	indtri  = [];

	% modifications of elt due to the renumerotation of model.Node
	nbnode = length(curnode);
	elt(2:end,1:nbnode) = reshape(out(elt(2:end,1:nbnode)),size(elt,1)-1,nbnode);

	for i=1:size(curface,1)	
		diffcurface = curface(i,:)-[curface(i,2:end) curface(i,1)];
		finddiffcur = find(diffcurface);
		if length(finddiffcur)==4; indquad(size(indquad,1)+1,1:4) = curface(i,finddiffcur); 
		elseif length(finddiffcur)==3
			indtri(size(indtri,1)+1,1:3) = curface(i,finddiffcur);
		elseif length(finddiffcur)==2
			indtri(size(indtri,1)+1,1:2) = curface(i,finddiffcur);
			indtri(size(indtri,1),3) = indtri(size(indtri,1),2);
		elseif length(finddiffcur)==1
			indtri(size(indtri,1)+1,1) = curface(i,finddiffcur);
			indtri(size(indtri,1),2) = indtri(size(indtri,1),1);
			indtri(size(indtri,1),3) = indtri(size(indtri,1),2);
		end
	end
	sizequadperelt = size(indquad,1);
	sizetriperelt  = size(indtri,1);
	
	for i=2:size(elt,1) % loop on elements
		addmatquad = reshape(elt(i,indquad),sizequadperelt,4);
		addmattri  = reshape(elt(i,indtri),sizetriperelt,3);
		addmatquad = [addmatquad zeros(size(addmatquad,1),1)];
		addmattri = [addmattri zeros(size(addmattri,1),1)];

		
		findrface = find(rface(EGroup(jGroup)+i-1,:));
		valrface = rface(EGroup(jGroup)+i-1,findrface);
		if findrface
			for j=1:length(findrface)
				goodface = elt(i,curfaceref(findrface(j),:)); % face referencee dans la numerotation "globale"
				goodface = unique(goodface); goodface = goodface(length(goodface):-1:1);

				%==============%
				% CASE P1 - Q1 %
				%==============%
				if indic_pb1
					% quads
					if length(goodface)==4
					tmpaddmatquad = sort(addmatquad(:,1:4),2); tmpaddmatquad = tmpaddmatquad(:,4:-1:1);
					tmpgood = sum(tmpaddmatquad(:,1:4)==ones(size(addmatquad,1),1)*goodface,2);
					tmpgood2 = find(tmpgood==4);
					addmatquad(tmpgood2,5) = valrface(j);

					% triangles
					elseif length(goodface)==3
						tmpaddmattri = sort(addmattri(:,1:3),2); tmpaddmattri = tmpaddmattri(:,3:-1:1);
						tmpgood = sum(tmpaddmattri(:,1:3)==ones(size(addmattri,1),1)*goodface,2);
						tmpgood2 = find(tmpgood==3);
						addmattri(tmpgood2,4) = valrface(j);
					end

				else
				%==============%
				% CASE P2 - Q2 %
				%==============%
				% length goodface > 4 !
					for k=1:size(addmatquad,1)
						if length(intersect(addmatquad(k,1:4),goodface))==4; addmatquad(k,5) = valrface(j);
						end
					end
					for k=1:size(addmattri,1)
						if length(intersect(addmattri(k,1:3),goodface))==3; addmattri(k,4) = valrface(j);
						end
					end

				end % case p1-q1
			end
		end
		
		% matquad and mattri
		matquad(size(matquad,1)+1:size(matquad,1)+sizequadperelt,1:5) = addmatquad;
		mattri(size(mattri,1)+1:size(mattri,1)+sizetriperelt,1:4)  = addmattri;
	end % loop on elements
end % loop on groups 

% internal faces are removed
if ~isempty(matquad)
    tmpmatquad = sort(matquad(:,1:4),2);
    [tmpmatquad,indtmpmatquad] = sortrows(tmpmatquad);
    tmpmatquad = tmpmatquad - [tmpmatquad(2:end,:);tmpmatquad(1,:)];
    tmpmatquad = sum(abs(tmpmatquad(:,1:4)),2);
    tmpmatquad(find(~tmpmatquad)+1) = 0;
    indmatquad = indtmpmatquad(find(tmpmatquad));
    matquad = matquad(indmatquad,:);
end
    
if ~isempty(mattri)    
    tmpmattri = sort(mattri(:,1:3),2);
    [tmpmattri,indtmpmattri] = sortrows(tmpmattri);
    tmpmattri = tmpmattri - [tmpmattri(2:end,:);tmpmattri(1,:)];
    tmpmattri = sum(abs(tmpmattri(:,1:3)),2);
    tmpmattri(find(~tmpmattri)+1) = 0;
    indmattri = indtmpmattri(find(tmpmattri));
    mattri  = mattri(indmattri,:);
end
    

%==================%
% elements writing %
%==================%
if ~isempty(matquad)
    fprintf(fd,'Quadrilaterals %d\n',size(matquad,1));
    for i=1:size(matquad,1)
	    fprintf(fd,'%d %d %d %d %d\n',matquad(i,1),matquad(i,2),matquad(i,3),matquad(i,4),matquad(i,5));
    end
end
   
if ~isempty(mattri)
    fprintf(fd,'Triangles %d\n',size(mattri,1));
    for i=1:size(mattri,1)
	    fprintf(fd,'%d %d %d %d\n',mattri(i,1),mattri(i,2),mattri(i,3),mattri(i,4));
    end
end
    
fclose(fd);




%---------------------------------------------------------------------------------------------------------------------------------------------------%
%							writemesh_and_deformed									     %
%---------------------------------------------------------------------------------------------------------------------------------------------------%
function [out,scale] = writemesh_and_deformed(filename,model,def,scale)
%=======================================================================================================%
% displays the mesh defined by model in a window and the deformation definded by def in an other window %
%=======================================================================================================%
% filename : filename where information for Medit will be written (without extension)
% model : a structure defining a model. It must contain at least fields .Node and .Elt
% def : a structure defining deformations with at least fields .def and .DOF
% scale : scale 
% out : vector that gives the corresponding between nodes numerotation in model.Node and in medit
%       out(model.Node(:,1)) = numerotation in medit
% the vector that gives the corresponding between medit and model.Node is not given because it is obvious (model.Node(:,1))

out = [];
filename1 = strcat([filename 'original.mesh']);
filename2 = strcat([filename 'deformed.mesh']);

[EGroup,nGroup] = getegroup(model.Elt);

matquad = [];
mattri = [];

%========================%
% builds rface and redge %  for the moment, only rface
%========================%
% determination of rface and dimensions
dimrface = 0;
for i=1:nGroup
	ElemF{i} = feutil('getelemf',model.Elt(EGroup(i),:));
	tmpdimrface = size(fe_super('face',ElemF{i}),1);
	dimrface = max(dimrface,tmpdimrface);
end
rface = zeros(size(model.Elt,1),dimrface);
 
eltid = feutil('eltidfix',model.Elt);
coreltid = zeros(1,max(eltid));
tmpeltid = find(eltid);
coreltid(eltid(tmpeltid)) = tmpeltid;
% construction de rface
if isfield(model,'Stack')
	for i=1:size(model.Stack,1)
		if comstr(model.Stack{i,1},'set') & comstr(model.Stack{i,2},'Face')
			indexrface = (model.Stack{i,3}.data(:,2)-1)*size(rface,1)+coreltid(model.Stack{i,3}.data(:,1))';
			rface(indexrface) = model.Stack{i,3}.ID;
		end
	end
end

%==============================%
% renumerotation of model.Node %
%==============================%
tmpindnode = model.Node(:,1);
model.Node(:,1) = [1:size(model.Node,1)]';
out = zeros(max(tmpindnode),1);
out(tmpindnode) = model.Node(:,1);

addmatquad = [];
addmattri = [];
%================================%
% building of matquad and mattri %
%================================%
for jGroup=1:nGroup % loop on groups
	elt = model.Elt(EGroup(jGroup):EGroup(jGroup+1)-1,:);
	
	curface = fe_super('sci_face',ElemF{jGroup});
	curnode = fe_super('node',ElemF{jGroup});
	curfaceref = fe_super('face',ElemF{jGroup});
	if size(curfaceref,2)<=4; indic_pb1 = 1;
	else indic_pb1 = 0;
	end
	indquad = [];
	indtri  = [];

	% modifications of elt due to the renumerotation of model.Node
	nbnode = length(curnode);
	elt(2:end,1:nbnode) = reshape(out(elt(2:end,1:nbnode)),size(elt,1)-1,nbnode);

	for i=1:size(curface,1)	
		diffcurface = curface(i,:)-[curface(i,2:end) curface(i,1)];
		finddiffcur = find(diffcurface);
		if length(finddiffcur)==4; indquad(size(indquad,1)+1,1:4) = curface(i,finddiffcur); 
		elseif length(finddiffcur)==3
			indtri(size(indtri,1)+1,1:3) = curface(i,finddiffcur);
		elseif length(finddiffcur)==2
			indtri(size(indtri,1)+1,1:2) = curface(i,finddiffcur);
			indtri(size(indtri,1),3) = indtri(size(indtri,1),2);
		elseif length(finddiffcur)==1
			indtri(size(indtri,1)+1,1) = curface(i,finddiffcur);
			indtri(size(indtri,1),2) = indtri(size(indtri,1),1);
			indtri(size(indtri,1),3) = indtri(size(indtri,1),2);
		end
	end
	sizequadperelt = size(indquad,1);
	sizetriperelt  = size(indtri,1);
	
	for i=2:size(elt,1) % loop on elements
		addmatquad = reshape(elt(i,indquad),sizequadperelt,4);
		addmattri  = reshape(elt(i,indtri),sizetriperelt,3);
		addmatquad = [addmatquad zeros(size(addmatquad,1),1)];
		addmattri = [addmattri zeros(size(addmattri,1),1)];

		
		findrface = find(rface(EGroup(jGroup)+i-1,:));
		valrface = rface(EGroup(jGroup)+i-1,findrface);
		if findrface
			for j=1:length(findrface)
				goodface = elt(i,curfaceref(findrface(j),:)); 
				goodface = unique(goodface); goodface = goodface(length(goodface):-1:1);

				%==============%
				% CASE P1 - Q1 %
				%==============%
				if indic_pb1
					% quads
					if length(goodface)==4
					tmpaddmatquad = sort(addmatquad,2); tmpaddmatquad = tmpaddmatquad(:,4:-1:1);
					tmpgood = sum(tmpaddmatquad(:,1:4)==ones(size(addmatquad,1),1)*goodface,2);
					tmpgood2 = find(tmpgood==4);
					addmatquad(tmpgood2,5) = valrface(j);

					% triangles
					elseif length(goodface)==3
						tmpaddmattri = sort(addmattri,2); tmpaddmattri = tmpaddmattri(:,3:-1:1);
						tmpgood = sum(tmpaddmattri(:,1:3)==ones(size(addmattri,1),1)*goodface,2);
						tmpgood2 = find(tmpgood==3);
						addmattri(tmpgood2,4) = valrface(j);
					end

				else
				%==============%
				% CASE P2 - Q2 %
				%==============%
				% length goodface > 4 !
					for k=1:size(addmatquad,1)
						if length(intersect(addmatquad(k,1:4),goodface))==4; addmatquad(k,5) = valrface(j);
						end
					end
					for k=1:size(addmattri,1)
						if length(intersect(addmattri(k,1:3),goodface))==3; addmattri(k,4) = valrface(j);
						end
					end

				end % case p1-q1
			end
		end
		
		% matquad and mattri
		matquad(size(matquad,1)+1:size(matquad,1)+sizequadperelt,1:5) = addmatquad;
		mattri(size(mattri,1)+1:size(mattri,1)+sizetriperelt,1:4)  = addmattri;
	end % loop on elements
end % loop on groups 

% internal faces are removed
if ~isempty(matquad)
	tmpmatquad = sort(matquad(:,1:4),2);
	[tmpmatquad,indtmpmatquad] = sortrows(tmpmatquad);
	tmpmatquad = tmpmatquad - [tmpmatquad(2:end,:);tmpmatquad(1,:)];
	tmpmatquad = sum(abs(tmpmatquad(:,1:4)),2);
	tmpmatquad(find(~tmpmatquad)+1) = 0;
	indmatquad = indtmpmatquad(find(tmpmatquad));
	matquad = matquad(indmatquad,:);
end

if ~isempty(mattri)
    tmpmattri = sort(mattri(:,1:3),2);
    [tmpmattri,indtmpmattri] = sortrows(tmpmattri);
    tmpmattri = tmpmattri - [tmpmattri(2:end,:);tmpmattri(1,:)];
    tmpmattri = sum(abs(tmpmattri(:,1:3)),2);
    tmpmattri(find(~tmpmattri)+1) = 0;
    indmattri = indtmpmattri(find(tmpmattri));
    mattri  = mattri(indmattri,:);
end

%===============%
% original mesh %
%===============%

%===============%
% nodes writing %
%===============%
[fd,err] = fopen(filename1,'wb');
fprintf(fd,'MeshVersionFormatted 1\n');
fprintf(fd,'Dimension 3\n');

fprintf(fd,'Vertices %d\n',size(model.Node,1));
for i=1:size(model.Node,1)
	fprintf(fd,'%f %f %f %d\n',model.Node(i,5),model.Node(i,6),model.Node(i,7),model.Node(i,4));
end

%==================%
% elements writing %
%==================%
if ~isempty(matquad)
    fprintf(fd,'Quadrilaterals %d\n',size(matquad,1));
    for i=1:size(matquad,1)
	    fprintf(fd,'%d %d %d %d %d\n',matquad(i,1),matquad(i,2),matquad(i,3),matquad(i,4),matquad(i,5));
    end
end
if ~isempty(mattri)
    fprintf(fd,'Triangles %d\n',size(mattri,1));
    for i=1:size(mattri,1)
	    fprintf(fd,'%d %d %d %d\n',mattri(i,1),mattri(i,2),mattri(i,3),mattri(i,4));
    end
end
fclose(fd);


%======================%
% computes new nodes   %
%======================%

% DOFs modifications due to nodes renumbering %
tmpdof1 = floor(def.DOF);
tmpdof2 = mod(def.DOF,1);
def.DOF = out(tmpdof1)+tmpdof2;

dof = round(def.DOF.*100);
ind_01 = find(mod(dof-1,100)==0);
ind_01 = ind_01';
find_01 = round(def.DOF(ind_01));
ind_02 = find(mod(dof-2,100)==0);
ind_02 = ind_02';
find_02 = round(def.DOF(ind_02));
ind_03 = find(mod(dof-3,100)==0);
ind_03 = ind_03';
find_03 = round(def.DOF(ind_03));

% deformation along x
model.Node(find_01,5) = model.Node(find_01,5)+scale*def.def(ind_01); 
% deformation along y
model.Node(find_02,6) = model.Node(find_02,6)+scale*def.def(ind_02); 
% deformation along z
model.Node(find_03,7) = model.Node(find_03,7)+scale*def.def(ind_03); 



%===============%
% deformed mesh %
%===============%

%===============%
% nodes writing %
%===============%
[fd,err] = fopen(filename2,'wb');
fprintf(fd,'MeshVersionFormatted 1\n');
fprintf(fd,'Dimension 3\n');

fprintf(fd,'Vertices %d\n',size(model.Node,1));
for i=1:size(model.Node,1)
	fprintf(fd,'%f %f %f %d\n',model.Node(i,5),model.Node(i,6),model.Node(i,7),model.Node(i,4));
end

%==================%
% elements writing %
%==================%
if ~isempty(matquad)
    fprintf(fd,'Quadrilaterals %d\n',size(matquad,1));
    for i=1:size(matquad,1)
	    fprintf(fd,'%d %d %d %d %d\n',matquad(i,1),matquad(i,2),matquad(i,3),matquad(i,4),matquad(i,5));
    end
end
if ~isempty(mattri)
    fprintf(fd,'Triangles %d\n',size(mattri,1));
    for i=1:size(mattri,1)
	    fprintf(fd,'%d %d %d %d\n',mattri(i,1),mattri(i,2),mattri(i,3),mattri(i,4));
    end
end
fclose(fd);



%---------------------------------------------------------------------------------------------------------------------------------------------------%
%						    writemesh_and_constraint									     %
%---------------------------------------------------------------------------------------------------------------------------------------------------%
function out = writemesh_and_constraint(CAM,model,strain)
%==========================================================================%
% displays the mesh defined by model and colors it with the help of strain %
%==========================================================================%
% CAM : filename where information for Medit will be written (without extension)
% model : a structure defining a model. It must contain at least fields .Node and .Elt
% strain : a structure defining the constraints. It must contain at least fields .data and .DOF
% out : vector that gives the corresponding between nodes numerotation in model.Node and in medit
%       out(model.Node(:,1)) = numerotation in medit
% the vector that gives the corresponding between medit and model.Node is not given because it is obvious (model.Node(:,1))

out = [];

[EGroup,nGroup] = getegroup(model.Elt);

matquad = [];
mattri = [];

%========================%
% builds rface and redge %  for the moment, only rface
%========================%
% determination of rface and dimensions
dimrface = 0;
for i=1:nGroup
	ElemF{i} = feutil('getelemf',model.Elt(EGroup(i),:));
	tmpdimrface = size(fe_super('face',ElemF{i}),1);
	dimrface = max(dimrface,tmpdimrface);
end
rface = zeros(size(model.Elt,1),dimrface);
% construction vecteur de correspondance entre les numeros des elements et leurs positions dans model.Elt 
eltid = feutil('eltidfix',model.Elt);
coreltid = zeros(1,max(eltid));
tmpeltid = find(eltid);
coreltid(eltid(tmpeltid)) = tmpeltid;
% construction de rface
if isfield(model,'Stack')
	for i=1:size(model.Stack,1)
		if comstr(model.Stack{i,1},'set') & comstr(model.Stack{i,2},'Face')
			indexrface = (model.Stack{i,3}.data(:,2)-1)*size(rface,1)+coreltid(model.Stack{i,3}.data(:,1))';
			rface(indexrface) = model.Stack{i,3}.ID;
		end
	end
end

%==============================%
% renumerotation of model.Node %
%==============================%
tmpindnode = model.Node(:,1);
model.Node(:,1) = [1:size(model.Node,1)]';
out = zeros(max(tmpindnode),1);
out(tmpindnode) = model.Node(:,1);


%===============%
% nodes writing %
%===============%
filename = strcat([CAM '.mesh']);
[fd,err] = fopen(filename,'wb');
fprintf(fd,'MeshVersionFormatted 1\n');
fprintf(fd,'Dimension 3\n');

fprintf(fd,'Vertices %d\n',size(model.Node,1));
for i=1:size(model.Node,1)
	fprintf(fd,'%f %f %f %d\n',model.Node(i,5),model.Node(i,6),model.Node(i,7),model.Node(i,4));
end

addmatquad = [];
addmattri = [];
%================================%
% building of matquad and mattri %
%================================%
for jGroup=1:nGroup % loop on groups
	elt = model.Elt(EGroup(jGroup):EGroup(jGroup+1)-1,:);
	
	curface = fe_super('sci_face',ElemF{jGroup});
	curnode = fe_super('node',ElemF{jGroup});
	curfaceref = fe_super('face',ElemF{jGroup});
	if size(curfaceref,2)<=4; indic_pb1 = 1;
	else indic_pb1 = 0;
	end
	indquad = [];
	indtri  = [];

	% modifications of elt due to the renumerotation of model.Node
	nbnode = length(curnode);
	elt(2:end,1:nbnode) = reshape(out(elt(2:end,1:nbnode)),size(elt,1)-1,nbnode);

	for i=1:size(curface,1)	
		diffcurface = curface(i,:)-[curface(i,2:end) curface(i,1)];
		finddiffcur = find(diffcurface);
		if length(finddiffcur)==4; indquad(size(indquad,1)+1,1:4) = curface(i,finddiffcur); 
		elseif length(finddiffcur)==3
			indtri(size(indtri,1)+1,1:3) = curface(i,finddiffcur);
		elseif length(finddiffcur)==2
			indtri(size(indtri,1)+1,1:2) = curface(i,finddiffcur);
			indtri(size(indtri,1),3) = indtri(size(indtri,1),2);
		elseif length(finddiffcur)==1
			indtri(size(indtri,1)+1,1) = curface(i,finddiffcur);
			indtri(size(indtri,1),2) = indtri(size(indtri,1),1);
			indtri(size(indtri,1),3) = indtri(size(indtri,1),2);
		end
	end
	sizequadperelt = size(indquad,1);
	sizetriperelt  = size(indtri,1);
	
	for i=2:size(elt,1) % loop on elements
		addmatquad = reshape(elt(i,indquad),sizequadperelt,4);
		addmattri  = reshape(elt(i,indtri),sizetriperelt,3);
		addmatquad = [addmatquad zeros(size(addmatquad,1),1)];
		addmattri = [addmattri zeros(size(addmattri,1),1)];

		
		findrface = find(rface(EGroup(jGroup)+i-1,:));
		valrface = rface(EGroup(jGroup)+i-1,findrface);
		if findrface
			for j=1:length(findrface)
				goodface = elt(i,curfaceref(findrface(j),:)); 
				goodface = unique(goodface); goodface = goodface(length(goodface):-1:1);

				%==============%
				% CASE P1 - Q1 %
				%==============%
				if indic_pb1
					% gestion references des quadrangles
					if length(goodface)==4
					tmpaddmatquad = sort(addmatquad,2); tmpaddmatquad = tmpaddmatquad(:,4:-1:1);
					tmpgood = sum(tmpaddmatquad(:,1:4)==ones(size(addmatquad,1),1)*goodface,2);
					tmpgood2 = find(tmpgood==4);
					addmatquad(tmpgood2,5) = valrface(j);

					% gestion references des triangles
					elseif length(goodface)==3
						tmpaddmattri = sort(addmattri,2); tmpaddmattri = tmpaddmattri(:,3:-1:1);
						tmpgood = sum(tmpaddmattri(:,1:3)==ones(size(addmattri,1),1)*goodface,2);
						tmpgood2 = find(tmpgood==3);
						addmattri(tmpgood2,4) = valrface(j);
					end

				else
				%==============%
				% CASE P2 - Q2 %
				%==============%
				% length goodface > 4 !
					for k=1:size(addmatquad,1)
						if length(intersect(addmatquad(k,1:4),goodface))==4; addmatquad(k,5) = valrface(j);
						end
					end
					for k=1:size(addmattri,1)
						if length(intersect(addmattri(k,1:3),goodface))==3; addmattri(k,4) = valrface(j);
						end
					end

				end % cas pb1 ou pb2
			end
		end
		
		% matquad and mattri
		matquad(size(matquad,1)+1:size(matquad,1)+sizequadperelt,1:5) = addmatquad;
		mattri(size(mattri,1)+1:size(mattri,1)+sizetriperelt,1:4)  = addmattri;
	end % loop on elements
end % loop on groups 

% internal faces are removed
if ~isempty(matquad)
    tmpmatquad = sort(matquad(:,1:4),2);
    [tmpmatquad,indtmpmatquad] = sortrows(tmpmatquad);
    tmpmatquad = tmpmatquad - [tmpmatquad(2:end,:);tmpmatquad(1,:)];
    tmpmatquad = sum(abs(tmpmatquad(:,1:4)),2);
    tmpmatquad(find(~tmpmatquad)+1) = 0;
    indmatquad = indtmpmatquad(find(tmpmatquad));
    matquad = matquad(indmatquad,:);
end

if ~isempty(mattri)   
    tmpmattri = sort(mattri(:,1:3),2);
    [tmpmattri,indtmpmattri] = sortrows(tmpmattri);
    tmpmattri = tmpmattri - [tmpmattri(2:end,:);tmpmattri(1,:)];
    tmpmattri = sum(abs(tmpmattri(:,1:3)),2);
    tmpmattri(find(~tmpmattri)+1) = 0;
    indmattri = indtmpmattri(find(tmpmattri));
    mattri  = mattri(indmattri,:);
end

%==================%
% elements writing %
%==================%
if ~isempty(matquad)
    fprintf(fd,'Quadrilaterals %d\n',size(matquad,1));
    for i=1:size(matquad,1)
	    fprintf(fd,'%d %d %d %d %d\n',matquad(i,1),matquad(i,2),matquad(i,3),matquad(i,4),matquad(i,5));
    end
end
if ~isempty(mattri)
    fprintf(fd,'Triangles %d\n',size(mattri,1));
    for i=1:size(mattri,1)
	    fprintf(fd,'%d %d %d %d\n',mattri(i,1),mattri(i,2),mattri(i,3),mattri(i,4));
    end
end

fclose(fd);

%==========%
% file .bb %
%==========%
filenamebb = strcat([CAM '.bb']);
fd = fopen(filenamebb,'wb');
if isfield(strain,'DOF')
	tmpdof1 = floor(strain.DOF);
	tmpdof2 = (tmpdof1*ones(1,length(tmpdof1)))' - tmpindnode*ones(1,length(tmpdof1));
	[inddata,tmp] = find(~tmpdof2);
	strain.data = strain.data(inddata);
	dim = 3;
	nbmet = 1;
	nbval = length(strain.data);
	typ = 2;
	fprintf(fd,'%d %d %d %d\n',dim,nbmet,nbval,typ);
	for i=1:nbval
	fprintf(fd,'%f ',strain.data(i));
	end
elseif isfield(strain,'EltId')
	warning('medit : this sort of constraint is not treated');
end

fclose(fd);



%---------------------------------------------------------------------------------------------------------------------------------------------------%
%							writemesh_and_defconst									     %
%---------------------------------------------------------------------------------------------------------------------------------------------------%
function [out,scale] = writemesh_and_defconst(filename,model,def,scale,strain)
%================================================================================%
% displays the mesh defined model in a window and the deformation defined by def %
% with colors due to strain in an other window                                   %
%================================================================================%
% filename : filename where information for Medit will be written (without extension)
% model : a structure defining a model. It must contain at least fields .Node and .Elt
% def : a structure defining deformations with at least fields .def and .DOF
% scale : scale 
% strain : a struture defining constraints with at least fields .data and .DOF
% out : vector that gives the corresponding between nodes numerotation in model.Node and in medit
%       out(model.Node(:,1)) = numerotation in medit
% the vector that gives the corresponding between medit and model.Node is not given because it is obvious (model.Node(:,1))

out = [];
filename1 = strcat([filename 'original.mesh']);
filename2 = strcat([filename 'deformed.mesh']);

[EGroup,nGroup] = getegroup(model.Elt);

matquad = [];
mattri = [];

%========================%
% builds rface and redge %  for the moment, only rface
%========================%
% determination of rface and dimensions
dimrface = 0;
for i=1:nGroup
	ElemF{i} = feutil('getelemf',model.Elt(EGroup(i),:));
	tmpdimrface = size(fe_super('face',ElemF{i}),1);
	dimrface = max(dimrface,tmpdimrface);
end
rface = zeros(size(model.Elt,1),dimrface);

eltid = feutil('eltidfix',model.Elt);
coreltid = zeros(1,max(eltid));
tmpeltid = find(eltid);
coreltid(eltid(tmpeltid)) = tmpeltid;
% construction de rface
if isfield(model,'Stack')
	for i=1:size(model.Stack,1)
		if comstr(model.Stack{i,1},'set') & comstr(model.Stack{i,2},'Face')
			indexrface = (model.Stack{i,3}.data(:,2)-1)*size(rface,1)+coreltid(model.Stack{i,3}.data(:,1))';
			rface(indexrface) = model.Stack{i,3}.ID;
		end
	end
end

%==============================%
% renumerotation of model.Node %
%==============================%
tmpindnode = model.Node(:,1);
model.Node(:,1) = [1:size(model.Node,1)]';
out = zeros(max(tmpindnode),1);
out(tmpindnode) = model.Node(:,1);

addmatquad = [];
addmattri = [];
%================================%
% building of matquad and mattri %
%================================%
for jGroup=1:nGroup % loop on groups
	elt = model.Elt(EGroup(jGroup):EGroup(jGroup+1)-1,:);
	
	curface = fe_super('sci_face',ElemF{jGroup});
	curnode = fe_super('node',ElemF{jGroup});
	curfaceref = fe_super('face',ElemF{jGroup});
	if size(curfaceref,2)<=4; indic_pb1 = 1;
	else indic_pb1 = 0;
	end
	indquad = [];
	indtri  = [];

	% modifications of elt due to the renumerotation of model.Node
	nbnode = length(curnode);
	elt(2:end,1:nbnode) = reshape(out(elt(2:end,1:nbnode)),size(elt,1)-1,nbnode);

	for i=1:size(curface,1)	
		diffcurface = curface(i,:)-[curface(i,2:end) curface(i,1)];
		finddiffcur = find(diffcurface);
		if length(finddiffcur)==4; indquad(size(indquad,1)+1,1:4) = curface(i,finddiffcur); 
		elseif length(finddiffcur)==3
			indtri(size(indtri,1)+1,1:3) = curface(i,finddiffcur);
		elseif length(finddiffcur)==2
			indtri(size(indtri,1)+1,1:2) = curface(i,finddiffcur);
			indtri(size(indtri,1),3) = indtri(size(indtri,1),2);
		elseif length(finddiffcur)==1
			indtri(size(indtri,1)+1,1) = curface(i,finddiffcur);
			indtri(size(indtri,1),2) = indtri(size(indtri,1),1);
			indtri(size(indtri,1),3) = indtri(size(indtri,1),2);
		end
	end
	sizequadperelt = size(indquad,1);
	sizetriperelt  = size(indtri,1);
	
	for i=2:size(elt,1) % loop on elements
		addmatquad = reshape(elt(i,indquad),sizequadperelt,4);
		addmattri  = reshape(elt(i,indtri),sizetriperelt,3);
		addmatquad = [addmatquad zeros(size(addmatquad,1),1)];
		addmattri = [addmattri zeros(size(addmattri,1),1)];

		
		findrface = find(rface(EGroup(jGroup)+i-1,:));
		valrface = rface(EGroup(jGroup)+i-1,findrface);
		if findrface
			for j=1:length(findrface)
				goodface = elt(i,curfaceref(findrface(j),:)); 
				goodface = unique(goodface); goodface = goodface(length(goodface):-1:1);

				%==============%
				% CASE P1 - Q1 %
				%==============%
				if indic_pb1
					% quads
					if length(goodface)==4
					tmpaddmatquad = sort(addmatquad,2); tmpaddmatquad = tmpaddmatquad(:,4:-1:1);
					tmpgood = sum(tmpaddmatquad(:,1:4)==ones(size(addmatquad,1),1)*goodface,2);
					tmpgood2 = find(tmpgood==4);
					addmatquad(tmpgood2,5) = valrface(j);

					% triangles
					elseif length(goodface)==3
						tmpaddmattri = sort(addmattri,2); tmpaddmattri = tmpaddmattri(:,3:-1:1);
						tmpgood = sum(tmpaddmattri(:,1:3)==ones(size(addmattri,1),1)*goodface,2);
						tmpgood2 = find(tmpgood==3);
						addmattri(tmpgood2,4) = valrface(j);
					end

				else
				%==============%
				% CASE P2 - Q2 %
				%==============%
				% length goodface > 4 !
					for k=1:size(addmatquad,1)
						if length(intersect(addmatquad(k,1:4),goodface))==4; addmatquad(k,5) = valrface(j);
						end
					end
					for k=1:size(addmattri,1)
						if length(intersect(addmattri(k,1:3),goodface))==3; addmattri(k,4) = valrface(j);
						end
					end

				end % cas pb1 ou pb2
			end
		end
		
		% matquad and mattri
		matquad(size(matquad,1)+1:size(matquad,1)+sizequadperelt,1:5) = addmatquad;
		mattri(size(mattri,1)+1:size(mattri,1)+sizetriperelt,1:4)  = addmattri;
	end % loop on elements
end % loop on groups 

% internal faces are removed
if ~isempty(matquad)
    tmpmatquad = sort(matquad(:,1:4),2);
    [tmpmatquad,indtmpmatquad] = sortrows(tmpmatquad);
    tmpmatquad = tmpmatquad - [tmpmatquad(2:end,:);tmpmatquad(1,:)];
    tmpmatquad = sum(abs(tmpmatquad(:,1:4)),2);
    tmpmatquad(find(~tmpmatquad)+1) = 0;
    indmatquad = indtmpmatquad(find(tmpmatquad));
    matquad = matquad(indmatquad,:);
end

if ~isempty(mattri)
    tmpmattri = sort(mattri(:,1:3),2);
    [tmpmattri,indtmpmattri] = sortrows(tmpmattri);
    tmpmattri = tmpmattri - [tmpmattri(2:end,:);tmpmattri(1,:)];
    tmpmattri = sum(abs(tmpmattri(:,1:3)),2);
    tmpmattri(find(~tmpmattri)+1) = 0;
    indmattri = indtmpmattri(find(tmpmattri));
    mattri  = mattri(indmattri,:);
end

%===============%
% original mesh %
%===============%

%===============%
% nodes writing %
%===============%
[fd,err] = fopen(filename1,'wb');
fprintf(fd,'MeshVersionFormatted 1\n');
fprintf(fd,'Dimension 3\n');

fprintf(fd,'Vertices %d\n',size(model.Node,1));
for i=1:size(model.Node,1)
	fprintf(fd,'%f %f %f %d\n',model.Node(i,5),model.Node(i,6),model.Node(i,7),model.Node(i,4));
end

%==================%
% elements writing %
%==================%
if ~isempty(matquad)
    fprintf(fd,'Quadrilaterals %d\n',size(matquad,1));
    for i=1:size(matquad,1)
	    fprintf(fd,'%d %d %d %d %d\n',matquad(i,1),matquad(i,2),matquad(i,3),matquad(i,4),matquad(i,5));
    end
end
if ~isempty(mattri)
    fprintf(fd,'Triangles %d\n',size(mattri,1));
    for i=1:size(mattri,1)
	    fprintf(fd,'%d %d %d %d\n',mattri(i,1),mattri(i,2),mattri(i,3),mattri(i,4));
    end
end
fclose(fd);


%======================%
% computes new nodes   %
%======================%

% DOFs modifications due to nodes renumbering %
tmpdof1 = floor(def.DOF);
tmpdof2 = mod(def.DOF,1);
def.DOF = out(tmpdof1)+tmpdof2;

dof = round(def.DOF.*100);
ind_01 = find(mod(dof-1,100)==0);
ind_01 = ind_01';
find_01 = round(def.DOF(ind_01));
ind_02 = find(mod(dof-2,100)==0);
ind_02 = ind_02';
find_02 = round(def.DOF(ind_02));
ind_03 = find(mod(dof-3,100)==0);
ind_03 = ind_03';
find_03 = round(def.DOF(ind_03));

% deformation along x
model.Node(find_01,5) = model.Node(find_01,5)+scale*def.def(ind_01); 
% deformation along y
model.Node(find_02,6) = model.Node(find_02,6)+scale*def.def(ind_02); 
% deformation along z
model.Node(find_03,7) = model.Node(find_03,7)+scale*def.def(ind_03); 



%===============%
% deformed mesh %
%===============%

%===============%
% nodes writing %
%===============%
[fd,err] = fopen(filename2,'wb');
fprintf(fd,'MeshVersionFormatted 1\n');
fprintf(fd,'Dimension 3\n');

fprintf(fd,'Vertices %d\n',size(model.Node,1));
for i=1:size(model.Node,1)
	fprintf(fd,'%f %f %f %d\n',model.Node(i,5),model.Node(i,6),model.Node(i,7),model.Node(i,4));
end

%==================%
% elements writing %
%==================%
if ~isempty(matquad)
    fprintf(fd,'Quadrilaterals %d\n',size(matquad,1));
    for i=1:size(matquad,1)
	    fprintf(fd,'%d %d %d %d %d\n',matquad(i,1),matquad(i,2),matquad(i,3),matquad(i,4),matquad(i,5));
    end
end
if ~isempty(mattri)
    fprintf(fd,'Triangles %d\n',size(mattri,1));
    for i=1:size(mattri,1)
	    fprintf(fd,'%d %d %d %d\n',mattri(i,1),mattri(i,2),mattri(i,3),mattri(i,4));
    end
end
fclose(fd);


%==========%
% file .bb %
%==========%
filenamebb = strcat([filename 'deformed.bb']);
fd = fopen(filenamebb,'wb');
if isfield(strain,'DOF')
	tmpdof1 = floor(strain.DOF);
	tmpdof2 = (tmpdof1*ones(1,length(tmpdof1)))' - tmpindnode*ones(1,length(tmpdof1));
	[inddata,tmp] = find(~tmpdof2);
	strain.data = strain.data(inddata);
	dim = 3;
	nbmet = 1;
	nbval = length(strain.data);
	typ = 2;
	fprintf(fd,'%d %d %d %d\n',dim,nbmet,nbval,typ);
	for i=1:nbval
	fprintf(fd,'%f ',strain.data(i));
	end
elseif isfield(strain,'EltId')
	warning('medit : this sort of constraint is not treated');
end

fclose(fd);



%---------------------------------------------------------------------------------------------------------------------------------------------------%
%							     animate_deformed									     %
%---------------------------------------------------------------------------------------------------------------------------------------------------%
function [out,scale] = animate_deformed(filename,model,def,scale,nbimag)
%===================================================%
% animates the deformations defined by def on model %
%===================================================%
% filename : filename where information for Medit will be written (without extension)
% model : a structure defining a model. It must contain at least fields .Node and .Elt
% def : a structure defining deformations with at least fields .def and .DOF
% scale : scale 
% nbimag : number of files to create for animation
% out : vector that gives the corresponding between nodes numerotation in model.Node and in medit
%       out(model.Node(:,1)) = numerotation in medit
% the vector that gives the corresponding between medit and model.Node is not given because it is obvious (model.Node(:,1))

out = [];

[EGroup,nGroup] = getegroup(model.Elt);

matquad = [];
mattri = [];

%========================%
% builds rface and redge %  for the moment, only rface
%========================%
% determination of rface and dimensions
dimrface = 0;
for i=1:nGroup
	ElemF{i} = feutil('getelemf',model.Elt(EGroup(i),:));
	tmpdimrface = size(fe_super('face',ElemF{i}),1);
	dimrface = max(dimrface,tmpdimrface);
end
rface = zeros(size(model.Elt,1),dimrface);

eltid = feutil('eltidfix',model.Elt);
coreltid = zeros(1,max(eltid));
tmpeltid = find(eltid);
coreltid(eltid(tmpeltid)) = tmpeltid;
% construction de rface
if isfield(model,'Stack')
	for i=1:size(model.Stack,1)
		if comstr(model.Stack{i,1},'set') & comstr(model.Stack{i,2},'Face')
			indexrface = (model.Stack{i,3}.data(:,2)-1)*size(rface,1)+coreltid(model.Stack{i,3}.data(:,1))';
			rface(indexrface) = model.Stack{i,3}.ID;
		end
	end
end

%==============================%
% renumerotation of model.Node %
%==============================%
tmpindnode = model.Node(:,1);
model.Node(:,1) = [1:size(model.Node,1)]';
out = zeros(max(tmpindnode),1);
out(tmpindnode) = model.Node(:,1);

addmatquad = [];
addmattri = [];
%================================%
% building of matquad and mattri %
%================================%
for jGroup=1:nGroup % loop on groups
	elt = model.Elt(EGroup(jGroup):EGroup(jGroup+1)-1,:);
	
	curface = fe_super('sci_face',ElemF{jGroup});
	curnode = fe_super('node',ElemF{jGroup});
	curfaceref = fe_super('face',ElemF{jGroup});
	if size(curfaceref,2)<=4; indic_pb1 = 1;
	else indic_pb1 = 0;
	end
	indquad = [];
	indtri  = [];

	% modifications of elt due to the renumerotation of model.Node
	nbnode = length(curnode);
	elt(2:end,1:nbnode) = reshape(out(elt(2:end,1:nbnode)),size(elt,1)-1,nbnode);

	for i=1:size(curface,1)	
		diffcurface = curface(i,:)-[curface(i,2:end) curface(i,1)];
		finddiffcur = find(diffcurface);
		if length(finddiffcur)==4; indquad(size(indquad,1)+1,1:4) = curface(i,finddiffcur); 
		elseif length(finddiffcur)==3
			indtri(size(indtri,1)+1,1:3) = curface(i,finddiffcur);
		elseif length(finddiffcur)==2
			indtri(size(indtri,1)+1,1:2) = curface(i,finddiffcur);
			indtri(size(indtri,1),3) = indtri(size(indtri,1),2);
		elseif length(finddiffcur)==1
			indtri(size(indtri,1)+1,1) = curface(i,finddiffcur);
			indtri(size(indtri,1),2) = indtri(size(indtri,1),1);
			indtri(size(indtri,1),3) = indtri(size(indtri,1),2);
		end
	end
	sizequadperelt = size(indquad,1);
	sizetriperelt  = size(indtri,1);
	
	for i=2:size(elt,1) % loop on elements
		addmatquad = reshape(elt(i,indquad),sizequadperelt,4);
		addmattri  = reshape(elt(i,indtri),sizetriperelt,3);
		addmatquad = [addmatquad zeros(size(addmatquad,1),1)];
		addmattri = [addmattri zeros(size(addmattri,1),1)];

		
		findrface = find(rface(EGroup(jGroup)+i-1,:));
		valrface = rface(EGroup(jGroup)+i-1,findrface);
		if findrface
			for j=1:length(findrface)
				goodface = elt(i,curfaceref(findrface(j),:)); 
				goodface = unique(goodface); goodface = goodface(length(goodface):-1:1);

				%==============%
				% CASE P1 - Q1 %
				%==============%
				if indic_pb1
					% quads
					if length(goodface)==4
					tmpaddmatquad = sort(addmatquad,2); tmpaddmatquad = tmpaddmatquad(:,4:-1:1);
					tmpgood = sum(tmpaddmatquad(:,1:4)==ones(size(addmatquad,1),1)*goodface,2);
					tmpgood2 = find(tmpgood==4);
					addmatquad(tmpgood2,5) = valrface(j);

					% triangles
					elseif length(goodface)==3
						tmpaddmattri = sort(addmattri,2); tmpaddmattri = tmpaddmattri(:,3:-1:1);
						tmpgood = sum(tmpaddmattri(:,1:3)==ones(size(addmattri,1),1)*goodface,2);
						tmpgood2 = find(tmpgood==3);
						addmattri(tmpgood2,4) = valrface(j);
					end

				else
				%==============%
				% CASE P2 - Q2 %
				%==============%
				% length goodface > 4 !
					for k=1:size(addmatquad,1)
						if length(intersect(addmatquad(k,1:4),goodface))==4; addmatquad(k,5) = valrface(j);
						end
					end
					for k=1:size(addmattri,1)
						if length(intersect(addmattri(k,1:3),goodface))==3; addmattri(k,4) = valrface(j);
						end
					end

				end % case p1-q1
			end
		end
		
		% matquad and mattri
		matquad(size(matquad,1)+1:size(matquad,1)+sizequadperelt,1:5) = addmatquad;
		mattri(size(mattri,1)+1:size(mattri,1)+sizetriperelt,1:4)  = addmattri;
	end % loop on elements
end % loop on groups 

% internal faces are removed
if ~isempty(matquad)
    tmpmatquad = sort(matquad(:,1:4),2);
    [tmpmatquad,indtmpmatquad] = sortrows(tmpmatquad);
    tmpmatquad = tmpmatquad - [tmpmatquad(2:end,:);tmpmatquad(1,:)];
    tmpmatquad = sum(abs(tmpmatquad(:,1:4)),2);
    tmpmatquad(find(~tmpmatquad)+1) = 0;
    indmatquad = indtmpmatquad(find(tmpmatquad));
    matquad = matquad(indmatquad,:);
end
    
if ~isempty(mattri)
    tmpmattri = sort(mattri(:,1:3),2);
    [tmpmattri,indtmpmattri] = sortrows(tmpmattri);
    tmpmattri = tmpmattri - [tmpmattri(2:end,:);tmpmattri(1,:)];
    tmpmattri = sum(abs(tmpmattri(:,1:3)),2);
    tmpmattri(find(~tmpmattri)+1) = 0;
    indmattri = indtmpmattri(find(tmpmattri));
    mattri  = mattri(indmattri,:);
end

% DOFs modifications due to nodes renumbering %
tmpdof1 = floor(def.DOF);
tmpdof2 = mod(def.DOF,1);
def.DOF = out(tmpdof1)+tmpdof2;

dof = round(def.DOF.*100);
ind_01 = find(mod(dof-1,100)==0);
ind_01 = ind_01';
find_01 = round(def.DOF(ind_01));
ind_02 = find(mod(dof-2,100)==0);
ind_02 = ind_02';
find_02 = round(def.DOF(ind_02));
ind_03 = find(mod(dof-3,100)==0);
ind_03 = ind_03';
find_03 = round(def.DOF(ind_03));

%======================%
% write the mesh files %
%======================%
for i = 1:nbimag

t = 4*pi*i/nbimag;
sint = sin(t);

%======================%
% computes new nodes   %
%======================%
modelnew.Node = model.Node;

% deformation along x
modelnew.Node(find_01,5) = model.Node(find_01,5)+scale*def.def(ind_01)*sint; 
% deformation along y
modelnew.Node(find_02,6) = model.Node(find_02,6)+scale*def.def(ind_02)*sint; 
% deformation along z
modelnew.Node(find_03,7) = model.Node(find_03,7)+scale*def.def(ind_03)*sint; 


%===============%
% nodes writing %
%===============%
[fd,err] = fopen(NumName(filename,i,'.mesh'),'wb');
fprintf(fd,'MeshVersionFormatted 1\n');
fprintf(fd,'Dimension 3\n');

fprintf(fd,'Vertices %d\n',size(modelnew.Node,1));
for i=1:size(modelnew.Node,1)
	fprintf(fd,'%f %f %f %d\n',modelnew.Node(i,5),modelnew.Node(i,6),modelnew.Node(i,7),modelnew.Node(i,4));
end

%==================%
% elements writing %
%==================%
if ~isempty(matquad)
    fprintf(fd,'Quadrilaterals %d\n',size(matquad,1));
    for i=1:size(matquad,1)
	    fprintf(fd,'%d %d %d %d %d\n',matquad(i,1),matquad(i,2),matquad(i,3),matquad(i,4),matquad(i,5));
    end
end
if ~isempty(mattri)
    fprintf(fd,'Triangles %d\n',size(mattri,1));
    for i=1:size(mattri,1)
	    fprintf(fd,'%d %d %d %d\n',mattri(i,1),mattri(i,2),mattri(i,3),mattri(i,4));
    end
end
fclose(fd);

end % for i = 1:nbimag




%---------------------------------------------------------------------------------------------------------------------------------------------------%
%						     animate_deformed_constraint								     %
%---------------------------------------------------------------------------------------------------------------------------------------------------%

function [out,scale] = animate_deformed_constraint(filename,model,def,scale,nbimag,strain)
%===========================================================================================%
% animates the deformations defined by def on model and colors them with the help of strain %
%===========================================================================================%
% filename : filename where information for Medit will be written (without extension)
% model : a structure defining a model. It must contain at least fields .Node and .Elt
% def : a structure defining deformations with at least fields .def and .DOF
% scale : scale 
% nb_imag : number of files to create for animation
% strain : structure with at least fields .data and .DOF or .EltId
% out : vector that gives the corresponding between nodes numerotation in model.Node and in medit
%       out(model.Node(:,1)) = numerotation in medit
% the vector that gives the corresponding between medit and model.Node is not given because it is obvious (model.Node(:,1))

out = [];

[EGroup,nGroup] = getegroup(model.Elt);

matquad = [];
mattri = [];

%========================%
% builds rface and redge %  for the moment, only rface
%========================%
% determination of rface and dimensions
dimrface = 0;
for i=1:nGroup
	ElemF{i} = feutil('getelemf',model.Elt(EGroup(i),:));
	tmpdimrface = size(fe_super('face',ElemF{i}),1);
	dimrface = max(dimrface,tmpdimrface);
end
rface = zeros(size(model.Elt,1),dimrface);

eltid = feutil('eltidfix',model.Elt);
coreltid = zeros(1,max(eltid));
tmpeltid = find(eltid);
coreltid(eltid(tmpeltid)) = tmpeltid;

if isfield(model,'Stack')
	for i=1:size(model.Stack,1)
		if comstr(model.Stack{i,1},'set') & comstr(model.Stack{i,2},'Face')
			indexrface = (model.Stack{i,3}.data(:,2)-1)*size(rface,1)+coreltid(model.Stack{i,3}.data(:,1))';
			rface(indexrface) = model.Stack{i,3}.ID;
		end
	end
end

%==============================%
% renumerotation of model.Node %
%==============================%
tmpindnode = model.Node(:,1);
model.Node(:,1) = [1:size(model.Node,1)]';
out = zeros(max(tmpindnode),1);
out(tmpindnode) = model.Node(:,1);

addmatquad = [];
addmattri = [];
%================================%
% building of matquad and mattri %
%================================%
for jGroup=1:nGroup % loop on groups
	elt = model.Elt(EGroup(jGroup):EGroup(jGroup+1)-1,:);
	
	curface = fe_super('sci_face',ElemF{jGroup});
	curnode = fe_super('node',ElemF{jGroup});
	curfaceref = fe_super('face',ElemF{jGroup});
	if size(curfaceref,2)<=4; indic_pb1 = 1;
	else indic_pb1 = 0;
	end
	indquad = [];
	indtri  = [];

	% modifications of elt due to the renumerotation of model.Node
	nbnode = length(curnode);
	elt(2:end,1:nbnode) = reshape(out(elt(2:end,1:nbnode)),size(elt,1)-1,nbnode);

	for i=1:size(curface,1)	
		diffcurface = curface(i,:)-[curface(i,2:end) curface(i,1)];
		finddiffcur = find(diffcurface);
		if length(finddiffcur)==4; indquad(size(indquad,1)+1,1:4) = curface(i,finddiffcur); 
		elseif length(finddiffcur)==3
			indtri(size(indtri,1)+1,1:3) = curface(i,finddiffcur);
		elseif length(finddiffcur)==2
			indtri(size(indtri,1)+1,1:2) = curface(i,finddiffcur);
			indtri(size(indtri,1),3) = indtri(size(indtri,1),2);
		elseif length(finddiffcur)==1
			indtri(size(indtri,1)+1,1) = curface(i,finddiffcur);
			indtri(size(indtri,1),2) = indtri(size(indtri,1),1);
			indtri(size(indtri,1),3) = indtri(size(indtri,1),2);
		end
	end
	sizequadperelt = size(indquad,1);
	sizetriperelt  = size(indtri,1);
	
	for i=2:size(elt,1) % loop on elements
		addmatquad = reshape(elt(i,indquad),sizequadperelt,4);
		addmattri  = reshape(elt(i,indtri),sizetriperelt,3);
		addmatquad = [addmatquad zeros(size(addmatquad,1),1)];
		addmattri = [addmattri zeros(size(addmattri,1),1)];

		
		findrface = find(rface(EGroup(jGroup)+i-1,:));
		valrface = rface(EGroup(jGroup)+i-1,findrface);
		if findrface
			for j=1:length(findrface)
				goodface = elt(i,curfaceref(findrface(j),:)); 
				goodface = unique(goodface); goodface = goodface(length(goodface):-1:1);

				%==============%
				% CASE P1 - Q1 %
				%==============%
				if indic_pb1
					% quads
					if length(goodface)==4
					tmpaddmatquad = sort(addmatquad,2); tmpaddmatquad = tmpaddmatquad(:,4:-1:1);
					tmpgood = sum(tmpaddmatquad(:,1:4)==ones(size(addmatquad,1),1)*goodface,2);
					tmpgood2 = find(tmpgood==4);
					addmatquad(tmpgood2,5) = valrface(j);

					% triangles
					elseif length(goodface)==3
						tmpaddmattri = sort(addmattri,2); tmpaddmattri = tmpaddmattri(:,3:-1:1);
						tmpgood = sum(tmpaddmattri(:,1:3)==ones(size(addmattri,1),1)*goodface,2);
						tmpgood2 = find(tmpgood==3);
						addmattri(tmpgood2,4) = valrface(j);
					end

				else
				%==============%
				% CASE P2 - Q2 %
				%==============%
				% length goodface > 4 !
					for k=1:size(addmatquad,1)
						if length(intersect(addmatquad(k,1:4),goodface))==4; addmatquad(k,5) = valrface(j);
						end
					end
					for k=1:size(addmattri,1)
						if length(intersect(addmattri(k,1:3),goodface))==3; addmattri(k,4) = valrface(j);
						end
					end

				end % cas pb1 ou pb2
			end
		end
		
		% matquad and mattri
		matquad(size(matquad,1)+1:size(matquad,1)+sizequadperelt,1:5) = addmatquad;
		mattri(size(mattri,1)+1:size(mattri,1)+sizetriperelt,1:4)  = addmattri;
	end % loop on elements
end % loop on groups 

% internal faces are removed
if ~isempty(matquad)
    tmpmatquad = sort(matquad(:,1:4),2);
    [tmpmatquad,indtmpmatquad] = sortrows(tmpmatquad);
    tmpmatquad = tmpmatquad - [tmpmatquad(2:end,:);tmpmatquad(1,:)];
    tmpmatquad = sum(abs(tmpmatquad(:,1:4)),2);
    tmpmatquad(find(~tmpmatquad)+1) = 0;
    indmatquad = indtmpmatquad(find(tmpmatquad));
    matquad = matquad(indmatquad,:);
end

if ~isempty(mattri)
    tmpmattri = sort(mattri(:,1:3),2);
    [tmpmattri,indtmpmattri] = sortrows(tmpmattri);
    tmpmattri = tmpmattri - [tmpmattri(2:end,:);tmpmattri(1,:)];
    tmpmattri = sum(abs(tmpmattri(:,1:3)),2);
    tmpmattri(find(~tmpmattri)+1) = 0;
    indmattri = indtmpmattri(find(tmpmattri));
    mattri  = mattri(indmattri,:);
end

% DOFs modifications due to nodes renumbering %
tmpdof1 = floor(def.DOF);
tmpdof2 = mod(def.DOF,1);
def.DOF = out(tmpdof1)+tmpdof2;

dof = round(def.DOF.*100);
ind_01 = find(mod(dof-1,100)==0);
ind_01 = ind_01';
find_01 = round(def.DOF(ind_01));
ind_02 = find(mod(dof-2,100)==0);
ind_02 = ind_02';
find_02 = round(def.DOF(ind_02));
ind_03 = find(mod(dof-3,100)==0);
ind_03 = ind_03';
find_03 = round(def.DOF(ind_03));

%======================%
% write the mesh files %
%======================%
for j1 = 1:nbimag

t = 4*pi*j1/nbimag;
sint = sin(t);

%======================%
% computes new nodes   %
%======================%
modelnew.Node = model.Node;

% deformation along x
modelnew.Node(find_01,5) = model.Node(find_01,5)+scale*def.def(ind_01)*sint; 
% deformation along y
modelnew.Node(find_02,6) = model.Node(find_02,6)+scale*def.def(ind_02)*sint; 
% deformation along z
modelnew.Node(find_03,7) = model.Node(find_03,7)+scale*def.def(ind_03)*sint; 


%===============%
% nodes writing %
%===============%
[fd,err] = fopen(NumName(filename,j1,'.mesh'),'wb');
fprintf(fd,'MeshVersionFormatted 1\n');
fprintf(fd,'Dimension 3\n');

fprintf(fd,'Vertices %d\n',size(modelnew.Node,1));
fprintf(fd,'%f %f %f %d\n',modelnew.Node(:,[5 6 7 4])');


%==================%
% elements writing %
%==================%
if ~isempty(matquad)
    fprintf(fd,'Quadrilaterals %d\n',size(matquad,1));
	fprintf(fd,'%d %d %d %d %d\n',matquad(:,1:5)');
end
if ~isempty(mattri)
    fprintf(fd,'Triangles %d\n',size(mattri,1));
	fprintf(fd,'%d %d %d %d\n',mattri(:,1:4)');
end
fclose(fd);

%===============================%
% write the associated .bb file %
%===============================%
[fd,err] = fopen(NumName(filename,j1,'.bb'),'wb');
if isfield(strain,'DOF')
	tmpdof1 = floor(strain.DOF);
	tmpdof2 = (tmpdof1*ones(1,length(tmpdof1)))' - tmpindnode*ones(1,length(tmpdof1));
	[inddata,tmp] = find(~tmpdof2);
	strain.data = strain.data(inddata);
	dim = 3;
	nbmet = 1;
	nbval = length(strain.data);
	typ = 2;
	fprintf(fd,'%d %d %d %d\n',dim,nbmet,nbval,typ);
	fprintf(fd,'%f ',strain.data);
elseif isfield(strain,'EltId')
	warning('medit : this sort of constraint is not treated');
end

end % for i = 1:nbimag

%% Generate file name with index number
function filenamei=NumName(filename,i,ext)

if 1==2 % Old version
 if i<10; filenamei =  strcat([filename '.00' num2str(i) ext]);
 elseif i<100; filenamei =  strcat([filename '.0' num2str(i) ext]);
 else;filenamei =  strcat([filename '.' num2str(i) ext]); 
 end
else % medit 3 compatible version
    filenamei =  strcat([filename '.' num2str(i) ext]);
end

