function ks=mumps_utils(varargin);


% Gateway function for inclusion of the MUMPS solver into the
% OpenFEM ofact object
%
% Call syntax:
%	kinv=ofact('fact',k);	performs standard fact
%	kinv=ofact('zfact',k);	performs complex fact
%	kinv=ofact('symbfact',k);	symbolic fact only
%	kinv=ofact('zsymbfact',k);	complex symbolic
%	kinv=ofact('numfact',k,kinv);	numerical fact all cases
%	x=kinv\f;			solution all cases
%	ofact('clear',kinv);	IMPORTANT: clear ofact object
%	x=ofact(k,f); does fact / solve / clear
%



Cam=varargin{1};carg=2; param=[];
% -------------------------------------------------------------------------
if comstr(Cam,'method') % This is where known methods are defined

	ks.name='mumps';
	ks.header='MUMPS sparse solver';
	ks.SymRenumber='';
	ks.FactBuild='ks=mumps_utils(''fact'',k,ks);';
	ks.Solve='q=mumps_utils(''solve'',k,full(b));';
	ks.Clear='mumps_utils(''clear'',ks);';
	ks.Available=exist('dmumpsmex','file') == 3;
	ks.HandlesComplex=exist('zmumpsmex','file') == 3;
	ks.param=[];



	% This is where the factorization is actually performed
elseif comstr(Cam,'fact')||comstr(Cam,'zfact')...
		||comstr(Cam,'symbfact')||comstr(Cam,'zsymbfact')...
		||comstr(Cam,'numfact')

	k=varargin{carg};carg=carg+1;
	ks=varargin{carg};carg=carg+1;

	if carg<nargin
		opt=varargin{carg};carg=carg+1;
	end
	
	%
	% Init when needed
	%
	if ~comstr(Cam,'numfact')
		id = initmumps;
	end
	% Flag for symmetry should be set here when needed
	%id.SYM = 1;
	if comstr(Cam,'zfact')||comstr(Cam,'zsymbfact')
		id = zmumps(id);
	elseif comstr(Cam,'fact')||comstr(Cam,'symbfact')
		id = dmumps(id);
	end
	%
	% Define job for each case
	%
	if comstr(Cam,'fact')||comstr(Cam,'zfact')
		id.JOB = 4;
	elseif comstr(Cam,'symbfact')||comstr(Cam,'zsymbfact')
		id.JOB = 1;
	else % numfact
		id = ks.data.mumpsid;
		id.JOB = 2;
	end
	%
	% Perform actual job
	%
	if comstr(Cam,'zfact')||comstr(Cam,'zsymbfact')...
			||(comstr(Cam,'numfact') && id.TYPE ==2)
		id = zmumps(id,k);
	else
		id = dmumps(id,k);
	end
	%
	% Store what's needed in ks
	%
	if ~comstr(Cam,'numfact')
		ks.ty = [10, 0, size( k, 1 )];
		% Warning: for some reason k needs to be provided at the solve stage -> copied in data.k
		ks.data = struct('mumpsid',id,'k',k);
	else
		ks.data.mumpsid = id;
	end


elseif comstr(Cam,'solve')

	k=varargin{2}; b=varargin{3};
	if issparse(b); b=full(b);end
	k.data.mumpsid.JOB=3;
	k.data.mumpsid.RHS=b;
	if k.data.mumpsid.TYPE == 2 % complex matrix
		k.data.mumpsid=zmumps(k.data.mumpsid,k.data.k);
	else
		k.data.mumpsid=dmumps(k.data.mumpsid,k.data.k);
	end
	q=k.data.mumpsid.SOL;
	ks=q;

elseif comstr(Cam,'clear')  % - - - - - - - - - - - - - - - - - - - - - -

	k=varargin{2};
	k.data.mumpsid.JOB=-2;
	if k.data.mumpsid.TYPE == 2 % complex matrix
		k.data.mumpsid = zmumps(k.data.mumpsid);
	else
		k.data.mumpsid = dmumps(k.data.mumpsid);
	end


elseif comstr(Cam,'silent')  % - - - - - - - - - - - - - - - - - - - - - -

	warning('Nothing done');

end % - - - - - - - - - - - - - - - - - - - - - -

