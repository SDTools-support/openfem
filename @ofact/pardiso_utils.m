function ks=pardiso_utils(varargin);


% Gateway function for inclusion of the PARDISO solver into the 
% OpenFEM ofact object
%
% http://www.computational.unibas.ch/cs/scicomp/software/pardiso/

%       $Revision: 1.15 $  $Date: 2009/10/30 08:10:53 $




Cam=varargin{1};carg=2; param=[];
% -------------------------------------------------------------------------
if comstr(Cam,'method') % This is where known methods are defined

   ks.name='pardiso';
   ks.header='PARDISO sparse solver';
   ks.SymRenumber='';
                                               %msg CPU
   ks.FactBuild='ks=pardiso_utils(''fact'',k,ks,[0 1]);';
   ks.Solve='q=pardiso_utils(''solve'',k,full(b));';
   ks.Clear='pardiso_utils(''clear'',ks);';
   ks.Available=exist('pardiso','file')==3;
   ks.HandlesComplex=0;

   if sp_util('mwIndex')==8; fun=@int32; else; fun=@int32;end
   param=feval(fun,[1 2 2 0 0 0 16 2 0 13 1 ...
    0 0 0 0 0 0 ...
    -1 -1 0 ]);
   param(18:20)=0; % no MFlop computation 
   param(101:105)=[1    1    1    0  0]; %maxfct mnum msglvl error nonsym
   %maxfct = 1;     /* Maximum number of numerical factorizations. */
   %mnum = 1;       /* Which factorization to use. */
   %msglvl = 0;     /* Print statistical information in file */
   %error = 0;      /* Initialize error flag */
   ks.param=param;


 try;
  i1=comstr(getenv('OMP_NUM_THREADS'),[-1 1]);
  ks.FactBuild=sprintf('ks=pardiso_utils(''fact'',k,ks,[0 %i]);',i1);
 catch
  if ~isunix
   ks.FactBuild=sprintf('ks=pardiso_utils(''fact'',k,ks,[0 %i]);', ...
    getenv('NUMBER_OF_PROCESSORS'));
  end
 end

% This is where the factorization is actually - - - - - - - - - - - - -
elseif comstr(Cam,'numfact')||comstr(Cam,'symbfact')||comstr(Cam,'fact')

 k=varargin{carg};carg=carg+1;
 ks=varargin{carg};carg=carg+1;

 if carg<nargin; opt=varargin{carg};carg=carg+1;
 else; opt=[0 1]; % msglvl CPU
 end
 param=ks.method.param; param(103)=opt(1); param(3)=opt(2); % CPU
 if ~isempty(regexp(Cam,'nonsym','once'));
	 param(105)=1;
	 % IMPORTANT: Pardiso storage along lines, hence transpose nonsym matrices
	 k=transpose(k);
 end
 fun=@int32; if sp_util('mwIndex')==8; k=v_handle('mkls',k);end
 ks.data=feval(fun,param);
 if comstr(Cam,'fact') 
 % fact performs BOTH symbolic AND numeric factorizationS
	i1=pardiso('symbfact',k,ks.data);
	i1=pardiso('numfact',k,ks.data);
 elseif comstr(Cam,'symbfact')
 % symbolic only
	i1=pardiso('symbfact',k,ks.data);
 elseif comstr(Cam,'numfact')
 % numeric only
	i1=pardiso('numfact',k,ks.data);
 end
 ks.ty(1:3)=[5 i1 size(k,1)];ks.ind='pardiso';

%iparm[0] = 1;   /* No solver default */
%iparm[1] = 2;   /* Fill-in reordering from METIS */
%iparm[2] = omp_get_max_threads(); /* Numbers of processors, value of OMP_NUM_THREADS */
%iparm[3] = 0;   /* No iterative-direct algorithm */
%iparm[4] = 0;   /* No user fill-in reducing permutation */
%iparm[5] = 0;   /* Write solution into x */
%iparm[6] = 0;   /* Not in use */
%iparm[7] = 0;   /* Max numbers of iterative refinement steps */
%iparm[8] = 0;   /* Not in use */
%iparm[9] = 13;  /* Perturb the pivot elements with 1E-13 */
%iparm[10] = 1;  /* Use nonsymmetric permutation and scaling MPS */
%iparm[11] = 0;  /* Not in use */
%iparm[12] = 0;  /* Not in use */
%iparm[13] = 0;  /* Output: Number of perturbed pivots */
%iparm[14] = 0;  /* Not in use */
%iparm[15] = 0;  /* Not in use */
%iparm[16] = 0;  /* Not in use */
%iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
%iparm[18] = -1; /* Output: Mflops for LU factorization */
%iparm[19] = 0;  /* Output: Numbers of CG Iterations */
%maxfct = 1;     /* Maximum number of numerical factorizations. */
%mnum = 1;       /* Which factorization to use. */
%msglvl = 0;     /* Print statistical information in file */
%error = 0;      /* Initialize error flag */

elseif comstr(Cam,'solve')

     k=varargin{2}; b=varargin{3};
     if issparse(b); b=full(b);end
     param=k.data;
     q=pardiso('solve',b,k.ty(2),param);
     ks=q;

elseif comstr(Cam,'clear')  % - - - - - - - - - - - - - - - - - - - - - -

     k=varargin{2};
     if k.ty(2)<0; pardiso('clear',0);
     else; pardiso('clear',k.ty(2),k.data);
     end
     return;
elseif comstr(Cam,'silent')  % - - - - - - - - - - - - - - - - - - - - - -

 warning('Nothing done');

end % - - - - - - - - - - - - - - - - - - - - - -


if isempty(param); elseif param(104)<0;
 st1={-1 'input inconsistent';
     -2 'not enough memory';
     -3 'reordering problem';
     -4 'zero pivot, numerical factorization problem'
     -5 'unclassified (internal) error'
     -6 'preordering failed (matrix types 11, 13 only)'
     -7 'diagonal matrix problem'};
 try; 
   st1=st1{abs(double(param(24))),2};
 catch
   try; st1=st1{abs(double(param(104))),2}; 
   catch
    st1=sprintf('%i (%i) error code',double(param([24 104])));
   end
 end
 error(st1);
end
