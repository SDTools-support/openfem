function [mode,freq,kd]=fe_eig(m,k,opt,mdof,adof,mode)

%FE_EIG compute and normalize normal modes from system matrices
%
%	Syntax : [mode,freq] = fe_eig(m,k)
%		 [mode,freq] = fe_eig(m,k,opt,mdof,adof,imode)
%
%	[MODE,FREQ]=FE_EIG(M,K) computes the eigensolution of the undamped
%	system equations and mass-normalizes the normal modes MODE
%
%	A number of options can be specified using FE_EIG(M,K,OPT). The row
%	of options contains
%	  [method nm shift print thres MaxIte]
%	with
%	  method -full matrix methods (for small or reduced problems) are
%		  0 : SVD based method (obsolete), 2 : cleaned up results of 
%		  EIG (default)
%		 - partial solutions for large sparse matrices
%		  1: subspace iteration (this method is VERY BASIC, much better
%                    solvers are provided in SDT)
%	  nm      number of modes to be found
%	  shift   mass shift value needed for cases with RIGID BODY modes
%		  (the square of the first flexible frequency is a good number)
%		  For iterations without mass shift but rigid body modes use
%		  OPT(3) = Inf
%	  print   level of printout
%	  thres   threshold for convergence of modes (1e-5 default)
%	  MaxIte  maximum number of iterations
%
%	MDOF and ADOF can be used to restrict the active degrees of freedom.
%	  See FE_C.
%
%	IMODE can provide an initial guess of modeshapes for iterative
%	      procedures
%
%
%	see also FE_MK

%	Etienne Balmes
%       Copyright (c) 1990-2012 by SDTools
%       All Rights Reserved.

%#ok<*NASGU,*ASGLU,*NOSEM>

time=cputime;cost(1)=0; in2=[]; iDOF=[]; st1='';

% input sorting section - - - - - - - - - - - - - - - - - - -
if isstruct(m) % deals with call : fe_eig(model,...)

  if nargin<3 st='';   Case=[];
  elseif isstruct(opt) Case=opt;st='';
  else st=opt; Case=[]; end

  model=m; if nargin>1 opt=k; else opt=2; end

  % determine case
  if isfield(model,'m')&&isfield(model,'k')  m=model.m; k=model.k;
  else 
    model=fe_mk(model);
    m=model.K{1};k=model.K{2};
    Case=fe_case(model,'gett');
    mdof=Case.DOF;
  end

  [mode,freq] = fe_eig(m,k,opt,mdof);
  if nargout==1 mode=struct('def',mode,'DOF',mdof,'data',freq); end
  return
end


if nargin == 2      opt = [2 size(m,2) 0 0]; end
if length(opt) < 6  opt = [opt(:)' zeros(1,6-length(opt))]; end
if opt(5)==0	    opt(5) = 1e-6*size(k,1); end
if opt(6)<1	    opt(6) = 99; end

if isempty(m) % matrices on file
  if isstr(k)  st=k; elseif isstruct(k) st=[k.file '_mk']; end
  load(st,'m','k');
  if opt(4)>0 disp(['Loading m,k matrices from ' st]);end
end

if opt(2)>size(m,1) opt(2)=size(m,1); end
if opt(1)>100 opt(7)=1; opt(1)=rem(opt(1),100); else opt(7)=0; end
if nargin~=6 mode = []; end
if opt(2)>size(m,1)/4 && any(opt(1)==[1 3 4 5 6]) opt(1)=2;
  if opt(4)>9  disp('Using fe_eig with opt(1)==2 for a small problem');end
end

nDOF0 = size(m,1);

if opt(3)==-1
  opt(3)=sqrt(real(sum(diag(k))/sum(diag(m))))/prod(size(m))/10;
end
if opt(3)~=0 && ~isinf(opt(3)) rDOF=find(diag(k)+opt(3)*diag(m));
else                          rDOF=find(diag(k)); end

if nargin>4  % DOF selection 
   if ~isempty(mdof) && length(mdof)~=size(m,1) 
           error('M and MDOF are incompatible in size');
   end
   if ~isempty(adof) [adof,i1] = fe_c(mdof(rDOF),adof); rDOF = rDOF(i1);end
   if isempty(rDOF) error('No active DOF selected'); end
end

if length(rDOF)~=nDOF0 opt(8)=0;
 if opt(4)>1 disp('fe_eig: deleting some DOF with fixed/empty column/rows');end
else opt(8)=1; end

if ~isreal(k) warning('taking real part of k');k=real(k); end
if ~isreal(m) warning('taking real part of m');m=real(m); end

if any(opt(1)==[0 .5 2])  % DOF elimination and mass-shift (full matrix)

  if opt(3)~=0 && ~isinf(opt(3))
         k1 = k + opt(3) * m;
         if opt(4)>9 disp(sprintf('Mass shift %7.3e',opt(3)));end
  else   k1 = k; end
  if length(rDOF)~=nDOF0 m1 = m(rDOF,rDOF); k1=k1(rDOF,rDOF);
  else                   m1 = m; end

elseif any(opt(1)==[1 3 4 5 6])  % DOF elimination and mass-shift (sparse)

  if ~opt(8) k = k(rDOF,rDOF);end  % eliminate DOFs
  if opt(7)~=1 % bandwidth
    if opt(1)==5 && size(k,1)>1000
      warning('you should do the bandwidth optimization before');
    end
    iDOF = symrcm(k);
    k = k(iDOF,iDOF); m=m(rDOF(iDOF),rDOF(iDOF));
  else % no bandwith optimization
    if ~opt(8)
      m = m(rDOF,rDOF);
    end  % eliminate DOFs
    iDOF = [1:length(rDOF)];
  end
  if opt(3)~=0 && ~isinf(opt(3))
    k = k + opt(3) * m;
    if opt(4)>9 disp(sprintf('Mass shift %7.3e',opt(3)));end
  end

  if     opt(1)==3||(opt(1)==1&&size(k,1)<1e3&&isfinite(opt(3)))
                   kd = ofact(k,'lu');
  elseif opt(1)==4 kd = ofact(k,'ch');
  else             kd = ofact(k);  end
  if size(k,1)>5000&&opt(4)>1 display(kd,0); end
  cost(2)=cputime-time;time=cputime;

end


if opt(1) == 0
% ----------------------------------------------------------------------------
% Updated 04/21/98
% computation of all vectors
	
m1 = chol(m1); m1 = inv(m1); k1 = m1' * k1 * m1;

[cmode,cfreq,v] = svd(full(k1));
cfreq = flipud(diag(cfreq))-opt(3);
cmode = fliplr(m1 * cmode);

ind = find(cfreq<0); opt(2)=size(cmode,2);
if ~isempty(ind) cfreq(ind) = cfreq(ind)*0; end

elseif opt(1) == 0.5
% ----------------------------------------------------------------------------
% Updated 04/21/98
% svd based computation of all vectors  

%k1=triu(k1)+triu(k1,1)';m1=triu(m1)+triu(m1,1)';
k1=(k1+k1')/2;m1=(m1+m1')/2;

[k1,s,v]=svd(full(k1)); i1 = find(diag(s)<s(1)*eps*size(s,1));
if ~isempty(i1) % there are rigid body modes
     i2 = 1:size(k1,1); i2(i1)=zeros(1,length(i1));i2=find(i2);
     s=diag(s);
     if ~isempty(i2) k1(:,i2) = k1(:,i2)*diag(s(i2).^(-.5)); end
     [cmode,cfreq,v] = svd(full(k1' * m1 * k1));

     cfreq = diag(cfreq).^(-1);
     i2 = find(isfinite(cfreq)); cmode(:,i2)=cmode(:,i2)*diag(sqrt(cfreq(i2)));
     cmode = k1 * cmode;
else
     k1 = k1*diag(diag(s).^(-.5));
     [cmode,cfreq,v] = svd(full(k1' * m * k1));
     cfreq = diag(cfreq).^(-1); cmode = k1 * cmode * diag(sqrt(cfreq));
end

if ~isinf(opt(3)) && opt(3)~=0
   cfreq=cfreq-opt(3);ind = find(cfreq<0);
   if ~isempty(ind) cfreq(ind) = cfreq(ind)*0; end
end

elseif opt(1) == 2
% ----------------------------------------------------------------------------
% Updated 04/21/98
% eig based computation of all vectors with cleaning up

if nnz(m1)==size(m1,1) && ~any(m1-eye(size(m1)))
  [u,s]=eig(full(k1)); s=diag(s);
else 
  k1=(k1+k1')/2;m1=(m1+m1')/2;
  [u,s]=eig(full(k1),full(m1)); s=diag(s);
end

i1=find(isfinite(s)); [s,i2]=sort(abs(s(i1))); u=u(:,i1(i2));
% eliminate small imaginary components
i1=find(abs(imag(u))<size(u,1)^2*eps*abs(u));u(i1)=real(u(i1));
i1=find(abs(imag(u))<eps);u(i1)=real(u(i1));

% complex eigenvalues (reference is i3)
i1=[diag(k1) diag(m1)];i3=find(i1(:,2));i3=mean(abs(i1(i3,1)./i1(i3,2)));
i1=find(abs(imag(s))>sqrt(eps)*i3);
if ~isempty(i1)  i2=i1(find(abs(s(i1))>sqrt(eps)*i3));
 if ~isempty(i2) disp('FE_EIG WARNING: problem with non-real eigenvalues');end
end
s=real(s);
i2=find(s<0&abs(s)<sqrt(eps)*i3);s(i2)=i2*0; %negative zero eigenvalue
cfreq=s;

% make multiple eigenvectors real
i2=[-1;find(diff(abs(s))<eps*abs(s(length(s))));Inf];
i3=[find(diff(i2)~=1)];
for j1=1:length(i3)-1
 v=u(:,i2(i3(j1)+1)+[0:i3(j1+1)-i3(j1)]); v=[real(v) imag(v)];
 [u1,s1,v1]=svd(v'*m1*v);i1=1:size(v,2)/2;
 s1=diag(s1);u1=v1(:,i1)*diag(s1(i1).^(-.5));
 kr=u1'*v'*k1*v*u1; [v1,s1]=eig(kr+kr');
 v=v*u1*v1;
 u(:,i2(i3(j1)+1)+[0:i3(j1+1)-i3(j1)])=v;
end

for j1=1:size(u,2)  v=u(:,j1);
  if ~isreal(v)
    v=[real(v) imag(v)];%error('complex mode')
    [u1,s1,v1]=svd(v'*m1*v);s1=diag(s1);
    v=v*u1(:,1)/sqrt(s1(1,1));
  else v=v/sqrt(v'*m1*v); end % mass normalization
  u(:,j1)=v;
end

cmode=u;

if ~isinf(opt(3)) && opt(3)~=0
   cfreq=cfreq-opt(3);ind = find(cfreq<0);
   if ~isempty(ind) cfreq(ind) = cfreq(ind)*0; end
end
if nargout==3 kd = ofact(k1,'lu'); end

elseif any(opt(1) == [1 ])
% ----------------------------------------------------------------------------
% subspace iteration
% Updated 04/21/98

if ~isempty(mode) mode = mode(rDOF(iDOF),:); end

if isinf(opt(3)) % detection of rigid body modes

   [l,u] = lu(sparse(k),0);

   if imag(opt(3))==0
     in2 = find(abs(diag(u))<norm(u,'fro')*eps*size(u,1));
   else [i1,i2]=sort(abs(diag(u)));in2=i2(1:abs(imag(opt(3))));end
   aDOF =1:size(k,1);
   if ~isempty(in2)
	aDOF(in2)=zeros(1,length(in2));aDOF=find(aDOF);

        rb=spalloc(size(k,1),length(in2),1);
        rb(in2,:)=eye(length(in2)); rb(aDOF,:) = -u(aDOF,aDOF)\u(aDOF,in2);

	u1 = full(rb'*m*rb);i2=find(diag(u1)<=eps);
	if ~isempty(i2)
               i1 = find(diag(u1)>eps);
	       rb = rb(:,i1);u1=u1(i1,i1);
               aDOF=sort([aDOF(:);in2(i2)]);in2=in2(i1);
        end
        if ~isempty(rb)
		[u1,s,v]=svd(u1);  s = diag(s);
		if length(s)==1 i1=find(s>eps); else i1=find(s>s(1)*eps); end
		rb = rb*v(:,i1)*diag(s(i1).^(-.5));
	end

        if ~isempty(rb)
	  rbm = rb'*m; opt(2)=max(0,opt(2)-size(rb,2));
	  if opt(4)> 0
	    disp(sprintf('Detected %i rigid body modes',size(rb,2)));
	  end
 	  clear u; l = chol(k(aDOF,aDOF));
	elseif ~isempty(in2) error('RB detection');
	else  clear u;
         if size(k,1)>1e3 l = ofact(k); else [l,u] = lu(k,0); end
        end
   end
   opt(3)=0;

end % detection of RB modes

    % initialisation of subspace
	
	if opt(2)==0 opt(2) = size(k,1)/2; end
	nM = min([2*opt(2),opt(2)+12,size(m,1)]);

	if opt(4)> 0   disp(sprintf(...
	    'Subspace iteration with %i modes and %i DOFs',nM,size(m,1)));
	end
	i1=diag(m);i2=diag(k);ind=find(i1>mean(i1)*sqrt(eps));
	[i1,i2] = sort(i2(ind)./i1(ind));  ind = ind(i2);
        if nM>length(ind) ind=1:nM;warning('computing too many modes');end
	if isempty(mode)
		mode = zeros(size(m,1),nM);  mode(ind(1:nM),:) = eye(nM);
	elseif size(mode,2) < nM
		i2 = size(mode,2);
		mode(ind(1:nM-i2),i2+1:nM) = eye(nM-i2);
	else	nM = size(mode,2);
	end

	Ite = 1;  freq = zeros(nM,1);

	cmode=[];res=[];cfreq=[]; cind = 1:nM;

while Ite < opt(6)+1 && Ite>0  % iterations - - - - - - - - - - - - - - - - - -

   mode(:,size(mode,2)) = rand(size(mode,1),1);
   mode = m*mode;
   if ~isempty(res) mode(:,size(mode,2)-1) = res;res=[]; end

   if isempty(in2) % no RB
             mode = kd\mode;
   else
             mode = mode - rbm' * (rb'*mode);
	     mode(in2,:) = zeros(length(in2),size(mode,2));
	     mode(aDOF,:) = l \ (l' \ mode(aDOF,:));
             mode = mode - rb * (rbm*mode);
   end

   if ~isempty(cmode)   mode = mode - cmode * (cmode'*m*mode);  end
   mr = mode'*m*mode;

   [bb,aa] = chol(mr);
   if aa~=0	[u1,s,v]=svd(full(mr));  mode = mode*v*diag(diag(s).^(-.5));
   else		mode = mode/bb;   end

   % partial eigenvalue solution ----------------------------------------
   if (rem(Ite,5)==0 || (Ite==1&&nargin==6) || Ite==opt(6))

   oldF = freq;
   [mdr,freq,v] = svd(full(mode'*k*mode)); [freq,ind] = sort(sqrt(diag(freq)));
   mode = mode * mdr(:,ind);

   if max(oldF)>0
	ca=[];for j1 = 1:size(mode,2)
           res = m*(-freq(j1).^2*mode(:,j1))+k*mode(:,j1);
	   ca(j1,1)= norm(res)/norm(mode(:,j1));
	   if ca(j1,1)>=opt(5)||abs(oldF(j1)-freq(j1))>1e-4*abs(freq(j1))
             break;
           end
	end

	ind = 1:j1-1; cind = j1:size(mode,2);
	%ind = ind(1:min([find(diff([0;ind])>1);length(ind)+1])-1);
	%cind=ones(size(mode,2),1); cind(ind)=ind*0;cind=find(cind);

        ca(1:j1,2) = abs(oldF(1:j1)-freq(1:j1))./abs(freq(1:j1));

	if opt(4) > 10   disp(sprintf( ...
             '\n Frequency   converged   residual   deltaFreq \n'));
	     disp(sprintf(' %7.3e       %i       %7.1e     %7.1e\n', ...
	     [sqrt(abs(freq(1:j1).^2-opt(3))) ca(:,1)<opt(5) ca]'));
	end

        if isempty(ind)
	   %i1 = fix((size(mode,2)+size(cmode,2)-opt(2))/2);
	   %mode(:,size(mode,2)+[-i1+1:0])=rand(size(mode,1),i1);
        else
             cmode = [cmode mode(:,ind)];cfreq=[cfreq;freq(ind)];
	     mode = mode(:,cind);freq = freq(cind);
        end 
	if size(cmode,2)>=opt(2)
	   [cfreq,i2]=sort(cfreq);cmode=cmode(:,i2);
           if min(freq)>max(cfreq(1:opt(2)))
	     disp(['Converged after ' num2str(Ite) ' iterations.']);Ite=-Ite;
           end
        end

   end % of non first eigenvalue case
   if opt(4)> 0 
     disp(sprintf('%i iterations, %i converged modes',abs(Ite),size(cmode,2)));
   end
   end % of fifth iteration treatment

          oldF = freq;
          Ite = Ite+1;

   end % of iterations

   if Ite>0
	  warning('max iteration reached without convergence');
	  [cfreq,i2]=sort(cfreq);cmode=cmode(:,i2);
          if isempty(cmode) cmode = mode; end
   end

   if opt(3)~=0&&~isinf(opt(3)) cfreq = abs(cfreq.^2-opt(3));
   else cfreq=cfreq.^2; end

   if opt(4)> 0
     cost(3)=cputime-time;cost(2)=cost(2)+cost(3);
     disp(sprintf( ...
      'Total time %i min. %i sec., iterations %i min. %i sec', ...
      [fix([cost(2)/60 rem(cost(2),60) cost(3)/60 rem(cost(3),60)])]))
   end

elseif any([3 4 5 6]==opt(1))
% ---------------------------------------------------------------------------- 
% IRA/Sorensen solver 

%sta =randn('state');randn('state',0) 
if opt(4)>0 disp(sprintf(' IRA/Sorensen solver with %i DOFs',size(m,1)));end 
opts=struct('issym',0,'isreal',1,'tol',min(opt(5)*1e-6,eps/2), ...
      'maxit',opt(2)*2,'p',opt(2)*2,'disp',0); 
%opts.v0=rand(size(k,1),1); 

afun=str2func('arp_eig'); 
[cmode,d]=eigs(afun,size(k,1),opt(2),'lm',opts,m,kd); 

if ~isreal(cmode) || ~isreal(d)
 r1=[real(cmode) imag(cmode)];r1=r1(:,find(sum(abs(r1))));
else r1=cmode; end

[mdr,d1]=fe_eig(r1'*m*r1,r1'*k*r1,2);

cfreq=abs(d1.^2-opt(3));cmode=r1*mdr;

if opt(4)> 0 
   cost(3)=cputime-time;cost(2)=cost(2)+cost(3); 
   disp(sprintf( ... 
    '\nTotal time %i min. %i sec., iterations %i min. %i sec', ... 
      [fix([cost(2)/60 rem(cost(2),60) cost(3)/60 rem(cost(3),60)])])) 
end 
if nargout<3 ofact('clear',kd);end 
%randn('state',sta) 

% ----------------------------------------------------------------------------
else error('not a known type of eigenvalue computation');
end
% ----------------------------------------------------------------------------
% do some cleaning up
if opt(2)==0	   opt(2) = size(cmode,2); end
 
opt(2) = min(opt(2),size(cmode,2));
cmode = cmode(:,1:opt(2));  freq = sqrt(cfreq(1:opt(2)));

if ~isempty(in2)
	freq=[zeros(size(rb,2),1);freq]; cmode = [rb cmode]; 
end

mode=zeros(nDOF0,size(cmode,2));

if isempty(iDOF)	mode(rDOF,:) = cmode;
else			mode(rDOF(iDOF),:) = cmode; end

if nargout<2 mode=freq; end

% ---------------------------------------------------------------------
function y = arp_eig(x,m,kd);

if nargin==1
  assignin('base','x',x);
  evalin('base','x=kd\ (m*x);');
  y=evalin('base','x');
else
   y=kd\ (m*x); 
end

