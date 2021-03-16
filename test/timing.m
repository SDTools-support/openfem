function out =timing(index,level)

% timers[0] =  mxCreateDoubleMatrix(1, 1, mxREAL);
%timers[1] =  mxCreateDoubleMatrix(1, 1, mxREAL);
%*mxGetPr(timers[1])=(double)matlabmsglvl; 
% *mxGetPr(timers[0])=0.;mexCallMATLAB(0,(mxArray**)NULL,2,&timers, "timing");
%
%
%

persistent T TIMES

if isempty(T) T=cputime; TIMES=[]; end

if nargout==1 out=TIMES; return; end

if nargin==0
 if isempty(TIMES) 
   T=cputime;
 else disp(TIMES); TIMES=[]; end
elseif index==0
 T=cputime;TIMES=[];
else 
 if length(TIMES)<index TIMES(index)=0; end
 TIMES(index)=TIMES(index)+cputime-T; T=cputime;
end
