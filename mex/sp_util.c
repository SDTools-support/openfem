#include "mex.h" 
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/*  Etienne Balmes, Jean Michel Leclere, Robert Cimrman */
/*  Claire Delforge                                     */

/* ------------------------------------------------------------------*/
/*     SCILAB/Matlab MACROS                                          */
/* ------------------------------------------------------------------*/
/* #include "../mex/of_def.h" */
#ifdef MatlabVER

#include  "matrix.h"
#include "pmat.c"
void SetPrToNull(mxArray *prhs)
{  mxSetPr(prhs,NULL); }


#else

void SetPrToNull(mxArray *prhs)
{
  double * tmp_pr;
  tmp_pr = mxGetPr(prhs);
  tmp_pr = NULL;
}
#endif

/* obsolete with scilab 3
if !defined(mxArrayToString)
 char* mxArrayToString(mxArray *prhs) {

  int buflen, status;
  char* buf;

  buflen = (mxGetM (prhs) * mxGetN(prhs)) + 1;
  buf = mxCalloc (buflen, sizeof(char));
  status = mxGetString (prhs, buf, buflen);
  return(buf);
}
*/ 


/* ------------------------------------------------------------------*/
/*        FullRealloc                                                  */
/* ------------------------------------------------------------------*/
mxArray* FullRealloc ( mxArray *full, mwSize newM, mwSize newN ) {

	double *ptr;
	double   *newptr;
	int oldN, oldM, i, j;
	mxArray* out; 
	
	if (mxIsComplex(full)) {
#if MatlabVER >= 904
		out = mxCreateDoubleMatrix(newM, newN, mxCOMPLEX);
		memcpy(mxGetComplexDoubles(out), mxGetComplexDoubles(full),
			mxGetNumberOfElements(full) *2* sizeof(double));
#else
		out=mxCreateDoubleMatrix(newM,newN,mxCOMPLEX);
 	    memcpy(mxGetPr(out),mxGetPr(full), 
                mxGetNumberOfElements(full)*sizeof(double));
 	    memcpy(mxGetPi(out),mxGetPi(full), 
                mxGetNumberOfElements(full)*sizeof(double));
#endif
	} else { 
		out=mxCreateDoubleMatrix(newM,newN,mxREAL);
 	    memcpy(mxGetPr(out),mxGetPr(full), 
                mxGetNumberOfElements(full)*sizeof(double));
	}
    return(out);
}


/* ------------------------------------------------------------------*/
/*      C: function :       heapsort, two *double fields are sorted  */
/* ------------------------------------------------------------------*/
void hpsort(const int n, double ra[], double rb[]) {
  unsigned long  i, ir, j, l;
  double         rra, rrb;

  ra--; rb--;                                   /* just for shifting */

  if (n<2)  return; 
  l = (n >> 1)+1; 
  ir = n;

  for (;;) {
    if (l>1) { rra = ra[--l]; rrb = rb[l];}
    else {
      rra = ra[ir];    rrb = rb[ir];
      ra[ir] = ra[1];  rb[ir] = rb[1];
      if (--ir == 1) {
        ra[1]  = rra;  rb[1]  = rrb;
        break;
      }
    }
    i = l; 
    j = l+l;
    while (j<=ir) {
      if ( j<ir && ra[j]<ra[j+1] ) j++;
      if ( rra<ra[j] ) {
        ra[i]=ra[j]; rb[i]=rb[j];
        i=j;
        j <<= 1;
      }
      else j = ir+1;
    }
    ra[i] = rra;  rb[i] = rrb;
  }
  ra++; rb++;                                   /* just for shifting */
}

/* ------------------------------------------------------------------  */
/*      C: function :       heapsort, three *double fields are sorted  */
/* ------------------------------------------------------------------  */
void hpsortc(const int n, double ra[], double rb[], double rc[]) {
  unsigned long  i, ir, j, l;
  double         rra, rrb,rrc;

  ra--; rb--; rc--;         /* just for shifting */

  if (n<2)  return; 
  l = (n >> 1)+1; 
  ir = n;

  for (;;) {
    if (l>1) { rra = ra[--l]; rrb = rb[l]; rrc = rc[l];}
    else {
      rra = ra[ir];    rrb = rb[ir];   rrc = rc[ir];
      ra[ir] = ra[1];  rb[ir] = rb[1]; rc[ir] = rc[1];
      if (--ir == 1) {
        ra[1]  = rra;  rb[1]  = rrb; rc[1]  = rrc;
        break;
      }
    }
    i = l; 
    j = l+l;
    while (j<=ir) {
      if ( j<ir && ra[j]<ra[j+1] ) j++;
      if ( rra<ra[j] ) {
        ra[i]=ra[j]; rb[i]=rb[j]; rc[i]=rc[j];
        i=j;
        j <<= 1;
      }
      else j = ir+1;
    }
    ra[i] = rra;  rb[i] = rrb; rc[i] = rrc;
  }
  ra++; rb++;rc++;         /* just for shifting */
}






#ifdef SDT_ADD
#include "../../sdtdev/mex50/sp_util_pre.c"
#include "../../openfem/mex/sp_util_subs.c"
#else
#include "sp_util_subs.c"
#endif

static double EPSL=1.e-6;
static double OpenFEMDIAG=0;

/* ------------------------------------------------------------------------ */
/*             The gatewayfunction                                          */
/* ------------------------------------------------------------------------ */

void mexFunction (int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]) {
char            *buf,sbuf[32] ;
int             buflen;

  int   i, i1;



  if (nrhs==0)  {
  double *a;
  plhs[0] = mxCreateDoubleMatrix (1, 1, mxREAL);
  a = mxGetPr (plhs[0]);  
#ifdef SDT_ADD
  *a=5.0008;
#else
  *a=1.0002;
#endif
  return; 
  }

  /* if (mxIsChar(prhs[0])) { mxGetString(prhs[0],buf,16); } else {buf[0]='\0';}*/
  buf = mxArrayToString(prhs[0]);
if (mxIsStruct(prhs[0])) { buf=(char*)mxGetFieldNameByNumber(prhs[0],0); buflen=-1;
} else { 
   buflen = 31;// ((int)mxGetM (prhs[0]) * (int)mxGetN(prhs[0])) + 1;
   buf=sbuf; //buf = (char*)mxCalloc (buflen, sizeof(char));
   if (mxGetString(prhs[0],buf,buflen)) mexErrMsgTxt("String conversion problem");
}
 

/* if (strcmp("diag",buf)&&strcmp("epsl",buf)) mexEvalString("disp('sp_util');dbstack;disp('__');");*/

/*----------------------------------------------------------- ismex */
  if (!strcmp("ismex",buf))  {
  double *a;

  plhs[0] = mxCreateDoubleMatrix (1, 1, mxREAL);
  a = mxGetPr (plhs[0]);
  *a=1.0;
/* #setinput,In,Val,Offset --------------------------------------------- */
} else if (!strcmp("setinput",buf))  {

  double   *inPr, *val, *opt, *iin=NULL, *ival=NULL;
  int     jrow, jcol, ioff, voff, typ;
  size_t M;
  mxArray* ToSet;
  char    *buf1;
  mwIndex jrow2;

  typ=0;
  if (nrhs==1) { /* 2023 copy bytes from used in nas2up */
   int* dims;
   if (!strcmp(mxGetFieldNameByNumber(prhs[0],1),"from")) {
     char* to;
     char* from;
     // copy bytes from to 
     to=(char*)mxGetData(mxGetFieldByNumber(prhs[0],0,0));
     from=(char*)mxGetData(mxGetFieldByNumber(prhs[0],0,1));
     dims=(int*)mxGetData(mxGetFieldByNumber(prhs[0],0,2));
     to+=dims[1];
     from+=dims[2];
     memcpy(to,from,(mwSize)dims[0]);
     return;
   }
 }
  if (nrhs>5) {
   buf1 = mxArrayToString( prhs[5] );
   if      (!strcmp("+",buf1))   typ=1;
   else if (!strcmp("-",buf1))   typ=2;
   else if (!strcmp(".*",buf1))  typ=3;
   else if (!strcmp("./",buf1))  typ=4;
   else if (!strcmp("*",buf1))   typ=5;
   else if (!strcmp("/",buf1))   typ=6;
   else if (!strcmp("=",buf1))   typ=0; /* default if not given */
   mxFree(buf1);
  }
 ToSet=(mxArray*)prhs[1]; 
 if (nrhs==4) {
     /* 4 inputs fill data from starting position - (out,in,start) - - - - - - - - - - - - - - - - - - - - - -*/
	  mwSize i1, i2;
	  of_ptrdiff two = 2, one = 1, dM;
      size_t irow;
      opt=mxGetPr(prhs[3]);  irow=(size_t)opt[0];
      M=mxGetNumberOfElements(prhs[2]);
      pmatGetPrAndSize(ToSet,&inPr,&i1);
      
      if (irow==-1) { /* realloc rows */
		if (nlhs==0) mexErrMsgTxt("Realloc requires output");
        if (mxGetN(prhs[3])==1) { return;} else irow=(size_t)opt[1];
        while (irow+M>mxGetNumberOfElements(ToSet)) {
          i2=(irow+M)/mxGetN(ToSet)+1; i1=mxGetM(ToSet); if (i2>i1) i1=i2; 
		  if (i1<4000) i1=i1+1000; else i1=(mwSize)((double)i1*1.25);
          ToSet=FullRealloc(ToSet,i1,mxGetN(ToSet));
        }
        pmatGetPrAndSize(ToSet,&inPr,&i1);
		plhs[0]=ToSet;
      } else if (irow==-2) { /* realloc cols */
		if (nlhs==0) mexErrMsgTxt("Realloc requires output");
        if (mxGetN(prhs[3])==1) {return;} else irow=(size_t)opt[1];
        while (irow+M>mxGetNumberOfElements(ToSet)) {
          i2=(irow+M)/mxGetM(ToSet)+1; i1=mxGetN(ToSet); if (i2>i1) i1=i2; 
		  if (i1<4000) i1=i1+1000; else i1=(mwSize)((double)i1*1.25);
          ToSet=FullRealloc(ToSet,mxGetM(ToSet),i1);
		}
        pmatGetPrAndSize(ToSet,&inPr,&i1);
		plhs[0]=ToSet;
      } else if (irow==-3) { /* GetInput */
        int step; double *OUT,*IN;
        if (mxGetN(prhs[3])==1) {return;}
        step=(int)opt[2];IN=mxGetPr(ToSet);IN+=(size_t)opt[1];
        OUT=mxGetPr(prhs[2]);
        if (opt[1]+step*(mxGetNumberOfElements(prhs[2]))-1>
                mxGetNumberOfElements(ToSet)) 
              mexErrMsgTxt("GetInput overflow");
        for (jcol=0;jcol<mxGetNumberOfElements(prhs[2]);jcol++) {
          OUT[jcol]=*(IN+jcol*step);
        }
        return;
      } else {
		  if (mxGetN(prhs[3])==2) {irow=(size_t)opt[1];}
		  if (i1-M<irow|i1-M>42949667296) mexErrMsgTxt("Setinput overflow");
	  }
      if (mxIsInt32(prhs[2])) { /* memcopy(in,val,size) */
         if (!mxIsInt32(ToSet)) mexErrMsgTxt("Type Mismatch");
          memcpy((int*)mxGetData(ToSet)+irow,mxGetData(prhs[2]),M*sizeof(int32)); 
	  }
	  else if (mxIsDouble(prhs[2])) {
         #if MatlabVER >= 904
		  if (!mxIsComplex(prhs[2])) {
			  if (!mxIsComplex(ToSet)) memcpy(inPr + irow, mxGetData(prhs[2]), M * sizeof(double));
			  else {dM = M;of_dcopy(&dM, (double*)mxGetData(prhs[2]), &one, inPr + irow, &two);
			  }
		  } else {
			  double* r1;
			  r1 = (double*)mxGetData(prhs[2]);
			  if (mxIsComplex(ToSet)) memcpy(inPr + 2*irow, mxGetData(prhs[2]), 2 * M * sizeof(double));
			  else mexErrMsgTxt("data must be complex");
	     }
         #else
		  memcpy(inPr + irow, mxGetData(prhs[2]), M * sizeof(double));
         #endif
      } else if (mxIsClass(prhs[2],"pmat")) {
          memcpy(inPr+irow,pmatGetData(prhs[2]),1*sizeof(double)); 
      } else if (mxIsChar(prhs[2]) && mxIsChar(ToSet)) {
          if (M>mxGetNumberOfElements(ToSet)) 
                M=mxGetNumberOfElements(ToSet);
          memcpy((mxChar*)mxGetData(ToSet)+irow,
                  mxGetData(prhs[2]),M*sizeof(mxChar)); 
      } else if (mxIsCell(prhs[2]) && mxIsCell(ToSet)) {
/*          for (j1=0;j1<M;j1++) {
//              mxSetCell(ToSet,(mwIndex)(j1+irow),
//                      mxGetCell(prhs[2],(mwIndex)j1));
          }*/
      } else { 
          mexPrintf("a=%s b=%s",mxGetClassName(ToSet),
                  mxGetClassName(prhs[2]));
          mexErrMsgTxt("Not a supported type of input");
      }
      /* both complex -assign */
     #if MatlabVER < 904
	 if (mxIsComplex(prhs[2])) {
       if (!mxIsComplex(ToSet)) mexErrMsgTxt("data must be complex");
         memcpy(mxGetPi(ToSet)+irow,mxGetPi(prhs[2]),M*sizeof(double));
     }
     #endif
	  if (mxGetN(prhs[3])==2) {opt[1]+=(double)M;} else opt[0]+= (double)M;
      return;
    /*  -------------------------------------------------------------------------------
	   with 5 inputs sp_util('setinput',IN,val,[ii],'IN') ---------------------------*/ 
    }  else if ( nrhs==5 & mxIsChar(prhs[4]) ) {
      int i1, i2;
      size_t irow;
	  double* pToSet;
      opt=mxGetPr(prhs[3]);  irow=(size_t)opt[0];if (mxGetN(prhs[3])==1) {return;}
      M=mxGetNumberOfElements(prhs[2]); buf1=NULL;
      if (opt[0]==-1) { /* realloc rows */
        irow=(size_t)opt[1];
		if (irow+M<=mxGetNumberOfElements(ToSet)) { 
	    } else {
		 buf1 = mxArrayToString(prhs[4]); ToSet=mexGetVariable("caller",buf1);
         while (irow+M>mxGetNumberOfElements(ToSet)) {
	      i2=(irow+M)/mxGetN(ToSet)+1; i1=mxGetM(ToSet); if (i2>i1) i1=i2; 
		  if (i1<4000) i1=i1+1000; else i1=(int)((double)i1*1.25);
          ToSet=FullRealloc(ToSet,i1,mxGetN(ToSet));
         }
		}
      } else if (opt[0]==-2) { /* realloc cols */
        irow=(size_t)opt[1];
		if (irow+M<=mxGetNumberOfElements(ToSet)) { 
	    } else {
         /* if no realloc no need to reassign */ 
		 buf1 = mxArrayToString(prhs[4]); ToSet=mexGetVariable("caller",buf1);
         while (irow+M>mxGetNumberOfElements(ToSet)) {
          i2=(irow+M)/mxGetM(ToSet)+1; i1=mxGetN(ToSet); if (i2>i1) i1=i2; 
		  if (i1<4000) i1=i1+1000; else i1=(int)((double)i1*1.25);
          ToSet=FullRealloc(ToSet,mxGetM(ToSet),i1);
		 }
		}
      } else if (mxGetNumberOfElements(ToSet)-M<irow) mexErrMsgTxt("Setinput overflow");
      if (mxIsDouble(prhs[2])) {
#if MatlabVER >= 904
		  if (mxIsComplex(prhs[2])) memcpy(mxGetComplexDoubles(ToSet) + irow, mxGetComplexDoubles(prhs[2]), M*sizeof(double)*2);
		  else memcpy(mxGetPr(ToSet) + irow, mxGetData(prhs[2]), M * sizeof(double));
#else
		  if (mxIsComplex(prhs[2])) {memcpy(mxGetPi(ToSet)+irow,mxGetPi(prhs[2]),M*sizeof(double)); }
          memcpy(mxGetPr(ToSet)+irow,mxGetData(prhs[2]),M*sizeof(double)); 
#endif
	  } else mexErrMsgTxt("Not a supported type of input");
      if (mxGetN(prhs[3])==2) {opt[1]+=M;} else opt[0]+=M;
	  if (buf1!=NULL) {	mexPutVariable("caller",buf1,ToSet);mxFree(buf1);}
      return;
      
    /*  with 4 inputs (out,in,irow,icol) --------------------------------------------------- */ 
    }  else if ( mxIsDouble(ToSet) & mxIsDouble(prhs[2]) ) {
	int *irow,*icol;
    if (mxIsComplex(prhs[2])) {
       if (!mxIsComplex(ToSet)) mexErrMsgTxt("data must be complex");
       #if MatlabVER >= 904
	   inPr = (double*)mxGetComplexDoubles(ToSet); val = (double*)mxGetComplexDoubles(prhs[2]);
       #else
       inPr = mxGetPr(ToSet); val = mxGetPr(prhs[2]);
	   iin=mxGetPi(ToSet); ival=mxGetPi(prhs[2]);
       #endif
	}
	else { inPr = mxGetPr(ToSet); val = mxGetPr(prhs[2]); }
    if (!mxIsInt32(prhs[3]) || !mxIsInt32(prhs[4])) mexErrMsgTxt("Non int32 indices");
	irow=(int*)mxGetData(prhs[3]); icol= (int*)mxGetData(prhs[4]);
    
    /* check size for case where row and col indices given */
    if (     (irow[mxGetM(prhs[3])-1]>(int)mxGetM(ToSet))
         ||  (icol[mxGetN(prhs[4])-1]>(int)mxGetN(ToSet))    )
                                            mexErrMsgTxt("Wrong indices");
    M=mxGetM(prhs[2]);
    if (mxIsSparse(prhs[2])) mexErrMsgTxt("Setinput not implemented for sparse");
    if (typ==0) { /* set value */
		for (jcol = 0; jcol < mxGetN(prhs[2]); jcol++) {
			ioff = mxGetM(ToSet)*(icol[jcol] - 1); voff = M * jcol;
#if MatlabVER >= 904
	    if (mxIsComplex(prhs[2])) {
			for (jrow = 0; jrow < M; jrow++) {
				inPr[(irow[jrow] - 1 + ioff)*2] = val[(voff + jrow)*2];
				inPr[(irow[jrow] - 1 + ioff)*2+1] = val[(voff+jrow)*2+1];
			}
		} else {
				for (jrow = 0; jrow < M; jrow++) inPr[irow[jrow] - 1 + ioff] = val[voff + jrow];
		}

#else
	   for (jrow = 0; jrow<M; jrow++) inPr[irow[jrow] - 1 + ioff] = val[voff + jrow];
	   if (ival!=NULL){
        for (jrow=0;jrow<M;jrow++) iin[irow[jrow]-1+ioff] = ival[voff+jrow];
       }
#endif
      }
    } else if (typ==1) {  /* add to matrix */
      for (jcol=0;jcol<mxGetN(prhs[2]);jcol++) {
       ioff=mxGetM(ToSet)*(icol[jcol]-1); voff=M*jcol;
       for (jrow=0;jrow<M;jrow++) inPr[irow[jrow]-1+ioff] += val[voff+jrow];
      }
    } else if (typ==2) {  /* minus to matrix */
      for (jcol=0;jcol<mxGetN(prhs[2]);jcol++) {
       ioff=mxGetM(ToSet)*(icol[jcol]-1); voff=M*jcol;
       for (jrow=0;jrow<M;jrow++) inPr[irow[jrow]-1+ioff] -= val[voff+jrow];
      }
    } else if (typ==3) {  /* multiply matrix, by coordinates  */
      for (jcol=0;jcol<mxGetN(prhs[2]);jcol++) {
       ioff=mxGetM(ToSet)*(icol[jcol]-1); voff=M*jcol;
       for (jrow=0;jrow<M;jrow++) inPr[irow[jrow]-1+ioff] *= val[voff+jrow];
      }
    } else if (typ==4) {  /* divide matrix, by coordinates */
      for (jcol=0;jcol<mxGetN(prhs[2]);jcol++) {
       ioff=mxGetM(ToSet)*(icol[jcol]-1); voff=M*jcol;
       for (jrow=0;jrow<M;jrow++) inPr[irow[jrow]-1+ioff] /= val[voff+jrow];
      }
    } else if (typ==5) {  /* multiply  matrix, by scalar */
      M=mxGetN(prhs[3]);
      for (jcol=0;jcol<mxGetN(prhs[4]);jcol++) {
       ioff=mxGetM(ToSet)*(icol[jcol]-1);
       for (jrow=0;jrow<M;jrow++) inPr[irow[jrow]-1+ioff] *= *val; 
      }
    } else if (typ==6) {  /* divide matrix, by scalar */
      M=mxGetN(prhs[3]);
      for (jcol=0;jcol<mxGetN(prhs[4]);jcol++) {
       ioff=mxGetM(ToSet)*(icol[jcol]-1);
       for (jrow=0;jrow<M;jrow++)  inPr[irow[jrow]-1+ioff] /= *val; 
      }
    } /* type of double operations */

   } else if (mxIsCell(ToSet)) {
    opt=mxGetPr(prhs[3]);
    jrow2=(mwIndex)opt[0]-1;
    mxDestroyArray(mxGetCell(ToSet,jrow2));
    mxSetCell(ToSet,jrow2,(mxArray*)prhs[2]);
   } else mexErrMsgTxt("Not a supported setinput case");

/*----------------------------------------------------------- mwIndex */
}  else if (!strcmp("mwIndex",buf)) {

    if (nlhs==1)  {
      double *out;
      plhs[0] =  mxCreateDoubleMatrix(1,1,mxREAL);
      out = mxGetPr(plhs[0]);
      out[0]=(double)sizeof(mwIndex);
      return;
    } else if (nlhs==0) {
#if MatlabVER>=73 /* NEW MATLAB 7.3 SIZE POINTERS */
      mexPrintf("SINCE73 : NO\n");
#else 
      mexPrintf("NO SINCE73 : YES\n");
#endif
    }


/*----------------------------------------------------------------- profile */
}  else if (!strcmp("profile",buf)) {
  double        *prind;
  int           N, i;
  mwIndex       *ir, *jc;
  /* check for proper number of arguments */
  if (nrhs != 2) mexErrMsgTxt ("Two inputs required."); 
  if (!mxIsSparse(prhs[1])) mexErrMsgTxt ("Matrix must be sparse.");

  N = mxGetN (prhs[1]); 
  ir = mxGetIr(prhs[1]); jc = mxGetJc(prhs[1]);  
  /* create matrix for return arguments */
  plhs[0] = mxCreateDoubleMatrix (N, 1, mxREAL); prind = mxGetPr (plhs[0]);

  /* calculate the profile prind of sparse matrix  */
  for (i=0; i<N; i++, prind++)   *prind = (double)(i-ir[(jc[i])]+1);
  
    

}/* k = sp_util('spind',k,ind); ------------------------------*/
else if (!strcmp("spind",buf))  {

  int        i, j, k, Nk, Nind, Nc, nzmax, cj;
  double     *ind, *pr, *pI, *ir_tmp, N, nj, kR;
  int        *rind;
  mwIndex    *ir, *jc, *ip;

  /* check if input is sparse */
  if (!mxIsSparse (prhs[1])) mexErrMsgTxt ("Matrix must be sparse.");
  if (nrhs != 3) mexErrMsgTxt ("three inputs are required."); 
  if (nlhs != 1) mexErrMsgTxt ("one output is required."); 
 
  pr = mxGetPr (prhs[1]); ir = mxGetIr (prhs[1]); jc = mxGetJc (prhs[1]);  

  Nk = mxGetN (prhs[1]); N = (double)(Nk);  Nc = jc[Nk];

  /* Get the starting position of indexation vector and its size */
  ind = mxGetPr (prhs[2]); Nind = mxGetM(prhs[2])*mxGetN(prhs[2]);


  if (Nind==Nk) {                   /* all rows are kept ---------*/

  /* This is very dirty : explicit modification of input */ 
    rind = (int*)mxCalloc (Nk,sizeof(int));
    for (i=0;i<Nind;i++) {j = (int)(ind[i])-1;
     if (j>Nk) mexErrMsgTxt ("Index larger than sparse matrix size"); 
     rind[j] = i;
    }
    ir_tmp = (double*)mxCalloc (jc[Nk],sizeof(double));
    for (i=0; i<Nk; i++) { kR=(double)(rind[i])*N+.1;
      for (j=jc[i]; j<jc[i+1]; j++) {
        ir_tmp[j] = kR+(double)(rind[ir[j]]);
      }
    }
    if (mxIsComplex(prhs[1])) {
#if MatlabVER >= 904
		mexErrMsgTxt("Error need complex reimplement");
#else
		pI = mxGetPi (prhs[1]);
#endif
		hpsortc (jc[Nk], ir_tmp, pr,pI); }/* sort ir_tmp, rearrange pr,pI */ 
    else {
      hpsort (jc[Nk], ir_tmp, pr); }/* sort ir_tmp, rearrange pr */ 

    k=0; nj=0; for (i=0;i<jc[Nk];i++){
      j = (int)(ir_tmp[i]/N);
      while (k!=j) {k++;jc[k]=i;nj=N*(double)(j);} /* column starts */
      ir[i] = (int)(ir_tmp[i]-nj);
    }
    while (k<Nind) {k++;jc[k]=jc[Nk];}

    plhs[0]=(mxArray *)(prhs[1]);
    /*    nzmax = jc[Nk];
    plhs[0] = mxCreateSparse (Nind,Nind,nzmax,mxREAL);
    ip = mxGetJc (plhs[0]);  for (i=0;i<=Nind;i++) {ip[i]=jc[i];}
    ip = mxGetIr (plhs[0]);  ind = mxGetPr (plhs[0]);
    for (i=0;i<nzmax;i++) {ip[i]=ir[i];ind[i]=pr[i];}*/
  }
  else {                 /* some rows are deleted ------------------ */

  rind = (int*)mxCalloc (Nk,sizeof(int));
  for (i=0;i<Nk;i++) rind[i] = Nk;
  for (i=0;i<Nind;i++) {j = (int)(ind[i])-1;
   if (j>Nk) mexErrMsgTxt ("Index larger than sparse matrix size"); 
   rind[j] = i;
  }

  ir_tmp = (double*)mxCalloc (jc[Nk],sizeof(double)); nj = N*(N+1);
  for (i=0; i<Nk; i++) { 
      if (rind[i]==Nk) {for (j=jc[i]; j<jc[i+1]; j++) {ir_tmp[j]=nj;}}
      else {kR=(double)(rind[i])*N+.1;
        for (j=jc[i]; j<jc[i+1]; j++) {
          cj=rind[ir[j]]; if (cj==Nk) ir_tmp[j] = nj;
          else ir_tmp[j]=kR+(double)(cj);}
      }
  }
  if (mxIsComplex(prhs[1])) {
#if MatlabVER >= 904
	  mexErrMsgTxt("Error need complex reimplement");
#else
	  pI = mxGetPi (prhs[1]);
#endif
      hpsortc (jc[Nk], ir_tmp, pr,pI); }/* sort ir_tmp, rearrange pr,pI */ 
  else {
      hpsort (jc[Nk], ir_tmp, pr); }/* sort ir_tmp, rearrange pr */ 

  k=0;jc[0]=0;nzmax=0;nj=0; for (i=0;i<jc[Nk];i++){
      j = (int)(ir_tmp[i]/N);
      if (j>=Nk){if (nzmax==0) {nzmax=i;} i = jc[Nk];}
      else {j = (int)(ir_tmp[i]/N);
       while (k!=j) {k++;jc[k]=i;nj=N*(double)(j);} /* column starts */
       ir[i] = (int)(ir_tmp[i]-nj);
      }
      /*printf("(%i,%i,%i,%5.1f)",i,j,k,pr[i]);*/
  }
  if (nzmax==0) nzmax = jc[Nk];
  while (k<Nind) {k++;jc[k]=nzmax;}

  if (mxIsComplex(prhs[1])) {
    plhs[0] = mxCreateSparse (Nind,Nind,nzmax,mxCOMPLEX);
    ip = mxGetJc (plhs[0]);   for (i=0;i<=Nind;i++) {ip[i]=jc[i];}
#if MatlabVER >= 904
	mexErrMsgTxt("Error need complex reimplement");
#else
	ip=mxGetIr (plhs[0]); ind = mxGetPr (plhs[0]); ir_tmp = mxGetPi (plhs[0]);
#endif
	for (i=0;i<nzmax;i++) {ip[i]=ir[i];ind[i]=pr[i];ir_tmp[i]=pI[i];}}
  else {
    plhs[0] = mxCreateSparse (Nind,Nind,nzmax,mxREAL);
    ip = mxGetJc (plhs[0]);   for (i=0;i<=Nind;i++) {ip[i]=jc[i];}
    ip = mxGetIr (plhs[0]);  ind = mxGetPr (plhs[0]);
    for (i=0;i<nzmax;i++) {ip[i]=ir[i];ind[i]=pr[i];}}

  }

}/*--------------------------------------------------------- xkx  */
  else if (!strcmp("xkx",buf))  {

    double *offset;
    int    type,ym,yn;

  mexWarnMsgTxt("Use OF_MK call instead of SP_UTIL");

  /* check for proper number of arguments */
  if (nrhs < 3) {
    mexErrMsgTxt ("Three inputs required.");
  } 
  else if(nlhs!=1) {
    mexErrMsgTxt ("One output required.");
  }
  
  /* check the inputs sizes */
  if (mxGetM (prhs[1]) != 3 || mxGetN (prhs[1]) != 3) {
   mexErrMsgTxt ("First argument X must be 3 by 3 in c = xkx(X,K)");
  }

  if (mxIsSparse(prhs[2])) mexErrMsgTxt ("k must be full");
  /* create matrix for the return argument */
  ym=mxGetM(prhs[2]);yn=mxGetN(prhs[2]);
  if (mxIsComplex(prhs[2])) { 
    plhs[0] = mxCreateDoubleMatrix (ym, yn, mxCOMPLEX);
  } else  plhs[0] = mxCreateDoubleMatrix (ym, yn, mxREAL);
  if (nrhs>3) {
   if   (mxGetM(prhs[3])==0) offset=NULL; 
   else if (mxGetM(prhs[3])!=3) { mexErrMsgTxt ("offset must be 3xN");}
   else offset=mxGetPr(prhs[3]); 
  } else offset=NULL;

  if (nrhs==5) type=(int)*mxGetPr(prhs[4]); else type=0;

  if (mxIsComplex (prhs[1])) mexErrMsgTxt ("Matrix of basis change must not be complex.");
  if (mxIsComplex (prhs[2])) { 
#if MatlabVER >= 904
	  mexErrMsgTxt("Error need complex reimplement");
#else
	  x_k_x(mxGetPr(prhs[1]),mxGetPi(prhs[2]),mxGetPi(plhs[0]),offset,type,
     ym,yn);
#endif
  }
  /* call the subroutine xkx */
  x_k_x (mxGetPr (prhs[1]),mxGetPr(prhs[2]),mxGetPr(plhs[0]),offset,type,
   ym,yn);

}/* ----------------------------------------------------------------
out1 = sp_util('mind',ki, ke,length(Up.DOF),[mind]);
out1 = sp_util('mindsym',ki, ke,length(Up.DOF),[mind]); 
out1 = sp_util('mind','xxx',[],length(Up.DOF),[mind]);
out1 = sp_util('mindsym','xxx',[],length(Up.DOF),[mind]);
nees fid=fopen('xxx','wb','l');
*/
else if (!strcmp ("mind",buf) || !strcmp ("mindsym",buf))  {

  double  *mind, *ki, *ke, K,  r3, nj, cur, nex, *pr, *r1, *r2; 
  int     i, i2, j, k, l, NE,  n, nzmax, nzc,  ooc, i1;
  
  mwIndex  *jci,*ir, *jc, *jcj;

  char    *buf1;

  if (nrhs<4) mexErrMsgTxt ("Bad number of inputs.");
  if (nlhs>1) mexErrMsgTxt ("Bad number of outputs.");
  if (mxIsComplex(prhs[2])) mexErrMsgTxt ("ke must be real");
  if (mxIsSparse(prhs[1])) mexErrMsgTxt ("ki must be full");
  if (mxIsSparse(prhs[2])) mexErrMsgTxt ("ke must be full");
 
  /* deal with how ki ke is given */
  if (mxIsInt32(prhs[2])) { /* as a single matrix */
   int *pi;
   ki = mxGetPr (prhs[1]);  ke = ki+mxGetM(prhs[1]);
   pi=(int*)mxGetData(prhs[2]);NE=(int)pi[0];  ooc=1;
  } else if (mxIsChar(prhs[1])) { /* kie in a file */

   /*    int       fb, n[1]; */
   int n[1];
   FILE *fb;

    buf1 = mxArrayToString( prhs[1] ); fb=fopen(buf1,"rb"); mxFree(buf1);
    fread(n,sizeof(int),1,fb);NE=n[0];
    ki=(double*)mxMalloc(2*n[0]*sizeof(double));
    fread(ki,sizeof(double),2*n[0],fb);
    fclose(fb);
    ke=ki+NE;  ooc=0;

  } else { /* Ki, ke given as arguments */ 
    ke = mxGetPr (prhs[2]);  ki = mxGetPr (prhs[1]);   
    NE = mxGetM (prhs[1])*mxGetN (prhs[1]); 
    ooc=2;
    if (mxGetM (prhs[3]) == 2 || mxGetN (prhs[3]) == 2) {
     pr = mxGetPr (prhs[3]); n = (int)(pr[1]+.1); 
     if (n>NE){ mexErrMsgTxt ("NW must be <= length(ki)");}
     NE=n;
    } 
    if (mxGetM (prhs[1])!=mxGetM (prhs[2])) {
      mexErrMsgTxt ("data and index must have the same size");
    }
  }  /* format of input arguments - - - - - - - - - - - - - - - - - - -*/
  if (NE==0)  { /* safe return if empty matrix */
     K = *mxGetPr (prhs[3]); n = (mwSize)(K+.1);
     plhs[0] = mxCreateSparse (n,n,0,mxREAL);
     return;
  } 
  if (nrhs==5)  {/* mind is given - - - - - - - - - - - - - - - - - - - -*/

    int Nmi;

    r1 = mxGetPr (prhs[4]);     Nmi = mxGetM (prhs[4]); 
    mind = (double*)mxCalloc (3*Nmi,sizeof(double));
    for (i=0; i<3*Nmi; i++)  mind[i]=r1[i];
    r1 = mind+Nmi; r2 = mind+2*Nmi; 
    hpsortc (Nmi,mind,r1,r2); /* sort mind */ 

    i=0; 
    for (j=0;j<Nmi;j++) { /* loop on elements with non-zero coefficients */
     /* printf("%i-%i ",i,(int)(mind[j]-.99)); */
     /* indices and coefficient */
     if (mind[j]<=mind[j+Nmi]) {
      i1=(int)(mind[j]-.99);  i2=(int)(mind[j+Nmi]-.99); K = mind[j+2*Nmi];
      if (i1<0)  mexErrMsgTxt ("indices must be positive");
      if (i2>NE)  mexErrMsgTxt ("mind > size(ke)");
      if (K==1.0) { /* printf("km%4i %4i %4i %e\n",i,i1,i2,K); */
        if (i>i2) {} /* empty matrix */  
        else if (i!=i1) { /* translate if zero coef before */
          for (k=i1;k<=i2;k++) { ke[i] = ke[k];ki[i] = ki[k];i++; }
        } else if (i<i2+1) {i=i2+1;}
      } else if (K!=0.0) {
          for (k=i1;k<=i2;k++) {ke[i] = K*ke[k];ki[i] = ki[k];i++; }
      }
    }}; 
    if (i<=NE) {NE = i;}
    mxFree(mind);
  } /* end of mind multiplication  - - - - - - - - - - - - - - - - - - - */

  K = *mxGetPr (prhs[3]); n = (int)(K+.1);

  /* sort ind and rearrange data in the same way */ 
  hpsort (NE, ki, ke);
  if (ki[NE-1]>((double)n)*((double)n)) {
   mexErrMsgTxt ("index exceeds matrix dimension (ki>N^2)");
  }

  if (NE==0)  { /* safe return if empty matrix, available for mind and mindsym */
     plhs[0] = mxCreateSparse (n,n,0,mxREAL);
     return;
  } 

  /* - - - -  - - - - - - - - - - - - - - - - - - - - - - - - mind - */
  if (!strcmp ("mind",buf)) {

   /* eliminate repeated values  */
   jci = (mwIndex*)mxCalloc(n+1, sizeof(mwIndex));  for(k=0;k<=n;k++){jci[k]=0;}
   nzmax=0; cur = ki[nzmax]-1;j=(int)(cur/K);nj=K*(double)(j)-0.1;
   ki[nzmax]=cur-nj;nzc = 0;

   for(k=1;k<NE;k++) {
     nex = ki[k]-1;
     if(nex==cur) { ke[nzmax]+=ke[k]; }  /* sum repeated */
     else {
       if  (ke[nzmax]!=0.0) {/* (ke[nzmax]!=0.0)accepting non-0 new values */
         jci[j+1]++;nzmax++;ke[nzmax]=ke[k];}
       else {     /* removing zero new value */
          ke[nzmax]=ke[k];/*jci[j+1]--;printf("bad");*/
       }
       cur = nex; r3=cur-nj; ki[nzmax]=r3;
       if(r3>=K){ /* new column */
          j=(int)(cur/K);nj=K*(double)(j)-0.1;r3=cur-nj; ki[nzmax]=r3;
       }
     }
   } /* loop on k */
   if (ke[nzmax]!=0.0) { jci[j+1]++; } else {nzmax--;}

   nzmax++;
   /* printf("nzmax %i",nzmax); */
   for(k=1;k<=n;k++){jci[k]+=jci[k-1];} /* cumsum for column starts */
   /* create sparse output matrix */
   plhs[0] = mxCreateSparse ((mwSize)n,(mwSize)n,(mwSize)nzmax,mxREAL);
   pr = mxGetPr (plhs[0]); ir = mxGetIr (plhs[0]); jc = mxGetJc (plhs[0]);
   memcpy (pr, ke, nzmax*sizeof(double)); 
   memcpy (jc, jci, (n+1)*sizeof(mwIndex));
   mxFree(jci);
   for(k=0;k<nzmax;k++) { ir[k]=(mwIndex)(ki[k]); }

  /* - - - -  - - - - - - - - - - - - - - - - - - - - - - - - mindsym - */
  } else if (!strcmp ("mindsym",buf)) {

   /* eliminate repeated values */
   K = (double)(n); 
   jci = (mwIndex*)mxCalloc(n+1, sizeof(mwIndex));  for(k=0;k<=n;k++){jci[k]=0;}
   jcj = (mwIndex*)mxCalloc(n+1, sizeof(mwIndex));  for(k=0;k<=n;k++){jcj[k]=0;}
   nzmax=0; cur = ki[nzmax]-1;j=(int)(cur/K);nj=K*(double)(j)-0.1;
   r3=cur-nj; ki[nzmax]=r3; nzc = 0;


   for(k=1;k<NE;k++) {
     nex = ki[k]-1;
     if (nex==cur) { ke[nzmax]+=ke[k]; }  /* sum repeated */
     else {
       if  (ke[nzmax]!=0.0) {/* (ke[i]!=0.0)accepting non-zero new values */
         jci[j+1]++; jcj[j+1]++; 
	 if (r3>2e9||(int)(r3)!=j) {jcj[(int)(r3)+1]++;} else {nzc++;}
	 /*	 if ((int)(r3)!=j) {jcj[(int)(r3)+1]++;} else {nzc++;}*/
         nzmax++;ke[nzmax]=ke[k];}
       else {     /* removing zero new value */
          ke[nzmax]=ke[k];
       }
       cur = nex; r3=cur-nj; ki[nzmax]=r3;
       if(r3>=K){ /* new column */
         j=(int)((cur+.01)/K);nj=K*(double)(j)-0.1;r3=cur-nj; 
	 /* if (r3>K) {
	  mexPrintf("%.1f j=%i=%.1f r3=%.1f\n",cur,j,nj,r3);
          mexErrMsgTxt("bad");
	  }*/
         ki[nzmax]=r3;
       }
     }
   } /* loop on k */

   if (ke[nzmax]!=0.0) {
    jci[j+1]++;jcj[j+1]++; 
    if(r3>2e9||(int)(r3)!=j) {jcj[(int)(r3)+1]++;} else {nzc++;}
   } else {nzmax--;}

   for(k=1;k<=n;k++){jci[k]+=jci[k-1];jcj[k]+=jcj[k-1]; } /* cumsum c. starts */

   /* create sparse output matrix */
   plhs[0] = mxCreateSparse (n,n,(int)jcj[n],mxREAL);
   pr = mxGetPr (plhs[0]); ir = mxGetIr (plhs[0]); jc = mxGetJc (plhs[0]);
   memcpy (jc, jcj, (n+1)*sizeof(mwIndex));

   for(j=n-1;j>=0;j--) {
     /* mexPrintf("\n j %i %i,",j+1,jci[j+1]-1);*/
     for(k=(int)jci[j+1]-1;k>=(int)jci[j];k--) {
         i = j+1; l = (int)jcj[i]-1;pr[l] = ke[k];ir[l] = (mwIndex)(ki[k]);jcj[i]--;
	 /* mexPrintf(" %i:%i:%i (%.1f) ",i,l,ir[l],ki[k]); */
	 i = (int)ir[l]+1;
         if (i!=j+1){l=(int)jcj[i]-1; pr[l] = ke[k]; ir[l] = (mwIndex)j;jcj[i]--; }
     }
   } 
   mxFree(jci); mxFree(jcj);

  } /* if mindsym */
  if      (ooc==0) { ki[0]=-1;  /* mxFree(ki); */} 
  else if (ooc==1) { ki[0]=-1;  /* mxFree(ki);  SetPrToNull(prhs[1]); */}
  else             {ki[0]=-1; 
  /* mxFree(ki); SetPrToNull(prhs[1]);
     mxFree(ke); SetPrToNull(prhs[2]);*/
 }



}/* ----------------------------------------------------------------
sparse = sp_util('sp2st',k); 
sparse=sp_util('sp2st',double([N1 N2 Nzmax]),ir,jc,pr)
 */
else if (!strcmp ("sp2st",buf))  {

  double        *pr, *v;
  mxArray       *field_value; 
  int           n, j1; /* , dims[2] = {1,1};*/
  mwSize        dims[2], nz;

  const char    *field_names[] = {"pr", "ir", "jc","nzmax"};
  mwIndex       *ir, *jc; 

  if (nrhs==2) {
  if (!mxIsSparse(prhs[1])) mexErrMsgTxt ("Matrix must be sparse."); 
  pr = mxGetPr (prhs[1]); ir = mxGetIr (prhs[1]); jc = mxGetJc (prhs[1]);

  dims[0]=1; dims[1]=1;

  plhs[0] = mxCreateStructArray(2, dims, 4, field_names);  
  n = mxGetNzmax (prhs[1]);
  field_value = mxCreateDoubleMatrix(1,n,mxREAL);
  v = mxGetPr (field_value); 
  mxSetField(plhs[0], 0, "pr", field_value); 
  for (j1=0;j1<n;j1++) { v[j1] = pr[j1]; }

  field_value = mxCreateDoubleMatrix(1,n,mxREAL);
  v = mxGetPr (field_value); 
  mxSetField(plhs[0], 0, "ir", field_value); 
  for (j1=0;j1<n;j1++) { v[j1] = (double)(ir[j1]);}

  n = mxGetN (prhs[1]);
  field_value = mxCreateDoubleMatrix(1,n+1,mxREAL);
  v = mxGetPr (field_value); 
  mxSetField(plhs[0], 0, "jc", field_value); 
  for (j1=0;j1<=n;j1++) { v[j1] = (double)(jc[j1]); }

  field_value = mxCreateDoubleMatrix(1,1,mxREAL);
  v = mxGetPr(field_value);  mxSetField(plhs[0], 0, "nzmax", field_value); 
  v[0]=(double)mxGetNzmax(prhs[1]); 
  } else if (nrhs>=5) { /* create a sparse matrix */

   if (!mxIsDouble(prhs[1])) mexErrMsgTxt ("Size must be given as double");
   pr=mxGetPr(prhs[1]); 
   if (mxIsDouble(prhs[2])) {
    nz=(mwSize) pr[2];
	if (mxIsComplex(prhs[4])) plhs[0] = mxCreateSparse((mwSize) pr[0], (mwSize) pr[1],nz,mxCOMPLEX);  
	else plhs[0] = mxCreateSparse((mwSize)pr[0], (mwSize)pr[1], nz, mxREAL);

    pr=mxGetPr(prhs[2]);ir=mxGetIr(plhs[0]);  
    for (j1=0;j1<nz;j1++) {ir[j1]=(mwIndex)pr[j1];} 

    jc=mxGetJc(plhs[0]); pr=mxGetPr(prhs[3]);
    for (j1=0;j1<mxGetNumberOfElements(prhs[3]);j1++) {jc[j1] = (mwIndex)(pr[j1]);}
    pr=mxGetPr(plhs[0]);
    memcpy (pr,mxGetPr(prhs[4]),mxGetNumberOfElements(prhs[4])*sizeof(double));  
    
  }  else if (mxGetElementSize(prhs[2])==sizeof(mwIndex)) {
	if (mxIsComplex(prhs[4])) plhs[0] = mxCreateSparse((mwSize) pr[0], (mwSize) pr[1],(mwSize) pr[2],mxCOMPLEX);
	else plhs[0] = mxCreateSparse((mwSize)pr[0], (mwSize)pr[1], (mwSize)pr[2], mxREAL);

    ir=mxGetIr(plhs[0]);   jc=mxGetJc(plhs[0]);  pr=mxGetPr(plhs[0]);
    
    /*mxFree(ir);mxFree(jc);
    mxSetIr(plhs[0],mxGetData(prhs[2]));mxSetPr(prhs[2],NULL);
    mxSetJc(plhs[0],mxGetData(prhs[3]));mxSetPr(prhs[2],NULL); */
    memcpy (ir,mxGetData(prhs[2]),mxGetNumberOfElements(prhs[2])*sizeof(mwIndex));
    memcpy (jc,mxGetData(prhs[3]),mxGetNumberOfElements(prhs[3])*sizeof(mwIndex));
    memcpy (pr,mxGetData(prhs[4]),mxGetNumberOfElements(prhs[4])*sizeof(double));
    if (mxIsComplex(prhs[4])) {
#if MatlabVER >= 904
		mexErrMsgTxt("Error need complex reimplement");
#else
		pr=mxGetPi(plhs[0]);
      memcpy (pr,mxGetPi(prhs[4]),mxGetNumberOfElements(prhs[4])*sizeof(double));
#endif
	} 
   } else { mexErrMsgTxt ("Pointer class mismatch");}
  
   /* mxSetNzmax(plhs[0],mxGetNumberOfElements(prhs[4]));*/
   
  }

}/* ----------------------------------------------------------------
 sp_util('insertinkie',kelem,iddl,kie,ik,N); 
kelem=[1 2 3;4 5 6;7 8 9]';
iddl=int32([5 6 9]);
kie=zeros(20,2);ik=int32(9);
sp_util('insertinkie',sparse(kelem),iddl,kie,ik,100);
kie
ik

*/
else if (!strcmp ("insertinkie",buf))  {

  double        *kelem, *kie, N;
  int           NE, Nkie, *ik, j1, j2, *iddl, flag; 
  mwIndex       *ir, *jc;

  /* handle with sparse or full input */ 
  if (mxIsSparse (prhs[1])) {
    flag=1;
    kelem = mxGetPr (prhs[1]); ir = mxGetIr (prhs[1]); jc = mxGetJc (prhs[1]);
  }
  else {
    flag=0;
    kelem = mxGetPr (prhs[1]);
  }

  if (nrhs != 6) mexErrMsgTxt ("Not enough inputs."); 
#ifdef MatlabVER
  if (nlhs != 0) mexErrMsgTxt ("No output is required."); 
#else
  if (nlhs != 2) mexErrMsgTxt ("Two output is required."); 
  plhs[0] = prhs[3];plhs[1] = prhs[4];
#endif
  /* matrix entry kelem */ 
  NE = mxGetN (prhs[1]); 
  if (NE!=mxGetM (prhs[1])) mexErrMsgTxt ("Matrix must be square.");

  iddl = (int*)mxGetData (prhs[2]);  
  kie = mxGetPr (prhs[3]); Nkie = mxGetM (prhs[3]);
  ik = (int*)mxGetData (prhs[4]); N = *mxGetPr (prhs[5]); 


  if (flag) { /* is Sparse */
    for (j1=0;j1<NE;j1++) { /* column jc */

      for (j2=jc[j1];j2<jc[j1+1];j2++) { /* row ir */
        if (iddl[j1]>=0 && iddl[ir[j2]]>=0 && kelem[j2]) {
          kie[*ik]  = (double)(1+iddl[ir[j2]])+N*(double)(iddl[j1]); /* ind */
          kie[*ik+Nkie]=    kelem[j2]; /* value */
          ik[0]++;
          if (ik[0]>Nkie) { 
            Nkie=(int)((double)Nkie*1.25)+10;
			mexErrMsgTxt("Realloc requires output");
            /* FullRealloc(prhs[3],Nkie,2); kie = mxGetPr (prhs[3]);*/
          }
	}
      }
    }
  }
  else {
    for (j1=0;j1<NE;j1++) {/* column */
      for (j2=0;j2<NE;j2++) { /* row */
        if (iddl[j1]>=0 && iddl[j2]>=0 && kelem[j2+NE*j1]) {
          if (ik[0]>Nkie) { 
            Nkie=(int)((double)Nkie*1.25)+NE^2;
            /*FullRealloc(prhs[3],Nkie,2); kie = mxGetPr (prhs[3]);*/
	        mexPrintf("At %i\n",Nkie);mexErrMsgTxt("Kie buffer overflow");
          }
          kie[*ik]=(double)(1+iddl[j2])+N*(double)(iddl[j1]);/* indice */
          kie[*ik+Nkie]=kelem[j2+NE*j1]; /* value */
          ik[0]++;
        }
      }
    }
  }


}/* ---------------------------------------------------------------- */
else if (!strcmp ("basis",buf))  {

  double     *x, *y, *z, *out2, eps=1.e-16;

  
  if (nrhs == 3) { /* sp_util('basis',x,y) */

   if (mxIsSparse(prhs[1])) mexErrMsgTxt ("X must be full");
   if (mxIsSparse(prhs[2])) mexErrMsgTxt ("Y must be full");
   if ((mxGetM (prhs[1])!=3 || mxGetM (prhs[1])!=3) &&
       (mxGetN (prhs[1])!=3 || mxGetN (prhs[1])!=3)  ) 
       mexErrMsgTxt ("Input must be 3 by 1 vectors.");
   x = mxGetPr (prhs[1]); y = mxGetPr (prhs[2]); 
   plhs[0] = mxCreateDoubleMatrix(3,3,mxREAL);  
   z = mxGetPr (plhs[0]);
   basis(x,y,z);

  } else if (nrhs == 2) { 
       /*
       Node=[0 0 0; 0 1 0; 1 1 0; 1 0 0]'
       [b,r,t]=sp_util('basis',Node)      
     */
    if ((mxGetN (prhs[1])==4) || (mxGetN (prhs[1])==8)) { /* quad4 basis */
      double     *nodeE, *x1;
      int      Elt[4]={1,2,3,4}; 

      nodeE = mxGetPr (prhs[1]);              
      plhs[1] = mxCreateDoubleMatrix(4,3,mxREAL);  z=mxGetPr(plhs[1]);
      plhs[2] = mxCreateDoubleMatrix(1,3,mxREAL);  out2=mxGetPr(plhs[2]);
      /* plhs[0] = mxCreateDoubleMatrix(9,1,mxREAL);  x1=mxGetPr(plhs[0]); */
      plhs[0] = mxCreateDoubleMatrix(3,3,mxREAL);  x1=mxGetPr(plhs[0]);
      basisQ4(nodeE, 1, Elt, x1,z,out2,mxGetN (prhs[1]));
    } else if ((mxGetN (prhs[1])==3) || (mxGetN (prhs[1])==6)) {/* tria3 basis */
      double     *nodeE, *z, x[3], y[3];

      nodeE = mxGetPr (prhs[1]);              
      plhs[0] = mxCreateDoubleMatrix(3,3,mxREAL);  z=mxGetPr(plhs[0]);
      x[0]=nodeE[3]-nodeE[0];x[1]=nodeE[4]-nodeE[1];x[2]=nodeE[5]-nodeE[2];
      y[0]=nodeE[6]-nodeE[0];y[1]=nodeE[7]-nodeE[1];y[2]=nodeE[8]-nodeE[2];
      basis(x,y,z);
    } else        mexErrMsgTxt ("Not a know number of nodes");

  } /* if (nrhs == 3)*/
}/* ----------------------------------------------------------------
 % standard call
 Node=[0 0 0; 0 1 0; 1 1 0; 1 0 0;0 0 1; 0 1 1; 1 1 1; 1 0 1]';
 Elt=int32([1 2 3 4;5 6 7 8]');
 b=sp_util('basiselt',Node,Elt);
 reshape(b,3,6)

 */
else if (!strcmp ("basiselt",buf))  {

  double     *Node, *b, eps=1.e-16;
  int      *Elt;
  int        Nelt, nPerElt;

  if (nrhs != 3)  
    mexErrMsgTxt ("3 input arguments needed :  b=sp_util('basiselt',Node,Elt)");
  if (mxIsSparse(prhs[1])) mexErrMsgTxt ("Node must be full");
  if (mxIsSparse(prhs[2])) mexErrMsgTxt ("Elt must be full");
  if (mxGetM (prhs[1])!=3) mexErrMsgTxt ("Node must have 3 rows");
  if (!mxIsInt32(prhs[2])) mexErrMsgTxt ("Elt connectivity must be int32");

  Nelt=(int)mxGetN(prhs[2]);  Node = mxGetPr (prhs[1]);  Elt = (int*)mxGetData(prhs[2]);
  plhs[0] = mxCreateDoubleMatrix(9, Nelt, mxREAL);  b=mxGetPr(plhs[0]);
  nPerElt=mxGetM (prhs[2]);
  if (nPerElt==3) {basisT3(Node, Nelt, Elt, b,3);}
  else if (nPerElt==4) { basisQ4(Node, Nelt, Elt, b, NULL, NULL,4);}
  else if (nPerElt==8) { basisQ4(Node, Nelt, Elt, b, NULL, NULL,8);}
  else if (nPerElt==6) { basisT3(Node, Nelt, Elt, b,6);}
  else {
   mexErrMsgTxt ("sp_util('basiselt',Node,Elt) : unrecognized element");
  }
 

/*---------------------------------------------------------  ERROR  */



#ifdef SDT_ADD
#include "../../sdtdev/mex50/sp_util_post.c"
#endif


} else if (!strcmp("cvs",buf)) {

    mxArray *st;
    mxArray *rhs[1], *lhs[1];
    mwSize  *dims;
    st=mxCreateString("$Revision: 1.109 $  $Date: 2023/02/21 10:29:41 $");

#ifdef SDT_ADD

    rhs[0] = mxCreateString("post_cvs");

    dims= (mwSize*)mxCalloc(2,sizeof(mwSize));

    dims[0]=4; dims[1]=2;  plhs[0]=mxCreateCellArray(2,dims);

    mxSetCell(plhs[0],0,      mxCreateString("pre"));
    mxSetCell(plhs[0],1,      mxCreateString("sp_util"));
    mxSetCell(plhs[0],2,      mxCreateString("post"));
    mxSetCell(plhs[0],3,      mxCreateString("subs"));

    mxSetCell(plhs[0],4,pre_cvs ());         /* pre */
    mxSetCell(plhs[0],5,st);                 /* sp_util */
    mexCallMATLAB(1, lhs, 1, rhs, "sp_util");
    mxSetCell(plhs[0],6,lhs[0]);             /* post */
    mxSetCell(plhs[0],7,subs_cvs_sputil());         /* subs */
    mxFree(dims); 

    mxDestroyArray(rhs[0]);
#else
    plhs[0]=st; 
#endif

} else if (!strcmp("subscvs",buf)) {

  plhs[0]=subs_cvs_sputil ();

/*----------------------------------------------------------- epsl */
}  else if (!strcmp("epsl",buf))  {
  double *a;

  if (nrhs==1) {
   plhs[0] = mxCreateDoubleMatrix (1, 1, mxREAL);
   a = mxGetPr (plhs[0]); *a=EPSL; 
  }
  else if (nrhs==2) { a = mxGetPr (prhs[1]); EPSL = *a;}

/*----------------------------------------------------------- epsl */
}  else if (!strcmp("diag",buf))  {
  double *a;

  if (nrhs==1) {
   plhs[0] = mxCreateDoubleMatrix (1, 1, mxREAL);
   a = mxGetPr (plhs[0]); *a=OpenFEMDIAG; 
  }
  else if (nrhs==2) { a = mxGetPr (prhs[1]); OpenFEMDIAG = *a;}

/*----------------------------------------------------------- ismex */
} else  if (!strcmp("issdt",buf))  {
  double *a;

  plhs[0] = mxCreateDoubleMatrix (1, 1, mxREAL);
  a = mxGetPr (plhs[0]);
#ifdef SDT_ADD
  *a=1.0;
#else
  *a=0.0;
#endif


/* addcol ------------------------------------------------------------
   add end column to sparse matrix k1 of size nk1 x nk1
   matrix k2 is size nk1+ncol x nk1+ncol
   k2=sp_util('addcol',k1,ncol)

                 --------------
                 -        -   -
                 -  k1    -col-
             k2= -        -   -
                 --------------
                 -  col   -   -
                 --------------

   28/02/2006 ------------------------------------------------------*/
    } else if (!strcmp("addcol",buf))  {

      mxArray  *out;
      int      nNod, nnz, ii, jj, ncol;
      mwIndex  *icol, *prow, *nprow, *nicol, aux;
      double   *val;


      nNod = mxGetM( prhs[1] );
      prow = mxGetJc( prhs[1] );
      icol = mxGetIr( prhs[1] );
      ncol = (int)mxGetScalar( prhs[2] );

 
      aux=prow[nNod]-1;
      aux += ncol*(ncol+2*nNod);

      nprow = (mwIndex*)mxCalloc( nNod+ncol+1, sizeof( mwIndex ) );
      nicol = (mwIndex*)mxCalloc( aux, sizeof( mwIndex ) );

      nprow[0]=prow[0];
      nicol[0]=icol[0];


      for (ii=1;ii<nNod+1;ii++) {
	nprow[ii]=prow[ii]+ncol*ii;
	for (jj=0;jj<prow[ii]-prow[ii-1];jj++) {
	  nicol[nprow[ii-1]+jj]=icol[prow[ii-1]+jj];
	}
	for (jj=0;jj<ncol;jj++) {
	  nicol[nprow[ii]-ncol+jj]=nNod+jj;
	}
      }

      for (ii=0;ii<ncol;ii++) {
	nprow[nNod+1+ii]=nprow[nNod+ii]+nNod+ncol;
	for (jj=0;jj<nNod+ncol;jj++) {
	  nicol[jj+nprow[nNod+ii]]=jj;
	}
      }

      nnz=nprow[nNod+ncol];
      out = mxCreateSparse( nNod+ncol, nNod+ncol, nnz, mxREAL );
      val = mxGetPr( out);

      for (ii = 0; ii < nnz; ii++) {
	val[ii] = 1.0;
      }
      mxSetIr( out, nicol );
      mxSetJc( out, nprow );
      plhs[0] = out;
      

} else {
    mexErrMsgTxt ("Unknown command.");
}
  
}

