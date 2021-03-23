#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "mex.h" 
/*#include "../../include/mkl_blas.h"*/

extern int omp_get_max_threads();
/* PARDISO prototype. */
#if defined(_WIN32) || defined(_WIN64)
#define pardiso_ PARDISO
#define pardisoinit_ PARDISOINIT
#else
#define PARDISO pardiso_
#define PARDISOINIT pardisoinit_
#define intP mwSize
#endif
extern int PARDISO
	(void *, int *, int *, int *, int *, int *,
	double *, int *, int *, int *, int *, int *,
	int *, double *, double *, int *);

extern int PARDISOINIT
        (void *, int *, int *);
/*
ofact('pardiso');ofact(speye(10),rand(10,1))
*/
/*-------------------------------------------------------------*/
void sd_mkl_get_sparse(mxArray *field, void* *irp, int *dims) {

  mxArray *f2;
  double *r1;
  
  dims[2]=0; /* default not transposed */
  if (mxIsClass(field,"v_handle")) {
    field= mxGetField(field,0,"file");
  } 
  if (mxIsSparse(field)) {
    /*if (sizeof(mwIndex)==8) mexErrMsgTxt("ir must be int32");*/
   irp[2]= mxGetPr(field);irp[1]= mxGetJc(field); irp[0]= mxGetIr(field);
   dims[0]=mxGetM(field); dims[1]=mxGetN(field);
  } else if (mxIsStruct(field)) {
   f2=mxGetField(field,0,"pr");
   if (f2==NULL)   {irp[2]=NULL; irp[0]=NULL; irp[1]=NULL;}
   else {
     irp[2]=mxGetPr(f2);
     irp[0]=mxGetData(mxGetField(field,0,"ir"));
     f2=mxGetField(field,0,"trans");
     if (f2!=NULL&&*mxGetPr(f2)!=0) {dims[2]=1;}
     f2=mxGetField(field,0,"dims"); 
     field=mxGetField(field,0,"jc");   irp[1]=mxGetData(field);
     if (f2!=NULL) {
      r1=mxGetPr(f2);dims[0]=(int)r1[0];dims[1]=(int)r1[1];  
     } else {dims[0]=mxGetM(field)-1;dims[1]=dims[0];} /* default square */
     /* mexPrintf("%p %p %p",irp[2],irp[0],irp[1]); */
   }
  } else  { /* Full matrix */
    irp[2]=mxGetPr(field); irp[0]=NULL; irp[1]=NULL;
    dims[0]=mxGetM(field); dims[1]=mxGetN(field);
  }
}
        

/*static void ClearAllFactors(void) {
   mexPrintf("Cleanup of PARDISO factors not properly implemented\n");
}*/

/* Internal solver memory pointer pt, */
/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
/* or void *pt[64] should be OK on both architectures */
/* static void *PardisoPt[64]; */
/*--------------------------------------------------------------------*/
typedef struct _FactorInfo FactorInfo ;
struct _FactorInfo {
  void         *PardisoPt[64];
  int          nRows;
  int          *rows_1;
  int          *columns_1;
  double       *pr;
} ;
/*--------------------------------------------------------------------*/
void ClearMtx(FactorInfo *FI, int cF) {
  mexPrintf("Clear factor\n");
}
/*--------------------------------------------------------------------*/
static FactorInfo     FI[20];
static int            cF=-1, indF[20];

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

void ClearAllFactors(void) 
{

  int             k, phase, idum, nRhs=0,maxfct=1,mnum=1,mtype,msglvl=0,error=0;
  int iparm[64];
  double          ddum;

  phase = -1;

  mexPrintf("Cleanup of PARDISO factors\n");
  k=0;
  while (!(FI[k].rows_1 == NULL))
    {
      
      FI[k].nRows=0;
      pardiso_ (FI[k].PardisoPt, &maxfct, &mnum, &mtype, &phase,
		&(FI[k].nRows), FI[k].pr, FI[k].rows_1, FI[k].columns_1, &idum, &nRhs,
		iparm, &msglvl, &ddum, &ddum, &error);
      mxFree(FI[k].rows_1); 
      mxFree(FI[k].columns_1); 
      mxFree(FI[k].pr);
      mexPrintf("factor %i cleaned\n",k);
      k++; 
    }
  return;
}



void PardisoError(int error) {
 if (error==-2)   mexPrintf("error = %i (Out of memory)\n",error);
}

void ClearFactor(int cF, int maxfct, int mnum, int mtype, int* param,
   int msglvl) {

  int             phase, idum, nRhs=0;
  double          ddum;

 phase = -1; /* Release internal memory. */
 cF=0;
 if (msglvl) mexPrintf("Clearing factor %i",cF);

 if ((FI[cF].rows_1 == NULL) && (FI[cF].columns_1 == NULL)  && (FI[cF].pr == NULL)) {
   return;  
 }
 /*
 if (FI[cF].rows_1 == NULL) mexErrMsgTxt("Problem with rows_1 pointers");
 if (FI[cF].columns_1 == NULL) mexErrMsgTxt("Problem with columns_1 pointers");
 if (FI[cF].pr == NULL) mexErrMsgTxt("Problem with pr pointers");

*/
 PARDISO (FI[cF].PardisoPt, &param[100], &param[101], &mtype, &phase,
    &(FI[cF].nRows), FI[cF].pr, FI[cF].rows_1, FI[cF].columns_1, &idum, &nRhs,
    param, &param[102], &ddum, &ddum, &param[103]);

 FI[cF].nRows=0;
 mxFree(FI[cF].rows_1); mxFree(FI[cF].columns_1); mxFree(FI[cF].pr); 

 FI[cF].rows_1=NULL;

}


/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
void mexFunction (int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]) {

char            *buf ;
 int             buflen;
double          *a;
int             j1, j2,j3, nCols, nZmax, nRhs, *param, dims[3];
int             *rows, *columns;
double          *rhs, *q, *qi;
char            statIn[] = "determinant";
int             mtype = -2;   /* Real symmetric matrix */
int             maxfct, mnum, phase, error, msglvl,nonsym;
int             i;    /* Auxiliary variables. */
double          ddum; /* Double dummy */
int             idum; /* Integer dummy. */
void*    irp[4];

if (nrhs==0) mexErrMsgTxt("no argument");

buflen = (mxGetM (prhs[0]) * mxGetN(prhs[0])) + 1;
buf = mxCalloc (buflen, sizeof(char));
error = mxGetString (prhs[0], buf, buflen);
if (cF==-1) { for (cF=0;cF<20;cF++) {FI[cF].nRows=0;}; cF=0; }

 if ((!strcmp("symbfact",buf)) || (!strcmp("numfact",buf)) || (!strcmp("fact",buf)) || (!strcmp("clear",buf)) ) { 
  if (!mxIsInt32(prhs[2])) mexErrMsgTxt("bad param type, clear");
    param=(int*)mxGetData(prhs[2]); 
} else if (!strcmp("solve",buf)) {
  if (!mxIsInt32(prhs[3])) mexErrMsgTxt("bad param type,solve");
  param=(int*)mxGetData(prhs[3]);
} else if (!strcmp("cvs",buf)) { 
 plhs[0]=mxCreateString("$Revision: 1.8 $  $Date: 2009/10/30 08:10:53 $");
 return;
} else {
  mexErrMsgTxt("bad command\n");
}

maxfct = param[100]; mnum = param[101];  
msglvl = param[102]; error = param[103]; nonsym =param[104]; 


/*----------------------------------------------------------- ismex */
if (!strcmp("ismex",buf))  {

 plhs[0] = mxCreateDoubleMatrix (1, 1, mxREAL);
 a = mxGetPr (plhs[0]);  *a=1.0;

/* --------------------- symbolic factorization------------*/
} else if (!strcmp("fact",buf))  {

 mexErrMsgTxt("Use symbfact ans numfact instead");

/* --------------------- symbolic factorization------------*/
} else if (!strcmp("symbfact",buf))  {

   cF=0; while (cF<20&&FI[cF].nRows!=0) { cF++; }
   if (cF==20) {mexErrMsgTxt("20 factors is the max, use pardiso clear");}
   cF=0;
   /* xxx can we have more than one pardiso factor ? */
   if (mxIsComplex(prhs[1])) mtype=13; /* complex, unsymmetric */

   plhs[0] = mxCreateDoubleMatrix (1,1, mxREAL);
   *mxGetPr(plhs[0])=(double)cF;
   /* PARDISOINIT(FI[cF].PardisoPt, &mtype, param);*/

   sd_mkl_get_sparse(prhs[1],irp,dims);
   columns=irp[0];rows=irp[1];q=irp[2];qi=irp[3];
   nCols=dims[0];FI[cF].nRows = nCols;

   if (mtype==13) { /* complex case */

     nZmax=rows[nCols];     /* warning pardiso permutes rows/columns */
     FI[cF].columns_1 = mxCalloc(nZmax,sizeof(int)); /* unsymnetric */

     mexMakeMemoryPersistent(FI[cF].columns_1);
     FI[cF].rows_1  = mxCalloc(FI[cF].nRows+1,sizeof(int));
     mexMakeMemoryPersistent(FI[cF].rows_1);
     FI[cF].pr =  mxCalloc(2*nZmax,sizeof(double)); /* unsymnetric */
     mexMakeMemoryPersistent(FI[cF].pr);

     j3=0; FI[cF].rows_1[0]=1;
     for(j1=0; j1<FI[cF].nRows; j1++)  {
       for (j2=rows[j1]; j2<rows[j1+1]; j2++) {  
	 FI[cF].columns_1[j3/2]  = columns[j2]+1;
         FI[cF].pr[j3]=q[j2];  j3++;
         FI[cF].pr[j3]=qi[j2]; j3++;
       }  
       FI[cF].rows_1[j1+1]=j3/2+1;  
     }

   } else if (nonsym == 1) {

     nZmax=rows[nCols];  /* warning pardiso permutes rows/columns */
     FI[cF].columns_1 = mxCalloc(nZmax,sizeof(int)); /* unsymnetric */
     mexMakeMemoryPersistent(FI[cF].columns_1);
     FI[cF].rows_1  = mxCalloc(FI[cF].nRows+1,sizeof(int));
     mexMakeMemoryPersistent(FI[cF].rows_1);
     FI[cF].pr =  mxCalloc(nZmax,sizeof(double)); /* unsymnetric */
     mexMakeMemoryPersistent(FI[cF].pr);

     j3=0; FI[cF].rows_1[0]=1;for(j1=0; j1<FI[cF].nRows; j1++)  {
       for (j2=rows[j1]; j2<rows[j1+1]; j2++) {  
	 FI[cF].columns_1[j3]  = columns[j2]+1;FI[cF].pr[j3]=q[j2]; j3++;  
       }  
       FI[cF].rows_1[j1+1]=j3+1;  
     }     
     mtype = 11;
   } else {

     nZmax=0; 
     for(j1=0; j1<nCols; j1++)  {
      for(j2=rows[j1]; j2<rows[j1+1]; j2++) {
       if (columns[j2]>=j1) { nZmax++; }
       }
     }

     FI[cF].columns_1 = mxCalloc(nZmax*2,sizeof(int));
     mexMakeMemoryPersistent(FI[cF].columns_1);
     FI[cF].rows_1  = mxCalloc(FI[cF].nRows+1,sizeof(int));
     mexMakeMemoryPersistent(FI[cF].rows_1);
     FI[cF].pr =  mxCalloc(nZmax*2,sizeof(double));
     mexMakeMemoryPersistent(FI[cF].pr);

     j3=0; FI[cF].rows_1[0]=1;for(j1=0; j1<FI[cF].nRows; j1++)  {
       for (j2=rows[j1]; j2<rows[j1+1]; j2++) {   
	 if (columns[j2]>=j1) {   
	   FI[cF].columns_1[j3]  = columns[j2]+1;FI[cF].pr[j3]=q[j2]; j3++;   
	 }  
       }   
       FI[cF].rows_1[j1+1]=j3+1;   
     }
     mtype = -2;
   }
msglvl=1;
   if (msglvl>0) {
     mexPrintf("Factor %d : nNonZeros= %i  nRows=%i  \n", cF, nZmax, FI[cF].nRows);
   }

/* -------------------------------------------------------------------- */
/* .. Initialize the internal solver memory pointer. This is only */
/* necessary for the FIRST call of the PARDISO solver. */
/* -------------------------------------------------------------------- */
   for (i = 0; i < 64; i++) {  FI[cF].PardisoPt[i] = 0; }
/* -------------------------------------------------------------------- */
/* .. Reordering and Symbolic Factorization. This step also allocates */
/* all memory that is necessary for the factorization. */
/* -------------------------------------------------------------------- */
   phase = 11; nRhs=0;
  mexPrintf("zzza --------------- pt %p %i",FI[cF].PardisoPt,FI[cF].PardisoPt[64]);

   PARDISO (FI[cF].PardisoPt, &maxfct, &mnum, &mtype, &phase,
	    &(FI[cF].nRows), FI[cF].pr, FI[cF].rows_1, FI[cF].columns_1, &idum, &nRhs,
	    param, &msglvl, &ddum, &ddum, &param[103]);

   if (param[103] != 0) return;
   if (msglvl>0) {
     mexPrintf("Reordering and symbolic factorization completed ...\n");
     mexPrintf("Memory used in symbfact: %d ko\n",param[14]+param[15]);
   }   

mexPrintf("zzz --------------- pt %p",FI[cF].PardisoPt);

 }

/* --------------------- nemerical factorization----------------*/
 else if (!strcmp("numfact",buf))  {

/* -------------------------------------------------------------------- */
/* .. Numerical factorization. */
/* -------------------------------------------------------------------- */
 /*
clear functions
sd('pardiso')
k=speye(10);b=ones(size(k,1),1);
ofact('method pardiso')
kd=ofact('symbfact',k);
kd = ofact('numfact',k,kd); % at each factorization
q=kd\b;

k=sparse([2 0 0;0 3 0;0 0 1])+1e-20*i;
b=[1 1+1e-17*i 1]';

  */
  sd_mkl_get_sparse(prhs[1],irp,dims);
  columns=irp[0];rows=irp[1];q=irp[2]; qi=irp[3];
    
  phase = 22; cF=0;
  if (mxIsComplex(prhs[1])) mtype=13; /* complex, unsymmetric */
  if (mtype==13) { /* complex case */
   j3=0;
   for(j1=0; j1<FI[cF].nRows; j1++)  {
     for (j2= FI[cF].rows_1[j1]-1;j2<(FI[cF].rows_1[j1+1]-1);j2++) {
       FI[cF].pr[j3]=q[j2];  j3++;  
       FI[cF].pr[j3]=qi[j2]; j3++;  
     }
   }

  } else if (nonsym == 1) {
   mtype=11;
   j3=0;
   for(j1=0; j1<FI[cF].nRows; j1++)  {
     for (j2= FI[cF].rows_1[j1]-1;j2<(FI[cF].rows_1[j1+1]-1);j2++) {
       mexPrintf("%i %i %i %g\n",j3,j1,j2,q[j2]);
       FI[cF].pr[j3]=q[j2]; j3++;  
     }
   }
  } else {
   mtype = -2;
   j3=0;
   for(j1=0; j1<FI[cF].nRows; j1++)  {
     for (j2=rows[j1]; j2<rows[j1+1]; j2++) {   
       if (columns[j2]>=j1) {   
	 FI[cF].pr[j3]=q[j2]; j3++;   
       }  
     }   
   }
 }
mexPrintf("\nnumfact zzz --------------- pt %p",FI[cF].PardisoPt);

 PARDISO (FI[cF].PardisoPt, &maxfct, &mnum, &mtype, &phase,
	  &(FI[cF].nRows), FI[cF].pr, FI[cF].rows_1, FI[cF].columns_1, &idum, &nRhs,
	  param, &msglvl, &ddum, &ddum, &param[103]);
mexPrintf("\nnumfact zzz --------------- pt %p",FI[cF].PardisoPt);

 plhs[0] = mxCreateDoubleMatrix (1,1, mxREAL);
 if (param[103] != 0) {
   mexPrintf("\nERROR during numerical factorization with code= %d\n",param[103]);
   if (param[103] == -2){
     mexPrintf("memory error code returned, try continue...\n");
     mexPrintf("Memory used: %d ko\n",param[16]);
     param[103] = 0;
   }
   else {
     return;
   }
 }
 /* { mexErrMsgTxt("\nERROR during numerical factorization"); }*/
 if (msglvl>0) {
   mexPrintf("Factor %e :     nNonZeros= %i    nRows=%i \n", 
              FI[cF].pr[1],  FI[cF].rows_1[FI[cF].nRows]-1, FI[cF].nRows);
   mexPrintf("Factorization completed ...\n");
   mexPrintf("Memory used in numfact: %i ko\n",param[16]);
 }
 *mxGetPr(plhs[0])=(double)cF;
mexPrintf("zzz------------- pt %p",FI[cF].PardisoPt);

 }

/*--------------------------------------------------------- solve */
else if (!strcmp("solve",buf)) {

 double *b, *bi, *qr, *qri; 

 if (nrhs<3) mexErrMsgTxt("Symbolic factor must be given");
 cF= (int)(*mxGetPr (prhs[2])); /* factor number */
 if (FI[cF].nRows==0) mexErrMsgTxt("Factor is not defined");

 nRhs=mxGetN(prhs[1]);

 if (mxIsComplex(prhs[1])) { /* xxx only complex matrix and complex rhs at the moment */
   mtype=13;
   b=mxGetPr(prhs[1]); bi=mxGetPi(prhs[1]);  
   rhs=mxCalloc (2*mxGetM(prhs[1])*mxGetN(prhs[1]), sizeof(double));
   for(j1=0; j1<mxGetM(prhs[1])*mxGetN(prhs[1]); j1++) {
     rhs[2*j1]=b[j1]; rhs[2*j1+1]=bi[j1];
   }
   q=mxCalloc (2*mxGetM(prhs[1])*mxGetN(prhs[1]), sizeof(double));

   plhs[0] = mxCreateDoubleMatrix (mxGetM(prhs[1]),mxGetN(prhs[1]), mxCOMPLEX);
   qr=mxGetPr(plhs[0]); qri=mxGetPi(plhs[0]);

 } else {

   rhs = mxGetPr(prhs[1]);  
   plhs[0] = mxCreateDoubleMatrix (mxGetM(prhs[1]),mxGetN(prhs[1]), mxREAL);
   q=mxGetPr(plhs[0]);
 }
/* -------------------------------------------------------------------- */
/* .. Back substitution and iterative refinement. */
/* -------------------------------------------------------------------- */

 if ((mxIsComplex(prhs[1])) && (nonsym==1)) {
   mtype=13;
 }
 else if (nonsym==1){
   mtype=11;
 }
mexPrintf("zzz------------- pt %p",FI[cF].PardisoPt);
 phase = 33;
 PARDISO (FI[cF].PardisoPt, &maxfct, &mnum, &mtype, &phase,
          &(FI[cF].nRows), FI[cF].pr, FI[cF].rows_1, FI[cF].columns_1, &idum, &nRhs,
          param, &msglvl, rhs,q, &param[103]);
 if (param[103] != 0) { mexErrMsgTxt("\nERROR during back substitution"); }


 if (mxIsComplex(prhs[1])) { /* xxx only complex matrix and complex rhs at the moment */

   for(j1=0; j1<mxGetM(prhs[1])*mxGetN(prhs[1]); j1++) {
     qr[j1]=q[2*j1]; qri[j1]=q[2*j1+1];
   }
   mxFree(q); mxFree(rhs);
 }


}/*--------------------------------------------------------- clear */
else if (!strcmp("clear",buf)) {

  if (nrhs<2) {
   cF=-1; if (msglvl>0) mexPrintf("Clearing all pardiso factors");
  } else {
    cF= (int)(*mxGetPr (prhs[1]));
   if (cF<0) if (msglvl>0) mexPrintf("Clearing all pardiso factors");
   else { 
       if (FI[cF].nRows==0) mexErrMsgTxt("Cannot clear undefined factor");
       if (msglvl>0) mexPrintf("Clearing factor %i\n",cF);     
     }
  }
  ClearFactor(cF, maxfct,mnum,mtype,param,msglvl);


}/* end command selection ---------------------------------------- */


 mexAtExit(ClearAllFactors);

 return; 

} /* EOF */

