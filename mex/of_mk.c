#include <stdio.h>
#include <string.h>
#ifdef OFOMP
#include <omp.h>
#endif
#include "mex.h"
#include "math.h"

extern char* pre_cvs();
#define size_3 3	  /* dimension of matrix */

#include "of_EltConst.h"
#include "of_mk_pre.h"

#ifdef MatlabVER /* matlab compile */
#include "matrix.h"
#include "of_def.h"
#include "of_EltConst.h"
#include "hyper.h"

#else /* Scilab compile */
#include "stack-c.h"
#include "of_def.h"
#include "of_mk_pre.c"
#endif

int nMiss, nMissAlloc;
        
void mexFunction (int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]) 
/*-----------------------------------------------------------------------*/
{
  double *out, *out1, *out2;
  char CAM[15], ierr[100]={'\0'};
  int i1, typ;

  /*-----------------------------------------------------------------------*/
  /* INIT used by all                                                      */

  if (nrhs==0)
    {
      plhs[0] =  mxCreateDoubleMatrix(1,1,mxREAL);
      *mxGetPr(plhs[0])=5.0003;
      return;
    }

  i1 = (int)mxGetNumberOfElements(prhs[0]);
  mxGetString(prhs[0],CAM,15);
  /* CAM = mxCalloc(i1+1,sizeof(char)); mxGetString(prhs[0],CAM,i1+1);*/

  if ((nrhs==1) && (!strcmp("cvs",CAM))) {
    mwSize     dims[2];

    dims[0]=3; dims[1]=2; plhs[0]=mxCreateCellArray((mwSize)2,dims);
    mxSetCell(plhs[0],0,mxCreateString("of_mk"));
    mxSetCell(plhs[0],1,mxCreateString("of_mk_subs"));
    mxSetCell(plhs[0],2,mxCreateString("of_mk_pre"));
    mxSetCell(plhs[0],3,mxCreateString("$Revision: 1.245 $  $Date: 2020/05/27 15:23:29 $"));
    mxSetCell(plhs[0],4,mxCreateString(pre_cvs()));
    mxSetCell(plhs[0],5,pre_cvs2());

    return;
  } else if ((nlhs==1) && (nrhs==1) && (!strcmp("mwIndex",CAM))) {
      plhs[0] =  mxCreateDoubleMatrix(1,1,mxREAL);
      out = mxGetPr(plhs[0]);
      out[0]=(double)sizeof(mwIndex);
      return;
  }  else if ((nlhs==0) && (nrhs==1) && (!strcmp("mwIndex",CAM))) {
#if MatlabVER<=73 /* NEW MATLAB 7.3 SIZE POINTERS */
     mexPrintf("SINCE73 : NO\n");
#else 
     mexPrintf("NO SINCE73 : YES\n");
#endif
 
      return;
  }

  /* INIT and other non standard calls. TEST IS STRING WITH MORE THAN 8 char
   -------------------------------------------------------------- */
  if (i1>=8){

    if (!strcmp(CAM,"mitcinit")) {
	 mexErrMsgTxt ("Obsolete.");
    } else if (!strcmp("xkx_trans",CAM)) {
    /* k=of_mk('xkx_trans',x,k);  */
    double *offset;
    int    type;
	mwIndex ym,yn;

    if(nlhs!=1) mexErrMsgTxt ("One output required.");
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
    if (mxIsComplex (prhs[1])) 
          mexErrMsgTxt ("Matrix of basis change must not be complex.");
     if (mxIsComplex (prhs[2])) { 
      #if MatlabVER >= 904
		mexErrMsgTxt("Error need complex reimplement");
      #else
	if (mxIsComplex(prhs[1])) {
          x_k_x(mxGetPr(prhs[1]),mxGetPi(prhs[2]),mxGetPi(plhs[0]),
                                               offset,type,(int)ym,(int)yn);
	}
      #endif
    }
    /* call the subroutine xkx */
    x_k_x (mxGetPr (prhs[1]),mxGetPr(prhs[2]),mxGetPr(plhs[0]),
                                              offset,type,(int)ym,(int)yn);
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    }  else if (!strcmp(CAM,"princstress"))  {
        double    *A, *b, c[size_3]; /*DUMMY[1][1], WORK[3*size_3];*/
        int       j1, ok=1, c1=size_3, c2=3*size_3, c3=1;
        char      c4='N';
	/*
              of_mk('princstress',eye(3))
 	*/
        if (mxIsSparse(prhs[1])) mexErrMsgTxt("Full matrix needed.");
        A=mxGetPr(prhs[1]);
        if ((mxGetM(prhs[1])!=3) || (mxGetN(prhs[1])!=3)) 
                                 mexErrMsgTxt("Only 3x3 matrix");

        plhs[0]= mxCreateDoubleMatrix(size_3,1,mxREAL); b=mxGetPr(plhs[0]);

	/* of_dgeev(&c4,&c4,&c1,A,&c1,b,c,DUMMY,&c3,DUMMY,&c3, WORK, &c2, &ok);  */
        /*   'N','N',3,input, 3,eigen,VL,  1,   VR,    1  , out, 2*3, out,info)
              SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
                                $                  LDVR, WORK, LWORK, INFO ) */
	/*
ALTERNATIVE IS:

      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DSYEV computes all eigenvalues and, optionally, eigenvectors of a
*  real symmetric matrix A.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of A contains the
*          upper triangular part of the matrix A.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of A contains
*          the lower triangular part of the matrix A.
*          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
*          orthonormal eigenvectors of the matrix A.
*          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
*          or the upper triangle (if UPLO='U') of A, including the
*          diagonal, is destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= max(1,3*N-1).
*          For optimal efficiency, LWORK >= (NB+2)*N,
*          where NB is the blocksize for DSYTRD returned by ILAENV.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the algorithm failed to converge; i
*                off-diagonal elements of an intermediate tridiagonal
*                form did not converge to zero.
*
*  =====================================================================
*/

        for (j1=0; j1<size_3; j1++) { 
          if (c[j1]!=0 ) mexErrMsgTxt("Complex eigenvalue exist");
        }
        if (ok!=0) mexErrMsgTxt("Problem to compute eigenvalues");
        /* for (j1=0; j1<size_3; j1++) mexPrintf("%f  %f \n", b[j1],c[j1]);*/


/* -----------------------------------------------------------------------*/
    }  else if (!strcmp(CAM,"buildndn"))  {

        /* rules=of_mk('buildndn',13,integrules,nodeE);
           rules=of_mk('buildndn',2 or 3 or 23,integrules,nodeE);
           rules=of_mk('buildndn',2 or 3 or 23,out,x);               */
      mxArray       *mEC,*InfoAtNode; 
      int           Nw, Nnode, NfieldE, t1, *ti, GSize[20];
      struct EltConst EC;
 
      struct GroupFields GF;

      double*     NeedFree[6]={NULL,NULL,NULL,NULL,NULL,NULL};

        if (mxIsInt32(prhs[1])) { ti = (int*)mxGetData(prhs[1]);t1=ti[0];}
        else                    { t1 = (int)(*mxGetPr(prhs[1])); }

       mEC=(mxArray*)prhs[2]; InfoAtNode=NULL;  
   
        GF.N=NULL; GF=initGroup(GF,prhs[3],mEC,GSize,NeedFree);
	   if (nrhs<4) {EC.nodeE=NULL; GSize[4]=0; GSize[0]=0;
	   } else { /* nodeE assembled externally, should dissapear for EC.nodeE */
	    mexErrMsgTxt("Should use EC.nodeE and EC.nodeEt ");
        EC.nodeE=mxGetPr(prhs[3]); GSize[0]=(int)mxGetM(prhs[3]);GSize[4]=(int)mxGetN(prhs[3])-4;
	   }
       EC=initInfoAtNode(InfoAtNode,mEC,EC,GF,GSize,NeedFree); 
       Nw =GSize[2];

       Nnode=GSize[0];NfieldE=GSize[4]+4;
        if (    mxGetM(mxGetField(prhs[2], 0,"NDN"))
             > mxGetN(mxGetField(prhs[2], 0,"N"))  )
          mexErrMsgTxt("Problem of row allocation in NDN ");

	if (t1==3 || t1==32|| t1==31) {
          if (Nw*4!= mxGetN(mxGetField(prhs[2], 0,"NDN"))) 
            mexErrMsgTxt("Problem of column allocation in NDN");
    } else if (t1==2 || t1==23) {
          if (Nw*3!= mxGetN(mxGetField(prhs[2], 0,"NDN"))) 
            mexErrMsgTxt("Problem of column allocation in NDN");
    } else if (t1==13 || t1==12) {
          if (Nw*2> mxGetN(mxGetField(prhs[2], 0,"NDN"))) 
            mexErrMsgTxt("Problem of column allocation in NDN");
    } else if (t1==0) { return;
	} else mexErrMsgTxt("Not a valid NDN mode");


    if (nrhs>4) mexErrMsgTxt("Expecting EltConst.nodeEt "); 
    NDNSwitch(t1,GF,&EC,Nw,Nnode,NfieldE);
	/* does not account for variability withing element */
    if (GF.CTable!=NULL) { int jw=0;
	 constitInterp(GF,&EC,jw,Nnode,GSize);
	} 

/* ------------------------------------------------------------------------*/
	/* CONDENSE :

 ke=t3p('testmat'); me=ke{2}; ke=ke{1}; [a,b]=elem0('condense',5,ke,me);
 [kr,mr]=of_mk('condense',5,ke,me);
 norm(kr-a)
 norm(mr-b)

*/
}  else if (!strcmp(CAM,"condense"))  {

  int i1, M;
  double *ke, *me, *kr, *mr;

  if (nrhs!=4) mexErrMsgTxt("4 input arguments needed");
  i1=(int)(*mxGetPr(prhs[1]))-1;
  M=(int)mxGetM(prhs[2]);   /* xxx check sizes */
  ke=mxGetPr(prhs[2]);
  me=mxGetPr(prhs[3]); /* xxx add check */

  plhs[0]= mxCreateDoubleMatrix(i1,i1,mxREAL); kr=mxGetPr(plhs[0]);
  if (mxGetM(prhs[3])>0) {
    plhs[1]= mxCreateDoubleMatrix(i1,i1,mxREAL); mr=mxGetPr(plhs[1]);
  } else { mr=NULL;  }

 
  elt_condense(i1, ke, kr, me, mr, M);  


/* ------------------------------------------------------------------------*/
      /* applies the matrix integration rule defined in R1     
        ke=of_mk('matrixintegration',jElt[1],NodePos[2],Case.Node[3], ...
            pointers[4],integ[5],constit[6],gstate[7],elmap[8],
            InfoAtNode[9],EltConst[10],def.def)
           of_mk('matrixintegration',DofPos[1],NodePos[2],Case.Node[3], ...
            pointers[4],integ[5],constit[6],gstate[7],elmap[8],
            InfoAtNode[9],EltConst[10],def.def[11],k[12],opt[13])

      */
}  
else if (!strcmp(CAM,"setomppro"))  {
#ifdef OFOMP 
  if (nrhs>1) {
   omp_set_num_threads((int)*mxGetPr(prhs[1]));
  }
  mexPrintf("Max_Threads = %i\n",omp_get_max_threads());
#else
  mexErrMsgTxt("of_mk was compiled without OFOMP");
#endif
}
 else if (!strcmp(CAM,"matrixintegrat"))  {

        /*  k=zeros(8); of_mk('matrixintegration',rules,r1,constit,k); */
        mxArray     *field,  *mEC;
        const mxArray *CallRHS[11];
        double      *node, *o_pr, zero=0, *Ener;
        double*     NeedFree[6]={NULL,NULL,NULL,NULL,NULL,NULL};
        int*       pEC=NULL;
        int         j1, j2, j3, Nterms, Nw, Mk, 
                    Nelt, jMat,  *opt, *CurDofPos,
                    *elmap, *DofPos, DofPerElt, StrategyType,ThreadiPos[10];
        int         jElt, jEl0, *NodePos, Ndof, Nnode, NfieldE,
                    Nmnode,*pointg, 
            	    Mdef=0, Ndef=0, GSize_ref[20],*GSize,
                    *cEGI, Mener, *rrule=NULL, Nrule=0, Mrule=0,*ind, Nk,nul=0,unit=1,ist=1;
        mwIndex     *o_ir, *o_jc;
        struct GroupFields GF;

        struct EltConst  EC;

	/* constitutive energy function handle (Integ,Constit,I,dWdI,d2WdI2) */


#define of_inc_fs "../mex/of_mk_fsmatrix.c"
#define MatrixIntegrationStep 0
#include of_inc_fs
#ifdef of_inc_user
#  include of_inc_user
#endif

     /* inputs */
     NodePos=(int*)mxGetData(prhs[2]); jEl0=0;Ndof=0; Ener=NULL;
     if (!mxIsInt32(prhs[2])) mexErrMsgTxt("NodePos must be int32");

     Nmnode=(int)mxGetM(prhs[3]); /* OffSet for rows in node matrix */ 
     node=mxGetPr(prhs[3]); 
     for (j1=0;j1<6;j1++) { if (NeedFree[j1]!=NULL) mexErrMsgTxt("alloc error");}

     mEC=(mxArray*)prhs[10];
     if (!mxIsStruct(mEC)) mexErrMsgTxt("EltConst must be a structure");
     if (mxIsEmpty(mEC)) mexErrMsgTxt("EltConst must not be empty");

     field = mxGetField(mEC, 0,"GSize");
     if (field!=NULL) {GSize=(int*)mxGetData(field);} else {GSize=GSize_ref;}
     GSize[0]=(int)mxGetM(prhs[2]);Nnode=GSize[0];GSize[5]=0;

     GF.N=NULL; GF=initGroup(GF,prhs[3],mEC,GSize,NeedFree);     
	 if (GF.Ncondense) {
       field = mxGetField(mEC, 0,"Ncondense");
       ind=(int*)mxGetData(mxGetField(mEC, 0,"CondenseInd"));
       Nk=(int)mxGetM(mxGetField(mEC, 0,"CondenseMat"));
	 }

     EC.nodeE=NULL; GSize[4]=0; /* nodeE possibly handled later */
     EC=initInfoAtNode(prhs[9],mEC,EC,GF,GSize,NeedFree); 

     field=mxGetField(mEC, 0,"RhsDefinition");
     if (field!=NULL) { 
       rrule=(int*)mxGetData(field);Mrule=(int)mxGetM(field);Nrule=(int)mxGetN(field);
     }
     if (nrhs>12) { /* direct assembly of a complete group  - - - - - - - - */

       Nelt=(int)mxGetN(prhs[1]); 
       opt=(int*)mxGetData(prhs[13]); DofPerElt=opt[0]; 
       elmap=(int*)mxGetData(prhs[8]);DofPos=(int*)mxGetData(prhs[1]);Ndof=(int)mxGetM(prhs[1]);
       if (opt[1]==-1) { /* xxx ENER GetField(cEGI), GetField(Ener) */
         cEGI=(int*)mxGetData(mxGetField(prhs[12], 0,"cEGI"));
         Ener=mxGetPr(mxGetField(prhs[12], 0,"Ener"));
         Mener=(int)mxGetM(mxGetField(prhs[12], 0,"Ener"));
         rrule=NULL;
       } else {
        o_pr=mxGetPr(prhs[12]); o_ir=mxGetIr(prhs[12]); o_jc=mxGetJc(prhs[12]); 
       }
       #if MatlabVER >= 904
       if (mxIsComplex(prhs[11])) {
	    GF.def=(double*)mxGetData(prhs[11]); GF.defi=GF.def; ist=(int)2;
	    /*{mexErrMsgTxt("Error need complex reimplement");}*/
       }  else {GF.def=(double*)mxGetDoubles(prhs[11]);GF.defi=NULL;}
       #else
	    GF.def=mxGetPr(prhs[11]); GF.defi=mxGetPi(prhs[11]);
       #endif
       if (mxIsSparse(prhs[11])) mexErrMsgTxt("def.def must be full");
       Mdef=(int)mxGetM(prhs[11]);Ndef=(int)mxGetN(prhs[11]);/* mexPrintf("%i %i %p\n",Mdef,Ndef,def);*/
       if (GF.def!=NULL&&Ndef>1) {GF.RHS=GF.def+ist*Mdef;} else {GF.RHS=NULL;}
     } else { /* assembly of a single element matrix */
       jEl0=(int)*mxGetPr(prhs[1])-1; Nelt=jEl0+1;
       opt=NULL; rrule=NULL;GF.RHS=NULL;
     } /*  assembly of a single matrix or a full group - - - - - - - - - - - */
     /* initialize element output matrices - - - - - - - - - - - - - - */

     pointg   = (int*)mxGetData(prhs[4]);
     /* always use first col. for matrix allocation */
     EC.integ   = (int*)mxGetData(prhs[5]);Mk=EC.integ[2]+GF.Ncondense; 
     EC.integ+=pointg[5]; if (Mk==0) mexErrMsgTxt("Mk cannot be 0");
     field = mxGetField(mEC,0,"ke"); if (field!=NULL) { EC.ke=mxGetPr(field); } 
	 field = mxGetField(mEC,0,"me");    
     if (opt==NULL && pointg[4]==0) { /* two matrices M,K return as outputs */
       if (nlhs==1) mexErrMsgTxt("m and k should be returned for MatDes 0");
     } 
 /*general inits that may be modified depending on assembly stratety*/
 jMat=0; StrategyType=0;   GF.rule_terms=NULL;Nterms=0;
 GF.ke=(double*)ofMalloc(Mk*Mk*sizeof(double)); 
 plhs[0]=NULL; plhs[1]=NULL; 

 while (jMat<4) { /* matrices to generate allow for multiple */

   if      (pointg[4]==0 &&jMat==0) { pointg[4]=2;jMat=2;  /* mass */
   } else if (pointg[4]==2 && jMat==2) { /* stiffness assembled in mass+stiff*/
       pointg[4]=1;        jEl0=(int)*mxGetPr(prhs[1])-1;jMat=3;
   } else if (pointg[4]==1 &&jMat==3) { /* exit when point[4]==0 and m and k done */ 
      pointg[4]=0;jMat=4;break; /*exit*/
   } else { jMat=4; } /* standard case one generates one element matrix */	  
   /* prepare for different matrix assembly strategies - - - - - - - - - - */
   /* generic multiphysic linear element unless no match in MatrixInt ...*/
   field = mxGetField(mEC,0,"MatrixIntegrationRule");
   if ((field==NULL)||(mxGetM(field)==0)||pointg[4]>mxGetN(field)) { 
      /* no match for matrule */
      field = mxGetField(mEC,0,"material");
      if (field==NULL) mexErrMsgTxt("Not a supported MatrixIntegrationRule 1");
   } else {
    field=mxGetCell(field,pointg[4]-1); /* rule for given matrix type */     
    /* MatrixRule{j1} defined something to be done */
    if (field!=NULL&&mxGetM(field)!=0) {
       GF.rule_terms   = (int*)mxGetData(field); Nterms=(int)mxGetM(field);
       if (pointg[3]==31) { /* this rule defined by INRIA only for hyperelastic */
         for (j1=0;j1<Nelt;j1++){ pointg[3+j1*mxGetM(prhs[4])]=3;}
       }
       StrategyType=1;field=mxGetField(mEC,0,"fHandle");
       if (field!=NULL) field=mxGetField(mEC,0,"material"); /* if fHandle force callback */
    } else  {field = mxGetField(mEC,0,"material");Nterms=0;}
   }
   if (GF.RHS!=NULL) {
	    Ndof=Mk; GF.NBe=Ndof;
		if (rrule!=NULL && StrategyType==0) StrategyType=1; /* rhs given */
       }
   if (field!=NULL) { /* a material field is defined */
     mxGetString(field,CAM,15);
     if (!strcmp("Elastic3DNL",CAM)) { /* geometric non linear elastic 3D */
       if (pointg[4]!=5) { StrategyType=1; /* lin_multi */
        if (Nterms==0) jEl0=Nelt; /* empty rule skip */
       } else {
        StrategyType=2;  
        Ndof = 3*GSize[0];GF.NBe=Ndof; GF.NdefE=Ndof*Ndef;
       }
     } else if (!strcmp("callback",CAM)) { /* comes back to matlab */
       StrategyType=5; 
       field = mxGetField(mEC,0,"defe"); GF.NdefE=0; /* no alloc possible for callback */
	   if (field!=NULL) {
		   EC.defE=mxGetPr(field);Ndef=(int)mxGetN(field);
	   }
       CallRHS[0]=mxGetField(mEC,0,"fHandle");
       if (CallRHS[0]==NULL) mexErrMsgTxt("fHandle not defined");
       CallRHS[1]=prhs[4];/* [0] function_pointer [1] pointers, integ, constit, gstate elmap */
       CallRHS[2]=prhs[5];CallRHS[3]=prhs[6];
       CallRHS[4]=prhs[7];CallRHS[5]=prhs[8];
       CallRHS[6]=prhs[9]; /* [6] InfoAtNode, EltConst */
       CallRHS[9]= mxCreateDoubleMatrix(1,2,mxREAL);
       if (nrhs<11)  CallRHS[7]=mxCreateDoubleMatrix(0,0,mxREAL);
       else CallRHS[7]=mEC; /* EltConst */
       field = mxGetField(mEC,0,"constit"); 
       if (field!=NULL) EC.constit=mxGetPr(field);

       if (nrhs<12)  CallRHS[8]=mxCreateDoubleMatrix(0,0,mxREAL);
       else CallRHS[8]=prhs[11];
       CallRHS[10]= prhs[1]; /* DofPos */
     } else if (!strcmp("hyperelastic",CAM)) { /* hyperelastic */
       StrategyType=3; /* MatrixIntegrationRule must be empty for desired mat */
       Ndof = 3*GSize[0];
       GF.NBe=7*Ndof;
       GF.NdefE=Ndof*Ndef;
     } 
#undef MatrixIntegrationStep 
#define MatrixIntegrationStep 1
#include of_inc_fs
#ifdef of_inc_user
#   include of_inc_user
#endif
     else if (StrategyType==0) {
       mexPrintf(".material='%s'",CAM);mexErrMsgTxt(" not supported");
     }
   }
   if (mxIsStruct(prhs[7])) {
       field=mxGetField(prhs[7],0,"Y"); 
       if (field==NULL) mexErrMsgTxt(" missing gstate.Y field");
       EC.gstate=mxGetPr(field);GSize[5]=(int)mxGetM(field); 
   } else { EC.gstate=mxGetPr(prhs[7]);GSize[5]=(int)mxGetM(prhs[7]); }
   
   Nnode=GSize[0];
   GSize[6]=(int)mxGetM(prhs[6]); Nw=GSize[2];
   InitThreadiPos(ThreadiPos,GSize,Mk);
/* OMP execution of the assembly loop - - - - - - - - - - - - - - - - -
 Private variables
 jElt, j1, j2, i1, i2  % loop indices
 point, integ, constit, CurDofPos  % pointers to read only data
 Variables that need to be reallocated for each thread
   nodeEj[GSize[0]*(4+NInfoAtNode)]
   NDNj[NdofPerField*(Nw*Ndim)]
   jdetj[Nw], basj[4*Nw]
   kej, Bej, defej
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

#undef MatrixIntegrationStep 
#define MatrixIntegrationStep 11
      #ifdef of_inc_user
      #   include of_inc_user
      #endif

pEC=(int*)FieldPointer(mEC,(char*)"pEC"); 

/*callback -> one thread not using critical to allow for EC.NodeE */
if (StrategyType==5) { /* matlab one thread, uses EC */
  double *ke, *constitE=NULL;
  if (GF.CTable!=NULL&&EC.constit==NULL) {
	  constitE=(double*)ofMalloc(GSize[6]*sizeof(double)); EC.constit=constitE;
  }
/* -------------------------------------------------------------------------------*/
for (jElt=jEl0;jElt<Nelt;jElt++) { /* loop on elements to assemble */
/* warning jElt MUST REMAIN SET BEFORE HERE FOR SINGLE MATRIX ASSEMBLY */

  double *pjElt;
  struct EltConst ECe;
  int* point;
   mxArray *BackEC[1];

  CurDofPos=DofPos+jElt*DofPerElt;
  point   = (int*)mxGetData(prhs[4])+jElt*mxGetM(prhs[4]);/*pointers(:,jElt)*/
  point[4]= *((int*)mxGetData(prhs[4])+4); 
  if ((point[2]==-1&&point[4]!=2)||point[2]==-2) continue; /* -1 only mass, -2 deactivated*/
  ECe=EC; NfieldE=GSize[4]+4;
  if (GF.NBe!=0) {mxArray* f2;
    f2=mxGetField(mEC,0,"Be");if (f2!=NULL) ECe.Be=mxGetPr(f2); 
  }  else ECe.Be=NULL;
  ECe.integ   = (int*)mxGetData(prhs[5])+point[5];
  if (GF.CTable!=NULL) { memcpy(ECe.constit,mxGetPr(prhs[6])+point[6] , GSize[6] * sizeof( double ) );
  } else { ECe.constit = mxGetPr(prhs[6])+point[6];}

   if (ECe.defE!=NULL) { /* element deformation */
      for (j2=0;j2<Ndef;j2++) { for (j1=0;j1<Ndof;j1++) { 
	   ECe.defE[j1+Ndof*j2]=GF.def[CurDofPos[j1]+Mdef*j2];
	 }}
   } 
   if (ECe.ke==NULL) ECe.ke=(double*)ofMalloc(Mk*Mk*sizeof(double));
   for (j2=0;j2<Mk*Mk;j2++) ECe.ke[j2]=0.;
   if (ECe.Be!=NULL) {
       of_ptrdiff pMk, pnul=0, punit=1; 
       pMk=Mk;for (j2=0;j2<Mk;j2++) ECe.Be[j2]=0;/*of_dcopy(&pMk,&zero,&pnul,Bee,&punit);*/
   }
   if (pEC!=NULL) memcpy(pEC,&EC,40);/* sizeof(struct EltConst)); */

   buildNodeE(&ECe,node,GSize,NodePos,jElt,Nmnode);
   NDNSwitch(point[3],GF,&ECe,Nw,Nnode,NfieldE);

   pjElt=mxGetPr(CallRHS[9]);pjElt[0]=(double)jElt;pjElt[1]=(double)jMat;
   /* mxSetField(BackEC[0],0,"ke",mxCreateDoubleMatrix(Mk,Mk,mxREAL)); */
   mexCallMATLAB(1,BackEC,11,(mxArray**)CallRHS,"feval");
   ke=mxGetPr(mxGetField(BackEC[0],0,"ke"));
   if (ke!=ECe.ke) memcpy(ECe.ke,ke,ThreadiPos[5]*sizeof(double));/* copy to threadMem */
   if (pjElt[0]==-1) { 
	    lin_multi(GF,&ECe,Nterms,Mrule,Ndef,jElt,GSize,Mk,point,rrule);
	  }
    AssembleMatVec(&ECe,GF,opt,StrategyType,Mk,CurDofPos,o_ir,o_jc,o_pr,
             elmap,DofPerElt,Mdef,Ndef,Ener,cEGI,jElt,Mener);

} /* end of jElt loop on elements - - - - - - - - - - - - - - - - -*/ 
memcpy(GF.ke,ke,ThreadiPos[5]*sizeof(double));
if (constitE!=NULL) ofFree(constitE);


/* -------------------------------------------------------------------------------*/
} else { /* strategy type !=5 */


#pragma omp parallel /* #pragma omp parallel */
{ /* omp parallel */
  /* ECe=threadAlloc(GF,ECe,ThreadiPos,GSize); */
   int it;
   double* ECv; /* vector for alloc */
   struct EltConst* ECp;
  it=(sizeof(struct EltConst)/sizeof(double))+1;
  ECv=(double*)ofMalloc( (it + /* memory for ECe */
	   ThreadiPos[0]+ThreadiPos[1]+ThreadiPos[2]+ /* nodeE, NDN, jdet */
       GSize[6]+ThreadiPos[3]+ThreadiPos[4]+ThreadiPos[5]+GF.NBe+GF.NdefE)
	    *sizeof(double));
   ECp=(struct EltConst*)ECv;
   memcpy(ECp,&EC,sizeof(struct EltConst));
   ECp[0].nodeE=ECv+it;it+=ThreadiPos[0];
   ECp[0].NDN  = ECv+it;it+=ThreadiPos[1];
   ECp[0].jdet = ECv+it;it+=ThreadiPos[2];
   ECp[0].constit = ECv+it;it+=GSize[6];/* (GF.CTable!=NULL) material interpolation */
   if (ECp[0].bas!=NULL) {ECp[0].bas=ECv+it;it+=ThreadiPos[3]; } else ECp[0].bas=NULL;
   if (ECp[0].J!=NULL) {ECp[0].J=ECv+it;it+=ThreadiPos[3];} else ECp[0].J=NULL;
   ECp[0].ke=ECv+it;it+=ThreadiPos[5]; ECp[0].ke[0]=100e100;
   if (GF.NBe!=0) {ECp[0].Be=ECv+it;it+=GF.NBe;} else ECp[0].Be=NULL;
   if (GF.NdefE!=0) {ECp[0].defE=ECv+it;it+=GF.NdefE;} else ECp[0].defE=NULL;

 /* pragma omp for */
#pragma omp for private(j1,j2,j3,i1,CurDofPos) 
for (jElt=jEl0;jElt<Nelt;jElt++) { /* loop on elements to assemble */
/* warning jElt MUST REMAIN SET BEFORE HERE FOR SINGLE MATRIX ASSEMBLY */
 #ifdef verbose
   /* mexPrintf("jElt %i\n",jElt);*/
   if (jElt<0) mexErrMsgTxt("jElt<0");
   if (jElt>mxGetM(prhs[4])) mexErrMsgTxt("jElt>nElt");
 #endif

  int* point;
     CurDofPos=DofPos+jElt*DofPerElt;
     point   = (int*)mxGetData(prhs[4])+jElt*mxGetM(prhs[4]);/*pointers(:,jElt)*/
     if ((point[2]==-1&&point[4]!=2)||point[2]==-2) continue; /* -1 only mass, -2 deactivated*/
     point[4]= *((int*)mxGetData(prhs[4])+4);
     ECp[0].constit0=mxGetPr(prhs[6])+point[6];
     ECp[0].integ   = (int*)mxGetData(prhs[5])+point[5];
	 if (GF.CTable!=NULL) { memcpy(ECp[0].constit,mxGetPr(prhs[6])+point[6] , GSize[6] * sizeof( double ) );
	 } else { ECp[0].constit = mxGetPr(prhs[6])+point[6];}

  /* mexPrintf("num%i/%i pos %p\n",omp_get_thread_num(),MaxThread0,nodeEe); */
   if (ECp[0].defE!=NULL) { /* element deformation */
      for (j2=0;j2<Ndef;j2++) { for (j1=0;j1<Ndof;j1++) { 
	   ECp[0].defE[j1+Ndof*j2]=GF.def[CurDofPos[j1]+Mdef*j2];
	 }}
   } 
   for (j2=0;j2<Mk*Mk;j2++) ECp[0].ke[j2]=0.;
   if (ECp[0].Be!=NULL) {
       of_ptrdiff pMk, pnul=0, punit=1; 
       pMk=Mk;for (j2=0;j2<Mk;j2++) ECp[0].Be[j2]=0;/*of_dcopy(&pMk,&zero,&pnul,Bee,&punit);*/
   }

   buildNodeE(ECp,node,GSize,NodePos,jElt,Nmnode);/* xxx need Infinite handling xxx */
   /* fields at node of current element */
   if (EC.giNodePos==NULL) { /* cols gives field at element nodes, OBSOLETE */
    int MInfoAtNode;
	MInfoAtNode=GSize[3];
    for (j1=0;j1<MInfoAtNode/Nnode;j1++) { for (j2=0;j2<Nnode;j2++) {
     i1=j2+Nnode*jElt;
     ECp[0].nodeE[j2+(4+j1)*Nnode]=EC.InfoAtNode[i1+MInfoAtNode*j1]; 
    }}
   } 
   NfieldE=4+GSize[4];
       NDNSwitch(point[3],GF,ECp,Nw,Nnode,NfieldE);
#undef MatrixIntegrationStep 
#define MatrixIntegrationStep 2

     switch (StrategyType) {
     case 0: 
       EC.Be=NULL; /* robust energy computation */
	   break;
     case 1: { 
       /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
       GENERIC LINEAR MULTIPHYSIC ELEMENTS mat_og with MatrixIntegrationRule
       see implementation in of_mk_pre.c
      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */   
       lin_multi(GF,ECp,Nterms,Mrule,Ndef,jElt,GSize,Mk,point,rrule);
    }  break;
     case 2:  { 
       /* case 2 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       geometric non linear elastic 3D (see elem0.m for m file implementation)
       see implementation in of_mk_pre.c
      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */
       nonlin_elas(GF,ECp,GSize,jElt);       
     }  
     break;
     case 3: {
       /* elementary hyperelastic matrix computation now in hyper.c */
       Hyper_main(GF,ECp,Nnode,Ndof,Nw,Mdef,Ndef,CurDofPos);	  
     }
     break;
#include of_inc_fs
#ifdef   of_inc_user
#   include of_inc_user
#endif
     case -1: /* Early return */
	   break;

 case 999:   
       /* When introducing a new element family you need to 
        use standard variables : integ,constit,nodeEe,NDN
        to fill in
         ke : the element matrix (which will be assemble below)
        and possibly a current RHS as illustrated in the Hyperelastic
        strategy. You will then need to use the DofPositions of the current
        element
         CurDofPos=DofPos+jElt*DofPerElt;
        and add the elemental RHS to the current force in def.def(:,2) 
         F=def+Mdef; 
       */
       mexPrintf("jMat %i, jElt %i (strategy %i)",jMat,jElt,StrategyType);
       mexErrMsgTxt( "User element family place holder" );
       break; 
 default:
       mexPrintf("jMat %i, jElt %i (strategy %i)",jMat,jElt,StrategyType);
       mexErrMsgTxt( "Not a supported matrix strategy" );
    
 } /* swith on mat_og matrix assembly */

   if (GF.Ncondense) {
       /* xxx JM renumeroter la matrice, condenser les DDLs
          voir elem0.m ligne 402
          pour finaliser les inits cherche Ncondense ci-dessus
          Ne faire aucune alloc.
       */
	   of_ptrdiff   opt1[5];
	   double *ke_tmp, r1;

	   ke_tmp=(double*)malloc(Nk*Nk*sizeof(double));
       opt1[2]=1; opt1[3]=1; /* no increment */
        for (j1=0; j1<Nk; j1++) { for (j2=0; j2<Nk; j2++) {
        ke_tmp[j1+Nk*j2]=ECp[0].ke[ind[j1]+Nk*ind[j2]]; 
       }}  
       /* condense step by step */
       for (j1=0; j1<GF.Ncondense; j1++) {
         opt1[0]=Nk-j1-1; opt1[1]=Nk-j1-1; opt1[4]=Nk;
         r1=-1./ke_tmp[(Nk-j1-1)*Nk+Nk-j1-1];
         of_dger(&opt1[0], &opt1[1],  &r1,
                 &ke_tmp[(Nk-j1-1)*Nk], &opt1[2],
                 &ke_tmp[(Nk-j1-1)*Nk], &opt1[3], /* valid if  symmetric */ 
                 &ke_tmp[0],            &opt1[4]);
       }
       /*  Par rapport a ton condense_matrix remettre 
            une phase de renumerotation de ke_tmp -> ke  */
       j3=0;
       for (j1=0; j1<Nk-GF.Ncondense; j1++) { for (j2=0; j2<Nk-GF.Ncondense; j2++) {
        ECp[0].ke[j3]=ke_tmp[j1+Nk*j2];
        j3++;
       }} 
	   free(ke_tmp);
     }
    AssembleMatVec(ECp,GF,opt,StrategyType,Mk,CurDofPos,o_ir,o_jc,o_pr,
             elmap,DofPerElt,Mdef,Ndef,Ener,cEGI,jElt,Mener);

} /* end of jElt loop on elements - - - - - - - - - - - - - - - - -*/ 
  if (ECp[0].ke[0]!=100e100) memcpy(GF.ke,ECp[0].ke,ThreadiPos[5]*sizeof(double)); /* must be on a used thread */
  ofFree(ECv);
} /* end parallel region  last private gives ECe */
} /* strategy type */

 if (Ener==NULL) {
  if (nlhs>0 && jMat!=2) { /* return stiffness */
	 plhs[0]= mxCreateDoubleMatrix(Mk,Mk,mxREAL);
	 memcpy(mxGetPr(plhs[0]),GF.ke,ThreadiPos[5]*sizeof(double));
  } 
  if (nlhs>1 && jMat==2) { /* return mass */
	 plhs[1]= mxCreateDoubleMatrix(Mk,Mk,mxREAL);
	 memcpy(mxGetPr(plhs[1]),GF.ke,ThreadiPos[5]*sizeof(double));
  }
 }

} /* jMat<4 end multiple matrix case - - - - - - - - - - - - - - - - - */
ofFree(GF.ke);
if (nlhs>1&&plhs[1]==NULL)   plhs[1]= mxCreateDoubleMatrix(0,0,mxREAL);
/*strategy dependent clean up*/
#undef MatrixIntegrationStep 
#define MatrixIntegrationStep 3
switch (StrategyType) {
#include of_inc_fs
#ifdef   of_inc_user
#   include of_inc_user
#endif
}

NeedFreeAlloc(NeedFree,-6,-1);
EC.defE=NULL; 
    
/* ------------------------------------------------------------------------*/
} else if (!strcmp(CAM,"StressObserve"))  {

	/* see elem0 stress_og for m file equivalent
	*/

        mxArray       *field,*mEC,*mnode,*InfoAtNode; 
        double        *constit, coef,GfCtable4;
        int           Nw, Nwe, Nnode, NfieldE,  *point, *rule, nNDN, GSize[20],
                      jConst, jN, jw, Mrule, i1, i2, Ndof, Nstress;
        double*     NeedFree[6]={NULL,NULL,NULL,NULL,NULL,NULL};
        struct EltConst EC; 
        struct GroupFields GF;
        void*       pEC=NULL;

        /* inputs */
       mEC=(mxArray*)prhs[1]; mnode=(mxArray*)prhs[4]; InfoAtNode=NULL; 
       GSize[0]=(int)mxGetM(prhs[4]);Nnode = GSize[0];NfieldE=(int)mxGetN(prhs[4]);
       GF.N=NULL;GF=initGroup(GF,mnode,mEC,GSize,NeedFree);
       EC.nodeE=mxGetPr(prhs[4]); GSize[4]=(int)mxGetN(prhs[4])-4;/* nodeE assembled externally */
       EC=initInfoAtNode(InfoAtNode,mEC,EC,GF,GSize,NeedFree); 

        if (nrhs<6) { mexErrMsgTxt("6 input arguments needed."); }
        Nwe = (int)*mxGetPr(mxGetField(mEC,0,"Nw")); Nw=GSize[2];/* number of weight points, nodes*/
        nNDN        = (int)mxGetM (mxGetField(mEC, 0,"NDN"));
        

        if (!mxIsInt32(prhs[5])) { mexErrMsgTxt("point must be int32"); }
        point =  (int*)mxGetData(prhs[5]); 
        rule =  (int*)mxGetData(prhs[2]); 
        Mrule = (int)mxGetM(prhs[2]); /* size(rule,1) */
        Nstress=rule[Mrule*4];
        
        NDNSwitch(point[3],GF,&EC,Nw,Nnode,NfieldE);
        field=mxGetField(prhs[1], 0,"VectMap");
        if (field==NULL) mexErrMsgTxt("VectMap must be defined");
        Ndof=(int)(mxGetNumberOfElements(field)/Nnode); 

        constit = mxGetPr(prhs[3]);GSize[6]=(int)mxGetM(prhs[3]);
        pEC=FieldPointer(prhs[1],(char*)"pEC"); 
        if (pEC!=NULL) memcpy(pEC,&EC,40);/* sizeof(struct EltConst)); */

		EC.constit=constit;/* +point[6];*/
	   if (GF.CTable!=NULL) { 
	     double *constite;
         GfCtable4=GF.CTable[4];
		 if (GF.CTable[4]==-2) GF.CTable[4]=-2.1;
 	     constite=(double*)malloc(GSize[6]*sizeof(double));
	     memcpy(constite,EC.constit , GSize[6] * sizeof( double ) );
	     EC.constit0=EC.constit; EC.constit = constite;
   	    }
	    if (nrhs>6) EC.defE=mxGetPr(prhs[6]);

		{ /* prepare for // */
   int it, ThreadiPos[10];
   double* ECv; /* vector for alloc */
   struct EltConst* ECp;
   InitThreadiPos(ThreadiPos,GSize,5); /* xxx Mk */

   it=(sizeof(struct EltConst)/sizeof(double))+1;
   ECv=(double*)ofMalloc( (it + /* memory for ECe */
	   ThreadiPos[0]+ThreadiPos[1]+ThreadiPos[2]+ /* nodeE, NDN, jdet */
       GSize[6]+ThreadiPos[3]+ThreadiPos[4]+ThreadiPos[5]+GF.NBe+GF.NdefE)
	    *sizeof(double));
   ECp=(struct EltConst*)ECv;
   memcpy(ECp,&EC,sizeof(struct EltConst));
   /* missing shift of buffers for .nodeE ... threadipos should be in GF */
   if (point[4]==5) {
   /* first try at non linear stress evaluation - - - - - - - - - - - - - - -*/ 
        field = mxGetField(mEC,0,"material");
        if (field!=NULL) {
         mxGetString(field,CAM,15); 
         if (!strcmp("Elastic3DNL",CAM)) { /* geometric non linear elastic 3D */
	       int jElt=0; /* single element call for now */
		   plhs[0] = mxCreateDoubleMatrix(6*Nwe,mxGetN(prhs[6]),mxREAL); 
		   if (mxGetN(prhs[6])>1) mexPrintf("Multi shapes not implemented");
           EC.gstate = mxGetPr(plhs[0]); 
		   EC.ke=NULL;
           nonlin_elas(GF,&EC,GSize,jElt);
	 }
	}
	} else {
	 /* classical linear interation on all Defs of element */
		plhs[0] = mxCreateDoubleMatrix(Nstress*Nwe,Nnode*Ndof,mxREAL); 
        EC.ke = mxGetPr(plhs[0]); 
        for (jConst=0;jConst<Nstress*Nwe*Nnode*Ndof;jConst++) EC.ke[jConst]=0.;

        for (jw=0;jw<Nwe ;jw++) {
          if (GF.CTable!=NULL) { constitInterp(GF,ECp,Nnode*jw,Nnode,GSize);}
        for (jN=0;jN<Nnode ;jN++) {
        for (jConst=0;jConst<Mrule ;jConst++) {
          coef=ECp[0].constit[rule[jConst+2*Mrule]+point[6]];
		  if (rule[jConst+6*Mrule]<0) coef=-coef;
          i1=jw*rule[jConst+4*Mrule]+rule[jConst];
          i2=rule[jConst+3*Mrule]+jN*Ndof;
          EC.ke[i1+Nwe*Nstress*i2] += 
                EC.NDN[jN+ nNDN * (Nw*rule[jConst+Mrule]+jw+rule[jConst+5*Mrule]) ]*coef;     
             /*if (i1+Nw*Nstress*i2>mxGetNumberOfElements(plhs[0])) mexErrMsgTxt("Overflow");*/
        } /* jConst */
        } /* jN */
        } /* jw */
	}
	if (GF.CTable!=NULL) { GF.CTable[4]=GfCtable4;free(EC.constit); }
		}
 /* ------------------------------------------------------------------------*/
} else if (!strcmp(CAM,"StressCrit"))  {       
 /* compute equivalent stresses according to specific criteria
  * input is a 6xNe matrix with lines s1,s2,s3,s23,s31,s12
  *only VonMises available for the moment, as no need to compute principal stesses
  *sig_e=sqrt(.5*((r1(1,:)-r1(2,:)).^2+(r1(2,:)-r1(3,:)).^2+(r1(3,:)-r1(1,:)).^2+6*sum(r1(4:end,:).^2,1)));
	*/
  char     *crit;
  mwSize      Ne;
  int j1;
  double   *r1, *out;
          
   r1=mxGetPr(prhs[1]); crit=mxArrayToString(prhs[2]); Ne=mxGetNumberOfElements(prhs[1])/mxGetM(prhs[1]);
   plhs[0]=mxCreateDoubleMatrix(1,Ne,mxREAL);
   out = mxGetPr (plhs[0]);
   if (!strcmp("VonMises",crit))  {  
	#pragma omp parallel for 
    for (j1=0; j1<Ne; j1++) { 
     double r2;
     r2 = sqrt( .5*(  pow(r1[6*j1] - r1[6*j1+1],2)
             + pow(r1[6*j1+1] - r1[6*j1+2],2)
             + pow(r1[6*j1+2] - r1[6*j1],2)
             + 6*(pow(r1[6*j1+3],2) + pow(r1[6*j1+4],2) + pow(r1[6*j1+5],2)) ));
	 out[j1]=r2; 
    }
    
   } else { mexErrMsgTxt("of_mk StressCrit : Unknown criterium"); }     
 
/* ------------------------------------------------------------------------*/
	/*
addpath c:/balmes/sdt.cur/6.5
ofutil of_mk
k=zeros(10); x=[1 2 3]; of_mk('k<-k+a*x*y',k,x,x,1,int32([0 0 0 size(k,1)  3 3 1 1]));k
k=zeros(10); x=[1 2 3;4 5 6]'; 
of_mk('k<-k+a*x*y',k,x,x,1,int32([0 0 0 size(k,1)  3 3 1 1]));k

	 */
} else if (!strcmp(CAM,"k<-k+a*x*y"))  {

   double        *a, *x, *y, *k;
   of_ptrdiff     *opt;

   opt=(of_ptrdiff*)mxGetData(prhs[5]);
   k=mxGetPr(prhs[1])+opt[0];  x=mxGetPr(prhs[2])+opt[1];  
   y=mxGetPr(prhs[3])+opt[2];  a=mxGetPr(prhs[4]);
   /* mexPrintf("a=%g %i %i",a[0],opt[6],opt[7]); */
   /* opt=[offk[0] offx[1] offy[2] size(k,1)[3] m[4] n[5]  incx[6] incy[7]] */
   if (mxGetM(prhs[4])>1) {
     int ji,jj;
     for (ji=0;ji<mxGetM(prhs[4]);ji++) {for (jj=0;jj<mxGetN(prhs[4]);jj++) {
       if (a[0]!=0) {
        of_dger(&opt[4],&opt[5],a,
          x+opt[4]*ji,&opt[6],
          y+opt[5]*jj,&opt[7],k,&opt[3]);
       }
       a++;
     }}
   } else { of_dger(&opt[4],&opt[5],a,x,&opt[6],y,&opt[7],k,&opt[3]);}


/* ------------------------------------------------------------------------*/
} else if (!strcmp(CAM,"Mecha3DInteg"))  {

  int* opt;
  double *Be=NULL;
  opt=(int*)mxGetData(prhs[9]); 

 /* of_mk('Mecha3DInteg',k,Be,F,d2wde2,Sigma,w,jdet,NDN,int32([Nnode Ndof Nw jW])); */
  mexErrMsgTxt("not reimplemented missing EC");
 /* Mecha3DInteg(mxGetPr(prhs[1]),mxGetPr(prhs[2]),mxGetPr(prhs[3]),
              mxGetPr(prhs[4]),
              mxGetPr(prhs[5]),mxGetPr(prhs[6]),mxGetPr(prhs[7]),
              mxGetPr(prhs[8]),
              opt[0],opt[1],opt[2],opt[3]);*/

/*---------------------------------------------------------  */
/* insert a matrix in a sparse one                           */
/*---------------------------------------------------------  */
}  else if (!strcmp(CAM,"asmsparse"))  {

      /*
k1=k+0.;
      sp_util('asmsparse',k,int32(i1),ke,elmap);
      of_mk('asmsparse',k1,int32(i1),ke,elmap);
size(find(k-k1))

       */      
    double      *pr; 
    int         nVal, *elmap, 
                NDDL, IsSymK, IsSymVal, *opt; 
    mwIndex     *vir, *vjc, *ir, *jc;

    if (nrhs < 4 ) { mexErrMsgTxt( "At least 4 args required!" );  }

    ir = mxGetIr( prhs[1] ); jc = mxGetJc( prhs[1] ); pr = mxGetPr( prhs[1] ); 
    opt = (int*)mxGetData( prhs[2] ); 
    NDDL = opt[0]; IsSymVal = opt[1]; IsSymK = opt[2];  /* get options */
    nVal = (int)mxGetNumberOfElements( prhs[3] );       /* length of values */

    if ((IsSymVal && (nVal != (NDDL*(NDDL+1)/2)))) {
      mexPrintf("%i %i %i",nVal,NDDL,IsSymVal);
      mexErrMsgTxt( "Number of values does not match triangle!" );
    }
    if (((!IsSymVal) && (nVal != (NDDL*NDDL)))) {
      mexPrintf("%i %i %i",nVal,NDDL,IsSymVal);
      mexErrMsgTxt( "Number of values does not match full matrix!" );
    }

    vir=NULL; vjc=NULL;
    if (mxIsSparse(prhs[3])) {
      vir = mxGetIr( prhs[3] ); vjc = mxGetJc( prhs[3] ); 
    }    
 
    if (nrhs<5) {elmap=NULL;} else {elmap = (int*)mxGetData( prhs[4] );}
    AssembleSparse(ir,jc,pr,opt[0],opt[1],opt[2],opt+3,elmap,
		    vir,vjc,mxGetPr(prhs[3]));

    /* other input checking should be here !!! */

    /* Build output args : problem with scope should not be used */
    if (nlhs>0) plhs[0] = (mxArray*)prhs[1];


  /* fillval -----------------------------------------------------------
     set all values of a (sparse) matrix or a vector to a constant
     k=sparse(ones(10,10));of_mk('fillval',k,0);
     ------------------------------------------------------------------ */
}  else if (!strcmp(CAM,"fillvalue"))  {

    int nVal;
    double *pr, *pr0, val;

    pr0 = pr = mxGetPr( prhs[1] ); 
    if (mxIsSparse( prhs[1] )) { nVal = (int)mxGetNzmax( prhs[1] );
    } else if (mxIsDouble(prhs[1])) {nVal=(int)mxGetNumberOfElements(prhs[1]);
    } else {
      mexErrMsgTxt( "Unsupported input type!" );
    }
    val = mxGetScalar( prhs[2] ); /*      mexPrintf( "%f %d\n", val, nVal ); */

    while (pr < (pr0 + nVal)) {
      *pr = val; pr++;
    }
#ifndef MatlabVER
    plhs[0] = prhs[1];
#endif  
 
  /* asmFull ----------------------------------------------------------
     Dangerous - no bound checks!
     ------------------------------------------------------------------ */
}  else if (!strcmp(CAM,"asm_full"))  {

    int nVal, ii, ir, *pm;
    double *pr, *pv;

    pr = mxGetPr( prhs[1] ); 
    pm = (int*)mxGetData( prhs[2] ); 
    pv = mxGetPr( prhs[3] ); 
    
    nVal = (int)mxGetNumberOfElements( prhs[2] );
    for (ii = 0; ii < nVal; ii++) {
      ir = pm[ii];
      if (ir == -1) continue;
      pr[ir] += pv[ii];
    }

  /* insert ----------------------------------------------------------       
     Offsets in C sense!
     ------------------------------------------------------------------ */
}  else if (!strcmp(CAM,"insertval"))  {

    int offR = (int) (*mxGetPr( prhs[2] ));
    int offA = (int) (*mxGetPr( prhs[4] ));
    int nVal = (int) (*mxGetPr( prhs[5] ));
    double *pr, *pa;

/*      mexPrintf( "offR: %d, offA: %d, nVal: %d\n", offR, offA, nVal ); */

    pr = mxGetPr( prhs[1] ) + offR;
    pa = mxGetPr( prhs[3] ) + offA;
    
    if ((mxGetNumberOfElements( prhs[1] ) - offR) < nVal) {
      mexErrMsgTxt( "insert(): result vector too small!" );
    }
    if ((mxGetNumberOfElements( prhs[3] ) - offA) < nVal) {
      mexErrMsgTxt( "insert(): source vector too small!" );
    }
    memcpy( pr, pa, nVal * sizeof( double ) );

  /* mesh_meshGraph ------------------------------------------------------ 
     ------------------------------------------------------------------ */
}  else if (!strcmp(CAM,"meshGraph"))  {

    mxArray *out;
    int nNod, nGr, mConn, ii;
    mwIndex *conn, *ptrGr, *nEl, *nEP, *jc, *icol, *aux;
    mwSize nnz;
    double *val;

    nrhs--; prhs++; /* !!! */
    nNod = (int) mxGetScalar( prhs[0] );
    nGr = (int) mxGetScalar( prhs[1] );
    mConn = (int) mxGetScalar( prhs[2] );
    ptrGr = (mwIndex *) mxGetData( prhs[3] ); /* group pointer */
    nEl = (mwIndex *) mxGetData(prhs[4] );
    nEP = (mwIndex *) mxGetData(prhs[5] );
    conn = (mwIndex *) mxGetData( prhs[6] );

    mesh_meshGraph( &nnz, &jc, &icol, nNod, nGr, mConn,
		    ptrGr, nEl, nEP, conn );

    out = mxCreateSparse( nNod, nNod, nnz, mxREAL );

    val = mxGetPr(out);
    for (ii = 0; ii < nnz; ii++) { val[ii] = 1.0; }
    aux = mxGetIr( out ); if (aux) mxFree( aux );
    mxSetIr(out,icol);

    aux = mxGetJc( out ); if (aux) mxFree( aux );
    mxSetJc(out,jc);

    plhs[0] = out;

}  else if (!strcmp(CAM,"dofGraph"))  {

    mxArray *out;
    int32 nNod, nEq, *dpn, *dofOffset, *eq;
    int32 ii;
    mwIndex *dprow, *dicol, *nprow, *nicol, *aux;
    double *val; 
    mwSize dnnz;

    nrhs--; prhs++; /* !!! */
    nNod = (int)mxGetM( prhs[0] );
    nprow = mxGetJc( prhs[0] );
    nicol = mxGetIr( prhs[0] );

    nEq = (int32) mxGetScalar( prhs[1] );
    dpn = (int32 *) mxGetPr( prhs[2] );
    dofOffset = (int32 *) mxGetPr( prhs[3] );
    eq = (int32 *) mxGetPr( prhs[4] );

    mesh_dofGraph( &dnnz, &dprow, &dicol,
		   nNod, nprow, nicol, nEq, dpn, dofOffset, eq );

    out = mxCreateSparse( nEq, nEq, dnnz, mxREAL );

    val = mxGetPr( out );
    for (ii = 0; ii < dnnz; ii++) {
      val[ii] = 1.0;
    }
    aux = mxGetIr( out );
    if (aux) mxFree( aux );
    mxSetIr( out, dicol );

    aux = mxGetJc( out );
    if (aux) mxFree( aux );
    mxSetJc( out, dprow );

    plhs[0] = out;
} else if (!strcmp(CAM,"pecgetcall"))  {
    if (1) { /*(nrhs==1) {*/
      plhs[0] =  mxCreateDoubleMatrix(1,1,mxREAL);
      out = mxGetPr(plhs[0]);
      out[0]=((double)sizeof(struct EltConst))/4;
      return;
    }
}  else if (!strcmp(CAM,"translam"))  {
    
    plhs[0] =  mxCreateDoubleMatrix(6,6,mxREAL);
    TransformLambda(mxGetPr(prhs[1]),mxGetPr(prhs[2]),mxGetPr(plhs[0]));
    
/*end of commands */
/* ------------------------------------------------------------------------*/
}  else {mexPrintf("%s\n",CAM);mexErrMsgTxt("Not a valid of_mk call");}

/* -----------------------------------------------------------------------*/
/* calls with commands with <8 characters MODULEF family of elements      */
/* -----------------------------------------------------------------------*/
} else {

    
  int *point, *integ;
  double *constit, *node, *estate,  *InfoAtNode, *defe, *eltconst;
  
  if (!mxIsInt32(prhs[2])) {mexErrMsgTxt("integ must be int32");}
  if (!mxIsInt32(prhs[1])) {mexErrMsgTxt("point must be int32");}
  
  point = (int*)mxGetData(prhs[1]); 
  integ = (int*)mxGetData(prhs[2]);
  constit = mxGetPr(prhs[3]);
  node = mxGetPr(prhs[4]);

  if (nrhs>5)  estate = mxGetPr(prhs[5]);     else         estate = node;
  if (nrhs>6)  defe = mxGetPr(prhs[6]);       else         defe = node;
  if (nrhs>7)  eltconst = mxGetPr(prhs[7]);   else         eltconst = node;
  if (nrhs>8)  InfoAtNode = mxGetPr(prhs[8]); else         InfoAtNode = node;

  integ[point[5]+4]=point[4]; /* output type in integ */
  typ=point[4];

  if (point[0]==0 ) { /*     DEAL WITH DEFAULT OUTPUT SIZES */
     if (typ == 0) {
	    point[0] = integ[point[5] +2]*(integ[point[5]+2]+1)/2;
	    point[1] = point[0];
    } else if (typ == 1 || typ == 2 || typ == 3) { /* symmetric one matrix only*/
	    point[0] = integ[point[5] +2]*(integ[point[5]+2]+1)/2;
    } else if (typ == 100) {/*  RHS is given*/
	    point[0] = integ[point[5] +2];
    } else { mexErrMsgTxt("Assembly option not supported"); }
  }

  plhs[0] = mxCreateDoubleMatrix(point[0],1,mxREAL); out = mxGetPr(plhs[0]);
  
  if (nlhs>1)  {
      plhs[1] = mxCreateDoubleMatrix(point[1],1,mxREAL);out1=mxGetPr(plhs[1]);
  } else out1 = out;

  out2 = out;
  
  ierr[0]=' '; ierr[1]=0;
  io_switch(CAM,ierr,point,integ,constit,node,estate,defe,eltconst,out,
      out1,out2,InfoAtNode);

  if (ierr[0]!=' ') mexPrintf("%s\n",ierr);

}}


