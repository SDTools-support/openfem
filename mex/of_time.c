#include "mex.h"
#include "of_def.h"
#include "of_time_interp.c"
#include "pmat.c"

#include "math.h"
#include <string.h>

/* #include "../../sdtdev/mex50/sdt_getFieldByNumber.c"*/
/*--------------------------------------------------------------------*/
/* LININTERP   */
/*--------------------------------------------------------------------*/

/* the mex function -----------------------------------------------*/
void mexFunction (int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]) {
char           *buf ;


  if (nrhs==0)  {
  double *a;
  plhs[0] = mxCreateDoubleMatrix (1, 1, mxREAL);
  a = mxGetPr (plhs[0]);   *a=5.0001;
  return; 
  }

  if (nrhs>1 && mxIsStruct(prhs[0])) { buf=mxGetFieldNameByNumber(prhs[0],0);
  } else if (nrhs>1 && !mxIsChar(prhs[0])) {
  /* setinput for double full vectors of_time([-1 offset],target,input) */
    
  mwSize len,offset,i1;
  double* command;
  
  if (*mxGetPr(prhs[0])==-1) { /* do usual setinput */
   mxArray    *field;

   if (mxGetN(prhs[0])>1) { command=mxGetPr(prhs[0]);offset=(mwSize)command[1]; }
   else { offset=0; }
   len=mxGetNumberOfElements(prhs[2]);
   if (mxIsDouble(prhs[1])) {field=prhs[1]; i1=mxGetNumberOfElements(field);
   } else { /* scalar for now */
    double* vect;
	vect=(double*)pmatGetData(prhs[1]);
    if (len>1) mexErrMsgIdAndTxt("SDT:SetInput","pMat only implemented for scalars");
	memcpy(vect,pmatGetData(prhs[2]),1*sizeof(double)); 
 if (nlhs>0) {
  double *a; /* xxx need to cleanup this duplicated part  */
       plhs[0] = mxCreateDoubleMatrix (1, 1, mxREAL);
       a = mxGetPr (plhs[0]);   *a=(double)(offset+len);
   } 
	return;
    /* mexPrintf("%i %p %g add\n",i1,vect,vect[0]); */
   }
   if (i1==1) mexErrMsgIdAndTxt("SDT:SetInput","SetInput only works on vectors");
   if (i1<len+offset||offset<0) mexErrMsgIdAndTxt("OpenFEM:of_time","overflow"); /* avoid diff with unsigned */
   memcpy(mxGetPr(field)+offset,pmatGetData(prhs[2]),len*sizeof(double)); 
   if (nlhs>0) {
       double *a;
       plhs[0] = mxCreateDoubleMatrix (1, 1, mxREAL);
       a = mxGetPr (plhs[0]);   *a=(double)(offset+len);
   } 
  } else if (*mxGetPr(prhs[0])==-2) { /* allocate memory and return address */
    double *vect, **pvect;
	vect=(double*)malloc(sizeof(double));vect[0]=-111; 

	pvect=&vect;
	/* mexPrintf("pInit  %p %p",vect,pvect[0]);*/
	plhs[0] = mxCreateDoubleMatrix (1, 1, mxREAL);
	memcpy(mxGetData(plhs[0]),pvect,sizeof(double*)); 

	/* pvect=(double**)mxGetData(plhs[0]); vect=pvect[0]; mexPrintf("  pend  %p \n ",vect); */

  } else if (*mxGetPr(prhs[0])==-2.1) { /* copy data to matlab */
    double* vect;
	vect=(double*)pmatGetData(prhs[1]);
	/* mexPrintf("okb  %p = %g\n",vect,vect[0]);*/
	plhs[0] = mxCreateDoubleMatrix (1, 1, mxREAL); 
	memcpy(mxGetData(plhs[0]),vect,sizeof(double));

  } else if (*mxGetPr(prhs[0])==-2.2) { /* free scalar memory */
    double* vect;
	vect=(double*)pmatGetData(prhs[1]);
    /* mexPrintf("ok %p",vect);mexEvalString("dbstack;");*/
    free(vect);

  } else if (*mxGetPr(prhs[0])==-2.3) { /* return isReal double 0 or 1 */
   double *a, *vect;
   plhs[0] = mxCreateDoubleMatrix (1, 1, mxREAL);
   a = mxGetPr(plhs[0]);
   vect=(double*)pmatGetData(prhs[1]);
   if (mxIsComplex(prhs[1])) {*a=(double)(0);} else; {*a=(double)(1);}
   
  } else { mexErrMsgIdAndTxt("OpenFEM:of_time","Not a valid command");
  }
  return;
} else buf = mxArrayToString(prhs[0]); /* xxx mxGetString for fixed buf */

/*----------------------------------------------------------- LinInterp */
if (!strcmp("lininterp",buf))  {

  /* out = of_time ('lininterp',table,val,last)
   if table is a row, constant value is used
   */

  double    *out, *outi, *table, *val, *last, *tablei;
  mwSize       Nc, M, Npoints;
  bool ti; 

/*
table=[   1     1     2; 2     2     4; 3     3     6];
val=2.5;last=[2 1 .5]; % index starts at 0
out = of_time ('lininterp',table,val,last);last
of_time('lininterp',[0 0;1 1],linspace(0,1,10)',[1 0 0])
out=of_time('lininterp',[0;1;2],linspace(0,2,10)',[0 0 0])
*/

  if (nrhs!=4)  mexErrMsgIdAndTxt("OpenFEM:of_time","4 inputs required.");
  if (!mxIsDouble(prhs[2]))  mexErrMsgIdAndTxt("OpenFEM:of_time","val must be double");

  /* inputs */
  M     = mxGetM(prhs[1]);   /* number of points */
  Nc    = mxGetN(prhs[1])-1; /* number of curves */
 #if MatlabVER >= 904
  /*table = mxGetPr(prhs[1]); ti =NULL; /* [xi yi ...] */
  ti=mxIsComplex(prhs[1]);tablei=NULL; 
  if (ti) table = (double*)mxGetComplexDoubles(prhs[1]);
  else table = mxGetDoubles(prhs[1]);
 #else
  ti=false; 
  table = mxGetPr(prhs[1]); tablei = mxGetPi(prhs[1]); /* [xi yi ...] */
 #endif 
  if ((mxGetM(prhs[2])==0) || (mxGetN(prhs[2])==0)) {
   plhs[0] = mxCreateDoubleMatrix (0, 0, mxREAL); return;
  }
  
  val   = mxGetPr(prhs[2]);  /*  xs         */
  last  = mxGetPr(prhs[3]);  /* [i1 xi si  ] */
  Npoints= mxGetM(prhs[2]);
  if (mxGetN(prhs[2])>Npoints) Npoints=mxGetN(prhs[2]);
  
  if (nlhs>0) {
   mwSize i1=Nc;
   if (Nc<1) i1=1;
    #if MatlabVER >= 904
    if (!ti) { 
     plhs[0] = mxCreateDoubleMatrix (Npoints,i1, mxREAL);
     out=(double*)mxGetDoubles(plhs[0]);outi=NULL; 
    }
    else {
     plhs[0] = mxCreateDoubleMatrix (Npoints,i1, mxCOMPLEX);
     out=(double*)mxGetComplexDoubles(plhs[0]);outi=NULL;
    }
    #else
   if (tablei!=NULL) plhs[0] = mxCreateDoubleMatrix (Npoints,i1, mxCOMPLEX);
   else plhs[0] = mxCreateDoubleMatrix (Npoints,i1, mxREAL);
   out= mxGetPr(plhs[0]);outi=mxGetPi(plhs[0]);
    #endif 
  } else {out=NULL;}
  of_time_LinInterp(table,val,last,out,M,Nc,Npoints,nlhs,ti,tablei,outi);
}
/*----------------------------------------------------------- storelaststep */
else if (!strcmp("storelaststep",buf))  {

#ifdef SDT_ADD
#define AddStep 101
#include "../../sdtdev/mex50/OpenFEMAdd.c"
#else
  /*
   of_time('storelaststep',uva,u,v,a);
  */
double     *uva, *u, *v, *a;
int        N;
 N=mxGetM(prhs[1]); 
 if (mxIsStruct(prhs[1])) {uva=mxGetPr(mxGetField(prhs[1],0,"uva"));} 
 else uva=mxGetPr(prhs[1]);
 u=mxGetPr(prhs[2]); memcpy(uva,u,N*sizeof(double)); 
 v=mxGetPr(prhs[3]); memcpy(uva+N,v,N*sizeof(double)); 
 a=mxGetPr(prhs[4]); memcpy(uva+2*N,a,N*sizeof(double)); 
#endif

/*----------------------------------------------------------- interp */
#ifdef SDT_ADD
#pragma warning( disable : 4005 )
#define AddStep 103
#include "../../sdtdev/mex50/OpenFEMAdd.c"
#endif

/*----------------------------------------------------------- newmarkinterp */
} else if (!strcmp("newmarkinterp",buf))  {
  /* of_time ('newmarkinterp', out, beta,gamma,uva,a1, t0,t1,other)
     %Geradin p.371 equation 7.3.9   */
double     beta, gamma, *uva, *a1,t0, t1, *def=NULL, *v=NULL, *a=NULL, *cur, *tout, dt, *other;
int         j1, j2, M, *ind, Nind, Ntout, Ndef, Nother;
mxArray    *field;

/* inputs */
if (nrhs<8)  mexErrMsgIdAndTxt("OpenFEM:of_time","8 inputs required.");
beta = *mxGetPr(prhs[2]); gamma = *mxGetPr(prhs[3]);
uva = mxGetPr(prhs[4]); a1 = mxGetPr(prhs[5]);
t0 = *mxGetPr(prhs[6]);  t1 = *mxGetPr(prhs[7]); 

/* out structure (input and output) */
field=mxGetField(prhs[1],0,"cur");
if (field==NULL) mexErrMsgIdAndTxt("OpenFEM:of_time","def.cur field required");
cur = mxGetPr (field);    /* [last_ind  last_t] */
field=mxGetField(prhs[1],0,"data");
if (field==NULL) mexErrMsgIdAndTxt("OpenFEM:of_time","def.data field required");
tout = mxGetPr(field);  Ntout=(int)mxGetM(field);

field=mxGetField(prhs[1],0,"OutInd");
if (field==NULL) mexErrMsgIdAndTxt("OpenFEM:of_time","out.OutInd must be defined");
else if (!mxIsInt32(field)) mexErrMsgIdAndTxt("OpenFEM:of_time","out.OutInd must be int32");
ind = mxGetData(field); Nind=(int)mxGetN(field);   /* length of ind  */
if ((int)mxGetM(field)>Nind) {Nind=(int)mxGetM(field);}

field=mxGetField(prhs[1],0,"def");
if (field!=NULL) {def = mxGetPr (field);Ndef=(int)mxGetM(field);} else def=NULL;
if (mxGetField(prhs[1],0,"v")!=NULL)   v   = mxGetPr (mxGetField(prhs[1],0,"v"));  
if (mxGetField(prhs[1],0,"a")!=NULL)   a   = mxGetPr (mxGetField(prhs[1],0,"a"));  
if (nrhs==9) {other=mxGetPr(prhs[8]);Nother=(int)mxGetM(prhs[8]);} else other=NULL;


M=(int)mxGetM(prhs[4]);   /* Nddl        */
j1=(int)cur[0];  
if (j1>Ntout-1) {j1=0;} else if (tout[j1]>t1) {j1=0;}
 while(j1<Ntout) {
  if (tout[j1]>=t0 && tout[j1]<=t1) {/* the current step contains a SaveTime */
  cur[0]=(double)(j1+1); cur[1]=tout[j1];   dt=tout[j1]-t0;
  /* mexPrintf(" j1=%i [%g : %g : %g ]\n",j1,t0,tout[j1],t1); */

  for (j2=0; j2<Nind; j2++) { 
   if (def!=NULL) def[j2+Ndef*j1]= uva[ind[j2]-1] + dt*uva[ind[j2]-1+M] 
                   + dt*dt*(.5-beta)*uva[ind[j2]-1+2*M] + dt*dt*beta*a1[ind[j2]-1];
   if (v!=NULL)  v[j2+Nind*j1]=uva[ind[j2]-1+M]+(1-gamma)*dt*uva[ind[j2]-1+2*M]
                                                           +gamma*dt*a1[ind[j2]-1];
   if (a!=NULL)  a[j2+Nind*j1]= a1[ind[j2]-1];
  
  }
#ifdef SDT_ADD
#define AddStep 102
#include "../../sdtdev/mex50/OpenFEMAdd.c"
#else
  /* if (def!=NULL&&other!=NULL) {
      for (j2=0; j2<Nother; j2++) {def[j2+Nind+Ndef*j1]=other[j2];}
  } */
  mexPrintf("Other FNL DOF not supported in OpenFEM (see SDT)\n");

#endif
  } /* current step contains save time */
 if (tout[j1]>t1)  {  break; }
 j1++;
}


} else if (!strcmp("cinterp",buf))  {
/*----------------------------------------------------------- TableInterp 
of_time('cinterp',nodeE,EC.NDN,constit,EltConst.CTable,double(jW))
*/

double     *nodeE, *NDN, *constit, *r1, *r2,*X,val[1];
int         j1, j2, i2, jw, M,n1, Nnode;

nodeE=mxGetPr(prhs[1]); Nnode=(int)mxGetM(prhs[1]); 
r1=mxGetPr(prhs[5]);jw=(int)r1[0];
NDN=mxGetPr(prhs[2])+Nnode*jw;constit=mxGetPr(prhs[3]);r1=mxGetPr(prhs[4]);
n1=(int)r1[0];
for (j1=0;j1<n1;j1++) {
  r2=r1+j1*7+1;M=(int)r2[4];X=r1+(int)r2[3]; i2=Nnode*((int)r2[5]-1);
    val[0]=0;for (j2=0;j2<Nnode;j2++) {val[0]+=NDN[j2]*nodeE[j2+i2];}
    if (M<0) constit[(int)r2[6]-1]=val[0]; /* use FieldAtNode Value */
    else of_time_LinInterp(X,val,r2,constit+(int)r2[6]-1,M,1,1,1,false,NULL,NULL);
}
} else if (!strcmp("issdt",buf))  {
#ifdef SDT_ADD
  plhs[0]=mxCreateDoubleScalar(1);
#else
  plhs[0]=mxCreateDoubleScalar(0);
#endif
/*----------------------------------------------------------- v_mkl 
} else if (!strcmp("v_mat",buf))  {
 
 mxArray *lhs[1];
 mxArray *mx, *PropertyPtr = NULL;
 mxArray *lhs[1];

 double *r1;
 if (nrhs<3) mexErrMsgTxt("3 args needed");
 r1=mxGetPr(prhs[2]);
 mx = mexCallMATLABWithTrap(1,lhs,2,prhs+1,"subsref");
 mxPrintf("ok %p ",mx);
 //plhs[0]=sdt_getFieldByNumber(mx,0,(int)r1[0]);

 */
  
} else if (!strcmp("cvs",buf))  {
/*----------------------------------------------------------- cvs */
mxArray    *field, *field2;
mwSize     dims[2];

field=mxCreateString("$Revision: 1.62 $  $Date: 2023/08/21 10:37:44 $");
#ifdef SDT_ADD
#define AddStep 104
#include "../../sdtdev/mex50/OpenFEMAdd.c"
    dims[0]=2; dims[1]=2; plhs[0]=mxCreateCellArray((mwSize)2,dims);
    mxSetCell(plhs[0],0,mxCreateString("of_time"));
    mxSetCell(plhs[0],1,mxCreateString("OpenFEMAdd"));
    mxSetCell(plhs[0],2,field);
    mxSetCell(plhs[0],3,field2);
#else
  plhs[0]=field;
#endif

}
/*----------------------------------------------------------- end command */
if (!mxIsStruct(prhs[0])) mxFree(buf);
} /* end mexFunction */    
