/* ------------------------------------------------------------------*/
/*            calculate X'*K*X                                       */
/*            X is a 3x3 block                                       */
/* ------------------------------------------------------------------*/

#include "mex.h"
#include <string.h>

void x_k_x (const double* x, const double* y, double* result, 
  double* offset, int type, int YM,int YN) { 

  int             j1, j2, jOff;   
  double          *K1, *cof;
    
  K1 = (double*)mxCalloc (YM*YN, sizeof(double));/* memory allocation for z1 */

  if (type!=2) { /* Right multiplication */
    jOff=0;cof=offset; j2=0;while (j2+2<YN) { /* loop on columns */
    for (j1=0;j1<YM;j1++) {
       K1[j1+j2*YM]=
         y[j1+j2*YM]*x[0]+y[j1+(j2+1)*YM]*x[1]+y[j1+(j2+2)*YM]*x[2];
       K1[j1+(j2+1)*YM]=
         y[j1+j2*YM]*x[3]+y[j1+(j2+1)*YM]*x[4]+y[j1+(j2+2)*YM]*x[5];
       K1[j1+(j2+2)*YM]=
         y[j1+j2*YM]*x[6]+y[j1+(j2+1)*YM]*x[7]+y[j1+(j2+2)*YM]*x[8];
    }
    j2+=3;jOff+=3;
    if (jOff==6 && offset!=NULL) {/*account for rigid link defined by offset*/
     for (j1=0;j1<YM;j1++) {
       K1[j1+(j2-3)*YM]-=cof[2]*K1[j1+(j2-5)*YM];
       K1[j1+(j2-3)*YM]+=cof[1]*K1[j1+(j2-4)*YM];
       K1[j1+(j2-2)*YM]+=cof[2]*K1[j1+(j2-6)*YM];
       K1[j1+(j2-2)*YM]-=cof[0]*K1[j1+(j2-4)*YM];
       K1[j1+(j2-1)*YM]-=cof[1]*K1[j1+(j2-6)*YM];
       K1[j1+(j2-1)*YM]+=cof[0]*K1[j1+(j2-5)*YM];
     }
      jOff=0;cof+=3;
    } /*offset*/
    }
  } else {memcpy(K1,y,YM*YN*sizeof(double));} /* right multiplication */

  if (type!=1) { /* Left multiplication */
    jOff=0;cof=offset; j1=0;while (j1+2<YM) { /* loop on columns */
    for (j2=0;j2<YN;j2++) {
       result[j1+j2*YM]=
         K1[j1+j2*YM]*x[0]+K1[j1+1+j2*YM]*x[1]+K1[j1+2+j2*YM]*x[2];
       result[j1+1+j2*YM]=
         K1[j1+j2*YM]*x[3]+K1[j1+1+j2*YM]*x[4]+K1[j1+2+j2*YM]*x[5];
       result[j1+2+j2*YM]=
         K1[j1+j2*YM]*x[6]+K1[j1+1+j2*YM]*x[7]+K1[j1+2+j2*YM]*x[8];
    }
    j1+=3;jOff+=3;
    if (jOff==6 && offset!=NULL) {/*account for rigid link defined by offset*/
     for (j2=0;j2<YN;j2++) {
       result[j1-3+j2*YM]-=cof[2]*result[j1-5+j2*YM];
       result[j1-3+j2*YM]+=cof[1]*result[j1-4+j2*YM];
       result[j1-2+j2*YM]+=cof[2]*result[j1-6+j2*YM];
       result[j1-2+j2*YM]-=cof[0]*result[j1-4+j2*YM];
       result[j1-1+j2*YM]-=cof[1]*result[j1-6+j2*YM];
       result[j1-1+j2*YM]+=cof[0]*result[j1-5+j2*YM];
     }
      jOff=0;cof+=3;
    } /*offset*/
    }
  } else {memcpy(result,K1,YM*YN*sizeof(double));} /* left multiplication */

  mxFree (K1);                                    /* desallocate k1 */
}
/* ------------------------------------------------------------------*/
