
#include "mex.h"
#include "of_def.h"
#include "math.h"
#include <string.h>

/* $Revision: 1.4 $  $Date: 2018/06/12 15:44:34 $ */
        
void* pmatGetData(const mxArray *cfield) {
	if (mxIsDouble(cfield)) { 
    void* v1;
	 v1=mxGetData(cfield);
 	 return v1; 
	} else {
     mxArray    *field;
   	 double *vect, **pvect;
 	 field=mxGetProperty(cfield,0,"c");
	 pvect=(double**)mxGetData(field); vect=pvect[0];
 	 return vect; 
	}
}

void pmatGetPrAndSize(const mxArray *cfield, double** pr, mwSize *i1) {
	if (mxIsDouble(cfield)) { 
	 pr[0]=(double*)mxGetData(cfield);
     i1[0]=mxGetNumberOfElements(cfield);
     if (i1[0]==1) mexErrMsgIdAndTxt("SDT:SetInput","SetInput only works on vectors");
	} else if ( mxIsClass(cfield,"pmat")) {
     mxArray    *field;
   	 double **pvect;
 	 field=mxGetProperty(cfield,0,"c");
     if (field!=NULL) {
	   i1[0]=1;
        pvect=(double**)mxGetData(field); pr[0]=(double*)pvect[0];
	 } 
	} else {
	 pr[0]=(double*)mxGetData(cfield);
     i1[0]=mxGetNumberOfElements(cfield);
	}

}
