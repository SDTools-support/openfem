/*
 Etienne Balmes, Jean Michel Leclere, Guillaume Vermot des Roches
 $Revision: 1.15 $  $Date: 2019/02/05 14:00:55 $
*/
/*----------------------------------------------------------- LinInterp 
  of_time_LinInterp(table,val,cur,out,Mtable,Ntable,Nval,nlhs,tablei,outi)
*/

#include "mex.h" 
#include "../mex/of_def.h"

#include "stddef.h" /* 4 NULL */

OF_EXPORT void of_time_LinInterp(double* table, double * val, double* last,
		       double* out, int M, int Nc, int Npoints, int nlhs,
               double* ti, double* outi) {

     int i1, j1, j2, j3, flag, ist;
     double s;

   #if MatlabVER >= 904
     if (ti==NULL) {ist=(int)1;} else {ist=(int)2; }
   #else
     ist=(int)1;
   #endif
  for (j2=0; j2<Npoints; j2++) { /* loop on points */
    /* search */
   i1=(int)(last[0]); if (i1<0){ i1=0;}
   flag=0; /* mexPrintf("M=%i %.5g , %.5g \n",M,table[i1+1],val[j2]); */
   while (1) {
    /* mexPrintf("i1=%i\n",i1);*/
      if (i1>=ist*(M-1)) { /* in the last interval */
        i1=ist*(M-1); 
        if ((i1>0)&&(val[j2]<table[i1-ist])) i1-=ist; 
        else break; 
      } else if (val[j2]>table[i1+ist]) {i1+=ist; /* need to move right*/
      } else if  (i1<=0) { /* in the first interval */
          if (val[j2]>table[ist]) i1=ist; else {i1=0;break;} 
      } else if (val[j2]<table[i1]) {i1-=ist;}  /* need to move left*/
      else break;
   }
   last[0]=i1; last[1]=val[0];
   /* mexPrintf("i1=%i\n",i1); */
   if (i1==ist*(M-1)&&i1>0) {i1-=ist;} /*If i1 is last table value move one backward to get slope for extrap*/
   if (i1<0) i1=0;
   /* mexPrintf("i1=%i\n",i1);*/
   if (Nc==0&&nlhs) out[j2]=i1+ist;
   /* linear interpolation */
   /* outputs */
   if (nlhs==0 || M==0) { 
	    s=(table[i1+ist]-table[i1]);
     if (M==1) {s=1;} else if (s==0) {s=.5;} else { s=(val[j2]-table[i1])/s;}
     last[2]=s;
   } else {
      /* if two equal X move interval to that x */
      if (i1<ist*(M-2) && val[j2]==table[i1+ist]&&table[i1+2*ist]==table[i1+ist]) {i1+=ist;}
      else if (i1>ist-1 && val[j2]==table[i1] && table[i1-ist]==table[i1]) {i1-=ist;}
      /*mexPrintf("M=%i %.5g , %.5g %.5g\n",M,table[i1+1],table[i1],val[j2]);*/
      s=(table[i1+ist]-table[i1]);
      if (M==1) {s=0;} /* single point, use first value */
	     else if (s==0) {s=.5;} 
      else { s=(val[j2]-table[i1])/s;}
      last[2]=s;
      /*mexPrintf("M=%i %.5g , %.5g %.5g  s=%.5g\n",M,table[i1+ist],table[i1],val[j2],s);*/
      #if MatlabVER >= 904
       for (j1=0; j1<Nc; j1++) { /* loop on curves */
        for (j3=0; j3<ist; j3++) { /* loop on interleaved */
out[ist*(j1*Npoints+j2)+j3] = (1-s)*table[i1+ist*(j1+1)*M+j3] + s*table[i1+ist+ist*(j1+1)*M+j3];
        }
       } 
      #else
       for (j1=0; j1<Nc; j1++) { /* loop on curves */
       out[j1*Npoints+j2] = (1-s)*table[i1+(j1+1)*M] + s*table[i1+1+(j1+1)*M];
       } 
       if (ti!=NULL && outi!=NULL) {
        for (j1=0; j1<Nc; j1++) { /* loop on curves */
         outi[j1*Npoints+j2] = (1-s)*ti[i1+(j1+1)*M] + s*ti[i1+1+(j1+1)*M];
        } 
       }
      #endif
   } /* trivial or not */

  } /* j2 points to interpolate */
}
