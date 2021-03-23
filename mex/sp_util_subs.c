
#include <math.h>

/*  Etienne Balmes, Jean Michel Leclere, Robert Cimrman */
/*  Claire Delforge                                     */

#ifndef int32
#define int32 int
#endif

mxArray* subs_cvs_sputil () {
 mxArray *st;
 st= mxCreateString("$Revision: 1.15 $  $Date: 2007/03/29 15:30:38 $");
 return(st);
}


#include "of_xkx.c"
#include "of_basis.c"


/* ------------------------------------------------------------------*/
/*            basis for T3                                           */
/*                                                                   */
/* ------------------------------------------------------------------*/
void basisT3 (const double* x, const int Nelt, int* Elt, double* x1, 
   int nPerElt ) { 

 double     r1, r2, node[9], eps=1.e-16;
 double     *y1, *z1;
 int        jElt, jNode, j1, j2, j3;

 y1=x1+3; z1=x1+6;
 
 for (jElt=0;jElt<Nelt;jElt++) {

   for (jNode=0;jNode<3;jNode++) {  /* node => tmp */
     for (j1=0;j1<3;j1++) node[3*jNode+j1] = x[ 3*(Elt[nPerElt*jElt+jNode]-1)+j1 ];
   }

   r1=0.; for (j1=0;j1<3;j1++) r1+=node[j1]*node[j1]; 
   if (r1<eps) { x1[9*jElt+0]=1.; x1[9*jElt+1]=0.; x1[9*jElt+2]=0.;  }
   else { /* x=x/sqrt(x'*x); */
    r1=sqrt(r1); for (j1=0;j1<3;j1++) { x1[9*jElt+j1]=node[j1]/r1; } 
   }

   r1=0.; for (j2=0;j2<3;j2++) r1+=node[j2]*node[3+j2];  /* y=y-(x'*y)*x; */
   for (j1=0;j1<3;j1++) y1[9*jElt+j1]=node[3+j1]-r1*x1[9*jElt+j1];
   r2=0.; for (j1=0;j1<3;j1++) r2+=y1[9*jElt+j1]*y1[9*jElt+j1]; /* r2=norm(y) */

   if (r2<eps) {
     j1=0; j3=0;  r1=fabs(z1[9*jElt]);
     while (j1<3) {                                   /* [y,i1]=min(abs(x)); */
       j1++; r2=fabs(x1[9*jElt+j1]);
       if (r2<r1) { r1=r2; j3=j1; }
     } /* determine j3 */
     for (j1=0;j1<3;j1++) y1[9*jElt+j1]=0.;            /* y=zeros(3,1);    */
     y1[9*jElt+j3]=1.;                                 /* y(i1(1))=1;      */
     r1=0.; for (j2=0;j2<3;j2++) r1+=node[j2]*node[3+j2];  /* y=y-(x'*y)*x; */
     for (j1=0;j1<3;j1++) y1[9*jElt+j1]=node[3+j1]-r1*x1[9*jElt+j1];
     r2=0.; for (j1=0;j1<3;j1++) r2+=y1[9*jElt+j1]*y1[9*jElt+j1]; /* r2=norm(y) */
   } /* if */
   r2=sqrt(r2);   
   for (j1=0;j1<3;j1++) { y1[9*jElt+j1]=y1[9*jElt+j1]/r2; }/* y=y/sqrt(y'*y); */
   /* z=-[y(2)*x(3)-y(3)*x(2);y(3)*x(1)-y(1)*x(3);y(1)*x(2)-y(2)*x(1)]; */
   z1[9*jElt+0]= -x1[9*jElt+4]*x1[9*jElt+2] + x1[9*jElt+5]*x1[9*jElt+1];
   z1[9*jElt+1]= -x1[9*jElt+5]*x1[9*jElt+0] + x1[9*jElt+3]*x1[9*jElt+2];
   z1[9*jElt+2]= -x1[9*jElt+3]*x1[9*jElt+1] + x1[9*jElt+4]*x1[9*jElt+0];
 } /* jElt */
}
/* ------------------------------------------------------------------*/
/*            basis for Q4                                           */
/*                                                                   */
/* ------------------------------------------------------------------*/
void basisQ4 (const double* x, const int Nelt, int32 *Elt, double *x1,
              double *X, double *out3, int nPerElt  ) { 

 double     r1, node[12];
 double     *y1, *z1, out2[3];
 int        jElt, jNode, j1, j2, j3;

 y1=x1+3; z1=x1+6;
 
 for (jElt=0;jElt<Nelt;jElt++) {

   for (jNode=0;jNode<4;jNode++) {  /* node => tmp */
     for (j1=0;j1<3;j1++) node[3*jNode+j1] = x[ 3*(Elt[nPerElt*jElt+jNode]-1)+j1 ];
   }
   for (j1=0;j1<3;j1++) { /* out2 = mean(x); */
     out2[j1]=0.; for (j2=0;j2<4;j2++) out2[j1]+=node[3*j2+j1]/4.;
   }
   for (j1=0;j1<3;j1++) { /* node = x(1:4,:)-out2(ones(4,1),:); */
     r1=out2[j1]; 
     for (j2=0;j2<4;j2++)  node[3*j2+j1]=x[3*(Elt[nPerElt*jElt+j2]-1)+j1]-r1;
   }
   /* x = node(3,:)-node(1,:); x=x/sqrt(x*x'); */
   for (j1=0;j1<3;j1++) { x1[9*jElt+j1]=node[6+j1]-node[j1];  }
   r1=0.; for (j1=0;j1<3;j1++) r1+=x1[9*jElt+j1]*x1[9*jElt+j1]; 
   for (j1=0;j1<3;j1++) { x1[9*jElt+j1]/=sqrt(r1); }
   /* y = node(4,:)-node(2,:); y=y/sqrt(y*y'); */
   for (j1=0;j1<3;j1++) { y1[9*jElt+j1]=node[9+j1]-node[3+j1]; }
   r1=0.; for (j1=0;j1<3;j1++) r1+=y1[9*jElt+j1]*y1[9*jElt+j1]; 
   for (j1=0;j1<3;j1++) { y1[9*jElt+j1]/=sqrt(r1); }
   /* x = x-y; x=x'/sqrt(x*x'); */
   for (j1=0;j1<3;j1++) x1[9*jElt+j1]-=y1[9*jElt+j1];
   r1=0.; for (j1=0;j1<3;j1++) r1+=x1[9*jElt+j1]*x1[9*jElt+j1]; 
   for (j1=0;j1<3;j1++) { x1[9*jElt+j1]/=sqrt(r1); }
   /* y=y';y=y-(x'*y)*x; y=y/sqrt(y'*y); */
   r1=0.; for (j1=0;j1<3;j1++) r1+=x1[9*jElt+j1]*y1[9*jElt+j1];
   for (j1=0;j1<3;j1++) y1[9*jElt+j1]-=r1*x1[9*jElt+j1];
   r1=0.; for (j1=0;j1<3;j1++) r1+=y1[9*jElt+j1]*y1[9*jElt+j1]; 
   for (j1=0;j1<3;j1++) { y1[9*jElt+j1]/=sqrt(r1);  }
   /* z=-[y(2)*x(3)-y(3)*x(2);  y(3)*x(1)-y(1)*x(3);  y(1)*x(2)-y(2)*x(1)]; */
   z1[9*jElt+0] = -y1[9*jElt+1]*x1[9*jElt+2] + y1[9*jElt+2]*x1[9*jElt+1];
   z1[9*jElt+1] = -y1[9*jElt+2]*x1[9*jElt+0] + y1[9*jElt+0]*x1[9*jElt+2];
   z1[9*jElt+2] = -y1[9*jElt+0]*x1[9*jElt+1] + y1[9*jElt+1]*x1[9*jElt+0];

 } /* jElt */

  if (X!=NULL) { /* X=node*p; */
   for (j1=0;j1<12;j1++) X[j1]=0.;

   for (j2=0;j2<4;j2++) {   
    for (j3=0;j3<3;j3++) {
      X[j2]   += node[3*j2+j3]*x1[j3]; 
      X[4+j2] += node[3*j2+j3]*y1[j3]; 
      X[8+j2] += node[3*j2+j3]*z1[j3];    
      /*      X[j2]   += node[j2+4*j3]*x1[j3]; 
      X[4+j2] += node[j2+4*j3]*y1[j3]; 
      X[8+j2] += node[j2+4*j3]*z1[j3];    */
    }   
   } 

  } /* X */
  if (out3!=NULL) {  for (j1=0;j1<3;j1++)  out3[j1]=out2[j1];   }


} /* end basisQ4 */
