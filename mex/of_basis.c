#include <math.h>

void basis_clean(double* z) { 
  int j1;
  for (j1=0;j1<9;j1++) {
     if (fabs(z[j1])<1e-13) z[j1]=0;
  }
}
/* ------------------------------------------------------------------*/
/*            calculate basis z, input are 2 vectors x and y         */
/*            z is a 3x3 matrix                                      */
/* ------------------------------------------------------------------*/
/* ------------------------------------------------------------------*/
void basis (const double* x, const double* y, double* z) { 

  double r1, r2, eps=1.e-16, r3;
  int        j1, j2, j3;

  r1=0.; for (j1=0;j1<3;j1++) r1+=x[j1]*x[j1];       /* x=x/sqrt(x'*x); */
  r3=r1;
  if (r1<eps) { z[0]=1.; z[1]=0.; z[2]=0.; r3=1.; }
  else {r1=sqrt(r1); z[0]=x[0]/r1; z[1]=x[1]/r1; z[2]=x[2]/r1;}
  r1=0.; for (j2=0;j2<3;j2++) r1+=z[j2]*y[j2];       /* y=y-(x'*y)*x; */
  for (j1=0;j1<3;j1++) z[j1+3]=y[j1]-r1*z[j1];

  r2=0.; for (j1=3;j1<6;j1++) r2+=z[j1]*z[j1]; /* r2=norm(y) */
  if (r2<eps*r3) {
    j1=0; j3=0;  r1=fabs(z[0]);r2=fabs(z[0]);
    while (j1<3) {                                   /* [y,i1]=min(abs(x)); */
      if (r2<r1) { r1=r2; j3=j1; }
      j1++;r2=fabs(z[j1]);
    } /* determine j3 */

    for (j1=3;j1<6;j1++) z[j1]=0.;                   /* y=zeros(3,1);    */
    z[j3+3]=1.;                                      /* y(i1(1))=1;      */

    r1=0.; for (j2=0;j2<3;j2++) r1+=z[j2]*z[j2+3];   /* y=y-(x'*y)*x;    */
    for (j1=0;j1<3;j1++)  z[j1+3]-=r1*z[j1];
    r2=0.; for (j1=3;j1<6;j1++) r2+=z[j1]*z[j1]; 

  } /* if */


  r2=sqrt(r2);   
  z[3]=z[3]/r2;z[4]=z[4]/r2;z[5]=z[5]/r2; /* y=y/sqrt(y'*y); */
  /* z=-[y(2)*x(3)-y(3)*x(2);y(3)*x(1)-y(1)*x(3);y(1)*x(2)-y(2)*x(1)]; */
  z[6]= -z[4]*z[2]+z[5]*z[1];
  z[7]= -z[5]*z[0]+z[3]*z[2];
  z[8]= -z[3]*z[1]+z[4]*z[0];
  basis_clean(z);
}
