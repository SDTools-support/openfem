#include "../mex/hyper.h"
#include "../mex/of_def.h"
#include "../mex/of_EltConst.h"

extern void Mecha3DIntegH(struct EltConst* ECp, double *F, double *F1, double* d2wde2, 
          double* d2wvde2, double* Sigma,
	      double* w, int Nnode,int Ndof,int Nw, int jw);

/*U_ij, F_ij,C,I1,I2,I3, et derivees--------------------------------------------------*/
void elemcalc(struct EltConst* ECp,int jw, int Nw, int Nnode, double *e0, double *e, double *e1, double *F_ij, double *F1,
	     double *I, double *dI2dc, double *dI3dc, double *d2I3dcdc, int Mdef)
{   
  double       C[9],U_ij0[9], U_ij[9], U_ij1[9],*def0,zero=0.,deux=2.,neg = -1.;
  int          j1,j2,j3,ji,jj,jk; 
  of_ptrdiff    neuf=9,trois=3,nul=0,unit=1; 
  
  def0=ECp[0].defE+2*Mdef;
  /* U */
  for (ji=0;ji<9;ji++)  {U_ij[ji]=0; U_ij0[ji]=0; U_ij1[ji]=0;}
  for (ji=0;ji<3;ji++){ for (jj=0;jj<3;jj++) {
      for (jk=0;jk<Nnode;jk++) {
 /* mexPrintf("U :  %g %g(%i)\n",ECp[0].defE[ji+3*jk],ECp[0].NDN[Nnode*jw*4+4*jk+jj+1],Nnode*jw*4+4*jk+jj+1);*/
	U_ij[ji+3*jj]+=ECp[0].defE[ji+3*jk]*ECp[0].NDN[Nnode*jw*4+4*jk+jj+1];
	U_ij0[ji+3*jj]+=def0[ji+3*jk]*ECp[0].NDN[Nnode*jw*4+4*jk+jj+1];
  }}}
  of_dcopy(&neuf,&zero,&nul,U_ij1,&unit);
  of_daxpy(&neuf,&deux,U_ij,&unit,U_ij1,&unit);
  of_daxpy(&neuf,&neg,U_ij0,&unit,U_ij1,&unit);
 /* E */

  for (j1=0;j1<3;j1++) {
    e0[j1]=U_ij0[4*j1] + 0.5*of_ddot(&trois,U_ij0+3*j1,&unit,U_ij0+3*j1,&unit);
    e[j1]=U_ij[4*j1] + 0.5*of_ddot(&trois,U_ij+3*j1,&unit,U_ij+3*j1,&unit);
    e1[j1]=U_ij1[4*j1] + 0.5*of_ddot(&trois,U_ij1+3*j1,&unit,U_ij1+3*j1,&unit);
    for (j2=j1+1;j2<3;j2++) {
      e0[6-(j1+j2)]=(U_ij0[j1+3*j2]+U_ij0[j2+3*j1] + of_ddot(&trois,U_ij0+3*j1,&unit,U_ij0+3*j2,&unit));
      e[6-(j1+j2)]=(U_ij[j1+3*j2]+U_ij[j2+3*j1] + of_ddot(&trois,U_ij+3*j1,&unit,U_ij+3*j2,&unit));
      e1[6-(j1+j2)]=(U_ij1[j1+3*j2]+U_ij1[j2+3*j1] + of_ddot(&trois,U_ij1+3*j1,&unit,U_ij1+3*j2,&unit));
    }
  }




  /* F */
  for (ji=0;ji<9;ji++)  {   
    F_ij[ji]=U_ij[ji];
    F1[ji]=U_ij1[ji];
  }
  F_ij[0]+=1;F_ij[4]+=1;F_ij[8]+=1;
  F1[0]+=1;F1[4]+=1;F1[8]+=1;
  


  /* C */
  for (j2=0;j2<9;j2++) C[j2]=0.;
  for (j1=0;j1<3;j1++) { 
    for (j2=0;j2<3;j2++) {
      for (j3=0;j3<3;j3++)  C[j1+3*j2]+=F_ij[j3+3*j1]*F_ij[j3+3*j2]; 
    }
  }
  
  /* I1 =tr(C)*/ 
  I[0]=C[0]+C[4]+C[8];
 
  /* I2 =1/2 (tr(C)^2-tr(C^2)) = tr(Cof(C))*/
  I[1]= C[4]*C[8]-C[5]*C[7]+C[0]*C[8]-C[2]*C[6]+C[0]*C[4]-C[1]*C[3];

  /* I3 = det(C)*/
  I[2] = C[0]*( C[4]*C[8]-C[5]*C[7] )
    +C[1]*( C[5]*C[6]-C[3]*C[8] )
    +C[2]*( C[3]*C[7]-C[4]*C[6] );

 
  /* d^2(I3)/dc^2 */
  for (j2=0;j2<36;j2++)  d2I3dcdc[j2]=0.;
  d2I3dcdc[1] =C[8];    d2I3dcdc[2] =C[4];      d2I3dcdc[3] =-C[7];
  d2I3dcdc[6] =C[8];    d2I3dcdc[8] =C[0];      d2I3dcdc[10]=-C[6];
  d2I3dcdc[12]=C[4];    d2I3dcdc[13]=C[0];      d2I3dcdc[17]=-C[3];
  d2I3dcdc[18]=-C[7];   d2I3dcdc[21]=-C[0]/2.;  d2I3dcdc[22]=C[1]/2.;
  d2I3dcdc[23]=C[2]/2.; d2I3dcdc[25]=-C[6];     d2I3dcdc[27]=C[3]/2.;
  d2I3dcdc[28]=-C[4]/2; d2I3dcdc[29]=C[5]/2.;   d2I3dcdc[32]=-C[3];
  d2I3dcdc[33]=C[6]/2.; d2I3dcdc[34]=C[7]/2.;   d2I3dcdc[35]=-C[8]/2.;
  
  /* dI2/dc = i1 d_ij-c_ij*/
  for (j1=0;j1<3;j1++)       dI2dc[j1]=I[0]-C[j1+3*j1];
      dI2dc[3] =-C[7]; 
      dI2dc[4] =-C[6]; 
      dI2dc[5] =-C[3]; 

  /* dI3/dc =i3*cij^-1=cof(C)*/

  dI3dc[0]= C[4]*C[8]-C[5]*C[7];
  dI3dc[1]= C[0]*C[8]-C[2]*C[6];
  dI3dc[2]= C[0]*C[4]-C[1]*C[3];
  dI3dc[3]= C[1]*C[6]-C[0]*C[7];
  dI3dc[4]= C[3]*C[7]-C[4]*C[6]; 
  dI3dc[5]= C[5]*C[6]-C[3]*C[8]; 
  
}



/*Energie branche passive---------------------------------------------------------*/
void EnPassiv(int *integ,double *constit,double *I,double *dWdI,double *d2WdI2)
{
  
  constit ++;
  
  if ((int) constit[0] == 0) {
    /*************************************************
     *energie du type C1*(J1-3)+C2*(J2-3)+K*(J3-1)^2
     *************************************************/
    dWdI[0] = constit[1]*pow(I[2],-1./3.);
    dWdI[1] = constit[2]*pow(I[2],-2./3.);
    dWdI[2] = -1./3.*constit[1]*I[0]*pow(I[2],-4./3.) 
      -2./3.*constit[2]*I[1]*pow(I[2],-5./3.) 
      + constit[3]*(1.-pow(I[2],-.5));
    d2WdI2[0]=0.; d2WdI2[3]=0.; d2WdI2[6]=-1./3.*constit[1]*pow(I[2],-4./3.)*4; 
    d2WdI2[1]=0.; d2WdI2[4]=0.; d2WdI2[7]=-2./3.*constit[2]*pow(I[2],-5./3.)*4;
    d2WdI2[2]= -1./3.*constit[1]*pow(I[2],-4./3.)*4;
    d2WdI2[5]= -2./3.*constit[2]*pow(I[2],-5./3.)*4;
    d2WdI2[8]=  4./9.*constit[1]*I[0]*pow(I[2],-7./3.)*4
      + 10./9.*constit[2]*I[1]*pow(I[2],-8./3.)*4
      + .5*constit[3]*pow(I[2],-3./2.)*4;
  }
  
  else if ((int) constit[0] == 1) {
    /*****************************************************************
    *energie du type C1*(J1-3)+C2*(J2-3)+K*(J3-1)-K*ln(J3)
    ******************************************************************/

    dWdI[0] = constit[1]*pow(I[2],-1./3.);
    dWdI[1] = constit[2]*pow(I[2],-2./3.);
    dWdI[2] = -1./3.*constit[1]*I[0]*pow(I[2],-4./3.) 
      -2./3.*constit[2]*I[1]*pow(I[2],-5./3.)
      +0.5*constit[3]*pow(I[2],-0.5)
      -0.5*constit[3]*pow(I[2],-1);

    d2WdI2[0]=0.; d2WdI2[3]=0.; d2WdI2[6]=-1./3.*constit[1]*pow(I[2],-4./3.)*4; 
    d2WdI2[1]=0.; d2WdI2[4]=0.; d2WdI2[7]=-2./3.*constit[2]*pow(I[2],-5./3.)*4;
    d2WdI2[2]= -1./3.*constit[1]*pow(I[2],-4./3.)*4;
    d2WdI2[5]= -2./3.*constit[2]*pow(I[2],-5./3.)*4;
    d2WdI2[8]=  4./9.*constit[1]*I[0]*pow(I[2],-7./3.)*4
      + 10./9.*constit[2]*I[1]*pow(I[2],-8./3.)*4
      -1./4.*constit[3]*pow(I[2],-3./2.)*4
      +0.5*constit[3]*pow(I[2],-2)*4;
  } 
}


/*calcul d'un element la matrice tangente au point d'integration-------------------------*/
void matPassiv( double *Sigma, double *dWdI, double *d2WdI2,
	       double *dI2dc, double *dI3dc, double *d2I3dcdc,double *d2wde2)
{
  int j1;
  of_ptrdiff  opt1[5];
        
    for (j1=0;j1<6;j1++) { /* compute Sigma */
         Sigma[j1]+=2.*( dWdI[0]*dI1dc[j1] + dWdI[1]*dI2dc[j1] + dWdI[2]*dI3dc[j1] );
       }
	/*  Equation (16) (second order derivative with respect to e */

       opt1[0]=6;opt1[1]=6;opt1[2]=1;opt1[3]=1;opt1[4]=6; /* for blas call */
       for (j1=0;j1<36;j1++) d2wde2[j1] += 4*(dWdI[1]*d2I2dcdc[j1] + dWdI[2]*d2I3dcdc[j1]);


       of_dger(&opt1[0],&opt1[1],&d2WdI2[0],dI1dc,&opt1[2],dI1dc,&opt1[3],d2wde2,&opt1[4]);
       of_dger(&opt1[0],&opt1[1],&d2WdI2[1],dI2dc,&opt1[2],dI1dc,&opt1[3],d2wde2,&opt1[4]);
       of_dger(&opt1[0],&opt1[1],&d2WdI2[2],dI3dc,&opt1[2],dI1dc,&opt1[3],d2wde2,&opt1[4]);
       of_dger(&opt1[0],&opt1[1],&d2WdI2[3],dI1dc,&opt1[2],dI2dc,&opt1[3],d2wde2,&opt1[4]);
       of_dger(&opt1[0],&opt1[1],&d2WdI2[4],dI2dc,&opt1[2],dI2dc,&opt1[3],d2wde2,&opt1[4]);
       of_dger(&opt1[0],&opt1[1],&d2WdI2[5],dI3dc,&opt1[2],dI2dc,&opt1[3],d2wde2,&opt1[4]);
       of_dger(&opt1[0],&opt1[1],&d2WdI2[6],dI1dc,&opt1[2],dI3dc,&opt1[3],d2wde2,&opt1[4]);
       of_dger(&opt1[0],&opt1[1],&d2WdI2[7],dI2dc,&opt1[2],dI3dc,&opt1[3],d2wde2,&opt1[4]);
       of_dger(&opt1[0],&opt1[1],&d2WdI2[8],dI3dc,&opt1[2],dI3dc,&opt1[3],d2wde2,&opt1[4]);



}

/* viscous term for dynamic case*/
void EnVisco(double *constit, double *e0, double *e, double *Sigma, double *d2wdep2)
{
  int j1;

  for (j1=0;j1<6;j1++) { 
    Sigma[j1]+=constit[0]*(e[j1]-e0[j1])/constit[1];
  }
  /*compute d2wdep2 = eta.I36 */
  for (j1=0;j1<6;j1++) {
    d2wdep2[j1+6*j1]+= constit[0]/constit[1];
  }


}



/* hyperelastic_main */

void Hyper_main(struct GroupFields GF,struct EltConst* ECp, int Nnode,int Ndof,int Nw,int Mdef, int Ndef,int *CurDofPos)
{
  int j1,j2,jw;
  double e0[6], e[6], e1[6], F_ij[9], F1[9], Sigma[6], d2wde2[36], d2wvde2[36];
  double d2I3dcdc[36], dI2dc[6], dI3dc[6], I[3], dWdI[3], d2WdI2[9];


  for (jw=0;jw<Nw;jw++) {/* loop on integ. points --------------------- */
    
    for (j1=0;j1<6;j1++) { 
      for (j2=0;j2<6;j2++) { d2wde2[j1+6*j2]=0; d2wvde2[j1+6*j2]=0;}	   
      Sigma[j1]=0;
    }

    elemcalc(ECp,jw,Nw,Nnode,e0,e,e1,F_ij,F1,I,dI2dc,dI3dc,d2I3dcdc,Ndof);
    EnPassiv(ECp[0].integ,ECp[0].constit,I,dWdI,d2WdI2);
    matPassiv(Sigma,dWdI,d2WdI2,dI2dc,dI3dc,d2I3dcdc,d2wde2);
    EnVisco(ECp[0].constit+5, e0, e1, Sigma, d2wvde2);
    Mecha3DIntegH(ECp,F_ij,F1,d2wde2,d2wvde2,Sigma,GF.w,Nnode,Ndof,Nw,jw);	

  }/*end loop on integ. points ------*/


}
