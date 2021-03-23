#include <stdio.h>
#include <string.h>
#ifdef OFOMP
#include <omp.h>
#endif
#include "mex.h"
#include "math.h"
#include "matrix.h"
#include "of_def.h"
#define mwSize int
#define mwIndex int
#define int32 int
#include "of_mk_pre.c"


void mexFunction (int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]) 
{
  double *out, *temp, *in;
  int i,j,k, m, n;
  char CAM[15];
  

  if (nrhs == 0)
    {
      plhs[0] =  mxCreateDoubleMatrix(1,1,mxREAL);
      out = mxGetPr(plhs[0]);
      out[0]=666;
      return;
    }

  mxGetString(prhs[0],CAM,15);

  if ((nrhs > 1) && (!strcmp("test",CAM))) {

    mxArray     *field;
    double   *ndnE, *jdetE, *w, *ke, *constit, *nodeE, *node, *N, *Nr, *Ns, *Nt, 
      *o_pr, *def, *defi, *defelem, *Be, *RHS, *bas, *J;
    int jcol, jw, j1, j2, j3, j4, Mk, Nelt, i1,i2, Nw, Nshape, *o_ir, *o_jc, *opt, *elmap, *DofPos, DofPerElt,nb;
    int lndn, jElt, *NodePos, Nnode, Ndof, Nmnode,Nfield,*point,*integ, Mdef, Ndef, opt1[5], *CurDofPos;
    static double dI1dc[6]={1.,1.,1.,0.,0.,0.};
    static double d2I2dcdc[36]={0.0,1.0,1.0,0.0,0.0,0.0,1.0,0.0,1.0,0.0,
				0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.5,0.0,0.0,
				0.0,0.0,0.0,0.0,-0.5,0.0,0.0,0.0,0.0,0.0,0.0,-0.5};   



    NodePos=mxGetData(prhs[2]); 
    if (!mxIsInt32(prhs[2])) mexErrMsgTxt("NodePos must be int32");

        
    Nnode=mxGetM(prhs[2]); 
    Nmnode=mxGetM(prhs[3]); node=mxGetPr(prhs[3]);
    Ndof = Nnode*3;

    if (!mxIsStruct(prhs[10])) mexErrMsgTxt("EltConst must be a structure");
    if (mxIsEmpty(prhs[10])) mexErrMsgTxt("EltConst must not be empty");

    field = mxGetField(prhs[10], 0,"NDN");
    if (field!=NULL) {
      field   = mxGetField(prhs[10], 0,"w");
      w  = mxGetPr (field); Nw=mxGetM(field); w = w+3*Nw;
      Nshape=mxGetN(mxGetField(prhs[10], 0,"Nr"));
     }
    N=mxGetPr(mxGetField(prhs[10], 0,"N"));
    Nr=mxGetPr(mxGetField(prhs[10], 0,"Nr"));
    Ns=mxGetPr(mxGetField(prhs[10], 0,"Ns"));
    
    field = mxGetField(prhs[10], 0,"Nt");
    if (field!=NULL) Nt=mxGetPr(field); else Nt=NULL;    
    if (field!=NULL) Nt=mxGetPr(field); else Nt=NULL;   
    field=mxGetField(prhs[10], 0,"bas"); 
    if (field==NULL||mxGetM(field)!=9||mxGetN(field)!=Nw) bas=NULL; 
    else  bas = mxGetPr (field);
    field=mxGetField(prhs[10], 0,"J"); 
    if (field==NULL||mxGetM(field)!=4||mxGetN(field)!=Nw) J=NULL; 
    else  J = mxGetPr (field);


    Nelt=mxGetN(prhs[1]); jElt=0;
    o_pr=mxGetPr(prhs[12]); o_ir=mxGetIr(prhs[12]); o_jc=mxGetJc(prhs[12]); 
    opt=mxGetData(prhs[13]); DofPerElt=opt[0]; 
    elmap=mxGetData(prhs[8]);DofPos=mxGetData(prhs[1]);
    def=mxGetPr(prhs[11]);defi=mxGetPi(prhs[11]);
    Mdef=mxGetM(prhs[11]);Ndef=mxGetN(prhs[11]);

    point   = (int*)mxGetData(prhs[4]);
    integ   = (int*)mxGetData(prhs[5])+point[5];
    Mk=integ[2]; 
     
    plhs[0] =  mxCreateDoubleMatrix(Mk,Mk,mxREAL);
    
    Mk=mxGetM(plhs[0]);
    nb = Mk*Mk;
    lndn = 4*Nw*Nnode;

#ifdef TESTALLOCFAIL
      nodeE=calloc(Nnode*4,sizeof(double));
      ke=calloc(nb,sizeof(double));
      Be=calloc(3*Nnode,sizeof(double));
      defelem=calloc(3*Nnode*Ndef,sizeof(double));
      ndnE = calloc(lndn,sizeof(double));
      jdetE = calloc(Nw,sizeof(double));
#endif



#pragma omp parallel for private(jElt,nodeE,ke,Be,defelem,point,constit,i1,ndnE,jdetE,j1,j2,jw,opt1)

    for (jElt=0;jElt<Nelt;jElt++) { /* loop on elements to assemble */

      int *CurDofPos;
      double F_ij[9], Sigma[6], d2wde2[36],zero=0,un=1;
      double d2I3dcdc[36], dI2dc[6], dI3dc[6], I[3], dWdI[3], d2WdI2[9];
      int ji,jk,jj,jl,jm,jn,nul=0,unit=1,deux=2,trois=3,quatre=4,six=6,neuf=9;
      double C[9], U_ij[9];
      double coef[1], DD[81] ,D[54], AUX[54], TEMP[81], *B, alpha, beta,zcoef;
      char norm = 'N', trans= 'T';

      CurDofPos=DofPos+jElt*DofPerElt;
      point   = (int*)mxGetData(prhs[4])+jElt*mxGetM(prhs[4]);
      point[4]= *((int*)mxGetData(prhs[4])+4);

      constit =  mxGetPr(prhs[6])+point[6];
#ifdef TESTALLOCWORKS
      nodeE=calloc(Nnode*4,sizeof(double));
      ke=calloc(nb,sizeof(double));
      Be=calloc(3*Nnode,sizeof(double));
      defelem=calloc(3*Nnode*Ndef,sizeof(double));
      ndnE = calloc(lndn,sizeof(double));
      jdetE = calloc(Nw,sizeof(double));
#endif

      for (j2=0;j2<Nnode;j2++) {/* nodeE nodes of current element*/
	i1=NodePos[j2+Nnode*jElt]-1;  
	nodeE[j2]=node[i1+Nmnode*4];   
	nodeE[j2+Nnode]=node[i1+Nmnode*5]; 
	nodeE[j2+2*Nnode]=node[i1+Nmnode*6]; 
      }

      NDNSwitch(point[3],w,N,Nr,Ns,Nt,nodeE,jdetE,ndnE,Nw,Nnode,Nfield,Nshape,bas,J);

      of_dcopy(&nb,&zero,&nul,ke,&unit);
      of_dcopy(&Mk,&zero,&nul,Be,&unit);

      for (j1=0;j1<Mk;j1++) { for (j2=0;j2<Ndef;j2++) {
	  defelem[j1+Mk*j2]=def[CurDofPos[j1]+Mdef*j2];
	}}
      for (jw=0;jw<Nw;jw++) {/* loop on integ. points --------------------- */
	
	for (j1=0;j1<6;j1++) { 
	  for (j2=0;j2<6;j2++) {
	    d2wde2[j1+6*j2]=0;
	  }	   
	  Sigma[j1]=0;
	}
	
	
	
	/* U */
	for (ji=0;ji<9;ji++)  {U_ij[ji]=0;}
	for (ji=0;ji<3;ji++){ for (jj=0;jj<3;jj++) {for (jk=0;jk<Nnode;jk++) {
		  U_ij[ji+3*jj]+=defelem[ji+3*jk]*ndnE[Nnode*jw*4+4*jk+jj+1];
	}}}

	/* F */
	for (ji=0;ji<9;ji++)  {   F_ij[ji]=U_ij[ji];}
	F_ij[0]+=1;F_ij[4]+=1;F_ij[8]+=1;
  
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




        
	/*EnPassiv(integ,constit,I,dWdI,d2WdI2);*/


	constit ++;
	
	if (constit[0] == 0) {
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
	
	else if (constit[0] == 1) {
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



	/*matPassiv(opt1,Sigma,dWdI,d2WdI2,dI2dc,dI3dc,d2I3dcdc,d2wde2);*/


	for (ji=0;ji<6;ji++) { /* compute Sigma */
	  Sigma[ji]+=2.*( dWdI[0]*dI1dc[ji] + dWdI[1]*dI2dc[ji] + dWdI[2]*dI3dc[ji] );
	}
	/*  Equation (16) (second order derivative with respect to e */

	opt1[0]=6;opt1[1]=6;opt1[2]=1;opt1[3]=1;opt1[4]=6; /* for blas call */
	for (ji=0;ji<36;ji++) d2wde2[ji] += 4*(dWdI[1]*d2I2dcdc[ji] + dWdI[2]*d2I3dcdc[ji]);

	of_dger(&opt1[0],&opt1[1],&d2WdI2[0],dI1dc,&opt1[2],dI1dc,&opt1[3],d2wde2,&opt1[4]);
	of_dger(&opt1[0],&opt1[1],&d2WdI2[1],dI2dc,&opt1[2],dI1dc,&opt1[3],d2wde2,&opt1[4]);
	of_dger(&opt1[0],&opt1[1],&d2WdI2[2],dI3dc,&opt1[2],dI1dc,&opt1[3],d2wde2,&opt1[4]);
	of_dger(&opt1[0],&opt1[1],&d2WdI2[3],dI1dc,&opt1[2],dI2dc,&opt1[3],d2wde2,&opt1[4]);
	of_dger(&opt1[0],&opt1[1],&d2WdI2[4],dI2dc,&opt1[2],dI2dc,&opt1[3],d2wde2,&opt1[4]);
	of_dger(&opt1[0],&opt1[1],&d2WdI2[5],dI3dc,&opt1[2],dI2dc,&opt1[3],d2wde2,&opt1[4]);
	of_dger(&opt1[0],&opt1[1],&d2WdI2[6],dI1dc,&opt1[2],dI3dc,&opt1[3],d2wde2,&opt1[4]);
	of_dger(&opt1[0],&opt1[1],&d2WdI2[7],dI2dc,&opt1[2],dI3dc,&opt1[3],d2wde2,&opt1[4]);
	of_dger(&opt1[0],&opt1[1],&d2WdI2[8],dI3dc,&opt1[2],dI3dc,&opt1[3],d2wde2,&opt1[4]);



		  
	/*Mecha3DInteg*/

	B=calloc(18*Nnode,sizeof(double));


	for (ji=0;ji<6;ji++) {for (jj=0;jj<9;jj++) D[9*ji+jj]=0;}

	D[0]=F_ij[0];  D[3]=F_ij[1];  D[6]=F_ij[2]; 
	D[10]=F_ij[3];  D[13]=F_ij[4];  D[16]=F_ij[5]; 
	D[20]=F_ij[6];  D[23]=F_ij[7];  D[26]=F_ij[8];
	
	D[28]=F_ij[6];  D[29]=F_ij[3];  D[31]=F_ij[7]; 
	D[32]=F_ij[4];  D[34]=F_ij[8];  D[35]=F_ij[5];

	D[36]=F_ij[6];  D[38]=F_ij[0];  D[39]=F_ij[7]; 
	D[41]=F_ij[1];  D[42]=F_ij[8];  D[44]=F_ij[2];

	D[45]=F_ij[3];  D[46]=F_ij[0];  D[48]=F_ij[4]; 
	D[49]=F_ij[1];  D[51]=F_ij[5];  D[52]=F_ij[2];

	/* aux=d2wde2 *d */
	for (ji=0;ji<6;ji++) {
	  for (jj=0;jj<9;jj++) {
	    
	    AUX[9*ji+jj] =0;
	    
	    for (jn=0;jn<6;jn++) {
	      AUX[9*ji+jj] +=  d2wde2[6*jn+ji]*D[9*jn+jj];
	    }
	  }
	}
	
	/* dd=d^t aux */

	for (ji=0;ji<9;ji++) {
	  for (jj=0;jj<9;jj++) {

	    DD[9*ji+jj] =0;

	    for (jn=0;jn<6;jn++) {
	      DD[9*ji+jj] +=  D[9*jn+ji]*AUX[9*jn+jj];
	    }}}

	/* rajouter k_nl */
	for (jj=0;jj<3;jj++) {
	  DD[9*jj+jj]     += Sigma[jj];
	  DD[9*(jj+3)+jj+3] += Sigma[jj];
	  DD[9*(jj+6)+jj+6] += Sigma[jj];
	  for (ji=0;ji<3;ji++) {
	    if(ji!=jj) {
	      DD[9*ji+jj]     += Sigma[6-jj-ji];
	      DD[9*(ji+3)+jj+3] += Sigma[6-jj-ji];
	      DD[9*(ji+6)+jj+6] += Sigma[6-jj-ji];
	    }
	  }
	}



	/*for (ji=0;ji<9;ji++) { for (jj=0;jj<9;jj++) { 
	  mexPrintf("%10.5g",DD[ji+9*jj]);} mexPrintf("\n");};mexPrintf("\n");*/
  
  
	alpha = jdetE[jw]*w[jw]; 

	/*              Multiplication par blocs 

              K_ij        = DP^T  DD_ij      DP 
         (Nnode x Nnode)        (3 x 3)  (3x Nnode)      */

	for (ji=0;ji<3;ji++) {
	  for (jj=0;jj<3;jj++) { 

          /* temp = dd_ij*Dp */
	    of_dgemm(&trans,&norm,&trois,&Nnode,&trois,&un,DD+9*(3*ji)+3*jj,&neuf,
		     ndnE+Nnode*4*jw+1,&quatre,&zero,TEMP,&trois); 

          /* K_ij = K_ij+jdet[jw]*w[jw] Dp^t Temp*/

	    of_dgemm(&trans,&norm,&Nnode,&Nnode,&trois,&alpha, ndnE+Nnode*4*jw+1,
		     &quatre,TEMP,&trois,&un, ke+Nnode*ji+Mk*Nnode*jj,&Mk);
	  }
      
      /* B = D*Dp (rhs)*/
      
	  of_dgemm(&trans,&norm,&six,&Nnode,&trois,&un,D+3*ji,&neuf,
		   ndnE+Nnode*4*jw+1,&quatre,&zero,B+6*ji*Nnode,&six); 
     
	}  

    /* Be = Sigma_i*B[i} (rhs)*/
	for (ji=0;ji<6;ji++) {
          beta =alpha*Sigma[ji];
	  of_daxpy(&Mk,&beta,B+ji,&six,Be,&unit);
	}


	free(B);

	
      }

      /*tabel(nodeE, constit, ke, Nnode, Ndof);*/
      #pragma omp critical
      {
	RHS=def+Mdef; 
	for (jk=0;jk<3;jk++) {
	  for (jj=0;jj<Nnode;jj++) {
	    /*#pragma omp atomic*/
	    RHS[CurDofPos[jk+3*jj]]+=Be[Nnode*jk+jj];
	  }
	}
      

      /*#pragma omp critical
      {*/
	AssembleSparse(o_ir,o_jc,o_pr,opt[0],opt[1],opt[2],
		       CurDofPos,elmap,NULL,NULL,ke);
      }

#ifdef TESTALLOCWORKS      
      free(ke);
      free(Be);
      free(nodeE);
      free(defelem);
      free(ndnE);
      free(jdetE);
#endif

    } /* end loop on elements */
#ifdef TESTALLOCFAIL      
      free(ke);
      free(Be);
      free(nodeE);
      free(defelem);
      free(ndnE);
      free(jdetE);
#endif


     
  }
    
  
  return;

}
