/* of_mk_subs.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "mex.h"
#include "../mex/of_mk_pre.h"
#if UseModulef>=1
 #include "f2c.h"
#endif

/* Table of constant values */
extern void compute_matrix(char*,int* iopt,double* constit, double* coor, 
        double* out, double* estate, double* defe, double* eltconst, double* infoatnode);
extern  void compute_b(char*, int* point,int* iopt,double* constit,
        double* coor, double* est,double* defe,double* eltconst,
        double* infoatnode, double* out);
extern void compute_stress(char*,int* iopt,double* constit,double* coor, 
        double* estate,double* defe,double* eltconst,double* infoatnode,
        double* out);

int c__1 = 1;
int c__4 = 4;
int c__0 = 0;
int c__2 = 2;

/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
char* pre_cvs () {
 char *st=(char*)"$Revision: 1.14 $  $Date: 2017/05/09 10:54:15 $";
 return st;
}

int io_switch(char* cam,char* ierr,int* pointers, int* integ, double* constit, 
  double* node,double*  estate,double* defe, double* eltconst,double*  out,double*  out1, double* out2,double*  infoatnode) {
    /* System generated locals */
  int *iopt;

    /* Local variables */
    int i1, typ;

    /* Parameter adjustments */
    --eltconst;
    --estate;

    /* Function Body */
    constit=constit+pointers[6];
    iopt=integ+pointers[5]+2;
    typ = iopt[2];
#if UseModulef>=1
    if (typ == 0) {
	iopt[2] = 1;
	compute_matrix(cam,iopt,constit, node, out, &estate[1],defe, &eltconst[1], infoatnode);
	iopt[2] = 2;
	compute_matrix(cam,iopt,constit, node, out1, &estate[1],defe, &eltconst[1], infoatnode);
	iopt[2] = 0;
    } else if (typ == 1 || typ == 2) {
	compute_matrix(cam,iopt,constit, node, out, &estate[1], defe, &eltconst[1], infoatnode);

/*     hysteretic damping matrix */
    } else if (typ == 4) {
	iopt[2] = 1;
	compute_matrix(cam,iopt,constit, node, out, &estate[1], defe, &eltconst[1], infoatnode);
	for (i1 = 0; i1 < pointers[1]; ++i1) { out[i1] *= constit[1];}
	iopt[2] = 4;
/*     volume forces */
    } else if (typ == 100) {
	compute_b(cam,pointers,iopt,constit, node, &estate[1], defe, &eltconst[1], infoatnode, out);
    } else if (typ == 200) {
	compute_stress(cam,iopt,constit, node, &estate[1], defe, &eltconst[1], infoatnode, out);
    } else {
      /* return a null matrix */
    }
    return 0;
#else
    return 1;
#endif
} /* io_switch */

/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* Subroutine */
void elt_condense(int i1, double* ke, double* kr, 
                          double* me, double* mr, int M) {

  int    j1, j2, j3; /* ind is the last index */
  double *vect, r1;
  
  for (j1=M-1;j1>=i1;j1--) { /* loop on columns to suppress */
    ke[j1+M*j1]=1./ke[j1+M*j1]; r1=ke[j1+M*j1];
    for (j2=0;j2<j1;j2++) { ke[j2+M*j1]*=-r1; }
    vect=ke+M*j1; /* ke(1:j1-1,j1) */
    for (j2=0;j2<j1;j2++) { /* change into BLAS call */
      for (j3=0;j3<j1;j3++) { 
        ke[j3+M*j2]-= 1./r1 * *(vect+j3) * *(vect+j2);
        if (me!=0) {
	  me[j3+M*j2]+=  - *(vect+j2) * me[j1+M*j3] - *(vect+j3) * me[M*j2+j1] 
                        + *(vect+j3) * me[j1+M*j1] * *(vect+j2);
	}
      } /* j3 */
    } /* j2 */
  } /* j1 */

  if (ke!=kr) {
    for (j1=0;j1<i1;j1++) for (j2=0;j2<i1;j2++) { kr[j1+i1*j2]=ke[j1+M*j2];  }
  } /* if */
  if ((me!=0) && (me!=mr)) {
    for (j1=0;j1<i1;j1++) for (j2=0;j2<i1;j2++) { mr[j1+i1*j2]=me[j1+M*j2];  }
  } /* if */


}
#if UseModulef==1
/* ----------------------------------------------------------------------- */
/* Subroutine */ 
void compute_matrix(char* cam,int* iopt, doublereal* constit, 
        doublereal* coor, doublereal* out, doublereal* estate, 
        doublereal* defe, doublereal* eltconst, doublereal* infoatnode) {

    /* Local variables */
    extern /* Subroutine */ int emaq2c_(), eraq2c_(), etm2p1d_(), etm3p1d_(), etm3q1d_(), etm3r1d_(), etm3p2c_(), etm3q2c_(), etr3p1d_(), etr3q1d_(), etr3r1d_(), etr3p2c_(), etr3q2c_(), etr3r2c_(), etm3r2c_(), etr2q1d_(), etm2q1d_(), etr2q2c_(), etm2q2c_(), etr2p1d_(), etr2p2c_(), etm2p2c_(), etmap1d_(), etmaq1d_(), etmaq2c_(), etmap2c_(), etrap1d_(), etraq1d_(), etraq2c_(), etrap2c_(), etm5noe_(), etr5noe_(), etrmit4_(), etmdktp_(), etrdktp_();

/* ----------------------------------------------------------------------- */
/* 3D ELEMENTS */
/* ----------------------------------------------------------------------- */

/*     CONSTIT  -> Car : caracteristic */
/*       StIFfness */
/*         IOPT=2  :   E, nu */
/*         IOPT=9  :   orthotropic material constants */
/*       Mass */
/*         IOPT=1  :   rho */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --infoatnode;
    --eltconst;
    --defe;
    --estate;
    --out;
    --coor;
    --constit;
    --iopt;

    /*  mexPrintf(">%i,%i,%i,%i,%i<",iopt[1],iopt[2],iopt[3],iopt[4],iopt[5]);*/

    /* Function Body */
    if (!strcmp(cam, "tetra4")) {
	if (iopt[3] == 1) {
	    etr3p1d_(&coor[1], &constit[3], &iopt[6], &out[1]);
	} else if (iopt[3] == 2) {
	    etm3p1d_(&coor[1], &constit[1], &iopt[6], &out[1]);
	}
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    } else if (!strcmp(cam, "hexa8")) {
	if (iopt[3] == 1) {
	    etr3q1d_(&coor[1], &constit[3], &iopt[6], &out[1]);
	} else if (iopt[3] == 2) {
	    etm3q1d_(&coor[1], &constit[1], &iopt[6], &out[1]);
	}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    } else if (!strcmp(cam, "penta6")) {
	if (iopt[3] == 1) {
	    etr3r1d_(&coor[1], &constit[3], &iopt[6], &out[1]);
	} else if (iopt[3] == 2) {
	    etm3r1d_(&coor[1], &constit[1], &iopt[6], &out[1]);
	}
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    } else if (!strcmp(cam, "tetra10")) {
	if (iopt[3] == 1) {
	    etr3p2c_(&coor[1], &constit[3], &iopt[6], &out[1]);
	} else if (iopt[3] == 2) {
	    etm3p2c_(&coor[1], &constit[1], &iopt[6], &out[1]);
	}
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    } else if (!strcmp(cam, "hexa20")) {
	if (iopt[3] == 1) {
	    etr3q2c_(&coor[1], &constit[3], &iopt[6], &out[1]);
	} else if (iopt[3] == 2) {
	    etm3q2c_(&coor[1], &constit[1], &iopt[6], &out[1]);
	}
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    } else if (!strcmp(cam, "penta15")) {
	if (iopt[3] == 1) {
	    etr3r2c_(&coor[1], &constit[3], &iopt[6], &out[1]);
	} else if (iopt[3] == 2) {
	    etm3r2c_(&coor[1], &constit[1], &iopt[6], &out[1]);
	}
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    } else if (!strcmp(cam, "dktp")) {
	if (iopt[3] == 1) {
	    etrdktp_(&coor[1], &constit[3], &iopt[6], &out[1]);
	} else if (iopt[3] == 2) {
	    etmdktp_(&coor[1], &constit[1], &iopt[6], &out[1]);
	}
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    } else if (!strcmp(cam, "mitc4")) {
	if (iopt[3] == 1) {
/*                    (COOR,CAR, IOPT,          POL,  XNORM,  ae) */
	    etrmit4_(&coor[1], &constit[3], &iopt[6], &eltconst[1], &infoatnode[1], &out[1]);
	} else if (iopt[3] == 2) {
	}
/* ----------------------------------------------------------------------- */
/* 2D PLANE ELEMENTS */
/*
   N=8;
   i1=find(triu(ones(N,N))); r1=zeros(N,N);r1(i1)=1:length(i1);
   out2=r1+tril(r1',-1); ind=diag(out2);
   fprintf('\n                ');for j1=length(ind):-1:2
    fprintf('out[%i]=out[%i];',ind(j1),j1);
    if ~any(ind==j1) fprintf('out[%i]=0;',j1);end
   end
   fprintf('\n');

 ----------------------------------------------------------------------- */
    } else if (!strcmp(cam, "q4p")) {
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - 2D STRESS/STRAIN */


	if (iopt[5] == 1 || iopt[5] == 2) {
	    if (iopt[3] == 1) {
		etr2q1d_(&coor[1], &constit[3], &iopt[6], &out[1]);
	    } else if (iopt[3] == 2) {
		etm2q1d_(&coor[1], &constit[1], &iopt[7], &out[1]);
         out[36]=out[8];out[8]=0;
         out[28]=out[7];out[7]=0;
         out[21]=out[6];
         out[15]=out[5];out[5]=0;
         out[10]=out[4];out[4]=0;
         out[6]=out[3];
         out[3]=out[2];out[2]=0;
	    }
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - 2D AXI ELEMENTS */
	} else if (iopt[5] == 3) {
	    if (iopt[3] == 1) {
		etraq1d_(&coor[1], &constit[3], &iopt[6], &out[1]);
	    } else if (iopt[3] == 2) {
		etmaq1d_(&coor[1], &constit[1], &iopt[7], &out[1]);
	    }
	}
    } else if (!strcmp(cam, "q8p")) {
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - 2D STRESS/STRAIN */
	if (iopt[5] == 1 || iopt[5] == 2) {
	    if (iopt[3] == 1) {
		etr2q2c_(&coor[1], &constit[3], &iopt[6], &out[1]);
	    } else if (iopt[3] == 2) {
		etm2q2c_(&coor[1], &constit[1], &iopt[7], &out[1]);
	    }
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - 2D AXI ELEMENTS */
	} else if (iopt[5] == 3) {
	    if (iopt[3] == 1) {
		etraq2c_(&coor[1], &constit[3], &iopt[6], &out[1]);
	    } else if (iopt[3] == 2) {
		etmaq2c_(&coor[1], &constit[1], &iopt[7], &out[1]);
	    }
	}
    } else if (!strcmp(cam, "q9a")) {
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - 2D AXI FOURIER */
	if (iopt[5] == 3) {
	    if (iopt[3] == 1) {
		eraq2c_(&iopt[6], &coor[1], &coor[10], &constit[3], &constit[4], &out[1]);
	    } else if (iopt[3] == 2) {
		emaq2c_(&iopt[6], &coor[1], &coor[10], &constit[1], &out[1]);
	    }
	}
    } else if (!strcmp(cam, "t3p")) {
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - 2D STRESS/STRAIN */

	if (iopt[5] == 1 || iopt[5] == 2) {
	    if (iopt[3] == 1) {
		etr2p1d_(&coor[1], &constit[3], &iopt[6], &out[1]);
	    } else if (iopt[3] == 2) {
		etm2p1d_(&coor[1], &constit[1], &iopt[7], &out[1]);
	    }
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - 2D AXI ELEMENTS */
	} else if (iopt[5] == 3) {
	    if (iopt[3] == 1) {
		etrap1d_(&coor[1], &constit[3], &iopt[6], &out[1]);
	    } else if (iopt[3] == 2) {
		etmap1d_(&coor[1], &constit[1], &iopt[7], &out[1]);
	    }
	}
    } else if (!strcmp(cam, "t6p")) {
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - 2D STRESS/STRAIN */
	if (iopt[5] == 1 || iopt[5] == 2) {
	    if (iopt[3] == 1) {
		etr2p2c_(&coor[1], &constit[3], &iopt[6], &out[1]);
	    } else if (iopt[3] == 2) {
		etm2p2c_(&coor[1], &constit[1], &iopt[7], &out[1]);
	    }
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - 2D AXI ELEMENTS */
	} else if (iopt[5] == 3) {
	    if (iopt[3] == 1) {
		etrap2c_(&coor[1], &constit[3], &iopt[6], &out[1]);
	    } else if (iopt[3] == 2) {
		etmap2c_(&coor[1], &constit[1], &iopt[7], &out[1]);
	    }
	}
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - 2D Q5P */
    } else if (!strcmp(cam, "q5p")) {
	if (iopt[3] == 1) {
	    etr5noe_(&coor[1], &constit[3], &iopt[6], &out[1]);
	} else if (iopt[3] == 2) {
	    etm5noe_(&coor[1], &constit[1], &iopt[7], &out[1]);
	}
/* ----------------------------------------------------------------------- */
/* OTHERWISE ERROR */
/* ----------------------------------------------------------------------- */
    } else {
      mexPrintf("'%s'\n",cam);
      mexErrMsgTxt("Error computing matrix : not implemented");
    }
} /* compute_matrix */

/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* Subroutine */ void compute_b(char* cam,int* point, int* iopt, 
    doublereal* constit, doublereal* coor, doublereal* est, 
   doublereal* defe, doublereal* eltconst,doublereal*  infoatnode, doublereal* out) {

    /* Local variables */
    extern /* Subroutine */ int ets2p1d_(), ets3p1d_(), ets3q1d_(), ets3r1d_(), ets3p2c_(), ets3q2c_(), ets3r2c_(), ets2q1d_(), ets2p2c_(), ets2q2c_(), etsap1d_(), etsaq1d_(), etsap2c_(), etsaq2c_(), ets5noe_(), etsmit4_(), etsdktp_();

/* ----------------------------------------------------------------------- */
/* 3D ELEMENTS */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --out;
    --infoatnode;
    --eltconst;
    --defe;
    --est;
    --coor;
    --constit;
    --iopt;
    --point;

    /* Function Body */
    if (!strcmp(cam, "tetra4")) {
/*                      [vol sur pre noref]  [3*3*4  3*4];cumsum(ans)+1 */
	ets3p1d_(&coor[1], &defe[1], &est[1], &est[37], &point[9], &out[1]);
    } else if (!strcmp(cam, "hexa8")) {
/*                 [vol sur pre noref] [3*4*6  3*6];cumsum(ans)+1 */
	ets3q1d_(&coor[1], &defe[1], &est[1], &est[73], &point[9], &out[1]);
    } else if (!strcmp(cam, "penta6")) {
/*                 [vol sur pre noref] [3*(3*2+4*3)  3*5];cumsum(ans)+1 */
	ets3r1d_(&coor[1], &defe[1], &est[1], &est[55], &point[9], &out[1]);
    } else if (!strcmp(cam, "tetra10")) {
/*                      [vol sur pre noref]  [3*6*4  3*4];cumsum(ans)+1 */
	ets3p2c_(&coor[1], &defe[1], &est[1], &est[73], &point[9], &out[1]);
    } else if (!strcmp(cam, "hexa20")) {
/*                 [vol sur pre noref] [3*8*6  3*6];cumsum(ans)+1 */
	ets3q2c_(&coor[1], &defe[1], &est[1], &est[145], &point[9], &out[1]);
    } else if (!strcmp(cam, "penta15")) {
/*                 [vol sur pre noref] [3*(6*2+8*3)  3*5];cumsum(ans)+1 */
	ets3r2c_(&coor[1], &defe[1], &est[1], &est[109], &point[9], &out[1]);
    } else if (!strcmp(cam, "dktp")) {
	etsdktp_(&coor[1], &defe[1], &out[1]);
    } else if (!strcmp(cam, "mitc4")) {
	etsmit4_(&coor[1], &defe[1], &iopt[6], &eltconst[1], &infoatnode[1], &out[1]);
/* ----------------------------------------------------------------------- */
/* 2D ELEMENTS */
/* ----------------------------------------------------------------------- */
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - 2D AXI ELEMENTS */
    } else if (!strcmp(cam, "q4p") && iopt[5] == 3) {
/*          [vol sur pre noref] [2*8];cumsum(ans)+1 */

	etsaq1d_(&coor[1], &defe[1], &est[1], &est[17], &point[9], &constit[9], &est[25], &constit[3], &iopt[6], &out[1]);
    } else if (!strcmp(cam, "t3p") && iopt[5] == 3) {
/*          [vol sur pre noref] [2*6];cumsum(ans)+1 */
	etsap1d_(&coor[1], &defe[1], &est[1], &est[13], &point[9], &constit[9], &est[19], &constit[3], &iopt[6], &out[1]);
    } else if (!strcmp(cam, "t6p") && iopt[5] == 3) {
/*          [vol sur pre noref] [2*9];cumsum(ans)+1 */
	etsap2c_(&coor[1], &defe[1], &est[1], &est[19], &point[9], &constit[9], &est[25], &constit[3], &iopt[6], &out[1]);
    } else if (!strcmp(cam, "q8p") && iopt[5] == 3) {
/*          [vol sur pre noref] [2*12];cumsum(ans)+1 */
	etsaq2c_(&coor[1], &defe[1], &est[1], &est[25], &point[9], &constit[9], &est[33], &constit[3], &iopt[6], &out[1]);
    } else if (!strcmp(cam, "q5p") && iopt[5] == 3) {
/*          CALL ets5noe( Coor,VOL,SUR,PRE,POINT(9),OUT) */
      mexErrMsgTxt("Error computing RHS : q5a not implemented");
      
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - 2D STRESS/STRAIN */
    } else if (!strcmp(cam, "q4p")) {
/*                      coor fv  fs  p    noref */
/*          [vol sur pre noref] [2*8];cumsum(ans)+1 */
	ets2q1d_(&coor[1], &defe[1], &est[1], &est[17], &point[9], &constit[9], &est[25], &constit[3], &iopt[6], &out[1]);
/*      OUT(1)=CONSTIT(3) */
    } else if (!strcmp(cam, "t3p")) {
	ets2p1d_(&coor[1], &defe[1], &est[1], &est[13], &point[9], &constit[9], &est[19], &constit[3], &iopt[6], &out[1]);
    } else if (!strcmp(cam, "t6p")) {
	ets2p2c_(&coor[1], &defe[1], &est[1], &est[19], &point[9], &constit[9], &est[25], &constit[3], &iopt[6], &out[1]);
    } else if (!strcmp(cam, "q8p")) {
	ets2q2c_(&coor[1], &defe[1], &est[1], &est[25], &point[9], &out[1]);
    } else if (!strcmp(cam, "q5p")) {
	ets5noe_(&coor[1], &defe[1], &est[1], &est[17], &point[9], &out[1]);
/* ----------------------------------------------------------------------- */
/* 2D AXI Fourier ELEMENTS */
/* ----------------------------------------------------------------------- */
    } else if (!strcmp(cam, "q9a")) {
      mexErrMsgTxt("Error computing RHS : q9a not implemented");

/* ----------------------------------------------------------------------- */
/* OTHERWISE ERROR */
/* ----------------------------------------------------------------------- */
    } else {
      mexErrMsgTxt("Error computing RHS : not implemented");
    }
} /* compute_b */

/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* Subroutine */ 
void compute_stress(char* cam,int* iopt,doublereal* constit,doublereal* coor, 
        doublereal* estate, doublereal* defe, doublereal* eltconst, 
        doublereal* infoatnode, doublereal* out) {
    /* System generated locals */
  /* real r__1, r__2; */

    /* Local variables */
    extern /* Subroutine */ int etc2p1d_(), etc3p1d_(), etc3q1d_(), etc3r1d_(), etc3p2c_(), etc3q2c_(), etc3r2c_(), etc2q1d_(), etc2p2c_(), etc2q2c_(), etcap1d_(), etcaq1d_(), etcap2c_(), etcaq2c_(), etc5noe_(), etcdktp_();

/*      real*8        IN(*), Coor(*), OUT(*), DEF(*), TEM(*) */
/* ----------------------------------------------------------------------- */
/* 3D ELEMENTS */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --out;
    --coor;
    --constit;
    --iopt;
    --estate; 

    /* Function Body */
    if (!strcmp(cam, "tetra4")) {
	etc3p1d_(&coor[1], &constit[3], &iopt[6], &out[1], defe);
    } else if (!strcmp(cam, "hexa8")) {
	etc3q1d_(&coor[1], &constit[3], &iopt[6], &out[1], defe);
    } else if (!strcmp(cam, "penta6")) {
	etc3r1d_(&coor[1], &constit[3], &iopt[6], &out[1], defe);
    } else if (!strcmp(cam, "tetra10")) {
	etc3p2c_(&coor[1], &constit[3], &iopt[6], &out[1], defe);
    } else if (!strcmp(cam, "hexa20")) {
	etc3q2c_(&coor[1], &constit[3], &iopt[6], &out[1], defe);
    } else if (!strcmp(cam, "penta15")) {
	etc3r2c_(&coor[1], &constit[3], &iopt[6], &out[1], defe);
    } else if (!strcmp(cam, "dktp")) {
	etcdktp_(&coor[1], &constit[3], &iopt[6], defe, &out[1]);
/* ----------------------------------------------------------------------- */
/* 2D AXI ELEMENTS */
/* ----------------------------------------------------------------------- */
    } else if (iopt[5] == 3) {
	if (!strcmp(cam, "t3p")) {
	    etcap1d_(&coor[1], &constit[3], &iopt[6], defe, &estate[1], &estate[4], &out[1]);
	} else if (!strcmp(cam, "q4p")) {
	    etcaq1d_(&coor[1], &constit[3], &iopt[6], defe, &estate[1], &estate[4], &out[1]);
	} else if (!strcmp(cam, "t6p")) {
	    etcap2c_(&coor[1], &constit[3], &iopt[6], defe, &estate[1], &estate[4], &out[1]);
	} else if (!strcmp(cam, "q8p")) {
	    etcaq2c_(&coor[1], &constit[3], &iopt[6], defe, &estate[1], &estate[4], &out[1]);
	}
/* ----------------------------------------------------------------------- */
/* 2D ELEMENTS */
/* ----------------------------------------------------------------------- */
    } else if (!strcmp(cam, "t3p")) {
/* selection of element type based on IOPT(5) (see matrix) */
	etc2p1d_(&coor[1], &constit[3], &iopt[6], defe, &estate[1], &estate[4], &out[1]);
    } else if (!strcmp(cam, "q4p")) {
	etc2q1d_(&coor[1], &constit[3], &iopt[6], defe, &estate[1], &estate[4], &out[1]);
    } else if (!strcmp(cam, "t6p")) {
	etc2p2c_(&coor[1], &constit[3], &iopt[6], defe, &estate[1], &estate[4], &out[1]);
    } else if (!strcmp(cam, "q8p")) {
	etc2q2c_(&coor[1], &constit[3], &iopt[6], defe, &out[1]);
    } else if (!strcmp(cam, "q5p")) {
	etc5noe_(&coor[1], &constit[3], &iopt[6], defe, &out[1]);
    } else {
      mexErrMsgTxt("Error computing stresses : not implemented");
    }
} /* compute_stress */

/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* Subroutine */ int idof_(char* cam, int* iopt, doublereal* coor, doublereal* out) {
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    int i1, j1, j2, j3, i2;

/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --out;
    --coor;
    --iopt;

    /* Function Body */
    if (!strcmp(cam, "tetra4") || !strcmp(cam, "hexa8") || !strcmp(cam, "penta6") || !strcmp(cam, "tetra10") || !strcmp(cam, "hexa20") || !strcmp(cam, "dktp") || !strcmp(cam, "penta15") || !strcmp(cam, "q9a")) {
	i2 = 3;
    } else if (!strcmp(cam, "q4p") || !strcmp(cam, "q8p") || !strcmp(cam, "t3p") || !strcmp(cam, "t6p") || !strcmp(cam, "q5p")) {
	i2 = 2;
    } else if (!strcmp(cam, "mitc4")) {
	i2 = 5;
    } else {
      mexErrMsgTxt("Error computing DOFs : not implemented");
    }
    i1 = iopt[2];
    if (!strcmp(cam, "dktp")) {
	j3 = 1;
	i__1 = i1;
	for (j1 = 1; j1 <= i__1; ++j1) {
	    i__2 = i2;
	    for (j2 = 1; j2 <= i__2; ++j2) {
		out[j3] = coor[i1 * 3 + j1] + (float).02 + j2 * (float).01;
		++j3;
/* L12: */
	    }
/* L11: */
	}
    } else if (!strcmp(cam, "mitc4")) {
	j3 = 1;
	i__1 = i1;
	for (j1 = 1; j1 <= i__1; ++j1) {
	    i__2 = i2;
	    for (j2 = 1; j2 <= i__2; ++j2) {
		out[j3] = coor[i1 * 3 + j1] + j2 * (float).01;
		++j3;
/* L14: */
	    }
/* L13: */
	}
    } else {
	j3 = 1;
	i__1 = i1;
	for (j1 = 1; j1 <= i__1; ++j1) {
	    i__2 = i2;
	    for (j2 = 1; j2 <= i__2; ++j2) {
		out[j3] = coor[i1 * 3 + j1] + j2 * (float).01;
		++j3;
/* L10: */
	    }
/* L9: */
	}
    }
    return 0;
} /* idof_ */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* Subroutine */ int edemit4_(doublereal* pol) {
    /* Initialized data */

    static doublereal yty[4] = { -1.,1.,0.,0. };
    static doublereal xint[4] = { -.5773502691896258,.5773502691896258,.5773502691896258,-.5773502691896258 };
    static doublereal yint[4] = { -.5773502691896258,-.5773502691896258,.5773502691896258,.5773502691896258 };
    static doublereal xty[4] = { 0.,0.,-1.,1. };

    int l, n;
    extern doublereal pcq12d_(int*,int*,doublereal*,doublereal*);
    int lt;

/* -----*--*---------*---------*---------*---------*---------*---------*-- */

/* BUT : CALCUL DES POLYNOMES DE BASE ET DERIVEES POUR MIT4 */
/* --- */
/*       Le tableau pol de 96 variables doit etre declare ds matlab et rempli */
/*       si l'elem mitc4 est utilise */
/* -----*--*---------*---------*---------*---------*---------*---------*-- */
/*     pol contient en sequence vp2(nno,npi+npt)   cad    vp2(4,4+4) */
/*                              vdpq2(2,nno,npi)   cad    vdpq2(2,4,4) */
/*                              vdtq2(2,nno,npt)   cad    vdtq2(2,4,4) */



/*     tying points */

    /* Parameter adjustments */
    --pol;

    /* Function Body */

/*     integration points en r,s */


/*        remplir les tableaux vdpq2(2,nno,npi),vdtq2(2,nno,npt), */
/*                             vp2(nno,npi+npt) */

    for (l = 1; l <= 4; ++l) {
	for (n = 1; n <= 4; ++n) {
/*           -- > vp2(n,l)     = pcq12d(0,n,xint(l),yint(l)) */
	    pol[((l - 1) << 2) + n] = pcq12d_(&c__0, &n, &xint[l - 1], &yint[l - 1]);
/*           -- > vdpq2(1,n,l) = pcq12d(1,n,xint(l),yint(l)) */
	    pol[((l - 1) << 3) + 33 + ((n - 1) << 1)] = pcq12d_(&c__1, &n, &xint[l - 1], &yint[l - 1]);
/*           -- > vdpq2(2,n,l) = pcq12d(2,n,xint(l),yint(l)) */
	    pol[((l - 1) << 3) + 34 + ((n - 1) << 1)] = pcq12d_(&c__2, &n, &xint[l - 1], &yint[l - 1]);
/* L1: */
	}
    }
    for (lt = 1; lt <= 4; ++lt) {
	for (n = 1; n <= 4; ++n) {
/*           -- >  vp2(n,npi+lt) = pcq12d(0,n,xty(lt),yty(lt)) */
	    pol[((lt + 3) << 2) + n] = pcq12d_(&c__0, &n, &xty[lt - 1], &yty[lt - 1]);
/*           -- >  vdtq2(1,n,lt) = pcq12d(1,n,xty(lt),yty(lt)) */
	    pol[((lt - 1) << 3) + 65 + ((n - 1) << 1)] = pcq12d_(&c__1, &n, &xty[lt - 1], &yty[lt - 1]);
/*           -- >  vdtq2(2,n,lt) = pcq12d(2,n,xty(lt),yty(lt)) */
	    pol[((lt - 1) << 3) + 66 + ((n - 1) << 1)] = pcq12d_(&c__2, &n, &xty[lt - 1], &yty[lt - 1]);
/* L2: */
	}
    }
} /* edemit4_ */

#endif

