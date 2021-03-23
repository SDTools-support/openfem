/* chan57mod.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int chan57mod_(nno, ae5, lae5, ae6, lae6, t, e, coor)
int32 *nno;
doublereal *ae5;
int32 *lae5;
doublereal *ae6;
int32 *lae6;
doublereal *t, *e, *coor;
{
    /* System generated locals */
    int32 coor_dim1, coor_offset, i__1, i__2;
    real r__1;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static int32 ncar, nloc;
    static doublereal mesk;
    static int32 i__, j, n;
    static real t1;
    static doublereal t2, t3, t4;
    static int32 k51, k61, k62, k52;
    static doublereal coef11x, coef12x, coef22x, coef14x, coef23x, coef33x, coef34x, coef12y, coef14y, coef22y, coef23y, coef33y, coef34y;
    static real coef11y, coef44x, coef44y;

/*     ............................................................ */

/*     transformer la matrice de rigidite de 5 a 6 dl/noeud */
/*     pour cas Coventor plaque mince */

/*     ............................................................. */
/*     Programmer: Matthieu Sors Coventor/INRIA      2004 */
/*     see Mathematica file for additional informations */
/*     ............................................................. */


/*     ............................................................. */
/*     Variables */

/*     ae5: Stiffness matrix without drilling */
/*     ae6: Stiffness matrix with drilling */
/*     Mesk: Area of the element */
/*     penal: coefficient of penalization */
/*     COOR: Coordinates of the nodes */
/*     ............................................................. */

/*    --------------- boucle sur les noeuds   ----------------- */
/*        ----- lignes */

/*     .................................................................. */
/*     Initialize ae6 */
/*     .................................................................. */

    /* Parameter adjustments */
    coor_dim1 = *nno;
    coor_offset = coor_dim1 + 1;
    coor -= coor_offset;
    --t;
    --ae5;
    --ae6;

    /* Function Body */
    for (n = 1; n <= 300; ++n) {
	ae6[n] = (float)0.;
/* L99: */
    }

/*     .................................................................. */
/*     Calculate area of the element */
/*     .................................................................. */

    mesk = (coor[coor_dim1 + 1] * coor[(coor_dim1 << 1) + 2] - coor[coor_dim1 + 2] * coor[(coor_dim1 << 1) + 1] + coor[coor_dim1 + 2] * coor[(coor_dim1 << 1) + 3] - coor[(coor_dim1 << 1) + 2] * coor[coor_dim1 + 3] + coor[coor_dim1 + 3] * coor[(coor_dim1 << 1) + 4] - coor[(coor_dim1 << 1) + 3] * coor[coor_dim1 + 4] + coor[coor_dim1 + 4] * coor[(coor_dim1 << 1) + 1] - coor[(coor_dim1 << 1) + 4] * coor[coor_dim1 + 1]) / 2;
/*     ................................................................... */
/*     Evaluate Coefficients of the Energy of penalization */
/*     ................................................................... */

/* Computing 2nd power */
    d__1 = -coor[coor_dim1 + 1] + coor[coor_dim1 + 2];
/* Computing 2nd power */
    d__2 = -coor[(coor_dim1 << 1) + 1] + coor[(coor_dim1 << 1) + 2];
    coef12x = (coor[(coor_dim1 << 1) + 1] - coor[(coor_dim1 << 1) + 2]) / ((d__1 * d__1 + d__2 * d__2) * (float)2.);
/* Computing 2nd power */
    d__1 = -coor[coor_dim1 + 1] + coor[coor_dim1 + 2];
/* Computing 2nd power */
    d__2 = -coor[(coor_dim1 << 1) + 1] + coor[(coor_dim1 << 1) + 2];
    coef12y = (-coor[coor_dim1 + 1] + coor[coor_dim1 + 2]) / ((d__1 * d__1 + d__2 * d__2) * (float)2.);
/* Computing 2nd power */
    d__1 = coor[coor_dim1 + 1] - coor[coor_dim1 + 4];
/* Computing 2nd power */
    d__2 = coor[(coor_dim1 << 1) + 1] - coor[(coor_dim1 << 1) + 4];
    coef14x = -(-coor[(coor_dim1 << 1) + 1] + coor[(coor_dim1 << 1) + 4]) / ((d__1 * d__1 + d__2 * d__2) * (float)2.);
/* Computing 2nd power */
    d__1 = coor[coor_dim1 + 1] - coor[coor_dim1 + 4];
/* Computing 2nd power */
    d__2 = coor[(coor_dim1 << 1) + 1] - coor[(coor_dim1 << 1) + 4];
    coef14y = (-coor[coor_dim1 + 1] + coor[coor_dim1 + 4]) / ((d__1 * d__1 + d__2 * d__2) * (float)2.);
/* Computing 2nd power */
    d__1 = -coor[coor_dim1 + 2] + coor[coor_dim1 + 3];
/* Computing 2nd power */
    d__2 = -coor[(coor_dim1 << 1) + 2] + coor[(coor_dim1 << 1) + 3];
    coef23x = (coor[(coor_dim1 << 1) + 2] - coor[(coor_dim1 << 1) + 3]) / ((d__1 * d__1 + d__2 * d__2) * (float)2.);
/* Computing 2nd power */
    d__1 = -coor[coor_dim1 + 2] + coor[coor_dim1 + 3];
/* Computing 2nd power */
    d__2 = -coor[(coor_dim1 << 1) + 2] + coor[(coor_dim1 << 1) + 3];
    coef23y = (-coor[coor_dim1 + 2] + coor[coor_dim1 + 3]) / ((d__1 * d__1 + d__2 * d__2) * (float)2.);
/* Computing 2nd power */
    d__1 = -coor[coor_dim1 + 3] + coor[coor_dim1 + 4];
/* Computing 2nd power */
    d__2 = -coor[(coor_dim1 << 1) + 3] + coor[(coor_dim1 << 1) + 4];
    coef34x = (coor[(coor_dim1 << 1) + 3] - coor[(coor_dim1 << 1) + 4]) / ((d__1 * d__1 + d__2 * d__2) * (float)2.);
/* Computing 2nd power */
    d__1 = -coor[coor_dim1 + 3] + coor[coor_dim1 + 4];
/* Computing 2nd power */
    d__2 = -coor[(coor_dim1 << 1) + 3] + coor[(coor_dim1 << 1) + 4];
    coef34y = (-coor[coor_dim1 + 3] + coor[coor_dim1 + 4]) / ((d__1 * d__1 + d__2 * d__2) * (float)2.);

    coef11x = -coef12x - coef14x;
    coef11y = -coef12y - coef14y;
    coef22x = coef12x - coef23x;
    coef22y = coef12y - coef23y;
    coef33x = coef23x - coef34x;
    coef33y = coef23y - coef34y;
    coef44x = coef14x + coef34x;
    coef44y = coef14y + coef34y;


    t1 = t[1];
    t2 = t[2];
    t3 = t[3];
    t4 = t[4];

/*      .......................................................... */
/*      Test */
/*      .......................................................... */

/*      if(COOR(3,2) .eq. 1) then */
/*          if(COOR(3,1) .eq. 2) then */
/*      temp=MesK */
/*      mp = matOpen('inputm.mat', 'w') */

/*      arN = mxCreateDoubleMatrix(1,1,0) */
/*      call mxCopyReal8ToPtr(MesK,mxGetPr(arN),1) */
/*      call mxSetName(arN,'MesK') */
/*      status = matPutMatrix(mp,arN) */
/*      status = matClose(mp) */
/*          endif */
/*      endif */

/*      ................................................................. */
/*      Copy ae5 in ae6 */
/*      ................................................................. */

    i__1 = *nno;
    for (n = 1; n <= i__1; ++n) {
	if (n != 1) {
	    nloc = n - 1;
	    i__2 = nloc;
	    for (ncar = 1; ncar <= i__2; ++ncar) {
		for (i__ = 1; i__ <= 5; ++i__) {
		    k61 = nloc * 6 + i__;
		    k51 = nloc * 5 + i__;
		    for (j = 1; j <= 5; ++j) {
			k62 = (ncar - 1) * 6 + j;
			k52 = (ncar - 1) * 5 + j;

			ae6[k61 * (k61 - 1) / 2 + k62] = ae5[k51 * (k51 - 1) / 2 + k52];
/* L4: */
		    }
/* L3: */
		}
/* L2: */
	    }
	}


	for (i__ = 1; i__ <= 5; ++i__) {
	    k61 = (n - 1) * 6 + i__;
	    k51 = (n - 1) * 5 + i__;
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		k62 = (n - 1) * 6 + j;
		k52 = (n - 1) * 5 + j;

		ae6[k61 * (k61 - 1) / 2 + k62] = ae5[k51 * (k51 - 1) / 2 + k52];
/* L31: */
	    }
/* L21: */
	}
/* L1: */
    }

/*     ................................................................... */
/*     Adding Energy of compensation for drilling degrees of freedom */
/*     ................................................................... */

/* Computing 2nd power */
    d__1 = coef11x;
/* Computing 2nd power */
    d__2 = coef12x;
/* Computing 2nd power */
    d__3 = coef14x;
    ae6[1] += *e * mesk * 100. * (d__1 * d__1 * t1 + d__2 * d__2 * t2 + d__3 * d__3 * t4);
    ae6[2] += *e * mesk * 100. * (coef11x * coef11y * t1 + coef12x * coef12y * t2 + coef14x * coef14y * t4);
/* Computing 2nd power */
    r__1 = coef11y;
/* Computing 2nd power */
    d__1 = coef12y;
/* Computing 2nd power */
    d__2 = coef14y;
    ae6[3] += *e * mesk * 100. * (r__1 * r__1 * t1 + d__1 * d__1 * t2 + d__2 * d__2 * t4);
    ae6[16] -= coef11x * *e * mesk * 100. * t1;
    ae6[17] -= coef11y * *e * mesk * 100. * t1;
    ae6[21] += *e * mesk * 100. * t1;
    ae6[22] += *e * mesk * 100. * (coef11x * coef12x * t1 - coef12x * coef22x * t2);
    ae6[23] += *e * mesk * 100. * (coef11y * coef12x * t1 - coef12y * coef22x * t2);
    ae6[27] -= coef12x * *e * mesk * 100. * t1;
/* Computing 2nd power */
    d__1 = coef12x;
/* Computing 2nd power */
    d__2 = coef22x;
/* Computing 2nd power */
    d__3 = coef23x;
    ae6[28] += *e * mesk * 100. * (d__1 * d__1 * t1 + d__2 * d__2 * t2 + d__3 * d__3 * t3);
    ae6[29] += *e * mesk * 100. * (coef11x * coef12y * t1 - coef12x * coef22y * t2);
    ae6[30] += *e * mesk * 100. * (coef11y * coef12y * t1 - coef12y * coef22y * t2);
    ae6[34] -= coef12y * *e * mesk * 100. * t1;
    ae6[35] += *e * mesk * 100. * (coef12x * coef12y * t1 + coef22x * coef22y * t2 + coef23x * coef23y * t3);
/* Computing 2nd power */
    d__1 = coef12y;
/* Computing 2nd power */
    d__2 = coef22y;
/* Computing 2nd power */
    d__3 = coef23y;
    ae6[36] += *e * mesk * 100. * (d__1 * d__1 * t1 + d__2 * d__2 * t2 + d__3 * d__3 * t3);
    ae6[67] += coef12x * *e * mesk * 100. * t2;
    ae6[68] += coef12y * *e * mesk * 100. * t2;
    ae6[73] -= coef22x * *e * mesk * 100. * t2;
    ae6[74] -= coef22y * *e * mesk * 100. * t2;
    ae6[78] += *e * mesk * 100. * t2;
    ae6[79] += *e * mesk * 100. * (-(coef12x * coef23x * t2) + coef14x * coef34x * t4);
    ae6[80] += *e * mesk * 100. * (-(coef12y * coef23x * t2) + coef14y * coef34x * t4);
    ae6[85] += *e * mesk * 100. * (coef22x * coef23x * t2 - coef23x * coef33x * t3);
    ae6[86] += *e * mesk * 100. * (coef22y * coef23x * t2 - coef23y * coef33x * t3);
    ae6[90] -= coef23x * *e * mesk * 100. * t2;
/* Computing 2nd power */
    d__1 = coef23x;
/* Computing 2nd power */
    d__2 = coef33x;
/* Computing 2nd power */
    d__3 = coef34x;
    ae6[91] += *e * mesk * 100. * (d__1 * d__1 * t2 + d__2 * d__2 * t3 + d__3 * d__3 * t4);
    ae6[92] += *e * mesk * 100. * (-(coef12x * coef23y * t2) + coef14x * coef34y * t4);
    ae6[93] += *e * mesk * 100. * (-(coef12y * coef23y * t2) + coef14y * coef34y * t4);
    ae6[98] += *e * mesk * 100. * (coef22x * coef23y * t2 - coef23x * coef33y * t3);
    ae6[99] += *e * mesk * 100. * (coef22y * coef23y * t2 - coef23y * coef33y * t3);
    ae6[103] -= coef23y * *e * mesk * 100. * t2;
    ae6[104] += *e * mesk * 100. * (coef23x * coef23y * t2 + coef33x * coef33y * t3 + coef34x * coef34y * t4);
/* Computing 2nd power */
    d__1 = coef23y;
/* Computing 2nd power */
    d__2 = coef33y;
/* Computing 2nd power */
    d__3 = coef34y;
    ae6[105] += *e * mesk * 100. * (d__1 * d__1 * t2 + d__2 * d__2 * t3 + d__3 * d__3 * t4);
    ae6[160] += coef23x * *e * mesk * 100. * t3;
    ae6[161] += coef23y * *e * mesk * 100. * t3;
    ae6[166] -= coef33x * *e * mesk * 100. * t3;
    ae6[167] -= coef33y * *e * mesk * 100. * t3;
    ae6[171] += *e * mesk * 100. * t3;
    ae6[172] += *e * mesk * 100. * (coef11x * coef14x * t1 - coef14x * coef44x * t4);
    ae6[173] += *e * mesk * 100. * (coef11y * coef14x * t1 - coef14y * coef44x * t4);
    ae6[177] -= coef14x * *e * mesk * 100. * t1;
    ae6[178] += *e * mesk * 100. * (coef12x * coef14x * t1 - coef23x * coef34x * t3);
    ae6[179] += *e * mesk * 100. * (coef12y * coef14x * t1 - coef23y * coef34x * t3);
    ae6[184] += *e * mesk * 100. * (coef33x * coef34x * t3 - coef34x * coef44x * t4);
    ae6[185] += *e * mesk * 100. * (coef33y * coef34x * t3 - coef34y * coef44x * t4);
    ae6[189] -= coef34x * *e * mesk * 100. * t3;
/* Computing 2nd power */
    d__1 = coef14x;
/* Computing 2nd power */
    d__2 = coef34x;
/* Computing 2nd power */
    r__1 = coef44x;
    ae6[190] += *e * mesk * 100. * (d__1 * d__1 * t1 + d__2 * d__2 * t3 + r__1 * r__1 * t4);
    ae6[191] += *e * mesk * 100. * (coef11x * coef14y * t1 - coef14x * coef44y * t4);
    ae6[192] += *e * mesk * 100. * (coef11y * coef14y * t1 - coef14y * coef44y * t4);
    ae6[196] -= coef14y * *e * mesk * 100. * t1;
    ae6[197] += *e * mesk * 100. * (coef12x * coef14y * t1 - coef23x * coef34y * t3);
    ae6[198] += *e * mesk * 100. * (coef12y * coef14y * t1 - coef23y * coef34y * t3);
    ae6[203] += *e * mesk * 100. * (coef33x * coef34y * t3 - coef34x * coef44y * t4);
    ae6[204] += *e * mesk * 100. * (coef33y * coef34y * t3 - coef34y * coef44y * t4);
    ae6[208] -= coef34y * *e * mesk * 100. * t3;
    ae6[209] += *e * mesk * 100. * (coef14x * coef14y * t1 + coef34x * coef34y * t3 + coef44x * coef44y * t4);
/* Computing 2nd power */
    d__1 = coef14y;
/* Computing 2nd power */
    d__2 = coef34y;
/* Computing 2nd power */
    r__1 = coef44y;
    ae6[210] += *e * mesk * 100. * (d__1 * d__1 * t1 + d__2 * d__2 * t3 + r__1 * r__1 * t4);
    ae6[277] += coef14x * *e * mesk * 100. * t4;
    ae6[278] += coef14y * *e * mesk * 100. * t4;
    ae6[289] += coef34x * *e * mesk * 100. * t4;
    ae6[290] += coef34y * *e * mesk * 100. * t4;
    ae6[295] -= coef44x * *e * mesk * 100. * t4;
    ae6[296] -= coef44y * *e * mesk * 100. * t4;
    ae6[300] += *e * mesk * 100. * t4;

/*      .......................................................... */
/*      Test */
/*      .......................................................... */

/*        if(COOR(3,2) .eq. 1) then */
/*          if(COOR(3,1) .eq. 1) then */
/*        mp = matOpen('inputm2.mat', 'w') */

/*        are6 = mxCreateDoubleMatrix(300,1,0) */
/*        call mxCopyReal8ToPtr(ae6,mxGetPr(are6),300) */
/*        call mxSetName(are6,'ae6ap') */
/*        status = matPutMatrix(mp,are6) */
/*        status = matClose(mp) */
/*          endif */
/*        endif */

} /* chan57mod_ */

