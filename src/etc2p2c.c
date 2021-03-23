/* etc2p2c.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int etc2p2c_(coor, car, iopt, u, alpha, theta, sigma)
doublereal *coor, *car;
int32 *iopt;
doublereal *u, *alpha, *theta, *sigma;
{
    /* Initialized data */

    static int32 ijt[12] = { 1,3,5,7,9,11,2,4,6,8,10,12 };
    static doublereal dp[48]	/* was [2][6][4] */ = { -1.,-1.,1.,0.,0.,1.,0.,-4.,4.,4.,-4.,0.,-3.,-3.,-1.,0.,0.,-1.,4.,0.,0.,0.,0.,4.,1.,1.,3.,0.,0.,-1.,-4.,-4.,0.,4.,0.,0.,1.,1.,-1.,0.,0.,3.,0.,0.,4.,0.,-4.,-4. };

    static doublereal edfp[18]	/* was [3][6] */, c__, e[6];
    static int32 i__, j, l;
    static doublereal s, dfidp[12]	/* was [2][6] */;
    static int32 ibloc;
    static doublereal delta, dfinv[16]	/* was [2][2][4] */, const__[3];
    static int32 i1, i2, i3, i4, i5, i6, j1;
    static doublereal young, unmnu, ed[6];
    static int32 kk;
    static doublereal x21, y21, x31, y31, x32, y32, x41, y41, dsigma[36], x42, y42, x54, y54, x61, y61, x63, y63, x65, y65, ed1[3], eap[18]	/* was [3][6] */, poisson;

/* *************************************************************** */
/* BUT: CALCUL DES CONTRAINTES DE L ELEMENT TRIA 2p2c */
/* --- */
/* in : coor(noe,ndim) : coor. 3 sommets + 3 milieux aretes */
/*      car, iopt      : caracteristiques des materiaux */
/*      U(ndim,noe): deplacements U_x et U_y aux 6 noeuds */
/*      alpha(3)   : tenseur de dilatation thermique */
/*      theta(6)   : temperature aux 6 noeuds. */
/* out: SIGMA(3)   :  S_xx, S_yy, S_xy elastiques */

/* programmeur : modulef */
/* ............................................................... */
    /* Parameter adjustments */
    --sigma;
    --theta;
    --alpha;
    u -= 3;
    --car;
    coor -= 7;

    /* Function Body */

    x21 = coor[8] - coor[7];
    y21 = coor[14] - coor[13];
    x31 = coor[9] - coor[7];
    y31 = coor[15] - coor[13];
    x32 = coor[9] - coor[8];
    y32 = coor[15] - coor[14];
    x41 = coor[10] - coor[7];
    y41 = coor[16] - coor[13];
    x42 = coor[10] - coor[8];
    y42 = coor[16] - coor[14];
    x54 = coor[11] - coor[10];
    y54 = coor[17] - coor[16];
    x61 = coor[12] - coor[7];
    y61 = coor[18] - coor[13];
    x63 = coor[12] - coor[9];
    y63 = coor[18] - coor[15];
    x65 = coor[12] - coor[11];
    y65 = coor[18] - coor[17];

    if (*iopt == 1) {
/*  --    CONTRAINTES PLANES     (ISOTROPE)     ----- */
	young = car[1];
	poisson = car[2];
	c__ = young / (1. - poisson * poisson);
	e[0] = c__;
	e[1] = c__ * poisson;
	e[2] = c__;
	e[3] = 0.;
	e[4] = 0.;
	e[5] = c__ * (1. - poisson) / 2.;
    } else if (*iopt == 2) {
/*  --    DEFORMATIONS PLANES (ISOTROPE)     ----- */
	young = car[1];
	poisson = car[2];
	unmnu = 1. - poisson;
	c__ = young * unmnu;
	c__ /= (poisson + 1.) * (1. - poisson * 2.);
	e[0] = c__;
	e[1] = poisson * c__ / unmnu;
	e[2] = c__;
	e[3] = 0.;
	e[4] = 0.;
	e[5] = c__ * (1. - poisson * 2.) / (unmnu * 2.);
    } else {
/*  --    CAS ANISOTROPE     ----- */
	for (i__ = 1; i__ <= 6; ++i__) {
	    e[i__ - 1] = car[i__];
/* L1: */
	}
    }

/* -- CALCUL DES CONTRAINTES ELEMENTAIRES AU BARYCENTRE */

/*     --  CALCUL DE DFINV  L INVERSE DE DF */

    dfinv[0] = y31 + y54 * 4.;
    dfinv[1] = -x31 - x54 * 4.;
    dfinv[2] = -y21 + y65 * 4.;
    dfinv[3] = x21 - x65 * 4.;
    delta = dfinv[0] * dfinv[3] - dfinv[2] * dfinv[1];
    for (j = 1; j <= 6; ++j) {
	e[j - 1] /= delta;
/* L2: */
    }

/*     -- CALCUL DE DFINV*DP */

/*                | dp1dx dp2dx dp3dx dp3dx dp4dx dp5dx dp6dx | */
/*   DFIDP(2,6) = |                                           | */
/*                | dp1dy dp2dy dp3dy dp3dy dp4dy dp5dy dp6dy | */
    for (i__ = 1; i__ <= 2; ++i__) {
	for (j = 1; j <= 6; ++j) {
	    s = 0.;
	    for (l = 1; l <= 2; ++l) {
		s += dfinv[i__ + (l + 2 << 1) - 7] * dp[l + (j + 6 << 1) - 15];
/* L3: */
	    }
	    dfidp[i__ + (j << 1) - 3] = s;
/* L4: */
	}
    }

/*     --CALCUL DE E*D */
    for (ibloc = 1; ibloc <= 2; ++ibloc) {
	if (ibloc == 1) {
	    i1 = 1;
	    i2 = 2;
	    i3 = 4;
	    i4 = 4;
	    i5 = 5;
	    i6 = 6;
	} else if (ibloc == 2) {
	    i1 = 4;
	    i2 = 5;
	    i3 = 6;
	    i4 = 2;
	    i5 = 3;
	    i6 = 5;
	}
	ed[0] = e[i1 - 1];
	ed[1] = e[i2 - 1];
	ed[2] = e[i3 - 1];
	ed[3] = e[i4 - 1];
	ed[4] = e[i5 - 1];
	ed[5] = e[i6 - 1];

/*       -- CALCUL DE ED * DFIDP PAR BLOC */

	for (j = 1; j <= 6; ++j) {
	    for (i__ = 1; i__ <= 3; ++i__) {
		s = 0.;
		for (l = 1; l <= 2; ++l) {
		    kk = (l - 1) * 3 + i__;
		    s += ed[kk - 1] * dfidp[l + (j << 1) - 3];
/* L5: */
		}
		edfp[i__ + j * 3 - 4] = s;
/* L6: */
	    }
	}

/*       -- CALCUL DE DSIGMA(36) */

	for (j = 1; j <= 6; ++j) {
	    for (i__ = 1; i__ <= 3; ++i__) {
		j1 = j;
		if (ibloc == 2) {
		    j1 = j + 6;
		}
		kk = (ijt[j1 - 1] - 1) * 3 + i__;
		dsigma[kk - 1] = edfp[i__ + j * 3 - 4];
/* L7: */
	    }
	}
/* L8: */
    }

/*  [      11]   [1 4 7 10 13 16 19 22 25 28 31 34]   [          ] */
/*  [SIGMA 22] = [2 5 8 11 14 17 20 23 26 29 32 35] * [u_solution] */
/*  [      12]   [3 6 9 12 15 18 21 24 27 30 33 36]   [          ] */
/*          3*1                                 3*12           12*1 */
    sigma[1] = dsigma[0] * u[3] + dsigma[3] * u[4] + dsigma[6] * u[5] + dsigma[9] * u[6] + dsigma[12] * u[7] + dsigma[15] * u[8] + dsigma[18] * u[9] + dsigma[21] * u[10] + dsigma[24] * u[11] + dsigma[27] * u[12] + dsigma[30] * u[13] + dsigma[33] * u[14];
    sigma[2] = dsigma[1] * u[3] + dsigma[4] * u[4] + dsigma[7] * u[5] + dsigma[10] * u[6] + dsigma[13] * u[7] + dsigma[16] * u[8] + dsigma[19] * u[9] + dsigma[22] * u[10] + dsigma[25] * u[11] + dsigma[28] * u[12] + dsigma[31] * u[13] + dsigma[34] * u[14];
    sigma[3] = dsigma[2] * u[3] + dsigma[5] * u[4] + dsigma[8] * u[5] + dsigma[11] * u[6] + dsigma[14] * u[7] + dsigma[17] * u[8] + dsigma[20] * u[9] + dsigma[23] * u[10] + dsigma[26] * u[11] + dsigma[29] * u[12] + dsigma[32] * u[13] + dsigma[35] * u[14];
/*      print *,'---------------- Verif avec impressions Modulef' */
/*      print *, (dsigma(i), i= 1, 6) */
/*      print *, (dsigma(i), i= 7, 12) */
/*      print *, (dsigma(i), i= 13, 18) */
/*      print *, (dsigma(i), i= 19, 24) */
/*      print *, (dsigma(i), i= 24, 30) */
/*      print *, (dsigma(i), i= 31, 36) */

/*     LES CONTRAINTES THERMIQUES */
/*     -------------------------- */

/*     SIGMA(TETA) = - (E) * (ALPHA) * (P) (X,Y) */
/*     ED1 = (E) * (ALPHA) */
/*            3*3       3*1 */
    ed1[0] = -(e[0] * alpha[1] + e[1] * alpha[3] + e[3] * alpha[2]);
    ed1[1] = -(e[1] * alpha[1] + e[2] * alpha[3] + e[4] * alpha[2]);
    ed1[2] = -(e[3] * alpha[1] + e[4] * alpha[3] + e[5] * alpha[2]);

/*     ED1  * P   => EAP */
/*      3*1  1*6      3*6 */

/*     CONTRAINTES AU BARYCENTRE */
    for (i__ = 1; i__ <= 3; ++i__) {
	const__[i__ - 1] = ed1[i__ - 1] / (float)9.;
/* L9: */
    }
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    eap[i__ + j * 3 - 4] = -const__[i__ - 1];
	    eap[i__ + (j + 3) * 3 - 4] = const__[i__ - 1] * (float)4.;
/* L10: */
	}
    }

    for (j = 1; j <= 6; ++j) {
	sigma[1] += eap[j * 3 - 3] * theta[j];
	sigma[2] += eap[j * 3 - 2] * theta[j];
	sigma[3] += eap[j * 3 - 1] * theta[j];
/* L11: */
    }
} /* etc2p2c_ */

