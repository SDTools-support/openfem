/* etc5noe.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int etc5noe_(coor, car, iopt, u, sigma)
doublereal *coor, *car;
int32 *iopt;
doublereal *u, *sigma;
{
    /* Initialized data */

    static int32 n[12]	/* was [3][4] */ = { 1,2,5,2,3,5,3,4,5,4,1,5 };

    static int32 i__, j, k;
    static doublereal s, x[5], y[5];
    extern /* Subroutine */ int mexErrMsgTxt();
    static doublereal e1;
    static int32 i1, i2;
    static doublereal young, b11, b12, d11[10], d12[10], d22[10];
    static int32 ii;
    static doublereal x31, y31;
    static int32 np;
    static doublereal x42, y42, rr, dsigma[30]	/* was [3][10] */, bi1[3], bi2[3], poison, div[10], div1;

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  BUT : CALCUL DU TABLEAU DES CONTRAINTES ELEMENTAIRES AU CENTRE DE */
/*  ---   L ELEMENT PAR LA FORMULE : */
/*             SIGMA(I)=2.*E1*(D(I)-.5*TR(D)*ID)+2.*RR*TR(D)*ID */
/*        POUR L ELEMENT QUADRILATERAL A 4 SOUS TRIANGLES EN CROIX */
/*              ---   VERSION PLANE  --- */
/*  N' EST ECRIT CI DESSOUS QUE LA PARTIE DE CALCUL DES CONTRAINTES */
/*  ELASTIQUES POUR UN MATERIAU HOMOGENE POUR L'OPTION NTHELA=0 */
/*  MANQUE: (1)   LE TRAITEMENT DES OPTIONS NTHELA DIFFERENTES DE ZERO */
/*          (2)   LE CALCUL POUR DES MATERIAUX NON HOMOGENES */

/* in : coor(noe,ndim) : coordones des 4 sommets. */
/*      car(2): caracteristiques des materiaux */
/*              car(1) = young */
/*              car(2) = poisson */
/*      iopt = 1 isotrope Deformations Planes */
/*           = 2 isotrope Contraintes  Planes */
/* out: sigma(3)  CONTRAINTES SXX,SYY,SXY */
/* .................................................................... */
    /* Parameter adjustments */
    --sigma;
    u -= 3;
    --car;
    coor -= 6;

    /* Function Body */

    young = car[1];
    poison = car[2];
    e1 = young * .5 / (poison + 1.);
    if (*iopt == 1) {
	rr = e1 * (float).5 / ((float)1. - poison * (float)2.);
    } else if (*iopt == 2) {
	rr = e1 * (float).5 * (poison + (float)1.) / ((float)1. - poison);
    } else {
	mexErrMsgTxt("CHOIX CONTRAINTES OU DEFORMATIONS PLANES ?", 42L);
	return 0;
    }

/*     CALCUL DES COORDONNEES DES NOEUDS */

    for (i__ = 1; i__ <= 4; ++i__) {
	x[i__ - 1] = coor[i__ + 5];
	y[i__ - 1] = coor[i__ + 10];
/* L1: */
    }
    x31 = x[2] - x[0];
    y31 = y[2] - y[0];
    x42 = x[3] - x[1];
    y42 = y[3] - y[1];
    div1 = (float)1. / (y31 * x42 - y42 * x31);
    x[4] = (x[0] * y31 * x42 - x[1] * y42 * x31 + (y[1] - y[0]) * x31 * x42) * div1;
    y[4] = (y[1] * y31 * x42 - y[0] * y42 * x31 - (x[1] - x[0]) * y31 * y42) * div1;

/*     INITIALISATION */

    div1 = (float)0.;
    for (i1 = 1; i1 <= 10; ++i1) {
	div[i1 - 1] = (float)0.;
	d11[i1 - 1] = (float)0.;
	d12[i1 - 1] = (float)0.;
	d22[i1 - 1] = (float)0.;
/* L100: */
    }

/*     BOUCLE SUR LES SOUS-ELEMENTS */

    for (np = 1; np <= 4; ++np) {
	i__ = n[np * 3 - 3];
	j = n[np * 3 - 2];
	k = n[np * 3 - 1];

/*        CALCUL DES GRADIENTS BIJ DES COORDONNEES BARYCENTRIQUES */
/*        ET DE S=2*AIRE DU SOUS ELEMENT */

	s = (x[i__ - 1] - x[j - 1]) * (y[i__ - 1] - y[k - 1]) - (x[i__ - 1] - x[k - 1]) * (y[i__ - 1] - y[j - 1]);
	bi1[0] = y[j - 1] - y[k - 1];
	bi1[1] = y[k - 1] - y[i__ - 1];
	bi1[2] = y[i__ - 1] - y[j - 1];
	bi2[0] = x[k - 1] - x[j - 1];
	bi2[1] = x[i__ - 1] - x[k - 1];
	bi2[2] = x[j - 1] - x[i__ - 1];

/*        ACTUALISATION DES DEFORMATIONS */
	div1 += s;
	for (ii = 1; ii <= 3; ++ii) {
	    i1 = (n[ii + np * 3 - 4] << 1) - 1;
	    i2 = i1 + 1;
	    b11 = bi1[ii - 1] * (float).5;
	    b12 = bi2[ii - 1] * (float).5;
	    d11[i1 - 1] += b11;
	    d12[i1 - 1] = b12 + d12[i1 - 1];
	    d22[i1 - 1] = -b11 + d22[i1 - 1];
	    div[i1 - 1] = b11 * (float)2. + div[i1 - 1];
	    d11[i2 - 1] = -b12 + d11[i2 - 1];
	    d12[i2 - 1] = b11 + d12[i2 - 1];
	    d22[i2 - 1] = b12 + d22[i2 - 1];
	    div[i2 - 1] += b12 * (float)2.;
/* L101: */
	}
/* L10: */
    }

/*     CALCUL DE SIGMA AU CENTRE (DANS L'ORDRE SXX,SYY,SXY) */
/*     SIGMA XX => DSIGMA(1,I) */
/*     SIGMA YY => DSIGMA(2,I) */
/*     SIGMA XY => DSIGMA(3,I) */

    div1 = (float)2. / div1;
    for (j = 1; j <= 10; ++j) {
	dsigma[j * 3 - 3] = (e1 * d11[j - 1] + rr * div[j - 1]) * div1;
	dsigma[j * 3 - 2] = e1 * d12[j - 1] * div1;
	dsigma[j * 3 - 1] = (e1 * d22[j - 1] + rr * div[j - 1]) * div1;
/* L102: */
    }
/*     [      11]   [1 4 7 10 13 16 19 22 25 28]   [          ] */
/*     [SIGMA 22] = [2 5 8 11 14 17 20 23 26 29] * [u_solution] */
/*     [      12]   [3 6 9 12 15 18 21 24 27 30]   [          ] */
/*             3*1                3*9                     10*1 */

    sigma[1] = dsigma[0] * u[3] + dsigma[3] * u[4] + dsigma[6] * u[5] + dsigma[9] * u[6] + dsigma[12] * u[7] + dsigma[15] * u[8] + dsigma[18] * u[9] + dsigma[21] * u[10] + dsigma[24] * u[11] + dsigma[27] * u[12];
    sigma[2] = dsigma[1] * u[3] + dsigma[4] * u[4] + dsigma[7] * u[5] + dsigma[10] * u[6] + dsigma[13] * u[7] + dsigma[16] * u[8] + dsigma[19] * u[9] + dsigma[22] * u[10] + dsigma[25] * u[11] + dsigma[28] * u[12];
    sigma[3] = dsigma[2] * u[3] + dsigma[5] * u[4] + dsigma[8] * u[5] + dsigma[11] * u[6] + dsigma[14] * u[7] + dsigma[17] * u[8] + dsigma[20] * u[9] + dsigma[23] * u[10] + dsigma[26] * u[11] + dsigma[29] * u[12];

} /* etc5noe_ */

