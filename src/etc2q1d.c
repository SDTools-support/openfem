/* etc2q1d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int etc2q1d_(coor, car, iopt, u, alpha, theta, sigma)
doublereal *coor, *car;
int32 *iopt;
doublereal *u, *alpha, *theta, *sigma;
{
    static doublereal c__, e[6];
    static int32 i__, j;
    static doublereal delta, young, unmnu, x21, y21, x41, y41, x32, y32, x42, y42, x43, y43, x31, y31, dsigma[24], ed1[3], eap[12]	/* was [3][4] */, poisson;

/* *************************************************************** */
/* BUT: CALCUL DES CONTRAINTES DE L ELEMENT QUAD 2q1d */
/* --- */
/* in : coor(noe,ndim) : coor. 3 sommets */
/*      car, iopt      : caracteristiques des materiaux */
/*      U(ndim,noe): deplacements U_x et U_y aux 3 sommets */
/*      alpha(3)   : tenseur de dilatation thermique */
/*      theta(4)   : temperature aux 3 somets */
/* out: SIGMA(3)   :  S_xx, S_yy, S_xy elastiques au barycentre */

/* programmeur : modulef */
/* ............................................................... */

    /* Parameter adjustments */
    --sigma;
    --theta;
    --alpha;
    u -= 3;
    --car;
    coor -= 5;

    /* Function Body */
    x21 = coor[6] - coor[5];
    y21 = coor[10] - coor[9];
    x31 = coor[7] - coor[5];
    y31 = coor[11] - coor[9];
    x32 = coor[7] - coor[6];
    y32 = coor[11] - coor[10];
    x41 = coor[8] - coor[5];
    y41 = coor[12] - coor[9];
    x42 = coor[8] - coor[6];
    y42 = coor[12] - coor[10];
    x43 = coor[8] - coor[7];
    y43 = coor[12] - coor[11];

    delta = ((x21 - x43) * (y41 + y32) - (y21 - y43) * (x41 + x32)) / 2.;

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

    for (j = 1; j <= 6; ++j) {
	e[j - 1] /= delta;
/* L2: */
    }

/*     [      11]   [E1 E2 E4]    [    ]    [          ] */
/*     [SIGMA 22] = [   E3 E5]  * [ DP ]  * [u_solution] */
/*     [      12]   [      E6]    [    ]    [          ] */
/*             3*1          3*3       3*8             8*1 */
/*     --- */
/*     [      11]   [1 4 7 10 13 16 19 22]     [          ] */
/*     [SIGMA 22] = [2 5 8 11 14 17 20 23]   * [u_solution] */
/*     [      12]   [3 6 9 12 15 18 21 24]     [          ] */
/*             3*1                      3*8              8*1 */

/*     ---    SIGMA11 */
    dsigma[0] = -e[0] * y42 + e[3] * x42;
    dsigma[6] = e[0] * y31 - e[3] * x31;
    dsigma[12] = e[0] * y42 - e[3] * x42;
    dsigma[18] = -e[0] * y31 + e[3] * x31;
    dsigma[3] = e[1] * x42 - e[3] * y42;
    dsigma[9] = -e[1] * x31 + e[3] * y31;
    dsigma[15] = -e[1] * x42 + e[3] * y42;
    dsigma[21] = e[1] * x31 - e[3] * y31;
    sigma[1] = dsigma[0] * u[3] + dsigma[3] * u[4] + dsigma[6] * u[5] + dsigma[9] * u[6] + dsigma[12] * u[7] + dsigma[15] * u[8] + dsigma[18] * u[9] + dsigma[21] * u[10];
/*     ---    SIGMA22 */
    dsigma[1] = -e[1] * y42 + e[4] * x42;
    dsigma[7] = e[1] * y31 - e[4] * x31;
    dsigma[13] = e[1] * y42 - e[4] * x42;
    dsigma[19] = -e[1] * y31 + e[4] * x31;
    dsigma[4] = e[2] * x42 - e[4] * y42;
    dsigma[10] = -e[2] * x31 + e[4] * y31;
    dsigma[16] = -e[2] * x42 + e[4] * y42;
    dsigma[22] = e[2] * x31 - e[4] * y31;
    sigma[2] = dsigma[1] * u[3] + dsigma[4] * u[4] + dsigma[7] * u[5] + dsigma[10] * u[6] + dsigma[13] * u[7] + dsigma[16] * u[8] + dsigma[19] * u[9] + dsigma[22] * u[10];
/*     ---    SIGMA12 */
    dsigma[2] = -e[3] * y42 + e[5] * x42;
    dsigma[8] = e[3] * y31 - e[5] * x31;
    dsigma[14] = e[3] * y42 - e[5] * x42;
    dsigma[20] = -e[3] * y31 + e[5] * x31;
    dsigma[5] = e[4] * x42 - e[5] * y42;
    dsigma[11] = -e[4] * x31 + e[5] * y31;
    dsigma[17] = -e[4] * x42 + e[5] * y42;
    dsigma[23] = e[4] * x31 - e[5] * y31;
    sigma[3] = dsigma[2] * u[3] + dsigma[5] * u[4] + dsigma[8] * u[5] + dsigma[11] * u[6] + dsigma[14] * u[7] + dsigma[17] * u[8] + dsigma[20] * u[9] + dsigma[23] * u[10];
/*      print *,'---------------- Verif avec impressions Modulef' */
/*      print *, (dsigma(i), i= 1, 6) */
/*      print *, (dsigma(i), i= 7, 12) */
/*      print *, (dsigma(i), i= 13, 18) */
/*      print *, (dsigma(i), i= 19, 24) */

/*     CONTRAINTES THERMIQUES = HOOKE * ALPHA * TETA */


/* [s_11]   [E1 E2 E4]    [alpha_x ]                   [teta1] */
/* [s_22] = [   E3 E5]  * [alpha_y ]  *[ p1, p2, p3] * [teta2] */
/* [s 12]   [      E6]    [alpha_xy]              1*4  [teta3] */
/*     3*1          3*3           3*1                  [teta4] */

/*     SIGMA(TETA) = - (E) * (ALPHA) * (P) (X,Y) */

/*     [ED1] = [E] * [ALPHA] */
/*              3*3       3*1 */

    ed1[0] = -(e[0] * alpha[1] + e[1] * alpha[3] + e[3] * alpha[2]);
    ed1[1] = -(e[1] * alpha[1] + e[2] * alpha[3] + e[4] * alpha[2]);
    ed1[2] = -(e[3] * alpha[1] + e[4] * alpha[3] + e[5] * alpha[2]);

/*     [ED1] * [P] => [EAP] */
/*        3*1   1*4      3*4 */


/*     ED1 * P => EAP */

    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 4; ++j) {
	    eap[i__ + j * 3 - 4] = ed1[i__ - 1] / (float)4.;
/* L3: */
	}
    }

    for (j = 1; j <= 4; ++j) {
	sigma[1] += eap[j * 3 - 3] * theta[j];
	sigma[2] += eap[j * 3 - 2] * theta[j];
	sigma[3] += eap[j * 3 - 1] * theta[j];
/* L4: */
    }

    return 0;
} /* etc2q1d_ */

