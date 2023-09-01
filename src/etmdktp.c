/* etmdktp.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int etmdktp_(coor, car, iopt, ae)
doublereal *coor, *car;
int32 *iopt;
doublereal *ae;
{
    /* Initialized data */

    static doublereal q[45] = { .0960317,.0103175,.0015377,.0103175,9.4246e-4,.0015377,.0353175,.00674603,.00376984,.0960317,.00376984,6.44841e-4,5.45635e-4,.0103175,.0015377,.00674603,.00124008,6.44841e-4,.0103175,9.4246e-4,.0015377,.0353175,.00376984,.00674603,.0353175,.00674603,.00376984,.0960317,.00674603,6.44841e-4,.00124008,.00376984,6.44841e-4,5.45635e-4,.0103175,.0015377,.00376984,5.45635e-4,6.44841e-4,.00674603,.00124008,6.44841e-4,.0103175,9.4246e-4,.0015377 };

    /* System generated locals */
    int32 i__1;

    /* Local variables */
    static int32 i__, j, k, l, n;
    static doublereal r__[81]	/* was [9][9] */, delta, epais;
    static int32 kk;
    static doublereal x12, x31, x23, y23, y31, y12, rho, roh;

/* *************************************************************** */
/* but: calcul de la matrice de masse de l element de plaque */
/*      tria DKTP (DISCRETE KIRCHHOFF THEORY) */
/* in : coor(noe,ndim) : coordonees des 3 sommets. */
/*      car            : caracteristiques des materiaux */
/*                       car(1) = masse volumique */
/*                       car(2) = epaisseur de la plaque */
/*      iopt           : ouvert si masse lumping ou autre ds futur */

/* out: ae            : matrice triangulaire sup */

/* programmeur : modulef */
/* ............................................................... */
    /* Parameter adjustments */
    --ae;
    --car;
    coor -= 4;

    /* Function Body */


    x23 = coor[5] - coor[6];
    y23 = coor[8] - coor[9];
    x31 = coor[6] - coor[4];
    y31 = coor[9] - coor[7];
    x12 = coor[4] - coor[5];
    y12 = coor[7] - coor[8];
    delta = x31 * y12 - x12 * y31;
    for (i__ = 1; i__ <= 9; ++i__) {
	for (j = 1; j <= 9; ++j) {
	    r__[i__ + j * 9 - 10] = (float)0.;
/* L1: */
	}
    }
    r__[0] = 1.;
    r__[30] = 1.;
    r__[60] = 1.;
    r__[10] = -y12;
    r__[19] = x12;
    r__[11] = y31;
    r__[20] = -x31;
    r__[40] = -y23;
    r__[49] = x23;
    r__[41] = y12;
    r__[50] = -x12;
    r__[70] = -y31;
    r__[79] = x31;
    r__[71] = y23;
    r__[80] = -x23;
    if (*iopt == 3) {
	roh = car[1] * abs(delta);
    } else {
	rho = car[1];
	epais = car[2];
	roh = rho * epais * abs(delta);
    }
    n = 0;
    for (i__ = 1; i__ <= 9; ++i__) {
	i__1 = i__;
	for (j = 1; j <= i__1; ++j) {
	    ++n;
	    ae[n] = (float)0.;
	    for (k = 1; k <= 9; ++k) {
		for (l = 1; l <= 9; ++l) {
		    if (k < l) {
			kk = (l - 1) * l / 2 + k;
		    } else {
			kk = (k - 1) * k / 2 + l;
		    }
		    ae[n] += r__[k + i__ * 9 - 10] * q[kk - 1] * r__[l + j * 9 - 10];
/* L3: */
		}
	    }
	    ae[n] *= roh;
/* L4: */
	}
    }
} /* etmdktp_ */

