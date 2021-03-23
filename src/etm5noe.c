/* etm5noe.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int etm5noe_(coor, car, iopt, ae)
doublereal *coor, *car;
int32 *iopt;
doublereal *ae;
{
    /* Initialized data */

    static int32 n[12]	/* was [3][4] */ = { 1,2,5,2,3,5,3,4,5,4,1,5 };

    static int32 i__, j, k;
    static doublereal s, x[5], y[5];
    static int32 ij, ni, nj;
    static doublereal x31, y31;
    static int32 np;
    static doublereal x42, y42, ss, rho, div1;

/* *************************************************************** */
/* but : calcul de la matrice de masse de l element QUAD 5NOE */
/* ---   4 SOUS TRIANGLES EN CROIX */
/* in : coor(noe,ndim) : coordones des 4 sommets. */
/*      car            : car(1) = rho masse volumique */
/*      iopt           : ouvert si masse lumping ou autre ds futur */

/* out: matrice de rigidite ae(55) DEMI-MATRICE SUPERIEURE */

/*  programmeur : modulef */
/* ................................................................. */
    /* Parameter adjustments */
    --ae;
    --car;
    coor -= 6;

    /* Function Body */
    rho = car[1];

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

    for (ij = 1; ij <= 55; ++ij) {
	ae[ij] = (float)0.;
/* L2: */
    }

/*     BOUCLE SUR LES SOUS-ELEMENTS */

    for (np = 1; np <= 4; ++np) {
	i__ = n[np * 3 - 3];
	j = n[np * 3 - 2];
	k = n[np * 3 - 1];
/*        ET DE S=2*AIRE DU SOUS ELEMENT */
	s = (x[i__ - 1] - x[j - 1]) * (y[i__ - 1] - y[k - 1]) - (x[i__ - 1] - x[k - 1]) * (y[i__ - 1] - y[j - 1]);
	ss = s * rho / (float)24.;
/*        ACTUALISATION DE LA MATRICE */
	for (k = 1; k <= 2; ++k) {
	    for (i__ = 1; i__ <= 3; ++i__) {
		ni = (n[i__ + np * 3 - 4] << 1) - 2 + k;
		for (j = 1; j <= 3; ++j) {
		    nj = (n[j + np * 3 - 4] << 1) - 2 + k;
		    if (nj <= ni) {
			ij = (ni - 1) * ni / 2 + nj;
			ae[ij] += ss;
			if (nj == ni) {
			    ae[ij] += ss;
			}
		    }
/* L103: */
		}
/* L102: */
	    }
/* L101: */
	}
/* L10: */
    }
} /* etm5noe_ */

