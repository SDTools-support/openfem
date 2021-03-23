/* pcq12d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

doublereal pcq12d_(j, i__, x, y)
int32 *j, *i__;
doublereal *x, *y;
{
    /* Initialized data */

    static doublereal xref[4] = { -1.,1.,1.,-1. };
    static doublereal yref[4] = { -1.,-1.,1.,1. };

    /* System generated locals */
    doublereal ret_val;

/* ......................................................... */
/*     polynomes de base du quadrilatere q2 a 9 noeuds */
/*            j = 0 : fonction */
/*            j = 1 : derivee par rapport a x */
/*            j = 2 : derivee par rapport a y */
/*     element de reference [-1,1]*[-1,1] */

/* ......................................................... */

/*     la fonction */
/*     ------------ */
    if (*j == 0) {
	ret_val = (xref[*i__ - 1] * *x + 1) * (float).25 * (yref[*i__ - 1] * *y + 1);

/*      la derivee par rapport a x */
/*      -------------------------- */
    } else if (*j == 1) {
	ret_val = xref[*i__ - 1] * (float).25 * (yref[*i__ - 1] * *y + 1);

/*      la derivee par rapport a y */
/*      -------------------------- */
    } else if (*j == 2) {
	ret_val = yref[*i__ - 1] * (float).25 * (xref[*i__ - 1] * *x + 1);
    }
    return ret_val;
} /* pcq12d_ */

