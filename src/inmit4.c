/* inmit4.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

doublereal inmit4_(j, i__, x, y)
int32 *j, *i__;
doublereal *x, *y;
{
    /* Initialized data */

    static doublereal yt1[2] = { -1.,1. };
    static doublereal xt2[2] = { -1.,1. };

    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    extern /* Subroutine */ int mexErrMsgTxt();

/* .............................................. */
/*     polynomes d'interpolation mit9 */

/*          j = 1  composantes  eps_rr, eps_rz */
/*          j = 2  composantes  eps_ss, eps_sz */
/*          j = 3  composante   eps_rs */
/* ............................................... */
/*     comp rr rz */
/*     comp ss  sz */
/*     comp rs */

/*     composantes eps_rr, eps_rz */
/*     -------------------------- */
    if (*j == 1) {
	ret_val = (*y / yt1[*i__ - 1] + 1) * (float).5;

/*     composantes eps_ss, eps_sz */
/*     -------------------------- */
    } else if (*j == 2) {
	ret_val = (*x / xt2[*i__ - 1] + 1) * (float).5;

/*     composante eps_rs */
/*     -------------------------- */
    } else if (*j == 3) {
	mexErrMsgTxt("eps rs : curieux!", 17L);
	return ret_val;
    }
    return ret_val;
} /* inmit4_ */

