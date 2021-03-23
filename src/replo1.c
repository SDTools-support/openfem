/* replo1.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int replo1_(npo, xnorm, v1, v2, ichoix, ierr, ierr_len)
int32 *npo;
doublereal *xnorm, *v1, *v2;
int32 *ichoix;
char *ierr;
ftnlen ierr_len;
{
    /* Initialized data */

    static int32 ip[3] = { 2,3,1 };
    static int32 ipas = 0;

    /* System generated locals */
    int32 i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal e[3];
    static int32 i__, n;
    static doublereal delta;
    extern /* Subroutine */ int mexPrintf();

/* ....................................................... */
/*     npo   : nombre de points a considerer */
/*     xnorm : tableau des normales pour  les points de l'element */
/*             xnorm(3,npo) */
/*     v1,v2 : repere local, rempli element/element */
/*             attention el isoparametrique noeuds et points coincident! */
/* ...................................................... */

    /* Parameter adjustments */
    v2 -= 4;
    v1 -= 4;
    xnorm -= 4;

    /* Function Body */
    i__1 = *npo;
    for (n = 1; n <= i__1; ++n) {
	if (*ichoix >= 0) {
/*            --------------------------------------------------- */

/*            choix regularise */

/*            --------------------------------------------------- */
/*            ^ vectoriel */
/*            e= x au depart */

	    e[0] = (float)1.;
	    e[1] = (float)0.;
	    e[2] = (float)0.;
	} else if (*ichoix == -1) {

/*            e = x */

	    e[0] = (float)1.;
	    e[1] = (float)0.;
	    e[2] = (float)0.;
	} else if (*ichoix == -2) {

/*            e = y */

	    e[0] = (float)0.;
	    e[1] = (float)1.;
	    e[2] = (float)0.;
	} else if (*ichoix == -3) {

/*            e = z */

	    e[0] = (float)0.;
	    e[1] = (float)0.;
	    e[2] = (float)1.;
	}

L102:
	v1[n * 3 + 1] = xnorm[n * 3 + 2] * e[2] - xnorm[n * 3 + 3] * e[1];
	v1[n * 3 + 2] = -xnorm[n * 3 + 1] * e[2] + xnorm[n * 3 + 3] * e[0];
	v1[n * 3 + 3] = xnorm[n * 3 + 1] * e[1] - xnorm[n * 3 + 2] * e[0];

/* Computing 2nd power */
	d__1 = v1[n * 3 + 1];
/* Computing 2nd power */
	d__2 = v1[n * 3 + 2];
/* Computing 2nd power */
	d__3 = v1[n * 3 + 3];
	delta = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);

/*         verifier que v1 n est pas nul (e//n) */

	if (delta <= (float).001) {
/*         mise en commentaire suite a l'utilisation de mexWarnMsgTxt */
/*             if (ipas.eq.0) print 100,n,e(1),e(2),e(3) */
/* 100         format( 'attention noeud local : ',i2,' direction ', */
/*     *       3(e8.2,1x), 'parallele a la normale ') */
/* L100: */
	    mexPrintf("attention direction parallele normale", 37L);
	    for (i__ = 1; i__ <= 3; ++i__) {
		if (e[i__ - 1] != (float)0.) {
		    e[i__ - 1] = 0.;
		    e[ip[i__ - 1] - 1] = 1.;
/*         mise en commentaire suite a l'utilisation de mexWarnMsgTxt */
/*                     if (ipas.eq.0) print 110, ip(i) */
/* 110                format(' nouvel axe privilegie e_',i1) */
/* L110: */
		    mexPrintf(" nouvel axe privilegie e_", 25L);
		    goto L102;
		}
/* L101: */
	    }
	    ++ipas;
	}
	v1[n * 3 + 1] /= delta;
	v1[n * 3 + 2] /= delta;
	v1[n * 3 + 3] /= delta;
/*         --------------------------------------------------- */

/*         calculer v2 = n^v1 */

	v2[n * 3 + 1] = xnorm[n * 3 + 2] * v1[n * 3 + 3] - xnorm[n * 3 + 3] * v1[n * 3 + 2];
	v2[n * 3 + 2] = xnorm[n * 3 + 3] * v1[n * 3 + 1] - xnorm[n * 3 + 1] * v1[n * 3 + 3];
	v2[n * 3 + 3] = xnorm[n * 3 + 1] * v1[n * 3 + 2] - xnorm[n * 3 + 2] * v1[n * 3 + 1];
/* Computing 2nd power */
	d__1 = v2[n * 3 + 1];
/* Computing 2nd power */
	d__2 = v2[n * 3 + 2];
/* Computing 2nd power */
	d__3 = v2[n * 3 + 3];
	delta = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);

	v2[n * 3 + 1] /= delta;
	v2[n * 3 + 2] /= delta;
	v2[n * 3 + 3] /= delta;
/* L1: */
    }
} /* replo1_ */

