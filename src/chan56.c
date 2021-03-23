/* chan56.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int chan56_(nno, ae5, lae5, ae6, lae6, t, v1, v2, xnorm, e)
int32 *nno;
doublereal *ae5;
int32 *lae5;
doublereal *ae6;
int32 *lae6;
doublereal *t, *v1, *v2, *xnorm, *e;
{
    /* System generated locals */
    int32 i__1, i__2, i__3;

    /* Local variables */
    static int32 i__, j;
    static doublereal k;
    static int32 m, n, jlige, jcole, jligs, jcols, je, js, je1, je2, js1, js2, js3;
    static doublereal ae2221;

/*     ............................................................ */

/*     transformer la matrice de rigidite de 5 a 6 dl/noeud */
/*     transformer bloc par bloc, pour chacun */

/*     K6 = P_t K5 P */

/*     avec */

/*           K11 K12 0                       I   0   0 */
/*     k5 =  K21 K22 0                       0   P22 P23 */
/*           0    0  k                       0   P32 P33 */

/*     ............................................................. */
/*    --------------- boucle sur les noeuds   ----------------- */
/*        ----- lignes */
    /* Parameter adjustments */
    xnorm -= 4;
    v2 -= 4;
    v1 -= 4;
    --t;
    --ae5;
    --ae6;

    /* Function Body */
    i__1 = *nno;
    for (n = 1; n <= i__1; ++n) {
	jlige = (n - 1) * 5;
	jligs = (n - 1) * 6;
	i__2 = *nno;
	for (m = n; m <= i__2; ++m) {
	    jcole = (m - 1) * 5;
	    jcols = (m - 1) * 6;

/*        multiplication bloc / bloc */

/*        sous-blocs diagonaux */
	    if (m == n) {
/*           -- > K6_11  bloc 3*3  sym */
/*                K6_11  = K5_11 */
		for (j = 1; j <= 3; ++j) {
		    js = (jcols + j) * (jcols + j - 1) / 2;
		    je = (jcole + j) * (jcole + j - 1) / 2;
		    i__3 = j;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			ae6[js + jligs + i__] = ae5[je + jlige + i__];
/* L2: */
		    }
		}
	    } else {
/*           -- > K6_11  bloc 3*3  nsym */
		for (j = 1; j <= 3; ++j) {
		    js = (jcols + j) * (jcols + j - 1) / 2;
		    je = (jcole + j) * (jcole + j - 1) / 2;
		    for (i__ = 1; i__ <= 3; ++i__) {
			ae6[js + jligs + i__] = ae5[je + jlige + i__];
/* L20: */
		    }
		}
	    }
/*        -- > K6_22  bloc 2,2 */

	    je1 = (jcole + 4) * (jcole + 3) / 2;
	    je2 = (jcole + 5) * (jcole + 4) / 2;
	    js1 = (jcols + 4) * (jcols + 3) / 2;
	    js2 = (jcols + 5) * (jcols + 4) / 2;
	    if (m != n) {
		ae2221 = ae5[je1 + jlige + 5];
		k = 0.;
	    } else {
		ae2221 = ae5[je2 + jlige + 4];
		k = *e * 1e5;
	    }
	    ae6[js1 + jligs + 4] = v1[n * 3 + 1] * ae5[je1 + jlige + 4] * v1[m * 3 + 1] + v1[n * 3 + 1] * ae5[je2 + jlige + 4] * v2[m * 3 + 1] + v2[n * 3 + 1] * ae2221 * v1[m * 3 + 1] + v2[n * 3 + 1] * ae5[je2 + jlige + 5] * v2[m * 3 + 1] + xnorm[n * 3 + 1] * k * xnorm[m * 3 + 1];
	    if (m != n) {
		ae6[js1 + jligs + 5] = v1[n * 3 + 2] * ae5[je1 + jlige + 4] * v1[m * 3 + 1] + v1[n * 3 + 2] * ae5[je2 + jlige + 4] * v2[m * 3 + 1] + v2[n * 3 + 2] * ae2221 * v1[m * 3 + 1] + v2[n * 3 + 2] * ae5[je2 + jlige + 5] * v2[m * 3 + 1] + xnorm[n * 3 + 2] * k * xnorm[m * 3 + 1];
	    }
	    ae6[js2 + jligs + 4] = v1[n * 3 + 1] * ae5[je1 + jlige + 4] * v1[m * 3 + 2] + v1[n * 3 + 1] * ae5[je2 + jlige + 4] * v2[m * 3 + 2] + v2[n * 3 + 1] * ae2221 * v1[m * 3 + 2] + v2[n * 3 + 1] * ae5[je2 + jlige + 5] * v2[m * 3 + 2] + xnorm[n * 3 + 1] * k * xnorm[m * 3 + 2];
	    ae6[js2 + jligs + 5] = v1[n * 3 + 2] * ae5[je1 + jlige + 4] * v1[m * 3 + 2] + v1[n * 3 + 2] * ae5[je2 + jlige + 4] * v2[m * 3 + 2] + v2[n * 3 + 2] * ae2221 * v1[m * 3 + 2] + v2[n * 3 + 2] * ae5[je2 + jlige + 5] * v2[m * 3 + 2] + xnorm[n * 3 + 2] * k * xnorm[m * 3 + 2];
/*        -- > K6_33  bloc 1,1 */

	    je1 = (jcole + 4) * (jcole + 3) / 2;
	    je2 = (jcole + 5) * (jcole + 4) / 2;
	    js = (jcols + 6) * (jcols + 5) / 2;
	    ae6[js + jligs + 6] = v1[n * 3 + 3] * ae5[je1 + jlige + 4] * v1[m * 3 + 3] + v1[n * 3 + 3] * ae5[je2 + jlige + 4] * v2[m * 3 + 3] + v2[n * 3 + 3] * ae2221 * v1[m * 3 + 3] + v2[n * 3 + 3] * ae5[je2 + jlige + 5] * v2[m * 3 + 3] + xnorm[n * 3 + 3] * k * xnorm[m * 3 + 3];

/*        sous-blocs non-diagonaux */

/*        -- >K6_12 bloc 3,2 et K6_13 bloc 3,1 */
	    je1 = (jcole + 4) * (jcole + 3) / 2;
	    je2 = (jcole + 5) * (jcole + 4) / 2;
	    js1 = (jcols + 4) * (jcols + 3) / 2;
	    js2 = (jcols + 5) * (jcols + 4) / 2;
	    js3 = (jcols + 6) * (jcols + 5) / 2;
	    for (i__ = 1; i__ <= 3; ++i__) {
		ae6[js1 + jligs + i__] = ae5[je1 + jlige + i__] * v1[m * 3 + 1] + ae5[je2 + jlige + i__] * v2[m * 3 + 1];
		ae6[js2 + jligs + i__] = ae5[je1 + jlige + i__] * v1[m * 3 + 2] + ae5[je2 + jlige + i__] * v2[m * 3 + 2];
		ae6[js3 + jligs + i__] = ae5[je1 + jlige + i__] * v1[m * 3 + 3] + ae5[je2 + jlige + i__] * v2[m * 3 + 3];
/* L3: */
	    }
/*        -- >K6_23 bloc 2,1 */
	    je1 = (jcole + 4) * (jcole + 3) / 2;
	    je2 = (jcole + 5) * (jcole + 4) / 2;
	    js = (jcols + 6) * (jcols + 5) / 2;
	    ae6[js + jligs + 4] = v1[n * 3 + 1] * ae5[je1 + jlige + 4] * v1[m * 3 + 3] + v2[n * 3 + 1] * ae2221 * v1[m * 3 + 3] + v1[n * 3 + 1] * ae5[je2 + jlige + 4] * v2[m * 3 + 3] + v2[n * 3 + 1] * ae5[je2 + jlige + 5] * v2[m * 3 + 3] + xnorm[n * 3 + 1] * k * xnorm[m * 3 + 3];
	    ae6[js + jligs + 5] = v1[n * 3 + 2] * ae5[je1 + jlige + 4] * v1[m * 3 + 3] + v2[n * 3 + 2] * ae2221 * v1[m * 3 + 3] + v1[n * 3 + 2] * ae5[je2 + jlige + 4] * v2[m * 3 + 3] + v2[n * 3 + 2] * ae5[je2 + jlige + 5] * v2[m * 3 + 3] + xnorm[n * 3 + 2] * k * xnorm[m * 3 + 3];
/*        partie triangulaire inferieure si m different de n */
	    if (m != n) {
/*           -- >K6_21 bloc 2,3 */
/*           -- >K6_31 bloc 1,3 */
		for (j = 1; j <= 3; ++j) {
		    js1 = (jcols + j) * (jcols + j - 1) / 2;
		    je1 = (jcole + j) * (jcole + j - 1) / 2;
		    ae6[js1 + jligs + 4] = ae5[je1 + jlige + 4] * v1[n * 3 + 1] + ae5[je1 + jlige + 5] * v2[n * 3 + 1];
		    ae6[js1 + jligs + 5] = ae5[je1 + jlige + 4] * v1[n * 3 + 2] + ae5[je1 + jlige + 5] * v2[n * 3 + 2];
		    ae6[js1 + jligs + 6] = ae5[je1 + jlige + 4] * v1[n * 3 + 3] + ae5[je1 + jlige + 5] * v2[n * 3 + 3];
/* L4: */
		}
/*           -- >K6_32 bloc 1,2 */
		je1 = (jcole + 4) * (jcole + 3) / 2;
		je2 = (jcole + 5) * (jcole + 4) / 2;
		js1 = (jcols + 4) * (jcols + 3) / 2;
		js2 = (jcols + 5) * (jcols + 4) / 2;
		ae6[js1 + jligs + 6] = v1[n * 3 + 3] * ae5[je1 + jlige + 4] * v1[m * 3 + 1] + v1[n * 3 + 3] * ae5[je2 + jlige + 4] * v2[m * 3 + 1] + v2[n * 3 + 3] * ae5[je1 + jlige + 5] * v1[m * 3 + 1] + v2[n * 3 + 3] * ae5[je2 + jlige + 5] * v2[m * 3 + 1] + xnorm[n * 3 + 3] * k * xnorm[m * 3 + 1];
		ae6[js2 + jligs + 6] = v1[n * 3 + 3] * ae5[je1 + jlige + 4] * v1[m * 3 + 2] + v1[n * 3 + 3] * ae5[je2 + jlige + 4] * v2[m * 3 + 2] + v2[n * 3 + 3] * ae5[je1 + jlige + 5] * v1[m * 3 + 2] + v2[n * 3 + 3] * ae5[je2 + jlige + 5] * v2[m * 3 + 2] + xnorm[n * 3 + 3] * k * xnorm[m * 3 + 2];
	    }
/* L1: */
	}
    }
} /* chan56_ */

