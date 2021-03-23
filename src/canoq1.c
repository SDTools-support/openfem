/* canoq1.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int canoq1_(ndim, npo, nopoi, coorp, xnorm, np)
int32 *ndim, *npo, *nopoi;
doublereal *coorp, *xnorm;
int32 *np;
{
    /* Initialized data */

    static int32 iq[8]	/* was [2][4] */ = { 2,4,3,1,4,2,1,3 };

    /* System generated locals */
    int32 coorp_dim1, coorp_offset, i__1;

    /* Local variables */
    static int32 i__, n;
    static doublereal x21, y21, x31, y31, z21, z31, nx, ny, nz;

/* ..................................................... */
/*     IN : */
/*     ---- */
/*      ndim : dimension */
/*      npo  : nombre de points de l'element (3 pour face triangulaire */
/*                                            4 face quadrangulaire) */
/*             --- > pour l'instant seul le cas 4 est traite */
/*             ......................................................... */
/*             pout triangle le calcul de nx,ny,nz est sorti de la boucle */
/*             et i=1 iq(1,i) = 2 iq(2,i) = 3 */
/*             .......................................................... */
/*      nopoi(npo)   : numeros des points de l'element (comme ds maillage) */
/*      coorp(npo,3) : coordonnees */
/*      np           : nombre total de points du maillage */
/*    OUT : */
/*    ----- */
/*     xnorm(1:3,n)  : somme(n(i)) au point n * aire */
/* ...................................................... */
    /* Parameter adjustments */
    coorp_dim1 = *npo;
    coorp_offset = coorp_dim1 + 1;
    coorp -= coorp_offset;
    --nopoi;
    xnorm -= 4;

    /* Function Body */
    i__1 = *npo;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*        n numero global du point */
	n = nopoi[i__];
/*        produit vectoriel des aretes adjacentes de i ds le bon sens! */
	x21 = coorp[iq[(i__ << 1) - 2] + coorp_dim1] - coorp[i__ + coorp_dim1];
	y21 = coorp[iq[(i__ << 1) - 2] + (coorp_dim1 << 1)] - coorp[i__ + (coorp_dim1 << 1)];
	x31 = coorp[iq[(i__ << 1) - 1] + coorp_dim1] - coorp[i__ + coorp_dim1];
	y31 = coorp[iq[(i__ << 1) - 1] + (coorp_dim1 << 1)] - coorp[i__ + (coorp_dim1 << 1)];
	z21 = coorp[iq[(i__ << 1) - 2] + coorp_dim1 * 3] - coorp[i__ + coorp_dim1 * 3];
	z31 = coorp[iq[(i__ << 1) - 1] + coorp_dim1 * 3] - coorp[i__ + coorp_dim1 * 3];
	nx = y21 * z31 - z21 * y31;
	ny = z21 * x31 - x21 * z31;
	nz = x21 * y31 - x31 * y21;
	xnorm[n * 3 + 1] += nx;
	xnorm[n * 3 + 2] += ny;
	xnorm[n * 3 + 3] += nz;
/* L1: */
    }
} /* canoq1_ */

