/* ets5noe.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int ets5noe_(coor, fomega, fgamma, pressi, norefs, be)
doublereal *coor, *fomega, *fgamma, *pressi;
int32 *norefs;
doublereal *be;
{
    /* Initialized data */

    static int32 ind[4] = { 2,3,4,1 };
    static int32 np[12]	/* was [3][4] */ = { 1,2,5,2,3,5,3,4,5,4,1,5 };

    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal aret[8]	/* was [2][4] */, surf, f[10];
    static int32 iface, i__, j, k;
    static doublereal x[5], y[5];
    static int32 i1, j1, j2, j3, k1, n1, n2, n3, k2, i2;
    static doublereal x31, y31, x42, y42;
    static int32 nt;
    static doublereal arelon, xji, yji, xnu, ynu, div1;

/* ................................................................... */
/* but: CALCULE LE SECOND MEMBRE LOCAL SUR UN ELEMENT L PAR */
/*      BE(I)=SOMME(F(X)*PHII) */
/*      PUIS AJOUTE LES INTEGRALES DE surface */

/* in : coor(4,ndim) : coordonees des 4 sommets. */
/*      fomega(2,4)  : fx, fy volumiques AUX 4 sommets. */
/*      fgamma(ndim,2*nbarete): fx, fy aux noeuds de chaque arete. */
/*      pressi(2*nbarete)     : pression aux noeuds de chaque arete */
/*      norefs(nbarete,2): norefs(i,1) = 0 si fgamma   = 0 sur arete_i */
/*                         norefs(i,2) = 0 si pression = 0 sur arete_i */

/* out: BE(10)       : SECOND MEMBRE */
/* ------------------------------------------------------------------- */
    /* Parameter adjustments */
    --be;
    norefs -= 5;
    --pressi;
    fgamma -= 3;
    fomega -= 3;
    coor -= 6;

    /* Function Body */

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

/*     -- FORCES VOLUMIQUES */

    for (i__ = 1; i__ <= 10; ++i__) {
	be[i__] = (float)0.;
/* L2: */
    }
    for (i__ = 1; i__ <= 5; ++i__) {
	for (k = 1; k <= 2; ++k) {
	    f[(i__ << 1) - 2 + k - 1] = fomega[k + (i__ << 1)];
/* L3: */
	}
    }

/* --  FOMEGA: INTEGRATION AUX 4 SOMMETS */

    for (nt = 1; nt <= 4; ++nt) {
	n1 = np[nt * 3 - 3];
	n2 = np[nt * 3 - 2];
	n3 = np[nt * 3 - 1];
/*       -- CALCUL DE SURF */

	surf = ((x[n1 - 1] - x[n2 - 1]) * (y[n1 - 1] - y[n3 - 1]) - (x[n1 - 1] - x[n3 - 1]) * (y[n1 - 1] - y[n2 - 1])) / (float)24.;
/*       -- BOUCLE SUR LES DIRECTIONS */

	for (k = 1; k <= 2; ++k) {
	    j1 = (n1 - 1 << 1) + k;
	    j2 = (n2 - 1 << 1) + k;
	    j3 = (n3 - 1 << 1) + k;
	    be[j1] = surf * (f[j1 - 1] * (float)2. + f[j2 - 1] + f[j3 - 1]) + be[j1];
	    be[j2] = surf * (f[j2 - 1] * (float)2. + f[j3 - 1] + f[j1 - 1]) + be[j2];
	    be[j3] = surf * (f[j3 - 1] * (float)2. + f[j1 - 1] + f[j2 - 1]) + be[j3];
/* L4: */
	}
/* L5: */
    }

/* --  FGAMMA: INTEGRATION AU milieu des aretes. */
/*     fgamma(2,2*nbarete) fx, fy aux extremites de chaque arete */
/*     pressi(2*nbarete) pression aux extremites de chaque arete */
/*     aret(2,nbarete) fx, fy au milieu de chaque arete */

    for (k = 1; k <= 4; ++k) {
	xnu = 0.;
	ynu = 0.;
	if (norefs[k + 8] != 0) {
/*         -- LONGEUR ARETE, COSINUS DIRECTEURS DE LA NORMALE ext. */
	    j = k % 4 + 1;
	    xji = coor[j + 5] - coor[k + 5];
	    yji = coor[j + 10] - coor[k + 10];
/* Computing 2nd power */
	    d__1 = xji;
/* Computing 2nd power */
	    d__2 = yji;
	    arelon = sqrt(d__1 * d__1 + d__2 * d__2);
	    xnu = yji / arelon;
	    ynu = -xji / arelon;
	}
	aret[(k << 1) - 2] = (-(pressi[(k - 1 << 1) + 1] + pressi[(k - 1 << 1) + 2]) * xnu + (fgamma[((k - 1 << 1) + 1 << 1) + 1] + fgamma[((k - 1 << 1) + 2 << 1) + 1])) / 2.;
	aret[(k << 1) - 1] = (-(pressi[(k - 1 << 1) + 1] + pressi[(k - 1 << 1) + 2]) * ynu + (fgamma[((k - 1 << 1) + 1 << 1) + 2] + fgamma[((k - 1 << 1) + 2 << 1) + 2])) / 2.;
/* L6: */
    }
    for (iface = 1; iface <= 4; ++iface) {
	i1 = iface;
	i2 = ind[iface - 1];
/* Computing 2nd power */
	d__1 = x[i1 - 1] - x[i2 - 1];
/* Computing 2nd power */
	d__2 = y[i1 - 1] - y[i2 - 1];
	surf = sqrt(d__1 * d__1 + d__2 * d__2) * (float).5;
	j1 = (i1 << 1) - 1;
	j2 = (i2 << 1) - 1;
	k1 = i1 << 1;
	k2 = i2 << 1;
	be[j1] += surf * aret[(iface << 1) - 2];
	be[j2] += surf * aret[(iface << 1) - 2];
	be[k1] += surf * aret[(iface << 1) - 1];
	be[k2] += surf * aret[(iface << 1) - 1];
/* L7: */
    }
} /* ets5noe_ */

