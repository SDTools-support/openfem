/* em3c2c.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static doublereal c_b2 = 0.;
static int32 c__0 = 0;
static int32 c__1 = 1;
static int32 c__3 = 3;

/* Subroutine */ int em3c2c_(nno, npo, x, y, z__, npi, ijt, poids, vp, vdpq1, ro, ae, delta)
int32 *nno, *npo;
doublereal *x, *y, *z__;
int32 *npi, *ijt;
doublereal *poids, *vp, *vdpq1, *ro, *ae, *delta;
{
    /* System generated locals */
    int32 vdpq1_dim2, vdpq1_offset, vp_dim1, vp_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static int32 iopt;
    static doublereal xint[1], yint[1], zint[1];
    static int32 i__, j, l;
    static doublereal dfinv[9]	/* was [3][3] */;
    extern /* Subroutine */ int dcopy_();
    static doublereal a2[1];
    static int32 k1, k2, k3;
    static doublereal df[9]	/* was [3][3] */, fi[1];
    static int32 ii, jj, nm, indice;
    extern /* Subroutine */ int fobase_();

/*  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  BUT :   CALCUL DE LA MATRICE ELEMENTAIRE DE MASSE */
/*  --- */
/*  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  PARAMETRES D'ENTREE : */
/*  ------------------- */
/*  NNO          : nombre de points (interpolation des fonctions) */
/*  NPO          : nombre de points (interpolation de la geometrie) */
/*  X,Y,Z       : COORDONNEES DES NOEUDS DE L ELEMENT */
/*  RO          : MASSE VOLUMIQUE CONSTANTE SUR L ELEMENT */
/*  IJT         : PERMUTATION */
/*  NPI         : NOMBRE DE POINTS DE LA FORMULE D INTEGRATION */
/*  POIDS       : POIDS DE LA FORMULE D INTEGRATION SUR LE CUBE */
/*  VDPQ1       : VALEURS DES DERIVES DES POLYNOMES DE BASE */
/*                AUX POINTS D INTEGRATION (geometrie) */
/*  VP          : VALEURS DES POLYNOMES DE BASE (Fonctions) */
/*                AUX POINTS D INTEGRATION */
/*  PARAMETRES RESULTATS : */
/*  -------------------- */
/*  AE          : MATRICE ELEMENTAIRE DE MASSE STOCKEE SYMETRIQUE */
/*  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  PROGRAMMEURS  : Marina Vidrascu INRIA septembre 2001 */
/*  .................................................................. */
/*     --- non utilises --- */

    /* Parameter adjustments */
    --z__;
    --y;
    --x;
    --delta;
    vdpq1_dim2 = *npo;
    vdpq1_offset = (vdpq1_dim2 + 1) * 3 + 1;
    vdpq1 -= vdpq1_offset;
    vp_dim1 = *nno;
    vp_offset = vp_dim1 + 1;
    vp -= vp_offset;
    --poids;
    --ijt;
    --ae;

    /* Function Body */
    nm = *nno * 3 * (*nno * 3 + 1) / 2;
    dcopy_(&nm, &c_b2, &c__0, &ae[1], &c__1);
    indice = 0;
    fobase_(&c__3, &c__3, nno, npo, npi, &vp[vp_offset], &vp[vp_offset], &vdpq1[vdpq1_offset], &vdpq1[vdpq1_offset], &x[1], &y[1], &z__[1], a2, xint, yint, zint, &delta[1], dfinv, df, &indice, fi, fi, fi, &iopt);

    i__1 = *nno;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {

/*          BLOC 11 */
/*          ------- */
/* Computing MIN */
	    i__3 = ijt[i__], i__4 = ijt[j];
	    ii = min(i__3,i__4);
/* Computing MAX */
	    i__3 = ijt[i__], i__4 = ijt[j];
	    jj = max(i__3,i__4);
	    k1 = jj * (jj - 1) / 2 + ii;

/*          BLOC 22 */
/*          ------- */
/* Computing MIN */
	    i__3 = ijt[*nno + i__], i__4 = ijt[*nno + j];
	    ii = min(i__3,i__4);
/* Computing MAX */
	    i__3 = ijt[*nno + i__], i__4 = ijt[*nno + j];
	    jj = max(i__3,i__4);
	    k2 = jj * (jj - 1) / 2 + ii;

/*          BLOC 33 */
/*          ------- */
/* Computing MIN */
	    i__3 = ijt[(*nno << 1) + i__], i__4 = ijt[(*nno << 1) + j];
	    ii = min(i__3,i__4);
/* Computing MAX */
	    i__3 = ijt[(*nno << 1) + i__], i__4 = ijt[(*nno << 1) + j];
	    jj = max(i__3,i__4);
	    k3 = jj * (jj - 1) / 2 + ii;
	    i__3 = *npi;
	    for (l = 1; l <= i__3; ++l) {
		ae[k1] += vp[j + l * vp_dim1] * vp[i__ + l * vp_dim1] * poids[l] * delta[l] * *ro;
		ae[k2] += vp[j + l * vp_dim1] * vp[i__ + l * vp_dim1] * poids[l] * delta[l] * *ro;
		ae[k3] += vp[j + l * vp_dim1] * vp[i__ + l * vp_dim1] * poids[l] * delta[l] * *ro;
/* L11: */
	    }
	}
    }
} /* em3c2c_ */

