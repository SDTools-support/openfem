/* fobase.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int fobase_(ndim, ndims, nno, npo, npi, vp, vp1, vdpq2, vdpq1, x, y, z__, vdpelc, xint, yint, zint, delta, dfinv, df, indic, fi, fn, pn, iopt)
int32 *ndim, *ndims, *nno, *npo, *npi;
doublereal *vp, *vp1, *vdpq2, *vdpq1, *x, *y, *z__, *vdpelc, *xint, *yint, *zint, *delta, *dfinv, *df;
int32 *indic;
doublereal *fi, *fn, *pn;
int32 *iopt;
{
    /* System generated locals */
    int32 vp_dim1, vp_offset, vp1_dim1, vp1_offset, vdpq2_dim1, vdpq2_dim2, vdpq2_offset, vdpq1_dim1, vdpq1_dim2, vdpq1_offset, vdpelc_dim1, vdpelc_dim2, vdpelc_offset, fi_dim1, fi_offset, fn_dim1, fn_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal norm[3];
    static int32 i__, j, l, n;
    static doublereal xnorm, pi;

/* ..................................................................... */
/* BUT : CALCULER  SUR L'ELEMENT COURANT, AU POINTS D'INTEGRATION */
/* ---      -  LES VALEURS DES GRADIENTS DES FONCTIONS DE BASE  (VDPELC) */
/*          -  LES VALEURS DES JACOBIENS                        (DELTA) */
/*          -  LES COORDONNEES DES POINTS D'INTEGRATION */
/* ..................................................................... */
/* PROGRAMMEUR : MARINA VIDRASCU  INRIA 2001 */
/* ..................................................................... */
/* PARAMETRES D'ENTREE : */
/* ------------------- */
/*  NDIM           :  DIMENSION */
/*  NDIMS          :  3 VOLUME */
/*                    2 SURFACE */
/*  NNO            :  NOMBRE DE NOEUDS */
/*  NPO            :  NOMBRE DE POINTS */
/*  NPI            :  NOMBRE DE POINTS D'INTEGRATION */
/*  VP             :  VALEUR DES POLYNOMES DE BASE (FONCTIONS) */
/*  VP1            :  VALEUR DES POLYNOMES DE BASE (GEOMETRIE - FT) */
/*  VDPQ2          :  DERIVEES DES POLYNOMES DE BASE (FONCTIONS) */
/*  VDPQ1          :  DERIVEES DES DES POLYNOMES DE BASE (GEOMETRIE - FT) */
/*  X,Y,Z          :  COORDONNEES DES POINTS */
/*  VDPELC         :  DERIVEES DES POLS DE BASE SUR ELEMENT  COURANT */
/*  INDIC          :  >0 et <= 2 CALCUL DES DERIVEES DES POLYNOMES */
/*                    = 2  CALCUL DES COORDONNEES DES POINTS D'INTEGRATION */
/*                    CALCUL DE DELTA toujours! */
/*                    = 3 Calcul d'efforts ou pressions aux */
/*                         points d'integration */
/*  IOPT           :  1 forces mortes */
/*                    2 pressions */
/*                    3 les deux */
/* PARAMETRES DE SORTIE : */
/* -------------------- */
/*  XINT,YINT,ZINT :  COORDONNEES DES POINTS D'INTEGRATION */
/*  DELTA          :  JACOBIEN */
/* --------------------------------------------------------------------- */
    /* Parameter adjustments */
    --pn;
    fn_dim1 = *ndim;
    fn_offset = fn_dim1 + 1;
    fn -= fn_offset;
    --z__;
    --y;
    --x;
    fi_dim1 = *ndim;
    fi_offset = fi_dim1 + 1;
    fi -= fi_offset;
    --delta;
    --zint;
    --yint;
    --xint;
    vdpelc_dim1 = *ndim;
    vdpelc_dim2 = *nno;
    vdpelc_offset = vdpelc_dim1 * (vdpelc_dim2 + 1) + 1;
    vdpelc -= vdpelc_offset;
    vdpq1_dim1 = *ndims;
    vdpq1_dim2 = *npo;
    vdpq1_offset = vdpq1_dim1 * (vdpq1_dim2 + 1) + 1;
    vdpq1 -= vdpq1_offset;
    vdpq2_dim1 = *ndim;
    vdpq2_dim2 = *nno;
    vdpq2_offset = vdpq2_dim1 * (vdpq2_dim2 + 1) + 1;
    vdpq2 -= vdpq2_offset;
    vp1_dim1 = *npo;
    vp1_offset = vp1_dim1 + 1;
    vp1 -= vp1_offset;
    vp_dim1 = *nno;
    vp_offset = vp_dim1 + 1;
    vp -= vp_offset;
    dfinv -= 4;
    df -= 4;

    /* Function Body */
    i__1 = *npi;
    for (l = 1; l <= i__1; ++l) {

	i__2 = *ndims;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = *ndim;
	    for (j = 1; j <= i__3; ++j) {
		df[i__ + j * 3] = 0.;
/* L21: */
	    }
	}

	i__3 = *ndims;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__2 = *npo;
	    for (n = 1; n <= i__2; ++n) {
		df[i__ + 3] += vdpq1[i__ + (n + l * vdpq1_dim2) * vdpq1_dim1] * x[n];
		df[i__ + 6] += vdpq1[i__ + (n + l * vdpq1_dim2) * vdpq1_dim1] * y[n];
		df[i__ + 9] += vdpq1[i__ + (n + l * vdpq1_dim2) * vdpq1_dim1] * z__[n];
/* L1: */
	    }
	}
	if (*ndims == 3) {
	    dfinv[4] = df[8] * df[12] - df[11] * df[9];
	    dfinv[8] = df[4] * df[12] - df[6] * df[10];
	    dfinv[12] = df[4] * df[8] - df[5] * df[7];
	    dfinv[7] = df[10] * df[9] - df[7] * df[12];
	    dfinv[5] = df[11] * df[6] - df[5] * df[12];
	    dfinv[10] = df[7] * df[11] - df[10] * df[8];
	    dfinv[6] = df[5] * df[9] - df[8] * df[6];
	    dfinv[11] = df[5] * df[10] - df[11] * df[4];
	    dfinv[9] = df[6] * df[7] - df[9] * df[4];

	    delta[l] = df[4] * dfinv[4] + df[5] * dfinv[7] + df[6] * dfinv[10];
	} else {

/* Computing 2nd power */
	    d__1 = df[4] * df[8] - df[7] * df[5];
/* Computing 2nd power */
	    d__2 = df[7] * df[11] - df[10] * df[8];
/* Computing 2nd power */
	    d__3 = df[10] * df[5] - df[4] * df[11];
	    delta[l] = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	}
	if (*indic > 0 && *indic <= 2) {

	    for (i__ = 1; i__ <= 3; ++i__) {
		i__2 = *nno;
		for (j = 1; j <= i__2; ++j) {
		    vdpelc[i__ + (j + l * vdpelc_dim2) * vdpelc_dim1] = dfinv[i__ + 3] * vdpq2[(j + l * vdpq2_dim2) * vdpq2_dim1 + 1] + dfinv[i__ + 6] * vdpq2[(j + l * vdpq2_dim2) * vdpq2_dim1 + 2] + dfinv[i__ + 9] * vdpq2[(j + l * vdpq2_dim2) * vdpq2_dim1 + 3];
/* L2: */
		}
	    }
	}

	if (*indic >= 3) {
	    if (*iopt == 1 || *iopt == 3) {

/*           forces mortes (interpoler) */

		i__2 = *nno;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    fi[l * fi_dim1 + 1] += vp[i__ + l * vp_dim1] * fn[i__ * fn_dim1 + 1];
		    fi[l * fi_dim1 + 2] += vp[i__ + l * vp_dim1] * fn[i__ * fn_dim1 + 2];
		    fi[l * fi_dim1 + 3] += vp[i__ + l * vp_dim1] * fn[i__ * fn_dim1 + 3];
/* L4: */
		}
	    }
	    if (*iopt == 2 || *iopt == 3) {

/*           forces de pression (interpoler) */

		pi = (float)0.;
		i__2 = *nno;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    pi += vp[i__ + l * vp_dim1] * pn[i__];
/* L5: */
		}
		norm[0] = df[7] * df[11] - df[10] * df[8];
		norm[1] = df[10] * df[5] - df[4] * df[11];
		norm[2] = df[4] * df[8] - df[7] * df[5];
/* Computing 2nd power */
		d__1 = norm[0];
/* Computing 2nd power */
		d__2 = norm[1];
/* Computing 2nd power */
		d__3 = norm[2];
		xnorm = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
		norm[0] /= xnorm;
		norm[1] /= xnorm;
		norm[2] /= xnorm;
		fi[l * fi_dim1 + 1] += pi * norm[0];
		fi[l * fi_dim1 + 2] += pi * norm[1];
		fi[l * fi_dim1 + 3] += pi * norm[2];
	    }
	}

	if (*indic == 2) {
	    xint[l] = 0.;
	    yint[l] = 0.;
	    zint[l] = 0.;
	    i__2 = *npo;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		xint[l] += vp1[i__ + l * vp1_dim1] * x[i__];
		yint[l] += vp1[i__ + l * vp1_dim1] * y[i__];
		zint[l] += vp1[i__ + l * vp1_dim1] * z__[i__];
/* L3: */
	    }
	}
/* L25: */
    }
} /* fobase_ */

