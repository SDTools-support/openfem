/* fonmit.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static doublereal c_b2 = 0.;
static int32 c__0 = 0;
static int32 c__1 = 1;

/* Subroutine */ int fonmit_(npo, npt, x, y, z__, xnorm, t, vdtq2, vp2, f, fz, ind, gzz, tt)
int32 *npo, *npt;
doublereal *x, *y, *z__, *xnorm, *t, *vdtq2, *vp2, *f, *fz;
int32 *ind;
doublereal *gzz, *tt;
{
    /* System generated locals */
    int32 vdtq2_dim2, vdtq2_offset, vp2_dim1, vp2_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static int32 i__, l;
    extern /* Subroutine */ int dcopy_();

/* ............................................................ */

/*   calcul de F = g_*,* aux tying points F= F+fz*z */
/*                       ou aux points d'integration */
/*   npo        : nombre de points */
/*   npt        : somme sur les points d'interpolation */
/*   x,y,z      : coordonnes des noeuds de l'element courant */
/*   xnorm      : normale aux noeuds */
/*   t          : epaisseur de la coque (donee aux noeuds) */
/*   vdtq2(1,*) : derivee par rapport a r aux pt d'interpolation (tying) */
/*   vdtq2(2,*) : derivee par rapport a s ux pt d'interpolation */
/*   vp2        : valeur des polynomes de base aux pt d'interpolation */

/*  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     PROGRAMMEUR  :  Marina Vidrascu INRIA 2000 */
/*  .................................................................. */
    /* Parameter adjustments */
    --t;
    xnorm -= 4;
    --z__;
    --y;
    --x;
    --tt;
    --gzz;
    fz -= 13;
    f -= 13;
    vp2_dim1 = *npo;
    vp2_offset = vp2_dim1 + 1;
    vp2 -= vp2_offset;
    vdtq2_dim2 = *npo;
    vdtq2_offset = (vdtq2_dim2 + 1 << 1) + 1;
    vdtq2 -= vdtq2_offset;

    /* Function Body */
    i__1 = *npt * 9;
    dcopy_(&i__1, &c_b2, &c__0, &f[13], &c__1);
    i__1 = *npt * 9;
    dcopy_(&i__1, &c_b2, &c__0, &fz[13], &c__1);
    if (*ind == 1) {
	dcopy_(npt, &c_b2, &c__0, &gzz[1], &c__1);
	dcopy_(npt, &c_b2, &c__0, &tt[1], &c__1);
    }
    i__1 = *npt;
    for (l = 1; l <= i__1; ++l) {
	i__2 = *npo;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    f[(l * 3 + 1) * 3 + 1] += vdtq2[(i__ + l * vdtq2_dim2 << 1) + 1] * x[i__];
	    fz[(l * 3 + 1) * 3 + 1] += t[i__] * (float).5 * vdtq2[(i__ + l * vdtq2_dim2 << 1) + 1] * xnorm[i__ * 3 + 1];
	    f[(l * 3 + 2) * 3 + 1] += vdtq2[(i__ + l * vdtq2_dim2 << 1) + 1] * y[i__];
	    fz[(l * 3 + 2) * 3 + 1] += t[i__] * (float).5 * vdtq2[(i__ + l * vdtq2_dim2 << 1) + 1] * xnorm[i__ * 3 + 2];
	    f[(l * 3 + 3) * 3 + 1] += vdtq2[(i__ + l * vdtq2_dim2 << 1) + 1] * z__[i__];
	    fz[(l * 3 + 3) * 3 + 1] += t[i__] * (float).5 * vdtq2[(i__ + l * vdtq2_dim2 << 1) + 1] * xnorm[i__ * 3 + 3];
	    f[(l * 3 + 1) * 3 + 2] += vdtq2[(i__ + l * vdtq2_dim2 << 1) + 2] * x[i__];
	    fz[(l * 3 + 1) * 3 + 2] += t[i__] * (float).5 * vdtq2[(i__ + l * vdtq2_dim2 << 1) + 2] * xnorm[i__ * 3 + 1];
	    f[(l * 3 + 2) * 3 + 2] += vdtq2[(i__ + l * vdtq2_dim2 << 1) + 2] * y[i__];
	    fz[(l * 3 + 2) * 3 + 2] += t[i__] * (float).5 * vdtq2[(i__ + l * vdtq2_dim2 << 1) + 2] * xnorm[i__ * 3 + 2];
	    f[(l * 3 + 3) * 3 + 2] += vdtq2[(i__ + l * vdtq2_dim2 << 1) + 2] * z__[i__];
	    fz[(l * 3 + 3) * 3 + 2] += t[i__] * (float).5 * vdtq2[(i__ + l * vdtq2_dim2 << 1) + 2] * xnorm[i__ * 3 + 3];
	    f[(l * 3 + 1) * 3 + 3] += t[i__] * (float).5 * vp2[i__ + l * vp2_dim1] * xnorm[i__ * 3 + 1];
	    f[(l * 3 + 2) * 3 + 3] += t[i__] * (float).5 * vp2[i__ + l * vp2_dim1] * xnorm[i__ * 3 + 2];
	    f[(l * 3 + 3) * 3 + 3] += t[i__] * (float).5 * vp2[i__ + l * vp2_dim1] * xnorm[i__ * 3 + 3];
	    if (*ind == 1) {
/* Computing 2nd power */
		d__1 = t[i__];
		gzz[l] += d__1 * d__1 * (float).25 * vp2[i__ + l * vp2_dim1];
		tt[l] += t[i__] * vp2[i__ + l * vp2_dim1];
	    }
/* L1: */
	}
    }
} /* fonmit_ */

