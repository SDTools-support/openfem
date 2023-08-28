/* etsmitx.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static doublereal c_b2 = 0.;
static int32 c__0 = 0;
static int32 c__1 = 1;

/* Subroutine */ int etsmitx_(gr, gs, sur, suri, nno, npi, vdpq2, vp2, x, y, z__)
doublereal *gr, *gs, *sur, *suri;
int32 *nno, *npi;
doublereal *vdpq2, *vp2, *x, *y, *z__;
{
    /* System generated locals */
    int32 vdpq2_dim2, vdpq2_offset, vp2_dim1, vp2_offset, i__1, i__2;

    /* Local variables */
    static int32 i__, l;
    extern /* Subroutine */ int dcopy_();

/*  .................................................................. */
/*  .................................................................. */
    /* Parameter adjustments */
    --z__;
    --y;
    --x;
    sur -= 6;
    vp2_dim1 = *nno;
    vp2_offset = vp2_dim1 + 1;
    vp2 -= vp2_offset;
    vdpq2_dim2 = *nno;
    vdpq2_offset = (vdpq2_dim2 + 1 << 1) + 1;
    vdpq2 -= vdpq2_offset;
    suri -= 6;
    gs -= 4;
    gr -= 4;

    /* Function Body */
    i__1 = *npi * 3;
    dcopy_(&i__1, &c_b2, &c__0, &gr[4], &c__1);
    i__1 = *npi * 3;
    dcopy_(&i__1, &c_b2, &c__0, &gs[4], &c__1);
    i__1 = *npi * 5;
    dcopy_(&i__1, &c_b2, &c__0, &suri[6], &c__1);
    i__1 = *npi;
    for (l = 1; l <= i__1; ++l) {
	i__2 = *nno;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    gr[l * 3 + 1] += vdpq2[(i__ + l * vdpq2_dim2 << 1) + 1] * x[i__];
	    gr[l * 3 + 2] += vdpq2[(i__ + l * vdpq2_dim2 << 1) + 1] * y[i__];
	    gr[l * 3 + 3] += vdpq2[(i__ + l * vdpq2_dim2 << 1) + 1] * z__[i__];
	    gs[l * 3 + 1] += vdpq2[(i__ + l * vdpq2_dim2 << 1) + 2] * x[i__];
	    gs[l * 3 + 2] += vdpq2[(i__ + l * vdpq2_dim2 << 1) + 2] * y[i__];
	    gs[l * 3 + 3] += vdpq2[(i__ + l * vdpq2_dim2 << 1) + 2] * z__[i__];
	    suri[l * 5 + 1] += vp2[i__ + l * vp2_dim1] * sur[i__ * 5 + 1];
	    suri[l * 5 + 2] += vp2[i__ + l * vp2_dim1] * sur[i__ * 5 + 2];
	    suri[l * 5 + 3] += vp2[i__ + l * vp2_dim1] * sur[i__ * 5 + 3];
	    suri[l * 5 + 4] += vp2[i__ + l * vp2_dim1] * sur[i__ * 5 + 4];
	    suri[l * 5 + 5] += vp2[i__ + l * vp2_dim1] * sur[i__ * 5 + 5];
/* L1: */
	}
    }
} /* etsmitx_ */

