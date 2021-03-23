/* bremit.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static doublereal c_b2 = 0.;
static int32 c__0 = 0;
static int32 c__1 = 1;
static int32 c__2 = 2;
static int32 c__3 = 3;

/* Subroutine */ int bremit_(nno, npt1, npt2, npt3, npt, xnorm, v1, v2, t, npi, xint, yint, vdtq2, vp2, f, fz, inmitc, b, bz, bzz, irint, fi, fiz, vdpq2, vp2i, l)
int32 *nno, *npt1, *npt2, *npt3, *npt;
doublereal *xnorm, *v1, *v2, *t;
int32 *npi;
doublereal *xint, *yint, *vdtq2, *vp2, *f, *fz;
doublereal (*inmitc) ();
doublereal *b, *bz, *bzz;
int32 *irint;
doublereal *fi, *fiz, *vdpq2, *vp2i;
int32 *l;
{
    /* System generated locals */
    int32 vdtq2_dim2, vdtq2_offset, vp2_dim1, vp2_offset, vdpq2_dim2, vdpq2_offset, vp2i_dim1, vp2i_offset, i__1, i__2, i__3;

    /* Local variables */
    static doublereal valx, valy;
    static int32 n;
    static doublereal valdx, valdy;
    extern /* Subroutine */ int dcopy_();
    static int32 lt;

/* ................................................................... */

/*   calcul de B aux points d'integration comme polynome de degre 2 */
/*       B(i) = b(i)*z*bz(i)+z**2*bzz(i) */

/*  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*  PARAMETRES D'ENTREE : */
/*  ------------------- */
/*  NNO         : NOMBRE DE POINTS DE L ELEMENT */
/*  npt         : nb de points d'interpolation (tying) */
/*                npt = npt1 (rr,rz) +npt2 (ss,sz) + npt3 (rs) */
/*  xnorm       : normale aux noeuds */
/*  v1,v2       : forment un repere local avec norm (calc reploc) */
/*  npi         : nombre de points d'integration (en r,s) */
/*  xint,yint   : coord des pt d'integration en r,s */
/*  inmitc      : nom generique de la fonction d'interpolation */
/*                inmitc(j,i,x,y) */
/*  vdtq2(1,*)  : derivee par rapport a r aux pt interpolation */
/*  vdtq2(2,*)  : derivee par rapport a s aux pt interpolation */
/*  f           : matrice f=f+fz*z */

/* ............................................................ */
/*  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     PROGRAMMEUR  :  Marina Vidrascu INRIA 2000 */
/*  .................................................................. */
    /* Parameter adjustments */
    bzz -= 6;
    bz -= 6;
    b -= 6;
    --t;
    fz -= 13;
    f -= 13;
    vp2_dim1 = *nno;
    vp2_offset = vp2_dim1 + 1;
    vp2 -= vp2_offset;
    vdtq2_dim2 = *nno;
    vdtq2_offset = (vdtq2_dim2 + 1 << 1) + 1;
    vdtq2 -= vdtq2_offset;
    xnorm -= 4;
    v1 -= 4;
    v2 -= 4;
    vp2i_dim1 = *nno;
    vp2i_offset = vp2i_dim1 + 1;
    vp2i -= vp2i_offset;
    vdpq2_dim2 = *nno;
    vdpq2_offset = (vdpq2_dim2 + 1 << 1) + 1;
    vdpq2 -= vdpq2_offset;
    --yint;
    --xint;
    --irint;
    fi -= 4;
    fiz -= 4;

    /* Function Body */
    i__1 = *nno * 25;
    dcopy_(&i__1, &c_b2, &c__0, &b[6], &c__1);
    i__1 = *nno * 25;
    dcopy_(&i__1, &c_b2, &c__0, &bz[6], &c__1);
    i__1 = *nno * 25;
    dcopy_(&i__1, &c_b2, &c__0, &bzz[6], &c__1);
/*      ------------------------------------------------------ */

/*            COMPOSANTE 1 rr */

/*      ------------------------------------------------------ */
    if (irint[1] == 1) {
/*      ---- > reinterpolation */
	i__1 = *npt1;
	for (lt = 1; lt <= i__1; ++lt) {
	    i__2 = *nno;
	    for (n = 1; n <= i__2; ++n) {
		valdx = vdtq2[(n + lt * vdtq2_dim2 << 1) + 1] * (*inmitc)(&c__1, &lt, &xint[*l], &yint[*l]);
		valx = vp2[n + lt * vp2_dim1] * (*inmitc)(&c__1, &lt, &xint[*l], &yint[*l]);

/*          calcul de B_rr (B(1,...) */

/*          (B_rr),1 */
		b[n * 5 + 1] += f[(lt * 3 + 1) * 3 + 1] * valdx;
		bz[n * 5 + 1] += fz[(lt * 3 + 1) * 3 + 1] * valdx;
/*          (B_rr),2 */
		b[(n + *nno) * 5 + 1] += f[(lt * 3 + 2) * 3 + 1] * valdx;
		bz[(n + *nno) * 5 + 1] += fz[(lt * 3 + 2) * 3 + 1] * valdx;
/*          (B_rr),3 */
		b[(n + (*nno << 1)) * 5 + 1] += f[(lt * 3 + 3) * 3 + 1] * valdx;
		bz[(n + (*nno << 1)) * 5 + 1] += fz[(lt * 3 + 3) * 3 + 1] * valdx;
/*          (B_rr),4 */
		bz[(n + *nno * 3) * 5 + 1] += (f[(lt * 3 + 1) * 3 + 1] * v1[n * 3 + 1] + f[(lt * 3 + 2) * 3 + 1] * v1[n * 3 + 2] + f[(lt * 3 + 3) * 3 + 1] * v1[n * 3 + 3]) * valdx * (float).5 * t[n];
		bzz[(n + *nno * 3) * 5 + 1] += (fz[(lt * 3 + 1) * 3 + 1] * v1[n * 3 + 1] + fz[(lt * 3 + 2) * 3 + 1] * v1[n * 3 + 2] + fz[(lt * 3 + 3) * 3 + 1] * v1[n * 3 + 3]) * valdx * (float).5 * t[n];
/*          (B_rr),5 */
		bz[(n + (*nno << 2)) * 5 + 1] += (f[(lt * 3 + 1) * 3 + 1] * v2[n * 3 + 1] + f[(lt * 3 + 2) * 3 + 1] * v2[n * 3 + 2] + f[(lt * 3 + 3) * 3 + 1] * v2[n * 3 + 3]) * valdx * (float).5 * t[n];
		bzz[(n + (*nno << 2)) * 5 + 1] += (fz[(lt * 3 + 1) * 3 + 1] * v2[n * 3 + 1] + fz[(lt * 3 + 2) * 3 + 1] * v2[n * 3 + 2] + fz[(lt * 3 + 3) * 3 + 1] * v2[n * 3 + 3]) * valdx * (float).5 * t[n];
/* L1: */
	    }
	}
    } else {
/*      ---- > pas de reinterpolation */
	i__2 = *nno;
	for (n = 1; n <= i__2; ++n) {

/*          calcul de B_rr (B(1,...) */

/*          (B_rr),1 */
	    b[n * 5 + 1] = fi[4] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 1];
	    bz[n * 5 + 1] = fiz[4] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 1];
/*          (B_rr),2 */
	    b[(n + *nno) * 5 + 1] = fi[7] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 1];
	    bz[(n + *nno) * 5 + 1] = fiz[7] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 1];
/*          (B_rr),3 */
	    b[(n + (*nno << 1)) * 5 + 1] = fi[10] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 1];
	    bz[(n + (*nno << 1)) * 5 + 1] = fiz[10] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 1];
/*          (B_rr),4 */
	    bz[(n + *nno * 3) * 5 + 1] = (fi[4] * v1[n * 3 + 1] + fi[7] * v1[n * 3 + 2] + fi[10] * v1[n * 3 + 3]) * vdpq2[(n + *l * vdpq2_dim2 << 1) + 1] * (float).5 * t[n];
	    bzz[(n + *nno * 3) * 5 + 1] = (fiz[4] * v1[n * 3 + 1] + fiz[7] * v1[n * 3 + 2] + fiz[10] * v1[n * 3 + 3]) * vdpq2[(n + *l * vdpq2_dim2 << 1) + 1] * (float).5 * t[n];
/*          (B_rr),5 */
	    bz[(n + (*nno << 2)) * 5 + 1] = (fi[4] * v2[n * 3 + 1] + fi[7] * v2[n * 3 + 2] + fi[10] * v2[n * 3 + 3]) * vdpq2[(n + *l * vdpq2_dim2 << 1) + 1] * (float).5 * t[n];
	    bzz[(n + (*nno << 2)) * 5 + 1] = (fiz[4] * v2[n * 3 + 1] + fiz[7] * v2[n * 3 + 2] + fiz[10] * v2[n * 3 + 3]) * vdpq2[(n + *l * vdpq2_dim2 << 1) + 1] * (float).5 * t[n];
/* L10: */
	}
    }
/*      ------------------------------------------------------ */

/*            COMPOSANTE 4 rz */

/*      ------------------------------------------------------ */
    if (irint[4] == 1) {
/*      ---- > reinterpolation */
	i__2 = *npt1;
	for (lt = 1; lt <= i__2; ++lt) {
	    i__1 = *nno;
	    for (n = 1; n <= i__1; ++n) {
		valdx = vdtq2[(n + lt * vdtq2_dim2 << 1) + 1] * (*inmitc)(&c__1, &lt, &xint[*l], &yint[*l]);
		valx = vp2[n + lt * vp2_dim1] * (*inmitc)(&c__1, &lt, &xint[*l], &yint[*l]);

/*          calcul de B_rz (B(4,...) */

/*                  attention comme Eps_rz(r,s,z) =  Eps_rz(r,s,0) */
/*                                  Bz =Bzz = 0 */
/*          (B_rz),1 */
		b[n * 5 + 4] += f[(lt * 3 + 1) * 3 + 3] * valdx;
/*          (B_rz),2 */
		b[(n + *nno) * 5 + 4] += f[(lt * 3 + 2) * 3 + 3] * valdx;
/*          (B_rz),3 */
		b[(n + (*nno << 1)) * 5 + 4] += f[(lt * 3 + 3) * 3 + 3] * valdx;
/*          (B_rz),4 */
		b[(n + *nno * 3) * 5 + 4] += (f[(lt * 3 + 1) * 3 + 1] * v1[n * 3 + 1] + f[(lt * 3 + 2) * 3 + 1] * v1[n * 3 + 2] + f[(lt * 3 + 3) * 3 + 1] * v1[n * 3 + 3]) * valx * (float).5 * t[n];
/*          (B_rz),5 */
		b[(n + (*nno << 2)) * 5 + 4] += (f[(lt * 3 + 1) * 3 + 1] * v2[n * 3 + 1] + f[(lt * 3 + 2) * 3 + 1] * v2[n * 3 + 2] + f[(lt * 3 + 3) * 3 + 1] * v2[n * 3 + 3]) * valx * (float).5 * t[n];
/* L4: */
	    }
	}
    } else {
/*      ---- > pas de reinterpolation */
	i__1 = *nno;
	for (n = 1; n <= i__1; ++n) {

/*          calcul de B_rz (B(4,...) */

/*                  attention comme Eps_rz(r,s,z) =  Eps_rz(r,s,0) */
/*                                  Bz =Bzz = 0 */
/*          (B_rz),1 */
	    b[n * 5 + 4] = fi[6] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 1];
/*          (B_rz),2 */
	    b[(n + *nno) * 5 + 4] = fi[9] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 1];
/*          (B_rz),3 */
	    b[(n + (*nno << 1)) * 5 + 4] = fi[12] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 1];
/*          (B_rz),4 */
	    b[(n + *nno * 3) * 5 + 4] = (fi[4] * v1[n * 3 + 1] + fi[7] * v1[n * 3 + 2] + fi[10] * v1[n * 3 + 3]) * vp2i[n + *l * vp2i_dim1] * (float).5 * t[n];
/*          (B_rz),5 */
	    b[(n + (*nno << 2)) * 5 + 4] = (fi[4] * v2[n * 3 + 1] + fi[7] * v2[n * 3 + 2] + fi[10] * v2[n * 3 + 3]) * vp2i[n + *l * vp2i_dim1] * (float).5 * t[n];
/* L40: */
	}
    }
/*      ------------------------------------------------------ */

/*            COMPOSANTE 2 ss */

/*      ------------------------------------------------------ */
    if (irint[2] == 1) {
/*      ---- > reinterpolation */
	i__1 = *npt1 + *npt2;
	for (lt = *npt1 + 1; lt <= i__1; ++lt) {
	    i__2 = *nno;
	    for (n = 1; n <= i__2; ++n) {
		i__3 = lt - *npt1;
		valdy = vdtq2[(n + lt * vdtq2_dim2 << 1) + 2] * (*inmitc)(&c__2, &i__3, &xint[*l], &yint[*l]);
		i__3 = lt - *npt1;
		valy = vp2[n + lt * vp2_dim1] * (*inmitc)(&c__2, &i__3, &xint[*l], &yint[*l]);
/*          calcul de B_ss (B(2,...) */

/*          (B_ss),1 */
		b[n * 5 + 2] += f[(lt * 3 + 1) * 3 + 2] * valdy;
		bz[n * 5 + 2] += fz[(lt * 3 + 1) * 3 + 2] * valdy;
/*          (B_ss),2 */
		b[(n + *nno) * 5 + 2] += f[(lt * 3 + 2) * 3 + 2] * valdy;
		bz[(n + *nno) * 5 + 2] += fz[(lt * 3 + 2) * 3 + 2] * valdy;
/*          (B_ss),3 */
		b[(n + (*nno << 1)) * 5 + 2] += f[(lt * 3 + 3) * 3 + 2] * valdy;
		bz[(n + (*nno << 1)) * 5 + 2] += fz[(lt * 3 + 3) * 3 + 2] * valdy;
/*          (B_ss),4 */
		bz[(n + *nno * 3) * 5 + 2] += (f[(lt * 3 + 1) * 3 + 2] * v1[n * 3 + 1] + f[(lt * 3 + 2) * 3 + 2] * v1[n * 3 + 2] + f[(lt * 3 + 3) * 3 + 2] * v1[n * 3 + 3]) * valdy * (float).5 * t[n];
		bzz[(n + *nno * 3) * 5 + 2] += (fz[(lt * 3 + 1) * 3 + 2] * v1[n * 3 + 1] + fz[(lt * 3 + 2) * 3 + 2] * v1[n * 3 + 2] + fz[(lt * 3 + 3) * 3 + 2] * v1[n * 3 + 3]) * valdy * (float).5 * t[n];
/*          (B_ss),5 */
		bz[(n + (*nno << 2)) * 5 + 2] += (f[(lt * 3 + 1) * 3 + 2] * v2[n * 3 + 1] + f[(lt * 3 + 2) * 3 + 2] * v2[n * 3 + 2] + f[(lt * 3 + 3) * 3 + 2] * v2[n * 3 + 3]) * valdy * (float).5 * t[n];
		bzz[(n + (*nno << 2)) * 5 + 2] += (fz[(lt * 3 + 1) * 3 + 2] * v2[n * 3 + 1] + fz[(lt * 3 + 2) * 3 + 2] * v2[n * 3 + 2] + fz[(lt * 3 + 3) * 3 + 2] * v2[n * 3 + 3]) * valdy * (float).5 * t[n];
/* L2: */
	    }
	}
    } else {
/*      ---- > pas de reinterpolation */
	i__2 = *nno;
	for (n = 1; n <= i__2; ++n) {
/*          calcul de B_ss (B(2,...) */

/*          (B_ss),1 */
	    b[n * 5 + 2] = fi[5] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 2];
	    bz[n * 5 + 2] = fiz[5] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 2];
/*          (B_ss),2 */
	    b[(n + *nno) * 5 + 2] = fi[8] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 2];
	    bz[(n + *nno) * 5 + 2] = fiz[8] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 2];
/*          (B_ss),3 */
	    b[(n + (*nno << 1)) * 5 + 2] = fi[11] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 2];
	    bz[(n + (*nno << 1)) * 5 + 2] = fiz[11] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 2];
/*          (B_ss),4 */
	    bz[(n + *nno * 3) * 5 + 2] = (fi[5] * v1[n * 3 + 1] + fi[8] * v1[n * 3 + 2] + fi[11] * v1[n * 3 + 3]) * vdpq2[(n + *l * vdpq2_dim2 << 1) + 2] * (float).5 * t[n];
	    bzz[(n + *nno * 3) * 5 + 2] = (fiz[5] * v1[n * 3 + 1] + fiz[8] * v1[n * 3 + 2] + fiz[11] * v1[n * 3 + 3]) * vdpq2[(n + *l * vdpq2_dim2 << 1) + 2] * (float).5 * t[n];
/*          (B_ss),5 */
	    bz[(n + (*nno << 2)) * 5 + 2] = (fi[5] * v2[n * 3 + 1] + fi[8] * v2[n * 3 + 2] + fi[11] * v2[n * 3 + 3]) * vdpq2[(n + *l * vdpq2_dim2 << 1) + 2] * (float).5 * t[n];
	    bzz[(n + (*nno << 2)) * 5 + 2] = (fiz[5] * v2[n * 3 + 1] + fiz[8] * v2[n * 3 + 2] + fiz[11] * v2[n * 3 + 3]) * vdpq2[(n + *l * vdpq2_dim2 << 1) + 2] * (float).5 * t[n];
/* L20: */
	}
    }
/*      ------------------------------------------------------ */

/*            COMPOSANTE 5 sz */

/*      ------------------------------------------------------ */
    if (irint[5] == 1) {
/*      ---- > reinterpolation */
	i__2 = *npt1 + *npt2;
	for (lt = *npt1 + 1; lt <= i__2; ++lt) {
	    i__1 = *nno;
	    for (n = 1; n <= i__1; ++n) {
		i__3 = lt - *npt1;
		valdy = vdtq2[(n + lt * vdtq2_dim2 << 1) + 2] * (*inmitc)(&c__2, &i__3, &xint[*l], &yint[*l]);
		i__3 = lt - *npt1;
		valy = vp2[n + lt * vp2_dim1] * (*inmitc)(&c__2, &i__3, &xint[*l], &yint[*l]);

/*          calcul de B_sz (B(5,...) */


/*                  attention comme Eps_sz(r,s,z) =  Eps_sz(r,s,0) */
/*                                  Bz =Bzz = 0 */
/*          (B_sz),1 */
		b[n * 5 + 5] += f[(lt * 3 + 1) * 3 + 3] * valdy;
/*          (B_sz),2 */
		b[(n + *nno) * 5 + 5] += f[(lt * 3 + 2) * 3 + 3] * valdy;
/*          (B_sz),3 */
		b[(n + (*nno << 1)) * 5 + 5] += f[(lt * 3 + 3) * 3 + 3] * valdy;
/*          (B_sz),4 */
		b[(n + *nno * 3) * 5 + 5] += (f[(lt * 3 + 1) * 3 + 2] * v1[n * 3 + 1] + f[(lt * 3 + 2) * 3 + 2] * v1[n * 3 + 2] + f[(lt * 3 + 3) * 3 + 2] * v1[n * 3 + 3]) * valy * (float).5 * t[n];
/*          (B_sz),5 */
		b[(n + (*nno << 2)) * 5 + 5] += (f[(lt * 3 + 1) * 3 + 2] * v2[n * 3 + 1] + f[(lt * 3 + 2) * 3 + 2] * v2[n * 3 + 2] + f[(lt * 3 + 3) * 3 + 2] * v2[n * 3 + 3]) * valy * (float).5 * t[n];
/* L5: */
	    }
	}
    } else {
/*      ---- > pas de reinterpolation */
	i__1 = *nno;
	for (n = 1; n <= i__1; ++n) {

/*          calcul de B_sz (B(5,...) */


/*                  attention comme Eps_sz(r,s,z) =  Eps_sz(r,s,0) */
/*                                  Bz =Bzz = 0 */
/*          (B_sz),1 */
	    b[n * 5 + 5] = fi[6] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 2];
/*          (B_sz),2 */
	    b[(n + *nno) * 5 + 5] = fi[9] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 2];
/*          (B_sz),3 */
	    b[(n + (*nno << 1)) * 5 + 5] = fi[12] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 2];
/*          (B_sz),4 */
	    b[(n + *nno * 3) * 5 + 5] = (fi[5] * v1[n * 3 + 1] + fi[8] * v1[n * 3 + 2] + fi[11] * v1[n * 3 + 3]) * vp2i[n + *l * vp2i_dim1] * (float).5 * t[n];
/*          (B_sz),5 */
	    b[(n + (*nno << 2)) * 5 + 5] = (fi[5] * v2[n * 3 + 1] + fi[8] * v2[n * 3 + 2] + fi[11] * v2[n * 3 + 3]) * vp2i[n + *l * vp2i_dim1] * (float).5 * t[n];
/* L50: */
	}
    }
/*      ------------------------------------------------------ */

/*            COMPOSANTE 3 rs */

/*      ------------------------------------------------------ */
    if (irint[3] == 1) {
/*      ---- > reinterpolation */
	i__1 = *npt;
	for (lt = *npt1 + *npt2 + 1; lt <= i__1; ++lt) {
	    i__2 = *nno;
	    for (n = 1; n <= i__2; ++n) {
		i__3 = lt - *npt1 - *npt2;
		valdx = vdtq2[(n + lt * vdtq2_dim2 << 1) + 1] * (*inmitc)(&c__3, &i__3, &xint[*l], &yint[*l]);
		i__3 = lt - *npt1 - *npt2;
		valdy = vdtq2[(n + lt * vdtq2_dim2 << 1) + 2] * (*inmitc)(&c__3, &i__3, &xint[*l], &yint[*l]);

/*          calcul de B_rs (B(3,...) */

/*          (B_rs),1 */
		b[n * 5 + 3] = b[n * 5 + 3] + f[(lt * 3 + 1) * 3 + 2] * valdx + f[(lt * 3 + 1) * 3 + 1] * valdy;
		bz[n * 5 + 3] = bz[n * 5 + 3] + fz[(lt * 3 + 1) * 3 + 2] * valdx + fz[(lt * 3 + 1) * 3 + 1] * valdy;
/*          (B_rs),2 */
		b[(n + *nno) * 5 + 3] = b[(n + *nno) * 5 + 3] + f[(lt * 3 + 2) * 3 + 2] * valdx + f[(lt * 3 + 2) * 3 + 1] * valdy;
		bz[(n + *nno) * 5 + 3] = bz[(n + *nno) * 5 + 3] + fz[(lt * 3 + 2) * 3 + 2] * valdx + fz[(lt * 3 + 2) * 3 + 1] * valdy;
/*          (B_rs),3 */
		b[(n + (*nno << 1)) * 5 + 3] = b[(n + (*nno << 1)) * 5 + 3] + f[(lt * 3 + 3) * 3 + 2] * valdx + f[(lt * 3 + 3) * 3 + 1] * valdy;
		bz[(n + (*nno << 1)) * 5 + 3] = bz[(n + (*nno << 1)) * 5 + 3] + fz[(lt * 3 + 3) * 3 + 2] * valdx + fz[(lt * 3 + 3) * 3 + 1] * valdy;
/*          (B_rs),4 */
		bz[(n + *nno * 3) * 5 + 3] = bz[(n + *nno * 3) * 5 + 3] + (f[(lt * 3 + 1) * 3 + 2] * v1[n * 3 + 1] + f[(lt * 3 + 2) * 3 + 2] * v1[n * 3 + 2] + f[(lt * 3 + 3) * 3 + 2] * v1[n * 3 + 3]) * valdx * (float).5 * t[n] + (f[(lt * 3 + 1) * 3 + 1] * v1[n * 3 + 1] + f[(lt * 3 + 2) * 3 + 1] * v1[n * 3 + 2] + f[(lt * 3 + 3) * 3 + 1] * v1[n * 3 + 3]) * valdy * (float).5 * t[n];
		bzz[(n + *nno * 3) * 5 + 3] = bzz[(n + *nno * 3) * 5 + 3] + (fz[(lt * 3 + 1) * 3 + 2] * v1[n * 3 + 1] + fz[(lt * 3 + 2) * 3 + 2] * v1[n * 3 + 2] + fz[(lt * 3 + 3) * 3 + 2] * v1[n * 3 + 3]) * valdx * (float).5 * t[n] + (fz[(lt * 3 + 1) * 3 + 1] * v1[n * 3 + 1] + fz[(lt * 3 + 2) * 3 + 1] * v1[n * 3 + 2] + fz[(lt * 3 + 3) * 3 + 1] * v1[n * 3 + 3]) * valdy * (float).5 * t[n];
/*          (B_rs),5 */
		bz[(n + (*nno << 2)) * 5 + 3] = bz[(n + (*nno << 2)) * 5 + 3] + (f[(lt * 3 + 1) * 3 + 2] * v2[n * 3 + 1] + f[(lt * 3 + 2) * 3 + 2] * v2[n * 3 + 2] + f[(lt * 3 + 3) * 3 + 2] * v2[n * 3 + 3]) * valdx * (float).5 * t[n] + (f[(lt * 3 + 1) * 3 + 1] * v2[n * 3 + 1] + f[(lt * 3 + 2) * 3 + 1] * v2[n * 3 + 2] + f[(lt * 3 + 3) * 3 + 1] * v2[n * 3 + 3]) * valdy * (float).5 * t[n];
		bzz[(n + (*nno << 2)) * 5 + 3] = bzz[(n + (*nno << 2)) * 5 + 3] + (fz[(lt * 3 + 1) * 3 + 2] * v2[n * 3 + 1] + fz[(lt * 3 + 2) * 3 + 2] * v2[n * 3 + 2] + fz[(lt * 3 + 3) * 3 + 2] * v2[n * 3 + 3]) * valdx * (float).5 * t[n] + (fz[(lt * 3 + 1) * 3 + 1] * v2[n * 3 + 1] + fz[(lt * 3 + 2) * 3 + 1] * v2[n * 3 + 2] + fz[(lt * 3 + 3) * 3 + 1] * v2[n * 3 + 3]) * valdy * (float).5 * t[n];
/* L3: */
	    }
	}
    } else {
/*      ---- > pas de reinterpolation */
	i__2 = *nno;
	for (n = 1; n <= i__2; ++n) {

/*          calcul de B_rs (B(3,...) */

/*          (B_rs),1 */
	    b[n * 5 + 3] = fi[5] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 1] + fi[4] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 2];
	    bz[n * 5 + 3] = fiz[5] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 1] + fiz[4] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 2];
/*          (B_rs),2 */
	    b[(n + *nno) * 5 + 3] = fi[8] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 1] + fi[7] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 2];
	    bz[(n + *nno) * 5 + 3] = fiz[8] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 1] + fiz[7] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 2];
/*          (B_rs),3 */
	    b[(n + (*nno << 1)) * 5 + 3] = fi[11] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 1] + fi[10] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 2];
	    bz[(n + (*nno << 1)) * 5 + 3] = fiz[11] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 1] + fiz[10] * vdpq2[(n + *l * vdpq2_dim2 << 1) + 2];
/*          (B_rs),4 */
	    bz[(n + *nno * 3) * 5 + 3] = (fi[5] * v1[n * 3 + 1] + fi[8] * v1[n * 3 + 2] + fi[11] * v1[n * 3 + 3]) * vdpq2[(n + *l * vdpq2_dim2 << 1) + 1] * (float).5 * t[n] + (fi[4] * v1[n * 3 + 1] + fi[7] * v1[n * 3 + 2] + fi[10] * v1[n * 3 + 3]) * vdpq2[(n + *l * vdpq2_dim2 << 1) + 2] * (float).5 * t[n];
	    bzz[(n + *nno * 3) * 5 + 3] = (fiz[5] * v1[n * 3 + 1] + fiz[8] * v1[n * 3 + 2] + fiz[11] * v1[n * 3 + 3]) * vdpq2[(n + *l * vdpq2_dim2 << 1) + 1] * (float).5 * t[n] + (fiz[4] * v1[n * 3 + 1] + fiz[7] * v1[n * 3 + 2] + fiz[10] * v1[n * 3 + 3]) * vdpq2[(n + *l * vdpq2_dim2 << 1) + 2] * (float).5 * t[n];
/*          (B_rs),5 */
	    bz[(n + (*nno << 2)) * 5 + 3] = (fi[5] * v2[n * 3 + 1] + fi[8] * v2[n * 3 + 2] + fi[11] * v2[n * 3 + 3]) * vdpq2[(n + *l * vdpq2_dim2 << 1) + 1] * (float).5 * t[n] + (fi[4] * v2[n * 3 + 1] + fi[7] * v2[n * 3 + 2] + fi[10] * v2[n * 3 + 3]) * vdpq2[(n + *l * vdpq2_dim2 << 1) + 2] * (float).5 * t[n];
	    bzz[(n + (*nno << 2)) * 5 + 3] = (fiz[5] * v2[n * 3 + 1] + fiz[8] * v2[n * 3 + 2] + fiz[11] * v2[n * 3 + 3]) * vdpq2[(n + *l * vdpq2_dim2 << 1) + 1] * (float).5 * t[n] + (fiz[4] * v2[n * 3 + 1] + fiz[7] * v2[n * 3 + 2] + fiz[10] * v2[n * 3 + 3]) * vdpq2[(n + *l * vdpq2_dim2 << 1) + 2] * (float).5 * t[n];
/* L30: */
	}
    }
} /* bremit_ */

