/* mitlin.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__3 = 3;

/* Subroutine */ int mitlin_(nno, npt1, npt2, npt3, npt, irint, npi, poids, npiz, zint, poiz, xint, yint, vdtq2, vt2, f, fz, inmitc, e, enu, t, x, y, z__, vdpq2, vp2, b, bz, bzz, xnorm, v1, v2, fi, fiz, gzz, tt, ijt, ae, nmat)
int32 *nno, *npt1, *npt2, *npt3, *npt, *irint, *npi;
doublereal *poids;
int32 *npiz;
doublereal *zint, *poiz, *xint, *yint, *vdtq2, *vt2, *f, *fz;
/* Subroutine */ int (*inmitc) ();
doublereal *e, *enu, *t, *x, *y, *z__, *vdpq2, *vp2, *b, *bz, *bzz, *xnorm, *v1, *v2, *fi, *fiz, *gzz, *tt;
int32 *ijt;
doublereal *ae;
int32 *nmat;
{
    /* System generated locals */
    int32 vdpq2_dim2, vdpq2_offset, vdtq2_dim2, vdtq2_offset, vp2_dim1, vp2_offset, vt2_dim1, vt2_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal elas[9];
    extern doublereal ddot_();
    static doublereal girr, girs, giss, grrz, grsz, gssz;
    static int32 i__, j, k, l, m;
    static doublereal s;
    static int32 k1;
    static doublereal grrzz, grszz, gsszz;
    static int32 ie, ii, jj, lz;
    static doublereal deltat;
    extern /* Subroutine */ int bremit_(), melmit_();
    static doublereal grr, grs, gss;

/* ............................................................ */
/*  BUT :   CALCUL DE LA MATRICE ELEMENTAIRE DE RAIDEUR  pour */
/*  ---     une coque mitc */

/*   calcul de la matrice d'elasticite elas au points d'integration */
/*   (polynome de degre 2 en z) */
/*   calcul de tB*elas*B ou B= */
/* ............................................................ */

/*   nno           : nombre de noeuds de l'element */
/*   npi           : nombre de points d'integration (en r,s) */
/*   poids         : poids de la formule d'integration en r,s */
/*   npiz          : nombre de points d'integration en z */
/*   zint(npiz)    : coord des nzpi points d'integration en z sur [-1,1] */
/*   poiz          : poids de la formule d'integration en z */
/*   e,enu         : coeff de Young, Poisson  aux points d'integration */
/*   t             : epaisseur de la coque (donee aux noeuds) */
/*   x,y,z         : coordonnes des noeuds de l'element courant */
/*   vdpq2(1,*)    : derivee par rapport a r aux pt integration */
/*   vdpq2(2,*)    : derivee par rapport a s aux pt integration */
/*   vp2           : valeur des polynomes de base aux pt d'integration */
/*   b,bz,bzz      : B = b+ bz*z+ bzz*z**2 */
/*   xnorm         : normale aux noeuds */
/*   ijt           : permutation composante / noeud */
/*   ae            : matrice elementaire triangulaire sup */
/*  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     PROGRAMMEUR  :  Marina Vidrascu INRIA 2002 */
/*  .................................................................. */

/*     calcul de la matrice de rigidite */
/*     -------------------------------- */
/*     Somme(omega) Somme(-1,1) {BT * C * B} */


/*     boucle sur les points d'integration (r,s) */

    /* Parameter adjustments */
    --ijt;
    v2 -= 4;
    v1 -= 4;
    xnorm -= 4;
    bzz -= 6;
    bz -= 6;
    b -= 6;
    --z__;
    --y;
    --x;
    --t;
    fz -= 13;
    f -= 13;
    vt2_dim1 = *nno;
    vt2_offset = vt2_dim1 + 1;
    vt2 -= vt2_offset;
    vdtq2_dim2 = *nno;
    vdtq2_offset = (vdtq2_dim2 + 1 << 1) + 1;
    vdtq2 -= vdtq2_offset;
    --tt;
    --gzz;
    fiz -= 13;
    fi -= 13;
    vp2_dim1 = *nno;
    vp2_offset = vp2_dim1 + 1;
    vp2 -= vp2_offset;
    vdpq2_dim2 = *nno;
    vdpq2_offset = (vdpq2_dim2 + 1 << 1) + 1;
    vdpq2 -= vdpq2_offset;
    --yint;
    --xint;
    --poids;
    --poiz;
    --zint;
    --ae;

    /* Function Body */
    i__1 = *npi;
    for (l = 1; l <= i__1; ++l) {

/*     calcul de B au points d'integration (utiliser B au tying points) */

	bremit_(nno, npt1, npt2, npt3, npt, &xnorm[4], &v1[4], &v2[4], &t[1], npi, &xint[1], &yint[1], &vdtq2[vdtq2_offset], &vt2[vt2_offset], &f[13], &fz[13], inmitc, &b[6], &bz[6], &bzz[6], irint, &fi[(l * 3 + 1) * 3 + 1], &fiz[(l * 3 + 1) * 3 + 1], &vdpq2[vdpq2_offset], &vp2[vp2_offset], &l);

/*        calcul de grr gss grs comme polynomes de degre 2 en z */
/*        gr (lig 1 de fi) gs (lig 2 de fi) */
	grr = ddot_(&c__3, &fi[(l * 3 + 1) * 3 + 1], &c__3, &fi[(l * 3 + 1) * 3 + 1], &c__3);
	grrz = ddot_(&c__3, &fi[(l * 3 + 1) * 3 + 1], &c__3, &fiz[(l * 3 + 1) * 3 + 1], &c__3) * 2;
	grrzz = ddot_(&c__3, &fiz[(l * 3 + 1) * 3 + 1], &c__3, &fiz[(l * 3 + 1) * 3 + 1], &c__3);
	gss = ddot_(&c__3, &fi[(l * 3 + 1) * 3 + 2], &c__3, &fi[(l * 3 + 1) * 3 + 2], &c__3);
	gssz = ddot_(&c__3, &fi[(l * 3 + 1) * 3 + 2], &c__3, &fiz[(l * 3 + 1) * 3 + 2], &c__3) * 2;
	gsszz = ddot_(&c__3, &fiz[(l * 3 + 1) * 3 + 2], &c__3, &fiz[(l * 3 + 1) * 3 + 2], &c__3);
	grs = ddot_(&c__3, &fi[(l * 3 + 1) * 3 + 1], &c__3, &fi[(l * 3 + 1) * 3 + 2], &c__3);
	grsz = ddot_(&c__3, &fi[(l * 3 + 1) * 3 + 1], &c__3, &fiz[(l * 3 + 1) * 3 + 2], &c__3) + ddot_(&c__3, &fi[(l * 3 + 1) * 3 + 2], &c__3, &fiz[(l * 3 + 1) * 3 + 1], &c__3);
	grszz = ddot_(&c__3, &fiz[(l * 3 + 1) * 3 + 1], &c__3, &fiz[(l * 3 + 1) * 3 + 2], &c__3);

/*        boucle sur les points d'integration z */

	i__2 = *npiz;
	for (lz = 1; lz <= i__2; ++lz) {

/*           calcul de la matrice d'elasticite */

	    melmit_(&grr, &grrz, &grrzz, &grs, &grsz, &grszz, &gss, &gssz, &gsszz, &gzz[1], &zint[1], npiz, npi, enu, e, &tt[1], &l, &lz, &deltat, elas, &girr, &giss, &girs);

/*           calcul de TB*elas*B */

	    i__3 = *nno * 5;
	    for (j = 1; j <= i__3; ++j) {
		i__4 = j;
		for (i__ = 1; i__ <= i__4; ++i__) {

/*              permutation ds la matrice */

/* Computing MIN */
		    i__5 = ijt[i__], i__6 = ijt[j];
		    ii = min(i__5,i__6);
/* Computing MAX */
		    i__5 = ijt[i__], i__6 = ijt[j];
		    jj = max(i__5,i__6);
		    k1 = jj * (jj - 1) / 2 + ii;
		    s = (float)0.;
/*              2 blocs 1 3*3 et un 3*2 */
		    for (m = 1; m <= 3; ++m) {
			for (k = 1; k <= 3; ++k) {
			    if (k <= m) {
				ie = (m - 1) * m / 2 + k;
			    } else {
				ie = (k - 1) * k / 2 + m;
			    }
/*                 s = s+ somme (k=1,3) B(m,i)*elas(m,k)*B(k,j) */
/* Computing 2nd power */
			    d__1 = zint[lz];
/* Computing 2nd power */
			    d__2 = zint[lz];
			    s += elas[ie - 1] * (b[m + i__ * 5] + zint[lz] * bz[m + i__ * 5] + d__1 * d__1 * bzz[m + i__ * 5]) * (b[k + j * 5] + zint[lz] * bz[k + j * 5] + d__2 * d__2 * bzz[k + j * 5]);
/* L4: */
			}
		    }
		    for (m = 4; m <= 5; ++m) {
			for (k = 4; k <= 5; ++k) {
			    if (k <= m) {
				ie = (m - 4) * (m - 3) / 2 + 6 + k - 3;
			    } else {
				ie = (k - 4) * (k - 3) / 2 + 6 + m - 3;
			    }
/*                 s = s+ somme (k=1,3) B(m,i)*elas(m,k)*B(k,j) */
			    s += elas[ie - 1] * b[m + i__ * 5] * b[k + j * 5];
/* L5: */
			}
		    }
		    ae[k1] += s * deltat * poids[l] * poiz[lz];
/* L3: */
		}
	    }
/* L10: */
	}
/* L2: */
    }
} /* mitlin_ */

