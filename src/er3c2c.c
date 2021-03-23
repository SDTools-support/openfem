/* er3c2c.f -- translated by f2c (version 19961017).
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

/* Subroutine */ int er3c2c_(nno, npo, x, y, z__, npi, ijt, poids, vp1, vdpq2, vdpq1, iopt, e, xnu, elas, ae, delta, xint, yint, zint, a2)
int32 *nno, *npo;
doublereal *x, *y, *z__;
int32 *npi, *ijt;
doublereal *poids, *vp1, *vdpq2, *vdpq1;
int32 *iopt;
doublereal *e, *xnu, *elas, *ae, *delta, *xint, *yint, *zint, *a2;
{
    /* System generated locals */
    int32 vp1_dim1, vp1_offset, vdpq1_dim2, vdpq1_offset, vdpq2_dim2, vdpq2_offset, a2_dim2, a2_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static doublereal xlam;
    static int32 i__, j, l;
    static doublereal dfinv[9]	/* was [3][3] */;
    extern /* Subroutine */ int dcopy_();
    static int32 k1, k2, k3;
    static doublereal df[9]	/* was [3][3] */, fi[1];
    static int32 ii, jj, nm, indice;
    extern /* Subroutine */ int fobase_();
    static doublereal xmu;

/*  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  BUT :   CALCUL DE LA MATRICE ELEMENTAIRE DE RAIDEUR */
/*  ---     D UN ELEMENT  HEXA3Q2D OU HEXA3Q2C ISOTROPE OU ORTHOTROPE */
/*  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  PARAMETRES D'ENTREE : */
/*  ------------------- */
/*  NNO         : NOMBRE DE NOEUDS DE L ELEMENT */
/*  NPO         : NOMBRE DE POINTS DE L ELEMENT */
/*  X,Y,Z       : COORDONNEES DES NOEUDS DE L ELEMENT */
/*  NPI         : NOMBRE DE POINTS  D INTEGRATION  SUR LE VOLUME */
/*  IJT         : PERMUTATION */
/*  POIDS       : POIDS DE LA FORMULE D INTEGRATION SUR LE CUBE */
/*  VP1         : VALEURS DES  DES POL DE BASE */
/*  NBAR        : NUMERO DU BARYCENTRE */
/*  VDPQ2       : VALEURS DES DERIVES DES POL DE BASE DE Q'2(R3) */
/*                AUX POINTS D INTEGRATION */
/*  VDPQ1       : VALEURS DES DERIVES DES POLYNOMES DE BASE DE Q1(R3) */
/*                AUX POINTS D INTEGRATION */
/*  ELAS        : MATRICE DE L ELASTICITE */
/*  DELT */
/*  DF */
/*  DFINV */
/*  XINT,YINT,ZINT */
/*  A2 */
/*  PARAMETRES RESULTATS : */
/*  -------------------- */
/*  AE          : MATRICE ELEMENTAIRE DE RAIDEUR STOCKEE SYMETRIQUE */
/*  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     PROGRAMMEUR  : MARINA VIDRASCU INRIA 2001 */
/*  .................................................................. */

/*     ---- non utilises --- */


    /* Parameter adjustments */
    --z__;
    --y;
    --x;
    a2_dim2 = *nno;
    a2_offset = (a2_dim2 + 1) * 3 + 1;
    a2 -= a2_offset;
    --zint;
    --yint;
    --xint;
    --delta;
    vdpq1_dim2 = *npo;
    vdpq1_offset = (vdpq1_dim2 + 1) * 3 + 1;
    vdpq1 -= vdpq1_offset;
    vdpq2_dim2 = *nno;
    vdpq2_offset = (vdpq2_dim2 + 1) * 3 + 1;
    vdpq2 -= vdpq2_offset;
    vp1_dim1 = *nno;
    vp1_offset = vp1_dim1 + 1;
    vp1 -= vp1_offset;
    --poids;
    --ijt;
    --elas;
    --ae;

    /* Function Body */
    nm = *nno * 3 * (*nno * 3 + 1) / 2;
    dcopy_(&nm, &c_b2, &c__0, &ae[1], &c__1);
    indice = 1;
    fobase_(&c__3, &c__3, nno, npo, npi, &vp1[vp1_offset], &vp1[vp1_offset], &vdpq2[vdpq2_offset], &vdpq1[vdpq1_offset], &x[1], &y[1], &z__[1], &a2[a2_offset], &xint[1], &yint[1], &zint[1], &delta[1], dfinv, df, &indice, fi, fi, fi, iopt);
/*           -------------- */
    if (*iopt == 2) {
	xmu = *e / ((*xnu + (float)1.) * (float)2.);
	xlam = *e * *xnu / ((*xnu + (float)1.) * ((float)1. - *xnu * (float)2.));
	elas[1] = xlam + xmu * (float)2.;
	elas[3] = xlam + xmu * (float)2.;
	elas[6] = xlam + xmu * (float)2.;
	elas[2] = xlam;
	elas[4] = xlam;
	elas[5] = xlam;
	elas[7] = xmu;
	elas[8] = xmu;
	elas[9] = xmu;
    }

/*       pour traiter l'orthotrope donner directement elas! */
/*       LES BLOCS DIAGONAUX */
/*       ------------------- */
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
		ae[k1] += poids[l] / delta[l] * (elas[1] * a2[(i__ + l * a2_dim2) * 3 + 1] * a2[(j + l * a2_dim2) * 3 + 1] + elas[7] * a2[(i__ + l * a2_dim2) * 3 + 2] * a2[(j + l * a2_dim2) * 3 + 2] + elas[9] * a2[(i__ + l * a2_dim2) * 3 + 3] * a2[(j + l * a2_dim2) * 3 + 3]);
		ae[k2] += poids[l] / delta[l] * (elas[7] * a2[(i__ + l * a2_dim2) * 3 + 1] * a2[(j + l * a2_dim2) * 3 + 1] + elas[3] * a2[(i__ + l * a2_dim2) * 3 + 2] * a2[(j + l * a2_dim2) * 3 + 2] + elas[8] * a2[(i__ + l * a2_dim2) * 3 + 3] * a2[(j + l * a2_dim2) * 3 + 3]);
		ae[k3] += poids[l] / delta[l] * (elas[9] * a2[(i__ + l * a2_dim2) * 3 + 1] * a2[(j + l * a2_dim2) * 3 + 1] + elas[8] * a2[(i__ + l * a2_dim2) * 3 + 2] * a2[(j + l * a2_dim2) * 3 + 2] + elas[6] * a2[(i__ + l * a2_dim2) * 3 + 3] * a2[(j + l * a2_dim2) * 3 + 3]);
/* L1: */
	    }
/* L9: */
	}
    }

/*       LES BLOCS EXTRA-DIAGONAUX */
/*       ------------------------- */
    i__2 = *nno;
    for (j = 1; j <= i__2; ++j) {
	i__1 = *nno;
	for (i__ = 1; i__ <= i__1; ++i__) {

/*          BLOC 12 */
/*          ------- */
/* Computing MIN */
	    i__3 = ijt[i__], i__4 = ijt[*nno + j];
	    ii = min(i__3,i__4);
/* Computing MAX */
	    i__3 = ijt[i__], i__4 = ijt[*nno + j];
	    jj = max(i__3,i__4);
	    k1 = jj * (jj - 1) / 2 + ii;

/*          BLOC 13 */
/*          ------- */
/* Computing MIN */
	    i__3 = ijt[i__], i__4 = ijt[(*nno << 1) + j];
	    ii = min(i__3,i__4);
/* Computing MAX */
	    i__3 = ijt[i__], i__4 = ijt[(*nno << 1) + j];
	    jj = max(i__3,i__4);
	    k2 = jj * (jj - 1) / 2 + ii;

/*          BLOC 23 */
/*          ------- */
/* Computing MIN */
	    i__3 = ijt[*nno + i__], i__4 = ijt[(*nno << 1) + j];
	    ii = min(i__3,i__4);
/* Computing MAX */
	    i__3 = ijt[*nno + i__], i__4 = ijt[(*nno << 1) + j];
	    jj = max(i__3,i__4);
	    k3 = jj * (jj - 1) / 2 + ii;
	    i__3 = *npi;
	    for (l = 1; l <= i__3; ++l) {
		ae[k1] += poids[l] / delta[l] * (elas[7] * a2[(i__ + l * a2_dim2) * 3 + 2] * a2[(j + l * a2_dim2) * 3 + 1] + elas[2] * a2[(i__ + l * a2_dim2) * 3 + 1] * a2[(j + l * a2_dim2) * 3 + 2]);
		ae[k2] += poids[l] / delta[l] * (elas[9] * a2[(i__ + l * a2_dim2) * 3 + 3] * a2[(j + l * a2_dim2) * 3 + 1] + elas[4] * a2[(i__ + l * a2_dim2) * 3 + 1] * a2[(j + l * a2_dim2) * 3 + 3]);
		ae[k3] += poids[l] / delta[l] * (elas[8] * a2[(i__ + l * a2_dim2) * 3 + 3] * a2[(j + l * a2_dim2) * 3 + 2] + elas[5] * a2[(i__ + l * a2_dim2) * 3 + 2] * a2[(j + l * a2_dim2) * 3 + 3]);
/* L4: */
	    }
/* L10: */
	}
    }
} /* er3c2c_ */

