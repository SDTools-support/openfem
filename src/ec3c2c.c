/* ec3c2c.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__3 = 3;

/* Subroutine */ int ec3c2c_(nno, npo, x, y, z__, ijt, vdpq2, vdpq1, iopt, e, xnu, elas, sigmae, sigma, u, delta, a2)
int32 *nno, *npo;
doublereal *x, *y, *z__;
int32 *ijt;
doublereal *vdpq2, *vdpq1;
int32 *iopt;
doublereal *e, *xnu, *elas, *sigmae, *sigma, *u, *delta, *a2;
{
    /* System generated locals */
    int32 vdpq1_dim2, vdpq1_offset, vdpq2_dim2, vdpq2_offset, a2_dim2, a2_offset, sigmae_dim2, sigmae_offset, i__1, i__2;

    /* Local variables */
    static doublereal xlam, xint[1], yint[1], zint[1];
    static int32 i__, j, l;
    static doublereal dfinv[9]	/* was [3][3] */;
    static int32 j1;
    static doublereal df[9]	/* was [3][3] */, fi[1];
    static int32 indice;
    extern /* Subroutine */ int fobase_();
    static doublereal vp1[1], xmu;

/*  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  BUT : CALCUL DES TABLEAUX SIGMAE ET SIGMA */
/*  ---        CONTRAINTES ELEMENTAIRES D UN ELEMENT */
/*  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*   PARAMETRES D'ENTREE : */
/*  ------------------- */
/*  NNO          : nombre de points (interpolation des fonctions) */
/*  NPO          : nombre de points (interpolation de la geometrie) */
/*  IJT          : PERMUTATION */
/*  X,Y,Z        : COORDONNEES DES NOEUDS DE L ELEMENT */
/*  ELAS         : MATRICE D ELASTICITE DE L ELEMENT */
/*  VDPQ2        : VALEURS DES DERIVES DES POL DE BASE (Fonctions) */
/*  VDPQ1        : VALEURS DES DERIVES DES POLYNOMES DE BASE (geometrie) */
/*  VP1          : VALEURS DES POLYNOMES DE BASE DE (geometrie) */
/*  POUR CES POL. OU DP. CES VALEURS  SONT AUX NOEUDS */
/*  U            : solution dans l element */

/*   TAB DE TRAVAIL */
/*  ---------------- */
/*  SIGMAE         : VALEURS DES CONTRAINTES ELEMENTAIRES */

/*  PARAMETRES RESULTATS : */
/*  -------------------- */
/*  SIGMA          : contraintes */
/*  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  PROGRAMMEURS  : Marina Vidrascu 2001 */
/*  ................................................................... */

/*     ---- non utilises --- */


/*   -----  CALCUL DE DELTA , DERIVEES  ----- */

    /* Parameter adjustments */
    a2_dim2 = *nno;
    a2_offset = (a2_dim2 + 1) * 3 + 1;
    a2 -= a2_offset;
    --delta;
    --u;
    sigma -= 7;
    sigmae_dim2 = *nno * 3;
    sigmae_offset = (sigmae_dim2 + 1) * 6 + 1;
    sigmae -= sigmae_offset;
    vdpq2_dim2 = *nno;
    vdpq2_offset = (vdpq2_dim2 + 1) * 3 + 1;
    vdpq2 -= vdpq2_offset;
    --ijt;
    vdpq1_dim2 = *npo;
    vdpq1_offset = (vdpq1_dim2 + 1) * 3 + 1;
    vdpq1 -= vdpq1_offset;
    --z__;
    --y;
    --x;
    --elas;

    /* Function Body */
    indice = 1;
    fobase_(&c__3, &c__3, nno, npo, nno, vp1, vp1, &vdpq2[vdpq2_offset], &vdpq1[vdpq1_offset], &x[1], &y[1], &z__[1], &a2[a2_offset], xint, yint, zint, &delta[1], dfinv, df, &indice, fi, fi, fi, iopt);


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

    i__1 = *nno;
    for (l = 1; l <= i__1; ++l) {

/*        B = DELTA * (DF)-1 * VDPQ2(L) */

	i__2 = *nno;
	for (j = 1; j <= i__2; ++j) {

/*   --- BLOC1 */

	    j1 = ijt[j];
	    sigmae[(j1 + l * sigmae_dim2) * 6 + 1] = elas[1] * a2[(j + l * a2_dim2) * 3 + 1] / delta[l];
	    sigmae[(j1 + l * sigmae_dim2) * 6 + 2] = elas[2] * a2[(j + l * a2_dim2) * 3 + 1] / delta[l];
	    sigmae[(j1 + l * sigmae_dim2) * 6 + 3] = elas[4] * a2[(j + l * a2_dim2) * 3 + 1] / delta[l];
	    sigmae[(j1 + l * sigmae_dim2) * 6 + 4] = elas[7] * a2[(j + l * a2_dim2) * 3 + 2] / delta[l];
	    sigmae[(j1 + l * sigmae_dim2) * 6 + 5] = 0.;
	    sigmae[(j1 + l * sigmae_dim2) * 6 + 6] = elas[9] * a2[(j + l * a2_dim2) * 3 + 3] / delta[l];

/*   --- BLOC2 */

	    j1 = ijt[j + *nno];
	    sigmae[(j1 + l * sigmae_dim2) * 6 + 1] = elas[2] * a2[(j + l * a2_dim2) * 3 + 2] / delta[l];
	    sigmae[(j1 + l * sigmae_dim2) * 6 + 2] = elas[3] * a2[(j + l * a2_dim2) * 3 + 2] / delta[l];
	    sigmae[(j1 + l * sigmae_dim2) * 6 + 3] = elas[5] * a2[(j + l * a2_dim2) * 3 + 2] / delta[l];
	    sigmae[(j1 + l * sigmae_dim2) * 6 + 4] = elas[7] * a2[(j + l * a2_dim2) * 3 + 1] / delta[l];
	    sigmae[(j1 + l * sigmae_dim2) * 6 + 5] = elas[8] * a2[(j + l * a2_dim2) * 3 + 3] / delta[l];
	    sigmae[(j1 + l * sigmae_dim2) * 6 + 6] = 0.;

/*   --- BLOC3 */

	    j1 = ijt[j + (*nno << 1)];
	    sigmae[(j1 + l * sigmae_dim2) * 6 + 1] = elas[4] * a2[(j + l * a2_dim2) * 3 + 3] / delta[l];
	    sigmae[(j1 + l * sigmae_dim2) * 6 + 2] = elas[5] * a2[(j + l * a2_dim2) * 3 + 3] / delta[l];
	    sigmae[(j1 + l * sigmae_dim2) * 6 + 3] = elas[6] * a2[(j + l * a2_dim2) * 3 + 3] / delta[l];
	    sigmae[(j1 + l * sigmae_dim2) * 6 + 4] = 0.;
	    sigmae[(j1 + l * sigmae_dim2) * 6 + 5] = elas[8] * a2[(j + l * a2_dim2) * 3 + 2] / delta[l];
	    sigmae[(j1 + l * sigmae_dim2) * 6 + 6] = elas[9] * a2[(j + l * a2_dim2) * 3 + 1] / delta[l];
/* L11: */
	}

/*       multiplier la contraine elementaire par le deplacement */

	for (i__ = 1; i__ <= 6; ++i__) {
	    sigma[i__ + l * 6] = (float)0.;
	    i__2 = *nno * 3;
	    for (j = 1; j <= i__2; ++j) {
		sigma[i__ + l * 6] += sigmae[i__ + (j + l * sigmae_dim2) * 6] * u[j];
/* L2: */
	    }
	}
/* L10: */
    }
} /* ec3c2c_ */

