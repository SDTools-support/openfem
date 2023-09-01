/* etsaq1d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__4 = 4;
static int32 c__1 = 1;
static int32 c__2 = 2;
static int32 c__8 = 8;

/* Subroutine */ int etsaq1d_(coor, fomega, fgamma, pressi, norefs, alpha, theta, car, iopt, be)
doublereal *coor, *fomega, *fgamma, *pressi;
int32 *norefs;
doublereal *alpha, *theta, *car;
int32 *iopt;
doublereal *be;
{
    /* Initialized data */

    static doublereal d2pi = 6.2831853071795862;
    static doublereal poids[4] = { .25,.25,.25,.25 };
    static doublereal q13[16]	/* was [4][4] */ = { .62200846792814612,.16666666666666662,.044658198738520435,.16666666666666662,.16666666666666662,.62200846792814623,.16666666666666662,.044658198738520449,.044658198738520504,.16666666666666662,.62200846792814623,.16666666666666662,.16666666666666668,.044658198738520449,.16666666666666662,.62200846792814623 };
    static doublereal dq13[32]	/* was [2][4][4] */ = { -.78867513459481286,-.78867513459481286,.78867513459481286,-.21132486540518707,.21132486540518707,.21132486540518707,-.21132486540518707,.78867513459481286,-.78867513459481286,-.21132486540518713,.78867513459481286,-.78867513459481286,.21132486540518707,.78867513459481286,-.21132486540518707,.21132486540518713,-.21132486540518713,-.21132486540518713,.21132486540518713,-.78867513459481286,.78867513459481286,.78867513459481286,
	    -.78867513459481286,.21132486540518713,-.21132486540518713,-.78867513459481286,.21132486540518713,-.21132486540518707,.78867513459481286,.21132486540518707,-.78867513459481286,.78867513459481286 };
    static int32 npia = 2;
    static doublereal poidsa[2] = { .5,.5 };
    static doublereal xa[2] = { .211324865404,.788675134593 };

    /* System generated locals */
    int32 i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal elas[10];
    extern /* Subroutine */ int tab0d_(), tab1d_();
    static doublereal a[4], fface[8]	/* was [2][4] */, d__;
    static int32 i__, j, k, l;
    extern /* Subroutine */ int e1aq1c_();
    static doublereal d1, f1[5], f2[5], g1[8];
    static int32 i1;
    static doublereal g2[8];
    static int32 k1, l1;
    static doublereal g4[32]	/* was [8][4] */, dfm1dp[32]	/* was [2][4][4] */;
    static int32 il, ip[8];
    static doublereal fx, fy, rr, farete[16]	/* was [2][8] */, poidel[4], arelon[4];
    extern /* Subroutine */ int hookax_();
    static doublereal q13a[4]	/* was [2][2] */, alpha_r__, alpha_t__, xji[4], yji[4], alpha_z__;
    extern /* Subroutine */ int ab1d_();
    static doublereal xnu, ynu, dfm1[16]	/* was [4][4] */;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* BUT: CALCUL DE SECONDS MEMBRES DE L ELEMENT TRIA AQ1D */
/* --- */
/* in : */
/*   coor(noe,ndim) : coordonnees R(4), Z(4) des 4 sommets. */
/*   fomega(ndim,noe)        : fx, fy aux noeuds => fomega(ndim,npi) */
/*                               => fface(ndim=2 , npi=4) */
/*   fgamma(ndim,noe*nbarete): fx, fy aux noeuds de chaque arete. */
/*   pressi(nnof*nbarete)    : pression aux noeuds de chaque arete */
/*                               => farete(ndim=2  ,npia*nbarete=6) */
/*   norefs(nbarete,2): norefs(i,1) = 0 si fgamma   = 0 sur arete_i */
/*                      norefs(i,2) = 0 si pression = 0 sur arete_i */
/*   theta(noe) : theta aux 4 sommets */
/*   alpha    : coef. dilatation thermique: radial, axial, tangentiel */
/*   car        : caracteristiques des materiaux */
/*        if(iopt .eq. 1) then */
/*          car(1) = young */
/*          car(2) = poisson */
/*        else */
/*          car(1) = E_1  (Young radial     -> E_r     ) */
/*          car(2) = nu_1 (poisson          -> Nu_r    ) */
/*          car(3) = E_2  (Young axial      -> E_z     ) */
/*          car(4) = nu_2 (poisson          -> Nu_z    ) */
/*          car(5) = E_3  (Young Tangentiel -> E_theta ) */
/*        end if */
/* out: BE(8): second membre. */
/*  .................................................................. */
    /* Parameter adjustments */
    --be;
    --car;
    --theta;
    --alpha;
    norefs -= 5;
    --pressi;
    fgamma -= 3;
    fomega -= 3;
    coor -= 5;

    /* Function Body */
/* 2Q13 -- POIDS: poids du schema d'integration numerique. */
/*     -- XYNPI: coordonnees pt. int. numeriques (element reference) */
/*     -- Valeurs des Polynomes de base aux pt. int. numerique. */
/*     -- Valeurs Derive'es des Poly. de base aux pt int. numerique. */
/*     -- integration / aretes */
/*   -- Valeurs aux pt d'int. num. a partir valeurs aux noeuds */
/*      Efforts volumiques  fomega(ndim,noe) -> fface(ndim,npi) */
    for (i__ = 1; i__ <= 4; ++i__) {
	fface[(i__ << 1) - 2] = 0.;
	fface[(i__ << 1) - 1] = 0.;
	for (j = 1; j <= 4; ++j) {
	    fface[(i__ << 1) - 2] += q13[j + (i__ << 2) - 5] * fomega[(j << 1) + 1];
	    fface[(i__ << 1) - 1] += q13[j + (i__ << 2) - 5] * fomega[(j << 1) + 2];
/* L1: */
	}
/* L2: */
    }

/*   -- Valeurs aux pt d'int. num. a partir valeurs aux noeuds */
/*      fgamma(ndim,(nnof*nbarete) -> farete(ndim, npia*nbarete) */
/*      pressi(nnof*nbarete)       -> idem */

    for (k = 1; k <= 4; ++k) {
/*       -- LONGEUR DE L ARETE, COSINUS DIRECTEURS DE LA NORMALE ext. */
	j = k % 4 + 1;
	xji[k - 1] = coor[j + 4] - coor[k + 4];
	yji[k - 1] = coor[j + 8] - coor[k + 8];
/* Computing 2nd power */
	d__1 = xji[k - 1];
/* Computing 2nd power */
	d__2 = yji[k - 1];
	arelon[k - 1] = sqrt(d__1 * d__1 + d__2 * d__2);
	xnu = yji[k - 1] / arelon[k - 1];
	ynu = -xji[k - 1] / arelon[k - 1];
	il = npia * (k - 1);
	i__1 = npia;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    fx = 0.;
	    fy = 0.;
	    q13a[(i__ << 1) - 2] = 1. - xa[i__ - 1];
	    q13a[(i__ << 1) - 1] = xa[i__ - 1];
	    for (j = 1; j <= 2; ++j) {
		fx += q13a[j + (i__ << 1) - 3] * (fgamma[(j + il << 1) + 1] - pressi[j + il] * xnu);
		fy += q13a[j + (i__ << 1) - 3] * (fgamma[(j + il << 1) + 2] - pressi[j + il] * ynu);
/* L3: */
	    }
	    farete[(il + i__ << 1) - 2] = fx;
	    farete[(il + i__ << 1) - 1] = fy;
/* L4: */
	}
/* L5: */
    }

/*     -- CALCUL DE F1,F2,FFM1,DFM1DP */

    e1aq1c_(&c__4, &c__4, poids, q13, dq13, ip, f1, f2, dfm1dp, poidel, dfm1, &coor[5]);

    for (i__ = 1; i__ <= 8; ++i__) {
	be[i__] = 0.;
/* L6: */
    }

/*     CONTRIBUTION DES EFFORTS SURFACIQUES */
/*             FOMEGA(2,npi=4) */
    for (i__ = 1; i__ <= 4; ++i__) {
	for (j = 1; j <= 2; ++j) {
	    i1 = ip[i__ + (j - 1 << 2) - 1];
	    d1 = 0.;
	    for (l = 1; l <= 4; ++l) {
		d1 += poidel[l - 1] * q13[i__ + (l << 2) - 5] * fface[j + (l << 1) - 3];
/* L7: */
	    }
	    be[i1] += d1;
/* L8: */
	}
/* L9: */
    }

/*     -- CONTRIBUTIONS DES EFFORTS SUR LES ARETES */

    for (k = 1; k <= 4; ++k) {
	k1 = k % 4 + 1;
	i__1 = npia;
	for (l = 1; l <= i__1; ++l) {
	    d1 = d2pi * arelon[k - 1];
	    rr = coor[k + 4] * (1 - xa[l - 1]) + coor[k1 + 4] * xa[l - 1];
	    l1 = l + npia * (k - 1);
	    for (j = 1; j <= 2; ++j) {
		i__ = ip[k + (j - 1 << 2) - 1];
		i1 = ip[k1 + (j - 1 << 2) - 1];
		d__ = rr * d1 * poidsa[l - 1] * farete[j + (l1 << 1) - 3];
		be[i__] += d__ * (1 - xa[l - 1]);
		be[i1] += d__ * xa[l - 1];
/* L10: */
	    }
/* L11: */
	}
/* L12: */
    }

/*     CONTRIBUTIONS DES CONTRAINTES THERMIQUES  E(10) */

    hookax_(iopt, &car[1], elas);

/*  [  Sig_zz  ]    [      ] [alpha_z] */
/*  [          ]    [      ] [       ] */
/*  [  Sig_rr  ]    [      ] [alpha_r]                [ theta(1) ] */
/*  [          ] := [ elas ]*[       ]*[ p1 p2 p3 ] * [ theta(2) ] */
/*  [  Sig_tt  ]    [      ] [alpha_t]           1*3  [ theta(3) ] */
/*  [          ]    [      ] [       ]                [ theta(4) ] */
/*  [  Sig_rz  ]    [      ] [  0    ]                          4*1 */
/*            4*1         4*4       4*1 */

/* -- [ A ]  = [ELAS] * [ALPHA] */
/*       4*1       4*4       4*1 */
    alpha_r__ = alpha[1];
    alpha_z__ = alpha[2];
    alpha_t__ = alpha[3];
    a[0] = alpha_z__ * elas[3] + alpha_r__ * elas[4] + alpha_t__ * elas[5];
    a[1] = alpha_z__ * elas[1] + alpha_r__ * elas[2] + alpha_t__ * elas[4];
    a[2] = alpha_z__ * elas[6] + alpha_r__ * elas[7] + alpha_t__ * elas[8];
    a[3] = alpha_z__ * elas[0] + alpha_r__ * elas[1] + alpha_t__ * elas[3];

    for (l = 1; l <= 8; ++l) {
	for (i__ = 1; i__ <= 4; ++i__) {
	    g4[l + (i__ << 3) - 9] = 0.;
/* L14: */
	}
    }

    for (l = 1; l <= 4; ++l) {
	d__ = poidel[l - 1];
	for (i__ = 2; i__ <= 4; ++i__) {
	    g1[i__ - 1] = d__ * a[i__ - 1];
/* L15: */
	}
	d__ /= f1[l - 1];
	g1[0] = d__ * a[0];

/*       -- TP * G1(1) + TDDFM1P * (G1(2),G1(3)) */

	tab0d_(&c__4, &c__1, &c__1, &q13[(l << 2) - 4], g1, g2);
	tab1d_(&c__4, &c__2, &c__1, &dfm1dp[((l << 2) + 1 << 1) - 10], &g1[1], g2);

/*       -- TDFM1DP * (G1(2),G1(4)) */

	tab0d_(&c__4, &c__2, &c__1, &dfm1dp[((l << 2) + 1 << 1) - 10], &g1[2], &g2[4]);
/*       -- G2 * P */
	ab1d_(&c__8, &c__1, &c__4, g2, &q13[(l << 2) - 4], g4);
/* L16: */
    }

/*     THETA * G4  => BE = BE + G4(IP) */

    for (i__ = 1; i__ <= 8; ++i__) {
	l = ip[i__ - 1];
	d__ = 0.;
	for (j = 1; j <= 4; ++j) {
	    d__ += theta[j] * g4[i__ + (j << 3) - 9];
/* L17: */
	}
	be[l] += d__;
/* L18: */
    }
} /* etsaq1d_ */

