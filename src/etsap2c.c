/* etsap2c.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__6 = 6;
static int32 c__7 = 7;
static int32 c__1 = 1;
static int32 c__2 = 2;
static int32 c__12 = 12;

/* Subroutine */ int etsap2c_(coor, fomega, fgamma, pressi, norefs, alpha, theta, car, iopt, be)
doublereal *coor, *fomega, *fgamma, *pressi;
int32 *norefs;
doublereal *alpha, *theta, *car;
int32 *iopt;
doublereal *be;
{
    /* Initialized data */

    static doublereal d2pi = 6.2831853071795862;
    static doublereal poids[7] = { .062969590272413583,.062969590272413583,.062969590272413583,.066197076394253095,.066197076394253095,.066197076394253095,.11249999701976776 };
    static doublereal p25[42]	/* was [6][7] */ = { .47435260858553857,-.080768594191887185,-.080768594191887185,.32307437676754874,.041035826263138293,.32307437676754874,-.080768594191886977,.4743526085855384,-.080768594191887185,.32307437676754879,.32307437676754868,.041035826263138341,-.080768594191887033,-.080768594191887185,.4743526085855384,.041035826263138334,.32307437676754868,.3230743767675489,-.028074943223078796,-.028074943223078852,-.052583901102545349,.88413424176407262,
	    .11229977289231514,.11229977289231517,-.052583901102545571,-.028074943223078852,-.028074943223078852,.1122997728923154,.8841342417640724,.1122997728923154,-.028074943223078765,-.052583901102545349,-.028074943223078852,.11229977289231514,.11229977289231514,.88413424176407262,-.11111111773384862,-.11111110779974175,-.11111110779974175,.44444443119896703,.44444447093539807,.44444443119896703 };
    static doublereal dp25[84]	/* was [2][6][7] */ = { -2.1897079414123492,-2.1897079414123492,-.59485397070617462,0.,0.,-.59485397070617462,2.7845619121185238,-.40514602929382531,.40514602929382531,.40514602929382531,-.40514602929382531,2.7845619121185238,.59485397070617418,.59485397070617418,2.1897079414123488,0.,0.,-.59485397070617462,-2.7845619121185229,-3.1897079414123488,.40514602929382531,3.1897079414123488,-.40514602929382531,0.,.59485397070617418,.59485397070617418,-.59485397070617462,
	    0.,0.,2.1897079414123488,0.,-.40514602929382531,3.1897079414123488,.40514602929382531,-3.1897079414123488,-2.7845619121185229,-.88056825642046065,-.88056825642046065,.88056825642046021,0.,0.,-.76113651284092076,0.,-1.8805682564204602,.2388634871590792,1.8805682564204602,-.2388634871590792,1.6417047692613815,.76113651284092043,.76113651284092043,.88056825642046021,0.,0.,.88056825642046021,-1.6417047692613806,-1.8805682564204602,1.8805682564204602,1.8805682564204602,-1.8805682564204602,
	    -1.6417047692613806,-.88056825642046054,-.88056825642046054,-.76113651284092076,0.,0.,.88056825642046021,1.6417047692613813,-.2388634871590792,1.8805682564204602,.2388634871590792,-1.8805682564204602,0.,-.33333325386047363,-.33333325386047363,.33333337306976318,0.,0.,.33333337306976318,-1.1920928955078125e-7,-1.3333333730697631,1.3333333730697631,1.3333333730697631,-1.3333333730697631,-1.1920928955078125e-7 };
    static int32 npia = 3;
    static doublereal poidsa[3] = { .277777777777,.444444444444,.277777777777 };
    static doublereal xa[3] = { .1127015,.5,.8872985 };

    /* System generated locals */
    int32 i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal elas[10];
    extern /* Subroutine */ int tab0d_(), tab1d_();
    static doublereal a[4], fface[14]	/* was [2][7] */, d__;
    static int32 i__, j, k, l, m;
    extern /* Subroutine */ int e2ap2c_();
    static doublereal d1, d2, f1[7], f2[7], g1[12];
    static int32 i1, i2, k1, l1;
    static doublereal g2[12];
    static int32 k4;
    static doublereal g4[72]	/* was [12][6] */, dfm1dp[84]	/* was [2][6][7] */;
    static int32 il, ip[12];
    static doublereal fx, fy, rr, farete[18]	/* was [2][9] */, poidel[7], arelon;
    extern /* Subroutine */ int hookax_();
    static doublereal p25a[9]	/* was [3][3] */, alpha_r__, alpha_t__, alpha_z__, xmi, ymi, xjm, yjm;
    extern /* Subroutine */ int ab1d_();
    static doublereal xnu[3], ynu[3], dfm1[28]	/* was [4][7] */, p1pa, p2pa, p3pa;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* BUT: CALCUL DE SECONDS MEMBRES DE L ELEMENT TRIA AP2C */
/* --- */
/* in : */
/*   coor(noe,ndim) : coordonnees R(6), Z(6) des 6 noeuds. */
/*   fomega(ndim,noe)        : fx, fy aux noeuds => fomega(ndim,npi) */
/*                               => fface(ndim=2 , npi=1) */
/*   fgamma(ndim,nnof*nbarete): fx, fy aux noeuds de chaque arete. */
/*   pressi(nnof*nbarete)     : pression aux noeuds de chaque arete */
/*                               => farete(ndim=2  ,npia*nbarete=6) */
/*   norefs(nbarete,2): norefs(i,1) = 0 si fgamma   = 0 sur arete_i */
/*                      norefs(i,2) = 0 si pression = 0 sur arete_i */
/*   theta(6)    : theta aux 6 noeuds */
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
/* out: BE(12): second membre. */
/*  .................................................................. */
    /* Parameter adjustments */
    --be;
    --car;
    --theta;
    --alpha;
    norefs -= 4;
    --pressi;
    fgamma -= 3;
    fomega -= 3;
    coor -= 7;

    /* Function Body */
/* 2P25 -- XYNPI: coordonnees pt. int. numeriques (element reference) */
/*     -- POIDS: poids du schema d'integration numerique. */
/*     -- Valeurs des Polynomes de base aux pt. int. numerique. */
/*     -- Valeurs Derive'es des Poly. de base aux pt int. numerique. */

/*   -- Valeurs aux pt d'int. num. a partir valeurs aux noeuds */
/*      Efforts volumiques  fomega(ndim,noe) -> fface(ndim,npi) */
    for (i__ = 1; i__ <= 7; ++i__) {
	fface[(i__ << 1) - 2] = 0.;
	fface[(i__ << 1) - 1] = 0.;
	for (j = 1; j <= 6; ++j) {
	    fface[(i__ << 1) - 2] += p25[j + i__ * 6 - 7] * fomega[(j << 1) + 1];
	    fface[(i__ << 1) - 1] += p25[j + i__ * 6 - 7] * fomega[(j << 1) + 2];
/* L1: */
	}
/* L2: */
    }
/*   -- Valeurs aux pt d'int. num. a partir valeurs aux noeuds */
/*      fgamma(ndim,(nnof*nbarete) -> farete(ndim, npia*nbarete) */
/*      pressi(nnof*nbarete)       -> idem */
    for (k = 1; k <= 3; ++k) {
/*       -- LONGEUR DE L ARETE, COSINUS DIRECTEURS DE LA NORMALE ext. */
	j = k % 3 + 1;
	m = k + 3;
	xmi = coor[m + 6] - coor[k + 6];
	ymi = coor[m + 12] - coor[k + 12];
/* Computing 2nd power */
	d__1 = xmi;
/* Computing 2nd power */
	d__2 = ymi;
	arelon = sqrt(d__1 * d__1 + d__2 * d__2);
	xnu[0] = ymi / arelon;
	ynu[0] = -xmi / arelon;
	xjm = coor[j + 6] - coor[m + 6];
	yjm = coor[j + 12] - coor[m + 12];
/* Computing 2nd power */
	d__1 = xjm;
/* Computing 2nd power */
	d__2 = yjm;
	arelon = sqrt(d__1 * d__1 + d__2 * d__2);
	xnu[1] = yjm / arelon;
	ynu[1] = -xjm / arelon;
	xnu[2] = (xnu[0] + xnu[1]) * .5;
	ynu[2] = (ynu[0] + ynu[1]) * .5;
	il = npia * (k - 1);
	i__1 = npia;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    fx = 0.;
	    fy = 0.;
	    p25a[i__ * 3 - 3] = (1. - xa[i__ - 1]) * (1. - xa[i__ - 1] * 2.);
	    p25a[i__ * 3 - 2] = xa[i__ - 1] * (xa[i__ - 1] * 2. - 1.);
	    p25a[i__ * 3 - 1] = (1. - xa[i__ - 1]) * 4. * xa[i__ - 1];
	    for (j = 1; j <= 3; ++j) {
		fx += p25a[j + i__ * 3 - 4] * (fgamma[(j + il << 1) + 1] - pressi[j + il] * xnu[j - 1]);
		fy += p25a[j + i__ * 3 - 4] * (fgamma[(j + il << 1) + 2] - pressi[j + il] * ynu[j - 1]);
/* L3: */
	    }
	    farete[(il + i__ << 1) - 2] = fx;
	    farete[(il + i__ << 1) - 1] = fy;
/* L4: */
	}
/* L5: */
    }

/*     -- CALCUL DE F1,F2,FFM1,DFM1DP */

    e2ap2c_(&c__6, &c__7, poids, p25, dp25, ip, f1, f2, dfm1dp, poidel, dfm1, &coor[7]);

    for (i__ = 1; i__ <= 12; ++i__) {
	be[i__] = 0.;
/* L6: */
    }

/*     -- CONTRIBUTION DES EFFORTS SURFACIQUES */
/*             FOMEGA(2,npi) */
    for (i__ = 1; i__ <= 6; ++i__) {
	for (j = 1; j <= 2; ++j) {
	    i1 = ip[i__ + (j - 1) * 6 - 1];
	    d1 = 0.;
	    for (l = 1; l <= 7; ++l) {
		d1 += poidel[l - 1] * p25[i__ + l * 6 - 7] * fface[j + (l << 1) - 3];
/* L7: */
	    }
	    be[i1] += d1;
/* L8: */
	}
/* L9: */
    }

/*     -- CONTRIBUTIONS DES EFFORTS SUR LES ARETES */

    for (k = 1; k <= 3; ++k) {
	k1 = k % 3 + 1;
	k4 = k + 3;
	i__1 = npia;
	for (l = 1; l <= i__1; ++l) {
	    rr = coor[k + 6] * p25a[l * 3 - 3] + coor[k1 + 6] * p25a[l * 3 - 2] + coor[k4 + 6] * p25a[l * 3 - 1];
	    p1pa = xa[l - 1] * 4. - 3.;
	    p2pa = xa[l - 1] * 4. - 1.;
	    p3pa = 4. - xa[l - 1] * 8.;
	    d1 = coor[k + 6] * p1pa + coor[k1 + 6] * p2pa + coor[k4 + 6] * p3pa;
	    d2 = coor[k + 12] * p1pa + coor[k1 + 12] * p2pa + coor[k4 + 12] * p3pa;
	    d1 = d1 * d1 + d2 * d2;
	    d1 = d2pi * sqrt(d1);
	    l1 = l + npia * (k - 1);
	    for (j = 1; j <= 2; ++j) {
		i__ = ip[k + (j - 1) * 6 - 1];
		i1 = ip[k1 + (j - 1) * 6 - 1];
		i2 = ip[k + 3 + (j - 1) * 6 - 1];
		d__ = rr * d1 * poidsa[l - 1] * farete[j + (l1 << 1) - 3];
		be[i__] += d__ * p25a[l * 3 - 3];
		be[i1] += d__ * p25a[l * 3 - 2];
		be[i2] += d__ * p25a[l * 3 - 1];
/* L10: */
	    }
/* L11: */
	}
/* L12: */
    }

/*     -- CONTRIBUTIONS DES CONTRAINTES THERMIQUES */

    hookax_(iopt, &car[1], elas);

/*  [  Sig_zz  ]    [      ] [alpha_z] */
/*  [          ]    [      ] [       ] */
/*  [  Sig_rr  ]    [      ] [alpha_r]                [ theta(1) ] */
/*  [          ] := [ elas ]*[       ]*[ p1 p2 p3 ] * [ theta(2) ] */
/*  [  Sig_tt  ]    [      ] [alpha_t]           1*3  [ theta(3) ] */
/*  [          ]    [      ] [       ]                          3*1 */
/*  [  Sig_rz  ]    [      ] [  0    ] */
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

    for (l = 1; l <= 6; ++l) {
	for (i__ = 1; i__ <= 12; ++i__) {
	    g4[i__ + l * 12 - 13] = 0.;
/* L15: */
	}
    }

    for (l = 1; l <= 7; ++l) {
	d__ = poidel[l - 1];
	for (i__ = 2; i__ <= 4; ++i__) {
	    g1[i__ - 1] = d__ * a[i__ - 1];
/* L16: */
	}
	d__ /= f1[l - 1];
	g1[0] = d__ * a[0];

/*       -- TP * G1(1) + TDFM1DP * (G1(2),G1(3)) */

	tab0d_(&c__6, &c__1, &c__1, &p25[l * 6 - 6], g1, g2);
	tab1d_(&c__6, &c__2, &c__1, &dfm1dp[(l * 6 + 1 << 1) - 14], &g1[1], g2);

/*       -- TDFM1DP * (G1(2),G1(4)) */

	tab0d_(&c__6, &c__2, &c__1, &dfm1dp[(l * 6 + 1 << 1) - 14], &g1[2], &g2[6]);
/*        G2 * P */
	ab1d_(&c__12, &c__1, &c__6, g2, &p25[l * 6 - 6], g4);
/* L17: */
    }

/*     THETA * G4  => BE = BE + G4(IP) */

    for (i__ = 1; i__ <= 12; ++i__) {
	l = ip[i__ - 1];
	d__ = 0.;
	for (j = 1; j <= 6; ++j) {
	    d__ += theta[j] * g4[i__ + j * 12 - 13];
/* L18: */
	}
	be[l] += d__;
/* L20: */
    }
} /* etsap2c_ */

