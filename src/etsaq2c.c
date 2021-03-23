/* etsaq2c.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__8 = 8;
static int32 c__9 = 9;
static int32 c__1 = 1;
static int32 c__2 = 2;
static int32 c__16 = 16;

/* Subroutine */ int etsaq2c_(coor, fomega, fgamma, pressi, norefs, alpha, theta, car, iopt, be)
doublereal *coor, *fomega, *fgamma, *pressi;
int32 *norefs;
doublereal *alpha, *theta, *car;
int32 *iopt;
doublereal *be;
{
    /* Initialized data */

    static doublereal d2pi = 6.2831853071795862;
    static doublereal poids[9] = { .077160493827160503,.12345679012345678,.077160493827160503,.12345679012345678,.19753086419753085,.12345679012345678,.077160493827160503,.12345679012345678,.077160493827160503 };
    static doublereal p25[72]	/* was [8][9] */ = { .43237900077244512,-.1,-.032379000772445008,-.099999999999999991,.35491933384829665,.045080666151703314,.045080666151703314,.35491933384829665,-.099999999999999977,-.099999999999999977,-.099999999999999977,-.1,.8872983346207417,.19999999999999998,.11270166537925829,.19999999999999998,-.10000000000000003,.43237900077244512,-.099999999999999977,-.032379000772445015,.35491933384829665,.35491933384829665,.045080666151703308,.045080666151703314,
	    -.099999999999999991,-.099999999999999991,-.099999999999999991,-.099999999999999991,.19999999999999998,.11270166537925829,.19999999999999998,.8872983346207417,-.25,-.25,-.25,-.25,.5,.5,.5,.5,-.10000000000000019,-.099999999999999977,-.099999999999999977,-.099999999999999977,.19999999999999995,.8872983346207417,.19999999999999995,.11270166537925829,-.099999999999999977,-.032379000772444987,-.1,.43237900077244512,.045080666151703308,.045080666151703308,.35491933384829665,
	    .3549193338482966,-.10000000000000019,-.099999999999999977,-.099999999999999977,-.099999999999999977,.11270166537925829,.19999999999999995,.8872983346207417,.19999999999999995,-.032379000772445598,-.099999999999999866,.43237900077244484,-.10000000000000008,.045080666151703141,.35491933384829676,.35491933384829676,.045080666151703141 };
    static doublereal dp25[144]	/* was [2][8][9] */ = { -2.061895003862225,-2.061895003862225,-.68729833462074174,-.087298334620741685,-.26189500386222502,-.26189500386222502,-.087298334620741685,-.68729833462074174,2.7491933384829669,-.39999999999999996,.39999999999999996,.34919333848296674,.34919333848296674,.39999999999999996,-.39999999999999996,2.7491933384829669,-.68729833462074152,-.7745966692414834,.68729833462074174,-.7745966692414834,-.087298334620741657,-.7745966692414834,
	    .087298334620741685,-.7745966692414834,0.,-1.,.39999999999999996,1.5491933384829668,0.,1.,-.39999999999999996,1.5491933384829668,.68729833462074196,-.087298334620742101,2.0618950038622254,-2.061895003862225,.087298334620741713,-.68729833462074174,.26189500386222508,-.26189500386222497,-2.7491933384829669,-.39999999999999991,.39999999999999996,2.7491933384829669,-.34919333848296674,.39999999999999991,-.39999999999999996,.34919333848296674,-.7745966692414834,-.68729833462074174,
	    -.7745966692414834,.087298334620741685,-.7745966692414834,-.087298334620741685,-.7745966692414834,.68729833462074174,1.5491933384829668,-.39999999999999996,1.,0.,1.5491933384829668,.39999999999999996,-1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-1.,1.,0.,0.,1.,-1.,0.,.7745966692414834,.087298334620741213,.7745966692414834,-.68729833462074174,.7745966692414834,.68729833462074174,.7745966692414834,-.087298334620741657,-1.5491933384829668,-.39999999999999991,1.,0.,-1.5491933384829668,
	    .39999999999999991,-1.,0.,-.087298334620742157,.68729833462074174,-.26189500386222502,.26189500386222502,-.68729833462074174,.087298334620741685,-2.0618950038622254,2.061895003862225,.34919333848296674,-.39999999999999996,.39999999999999991,-.34919333848296674,2.7491933384829669,.39999999999999996,-.39999999999999991,-2.7491933384829669,.087298334620741213,.7745966692414834,-.087298334620741657,.7745966692414834,.68729833462074174,.7745966692414834,-.68729833462074196,
	    .7745966692414834,0.,-1.,.39999999999999991,-1.5491933384829668,0.,1.,-.39999999999999991,-1.5491933384829668,.26189500386222475,.26189500386222452,.087298334620741879,.68729833462074174,2.061895003862225,2.061895003862225,.68729833462074152,.087298334620741657,-.34919333848296663,-.39999999999999991,.39999999999999991,-2.7491933384829669,-2.7491933384829669,.39999999999999991,-.39999999999999991,-.34919333848296663 };
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
    static doublereal a[4], fface[18]	/* was [2][9] */, d__;
    static int32 i__, j, k, l, m;
    extern /* Subroutine */ int e2aq2c_();
    static doublereal d1, d2, f1[9], f2[9], g1[16];
    static int32 i1, i2, k1, l1;
    static doublereal g2[16];
    static int32 k4;
    static doublereal g4[128]	/* was [16][8] */, dfm1dp[144]	/* was [2][8][9] */;
    static int32 il, ip[16];
    static doublereal fx, fy, rr, farete[24]	/* was [2][12] */, poidel[9], arelon;
    extern /* Subroutine */ int hookax_();
    static doublereal p25a[9]	/* was [3][3] */, alpha_r__, alpha_t__, alpha_z__, xmi, ymi, xjm, yjm;
    extern /* Subroutine */ int ab1d_();
    static doublereal xnu[3], ynu[3], dfm1[36]	/* was [4][9] */, p1pa, p2pa, p3pa;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* BUT: CALCUL DE SECONDS MEMBRES DE L ELEMENT QUAD AQ2C */
/* --- */
/* in : */
/*   coor(noe,ndim) : coordonnees R(6), Z(6) des 8 noeuds. */
/*   fomega(ndim,noe)        : fx, fy aux noeuds => fomega(ndim,npi) */
/*                               => fface(ndim=2 , npi=1) */
/*   fgamma(ndim,nnof*nbarete): fx, fy aux noeuds de chaque arete. */
/*   pressi(nnof*nbarete)     : pression aux noeuds de chaque arete */
/*                               => farete(ndim=2  ,npia*nbarete=6) */
/*   norefs(nbarete,2): norefs(i,1) = 0 si fgamma   = 0 sur arete_i */
/*                      norefs(i,2) = 0 si pression = 0 sur arete_i */
/*   theta(8) : theta aux noeuds */
/*   alpha    : coef. dilatation thermique: radial, axial, tangentiel */
/*   car      : caracteristiques des materiaux */
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
/* out: BE(16): second membre. */
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
    coor -= 9;

    /* Function Body */
/* 2Q25 -- XYNPI: coordonnees pt. int. numeriques (element reference) */
/*     -- POIDS: poids du schema d'integration numerique. */
/*     -- Valeurs des Polynomes de base aux pt. int. numerique. */
/*     -- Valeurs Derive'es des Poly. de base aux pt int. numerique. */

/*   -- Valeurs aux pt d'int. num. a partir valeurs aux noeuds */
/*      Efforts volumiques  fomega(ndim,noe) -> fface(ndim,npi) */
    for (i__ = 1; i__ <= 9; ++i__) {
	fface[(i__ << 1) - 2] = 0.;
	fface[(i__ << 1) - 1] = 0.;
	for (j = 1; j <= 8; ++j) {
	    fface[(i__ << 1) - 2] += p25[j + (i__ << 3) - 9] * fomega[(j << 1) + 1];
	    fface[(i__ << 1) - 1] += p25[j + (i__ << 3) - 9] * fomega[(j << 1) + 2];
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
	m = k + 4;
	xmi = coor[m + 8] - coor[k + 8];
	ymi = coor[m + 16] - coor[k + 16];
/* Computing 2nd power */
	d__1 = xmi;
/* Computing 2nd power */
	d__2 = ymi;
	arelon = sqrt(d__1 * d__1 + d__2 * d__2);
	xnu[0] = ymi / arelon;
	ynu[0] = -xmi / arelon;
	xjm = coor[j + 8] - coor[m + 8];
	yjm = coor[j + 16] - coor[m + 16];
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
	    p25a[i__ * 3 - 1] = xa[i__ - 1] * 4. * (1. - xa[i__ - 1]);
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

/*     F1 , F2 , DFM1 , POIDEL , IP , DP PRETS A L EMPLOI */

    e2aq2c_(&c__8, &c__9, poids, p25, dp25, ip, f1, f2, dfm1dp, poidel, dfm1, &coor[9]);
    for (i__ = 1; i__ <= 16; ++i__) {
	be[i__] = 0.;
/* L6: */
    }

/*     CONTRIBUTION DES EFFORTS SURFACIQUES */

    for (i__ = 1; i__ <= 8; ++i__) {
	for (j = 1; j <= 2; ++j) {
	    i1 = ip[i__ + (j - 1 << 3) - 1];
	    d1 = 0.;
	    for (l = 1; l <= 9; ++l) {
		d1 += poidel[l - 1] * p25[i__ + (l << 3) - 9] * fface[j + (l << 1) - 3];
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
	k4 = k + 4;
	i__1 = npia;
	for (l = 1; l <= i__1; ++l) {
	    rr = coor[k + 8] * p25a[l * 3 - 3] + coor[k1 + 8] * p25a[l * 3 - 2] + coor[k4 + 8] * p25a[l * 3 - 1];
	    p1pa = xa[l - 1] * 4. - 3.;
	    p2pa = xa[l - 1] * 4. - 1.;
	    p3pa = 4. - xa[l - 1] * 8.;
	    d1 = coor[k + 8] * p1pa + coor[k1 + 8] * p2pa + coor[k4 + 8] * p3pa;
	    d2 = coor[k + 16] * p1pa + coor[k1 + 16] * p2pa + coor[k4 + 16] * p3pa;
	    d1 = d1 * d1 + d2 * d2;
	    d1 = d2pi * sqrt(d1);
	    l1 = l + npia * (k - 1);
	    for (j = 1; j <= 2; ++j) {
		i__ = ip[k + (j - 1 << 3) - 1];
		i1 = ip[k1 + (j - 1 << 3) - 1];
		i2 = ip[k + 4 + (j - 1 << 3) - 1];
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

/* [Sig_zz]   [    ] [alpha_z]                               [theta(1)] */
/* [      ]   [    ] [       ]                               [theta(2)] */
/* [Sig_rr]   [    ] [alpha_r]                               [theta(3)] */
/* [      ] = [elas]*[       ]*[ q1 q2 q3 q4 q5 q6 q7 q8 ] * [theta(4)] */
/* [Sig_tt]   [    ] [alpha_t]                          1*8  [theta(5)] */
/* [      ]   [    ] [       ]                               [theta(6)] */
/* [Sig_rz]   [    ] [  0    ]                               [theta(7)] */
/*       4*1      4*4       4*1                              [theta(8)] */
/*                                                                  8*1 */
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
	for (i__ = 1; i__ <= 16; ++i__) {
	    g4[i__ + (l << 4) - 17] = 0.;
/* L20: */
	}
    }

    for (l = 1; l <= 9; ++l) {
	d__ = poidel[l - 1];
	for (i__ = 2; i__ <= 4; ++i__) {
	    g1[i__ - 1] = d__ * a[i__ - 1];
/* L13: */
	}
	d__ /= f1[l - 1];
	g1[0] = d__ * a[0];

/*       -- TP * G1(1) + TDFM1DP * (G1(2),G1(3)) */

	tab0d_(&c__8, &c__1, &c__1, &p25[(l << 3) - 8], g1, g2);
	tab1d_(&c__8, &c__2, &c__1, &dfm1dp[((l << 3) + 1 << 1) - 18], &g1[1], g2);

/*       -- TDFM1DP * (G1(3),G1(4)) */

	tab0d_(&c__8, &c__2, &c__1, &dfm1dp[((l << 3) + 1 << 1) - 18], &g1[2], &g2[8]);
/*        G2 * P */
	ab1d_(&c__16, &c__1, &c__8, g2, &p25[(l << 3) - 8], g4);
/* L14: */
    }

/*     Theta * G4  => BE = BE + G4(IP) */

    for (i__ = 1; i__ <= 16; ++i__) {
	l = ip[i__ - 1];
	d__ = 0.;
	for (j = 1; j <= 8; ++j) {
	    d__ += theta[j] * g4[i__ + (j << 4) - 17];
/* L15: */
	}
	be[l] += d__;
/* L16: */
    }
} /* etsaq2c_ */

