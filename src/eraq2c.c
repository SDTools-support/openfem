/* eraq2c.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__9 = 9;

/* Subroutine */ int eraq2c_(n, r__, z__, young, poisson, ae)
int32 *n;
doublereal *r__, *z__, *young, *poisson, *ae;
{
    /* Initialized data */

    static doublereal poidsg[9] = { .30864197530864201,.49382716049382713,.30864197530864201,.49382716049382713,.79012345679012341,.49382716049382713,.30864197530864201,.49382716049382713,.30864197530864201 };
    static doublereal p[81]	/* was [9][9] */ = { .4723790007724451,-.060000000000000004,.0076209992275549851,-.060000000000000004,.27491933384829669,-.034919333848296665,-.034919333848296665,.27491933384829669,.15999999999999997,0.,0.,0.,0.,.68729833462074174,0.,-.087298334620741685,0.,.39999999999999996,-.060000000000000004,.4723790007724451,-.060000000000000004,.0076209992275549851,.27491933384829669,.27491933384829669,-.034919333848296665,-.034919333848296665,.15999999999999997,0.,0.,0.,
	    0.,0.,-.087298334620741685,0.,.68729833462074174,.39999999999999996,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,.68729833462074174,0.,-.087298334620741685,.39999999999999996,-.060000000000000004,.0076209992275549851,-.060000000000000004,.4723790007724451,-.034919333848296665,-.034919333848296665,.27491933384829669,.27491933384829669,.15999999999999997,0.,0.,0.,0.,-.087298334620741685,0.,.68729833462074174,0.,.39999999999999996,.0076209992275549851,-.060000000000000004,.4723790007724451,
	    -.060000000000000004,-.034919333848296665,.27491933384829669,.27491933384829669,-.034919333848296665,.15999999999999997 };
    static doublereal dp[162]	/* was [2][9][9] */ = { -.87602816808281591,-.87602816808281591,-.18872983346207417,.11127016653792583,.023971831917184147,.023971831917184147,.11127016653792583,-.18872983346207417,1.0647580015448903,-.50983866769659336,-.10983866769659336,-.13524199845510997,-.13524199845510997,-.10983866769659336,-.50983866769659336,1.0647580015448903,.6196773353931867,.6196773353931867,-.34364916731037087,0.,.34364916731037087,0.,-.043649167310370842,0.,.043649167310370842,0.,
	    0.,-1.2745966692414834,.19999999999999998,0.,0.,-.2745966692414834,-.19999999999999998,0.,0.,1.5491933384829668,.18872983346207417,.11127016653792583,.87602816808281591,-.87602816808281591,-.11127016653792583,-.18872983346207417,-.023971831917184147,.023971831917184147,-1.0647580015448903,-.50983866769659336,.50983866769659336,1.0647580015448903,.13524199845510997,-.10983866769659336,.10983866769659336,-.13524199845510997,-.6196773353931867,.6196773353931867,0.,-.34364916731037087,0.,
	    .043649167310370842,0.,-.043649167310370842,0.,.34364916731037087,0.,-.19999999999999998,-.2745966692414834,0.,0.,.19999999999999998,-1.2745966692414834,0.,1.5491933384829668,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-.5,.5,0.,0.,.5,-.5,0.,0.,0.,0.,.043649167310370842,0.,-.34364916731037087,0.,.34364916731037087,0.,-.043649167310370842,0.,-.19999999999999998,1.2745966692414834,0.,0.,.19999999999999998,.2745966692414834,0.,-1.5491933384829668,0.,.11127016653792583,.18872983346207417,
	    .023971831917184147,-.023971831917184147,-.18872983346207417,-.11127016653792583,-.87602816808281591,.87602816808281591,-.13524199845510997,.10983866769659336,-.10983866769659336,.13524199845510997,1.0647580015448903,.50983866769659336,-.50983866769659336,-1.0647580015448903,.6196773353931867,-.6196773353931867,.043649167310370842,0.,-.043649167310370842,0.,.34364916731037087,0.,-.34364916731037087,0.,0.,.2745966692414834,.19999999999999998,0.,0.,1.2745966692414834,
	    -.19999999999999998,0.,0.,-1.5491933384829668,-.023971831917184147,-.023971831917184147,-.11127016653792583,.18872983346207417,.87602816808281591,.87602816808281591,.18872983346207417,-.11127016653792583,.13524199845510997,.10983866769659336,.50983866769659336,-1.0647580015448903,-1.0647580015448903,.50983866769659336,.10983866769659336,.13524199845510997,-.6196773353931867,-.6196773353931867 };
    static doublereal pi = 3.14159265358979;

    /* System generated locals */
    int32 i__1;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7;

    /* Local variables */
    static doublereal coef, e[9], f;
    static int32 i__;
    static doublereal poids, f1[9], f2[9], p1, p2, p3, p4, p5, p6, p7, p8, p9, dfm1dp[162]	/* was [2][9][9] */;
    extern /* Subroutine */ int dpaq2c_();
    static doublereal lambda, mu, poidel[9], dpr1, dpr2, dpr3, dpr4, dpr5, dpr6, dpr7, dpr8, dpr9, dpz1, dpz2, dpz3, dpz4, dpz5, dpz6, dpz7, dpz8, dpz9;

/* *************************************************************** */
/* BUT: CALCUL DE LA MATRICE DE RIGIDITE POUR L'ELEMENT QUAD AQ2C */
/* ---  QUADRANGLE Q2 AXISYMETRIQUE */
/* IN : */
/*      N           : Nb ondes circonferentielles */
/*      R(9) , Z(9) : Coordonnees des sommets de l'element courant */
/*      YOUNG et POISSON */
/* OUT: AE: MATRICE DE RIGIDITE DE L'ELEMENT. */

/*    [ AE(1) AE(2) AE(4) AE(7)  .......... ]  { Ur    (1) } */
/*    [       AE(3) AE(5)  :     .......... ]  { Utheta(1) } */
/*    [             AE(6)  :     .......... ]  { Uz    (1) } */
/*    [                    :     .......... ]  {   .....   } */
/*    [                          .......... ]  {   .....   } */
/*    [                          .......... ]  { Ur    (9) ] */
/*    [                             AE(377) ]  { Utheta(9) } */
/*    [                             AE(378) ]  { Uz    (9) } */
/*                                        27*27 */
/*    SEULE LA PARTIE SUPERIEURE DE HAUT EN BAS ET DE GAUCHE VERS LA */
/*    DROITE EST STOCKEE. Soit AE(27 * (27+1)/2) = AE(378) */

/*      4----- 7 -------3 */
/*      |               | */
/*      |               |  3 degrees de liberte' par noeud : */
/*      8      9        6  ( Ur(i), Utheta(i), Uz (i) ; i=1,9) */
/*      |               | */
/*      |               | */
/*      1------5--------2 */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* PROGRAMMEUR A. Hassim */
/*  .............................................................. */
/*     -- XYNPI: coordonnees pt. int. numeriques (element reference) */
    /* Parameter adjustments */
    --ae;
    --z__;
    --r__;

    /* Function Body */
/*     -- POIDSG: poids du schema d'integration numerique. */
/*     -- Valeurs des Polynomes de base aux pt. int. numerique. */
/*     -- Valeurs Derive'es des Poly. de base aux pt int. numerique. */


    lambda = *young * *poisson / ((*poisson + 1.) * (1. - *poisson * 2.));
    mu = *young / ((*poisson + 1.) * 2.);
    e[0] = lambda + mu * 2.;
    e[1] = lambda;
    e[2] = e[0];
    e[3] = lambda;
    e[4] = lambda;
    e[5] = e[0];
    e[6] = mu;
    e[7] = mu;
    e[8] = mu;

/*     -- Calcul de [DF]  ; [DF] * [DP] ; ... */

    dpaq2c_(&r__[1], &z__[1], &c__9, &c__9, poidsg, p, dp, f1, f2, dfm1dp, poidel);

/*     -- INITIALISATION DE AE: */

    for (i__ = 1; i__ <= 378; ++i__) {
	ae[i__] = 0.;
/* L1: */
    }

    if (*n == 0) {
	coef = pi * 2.;
    } else {
	coef = pi;
    }

    for (i__ = 1; i__ <= 9; ++i__) {
	poids = coef * poidel[i__ - 1];
	f = f1[i__ - 1];
	p1 = p[i__ * 9 - 9];
	p2 = p[i__ * 9 - 8];
	p3 = p[i__ * 9 - 7];
	p4 = p[i__ * 9 - 6];
	p5 = p[i__ * 9 - 5];
	p6 = p[i__ * 9 - 4];
	p7 = p[i__ * 9 - 3];
	p8 = p[i__ * 9 - 2];
	p9 = p[i__ * 9 - 1];

	dpr1 = dfm1dp[(i__ * 9 + 1 << 1) - 20];
	dpr2 = dfm1dp[(i__ * 9 + 2 << 1) - 20];
	dpr3 = dfm1dp[(i__ * 9 + 3 << 1) - 20];
	dpr4 = dfm1dp[(i__ * 9 + 4 << 1) - 20];
	dpr5 = dfm1dp[(i__ * 9 + 5 << 1) - 20];
	dpr6 = dfm1dp[(i__ * 9 + 6 << 1) - 20];
	dpr7 = dfm1dp[(i__ * 9 + 7 << 1) - 20];
	dpr8 = dfm1dp[(i__ * 9 + 8 << 1) - 20];
	dpr9 = dfm1dp[(i__ * 9 + 9 << 1) - 20];

	dpz1 = dfm1dp[(i__ * 9 + 1 << 1) - 19];
	dpz2 = dfm1dp[(i__ * 9 + 2 << 1) - 19];
	dpz3 = dfm1dp[(i__ * 9 + 3 << 1) - 19];
	dpz4 = dfm1dp[(i__ * 9 + 4 << 1) - 19];
	dpz5 = dfm1dp[(i__ * 9 + 5 << 1) - 19];
	dpz6 = dfm1dp[(i__ * 9 + 6 << 1) - 19];
	dpz7 = dfm1dp[(i__ * 9 + 7 << 1) - 19];
	dpz8 = dfm1dp[(i__ * 9 + 8 << 1) - 19];
	dpz9 = dfm1dp[(i__ * 9 + 9 << 1) - 19];
/* Computing 2nd power */
	d__1 = p1;
/* Computing 2nd power */
	d__2 = p1;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__3 = dpr1;
/* Computing 2nd power */
	d__4 = f;
/* Computing 2nd power */
	d__5 = dpz1;
/* Computing 2nd power */
	d__6 = f;
/* Computing 2nd power */
	d__7 = f;
	ae[1] += poids * (d__1 * d__1 * e[2] + d__2 * d__2 * (i__1 * i__1) * e[7] + p1 * 2. * e[1] * dpr1 * f + d__3 * d__3 * (d__4 * d__4) * e[0] + e[6] * (d__5 * d__5) * (d__6 * d__6)) / (d__7 * d__7);
/* Computing 2nd power */
	d__1 = f;
	ae[2] += poids * p1 * *n * (p1 * e[2] + p1 * e[7] + e[1] * dpr1 * f - e[7] * dpr1 * f) / (d__1 * d__1);
/* Computing 2nd power */
	d__1 = p1;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__2 = p1;
/* Computing 2nd power */
	d__3 = dpr1;
/* Computing 2nd power */
	d__4 = f;
/* Computing 2nd power */
	d__5 = dpz1;
/* Computing 2nd power */
	d__6 = f;
/* Computing 2nd power */
	d__7 = f;
	ae[3] += poids * (d__1 * d__1 * (i__1 * i__1) * e[2] + d__2 * d__2 * e[7] - p1 * 2. * e[7] * dpr1 * f + e[7] * (d__3 * d__3) * (d__4 * d__4) + e[8] * (d__5 * d__5) * (d__6 * d__6)) / (d__7 * d__7);
	ae[4] += poids * dpz1 * (e[6] * dpr1 * f + e[4] * p1 + e[3] * dpr1 * f) / f;
	ae[5] += poids * p1 * *n * dpz1 * (e[4] - e[8]) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = p1;
/* Computing 2nd power */
	d__2 = dpr1;
/* Computing 2nd power */
	d__3 = f;
/* Computing 2nd power */
	d__4 = dpz1;
/* Computing 2nd power */
	d__5 = f;
/* Computing 2nd power */
	d__6 = f;
	ae[6] += poids * (i__1 * i__1 * e[8] * (d__1 * d__1) + e[6] * (d__2 * d__2) * (d__3 * d__3) + e[5] * (d__4 * d__4) * (d__5 * d__5)) / (d__6 * d__6);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[7] += poids * (p1 * p2 * e[2] + p1 * p2 * (i__1 * i__1) * e[7] + p1 * dpr2 * e[1] * f + dpr1 * f * e[1] * p2 + dpr1 * (d__1 * d__1) * e[0] * dpr2 + dpz1 * e[6] * dpz2 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[8] -= poids * *n * (-p1 * p2 * e[2] - p1 * p2 * e[7] + p2 * e[7] * dpr1 * f - p1 * dpr2 * e[1] * f) / (d__1 * d__1);
	ae[9] += poids * (dpr1 * e[6] * dpz2 * f + dpz1 * e[4] * p2 + dpz1 * e[3] * dpr2 * f) / f;
/* Computing 2nd power */
	d__1 = p2;
/* Computing 2nd power */
	d__2 = p2;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__3 = dpr2;
/* Computing 2nd power */
	d__4 = f;
/* Computing 2nd power */
	d__5 = dpz2;
/* Computing 2nd power */
	d__6 = f;
/* Computing 2nd power */
	d__7 = f;
	ae[10] += poids * (d__1 * d__1 * e[2] + d__2 * d__2 * (i__1 * i__1) * e[7] + p2 * 2. * dpr2 * e[1] * f + d__3 * d__3 * (d__4 * d__4) * e[0] + e[6] * (d__5 * d__5) * (d__6 * d__6)) / (d__7 * d__7);
/* Computing 2nd power */
	d__1 = f;
	ae[11] += poids * *n * (p1 * p2 * e[2] + p1 * p2 * e[7] + dpr1 * f * e[1] * p2 - dpr2 * e[7] * p1 * f) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[12] += poids * (p2 * p1 * (i__1 * i__1) * e[2] + p1 * p2 * e[7] - p2 * e[7] * dpr1 * f - dpr2 * e[7] * p1 * f + e[7] * dpr2 * (d__1 * d__1) * dpr1 + dpz1 * e[8] * dpz2 * (d__2 * d__2)) / (d__3 * d__3);
	ae[13] += poids * *n * (-p1 * e[8] * dpz2 + dpz1 * e[4] * p2) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[14] += poids * p2 * *n * (p2 * e[2] + p2 * e[7] - e[7] * dpr2 * f + dpr2 * e[1] * f) / (d__1 * d__1);
/* Computing 2nd power */
	d__1 = p2;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__2 = p2;
/* Computing 2nd power */
	d__3 = dpr2;
/* Computing 2nd power */
	d__4 = f;
/* Computing 2nd power */
	d__5 = dpz2;
/* Computing 2nd power */
	d__6 = f;
/* Computing 2nd power */
	d__7 = f;
	ae[15] += poids * (d__1 * d__1 * (i__1 * i__1) * e[2] + d__2 * d__2 * e[7] - p2 * 2. * e[7] * dpr2 * f + e[7] * (d__3 * d__3) * (d__4 * d__4) + e[8] * (d__5 * d__5) * (d__6 * d__6)) / (d__7 * d__7);
	ae[16] += poids * (p1 * e[4] * dpz2 + dpr1 * e[3] * dpz2 * f + dpz1 * e[6] * dpr2 * f) / f;
	ae[17] -= poids * *n * (-p1 * e[4] * dpz2 + dpz1 * e[8] * p2) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[18] += poids * (p1 * (i__1 * i__1) * e[8] * p2 + dpr1 * e[6] * dpr2 * (d__1 * d__1) + dpz1 * e[5] * dpz2 * (d__2 * d__2)) / (d__3 * d__3);
	ae[19] += poids * dpz2 * (e[6] * dpr2 * f + e[4] * p2 + e[3] * dpr2 * f) / f;
	ae[20] += poids * p2 * *n * dpz2 * (e[4] - e[8]) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = p2;
/* Computing 2nd power */
	d__2 = dpr2;
/* Computing 2nd power */
	d__3 = f;
/* Computing 2nd power */
	d__4 = dpz2;
/* Computing 2nd power */
	d__5 = f;
/* Computing 2nd power */
	d__6 = f;
	ae[21] += poids * (i__1 * i__1 * e[8] * (d__1 * d__1) + e[6] * (d__2 * d__2) * (d__3 * d__3) + e[5] * (d__4 * d__4) * (d__5 * d__5)) / (d__6 * d__6);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[22] += poids * (p1 * p3 * e[2] + p3 * p1 * (i__1 * i__1) * e[7] + dpr1 * e[1] * p3 * f + dpr3 * f * e[1] * p1 + dpr3 * (d__1 * d__1) * e[0] * dpr1 + dpz1 * e[6] * dpz3 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[23] -= poids * *n * (-p1 * p3 * e[2] - p1 * p3 * e[7] + p3 * e[7] * dpr1 * f - dpr3 * f * e[1] * p1) / (d__1 * d__1);
	ae[24] += poids * (p3 * e[4] * dpz1 + dpr3 * e[3] * dpz1 * f + dpr1 * e[6] * dpz3 * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[25] += poids * (p2 * p3 * e[2] + p2 * p3 * (i__1 * i__1) * e[7] + p2 * e[1] * dpr3 * f + p3 * dpr2 * e[1] * f + dpr2 * (d__1 * d__1) * e[0] * dpr3 + dpz2 * e[6] * dpz3 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[26] += poids * *n * (p2 * p3 * e[2] + p2 * p3 * e[7] + p2 * e[1] * dpr3 * f - dpr2 * e[7] * p3 * f) / (d__1 * d__1);
	ae[27] += poids * (dpz2 * e[4] * p3 + dpr3 * e[3] * dpz2 * f + dpr2 * e[6] * dpz3 * f) / f;
/* Computing 2nd power */
	d__1 = p3;
/* Computing 2nd power */
	d__2 = p3;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__3 = dpr3;
/* Computing 2nd power */
	d__4 = f;
/* Computing 2nd power */
	d__5 = dpz3;
/* Computing 2nd power */
	d__6 = f;
/* Computing 2nd power */
	d__7 = f;
	ae[28] += poids * (d__1 * d__1 * e[2] + d__2 * d__2 * (i__1 * i__1) * e[7] + p3 * 2. * e[1] * dpr3 * f + d__3 * d__3 * (d__4 * d__4) * e[0] + e[6] * (d__5 * d__5) * (d__6 * d__6)) / (d__7 * d__7);
/* Computing 2nd power */
	d__1 = f;
	ae[29] += poids * *n * (p1 * p3 * e[2] + p1 * p3 * e[7] - p1 * e[7] * dpr3 * f + dpr1 * e[1] * p3 * f) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[30] += poids * (p3 * p1 * (i__1 * i__1) * e[2] + p1 * p3 * e[7] - p3 * e[7] * dpr1 * f - p1 * e[7] * dpr3 * f + e[7] * dpr3 * (d__1 * d__1) * dpr1 + dpz1 * e[8] * dpz3 * (d__2 * d__2)) / (d__3 * d__3);
	ae[31] += poids * *n * (-p1 * e[8] * dpz3 + p3 * e[4] * dpz1) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[32] -= poids * *n * (-p2 * p3 * e[2] - p2 * p3 * e[7] - p3 * dpr2 * e[1] * f + dpr3 * e[7] * p2 * f) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[33] += poids * (p2 * p3 * (i__1 * i__1) * e[2] + p2 * p3 * e[7] - dpr3 * e[7] * p2 * f - dpr2 * e[7] * p3 * f + e[7] * dpr2 * (d__1 * d__1) * dpr3 + dpz2 * e[8] * dpz3 * (d__2 * d__2)) / (d__3 * d__3);
	ae[34] -= poids * *n * (p2 * e[8] * dpz3 - dpz2 * e[4] * p3) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[35] += poids * p3 * *n * (p3 * e[2] + p3 * e[7] - e[7] * dpr3 * f + e[1] * dpr3 * f) / (d__1 * d__1);
/* Computing 2nd power */
	d__1 = p3;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__2 = p3;
/* Computing 2nd power */
	d__3 = dpr3;
/* Computing 2nd power */
	d__4 = f;
/* Computing 2nd power */
	d__5 = dpz3;
/* Computing 2nd power */
	d__6 = f;
/* Computing 2nd power */
	d__7 = f;
	ae[36] += poids * (d__1 * d__1 * (i__1 * i__1) * e[2] + d__2 * d__2 * e[7] - p3 * 2. * e[7] * dpr3 * f + e[7] * (d__3 * d__3) * (d__4 * d__4) + e[8] * (d__5 * d__5) * (d__6 * d__6)) / (d__7 * d__7);
	ae[37] += poids * (p1 * e[4] * dpz3 + dpr1 * e[3] * dpz3 * f + dpz1 * e[6] * dpr3 * f) / f;
	ae[38] -= poids * *n * (-p1 * e[4] * dpz3 + dpz1 * e[8] * p3) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[39] += poids * (p1 * (i__1 * i__1) * e[8] * p3 + dpr1 * e[6] * dpr3 * (d__1 * d__1) + dpz1 * e[5] * dpz3 * (d__2 * d__2)) / (d__3 * d__3);
	ae[40] += poids * (dpz2 * e[6] * dpr3 * f + dpz3 * e[4] * p2 + dpz3 * e[3] * dpr2 * f) / f;
	ae[41] += poids * *n * (dpz3 * e[4] * p2 - dpz2 * e[8] * p3) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[42] += poids * (p2 * (i__1 * i__1) * e[8] * p3 + dpr2 * e[6] * dpr3 * (d__1 * d__1) + dpz2 * e[5] * dpz3 * (d__2 * d__2)) / (d__3 * d__3);
	ae[43] += poids * dpz3 * (e[4] * p3 + e[3] * dpr3 * f + e[6] * dpr3 * f) / f;
	ae[44] += poids * p3 * *n * dpz3 * (e[4] - e[8]) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = p3;
/* Computing 2nd power */
	d__2 = dpr3;
/* Computing 2nd power */
	d__3 = f;
/* Computing 2nd power */
	d__4 = dpz3;
/* Computing 2nd power */
	d__5 = f;
/* Computing 2nd power */
	d__6 = f;
	ae[45] += poids * (i__1 * i__1 * e[8] * (d__1 * d__1) + e[6] * (d__2 * d__2) * (d__3 * d__3) + e[5] * (d__4 * d__4) * (d__5 * d__5)) / (d__6 * d__6);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[46] += poids * (p4 * p1 * e[2] + p4 * p1 * (i__1 * i__1) * e[7] + p4 * e[1] * dpr1 * f + p1 * e[1] * dpr4 * f + dpr4 * (d__1 * d__1) * e[0] * dpr1 + dpz1 * e[6] * dpz4 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[47] -= poids * *n * (-p4 * p1 * e[2] - p4 * p1 * e[7] - p1 * e[1] * dpr4 * f + p4 * e[7] * dpr1 * f) / (d__1 * d__1);
	ae[48] += poids * (p4 * e[4] * dpz1 + dpr4 * e[3] * dpz1 * f + dpr1 * e[6] * dpz4 * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[49] += poids * (p2 * p4 * e[2] + p4 * p2 * (i__1 * i__1) * e[7] + dpr2 * e[1] * p4 * f + dpr4 * e[1] * p2 * f + dpr4 * (d__1 * d__1) * e[0] * dpr2 + dpz2 * e[6] * dpz4 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[50] -= poids * *n * (-p2 * p4 * e[2] - p2 * p4 * e[7] + p4 * e[7] * dpr2 * f - dpr4 * e[1] * p2 * f) / (d__1 * d__1);
	ae[51] += poids * (dpr2 * e[6] * dpz4 * f + dpz2 * e[4] * p4 + dpz2 * e[3] * dpr4 * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[52] += poids * (p4 * p3 * e[2] + p4 * p3 * (i__1 * i__1) * e[7] + p4 * e[1] * dpr3 * f + dpr4 * f * e[1] * p3 + dpr4 * (d__1 * d__1) * e[0] * dpr3 + dpz3 * e[6] * dpz4 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[53] -= poids * *n * (-p4 * p3 * e[2] - p4 * p3 * e[7] - dpr4 * f * e[1] * p3 + dpr3 * e[7] * p4 * f) / (d__1 * d__1);
	ae[54] += poids * (dpr3 * e[6] * dpz4 * f + dpz3 * e[4] * p4 + dpz3 * e[3] * dpr4 * f) / f;
/* Computing 2nd power */
	d__1 = p4;
/* Computing 2nd power */
	d__2 = p4;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__3 = dpr4;
/* Computing 2nd power */
	d__4 = f;
/* Computing 2nd power */
	d__5 = dpz4;
/* Computing 2nd power */
	d__6 = f;
/* Computing 2nd power */
	d__7 = f;
	ae[55] += poids * (d__1 * d__1 * e[2] + d__2 * d__2 * (i__1 * i__1) * e[7] + p4 * 2. * e[1] * dpr4 * f + d__3 * d__3 * (d__4 * d__4) * e[0] + e[6] * (d__5 * d__5) * (d__6 * d__6)) / (d__7 * d__7);
/* Computing 2nd power */
	d__1 = f;
	ae[56] += poids * *n * (p4 * p1 * e[2] + p4 * p1 * e[7] + p4 * e[1] * dpr1 * f - e[7] * dpr4 * f * p1) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[57] += poids * (p4 * p1 * (i__1 * i__1) * e[2] + p4 * p1 * e[7] - p4 * e[7] * dpr1 * f - e[7] * dpr4 * f * p1 + e[7] * dpr4 * (d__1 * d__1) * dpr1 + dpz1 * e[8] * dpz4 * (d__2 * d__2)) / (d__3 * d__3);
	ae[58] += poids * *n * (-p1 * e[8] * dpz4 + p4 * e[4] * dpz1) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[59] += poids * *n * (p2 * p4 * e[2] + p2 * p4 * e[7] - p2 * e[7] * dpr4 * f + dpr2 * e[1] * p4 * f) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[60] += poids * (p4 * p2 * (i__1 * i__1) * e[2] + p2 * p4 * e[7] - p4 * e[7] * dpr2 * f - p2 * e[7] * dpr4 * f + e[7] * dpr4 * (d__1 * d__1) * dpr2 + dpz2 * e[8] * dpz4 * (d__2 * d__2)) / (d__3 * d__3);
	ae[61] += poids * *n * (-p2 * e[8] * dpz4 + dpz2 * e[4] * p4) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[62] += poids * *n * (p4 * p3 * e[2] + p4 * p3 * e[7] + p4 * e[1] * dpr3 * f - p3 * e[7] * dpr4 * f) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[63] += poids * (p3 * p4 * (i__1 * i__1) * e[2] + p4 * p3 * e[7] - p3 * e[7] * dpr4 * f - dpr3 * e[7] * p4 * f + e[7] * dpr3 * (d__1 * d__1) * dpr4 + dpz3 * e[8] * dpz4 * (d__2 * d__2)) / (d__3 * d__3);
	ae[64] += poids * *n * (-p3 * e[8] * dpz4 + dpz3 * e[4] * p4) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[65] += poids * p4 * *n * (p4 * e[2] + p4 * e[7] - e[7] * dpr4 * f + e[1] * dpr4 * f) / (d__1 * d__1);
/* Computing 2nd power */
	d__1 = p4;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__2 = p4;
/* Computing 2nd power */
	d__3 = dpr4;
/* Computing 2nd power */
	d__4 = f;
/* Computing 2nd power */
	d__5 = dpz4;
/* Computing 2nd power */
	d__6 = f;
/* Computing 2nd power */
	d__7 = f;
	ae[66] += poids * (d__1 * d__1 * (i__1 * i__1) * e[2] + d__2 * d__2 * e[7] - p4 * 2. * e[7] * dpr4 * f + e[7] * (d__3 * d__3) * (d__4 * d__4) + e[8] * (d__5 * d__5) * (d__6 * d__6)) / (d__7 * d__7);
	ae[67] += poids * (dpz1 * e[6] * dpr4 * f + dpz4 * e[4] * p1 + dpz4 * e[3] * dpr1 * f) / f;
	ae[68] -= poids * *n * (-dpz4 * e[4] * p1 + dpz1 * e[8] * p4) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[69] += poids * (p1 * (i__1 * i__1) * e[8] * p4 + dpr1 * e[6] * dpr4 * (d__1 * d__1) + dpz1 * e[5] * dpz4 * (d__2 * d__2)) / (d__3 * d__3);
	ae[70] += poids * (p2 * e[4] * dpz4 + dpr2 * e[3] * dpz4 * f + dpz2 * e[6] * dpr4 * f) / f;
	ae[71] -= poids * *n * (-p2 * e[4] * dpz4 + dpz2 * e[8] * p4) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[72] += poids * (p2 * (i__1 * i__1) * e[8] * p4 + dpr2 * e[6] * dpr4 * (d__1 * d__1) + dpz2 * e[5] * dpz4 * (d__2 * d__2)) / (d__3 * d__3);
	ae[73] += poids * (p3 * e[4] * dpz4 + dpr3 * e[3] * dpz4 * f + dpz3 * e[6] * dpr4 * f) / f;
	ae[74] -= poids * *n * (-p3 * e[4] * dpz4 + dpz3 * e[8] * p4) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[75] += poids * (p3 * (i__1 * i__1) * e[8] * p4 + dpr3 * e[6] * dpr4 * (d__1 * d__1) + dpz3 * e[5] * dpz4 * (d__2 * d__2)) / (d__3 * d__3);
	ae[76] += poids * dpz4 * (e[4] * p4 + e[3] * dpr4 * f + e[6] * dpr4 * f) / f;
	ae[77] += poids * p4 * *n * dpz4 * (e[4] - e[8]) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = p4;
/* Computing 2nd power */
	d__2 = dpr4;
/* Computing 2nd power */
	d__3 = f;
/* Computing 2nd power */
	d__4 = dpz4;
/* Computing 2nd power */
	d__5 = f;
/* Computing 2nd power */
	d__6 = f;
	ae[78] += poids * (i__1 * i__1 * e[8] * (d__1 * d__1) + e[6] * (d__2 * d__2) * (d__3 * d__3) + e[5] * (d__4 * d__4) * (d__5 * d__5)) / (d__6 * d__6);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[79] += poids * (p5 * p1 * e[2] + p1 * p5 * (i__1 * i__1) * e[7] + p1 * dpr5 * e[1] * f + dpr1 * f * e[1] * p5 + dpr1 * (d__1 * d__1) * e[0] * dpr5 + dpz1 * e[6] * dpz5 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[80] += poids * *n * (p5 * p1 * e[2] + p5 * p1 * e[7] - p5 * e[7] * dpr1 * f + p1 * dpr5 * e[1] * f) / (d__1 * d__1);
	ae[81] += poids * (dpr1 * e[6] * dpz5 * f + dpz1 * e[4] * p5 + dpz1 * dpr5 * e[3] * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[82] += poids * (p5 * p2 * e[2] + p5 * p2 * (i__1 * i__1) * e[7] + p5 * dpr2 * e[1] * f + dpr5 * f * e[1] * p2 + dpr5 * (d__1 * d__1) * e[0] * dpr2 + dpz2 * e[6] * dpz5 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[83] += poids * *n * (p5 * p2 * e[2] + p5 * p2 * e[7] - p5 * e[7] * dpr2 * f + dpr5 * f * e[1] * p2) / (d__1 * d__1);
	ae[84] += poids * (p5 * e[4] * dpz2 + dpr5 * e[3] * dpz2 * f + dpr2 * e[6] * dpz5 * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[85] += poids * (p3 * p5 * e[2] + p3 * p5 * (i__1 * i__1) * e[7] + p3 * dpr5 * e[1] * f + dpr3 * f * e[1] * p5 + dpr3 * (d__1 * d__1) * e[0] * dpr5 + dpz3 * e[6] * dpz5 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[86] += poids * *n * (p3 * p5 * e[2] + p3 * p5 * e[7] - p5 * e[7] * dpr3 * f + p3 * dpr5 * e[1] * f) / (d__1 * d__1);
	ae[87] += poids * (dpr3 * e[6] * dpz5 * f + dpz3 * e[4] * p5 + dpz3 * dpr5 * e[3] * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[88] += poids * (p4 * p5 * e[2] + p4 * p5 * (i__1 * i__1) * e[7] + p4 * dpr5 * e[1] * f + p5 * e[1] * dpr4 * f + dpr4 * (d__1 * d__1) * e[0] * dpr5 + dpz4 * e[6] * dpz5 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[89] += poids * *n * (p4 * p5 * e[2] + p4 * p5 * e[7] + p4 * dpr5 * e[1] * f - dpr4 * e[7] * p5 * f) / (d__1 * d__1);
	ae[90] += poids * (dpz4 * e[4] * p5 + dpr5 * e[3] * dpz4 * f + dpr4 * e[6] * dpz5 * f) / f;
/* Computing 2nd power */
	d__1 = p5;
/* Computing 2nd power */
	d__2 = p5;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__3 = dpr5;
/* Computing 2nd power */
	d__4 = f;
/* Computing 2nd power */
	d__5 = dpz5;
/* Computing 2nd power */
	d__6 = f;
/* Computing 2nd power */
	d__7 = f;
	ae[91] += poids * (d__1 * d__1 * e[2] + d__2 * d__2 * (i__1 * i__1) * e[7] + p5 * 2. * dpr5 * e[1] * f + d__3 * d__3 * (d__4 * d__4) * e[0] + e[6] * (d__5 * d__5) * (d__6 * d__6)) / (d__7 * d__7);
/* Computing 2nd power */
	d__1 = f;
	ae[92] += poids * *n * (p5 * p1 * e[2] + p5 * p1 * e[7] - dpr5 * e[7] * f * p1 + dpr1 * f * e[1] * p5) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[93] += poids * (p5 * p1 * (i__1 * i__1) * e[2] + p5 * p1 * e[7] - p5 * e[7] * dpr1 * f - dpr5 * e[7] * f * p1 + dpr5 * e[7] * (d__1 * d__1) * dpr1 + dpz1 * e[8] * dpz5 * (d__2 * d__2)) / (d__3 * d__3);
	ae[94] += poids * *n * (-p1 * e[8] * dpz5 + dpz1 * e[4] * p5) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[95] -= poids * *n * (-p5 * p2 * e[2] - p5 * p2 * e[7] + p2 * dpr5 * e[7] * f - p5 * dpr2 * e[1] * f) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[96] += poids * (p5 * p2 * (i__1 * i__1) * e[2] + p5 * p2 * e[7] - p5 * e[7] * dpr2 * f - p2 * dpr5 * e[7] * f + dpr5 * e[7] * (d__1 * d__1) * dpr2 + dpz2 * e[8] * dpz5 * (d__2 * d__2)) / (d__3 * d__3);
	ae[97] -= poids * *n * (p2 * e[8] * dpz5 - p5 * e[4] * dpz2) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[98] -= poids * *n * (-p3 * p5 * e[2] - p3 * p5 * e[7] + p3 * dpr5 * e[7] * f - dpr3 * f * e[1] * p5) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[99] += poids * (p3 * p5 * (i__1 * i__1) * e[2] + p3 * p5 * e[7] - p3 * dpr5 * e[7] * f - p5 * e[7] * dpr3 * f + e[7] * dpr3 * (d__1 * d__1) * dpr5 + dpz3 * e[8] * dpz5 * (d__2 * d__2)) / (d__3 * d__3);
	ae[100] -= poids * *n * (p3 * e[8] * dpz5 - dpz3 * e[4] * p5) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[101] -= poids * *n * (-p4 * p5 * e[2] - p4 * p5 * e[7] - p5 * e[1] * dpr4 * f + dpr5 * e[7] * p4 * f) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[102] += poids * (p4 * p5 * (i__1 * i__1) * e[2] + p4 * p5 * e[7] - dpr5 * e[7] * p4 * f - dpr4 * e[7] * p5 * f + e[7] * dpr4 * (d__1 * d__1) * dpr5 + dpz4 * e[8] * dpz5 * (d__2 * d__2)) / (d__3 * d__3);
	ae[103] -= poids * *n * (p4 * e[8] * dpz5 - dpz4 * e[4] * p5) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[104] += poids * p5 * *n * (p5 * e[2] + p5 * e[7] - dpr5 * e[7] * f + dpr5 * e[1] * f) / (d__1 * d__1);
/* Computing 2nd power */
	d__1 = p5;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__2 = p5;
/* Computing 2nd power */
	d__3 = dpr5;
/* Computing 2nd power */
	d__4 = f;
/* Computing 2nd power */
	d__5 = dpz5;
/* Computing 2nd power */
	d__6 = f;
/* Computing 2nd power */
	d__7 = f;
	ae[105] += poids * (d__1 * d__1 * (i__1 * i__1) * e[2] + d__2 * d__2 * e[7] - p5 * 2. * dpr5 * e[7] * f + d__3 * d__3 * e[7] * (d__4 * d__4) + e[8] * (d__5 * d__5) * (d__6 * d__6)) / (d__7 * d__7);
	ae[106] += poids * (dpz1 * e[6] * dpr5 * f + p1 * e[4] * dpz5 + dpz5 * e[3] * dpr1 * f) / f;
	ae[107] += poids * *n * (p1 * e[4] * dpz5 - dpz1 * e[8] * p5) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[108] += poids * (p1 * (i__1 * i__1) * e[8] * p5 + dpr1 * e[6] * dpr5 * (d__1 * d__1) + dpz1 * e[5] * dpz5 * (d__2 * d__2)) / (d__3 * d__3);
	ae[109] += poids * (dpz2 * e[6] * dpr5 * f + dpz5 * e[4] * p2 + dpz5 * e[3] * dpr2 * f) / f;
	ae[110] += poids * *n * (dpz5 * e[4] * p2 - dpz2 * e[8] * p5) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[111] += poids * (p2 * (i__1 * i__1) * e[8] * p5 + dpr2 * e[6] * dpr5 * (d__1 * d__1) + dpz2 * e[5] * dpz5 * (d__2 * d__2)) / (d__3 * d__3);
	ae[112] += poids * (p3 * e[4] * dpz5 + dpr3 * e[3] * dpz5 * f + dpz3 * e[6] * dpr5 * f) / f;
	ae[113] += poids * *n * (p3 * e[4] * dpz5 - dpz3 * e[8] * p5) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[114] += poids * (p3 * (i__1 * i__1) * e[8] * p5 + dpr3 * e[6] * dpr5 * (d__1 * d__1) + dpz3 * e[5] * dpz5 * (d__2 * d__2)) / (d__3 * d__3);
	ae[115] += poids * (dpz4 * e[6] * dpr5 * f + dpz5 * e[4] * p4 + dpz5 * e[3] * dpr4 * f) / f;
	ae[116] += poids * *n * (dpz5 * e[4] * p4 - dpz4 * e[8] * p5) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[117] += poids * (p4 * (i__1 * i__1) * e[8] * p5 + dpr4 * e[6] * dpr5 * (d__1 * d__1) + dpz4 * e[5] * dpz5 * (d__2 * d__2)) / (d__3 * d__3);
	ae[118] += poids * dpz5 * (p5 * e[4] + dpr5 * e[3] * f + e[6] * dpr5 * f) / f;
	ae[119] += poids * p5 * *n * dpz5 * (e[4] - e[8]) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = p5;
/* Computing 2nd power */
	d__2 = dpr5;
/* Computing 2nd power */
	d__3 = f;
/* Computing 2nd power */
	d__4 = dpz5;
/* Computing 2nd power */
	d__5 = f;
/* Computing 2nd power */
	d__6 = f;
	ae[120] += poids * (i__1 * i__1 * e[8] * (d__1 * d__1) + e[6] * (d__2 * d__2) * (d__3 * d__3) + e[5] * (d__4 * d__4) * (d__5 * d__5)) / (d__6 * d__6);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[121] += poids * (p6 * p1 * e[2] + p1 * p6 * (i__1 * i__1) * e[7] + dpr6 * e[1] * p1 * f + p6 * e[1] * dpr1 * f + dpr1 * (d__1 * d__1) * e[0] * dpr6 + dpz1 * e[6] * dpz6 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[122] += poids * *n * (p6 * p1 * e[2] + p6 * p1 * e[7] - p6 * e[7] * dpr1 * f + dpr6 * e[1] * p1 * f) / (d__1 * d__1);
	ae[123] += poids * (p6 * e[4] * dpz1 + dpr6 * e[3] * dpz1 * f + dpr1 * e[6] * dpz6 * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[124] += poids * (p2 * p6 * e[2] + p6 * p2 * (i__1 * i__1) * e[7] + p6 * dpr2 * e[1] * f + p2 * e[1] * dpr6 * f + dpr6 * (d__1 * d__1) * e[0] * dpr2 + dpz2 * e[6] * dpz6 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[125] += poids * *n * (p2 * p6 * e[2] + p2 * p6 * e[7] + p2 * e[1] * dpr6 * f - dpr2 * e[7] * p6 * f) / (d__1 * d__1);
	ae[126] += poids * (dpz2 * e[4] * p6 + dpr6 * e[3] * dpz2 * f + dpr2 * e[6] * dpz6 * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[127] += poids * (p3 * p6 * e[2] + p3 * p6 * (i__1 * i__1) * e[7] + p3 * e[1] * dpr6 * f + p6 * e[1] * dpr3 * f + dpr3 * (d__1 * d__1) * e[0] * dpr6 + dpz3 * e[6] * dpz6 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[128] += poids * *n * (p3 * p6 * e[2] + p3 * p6 * e[7] + p3 * e[1] * dpr6 * f - dpr3 * e[7] * p6 * f) / (d__1 * d__1);
	ae[129] += poids * (p6 * e[4] * dpz3 + dpr6 * e[3] * dpz3 * f + dpr3 * e[6] * dpz6 * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[130] += poids * (p6 * p4 * e[2] + p6 * p4 * (i__1 * i__1) * e[7] + p6 * e[1] * dpr4 * f + dpr6 * f * e[1] * p4 + dpr6 * (d__1 * d__1) * e[0] * dpr4 + dpz4 * e[6] * dpz6 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[131] += poids * *n * (p6 * p4 * e[2] + p6 * p4 * e[7] - p6 * e[7] * dpr4 * f + dpr6 * f * e[1] * p4) / (d__1 * d__1);
	ae[132] += poids * (p6 * e[4] * dpz4 + dpr6 * e[3] * dpz4 * f + dpr4 * e[6] * dpz6 * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[133] += poids * (p5 * p6 * e[2] + p5 * p6 * (i__1 * i__1) * e[7] + p5 * e[1] * dpr6 * f + dpr5 * f * e[1] * p6 + dpr5 * (d__1 * d__1) * e[0] * dpr6 + dpz5 * e[6] * dpz6 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[134] += poids * *n * (p5 * p6 * e[2] + p5 * p6 * e[7] - dpr5 * e[7] * f * p6 + p5 * e[1] * dpr6 * f) / (d__1 * d__1);
	ae[135] += poids * (dpz5 * e[4] * p6 + dpr6 * e[3] * dpz5 * f + dpr5 * e[6] * dpz6 * f) / f;
/* Computing 2nd power */
	d__1 = p6;
/* Computing 2nd power */
	d__2 = p6;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__3 = dpr6;
/* Computing 2nd power */
	d__4 = f;
/* Computing 2nd power */
	d__5 = dpz6;
/* Computing 2nd power */
	d__6 = f;
/* Computing 2nd power */
	d__7 = f;
	ae[136] += poids * (d__1 * d__1 * e[2] + d__2 * d__2 * (i__1 * i__1) * e[7] + p6 * 2. * e[1] * dpr6 * f + d__3 * d__3 * (d__4 * d__4) * e[0] + e[6] * (d__5 * d__5) * (d__6 * d__6)) / (d__7 * d__7);
/* Computing 2nd power */
	d__1 = f;
	ae[137] += poids * *n * (p6 * p1 * e[2] + p6 * p1 * e[7] + p6 * e[1] * dpr1 * f - p1 * e[7] * dpr6 * f) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[138] += poids * (p1 * p6 * (i__1 * i__1) * e[2] + p6 * p1 * e[7] - p1 * e[7] * dpr6 * f - p6 * e[7] * dpr1 * f + e[7] * dpr1 * (d__1 * d__1) * dpr6 + dpz1 * e[8] * dpz6 * (d__2 * d__2)) / (d__3 * d__3);
	ae[139] += poids * *n * (-p1 * e[8] * dpz6 + p6 * e[4] * dpz1) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[140] -= poids * *n * (-p2 * p6 * e[2] - p2 * p6 * e[7] + p2 * e[7] * dpr6 * f - p6 * dpr2 * e[1] * f) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[141] += poids * (p2 * p6 * (i__1 * i__1) * e[2] + p2 * p6 * e[7] - p2 * e[7] * dpr6 * f - dpr2 * e[7] * p6 * f + e[7] * dpr2 * (d__1 * d__1) * dpr6 + dpz2 * e[8] * dpz6 * (d__2 * d__2)) / (d__3 * d__3);
	ae[142] -= poids * *n * (p2 * e[8] * dpz6 - dpz2 * e[4] * p6) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[143] -= poids * *n * (-p3 * p6 * e[2] - p3 * p6 * e[7] - p6 * e[1] * dpr3 * f + dpr6 * e[7] * p3 * f) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[144] += poids * (p3 * p6 * (i__1 * i__1) * e[2] + p3 * p6 * e[7] - dpr6 * e[7] * p3 * f - dpr3 * e[7] * p6 * f + e[7] * dpr3 * (d__1 * d__1) * dpr6 + dpz3 * e[8] * dpz6 * (d__2 * d__2)) / (d__3 * d__3);
	ae[145] -= poids * *n * (p3 * e[8] * dpz6 - p6 * e[4] * dpz3) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[146] -= poids * *n * (-p6 * p4 * e[2] - p6 * p4 * e[7] - p6 * e[1] * dpr4 * f + p4 * e[7] * dpr6 * f) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[147] += poids * (p4 * p6 * (i__1 * i__1) * e[2] + p6 * p4 * e[7] - p4 * e[7] * dpr6 * f - p6 * e[7] * dpr4 * f + e[7] * dpr4 * (d__1 * d__1) * dpr6 + dpz4 * e[8] * dpz6 * (d__2 * d__2)) / (d__3 * d__3);
	ae[148] -= poids * *n * (p4 * e[8] * dpz6 - p6 * e[4] * dpz4) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[149] += poids * *n * (p5 * p6 * e[2] + p5 * p6 * e[7] - p5 * e[7] * dpr6 * f + dpr5 * f * e[1] * p6) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[150] += poids * (p5 * p6 * (i__1 * i__1) * e[2] + p5 * p6 * e[7] - p5 * e[7] * dpr6 * f - dpr5 * e[7] * f * p6 + dpr5 * e[7] * (d__1 * d__1) * dpr6 + dpz5 * e[8] * dpz6 * (d__2 * d__2)) / (d__3 * d__3);
	ae[151] += poids * *n * (-p5 * e[8] * dpz6 + dpz5 * e[4] * p6) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[152] += poids * p6 * *n * (p6 * e[2] + p6 * e[7] + e[1] * dpr6 * f - e[7] * dpr6 * f) / (d__1 * d__1);
/* Computing 2nd power */
	d__1 = p6;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__2 = p6;
/* Computing 2nd power */
	d__3 = dpr6;
/* Computing 2nd power */
	d__4 = f;
/* Computing 2nd power */
	d__5 = dpz6;
/* Computing 2nd power */
	d__6 = f;
/* Computing 2nd power */
	d__7 = f;
	ae[153] += poids * (d__1 * d__1 * (i__1 * i__1) * e[2] + d__2 * d__2 * e[7] - p6 * 2. * e[7] * dpr6 * f + e[7] * (d__3 * d__3) * (d__4 * d__4) + e[8] * (d__5 * d__5) * (d__6 * d__6)) / (d__7 * d__7);
	ae[154] += poids * (p1 * e[4] * dpz6 + dpr1 * e[3] * dpz6 * f + dpz1 * e[6] * dpr6 * f) / f;
	ae[155] += poids * *n * (p1 * e[4] * dpz6 - dpz1 * e[8] * p6) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[156] += poids * (p1 * (i__1 * i__1) * e[8] * p6 + dpr1 * e[6] * dpr6 * (d__1 * d__1) + dpz1 * e[5] * dpz6 * (d__2 * d__2)) / (d__3 * d__3);
	ae[157] += poids * (dpz2 * e[6] * dpr6 * f + dpz6 * e[4] * p2 + dpz6 * e[3] * dpr2 * f) / f;
	ae[158] += poids * *n * (dpz6 * e[4] * p2 - dpz2 * e[8] * p6) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[159] += poids * (p2 * (i__1 * i__1) * e[8] * p6 + dpr2 * e[6] * dpr6 * (d__1 * d__1) + dpz2 * e[5] * dpz6 * (d__2 * d__2)) / (d__3 * d__3);
	ae[160] += poids * (dpz3 * e[6] * dpr6 * f + p3 * e[4] * dpz6 + dpz6 * e[3] * dpr3 * f) / f;
	ae[161] += poids * *n * (p3 * e[4] * dpz6 - dpz3 * e[8] * p6) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[162] += poids * (p3 * (i__1 * i__1) * e[8] * p6 + dpr3 * e[6] * dpr6 * (d__1 * d__1) + dpz3 * e[5] * dpz6 * (d__2 * d__2)) / (d__3 * d__3);
	ae[163] += poids * (p4 * e[4] * dpz6 + dpr4 * e[3] * dpz6 * f + dpz4 * e[6] * dpr6 * f) / f;
	ae[164] += poids * *n * (p4 * e[4] * dpz6 - dpz4 * e[8] * p6) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[165] += poids * (p4 * (i__1 * i__1) * e[8] * p6 + dpr4 * e[6] * dpr6 * (d__1 * d__1) + dpz4 * e[5] * dpz6 * (d__2 * d__2)) / (d__3 * d__3);
	ae[166] += poids * (p5 * e[4] * dpz6 + dpr5 * e[3] * dpz6 * f + dpz5 * e[6] * dpr6 * f) / f;
	ae[167] += poids * *n * (p5 * e[4] * dpz6 - dpz5 * e[8] * p6) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[168] += poids * (p5 * (i__1 * i__1) * e[8] * p6 + dpr5 * e[6] * dpr6 * (d__1 * d__1) + dpz5 * e[5] * dpz6 * (d__2 * d__2)) / (d__3 * d__3);
	ae[169] += poids * dpz6 * (e[4] * p6 + e[3] * dpr6 * f + e[6] * dpr6 * f) / f;
	ae[170] += poids * p6 * *n * dpz6 * (e[4] - e[8]) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = p6;
/* Computing 2nd power */
	d__2 = dpr6;
/* Computing 2nd power */
	d__3 = f;
/* Computing 2nd power */
	d__4 = dpz6;
/* Computing 2nd power */
	d__5 = f;
/* Computing 2nd power */
	d__6 = f;
	ae[171] += poids * (i__1 * i__1 * e[8] * (d__1 * d__1) + e[6] * (d__2 * d__2) * (d__3 * d__3) + e[5] * (d__4 * d__4) * (d__5 * d__5)) / (d__6 * d__6);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[172] += poids * (p7 * p1 * e[2] + p1 * p7 * (i__1 * i__1) * e[7] + dpr7 * e[1] * p1 * f + p7 * e[1] * dpr1 * f + dpr1 * (d__1 * d__1) * e[0] * dpr7 + dpz1 * e[6] * dpz7 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[173] += poids * *n * (p7 * p1 * e[2] + p7 * p1 * e[7] - p7 * e[7] * dpr1 * f + dpr7 * e[1] * p1 * f) / (d__1 * d__1);
	ae[174] += poids * (dpr1 * e[6] * dpz7 * f + dpz1 * e[4] * p7 + dpz1 * e[3] * dpr7 * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[175] += poids * (p2 * p7 * e[2] + p7 * p2 * (i__1 * i__1) * e[7] + p7 * dpr2 * e[1] * f + p2 * e[1] * dpr7 * f + dpr7 * (d__1 * d__1) * e[0] * dpr2 + dpz2 * e[6] * dpz7 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[176] += poids * *n * (p2 * p7 * e[2] + p2 * p7 * e[7] + p2 * e[1] * dpr7 * f - dpr2 * e[7] * p7 * f) / (d__1 * d__1);
	ae[177] += poids * (dpr2 * e[6] * dpz7 * f + dpz2 * e[4] * p7 + dpz2 * e[3] * dpr7 * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[178] += poids * (p7 * p3 * e[2] + p3 * p7 * (i__1 * i__1) * e[7] + dpr7 * e[1] * p3 * f + dpr3 * f * e[1] * p7 + dpr3 * (d__1 * d__1) * e[0] * dpr7 + dpz3 * e[6] * dpz7 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[179] += poids * *n * (p7 * p3 * e[2] + p7 * p3 * e[7] - p7 * e[7] * dpr3 * f + dpr7 * e[1] * p3 * f) / (d__1 * d__1);
	ae[180] += poids * (p7 * e[4] * dpz3 + dpr7 * e[3] * dpz3 * f + dpr3 * e[6] * dpz7 * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[181] += poids * (p7 * p4 * e[2] + p7 * p4 * (i__1 * i__1) * e[7] + p7 * e[1] * dpr4 * f + dpr7 * f * e[1] * p4 + dpr7 * (d__1 * d__1) * e[0] * dpr4 + dpz4 * e[6] * dpz7 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[182] += poids * *n * (p7 * p4 * e[2] + p7 * p4 * e[7] - p7 * e[7] * dpr4 * f + dpr7 * f * e[1] * p4) / (d__1 * d__1);
	ae[183] += poids * (dpz4 * e[4] * p7 + dpr7 * e[3] * dpz4 * f + dpr4 * e[6] * dpz7 * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[184] += poids * (p7 * p5 * e[2] + p7 * p5 * (i__1 * i__1) * e[7] + p7 * dpr5 * e[1] * f + p5 * e[1] * dpr7 * f + dpr7 * (d__1 * d__1) * e[0] * dpr5 + dpz5 * e[6] * dpz7 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[185] += poids * *n * (p7 * p5 * e[2] + p7 * p5 * e[7] + p5 * e[1] * dpr7 * f - dpr5 * e[7] * p7 * f) / (d__1 * d__1);
	ae[186] += poids * (dpr5 * e[6] * dpz7 * f + dpz5 * e[4] * p7 + dpz5 * e[3] * dpr7 * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[187] += poids * (p6 * p7 * e[2] + p7 * p6 * (i__1 * i__1) * e[7] + dpr6 * e[1] * p7 * f + p6 * e[1] * dpr7 * f + dpr7 * (d__1 * d__1) * e[0] * dpr6 + dpz6 * e[6] * dpz7 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[188] += poids * *n * (p6 * p7 * e[2] + p6 * p7 * e[7] + p6 * e[1] * dpr7 * f - dpr6 * e[7] * p7 * f) / (d__1 * d__1);
	ae[189] += poids * (dpr6 * e[6] * dpz7 * f + dpz6 * e[4] * p7 + dpz6 * e[3] * dpr7 * f) / f;
/* Computing 2nd power */
	d__1 = p7;
/* Computing 2nd power */
	d__2 = p7;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__3 = dpr7;
/* Computing 2nd power */
	d__4 = f;
/* Computing 2nd power */
	d__5 = dpz7;
/* Computing 2nd power */
	d__6 = f;
/* Computing 2nd power */
	d__7 = f;
	ae[190] += poids * (d__1 * d__1 * e[2] + d__2 * d__2 * (i__1 * i__1) * e[7] + p7 * 2. * e[1] * dpr7 * f + d__3 * d__3 * (d__4 * d__4) * e[0] + e[6] * (d__5 * d__5) * (d__6 * d__6)) / (d__7 * d__7);
/* Computing 2nd power */
	d__1 = f;
	ae[191] += poids * *n * (p7 * p1 * e[2] + p7 * p1 * e[7] + p7 * e[1] * dpr1 * f - e[7] * dpr7 * f * p1) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[192] += poids * (p7 * p1 * (i__1 * i__1) * e[2] + p7 * p1 * e[7] - p7 * e[7] * dpr1 * f - e[7] * dpr7 * f * p1 + e[7] * dpr7 * (d__1 * d__1) * dpr1 + dpz1 * e[8] * dpz7 * (d__2 * d__2)) / (d__3 * d__3);
	ae[193] += poids * *n * (-p1 * e[8] * dpz7 + dpz1 * e[4] * p7) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[194] -= poids * *n * (-p2 * p7 * e[2] - p2 * p7 * e[7] - p7 * dpr2 * e[1] * f + p2 * e[7] * dpr7 * f) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[195] += poids * (p2 * p7 * (i__1 * i__1) * e[2] + p2 * p7 * e[7] - p2 * e[7] * dpr7 * f - dpr2 * e[7] * p7 * f + e[7] * dpr2 * (d__1 * d__1) * dpr7 + dpz2 * e[8] * dpz7 * (d__2 * d__2)) / (d__3 * d__3);
	ae[196] -= poids * *n * (p2 * e[8] * dpz7 - dpz2 * e[4] * p7) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[197] -= poids * *n * (-p7 * p3 * e[2] - p7 * p3 * e[7] + e[7] * dpr7 * f * p3 - dpr3 * f * e[1] * p7) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[198] += poids * (p7 * p3 * (i__1 * i__1) * e[2] + p7 * p3 * e[7] - p7 * e[7] * dpr3 * f - e[7] * dpr7 * f * p3 + e[7] * dpr7 * (d__1 * d__1) * dpr3 + dpz3 * e[8] * dpz7 * (d__2 * d__2)) / (d__3 * d__3);
	ae[199] -= poids * *n * (p3 * e[8] * dpz7 - p7 * e[4] * dpz3) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[200] -= poids * *n * (-p7 * p4 * e[2] - p7 * p4 * e[7] + e[7] * dpr7 * f * p4 - p7 * e[1] * dpr4 * f) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[201] += poids * (p7 * p4 * (i__1 * i__1) * e[2] + p7 * p4 * e[7] - p7 * e[7] * dpr4 * f - e[7] * dpr7 * f * p4 + e[7] * dpr7 * (d__1 * d__1) * dpr4 + dpz4 * e[8] * dpz7 * (d__2 * d__2)) / (d__3 * d__3);
	ae[202] -= poids * *n * (p4 * e[8] * dpz7 - dpz4 * e[4] * p7) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[203] += poids * *n * (p7 * p5 * e[2] + p7 * p5 * e[7] + p7 * dpr5 * e[1] * f - dpr7 * e[7] * p5 * f) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[204] += poids * (p7 * p5 * (i__1 * i__1) * e[2] + p7 * p5 * e[7] - dpr5 * e[7] * p7 * f - dpr7 * e[7] * p5 * f + e[7] * dpr7 * (d__1 * d__1) * dpr5 + dpz5 * e[8] * dpz7 * (d__2 * d__2)) / (d__3 * d__3);
	ae[205] += poids * *n * (-p5 * e[8] * dpz7 + dpz5 * e[4] * p7) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[206] += poids * *n * (p6 * p7 * e[2] + p6 * p7 * e[7] - p6 * e[7] * dpr7 * f + dpr6 * e[1] * p7 * f) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[207] += poids * (p7 * p6 * (i__1 * i__1) * e[2] + p6 * p7 * e[7] - dpr6 * e[7] * p7 * f - p6 * e[7] * dpr7 * f + e[7] * dpr7 * (d__1 * d__1) * dpr6 + dpz6 * e[8] * dpz7 * (d__2 * d__2)) / (d__3 * d__3);
	ae[208] += poids * *n * (-p6 * e[8] * dpz7 + dpz6 * e[4] * p7) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[209] += poids * p7 * *n * (p7 * e[2] + p7 * e[7] - e[7] * dpr7 * f + e[1] * dpr7 * f) / (d__1 * d__1);
/* Computing 2nd power */
	d__1 = p7;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__2 = p7;
/* Computing 2nd power */
	d__3 = dpr7;
/* Computing 2nd power */
	d__4 = f;
/* Computing 2nd power */
	d__5 = dpz7;
/* Computing 2nd power */
	d__6 = f;
/* Computing 2nd power */
	d__7 = f;
	ae[210] += poids * (d__1 * d__1 * (i__1 * i__1) * e[2] + d__2 * d__2 * e[7] - p7 * 2. * e[7] * dpr7 * f + e[7] * (d__3 * d__3) * (d__4 * d__4) + e[8] * (d__5 * d__5) * (d__6 * d__6)) / (d__7 * d__7);
	ae[211] += poids * (dpz1 * e[6] * dpr7 * f + dpz7 * e[4] * p1 + dpz7 * e[3] * dpr1 * f) / f;
	ae[212] += poids * *n * (dpz7 * e[4] * p1 - dpz1 * e[8] * p7) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[213] += poids * (p1 * (i__1 * i__1) * e[8] * p7 + dpr1 * e[6] * dpr7 * (d__1 * d__1) + dpz1 * e[5] * dpz7 * (d__2 * d__2)) / (d__3 * d__3);
	ae[214] += poids * (p2 * e[4] * dpz7 + dpr2 * e[3] * dpz7 * f + dpz2 * e[6] * dpr7 * f) / f;
	ae[215] += poids * *n * (p2 * e[4] * dpz7 - dpz2 * e[8] * p7) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[216] += poids * (p2 * (i__1 * i__1) * e[8] * p7 + dpr2 * e[6] * dpr7 * (d__1 * d__1) + dpz2 * e[5] * dpz7 * (d__2 * d__2)) / (d__3 * d__3);
	ae[217] += poids * (dpz3 * e[6] * dpr7 * f + dpz7 * e[4] * p3 + dpz7 * e[3] * dpr3 * f) / f;
	ae[218] += poids * *n * (dpz7 * e[4] * p3 - dpz3 * e[8] * p7) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[219] += poids * (p3 * (i__1 * i__1) * e[8] * p7 + dpr3 * e[6] * dpr7 * (d__1 * d__1) + dpz3 * e[5] * dpz7 * (d__2 * d__2)) / (d__3 * d__3);
	ae[220] += poids * (p4 * e[4] * dpz7 + dpr4 * e[3] * dpz7 * f + dpz4 * e[6] * dpr7 * f) / f;
	ae[221] += poids * *n * (p4 * e[4] * dpz7 - dpz4 * e[8] * p7) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[222] += poids * (p4 * (i__1 * i__1) * e[8] * p7 + dpr4 * e[6] * dpr7 * (d__1 * d__1) + dpz4 * e[5] * dpz7 * (d__2 * d__2)) / (d__3 * d__3);
	ae[223] += poids * (p5 * e[4] * dpz7 + dpr5 * e[3] * dpz7 * f + dpz5 * e[6] * dpr7 * f) / f;
	ae[224] += poids * *n * (p5 * e[4] * dpz7 - dpz5 * e[8] * p7) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[225] += poids * (p5 * (i__1 * i__1) * e[8] * p7 + dpr5 * e[6] * dpr7 * (d__1 * d__1) + dpz5 * e[5] * dpz7 * (d__2 * d__2)) / (d__3 * d__3);
	ae[226] += poids * (p6 * e[4] * dpz7 + dpr6 * e[3] * dpz7 * f + dpz6 * e[6] * dpr7 * f) / f;
	ae[227] += poids * *n * (p6 * e[4] * dpz7 - dpz6 * e[8] * p7) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[228] += poids * (p6 * (i__1 * i__1) * e[8] * p7 + dpr6 * e[6] * dpr7 * (d__1 * d__1) + dpz6 * e[5] * dpz7 * (d__2 * d__2)) / (d__3 * d__3);
	ae[229] += poids * dpz7 * (e[4] * p7 + e[3] * dpr7 * f + e[6] * dpr7 * f) / f;
	ae[230] += poids * p7 * *n * dpz7 * (e[4] - e[8]) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = p7;
/* Computing 2nd power */
	d__2 = dpr7;
/* Computing 2nd power */
	d__3 = f;
/* Computing 2nd power */
	d__4 = dpz7;
/* Computing 2nd power */
	d__5 = f;
/* Computing 2nd power */
	d__6 = f;
	ae[231] += poids * (i__1 * i__1 * e[8] * (d__1 * d__1) + e[6] * (d__2 * d__2) * (d__3 * d__3) + e[5] * (d__4 * d__4) * (d__5 * d__5)) / (d__6 * d__6);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[232] += poids * (p1 * p8 * e[2] + p8 * p1 * (i__1 * i__1) * e[7] + dpr1 * e[1] * p8 * f + p1 * e[1] * dpr8 * f + dpr8 * (d__1 * d__1) * e[0] * dpr1 + dpz1 * e[6] * dpz8 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[233] -= poids * *n * (-p1 * p8 * e[2] - p1 * p8 * e[7] - p1 * e[1] * dpr8 * f + p8 * e[7] * dpr1 * f) / (d__1 * d__1);
	ae[234] += poids * (p8 * e[4] * dpz1 + dpr8 * e[3] * dpz1 * f + dpr1 * e[6] * dpz8 * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[235] += poids * (p2 * p8 * e[2] + p8 * p2 * (i__1 * i__1) * e[7] + dpr2 * e[1] * p8 * f + dpr8 * f * e[1] * p2 + dpr8 * (d__1 * d__1) * e[0] * dpr2 + dpz2 * e[6] * dpz8 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[236] += poids * *n * (p2 * p8 * e[2] + p2 * p8 * e[7] - p8 * e[7] * dpr2 * f + dpr8 * f * e[1] * p2) / (d__1 * d__1);
	ae[237] += poids * (dpr2 * e[6] * dpz8 * f + dpz2 * e[4] * p8 + dpz2 * e[3] * dpr8 * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[238] += poids * (p8 * p3 * e[2] + p3 * p8 * (i__1 * i__1) * e[7] + p3 * e[1] * dpr8 * f + p8 * e[1] * dpr3 * f + dpr3 * (d__1 * d__1) * e[0] * dpr8 + dpz3 * e[6] * dpz8 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[239] += poids * *n * (p8 * p3 * e[2] + p8 * p3 * e[7] + p3 * e[1] * dpr8 * f - p8 * e[7] * dpr3 * f) / (d__1 * d__1);
	ae[240] += poids * (dpr3 * e[6] * dpz8 * f + dpz3 * e[4] * p8 + dpz3 * e[3] * dpr8 * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[241] += poids * (p8 * p4 * e[2] + p8 * p4 * (i__1 * i__1) * e[7] + p8 * e[1] * dpr4 * f + dpr8 * f * e[1] * p4 + dpr8 * (d__1 * d__1) * e[0] * dpr4 + dpz4 * e[6] * dpz8 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[242] += poids * *n * (p8 * p4 * e[2] + p8 * p4 * e[7] + dpr8 * f * e[1] * p4 - dpr4 * e[7] * p8 * f) / (d__1 * d__1);
	ae[243] += poids * (dpr4 * e[6] * dpz8 * f + dpz4 * e[4] * p8 + dpz4 * e[3] * dpr8 * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[244] += poids * (p8 * p5 * e[2] + p8 * p5 * (i__1 * i__1) * e[7] + p8 * dpr5 * e[1] * f + dpr8 * f * e[1] * p5 + dpr8 * (d__1 * d__1) * e[0] * dpr5 + dpz5 * e[6] * dpz8 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[245] -= poids * *n * (-p8 * p5 * e[2] - p8 * p5 * e[7] - dpr8 * f * e[1] * p5 + p8 * dpr5 * e[7] * f) / (d__1 * d__1);
	ae[246] += poids * (dpr5 * e[6] * dpz8 * f + dpz5 * e[4] * p8 + dpz5 * e[3] * dpr8 * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[247] += poids * (p6 * p8 * e[2] + p6 * p8 * (i__1 * i__1) * e[7] + p6 * e[1] * dpr8 * f + dpr6 * f * e[1] * p8 + dpr6 * (d__1 * d__1) * e[0] * dpr8 + dpz6 * e[6] * dpz8 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[248] -= poids * *n * (-p6 * p8 * e[2] - p6 * p8 * e[7] + p8 * e[7] * dpr6 * f - p6 * e[1] * dpr8 * f) / (d__1 * d__1);
	ae[249] += poids * (dpr6 * e[6] * dpz8 * f + dpz6 * e[4] * p8 + dpz6 * e[3] * dpr8 * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[250] += poids * (p8 * p7 * e[2] + p8 * p7 * (i__1 * i__1) * e[7] + p8 * e[1] * dpr7 * f + dpr8 * f * e[1] * p7 + dpr8 * (d__1 * d__1) * e[0] * dpr7 + dpz7 * e[6] * dpz8 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[251] -= poids * *n * (-p8 * p7 * e[2] - p8 * p7 * e[7] - dpr8 * f * e[1] * p7 + dpr7 * e[7] * p8 * f) / (d__1 * d__1);
	ae[252] += poids * (p8 * e[4] * dpz7 + dpr8 * e[3] * dpz7 * f + dpr7 * e[6] * dpz8 * f) / f;
/* Computing 2nd power */
	d__1 = p8;
/* Computing 2nd power */
	d__2 = p8;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__3 = dpr8;
/* Computing 2nd power */
	d__4 = f;
/* Computing 2nd power */
	d__5 = dpz8;
/* Computing 2nd power */
	d__6 = f;
/* Computing 2nd power */
	d__7 = f;
	ae[253] += poids * (d__1 * d__1 * e[2] + d__2 * d__2 * (i__1 * i__1) * e[7] + p8 * 2. * e[1] * dpr8 * f + d__3 * d__3 * (d__4 * d__4) * e[0] + e[6] * (d__5 * d__5) * (d__6 * d__6)) / (d__7 * d__7);
/* Computing 2nd power */
	d__1 = f;
	ae[254] += poids * *n * (p1 * p8 * e[2] + p1 * p8 * e[7] - p1 * dpr8 * e[7] * f + dpr1 * e[1] * p8 * f) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[255] += poids * (p8 * p1 * (i__1 * i__1) * e[2] + p1 * p8 * e[7] - p8 * e[7] * dpr1 * f - p1 * dpr8 * e[7] * f + dpr8 * e[7] * (d__1 * d__1) * dpr1 + dpz1 * e[8] * dpz8 * (d__2 * d__2)) / (d__3 * d__3);
	ae[256] += poids * *n * (-p1 * e[8] * dpz8 + p8 * e[4] * dpz1) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[257] -= poids * *n * (-p2 * p8 * e[2] - p2 * p8 * e[7] + p2 * dpr8 * e[7] * f - dpr2 * e[1] * p8 * f) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[258] += poids * (p2 * p8 * (i__1 * i__1) * e[2] + p2 * p8 * e[7] - p2 * dpr8 * e[7] * f - p8 * e[7] * dpr2 * f + e[7] * dpr2 * (d__1 * d__1) * dpr8 + dpz2 * e[8] * dpz8 * (d__2 * d__2)) / (d__3 * d__3);
	ae[259] -= poids * *n * (p2 * e[8] * dpz8 - dpz2 * e[4] * p8) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[260] -= poids * *n * (-p8 * p3 * e[2] - p8 * p3 * e[7] - p8 * e[1] * dpr3 * f + dpr8 * e[7] * p3 * f) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[261] += poids * (p8 * p3 * (i__1 * i__1) * e[2] + p8 * p3 * e[7] - p8 * e[7] * dpr3 * f - dpr8 * e[7] * p3 * f + dpr8 * e[7] * (d__1 * d__1) * dpr3 + dpz3 * e[8] * dpz8 * (d__2 * d__2)) / (d__3 * d__3);
	ae[262] -= poids * *n * (p3 * e[8] * dpz8 - dpz3 * e[4] * p8) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[263] -= poids * *n * (-p8 * p4 * e[2] - p8 * p4 * e[7] + p4 * dpr8 * e[7] * f - p8 * e[1] * dpr4 * f) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[264] += poids * (p4 * p8 * (i__1 * i__1) * e[2] + p8 * p4 * e[7] - p4 * dpr8 * e[7] * f - dpr4 * e[7] * p8 * f + e[7] * dpr4 * (d__1 * d__1) * dpr8 + dpz4 * e[8] * dpz8 * (d__2 * d__2)) / (d__3 * d__3);
	ae[265] -= poids * *n * (p4 * e[8] * dpz8 - dpz4 * e[4] * p8) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[266] += poids * *n * (p8 * p5 * e[2] + p8 * p5 * e[7] + p8 * dpr5 * e[1] * f - dpr8 * e[7] * p5 * f) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[267] += poids * (p8 * p5 * (i__1 * i__1) * e[2] + p8 * p5 * e[7] - p8 * dpr5 * e[7] * f - dpr8 * e[7] * p5 * f + dpr8 * e[7] * (d__1 * d__1) * dpr5 + dpz5 * e[8] * dpz8 * (d__2 * d__2)) / (d__3 * d__3);
	ae[268] += poids * *n * (-p5 * e[8] * dpz8 + dpz5 * e[4] * p8) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[269] += poids * *n * (p6 * p8 * e[2] + p6 * p8 * e[7] - p6 * dpr8 * e[7] * f + dpr6 * f * e[1] * p8) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[270] += poids * (p8 * p6 * (i__1 * i__1) * e[2] + p6 * p8 * e[7] - p8 * e[7] * dpr6 * f - p6 * dpr8 * e[7] * f + dpr8 * e[7] * (d__1 * d__1) * dpr6 + dpz6 * e[8] * dpz8 * (d__2 * d__2)) / (d__3 * d__3);
	ae[271] += poids * *n * (-p6 * e[8] * dpz8 + dpz6 * e[4] * p8) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[272] += poids * *n * (p8 * p7 * e[2] + p8 * p7 * e[7] - p7 * dpr8 * e[7] * f + p8 * e[1] * dpr7 * f) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[273] += poids * (p7 * p8 * (i__1 * i__1) * e[2] + p8 * p7 * e[7] - p7 * dpr8 * e[7] * f - dpr7 * e[7] * p8 * f + e[7] * dpr7 * (d__1 * d__1) * dpr8 + dpz7 * e[8] * dpz8 * (d__2 * d__2)) / (d__3 * d__3);
	ae[274] += poids * *n * (-p7 * e[8] * dpz8 + p8 * e[4] * dpz7) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[275] += poids * p8 * *n * (p8 * e[2] + p8 * e[7] + e[1] * dpr8 * f - dpr8 * e[7] * f) / (d__1 * d__1);
/* Computing 2nd power */
	d__1 = p8;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__2 = p8;
/* Computing 2nd power */
	d__3 = dpr8;
/* Computing 2nd power */
	d__4 = f;
/* Computing 2nd power */
	d__5 = dpz8;
/* Computing 2nd power */
	d__6 = f;
/* Computing 2nd power */
	d__7 = f;
	ae[276] += poids * (d__1 * d__1 * (i__1 * i__1) * e[2] + d__2 * d__2 * e[7] - p8 * 2. * dpr8 * e[7] * f + d__3 * d__3 * e[7] * (d__4 * d__4) + e[8] * (d__5 * d__5) * (d__6 * d__6)) / (d__7 * d__7);
	ae[277] += poids * (p1 * e[4] * dpz8 + dpr1 * e[3] * dpz8 * f + dpz1 * e[6] * dpr8 * f) / f;
	ae[278] -= poids * *n * (-p1 * e[4] * dpz8 + dpz1 * e[8] * p8) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[279] += poids * (p1 * (i__1 * i__1) * e[8] * p8 + dpr1 * e[6] * dpr8 * (d__1 * d__1) + dpz1 * e[5] * dpz8 * (d__2 * d__2)) / (d__3 * d__3);
	ae[280] += poids * (p2 * e[4] * dpz8 + dpr2 * e[3] * dpz8 * f + dpz2 * e[6] * dpr8 * f) / f;
	ae[281] += poids * *n * (p2 * e[4] * dpz8 - dpz2 * e[8] * p8) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[282] += poids * (p2 * (i__1 * i__1) * e[8] * p8 + dpr2 * e[6] * dpr8 * (d__1 * d__1) + dpz2 * e[5] * dpz8 * (d__2 * d__2)) / (d__3 * d__3);
	ae[283] += poids * (p3 * e[4] * dpz8 + dpr3 * e[3] * dpz8 * f + dpz3 * e[6] * dpr8 * f) / f;
	ae[284] += poids * *n * (p3 * e[4] * dpz8 - dpz3 * e[8] * p8) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[285] += poids * (p3 * (i__1 * i__1) * e[8] * p8 + dpr3 * e[6] * dpr8 * (d__1 * d__1) + dpz3 * e[5] * dpz8 * (d__2 * d__2)) / (d__3 * d__3);
	ae[286] += poids * (dpz4 * e[6] * dpr8 * f + p4 * e[4] * dpz8 + dpz8 * e[3] * dpr4 * f) / f;
	ae[287] += poids * *n * (p4 * e[4] * dpz8 - dpz4 * e[8] * p8) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[288] += poids * (p4 * (i__1 * i__1) * e[8] * p8 + dpr4 * e[6] * dpr8 * (d__1 * d__1) + dpz4 * e[5] * dpz8 * (d__2 * d__2)) / (d__3 * d__3);
	ae[289] += poids * (p5 * e[4] * dpz8 + dpr5 * e[3] * dpz8 * f + dpz5 * e[6] * dpr8 * f) / f;
	ae[290] -= poids * *n * (-p5 * e[4] * dpz8 + dpz5 * e[8] * p8) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[291] += poids * (p5 * (i__1 * i__1) * e[8] * p8 + dpr5 * e[6] * dpr8 * (d__1 * d__1) + dpz5 * e[5] * dpz8 * (d__2 * d__2)) / (d__3 * d__3);
	ae[292] += poids * (dpz6 * e[6] * dpr8 * f + p6 * e[4] * dpz8 + dpz8 * e[3] * dpr6 * f) / f;
	ae[293] -= poids * *n * (-p6 * e[4] * dpz8 + dpz6 * e[8] * p8) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[294] += poids * (p6 * (i__1 * i__1) * e[8] * p8 + dpr6 * e[6] * dpr8 * (d__1 * d__1) + dpz6 * e[5] * dpz8 * (d__2 * d__2)) / (d__3 * d__3);
	ae[295] += poids * (p7 * e[4] * dpz8 + dpr7 * e[3] * dpz8 * f + dpz7 * e[6] * dpr8 * f) / f;
	ae[296] -= poids * *n * (-p7 * e[4] * dpz8 + dpz7 * e[8] * p8) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[297] += poids * (p7 * (i__1 * i__1) * e[8] * p8 + dpr7 * e[6] * dpr8 * (d__1 * d__1) + dpz7 * e[5] * dpz8 * (d__2 * d__2)) / (d__3 * d__3);
	ae[298] += poids * dpz8 * (e[6] * dpr8 * f + e[4] * p8 + e[3] * dpr8 * f) / f;
	ae[299] += poids * p8 * *n * dpz8 * (e[4] - e[8]) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = p8;
/* Computing 2nd power */
	d__2 = dpr8;
/* Computing 2nd power */
	d__3 = f;
/* Computing 2nd power */
	d__4 = dpz8;
/* Computing 2nd power */
	d__5 = f;
/* Computing 2nd power */
	d__6 = f;
	ae[300] += poids * (i__1 * i__1 * e[8] * (d__1 * d__1) + e[6] * (d__2 * d__2) * (d__3 * d__3) + e[5] * (d__4 * d__4) * (d__5 * d__5)) / (d__6 * d__6);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[301] += poids * (p9 * p1 * e[2] + p1 * p9 * (i__1 * i__1) * e[7] + dpr9 * e[1] * p1 * f + p9 * e[1] * dpr1 * f + dpr1 * (d__1 * d__1) * e[0] * dpr9 + dpz1 * e[6] * dpz9 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[302] += poids * *n * (p9 * p1 * e[2] + p9 * p1 * e[7] - p9 * e[7] * dpr1 * f + dpr9 * e[1] * p1 * f) / (d__1 * d__1);
	ae[303] += poids * (p9 * e[4] * dpz1 + dpr9 * e[3] * dpz1 * f + dpr1 * e[6] * dpz9 * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[304] += poids * (p2 * p9 * e[2] + p2 * p9 * (i__1 * i__1) * e[7] + p2 * e[1] * dpr9 * f + dpr2 * f * e[1] * p9 + dpr2 * (d__1 * d__1) * e[0] * dpr9 + dpz2 * e[6] * dpz9 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[305] += poids * *n * (p2 * p9 * e[2] + p2 * p9 * e[7] - e[7] * dpr2 * f * p9 + p2 * e[1] * dpr9 * f) / (d__1 * d__1);
	ae[306] += poids * (dpr2 * e[6] * dpz9 * f + dpz2 * e[4] * p9 + dpz2 * e[3] * dpr9 * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[307] += poids * (p3 * p9 * e[2] + p3 * p9 * (i__1 * i__1) * e[7] + p3 * e[1] * dpr9 * f + p9 * e[1] * dpr3 * f + dpr3 * (d__1 * d__1) * e[0] * dpr9 + dpz3 * e[6] * dpz9 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[308] += poids * *n * (p3 * p9 * e[2] + p3 * p9 * e[7] + p3 * e[1] * dpr9 * f - dpr3 * e[7] * p9 * f) / (d__1 * d__1);
	ae[309] += poids * (p9 * e[4] * dpz3 + dpr9 * e[3] * dpz3 * f + dpr3 * e[6] * dpz9 * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[310] += poids * (p9 * p4 * e[2] + p9 * p4 * (i__1 * i__1) * e[7] + p9 * e[1] * dpr4 * f + p4 * e[1] * dpr9 * f + dpr9 * (d__1 * d__1) * e[0] * dpr4 + dpz4 * e[6] * dpz9 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[311] += poids * *n * (p9 * p4 * e[2] + p9 * p4 * e[7] + p4 * e[1] * dpr9 * f - p9 * e[7] * dpr4 * f) / (d__1 * d__1);
	ae[312] += poids * (dpr4 * e[6] * dpz9 * f + dpz4 * e[4] * p9 + dpz4 * e[3] * dpr9 * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[313] += poids * (p9 * p5 * e[2] + p5 * p9 * (i__1 * i__1) * e[7] + dpr9 * e[1] * p5 * f + dpr5 * f * e[1] * p9 + dpr5 * (d__1 * d__1) * e[0] * dpr9 + dpz5 * e[6] * dpz9 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[314] += poids * *n * (p9 * p5 * e[2] + p9 * p5 * e[7] - p9 * dpr5 * e[7] * f + dpr9 * e[1] * p5 * f) / (d__1 * d__1);
	ae[315] += poids * (p9 * e[4] * dpz5 + dpr9 * e[3] * dpz5 * f + dpr5 * e[6] * dpz9 * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[316] += poids * (p9 * p6 * e[2] + p9 * p6 * (i__1 * i__1) * e[7] + p9 * e[1] * dpr6 * f + p6 * e[1] * dpr9 * f + dpr9 * (d__1 * d__1) * e[0] * dpr6 + dpz6 * e[6] * dpz9 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[317] += poids * *n * (p9 * p6 * e[2] + p9 * p6 * e[7] + p6 * e[1] * dpr9 * f - p9 * e[7] * dpr6 * f) / (d__1 * d__1);
	ae[318] += poids * (dpz6 * e[4] * p9 + dpr9 * e[3] * dpz6 * f + dpr6 * e[6] * dpz9 * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[319] += poids * (p9 * p7 * e[2] + p9 * p7 * (i__1 * i__1) * e[7] + p9 * e[1] * dpr7 * f + p7 * e[1] * dpr9 * f + dpr9 * (d__1 * d__1) * e[0] * dpr7 + dpz7 * e[6] * dpz9 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[320] += poids * *n * (p9 * p7 * e[2] + p9 * p7 * e[7] + p7 * e[1] * dpr9 * f - p9 * e[7] * dpr7 * f) / (d__1 * d__1);
	ae[321] += poids * (dpr7 * e[6] * dpz9 * f + dpz7 * e[4] * p9 + dpz7 * e[3] * dpr9 * f) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[322] += poids * (p8 * p9 * e[2] + p8 * p9 * (i__1 * i__1) * e[7] + p8 * e[1] * dpr9 * f + p9 * e[1] * dpr8 * f + dpr8 * (d__1 * d__1) * e[0] * dpr9 + dpz8 * e[6] * dpz9 * (d__2 * d__2)) / (d__3 * d__3);
/* Computing 2nd power */
	d__1 = f;
	ae[323] += poids * *n * (p8 * p9 * e[2] + p8 * p9 * e[7] + p8 * e[1] * dpr9 * f - dpr8 * e[7] * p9 * f) / (d__1 * d__1);
	ae[324] += poids * (p9 * e[4] * dpz8 + dpr9 * e[3] * dpz8 * f + dpr8 * e[6] * dpz9 * f) / f;
/* Computing 2nd power */
	d__1 = p9;
/* Computing 2nd power */
	d__2 = p9;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__3 = dpr9;
/* Computing 2nd power */
	d__4 = f;
/* Computing 2nd power */
	d__5 = dpz9;
/* Computing 2nd power */
	d__6 = f;
/* Computing 2nd power */
	d__7 = f;
	ae[325] += poids * (d__1 * d__1 * e[2] + d__2 * d__2 * (i__1 * i__1) * e[7] + p9 * 2. * e[1] * dpr9 * f + d__3 * d__3 * (d__4 * d__4) * e[0] + e[6] * (d__5 * d__5) * (d__6 * d__6)) / (d__7 * d__7);
/* Computing 2nd power */
	d__1 = f;
	ae[326] += poids * *n * (p9 * p1 * e[2] + p9 * p1 * e[7] + p9 * e[1] * dpr1 * f - p1 * e[7] * dpr9 * f) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[327] += poids * (p1 * p9 * (i__1 * i__1) * e[2] + p9 * p1 * e[7] - p1 * e[7] * dpr9 * f - p9 * e[7] * dpr1 * f + e[7] * dpr1 * (d__1 * d__1) * dpr9 + dpz1 * e[8] * dpz9 * (d__2 * d__2)) / (d__3 * d__3);
	ae[328] += poids * *n * (-p1 * e[8] * dpz9 + p9 * e[4] * dpz1) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[329] -= poids * *n * (-p2 * p9 * e[2] - p2 * p9 * e[7] + p2 * e[7] * dpr9 * f - dpr2 * f * e[1] * p9) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[330] += poids * (p2 * p9 * (i__1 * i__1) * e[2] + p2 * p9 * e[7] - p2 * e[7] * dpr9 * f - e[7] * dpr2 * f * p9 + e[7] * dpr2 * (d__1 * d__1) * dpr9 + dpz2 * e[8] * dpz9 * (d__2 * d__2)) / (d__3 * d__3);
	ae[331] -= poids * *n * (p2 * e[8] * dpz9 - dpz2 * e[4] * p9) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[332] -= poids * *n * (-p3 * p9 * e[2] - p3 * p9 * e[7] - p9 * e[1] * dpr3 * f + p3 * e[7] * dpr9 * f) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[333] += poids * (p3 * p9 * (i__1 * i__1) * e[2] + p3 * p9 * e[7] - p3 * e[7] * dpr9 * f - dpr3 * e[7] * p9 * f + e[7] * dpr3 * (d__1 * d__1) * dpr9 + dpz3 * e[8] * dpz9 * (d__2 * d__2)) / (d__3 * d__3);
	ae[334] -= poids * *n * (p3 * e[8] * dpz9 - p9 * e[4] * dpz3) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[335] -= poids * *n * (-p9 * p4 * e[2] - p9 * p4 * e[7] - p9 * e[1] * dpr4 * f + dpr9 * e[7] * p4 * f) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[336] += poids * (p9 * p4 * (i__1 * i__1) * e[2] + p9 * p4 * e[7] - p9 * e[7] * dpr4 * f - dpr9 * e[7] * p4 * f + e[7] * dpr9 * (d__1 * d__1) * dpr4 + dpz4 * e[8] * dpz9 * (d__2 * d__2)) / (d__3 * d__3);
	ae[337] -= poids * *n * (p4 * e[8] * dpz9 - dpz4 * e[4] * p9) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[338] += poids * *n * (p9 * p5 * e[2] + p9 * p5 * e[7] - p5 * e[7] * dpr9 * f + dpr5 * f * e[1] * p9) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[339] += poids * (p9 * p5 * (i__1 * i__1) * e[2] + p9 * p5 * e[7] - p9 * dpr5 * e[7] * f - p5 * e[7] * dpr9 * f + e[7] * dpr9 * (d__1 * d__1) * dpr5 + dpz5 * e[8] * dpz9 * (d__2 * d__2)) / (d__3 * d__3);
	ae[340] += poids * *n * (-p5 * e[8] * dpz9 + p9 * e[4] * dpz5) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[341] += poids * *n * (p9 * p6 * e[2] + p9 * p6 * e[7] - e[7] * dpr9 * f * p6 + p9 * e[1] * dpr6 * f) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[342] += poids * (p9 * p6 * (i__1 * i__1) * e[2] + p9 * p6 * e[7] - p9 * e[7] * dpr6 * f - e[7] * dpr9 * f * p6 + e[7] * dpr9 * (d__1 * d__1) * dpr6 + dpz6 * e[8] * dpz9 * (d__2 * d__2)) / (d__3 * d__3);
	ae[343] += poids * *n * (-p6 * e[8] * dpz9 + dpz6 * e[4] * p9) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[344] += poids * *n * (p9 * p7 * e[2] + p9 * p7 * e[7] + p9 * e[1] * dpr7 * f - dpr9 * e[7] * p7 * f) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[345] += poids * (p9 * p7 * (i__1 * i__1) * e[2] + p9 * p7 * e[7] - p9 * e[7] * dpr7 * f - dpr9 * e[7] * p7 * f + e[7] * dpr9 * (d__1 * d__1) * dpr7 + dpz7 * e[8] * dpz9 * (d__2 * d__2)) / (d__3 * d__3);
	ae[346] += poids * *n * (-p7 * e[8] * dpz9 + dpz7 * e[4] * p9) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[347] -= poids * *n * (-p8 * p9 * e[2] - p8 * p9 * e[7] - p9 * e[1] * dpr8 * f + p8 * e[7] * dpr9 * f) / (d__1 * d__1);
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[348] += poids * (p8 * p9 * (i__1 * i__1) * e[2] + p8 * p9 * e[7] - p8 * e[7] * dpr9 * f - dpr8 * e[7] * p9 * f + dpr8 * e[7] * (d__1 * d__1) * dpr9 + dpz8 * e[8] * dpz9 * (d__2 * d__2)) / (d__3 * d__3);
	ae[349] -= poids * *n * (p8 * e[8] * dpz9 - p9 * e[4] * dpz8) / f;
/* Computing 2nd power */
	d__1 = f;
	ae[350] += poids * p9 * *n * (p9 * e[2] + p9 * e[7] - e[7] * dpr9 * f + e[1] * dpr9 * f) / (d__1 * d__1);
/* Computing 2nd power */
	d__1 = p9;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__2 = p9;
/* Computing 2nd power */
	d__3 = dpr9;
/* Computing 2nd power */
	d__4 = f;
/* Computing 2nd power */
	d__5 = dpz9;
/* Computing 2nd power */
	d__6 = f;
/* Computing 2nd power */
	d__7 = f;
	ae[351] += poids * (d__1 * d__1 * (i__1 * i__1) * e[2] + d__2 * d__2 * e[7] - p9 * 2. * e[7] * dpr9 * f + e[7] * (d__3 * d__3) * (d__4 * d__4) + e[8] * (d__5 * d__5) * (d__6 * d__6)) / (d__7 * d__7);
	ae[352] += poids * (p1 * e[4] * dpz9 + dpr1 * e[3] * dpz9 * f + dpz1 * e[6] * dpr9 * f) / f;
	ae[353] += poids * *n * (p1 * e[4] * dpz9 - dpz1 * e[8] * p9) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[354] += poids * (p1 * (i__1 * i__1) * e[8] * p9 + dpr1 * e[6] * dpr9 * (d__1 * d__1) + dpz1 * e[5] * dpz9 * (d__2 * d__2)) / (d__3 * d__3);
	ae[355] += poids * (dpz2 * e[6] * dpr9 * f + dpz9 * e[4] * p2 + dpz9 * e[3] * dpr2 * f) / f;
	ae[356] += poids * *n * (dpz9 * e[4] * p2 - dpz2 * e[8] * p9) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[357] += poids * (p2 * (i__1 * i__1) * e[8] * p9 + dpr2 * e[6] * dpr9 * (d__1 * d__1) + dpz2 * e[5] * dpz9 * (d__2 * d__2)) / (d__3 * d__3);
	ae[358] += poids * (dpz3 * e[6] * dpr9 * f + p3 * e[4] * dpz9 + dpz9 * e[3] * dpr3 * f) / f;
	ae[359] += poids * *n * (p3 * e[4] * dpz9 - dpz3 * e[8] * p9) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[360] += poids * (p3 * (i__1 * i__1) * e[8] * p9 + dpr3 * e[6] * dpr9 * (d__1 * d__1) + dpz3 * e[5] * dpz9 * (d__2 * d__2)) / (d__3 * d__3);
	ae[361] += poids * (dpz4 * e[6] * dpr9 * f + p4 * e[4] * dpz9 + dpz9 * e[3] * dpr4 * f) / f;
	ae[362] += poids * *n * (p4 * e[4] * dpz9 - dpz4 * e[8] * p9) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[363] += poids * (p4 * (i__1 * i__1) * e[8] * p9 + dpr4 * e[6] * dpr9 * (d__1 * d__1) + dpz4 * e[5] * dpz9 * (d__2 * d__2)) / (d__3 * d__3);
	ae[364] += poids * (dpz5 * e[6] * dpr9 * f + p5 * e[4] * dpz9 + dpz9 * dpr5 * e[3] * f) / f;
	ae[365] += poids * *n * (p5 * e[4] * dpz9 - dpz5 * e[8] * p9) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[366] += poids * (p5 * (i__1 * i__1) * e[8] * p9 + dpr5 * e[6] * dpr9 * (d__1 * d__1) + dpz5 * e[5] * dpz9 * (d__2 * d__2)) / (d__3 * d__3);
	ae[367] += poids * (dpz6 * e[6] * dpr9 * f + p6 * e[4] * dpz9 + dpz9 * e[3] * dpr6 * f) / f;
	ae[368] += poids * *n * (p6 * e[4] * dpz9 - dpz6 * e[8] * p9) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[369] += poids * (p6 * (i__1 * i__1) * e[8] * p9 + dpr6 * e[6] * dpr9 * (d__1 * d__1) + dpz6 * e[5] * dpz9 * (d__2 * d__2)) / (d__3 * d__3);
	ae[370] += poids * (dpz7 * e[6] * dpr9 * f + p7 * e[4] * dpz9 + dpz9 * e[3] * dpr7 * f) / f;
	ae[371] += poids * *n * (p7 * e[4] * dpz9 - dpz7 * e[8] * p9) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[372] += poids * (p7 * (i__1 * i__1) * e[8] * p9 + dpr7 * e[6] * dpr9 * (d__1 * d__1) + dpz7 * e[5] * dpz9 * (d__2 * d__2)) / (d__3 * d__3);
	ae[373] += poids * (dpz8 * e[6] * dpr9 * f + dpz9 * e[4] * p8 + dpz9 * e[3] * dpr8 * f) / f;
	ae[374] += poids * *n * (dpz9 * e[4] * p8 - dpz8 * e[8] * p9) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = f;
/* Computing 2nd power */
	d__2 = f;
/* Computing 2nd power */
	d__3 = f;
	ae[375] += poids * (p8 * (i__1 * i__1) * e[8] * p9 + dpr8 * e[6] * dpr9 * (d__1 * d__1) + dpz8 * e[5] * dpz9 * (d__2 * d__2)) / (d__3 * d__3);
	ae[376] += poids * dpz9 * (e[4] * p9 + e[3] * dpr9 * f + e[6] * dpr9 * f) / f;
	ae[377] += poids * p9 * *n * dpz9 * (e[4] - e[8]) / f;
/* Computing 2nd power */
	i__1 = *n;
/* Computing 2nd power */
	d__1 = p9;
/* Computing 2nd power */
	d__2 = dpr9;
/* Computing 2nd power */
	d__3 = f;
/* Computing 2nd power */
	d__4 = dpz9;
/* Computing 2nd power */
	d__5 = f;
/* Computing 2nd power */
	d__6 = f;
	ae[378] += poids * (i__1 * i__1 * e[8] * (d__1 * d__1) + e[6] * (d__2 * d__2) * (d__3 * d__3) + e[5] * (d__4 * d__4) * (d__5 * d__5)) / (d__6 * d__6);
/* L2: */
    }
} /* eraq2c_ */

