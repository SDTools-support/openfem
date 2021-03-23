/* etmap1d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__3 = 3;
static int32 c__1 = 1;
static int32 c__2 = 2;

/* Subroutine */ int etmap1d_(coor, car, iopt, ae)
doublereal *coor, *car;
int32 *iopt;
doublereal *ae;
{
    /* Initialized data */

    static doublereal poids[1] = { .5 };
    static doublereal p13[3]	/* was [3][1] */ = { .33333333333333333,.33333333333333333,.33333333333333333 };
    static doublereal dp13[6]	/* was [2][3][1] */ = { -1.,-1.,1.,0.,0.,1. };

    /* System generated locals */
    int32 i__1;

    /* Local variables */
    static int32 i__, j;
    extern /* Subroutine */ int e1ap1d_();
    static int32 l, n;
    static doublereal s, f1[1], f2[1], dfm1dp[18]	/* was [2][3][3] */;
    static int32 ip[6];
    static doublereal ss[6], poidel[1];
    extern /* Subroutine */ int plmasd_();
    static doublereal rho, dfm1[16]	/* was [4][4] */;

/* *************************************************************** */
/* BUT : CALCUL DE LA MATRICE DE MASSE DE L ELEMENT AXI: TRIA AP1D */
/* --- */
/* in : coor(noe,ndim) : coordones R(3), Z(3) des 3 sommets. */
/*      car            : caracteristiques des materiaux */
/*                       car(1) = rho masse volumique */
/*      iopt           : ouvert si masse lumping ou autre ds futur */

/* out: ae            : matrice triangulaire sup */

/* programmeur : modulef */
/* ............................................................... */
/*     -- XYNPI: coordonnees pt. int. numeriques (element reference) */
    /* Parameter adjustments */
    --ae;
    --car;
    coor -= 4;

    /* Function Body */
/*     -- POIDS: poids du schema d'integration numerique. */
/*     -- Valeurs des Polynomes de base aux pt. int. numerique. */
/*     -- Valeurs Derive'es des Poly. de base aux pt int. numerique. */

    e1ap1d_(&c__3, &c__1, poids, p13, dp13, ip, f1, f2, dfm1dp, poidel, dfm1, &coor[4]);

/*     F1 , F2 , DFM1 , POIDEL , IP ,  DFM1DP PRETS A L EMPLOI */
/*     -------------------------------------------------- */
    rho = car[1];
    n = 0;
    for (j = 1; j <= 3; ++j) {
	i__1 = j;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s = 0.;
	    for (l = 1; l <= 1; ++l) {
		s += poidel[l - 1] * p13[j + l * 3 - 4] * p13[i__ + l * 3 - 4];
/* L5: */
	    }
	    ++n;
	    ss[n - 1] = s * rho;
/* L4: */
	}
/* L3: */
    }

/*     PLONGER SS DANS AE */
/*     ------------------ */
    plmasd_(&c__3, &c__2, ss, ip, &ae[1]);

} /* etmap1d_ */

