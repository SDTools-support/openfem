/* etmaq1d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__4 = 4;
static int32 c__2 = 2;

/* Subroutine */ int etmaq1d_(coor, car, iopt, ae)
doublereal *coor, *car;
int32 *iopt;
doublereal *ae;
{
    /* Initialized data */

    static doublereal poids[4] = { .25,.25,.25,.25 };
    static doublereal q13[16]	/* was [4][4] */ = { .62200846792814612,.16666666666666662,.044658198738520435,.16666666666666662,.16666666666666662,.62200846792814623,.16666666666666662,.044658198738520449,.044658198738520504,.16666666666666662,.62200846792814623,.16666666666666662,.16666666666666668,.044658198738520449,.16666666666666662,.62200846792814623 };
    static doublereal dq13[32]	/* was [2][4][4] */ = { -.78867513459481286,-.78867513459481286,.78867513459481286,-.21132486540518707,.21132486540518707,.21132486540518707,-.21132486540518707,.78867513459481286,-.78867513459481286,-.21132486540518713,.78867513459481286,-.78867513459481286,.21132486540518707,.78867513459481286,-.21132486540518707,.21132486540518713,-.21132486540518713,-.21132486540518713,.21132486540518713,-.78867513459481286,.78867513459481286,.78867513459481286,
	    -.78867513459481286,.21132486540518713,-.21132486540518713,-.78867513459481286,.21132486540518713,-.21132486540518707,.78867513459481286,.21132486540518707,-.78867513459481286,.78867513459481286 };

    /* System generated locals */
    int32 i__1;

    /* Local variables */
    static int32 i__, j;
    extern /* Subroutine */ int e1aq1c_();
    static int32 l, n;
    static doublereal s, f1[5], f2[5], dfm1dp[32]	/* was [2][4][4] */;
    static int32 ip[8];
    static doublereal ss[10], poidel[4];
    extern /* Subroutine */ int plmasd_();
    static doublereal rho, dfm1[16]	/* was [4][4] */;

/* *************************************************************** */
/* BUT : CALCUL DE LA MATRICE DE MASSE DE L ELEMENT AXI: quad aq1d */
/* --- */
/* in : coor(noe,ndim) : coordones R(4), Z(4) des 4 sommets. */
/*      car            : caracteristiques des materiaux */
/*                       car(1) = rho masse volumique */
/*      iopt           : ouvert si masse lumping ou autre ds futur */

/* out: ae            : matrice triangulaire sup */

/* programmeur : modulef */
/* ............................................................... */
/* 2Q13 -- POIDS: poids du schema d'integration numerique. */
    /* Parameter adjustments */
    --ae;
    --car;
    coor -= 5;

    /* Function Body */
/*     -- XYNPI: coordonnees pt. int. numeriques (element reference) */
/*     -- Valeurs des Polynomes de base aux pt. int. numerique. */
/*     -- Valeurs Derive'es des Poly. de base aux pt int. numerique. */

/*     F1 , F2 , DFM1 , POIDEL , IP , DFM1DP PRETS A L EMPLOI */

    e1aq1c_(&c__4, &c__4, poids, q13, dq13, ip, f1, f2, dfm1dp, poidel, dfm1, &coor[5]);

    rho = car[1];
    n = 0;
    for (j = 1; j <= 4; ++j) {
	i__1 = j;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s = 0.;
	    for (l = 1; l <= 4; ++l) {
		s += poidel[l - 1] * q13[j + (l << 2) - 5] * q13[i__ + (l << 2) - 5];
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
    plmasd_(&c__4, &c__2, ss, ip, &ae[1]);

} /* etmaq1d_ */

