/* etmap2c.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__6 = 6;
static int32 c__7 = 7;
static int32 c__2 = 2;

/* Subroutine */ int etmap2c_(coor, car, iopt, ae)
doublereal *coor, *car;
int32 *iopt;
doublereal *ae;
{
    /* Initialized data */

    static doublereal poids[7] = { .062969590272413583,.062969590272413583,.062969590272413583,.066197076394253095,.066197076394253095,.066197076394253095,.11249999701976776 };
    static doublereal p25[42]	/* was [6][7] */ = { .47435260858553857,-.080768594191887185,-.080768594191887185,.32307437676754874,.041035826263138293,.32307437676754874,-.080768594191886977,.4743526085855384,-.080768594191887185,.32307437676754879,.32307437676754868,.041035826263138341,-.080768594191887033,-.080768594191887185,.4743526085855384,.041035826263138334,.32307437676754868,.3230743767675489,-.028074943223078796,-.028074943223078852,-.052583901102545349,.88413424176407262,
	    .11229977289231514,.11229977289231517,-.052583901102545571,-.028074943223078852,-.028074943223078852,.1122997728923154,.8841342417640724,.1122997728923154,-.028074943223078765,-.052583901102545349,-.028074943223078852,.11229977289231514,.11229977289231514,.88413424176407262,-.11111111773384862,-.11111110779974175,-.11111110779974175,.44444443119896703,.44444447093539807,.44444443119896703 };
    static doublereal dp25[84]	/* was [2][6][7] */ = { -2.1897079414123492,-2.1897079414123492,-.59485397070617462,0.,0.,-.59485397070617462,2.7845619121185238,-.40514602929382531,.40514602929382531,.40514602929382531,-.40514602929382531,2.7845619121185238,.59485397070617418,.59485397070617418,2.1897079414123488,0.,0.,-.59485397070617462,-2.7845619121185229,-3.1897079414123488,.40514602929382531,3.1897079414123488,-.40514602929382531,0.,.59485397070617418,.59485397070617418,-.59485397070617462,
	    0.,0.,2.1897079414123488,0.,-.40514602929382531,3.1897079414123488,.40514602929382531,-3.1897079414123488,-2.7845619121185229,-.88056825642046065,-.88056825642046065,.88056825642046021,0.,0.,-.76113651284092076,0.,-1.8805682564204602,.2388634871590792,1.8805682564204602,-.2388634871590792,1.6417047692613815,.76113651284092043,.76113651284092043,.88056825642046021,0.,0.,.88056825642046021,-1.6417047692613806,-1.8805682564204602,1.8805682564204602,1.8805682564204602,-1.8805682564204602,
	    -1.6417047692613806,-.88056825642046054,-.88056825642046054,-.76113651284092076,0.,0.,.88056825642046021,1.6417047692613813,-.2388634871590792,1.8805682564204602,.2388634871590792,-1.8805682564204602,0.,-.33333325386047363,-.33333325386047363,.33333337306976318,0.,0.,.33333337306976318,-1.1920928955078125e-7,-1.3333333730697631,1.3333333730697631,1.3333333730697631,-1.3333333730697631,-1.1920928955078125e-7 };

    /* System generated locals */
    int32 i__1;

    /* Local variables */
    static int32 i__, j, l;
    extern /* Subroutine */ int e2ap2c_();
    static int32 n;
    static doublereal s, f1[7], f2[7], dfm1dp[84]	/* was [2][6][7] */;
    static int32 ip[12];
    static doublereal ss[21], poidel[7];
    extern /* Subroutine */ int plmasd_();
    static doublereal rho, dfm1[28]	/* was [4][7] */;

/* *************************************************************** */
/* BUT : CALCUL DE LA MATRICE DE MASSE DE L ELEMENT AXI: TRIA AP2C */
/* --- */
/* in : coor(noe,ndim) : coor. 3 sommets + 3 milieux aretes */
/*      car            : caracteristiques des materiaux */
/*                       car(1) = rho masse volumique */
/*      iopt           : ouvert si masse lumping ou autre ds futur */

/* out: ae            : matrice triangulaire sup */

/* programmeur : modulef */
/* ............................................................... */
/* 2P25 -- XYNPI: coordonnees pt. int. numeriques (element reference) */
    /* Parameter adjustments */
    --ae;
    --car;
    coor -= 7;

    /* Function Body */
/*     -- POIDS: poids du schema d'integration numerique. */
/*     -- Valeurs des Polynomes de base aux pt. int. numerique. */
/*     -- Valeurs Derive'es des Poly. de base aux pt int. numerique. */

/*     F1 , F2 , DFM1 , POIDEL , IP , DP PRETS A L EMPLOI */
/*     -------------------------------------------------- */
    e2ap2c_(&c__6, &c__7, poids, p25, dp25, ip, f1, f2, dfm1dp, poidel, dfm1, &coor[7]);

    rho = car[1];
    n = 0;
    for (j = 1; j <= 6; ++j) {
	i__1 = j;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s = 0.;
	    for (l = 1; l <= 7; ++l) {
		s += poidel[l - 1] * p25[j + l * 6 - 7] * p25[i__ + l * 6 - 7];
/* L2: */
	    }
	    ++n;
	    ss[n - 1] = s * rho;
/* L3: */
	}
/* L4: */
    }

/*     PLONGER SS DANS AE */
/*     ------------------ */
    plmasd_(&c__6, &c__2, ss, ip, &ae[1]);

} /* etmap2c_ */

