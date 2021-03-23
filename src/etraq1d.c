/* etraq1d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__4 = 4;
static int32 c__1 = 1;
static int32 c__2 = 2;
static int32 c_n1 = -1;

/* Subroutine */ int etraq1d_(coor, car, iopt, ae)
doublereal *coor, *car;
int32 *iopt;
doublereal *ae;
{
    /* Initialized data */

    static doublereal poids[4] = { .25,.25,.25,.25 };
    static doublereal q13[16]	/* was [4][4] */ = { .62200846792814612,.16666666666666662,.044658198738520435,.16666666666666662,.16666666666666662,.62200846792814623,.16666666666666662,.044658198738520449,.044658198738520504,.16666666666666662,.62200846792814623,.16666666666666662,.16666666666666668,.044658198738520449,.16666666666666662,.62200846792814623 };
    static doublereal dq13[32]	/* was [2][4][4] */ = { -.78867513459481286,-.78867513459481286,.78867513459481286,-.21132486540518707,.21132486540518707,.21132486540518707,-.21132486540518707,.78867513459481286,-.78867513459481286,-.21132486540518713,.78867513459481286,-.78867513459481286,.21132486540518707,.78867513459481286,-.21132486540518707,.21132486540518713,-.21132486540518713,-.21132486540518713,.21132486540518713,-.78867513459481286,.78867513459481286,.78867513459481286,
	    -.78867513459481286,.21132486540518713,-.21132486540518713,-.78867513459481286,.21132486540518713,-.21132486540518707,.78867513459481286,.21132486540518707,-.78867513459481286,.78867513459481286 };

    static doublereal elas[10];
    extern /* Subroutine */ int tab0d_(), tab1d_();
    static int32 i__;
    extern /* Subroutine */ int e1aq1c_();
    static int32 l;
    static doublereal f1[5], f2[5], g1[16]	/* was [4][4] */, g2[4], g3[8]	/* was [4][2] */, g5[8]	/* was [2][4] */, g4[8]	/* was [4][2] */;
    extern /* Subroutine */ int taba8d_();
    static doublereal dfm1dp[32]	/* was [2][4][4] */, a11[1], a12[2], a13[2], a22[4]	/* was [2][2] */, a23[4]	/* was [2][2] */, a33[3];
    static int32 ip[8];
    static doublereal poidel[4];
    extern /* Subroutine */ int plonad_(), hookax_();
    static doublereal aux;
    extern /* Subroutine */ int ab0d_(), ab1d_();
    static doublereal dfm1[16]	/* was [4][4] */;

/*  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  BUT  : CALCUL DE LA MATRICE DE RAIDEUR DE L ELEMENT TRIA AP1D */
/*  --- */
/* in : coor(noe,ndim) : coordonnees R(3), Z(3) des 3 sommets. */
/*      car            : caracteristiques des materiaux */
/*           if(iopt .lt. 3) then */
/*             car(1) = young */
/*             car(2) = poisson */
/*           else */
/*             car(1) = E_1  (Young radial     -> E_r     ) */
/*             car(2) = nu_1 (poisson          -> Nu_r    ) */
/*             car(3) = E_2  (Young axial      -> E_z     ) */
/*             car(4) = nu_2 (poisson          -> Nu_z    ) */
/*             car(5) = E_3  (Young Tangentiel -> E_theta ) */
/*           end if */
/* out: ae            : matrice triangulaire sup */
/* ............................................................... */
/*  programmeur : modulef */
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

    e1aq1c_(&c__4, &c__4, poids, q13, dq13, ip, f1, f2, dfm1dp, poidel, dfm1, &coor[5]);

    hookax_(iopt, &car[1], elas);
/*     -- INITIALISATION DE AE */
    for (i__ = 1; i__ <= 36; ++i__) {
	ae[i__] = 0.;
/* L1: */
    }

/*     CONSTRUCTION DES MATRICES   (A11 , A12 , A13) */
/*     -------------------------   (      A22 , A23)  =  (TD * E * D) */
/*                                 (            A33) */
    for (l = 1; l <= 4; ++l) {
	aux = poidel[l - 1] / f1[l - 1];
	a11[0] = aux * elas[5] / f1[l - 1];
	a12[0] = aux * elas[4];
	a12[1] = aux * elas[8];
	a13[0] = aux * elas[8];
	a13[1] = aux * elas[3];
	aux = poidel[l - 1];
	a22[0] = aux * elas[2];
	a22[1] = aux * elas[7];
	a22[2] = a22[1];
	a22[3] = aux * elas[9];
	a23[0] = a22[1];
	a23[1] = a22[3];
	a23[2] = aux * elas[1];
	a23[3] = aux * elas[6];
	a33[0] = a22[3];
	a33[1] = a23[3];
	a33[2] = aux * elas[0];

/*        CALCUL DE G1 */
/*        ------------ */
	ab0d_(&c__1, &c__1, &c__4, a11, &q13[(l << 2) - 4], g2);
	ab1d_(&c__1, &c__2, &c__4, a12, &dfm1dp[((l << 2) + 1 << 1) - 10], g2);

	tab0d_(&c__4, &c__1, &c__4, &q13[(l << 2) - 4], g2, g1);

	tab0d_(&c__2, &c__1, &c__4, a12, &q13[(l << 2) - 4], g5);
	ab1d_(&c__2, &c__2, &c__4, a22, &dfm1dp[((l << 2) + 1 << 1) - 10], g5);

	tab1d_(&c__4, &c__2, &c__4, &dfm1dp[((l << 2) + 1 << 1) - 10], g5, g1);
/*        ON PLONGE G1 DANS AE */
	plonad_(&c_n1, &c__4, &c__4, &c__1, &c__1, ip, g1, g1, &ae[1]);

/*        CALCUL DE G2 */
/*        ------------ */
	tab0d_(&c__4, &c__1, &c__2, &q13[(l << 2) - 4], a13, g3);
	tab1d_(&c__4, &c__2, &c__2, &dfm1dp[((l << 2) + 1 << 1) - 10], a23, g3);

	ab0d_(&c__4, &c__2, &c__4, g3, &dfm1dp[((l << 2) + 1 << 1) - 10], g1);
/*        ON PLONGE G2 DANS AE */
	plonad_(&c_n1, &c__4, &c__4, &c__1, &c__2, ip, g1, g1, &ae[1]);

/*        CALCUL DE G3 */
/*        ------------ */
	taba8d_(&c__4, &c__2, &dfm1dp[((l << 2) + 1 << 1) - 10], a33, g1, g4);
/*        ON PLONGE G3 DANS AE */
	plonad_(&c__1, &c__4, &c__4, &c__2, &c__2, ip, g1, g1, &ae[1]);
/* L2: */
    }
} /* etraq1d_ */

