/* etrap1d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__3 = 3;
static int32 c__1 = 1;
static int32 c__2 = 2;
static int32 c_n1 = -1;

/* Subroutine */ int etrap1d_(coor, car, iopt, ae)
doublereal *coor, *car;
int32 *iopt;
doublereal *ae;
{
    /* Initialized data */

    static doublereal poids[1] = { .5 };
    static doublereal p13[3]	/* was [3][1] */ = { .33333333333333333,.33333333333333333,.33333333333333333 };
    static doublereal dp13[6]	/* was [2][3][1] */ = { -1.,-1.,1.,0.,0.,1. };

    static doublereal elas[10];
    extern /* Subroutine */ int tab0d_(), tab1d_();
    static int32 i__;
    extern /* Subroutine */ int e1ap1d_();
    static int32 l;
    static doublereal f1[1], f2[1], g1[9]	/* was [3][3] */, g2[3], g3[6]	/* was [3][2] */, g5[6]	/* was [2][3] */, g4[6]	/* was [3][2] */, dfm1dp[18]	/* was [2][3][3] */, a11[1], a12[2], a13[2], a22[4]	/* was [2][2] */, a23[4]	/* was [2][2] */, a33[3];
    static int32 ip[6];
    extern /* Subroutine */ int tabaxd_();
    static doublereal poidel[1];
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
/*     -- XYNPI: coordonnees pt. int. numeriques (element reference) */
    /* Parameter adjustments */
    --ae;
    --car;
    coor -= 4;

    /* Function Body */
/*     -- POIDSG: poids du schema d'integration numerique. */
/*     -- Valeurs des Polynomes de base aux pt. int. numerique. */
/*     -- Valeurs Derive'es des Poly. de base aux pt int. numerique. */

    e1ap1d_(&c__3, &c__1, poids, p13, dp13, ip, f1, f2, dfm1dp, poidel, dfm1, &coor[4]);

    hookax_(iopt, &car[1], elas);

/*     -- INITIALISATION DE AE */
    for (i__ = 1; i__ <= 21; ++i__) {
	ae[i__] = 0.;
/* L1: */
    }

/*     CONSTRUCTION DES MATRICES   (A11 , A12 , A13) */
/*     -------------------------   (      A22 , A23)  =  (TD * E * D) */
/*                                 (            A33) */
    for (l = 1; l <= 1; ++l) {
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
	ab0d_(&c__1, &c__1, &c__3, a11, &p13[l * 3 - 3], g2);
	ab1d_(&c__1, &c__2, &c__3, a12, &dfm1dp[(l * 3 + 1 << 1) - 8], g2);

	tab0d_(&c__3, &c__1, &c__3, &p13[l * 3 - 3], g2, g1);

	tab0d_(&c__2, &c__1, &c__3, a12, &p13[l * 3 - 3], g5);
	ab1d_(&c__2, &c__2, &c__3, a22, &dfm1dp[(l * 3 + 1 << 1) - 8], g5);

	tab1d_(&c__3, &c__2, &c__3, &dfm1dp[(l * 3 + 1 << 1) - 8], g5, g1);
/*        ON PLONGE G1 DANS AE */
	plonad_(&c_n1, &c__3, &c__3, &c__1, &c__1, ip, g1, g1, &ae[1]);

/*        CALCUL DE G2 */
/*        ------------ */
	tab0d_(&c__3, &c__1, &c__2, &p13[l * 3 - 3], a13, g3);
	tab1d_(&c__3, &c__2, &c__2, &dfm1dp[(l * 3 + 1 << 1) - 8], a23, g3);

	ab0d_(&c__3, &c__2, &c__3, g3, &dfm1dp[(l * 3 + 1 << 1) - 8], g1);
/*        ON PLONGE G2 DANS AE */
	plonad_(&c_n1, &c__3, &c__3, &c__1, &c__2, ip, g1, g1, &ae[1]);

/*        CALCUL DE G3 */
/*        ------------ */
	tabaxd_(&c__3, &c__2, &dfm1dp[(l * 3 + 1 << 1) - 8], a33, g1, g4, &c__1);
/*        ON PLONGE G3 DANS AE */
	plonad_(&c__1, &c__3, &c__3, &c__2, &c__2, ip, g1, g1, &ae[1]);
/* L2: */
    }
} /* etrap1d_ */

