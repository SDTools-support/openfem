/* etraq2c.f -- translated by f2c (version 19961017).
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
static int32 c_n1 = -1;

/* Subroutine */ int etraq2c_(coor, car, iopt, ae)
doublereal *coor, *car;
int32 *iopt;
doublereal *ae;
{
    /* Initialized data */

    static doublereal poids[9] = { .077160493827160503,.12345679012345678,.077160493827160503,.12345679012345678,.19753086419753085,.12345679012345678,.077160493827160503,.12345679012345678,.077160493827160503 };
    static doublereal p25[72]	/* was [8][9] */ = { .43237900077244512,-.1,-.032379000772445008,-.099999999999999991,.35491933384829665,.045080666151703314,.045080666151703314,.35491933384829665,-.099999999999999977,-.099999999999999977,-.099999999999999977,-.1,.8872983346207417,.19999999999999998,.11270166537925829,.19999999999999998,-.10000000000000003,.43237900077244512,-.099999999999999977,-.032379000772445015,.35491933384829665,.35491933384829665,.045080666151703308,.045080666151703314,
	    -.099999999999999991,-.099999999999999991,-.099999999999999991,-.099999999999999991,.19999999999999998,.11270166537925829,.19999999999999998,.8872983346207417,-.25,-.25,-.25,-.25,.5,.5,.5,.5,-.10000000000000019,-.099999999999999977,-.099999999999999977,-.099999999999999977,.19999999999999995,.8872983346207417,.19999999999999995,.11270166537925829,-.099999999999999977,-.032379000772444987,-.1,.43237900077244512,.045080666151703308,.045080666151703308,.35491933384829665,
	    .3549193338482966,-.10000000000000019,-.099999999999999977,-.099999999999999977,-.099999999999999977,.11270166537925829,.19999999999999995,.8872983346207417,.19999999999999995,-.032379000772445598,-.099999999999999866,.43237900077244484,-.10000000000000008,.045080666151703141,.35491933384829676,.35491933384829676,.045080666151703141 };
    static doublereal dp25[144]	/* was [2][8][9] */ = { -2.061895003862225,-2.061895003862225,-.68729833462074174,-.087298334620741685,-.26189500386222502,-.26189500386222502,-.087298334620741685,-.68729833462074174,2.7491933384829669,-.39999999999999996,.39999999999999996,.34919333848296674,.34919333848296674,.39999999999999996,-.39999999999999996,2.7491933384829669,-.68729833462074152,-.7745966692414834,.68729833462074174,-.7745966692414834,-.087298334620741657,-.7745966692414834,
	    .087298334620741685,-.7745966692414834,0.,-1.,.39999999999999996,1.5491933384829668,0.,1.,-.39999999999999996,1.5491933384829668,.68729833462074196,-.087298334620742101,2.0618950038622254,-2.061895003862225,.087298334620741713,-.68729833462074174,.26189500386222508,-.26189500386222497,-2.7491933384829669,-.39999999999999991,.39999999999999996,2.7491933384829669,-.34919333848296674,.39999999999999991,-.39999999999999996,.34919333848296674,-.7745966692414834,-.68729833462074174,
	    -.7745966692414834,.087298334620741685,-.7745966692414834,-.087298334620741685,-.7745966692414834,.68729833462074174,1.5491933384829668,-.39999999999999996,1.,0.,1.5491933384829668,.39999999999999996,-1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-1.,1.,0.,0.,1.,-1.,0.,.7745966692414834,.087298334620741213,.7745966692414834,-.68729833462074174,.7745966692414834,.68729833462074174,.7745966692414834,-.087298334620741657,-1.5491933384829668,-.39999999999999991,1.,0.,-1.5491933384829668,
	    .39999999999999991,-1.,0.,-.087298334620742157,.68729833462074174,-.26189500386222502,.26189500386222502,-.68729833462074174,.087298334620741685,-2.0618950038622254,2.061895003862225,.34919333848296674,-.39999999999999996,.39999999999999991,-.34919333848296674,2.7491933384829669,.39999999999999996,-.39999999999999991,-2.7491933384829669,.087298334620741213,.7745966692414834,-.087298334620741657,.7745966692414834,.68729833462074174,.7745966692414834,-.68729833462074196,
	    .7745966692414834,0.,-1.,.39999999999999991,-1.5491933384829668,0.,1.,-.39999999999999991,-1.5491933384829668,.26189500386222475,.26189500386222452,.087298334620741879,.68729833462074174,2.061895003862225,2.061895003862225,.68729833462074152,.087298334620741657,-.34919333848296663,-.39999999999999991,.39999999999999991,-2.7491933384829669,-2.7491933384829669,.39999999999999991,-.39999999999999991,-.34919333848296663 };

    static doublereal elas[10];
    extern /* Subroutine */ int tab0d_(), tab1d_();
    static int32 i__, l;
    extern /* Subroutine */ int e2aq2c_();
    static doublereal f1[9], f2[9], g1[64]	/* was [8][8] */, g2[8], g3[16]	/* was [8][2] */, g5[16]	/* was [2][8] */, g4[16]	/* was [8][2] */, dfm1dp[144]	/* was [2][8][9] */, a11[1], a12[2], a13[2], a22[4]	/* was [2][2] */, a23[4]	/* was [2][2] */, a33[3];
    static int32 ip[16];
    extern /* Subroutine */ int tabaxd_();
    static doublereal poidel[9];
    extern /* Subroutine */ int plonad_(), hookax_();
    static doublereal aux;
    extern /* Subroutine */ int ab0d_(), ab1d_();
    static doublereal dfm1[36]	/* was [4][9] */;

/*  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  BUT  : MATRICE DE RAIDEUR DE L ELEMENT AXI: QUAD AQ2C */
/*  --- */
/* in : coor(noe,ndim) : coordonnees R(3), Z(3) des 3 sommets. */
/*      car            : caracteristiques des materiaux */
/*           if(iopt .eq. 1) then */
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
/* 2Q25 -- XYNPI: coordonnees pt. int. numeriques (element reference) */
    /* Parameter adjustments */
    --ae;
    --car;
    coor -= 9;

    /* Function Body */
/*     -- POIDS: poids du schema d'integration numerique. */
/*     -- Valeurs des Polynomes de base aux pt. int. numerique. */
/*     -- Valeurs Derive'es des Poly. de base aux pt int. numerique. */

/*     F1 , F2 , DFM1 , POIDEL , IP , DP PRETS A L EMPLOI */

    e2aq2c_(&c__8, &c__9, poids, p25, dp25, ip, f1, f2, dfm1dp, poidel, dfm1, &coor[9]);

    hookax_(iopt, &car[1], elas);

    for (i__ = 1; i__ <= 136; ++i__) {
	ae[i__] = 0.;
/* L1: */
    }

/*     CONSTRUCTION DES MATRICES   (A11 , A12 , A13) */
/*     -------------------------   (      A22 , A23)  =  (TD * E * D) */
/*                                 (            A33) */
    for (l = 1; l <= 9; ++l) {
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
	ab0d_(&c__1, &c__1, &c__8, a11, &p25[(l << 3) - 8], g2);
	ab1d_(&c__1, &c__2, &c__8, a12, &dfm1dp[((l << 3) + 1 << 1) - 18], g2);

	tab0d_(&c__8, &c__1, &c__8, &p25[(l << 3) - 8], g2, g1);

	tab0d_(&c__2, &c__1, &c__8, a12, &p25[(l << 3) - 8], g5);
	ab1d_(&c__2, &c__2, &c__8, a22, &dfm1dp[((l << 3) + 1 << 1) - 18], g5);

	tab1d_(&c__8, &c__2, &c__8, &dfm1dp[((l << 3) + 1 << 1) - 18], g5, g1);
/*        ON PLONGE G1 DANS AE */
	plonad_(&c_n1, &c__8, &c__8, &c__1, &c__1, ip, g1, g1, &ae[1]);

/*        CALCUL DE G2 */
/*        ------------ */
	tab0d_(&c__8, &c__1, &c__2, &p25[(l << 3) - 8], a13, g3);
	tab1d_(&c__8, &c__2, &c__2, &dfm1dp[((l << 3) + 1 << 1) - 18], a23, g3);

	ab0d_(&c__8, &c__2, &c__8, g3, &dfm1dp[((l << 3) + 1 << 1) - 18], g1);
/*        ON PLONGE G2 DANS AE */
	plonad_(&c_n1, &c__8, &c__8, &c__1, &c__2, ip, g1, g1, &ae[1]);

/*        CALCUL DE G3 */
/*        ------------ */
	tabaxd_(&c__8, &c__2, &dfm1dp[((l << 3) + 1 << 1) - 18], a33, g1, g4, &c__1);
/*        ON PLONGE G3 DANS AE */
	plonad_(&c__1, &c__8, &c__8, &c__2, &c__2, ip, g1, g1, &ae[1]);
/* L5: */
    }
} /* etraq2c_ */

