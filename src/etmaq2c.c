/* etmaq2c.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__8 = 8;
static int32 c__9 = 9;
static int32 c__2 = 2;

/* Subroutine */ int etmaq2c_(coor, car, iopt, ae)
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

    /* System generated locals */
    int32 i__1;

    /* Local variables */
    static int32 i__, j, l;
    extern /* Subroutine */ int e2aq2c_();
    static int32 n;
    static doublereal s, f1[9], f2[9], dfm1dp[144]	/* was [2][8][9] */;
    static int32 ip[16];
    static doublereal ss[36], poidel[9];
    extern /* Subroutine */ int plmasd_();
    static doublereal rho, dfm1[36]	/* was [4][9] */;

/* *************************************************************** */
/* BUT : MATRICE DE MASSE DE L ELEMENT AXI: QUAD AQ2C */
/* --- */
/* in : coor(noe,ndim) : coordonees R,Z  des 8 noeuds */
/*      car            : caracteristiques des materiaux */
/*                       car(1) = rho masse volumique */
/*      iopt           : ouvert si masse lumping ou autre ds futur */

/* out: ae            : matrice triangulaire sup */

/* programmeur : modulef */
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

    rho = car[1];
    n = 0;
    for (j = 1; j <= 8; ++j) {
	i__1 = j;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s = 0.;
	    for (l = 1; l <= 9; ++l) {
		s += poidel[l - 1] * p25[j + (l << 3) - 9] * p25[i__ + (l << 3) - 9];
/* L1: */
	    }
	    ++n;
	    ss[n - 1] = s * rho;
/* L2: */
	}
/* L3: */
    }

/*     PLONGER SS DANS AE */
/*     ------------------ */
    plmasd_(&c__8, &c__2, ss, ip, &ae[1]);

} /* etmaq2c_ */

