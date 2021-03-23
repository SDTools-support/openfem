/* etr2q2c.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__8 = 8;
static int32 c__2 = 2;

/* Subroutine */ int etr2q2c_(coor, car, iopt, ae)
doublereal *coor, *car;
int32 *iopt;
doublereal *ae;
{
    /* Initialized data */

    static int32 ijt[16] = { 1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16 };
    static int32 iblt[3] = { 0,0,8 };
    static int32 jblt[3] = { 0,8,8 };
    static doublereal poids[9] = { .077160493827160503,.12345679012345678,.077160493827160503,.12345679012345678,.19753086419753085,.12345679012345678,.077160493827160503,.12345679012345678,.077160493827160503 };
    static doublereal dp25[144]	/* was [2][8][9] */ = { -2.061895003862225,-2.061895003862225,-.68729833462074174,-.087298334620741685,-.26189500386222502,-.26189500386222502,-.087298334620741685,-.68729833462074174,2.7491933384829669,-.39999999999999996,.39999999999999996,.34919333848296674,.34919333848296674,.39999999999999996,-.39999999999999996,2.7491933384829669,-.68729833462074152,-.7745966692414834,.68729833462074174,-.7745966692414834,-.087298334620741657,-.7745966692414834,
	    .087298334620741685,-.7745966692414834,0.,-1.,.39999999999999996,1.5491933384829668,0.,1.,-.39999999999999996,1.5491933384829668,.68729833462074196,-.087298334620742101,2.0618950038622254,-2.061895003862225,.087298334620741713,-.68729833462074174,.26189500386222508,-.26189500386222497,-2.7491933384829669,-.39999999999999991,.39999999999999996,2.7491933384829669,-.34919333848296674,.39999999999999991,-.39999999999999996,.34919333848296674,-.7745966692414834,-.68729833462074174,
	    -.7745966692414834,.087298334620741685,-.7745966692414834,-.087298334620741685,-.7745966692414834,.68729833462074174,1.5491933384829668,-.39999999999999996,1.,0.,1.5491933384829668,.39999999999999996,-1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-1.,1.,0.,0.,1.,-1.,0.,.7745966692414834,.087298334620741213,.7745966692414834,-.68729833462074174,.7745966692414834,.68729833462074174,.7745966692414834,-.087298334620741657,-1.5491933384829668,-.39999999999999991,1.,0.,-1.5491933384829668,
	    .39999999999999991,-1.,0.,-.087298334620742157,.68729833462074174,-.26189500386222502,.26189500386222502,-.68729833462074174,.087298334620741685,-2.0618950038622254,2.061895003862225,.34919333848296674,-.39999999999999996,.39999999999999991,-.34919333848296674,2.7491933384829669,.39999999999999996,-.39999999999999991,-2.7491933384829669,.087298334620741213,.7745966692414834,-.087298334620741657,.7745966692414834,.68729833462074174,.7745966692414834,-.68729833462074196,
	    .7745966692414834,0.,-1.,.39999999999999991,-1.5491933384829668,0.,1.,-.39999999999999991,-1.5491933384829668,.26189500386222475,.26189500386222452,.087298334620741879,.68729833462074174,2.061895003862225,2.061895003862225,.68729833462074152,.087298334620741657,-.34919333848296663,-.39999999999999991,.39999999999999991,-2.7491933384829669,-2.7491933384829669,.39999999999999991,-.39999999999999991,-.34919333848296663 };

    /* System generated locals */
    int32 i__1;

    /* Local variables */
    static int32 ifin;
    static doublereal c__, e[6], fdedf[36]	/* was [4][9] */;
    static int32 i__, j, k, ibloc;
    static doublereal delta[9], dfinv[36]	/* was [4][9] */;
    static int32 i1, i2, i3, i4, j1;
    extern /* Subroutine */ int taba2d_();
    static doublereal young, unmnu;
    extern /* Subroutine */ int taba6d_();
    static doublereal f11[9], f12[9], f21[9], f22[9];
    static int32 ij, kk, kin;
    static doublereal poisson, aux1[64];

/*  .................................................................... */
/* but : matrice de rigidite de l element membrane:  QUAD 2Q2C */
/* --- */
/* in : coor(noe,ndim) : coordonees 8 noeuds */

/*      iopt = 1 isotrope Deformations Planes */
/*           = 2 isotrope Contraintes  Planes */
/*           = sinon anisotrope */
/*      car(6): caracteristiques des materiaux */
/*              if(iopt .eq. 1 .or. iopt.eq. 2) then */
/*                car(1) = young */
/*                car(2) = poisson */
/*              else */
/*                car: E11, E12, E22, E13, E23, E33 avec */

/*                     E11   E12   E13 */
/*                           E22   E23 */
/*                                 E33 */
/*              end if */
/* out: ae(136)        : matrice triangulaire sup */

/* programmeur : modulef */
/* ............................................................... */
    /* Parameter adjustments */
    --ae;
    --car;
    coor -= 9;

    /* Function Body */
/* 2Q25 -- XYNPI: coordonnees pt. int. numeriques (element reference) */
/*     -- POIDS: poids du schema d'integration numerique. */
/*     -- Valeurs des Polynomes de base aux pt. int. numerique. */
/*     -- Valeurs Derive'es des Poly. de base aux pt int. numerique. */

    for (j = 1; j <= 136; ++j) {
	ae[j] = 0.;
/* L20: */
    }

    if (*iopt == 1) {
/*  --    CONTRAINTES PLANES     (ISOTROPE)     ----- */
	young = car[1];
	poisson = car[2];
	c__ = young / (1. - poisson * poisson);
	e[0] = c__;
	e[1] = c__ * poisson;
	e[2] = c__;
	e[3] = 0.;
	e[4] = 0.;
	e[5] = c__ * (1. - poisson) / 2.;
    } else if (*iopt == 2) {
/*  --    DEFORMATIONS PLANES (ISOTROPE)     ----- */
	young = car[1];
	poisson = car[2];
	unmnu = 1. - poisson;
	c__ = young * unmnu;
	c__ /= (poisson + 1.) * (1. - poisson * 2.);
	e[0] = c__;
	e[1] = poisson * c__ / unmnu;
	e[2] = c__;
	e[3] = 0.;
	e[4] = 0.;
	e[5] = c__ * (1. - poisson * 2.) / (unmnu * 2.);
    } else {
/*  --    CAS ANISOTROPE     ----- */
	for (i__ = 1; i__ <= 6; ++i__) {
	    e[i__ - 1] = car[i__];
/* L1: */
	}
    }

/*    -----   CALCUL DE DELTA AUX NPI NOEUDS  ----- */

    for (k = 1; k <= 9; ++k) {
	f11[k - 1] = 0.;
	f12[k - 1] = 0.;
	f21[k - 1] = 0.;
	f22[k - 1] = 0.;
/* L2: */
    }
    for (k = 1; k <= 9; ++k) {
	for (i__ = 1; i__ <= 8; ++i__) {
	    f11[k - 1] += dp25[(i__ + (k << 3) << 1) - 18] * coor[i__ + 8];
	    f12[k - 1] += dp25[(i__ + (k << 3) << 1) - 17] * coor[i__ + 8];
	    f21[k - 1] += dp25[(i__ + (k << 3) << 1) - 18] * coor[i__ + 16];
	    f22[k - 1] += dp25[(i__ + (k << 3) << 1) - 17] * coor[i__ + 16];
/* L3: */
	}
	delta[k - 1] = f11[k - 1] * f22[k - 1] - f12[k - 1] * f21[k - 1];
/* L4: */
    }

/*  ----       CALCUL DE DFINV :  INVERSE DE DF    ---- */

    for (k = 1; k <= 9; ++k) {
	dfinv[(k << 2) - 4] = f22[k - 1];
	dfinv[(k << 2) - 3] = -f12[k - 1];
	dfinv[(k << 2) - 2] = -f21[k - 1];
	dfinv[(k << 2) - 1] = f11[k - 1];
/* L5: */
    }

/*  -----      COEFFICIENT DE LA MATRICE ELEMENTAIRE     ----- */

    for (ibloc = 1; ibloc <= 3; ++ibloc) {
	if (*iopt == 1 || *iopt == 2) {
/*         -- isotrope Deformations planes ou Contraintes planes */
	    if (ibloc == 1) {
/*           -- BLOC  1,1 DE    TDF * TD * E * D * DF   (ISOTROPE) */
		for (j = 1; j <= 9; ++j) {
		    fdedf[(j << 2) - 4] = e[0] * dfinv[(j << 2) - 4] * dfinv[(j << 2) - 4] + e[5] * dfinv[(j << 2) - 3] * dfinv[(j << 2) - 3];
		    fdedf[(j << 2) - 3] = e[0] * dfinv[(j << 2) - 2] * dfinv[(j << 2) - 4] + e[5] * dfinv[(j << 2) - 1] * dfinv[(j << 2) - 3];
		    fdedf[(j << 2) - 2] = fdedf[(j << 2) - 3];
		    fdedf[(j << 2) - 1] = e[0] * dfinv[(j << 2) - 2] * dfinv[(j << 2) - 2] + e[5] * dfinv[(j << 2) - 1] * dfinv[(j << 2) - 1];
/* L6: */
		}
	    } else if (ibloc == 2) {
/*           --- BLOC  1,2 */
		for (j = 1; j <= 9; ++j) {
		    fdedf[(j << 2) - 4] = (e[1] + e[5]) * dfinv[(j << 2) - 3] * dfinv[(j << 2) - 4];
		    fdedf[(j << 2) - 3] = e[5] * dfinv[(j << 2) - 1] * dfinv[(j << 2) - 4] + e[1] * dfinv[(j << 2) - 3] * dfinv[(j << 2) - 2];
		    fdedf[(j << 2) - 2] = e[5] * dfinv[(j << 2) - 3] * dfinv[(j << 2) - 2] + e[1] * dfinv[(j << 2) - 4] * dfinv[(j << 2) - 1];
		    fdedf[(j << 2) - 1] = (e[1] + e[5]) * dfinv[(j << 2) - 1] * dfinv[(j << 2) - 2];
/* L7: */
		}
	    } else if (ibloc == 3) {
/*           --- BLOC  2,2 */
		for (j = 1; j <= 9; ++j) {
		    fdedf[(j << 2) - 4] = e[2] * dfinv[(j << 2) - 3] * dfinv[(j << 2) - 3] + e[5] * dfinv[(j << 2) - 4] * dfinv[(j << 2) - 4];
		    fdedf[(j << 2) - 3] = e[2] * dfinv[(j << 2) - 1] * dfinv[(j << 2) - 3] + e[5] * dfinv[(j << 2) - 2] * dfinv[(j << 2) - 4];
		    fdedf[(j << 2) - 2] = fdedf[(j << 2) - 3];
		    fdedf[(j << 2) - 1] = e[2] * dfinv[(j << 2) - 1] * dfinv[(j << 2) - 1] + e[5] * dfinv[(j << 2) - 2] * dfinv[(j << 2) - 2];
/* L8: */
		}
	    }
	} else {
/*         -- anisotrope */
	    if (ibloc == 1) {
/*           ---  BLOC  1,1        CAS   ANISOTROPE */
		i1 = 1;
		i2 = 4;
		i3 = 4;
		i4 = 6;
	    } else if (ibloc == 2) {
/*           --- BLOC  1,2 */
		i1 = 4;
		i2 = 6;
		i3 = 2;
		i4 = 5;
	    } else if (ibloc == 3) {
/*           --- BLOC  2,2 */
		i1 = 6;
		i2 = 5;
		i3 = 5;
		i4 = 3;
	    }

	    for (j = 1; j <= 9; ++j) {
		fdedf[(j << 2) - 4] = e[i1 - 1] * dfinv[(j << 2) - 4] * dfinv[(j << 2) - 4] + e[i4 - 1] * dfinv[(j << 2) - 3] * dfinv[(j << 2) - 3] + (e[i2 - 1] + e[i3 - 1]) * dfinv[(j << 2) - 3] * dfinv[(j << 2) - 4];
		fdedf[(j << 2) - 3] = e[i1 - 1] * dfinv[(j << 2) - 2] * dfinv[(j << 2) - 4] + e[i4 - 1] * dfinv[(j << 2) - 1] * dfinv[(j << 2) - 3] + e[i2 - 1] * dfinv[(j << 2) - 1] * dfinv[(j << 2) - 4] + e[i3 - 1] * dfinv[(j << 2) - 3] * dfinv[(j << 2) - 2];
		fdedf[(j << 2) - 2] = e[i1 - 1] * dfinv[(j << 2) - 2] * dfinv[(j << 2) - 4] + e[i4 - 1] * dfinv[(j << 2) - 1] * dfinv[(j << 2) - 3] + e[i3 - 1] * dfinv[(j << 2) - 1] * dfinv[(j << 2) - 4] + e[i2 - 1] * dfinv[(j << 2) - 3] * dfinv[(j << 2) - 2];
		fdedf[(j << 2) - 1] = e[i4 - 1] * dfinv[(j << 2) - 1] * dfinv[(j << 2) - 1] + e[i1 - 1] * dfinv[(j << 2) - 2] * dfinv[(j << 2) - 2] + (e[i2 - 1] + e[i3 - 1]) * dfinv[(j << 2) - 1] * dfinv[(j << 2) - 2];
/* L9: */
	    }
	}

/*       -- CALCUL DE AE:  BLOCS (IBLOC) et  PASSAGE AU PAR D.L. (IJT) */
	for (kin = 1; kin <= 9; ++kin) {
	    if (ibloc != 2) {
		taba2d_(&c__8, &c__2, &dp25[((kin << 3) + 1 << 1) - 18], &fdedf[(kin << 2) - 4], aux1);
	    } else {
		taba6d_(&c__8, &c__2, &dp25[((kin << 3) + 1 << 1) - 18], &fdedf[(kin << 2) - 4], aux1);
	    }

	    for (j = 1; j <= 8; ++j) {
		ifin = j;
		if (ibloc == 2) {
		    ifin = 8;
		}
		i__1 = ifin;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    i1 = iblt[ibloc - 1] + i__;
		    j1 = jblt[ibloc - 1] + j;
		    if (ijt[i1 - 1] <= ijt[j1 - 1]) {
			kk = ijt[j1 - 1] * (ijt[j1 - 1] - 1) / 2 + ijt[i1 - 1];
		    } else {
			kk = ijt[i1 - 1] * (ijt[i1 - 1] - 1) / 2 + ijt[j1 - 1];
		    }
		    ij = j * (j - 1) / 2 + i__;
		    if (ibloc == 2) {
			ij = (j - 1 << 3) + i__;
		    }
		    ae[kk] += aux1[ij - 1] * poids[kin - 1] / delta[kin - 1];
/* L10: */
		}
/* L11: */
	    }
/* L12: */
	}
/* L13: */
    }

} /* etr2q2c_ */

