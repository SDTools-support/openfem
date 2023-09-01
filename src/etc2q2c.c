/* etc2q2c.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int etc2q2c_(coor, car, iopt, u, sigma)
doublereal *coor, *car;
int32 *iopt;
doublereal *u, *sigma;
{
    /* Initialized data */

    static int32 ijt[16] = { 1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16 };
    static doublereal dp25[144]	/* was [2][8][9] */ = { -2.061895003862225,-2.061895003862225,-.68729833462074174,-.087298334620741685,-.26189500386222502,-.26189500386222502,-.087298334620741685,-.68729833462074174,2.7491933384829669,-.39999999999999996,.39999999999999996,.34919333848296674,.34919333848296674,.39999999999999996,-.39999999999999996,2.7491933384829669,-.68729833462074152,-.7745966692414834,.68729833462074174,-.7745966692414834,-.087298334620741657,-.7745966692414834,
	    .087298334620741685,-.7745966692414834,0.,-1.,.39999999999999996,1.5491933384829668,0.,1.,-.39999999999999996,1.5491933384829668,.68729833462074196,-.087298334620742101,2.0618950038622254,-2.061895003862225,.087298334620741713,-.68729833462074174,.26189500386222508,-.26189500386222497,-2.7491933384829669,-.39999999999999991,.39999999999999996,2.7491933384829669,-.34919333848296674,.39999999999999991,-.39999999999999996,.34919333848296674,-.7745966692414834,-.68729833462074174,
	    -.7745966692414834,.087298334620741685,-.7745966692414834,-.087298334620741685,-.7745966692414834,.68729833462074174,1.5491933384829668,-.39999999999999996,1.,0.,1.5491933384829668,.39999999999999996,-1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-1.,1.,0.,0.,1.,-1.,0.,.7745966692414834,.087298334620741213,.7745966692414834,-.68729833462074174,.7745966692414834,.68729833462074174,.7745966692414834,-.087298334620741657,-1.5491933384829668,-.39999999999999991,1.,0.,-1.5491933384829668,
	    .39999999999999991,-1.,0.,-.087298334620742157,.68729833462074174,-.26189500386222502,.26189500386222502,-.68729833462074174,.087298334620741685,-2.0618950038622254,2.061895003862225,.34919333848296674,-.39999999999999996,.39999999999999991,-.34919333848296674,2.7491933384829669,.39999999999999996,-.39999999999999991,-2.7491933384829669,.087298334620741213,.7745966692414834,-.087298334620741657,.7745966692414834,.68729833462074174,.7745966692414834,-.68729833462074196,
	    .7745966692414834,0.,-1.,.39999999999999991,-1.5491933384829668,0.,1.,-.39999999999999991,-1.5491933384829668,.26189500386222475,.26189500386222452,.087298334620741879,.68729833462074174,2.061895003862225,2.061895003862225,.68729833462074152,.087298334620741657,-.34919333848296663,-.39999999999999991,.39999999999999991,-2.7491933384829669,-2.7491933384829669,.39999999999999991,-.39999999999999991,-.34919333848296663 };

    static doublereal edfp[24]	/* was [3][8] */, delt5, c__, e[6];
    static int32 i__, j, l;
    static doublereal s, dfidp[16]	/* was [2][8] */;
    static int32 ibloc;
    static doublereal dfinv[4]	/* was [2][2] */;
    static int32 i1, i2, i3, i4, i5, i6, j1;
    static doublereal young, unmnu, f11[9], f12[9], f21[9], f22[9], ed[6];
    static int32 kk;
    static doublereal dsigma[48], poisson;

/* *************************************************************** */
/* BUT: CALCUL DES CONTRAINTES DE L ELEMENT QUAD 2Q2C */
/* --- */
/* in : coor(noe,ndim) : coordonnees des 8 noeuds. */
/*      car, iopt      : caracteristiques des materiaux */
/*      U(ndim,noe): deplacements U_x et U_y aux 8 noeuds */
/* out: SIGMA(3)   :  S_xx, S_yy, S_xy elastiques */
/* programmeur : modulef */
/* ............................................................... */
    /* Parameter adjustments */
    --sigma;
    u -= 9;
    --car;
    coor -= 9;

    /* Function Body */
/* 2Q25 -- XYNPI: coordonnees pt. int. numeriques (element reference) */
/*     -- POIDS: poids du schema d'integration numerique. */
/*     -- Valeurs des Polynomes de base aux pt. int. numerique. */
/*     -- Valeurs Derive'es des Poly. de base aux pt int. numerique. */


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

/* -- CALCUL DES CONTRAINTES ELEMENTAIRES AU BARYCENTRE */
/*            CALCUL DE DELT5 = DELTA(5) */

    f11[4] = 0.;
    f12[4] = 0.;
    f21[4] = 0.;
    f22[4] = 0.;
    for (i__ = 1; i__ <= 8; ++i__) {
	f11[4] += dp25[(i__ + 40 << 1) - 18] * coor[i__ + 8];
	f12[4] += dp25[(i__ + 40 << 1) - 17] * coor[i__ + 8];
	f21[4] += dp25[(i__ + 40 << 1) - 18] * coor[i__ + 16];
	f22[4] += dp25[(i__ + 40 << 1) - 17] * coor[i__ + 16];
/* L2: */
    }
    delt5 = f11[4] * f22[4] - f12[4] * f21[4];

/*     ----  CALCUL DE DFINV  L INVERSE DE DF  ----- */

    dfinv[0] = f22[4];
    dfinv[1] = -f12[4];
    dfinv[2] = -f21[4];
    dfinv[3] = f11[4];

/*     ----  CALCUL DE DFINV*DP */

    for (i__ = 1; i__ <= 2; ++i__) {
	for (j = 1; j <= 8; ++j) {
	    s = 0.;
	    for (l = 1; l <= 2; ++l) {
		s += dfinv[i__ + (l << 1) - 3] * dp25[l + (j + 40 << 1) - 19];
/* L3: */
	    }
	    dfidp[i__ + (j << 1) - 3] = s / delt5;
/* L4: */
	}
    }

/*     --CALCUL DE E*D */
    for (ibloc = 1; ibloc <= 2; ++ibloc) {
	if (ibloc == 1) {
	    i1 = 1;
	    i2 = 2;
	    i3 = 4;
	    i4 = 4;
	    i5 = 5;
	    i6 = 6;
	} else if (ibloc == 2) {
	    i1 = 4;
	    i2 = 5;
	    i3 = 6;
	    i4 = 2;
	    i5 = 3;
	    i6 = 5;
	}
	ed[0] = e[i1 - 1];
	ed[1] = e[i2 - 1];
	ed[2] = e[i3 - 1];
	ed[3] = e[i4 - 1];
	ed[4] = e[i5 - 1];
	ed[5] = e[i6 - 1];

/*       -- CALCUL DE ED * DFIDP PAR BLOC */

	for (j = 1; j <= 8; ++j) {
	    for (i__ = 1; i__ <= 3; ++i__) {
		s = 0.;
		for (l = 1; l <= 2; ++l) {
		    kk = (l - 1) * 3 + i__;
		    s += ed[kk - 1] * dfidp[l + (j << 1) - 3];
/* L5: */
		}
		edfp[i__ + j * 3 - 4] = s;
/* L6: */
	    }
	}

/*       -- CALCUL DE DSIGMA(48) */

	for (j = 1; j <= 8; ++j) {
	    for (i__ = 1; i__ <= 3; ++i__) {
		j1 = j;
		if (ibloc == 2) {
		    j1 = j + 8;
		}
		kk = (ijt[j1 - 1] - 1) * 3 + i__;
		dsigma[kk - 1] = edfp[i__ + j * 3 - 4];
/* L7: */
	    }
	}
/* L8: */
    }


/*  [  11] [1 4 7 10 13 16 19 22 25 28 31 34 37 40 43 46] [     ] */
/*  [S 22]=[2 5 8 11 14 17 20 23 26 29 32 35 38 41 44 47]*[u_sol] */
/*  [  12] [3 6 9 12 15 18 21 24 27 30 33 36 39 42 45 48] [     ] */
/*      3*1                                          3*48    16*1 */
    sigma[1] = dsigma[0] * u[9] + dsigma[3] * u[17] + dsigma[6] * u[10] + dsigma[9] * u[18] + dsigma[12] * u[11] + dsigma[15] * u[19] + dsigma[18] * u[12] + dsigma[21] * u[20] + dsigma[24] * u[13] + dsigma[27] * u[21] + dsigma[30] * u[14] + dsigma[33] * u[22] + dsigma[36] * u[15] + dsigma[39] * u[23] + dsigma[42] * u[16] + dsigma[45] * u[24];
    sigma[2] = dsigma[1] * u[9] + dsigma[4] * u[17] + dsigma[7] * u[10] + dsigma[10] * u[18] + dsigma[13] * u[11] + dsigma[16] * u[19] + dsigma[19] * u[12] + dsigma[22] * u[20] + dsigma[25] * u[13] + dsigma[28] * u[21] + dsigma[31] * u[14] + dsigma[34] * u[22] + dsigma[37] * u[15] + dsigma[40] * u[23] + dsigma[43] * u[16] + dsigma[46] * u[24];
    sigma[3] = dsigma[2] * u[9] + dsigma[5] * u[17] + dsigma[8] * u[10] + dsigma[11] * u[18] + dsigma[14] * u[11] + dsigma[17] * u[19] + dsigma[20] * u[12] + dsigma[23] * u[20] + dsigma[26] * u[13] + dsigma[29] * u[21] + dsigma[32] * u[14] + dsigma[35] * u[22] + dsigma[38] * u[15] + dsigma[41] * u[23] + dsigma[44] * u[16] + dsigma[47] * u[24];

/*      print *,'---------------- Verif avec impressions Modulef' */
/*      print *, (dsigma(i), i= 1, 6) */
/*      print *, (dsigma(i), i= 7, 12) */
/*      print *, (dsigma(i), i= 13, 18) */
/*      print *, (dsigma(i), i= 19, 24) */
/*      print *, (dsigma(i), i= 25, 30) */
/*      print *, (dsigma(i), i= 31, 36) */
/*      print *, (dsigma(i), i= 37, 42) */
/*      print *, (dsigma(i), i= 43, 48) */

} /* etc2q2c_ */

