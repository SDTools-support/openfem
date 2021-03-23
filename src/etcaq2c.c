/* etcaq2c.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__8 = 8;
static int32 c__9 = 9;

/* Subroutine */ int etcaq2c_(coor, car, iopt, u, alpha, theta, sigma)
doublereal *coor, *car;
int32 *iopt;
doublereal *u, *alpha, *theta, *sigma;
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

    static doublereal elas[10], dq1dr, dq2dr, dq3dr, dq4dr, dq5dr, dq6dr, dq7dr, dq8dr, dq1dz, dq2dz, dq3dz, dq4dz, dq5dz, dq6dz, dq7dz, dq8dz;
    static int32 i__;
    extern /* Subroutine */ int e2aq2c_();
    static doublereal r__, f1[9], f2[9], q1, q2, q3, q4, q5, q6, q7, q8, dfm1dp[144]	/* was [2][8][9] */;
    static int32 ip[12];
    static doublereal poidel[9];
    extern /* Subroutine */ int hookax_();
    static doublereal estrain_rr__, estrain_tt__, estrain_rz__, tstrain_rr__, estrain_zz__, tstrain_tt__, tstrain_zz__, alpha_r__, alpha_t__, alpha_z__, dfm1[36]	/* was [4][9] */;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* BUT: CALCUL DES CONTRAINTES DE L ELEMENT TRIA AP1D */
/* --- */
/* in : coor(noe,ndim) : coordonnees R(8), Z(8) des 8 noeuds */
/*      car        : caracteristiques des materiaux */
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
/*      U(ndim,noe): deplacements U_r et U_z aux 8 noeuds */
/*      alpha(3) : coef. dilatation thermique: radial, axial, tangentiel */
/*      theta(8) : theta aux 8 noeuds */
/* out: */
/*      SIGMA(4) :  S_rr, S_rz, S_zz, S_tt elastiques */
/* ..................................................................... */
/* 2Q25 -- XYNPI: coordonnees pt. int. numeriques (element reference) */
    /* Parameter adjustments */
    --sigma;
    --theta;
    --alpha;
    u -= 3;
    --car;
    coor -= 9;

    /* Function Body */
/*     -- POIDS: poids du schema d'integration numerique. */
/*     -- Valeurs des Polynomes de base aux pt. int. numerique. */
/*     -- Valeurs Derive'es des Poly. de base aux pt int. numerique. */

/*     -- CALCUL DE F1,F2,FFM1,DFM1DP */

    e2aq2c_(&c__8, &c__9, poids, p25, dp25, ip, f1, f2, dfm1dp, poidel, dfm1, &coor[9]);

    hookax_(iopt, &car[1], elas);


/*     Contriantes [SIG] = [Elas] * [STRAIN_elas - STRAIN_ther] + SIG_0 */


/*     1) --- STRAIN Elastiques */
/*                                                             [U_r(1)] */
/* [Strain_zz]   [  0     0   ...   0   dq1dz dq2dz ... dq8dz] [  .   ] */
/* [         ]   [                                           ] [  .   ] */
/* [Strain_rr]   [dq1dr dq2dr ... dq8dr   0      0  ...   0  ] [  .   ] */
/* [         ]   [                                           ] [U_r(8)] */
/* [Strain_tt] = [ p1    p2        p8                        ] [      ] */
/* [         ]   [---   ---  ...  ---     0      0  ...   0  ] [U_z(1)] */
/* [2*Stra_rz]   [ r     r         r                         ] [  .   ] */
/*               [                                           ] [  .   ] */
/*               [ dq1dz .......  dq8dz dq1dr ........  dq8dr] [  .   ] */
/*                                                             [U_z(8)] */

/*  barycentre F1(5) etr F2(5) */
    i__ = 5;
    r__ = f1[i__ - 1];
    q1 = p25[(i__ << 3) - 8];
    q2 = p25[(i__ << 3) - 7];
    q3 = p25[(i__ << 3) - 6];
    q4 = p25[(i__ << 3) - 5];
    q5 = p25[(i__ << 3) - 4];
    q6 = p25[(i__ << 3) - 3];
    q7 = p25[(i__ << 3) - 2];
    q8 = p25[(i__ << 3) - 1];
    dq1dr = dfm1dp[((i__ << 3) + 1 << 1) - 18];
    dq2dr = dfm1dp[((i__ << 3) + 2 << 1) - 18];
    dq3dr = dfm1dp[((i__ << 3) + 3 << 1) - 18];
    dq4dr = dfm1dp[((i__ << 3) + 4 << 1) - 18];
    dq5dr = dfm1dp[((i__ << 3) + 5 << 1) - 18];
    dq6dr = dfm1dp[((i__ << 3) + 6 << 1) - 18];
    dq7dr = dfm1dp[((i__ << 3) + 7 << 1) - 18];
    dq8dr = dfm1dp[((i__ << 3) + 8 << 1) - 18];
    dq1dz = dfm1dp[((i__ << 3) + 1 << 1) - 17];
    dq2dz = dfm1dp[((i__ << 3) + 2 << 1) - 17];
    dq3dz = dfm1dp[((i__ << 3) + 3 << 1) - 17];
    dq4dz = dfm1dp[((i__ << 3) + 4 << 1) - 17];
    dq5dz = dfm1dp[((i__ << 3) + 5 << 1) - 17];
    dq6dz = dfm1dp[((i__ << 3) + 6 << 1) - 17];
    dq7dz = dfm1dp[((i__ << 3) + 7 << 1) - 17];
    dq8dz = dfm1dp[((i__ << 3) + 8 << 1) - 17];
    estrain_zz__ = dq1dz * u[4] + dq2dz * u[6] + dq3dz * u[8] + dq4dz * u[10] + dq5dz * u[12] + dq6dz * u[14] + dq7dz * u[16] + dq8dz * u[18];
    estrain_rr__ = dq1dr * u[3] + dq2dr * u[5] + dq3dr * u[7] + dq4dr * u[9] + dq5dr * u[11] + dq6dr * u[13] + dq7dr * u[15] + dq8dr * u[17];
    estrain_tt__ = (q1 * u[3] + q2 * u[5] + q3 * u[7] + q4 * u[9] + q5 * u[11] + q6 * u[13] + q7 * u[15] + q8 * u[17]) / r__;
    estrain_rz__ = (dq1dz * u[3] + dq1dr * u[4] + dq2dz * u[5] + dq2dr * u[6] + dq3dz * u[7] + dq3dr * u[8] + dq4dz * u[9] + dq4dr * u[10] + dq5dz * u[11] + dq5dr * u[12] + dq6dz * u[13] + dq6dr * u[14] + dq7dz * u[15] + dq7dr * u[16] + dq8dz * u[17] + dq8dr * u[18]) * (float).5;

/*     2) --- STRAIN THermiques =  (ALPHA) * (P) * THETA */

/* [ Strain_zz ] [alpha_z]                            [theta(1)] */
/* [           ] [       ]                            [theta(2)] */
/* [ Strain_rr ] [alpha_r]                            [theta(3)] */
/* [           ]=[       ]*[q1 q2 q3  q4 q5 q6 q7 q8]*[theta(4)] */
/* [ Strain_tt ] [alpha_t]                            [theta(5)] */
/* [           ] [       ]                            [theta(6)] */
/* [2*Strain_rz] [  0    ]                            [theta(7)] */
/*                                                    [theta(8)] */
    alpha_r__ = alpha[1];
    alpha_z__ = alpha[2];
    alpha_t__ = alpha[3];
    q1 = p25[(i__ << 3) - 8];
    q2 = p25[(i__ << 3) - 7];
    q3 = p25[(i__ << 3) - 6];
    q4 = p25[(i__ << 3) - 5];
    q5 = p25[(i__ << 3) - 4];
    q6 = p25[(i__ << 3) - 3];
    q7 = p25[(i__ << 3) - 2];
    q8 = p25[(i__ << 3) - 1];
    tstrain_zz__ = alpha_z__ * (q1 * theta[1] + q2 * theta[2] + q3 * theta[3] + q4 * theta[4] + q5 * theta[5] + q6 * theta[6] + q7 * theta[7] + q8 * theta[8]);
    tstrain_rr__ = alpha_r__ * (q1 * theta[1] + q2 * theta[2] + q3 * theta[3] + q4 * theta[4] + q5 * theta[5] + q6 * theta[6] + q7 * theta[7] + q8 * theta[8]);
    tstrain_tt__ = alpha_t__ * (q1 * theta[1] + q2 * theta[2] + q3 * theta[3] + q4 * theta[4] + q5 * theta[5] + q6 * theta[6] + q7 * theta[7] + q8 * theta[8]);

/*     3) --- Contraintes */

/*  [  Sig_zz  ]    [ elas(1)  elas(2)  elas(4)   elas(7) ] [  Strain_zz  ] */
/*  [          ]    [                                     ] [             ] */
/*  [  Sig_rr  ]    [ elas(2)  elas(3)  elas(5)   elas(8) ] [  Strain_rr  ] */
/*  [          ] := [                                     ] [             ] */
/*  [  Sig_tt  ]    [ elas(4)  elas(5)  elas(6)   elas(9) ] [  Strain_tt  ] */
/*  [          ]    [                                     ] [             ] */
/*  [  Sig_rz  ]    [ elas(7)  elas(8)  elas(9)  elas(10) ] [ 2*Strain_rz ] */

/*    Rangement SIGMA(4): S_rr, S_rz, S_zz, S_tt */

    sigma[3] = elas[0] * (estrain_zz__ - tstrain_zz__) + elas[1] * (estrain_rr__ - tstrain_rr__) + elas[3] * (estrain_tt__ - tstrain_tt__) + elas[6] * 2 * estrain_rz__;
    sigma[1] = elas[1] * (estrain_zz__ - tstrain_zz__) + elas[2] * (estrain_rr__ - tstrain_rr__) + elas[4] * (estrain_tt__ - tstrain_tt__) + elas[7] * 2 * estrain_rz__;
    sigma[4] = elas[3] * (estrain_zz__ - tstrain_zz__) + elas[4] * (estrain_rr__ - tstrain_rr__) + elas[5] * (estrain_tt__ - tstrain_tt__) + elas[8] * 2 * estrain_rz__;
    sigma[2] = elas[6] * (estrain_zz__ - tstrain_zz__) + elas[7] * (estrain_rr__ - tstrain_rr__) + elas[8] * (estrain_tt__ - tstrain_tt__) + elas[9] * 2 * estrain_rz__;
    return 0;
} /* etcaq2c_ */

