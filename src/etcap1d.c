/* etcap1d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__3 = 3;
static int32 c__1 = 1;

/* Subroutine */ int etcap1d_(coor, car, iopt, u, alpha, theta, sigma)
doublereal *coor, *car;
int32 *iopt;
doublereal *u, *alpha, *theta, *sigma;
{
    /* Initialized data */

    static doublereal poids[1] = { .5 };
    static doublereal p13[3]	/* was [3][1] */ = { .33333333333333333,.33333333333333333,.33333333333333333 };
    static doublereal dp13[6]	/* was [2][3][1] */ = { -1.,-1.,1.,0.,0.,1. };

    static doublereal elas[10], dp1dr, dp2dr, dp3dr, dp1dz, dp2dz, dp3dz;
    static int32 i__;
    extern /* Subroutine */ int e1ap1d_();
    static doublereal r__, f1[1], f2[1], p1, p2, p3, dfm1dp[18]	/* was [2][3][3] */;
    static int32 ip[6];
    static doublereal poidel[1];
    extern /* Subroutine */ int hookax_();
    static doublereal estrain_rr__, estrain_tt__, estrain_rz__, tstrain_rr__, estrain_zz__, tstrain_tt__, tstrain_zz__, alpha_r__, alpha_t__, alpha_z__, dfm1[16]	/* was [4][4] */;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* BUT: CALCUL DES CONTRAINTES DE L ELEMENT TRIA AP1D */
/* --- */
/* in : coor(noe,ndim) : coordonnees R(3), Z(3) des 3 sommets */
/*      car, iopt      : caracteristiques des materiaux */
/*      U(ndim,noe): deplacements U_r et U_z aux 3 noeuds (sommets) */
/*      alpha(3) : coef. dilatation thermique: radial, axial, tangentiel */
/*      theta(3)    : theta aux 3 sommets */
/* out: */
/*      SIGMA(4) :  S_rr, S_rz, S_zz, S_tt elastiques */
/* ..................................................................... */
/*     -- XYNPI: coordonnees pt. int. numeriques (element reference) */
    /* Parameter adjustments */
    --sigma;
    --theta;
    --alpha;
    u -= 3;
    --car;
    coor -= 4;

    /* Function Body */
/*     -- POIDS: poids du schema d'integration numerique. */
/*     -- Valeurs des Polynomes de base aux pt. int. numerique. */
/*     -- Valeurs Derive'es des Poly. de base aux pt int. numerique. */

/*     CALCUL DE F1,F2,DFM1,DFM1DP */
/*     ----------------------- */
    e1ap1d_(&c__3, &c__1, poids, p13, dp13, ip, f1, f2, dfm1dp, poidel, dfm1, &coor[4]);

    hookax_(iopt, &car[1], elas);

/*     Contriantes [SIG] = [Elas] * [STRAIN_elas - STRAIN_ther] + SIG_0 */

/*     1) --- Deformations Elastiques */
/*                                                                  [ U_r(1) ] */
/*  [  Strain_zz  ]    [   0      0      0    dp1dz  dp2dz  dp3dz ] [        ] */
/*  [             ]    [                                          ] [ U_r(2) ] */
/*  [  Strain_rr  ]    [ dp1dr  dp2dr  dp3dr    0      0      0   ] [        ] */
/*  [             ]    [                                          ] [ U_r(3) ] */
/*  [  Strain_tt  ] := [   p1     p2     p3                       ] [        ] */
/*  [             ]    [  ----   ----   ----    0      0      0   ] [ U_z(1) ] */
/*  [ 2*Strain_rz ]    [    r      r      r                       ] [        ] */
/*                     [                                          ] [ U_z(2) ] */
/*                     [ dp1dz  dp2dz  dp3dz  dp1dr  dp2dr  dp3dr ] [        ] */
/*                                                                  [ U_z(3) ] */

    i__ = 1;
    r__ = f1[i__ - 1];
    p1 = p13[i__ * 3 - 3];
    p2 = p13[i__ * 3 - 2];
    p3 = p13[i__ * 3 - 1];
    dp1dr = dfm1dp[(i__ * 3 + 1 << 1) - 8];
    dp2dr = dfm1dp[(i__ * 3 + 2 << 1) - 8];
    dp3dr = dfm1dp[(i__ * 3 + 3 << 1) - 8];
    dp1dz = dfm1dp[(i__ * 3 + 1 << 1) - 7];
    dp2dz = dfm1dp[(i__ * 3 + 2 << 1) - 7];
    dp3dz = dfm1dp[(i__ * 3 + 3 << 1) - 7];
    estrain_zz__ = dp1dz * u[4] + dp2dz * u[6] + dp3dz * u[8];
    estrain_rr__ = dp1dr * u[3] + dp2dr * u[5] + dp3dr * u[7];
    estrain_tt__ = (p1 * u[3] + p2 * u[5] + p3 * u[7]) / r__;
    estrain_rz__ = (dp1dz * u[3] + dp2dz * u[5] + dp3dz * u[7] + dp1dr * u[4] + dp2dr * u[6] + dp3dr * u[8]) * (float).5;

/*     2) --- Deformations THermiques =  (ALPHA) * (P) * THETA */

/*  [  Strain_zz ]  [alpha_z] */
/*  [            ]  [       ] */
/*  [  Strain_rr ]  [alpha_r]                  [ theta(1) ] */
/*  [            ]:=[       ] * [ p1 p2 p3 ] * [ theta(2) ] */
/*  [  Strain_tt ]  [alpha_t]                  [ theta(3) ] */
/*  [            ]  [       ] */
/*  [ 2*Strain_rz]  [  0    ] */

    alpha_r__ = alpha[1];
    alpha_z__ = alpha[2];
    alpha_t__ = alpha[3];
    tstrain_zz__ = alpha_z__ * (p1 * theta[1] + p2 * theta[2] + p3 * theta[3]);
    tstrain_rr__ = alpha_r__ * (p1 * theta[1] + p2 * theta[2] + p3 * theta[3]);
    tstrain_tt__ = alpha_t__ * (p1 * theta[1] + p2 * theta[2] + p3 * theta[3]);
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
} /* etcap1d_ */

