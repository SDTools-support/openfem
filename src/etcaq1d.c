/* etcaq1d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__4 = 4;

/* Subroutine */ int etcaq1d_(coor, car, iopt, u, alpha, theta, sigma)
doublereal *coor, *car;
int32 *iopt;
doublereal *u, *alpha, *theta, *sigma;
{
    /* Initialized data */

    static doublereal poids[4] = { .25,.25,.25,.25 };
    static doublereal q13[16]	/* was [4][4] */ = { .62200846792814612,.16666666666666662,.044658198738520435,.16666666666666662,.16666666666666662,.62200846792814623,.16666666666666662,.044658198738520449,.044658198738520504,.16666666666666662,.62200846792814623,.16666666666666662,.16666666666666668,.044658198738520449,.16666666666666662,.62200846792814623 };
    static doublereal dq13[32]	/* was [2][4][4] */ = { -.78867513459481286,-.78867513459481286,.78867513459481286,-.21132486540518707,.21132486540518707,.21132486540518707,-.21132486540518707,.78867513459481286,-.78867513459481286,-.21132486540518713,.78867513459481286,-.78867513459481286,.21132486540518707,.78867513459481286,-.21132486540518707,.21132486540518713,-.21132486540518713,-.21132486540518713,.21132486540518713,-.78867513459481286,.78867513459481286,.78867513459481286,
	    -.78867513459481286,.21132486540518713,-.21132486540518713,-.78867513459481286,.21132486540518713,-.21132486540518707,.78867513459481286,.21132486540518707,-.78867513459481286,.78867513459481286 };

    static doublereal elas[10], dq1dr, dq2dr, dq3dr, dq4dr, dq1dz, dq2dz, dq3dz, dq4dz;
    static int32 i__;
    extern /* Subroutine */ int e1aq1c_();
    static doublereal r__, f1[5], f2[5], q1, q2, q3, q4, dfm1dp[32]	/* was [2][4][4] */;
    static int32 ip[8];
    static doublereal poidel[4];
    extern /* Subroutine */ int hookax_();
    static doublereal estrain_rr__, estrain_tt__, estrain_rz__, tstrain_rr__, estrain_zz__, tstrain_tt__, tstrain_zz__, alpha_r__, alpha_t__, alpha_z__, dfm1[16]	/* was [4][4] */;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* BUT: CALCUL DES CONTRAINTES DE L ELEMENT TRIA AP1D */
/* --- */
/* in : coor(noe,ndim) : coordonnees R(3), Z(3) des 3 sommets. */

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
/*      U(ndim,noe): deplacements U_r et U_z aux 3 noeuds (sommets) */
/*      alpha(3) : coef. dilatation thermique: radial, axial, tangentiel */
/*      theta(4)    : theta aux 4 sommets */
/* out: */
/*      SIGMA(4) :  S_rr, S_rz, S_zz, S_tt elastiques */
/* ..................................................................... */
/* 2Q13 -- POIDS: poids du schema d'integration numerique. */
    /* Parameter adjustments */
    --sigma;
    --theta;
    --alpha;
    u -= 3;
    --car;
    coor -= 5;

    /* Function Body */
/*     -- XYNPI: coordonnees pt. int. numeriques (element reference) */
/*     -- Valeurs des Polynomes de base aux pt. int. numerique. */
/*     -- Valeurs Derive'es des Poly. de base aux pt int. numerique. */

/*     CALCUL DE F1,F2,DFM1,DFM1DP */
/*     ----------------------- */
    e1aq1c_(&c__4, &c__4, poids, q13, dq13, ip, f1, f2, dfm1dp, poidel, dfm1, &coor[5]);

    hookax_(iopt, &car[1], elas);

/*     1) --- STRAIN Elastiques                                         [ U_r(1) ] */
/*                                                                      [ U_r(2) ] */
/*  [ Estrain_zz ]    [  0     0     0     0   dq1dz dq2dz dq3dz dq4dz] [        ] */
/*  [            ]    [                                               ] [ U_r(3) ] */
/*  [ Estrain_rr ]    [dq1dr dq2dr dq3dr dq4dr   0      0      0   0  ] [        ] */
/*  [            ]    [                                               ] [ U_r(4) ] */
/*  [ Estrain_tt ] := [  p1    p2    p3    p4                         ] [        ] */
/*  [            ]    [ ----  ----  ----  ----   0      0      0   0  ] [ U_z(1) ] */
/*  [2*Estrain_rz]    [   r     r     r     r                         ] [        ] */
/*                    [                                               ] [ U_z(2) ] */
/*                    [dq1dz dq2dz dq3dz dq4dz dq1dr dq2dr dq3dr dq4dr] [        ] */
/*                                                                      [ U_z(3) ] */
/*                                                                      [ U_z(4) ] */
    sigma[1] = 0.;
    sigma[2] = 0.;
    sigma[3] = 0.;
    sigma[4] = 0.;
    for (i__ = 1; i__ <= 4; ++i__) {
	q1 = q13[(i__ << 2) - 4];
	q2 = q13[(i__ << 2) - 3];
	q3 = q13[(i__ << 2) - 2];
	q4 = q13[(i__ << 2) - 1];
	dq1dr = dfm1dp[((i__ << 2) + 1 << 1) - 10];
	dq2dr = dfm1dp[((i__ << 2) + 2 << 1) - 10];
	dq3dr = dfm1dp[((i__ << 2) + 3 << 1) - 10];
	dq4dr = dfm1dp[((i__ << 2) + 4 << 1) - 10];
	dq1dz = dfm1dp[((i__ << 2) + 1 << 1) - 9];
	dq2dz = dfm1dp[((i__ << 2) + 2 << 1) - 9];
	dq3dz = dfm1dp[((i__ << 2) + 3 << 1) - 9];
	dq4dz = dfm1dp[((i__ << 2) + 4 << 1) - 9];
	r__ = f1[i__ - 1];
	estrain_zz__ = dq1dz * u[4] + dq2dz * u[6] + dq3dz * u[8] + dq4dz * u[10];
	estrain_rr__ = dq1dr * u[3] + dq2dr * u[5] + dq3dr * u[7] + dq4dr * u[9];
	estrain_tt__ = (q1 * u[3] + q2 * u[5] + q3 * u[7] + q4 * u[9]) / r__;
	estrain_rz__ = (dq1dz * u[3] + dq2dz * u[5] + dq3dz * u[7] + dq4dz * u[9] + dq1dr * u[4] + dq2dr * u[6] + dq3dr * u[8] + dq4dr * u[10]) * (float).5;
/*     2) --- Thermal STRAIN */

/*  [ Tstrain_zz ]    [ alpha[z] ]                    [ theta[1] ] */
/*  [            ]    [          ]                    [          ] */
/*  [ Tstrain_rr ]    [ alpha[r] ]                    [ theta[2] ] */
/*  [            ] := [          ] [ q1  q2  q3  q4 ] [          ] */
/*  [ Tstrain_tt ]    [ alpha[t] ]                    [ theta[3] ] */
/*  [            ]    [          ]                    [          ] */
/*  [2*Tstrain_rz]    [   0      ]                    [ theta[4] ] */


	alpha_r__ = alpha[1];
	alpha_z__ = alpha[2];
	alpha_t__ = alpha[3];
	tstrain_zz__ = alpha_z__ * (q1 * theta[1] + q2 * theta[2] + q3 * theta[3] + q4 * theta[4]);
	tstrain_rr__ = alpha_r__ * (q1 * theta[1] + q2 * theta[2] + q3 * theta[3] + q4 * theta[4]);
	tstrain_tt__ = alpha_t__ * (q1 * theta[1] + q2 * theta[2] + q3 * theta[3] + q4 * theta[4]);

/*     3) --- STRESS */
/*                                            [                        ] */
/*                                            [ Estrain_zz - Tstrain_zz] */
/* [ sigma_zz ]  [ elas1 elas2 elas4  elas7 ] [                        ] */
/* [          ]  [                          ] [                        ] */
/* [ sigma_rr ]  [ elas2 elas3 elas5  elas8 ] [ Estrain_rr - Tstrain_rr] */
/* [          ]= [                          ] [                        ] */
/* [ sigma_tt ]  [ elas4 elas5 elas6  elas9 ] [                        ] */
/* [          ]  [                          ] [ Estrain_tt - Tstrain_tt] */
/* [ sigma_rz ]  [ elas7 elas8 elas9 elas10 ] [                        ] */
/*                                            [                        ] */
/*                                            [2Estrain_rz             ] */
/*    Rangement SIGMA(4): S_rr, S_rz, S_zz, S_tt */
	sigma[3] = sigma[3] + elas[0] * (estrain_zz__ - tstrain_zz__) + elas[1] * (estrain_rr__ - tstrain_rr__) + elas[3] * (estrain_tt__ - tstrain_tt__) + elas[6] * 2 * estrain_rz__;
	sigma[1] = sigma[1] + elas[1] * (estrain_zz__ - tstrain_zz__) + elas[2] * (estrain_rr__ - tstrain_rr__) + elas[4] * (estrain_tt__ - tstrain_tt__) + elas[7] * 2 * estrain_rz__;
	sigma[4] = sigma[4] + elas[3] * (estrain_zz__ - tstrain_zz__) + elas[4] * (estrain_rr__ - tstrain_rr__) + elas[5] * (estrain_tt__ - tstrain_tt__) + elas[8] * 2 * estrain_rz__;
	sigma[2] = sigma[2] + elas[6] * (estrain_zz__ - tstrain_zz__) + elas[7] * (estrain_rr__ - tstrain_rr__) + elas[8] * (estrain_tt__ - tstrain_tt__) + elas[9] * 2 * estrain_rz__;

/* L1: */
    }
    sigma[1] /= 4;
    sigma[2] /= 4;
    sigma[3] /= 4;
    sigma[4] /= 4;

} /* etcaq1d_ */

