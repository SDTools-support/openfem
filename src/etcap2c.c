/* etcap2c.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__6 = 6;
static int32 c__7 = 7;

/* Subroutine */ int etcap2c_(coor, car, iopt, u, alpha, theta, sigma)
doublereal *coor, *car;
int32 *iopt;
doublereal *u, *alpha, *theta, *sigma;
{
    /* Initialized data */

    static doublereal poids[7] = { .062969590272413583,.062969590272413583,.062969590272413583,.066197076394253095,.066197076394253095,.066197076394253095,.11249999701976776 };
    static doublereal p25[42]	/* was [6][7] */ = { .47435260858553857,-.080768594191887185,-.080768594191887185,.32307437676754874,.041035826263138293,.32307437676754874,-.080768594191886977,.4743526085855384,-.080768594191887185,.32307437676754879,.32307437676754868,.041035826263138341,-.080768594191887033,-.080768594191887185,.4743526085855384,.041035826263138334,.32307437676754868,.3230743767675489,-.028074943223078796,-.028074943223078852,-.052583901102545349,.88413424176407262,
	    .11229977289231514,.11229977289231517,-.052583901102545571,-.028074943223078852,-.028074943223078852,.1122997728923154,.8841342417640724,.1122997728923154,-.028074943223078765,-.052583901102545349,-.028074943223078852,.11229977289231514,.11229977289231514,.88413424176407262,-.11111111773384862,-.11111110779974175,-.11111110779974175,.44444443119896703,.44444447093539807,.44444443119896703 };
    static doublereal dp25[84]	/* was [2][6][7] */ = { -2.1897079414123492,-2.1897079414123492,-.59485397070617462,0.,0.,-.59485397070617462,2.7845619121185238,-.40514602929382531,.40514602929382531,.40514602929382531,-.40514602929382531,2.7845619121185238,.59485397070617418,.59485397070617418,2.1897079414123488,0.,0.,-.59485397070617462,-2.7845619121185229,-3.1897079414123488,.40514602929382531,3.1897079414123488,-.40514602929382531,0.,.59485397070617418,.59485397070617418,-.59485397070617462,
	    0.,0.,2.1897079414123488,0.,-.40514602929382531,3.1897079414123488,.40514602929382531,-3.1897079414123488,-2.7845619121185229,-.88056825642046065,-.88056825642046065,.88056825642046021,0.,0.,-.76113651284092076,0.,-1.8805682564204602,.2388634871590792,1.8805682564204602,-.2388634871590792,1.6417047692613815,.76113651284092043,.76113651284092043,.88056825642046021,0.,0.,.88056825642046021,-1.6417047692613806,-1.8805682564204602,1.8805682564204602,1.8805682564204602,-1.8805682564204602,
	    -1.6417047692613806,-.88056825642046054,-.88056825642046054,-.76113651284092076,0.,0.,.88056825642046021,1.6417047692613813,-.2388634871590792,1.8805682564204602,.2388634871590792,-1.8805682564204602,0.,-.33333325386047363,-.33333325386047363,.33333337306976318,0.,0.,.33333337306976318,-1.1920928955078125e-7,-1.3333333730697631,1.3333333730697631,1.3333333730697631,-1.3333333730697631,-1.1920928955078125e-7 };

    static doublereal elas[10], dp1dr, dp2dr, dp3dr, dp4dr, dp5dr, dp6dr, dp1dz, dp2dz, dp3dz, dp4dz, dp5dz, dp6dz;
    static int32 i__;
    extern /* Subroutine */ int e2ap2c_();
    static doublereal r__, f1[7], f2[7], p1, p2, p3, p4, p5, p6, dfm1dp[84]	/* was [2][6][7] */;
    static int32 ip[12];
    static doublereal poidel[7];
    extern /* Subroutine */ int hookax_();
    static doublereal estrain_rr__, estrain_tt__, estrain_rz__, tstrain_rr__, estrain_zz__, tstrain_tt__, tstrain_zz__, alpha_r__, alpha_t__, alpha_z__, dfm1[28]	/* was [4][7] */;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* BUT: CALCUL DES CONTRAINTES DE L ELEMENT TRIA AP2C */
/* --- */
/* in : coor(noe,ndim) : coordonnees R(6), Z(6) des 6 noeuds. */
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
/*      U(ndim,noe): deplacements U_r et U_z aux 6 noeuds */
/*      alpha(3) : coef. dilatation thermique: radial, axial, tangentiel */
/*      theta(6)   : theta aux 6 noeuds */
/* out: */
/*      SIGMA(4)   :  S_rr, S_rz, S_zz, S_tt elastiques */
/* ..................................................................... */
/* 2P25 -- XYNPI: coordonnees pt. int. numeriques (element reference) */
    /* Parameter adjustments */
    --sigma;
    --theta;
    --alpha;
    u -= 3;
    --car;
    coor -= 7;

    /* Function Body */
/*     -- POIDS: poids du schema d'integration numerique. */
/*     -- Valeurs des Polynomes de base aux pt. int. numerique. */
/*     -- Valeurs Derive'es des Poly. de base aux pt int. numerique. */

/*     -- CALCUL DE F1,F2,FFM1,DFM1DP */

    e2ap2c_(&c__6, &c__7, poids, p25, dp25, ip, f1, f2, dfm1dp, poidel, dfm1, &coor[7]);

    hookax_(iopt, &car[1], elas);

/*     Contriantes [SIG] = [Elas] * [STRAIN_elas - STRAIN_ther] + SIG_0 */


/*     1) --- STRAIN Elastiques */
/*                                                             [U_r(1)] */
/* [Strain_zz]   [  0     0   ...   0   dp1dz dp2dz ... dp6dz] [  .   ] */
/* [         ]   [                                           ] [  .   ] */
/* [Strain_rr]   [dp1dr dp2dr ... dp6dr   0      0  ...   0  ] [  .   ] */
/* [         ]   [                                           ] [U_r(6)] */
/* [Strain_tt] = [ p1    p2        p6                        ] [      ] */
/* [         ]   [---   ---  ...  ---     0      0  ...   0  ] [U_z(1)] */
/* [2*Stra_rz]   [ r     r         r                         ] [  .   ] */
/*               [                                           ] [  .   ] */
/*               [ dp1dz .......  dp86z dp1dr ........  dp6dr] [  .   ] */
/*                                                             [U_z(6)] */

/*  barycentre F1(1) etr F2(1) */
    i__ = 5;
    r__ = f1[i__ - 1];
    p1 = p25[i__ * 6 - 6];
    p2 = p25[i__ * 6 - 5];
    p3 = p25[i__ * 6 - 4];
    p4 = p25[i__ * 6 - 3];
    p5 = p25[i__ * 6 - 2];
    p6 = p25[i__ * 6 - 1];
    dp1dr = dfm1dp[(i__ * 6 + 1 << 1) - 14];
    dp2dr = dfm1dp[(i__ * 6 + 2 << 1) - 14];
    dp3dr = dfm1dp[(i__ * 6 + 3 << 1) - 14];
    dp4dr = dfm1dp[(i__ * 6 + 4 << 1) - 14];
    dp5dr = dfm1dp[(i__ * 6 + 5 << 1) - 14];
    dp6dr = dfm1dp[(i__ * 6 + 6 << 1) - 14];
    dp1dz = dfm1dp[(i__ * 6 + 1 << 1) - 13];
    dp2dz = dfm1dp[(i__ * 6 + 2 << 1) - 13];
    dp3dz = dfm1dp[(i__ * 6 + 3 << 1) - 13];
    dp4dz = dfm1dp[(i__ * 6 + 4 << 1) - 13];
    dp5dz = dfm1dp[(i__ * 6 + 5 << 1) - 13];
    dp6dz = dfm1dp[(i__ * 6 + 6 << 1) - 13];
    estrain_zz__ = dp1dz * u[4] + dp2dz * u[6] + dp3dz * u[8] + dp4dz * u[10] + dp5dz * u[12] + dp6dz * u[14];
    estrain_rr__ = dp1dr * u[3] + dp2dr * u[5] + dp3dr * u[7] + dp4dr * u[9] + dp5dr * u[11] + dp6dr * u[13];
    estrain_tt__ = (p1 * u[3] + p2 * u[5] + p3 * u[7] + p4 * u[9] + p5 * u[11] + p6 * u[13]) / r__;
    estrain_rz__ = (dp1dz * u[3] + dp1dr * u[4] + dp2dz * u[5] + dp2dr * u[6] + dp3dz * u[7] + dp3dr * u[8] + dp4dz * u[9] + dp4dr * u[10] + dp5dz * u[11] + dp5dr * u[12] + dp6dz * u[13] + dp6dr * u[14]) * (float).5;

/*     2) --- STRAIN THermiques =  (ALPHA) * (P) * THETA */

/* [ Strain_zz ] [alpha_z] */
/* [           ] [       ]                         [ theta(1)] */
/* [ Strain_rr ] [alpha_r]                         [ theta(2)] */
/* [           ]=[       ]*[ p1 p2 p3  p4 p5 p6] * [ theta(3)] */
/* [ Strain_tt ] [alpha_t]                         [ theta(4)] */
/* [           ] [       ]                         [ theta(5)] */
/* [2*Strain_rz] [  0    ]                         [ theta(6)] */

    alpha_r__ = alpha[1];
    alpha_z__ = alpha[2];
    alpha_t__ = alpha[3];
    tstrain_zz__ = alpha_z__ * (p1 * theta[1] + p2 * theta[2] + p3 * theta[3] + p4 * theta[4] + p5 * theta[5] + p6 * theta[6]);
    tstrain_rr__ = alpha_r__ * (p1 * theta[1] + p2 * theta[2] + p3 * theta[3] + p4 * theta[4] + p5 * theta[5] + p6 * theta[6]);
    tstrain_tt__ = alpha_t__ * (p1 * theta[1] + p2 * theta[2] + p3 * theta[3] + p4 * theta[4] + p5 * theta[5] + p6 * theta[6]);

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
} /* etcap2c_ */

