/* etcdktp.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int etcdktp_(coor, car, iopt, u, sigma)
doublereal *coor, *car;
int32 *iopt;
doublereal *u, *sigma;
{
    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static doublereal dxhx[9], dyhx[9], dxhy[9], dyhy[9], xksi, b[27]	/* was [3][9] */, d__[6];
    static int32 j, k;
    static doublereal coeff, delta, epais, x12sde, x31sde, y31sde, y12sde, yneta, young, p4, p5, p6, q4, t4, t5, t6, q5, q6, r4, r5, r6, x12, x31, x23, y23, y31, y12, dsigma[27], al4, al5, al6, poison;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* BUT: CALCUL DE LA MATRICE ELEMENTAIRE ELASTIQUE DES CONTRAINTES */
/*      AU BARYCENTRE DE L ELEMENT DE PLAQUE TRIA_DKTP . */
/* in : coor(noe,ndim) : coordonees des 3 sommets. */
/*      car            : caracteristiques des materiaux */
/*                       car(1) = masse volumique */
/*                       car(2) = epaisseur de la plaque */
/*      iopt           : ouvert */
/*      U(3,noe)       : U_z, Rotx, Roty */
/* out: SIGMA(3)      : contraintes elastiques */

/*  programmeur : modulef */
/* ............................................................... */

    /* Parameter adjustments */
    --sigma;
    u -= 4;
    --car;
    coor -= 4;

    /* Function Body */
    /* young = car[1];
    poison = car[2];
    epais = car[3]; */

    x23 = coor[5] - coor[6];
    y23 = coor[8] - coor[9];
    x31 = coor[6] - coor[4];
    y31 = coor[9] - coor[7];
    x12 = coor[4] - coor[5];
    y12 = coor[7] - coor[8];
    delta = x31 * y12 - x12 * y31;

    if (*iopt == 1) {
    young = car[1];
    poison = car[2];
    epais = car[3];
/*  MATRICE [ D ] (SYMETRIQUE STOKEE SOUS FORME TRIANGULAIRE => D(6)) */
/*             3*3 */
/* Computing 3rd power */
    d__1 = epais, d__2 = d__1;
/* Computing 2nd power */
    d__3 = poison;
    coeff = young * (d__2 * (d__1 * d__1)) / ((1. - d__3 * d__3) * 12.);
    d__[0] = coeff;
    d__[1] = coeff * poison;
    d__[2] = coeff;
    d__[3] = 0.;
    d__[4] = 0.;
    d__[5] = coeff * (1. - poison) / 2.; 
    } else  {
    epais  = car[7];
    coeff  = (epais * epais * epais) / 12.;
    d__[0] = coeff * car[1] ;
    d__[1] = coeff * car[2] ;
    d__[2] = coeff * car[3] ;
    d__[3] = coeff * car[4] ;
    d__[4] = coeff * car[5] ;
    d__[5] = coeff * car[6] ;
    }

/*  CALCUL DES LK ,PK ,TK ,QK ,RK */
/*                 (CF. J.L. BATOZ) */
/* Computing 2nd power */
    d__1 = x23;
/* Computing 2nd power */
    d__2 = y23;
    al4 = d__1 * d__1 + d__2 * d__2;
/* Computing 2nd power */
    d__1 = x31;
/* Computing 2nd power */
    d__2 = y31;
    al5 = d__1 * d__1 + d__2 * d__2;
/* Computing 2nd power */
    d__1 = x12;
/* Computing 2nd power */
    d__2 = y12;
    al6 = d__1 * d__1 + d__2 * d__2;
    p4 = x23 * -6. / al4;
    p5 = x31 * -6. / al5;
    p6 = x12 * -6. / al6;
    t4 = y23 * -6. / al4;
    t5 = y31 * -6. / al5;
    t6 = y12 * -6. / al6;
    q4 = x23 * 3. * y23 / al4;
    q5 = x31 * 3. * y31 / al5;
    q6 = x12 * 3. * y12 / al6;
    r4 = y23 * 3. * y23 / al4;
    r5 = y31 * 3. * y31 / al5;
    r6 = y12 * 3. * y12 / al6;

/*  CALCUL DES MATRICES [DX_HX] , [DY_HX] , [DX_HY] , [DY_HY] */
/*                           1*9       1*9       1*9       1*9 */

    xksi = .33333333333333331;
    yneta = xksi;
    dxhx[0] = p6 * (1. - xksi * 2.) + (p5 - p6) * yneta;
    dxhx[1] = q6 * (1. - xksi * 2.) - (q5 + q6) * yneta;
    dxhx[2] = (xksi + yneta) * 6. - 4. + r6 * (1. - xksi * 2.) - yneta * (r5 + r6);
    dxhx[3] = -p6 * (1. - xksi * 2.) + (p4 + p6) * yneta;
    dxhx[4] = q6 * (1. - xksi * 2.) - (q6 - q4) * yneta;
    dxhx[5] = xksi * 6. - 2. + r6 * (1. - xksi * 2.) + (r4 - r6) * yneta;
    dxhx[6] = -(p5 + p4) * yneta;
    dxhx[7] = (q4 - q5) * yneta;
    dxhx[8] = -(r5 - r4) * yneta;

    dyhx[0] = -p5 * (1. - yneta * 2.) - (p6 - p5) * xksi;
    dyhx[1] = q5 * (1. - yneta * 2.) - (q5 + q6) * xksi;
    dyhx[2] = (xksi + yneta) * 6. - 4. + r5 * (1. - yneta * 2.) - xksi * (r5 + r6);
    dyhx[3] = (p4 + p6) * xksi;
    dyhx[4] = (q4 - q6) * xksi;
    dyhx[5] = -(r6 - r4) * xksi;
    dyhx[6] = p5 * (1. - yneta * 2.) - (p4 + p5) * xksi;
    dyhx[7] = q5 * (1. - yneta * 2.) + (q4 - q5) * xksi;
    dyhx[8] = yneta * 6. - 2. + (r4 - r5) * xksi + r5 * (1. - yneta * 2.);

    dxhy[0] = t6 * (1. - xksi * 2.) + (t5 - t6) * yneta;
    dxhy[1] = r6 * (1. - xksi * 2.) + 1. - (r5 + r6) * yneta;
    dxhy[2] = -q6 * (1. - xksi * 2.) + yneta * (q5 + q6);
    dxhy[3] = -t6 * (1. - xksi * 2.) + (t4 + t6) * yneta;
    dxhy[4] = r6 * (1. - xksi * 2.) - 1. + (r4 - r6) * yneta;
    dxhy[5] = -q6 * (1. - xksi * 2.) - (q4 - q6) * yneta;
    dxhy[6] = -(t4 + t5) * yneta;
    dxhy[7] = (r4 - r5) * yneta;
    dxhy[8] = -(q4 - q5) * yneta;

    dyhy[0] = -t5 * (1. - yneta * 2.) - (t6 - t5) * xksi;
    dyhy[1] = r5 * (1. - yneta * 2.) + 1. - (r5 + r6) * xksi;
    dyhy[2] = -q5 * (1. - yneta * 2.) + xksi * (q5 + q6);
    dyhy[3] = (t4 + t6) * xksi;
    dyhy[4] = (r4 - r6) * xksi;
    dyhy[5] = -(q4 - q6) * xksi;
    dyhy[6] = t5 * (1. - yneta * 2.) - (t4 + t5) * xksi;
    dyhy[7] = r5 * (1. - yneta * 2.) - 1. + (r4 - r5) * xksi;
    dyhy[8] = -q5 * (1. - yneta * 2.) - (q4 - q5) * xksi;

/*  CALCUL DE LA MATRICE [ B ] (CF. J.L. BATOZ) */
/*                          3*9 */
    y31sde = y31 / delta;
    x31sde = x31 / delta;
    y12sde = y12 / delta;
    x12sde = x12 / delta;
    for (j = 1; j <= 9; ++j) {
	b[j * 3 - 3] = y31sde * dxhx[j - 1] + y12sde * dyhx[j - 1];
	b[j * 3 - 2] = -x31sde * dxhy[j - 1] - x12sde * dyhy[j - 1];
	b[j * 3 - 1] = -x31sde * dxhx[j - 1] - x12sde * dyhx[j - 1] + y31sde * dxhy[j - 1] + y12sde * dyhy[j - 1];
/* L1: */
    }

/*  CALCUL DE DSIGMA ( PASSAGE AU RANGEMENT PAR D.L. ) */

    k = 1;
    for (j = 1; j <= 9; ++j) {
	/*dsigma[k - 1] = d__[0] * b[j * 3 - 3] + d__[1] * b[j * 3 - 2];
	dsigma[k] = d__[1] * b[j * 3 - 3] + d__[2] * b[j * 3 - 2];
	dsigma[k + 1] = -d__[5] * b[j * 3 - 1];*/
    dsigma[k - 1] = d__[0] * b[j * 3 - 3] + d__[1] * b[j * 3 - 2] + d__[3] * b[j * 3 - 1];
    dsigma[k]     = d__[1] * b[j * 3 - 3] + d__[2] * b[j * 3 - 2] + d__[4] * b[j * 3 - 1];
    dsigma[k + 1] = d__[3] * b[j * 3 - 3] + d__[4] * b[j * 3 - 2] + d__[5] * b[j * 3 - 1];
  	k += 3;
/* L2: */
    }

/*     [      11]   [1 4 7 10 13 16 19 22 25]   [          ] */
/*     [SIGMA 22] = [2 5 8 11 14 17 20 23 26] * [u_solution] */
/*     [      12]   [3 6 9 12 15 18 21 24 27]   [          ] */
/*             3*1                3*9                     9*1 */

    sigma[1] = dsigma[0] * u[4] + dsigma[3] * u[5] + dsigma[6] * u[6] + dsigma[9] * u[7] + dsigma[12] * u[8] + dsigma[15] * u[9] + dsigma[18] * u[10] + dsigma[21] * u[11] + dsigma[24] * u[12];
    sigma[2] = dsigma[1] * u[4] + dsigma[4] * u[5] + dsigma[7] * u[6] + dsigma[10] * u[7] + dsigma[13] * u[8] + dsigma[16] * u[9] + dsigma[19] * u[10] + dsigma[22] * u[11] + dsigma[25] * u[12];
    sigma[3] = dsigma[2] * u[4] + dsigma[5] * u[5] + dsigma[8] * u[6] + dsigma[11] * u[7] + dsigma[14] * u[8] + dsigma[17] * u[9] + dsigma[20] * u[10] + dsigma[23] * u[11] + dsigma[26] * u[12];

} /* etcdktp_ */

