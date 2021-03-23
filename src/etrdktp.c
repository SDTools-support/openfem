/* etrdktp.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__9 = 9;
static int32 c__3 = 3;

/* Subroutine */ int etrdktp_(coor, car, iopt, ae)
doublereal *coor, *car;
int32 *iopt;
doublereal *ae;
{
    /* Initialized data */

    static doublereal xksi[3] = { .5,.5,0. };
    static doublereal yneta[3] = { 0.,.5,.5 };

    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static doublereal dxhx[9], dyhx[9], dxhy[9], dyhy[9], b[27]	/* was [3][9] */, d__[6];
    static int32 i__, j;
    static doublereal coeff, delta, epais, young, p4, p5, p6, q4, q5, q6, t4, t5, t6, r4, r5, r6, x12, x31, x23, y23, y31, y12;
    extern /* Subroutine */ int tabaxd_();
    static doublereal al4, al5, al6, poison, aux[27]	/* was [9][3] */;

/* *************************************************************** */
/* but: calcul de la matrice de rigidite de l element de plaque */
/*      tria DKTP (DISCRETE KIRCHHOFF THEORY) */
/* in : coor(noe,ndim) : coordonees des 3 sommets. */
/*      car(3)         : caracteristiques des materiaux */
/*                       young, poisson, epaisseur */
/*      iopt           : ouvert si masse lumping ou autre ds futur */

/* out: ae            : matrice triangulaire sup */

/*  programmeur : modulef */
/* ............................................................... */
    /* Parameter adjustments */
    --ae;
    --car;
    coor -= 4;

    /* Function Body */


    x23 = coor[5] - coor[6];
    y23 = coor[8] - coor[9];
    x31 = coor[6] - coor[4];
    y31 = coor[9] - coor[7];
    x12 = coor[4] - coor[5];
    y12 = coor[7] - coor[8];
    delta = x31 * y12 - x12 * y31;


/*  CALCUL DES LK ,PK ,QK ,TK ,RK */
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
/* Computing 2nd power */
    d__1 = y23;
    r4 = d__1 * d__1 * 3. / al4;
/* Computing 2nd power */
    d__1 = y31;
    r5 = d__1 * d__1 * 3. / al5;
/* Computing 2nd power */
    d__1 = y12;
    r6 = d__1 * d__1 * 3. / al6;
/*     CAS OU DB EST FOURNI */
    if (*iopt == 3) {
/*	d__[0] = car[1];
	d__[1] = car[2];
	d__[2] = car[3];
	d__[3] = car[4];
	d__[4] = car[5];
	d__[5] = car[6]; */
    epais = car[7];
    coeff = (epais * epais * epais) / 12.;
    d__[0] = coeff * car[1] / (delta * 6.);
    d__[1] = coeff * car[2] / (delta * 6.);
    d__[2] = coeff * car[3] / (delta * 6.);
    d__[3] = coeff * car[4] / (delta * 6.);
    d__[4] = coeff * car[5] / (delta * 6.);
    d__[5] = coeff * car[6] / (delta * 6.);

/*  MATRICE [ D ] (SYMETRIQUE STOKEE SOUS FORME TRIANGULAIRE => D(6)) */
/*             3*3 */
    } else {
	young = car[1];
	poison = car[2];
	epais = car[3];
/* Computing 3rd power */
	d__1 = epais, d__2 = d__1;
/* Computing 2nd power */
	d__3 = poison;
	coeff = young * (d__2 * (d__1 * d__1)) / ((1. - d__3 * d__3) * 12.);
	d__[0] = coeff / (delta * 6.);
	d__[1] = d__[0] * poison;
	d__[2] = d__[0];
	d__[3] = 0.;
	d__[4] = 0.;
	d__[5] = coeff * (1. - poison) / (delta * 12.);
    }

/*  BOUCLE SUR LES 3 POINTS D'INTEGRATION NUMERIQUES */
/*   ( MILIEU DES 3 COTES DU TRIANGLE DE REFERENCE ) */

    for (i__ = 1; i__ <= 3; ++i__) {

/*  CALCUL DES MATRICES [DX_HX] , [DY_HX] , [DX_HY] , [DY_HY] */
/*                           1*9       1*9       1*9       1*9 */
	dxhx[0] = p6 * (1. - xksi[i__ - 1] * 2.) + (p5 - p6) * yneta[i__ - 1];
	dxhx[1] = q6 * (1. - xksi[i__ - 1] * 2.) - (q5 + q6) * yneta[i__ - 1];
	dxhx[2] = (xksi[i__ - 1] + yneta[i__ - 1]) * 6. - 4. + r6 * (1. - xksi[i__ - 1] * 2.) - yneta[i__ - 1] * (r5 + r6);
	dxhx[3] = -p6 * (1. - xksi[i__ - 1] * 2.) + (p4 + p6) * yneta[i__ - 1];
	dxhx[4] = q6 * (1. - xksi[i__ - 1] * 2.) - (q6 - q4) * yneta[i__ - 1];
	dxhx[5] = xksi[i__ - 1] * 6. - 2. + r6 * (1. - xksi[i__ - 1] * 2.) + (r4 - r6) * yneta[i__ - 1];
	dxhx[6] = -(p5 + p4) * yneta[i__ - 1];
	dxhx[7] = (q4 - q5) * yneta[i__ - 1];
	dxhx[8] = -(r5 - r4) * yneta[i__ - 1];

	dyhx[0] = -p5 * (1. - yneta[i__ - 1] * 2.) - (p6 - p5) * xksi[i__ - 1];
	dyhx[1] = q5 * (1. - yneta[i__ - 1] * 2.) - (q5 + q6) * xksi[i__ - 1];
	dyhx[2] = (xksi[i__ - 1] + yneta[i__ - 1]) * 6. - 4. + r5 * (1. - yneta[i__ - 1] * 2.) - xksi[i__ - 1] * (r5 + r6);
	dyhx[3] = (p4 + p6) * xksi[i__ - 1];
	dyhx[4] = (q4 - q6) * xksi[i__ - 1];
	dyhx[5] = -(r6 - r4) * xksi[i__ - 1];
	dyhx[6] = p5 * (1. - yneta[i__ - 1] * 2.) - (p4 + p5) * xksi[i__ - 1];
	dyhx[7] = q5 * (1. - yneta[i__ - 1] * 2.) + (q4 - q5) * xksi[i__ - 1];
	dyhx[8] = yneta[i__ - 1] * 6. - 2. + (r4 - r5) * xksi[i__ - 1] + r5 * (1. - yneta[i__ - 1] * 2.);

	dxhy[0] = t6 * (1. - xksi[i__ - 1] * 2.) + (t5 - t6) * yneta[i__ - 1];
	dxhy[1] = r6 * (1. - xksi[i__ - 1] * 2.) + 1. - (r5 + r6) * yneta[i__ - 1];
	dxhy[2] = -q6 * (1. - xksi[i__ - 1] * 2.) + yneta[i__ - 1] * (q5 + q6);
	dxhy[3] = -t6 * (1. - xksi[i__ - 1] * 2.) + (t4 + t6) * yneta[i__ - 1];
	dxhy[4] = r6 * (1. - xksi[i__ - 1] * 2.) - 1. + (r4 - r6) * yneta[i__ - 1];
	dxhy[5] = -q6 * (1. - xksi[i__ - 1] * 2.) - (q4 - q6) * yneta[i__ - 1];
	dxhy[6] = -(t4 + t5) * yneta[i__ - 1];
	dxhy[7] = (r4 - r5) * yneta[i__ - 1];
	dxhy[8] = -(q4 - q5) * yneta[i__ - 1];

	dyhy[0] = -t5 * (1. - yneta[i__ - 1] * 2.) - (t6 - t5) * xksi[i__ - 1];
	dyhy[1] = r5 * (1. - yneta[i__ - 1] * 2.) + 1. - (r5 + r6) * xksi[i__ - 1];
	dyhy[2] = -q5 * (1. - yneta[i__ - 1] * 2.) + xksi[i__ - 1] * (q5 + q6);
	dyhy[3] = (t4 + t6) * xksi[i__ - 1];
	dyhy[4] = (r4 - r6) * xksi[i__ - 1];
	dyhy[5] = -(q4 - q6) * xksi[i__ - 1];
	dyhy[6] = t5 * (1. - yneta[i__ - 1] * 2.) - (t4 + t5) * xksi[i__ - 1];
	dyhy[7] = r5 * (1. - yneta[i__ - 1] * 2.) - 1. + (r4 - r5) * xksi[i__ - 1];
	dyhy[8] = -q5 * (1. - yneta[i__ - 1] * 2.) - (q4 - q5) * xksi[i__ - 1];

/*     --CALCUL DE LA MATRICE [ B ] (CF. J.L. BATOZ) */
/*                          3*9 */
	for (j = 1; j <= 9; ++j) {
	    b[j * 3 - 3] = y31 * dxhx[j - 1] + y12 * dyhx[j - 1];
	    b[j * 3 - 2] = -x31 * dxhy[j - 1] - x12 * dyhy[j - 1];
	    b[j * 3 - 1] = -x31 * dxhx[j - 1] - x12 * dyhx[j - 1] + y31 * dxhy[j - 1] + y12 * dyhy[j - 1];
/* L1: */
	}
/*                                 T */
/*     --CALCUL DU PRODUIT [ AE ] = [ B ] * [ D ] * [ B ] */
/*                             9*9     9*3     3*3     3*9 */
/*     ( [ AE ] SYMETRIQUE STOKEE SOUS FORME PROFIL => AE(45) ) */

	tabaxd_(&c__9, &c__3, b, d__, &ae[1], aux, &i__);

/* L2: */
    }
} /* etrdktp_ */

