/* etm2q1d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int etm2q1d_(coor, car, iopt, ae)
doublereal *coor, *car;
int32 *iopt;
doublereal *ae;
{
    static doublereal beta4, delt1, delt2, delt3, delt4;
    static int32 i__;
    static doublereal x21, y21, x31, y31, x32, y32, x41, y41, x43, y43, x42, y42, rho;

/* *************************************************************** */
/* but : calcul de la matrice de masse de l element quad 2q1d */
/* --- */
/* in : coor(noe,ndim) : coordones des 4 sommets. */
/*      car            : caracteristiques des materiaux */
/*                       car(1) = rho masse volumique */
/*      iopt           : ouvert si masse lumping ou autre ds futur */

/* out: ae(8)          : matrice DIAGONALE !! */

/*  programmeur : modulef */
/* ................................................................. */

    /* Parameter adjustments */
    --ae;
    --car;
    coor -= 5;

    /* Function Body */
    rho = car[1];
    x21 = coor[6] - coor[5];
    y21 = coor[10] - coor[9];
    x31 = coor[7] - coor[5];
    y31 = coor[11] - coor[9];
    x32 = coor[7] - coor[6];
    y32 = coor[11] - coor[10];
    x41 = coor[8] - coor[5];
    y41 = coor[12] - coor[9];
    x42 = coor[8] - coor[6];
    y42 = coor[12] - coor[10];
    x43 = coor[8] - coor[7];
    y43 = coor[12] - coor[11];
    delt1 = x21 * y41 - x41 * y21;
    delt2 = x21 * y32 - x32 * y21;
    delt3 = x32 * y43 - x43 * y32;
    delt4 = x41 * y43 - x43 * y41;

    for (i__ = 1; i__ <= 8; ++i__) {
	ae[i__] = 0.;
/* L1: */
    }
    beta4 = rho / 4.;

    ae[1] = delt1 * beta4;
    ae[3] = delt2 * beta4;
    ae[5] = delt3 * beta4;
    ae[7] = delt4 * beta4;
    ae[2] = delt1 * beta4;
    ae[4] = delt2 * beta4;
    ae[6] = delt3 * beta4;
    ae[8] = delt4 * beta4;
} /* etm2q1d_ */

