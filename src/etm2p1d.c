/* etm2p1d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int etm2p1d_(coor, car, iopt, ae)
doublereal *coor, *car;
int32 *iopt;
doublereal *ae;
{
    static doublereal beta, beta2;
    static int32 i__;
    static doublereal delta, x21, y21, x31, y31, x32, y32, rho;

/* *************************************************************** */
/* but : calcul de la matrice de masse de l element tria 2p1d */
/* --- */
/* in : coor(noe,ndim) : coordones des 3 sommets. */
/*      car            : caracteristiques des materiaux */
/*                       car(1) = rho masse volumique */
/*      iopt           : ouvert si masse lumping ou autre ds futur */

/* out: ae            : matrice triangulaire sup */

/*  programmeur : modulef */
/* ................................................................. */

    /* Parameter adjustments */
    --ae;
    --car;
    coor -= 4;

    /* Function Body */
    rho = car[1];
    x21 = coor[5] - coor[4];
    y21 = coor[8] - coor[7];
    x31 = coor[6] - coor[4];
    y31 = coor[9] - coor[7];
    x32 = coor[6] - coor[5];
    y32 = coor[9] - coor[8];
    delta = x21 * y31 - x31 * y21;
    beta = rho * delta / 24.;
    beta2 = beta * 2.;

/*    -----   calcul de ae(ck)  ----- */

    for (i__ = 1; i__ <= 21; ++i__) {
	ae[i__] = 0.;
/* L1: */
    }
    ae[1] = beta2;
    ae[4] = beta;
    ae[6] = beta2;
    ae[11] = beta;
    ae[13] = beta;
    ae[15] = beta2;
    ae[3] = beta2;
    ae[8] = beta;
    ae[10] = beta2;
    ae[17] = beta;
    ae[19] = beta;
    ae[21] = beta2;

} /* etm2p1d_ */

