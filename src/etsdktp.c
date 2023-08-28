/* etsdktp.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int etsdktp_(coor, fomega, be)
doublereal *coor, *fomega, *be;
{
    static doublereal d__, delta, x12, x31, x23, y23, y31, y12;

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* BUT   : CALCUL DES SECONDS MEMBRES DE L ELEMENT DE PLAQUE */
/* ---     TRIA-DKTP (J.L. BATOZ). */
/* in : coor(noe,ndim) : coordonees des 3 sommets. */
/*      FOMEGA(3): Fz  surfaciques aux sommets 1 , 2, 3. */
/* out: BE(9) LES  SECONDS MEMBRES ELEMENTAIRES. */

/*  programmeur : modulef */
/* ............................................................... */

    /* Parameter adjustments */
    --be;
    --fomega;
    coor -= 4;

    /* Function Body */
    x23 = coor[5] - coor[6];
    y23 = coor[8] - coor[9];
    x31 = coor[6] - coor[4];
    y31 = coor[9] - coor[7];
    x12 = coor[4] - coor[5];
    y12 = coor[7] - coor[8];
    delta = x31 * y12 - x12 * y31;
    d__ = delta / 6.;
    be[1] = fomega[1] * d__;
    be[2] = 0.;
    be[3] = 0.;
    be[4] = fomega[2] * d__;
    be[5] = 0.;
    be[6] = 0.;
    be[7] = fomega[3] * d__;
    be[8] = 0.;
    be[9] = 0.;
} /* etsdktp_ */

