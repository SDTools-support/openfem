/* etr2p1d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int etr2p1d_(coor, car, iopt, ae)
doublereal *coor, *car;
int32 *iopt;
doublereal *ae;
{
    static doublereal c__, e[6];
    static int32 i__;
    static doublereal delta, young, unmnu, x21, y21, x31, y31, x32, y32, poisson;

/*  .................................................................... */
/* but : calcul de la matrice de rigidite de l element tria 2p1d */
/* --- */
/* in : coor(noe,ndim) : coordonnees des 3 sommets. */
/*      car            : caracteristiques des materiaux */
/*      iopt = 1 isotrope Contraintes  Planes */
/*           = 2 isotrope Deformations Planes */
/*           = sinon anisotrope */
/*      car(6): caracteristiques des materiaux */
/*              if(iopt .eq. 1 .or. iopt.eq. 2) then */
/*                car(1) = young */
/*                car(2) = poisson */
/*              else */
/*                car: E11, E12, E22, E13, E23, E33 avec */

/*                     E11   E12   E13 */
/*                           E22   E23 */
/*                                 E33 */
/*              end if */

/* out: ae(21): matrice triangulaire sup */
/* ..................................................................... */

    /* Parameter adjustments */
    --ae;
    --car;
    coor -= 4;

    /* Function Body */
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

    x21 = coor[5] - coor[4];
    y21 = coor[8] - coor[7];
    x31 = coor[6] - coor[4];
    y31 = coor[9] - coor[7];
    x32 = coor[6] - coor[5];
    y32 = coor[9] - coor[8];
    delta = x21 * y31 - x31 * y21;

    for (i__ = 1; i__ <= 6; ++i__) {
	e[i__ - 1] = e[i__ - 1] * .5 / delta;
/* L2: */
    }

/*  -----      COEFFICIENT DE LA MATRICE ELEMENTAIRE     ----- */

    ae[1] = e[0] * y32 * y32 - e[3] * 2. * x32 * y32 + e[5] * x32 * x32;
    ae[4] = -e[0] * y31 * y32 + e[3] * (x31 * y32 + x32 * y31) - e[5] * x31 * x32;
    ae[6] = e[0] * y31 * y31 - e[3] * 2. * x31 * y31 + e[5] * x31 * x31;
    ae[11] = e[0] * y21 * y32 - e[3] * (x21 * y32 + x32 * y21) + e[5] * x21 * x32;
    ae[13] = -e[0] * y21 * y31 + e[3] * (x21 * y31 + x31 * y21) - e[5] * x21 * x31;
    ae[15] = e[0] * y21 * y21 - e[3] * 2. * x21 * y21 + e[5] * x21 * x21;
    ae[2] = -(e[1] + e[5]) * x32 * y32 + e[3] * y32 * y32 + e[4] * x32 * x32;
    ae[5] = e[1] * x32 * y31 - e[3] * y32 * y31 - e[4] * x32 * x31 + e[5] * x31 * y32;
    ae[12] = -e[1] * x32 * y21 + e[3] * y32 * y21 + e[4] * x32 * x21 - e[5] * x21 * y32;
    ae[3] = e[2] * x32 * x32 - e[4] * 2. * x32 * y32 + e[5] * y32 * y32;
    ae[7] = e[1] * x31 * y32 - e[3] * y31 * y32 - e[4] * x31 * x32 + e[5] * x32 * y31;
    ae[9] = -(e[1] + e[5]) * x31 * y31 + e[3] * y31 * y31 + e[4] * x31 * x31;
    ae[14] = e[1] * x31 * y21 - e[3] * y31 * y21 - e[4] * x31 * x21 + e[5] * x21 * y31;
    ae[8] = -e[2] * x31 * x32 + e[4] * (x32 * y31 + x31 * y32) - e[5] * y31 * y32;
    ae[10] = e[2] * x31 * x31 - e[4] * 2. * x31 * y31 + e[5] * y31 * y31;
    ae[16] = -e[1] * x21 * y32 + e[3] * y21 * y32 + e[4] * x21 * x32 - e[5] * x32 * y21;
    ae[18] = e[1] * x21 * y31 - e[3] * y21 * y31 - e[4] * x21 * x31 + e[5] * x31 * y21;
    ae[20] = -(e[1] + e[5]) * x21 * y21 + e[3] * y21 * y21 + e[4] * x21 * x21;
    ae[17] = e[2] * x21 * x32 - e[4] * (x32 * y21 + x21 * y32) + e[5] * y21 * y32;
    ae[19] = -e[2] * x21 * x31 + e[4] * (x31 * y21 + x21 * y31) - e[5] * y21 * y31;
    ae[21] = e[2] * x21 * x21 - e[4] * 2. * x21 * y21 + e[5] * y21 * y21;

} /* etr2p1d_ */

