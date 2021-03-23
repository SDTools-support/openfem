/* etr2q1d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int etr2q1d_(coor, car, iopt, ae)
doublereal *coor, *car;
int32 *iopt;
doublereal *ae;
{
    static doublereal delt1, delt2, delt3, delt4, c__, e[6];
    static int32 i__;
    static doublereal young, unmnu, x21, y21, x31, y31, x32, y32, x41, y41, x43, y43, x42, y42, x21c, y21c, x41c, x42c, y42c, y41c, x32c, y32c, x31c, y31c, x43c, y43c, del1, del2, del3, del4, poisson;

/*  .................................................................... */
/* but : calcul de la matrice de rigidite de l element quad 2q1d */
/* --- */
/* in : coor(noe,ndim) : coordones des 4 sommets. */
/*      car            : caracteristiques des materiaux */
/*      iopt = 1 isotrope Deformations Planes */
/*           = 2 isotrope Contraintes  Planes */
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

/* out: ae(36): matrice triangulaire sup */
/* ..................................................................... */

    /* Parameter adjustments */
    --ae;
    --car;
    coor -= 5;

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

    y42c = y42 * y42;
    x42c = x42 * x42;
    y41c = y41 * y41;
    x41c = x41 * x41;
    y21c = y21 * y21;
    x21c = x21 * x21;
    y32c = y32 * y32;
    x32c = x32 * x32;
    y31c = y31 * y31;
    x31c = x31 * x31;
    y43c = y43 * y43;
    x43c = x43 * x43;
    del1 = .25 / delt1;
    del2 = .25 / delt2;
    del3 = .25 / delt3;
    del4 = .25 / delt4;

    ae[1] = (e[0] * y42c + e[5] * x42c - e[3] * 2. * x42 * y42) * del1 + (e[0] * y43c + e[5] * x43c - e[3] * 2. * x43 * y43) * del4 + (e[0] * y32c + e[5] * x32c - e[3] * 2. * x32 * y32) * del2;
    ae[4] = -(e[0] * y42 * y41 + e[5] * x41 * x42 - e[3] * (x41 * y42 + x42 * y41)) * del1 - (e[0] * y31 * y32 + e[5] * x31 * x32 - e[3] * (x31 * y32 + x32 * y31)) * del2;
    ae[6] = (e[0] * y41c + e[5] * x41c - e[3] * 2. * x41 * y41) * del1 + (e[0] * y31c + e[5] * x31c - e[3] * 2. * x31 * y31) * del2 + (e[0] * y43c + e[5] * x43c - e[3] * 2. * x43 * y43) * del3;
    ae[11] = (e[0] * y21 * y32 + e[5] * x21 * x32 - e[3] * (x21 * y32 + x32 * y21)) * del2 - (e[0] * y41 * y43 + e[5] * x41 * x43 - e[3] * (x41 * y43 + x43 * y41)) * del4;
    ae[13] = -(e[0] * y21 * y31 + e[5] * x21 * x31 - e[3] * (x21 * y31 + x31 * y21)) * del2 - (e[0] * y42 * y43 + e[5] * x43 * x42 - e[3] * (x42 * y43 + x43 * y42)) * del3;
    ae[15] = (e[0] * y21c + e[5] * x21c - e[3] * 2. * x21 * y21) * del2 + (e[0] * y42c + e[5] * x42c - e[3] * 2. * x42 * y42) * del3 + (e[0] * y41c + e[5] * x41c - e[3] * 2. * x41 * y41) * del4;
    ae[22] = (e[0] * y42 * y21 + e[5] * x42 * x21 - e[3] * (x21 * y42 + x42 * y21)) * del1 + (e[0] * y31 * y43 + e[5] * x31 * x43 - e[3] * (x31 * y43 + x43 * y31)) * del4;
    ae[24] = -(e[0] * y41 * y21 + e[5] * x41 * x21 - e[3] * (x21 * y41 + x41 * y21)) * del1 + (e[0] * y32 * y43 + e[5] * x32 * x43 - e[3] * (x32 * y43 + x43 * y32)) * del3;
    ae[26] = -(e[0] * y32 * y42 + e[5] * x32 * x42 - e[3] * (x32 * y42 + x42 * y32)) * del3 - (e[0] * y31 * y41 + e[5] * x31 * x41 - e[3] * (x31 * y41 + x41 * y31)) * del4;
    ae[28] = (e[0] * y21c + e[5] * x21c - e[3] * 2. * x21 * y21) * del1 + (e[0] * y32c + e[5] * x32c - e[3] * 2. * x32 * y32) * del3 + (e[0] * y31c + e[5] * x31c - e[3] * 2. * x31 * y31) * del4;
    ae[2] = -((e[1] + e[5]) * x42 * y42 - e[3] * y42c - e[4] * x42c) * del1 - ((e[1] + e[5]) * x32 * y32 - e[3] * y32c - e[4] * x32c) * del2 - ((e[1] + e[5]) * y43 * x43 - e[3] * y43c - e[4] * x43c) * del4;
    ae[5] = (e[1] * x42 * y41 + e[5] * y42 * x41 - e[3] * y42 * y41 - e[4] * x42 * x41) * del1 + (e[1] * x32 * y31 + e[5] * y32 * x31 - e[3] * y32 * y31 - e[4] * x32 * x31) * del2;
    ae[12] = -(e[1] * x32 * y21 + e[5] * x21 * y32 - e[3] * y32 * y21 - e[4] * x32 * x21) * del2 + (e[1] * x43 * y41 + e[5] * y43 * x41 - e[3] * y43 * y41 - e[4] * x43 * x41) * del4;
    ae[23] = -(e[1] * x42 * y21 + e[5] * y42 * x21 - e[3] * y42 * y21 - e[4] * x42 * x21) * del1 - (e[1] * x43 * y31 + e[5] * y43 * x31 - e[3] * y43 * y31 - e[4] * x43 * x31) * del4;
    ae[3] = (e[5] * y42c + e[2] * x42c - e[4] * 2. * x42 * y42) * del1 + (e[5] * y43c + e[2] * x43c - e[4] * 2. * x43 * y43) * del4 + (e[5] * y32c + e[2] * x32c - e[4] * 2. * x32 * y32) * del2;
    ae[7] = (e[1] * x41 * y42 + e[5] * y41 * x42 - e[3] * y41 * y42 - e[4] * x41 * x42) * del1 + (e[1] * x31 * y32 + e[5] * x32 * y31 - e[3] * y31 * y32 - e[4] * x31 * x32) * del2;
    ae[9] = -((e[1] + e[5]) * x41 * y41 - e[3] * y41c - e[4] * x41c) * del1 - ((e[1] + e[5]) * x31 * y31 - e[3] * y31c - e[4] * x31c) * del2 - ((e[1] + e[5]) * x43 * y43 - e[3] * y43c - e[4] * x43c) * del3;
    ae[14] = (e[1] * x31 * y21 + e[5] * x21 * y31 - e[3] * y31 * y21 - e[4] * x31 * x21) * del2 + (e[1] * x43 * y42 + e[5] * x42 * y43 - e[3] * y43 * y42 - e[4] * x43 * x42) * del3;
    ae[25] = (e[1] * x41 * y21 + e[5] * x21 * y41 - e[3] * y41 * y21 - e[4] * x41 * x21) * del1 - (e[1] * x43 * y32 + e[5] * y43 * x32 - e[3] * y43 * y32 - e[4] * x43 * x32) * del3;
    ae[8] = -(e[5] * y42 * y41 + e[2] * x41 * x42 - e[4] * (x42 * y41 + x41 * y42)) * del1 - (e[5] * y31 * y32 + e[2] * x31 * x32 - e[4] * (x32 * y31 + x31 * y32)) * del2;
    ae[10] = (e[5] * y41c + e[2] * x41c - e[4] * 2. * x41 * y41) * del1 + (e[5] * y31c + e[2] * x31c - e[4] * 2. * x31 * y31) * del2 + (e[5] * y43c + e[2] * x43c - e[4] * 2. * x43 * y43) * del3;
    ae[16] = -(e[1] * x21 * y32 + e[5] * y21 * x32 - e[3] * y21 * y32 - e[4] * x21 * x32) * del2 + (e[1] * x41 * y43 + e[5] * x43 * y41 - e[3] * y41 * y43 - e[4] * x41 * x43) * del4;
    ae[18] = (e[1] * x21 * y31 + e[5] * x31 * y21 - e[3] * y21 * y31 - e[4] * x21 * x31) * del2 + (e[1] * x42 * y43 + e[5] * y42 * x43 - e[3] * y42 * y43 - e[4] * x42 * x43) * del3;
    ae[20] = -((e[1] + e[5]) * y21 * x21 - e[3] * y21c - e[4] * x21c) * del2 - ((e[1] + e[5]) * x42 * y42 - e[3] * y42c - e[4] * x42c) * del3 - ((e[1] + e[5]) * x41 * y41 - e[3] * y41c - e[4] * x41c) * del4;
    ae[27] = (e[1] * x42 * y32 + e[5] * y42 * x32 - e[3] * y42 * y32 - e[4] * x42 * x32) * del3 + (e[1] * x41 * y31 + e[5] * y41 * x31 - e[3] * y41 * y31 - e[4] * x41 * x31) * del4;
    ae[17] = (e[5] * y21 * y32 + e[2] * x21 * x32 - e[4] * (y21 * x32 + y32 * x21)) * del2 - (e[5] * y41 * y43 + e[2] * x41 * x43 - e[4] * (y41 * x43 + y43 * x41)) * del4;
    ae[19] = -(e[5] * y21 * y31 + e[2] * x21 * x31 - e[4] * (y21 * x31 + y31 * x21)) * del2 - (e[5] * y42 * y43 + e[2] * x43 * x42 - e[4] * (y42 * x43 + y43 * x42)) * del3;
    ae[21] = (e[5] * y21c + e[2] * x21c - e[4] * 2. * x21 * y21) * del2 + (e[5] * y42c + e[2] * x42c - e[4] * 2. * x42 * y42) * del3 + (e[5] * y41c + e[2] * x41c - e[4] * 2. * x41 * y41) * del4;
    ae[29] = -(e[1] * x21 * y42 + e[5] * x42 * y21 - e[3] * y21 * y42 - e[4] * x21 * x42) * del1 - (e[1] * x31 * y43 + e[5] * y31 * x43 - e[3] * y31 * y43 - e[4] * x31 * x43) * del4;
    ae[31] = (e[1] * x21 * y41 + e[5] * x41 * y21 - e[3] * y21 * y41 - e[4] * x21 * x41) * del1 - (e[1] * x32 * y43 + e[5] * x43 * y32 - e[3] * y32 * y43 - e[4] * x32 * x43) * del3;
    ae[33] = (e[1] * x32 * y42 + e[5] * y32 * x42 - e[3] * y32 * y42 - e[4] * x32 * x42) * del3 + (e[1] * x31 * y41 + e[5] * x41 * y31 - e[3] * y31 * y41 - e[4] * x31 * x41) * del4;
    ae[35] = -((e[5] + e[1]) * x21 * y21 - e[3] * y21c - e[4] * x21c) * del1 - ((e[1] + e[5]) * x32 * y32 - e[3] * y32c - e[4] * x32c) * del3 - ((e[1] + e[5]) * x31 * y31 - e[3] * y31c - e[4] * x31c) * del4;
    ae[30] = (e[5] * y42 * y21 + e[2] * x42 * x21 - e[4] * (x42 * y21 + x21 * y42)) * del1 + (e[5] * y31 * y43 + e[2] * x31 * x43 - e[4] * (y31 * x43 + y43 * x31)) * del4;
    ae[32] = -(e[5] * y41 * y21 + e[2] * x41 * x21 - e[4] * (y21 * x41 + y41 * x21)) * del1 + (e[5] * y32 * y43 + e[2] * x32 * x43 - e[4] * (x43 * y32 + x32 * y43)) * del3;
    ae[34] = -(e[5] * y32 * y42 + e[2] * x32 * x42 - e[4] * (y32 * x42 + x32 * y42)) * del3 - (e[5] * y31 * y41 + e[2] * x31 * x41 - e[4] * (y31 * x41 + y41 * x31)) * del4;
    ae[36] = (e[5] * y21c + e[2] * x21c - e[4] * 2. * x21 * y21) * del1 + (e[5] * y32c + e[2] * x32c - e[4] * 2. * x32 * y32) * del3 + (e[5] * y31c + e[2] * x31c - e[4] * 2. * x31 * y31) * del4;

} /* etr2q1d_ */

