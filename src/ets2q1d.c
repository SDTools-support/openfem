/* ets2q1d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int ets2q1d_(coor, fomega, fgamma, pressi, norefs, alpha, theta, car, iopt, be)
doublereal *coor, *fomega, *fgamma, *pressi;
int32 *norefs;
doublereal *alpha, *theta, *car;
int32 *iopt;
doublereal *be;
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal delt1, delt2, delt3, delt4, a[3], c__, e[6];
    static int32 i__, j;
    static doublereal delta, young, unmnu;
    static int32 ij;
    static doublereal x21, y21, x41, y41, x32, y32, x42, y42, x43, y43, x31, y31, arelon, del[4], xji, yji, xnu, ynu, poisson;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* BUT: CALCUL DE SECONDS MEMBRES DE L ELEMENT QUAD 2Q1D */
/* --- */
/* in : */
/*   coor(noe,ndim)        : coordones des 4 sommets. */
/*   fomega(ndim,noe)      : fx, fy aux noeuds */
/*   fgamma(ndim,2*nbarete): fx, fy aux noeuds de chaque arete. */
/*   pressi(2*nbarete)     : pression aux noeuds de chaque arete */
/*   norefs(nbarete,2): norefs(i,1) = 0 si fgamma   = 0 sur arete_i */
/*                      norefs(i,2) = 0 si pression = 0 sur arete_i */
/*   theta(4): temperature aux 4 sommets */
/*   alpha(3): tenseur de dilatation thermique (xx, xy, yy) */

/*   car            : caracteristiques des materiaux */
/*   iopt = 1 isotrope Deformations Planes */
/*        = 2 isotrope Contraintes  Planes */
/*        = sinon anisotrope */
/*   car(6): caracteristiques des materiaux */
/*           if(iopt .eq. 1 .or. iopt.eq. 2) then */
/*             car(1) = young */
/*             car(2) = poisson */
/*           else */
/*             car: E11, E12, E22, E13, E23, E33 avec */

/*                     E11   E12   E13 */
/*                           E22   E23 */
/*                                 E33 */
/*            end if */

/* out: BE(8): second membre. */
/*  .................................................................... */

    /* Parameter adjustments */
    --be;
    --car;
    --theta;
    --alpha;
    norefs -= 5;
    --pressi;
    fgamma -= 3;
    fomega -= 3;
    coor -= 5;

    /* Function Body */
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

    del[0] = delt1 * .25;
    del[1] = delt2 * .25;
    del[2] = delt3 * .25;
    del[3] = delt4 * .25;


/* --  FOMEGA: INTEGRATION AUX 4 SOMMETS (F(0,0)+F(1,0)+F(0,1)+F(1,1) */

    for (i__ = 1; i__ <= 2; ++i__) {
	for (j = 1; j <= 4; ++j) {
	    ij = (j - 1 << 1) + i__;
	    be[ij] = fomega[i__ + (j << 1)] * del[j - 1];
/* L1: */
	}
    }

/* --  FGAMMA: INTEGRATION AUX 2 EXTREMITES DES COTES (F(0)+F(1))/2 */

    if (norefs[9] != 0) {
	xji = coor[6] - coor[5];
	yji = coor[10] - coor[9];
/* Computing 2nd power */
	d__1 = xji;
/* Computing 2nd power */
	d__2 = yji;
	arelon = sqrt(d__1 * d__1 + d__2 * d__2);
	xnu = yji / arelon;
	ynu = -xji / arelon;
    }
    delta = sqrt(x21 * x21 + y21 * y21) / 2.;
    be[1] += delta * (fgamma[3] - pressi[1] * xnu);
    be[3] += delta * (fgamma[5] - pressi[2] * xnu);
    be[2] += delta * (fgamma[4] - pressi[1] * ynu);
    be[4] += delta * (fgamma[6] - pressi[2] * ynu);

    if (norefs[10] != 0) {
	xji = coor[7] - coor[6];
	yji = coor[11] - coor[10];
/* Computing 2nd power */
	d__1 = xji;
/* Computing 2nd power */
	d__2 = yji;
	arelon = sqrt(d__1 * d__1 + d__2 * d__2);
	xnu = yji / arelon;
	ynu = -xji / arelon;
    }
    delta = sqrt(x32 * x32 + y32 * y32) / 2.;
    be[3] += delta * (fgamma[7] - pressi[3] * xnu);
    be[5] += delta * (fgamma[9] - pressi[4] * xnu);
    be[4] += delta * (fgamma[8] - pressi[3] * ynu);
    be[6] += delta * (fgamma[10] - pressi[4] * ynu);

    if (norefs[11] != 0) {
	xji = coor[8] - coor[7];
	yji = coor[12] - coor[11];
/* Computing 2nd power */
	d__1 = xji;
/* Computing 2nd power */
	d__2 = yji;
	arelon = sqrt(d__1 * d__1 + d__2 * d__2);
	xnu = yji / arelon;
	ynu = -xji / arelon;
    }
    delta = sqrt(x43 * x43 + y43 * y43) / 2.;
    be[5] += delta * (fgamma[11] - pressi[5] * xnu);
    be[7] += delta * (fgamma[13] - pressi[6] * xnu);
    be[6] += delta * (fgamma[12] - pressi[5] * ynu);
    be[8] += delta * (fgamma[14] - pressi[6] * ynu);

    if (norefs[12] != 0) {
	xji = coor[5] - coor[8];
	yji = coor[9] - coor[12];
/* Computing 2nd power */
	d__1 = xji;
/* Computing 2nd power */
	d__2 = yji;
	arelon = sqrt(d__1 * d__1 + d__2 * d__2);
	xnu = yji / arelon;
	ynu = -xji / arelon;
    }
    delta = sqrt(x41 * x41 + y41 * y41) / 2.;
    be[1] += delta * (fgamma[17] - pressi[8] * xnu);
    be[7] += delta * (fgamma[15] - pressi[7] * xnu);
    be[2] += delta * (fgamma[18] - pressi[8] * ynu);
    be[8] += delta * (fgamma[16] - pressi[7] * ynu);
/* L2: */

/*     CONTRAINTES THERMIQUES = HOOKE * ALPHA * THETA */
/*                         3,1     3,3   3,3       3,1 */
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
/* L3: */
	}
    }

/*     TENSEUR D'ELASTICITE * TENSEUR THERMIQUE */

    a[0] = alpha[1] * e[0] + alpha[3] * e[1] + alpha[2] * e[3];
    a[1] = alpha[1] * e[1] + alpha[3] * e[2] + alpha[2] * e[4];
    a[2] = alpha[1] * e[3] + alpha[3] * e[4] + alpha[2] * e[5];

    be[1] = be[1] + (x42 * a[2] * theta[1] - y42 * a[0] * theta[1]) / (float)4. + (x32 * a[2] * theta[2] - y32 * a[0] * theta[2]) / (float)4. + (x43 * a[2] * theta[4] - y43 * a[0] * theta[4]) / (float)4.;
    be[2] = be[2] + (x42 * a[1] * theta[1] - y42 * a[2] * theta[1]) / (float)4. + (x32 * a[1] * theta[2] - y32 * a[2] * theta[2]) / (float)4. + (x43 * a[1] * theta[4] - y43 * a[2] * theta[4]) / (float)4.;
    be[3] = be[3] + (-x41 * a[2] * theta[1] + y41 * a[0] * theta[1]) / (float)4. + (-x31 * a[2] * theta[2] + y31 * a[0] * theta[2]) / (float)4. + (x43 * a[2] * theta[3] - y43 * a[0] * theta[3]) / (float)4.;
    be[4] = be[4] + (-x41 * a[1] * theta[1] + y41 * a[2] * theta[1]) / (float)4. + (-x31 * a[1] * theta[2] + y31 * a[2] * theta[2]) / (float)4. + (x43 * a[1] * theta[3] - y43 * a[2] * theta[3]) / (float)4.;
    be[5] = be[5] + (x21 * a[2] * theta[2] - y21 * a[0] * theta[2]) / (float)4. + (-x42 * a[2] * theta[3] + y42 * a[0] * theta[3]) / (float)4. + (-x41 * a[2] * theta[4] + y41 * a[0] * theta[4]) / (float)4.;
    be[6] = be[6] + (x21 * a[1] * theta[2] - y21 * a[2] * theta[2]) / (float)4. + (-x42 * a[1] * theta[3] + y42 * a[2] * theta[3]) / (float)4. + (-x41 * a[1] * theta[4] + y41 * a[2] * theta[4]) / (float)4.;
    be[7] = be[7] + (x21 * a[2] * theta[1] - y21 * a[0] * theta[1]) / (float)4. + (x32 * a[2] * theta[3] - y32 * a[0] * theta[3]) / (float)4. + (x31 * a[2] * theta[4] - y31 * a[0] * theta[4]) / (float)4.;
    be[8] = be[8] + (x21 * a[1] * theta[1] - y21 * a[2] * theta[1]) / (float)4. + (x32 * a[1] * theta[3] - y32 * a[2] * theta[3]) / (float)4. + (x31 * a[1] * theta[4] - y31 * a[2] * theta[4]) / (float)4.;
} /* ets2q1d_ */

