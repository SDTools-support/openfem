/* ets2p1d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int ets2p1d_(coor, fomega, fgamma, pressi, norefs, alpha, theta, car, iopt, be)
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
    static doublereal delt1, a[3], c__, e[6];
    static int32 i__, j;
    static doublereal delta, somme, young, unmnu;
    static int32 ij;
    static doublereal x21, y21, x31, y31, x32, y32, arelon, xji, yji, xnu, ynu, poisson;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* BUT: CALCUL DE SECONDS MEMBRES DE L ELEMENT TRIA 2P1D */
/* --- */
/* in : */
/*   coor(noe,ndim)        : coordonnees des 3 sommets. */
/*   fomega(ndim,noe)      : fx, fy aux noeuds */
/*   fgamma(ndim,2*nbarete): fx, fy aux noeuds de chaque arete. */
/*   pressi(2*nbarete)     : pression aux noeuds de chaque arete */
/*   norefs(nbarete,2): norefs(i,1) = 0 si fgamma   = 0 sur arete_i */
/*                      norefs(i,2) = 0 si pression = 0 sur arete_i */
/*   theta(3): temperature aux 3 sommets */
/*   alpha(3): tenseur de dilatation thermique (k11, k12, k22) */

/*   car            : caracteristiques des materiaux */
/*   iopt = 1 isotrope Contraintes  Planes */
/*        = 2 isotrope Deformations Planes */
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

/* out: BE(6): second membre. */
/*  .................................................................... */

    /* Parameter adjustments */
    --be;
    --car;
    --theta;
    --alpha;
    norefs -= 4;
    --pressi;
    fgamma -= 3;
    fomega -= 3;
    coor -= 4;

    /* Function Body */
    x21 = coor[5] - coor[4];
    y21 = coor[8] - coor[7];
    x31 = coor[6] - coor[4];
    y31 = coor[9] - coor[7];
    x32 = coor[6] - coor[5];
    y32 = coor[9] - coor[8];
    delta = x21 * y31 - x31 * y21;
    delt1 = delta / 6.;

/* --  FOMEGA: INTEGRATION AUX 3 SOMMETS (F(0,0)+F(1,0)+F(0,1)) */

    for (i__ = 1; i__ <= 2; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    ij = (j - 1 << 1) + i__;
	    be[ij] = fomega[i__ + (j << 1)] * delt1;
/* L1: */
	}
    }

/* --  FGAMMA: INTEGRATION AUX 2 EXTREMITES DES COTES (F(0)+F(1))/2. */

    if (norefs[7] != 0) {
	xji = coor[5] - coor[4];
	yji = coor[8] - coor[7];
/* Computing 2nd power */
	d__1 = xji;
/* Computing 2nd power */
	d__2 = yji;
	arelon = sqrt(d__1 * d__1 + d__2 * d__2);
	xnu = yji / arelon;
	ynu = -xji / arelon;
    }
    delt1 = sqrt(x21 * x21 + y21 * y21) / 2.;
    be[1] += delt1 * (fgamma[3] - pressi[1] * xnu);
    be[3] += delt1 * (fgamma[5] - pressi[2] * xnu);
    be[2] += delt1 * (fgamma[4] - pressi[1] * ynu);
    be[4] += delt1 * (fgamma[6] - pressi[2] * ynu);
    if (norefs[8] != 0) {
	xji = coor[6] - coor[5];
	yji = coor[9] - coor[8];
/* Computing 2nd power */
	d__1 = xji;
/* Computing 2nd power */
	d__2 = yji;
	arelon = sqrt(d__1 * d__1 + d__2 * d__2);
	xnu = yji / arelon;
	ynu = -xji / arelon;
    }
    delt1 = sqrt(x32 * x32 + y32 * y32) / 2.;
    be[3] += delt1 * (fgamma[7] - pressi[3] * xnu);
    be[5] += delt1 * (fgamma[9] - pressi[4] * xnu);
    be[4] += delt1 * (fgamma[8] - pressi[3] * ynu);
    be[6] += delt1 * (fgamma[10] - pressi[4] * ynu);
    if (norefs[9] != 0) {
	xji = coor[4] - coor[6];
	yji = coor[7] - coor[9];
/* Computing 2nd power */
	d__1 = xji;
/* Computing 2nd power */
	d__2 = yji;
	arelon = sqrt(d__1 * d__1 + d__2 * d__2);
	xnu = yji / arelon;
	ynu = -xji / arelon;
    }
    delt1 = sqrt(x31 * x31 + y31 * y31) / 2.;
    be[1] += delt1 * (fgamma[13] - pressi[6] * xnu);
    be[5] += delt1 * (fgamma[11] - pressi[5] * xnu);
    be[2] += delt1 * (fgamma[14] - pressi[6] * ynu);
    be[6] += delt1 * (fgamma[12] - pressi[5] * ynu);

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

/*     CALCUL DE LA MATRICE TENSEUR D'ELASTICITE * TENSEUR THERMIQUE */

    a[0] = alpha[1] * e[0] + alpha[3] * e[1] + alpha[2] * e[3];
    a[1] = alpha[1] * e[1] + alpha[3] * e[2] + alpha[2] * e[4];
    a[2] = alpha[1] * e[3] + alpha[3] * e[4] + alpha[2] * e[5];

    somme = (theta[1] + theta[2] + theta[3]) / (float)6.;
    be[1] += (x32 * a[2] - y32 * a[0]) * somme;
    be[2] += (x32 * a[1] - y32 * a[2]) * somme;
    be[3] += (-x31 * a[2] + y31 * a[0]) * somme;
    be[4] += (-x31 * a[1] + y31 * a[2]) * somme;
    be[5] += (x21 * a[2] - y21 * a[0]) * somme;
    be[6] += (x21 * a[1] - y21 * a[2]) * somme;

} /* ets2p1d_ */

