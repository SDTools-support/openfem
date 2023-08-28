/* ets2p2c.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int ets2p2c_(coor, fomega, fgamma, pressi, norefs, alpha, theta, car, iopt, be)
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
    static doublereal aret[18]	/* was [2][9] */, a[3], c__, e[6];
    static int32 i__, j, k, m;
    static doublereal delta[6], desur[3], young, unmnu, x21, y21, x31, y31, x32, y32, x41, y41, x42, y42, x54, y54, x61, y61, x63, y63, x65, y65, x52, y52, arelon, x53, y53, xmi, ymi, xjm, yjm, xnu[3], ynu[3], poisson;

/* ................................................................... */
/* but: SECONDS MEMBRES ELEMENTAIRES TRIA 2P2C */
/* ---  INTERPOLATION LAGRANGE P2 SUR UN TRIANGLE ISOPARAMETRIQUE */

/* in : coor(6,ndim) : coordonnees des 6 noeuds */
/*      fomega(ndim,noe)  : fx, fy volumiques aux 6 noeuds */
/*      fgamma(ndim,3*nbarete): fx, fy aux 3 noeuds de chaque arete. */
/*      pressi(3*nbarete)     : pression aux 3 noeuds de chaque arete */
/*      norefs(nbarete,2): norefs(i,1) = 0 si fgamma   = 0 sur arete_i */
/*                         norefs(i,2) = 0 si pression = 0 sur arete_i */
/*   theta(6): temperature aux 6 noeuds */
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

/* out: BE(12) */

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
    coor -= 7;

    /* Function Body */
    x21 = coor[8] - coor[7];
    y21 = coor[14] - coor[13];
    x31 = coor[9] - coor[7];
    y31 = coor[15] - coor[13];
    x32 = coor[9] - coor[8];
    y32 = coor[15] - coor[14];
    x41 = coor[10] - coor[7];
    y41 = coor[16] - coor[13];
    x42 = coor[10] - coor[8];
    y42 = coor[16] - coor[14];
    x54 = coor[11] - coor[10];
    y54 = coor[17] - coor[16];
    x61 = coor[12] - coor[7];
    y61 = coor[18] - coor[13];
    x63 = coor[12] - coor[9];
    y63 = coor[18] - coor[15];
    x65 = coor[12] - coor[11];
    y65 = coor[18] - coor[17];

/* --  INTEGRATION : TERMES SURFACIQUES  ---- */

/* --  CALCUL DE DELTA AUX NOEUDS  -- */

    delta[0] = (x41 * 3. + x42) * (y61 * 3. + y63) - (x61 * 3. + x63) * (y41 * 3. + y42);
    delta[1] = (x42 * 3. + x41) * (y54 * -4. + y31) + (x54 * 4. - x31) * (y42 * 3. + y41);
    delta[2] = (x65 * 4. + x21) * (y63 * 3. + y61) - (x63 * 3. + x61) * (y65 * 4. + y21);
    delta[3] = x21 * (y54 * 2. + y61 + y63) - y21 * (x54 * 2. + x61 + x63);
    delta[4] = -(x65 * 2. + x41 + x42) * (y54 * 2. - y61 - y63) + (y65 * 2. + y41 + y42) * (x54 * 2. - x61 - x63);
    delta[5] = y31 * (x65 * -2. + x41 + x42) + x31 * (y65 * 2. - y41 - y42);
    for (j = 1; j <= 6; ++j) {
	delta[j - 1] /= 6.;
/* L1: */
    }

    for (j = 1; j <= 2; ++j) {
	be[j] = delta[0] * fomega[j + 2] / 10. - (delta[1] * fomega[j + 4] + delta[2] * fomega[j + 6]) / 60. - delta[4] * fomega[j + 10] / 15.;
	be[j + 2] = delta[1] * fomega[j + 4] / 10. - (delta[0] * fomega[j + 2] + delta[2] * fomega[j + 6]) / 60. - delta[5] * fomega[j + 12] / 15.;
	be[j + 4] = delta[2] * fomega[j + 6] / 10. - (delta[0] * fomega[j + 2] + delta[1] * fomega[j + 4]) / 60. - delta[3] * fomega[j + 8] / 15.;
	be[j + 6] = delta[3] * 8. * fomega[j + 8] / 15. + (delta[4] * fomega[j + 10] + delta[5] * fomega[j + 12]) * 4. / 15. - delta[2] * fomega[j + 6] / 15.;
	be[j + 8] = delta[4] * 8. * fomega[j + 10] / 15. + (delta[5] * fomega[j + 12] + delta[3] * fomega[j + 8]) * 4. / 15. - delta[0] * fomega[j + 2] / 15.;
	be[j + 10] = delta[5] * 8. * fomega[j + 12] / 15. + (delta[3] * fomega[j + 8] + delta[4] * fomega[j + 10]) * 4. / 15. - delta[1] * fomega[j + 4] / 15.;
/* L2: */
    }

/* --  FGAMMA: INTEGRATION aux 3 noeuds de l'arete. */
/*     fgamma(2,3*nbarete) fx, fy aux 3 noeuds de chaque arete */
/*     pressi(3*nbarete) pression aux 3 noeuds de chaque arete */
/*     aret(2,3*nbarete) fx, fy aux 3 noeuds de chaque arete */

/*     -- arete 1 */
    if (norefs[7] != 0) {
	k = 1;
	j = 2;
	m = k + 3;
	xmi = coor[m + 6] - coor[k + 6];
	ymi = coor[m + 12] - coor[k + 12];
/* Computing 2nd power */
	d__1 = xmi;
/* Computing 2nd power */
	d__2 = ymi;
	arelon = sqrt(d__1 * d__1 + d__2 * d__2);
	xnu[0] = ymi / arelon;
	ynu[0] = -xmi / arelon;
	xjm = coor[j + 6] - coor[m + 6];
	yjm = coor[j + 12] - coor[m + 12];
/* Computing 2nd power */
	d__1 = xjm;
/* Computing 2nd power */
	d__2 = yjm;
	arelon = sqrt(d__1 * d__1 + d__2 * d__2);
	xnu[1] = yjm / arelon;
	ynu[1] = -xjm / arelon;
	xnu[2] = (xnu[0] + xnu[1]) * .5;
	ynu[2] = (ynu[0] + ynu[1]) * .5;
	aret[0] = fgamma[3] - pressi[1] * xnu[0];
	aret[1] = fgamma[4] - pressi[1] * ynu[0];
	aret[2] = fgamma[5] - pressi[2] * xnu[1];
	aret[3] = fgamma[6] - pressi[2] * ynu[1];
	aret[4] = fgamma[7] - pressi[3] * xnu[2];
	aret[5] = fgamma[8] - pressi[3] * ynu[2];
    } else {
	aret[0] = fgamma[3];
	aret[1] = fgamma[4];
	aret[2] = fgamma[5];
	aret[3] = fgamma[6];
	aret[4] = fgamma[7];
	aret[5] = fgamma[8];
    }
    desur[0] = sqrt((x41 * 3. + x42) * (x41 * 3. + x42) + (y41 * 3. + y42) * (y41 * 3. + y42)) / 6.;
    desur[1] = sqrt((x42 * 3. + x41) * (x42 * 3. + x41) + (y42 * 3. + y41) * (y42 * 3. + y41)) / 6.;
    desur[2] = sqrt(x21 * x21 + y21 * y21) * 4. / 6.;

    be[1] += desur[0] * aret[0];
    be[2] += desur[0] * aret[1];
    be[3] += desur[1] * aret[2];
    be[4] += desur[1] * aret[3];
    be[7] += desur[2] * aret[4];
    be[8] += desur[2] * aret[5];
/*     -- arete 2 */
    if (norefs[8] != 0) {
	k = 2;
	j = 3;
	m = k + 3;
	xmi = coor[m + 6] - coor[k + 6];
	ymi = coor[m + 12] - coor[k + 12];
/* Computing 2nd power */
	d__1 = xmi;
/* Computing 2nd power */
	d__2 = ymi;
	arelon = sqrt(d__1 * d__1 + d__2 * d__2);
	xnu[0] = ymi / arelon;
	ynu[0] = -xmi / arelon;
	xjm = coor[j + 6] - coor[m + 6];
	yjm = coor[j + 12] - coor[m + 12];
/* Computing 2nd power */
	d__1 = xjm;
/* Computing 2nd power */
	d__2 = yjm;
	arelon = sqrt(d__1 * d__1 + d__2 * d__2);
	xnu[1] = yjm / arelon;
	ynu[1] = -xjm / arelon;
	xnu[2] = (xnu[0] + xnu[1]) * .5;
	ynu[2] = (ynu[0] + ynu[1]) * .5;
	aret[6] = fgamma[9] - pressi[4] * xnu[0];
	aret[7] = fgamma[10] - pressi[4] * ynu[0];
	aret[8] = fgamma[11] - pressi[5] * xnu[1];
	aret[9] = fgamma[12] - pressi[5] * ynu[1];
	aret[10] = fgamma[13] - pressi[6] * xnu[2];
	aret[11] = fgamma[14] - pressi[6] * ynu[2];
    } else {
	aret[6] = fgamma[9];
	aret[7] = fgamma[10];
	aret[8] = fgamma[11];
	aret[9] = fgamma[12];
	aret[10] = fgamma[13];
	aret[11] = fgamma[14];
    }
    x52 = coor[11] - coor[8];
    y52 = coor[17] - coor[14];
    x53 = coor[11] - coor[9];
    y53 = coor[17] - coor[15];
    desur[0] = sqrt((x52 * 3. + x53) * (x52 * 3. + x53) + (y52 * 3. + y53) * (y52 * 3. + y53)) / 6.;
    desur[1] = sqrt((x53 * 3. + x52) * (x53 * 3. + x52) + (y53 * 3. + y52) * (y53 * 3. + y52)) / 6.;
    desur[2] = sqrt(x32 * x32 + y32 * y32) * 4. / 6.;

    be[3] += desur[0] * aret[6];
    be[4] += desur[0] * aret[7];
    be[5] += desur[1] * aret[8];
    be[6] += desur[1] * aret[9];
    be[9] += desur[2] * aret[10];
    be[10] += desur[2] * aret[11];
/*     -- arete 3 */
    if (norefs[9] != 0) {
	k = 3;
	j = 1;
	m = k + 3;
	xmi = coor[m + 6] - coor[k + 6];
	ymi = coor[m + 12] - coor[k + 12];
/* Computing 2nd power */
	d__1 = xmi;
/* Computing 2nd power */
	d__2 = ymi;
	arelon = sqrt(d__1 * d__1 + d__2 * d__2);
	xnu[0] = ymi / arelon;
	ynu[0] = -xmi / arelon;
	xjm = coor[j + 6] - coor[m + 6];
	yjm = coor[j + 12] - coor[m + 12];
/* Computing 2nd power */
	d__1 = xjm;
/* Computing 2nd power */
	d__2 = yjm;
	arelon = sqrt(d__1 * d__1 + d__2 * d__2);
	xnu[1] = yjm / arelon;
	ynu[1] = -xjm / arelon;
	xnu[2] = (xnu[0] + xnu[1]) * .5;
	ynu[2] = (ynu[0] + ynu[1]) * .5;
	aret[12] = fgamma[15] - pressi[7] * xnu[0];
	aret[13] = fgamma[16] - pressi[7] * ynu[0];
	aret[14] = fgamma[17] - pressi[8] * xnu[1];
	aret[15] = fgamma[18] - pressi[8] * ynu[1];
	aret[16] = fgamma[19] - pressi[9] * xnu[2];
	aret[17] = fgamma[20] - pressi[9] * ynu[2];
    } else {
	aret[12] = fgamma[15];
	aret[13] = fgamma[16];
	aret[14] = fgamma[17];
	aret[15] = fgamma[18];
	aret[16] = fgamma[19];
	aret[17] = fgamma[20];
    }

    desur[0] = sqrt((x63 * 3. + x61) * (x63 * 3. + x61) + (y63 * 3. + y61) * (y63 * 3. + y61)) / 6.;
    desur[1] = sqrt((x61 * 3. + x63) * (x61 * 3. + x63) + (y61 * 3. + y63) * (y61 * 3. + y63)) / 6.;
    desur[2] = sqrt(x31 * x31 + y31 * y31) * 4. / 6.;

    be[5] += desur[0] * aret[12];
    be[6] += desur[0] * aret[13];
    be[1] += desur[1] * aret[14];
    be[2] += desur[1] * aret[15];
    be[11] += desur[2] * aret[16];
    be[12] += desur[2] * aret[17];

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

    be[1] += ((-a[0] * y32 + a[2] * x32) * theta[4] + (a[0] * y32 - a[2] * x32) * theta[5] + (-a[0] * y32 + a[2] * x32) * theta[6]) / (float)6.;
    be[2] += ((-a[2] * y32 + a[1] * x32) * theta[4] + (a[2] * y32 - a[1] * x32) * theta[5] + (-a[2] * y32 + a[1] * x32) * theta[6]) / (float)6.;
    be[3] += ((a[0] * y31 - a[2] * x31) * theta[4] + (a[0] * y31 - a[2] * x31) * theta[5] + (-a[0] * y31 + a[2] * x31) * theta[6]) / (float)6.;
    be[4] += ((a[2] * y31 - a[1] * x31) * theta[4] + (a[2] * y31 - a[1] * x31) * theta[5] + (-a[2] * y31 + a[1] * x31) * theta[6]) / (float)6.;
    be[5] += ((a[0] * y21 - a[2] * x21) * theta[4] + (-a[0] * y21 + a[2] * x21) * theta[5] + (-a[0] * y21 + a[2] * x21) * theta[6]) / (float)6.;
    be[6] += ((a[2] * y21 - a[1] * x21) * theta[4] + (-a[2] * y21 + a[1] * x21) * theta[5] + (-a[2] * y21 + a[1] * x21) * theta[6]) / (float)6.;
    be[7] += ((a[0] * y21 - a[2] * x21) * (float)2. * theta[4] + (-a[0] * y32 + a[2] * x32) * (float)2. * theta[5] + (a[0] * y31 - a[2] * x31) * (float)2. * theta[6]) / (float)6.;
    be[8] += ((a[2] * y21 - a[1] * x21) * (float)2. * theta[4] + (-a[2] * y32 + a[1] * x32) * (float)2. * theta[5] + (a[2] * y31 - a[1] * x31) * (float)2. * theta[6]) / (float)6.;
    be[9] += ((-a[0] * y21 + a[2] * x21) * (float)2. * theta[4] + (a[0] * y32 - a[2] * x32) * (float)2. * theta[5] + (a[0] * y31 - a[2] * x31) * (float)2. * theta[6]) / (float)6.;
    be[10] += ((-a[2] * y21 + a[1] * x21) * (float)2. * theta[4] + (a[2] * y32 - a[1] * x32) * (float)2. * theta[5] + (a[2] * y31 - a[1] * x31) * (float)2. * theta[6]) / (float)6.;
    be[11] += ((-a[0] * y21 + a[2] * x21) * (float)2. * theta[4] + (-a[0] * y32 + a[2] * x32) * (float)2. * theta[5] + (-a[0] * y31 + a[2] * x31) * (float)2. * theta[6]) / (float)6.;
    be[12] += ((-a[2] * y21 + a[1] * x21) * (float)2. * theta[4] + (-a[2] * y32 + a[1] * x32) * (float)2. * theta[5] + (-a[2] * y31 + a[1] * x31) * (float)2. * theta[6]) / (float)6.;
} /* ets2p2c_ */

