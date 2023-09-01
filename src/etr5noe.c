/* etr5noe.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int etr5noe_(coor, car, iopt, ae)
doublereal *coor, *car;
int32 *iopt;
doublereal *ae;
{
    /* Initialized data */

    static int32 n[12]	/* was [3][4] */ = { 1,2,5,2,3,5,3,4,5,4,1,5 };

    /* System generated locals */
    int32 i__1;

    /* Local variables */
    static int32 i__, j, k;
    static doublereal s, x[5], y[5];
    static logical icomp;
    extern /* Subroutine */ int mexErrMsgTxt();
    static doublereal e1;
    static int32 i1, i2;
    static doublereal young, b11, b12, d11[10], d22[10], d12[10];
    static int32 ii, ij;
    static doublereal x31, y31;
    static int32 np;
    static doublereal x42, y42, rr, ss, bi1[3], bi2[3], poison, div[10], div1;

/*  .................................................................... */
/* but : CREATION DES MATRICES ELEMENTAIRES RAIDEUR POUR L ELEMENT */
/* ---   QUADRILATERAL A 4 SOUS TRIANGLES EN CROIX */
/*              ---   VERSION PLANE  --- */
/* ON DISTINGUE DEUX CAS SUIVANT QUE poison EST PROCHE DE .5 EN */
/* DEFORMATIONS PLANES (CAS INCOMPRESSIBLE) OU NON (CAS COMPRESSIBLE) */
/* 1) CAS COMPRESSIBLE */
/*    A(I,J) = 2.*E1 * D(I).D(J) */
/*           + RR*2.* DIV(D(I))*DIV(D(J)) */
/* 2) CAS INCOMPRESSIBLE */
/*    A(I,J) = 2.*E1 *( (D(I)-.5*TR(D(I))*ID) * (D(J)-.5*TR(D(J))*ID ) */
/*           + RR*2. * DIV(D(I))*DIV(D(J)) */

/* in : coor(noe,ndim) : coordones des 4 sommets. */
/*      car(2): caracteristiques des materiaux */
/*              car(1) = young */
/*              car(2) = poisson */
/*      iopt = 1 isotrope Contraintes  Planes */
/*           = 2 isotrope Deformations Planes */
/* out: ae(55)    : DEMI-MATRICE SUPERIEURE */
/* .................................................................... */
    /* Parameter adjustments */
    --ae;
    --car;
    coor -= 6;

    /* Function Body */

    young = car[1];
    poison = car[2];
    e1 = young * .5 / (poison + 1.);
    icomp = TRUE_;
    if (*iopt == 2) {
	ss = 1. - poison * 2.;
	if (abs(ss) < (float).1) {
	    icomp = FALSE_;
	}
	rr = e1 * .5 / ss;
	if (icomp) {
	    rr = rr * 2. * poison;
	}
    } else if (*iopt == 1) {
	rr = e1 * poison / (1. - poison);
    } else {
	mexErrMsgTxt("CHOIX CONTRAINTES OU DEFORMATIONS PLANES ?", 42L);
	return 0;
    }

/*     CALCUL DES COORDONNEES DES NOEUDS */

    for (i__ = 1; i__ <= 4; ++i__) {
	x[i__ - 1] = coor[i__ + 5];
	y[i__ - 1] = coor[i__ + 10];
/* L1: */
    }
    x31 = x[2] - x[0];
    y31 = y[2] - y[0];
    x42 = x[3] - x[1];
    y42 = y[3] - y[1];
    div1 = (float)1. / (y31 * x42 - y42 * x31);
    x[4] = (x[0] * y31 * x42 - x[1] * y42 * x31 + (y[1] - y[0]) * x31 * x42) * div1;
    y[4] = (y[1] * y31 * x42 - y[0] * y42 * x31 - (x[1] - x[0]) * y31 * y42) * div1;

/*     INITIALISATION */

    for (ij = 1; ij <= 55; ++ij) {
	ae[ij] = (float)0.;
/* L2: */
    }

/*     BOUCLE SUR LES SOUS-ELEMENTS */

    for (np = 1; np <= 4; ++np) {
	for (i1 = 1; i1 <= 10; ++i1) {
	    div[i1 - 1] = 0.;
	    d11[i1 - 1] = 0.;
	    d12[i1 - 1] = 0.;
	    d22[i1 - 1] = 0.;
/* L100: */
	}
	i__ = n[np * 3 - 3];
	j = n[np * 3 - 2];
	k = n[np * 3 - 1];

/*        CALCUL DES GRADIENTS BIJ DES COORDONNEES BARYCENTRIQUES */
/*        ET DE S=2*AIRE DU SOUS ELEMENT */

	s = (x[i__ - 1] - x[j - 1]) * (y[i__ - 1] - y[k - 1]) - (x[i__ - 1] - x[k - 1]) * (y[i__ - 1] - y[j - 1]);
	ss = 1. / s;
	bi1[0] = (y[j - 1] - y[k - 1]) * ss;
	bi1[1] = (y[k - 1] - y[i__ - 1]) * ss;
	bi1[2] = (y[i__ - 1] - y[j - 1]) * ss;
	bi2[0] = (x[k - 1] - x[j - 1]) * ss;
	bi2[1] = (x[i__ - 1] - x[k - 1]) * ss;
	bi2[2] = (x[j - 1] - x[i__ - 1]) * ss;

/*        ACTUALISATION DE LA MATRICE */
	for (ii = 1; ii <= 3; ++ii) {
	    i1 = (n[ii + np * 3 - 4] << 1) - 1;
	    i2 = i1 + 1;
	    b11 = bi1[ii - 1] * .5;
	    b12 = bi2[ii - 1] * .5;
/*             D11(I1)=B11 */
	    d12[i1 - 1] = b12;
/*             D22(I1)=-B11 */
	    div[i1 - 1] = b11 * 2.;
/*             D11(I2)=-B12 */
	    d12[i2 - 1] = b11;
/*             D22(I2)=B12 */
	    div[i2 - 1] = b12 * 2.;
	    if (icomp) {
		d11[i1 - 1] = b11 * 2.;
		d22[i2 - 1] = b12 * 2.;
	    } else {
		d11[i1 - 1] = b11;
		d22[i1 - 1] = -b11;
		d11[i2 - 1] = -b12;
		d22[i2 - 1] = b12;
	    }
/* L101: */
	}
	for (i__ = 1; i__ <= 10; ++i__) {
	    i__1 = i__;
	    for (j = 1; j <= i__1; ++j) {
		ij = (i__ - 1) * i__ / 2 + j;
		ae[ij] = ae[ij] + e1 * s * (d11[i__ - 1] * d11[j - 1] + d12[i__ - 1] * 2. * d12[j - 1] + d22[i__ - 1] * d22[j - 1]) + rr * s * div[i__ - 1] * div[j - 1];
/* L102: */
	    }
	}
/* L10: */
    }
} /* etr5noe_ */

