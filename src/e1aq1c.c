/* e1aq1c.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int e1aq1c_(nbpoly, npi, poids, q13, dq13, ip, f1, f2, dfm1dp, poidel, dfm1, rz)
int32 *nbpoly, *npi;
doublereal *poids, *q13, *dq13;
int32 *ip;
doublereal *f1, *f2, *dfm1dp, *poidel, *dfm1, *rz;
{
    /* Initialized data */

    static doublereal d2pi = 6.2831853071795862;

    /* System generated locals */
    int32 q13_dim1, q13_offset, dq13_dim2, dq13_offset, rz_dim1, rz_offset, i__1, i__2;

    /* Local variables */
    static doublereal d__;
    static int32 i__, j, l, i2;
    static doublereal rr, zz, df1, df2, df3, df4;

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                    S.P. E1AP1D */
/*                   ------------ */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* BUT: IP     : PASSAGE DES D.L. PAR COMPOSANTE AUX NOEUDS */
/*      F1,F2  : COOR DES POINTS D INT.  NUM. DE L ELEMENT COURANT */
/*      DFM1   = (DF) -1 AUX POINTS D INTEGRATION */
/*      DFM1DP = DERIVEES DES POLYNOMES AUX POINTS D INTEGRATION */
/*      POIDEL : POIDS * DELTA * R  AUX POINTS D'INTEGRATION */

/* PARAMETRES D ENTREE : */
/* --------------------- */
/* NBPOLY : NBRE DE POLYNOMES DE BASE */
/* NPI    : NBRE DE POINTS D INTEGRATION SUR L ELEMENT */
/* POIDS  : VALEUR DES NPI POIDS */
/* Q13    : Q13(I,J)=VALEUR DE PI(POINT J D INTEGRATION) */
/* DQ13   : DQ13(I,J,L)=DERIVEE DPJ/DXI(POINT L) */
/* RZ(NOE,NDIM): RAYON ET COTE DES 3 POINTS DE L ELEMENT */
/* ...................................................................... */
/*  programmeur : modulef */
/* ...................................................................... */
    /* Parameter adjustments */
    rz_dim1 = *nbpoly;
    rz_offset = rz_dim1 + 1;
    rz -= rz_offset;
    dq13_dim2 = *nbpoly;
    dq13_offset = (dq13_dim2 + 1 << 1) + 1;
    dq13 -= dq13_offset;
    q13_dim1 = *nbpoly;
    q13_offset = q13_dim1 + 1;
    q13 -= q13_offset;
    --poids;
    --ip;
    --f1;
    --f2;
    dfm1dp -= 11;
    --poidel;
    dfm1 -= 5;

    /* Function Body */

/*     IP(I) = POSITION DU I-EME D.L. PAR COMPOSANTE DANS L ORDRE */
/*             NOEUD PAR NOEUD */

    ip[1] = 1;
    ip[8] = 8;
    j = 1;
    i2 = (*nbpoly << 1) - 2;
    i__1 = i2;
    for (i__ = 2; i__ <= i__1; i__ += 2) {
	ip[j + *nbpoly] = i__;
	++j;
	ip[j] = i__ + 1;
/* L1: */
    }

/*      F1,DFM1,DELTA */

    i__1 = *npi;
    for (l = 1; l <= i__1; ++l) {
	rr = 0.;
	zz = 0.;
	df1 = 0.;
	df2 = 0.;
	df3 = 0.;
	df4 = 0.;
	i__2 = *nbpoly;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    rr += q13[i__ + l * q13_dim1] * rz[i__ + rz_dim1];
	    zz += q13[i__ + l * q13_dim1] * rz[i__ + (rz_dim1 << 1)];
	    df1 += dq13[(i__ + l * dq13_dim2 << 1) + 1] * rz[i__ + rz_dim1];
	    df2 += dq13[(i__ + l * dq13_dim2 << 1) + 2] * rz[i__ + rz_dim1];
	    df3 += dq13[(i__ + l * dq13_dim2 << 1) + 1] * rz[i__ + (rz_dim1 << 1)];
	    df4 += dq13[(i__ + l * dq13_dim2 << 1) + 2] * rz[i__ + (rz_dim1 << 1)];
/* L3: */
	}
	f1[l] = rr;
	f2[l] = zz;
	d__ = df1 * df4 - df2 * df3;
	dfm1[(l << 2) + 1] = df4 / d__;
	dfm1[(l << 2) + 2] = -df2 / d__;
	dfm1[(l << 2) + 3] = -df3 / d__;
	dfm1[(l << 2) + 4] = df1 / d__;
	poidel[l] = abs(d__) * poids[l] * d2pi * rr;

/*        DFM1DP = DF-1 * DP SUR T UNITE */

	i__2 = *nbpoly;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dfm1dp[(i__ + (l << 2) << 1) + 1] = dfm1[(l << 2) + 1] * dq13[(i__ + l * dq13_dim2 << 1) + 1] + dfm1[(l << 2) + 3] * dq13[(i__ + l * dq13_dim2 << 1) + 2];
	    dfm1dp[(i__ + (l << 2) << 1) + 2] = dfm1[(l << 2) + 2] * dq13[(i__ + l * dq13_dim2 << 1) + 1] + dfm1[(l << 2) + 4] * dq13[(i__ + l * dq13_dim2 << 1) + 2];
/* L4: */
	}
/* L2: */
    }

/*     COORDONNEES DU BARYCENTRE DE L'ELEMENT COURANT */

    rr = 0.;
    zz = 0.;
    l = *npi + 1;
    i__1 = *nbpoly;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rr += rz[i__ + rz_dim1];
	zz += rz[i__ + (rz_dim1 << 1)];
/* L5: */
    }
    f1[l] = rr / 4.;
    f2[l] = zz / 4.;
} /* e1aq1c_ */

