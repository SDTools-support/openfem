/* dpaq2c.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int dpaq2c_(r__, z__, nbpoly, npi, poidsg, p, dp, f1, f2, dfm1dp, poidel)
doublereal *r__, *z__;
int32 *nbpoly, *npi;
doublereal *poidsg, *p, *dp, *f1, *f2, *dfm1dp, *poidel;
{
    /* System generated locals */
    int32 p_dim1, p_offset, dp_dim2, dp_offset, i__1, i__2;

    /* Local variables */
    static doublereal d__;
    static int32 i__, l;
    static doublereal rr, zz, df1, df2, df3, df4, dfm1[36]	/* was [4][9] */;

/*  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* IN:  R et Z des 9 points de l'element + (Cf. include polbas.ins) */
/*      NBPOLY : NBRE DE POLYNOMES DE BASE = 9 */
/*      NPI    : NBRE DE POINTS D INTEGRATION SUR L ELEMENT = 9 */
/*      POIDSG(NPI): VALEUR DES NPI POIDS GAUSS */
/*      P(NBPOLY,NPI): VALEUR des polynomes de base aux points d'int. num */
/*      DP(2,NBPOLY,NPI): Valeur DERIVEE DP/DX,  DP/DY */
/* OUT: F1,F2  COORDONNEES DES PTS D INT. NUMERIQUES DE L ELEMENT COURANT */
/*      DFM1DP = PRODUIT DF-1 * DP */
/*      POIDEL : PRODUIT DU POIDSG*DELTA*R AUX POINTS D'INTEGRATION */
/* ...................................................................... */

/*      F1,DFM1,DELTA */

    /* Parameter adjustments */
    --r__;
    --z__;
    dp_dim2 = *nbpoly;
    dp_offset = (dp_dim2 + 1 << 1) + 1;
    dp -= dp_offset;
    p_dim1 = *nbpoly;
    p_offset = p_dim1 + 1;
    p -= p_offset;
    --poidsg;
    --f1;
    --f2;
    dfm1dp -= 21;
    --poidel;

    /* Function Body */
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
	    rr += p[i__ + l * p_dim1] * r__[i__];
	    zz += p[i__ + l * p_dim1] * z__[i__];
	    df1 += dp[(i__ + l * dp_dim2 << 1) + 1] * r__[i__];
	    df2 += dp[(i__ + l * dp_dim2 << 1) + 2] * r__[i__];
	    df3 += dp[(i__ + l * dp_dim2 << 1) + 1] * z__[i__];
	    df4 += dp[(i__ + l * dp_dim2 << 1) + 2] * z__[i__];
/* L3: */
	}
	f1[l] = rr;
	f2[l] = zz;
	d__ = df1 * df4 - df2 * df3;
	dfm1[(l << 2) - 4] = df4 / d__;
	dfm1[(l << 2) - 3] = -df2 / d__;
	dfm1[(l << 2) - 2] = -df3 / d__;
	dfm1[(l << 2) - 1] = df1 / d__;
	poidel[l] = abs(d__) * poidsg[l] * rr;

/*        DFM1DP = DF-1 * DP SUR T UNITE */

	i__2 = *nbpoly;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dfm1dp[(i__ + l * 9 << 1) + 1] = dfm1[(l << 2) - 4] * dp[(i__ + l * dp_dim2 << 1) + 1] + dfm1[(l << 2) - 2] * dp[(i__ + l * dp_dim2 << 1) + 2];
	    dfm1dp[(i__ + l * 9 << 1) + 2] = dfm1[(l << 2) - 3] * dp[(i__ + l * dp_dim2 << 1) + 1] + dfm1[(l << 2) - 1] * dp[(i__ + l * dp_dim2 << 1) + 2];
/* L4: */
	}
/* L2: */
    }

/*     COORDONNEES DU BARYCENTRE DE L'ELEMENT COURANT:  F1(5) ET F2(5) */
} /* dpaq2c_ */

