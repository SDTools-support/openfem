/* tab0d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int tab0d_(i1, i2, i3, a, b, c__)
int32 *i1, *i2, *i3;
doublereal *a, *b, *c__;
{
    /* System generated locals */
    int32 a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, i__3;

    /* Local variables */
    static int32 i__, j, k;
    static doublereal s;

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                        S.P. TAB0D */
/*                        ---------- */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* BUT : C = TA * B OU  C(I1,I3),A(I2,I1),B(I2,I3) */
/* ---- */
/* PARAMETRES D ENTREE : */
/* --------------------- */
/* I1  : NOMBRE DE LIGNES DE C ET DE COLONNES DEA */
/* I2  : NOMBRE DE LIGNES DE A ET DE B */
/* I3  : NOMBRE DE COLONNES DE C ET DE B */
/* A   : MATRICE A */
/* B   : MATRICE B */

/* PARAMETRE RESULTAT : */
/* -------------------- */
/* C   : MATRICE TA * B */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* PROGRAMMEUR : A.PERRONNET LAN189 ET INRIA PARIS  FEVRIER 1980 */
/* ...................................................................... */

    /* Parameter adjustments */
    a_dim1 = *i2;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    c_dim1 = *i1;
    c_offset = c_dim1 + 1;
    c__ -= c_offset;
    b_dim1 = *i2;
    b_offset = b_dim1 + 1;
    b -= b_offset;

    /* Function Body */
    i__1 = *i3;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *i1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    s = 0.;
	    i__3 = *i2;
	    for (k = 1; k <= i__3; ++k) {
		s += a[k + i__ * a_dim1] * b[k + j * b_dim1];
/* L3: */
	    }
	    c__[i__ + j * c_dim1] = s;
/* L2: */
	}
/* L1: */
    }
} /* tab0d_ */

