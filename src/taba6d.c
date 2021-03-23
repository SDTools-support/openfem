/* taba6d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int taba6d_(i1, i2, a, b, c__)
int32 *i1, *i2;
doublereal *a, *b, *c__;
{
    /* System generated locals */
    int32 a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static int32 i__, j, k, l;
    static doublereal s;

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                  S.P. TABA6D */
/*                  ---------------- */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* BUT : CALCUL DU PRODUIT MATRICIEL TRANSPOSEE(A)*B*A */
/* ----  C = TA * B * A */
/*       C(I1,I1) B(I2,I2) A(I2,I1) */

/* PARAMETRES D ENTREE : */
/* ----------------------- */
/* I1  : NBRE DE COLONNES DE A */
/* I2  : NBRE DE LIGNES DE A */
/* A   : MATRICE I2,I1 */
/* B   : MATRICE I2,I2 */

/* PARAMETRE RESULTAT : */
/* ----------------------- */
/* C   : C=TRANSPOSEE(A)*B*A */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* PROGRAMMEUR T. CAMPOS LAN 189 PARIS SEPTEMBRE 1979 */
/* ...................................................................... */

    /* Parameter adjustments */
    c_dim1 = *i1;
    c_offset = c_dim1 + 1;
    c__ -= c_offset;
    b_dim1 = *i2;
    b_offset = b_dim1 + 1;
    b -= b_offset;
    a_dim1 = *i2;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    i__1 = *i1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *i1;
	for (j = 1; j <= i__2; ++j) {
	    s = 0.;
	    i__3 = *i2;
	    for (l = 1; l <= i__3; ++l) {
		i__4 = *i2;
		for (k = 1; k <= i__4; ++k) {
		    s += a[k + i__ * a_dim1] * b[k + l * b_dim1] * a[l + j * a_dim1];
/* L4: */
		}
/* L3: */
	    }
	    c__[i__ + j * c_dim1] = s;
/* L2: */
	}
/* L1: */
    }
} /* taba6d_ */

