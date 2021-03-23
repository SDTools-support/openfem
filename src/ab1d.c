/* ab1d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int ab1d_(i__, k, j, a, b, c__)
int32 *i__, *k, *j;
doublereal *a, *b, *c__;
{
    /* System generated locals */
    int32 a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, i__3;

    /* Local variables */
    static int32 l, m, n;
    static doublereal s;

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                   S.P. AB1D */
/*                  ---------------- */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* BUT : CALCUL DU PRODUIT MATRICIEL C = C + A * B */
/* ----- */
/* PARAMETRES D ENTREE : */
/* ---------------------- */
/* I   : NBRE DE LIGNES DE A */
/* J   : NBRE DE COLONNES DE B */
/* K   : NBRE DE COLONNES DE A,ET DE LIGNES DE B */
/* A   : MATRICE I,K */
/* B   : MATRICE K,J */

/* PARAMETRE RESULTAT : */
/* ---------------------- */
/* C   : C = C + A * B */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* PROGRAMMEUR T. CAMPOS LAN 189 PARIS SEPTEMBRE 1979 */
/* ...................................................................... */

    /* Parameter adjustments */
    a_dim1 = *i__;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    c_dim1 = *i__;
    c_offset = c_dim1 + 1;
    c__ -= c_offset;
    b_dim1 = *k;
    b_offset = b_dim1 + 1;
    b -= b_offset;

    /* Function Body */
    i__1 = *i__;
    for (l = 1; l <= i__1; ++l) {
	i__2 = *j;
	for (m = 1; m <= i__2; ++m) {
	    s = c__[l + m * c_dim1];
	    i__3 = *k;
	    for (n = 1; n <= i__3; ++n) {
		s += a[l + n * a_dim1] * b[n + m * b_dim1];
/* L3: */
	    }
	    c__[l + m * c_dim1] = s;
/* L2: */
	}
/* L1: */
    }
} /* ab1d_ */

