/* ab4d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int ab4d_(i1, i2, a, b, c__)
int32 *i1, *i2;
doublereal *a, *b, *c__;
{
    /* System generated locals */
    int32 a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;

    /* Local variables */
    static int32 j, l, m, n;
    static doublereal s;

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                   S.P. AB4D */
/*                   --------- */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* BUT : CALCUL DU PRODUIT MATRICIEL [ C ] = [ A ] * [ B ] */
/* -----                                 S     I1*I2   I2*I1 */
/* PARAMETRES D ENTREE : */
/* ---------------------- */
/* I1     : NBRE DE LIGNES DE A, ET DE COLONNES DE B */
/* I2     : NBRE DE COLONNES DE A,ET DE LIGNES DE B */
/* A , B  : MATRICES [ A ] , [ B ] */
/*                     I1*I2   I2*I1 */
/* PARAMETRE RESULTAT : */
/* ------------------ */
/* C   : MATRICE [ C ] = [ A ] * [ B ] . [ C ] EST SYMETRIQUE STOCKEE */
/*       DE HAUT EN BAS ET DE LA GAUCHE VERS LA DROITE. SEULE LA PARTIE */
/*       TRIANGULAIRE SUPERIEURE EST STOCKEE. */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* PROGRAMMEUR : A. HASSIM INRIA */
/* .................................................................... */

    /* Parameter adjustments */
    b_dim1 = *i2;
    b_offset = b_dim1 + 1;
    b -= b_offset;
    a_dim1 = *i1;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --c__;

    /* Function Body */
    l = 0;
    i__1 = *i1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	for (m = 1; m <= i__2; ++m) {
	    s = 0.;
	    i__3 = *i2;
	    for (n = 1; n <= i__3; ++n) {
		s += a[m + n * a_dim1] * b[n + j * b_dim1];
/* L2: */
	    }
	    ++l;
	    c__[l] = s;
/* L1: */
	}
    }
} /* ab4d_ */

