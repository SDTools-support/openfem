/* tab2d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int tab2d_(i1, i2, a, b, c__)
int32 *i1, *i2;
doublereal *a, *b, *c__;
{
    /* System generated locals */
    int32 a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3;

    /* Local variables */
    static int32 i__, j, k;
    static doublereal s;
    static int32 j1, jc, jk, jjj;

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*             S.P  TAB2D */
/*             ---------- */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* BUT:CALCUL DU PRODUIT MATRICIEL TRANSPOSEE(A)*B DANS LE CAS OU */
/* ---  B EST SYMETRIQUE STOCKEE SOUS FORME TRIANGULAIRE */
/*      [ C ] = 0 + [ A ] * [ B ] */
/*        I1*I2       I1*I2     S */
/* PARAMETRES D ENTREE */
/* -------------------- */
/*  I1   : NOMBRE DE COLONNES  DE A */
/*  I2   : NOMBRE DE LIGNES DE A ,DE LIGNES ET COLONNES DE B */
/*  A    : MATRICE I1,I2 */
/*  B    : MATRICE SYMETRIQUE STOCKEE DE HAUT EN BAS ET */
/*        DE LA GAUCHE VERS LA DROITE.LA PARTIE TRIANGULAIRE */
/*        SUPERIEURE EST STOCKEE */

/* PARAMETRES DE SORTIE */
/* -------------------- */
/*  C    : PRODUIT MATRICIEL A*B  C(I1,I2) */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  PROGRAMMEUR:MARINA VIDRASCU  INRIA  FEVRIER 1980 */
/*              ALAIN  PERRONNET LAN189 PARIS ET INRIA  MARS 1980 */
/* ...................................................................... */

    /* Parameter adjustments */
    c_dim1 = *i1;
    c_offset = c_dim1 + 1;
    c__ -= c_offset;
    a_dim1 = *i2;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --b;

    /* Function Body */
    jc = 0;
    i__1 = *i2;
    for (j = 1; j <= i__1; ++j) {
	jjj = jc + j;
	j1 = j + 1;
	i__2 = *i1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    s = 0.;
	    i__3 = j;
	    for (k = 1; k <= i__3; ++k) {
		s += b[jc + k] * a[k + i__ * a_dim1];
/* L3: */
	    }

	    if (j == *i2) {
		goto L5;
	    }

	    jk = jjj;
	    i__3 = *i2;
	    for (k = j1; k <= i__3; ++k) {
		jk = jk + k - 1;
		s += b[jk] * a[k + i__ * a_dim1];
/* L4: */
	    }
L5:
	    c__[i__ + j * c_dim1] = s;
/* L2: */
	}
	jc = jjj;
/* L1: */
    }

} /* tab2d_ */

