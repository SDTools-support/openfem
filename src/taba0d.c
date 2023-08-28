/* taba0d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int taba0d_(i1, i2, a, b, c__)
int32 *i1, *i2;
doublereal *a, *b, *c__;
{
    /* System generated locals */
    int32 a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static int32 i__, j, k, l, n;
    static doublereal s;
    static int32 l1, n1;
    static doublereal s1, s2, s3;

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                  S.P. TABA0D */
/*                  ----------- */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* BUT : CALCUL DU PRODUIT MATRICIEL TRANSPOSEE(A)*B*A DANS LE CAS */
/* ----- OU B EST SYMETRIQUE STOCKEE SOUS FORME TRIANGULAIRE */
/*       OU C EST SYMETRIQUE STOCKEE SOUS FORME TRIANGULAIRE */
/*       OU C = TA * B * A */

/* PARAMETRES D ENTREE : */
/* ----------------------- */
/* I1  : NBRE DE COLONNES DE A */
/* I2  : NBRE DE LIGNES DE A */
/* A   : MATRICE I2,I1 */
/* B   : MATRICE SYMETRIQUE STOCKEE DE HAUT EN BAS ET DE LA GAUCHE VE */
/*       VERS LA DROITE SEULE LA PARTIE TRIANGULAIRE SUPERIEURE EST */
/*       STOCKEE */

/* PARAMETRE RESULTAT : */
/* ----------------------- */
/* C   : PRODUIT MATRICIEL TRANSPOSEE(A)*B*A */
/*       DANS LE TABLEAU C SEULE LA PARTIE TRIANGULAIRE SUPERIEURE */
/*       DE HAUT EN BAS ET DE LA GAUCHE VERS LA DROITE EST STOCKEE */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* PROGRAMMEUR : A.PERRONNET LAN 189 PARIS ET IRIA OCTOBRE 1979 */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

    /* Parameter adjustments */
    a_dim1 = *i2;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --b;
    --c__;

    /* Function Body */
    n = 0;
    i__1 = *i1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ++n;
	    n1 = 0;
	    s = 0.;
	    i__3 = *i2;
	    for (l = 1; l <= i__3; ++l) {
		l1 = l - 1;
		s1 = a[l + j * a_dim1];
		s2 = a[l + i__ * a_dim1];
		if (l1 == 0) {
		    goto L5;
		}
		i__4 = l1;
		for (k = 1; k <= i__4; ++k) {
		    ++n1;
		    s3 = b[n1];
		    s += s3 * (a[k + i__ * a_dim1] * s1 + s2 * a[k + j * a_dim1]);
/* L4: */
		}
L5:
		++n1;
		s += s2 * b[n1] * s1;
/* L3: */
	    }
	    c__[n] = s;
/* L2: */
	}
/* L1: */
    }
} /* taba0d_ */

