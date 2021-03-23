/* taba8d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int taba8d_(i1, i2, a, b, c__, aux)
int32 *i1, *i2;
doublereal *a, *b, *c__, *aux;
{
    /* System generated locals */
    int32 a_dim1, a_offset, aux_dim1, aux_offset, i__1, i__2, i__3;

    /* Local variables */
    static int32 i__, j, k, n;
    static doublereal s;
    extern /* Subroutine */ int tabszd_();

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*             S.P  TABA8D */
/*             ---------- */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* BUT:CALCUL DU PRODUIT MATRICIEL TRANSPOSE(A)*B*A DANS LE CAS OU */
/* ---  B EST SYMETRIQUE STOCKEE SOUS FORME TRIANGULAIRE */
/*      C=TRANSPOSEE(A)*B*A */

/* PARAMETRES D ENTREE: */
/* -------------------- */
/*      I2   :NOMBRE DE LIGNES  DE A , DE LIGNES ET COLONNES DE B */
/*      I1   :NOMBRE DE COLONNES DE A */
/*      A    :MATRICE I2,I1 */
/*      B    :MATRICE (I2,I2) SYMETRIQUE STOCKEE DE HAUT EN BAS ET */
/*            DE LA GAUCHE VERS LA DROITE.LA PARTIE TRIANGULAIRE */
/*            SUPERIEURE EST STOCKEE */
/*      AUX  :MATRICE AUXILIAIRE I1 I2 .CONTIENT LE PRODUIT */
/*             TRANSPOSEE(A)*B.PERMET DE CALCULER LE PRODUIT FINAL */
/*             EN 2*O(N**3) OPERATIONS */

/* PARAMETRES DE SORTIE */
/* -------------------- */
/*       C   : PRODUIT MATRICIEL TRANSPOSE(A)*B*A  C(I1,I1) */
/*             C EST SYMETRIQUE STOCKEE DE HAUT EN BAS ET DE LA GAUCHE */
/*             VERS LA DROITE . SEULE LA PARTIE TRIANGULAIRE */
/*             SUPERIEURE EST STOCKEE */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  PROGRAMMEUR:MARINA VIDRASCU  INRIA  FEVRIER 1980 */
/* ...................................................................... */
/* ...................................................................... */


/*     AUX=T(A)*B */
/*     ---------- */
    /* Parameter adjustments */
    aux_dim1 = *i1;
    aux_offset = aux_dim1 + 1;
    aux -= aux_offset;
    a_dim1 = *i2;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --b;
    --c__;

    /* Function Body */
    tabszd_(i1, i2, &a[a_offset], &b[1], &aux[aux_offset]);

/*     C=AUX*A */
/*     ------ */
    n = 0;
    i__1 = *i1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    s = 0.;
	    i__3 = *i2;
	    for (k = 1; k <= i__3; ++k) {
		s += aux[i__ + k * aux_dim1] * a[k + j * a_dim1];
/* L2: */
	    }
	    ++n;
	    c__[n] = s;
/* L1: */
	}
    }
} /* taba8d_ */

