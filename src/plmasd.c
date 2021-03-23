/* plmasd.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int plmasd_(i1, nbloc, gs, ijt, ae)
int32 *i1, *nbloc;
doublereal *gs;
int32 *ijt;
doublereal *ae;
{
    /* System generated locals */
    int32 i__1, i__2, i__3, i__4;

    /* Local variables */
    static int32 i__, j, k, n, j1, k1, k2, j2, ii, jj, jp, ii1, ii2, jj1, jj2, iii;

/*  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  BUT : PLONGER LA MATRICE TRIANGULAIRE SUPERIEURE GS DANS LES */
/*  ---   NBLOC BLOCS DIAGONAUX DE LA MATRICE AE ET INITIALISER LES */
/*        AUTRES A ZERO */
/*  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  PARAMETRES D'ENTREE : */
/*  ------------------- */
/*  I1    : NOMBRE DE LIGNES ET COLONNES DU BLOC DIAGONAL */
/*  NBLOC : NOMBRE DE BLOCS DIAGONAUX DE I1 LIGNES A PORTER DANS AE */
/*  GS    : LE BLOC DIAGONAL GS STOCKE 1 2 4  ... */
/*                                       3 5  ... */
/*                                         6  ... */
/*  IJT   : PERMUTATION DES D.L. PAR COMPOSANTE => PAR NOEUD */
/*  PARAMETRE RESULTAT : */
/*  ------------------ */
/*  AE     : MATRICE TRIANGULAIRE SUPERIEURE */
/* ....................................................................... */
/*  programmeur : modulef */
/* ............................................................... */
/*     ------   LE BLOC DE COLONNES JJ   ------ */
    /* Parameter adjustments */
    --ae;
    --ijt;
    --gs;

    /* Function Body */
    i__1 = *nbloc;
    for (jj = 1; jj <= i__1; ++jj) {
	jj1 = jj - 1;
	jj2 = jj1 * *i1;
/*        ------   LE BLOC DE LIGNES II   ------ */
	i__2 = jj;
	for (ii = 1; ii <= i__2; ++ii) {
	    n = 0;
	    ii1 = ii - 1;
	    ii2 = ii1 * *i1;
	    iii = *i1;
/*           ------   LA COLONNE J DU BLOC JJ   ------ */
	    i__3 = *i1;
	    for (j = 1; j <= i__3; ++j) {
		if (jj == ii) {
		    iii = j;
		}
		j1 = j + jj2;
		jp = ijt[j1];
/*              ------   LA LIGNE I DU BLOC II   ------ */
		i__4 = iii;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    k1 = i__ + ii2;
		    k2 = ijt[k1];
		    j2 = jp;
		    if (k2 <= j2) {
			goto L1;
		    }
		    k = k2;
		    k2 = j2;
		    j2 = k;
L1:
		    k = j2 * (j2 - 1) / 2 + k2;
		    if (jj == ii) {
			goto L2;
		    }
/*                 ------   BLOC NON DIAGONAL A ANNULER   ------ */
		    ae[k] = 0.;
		    goto L3;
/*                 ------   BLOC DIAGONAL   ------ */
L2:
		    ++n;
		    ae[k] = gs[n];
L3:
		    ;
		}
/* L4: */
	    }
/* L6: */
	}
    }
} /* plmasd_ */

