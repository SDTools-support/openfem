/* plonad.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int plonad_(nstoc, i1, i2, nbloci, nblocj, ijt, gs, gn, ae)
int32 *nstoc, *i1, *i2, *nbloci, *nblocj, *ijt;
doublereal *gs, *gn, *ae;
{
    /* System generated locals */
    int32 gn_dim1, gn_offset, i__1, i__2;

    /* Local variables */
    static int32 i__, j, k, ig, ii, jj, ni, nj, jjj;

/*  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  BUT : PLONGER DANS LA MATRICE AE SYMETRIQUE STOCKEE TRIANGULAIREMENT */
/*  ---   UNE MATRICE GN(I1,I2) SI NSTOC = -1 */
/*        UNE MATRICE GS(I1,I1) SI NSTOC =  1   STOCKEE TRIANGULAIREMENT */
/*        DE HAUT EN BAS ET DE LA GAUCHE VERS LA DROITE  1  2  4  ... */
/*                                                          3  5  ... */
/*  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  PARAMETRES D'ENTREE : */
/*  ------------------- */
/*  NSTOC  : 1 SI GS(1 A I1*(I1+1)/2) , -1 SI GN(I1,I2) */
/*  I1     : NOMBRE DE LIGNES DE GS(N) */
/*  I2     : NOMBRE DE COLONNES DE GN SI NSTOC = -1 */
/*  NBLOCI : NUMERO DU BLOC DE LIGNES I DANS LA MATRICE AE CORRESPONDANT */
/*           AUX I1 LIGNES DE GS(N) */
/*  NBLOCJ : NUMERO DU BLOC DE COLONNES DANS LA MATRICE AE CORRESPONDANT */
/*           AUX I2 COLONNES DE GN */
/*  IJT    : PERMUTATION DES D.L. PAR COMPOSANTE => PAR NOEUD */
/*  GS     : MATRICE POUR NSTOC = 1 */
/*  GN     : MATRICE POUR NSTOC =-1 */
/*  PARAMETRES DE SORTIE : */
/*  -------------------- */
/*  AE     : MATRICE TRIANGULAIRE SUPERIEURE SOMMEE A GS(N) */
/*  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  PROGRAMMEUR : A.PERRONNET LAN189 PARIS ET INRIA  FEVRIER 1980 */
/*  .................................................................... */

    /* Parameter adjustments */
    gn_dim1 = *i1;
    gn_offset = gn_dim1 + 1;
    gn -= gn_offset;
    --ijt;
    --gs;
    --ae;

    /* Function Body */
    ni = *i1 * (*nbloci - 1);
    nj = *i2 * (*nblocj - 1);
    if (*nstoc <= 0) {
	goto L5;
    } else {
	goto L1;
    }
/*     ------   GS SYMETRIQUE STOCKEE TRIANGULAIREMENT   ------ */
L1:
    ig = 0;
    i__1 = *i1;
    for (j = 1; j <= i__1; ++j) {
	jjj = ijt[j + nj];
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ii = ijt[i__ + ni];
	    jj = jjj;
	    if (ii <= jj) {
		goto L2;
	    }
	    k = ii;
	    ii = jj;
	    jj = k;
L2:
	    k = jj * (jj - 1) / 2 + ii;
	    ++ig;
	    ae[k] += gs[ig];
/* L3: */
	}
/* L4: */
    }
    return 0;
/*     ------   MATRICE NON SYMETRIQUE  GN(II,I2)   ------ */
L5:
    ig = *i1;
    i__1 = *i2;
    for (j = 1; j <= i__1; ++j) {
	if (*nbloci == *nblocj) {
	    ig = j;
	}
	jjj = ijt[j + nj];
	i__2 = ig;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ii = ijt[i__ + ni];
	    jj = jjj;
	    if (ii <= jj) {
		goto L6;
	    }
	    k = ii;
	    ii = jj;
	    jj = k;
L6:
	    k = jj * (jj - 1) / 2 + ii;
	    ae[k] += gn[i__ + j * gn_dim1];
/* L7: */
	}
/* L8: */
    }
} /* plonad_ */

