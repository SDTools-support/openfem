/* es3c2c.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static doublereal c_b2 = 0.;
static int32 c__0 = 0;
static int32 c__1 = 1;
static int32 c__3 = 3;
static int32 c__2 = 2;

/* Subroutine */ int es3c2c_(nno, npo, nbface, npof, nnof, noref, npi, poids, npis, poisms, nloc, x, y, z__, fom, fgam, pn, vp, vp1, vdpq2, vdpq1, vdpsq1, vps, vps1, be, delts, delta)
int32 *nno, *npo, *nbface, *npof, *nnof, *noref, *npi;
doublereal *poids;
int32 *npis;
doublereal *poisms;
int32 *nloc;
doublereal *x, *y, *z__, *fom, *fgam, *pn, *vp, *vp1, *vdpq2, *vdpq1, *vdpsq1, *vps, *vps1, *be, *delts, *delta;
{
    /* System generated locals */
    int32 fgam_dim2, fgam_offset, pn_dim1, pn_offset, vp_dim1, vp_offset, vp1_dim1, vp1_offset, vdpq2_dim2, vdpq2_offset, vdpq1_dim2, vdpq1_offset, vdpsq1_dim2, vdpsq1_offset, vps_dim1, vps_offset, vps1_dim1, vps1_offset, noref_dim1, noref_offset, nloc_dim1, nloc_offset, i__1, i__2, i__3;

    /* Local variables */
    static int32 iopt;
    static doublereal xint[1], yint[1], zint[1];
    static int32 i__, j, l;
    static doublereal dfinv[9]	/* was [3][3] */;
    extern /* Subroutine */ int dcopy_();
    static doublereal a2[1], df[9]	/* was [3][3] */;
    static int32 nf, indice;
    extern /* Subroutine */ int fobase_();
    static doublereal xs[9], ys[9], zs[9], ff1[1], ff2[1], ff3[1], vol[81]	/* was [3][27] */, sur[27]	/* was [3][9] */;

/*  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     BUT :   CALCUL DU SECOND MEMBRE ELEMENTAIRE */
/*     ---     D UN ELEMENT  HEXAEDRIQUE OU TETRAEDRIQUE */
/*  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     PARAMETRES D ENTREE */
/*     ------------------- */
/*     NNO           : NOMBRE DE NOEUDS DE L ELEMENT */
/*     NPO           : NOMBRE DE POINTS DE L ELEMENT */
/*     NBFACE        : NOMBRE DE FACES DE L ELEMENT */
/*     NNOF,NPOF     : NOMBRE DE NOEUDS (POINTS) PAR FACE */
/*     NPI           : NOMBRE DE POINTS DE LA FORMULE D INTEGRATION */
/*                     SUR LE CUBE DE REFERENCE */
/*     NLOC          : NLOC(NNOF,NBFACE) */
/*     POIDS         : POIDS DE LA FORMULE D INTEGRATION */
/*     NPIS POISMS   : MEME CHOSE POUR LES FACES */
/*     X,Y,Z         : COORDONNEES DES NOEUDS DE L ELEMENT */
/*     VOL           : DONNEES VOLUMIQUES */
/*     VP            : VALEURS DES POLYNOMES DE BASE FONCTIONS */
/*                     AUX POINTS D INTEGRATION */
/*     VP1           : VALEUR DES POLYNOMES DE BASE (GEOMETRIE - FT) */
/*     VDPQ2         : VALEURS DES DERIVES DES POL DE BASE FONCTIONS */
/*                     AUX POINTS D INTEGRATION */
/*     VDPQ1         : VALEURS DES DERIVES DES POLYNOMES DE BASE GEOM */
/*                     AUX POINTS D INTEGRATION */
/*     DELTA         : JACOBIEN  VOLUME */
/*     DELTS         : JACOBIEN FACES */
/*     FOM           : efforts volumique aux noeuds */
/*     FGAM          : esfforts surfaciques aux noeuds */
/*     PN            : pression aux noeuds */
/*     NOREF         : noref(nbface,2) 1 =/ 0 si effort surfacique */
/*                                     2 =/ 0 si pression */
/*     PARAMETRES RESULTATS */
/*     -------------------- */
/*     BE            : SECOND MEMBRE ELEMENTAIRE */
/*  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     PROGRAMMEUR  :  M VIDRASCU 2001 */
/*  ................................................................... */
/*     --- non utilises --- */

    /* Parameter adjustments */
    be -= 4;
    fom -= 4;
    --z__;
    --y;
    --x;
    noref_dim1 = *nbface;
    noref_offset = noref_dim1 + 1;
    noref -= noref_offset;
    pn_dim1 = *nnof;
    pn_offset = pn_dim1 + 1;
    pn -= pn_offset;
    fgam_dim2 = *nnof;
    fgam_offset = (fgam_dim2 + 1) * 3 + 1;
    fgam -= fgam_offset;
    nloc_dim1 = *nnof;
    nloc_offset = nloc_dim1 + 1;
    nloc -= nloc_offset;
    --delta;
    vdpq1_dim2 = *npo;
    vdpq1_offset = (vdpq1_dim2 + 1) * 3 + 1;
    vdpq1 -= vdpq1_offset;
    vdpq2_dim2 = *nno;
    vdpq2_offset = (vdpq2_dim2 + 1) * 3 + 1;
    vdpq2 -= vdpq2_offset;
    vp1_dim1 = *npo;
    vp1_offset = vp1_dim1 + 1;
    vp1 -= vp1_offset;
    vp_dim1 = *nno;
    vp_offset = vp_dim1 + 1;
    vp -= vp_offset;
    --poids;
    --delts;
    vps1_dim1 = *npof;
    vps1_offset = vps1_dim1 + 1;
    vps1 -= vps1_offset;
    vps_dim1 = *nnof;
    vps_offset = vps_dim1 + 1;
    vps -= vps_offset;
    vdpsq1_dim2 = *npof;
    vdpsq1_offset = (vdpsq1_dim2 + 1 << 1) + 1;
    vdpsq1 -= vdpsq1_offset;
    --poisms;

    /* Function Body */
    indice = 3;
    i__1 = *npi * 3;
    dcopy_(&i__1, &c_b2, &c__0, vol, &c__1);
    fobase_(&c__3, &c__3, nno, npo, npi, &vp[vp_offset], &vp1[vp1_offset], &vdpq2[vdpq2_offset], &vdpq1[vdpq1_offset], &x[1], &y[1], &z__[1], a2, xint, yint, zint, &delta[1], dfinv, df, &indice, vol, &fom[4], &fom[4], &c__1);
/*     TERMES DE VOLUME */
/*     ---------------- */

    i__1 = *nno * 3;
    dcopy_(&i__1, &c_b2, &c__0, &be[4], &c__1);
    i__1 = *nno;
    for (j = 1; j <= i__1; ++j) {
	for (i__ = 1; i__ <= 3; ++i__) {
	    i__2 = *npi;
	    for (l = 1; l <= i__2; ++l) {
		be[i__ + j * 3] += vp[j + l * vp_dim1] * poids[l] * delta[l] * vol[i__ + l * 3 - 4];
/* L12: */
	    }
	}
    }

/*    TERMES DE BORDS */

    i__2 = *nbface;
    for (nf = 1; nf <= i__2; ++nf) {
	if (noref[nf + noref_dim1] > 0 || noref[nf + (noref_dim1 << 1)] > 0) {

/*       calcul de iopt */

	    if (noref[nf + noref_dim1] != 0) {
		if (noref[nf + (noref_dim1 << 1)] != 0) {
		    iopt = 3;
		} else {
		    iopt = 1;
		}
	    } else {
		iopt = 2;
	    }
	    i__1 = *npof;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		xs[i__ - 1] = x[nloc[i__ + nf * nloc_dim1]];
		ys[i__ - 1] = y[nloc[i__ + nf * nloc_dim1]];
		zs[i__ - 1] = z__[nloc[i__ + nf * nloc_dim1]];
/* L115: */
	    }
	    indice = 3;
	    i__1 = *npis * 3;
	    dcopy_(&i__1, &c_b2, &c__0, sur, &c__1);
	    fobase_(&c__3, &c__2, nnof, npof, npis, &vps[vps_offset], &vps1[vps1_offset], &vdpsq1[vdpsq1_offset], &vdpsq1[vdpsq1_offset], xs, ys, zs, a2, ff1, ff2, ff3, &delts[1], dfinv, df, &indice, sur, &fgam[(nf * fgam_dim2 + 1) * 3 + 1], &pn[nf * pn_dim1 + 1], &iopt);
	    i__1 = *nnof;
	    for (j = 1; j <= i__1; ++j) {
		for (i__ = 1; i__ <= 3; ++i__) {
		    i__3 = *npis;
		    for (l = 1; l <= i__3; ++l) {
			be[i__ + nloc[j + nf * nloc_dim1] * 3] += vps[j + l * vps_dim1] * delts[l] * poisms[l] * sur[i__ + l * 3 - 4];
/* L23: */
		    }
		}
	    }
	}
/* L26: */
    }
} /* es3c2c_ */

