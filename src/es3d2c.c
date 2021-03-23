/* es3d2c.f -- translated by f2c (version 19961017).
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

/* Subroutine */ int es3d2c_(nno, npo, nbface, npoq, nnoq, npop, nnop, noref, npi, poir, npisq, npisp, poiqs, poips, nloc, x, y, z__, fom, fgam, pn, pr, pr1, dpr, dpr1, dpqs1, dpps1, pqs, pqs1, pps, pps1, delts, delta, be)
int32 *nno, *npo, *nbface, *npoq, *nnoq, *npop, *nnop, *noref, *npi;
doublereal *poir;
int32 *npisq, *npisp;
doublereal *poiqs, *poips;
int32 *nloc;
doublereal *x, *y, *z__, *fom, *fgam, *pn, *pr, *pr1, *dpr, *dpr1, *dpqs1, *dpps1, *pqs, *pqs1, *pps, *pps1, *delts, *delta, *be;
{
    /* System generated locals */
    int32 fgam_dim2, fgam_offset, pn_dim1, pn_offset, pr_dim1, pr_offset, dpr_dim2, dpr_offset, pr1_dim1, pr1_offset, dpr1_dim2, dpr1_offset, pqs_dim1, pqs_offset, pqs1_dim1, pqs1_offset, dpqs1_dim2, dpqs1_offset, pps_dim1, pps_offset, pps1_dim1, pps1_offset, dpps1_dim2, dpps1_offset, noref_dim1, noref_offset, nloc_dim1, nloc_offset, i__1, i__2, i__3;

    /* Local variables */
    static doublereal dpps[1], dpqs[1];
    static int32 iopt;
    static doublereal surf[27]	/* was [3][9] */, xint[1], yint[1], zint[1];
    static int32 i__, j, l;
    static doublereal dfinv[9]	/* was [3][3] */;
    extern /* Subroutine */ int dcopy_();
    static doublereal a2[1], df[9]	/* was [3][3] */;
    static int32 nf, indice;
    extern /* Subroutine */ int fobase_();
    static doublereal xs[9], ys[9], zs[9], ff1[1], ff2[1], ff3[1], vol[81]	/* was [3][27] */;

/*  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     BUT :   CALCUL DU SECOND MEMBRE ELEMENTAIRE */
/*     ---     D UN ELEMENT  PENTAEDRE */
/*  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     PARAMETRES D ENTREE */
/*     ------------------- */

/*     NNO           : NOMBRE DE NOEUDS DE L ELEMENT */
/*     NPO           : NOMBRE DE POINTS DE L ELEMENT */
/*     NBFACE        : NOMBRE DE FACES DE L ELEMENT */
/*     nloc : numerotation locale des faces */

/*     --- volume -- */
/*     nno ,npo           : nombre de noeuds , points */
/*     nbface             : nombre de faces */
/*     poir(npi)          : poids */
/*     pr,dpr             : interpolation */
/*     pr1,dpr1           : geometrie */
/*     --- face triangulaire -- */
/*     nnop ,npop         : nombre de noeuds , points */
/*     poips(npisp)       : poids */
/*     pps,dpps           : interpolation */
/*     pps1,dpps1         : geometrie */

/*     --- face quadrangulaire -- */
/*     nnop ,npop         : nombre de noeuds , points */
/*     poiqs(npisq)       : poids */
/*     pqs,dpqs           : interpolation */
/*     pqs1,dpqs1         : geometrie */
/*     FOM                : efforts volumique aux noeuds */
/*     FGAM               : esfforts surfaciques aux noeuds */
/*     PN                 : pression aux noeuds */
/*     NOREF              : noref(nbface,2) 1 =/ 0 si effort surfacique */
/*                                          2 =/ 0 si pression */
/*     PARAMETRES RESULTATS */
/*     -------------------- */
/*     BE            : SECOND MEMBRE ELEMENTAIRE */
/*  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     PROGRAMMEUR   : Marina Vidrascu INRIA 2001 */
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
    pn_dim1 = *nnoq;
    pn_offset = pn_dim1 + 1;
    pn -= pn_offset;
    fgam_dim2 = *nnoq;
    fgam_offset = (fgam_dim2 + 1) * 3 + 1;
    fgam -= fgam_offset;
    nloc_dim1 = *nnoq;
    nloc_offset = nloc_dim1 + 1;
    nloc -= nloc_offset;
    --delta;
    dpr1_dim2 = *npo;
    dpr1_offset = (dpr1_dim2 + 1) * 3 + 1;
    dpr1 -= dpr1_offset;
    dpr_dim2 = *nno;
    dpr_offset = (dpr_dim2 + 1) * 3 + 1;
    dpr -= dpr_offset;
    pr1_dim1 = *npo;
    pr1_offset = pr1_dim1 + 1;
    pr1 -= pr1_offset;
    pr_dim1 = *nno;
    pr_offset = pr_dim1 + 1;
    pr -= pr_offset;
    --poir;
    --delts;
    pqs1_dim1 = *npoq;
    pqs1_offset = pqs1_dim1 + 1;
    pqs1 -= pqs1_offset;
    pqs_dim1 = *nnoq;
    pqs_offset = pqs_dim1 + 1;
    pqs -= pqs_offset;
    dpqs1_dim2 = *npoq;
    dpqs1_offset = (dpqs1_dim2 + 1 << 1) + 1;
    dpqs1 -= dpqs1_offset;
    --poiqs;
    pps1_dim1 = *npop;
    pps1_offset = pps1_dim1 + 1;
    pps1 -= pps1_offset;
    pps_dim1 = *nnop;
    pps_offset = pps_dim1 + 1;
    pps -= pps_offset;
    dpps1_dim2 = *npop;
    dpps1_offset = (dpps1_dim2 + 1 << 1) + 1;
    dpps1 -= dpps1_offset;
    --poips;

    /* Function Body */
    indice = 3;
    i__1 = *npi * 3;
    dcopy_(&i__1, &c_b2, &c__0, vol, &c__1);
    fobase_(&c__3, &c__3, nno, npo, npi, &pr[pr_offset], &pr1[pr1_offset], &dpr[dpr_offset], &dpr1[dpr1_offset], &x[1], &y[1], &z__[1], a2, xint, yint, zint, &delta[1], dfinv, df, &indice, vol, &fom[4], &fom[4], &c__1);
/*     TERMES DE VOLUME */
/*     ---------------- */
    i__1 = *nno * 3;
    dcopy_(&i__1, &c_b2, &c__0, &be[4], &c__1);
    i__1 = *nno;
    for (j = 1; j <= i__1; ++j) {
	for (i__ = 1; i__ <= 3; ++i__) {
	    i__2 = *npi;
	    for (l = 1; l <= i__2; ++l) {
		be[i__ + j * 3] += pr[j + l * pr_dim1] * poir[l] * delta[l] * vol[i__ + l * 3 - 4];
/* L13: */
	    }
	}
    }

/*     TERMES DE BORD */
/*    ----------------- */
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

	    if (nf == 1 || nf == 4) {
/*              FACE TRIANGULAIRE */
		i__1 = *npop;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    xs[i__ - 1] = x[nloc[i__ + nf * nloc_dim1]];
		    ys[i__ - 1] = y[nloc[i__ + nf * nloc_dim1]];
		    zs[i__ - 1] = z__[nloc[i__ + nf * nloc_dim1]];
/* L116: */
		}
		indice = 3;
		i__1 = *npisp * 3;
		dcopy_(&i__1, &c_b2, &c__0, surf, &c__1);
		fobase_(&c__3, &c__2, nnop, npop, npisp, &pps[pps_offset], &pps1[pps1_offset], dpps, &dpps1[dpps1_offset], xs, ys, zs, a2, ff1, ff2, ff3, &delts[1], dfinv, df, &indice, surf, &fgam[(nf * fgam_dim2 + 1) * 3 + 1], &pn[nf * pn_dim1 + 1], &iopt);

		i__1 = *nnop;
		for (j = 1; j <= i__1; ++j) {
		    for (i__ = 1; i__ <= 3; ++i__) {
			i__3 = *npisp;
			for (l = 1; l <= i__3; ++l) {
			    be[i__ + nloc[j + nf * nloc_dim1] * 3] += pps[j + l * pps_dim1] * delts[l] * poips[l] * surf[i__ + l * 3 - 4];
/* L125: */
			}
		    }
		}
	    } else {
/*              FACE QUADRANGULAIRE */
		i__3 = *npoq;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    xs[i__ - 1] = x[nloc[i__ + nf * nloc_dim1]];
		    ys[i__ - 1] = y[nloc[i__ + nf * nloc_dim1]];
		    zs[i__ - 1] = z__[nloc[i__ + nf * nloc_dim1]];
/* L115: */
		}
		indice = 3;
		i__3 = *npisq * 3;
		dcopy_(&i__3, &c_b2, &c__0, surf, &c__1);
		fobase_(&c__3, &c__2, nnoq, npoq, npisq, &pqs[pqs_offset], &pqs1[pqs1_offset], dpqs, &dpqs1[dpqs1_offset], xs, ys, zs, a2, ff1, ff2, ff3, &delts[1], dfinv, df, &indice, surf, &fgam[(nf * fgam_dim2 + 1) * 3 + 1], &pn[nf * pn_dim1 + 1], &iopt);

		i__3 = *nnoq;
		for (j = 1; j <= i__3; ++j) {
		    for (i__ = 1; i__ <= 3; ++i__) {
			i__1 = *npisq;
			for (l = 1; l <= i__1; ++l) {
			    be[i__ + nloc[j + nf * nloc_dim1] * 3] += pqs[j + l * pqs_dim1] * delts[l] * poiqs[l] * surf[i__ + l * 3 - 4];
/* L25: */
			}
		    }
		}
	    }
	}
/* L26: */
    }
} /* es3d2c_ */

