/* etsmit4.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__20 = 20;
static int32 c__1 = 1;
static doublereal c_b6 = 0.;
static int32 c__0 = 0;
static int32 c__4 = 4;
static int32 c__3 = 3;
static int32 c__24 = 24;

/* Subroutine */ int etsmit4_(coor, sur, iopt, pol, xnorm, be)
doublereal *coor, *sur;
int32 *iopt;
doublereal *pol, *xnorm, *be;
{
    /* Initialized data */

    static doublereal poids[4] = { 1.,1.,1.,1. };

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    extern doublereal ddot_();
    static doublereal suri[20]	/* was [5][4] */;
    static int32 i__, k, l, n;
    static doublereal delta;
    extern /* Subroutine */ int dcopy_();
    static doublereal v1[12]	/* was [3][4] */, v2[12]	/* was [3][4] */;
    extern /* Subroutine */ int replo1_();
    static doublereal gr[12]	/* was [3][4] */, gs[12]	/* was [3][4] */;
    static int32 ichoix;
    static doublereal be5[20]	/* was [5][4] */, be6[24]	/* was [6][4] */, grr, grs, gss;
    extern /* Subroutine */ int etsmitx_();

/* -----*--*---------*---------*---------*---------*---------*---------*-- */

/* BUT : CALCUL DU SECOND MEMBRE DE L ELEMENT MIT4 */
/* --- */
/* -----*--*---------*---------*---------*---------*---------*---------*-- */
/* parametres d entree : */
/* ------------------- */
/* coor   : coor(npo,ndim) coordonnees */
/* SUR    : sur(5,nno) efforts volumiques aux noeuds */
/* iopt   : options */
/*          iopt(1) ; pas utilise */
/*          iopt(2) = ichoix >0 regulier */
/*                    -i e=x_i */
/*          iopt(3) = 5 si 5ddl par noud */
/*                    6 si 6 ddl par noud */

/* parametres resultats : */
/* -------------------- */
/* BE     : second memebre */
/* -----*--*---------*---------*---------*---------*---------*---------*-- */


/*     integration points en r,s */

    /* Parameter adjustments */
    --be;
    xnorm -= 4;
    --pol;
    --iopt;
    sur -= 6;
    coor -= 5;

    /* Function Body */
    if (FALSE_) {
	dcopy_(&c__20, &sur[6], &c__1, &be[1], &c__1);
	return 0;
    }
    dcopy_(&c__20, &c_b6, &c__0, be5, &c__1);
/*     *---------*---------*---------*---------*---------*---------*-- */
    etsmitx_(gr, gs, &sur[6], suri, &c__4, &c__4, &pol[33], &pol[1], &coor[5], &coor[9], &coor[13]);
/*     *---------*---------*---------*---------*---------*---------*-- */

/*     calcul du second membre */
/*     -------------------------------- */
/*     Somme(omega) f*vds */


/*     boucle sur les points d'integration (r,s) */
    for (l = 1; l <= 4; ++l) {

/*        calcul de grr gss grs comme polynomes de degre 2 en z */

	grr = ddot_(&c__3, &gr[l * 3 - 3], &c__1, &gr[l * 3 - 3], &c__1);
	gss = ddot_(&c__3, &gs[l * 3 - 3], &c__1, &gs[l * 3 - 3], &c__1);
	grs = ddot_(&c__3, &gr[l * 3 - 3], &c__1, &gs[l * 3 - 3], &c__1);
	delta = sqrt(grr * gss - grs * grs);
	for (n = 1; n <= 4; ++n) {
	    for (k = 1; k <= 5; ++k) {
/*          --- > be(k,n) = +vp2(n,l)*delta*suri(k,l)*poids(l) */
		be5[k + n * 5 - 6] += pol[(l - 1 << 2) + n] * delta * suri[k + l * 5 - 6] * poids[l - 1];
/* L2: */
	    }
	}
    }
    if (iopt[3] == 5) {

/*        --- > 5 ddl/noeud */
	dcopy_(&c__20, be5, &c__1, &be[1], &c__1);
    } else {

/*        --- > 6 ddl/noeud */
	ichoix = iopt[2];
	replo1_(&c__4, &xnorm[4], v1, v2, &ichoix);
	for (i__ = 1; i__ <= 4; ++i__) {
	    be6[i__ * 6 - 6] = be5[i__ * 5 - 5];
	    be6[i__ * 6 - 5] = be5[i__ * 5 - 4];
	    be6[i__ * 6 - 4] = be5[i__ * 5 - 3];
	    be6[i__ * 6 - 3] = v1[i__ * 3 - 3] * be5[i__ * 5 - 2] + v2[i__ * 3 - 3] * be5[i__ * 5 - 1];
	    be6[i__ * 6 - 2] = v1[i__ * 3 - 2] * be5[i__ * 5 - 2] + v2[i__ * 3 - 2] * be5[i__ * 5 - 1];
	    be6[i__ * 6 - 1] = v1[i__ * 3 - 1] * be5[i__ * 5 - 2] + v2[i__ * 3 - 1] * be5[i__ * 5 - 1];
/* L1: */
	}
	dcopy_(&c__24, be6, &c__1, &be[1], &c__1);
    }
} /* etsmit4_ */

