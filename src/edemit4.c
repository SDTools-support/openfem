/* edemit4.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__0 = 0;
static int32 c__1 = 1;
static int32 c__2 = 2;

/* Subroutine */ int edemit4_(pol)
doublereal *pol;
{
    /* Initialized data */

    static doublereal xty[4] = { 0.,0.,-1.,1. };
    static doublereal yty[4] = { -1.,1.,0.,0. };
    static doublereal xint[4] = { -.5773502691896258,.5773502691896258,.5773502691896258,-.5773502691896258 };
    static doublereal yint[4] = { -.5773502691896258,-.5773502691896258,.5773502691896258,.5773502691896258 };

    static int32 l, n;
    extern doublereal pcq12d_();
    static int32 lt;

/* -----*--*---------*---------*---------*---------*---------*---------*-- */

/* BUT : CALCUL DES POLYNOMES DE BASE ET DERIVEES POUR MIT4 */
/* --- */
/*       Le tableau pol de 96 variables doit etre declare ds matlab et rempli */
/*       si l'elem mitc4 est utilise */
/* -----*--*---------*---------*---------*---------*---------*---------*-- */
/*     pol contient en sequence vp2(nno,npi+npt)   cad    vp2(4,4+4) */
/*                              vdpq2(2,nno,npi)   cad    vdpq2(2,4,4) */
/*                              vdtq2(2,nno,npt)   cad    vdtq2(2,4,4) */



/*     tying points */

    /* Parameter adjustments */
    --pol;

    /* Function Body */

/*     integration points en r,s */


/*        remplir les tableaux vdpq2(2,nno,npi),vdtq2(2,nno,npt), */
/*                             vp2(nno,npi+npt) */

    for (l = 1; l <= 4; ++l) {
	for (n = 1; n <= 4; ++n) {
/*           -- > vp2(n,l)     = pcq12d(0,n,xint(l),yint(l)) */
	    pol[(l - 1 << 2) + n] = pcq12d_(&c__0, &n, &xint[l - 1], &yint[l - 1]);
/*           -- > vdpq2(1,n,l) = pcq12d(1,n,xint(l),yint(l)) */
	    pol[(l - 1 << 3) + 33 + (n - 1 << 1)] = pcq12d_(&c__1, &n, &xint[l - 1], &yint[l - 1]);
/*           -- > vdpq2(2,n,l) = pcq12d(2,n,xint(l),yint(l)) */
	    pol[(l - 1 << 3) + 34 + (n - 1 << 1)] = pcq12d_(&c__2, &n, &xint[l - 1], &yint[l - 1]);
/* L1: */
	}
    }
    for (lt = 1; lt <= 4; ++lt) {
	for (n = 1; n <= 4; ++n) {
/*           -- >  vp2(n,npi+lt) = pcq12d(0,n,xty(lt),yty(lt)) */
	    pol[(lt + 3 << 2) + n] = pcq12d_(&c__0, &n, &xty[lt - 1], &yty[lt - 1]);
/*           -- >  vdtq2(1,n,lt) = pcq12d(1,n,xty(lt),yty(lt)) */
	    pol[(lt - 1 << 3) + 65 + (n - 1 << 1)] = pcq12d_(&c__1, &n, &xty[lt - 1], &yty[lt - 1]);
/*           -- >  vdtq2(2,n,lt) = pcq12d(2,n,xty(lt),yty(lt)) */
	    pol[(lt - 1 << 3) + 66 + (n - 1 << 1)] = pcq12d_(&c__2, &n, &xty[lt - 1], &yty[lt - 1]);
/* L2: */
	}
    }
} /* edemit4_ */

