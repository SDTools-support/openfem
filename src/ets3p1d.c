/* ets3p1d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__4 = 4;
static int32 c__3 = 3;

/* Subroutine */ int ets3p1d_(coor, volf, sur, pres, noref, be)
doublereal *coor, *volf, *sur, *pres;
int32 *noref;
doublereal *be;
{
    /* Initialized data */

    static int32 nloc[12]	/* was [3][4] */ = { 1,3,2,1,4,3,1,2,4,2,3,4 };
    static doublereal dp[48]	/* was [3][4][4] */ = { -1.,-1.,-1.,1.,0.,0.,0.,1.,0.,0.,0.,1.,-1.,-1.,-1.,1.,0.,0.,0.,1.,0.,0.,0.,1.,-1.,-1.,-1.,1.,0.,0.,0.,1.,0.,0.,0.,1.,-1.,-1.,-1.,1.,0.,0.,0.,1.,0.,0.,0.,1. };
    static doublereal vp1[16]	/* was [4][4] */ = { 1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1. };
    static doublereal dps1[18]	/* was [2][3][3] */ = { -1.,-1.,1.,0.,0.,1.,-1.,-1.,1.,0.,0.,1.,-1.,-1.,1.,0.,0.,1. };
    static doublereal vps1[9]	/* was [3][3] */ = { 1.,0.,0.,0.,1.,0.,0.,0.,1. };
    static doublereal poids[4] = { .04166666666666666,.04166666666666666,.04166666666666666,.04166666666666666 };
    static doublereal poisms[3] = { .1666666666666667,.1666666666666667,.1666666666666667 };

    extern /* Subroutine */ int es3c2c_();
    static doublereal delta[4], delts[3];

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                  S.P. ES3P1D */
/*                  ------------ */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  BUT   : CALCUL DU SECOND MEMBRE ELEMENTAIRE ELASTIQUE */
/*  ---      TETRAEDRE P1 DROIT  . */

/*  PARAMETRES D ENTREE  : */
/*  ------------------- */
/*   X,Y,Z   : COORDONNEES DES POINTS DE L ELEMENT (SIMPLE PRECISION) */
/*   VOLF    : FORCES VOLUMIQUES */
/*   NOREF   : noref(nbface,2) avec nbface = nombre de faces */
/*             noref(i,1) =/ 0 si face i soumise force surfacique */
/*             noref(2,1) =/ 0 si face i soumise a pression */
/*   SUR     : efforts surfaciques */
/*   pres    : pression */
/*  PARAMETRE DE SORTIE  : */
/*  -------------------- */
/*   BE      : TABLEAU DES SECONDS MEMBRES */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  PROGRAMMEUR  : Marina Vidrascu 2001 */
/* ................................................................... */

    /* Parameter adjustments */
    --be;
    noref -= 5;
    pres -= 4;
    sur -= 13;
    volf -= 4;
    coor -= 5;

    /* Function Body */

/*         poids = .25/6. */
/*         poisms = .5/3 */

    es3c2c_(&c__4, &c__4, &c__4, &c__3, &c__3, &noref[5], &c__4, poids, &c__3, poisms, nloc, &coor[5], &coor[9], &coor[13], &volf[4], &sur[13], &pres[4], vp1, vp1, dp, dp, dps1, vps1, vps1, &be[1], delts, delta);
} /* ets3p1d_ */

