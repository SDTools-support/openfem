/* etr3r1d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__6 = 6;

/* Subroutine */ int etr3r1d_(coor, car, iopt, ae)
doublereal *coor, *car;
int32 *iopt;
doublereal *ae;
{
    /* Initialized data */

    static int32 ijt[18] = { 1,4,7,10,13,16,2,5,8,11,14,17,3,6,9,12,15,18 };
    static doublereal dp[108]	/* was [3][6][6] */ = { -.78867513459481,-.78867513459481,-.66666666666667,.78867513459481,0.,-.16666666666667,0.,.78867513459481,-.16666666666667,-.21132486540519,-.21132486540519,.66666666666667,.21132486540519,0.,.16666666666667,0.,.21132486540519,.16666666666667,-.78867513459481,-.78867513459481,-.16666666666667,.78867513459481,0.,-.66666666666667,0.,.78867513459481,-.16666666666667,-.21132486540519,-.21132486540519,.16666666666667,.21132486540519,0.,
	    .66666666666667,0.,.21132486540519,.16666666666667,-.78867513459481,-.78867513459481,-.16666666666667,.78867513459481,0.,-.16666666666667,0.,.78867513459481,-.66666666666667,-.21132486540519,-.21132486540519,.16666666666667,.21132486540519,0.,.16666666666667,0.,.21132486540519,.66666666666667,-.21132486540519,-.21132486540519,-.66666666666667,.21132486540519,0.,-.16666666666667,0.,.21132486540519,-.16666666666667,-.78867513459481,-.78867513459481,.66666666666667,.78867513459481,0.,
	    .16666666666667,0.,.78867513459481,.16666666666667,-.21132486540519,-.21132486540519,-.16666666666667,.21132486540519,0.,-.66666666666667,0.,.21132486540519,-.16666666666667,-.78867513459481,-.78867513459481,.16666666666667,.78867513459481,0.,.66666666666667,0.,.78867513459481,.16666666666667,-.21132486540519,-.21132486540519,-.16666666666667,.21132486540519,0.,-.16666666666667,0.,.21132486540519,-.66666666666667,-.78867513459481,-.78867513459481,.16666666666667,.78867513459481,0.,
	    .16666666666667,0.,.78867513459481,.66666666666667 };
    static doublereal vp1[36]	/* was [6][6] */ = { .52578342306321,.1314458557658,.1314458557658,.14088324360346,.035220810900865,.035220810900865,.1314458557658,.52578342306321,.1314458557658,.035220810900865,.14088324360346,.035220810900865,.1314458557658,.1314458557658,.52578342306321,.035220810900865,.035220810900865,.14088324360346,.14088324360346,.035220810900865,.035220810900865,.52578342306321,.1314458557658,.1314458557658,.035220810900865,.14088324360346,.035220810900865,
	    .1314458557658,.52578342306321,.1314458557658,.035220810900865,.035220810900865,.14088324360346,.1314458557658,.1314458557658,.52578342306321 };
    static doublereal poids[6] = { .083333333333333,.083333333333333,.083333333333333,.083333333333333,.083333333333333,.083333333333333 };

    static doublereal elas[9], xint[6], yint[6], zint[6];
    extern /* Subroutine */ int er3c2c_();
    static doublereal delta[6], a2[108]	/* was [3][6][6] */;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  BUT   : CALCUL DE LA MATRICE ELEMENTAIRE ELASTIQUE DE RIGIDITE */
/*  ---     pentaedre r1 DROIT */
/*          formule d'integration a 6 points (gauss) 3*2 */

/*  PARAMETRES D ENTREE  : */
/*  ------------------- */
/*   X,Y,Z   : TABLEAUX DES COORDONNEES DES POINTS DE L ELEMENT */
/*   ijt     : permutation pour oasser de la numerotation par inconnues */
/*             a celle par noeuds */
/*   nno     : nombre de noeuds de l'element */
/*   npo     : nombre de points */
/*   npi     : nombre de points d'integration */
/*   dp      : valeur des derivees des polynomes de base aux points */
/*             d'integration sur l'element de reference */
/*   vp1     : valeur des polynomes de base aux points */
/*             d'integration sur l'element de reference */
/*   poids   : poids de la formule d'integration */

/*  tableaux de travail : */
/*  ------------------- */
/*   delta   : jacobien aux points d'integration */
/*   (x y z)int : coordonnees des points d'integration sur */
/*              l'element courant */

/*  PARAMETRE DE SORTIE  : */
/*  -------------------- */
/*   AE     : MATRICE DE RIGIDITE . SYMETRIQUE */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  PROGRAMMEURS  : Marina Vidrascu 2001 */
/* ................................................................... */


    /* Parameter adjustments */
    coor -= 7;
    --car;
    --ae;

    /* Function Body */



    if (*iopt == 2) {

/*         CAS ISOTROPE */
/*         ------------ */
	er3c2c_(&c__6, &c__6, &coor[7], &coor[13], &coor[19], &c__6, ijt, poids, vp1, dp, dp, iopt, &car[1], &car[2], elas, &ae[1], delta, xint, yint, zint, a2);
    } else {

/*         CAS ORTHOTROPE */
/*         -------------- */
	er3c2c_(&c__6, &c__6, &coor[7], &coor[13], &coor[19], &c__6, ijt, poids, vp1, dp, dp, iopt, &car[1], &car[2], &car[1], &ae[1], delta, xint, yint, zint, a2);
    }
} /* etr3r1d_ */

