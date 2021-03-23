/* etm3p1d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__4 = 4;

/* Subroutine */ int etm3p1d_(coor, ro, iopt, ae)
doublereal *coor, *ro;
int32 *iopt;
doublereal *ae;
{
    /* Initialized data */

    static doublereal dp[48]	/* was [3][4][4] */ = { -1.,-1.,-1.,1.,0.,0.,0.,1.,0.,0.,0.,1.,-1.,-1.,-1.,1.,0.,0.,0.,1.,0.,0.,0.,1.,-1.,-1.,-1.,1.,0.,0.,0.,1.,0.,0.,0.,1.,-1.,-1.,-1.,1.,0.,0.,0.,1.,0.,0.,0.,1. };
    static int32 ijt[12] = { 1,4,7,10,2,5,8,11,3,6,9,12 };
    static doublereal vp1[16]	/* was [4][4] */ = { 1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1. };
    static doublereal poids[4] = { .04166666666666666,.04166666666666666,.04166666666666666,.04166666666666666 };

    extern /* Subroutine */ int em3c2c_();
    static doublereal delta[4];

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  BUT   : CALCUL DE LA MATRICE ELEMENTAIRE de masse */
/*  ---     TETRAEDRE P1 DROIT */
/*          formule d'integration a 4 points, les sommets */

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
/*   AE     : MATRICE DE Masse. SYMETRIQUE */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  PROGRAMMEURS  : Marina Vidrascu 2001 */
/* ................................................................... */


    /* Parameter adjustments */
    --ae;
    coor -= 5;

    /* Function Body */

/*     poids = .25/6. */

    em3c2c_(&c__4, &c__4, &coor[5], &coor[9], &coor[13], &c__4, ijt, poids, vp1, dp, ro, &ae[1], delta);
} /* etm3p1d_ */

