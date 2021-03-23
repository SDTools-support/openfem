/* etr3p1d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__4 = 4;

/* Subroutine */ int etr3p1d_(coor, car, iopt, ae)
doublereal *coor, *car;
int32 *iopt;
doublereal *ae;
{
    /* Initialized data */

    static doublereal dp[48]	/* was [3][4][4] */ = { -1.,-1.,-1.,1.,0.,0.,0.,1.,0.,0.,0.,1.,-1.,-1.,-1.,1.,0.,0.,0.,1.,0.,0.,0.,1.,-1.,-1.,-1.,1.,0.,0.,0.,1.,0.,0.,0.,1.,-1.,-1.,-1.,1.,0.,0.,0.,1.,0.,0.,0.,1. };
    static int32 ijt[12] = { 1,4,7,10,2,5,8,11,3,6,9,12 };
    static doublereal vp1[16]	/* was [4][4] */ = { 1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1. };
    static doublereal poids[4] = { .04166666666666666,.04166666666666666,.04166666666666666,.04166666666666666 };

    static doublereal elas[9], xint[4], yint[4], zint[4];
    extern /* Subroutine */ int er3c2c_();
    static doublereal delta[4], a2[48]	/* was [3][4][4] */;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  BUT   : CALCUL DE LA MATRICE ELEMENTAIRE ELASTIQUE DE RIGIDITE */
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
/*   AE     : MATRICE DE RIGIDITE . SYMETRIQUE */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  PROGRAMMEURS  : Marina Vidrascu 2001 */
/* ................................................................... */


    /* Parameter adjustments */
    coor -= 5;
    --car;
    --ae;

    /* Function Body */

/*     poids = .25/6. */

    if (*iopt == 2) {

/*         CAS ISOTROPE */
/*         ------------ */
	er3c2c_(&c__4, &c__4, &coor[5], &coor[9], &coor[13], &c__4, ijt, poids, vp1, dp, dp, iopt, &car[1], &car[2], elas, &ae[1], delta, xint, yint, zint, a2);
    } else {

/*         CAS ORTHOTROPE */
/*         -------------- */
	er3c2c_(&c__4, &c__4, &coor[5], &coor[9], &coor[13], &c__4, ijt, poids, vp1, dp, dp, iopt, &car[1], &car[2], &car[1], &ae[1], delta, xint, yint, zint, a2);
    }
} /* etr3p1d_ */

