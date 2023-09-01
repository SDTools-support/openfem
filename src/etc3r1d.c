/* etc3r1d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__6 = 6;

/* Subroutine */ int etc3r1d_(coor, car, iopt, sigma, u)
doublereal *coor, *car;
int32 *iopt;
doublereal *sigma, *u;
{
    /* Initialized data */

    static int32 ijt[18] = { 1,4,7,10,13,16,2,5,8,11,14,17,3,6,9,12,15,18 };
    static doublereal dp[108]	/* was [3][6][6] */ = { -1.,-1.,-1.,1.,0.,0.,0.,1.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,-1.,-1.,0.,1.,0.,-1.,0.,1.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,-1.,-1.,0.,1.,0.,0.,0.,1.,-1.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,-1.,0.,0.,0.,0.,0.,0.,-1.,-1.,1.,1.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,-1.,0.,0.,0.,-1.,-1.,0.,1.,0.,1.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-1.,-1.,-1.,0.,1.,0.,0.,0.,1.,1. };

    static doublereal elas[9];
    extern /* Subroutine */ int ec3c2c_();
    static doublereal delta[6], a2[108]	/* was [3][6][6] */, sigmae[648]	/* was [6][18][6] */;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  BUT   : CALCUL DE LA CONTRAINTE */
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

/*  PARAMETRE DE SORTIE  : */
/*  -------------------- */
/*   sigma     : */


/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  PROGRAMMEURS  : Marina Vidrascu 2003 */
/* ................................................................... */


    /* Parameter adjustments */
    coor -= 7;
    --car;
    sigma -= 7;
    --u;

    /* Function Body */

/*     si isotrope e=car(1), nu=car(2) */
/*     si orthotrope elas=car */
    if (*iopt == 2) {

/*         CAS ISOTROPE */
/*         ------------ */
	ec3c2c_(&c__6, &c__6, &coor[7], &coor[13], &coor[19], ijt, dp, dp, iopt, &car[1], &car[2], elas, sigmae, &sigma[7], &u[1], delta, a2);
    } else {

/*         CAS ORTHOTROPE */
/*         -------------- */
	ec3c2c_(&c__6, &c__6, &coor[7], &coor[13], &coor[19], ijt, dp, dp, iopt, &car[1], &car[2], &car[1], sigmae, &sigma[7], &u[1], delta, a2);
    }
} /* etc3r1d_ */

