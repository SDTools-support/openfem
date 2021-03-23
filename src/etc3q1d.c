/* etc3q1d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__8 = 8;

/* Subroutine */ int etc3q1d_(coor, car, iopt, sigma, u)
doublereal *coor, *car;
int32 *iopt;
doublereal *sigma, *u;
{
    /* Initialized data */

    static int32 ijt[24] = { 1,4,7,10,13,16,19,22,2,5,8,11,14,17,20,23,3,6,9,12,15,18,21,24 };
    static doublereal dp[192]	/* was [3][8][8] */ = { -.5,-.5,-.5,.5,0.,0.,0.,0.,0.,0.,.5,0.,0.,0.,.5,0.,0.,0.,0.,0.,0.,0.,0.,0.,-.5,0.,0.,.5,-.5,-.5,0.,.5,0.,0.,0.,0.,0.,0.,0.,0.,0.,.5,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-.5,0.,.5,.5,-.5,-.5,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,.5,0.,0.,0.,0.,-.5,0.,0.,0.,0.,.5,0.,0.,-.5,.5,-.5,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,.5,0.,0.,-.5,0.,0.,0.,0.,0.,0.,0.,0.,0.,-.5,-.5,.5,.5,0.,0.,0.,0.,0.,0.,.5,0.,0.,0.,0.,0.,0.,-.5,0.,0.,0.,0.,0.,0.,-.5,0.,0.,.5,-.5,.5,0.,.5,0.,0.,
	    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-.5,0.,0.,0.,0.,0.,0.,0.,-.5,0.,.5,.5,.5,-.5,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-.5,0.,-.5,0.,0.,0.,0.,.5,0.,0.,-.5,.5,.5 };

    static doublereal elas[9];
    extern /* Subroutine */ int ec3c2c_();
    static doublereal delta[8], a2[192]	/* was [3][8][8] */, sigmae[1152]	/* was [6][24][8] */;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  BUT   : CALCUL DES CONTRAINTES */
/*  ---     hexaedre q1 DROIT */
/*          calcul aux noeuds */

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

/*  tableaux de travail : */
/*  ------------------- */
/*   delta   : jacobien aux points d'integration */

/*  PARAMETRE DE SORTIE  : */
/*  -------------------- */
/*   SIGMA     : Tenseur des contraintes, aux noeuds SIGMA(6,nno) */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  PROGRAMMEURS  : Marina Vidrascu 2001 */
/* ................................................................... */


    /* Parameter adjustments */
    coor -= 9;
    --car;
    sigma -= 7;
    --u;

    /* Function Body */

/*     si isotrope e=car(1), nu=car(2) */
/*     si orthotrope elas=car */
    if (*iopt == 2) {

/*         CAS ISOTROPE */
/*         ------------ */
	ec3c2c_(&c__8, &c__8, &coor[9], &coor[17], &coor[25], ijt, dp, dp, iopt, &car[1], &car[2], elas, sigmae, &sigma[7], &u[1], delta, a2);
    } else {

/*         CAS ORTHOTROPE */
/*         -------------- */
	ec3c2c_(&c__8, &c__8, &coor[9], &coor[17], &coor[25], ijt, dp, dp, iopt, &car[1], &car[2], &car[1], sigmae, &sigma[7], &u[1], delta, a2);
    }
} /* etc3q1d_ */

