/* etc3p1d.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__4 = 4;

/* Subroutine */ int etc3p1d_(coor, car, iopt, sigma, u)
doublereal *coor, *car;
int32 *iopt;
doublereal *sigma, *u;
{
    /* Initialized data */

    static int32 ijt[12] = { 1,4,7,10,2,5,8,11,3,6,9,12 };
    static doublereal dp[48]	/* was [3][4][4] */ = { -1.,-1.,-1.,1.,0.,0.,0.,1.,0.,0.,0.,1.,-1.,-1.,-1.,1.,0.,0.,0.,1.,0.,0.,0.,1.,-1.,-1.,-1.,1.,0.,0.,0.,1.,0.,0.,0.,1.,-1.,-1.,-1.,1.,0.,0.,0.,1.,0.,0.,0.,1. };

    static doublereal elas[9];
    extern /* Subroutine */ int ec3c2c_();
    static doublereal delta[4], a2[48]	/* was [3][4][4] */, sigmae[288]	/* was [6][12][4] */;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  BUT   : CALCUL DES CONTRAINTES TETRAEDRE P1 DROIT */
/*  --- */
/*          aux points d'integration, les 4 sommets */

/*  PARAMETRES D ENTREE  : */
/*  ------------------- */
/*   coor  : TABLEAUX DES COORDONNEES DES POINTS DE L ELEMENT */
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
/*            STockee dans ae(nmax,nmax) */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  PROGRAMMEURS  : Marina Vidrascu 2003 */
/* ................................................................... */


    /* Parameter adjustments */
    coor -= 5;
    --car;
    sigma -= 7;
    --u;

    /* Function Body */


/*     si isotrope e=car(1), nu=car(2) */
/*     si orthotrope elas=car */
    if (*iopt == 2) {

/*         CAS ISOTROPE */
/*         ------------ */
	ec3c2c_(&c__4, &c__4, &coor[5], &coor[9], &coor[13], ijt, dp, dp, iopt, &car[1], &car[2], elas, sigmae, &sigma[7], &u[1], delta, a2);
    } else {

/*         CAS ORTHOTROPE */
/*         -------------- */
	ec3c2c_(&c__4, &c__4, &coor[5], &coor[9], &coor[13], ijt, dp, dp, iopt, &car[1], &car[2], &car[1], sigmae, &sigma[7], &u[1], delta, a2);
    }
} /* etc3p1d_ */

