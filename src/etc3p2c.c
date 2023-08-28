/* etc3p2c.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__10 = 10;

/* Subroutine */ int etc3p2c_(coor, car, iopt, sigma, u)
doublereal *coor, *car;
int32 *iopt;
doublereal *sigma, *u;
{
    /* Initialized data */

    static int32 ijt[30] = { 1,4,7,10,13,16,19,22,25,28,2,5,8,11,14,17,20,23,26,29,3,6,9,12,15,18,21,24,27,30 };
    static doublereal dp[300]	/* was [3][10][10] */ = { -3.,-3.,-3.,-1.,0.,0.,0.,-1.,0.,0.,0.,-1.,4.,0.,0.,0.,0.,0.,0.,4.,0.,0.,0.,4.,0.,0.,0.,0.,0.,0.,1.,1.,1.,3.,0.,0.,0.,-1.,0.,0.,0.,-1.,-4.,-4.,-4.,0.,4.,0.,0.,0.,0.,0.,0.,0.,0.,0.,4.,0.,0.,0.,1.,1.,1.,-1.,0.,0.,0.,3.,0.,0.,0.,-1.,0.,0.,0.,4.,0.,0.,-4.,-4.,-4.,0.,0.,0.,0.,0.,0.,0.,0.,4.,1.,1.,1.,-1.,0.,0.,0.,-1.,0.,0.,0.,3.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-4.,-4.,-4.,4.,0.,0.,0.,4.,0.,-1.,-1.,-1.,1.,0.,0.,0.,-1.,0.,0.,0.,-1.,0.,-2.,-2.,0.,2.,0.,
	    0.,2.,0.,0.,0.,2.,0.,0.,2.,0.,0.,0.,1.,1.,1.,1.,0.,0.,0.,1.,0.,0.,0.,-1.,-2.,-2.,-2.,2.,2.,0.,-2.,-2.,-2.,0.,0.,0.,0.,0.,2.,0.,0.,2.,-1.,-1.,-1.,-1.,0.,0.,0.,1.,0.,0.,0.,-1.,2.,0.,0.,2.,0.,0.,-2.,0.,-2.,0.,0.,2.,0.,0.,0.,0.,0.,2.,-1.,-1.,-1.,-1.,0.,0.,0.,-1.,0.,0.,0.,1.,2.,0.,0.,0.,0.,0.,0.,2.,0.,-2.,-2.,0.,2.,0.,0.,0.,2.,0.,1.,1.,1.,1.,0.,0.,0.,-1.,0.,0.,0.,1.,-2.,-2.,-2.,0.,2.,0.,0.,0.,0.,-2.,-2.,-2.,2.,0.,2.,0.,2.,0.,1.,1.,1.,-1.,0.,0.,0.,1.,0.,0.,0.,1.,0.,0.,0.,2.,0.,0.,-2.,-2.,
	    -2.,-2.,-2.,-2.,2.,0.,0.,0.,2.,2. };

    static doublereal elas[9];
    extern /* Subroutine */ int ec3c2c_();
    static doublereal delta[10], a2[300]	/* was [3][10][10] */, sigmae[1800]	/* was [6][30][10] */;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  BUT   : CALCUL DES CONTRAINTES */
/*  ---     TETRAEDRE P2 ISOPARAMETRIQUE */
/*          formule d'integration a */

/*  PARAMETRES D ENTREE  : */
/*  ------------------- */
/*   coor   : TABLEAUX DES COORDONNEES DES POINTS DE L ELEMENT */
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
    coor -= 11;
    --car;
    sigma -= 7;
    --u;

    /* Function Body */


/*     si isotrope e=car(1), nu=car(2) */
/*     si orthotrope elas=car */
    if (*iopt == 2) {

/*         CAS ISOTROPE */
/*         ------------ */
	ec3c2c_(&c__10, &c__10, &coor[11], &coor[21], &coor[31], ijt, dp, dp, iopt, &car[1], &car[2], elas, sigmae, &sigma[7], &u[1], delta, a2);
    } else {

/*         CAS ORTHOTROPE */
/*         -------------- */
	ec3c2c_(&c__10, &c__10, &coor[11], &coor[21], &coor[31], ijt, dp, dp, iopt, &car[1], &car[2], &car[1], sigmae, &sigma[7], &u[1], delta, a2);
    }
} /* etc3p2c_ */

