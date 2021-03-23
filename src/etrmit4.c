/* etrmit4.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__4 = 4;
static int32 c__210 = 210;
static doublereal c_b9 = 0.;
static int32 c__0 = 0;
static int32 c__1 = 1;
static int32 c__2 = 2;
static int32 c__300 = 300;

/* Subroutine */ int etrmit4_(coor, car, iopt, pol, xnorm, ae)
doublereal *coor, *car;
int32 *iopt;
doublereal *pol, *xnorm, *ae;
{
    /* Initialized data */

    static int32 ijt[20] = { 1,6,11,16,2,7,12,17,3,8,13,18,4,9,14,19,5,10,15,20 };
    static int32 irint[5] = { 0,0,0,1,1 };
    static doublereal poids[4] = { 1.,1.,1.,1. };
    static doublereal xint[4] = { -.5773502691896258,.5773502691896258,.5773502691896258,-.5773502691896258 };
    static doublereal yint[4] = { -.5773502691896258,-.5773502691896258,.5773502691896258,.5773502691896258 };
    static doublereal poiz[2] = { 1.,1. };
    static doublereal zint[2] = { -.5773502691896258,.5773502691896258 };

    extern /* Subroutine */ int chan57mod_();
    static doublereal b[100]	/* was [5][20] */, f[36]	/* was [3][3][4] */;
    static int32 i__;
    extern /* Subroutine */ int chan56_(), dcopy_();
    static doublereal v1[12]	/* was [3][4] */, v2[12]	/* was [3][4] */;
    extern /* Subroutine */ int replo1_();
    extern doublereal inmit4_();
    static doublereal fi[36]	/* was [3][3][4] */, bz[100]	/* was [5][20] */, fz[36]	/* was [3][3][4] */, tt[4];
    static int32 ichoix;
    static doublereal ae5[210];
    extern /* Subroutine */ int fonmit_(), mitlin_();
    static int32 ind;
    static doublereal fiz[36]	/* was [3][3][4] */, bzz[100]	/* was [5][20] */, gzz[4];

/* -----*--*---------*---------*---------*---------*---------*---------*-- */

/* BUT : CALCUL DE LA MATRICE DE RIGIDITE DE L ELEMENT MIT4 */
/* --- */
/* -----*--*---------*---------*---------*---------*---------*---------*-- */
/* parametres d entree : */
/* ------------------- */
/* coor   : coor(npo,ndim) coordonnees */
/* car    : caracteristiques des materiaux */
/*          e, nu, t(nno) */
/* iopt   : options */
/*          iopt(1) = notel = 1 cas stadard, */
/*                           10 si pas interpolatio */
/*          iopt(2) = ichoix >0 regulier */
/*                    -i e=x_i */
/*          iopt(3) = 5 si 5ddl par noud */
/*                    6 si 6 ddl par noud */
/* pol    : pol(96) contient vp2,vdpq2,vdtq2 */
/*          1   : vp2(nno,npi+npt) */
/*          17  : vp2(1,npi+1) */
/*          33  : vdpq2(2,nno,npi) */
/*          65  : vdpt2(2,nno,npt) */

/* parametres resultats : */
/* -------------------- */
/* AE     : MATRICE DE RIGIDITE DE L ELEMENT */
/*          SEULE LA PARTIE SUPERIEURE DE HAUT EN BAS ET DE LA GAUCHE */
/*          VERS LA DROITE EST STOCKEE */

/* -----*--*---------*---------*---------*---------*---------*---------*-- */
/*      parameter (xg3 = sqtr(3./5.), xg2 = 1/sqrt(3)) */



    /* Parameter adjustments */
    --ae;
    xnorm -= 4;
    --pol;
    --iopt;
    --car;
    coor -= 5;

    /* Function Body */

/*     composantes de eps reinterpolees */

/*     integration points en r,s */


/*     integration points en z */


/*     calcul des normales (tying+integration) et V1,V2 (tying) */

    ichoix = iopt[2];
    replo1_(&c__4, &xnorm[4], v1, v2, &ichoix);

/*     calcul de F  f1,i=g,i au tying points */

    ind = 0;
    fonmit_(&c__4, &c__4, &coor[5], &coor[9], &coor[13], &xnorm[4], &car[3], &pol[65], &pol[17], f, fz, &ind, gzz, tt);

/*     F au pt integration */

    ind = 1;
    fonmit_(&c__4, &c__4, &coor[5], &coor[9], &coor[13], &xnorm[4], &car[3], &pol[33], &pol[1], fi, fiz, &ind, gzz, tt);

/*     calcul de la matrice d'elasticite */
/*     --------------------------------- */
/*     ds car il y a ds l ordre E,enu,t(nno) */
/*     t(nno) -- > car(3) */

/*     si notel = 10 pas de reinterpolation */

    if (iopt[1] == 10) {
	for (i__ = 1; i__ <= 5; ++i__) {
	    irint[i__ - 1] = 0;
/* L1: */
	}
    }
    if (iopt[3] == 5) {

/*        --- > 5 ddl/noeud */
	dcopy_(&c__210, &c_b9, &c__0, &ae[1], &c__1);
	mitlin_(&c__4, &c__2, &c__2, &c__0, &c__4, irint, &c__4, poids, &c__2, zint, poiz, xint, yint, &pol[65], &pol[17], f, fz, inmit4_, &car[1], &car[2], &car[3], &coor[5], &coor[9], &coor[13], &pol[33], &pol[1], b, bz, bzz, &xnorm[4], v1, v2, fi, fiz, gzz, tt, ijt, &ae[1], &c__210);
    } else if (iopt[3] == 7) {
/*        --- > 6 ddl/noeud Cas Coventor Plaques minces */

	dcopy_(&c__210, &c_b9, &c__0, ae5, &c__1);
	mitlin_(&c__4, &c__2, &c__2, &c__0, &c__4, irint, &c__4, poids, &c__2, zint, poiz, xint, yint, &pol[65], &pol[17], f, fz, inmit4_, &car[1], &car[2], &car[3], &coor[5], &coor[9], &coor[13], &pol[33], &pol[1], b, bz, bzz, &xnorm[4], v1, v2, fi, fiz, gzz, tt, ijt, ae5, &c__210);
	chan57mod_(&c__4, ae5, &c__210, &ae[1], &c__300, &car[3], &car[1], &coor[5]);
    } else {

/*        --- > 6 ddl/noeud */
	dcopy_(&c__210, &c_b9, &c__0, ae5, &c__1);
	mitlin_(&c__4, &c__2, &c__2, &c__0, &c__4, irint, &c__4, poids, &c__2, zint, poiz, xint, yint, &pol[65], &pol[17], f, fz, inmit4_, &car[1], &car[2], &car[3], &coor[5], &coor[9], &coor[13], &pol[33], &pol[1], b, bz, bzz, &xnorm[4], v1, v2, fi, fiz, gzz, tt, ijt, ae5, &c__210);
	chan56_(&c__4, ae5, &c__210, &ae[1], &c__300, &car[3], v1, v2, &xnorm[4], &car[1]);
    }
} /* etrmit4_ */

