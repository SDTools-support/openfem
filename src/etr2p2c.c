/* etr2p2c.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int etr2p2c_(coor, car, iopt, ae)
doublereal *coor, *car;
int32 *iopt;
doublereal *ae;
{
    /* Initialized data */

    static int32 ijt[12] = { 1,3,5,7,9,11,2,4,6,8,10,12 };
    static int32 iblt[3] = { 0,0,6 };
    static int32 jblt[3] = { 0,6,6 };
    static doublereal dpx[18]	/* was [3][6] */ = { -1.,1.,-1.,1.,1.,-1.,0.,0.,0.,0.,-2.,2.,0.,2.,2.,0.,-2.,-2. };
    static doublereal dpy[18]	/* was [3][6] */ = { -1.,1.,-1.,0.,0.,0.,-1.,1.,1.,-2.,-2.,0.,2.,2.,0.,2.,-2.,0. };

    /* System generated locals */
    int32 i__1;

    /* Local variables */
    static int32 ifin;
    static doublereal c__, e[6], fdedf[12]	/* was [4][3] */;
    static int32 i__, j, ibloc;
    static doublereal delta[3], dfinv[12]	/* was [4][3] */, c1;
    static int32 i1, i2, i3, i4, j1;
    static doublereal young, unmnu;
    static int32 kk;
    static doublereal x21, y21, x31, y31, x32, y32, x41, y41, x42, y42, x54, y54, x61, y61, x63, y63, x65, y65;
    static int32 kin;
    static doublereal poisson;

/*  .................................................................... */
/* but : matrice de rigidite de l element membrane: TRIA 2P2C */
/* --- */
/* in : coor(noe,ndim) : coor. 3 sommets + 3 milieux aretes */
/*      iopt = 1 isotrope Contraintes  Planes */
/*           = 2 isotrope Deformations Planes */
/*           = sinon anisotrope */
/*      car(6): caracteristiques des materiaux */
/*              if(iopt .eq. 1 .or. iopt.eq. 2) then */
/*                car(1) = young */
/*                car(2) = poisson */
/*              else */
/*                car: E11, E12, E22, E13, E23, E33 avec */

/*                     E11   E12   E13 */
/*                           E22   E23 */
/*                                 E33 */
/*              end if */
/* out: ae(78)        : matrice triangulaire sup */
/* ..................................................................... */

    /* Parameter adjustments */
    --ae;
    --car;
    coor -= 7;

    /* Function Body */

    x21 = coor[8] - coor[7];
    y21 = coor[14] - coor[13];
    x31 = coor[9] - coor[7];
    y31 = coor[15] - coor[13];
    x32 = coor[9] - coor[8];
    y32 = coor[15] - coor[14];
    x41 = coor[10] - coor[7];
    y41 = coor[16] - coor[13];
    x42 = coor[10] - coor[8];
    y42 = coor[16] - coor[14];
    x54 = coor[11] - coor[10];
    y54 = coor[17] - coor[16];
    x61 = coor[12] - coor[7];
    y61 = coor[18] - coor[13];
    x63 = coor[12] - coor[9];
    y63 = coor[18] - coor[15];
    x65 = coor[12] - coor[11];
    y65 = coor[18] - coor[17];

    if (*iopt == 1) {
/*  --    CONTRAINTES PLANES     (ISOTROPE)     ----- */
	young = car[1];
	poisson = car[2];
	c__ = young / (1. - poisson * poisson);
	e[0] = c__;
	e[1] = c__ * poisson;
	e[2] = c__;
	e[3] = 0.;
	e[4] = 0.;
	e[5] = c__ * (1. - poisson) / 2.;
    } else if (*iopt == 2) {
/*  --    DEFORMATIONS PLANES (ISOTROPE)     ----- */
	young = car[1];
	poisson = car[2];
	unmnu = 1. - poisson;
	c__ = young * unmnu;
	c__ /= (poisson + 1.) * (1. - poisson * 2.);
	e[0] = c__;
	e[1] = poisson * c__ / unmnu;
	e[2] = c__;
	e[3] = 0.;
	e[4] = 0.;
	e[5] = c__ * (1. - poisson * 2.) / (unmnu * 2.);
    } else {
/*  --    CAS ANISOTROPE     ----- */
	for (i__ = 1; i__ <= 6; ++i__) {
	    e[i__ - 1] = car[i__];
/* L1: */
	}
    }

/*  ----       CALCUL DE DFINV :  INVERSE DE DF    ---- */

    dfinv[0] = y54 * 2. + y61 + y63;
    dfinv[4] = y54 * 2. - y61 - y63;
    dfinv[8] = y31;

    dfinv[1] = -(x54 * 2. + x61 + x63);
    dfinv[5] = -(x54 * 2. - x61 - x63);
    dfinv[9] = -x31;

    dfinv[2] = -y21;
    dfinv[6] = y65 * 2. + y41 + y42;
    dfinv[10] = y65 * 2. - y41 - y42;

    dfinv[3] = x21;
    dfinv[7] = -(x65 * 2. + x41 + x42);
    dfinv[11] = -(x65 * 2. - x41 - x42);

    for (j = 1; j <= 3; ++j) {
	delta[j - 1] = dfinv[(j << 2) - 1] * dfinv[(j << 2) - 4] - dfinv[(j << 2) - 2] * dfinv[(j << 2) - 3];
/* L2: */
    }

/*  -----      COEFFICIENT DE LA MATRICE ELEMENTAIRE     ----- */

    for (ibloc = 1; ibloc <= 3; ++ibloc) {
	if (*iopt == 1 || *iopt == 2) {
/*         -- isotrope Deformations planes ou Contraintes planes */
	    if (ibloc == 1) {
/*           -- BLOC  1,1 DE    TDF * TD * E * D * DF   (ISOTROPE) */
		for (j = 1; j <= 3; ++j) {
		    fdedf[(j << 2) - 4] = e[0] * dfinv[(j << 2) - 4] * dfinv[(j << 2) - 4] + e[5] * dfinv[(j << 2) - 3] * dfinv[(j << 2) - 3];
		    fdedf[(j << 2) - 3] = e[0] * dfinv[(j << 2) - 2] * dfinv[(j << 2) - 4] + e[5] * dfinv[(j << 2) - 1] * dfinv[(j << 2) - 3];
		    fdedf[(j << 2) - 2] = fdedf[(j << 2) - 3];
		    fdedf[(j << 2) - 1] = e[0] * dfinv[(j << 2) - 2] * dfinv[(j << 2) - 2] + e[5] * dfinv[(j << 2) - 1] * dfinv[(j << 2) - 1];
/* L3: */
		}
	    } else if (ibloc == 2) {
/*           --- BLOC  1,2 */
		for (j = 1; j <= 3; ++j) {
		    fdedf[(j << 2) - 4] = (e[1] + e[5]) * dfinv[(j << 2) - 3] * dfinv[(j << 2) - 4];
		    fdedf[(j << 2) - 3] = e[5] * dfinv[(j << 2) - 1] * dfinv[(j << 2) - 4] + e[1] * dfinv[(j << 2) - 3] * dfinv[(j << 2) - 2];
		    fdedf[(j << 2) - 2] = e[5] * dfinv[(j << 2) - 3] * dfinv[(j << 2) - 2] + e[1] * dfinv[(j << 2) - 4] * dfinv[(j << 2) - 1];
		    fdedf[(j << 2) - 1] = (e[1] + e[5]) * dfinv[(j << 2) - 1] * dfinv[(j << 2) - 2];
/* L4: */
		}
	    } else if (ibloc == 3) {
/*           --- BLOC  2,2 */
		for (j = 1; j <= 3; ++j) {
		    fdedf[(j << 2) - 4] = e[2] * dfinv[(j << 2) - 3] * dfinv[(j << 2) - 3] + e[5] * dfinv[(j << 2) - 4] * dfinv[(j << 2) - 4];
		    fdedf[(j << 2) - 3] = e[2] * dfinv[(j << 2) - 1] * dfinv[(j << 2) - 3] + e[5] * dfinv[(j << 2) - 2] * dfinv[(j << 2) - 4];
		    fdedf[(j << 2) - 2] = fdedf[(j << 2) - 3];
		    fdedf[(j << 2) - 1] = e[2] * dfinv[(j << 2) - 1] * dfinv[(j << 2) - 1] + e[5] * dfinv[(j << 2) - 2] * dfinv[(j << 2) - 2];
/* L5: */
		}
	    }
	} else {
/*         -- anisotrope */
	    if (ibloc == 1) {
/*           ---  BLOC  1,1        CAS   ANISOTROPE */
		i1 = 1;
		i2 = 4;
		i3 = 4;
		i4 = 6;
	    } else if (ibloc == 2) {
/*           --- BLOC  1,2 */
		i1 = 4;
		i2 = 6;
		i3 = 2;
		i4 = 5;
	    } else if (ibloc == 3) {
/*           --- BLOC  2,2 */
		i1 = 6;
		i2 = 5;
		i3 = 5;
		i4 = 3;
	    }

	    for (j = 1; j <= 3; ++j) {
		fdedf[(j << 2) - 4] = e[i1 - 1] * dfinv[(j << 2) - 4] * dfinv[(j << 2) - 4] + e[i4 - 1] * dfinv[(j << 2) - 3] * dfinv[(j << 2) - 3] + (e[i2 - 1] + e[i3 - 1]) * dfinv[(j << 2) - 3] * dfinv[(j << 2) - 4];
		fdedf[(j << 2) - 3] = e[i1 - 1] * dfinv[(j << 2) - 2] * dfinv[(j << 2) - 4] + e[i4 - 1] * dfinv[(j << 2) - 1] * dfinv[(j << 2) - 3] + e[i2 - 1] * dfinv[(j << 2) - 1] * dfinv[(j << 2) - 4] + e[i3 - 1] * dfinv[(j << 2) - 3] * dfinv[(j << 2) - 2];
		fdedf[(j << 2) - 2] = e[i1 - 1] * dfinv[(j << 2) - 2] * dfinv[(j << 2) - 4] + e[i4 - 1] * dfinv[(j << 2) - 1] * dfinv[(j << 2) - 3] + e[i3 - 1] * dfinv[(j << 2) - 1] * dfinv[(j << 2) - 4] + e[i2 - 1] * dfinv[(j << 2) - 3] * dfinv[(j << 2) - 2];
		fdedf[(j << 2) - 1] = e[i4 - 1] * dfinv[(j << 2) - 1] * dfinv[(j << 2) - 1] + e[i1 - 1] * dfinv[(j << 2) - 2] * dfinv[(j << 2) - 2] + (e[i2 - 1] + e[i3 - 1]) * dfinv[(j << 2) - 1] * dfinv[(j << 2) - 2];
/* L6: */
	    }
	}

/*       -- CALCUL DE AE  ( CONTRIBUTION DU BLOC ) */
/*       --   1.  GESTION DES BLOCS ( IBLOC ) */
/*       --   2.  PASSAGE AU RANGEMENT PAR D.L. ( IJT ) */

	for (j = 1; j <= 6; ++j) {
	    ifin = j;
	    if (ibloc == 2) {
		ifin = 6;
	    }
	    i__1 = ifin;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i1 = iblt[ibloc - 1] + i__;
		j1 = jblt[ibloc - 1] + j;
		if (ijt[i1 - 1] <= ijt[j1 - 1]) {
		    kk = ijt[j1 - 1] * (ijt[j1 - 1] - 1) / 2 + ijt[i1 - 1];
		} else {
		    kk = ijt[i1 - 1] * (ijt[i1 - 1] - 1) / 2 + ijt[j1 - 1];
		}
		c__ = 0.;
		for (kin = 1; kin <= 3; ++kin) {
		    c1 = fdedf[(kin << 2) - 4] * dpx[kin + j * 3 - 4] * dpx[kin + i__ * 3 - 4] + fdedf[(kin << 2) - 3] * dpx[kin + j * 3 - 4] * dpy[kin + i__ * 3 - 4] + fdedf[(kin << 2) - 2] * dpx[kin + i__ * 3 - 4] * dpy[kin + j * 3 - 4] + fdedf[(kin << 2) - 1] * dpy[kin + i__ * 3 - 4] * dpy[kin + j * 3 - 4];
		    c__ += c1 / delta[kin - 1];
/* L7: */
		}
		ae[kk] = c__ / 6.;
/* L8: */
	    }
/* L9: */
	}

/* L10: */
    }

} /* etr2p2c_ */

