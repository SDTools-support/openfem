/* melmit.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int melmit_(grr, grrz, grrzz, grs, grsz, grszz, gss, gssz, gsszz, gzz, zint, npiz, npi, enu, e, tt, l, lz, deltat, elas, girr, giss, girs)
doublereal *grr, *grrz, *grrzz, *grs, *grsz, *grszz, *gss, *gssz, *gsszz, *gzz, *zint;
int32 *npiz, *npi;
doublereal *enu, *e, *tt;
int32 *l, *lz;
doublereal *deltat, *elas, *girr, *giss, *girs;
{
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal delta;

/*    ............................................................... */
/*    calcul de la matrice d'elasticite */
/*    .............................................................. */

    /* Parameter adjustments */
    --zint;
    --tt;
    --gzz;
    --elas;

    /* Function Body */
/* Computing 2nd power */
    d__1 = zint[*lz];
/* Computing 2nd power */
    d__2 = zint[*lz];
/* Computing 2nd power */
    d__4 = zint[*lz];
/* Computing 2nd power */
    d__3 = *grs + *grsz * zint[*lz] + *grszz * (d__4 * d__4);
    delta = (*grr + *grrz * zint[*lz] + *grrzz * (d__1 * d__1)) * (*gss + *gssz * zint[*lz] + *gsszz * (d__2 * d__2)) - d__3 * d__3;
/* Computing 2nd power */
    d__1 = zint[*lz];
    *girr = (*gss + *gssz * zint[*lz] + *gsszz * (d__1 * d__1)) / delta;
/* Computing 2nd power */
    d__1 = zint[*lz];
    *giss = (*grr + *grrz * zint[*lz] + *grrzz * (d__1 * d__1)) / delta;
/* Computing 2nd power */
    d__1 = zint[*lz];
    *girs = -(*grs + *grsz * zint[*lz] + *grszz * (d__1 * d__1)) / delta;
    *deltat = sqrt(delta * gzz[*l]);

/*           matrice d'elasticite 2 blocs diagonaux */
/*           -------------------------------------- */
/*           c11 c12 c13                el(1)  el(2)  el(4) */
/*           c12 c22 c23                       el(3)  el(5) */
/*           c13 c23 c33          =                   el(6) */
/*                      c44 c45                            el(7) el(8) */
/*                      c45 c55                                  el(9) */
    elas[1] = *e / ((*enu + 1) * 2) * (*girr * 2 * *girr + *enu * 2 / (1 - *enu) * *girr * *girr);
    elas[2] = *e / ((*enu + 1) * 2) * (*girs * 2 * *girs + *enu * 2 / (1 - *enu) * *girr * *giss);
    elas[4] = *e / ((*enu + 1) * 2) * (*girr * *girs + *girs * *girr + *enu * 2 / (1 - *enu) * *girr * *girs);
    elas[3] = *e / ((*enu + 1) * 2) * (*giss * *giss + *giss * *giss + *enu * 2 / (1 - *enu) * *giss * *giss);
    elas[5] = *e / ((*enu + 1) * 2) * (*girs * *giss + *giss * *girs + *enu * 2 / (1 - *enu) * *giss * *girs);
    elas[6] = *e / ((*enu + 1) * 2) * (*girr * *giss + *girs * *girs + *enu * 2 / (1 - *enu) * *girs * *girs);
/* Computing 2nd power */
    d__1 = tt[*l];
    elas[7] = 4 / (d__1 * d__1) * *e / ((*enu + 1) * 2) * *girr;
/* Computing 2nd power */
    d__1 = tt[*l];
    elas[8] = 4 / (d__1 * d__1) * *e / ((*enu + 1) * 2) * *girs;
/* Computing 2nd power */
    d__1 = tt[*l];
    elas[9] = 4 / (d__1 * d__1) * *e / ((*enu + 1) * 2) * *giss;
} /* melmit_ */

