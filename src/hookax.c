/* hookax.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int hookax_(iopt, car, elas)
int32 *iopt;
doublereal *car, *elas;
{
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7;

    /* Local variables */
    static doublereal e, e1, e2, e3, nu, nu1, nu2;

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* BUT: CALCUL DE LA MATRICE D ELASTICITE ( 4*4 ) AXISYMETRIQUE . */

/* in : car: caracteristiques des materiaux */
/*           if(iopt .lt. 3) then */
/*             car(1) = young */
/*             car(2) = poisson */
/*           else */
/*             car(1) = E_1  (Young radial     -> E_r     ) */
/*             car(2) = nu_1 (poisson          -> Nu_r    ) */
/*             car(3) = E_2  (Young axial      -> E_z     ) */
/*             car(4) = nu_2 (poisson          -> Nu_z    ) */
/*             car(5) = E_3  (Young Tangentiel -> E_theta ) */
/*           end if */
/* out: elas(10)                        -1 */
/*      Deformation             [ Elas ]                        Stress */
/*         ||                      ||                             || */
/*      | s_zz |   |   1/E_2    -Nu_2/E_2  -Nu_2/E_2   0   |   | s_zz | */
/*      | s_rr | - | -Nu_1/E_1    1/E_1    -Nu_1/E_1   0   |   | s_rr | */
/*      | s_tt | - | -Nu_3/E3   -Nu_3/E_3    1/E_3     0   | * | s_tt | */
/*      | s_rz |   |     0          0          0     1/G_2 |   | s_rz | */

/* isotrope  : E_1 = E_2 = E_3     ;    Nu_1 = Nu_2 = Nu_3 */

/* orthotrope: -Nu_2/E_2 = -Nu_1/E_1 ; -Nu_3/E3 = -Nu_2/E_2 ; */
/*             -Nu_3/E_3 = -Nu_1/E_1 ; G_2 = E_2 / (2*(1+Nu_2)) */
/* .................................................................... */
/* ............................................................... */
/*  programmeur : modulef */
/* ............................................................... */

    /* Parameter adjustments */
    --elas;
    --car;

    /* Function Body */
    if (*iopt < 3) {
	e = car[1];
	nu = car[2];
/* Computing 2nd power */
	d__1 = nu;
	elas[1] = e * (nu - 1) / (d__1 * d__1 * 2 + nu - 1);
/* Computing 2nd power */
	d__1 = nu;
	elas[2] = -e * nu / (d__1 * d__1 * 2 + nu - 1);
/* Computing 2nd power */
	d__1 = nu;
	elas[3] = e * (nu - 1) / (d__1 * d__1 * 2 + nu - 1);
/* Computing 2nd power */
	d__1 = nu;
	elas[4] = -e * nu / (d__1 * d__1 * 2 + nu - 1);
/* Computing 2nd power */
	d__1 = nu;
	elas[5] = -e * nu / (d__1 * d__1 * 2 + nu - 1);
/* Computing 2nd power */
	d__1 = nu;
	elas[6] = e * (nu - 1) / (d__1 * d__1 * 2 + nu - 1);
	elas[7] = (float)0.;
	elas[8] = (float)0.;
	elas[9] = (float)0.;
	elas[10] = e / ((nu + 1) * 2);
    } else {
	e1 = car[1];
	nu1 = car[2];
	e2 = car[3];
	nu2 = car[4];
	e3 = car[5];
/* Computing 2nd power */
	d__1 = nu1;
/* Computing 2nd power */
	d__2 = nu1;
/* Computing 2nd power */
	d__3 = nu2;
/* Computing 2nd power */
	d__4 = e1;
/* Computing 2nd power */
	d__5 = nu2;
/* Computing 2nd power */
	d__6 = nu2;
/* Computing 2nd power */
	d__7 = e2;
	elas[1] = (e1 - d__1 * d__1 * e3) / (e2 * e1 - e2 * (d__2 * d__2) * e3 - d__3 * d__3 * (d__4 * d__4) - d__5 * d__5 * 2 * e1 * nu1 * e3 - d__6 * d__6 * e1 * e3) * (d__7 * d__7);
/* Computing 2nd power */
	d__1 = nu1;
/* Computing 2nd power */
	d__2 = nu2;
/* Computing 2nd power */
	d__3 = e1;
/* Computing 2nd power */
	d__4 = nu2;
/* Computing 2nd power */
	d__5 = nu2;
	elas[2] = nu2 * (e1 + nu1 * e3) * e2 * e1 / (e2 * e1 - e2 * (d__1 * d__1) * e3 - d__2 * d__2 * (d__3 * d__3) - d__4 * d__4 * 2 * e1 * nu1 * e3 - d__5 * d__5 * e1 * e3);
/* Computing 2nd power */
	d__1 = nu2;
/* Computing 2nd power */
	d__2 = nu1;
/* Computing 2nd power */
	d__3 = nu2;
/* Computing 2nd power */
	d__4 = e1;
/* Computing 2nd power */
	d__5 = nu2;
/* Computing 2nd power */
	d__6 = nu2;
/* Computing 2nd power */
	d__7 = e1;
	elas[3] = (e2 - d__1 * d__1 * e3) / (e2 * e1 - e2 * (d__2 * d__2) * e3 - d__3 * d__3 * (d__4 * d__4) - d__5 * d__5 * 2 * e1 * nu1 * e3 - d__6 * d__6 * e1 * e3) * (d__7 * d__7);
/* Computing 2nd power */
	d__1 = nu1;
/* Computing 2nd power */
	d__2 = nu2;
/* Computing 2nd power */
	d__3 = e1;
/* Computing 2nd power */
	d__4 = nu2;
/* Computing 2nd power */
	d__5 = nu2;
	elas[4] = nu2 * (nu1 + 1) * e2 * e1 / (e2 * e1 - e2 * (d__1 * d__1) * e3 - d__2 * d__2 * (d__3 * d__3) - d__4 * d__4 * 2 * e1 * nu1 * e3 - d__5 * d__5 * e1 * e3) * e3;
/* Computing 2nd power */
	d__1 = nu2;
/* Computing 2nd power */
	d__2 = nu1;
/* Computing 2nd power */
	d__3 = nu2;
/* Computing 2nd power */
	d__4 = e1;
/* Computing 2nd power */
	d__5 = nu2;
/* Computing 2nd power */
	d__6 = nu2;
	elas[5] = (nu1 * e2 + d__1 * d__1 * e1) * e1 / (e2 * e1 - e2 * (d__2 * d__2) * e3 - d__3 * d__3 * (d__4 * d__4) - d__5 * d__5 * 2 * e1 * nu1 * e3 - d__6 * d__6 * e1 * e3) * e3;
/* Computing 2nd power */
	d__1 = nu2;
/* Computing 2nd power */
	d__2 = nu1;
/* Computing 2nd power */
	d__3 = nu2;
/* Computing 2nd power */
	d__4 = e1;
/* Computing 2nd power */
	d__5 = nu2;
/* Computing 2nd power */
	d__6 = nu2;
	elas[6] = (e2 - d__1 * d__1 * e1) * e1 / (e2 * e1 - e2 * (d__2 * d__2) * e3 - d__3 * d__3 * (d__4 * d__4) - d__5 * d__5 * 2 * e1 * nu1 * e3 - d__6 * d__6 * e1 * e3) * e3;
	elas[7] = (float)0.;
	elas[8] = (float)0.;
	elas[9] = (float)0.;
	elas[10] = e2 / ((nu2 + 1) * 2);
    }

    return 0;
} /* hookax_ */

