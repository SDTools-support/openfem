/* mainf.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__1 = 1;
static int32 c__50 = 50;

/* Subroutine */ int mexfunction_(nlhs, plhs, nrhs, prhs)
int32 *nlhs, *plhs, *nrhs, *prhs;
{
    /* Builtin functions */
    int32 s_wsfi(), do_fio(), e_wsfi();

    /* Local variables */
    static doublereal x, y;
    extern /* Subroutine */ int mexprintf_();
    static char string[100];

    /* Fortran I/O blocks */
    static icilist io___4 = { 0, string, 0, "(\"i want this \",F,\" and \",F,\"and furthermore \",I)", 100, 1 };


/* ----------------------------------------------------------------------- */
/*     (integer) Replace int32 by integer*8 on the DEC Alpha and the */
/*     SGI 64-bit platforms */

/* ----------------------------------------------------------------------- */


/*      WRITE(6,*) 'hello' */
    /* Parameter adjustments */
    --prhs;
    --plhs;

    /* Function Body */
    mexprintf_("hello", 5L);
/*     PRINT *, 'hello' */

    x = (float)5.;
    y = (float)10.;
    s_wsfi(&io___4);
    do_fio(&c__1, (char *)&x, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&y, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&c__50, (ftnlen)sizeof(integer));
    e_wsfi();
    mexprintf_(string, 100L);
    return 0;
} /* mexfunction_ */

