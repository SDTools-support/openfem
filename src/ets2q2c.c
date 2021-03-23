/* ets2q2c.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int ets2q2c_(coor, fomega, fgamma, pressi, norefs, be)
doublereal *coor, *fomega, *fgamma, *pressi;
int32 *norefs;
doublereal *be;
{
    /* Initialized data */

    static doublereal poids[9] = { .077160493827160503,.12345679012345678,.077160493827160503,.12345679012345678,.19753086419753085,.12345679012345678,.077160493827160503,.12345679012345678,.077160493827160503 };
    static doublereal p25[72]	/* was [8][9] */ = { .43237900077244512,-.1,-.032379000772445008,-.099999999999999991,.35491933384829665,.045080666151703314,.045080666151703314,.35491933384829665,-.099999999999999977,-.099999999999999977,-.099999999999999977,-.1,.8872983346207417,.19999999999999998,.11270166537925829,.19999999999999998,-.10000000000000003,.43237900077244512,-.099999999999999977,-.032379000772445015,.35491933384829665,.35491933384829665,.045080666151703308,.045080666151703314,
	    -.099999999999999991,-.099999999999999991,-.099999999999999991,-.099999999999999991,.19999999999999998,.11270166537925829,.19999999999999998,.8872983346207417,-.25,-.25,-.25,-.25,.5,.5,.5,.5,-.10000000000000019,-.099999999999999977,-.099999999999999977,-.099999999999999977,.19999999999999995,.8872983346207417,.19999999999999995,.11270166537925829,-.099999999999999977,-.032379000772444987,-.1,.43237900077244512,.045080666151703308,.045080666151703308,.35491933384829665,
	    .3549193338482966,-.10000000000000019,-.099999999999999977,-.099999999999999977,-.099999999999999977,.11270166537925829,.19999999999999995,.8872983346207417,.19999999999999995,-.032379000772445598,-.099999999999999866,.43237900077244484,-.10000000000000008,.045080666151703141,.35491933384829676,.35491933384829676,.045080666151703141 };
    static doublereal dp25[144]	/* was [2][8][9] */ = { -2.061895003862225,-2.061895003862225,-.68729833462074174,-.087298334620741685,-.26189500386222502,-.26189500386222502,-.087298334620741685,-.68729833462074174,2.7491933384829669,-.39999999999999996,.39999999999999996,.34919333848296674,.34919333848296674,.39999999999999996,-.39999999999999996,2.7491933384829669,-.68729833462074152,-.7745966692414834,.68729833462074174,-.7745966692414834,-.087298334620741657,-.7745966692414834,
	    .087298334620741685,-.7745966692414834,0.,-1.,.39999999999999996,1.5491933384829668,0.,1.,-.39999999999999996,1.5491933384829668,.68729833462074196,-.087298334620742101,2.0618950038622254,-2.061895003862225,.087298334620741713,-.68729833462074174,.26189500386222508,-.26189500386222497,-2.7491933384829669,-.39999999999999991,.39999999999999996,2.7491933384829669,-.34919333848296674,.39999999999999991,-.39999999999999996,.34919333848296674,-.7745966692414834,-.68729833462074174,
	    -.7745966692414834,.087298334620741685,-.7745966692414834,-.087298334620741685,-.7745966692414834,.68729833462074174,1.5491933384829668,-.39999999999999996,1.,0.,1.5491933384829668,.39999999999999996,-1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-1.,1.,0.,0.,1.,-1.,0.,.7745966692414834,.087298334620741213,.7745966692414834,-.68729833462074174,.7745966692414834,.68729833462074174,.7745966692414834,-.087298334620741657,-1.5491933384829668,-.39999999999999991,1.,0.,-1.5491933384829668,
	    .39999999999999991,-1.,0.,-.087298334620742157,.68729833462074174,-.26189500386222502,.26189500386222502,-.68729833462074174,.087298334620741685,-2.0618950038622254,2.061895003862225,.34919333848296674,-.39999999999999996,.39999999999999991,-.34919333848296674,2.7491933384829669,.39999999999999996,-.39999999999999991,-2.7491933384829669,.087298334620741213,.7745966692414834,-.087298334620741657,.7745966692414834,.68729833462074174,.7745966692414834,-.68729833462074196,
	    .7745966692414834,0.,-1.,.39999999999999991,-1.5491933384829668,0.,1.,-.39999999999999991,-1.5491933384829668,.26189500386222475,.26189500386222452,.087298334620741879,.68729833462074174,2.061895003862225,2.061895003862225,.68729833462074152,.087298334620741657,-.34919333848296663,-.39999999999999991,.39999999999999991,-2.7491933384829669,-2.7491933384829669,.39999999999999991,-.39999999999999991,-.34919333848296663 };

    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal aret[24]	/* was [2][12] */, fface[18]	/* was [2][9] */;
    static int32 i__, j, k, m;
    static doublereal delta[9], desur[3], f11[9], f12[9], f21[9], f22[9], x21, y21, x32, x51, y51, x52, y52, x62, y62, x63, y63, y32, x73, y73, x74, sx, sy, y74, x43, y43, arelon, x81, y81, x84, y84, x41, y41, xmi, ymi, xjm, yjm, xnu[3], ynu[3];

/*  .................................................................... */
/* but : second membre de l element membrane:  QUAD 2Q2C */
/* --- */
/* in : coor(noe,ndim) : coordonees 8 noeuds */
/*      fomega(2,8)  : fx, fy volumiques aux 8 noeuds */
/*      fgamma(ndim,3*nbarete): fx, fy aux 3 noeuds de chaque arete. */
/*      pressi(3*nbarete)     : pression aux 3 noeuds de chaque arete */
/*      norefs(nbarete,2): norefs(i,1) = 0 si fgamma   = 0 sur arete_i */
/*                         norefs(i,2) = 0 si pression = 0 sur arete_i */

/* out: BE(16) */

/* programmeur : modulef */
/* ............................................................... */
/* 2Q25 -- XYNPI: coordonnees pt. int. numeriques (element reference) */
    /* Parameter adjustments */
    --be;
    norefs -= 5;
    --pressi;
    fgamma -= 3;
    fomega -= 3;
    coor -= 9;

    /* Function Body */
/*     -- POIDS: poids du schema d'integration numerique. */
/*     -- Valeurs des Polynomes de base aux pt. int. numerique. */
/*     -- Valeurs Derive'es des Poly. de base aux pt int. numerique. */

/* --  INTEGRATION : TERMES SURFACIQUES  ---- */
/*     -- Valeurs aux pt d'int. num. a partir valeurs aux noeuds */
/*        Efforts volumiques  fomega(ndim,noe) -> fface(ndim,npi) */
    for (i__ = 1; i__ <= 9; ++i__) {
	fface[(i__ << 1) - 2] = 0.;
	fface[(i__ << 1) - 1] = 0.;
	for (j = 1; j <= 8; ++j) {
	    fface[(i__ << 1) - 2] += p25[j + (i__ << 3) - 9] * fomega[(j << 1) + 1];
	    fface[(i__ << 1) - 1] += p25[j + (i__ << 3) - 9] * fomega[(j << 1) + 2];
/* L1: */
	}
/* L2: */
    }
/*     -----   CALCUL DE DELTA AUX NPI NOEUDS  ----- */
    for (k = 1; k <= 9; ++k) {
	f11[k - 1] = 0.;
	f12[k - 1] = 0.;
	f21[k - 1] = 0.;
	f22[k - 1] = 0.;
/* L3: */
    }
    for (k = 1; k <= 9; ++k) {
	for (i__ = 1; i__ <= 8; ++i__) {
	    f11[k - 1] += dp25[(i__ + (k << 3) << 1) - 18] * coor[i__ + 8];
	    f12[k - 1] += dp25[(i__ + (k << 3) << 1) - 17] * coor[i__ + 8];
	    f21[k - 1] += dp25[(i__ + (k << 3) << 1) - 18] * coor[i__ + 16];
	    f22[k - 1] += dp25[(i__ + (k << 3) << 1) - 17] * coor[i__ + 16];
/* L4: */
	}
	delta[k - 1] = f11[k - 1] * f22[k - 1] - f12[k - 1] * f21[k - 1];
/* L5: */
    }
    for (i__ = 1; i__ <= 16; ++i__) {
	be[i__] = 0.;
/* L6: */
    }

    for (j = 1; j <= 8; ++j) {
	sx = 0.;
	sy = 0.;
	for (k = 1; k <= 9; ++k) {
	    sx += poids[k - 1] * delta[k - 1] * p25[j + (k << 3) - 9] * fface[(k << 1) - 2];
	    sy += poids[k - 1] * delta[k - 1] * p25[j + (k << 3) - 9] * fface[(k << 1) - 1];
/* L7: */
	}
	be[(j << 1) - 1] = sx;
	be[j * 2] = sy;
/* L8: */
    }

/*  -- INTEGRATION : TERMES DE BORD  ---- */

/*     -- arete 1 */
    if (norefs[9] != 0) {
	k = 1;
	j = 2;
	m = k + 4;
	xmi = coor[m + 8] - coor[k + 8];
	ymi = coor[m + 16] - coor[k + 16];
/* Computing 2nd power */
	d__1 = xmi;
/* Computing 2nd power */
	d__2 = ymi;
	arelon = sqrt(d__1 * d__1 + d__2 * d__2);
	xnu[0] = ymi / arelon;
	ynu[0] = -xmi / arelon;
	xjm = coor[j + 8] - coor[m + 8];
	yjm = coor[j + 16] - coor[m + 16];
/* Computing 2nd power */
	d__1 = xjm;
/* Computing 2nd power */
	d__2 = yjm;
	arelon = sqrt(d__1 * d__1 + d__2 * d__2);
	xnu[1] = yjm / arelon;
	ynu[1] = -xjm / arelon;
	xnu[2] = (xnu[0] + xnu[1]) * .5;
	ynu[2] = (ynu[0] + ynu[1]) * .5;
	aret[0] = fgamma[3] - pressi[1] * xnu[0];
	aret[1] = fgamma[4] - pressi[1] * ynu[0];
	aret[2] = fgamma[5] - pressi[2] * xnu[1];
	aret[3] = fgamma[6] - pressi[2] * ynu[1];
	aret[4] = fgamma[7] - pressi[3] * xnu[2];
	aret[5] = fgamma[8] - pressi[3] * ynu[2];
    } else {
	aret[0] = fgamma[3];
	aret[1] = fgamma[4];
	aret[2] = fgamma[5];
	aret[3] = fgamma[6];
	aret[4] = fgamma[7];
	aret[5] = fgamma[8];
    }
    x51 = coor[13] - coor[9];
    y51 = coor[21] - coor[17];
    x52 = coor[13] - coor[10];
    y52 = coor[21] - coor[18];
    x21 = coor[10] - coor[9];
    y21 = coor[18] - coor[17];
    desur[0] = sqrt((x51 * 3. + x52) * (x51 * 3. + x52) + (y51 * 3. + y52) * (y51 * 3. + y52)) / 6.;
    desur[1] = sqrt((x52 * 3. + x51) * (x52 * 3. + x51) + (y52 * 3. + y51) * (y52 * 3. + y51)) / 6.;
    desur[2] = sqrt(x21 * x21 + y21 * y21) * 4. / 6.;
    be[1] += desur[0] * aret[0];
    be[2] += desur[0] * aret[1];
    be[3] += desur[1] * aret[2];
    be[4] += desur[1] * aret[3];
    be[9] += desur[2] * aret[4];
    be[10] += desur[2] * aret[5];

/*     -- arete 2 */
    if (norefs[10] != 0) {
	k = 2;
	j = 3;
	m = k + 4;
	xmi = coor[m + 8] - coor[k + 8];
	ymi = coor[m + 16] - coor[k + 16];
/* Computing 2nd power */
	d__1 = xmi;
/* Computing 2nd power */
	d__2 = ymi;
	arelon = sqrt(d__1 * d__1 + d__2 * d__2);
	xnu[0] = ymi / arelon;
	ynu[0] = -xmi / arelon;
	xjm = coor[j + 8] - coor[m + 8];
	yjm = coor[j + 16] - coor[m + 16];
/* Computing 2nd power */
	d__1 = xjm;
/* Computing 2nd power */
	d__2 = yjm;
	arelon = sqrt(d__1 * d__1 + d__2 * d__2);
	xnu[1] = yjm / arelon;
	ynu[1] = -xjm / arelon;
	xnu[2] = (xnu[0] + xnu[1]) * .5;
	ynu[2] = (ynu[0] + ynu[1]) * .5;
	aret[6] = fgamma[9] - pressi[4] * xnu[0];
	aret[7] = fgamma[10] - pressi[4] * ynu[0];
	aret[8] = fgamma[11] - pressi[5] * xnu[1];
	aret[9] = fgamma[12] - pressi[5] * ynu[1];
	aret[10] = fgamma[13] - pressi[6] * xnu[2];
	aret[11] = fgamma[14] - pressi[6] * ynu[2];
    } else {
	aret[6] = fgamma[9];
	aret[7] = fgamma[10];
	aret[8] = fgamma[11];
	aret[9] = fgamma[12];
	aret[10] = fgamma[13];
	aret[11] = fgamma[14];
    }
    x62 = coor[14] - coor[10];
    y62 = coor[22] - coor[18];
    x63 = coor[14] - coor[11];
    y63 = coor[22] - coor[19];
    x32 = coor[11] - coor[10];
    y32 = coor[19] - coor[18];
    desur[0] = sqrt((x62 * 3. + x63) * (x62 * 3. + x63) + (y62 * 3. + y63) * (y62 * 3. + y63)) / 6.;
    desur[1] = sqrt((x63 * 3. + x62) * (x63 * 3. + x62) + (y63 * 3. + y62) * (y63 * 3. + y62)) / 6.;
    desur[2] = sqrt(x32 * x32 + y32 * y32) * 4. / 6.;
    be[3] += desur[0] * aret[6];
    be[4] += desur[0] * aret[7];
    be[5] += desur[1] * aret[8];
    be[6] += desur[1] * aret[9];
    be[11] += desur[2] * aret[10];
    be[12] += desur[2] * aret[11];

/*     -- arete 3 */
    if (norefs[11] != 0) {
	k = 3;
	j = 4;
	m = k + 4;
	xmi = coor[m + 8] - coor[k + 8];
	ymi = coor[m + 16] - coor[k + 16];
/* Computing 2nd power */
	d__1 = xmi;
/* Computing 2nd power */
	d__2 = ymi;
	arelon = sqrt(d__1 * d__1 + d__2 * d__2);
	xnu[0] = ymi / arelon;
	ynu[0] = -xmi / arelon;
	xjm = coor[j + 8] - coor[m + 8];
	yjm = coor[j + 16] - coor[m + 16];
/* Computing 2nd power */
	d__1 = xjm;
/* Computing 2nd power */
	d__2 = yjm;
	arelon = sqrt(d__1 * d__1 + d__2 * d__2);
	xnu[1] = yjm / arelon;
	ynu[1] = -xjm / arelon;
	xnu[2] = (xnu[0] + xnu[1]) * .5;
	ynu[2] = (ynu[0] + ynu[1]) * .5;
	aret[12] = fgamma[15] - pressi[7] * xnu[0];
	aret[13] = fgamma[16] - pressi[7] * ynu[0];
	aret[14] = fgamma[17] - pressi[8] * xnu[1];
	aret[15] = fgamma[18] - pressi[8] * ynu[1];
	aret[16] = fgamma[19] - pressi[9] * xnu[2];
	aret[17] = fgamma[20] - pressi[9] * ynu[2];
    } else {
	aret[12] = fgamma[15];
	aret[13] = fgamma[16];
	aret[14] = fgamma[17];
	aret[15] = fgamma[18];
	aret[16] = fgamma[19];
	aret[17] = fgamma[20];
    }

    x73 = coor[15] - coor[11];
    y73 = coor[23] - coor[19];
    x74 = coor[15] - coor[12];
    y74 = coor[23] - coor[20];
    x43 = coor[12] - coor[11];
    y43 = coor[20] - coor[19];
    desur[0] = sqrt((x73 * 3. + x74) * (x73 * 3. + x74) + (y73 * 3. + y74) * (y73 * 3. + y74)) / 6.;
    desur[1] = sqrt((x74 * 3. + x73) * (x74 * 3. + x73) + (y74 * 3. + y73) * (y74 * 3. + y73)) / 6.;
    desur[2] = sqrt(x43 * x43 + y43 * y43) * 4. / 6.;
    be[5] += desur[0] * aret[12];
    be[6] += desur[0] * aret[13];
    be[7] += desur[1] * aret[14];
    be[8] += desur[1] * aret[15];
    be[13] += desur[2] * aret[16];
    be[14] += desur[2] * aret[17];

/*     -- arete 4 */
    if (norefs[12] != 0) {
	k = 4;
	j = 1;
	m = k + 4;
	xmi = coor[m + 8] - coor[k + 8];
	ymi = coor[m + 16] - coor[k + 16];
/* Computing 2nd power */
	d__1 = xmi;
/* Computing 2nd power */
	d__2 = ymi;
	arelon = sqrt(d__1 * d__1 + d__2 * d__2);
	xnu[0] = ymi / arelon;
	ynu[0] = -xmi / arelon;
	xjm = coor[j + 8] - coor[m + 8];
	yjm = coor[j + 16] - coor[m + 16];
/* Computing 2nd power */
	d__1 = xjm;
/* Computing 2nd power */
	d__2 = yjm;
	arelon = sqrt(d__1 * d__1 + d__2 * d__2);
	xnu[1] = yjm / arelon;
	ynu[1] = -xjm / arelon;
	xnu[2] = (xnu[0] + xnu[1]) * .5;
	ynu[2] = (ynu[0] + ynu[1]) * .5;
	aret[18] = fgamma[21] - pressi[10] * xnu[0];
	aret[19] = fgamma[22] - pressi[10] * ynu[0];
	aret[20] = fgamma[23] - pressi[11] * xnu[1];
	aret[21] = fgamma[24] - pressi[11] * ynu[1];
	aret[22] = fgamma[25] - pressi[12] * xnu[2];
	aret[23] = fgamma[26] - pressi[12] * ynu[2];
    } else {
	aret[18] = fgamma[21];
	aret[19] = fgamma[22];
	aret[20] = fgamma[23];
	aret[21] = fgamma[24];
	aret[22] = fgamma[25];
	aret[23] = fgamma[26];
    }

    x81 = coor[16] - coor[9];
    y81 = coor[24] - coor[17];
    x84 = coor[16] - coor[12];
    y84 = coor[24] - coor[20];
    x41 = coor[12] - coor[9];
    y41 = coor[20] - coor[17];
    desur[0] = sqrt((x84 * 3. + x81) * (x84 * 3. + x81) + (y84 * 3. + y81) * (y84 * 3. + y81)) / 6.;
    desur[1] = sqrt((x81 * 3. + x84) * (x81 * 3. + x84) + (y81 * 3. + y84) * (y81 * 3. + y84)) / 6.;
    desur[2] = sqrt(x41 * x41 + y41 * y41) * 4. / 6.;
    be[7] += desur[0] * aret[18];
    be[8] += desur[0] * aret[19];
    be[1] += desur[1] * aret[20];
    be[2] += desur[1] * aret[21];
    be[15] += desur[2] * aret[22];
    be[16] += desur[2] * aret[23];

} /* ets2q2c_ */

