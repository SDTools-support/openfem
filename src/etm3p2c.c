/* etm3p2c.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__10 = 10;
static int32 c__15 = 15;

/* Subroutine */ int etm3p2c_(coor, ro, iopt, ae)
doublereal *coor, *ro;
int32 *iopt;
doublereal *ae;
{
    /* Initialized data */

    static int32 ijt[30] = { 1,4,7,10,13,16,19,22,25,28,2,5,8,11,14,17,20,23,26,29,3,6,9,12,15,18,21,24,27,30 };
    static doublereal vp1[150]	/* was [10][15] */ = { -.125,-.125,-.125,-.125,.25,.25,.25,.25,.25,.25,.324516523092734,-.0750537196563629,-.0750537196563629,-.0750537196563629,.266380161832731,.0338347167927203,.266380161832731,.266380161832731,.0338347167927203,.0338347167927203,-.075053719656363,.324516523092734,-.0750537196563629,-.0750537196563629,.266380161832731,.266380161832731,.0338347167927203,.0338347167927203,.266380161832731,.0338347167927203,-.0750537196563627,-.0750537196563629,
	    .324516523092734,-.0750537196563629,.0338347167927203,.266380161832731,.266380161832731,.0338347167927203,.0338347167927203,.266380161832731,-.0750537196563627,-.0750537196563629,-.0750537196563629,.324516523092734,.0338347167927203,.0338347167927203,.0338347167927203,.266380161832731,.266380161832731,.266380161832731,-.0373192912588241,-.115257699028758,-.115257699028758,-.115257699028758,.0519589385132896,.409071857601743,.0519589385132896,.0519589385132896,.409071857601743,
	    .409071857601743,-.115257699028758,-.0373192912588239,-.115257699028758,-.115257699028758,.0519589385132895,.0519589385132894,.409071857601743,.409071857601743,.0519589385132894,.409071857601743,-.115257699028758,-.115257699028758,-.0373192912588239,-.115257699028758,.409071857601743,.0519589385132894,.0519589385132895,.409071857601743,.409071857601743,.0519589385132894,-.115257699028758,-.115257699028758,-.115257699028758,-.0373192912588239,.409071857601743,.409071857601743,
	    .409071857601743,.0519589385132895,.0519589385132894,.0519589385132894,-.05,-.05,-.05,-.05,.0127016653792583,.1,.1,.1,.1,.787298334620742,-.05,-.05,-.05,-.05,.1,.0127016653792583,.1,.787298334620742,.1,.1,-.0500000000000001,-.05,-.05,-.05,.787298334620742,.1,.1,.1,.1,.0127016653792583,-.0499999999999999,-.05,-.05,-.05,.0999999999999999,.787298334620742,.0999999999999999,.0127016653792583,.1,.1,-.05,-.05,-.05,-.05,.1,.1,.0127016653792583,.1,.787298334620742,.1,-.05,-.05,-.05,-.05,.1,
	    .1,.787298334620742,.1,.0127016653792583,.1 };
    static doublereal poids[15] = { .0197530873119831,.0119895139631698,.0119895139631698,.0119895139631698,.0119895139631698,.0115113678710454,.0115113678710454,.0115113678710454,.0115113678710454,.00881834235042334,.00881834235042334,.00881834235042334,.00881834235042334,.00881834235042334,.00881834235042334 };
    static struct {
	doublereal e_1[450];
	} equiv_2 = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1., -1., 1., 1., 0., -1., 0., -1., -1., -1., 0., 1., 0., 1., 0., 1., 1., -1.89634706336732, -1.89634706336732, -1.89634706336732, -.632115687789108, 0., 0., 0., -.632115687789108, 0., 0., 0., -.632115687789108, 2.52846275115643, -.367884312210892, -.367884312210892, .367884312210892, .367884312210892, 0., -.367884312210892, 2.52846275115643, -.367884312210892, -.367884312210892, -.367884312210892, 2.52846275115643, 
		.367884312210892, 0., .367884312210892, 0., .367884312210892, .367884312210892, .632115687789108, .632115687789108, .632115687789108, 1.89634706336732, 0., 0., 0., -.632115687789108, 0., 0., 0., -.632115687789108, -2.52846275115643, -2.89634706336732, -2.89634706336732, .367884312210892, 2.89634706336732, 0., -.367884312210892, 0., -.367884312210892, -.367884312210892, -.367884312210892, 0., .367884312210892, 0., 2.89634706336732, 0., .367884312210892, .367884312210892, 
		.632115687789108, .632115687789108, .632115687789108, -.632115687789108, 0., 0., 0., 1.89634706336732, 0., 0., 0., -.632115687789108, 0., -.367884312210892, -.367884312210892, 2.89634706336732, .367884312210892, 0., -2.89634706336732, -2.52846275115643, -2.89634706336732, -.367884312210892, -.367884312210892, 0., .367884312210892, 0., .367884312210892, 0., .367884312210892, 2.89634706336732, .632115687789108, .632115687789108, .632115687789108, -.632115687789108, 0., 0., 0., 
		-.632115687789108, 0., 0., 0., 1.89634706336732, 0., -.367884312210892, -.367884312210892, .367884312210892, .367884312210892, 0., -.367884312210892, 0., -.367884312210892, -2.89634706336732, -2.89634706336732, -2.52846275115643, 2.89634706336732, 0., .367884312210892, 0., 2.89634706336732, .367884312210892, .837523533955559, .837523533955559, .837523533955559, .27917451131852, 0., 0., 0., .27917451131852, 0., 0., 0., .27917451131852, -1.11669804527408, -1.27917451131852, 
		-1.27917451131852, 1.27917451131852, 1.27917451131852, 0., -1.27917451131852, -1.11669804527408, -1.27917451131852, -1.27917451131852, -1.27917451131852, -1.11669804527408, 1.27917451131852, 0., 1.27917451131852, 0., 1.27917451131852, 1.27917451131852, -.27917451131852, -.27917451131852, -.27917451131852, -.837523533955559, 0., 0., 0., .27917451131852, 0., 0., 0., .27917451131852, 1.11669804527408, -.162476466044441, -.162476466044441, 1.27917451131852, .162476466044441, 0., 
		-1.27917451131852, 0., -1.27917451131852, -1.27917451131852, -1.27917451131852, 0., 1.27917451131852, 0., .162476466044441, 0., 1.27917451131852, 1.27917451131852, -.27917451131852, -.27917451131852, -.27917451131852, .27917451131852, 0., 0., 0., -.837523533955559, 0., 0., 0., .27917451131852, 0., -1.27917451131852, -1.27917451131852, .162476466044441, 1.27917451131852, 0., -.162476466044441, 1.11669804527408, -.162476466044441, -1.27917451131852, -1.27917451131852, 0., 
		1.27917451131852, 0., 1.27917451131852, 0., 1.27917451131852, .162476466044441, -.27917451131852, -.27917451131852, -.27917451131852, .27917451131852, 0., 0., 0., .27917451131852, 0., 0., 0., -.837523533955559, 0., -1.27917451131852, -1.27917451131852, 1.27917451131852, 1.27917451131852, 0., -1.27917451131852, 0., -1.27917451131852, -.162476466044441, -.162476466044441, 1.11669804527408, .162476466044441, 0., 1.27917451131852, 0., .162476466044441, 1.27917451131852, 
		.774596669241483, .774596669241483, .774596669241483, -.774596669241483, 0., 0., 0., .774596669241483, 0., 0., 0., .774596669241483, 0., -.225403330758517, -.225403330758517, 1.77459666924148, .225403330758517, 0., -1.77459666924148, -1.54919333848297, -1.77459666924148, -1.77459666924148, -1.77459666924148, -1.54919333848297, 1.77459666924148, 0., .225403330758517, 0., 1.77459666924148, 1.77459666924148, -.774596669241483, -.774596669241483, -.774596669241483, -.774596669241483, 
		0., 0., 0., -.774596669241483, 0., 0., 0., .774596669241483, 1.54919333848297, -.225403330758517, -.225403330758517, .225403330758517, .225403330758517, 0., -.225403330758517, 1.54919333848297, -.225403330758517, -1.77459666924148, -1.77459666924148, 0., 1.77459666924148, 0., .225403330758517, 0., 1.77459666924148, .225403330758517, -.774596669241483, -.774596669241483, -.774596669241483, .774596669241483, 0., 0., 0., -.774596669241483, 0., 0., 0., -.774596669241483, 0., 
		-1.77459666924148, -1.77459666924148, .225403330758517, 1.77459666924148, 0., -.225403330758517, 1.54919333848297, -.225403330758517, -.225403330758517, -.225403330758517, 1.54919333848297, .225403330758517, 0., 1.77459666924148, 0., .225403330758517, .225403330758517, .774596669241483, .774596669241483, .774596669241483, .774596669241483, 0., 0., 0., .774596669241483, 0., 0., 0., -.774596669241483, -1.54919333848297, -1.77459666924148, -1.77459666924148, 1.77459666924148, 
		1.77459666924148, 0., -1.77459666924148, -1.54919333848297, -1.77459666924148, -.225403330758517, -.225403330758517, 0., .225403330758517, 0., 1.77459666924148, 0., .225403330758517, 1.77459666924148, .774596669241483, .774596669241483, .774596669241483, .774596669241483, 0., 0., 0., -.774596669241483, 0., 0., 0., .774596669241483, -1.54919333848297, -1.77459666924148, -1.77459666924148, .225403330758517, 1.77459666924148, 0., -.225403330758517, 0., -.225403330758517, 
		-1.77459666924148, -1.77459666924148, -1.54919333848297, 1.77459666924148, 0., 1.77459666924148, 0., 1.77459666924148, .225403330758517, -.774596669241483, -.774596669241483, -.774596669241483, -.774596669241483, 0., 0., 0., .774596669241483, 0., 0., 0., -.774596669241483, 1.54919333848297, -.225403330758517, -.225403330758517, 1.77459666924148, .225403330758517, 0., -1.77459666924148, 0., -1.77459666924148, -.225403330758517, -.225403330758517, 1.54919333848297, .225403330758517,
		 0., .225403330758517, 0., .225403330758517, 1.77459666924148 };


    /* Local variables */
#define vdpq2 ((doublereal *)&equiv_2)
    extern /* Subroutine */ int em3c2c_();
    static doublereal delta[15];
#define dp1 ((doublereal *)&equiv_2)
#define dp2 ((doublereal *)&equiv_2 + 150)
#define dp3 ((doublereal *)&equiv_2 + 300)

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  BUT   : CALCUL DE LA MATRICE ELEMENTAIRE ELASTIQUE DE RIGIDITE */
/*  ---     TETRAEDRE P2 ISOPARAMETRIQUE */
/*          formule d'integration a 15 points (Gauss) */

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
/*   (x y z)int : coordonnees des points d'integration sur */
/*              l'element courant */

/*  PARAMETRE DE SORTIE  : */
/*  -------------------- */
/*   AE     : MATRICE DE RIGIDITE . SYMETRIQUE */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  PROGRAMMEURS  : Marina Vidrascu 2001 */
/* ................................................................... */


    /* Parameter adjustments */
    --ae;
    coor -= 11;

    /* Function Body */

    em3c2c_(&c__10, &c__10, &coor[11], &coor[21], &coor[31], &c__15, ijt, poids, vp1, vdpq2, ro, &ae[1], delta);
} /* etm3p2c_ */

#undef dp3
#undef dp2
#undef dp1
#undef vdpq2

