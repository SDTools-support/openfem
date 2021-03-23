/* tabaxd.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Subroutine */ int tabaxd_(i1, i2, a, b, c__, aux, iopt)
int32 *i1, *i2;
doublereal *a, *b, *c__, *aux;
int32 *iopt;
{
    /* System generated locals */
    int32 a_dim1, a_offset, aux_dim1, aux_offset;

    /* Local variables */
    extern /* Subroutine */ int tab2d_(), ab4d_(), ab5d_();

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  BUT : CALCUL DU PRODUIT MATRICIEL TRANSPOSEE [ A ] * [ B ] * [ A ] */
/*  ---   SI IOPT = 1 : [ C ]  =   0   + [ TA ] * [ B ] * [ A ] */
/*                         S               I1*I2     S     I2*I1 */
/*        SI IOPT > 1 : [ C ]  = [ C ] + [ TA ] * [ B ] * [ A ] */
/*                         S        S      I1*I2     S     I2*I1 */
/*       AVEC B ET C SYMETRIQUES , STOCKEES DE HAUT EN BAS ET DE LA */
/*       GAUCHE VERS LA DROITE. SEULES LES PARTIES TRIANGULAIRES */
/*       SUPERIEURES SONT STOCKEES. */
/*  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  PARAMETRES D'ENTREE : */
/*  ------------------- */
/*  I1    : NOMBRE DE COLONNES DE A */
/*  I2    : NOMBRE DE LIGNES DE A */
/*  A , B : MATRICES [ A ] I2*I1  ET [ B ] S */
/*  IOPT  : OPTION DE CALCUL */
/*  PARAMETRES RESULTATS : */
/*  -------------------- */
/*  C   : MATRICE [ C ] SYMETRIQUE */
/*  AUX : MATRICE AUXILLAIRE [ AUX ] I1*I2 */
/*  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  PROGRAMMEUR : A. HASSIM INRIA */
/*  ................................................................... */

/*  CALCUL DE [ AUX ] = TRANSP[ A ] * [ B ] */
/*                I1*I2         I1*I2     S */
    /* Parameter adjustments */
    aux_dim1 = *i1;
    aux_offset = aux_dim1 + 1;
    aux -= aux_offset;
    a_dim1 = *i2;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --b;
    --c__;

    /* Function Body */
    tab2d_(i1, i2, &a[a_offset], &b[1], &aux[aux_offset]);

    if (*iopt == 1) {

/*            CALCUL DE [ C ] = [ AUX ] * [ A ] */
/*                          S       I1*I2   I2*I1 */
	ab4d_(i1, i2, &aux[aux_offset], &a[a_offset], &c__[1]);

    } else {

/*            CALCUL DE [ C ] = [ C ] + [ AUX ] * [ A ] */
/*                          S       S       I1*I2   I2*I1 */
	ab5d_(i1, i2, &aux[aux_offset], &a[a_offset], &c__[1]);

    }
} /* tabaxd_ */

