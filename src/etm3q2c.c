/* etm3q2c.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "../mex/of_def.h"

/* Table of constant values */

static int32 c__20 = 20;
static int32 c__27 = 27;

/* Subroutine */ int etm3q2c_(coor, ro, iopt, ae)
doublereal *coor, *ro;
int32 *iopt;
doublereal *ae;
{
    /* Initialized data */

    static int32 ijt[60] = { 1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60 };
    static doublereal poids[27] = { .0214334705075446,.0342935528120713,.0214334705075446,.0342935528120713,.0548696844993141,.0342935528120713,.0214334705075446,.0342935528120713,.0214334705075446,.0342935528120713,.0548696844993141,.0342935528120713,.0548696844993141,.0877914951989026,.0548696844993141,.0342935528120713,.0548696844993141,.0342935528120713,.0214334705075446,.0342935528120713,.0214334705075446,.0342935528120713,.0548696844993141,.0342935528120713,.0214334705075446,
	    .0342935528120713,.0214334705075446 };
    static struct {
	doublereal e_1[1620];
	} equiv_6 = { -1.65205633616563, -1.65205633616563, -1.65205633616563, -.787298334620742, -.0549193338482967, -.0549193338482967, -.254919333848297, -.254919333848297, .0127016653792583, -.0549193338482967, -.787298334620742, -.0549193338482967, -.0549193338482967, -.0549193338482967, -.787298334620742, -.254919333848297, .0127016653792583, -.254919333848297, -.0520563361656317, -.0520563361656317, -.0520563361656317, .0127016653792583, -.254919333848297, -.254919333848297, 
		2.43935467078637, -.354919333848297, -.354919333848297, .354919333848297, .309838667696593, -.0450806661517033, .309838667696593, .354919333848297, -.0450806661517033, -.354919333848297, 2.43935467078637, -.354919333848297, -.354919333848297, -.354919333848297, 2.43935467078637, .354919333848297, -.0450806661517033, .309838667696593, .0450806661517033, .0450806661517033, .0393546707863734, -.0450806661517033, .354919333848297, .309838667696593, .309838667696593, 
		-.0450806661517033, .354919333848297, .0450806661517033, .0393546707863734, .0450806661517033, .0393546707863734, .0450806661517033, .0450806661517033, -.0450806661517033, .309838667696593, .354919333848297, -.432379000772445, -.587298334620742, -.587298334620742, .432379000772445, -.587298334620742, -.587298334620742, -.1, -.787298334620742, .0127016653792583, .1, -.787298334620742, .0127016653792583, .1, .0127016653792583, -.787298334620742, -.1, .0127016653792583, 
		-.787298334620742, -.032379000772445, -.187298334620742, -.187298334620742, .032379000772445, -.187298334620742, -.187298334620742, 0., -.887298334620742, -.887298334620742, .354919333848297, 1.37459666924148, -.2, 0., .887298334620742, -.112701665379258, -.354919333848297, 1.37459666924148, -.2, -.354919333848297, -.2, 1.37459666924148, .354919333848297, -.2, 1.37459666924148, .0450806661517033, .2, .174596669241483, -.0450806661517033, .2, .174596669241483, 0., 
		-.112701665379258, .887298334620742, .0450806661517033, .174596669241483, .2, 0., .112701665379258, .112701665379258, -.0450806661517033, .174596669241483, .2, .787298334620742, -.054919333848297, -.054919333848297, 1.65205633616563, -1.65205633616563, -1.65205633616563, .0549193338482967, -.787298334620742, -.0549193338482967, .254919333848297, -.254919333848296, .0127016653792583, .254919333848297, .0127016653792583, -.254919333848296, .0549193338482967, -.0549193338482967, 
		-.787298334620742, -.0127016653792583, -.254919333848297, -.254919333848297, .0520563361656317, -.0520563361656317, -.0520563361656317, -2.43935467078637, -.354919333848297, -.354919333848297, .354919333848297, 2.43935467078637, -.354919333848297, -.309838667696593, .354919333848297, -.0450806661517033, -.354919333848297, .309838667696593, -.0450806661517033, -.354919333848297, -.0450806661517033, .309838667696593, .354919333848297, -.354919333848297, 2.43935467078637, 
		.0450806661517033, .354919333848297, .309838667696593, -.0450806661517033, .0450806661517033, .0393546707863734, -.309838667696593, -.0450806661517033, .354919333848297, .0450806661517033, .309838667696593, .354919333848297, -.0393546707863734, .0450806661517033, .0450806661517033, -.0450806661517033, .0393546707863734, .0450806661517033, -.587298334620741, -.432379000772445, -.587298334620742, -.787298334620742, .1, .0127016653792583, -.787298334620742, -.1, .0127016653792583, 
		-.587298334620742, .432379000772445, -.587298334620742, .0127016653792583, .1, -.787298334620742, -.187298334620742, .032379000772445, -.187298334620742, -.187298334620742, -.032379000772445, -.187298334620742, .0127016653792583, -.1, -.787298334620742, 1.37459666924148, -.354919333848297, -.2, .887298334620742, 0., -.112701665379258, 1.37459666924148, .354919333848297, -.2, -.887298334620742, 0., -.887298334620742, -.2, -.354919333848297, 1.37459666924148, .2, -.0450806661517033,
		 .174596669241483, .2, .0450806661517033, .174596669241483, -.2, .354919333848297, 1.37459666924148, .174596669241483, -.0450806661517033, .2, .112701665379258, 0., .112701665379258, .174596669241483, .0450806661517033, .2, -.112701665379258, 0., .887298334620742, .1, .1, -.137298334620742, -.0999999999999999, .1, -.137298334620742, -.1, -.1, -.137298334620742, .1, -.1, -.137298334620742, .1, .1, -.637298334620742, -.1, .1, -.637298334620742, -.1, -.1, -.637298334620742, .1, -.1, 
		-.637298334620742, 0., -.887298334620742, -.5, .887298334620742, 0., -.5, 0., .887298334620742, -.5, -.887298334620742, 0., -.5, -.2, -.2, .774596669241483, .2, -.2, .774596669241483, .2, .2, .774596669241483, -.2, .2, .774596669241483, 0., -.112701665379258, .5, .112701665379258, 0., .5, 0., .112701665379258, .5, -.112701665379258, 0., .5, .787298334620742, .0999999999999998, .0127016653792577, .587298334620742, -.432379000772445, -.587298334620742, .587298334620742, 
		.432379000772445, -.587298334620742, .787298334620742, -.1, .0127016653792583, .187298334620742, .0323790007724451, -.187298334620742, -.0127016653792583, .1, -.787298334620742, -.0127016653792583, -.0999999999999999, -.787298334620741, .187298334620742, -.032379000772445, -.187298334620742, -1.37459666924148, -.354919333848297, -.2, .887298334620742, 0., -.887298334620742, -1.37459666924148, .354919333848297, -.2, -.887298334620742, 0., -.112701665379258, -.2, -.0450806661517033,
		 .174596669241483, .2, -.354919333848297, 1.37459666924148, .2, .354919333848297, 1.37459666924148, -.2, .0450806661517033, .174596669241483, -.174596669241483, -.0450806661517033, .2, .112701665379258, 0., .887298334620742, -.174596669241483, .0450806661517033, .2, -.112701665379258, 0., .112701665379258, -.0549193338482971, .787298334620742, -.0549193338482971, -.254919333848296, .254919333848297, .0127016653792583, -.787298334620742, .0549193338482967, -.0549193338482967, 
		-1.65205633616563, 1.65205633616563, -1.65205633616563, .0127016653792583, .254919333848297, -.254919333848297, -.0520563361656317, .0520563361656317, -.0520563361656317, -.254919333848297, -.0127016653792583, -.254919333848297, -.0549193338482967, .0549193338482967, -.787298334620742, .309838667696593, -.354919333848297, -.0450806661517033, .354919333848297, -.309838667696593, -.0450806661517033, 2.43935467078637, .354919333848297, -.354919333848297, -.354919333848297, 
		-2.43935467078637, -.354919333848297, -.0450806661517033, -.354919333848297, .309838667696593, .0450806661517033, -.0450806661517033, .0393546707863734, .354919333848297, .0450806661517033, .309838667696593, -.354919333848297, .354919333848297, 2.43935467078637, .0393546707863734, -.0450806661517033, .0450806661517033, .0450806661517033, -.0393546707863734, .0450806661517033, .309838667696593, .0450806661517033, .354919333848297, -.0450806661517033, -.309838667696593, 
		.354919333848297, .0999999999999997, .787298334620742, .0127016653792574, -.0999999999999998, .787298334620742, .0127016653792583, .432379000772445, .587298334620742, -.587298334620742, -.432379000772445, .587298334620742, -.587298334620742, .032379000772445, .187298334620742, -.187298334620742, -.032379000772445, .187298334620742, -.187298334620742, -.0999999999999999, -.0127016653792583, -.787298334620741, .1, -.0127016653792583, -.787298334620742, 0., -.887298334620742, 
		-.112701665379258, .354919333848297, -1.37459666924148, -.2, 0., .887298334620742, -.887298334620742, -.354919333848297, -1.37459666924148, -.2, -.0450806661517033, -.2, .174596669241483, .0450806661517033, -.2, .174596669241483, .354919333848297, .2, 1.37459666924148, -.354919333848297, .2, 1.37459666924148, 0., -.112701665379258, .112701665379258, .0450806661517033, -.174596669241483, .2, 0., .112701665379258, .887298334620742, -.0450806661517033, -.174596669241483, .2, 
		.254919333848296, .254919333848296, .0127016653792571, .0549193338482968, .787298334620742, -.0549193338482967, 1.65205633616563, 1.65205633616563, -1.65205633616563, .787298334620742, .054919333848297, -.0549193338482965, .0520563361656317, .0520563361656317, -.0520563361656314, -.0127016653792583, .254919333848297, -.254919333848297, .0549193338482968, .0549193338482968, -.787298334620742, .254919333848297, -.0127016653792583, -.254919333848297, -.309838667696593, 
		-.354919333848297, -.0450806661517031, .354919333848297, -2.43935467078637, -.354919333848297, -2.43935467078637, .354919333848297, -.354919333848297, -.354919333848297, -.309838667696593, -.0450806661517031, -.0450806661517033, -.0450806661517033, .0393546707863732, .0450806661517033, -.354919333848297, .309838667696593, .354919333848297, .354919333848297, 2.43935467078637, -.354919333848297, .0450806661517033, .309838667696593, -.0393546707863734, -.0450806661517033, 
		.0450806661517031, .0450806661517033, -.309838667696593, .354919333848297, -.309838667696593, .0450806661517033, .354919333848297, -.0450806661517033, -.0393546707863734, .0450806661517031, -.587298334620742, -.587298334620742, -.432379000772445, -.787298334620742, .0127016653792583, .1, -.187298334620742, -.187298334620742, .032379000772445, .0127016653792583, -.787298334620742, .1, -.587298334620742, -.587298334620742, .432379000772445, -.787298334620742, .0127016653792583, -.1,
		 -.187298334620742, -.187298334620742, -.032379000772445, .0127016653792583, -.787298334620742, -.1, 1.37459666924148, -.2, -.354919333848297, .2, .174596669241483, -.0450806661517033, .174596669241483, .2, -.0450806661517033, -.2, 1.37459666924148, -.354919333848297, -.887298334620742, -.887298334620742, 0., .887298334620742, -.112701665379258, 0., .112701665379258, .112701665379258, 0., -.112701665379258, .887298334620742, 0., 1.37459666924148, -.2, .354919333848297, .2, 
		.174596669241483, .0450806661517033, .174596669241483, .2, .0450806661517033, -.2, 1.37459666924148, .354919333848297, .1, -.137298334620742, .1, -.1, -.137298334620742, .1, -.1, -.637298334620742, .1, .1, -.637298334620742, .1, .1, -.137298334620742, -.1, -.1, -.137298334620742, -.1, -.1, -.637298334620742, -.1, .1, -.637298334620742, -.1, 0., -.5, -.887298334620742, .2, .774596669241483, -.2, 0., .5, -.112701665379258, -.2, .774596669241483, -.2, -.887298334620742, -.5, 0., 
		.887298334620742, -.5, 0., .112701665379258, .5, 0., -.112701665379258, .5, 0., 0., -.5, .887298334620742, .2, .774596669241483, .2, 0., .5, .112701665379258, -.2, .774596669241483, .2, .787298334620742, .0127016653792577, .0999999999999991, .587298334620742, -.587298334620742, -.432379000772445, -.0127016653792583, -.787298334620742, .1, .187298334620742, -.187298334620742, .032379000772445, .787298334620742, .0127016653792583, -.0999999999999998, .587298334620742, 
		-.587298334620742, .432379000772445, -.0127016653792583, -.787298334620742, -.0999999999999999, .187298334620742, -.187298334620742, -.032379000772445, -1.37459666924148, -.2, -.354919333848297, .2, 1.37459666924148, -.354919333848297, -.174596669241483, .2, -.0450806661517033, -.2, .174596669241483, -.0450806661517033, -.887298334620742, -.112701665379258, 0., .887298334620742, -.887298334620742, 0., .112701665379258, .887298334620742, 0., -.112701665379258, .112701665379258, 0.,
		 -1.37459666924148, -.2, .354919333848297, .2, 1.37459666924148, .354919333848297, -.174596669241483, .2, .0450806661517033, -.2, .174596669241483, .0450806661517033, -.137298334620742, .1, .1, -.637298334620742, .1, .1, -.637298334620742, -.1, .1, -.137298334620742, -.1, .1, -.137298334620742, .1, -.1, -.637298334620742, .1, -.1, -.637298334620742, -.1, -.1, -.137298334620742, -.1, -.1, .774596669241483, -.2, -.2, .5, 0., -.112701665379258, .774596669241483, .2, -.2, -.5, 0., 
		-.887298334620742, -.5, -.887298334620742, 0., .5, -.112701665379258, 0., .5, .112701665379258, 0., -.5, .887298334620742, 0., .774596669241483, -.2, .2, .5, 0., .112701665379258, .774596669241483, .2, .2, -.5, 0., .887298334620742, .25, .25, .25, -.25, .25, .25, -.25, -.25, .25, .25, -.25, .25, .25, .25, -.25, -.25, .25, -.25, -.25, -.25, -.25, .25, -.25, -.25, 0., -.5, -.5, .5, 0., -.5, 0., .5, -.5, -.5, 0., -.5, -.5, -.5, 0., .5, -.5, 0., .5, .5, 0., -.5, .5, 0., 0., -.5, .5, 
		.5, 0., .5, 0., .5, .5, -.5, 0., .5, .637298334620742, .0999999999999993, .0999999999999993, .137298334620742, .1, .0999999999999998, .137298334620742, -.0999999999999998, .1, .637298334620742, -.1, .1, .637298334620742, .1, -.1, .137298334620742, .1, -.1, .137298334620742, -.0999999999999998, -.0999999999999998, .637298334620742, -.1, -.1, -.774596669241483, -.2, -.2, .5, 0., -.887298334620742, -.774596669241483, .2, -.2, -.5, 0., -.112701665379258, -.5, -.112701665379258, 0., 
		.5, -.887298334620742, 0., .5, .887298334620742, 0., -.5, .112701665379258, 0., -.774596669241483, -.2, .2, .5, 0., .887298334620742, -.774596669241483, .2, .2, -.5, 0., .112701665379258, .0127016653792576, .787298334620742, .0999999999999995, -.187298334620742, .187298334620742, .032379000772445, -.787298334620741, -.0127016653792583, .1, -.587298334620742, .587298334620742, -.432379000772445, .0127016653792583, .787298334620742, -.1, -.187298334620742, .187298334620742, 
		-.032379000772445, -.787298334620742, -.0127016653792583, -.1, -.587298334620742, .587298334620742, .432379000772445, .174596669241483, -.2, -.0450806661517033, .2, -.174596669241483, -.0450806661517033, 1.37459666924148, .2, -.354919333848297, -.2, -1.37459666924148, -.354919333848297, -.112701665379258, -.887298334620742, 0., .112701665379258, -.112701665379258, 0., .887298334620742, .112701665379258, 0., -.887298334620742, .887298334620742, 0., .174596669241483, -.2, 
		.0450806661517033, .2, -.174596669241483, .0450806661517033, 1.37459666924148, .2, .354919333848297, -.2, -1.37459666924148, .354919333848297, .0999999999999992, .637298334620742, .0999999999999995, -.1, .637298334620742, .1, -.0999999999999998, .137298334620742, .1, .1, .137298334620742, .1, .1, .637298334620742, -.1, -.1, .637298334620742, -.1, -.0999999999999998, .137298334620742, -.0999999999999998, .1, .137298334620742, -.1, 0., -.5, -.112701665379258, .2, -.774596669241483, 
		-.2, 0., .5, -.887298334620742, -.2, -.774596669241483, -.2, -.112701665379258, -.5, 0., .112701665379258, -.5, 0., .887298334620742, .5, 0., -.887298334620742, .5, 0., 0., -.5, .112701665379258, .2, -.774596669241483, .2, 0., .5, .887298334620742, -.2, -.774596669241483, .2, .187298334620741, .187298334620741, .032379000772444, -.0127016653792584, .787298334620742, .1, .587298334620742, .587298334620742, -.432379000772445, .787298334620742, -.0127016653792583, .1, 
		.187298334620742, .187298334620742, -.0323790007724449, -.0127016653792584, .787298334620742, -.0999999999999999, .587298334620742, .587298334620742, .432379000772445, .787298334620742, -.0127016653792583, -.0999999999999999, -.174596669241483, -.2, -.0450806661517031, .2, -1.37459666924148, -.354919333848297, -1.37459666924148, .2, -.354919333848297, -.2, -.174596669241483, -.0450806661517031, -.112701665379258, -.112701665379258, 0., .112701665379258, -.887298334620742, 0., 
		.887298334620742, .887298334620742, 0., -.887298334620742, .112701665379258, 0., -.174596669241483, -.2, .0450806661517031, .2, -1.37459666924148, .354919333848297, -1.37459666924148, .2, .354919333848297, -.2, -.174596669241483, .0450806661517031, -.0549193338482972, -.0549193338482972, .787298334620742, -.254919333848297, .0127016653792583, .254919333848297, -.0520563361656317, -.0520563361656317, .0520563361656317, .0127016653792583, -.254919333848297, .254919333848297, 
		-1.65205633616563, -1.65205633616563, 1.65205633616563, -.787298334620742, -.0549193338482967, .0549193338482967, -.254919333848297, -.254919333848297, -.0127016653792583, -.0549193338482967, -.787298334620742, .0549193338482967, .309838667696593, -.0450806661517033, -.354919333848297, .0450806661517033, .0393546707863734, -.0450806661517033, .0393546707863734, .0450806661517033, -.0450806661517033, -.0450806661517033, .309838667696593, -.354919333848297, -.354919333848297, 
		-.354919333848297, -2.43935467078637, .354919333848297, -.0450806661517033, -.309838667696593, .0450806661517033, .0450806661517033, -.0393546707863734, -.0450806661517033, .354919333848297, -.309838667696593, 2.43935467078637, -.354919333848297, .354919333848297, .354919333848297, .309838667696593, .0450806661517033, .309838667696593, .354919333848297, .0450806661517033, -.354919333848297, 2.43935467078637, .354919333848297, .0999999999999995, .0127016653792574, .787298334620742,
		 -.0999999999999998, .0127016653792583, .787298334620742, -.032379000772445, -.187298334620742, .187298334620742, .032379000772445, -.187298334620742, .187298334620742, -.432379000772446, -.587298334620742, .587298334620742, .432379000772445, -.587298334620742, .587298334620742, -.1, -.787298334620742, -.0127016653792583, .1, -.787298334620742, -.0127016653792583, 0., -.112701665379258, -.887298334620742, .0450806661517033, .174596669241483, -.2, 0., .112701665379258, 
		-.112701665379258, -.0450806661517033, .174596669241483, -.2, -.354919333848297, -.2, -1.37459666924148, .354919333848297, -.2, -1.37459666924148, .0450806661517033, .2, -.174596669241483, -.0450806661517033, .2, -.174596669241483, 0., -.887298334620742, .887298334620742, .354919333848297, 1.37459666924148, .2, 0., .887298334620742, .112701665379258, -.354919333848297, 1.37459666924148, .2, .254919333848296, .0127016653792569, .254919333848296, .0549193338482968, 
		-.0549193338482965, .787298334620742, -.0127016653792583, -.254919333848297, .254919333848297, .0520563361656317, -.0520563361656312, .0520563361656317, .787298334620741, -.0549193338482969, .0549193338482966, 1.65205633616563, -1.65205633616563, 1.65205633616563, .0549193338482967, -.787298334620741, .0549193338482967, .254919333848297, -.254919333848297, -.0127016653792583, -.309838667696593, -.0450806661517031, -.354919333848297, .0450806661517033, .309838667696593, 
		-.354919333848297, -.0393546707863734, .0450806661517031, -.0450806661517033, -.0450806661517033, .0393546707863733, -.0450806661517033, -.354919333848297, -.0450806661517031, -.309838667696593, .354919333848297, -.354919333848297, -2.43935467078637, .0450806661517033, .354919333848297, -.309838667696593, -.0450806661517033, .0450806661517031, -.0393546707863734, -2.43935467078637, -.354919333848297, .354919333848297, .354919333848297, 2.43935467078637, .354919333848297, 
		-.309838667696593, .354919333848297, .0450806661517033, -.354919333848297, .309838667696593, .0450806661517033, .0127016653792578, .0999999999999995, .787298334620742, -.187298334620742, .032379000772445, .187298334620742, -.187298334620742, -.032379000772445, .187298334620742, .0127016653792583, -.1, .787298334620742, -.587298334620742, -.432379000772445, .587298334620742, -.787298334620742, .1, -.0127016653792583, -.787298334620742, -.1, -.0127016653792583, -.587298334620742, 
		.432379000772445, .587298334620742, .174596669241483, -.0450806661517033, -.2, .112701665379258, 0., -.112701665379258, .174596669241483, .0450806661517033, -.2, -.112701665379258, 0., -.887298334620742, -.2, -.354919333848297, -1.37459666924148, .2, -.0450806661517033, -.174596669241483, .2, .0450806661517033, -.174596669241483, -.2, .354919333848297, -1.37459666924148, 1.37459666924148, -.354919333848297, .2, .887298334620742, 0., .112701665379258, 1.37459666924148, 
		.354919333848297, .2, -.887298334620742, 0., .887298334620742, .0999999999999996, .0999999999999995, .637298334620742, -.1, .0999999999999999, .637298334620742, -.1, -.1, .637298334620742, .1, -.0999999999999999, .637298334620742, .0999999999999998, .0999999999999999, .137298334620742, -.0999999999999999, .1, .137298334620742, -.0999999999999999, -.0999999999999999, .137298334620742, .1, -.1, .137298334620742, 0., -.112701665379258, -.5, .112701665379258, 0., -.5, 0., 
		.112701665379258, -.5, -.112701665379258, 0., -.5, -.2, -.2, -.774596669241483, .2, -.2, -.774596669241483, .2, .2, -.774596669241483, -.2, .2, -.774596669241483, 0., -.887298334620742, .5, .887298334620742, 0., .5, 0., .887298334620742, .5, -.887298334620742, 0., .5, .187298334620741, .032379000772444, .187298334620741, -.0127016653792584, .1, .787298334620742, -.0127016653792584, -.1, .787298334620742, .187298334620742, -.0323790007724447, .187298334620742, .787298334620741, 
		.0999999999999999, -.0127016653792583, .587298334620742, -.432379000772445, .587298334620742, .587298334620742, .432379000772446, .587298334620742, .787298334620742, -.0999999999999999, -.0127016653792583, -.174596669241483, -.0450806661517031, -.2, .112701665379258, 0., -.887298334620742, -.174596669241483, .0450806661517031, -.2, -.112701665379258, 0., -.112701665379258, -.2, -.0450806661517031, -.174596669241483, .2, -.354919333848297, -1.37459666924148, .2, .354919333848297, 
		-1.37459666924148, -.2, .0450806661517031, -.174596669241483, -1.37459666924148, -.354919333848297, .2, .887298334620742, 0., .887298334620742, -1.37459666924148, .354919333848297, .2, -.887298334620742, 0., .112701665379258, .0127016653792573, .254919333848296, .254919333848296, -.0520563361656314, .0520563361656317, .0520563361656317, -.254919333848297, -.0127016653792583, .254919333848297, -.0549193338482966, .0549193338482969, .787298334620742, -.0549193338482968, 
		.787298334620742, .0549193338482967, -.254919333848297, .254919333848297, -.0127016653792583, -.787298334620741, .0549193338482967, .0549193338482967, -1.65205633616563, 1.65205633616563, 1.65205633616563, .0393546707863734, -.0450806661517033, -.0450806661517033, .0450806661517031, -.0393546707863734, -.0450806661517033, .309838667696593, .0450806661517033, -.354919333848297, -.0450806661517031, -.309838667696593, -.354919333848297, -.0450806661517031, -.354919333848297, 
		-.309838667696593, .0450806661517031, -.0450806661517033, -.0393546707863734, .354919333848297, .0450806661517033, -.309838667696593, -.354919333848297, .354919333848297, -2.43935467078637, .309838667696593, -.354919333848297, .0450806661517033, .354919333848297, -.309838667696593, .0450806661517033, 2.43935467078637, .354919333848297, .354919333848297, -.354919333848297, -2.43935467078637, .354919333848297, .0323790007724443, .187298334620741, .187298334620741, 
		-.0323790007724447, .187298334620742, .187298334620742, -.1, -.0127016653792582, .787298334620742, .1, -.0127016653792582, .787298334620742, .0999999999999999, .787298334620741, -.0127016653792585, -.0999999999999999, .787298334620742, -.0127016653792583, .432379000772446, .587298334620742, .587298334620742, -.432379000772445, .587298334620741, .587298334620741, 0., -.112701665379258, -.112701665379258, .0450806661517031, -.174596669241483, -.2, 0., .112701665379258, 
		-.887298334620742, -.0450806661517031, -.174596669241483, -.2, -.0450806661517031, -.2, -.174596669241483, .0450806661517031, -.2, -.174596669241483, .354919333848297, .2, -1.37459666924148, -.354919333848297, .2, -1.37459666924148, 0., -.887298334620742, .112701665379258, .354919333848297, -1.37459666924148, .2, 0., .887298334620742, .887298334620742, -.354919333848297, -1.37459666924148, .2, .0520563361656308, .052056336165631, .0520563361656301, -.0127016653792582, 
		.254919333848297, .254919333848296, .0549193338482965, .0549193338482967, .787298334620742, .254919333848297, -.012701665379258, .254919333848297, .254919333848297, .254919333848296, -.0127016653792584, .0549193338482969, .787298334620742, .0549193338482967, 1.65205633616563, 1.65205633616563, 1.65205633616563, .787298334620742, .0549193338482965, .0549193338482965, -.0393546707863734, -.0450806661517031, -.0450806661517031, .0450806661517031, -.309838667696594, -.354919333848297,
		 -.309838667696594, .0450806661517031, -.354919333848297, -.0450806661517031, -.0393546707863734, -.0450806661517031, -.0450806661517031, -.0450806661517031, -.0393546707863734, .0450806661517031, -.354919333848297, -.309838667696594, .354919333848297, .354919333848297, -2.43935467078637, -.354919333848297, .0450806661517031, -.309838667696594, -.309838667696594, -.354919333848297, .0450806661517031, .354919333848297, -2.43935467078637, .354919333848297, -2.43935467078637, 
		.354919333848297, .354919333848297, -.354919333848297, -.309838667696594, .0450806661517031 };

    static struct {
	doublereal e_1[540];
	} equiv_9 = { .226189500386223, -.108729833462074, -.0312701665379258, -.108729833462074, -.108729833462074, -.0312701665379258, -.00618950038622251, -.0312701665379258, .314919333848297, .04, .04, .314919333848297, .314919333848297, .04, .00508066615170332, .04, .04, .00508066615170332, .00508066615170332, .04, -.177459666924148, -.177459666924148, -.1, -.1, -.1, -.1, -.0225403330758517, -.0225403330758517, .787298334620742, .177459666924148, .1, .177459666924148, .177459666924148, 
		.177459666924148, .0225403330758517, .0225403330758517, .1, .0225403330758517, .0127016653792583, .0225403330758517, -.108729833462074, .226189500386223, -.108729833462074, -.0312701665379258, -.0312701665379258, -.108729833462074, -.0312701665379258, -.00618950038622251, .314919333848297, .314919333848297, .04, .04, .04, .314919333848297, .04, .00508066615170332, .04, .04, .00508066615170332, .00508066615170332, -.177459666924148, -.1, -.1, -.177459666924148, -.1, 
		-.0225403330758517, -.0225403330758517, -.1, .177459666924148, .1, .177459666924148, .787298334620742, .177459666924148, .0225403330758517, .0225403330758517, .177459666924148, .0225403330758517, .0127016653792583, .0225403330758517, .1, -.271824583655185, -.271824583655185, -.271824583655185, -.271824583655185, -.0781754163448146, -.0781754163448146, -.0781754163448146, -.0781754163448146, .443649167310371, .443649167310371, .443649167310371, .443649167310371, .1, .1, .1, .1, 
		.0563508326896291, .0563508326896291, .0563508326896291, .0563508326896291, -.1, -.177459666924148, -.177459666924148, -.0999999999999998, -.0225403330758517, -.1, -.1, -.0225403330758516, .177459666924148, .787298334620742, .177459666924148, .1, .0225403330758517, .177459666924148, .177459666924148, .0225403330758517, .0225403330758517, .1, .0225403330758517, .0127016653792583, -.108729833462074, -.0312701665379258, -.108729833462074, .226189500386223, -.0312701665379258, 
		-.00618950038622251, -.0312701665379258, -.108729833462074, .04, .04, .314919333848297, .314919333848297, .04, .00508066615170332, .04, .314919333848297, .00508066615170332, .00508066615170332, .04, .04, -.1, -.0999999999999999, -.177459666924148, -.177459666924148, -.0225403330758517, -.0225403330758516, -.1, -.1, .1, .177459666924148, .787298334620742, .177459666924148, .0225403330758517, .0225403330758517, .177459666924148, .177459666924148, .0127016653792583, 
		.0225403330758517, .1, .0225403330758517, -.0312701665379259, -.108729833462074, .226189500386222, -.108729833462074, -.00618950038622246, -.0312701665379258, -.108729833462074, -.0312701665379258, .0400000000000001, .314919333848297, .314919333848297, .0400000000000001, .00508066615170331, .04, .314919333848297, .04, .00508066615170327, .04, .04, .00508066615170327, -.177459666924148, -.1, -.0225403330758517, -.1, -.177459666924148, -.1, -.0225403330758517, -.1, .177459666924148,
		 .0225403330758517, .0225403330758517, .177459666924148, .787298334620742, .1, .0127016653792583, .1, .177459666924148, .0225403330758517, .0225403330758517, .177459666924148, -.271824583655186, -.271824583655185, -.0781754163448146, -.0781754163448145, -.271824583655185, -.271824583655185, -.0781754163448146, -.0781754163448146, .443649167310371, .1, .0563508326896291, .1, .443649167310371, .443649167310371, .0563508326896291, .0563508326896291, .443649167310371, .1, 
		.0563508326896291, .1, -.1, -.177459666924148, -.1, -.0225403330758516, -.0999999999999998, -.177459666924148, -.1, -.0225403330758517, .177459666924148, .177459666924148, .0225403330758517, .0225403330758517, .1, .787298334620742, .1, .0127016653792583, .177459666924148, .177459666924148, .0225403330758517, .0225403330758517, -.271824583655185, -.0781754163448146, -.0781754163448146, -.271824583655185, -.271824583655185, -.0781754163448146, -.0781754163448146, -.271824583655185, 
		.1, .0563508326896291, .1, .443649167310371, .443649167310371, .0563508326896291, .0563508326896291, .443649167310371, .1, .0563508326896291, .1, .443649167310371, -.25, -.25, -.25, -.25, -.25, -.25, -.25, -.25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, .25, -.0781754163448151, -.271824583655185, -.271824583655185, -.0781754163448146, -.0781754163448146, -.271824583655185, -.271824583655185, -.0781754163448146, .1, .443649167310371, .1, .0563508326896291, 
		.0563508326896291, .443649167310371, .443649167310371, .0563508326896291, .1, .443649167310371, .1, .0563508326896291, -.1, -.0225403330758516, -.1, -.177459666924148, -.1, -.0225403330758517, -.1, -.177459666924148, .0225403330758517, .0225403330758517, .177459666924148, .177459666924148, .1, .0127016653792583, .1, .787298334620742, .0225403330758517, .0225403330758517, .177459666924148, .177459666924148, -.078175416344815, -.0781754163448146, -.271824583655185, 
		-.271824583655185, -.0781754163448146, -.0781754163448146, -.271824583655185, -.271824583655185, .0563508326896291, .1, .443649167310371, .1, .0563508326896291, .0563508326896291, .443649167310371, .443649167310371, .0563508326896291, .1, .443649167310371, .1, -.0225403330758523, -.0999999999999998, -.177459666924149, -.0999999999999999, -.0225403330758515, -.1, -.177459666924148, -.1, .0225403330758513, .177459666924148, .177459666924148, .0225403330758513, .0127016653792583, .1,
		 .787298334620742, .1, .0225403330758516, .177459666924148, .177459666924148, .0225403330758516, -.108729833462074, -.0312701665379258, -.0061895003862225, -.0312701665379258, .226189500386223, -.108729833462074, -.0312701665379258, -.108729833462074, .04, .00508066615170332, .00508066615170332, .04, .314919333848297, .04, .00508066615170332, .04, .314919333848297, .04, .04, .314919333848297, -.1, -.1, -.0225403330758517, -.0225403330758516, -.177459666924148, -.177459666924148, 
		-.1, -.1, .1, .0225403330758517, .0127016653792583, .0225403330758517, .177459666924148, .177459666924148, .0225403330758517, .0225403330758517, .787298334620742, .177459666924148, .1, .177459666924148, -.0312701665379262, -.108729833462074, -.0312701665379258, -.00618950038622246, -.108729833462074, .226189500386222, -.108729833462074, -.0312701665379258, .04, .04, .00508066615170333, .00508066615170332, .04, .314919333848297, .04, .00508066615170333, .314919333848297, 
		.314919333848297, .04, .04, -.1, -.0225403330758516, -.0225403330758517, -.1, -.177459666924148, -.1, -.1, -.177459666924148, .0225403330758517, .0127016653792583, .0225403330758517, .1, .177459666924148, .0225403330758517, .0225403330758517, .177459666924148, .177459666924148, .1, .177459666924148, .787298334620742, -.0781754163448149, -.0781754163448146, -.0781754163448146, -.0781754163448146, -.271824583655185, -.271824583655185, -.271824583655185, -.271824583655185, 
		.0563508326896291, .0563508326896291, .0563508326896291, .0563508326896291, .1, .1, .1, .1, .443649167310371, .443649167310371, .443649167310371, .443649167310371, -.0225403330758526, -.0999999999999998, -.1, -.0225403330758516, -.0999999999999999, -.177459666924149, -.177459666924148, -.1, .0225403330758513, .1, .0225403330758516, .0127016653792583, .0225403330758513, .177459666924148, .177459666924148, .0225403330758516, .177459666924148, .787298334620742, .177459666924148, .1, 
		-.0312701665379264, -.00618950038622245, -.0312701665379258, -.108729833462074, -.108729833462074, -.0312701665379258, -.108729833462074, .226189500386222, .00508066615170332, .00508066615170333, .04, .0399999999999998, .0399999999999998, .00508066615170333, .04, .314919333848297, .04, .04, .314919333848297, .314919333848297, -.022540333075853, -.0225403330758515, -.1, -.0999999999999998, -.1, -.0999999999999999, -.177459666924148, -.177459666924149, .0127016653792583, 
		.0225403330758516, .1, .0225403330758516, .0225403330758516, .0225403330758516, .177459666924148, .177459666924148, .1, .177459666924148, .787298334620742, .177459666924148, -.00618950038622335, -.0312701665379256, -.108729833462074, -.0312701665379258, -.0312701665379254, -.108729833462075, .226189500386222, -.108729833462075, .00508066615170311, .04, .04, .00508066615170311, .00508066615170266, .0400000000000005, .314919333848296, .0400000000000005, .04, .314919333848297, 
		.314919333848297, .04 };


    /* Local variables */
#define vdpq2 ((doublereal *)&equiv_6)
    extern /* Subroutine */ int em3c2c_();
    static doublereal delta[27];
#define dp1 ((doublereal *)&equiv_6)
#define dp2 ((doublereal *)&equiv_6 + 240)
#define dp3 ((doublereal *)&equiv_6 + 480)
#define dp4 ((doublereal *)&equiv_6 + 720)
#define dp5 ((doublereal *)&equiv_6 + 960)
#define dp6 ((doublereal *)&equiv_6 + 1200)
#define dp7 ((doublereal *)&equiv_6 + 1440)
#define va1 ((doublereal *)&equiv_9)
#define va2 ((doublereal *)&equiv_9 + 180)
#define va3 ((doublereal *)&equiv_9 + 360)
#define vp1 ((doublereal *)&equiv_9)

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  BUT   : CALCUL DE LA MATRICE ELEMENTAIRE ELASTIQUE DE MASSE */
/*  ---     Hexaedre Q2 ISOPARAMETRIQUE */
/*          formule d'integration a 27 points, (Gauss) (3*3*3) */

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
/*   AE     : MATRICE DE Masse . SYMETRIQUE */


/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  PROGRAMMEURS  : Marina Vidrascu 2001 */
/* ................................................................... */


    /* Parameter adjustments */
    --ae;
    coor -= 21;

    /* Function Body */

    em3c2c_(&c__20, &c__20, &coor[21], &coor[41], &coor[61], &c__27, ijt, poids, vp1, vdpq2, ro, &ae[1], delta);
} /* etm3q2c_ */

#undef vp1
#undef va3
#undef va2
#undef va1
#undef dp7
#undef dp6
#undef dp5
#undef dp4
#undef dp3
#undef dp2
#undef dp1
#undef vdpq2


