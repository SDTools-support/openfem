#pragma once

#ifndef HYPER_H
#define HYPER_H
#include <stdio.h>
#include <string.h>
#include "mex.h"
#include "math.h"

#ifdef MatlabVER
#include "matrix.h"
#else
#include "stack-c.h"
#endif
#include "../mex/of_EltConst.h"



static double dI1dc[6]={1.,1.,1.,0.,0.,0.};
static double d2I2dcdc[36]={0.0,1.0,1.0,0.0,0.0,0.0,1.0,0.0,1.0,0.0,
            0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.5,0.0,0.0,
            0.0,0.0,0.0,0.0,-0.5,0.0,0.0,0.0,0.0,0.0,0.0,-0.5};

void elemcalc(struct EltConst* ECp,int jw, int Nw, int Nnode, double *e0, double *e, double *e1, double *F_ij, double *F1,
	     double *I, double *dI2dc, double *dI3dc, double *d2I3dcdc, int Mdef);

void EnPassiv(int *integ,double *constit,double *I,double *dWdI,double *d2WdI2);

void matPassiv( double *Sigma, double *dWdI, double *d2WdI2,
	       double *dI2dc, double *dI3dc, double *d2I3dcdc,double *d2wde2);

void EnVisco(double *constit,double *e0, double *e, double *Sigma, double *d2wvdep2);


#endif
