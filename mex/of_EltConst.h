
#pragma once
#include "../mex/of_def.h"

/* int values should be first for easier access in Matlab Callbacks
   where pointers may by 32 or 64 bits */

struct EltConst {
    int v1x; /* column of v1x in NodeE */
    int v3x; /* column of v3x in NodeE */
    int t;   /* column of t in NodeE (temperature)*/
    int h;   /* column of h  thickness */
    double* nodeE;
    double* defE;
    double* Be;
    double* ke;
    double* bas;
    double* J;
    double* jdet;
    double* constit; /* material data after possible manimulation */
    double* constit0;/* read only material data */
    double* gstate;
    double* NDN;
    double* InfoAtNode;
    int* integ;
    int* giNodePos; /* group iNodePos for infoAtNode */
    int* NodeET;
};

struct GroupFields {
    double* N;
    double* Nr;
    double* Ns;
    double* Nt; 
    double* w;
    double* def;
    double* defi;
    double* RHS;
    double* ke;
    double* CTable;  /* interpolation table */
    int* rule_terms;
    int* topo;
    int* VectMap;
    int Nshape;
	int NBe;
	int NdefE;
	int Ncondense;
};


OF_EXPORT void GetpEC(struct EltConst, int *);
OF_EXPORT void SetpEC(struct EltConst *, int *);

