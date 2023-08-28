#pragma once
#define UseModulef   1

struct GroupFields initGroup(struct GroupFields GF, const mxArray *mNode, const mxArray *mEC, int *GSize, double* *NeedFree);
struct EltConst  initInfoAtNode(const mxArray *infoAtNode,const mxArray *mEC,struct EltConst EC, struct GroupFields GF, int *GSize, double* *NeedFree);
void AssembleSparse(mwIndex* ir, mwIndex* jc, double* pr, 
		   int NDDL, int IsSymVal, int IsSymK, int* keind, int* elmap, 
	           mwIndex* vir, mwIndex* vjc, double* val);

extern void x_k_x (const double* x, const double* y, double* result, 
	double* offset, int type, int YM,int YN);

extern void basis_clean(double* z);
extern void basis (const double* x, const double* y, double* z);

void buildNodeE(struct EltConst* ECp, double *node, int *GSize,
	   int *NodePos,  int jElt, int Nmnode);
void Transform33(double *IN, double *T, double *S);
void Transform33s(double *IN, double *T, double *S);
void TransformLambda(double *IN, double *T, double *S);

char* pre_cvs ();
mxArray* pre_cvs2();
void* ofMalloc(int len);
void ofFree(void* in);
void* NeedFreeAlloc(double** NeedFree,int pos, int len);
void* FieldPointer(const mxArray *field, char* name);

void AssembleMatVec(struct EltConst* EC,struct GroupFields GF, int* opt,
        int StrategyType, int Mk, int* CurDofPos,
        mwIndex* o_ir, mwIndex* o_jc, double* o_pr,
        int* elmap, int DofPerElt, int Mdef,int Ndef,double *Ener,
        int* cEGI, int jElt,int Mener);

void NDNSwitch(int type, struct GroupFields GF, struct EltConst* ECp, 
           int Nw, int Nnode, int Nfield);
void constitInterp(struct GroupFields GF,struct EltConst* ECp, int NDNoff, int Nnode, int *GSize);
void elt_condense(int i1, double* ke, double* kr, 
                          double* me, double* mr, int M);
#if UseModulef==1
int edemit4_(double*);
#endif
void   InitThreadiPos(int* ThreadiPos, int* GSize, int Mk);
void lin_multi(struct GroupFields GF,struct EltConst* ECp, int Nterms,int Mrule,int Ndef,int jElt,
	  int *GSize,int Mk,int *point,int *rrule);
void nonlin_elas(struct GroupFields GF,struct EltConst* ECp, int *GSize,int jElt);

void Mecha3DInteg( struct EltConst* ECp, double* F, double* d2wde2, double* Sigma, 
             double* w,int Nnode,int Ndof,int Nw, int jw);
void Mecha3DIntegH(struct EltConst* ECp, double *F, double *F1, double* d2wde2, 
          double* d2wvde2, double* Sigma,
	      double* w, int Nnode,int Ndof,int Nw, int jw);

mwIndex mesh_meshGraph( mwSize *p_nnz, mwIndex **p_prow, mwIndex **p_icol,
		      mwSize nNod, int nGr, int mConn,
		      mwIndex *EGroup, mwIndex *nEl,
		      mwIndex *DofPerElt, mwIndex *conn );
int mesh_dofGraph( mwSize *p_dnnz, mwIndex **p_dprow, mwIndex **p_dicol,
		     int nNod, mwIndex *nprow, mwIndex *nicol,
		     int nEq, int *dpn, int *dofOffset, int *eq );

void Hyper_main(struct GroupFields GF,struct EltConst* ECp, int Nnode,
         int Ndof,int Nw,int Mdef, int Ndef,int *CurDofPos);
int io_switch(char* cam,char* ierr,int* pointers, int* integ, double* constit, 
  double* node,double*  estate,double* defe, 
        double* eltconst,double*  out,double*  out1, double* out2,double*  infoatnode);

void F_pressure(struct GroupFields GF,struct EltConst* ECp, int Nnode,int Nw,int *point);
void of_time_LinInterp(double* table, double * val, double* last,
		       double* out, mwSize M, mwSize Nc, mwSize Npoints, mwSize nlhs,
               bool ti, double* tablei, double* outi);

