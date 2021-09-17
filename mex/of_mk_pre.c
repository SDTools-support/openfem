#include "mex.h" 
#include "math.h" 
#include <string.h>
#ifdef OFOMP
#include <omp.h>
#endif

#include "../mex/of_EltConst.h"
#include "../mex/of_mk_pre.h"


void ofFree(void* in) {
/* mexPrintf("Free %p\n",in);*/
   free(in);
 }

void* ofMalloc(int len) {
	void* out;
	out=malloc(len);
	return(out);
}

/*----------------------------------------------------------------------- */
void* NeedFreeAlloc(double** NeedFree,int pos, int len) {
 int j1;
#if 0
 if (pos<0) { /* -neg do free */
   pos=-pos;
   for (j1=0;j1<-pos;j1++) {
    if (NeedFree[j1]!=NULL) {
        mxFree(NeedFree[j1]); NeedFree[j1]=NULL;
    }
   }
  } else {
    if (NeedFree[pos]!=NULL) mxFree(NeedFree[pos]);
    NeedFree[pos]=mxCalloc(len,sizeof(double));
    return NeedFree[pos];
  }
#else 
 if (pos<0) { /* -neg do free */
   for (j1=0;j1<-pos;j1++) {
    if (NeedFree[j1]!=NULL) {
        ofFree(NeedFree[j1]); NeedFree[j1]=NULL; 
    }
   }
   return NULL;
  } else {
    if (NeedFree[pos]!=NULL) free(NeedFree[pos]);
    NeedFree[pos]=(double*)ofMalloc(len*sizeof(double));
    return NeedFree[pos];
  }
#endif  
}

/*----------------------------------------------------------------------- */
   /* GSize : [0] Nnode (nodes per element)
      [1] NDofPerField (rows in NDN)
      [2] Nw (number of integrations points in .w)
      [3] MInfoAtNode (number of nodes in InfoAtNode)
      [4] NInfoAtNode (number of fields)
      [5] NState (number of states at each integration point for gstate)
      [6] M_constit (number of rows in constit)
     DTable[0] gstate
     DTable[1] CTable (constitutive law interpolation tables
  */
void   InitThreadiPos(int* ThreadiPos, int* GSize, int Mk) {
  int Nnode, Nw;
  Nnode=GSize[0];
  Nw=GSize[2];
     ThreadiPos[0]=Nnode*(4+GSize[4]); /* Nnode*(4+NInfoAtNode) */
     ThreadiPos[1]=GSize[1]*Nw*4;
     ThreadiPos[2]=Nw;
     ThreadiPos[3]=9*Nw;
     ThreadiPos[4]=9*Nw;
     ThreadiPos[5]=Mk*Mk;
	 ThreadiPos[6]=0; /* no longer used */
	 ThreadiPos[7]=0; /* no longer used */
}

/*----------------------------------------------------------------------- */
/* #constitInterp Interpolation of a constit field
 * sdtweb m_elastic ctable 
 */
void constitInterp(struct GroupFields GF,struct EltConst* ECp, int NDNoff, int Nnode, int *GSize) 
{
	/* interpolate values */
      double     *t1,*t2,*X,val[1];
      int        M,n1,i2, j1,j2,offset;

      t1=GF.CTable; n1=(int)t1[0];if (ECp[0].constit==NULL) return;
      for (j1=0;j1<n1;j1++) {
       t2=t1+j1*7+1;M=(int)t2[4];i2=Nnode*((int)t2[5]-1);
       val[0]=0;for (j2=0;j2<Nnode;j2++) {val[0]+=ECp[0].NDN[j2+NDNoff]*ECp[0].nodeE[j2+i2];}
       offset=(int)t2[6]-1;
       if ((int)t2[3]==-1) { /* constit interpolated with InfoAtNode
          initialized in p_solid ctableGen, tested in t_thermal
        */ 
         ECp[0].constit[offset]=val[0]; 
		 /* if (offset>2) {mexPrintf("%g %g %g %g\n",ECp[0].nodeE[0],ECp[0].nodeE[8],ECp[0].nodeE[16],val[0]/210e9);}*/
	   } else if (t2[3]==-2.1) { /* local orient strain only 3x3 */
         offset=(int)t2[0]-1;
		 Transform33s(ECp[0].constit0+offset,ECp[0].bas+9*NDNoff/Nnode,ECp[0].constit+offset);
	   } else if ((int)t2[3]==-2) { /* local orient 3x3 used by p_zt */
         offset=(int)t2[0]-1;
		 Transform33(ECp[0].constit0+offset,ECp[0].bas+9*NDNoff/Nnode,ECp[0].constit+offset);
	   } else if ((int)t2[3]==-3) { /* local mecha tensor transform */
         offset=(int)t2[0]-1;
		 TransformLambda(ECp[0].constit0+offset,ECp[0].bas+9*NDNoff/Nnode,ECp[0].constit+offset);
       } else { /* val is field at node : temp, ... */
           X=t1+(int)t2[3]; /* table [X;Y] */
           of_time_LinInterp(X,val,t2,ECp[0].constit+offset,M,1,1,1,NULL,NULL);
       }
       /* mexPrintf("%g %g %g %g (%i)\n",X[0],X[1],X[2],X[3],M);
       mexPrintf(" %i %.2f i2=%i (%.3f)=%.3f\n",j1,val[0],i2,nodeE[i2],constit[(int)t2[6]-1]);
       */
	  }

}

/*----------------------------------------------------------------------- */
void buildNodeE(struct EltConst* ECp, double *node, int *GSize,
	   int *NodePos,  int jElt, int Nmnode) {

  int i1, j1, j2, Nnode, NInfoAtNode;

  Nnode=GSize[0];
  NInfoAtNode=GSize[4];
   for (j2=0;j2<Nnode;j2++) {/* nodeEe nodes of current element, x,y,z */
         i1=NodePos[j2+Nnode*jElt]-1;
         ECp[0].nodeE[j2]=node[i1+Nmnode*4]; 
         ECp[0].nodeE[j2+Nnode]=node[i1+Nmnode*5]; 
         ECp[0].nodeE[j2+2*Nnode]=node[i1+Nmnode*6]; 
   }
   if (ECp[0].giNodePos!=NULL) { /* rows correspond to fields, NodePos gives cols */
    for (j2=0;j2<Nnode;j2++) { 
     double *ccoef;
     ccoef=ECp[0].InfoAtNode+NInfoAtNode*(ECp[0].giNodePos[j2+Nnode*jElt]-1);
     for (j1=0;j1<NInfoAtNode;j1++) { 
       ECp[0].nodeE[j2+(4+j1)*Nnode]=ccoef[0];ccoef++; 
     }
    }
   }
}

/*----------------------------------------------------------------------- */
OF_EXPORT void GetpEC(struct EltConst EC, int *value) {
    void* pointer[1];
    pointer[0]=&EC; if (pointer[0]=NULL) mexErrMsgTxt("incoherent EC pointer");
    memcpy(value,&pointer,sizeof(mwIndex));
    value[3]=1;
}
/*----------------------------------------------------------------------- */
OF_EXPORT void SetpEC(struct EltConst *EC, int *value) {
   memcpy(EC,value,sizeof(struct EltConst));
}

/* #FieldPointer -----------------------------------------------------------*/
void* FieldPointer(const mxArray *field, char* name) {
    void* pointer;
    
    field=mxGetField(field,0,name);
    if (field !=NULL) {pointer =  mxGetData(field);} else { pointer=NULL;}
    return(pointer); 
}
/* #CheckFieldTypes --------------------------------------------------------*/
struct EltConst CheckFieldTypes(struct EltConst EC, mxArray *field, 
   int i0, int NfieldE) {

int j1,j2,j3;

EC.v1x=0; EC.v3x=0;EC.t=0;  EC.t=0;EC.h=0;/* EC.v3x=0; EC.h=0; EC.T=0; */
if (EC.NodeET==NULL) return(EC);
for (j1=0;j1<NfieldE-i0;j1++) { 
    mxArray *aSt; 
    mxChar *pSt;
    if (field!=NULL) {      
     aSt=mxGetCell(field,j1); pSt=(mxChar*)mxGetData(aSt); EC.NodeET[j1+i0]=0;j3=1;
     for (j2=0;j2<mxGetN(aSt);j2++) {
         EC.NodeET[j1+i0]+=(int)pSt[j2]*j3;j3*=256;
     }
    }
    if (EC.NodeET[j1+i0]==7876982) {EC.v1x=j1+i0;} /*comstr('v1x',-32)*/
    if (EC.NodeET[j1+i0]==7877494) {EC.v3x=j1+i0;} /*comstr('v3x',-32)*/
    if ((EC.NodeET[j1+i0]==84)||(EC.NodeET[j1+i0]==116)) {EC.t=j1+i0;} /*comstr('T',-32)*/
    if ((EC.NodeET[j1+i0]==72)||(EC.NodeET[j1+i0]==114)) {EC.h=j1+i0;} /*comstr('h',-32)*/
 }
 return(EC); 
}
/*----------------------------------------------------------------------- */
struct GroupFields initGroup(struct GroupFields GF, const mxArray *mNode, const mxArray *mEC, int *GSize, double* *NeedFree) {

  mxArray *field;
  int Nw=0;
  GF.RHS=NULL;GF.def=NULL;GF.defi=NULL;  GF.ke=NULL;
  GF.NBe=0;

  if (!mxIsStruct(mEC)) mexErrMsgTxt("EltConst must be a structure");
  if (mxIsEmpty(mEC)) mexErrMsgTxt("EltConst must not be empty");

  GF.N=mxGetPr(mxGetField(mEC, 0,"N")); 
  field = mxGetField(mEC, 0,"Nr");
  GF.Nshape=(int)mxGetN(field); GF.Nr=mxGetPr(field); GF.NdefE=0;

  field = mxGetField(mEC, 0,"Ns");
  if (field!=NULL) GF.Ns=mxGetPr(field); else GF.Ns=NULL;
  field = mxGetField(mEC, 0,"Nt");
  if (field!=NULL) GF.Nt=mxGetPr(field); else GF.Nt=NULL;

  field = mxGetField(mEC, 0,"NDN");
  if (field!=NULL) {
       GSize[1]=(int)mxGetM(field); /* NdofPerField */
       field = mxGetField(mEC, 0,"w");
       Nw=(int)mxGetM(field); GF.w = mxGetPr(field)+3*Nw; GSize[2]=Nw;
     }
  if (Nw==0) {Nw=(int)mxGetN(mxGetField(mEC, 0,"N"));}
  GSize[2]=Nw; 

  field = mxGetField(mEC, 0,"Ncondense");
  if (field!=NULL) {
       GF.Ncondense=(int)*mxGetPr (field);
   } else {GF.Ncondense=0;}

  field=mxGetField(mEC, 0,"VectMap");
  if (field==NULL) {
    GF.VectMap=NULL;
  } else {
   if (!mxIsInt32(field)) mexErrMsgTxt("VectMap must be int32");
   GF.VectMap=(int*)mxGetData(field);
  }
  field=mxGetField(mEC,0,"ConstitTopology");
  if (field!=NULL) { /* currently only first matrix topology */
   field = mxGetCell(field,0);
   if (field==NULL) {GF.topo=NULL;
   } else {
    if (mxIsInt32(field)) { GF.topo=(int*)mxGetData(field);
    }  else if (!mxIsChar(field)) mexErrMsgTxt("ConstitTopology must be int32");
   }
  } else {GF.topo=NULL; }
  GF.CTable=(double*)FieldPointer(mEC,(char*)"CTable");
  return(GF);
}

/*----------------------------------------------------------------------- */
struct EltConst  initInfoAtNode(const mxArray *infoAtNode,const mxArray *mEC,
              struct EltConst EC, struct GroupFields GF, int *GSize, double* *NeedFree) {

  int  NfieldE, Nw;
  mxArray     *field,*InfoAtNode_lab=NULL;

  /* NfieldE=4+GSize[4];*/
  Nw=GSize[2];

  EC.v1x=0; EC.ke=NULL;EC.Be=NULL; EC.defE=NULL; EC.NodeET=NULL;EC.gstate=NULL;

     if (infoAtNode==NULL) {
     } else if (mxIsStruct(infoAtNode)) {
       field=mxGetField(infoAtNode, 0,"data");
       if (field==NULL||mxIsSparse(field)) mexErrMsgTxt("IntoAtNode.data must be full");
       EC.InfoAtNode=mxGetPr(field);
       GSize[3]=(int)mxGetN(field);GSize[4]=(int)mxGetM(field);
       field=mxGetField(infoAtNode, 0,"NodePos");
       if (field==NULL||!mxIsInt32(field)) mexErrMsgTxt("IntoAtNode.NodePos must be int32");
       EC.giNodePos=(int*)mxGetData(field);
       if (EC.giNodePos[0]>GSize[3]) mexErrMsgTxt(".NodePos points to columns");
       field=mxGetField(infoAtNode, 0,"lab");InfoAtNode_lab=field;
       if (field==NULL||!mxIsCell(field)|| mxGetNumberOfElements(field)!=GSize[4]) {
           mexErrMsgTxt("IntoAtNode.lab must be a cell of length nField");
       }
    }   else { /* format with data at all nodes */
       EC.InfoAtNode=mxGetPr(infoAtNode);  EC.giNodePos=NULL; 
       GSize[3]=(int)mxGetM(infoAtNode); GSize[4]=(int)mxGetN(infoAtNode);field=NULL;
       if (GSize[3]>0) {
       mexWarnMsgTxt("ErrorTODO : matrix infoAtNode will disappear shortly");
       }
     } 
	 field=mxGetField(mEC, 0,"nodeEt");
     if (field==NULL) EC.NodeET=(int*)NeedFreeAlloc(NeedFree,5,GSize[4]+4);
	 else if (mxGetNumberOfElements(field)<GSize[4]+4)  mexErrMsgTxt("EC.nodeEt too small");
	 else EC.NodeET=(int*)mxGetData(field);
	 /* FieldPointer(mEC,"nodeEt");
      if (EC.NodeET==NULL) EC.NodeET=NeedFreeAlloc(NeedFree,5,GSize[4]+4);*/
     field=mxGetField(mEC, 0,"bas"); /* bas can be 4 in 2D and 9 in 3D */
     if (field==NULL||mxGetN(field)<Nw) EC.bas=NULL; 
     else  EC.bas = mxGetPr (field);
     field = mxGetField(mEC, 0,"NDN");
     if (field==NULL) EC.NDN=NULL; else EC.NDN= mxGetPr(field);
     field = mxGetField(mEC, 0,"jdet");
     if (field==NULL) EC.jdet=NULL; else {
       EC.jdet= mxGetPr(field); 
       if ( mxGetM(field)!=Nw ) mexErrMsgTxt("Problem of row allocation in jdet");
     }
     field=mxGetField(mEC, 0,"J"); 
     if (EC.jdet==NULL||field==NULL) {EC.J=NULL;
     } else if (mxGetM(field)==9&&mxGetN(field)==Nw) { EC.J = mxGetPr(field);
     } else if (mxGetM(field)!=4||mxGetN(field)!=Nw||mxGetField(mEC,0,"Nt")!=NULL) {EC.J=NULL;
     } else  EC.J = mxGetPr (field);

     field = mxGetField(mEC, 0,"gstate");
	 if (field!=NULL) { EC.gstate=mxGetPr(field);};

     field = mxGetField(mEC, 0,"nodeE");
     if (field!=NULL) {
      if (GSize[0]==0) { GSize[0]=(int)mxGetM(field);
	   } else if (mxGetM(field)!=GSize[0]) {
         mexErrMsgTxt("EltConst.nodeE and NodePos are inconsistent");
       }
	   if (EC.nodeE!=NULL && mxGetPr(field)!=EC.nodeE) mexErrMsgTxt("EltConst.nodeE redefines provided nodeE");
       NfieldE=(int)mxGetN(field); if (GSize[4]==0) GSize[4]=NfieldE-4;
       EC.nodeE=mxGetPr(field); NeedFree[0]=NULL; 
     } 

     EC=CheckFieldTypes(EC,InfoAtNode_lab,4,GSize[4]+4); 

     /* allow storage of fields at node to be passed to elements in nodeE*/
     if (field!=NULL) {
      if (NfieldE<(4+GSize[4]))  mexErrMsgTxt("nodeE must have at least 4+NInfoAtNode columns");
     } else {
      NfieldE=4+GSize[4];
	  /* EC.nodeE can be external nodeE provided as argument for BuildNDN */
      if (EC.nodeE==NULL) EC.nodeE=(double*)NeedFreeAlloc(NeedFree,0,GSize[0]*NfieldE);
     }
	 if (GF.CTable!=NULL) {
      EC.constit=(double*)FieldPointer(mEC,(char*)"constit");
	  EC.constit0=(double*)FieldPointer(mEC,(char*)"constit0");
     }
     return(EC);

}

/* #Transform33 : T=cGL --------------------------------*/
void Transform33(double *IN, double *T, double *S) {

	int ji,jj,jk,jl;

	for (ji=0;ji<9;ji++) S[ji]=0;
	for (ji=0;ji<3;ji++) {
	for (jj=0;jj<3;jj++) {
	for (jk=0;jk<3;jk++) {
	for (jl=0;jl<3;jl++) {
	 S[ji+3*jl]+=T[ji+3*jj]*IN[jj+3*jk]*T[jl+3*jk];
	}
	}
	}
	}

}
/* #Transform33s transform RHS to global basis but out in local -*/
void Transform33s(double *IN, double *T, double *S) {
	int ji,jj,jk;
	for (ji=0;ji<9;ji++) S[ji]=0;
	for (ji=0;ji<3;ji++) {
	for (jj=0;jj<3;jj++) {
	for (jk=0;jk<3;jk++) {
	 S[ji+3*jk]+=IN[ji+3*jj]*T[jk+3*jj];
	}
	}
	}
}

/* #TransformLambda stress tensor transform --------------------------------*/
void TransformLambda(double *IN, double *T, double *S) {


S[0]=((T[0]*T[0])*(T[0]*T[0])*IN[0]+(T[0]*T[0])*(T[3]*T[3])*IN[6]+(T[0]*T[0])*(T[6]*T[6])*IN[12]+(T[0]*T[0])*(T[3]*T[6]+T[6]*T[3])*IN[18]+(T[0]*T[0])*(T[6]*T[0]+T[0]*T[6])*IN[24]+(T[0]*T[0])*(T[0]*T[3]+T[3]*T[0])*IN[30]+(T[3]*T[3])*(T[0]*T[0])*IN[1]+(T[3]*T[3])*(T[3]*T[3])*IN[7]+(T[3]*T[3])*(T[6]*T[6])*IN[13]+(T[3]*T[3])*(T[3]*T[6]+T[6]*T[3])*IN[19]+(T[3]*T[3])*(T[6]*T[0]+T[0]*T[6])*IN[25]+(T[3]*T[3])*(T[0]*T[3]+T[3]*T[0])*IN[31]+(T[6]*T[6])*(T[0]*T[0])*IN[2]+(T[6]*T[6])*(T[3]*T[3])*IN[8]+(T[6]*T[6])*(T[6]*T[6])*IN[14]+(T[6]*T[6])*(T[3]*T[6]+T[6]*T[3])*IN[20]+(T[6]*T[6])*(T[6]*T[0]+T[0]*T[6])*IN[26]+(T[6]*T[6])*(T[0]*T[3]+T[3]*T[0])*IN[32]+(T[3]*T[6]+T[6]*T[3])*(T[0]*T[0])*IN[3]+(T[3]*T[6]+T[6]*T[3])*(T[3]*T[3])*IN[9]+(T[3]*T[6]+T[6]*T[3])*(T[6]*T[6])*IN[15]+(T[3]*T[6]+T[6]*T[3])*(T[3]*T[6]+T[6]*T[3])*IN[21]+(T[3]*T[6]+T[6]*T[3])*(T[6]*T[0]+T[0]*T[6])*IN[27]+(T[3]*T[6]+T[6]*T[3])*(T[0]*T[3]+T[3]*T[0])*IN[33]+(T[6]*T[0]+T[0]*T[6])*(T[0]*T[0])*IN[4]+(T[6]*T[0]+T[0]*T[6])*(T[3]*T[3])*IN[10]+(T[6]*T[0]+T[0]*T[6])*(T[6]*T[6])*IN[16]+(T[6]*T[0]+T[0]*T[6])*(T[3]*T[6]+T[6]*T[3])*IN[22]+(T[6]*T[0]+T[0]*T[6])*(T[6]*T[0]+T[0]*T[6])*IN[28]+(T[6]*T[0]+T[0]*T[6])*(T[0]*T[3]+T[3]*T[0])*IN[34]+(T[0]*T[3]+T[3]*T[0])*(T[0]*T[0])*IN[5]+(T[0]*T[3]+T[3]*T[0])*(T[3]*T[3])*IN[11]+(T[0]*T[3]+T[3]*T[0])*(T[6]*T[6])*IN[17]+(T[0]*T[3]+T[3]*T[0])*(T[3]*T[6]+T[6]*T[3])*IN[23]+(T[0]*T[3]+T[3]*T[0])*(T[6]*T[0]+T[0]*T[6])*IN[29]+(T[0]*T[3]+T[3]*T[0])*(T[0]*T[3]+T[3]*T[0])*IN[35]);
S[6]=((T[0]*T[0])*(T[1]*T[1])*IN[0]+(T[0]*T[0])*(T[4]*T[4])*IN[6]+(T[0]*T[0])*(T[7]*T[7])*IN[12]+(T[0]*T[0])*(T[4]*T[7]+T[7]*T[4])*IN[18]+(T[0]*T[0])*(T[7]*T[1]+T[1]*T[7])*IN[24]+(T[0]*T[0])*(T[1]*T[4]+T[4]*T[1])*IN[30]+(T[3]*T[3])*(T[1]*T[1])*IN[1]+(T[3]*T[3])*(T[4]*T[4])*IN[7]+(T[3]*T[3])*(T[7]*T[7])*IN[13]+(T[3]*T[3])*(T[4]*T[7]+T[7]*T[4])*IN[19]+(T[3]*T[3])*(T[7]*T[1]+T[1]*T[7])*IN[25]+(T[3]*T[3])*(T[1]*T[4]+T[4]*T[1])*IN[31]+(T[6]*T[6])*(T[1]*T[1])*IN[2]+(T[6]*T[6])*(T[4]*T[4])*IN[8]+(T[6]*T[6])*(T[7]*T[7])*IN[14]+(T[6]*T[6])*(T[4]*T[7]+T[7]*T[4])*IN[20]+(T[6]*T[6])*(T[7]*T[1]+T[1]*T[7])*IN[26]+(T[6]*T[6])*(T[1]*T[4]+T[4]*T[1])*IN[32]+(T[3]*T[6]+T[6]*T[3])*(T[1]*T[1])*IN[3]+(T[3]*T[6]+T[6]*T[3])*(T[4]*T[4])*IN[9]+(T[3]*T[6]+T[6]*T[3])*(T[7]*T[7])*IN[15]+(T[3]*T[6]+T[6]*T[3])*(T[4]*T[7]+T[7]*T[4])*IN[21]+(T[3]*T[6]+T[6]*T[3])*(T[7]*T[1]+T[1]*T[7])*IN[27]+(T[3]*T[6]+T[6]*T[3])*(T[1]*T[4]+T[4]*T[1])*IN[33]+(T[6]*T[0]+T[0]*T[6])*(T[1]*T[1])*IN[4]+(T[6]*T[0]+T[0]*T[6])*(T[4]*T[4])*IN[10]+(T[6]*T[0]+T[0]*T[6])*(T[7]*T[7])*IN[16]+(T[6]*T[0]+T[0]*T[6])*(T[4]*T[7]+T[7]*T[4])*IN[22]+(T[6]*T[0]+T[0]*T[6])*(T[7]*T[1]+T[1]*T[7])*IN[28]+(T[6]*T[0]+T[0]*T[6])*(T[1]*T[4]+T[4]*T[1])*IN[34]+(T[0]*T[3]+T[3]*T[0])*(T[1]*T[1])*IN[5]+(T[0]*T[3]+T[3]*T[0])*(T[4]*T[4])*IN[11]+(T[0]*T[3]+T[3]*T[0])*(T[7]*T[7])*IN[17]+(T[0]*T[3]+T[3]*T[0])*(T[4]*T[7]+T[7]*T[4])*IN[23]+(T[0]*T[3]+T[3]*T[0])*(T[7]*T[1]+T[1]*T[7])*IN[29]+(T[0]*T[3]+T[3]*T[0])*(T[1]*T[4]+T[4]*T[1])*IN[35]);
S[12]=((T[0]*T[0])*(T[2]*T[2])*IN[0]+(T[0]*T[0])*(T[5]*T[5])*IN[6]+(T[0]*T[0])*(T[8]*T[8])*IN[12]+(T[0]*T[0])*(T[5]*T[8]+T[8]*T[5])*IN[18]+(T[0]*T[0])*(T[8]*T[2]+T[2]*T[8])*IN[24]+(T[0]*T[0])*(T[2]*T[5]+T[5]*T[2])*IN[30]+(T[3]*T[3])*(T[2]*T[2])*IN[1]+(T[3]*T[3])*(T[5]*T[5])*IN[7]+(T[3]*T[3])*(T[8]*T[8])*IN[13]+(T[3]*T[3])*(T[5]*T[8]+T[8]*T[5])*IN[19]+(T[3]*T[3])*(T[8]*T[2]+T[2]*T[8])*IN[25]+(T[3]*T[3])*(T[2]*T[5]+T[5]*T[2])*IN[31]+(T[6]*T[6])*(T[2]*T[2])*IN[2]+(T[6]*T[6])*(T[5]*T[5])*IN[8]+(T[6]*T[6])*(T[8]*T[8])*IN[14]+(T[6]*T[6])*(T[5]*T[8]+T[8]*T[5])*IN[20]+(T[6]*T[6])*(T[8]*T[2]+T[2]*T[8])*IN[26]+(T[6]*T[6])*(T[2]*T[5]+T[5]*T[2])*IN[32]+(T[3]*T[6]+T[6]*T[3])*(T[2]*T[2])*IN[3]+(T[3]*T[6]+T[6]*T[3])*(T[5]*T[5])*IN[9]+(T[3]*T[6]+T[6]*T[3])*(T[8]*T[8])*IN[15]+(T[3]*T[6]+T[6]*T[3])*(T[5]*T[8]+T[8]*T[5])*IN[21]+(T[3]*T[6]+T[6]*T[3])*(T[8]*T[2]+T[2]*T[8])*IN[27]+(T[3]*T[6]+T[6]*T[3])*(T[2]*T[5]+T[5]*T[2])*IN[33]+(T[6]*T[0]+T[0]*T[6])*(T[2]*T[2])*IN[4]+(T[6]*T[0]+T[0]*T[6])*(T[5]*T[5])*IN[10]+(T[6]*T[0]+T[0]*T[6])*(T[8]*T[8])*IN[16]+(T[6]*T[0]+T[0]*T[6])*(T[5]*T[8]+T[8]*T[5])*IN[22]+(T[6]*T[0]+T[0]*T[6])*(T[8]*T[2]+T[2]*T[8])*IN[28]+(T[6]*T[0]+T[0]*T[6])*(T[2]*T[5]+T[5]*T[2])*IN[34]+(T[0]*T[3]+T[3]*T[0])*(T[2]*T[2])*IN[5]+(T[0]*T[3]+T[3]*T[0])*(T[5]*T[5])*IN[11]+(T[0]*T[3]+T[3]*T[0])*(T[8]*T[8])*IN[17]+(T[0]*T[3]+T[3]*T[0])*(T[5]*T[8]+T[8]*T[5])*IN[23]+(T[0]*T[3]+T[3]*T[0])*(T[8]*T[2]+T[2]*T[8])*IN[29]+(T[0]*T[3]+T[3]*T[0])*(T[2]*T[5]+T[5]*T[2])*IN[35]);
S[18]=((T[0]*T[0])*(T[1]*T[2])*IN[0]+(T[0]*T[0])*(T[4]*T[5])*IN[6]+(T[0]*T[0])*(T[7]*T[8])*IN[12]+(T[0]*T[0])*(T[4]*T[8]+T[7]*T[5])*IN[18]+(T[0]*T[0])*(T[7]*T[2]+T[1]*T[8])*IN[24]+(T[0]*T[0])*(T[1]*T[5]+T[4]*T[2])*IN[30]+(T[3]*T[3])*(T[1]*T[2])*IN[1]+(T[3]*T[3])*(T[4]*T[5])*IN[7]+(T[3]*T[3])*(T[7]*T[8])*IN[13]+(T[3]*T[3])*(T[4]*T[8]+T[7]*T[5])*IN[19]+(T[3]*T[3])*(T[7]*T[2]+T[1]*T[8])*IN[25]+(T[3]*T[3])*(T[1]*T[5]+T[4]*T[2])*IN[31]+(T[6]*T[6])*(T[1]*T[2])*IN[2]+(T[6]*T[6])*(T[4]*T[5])*IN[8]+(T[6]*T[6])*(T[7]*T[8])*IN[14]+(T[6]*T[6])*(T[4]*T[8]+T[7]*T[5])*IN[20]+(T[6]*T[6])*(T[7]*T[2]+T[1]*T[8])*IN[26]+(T[6]*T[6])*(T[1]*T[5]+T[4]*T[2])*IN[32]+(T[3]*T[6]+T[6]*T[3])*(T[1]*T[2])*IN[3]+(T[3]*T[6]+T[6]*T[3])*(T[4]*T[5])*IN[9]+(T[3]*T[6]+T[6]*T[3])*(T[7]*T[8])*IN[15]+(T[3]*T[6]+T[6]*T[3])*(T[4]*T[8]+T[7]*T[5])*IN[21]+(T[3]*T[6]+T[6]*T[3])*(T[7]*T[2]+T[1]*T[8])*IN[27]+(T[3]*T[6]+T[6]*T[3])*(T[1]*T[5]+T[4]*T[2])*IN[33]+(T[6]*T[0]+T[0]*T[6])*(T[1]*T[2])*IN[4]+(T[6]*T[0]+T[0]*T[6])*(T[4]*T[5])*IN[10]+(T[6]*T[0]+T[0]*T[6])*(T[7]*T[8])*IN[16]+(T[6]*T[0]+T[0]*T[6])*(T[4]*T[8]+T[7]*T[5])*IN[22]+(T[6]*T[0]+T[0]*T[6])*(T[7]*T[2]+T[1]*T[8])*IN[28]+(T[6]*T[0]+T[0]*T[6])*(T[1]*T[5]+T[4]*T[2])*IN[34]+(T[0]*T[3]+T[3]*T[0])*(T[1]*T[2])*IN[5]+(T[0]*T[3]+T[3]*T[0])*(T[4]*T[5])*IN[11]+(T[0]*T[3]+T[3]*T[0])*(T[7]*T[8])*IN[17]+(T[0]*T[3]+T[3]*T[0])*(T[4]*T[8]+T[7]*T[5])*IN[23]+(T[0]*T[3]+T[3]*T[0])*(T[7]*T[2]+T[1]*T[8])*IN[29]+(T[0]*T[3]+T[3]*T[0])*(T[1]*T[5]+T[4]*T[2])*IN[35]);
S[24]=((T[0]*T[0])*(T[2]*T[0])*IN[0]+(T[0]*T[0])*(T[5]*T[3])*IN[6]+(T[0]*T[0])*(T[8]*T[6])*IN[12]+(T[0]*T[0])*(T[5]*T[6]+T[8]*T[3])*IN[18]+(T[0]*T[0])*(T[8]*T[0]+T[2]*T[6])*IN[24]+(T[0]*T[0])*(T[2]*T[3]+T[5]*T[0])*IN[30]+(T[3]*T[3])*(T[2]*T[0])*IN[1]+(T[3]*T[3])*(T[5]*T[3])*IN[7]+(T[3]*T[3])*(T[8]*T[6])*IN[13]+(T[3]*T[3])*(T[5]*T[6]+T[8]*T[3])*IN[19]+(T[3]*T[3])*(T[8]*T[0]+T[2]*T[6])*IN[25]+(T[3]*T[3])*(T[2]*T[3]+T[5]*T[0])*IN[31]+(T[6]*T[6])*(T[2]*T[0])*IN[2]+(T[6]*T[6])*(T[5]*T[3])*IN[8]+(T[6]*T[6])*(T[8]*T[6])*IN[14]+(T[6]*T[6])*(T[5]*T[6]+T[8]*T[3])*IN[20]+(T[6]*T[6])*(T[8]*T[0]+T[2]*T[6])*IN[26]+(T[6]*T[6])*(T[2]*T[3]+T[5]*T[0])*IN[32]+(T[3]*T[6]+T[6]*T[3])*(T[2]*T[0])*IN[3]+(T[3]*T[6]+T[6]*T[3])*(T[5]*T[3])*IN[9]+(T[3]*T[6]+T[6]*T[3])*(T[8]*T[6])*IN[15]+(T[3]*T[6]+T[6]*T[3])*(T[5]*T[6]+T[8]*T[3])*IN[21]+(T[3]*T[6]+T[6]*T[3])*(T[8]*T[0]+T[2]*T[6])*IN[27]+(T[3]*T[6]+T[6]*T[3])*(T[2]*T[3]+T[5]*T[0])*IN[33]+(T[6]*T[0]+T[0]*T[6])*(T[2]*T[0])*IN[4]+(T[6]*T[0]+T[0]*T[6])*(T[5]*T[3])*IN[10]+(T[6]*T[0]+T[0]*T[6])*(T[8]*T[6])*IN[16]+(T[6]*T[0]+T[0]*T[6])*(T[5]*T[6]+T[8]*T[3])*IN[22]+(T[6]*T[0]+T[0]*T[6])*(T[8]*T[0]+T[2]*T[6])*IN[28]+(T[6]*T[0]+T[0]*T[6])*(T[2]*T[3]+T[5]*T[0])*IN[34]+(T[0]*T[3]+T[3]*T[0])*(T[2]*T[0])*IN[5]+(T[0]*T[3]+T[3]*T[0])*(T[5]*T[3])*IN[11]+(T[0]*T[3]+T[3]*T[0])*(T[8]*T[6])*IN[17]+(T[0]*T[3]+T[3]*T[0])*(T[5]*T[6]+T[8]*T[3])*IN[23]+(T[0]*T[3]+T[3]*T[0])*(T[8]*T[0]+T[2]*T[6])*IN[29]+(T[0]*T[3]+T[3]*T[0])*(T[2]*T[3]+T[5]*T[0])*IN[35]);
S[30]=((T[0]*T[0])*(T[0]*T[1])*IN[0]+(T[0]*T[0])*(T[3]*T[4])*IN[6]+(T[0]*T[0])*(T[6]*T[7])*IN[12]+(T[0]*T[0])*(T[3]*T[7]+T[6]*T[4])*IN[18]+(T[0]*T[0])*(T[6]*T[1]+T[0]*T[7])*IN[24]+(T[0]*T[0])*(T[0]*T[4]+T[3]*T[1])*IN[30]+(T[3]*T[3])*(T[0]*T[1])*IN[1]+(T[3]*T[3])*(T[3]*T[4])*IN[7]+(T[3]*T[3])*(T[6]*T[7])*IN[13]+(T[3]*T[3])*(T[3]*T[7]+T[6]*T[4])*IN[19]+(T[3]*T[3])*(T[6]*T[1]+T[0]*T[7])*IN[25]+(T[3]*T[3])*(T[0]*T[4]+T[3]*T[1])*IN[31]+(T[6]*T[6])*(T[0]*T[1])*IN[2]+(T[6]*T[6])*(T[3]*T[4])*IN[8]+(T[6]*T[6])*(T[6]*T[7])*IN[14]+(T[6]*T[6])*(T[3]*T[7]+T[6]*T[4])*IN[20]+(T[6]*T[6])*(T[6]*T[1]+T[0]*T[7])*IN[26]+(T[6]*T[6])*(T[0]*T[4]+T[3]*T[1])*IN[32]+(T[3]*T[6]+T[6]*T[3])*(T[0]*T[1])*IN[3]+(T[3]*T[6]+T[6]*T[3])*(T[3]*T[4])*IN[9]+(T[3]*T[6]+T[6]*T[3])*(T[6]*T[7])*IN[15]+(T[3]*T[6]+T[6]*T[3])*(T[3]*T[7]+T[6]*T[4])*IN[21]+(T[3]*T[6]+T[6]*T[3])*(T[6]*T[1]+T[0]*T[7])*IN[27]+(T[3]*T[6]+T[6]*T[3])*(T[0]*T[4]+T[3]*T[1])*IN[33]+(T[6]*T[0]+T[0]*T[6])*(T[0]*T[1])*IN[4]+(T[6]*T[0]+T[0]*T[6])*(T[3]*T[4])*IN[10]+(T[6]*T[0]+T[0]*T[6])*(T[6]*T[7])*IN[16]+(T[6]*T[0]+T[0]*T[6])*(T[3]*T[7]+T[6]*T[4])*IN[22]+(T[6]*T[0]+T[0]*T[6])*(T[6]*T[1]+T[0]*T[7])*IN[28]+(T[6]*T[0]+T[0]*T[6])*(T[0]*T[4]+T[3]*T[1])*IN[34]+(T[0]*T[3]+T[3]*T[0])*(T[0]*T[1])*IN[5]+(T[0]*T[3]+T[3]*T[0])*(T[3]*T[4])*IN[11]+(T[0]*T[3]+T[3]*T[0])*(T[6]*T[7])*IN[17]+(T[0]*T[3]+T[3]*T[0])*(T[3]*T[7]+T[6]*T[4])*IN[23]+(T[0]*T[3]+T[3]*T[0])*(T[6]*T[1]+T[0]*T[7])*IN[29]+(T[0]*T[3]+T[3]*T[0])*(T[0]*T[4]+T[3]*T[1])*IN[35]);
S[1]=((T[1]*T[1])*(T[0]*T[0])*IN[0]+(T[1]*T[1])*(T[3]*T[3])*IN[6]+(T[1]*T[1])*(T[6]*T[6])*IN[12]+(T[1]*T[1])*(T[3]*T[6]+T[6]*T[3])*IN[18]+(T[1]*T[1])*(T[6]*T[0]+T[0]*T[6])*IN[24]+(T[1]*T[1])*(T[0]*T[3]+T[3]*T[0])*IN[30]+(T[4]*T[4])*(T[0]*T[0])*IN[1]+(T[4]*T[4])*(T[3]*T[3])*IN[7]+(T[4]*T[4])*(T[6]*T[6])*IN[13]+(T[4]*T[4])*(T[3]*T[6]+T[6]*T[3])*IN[19]+(T[4]*T[4])*(T[6]*T[0]+T[0]*T[6])*IN[25]+(T[4]*T[4])*(T[0]*T[3]+T[3]*T[0])*IN[31]+(T[7]*T[7])*(T[0]*T[0])*IN[2]+(T[7]*T[7])*(T[3]*T[3])*IN[8]+(T[7]*T[7])*(T[6]*T[6])*IN[14]+(T[7]*T[7])*(T[3]*T[6]+T[6]*T[3])*IN[20]+(T[7]*T[7])*(T[6]*T[0]+T[0]*T[6])*IN[26]+(T[7]*T[7])*(T[0]*T[3]+T[3]*T[0])*IN[32]+(T[4]*T[7]+T[7]*T[4])*(T[0]*T[0])*IN[3]+(T[4]*T[7]+T[7]*T[4])*(T[3]*T[3])*IN[9]+(T[4]*T[7]+T[7]*T[4])*(T[6]*T[6])*IN[15]+(T[4]*T[7]+T[7]*T[4])*(T[3]*T[6]+T[6]*T[3])*IN[21]+(T[4]*T[7]+T[7]*T[4])*(T[6]*T[0]+T[0]*T[6])*IN[27]+(T[4]*T[7]+T[7]*T[4])*(T[0]*T[3]+T[3]*T[0])*IN[33]+(T[7]*T[1]+T[1]*T[7])*(T[0]*T[0])*IN[4]+(T[7]*T[1]+T[1]*T[7])*(T[3]*T[3])*IN[10]+(T[7]*T[1]+T[1]*T[7])*(T[6]*T[6])*IN[16]+(T[7]*T[1]+T[1]*T[7])*(T[3]*T[6]+T[6]*T[3])*IN[22]+(T[7]*T[1]+T[1]*T[7])*(T[6]*T[0]+T[0]*T[6])*IN[28]+(T[7]*T[1]+T[1]*T[7])*(T[0]*T[3]+T[3]*T[0])*IN[34]+(T[1]*T[4]+T[4]*T[1])*(T[0]*T[0])*IN[5]+(T[1]*T[4]+T[4]*T[1])*(T[3]*T[3])*IN[11]+(T[1]*T[4]+T[4]*T[1])*(T[6]*T[6])*IN[17]+(T[1]*T[4]+T[4]*T[1])*(T[3]*T[6]+T[6]*T[3])*IN[23]+(T[1]*T[4]+T[4]*T[1])*(T[6]*T[0]+T[0]*T[6])*IN[29]+(T[1]*T[4]+T[4]*T[1])*(T[0]*T[3]+T[3]*T[0])*IN[35]);
S[7]=((T[1]*T[1])*(T[1]*T[1])*IN[0]+(T[1]*T[1])*(T[4]*T[4])*IN[6]+(T[1]*T[1])*(T[7]*T[7])*IN[12]+(T[1]*T[1])*(T[4]*T[7]+T[7]*T[4])*IN[18]+(T[1]*T[1])*(T[7]*T[1]+T[1]*T[7])*IN[24]+(T[1]*T[1])*(T[1]*T[4]+T[4]*T[1])*IN[30]+(T[4]*T[4])*(T[1]*T[1])*IN[1]+(T[4]*T[4])*(T[4]*T[4])*IN[7]+(T[4]*T[4])*(T[7]*T[7])*IN[13]+(T[4]*T[4])*(T[4]*T[7]+T[7]*T[4])*IN[19]+(T[4]*T[4])*(T[7]*T[1]+T[1]*T[7])*IN[25]+(T[4]*T[4])*(T[1]*T[4]+T[4]*T[1])*IN[31]+(T[7]*T[7])*(T[1]*T[1])*IN[2]+(T[7]*T[7])*(T[4]*T[4])*IN[8]+(T[7]*T[7])*(T[7]*T[7])*IN[14]+(T[7]*T[7])*(T[4]*T[7]+T[7]*T[4])*IN[20]+(T[7]*T[7])*(T[7]*T[1]+T[1]*T[7])*IN[26]+(T[7]*T[7])*(T[1]*T[4]+T[4]*T[1])*IN[32]+(T[4]*T[7]+T[7]*T[4])*(T[1]*T[1])*IN[3]+(T[4]*T[7]+T[7]*T[4])*(T[4]*T[4])*IN[9]+(T[4]*T[7]+T[7]*T[4])*(T[7]*T[7])*IN[15]+(T[4]*T[7]+T[7]*T[4])*(T[4]*T[7]+T[7]*T[4])*IN[21]+(T[4]*T[7]+T[7]*T[4])*(T[7]*T[1]+T[1]*T[7])*IN[27]+(T[4]*T[7]+T[7]*T[4])*(T[1]*T[4]+T[4]*T[1])*IN[33]+(T[7]*T[1]+T[1]*T[7])*(T[1]*T[1])*IN[4]+(T[7]*T[1]+T[1]*T[7])*(T[4]*T[4])*IN[10]+(T[7]*T[1]+T[1]*T[7])*(T[7]*T[7])*IN[16]+(T[7]*T[1]+T[1]*T[7])*(T[4]*T[7]+T[7]*T[4])*IN[22]+(T[7]*T[1]+T[1]*T[7])*(T[7]*T[1]+T[1]*T[7])*IN[28]+(T[7]*T[1]+T[1]*T[7])*(T[1]*T[4]+T[4]*T[1])*IN[34]+(T[1]*T[4]+T[4]*T[1])*(T[1]*T[1])*IN[5]+(T[1]*T[4]+T[4]*T[1])*(T[4]*T[4])*IN[11]+(T[1]*T[4]+T[4]*T[1])*(T[7]*T[7])*IN[17]+(T[1]*T[4]+T[4]*T[1])*(T[4]*T[7]+T[7]*T[4])*IN[23]+(T[1]*T[4]+T[4]*T[1])*(T[7]*T[1]+T[1]*T[7])*IN[29]+(T[1]*T[4]+T[4]*T[1])*(T[1]*T[4]+T[4]*T[1])*IN[35]);
S[13]=((T[1]*T[1])*(T[2]*T[2])*IN[0]+(T[1]*T[1])*(T[5]*T[5])*IN[6]+(T[1]*T[1])*(T[8]*T[8])*IN[12]+(T[1]*T[1])*(T[5]*T[8]+T[8]*T[5])*IN[18]+(T[1]*T[1])*(T[8]*T[2]+T[2]*T[8])*IN[24]+(T[1]*T[1])*(T[2]*T[5]+T[5]*T[2])*IN[30]+(T[4]*T[4])*(T[2]*T[2])*IN[1]+(T[4]*T[4])*(T[5]*T[5])*IN[7]+(T[4]*T[4])*(T[8]*T[8])*IN[13]+(T[4]*T[4])*(T[5]*T[8]+T[8]*T[5])*IN[19]+(T[4]*T[4])*(T[8]*T[2]+T[2]*T[8])*IN[25]+(T[4]*T[4])*(T[2]*T[5]+T[5]*T[2])*IN[31]+(T[7]*T[7])*(T[2]*T[2])*IN[2]+(T[7]*T[7])*(T[5]*T[5])*IN[8]+(T[7]*T[7])*(T[8]*T[8])*IN[14]+(T[7]*T[7])*(T[5]*T[8]+T[8]*T[5])*IN[20]+(T[7]*T[7])*(T[8]*T[2]+T[2]*T[8])*IN[26]+(T[7]*T[7])*(T[2]*T[5]+T[5]*T[2])*IN[32]+(T[4]*T[7]+T[7]*T[4])*(T[2]*T[2])*IN[3]+(T[4]*T[7]+T[7]*T[4])*(T[5]*T[5])*IN[9]+(T[4]*T[7]+T[7]*T[4])*(T[8]*T[8])*IN[15]+(T[4]*T[7]+T[7]*T[4])*(T[5]*T[8]+T[8]*T[5])*IN[21]+(T[4]*T[7]+T[7]*T[4])*(T[8]*T[2]+T[2]*T[8])*IN[27]+(T[4]*T[7]+T[7]*T[4])*(T[2]*T[5]+T[5]*T[2])*IN[33]+(T[7]*T[1]+T[1]*T[7])*(T[2]*T[2])*IN[4]+(T[7]*T[1]+T[1]*T[7])*(T[5]*T[5])*IN[10]+(T[7]*T[1]+T[1]*T[7])*(T[8]*T[8])*IN[16]+(T[7]*T[1]+T[1]*T[7])*(T[5]*T[8]+T[8]*T[5])*IN[22]+(T[7]*T[1]+T[1]*T[7])*(T[8]*T[2]+T[2]*T[8])*IN[28]+(T[7]*T[1]+T[1]*T[7])*(T[2]*T[5]+T[5]*T[2])*IN[34]+(T[1]*T[4]+T[4]*T[1])*(T[2]*T[2])*IN[5]+(T[1]*T[4]+T[4]*T[1])*(T[5]*T[5])*IN[11]+(T[1]*T[4]+T[4]*T[1])*(T[8]*T[8])*IN[17]+(T[1]*T[4]+T[4]*T[1])*(T[5]*T[8]+T[8]*T[5])*IN[23]+(T[1]*T[4]+T[4]*T[1])*(T[8]*T[2]+T[2]*T[8])*IN[29]+(T[1]*T[4]+T[4]*T[1])*(T[2]*T[5]+T[5]*T[2])*IN[35]);
S[19]=((T[1]*T[1])*(T[1]*T[2])*IN[0]+(T[1]*T[1])*(T[4]*T[5])*IN[6]+(T[1]*T[1])*(T[7]*T[8])*IN[12]+(T[1]*T[1])*(T[4]*T[8]+T[7]*T[5])*IN[18]+(T[1]*T[1])*(T[7]*T[2]+T[1]*T[8])*IN[24]+(T[1]*T[1])*(T[1]*T[5]+T[4]*T[2])*IN[30]+(T[4]*T[4])*(T[1]*T[2])*IN[1]+(T[4]*T[4])*(T[4]*T[5])*IN[7]+(T[4]*T[4])*(T[7]*T[8])*IN[13]+(T[4]*T[4])*(T[4]*T[8]+T[7]*T[5])*IN[19]+(T[4]*T[4])*(T[7]*T[2]+T[1]*T[8])*IN[25]+(T[4]*T[4])*(T[1]*T[5]+T[4]*T[2])*IN[31]+(T[7]*T[7])*(T[1]*T[2])*IN[2]+(T[7]*T[7])*(T[4]*T[5])*IN[8]+(T[7]*T[7])*(T[7]*T[8])*IN[14]+(T[7]*T[7])*(T[4]*T[8]+T[7]*T[5])*IN[20]+(T[7]*T[7])*(T[7]*T[2]+T[1]*T[8])*IN[26]+(T[7]*T[7])*(T[1]*T[5]+T[4]*T[2])*IN[32]+(T[4]*T[7]+T[7]*T[4])*(T[1]*T[2])*IN[3]+(T[4]*T[7]+T[7]*T[4])*(T[4]*T[5])*IN[9]+(T[4]*T[7]+T[7]*T[4])*(T[7]*T[8])*IN[15]+(T[4]*T[7]+T[7]*T[4])*(T[4]*T[8]+T[7]*T[5])*IN[21]+(T[4]*T[7]+T[7]*T[4])*(T[7]*T[2]+T[1]*T[8])*IN[27]+(T[4]*T[7]+T[7]*T[4])*(T[1]*T[5]+T[4]*T[2])*IN[33]+(T[7]*T[1]+T[1]*T[7])*(T[1]*T[2])*IN[4]+(T[7]*T[1]+T[1]*T[7])*(T[4]*T[5])*IN[10]+(T[7]*T[1]+T[1]*T[7])*(T[7]*T[8])*IN[16]+(T[7]*T[1]+T[1]*T[7])*(T[4]*T[8]+T[7]*T[5])*IN[22]+(T[7]*T[1]+T[1]*T[7])*(T[7]*T[2]+T[1]*T[8])*IN[28]+(T[7]*T[1]+T[1]*T[7])*(T[1]*T[5]+T[4]*T[2])*IN[34]+(T[1]*T[4]+T[4]*T[1])*(T[1]*T[2])*IN[5]+(T[1]*T[4]+T[4]*T[1])*(T[4]*T[5])*IN[11]+(T[1]*T[4]+T[4]*T[1])*(T[7]*T[8])*IN[17]+(T[1]*T[4]+T[4]*T[1])*(T[4]*T[8]+T[7]*T[5])*IN[23]+(T[1]*T[4]+T[4]*T[1])*(T[7]*T[2]+T[1]*T[8])*IN[29]+(T[1]*T[4]+T[4]*T[1])*(T[1]*T[5]+T[4]*T[2])*IN[35]);
S[25]=((T[1]*T[1])*(T[2]*T[0])*IN[0]+(T[1]*T[1])*(T[5]*T[3])*IN[6]+(T[1]*T[1])*(T[8]*T[6])*IN[12]+(T[1]*T[1])*(T[5]*T[6]+T[8]*T[3])*IN[18]+(T[1]*T[1])*(T[8]*T[0]+T[2]*T[6])*IN[24]+(T[1]*T[1])*(T[2]*T[3]+T[5]*T[0])*IN[30]+(T[4]*T[4])*(T[2]*T[0])*IN[1]+(T[4]*T[4])*(T[5]*T[3])*IN[7]+(T[4]*T[4])*(T[8]*T[6])*IN[13]+(T[4]*T[4])*(T[5]*T[6]+T[8]*T[3])*IN[19]+(T[4]*T[4])*(T[8]*T[0]+T[2]*T[6])*IN[25]+(T[4]*T[4])*(T[2]*T[3]+T[5]*T[0])*IN[31]+(T[7]*T[7])*(T[2]*T[0])*IN[2]+(T[7]*T[7])*(T[5]*T[3])*IN[8]+(T[7]*T[7])*(T[8]*T[6])*IN[14]+(T[7]*T[7])*(T[5]*T[6]+T[8]*T[3])*IN[20]+(T[7]*T[7])*(T[8]*T[0]+T[2]*T[6])*IN[26]+(T[7]*T[7])*(T[2]*T[3]+T[5]*T[0])*IN[32]+(T[4]*T[7]+T[7]*T[4])*(T[2]*T[0])*IN[3]+(T[4]*T[7]+T[7]*T[4])*(T[5]*T[3])*IN[9]+(T[4]*T[7]+T[7]*T[4])*(T[8]*T[6])*IN[15]+(T[4]*T[7]+T[7]*T[4])*(T[5]*T[6]+T[8]*T[3])*IN[21]+(T[4]*T[7]+T[7]*T[4])*(T[8]*T[0]+T[2]*T[6])*IN[27]+(T[4]*T[7]+T[7]*T[4])*(T[2]*T[3]+T[5]*T[0])*IN[33]+(T[7]*T[1]+T[1]*T[7])*(T[2]*T[0])*IN[4]+(T[7]*T[1]+T[1]*T[7])*(T[5]*T[3])*IN[10]+(T[7]*T[1]+T[1]*T[7])*(T[8]*T[6])*IN[16]+(T[7]*T[1]+T[1]*T[7])*(T[5]*T[6]+T[8]*T[3])*IN[22]+(T[7]*T[1]+T[1]*T[7])*(T[8]*T[0]+T[2]*T[6])*IN[28]+(T[7]*T[1]+T[1]*T[7])*(T[2]*T[3]+T[5]*T[0])*IN[34]+(T[1]*T[4]+T[4]*T[1])*(T[2]*T[0])*IN[5]+(T[1]*T[4]+T[4]*T[1])*(T[5]*T[3])*IN[11]+(T[1]*T[4]+T[4]*T[1])*(T[8]*T[6])*IN[17]+(T[1]*T[4]+T[4]*T[1])*(T[5]*T[6]+T[8]*T[3])*IN[23]+(T[1]*T[4]+T[4]*T[1])*(T[8]*T[0]+T[2]*T[6])*IN[29]+(T[1]*T[4]+T[4]*T[1])*(T[2]*T[3]+T[5]*T[0])*IN[35]);
S[31]=((T[1]*T[1])*(T[0]*T[1])*IN[0]+(T[1]*T[1])*(T[3]*T[4])*IN[6]+(T[1]*T[1])*(T[6]*T[7])*IN[12]+(T[1]*T[1])*(T[3]*T[7]+T[6]*T[4])*IN[18]+(T[1]*T[1])*(T[6]*T[1]+T[0]*T[7])*IN[24]+(T[1]*T[1])*(T[0]*T[4]+T[3]*T[1])*IN[30]+(T[4]*T[4])*(T[0]*T[1])*IN[1]+(T[4]*T[4])*(T[3]*T[4])*IN[7]+(T[4]*T[4])*(T[6]*T[7])*IN[13]+(T[4]*T[4])*(T[3]*T[7]+T[6]*T[4])*IN[19]+(T[4]*T[4])*(T[6]*T[1]+T[0]*T[7])*IN[25]+(T[4]*T[4])*(T[0]*T[4]+T[3]*T[1])*IN[31]+(T[7]*T[7])*(T[0]*T[1])*IN[2]+(T[7]*T[7])*(T[3]*T[4])*IN[8]+(T[7]*T[7])*(T[6]*T[7])*IN[14]+(T[7]*T[7])*(T[3]*T[7]+T[6]*T[4])*IN[20]+(T[7]*T[7])*(T[6]*T[1]+T[0]*T[7])*IN[26]+(T[7]*T[7])*(T[0]*T[4]+T[3]*T[1])*IN[32]+(T[4]*T[7]+T[7]*T[4])*(T[0]*T[1])*IN[3]+(T[4]*T[7]+T[7]*T[4])*(T[3]*T[4])*IN[9]+(T[4]*T[7]+T[7]*T[4])*(T[6]*T[7])*IN[15]+(T[4]*T[7]+T[7]*T[4])*(T[3]*T[7]+T[6]*T[4])*IN[21]+(T[4]*T[7]+T[7]*T[4])*(T[6]*T[1]+T[0]*T[7])*IN[27]+(T[4]*T[7]+T[7]*T[4])*(T[0]*T[4]+T[3]*T[1])*IN[33]+(T[7]*T[1]+T[1]*T[7])*(T[0]*T[1])*IN[4]+(T[7]*T[1]+T[1]*T[7])*(T[3]*T[4])*IN[10]+(T[7]*T[1]+T[1]*T[7])*(T[6]*T[7])*IN[16]+(T[7]*T[1]+T[1]*T[7])*(T[3]*T[7]+T[6]*T[4])*IN[22]+(T[7]*T[1]+T[1]*T[7])*(T[6]*T[1]+T[0]*T[7])*IN[28]+(T[7]*T[1]+T[1]*T[7])*(T[0]*T[4]+T[3]*T[1])*IN[34]+(T[1]*T[4]+T[4]*T[1])*(T[0]*T[1])*IN[5]+(T[1]*T[4]+T[4]*T[1])*(T[3]*T[4])*IN[11]+(T[1]*T[4]+T[4]*T[1])*(T[6]*T[7])*IN[17]+(T[1]*T[4]+T[4]*T[1])*(T[3]*T[7]+T[6]*T[4])*IN[23]+(T[1]*T[4]+T[4]*T[1])*(T[6]*T[1]+T[0]*T[7])*IN[29]+(T[1]*T[4]+T[4]*T[1])*(T[0]*T[4]+T[3]*T[1])*IN[35]);
S[2]=((T[2]*T[2])*(T[0]*T[0])*IN[0]+(T[2]*T[2])*(T[3]*T[3])*IN[6]+(T[2]*T[2])*(T[6]*T[6])*IN[12]+(T[2]*T[2])*(T[3]*T[6]+T[6]*T[3])*IN[18]+(T[2]*T[2])*(T[6]*T[0]+T[0]*T[6])*IN[24]+(T[2]*T[2])*(T[0]*T[3]+T[3]*T[0])*IN[30]+(T[5]*T[5])*(T[0]*T[0])*IN[1]+(T[5]*T[5])*(T[3]*T[3])*IN[7]+(T[5]*T[5])*(T[6]*T[6])*IN[13]+(T[5]*T[5])*(T[3]*T[6]+T[6]*T[3])*IN[19]+(T[5]*T[5])*(T[6]*T[0]+T[0]*T[6])*IN[25]+(T[5]*T[5])*(T[0]*T[3]+T[3]*T[0])*IN[31]+(T[8]*T[8])*(T[0]*T[0])*IN[2]+(T[8]*T[8])*(T[3]*T[3])*IN[8]+(T[8]*T[8])*(T[6]*T[6])*IN[14]+(T[8]*T[8])*(T[3]*T[6]+T[6]*T[3])*IN[20]+(T[8]*T[8])*(T[6]*T[0]+T[0]*T[6])*IN[26]+(T[8]*T[8])*(T[0]*T[3]+T[3]*T[0])*IN[32]+(T[5]*T[8]+T[8]*T[5])*(T[0]*T[0])*IN[3]+(T[5]*T[8]+T[8]*T[5])*(T[3]*T[3])*IN[9]+(T[5]*T[8]+T[8]*T[5])*(T[6]*T[6])*IN[15]+(T[5]*T[8]+T[8]*T[5])*(T[3]*T[6]+T[6]*T[3])*IN[21]+(T[5]*T[8]+T[8]*T[5])*(T[6]*T[0]+T[0]*T[6])*IN[27]+(T[5]*T[8]+T[8]*T[5])*(T[0]*T[3]+T[3]*T[0])*IN[33]+(T[8]*T[2]+T[2]*T[8])*(T[0]*T[0])*IN[4]+(T[8]*T[2]+T[2]*T[8])*(T[3]*T[3])*IN[10]+(T[8]*T[2]+T[2]*T[8])*(T[6]*T[6])*IN[16]+(T[8]*T[2]+T[2]*T[8])*(T[3]*T[6]+T[6]*T[3])*IN[22]+(T[8]*T[2]+T[2]*T[8])*(T[6]*T[0]+T[0]*T[6])*IN[28]+(T[8]*T[2]+T[2]*T[8])*(T[0]*T[3]+T[3]*T[0])*IN[34]+(T[2]*T[5]+T[5]*T[2])*(T[0]*T[0])*IN[5]+(T[2]*T[5]+T[5]*T[2])*(T[3]*T[3])*IN[11]+(T[2]*T[5]+T[5]*T[2])*(T[6]*T[6])*IN[17]+(T[2]*T[5]+T[5]*T[2])*(T[3]*T[6]+T[6]*T[3])*IN[23]+(T[2]*T[5]+T[5]*T[2])*(T[6]*T[0]+T[0]*T[6])*IN[29]+(T[2]*T[5]+T[5]*T[2])*(T[0]*T[3]+T[3]*T[0])*IN[35]);
S[8]=((T[2]*T[2])*(T[1]*T[1])*IN[0]+(T[2]*T[2])*(T[4]*T[4])*IN[6]+(T[2]*T[2])*(T[7]*T[7])*IN[12]+(T[2]*T[2])*(T[4]*T[7]+T[7]*T[4])*IN[18]+(T[2]*T[2])*(T[7]*T[1]+T[1]*T[7])*IN[24]+(T[2]*T[2])*(T[1]*T[4]+T[4]*T[1])*IN[30]+(T[5]*T[5])*(T[1]*T[1])*IN[1]+(T[5]*T[5])*(T[4]*T[4])*IN[7]+(T[5]*T[5])*(T[7]*T[7])*IN[13]+(T[5]*T[5])*(T[4]*T[7]+T[7]*T[4])*IN[19]+(T[5]*T[5])*(T[7]*T[1]+T[1]*T[7])*IN[25]+(T[5]*T[5])*(T[1]*T[4]+T[4]*T[1])*IN[31]+(T[8]*T[8])*(T[1]*T[1])*IN[2]+(T[8]*T[8])*(T[4]*T[4])*IN[8]+(T[8]*T[8])*(T[7]*T[7])*IN[14]+(T[8]*T[8])*(T[4]*T[7]+T[7]*T[4])*IN[20]+(T[8]*T[8])*(T[7]*T[1]+T[1]*T[7])*IN[26]+(T[8]*T[8])*(T[1]*T[4]+T[4]*T[1])*IN[32]+(T[5]*T[8]+T[8]*T[5])*(T[1]*T[1])*IN[3]+(T[5]*T[8]+T[8]*T[5])*(T[4]*T[4])*IN[9]+(T[5]*T[8]+T[8]*T[5])*(T[7]*T[7])*IN[15]+(T[5]*T[8]+T[8]*T[5])*(T[4]*T[7]+T[7]*T[4])*IN[21]+(T[5]*T[8]+T[8]*T[5])*(T[7]*T[1]+T[1]*T[7])*IN[27]+(T[5]*T[8]+T[8]*T[5])*(T[1]*T[4]+T[4]*T[1])*IN[33]+(T[8]*T[2]+T[2]*T[8])*(T[1]*T[1])*IN[4]+(T[8]*T[2]+T[2]*T[8])*(T[4]*T[4])*IN[10]+(T[8]*T[2]+T[2]*T[8])*(T[7]*T[7])*IN[16]+(T[8]*T[2]+T[2]*T[8])*(T[4]*T[7]+T[7]*T[4])*IN[22]+(T[8]*T[2]+T[2]*T[8])*(T[7]*T[1]+T[1]*T[7])*IN[28]+(T[8]*T[2]+T[2]*T[8])*(T[1]*T[4]+T[4]*T[1])*IN[34]+(T[2]*T[5]+T[5]*T[2])*(T[1]*T[1])*IN[5]+(T[2]*T[5]+T[5]*T[2])*(T[4]*T[4])*IN[11]+(T[2]*T[5]+T[5]*T[2])*(T[7]*T[7])*IN[17]+(T[2]*T[5]+T[5]*T[2])*(T[4]*T[7]+T[7]*T[4])*IN[23]+(T[2]*T[5]+T[5]*T[2])*(T[7]*T[1]+T[1]*T[7])*IN[29]+(T[2]*T[5]+T[5]*T[2])*(T[1]*T[4]+T[4]*T[1])*IN[35]);
S[14]=((T[2]*T[2])*(T[2]*T[2])*IN[0]+(T[2]*T[2])*(T[5]*T[5])*IN[6]+(T[2]*T[2])*(T[8]*T[8])*IN[12]+(T[2]*T[2])*(T[5]*T[8]+T[8]*T[5])*IN[18]+(T[2]*T[2])*(T[8]*T[2]+T[2]*T[8])*IN[24]+(T[2]*T[2])*(T[2]*T[5]+T[5]*T[2])*IN[30]+(T[5]*T[5])*(T[2]*T[2])*IN[1]+(T[5]*T[5])*(T[5]*T[5])*IN[7]+(T[5]*T[5])*(T[8]*T[8])*IN[13]+(T[5]*T[5])*(T[5]*T[8]+T[8]*T[5])*IN[19]+(T[5]*T[5])*(T[8]*T[2]+T[2]*T[8])*IN[25]+(T[5]*T[5])*(T[2]*T[5]+T[5]*T[2])*IN[31]+(T[8]*T[8])*(T[2]*T[2])*IN[2]+(T[8]*T[8])*(T[5]*T[5])*IN[8]+(T[8]*T[8])*(T[8]*T[8])*IN[14]+(T[8]*T[8])*(T[5]*T[8]+T[8]*T[5])*IN[20]+(T[8]*T[8])*(T[8]*T[2]+T[2]*T[8])*IN[26]+(T[8]*T[8])*(T[2]*T[5]+T[5]*T[2])*IN[32]+(T[5]*T[8]+T[8]*T[5])*(T[2]*T[2])*IN[3]+(T[5]*T[8]+T[8]*T[5])*(T[5]*T[5])*IN[9]+(T[5]*T[8]+T[8]*T[5])*(T[8]*T[8])*IN[15]+(T[5]*T[8]+T[8]*T[5])*(T[5]*T[8]+T[8]*T[5])*IN[21]+(T[5]*T[8]+T[8]*T[5])*(T[8]*T[2]+T[2]*T[8])*IN[27]+(T[5]*T[8]+T[8]*T[5])*(T[2]*T[5]+T[5]*T[2])*IN[33]+(T[8]*T[2]+T[2]*T[8])*(T[2]*T[2])*IN[4]+(T[8]*T[2]+T[2]*T[8])*(T[5]*T[5])*IN[10]+(T[8]*T[2]+T[2]*T[8])*(T[8]*T[8])*IN[16]+(T[8]*T[2]+T[2]*T[8])*(T[5]*T[8]+T[8]*T[5])*IN[22]+(T[8]*T[2]+T[2]*T[8])*(T[8]*T[2]+T[2]*T[8])*IN[28]+(T[8]*T[2]+T[2]*T[8])*(T[2]*T[5]+T[5]*T[2])*IN[34]+(T[2]*T[5]+T[5]*T[2])*(T[2]*T[2])*IN[5]+(T[2]*T[5]+T[5]*T[2])*(T[5]*T[5])*IN[11]+(T[2]*T[5]+T[5]*T[2])*(T[8]*T[8])*IN[17]+(T[2]*T[5]+T[5]*T[2])*(T[5]*T[8]+T[8]*T[5])*IN[23]+(T[2]*T[5]+T[5]*T[2])*(T[8]*T[2]+T[2]*T[8])*IN[29]+(T[2]*T[5]+T[5]*T[2])*(T[2]*T[5]+T[5]*T[2])*IN[35]);
S[20]=((T[2]*T[2])*(T[1]*T[2])*IN[0]+(T[2]*T[2])*(T[4]*T[5])*IN[6]+(T[2]*T[2])*(T[7]*T[8])*IN[12]+(T[2]*T[2])*(T[4]*T[8]+T[7]*T[5])*IN[18]+(T[2]*T[2])*(T[7]*T[2]+T[1]*T[8])*IN[24]+(T[2]*T[2])*(T[1]*T[5]+T[4]*T[2])*IN[30]+(T[5]*T[5])*(T[1]*T[2])*IN[1]+(T[5]*T[5])*(T[4]*T[5])*IN[7]+(T[5]*T[5])*(T[7]*T[8])*IN[13]+(T[5]*T[5])*(T[4]*T[8]+T[7]*T[5])*IN[19]+(T[5]*T[5])*(T[7]*T[2]+T[1]*T[8])*IN[25]+(T[5]*T[5])*(T[1]*T[5]+T[4]*T[2])*IN[31]+(T[8]*T[8])*(T[1]*T[2])*IN[2]+(T[8]*T[8])*(T[4]*T[5])*IN[8]+(T[8]*T[8])*(T[7]*T[8])*IN[14]+(T[8]*T[8])*(T[4]*T[8]+T[7]*T[5])*IN[20]+(T[8]*T[8])*(T[7]*T[2]+T[1]*T[8])*IN[26]+(T[8]*T[8])*(T[1]*T[5]+T[4]*T[2])*IN[32]+(T[5]*T[8]+T[8]*T[5])*(T[1]*T[2])*IN[3]+(T[5]*T[8]+T[8]*T[5])*(T[4]*T[5])*IN[9]+(T[5]*T[8]+T[8]*T[5])*(T[7]*T[8])*IN[15]+(T[5]*T[8]+T[8]*T[5])*(T[4]*T[8]+T[7]*T[5])*IN[21]+(T[5]*T[8]+T[8]*T[5])*(T[7]*T[2]+T[1]*T[8])*IN[27]+(T[5]*T[8]+T[8]*T[5])*(T[1]*T[5]+T[4]*T[2])*IN[33]+(T[8]*T[2]+T[2]*T[8])*(T[1]*T[2])*IN[4]+(T[8]*T[2]+T[2]*T[8])*(T[4]*T[5])*IN[10]+(T[8]*T[2]+T[2]*T[8])*(T[7]*T[8])*IN[16]+(T[8]*T[2]+T[2]*T[8])*(T[4]*T[8]+T[7]*T[5])*IN[22]+(T[8]*T[2]+T[2]*T[8])*(T[7]*T[2]+T[1]*T[8])*IN[28]+(T[8]*T[2]+T[2]*T[8])*(T[1]*T[5]+T[4]*T[2])*IN[34]+(T[2]*T[5]+T[5]*T[2])*(T[1]*T[2])*IN[5]+(T[2]*T[5]+T[5]*T[2])*(T[4]*T[5])*IN[11]+(T[2]*T[5]+T[5]*T[2])*(T[7]*T[8])*IN[17]+(T[2]*T[5]+T[5]*T[2])*(T[4]*T[8]+T[7]*T[5])*IN[23]+(T[2]*T[5]+T[5]*T[2])*(T[7]*T[2]+T[1]*T[8])*IN[29]+(T[2]*T[5]+T[5]*T[2])*(T[1]*T[5]+T[4]*T[2])*IN[35]);
S[26]=((T[2]*T[2])*(T[2]*T[0])*IN[0]+(T[2]*T[2])*(T[5]*T[3])*IN[6]+(T[2]*T[2])*(T[8]*T[6])*IN[12]+(T[2]*T[2])*(T[5]*T[6]+T[8]*T[3])*IN[18]+(T[2]*T[2])*(T[8]*T[0]+T[2]*T[6])*IN[24]+(T[2]*T[2])*(T[2]*T[3]+T[5]*T[0])*IN[30]+(T[5]*T[5])*(T[2]*T[0])*IN[1]+(T[5]*T[5])*(T[5]*T[3])*IN[7]+(T[5]*T[5])*(T[8]*T[6])*IN[13]+(T[5]*T[5])*(T[5]*T[6]+T[8]*T[3])*IN[19]+(T[5]*T[5])*(T[8]*T[0]+T[2]*T[6])*IN[25]+(T[5]*T[5])*(T[2]*T[3]+T[5]*T[0])*IN[31]+(T[8]*T[8])*(T[2]*T[0])*IN[2]+(T[8]*T[8])*(T[5]*T[3])*IN[8]+(T[8]*T[8])*(T[8]*T[6])*IN[14]+(T[8]*T[8])*(T[5]*T[6]+T[8]*T[3])*IN[20]+(T[8]*T[8])*(T[8]*T[0]+T[2]*T[6])*IN[26]+(T[8]*T[8])*(T[2]*T[3]+T[5]*T[0])*IN[32]+(T[5]*T[8]+T[8]*T[5])*(T[2]*T[0])*IN[3]+(T[5]*T[8]+T[8]*T[5])*(T[5]*T[3])*IN[9]+(T[5]*T[8]+T[8]*T[5])*(T[8]*T[6])*IN[15]+(T[5]*T[8]+T[8]*T[5])*(T[5]*T[6]+T[8]*T[3])*IN[21]+(T[5]*T[8]+T[8]*T[5])*(T[8]*T[0]+T[2]*T[6])*IN[27]+(T[5]*T[8]+T[8]*T[5])*(T[2]*T[3]+T[5]*T[0])*IN[33]+(T[8]*T[2]+T[2]*T[8])*(T[2]*T[0])*IN[4]+(T[8]*T[2]+T[2]*T[8])*(T[5]*T[3])*IN[10]+(T[8]*T[2]+T[2]*T[8])*(T[8]*T[6])*IN[16]+(T[8]*T[2]+T[2]*T[8])*(T[5]*T[6]+T[8]*T[3])*IN[22]+(T[8]*T[2]+T[2]*T[8])*(T[8]*T[0]+T[2]*T[6])*IN[28]+(T[8]*T[2]+T[2]*T[8])*(T[2]*T[3]+T[5]*T[0])*IN[34]+(T[2]*T[5]+T[5]*T[2])*(T[2]*T[0])*IN[5]+(T[2]*T[5]+T[5]*T[2])*(T[5]*T[3])*IN[11]+(T[2]*T[5]+T[5]*T[2])*(T[8]*T[6])*IN[17]+(T[2]*T[5]+T[5]*T[2])*(T[5]*T[6]+T[8]*T[3])*IN[23]+(T[2]*T[5]+T[5]*T[2])*(T[8]*T[0]+T[2]*T[6])*IN[29]+(T[2]*T[5]+T[5]*T[2])*(T[2]*T[3]+T[5]*T[0])*IN[35]);
S[32]=((T[2]*T[2])*(T[0]*T[1])*IN[0]+(T[2]*T[2])*(T[3]*T[4])*IN[6]+(T[2]*T[2])*(T[6]*T[7])*IN[12]+(T[2]*T[2])*(T[3]*T[7]+T[6]*T[4])*IN[18]+(T[2]*T[2])*(T[6]*T[1]+T[0]*T[7])*IN[24]+(T[2]*T[2])*(T[0]*T[4]+T[3]*T[1])*IN[30]+(T[5]*T[5])*(T[0]*T[1])*IN[1]+(T[5]*T[5])*(T[3]*T[4])*IN[7]+(T[5]*T[5])*(T[6]*T[7])*IN[13]+(T[5]*T[5])*(T[3]*T[7]+T[6]*T[4])*IN[19]+(T[5]*T[5])*(T[6]*T[1]+T[0]*T[7])*IN[25]+(T[5]*T[5])*(T[0]*T[4]+T[3]*T[1])*IN[31]+(T[8]*T[8])*(T[0]*T[1])*IN[2]+(T[8]*T[8])*(T[3]*T[4])*IN[8]+(T[8]*T[8])*(T[6]*T[7])*IN[14]+(T[8]*T[8])*(T[3]*T[7]+T[6]*T[4])*IN[20]+(T[8]*T[8])*(T[6]*T[1]+T[0]*T[7])*IN[26]+(T[8]*T[8])*(T[0]*T[4]+T[3]*T[1])*IN[32]+(T[5]*T[8]+T[8]*T[5])*(T[0]*T[1])*IN[3]+(T[5]*T[8]+T[8]*T[5])*(T[3]*T[4])*IN[9]+(T[5]*T[8]+T[8]*T[5])*(T[6]*T[7])*IN[15]+(T[5]*T[8]+T[8]*T[5])*(T[3]*T[7]+T[6]*T[4])*IN[21]+(T[5]*T[8]+T[8]*T[5])*(T[6]*T[1]+T[0]*T[7])*IN[27]+(T[5]*T[8]+T[8]*T[5])*(T[0]*T[4]+T[3]*T[1])*IN[33]+(T[8]*T[2]+T[2]*T[8])*(T[0]*T[1])*IN[4]+(T[8]*T[2]+T[2]*T[8])*(T[3]*T[4])*IN[10]+(T[8]*T[2]+T[2]*T[8])*(T[6]*T[7])*IN[16]+(T[8]*T[2]+T[2]*T[8])*(T[3]*T[7]+T[6]*T[4])*IN[22]+(T[8]*T[2]+T[2]*T[8])*(T[6]*T[1]+T[0]*T[7])*IN[28]+(T[8]*T[2]+T[2]*T[8])*(T[0]*T[4]+T[3]*T[1])*IN[34]+(T[2]*T[5]+T[5]*T[2])*(T[0]*T[1])*IN[5]+(T[2]*T[5]+T[5]*T[2])*(T[3]*T[4])*IN[11]+(T[2]*T[5]+T[5]*T[2])*(T[6]*T[7])*IN[17]+(T[2]*T[5]+T[5]*T[2])*(T[3]*T[7]+T[6]*T[4])*IN[23]+(T[2]*T[5]+T[5]*T[2])*(T[6]*T[1]+T[0]*T[7])*IN[29]+(T[2]*T[5]+T[5]*T[2])*(T[0]*T[4]+T[3]*T[1])*IN[35]);
S[3]=((T[1]*T[2])*(T[0]*T[0])*IN[0]+(T[1]*T[2])*(T[3]*T[3])*IN[6]+(T[1]*T[2])*(T[6]*T[6])*IN[12]+(T[1]*T[2])*(T[3]*T[6]+T[6]*T[3])*IN[18]+(T[1]*T[2])*(T[6]*T[0]+T[0]*T[6])*IN[24]+(T[1]*T[2])*(T[0]*T[3]+T[3]*T[0])*IN[30]+(T[4]*T[5])*(T[0]*T[0])*IN[1]+(T[4]*T[5])*(T[3]*T[3])*IN[7]+(T[4]*T[5])*(T[6]*T[6])*IN[13]+(T[4]*T[5])*(T[3]*T[6]+T[6]*T[3])*IN[19]+(T[4]*T[5])*(T[6]*T[0]+T[0]*T[6])*IN[25]+(T[4]*T[5])*(T[0]*T[3]+T[3]*T[0])*IN[31]+(T[7]*T[8])*(T[0]*T[0])*IN[2]+(T[7]*T[8])*(T[3]*T[3])*IN[8]+(T[7]*T[8])*(T[6]*T[6])*IN[14]+(T[7]*T[8])*(T[3]*T[6]+T[6]*T[3])*IN[20]+(T[7]*T[8])*(T[6]*T[0]+T[0]*T[6])*IN[26]+(T[7]*T[8])*(T[0]*T[3]+T[3]*T[0])*IN[32]+(T[4]*T[8]+T[7]*T[5])*(T[0]*T[0])*IN[3]+(T[4]*T[8]+T[7]*T[5])*(T[3]*T[3])*IN[9]+(T[4]*T[8]+T[7]*T[5])*(T[6]*T[6])*IN[15]+(T[4]*T[8]+T[7]*T[5])*(T[3]*T[6]+T[6]*T[3])*IN[21]+(T[4]*T[8]+T[7]*T[5])*(T[6]*T[0]+T[0]*T[6])*IN[27]+(T[4]*T[8]+T[7]*T[5])*(T[0]*T[3]+T[3]*T[0])*IN[33]+(T[7]*T[2]+T[1]*T[8])*(T[0]*T[0])*IN[4]+(T[7]*T[2]+T[1]*T[8])*(T[3]*T[3])*IN[10]+(T[7]*T[2]+T[1]*T[8])*(T[6]*T[6])*IN[16]+(T[7]*T[2]+T[1]*T[8])*(T[3]*T[6]+T[6]*T[3])*IN[22]+(T[7]*T[2]+T[1]*T[8])*(T[6]*T[0]+T[0]*T[6])*IN[28]+(T[7]*T[2]+T[1]*T[8])*(T[0]*T[3]+T[3]*T[0])*IN[34]+(T[1]*T[5]+T[4]*T[2])*(T[0]*T[0])*IN[5]+(T[1]*T[5]+T[4]*T[2])*(T[3]*T[3])*IN[11]+(T[1]*T[5]+T[4]*T[2])*(T[6]*T[6])*IN[17]+(T[1]*T[5]+T[4]*T[2])*(T[3]*T[6]+T[6]*T[3])*IN[23]+(T[1]*T[5]+T[4]*T[2])*(T[6]*T[0]+T[0]*T[6])*IN[29]+(T[1]*T[5]+T[4]*T[2])*(T[0]*T[3]+T[3]*T[0])*IN[35]);
S[9]=((T[1]*T[2])*(T[1]*T[1])*IN[0]+(T[1]*T[2])*(T[4]*T[4])*IN[6]+(T[1]*T[2])*(T[7]*T[7])*IN[12]+(T[1]*T[2])*(T[4]*T[7]+T[7]*T[4])*IN[18]+(T[1]*T[2])*(T[7]*T[1]+T[1]*T[7])*IN[24]+(T[1]*T[2])*(T[1]*T[4]+T[4]*T[1])*IN[30]+(T[4]*T[5])*(T[1]*T[1])*IN[1]+(T[4]*T[5])*(T[4]*T[4])*IN[7]+(T[4]*T[5])*(T[7]*T[7])*IN[13]+(T[4]*T[5])*(T[4]*T[7]+T[7]*T[4])*IN[19]+(T[4]*T[5])*(T[7]*T[1]+T[1]*T[7])*IN[25]+(T[4]*T[5])*(T[1]*T[4]+T[4]*T[1])*IN[31]+(T[7]*T[8])*(T[1]*T[1])*IN[2]+(T[7]*T[8])*(T[4]*T[4])*IN[8]+(T[7]*T[8])*(T[7]*T[7])*IN[14]+(T[7]*T[8])*(T[4]*T[7]+T[7]*T[4])*IN[20]+(T[7]*T[8])*(T[7]*T[1]+T[1]*T[7])*IN[26]+(T[7]*T[8])*(T[1]*T[4]+T[4]*T[1])*IN[32]+(T[4]*T[8]+T[7]*T[5])*(T[1]*T[1])*IN[3]+(T[4]*T[8]+T[7]*T[5])*(T[4]*T[4])*IN[9]+(T[4]*T[8]+T[7]*T[5])*(T[7]*T[7])*IN[15]+(T[4]*T[8]+T[7]*T[5])*(T[4]*T[7]+T[7]*T[4])*IN[21]+(T[4]*T[8]+T[7]*T[5])*(T[7]*T[1]+T[1]*T[7])*IN[27]+(T[4]*T[8]+T[7]*T[5])*(T[1]*T[4]+T[4]*T[1])*IN[33]+(T[7]*T[2]+T[1]*T[8])*(T[1]*T[1])*IN[4]+(T[7]*T[2]+T[1]*T[8])*(T[4]*T[4])*IN[10]+(T[7]*T[2]+T[1]*T[8])*(T[7]*T[7])*IN[16]+(T[7]*T[2]+T[1]*T[8])*(T[4]*T[7]+T[7]*T[4])*IN[22]+(T[7]*T[2]+T[1]*T[8])*(T[7]*T[1]+T[1]*T[7])*IN[28]+(T[7]*T[2]+T[1]*T[8])*(T[1]*T[4]+T[4]*T[1])*IN[34]+(T[1]*T[5]+T[4]*T[2])*(T[1]*T[1])*IN[5]+(T[1]*T[5]+T[4]*T[2])*(T[4]*T[4])*IN[11]+(T[1]*T[5]+T[4]*T[2])*(T[7]*T[7])*IN[17]+(T[1]*T[5]+T[4]*T[2])*(T[4]*T[7]+T[7]*T[4])*IN[23]+(T[1]*T[5]+T[4]*T[2])*(T[7]*T[1]+T[1]*T[7])*IN[29]+(T[1]*T[5]+T[4]*T[2])*(T[1]*T[4]+T[4]*T[1])*IN[35]);
S[15]=((T[1]*T[2])*(T[2]*T[2])*IN[0]+(T[1]*T[2])*(T[5]*T[5])*IN[6]+(T[1]*T[2])*(T[8]*T[8])*IN[12]+(T[1]*T[2])*(T[5]*T[8]+T[8]*T[5])*IN[18]+(T[1]*T[2])*(T[8]*T[2]+T[2]*T[8])*IN[24]+(T[1]*T[2])*(T[2]*T[5]+T[5]*T[2])*IN[30]+(T[4]*T[5])*(T[2]*T[2])*IN[1]+(T[4]*T[5])*(T[5]*T[5])*IN[7]+(T[4]*T[5])*(T[8]*T[8])*IN[13]+(T[4]*T[5])*(T[5]*T[8]+T[8]*T[5])*IN[19]+(T[4]*T[5])*(T[8]*T[2]+T[2]*T[8])*IN[25]+(T[4]*T[5])*(T[2]*T[5]+T[5]*T[2])*IN[31]+(T[7]*T[8])*(T[2]*T[2])*IN[2]+(T[7]*T[8])*(T[5]*T[5])*IN[8]+(T[7]*T[8])*(T[8]*T[8])*IN[14]+(T[7]*T[8])*(T[5]*T[8]+T[8]*T[5])*IN[20]+(T[7]*T[8])*(T[8]*T[2]+T[2]*T[8])*IN[26]+(T[7]*T[8])*(T[2]*T[5]+T[5]*T[2])*IN[32]+(T[4]*T[8]+T[7]*T[5])*(T[2]*T[2])*IN[3]+(T[4]*T[8]+T[7]*T[5])*(T[5]*T[5])*IN[9]+(T[4]*T[8]+T[7]*T[5])*(T[8]*T[8])*IN[15]+(T[4]*T[8]+T[7]*T[5])*(T[5]*T[8]+T[8]*T[5])*IN[21]+(T[4]*T[8]+T[7]*T[5])*(T[8]*T[2]+T[2]*T[8])*IN[27]+(T[4]*T[8]+T[7]*T[5])*(T[2]*T[5]+T[5]*T[2])*IN[33]+(T[7]*T[2]+T[1]*T[8])*(T[2]*T[2])*IN[4]+(T[7]*T[2]+T[1]*T[8])*(T[5]*T[5])*IN[10]+(T[7]*T[2]+T[1]*T[8])*(T[8]*T[8])*IN[16]+(T[7]*T[2]+T[1]*T[8])*(T[5]*T[8]+T[8]*T[5])*IN[22]+(T[7]*T[2]+T[1]*T[8])*(T[8]*T[2]+T[2]*T[8])*IN[28]+(T[7]*T[2]+T[1]*T[8])*(T[2]*T[5]+T[5]*T[2])*IN[34]+(T[1]*T[5]+T[4]*T[2])*(T[2]*T[2])*IN[5]+(T[1]*T[5]+T[4]*T[2])*(T[5]*T[5])*IN[11]+(T[1]*T[5]+T[4]*T[2])*(T[8]*T[8])*IN[17]+(T[1]*T[5]+T[4]*T[2])*(T[5]*T[8]+T[8]*T[5])*IN[23]+(T[1]*T[5]+T[4]*T[2])*(T[8]*T[2]+T[2]*T[8])*IN[29]+(T[1]*T[5]+T[4]*T[2])*(T[2]*T[5]+T[5]*T[2])*IN[35]);
S[21]=((T[1]*T[2])*(T[1]*T[2])*IN[0]+(T[1]*T[2])*(T[4]*T[5])*IN[6]+(T[1]*T[2])*(T[7]*T[8])*IN[12]+(T[1]*T[2])*(T[4]*T[8]+T[7]*T[5])*IN[18]+(T[1]*T[2])*(T[7]*T[2]+T[1]*T[8])*IN[24]+(T[1]*T[2])*(T[1]*T[5]+T[4]*T[2])*IN[30]+(T[4]*T[5])*(T[1]*T[2])*IN[1]+(T[4]*T[5])*(T[4]*T[5])*IN[7]+(T[4]*T[5])*(T[7]*T[8])*IN[13]+(T[4]*T[5])*(T[4]*T[8]+T[7]*T[5])*IN[19]+(T[4]*T[5])*(T[7]*T[2]+T[1]*T[8])*IN[25]+(T[4]*T[5])*(T[1]*T[5]+T[4]*T[2])*IN[31]+(T[7]*T[8])*(T[1]*T[2])*IN[2]+(T[7]*T[8])*(T[4]*T[5])*IN[8]+(T[7]*T[8])*(T[7]*T[8])*IN[14]+(T[7]*T[8])*(T[4]*T[8]+T[7]*T[5])*IN[20]+(T[7]*T[8])*(T[7]*T[2]+T[1]*T[8])*IN[26]+(T[7]*T[8])*(T[1]*T[5]+T[4]*T[2])*IN[32]+(T[4]*T[8]+T[7]*T[5])*(T[1]*T[2])*IN[3]+(T[4]*T[8]+T[7]*T[5])*(T[4]*T[5])*IN[9]+(T[4]*T[8]+T[7]*T[5])*(T[7]*T[8])*IN[15]+(T[4]*T[8]+T[7]*T[5])*(T[4]*T[8]+T[7]*T[5])*IN[21]+(T[4]*T[8]+T[7]*T[5])*(T[7]*T[2]+T[1]*T[8])*IN[27]+(T[4]*T[8]+T[7]*T[5])*(T[1]*T[5]+T[4]*T[2])*IN[33]+(T[7]*T[2]+T[1]*T[8])*(T[1]*T[2])*IN[4]+(T[7]*T[2]+T[1]*T[8])*(T[4]*T[5])*IN[10]+(T[7]*T[2]+T[1]*T[8])*(T[7]*T[8])*IN[16]+(T[7]*T[2]+T[1]*T[8])*(T[4]*T[8]+T[7]*T[5])*IN[22]+(T[7]*T[2]+T[1]*T[8])*(T[7]*T[2]+T[1]*T[8])*IN[28]+(T[7]*T[2]+T[1]*T[8])*(T[1]*T[5]+T[4]*T[2])*IN[34]+(T[1]*T[5]+T[4]*T[2])*(T[1]*T[2])*IN[5]+(T[1]*T[5]+T[4]*T[2])*(T[4]*T[5])*IN[11]+(T[1]*T[5]+T[4]*T[2])*(T[7]*T[8])*IN[17]+(T[1]*T[5]+T[4]*T[2])*(T[4]*T[8]+T[7]*T[5])*IN[23]+(T[1]*T[5]+T[4]*T[2])*(T[7]*T[2]+T[1]*T[8])*IN[29]+(T[1]*T[5]+T[4]*T[2])*(T[1]*T[5]+T[4]*T[2])*IN[35]);
S[27]=((T[1]*T[2])*(T[2]*T[0])*IN[0]+(T[1]*T[2])*(T[5]*T[3])*IN[6]+(T[1]*T[2])*(T[8]*T[6])*IN[12]+(T[1]*T[2])*(T[5]*T[6]+T[8]*T[3])*IN[18]+(T[1]*T[2])*(T[8]*T[0]+T[2]*T[6])*IN[24]+(T[1]*T[2])*(T[2]*T[3]+T[5]*T[0])*IN[30]+(T[4]*T[5])*(T[2]*T[0])*IN[1]+(T[4]*T[5])*(T[5]*T[3])*IN[7]+(T[4]*T[5])*(T[8]*T[6])*IN[13]+(T[4]*T[5])*(T[5]*T[6]+T[8]*T[3])*IN[19]+(T[4]*T[5])*(T[8]*T[0]+T[2]*T[6])*IN[25]+(T[4]*T[5])*(T[2]*T[3]+T[5]*T[0])*IN[31]+(T[7]*T[8])*(T[2]*T[0])*IN[2]+(T[7]*T[8])*(T[5]*T[3])*IN[8]+(T[7]*T[8])*(T[8]*T[6])*IN[14]+(T[7]*T[8])*(T[5]*T[6]+T[8]*T[3])*IN[20]+(T[7]*T[8])*(T[8]*T[0]+T[2]*T[6])*IN[26]+(T[7]*T[8])*(T[2]*T[3]+T[5]*T[0])*IN[32]+(T[4]*T[8]+T[7]*T[5])*(T[2]*T[0])*IN[3]+(T[4]*T[8]+T[7]*T[5])*(T[5]*T[3])*IN[9]+(T[4]*T[8]+T[7]*T[5])*(T[8]*T[6])*IN[15]+(T[4]*T[8]+T[7]*T[5])*(T[5]*T[6]+T[8]*T[3])*IN[21]+(T[4]*T[8]+T[7]*T[5])*(T[8]*T[0]+T[2]*T[6])*IN[27]+(T[4]*T[8]+T[7]*T[5])*(T[2]*T[3]+T[5]*T[0])*IN[33]+(T[7]*T[2]+T[1]*T[8])*(T[2]*T[0])*IN[4]+(T[7]*T[2]+T[1]*T[8])*(T[5]*T[3])*IN[10]+(T[7]*T[2]+T[1]*T[8])*(T[8]*T[6])*IN[16]+(T[7]*T[2]+T[1]*T[8])*(T[5]*T[6]+T[8]*T[3])*IN[22]+(T[7]*T[2]+T[1]*T[8])*(T[8]*T[0]+T[2]*T[6])*IN[28]+(T[7]*T[2]+T[1]*T[8])*(T[2]*T[3]+T[5]*T[0])*IN[34]+(T[1]*T[5]+T[4]*T[2])*(T[2]*T[0])*IN[5]+(T[1]*T[5]+T[4]*T[2])*(T[5]*T[3])*IN[11]+(T[1]*T[5]+T[4]*T[2])*(T[8]*T[6])*IN[17]+(T[1]*T[5]+T[4]*T[2])*(T[5]*T[6]+T[8]*T[3])*IN[23]+(T[1]*T[5]+T[4]*T[2])*(T[8]*T[0]+T[2]*T[6])*IN[29]+(T[1]*T[5]+T[4]*T[2])*(T[2]*T[3]+T[5]*T[0])*IN[35]);
S[33]=((T[1]*T[2])*(T[0]*T[1])*IN[0]+(T[1]*T[2])*(T[3]*T[4])*IN[6]+(T[1]*T[2])*(T[6]*T[7])*IN[12]+(T[1]*T[2])*(T[3]*T[7]+T[6]*T[4])*IN[18]+(T[1]*T[2])*(T[6]*T[1]+T[0]*T[7])*IN[24]+(T[1]*T[2])*(T[0]*T[4]+T[3]*T[1])*IN[30]+(T[4]*T[5])*(T[0]*T[1])*IN[1]+(T[4]*T[5])*(T[3]*T[4])*IN[7]+(T[4]*T[5])*(T[6]*T[7])*IN[13]+(T[4]*T[5])*(T[3]*T[7]+T[6]*T[4])*IN[19]+(T[4]*T[5])*(T[6]*T[1]+T[0]*T[7])*IN[25]+(T[4]*T[5])*(T[0]*T[4]+T[3]*T[1])*IN[31]+(T[7]*T[8])*(T[0]*T[1])*IN[2]+(T[7]*T[8])*(T[3]*T[4])*IN[8]+(T[7]*T[8])*(T[6]*T[7])*IN[14]+(T[7]*T[8])*(T[3]*T[7]+T[6]*T[4])*IN[20]+(T[7]*T[8])*(T[6]*T[1]+T[0]*T[7])*IN[26]+(T[7]*T[8])*(T[0]*T[4]+T[3]*T[1])*IN[32]+(T[4]*T[8]+T[7]*T[5])*(T[0]*T[1])*IN[3]+(T[4]*T[8]+T[7]*T[5])*(T[3]*T[4])*IN[9]+(T[4]*T[8]+T[7]*T[5])*(T[6]*T[7])*IN[15]+(T[4]*T[8]+T[7]*T[5])*(T[3]*T[7]+T[6]*T[4])*IN[21]+(T[4]*T[8]+T[7]*T[5])*(T[6]*T[1]+T[0]*T[7])*IN[27]+(T[4]*T[8]+T[7]*T[5])*(T[0]*T[4]+T[3]*T[1])*IN[33]+(T[7]*T[2]+T[1]*T[8])*(T[0]*T[1])*IN[4]+(T[7]*T[2]+T[1]*T[8])*(T[3]*T[4])*IN[10]+(T[7]*T[2]+T[1]*T[8])*(T[6]*T[7])*IN[16]+(T[7]*T[2]+T[1]*T[8])*(T[3]*T[7]+T[6]*T[4])*IN[22]+(T[7]*T[2]+T[1]*T[8])*(T[6]*T[1]+T[0]*T[7])*IN[28]+(T[7]*T[2]+T[1]*T[8])*(T[0]*T[4]+T[3]*T[1])*IN[34]+(T[1]*T[5]+T[4]*T[2])*(T[0]*T[1])*IN[5]+(T[1]*T[5]+T[4]*T[2])*(T[3]*T[4])*IN[11]+(T[1]*T[5]+T[4]*T[2])*(T[6]*T[7])*IN[17]+(T[1]*T[5]+T[4]*T[2])*(T[3]*T[7]+T[6]*T[4])*IN[23]+(T[1]*T[5]+T[4]*T[2])*(T[6]*T[1]+T[0]*T[7])*IN[29]+(T[1]*T[5]+T[4]*T[2])*(T[0]*T[4]+T[3]*T[1])*IN[35]);
S[4]=((T[2]*T[0])*(T[0]*T[0])*IN[0]+(T[2]*T[0])*(T[3]*T[3])*IN[6]+(T[2]*T[0])*(T[6]*T[6])*IN[12]+(T[2]*T[0])*(T[3]*T[6]+T[6]*T[3])*IN[18]+(T[2]*T[0])*(T[6]*T[0]+T[0]*T[6])*IN[24]+(T[2]*T[0])*(T[0]*T[3]+T[3]*T[0])*IN[30]+(T[5]*T[3])*(T[0]*T[0])*IN[1]+(T[5]*T[3])*(T[3]*T[3])*IN[7]+(T[5]*T[3])*(T[6]*T[6])*IN[13]+(T[5]*T[3])*(T[3]*T[6]+T[6]*T[3])*IN[19]+(T[5]*T[3])*(T[6]*T[0]+T[0]*T[6])*IN[25]+(T[5]*T[3])*(T[0]*T[3]+T[3]*T[0])*IN[31]+(T[8]*T[6])*(T[0]*T[0])*IN[2]+(T[8]*T[6])*(T[3]*T[3])*IN[8]+(T[8]*T[6])*(T[6]*T[6])*IN[14]+(T[8]*T[6])*(T[3]*T[6]+T[6]*T[3])*IN[20]+(T[8]*T[6])*(T[6]*T[0]+T[0]*T[6])*IN[26]+(T[8]*T[6])*(T[0]*T[3]+T[3]*T[0])*IN[32]+(T[5]*T[6]+T[8]*T[3])*(T[0]*T[0])*IN[3]+(T[5]*T[6]+T[8]*T[3])*(T[3]*T[3])*IN[9]+(T[5]*T[6]+T[8]*T[3])*(T[6]*T[6])*IN[15]+(T[5]*T[6]+T[8]*T[3])*(T[3]*T[6]+T[6]*T[3])*IN[21]+(T[5]*T[6]+T[8]*T[3])*(T[6]*T[0]+T[0]*T[6])*IN[27]+(T[5]*T[6]+T[8]*T[3])*(T[0]*T[3]+T[3]*T[0])*IN[33]+(T[8]*T[0]+T[2]*T[6])*(T[0]*T[0])*IN[4]+(T[8]*T[0]+T[2]*T[6])*(T[3]*T[3])*IN[10]+(T[8]*T[0]+T[2]*T[6])*(T[6]*T[6])*IN[16]+(T[8]*T[0]+T[2]*T[6])*(T[3]*T[6]+T[6]*T[3])*IN[22]+(T[8]*T[0]+T[2]*T[6])*(T[6]*T[0]+T[0]*T[6])*IN[28]+(T[8]*T[0]+T[2]*T[6])*(T[0]*T[3]+T[3]*T[0])*IN[34]+(T[2]*T[3]+T[5]*T[0])*(T[0]*T[0])*IN[5]+(T[2]*T[3]+T[5]*T[0])*(T[3]*T[3])*IN[11]+(T[2]*T[3]+T[5]*T[0])*(T[6]*T[6])*IN[17]+(T[2]*T[3]+T[5]*T[0])*(T[3]*T[6]+T[6]*T[3])*IN[23]+(T[2]*T[3]+T[5]*T[0])*(T[6]*T[0]+T[0]*T[6])*IN[29]+(T[2]*T[3]+T[5]*T[0])*(T[0]*T[3]+T[3]*T[0])*IN[35]);
S[10]=((T[2]*T[0])*(T[1]*T[1])*IN[0]+(T[2]*T[0])*(T[4]*T[4])*IN[6]+(T[2]*T[0])*(T[7]*T[7])*IN[12]+(T[2]*T[0])*(T[4]*T[7]+T[7]*T[4])*IN[18]+(T[2]*T[0])*(T[7]*T[1]+T[1]*T[7])*IN[24]+(T[2]*T[0])*(T[1]*T[4]+T[4]*T[1])*IN[30]+(T[5]*T[3])*(T[1]*T[1])*IN[1]+(T[5]*T[3])*(T[4]*T[4])*IN[7]+(T[5]*T[3])*(T[7]*T[7])*IN[13]+(T[5]*T[3])*(T[4]*T[7]+T[7]*T[4])*IN[19]+(T[5]*T[3])*(T[7]*T[1]+T[1]*T[7])*IN[25]+(T[5]*T[3])*(T[1]*T[4]+T[4]*T[1])*IN[31]+(T[8]*T[6])*(T[1]*T[1])*IN[2]+(T[8]*T[6])*(T[4]*T[4])*IN[8]+(T[8]*T[6])*(T[7]*T[7])*IN[14]+(T[8]*T[6])*(T[4]*T[7]+T[7]*T[4])*IN[20]+(T[8]*T[6])*(T[7]*T[1]+T[1]*T[7])*IN[26]+(T[8]*T[6])*(T[1]*T[4]+T[4]*T[1])*IN[32]+(T[5]*T[6]+T[8]*T[3])*(T[1]*T[1])*IN[3]+(T[5]*T[6]+T[8]*T[3])*(T[4]*T[4])*IN[9]+(T[5]*T[6]+T[8]*T[3])*(T[7]*T[7])*IN[15]+(T[5]*T[6]+T[8]*T[3])*(T[4]*T[7]+T[7]*T[4])*IN[21]+(T[5]*T[6]+T[8]*T[3])*(T[7]*T[1]+T[1]*T[7])*IN[27]+(T[5]*T[6]+T[8]*T[3])*(T[1]*T[4]+T[4]*T[1])*IN[33]+(T[8]*T[0]+T[2]*T[6])*(T[1]*T[1])*IN[4]+(T[8]*T[0]+T[2]*T[6])*(T[4]*T[4])*IN[10]+(T[8]*T[0]+T[2]*T[6])*(T[7]*T[7])*IN[16]+(T[8]*T[0]+T[2]*T[6])*(T[4]*T[7]+T[7]*T[4])*IN[22]+(T[8]*T[0]+T[2]*T[6])*(T[7]*T[1]+T[1]*T[7])*IN[28]+(T[8]*T[0]+T[2]*T[6])*(T[1]*T[4]+T[4]*T[1])*IN[34]+(T[2]*T[3]+T[5]*T[0])*(T[1]*T[1])*IN[5]+(T[2]*T[3]+T[5]*T[0])*(T[4]*T[4])*IN[11]+(T[2]*T[3]+T[5]*T[0])*(T[7]*T[7])*IN[17]+(T[2]*T[3]+T[5]*T[0])*(T[4]*T[7]+T[7]*T[4])*IN[23]+(T[2]*T[3]+T[5]*T[0])*(T[7]*T[1]+T[1]*T[7])*IN[29]+(T[2]*T[3]+T[5]*T[0])*(T[1]*T[4]+T[4]*T[1])*IN[35]);
S[16]=((T[2]*T[0])*(T[2]*T[2])*IN[0]+(T[2]*T[0])*(T[5]*T[5])*IN[6]+(T[2]*T[0])*(T[8]*T[8])*IN[12]+(T[2]*T[0])*(T[5]*T[8]+T[8]*T[5])*IN[18]+(T[2]*T[0])*(T[8]*T[2]+T[2]*T[8])*IN[24]+(T[2]*T[0])*(T[2]*T[5]+T[5]*T[2])*IN[30]+(T[5]*T[3])*(T[2]*T[2])*IN[1]+(T[5]*T[3])*(T[5]*T[5])*IN[7]+(T[5]*T[3])*(T[8]*T[8])*IN[13]+(T[5]*T[3])*(T[5]*T[8]+T[8]*T[5])*IN[19]+(T[5]*T[3])*(T[8]*T[2]+T[2]*T[8])*IN[25]+(T[5]*T[3])*(T[2]*T[5]+T[5]*T[2])*IN[31]+(T[8]*T[6])*(T[2]*T[2])*IN[2]+(T[8]*T[6])*(T[5]*T[5])*IN[8]+(T[8]*T[6])*(T[8]*T[8])*IN[14]+(T[8]*T[6])*(T[5]*T[8]+T[8]*T[5])*IN[20]+(T[8]*T[6])*(T[8]*T[2]+T[2]*T[8])*IN[26]+(T[8]*T[6])*(T[2]*T[5]+T[5]*T[2])*IN[32]+(T[5]*T[6]+T[8]*T[3])*(T[2]*T[2])*IN[3]+(T[5]*T[6]+T[8]*T[3])*(T[5]*T[5])*IN[9]+(T[5]*T[6]+T[8]*T[3])*(T[8]*T[8])*IN[15]+(T[5]*T[6]+T[8]*T[3])*(T[5]*T[8]+T[8]*T[5])*IN[21]+(T[5]*T[6]+T[8]*T[3])*(T[8]*T[2]+T[2]*T[8])*IN[27]+(T[5]*T[6]+T[8]*T[3])*(T[2]*T[5]+T[5]*T[2])*IN[33]+(T[8]*T[0]+T[2]*T[6])*(T[2]*T[2])*IN[4]+(T[8]*T[0]+T[2]*T[6])*(T[5]*T[5])*IN[10]+(T[8]*T[0]+T[2]*T[6])*(T[8]*T[8])*IN[16]+(T[8]*T[0]+T[2]*T[6])*(T[5]*T[8]+T[8]*T[5])*IN[22]+(T[8]*T[0]+T[2]*T[6])*(T[8]*T[2]+T[2]*T[8])*IN[28]+(T[8]*T[0]+T[2]*T[6])*(T[2]*T[5]+T[5]*T[2])*IN[34]+(T[2]*T[3]+T[5]*T[0])*(T[2]*T[2])*IN[5]+(T[2]*T[3]+T[5]*T[0])*(T[5]*T[5])*IN[11]+(T[2]*T[3]+T[5]*T[0])*(T[8]*T[8])*IN[17]+(T[2]*T[3]+T[5]*T[0])*(T[5]*T[8]+T[8]*T[5])*IN[23]+(T[2]*T[3]+T[5]*T[0])*(T[8]*T[2]+T[2]*T[8])*IN[29]+(T[2]*T[3]+T[5]*T[0])*(T[2]*T[5]+T[5]*T[2])*IN[35]);
S[22]=((T[2]*T[0])*(T[1]*T[2])*IN[0]+(T[2]*T[0])*(T[4]*T[5])*IN[6]+(T[2]*T[0])*(T[7]*T[8])*IN[12]+(T[2]*T[0])*(T[4]*T[8]+T[7]*T[5])*IN[18]+(T[2]*T[0])*(T[7]*T[2]+T[1]*T[8])*IN[24]+(T[2]*T[0])*(T[1]*T[5]+T[4]*T[2])*IN[30]+(T[5]*T[3])*(T[1]*T[2])*IN[1]+(T[5]*T[3])*(T[4]*T[5])*IN[7]+(T[5]*T[3])*(T[7]*T[8])*IN[13]+(T[5]*T[3])*(T[4]*T[8]+T[7]*T[5])*IN[19]+(T[5]*T[3])*(T[7]*T[2]+T[1]*T[8])*IN[25]+(T[5]*T[3])*(T[1]*T[5]+T[4]*T[2])*IN[31]+(T[8]*T[6])*(T[1]*T[2])*IN[2]+(T[8]*T[6])*(T[4]*T[5])*IN[8]+(T[8]*T[6])*(T[7]*T[8])*IN[14]+(T[8]*T[6])*(T[4]*T[8]+T[7]*T[5])*IN[20]+(T[8]*T[6])*(T[7]*T[2]+T[1]*T[8])*IN[26]+(T[8]*T[6])*(T[1]*T[5]+T[4]*T[2])*IN[32]+(T[5]*T[6]+T[8]*T[3])*(T[1]*T[2])*IN[3]+(T[5]*T[6]+T[8]*T[3])*(T[4]*T[5])*IN[9]+(T[5]*T[6]+T[8]*T[3])*(T[7]*T[8])*IN[15]+(T[5]*T[6]+T[8]*T[3])*(T[4]*T[8]+T[7]*T[5])*IN[21]+(T[5]*T[6]+T[8]*T[3])*(T[7]*T[2]+T[1]*T[8])*IN[27]+(T[5]*T[6]+T[8]*T[3])*(T[1]*T[5]+T[4]*T[2])*IN[33]+(T[8]*T[0]+T[2]*T[6])*(T[1]*T[2])*IN[4]+(T[8]*T[0]+T[2]*T[6])*(T[4]*T[5])*IN[10]+(T[8]*T[0]+T[2]*T[6])*(T[7]*T[8])*IN[16]+(T[8]*T[0]+T[2]*T[6])*(T[4]*T[8]+T[7]*T[5])*IN[22]+(T[8]*T[0]+T[2]*T[6])*(T[7]*T[2]+T[1]*T[8])*IN[28]+(T[8]*T[0]+T[2]*T[6])*(T[1]*T[5]+T[4]*T[2])*IN[34]+(T[2]*T[3]+T[5]*T[0])*(T[1]*T[2])*IN[5]+(T[2]*T[3]+T[5]*T[0])*(T[4]*T[5])*IN[11]+(T[2]*T[3]+T[5]*T[0])*(T[7]*T[8])*IN[17]+(T[2]*T[3]+T[5]*T[0])*(T[4]*T[8]+T[7]*T[5])*IN[23]+(T[2]*T[3]+T[5]*T[0])*(T[7]*T[2]+T[1]*T[8])*IN[29]+(T[2]*T[3]+T[5]*T[0])*(T[1]*T[5]+T[4]*T[2])*IN[35]);
S[28]=((T[2]*T[0])*(T[2]*T[0])*IN[0]+(T[2]*T[0])*(T[5]*T[3])*IN[6]+(T[2]*T[0])*(T[8]*T[6])*IN[12]+(T[2]*T[0])*(T[5]*T[6]+T[8]*T[3])*IN[18]+(T[2]*T[0])*(T[8]*T[0]+T[2]*T[6])*IN[24]+(T[2]*T[0])*(T[2]*T[3]+T[5]*T[0])*IN[30]+(T[5]*T[3])*(T[2]*T[0])*IN[1]+(T[5]*T[3])*(T[5]*T[3])*IN[7]+(T[5]*T[3])*(T[8]*T[6])*IN[13]+(T[5]*T[3])*(T[5]*T[6]+T[8]*T[3])*IN[19]+(T[5]*T[3])*(T[8]*T[0]+T[2]*T[6])*IN[25]+(T[5]*T[3])*(T[2]*T[3]+T[5]*T[0])*IN[31]+(T[8]*T[6])*(T[2]*T[0])*IN[2]+(T[8]*T[6])*(T[5]*T[3])*IN[8]+(T[8]*T[6])*(T[8]*T[6])*IN[14]+(T[8]*T[6])*(T[5]*T[6]+T[8]*T[3])*IN[20]+(T[8]*T[6])*(T[8]*T[0]+T[2]*T[6])*IN[26]+(T[8]*T[6])*(T[2]*T[3]+T[5]*T[0])*IN[32]+(T[5]*T[6]+T[8]*T[3])*(T[2]*T[0])*IN[3]+(T[5]*T[6]+T[8]*T[3])*(T[5]*T[3])*IN[9]+(T[5]*T[6]+T[8]*T[3])*(T[8]*T[6])*IN[15]+(T[5]*T[6]+T[8]*T[3])*(T[5]*T[6]+T[8]*T[3])*IN[21]+(T[5]*T[6]+T[8]*T[3])*(T[8]*T[0]+T[2]*T[6])*IN[27]+(T[5]*T[6]+T[8]*T[3])*(T[2]*T[3]+T[5]*T[0])*IN[33]+(T[8]*T[0]+T[2]*T[6])*(T[2]*T[0])*IN[4]+(T[8]*T[0]+T[2]*T[6])*(T[5]*T[3])*IN[10]+(T[8]*T[0]+T[2]*T[6])*(T[8]*T[6])*IN[16]+(T[8]*T[0]+T[2]*T[6])*(T[5]*T[6]+T[8]*T[3])*IN[22]+(T[8]*T[0]+T[2]*T[6])*(T[8]*T[0]+T[2]*T[6])*IN[28]+(T[8]*T[0]+T[2]*T[6])*(T[2]*T[3]+T[5]*T[0])*IN[34]+(T[2]*T[3]+T[5]*T[0])*(T[2]*T[0])*IN[5]+(T[2]*T[3]+T[5]*T[0])*(T[5]*T[3])*IN[11]+(T[2]*T[3]+T[5]*T[0])*(T[8]*T[6])*IN[17]+(T[2]*T[3]+T[5]*T[0])*(T[5]*T[6]+T[8]*T[3])*IN[23]+(T[2]*T[3]+T[5]*T[0])*(T[8]*T[0]+T[2]*T[6])*IN[29]+(T[2]*T[3]+T[5]*T[0])*(T[2]*T[3]+T[5]*T[0])*IN[35]);
S[34]=((T[2]*T[0])*(T[0]*T[1])*IN[0]+(T[2]*T[0])*(T[3]*T[4])*IN[6]+(T[2]*T[0])*(T[6]*T[7])*IN[12]+(T[2]*T[0])*(T[3]*T[7]+T[6]*T[4])*IN[18]+(T[2]*T[0])*(T[6]*T[1]+T[0]*T[7])*IN[24]+(T[2]*T[0])*(T[0]*T[4]+T[3]*T[1])*IN[30]+(T[5]*T[3])*(T[0]*T[1])*IN[1]+(T[5]*T[3])*(T[3]*T[4])*IN[7]+(T[5]*T[3])*(T[6]*T[7])*IN[13]+(T[5]*T[3])*(T[3]*T[7]+T[6]*T[4])*IN[19]+(T[5]*T[3])*(T[6]*T[1]+T[0]*T[7])*IN[25]+(T[5]*T[3])*(T[0]*T[4]+T[3]*T[1])*IN[31]+(T[8]*T[6])*(T[0]*T[1])*IN[2]+(T[8]*T[6])*(T[3]*T[4])*IN[8]+(T[8]*T[6])*(T[6]*T[7])*IN[14]+(T[8]*T[6])*(T[3]*T[7]+T[6]*T[4])*IN[20]+(T[8]*T[6])*(T[6]*T[1]+T[0]*T[7])*IN[26]+(T[8]*T[6])*(T[0]*T[4]+T[3]*T[1])*IN[32]+(T[5]*T[6]+T[8]*T[3])*(T[0]*T[1])*IN[3]+(T[5]*T[6]+T[8]*T[3])*(T[3]*T[4])*IN[9]+(T[5]*T[6]+T[8]*T[3])*(T[6]*T[7])*IN[15]+(T[5]*T[6]+T[8]*T[3])*(T[3]*T[7]+T[6]*T[4])*IN[21]+(T[5]*T[6]+T[8]*T[3])*(T[6]*T[1]+T[0]*T[7])*IN[27]+(T[5]*T[6]+T[8]*T[3])*(T[0]*T[4]+T[3]*T[1])*IN[33]+(T[8]*T[0]+T[2]*T[6])*(T[0]*T[1])*IN[4]+(T[8]*T[0]+T[2]*T[6])*(T[3]*T[4])*IN[10]+(T[8]*T[0]+T[2]*T[6])*(T[6]*T[7])*IN[16]+(T[8]*T[0]+T[2]*T[6])*(T[3]*T[7]+T[6]*T[4])*IN[22]+(T[8]*T[0]+T[2]*T[6])*(T[6]*T[1]+T[0]*T[7])*IN[28]+(T[8]*T[0]+T[2]*T[6])*(T[0]*T[4]+T[3]*T[1])*IN[34]+(T[2]*T[3]+T[5]*T[0])*(T[0]*T[1])*IN[5]+(T[2]*T[3]+T[5]*T[0])*(T[3]*T[4])*IN[11]+(T[2]*T[3]+T[5]*T[0])*(T[6]*T[7])*IN[17]+(T[2]*T[3]+T[5]*T[0])*(T[3]*T[7]+T[6]*T[4])*IN[23]+(T[2]*T[3]+T[5]*T[0])*(T[6]*T[1]+T[0]*T[7])*IN[29]+(T[2]*T[3]+T[5]*T[0])*(T[0]*T[4]+T[3]*T[1])*IN[35]);
S[5]=((T[0]*T[1])*(T[0]*T[0])*IN[0]+(T[0]*T[1])*(T[3]*T[3])*IN[6]+(T[0]*T[1])*(T[6]*T[6])*IN[12]+(T[0]*T[1])*(T[3]*T[6]+T[6]*T[3])*IN[18]+(T[0]*T[1])*(T[6]*T[0]+T[0]*T[6])*IN[24]+(T[0]*T[1])*(T[0]*T[3]+T[3]*T[0])*IN[30]+(T[3]*T[4])*(T[0]*T[0])*IN[1]+(T[3]*T[4])*(T[3]*T[3])*IN[7]+(T[3]*T[4])*(T[6]*T[6])*IN[13]+(T[3]*T[4])*(T[3]*T[6]+T[6]*T[3])*IN[19]+(T[3]*T[4])*(T[6]*T[0]+T[0]*T[6])*IN[25]+(T[3]*T[4])*(T[0]*T[3]+T[3]*T[0])*IN[31]+(T[6]*T[7])*(T[0]*T[0])*IN[2]+(T[6]*T[7])*(T[3]*T[3])*IN[8]+(T[6]*T[7])*(T[6]*T[6])*IN[14]+(T[6]*T[7])*(T[3]*T[6]+T[6]*T[3])*IN[20]+(T[6]*T[7])*(T[6]*T[0]+T[0]*T[6])*IN[26]+(T[6]*T[7])*(T[0]*T[3]+T[3]*T[0])*IN[32]+(T[3]*T[7]+T[6]*T[4])*(T[0]*T[0])*IN[3]+(T[3]*T[7]+T[6]*T[4])*(T[3]*T[3])*IN[9]+(T[3]*T[7]+T[6]*T[4])*(T[6]*T[6])*IN[15]+(T[3]*T[7]+T[6]*T[4])*(T[3]*T[6]+T[6]*T[3])*IN[21]+(T[3]*T[7]+T[6]*T[4])*(T[6]*T[0]+T[0]*T[6])*IN[27]+(T[3]*T[7]+T[6]*T[4])*(T[0]*T[3]+T[3]*T[0])*IN[33]+(T[6]*T[1]+T[0]*T[7])*(T[0]*T[0])*IN[4]+(T[6]*T[1]+T[0]*T[7])*(T[3]*T[3])*IN[10]+(T[6]*T[1]+T[0]*T[7])*(T[6]*T[6])*IN[16]+(T[6]*T[1]+T[0]*T[7])*(T[3]*T[6]+T[6]*T[3])*IN[22]+(T[6]*T[1]+T[0]*T[7])*(T[6]*T[0]+T[0]*T[6])*IN[28]+(T[6]*T[1]+T[0]*T[7])*(T[0]*T[3]+T[3]*T[0])*IN[34]+(T[0]*T[4]+T[3]*T[1])*(T[0]*T[0])*IN[5]+(T[0]*T[4]+T[3]*T[1])*(T[3]*T[3])*IN[11]+(T[0]*T[4]+T[3]*T[1])*(T[6]*T[6])*IN[17]+(T[0]*T[4]+T[3]*T[1])*(T[3]*T[6]+T[6]*T[3])*IN[23]+(T[0]*T[4]+T[3]*T[1])*(T[6]*T[0]+T[0]*T[6])*IN[29]+(T[0]*T[4]+T[3]*T[1])*(T[0]*T[3]+T[3]*T[0])*IN[35]);
S[11]=((T[0]*T[1])*(T[1]*T[1])*IN[0]+(T[0]*T[1])*(T[4]*T[4])*IN[6]+(T[0]*T[1])*(T[7]*T[7])*IN[12]+(T[0]*T[1])*(T[4]*T[7]+T[7]*T[4])*IN[18]+(T[0]*T[1])*(T[7]*T[1]+T[1]*T[7])*IN[24]+(T[0]*T[1])*(T[1]*T[4]+T[4]*T[1])*IN[30]+(T[3]*T[4])*(T[1]*T[1])*IN[1]+(T[3]*T[4])*(T[4]*T[4])*IN[7]+(T[3]*T[4])*(T[7]*T[7])*IN[13]+(T[3]*T[4])*(T[4]*T[7]+T[7]*T[4])*IN[19]+(T[3]*T[4])*(T[7]*T[1]+T[1]*T[7])*IN[25]+(T[3]*T[4])*(T[1]*T[4]+T[4]*T[1])*IN[31]+(T[6]*T[7])*(T[1]*T[1])*IN[2]+(T[6]*T[7])*(T[4]*T[4])*IN[8]+(T[6]*T[7])*(T[7]*T[7])*IN[14]+(T[6]*T[7])*(T[4]*T[7]+T[7]*T[4])*IN[20]+(T[6]*T[7])*(T[7]*T[1]+T[1]*T[7])*IN[26]+(T[6]*T[7])*(T[1]*T[4]+T[4]*T[1])*IN[32]+(T[3]*T[7]+T[6]*T[4])*(T[1]*T[1])*IN[3]+(T[3]*T[7]+T[6]*T[4])*(T[4]*T[4])*IN[9]+(T[3]*T[7]+T[6]*T[4])*(T[7]*T[7])*IN[15]+(T[3]*T[7]+T[6]*T[4])*(T[4]*T[7]+T[7]*T[4])*IN[21]+(T[3]*T[7]+T[6]*T[4])*(T[7]*T[1]+T[1]*T[7])*IN[27]+(T[3]*T[7]+T[6]*T[4])*(T[1]*T[4]+T[4]*T[1])*IN[33]+(T[6]*T[1]+T[0]*T[7])*(T[1]*T[1])*IN[4]+(T[6]*T[1]+T[0]*T[7])*(T[4]*T[4])*IN[10]+(T[6]*T[1]+T[0]*T[7])*(T[7]*T[7])*IN[16]+(T[6]*T[1]+T[0]*T[7])*(T[4]*T[7]+T[7]*T[4])*IN[22]+(T[6]*T[1]+T[0]*T[7])*(T[7]*T[1]+T[1]*T[7])*IN[28]+(T[6]*T[1]+T[0]*T[7])*(T[1]*T[4]+T[4]*T[1])*IN[34]+(T[0]*T[4]+T[3]*T[1])*(T[1]*T[1])*IN[5]+(T[0]*T[4]+T[3]*T[1])*(T[4]*T[4])*IN[11]+(T[0]*T[4]+T[3]*T[1])*(T[7]*T[7])*IN[17]+(T[0]*T[4]+T[3]*T[1])*(T[4]*T[7]+T[7]*T[4])*IN[23]+(T[0]*T[4]+T[3]*T[1])*(T[7]*T[1]+T[1]*T[7])*IN[29]+(T[0]*T[4]+T[3]*T[1])*(T[1]*T[4]+T[4]*T[1])*IN[35]);
S[17]=((T[0]*T[1])*(T[2]*T[2])*IN[0]+(T[0]*T[1])*(T[5]*T[5])*IN[6]+(T[0]*T[1])*(T[8]*T[8])*IN[12]+(T[0]*T[1])*(T[5]*T[8]+T[8]*T[5])*IN[18]+(T[0]*T[1])*(T[8]*T[2]+T[2]*T[8])*IN[24]+(T[0]*T[1])*(T[2]*T[5]+T[5]*T[2])*IN[30]+(T[3]*T[4])*(T[2]*T[2])*IN[1]+(T[3]*T[4])*(T[5]*T[5])*IN[7]+(T[3]*T[4])*(T[8]*T[8])*IN[13]+(T[3]*T[4])*(T[5]*T[8]+T[8]*T[5])*IN[19]+(T[3]*T[4])*(T[8]*T[2]+T[2]*T[8])*IN[25]+(T[3]*T[4])*(T[2]*T[5]+T[5]*T[2])*IN[31]+(T[6]*T[7])*(T[2]*T[2])*IN[2]+(T[6]*T[7])*(T[5]*T[5])*IN[8]+(T[6]*T[7])*(T[8]*T[8])*IN[14]+(T[6]*T[7])*(T[5]*T[8]+T[8]*T[5])*IN[20]+(T[6]*T[7])*(T[8]*T[2]+T[2]*T[8])*IN[26]+(T[6]*T[7])*(T[2]*T[5]+T[5]*T[2])*IN[32]+(T[3]*T[7]+T[6]*T[4])*(T[2]*T[2])*IN[3]+(T[3]*T[7]+T[6]*T[4])*(T[5]*T[5])*IN[9]+(T[3]*T[7]+T[6]*T[4])*(T[8]*T[8])*IN[15]+(T[3]*T[7]+T[6]*T[4])*(T[5]*T[8]+T[8]*T[5])*IN[21]+(T[3]*T[7]+T[6]*T[4])*(T[8]*T[2]+T[2]*T[8])*IN[27]+(T[3]*T[7]+T[6]*T[4])*(T[2]*T[5]+T[5]*T[2])*IN[33]+(T[6]*T[1]+T[0]*T[7])*(T[2]*T[2])*IN[4]+(T[6]*T[1]+T[0]*T[7])*(T[5]*T[5])*IN[10]+(T[6]*T[1]+T[0]*T[7])*(T[8]*T[8])*IN[16]+(T[6]*T[1]+T[0]*T[7])*(T[5]*T[8]+T[8]*T[5])*IN[22]+(T[6]*T[1]+T[0]*T[7])*(T[8]*T[2]+T[2]*T[8])*IN[28]+(T[6]*T[1]+T[0]*T[7])*(T[2]*T[5]+T[5]*T[2])*IN[34]+(T[0]*T[4]+T[3]*T[1])*(T[2]*T[2])*IN[5]+(T[0]*T[4]+T[3]*T[1])*(T[5]*T[5])*IN[11]+(T[0]*T[4]+T[3]*T[1])*(T[8]*T[8])*IN[17]+(T[0]*T[4]+T[3]*T[1])*(T[5]*T[8]+T[8]*T[5])*IN[23]+(T[0]*T[4]+T[3]*T[1])*(T[8]*T[2]+T[2]*T[8])*IN[29]+(T[0]*T[4]+T[3]*T[1])*(T[2]*T[5]+T[5]*T[2])*IN[35]);
S[23]=((T[0]*T[1])*(T[1]*T[2])*IN[0]+(T[0]*T[1])*(T[4]*T[5])*IN[6]+(T[0]*T[1])*(T[7]*T[8])*IN[12]+(T[0]*T[1])*(T[4]*T[8]+T[7]*T[5])*IN[18]+(T[0]*T[1])*(T[7]*T[2]+T[1]*T[8])*IN[24]+(T[0]*T[1])*(T[1]*T[5]+T[4]*T[2])*IN[30]+(T[3]*T[4])*(T[1]*T[2])*IN[1]+(T[3]*T[4])*(T[4]*T[5])*IN[7]+(T[3]*T[4])*(T[7]*T[8])*IN[13]+(T[3]*T[4])*(T[4]*T[8]+T[7]*T[5])*IN[19]+(T[3]*T[4])*(T[7]*T[2]+T[1]*T[8])*IN[25]+(T[3]*T[4])*(T[1]*T[5]+T[4]*T[2])*IN[31]+(T[6]*T[7])*(T[1]*T[2])*IN[2]+(T[6]*T[7])*(T[4]*T[5])*IN[8]+(T[6]*T[7])*(T[7]*T[8])*IN[14]+(T[6]*T[7])*(T[4]*T[8]+T[7]*T[5])*IN[20]+(T[6]*T[7])*(T[7]*T[2]+T[1]*T[8])*IN[26]+(T[6]*T[7])*(T[1]*T[5]+T[4]*T[2])*IN[32]+(T[3]*T[7]+T[6]*T[4])*(T[1]*T[2])*IN[3]+(T[3]*T[7]+T[6]*T[4])*(T[4]*T[5])*IN[9]+(T[3]*T[7]+T[6]*T[4])*(T[7]*T[8])*IN[15]+(T[3]*T[7]+T[6]*T[4])*(T[4]*T[8]+T[7]*T[5])*IN[21]+(T[3]*T[7]+T[6]*T[4])*(T[7]*T[2]+T[1]*T[8])*IN[27]+(T[3]*T[7]+T[6]*T[4])*(T[1]*T[5]+T[4]*T[2])*IN[33]+(T[6]*T[1]+T[0]*T[7])*(T[1]*T[2])*IN[4]+(T[6]*T[1]+T[0]*T[7])*(T[4]*T[5])*IN[10]+(T[6]*T[1]+T[0]*T[7])*(T[7]*T[8])*IN[16]+(T[6]*T[1]+T[0]*T[7])*(T[4]*T[8]+T[7]*T[5])*IN[22]+(T[6]*T[1]+T[0]*T[7])*(T[7]*T[2]+T[1]*T[8])*IN[28]+(T[6]*T[1]+T[0]*T[7])*(T[1]*T[5]+T[4]*T[2])*IN[34]+(T[0]*T[4]+T[3]*T[1])*(T[1]*T[2])*IN[5]+(T[0]*T[4]+T[3]*T[1])*(T[4]*T[5])*IN[11]+(T[0]*T[4]+T[3]*T[1])*(T[7]*T[8])*IN[17]+(T[0]*T[4]+T[3]*T[1])*(T[4]*T[8]+T[7]*T[5])*IN[23]+(T[0]*T[4]+T[3]*T[1])*(T[7]*T[2]+T[1]*T[8])*IN[29]+(T[0]*T[4]+T[3]*T[1])*(T[1]*T[5]+T[4]*T[2])*IN[35]);
S[29]=((T[0]*T[1])*(T[2]*T[0])*IN[0]+(T[0]*T[1])*(T[5]*T[3])*IN[6]+(T[0]*T[1])*(T[8]*T[6])*IN[12]+(T[0]*T[1])*(T[5]*T[6]+T[8]*T[3])*IN[18]+(T[0]*T[1])*(T[8]*T[0]+T[2]*T[6])*IN[24]+(T[0]*T[1])*(T[2]*T[3]+T[5]*T[0])*IN[30]+(T[3]*T[4])*(T[2]*T[0])*IN[1]+(T[3]*T[4])*(T[5]*T[3])*IN[7]+(T[3]*T[4])*(T[8]*T[6])*IN[13]+(T[3]*T[4])*(T[5]*T[6]+T[8]*T[3])*IN[19]+(T[3]*T[4])*(T[8]*T[0]+T[2]*T[6])*IN[25]+(T[3]*T[4])*(T[2]*T[3]+T[5]*T[0])*IN[31]+(T[6]*T[7])*(T[2]*T[0])*IN[2]+(T[6]*T[7])*(T[5]*T[3])*IN[8]+(T[6]*T[7])*(T[8]*T[6])*IN[14]+(T[6]*T[7])*(T[5]*T[6]+T[8]*T[3])*IN[20]+(T[6]*T[7])*(T[8]*T[0]+T[2]*T[6])*IN[26]+(T[6]*T[7])*(T[2]*T[3]+T[5]*T[0])*IN[32]+(T[3]*T[7]+T[6]*T[4])*(T[2]*T[0])*IN[3]+(T[3]*T[7]+T[6]*T[4])*(T[5]*T[3])*IN[9]+(T[3]*T[7]+T[6]*T[4])*(T[8]*T[6])*IN[15]+(T[3]*T[7]+T[6]*T[4])*(T[5]*T[6]+T[8]*T[3])*IN[21]+(T[3]*T[7]+T[6]*T[4])*(T[8]*T[0]+T[2]*T[6])*IN[27]+(T[3]*T[7]+T[6]*T[4])*(T[2]*T[3]+T[5]*T[0])*IN[33]+(T[6]*T[1]+T[0]*T[7])*(T[2]*T[0])*IN[4]+(T[6]*T[1]+T[0]*T[7])*(T[5]*T[3])*IN[10]+(T[6]*T[1]+T[0]*T[7])*(T[8]*T[6])*IN[16]+(T[6]*T[1]+T[0]*T[7])*(T[5]*T[6]+T[8]*T[3])*IN[22]+(T[6]*T[1]+T[0]*T[7])*(T[8]*T[0]+T[2]*T[6])*IN[28]+(T[6]*T[1]+T[0]*T[7])*(T[2]*T[3]+T[5]*T[0])*IN[34]+(T[0]*T[4]+T[3]*T[1])*(T[2]*T[0])*IN[5]+(T[0]*T[4]+T[3]*T[1])*(T[5]*T[3])*IN[11]+(T[0]*T[4]+T[3]*T[1])*(T[8]*T[6])*IN[17]+(T[0]*T[4]+T[3]*T[1])*(T[5]*T[6]+T[8]*T[3])*IN[23]+(T[0]*T[4]+T[3]*T[1])*(T[8]*T[0]+T[2]*T[6])*IN[29]+(T[0]*T[4]+T[3]*T[1])*(T[2]*T[3]+T[5]*T[0])*IN[35]);
S[35]=((T[0]*T[1])*(T[0]*T[1])*IN[0]+(T[0]*T[1])*(T[3]*T[4])*IN[6]+(T[0]*T[1])*(T[6]*T[7])*IN[12]+(T[0]*T[1])*(T[3]*T[7]+T[6]*T[4])*IN[18]+(T[0]*T[1])*(T[6]*T[1]+T[0]*T[7])*IN[24]+(T[0]*T[1])*(T[0]*T[4]+T[3]*T[1])*IN[30]+(T[3]*T[4])*(T[0]*T[1])*IN[1]+(T[3]*T[4])*(T[3]*T[4])*IN[7]+(T[3]*T[4])*(T[6]*T[7])*IN[13]+(T[3]*T[4])*(T[3]*T[7]+T[6]*T[4])*IN[19]+(T[3]*T[4])*(T[6]*T[1]+T[0]*T[7])*IN[25]+(T[3]*T[4])*(T[0]*T[4]+T[3]*T[1])*IN[31]+(T[6]*T[7])*(T[0]*T[1])*IN[2]+(T[6]*T[7])*(T[3]*T[4])*IN[8]+(T[6]*T[7])*(T[6]*T[7])*IN[14]+(T[6]*T[7])*(T[3]*T[7]+T[6]*T[4])*IN[20]+(T[6]*T[7])*(T[6]*T[1]+T[0]*T[7])*IN[26]+(T[6]*T[7])*(T[0]*T[4]+T[3]*T[1])*IN[32]+(T[3]*T[7]+T[6]*T[4])*(T[0]*T[1])*IN[3]+(T[3]*T[7]+T[6]*T[4])*(T[3]*T[4])*IN[9]+(T[3]*T[7]+T[6]*T[4])*(T[6]*T[7])*IN[15]+(T[3]*T[7]+T[6]*T[4])*(T[3]*T[7]+T[6]*T[4])*IN[21]+(T[3]*T[7]+T[6]*T[4])*(T[6]*T[1]+T[0]*T[7])*IN[27]+(T[3]*T[7]+T[6]*T[4])*(T[0]*T[4]+T[3]*T[1])*IN[33]+(T[6]*T[1]+T[0]*T[7])*(T[0]*T[1])*IN[4]+(T[6]*T[1]+T[0]*T[7])*(T[3]*T[4])*IN[10]+(T[6]*T[1]+T[0]*T[7])*(T[6]*T[7])*IN[16]+(T[6]*T[1]+T[0]*T[7])*(T[3]*T[7]+T[6]*T[4])*IN[22]+(T[6]*T[1]+T[0]*T[7])*(T[6]*T[1]+T[0]*T[7])*IN[28]+(T[6]*T[1]+T[0]*T[7])*(T[0]*T[4]+T[3]*T[1])*IN[34]+(T[0]*T[4]+T[3]*T[1])*(T[0]*T[1])*IN[5]+(T[0]*T[4]+T[3]*T[1])*(T[3]*T[4])*IN[11]+(T[0]*T[4]+T[3]*T[1])*(T[6]*T[7])*IN[17]+(T[0]*T[4]+T[3]*T[1])*(T[3]*T[7]+T[6]*T[4])*IN[23]+(T[0]*T[4]+T[3]*T[1])*(T[6]*T[1]+T[0]*T[7])*IN[29]+(T[0]*T[4]+T[3]*T[1])*(T[0]*T[4]+T[3]*T[1])*IN[35]);


}

/*-----------------------------------------------------------------------*/
/* #Inv3  inverse of 3x3 matrix ------------------------------------------- */
void Inv3(double *C, double *CI) {

 double jdet;
 int    j1;
 /*
 cof11 = ys*zt - yt*zs; cof12 = yt*zr - yr*zt; cof13 = yr*zs - ys*zr;
 cof21 = zs*xt - zt*xs; cof22 = zt*xr - zr*xt; cof23 = zr*xs - zs*xr;
 cof31 = xs*yt - xt*ys; cof32 = xt*yr - xr*yt; cof33 = xr*ys - xs*yr;
 jdet  = xr*cof11+xs*cof12+xt*cof13; 
 */
 

 CI[0] = C[4]*C[8] - C[5]*C[7]; CI[3] = C[5]*C[6] - C[3]*C[8]; CI[6] = C[3]*C[7] - C[4]*C[6];
 CI[1] = C[7]*C[2] - C[8]*C[1]; CI[4] = C[8]*C[0] - C[6]*C[2]; CI[7] = C[6]*C[1] - C[7]*C[0];
 CI[2] = C[1]*C[5] - C[2]*C[4]; CI[5] = C[2]*C[3] - C[0]*C[5]; CI[8] = C[0]*C[4] - C[1]*C[3];

 jdet = C[0]*CI[0] + C[1]*CI[3] + C[2]*CI[6];
 jdet = 1./jdet;
 if (!mxIsFinite(jdet)) mexErrMsgTxt("Singular matrix"); 
 for (j1=0;j1<9;j1++) CI[j1]*= jdet;

}
/* #cross ------------------------------------------------------------------*/
void cross(double *u, double *v, double *w){

  w[0] = u[1]*v[2]-u[2]*v[1];
  w[1] = u[2]*v[0]-u[0]*v[2];
  w[2] = u[0]*v[1]-u[1]*v[0];

}

/* #EnHeart ----------------------------------------------------------------*/
void EnHeart(int *integ,double *constit,double *I,double *dWdI,double *d2WdI2) {
  /* %C1=0.3MPa, C2=0.2MPa, K=0.3MPa
constit=[.3e6 .2e6 .3e6]/1e6;  % Cenerg : C1 C2 K

dWdI(1) = constit(1)*I(3)^(-1./3.);
dWdI(2) = constit(2)*I(3)^(-2./3.);
dWdI(3) = - 1./3.* constit(1)*I(1)*I(3)^(-4./3.) ...
          - 2./3.* constit(2)*I(2)*I(3)^(-5./3.) ...
              +   constit(3)*(1-I(3)^(-1/2)) ;

d2WdI2=[0 0 -1./3.*constit(1)*I(3)^(-4./3.) ;
        0 0  -2./3.*constit(2)*I(3)^(-5./3.);
       -1./3.*constit(1)*I(3)^(-4./3.)  -2./3.*constit(2)*I(3)^(-5./3.) ...
        (4/9*constit(1)*I(1)*I(3)^(-7./3.) ...
              + 10./9.* constit(2)*I(2)*I(3)^(-8./3.) ...
              + 1/2 * constit(3)*I(3)^(-3/2))]; */
  /* constit[0]=.3e6/1.e6; 
  constit[1]=.2e6/1.e6; 
   constit[2]=.3e6/1.e6; */
  dWdI[0] = constit[0]*pow(I[2],-1./3.);
  dWdI[1] = constit[1]*pow(I[2],-2./3.);
  dWdI[2] = -1./3.*constit[0]*I[0]*pow(I[2],-4./3.) 
            -2./3.*constit[1]*I[1]*pow(I[2],-5./3.) 
                 + constit[2]*(1.-pow(I[2],-.5));
  d2WdI2[0]=0.; d2WdI2[3]=0.; d2WdI2[6]=-1./3.*constit[0]*pow(I[2],-4./3.); 
  d2WdI2[1]=0.; d2WdI2[4]=0.; d2WdI2[7]=-2./3.*constit[1]*pow(I[2],-5./3.);
  d2WdI2[2]= -1./3.*constit[0]*pow(I[2],-4./3.);
  d2WdI2[5]= -2./3.*constit[1]*pow(I[2],-5./3.);
  d2WdI2[8]=  4./9.*constit[0]*I[0]*pow(I[2],-7./3.)
            + 10./9.*constit[1]*I[1]*pow(I[2],-8./3.)
    + .5*constit[2]*pow(I[2],-3./2.);

  /* mexPrintf("%10.5g %10.5g %10.5g \n%10.5g %10.5g %10.5g \n%10.5g %10.5g %10.5g \n",d2WdI2[0],d2WdI2[1],d2WdI2[2],d2WdI2[3],d2WdI2[4],d2WdI2[5],d2WdI2[6],d2WdI2[7],d2WdI2[8]); */
}
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
mxArray* pre_cvs2 () {
 mxArray *st;
 st= mxCreateString("$Revision: 1.135 $  $Date: 2021/09/15 15:21:26 $");
 return(st);
}

#ifndef NULL
#define NULL 0
#endif

/*-----------------------------------------------------------------------*/
int findPosInSparse( int row, int col, mwIndex *jc, mwIndex *ir )
{
  int ii;


#if 1
  for (ii = (int)jc[col]; ii < ((int)jc[col+1]); ii++) {
    if (ir[ii] == row) {
      return( ii );
    }
  }
#else
  RO[0]=row;
  bsearch(RO,ir+jc[col],jc[col+1]-jc[col],sizeof( int32 ), &compareI32 );
#endif

  return( -1 );
}

/* ------------------------------------------------------------------*/
/*      Starting information for meshGraph DOFGraph ...   */
/* ------------------------------------------------------------------*/

/* #ifdef OSTYPEmexrs6
#else
typedef int32 int;
#endif*/
#define float64 double;
/* Globals. */
#define allocMem( Type, num ) \
 (Type *) mxCalloc( 1, (num) * sizeof( Type ) )
#define freeMem( p ) do {\
  mxFree( p ); p = 0; } while (0)
typedef enum ReturnStatus {
  RET_OK,
  RET_Fail
} ReturnStatus;
#if !defined(Max)
#define Max(a,b) (((a) > (b)) ? (a) : (b))
#endif

int32 compareMWI( const void *a, const void *b )
{
  mwIndex i1, i2;

  i1 = *((mwIndex *) a);
  i2 = *((mwIndex *) b);
  return( (int)(i1 - i2) );
}


int32 compareI32( const void *a, const void *b )
{
  int32 i1, i2;

  i1 = *((int32 *) a);
  i2 = *((int32 *) b);
  return( i1 - i2 );
}

/* #mesh_meshGraph -------  */
#undef __FUNC__
#define __FUNC__ "mesh_meshGraph"
mwIndex mesh_meshGraph( mwSize *p_nnz, mwIndex **p_prow, mwIndex **p_icol,
		      mwSize nNod, int nGr, int mConn,
		      mwIndex *EGroup, mwIndex *nEl,
		      mwIndex *DofPerElt, mwIndex *conn )
{

  int     iir, iic, found, in, ii, ip, ig, iel, iep, ir, ic, nn, np, pr, niecMax, nUnique;
  mwIndex   *niec,*pconn, *eonlist, *nir, *nods, *icol;

  niec = allocMem( mwIndex,nNod+1);
  memset(niec,0,(nNod+1)*sizeof(mwIndex)); /* fills with 0 */
  /* nNod number of DOFs, nGr number of groups, mConn=size(conn,1)
     EGroup=EGroup
  */
  pconn=NULL;
  nn = 0;
  for (ig = 0; ig < nGr; ig++) {
    nn += DofPerElt[ig] * nEl[ig]; /* max matrix size */
     for (iel = EGroup[ig]; iel < (EGroup[ig] + nEl[ig]); iel++) {
         pconn = conn + (mConn * iel);if (DofPerElt[ig]>mConn) DofPerElt[ig]=0;
      for (iep = 0; iep < DofPerElt[ig]; iep++) {
	    niec[1+pconn[iep]]++;/* count repetitions of each dof */
      }
    }
  }
  /* mexPrintf( "0b\n"); */

  niec[0] = 0;  niecMax = 1; /* max repetition of each dof */
  for (in = 0; in <= nNod; in++) { niecMax = Max(niecMax,niec[in]);}

  /* mexPrintf( "%d dofs not used, niecMax %i\n", niec[0],niecMax); */

  for (in = 0; in < nNod; in++) {
    niec[in+1] += niec[in];
  }

  eonlist = allocMem(mwIndex,2*nn);

  nir = allocMem( mwIndex, nNod + 1 );
  memset( nir, 0, (nNod + 1) * sizeof( mwIndex ) );

  /* mexPrintf( "1b\n" ); */

  for (ig = 0; ig < nGr; ig++) {
    for (iel = EGroup[ig]; iel < (EGroup[ig] + nEl[ig]); iel++) {
      pconn = conn + mConn * iel;
      for (iep = 0; iep < DofPerElt[ig]; iep++) {
	np = (int)(pconn[iep]);
        if (np==-1) {  }
	else if (np>=0) {
	  eonlist[2*(niec[np]+nir[np])+0] = iel;
	  eonlist[2*(niec[np]+nir[np])+1] = ig;
      /* mexPrintf( "  %d %d %d %d\n", np, eonlist[2*(niec[np]+nir[np])+0], */
      /* eonlist[2*(niec[np]+nir[np])+1], nir[np] ); */
	  nir[np]++;
	}
      }
    }
  }

  /*mexPrintf( "2==================================================\n" );*/ 
 
  memset( nir, 0, (nNod) * sizeof( mwIndex ) );
  nods = allocMem( mwIndex, niecMax * mConn );

  nn = 0;
  for (in = 0; in < nNod; in++) {
    ii = 0;/* mexPrintf( "in=%d/%d\n",in,nNod); */
    for (ip = niec[in]; ip < niec[in+1]; ip++) {
      iel = eonlist[2*(ip)+0];
      ig = eonlist[2*(ip)+1];
      /* mexPrintf( " %d ddl[%d] %d\n", ip, DofPerElt[ig], iel ); */
      for (iep = 0; iep < DofPerElt[ig]; iep++) {
       np = (int)(conn[mConn*iel+iep]);
       if (np==-1) {  }
       else	if (np >= 0) {nods[ii] = conn[mConn*iel+iep];ii++;}
      }
    }
    qsort( nods, ii, sizeof( mwIndex ), &compareMWI );

    nUnique = 1;
    for (ir = 0; ir < (ii - 1); ir++) {
      if (nods[ir] != nods[ir+1]) {
	nUnique++;
      }
    }
    nn += nUnique;/* mexPrintf( " -> %d\n", nUnique ); */
    nir[in] = nUnique;
  }
  
/* mexPrintf( "3==================================================\n" );*/ 

  *p_nnz = nn;
  *p_prow = niec;
  icol = *p_icol = allocMem( mwIndex, nn );

  niec[0] = 0;
  for (in = 0; in < nNod; in++) {
    niec[in+1] = niec[in] + nir[in]; /*      mexPrintf( " %d\n", niec[in+1] ); */
  }
  memset( nir, 0, (nNod + 1) * sizeof( mwIndex ) );

  /* mexPrintf( "4==================================================\n" ); */

  for (ig = 0; ig < nGr; ig++) {
    for (iel = EGroup[ig]; iel < (EGroup[ig] + nEl[ig]); iel++) {
      pconn = conn + mConn * iel;
      for (ir = 0; ir < DofPerElt[ig]; ir++) {
	iir = (int)(pconn[ir]);
	if (iir < 0) continue;
	pr = niec[iir];
  	/* mexPrintf( " %d %d %d\n", iir, pr, niec[iir+1] - pr ); */
	for (ic = 0; ic < DofPerElt[ig]; ic++) {
	  iic = (int)(pconn[ic]);
	  if (iic < 0) continue;
  	  /* mexPrintf( "   %d %d\n", iic, nir[iir] ); */
	  found = 0;
	  for (ii = pr; ii < (pr + nir[iir]); ii++) {
	    if (icol[ii] == iic) {
	      found = 1;
	      break;
	    }
	  }
  	  /* mexPrintf( "  ? %d\n", found ); */
	  if (!found) {
	    if ((int)(nir[iir]) < (int)(niec[iir+1] - pr)) {
	      icol[pr+nir[iir]] = iic;
	      nir[iir]++;
  	      /* mexPrintf( "  + %d %d\n", nir[iir], niec[iir+1] - pr ); */
	    } else {
	      mexPrintf( "  %d %d\n", nir[iir], niec[iir+1] - pr );
	      mexErrMsgTxt( "ERR_VerificationFail\n" );
	    }
	  }
	}
	qsort( icol + pr, nir[iir], sizeof( mwIndex ), &compareMWI );
      }
    }
  }

/*    mexCallMATLAB( 0, 0, 0, 0, "pause" ); */

  freeMem( nods );
  freeMem( nir );
  freeMem( eonlist );

  return( RET_OK );
}
/*
  #mesh_dofGraph
  mex_mesh( 'mesh_meshGraph', size( m.Node, 1 ), m.eli.nGr, size( m.Elt, 2 ), int( m.eli.EGroup ), int( m.eli.nEl ), int( m.eli.nEP ), int( m.Elt' - 1 ) )
*/

#undef __FUNC__
#define __FUNC__ "mesh_dofGraph"
int32 mesh_dofGraph( mwSize *p_dnnz, mwIndex **p_dprow, mwIndex **p_dicol,
		     int32 nNod, mwIndex *nprow, mwIndex *nicol,
		     int32 nEq, int32 *dpn, int32 *dofOffset, int32 *eq )
{
  int32 in, ir, ic, ird, icd, er, ec;
  mwSize dnnz;
  mwIndex *dir, *dprow, *dicol;

  *p_dprow = dprow = allocMem( mwIndex, nEq + 1 );
  memset( dprow, 0, (nEq + 1) * sizeof(mwIndex) );

/*    mexPrintf( "1\n" ); 
  mexPrintf( "%d %d\n", nNod, nEq );
*/
  dnnz = 0;

  /* Row pointers and total number of nonzeros. */
  for (ir = 0; ir < nNod; ir++) {
/*      mexPrintf( "%d\n", dpn[ir] ); */
    for (ird = 0; ird < dpn[ir]; ird++) {
/*        mexPrintf( "%d\n", dofOffset[ir] ); */
      er = eq[dofOffset[ir]+ird];
/*        mexPrintf( "%d %d %d\n", ir, ird, er ); */
      if (er >= 0) {
	for (ic = nprow[ir]; ic < nprow[ir+1]; ic++) {
	  in = nicol[ic];
	  for (icd = 0; icd < dpn[in]; icd++) {
	    ec = eq[dofOffset[in]+icd];
/*  	    mexPrintf( "  %d %d %d\n", ic, icd, ec ); */
	    if (ec >= 0) dprow[er+1]++;
	  }
	}
	dnnz += dprow[er+1];
      }
    }
  }

/*    mexPrintf( "2\n" ); */

  for (ir = 0; ir < nEq; ir++) {
    dprow[ir+1] += dprow[ir];
  }

  if (dnnz != dprow[nEq]) {
    mexPrintf( "(%d == %d)\n", dnnz, dprow[nEq] );
    mexErrMsgTxt( "ERR_VerificationFail!" );
  }

  *p_dnnz = dnnz;
  *p_dicol = dicol = allocMem( mwIndex , dnnz );

  dir = allocMem( mwIndex, nEq + 1 );
  memset( dir, 0, (nEq + 1) * sizeof( mwIndex ) );
  
  /* Sorted column numbers. */
  for (ir = 0; ir < nNod; ir++) {
    for (ird = 0; ird < dpn[ir]; ird++) {
      er = eq[dofOffset[ir]+ird];
      if (er >= 0) {
	for (ic = nprow[ir]; ic < nprow[ir+1]; ic++) {
	  in = nicol[ic];
	  for (icd = 0; icd < dpn[in]; icd++) {
	    ec = eq[dofOffset[in]+icd];
	    if (ec >= 0) {
	      dicol[dprow[er]+dir[er]] = ec;
	      dir[er]++;
	    }
	  }
	}
	qsort( dicol + dprow[er], dir[er], sizeof( int32 ), &compareI32 );
      }
    }
  }
/*    mexPrintf( "3\n" ); */

#ifdef DEBUG_MESH
  for (er = 0; er < nEq; er++) {
    if (dir[er] != (dprow[er+1] - dprow[er])) {
      mexPrintf( "(%d == %d)\n", dir[er], dprow[er+1] - dprow[er] );
      mexErrMsgTxt( "ERR_VerificationFail!" );
    }
  }
#endif

  freeMem( dir );

  return( RET_OK );
}


/* #AssembleSparse ---------------------------------------------------------*/
void AssembleSparse(mwIndex* ir, mwIndex* jc, double* pr, 
		   int NDDL, int IsSymVal, int IsSymK, int* keind, int* elmap, 
	           mwIndex* vir, mwIndex* vjc, double* val) {
    int         ie, lr, lc, gr, gc, ii, jj, j_row_k1; 

	if (pr==NULL) { return; } /* empty */
    else if (vjc==NULL) { /* full matrix */

      for (lc = 0; lc < NDDL; lc++) {
	gc = keind[lc];
	if (gc == -1) continue;
	for (lr = 0; lr < NDDL; lr++) {
	  gr = keind[lr];
	  if (gr == -1) continue;
	  if ((IsSymK) && (gr < gc)) continue;
	  ii = findPosInSparse( gr, gc, jc, ir );
	  if (elmap!=NULL) {ie = elmap[NDDL*lc+lr] - 1;} 
          else {ie = NDDL*lc+lr ;}
	  if (ii < 0) {
	    if (val[ie]==0) {}
	    else {
              mexPrintf("%i (C-row %i, C-col %i) %f",ie+1,gr,gc,val[ie]);
              mexErrMsgTxt( "Non-existent matrix item!" );
            }
	  } else { 
         #pragma omp atomic
         pr[ii] += val[ie];
      }
	}
      }

    } else { /* sparse matrix */
      if (IsSymVal) {
	mexErrMsgTxt( "Symmetric sparse storage of element matrix"
		      " not supported!" );
      }

      for (lc = 0; lc < NDDL; lc++) {
	gc = keind[lc];
	if (gc == -1) continue;
	for (j_row_k1=vjc[lc]; j_row_k1< vjc[lc+1]; j_row_k1++) {
	  lr=vir[j_row_k1]; gr = keind[lr];
	  if (gr == -1) continue;
	  if ((IsSymK) && (gr < gc)) continue;
	  ii = findPosInSparse( gr, gc, jc, ir );
	  if (ii < 0) {
	    if (val[j_row_k1]!=0.0)
             { mexPrintf("%i (%i,%i) %.15g", 
                  j_row_k1+1,gr+1,gc+1,val[j_row_k1]);
	       mexErrMsgTxt( "Non-existent matrix item!" );}
	  } else {
	    jj = findPosInSparse( lr, lc, vjc, vir );
	    if (jj >= 0) {
         #pragma omp atomic
         pr[ii] += val[jj];
        }
	  }
	}
      }
   }

}

/* #AssembleMatVec ---------------------------------------------------------*/
void AssembleMatVec(struct EltConst* ECp,struct GroupFields GF, int* opt,
        int StrategyType, int Mk, int* CurDofPos,
        mwIndex* o_ir, mwIndex* o_jc, double* o_pr,
        int* elmap, int DofPerElt, int Mdef,int Ndef,double *Ener,
        int* cEGI, int jElt,int Mener) {
	 /* assemble elementwise RHS if needed but not ener - - - - - - - - - - - - - - - - - - */
	int i1,j1,j2,j3,jk;
    double *RHS=GF.RHS;
    
     if ((ECp[0].Be!=NULL) && (GF.RHS!=NULL) && ((opt!=NULL) && (opt[1]!=-1))){
      if (GF.VectMap==NULL) mexErrMsgTxt("EltConst.VectMap missing");
      for (jk=0;jk<Mk;jk++) {
		   #ifdef verbose
		    int iz;
		    iz=CurDofPos[GF.VectMap[jk]-1];
            if (iz<0||iz>Mdef) mexErrMsgTxt("position RHS");
          #endif
	/*		{
          #pragma omp critical
				if (ECp[0].Be[jk]!=0) mexPrintf("%i (%i/%i) %g\n",jElt,jk,Mk,ECp[0].Be[jk]);
            }*/
          #pragma omp atomic
		  RHS[CurDofPos[GF.VectMap[jk]-1]]+=ECp[0].Be[jk];
	  }
	 }
     if (StrategyType<=0) {}
     else if ((opt!=NULL) && (opt[1]!=-1)){ 
       /* assemble current element matrix if appropriate */
       AssembleSparse(o_ir,o_jc,o_pr,opt[0],opt[1],opt[2],
         CurDofPos,elmap,NULL,NULL,ECp[0].ke);
     } else if (opt==NULL) {}
     else {/* energy output  - - - - - - - - - - - - - - - - - - -*/ 
        double* vecte, vtmp;
	int ist, j4;
	#if MatlabVER >= 904
	  if (GF.defi==NULL) {ist=(int)1;} else {ist=(int)2;}
    #else
         ist=(int)1;
    #endif
        vecte=(double*)ofMalloc(ist*DofPerElt*sizeof(double));
        if (GF.VectMap==NULL) mexErrMsgTxt("EltConst.VectMap missing");
        for (j2=0;j2<Ndef;j2++) {
          for (j1=0;j1<DofPerElt;j1++) {
            if (elmap!=NULL) i1=GF.VectMap[j1]-1; else i1=j1;
            i1 = CurDofPos[i1];
	    for (j3=0;j3<ist;j3++) {
             if (i1<0) vecte[ist*j1+j3]=0.; else { vecte[ist*j1+j3]=GF.def[ist*i1+j3+j2*ist*Mdef]; }
	    }
          }
          vtmp=0.;
	  for (j4=0;j4<ist;j4++) {
          for  (j1=0;j1<DofPerElt;j1++) {
	    for (j3=0;j3<DofPerElt;j3++) {
	      vtmp+=vecte[ist*j3+j4]*ECp[0].ke[j3+j1*DofPerElt]*vecte[ist*j1+j4];
	    } } }
	  #if MatlabVER < 904
           if (GF.defi!=NULL) {
            for (j1=0;j1<DofPerElt;j1++) {
             if (elmap!=NULL) i1=GF.VectMap[j1]-1; else i1=j1;
             i1 = CurDofPos[i1];
             if (i1<0) vecte[j1]=0.; else { vecte[j1]=GF.defi[i1+j2*Mdef]; }
            }
            for (j1=0;j1<DofPerElt;j1++) { for(j3=0;j3<DofPerElt;j3++) {
             vtmp+=vecte[j3]*ECp[0].ke[j3+j1*DofPerElt]*vecte[j1];
            }}
          } /* if imaginary part consider it */
	  #endif
          Ener[cEGI[jElt]-1 + j2*Mener]=vtmp;
        } /* j2 */ 
        ofFree(vecte);
     }  /* (opt!=NULL) - - - - - - - - - - - - - - - - - - - - - - - */

}



void threadToEC(struct GroupFields GF, struct EltConst EC,int* ThreadiPos,struct EltConst ECe) {

	  memcpy(EC.nodeE,ECe.nodeE,ThreadiPos[0]*sizeof(double));
	  memcpy(EC.NDN,ECe.NDN,ThreadiPos[1]*sizeof(double));
	  memcpy(EC.jdet,ECe.jdet,ThreadiPos[2]*sizeof(double));
      if (EC.bas!=NULL) memcpy(EC.bas,ECe.bas,ThreadiPos[3]*sizeof(double));
      if (EC.J!=NULL) memcpy(EC.J,ECe.J,ThreadiPos[4]*sizeof(double));
      if (EC.ke!=NULL) memcpy(EC.ke,ECe.ke,ThreadiPos[5]*sizeof(double));
      if (GF.NBe!=0) memcpy(EC.Be,ECe.Be,GF.NBe*sizeof(double)); 
      if (GF.NdefE!=0) memcpy(EC.defE,ECe.defE,GF.NdefE*sizeof(double)); 

}

/* #buildNDN #NDNSwitch ----------------------------------------------------*/
void NDNSwitch(int type, struct GroupFields GF, struct EltConst* ECp, 
           int Nw, int Nnode, int Nfield)
{
  int j1, j2 ,jw, Nshape;
 double *N,*Nr,*Ns,*Nt;
 
 N=GF.N; Nr=GF.Nr;  Ns=GF.Ns; Nt=GF.Nt; 

 Nshape=GF.Nshape;
    switch (type){
    case 0: { break; }/* do nothing */
/*begin 3D volume integration */
    case 3: {
        double  xr, yr, zr, xs, ys, zs, xt, yt, zt,
	        cof11, cof12, cof13, cof21, cof22, cof23, cof31, cof32, cof33;
        int           jw, j2;
	
	if (Nt==NULL) mexErrMsgTxt("3D integration rules called for 2D element");
    for (jw=0;jw<Nw;jw++) { /* loop on integration points */

          xr=0.;xs=0.;xt=0.;yr=0.;ys=0.;yt=0.;zr=0.;zs=0.;zt=0.;
          for (j2=0;j2<Nnode;j2++) {
            xr += Nr[j2*Nw+jw] * ECp[0].nodeE[j2]; 
            xs += Ns[j2*Nw+jw] * ECp[0].nodeE[j2]; 
            xt += Nt[j2*Nw+jw] * ECp[0].nodeE[j2]; 
            yr += Nr[j2*Nw+jw] * ECp[0].nodeE[j2+Nnode]; 
            ys += Ns[j2*Nw+jw] * ECp[0].nodeE[j2+Nnode]; 
            yt += Nt[j2*Nw+jw] * ECp[0].nodeE[j2+Nnode]; 
            zr += Nr[j2*Nw+jw] * ECp[0].nodeE[j2+2*Nnode]; 
            zs += Ns[j2*Nw+jw] * ECp[0].nodeE[j2+2*Nnode]; 
            zt += Nt[j2*Nw+jw] * ECp[0].nodeE[j2+2*Nnode]; 
          } /* j2 */

         cof11 = ys*zt - yt*zs; cof12 = yt*zr - yr*zt; cof13 = yr*zs - ys*zr;
         cof21 = zs*xt - zt*xs; cof22 = zt*xr - zr*xt; cof23 = zr*xs - zs*xr;
         cof31 = xs*yt - xt*ys; cof32 = xt*yr - xr*yt; cof33 = xr*ys - xs*yr;
		 if (ECp[0].jdet==NULL) mexErrMsgTxt("jdet init problem");
         ECp[0].jdet[jw] = xr*cof11+xs*cof12+xt*cof13; 

         if (ECp[0].J!=NULL) { /* needed for blockstrain, check init in initInfoAtNode */
             /* mexEvalString("dbstack;"); */
         	 ECp[0].J[0 + 9 * jw] = cof11 / ECp[0].jdet[jw];
			 ECp[0].J[1 + 9 * jw] = cof12 / ECp[0].jdet[jw];
			 ECp[0].J[2 + 9 * jw] = cof13 / ECp[0].jdet[jw];
			 ECp[0].J[3 + 9 * jw] = cof21 / ECp[0].jdet[jw];
			 ECp[0].J[4 + 9 * jw] = cof22 / ECp[0].jdet[jw];
			 ECp[0].J[5 + 9 * jw] = cof23 / ECp[0].jdet[jw];
			 ECp[0].J[6 + 9 * jw] = cof31 / ECp[0].jdet[jw];
			 ECp[0].J[7 + 9 * jw] = cof32 / ECp[0].jdet[jw];
			 ECp[0].J[8 + 9 * jw] = cof33 / ECp[0].jdet[jw];
		 } 
		 
         for (j2=0;j2<Nshape;j2++) {
          ECp[0].NDN[j2+jw*Nshape]=N[jw+j2*Nw];                              /*N*/
          ECp[0].NDN[j2+(jw+Nw)*Nshape] = ( cof11*Nr[jw+j2*Nw]
                                  + cof12*Ns[jw+j2*Nw]
                                  + cof13*Nt[jw+j2*Nw] )/ECp[0].jdet[jw];   /*Nx*/
          ECp[0].NDN[j2+(jw+2*Nw)*Nshape] = ( cof21*Nr[jw+j2*Nw]
                                    + cof22*Ns[jw+j2*Nw]
             		            + cof23*Nt[jw+j2*Nw]) /ECp[0].jdet[jw]; /*Ny*/ 
          ECp[0].NDN[j2+(jw+3*Nw)*Nshape] = ( cof31*Nr[jw+j2*Nw]
                                    + cof32*Ns[jw+j2*Nw]
                                    + cof33*Nt[jw+j2*Nw]) /ECp[0].jdet[jw]; /*Nz*/
         } /* j2 NNode*/
     } /* jw Nw */

    break;
    }


/* begin 3D volume integration specific to mecha3DintegH - - - - - */
/* NDN is now NDN(ndim+1,nNode,Nw) to enable block multiplication  */
/* using blas functions in mecha3DintegH                           */
/*      | N_1   ...   N_Nw  |                                      */
/* NDN =| Nx_1  ...   Nx_Nw |                                      */
/*      | Ny_1  ...   Ny_Nw |                                      */
/*      | Nz_1  ...   Nz_Nw |                                      */

    case 31: {
        double        xr, yr, zr, xs, ys, zs, xt, yt, zt,
	              cof11, cof12, cof13, cof21, cof22, cof23, 
                      cof31, cof32, cof33;
        int           jw, j2;
        of_ptrdiff pNnode,pNw, unit=1;

	if (Nt==NULL) mexErrMsgTxt("3D integration rules called for 2D element");
    pNnode=Nnode;pNw=Nw;
    
	for (jw=0;jw<Nw;jw++) { /* loop on integration points */

          xr=0.;xs=0.;xt=0.;yr=0.;ys=0.;yt=0.;zr=0.;zs=0.;zt=0.;
          for (j2=0;j2<Nnode;j2++) {
            xr += Nr[j2*Nw+jw] * ECp[0].nodeE[j2]; 
            xs += Ns[j2*Nw+jw] * ECp[0].nodeE[j2]; 
            xt += Nt[j2*Nw+jw] * ECp[0].nodeE[j2]; 
            yr += Nr[j2*Nw+jw] * ECp[0].nodeE[j2+Nnode]; 
            ys += Ns[j2*Nw+jw] * ECp[0].nodeE[j2+Nnode]; 
            yt += Nt[j2*Nw+jw] * ECp[0].nodeE[j2+Nnode]; 
            zr += Nr[j2*Nw+jw] * ECp[0].nodeE[j2+2*Nnode]; 
            zs += Ns[j2*Nw+jw] * ECp[0].nodeE[j2+2*Nnode]; 
            zt += Nt[j2*Nw+jw] * ECp[0].nodeE[j2+2*Nnode]; 
	  } 

	    /*Riv
	    xr=of_ddot(&pNnode,Nr+jw,&pNw,ECp[0].nodeE,&unit);
	    xs=of_ddot(&pNnode,Ns+jw,&pNw,ECp[0].nodeE,&unit);
	    xt=of_ddot(&pNnode,Nt+jw,&pNw,ECp[0].nodeE,&unit);
	    yr=of_ddot(&pNnode,Nr+jw,&pNw,ECp[0].nodeE+Nnode,&unit);
	    ys=of_ddot(&pNnode,Ns+jw,&pNw,ECp[0].nodeE+Nnode,&unit);
	    yt=of_ddot(&pNnode,Nt+jw,&pNw,ECp[0].nodeE+Nnode,&unit);
	    zr=of_ddot(&pNnode,Nr+jw,&pNw,ECp[0].nodeE+2*Nnode,&unit);
	    zs=of_ddot(&pNnode,Ns+jw,&pNw,ECp[0].nodeE+2*Nnode,&unit);
	    zt=of_ddot(&pNnode,Nt+jw,&pNw,ECp[0].nodeE+2*Nnode,&unit);
            */


         cof11 = ys*zt - yt*zs; cof12 = yt*zr - yr*zt; cof13 = yr*zs - ys*zr;
         cof21 = zs*xt - zt*xs; cof22 = zt*xr - zr*xt; cof23 = zr*xs - zs*xr;
         cof31 = xs*yt - xt*ys; cof32 = xt*yr - xr*yt; cof33 = xr*ys - xs*yr;
         ECp[0].jdet[jw] = xr*cof11+xs*cof12+xt*cof13; 
         for (j2=0;j2<Nshape;j2++) {
          ECp[0].NDN[j2*4+jw*Nshape*4]   =N[jw+j2*Nw];                      /*N*/
          ECp[0].NDN[j2*4+jw*Nshape*4+1] = ( cof11*Nr[jw+j2*Nw]
                                  + cof12*Ns[jw+j2*Nw]
                                  + cof13*Nt[jw+j2*Nw] )/ECp[0].jdet[jw];   /*Nx*/
          ECp[0].NDN[j2*4+jw*Nshape*4+2] = ( cof21*Nr[jw+j2*Nw]
                                    + cof22*Ns[jw+j2*Nw]
             		            + cof23*Nt[jw+j2*Nw]) /ECp[0].jdet[jw]; /*Ny*/ 
          ECp[0].NDN[j2*4+jw*Nshape*4+3] = ( cof31*Nr[jw+j2*Nw]
                                    + cof32*Ns[jw+j2*Nw]
                                    + cof33*Nt[jw+j2*Nw]) /ECp[0].jdet[jw]; /*Nz*/
	 } /* j2 NNode*/
        } /* jw Nw */

    break;

/* end 3D volume integration - - - - - - - - - - - - - - - - */
/* begin 2D volume integration - - - - - - - - - - - - - - - - */
}  case 2: {
        double        xr, yr, xs, ys; 
        int           jw, j2;

    if (ECp[0].v1x) {
        mexPrintf("Nfield=%i,v1x_col=%i\n",Nfield,ECp[0].v1x);mexErrMsgTxt("Error use rule 23 for local basis support");}
    for (jw=0;jw<Nw;jw++) { /* loop on integration points */

          xr=0.;xs=0.;yr=0.;ys=0.;
          for (j2=0;j2<Nnode;j2++) {
            xr += Nr[j2*Nw+jw] * ECp[0].nodeE[j2];        xs += Ns[j2*Nw+jw] * ECp[0].nodeE[j2]; 
            yr += Nr[j2*Nw+jw] * ECp[0].nodeE[j2+Nnode];  ys += Ns[j2*Nw+jw] * ECp[0].nodeE[j2+Nnode]; 
	      } /* j2 */

          ECp[0].jdet[jw] = xr*ys-xs*yr; 

          for (j2=0;j2<Nshape;j2++) {   
            ECp[0].NDN[j2+jw*Nshape]        = N[jw+j2*Nw];                  /* N */
            ECp[0].NDN[j2+(jw+Nw)*Nshape]   = ( ys*Nr[jw+j2*Nw]
			              - yr*Ns[jw+j2*Nw] )/ECp[0].jdet[jw]; /* Nx */
            ECp[0].NDN[j2+(jw+2*Nw)*Nshape] = (- xs*Nr[jw+j2*Nw]
                                        + xr*Ns[jw+j2*Nw] )/ECp[0].jdet[jw];/* Ny */
          } /* j2 */
    }/* jw */
    break;
/* end 2D volume integration - - - - - - - - - - - - - - - - */
/* begin 3D surface integration - - - - - - - - - - - - - - - - */
} case 32: {
  /*
model=femesh('testhexa8');  EC.nodeE=model.Node(:,5:7);
opt=integrules('hexa8',-1);
EC.nodeE(:,5:10)=0; EC.nodeE(:,7)=1;  EC.nodeE(:,8)=1; % xe=z and ye=y
integrules('buildndn',32,opt,EC.nodeE)

   */
        double        xr, yr, zr, xs, ys, zs, xt, yt, zt,
	              cof11, cof12, cof13, cof21, cof22, cof23, 
                      cof31, cof32, cof33;
        double        xe[3], ye[3], bas1[9], cof_bas[9];
        int           jw, j2;
	if (Nt==NULL) mexErrMsgTxt("3D integration rules called for 2D element");

        for (jw=0;jw<Nw;jw++) { /* loop on integration points */

          xr=0.;xs=0.;xt=0.;yr=0.;ys=0.;yt=0.;zr=0.;zs=0.;zt=0.;
          xe[0]=0.;xe[1]=0.;xe[2]=0.;
          ye[0]=0.;ye[1]=0.;ye[2]=0.;
          for (j2=0;j2<Nnode;j2++) {

            ECp[0].NDN[j2+jw*Nnode]=N[jw];
            xr += Nr[j2*Nw+jw] * ECp[0].nodeE[j2]; 
            xs += Ns[j2*Nw+jw] * ECp[0].nodeE[j2]; 
            xt += Nt[j2*Nw+jw] * ECp[0].nodeE[j2]; 
            yr += Nr[j2*Nw+jw] * ECp[0].nodeE[j2+Nnode]; 
            ys += Ns[j2*Nw+jw] * ECp[0].nodeE[j2+Nnode]; 
            yt += Nt[j2*Nw+jw] * ECp[0].nodeE[j2+Nnode]; 
            zr += Nr[j2*Nw+jw] * ECp[0].nodeE[j2+2*Nnode]; 
            zs += Ns[j2*Nw+jw] * ECp[0].nodeE[j2+2*Nnode]; 
            zt += Nt[j2*Nw+jw] * ECp[0].nodeE[j2+2*Nnode]; 
	    /* xe=opt.NDN(:,jw)'*EC.nodeE(:,5:7); ye=opt.NDN(:,jw)'*EC.nodeE(:,8:10); */
            xe[0] += ECp[0].NDN[j2] * ECp[0].nodeE[j2+4*Nnode]; 
            xe[1] += ECp[0].NDN[j2] * ECp[0].nodeE[j2+5*Nnode]; 
            xe[2] += ECp[0].NDN[j2] * ECp[0].nodeE[j2+6*Nnode]; 
            ye[0] += ECp[0].NDN[j2] * ECp[0].nodeE[j2+7*Nnode]; 
            ye[1] += ECp[0].NDN[j2] * ECp[0].nodeE[j2+8*Nnode]; 
            ye[2] += ECp[0].NDN[j2] * ECp[0].nodeE[j2+9*Nnode]; 
	   } /* j2 */
           cof11 = ys*zt - yt*zs; cof12 = yt*zr - yr*zt; cof13 = yr*zs - ys*zr;
           cof21 = zs*xt - zt*xs; cof22 = zt*xr - zr*xt; cof23 = zr*xs - zs*xr;
           cof31 = xs*yt - xt*ys; cof32 = xt*yr - xr*yt; cof33 = xr*ys - xs*yr;
           ECp[0].jdet[jw] = xr*cof11+xs*cof12+xt*cof13; 
           basis(xe,ye,bas1);

	   /*
              Nrst=[opt.Nr(jw,:)' opt.Ns(jw,:)' opt.Nt(jw,:)'];
              Nxyze=Nrst*(cof*basis(xe,ye)/EC.jdet); % dxe/dx^T=basis(xe,ye);
              opt.NDN(:,jw+[Nw 2*Nw 3*Nw])=Nxyze;
	   */
           cof_bas[0]= cof11*bas1[0] + cof12*bas1[1] + cof13*bas1[2];
           cof_bas[1]= cof21*bas1[0] + cof22*bas1[1] + cof23*bas1[2];
           cof_bas[2]= cof31*bas1[0] + cof32*bas1[1] + cof33*bas1[2];
           cof_bas[3]= cof11*bas1[3] + cof12*bas1[4] + cof13*bas1[5];
           cof_bas[4]= cof21*bas1[3] + cof22*bas1[4] + cof23*bas1[5];
           cof_bas[5]= cof31*bas1[3] + cof32*bas1[4] + cof33*bas1[5];
           cof_bas[6]= cof11*bas1[6] + cof12*bas1[7] + cof13*bas1[8];
           cof_bas[7]= cof21*bas1[6] + cof22*bas1[7] + cof23*bas1[8];
           cof_bas[8]= cof31*bas1[6] + cof32*bas1[7] + cof33*bas1[8];

           for (j2=0;j2<Nnode;j2++) {
	     ECp[0].NDN[j2+(jw+Nw)*Nnode]=0.;ECp[0].NDN[j2+(jw+2*Nw)*Nnode]=0.;ECp[0].NDN[j2+(jw+3*Nw)*Nnode]=0.;
	   }
           for (j2=0;j2<Nnode;j2++) {
             ECp[0].NDN[j2+(jw+Nw)*Nnode] += Nr[j2*Nw+jw]*cof_bas[0] + Ns[j2*Nw+jw]*cof_bas[1]
                                                        + Nt[j2*Nw+jw]*cof_bas[2];
             ECp[0].NDN[j2+(jw+2*Nw)*Nnode] += Nr[j2*Nw+jw]*cof_bas[3] + Ns[j2*Nw+jw]*cof_bas[4]
                                                        + Nt[j2*Nw+jw]*cof_bas[5];
             ECp[0].NDN[j2+(jw+3*Nw)*Nnode] += Nr[j2*Nw+jw]*cof_bas[6] + Ns[j2*Nw+jw]*cof_bas[7]
                                                        + Nt[j2*Nw+jw]*cof_bas[8];
	   }
           for (j2=0;j2<Nnode;j2++) {
             ECp[0].NDN[j2+(jw+Nw)*Nnode]  /=ECp[0].jdet[jw];
             ECp[0].NDN[j2+(jw+2*Nw)*Nnode]/=ECp[0].jdet[jw];
             ECp[0].NDN[j2+(jw+3*Nw)*Nnode]/=ECp[0].jdet[jw];
           }

        } /* jw Nw */
	break;

    }
/* begin 3D surface integration - - - - - - - - - - - - - - - - */
    case 23: case 231: {

  double        bas1[9],*x, *y, *z, a, b,c,xr, xs, ys, yr, *J1, Je[4],Jr[4],xe[3];

  /* opt.NDN(:,1:size(opt.N,1))=opt.N'; */
  x=bas1;y=bas1+3;z=bas1+6;
  for (j1=0;j1<Nshape;j1++) for (jw=0;jw<Nw;jw++) ECp[0].NDN[j1+Nshape*jw]=N[jw+Nw*j1];
  for (jw=0;jw<Nw;jw++) { /* loop on integration points */

    for (j1=0;j1<3;j1++) x[j1] = 0.; 
    for (j1=0;j1<3;j1++) for (j2=0;j2<Nnode;j2++) x[j1] += Nr[jw+j2*Nw]*ECp[0].nodeE[j2+j1*Nnode];
    xr = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    for (j1=0;j1<3;j1++) x[j1] /=xr;    
    xs = 0.;
    
    for (j1=0;j1<3;j1++) y[j1] = 0.; 
    for (j1=0;j1<3;j1++) for (j2=0;j2<Nnode;j2++) {y[j1] += Ns[jw+j2*Nw]*ECp[0].nodeE[j2+j1*Nnode];}

    a = 0.;    for (j1=0;j1<3;j1++) a -= x[j1] * y[j1];
    for (j1=0;j1<3;j1++) y[j1] += a * x[j1];
    ys = sqrt(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]);  yr=-a;
    for (j1=0;j1<3;j1++) y[j1] /=ys;
    /*    z=-[y(2)*x(3)-y(3)*x(2);y(3)*x(1)-y(1)*x(3);y(1)*x(2)-y(2)*x(1)]    */
    z[0] = -y[1]*x[2] + y[2]*x[1];
    z[1] = -y[2]*x[0] + y[0]*x[2];
    z[2] = -y[0]*x[1] + y[1]*x[0];
    ECp[0].jdet[jw] = xr * ys;
 #if 0 /* earlier error  was compatible with t_rigid('offset') */
    Je[0]= ys/ECp[0].jdet[jw]; Je[2]=-xs/ECp[0].jdet[jw];
    Je[1]=-yr/ECp[0].jdet[jw]; Je[3]= xr/ECp[0].jdet[jw]; 
 #else   /* sdt 6.3 */
    Je[0]= ys/ECp[0].jdet[jw]; Je[1]=-xs/ECp[0].jdet[jw]; 
    Je[2]=-yr/ECp[0].jdet[jw]; Je[3]= xr/ECp[0].jdet[jw]; 
 #endif 
 /*mexPrintf("jxe %i %.3e %.3e %.3e %.3e\n",jxe,Je[0],Je[1],Je[2],Je[3]);*/
 
    if (ECp[0].v1x) { /* a local xe is specified */
      for (j1=0;j1<3;j1++) xe[j1] = 0.; 
      for (j1=0;j1<3;j1++) for (j2=0;j2<Nnode;j2++) {xe[j1]+=N[jw+j2*Nw]*ECp[0].nodeE[j2+(j1+ECp[0].v1x)*Nnode];}
      b=xe[0]*x[0]+xe[1]*x[1]+xe[2]*x[2];
      c=xe[0]*y[0]+xe[1]*y[1]+xe[2]*y[2];
      a=sqrt(b*b+c*c);
      if (a>1e-10) {b/=a;c/=a;} else {b=1;c=0;}/* in-plane rotation is [b -c;c b] */
      /* mexPrintf("[%g %g;%g %g]",Je[0],Je[2],Je[1],Je[3]);
       mexPrintf("*[%g %g;%g %g]\n",b,-c,c,b); */
      Jr[0]=Je[0]*b+Je[2]*c; Jr[2]=-Je[0]*c+Je[2]*b;
      Jr[1]=Je[1]*b+Je[3]*c; Jr[3]=-Je[1]*c+Je[3]*b; 
       /*mexPrintf("[%g %g;%g %g]'*[%g %g;%g %g]\n",b,-c,c,b,Jr[0],Jr[1],Jr[2],Jr[3]);*/
      /*
       tr={'b','c','-c','b'}; out=cell(4,1);
     for jjk=0:3; jj=fix(jjk/2); jk=rem(jjk,2); 
      for ji=0:1
      for jl=0:1
         out{ji+2*jl+1}=[out{ji+2*jl+1} sprintf('%s*%s*Je[%i]+',tr{jj+2*ji+1},tr{jk+jl*2+1},jj+jk*2)]
      end
      end  
     end */
       /*
      Jr[0]= b*b*Je[0]+b*c*Je[2]+c*b*Je[1]+c*c*Je[3];
      Jr[1]=-c*b*Je[0]-c*c*Je[2]+b*b*Je[1]+b*c*Je[3];
      Jr[2]=-b*c*Je[0]+b*b*Je[2]-c*c*Je[1]+c*b*Je[3];
      Jr[3]= c*c*Je[0]-c*b*Je[2]-b*c*Je[1]+b*b*Je[3];*/
      if (ECp[0].bas!=NULL) { 
       xe[0]=x[0]*b+y[0]*c; xe[1]=x[1]*b+y[1]*c; xe[2]=x[2]*b+y[2]*c; 
       y[0]=-x[0]*c+y[0]*b; y[1]=-x[1]*c+y[1]*b; y[2]=-x[2]*c+y[2]*b;
       x[0]=xe[0];x[1]=xe[1];x[2]=xe[2];
       }
       basis_clean(bas1);
     } else { /* standard approach where xe is along Nr */ 
       Jr[0]=Je[0]; Jr[1]=Je[1]; Jr[2]=Je[2]; Jr[3]=Je[3]; 
       basis_clean(bas1);
     }
     if (type==231) {xs=Jr[2];Jr[2]=Jr[1];Jr[1]=xs; }/* Transpose J*/
     /* opt.NDN(:,nw+jw)=opt.Nr(jw,:)'/xr; % N,x(jw)
       opt.NDN(:,2*nw+jw)=opt.Ns(jw,:)'/ys-opt.Nr(jw,:)'*(yr/xr/ys); % N,x(jw) */    
     for (j1=0;j1<Nshape;j1++) ECp[0].NDN[ j1+Nshape*(Nw+jw) ] = Nr[jw+Nw*j1]*Jr[0]+Ns[jw+Nw*j1]*Jr[2];
     for (j1=0;j1<Nshape;j1++) ECp[0].NDN[j1+Nshape*(2*Nw+jw)] = Nr[jw+Nw*j1]*Jr[1]+Ns[jw+Nw*j1]*Jr[3];

    /* opt.bas(:,jw)=[x';y';z]; */
    if (ECp[0].J!=NULL) { 
     J1=ECp[0].J+4*jw; J1[0]=Jr[0];J1[1]=Jr[1];J1[2]=Jr[2];J1[3]=Jr[3];
     }
    if (ECp[0].bas!=NULL) {  for (j1=0;j1<9;j1++) ECp[0].bas[j1+9*jw]=bas1[j1];  }
  }   /* for jW */
    break;
/* end 3D surface integration - - - - - - - - - - - - - - - - */
/* #rule13 1D line - - - - - - - - - - - - - - - - - - - - - - - - - -*/
} case 13: {
 double        *r, xr, yr, zr;
 int           jNode;
 if (Nshape>Nnode) {
  r=GF.w-3*Nw; 
  if (Nnode!=2) mexErrMsgTxt("Quadratic rule only implemented for beam1");
  if (GF.w==NULL) mexErrMsgTxt("Error GF.w is needed");
  for (jw=0;jw<Nw;jw++) { /* loop on integration points */
      xr=ECp[0].nodeE[1]-ECp[0].nodeE[0]; yr=ECp[0].nodeE[3]-ECp[0].nodeE[2]; zr=ECp[0].nodeE[5]-ECp[0].nodeE[4]; 
      ECp[0].jdet[jw] = sqrt(xr*xr+yr*yr+zr*zr); 
      /* mexPrintf("%g %g %g %p\n",ECp[0].jdet[jw],r[jw],GF.w[jw],r);*/
      ECp[0].NDN[0+jw*Nshape]=1-r[jw];          ECp[0].NDN[1+jw*Nshape]=r[jw]; /* N */
      ECp[0].NDN[0+(Nw+jw)*Nshape]=-1/ECp[0].jdet[jw]; ECp[0].NDN[1+(Nw+jw)*Nshape]=1/ECp[0].jdet[jw]; /* Nx */
   } /* jw Nw */
 } else {
  for (jw=0;jw<Nw;jw++) { /* loop on integration points */
       xr=0.;yr=0.;zr=0.;
      for (jNode=0;jNode<Nnode;jNode++) {
          xr += Nr[jNode*Nw+jw] * ECp[0].nodeE[jNode]; 
          yr += Nr[jNode*Nw+jw] * ECp[0].nodeE[jNode+Nnode]; 
          zr += Nr[jNode*Nw+jw] * ECp[0].nodeE[jNode+2*Nnode]; 
	  } /* jNode */
      ECp[0].jdet[jw] = sqrt(xr*xr+yr*yr+zr*zr); 
      for (j2=0;j2<Nshape;j2++) {   
            ECp[0].NDN[j2+jw*Nshape]        = N[jw+j2*Nw]; /* N */
            ECp[0].NDN[j2+(jw+Nw)*Nshape]   = Nr[jw+j2*Nw]/ECp[0].jdet[jw]; /* Nx */
          } /* j2 */
   } /* jw Nw */
 }
   break;
/* 1D line - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/* #rule12 line in 2D - - - - - - - - - - - - - - - - - - - - - - - - - -*/
} case 12: {
      double        xr, yr;
      int           jNode;

  for (j1=0;j1<Nshape;j1++) for (jw=0;jw<Nw;jw++) ECp[0].NDN[j1+Nshape*jw]=N[jw+Nw*j1];
        for (jw=0;jw<Nw;jw++) { /* loop on integration points */

          xr=0.;yr=0.;
          for (jNode=0;jNode<Nnode;jNode++) {
            xr += Nr[jNode*Nw+jw] * ECp[0].nodeE[jNode]; 
            yr += Nr[jNode*Nw+jw] * ECp[0].nodeE[jNode+Nnode]; 
	  } /* jNode */

         ECp[0].jdet[jw] = sqrt(xr*xr+yr*yr); 
         if (ECp[0].bas!=NULL) { 
          xr=xr/ECp[0].jdet[jw]; yr=yr/ECp[0].jdet[jw]; 
          ECp[0].bas[4*jw]=xr; ECp[0].bas[4*jw+1]=yr; ECp[0].bas[4*jw+2]=-yr; ECp[0].bas[4*jw+3]=xr; 
         }

        } /* jw Nw */
    break;
/* 1D line - - - - - - - - - - - - - - - - - - - - - - - - - -*/
      } default:
       mexPrintf("type=%i",type);
       mexErrMsgTxt( "Not a supported NDN method" );
  }
}



/* #Mecha3DInteg -------------------------------------------------------------
function Mecha3DInteg(k,Be,DD,w,jdet,NDN,Nnode,Ndof,Nw,jW)
   DD=zeros(9); 
   for jn=0:2; for jj=0:2; for jm=0:2; for jl=0:2; 
     for ji=0:2; for jk=0:2
      DD(jn+3*jj+9*jm+27*jl+1)=DD(jn+3*jj+9*jm+27*jl+1)+ ...
       F_ij(jn+ji*3+1)* ...
       d2wde2(ind_ts_eg(ji+3*jj+1)+6*(ind_ts_eg(jk+3*jl+1)-1))* ...
       F_ij(jm+jk*3+1);
     end;end
     if jn==jm; 
       DD(jn+3*jj+9*jm+27*jl+1)=DD(jn+3*jj+9*jm+27*jl+1)+ ...
        Sigma(jj+3*jl+1);
     end
   end;end;end;end

   for ji=0:2;for jj=0:2;for jk=0:2;for jl=0:2;
    coef=DD(ji+jj*3+jk*9+jl*27+1)*jdet(jW+1)*w(jW+1);
    of_mk('k<-k+a*x*y',k,NDN,NDN,coef, ...
     int32([Nnode*ji+Ndof*Nnode*jk ... % block in stiffness matrix
            Nnode*Nw*[1+jj 1+jl]   ... % columns in NDN
            Ndof Nnode Nnode 1 1]));
   end;end;end;end
*/

void Mecha3DInteg( struct EltConst* ECp, double* F, double* d2wde2, double* Sigma, 
             double* w,int Nnode,int Ndof,int Nw, int jw) {
  double coef[1], DD[81], val;
  int    ji,jj,jk,jl,jm,jn,ind_ts_eg[9]={0,5,4,5,1,3,4,3,2};
  of_ptrdiff pNnode=Nnode, unit=1, pNdof=Ndof;
  
  memset(DD,0,81*sizeof(double)); 
  for (jn=0;jn<3;jn++) {
  for (jj=0;jj<3;jj++) {
  for (jm=0;jm<3;jm++) {
  for (jl=0;jl<3;jl++) {
  for (ji=0;ji<3;ji++) {
  for (jk=0;jk<3;jk++) {
      DD[jn+3*jj+9*jm+27*jl]+=
       F[jn+ji*3]*
       d2wde2[ind_ts_eg[ji+3*jj]+6*ind_ts_eg[jk+3*jl]]*
       F[jm+jk*3];
  }}
  if (jn==jm) {
       DD[jn+3*jj+9*jm+27*jl]+=Sigma[jj+3*jl];
  }
  }}}}

/*
#pragma omp critical
		{
			int it;
		 it=omp_get_thread_num();if (it==0) mexPrintf(" (%i) -> %p\n",it,EC.NDN);
		}
		*/
  /*for (ji=0;ji<9;ji++) { for (jj=0;jj<9;jj++) {
     mexPrintf("%10.5g",DD[ji+9*jj]);} mexPrintf("\n");};mexPrintf("\n");
  */
   for (ji=0;ji<3;ji++) {for (jj=0;jj<3;jj++) { for (jk=0;jk<3;jk++) { for (jl=0;jl<3;jl++) {
    coef[0]=DD[ji+jj*3+jk*9+jl*27]*ECp[0].jdet[jw]*w[jw];
   /* opt=[offk[0] offx[1] offy[2] size(k,1)[3] m[4] n[5]  incx[6] incy[7]] */
    of_dger(&pNnode,&pNnode,coef,ECp[0].NDN+Nnode*(Nw*(jj+1)+jw),&unit,
      ECp[0].NDN+Nnode*(Nw*(jl+1)+jw),&unit,
      ECp[0].ke+Nnode*ji+Ndof*Nnode*jk,&pNdof);
   }}}}
   /* Be{k} = {Ni} sigma_ij F_kj */
   for (ji=0;ji<3;ji++) { 
     for (jk=0;jk<3;jk++) {
     val=(Sigma[ji]*F[jk]+Sigma[ji+3]*F[jk+3]+Sigma[ji+6]*F[jk+6])*ECp[0].jdet[jw]*w[jw];
     for (jl=0;jl<Nnode;jl++) ECp[0].Be[jl+jk*Nnode]+=ECp[0].NDN[jl+Nnode*(Nw*(ji+1)+jw)]*val;
   }}
   /* mexPrintf("Be=\n");for (ji=0;ji<3*Nnode;ji++) mexPrintf("%g\n",Be[ji]);*/

}


/* #Mecha3DIntegH hyperelastic mecha 3D ------------------------------------*/
void Mecha3DIntegH(struct EltConst* ECp, double *F, double *F1, double* d2wde2, 
          double* d2wvde2, double* Sigma,
	      double* w, int Nnode,int Ndof,int Nw, int jw)
{
    double DD[81] ,D[54], D1[54], AUX[54], TEMP[81],*B, alpha, un=1.,zero=0. ;
    /*DD[9][9] D[6][9] D1[6][9] AUX[9][6]*/
    int    ji,jj,jk,jl,jn,Mk;
    int  unit=1,deux=2, trois=3 ,quatre=4, six=6, neuf=9,nul=0;
    char norm = 'N', trans= 'T';

    Mk = 6*3*Nnode;
    B=ECp[0].Be+3*Nnode;
    /*  of_dcopy(&Mk,&zero,&nul,B,&unit); */
    for (ji=0;ji<18*Nnode-1;ji++) B[ji]=0;

    /*B=calloc(18*Nnode,sizeof(double));*/

    for (ji=0;ji<6;ji++) {for (jj=0;jj<9;jj++) {D[9*ji+jj]=0; D1[9*ji+jj]=0;}}


    D[0]=F[0];  D[3]=F[1];  D[6]=F[2]; 
    D[10]=F[3];  D[13]=F[4];  D[16]=F[5]; 
    D[20]=F[6];  D[23]=F[7];  D[26]=F[8];

    D[28]=F[6];  D[29]=F[3];  D[31]=F[7]; 
    D[32]=F[4];  D[34]=F[8];  D[35]=F[5];

    D[36]=F[6];  D[38]=F[0];  D[39]=F[7]; 
    D[41]=F[1];  D[42]=F[8];  D[44]=F[2];

    D[45]=F[3];  D[46]=F[0];  D[48]=F[4]; 
    D[49]=F[1];  D[51]=F[5];  D[52]=F[2];


    D1[0]=F1[0];  D1[3]=F1[1];  D1[6]=F1[2]; 
    D1[10]=F1[3];  D1[13]=F1[4];  D1[16]=F1[5]; 
    D1[20]=F1[6];  D1[23]=F1[7];  D1[26]=F1[8];

    D1[28]=F1[6];  D1[29]=F1[3];  D1[31]=F1[7]; 
    D1[32]=F1[4];  D1[34]=F1[8];  D1[35]=F1[5];

    D1[36]=F1[6];  D1[38]=F1[0];  D1[39]=F1[7]; 
    D1[41]=F1[1];  D1[42]=F1[8];  D1[44]=F1[2];

    D1[45]=F1[3];  D1[46]=F1[0];  D1[48]=F1[4]; 
    D1[49]=F1[1];  D1[51]=F1[5];  D1[52]=F1[2];


    /* d2wde2 for hyperelastic energy
    d2wvde2 for viscoelastic energy */

    /* aux = d^t * d2wde2 + d1^t * d2wvde2*/
    for (ji=0;ji<9;ji++) 
    {
        for (jj=0;jj<6;jj++) 
        {  
            AUX[6*ji+jj] =0;  

            for (jn=0;jn<6;jn++) 
            {  
                AUX[6*ji+jj] +=  (D[9*jn+ji]*d2wde2[6*jj+jn] + 2*D1[9*jn+ji]*d2wvde2[6*jj+jn]);  
            }
        }
    }


    /* dd = aux * d */

    for (ji=0;ji<9;ji++) 
    { for (jj=0;jj<9;jj++) 
       { DD[9*ji+jj] =0;
          for (jn=0;jn<6;jn++) { DD[9*ji+jj] +=  AUX[6*ji+jn]*D[9*jn+jj];}
       }
    }


    /* rajouter k_nl */

    for (jj=0;jj<3;jj++) 
    {
        DD[9*jj+jj]     += Sigma[jj];
        DD[9*(jj+3)+jj+3] += Sigma[jj];
        DD[9*(jj+6)+jj+6] += Sigma[jj];
        for (ji=0;ji<3;ji++) 
        {
            if(ji!=jj) 
            {
                DD[9*ji+jj]     += Sigma[6-jj-ji];
                DD[9*(ji+3)+jj+3] += Sigma[6-jj-ji];
                DD[9*(ji+6)+jj+6] += Sigma[6-jj-ji];
            }
        }
    }
     alpha = ECp[0].jdet[jw]*w[jw]; 

    /*              Multiplication par blocs 

              K_ij        = DP^T  DD_ij      DP 
         (Nnode x Nnode)        (3 x 3)  (3x Nnode)      */

    for (ji=0;ji<3;ji++) 
    {
        for (jj=0;jj<3;jj++) 
        { 

          /* temp = dd_ij*Dp */

            for (jl=0;jl<3;jl++)
            {
                for (jk=0;jk<Nnode;jk++) 
                {
                    TEMP[jk*3+jl] = DD[9*(3*jj)+3*ji+jl]  *ECp[0].NDN[jw*4*Nnode+4*jk+1]+
                                    DD[9*(3*jj+1)+3*ji+jl]*ECp[0].NDN[jw*4*Nnode+4*jk+2]+
                                    DD[9*(3*jj+2)+3*ji+jl]*ECp[0].NDN[jw*4*Nnode+4*jk+3];
                }
            }  
            /* K_ij = K_ij+jdet[jw]*w[jw] Dp^t Temp*/

            for (jl=0;jl<Nnode;jl++)
            {
                for (jk=0;jk<Nnode;jk++) 
                {
                    ECp[0].ke[Ndof*(Nnode*jj+jk)+ji*Nnode+jl] += alpha*(ECp[0].NDN[jw*4*Nnode+4*jl+1]*TEMP[jk*3]+
                                                                ECp[0].NDN[jw*4*Nnode+4*jl+2]*TEMP[jk*3+1]+
                                                                ECp[0].NDN[jw*4*Nnode+4*jl+3]*TEMP[jk*3+2]);
                }
            }  
        }

        /* B = D*Dp (rhs)*/
        for (jl=0;jl<6;jl++)
        {
            for (jk=0;jk<Nnode;jk++) 
            {
                B[6*(Nnode*ji+jk)+jl]=  D[3*ji+9*jl]  *ECp[0].NDN[jw*4*Nnode+4*jk+1]+
                                        D[3*ji+1+9*jl]*ECp[0].NDN[jw*4*Nnode+4*jk+2]+ 
                                        D[3*ji+2+9*jl]*ECp[0].NDN[jw*4*Nnode+4*jk+3];          
            }
        }      

    }  

    /* Be = Sigma_i*B[i} (rhs)*/
    
    for (ji=0;ji<Ndof;ji++) 
    {
        ECp[0].Be[ji] += alpha*( Sigma[0]*B[6*ji]+Sigma[1]*B[6*ji+1]+Sigma[2]*B[6*ji+2]
            +Sigma[3]*B[6*ji+3]+Sigma[4]*B[6*ji+4]+Sigma[5]*B[6*ji+5] );
    }
     
    /*      if (B!=NULL) free(B);*/

}


/*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   GENERIC LINEAR MULTIPHYSIC ELEMENTS
   standard mat_og matrix assembly with MatrixIntegrationRule
 */

void lin_multi(struct GroupFields GF,struct EltConst* ECp, int Nterms,int Mrule,int Ndef,int jElt,
	  int *GSize,int Mk,int *point,int *rrule)
{

  int i1,i2,j2,j3,j4,jrule,jk,Ndof,jw,wstep,NdofPerField,Nnode;
  double coef,r1;
  int NInfoAtNode,MInfoAtNode;

  MInfoAtNode=GSize[3];NInfoAtNode=GSize[4];
#ifdef verbose
 mexPrintf("LinMulti rule_terms %p constit %p Nterms %i\n", 
        GF.rule_terms,EC.constit,Nterms);
#endif
  
  NdofPerField=GSize[1]; Nnode=GSize[0];
 /* {
  #pragma omp critical
				  for (jk=0;jk<Nnode*3;jk++) {
					  if (Be[jk]!=0)
					  {			  mexPrintf("xxx %i (%i/%i) %g\n",jElt,jk,Nnode*3,Be[jk]);}
				  }
	  } */


  for (j2=0;j2<Nterms;j2++) { /* loop on terms of integration rule */  
  wstep=GF.rule_terms[j2+Nterms*5]; 
  for (jw=0;jw<GF.rule_terms[j2+Nterms*6];jw++) {/* loop on integ. points */ 

      /* reinterp at each point stupid but needed j2 */
      if (GF.CTable!=NULL) { constitInterp(GF,ECp,Nnode*jw,Nnode,GSize);} 
      i1=GF.rule_terms[j2+4*Nterms];
      if (i1<0) {coef=-1;} else coef=1;
      i1=abs(i1); i2=jw+GF.rule_terms[j2+Nterms*7];
      coef = coef* GF.w[i2]*ECp[0].jdet[i2] * ECp[0].constit[ jw*wstep+i1 ];
      if (coef!=0) { 
	/*    mexPrintf("consit(%i)*%.2f(%i)*NDN(i,%i)\n", 
	      jw*wstep+abs(i1)+1,w[i2],i2,(1+jw+rule_terms[j2+2*Nterms]));  */
	for (j4=0;j4<NdofPerField;j4++) { 
	  mwSize ioff;
	  ioff = (GF.rule_terms[j2+Nterms]+j4)*Mk + GF.rule_terms[j2]; 
	  for (j3=0;j3<NdofPerField;j3++) {
	    ECp[0].ke[ioff] += coef*ECp[0].NDN[j3+(jw+GF.rule_terms[j2+2*Nterms])*NdofPerField] *
	      ECp[0].NDN[j4+(jw+GF.rule_terms[j2+3*Nterms])*NdofPerField ];
	    ioff++; 
	  } /* j3 */
	} /* j4 */
	
      } /* coef!= 0*/
  } /* jw point */
  } /* terms of constitutive relation */
#ifdef verbose
 mexPrintf("Be %p rrule %p Ndef %i point[4] %i\n",Be,rrule,Ndef,point[4]);
#endif
  if (ECp[0].Be !=NULL && rrule!=NULL && Ndef>1 && 
          (point[4]==1 || point[4]==5 || point[4]>=100) ) { 
#ifdef verbose
 mexPrintf("Nnode %i InfoAtNode %ix%i->%p Be %p",Nnode,
        MInfoAtNode,NInfoAtNode,Be,Be);
#endif
    Ndof=(int)(MInfoAtNode/Nnode); /* 1 = 8/8*/
for (jrule=0;jrule<Mrule;jrule++) {
      switch (rrule[jrule]) {
      case 101 :  /* #rhs_og Volume load defined in InfoAtNode */
	{
	  if (!(NInfoAtNode*MInfoAtNode)) { } else {/* skipping */
       double *FieldAtNode;
       FieldAtNode=ECp[0].nodeE+Nnode*4;
	  for (jw=rrule[7*Mrule+jrule];
	       jw< rrule[7*Mrule+jrule]+rrule[8*Mrule+jrule];jw++) { 
	    r1=0.;
	    /* r1=EltConst.NDN(:,jw+r3(4)+1)'
	     *InfoAtNode(  r3(2)+[0:Nnode-1]*r3(3)+1  ,jElt);  */
	    for (jk=0;jk<Nnode;jk++) {
	      r1 += ECp[0].NDN[jk+Nnode*(jw+rrule[3*Mrule+jrule])] * 
             FieldAtNode[ rrule[1*Mrule+jrule]+jk*rrule[2*Mrule+jrule]]; 
	    }
		if (r1==0) continue;
	    r1*=ECp[0].jdet[jw]*GF.w[jw]; /* r1=r1*jdet(jw+1)*w(jw+1); */
	    /* if r3(7)>=0; r1=r1*EltConst.bas(r3(7)+1,jw+1);end  */
	    if (rrule[6*Mrule+jrule]>=0) mexErrMsgTxt(".bas not implemented"); 
	    /* in1=double(DofPos(EltConst.VectMap(Nnode*r3(5)+[1:Nnode]),jElt));
	       F(in1)=F(in1)+  EltConst.NDN(:,jw+r3(6)+1)*r1; */
	    for (jk=0;jk<Nnode;jk++) { 
	      /*RHS[CurDofPos[ VectMap[ Nnode*(rrule[4*Mrule+jrule])+jk ]-1 ]]*/ 
	      ECp[0].Be[ Nnode*(rrule[4*Mrule+jrule])+jk ]   += ECp[0].NDN[ jk+Nnode*(jw+rrule[5*Mrule+jrule]) ]*r1;
	    }      

	  } /* jw */
/*	  		  { 
#pragma omp critical
				  for (jk=0;jk<Nnode*3;jk++) {
					  if (Be[jk]!=0)
					  {			  mexPrintf("a %i (%i/%i) %g\n",jElt,jk,Nnode*3,Be[jk]);}
				  }
	  }*/
	  } break;
	} /* case */
      default : 
	{ /* mexErrMsgTxt("Not a supported RHS method"); */ }
      } /* switch */
    } /* jrule */
  }

}



/* #nonlin_elas - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 geometric non linear elastic 3D (see elem0.m for m file implementation)
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */
void nonlin_elas(struct GroupFields GF,struct EltConst* ECp, int *GSize,int jElt)
{

  int j2,jw,ji,jj,jk, Nnode, Nw, Nstate,IsVar=0,*topo=GF.topo;
  int ind_ts_eg[9]={0,5,4,5,1,3,4,3,2};
  double U_ij[9], F_ij[9], d2wde2G[36], d2wde2L[36], *d2wde2, 
	       Sigma[9], e_ij[9], Alpha[9], *pAlpha;

  
  Nnode=GSize[0];Nw=GSize[2];Nstate=GSize[5];
  if (GSize[6]<10) { /* thermoelastic isotropic material */
     IsVar=1;for (ji=0;ji<36;ji++) d2wde2G[ji]=0; pAlpha=Alpha;
  } else { /* general anisotropic elastic material */
   IsVar=0;pAlpha=ECp[0].constit+38;
   for (ji=0;ji<36;ji++) { 
    if (topo[ji]==0) d2wde2G[ji]=0; else; d2wde2G[ji]=ECp[0].constit[topo[ji]-1];
   }
  }
  
  for (jw=0;jw<Nw;jw++) {/* loop on integ. points*/

    /* build U and F */
    for (ji=0;ji<9;ji++)  U_ij[ji]=0;
    if (ECp[0].defE!=NULL) { 
    for (ji=0;ji<3;ji++){ for (jj=0;jj<3;jj++) {for (jk=0;jk<Nnode;jk++) {
	 /*U_ij[ji+3*jj]+=def[CurDofPos[ji+3*jk]]*NDN[jk+Nnode*(Nw*(jj+1)+jw)];*/
	 U_ij[ji+3*jj]+=ECp[0].defE[ji+3*jk]*ECp[0].NDN[jk+Nnode*(Nw*(jj+1)+jw)];
	 }}}
	}
    for (ji=0;ji<9;ji++)  { F_ij[ji]=U_ij[ji];e_ij[ji]=0;}
    F_ij[0]+=1;F_ij[4]+=1;F_ij[8]+=1;
    for (ji=0;ji<3;ji++){ for (jj=0;jj<3;jj++) {for (jk=0;jk<3;jk++) {
	  e_ij[ji+3*jj]+=F_ij[jk+3*ji]*F_ij[jk+3*jj]/2;
	 }
	}}
    e_ij[0]-=.5;e_ij[4]-=.5;e_ij[8]-=.5;

    if (GF.CTable!=NULL) { 
     constitInterp(GF,ECp,Nnode*jw,Nnode,GSize);
     if (IsVar) { /* interpolate isotropic constitutive parameters */
      double nu,E,G,r1,r2;     /* table for rho eta E nu G alpha */
      nu=ECp[0].constit[3]; E=ECp[0].constit[2]; 
      G=ECp[0].constit[5]; if (G==0) G=E/2/(1+nu);
      /* mexPrintf("%15.8e %15.8e %15.8e %15.8e\n",E,nu,G,constit[7]);*/
      r1=nu/(1-nu); /* n/(1-n) */ 
      r2=E*(1-nu)/(1+nu)/(1-2*nu); /* E(1-n)/(1+n)(1-2*n) */
      d2wde2G[0]=r2;d2wde2G[7]=r2;d2wde2G[14]=r2;
      r2=r2*r1;d2wde2G[1]=r2;d2wde2G[2]=r2;d2wde2G[6]=r2;
      d2wde2G[8]=r2;d2wde2G[12]=r2;d2wde2G[13]=r2;
      d2wde2G[21]=G;d2wde2G[28]=G;d2wde2G[35]=G;
      /* missing interpolation */
      Alpha[1]=0;Alpha[2]=0;Alpha[3]=0;
      Alpha[5]=0;Alpha[6]=0;Alpha[7]=0;
	  if (GSize[6]>7) {
	   Alpha[0]=ECp[0].constit[7];Alpha[8]=ECp[0].constit[7];Alpha[4]=ECp[0].constit[7];
	  } else { Alpha[0]=0;Alpha[8]=0;Alpha[4]=0;}
	 } else { /* update interpolated law */
      for (ji=0;ji<36;ji++) { 
       if (topo[ji]==0) d2wde2G[ji]=0; else; d2wde2G[ji]=ECp[0].constit[topo[ji]-1];
      }
	 }
	}

    if (ECp[0].v1x) {    /* here possibly do local coordinate transforms */

     double v1[6], *NDNjw, be[9];
     int i2, j3;

     if (ECp[0].bas==NULL) ECp[0].bas=be; 
     i2=Nnode*ECp[0].v1x;NDNjw=ECp[0].NDN+Nnode*jw;     
     for (j3=0;j3<6;j3++) {v1[j3]=0; /* Interpolate v1x ... */
      for (j2=0;j2<Nnode;j2++) {v1[j3]+=NDNjw[j2]*ECp[0].nodeE[j2+i2];}
      i2+=Nnode;
     }
     basis(v1,v1+3,be);
     TransformLambda(d2wde2G,be,d2wde2L); /* build local dd from global */
	 d2wde2=d2wde2L;
	} else { d2wde2=d2wde2G; }

    /* compute Sigma */
	for (ji=0;ji<9;ji++)  { Sigma[ji]=0; }
	if (ECp[0].integ[0] >= 0) { /* with integ[0]==-1 transport no prestress (visco) */
		for (ji = 0; ji < 9; ji++) {
			for (jj = 0; jj < 9; jj++) {
				Sigma[ji] += d2wde2[ind_ts_eg[ji] + 6 * ind_ts_eg[jj]] * e_ij[jj];
			}
		}
	}
    if (ECp[0].gstate!=NULL&&ECp[0].ke!=NULL) { /* allow for non elastic stress */
      double *gstate;
      jj=Nstate*jElt+jw*6;gstate=ECp[0].gstate;
	  for (ji = 0; ji < 9; ji++) { Sigma[ji] += gstate[jj + ind_ts_eg[ji]]; }
    }
    if (ECp[0].t) { /* thermal pre-stress  - - - - - - - - - - -*/
     double *T, dT, eT[9];
     if (IsVar) dT=-ECp[0].constit[8]; else dT=-ECp[0].constit[47];
     T=ECp[0].nodeE+Nnode*4; 
     for (jk=0;jk<Nnode;jk++) {dT+=ECp[0].NDN[Nnode*jw+jk]*T[jk];}
     for (jk=0;jk<9;jk++) {eT[jk]=dT*pAlpha[jk];} /* thermal strain */
     for (ji=0;ji<9;ji++) { for (jj=0;jj<9;jj++) {
       Sigma[ji]-=d2wde2[ind_ts_eg[ji]+6*ind_ts_eg[jj]]*eT[jj];
      }} 
    } /* thermal pre-stress  - - - - - - - - - - -*/
    
    if (ECp[0].ke!=NULL) { /* assemble matrix */
     Mecha3DInteg(ECp,F_ij,d2wde2,Sigma,GF.w,Nnode,3*Nnode,Nw,jw);
    } else if (ECp[0].gstate!=NULL) { /* return Pioloa Kirchhoff stress */
      double *gstate;
      gstate=ECp[0].gstate+6*jw;
      for (ji=0;ji<9;ji++) {gstate[ind_ts_eg[ji]]=Sigma[ji];}
    }
	/* if (ECp[0].ke==NULL) {} else if (mxIsNaN(ECp[0].ke[0])) { mexErrMsgTxt("Non finite matrix");} */
    /* mexPrintf("sig %10.3e %10.3e %10.3e\n",Sigma[0],Sigma[1],Sigma[2]);
     mexPrintf("sig %10.3e %10.3e %10.3e\n",Sigma[3],Sigma[4],Sigma[5]);
     mexPrintf("sig %10.3e %10.3e %10.3e\n",Sigma[6],Sigma[7],Sigma[8]);*/
    /* allow for non elastic stress */
 } /* loop on gauss points */

}


/* --- Follower Pressure implementation, see implementation in elem0.m---*/

void F_pressure(struct GroupFields GF,struct EltConst* ECp, int Nnode,int Nw,int *point)
{

  int ji,jj,jk,jl,j1,jw,jq,jqq;
  double *w,*N,*Nr,*Ns;

  w=GF.w; N=GF.N; Nr=GF.Nr; Ns=GF.Ns; 
  if (point[4]==5) {

    double Axr[9], Axs[9], p, axr[3], axs[3], res[3];
  
    for (j1=0;j1<9;j1++)  {Axr[j1] = 0.; Axs[j1] = 0.; } 
    for (j1=0;j1<3*Nnode;j1++) ECp[0].Be[j1] = 0.;
    for (j1=0;j1<pow(4*Nnode,2);j1++) ECp[0].ke[j1] = 0.;
       
    for (ji=0;ji<3;ji++) { for(jq=0;jq<Nnode;jq++) {
	ECp[0].defE[ji+4*jq] += ECp[0].nodeE[jq+Nnode*ji];
    }}
  
    for (jw=0;jw<Nw;jw++) {/* LOOP ON INTEGRATION POINTS */
      p=0.;
      for (jq=0;jq<Nnode;jq++) p+=N[jw+Nw*jq]*ECp[0].defE[3+4*jq];
      for (j1=0;j1<3;j1++) {axr[j1] = 0; axs[j1] = 0;}
      for (ji=0;ji<3;ji++) { for(jq=0;jq<Nnode;jq++) {
	  axr[ji]+= Nr[jw+Nw*jq]*ECp[0].defE[ji+4*jq];
	  axs[ji]+= Ns[jw+Nw*jq]*ECp[0].defE[ji+4*jq];
	}}
      Axr[5] = axr[0]; Axr[6] = axr[1]; Axr[1] = axr[2];
      Axr[7] = -axr[0]; Axr[2] = -axr[1]; Axr[3] = -axr[2];
      Axs[5] = axs[0]; Axs[6] = axs[1]; Axs[1] = axs[2];
      Axs[7] = -axs[0]; Axs[2] = -axs[1]; Axs[3] = -axs[2]; 
    
      cross(axr,axs,res);
      
      for (jq=0;jq<Nnode;jq++) { for (jk=0;jk<3;jk++) {
	  ECp[0].Be[jq+Nnode*jk] += p*w[jw]*res[jk]*N[jw+Nw*jq];
	  for (jqq=0;jqq<Nnode;jqq++) { for (jl=0;jl<3;jl++) {
	      ECp[0].ke[(jq+Nnode*jk)+4*Nnode*(jqq+Nnode*jl)] +=
		p*w[jw]*N[jw+Nw*jq]*(Axr[3*jl+jk]*Ns[jw+Nw*jqq]
				     - Axs[3*jl+jk]*Nr[jw+Nw*jqq]);
	    }}
	}}
    }/* end of the loop on integration points */

       
  } else {
  
  /* this needs to be compiled
     for jW=0:Nw-1; % HERE IS THE LOOP ON INTEGRATION POINTS 
     % ({u}.{n}) p ds
     for ji=0:2; 
     coef=EltConst.bas(7+ji,jW+1)*EltConst.jdet(jW+1)*EltConst.w(jW+1,4);
     
     if coef; jj=Nnode*ji+[1:Nnode]; jk=Nnode*3+[1:Nnode];
     
     if point(5)==1;
     ke(jj,jk)=ke(jj,jk)+coef*EltConst.N(jW+1,:)'*EltConst.N(jW+1,:);
     elseif point(5)==2
       ke(jk,jj)=ke(jk,jj)+coef*EltConst.N(jW+1,:)'*EltConst.N(jW+1,:);
       end
       
       end
       end
       end
       mexErrMsgTxt(" Linear fs_coupling still in elem0.m");*/
  
    double coef;

    for (j1=0;j1<pow(4*Nnode,2);j1++) ECp[0].ke[j1] = 0.;
  
    for (jw=0;jw<Nw;jw++) {
      for (ji=0;ji<3;ji++) {
	coef=ECp[0].bas[6+ji+jw*9]*ECp[0].jdet[jw]*w[jw];
	if (coef!=0.) {
	  for (jj=0; jj<Nnode;jj++) { for (jk=0; jk<Nnode;jk++){ 
	      if (point[4]==1) {
		ECp[0].ke[Nnode*ji+jj+4*Nnode*(jk+Nnode*3)] -= coef*N[jw+Nw*jj]*N[jw+Nw*jk];
	      } else if (point[4]==2) {
		ECp[0].ke[jk+Nnode*3+4*Nnode*(Nnode*ji+jj)] += coef*N[jw+Nw*jj]*N[jw+Nw*jk];
	      } /* if */
	    }} /* jk, jj */
	} /* if coef */
      
      }} /* ji, jw */
  


  }
}

/* -----------------------------------------------------------------------*/

void PrestressLaw(double *k, double *Be, double *F, double* Sigma,
	      double* w, double* jdet, double* NDN, int Nnode,int Ndof,int Nw, int jw, int isNonLinear)
{
   double DD[81] ,D[54], TEMP[81],*B, alpha, un=1.,zero=0. ;
/*DD[9][9] D[6][9] D1[6][9] AUX[9][6]*/
  int    ji,jj,jk,jl,Mk;
  int  unit=1,deux=2, trois=3 ,quatre=4, six=6, neuf=9,nul=0;
  char norm = 'N', trans= 'T';


  Mk = 6*3*Nnode;
  B=Be+3*Nnode;
  /*of_dcopy(&Mk,&zero,&nul,B,&unit);*/
  /*B=mxCalloc(18*Nnode,sizeof(double));*/
  for (ji=0;ji<18*Nnode-1;ji++) B[ji]=0;	


  
  for (ji=0;ji<6;ji++) {for (jj=0;jj<9;jj++) {D[9*ji+jj]=0;}}


  D[0]=F[0];  D[3]=F[1];  D[6]=F[2]; 
  D[10]=F[3];  D[13]=F[4];  D[16]=F[5]; 
  D[20]=F[6];  D[23]=F[7];  D[26]=F[8];
 
  D[28]=F[6];  D[29]=F[3];  D[31]=F[7]; 
  D[32]=F[4];  D[34]=F[8];  D[35]=F[5];

  D[36]=F[6];  D[38]=F[0];  D[39]=F[7]; 
  D[41]=F[1];  D[42]=F[8];  D[44]=F[2];

  D[45]=F[3];  D[46]=F[0];  D[48]=F[4]; 
  D[49]=F[1];  D[51]=F[5];  D[52]=F[2];

  alpha = jdet[jw]*w[jw];
  
  if (isNonLinear==1) {
      
	/* dd initialisation */
	for (ji=0;ji<9;ji++) {
		for (jj=0;jj<9;jj++) {
			DD[9*ji+jj] =0;
		}
	}
	
	/* rajouter k_nl */	
	for (jj=0;jj<3;jj++) {
		DD[9*jj+jj]     += Sigma[jj];
		DD[9*(jj+3)+jj+3] += Sigma[jj];
		DD[9*(jj+6)+jj+6] += Sigma[jj];
		for (ji=0;ji<3;ji++) {
			if(ji!=jj) {
				DD[9*ji+jj]     += Sigma[6-jj-ji];
				DD[9*(ji+3)+jj+3] += Sigma[6-jj-ji];
				DD[9*(ji+6)+jj+6] += Sigma[6-jj-ji];
			}
		}
	}
	
	/* Multiplication par blocs 
		K_ij        = DP^T  DD_ij      DP 
	   (Nnode x Nnode)        (3 x 3)  (3x Nnode)
	*/	
	for (ji=0;ji<3;ji++) 
	{
		for (jj=0;jj<3;jj++) 
		{ 
	
			/* temp = dd_ij*Dp */			
			for (jl=0;jl<3;jl++)
            	{
                		for (jk=0;jk<Nnode;jk++) 
                		{
				  TEMP[jk*3+jl] = DD[9*(3*jj)+3*ji+jl]  *NDN[jw*4*Nnode+4*jk+1]+
							DD[9*(3*jj+1)+3*ji+jl]*NDN[jw*4*Nnode+4*jk+2]+
							DD[9*(3*jj+2)+3*ji+jl]*NDN[jw*4*Nnode+4*jk+3];
                		}
            	}  
				 
			/* K_ij = K_ij+jdet[jw]*w[jw] Dp^t Temp*/	
			for (jl=0;jl<Nnode;jl++)
            	{
                		for (jk=0;jk<Nnode;jk++) 
                		{
				  k[Ndof*(Nnode*jj+jk)+ji*Nnode+jl] += alpha*(NDN[jw*4*Nnode+4*jl+1]*TEMP[jk*3]+
											    NDN[jw*4*Nnode+4*jl+2]*TEMP[jk*3+1]+
											    NDN[jw*4*Nnode+4*jl+3]*TEMP[jk*3+2]);
               		 }
            	}  
		}
	}  
  }
  
  for (ji=0;ji<3;ji++) {
	/* B = D*Dp (rhs)*/
	for (jl=0;jl<6;jl++)
	{
		for (jk=0;jk<Nnode;jk++) 
		{
		    B[6*(Nnode*ji+jk)+jl]=  D[3*ji+9*jl]  *NDN[jw*4*Nnode+4*jk+1]+
						    D[3*ji+1+9*jl]*NDN[jw*4*Nnode+4*jk+2]+ 
						    D[3*ji+2+9*jl]*NDN[jw*4*Nnode+4*jk+3];          
		}
	}  

  }  
   
  /* Be = Sigma_i*B[i} (rhs)*/
   for (ji=0;ji<Ndof;ji++) 
   {
        Be[ji] += alpha*( Sigma[0]*B[6*ji]+Sigma[1]*B[6*ji+1]+Sigma[2]*B[6*ji+2]
            +Sigma[3]*B[6*ji+3]+Sigma[4]*B[6*ji+4]+Sigma[5]*B[6*ji+5] );
   }
  
  /*if (B!=NULL) mxFree(B); */

}


