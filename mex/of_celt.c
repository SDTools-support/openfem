#include "mex.h"

/*--------------------------------------------------------------------*/
/*  COMPILED ELEMENTS  */
/*--------------------------------------------------------------------*/

#include "sp_util_subs.c"

void bar1 (double *point, double *integ, double *constit, double *node,
           double *elt, double *k, double *m) ;

void beam1 (double *point, double *integ, double *constit, double *node,
            double *elt, int N[], double *k, double *m) ;




/* ------------------------------------------------------------------------ */
/* the mex function          -----------------------------------------------*/
/* ------------------------------------------------------------------------ */
void mexFunction (int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]) {

 char     *buf ;
 double   *point, *integ, *constit, *node, *elt,
          *k, *m;
 int      N[2], Ndof;

 if (nrhs==0)  {
   double *a;
   plhs[0] = mxCreateDoubleMatrix (1, 1, mxREAL);
   a = mxGetPr (plhs[0]);  *a=5.01;  return; 
   }
 else if (nrhs==1) {
    mxArray *st;
    mxArray *rhs[1], *lhs[1];
    int     *dims;
    plhs[0]=mxCreateString("$Revision: 1.4 $  $Date: 2005/03/29 08:47:32 $");
    return;
 }

 if (nrhs!=6)  mexErrMsgTxt ("6 inputs required");

 /* k1=of_celt(ElemF,point,integ,constit,nodeE); */
 buf = mxArrayToString( prhs[0] ); point   = mxGetPr(prhs[1]); 
 integ   = mxGetPr(prhs[2]);       constit = mxGetPr(prhs[3]); 
 node    = mxGetPr(prhs[4]);       elt     = mxGetPr(prhs[5]); 

 if (mxGetN(prhs[4])!=4) { mexErrMsgTxt ("node : 4 column format needed"); }

 if (!strcmp("bar1",buf)) Ndof=6;
 else if (!strcmp("beam1",buf)) Ndof=12;

 plhs[0] = mxCreateDoubleMatrix (Ndof, Ndof, mxREAL);
 if      (nlhs==1 && point[4]==1.) k = mxGetPr (plhs[0]); 
 else if (nlhs==1 && point[4]==2.) m = mxGetPr (plhs[0]); 
 
 if (nlhs==2) {
  plhs[1] = mxCreateDoubleMatrix (Ndof, Ndof, mxREAL);
  m = mxGetPr (plhs[1]); 
 }





/*----------------------------------------------------------- ismex */
  if (!strcmp("ismex",buf))  {

   double *a;
   plhs[0] = mxCreateDoubleMatrix (1, 1, mxREAL);
   a = mxGetPr (plhs[0]);
   *a=1.0;

/*----------------------------------------------------------------- bar1 */
}  else if (!strcmp("bar1",buf)) {

  bar1 (point, integ, constit, node, elt, k, m) ;


/*----------------------------------------------------------------- beam1 */
}  else if (!strcmp("beam1",buf)) {

  N[0]=mxGetN(prhs[5]);
  N[1]=(int)point[4];
  beam1 (point, integ, constit, node, elt, N, k, m) ;

}


} /* end of mexFunction */







/* ------------------------------------------------------------------------ */
void bar1 (double *point, double *integ, double *constit, 
           double *node, double *elt, double *k, double *m) {

 double   x[6], *x1, *x2, x0[3], l, l2, r1, Coo[36], *z, *kt;
 int      j1, j2, j3;

 if (elt[4]==0.) { x0[0]=1.5; x0[1]=1.5; x0[2]=1.5;}
 else {mexErrMsgTxt ("not allowed at the moment"); }

 for (j1=0;j1<6;j1++) x[j1]=node[j1];
 /* if (node[6]==elt[0]) {    for (j1=0;j1<6;j1++) x[j1]=node[j1]; }
    else {   x[0]=node[1]; x[1]=node[0]; x[2]=node[3]; 
          x[3]=node[2]; x[4]=node[5]; x[5]=node[4];             } */
 l2=0.; for (j1=0;j1<5;j1+=2) l2 += (x[j1]-x[j1+1])*(x[j1]-x[j1+1]);
 l=sqrt(l2);
 /* x(1,:) = x(2,:)-x(1,:); x(2,:) = x0-x(2,:); x=basis(x(1,:),x(2,:),1); */
 x1=mxCalloc(3,sizeof(double));  x2=mxCalloc(3,sizeof(double)); 
 x1[0]=x[1]-x[0];  x1[1]=x[3]-x[2];  x1[2]=x[5]-x[4];
 x2[0]=x0[0]-x[1]; x2[1]=x0[1]-x[3]; x2[2]=x0[2]-x[5];
 z=mxCalloc(9,sizeof(double)); 
 basis(x1,x2,z);
 mxFree(x1);  mxFree(x2);

 if (point[4]==0. || point[4]==1.) { /* stiffness */

  for (j1=0;j1<36;j1++) k[j1]=0.;
  r1= constit[(int)point[6]];/*k(ind,ind)=constit(1+point(7))/l*[1 -1;-1 1];*/
  k[0] = r1/l; k[3] =-r1/l;
  k[18]=-r1/l; k[21]= r1/l;
  /*  x=x'; Coo=zeros(6,6); for j1 = 1:3:6; Coo(j1+[0:2],j1+[0:2]) = x; end */
  for (j1=0;j1<36;j1++)  Coo[j1]=0.;
  for (j1=0;j1<3;j1++)   Coo[j1]    = z[3*j1];
  for (j1=0;j1<3;j1++)   Coo[6+j1]  = z[1+3*j1];
  for (j1=0;j1<3;j1++)   Coo[12+j1] = z[2+3*j1];
  for (j1=0;j1<3;j1++)   Coo[21+j1] = z[3*j1];
  for (j1=0;j1<3;j1++)   Coo[27+j1] = z[1+3*j1];
  for (j1=0;j1<3;j1++)   Coo[33+j1] = z[2+3*j1];
  mxFree(z);
  /*   out = Coo'*k*Coo;   */ 
  kt=mxCalloc(36,sizeof(double)); 
  for (j1=0;j1<6;j1++) {
   for (j2=0;j2<6;j2++) {
    r1=0.; for (j3=0;j3<6;j3++)  r1+=k[j1+6*j3]*Coo[6*j2+j3];
    kt[j1+6*j2] = r1;
   }
  }
  for (j1=0;j1<6;j1++) {
   for (j2=0;j2<6;j2++) {
    r1=0.; for (j3=0;j3<6;j3++) r1+=Coo[6*j1+j3]*kt[6*j2+j3];
    k[j1+6*j2] = r1;
   }
  }
  mxFree(kt);

 } /* end stiffness */


 if (point[4]==0. || point[4]==2.) { /* mass */

  for (j1=0;j1<36;j1++) m[j1]=0.;
  r1= constit[2+(int)point[6]]*l/3;
  m[0]=r1;  m[7]=r1;  m[14]=r1;
  m[21]=r1; m[28]=r1; m[35]=r1;
  r1*=.5;
  m[3]=r1;  m[10]=r1; m[17]=r1;
  m[18]=r1; m[25]=r1; m[32]=r1;

 } /* end mass */


} /* end bar1 */



/* ------------------------------------------------------------------------ */
void beam1 (double *point, double *integ, double *constit, 
           double *node, double *elt, int N[], double *k, double *m) {

 double   x[6], x1[3], x2[3], x0[3], l, l2, r1, *z,
          tr[144], pe[6], ie[6], *z1;
 int      j1, j2, off=0;

 for (j1=0;j1<6;j1++) x[j1]=node[j1];
 
 if (N[0]>6) { /* if normal is given */
 
   /* 
  if (any(elt(6:7))|rem(elt(5),1))  % normal given
   x0 = x(1,:)+elt(5:7);
  else
   x0 = find(node(:,4)==elt(1,5));
   if isempty(x0); x0=[1.5 1.5 1.5];
   elseif size(node,2)==7; x0 = node(x0,5:7); 
   else x0=node(x0,1:3);end
  end
  */
 }
 else if (N[0]<5 || elt[6]==0.) { 
   /*  default if nothing  */
   x0[0]=1.5; x0[1]=1.5; x0[2]=1.5;
 }
 else {
   /* 
  x0 = find(node(:,1)==elt(1,5));
  if isempty(x0); x0=[1.5 1.5 1.5];
  else x0 = node(x0,5:7); end
   */
 }
 for (j1=0;j1<144;j1++) tr[j1]=0.;
 for (j1=0;j1<12;j1++) tr[13*j1]=1.; 
 
 /*
 % off-set if any - - - - - - - - - - - - - - - - - - - - - - - - -
tr = eye(12,12);off=0;

if size(elt,2)>12
if any(elt(11:13))
  tr([51 61 38])= elt(11:13);%[3 5;1 6;2 4];ans(:,1)+12*(ans(:,2)-1)
  tr([62 39 49])=-elt(11:13);%[2 6;3 4;1 5];ans(:,1)+12*(ans(:,2)-1);ans'
  off=1;
end
end
if size(elt,2)>15 
if any(elt(14:16))
for (j1=0;j1<8;j1++) mexPrintf("%i  %i \n",j1,(int)point[j1]);
for (j1=0;j1<3;j1++) mexPrintf("elt %i  %i \n",j1,elt[j1]);
mexErrMsgTxt ("not allowed at the moment");

  tr([129 139 116])= elt(14:16);%[3 5;1 6;2 4]+6;ans(:,1)+12*(ans(:,2)-1)
  tr([140 117 127])=-elt(14:16);%[2 6;3 4;1 5]+6;ans(:,1)+12*(ans(:,2)-1);ans'
  off=1;
end
end

% basis function - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 */
 l2=0.; for (j2=0;j2<5;j2+=2) l2 += (x[j2+1]-x[j2])*(x[j2+1]-x[j2]);
 l=sqrt(l2);
 if (l==0.) mexErrMsgTxt ("beam with zero length");

 /*
 x=sp_util('basis',x(2,:)-x(1,:),x0-x(2,:))';
 */


 x1[0]=x[1]-x[0];  x1[1]=x[3]-x[2];  x1[2]=x[5]-x[4];
 x2[0]=x0[0]-x[1]; x2[1]=x0[1]-x[3]; x2[2]=x0[2]-x[5];

 z=mxCalloc(9,sizeof(double)); 
 basis(x1,x2,z);
 mxFree(x1); mxFree(x2);

 /*
pe=constit(point(7)+[1:6]);pe=pe(:)';
ie=constit(point(7)+[7:12]);ie=ie(:)';
 */

 for (j1=0;j1<6;j1++) pe[j1]=constit[(int)point[6]+j1];
 for (j1=0;j1<6;j1++) ie[j1]=constit[(int)point[6]+6+j1]; /* xxx */

 /*
if ~any([0 1 2 3 100]==typ)
 warning(sprintf('Matrix type %i not supported by beam1',typ))
end
Tpin=[];
 */

 z1 = mxCalloc(144,sizeof(double)); 
 for (j1=0;j1<144;j1++) z1[j1]=0.;

 if (point[4]==0. || point[4]==1.) { /* stiffness */

  r1=pe[0]*ie[3]/l; 
  z1[0]=r1;  z1[6]=-r1; z1[72]=-r1; z1[78]=r1;
  r1=pe[3]*ie[0]/l; /* torsionnal stiffness GJ/l */
  z1[39]=r1; z1[45]=-r1; z1[111]=-r1; z1[117]=r1;
  /* ind = [14 18 20 24 62 66 68 72 86 90 92 96 134 138 140 144];  */
  if (ie[4]==0.) { /* bernoulli */
    /*    k(ind) = (pe(1,1)*ie(1,2)/l)* ...
	           [12/l2    6/l   -12/l2   6/l     6/l   4  -6/l  2 ...
	           -12/l2  -6/l    12/l2  -6/l     6/l   2  -6/l  4];   */
    r1=pe[0]*ie[1]/l;
    z1[13] = 12/l2 * r1; z1[17] = 6/l * r1;   z1[19] = -12/l2 * r1;
    z1[23] = 6/l * r1;   z1[61] = 6/l * r1;   z1[65] = 4 * r1;
    z1[67] = -6/l * r1;  z1[71] = 2 * r1;     z1[85] = -12/l2 * r1;
    z1[89] = -6/l * r1;  z1[91] = 12/l2 * r1; z1[95] = -6/l * r1;
    z1[133]= 6/l * r1;   z1[137]= 2 * r1;     z1[139]= -6/l * r1;   
    z1[143]= 4 * r1;
  }
  else { /* timoshenko with reduced integration */
   mexErrMsgTxt ("xxx timoshenko with reduced integration");
  }
  /*
   ind = [3 5 9 11];
    if ie(6)==0
      k(ind,ind) = (pe(1,1)*ie(1,3)/l)* ...
	 [12/l2   -6/l   -12/l2  -6/l;-6/l   4   6/l  2;
	  -12/l2   6/l    12/l2   6/l;-6/l   2   6/l  4];
    else
       k(ind,ind) = ...
         pe(1,1)*ie(1,3)/l*[0 0 0 0;0 1 0 -1;0 0 0 0;0 -1 0 1] + ... % flexion
         pe(4)*ie(4)*ie(6)/l*[1 -l/2 -1 -l/2;-l/2 l2/4 l/2 l2/4;         % shear
                            -1 l/2 1 l/2;-l/2 l2/4 l/2 l2/4];
    end
  */
  if (ie[5]==0.) {
    r1=pe[0]*ie[2]/l;
    z1[24+2]  = 12/l2 * r1; z1[24+4]  = -6/l * r1;   z1[24+8]  = -12/l2 * r1;
    z1[24+10] = -6/l * r1;  z1[48+2]  = -6/l * r1;   z1[48+4]  = 4 * r1;
    z1[48+8]  = 6/l * r1;   z1[48+10] = 2 * r1;      z1[96+2]  = -12/l2 * r1;
    z1[96+4]  = 6/l * r1;   z1[96+8]  = 12/l2 * r1;  z1[96+10] = 6/l * r1;
    z1[120+2] = -6/l * r1;  z1[120+4] = 2 * r1;      z1[120+8] = 6/l * r1;
    z1[120+10]= 4 * r1;
  }
  else {
   mexErrMsgTxt ("xxx timoshenko with reduced integration");
  }
 }/* end stiffness */


 /*
% pin flag handling (condense DOF a put zero stiffness)
   if size(elt,2)>9
    i2 = elt(1,9);i2=[rem(i2,10) rem(fix(i2/10),10) rem(fix(i2/100),10) ...
         rem(fix(i2/1e3),10) rem(fix(i2/1e4),10) rem(fix(i2/1e5),10)];
    i1=i2(find(i2)); 
    i2 = elt(1,10);i2=[rem(i2,10) rem(fix(i2/10),10) rem(fix(i2/100),10) ...
         rem(fix(i2/1e3),10) rem(fix(i2/1e4),10) rem(fix(i2/1e5),10)];
    i1=[i1 i2(find(i2))+6];
    if ~isempty(i1)
      i2=1:12;i2(i1)=0;i2=find(i2);
      Tpin([i2(:);i1(:)],1:length(i2))=[eye(length(i2),length(i2));-k(i1,i1)\k(i1,i2)];
      k1=zeros(12,12); k1(i2,i2)=Tpin'*k*Tpin; k=k1;
    else Tpin=[]; end
   else Tpin=[]; end
   %if size(elt,2)>9 & any(elt(9:10)) warning('pin flags not implemented');end
 */

 if   (N[1]==0 || N[1]==1) {
  x_k_x (z, z1, k, 12);
  mxFree(z1); mxFree(z);
 }

 if (m!=NULL) { /* mass */

   /* assemble mass matrix in m */

 }/* end mass */
 


 /*
   k=sp_util('xkx',x,k);  if off; k = tr'*k*tr;end % coordinate transformation
   out=k; out1=[];
 */

}
