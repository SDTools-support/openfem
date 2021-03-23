#pragma once

/* Etienne Balmes                                                     */
/*  $Revision: 1.29 $  $Date: 2017/05/09 10:55:00 $                    */
/* This is here to allow platform dependent definition of blass calls */

#ifndef int32
#define int32 int
#endif
#ifndef OF_SMALL
#if MatlabVER<73  /* I don't know why mwSize is not defined yet */
#ifndef mwSize
#define mwSize int
#endif
#ifndef mwIndex
#define mwIndex int
#endif
#endif


#if defined(OSTYPEmexhp7)||defined(OSTYPEmexrs6)||defined(OSTYPEmexglx)||defined(OSTYPEmexmac)||defined(OSTYPEmexmaci)||defined(OSTYPEmexmaci64)||defined(OSTYPEmexsol)||defined(OSTYPEmexa64)

#define of_dger   dger_
#define of_dgeev  dgeev_
#define of_daxpy  daxpy_
#define of_ddot   ddot_
#define of_dnrm2  dnrm2_
/*#if !defined(FORMATLAB)*/
#define of_dcopy  dcopy_
/*#endif*/

#define of_dgemm  dgemm_
#define of_dgemv  dgemv_
#define of_dpotrf dpotrf_
#define of_dsyrk  dsyrk_
#define of_dtrsm  dtrsm_

#define of_zgemm  zgemm_
#define of_zpotrf zpotrf_
#define of_zsyrk  zsyrk_
#define of_ztrsm  ztrsm_

#else /* if defined(OSTYPEdll) */

#define of_dger   dger
#define of_dgeev  dgeev
#define of_daxpy  daxpy
#define of_ddot   ddot
#define of_dnrm2  dnrm2
#define of_dcopy  dcopy

#define of_dgemm  dgemm
#define of_dgemv  dgemv
#define of_dpotrf dpotrf
#define of_dsyrk  dsyrk
#define of_dtrsm  dtrsm

#define of_zgemm  zgemm
#define of_zpotrf zpotrf
#define of_zherk  zherk
#define of_ztrsm  ztrsm

#endif


#if MatlabVER>=78  /* I don't know why mwSize is not defined yet */
 #ifndef MKL
#if defined(OSTYPEmexw64)
 #include "crtdefs.h" 
#endif
 #include "blas.h"
 #endif
 #ifndef of_ptrdiff
 #define of_ptrdiff ptrdiff_t
 #endif
#else
#define of_ptrdiff int
#ifndef of_dger
 extern void of_dger(int *,int *, double *, double *, int *, double *, int *, double *,int *);
 extern void of_dgeev(char *, char*, char*, int*, double *, int*, double *, double*, double **, int*, double **, int *,int*);
#endif
#ifndef of_dcopy
 extern void of_dcopy(int *, double *, int *, double *, int *);
#endif
#ifndef of_daxpy
 extern void of_daxpy(int *, double *, double *, int *, double *, int *);
 extern double of_dnrm2(int *, double *, int *);
 extern double of_ddot(int *, double *, int *, double *, int *);  
 extern void of_dgemm(char *, char *, int *,int *,int *, double *,double *,int *, double *, int *,double *,double *, int *);
 extern void of_dgemv(char *, int *, int *, double *,double *,int *, double *, int *,double *,double *, int *);
#endif
#endif

#ifndef sNDN
struct sNDN {
    int Nnode; 
    int Nw; 
    double* N;
    double* Nr;
    double* Ns;
    double* Nt;
};
#endif
#ifdef _WIN32
#define OF_EXPORT __declspec(dllexport)
#else
#define OF_EXPORT
#endif

#endif 
/*end of_small for src files */
