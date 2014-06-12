#ifndef DS_COMMON_H
#define DS_COMMON_H

#include "machine.h" /* C2F */

#ifdef MSCVER
#ifdef COMMON_EXPORT
#define my_export __declspec (dllexport)
#else
#define my_export __declspec (dllimport)
#endif
#else
#define my_export extern
#endif

typedef void (*fobj_filtersd_ptr)(int *, int *, double *, double *, double *, double *, int *);
typedef void (*grad_filtersd_ptr)(int *, int *, double *, double *, double *, int *);

typedef void (*fobj_glcpd_ptr)(int *, double *, double *, double *, int *, char *);
typedef void (*grad_glcpd_ptr)(int *, double *, double *, double *, int *, char *);

//                              n,    m,    x,       al,      f,       fmin,    cstype,bl,      bu,      ws,      lws,  v,       nv,
my_export void C2F(ds_filtersd)(int *,int *,double *,double *,double *,double *,char *,double *,double *,double *,int *,double *,int *,
				int *,int *,int *,int *,int *,int *,double *,double *,double *,int *,int *, int *,int *, fobj_filtersd_ptr,  grad_filtersd_ptr);
//                              maxa, maxla,maxu, maxiu,kmax, maxg, rho,     htol,    rgtol,   maxit,iprint,nout, ifail, functions,          gradients

//                           n,     m,     k,     kmax,  maxg,  a,        la,    x,        bl,       bu,       f,        fmin,     g,         r,        w,        e,
my_export void C2F(ds_glcpd)(int *, int *, int *, int *, int * ,double *, int *, double *, double *, double *, double *, double *, double * , double *, double *, double *,
			     int *, double *, int *, int *, int *, double *, int *, char *, double *, int *, double *, int *, int *, int *, int *,  int *, fobj_glcpd_ptr,  grad_glcpd_ptr);
//                           ls,    alp,      lp,    mlp,   peq,   ws,       lws,   cws,    v,        nv,    rgtol,    mode,  ifail, mxgr,  iprint, nout,  functions,       gradients

//                               n,     m,     k,     kmax,  maxg,  a,        la,    x,        bl,       bu,       f,        fmin,     g,         r,        w,        e,
my_export void C2F(sci_ds_glcpd)(int *, int *, int *, int *, int * ,double *, int *, double *, double *, double *, double *, double *, double * , double *, double *, double *,
				 int *, double *, int *, int *, int *, double *, int *, char *, double *, int *, double *, int *, int *, int *, int *,  int *, fobj_glcpd_ptr,  grad_glcpd_ptr);
//                               ls,    alp,      lp,    mlp,   peq,   ws,       lws,   cws,    v,        nv,    rgtol,    mode,  ifail, mxgr,  iprint, nout,  functions,       gradients

typedef struct
{
  double ainfty;
  double ubd;
  int    mlp;
  int    mxf;
} FILTERSD_defaultc_struct;

my_export FILTERSD_defaultc_struct C2F(defaultc);

// C2F(defaultc).ainfty
// C2F(defaultc).ubd
// C2F(defaultc).mlp
// C2F(defaultc).mxf

typedef struct
{
  int kk;
  int ll;
  int kkk;
  int lll;
  int mxws;
  int mxlws;
} FILTERSD_wsc_struct;

my_export FILTERSD_wsc_struct C2F(wsc);

// C2F(wsc).kk
// C2F(wsc).ll
// C2F(wsc).kkk
// C2F(wsc).lll
// C2F(wsc).mxws
// C2F(wsc).mxlws

typedef struct
{
  double dnorm;
  double h;
  double hJt;
  double hJ;
  int ipeq;
  int k;
  int itn;
  int nft;
  int ngt;
} FILTERSD_statsc_struct;

my_export FILTERSD_statsc_struct C2F(statsc);

// C2F(statsc).dnorm
// C2F(statsc).h
// C2F(statsc).hJt
// C2F(statsc).hJ
// C2F(statsc).ipeq
// C2F(statsc).k
// C2F(statsc).itn
// C2F(statsc).nft
// C2F(statsc).ngt

typedef struct
{
  int mxgr;
} FILTERSD_ngrc_struct;

my_export FILTERSD_ngrc_struct C2F(ngrc);

// C2F(ngrc).mxgr

typedef struct
{
  int mxm1;
} FILTERSD_mxm1c_struct;

my_export FILTERSD_mxm1c_struct C2F(mxm1c);

// C2F(mxm1c).mxm1

typedef struct
{
  double eps;
  double tol;
  double emin;
} FILTERSD_epsc_struct;

my_export FILTERSD_epsc_struct C2F(epsc);

// C2F(epsc).eps
// C2F(epsc).tol
// C2F(epsc).emin

typedef struct
{
  double sgnf;
  int nrep;
  int npiv;
  int nres;
} FILTERSD_repc_struct;

my_export FILTERSD_repc_struct C2F(repc);

// C2F(repc).sgnf
// C2F(repc).nrep
// C2F(repc).npiv
// C2F(repc).nres

typedef struct
{
  double rgnorm;
  double vstep;
  int    iter;
  int    npv;
  int    nfn;
  int    ngr;
} FILTERSD_infoc_struct;

my_export FILTERSD_infoc_struct C2F(infoc);

// C2F(infoc).rgnorm
// C2F(infoc).vstep
// C2F(infoc).iter
// C2F(infoc).npv
// C2F(infoc).nfn
// C2F(infoc).ngr

#endif
