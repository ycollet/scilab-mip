// Fortran rule:
// IMPLICIT: default IMPLICIT REAL(A-H, O-Z), INTEGER(I-N)
// here, we wil have IMPLICIT DOUBLE PRECISION(a-h,o-z), INTEGER(i-n)

#include <stack-c.h>
#include <sciprint.h>
#include <Scierror.h>
#include <MALLOC.h>
#include <api_parameters.h>

#include <api_scilab.h>

#include <stdlib.h>
#include <string.h>

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

#define CONMIN_ERROR if(_SciErr.iErr)		\
    {                                           \
      printError(&_SciErr, 0);                  \
      return _SciErr.iErr;                      \
    }

//#define DEBUG 1
#define DBGPRINTF sciprint
//#define DBGPRINTF printf

#ifdef WIN32
#define my_import __declspec(dllimport)
#define my_export __declspec(dllexport)
#else
#define my_import extern
#define my_export extern
#endif

// Order of the parameters:
#define X0_IN      1
#define FOBJ_IN    2
#define GOBJ_IN    3
#define NCON_IN    4
#define VUB_IN     5
#define VLB_IN     6
#define ITMAX_IN   7
#define PARAM_IN   8
#define LAST_PARAM 8

#define X_OPT_OUT  Rhs+1
#define F_OPT_OUT  Rhs+2
#define DF_OPT_OUT Rhs+3
#define G_OPT_OUT  Rhs+4
#define DG_OPT_OUT Rhs+5
#define IC_OPT_OUT Rhs+6

extern void C2F(get_ct)(double *);
extern void C2F(get_obj)(double *);
extern void C2F(get_igoto)(int *);
extern void C2F(get_nac)(int *);
extern void C2F(set_delfun)(double *);
extern void C2F(set_dabfun)(double *);
extern void C2F(set_fdch)(double *);
extern void C2F(set_fdchm)(double *);
extern void C2F(set_ct)(double *);
extern void C2F(set_ctmin)(double *);
extern void C2F(set_ctl)(double *);
extern void C2F(set_ctlmin)(double *);
extern void C2F(set_alphax)(double *);
extern void C2F(set_abobj1)(double *);
extern void C2F(set_theta)(double *);
extern void C2F(set_obj)(double *);
extern void C2F(set_ndv)(int *);
extern void C2F(set_ncon)(int *);
extern void C2F(set_nside)(int *);
extern void C2F(set_iprint)(int *);
extern void C2F(set_nfdg)(int *);
extern void C2F(set_nscal)(int *);
extern void C2F(set_linobj)(int *);
extern void C2F(set_itmax)(int *);
extern void C2F(set_itrm)(int *);
extern void C2F(set_icndir)(int *);
extern void C2F(set_igoto)(int *);
extern void C2F(set_nac)(int *);
extern void C2F(set_info)(int *);
extern void C2F(set_infog)(int *);
extern void C2F(set_iter)(int *);

typedef struct {
  double delfun;
  double dabfun;
  double fdch;
  double fdchm;
  double ct;
  double ctmin;
  double ctl;
  double ctlmin;
  double alphax;
  double abobj1;
  double theta;
  double obj;
  int    ndv;
  int    ncon;
  int    nside;
  int    iprint;
  int    nfdg;
  int    nscal;
  int    linobj;
  int    itmax;
  int    itrm;
  int    icndir;
  int    igoto;
  int    nac;
  int    info;
  int    infog;
  int    iter;
} MY_CNMN1;

my_export MY_CNMN1 C2F(cnmn1);

typedef struct {
  double dm1;
  double dm2;
  double dm3;
  double dm4;
  double dm5;
  double dm6;
  double dm7;
  double dm8;
  double dm9;
  double dm10;
  double dm11;
  double dm12;
  double dct;
  double dctl;
  double phi;
  double abobj;
  double cta;
  double ctam;
  double ctbm;
  double obj1;
  double slope;
  double dx;
  double dx1;
  double fi;
  double xi;
  double dftdf1;
  double alp;
  double fff;
  double a1;
  double a2;
  double a3;
  double a4;
  double f1;
  double f2;
  double f3;
  double f4;
  double cv1;
  double cv2;
  double cv3;
  double cv4;
  double app;
  double alpca;
  double alpfes;
  double alpln;
  double alpmin;
  double alpnc;
  double alpsav;
  double alpsid;
  double alptot;
  double rspace;
  int idm1;
  int idm2;
  int idm3;
  int jdir;
  int iobj;
  int kobj;
  int kcount;
  int ncal[2];
  int nfeas;
  int mscal;
  int ncobj;
  int nvc;
  int kount;
  int icount;
  int igood1;
  int igood2;
  int igood3;
  int igood4;
  int ibest;
  int iii;
  int nlnc;
  int jgoto;
  int ispace[2];
} MY_CONSAV;

my_export MY_CONSAV C2F(consav);

//                       X         VLB       VUB       G         SCAL      DF        A         S         G1        G2        B          C
extern  void C2F(conmin)(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *,
                         int *, int *, int *, int *, int *, int *, int *, int *);
//                       ISC    IC     MS1    N1     N2     N3     N4     N5

// Example of settings from: example1.f       example2.f       example3.f
// INFOG  = 0                                 Idem             Idem
// INFO   = 0                                 Idem             Idem
// NFDG   = 0                                 NFDG   = 1       NFDG  = 2
// IPRINT = 2                                 IPRINT = 1       Idem
// NDV    = 4                                 Idem             NDV   = 2
// ITMAX  = 40                                Idem             Idem
// NCON   = 3                                 Idem             NCON  = 6
// NSIDE  = 0                                 Idem             NSIDE = 1
// ICNDIR = 0                                 Idem             Idem
// NSCAL  = 0                                 Idem             Idem
// FDCH   = 0.0                               Idem             Idem
// FDCHM  = 0.0                               Idem             Idem
// CT     = 0.0                               Idem             Idem
// CTMIN  = 0.0                               Idem             Idem
// CTL    = 0.0                               Idem             Idem
// CTLMIN = 0.0                               Idem             Idem
// THETA  = 0.0                               Idem             Idem
// DELFUN = 0.0                               Idem             Idem
// DABFUN = 0.0                               Idem             Idem
// LINOBJ = 0.0                               Idem             LINOBJ = 1
// ITRM   = 0                                 Idem             Idem
// N1     = 6                                 Idem             N1     = 4
// N2     = 11                                Idem             N2     = 10
// N3     = 11                                Idem             N3     = 10
// N4     = 11                                Idem             N4     = 10
// N5     = 22                                Idem             N5     = 20
// ALPHAX = 0.0                               Idem             Idem
// ABOBJ1 = 0.0                               Idem             Idem
// CTL    = 0.0                               Idem             Idem
// DO 5 I=1,NDV
//   X(I)   = 1.0                             Idem             Idem
//   VLB(I) = -99999.                         Idem             VLB(I) = 0.001
//   VUB(I) =  99999.                         Idem             VUB(I) = 1.0E+10
// 5 CONTINUE
//
// DO 6 J=1,NCON
//   ISC(J) = 0                               Idem             Idem
// 6 CONTINUE

// Prototypes for objective function and constraints:
// Objective function: function [f,df] = fobj(x)
// x  = [nx][1]
// f  = [1][1]
// df = [nx][1]
// Constraints:        function [g,dg,ic] = constr(x,ct)
// x  = [nx][1]
// ct = [1][1]
// g  = [nac][1]
// dg = [nx][nac] 
// ic = [nac][1]  - contains the index of active constraints
// The constraints tested must be
// g_i(x)<=ct

int sci_conmin(char * fname)
{
  int n_x,      m_x,      * pi_x     = NULL; double * x     = NULL; // Input - X + ndv
  int fobj_lhs, fobj_rhs, l_fobj; // Input - f_obj
  int gobj_lhs, gobj_rhs, l_gobj;  // Input - g_obj + ncon
  int n_vub,    m_vub,    * pi_vub   = NULL; double * vub   = NULL; // Input - VUB + NSIDE = 1 if vub and / or vlb are present
  int n_vlb,    m_vlb,    * pi_vlb   = NULL; double * vlb   = NULL; // Input - VLB
  int n_itmax,  m_itmax,  * pi_itmax = NULL; double * itmax = NULL; // Input - ITMAX
  int n_ncon,   m_ncon,   * pi_ncon  = NULL; double * ncon  = NULL; // Input - NCON

  // Temporary variables for the call to objective function and constraints
  int n_x1_tmp    = 0, m_x1_tmp    = 0; double * x1_tmp = NULL;
  int n_x2_tmp    = 0, m_x2_tmp    = 0; double * x2_tmp = NULL;
  int n_x3_tmp    = 0, m_x3_tmp    = 0; double * x3_tmp = NULL;
  int n_fobj_tmp  = 0, m_fobj_tmp  = 0, * pi_fobj_tmp  = NULL; double * fobj_tmp = NULL;
  int n_dfobj_tmp = 0, m_dfobj_tmp = 0, * pi_dfobj_tmp = NULL; double * dfobj_tmp = NULL;
  int n_gobj_tmp  = 0, m_gobj_tmp  = 0, * pi_gobj_tmp  = NULL; double * gobj_tmp = NULL;
  int n_dgobj_tmp = 0, m_dgobj_tmp = 0, * pi_dgobj_tmp = NULL; double * dgobj_tmp = NULL;
  int n_ic        = 0, m_ic        = 0, * pi_ic        = NULL; double * ic = NULL; int * ic_int = NULL;
  // Output variables
  int n_x_opt  = 0, m_x_opt  = 0; double * x_opt  = NULL; // Output
  int n_f_opt  = 0, m_f_opt  = 0; double * f_opt  = NULL; // Output - optional
  int n_df_opt = 0, m_df_opt = 0; double * df_opt = NULL; // Output - optional
  int n_g_opt  = 0, m_g_opt  = 0; double * g_opt  = NULL; // Output - optional
  int n_dg_opt = 0, m_dg_opt = 0; double * dg_opt = NULL; // Output - optional
  int n_ic_res = 0, m_ic_res = 0; double * ic_opt = NULL; // Output - optional

  int n1,n2,n3,n4,n5;
  int ibegin; // To be used in SciFunction
  int i,j;
  int tmp_int, tmp_int_2, tmp_res;
  double tmp_dbl;

  int    * pi_param = NULL;
  int    * pIsc  = NULL;
  double * pScal = NULL;
  int szIsc = 0, szScal = 0;
  int Isc_mallocated = 0, Scal_mallocated = 0;
  int nbvars_old = Nbvars;

  double * G   = NULL;
  double * A   = NULL;
  double * S   = NULL;
  double * G1  = NULL;
  double * G2  = NULL;
  double * B   = NULL;
  double * C   = NULL;
  double * DF  = NULL;
  int    * MS1 = NULL;

  SciErr _SciErr;

  // Order of the parameters:
  // x0
  // fobj
  // gobj
  // ncon
  // nacmx
  // vub
  // vlb
  // itmax
  // isc
  // scal
  // nscal
  // nfdg
  // icndir
  // fdch
  // fdchm
  // ct
  // ctmin
  // ctl
  // ctlmin
  // theta
  // delfun
  // dabfun
  // linobj
  // itrm
  // alphax
  // abobj1
  // infog
  // info
  // iprint

#ifdef DEBUG
  DBGPRINTF("conmin_optim: Lhs = %d Rhs = %d\n",Lhs,Rhs);
#endif

  if (Rhs<LAST_PARAM)
    {
      Scierror(999,"%s: %d parameters are required\n", fname, LAST_PARAM);
      return 0;
    } /* End If */

  ////////////////////////
  // Get the parameters //
  ////////////////////////

  _SciErr = getVarAddressFromPosition(pvApiCtx, X0_IN, &pi_x); CONMIN_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, pi_x, &n_x, &m_x, &x); CONMIN_ERROR;

  tmp_int = n_x*m_x;
  C2F(set_ndv)(&tmp_int);

  GetRhsVar(FOBJ_IN, EXTERNAL_DATATYPE, &fobj_lhs, &fobj_rhs, &l_fobj);
  GetRhsVar(GOBJ_IN, EXTERNAL_DATATYPE, &gobj_lhs, &gobj_rhs, &l_gobj);

  /////////////////////////////////////////////
  // Initialization of some "size" variables //
  /////////////////////////////////////////////

  _SciErr = getVarAddressFromPosition(pvApiCtx, NCON_IN, &pi_ncon); CONMIN_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, pi_ncon, &n_ncon, &m_ncon, &ncon); CONMIN_ERROR;
  tmp_int = (int)ncon[0];
  C2F(set_ncon)(&tmp_int);

  _SciErr = getVarAddressFromPosition(pvApiCtx, VUB_IN, &pi_vub); CONMIN_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, pi_vub, &n_vub, &m_vub, &vub); CONMIN_ERROR;
  if (n_vub*m_vub!=n_x*m_x)
    {
      Scierror(999,"%s: upper must have the same length as x\n", fname);
      return 0;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, VLB_IN, &pi_vlb); CONMIN_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, pi_vlb, &n_vlb, &m_vlb, &vlb); CONMIN_ERROR;
  if (n_vlb*m_vlb!=n_x*m_x)
    {
      Scierror(999,"%s: upper must have the same length as x\n", fname);
      return 0;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, ITMAX_IN, &pi_itmax); CONMIN_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, pi_itmax, &n_itmax, &m_itmax, &itmax); CONMIN_ERROR;
  tmp_int = (int)itmax[0];
  C2F(set_itmax)(&tmp_int);

  // YC: si nvub != nvlb mettre le plus petit à +inf ou -inf
  // YC: si nvub != nvlb != ndv mettre tout à +inf et -inf et taille de ndv
  // C2F(cnmn1).nside = n_vub * m_vub;
  tmp_int = 2 * n_x * m_x;
  C2F(set_nside)(&tmp_int);
  //if (C2F(cnmn1).nside < (n_vlb * m_vlb)) C2F(cnmn1).nside = n_vlb * m_vlb;

  ///////////////////////////////////////
  // Get the parameters from the plist //
  ///////////////////////////////////////

  _SciErr = initPList(pvApiCtx, PARAM_IN, &pi_param); CONMIN_ERROR;
  if (!checkPList(pvApiCtx, pi_param))
    {
      Scierror(999, "%s: argument n° %d is not a plist\n", fname, PARAM_IN);
      
      return 0;
    }

  tmp_int_2 = ncon[0] + 2 * n_x * m_x + 1;

  _SciErr = getIntInPList(pvApiCtx, pi_param, "nacmx", &tmp_int, &tmp_res, 0, 0, CHECK_NONE); CONMIN_ERROR;
  n3 = tmp_int;

  pIsc = (int *)MALLOC(10*sizeof(int)); // Bug in scilab: default value vector not mallocated -> workaround with a memory leak ...
  _SciErr = getColVectorOfIntInPList(pvApiCtx, pi_param, "isc", pIsc, &tmp_res, 0, 1, &szIsc, 0, CHECK_NONE); CONMIN_ERROR;

  if (!((tmp_res!=-1) && (szIsc==n_x*m_x)))
    {
      if (pIsc) FREE(pIsc);
      // We must allocate our own pIsc
      pIsc = (int *)MALLOC(n_x*m_x*sizeof(int));
      for(i=0;i<n_x*m_x;i++) pIsc[i] = 0;
      Isc_mallocated = 1;
    }

  pScal = (double *)MALLOC(10*sizeof(double)); // Bug in scilab: default value vector not mallocated -> workaround with a memory leak ...
  _SciErr = getColVectorOfDoubleInPList(pvApiCtx, pi_param, "scal", pScal, &tmp_res, 0, 1, &szScal, 0, CHECK_NONE); CONMIN_ERROR;

  if (!((tmp_res!=-1) && (szScal==n_x*m_x)))
    {
      if (pScal) FREE(pScal);
      // We must allocate our own pScal
      pScal = (double *)MALLOC(n_x*m_x*sizeof(double));
      for(i=0;i<n_x*m_x;i++) pScal[i] = 1.0;
      Scal_mallocated = 1;
    }

  _SciErr = getIntInPList(pvApiCtx, pi_param, "nscal", &tmp_int, &tmp_res, 0, 0, CHECK_NONE); CONMIN_ERROR;
  if (tmp_res!=-1) C2F(set_nscal)(&tmp_int);

  _SciErr = getIntInPList(pvApiCtx, pi_param, "nfdg", &tmp_int, &tmp_res, 0, 0, CHECK_NONE); CONMIN_ERROR;
  if (tmp_res!=-1) C2F(set_nfdg)(&tmp_int);

  _SciErr = getIntInPList(pvApiCtx, pi_param, "icndir", &tmp_int, &tmp_res, n_x*m_x+1, 0, CHECK_NONE); CONMIN_ERROR;
  if (tmp_res!=-1) C2F(set_icndir)(&tmp_int);

  _SciErr = getDoubleInPList(pvApiCtx, pi_param, "fdch", &tmp_dbl, &tmp_res, 0.01, 0, CHECK_NONE); CONMIN_ERROR;
  if (tmp_res!=-1) C2F(set_fdch)(&tmp_dbl);

  _SciErr = getDoubleInPList(pvApiCtx, pi_param, "fdchm", &tmp_dbl, &tmp_res, 0.01, 0, CHECK_NONE); CONMIN_ERROR;
  if (tmp_res!=-1) C2F(set_fdchm)(&tmp_dbl);

  _SciErr = getDoubleInPList(pvApiCtx, pi_param, "ct", &tmp_dbl, &tmp_res, -0.1, 0, CHECK_NONE); CONMIN_ERROR;
  if (tmp_res!=-1) C2F(set_ct)(&tmp_dbl);

  _SciErr = getDoubleInPList(pvApiCtx, pi_param, "ctmin", &tmp_dbl, &tmp_res, 0.004, 0, CHECK_NONE); CONMIN_ERROR;
  if (tmp_res!=-1) C2F(set_ctmin)(&tmp_dbl);

  _SciErr = getDoubleInPList(pvApiCtx, pi_param, "ctl", &tmp_dbl, &tmp_res, -0.01, 0, CHECK_NONE); CONMIN_ERROR;
  if (tmp_res!=-1) C2F(set_ctl)(&tmp_dbl);

  _SciErr = getDoubleInPList(pvApiCtx, pi_param, "ctlmin", &tmp_dbl, &tmp_res, 0.001, 0, CHECK_NONE); CONMIN_ERROR;
  if (tmp_res!=-1) C2F(set_ctlmin)(&tmp_dbl);

  _SciErr = getDoubleInPList(pvApiCtx, pi_param, "theta", &tmp_dbl, &tmp_res, 1.0, 0, CHECK_NONE); CONMIN_ERROR;
  if (tmp_res!=-1) C2F(set_theta)(&tmp_dbl);

  _SciErr = getDoubleInPList(pvApiCtx, pi_param, "delfun", &tmp_dbl, &tmp_res, 0.001, 0, CHECK_NONE); CONMIN_ERROR;
  if (tmp_res!=-1) C2F(set_delfun)(&tmp_dbl);

  _SciErr = getDoubleInPList(pvApiCtx, pi_param, "dabfun", &tmp_dbl, &tmp_res, 0.001, 0, CHECK_NONE); CONMIN_ERROR;
  if (tmp_res!=-1) C2F(set_dabfun)(&tmp_dbl);

  _SciErr = getIntInPList(pvApiCtx, pi_param, "linobj", &tmp_int, &tmp_res, 0, 0, CHECK_NONE); CONMIN_ERROR;
  if (tmp_res!=-1) C2F(set_linobj)(&tmp_int);

  _SciErr = getIntInPList(pvApiCtx, pi_param, "itrm", &tmp_int, &tmp_res, 3, 0, CHECK_NONE); CONMIN_ERROR;
  if (tmp_res!=-1) C2F(set_itrm)(&tmp_int);

  _SciErr = getDoubleInPList(pvApiCtx, pi_param, "alphax", &tmp_dbl, &tmp_res, 0.3, 0, CHECK_NONE); CONMIN_ERROR;
  if (tmp_res!=-1) C2F(set_alphax)(&tmp_dbl);

  _SciErr = getDoubleInPList(pvApiCtx, pi_param, "abobj1", &tmp_dbl, &tmp_res, 0.2, 0, CHECK_NONE); CONMIN_ERROR;
  if (tmp_res!=-1) C2F(set_abobj1)(&tmp_dbl);

  _SciErr = getIntInPList(pvApiCtx, pi_param, "infog", &tmp_int, &tmp_res, 0, 0, CHECK_NONE); CONMIN_ERROR;
  if (tmp_res!=-1) C2F(set_infog)(&tmp_int);

  _SciErr = getIntInPList(pvApiCtx, pi_param, "info", &tmp_int, &tmp_res, 1, 0, CHECK_NONE); CONMIN_ERROR;
  if (tmp_res!=-1) C2F(set_info)(&tmp_int);

  _SciErr = getIntInPList(pvApiCtx, pi_param, "iprint", &tmp_int, &tmp_res, 0, 0, CHECK_NONE); CONMIN_ERROR;
  if (tmp_res!=-1) C2F(set_iprint)(&tmp_int);

  // Verification: 
  // size(x) = size(vlb) = size(vub) 
  // size(x) = size(scal) ou (size(scal)==0 et nscal = 0)
  // size(isc) = ncon ou 0

  n1 = n_x*m_x + 2;
  n2 = ncon[0] + 2*n_x*m_x;
  n4 = (n3>n_x*m_x)?n3:n_x*m_x;
  n5 = 2*n4;

  /////////////////////////////////////////////////////
  // First call to objective function and constraint //
  /////////////////////////////////////////////////////

  ////////////////////////////////
  // Call to objective function //
  ////////////////////////////////

  ibegin = Rhs + 1;
  Nbvars = ibegin + MAX(fobj_rhs,fobj_lhs);

  // Create the first variable: x
  n_x1_tmp = n_x; m_x1_tmp = m_x;
  _SciErr = allocMatrixOfDouble(pvApiCtx, ibegin, n_x1_tmp, m_x1_tmp, &x1_tmp); CONMIN_ERROR;
  memcpy(x1_tmp,x,n_x*m_x*sizeof(double));

  // Create a fake 2nd variable on the stack. This will allow to call GetRhsVar to retrieve the 2nd output argument
  n_x2_tmp = 1; m_x2_tmp = 1;
  _SciErr = allocMatrixOfDouble(pvApiCtx, ibegin + 1, n_x2_tmp, m_x2_tmp, &x2_tmp); CONMIN_ERROR;

  // Call to the scilab function
  SciFunction(&ibegin,&l_fobj,&fobj_lhs,&fobj_rhs);
  if (Err>0) 
    {
      sciprint("conmin: error when calling objective function\n");
      return 0;
    } /* End If */
  
  // YC: we should add a check on fobj_lhs and fobj_rhs to check the scilab prototype of fobj

  // Get fobj
  _SciErr = getVarAddressFromPosition(pvApiCtx, ibegin, &pi_fobj_tmp); CONMIN_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, pi_fobj_tmp, &n_fobj_tmp, &m_fobj_tmp, &fobj_tmp); CONMIN_ERROR;
  // Get gradient of fobj
  _SciErr = getVarAddressFromPosition(pvApiCtx, ibegin + 1, &pi_dfobj_tmp); CONMIN_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, pi_dfobj_tmp, &n_dfobj_tmp, &m_dfobj_tmp, &dfobj_tmp); CONMIN_ERROR;

  Nbvars = nbvars_old;

  // Fill C2F(cnmn1).obj: objective function value
  C2F(set_obj)(fobj_tmp);
  
  // DF(N1) Analytic gradient of the objective function for the current decision variables, X(I).  DF(I) contains the partial derivative
  //        of OBJ with respect to X(I).  Calculate DF(I), I = 1, NDV if INFO = 3 or INFO = 4 and if NFDG = 0 or NFDG = 2.
  
  DF = (double *)MALLOC(n1*sizeof(double));
  
  // Fill DF(N1): gradient of objective function value
  memcpy(DF,dfobj_tmp,n_x*m_x*sizeof(double));

  ////////////////////////
  // Call to constraint //
  ////////////////////////

  ibegin = Rhs + 1; // YC: for SciFunction
  Nbvars = ibegin + MAX(gobj_rhs,gobj_lhs);

  // YC: test nombre de param en entree (2) / sortie (3)

  // First parameter: the vector x
  n_x1_tmp = n_x; m_x1_tmp = m_x;
  _SciErr = allocMatrixOfDouble(pvApiCtx, ibegin, n_x1_tmp, m_x1_tmp, &x1_tmp); CONMIN_ERROR;
  memcpy(x1_tmp,x,n_x*m_x*sizeof(double));

  // Second parameter: the scalar CT (Constraint Threshold)
  n_x2_tmp = 1; m_x2_tmp = 1;
  _SciErr = allocMatrixOfDouble(pvApiCtx, ibegin + 1, n_x2_tmp, m_x2_tmp, &x2_tmp); CONMIN_ERROR;
  C2F(get_ct)(x2_tmp);
  
  // Third parameter: a fake parameter which will receive IC (index of active constraints)
  n_x3_tmp = 1; m_x3_tmp = 1;
  _SciErr = allocMatrixOfDouble(pvApiCtx, ibegin + 2, n_x3_tmp, m_x3_tmp, &x3_tmp); CONMIN_ERROR;

  SciFunction(&ibegin,&l_gobj,&gobj_lhs,&gobj_rhs);
  if (Err>0) 
    {
      sciprint("conmin: error when calling constraints function\n");
      return 0;
    } /* End If */

  // YC: we should add a check on gobj_lhs and gobj_rhs to check the scilab prototype of gobj

  // After the first call to the constraints, we get the number of constraints via the first variable returned by constraint function
  // Get gobj
  _SciErr = getVarAddressFromPosition(pvApiCtx, ibegin, &pi_gobj_tmp); CONMIN_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, pi_gobj_tmp, &n_gobj_tmp, &m_gobj_tmp, &gobj_tmp); CONMIN_ERROR;
  // Get gradient of gobj
  _SciErr = getVarAddressFromPosition(pvApiCtx, ibegin + 1, &pi_dgobj_tmp); CONMIN_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, pi_dgobj_tmp, &n_dgobj_tmp, &m_dgobj_tmp, &dgobj_tmp); CONMIN_ERROR;
  // Get indexes of active constraints
  _SciErr = getVarAddressFromPosition(pvApiCtx, ibegin + 2, &pi_ic); CONMIN_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, pi_ic, &n_ic, &m_ic, &ic); CONMIN_ERROR;
  if (ic_int) FREE(ic_int);
  ic_int = (int *)MALLOC(n_ic*m_ic*sizeof(int));
  for(i=0; i<m_ic*m_ic; i++) ic_int[i] = (int)ic[i];

  Nbvars = nbvars_old;

  // dgobj_tmp is of size [ndv][nac]
  tmp_int = n_ic*m_ic;
  C2F(set_nac)(&tmp_int);
  
  if (tmp_int>0)
    {
      // In scilab, the index start from 1. We need an index which start from 0
      for(j=0;j<tmp_int;j++)
	{
	  ic[j] = ic[j] - 1;
	} /* End For */
    } /* End If */

  ///////////////////////
  // Memory allocation //
  ///////////////////////

  // S(N1)    Move direction in the NDV-dimensional optimization space. S(I) gives the rate at which variable X(I) changes with respect toALPHA.  
  S  = (double *)MALLOC(n1*sizeof(double));
  
  // G1(N2)   Not used if NCON = NSIDE = NSCAL = 0.  Used for temporary storage of constraint values G(J), J = 1, NCON and decision
  //          variables X(I), I = 1, NDV.
  G1 = (double *)MALLOC(n2*sizeof(double));
  
  // G2(N2)   Not used if NCON = NSIDE = 0.  Used for temporary storage of constraint values G(J), J = 1, NCON.
  G2 = (double *)MALLOC(n2*sizeof(double));
  
  // B(N3,N3) Not used if NCON = NSIDE = 0.  Used in determining direction vector S for constrained minimization problems.  Array B may
  //          be used for temporary storage in external routine SUB1.
  B = (double *)MALLOC(n3*n3*sizeof(double));

  // C(N4)    Not used in NCON = NSIDE = 0.  Used with array B in determining direction vector S for constrained minimization problems.  Used
  //          for temporary storage of vector X if NSCAL.NE.0. routine SUB1.
  C   = (double *)MALLOC(n4*sizeof(double));

  // MS1(N5)  Not used if NCON = NSIDE = 0.  Used with array B in determining direction vector S for constrained minimization problems.  Array
  //          MS1 may be used for temporary storage in external routine SUB1.
  MS1 = (int *)MALLOC(n5*sizeof(int));

  // G(N2)  Not used if NCON = NSIDE = 0.  Vector containing all constraint functions, G(J), J = 1, NCON for current decision variables, X.
  //        Calculate G(J), J = 1, NCON if INFO = 2.
  G  = (double *)MALLOC(n2*sizeof(double));
  // Fill G(N2): value of constraint
  for(i=0;i<ncon[0];i++) G[i] = gobj_tmp[i];

  // A(N1,N3)  Not used if NCON = NSIDE = 0.  Gradients of active or violated constraints, for current decision variables, X(I).
  //           A(J,I) contains the gradient of the Jth active or violated constraint, G(J), with respect to the Ith decision variable,
  //           X(I) for J = 1, NAC and I = 1, NDV.  Calculate if INFO = 4 and NFDG = 0.
  // YC: adress if 0 bytes after a block of size 48 alloc'd
  A = (double *)MALLOC(n1*n3*sizeof(double));

  // Fill A(N1,N3): matrix of gradient of active constraint
  if (tmp_int>0) // tmp_int is still equal to C2F(cnmn1).nac
    {
      memcpy(A,dgobj_tmp,n_dgobj_tmp*m_dgobj_tmp*sizeof(double));
    } /* End If */

  tmp_int = 0;
  C2F(set_obj)(fobj_tmp);
  C2F(set_igoto)(&tmp_int);
  C2F(set_iter)(&tmp_int);

  // Start the optimization

  // Optimization loop
  for(i=0;i<itmax[0];i++)
    {
#ifdef DEBUG
      DBGPRINTF("Iteration %d / %d\n",i,itmax[0]);
#endif
      // Call to objective function and constraint function. 
      // Create a new variable: x1_tmp on top of the stack - outside of the loop
      // Call SciFunction on fobj with ibegin=29
      // Get fobj and dfobj on stack at position 30 and 31
      // Call SciFunction on gobj with ibegin=29
      // Get gobj and dgobj on stack at position 30 and 31 

      if (i>0) 
	{
	  // We skip the first iteration because we already have made a call to objective function and constraints

	  ////////////////////////////////
	  // Call to objective function //
	  ////////////////////////////////

	  ibegin = Rhs + 1; // YC: for SciFunction
	  Nbvars = ibegin + MAX(fobj_rhs,fobj_lhs);

	  // First parameter: the vector x
	  n_x1_tmp = n_x; m_x1_tmp = m_x;
	  _SciErr = allocMatrixOfDouble(pvApiCtx, ibegin, n_x1_tmp, m_x1_tmp, &x1_tmp); CONMIN_ERROR;
	  memcpy(x1_tmp,x,n_x*m_x*sizeof(double));

	  // Create a fake 2nd variable on the stack. This will allow to call GetRhsVar to retrieve the 2nd output argument
	  n_x2_tmp = 1; m_x2_tmp = 1;
	  _SciErr = allocMatrixOfDouble(pvApiCtx, ibegin + 1, n_x2_tmp, m_x2_tmp, &x2_tmp); CONMIN_ERROR;

	  SciFunction(&ibegin,&l_fobj,&fobj_lhs,&fobj_rhs);
	  if (Err>0) 
	    {
	      sciprint("conmin: error when calling objective function at iteration %d / %d\n",i,itmax[0]);
	      return 0;
	    } /* End If */

	  // Get fobj
	  _SciErr = getVarAddressFromPosition(pvApiCtx, ibegin, &pi_fobj_tmp); CONMIN_ERROR;
	  _SciErr = getMatrixOfDouble(pvApiCtx, pi_fobj_tmp, &n_fobj_tmp, &m_fobj_tmp, &fobj_tmp); CONMIN_ERROR;
	  // Get gradient of fobj
	  _SciErr = getVarAddressFromPosition(pvApiCtx, ibegin + 1, &pi_dfobj_tmp); CONMIN_ERROR;
	  _SciErr = getMatrixOfDouble(pvApiCtx, pi_dfobj_tmp, &n_dfobj_tmp, &m_dfobj_tmp, &dfobj_tmp); CONMIN_ERROR;

	  Nbvars = nbvars_old;

	  // Fill C2F(cnmn1).obj: objective function value
	  C2F(set_obj)(fobj_tmp);
	  // Fill DF(N1): gradient of objective function value
	  memcpy(DF,dfobj_tmp,n_x*m_x*sizeof(double));

	  //////////////////////////////////
	  // Call to constraints function //
	  //////////////////////////////////

	  ibegin = Rhs + 1; // YC: for SciFunction
	  Nbvars = ibegin + MAX(gobj_rhs,gobj_lhs);

	  // First parameter: the vector x
	  n_x1_tmp = n_x; m_x1_tmp = m_x;
	  _SciErr = allocMatrixOfDouble(pvApiCtx, ibegin, n_x1_tmp, m_x1_tmp, &x1_tmp); CONMIN_ERROR;
	  memcpy(x1_tmp,x,n_x*m_x*sizeof(double));
	
	  // Second parameter: the scalar CT (Constraint Threshold)
	  n_x2_tmp = 1; m_x2_tmp = 1;
	  _SciErr = allocMatrixOfDouble(pvApiCtx, ibegin + 1, n_x2_tmp, m_x2_tmp, &x2_tmp); CONMIN_ERROR;
	  C2F(get_ct)(x2_tmp);
	  
	  // Third parameter: a fake parameter which will receive the IC (index of active constraint) list
	  n_x3_tmp = 1; m_x3_tmp = 1;
	  _SciErr = allocMatrixOfDouble(pvApiCtx, ibegin + 2, n_x3_tmp, m_x3_tmp, &x3_tmp); CONMIN_ERROR;

	  SciFunction(&ibegin,&l_gobj,&gobj_lhs,&gobj_rhs);
	  if (Err>0) 
	    {
	      sciprint("conmin: error when calling constraints function at iteration %d / %d\n",i,itmax[0]);
	      return 0;
	    } /* End If */

	  // Now, we get the constraint values and the derivate of the constraints
	  // Get gobj
	  _SciErr = getVarAddressFromPosition(pvApiCtx, ibegin, &pi_gobj_tmp); CONMIN_ERROR;
	  _SciErr = getMatrixOfDouble(pvApiCtx, pi_gobj_tmp, &n_gobj_tmp, &m_gobj_tmp, &gobj_tmp); CONMIN_ERROR;
	  // Get gradient of gobj
	  _SciErr = getVarAddressFromPosition(pvApiCtx, ibegin + 1, &pi_dgobj_tmp); CONMIN_ERROR;
	  _SciErr = getMatrixOfDouble(pvApiCtx, pi_dgobj_tmp, &n_dgobj_tmp, &m_dgobj_tmp, &dgobj_tmp); CONMIN_ERROR;
	  // Get indexes of active constraints
	  _SciErr = getVarAddressFromPosition(pvApiCtx, ibegin + 2, &pi_ic); CONMIN_ERROR;
	  _SciErr = getMatrixOfDouble(pvApiCtx, pi_ic, &n_ic, &m_ic, &ic); CONMIN_ERROR;
	  if (ic_int) FREE(ic_int);
	  ic_int = (int *)MALLOC(n_ic*m_ic*sizeof(int));
	  for(i=0; i<m_ic*m_ic; i++) ic_int[i] = (int)ic[i];

	  Nbvars = nbvars_old;

	  // dgobj_tmp is of size [ndv][nac]
	  tmp_int = n_ic*m_ic;
	  C2F(set_nac)(&tmp_int);

	  // Fill G(N2): value of constraints
	  memcpy(G,gobj_tmp,ncon[0]*sizeof(double));

	  if (tmp_int>0) // tmp_int is still equal to nac
	    {
	      // In scilab, the index start from 1. We need an index which start from 0
	      for(j=0;j<tmp_int;j++)
		{
		  ic[j] = ic[j] - 1;
		} /* End For */

	      // Fill A(N1,N3): matrix of gradient of active constraint
	      memcpy(A,dgobj_tmp,n_dgobj_tmp*m_dgobj_tmp*sizeof(double));
	    } /* End If */
	} /* End If */
      
      //                        X         VLB       VUB       G         SCAL      DF        A         S         G1        G2        B         C
      //extern void C2F(conmin)(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *,
      //                        int *, int *, int *, int *, int *, int *, int *, int *);
      //                        ISC    IC     MS1    N1     N2     N3     N4     N5

      // YC: revoir ic: mauvais type
      C2F(conmin)(x,vlb,vub,G,pScal,DF,A,S,G1,G2,B,C,pIsc,ic_int,MS1,&n1,&n2,&n3,&n4,&n5);

      C2F(get_igoto)(&tmp_int);
      if (tmp_int == 0) break; // Optimization has converged
    } /* End For */

  //////////////////////////
  // return the variables //
  //////////////////////////

  // x_opt
  n_x_opt = n_x; m_x_opt = m_x;
  _SciErr = allocMatrixOfDouble(pvApiCtx, X_OPT_OUT, n_x_opt, m_x_opt, &x_opt); CONMIN_ERROR;
  memcpy(x_opt,x,n_x*m_x*sizeof(double));

  // f_opt
  n_f_opt = 1; m_f_opt = 1;
  _SciErr = allocMatrixOfDouble(pvApiCtx, F_OPT_OUT, n_f_opt, m_f_opt, &f_opt); CONMIN_ERROR;
  C2F(get_obj)(f_opt);

  // df_opt
  n_df_opt = n_x; m_df_opt = m_x;
  _SciErr = allocMatrixOfDouble(pvApiCtx, DF_OPT_OUT, n_df_opt, m_df_opt, &df_opt); CONMIN_ERROR;
  memcpy(df_opt,DF,n_x*m_x*sizeof(double));
  
  // g_opt
  n_g_opt = ncon[0]; m_g_opt = 1;
  _SciErr = allocMatrixOfDouble(pvApiCtx, G_OPT_OUT, n_g_opt, m_g_opt, &g_opt); CONMIN_ERROR;
  memcpy(g_opt,G,ncon[0]*sizeof(double));
  
  // dg_opt
  n_dg_opt = n_dgobj_tmp; m_dg_opt = m_dgobj_tmp;
  _SciErr = allocMatrixOfDouble(pvApiCtx, DG_OPT_OUT, n_dg_opt, m_dg_opt, &dg_opt); CONMIN_ERROR;
  // Now, we stored the gradient of active constraint in the matrix
  memcpy(dg_opt,A,n_dgobj_tmp*m_dgobj_tmp*sizeof(double));

  // IC
  C2F(get_nac)(&tmp_int);
  if (tmp_int!=0)
    {
      n_ic_res = tmp_int; m_ic_res = 1;
      _SciErr = allocMatrixOfDouble(pvApiCtx, IC_OPT_OUT, n_ic_res, m_ic_res, &ic_opt); CONMIN_ERROR;
      // Now, we stored the gradient of active constraint in the matrix
      memcpy(ic_opt,ic_int,tmp_int*sizeof(int));
    }
  else
    {
      n_ic_res = 1; m_ic_res = 1;
      _SciErr = allocMatrixOfDouble(pvApiCtx, IC_OPT_OUT, n_ic_res, m_ic_res, &ic_opt); CONMIN_ERROR;
      ic_opt[0] = -1;
    }

  LhsVar(1) = X_OPT_OUT;
  LhsVar(2) = F_OPT_OUT;
  LhsVar(3) = DF_OPT_OUT;
  LhsVar(4) = G_OPT_OUT;
  LhsVar(5) = DG_OPT_OUT;
  LhsVar(6) = IC_OPT_OUT;

  ///////////////////////////
  // Free allocated memory //
  ///////////////////////////

  if (G)   FREE(G);
  if (A)   FREE(A);
  if (DF)  FREE(DF);
  if (S)   FREE(S);
  if (G1)  FREE(G1);
  if (G2)  FREE(G2);
  if (B)   FREE(B);
  if (C)   FREE(C);
  if (MS1) FREE(MS1);

  if (Isc_mallocated)  FREE(pIsc);
  if (Scal_mallocated) FREE(pScal);

  return 0;
}
