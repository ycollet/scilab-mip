/* Sample scilab interface functio for getting functions, gradients, and dense Hessians from an AMPL .nl file.
   Start with:

	[asl,x,bl,bu,v,cl,cu] = ampl_init(<filename>)

   to read in a problem (discarding the previous problem, if any).
   The return values are:

	x = primal initial guess
	v = dual initial guess
	bl, bu = lower and upper bounds on x
	cl, cu = lower and upper bounds on c (the constraint bodies).

   Then

	[f,c] = ampl_evalf(asl,x)

   gives the function (f) and constraint bodies (c) at x;

	[g,Jac] = ampl_evalg(asl,x)

   gives the gradient g of f, the Jacobian matrix J of c, and

	W = ampl_evalw(asl,v)

   gives the Hessian W of the Lagrangian function L = f + v*c
   (at the last x at which amplfunc(x,1) was called).

   After finding optimal values for x and v,

	ampl_write_sol(asl,'solution message',x,v)

   to write a stub.sol file.
*/

#ifdef WIN32
#include <windows.h>
#endif

#include <string.h>

#include <stack-c.h>
#include <sciprint.h>
#include <Scierror.h>
#include <freeArrayOfString.h>
#include <MALLOC.h>

#include <api_common.h>
#include <api_pointer.h>
#include <api_string.h>
#include <api_double.h>
#include <api_list.h>
#include <api_sparse.h>

//#include <asl.h>
#include <asl_pfgh.h>
#include <ampl_ops.h>

#define DEBUG 1
#define DBGPRINTF sciprint
//#define DBGPRINTF printf

#ifdef _MSC_VER
#define strdup _strdup
/* Omit sw "signal" catching and x86 precision adjustment. */
#define ASL_NO_FP_INIT
#include <fpinit.c>
#endif /* _MSC_VER */

#ifdef DEBUG
#define CHECK_ERROR(MSG)				\
  DBGPRINTF("%s\n",MSG);				\
  if(_SciErr.iErr)					\
    {							\
      printError(&_SciErr, 0);				\
      return 0;						\
    }
#else
#define CHECK_ERROR(MSG)			\
  if(_SciErr.iErr)				\
    {						\
      printError(&_SciErr, 0);			\
      return 0;					\
    }
#endif

///////////////////////////////////////////////////////////////
// Load and get the initial values of the given AMPL problem //
///////////////////////////////////////////////////////////////

struct my_asl
{
  ASL * asl_pointer;
  real * J;
};

struct my_asl_fgh
{
  ASL_fgh * asl_pointer;
  real * J;
};

int sci_ampl_init(char * fname)
{
  int m_string, n_string, * pstring = NULL, * string_len = NULL;
  int m_X0,     n_X0;
  int m_LUv,    n_LUv;
  int m_Uvx,    n_Uvx;
  int m_pi0,    n_pi0;
  int m_LUrhs,  n_LUrhs;
  int m_Urhsx,  n_Urhsx;
  int i, n, nc, nz, Type;
  char ** String = NULL;

  struct my_asl * my_asl_pointer = NULL;
  FILE * nl;

  SciErr _SciErr;
  _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &pstring);
  if(_SciErr.iErr)
    {
      printError(&_SciErr, 0);
      return 0;
    }

  _SciErr = getVarType(pvApiCtx, pstring, &Type); CHECK_ERROR("ampl_init - 1");

  // The first parameter is a string, then, we open the given file
  if (((Rhs!=1)&&(Type!=sci_strings))||(Lhs!=7))
    {
      sciprint("%s: [asl,x,bl,bu,v,cl,cu] = ampl_init(filename)\n", fname);
      return 0;
    }

  _SciErr = getMatrixOfString(pvApiCtx, pstring, &m_string, &n_string, NULL, NULL);  CHECK_ERROR("ampl_init - 2");

  string_len = (int *)MALLOC(m_string*n_string*sizeof(int));
  _SciErr = getMatrixOfString(pvApiCtx, pstring, &m_string, &n_string, string_len, NULL); CHECK_ERROR("ampl_init - 3");

  String = (char **)MALLOC(m_string*n_string*sizeof(char *));
  for(i=0; i<m_string*n_string; i++)
    {
      String[i] = (char *)MALLOC((string_len[i]+1)*sizeof(char));
    }
  _SciErr = getMatrixOfString(pvApiCtx, pstring, &m_string, &n_string, string_len, String); CHECK_ERROR("ampl_init - 4");

#ifdef DEBUG
  sciprint("DEBUG: |%s| - %d - %d\n", String[0], m_string, n_string);
#endif

  my_asl_pointer = (struct my_asl *)MALLOC(1*sizeof(struct my_asl));

  my_asl_pointer->asl_pointer = ASL_alloc(ASL_read_pfgh);

  /* 0 ==> jacdim0 should exit if filename.nl does not exist; 
     1 ==> return 0 */
  my_asl_pointer->asl_pointer->i.return_nofile_ = 1;
  
  if (!(nl = jac0dim_ASL(my_asl_pointer->asl_pointer, String[0],string_len[0]))) 
    {
      Scierror(999,"%s: can't open %s\n", fname, String[0]);
      ASL_free(&(my_asl_pointer->asl_pointer));
      if (string_len) FREE(string_len);
      for(i=0;i<m_string*n_string; i++)
	{
	  if (String[i]) FREE(String[i]);
	}
      if (String) FREE(String);

      return 0;
    }

  if (my_asl_pointer->asl_pointer->i.n_obj_ <= 0) sciprint("%s: warning: objective == 0\n", fname);
  
  n  = my_asl_pointer->asl_pointer->i.n_var_;
  nc = my_asl_pointer->asl_pointer->i.n_con_;
  nz = my_asl_pointer->asl_pointer->i.nzc_;

  my_asl_pointer->J = (real *)M1alloc_ASL(&(my_asl_pointer->asl_pointer->i),nz*sizeof(real));

  //  nhnz = sphsetup(0, 0, nc > 0, 0); // For the sparse things.

  m_X0    = n;  n_X0    = 1;
  m_LUv   = n;  n_LUv   = 1;
  m_Uvx   = n;  n_Uvx   = 1;
  m_pi0   = nc; n_pi0   = 1;
  m_LUrhs = nc; n_LUrhs = 1;
  m_Urhsx = nc; n_Urhsx = 1;
  
  my_asl_pointer->asl_pointer->i.X0_    = (double *)MALLOC(m_X0*n_X0*sizeof(double));
  my_asl_pointer->asl_pointer->i.LUv_   = (double *)MALLOC(m_LUv*n_LUv*sizeof(double));
  my_asl_pointer->asl_pointer->i.Uvx_   = (double *)MALLOC(m_Uvx*n_Uvx*sizeof(double));
  my_asl_pointer->asl_pointer->i.pi0_   = (double *)MALLOC(m_pi0*n_pi0*sizeof(double));
  my_asl_pointer->asl_pointer->i.LUrhs_ = (double *)MALLOC(m_LUrhs*n_LUrhs*sizeof(double));
  my_asl_pointer->asl_pointer->i.Urhsx_ = (double *)MALLOC(m_Urhsx*n_Urhsx*sizeof(double));

  pfgh_read_ASL(my_asl_pointer->asl_pointer, nl, ASL_findgroups | ASL_find_co_class | ASL_return_read_err);

  // ASL_readerr_none    = 0, all went well
  // ASL_readerr_nofile  = 1, cannot open .nl file
  // ASL_readerr_nonlin  = 2, model involves nonlinearities (ed0read)
  // ASL_readerr_argerr  = 3, user-defined function with bad args
  // ASL_readerr_unavail = 4, user-defined function not available
  // ASL_readerr_corrupt = 5, corrupt .nl file
  // ASL_readerr_bug	 = 6, bug in .nl reader
  // ASL_readerr_CLP     = 7  solver cannot handle CLP extensions

  _SciErr = createPointer(pvApiCtx, Rhs+1, (void*)my_asl_pointer); CHECK_ERROR("ampl_init - 5");
  _SciErr = createMatrixOfDouble(pvApiCtx, Rhs+2, m_X0, n_X0, my_asl_pointer->asl_pointer->i.X0_); CHECK_ERROR("ampl_init - 6");
  _SciErr = createMatrixOfDouble(pvApiCtx, Rhs+3, m_LUv, n_LUv, my_asl_pointer->asl_pointer->i.LUv_); CHECK_ERROR("ampl_init - 7");
  _SciErr = createMatrixOfDouble(pvApiCtx, Rhs+4, m_Uvx, n_Uvx, my_asl_pointer->asl_pointer->i.Uvx_); CHECK_ERROR("ampl_init - 8");
  _SciErr = createMatrixOfDouble(pvApiCtx, Rhs+5, m_pi0, n_pi0, my_asl_pointer->asl_pointer->i.pi0_); CHECK_ERROR("ampl_init - 9");
  _SciErr = createMatrixOfDouble(pvApiCtx, Rhs+6, m_LUrhs, n_LUrhs, my_asl_pointer->asl_pointer->i.LUrhs_); CHECK_ERROR("ampl_init - 10");
  _SciErr = createMatrixOfDouble(pvApiCtx, Rhs+7, m_Urhsx, n_Urhsx, my_asl_pointer->asl_pointer->i.Urhsx_); CHECK_ERROR("ampl_init - 11");

  LhsVar(1) = Rhs + 1;
  LhsVar(2) = Rhs + 2;
  LhsVar(3) = Rhs + 3;
  LhsVar(4) = Rhs + 4;
  LhsVar(5) = Rhs + 5;
  LhsVar(6) = Rhs + 6;
  LhsVar(7) = Rhs + 7;

#ifdef DEBUG
  sciprint("DEBUG: |%s| - %d - %d\n", String[0], m_string, n_string);
#endif

  if (string_len) FREE(string_len);
  for(i=0;i<m_string*n_string; i++)
    {
      if (String[i]) FREE(String[i]);
    }
  if (String) FREE(String);
  
  return 0;
}

//////////////////////////////////////////////////////
// Release the memory allocated to the AMPL problem //
//////////////////////////////////////////////////////

int sci_ampl_free(char * fname)
{
  int  * p_pointer   = NULL;
  void * asl_pointer = NULL;
  SciErr _SciErr;
  struct my_asl * my_asl_pointer = NULL;

  if (Rhs!=1)
    {
      sciprint("%s usage: %s(asl)\n", fname, fname);
      return 0;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_pointer); CHECK_ERROR("ampl_free - 1");
  _SciErr = getPointer(pvApiCtx, p_pointer, &asl_pointer); CHECK_ERROR("ampl_free - 2");

  my_asl_pointer = (struct my_asl *)asl_pointer;

  ASL_free(&(my_asl_pointer->asl_pointer));

  if (my_asl_pointer) FREE(my_asl_pointer);

  return 0;
}

///////////////////////////////////////////
// Computation of the objective function //
///////////////////////////////////////////

int sci_ampl_evalf(char * fname)
{
  int  * p_pointer   = NULL;
  void * asl_pointer = NULL;
  int m_x, n_x, * p_x = NULL;
  int m_f, n_f;
  int m_c, n_c;
  int n, nc;
  fint nerror;
  double * x = NULL, * f = NULL, * c = NULL;
  char * what;
  SciErr _SciErr;
  struct my_asl * my_asl_pointer = NULL;
  Jmp_buf err_jmp0;
  
  if (Rhs!=2) 
    {
      sciprint("%s usage: [f,c] = %s(asl, x)", fname, fname);
      return 0;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_pointer); CHECK_ERROR("ampl_evalf - 1");
  _SciErr = getPointer(pvApiCtx, p_pointer, &asl_pointer); CHECK_ERROR("ampl_evalf - 2");

  my_asl_pointer = (struct my_asl *)asl_pointer;

  nerror = -1;
  what = "i";
  my_asl_pointer->asl_pointer->i.err_jmp1_ = &err_jmp0;

  n  = my_asl_pointer->asl_pointer->i.n_var_;
  nc = my_asl_pointer->asl_pointer->i.n_con_;

  if (!my_asl_pointer->asl_pointer->i.filename_)
    {
      sciprint("%s: ampl_init(filename) has not been called\n", fname);
      return 0;
    }
      
  CheckVector(2,n,1);
  
  _SciErr = getVarAddressFromPosition(pvApiCtx, 2, &p_x); CHECK_ERROR("ampl_evalf - 3");
  _SciErr = getMatrixOfDouble(pvApiCtx, p_x, &m_x, &n_x, &x); CHECK_ERROR("ampl_evalf - 4");

  m_f = 1;  n_f = 1;
  m_c = nc; n_c = 1;
  
  f = (double *)MALLOC(m_f*n_f*sizeof(double));
  c = (double *)MALLOC(m_c*n_c*sizeof(double));

  if (setjmp(err_jmp0.jb))
    {
      sciprint("%s: trouble evaluating %s\n", fname, what);
      return 0;
    }

  what = "f";
  *f = (*((ASL*)(my_asl_pointer->asl_pointer))->p.Objval)((ASL*)my_asl_pointer->asl_pointer,0,x,&nerror);
  if (nerror!=-1)
    {
      sciprint("%s: objective function value - error %d\n", fname, nerror);
    }

  what = "c";
  (*((ASL*)(my_asl_pointer->asl_pointer))->p.Conval)((ASL*)(my_asl_pointer->asl_pointer),x,c,&nerror);
  if (nerror!=-1)
    {
      sciprint("%s: constraint function value - error %d\n", fname, nerror);
    }
  
  _SciErr = createMatrixOfDouble(pvApiCtx, Rhs+1, m_f, n_f, f); CHECK_ERROR("ampl_evalf - 5");
  _SciErr = createMatrixOfDouble(pvApiCtx, Rhs+2, m_c, n_c, c); CHECK_ERROR("ampl_evalf - 6");

  LhsVar(1) = Rhs + 1;
  LhsVar(2) = Rhs + 2;

  if (f) FREE(f);
  if (c) FREE(c);

  return 0;
}

///////////////////////////////////
// Dense version of the Jacobian //
///////////////////////////////////

int sci_ampl_evalg(char * fname)
{
  int  * p_pointer   = NULL;
  void * asl_pointer = NULL;
  int m_x,  n_x, * p_x  = NULL;
  int m_g,  n_g;
  int m_J1, n_J1;
  double * x = NULL, * g = NULL, * J1 = NULL;
  int n, nc, nz;
  fint nerror;
  char   * what = NULL;
  double * tmp_dbl = NULL;
  real   * J = NULL;
  SciErr _SciErr;
  struct my_asl * my_asl_pointer = NULL;
  cgrad  * cg = NULL, ** cgp = NULL, ** cgpe = NULL;
  Jmp_buf err_jmp0;

  if (Rhs!=2) 
    {
      sciprint("%s usage: [g,Jac] = %s(asl,x)\n", fname, fname);
      return 0;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_pointer); CHECK_ERROR("ampl_evalg - 1");
  _SciErr = getPointer(pvApiCtx, p_pointer, &asl_pointer); CHECK_ERROR("ampl_evalg - 2");

  my_asl_pointer = (struct my_asl *)asl_pointer;

  if (!my_asl_pointer->asl_pointer->i.filename_)
    {
      sciprint("%s: ampl_init(filename) has not been called\n", fname);
      return 0;
    }
      
  n  = my_asl_pointer->asl_pointer->i.n_var_;
  nc = my_asl_pointer->asl_pointer->i.n_con_;
  nz = my_asl_pointer->asl_pointer->i.nzc_;

  CheckVector(2,n,1);

  _SciErr = getVarAddressFromPosition(pvApiCtx, 2, &p_x); CHECK_ERROR("ampl_evalg - 3");
  _SciErr = getMatrixOfDouble(pvApiCtx, p_x, &m_x, &n_x, &x); CHECK_ERROR("ampl_evalg - 4");

  J = my_asl_pointer->J;

  nerror = -1;
  my_asl_pointer->asl_pointer->i.err_jmp1_ = &err_jmp0;

  m_g  = n;  n_g  = 1;
  m_J1 = nc; n_J1 = n;

  g  = (double *)MALLOC(m_g*n_g*sizeof(double));
  J1 = (double *)MALLOC(m_J1*n_J1*sizeof(double));

  if (setjmp(err_jmp0.jb))
    {
      sciprint("%s: trouble evaluating %s\n", fname, what);
      return 0;
    }

  what = "g";
  (*((ASL*)(my_asl_pointer->asl_pointer))->p.Objgrd)((ASL*)(my_asl_pointer->asl_pointer),0,x,g,&nerror);
  if (nerror!=-1)
    {
      sciprint("%s: gradient value - error %d\n", fname, nerror);
    }
  
  if (nc) 
    {
      what = "J";
      (*((ASL*)(my_asl_pointer->asl_pointer))->p.Jacval)((ASL*)(my_asl_pointer->asl_pointer),x,J,&nerror);
      if (nerror!=-1)
	{
	  sciprint("%s: Jacobian value - error %d\n", fname, nerror);
	}
      cgp = my_asl_pointer->asl_pointer->i.Cgrad_;

      memset(J1,0,n*nc*sizeof(double));
      tmp_dbl = J1;
      for(cgpe=cgp+nc; cgp<cgpe; tmp_dbl++)
	for(cg=*cgp++; cg; cg=cg->next)
	  tmp_dbl[nc*cg->varno] = J[cg->goff];
    }
  
  _SciErr = createMatrixOfDouble(pvApiCtx, Rhs+1, m_g,  n_g,  g); CHECK_ERROR("ampl_evalg - 5");
  _SciErr = createMatrixOfDouble(pvApiCtx, Rhs+2, m_J1, n_J1, J1); CHECK_ERROR("ampl_evalg - 6");

  LhsVar(1) = Rhs + 1;
  LhsVar(2) = Rhs + 2;

  if (g)  FREE(g);
  if (J1) FREE(J1);

  return 0;
}

//////////////////////////////////
// Dense version of the Hessian //
//////////////////////////////////

int sci_ampl_evalw(char * fname)
{
  int  * p_pointer   = NULL;
  void * asl_pointer = NULL;
  int m_v, n_v, * p_v = NULL;
  int m_W, n_W;
  int n, nc, nz;
  fint nerror;
  char * what;
  double * v = NULL, * W = NULL;
  SciErr _SciErr;
  struct my_asl * my_asl_pointer = NULL;
  Jmp_buf err_jmp0;
  
  if ((Lhs!=1)||(Rhs!=2))
    {
      sciprint("%s usage: W = %s(asl, v)\n", fname, fname);
      return 0;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_pointer); CHECK_ERROR("ampl_evalw - 1");
  _SciErr = getPointer(pvApiCtx, p_pointer, &asl_pointer); CHECK_ERROR("ampl_evalw - 2");

  my_asl_pointer = (struct my_asl *)asl_pointer;

  n  = my_asl_pointer->asl_pointer->i.n_var_;
  nc = my_asl_pointer->asl_pointer->i.n_con_;
  nz = my_asl_pointer->asl_pointer->i.nzc_;

  CheckVector(2,nc,1);
  
  nerror = -1;
  what = "i";
  my_asl_pointer->asl_pointer->i.err_jmp1_ = &err_jmp0;

  if (!my_asl_pointer->asl_pointer->i.filename_)
    {
      sciprint("%s: ampl_init(filename) has not been called\n", fname);
      return 0;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, 2, &p_v); CHECK_ERROR("ampl_evalw - 3");
  _SciErr = getMatrixOfDouble(pvApiCtx, p_v, &m_v, &n_v, &v); CHECK_ERROR("ampl_evalw - 4");

  m_W = n; n_W = n;

  W = (double *)MALLOC(n_W*m_W*sizeof(double));

  if (setjmp(err_jmp0.jb))
    {
      sciprint("%s: trouble evaluating %s\n", fname, what);
      return 0;
    }

  what = "W";
  (*((ASL*)(my_asl_pointer->asl_pointer))->p.Fulhes)((ASL*)(my_asl_pointer->asl_pointer),W,n,0,0,v);

  _SciErr = createMatrixOfDouble(pvApiCtx, Rhs+1, m_W,  n_W,  W); CHECK_ERROR("ampl_evalw - 5");

  LhsVar(1) = Rhs + 1;

  if (W) FREE(W);

  return 0;
}

////////////////////////////////////
// Sparse version of the Jacobian //
////////////////////////////////////

int sci_ampl_eval_sp_g(char * fname)
{
  int  * p_pointer   = NULL;
  void * asl_pointer = NULL;
  int m_x, n_x, * p_x = NULL;
  int m_g, n_g;
  int n, nc, nz, i;
  fint nerror;
  double * x = NULL, * g = NULL;
  char   * what        = NULL;
  int * mnel = NULL, * icol = NULL;
  double * R = NULL;
  SciErr _SciErr;
  struct my_asl * my_asl_pointer = NULL;
  cgrad  * cg = NULL, ** cgp = NULL;
  Jmp_buf err_jmp0;

  if (Rhs!=2) 
    {
      sciprint("%s usage: [g,spJac] = %s(asl,x)\n", fname, fname);
      return 0;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_pointer); CHECK_ERROR("ampl_eval_sp_g - 1");
  _SciErr = getPointer(pvApiCtx, p_pointer, &asl_pointer); CHECK_ERROR("ampl_eval_sp_g - 2");

  my_asl_pointer = (struct my_asl *)asl_pointer;

  if (!my_asl_pointer->asl_pointer->i.filename_)
    {
      sciprint("%s: ampl_init(filename) has not been called\n", fname);
      return 0;
    }
      
  n  = my_asl_pointer->asl_pointer->i.n_var_;
  nc = my_asl_pointer->asl_pointer->i.n_con_;
  nz = my_asl_pointer->asl_pointer->i.nzc_;

  CheckVector(2,n,1);

  _SciErr = getVarAddressFromPosition(pvApiCtx, 2, &p_x); CHECK_ERROR("ampl_eval_sp_g - 3");
  _SciErr = getMatrixOfDouble(pvApiCtx, p_x, &m_x, &n_x, &x); CHECK_ERROR("ampl_eval_sp_g - 4");

  nerror = -1;
  my_asl_pointer->asl_pointer->i.err_jmp1_ = &err_jmp0;

  m_g  = n;  n_g  = 1;
  g = (double *)MALLOC(m_g*n_g*sizeof(double));

  if (setjmp(err_jmp0.jb))
    {
      sciprint("%s: trouble evaluating %s\n", fname, what);
      return 0;
    }

  // Gradient of the objective function (dense)
  what = "g";
  (*((ASL*)(my_asl_pointer->asl_pointer))->p.Objgrd)((ASL*)(my_asl_pointer->asl_pointer),0,x,g,&nerror);
  if (nerror!=-1)
    {
      sciprint("%s: gradient value - error %d\n", fname, nerror);
    }

  // Jacobian of the constraints (sparse)

  mnel = (int *)MALLOC(n*sizeof(int));
  icol = (int *)MALLOC(nz*sizeof(int));
  R    = (double *)MALLOC(nz*sizeof(double));

  if (nc) 
    {
      what = "J";

      memset(R,0,nz*sizeof(double));

      (*((ASL*)(my_asl_pointer->asl_pointer))->p.Jacval)((ASL*)(my_asl_pointer->asl_pointer),x, R, &nerror); // R est de taille nz
      if (nerror!=-1)
	{
	  sciprint("%s: Jacobian value - error %d\n", fname, nerror);
	}

      for(i=0;i<n;i++)
	{
	  mnel[i] = my_asl_pointer->asl_pointer->i.A_colstarts_[i+1] - my_asl_pointer->asl_pointer->i.A_colstarts_[i];
	}

      cgp = my_asl_pointer->asl_pointer->i.Cgrad_;
      for(i=0;i<nc;i++)
	for(cg=*cgp++;cg;cg=cg->next)
	  {
	    icol[cg->goff] = i + 1;
	  }
    }
  
#ifdef DEBUG
  DBGPRINTF("DEBUG: n = %d nc = %d nz = %d\n", n, nc, nz);
#endif

  _SciErr = createMatrixOfDouble(pvApiCtx, Rhs+1, m_g,  n_g,  g); CHECK_ERROR("ampl_eval_sp_g - 5");
  _SciErr = createSparseMatrix(pvApiCtx, Rhs+2, n, nc, nz, mnel, icol, R); CHECK_ERROR("ampl_eval_sp_g - 6");

  LhsVar(1) = Rhs + 1;
  LhsVar(2) = Rhs + 2;
  
  if (mnel) FREE(mnel);
  if (icol) FREE(icol);
  if (R)    FREE(R);
  if (g)    FREE(g);

  return 0;
}

///////////////////////////////////
// Sparse version of the Hessian //
///////////////////////////////////

int sci_ampl_eval_sp_w(char * fname)
{
  int  * p_pointer   = NULL;
  void * asl_pointer = NULL;
  int m_v, n_v, * p_v = NULL;
  int n, nc, nz, nhnz, i;
  fint nerror;
  char * what;
  double * v = NULL;
  SciErr _SciErr;
  struct my_asl * my_asl_pointer = NULL;
  Jmp_buf err_jmp0;
  int * mnel = NULL, * icol = NULL;
  double * R = NULL;

  if ((Lhs!=1)||(Rhs!=2))
    {
      sciprint("%s usage: spW = %s(asl, v)\n", fname, fname);
      return 0;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_pointer); CHECK_ERROR("ampl_eval_sp_w - 1");
  _SciErr = getPointer(pvApiCtx, p_pointer, &asl_pointer); CHECK_ERROR("ampl_eval_sp_w - 2");

  my_asl_pointer = (struct my_asl *)asl_pointer;

  n  = my_asl_pointer->asl_pointer->i.n_var_;
  nc = my_asl_pointer->asl_pointer->i.n_con_;
  nz = my_asl_pointer->asl_pointer->i.nzc_;

  CheckVector(2,nc,1);
  
  nerror = -1;
  what = "i";
  my_asl_pointer->asl_pointer->i.err_jmp1_ = &err_jmp0;

  if (!my_asl_pointer->asl_pointer->i.filename_)
    {
      sciprint("%s: ampl_init(filename) has not been called\n", fname);
      return 0;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, 2, &p_v); CHECK_ERROR("ampl_eval_sp_w - 3");
  _SciErr = getMatrixOfDouble(pvApiCtx, p_v, &m_v, &n_v, &v); CHECK_ERROR("ampl_eval_sp_w - 4");

  if (setjmp(err_jmp0.jb))
    {
      sciprint("%s: trouble evaluating %s\n", fname, what);
      return 0;
    }

  nhnz = (*((ASL*)(my_asl_pointer->asl_pointer))->p.Sphset)((ASL*)(my_asl_pointer->asl_pointer), 0, 0, 0, nc > 0, 0);

  mnel = (int *)MALLOC(n*sizeof(int));
  icol = (int *)MALLOC(nhnz*sizeof(int));
  R    = (double *)MALLOC(nhnz*sizeof(double));

  my_asl_pointer->asl_pointer->i.sputinfo_->uptri = 0; // We want all the Hessian (upper and lower triangles)

  what = "W";
  if (nc > 0) (*((ASL*)(my_asl_pointer->asl_pointer))->p.Sphes)((ASL*)(my_asl_pointer->asl_pointer),0,R,0,0,v);
  else        (*((ASL*)(my_asl_pointer->asl_pointer))->p.Sphes)((ASL*)(my_asl_pointer->asl_pointer),0,R,0,0,0);

  for(i=0;i<n;i++) 
    mnel[i] = my_asl_pointer->asl_pointer->i.sputinfo_->hcolstarts[i+1] - my_asl_pointer->asl_pointer->i.sputinfo_->hcolstarts[i];

  // YC: tester sputinfo == NULL et hrownos == NULL
  for(i=0;i<nhnz;i++) 
    icol[i] = my_asl_pointer->asl_pointer->i.sputinfo_->hrownos[i]+1;

  _SciErr = createSparseMatrix(pvApiCtx, Rhs+1, n, n, nhnz, mnel, icol, R); CHECK_ERROR("ampl_eval_sp_w - 5");

  LhsVar(1) = Rhs + 1;

  if (mnel) FREE(mnel);
  if (icol) FREE(icol);
  if (R)    FREE(R);

  return 0;
}

////////////////////////////////////////////////////
// Returns the sparsity structure of the Jacobian //
////////////////////////////////////////////////////

int sci_ampl_eval_spst_g_rc(char * fname)
{
  int  * p_pointer   = NULL;
  void * asl_pointer = NULL;
  int m_x,    n_x, * p_x = NULL;
  int m_irow, n_irow;
  int m_jcol, n_jcol;
  int n, nc, nz, i, j, Index;
  fint nerror;
  char   * what        = NULL;
  double * tmp_dbl     = NULL;
  double * irow = NULL, * jcol = NULL, * x = NULL;
  struct my_asl * my_asl_pointer = NULL;
  cgrad  * cg = NULL, ** cgp = NULL;
  Jmp_buf err_jmp0;
  SciErr _SciErr;

  if (Rhs!=2) 
    {
      sciprint("%s usage: [irow,jcol] = %s(asl,x)\n", fname, fname);
      return 0;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_pointer); CHECK_ERROR("ampl_eval_spst_rc - 1");
  _SciErr = getPointer(pvApiCtx, p_pointer, &asl_pointer); CHECK_ERROR("ampl_eval_spst_rc - 2");

  my_asl_pointer = (struct my_asl *)asl_pointer;

  if (!my_asl_pointer->asl_pointer->i.filename_)
    {
      sciprint("%s: ampl_init(filename) has not been called\n", fname);
      return 0;
    }
      
  n  = my_asl_pointer->asl_pointer->i.n_var_;
  nc = my_asl_pointer->asl_pointer->i.n_con_;
  nz = my_asl_pointer->asl_pointer->i.nzc_;

  CheckVector(2,n,1);

  _SciErr = getVarAddressFromPosition(pvApiCtx, 2, &p_x); CHECK_ERROR("ampl_eval_spst_rc - 3");
  _SciErr = getMatrixOfDouble(pvApiCtx, p_x, &m_x, &n_x, &x); CHECK_ERROR("ampl_eval_spst_rc - 4");

  nerror = -1;
  my_asl_pointer->asl_pointer->i.err_jmp1_ = &err_jmp0;

  if (setjmp(err_jmp0.jb))
    {
      sciprint("%s: trouble evaluating %s\n", fname, what);
      return 0;
    }

  // Sparsity structure of the Jacobian

  m_irow = nz; n_irow = 1;
  m_jcol = nz; n_jcol = 1;
  
  irow = (double *)MALLOC(m_irow*n_irow*sizeof(double));
  jcol = (double *)MALLOC(m_jcol*n_jcol*sizeof(double));
  tmp_dbl = (double *)MALLOC(nz*sizeof(double));

  if (nc) 
    {
      what = "J";

      (*((ASL*)(my_asl_pointer->asl_pointer))->p.Jacval)((ASL*)(my_asl_pointer->asl_pointer), x, tmp_dbl, &nerror);
      if (nerror!=-1)
	{
	  sciprint("%s: Jacobian value - error %d\n", fname, nerror);
	}

      Index = 0;
      for(i=0;i<n;i++)
	{
	  for(j=0;j<my_asl_pointer->asl_pointer->i.A_colstarts_[i+1] - my_asl_pointer->asl_pointer->i.A_colstarts_[i];j++)
	    {
	      jcol[Index] = (double)(i + 1);
	      Index++;
	    }
	}

      cgp = my_asl_pointer->asl_pointer->i.Cgrad_;
      for(i=0;i<nc;i++)
	for(cg=*cgp++;cg;cg=cg->next)
	  {
	    irow[cg->goff] = (double)(i + 1);
	  }
    }
  
  _SciErr = createMatrixOfDouble(pvApiCtx, Rhs+1, m_irow,  n_irow,  irow); CHECK_ERROR("ampl_eval_spst_rc - 5");
  _SciErr = createMatrixOfDouble(pvApiCtx, Rhs+2, m_jcol,  n_jcol,  jcol); CHECK_ERROR("ampl_eval_spst_rc - 6");

  LhsVar(1) = Rhs + 1;
  LhsVar(2) = Rhs + 2;
  
  if (tmp_dbl) FREE(tmp_dbl);
  if (irow)    FREE(irow);
  if (jcol)    FREE(jcol);

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// Returns the values corresponding to the sparsity structure of the Jacobian //
////////////////////////////////////////////////////////////////////////////////

int sci_ampl_eval_spst_g_val(char * fname)
{
  int  * p_pointer   = NULL;
  void * asl_pointer = NULL;
  int m_x,   n_x, * p_x   = NULL;
  int m_val, n_val;
  int n, nc, nz;
  fint nerror;
  char   * what        = NULL;
  double * x = NULL, * val = NULL;
  struct my_asl * my_asl_pointer = NULL;
  Jmp_buf err_jmp0;
  SciErr _SciErr;

  if (Rhs!=2) 
    {
      sciprint("%s usage: val = (asl,x)\n", fname, fname);
      return 0;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_pointer); CHECK_ERROR("ampl_eval_spst_val - 1");
  _SciErr = getPointer(pvApiCtx, p_pointer, &asl_pointer); CHECK_ERROR("ampl_eval_spst_val - 2");

  my_asl_pointer = (struct my_asl *)asl_pointer;

  if (!my_asl_pointer->asl_pointer->i.filename_)
    {
      sciprint("%s: ampl_init(filename) has not been called\n", fname);
      return 0;
    }
      
  n  = my_asl_pointer->asl_pointer->i.n_var_;
  nc = my_asl_pointer->asl_pointer->i.n_con_;
  nz = my_asl_pointer->asl_pointer->i.nzc_;

  CheckVector(2,n,1);

  _SciErr = getVarAddressFromPosition(pvApiCtx, 2, &p_x); CHECK_ERROR("ampl_eval_spst_val - 3");
  _SciErr = getMatrixOfDouble(pvApiCtx, p_x, &m_x, &n_x, &x); CHECK_ERROR("ampl_eval_spst_val - 4");

  nerror = -1;
  my_asl_pointer->asl_pointer->i.err_jmp1_ = &err_jmp0;

  if (setjmp(err_jmp0.jb))
    {
      sciprint("%s: trouble evaluating %s\n", fname, what);
      return 0;
    }

  // Sparsity structure of the Jacobian

  m_val = nz; n_val = 1;
  val   = (double *)MALLOC(m_val*n_val*sizeof(double));
  memset(val,0,m_val*n_val*sizeof(double));

  if (nc) 
    {
      what = "J";

      (*((ASL*)(my_asl_pointer->asl_pointer))->p.Jacval)((ASL*)(my_asl_pointer->asl_pointer),x, val, &nerror);
      if (nerror!=-1)
	{
	  sciprint("%s: Jacobian value - error %d\n", fname, nerror);
	}
    }

  _SciErr = createMatrixOfDouble(pvApiCtx, Rhs+1, m_val, n_val, val); CHECK_ERROR("ampl_eval_spst_val - 5");

  LhsVar(1) = Rhs + 1;
  
  if (val) FREE(val);

  return 0;
}

///////////////////////////////////////////////////
// Returns the sparsity structure of the Hessian //
///////////////////////////////////////////////////

int sci_ampl_eval_spst_w_rc(char * fname)
{
  int * p_pointer = NULL;
  void * asl_pointer = NULL;
  int m_v,    n_v, * p_v = NULL;
  int m_irow, n_irow;
  int m_jcol, n_jcol;
  int n, nc, nz, nhnz = 0;
  int i, j, Index;
  fint nerror;
  char   * what        = NULL;
  double * tmp_dbl     = NULL;
  double * irow = NULL, * jcol = NULL, * v = NULL;
  struct my_asl * my_asl_pointer = NULL;
  Jmp_buf err_jmp0;
  SciErr _SciErr;

  if ((Lhs!=1)||(Rhs!=2))
    {
      sciprint("%s usage: [irow,jcol] = %s(asl, v)\n", fname, fname);
      return 0;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_pointer); CHECK_ERROR("ampl_eval_spst_w_rc - 1");
  _SciErr = getPointer(pvApiCtx, p_pointer, &asl_pointer); CHECK_ERROR("ampl_eval_spst_w_rc - 2");

  my_asl_pointer = (struct my_asl *)asl_pointer;

  n  = my_asl_pointer->asl_pointer->i.n_var_;
  nc = my_asl_pointer->asl_pointer->i.n_con_;
  nz = my_asl_pointer->asl_pointer->i.nzc_;

  CheckVector(2,nc,1);
  
  nerror = -1;
  my_asl_pointer->asl_pointer->i.err_jmp1_ = &err_jmp0;

  if (!my_asl_pointer->asl_pointer->i.filename_)
    {
      sciprint("%s: ampl_init(filename) has not been called\n", fname);
      return 0;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, 2, &p_v); CHECK_ERROR("ampl_eval_spst_w_rc - 3");
  _SciErr = getMatrixOfDouble(pvApiCtx, p_v, &m_v, &n_v, &v); CHECK_ERROR("ampl_eval_spst_w_rc - 4");

  if (setjmp(err_jmp0.jb))
    {
      sciprint("%s: trouble evaluating %s\n", fname, what);
      return 0;
    }

  // Sparsity structure of the Jacobian

  m_irow = nz; n_irow = 1;
  m_jcol = nz; n_jcol = 1;

  irow = (double *)MALLOC(m_irow*n_irow*sizeof(double));
  jcol = (double *)MALLOC(m_jcol*n_jcol*sizeof(double));
  tmp_dbl = (double *)MALLOC(nhnz*sizeof(double));

  nhnz = (*((ASL*)(my_asl_pointer->asl_pointer))->p.Sphset)((ASL*)(my_asl_pointer->asl_pointer), 0, 0, 0, nc > 0, 0);

  my_asl_pointer->asl_pointer->i.sputinfo_->uptri = 0; // We want all the Hessian (upper and lower triangles)

  what = "W";
  if (nc > 0) (*((ASL*)(my_asl_pointer->asl_pointer))->p.Sphes)((ASL*)(my_asl_pointer->asl_pointer),0,tmp_dbl,0,0,v);
  else        (*((ASL*)(my_asl_pointer->asl_pointer))->p.Sphes)((ASL*)(my_asl_pointer->asl_pointer),0,tmp_dbl,0,0,0);

  Index = 0;
  for(i=0;i<n;i++)
    {
      for(j=0;j<my_asl_pointer->asl_pointer->i.sputinfo_->hcolstarts[i+1] - my_asl_pointer->asl_pointer->i.sputinfo_->hcolstarts[i];j++)
	{
	  jcol[Index] = (double)(i + 1);
	  Index++;
	}
    }

  for(i=0;i<nhnz;i++) irow[i] = (double)(my_asl_pointer->asl_pointer->i.sputinfo_->hrownos[i] + 1);

  _SciErr = createMatrixOfDouble(pvApiCtx, Rhs+1, m_irow,  n_irow,  irow); CHECK_ERROR("ampl_eval_spst_w_rc - 5");
  _SciErr = createMatrixOfDouble(pvApiCtx, Rhs+2, m_jcol,  n_jcol,  jcol); CHECK_ERROR("ampl_eval_spst_w_rc - 6");

  LhsVar(1) = Rhs + 1;
  LhsVar(2) = Rhs + 2;
  
  if (tmp_dbl) FREE(tmp_dbl);
  if (irow)    FREE(irow);
  if (jcol)    FREE(jcol);

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
// Returns the values corresponding to the sparsity structure of the Hessian //
///////////////////////////////////////////////////////////////////////////////

int sci_ampl_eval_spst_w_val(char * fname)
{
  int * p_pointer = NULL;
  void * asl_pointer = NULL;
  int m_v,   n_v,   * p_v   = NULL;
  int m_val, n_val;
  int n, nc, nz, nhnz = 0;
  fint nerror;
  double * v = NULL, * val = NULL;
  char   * what        = NULL;
  struct my_asl * my_asl_pointer = NULL;
  Jmp_buf err_jmp0;
  SciErr _SciErr;

  if ((Lhs!=1)||(Rhs!=2))
    {
      sciprint("%s usage: val = %s(asl, v)\n", fname, fname);
      return 0;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_pointer); CHECK_ERROR("ampl_eval_spst_w_val - 1");
  _SciErr = getPointer(pvApiCtx, p_pointer, &asl_pointer); CHECK_ERROR("ampl_eval_spst_w_val - 2");

  my_asl_pointer = (struct my_asl *)asl_pointer;

  n  = my_asl_pointer->asl_pointer->i.n_var_;
  nc = my_asl_pointer->asl_pointer->i.n_con_;
  nz = my_asl_pointer->asl_pointer->i.nzc_;

  CheckVector(2,nc,1);
  
  nerror = -1;
  my_asl_pointer->asl_pointer->i.err_jmp1_ = &err_jmp0;

  if (!my_asl_pointer->asl_pointer->i.filename_)
    {
      sciprint("%s: ampl_init(filename) has not been called\n", fname);
      return 0;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, 2, &p_v); CHECK_ERROR("ampl_eval_spst_w_val - 3");
  _SciErr = getMatrixOfDouble(pvApiCtx, p_v, &m_v, &n_v, &v); CHECK_ERROR("ampl_eval_spst_w_val - 4");
  
  if (setjmp(err_jmp0.jb))
    {
      sciprint("%s: trouble evaluating %s\n", fname, what);
      return 0;
    }

  // Sparsity structure of the Jacobian

  m_val = nhnz; n_val = 1;

  val = (double *)MALLOC(m_val*n_val*sizeof(double));

  nhnz = (*((ASL*)(my_asl_pointer->asl_pointer))->p.Sphset)((ASL*)(my_asl_pointer->asl_pointer), 0, 0, 0, nc > 0, 0);

  my_asl_pointer->asl_pointer->i.sputinfo_->uptri = 0; // We want all the Hessian (upper and lower triangles)

  what = "W";
  if (nc > 0) (*((ASL*)(my_asl_pointer->asl_pointer))->p.Sphes)((ASL*)(my_asl_pointer->asl_pointer),0,val,0,0,v);
  else        (*((ASL*)(my_asl_pointer->asl_pointer))->p.Sphes)((ASL*)(my_asl_pointer->asl_pointer),0,val,0,0,0);

  _SciErr = createMatrixOfDouble(pvApiCtx, Rhs+1, m_val, n_val, val); CHECK_ERROR("ampl_eval_spst_w_val - 5");

  LhsVar(1) = Rhs + 1;

  if (val) FREE(val);

  return 0;
}

///////////////////////////
// Computation of hvcomp //
///////////////////////////

int sci_ampl_eval_hvcomp(char * fname)
{
  int  * p_pointer   = NULL;
  void * asl_pointer = NULL;
  int m_p,  n_p,  * p_p  = NULL;
  int m_ow, n_ow, * p_ow = NULL;
  int m_y,  n_y,  * p_y  = NULL;
  int m_hv, n_hv;
  fint nerror;
  char   * what = NULL;
  double * p = NULL, * hv = NULL, * OW = NULL, * Y = NULL;
  int    nv, nc, nf;
  SciErr _SciErr;
  struct my_asl * my_asl_pointer = NULL;
  Jmp_buf err_jmp0;
  
  if ((Rhs<2) && (Rhs>4))
    {
      sciprint("%s usage: [hv] = %s(asl,p[,OW[,Y]])", fname, fname);
      return 0;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_pointer); CHECK_ERROR("ampl_eval_hvcomp - 1");
  _SciErr = getPointer(pvApiCtx, p_pointer, &asl_pointer); CHECK_ERROR("ampl_eval_hvcomp - 2");

  my_asl_pointer = (struct my_asl *)asl_pointer;

  nerror = -1;
  what = "i";
  my_asl_pointer->asl_pointer->i.err_jmp1_ = &err_jmp0;

  nv = my_asl_pointer->asl_pointer->i.n_var_;
  nc = my_asl_pointer->asl_pointer->i.n_con_;
  nf = my_asl_pointer->asl_pointer->i.n_obj_;

  if (!my_asl_pointer->asl_pointer->i.filename_)
    {
      sciprint("%s: ampl_init(filename) has not been called\n", fname);
      return 0;
    }
      
  CheckVector(2,nv,1);
  
  _SciErr = getVarAddressFromPosition(pvApiCtx, 2, &p_p); CHECK_ERROR("ampl_eval_hvcomp - 3");
  _SciErr = getMatrixOfDouble(pvApiCtx, p_p, &m_p, &n_p, &p); CHECK_ERROR("ampl_eval_hvcomp - 4");

  if (Rhs>=3)
    {
      _SciErr = getVarAddressFromPosition(pvApiCtx, 3, &p_ow); CHECK_ERROR("ampl_eval_hvcomp - 5");
      _SciErr = getMatrixOfDouble(pvApiCtx, p_ow, &m_ow, &n_ow, &OW); CHECK_ERROR("ampl_eval_hvcomp - 6");
      if (m_ow*n_ow!=nf)
        {
          Scierror(999,"%s: the OW parameter must be of the size of the objective vector\n", fname);
          return 0;
        }
    }

  if (Rhs>=4)
    {
      _SciErr = getVarAddressFromPosition(pvApiCtx, 4, &p_y); CHECK_ERROR("ampl_eval_hvcomp - 7");
      _SciErr = getMatrixOfDouble(pvApiCtx, p_y, &m_y, &n_y, &Y); CHECK_ERROR("ampl_eval_hvcomp - 8");
      if (m_y*n_y!=nc)
        {
          Scierror(999,"%s: the Y parameter must be of the size of the constraint vector\n", fname);
          return 0;
        }
    }

  if (setjmp(err_jmp0.jb))
    {
      sciprint("%s: trouble evaluating %s\n", fname, what);
      return 0;
    }

  m_hv = nv;
  n_hv = 1;
  _SciErr = allocMatrixOfDouble(pvApiCtx, Rhs+1, m_hv, n_hv, &hv); CHECK_ERROR("ampl_eval_hvcomp - 9");

  (*((ASL*)(my_asl_pointer->asl_pointer))->p.Hvcomp)((ASL*)(my_asl_pointer->asl_pointer), hv, p, 0, OW, Y);

  LhsVar(1) = Rhs + 1;

  return 0;
}

///////////////////////
// Utility functions //
///////////////////////

int sci_ampl_write_sol(char * fname)
{
  int * p_pointer;
  void * asl_pointer = NULL;
  int m_string,  n_string,  * p_string;
  int m_x,       n_x,       * p_x;
  int m_v,       n_v,       * p_v;
  int Type;
  int * string_len = NULL;
  char ** String = NULL;
  double * x = NULL, * v = NULL;
  struct my_asl * my_asl_pointer = NULL;
  int i, n, nc;
  SciErr _SciErr;

  /* ampl_write_sol(asl, 'solution message', x, v): x = primal, v = dual */

  _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_pointer); CHECK_ERROR("ampl_write_sol - 1");
  _SciErr = getPointer(pvApiCtx, p_pointer, &asl_pointer); CHECK_ERROR("ampl_write_sol - 2");

  my_asl_pointer = (struct my_asl *)asl_pointer;

  n  = my_asl_pointer->asl_pointer->i.n_var_;
  nc = my_asl_pointer->asl_pointer->i.n_con_;

  CheckVector(3,n,1);
  CheckVector(4,nc,1);

  if (!my_asl_pointer->asl_pointer->i.filename_)
    {
      sciprint("%s: ampl_init(filename) has not been called\n", fname);
      return 0;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, 2, &p_string); CHECK_ERROR("ampl_write_sol - 3");
  _SciErr = getVarType(pvApiCtx, p_string, &Type); CHECK_ERROR("ampl_write_sol - 4");

  if ((Type!=sci_strings)||(Rhs!=4))
    {
      sciprint("%s usage: %s(asl, 'solution message', x, v)\n", fname, fname);
      return 0;
    }

  _SciErr = getMatrixOfString(pvApiCtx, p_string, &m_string, &n_string, NULL, NULL); CHECK_ERROR("ampl_write_sol - 5");

  string_len = (int *)MALLOC(m_string*n_string*sizeof(int));
  _SciErr = getMatrixOfString(pvApiCtx, p_string, &m_string, &n_string, string_len, NULL); CHECK_ERROR("ampl_write_sol - 6");
  
  String = (char **)MALLOC(m_string*n_string*sizeof(char *));
  for(i=0; i<m_string*n_string; i++)
    {
      String[i] = (char *)MALLOC((string_len[i]+1)*sizeof(char));
    }
  _SciErr = getMatrixOfString(pvApiCtx, p_string, &m_string, &n_string, string_len, String); CHECK_ERROR("ampl_write_sol - 7");

  _SciErr = getVarAddressFromPosition(pvApiCtx, 3, &p_x); CHECK_ERROR("ampl_write_sol - 8");
  _SciErr = getMatrixOfDouble(pvApiCtx, p_x, &m_x, &n_x, &x); CHECK_ERROR("ampl_write_sol - 9");
  _SciErr = getVarAddressFromPosition(pvApiCtx, 4, &p_v); CHECK_ERROR("ampl_write_sol - 10");
  _SciErr = getMatrixOfDouble(pvApiCtx, p_v, &m_v, &n_v, &v); CHECK_ERROR("ampl_write_sol - 11");

  write_sol_ASL((ASL*)(my_asl_pointer->asl_pointer),String[0],x,v,0);
  
  if (string_len) FREE(string_len);
  for(i=0;i<m_string*n_string; i++)
    {
      if (String[i]) FREE(String[i]);
    }
  if (String) FREE(String);
  
  return 0;
}

//////////////////////////////////////////////////
// Get informations related to the AMPL problem //
//////////////////////////////////////////////////

int sci_ampl_get_size(char * fname)
{
  int * p_extra = NULL; 
  int * p_pointer = NULL;
  void * asl_pointer = NULL;
  double dtmp;
  SciErr _SciErr;
  static char * ListLabels [] = {"plist",
				 "nbv", /* no. of linear binary variables */
				 "niv", /* no. of linear integer variables */
				 "nlc", /* total no. of nonlinear constraints */
				 "n_eqn", /* number of equality constraints or -1 if unknown (ampl prior to 19970627) */
				 "n_cc", /* total complementarity conditions */
				 "nlcc", /* nonlinear complementarity conditions */
				 "nlnc", /* no. of nonlinear network constraints */
				 "nlo", /* no. of nonlinear objectives */
				 "nlvb", /* no. of nonlinear variables in both constraints and objectives */
				 "nlvc", /* no. of nonlinear variables in constraints */
				 "nlvo", /* no. of nonlinear variables in objectives nlvc_ and nlvo_ include nlvb_ */
				 "nlvbi", /* integer nonlinear variables in both constraints and objectives */
				 "nlvci", /* integer nonlinear vars just in constraints */
				 "nlvoi", /* integer nonlinear vars just in objectives */
				 "nwv", /* no. of (linear) network variables (arcs) */
				 "nzc", /* no. of nonzeros in constraints' Jacobian */
				 "nzo", /* no. of nonzeros in all objective gradients */
				 "n_var", /* total no. of variables */
				 "n_con", /* total no. of constraints */
				 "n_obj", /* total no. of objectives */
				 "n_lcon"}; /* no. of logical constraints */
  struct my_asl * my_asl_pointer = NULL;

  _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_pointer); CHECK_ERROR("ampl_get_size - 1");
  _SciErr = getPointer(pvApiCtx, p_pointer, &asl_pointer); CHECK_ERROR("ampl_get_size - 2");

  my_asl_pointer = (struct my_asl *)asl_pointer;

  if (!my_asl_pointer->asl_pointer->i.filename_)
    {
      sciprint("%s: ampl_init(filename) has not been called\n", fname);
      return 0;
    }

  if ((Lhs!=1)||(Rhs!=1))
    {
      sciprint("%s usage: info = %s(asl)\n", fname, fname);
      return 0;
    }

  _SciErr = createMList(pvApiCtx, Rhs+1, 22, &p_extra); CHECK_ERROR("ampl_get_size - 3");
  _SciErr = createMatrixOfStringInList(pvApiCtx, Rhs + 1, p_extra, 1, 1, 22, ListLabels); CHECK_ERROR("ampl_get_size - 4");

  dtmp = (double)my_asl_pointer->asl_pointer->i.nbv_;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs + 1, p_extra, 2, 1, 1, &dtmp); CHECK_ERROR("ampl_get_size - 5");
  dtmp = (double)my_asl_pointer->asl_pointer->i.niv_;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs + 1, p_extra, 3,  1, 1, &dtmp); CHECK_ERROR("ampl_get_size - 6");
  dtmp = (double)my_asl_pointer->asl_pointer->i.nlc_;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs + 1, p_extra, 4,  1, 1, &dtmp); CHECK_ERROR("ampl_get_size - 7");
  dtmp = (double)my_asl_pointer->asl_pointer->i.n_eqn_;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs + 1, p_extra, 5,  1, 1, &dtmp); CHECK_ERROR("ampl_get_size - 8");
  dtmp = (double)my_asl_pointer->asl_pointer->i.n_cc_;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs + 1, p_extra, 6,  1, 1, &dtmp); CHECK_ERROR("ampl_get_size - 9");
  dtmp = (double)my_asl_pointer->asl_pointer->i.nlcc_;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs + 1, p_extra, 7,  1, 1, &dtmp); CHECK_ERROR("ampl_get_size - 10");
  dtmp = (double)my_asl_pointer->asl_pointer->i.nlnc_;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs + 1, p_extra, 8,  1, 1, &dtmp); CHECK_ERROR("ampl_get_size - 11");
  dtmp = (double)my_asl_pointer->asl_pointer->i.nlo_;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs + 1, p_extra, 9,  1, 1, &dtmp); CHECK_ERROR("ampl_get_size - 12");
  dtmp = (double)my_asl_pointer->asl_pointer->i.nlvb_;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs + 1, p_extra, 10, 1, 1, &dtmp); CHECK_ERROR("ampl_get_size - 13");
  dtmp = (double)my_asl_pointer->asl_pointer->i.nlvc_;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs + 1, p_extra, 11, 1, 1, &dtmp); CHECK_ERROR("ampl_get_size - 14");
  dtmp = (double)my_asl_pointer->asl_pointer->i.nlvo_;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs + 1, p_extra, 12, 1, 1, &dtmp); CHECK_ERROR("ampl_get_size - 15");
  dtmp = (double)my_asl_pointer->asl_pointer->i.nlvbi_;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs + 1, p_extra, 13, 1, 1, &dtmp); CHECK_ERROR("ampl_get_size - 16");
  dtmp = (double)my_asl_pointer->asl_pointer->i.nlvci_;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs + 1, p_extra, 14, 1, 1, &dtmp); CHECK_ERROR("ampl_get_size - 17");
  dtmp = (double)my_asl_pointer->asl_pointer->i.nlvoi_;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs + 1, p_extra, 15, 1, 1, &dtmp); CHECK_ERROR("ampl_get_size - 18");
  dtmp = (double)my_asl_pointer->asl_pointer->i.nwv_;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs + 1, p_extra, 16, 1, 1, &dtmp); CHECK_ERROR("ampl_get_size - 19");
  dtmp = (double)my_asl_pointer->asl_pointer->i.nzc_;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs + 1, p_extra, 17, 1, 1, &dtmp); CHECK_ERROR("ampl_get_size - 20");
  dtmp = (double)my_asl_pointer->asl_pointer->i.nzo_;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs + 1, p_extra, 18, 1, 1, &dtmp); CHECK_ERROR("ampl_get_size - 21");
  dtmp = (double)my_asl_pointer->asl_pointer->i.n_var_;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs + 1, p_extra, 19, 1, 1, &dtmp); CHECK_ERROR("ampl_get_size - 22");
  dtmp = (double)my_asl_pointer->asl_pointer->i.n_con_;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs + 1, p_extra, 20, 1, 1, &dtmp); CHECK_ERROR("ampl_get_size - 23");
  dtmp = (double)my_asl_pointer->asl_pointer->i.n_obj_;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs + 1, p_extra, 21, 1, 1, &dtmp); CHECK_ERROR("ampl_get_size - 24");
  dtmp = (double)my_asl_pointer->asl_pointer->i.n_lcon_;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs + 1, p_extra, 22, 1, 1, &dtmp); CHECK_ERROR("ampl_get_size - 25");

  LhsVar(1) = Rhs + 1;

  return 0;
}

/////////////////////////////////////
// Get the complementarity indexes //
/////////////////////////////////////

int sci_ampl_get_compl(char * fname)
{
  int * p_pointer = NULL;
  void * asl_pointer = NULL;
  int m_cvar, n_cvar;
  int nc, i;
  double * mycvar = NULL;
  struct my_asl * my_asl_pointer = NULL;
  SciErr _SciErr;

  if ((Lhs!=1)||(Rhs!=1))
    {
      sciprint("%s usage: cvar = %s(asl)\n", fname, fname);
      return 0;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_pointer); CHECK_ERROR("ampl_get_compl - 1");
  _SciErr = getPointer(pvApiCtx, p_pointer, &asl_pointer); CHECK_ERROR("ampl_get_compl - 2");

  my_asl_pointer = (struct my_asl *)asl_pointer;

  if (!my_asl_pointer->asl_pointer->i.filename_)
    {
      sciprint("%s: ampl_init(filename) has not been called\n", fname);
      return 0;
    }

  nc = my_asl_pointer->asl_pointer->i.n_con_;

  if (!my_asl_pointer->asl_pointer->i.n_cc_)
    {
      m_cvar = 1; n_cvar = 1;
      mycvar = (double *)MALLOC(m_cvar*n_cvar*sizeof(double));

      mycvar[0] = -1;

      _SciErr = createMatrixOfDouble(pvApiCtx, Rhs+1, m_cvar, n_cvar, mycvar); CHECK_ERROR("ampl_get_compl - 3");

      LhsVar(1) = Rhs + 1;

      if (mycvar) FREE(mycvar);

      return 0;
    }

  /* cvar[i] > 0 means constraint i complements variable cvar[i] - 1 */
  
  m_cvar = nc; n_cvar = 1;
  
  mycvar = (double *)MALLOC(m_cvar*n_cvar*sizeof(double));

  for(i=0;i<nc;i++) mycvar[i] = my_asl_pointer->asl_pointer->i.cvar_[i];

  _SciErr = createMatrixOfDouble(pvApiCtx, Rhs+1, m_cvar, n_cvar, mycvar); CHECK_ERROR("ampl_get_compl - 4");

  LhsVar(1) = Rhs + 1;

  if (mycvar) FREE(mycvar);

  return 0;
}

///////////////////////////////
// Get the type of variables //
///////////////////////////////

int sci_ampl_get_type(char * fname)
{
  int * p_pointer = NULL;
  void * asl_pointer = NULL;
  int m_type, n_type;
  int n, i, i_start, i_end, tmp_int;
  char ** type = NULL;
  struct my_asl * my_asl_pointer = NULL;
  SciErr _SciErr;

  if ((Lhs!=1)||(Rhs!=1))
    {
      sciprint("%s usage: var_type = %s(asl)\n", fname, fname);
      return 0;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_pointer); CHECK_ERROR("ampl_get_type - 1");
  _SciErr = getPointer(pvApiCtx, p_pointer, &asl_pointer); CHECK_ERROR("ampl_get_type - 2");

  my_asl_pointer = (struct my_asl *)asl_pointer;

  if (!my_asl_pointer->asl_pointer->i.filename_)
    {
      sciprint("%s: ampl_init(filename) has not been called\n", fname);
      return 0;
    }
      
  n  = my_asl_pointer->asl_pointer->i.n_var_;

  m_type = n; n_type = 1;

  type    = (char **)MALLOC(sizeof(char *));
  type[0] = (char *)MALLOC((m_type*n_type+1)*sizeof(char));

  tmp_int = (my_asl_pointer->asl_pointer->i.nlvc_<my_asl_pointer->asl_pointer->i.nlvo_)?my_asl_pointer->asl_pointer->i.nlvo_:my_asl_pointer->asl_pointer->i.nlvc_;

  i_start = 0;
  i_end   = i_start + tmp_int;
  
  for(i=i_start; i<i_end;i++) type[0][i] = 'n'; // n: nonlinear variable

  i_start = i_end;
  i_end   = i_start + my_asl_pointer->asl_pointer->i.nwv_;

  for(i=i_start; i<i_end;i++) type[0][i] = 'a'; // a: linear arcs

  i_start = i_end;
  i_end   = i_start + my_asl_pointer->asl_pointer->i.n_var_ - (tmp_int + 
							       my_asl_pointer->asl_pointer->i.niv_ + 
							       my_asl_pointer->asl_pointer->i.nbv_ + 
							       my_asl_pointer->asl_pointer->i.nwv_);

  for(i=i_start; i<i_end;i++) type[0][i] = 'o'; // o: other linear

  i_start = i_end;
  i_end   = i_start + my_asl_pointer->asl_pointer->i.nbv_;

  for(i=i_start; i<i_end;i++) type[0][i] = 'b'; // b: binary variable

  i_start = i_end;
  i_end   = i_start + my_asl_pointer->asl_pointer->i.niv_;

  for(i=i_start; i<i_end;i++) type[0][i] = 'i'; // i: integer variable

  type[0][i_end] = '\0';

#ifdef DEBUG
  DBGPRINTF("DEBUG: i_end = %d n = %d\n", i_end, n);
#endif

  _SciErr = createMatrixOfString(pvApiCtx, Rhs + 1, 1, 1, type); CHECK_ERROR("ampl_get_type - 3");

  LhsVar(1) = Rhs + 1;

  if (type[0]) FREE(type[0]);
  if (type)    FREE(type);

  return 0;
}


void count_tree_leaf_node(expr2 * E, int * nb_node, int * nb_leaf);
void proceed_tree(expr2 * E, double * node_start, double * node_end, double * node_value, char ** node_string, 
		  int * parent_node, int * current_node);

void count_expr_node(expr2_v * E, int * nb_node);
void proceed_expr(expr2_v * E, double * node_index, double * node_value, char ** node_string, int * Index); 
void proceed_var(expr2_v * E, double * node_index, double * node_value, char ** node_string, int * Index); 

int sci_ampl_get_dag(char * fname)
{
  
  int * p_pointer = NULL;
  int * p_extra   = NULL;
  int * p_list    = NULL;
  int * p_objcon  = NULL;
  void * asl_pointer = NULL;
  struct my_asl_fgh * my_asl_pointer = NULL;
  int _n_var, _n_con, _n_obj, nb_node = 0, nb_leaf = 0;
  int i, j, k, size_lin_part, parent_node, current_node, Index;
  int has_lin_part = 0, has_nonlin_part = 0;
  ograd * og;
  SciErr _SciErr;
  double * node_start = NULL, * node_end = NULL, * node_index = NULL;
  double * lin_part_var = NULL;
  double * node_value = NULL, * lin_part_coeff = NULL;
  char ** node_string = NULL;
  static char * ListLabels [] = {"ampl",
				 "obj", // no. of linear binary variables
				 "constr",  // no. of linear integer variables 
                                 "obj_type", // objective type:
						// 0 = constant
						// 1 = linear
						// 2 = quadratic
						// 3 = general nonlinear
				 "constr_type", // constraint type 
				 "var_type"}; // variable type 

  static char * ObjConLabels [] = {"objcon",
				   "start",
				   "end",
				   "value",
				   "name",
				   "lin_part_var",
				   "lin_part_coeff"};
				 
  if ((Lhs!=1)||(Rhs!=1))
    {
      sciprint("%s usage: dag = %s(asl)\n", fname, fname);
      return 0;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, 1, &p_pointer); CHECK_ERROR("ampl_get_dag - 1");
  _SciErr = getPointer(pvApiCtx, p_pointer, &asl_pointer); CHECK_ERROR("ampl_get_dag - 2");

  my_asl_pointer = (struct my_asl_fgh *)asl_pointer;

  if (!my_asl_pointer->asl_pointer->i.filename_)
    {
      sciprint("%s: ampl_init(filename) has not been called\n", fname);
      return 0;
    }
      
  _SciErr = createMList(pvApiCtx, Rhs+1, 6, &p_extra); CHECK_ERROR("ampl_get_dag - 3");
  _SciErr = createMatrixOfStringInList(pvApiCtx, Rhs + 1, p_extra, 1, 1, 6, ListLabels); CHECK_ERROR("ampl_get_dag - 4");

  _n_var = my_asl_pointer->asl_pointer->i.n_var_;
  _n_con = my_asl_pointer->asl_pointer->i.n_con_;
  _n_obj = my_asl_pointer->asl_pointer->i.n_obj_;

#ifdef DEBUG
  sciprint("_n_var = %d, _n_con = %d, _n_obj = %d\n", _n_var, _n_con, _n_obj);
  sciprint("con2_de  = %x\n", my_asl_pointer->asl_pointer->I.con2_de_);
  sciprint("lcon2_de = %x\n", my_asl_pointer->asl_pointer->I.lcon2_de_);
  sciprint("obj2_de  = %x\n", my_asl_pointer->asl_pointer->I.obj2_de_);
  sciprint("var2_e   = %x\n", my_asl_pointer->asl_pointer->I.var2_e_);
  sciprint("c_class  = %s\n", my_asl_pointer->asl_pointer->I.c_class);
  sciprint("o_class  = %s\n", my_asl_pointer->asl_pointer->I.o_class);
  sciprint("v_class  = %s\n", my_asl_pointer->asl_pointer->I.v_class);
#endif

  /////////////////////////////////////
  // Proceed the objective functions //
  /////////////////////////////////////

  // Add a list in the output structure
  _SciErr = createListInList(pvApiCtx, Rhs+1, p_extra, 2, _n_obj, &p_list); CHECK_ERROR("ampl_get_dag - 5");

  for(i=0;i<_n_obj;i++)
    {
      // Proceed the non linear part of the objective functions
      nb_node = 0; nb_leaf = 0;
      if (my_asl_pointer->asl_pointer->I.obj2_de_)
	{
	  if ((my_asl_pointer->asl_pointer->I.obj2_de_ + i)->e)
	    {
	      has_nonlin_part = TRUE;
              // sci_ampl.c(1599) : warning C4133: 'function' : incompatible types - from 'expr2 *' to 'expr2_v *'
              // sci_ampl.c(2217) : warning C4133: 'function' : incompatible types - from 'expr2 *' to 'expr2_v *'
              // sci_ampl.c(2218) : warning C4133: 'function' : incompatible types - from 'expr2 *' to 'expr2_v *'
              // sci_ampl.c(2223) : warning C4133: 'function' : incompatible types - from 'expr2 *' to 'expr2_v *'
              // sci_ampl.c(2253) : warning C4133: 'function' : incompatible types - from 'expr2 *' to 'expr2_v *'
              // sci_ampl.c(2254) : warning C4133: 'function' : incompatible types - from 'expr2 *' to 'expr2_v *'
              // sci_ampl.c(2264) : warning C4133: 'function' : incompatible types - from 'expr2 *' to 'expr2_v *'

	      count_expr_node((my_asl_pointer->asl_pointer->I.obj2_de_ + i)->e, &nb_node);
	      
	      node_start  = (double *)MALLOC((nb_node + nb_leaf)*sizeof(double));
	      node_end    = (double *)MALLOC((nb_node + nb_leaf)*sizeof(double));
	      node_value  = (double *)MALLOC((nb_node + nb_leaf)*sizeof(double));
	      node_string = (char **)MALLOC((nb_node + nb_leaf)*sizeof(char *));
	      
	      parent_node = -1;
	      current_node = 0;
	      proceed_tree((my_asl_pointer->asl_pointer->I.obj2_de_ + i)->e, 
			   node_start, node_end, node_value, node_string, 
			   &parent_node, &current_node);
	    }
	  else
	    {
	      has_nonlin_part = FALSE;
	      node_start   = NULL;
	      node_end     = NULL;
	      node_value   = NULL;
	      node_string  = NULL;
	    }
	}
      else
	{
	  has_nonlin_part = FALSE;
	  node_start   = NULL;
	  node_end     = NULL;
	  node_value   = NULL;
	  node_string  = NULL;
	}

      // Proceed the linear part of the objective functions
      size_lin_part = 0;
      for(og = my_asl_pointer->asl_pointer->i.Ograd_[i]; og; og = og->next) size_lin_part++;

      if (size_lin_part)
	{
	  has_lin_part   = TRUE;
	  lin_part_var   = (double *)MALLOC(size_lin_part*sizeof(double));
	  lin_part_coeff = (double *)MALLOC(size_lin_part*sizeof(double));
	  
	  j = 0;
	  for(og = my_asl_pointer->asl_pointer->i.Ograd_[i]; og; og = og->next)
	    {
	      lin_part_coeff[j] = og->coef;
	      lin_part_var[j]   = og->varno + 1;
	      j++;
	    }
	}
      else
	{
	  has_lin_part   = FALSE;
	  lin_part_var   = NULL;
	  lin_part_coeff = NULL;
	}

      _SciErr = createMListInList(pvApiCtx, Rhs+1, p_list, i+1, 7, &p_objcon); CHECK_ERROR("ampl_get_dag - 6");
      _SciErr = createMatrixOfStringInList(pvApiCtx, Rhs+1, p_objcon, 1, 1, 7, ObjConLabels); CHECK_ERROR("ampl_get_dag - 7");
      if (has_nonlin_part)
	{
	  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, p_objcon, 2, 1, nb_node + nb_leaf, node_start); CHECK_ERROR("ampl_get_dag - 8");
	  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, p_objcon, 3, 1, nb_node + nb_leaf, node_end); CHECK_ERROR("ampl_get_dag - 9");
	  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, p_objcon, 4, 1, nb_node + nb_leaf, node_value); CHECK_ERROR("ampl_get_dag - 10");
	  _SciErr = createMatrixOfStringInList(pvApiCtx, Rhs+1, p_objcon, 5, 1, nb_node + nb_leaf, node_string); CHECK_ERROR("ampl_get_dag - 11");
	}
      else
	{
	  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, p_objcon, 2, 0, 0, NULL); CHECK_ERROR("ampl_get_dag - 12");
	  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, p_objcon, 3, 0, 0, NULL); CHECK_ERROR("ampl_get_dag - 13");
	  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, p_objcon, 4, 0, 0, NULL); CHECK_ERROR("ampl_get_dag - 14");
	  // Empty matrixes are of type double. So, we return a double matrix instead of a string one
	  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, p_objcon, 5, 0, 0, NULL); CHECK_ERROR("ampl_get_dag - 15");
	  node_string = NULL;
	}

      if (has_lin_part)
	{
#ifdef DEBUG
	  sciprint("size_lin_part %d\n", size_lin_part);
	  for(k=0;k<size_lin_part;k++)
	    {
	      sciprint("lin_part_var[%d] = %f\n", k, lin_part_var[k]);
	    }
#endif
	  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, p_objcon, 6, 1, size_lin_part, lin_part_var); CHECK_ERROR("ampl_get_dag - 16");
	  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, p_objcon, 7, 1, size_lin_part, lin_part_coeff); CHECK_ERROR("ampl_get_dag - 17");
	}
      else
	{
	  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, p_objcon, 6, 0, 0, NULL); CHECK_ERROR("ampl_get_dag - 17");
	  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, p_objcon, 7, 0, 0, NULL); CHECK_ERROR("ampl_get_dag - 18");
	}

      if (node_start)     FREE(node_start);
      if (node_end)       FREE(node_end);
      if (node_value)     FREE(node_value);
      if (node_string)    freeArrayOfString(node_string, nb_node + nb_leaf);
      if (lin_part_var)   FREE(lin_part_var);
      if (lin_part_coeff) FREE(lin_part_coeff);
    }

  /////////////////////////////
  // Proceed the constraints //
  /////////////////////////////

  // Add a list in the output structure
  _SciErr = createListInList(pvApiCtx, Rhs+1, p_extra, 3, _n_con, &p_list); CHECK_ERROR("ampl_get_dag - 19");

  for(i=0;i<_n_con;i++)
    {
      // Proceed the nonlinear part of the constraints
      nb_node = 0; nb_leaf = 0;

      if (my_asl_pointer->asl_pointer->I.con2_de_)
	{
	  if ((my_asl_pointer->asl_pointer->I.con2_de_ + i)->e)
	    {
#ifdef DEBUG
	      sciprint("Generic constraints + non linear - %d\n",i);
	      sciprint("Generic constraints + non linear e=%x ee=%x ef=%x\n", 
		       (my_asl_pointer->asl_pointer->I.con2_de_ + i)->e, 
		       (my_asl_pointer->asl_pointer->I.con2_de_ + i)->ee, 
		       (my_asl_pointer->asl_pointer->I.con2_de_ + i)->ef);
#endif
	      has_nonlin_part = TRUE;
	      count_tree_leaf_node((my_asl_pointer->asl_pointer->I.con2_de_ + i)->e, &nb_node, &nb_leaf);
	      
	      node_start  = (double *)MALLOC((nb_node + nb_leaf)*sizeof(double));
	      node_end    = (double *)MALLOC((nb_node + nb_leaf)*sizeof(double));
	      node_value  = (double *)MALLOC((nb_node + nb_leaf)*sizeof(double));
	      node_string = (char **)MALLOC((nb_node + nb_leaf)*sizeof(char *));
	      
	      parent_node = -1;
	      current_node = 0;
	      proceed_tree((my_asl_pointer->asl_pointer->I.con2_de_ + i)->e, 
			   node_start, node_end, node_value, node_string, 
			   &parent_node, &current_node);
	    }
	  else
	    {
#ifdef DEBUG
	      sciprint("Generic constraints + linear- %d\n", i);
	      sciprint("Generic constraints + non linear e=%x ee=%x ef=%x\n", 
		       (my_asl_pointer->asl_pointer->I.con2_de_ + i)->e, 
		       (my_asl_pointer->asl_pointer->I.con2_de_ + i)->ee, 
		       (my_asl_pointer->asl_pointer->I.con2_de_ + i)->ef);
#endif
	      has_nonlin_part = FALSE;
	      node_start  = NULL;
	      node_end    = NULL;
	      node_value  = NULL;
	      node_string = NULL;
	    }
	}
      else if (my_asl_pointer->asl_pointer->I.lcon2_de_)
	{
	  if ((my_asl_pointer->asl_pointer->I.lcon2_de_ + i)->e)
	    {
#ifdef DEBUG
	      sciprint("Logical constraints + non linear - %d\n",i);
	      sciprint("Logical constraints + non linear e=%x ee=%x ef=%x\n", 
		       (my_asl_pointer->asl_pointer->I.lcon2_de_ + i)->e, 
		       (my_asl_pointer->asl_pointer->I.lcon2_de_ + i)->ee, 
		       (my_asl_pointer->asl_pointer->I.lcon2_de_ + i)->ef);
#endif
	      has_nonlin_part = TRUE;
	      count_tree_leaf_node((my_asl_pointer->asl_pointer->I.lcon2_de_ + i)->e, &nb_node, &nb_leaf);
	      
	      node_start  = (double *)MALLOC((nb_node + nb_leaf)*sizeof(double));
	      node_end    = (double *)MALLOC((nb_node + nb_leaf)*sizeof(double));
	      node_value  = (double *)MALLOC((nb_node + nb_leaf)*sizeof(double));
	      node_string = (char **)MALLOC((nb_node + nb_leaf)*sizeof(char *));
	      
	      parent_node = -1;
	      current_node = 0;
	      proceed_tree((my_asl_pointer->asl_pointer->I.lcon2_de_ + i)->e, 
			   node_start, node_end, node_value, node_string, 
			   &parent_node, &current_node);
	    }
	  else
	    {
#ifdef DEBUG
	      sciprint("Logical constraints + linear - %d\n",i);
	      sciprint("Logical constraints + non linear e=%x ee=%x ef=%x\n", 
		       (my_asl_pointer->asl_pointer->I.lcon2_de_ + i)->e, 
		       (my_asl_pointer->asl_pointer->I.lcon2_de_ + i)->ee, 
		       (my_asl_pointer->asl_pointer->I.lcon2_de_ + i)->ef);
#endif
	      has_nonlin_part = FALSE;
	      node_start  = NULL;
	      node_end    = NULL;
	      node_value  = NULL;
	      node_string = NULL;
	    }
	}
      else if (my_asl_pointer->asl_pointer->I.var2_e_)
	{
#ifdef DEBUG
	  sciprint("expressions - %d \n",i);
#endif
	  has_nonlin_part = TRUE;
	  count_expr_node((my_asl_pointer->asl_pointer->I.var2_e_ + i), &nb_node);
	  
	  node_index  = (double *)MALLOC((nb_node + nb_leaf)*sizeof(double));
	  node_value  = (double *)MALLOC((nb_node + nb_leaf)*sizeof(double));
	  node_string = (char **)MALLOC((nb_node + nb_leaf)*sizeof(char *));
	  
	  Index = 0;
	  proceed_var((my_asl_pointer->asl_pointer->I.var2_e_ + i), 
		      node_index, node_value, node_string, &Index);
#ifdef DEBUG
	  sciprint("nb_node = %d nb_leaf = %d\n", nb_node, nb_leaf);
#endif
	}
      else
	{
#ifdef DEBUG
	  sciprint("expression - %d\n", i);
#endif
	  has_nonlin_part = FALSE;
	  node_start  = NULL;
	  node_end    = NULL;
	  node_value  = NULL;
	  node_string = NULL;
	}

      // Proceed the linear part of the constraints
      size_lin_part = 0;
      if (my_asl_pointer->asl_pointer->i.A_vals_)
	{
#ifdef DEBUG
	  sciprint("Linear constraints - %d\n", i);
#endif
	  has_lin_part   = TRUE;

	  for(j=my_asl_pointer->asl_pointer->i.A_colstarts_[i]; j < my_asl_pointer->asl_pointer->i.A_colstarts_[i+1]; j++) size_lin_part++;

	  lin_part_var   = (double *)MALLOC(size_lin_part*sizeof(double));
	  lin_part_coeff = (double *)MALLOC(size_lin_part*sizeof(double));
	  
	  k = 0;
	  for(j=my_asl_pointer->asl_pointer->i.A_colstarts_[i]; j < my_asl_pointer->asl_pointer->i.A_colstarts_[i+1]; j++)
	    {
	      lin_part_var[k]   = my_asl_pointer->asl_pointer->i.A_rownos_[j];
	      lin_part_coeff[k] = my_asl_pointer->asl_pointer->i.A_vals_[j];
	      k++;
	    }
	}
      else
	{
#ifdef DEBUG
	  sciprint("Linear constraints - %d\n", i);
#endif
	  has_lin_part   = FALSE;
	  lin_part_var   = NULL;
	  lin_part_coeff = NULL;
	}

      _SciErr = createMListInList(pvApiCtx, Rhs+1, p_list, i+1, 7, &p_objcon); CHECK_ERROR("ampl_get_dag - 20bis");
      _SciErr = createMatrixOfStringInList(pvApiCtx, Rhs+1, p_objcon, 1, 1, 7, ObjConLabels); CHECK_ERROR("ampl_get_dag - 20");
      if (has_nonlin_part)
	{
	  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, p_objcon, 2, 1, nb_node + nb_leaf, node_start); CHECK_ERROR("ampl_get_dag - 21");
	  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, p_objcon, 3, 1, nb_node + nb_leaf, node_end); CHECK_ERROR("ampl_get_dag - 22");
	  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, p_objcon, 4, 1, nb_node + nb_leaf, node_value); CHECK_ERROR("ampl_get_dag - 23");
#ifdef DEBUG
	  sciprint("nb_node = %d, nd_leaf = %d\n", nb_node, nb_leaf);
	  for(k=0;k<nb_node+nb_leaf;k++)
	    {
	      sciprint("string[%d] = %s\n", k, node_string[k]);
	    }
#endif
	  _SciErr = createMatrixOfStringInList(pvApiCtx, Rhs+1, p_objcon, 5, 1, nb_node + nb_leaf, node_string); CHECK_ERROR("ampl_get_dag - 24");
	}
      else
	{
	  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, p_objcon, 2, 0, 0, NULL); CHECK_ERROR("ampl_get_dag - 25");
	  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, p_objcon, 3, 0, 0, NULL); CHECK_ERROR("ampl_get_dag - 26");
	  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, p_objcon, 4, 0, 0, NULL); CHECK_ERROR("ampl_get_dag - 27");
	  // Empty matrixes are of type double. So, we return a double matrix instead of a string one
	  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, p_objcon, 5, 0, 0, NULL); CHECK_ERROR("ampl_get_dag - 28");
	  node_string = NULL;
	}

      if (has_lin_part)
	{
	  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, p_objcon, 6, 1, size_lin_part, lin_part_var); CHECK_ERROR("ampl_get_dag - 29");
	  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, p_objcon, 7, 1, size_lin_part, lin_part_coeff); CHECK_ERROR("ampl_get_dag - 30");
	}
      else
	{
	  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, p_objcon, 6, 0, 0, NULL); CHECK_ERROR("ampl_get_dag - 31");
	  _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs+1, p_objcon, 7, 0, 0, NULL); CHECK_ERROR("ampl_get_dag - 32");
	}

      if (node_start)     FREE(node_start);
      if (node_end)       FREE(node_end);
      if (node_value)     FREE(node_value);
      if (node_string)    freeArrayOfString(node_string, nb_node + nb_leaf);
      if (lin_part_var)   FREE(lin_part_var);
      if (lin_part_coeff) FREE(lin_part_coeff);
    }

  ////////////////////////////////////////////////////////////////////
  // Proceed the variables, objectives and constraints informations //
  ////////////////////////////////////////////////////////////////////

  // YC: need to pass ASL_find_co_class to pfgh_read. But, once this flag is passed,
  // I still have no informations in *_class.

#ifdef DEBUG
  sciprint("con2_de = %ld\n", my_asl_pointer->asl_pointer->I.con2_de_);
  sciprint("lcon2_de = %ld\n", my_asl_pointer->asl_pointer->I.lcon2_de_);
  sciprint("obj2_de = %ld\n", my_asl_pointer->asl_pointer->I.obj2_de_);
  sciprint("var2_e = %ld\n", my_asl_pointer->asl_pointer->I.var2_e_);
#endif

  // Variables
  if (my_asl_pointer->asl_pointer->I.v_class)
    {
#ifdef DEBUG
      sciprint("v_class is not null.\n");
#endif
      node_string = (char **)MALLOC(1*sizeof(char *));
      node_string[0] = strdup(strcat(my_asl_pointer->asl_pointer->I.v_class,"\0"));
      
      _SciErr = createMatrixOfStringInList(pvApiCtx, Rhs + 1, p_extra, 4, 1, 1, node_string); CHECK_ERROR("ampl_get_dag - 33");
      if (node_string) freeArrayOfString(node_string, 1);
    }
  else
    {
#ifdef DEBUG
      sciprint("v_class is null.\n");
#endif
      _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs + 1, p_extra, 4, 0, 0, NULL); CHECK_ERROR("ampl_get_dag - 34");
    }

  // Objectives
  if (my_asl_pointer->asl_pointer->I.o_class)
    {
#ifdef DEBUG
      sciprint("o_class is not null. o_class_max = %d\n", my_asl_pointer->asl_pointer->I.o_class_max);
#endif
      node_string = (char **)MALLOC(1*sizeof(char *));
      node_string[0] = strdup(strcat(my_asl_pointer->asl_pointer->I.o_class,"\0"));
      
      _SciErr = createMatrixOfStringInList(pvApiCtx, Rhs + 1, p_extra, 5, 1, 1, node_string); CHECK_ERROR("ampl_get_dag - 35");
      if (node_string) freeArrayOfString(node_string, 1);
    }
  else
    {
#ifdef DEBUG
      sciprint("o_class is null. o_class_max = %d\n", my_asl_pointer->asl_pointer->I.o_class_max);
#endif
      _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs + 1, p_extra, 5, 0, 0, NULL); CHECK_ERROR("ampl_get_dag - 36");
    }

  // Constraints
  if (my_asl_pointer->asl_pointer->I.c_class)
    {
#ifdef DEBUG
      sciprint("c_class is not null. c_class_max = %d\n", my_asl_pointer->asl_pointer->I.c_class_max);
#endif
      node_string = (char **)MALLOC(1*sizeof(char *));
      node_string[0] = strdup(strcat(my_asl_pointer->asl_pointer->I.c_class,"\0"));
      
      _SciErr = createMatrixOfStringInList(pvApiCtx, Rhs + 1, p_extra, 6, 1, 1, node_string); CHECK_ERROR("ampl_get_dag - 37");
      if (node_string) freeArrayOfString(node_string, 1);
    }
  else
    {
#ifdef DEBUG
      sciprint("c_class is null. c_class_max = %d\n", my_asl_pointer->asl_pointer->I.c_class_max);
#endif
      _SciErr = createMatrixOfDoubleInList(pvApiCtx, Rhs + 1, p_extra, 6, 0, 0, NULL); CHECK_ERROR("ampl_get_dag - 38");
    }

  LhsVar(1) = Rhs + 1;
  
  return 0;
}

void count_tree_leaf_node(expr2 * E, int * nb_nodes, int * nb_leaf)
{
  long int opnum = (long int)E->op;
  if (opnum == (long int)AMPL_OPNUM)
    {
      (*nb_leaf)++;
    }
  else if (opnum == (long int)AMPL_OPVARVAL)
    {
      (*nb_leaf)++;
    }
  else if (opnum<=(long int)AMPL_OPPOW)
    {
      (*nb_nodes)++;
      if (E->L.e) count_tree_leaf_node(E->L.e, nb_nodes, nb_leaf);
      if (E->R.e) count_tree_leaf_node(E->R.e, nb_nodes, nb_leaf);
    }
  else
    {
      (*nb_nodes)++;
      if (E->L.e) count_tree_leaf_node(E->L.e, nb_nodes, nb_leaf);
    }
}

char * set_operator(long int opnum)
{
  char * string = NULL;

#ifdef DEBUG
  sciprint("AMPL_OPPLUS       = %ld - %ld\n", AMPL_OPPLUS, opnum);
  sciprint("AMPL_OPMINUS      = %ld - %ld\n", AMPL_OPMINUS, opnum);
  sciprint("AMPL_OPMULT       = %ld - %ld\n", AMPL_OPMULT, opnum);
  sciprint("AMPL_OPDIV        = %ld - %ld\n", AMPL_OPDIV, opnum);
  sciprint("AMPL_OPREM        = %ld - %ld\n", AMPL_OPREM, opnum);
  sciprint("AMPL_OPPOW        = %ld - %ld\n", AMPL_OPPOW, opnum);
  sciprint("AMPL_OPLESS       = %ld - %ld\n", AMPL_OPLESS, opnum);
  sciprint("AMPL_MINLIST      = %ld - %ld\n", AMPL_MINLIST, opnum);
  sciprint("AMPL_MAXLIST      = %ld - %ld\n", AMPL_MAXLIST, opnum);
  sciprint("AMPL_FLOOR        = %ld - %ld\n", AMPL_FLOOR, opnum);
  sciprint("AMPL_CEIL         = %ld - %ld\n", AMPL_CEIL, opnum);
  sciprint("AMPL_ABS          = %ld - %ld\n", AMPL_ABS, opnum);
  sciprint("AMPL_OPUMINUS     = %ld - %ld\n", AMPL_OPUMINUS, opnum);
  sciprint("AMPL_OPOR         = %ld - %ld\n", AMPL_OPOR, opnum);
  sciprint("AMPL_OPAND        = %ld - %ld\n", AMPL_OPAND, opnum);
  sciprint("AMPL_LT           = %ld - %ld\n", AMPL_LT, opnum);
  sciprint("AMPL_LE           = %ld - %ld\n", AMPL_LE, opnum);
  sciprint("AMPL_EQ           = %ld - %ld\n", AMPL_EQ, opnum);
  sciprint("AMPL_GE           = %ld - %ld\n", AMPL_GE, opnum);
  sciprint("AMPL_GT           = %ld - %ld\n", AMPL_GT, opnum);
  sciprint("AMPL_NE           = %ld - %ld\n", AMPL_NE, opnum);
  sciprint("AMPL_OPNOT        = %ld - %ld\n", AMPL_OPNOT, opnum);
  sciprint("AMPL_OPIFnl       = %ld - %ld\n", AMPL_OPIFnl, opnum);
  sciprint("AMPL_OP_tanh      = %ld - %ld\n", AMPL_OP_tanh, opnum);
  sciprint("AMPL_OP_tan       = %ld - %ld\n", AMPL_OP_tan, opnum);
  sciprint("AMPL_OP_sqrt      = %ld - %ld\n", AMPL_OP_sqrt, opnum);
  sciprint("AMPL_OP_sinh      = %ld - %ld\n", AMPL_OP_sinh, opnum);
  sciprint("AMPL_OP_sin       = %ld - %ld\n", AMPL_OP_sin, opnum);
  sciprint("AMPL_OP_log10     = %ld - %ld\n", AMPL_OP_log10, opnum);
  sciprint("AMPL_OP_log       = %ld - %ld\n", AMPL_OP_log, opnum);
  sciprint("AMPL_OP_exp       = %ld - %ld\n", AMPL_OP_exp, opnum);
  sciprint("AMPL_OP_cosh      = %ld - %ld\n", AMPL_OP_cosh, opnum);
  sciprint("AMPL_OP_cos       = %ld - %ld\n", AMPL_OP_cos, opnum);
  sciprint("AMPL_OP_atanh     = %ld - %ld\n", AMPL_OP_atanh, opnum);
  sciprint("AMPL_OP_atan2     = %ld - %ld\n", AMPL_OP_atan2, opnum);
  sciprint("AMPL_OP_atan      = %ld - %ld\n", AMPL_OP_atan, opnum);
  sciprint("AMPL_OP_asinh     = %ld - %ld\n", AMPL_OP_asinh, opnum);
  sciprint("AMPL_OP_asin      = %ld - %ld\n", AMPL_OP_asin, opnum);
  sciprint("AMPL_OP_acosh     = %ld - %ld\n", AMPL_OP_acosh, opnum);
  sciprint("AMPL_OP_acos      = %ld - %ld\n", AMPL_OP_acos, opnum);
  sciprint("AMPL_OPSUMLIST    = %ld - %ld\n", AMPL_OPSUMLIST, opnum);
  sciprint("AMPL_OPintDIV     = %ld - %ld\n", AMPL_OPintDIV, opnum);
  sciprint("AMPL_OPprecision  = %ld - %ld\n", AMPL_OPprecision, opnum);
  sciprint("AMPL_OPround      = %ld - %ld\n", AMPL_OPround, opnum);
  sciprint("AMPL_OPtrunc      = %ld - %ld\n", AMPL_OPtrunc, opnum);
  sciprint("AMPL_OPCOUNT      = %ld - %ld\n", AMPL_OPCOUNT, opnum);
  sciprint("AMPL_OPNUMBEROF   = %ld - %ld\n", AMPL_OPNUMBEROF, opnum);
  sciprint("AMPL_OPNUMBEROFs  = %ld - %ld\n", AMPL_OPNUMBEROFs, opnum);
  sciprint("AMPL_OPATLEAST    = %ld - %ld\n", AMPL_OPATLEAST, opnum);
  sciprint("AMPL_OPATMOST     = %ld - %ld\n", AMPL_OPATMOST, opnum);
  sciprint("AMPL_OPPLTERM     = %ld - %ld\n", AMPL_OPPLTERM, opnum);
  sciprint("AMPL_OPIFSYM      = %ld - %ld\n", AMPL_OPIFSYM, opnum);
  sciprint("AMPL_OPEXACTLY    = %ld - %ld\n", AMPL_OPEXACTLY, opnum);
  sciprint("AMPL_OPNOTATLEAST = %ld - %ld\n", AMPL_OPNOTATLEAST, opnum);
  sciprint("AMPL_OPNOTATMOST  = %ld - %ld\n", AMPL_OPNOTATMOST, opnum);
  sciprint("AMPL_OPNOTEXACTLY = %ld - %ld\n", AMPL_OPNOTEXACTLY, opnum);
  sciprint("AMPL_ANDLIST      = %ld - %ld\n", AMPL_ANDLIST, opnum);
  sciprint("AMPL_ORLIST       = %ld - %ld\n", AMPL_ORLIST, opnum);
  sciprint("AMPL_OPIMPELSE    = %ld - %ld\n", AMPL_OPIMPELSE, opnum);
  sciprint("AMPL_OP_IFF       = %ld - %ld\n", AMPL_OP_IFF, opnum);
  sciprint("AMPL_OPALLDIFF    = %ld - %ld\n", AMPL_OPALLDIFF, opnum);
  sciprint("AMPL_OP1POW       = %ld - %ld\n", AMPL_OP1POW, opnum);
  sciprint("AMPL_OP2POW       = %ld - %ld\n", AMPL_OP2POW, opnum);
  sciprint("AMPL_OPCPOW       = %ld - %ld\n", AMPL_OPCPOW, opnum);
  sciprint("AMPL_OPFUNCALL    = %ld - %ld\n", AMPL_OPFUNCALL, opnum);
  sciprint("AMPL_OPNUM        = %ld - %ld\n", AMPL_OPNUM, opnum);
  sciprint("AMPL_OPHOL        = %ld - %ld\n", AMPL_OPHOL, opnum);
  sciprint("AMPL_OPVARVAL     = %ld - %ld\n", AMPL_OPVARVAL, opnum);
#endif

  if      (opnum == (long int)AMPL_OPPLUS)       string = strdup("OPPLUS");
  else if (opnum == (long int)AMPL_OPMINUS)      string = strdup("OPMINUS");
  else if (opnum == (long int)AMPL_OPMULT)       string = strdup("OPMULT");
  else if (opnum == (long int)AMPL_OPDIV)        string = strdup("OPDIV");
  else if (opnum == (long int)AMPL_OPREM)        string = strdup("PREM");
  else if (opnum == (long int)AMPL_OPPOW)        string = strdup("OPPOW");
  else if (opnum == (long int)AMPL_OPLESS)       string = strdup("OPLESS");
  else if (opnum == (long int)AMPL_MINLIST)      string = strdup("MINLIST");
  else if (opnum == (long int)AMPL_MAXLIST)      string = strdup("MAXLIST");
  else if (opnum == (long int)AMPL_FLOOR)        string = strdup("FLOOR");
  else if (opnum == (long int)AMPL_CEIL)         string = strdup("CEIL");
  else if (opnum == (long int)AMPL_ABS)          string = strdup("ABS");
  else if (opnum == (long int)AMPL_OPUMINUS)     string = strdup("OPUMINUS");
  else if (opnum == (long int)AMPL_OPOR)         string = strdup("OPOR");
  else if (opnum == (long int)AMPL_OPAND)        string = strdup("OPAND");
  else if (opnum == (long int)AMPL_LT)           string = strdup("LT");
  else if (opnum == (long int)AMPL_LE)           string = strdup("LE");
  else if (opnum == (long int)AMPL_EQ)           string = strdup("EQ");
  else if (opnum == (long int)AMPL_GE)           string = strdup("GE");
  else if (opnum == (long int)AMPL_GT)           string = strdup("GT");
  else if (opnum == (long int)AMPL_NE)           string = strdup("NE");
  else if (opnum == (long int)AMPL_OPNOT)        string = strdup("OPNOT");
  else if (opnum == (long int)AMPL_OPIFnl)       string = strdup("OPIFnl");
  else if (opnum == (long int)AMPL_OP_tanh)      string = strdup("OP_tanh");
  else if (opnum == (long int)AMPL_OP_tan)       string = strdup("OP_tan");
  else if (opnum == (long int)AMPL_OP_sqrt)      string = strdup("OP_sqrt");
  else if (opnum == (long int)AMPL_OP_sinh)      string = strdup("OP_sinh");
  else if (opnum == (long int)AMPL_OP_sin)       string = strdup("OP_sin");
  else if (opnum == (long int)AMPL_OP_log10)     string = strdup("OP_log10");
  else if (opnum == (long int)AMPL_OP_log)       string = strdup("OP_log");
  else if (opnum == (long int)AMPL_OP_exp)       string = strdup("OP_exp");
  else if (opnum == (long int)AMPL_OP_cosh)      string = strdup("OP_cosh");
  else if (opnum == (long int)AMPL_OP_cos)       string = strdup("OP_cos");
  else if (opnum == (long int)AMPL_OP_atanh)     string = strdup("OP_atanh");
  else if (opnum == (long int)AMPL_OP_atan2)     string = strdup("OP_atan2");
  else if (opnum == (long int)AMPL_OP_atan)      string = strdup("OP_atan");
  else if (opnum == (long int)AMPL_OP_asinh)     string = strdup("OP_asinh");
  else if (opnum == (long int)AMPL_OP_asin)      string = strdup("OP_asin");
  else if (opnum == (long int)AMPL_OP_acosh)     string = strdup("OP_acosh");
  else if (opnum == (long int)AMPL_OP_acos)      string = strdup("OP_acos");
  else if (opnum == (long int)AMPL_OPSUMLIST)    string = strdup("OPSUMLIST");
  else if (opnum == (long int)AMPL_OPintDIV)     string = strdup("OPIntDIV");
  else if (opnum == (long int)AMPL_OPprecision)  string = strdup("OPprecision");
  else if (opnum == (long int)AMPL_OPround)      string = strdup("OPround");
  else if (opnum == (long int)AMPL_OPtrunc)      string = strdup("OPtrunc");
  else if (opnum == (long int)AMPL_OPCOUNT)      string = strdup("OPCOUNT");
  else if (opnum == (long int)AMPL_OPNUMBEROF)   string = strdup("OPNUMBEROF");
  else if (opnum == (long int)AMPL_OPNUMBEROFs)  string = strdup("OPNUMBEROFs");
  else if (opnum == (long int)AMPL_OPATLEAST)    string = strdup("OPATLEAST");
  else if (opnum == (long int)AMPL_OPATMOST)     string = strdup("OPATMOST");
  else if (opnum == (long int)AMPL_OPPLTERM)     string = strdup("OPPLTERM");
  else if (opnum == (long int)AMPL_OPIFSYM)      string = strdup("OPIFSYM");
  else if (opnum == (long int)AMPL_OPEXACTLY)    string = strdup("OPEXACTLY");
  else if (opnum == (long int)AMPL_OPNOTATLEAST) string = strdup("OPNOTATLEAST");
  else if (opnum == (long int)AMPL_OPNOTATMOST)  string = strdup("OPNOTATMOST");
  else if (opnum == (long int)AMPL_OPNOTEXACTLY) string = strdup("OPNOTEXACTLY");
  else if (opnum == (long int)AMPL_ANDLIST)      string = strdup("ANDLIST");
  else if (opnum == (long int)AMPL_ORLIST)       string = strdup("ORLIST");
  else if (opnum == (long int)AMPL_OPIMPELSE)    string = strdup("OPIMPELSE");
  else if (opnum == (long int)AMPL_OP_IFF)       string = strdup("OP_IFF");
  else if (opnum == (long int)AMPL_OPALLDIFF)    string = strdup("OPALLDIFF");
  else if (opnum == (long int)AMPL_OP1POW)       string = strdup("OP1POW");
  else if (opnum == (long int)AMPL_OP2POW)       string = strdup("OP2POW");
  else if (opnum == (long int)AMPL_OPCPOW)       string = strdup("OPCPOW");
  else if (opnum == (long int)AMPL_OPFUNCALL)    string = strdup("OPFUNCALL");
  else if (opnum == (long int)AMPL_OPNUM)        string = strdup("OPNUM");
  else if (opnum == (long int)AMPL_OPHOL)        string = strdup("OPHOL");
  else if (opnum == (long int)AMPL_OPVARVAL)     string = strdup("OPVARVAL");
  else
    string = strdup("OPUNKNOW");

  return string;
}


void proceed_tree(expr2 * E, double * node_start, double * node_end, double * node_value, char ** node_string, 
		  int * parent_node, int * current_node)
{
  long int opnum = (long int)E->op;
  
  if (opnum == (long int)AMPL_OPNUM)
    {
      node_start[*current_node]  = *parent_node;
      node_end[*current_node]    = *current_node;
      node_value[*current_node]  = ((expr_n *)E)->v;
      node_string[*current_node] = set_operator(opnum);
      *current_node = *current_node + 1;
    }
  else if (opnum == (long int)AMPL_OPVARVAL)
    {
      node_start[*current_node]  = *parent_node;
      node_end[*current_node]    = *current_node;
      node_value[*current_node]  = E->a; // the ath variable
      node_string[*current_node] = set_operator(opnum);
      *current_node = *current_node + 1;
    }
  else if (opnum<=(long int)AMPL_OPPOW)
    {
      // Binary operators
      node_start[*current_node]  = *parent_node;
      node_end[*current_node]    = *current_node;
      node_value[*current_node]  = -1;
      node_string[*current_node] = set_operator(opnum);
      *parent_node = *current_node;
      *current_node = *current_node + 1;

      proceed_tree(E->L.e, node_start, node_end, node_value, node_string, parent_node, current_node);
      proceed_tree(E->R.e, node_start, node_end, node_value, node_string, parent_node, current_node);
    }
  else
    {
      // Unary operators
      node_start[*current_node]  = *parent_node;
      node_end[*current_node]    = *current_node;
      node_value[*current_node]  = -1;
      node_string[*current_node] = set_operator(opnum);
      *parent_node = *current_node;
      *current_node = *current_node + 1;

      proceed_tree(E->L.e, node_start, node_end, node_value, node_string, parent_node, current_node);
    }
}

void count_expr_node(expr2_v * E, int * nb_node)
{
  long int opnum = (long int)E->op;
  if (opnum == (long int)AMPL_OPNUM)
    {
      (*nb_node)++;
    }
  else if (opnum == (long int)AMPL_OPVARVAL)
    {
      (*nb_node)++;
    }
  else if (opnum<=(long int)AMPL_OPPOW)
    {
      (*nb_node)++;
      if (E->bak) count_expr_node(E->bak, nb_node);
      if (E->fwd) count_expr_node(E->fwd, nb_node);
    }
  else
    {
      (*nb_node)++;
      if (E->fwd) count_expr_node(E->fwd, nb_node);
    }
}

void proceed_expr(expr2_v * E, double * node_index, double * node_value, char ** node_string, int * Index)
{
  long int opnum = (long int)E->op;
  
  if (opnum == (long int)AMPL_OPNUM)
    {
      node_index[*Index]  = E->a;
      node_value[*Index]  = E->v;
      node_string[*Index] = set_operator(opnum);
      *Index = *Index + 1;
    }
  else if (opnum == (long int)AMPL_OPVARVAL)
    {
      node_index[*Index]  = E->a;
      node_value[*Index]  = E->v;
      node_string[*Index] = set_operator(opnum);
      *Index = *Index + 1;
    }
  else if (opnum<=(long int)AMPL_OPPOW)
    {
      // Binary operators
      node_index[*Index]  = E->a;
      node_value[*Index]  = E->v;
      node_string[*Index] = set_operator(opnum);
      *Index = *Index + 1;

      proceed_expr(E->bak, node_index, node_value, node_string, Index);
      proceed_expr(E->fwd, node_index, node_value, node_string, Index);
    }
  else
    {
      // Unary operators
      node_index[*Index]  = E->a;
      node_value[*Index]  = E->v;
      node_string[*Index] = set_operator(opnum);
      *Index = *Index + 1;

      proceed_expr(E->fwd, node_index, node_value, node_string, Index);
    }
}

void proceed_var(expr2_v * E, double * node_index, double * node_value, char ** node_string, int * Index)
{
  long int opnum = (long int)E->op;
  
  node_index[*Index]  = E->a;
  node_value[*Index]  = E->v;
  node_string[*Index] = strdup("OPVAR");
  *Index = *Index + 1;
}
