//
//  Copyright (C) 2008-2010 Yann Collette.
//  Copyright (C) 2001-2007 Nicolo' Giorgetti.
//
//  SCIGLPK is free software; you can redistribute it and/or modify it
//  under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 2, or (at your option) any
//  later version.
//
//  This part of code is distributed with the FURTHER condition that it 
//  can be compiled and linked with the Matlab libraries and it can be 
//  used within the Matlab environment.
//
//  SCIGLPK is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with SciCoinOR; see the file COPYING.  If not, write to the Free
//  Software Foundation, 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//

#include <ctime>
#include <limits>
#include <cstdio>

extern "C"
{
#include <stack-c.h>
#include <sciprint.h>
#include <Scierror.h>
#include <MALLOC.h>
#include <lpx.h>
#include <glpk.h>
#include <api_parameters.h>
}

#include <string.h>

//#define DEBUG 1

#include <api_scilab.h>

#include <helper.hpp>

using namespace std;

//#ifdef WIN32
//template<typename T> bool isinf(T value)
//{
//  return numeric_limits<T>::has_infinity && value == numeric_limits<T>::infinity();
//}
//#endif

static int glpk_print_hook(void *info, const char *msg)
{
  sciprint("%s",msg);
  return 1;
}

// This callback allows glpk to terminate when ctrl+c is hit
void cb_func(glp_tree * tree, void * info)
{
  if (C2F(basbrk).iflag==-1)
    {
      C2F(basbrk).iflag=0;
      glp_ios_terminate(tree);
    }
}

// [xmin, fmin, status, extra] = sciglpk(c, a, lhs, rhs, lb, ub, ctype, vartype, param);

// Input arguments
#define	C_IN	     1
#define	A_IN	     2
#define	LHS_IN	     3
#define	RHS_IN	     4
#define LB_IN	     5
#define UB_IN        6
#define CTYPE_IN     7
#define VARTYPE_IN   8
#define PARAM_IN     9
#define LASTPARAM    9

// Output Arguments
#define	 XMIN_OUT     10
#define	 FMIN_OUT     11
#define	 STATUS_OUT   12
#define  EXTRA_OUT    13

using namespace std;

extern "C" int sciglpk(char * fname)
{
  int m_c,   n_c,   * c_addr = NULL;
  int m_a,   n_a,   * a_addr = NULL;
  int m_lhs, n_lhs, * lhs_addr = NULL;
  int m_rhs, n_rhs, * rhs_addr = NULL;
  int m_lb,  n_lb,  * lb_addr = NULL;
  int m_ub,  n_ub,  * ub_addr = NULL;
  int * ctype_addr = NULL;
  int * vartype_addr = NULL;
  int * extra_addr = NULL;
  double * c = NULL, * a = NULL, * lhs = NULL, * rhs = NULL; 
  double * lb = NULL, * ub = NULL;
  char * ctype = NULL, * vartype = NULL;
  int Log = 0;
  int i, j, n_var, n_constr;
  int size, isMIP = 0, msglev = 0, isRelaxed = 0;
  SciErr _SciErr;

  if (strcmp(glp_version(),"4.33")<0)
    {
      Scierror(999,"%s: This Scilab interface is compatible only with GLPK version 4.33 or higher.\n",fname);
      return 0;
    }
  
  if (Rhs != LASTPARAM)
    {
      sciprint("Scilab interface to GLPK Version %s\n",glp_version());
      sciprint("Internal interface for the GNU GLPK library.\n");
      sciprint("You should use the 'glpk' function instead.\n\n");
      sciprint("SYNTAX: [xopt, fmin, status, extra] = sciglpk(c, a, b, lb, ub, ctype, vartype, param)\n");
      return 0;
    }
  
  //
  //  1nd Input. A column array containing the objective function coefficients 
  //

  _SciErr = getVarAddressFromPosition(pvApiCtx, C_IN, &c_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, c_addr, &m_c, &n_c, &c); SCICOINOR_ERROR;
  
  n_var = n_c*m_c;

#ifdef DEBUG
  sciprint("DEBUG m_c = %d n_c = %d\n", m_c, n_c);
  for(i=0;i<m_c*n_c;i++) sciprint("c[%d] = %f\n",i,*(c+i)); 
#endif

  //
  //  2nd Input. A matrix containing the constraints coefficients. If matrix A is NOT a sparse matrix 
  //

  int    *glpk_rn = NULL;
  int    *glpk_cn = NULL;
  double *glpk_a  = NULL;
  int nz = 0, type;

  _SciErr = getVarAddressFromPosition(pvApiCtx, A_IN, &a_addr); SCICOINOR_ERROR;
  _SciErr = getVarType(pvApiCtx, a_addr, &type); SCICOINOR_ERROR;
  
  if(type!=sci_sparse)
    {
#ifdef DEBUG
      sciprint("DEBUG: A_IN is not sparse\n");
#endif
      
      _SciErr = getMatrixOfDouble(pvApiCtx, a_addr, &m_a, &n_a, &a); SCICOINOR_ERROR;
      
      if (a==NULL) 
	{
	  Scierror(999,"%s: invalid value of matrix a\n",fname);
	  return 0;
	}
      
      glpk_rn = new int[m_a*n_a+2];
      glpk_cn = new int[m_a*n_a+2];
      glpk_a  = new double[m_a*n_a+2];
      
      for(i=0; i<m_a; i++)
	{
	  for(j=0; j<n_a; j++)
	    {
	      if (*(a+i+j*m_a)!=0)
		{
		  nz++;
		  glpk_rn[nz] = i+1;
		  glpk_cn[nz] = j+1;
		  glpk_a[nz]  = *(a+i+j*m_a);
		}
	    }
	}
    }
  else
    {
#ifdef DEBUG
      sciprint("DEBUG: A_IN is sparse\n");
#endif
      
      SciSparse S;
      
      getAllocatedSparseMatrix(pvApiCtx, a_addr, &S.m, &S.n, &S.nel, &S.mnel, &S.icol, &S.R);
      
      nz = S.nel;
      
      glpk_rn = new int[nz+1];
      glpk_cn = new int[nz+1];
      glpk_a  = new double[nz+1];
      
      int count = 0;
      for(i=0;i<S.m;i++)
	{
	  if (S.mnel[i]!=0) 
	    {
	      for(j=0;j<S.mnel[i];j++)
		{
		  count++;
		  glpk_rn[count] = i+1;
		  glpk_cn[count] = S.icol[count-1];
		  glpk_a[count]  = S.R[count-1];
		}
	    }
	}
      
      freeAllocatedSparseMatrix(S.mnel, S.icol, S.R);
    }
  
#ifdef DEBUG
  sciprint("DEBUG A - nz = %d\n", nz);
  for(i=0;i<nz;i++) sciprint("a[%d] = %f - cn[%d] = %d - rn[%d] = %d\n",i+1,glpk_a[i+1],i+1,glpk_cn[i+1],i+1,glpk_rn[i+1]);
#endif

  //
  // 3rd Input. A column array containing the left-hand side value for each constraint in the constraint matrix 
  //

  _SciErr = getVarAddressFromPosition(pvApiCtx, LHS_IN, &lhs_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, lhs_addr, &m_lhs, &n_lhs, &lhs);  SCICOINOR_ERROR;

  n_constr = m_lhs*n_lhs;

  //
  // 4th Input. A column array containing the right-hand side value for each constraint in the constraint matrix 
  //
  
  _SciErr = getVarAddressFromPosition(pvApiCtx, RHS_IN, &rhs_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, rhs_addr, &m_rhs, &n_rhs, &rhs);  SCICOINOR_ERROR;

  if (n_constr==0) n_constr = m_rhs*n_rhs;

#ifdef DEBUG
  sciprint("DEBUG m_b = %d n_b = %d\n", m_rhs, n_rhs);
  for(i=0;i<n_constr;i++) sciprint("b[%d] = %f\n",i,*(rhs+i));
#endif

  //
  //  5th Input. An array of length n_var containing the lower bound on each of the variables 
  //

  _SciErr = getVarAddressFromPosition(pvApiCtx, LB_IN, &lb_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, lb_addr, &m_lb, &n_lb, &lb);  SCICOINOR_ERROR;

#ifdef DEBUG
  sciprint("DEBUG m_lb = %d n_lb = %d n_var = %d\n", m_lb, n_lb, n_var);
  for(i=0;i<n_var;i++) sciprint("lb[%d] = %f\n",i,*(lb+i)); 
#endif

  //
  // 6th Input. An array of at least length numcols containing the upper bound on each of the variables
  //

  _SciErr = getVarAddressFromPosition(pvApiCtx, UB_IN, &ub_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, ub_addr, &m_ub, &n_ub, &ub);  SCICOINOR_ERROR;

#ifdef DEBUG
  sciprint("DEBUG m_ub = %d n_ub = %d n_var = %d \n", m_ub, n_ub, n_var);
  for(i=0;i<n_var;i++) sciprint("ub[%d] = %f\n",i,*(ub+i));
#endif

  //
  // 7th Input. A column array containing the sense of each constraint in the constraint matrix 
  //

  _SciErr = getVarAddressFromPosition(pvApiCtx, CTYPE_IN, &ctype_addr); SCICOINOR_ERROR;
  getAllocatedSingleString(pvApiCtx, ctype_addr, &ctype);

#ifdef DEBUG
  sciprint("DEBUG: length = %d string = %s\n", strlen(ctype), ctype);
  sciprint("DEBUG: n_constr = %d\n", n_constr);
#endif

  //
  // 8th Input. A column array containing the types of the variables
  //

  _SciErr = getVarAddressFromPosition(pvApiCtx, VARTYPE_IN, &vartype_addr); SCICOINOR_ERROR;
  getAllocatedSingleString(pvApiCtx, vartype_addr, &vartype);

  size = strlen(vartype);
  if (size==0) 
    {
      Scierror(999,"%s: invalid value of vartype\n",fname);

      freeAllocatedSingleString(ctype);
      freeAllocatedSingleString(vartype);

      if (glpk_rn) delete [] glpk_rn;
      if (glpk_cn) delete [] glpk_cn;
      if (glpk_a)  delete [] glpk_a;
		     
      return 0;
    }

  //
  // 9th Input. A structure containing the control parameters 
  //
 
  int typx = 0;

  clock_t t_start = clock();

  // Create an empty LP/MILP object 
  glp_prob * lp = lpx_create_prob();

  int * param_in_addr = NULL;
  int tmp_res, tmp_int;
  char * tmp_char = NULL;
  double tmp_double;

  initPList(pvApiCtx, PARAM_IN, &param_in_addr);
  if (!checkPList(pvApiCtx, param_in_addr))
    {
      Scierror(999, "%s: argument n° %d is not a plist\n", fname, PARAM_IN);

      freeAllocatedSingleString(ctype);
      freeAllocatedSingleString(vartype);

      if (glpk_rn) delete [] glpk_rn;
      if (glpk_cn) delete [] glpk_cn;
      if (glpk_a)  delete [] glpk_a;

      return 0;
    }
  
  getIntInPList(pvApiCtx, param_in_addr, "msglev", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  lpx_set_int_parm(lp,  LPX_K_MSGLEV, tmp_int);
  msglev = tmp_int;

  // Redirect standard output
  if (msglev > 1) glp_term_hook(glpk_print_hook, NULL);
  else            glp_term_hook(NULL, NULL);
  
  getIntInPList(pvApiCtx, param_in_addr, "scale", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) lpx_set_int_parm(lp,  LPX_K_SCALE, tmp_int);

  getIntInPList(pvApiCtx, param_in_addr, "dual", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) lpx_set_int_parm(lp,  LPX_K_DUAL, tmp_int);

  getIntInPList(pvApiCtx, param_in_addr, "price", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) lpx_set_int_parm(lp,  LPX_K_PRICE, tmp_int);

  getIntInPList(pvApiCtx, param_in_addr, "round", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) lpx_set_int_parm(lp,  LPX_K_ROUND, tmp_int);

  getIntInPList(pvApiCtx, param_in_addr, "itlim", &tmp_int, &tmp_res, 10000, Log, CHECK_NONE);
  if (tmp_res!=-1) lpx_set_int_parm(lp,  LPX_K_ITLIM, tmp_int);

  getIntInPList(pvApiCtx, param_in_addr, "outfrq", &tmp_int, &tmp_res, 100, Log, CHECK_NONE);
  if (tmp_res!=-1) lpx_set_int_parm(lp,  LPX_K_OUTFRQ, tmp_int);

  getIntInPList(pvApiCtx, param_in_addr, "mpsinfo", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) lpx_set_int_parm(lp,  LPX_K_MPSINFO, tmp_int);

  getIntInPList(pvApiCtx, param_in_addr, "mpsobj", &tmp_int, &tmp_res, 2, Log, CHECK_NONE);
  if (tmp_res!=-1) lpx_set_int_parm(lp,  LPX_K_MPSOBJ, tmp_int);

  getIntInPList(pvApiCtx, param_in_addr, "mpsorig", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) lpx_set_int_parm(lp,  LPX_K_MPSORIG, tmp_int);

  getIntInPList(pvApiCtx, param_in_addr, "mpswide", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) lpx_set_int_parm(lp,  LPX_K_MPSWIDE, tmp_int);

  getIntInPList(pvApiCtx, param_in_addr, "mpsfree", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) lpx_set_int_parm(lp,  LPX_K_MPSFREE, tmp_int);

  getIntInPList(pvApiCtx, param_in_addr, "mpsskip", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) lpx_set_int_parm(lp,  LPX_K_MPSSKIP, tmp_int);

  getIntInPList(pvApiCtx, param_in_addr, "branch", &tmp_int, &tmp_res, 2, Log, CHECK_NONE);
  if (tmp_res!=-1) lpx_set_int_parm(lp,  LPX_K_BRANCH, tmp_int);

  getIntInPList(pvApiCtx, param_in_addr, "btrack", &tmp_int, &tmp_res, 2, Log, CHECK_NONE);
  if (tmp_res!=-1) lpx_set_int_parm(lp,  LPX_K_BTRACK, tmp_int);

  getIntInPList(pvApiCtx, param_in_addr, "presol", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) lpx_set_int_parm(lp,  LPX_K_PRESOL, tmp_int);

  getIntInPList(pvApiCtx, param_in_addr, "usecuts", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) lpx_set_int_parm(lp,  LPX_K_USECUTS, tmp_int);

  // Set the sense of optimization
  getIntInPList(pvApiCtx, param_in_addr, "sense", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_int == 1) glp_set_obj_dir(lp, GLP_MIN);
  else                 glp_set_obj_dir(lp, GLP_MAX);

  getDoubleInPList(pvApiCtx, param_in_addr, "relax", &tmp_double, &tmp_res, 0.07, Log, CHECK_NONE);
  if (tmp_res!=-1) lpx_set_real_parm(lp, LPX_K_RELAX, tmp_double);

  // Relative tolerance used to check if the current basic solution is primal feasible
  getDoubleInPList(pvApiCtx, param_in_addr, "tolbnd", &tmp_double, &tmp_res, 1e-7, Log, CHECK_NONE);
  if (tmp_res!=-1) lpx_set_real_parm(lp, LPX_K_TOLBND, tmp_double);

  // Absolute tolerance used to check if the current basic solution is dual feasible
  getDoubleInPList(pvApiCtx, param_in_addr, "toldj", &tmp_double, &tmp_res, 1e-7, Log, CHECK_NONE);
  if (tmp_res!=-1) lpx_set_real_parm(lp, LPX_K_TOLDJ, tmp_double);

  // Relative tolerance used to choose eligible pivotal elements of the simplex table in the ratio test
  getDoubleInPList(pvApiCtx, param_in_addr, "tolpiv", &tmp_double, &tmp_res, 1e-9, Log, CHECK_NONE);
  if (tmp_res!=-1) lpx_set_real_parm(lp, LPX_K_TOLPIV, tmp_double);

  getDoubleInPList(pvApiCtx, param_in_addr, "objll", &tmp_double, &tmp_res, numeric_limits<double>::min(), Log, CHECK_NONE);
  if (tmp_res!=-1) lpx_set_real_parm(lp, LPX_K_OBJLL, tmp_double);

  getDoubleInPList(pvApiCtx, param_in_addr, "objul", &tmp_double, &tmp_res, numeric_limits<double>::max(), Log, CHECK_NONE);
  if (tmp_res!=-1) lpx_set_real_parm(lp, LPX_K_OBJUL, tmp_double);

  getDoubleInPList(pvApiCtx, param_in_addr, "tmlim", &tmp_double, &tmp_res, 10000, Log, CHECK_NONE);
  lpx_set_real_parm(lp, LPX_K_TMLIM, tmp_double);

  getDoubleInPList(pvApiCtx, param_in_addr, "outdly", &tmp_double, &tmp_res, 0.0, Log, CHECK_NONE);
  if (tmp_res!=-1) lpx_set_real_parm(lp, LPX_K_OUTDLY, tmp_double);

  getDoubleInPList(pvApiCtx, param_in_addr, "tolint", &tmp_double, &tmp_res, 1e-5, Log, CHECK_NONE);
  if (tmp_res!=-1) lpx_set_real_parm(lp, LPX_K_TOLINT, tmp_double);

  getDoubleInPList(pvApiCtx, param_in_addr, "tolobj", &tmp_double, &tmp_res, 1e-7, Log, CHECK_NONE);
  if (tmp_res!=-1) lpx_set_real_parm(lp, LPX_K_TOLOBJ, tmp_double);

  int lpsolver = 1;
  getIntInPList(pvApiCtx, param_in_addr, "lpsolver", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);

  //
  //  Save option
  //

  int save_pb = 0;
  char *save_filename = NULL;

  getStringInPList(pvApiCtx, param_in_addr, "writemps", &tmp_char, &tmp_res, "test.mps", Log, CHECK_NONE);
  save_pb = (save_pb!=-1);
  save_filename = tmp_char;

  // Define the number of unknowns and their domains.
  glp_add_cols(lp, n_var);
  for (i=0; i<n_var; i++)
    {
      // Define type of the structural variables
      
      if (!std::isinf(*(lb+i)) && !std::isinf(*(ub+i))) glp_set_col_bnds(lp, i+1, GLP_DB, *(lb+i), *(ub+i));
      else
	{
	  if (!std::isinf(*(lb+i)) && std::isinf(*(ub+i))) glp_set_col_bnds(lp, i+1, GLP_LO, *(lb+i), *(ub+i));
	  else
	    {
	      if (std::isinf(*(lb+i)) && !std::isinf(*(ub+i))) glp_set_col_bnds(lp, i+1, GLP_UP, *(lb+i), *(ub+i));
	      else                                             glp_set_col_bnds(lp, i+1, GLP_FR, *(lb+i), *(ub+i));
	    }
	}
    
      // Set the objective coefficient of the corresponding
      // structural variable. No constant term is assumed.
      glp_set_obj_coef(lp,i+1,*(c+i));

      switch(*(vartype+i))
	{
	case 'i':
	case 'I':
	  glp_set_col_kind(lp, i+1, GLP_IV);
	  isMIP = 1;
	  break;
	case 'b':
	case 'B':
	  glp_set_col_kind(lp, i+1, GLP_BV);
	  isMIP = 1;
	  break;
	default:
	  glp_set_col_kind(lp, i+1, GLP_CV);
	  break;
	}

#ifdef DEBUG
      sciprint("i = %d lb = %f ub = %f c = %f vartype = %c\n", i, *(lb+i), *(ub+i), *(c+i), *(vartype+i));
#endif
    }

  glp_add_rows(lp, n_constr);

#ifdef DEBUG
  sciprint("DEBUG: glp->m = %d\n", glp_get_num_rows(lp));
  sciprint("DEBUG: glp->n = %d\n", glp_get_num_cols(lp));
#endif

  for (i=0; i<n_constr; i++)
    {
      // If the i-th row has no lower bound (types F,U), the correspondent parameter will be ignored.
      // If the i-th row has no upper bound (types F,L), the correspondent parameter will be ignored.
      // If the i-th row is of S type, the i-th LB is used, but the i-th UB is ignored.

      // 'L' - smaller than - <=
      // 'E' - equality     - =
      // 'G' - greater than - >=
      // 'R' - Range        - <= + >=
      // 'N' - Free         - no constraints

      switch (*(ctype+i))
	{
	case 'N': // free bound
	case 'n':
	  typx = GLP_FR;
	  break;
	case 'G': // upper bound
	case 'g':
	  typx = GLP_LO;
	  break;
	case 'L': // lower bound
	case 'l':
	  typx = GLP_UP;
	  break;
	case 'E': // fixed constraint
	case 'e':
	  typx = GLP_FX;
	  break;
	case 'R': // double-bounded variable
	case 'r':
	  typx = GLP_DB;
	  break;
	}
      
      glp_set_row_bnds(lp, i+1, typx, *(lhs+i), *(rhs+i));

#ifdef DEBUG
      sciprint("i = %d lhs = %f rhs = %f typex = %d\n", i, *(lhs+i), *(rhs+i), typx);
#endif
    }

  // Load constraint matrix A
  glp_load_matrix(lp, nz, glpk_rn, glpk_cn, glpk_a);

  // Save the problem in freemps format if required
  if (save_pb)
    {
      if (lpx_write_freemps(lp, save_filename) != 0)
	{
	  Scierror(999,"%s: unable to write the problem",fname);

	  freeAllocatedSingleString(ctype);
	  freeAllocatedSingleString(vartype);

	  if (save_filename) FREE(save_filename);
	  if (glpk_rn) delete [] glpk_rn;
	  if (glpk_cn) delete [] glpk_cn;
	  if (glpk_a)  delete [] glpk_a;

	  return 0;
	}
    }

  // scale the problem data (if required)
  // The parameter flags specifies scaling options used by the routine.
  // Options can be combined with the bitwise OR operator and may be the following:
  // GLP_SF_GM      perform geometric mean scaling;
  // GLP_SF_EQ      perform equilibration scaling;
  // GLP_SF_2N      round scale factors to nearest power of two;
  // GLP_SF_SKIP    skip scaling, if the problem is well scaled.
  // The parameter flags may be specified as GLP_SF_AUTO, in which case the routine chooses scaling options automatically.

  if (lpx_get_int_parm(lp, LPX_K_SCALE) && (!lpx_get_int_parm(lp, LPX_K_PRESOL) || lpsolver != 1))
    {
      getIntInPList(pvApiCtx, param_in_addr, "scale_flags", &tmp_int, &tmp_res, GLP_SF_AUTO, Log, CHECK_NONE);
      glp_scale_prob(lp, tmp_int);
    }

  // build advanced initial basis (if required)
  if (lpsolver == 1 && !lpx_get_int_parm(lp, LPX_K_PRESOL)) 
    {
      // basis type:
      // 0: standard
      // 1: advanced
      // 2: bixby basis
      getIntInPList(pvApiCtx, param_in_addr, "basis_type", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) 
	{
	  switch(tmp_int)
	    {
	    case 1:
	      // The parameter flags is reserved for use in the future and must be specified as zero.
	      glp_adv_basis(lp, 0);
	      break;
	    case 2:
	      glp_cpx_basis(lp);
	      break;
	    default:
	      glp_std_basis(lp);
	      break;
	    } // End Switch
	} // End If 
    } // End If 

  ////////////////////////////
  // set control parameters //
  ////////////////////////////

  glp_smcp smcp_param;
  glp_iocp iocp_param;

  /////////////////////
  // For the simplex //
  /////////////////////

  glp_init_smcp(&smcp_param);
  // remap of control parameters for simplex method
  smcp_param.msg_lev = lpx_get_int_parm(lp, LPX_K_MSGLEV); // message level
  
  // simplex method: primal/dual
  if (lpx_get_int_parm(lp, LPX_K_DUAL)==0) smcp_param.meth = GLP_PRIMAL;		
  else                                     smcp_param.meth = GLP_DUALP;
  
  // pricing technique 
  if (lpx_get_int_parm(lp, LPX_K_PRICE)==0) smcp_param.pricing = GLP_PT_STD;
  else                                      smcp_param.pricing = GLP_PT_PSE;
  
  // smcp_param.r_test not available 
  
  smcp_param.tol_bnd = lpx_get_real_parm(lp, LPX_K_TOLBND); // primal feasible tolerance
  smcp_param.tol_dj  = lpx_get_real_parm(lp, LPX_K_TOLDJ);  // dual feasible tolerance
  smcp_param.tol_piv = lpx_get_real_parm(lp, LPX_K_TOLPIV); // pivot tolerance
  smcp_param.obj_ll  = lpx_get_real_parm(lp, LPX_K_OBJLL);  // lower limit
  smcp_param.obj_ul  = lpx_get_real_parm(lp, LPX_K_OBJUL);  // upper limit
  
  // iteration limit
  if (lpx_get_int_parm(lp, LPX_K_ITLIM)==-1) smcp_param.it_lim = INT_MAX;
  else                                       smcp_param.it_lim = lpx_get_int_parm(lp, LPX_K_ITLIM);
  
  // time limit
  if (lpx_get_real_parm(lp, LPX_K_TMLIM)==-1) smcp_param.tm_lim = INT_MAX;
  else                                        smcp_param.tm_lim = (int)lpx_get_real_parm(lp, LPX_K_TMLIM);
  
  smcp_param.out_frq = lpx_get_int_parm(lp, LPX_K_OUTFRQ);  // output frequency
  smcp_param.out_dly = (int)lpx_get_real_parm(lp, LPX_K_OUTDLY); // output delay
  
  // presolver
  if (lpx_get_int_parm(lp, LPX_K_PRESOL)) smcp_param.presolve = GLP_ON;
  else                                    smcp_param.presolve = GLP_OFF;

  //////////////////////////////
  // For the branch and bound //
  //////////////////////////////

  glp_init_iocp(&iocp_param);

  // add the ctrl+c callback
  iocp_param.cb_func = cb_func;

  // remap of control parameters for simplex method
  iocp_param.msg_lev = lpx_get_int_parm(lp, LPX_K_MSGLEV); // message level
  
  // branching technique:
  // GLP_BR_FFV 1 - first fractional variable
  // GLP_BR_LFV 2 - last fractional variable
  // GLP_BR_MFV 3 - most fractional variable
  // GLP_BR_DTH 4 - heuristic by Driebeck and Tomlin
  
  getIntInPList(pvApiCtx, param_in_addr, "br_tech", &tmp_int, &tmp_res, GLP_BR_FFV, Log, CHECK_NONE);
  if (tmp_res!=-1) iocp_param.br_tech = tmp_int;
  
  // backtracking technique:
  // GLP_BT_DFS 1 - depth first search
  // GLP_BT_BFS 2 - breadth first search
  // GLP_BT_BLB 3 - best local bound
  // GLP_BT_BPH 4 - best projection heuristic
  
  getIntInPList(pvApiCtx, param_in_addr, "bt_tech", &tmp_int, &tmp_res, GLP_BT_DFS, Log, CHECK_NONE);
  if (tmp_res!=-1) iocp_param.bt_tech = tmp_int;
  
  getDoubleInPList(pvApiCtx, param_in_addr, "tolint", &tmp_double, &tmp_res, 1e-5, Log, CHECK_NONE);
  if (tmp_res!=-1) iocp_param.tol_int = tmp_double;
  
  getDoubleInPList(pvApiCtx, param_in_addr, "tolobj", &tmp_double, &tmp_res, 1e-7, Log, CHECK_NONE);
  if (tmp_res!=-1) iocp_param.tol_obj = tmp_double;
  
  getDoubleInPList(pvApiCtx, param_in_addr, "tmlim", &tmp_double, &tmp_res, -1, Log, CHECK_NONE);
  if (tmp_res!=-1) iocp_param.tm_lim = (int)tmp_double;

  // Already defined precedingly
  iocp_param.out_frq = lpx_get_int_parm(lp, LPX_K_OUTFRQ);
  
  getDoubleInPList(pvApiCtx, param_in_addr, "outdly", &tmp_double, &tmp_res, 0.0, Log, CHECK_NONE);
  if (tmp_res!=-1) iocp_param.out_dly = (int)tmp_double;
  
  // preprocessing technique:
  // GLP_PP_NONE 0 - disable preprocessing
  // GLP_PP_ROOT 1 - preprocessing only on root level
  // GLP_PP_ALL  2 - preprocessing on all levels
  getIntInPList(pvApiCtx, param_in_addr, "pp_tech", &tmp_int, &tmp_res, GLP_PP_ROOT, Log, CHECK_NONE);
  if (tmp_res!=-1) iocp_param.pp_tech = tmp_int;
  
  // relative MIP gap tolerance
  getDoubleInPList(pvApiCtx, param_in_addr, "mip_gap", &tmp_double, &tmp_res, 0.0, Log, CHECK_NONE);
  if (tmp_res!=-1) iocp_param.mip_gap = tmp_double;
  
  // MIR cuts (GLP_ON/GLP_OFF)
  getIntInPList(pvApiCtx, param_in_addr, "mir_cuts", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) iocp_param.mir_cuts = tmp_int;
  
  // GOMORY cuts (GLP_ON/GLP_OFF)
  getIntInPList(pvApiCtx, param_in_addr, "gmi_cuts", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) iocp_param.gmi_cuts = tmp_int;
	  
  // COVER cuts (GLP_ON/GLP_OFF)
  getIntInPList(pvApiCtx, param_in_addr, "cov_cuts", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) iocp_param.cov_cuts = tmp_int;
  
  // CLIQUE cuts (GLP_ON/GLP_OFF)
  getIntInPList(pvApiCtx, param_in_addr, "clq_cuts", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) iocp_param.clq_cuts = tmp_int;
  
  // enable/disable using MIP presolver
  getIntInPList(pvApiCtx, param_in_addr, "presol", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) iocp_param.presolve = tmp_int;
  // We enable the presol option because, when solving an integer programming problem, it is required to
  // first to solve a relaxed problem and the to perform the search.
  // Using the GLPK API, you can solve the problem using a simplex and then call the branch and bound part
  // We let glpk do the presolve thing for us
  // try to binarize integer variables 
  //iocp_param.presolve = 1;
  
  getIntInPList(pvApiCtx, param_in_addr, "binarize", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) iocp_param.binarize = tmp_int;

  // Sometimes we want to solve a relaxed problem
  getIntInPList(pvApiCtx, param_in_addr, "solve_relaxed", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) isRelaxed = tmp_int;

  if (isRelaxed) isMIP = 0;

  ///////////////////////
  // Solve the problem //
  ///////////////////////

  // Dual costs

  // Standard formulae are used:

  // Current values of simplex multipliers pi[i], i = 1, ..., m (which are values of Lagrange multipliers 
  // for equality constraints (7) also called shadow prices) are computed as follows:
  //
  // pi = inv(B') * cB,                                            (12)
  //
  // where B' is a matrix transposed to B, cB is a vector of objective coefficients at basic variables xB.
  //
  // Current values of reduced costs d[j], j = 1, ..., n, (which are values of Langrange multipliers for
  // active inequality constraints corresponding to non-basic variables) are computed as follows:
  //
  // d = cN - N' * pi,                                             (13)
  //
  // where N' is a matrix transposed to N, cN is a vector of objective coefficients at non-basic variables xN.
  //
  // -> The problem that this is easiest to examine this on is the "pk1.mps" from the MIPLib2003 problem set.
  // -> The very first node is returning 1 5zero-valued dual costs for the first node.
  // I see nothing unusual. Zero reduced costs just say that the problem is dual degenerative. 

  double * xmin = NULL, * redcosts = NULL, * lambda = NULL;
  double fmin[1], status[1], errnum[1], time[1], mem[1];

  xmin     = (double *)MALLOC(n_var*sizeof(double));
  redcosts = (double *)MALLOC(n_var*sizeof(double));
  lambda   = (double *)MALLOC(n_constr*sizeof(double));

  // give a default value to xmin, redcosts, lambda and fmin

  for(i=0; i<n_var; i++)
    {
      // Primal values
      *(xmin+i) = 0.0;
      // Reduced costs
      *(redcosts+i) = 0.0;
    }
  // Dual values
  for (i=0; i<n_constr; i++) *(lambda+i) = 0.0;
  fmin[0]   = 0.0;
  status[0] = 0;

  // Solve
  switch (lpsolver)
    {
    case 1: 
      {
	/////////////////////////////
	// Deal with a MIP problem //
	/////////////////////////////

	if (isMIP)
	  {
	    // glp_intopt - solve MIP problem with the branch-and-bound method
	    //
	    // DESCRIPTION
	    //
	    // The routine glp_intopt is a driver to the MIP solver based on the branch-and-bound method.
	    //
	    // On entry the problem object should contain optimal solution to LP relaxation (which can be
	    // obtained with the routine glp_simplex).
	    //
	    // The MIP solver has a set of control parameters. Values of the control
	    // parameters can be passed in a structure glp_iocp, which the parameter parm points to.
	    //
	    // The parameter parm can be specified as NULL, in which case the MIP solver uses default settings.
	    //
	    // RETURNS
	    //
	    // 0                    The MIP problem instance has been successfully solved. This code does not necessarily mean that the solver has found optimal
	    //                      solution. It only means that the solution process was successful.
	    // GLP_EBOUND  - 0x04 - Unable to start the search, because some double-bounded variables
	    //                      have incorrect bounds or some integer variables have non-integer (fractional) bounds.
	    // GLP_EROOT   - 0x0C - Unable to start the search, because optimal basis for initial LP relaxation is not provided.
	    // GLP_EFAIL   - 0x05 - The search was prematurely terminated due to the solver failure.
	    // GLP_EMIPGAP - 0x0E - The search was prematurely terminated, because the relative mip gap tolerance has been reached.
	    // GLP_ETMLIM  - 0x09 - The search was prematurely terminated, because the time limit has been exceeded.
	    // GLP_ENOPFS  - 0x0A - The MIP problem instance has no primal feasible solution (only if the MIP presolver is used).
	    // GLP_ENODFS  - 0x0B - LP relaxation of the MIP problem instance has no dual feasible solution (only if the MIP presolver is used).
	    // GLP_ESTOP   - 0x0D - The search was prematurely terminated by application.

#ifdef DEBUG
	    sciprint("MIP solver\n");
#endif
	    getIntInPList(pvApiCtx, param_in_addr, "mip_presolve", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
	    switch(tmp_int)
	      {
	      case 1: // simplex
		errnum[0] = (double)glp_simplex(lp, &smcp_param);
		iocp_param.presolve = 0;
		// If a problem occured, we switch to the default presolve method of glp_intopt
		if (errnum[0]) iocp_param.presolve = 1;
		break;
	      case 2: // interior point
		errnum[0] = (double)glp_interior(lp,NULL); 
		iocp_param.presolve = 0;
		// If a problem occured, we switch to the default presolve method of glp_intopt
		if (errnum[0]) iocp_param.presolve = 1;
		break;
	      case 3: // extact resolution
		errnum[0] = (double)glp_exact(lp,&smcp_param); 
		iocp_param.presolve = 0;
		// If a problem occured, we switch to the default presolve method of glp_intopt
		if (errnum[0]) iocp_param.presolve = 1;
		break;
	      default: // use the glp_intopt presolve
		iocp_param.presolve = 1;
	      }
	    // Now perform integer programming
	    errnum[0] = (double)glp_intopt(lp, &iocp_param);
	    status[0] = (double)glp_mip_status(lp);
#ifdef DEBUG
	    sciprint("status = %f\n", status);
#endif
 	    if ((status[0]==(double)GLP_OPT) || (status[0]==(double)GLP_FEAS))
 	      {
		fmin[0] = glp_mip_obj_val(lp);
		
		for(i=0; i<n_var; i++)
		  {
		    *(xmin+i)       = glp_mip_col_val(lp, i+1);
		    *(redcosts + i) = 0.0;
		  }
		for(i=0; i<n_constr; i++) *(lambda+i) = glp_mip_row_val(lp, i+1);
 	      }
	  }
	else
	  {
	    ////////////////////////////
	    // Deal with a LP problem //
	    ////////////////////////////

	    // glp_simplex - solve LP problem with the simplex method
	    //
	    // DESCRIPTION
	    //
	    // The routine glp_simplex is a driver to the LP solver based on the simplex method.
	    // This routine retrieves problem data from the specified problem object, calls the solver to solve the problem
	    // instance, and stores results of computations back into the problem object.
	    //
	    // The simplex solver has a set of control parameters. Values of the control parameters can be passed in a
	    // structure glp_smcp, which the parameter parm points to.
	    //
	    // The parameter parm can be specified as NULL, in which case the LP solver uses default settings.
	    //
	    // RETURNS
	    //
	    // 0                    The LP problem instance has been successfully solved. This code does not necessarily mean that
	    //                      the solver has found optimal solution. It only means that the solution process was successful.
	    // GLP_EBADB   - 0x01 - Unable to start the search, because the initial basis specified
	    //                      in the problem object is invalid--the number of basic (auxiliary and structural) variables is
	    //                      not the same as the number of rows in the problem object.
	    // GLP_ESING   - 0x02 - Unable to start the search, because the basis matrix correspodning
	    //                      to the initial basis is singular within the working precision.
	    // GLP_ECOND   - 0x03 - Unable to start the search, because the basis matrix correspodning
	    //                      to the initial basis is ill-conditioned, i.e. its condition number is too large.
	    // GLP_EBOUND  - 0x04 - Unable to start the search, because some double-bounded variables have incorrect bounds.
	    // GLP_EFAIL   - 0x05 - The search was prematurely terminated due to the solver failure.
	    // GLP_EOBJLL  - 0x06 - The search was prematurely terminated, because the objective function being maximized has
	    //                      reached its lower limit and continues decreasing (dual simplex only).
	    // GLP_EOBJUL  - 0x07 - The search was prematurely terminated, because the objective function being minimized has
	    //                      reached its upper limit and continues increasing (dual simplex only).
	    // GLP_EITLIM  - 0x08 - The search was prematurely terminated, because the simplex iteration limit has been exceeded.
	    // GLP_ETMLIM  - 0x09 - The search was prematurely terminated, because the time limit has been exceeded.
	    // GLP_ENOPFS  - 0x0A - The LP problem instance has no primal feasible solution (only if the LP presolver is used).
	    // GLP_ENODFS  - 0x0B - The LP problem instance has no dual feasible solution (only if the LP presolver is used).

#ifdef DEBUG
	    sciprint("SIMPLEX chosen\n");
#endif
	    errnum[0] = (double)glp_simplex(lp, &smcp_param);
	    status[0] = (double)glp_get_status(lp);
#ifdef DEBUG
	    sciprint("status = %f\n", status[0]);
#endif
	    
 	    if ((status[0]==(double)GLP_OPT) ||
		(status[0]==(double)GLP_FEAS) ||
		(status[0]==(double)GLP_INFEAS) ||
		(status[0]==(double)GLP_NOFEAS) ||
		(status[0]==(double)GLP_UNBND) ||
		(status[0]==(double)GLP_UNDEF))
 	      {
		fmin[0] = glp_get_obj_val(lp);
		
		for (i=0; i<n_var; i++)
		  {
		    // Primal values
		    *(xmin+i) = glp_get_col_prim(lp, i+1);
		    // Reduced costs
		    *(redcosts+i) = glp_get_col_dual(lp, i+1);
		  }
		// Dual values
		for (i=0; i<n_constr; i++)
		  {
		    *(lambda+i) = glp_get_row_dual(lp, i+1);
		  }
 	      }
	  }
      }
      break;

    case 2: 
      ////////////////////////////
      // Deal with a LP problem //
      ////////////////////////////
      
      // glp_exact - solve LP problem in exact arithmetic
      //
      // DESCRIPTION
      //
      // The routine glp_exact is a tentative implementation of the primal two-phase simplex method based on exact (rational) arithmetic. It is
      // similar to the routine glp_simplex, however, for all internal computations it uses arithmetic of rational numbers, which is exact
      // in mathematical sense, i.e. free of round-off errors unlike floating point arithmetic.
      //
      // Note that the routine glp_exact uses inly two control parameters passed in the structure glp_smcp, namely, it_lim and tm_lim.
      //
      // RETURNS
      //
      // 0          The LP problem instance has been successfully solved. This code does not necessarily mean that the solver has found optimal
      //            solution. It only means that the solution process was successful.
      // GLP_EBADB  Unable to start the search, because the initial basis specified in the problem object is invalid--the number of basic (auxiliary
      //            and structural) variables is not the same as the number of rows in the problem object.
      // GLP_ESING  Unable to start the search, because the basis matrix correspodning to the initial basis is exactly singular.
      // GLP_EBOUND Unable to start the search, because some double-bounded variables have incorrect bounds.
      // GLP_EFAIL  The problem has no rows/columns.
      // GLP_EITLIM The search was prematurely terminated, because the simplex iteration limit has been exceeded.
      // GLP_ETMLIM The search was prematurely terminated, because the time limit has been exceeded.

#ifdef DEBUG
      sciprint("exec solver chosen\n");
#endif

      errnum[0] = (double)glp_exact(lp, &smcp_param); 
      status[0] = (double)glp_get_status(lp);

#ifdef DEBUG
      sciprint("status = %f\n", status[0]);
#endif

      if ((status[0]==(double)GLP_OPT) ||
	  (status[0]==(double)GLP_FEAS) ||
	  (status[0]==(double)GLP_INFEAS) ||
	  (status[0]==(double)GLP_NOFEAS) ||
	  (status[0]==(double)GLP_UNBND) ||
	  (status[0]==(double)GLP_UNDEF))
 	{
	  fmin[0] = glp_get_obj_val(lp);
	  
	  for (i=0; i<n_var; i++)
	    {
	      // Primal values
	      *(xmin+i) = glp_get_col_prim(lp, i+1);
	      // Reduced costs
	      *(redcosts+i) = glp_get_col_dual(lp, i+1);
	    }
	  // Dual values
	  for (i=0; i<n_constr; i++)
	    {
	      *(lambda+i) = glp_get_row_dual(lp, i+1);
	    }
 	}
      break;

    default: 
      ////////////////////////////
      // Deal with a LP problem //
      ////////////////////////////
      
      // glp_interior - solve LP problem with the interior-point method
      //
      // DESCRIPTION
      //
      // The routine glp_interior is a driver to the LP solver based on the interior-point method.
      // The parameter parm is reserved for use in the future and must be specified as NULL.
      // Currently this routine implements an easy variant of the primal-dual interior-point method based on Mehrotra's technique.
      //
      // This routine transforms the original LP problem to an equivalent LP problem in the standard formulation (all constraints are equalities,
      // all variables are non-negative), calls the routine ipm_main to solve the transformed problem, and then transforms an obtained solution to
      // the solution of the original problem.
      //
      // RETURNS
      //
      // 0           The LP problem instance has been successfully solved. This code does not necessarily mean that the solver has found optimal
      //             solution. It only means that the solution process was successful.
      // GLP_EFAIL   The problem has no rows/columns.
      // GLP_ENOFEAS The problem has no feasible (primal/dual) solution.
      // GLP_ENOCVG  Very slow convergence or divergence.
      // GLP_EITLIM  Iteration limit exceeded.
      // GLP_EINSTAB Numerical instability on solving Newtonian system.

#ifdef DEBUG
      sciprint("interior point chosen\n");
#endif

      errnum[0] = (double)glp_interior(lp,NULL); 
      status[0] = (double)glp_ipt_status(lp);

#ifdef DEBUG
      sciprint("status = %f\n", status);
#endif

      if (status[0]==(double)GLP_OPT)
 	{
	  fmin[0] = glp_ipt_obj_val(lp);
	  
	  for (i=0; i<n_var; i++)
	    {
	      // Primal values
	      *(xmin+i) = glp_ipt_col_prim(lp, i+1);
	      // Reduced costs
	      *(redcosts+i) = glp_ipt_col_dual(lp, i+1);
	    }
	  // Dual values
	  for (i=0; i<n_constr; i++)
	    {
	      *(lambda+i) = glp_ipt_row_dual(lp, i+1);
	    }
 	}
    }

//   if (status!=GLP_OPT)
//     {
//       // A problem has been raised, the value of the solution is not defined
//       memset(xmin,     0, sizeof(double)*n_var);
//       memset(redcosts, 0, sizeof(double)*n_var);
//       memset(lambda,   0, sizeof(double)*n_constr);
//       fmin = 0.0;
//     }

  time[0] = (clock() - t_start) / CLOCKS_PER_SEC;

#ifdef DEBUG
  sciprint("DEBUG: status = %f\n", status);
  sciprint("DEBUG: time   = %f\n", time);
#endif

  size_t tpeak;
  glp_mem_usage(NULL, NULL, NULL, &tpeak);
  // Memory measured in ko ?
  mem[0] = tpeak / (1024);

  if (errnum[0])
    {
      if (msglev)
	{
	  sciprint("sciglpk: an error occured: errnum = %d\n", (int)errnum[0]);
	}
    }

  _SciErr = createMatrixOfDouble(pvApiCtx, XMIN_OUT,   n_var, 1, xmin); SCICOINOR_ERROR;
  _SciErr = createMatrixOfDouble(pvApiCtx, FMIN_OUT,   1,     1, fmin); SCICOINOR_ERROR;
  _SciErr = createMatrixOfDouble(pvApiCtx, STATUS_OUT, 1,     1, status); SCICOINOR_ERROR;

  char * ListLabels [] = {"lambda","redcosts","time","mem","errnum"};

  _SciErr = createPList(pvApiCtx, EXTRA_OUT, &extra_addr, (char **)ListLabels, 5); SCICOINOR_ERROR;

  _SciErr = createColVectorOfDoubleInPList(pvApiCtx, EXTRA_OUT, extra_addr, "lambda", n_constr, lambda); SCICOINOR_ERROR;
  _SciErr = createColVectorOfDoubleInPList(pvApiCtx, EXTRA_OUT, extra_addr, "redcosts", n_var, redcosts); SCICOINOR_ERROR;
  _SciErr = createDoubleInPList(pvApiCtx,            EXTRA_OUT, extra_addr, "time",   time[0]); SCICOINOR_ERROR;
  _SciErr = createDoubleInPList(pvApiCtx,            EXTRA_OUT, extra_addr, "mem",    mem[0]); SCICOINOR_ERROR;
  _SciErr = createIntInPList(pvApiCtx,               EXTRA_OUT, extra_addr, "errnum", errnum[0]); SCICOINOR_ERROR;

  LhsVar(1) = XMIN_OUT;
  LhsVar(2) = FMIN_OUT;
  LhsVar(3) = STATUS_OUT;
  LhsVar(4) = EXTRA_OUT;

  if (glpk_rn) delete [] glpk_rn;
  if (glpk_cn) delete [] glpk_cn;
  if (glpk_a)  delete [] glpk_a;

  glp_delete_prob(lp);

  freeAllocatedSingleString(ctype);
  freeAllocatedSingleString(vartype);
  if (save_filename) FREE(save_filename);

  if (xmin)     FREE(xmin);
  if (redcosts) FREE(redcosts);
  if (lambda)   FREE(lambda);

  return 0;
}
