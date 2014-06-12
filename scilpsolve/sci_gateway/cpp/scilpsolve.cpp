//////////////////////////////////////////////////////////////////////////////////
// scilpsolve: A scilab interface to the LPSOLVE library for linear programming //
//////////////////////////////////////////////////////////////////////////////////

//  Copyright (C) 2008-2010 Yann Collette.
//  Copyright (C) 2001-2007 Nicolo' Giorgetti.
//
//  SCILPSOLVE is free software; you can redistribute it and/or modify it
//  under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 2, or (at your option) any
//  later version.
//
//  This part of code is distributed with the FURTHER condition that it 
//  can be compiled and linked with the Matlab libraries and it can be 
//  used within the Matlab environment.
//
//  SCILPSOLVE is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with SciCoinOR; see the file COPYING.  If not, write to the Free
//  Software Foundation, 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#ifdef WIN32
#include <windows.h>
#endif

#include <string.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>

extern "C"
{
#include <stack-c.h>
#include <Scierror.h>
#include <sciprint.h>
#include <MALLOC.h>
#include <api_parameters.h>
}

#include <api_scilab.h>

#include <lp_lib.h>
#include <lp_simplex.h>

#define DEBUG 1

#include <helper.hpp>

int __WINAPI abortfunction(lprec *lp, void *userhandle)
{
  int doabort = (C2F(basbrk).iflag == -1);

  if (doabort)
    {
      C2F(basbrk).iflag = 0;
    }
#ifdef DEBUG
  sciprint("abortfunction: iflag = %d interruptible = %d\n", C2F(basbrk).iflag, C2F(basbrk).interruptible);
#endif

  return doabort;
}

void __WINAPI msgfunction(lprec *lp, void *userhandle, int msg)
{
  switch(msg)
    {
    case MSG_PRESOLVE:
      sciprint("scilpsolve: Presolve done\n");
      break;
    case MSG_LPFEASIBLE:
      sciprint("scilpsolve: Feasible solution found\n");
      break;
    case MSG_LPOPTIMAL:
      sciprint("scilpsolve: Real optimal solution found. Only fired when there are integer variables at the start of B&B\n");
      break;
    case MSG_MILPEQUAL:
      sciprint("scilpsolve: Equal MILP solution found. Only fired when there are integer variables during B&B\n");
      break;
    case MSG_MILPFEASIBLE:
      sciprint("scilpsolve: Integer feasible solution found\n");
      break;
    case MSG_MILPBETTER:
      sciprint("scilpsolve: Better integer feasible solution found\n");
      break;
    }
}

void __WINAPI logfunction(lprec *lp, void *userhandle, char *buf)
{
  sciprint("scilpsolve: %s", buf);
}
    
// Input arguments
#define	C_IN	               1
#define	A_IN	               2
#define	LHS_IN	               3
#define	RHS_IN	               4
#define LB_IN	               5
#define UB_IN                  6
#define CTYPE_IN               7
#define VARTYPE_IN             8
#define PARAM_IN               9
#define SPECIAL_CONSTR_IN      10
#define LAST_PARAM             10
// Output Arguments
#define	 XMIN_OUT              Rhs+1
#define	 FMIN_OUT              Rhs+2
#define	 STATUS_OUT            Rhs+3
#define  EXTRA_OUT             Rhs+4

#define SPECIAL_WHICH_IN       2
#define SPECIAL_WEIGHT_IN      3
#define SPECIAL_TYPE_IN        4
#define SPECIAL_COLUMN_IN      5
#define SPECIAL_RANGE_IN       6
#define SPECIAL_ID_OBJ_IN      7
#define SPECIAL_LENGTH_IN      8
#define SPECIAL_ID_IN          9
#define SPECIAL_CLIQUE_TYPE_IN 10

extern "C" int scilpsolve(char * fname)
{
  int m_c,   n_c,   * c_addr = NULL;
  int m_a,   n_a,   * a_addr = NULL;
  int m_lhs, n_lhs, * lhs_addr = NULL;
  int m_rhs, n_rhs, * rhs_addr = NULL;
  int m_lb,  n_lb,  * lb_addr = NULL;
  int m_ub,  n_ub,  * ub_addr = NULL;
  int * ctype_addr = NULL;
  int * vartype_addr = NULL;
  double * c = NULL, * a = NULL, * lhs = NULL;
  double * rhs = NULL, * lb = NULL, * ub = NULL;
  char * ctype = NULL, * vartype = NULL;
  int    * rn = NULL, * cn = NULL;
  double * an = NULL;
  lprec  * lp = NULL;
  char   * tmp_char = NULL;
  double  tmp_double;
  int     tmp_int, tmp_int_1, tmp_int_2, tmp_res, tmp_res_1, tmp_res_2, status = 0, nz = 0;
  int i, j, n_var, n_constr, type, Log = 0;
  int majorversion, minorversion, release, build;
  SciErr _SciErr;
  
  lp_solve_version(&majorversion, &minorversion, &release, &build);

  if(!(majorversion>=5) && !(minorversion>=5) && !(release>=0) && !(build>=12))
    {
      Scierror(999,"%s: this interface is compatible only with LPSOLVE version 5.5.0.12 or higher.\n",fname);
      return 0;
    }
  
  if (Rhs != LAST_PARAM)
    {
      sciprint("Interface to LPSOLVE Version %d.%d.%d.%d\n",majorversion,minorversion,release,build);
      sciprint("Internal interface for the LPSOLVE library.\n");
      sciprint("You should use the 'lpsolve' function instead.\n\n");
      sciprint("Syntax: [xopt, fmin, status, extra] = scilpsolve(c, a, lhs, rhs, lb, ub, ctype, vartype, param, special)\n");
      return 0;
    }
 
  //////////////////////////////////////////////////////////////////////////////
  // 1nd Input. A column array containing the objective function coefficients //
  //////////////////////////////////////////////////////////////////////////////

  _SciErr = getVarAddressFromPosition(pvApiCtx, C_IN, &c_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, c_addr, &n_c, &m_c, &c); SCICOINOR_ERROR;

  n_var = n_c*m_c;

  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // 2nd Input. A matrix containing the constraints coefficients. If matrix A is NOT a sparse matrix //
  /////////////////////////////////////////////////////////////////////////////////////////////////////

  _SciErr = getVarAddressFromPosition(pvApiCtx, A_IN, &a_addr); SCICOINOR_ERROR;
  _SciErr = getVarType(pvApiCtx, a_addr, &type); SCICOINOR_ERROR;
  if(type!=sci_sparse)
    {
#ifdef DEBUG
      sciprint("DEBUG: A_IN is not sparse\n");
#endif
      
      _SciErr = getMatrixOfDouble(pvApiCtx, a_addr, &m_a, &n_a, &a); SCICOINOR_ERROR;
      
      rn = (int *)MALLOC((m_a*n_a+1)*sizeof(int));
      cn = (int *)MALLOC((m_a*n_a+1)*sizeof(int));
      an = (double *)MALLOC((m_a*n_a+1)*sizeof(double));
      
      nz = 0;

      for(i=0; i<m_a; i++)
	{
	  for(j=0; j<n_a; j++)
	    {
	      if (*(a+i*n_a+j) != 0)
		{
		  nz++;
 		  rn[nz-1] = i+1;
 		  cn[nz-1] = j+1;
		  an[nz-1] = *(a+i+j*m_a);
		}
	    }
	}
    }
  else
    {
#ifdef DEBUG
      sciprint("DEBUG: A_IN is sparse\n");
#endif
      
      SciSparse S_A;
      
      getAllocatedSparseMatrix(pvApiCtx, a_addr, &S_A.m, &S_A.n, &S_A.nel, &S_A.mnel, &S_A.icol, &S_A.R);
      
      nz = S_A.nel;
      
      rn = (int *)MALLOC((nz+1)*sizeof(int));
      cn = (int *)MALLOC((nz+1)*sizeof(int));
      an = (double *)MALLOC((nz+1)*sizeof(double));
      
      int count = -1;
      for(i=0;i<S_A.m;i++)
	{
	  if (S_A.mnel[i]!=0) 
	    {
	      for(j=0;j<S_A.mnel[i];j++)
		{
		  count++;
		  rn[count] = i+1;
		  cn[count] = S_A.icol[count];
		  an[count] = S_A.R[count];
		}
	    }
	}
      nz = count + 1;
      freeAllocatedSparseMatrix(S_A.mnel, S_A.icol, S_A.R);
    }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // 3rd Input. A column array containing the left-hand side value for each constraint in the constraint matrix //
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  _SciErr = getVarAddressFromPosition(pvApiCtx, LHS_IN, &lhs_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, lhs_addr, &n_lhs, &m_lhs, &lhs); SCICOINOR_ERROR;

  n_constr = m_lhs*n_lhs;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // 4th Input. A column array containing the right-hand side value for each constraint in the constraint matrix //
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  _SciErr = getVarAddressFromPosition(pvApiCtx, RHS_IN, &rhs_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, rhs_addr, &n_rhs, &m_rhs, &rhs); SCICOINOR_ERROR;

  n_constr = n_constr<(m_rhs*n_rhs)?(m_rhs*n_rhs):n_constr;

  /////////////////////////////////////////////////////////////////////////////////////////////
  // 5th Input. An array of length n_var containing the lower bound on each of the variables //
  /////////////////////////////////////////////////////////////////////////////////////////////

  _SciErr = getVarAddressFromPosition(pvApiCtx, LB_IN, &lb_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, lb_addr, &n_lb, &m_lb, &lb); SCICOINOR_ERROR;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // 6th Input. An array of at least length numcols containing the upper bound on each of the variables //
  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  _SciErr = getVarAddressFromPosition(pvApiCtx, UB_IN, &ub_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, ub_addr, &n_ub, &m_ub, &ub); SCICOINOR_ERROR;

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // 7th Input. A column array containing the sense of each constraint in the constraint matrix //
  ////////////////////////////////////////////////////////////////////////////////////////////////

  _SciErr = getVarAddressFromPosition(pvApiCtx, CTYPE_IN, &ctype_addr); SCICOINOR_ERROR;
  getAllocatedSingleString(pvApiCtx, ctype_addr, &ctype);

  /////////////////////////////////////////////////////////////////////
  // 8th Input. A column array containing the types of the variables //
  /////////////////////////////////////////////////////////////////////

  _SciErr = getVarAddressFromPosition(pvApiCtx, VARTYPE_IN, &vartype_addr); SCICOINOR_ERROR;
  getAllocatedSingleString(pvApiCtx, vartype_addr, &vartype);

#ifdef DEBUG
  sciprint("DEBUG m_vartype = %d n_vartype = %d vartype = %s\n", m_vartype, n_vartype, vartype);
  sciprint("DEBUG m_lhs = %d n_lhs = %d\n", m_lhs, n_lhs);
  for(i=0;i<m_lhs*n_lhs;i++) sciprint("lhs[%d] = %f\n",i,lhs[i]);
  sciprint("DEBUG m_c = %d n_c = %d\n", m_c, n_c);
  for(i=0;i<m_c*n_c;i++) sciprint("c[%d] = %f\n",i,c[i]); 
  sciprint("DEBUG A - nz = %d\n", nz);
  for(i=0;i<nz;i++) sciprint("a[%d] = %f - cn[%d] = %d - rn[%d] = %d\n",i,an[i],i,rn[i],i,cn[i]);
  sciprint("DEBUG m_rhs = %d n_rhs = %d\n", m_rhs, n_rhs);
  for(i=0;i<m_rhs*n_rhs;i++) sciprint("rhs[%d] = %f\n",i,rhs[i]);
  sciprint("DEBUG m_lb = %d n_lb = %d n_var = %d\n", m_lb, n_lb, n_var);
  for(i=0;i<n_var;i++) sciprint("lb[%d] = %f\n",i,lb[i]); 
  sciprint("DEBUG m_ub = %d n_ub = %d n_var = %d \n", m_ub, n_ub, n_var);
  for(i=0;i<n_var;i++) sciprint("ub[%d] = %f\n",i,ub[i]);
  sciprint("DEBUG m_ctype = %d n_ctype = %d ctype = %s\n", m_ctype, n_ctype, ctype);
#endif

  ////////////////////////////////////////
  // Initialization of the lp structure //
  ////////////////////////////////////////

  lp = make_lp(n_constr, n_var);
  if (lp==NULL)
    {
      sciprint("scilpsolve: unable to initialize the lp structure\n");

      freeAllocatedSingleString(ctype);
      freeAllocatedSingleString(vartype);

      return 0;
    }

  //////////////////////////////////////////////////////////////
  // 9th Input. A structure containing the control parameters //
  //////////////////////////////////////////////////////////////
  //
  // Lpsolve parameters

  // 'anti_degen'      -> int    Specifies if special handling must be done to reduce degeneracy/cycling while solving.    
  // 'verbose'         -> int    message level
  // 'pivoting'        -> int    Sets the pivot rule and mode. PRICER_* and PRICE_* can be ORed
  // 'epsb'            -> double the value that is used as a tolerance for the Right Hand Side (RHS) to determine whether
  //                             a value should be considered as 0
  // 'epsd'            -> double the value that is used as a tolerance for the reduced costs to determine whether a value 
  //                             should be considered as 0.
  // 'epspivot'        -> double the value that is used as a tolerance for the pivot element to determine whether a value
  //                             should be considered as 0.
  // 'epsel'           -> double the value that is used as a tolerance for rounding values to zero.
  // 'epsint'          -> double the tolerance that is used to determine whether a floating-point number is in fact an integer.
  // 'epsperturb'      -> double the value that is used as perturbation scalar for degenerative problems.
  // 'infinite'        -> double Specifies the practical value for "infinite".
  // 'break_at_first'  -> int    Specifies if the branch-and-bound algorithm stops at first found solution.
  // 'break_at_value'  -> double Specifies if the branch-and-bound algorithm stops when the object value is better than a given value.
  // 'basiscrash'      -> int     
  // 'bb_depthlimit'   -> int    Sets the maximum branch-and-bound depth.
  // 'bb_floorfirst'   -> int    Specifies which branch to take first in branch-and-bound algorithm.
  // 'debug'           -> int    Sets a flag if all intermediate results and the branch-and-bound decisions must be printed while solving.
  // 'lag_trace'       -> int    Sets a flag if Lagrangian progression must be printed while solving.
  // 'maxpivot'        -> int    Sets the maximum number of pivots between a re-inversion of the matrix.
  // 'mip_gap'         -> int    Absolute gap if 1
  //                   -> double Specifies the MIP gap value. 17 -> TRUE / FALSE, 18 -> gap
  // 'bb_rule'         -> int     
  // 'preferdual'      -> int    Sets the desired combination of primal and dual simplex algorithms.
  // 'presolve'        -> int    Do presolve in 1
  //                   -> int    maxloops - Specifies if a presolve must be done before solving.
  // 'scalelimit'      -> double Sets the relative scaling convergence criterion for the active scaling mode; 
  //                             the integer part specifies the maximum number of iterations.
  // 'scaling'         -> int    Specifies which scaling algorithm must be used.
  // 'solutionlimit'   -> int    Sets the solution number that must be returned
  //                             (for problem with binary, integer or semicontinuous variables.
  // 'timeout'         -> int    set a timeout in second.
  // 'trace'           -> int    Sets a flag if pivot selection must be printed while solving.
  // 'dualize'         -> int    
  // 'method'          -> int    lagragian relaxation or classical
  //                   -> double Starting bound
  //                   -> int    Number of iterations
  
  int * param_in_addr = NULL;

  initPList(pvApiCtx, PARAM_IN, &param_in_addr);
  if (!checkPList(pvApiCtx, param_in_addr))
    {
      Scierror(999, "%s: argument n° %d is not a plist\n", fname, PARAM_IN);

      freeAllocatedSingleString(ctype);
      freeAllocatedSingleString(vartype);

      return 0;
    }
  
  //
  // First, we process the solver option to allocate the good solver
  //

  ////////////////////////
  // Integer parameters //
  ////////////////////////

  // 'anti_degen' - int - Specifies if special handling must be done to reduce degeneracy/cycling while solving.    
  // Strategy codes to avoid or recover from degenerate pivots, infeasibility or numeric errors via randomized bound relaxation 
  // ANTIDEGEN_NONE           0
  // ANTIDEGEN_FIXEDVARS      1
  // ANTIDEGEN_COLUMNCHECK    2
  // ANTIDEGEN_STALLING       4
  // ANTIDEGEN_NUMFAILURE     8
  // ANTIDEGEN_LOSTFEAS      16
  // ANTIDEGEN_INFEASIBLE    32
  // ANTIDEGEN_DYNAMIC       64
  // ANTIDEGEN_DURINGBB     128
  // ANTIDEGEN_RHSPERTURB   256
  // ANTIDEGEN_BOUNDFLIP    512
  // ANTIDEGEN_DEFAULT      (ANTIDEGEN_FIXEDVARS | ANTIDEGEN_STALLING | ANTIDEGEN_INFEASIBLE)
  getIntInPList(pvApiCtx, param_in_addr, "anti_degen", &tmp_int, &tmp_res, ANTIDEGEN_DEFAULT, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) set_anti_degen(lp, tmp_int);

  // 'verbose' - int - message level
  //
  // MSG_NONE                 0
  // MSG_PRESOLVE             1
  // MSG_ITERATION            2
  // MSG_INVERT               4
  // MSG_LPFEASIBLE           8
  // MSG_LPOPTIMAL           16
  // MSG_LPEQUAL             32
  // MSG_LPBETTER            64
  // MSG_MILPFEASIBLE       128
  // MSG_MILPEQUAL          256
  // MSG_MILPBETTER         512
  // MSG_MILPSTRATEGY      1024
  // MSG_MILPOPTIMAL       2048
  // MSG_PERFORMANCE       4096
  // MSG_INITPSEUDOCOST    8192
  getIntInPList(pvApiCtx, param_in_addr, "verbose", &tmp_int, &tmp_res, MSG_NONE, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) set_verbose(lp, tmp_int);

  // 'pivoting' -> int -> Sets the pivot rule and mode. PRICER_* and PRICE_* can be ORed
  // 
  // Pricing methods
  //   
  // PRICER_FIRSTINDEX        0
  // PRICER_DANTZIG           1
  // PRICER_DEVEX             2
  // PRICER_STEEPESTEDGE      3
  // PRICER_LASTOPTION        PRICER_STEEPESTEDGE
  //
  // Pricing strategies 
  //   
  // PRICE_PRIMALFALLBACK     4    // In case of Steepest Edge, fall back to DEVEX in primal
  // PRICE_MULTIPLE           8    // Enable multiple pricing (primal simplex)
  // PRICE_PARTIAL           16    // Enable partial pricing
  // PRICE_ADAPTIVE          32    // Temporarily use alternative strategy if cycling is detected 
  // PRICE_RANDOMIZE        128    // Adds a small randomization effect to the selected pricer 
  // PRICE_AUTOPARTIAL      256    // Detect and use data on the block structure of the model (primal) 
  // PRICE_AUTOMULTIPLE     512    // Automatically select multiple pricing (primal simplex) 
  // PRICE_LOOPLEFT        1024    // Scan entering/leaving columns left rather than right 
  // PRICE_LOOPALTERNATE   2048    // Scan entering/leaving columns alternatingly left/right 
  // PRICE_HARRISTWOPASS   4096    // Use Harris' primal pivot logic rather than the default 
  // PRICE_FORCEFULL       8192    // Non-user option to force full pricing 
  // PRICE_TRUENORMINIT   16384    // Use true norms for Devex and Steepest Edge initializations 
  getIntInPList(pvApiCtx, param_in_addr, "pivoting", &tmp_int, &tmp_res, PRICER_DEVEX | PRICE_ADAPTIVE, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) set_pivoting(lp, tmp_int);

  // 'epsb' -> double -> the value that is used as a tolerance for the Right Hand Side (RHS) to determine 
  //                     whether a value should be considered as 0
  getDoubleInPList(pvApiCtx, param_in_addr, "epsb", &tmp_double, &tmp_res, 1e-10, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) set_epsb(lp, tmp_double);

  // 'epsd' -> double -> the value that is used as a tolerance for the reduced costs to determine 
  //                     whether a value should be considered as 0.
  getDoubleInPList(pvApiCtx, param_in_addr, "epsd", &tmp_double, &tmp_res, 1e-9, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) set_epsd(lp, tmp_double);

  // 'epspivot' -> double -> the value that is used as a tolerance for the pivot element to determine
  //                         whether a value should be considered as 0.
  getDoubleInPList(pvApiCtx, param_in_addr, "epspivot", &tmp_double, &tmp_res, 2e-7, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) set_epspivot(lp, tmp_double);

  // 'epsel' -> double -> the value that is used as a tolerance for rounding values to zero.
  getDoubleInPList(pvApiCtx, param_in_addr, "epsel", &tmp_double, &tmp_res, 1e-12, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) set_epsel(lp, tmp_double);

  // 'epsint' -> double -> the tolerance that is used to determine whether a floating-point number is in fact an integer.
  getDoubleInPList(pvApiCtx, param_in_addr, "epsint", &tmp_double, &tmp_res, 1e-7, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) set_epsint(lp, tmp_double);

  // 'epsperturb' -> double -> the value that is used as perturbation scalar for degenerative problems.
  getDoubleInPList(pvApiCtx, param_in_addr, "epsperturb", &tmp_double, &tmp_res, 1e-5, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) set_epsperturb(lp, tmp_double);

  // 'epslevel' -> int -> This is a simplified way of specifying multiple eps thresholds that are "logically" consistent.
  // EPS_TIGHT  (0) Very tight epsilon values (default)
  // EPS_MEDIUM (1) Medium epsilon values
  // EPS_LOOSE  (2) Loose epsilon values
  // EPS_BAGGY  (3) Very loose epsilon values
  getIntInPList(pvApiCtx, param_in_addr, "epslevel", &tmp_int, &tmp_res, EPS_TIGHT, Log, CHECK_VALUES, 4, 0, 1, 2, 3);
  if (tmp_res!=-1) set_epslevel(lp, tmp_int);

  // 'infinite' -> double -> Specifies the practical value for "infinite".
  getDoubleInPList(pvApiCtx, param_in_addr, "infinite", &tmp_double, &tmp_res, 1e30, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) set_infinite(lp, tmp_double);

  // 'break_at_first' -> int -> Specifies if the branch-and-bound algorithm stops at first found solution.
  getIntInPList(pvApiCtx, param_in_addr, "break_at_first", &tmp_int, &tmp_res, FALSE, Log, CHECK_NONE);
  if (tmp_res!=-1) set_break_at_first(lp, tmp_int);

  // 'break_at_value' -> boolean -> Specifies if the branch-and-bound algorithm stops when the object value is better than a given value.
  getDoubleInPList(pvApiCtx, param_in_addr, "break_at_value", &tmp_double, &tmp_res, -1e30, Log, CHECK_NONE);
  if (tmp_res!=-1) set_break_at_value(lp, tmp_double);

  // 'basiscrash' -> int -> ?
  //
  // Basis crash options
  //
  // CRASH_NONE               0
  // CRASH_NONBASICBOUNDS     1
  // CRASH_MOSTFEASIBLE       2
  // CRASH_LEASTDEGENERATE    3
  getIntInPList(pvApiCtx, param_in_addr, "basiscrash", &tmp_int, &tmp_res, CRASH_NONE, Log, CHECK_VALUES, 4, 0, 1, 2, 3);
  if (tmp_res!=-1) set_basiscrash(lp, tmp_int);

  // 'bb_depthlimit' -> int -> Sets the maximum branch-and-bound depth.
  getIntInPList(pvApiCtx, param_in_addr, "bb_depthlimit", &tmp_int, &tmp_res, -50, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) set_bb_depthlimit(lp, tmp_int);

  // 'bb_floorfirst' -> int -> Specifies which branch to take first in branch-and-bound algorithm.
  //
  // BRANCH_CEILING           0
  // BRANCH_FLOOR             1
  // BRANCH_AUTOMATIC         2
  // BRANCH_DEFAULT           3
  getIntInPList(pvApiCtx, param_in_addr, "bb_floorfirst", &tmp_int, &tmp_res, BRANCH_AUTOMATIC, Log, CHECK_VALUES, 4, 0, 1, 2, 3);
  if (tmp_int!=-1) set_bb_floorfirst(lp, tmp_int);

  // 'bb_rule' -> int -> Specifies the branch-and-bound rule.
  //
  // NODE_FIRSTSELECT         0
  // NODE_GAPSELECT           1
  // NODE_RANGESELECT         2
  // NODE_FRACTIONSELECT      3
  // NODE_PSEUDOCOSTSELECT    4
  // NODE_PSEUDONONINTSELECT  5    // Kjell Eikland #1 - Minimize B&B depth
  // NODE_PSEUDOFEASSELECT   (NODE_PSEUDONONINTSELECT+NODE_WEIGHTREVERSEMODE)
  // NODE_PSEUDORATIOSELECT   6    // Kjell Eikland #2 - Minimize a "cost/benefit" ratio 
  // NODE_USERSELECT          7
  // NODE_STRATEGYMASK        (NODE_WEIGHTREVERSEMODE-1) // Mask for B&B strategies
  // NODE_WEIGHTREVERSEMODE   8
  // NODE_BRANCHREVERSEMODE  16
  // NODE_GREEDYMODE         32
  // NODE_PSEUDOCOSTMODE     64
  // NODE_DEPTHFIRSTMODE    128
  // NODE_RANDOMIZEMODE     256
  // NODE_GUBMODE           512
  // NODE_DYNAMICMODE      1024
  // NODE_RESTARTMODE      2048
  // NODE_BREADTHFIRSTMODE 4096
  // NODE_AUTOORDER        8192
  // NODE_RCOSTFIXING     16384
  // NODE_STRONGINIT      32768
  getIntInPList(pvApiCtx, param_in_addr, "bb_rule", &tmp_int, &tmp_res, NODE_PSEUDONONINTSELECT + NODE_GREEDYMODE + NODE_DYNAMICMODE + NODE_RCOSTFIXING, Log, CHECK_MIN, 0);
  if (tmp_int!=-1) set_bb_rule(lp, tmp_int);

  // 'debug' -> int -> Sets a flag if all intermediate results and the branch-and-bound decisions must be printed while solving.
  getIntInPList(pvApiCtx, param_in_addr, "debug", &tmp_int, &tmp_res, FALSE, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) set_debug(lp, tmp_res);

  // 'lag_trace' -> int -> Sets a flag if Lagrangian progression must be printed while solving.
  getIntInPList(pvApiCtx, param_in_addr, "lag_trace", &tmp_int, &tmp_res, FALSE, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) set_lag_trace(lp, tmp_int);

  // 'maxpivot' -> int -> Sets the maximum number of pivots between a re-inversion of the matrix.
  getIntInPList(pvApiCtx, param_in_addr, "maxpivot", &tmp_int, &tmp_res, 250, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) set_maxpivot(lp, tmp_int);

  // 'mip_gap_abs_[12]' -> int -> Absolute gap if 1
  getIntInPList(pvApiCtx, param_in_addr, "mip_gap_abs", &tmp_int, &tmp_res_1, 0, Log, CHECK_MIN, 0);
  getDoubleInPList(pvApiCtx, param_in_addr, "mip_gap_gap", &tmp_double, &tmp_res_2, 1e-11, Log, CHECK_MIN, 0);
  if ((tmp_res_1!=-1)&&(tmp_res_2!=-1)) set_mip_gap(lp, tmp_int, tmp_double);

  // 'improve' -> int -> Specifies the iterative improvement level.
  // IMPROVE_NONE      (0) improve none
  // IMPROVE_SOLUTION  (1) Running accuracy measurement of solved equations based on Bx=r (primal simplex), remedy is refactorization.
  // IMPROVE_DUALFEAS  (2) Improve initial dual feasibility by bound flips (highly recommended, and default)
  // IMPROVE_THETAGAP  (4) Low-cost accuracy monitoring in the dual, remedy is refactorization
  // IMPROVE_BBSIMPLEX (8) By default there is a check for primal/dual feasibility at optimum only for the relaxed problem, 
  //                       this also activates the test at the node level
  getIntInPList(pvApiCtx, param_in_addr, "improve", &tmp_int, &tmp_res, IMPROVE_DUALFEAS + IMPROVE_THETAGAP, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) set_improve(lp, tmp_int);

  // 'bounds_tighter' -> boolean -> Specifies if set bounds may only be tighter or also less restrictive.
  getIntInPList(pvApiCtx, param_in_addr, "bounds_tighter", &tmp_int, &tmp_res, 0, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) set_bounds_tighter(lp, tmp_int);

  // To Add:
  // - set_BFP - need a string parameter
  // - set_XLI

  // 'preferdual' -> int -> Sets the desired combination of primal and dual simplex algorithms.
  //
  // SIMPLEX_UNDEFINED        0
  // SIMPLEX_Phase1_PRIMAL    1
  // SIMPLEX_Phase1_DUAL      2
  // SIMPLEX_Phase2_PRIMAL    4
  // SIMPLEX_Phase2_DUAL      8
  // SIMPLEX_DYNAMIC         16
  // SIMPLEX_AUTODUALIZE     32
  //
  // SIMPLEX_PRIMAL_PRIMAL   (SIMPLEX_Phase1_PRIMAL + SIMPLEX_Phase2_PRIMAL)
  // SIMPLEX_DUAL_PRIMAL     (SIMPLEX_Phase1_DUAL   + SIMPLEX_Phase2_PRIMAL)
  // SIMPLEX_PRIMAL_DUAL     (SIMPLEX_Phase1_PRIMAL + SIMPLEX_Phase2_DUAL)
  // SIMPLEX_DUAL_DUAL       (SIMPLEX_Phase1_DUAL   + SIMPLEX_Phase2_DUAL)
  // SIMPLEX_DEFAULT         (SIMPLEX_DUAL_PRIMAL)
  getIntInPList(pvApiCtx, param_in_addr, "preferdual", &tmp_int, &tmp_res, SIMPLEX_DEFAULT, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) set_preferdual(lp, tmp_int);

  // 'presolve' -> int -> Do presolve in 1
  //
  // PRESOLVE_NONE            0
  // PRESOLVE_ROWS            1
  // PRESOLVE_COLS            2
  // PRESOLVE_LINDEP          4
  // PRESOLVE_SOS            32
  // PRESOLVE_REDUCEMIP      64
  // PRESOLVE_KNAPSACK      128  // Implementation not tested completely
  // PRESOLVE_ELIMEQ2       256
  // PRESOLVE_IMPLIEDFREE   512
  // PRESOLVE_REDUCEGCD    1024
  // PRESOLVE_PROBEFIX     2048
  // PRESOLVE_PROBEREDUCE  4096
  // PRESOLVE_ROWDOMINATE  8192
  // PRESOLVE_COLDOMINATE 16384  // Reduced functionality, should be expanded 
  // PRESOLVE_MERGEROWS   32768
  // PRESOLVE_IMPLIEDSLK  65536
  // PRESOLVE_COLFIXDUAL 131072
  // PRESOLVE_BOUNDS     262144
  // PRESOLVE_LASTMASKMODE    (PRESOLVE_DUALS - 1)
  // PRESOLVE_DUALS      524288
  // PRESOLVE_SENSDUALS 1048576
  //
  // Can be ORed
  getIntInPList(pvApiCtx, param_in_addr, "presolve_do", &tmp_int_1, &tmp_res_1, PRESOLVE_NONE, Log, CHECK_MIN, 0);
  // 'presolve_maxloops' -> int -> maxloops - Specifies if a presolve must be done before solving.
  getIntInPList(pvApiCtx, param_in_addr, "presolve_maxloops", &tmp_int_2, &tmp_res_2, 200, Log, CHECK_MIN, 0);
  if ((tmp_res_1!=-1)&&(tmp_res_2==-1)) set_presolve(lp, tmp_int_1, get_presolveloops(lp));
  if ((tmp_res_1!=-1)&&(tmp_res_2!=-1)) set_presolve(lp, tmp_int_1, tmp_int_2);


  // 'scalelimit'      -> double -> Sets the relative scaling convergence criterion for the active scaling mode; 
  //                                the integer part specifies the maximum number of iterations.
  getDoubleInPList(pvApiCtx, param_in_addr, "scalelimit", &tmp_double, &tmp_res, 5, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) set_scalelimit(lp, tmp_double);

  // 'scaling' -> int -> Specifies which scaling algorithm must be used.
  //
  // SCALE_NONE               0
  // SCALE_EXTREME            1
  // SCALE_RANGE              2
  // SCALE_MEAN               3
  // SCALE_GEOMETRIC          4
  // SCALE_FUTURE1            5
  // SCALE_FUTURE2            6
  // SCALE_CURTISREID         7   // Override to Curtis-Reid "optimal" scaling
  //  
  // Alternative scaling weights 
  //  
  // SCALE_LINEAR             0
  // SCALE_QUADRATIC          8
  // SCALE_LOGARITHMIC       16
  // SCALE_USERWEIGHT        31
  // SCALE_MAXTYPE            (SCALE_QUADRATIC-1)
  //  
  // Scaling modes 
  //  
  // SCALE_POWER2            32   // As is or rounded to power of 2 
  // SCALE_EQUILIBRATE       64   // Make sure that no scaled number is above 1 
  // SCALE_INTEGERS         128   // Apply to integer columns/variables 
  // SCALE_DYNUPDATE        256   // Apply incrementally every solve() 
  // SCALE_ROWSONLY         512   // Override any scaling to only scale the rows 
  // SCALE_COLSONLY        1024   // Override any scaling to only scale the rows 
  //
  // Can be ORed
  getIntInPList(pvApiCtx, param_in_addr, "scaling", &tmp_int, &tmp_res, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) set_scaling(lp, tmp_int);

  // 'solutionlimit' -> int -> Sets the solution number that must be returned (for problem with binary, integer or semicontinuous variables.
  getIntInPList(pvApiCtx, param_in_addr, "solutionlimit", &tmp_int, &tmp_res, 1, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) set_solutionlimit(lp, tmp_int);

  // 'timeout' -> int -> set a timeout in second.
  getIntInPList(pvApiCtx, param_in_addr, "timeout", &tmp_int, &tmp_res, 10000, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) set_timeout(lp, tmp_int);
  
  // 'trace' -> int -> Sets a flag if pivot selection must be printed while solving.
  getIntInPList(pvApiCtx, param_in_addr, "trace", &tmp_int, &tmp_res, FALSE, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) set_trace(lp, tmp_int);

  // 'negrange' -> double -> Set negative value below which variables are split into a negative and a positive part.
  getDoubleInPList(pvApiCtx, param_in_addr, "nerange", &tmp_double, &tmp_res, -1e-6, Log, CHECK_NONE);
  if (tmp_res!=-1) set_negrange(lp, tmp_double);

  // 'sense' -> int -> optimization direction
  getIntInPList(pvApiCtx, param_in_addr, "sense", &tmp_int, &tmp_res, 0, Log, CHECK_MIN, 0);
  if (tmp_res!=-1)
    {
      if (tmp_int==1) set_maxim(lp);
      else            set_minim(lp);
    }
  
  // 'simplextype' -> int -> Sets the desired combination of primal and dual simplex algorithms.
  // SIMPLEX_PRIMAL_PRIMAL (5) Phase1 Primal, Phase2 Primal
  // SIMPLEX_DUAL_PRIMAL (6)   Phase1 Dual, Phase2 Primal (default value)
  // SIMPLEX_PRIMAL_DUAL (9)   Phase1 Primal, Phase2 Dual
  // SIMPLEX_DUAL_DUAL (10)    Phase1 Dual, Phase2 Dual  
  getIntInPList(pvApiCtx, param_in_addr, "simplextype", &tmp_int, &tmp_res, SIMPLEX_DUAL_PRIMAL, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) set_simplextype(lp, tmp_int);

  /////////////////////////////
  // Process SOS constraints //
  /////////////////////////////

  int nb_constr_which        = 0, * special_constr_addr      = NULL;
  int nb_special_which       = 0, * special_which_addr       = NULL;
  int nb_special_weight      = 0, * special_weight_addr      = NULL;
  int nb_special_type        = 0, * special_type_addr        = NULL;
  int nb_special_column      = 0, * special_column_addr      = NULL;
  int nb_special_range       = 0, * special_range_addr       = NULL;
  int nb_special_id_obj      = 0, * special_id_obj_addr      = NULL;
  int nb_special_length      = 0, * special_length_addr      = NULL;
  int nb_special_id          = 0, * special_id_addr          = NULL;
  int nb_special_clique_type = 0, * special_clique_type_addr = NULL;
  int m_elem_which          = 0, n_elem_which          = 0;
  int m_elem_weight         = 0, n_elem_weight         = 0;
  int m_elem_type           = 0, n_elem_type           = 0;
  int m_elem_column         = 0, n_elem_column         = 0;
  int m_elem_range          = 0, n_elem_range          = 0;
  int m_elem_id_obj         = 0, n_elem_id_obj         = 0;
  int m_elem_length         = 0, n_elem_length         = 0;
  int m_elem_id             = 0, n_elem_id             = 0;
  int m_elem_clique_type    = 0, n_elem_clique_type    = 0;
  double * elem_which = NULL, * elem_weight = NULL, * elem_type = NULL;
  double * elem_column = NULL, * elem_range = NULL, * elem_id_obj = NULL;
  double * elem_length = NULL, * elem_id = NULL, * elem_clique_type = NULL;
  bool special_stored = false;
  int nb_elem = 0, var_type;

  // We must get something from the stack for these parameters so as to be able to CreateVar for the outputs

  _SciErr = getVarAddressFromPosition(pvApiCtx, SPECIAL_CONSTR_IN, &special_constr_addr); SCICOINOR_ERROR;

  if (!isEmptyMatrix(pvApiCtx, special_constr_addr))
    {
      _SciErr = getVarType(pvApiCtx, special_constr_addr, &var_type); SCICOINOR_ERROR;
      if (var_type!=sci_mlist) 
	{
	  Scierror(999,"%s: parameter %d must be a mlist\n", fname, SPECIAL_CONSTR_IN);
	  return 0;
	}
      
      _SciErr = getListInList(pvApiCtx, special_constr_addr, SPECIAL_WHICH_IN, &special_which_addr); SCICOINOR_ERROR;
      _SciErr = getVarType(pvApiCtx, special_which_addr, &var_type); SCICOINOR_ERROR;
      if (var_type!=sci_list) 
	{
	  Scierror(999,"%s: special element %d must be a list\n", fname, SPECIAL_WHICH_IN);
	  return 0;
	}
      
      _SciErr = getListInList(pvApiCtx, special_constr_addr, SPECIAL_WEIGHT_IN, &special_weight_addr); SCICOINOR_ERROR;
      _SciErr = getVarType(pvApiCtx, special_weight_addr, &var_type); SCICOINOR_ERROR;
      if (var_type!=sci_list) 
	{
	  Scierror(999,"%s: special element %d must be a list\n", fname, SPECIAL_WEIGHT_IN);
	  return 0;
	}
      
      _SciErr = getListInList(pvApiCtx, special_constr_addr, SPECIAL_TYPE_IN, &special_type_addr); SCICOINOR_ERROR;
      _SciErr = getVarType(pvApiCtx, special_type_addr, &var_type); SCICOINOR_ERROR;
      if (var_type!=sci_list) 
	{
	  Scierror(999,"%s: special element %d must be a list\n", fname, SPECIAL_TYPE_IN);
	  return 0;
	}
      
      _SciErr = getListInList(pvApiCtx, special_constr_addr, SPECIAL_COLUMN_IN, &special_column_addr); SCICOINOR_ERROR;
      _SciErr = getVarType(pvApiCtx, special_column_addr, &var_type); SCICOINOR_ERROR;
      if (var_type!=sci_list) 
	{
	  Scierror(999,"%s: special element %d must be a list\n", fname, SPECIAL_COLUMN_IN);
	  return 0;
	}
      
      _SciErr = getListInList(pvApiCtx, special_constr_addr, SPECIAL_RANGE_IN, &special_range_addr); SCICOINOR_ERROR;
      _SciErr = getVarType(pvApiCtx, special_range_addr, &var_type); SCICOINOR_ERROR;
      if (var_type!=sci_list) 
	{
	  Scierror(999,"%s: special element %d must be a list\n", fname, SPECIAL_RANGE_IN);
	  return 0;
	}
      
      _SciErr = getListInList(pvApiCtx, special_constr_addr, SPECIAL_ID_OBJ_IN, &special_id_obj_addr); SCICOINOR_ERROR;
      _SciErr = getVarType(pvApiCtx, special_id_obj_addr, &var_type); SCICOINOR_ERROR;
      if (var_type!=sci_list) 
	{
	  Scierror(999,"%s: special element %d must be a list\n", fname, SPECIAL_ID_OBJ_IN);
	  return 0;
	}
      
      _SciErr = getListInList(pvApiCtx, special_constr_addr, SPECIAL_LENGTH_IN, &special_length_addr); SCICOINOR_ERROR;
      _SciErr = getVarType(pvApiCtx, special_length_addr, &var_type); SCICOINOR_ERROR;
      if (var_type!=sci_list) 
	{
	  Scierror(999,"%s: special element %d must be a list\n", fname, SPECIAL_LENGTH_IN);
	  return 0;
	}
      
      _SciErr = getListInList(pvApiCtx, special_constr_addr, SPECIAL_ID_IN, &special_id_addr); SCICOINOR_ERROR;
      _SciErr = getVarType(pvApiCtx, special_id_addr, &var_type); SCICOINOR_ERROR;
      if (var_type!=sci_list) 
	{
	  Scierror(999,"%s: special element %d must be a list\n", fname, SPECIAL_ID_IN);
	  return 0;
	}
      
      _SciErr = getListInList(pvApiCtx, special_constr_addr, SPECIAL_CLIQUE_TYPE_IN, &special_clique_type_addr); SCICOINOR_ERROR;
      _SciErr = getVarType(pvApiCtx, special_clique_type_addr, &var_type); SCICOINOR_ERROR;
      if (var_type!=sci_list) 
	{
	  Scierror(999,"%s: special element %d must be a list\n", fname, SPECIAL_CLIQUE_TYPE_IN);
	  return 0;
	}
      
      _SciErr = getListItemNumber(pvApiCtx, special_which_addr,       &nb_special_which);
      _SciErr = getListItemNumber(pvApiCtx, special_weight_addr,      &nb_special_weight);
      _SciErr = getListItemNumber(pvApiCtx, special_type_addr,        &nb_special_type);
      _SciErr = getListItemNumber(pvApiCtx, special_column_addr,      &nb_special_column);
      _SciErr = getListItemNumber(pvApiCtx, special_range_addr,       &nb_special_range);
      _SciErr = getListItemNumber(pvApiCtx, special_id_obj_addr,      &nb_special_id_obj);
      _SciErr = getListItemNumber(pvApiCtx, special_length_addr,      &nb_special_length);
      _SciErr = getListItemNumber(pvApiCtx, special_id_addr,          &nb_special_id);
      _SciErr = getListItemNumber(pvApiCtx, special_clique_type_addr, &nb_special_clique_type);
      
#ifdef DEBUG
      sciprint("nb_special_which  = %d\n", nb_special_which);
      sciprint("nb_special_weight = %d\n", nb_special_weight);
      sciprint("nb_special_type   = %d\n", nb_special_type);
      sciprint("nb_special_column = %d\n", nb_special_column);
      sciprint("nb_special_range  = %d\n", nb_special_range);
      sciprint("nb_special_id_obj = %d\n", nb_special_id_obj);
      sciprint("nb_special_length = %d\n", nb_special_length);
      sciprint("nb_special_id     = %d\n", nb_special_id);
      sciprint("nb_special_clique_type = %d\n", nb_special_clique_type);
#endif
      
      if (nb_special_which != 0) special_stored = true;
    }
  else
    {
      special_stored = false;
    }

  if (special_stored)
    {
      nb_elem = nb_special_id;

      for(i=0;i<nb_elem; i++)
	{
	  // Get the element i from the parameters list
	  _SciErr = getMatrixOfDoubleInList(pvApiCtx, special_which_addr,  i+1, &m_elem_which,  &n_elem_which,  &elem_which);
	  _SciErr = getMatrixOfDoubleInList(pvApiCtx, special_weight_addr, i+1, &m_elem_weight, &n_elem_weight, &elem_weight);
	  _SciErr = getMatrixOfDoubleInList(pvApiCtx, special_type_addr,   i+1, &m_elem_type,   &n_elem_type,   &elem_type);
	  _SciErr = getMatrixOfDoubleInList(pvApiCtx, special_column_addr, i+1, &m_elem_column, &n_elem_column, &elem_column);
	  _SciErr = getMatrixOfDoubleInList(pvApiCtx, special_range_addr,  i+1, &m_elem_range,  &n_elem_range,  &elem_range);
	  _SciErr = getMatrixOfDoubleInList(pvApiCtx, special_id_obj_addr, i+1, &m_elem_id_obj, &n_elem_id_obj, &elem_id_obj);
	  _SciErr = getMatrixOfDoubleInList(pvApiCtx, special_length_addr, i+1, &m_elem_length, &n_elem_length, &elem_length);
	  _SciErr = getMatrixOfDoubleInList(pvApiCtx, special_id_addr,     i+1, &m_elem_id,     &n_elem_id,     &elem_id);
	  _SciErr = getMatrixOfDoubleInList(pvApiCtx, special_clique_type_addr, i+1, &m_elem_clique_type, &n_elem_clique_type, &elem_clique_type);
	  
	  if ((int)*elem_id_obj==0) // SOS constraint
	    {
	      int * int_elem_which = (int *)MALLOC(m_elem_which*n_elem_which*sizeof(int));
	      for(i=0;i<m_elem_which*n_elem_which;i++) int_elem_which[i] = elem_which[i];
	      add_SOS(lp, "SOS", (int)*elem_type, 1, (int)*elem_length, int_elem_which, elem_weight);
	      if (int_elem_which) FREE(int_elem_which);
	    }
	}
    }

  ///////////
  // TO DO //
  ///////////

  // YC: option to add ??
  // default_basis
  // get_basis
  // get_basiscrash
  // guess_basis
  // read_basis
  // reset_basis
  // set_basis
  // set_basiscrash
  // set_basisvar
  // write_basis

  // get_sensitivity_obj, get_ptr_sensitivity_obj, get_sensitivity_objex, get_ptr_sensitivity_objex
  // get_sensitivity_rhs, get_ptr_sensitivity_rhs, get_dual_solution, get_ptr_dual_solution, get_var_dualresult

  // is_feasible
  // add_lag_con, str_add_lag_con
  // add_SOS, is_SOS_var

  // TO DO: add a branch (int) parameter to lpsolve
  // TO DO: add a Parameters (int / double) parameter to lpsolve
  // TO DO: add a Weights (double) parameter to lpsolve

  //
  // BRANCH_CEILING           0
  // BRANCH_FLOOR             1
  // BRANCH_AUTOMATIC         2
  // BRANCH_DEFAULT           3
  //
  
  put_abortfunc(lp, abortfunction, NULL);
  put_msgfunc(lp, msgfunction, NULL, MSG_PRESOLVE | MSG_LPFEASIBLE | MSG_LPOPTIMAL | MSG_MILPEQUAL | MSG_MILPFEASIBLE | MSG_MILPBETTER);
  put_logfunc(lp, logfunction, NULL);

  // Define the constraints
  for (i=0; i<n_constr; i++)
    {
      switch(ctype[i])
	{
	case 'N': // FR (0) Free
	case 'n': // FR (0) Free
	  tmp_int = set_constr_type(lp, i+1, FR);
	  break;
	case 'G': // GE (2) Greater than or equal (>=)
	case 'g': // GE (2) Greater than or equal (>=)
	  tmp_int = set_constr_type(lp, i+1, GE);
	  set_rh(lp, i+1, lhs[i]);
	  break;
	case 'E': // EQ (3) Equal (=)
	case 'e': // EQ (3) Equal (=)
	  tmp_int = set_constr_type(lp, i+1, EQ);
	  set_rh(lp, i+1, rhs[i]);
	  break;
	case 'R': // Range
	case 'r': // Range
	  tmp_int = set_constr_type(lp, i+1, LE);
	  set_rh(lp, i+1, rhs[i]);
	  set_rh_range(lp, i+1, lhs[i] - rhs[i]);
	  break;
	default: // LE (1) Less than or equal (<=)
	  tmp_int = set_constr_type(lp, i+1, LE);
	  set_rh(lp, i+1, rhs[i]);
	  break;
	} 

#ifdef DEBUG
      sciprint("set_constr_type: constr %d - tmp_int = %d\n", i, tmp_int);
      sciprint("constr %d : rh = %f\n", i+1, get_rh(lp,i+1));
#endif
    }

#ifdef DEBUG
  printf("DEBUG: lp->rows   = %d / n_constr = %d\n", lp->rows,    n_constr);
  printf("DEBUG: lp->column = %d / n_var    = %d\n", lp->columns, n_var);
#endif

  for (i=0; i<nz; i++)
    {
      // Load constraint matrix A
      // row:    Row number of the matrix. Must be between 0 and number of rows in the lp. Row 0 is objective function.
      // column: Column number of the matrix. Must be between 1 and number of columns in the lp.
      set_mat(lp, rn[i], cn[i], an[i]);
    }

  set_constr_type(lp, 0, FR); // The objective function correspond to a FREE constraint

  for (i=0; i<n_var; i++)
    {
      set_obj(lp, i+1, c[i]);

      switch(vartype[i])
	{
	case 'B':
	case 'b':
	  set_binary(lp, i+1, 1);
	  break;
	case 'I':
	case 'i':
	  set_int(lp, i+1, 1);
	  break;
	case 'S':
	case 's':
	  set_semicont(lp, i+1, 1);
	  break;
	default:
	  set_binary(lp, i+1, 0);
	  set_int(lp, i+1, 0);
	  set_semicont(lp, i+1, 0);
	  break;
	}
 
      set_lowbo(lp, i+1, lb[i]);
      set_upbo(lp, i+1, ub[i]);
#ifdef DEBUG
      sciprint("var %d : lowbo = %f - upbo = %f\n", i+1, lb[i], ub[i]);
#endif
    }
  
#ifdef DEBUG
  sciprint("non zeros elements in the matrix: %d\n", get_nonzeros(lp));
#endif

  //
  // Dualisation of the problem
  //

  // 'dualize' -> int -> Create the dual of the current model.
  getIntInPList(pvApiCtx, param_in_addr, "dualize", &tmp_int, &tmp_res, FALSE, Log, CHECK_NONE);
  if (tmp_res!=-1) tmp_int = dualize_lp(lp);

  /////////////////
  // Save option //
  /////////////////

  getStringInPList(pvApiCtx, param_in_addr, "writemps", &tmp_char, &tmp_res, "test.mps", Log, CHECK_NONE);
  // YC: we have a problem with the writemps option.
  if (tmp_res!=-1) 
    {
      if (write_freemps(lp, tmp_char)==0)
	{
	  sciprint("scilpsolve: unable to write the problem in file %s\n",tmp_char);
	}
    }
  FREE(tmp_char);

  ///////////
  // Solve //
  ///////////
  // YC: Removed: lag_solve: we need to find an efficient way to deal with these kind of constraints
  // via a GET_PARAM_DOUBLE_VECTOR ?
  status = solve(lp);

  ///////////////////////////////////
  // Prepare the output parameters //
  ///////////////////////////////////

  double * xmin   = (double *)MALLOC(n_var*sizeof(double));
  double * lambda = (double *)MALLOC(n_constr*sizeof(double));
  double fmin = 0;
  int solutioncount, total_iter, total_nodes;
  double time, mem;

  // Solver status
  //
  // UNKNOWNERROR            -5
  // DATAIGNORED             -4
  // NOBFP                   -3
  // NOMEMORY                -2
  // NOTRUN                  -1
  // OPTIMAL                  0
  // SUBOPTIMAL               1
  // INFEASIBLE               2
  // UNBOUNDED                3
  // DEGENERATE               4
  // NUMFAILURE               5
  // USERABORT                6
  // TIMEOUT                  7
  // RUNNING                  8
  // PRESOLVED                9

  double * tmp_result_lpsolve = NULL;

  // if ((status==OPTIMAL)    || 
  //     (status==SUBOPTIMAL) ||
  //     (status==UNBOUNDED)  || 
  //     (status==DEGENERATE) || 
  //     (status==INFEASIBLE) ||
  //     (status==USERABORT)  ||
  //     (status==NUMFAILURE))
  if (status>=0)
    {
      // Primal values
      // An array that will contain the value of the objective function (element 0), values of the constraints (elements 1 till Nrows), 
      // and the values of the variables (elements Nrows+1 till Nrows+NColumns).
      get_ptr_primal_solution(lp, &tmp_result_lpsolve);
      for(i=0;i<n_var;i++)
	{
	  xmin[i] = get_var_primalresult(lp,n_constr + 1 + i);
	}

      // Dual values
      get_ptr_dual_solution(lp, &tmp_result_lpsolve);
      for(i=0;i<n_constr;i++)
	{
	  lambda[i] = get_var_dualresult(lp, i);
	}

      fmin = get_objective(lp);
#ifdef DEBUG
      sciprint("solution found, writing\n");
#endif
    }
  else
    {
      // Primal values
      memset(xmin, 0, n_var*sizeof(double));
      // Dual values
      memset(lambda, 0, n_constr*sizeof(double));

      fmin = 0.0;
#ifdef DEBUG
      sciprint("no solution found, initializing to 0\n");
#endif
    }
  
  //status        = status;
  solutioncount = get_solutioncount(lp);
  total_iter    = get_total_iter(lp);
  total_nodes   = get_total_nodes(lp);
  time          = time_elapsed(lp);
  mem           = lp->sum_alloc;

#ifdef DEBUG
  for(i=0;i<n_constr;i++) sciprint("lambda[%d] = %f\n",i,lambda[i]);
  for(i=0;i<n_var;i++)    sciprint("xmin[%d]   = %f\n",i,xmin[i]);
  sciprint("time          = %f\n", time);
  sciprint("mem           = %f\n", mem);
  sciprint("fmin          = %f\n", fmin);
  sciprint("status        = %d\n", status);
  sciprint("solutioncount = %d\n", solutioncount);
  sciprint("total_iter    = %d\n", total_iter);
  sciprint("total_nodes   = %d\n", total_nodes);
#endif
  
  // 
  // Store the informations related to the convergence in a mlist
  //

  int * extra_addr = NULL;
  double tmp_dbl[1];

  char * ListLabels [] = {"lambda","time","mem","solutioncount","total_iter","total_nodes"};

  _SciErr = createMatrixOfDouble(pvApiCtx, XMIN_OUT, n_var, 1, xmin); SCICOINOR_ERROR;
  createScalarDouble(pvApiCtx, FMIN_OUT, fmin);
  createScalarDouble(pvApiCtx, STATUS_OUT, (double)status);

  _SciErr = createPList(pvApiCtx, EXTRA_OUT, &extra_addr, (char **)ListLabels, 6); SCICOINOR_ERROR;

  _SciErr = createColVectorOfDoubleInPList(pvApiCtx, EXTRA_OUT, extra_addr, "lambda", n_constr, lambda); SCICOINOR_ERROR;
  _SciErr = createDoubleInPList(pvApiCtx,            EXTRA_OUT, extra_addr, "time",          time); SCICOINOR_ERROR;
  _SciErr = createDoubleInPList(pvApiCtx,            EXTRA_OUT, extra_addr, "mem",           mem); SCICOINOR_ERROR;
  _SciErr = createIntInPList(pvApiCtx,               EXTRA_OUT, extra_addr, "solutioncount", solutioncount); SCICOINOR_ERROR;
  _SciErr = createIntInPList(pvApiCtx,               EXTRA_OUT, extra_addr, "total_iter",    total_iter); SCICOINOR_ERROR;
  _SciErr = createIntInPList(pvApiCtx,               EXTRA_OUT, extra_addr, "total_nodes",   total_nodes); SCICOINOR_ERROR;

  LhsVar(1) = XMIN_OUT;
  LhsVar(2) = FMIN_OUT;
  LhsVar(3) = STATUS_OUT;
  LhsVar(4) = EXTRA_OUT;

  if (rn) FREE(rn);
  if (cn) FREE(cn);
  if (an) FREE(an);

  freeAllocatedSingleString(ctype);
  freeAllocatedSingleString(vartype);

  delete_lp(lp);

  return 0;
}
