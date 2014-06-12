/////////////////////////////////////////////////////////////////////////////////////////////////////////
// scisymphony: A scilab interface to Symphony library for bicriteria mixed integer linear programming //
/////////////////////////////////////////////////////////////////////////////////////////////////////////

//  Copyright (C) 2008-2010 Yann Collette.
//
//  SciSymphony is free software; you can redistribute it and/or modify it
//  under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 2, or (at your option) any
//  later version.
//
//  This part of code is distributed with the FURTHER condition that it 
//  can be compiled and linked with the Matlab libraries and it can be 
//  used within the Matlab environment.
//
//  SciSymphony is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with Octave; see the file COPYING.  If not, write to the Free
//  Software Foundation, 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#include <cstdio>
#include <cstring>
#include <exception>
#include <iostream>
#include <string.h>

#include <symphony.h>
#include <sym_constants.h>
#include "stream_redirect.h"

extern "C"
{
#include <stack-c.h>
#include <sciprint.h>
#include <Scierror.h>
#include <MALLOC.h>
#include <api_parameters.h>
}

#include <api_scilab.h>

//#define DEBUG 1

#include <helper.hpp>

using namespace std;

#define A_IN        1
#define C1_IN       2
#define C2_IN       3
#define LHS_IN      4
#define RHS_IN      5
#define UPPER_IN    6
#define LOWER_IN    7
#define CTYPE_IN    8
#define VTYPE_IN    9
#define PARAM_IN    10
#define LASTPARAM   10

// Output Arguments
#define  XMIN_OUT       11
#define  FMIN_OUT       12
#define  STATUS_OUT     13
#define  LAMBDA_OUT     14
#define  SYM_STAT_OUT   15
#define  ITERATIONS_OUT 16
#define  EXTRA_OUT      17

// To define the function scisymphony as a C function and allow this function to be loaded easily in scilab
extern "C" int scisymphony(char * fname)
{
  sym_environment * sci_environment = NULL;
  int numrows = -1, numcols = -1;
  int * start = NULL, * index = NULL;
  double * value = NULL;
  double *primal = NULL;
  double *dual   = NULL;
    
  int i, j, nz = 0, Log = 0;
  int count = 0, status = 0;
  int * param_in_addr = NULL;
  SciSparse S_A;
  ScilabStream sci_cout(std::cout);
  ScilabStream sci_cerr(std::cerr);

  if (Rhs<LASTPARAM) 
    {
      Scierror(999,"%s: %d inputs required in call to scisymphony. Bug in symphony.sci ?...\n",fname, LASTPARAM);
      return 0;
    }
		
  /* Get pointers to input */

  int n_a,     m_a,     * a_addr     = NULL;
  int n_c1,    m_c1,    * c1_addr    = NULL;
  int n_c2,    m_c2,    * c2_addr    = NULL;
  int n_lhs,   m_lhs,   * lhs_addr   = NULL;
  int n_rhs,   m_rhs,   * rhs_addr   = NULL;
  int n_upper, m_upper, * upper_addr = NULL;
  int n_lower, m_lower, * lower_addr = NULL;
  int * vtype_addr = NULL;
  int * ctype_addr = NULL;
  double * c1 = NULL, * c2 = NULL, * lhs = NULL, * rhs = NULL, * lower = NULL;
  double * upper = NULL, * a = NULL;
  char * ctype = NULL, * vtype = NULL;
  int type;
  SciErr _SciErr;

  _SciErr = getVarAddressFromPosition(pvApiCtx, C1_IN, &c1_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, c1_addr, &n_c1, &m_c1, &c1); SCICOINOR_ERROR;
  _SciErr = getVarAddressFromPosition(pvApiCtx, C2_IN, &c2_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, c2_addr, &n_c2, &m_c2, &c2); SCICOINOR_ERROR;
  _SciErr = getVarAddressFromPosition(pvApiCtx, LHS_IN, &lhs_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, lhs_addr, &n_lhs, &m_lhs, &lhs); SCICOINOR_ERROR;
  _SciErr = getVarAddressFromPosition(pvApiCtx, RHS_IN, &rhs_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, rhs_addr, &n_rhs, &m_rhs, &rhs); SCICOINOR_ERROR;
  _SciErr = getVarAddressFromPosition(pvApiCtx, LOWER_IN, &lower_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, lower_addr, &n_lower, &m_lower, &lower); SCICOINOR_ERROR;
  _SciErr = getVarAddressFromPosition(pvApiCtx, UPPER_IN, &upper_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, upper_addr, &n_upper, &m_upper, &upper); SCICOINOR_ERROR;
  _SciErr = getVarAddressFromPosition(pvApiCtx, CTYPE_IN, &ctype_addr); SCICOINOR_ERROR;
  getAllocatedSingleString(pvApiCtx, ctype_addr, &ctype); // in 'L' (<=),'E'(=),'G'(>=),'R'(<=+<=),'N'(free)
  _SciErr = getVarAddressFromPosition(pvApiCtx, VTYPE_IN, &vtype_addr); SCICOINOR_ERROR;
  getAllocatedSingleString(pvApiCtx, vtype_addr, &vtype);

  numcols = n_c1  * m_c1;  /* Length of c1 == number of columns */
  numrows = n_rhs * m_rhs; /* length of b == number of rows    */
  if (numrows==0) numrows = n_lhs * m_lhs;

#ifdef DEBUG
  sciprint("DEBUG: c1          size: m = %d n = %d\n", m_c1,    n_c1);
  sciprint("DEBUG: c2          size: m = %d n = %d\n", m_c2,    n_c2);
  sciprint("DEBUG: lhs         size: m = %d n = %d\n", m_lhs,   n_lhs);
  sciprint("DEBUG: rhs         size: m = %d n = %d\n", m_rhs,   n_rhs);
  sciprint("DEBUG: lower       size: m = %d n = %d\n", m_lower, n_lower);
  sciprint("DEBUG: upper       size: m = %d n = %d\n", m_upper, n_upper);
  sciprint("DEBUG: ctype       size: m = %d n = %d\n", m_ctype, n_ctype);
  sciprint("DEBUG: vtype       size: m = %d n = %d\n", m_vtype, n_vtype);
  sciprint("DEBUG: numrows = %d\n", numrows);
  sciprint("DEBUG: numcols = %d\n", numcols);

  sciprint("c1 :");
  for(i=0;i<numcols; i++) sciprint("%f ",*stk(l_c1+i));
  sciprint("\n");
  if ((m_c2!=0)&&(n_c2!=0))
    {
      sciprint("c2 :");
      for(i=0;i<numcols; i++) sciprint("%f ",*stk(l_c2+i));
      sciprint("\n");
    }
  sciprint("lhs :");
  for(i=0;i<numrows; i++) sciprint("%f ",*stk(l_lhs+i));
  sciprint("\n");
  sciprint("rhs :");
  for(i=0;i<numrows; i++) sciprint("%f ",*stk(l_rhs+i));
  sciprint("\n");
  sciprint("lb :");
  for(i=0;i<numcols; i++) sciprint("%f ",*stk(l_lower+i));
  sciprint("\n");
  sciprint("ub :");
  for(i=0;i<numcols; i++) sciprint("%f ",*stk(l_upper+i));
  sciprint("\n");
  sciprint("ctype = %s\n", cstk(l_ctype));
  sciprint("vtype = %s\n", cstk(l_vtype));
#endif

  ////////////////////
  // Initialization //
  ////////////////////

  sci_environment = sym_open_environment();

  /////////////////
  // Get options //
  /////////////////

  // Set default values
  status = sym_set_defaults(sci_environment);

#ifdef DEBUG
  sciprint("DEBUG: get options\n");
#endif

  int     tmp_int, tmp_res;
  double  tmp_double;
  char  * tmp_char;

  initPList(pvApiCtx, PARAM_IN, &param_in_addr);
  if (!checkPList(pvApiCtx, param_in_addr))
    {
      Scierror(999, "%s: argument n° %d is not a plist\n", fname, PARAM_IN);

      return 0;
    }

  ///////////////////////
  // Global parameters //
  ///////////////////////
  // 'verbosity'  integer (0). Sets the verbosity of all modules to the given value. In general, the greater this number the more 
  //                           verbose each module is.
  getIntInPList(pvApiCtx, param_in_addr, "verbosity", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "verbosity",tmp_int);
      
  // 'random_seed'  integer (17). A random seed. 
  getIntInPList(pvApiCtx, param_in_addr, "random_seed", &tmp_int, &tmp_res, 17, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "random_seed",tmp_int);
      
  // 'granularity'  double (1e-6). Should be set to "the minimum difference between two distinct objective function values" 
  //                               less the epsilon tolerance. E.g., if every variable is integral and the objective 
  //                               coefficients are integral then for any feasible solution the objective value is integer,
  //                               so granularity could be correctly set to .99999.
  getDoubleInPList(pvApiCtx, param_in_addr, "granularity", &tmp_double, &tmp_res, 1e-6, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_dbl_param(sci_environment, "granularity",tmp_double);
      
  // 'upper_bound'  double (none) . The value of the best known upper bound.
  getDoubleInPList(pvApiCtx, param_in_addr, "upper_bound", &tmp_double, &tmp_res, 1e100, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_dbl_param(sci_environment, "upper_bound",tmp_double);
      
  // 'lower_bound'  double (no lower bound). This parameter is used if the user wants to artificially impose a lower bound.
  getDoubleInPList(pvApiCtx, param_in_addr, "lower_bound", &tmp_double, &tmp_res, -1e100, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_dbl_param(sci_environment, "lower_bound",tmp_double);

  // 'probname'  string (empty string). The name of the problem name.
  getStringInPList(pvApiCtx, param_in_addr, "probname", &tmp_char, &tmp_res, "test", Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_str_param(sci_environment, "probname",tmp_char);
  FREE(tmp_char);

  // 'infile_name'  string (empty string). The name of the input file that was read by '-F' or the '-L' flag.
  getStringInPList(pvApiCtx, param_in_addr, "infile_name", &tmp_char, &tmp_res, "test.mps", Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_str_param(sci_environment, "infile_name",tmp_char);
  FREE(tmp_char);

  //////////////////////////////
  // Master module parameters //
  //////////////////////////////
  // 'M_verbosity'  integer (0).
  getIntInPList(pvApiCtx, param_in_addr, "M_verbosity", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "M_verbosity",tmp_int);

  // 'M_random_seed'  integer (17). A random seed just for the Master module.
  getIntInPList(pvApiCtx, param_in_addr, "M_random_seed", &tmp_int, &tmp_res, 17, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "M_random_seed",tmp_int);

  // 'upper_bound_estimate'  double (no estimate). This parameter is used if the user wants to provide an estimate
  //                                 of the optimal value which will help guide the search. This is used in
  //                                 conjunction with the diving strategy BEST ESTIMATE.
  getDoubleInPList(pvApiCtx, param_in_addr, "upper_bound_estimate", &tmp_double, &tmp_res, 1e100, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_dbl_param(sci_environment, "upper_bound_estimate",tmp_double);

  // 'tm_exe', 'dg_exe'  strings ('tm', 'dg'). The name of the executable files of the TM and DG modules. Note that 
  //                             the TM executable name may have extensions that depend on the configuration of the
  //                             modules, but the default is always set to the file name produced by the makefile.
  //                             If you change the name of the treemanager executable from the default, you must 
  //                             set this parameter to the new name.
  getStringInPList(pvApiCtx, param_in_addr, "tm_exe", &tmp_char, &tmp_res, "tm", Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_str_param(sci_environment, "tm_exe",tmp_char);
  FREE(tmp_char);
  getStringInPList(pvApiCtx, param_in_addr, "dg_exe", &tmp_char, &tmp_res, "dg", Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_str_param(sci_environment, "dg_exe",tmp_char);
  FREE(tmp_char);

  // 'tm_debug', 'dg_debug'  boolean (both FALSE). Whether these modules should be started under a debugger or not 
  //                             (see 5.6.2 for more details on this).
  getIntInPList(pvApiCtx, param_in_addr, "tm_debug", &tmp_int, &tmp_res, FALSE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "tm_debug",tmp_int);
  getIntInPList(pvApiCtx, param_in_addr, "dg_debug", &tmp_int, &tmp_res, FALSE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "dg_debug",tmp_int);

  // 'tm_machine'  string (empty string). On which processor of the virtual machine the TM should be run. Leaving 
  //                             this parameter as an empty string means arbitrary selection.
  getStringInPList(pvApiCtx, param_in_addr, "tm_machine", &tmp_char, &tmp_res, "", Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_str_param(sci_environment, "tm_machine",tmp_char);
  FREE(tmp_char);

  //// 'do_draw_graph'  boolean (FALSE). Whether to start up the DG module or not (see Section 5.6.4 for an introduction to this).
  // 'do_branch_and_cut'  boolean (TRUE). Whether to run the branch and cut algorithm or not. (Set this to FALSE
  //                             to run the userN"s heuristics only.)
  getIntInPList(pvApiCtx, param_in_addr, "do_branch_and_cut", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "do_branch_and_cut",tmp_int);

  // 'mc_search_order'  integer (MC_FIFO). Use the fifo (MC_FIFO) or lifo (MC_LIFO) searh order during the multi
  //                             criteria solution procedure.
  getIntInPList(pvApiCtx, param_in_addr, "mc_search_order", &tmp_int, &tmp_res, MC_FIFO, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "mc_search_order",tmp_int);

  // mc_warm_start  boolean (FALSE). Whether to solve the corresponding problem of each iteration from a warm start
  //                             loaded from a base iteration (which is the first iteration where gamma = 1.0 and
  //                             tau = 0.0) or from scratch. Currently, this option is supported if only the supported
  //                             solutions are desired to be found.
  getIntInPList(pvApiCtx, param_in_addr, "mc_warm_start", &tmp_int, &tmp_res, FALSE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "mc_warm_start",tmp_int);

  // 'trim_warm_tree'  boolean (FALSE). Whether to trim the warm start tree before re-solving. This consists of
  //                             locating nodes whose descendants are all likely to be pruned in the resolve and
  //                             eliminating those descendants in favor of processing the parent node itself.
  getIntInPList(pvApiCtx, param_in_addr, "trim_warm_tree", &tmp_int, &tmp_res, FALSE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "trim_warm_tree",tmp_int);

  // 'mc_compare_solution_tolerance'  double (0.001). If the difference between the objective values of two 
  //                             solutions to be compared, during the bicriteria solution procedure, are less
  //                             than this tolerance, then assume them to be equal.
  getDoubleInPList(pvApiCtx, param_in_addr, "mc_compare_solution_tolerance", &tmp_double, &tmp_res, 0.001, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_dbl_param(sci_environment, "mc_compare_solution_tolerance",tmp_double);

  // 'mc_binary_search_tolerance'  double (0). The tolerance to be used to differentiate the gamma values if 
  //                             binary search is used during the bicriteria solution procedure. A value greater
  //                             than zero will cause the binary search to be activated.
  getDoubleInPList(pvApiCtx, param_in_addr, "mc_binary_search_tolerance", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_dbl_param(sci_environment, "mc_binary_search_tolerance",tmp_double);

  /////////////////////////////
  // Tree manager parameters //
  /////////////////////////////
  // 'TM_verbosity'  integer (0). The verbosity of the TM module.
  getIntInPList(pvApiCtx, param_in_addr, "TM_verbosity", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "TM_verbosity",tmp_int);

  // 'lp_exe', 'cg_exe', 'cp_exe'  strings ('lp', 'cg', 'cp'). The name of the LP, CG, and CP module binaries. 
  //                             Note: when running in parallel using PVM, these executables (or links to them) 
  //                             must reside in the PVM ROOT/bin/PVM ARCH/ directory. Also, be sure to note that 
  //                             the executable names may have extensions that depend on the configuration of 
  //                             the modules, but the defaults will always be set to the name that the makefile 
  //                             produces.
  getStringInPList(pvApiCtx, param_in_addr, "lp_exe", &tmp_char, &tmp_res, "lp", Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_str_param(sci_environment, "lp_exe",tmp_char);
  FREE(tmp_char);
  getStringInPList(pvApiCtx, param_in_addr, "cg_exe", &tmp_char, &tmp_res, "cg", Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_str_param(sci_environment, "cg_exe",tmp_char);
  FREE(tmp_char);
  getStringInPList(pvApiCtx, param_in_addr, "cp_exe", &tmp_char, &tmp_res, "cp", Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_str_param(sci_environment, "cp_exe",tmp_char);
  FREE(tmp_char);

  // 'lp_debug', 'cg_debug', 'cp_debug'  boolean (all FALSE). Whether the modules should be started under a debugger or not.
  getIntInPList(pvApiCtx, param_in_addr, "lp_debug", &tmp_int, &tmp_res, FALSE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "lp_debug",tmp_int);
  getIntInPList(pvApiCtx, param_in_addr, "cg_debug", &tmp_int, &tmp_res, FALSE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "cg_debug",tmp_int);
  getIntInPList(pvApiCtx, param_in_addr, "cp_debug", &tmp_int, &tmp_res, FALSE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "cp_debug",tmp_int);

  // 'max_active_nodes'  integer (1). The maximum number of active search tree nodesequal to the number of LP and CG 
  //                               tandems to be started up.
  getIntInPList(pvApiCtx, param_in_addr, "max_active_nodes", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "max_active_nodes",tmp_int);

  // 'max_cp_num'  integer (0). The maximum number of cut pools to be used.
  getIntInPList(pvApiCtx, param_in_addr, "max_cp_num", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "max_cp_num",tmp_int);

  // 'lp_mach_num', 'cg_mach_num', 'cp_mach_num'  integers (all 0). The number of processors in the virtual machine 
  //                               to run LP (CG, CP) processes. If this value is 0 then the processes will be assigned
  //                               to processors in round-robin order. Otherwise the next xx mach num lines describe
  //                               the processors where the LP (CG, CP) modules must run. The keyword value pairs on
  //                               these lines must be TM xx machine and the name or IP address of a processor (the processor
  //                               names need not be distinct). In this case the actual processes are assigned in a
  //                               round robin fashion to the processors on this list. This feature is useful if a
  //                               specific software package is needed for some module, but that software is not
  //                               licensed for every node of the virtual machine or if a certain process must run on a
  //                               certain type of machine due to resource requirements.
  getIntInPList(pvApiCtx, param_in_addr, "lp_mach_num", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "lp_mach_num",tmp_int);
  getIntInPList(pvApiCtx, param_in_addr, "cg_mach_num", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "cg_mach_num",tmp_int);
  getIntInPList(pvApiCtx, param_in_addr, "cp_mach_num", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "cp_mach_num",tmp_int);

  // 'use_cg'  boolean (FALSE).     Whether to use a cut generator or not.
  getIntInPList(pvApiCtx, param_in_addr, "use_cg", &tmp_int, &tmp_res, FALSE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "use_cg",tmp_int);

  // 'TM_random_seed'  integer (17). The random seed used in the TM.
  getIntInPList(pvApiCtx, param_in_addr, "TM_random_seed", &tmp_int, &tmp_res, 17, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "TM_random_seed",tmp_int);

  // 'unconditional_dive_frac'  double (0.0). The fraction of the nodes on which SYMPHONY randomly dives 
  //                                unconditionally into one of the children.
  getDoubleInPList(pvApiCtx, param_in_addr, "unconditional_dive_frac", &tmp_double, &tmp_res, 0.0, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_dbl_param(sci_environment, "unconditional_dive_frac",tmp_double);

  // 'diving_strategy'  integer (BEST_ESTIMATE{0}). The strategy employed when deciding whether to dive or
  //                                not. The BEST_ESTIMATE{0} strategy continues to dive until the lower bound
  //                                in the child to be dived into exceeds the parameter upper_bound_estimate, 
  //                                which is given by the user. The COMP_BEST_K{1} strategy computes the average 
  //                                lower bound on the best diving_k search tree nodes and decides to dive if
  //                                the lower bound of the child to be dived into does not exceed this average by
  //                                more than the fraction diving_threshold. The COMP_BEST_K_GAP{2} strategy takes 
  //                                the size of the gap into account when deciding whether to dive. After the average
  //                                lower bound of the best diving_k nodes is computed, the gap between this average 
  //                                lower bound and the current upper bound is computed. Diving only occurs if the
  //                                difference between the computed average lower bound and the lower bound of the 
  //                                child to be dived into is at most the fraction diving_threshold of the gap.
  //                                Note that fractional diving settings can override these strategies. See below.
  getIntInPList(pvApiCtx, param_in_addr, "diving_strategy", &tmp_int, &tmp_res, BEST_ESTIMATE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "diving_strategy",tmp_int);

  // 'diving_k', 'diving_threshold'  integer, double (1, 0.05). See above.
  getIntInPList(pvApiCtx, param_in_addr, "diving_k", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "diving_k",tmp_int);
  getDoubleInPList(pvApiCtx, param_in_addr, "diving_threshold", &tmp_double, &tmp_res, 0.05, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_dbl_param(sci_environment, "diving_threshold",tmp_double);

  // 'fractional_diving_ratio', 'fractional_diving_num'  integer (0.02, 0). Diving occurs automatically if the 
  //                                number of fractional variables in the child to be dived into is less than 
  //                                fractional diving num or the fraction of total variables that are fractional is less
  //                                than fractional_diving_ratio. This overrides the other diving rules. Note that in 
  //                                order for this option to work, the code must be compiled with FRACTIONAL_BRANCHING defined. 
  //                                This is the default. See the makefile for more details.
  getDoubleInPList(pvApiCtx, param_in_addr, "fractional_diving_ratio", &tmp_double, &tmp_res, 0.02, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_dbl_param(sci_environment, "fractional_diving_ratio",tmp_double);
  getIntInPList(pvApiCtx, param_in_addr, "fractional_diving_num", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "fractional_diving_num",tmp_int);

  // 'node_selection_rule'  integer (LOWEST_LP_FIRST{0}). The rule for selecting the next search tree node to be processed.
  //                                This rule selects the one with lowest lower bound. Other possible values are: 
  //                                HIGHEST_LP_FIRST{1}, BREADTH_FIRST_SEARCH{2} and DEPTH_FIRST_SEARCH{3}.
  getIntInPList(pvApiCtx, param_in_addr, "node_selection_rule", &tmp_int, &tmp_res, LOWEST_LP_FIRST, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "node_selection_rule",tmp_int);

  // 'load_balance_level'  integer (-1) A naive attempt at load balancing on problems where significant time is spent in
  //                                the root node, contributing to a lack of parallel speed-up. Only a prescribed number 
  //                                of iterations (load balance iter) are performed in the root node (and in each subsequent
  //                                node on a level less than or equal to load_balance_level) before branching is forced in 
  //                                order to provide additional subproblems for the idle processors to work on. 
  //                                This doesn't work well in general. 
  getIntInPList(pvApiCtx, param_in_addr, "load_balance_level", &tmp_int, &tmp_res, -1, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "load_balance_level",tmp_int);

  // 'load_balance_iter'  integer (-1) Works in tandem with the load_balance_level to attempt some simple load balancing.
  //                                See the above description.
  getIntInPList(pvApiCtx, param_in_addr, "load_balance_iter", &tmp_int, &tmp_res, -1, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "load_balance_iter",tmp_int);

  // 'keep_description_of_pruned'  integer (DISCARD{0}). Whether to keep the description of pruned search tree nodes or not.
  //                                The reasons to do this are (1) if the user wants to write out a proof of optimality using
  //                                the logging function, (2) for debugging, or (3) to get a visual picture of the tree
  //                                using the software VBCTOOL. Otherwise, keeping the pruned nodes around just takes up memory.
  //                                There are three options if it is desired to keep some description of the pruned nodes 
  //                                around. First, their full description can be written out to disk and freed from memory 
  //                                (KEEP_ON_DISK_FULL{1}). There is not really too much you can do with this kind of file, 
  //                                but theoretically, it contains a full record of the solution process and could be used 
  //                                to provide a certificate of optimality (if we were using exact arithmetic) using an 
  //                                independent verifier.
  //                                In this case, the line following keep_description_of_pruned should be a line containing 
  //                                the keyword pruned_node_file_name with its corresponding value being the name of a file
  //                                to which a description of the pruned nodes can be written. The file does not need to exist
  //                                and will be over-written if it does exist.
  //                                If you have the software VBCTOOL, then you can alternatively just write out the information
  //                                VBCTOOL needs to display the tree (KEEP_ON_DISK_VBC_TOOL{2}).
  //                                Finally, the user can set the value to of this parameter to KEEP_IN_MEMORY{2}, 
  //                                in which case all pruned nodes will be kept in memory and written out to the regular log
  //                                file if that option is chosen. This is really only useful for debugging. 
  //                                Otherwise, pruned nodes should be flushed.
  getIntInPList(pvApiCtx, param_in_addr, "keep_description_of_pruned", &tmp_int, &tmp_res, DISCARD, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "keep_description_of_pruned",tmp_int);

  // 'keep_warm_start'  boolean (FALSE). Turning this parameter on will have exactly the same impact with setting the
  //                                keep_description_of_pruned to KEEP_IN_MEMORY{2}. This will allow SYMPHONY to keep all
  //                                the necessary information obtained from the branching tree of the original problem to 
  //                                be able to warm start after a parameter or problem data modification. Thus, if it is
  //                                intended to warm start later, the user should set this parameter before solving the 
  //                                original problem.
  getIntInPList(pvApiCtx, param_in_addr, "keep_warm_start", &tmp_int, &tmp_res, FALSE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "keep_warm_start",tmp_int);

  // 'warm_start_node_limit'  integer (SYM_INFINITY). Setting this parameter will start the warm start routine using
  //                                only the first warm_start_node_limit nodes generated during the previous solve procedure.
  //                                The rest of the tree will be trimmed.
  getIntInPList(pvApiCtx, param_in_addr, "warm_start_node_limit", &tmp_int, &tmp_res, FALSE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "warm_start_node_limit",tmp_int);

  // 'warm_start_node_ratio'  double (0.0). Setting this parameter will start the warm start routine using only the
  //                                first warm_start_node_ratio% of the nodes generated during the previous solve procedure.
  getDoubleInPList(pvApiCtx, param_in_addr, "warm_start_node_ratio", &tmp_double, &tmp_res, 0.0, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_dbl_param(sci_environment, "warm_start_node_ratio",tmp_double);

  // 'warm_start_node_level'  integer (SYM_INFINITY). Setting this parameter will start the warm start routine using
  //                                all the nodes above the level warm_start_node_level of the tree generated during
  //                                the previous solve procedure. The rest of the tree will be trimmed.
  // YC: a problem here: overflow in the implicit conversion of the constant
  getIntInPList(pvApiCtx, param_in_addr, "warm_start_node_level", &tmp_int, &tmp_res, SYM_INFINITY, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "warm_start_node_level",tmp_int);

  // 'warm_start_node_level_ratio'  double (0.0). Setting this parameter will start the warm start routine using all
  //                                the nodes above the level warm_start_node_level% of the warm start tree depth.
  //                                The rest of the tree will be trimmed 
  getDoubleInPList(pvApiCtx, param_in_addr, "warm_start_node_level_ratio", &tmp_double, &tmp_res, 0.0, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_dbl_param(sci_environment, "warm_start_node_level_ratio",tmp_double);

  // 'logging'  integer (NO_LOGGING{0}). Whether or not to write out the state of the search tree and all other necessary
  //                                data to disk periodically in order to allow a warm start in the case of a system crash
  //                                or to allow periodic viewing with VBCTOOL. If the value of this parameter is set to
  //                                FULL_LOGGING{1}, then all information needed to warm start the calculation will written 
  //                                out periodically. The next two lines of the parameter file following should contain
  //                                the keywords tree_log_file_name and cut_log_file_name along with corresponding file names
  //                                as values. These will be the files used to record the search tree and related data and
  //                                the list of cuts needed to reconstruct the tree.
  //                                If the value of the parameter is set to VBC_TOOL{2}, then only the information VBCTOOL 
  //                                needs to display the tree will be logged. 
  //                                This is not really a very useful option since a 'live' picture of the tree can be obtained
  //                                using the vbc_emulation parameter described below.
  getIntInPList(pvApiCtx, param_in_addr, "logging", &tmp_int, &tmp_res, NO_LOGGING, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "logging",tmp_int);

  // 'logging_interval'  integer (1800). Interval (in seconds) between writing out the above log files.
  getIntInPList(pvApiCtx, param_in_addr, "logging_interval", &tmp_int, &tmp_res, 1800, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "logging_interval",tmp_int);

  // 'warm_start' boolean (0). Used to allow the tree manager to make a warm start by reading in previously written log 
  //                                files. If this option is set, then the two line following must start with the 
  //                                keywords warm_start_tree_file_name and warm_start_cut_file_name and include the 
  //                                appropriate file names as the corresponding values.
  getIntInPList(pvApiCtx, param_in_addr, "warm_start", &tmp_int, &tmp_res, FALSE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "warm_start",tmp_int);

  // 'vbc_emulation'  integer (NO_VBC_EMULATION{0}) Determines whether or not to employ the VBCTOOL emulation mode.
  //                                If one of these modes is chosen, then the tree will be displayed in 'real time' 
  //                                using the VBCTOOL Software. When using the option VBC_EMULATION_LIVE{2} and piping
  //                                the output directly to VBCTOOL, the tree will be displayed as it is constructed,
  //                                with color coding indicating the status of each node. 
  //                                With VBC_EMULATION_FILE{1} selected, a log file will be produced which can later be
  //                                read into VBCTOOL to produce an emulation of the solution process at any desired speed.
  //                                If VBC_EMULATION_FILE is selected, the the following line should contain the keyword 
  //                                vbc_emulation_file_name along with the corresponding file name for a value.
  getIntInPList(pvApiCtx, param_in_addr, "vbc_emulation", &tmp_int, &tmp_res, NO_VBC_EMULATION, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "vbc_emulation",tmp_int);

  // 'price_in_root'  boolean (FALSE). Whether to price out variables in the root node before the second phase starts 
  //                                (called repricing the root).
  getIntInPList(pvApiCtx, param_in_addr, "price_in_root", &tmp_int, &tmp_res, NO_VBC_EMULATION, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "price_in_root",tmp_int);

  // 'trim_search_tree'  boolean (FALSE). Whether to trim the search tree before the second phase starts or not. Useful
  //                                only if there are two phases. (It is very useful then.)
  getIntInPList(pvApiCtx, param_in_addr, "trim_search_tree", &tmp_int, &tmp_res, FALSE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "trim_search_tree",tmp_int);

  // 'colgen_in_first_phase', 'colgen_in_second_phase'  integers (both 4). These parameters determine if and when to do
  //                                column generation in the first and second phase of the algorithm. The value of 
  //                                each parameter is obtained by setting the last four bits. The last two bits refer
  //                                to what to do when attempting to prune a node. If neither of the last two bits are
  //                                set, then we donN"t do anything we just prune it. If only the last bit is set,
  //                                then we simply save the node for the second phase without doing any column generation
  //                                (yet). If only the second to last bit is set, then we do column generation 
  //                                immediately and resolve if any new columns are found.
  //                                The next two higher bits determine whether or not to do column generation before 
  //                                branching. If only the third lowest bit is set, then no column generation occurs 
  //                                before branching. If only the fourth lowest bit is set, then column generation is
  //                                attempted before branching. The default is not to generate columns before branching 
  //                                or fathoming, which corresponds to only the third lowest bit being set, resulting in
  //                                a default value of 4.
  getIntInPList(pvApiCtx, param_in_addr, "colgen_in_first_phase", &tmp_int, &tmp_res, 4, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "colgen_in_first_phase",tmp_int);
  getIntInPList(pvApiCtx, param_in_addr, "colgen_in_second_phase", &tmp_int, &tmp_res, 4, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "colgen_in_second_phase",tmp_int);

  // 'time_limit'  double (-1.0). Number of seconds of wall-clock time allowed for solution. When this time limit is 
  //                                reached, the solution process will stop and the best solution found to that point,
  //                                along with other relevant data, will be output. A time limit less than 0.0 means
  //                                there is no limit.
  getDoubleInPList(pvApiCtx, param_in_addr, "time_limit", &tmp_double, &tmp_res, -1.0, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_dbl_param(sci_environment, "time_limit",tmp_double);

  // 'node_limit'  integer (-1). Number of nodes allowed to be analyzed during the solution. When this node limit is
  //                                reached, the solution process will stop and the best solution found to that point,
  //                                along with other relevant data, will be output. A node limit less than 0 means there
  //                                is no limit.
  getIntInPList(pvApiCtx, param_in_addr, "node_limit", &tmp_int, &tmp_res, -1, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "node_limit",tmp_int);

  // 'gap_limit'  double (-1.0). Target gap limit allowed for solution. When the gap between the lower and the upper
  //                                bound reaches this point, the solution process will stop and the best solution
  //                                found to that point, along with other relevant data, will be output. A gap limit
  //                                less than 0 means there is no limit.
  getDoubleInPList(pvApiCtx, param_in_addr, "gap_limit", &tmp_double, &tmp_res, -1.0, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_dbl_param(sci_environment, "gap_limit",tmp_double);

  // 'find_first_feasible'  boolean (FALSE). Whether to stop after finding the first feasible solution or not.
  getIntInPList(pvApiCtx, param_in_addr, "find_first_feasible", &tmp_int, &tmp_res, FALSE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "find_first_feasible",tmp_int);

  // 'sensitivity_analysis'  boolean (FALSE). If the user wants to do the rudimentary sensitivity analysis, which will 
  //                                give a lower bound for the problem modified by the right hand side, then, this 
  //                                parameter has to be set before solving the original problem. If it is set, SYMPHONY
  //                                will keep the necessary information from the solution processes of the original 
  //                                problem to be able to do the sensitivity analysis later.
  getIntInPList(pvApiCtx, param_in_addr, "sensitivity_analysis", &tmp_int, &tmp_res, FALSE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "sensitivity_analysis",tmp_int);

  ///////////////////
  // LP Parameters //
  ///////////////////
  // 'LP_verbosity'  integer (0). Verbosity level of the LP module.
  getIntInPList(pvApiCtx, param_in_addr, "LP_verbosity", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "LP_verbosity",tmp_int);

  // 'set_obj_upper_lim'  boolean (FALSE). Whether to stop solving the LP relaxation when itN"s optimal value is 
  //                                provably higher than the global upper bound. There are some advantages to continuing
  //                                the solution process anyway. For instance, this results in the highest possible lower
  //                                bound. On the other hand, if the matrix is full, this node will be pruned anyway and
  //                                the rest of the computation is pointless. This option should be set at FALSE for column
  //                                generation since the LP dual values may not be reliable otherwise.
  getIntInPList(pvApiCtx, param_in_addr, "set_obj_upper_lim", &tmp_int, &tmp_res, FALSE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "set_obj_upper_lim",tmp_int);

  // 'try_to_recover_from_error'  boolean (TRUE). Indicates what should be done in case the LP solver is unable to solve a
  //                                particular LP relaxation because of numerical problems. It is possible to recover from
  //                                this situation but further results may be suspect. On the other hand, the entire solution
  //                                process can be abandoned.
  getIntInPList(pvApiCtx, param_in_addr, "try_to_recover_from_error", &tmp_int, &tmp_res, TRUE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "try_to_recover_from_error",tmp_int);

  // 'problem_type'  integer (ZERO_ONE_PROBLEM{0}). The type of problem being solved. Other values are INTEGER_PROBLEM{1}
  //                                or MIXED_INTEGER_PROBLEM{2}. (Caution: The mixed-integer option is not well tested.)
  getIntInPList(pvApiCtx, param_in_addr, "problem_type", &tmp_int, &tmp_res, ZERO_ONE_PROBLEM, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "problem_type",tmp_int);

  // 'cut_pool_check_frequency'  integer (10). The number of iterations between sending LP solutions to the cut pool to 
  //                                find violated cuts. It is not advisable to check the cut pool too frequently as the
  //                                cut pool module can get bogged down and the LP solution generally do not change that 
  //                                drastically from one iteration to the next anyway.
  getIntInPList(pvApiCtx, param_in_addr, "cut_pool_check_frequency", &tmp_int, &tmp_res, 10, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "cut_pool_check_frequency",tmp_int);

  // 'not_fixed_storage_size'  integer (2048). The not fixed list is a partial list of indices of variables not in the
  //                                matrix that have not been fixed by reduced cost. Keeping this list allows SYMPHONY 
  //                                to avoid repricing variables (an expensive operation) that are not in the matrix 
  //                                because they have already been permanently fixed. When this array reaches its maximum
  //                                size, no more variable indices can be stored. It is therefore advisable to keep the
  //                                maximum size of this array as large as possible, given memory limitations.
  getIntInPList(pvApiCtx, param_in_addr, "not_fixed_storage_size", &tmp_int, &tmp_res, 2048, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "not_fixed_storage_size",tmp_int);

  // 'max_non_dual_feas_to_add_min',
  // 'max_non_dual_feas_to_add_max',
  // 'max_non_dual_feas_to_add_frac'  integer, integer, double (20, 200, .05). These three parameters determine the maximum
  //                                number of non-dual-feasible columns that can be added in any one iteration after pricing.
  //                                This maximum is set to the indicated fraction of the current number of active columns 
  //                                unless this numbers exceeds the given maximum or is less than the given minimum, in which
  //                                case, it is set to the max or min, respectively.
  getIntInPList(pvApiCtx, param_in_addr, "max_non_dual_feas_to_add_min", &tmp_int, &tmp_res, 20, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "max_non_dual_feas_to_add_min",tmp_int);
  getIntInPList(pvApiCtx, param_in_addr, "max_non_dual_feas_to_add_max", &tmp_int, &tmp_res, 200, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "max_non_dual_feas_to_add_max",tmp_int);
  getDoubleInPList(pvApiCtx, param_in_addr, "max_non_dual_feas_to_add_frac", &tmp_double, &tmp_res, 0.05, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_dbl_param(sci_environment, "max_non_dual_feas_to_add_frac",tmp_double);

  // 'max_not_fixable_to_add_min',
  // 'max_not_fixable_to_add_max',
  // 'max_not_fixable_to_add_frac'  integer, integer, double (100, 500, .1) As above, these three parameters determine the 
  //                                maximum number of new columns to be added to the problem because they cannot be priced out.
  //                                These variables are only added when trying to restore infeasibility and usually, this does 
  //                                not require many variables anyway.
  getIntInPList(pvApiCtx, param_in_addr, "max_not_fixable_to_add_min", &tmp_int, &tmp_res, 100, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "max_not_fixable_to_add_min",tmp_int);
  getIntInPList(pvApiCtx, param_in_addr, "max_not_fixable_to_add_max", &tmp_int, &tmp_res, 500, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "max_not_fixable_to_add_max",tmp_int);
  getDoubleInPList(pvApiCtx, param_in_addr, "max_not_fixable_to_add_frac", &tmp_double, &tmp_res, 0.1, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_dbl_param(sci_environment, "max_not_fixable_to_add_frac",tmp_double);

  // 'mat_col_compress_num', 'mat_col_compress_ratio'  integer, double (50, .05). Determines when the matrix should be 
  //                                physically compressed. This only happens when the number of columns is high enough to
  //                                make it "worthwhile." The matrix is physically compressed when the number of deleted
  //                                columns exceeds either an absolute number and a specified fraction of the current number
  //                                of active columns.
  getIntInPList(pvApiCtx, param_in_addr, "mat_col_compress_num", &tmp_int, &tmp_res, 50, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "mat_col_compress_num",tmp_int);
  getDoubleInPList(pvApiCtx, param_in_addr, "mat_col_compress_ratio", &tmp_double, &tmp_res, 0.05, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_dbl_param(sci_environment, "mat_col_compress_ratio",tmp_double);

  // 'mat_row_compress_num', 'mat_row_compress_ratio'  integer, double (20, .05). Same as above except for rows.
  getIntInPList(pvApiCtx, param_in_addr, "mat_row_compress_num", &tmp_int, &tmp_res, 20, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "mat_row_compress_num",tmp_int);
  getDoubleInPList(pvApiCtx, param_in_addr, "mat_row_compress_ratio", &tmp_double, &tmp_res, 0.05, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_dbl_param(sci_environment, "mat_row_compress_ratio",tmp_double);

  // 'tailoff_gap_backsteps', 'tailoff_gap_frac'  integer, double (2, .99). Determines when tailoff is detected in 
  //                                the LP module. Tailoff is reported if the average ratio of the current gap to the 
  //                                previous iteration's gap over the last tailoff_gap_backsteps iterations wasn't at least
  //                                tailoff_gap_frac.
  getIntInPList(pvApiCtx, param_in_addr, "tailoff_gap_backsteps", &tmp_int, &tmp_res, 2, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "tailoff_gap_backsteps",tmp_int);
  getDoubleInPList(pvApiCtx, param_in_addr, "tailoff_gap_frac", &tmp_double, &tmp_res, 0.99, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_dbl_param(sci_environment, "tailoff_gap_frac",tmp_double);

  // 'tailoff_obj_backsteps', 'tailoff_obj_frac'  integer, double (2, .99). Same as above, only the ratio is taken with
  //                                respect to the change in objective function values instead of the change in the gap.
  getIntInPList(pvApiCtx, param_in_addr, "tailoff_obj_backsteps", &tmp_int, &tmp_res, 2, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "tailoff_obj_backsteps",tmp_int);

  // 'ineff_cnt_to_delete'  integer (0). Determines after how many iterations of being deemed ineffective a constraint 
  //                                is removed from the current relaxation.
  getIntInPList(pvApiCtx, param_in_addr, "ineff_cnt_to_delete", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "ineff_cnt_to_delete",tmp_int);

  // 'eff_cnt_before_cutpool'  integer (3). Determines after how many iterations of being deemed effective each cut will
  //                                be sent to the global pool.
  getIntInPList(pvApiCtx, param_in_addr, "eff_cnt_before_cutpool", &tmp_int, &tmp_res, 3, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "eff_cnt_before_cutpool",tmp_int);

  // 'ineffective_constraints'  integer (BASIC_SLACKS_ARE_INEFFECTIVE{2}). Determines under what condition a constraint 
  //                                is deemed ineffective in the current relaxation. Other possible values are 
  //                                NO_CONSTRAINT_IS_INEFFECTIVE{0}, NONZERO_SLACKS_ARE_INEFFECTIVE{1}, and 
  //                                ZERO_DUAL_VALUES_ARE_INEFFECTIVE{3}.
  getIntInPList(pvApiCtx, param_in_addr, "ineffective_constraints", &tmp_int, &tmp_res, BASIC_SLACKS_ARE_INEFFECTIVE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "ineffective_constraints",tmp_int);

  // 'base_constraints_always_effective'  boolean (TRUE). Determines whether the base constraints can ever be removed from
  //                                the relaxation. In some case, removing the base constraints from the problem can be
  //                                disastrous depending on the assumptions made by the cut generator. 
  getIntInPList(pvApiCtx, param_in_addr, "base_constraints_always_effective", &tmp_int, &tmp_res, TRUE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "base_constraints_always_effective",tmp_int);

  // 'branch_on_cuts'  boolean (FALSE). This informs the framework whether the user plans on branching on cuts or not.
  //                                If so, there is additional bookkeeping to be done, such as maintaining a pool of
  //                                slack cuts to be used for branching. Therefore, the user should not set this flag
  //                                unless he actually plans on using this feature. 
  getIntInPList(pvApiCtx, param_in_addr, "branch_on_cuts", &tmp_int, &tmp_res, FALSE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "branch_on_cuts",tmp_int);

  // 'discard_slack_cuts'  integer (DISCARD_SLACKS_BEFORE_NEW_ITERATION{0}). Determines when the pool of slack cuts is
  //                                discarded. The other option is DISCARD_SLACKS_WHEN_STARTING_NEW_NODE{1}.
  getIntInPList(pvApiCtx, param_in_addr, "discard_slack_cuts", &tmp_int, &tmp_res, DISCARD_SLACKS_BEFORE_NEW_ITERATION, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "discard_slack_cuts",tmp_int);

  // 'first_lp_first_cut_time_out',
  // 'first_lp_all_cuts_time_out',
  // 'later_lp_first_cut_time_out',
  // 'later_lp_all_cuts_time_out'  double (0, 0, 5, 1). The next group of parameters determines when the LP should give 
  //                                up waiting for cuts from the cut generator and start to solve the relaxation in its
  //                                current form or possibly branch if necessary. There are two factors that contribute 
  //                                to determining this timeout. First is whether this is the first LP in the search node
  //                                of whether it is a later LP. Second is whether any cuts have been added already in this
  //                                iteration. The four timeout parameters correspond to the four possible combinations of 
  //                                these two variables.
  getDoubleInPList(pvApiCtx, param_in_addr, "first_lp_first_cut_time_out", &tmp_double, &tmp_res, 0.0, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_dbl_param(sci_environment, "first_lp_first_cut_time_out",tmp_double);
  getDoubleInPList(pvApiCtx, param_in_addr, "first_lp_all_cuts_time_out", &tmp_double, &tmp_res, 0.0, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_dbl_param(sci_environment, "first_lp_all_cuts_time_out",tmp_double);
  getDoubleInPList(pvApiCtx, param_in_addr, "later_lp_first_cut_time_out", &tmp_double, &tmp_res, 5.0, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_dbl_param(sci_environment, "later_lp_first_cut_time_out",tmp_double);
  getDoubleInPList(pvApiCtx, param_in_addr, "later_lp_all_cuts_time_out", &tmp_double, &tmp_res, 1.0, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_dbl_param(sci_environment, "later_lp_all_cuts_time_out",tmp_double);

  // 'no_cut_timeout'  This keyword does not have an associated value. If this keyword appears on a line by itself or with
  //                                a value, this tells the framework not to time out while waiting for cuts. This is useful
  //                                for debugging since it enables runs with a single LP module to be duplicated.
  getIntInPList(pvApiCtx, param_in_addr, "no_cut_timeout", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "no_cut_timeout",tmp_int);

  // 'all_cut_timeout'  double (no default). This keyword tells the framework to set all of the above timeout parameters 
  //                                to the value indicated.
  getDoubleInPList(pvApiCtx, param_in_addr, "all_cut_timeout", &tmp_double, &tmp_res, 1.0, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_dbl_param(sci_environment, "all_cut_timeout",tmp_double);

  // 'max_cut_num_per_iter'  integer (20). The maximum number of cuts that can be added to the LP in an iteration. The
  //                                remaining cuts stay in the local pool to be added in subsequent iterations, if they 
  //                                are strong enough.
  getIntInPList(pvApiCtx, param_in_addr, "max_cut_num_per_iter", &tmp_int, &tmp_res, 20, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "max_cut_num_per_iter",tmp_int);

  // 'do_reduced_cost_fixing'  boolean (FALSE). Whether or not to attempt to fix variables by reduced cost. This option
  //                                is highly recommended
  getIntInPList(pvApiCtx, param_in_addr, "do_reduce_cost_fixing", &tmp_int, &tmp_res, FALSE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "do_reduce_cost_fixing",tmp_int);

  // 'gap_as_ub_frac', 'gap_as_last_gap_frac'  double (.1, .7). Determines when reduced cost fixing should be attempted.
  //                                It is only done when the gap is within the fraction gap as ub frac of the upper bound
  //                                or when the gap has decreased by the fraction gap as last gap frac since the last
  //                                time variables were fixed.
  getDoubleInPList(pvApiCtx, param_in_addr, "gap_as_ub_frac", &tmp_double, &tmp_res, 0.1, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_dbl_param(sci_environment, "gap_as_ub_frac",tmp_double);
  getDoubleInPList(pvApiCtx, param_in_addr, "gap_as_last_gap_frac", &tmp_double, &tmp_res, 0.7, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_dbl_param(sci_environment, "gap_as_last_gap_frac",tmp_double);

  // 'do_logical_fixing'  boolean (FALSE). Determines whether the user's logical fixing routine should be used.
  getIntInPList(pvApiCtx, param_in_addr, "do_logical_fixing", &tmp_int, &tmp_res, FALSE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "do_logical_fixing",tmp_int);

  // 'fixed_to_ub_before_logical_fixing',
  // 'fixed_to_ub_frac_before_logical_fixing'  integer, double (1, .01). Determines when logical fixing should be
  //                                attempted. It will be called only when a certain absolute number and a certain number
  //                                of variables have been fixed to their upper bounds by reduced cost. This is because
  //                                it is typically only after fixing variables to their upper bound that other variables
  //                                can be logically fixed.
  getIntInPList(pvApiCtx, param_in_addr, "fixed_to_ub_before_logical_fixing", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "fixed_to_ub_before_logical_fixing",tmp_int);
  getDoubleInPList(pvApiCtx, param_in_addr, "fixed_to_ub_frac_before_logical_fixing", &tmp_double, &tmp_res, 0.01, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_dbl_param(sci_environment, "fixed_to_ub_frac_before_logical_fixing",tmp_double);

  // 'max_presolve_iter'  integer (10). Number of simplex iterations to be performed in the presolve for strong branching.
  getIntInPList(pvApiCtx, param_in_addr, "max_presolve_iter", &tmp_int, &tmp_res, 10, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "max_presolve_iter",tmp_int);

  // 'strong_branching_cand_num_max',
  // 'strong_branching_cand_num_min',
  // 'strong_branching_red_ratio'  integer (10, 5, 1). These three parameters together determine the number of strong branching
  //                               candidates to be used by default. In the root node, strong_branching_cand_num_max candidates
  //                               are used. On each succeeding level, this number is reduced by the number 
  //                               strong_branching_red_ratio multiplied by the square of the level. This continues until the number
  //                               of candidates is reduced to strong_branching_cand_num_min and then that number of candidates 
  //                               is used in all lower levels of the tree.
  getIntInPList(pvApiCtx, param_in_addr, "strong_branching_cand_num_max", &tmp_int, &tmp_res, 10, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "strong_branching_cand_num_max",tmp_int);
  getIntInPList(pvApiCtx, param_in_addr, "strong_branching_cand_num_min", &tmp_int, &tmp_res, 5, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "strong_branching_cand_num_min",tmp_int);
  getIntInPList(pvApiCtx, param_in_addr, "strong_branching_red_ratio", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "strong_branching_red_ratio",tmp_int);

  // 'is_feasible_default'  integer (TEST_INTEGRALITY{1}). Determines the default test to be used to determine feasibility.
  //                               This parameter is provided so that the user can change the default behavior without 
  //                               recompiling. The only other option is TEST_ZERO_ONE{0}. 
  getIntInPList(pvApiCtx, param_in_addr, "is_feasible_default", &tmp_int, &tmp_res, TEST_INTEGRALITY, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "is_feasible_default",tmp_int);

  // 'send_feasible_solution_default'  integer (SEND_NONZEROS{0}). Determines the form in which to send the feasible solution.
  //                               This parameter is provided so that the user can change the default behavior without 
  //                               recompiling. This is currently the only option.
  getIntInPList(pvApiCtx, param_in_addr, "send_feasible_solution_default", &tmp_int, &tmp_res, SEND_NONZEROS, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "send_feasible_solution_default",tmp_int);

  // 'send_lp_solution_default'  integer (SEND_NONZEROS{0}). Determines the default form in which to send the LP solution to
  //                               the cut generator and cut pool. This parameter is provided so that the user can change the
  //                               default behavior without recompiling. The other option is SEND_FRACTIONS{1}.
  getIntInPList(pvApiCtx, param_in_addr, "send_lp_solution_default", &tmp_int, &tmp_res, SEND_NONZEROS, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "send_lp_solution_default",tmp_int);

  // 'display_solution_default'  integer (DISP_NOTHING{0}). Determines how to display the current LP solution if desired. 
  //                               See the description of user display solution() for other possible values. This parameter
  //                               is provided so that the user can change the default behavior without recompiling.
  getIntInPList(pvApiCtx, param_in_addr, "display_solution_default", &tmp_int, &tmp_res, DISP_NOTHING, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "display_solution_default",tmp_int);

  // 'shall_we_branch_default'  integer (USER_BRANCH_IF_MUST{2}). Determines the default branching behavior. Other values 
  //                               are USER__DO_NOT_BRANCH{0} (not recommended as a default), USER__DO_BRANCH{1} 
  //                               (also not recommended as a default), and USER__BRANCH_IF_TAILOFF{3}. This parameter 
  //                               is provided so that the user can change the default behavior without recompiling.
  getIntInPList(pvApiCtx, param_in_addr, "shall_we_branch_default", &tmp_int, &tmp_res, USER__BRANCH_IF_MUST, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "shall_we_branch_default",tmp_int);

  // 'select_candidates_default'  integer (USER__CLOSE_TO_HALF_AND_EXPENSIVE{10}). Determines the default rule for selecting
  //                               strong branching candidates. Other values are USER__CLOSE_TO_HALF{10} and 
  //                               USER__CLOSE_TO_ONE_AND_CHEAP{12}. This parameter is provided so that the user can change the
  //                               default behavior without recompiling.
  getIntInPList(pvApiCtx, param_in_addr, "select_candidate_default", &tmp_int, &tmp_res, USER__CLOSE_TO_HALF_AND_EXPENSIVE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "select_candidate_default",tmp_int);

  // 'compare_candidates_default'  integer (HIGHEST_LOW_OBJ{2}). Determines the default rule for comparing candidates. See the
  //                               description of user_compare_candidates() for other values. This parameter is provided so that
  //                               the user can change the default behavior without recompiling.
  getIntInPList(pvApiCtx, param_in_addr, "compare_candidates_default", &tmp_int, &tmp_res, HIGHEST_LOW_OBJ, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "compare_candidates_default",tmp_int);

  // 'select_child_default'  integer (PREFER_LOWER_OBJ_VALUE{0}). Determines the default rule for selecting the child to be
  //                               processed next. For other possible values, see the description user_select_child(). 
  //                               This parameter is provided so that the user can change the default behavior without recompiling.
  getIntInPList(pvApiCtx, param_in_addr, "select_child_default", &tmp_int, &tmp_res, PREFER_LOWER_OBJ_VALUE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "select_child_default",tmp_int);

  // 'mc_find_supported_solutions'  boolean (FALSE). By default, sym_mc_solve routine will find all the non-dominated solutions
  //                               if the problem to be solved is a bicriteria problem. However, if the user plans to find only
  //                               the supported solutions, then, this parameter has to be set before calling sym_mc_solve routine.
  getIntInPList(pvApiCtx, param_in_addr, "mc_find_supported_solutions", &tmp_int, &tmp_res, FALSE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "mc_find_supported_solutions",tmp_int);

  // 'mc_rho'  double (0.00001). The value used in augmented Chebyshev norm during the bicriteria solution procedure.
  getDoubleInPList(pvApiCtx, param_in_addr, "mc_rho", &tmp_double, &tmp_res, 0.00001, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_dbl_param(sci_environment, "mc_rho",tmp_double);

  // 'generate_cgl_cuts'  boolean (TRUE). Whether or not to generate cuts using COIN's cut generation library. Note that,
  //                               to use CGL cuts, OSI interface has to be used and moreover the corresponding flags 
  //                               have to be set during installation. See the makefile for more details.
  getIntInPList(pvApiCtx, param_in_addr, "generate_cgl_cuts", &tmp_int, &tmp_res, TRUE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "generate_cgl_cuts",tmp_int);

  // 'generate_cgl_gomory_cuts'  boolean (TRUE). Whether or not to generate Gomory cuts using COIN's cut generation library.
  getIntInPList(pvApiCtx, param_in_addr, "generate_cgl_gomory_cuts", &tmp_int, &tmp_res, TRUE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "generate_cgl_gomory_cuts",tmp_int);

  // 'generate_cgl_knapsack_cuts'  boolean (TRUE). Whether or not to generate knapsack cover cuts using COIN's cut generation library.
  getIntInPList(pvApiCtx, param_in_addr, "generate_cgl_knapsack_cuts", &tmp_int, &tmp_res, TRUE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "generate_cgl_knapsack_cuts",tmp_int);

  // 'generate_cgl_oddhole_cuts'  boolean (TRUE). Whether or not to generate generalized odd hole cuts using COIN's cut generation library.
  getIntInPList(pvApiCtx, param_in_addr, "generate_cgl_oddhole_cuts", &tmp_int, &tmp_res, TRUE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "generate_oddhole_cgl_cuts",tmp_int);

  // 'generate_cgl_probing_cuts'  boolean (TRUE). Whether or not to generate probing cuts using COIN's cut generation library.
  getIntInPList(pvApiCtx, param_in_addr, "generate_cgl_probing_cuts", &tmp_int, &tmp_res, TRUE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "generate_cgl_probing_cuts",tmp_int);

  // 'generate_cgl_clique_cuts'  boolean (TRUE). Whether or not to generate clique cuts using COIN's cut generation library.
  getIntInPList(pvApiCtx, param_in_addr, "generate_cgl_clique_cuts", &tmp_int, &tmp_res, TRUE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "generate_cgl_clique_cuts",tmp_int);

  // 'generate_cgl_flow_and_cover_cuts'  boolean (FALSE). Whether or not to generate flow and cover cuts using COIN's cut
  //                               generation library.
  getIntInPList(pvApiCtx, param_in_addr, "generate_cgl_flow_and_cover_cuts", &tmp_int, &tmp_res, FALSE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "generate_cgl_flow_and_cover_cuts",tmp_int);

  // 'generate_cgl_rounding_cuts'  boolean (FALSE). Whether or not to generate simple rounding cuts using COIN's cut generation library.
  getIntInPList(pvApiCtx, param_in_addr, "generate_cgl_rounding_cuts", &tmp_int, &tmp_res, FALSE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "generate_cgl_rounding_cuts",tmp_int);

  // 'generate_cgl_lift_and_project_cuts'  boolean (FALSE). Whether or not to generate lift-andproject cuts using COIN's
  //                               cut generation library.
  getIntInPList(pvApiCtx, param_in_addr, "generate_cgl_lift_and_project_cuts", &tmp_int, &tmp_res, FALSE, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "generate_cgl_lift_and_project_cuts",tmp_int);

  //////////////////////////////
  // Cut generator parameters //
  //////////////////////////////
  // 'CG_verbosity'  integer (0). Verbosity level for the cut generator module.
  getIntInPList(pvApiCtx, param_in_addr, "CG_verbosity", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "CG_verbosity",tmp_int);

  /////////////////////////
  // Cut pool parameters //
  /////////////////////////
  // 'CP_verbosity'  integer (0). Verbosity of the cut pool module.
  getIntInPList(pvApiCtx, param_in_addr, "CP_verbosity", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "CP_verbosity",tmp_int);

  // 'cp_logging'  boolean (0). Determines whether the logging option is enabled. In this case, the entire contents of the
  //                               cut pool are written out periodically to disk (at the same interval as the tree manager
  //                               log files are written). If this option is set, then the line following must start with
  //                               the keyword cp_log_file_name and include the appropriate file name as the value.
  getIntInPList(pvApiCtx, param_in_addr, "cp_logging", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "cp_logging",tmp_int);

  // 'cp_warm_start'  boolean (0). Used to allow the cut pool to make a warm start by reading in a previously written log
  //                               file. If this option is set, then the line following must start with the keyword
  //                               cp_warm_start_file_name and include the appropriate file name as the value.
  getIntInPList(pvApiCtx, param_in_addr, "cp_warm_start", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "cp_warm_start",tmp_int);

  // 'block_size'  integer (5000). Indicates the size of the blocks to allocate when more space is needed in the cut list.
  getIntInPList(pvApiCtx, param_in_addr, "block_size", &tmp_int, &tmp_res, 5000, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "block_size",tmp_int);

  // 'max_size'  integer (2000000). Indicates the maximum size of the cut pool in bytes. This is the total memory taken up
  //                               by the cut list, including all data structures and the array of pointers itself.
  getIntInPList(pvApiCtx, param_in_addr, "max_size", &tmp_int, &tmp_res, 2000000, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "max_size",tmp_int);

  // 'max_number_of_cuts'  integer (10000). Indicates the maximum number of cuts allowed to be stored. When this max is
  //                               reached, cuts are forcibly purged, starting with duplicates and then those indicated 
  //                               by the parameter delete_which (see below), until the list is below the allowable size.
  getIntInPList(pvApiCtx, param_in_addr, "max_number_of_cuts", &tmp_int, &tmp_res, 10000, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "max_number_of_cuts",tmp_int);

  // 'min_to_delete'  integer (1000). Indicates the number of cuts required to be deleted when the pool reaches it's maximum size.
  getIntInPList(pvApiCtx, param_in_addr, "min_to_delete", &tmp_int, &tmp_res, 1000, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "min_to_delete",tmp_int);

  // 'touches_until_deletion'  integer (10). When using the number of touches a cut has as a measure of its quality, this
  //                               parameter indicates the number of touches a cut can have before being deleted from the
  //                               pool. The number of touches is the number of times in a row that a cut has been checked
  //                               without being found to be violated. It is a measure of a cut's relevance or effectiveness.
  getIntInPList(pvApiCtx, param_in_addr, "touches_until_deletion", &tmp_int, &tmp_res, 10, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "touches_until_deletion",tmp_int);

  // 'delete_which'  integer (DELETE_BY_TOUCHES{2}). Indicates which cuts to delete when purging the pool. DELETE_BY_TOUCHES
  //                               indicates that cuts whose number of touches is above the threshold (see 
  //                               touches_until_deletion above) should be purged if the pool gets too large. 
  //                               DELETE_BY_QUALITY{1} indicates that a user-defined measure of quality should be used
  //                               (see the function user_check_cuts in Section6.3.4).
  getIntInPList(pvApiCtx, param_in_addr, "delete_which", &tmp_int, &tmp_res, DELETE_BY_TOUCHES, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "delete_which",tmp_int);

  // 'check_which'  integer (CHECK_ALL_CUTS{0}). Indicates which cuts should be checked for violation. The choices are to check
  //                               all cuts (CHECK_ALL_CUTS{0}); only those that have number of touches below the threshold
  //                               (CHECK_TOUCHES{2}); only those that were generated at a level higher in the tree than the
  //                               current one (CHECK_LEVEL{1}); or both (CHECK_LEVEL_AND_TOUCHES{3}). Note that with 
  //                               CHECK_ALL_CUTS set, SYMPHONY will still only check the first cuts to check cuts in the list
  //                               ordered by quality (see the function user_check_cut).
  getIntInPList(pvApiCtx, param_in_addr, "check_which", &tmp_int, &tmp_res, CHECK_ALL_CUTS, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "check_which",tmp_int);

  // 'cuts_to_check'  integer (1000). Indicates how many cuts in the pool to actually check. The list is ordered by quality
  //                               and the first cuts to check cuts are checked for violation.
  getIntInPList(pvApiCtx, param_in_addr, "cuts_to_check", &tmp_int, &tmp_res, 1000, Log, CHECK_NONE);
  if (tmp_res!=-1) sym_set_int_param(sci_environment, "cuts_to_check",tmp_int);

  //////////////////
  // The A matrix //
  //////////////////

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
	  freeAllocatedSingleString(ctype);
	  freeAllocatedSingleString(vtype);
	  
	  return 0;
	}
      
      // The matrix has been transposed to ease to transfert between scilab and symphony
      numrows = n_a;
      numcols = m_a;
      nz = numrows*numcols;

      start = (int *)MALLOC((numcols+1)*sizeof(int));
      index = (int *)MALLOC(nz*sizeof(int));
      value = (double *)MALLOC(nz*sizeof(double));
      
      start[numcols] = nz;

      count = 0;
      for(i=0; i<numcols; i++)
	{
	  start[i] = count;
	  for(j=0; j<numrows; j++)
	    {
	      index[count] = j;
	      value[count] = *(a+i+j*numcols);
	      count++;
#ifdef DEBUG
	      sciprint("%f ",*(a+i+j*numcols));
#endif
	    }
#ifdef DEBUG
	  sciprint("\n");
#endif
	}
    }
  else
    {
#ifdef DEBUG
      sciprint("DEBUG: A_IN is sparse\n");
#endif
      
      getAllocatedSparseMatrix(pvApiCtx, a_addr, &S_A.m, &S_A.n, &S_A.nel, &S_A.mnel, &S_A.icol, &S_A.R);

#ifdef DEBUG
      sciprint("A' = [%d,%d]\n", m_a, n_a);
      int _dbg_max = -S_A.nel, _dbg_min = S_A.nel;
#endif

      // The matrix has been transposed to ease to transfert between scilab and symphony
      numrows = S_A.n;
      numcols = S_A.m;
      nz = S_A.nel;

      start = (int *)MALLOC((numcols+1)*sizeof(int));
      index = (int *)MALLOC(nz*sizeof(int));
      value = (double *)MALLOC(nz*sizeof(double));

      start[numcols] = nz;

      count = 0;
      for(i=0;i<S_A.m;i++)
	{
	  start[i] = count;

	  if (S_A.mnel[i]!=0) 
	    {
#ifdef DEBUG
	      sciprint("mnel[%d] = %d - start[%d] = %d - ",i, S_A.mnel[i], i, start[i]);
#endif
	      for(j=0;j<S_A.mnel[i];j++)
		{
		  index[count] = S_A.icol[count]-1;
		  value[count] = S_A.R[count];
		  count++;
#ifdef DEBUG
		  if (_dbg_max<index[count-1]) _dbg_max = index[count-1];
		  if (_dbg_min>index[count-1]) _dbg_min = index[count-1];

		  sciprint("[%d] = %f ", S_A.icol[count-1]-1, S_A.R[count-1]);
#endif
		}
#ifdef DEBUG
	      sciprint("\n");
#endif
	    }
	}
#ifdef DEBUG
      sciprint("start[%d] = %d - ",numcols, start[numcols]);
      sciprint("_dbg_max  = %d - ",_dbg_max);
      sciprint("_dbg_min  = %d - ",_dbg_min);
      sciprint("count     = %d - ",count);
#endif

      freeAllocatedSparseMatrix(S_A.mnel, S_A.icol, S_A.R);
    }

#ifdef DEBUG
  sciprint("S_A.m = %d S_A.n = %d S_A.nel = %d numrows = %d numcols = %d - A matrix loaded\n",S_A.m, S_A.n, S_A.nel, numrows, numcols);
#endif

  ////////////////////////////
  // Deal with the variable //
  ////////////////////////////

  char * vartype = (char *)MALLOC(numcols*sizeof(char));
  for(i=0;i<numcols;i++)
    {
      if (*(vtype+i)=='C') vartype[i] = FALSE;
      else                 vartype[i] = TRUE;
    }
#ifdef DEBUG
  sciprint("vartype = ");
  for(i=0;i<numcols;i++)
    {
      sciprint("%d ", (int)(vartype[i]));
    }
  sciprint("\n");
#endif

  //////////////////////
  // Load the problem //
  //////////////////////

#ifdef DEBUG
  sciprint("numcols = %d numrows = %d m_lower = %d n_lower = %d m_upper = %d n_upper = %d \n", numcols, numrows, m_lower, n_lower, m_upper, n_upper);
  sciprint("m_c1 = %d n_c1 = %d m_c2 = %d n_c2 = %d len ctype = %d len vtype = %d\n", m_c1, n_c1, m_c2, n_c2, strlen(ctype), strlen(vtype));
  sciprint("m_rhs = %d n_rhs = %d m_lhs = %d n_lhs = %d\n",m_rhs, n_rhs, m_lhs, n_lhs);
  sciprint("vtype = %s\n", vtype);
  sciprint("ctype = %s\n", ctype);
#endif

  double * in_upper = NULL, * in_lower = NULL;
  double * in_rhs = NULL, * in_lhs = NULL;

  in_upper = (double *)MALLOC(numcols*sizeof(double));
  in_lower = (double *)MALLOC(numcols*sizeof(double));
  in_rhs   = (double *)MALLOC(numrows*sizeof(double));
  in_lhs   = (double *)MALLOC(numrows*sizeof(double));

  for(i=0;i<numcols;i++)
    {
      in_upper[i] = upper[i];
      in_lower[i] = lower[i];
      if (in_upper[i]>SYM_INFINITY)  in_upper[i] = SYM_INFINITY;
      if (in_upper[i]<-SYM_INFINITY) in_upper[i] = -SYM_INFINITY;
      if (in_lower[i]>SYM_INFINITY)  in_lower[i] = SYM_INFINITY;
      if (in_lower[i]<-SYM_INFINITY) in_lower[i] = -SYM_INFINITY;
    }

  for(i=0;i<numrows;i++)
    {
      in_rhs[i] = rhs[i];
      in_lhs[i] = lhs[i];
      if (in_rhs[i]>SYM_INFINITY)  in_rhs[i] = SYM_INFINITY;
      if (in_rhs[i]<-SYM_INFINITY) in_rhs[i] = -SYM_INFINITY;
      if (in_lhs[i]>SYM_INFINITY)  in_lhs[i] = SYM_INFINITY;
      if (in_lhs[i]<-SYM_INFINITY) in_lhs[i] = -SYM_INFINITY;
    }

#ifdef DEBUG
  for(i=0;i<numrows;i++)
    {
      sciprint("row %d - %f - %f == %f - %f \n", i, in_lhs[i], in_rhs[i], lhs[i], rhs[i]);
    }

  for(i=0;i<numcols;i++)
    {
      sciprint("col %d - %f - %f == %f - %f\n", i, in_lower[i], in_upper[i], lower[i], upper[i]);
    }
#endif

  if ((m_c2==0)&&(n_c2==0))
    {
#ifdef DEBUG
      sciprint("DEBUG: single criteria problem loaded\n");
#endif

//       status = sym_explicit_load_problem(sci_environment, numcols, numrows, start, index, value, in_lower, in_upper, vartype, c1,
//                                          NULL, ctype, in_rhs, in_lhs, TRUE); // YC can be replaced by FALSE ?
      status = sym_explicit_load_problem(sci_environment, numcols, numrows, start, index, value, in_lower, in_upper, vartype, c1, 
	                                 NULL, ctype, in_rhs, in_lhs, TRUE); // YC can be replaced by FALSE ?
    }
  else
    {
#ifdef DEBUG
      sciprint("DEBUG: multi criteria problem loaded\n");
#endif

//       status = sym_explicit_load_problem(sci_environment, numcols, numrows, start, index, value, in_lower, in_upper, 
//                                          vartype, c1, c2, ctype, in_rhs), in_lhs, FALSE); // YC can be replaced by FALSE ?
      status = sym_explicit_load_problem(sci_environment, numcols, numrows, start, index, value, in_lower, in_upper,
					 vartype, c1, c2, ctype, in_rhs, in_lhs, FALSE); // YC can be replaced by FALSE ?
    }

#ifdef DEBUG
  sciprint("sym_explicit_load_problem: status = %d\n", status);
#endif

  // status return value:
  // - ERROR__USER
  // - FUNCTION_TERMINATED_ABNORMALLY
  // - FUNCTION_TERMINATED_NORMALLY


  ////////////////////////////////
  // Write problem if requested //
  ////////////////////////////////

  getStringInPList(pvApiCtx, param_in_addr, "writemps", &tmp_char, &tmp_res, "test.mps", Log, CHECK_NONE);
  if (tmp_res!=-1) sym_write_mps(sci_environment, tmp_char);
  FREE(tmp_char);

  ////////////////
  // Resolution //
  ////////////////

#ifdef DEBUG
  sciprint("DEBUG: resolution\n");
  sciprint("model number columns = %d\n", numcols);
  sciprint("model number rows    = %d\n", numrows);
#endif

  // YC: ajouter la possibilite d'un warm solve
  // status = sym_warm_solve(sci_environment);

  if ((n_c2==0)&&(m_c2==0))
    {
#ifdef DEBUG
      sciprint("DEBUG: single criteria resolution\n");
#endif
      // Single criteria resolution
      status = sym_solve(sci_environment);
      
      // status return value:
      // ERROR__USER
      // TM_OPTIMAL_SOLUTION_FOUND
      // TM_TIME_LIMIT_EXCEEDED
      // TM_NODE_LIMIT_EXCEEDED
      // TM_TARGET_GAP_ACHIEVED
      // TM_FOUND_FIRST_FEASIBLE
      // TM_ERROR__NO_BRANCHING_CANDIDATE
      // TM_ERROR__ILLEGAL_RETURN_CODE
      // TM_ERROR__NUMERICAL_INSTABILITY
      // TM_ERROR__COMM_ERROR
      // TM_ERROR__USER
    }
  else
    {
#ifdef DEBUG
      sciprint("DEBUG: multi criteria resolution\n");
#endif
      // Multi criteria resolution
      status = sym_mc_solve(sci_environment);
      
      // status return value:
      // ERROR__USER
      // TM_OPTIMAL_SOLUTION_FOUND
      // TM_ERROR__NO_BRANCHING_CANDIDATE
      // TM_ERROR__ILLEGAL_RETURN_CODE
      // TM_ERROR__NUMERICAL_INSTABILITY
      // TM_ERROR__COMM_ERROR
      // TM_ERROR__USER
      // FUNCTION_TERMINATED_ABNORMALLY
    }

#ifdef DEBUG
  sciprint("\nOptimization status = %d \n", status);
#endif

  int sym_status = 0;

  sym_status = sym_get_status(sci_environment);
  
  // sym_status return value:
  // ERROR__USER
  // TM_OPTIMAL_SOLUTION_FOUND
  // TM_TIME_LIMIT_EXCEEDED
  // TM_NODE_LIMIT_EXCEEDED
  // TM_TARGET_GAP_ACHIEVED
  // TM_FOUND_FIRST_FEASIBLE
  // TM_ERROR__NO_BRANCHING_CANDIDATE
  // TM_ERROR__ILLEGAL_RETURN_CODE
  // TM_ERROR__NUMERICAL_INSTABILITY
  // TM_ERROR__COMM_ERROR
  // TM_ERROR__USER

  int sym_status_2 = 0;

  sym_status_2  = (int)(pow(2.0,0.0)*sym_is_proven_optimal(sci_environment));
  sym_status_2 += (int)(pow(2.0,1.0)*sym_is_proven_primal_infeasible(sci_environment));
  sym_status_2 += (int)(pow(2.0,2.0)*sym_is_target_gap_achieved(sci_environment));
  sym_status_2 += (int)(pow(2.0,3.0)*sym_is_abandoned(sci_environment));

  //////////////////////////////
  // Allocate for return data //
  //////////////////////////////

#ifdef DEBUG
  sciprint("DEBUG: allocating data\n");
#endif

  double fmin = 0;
  int iterations = 0;

#ifdef DEBUG
  sciprint("DEBUG: getting solution\n");
#endif

  // YC: si pas de solution alors ces fonctions ne doivent pas etre appelees !!

  if ((status==TM_OPTIMAL_SOLUTION_FOUND) ||
      (status==TM_TIME_LIMIT_EXCEEDED) ||
      (status==TM_NODE_LIMIT_EXCEEDED) ||
      (status==TM_TARGET_GAP_ACHIEVED) ||
      (status==TM_FOUND_FIRST_FEASIBLE))
    {  
      primal = (double *)MALLOC(numcols*sizeof(double));
      dual   = (double *)MALLOC(numrows*sizeof(double));
      sym_get_col_solution(sci_environment,primal);
      sym_get_row_activity(sci_environment,dual);
    }

  sym_get_iteration_count(sci_environment,&iterations);

  /////////////////////////////////
  // Copy solutions if available //
  /////////////////////////////////

#ifdef DEBUG
  sciprint("DEBUG: returning data\n");
#endif

  if ((status==TM_OPTIMAL_SOLUTION_FOUND) ||
      (status==TM_TIME_LIMIT_EXCEEDED) ||
      (status==TM_NODE_LIMIT_EXCEEDED) ||
      (status==TM_TARGET_GAP_ACHIEVED) ||
      (status==TM_FOUND_FIRST_FEASIBLE))
    {  
      sym_get_obj_val(sci_environment,&fmin);
    }
  else
    {
      // An error has occured
      primal = (double *)MALLOC(numcols*sizeof(double));
      dual   = (double *)MALLOC(numrows*sizeof(double));
      memset(primal, 0, numcols*sizeof(double));
      memset(dual,   0, numrows*sizeof(double));
      fmin = 0.0;
    }

  int * extra_addr = NULL;
  double tmp_dbl[1];

  char * ListLabels [] = {"lambda","sym_status", "iterations"};

  _SciErr = createMatrixOfDouble(pvApiCtx, XMIN_OUT, numcols, 1, primal); SCICOINOR_ERROR;
  createScalarDouble(pvApiCtx, FMIN_OUT, fmin);
  createScalarDouble(pvApiCtx, STATUS_OUT, (double)sym_status);

  _SciErr = createPList(pvApiCtx, EXTRA_OUT, &extra_addr, (char **)ListLabels, 3); SCICOINOR_ERROR;

  _SciErr = createColVectorOfDoubleInPList(pvApiCtx, EXTRA_OUT, extra_addr, "lambda", numrows, dual); SCICOINOR_ERROR;
  _SciErr = createIntInPList(pvApiCtx,               EXTRA_OUT, extra_addr, "sym_status", sym_status_2); SCICOINOR_ERROR;
  _SciErr = createIntInPList(pvApiCtx,               EXTRA_OUT, extra_addr, "iterations", iterations); SCICOINOR_ERROR;

  LhsVar(1) = XMIN_OUT;
  LhsVar(2) = FMIN_OUT;
  LhsVar(3) = STATUS_OUT;
  LhsVar(4) = EXTRA_OUT;

  //////////////////////////////
  // Delete allocated objects //
  //////////////////////////////

  FREE(primal);
  FREE(dual);
  FREE(start);
  FREE(index);
  FREE(value);
  FREE(vartype);

  FREE(in_lower);
  FREE(in_upper);
  FREE(in_rhs);
  FREE(in_lhs);

  freeAllocatedSingleString(ctype);
  freeAllocatedSingleString(vtype);

  status = sym_close_environment(sci_environment);

  return 0;
}
