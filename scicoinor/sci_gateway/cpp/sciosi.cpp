//////////////////////////////////////////////////////////////////////////
// sciosi: A scilab interface to the OSI library for linear programming //
//////////////////////////////////////////////////////////////////////////

//  Copyright (C) 2010 Yann Collette.
//
//  SCIOSI is free software; you can redistribute it and/or modify it
//  under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 2, or (at your option) any
//  later version.
//
//  This part of code is distributed with the FURTHER condition that it 
//  can be compiled and linked with the Matlab libraries and it can be 
//  used within the Matlab environment.
//
//  SCIOSI is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with SciCoinOR; see the file COPYING.  If not, write to the Free
//  Software Foundation, 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#define HAS_CBC_OSI  1
//#define HAS_CPX_OSI  1
//#define HAS_DYLP_OSI 1
//#define HAS_FMP_OSI  1
//#define HAS_GLPK_OSI 1
//#define HAS_GRB_OSI  1
//#define HAS_MSK_OSI  1
//#define HAS_OSL_OSI  1
//#define HAS_SPX_OSI  1
//#define HAS_SYM_OSI  1
//#define HAS_VOL_OSI  1
//#define HAS_XPR_OSI  1

#include <cstdio>
#include <cstring>
#include <exception>

#undef min
#undef max

#include <Coin_C_defines.h>
#include <CoinPackedMatrix.hpp>
#include <CoinMessageHandler.hpp>
#include <CoinError.hpp>

#include <OsiSolverInterface.hpp>
#include <OsiClpSolverInterface.hpp>

#ifdef HAS_CBC_OSI
#include <OsiCbcSolverInterface.hpp>
#endif
#ifdef HAS_CPX_OSI
#include <OsiCpxSolverInterface.hpp>
#endif
#ifdef HAS_DYLP_OSI
#include <OsiDylpSolverInterface.hpp>
#endif
#ifdef HAS_FMP_OSI
#include <OsiFmpSolverInterface.hpp>
#endif
#ifdef HAS_GLPK_OSI
#include <OsiGlpkSolverInterface.hpp>
#endif
#ifdef HAS_GRB_OSI
#include <OsiGrbSolverInterface.hpp>
#endif
#ifdef HAS_MSK_OSI
#include <OsiMskSolverInterface.hpp>
#endif
#ifdef HAS_OSL_OSI
#include <OsiOslSolverInterface.hpp>
#endif
#ifdef HAS_SPX_OSI
#include <OsiSpxSolverInterface.hpp>
#endif
#ifdef HAS_SYM_OSI
#include <OsiSymSolverInterface.hpp>
#endif
#ifdef HAS_VOL_OSI
#include <OsiVolSolverInterface.hpp>
#endif
#ifdef HAS_XPR_OSI
#include <OsiXprSolverInterface.hpp>
#endif

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
#define C_IN        2
#define LHS_IN      3
#define RHS_IN      4
#define UPPER_IN    5
#define LOWER_IN    6
#define CTYPE_IN    7
#define VTYPE_IN    8
#define PARAM_IN    9
#define LASTPARAM   9

// Output Arguments
#define  XMIN_OUT     Rhs+1
#define  FMIN_OUT     Rhs+2
#define  STATUS_OUT   Rhs+3
#define  EXTRA_OUT    Rhs+4

#define TRYCATCH(FUNCTION) FUNCTION;

// // In the catch, we create 4 empty return variables
// #define TRYCATCH(FUNCTION)						\
//   try									\
//     {									\
//       FUNCTION;								\
//     }									\
//   catch(CoinError & e)							\
//     {									\
//       sciprint("sciosi: error %s line %d %s::%s - %s\n", e.fileName().c_str(), \
// 	       e.lineNumber(),						\
// 	       e.className().c_str(),					\
// 	       e.methodName().c_str(),					\
// 	       e.message());						\
//       									\
//       if (printer)       delete printer;				\
//       freeAllocatedSingleString(ctype);					\
//       freeAllocatedSingleString(vtype);					\
//       if (writemps_filename) FREE(writemps_filename);			\
// 									\
//       createScalarDouble(pvApiCtx, Rhs+1, 0);				\
//       createScalarDouble(pvApiCtx, Rhs+2, 0);				\
//       createScalarDouble(pvApiCtx, Rhs+3, 0);				\
//       createScalarDouble(pvApiCtx, Rhs+4, 0);				\
//       									\
//       LhsVar(1) = Rhs+1;						\
//       LhsVar(2) = Rhs+2;						\
//       LhsVar(3) = Rhs+3;						\
//       LhsVar(4) = Rhs+4;						\
//       									\
//       return 0;								\
//     }									\
//   catch(std::bad_alloc & e)						\
//     {									\
//       sciprint("sciosi: exception raised - %s\n", e.what());		\
//       									\
//       if (printer)       delete printer;				\
//       freeAllocatedSingleString(ctype);					\
//       freeAllocatedSingleString(vtype);					\
//       if (writemps_filename) FREE(writemps_filename);			\
//       									\
//       createScalarDouble(pvApiCtx, Rhs+1, 0);				\
//       createScalarDouble(pvApiCtx, Rhs+2, 0);				\
//       createScalarDouble(pvApiCtx, Rhs+3, 0);				\
//       createScalarDouble(pvApiCtx, Rhs+4, 0);				\
//       									\
//       LhsVar(1) = Rhs+1;						\
//       LhsVar(2) = Rhs+2;						\
//       LhsVar(3) = Rhs+3;						\
//       LhsVar(4) = Rhs+4;						\
//       									\
//       return 0;								\
//     }									\
//   catch(std::bad_exception & e)						\
//     {									\
//       sciprint("sciosi: exception raise - %s\n", e.what());		\
//       									\
//       if (printer)       delete printer;				\
//       freeAllocatedSingleString(ctype);					\
//       freeAllocatedSingleString(vtype);					\
//       if (writemps_filename) FREE(writemps_filename);			\
//       									\
//       createScalarDouble(pvApiCtx, Rhs+1, 0);				\
//       createScalarDouble(pvApiCtx, Rhs+2, 0);				\
//       createScalarDouble(pvApiCtx, Rhs+3, 0);				\
//       createScalarDouble(pvApiCtx, Rhs+4, 0);				\
//       									\
//       LhsVar(1) = Rhs+1;						\
//       LhsVar(2) = Rhs+2;						\
//       LhsVar(3) = Rhs+3;						\
//       LhsVar(4) = Rhs+4;						\
//       									\
//       return 0;								\
//     }									\
//   catch(...)								\
//     {									\
//       sciprint("sciosi: unknow exception raised\n");			\
//       throw;								\
//     }

// To define the function sciosi as a C function and allow this function to be loaded easily in scilab
extern "C" int sciosi(char * fname)
{
  OsiSolverInterface * modelOSI = NULL;
  CoinMessageHandler * printer  = NULL;
  CoinPackedMatrix  A_matrix;
  int Log = 0;
  int ncols = 0, nrows = 0, i, j, nz = 0, writemps = 0;
  int count = 0, status = 0, type;
  char  * writemps_filename = NULL;
  SciSparse S_A;
  
  if (Rhs<LASTPARAM) 
    {
      Scierror(999,"%s: 12 inputs required in call to %s. Bug in osi.sci ?...\n", fname, fname);
      return 0;
    }
		
  /* Get pointers to input */

  int n_a,     m_a,     * a_addr     = NULL;
  int n_c,     m_c,     * c_addr     = NULL;
  int n_lhs,   m_lhs,   * lhs_addr   = NULL;
  int n_rhs,   m_rhs,   * rhs_addr   = NULL;
  int n_upper, m_upper, * upper_addr = NULL;
  int n_lower, m_lower, * lower_addr = NULL;
  int * vtype_addr = NULL;
  int * ctype_addr = NULL;
  double * c = NULL, * lhs = NULL, * rhs = NULL, * lower = NULL, * upper = NULL, * a = NULL;
  char * ctype = NULL, * vtype = NULL;
  SciErr _SciErr;
  
  _SciErr = getVarAddressFromPosition(pvApiCtx, C_IN, &c_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, c_addr, &n_c, &m_c, &c); SCICOINOR_ERROR;
  _SciErr = getVarAddressFromPosition(pvApiCtx, LHS_IN, &lhs_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, lhs_addr, &n_lhs, &m_lhs, &lhs); SCICOINOR_ERROR;
  _SciErr = getVarAddressFromPosition(pvApiCtx, RHS_IN, &rhs_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, rhs_addr, &n_rhs, &m_rhs, &rhs); SCICOINOR_ERROR;
  _SciErr = getVarAddressFromPosition(pvApiCtx, LOWER_IN, &lower_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, lower_addr, &n_lower, &m_lower, &lower); SCICOINOR_ERROR;
  _SciErr = getVarAddressFromPosition(pvApiCtx, UPPER_IN, &upper_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, upper_addr, &n_upper, &m_upper, &upper); SCICOINOR_ERROR;
  _SciErr = getVarAddressFromPosition(pvApiCtx, CTYPE_IN, &ctype_addr); SCICOINOR_ERROR;
  getAllocatedSingleString(pvApiCtx, ctype_addr, &ctype);
  _SciErr = getVarAddressFromPosition(pvApiCtx, VTYPE_IN, &vtype_addr); SCICOINOR_ERROR;
  getAllocatedSingleString(pvApiCtx, vtype_addr, &vtype);

  ncols = n_c   * m_c;   /* Length of c == number of columns */
  nrows = n_rhs * m_rhs; /* length of b == number of rows    */
  if (nrows==0) nrows = n_lhs * m_lhs;

#ifdef DEBUG
  sciprint("DEBUG: c           size: m = %d n = %d\n", m_c,     n_c);
  sciprint("DEBUG: lhs         size: m = %d n = %d\n", m_lhs,   n_lhs);
  sciprint("DEBUG: rhs         size: m = %d n = %d\n", m_rhs,   n_rhs);
  sciprint("DEBUG: lower       size: m = %d n = %d\n", m_lower, n_lower);
  sciprint("DEBUG: upper       size: m = %d n = %d\n", m_upper, n_upper);
  sciprint("DEBUG: ctype       size: m = %d n = %d\n", m_ctype, n_ctype);
  sciprint("DEBUG: vtype       size: m = %d n = %d\n", m_vtype, n_vtype);
  sciprint("DEBUG: nrows = %d\n", nrows);
  sciprint("DEBUG: ncols = %d\n", ncols);

  sciprint("c :");
  for(i=0;i<ncols; i++) sciprint("%f ",*(c+i));
  sciprint("\n");
  sciprint("lhs :");
  for(i=0;i<nrows; i++) sciprint("%f ",*(lhs+i));
  sciprint("\n");
  sciprint("rhs :");
  for(i=0;i<nrows; i++) sciprint("%f ",*(rhs+i));
  sciprint("\n");
  sciprint("lb :");
  for(i=0;i<ncols; i++) sciprint("%f ",*(lower+i));
  sciprint("\n");
  sciprint("ub :");
  for(i=0;i<ncols; i++) sciprint("%f ",*(upper+i));
  sciprint("\n");
  sciprint("ctype = %s\n", ctype);
  sciprint("vtype = %s\n", vtype);
#endif

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
      A_matrix.setDimensions(nrows,ncols);
      
      if (a==NULL) 
	{
	  Scierror(999,"%s: invalid value of matrix a\n",fname);
	  
	  freeAllocatedSingleString(ctype);
	  freeAllocatedSingleString(vtype);
	  
	  return 0;
	}
      
      for(i=0; i<m_a; i++)
	{
	  for(j=0; j<n_a; j++)
	    {
	      if (*(a+i+j*m_a) != 0) A_matrix.modifyCoefficient(i,j,*(a+i+j*m_a));
#ifdef DEBUG
	      sciprint("%f ",*(a+i+j*m_a));
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
      A_matrix.setDimensions(nrows,ncols);
      
#ifdef DEBUG
      sciprint("A = [%d,%d]\n", m_a, n_a);
#endif
      
      nz = S_A.nel;
      
      count = 0;
      for(i=0;i<S_A.m;i++)
	{
	  if (S_A.mnel[i]!=0) 
	    {
#ifdef DEBUG
	      sciprint("mnel[%d] = %d - ",i, S_A.mnel[i]);
#endif
	      for(j=0;j<S_A.mnel[i];j++)
		{
		  count++;
		  A_matrix.modifyCoefficient(i,S_A.icol[count-1]-1,S_A.R[count-1]);
#ifdef DEBUG
		  sciprint("[%d] = %f ", S_A.icol[count-1]-1, S_A.R[count-1]);
#endif
		}
#ifdef DEBUG
	      sciprint("\n");
#endif
	    }
	}
      
      freeAllocatedSparseMatrix(S_A.mnel, S_A.icol, S_A.R);
    }
  
  /////////////////
  // Get options //
  /////////////////

#ifdef DEBUG
  sciprint("DEBUG: get options\n");
#endif

  // Default settings
  int * param_in_addr = NULL;
  int     loglevel = 0;
  int     tmp_int, tmp_res;
  double  tmp_double;
  char *  tmp_char;

  initPList(pvApiCtx, PARAM_IN, &param_in_addr);
  if (!checkPList(pvApiCtx, param_in_addr))
    {
      Scierror(999, "%s: argument n° %d is not a plist\n", fname, PARAM_IN);

      return 0;
    }

  // writemps option
  getStringInPList(pvApiCtx, param_in_addr, "solvername", &tmp_char, &tmp_res, "clp", Log, CHECK_NONE);
  if (strcmp(tmp_char,"clp")==0)
    {
      modelOSI = new OsiClpSolverInterface();
    }
#ifdef HAS_CBC_OSI
  else if (strcmp(tmp_char,"cbc")==0)
    {
      modelOSI = new OsiCbcSolverInterface();
    }
#endif
#ifdef HAS_CPX_OSI
  else if (strcmp(tmp_char,"cpx")==0)
    {
      modelOSI = new OsiCpxSolverInterface();
    }
#endif
#ifdef HAS_DYLP_OSI
  else if (strcmp(tmp_char,"dylp")==0)
    {
      modelOSI = new OsiDylpSolverInterface();
    }
#endif
#ifdef HAS_FMP_OSI
  else if (strcmp(tmp_char,"fmp")==0)
    {
      modelOSI = new OsiFmpSolverInterface();
    }
#endif
#ifdef HAS_GLPK_OSI
  else if (strcmp(tmp_char,"glpk")==0)
    {
      modelOSI = new OsiGlpkSolverInterface();
    }
#endif
#ifdef HAS_GRB_OSI
  else if (strcmp(tmp_char,"grb")==0)
    {
      modelOSI = new OsiGrbSolverInterface();
    }
#endif
#ifdef HAS_MSK_OSI
  else if (strcmp(tmp_char,"msk")==0)
    {
      modelOSI = new OsiMskSolverInterface();
    }
#endif
#ifdef HAS_OSL_OSI
  else if (strcmp(tmp_char,"osl")==0)
    {
      modelOSI = new OsiOslSolverInterface();
    }
#endif
#ifdef HAS_SPX_OSI
  else if (strcmp(tmp_char,"spx")==0)
    {
      modelOSI = new OsiSpxSolverInterface();
    }
#endif
#ifdef HAS_SYM_OSI
  else if (strcmp(tmp_char,"sym")==0)
    {
      modelOSI = new OsiSymSolverInterface();
    }
#endif
#ifdef HAS_VOL_OSI
  else if (strcmp(tmp_char,"vol")==0)
    {
      modelOSI = new OsiVolSolverInterface();
    }
#endif
#ifdef HAS_XPR_OSI
  else if (strcmp(tmp_char,"xpr")==0)
    {
      modelOSI = new OsiXprSolverInterface();
    }
#endif
  else
    {
      Scierror(999,"%s: wrong solver name.\n", fname);
      return 0;
    }

  // verbose option
  getIntInPList(pvApiCtx, param_in_addr, "verbose", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) 
    {
      ///////////////////////////////
      // Enable printing in Scilab //
      // and set parameters of osi //
      ///////////////////////////////
	  
      loglevel = tmp_int;

      printer = new CoinMessageHandler(); // assumed open	
      printer->setLogLevel(loglevel);
      modelOSI->passInMessageHandler(printer);
    }

  // writemps option
  getStringInPList(pvApiCtx, param_in_addr, "writemps", &tmp_char, &tmp_res, "test.mps", Log, CHECK_NONE);
  if (tmp_res!=-1) 
    {
      writemps_filename = tmp_char;
      writemps = 1;
#ifdef DEBUG
      sciprint("DEBUG: writemps_filename = %s\n", writemps_filename);
#endif
    }

  // optimization direction
  // 1: minimization, -1 maximization
  getDoubleInPList(pvApiCtx, param_in_addr, "optim_dir", &tmp_double, &tmp_res, 1, Log, CHECK_VALUES, 2, 1, -1);
  if (tmp_res!=-1) modelOSI->setObjSense(tmp_double);

  // maxnumiteration
  getIntInPList(pvApiCtx, param_in_addr, "maxnumiteration", &tmp_int, &tmp_res, 1000, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) modelOSI->setIntParam(OsiMaxNumIteration, tmp_int);

  // maxnumiterationshotstart
  getIntInPList(pvApiCtx, param_in_addr, "maxnumiterationshotstart", &tmp_int, &tmp_res, 1000, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) modelOSI->setIntParam(OsiMaxNumIterationHotStart, tmp_int);

  // name discipline
  getIntInPList(pvApiCtx, param_in_addr, "namediscipline", &tmp_int, &tmp_res, 0, Log, CHECK_VALUES, 3, 0, 1, 2);
  if (tmp_res!=-1) modelOSI->setIntParam(OsiNameDiscipline, tmp_int);

  // dualobjectivelimit
  getDoubleInPList(pvApiCtx, param_in_addr, "dualobjectivelimit", &tmp_double, &tmp_res, 1e-7, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) modelOSI->setDblParam(OsiDualObjectiveLimit, tmp_double);

  // primalobjectivelimit
  getDoubleInPList(pvApiCtx, param_in_addr, "primalobjectivelimit", &tmp_double, &tmp_res, 1e-7, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) modelOSI->setDblParam(OsiPrimalObjectiveLimit, tmp_double);

  // obj offset
  getDoubleInPList(pvApiCtx, param_in_addr, "objoffset", &tmp_double, &tmp_res, 0.0, Log, CHECK_NONE);
  if (tmp_res!=-1) modelOSI->setDblParam(OsiObjOffset, tmp_double);

  // dual tolerance
  getDoubleInPList(pvApiCtx, param_in_addr, "dualtolerance", &tmp_double, &tmp_res, 1e-7, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) modelOSI->setDblParam(OsiDualTolerance, tmp_double);
      
  // primal tolerance
  getDoubleInPList(pvApiCtx, param_in_addr, "primaltolerance", &tmp_double, &tmp_res, 1e-7, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) modelOSI->setDblParam(OsiPrimalTolerance, tmp_double);

  // probname option
  getStringInPList(pvApiCtx, param_in_addr, "probname", &tmp_char, &tmp_res, "test", Log, CHECK_NONE);
  if (tmp_res!=-1) modelOSI->setStrParam(OsiProbName, tmp_char);

  // OsiDoPresolveInInitial     Whether to do a presolve in initialSolve.
  getIntInPList(pvApiCtx, param_in_addr, "dopresolveininitial", &tmp_int, &tmp_res, 0, Log, CHECK_VALUES, 2, 0, 1);
  if (tmp_res!=-1) modelOSI->setHintParam(OsiDoPresolveInInitial, tmp_int);

  // OsiDoDualInInitial         Whether to use a dual algorithm in initialSolve.
  getIntInPList(pvApiCtx, param_in_addr, "dodualininitial", &tmp_int, &tmp_res, 0, Log, CHECK_VALUES, 2, 0, 1);
  if (tmp_res!=-1) modelOSI->setHintParam(OsiDoDualInInitial, tmp_int);

  // The reverse is to use a primal algorithm
  // OsiDoPresolveInResolve     Whether to do a presolve in resolve.
  getIntInPList(pvApiCtx, param_in_addr, "dopresolveinresolve", &tmp_int, &tmp_res, 0, Log, CHECK_VALUES, 2, 0, 1);
  if (tmp_res!=-1) modelOSI->setHintParam(OsiDoPresolveInResolve, tmp_int);

  // OsiDoDualInResolve         Whether to use a dual algorithm in resolve.
  getIntInPList(pvApiCtx, param_in_addr, "dodualinresolve", &tmp_int, &tmp_res, 0, Log, CHECK_VALUES, 2, 0, 1);
  if (tmp_res!=-1) modelOSI->setHintParam(OsiDoDualInResolve, tmp_int);

  // The reverse is to use a primal algorithm
  // OsiDoScale                 Whether to scale problem.
  getIntInPList(pvApiCtx, param_in_addr, "doscale", &tmp_int, &tmp_res, 0, Log, CHECK_VALUES, 2, 0, 1);
  if (tmp_res!=-1) modelOSI->setHintParam(OsiDoScale, tmp_int);

  // OsiDoCrash                 Whether to create a non-slack basis (only in initialSolve).
  getIntInPList(pvApiCtx, param_in_addr, "docrash", &tmp_int, &tmp_res, 0, Log, CHECK_VALUES, 2, 0, 1);
  if (tmp_res!=-1) modelOSI->setHintParam(OsiDoCrash, tmp_int);

  // OsiDoReducePrint           Whether to reduce amount of printout, e.g., for branch and cut.
  getIntInPList(pvApiCtx, param_in_addr, "doreduceprint", &tmp_int, &tmp_res, 0, Log, CHECK_VALUES, 2, 0, 1);
  if (tmp_res!=-1) modelOSI->setHintParam(OsiDoReducePrint, tmp_int);

  // OsiDoInBranchAndCut        Whether we are in branch and cut - so can modify behavior. 
  getIntInPList(pvApiCtx, param_in_addr, "doinbranchandcut", &tmp_int, &tmp_res, 0, Log, CHECK_VALUES, 2, 0, 1);
  if (tmp_res!=-1) modelOSI->setHintParam(OsiDoInBranchAndCut, tmp_int);

#ifdef HAS_OSI_CBC
  //////////////////////////
  // CBC Specific options //
  //////////////////////////

  getDoubleInPList(pvApiCtx, param_in_addr, "cutoff", &tmp_double, &tmp_res, 1e50, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) modelOSI->setCutOff(tmp_double);

  getIntInPList(pvApiCtx, param_in_addr, "maximumnodes", &tmp_int, &tmp_res, 50000, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) modelOSI->setMaximumNodes(tmp_int);

  getIntInPList(pvApiCtx, param_in_addr, "maximumsolutions", &tmp_int, &tmp_res, 1, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) modelOSI->setMaximumSolutions(tmp_int);

  getDoubleInPList(pvApiCtx, param_in_addr, "maximumseconds", &tmp_double, &tmp_res, 1e100, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) modelOSI->setMaximumSeconds(tmp_double);
#endif

#ifdef HAS_OSI_CLP
  // void setSpecialOptions (unsigned int value)
  // Get pointer to Clp model.
  // int lastAlgorithm () const
  // Last algorithm used , 1 = primal, 2 = dual.
  // int cleanupScaling () const
  // Get scaling action option.
  // void setCleanupScaling (int value)
  // Set Scaling option When scaling is on it is possible that the scaled problem is feasible but the unscaled is not.
  // double smallestElementInCut () const
  // Get smallest allowed element in cut.
  // void setSmallestElementInCut (double value)
  // Set smallest allowed element in cut.
  // double smallestChangeInCut () const
  // Get smallest change in cut.
  // void setSmallestChangeInCut (double value)
  // Set smallest change in cut. 
#endif

#ifdef HAS_OSI_CPX
  // Access to a pointer to a cplex instance - need cplex documentation
#endif

#ifdef HAS_OSI_DYLP
  // void dylp_controlfile (const char *name, const bool silent, const bool mustexist=true)
  // Process an options (.spc) file.
  // void dylp_logfile (const char *name, bool echo=false)
  // Establish a log file.
  // void dylp_outfile (const char *name)
  // Establish an output (solution and/or statistics) file.
  // void dylp_printsoln (bool wantSoln, bool wantStats)
  // Print the solution and/or statistics to the output file. #endif
#endif

#ifdef HAS_OSI_FMP
  // No supplementary options available
#endif

#ifdef HAS_OSI_GLPK
  // No supplementary options available
#endif

#ifdef HAS_OSI_GRB
  // Access to a pointer to a gurobi instance - need a gurobi documentation
#endif

#ifdef HAS_OSI_MSK
  // Access to a pointer to a mosek instance - need a mosek documentation
#endif

#ifdef HAS_OSI_OSL
  // Access to a pointer to a OSL instance - need a OSL documentation
#endif

#ifdef HAS_OSI_SPX
  // No supplementary options available
#endif

#ifdef HAS_OSI_SYM
  // No supplementary options available
#endif

#ifdef HAS_OSI_VOL
  // Access to a pointer to a VOL instance
#endif

#ifdef HAS_OSI_XPR
  // No supplementary options available
#endif

  /////////////////////////////////////
  // Set the bounds on the variables //
  /////////////////////////////////////

  modelOSI->loadProblem(A_matrix,lhs,rhs,c,lower,upper);

  for(i=0;i<ncols; i++)
    {
      modelOSI->setColUpper(i, *(upper+i));
      modelOSI->setColLower(i, *(lower+i));
      modelOSI->setObjCoeff(i,*(c+i));
      if ((*(vtype+i)=='I')||(*(vtype+i)=='i'))
	{
	  modelOSI->setInteger(i);
	}
      else
	{
	  modelOSI->setContinuous(i);
	}
    }
  
  /////////////////////////////////////////
  // Set the boundary of the constraints //
  /////////////////////////////////////////

  // 'L' - smaller than - <=
  // 'E' - equality     - =
  // 'G' - greater than - >=
  // 'R' - Range        - <= + >=
  // 'N' - Free         - no constraints

#ifdef DEBUG
  sciprint("DEBUG: dealing with btype\n");
#endif
  for(i=0;i<nrows; i++)
    {
      switch(*(ctype+i))
	{
	case 'l':
	case 'L':
	  modelOSI->setRowUpper(i, *(rhs+i));
	  modelOSI->setRowLower(i, -COIN_DBL_MAX);
	  break;
	case 'e':
	case 'E':
	  modelOSI->setRowUpper(i, *(rhs+i));
	  modelOSI->setRowLower(i, *(rhs+i));
	  break;
	case 'n':
	case 'N':
	  modelOSI->setRowUpper(i,  COIN_DBL_MAX);
	  modelOSI->setRowLower(i, -COIN_DBL_MAX);
	  break;
	case 'r':
	case 'R':
	  modelOSI->setRowUpper(i, *(rhs+i));
	  modelOSI->setRowLower(i, *(lhs+i));
	  break;
	case 'g':
	case 'G':
	default:
	  modelOSI->setRowUpper(i, COIN_DBL_MAX);
	  modelOSI->setRowLower(i, *(lhs+i));
	  break;
	}
#ifdef DEBUG
      sciprint("row lower[%d] = %f row upper[%d] = %f\n", i, modelOSI->getRowLower()[i], i, modelOSI->getRowUpper()[i]);
#endif
    }
  
  ////////////////////////////////////////////
  // If needed, write the problem in a file //
  ////////////////////////////////////////////

#ifdef DEBUG
  sciprint("DEBUG: dealing with writemps\n");
#endif
  if (writemps)
    {
      modelOSI->writeMps(writemps_filename);
      if (loglevel) sciprint("sciosi: writing %s mps file\n",writemps_filename);
    }

  ////////////////
  // Resolution //
  ////////////////

#ifdef DEBUG
  sciprint("DEBUG: resolution\n");
  sciprint("model number columns = %d\n", modelOSI->numberColumns());
  sciprint("model number rows    = %d\n", modelOSI->numberRows());
#endif

  TRYCATCH(modelOSI->initialSolve())

#ifdef DEBUG
  if (status && loglevel) sciprint("Optimization failed\n");
#endif

  int osi_status = 0;

  osi_status  = (int)(pow(2.0,0.0)*modelOSI->isAbandoned());
  osi_status += (int)(pow(2.0,1.0)*modelOSI->isProvenOptimal());
  osi_status += (int)(pow(2.0,2.0)*modelOSI->isProvenPrimalInfeasible());
  osi_status += (int)(pow(2.0,3.0)*modelOSI->isProvenDualInfeasible());
  osi_status += (int)(pow(2.0,4.0)*modelOSI->isPrimalObjectiveLimitReached());
  osi_status += (int)(pow(2.0,5.0)*modelOSI->isDualObjectiveLimitReached());  
  osi_status += (int)(pow(2.0,6.0)*modelOSI->isIterationLimitReached());      

  ////////////////////////////
  // Specific return values //
  ////////////////////////////

#ifdef HAS_OSI_CBC
  // bool isNodeLimitReached () const
  // Node limit reached?
  // bool isSolutionLimitReached () const
  // Solution limit reached?
  // int getNodeCount () const
  // Get how many Nodes it took to solve the problem.
  // int status () const
  // Final status of problem - 0 finished, 1 stopped, 2 difficulties. 
#endif

  //////////////////////////////
  // Allocate for return data //
  //////////////////////////////

#ifdef DEBUG
  sciprint("DEBUG: allocating data\n");
#endif

  int m_xmin   = ncols, n_xmin   = 1;
  int m_lambda   = 1, n_lambda   = nrows;
  int * extra_addr = NULL;
  double tmp_dbl[1];
  char * ListLabels [] = {"lambda"};

  _SciErr = createMatrixOfDouble(pvApiCtx, XMIN_OUT, m_xmin, n_xmin, (double *)modelOSI->getColSolution()); SCICOINOR_ERROR;
  createScalarDouble(pvApiCtx, FMIN_OUT, modelOSI->getObjValue());
  createScalarDouble(pvApiCtx, STATUS_OUT, (double)osi_status);

  _SciErr = createPList(pvApiCtx, EXTRA_OUT, &extra_addr, (char **)ListLabels, 1); SCICOINOR_ERROR;

  _SciErr = createColVectorOfDoubleInPList(pvApiCtx, EXTRA_OUT, extra_addr, "lambda", m_lambda*n_lambda, (double *)modelOSI->getRowPrice()); SCICOINOR_ERROR;

#ifdef DEBUG
  sciprint("DEBUG: getting solution\n");
#endif

  //
  // status of problem:
  //  -1 - unknown e.g. before solve or if postSolve says not optimal
  //   0 - optimal
  //   1 - primal infeasible
  //   2 - dual infeasible
  //   3 - stopped on iterations or time
  //   4 - stopped due to errors
  //   5 - stopped by event handler (virtual int ClpEventHandler::event())
  //

  // Secondary status of problem - may get extended 
  // - 0 - none 
  // - 1 - primal infeasible because dual limit reached OR probably primal infeasible but can't prove it (main status 4)
  // - 2 - scaled problem optimal - unscaled problem has primal infeasibilities
  // - 3 - scaled problem optimal - unscaled problem has dual infeasibilities
  // - 4 - scaled problem optimal - unscaled problem has primal and dual infeasibilities
  // - 5 - giving up in primal with flagged variables
  // - 6 - failed due to empty problem check
  // - 7 - postSolve says not optimal
  // - 8 - failed due to bad element check
  // - 9 - status was 3 and stopped on time 
  // - 100 up - translation of enum from ClpEventHandler.

  /////////////////////////////////
  // Copy solutions if available //
  /////////////////////////////////

#ifdef DEBUG
  sciprint("DEBUG: returning data\n");
#endif

  LhsVar(1) = XMIN_OUT;
  LhsVar(2) = FMIN_OUT;
  LhsVar(3) = STATUS_OUT;
  LhsVar(4) = EXTRA_OUT;

  //////////////////////////////
  // Delete allocated objects //
  //////////////////////////////

  if (printer)       delete printer;

  freeAllocatedSingleString(ctype);
  freeAllocatedSingleString(vtype);

  if (writemps_filename) FREE(writemps_filename);

  return 0;
}
