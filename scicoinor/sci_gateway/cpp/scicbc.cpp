//////////////////////////////////////////////////////////////////////////////////////////////////////
// scicbc: A scilab interface to the Bcp library for linear and quadratic mixed integer programming //
//////////////////////////////////////////////////////////////////////////////////////////////////////

//  Copyright (C) 2008-2010 Yann Collette.
//
//  SCICBC is free software; you can redistribute it and/or modify it
//  under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 2, or (at your option) any
//  later version.
//
//  This part of code is distributed with the FURTHER condition that it 
//  can be compiled and linked with the Matlab libraries and it can be 
//  used within the Matlab environment.
//
//  SCICLP is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with SciCoinOR; see the file COPYING.  If not, write to the Free
//  Software Foundation, 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#include <signal.h>

#include <new>
#include <exception>

#undef min
#undef max

#include <Coin_C_defines.h>
#include <CoinPackedMatrix.hpp>
#include <CoinMessageHandler.hpp>
#include <CoinError.hpp>
#include <ClpSimplex.hpp>
#include <ClpInterior.hpp>
#include <ClpCholeskyBase.hpp>

#include <CbcConfig.h>
#include <OsiSolverInterface.hpp>
#include <CbcModel.hpp>
#include <CbcEventHandler.hpp>

#include <CglAllDifferent.hpp>
#include <CglClique.hpp>
#include <CglCutGenerator.hpp>
#include <CglDuplicateRow.hpp>
#include <CglFlowCover.hpp>
#include <CglGomory.hpp>
#include <CglKnapsackCover.hpp>
#include <CglLandP.hpp>
#include <CglLiftAndProject.hpp>
#include <CglMixedIntegerRounding2.hpp>
#include <CglMixedIntegerRounding.hpp>
#include <CglOddHole.hpp>
#include <CglProbing.hpp>
#include <CglRedSplit.hpp>
#include <CglResidualCapacity.hpp>
#include <CglSimpleRounding.hpp>
#include <CglStored.hpp>
#include <CglTreeInfo.hpp>
#include <CglTwomir.hpp>
#include <CglPreProcess.hpp>

#include <CbcHeuristic.hpp>
#include <CbcHeuristicDive.hpp>
#include <CbcHeuristicDiveCoefficient.hpp>
#include <CbcHeuristicDiveFractional.hpp>
#include <CbcHeuristicDiveGuided.hpp>
#include <CbcHeuristicDiveVectorLength.hpp>
#define COIN_HAS_CLP 1
#include <CbcLinked.hpp>
#include <CbcHeuristicFPump.hpp>
#include <CbcHeuristicGreedy.hpp>
#include <CbcHeuristicLocal.hpp>
#include <CbcHeuristic.hpp>
#include <CbcHeuristicRINS.hpp>

#include <CbcCompareActual.hpp>
#include <CbcCompareEstimate.hpp>
#include <CbcCompareObjective.hpp>
#include <CbcStrategy.hpp>
#include <ClpPresolve.hpp>
#include <CbcBranchActual.hpp>
#include <CbcBranchDynamic.hpp>
#include <CbcBranchLotsize.hpp>

#include <CglPreProcess.hpp>

#include <CbcSolver.hpp>
#include <OsiSolverInterface.hpp>
#include <OsiClpSolverInterface.hpp>
#include <OsiCbcSolverInterface.hpp>

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
//#define DBGPRINTF printf
#define DBGPRINTF sciprint

#include <helper.hpp>

using namespace std;

//   [xmin,fmin,status,extra] = scicbc(A,c,lhs,rhs,ub,lb,btype,vartype,Q,options,...
//                                     special('which'), ...
//                                     special('weight'), ...
//                                     special('type'), ...
//                                     special('column'), ...
//                                     special('range'), ...
//                                     special('id_obj'), ...
//                                     special('length'), ...
//                                     special('id'), ...
//                                     special('clique_type'));

#define A_IN                   1
#define C_IN                   2
#define LHS_IN                 3
#define RHS_IN                 4
#define UPPER_IN               5
#define LOWER_IN               6
#define CTYPE_IN               7
#define VTYPE_IN               8
#define Q_IN                   9
#define PARAM_IN               10
#define SPECIAL_CONSTR_IN      11
#define LASTPARAM              11
// Output Arguments
#define  XMIN_OUT              Rhs+1
#define  FMIN_OUT              Rhs+2
#define  STATUS_OUT            Rhs+3
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


// scilab master declares the basbrk symbol which allows to handle the ctrl+c event in scilab
#define HANDLE_CTRLC 1

#define TRYCATCH(FUNCTION) FUNCTION;

// // In the catch, we create 4 empty return variables
// #define TRYCATCH(FUNCTION)						\
//   try									\
//     {									\
//       FUNCTION;								\
//     }									\
//   catch(CoinError & e)							\
//     {									\
//       sciprint("scicbc: error %s line %d %s::%s - %s\n", e.fileName().c_str(), \
// 	       e.lineNumber(),						\
// 	       e.className().c_str(),					\
// 	       e.methodName().c_str(),					\
// 	       e.message());						\
//       									\
//       if (printer)      delete printer;					\
//       if (clpprinter)   delete clpprinter;				\
//       if (A_matrix)     delete A_matrix;				\
//       if (Q_matrix)     delete Q_matrix;				\
//       									\
//       createEmptyMatrix(pvApiCtx,Rhs+1);				\
//       createEmptyMatrix(pvApiCtx,Rhs+2);				\
//       createEmptyMatrix(pvApiCtx,Rhs+3);				\
//       createEmptyMatrix(pvApiCtx,Rhs+4);				\
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
//       sciprint("scicbc: %s\n", e.what());				\
//       									\
//       if (printer)      delete printer;					\
//       if (clpprinter)   delete clpprinter;				\
//       if (A_matrix)     delete A_matrix;				\
//       if (Q_matrix)     delete Q_matrix;				\
//       									\
//       createEmptyMatrix(pvApiCtx,Rhs+1);				\
//       createEmptyMatrix(pvApiCtx,Rhs+2);				\
//       createEmptyMatrix(pvApiCtx,Rhs+3);				\
//       createEmptyMatrix(pvApiCtx,Rhs+4);				\
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
//       sciprint("scicbc: %s\n", e.what());				\
//       									\
//       if (printer)      delete printer;					\
//       if (clpprinter)   delete clpprinter;				\
//       if (A_matrix)     delete A_matrix;				\
//       if (Q_matrix)     delete Q_matrix;				\
//       									\
//       createEmptyMatrix(pvApiCtx,Rhs+1);				\
//       createEmptyMatrix(pvApiCtx,Rhs+2);				\
//       createEmptyMatrix(pvApiCtx,Rhs+3);				\
//       createEmptyMatrix(pvApiCtx,Rhs+4);				\
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
//       sciprint("scicbc: unknow exception raised\n");			\
//       throw;								\
//     }

// from bcp / sample1.cpp
double cpuTime()
{
  double cpu_temp;
#if defined(_MSC_VER)
  unsigned int ticksnow;
  ticksnow = (unsigned int)clock();
  cpu_temp = (double)((double)ticksnow/ CLOCKS_PER_SEC);
#else
  struct rusage usage;
  getrusage(RUSAGE_SELF,&usage);
  cpu_temp = usage.ru_utime.tv_sec;
  cpu_temp += 1.0e-6*((double)usage.ru_utime.tv_usec);
#endif
  return cpu_temp;
}

// If cbc has been compiled with thread support (--enable-cbc-parallel) then uncomment this #define
//#define THREAD_SUPPORT

// To define the function sciclp as a C function and allow this function to be loaded easily in scilab
extern "C" int scicbc(char * fname)
{
  ClpSimplex     * modelSimplex     = NULL;
  ClpSimplex     * modelSimplexCopy = NULL; // To keep track of the allocated modelSimplex pointer
  DerivedHandler * printer          = NULL;
  DerivedHandler * clpprinter       = NULL;
  CbcObject     ** specialObjects   = NULL;
  CoinPackedMatrix * A_matrix       = NULL;
  CoinPackedMatrix * Q_matrix       = NULL;
#ifdef HANDLE_CTRLC
  MyClpEventHandler SciClpEventHandler;
  SciClpEventHandler.setInCbc(true);
  MyCbcEventHandler SciCbcEventHandler;
#endif
  double time0, tmp_double;
  int Log = 0;
  int ncols, nrows, i, j, nz;
  int count = 0, tmp_int;
  int IsMIP = 0;
  double status = -1;

  SciSparse S_A, S_Q;

  time0 = cpuTime();

  if (Rhs<LASTPARAM) 
    {
      Scierror(999,"%s: %d inputs required in call to %s. Bug in cbc.sci ?...\n",fname, fname, LASTPARAM);
      return 0;
    }
		
  // Get pointers to input

  int m_a                   = 0, n_a                   = 0, * a_addr              = NULL;
  int m_c                   = 0, n_c                   = 0, * c_addr              = NULL;
  int m_rhs                 = 0, n_rhs                 = 0, * rhs_addr            = NULL;
  int m_lhs                 = 0, n_lhs                 = 0, * lhs_addr            = NULL;
  int m_upper               = 0, n_upper               = 0, * upper_addr          = NULL;
  int m_lower               = 0, n_lower               = 0, * lower_addr          = NULL;
  int * vtype_addr          = NULL;
  int * ctype_addr          = NULL;
  int m_q                   = 0, n_q                   = 0, * q_addr              = NULL;
  int var_type;
  double * a = NULL, * c = NULL, * lhs = NULL, * rhs = NULL, * lower = NULL, * upper = NULL, * q = NULL;
  char * ctype = NULL, * vtype = NULL;
  bool special_stored = false;
  SciErr _SciErr;

  _SciErr = getVarAddressFromPosition(pvApiCtx, C_IN, &c_addr);
  _SciErr = getMatrixOfDouble(pvApiCtx, c_addr, &n_c, &m_c, &c);

  _SciErr = getVarAddressFromPosition(pvApiCtx, LHS_IN, &lhs_addr);
  _SciErr = getMatrixOfDouble(pvApiCtx, lhs_addr, &n_lhs, &m_lhs, &lhs);

  _SciErr = getVarAddressFromPosition(pvApiCtx, RHS_IN, &rhs_addr);
  _SciErr = getMatrixOfDouble(pvApiCtx, rhs_addr, &n_rhs, &m_rhs, &rhs);

  _SciErr = getVarAddressFromPosition(pvApiCtx, UPPER_IN, &upper_addr);
  _SciErr = getMatrixOfDouble(pvApiCtx, upper_addr, &n_upper, &m_upper, &upper);

  _SciErr = getVarAddressFromPosition(pvApiCtx, LOWER_IN, &lower_addr);
  _SciErr = getMatrixOfDouble(pvApiCtx, lower_addr, &n_lower, &m_lower, &lower);

  _SciErr = getVarAddressFromPosition(pvApiCtx, CTYPE_IN, &ctype_addr);
  getAllocatedSingleString(pvApiCtx, ctype_addr, &ctype);

  _SciErr = getVarAddressFromPosition(pvApiCtx, VTYPE_IN, &vtype_addr);
  getAllocatedSingleString(pvApiCtx, vtype_addr, &vtype);

  ncols = n_c   * m_c;   // Length of c == number of columns
  nrows = n_rhs * m_rhs; // length of b == number of rows 
  if (nrows==0) nrows = n_lhs * m_lhs;

#ifdef DEBUG
  sciprint("DEBUG: c           size: m = %d n = %d\n", m_c,     n_c);
  sciprint("DEBUG: lhs         size: m = %d n = %d\n", m_lhs,   n_lhs);
  sciprint("DEBUG: rhs         size: m = %d n = %d\n", m_rhs,   n_rhs);
  sciprint("DEBUG: lower       size: m = %d n = %d\n", m_lower, n_lower);
  sciprint("DEBUG: upper       size: m = %d n = %d\n", m_upper, n_upper);
  sciprint("DEBUG: nrows = %d\n", nrows);
  sciprint("DEBUG: ncols = %d\n", ncols);
#endif

  _SciErr = getVarAddressFromPosition(pvApiCtx, A_IN, &a_addr);
  _SciErr = getVarType(pvApiCtx, a_addr, &var_type);
  if (var_type!=sci_sparse)
    {
#ifdef DEBUG
      sciprint("DEBUG: A_IN is not sparse\n");
#endif

      _SciErr = getMatrixOfDouble(pvApiCtx, a_addr, &m_a, &n_a, &a); // get the matrix
      
      A_matrix = new CoinPackedMatrix();
      A_matrix->setDimensions(nrows,ncols);
      
      if (a==NULL) 
	{
	  Scierror(999,"%s: invalid value of matrix a\n",fname);
	  return 0;
	}
      
      for(i=0; i<m_a; i++)
	{
	  for(j=0; j<n_a; j++)
	    {
	      if (*(a+i+j*m_a) != 0) A_matrix->modifyCoefficient(i,j,*(a+i+j*m_a));
	    }
	}
    }
  else
    {
#ifdef DEBUG
      sciprint("DEBUG: A_IN is sparse\n");
#endif
      
      getAllocatedSparseMatrix(pvApiCtx, a_addr, &S_A.m, &S_A.n, &S_A.nel, &S_A.mnel, &S_A.icol, &S_A.R);
      
      A_matrix = new CoinPackedMatrix();
      A_matrix->setDimensions(nrows,ncols);
      
      nz = S_A.nel;
      
      count = 0;
      for(i=0;i<S_A.m;i++)
	{
	  if (S_A.mnel[i]!=0) 
	    {
	      for(j=0;j<S_A.mnel[i];j++)
		{
		  count++;
		  A_matrix->modifyCoefficient(i,S_A.icol[count-1]-1,S_A.R[count-1]);
		}
	    }
	}

      freeAllocatedSparseMatrix(S_A.mnel, S_A.icol, S_A.R);
    }
  
  /////////////////////
  // Set the problem //
  /////////////////////

  modelSimplex = new ClpSimplex();
  modelSimplexCopy = modelSimplex;

  modelSimplex->loadProblem(*A_matrix, lhs, rhs, c, lower, upper);

  for(i=0;i<ncols; i++)
    {
      modelSimplex->setColUpper(i, upper[i]);
      modelSimplex->setColLower(i, lower[i]);
      modelSimplex->setObjectiveCoefficient(i, c[i]);

      if ((vtype[i]=='I')||(vtype[i]=='i'))
	{
	  modelSimplex->setInteger(i);
	  IsMIP = 1;
	}
      else
	{
	  modelSimplex->setContinuous(i);
	}
    }
  
  freeAllocatedSingleString(vtype);

  // Set the boundary of the constraints
  
  // type of auxiliary/structural variable:
  // 'L' - smaller than - <=
  // 'E' - equality     - =
  // 'G' - greater than - >=
  // 'R' - Range        - <= + >=
  // 'N' - Free         - no constraints

  for(i=0;i<nrows; i++)
    {
      switch(ctype[i])
	{
	case 'l':
	case 'L':
	  modelSimplex->setRowUpper(i, rhs[i]);
	  modelSimplex->setRowLower(i, -COIN_DBL_MAX);
	  break;
	case 'e':
	case 'E':
	  modelSimplex->setRowUpper(i, rhs[i]);
	  modelSimplex->setRowLower(i, rhs[i]);
	  break;
	case 'n':
	case 'N':
	  modelSimplex->setRowUpper(i,  COIN_DBL_MAX);
	  modelSimplex->setRowLower(i, -COIN_DBL_MAX);
	  break;
	case 'r':
	case 'R':
	  modelSimplex->setRowUpper(i, rhs[i]);
	  modelSimplex->setRowLower(i, lhs[i]);
	case 'g':
	case 'G':
	default:
	  modelSimplex->setRowUpper(i, COIN_DBL_MAX);
	  modelSimplex->setRowLower(i, lhs[i]);
	  break;
	}
    }

  freeAllocatedSingleString(ctype);

  ///////////////////////////////////
  // Affect names to rows and cols //
  ///////////////////////////////////

  modelSimplex->setIntParam(ClpNameDiscipline,0);

#ifdef HANDLE_CTRLC
  ///////////////////
  // Handle Ctrl+C //
  ///////////////////

  modelSimplex->passInEventHandler(&SciClpEventHandler);
#endif

  ////////////////////////
  // Any quadratic part //
  ////////////////////////
  
  _SciErr = getVarAddressFromPosition(pvApiCtx, Q_IN, &q_addr);
  _SciErr = getVarType(pvApiCtx, q_addr, &var_type);

  if(var_type!=sci_sparse)
    {
      _SciErr = getMatrixOfDouble(pvApiCtx, q_addr, &m_q, &n_q, &q);
      if (n_q * m_q !=0) // the matrix can be empty
	{
#ifdef DEBUG
	  sciprint("DEBUG: Q_IN is not sparse\n");
#endif
	  Q_matrix = new CoinPackedMatrix();
	  Q_matrix->setDimensions(nrows,ncols);
	  
	  nz = n_q * m_q;

	  if (q==NULL)
	    {
	      Scierror(999,"%s: invalid value of matrix q",fname);
	      return 0;
	    }
	  
	  for(i=0; i<m_q; i++)
	    {
	      for(j=0; j<n_q; j++)
		{
		  if (*(q+i+j*m_q) != 0) Q_matrix->modifyCoefficient(i,j,*(q+i+j*m_q));
		}
	    }
	  modelSimplex->loadQuadraticObjective(*Q_matrix);
	}
    }
  else
    {
      getAllocatedSparseMatrix(pvApiCtx, q_addr, &S_Q.m, &S_Q.n, &S_Q.nel, &S_Q.mnel, &S_Q.icol, &S_Q.R);

      if (S_Q.m * S_Q.n!=0) // the matrix can be empty
	{
#ifdef DEBUG
	  sciprint("DEBUG: Q_IN is sparse\n");
#endif
	  Q_matrix = new CoinPackedMatrix();
	  Q_matrix->setDimensions(nrows,ncols);
	  
	  nz = S_Q.nel;
	  
	  count = 0;
	  for(i=0;i<S_Q.m;i++)
	    {
	      if (S_Q.mnel[i]!=0) 
		{
		  for(j=0;j<S_Q.mnel[i];j++)
		    {
		      count++;
		      Q_matrix->modifyCoefficient(i,S_Q.icol[count-1]-1,S_Q.R[count-1]);
		    }
		}
	    }
	  modelSimplex->loadQuadraticObjective(*Q_matrix);
	}

      freeAllocatedSparseMatrix(S_Q.mnel, S_Q.icol, S_Q.R);
    }

  /////////////////
  // Get options //
  /////////////////

#ifdef DEBUG
  sciprint("DEBUG: get options\n");
#endif

  // Default settings
  // optimisation direction: 1 - minimize, -1 - maximize, 0 - ignore
  int     * param_in_addr = NULL;
  int     maxnumiterations = 99999999, loglevel = 0, clp_loglevel = 0, optim_dir = 1;
  int     stoponfirstsol = 0;
  double  maxnumseconds = 3600.0, primaltolerance = 1e-7, dualtolerance = 1e-7;
  int     tmp_res;

  initPList(pvApiCtx, PARAM_IN, &param_in_addr);
  if (!checkPList(pvApiCtx, param_in_addr))
    {
      Scierror(999, "%s: argument n° %d is not a plist\n", fname, PARAM_IN);

      return 0;
    }

  //////////////////////////////
  // Set some simplex options //
  //////////////////////////////

  // We set some global parameters
  getIntInPList(pvApiCtx, param_in_addr, "maxnumiterations", &tmp_int, &tmp_res, 1000, Log, CHECK_NONE);
  if (tmp_res!=-1) modelSimplex->setMaximumIterations(maxnumiterations);

  getDoubleInPList(pvApiCtx, param_in_addr, "maxnumseconds", &tmp_double, &tmp_res, 1000, Log, CHECK_NONE);
  if (tmp_res!=-1) modelSimplex->setMaximumSeconds(maxnumseconds);

  getDoubleInPList(pvApiCtx, param_in_addr, "primaltolerance", &tmp_double, &tmp_res, 1e-8, Log, CHECK_NONE);
  if (tmp_res!=-1) modelSimplex->setPrimalTolerance(primaltolerance);

  getDoubleInPList(pvApiCtx, param_in_addr, "dualtolerance", &tmp_double, &tmp_res, 1e-8, Log, CHECK_NONE);
  if (tmp_res!=-1) modelSimplex->setDualTolerance(dualtolerance);

  getIntInPList(pvApiCtx, param_in_addr, "stoponfirstsol", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) SciCbcEventHandler.setStopOnFirstSol((bool)stoponfirstsol);

  // There are several log levels. Setting the log level to be i produces the log messages for level i and all levels less than i.
  // * Log Level 0: Switches off all CBC messages, but one.
  // * Log Level 1: The default.
  // * Log Level 2: Substantial amount of information, e.g., message 15 is generated once per node. 
  //                Can be useful when the evaluation at each node is slow.
  // * Log Level 3: Tremendous amount of information, e.g., multiple messages per node. 

  getIntInPList(pvApiCtx, param_in_addr, "verbose", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  printer = new DerivedHandler();
  if (tmp_res!=-1) printer->setLogLevel(loglevel);

  getIntInPList(pvApiCtx, param_in_addr, "clpverbose", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  clpprinter = new DerivedHandler();
  if (tmp_res!=-1) clpprinter->setLogLevel(clp_loglevel);
  modelSimplex->passInMessageHandler(clpprinter);	

  getIntInPList(pvApiCtx, param_in_addr, "optim_dir", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) modelSimplex->setOptimizationDirection(optim_dir);

  ////////////////////////////////////
  // Do we have special constraints //
  ////////////////////////////////////

  int nb_special_constr      = 0, * special_constr_addr      = NULL;
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
  double * elem_which = NULL, * elem_weight = NULL, * elem_type = NULL, * elem_column = NULL;
  double * elem_range = NULL, * elem_id_obj = NULL, * elem_length = NULL, * elem_id = NULL;
  double * elem_clique_type = NULL;

  _SciErr = getVarAddressFromPosition(pvApiCtx, SPECIAL_CONSTR_IN, &special_constr_addr);
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
      
      if (nb_special_which != 0) special_stored = true;
    }
  else
    {
      special_stored = false;
    }

  ///////////////////
  // Set presolver //
  ///////////////////

  // Extract of the CoinPresolveAction.hpp doxygen documentation:

  // Abstract base class of all presolve routines.

  // The details will make more sense after a quick overview of the grand plan:
  // A presolve object is handed a problem object, which it is expected to modify in some useful way.
  // Assuming that it succeeds, the presolve object should create a postsolve object, i.e., an object that contains 
  // instructions for backing out the presolve transform to recover the original problem.
  // These postsolve objects are accumulated in a linked list, with each successive presolve action adding its postsolve action
  //  to the head of the list.
  // The end result of all this is a presolved problem object, and a list of postsolve objects. The presolved problem object is 
  // then handed to a solver for optimization, and the problem object augmented with the results.
  // The list of postsolve objects is then traversed. Each of them (un)modifies the problem object, with the end result being
  // the original problem, augmented with solution information.
  //
  // The problem object representation is CoinPrePostsolveMatrix and subclasses. Check there for details. 
  // The CoinPresolveAction class and subclasses represent the presolve and postsolve objects.
  //
  // In spite of the name, the only information held in a CoinPresolveAction object is the information needed to postsolve
  // (i.e., the information needed to back out the presolve transformation). 
  // This information is not expected to change, so the fields are all const.
  //
  // A subclass of CoinPresolveAction, implementing a specific pre/postsolve action, is expected to declare a static function
  // that attempts to perform a presolve transformation. 
  // This function will be handed a CoinPresolveMatrix to transform, and a pointer to the head of the list of postsolve objects.
  // If the transform is successful, the function will create a new CoinPresolveAction object, link it at the head of the list
  // of postsolve objects, and return a pointer to the postsolve object it has just created.
  // Otherwise, it should return 0. It is expected that these static functions will be the only things that can create new
  // CoinPresolveAction objects; this is expressed by making each subclass' constructor(s) private.
  //
  // Every subclass must also define a postsolve method. This function will be handed a CoinPostsolveMatrix to transform.
  //
  // It is the client's responsibility to implement presolve and postsolve driver routines. See OsiPresolve for examples.
  //
  // Note:
  // Since the only fields in a CoinPresolveAction are const, anything one can do with a variable declared CoinPresolveAction* 
  // can also be done with a variable declared const CoinPresolveAction* It is expected that all derived subclasses
  // of CoinPresolveAction also have this property. 

  // Actions performed by ClpPresolve:
  //
  // 1 - doubleton_action:
  // Solve ax+by=c for y and substitute y out of the problem.
  // This moves the bounds information for y onto x, making y free and allowing us to substitute it away.
  //	   a x + b y = c
  //	   l1 <= x <= u1
  //	   l2 <= y <= u2	==>
  //	  
  //	   l2 <= (c - a x) / b <= u2
  //	   b/-a > 0 ==> (b l2 - c) / -a <= x <= (b u2 - c) / -a
  //	   b/-a < 0 ==> (b u2 - c) / -a <= x <= (b l2 - c) / -a
  
  // 2 - dupcol_action:
  // Detect and remove duplicate columns.
  // The general technique is to sum the coefficients a_(*,j) of each column. Columns with identical sums are duplicates. 
  // The obvious problem is that, e.g., [1 0 1 0] and [0 1 0 1] both add to 2. To minimize the chances of false positives, 
  // the coefficients of each row are multipled by a random number r_i, so that we sum r_i*a_ij.
  // Candidate columns are checked to confirm they are identical. Where the columns have the same objective coefficient,
  // the two are combined. If the columns have different objective coefficients, complications ensue.
  // In order to remove the duplicate, it must be possible to fix the variable at a bound. 

  // 3 - duprow_action:
  // Detect and remove duplicate rows.
  // The algorithm to detect duplicate rows is as outlined for dupcol_action.
  // If the feasible interval for one constraint is strictly contained in the other, the tighter (contained) constraint is kept.
  // If the feasible intervals are disjoint, the problem is infeasible. 
  // If the feasible intervals overlap, both constraints are kept.
  // duprow_action is definitely a work in progress; postsolve is unimplemented. 

  // 4 - drop_empty_cols_action:
  // Physically removes empty columns in presolve, and reinserts empty columns in postsolve.
  // Physical removal of rows and columns should be the last activities performed during presolve. 
  // Do them exactly once. 
  // The row-major matrix is not maintained by this transform.
  // To physically drop the columns, CoinPrePostsolveMatrix::mcstrt_ and CoinPrePostsolveMatrix::hincol_ are compressed,
  // along with column bounds, objective, and (if present) the column portions of the solution.
  // This renumbers the columns. drop_empty_cols_action::presolve will reconstruct CoinPresolveMatrix::clink_.
  // Todo: Confirm correct behaviour with solution in presolve. 

  // 5 - drop_empty_rows_action:
  // Physically removes empty rows in presolve, and reinserts empty rows in postsolve.
  // Physical removal of rows and columns should be the last activities performed during presolve. 
  // Do them exactly once. The row-major matrix is not maintained by this transform.
  // To physically drop the rows, the rows are renumbered, excluding empty rows.
  // This involves rewriting CoinPrePostsolveMatrix::hrow_ and compressing the row bounds and (if present)
  // the row portions of the solution.
  // Todo: Confirm behaviour when a solution is present in presolve. 

  // 6 - remove_fixed_action:
  // Excise fixed variables from the model.
  // Implements the action of removing one or more fixed variables x_j from the model by substituting the value sol_j in 
  // each constraint. Specifically, for each constraint i where a_ij != 0, rlo_i and rup_i are adjusted by -a_ij*sol_j
  // and a_ij is set to 0.
  // There is an implicit assumption that the variable already has the correct value. If this isn't true,
  // corrections to row activity may be incorrect. If you want to guard against this possibility, consider make_fixed_action.
  // Actual removal of the column from the matrix is handled by drop_empty_cols_action. 
  // Correction of the objective function is done there. 

  // 7 - make_fixed_action:
  // Fix a variable at a specified bound.
  // Implements the action of fixing a variable by forcing both bounds to the same value and forcing the value of the variable to match.
  // If the bounds are already equal, and the value of the variable is already correct, consider remove_fixed_action. 

  // 8 - forcing_constraint_action:
  // Detect and process forcing constraints and useless constraints.
  // A constraint is useless if the bounds on the variables prevent the constraint from ever being violated.
  // A constraint is a forcing constraint if the bounds on the constraint force the value of an involved variable to one
  // of its bounds. A constraint can force more than one variable. 

  // 8 - implied_free_action:
  // Detect and process implied free variables.
  // Consider a singleton variable x (i.e., a variable involved in only one constraint). Suppose that the bounds on that constraint, 
  // combined with the bounds on the other variables involved in the constraint, are such that even the worst case values of
  // the other variables still imply bounds for x which are tighter than the variable's original bounds. 
  // Since x can never reach its upper or lower bounds, it is an implied free variable. 
  // Both x and the constraint can be deleted from the problem.
  // The transform also handles more complicated variations, where x is not a singleton. 

  // 10 - isolated_constraint_action:

  // 11 - slack_doubleton_action:
  // Convert an explicit bound constraint to a column bound.
  // This transform looks for explicit bound constraints for a variable and transfers the bound to the appropriate
  // column bound array. The constraint is removed from the constraint system. 

  // 12 - slack_singleton_action:
  // For variables with one entry.
  // If we have a variable with one entry and no cost then we can transform the row from E to G etc. 
  // If there is a row objective region then we may be able to do this even with a cost. 

  // 13 - subst_constraint_action:

  // 14 - do_tighten_action:
  // it decides which columns can be made fixed and calls make_fixed_action::presolve

  // 15 - tripleton_action:
  // We are only going to do this if it does not increase number of elements?.
  // It could be generalized to more than three but it seems unlikely it would help.
  // As it is adapted from doubleton icoly is one dropped. 

  // 16 - useless_constraint_action:

  // 17 - drop_zero_coefficients_action:
  // Removal of explicit zeros.
  // The presolve action for this class removes explicit zeros from the constraint matrix. The postsolve action puts them back. 

  // 18 - remove_dual_action:
  // Attempt to fix variables by bounding reduced costs.
  // The reduced cost of x_j is d_j = c_j - y*a_j (1). Assume minimization, so that at optimality d_j >= 0 for x_j nonbasic
  // at lower bound, and d_j <= 0 for x_j nonbasic at upper bound.
  // For a slack variable s_i, c_(n+i) = 0 and a_(n+i) is a unit vector, hence d_(n+i) = -y_i. 
  // If s_i has a finite lower bound and no upper bound, we must have y_i <= 0 at optimality. 
  // Similarly, if s_i has no lower bound and a finite upper bound, we must have y_i >= 0.
  // For a singleton variable x_j, d_j = c_j - y_i*a_ij. Given x_j with a single finite bound, we can bound d_j greater
  // or less than 0 at optimality, and that allows us to calculate an upper or lower bound on y_i (depending on the bound
  // on d_j and the sign of a_ij).
  // Now we have bounds on some subset of the y_i, and we can use these to calculate upper and lower bounds on the d_j,
  // using bound propagation on (1). If we can manage to bound some d_j as strictly positive or strictly negative, 
  // then at optimality the corresponding variable must be nonbasic at its lower or upper bound, 
  // respectively. If the required bound is lacking, the problem is unbounded.
  // There is no postsolve object specific to remove_dual_action, but execution will queue postsolve actions for any variables
  // that are fixed. 

#ifdef DEBUG
  sciprint("Set presolver\n");
#endif

  // clppresolve is desactivated when we use special constraints because clppresolve can change the size of the matrix A
  // (removing constraints, removing unused constraints, etc ...)
  ClpPresolve  pinfo;
  ClpSimplex * modelSimplex2 = NULL;
  
  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "clppresolve");
  if ((tmp_int!=-1) && !special_stored)
    {
      TRYCATCH(modelSimplex2 = pinfo.presolvedModel(*modelSimplex))
    }
  else
    {
#ifdef DEBUG
      sciprint("presolver not set\n");
#endif
      modelSimplex2 = modelSimplex;
    }

  // The presolver can fail
  if (modelSimplex2==NULL)
    {
      if (clp_loglevel) 
	{
	  sciprint("%s: presolver says model is infeasible or preprocess is not selected. Bypassing presolver\n", fname);
	}
      modelSimplex2 = modelSimplex;
    }

  clpprinter->setLogLevel(clp_loglevel);
  modelSimplex2->passInMessageHandler(clpprinter);

  ////////////////
  // PreProcess //
  ////////////////

  // Class for preProcessing and postProcessing.
  //
  // While cuts can be added at any time in the tree, some cuts are actually just stronger versions of existing constraints. 
  // In this case they can replace those constraints rather than being added as new constraints. 
  // This is awkward in the tree but reasonable at the root node.
  // This is a general process class which uses other cut generators to strengthen constraints, establish that constraints
  // are redundant, fix variables and find relationships such as x + y == 1.
  // Presolve will also be done.
  // If row names existed they may be replaced by R0000000 etc 

  // YC: Add cut generators.Maybe we can set the cuts before the call to CglPreProcess and then get the generators stored
  // in model and feed these generators in process via addCutGenerator.
  // Seems to be difficults ...

#ifdef DEBUG
  sciprint("Set preprocess\n");
#endif

  int  cgl_nbpass = 5, cgl_type = 2;
  CglPreProcess process;
  OsiSolverInterface * solver4 = NULL;
  OsiClpSolverInterface solver3(modelSimplex2);

  clpprinter->setLogLevel(clp_loglevel);
  solver3.passInMessageHandler(clpprinter);

#ifdef HANDLE_CTRLC
  ///////////////////
  // Handle Ctrl+C //
  ///////////////////

  solver3.getModelPtr()->passInEventHandler(&SciClpEventHandler);
#endif

  // This is a class which uses other cut generators to strengthen cuts, establish that some cuts are redundant, 
  // fix variables and find relationships such as x + y == 1.
  //
  // While cuts can be added at any time in the tree, some cuts are actually just stronger versions of existing ones. 
  // They can thus replace the existing cuts rather than being added as new cuts. 
  // This is awkward in the tree but reasonable at the root node. 

  // cglpreprocess is desactivated when we use special constraints because cglppreprocess can change the size of the matrix A
  // (removing constraints, removing unused constraints, etc ...)

  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cglpreprocess");
  if ((tmp_int!=-1) && !special_stored)
    {
      getIntInPList(pvApiCtx, param_in_addr, "cglpreprocess_other", &cgl_nbpass, &tmp_res, 5, Log, CHECK_NONE);
      if (tmp_res!=-1) solver4 = process.preProcess(solver3,false,cgl_nbpass);
      getIntInPList(pvApiCtx, param_in_addr, "cglpreprocess_nondefault", &cgl_nbpass, &tmp_res, 5, Log, CHECK_NONE);
      if (tmp_res!=-1)
	{
	  getIntInPList(pvApiCtx, param_in_addr, "cglpreprocess_nondefault_type", &cgl_type, &tmp_res, 2, Log, CHECK_NONE);
	  CglGomory                generator_preprocess_gomory;
	  CglRedSplit              generator_preprocess_redsplit;
	  CglMixedIntegerRounding  generator_preprocess_mixgen;
	  CglMixedIntegerRounding2 generator_preprocess_mixgen2;
	  CglClique                generator_preprocess_clique;
	  CglTwomir                generator_preprocess_twomir;
	  
	  process.addCutGenerator(&generator_preprocess_gomory);
	  process.addCutGenerator(&generator_preprocess_redsplit);
	  process.addCutGenerator(&generator_preprocess_mixgen);
	  process.addCutGenerator(&generator_preprocess_mixgen2);
	  process.addCutGenerator(&generator_preprocess_clique);
	  process.addCutGenerator(&generator_preprocess_twomir);

	  // 1 -> clique, 2 -> SOS
	  TRYCATCH(solver4 = process.preProcessNonDefault(solver3,cgl_type,cgl_nbpass))
	  
	  // YC:
	  // Once the particular constraints has been detected by process, we must add objects for these
	  // constraints manually
	}

      if (clp_loglevel) sciprint("numberSOS = %d\n",process.numberSOS());

      if (!solver4)
	{
	  if (clp_loglevel) 
	    {
	      sciprint("%s: preprocessing says model is infeasible or preprocess is not selected. Bypassing preprocessing\n", fname);
	    }
	  solver4 = &solver3;
	}
      else
	{
	  for(i=0;i<process.numberSOS();i++)
	    {
	      if (clp_loglevel) 
		{
		  sciprint("%s: %d - whichSOS = %d typeSOS = %d\n", fname, i, process.whichSOS()[i], process.typeSOS()[i]);
		}
	    }
	}

      clpprinter->setLogLevel(clp_loglevel);
      solver4->passInMessageHandler(clpprinter);	

      TRYCATCH(solver4->resolve());
    }
  else
    {
      // The option has not been selected
      solver4 = &solver3;
      clpprinter->setLogLevel(clp_loglevel);
      solver4->passInMessageHandler(clpprinter);	

#ifdef DEBUG
      sciprint("presolver not set\n");
#endif
    }

  /////////////////////////////
  // Launch first resolution //
  // before branch and bound //
  /////////////////////////////

  solver4->setHintParam(OsiDoReducePrint,true,OsiHintTry); // from sample2.cpp
  
  TRYCATCH(solver4->initialSolve())

  status = solver4->isProvenPrimalInfeasible() || solver4->isProvenDualInfeasible() ||
    solver4->isPrimalObjectiveLimitReached() || solver4->isDualObjectiveLimitReached() ||
    solver4->isIterationLimitReached();

  if ((status!=0) && clp_loglevel)
      {
	sciprint("%s: initialSolve fails. Problem is not feasible.\n", fname);
      }
  
  CbcModel model(*solver4);

  // CbcMain0 allows to add some "generic" options to model
  getIntInPList(pvApiCtx, param_in_addr, "cbcmaininit", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) CbcMain0(model);

  // CbcMain0 doesn't deal with passInMessageHandler (to be verified). So, we now set the printer handler.
  printer->setLogLevel(loglevel);
  model.passInMessageHandler(printer);

#ifdef HANDLE_CTRLC
  ////////////////////////////////////////////////
  // Pass in the event handler to handle ctrl+c //
  ////////////////////////////////////////////////

  model.passInEventHandler(&SciCbcEventHandler);
#endif

  /////////////////////////////////
  // Process Special constraints //
  /////////////////////////////////

  //////////////////////////////
  // Read Special Constraints //
  // SOS, Clique, nWay        //
  //////////////////////////////

  // CbcClique: Le choix fait sur une variable implique un choix sur les autres variables
  // x1 + x2 + x3 = 1
  // si x1 = 1 alors x2 = x3 = 0
  // CbcClique(CbcModel *, int CliqueType, int nbMember, const int * which, const char * type, int id, int slack = -1)
  // CliqueType: 0 alors <= , 1 alors ==
  // Type: 0 alors non SOS, 1 à nbMember - 1 SOS
  
  // CbcLotSize: une variable peut prendre des valeurs du type 0, 1.2 à 2.1, 3, 4.0 à 4.2
  // CbcLotSize(CbcModel *, int Column, int NbPts, const double * Pts, bool range = false)
  // Column: sur quelle variable on applique cette méthode
  // NbPts: nombre de valeurs de la variables. Si range = true alors la liste doit être du type (min1 max1 min2 max2 ...)
  
  // CbcSOS: contraintes de type Special Ordered Sets
  // On définit un ensemble de variables parmis lesquels un certain nombre peuvent être non nulles
  // SOS de type 1 impliquant x1, x2, x3 dans un problème à 5 variables
  // si on affecte x1 alors x2 et x3 seront = 0
  // si on affecte x2 alors x1 et x3 seront = 0
  // si on affecte x3 alors x1 et x2 seront = 0
  // SOS de type 2 impliquant x1, x2, x3 dans un problème à 5 variables
  // si on affecte x1 et x2 alors x3 sera = 0
  // si on affecte x1 et x3 alors x2 sera = 0
  // si on affecte x3 et x2 alors x1 sera = 0

  /////////////////////////////////////////
  // Get the list of special constraints //
  /////////////////////////////////////////

#ifdef DEBUG
  sciprint("Set specials\n");
#endif
  
  int type, length;
  char * tmp_char = NULL;
  int * which = NULL;
  // We must get something from the stack for these parameters so as to be able to CreateVar for the outputs

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
  int nb_elem = nb_special_id;

  if (special_stored)
    {
      specialObjects = new CbcObject * [nb_elem]; 
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

	  which = (int *)MALLOC(((int)*elem_length)*sizeof(int));
	  for(j=0;j<(int)*elem_length;j++) which[j] = (int)elem_which[j];

#ifdef DEBUG
	  sciprint("%d - type = %d - column = %d range = %d id_obj = %d length = %d id = %d clique_type = %d\n", 
		   i,
		   (int)*elem_type,
		   (int)*elem_column,
		   (int)*elem_range,
		   (int)*elem_id_obj,
		   (int)*elem_length, 
		   (int)*elem_id, 
		   (int)*elem_clique_type);

	  for(j=0;j<(int)*elem_length;j++) sciprint("element %d: which = %d weight = %f\n", j, (int)elem_which[j], elem_weight[j]);
#endif

	  // Process element i
	  switch((int)*elem_id_obj)
	    {
	    case 0: // CbcSOS
#ifdef DEBUG
	      sciprint("adding CbcSOS object %d\n",i);
#endif
 	      // Take off integers if sos type is above 1.
 	      if ((int)*elem_type==2)
 		{
 		  for(j=0;j<(int)*elem_length;j++) model.solver()->setContinuous(which[j]);
 		}

	      specialObjects[i] = new CbcSOS(&model,(int)*elem_length,(const int *)which,(const double *)elem_weight,(int)*elem_id,(int)*elem_type);
	      break;
	    case 1: // CbcClique
#ifdef DEBUG
	      sciprint("adding CbcClique object %d\n",i);
#endif
	      // CbcClique(CbcModel * model, int cliqueType, int numberMembers, const int * which, const char * type, int identifier,int slack=-1);
	      // members_: Members (indices in range 0 ... numberIntegers_-1)
	      // type_:    Type of each member 0=SOS, 1 =clique
	      // cliqueType_: Clique type - 0 <=, 1 ==
	      // slack_: Which one is slack (if any) sequence within this set

	      tmp_char = new char[length];
	      
	      if (type==0) memset(tmp_char,0,length*sizeof(char)); // Checked in CbcBranchActual.cpp - class CbcClique
	      else         memset(tmp_char,1,length*sizeof(char));

	      specialObjects[i] = new CbcClique(&model,(int)*elem_clique_type,(int)*elem_length,(const int *)which,tmp_char,(int)*elem_id);

	      if (tmp_char) delete [] tmp_char;

	      break;
	    case 2: // CbcLotsize
#ifdef DEBUG
	      sciprint("adding CbcLotsize object %d\n",i);
#endif
	      specialObjects[i] = new CbcLotsize(&model,(int)*elem_column,(int)*elem_length,(const double *)elem_weight,(int)*elem_range);
	      break;
	    case 3: // CbcNWay
#ifdef DEBUG
	      sciprint("adding CbcNWay object %d\n",i);
#endif
	      specialObjects[i] = new CbcNWay(&model,(int)*elem_length,(const int *)which,(int)*elem_id);
	      break;
	    default:
	      Scierror(999,"%s: Wrong id_obj\n",fname);
	      return 0;
	    }
	}

      if (which) FREE(which);

      // First we clean model
#ifdef DEBUG
      sciprint("Number of pre-existing CbcObjects: %d - nb_elem = %d\n",model.numberObjects(), nb_elem);
#endif
      // addObjects will add a CbcSimpleInteger for each integer variables.
      // so, we will have nb_elem + nb_integer objects
      model.addObjects(nb_elem,specialObjects);
      // Ensure all attached objects (OsiObjects, heuristics, and cut
      // generators) point to this model.
      model.synchronizeModel();

      for(i=0;i<nb_elem;i++) delete specialObjects[i];
      delete [] specialObjects;

//       int * priority = new int[nb_elem];
//       // Set Special object priorities high
//       CoinFillN(priority,nb_elem,1);

//       model.passInPriorities(priority,true);
//       delete [] priority;


      getIntInPList(pvApiCtx, param_in_addr, "cglpreprocess_nondefault", &cgl_nbpass, &tmp_res, 5, Log, CHECK_NONE);
      if (tmp_res!=-1)
	{
	  getIntInPList(pvApiCtx, param_in_addr, "cglpreprocess_nondefault_type", &cgl_type, &tmp_res, 2, Log, CHECK_NONE);
	  // Now, add SOS detected by CglPreProcess
	  if (cgl_type==2)
	    {
	      specialObjects = new CbcObject * [process.numberSOS()]; 

#ifdef DEBUG
	      sciprint("number of SOS detected by preprocess: %d\n", process.numberSOS());
#endif

	      for(i=0;i<process.numberSOS();i++)
		{
		  int      sos_length = process.startSOS()[i+1] - process.startSOS()[i] - 1;
		  int    * sos_which  = new int[sos_length];
		  double * sos_weight = new double[sos_length];
		  int      sos_type   = process.typeSOS()[i];

		  memcpy(sos_which, &process.whichSOS()[process.startSOS()[i]], sos_length);
		  memcpy(sos_weight,&process.weightSOS()[process.startSOS()[i]],sos_length);
		  specialObjects[i] = new CbcSOS(&model,sos_length,sos_which,sos_weight,i,sos_type);

#ifdef DEBUG
		  sciprint("%d - typeSOS = %d - startSOS = %d\n", i, process.typeSOS()[i], process.startSOS()[i]);
#endif

		  delete [] sos_which;
		  delete [] sos_weight;
		}
	      model.addObjects(process.numberSOS(),specialObjects);
	      model.synchronizeModel();

	      for(i=0;i<process.numberSOS();i++) delete specialObjects[i];
	      delete [] specialObjects;
	    }
#ifdef DEBUG
	  sciprint("Number of CbcObjects: %d\n",model.numberObjects());
#endif
	}

    } // end special stored

  //
  // Now we set some Cbc classes to fine tune the behavior of CBC
  // CbcBranch...     These classes define the nature of MIP's discontinuity. The simplest discontinuity is a 
  //                  variable which must take an integral value. Other types of discontinuities exist, e.g., lot-sizing variables.
  //
  // CbcNode          This class decides which variable/entity to branch on next. Even advanced users will probably 
  //                  only interact with this class by setting CbcModel parameters ( e.g., priorities).
  //
  // CbcTree          All unsolved models can be thought of as being nodes on a tree where each node (model) can
  //                  branch two or more times. The interface with this class is helpful to know, but the user can pretty
  //                  safely ignore the inner workings of this class.
  //
  // CbcCompare...    These classes are used in determine which of the unexplored nodes in the tree to consider next.
  //                  These classes are very small simple classes that can be tailored to suit the problem.
  //
  // CglCutGenerators Any cut generator from CGL can be used in CBC. The cut generators are passed to CBC with parameters
  //                  which modify when each generator will be tried. All cut generators should be tried to determine which are
  //                  effective. Few users will write their own cut generators.
  //
  // CbcHeuristics    Heuristics are very important for obtaining valid solutions quickly. Some heuristics are available,
  //                  but this is an area where it is useful and interesting to write specialized ones. 

  // CbcCompareBase   Controls which node on the tree is selected.
  //                  The default is CbcCompareDefault. Other comparison classes in CbcCompareActual.hpp include CbcCompareDepth
  //                  and CbcCompareObjective. Experimenting with these classes and creating new compare classes is easy.
  // CbcCutGenerator  A wrapper for CglCutGenerator with additional data to control when the cut generator is invoked during the 
  //                  tree search.
  //                  Other than knowing how to add a cut generator to CbcModel, there is not much the average user needs to know
  //                  about this class. However, sophisticated users can implement their own cut generators.
  // CbcHeuristic     Heuristic that attempts to generate valid MIP-solutions leading to good upper bounds.
  //                  Specialized heuristics can dramatically improve branch-and-cut performance. As many different heuristics as
  //                  desired can be used in CBC. Advanced users should consider implementing custom heuristics when tackling 
  //                  difficult problems.
  // CbcObject        Defines what it means for a variable to be satisfied. Used in branching.
  //                  Virtual class. CBC's concept of branching is based on the idea of an "object". An object has 
  //                  (i) a feasible region,
  //                  (ii) can be evaluated for infeasibility,
  //                  (iii) can be branched on, e.g., a method of generating a branching object, which defines an up branch and
  //                        a down branch, and 
  //                  (iv) allows comparison of the effect of branching.
  //                  Instances of objects include CbcSimpleInteger, CbcSimpleIntegerPseudoCosts, CbcClique, CbcSOS (type 1 and 2),
  //                  CbcFollowOn, and CbcLotsize. 

  ///////////////////////////
  // Set the compare class //
  ///////////////////////////

  // CbcCompareDepth     This will always choose the node deepest in tree.
  //                     It gives minimum tree size but may take a long time to find the best solution.
  // CbcCompareObjective This will always choose the node with the best objective value. This may give a very large tree.
  //                     It is likely that the first solution found will be the best and the search should finish soon
  //                     after the first solution is found.
  // CbcCompareDefault   This is designed to do a mostly depth-first search until a solution has been found. 
  //                     It then use estimates that are designed to give a slightly better solution. If a reasonable number of
  //                     nodes have been explored (or a reasonable number of solutions found), then this class will adopt a
  //                     breadth-first search (i.e., making a comparison based strictly on objective function values) 
  //                     unless the tree is very large, in which case it will revert to depth-first search. 
  //                     A better description of CbcCompareUser is given below.
  // CbcCompareEstimate  When pseudo costs are invoked, CBC uses the psuedo costs to guess a solution.
  //                     This class uses the guessed solution. 
      
#ifdef DEBUG
  sciprint("Set options\n");
#endif

  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cbccompareuser");
  if (tmp_int!=-1)
    {
      CbcCompareDefault compare1;
      getDoubleInPList(pvApiCtx, param_in_addr, "cbccompareuser_weight", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) compare1.setWeight(tmp_double);
      model.setNodeComparison(compare1);
    }

  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cbccomparedefault");
  if (tmp_int!=-1)
    {
      CbcCompareDefault compare2;
      getDoubleInPList(pvApiCtx, param_in_addr, "cbccomparedefault_weight", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) compare2.setWeight(tmp_double);
      getIntInPList(pvApiCtx, param_in_addr, "cbccomparedefault_breadthdepth", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) compare2.setBreadthDepth(tmp_int);
      model.setNodeComparison(compare2);
    }

  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cbccomparedepth");
  if (tmp_int!=-1)
    {
      CbcCompareDepth compare3;
      model.setNodeComparison(compare3);
    }

  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cbccompareestimate");
  if (tmp_int!=-1)
    {
      CbcCompareEstimate compare4;
      model.setNodeComparison(compare4);
    }

  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cbccompareobjective");
  if (tmp_int!=-1)
    {
      CbcCompareObjective compare5;
      model.setNodeComparison(compare5);
    }

  ///////////////////////////////
  // Set the Cgl cut generator //
  ///////////////////////////////

  // This is a simplification of probing ideas put into OSL about ten years ago.
  // The only known documentation is a copy of a talk handout - we think Robin Lougee-Heimer has a copy!
  //
  // For selected integer variables (e.g. unsatisfied ones) the effect of setting them up or down is investigated. 
  // Setting a variable up may in turn set other variables (continuous as well as integer). There are various possible results:
  //
  // 1) It is shown that the problem is infeasible (this may also be because the objective function or reduced costs 
  //    show worse than best solution value). If the other way is feasible we can generate a column cut (and continue probing), 
  //    otherwise we can say that the problem is infeasible. 
  //
  // 2) If both ways are feasible, it can happen that
  //    2.1) setting x to 0 implies that y must be set to 1 and
  //    2.2) setting x to 1 implies that y must be set to 1
  //         yielding again a column cut. (2.2 is not done in this code as there is no mechanism for returning the information.) 
  //
  // More common is that
  //
  //    2.3) setting x to 0 implies that y must be set to 1 and
  //    2.4) setting x to 1 implies that y must be set to 0
  //
  // so we can substitute for y which might lead later to more powerful cuts.
  //
  // 3) When setting x to 1, a constraint went slack by c. We can tighten the constraint ax + .... <= b (where a may be zero)
  //    to (a+c)x + .... <= b. If this cut is violated then is generated. 
  //
  // 4) Similarly we can generate implied disaggregation cuts 
  //
  // Note - differences to cuts in OSL.
  //
  // a) OSL had structures intended to make this faster.
  // b) The "chaining" in 2) was done
  // c) Row cuts modified original constraint rather than adding cut
  // d) This code can cope with general integer variables. 
  //
  // Parameters:
  //
  // The mode options are:
  //
  // 0) Only unsatisfied integer variables will be looked at. If no information exists for that variable then probing
  //    will be done so as a by-product you "may" get a fixing or infeasibility. This will be fast and is only available 
  //    if a snapshot exists (otherwise as 1). The bounds in the snapshot are the ones used.
  // 1) Look at unsatisfied integer variables, using current bounds. Probing will be done on all looked at.
  // 2) Look at all integer variables, using current bounds. Probing will be done on all 
  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cglprobing");
  if (tmp_int!=-1)
    {
      CglProbing generator1;

      getIntInPList(pvApiCtx, param_in_addr, "cglprobing_mode", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator1.setMode(tmp_int);
      // Set maximum number of passes per node.
      getIntInPList(pvApiCtx, param_in_addr, "cglprobing_maxpass", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator1.setMaxPass(tmp_int);
      // Set log level - 0 none, 1 - a bit, 2 - more details.
      getIntInPList(pvApiCtx, param_in_addr, "cglprobing_loglevel", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator1.setLogLevel(tmp_int);
      // Set maximum number of unsatisfied variables to look at.
      getIntInPList(pvApiCtx, param_in_addr, "cglprobing_maxprobe", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator1.setMaxProbe(tmp_int);
      // Set maximum number of variables to look at in one probe.
      getIntInPList(pvApiCtx, param_in_addr, "cglprobing_maxlook", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator1.setMaxLook(tmp_int);
      // Set maximum number of elements in row for it to be considered.
      getIntInPList(pvApiCtx, param_in_addr, "cglprobing_maxelements", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator1.setMaxElements(tmp_int);
      // Set maximum number of passes per node (root node).
      getIntInPList(pvApiCtx, param_in_addr, "cglprobing_maxpassroot", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator1.setMaxPassRoot(tmp_int);
      // Set maximum number of unsatisfied variables to look at (root node).
      getIntInPList(pvApiCtx, param_in_addr, "cglprobing_maxproberoot", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator1.setMaxProbeRoot(tmp_int);
      // Set maximum number of variables to look at in one probe (root node).
      getIntInPList(pvApiCtx, param_in_addr, "cglprobing_maxlookroot", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator1.setMaxLookRoot(tmp_int);
      // Set maximum number of elements in row for it to be considered (root node).
      getIntInPList(pvApiCtx, param_in_addr, "cglprobing_maxelementsroot", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator1.setMaxElementsRoot(tmp_int);
      // Stop or restart row cuts (otherwise just fixing from probing)
      // Set 0 no cuts, 1 just disaggregation type, 2 coefficient ( 3 both).
      getIntInPList(pvApiCtx, param_in_addr, "cglprobing_rowcuts", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator1.setRowCuts(tmp_int);
      // Whether use objective as constraint
      // Set 0 don't 1 do -1 don't even think about it.
      getIntInPList(pvApiCtx, param_in_addr, "cglprobing_usingobjective", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator1.setUsingObjective(tmp_int);

      getIntInPList(pvApiCtx, param_in_addr, "cglprobing_priority", &tmp_int, &tmp_res, -1, Log, CHECK_NONE);
      model.addCutGenerator(&generator1,tmp_int,"Probing");
    }

  // Generates mixed integer Gomory Cuts.
  // Notes: By default, does not generate cuts with more than 50 non zero coefficients. 
  // To get more dense cuts, modify the parameter Limit. See the method setLimit(). 
  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cglgomory");
  if (tmp_int!=-1)
    {
      CglGomory generator2;

      // Change limit on how many variables in cut (default 50)
      getIntInPList(pvApiCtx, param_in_addr, "cglgomory_limit", &tmp_int, &tmp_res, 50, Log, CHECK_NONE);
      if (tmp_res!=-1) generator2.setLimit(tmp_int); 
      // Change limit on how many variables in cut (default 50)
      getIntInPList(pvApiCtx, param_in_addr, "cglgomory_limitatroot", &tmp_int, &tmp_res, 50, Log, CHECK_NONE);
      if (tmp_res!=-1) generator2.setLimitAtRoot(tmp_int);
      // Change criterion on which variables to look at. All ones more than "away" away from integrality will be investigated (default 0.05)
      getDoubleInPList(pvApiCtx, param_in_addr, "cglgomory_away", &tmp_double, &tmp_res, 0.05, Log, CHECK_NONE);
      if (tmp_res!=-1) generator2.setAway(tmp_double);

      getIntInPList(pvApiCtx, param_in_addr, "cglgomory_priority", &tmp_int, &tmp_res, -1, Log, CHECK_NONE);
      model.addCutGenerator(&generator2,tmp_int,"Gomory");
    }

  // CglKnapsackCover generates "knapsack cover cuts". It looks for a series of different types of minimal covers. 
  // If a minimal cover is found, it lifts the associated minimal cover inequality and adds the lifted cut to the cut set.
  // CglKnapsackCover has a variety of methods for finding and lifting covers.
  // You can re-use the various cover-finding methods and cover-lifting methods to build your own variations of this classic cut.
  //
  // Cover-finding methods:
  //
  // * findGreedyCover: Try to generate a violated minimal cover greedily from fractional vars. 
  // * findJohnAndEllisCover: Try to generated a violated minimal cover using "John and Ellis" logic 
  //                          (i.e., my understanding of some of what John Forrest and Ellis Johnson used in OSL). 
  // * findPseudoJohnAndEllisCover: A variation on findJohnAndEllisCover. 
  // * findExactMostViolatedMinCover: Use an exact algorithm to find the most violated (minimal) cover. 
  // * findLPMostViolatedMinCover: Use an lp-relaxation to find the approximately most violated (minimal) cover. 
  //
  // Cover-lifting methods:
  //
  // * liftUpDownAndUncomplementAndAdd
  // * seqLiftAndUncomplementAndAdd: Sequence-dependent lifting 
  // * liftCoverCut: Sequence-independent lifting 
  //
  // Exact solvers for Knapsack Problems:
  // 
  // * exactSolveKnapsack: A goto-less implementation of the Horowitz-Sahni exact solution procedure for solving knapsack problem.
  //                       Reference: Martello and Toth, Knapsack Problems, Wiley, 1990, p30-31. 
  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cglknapsackcover");
  if (tmp_int!=-1)
    {
      CglKnapsackCover generator3;

      // Set limit on number in knapsack.
      getIntInPList(pvApiCtx, param_in_addr, "cglknapsackcover_maxinknapsack", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator3.setMaxInKnapsack(tmp_int);
      // - void switchOffExpensive()        - Switch off expensive cuts.
      getIntInPList(pvApiCtx, param_in_addr, "cglknapsackcover_switchoffexpensive", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) 
	if (tmp_int) generator3.switchOffExpensive();
      // - void switchOnExpensive()         - Switch on expensive cuts.
      getIntInPList(pvApiCtx, param_in_addr, "cglknapsackcover_switchonexpensive", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) 
	if (tmp_int) generator3.switchOnExpensive();

      getIntInPList(pvApiCtx, param_in_addr, "cglknapsackcover_priority", &tmp_int, &tmp_res, -1, Log, CHECK_NONE);
      model.addCutGenerator(&generator3,tmp_int,"Knapsack");
    }

  // Reduce-and-Split cuts are variants of Gomory cuts:
  // Starting from the current optimal tableau, linear combinations of the rows of the current optimal simplex tableau
  // are used for generating Gomory cuts. The choice of the linear combinations is driven by the objective of reducing the
  // coefficients of the non basic continuous variables in the resulting row.
  // Publication:
  // "Reduce-and-Split Cuts: Improving the Performance of Mixed Integer Gomory Cuts", Management Science 51 (2005)
  // by K. Anderson, G. Cornuejols, Yanjun Li.
  //
  // Warning: This generator currently works only with the Lp solvers Clp or Cplex9.0 or higher. 
  // It requires access to the optimal tableau and optimal basis inverse and makes assumptions on the way slack variables 
  // are added by the solver. The Osi implementations for Clp and Cplex verify these assumptions.
  //
  // Note that this generator might not be able to generate cuts for some solutions violating integrality constraints.
  //
  // Note also that when calling the generator, the solver interface si must contain an optimized problem and information related
  // to the optimal basis must be available through the OsiSolverInterface methods (si->optimalBasisIsAvailable() must return 'true'). 
  // It is also essential that the integrality of structural variable i can be obtained using si->isInteger(i).
  // If the first condition is not met, no cuts are generated. If the second condition is not met, no cuts or invalid cuts may be generated.
  //
  // Parameters:
  // Parameters of the generator are listed below. Modifying the default values for parameters other than the last five
  // might result in invalid cuts.
  //
  // * LUB: Value considered large for the absolute value of a lower or upper bound on a variable.
  //        Default value: 1000. See method setLUB().
  // * EPS: Precision of double computations.
  //        Default value: 1e-7. See method setEPS(). 
  // * EPS_COEFF: Precision for deciding if a coefficient of a generated cut is zero.
  //        Default value: 1e-8. See method setEPS_COEFF(). 
  // * EPS_COEFF_LUB: Precision for deciding if a coefficient of a generated cut is zero when the corresponding variable
  //                  has a lower or upper bound larger than LUB in absolute value.
  //        Default value: 1e-13. See method setEPS_COEFF_LUB(). 
  // * EPS_RELAX: Value used to relax slightly the right hand side of each generated cut.
  //        Default value: 1e-8. See method setEPS_RELAX(). 
  // * normIsZero: Norm of a vector is considered zero if smaller than this value.
  //        Default value: 1e-5. See method setNormIsZero(). 
  // * minReduc: Reduction is performed only if the norm of the vector is reduced by this fraction.
  //        Default value: 0.05. See method setMinReduc(). 
  // * limit: Generate cuts with at most this number of nonzero entries.
  //        Default value: 50. See method setLimit(). 
  // * away: Look only at basic integer variables whose current value is at least this value from being integer.
  //        Default value: 0.05. See method setAway(). 
  // * maxTab: Controls the number of rows selected for the generation.
  //        Default value: 1e7. See method setMaxTab(). 
  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cglredsplit");
  if (tmp_int!=-1)
    {
      CglRedSplit generator4;

      // Set limit, the maximum number of non zero coefficients in generated cut; Default: 50.
      getIntInPList(pvApiCtx, param_in_addr, "cglredsplit_limit", &tmp_int, &tmp_res, 50, Log, CHECK_NONE);
      if (tmp_res!=-1) generator4.setLimit(tmp_int);
      // Set away, the minimum distance from being integer used for selecting rows for cut generation;
      // all rows whose pivot variable should be integer but is more than away from integrality will be selected;
      // Default: 0.05.
      getDoubleInPList(pvApiCtx, param_in_addr, "cglredsplit_away", &tmp_double, &tmp_res, 0.05, Log, CHECK_NONE);
      if (tmp_res!=-1) generator4.setAway(tmp_double);
      // Set the value of LUB, value considered large for the absolute value of a lower or upper bound on a variable; Default: 1000.
      getDoubleInPList(pvApiCtx, param_in_addr, "cglredsplit_lub", &tmp_double, &tmp_res, 1000, Log, CHECK_NONE);
      if (tmp_res!=-1) generator4.setLUB(tmp_double);
      // Set the value of EPS, epsilon for double computations; Default: 1e-7.
      getDoubleInPList(pvApiCtx, param_in_addr, "cglredsplit_eps", &tmp_double, &tmp_res, 1e-7, Log, CHECK_NONE);
      if (tmp_res!=-1) generator4.setEPS(tmp_double);
      // Set the value of EPS_COEFF, epsilon for values of coefficients; Default: 1e-8.
      getDoubleInPList(pvApiCtx, param_in_addr, "cglredsplit_eps_coeff", &tmp_double, &tmp_res, 1e-8, Log, CHECK_NONE);
      if (tmp_res!=-1) generator4.setEPS_COEFF(tmp_double);
      // Set the value of EPS_COEFF_LUB, epsilon for values of coefficients for variables with absolute value of lower
      // or upper bound larger than LUB; Default: 1e-13.
      getDoubleInPList(pvApiCtx, param_in_addr, "cglredsplit_eps_coeff_lub", &tmp_double, &tmp_res, 1e-13, Log, CHECK_NONE);
      if (tmp_res!=-1) generator4.setEPS_COEFF_LUB(tmp_double);
      // Set the value of EPS_RELAX, value used for relaxing the right hand side of each generated cut; Default: 1e-8.
      getDoubleInPList(pvApiCtx, param_in_addr, "cglredsplit_eps_relax", &tmp_double, &tmp_res, 1e-8, Log, CHECK_NONE);
      if (tmp_res!=-1) generator4.setEPS_RELAX(tmp_double);
      // Set the value of normIsZero, the threshold for considering a norm to be 0; Default: 1e-5.
      getDoubleInPList(pvApiCtx, param_in_addr, "cglredsplit_normiszero", &tmp_double, &tmp_res, 1e-5, Log, CHECK_NONE);
      if (tmp_res!=-1) generator4.setNormIsZero(tmp_double);
      // Set the value of minReduc, threshold for relative norm improvement for performing a reduction; Default: 0.05.
      getDoubleInPList(pvApiCtx, param_in_addr, "cglredsplit_minreduc", &tmp_double, &tmp_res, 0.05, Log, CHECK_NONE);
      if (tmp_res!=-1) generator4.setMinReduc(tmp_double);
      // Set the maximum allowed value for (mTab * mTab * max(mTab, nTab)) where mTab is the number of rows 
      // used in the combinations and nTab is the number of continuous non basic variables.
      getDoubleInPList(pvApiCtx, param_in_addr, "cglredsplit_maxtab", &tmp_double, &tmp_res, 100, Log, CHECK_NONE);
      if (tmp_res!=-1) generator4.setMaxTab(tmp_double);

      getIntInPList(pvApiCtx, param_in_addr, "cglredsplit_priority", &tmp_int, &tmp_res, -1, Log, CHECK_NONE);
      model.addCutGenerator(&generator4,tmp_int,"RedSplit");
    }

  // Generates cuts of the form sum of a set of variables <= 1.
  // Note: This implementation is very fast, but design for set partitioning problems.
  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cglclique");
  if (tmp_int!=-1)
    {
      CglClique generator5;

      // possible choices for selecting the next node in the star clique search
      getIntInPList(pvApiCtx, param_in_addr, "cglclique_starcliquecandidatelengththreshold", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator5.setStarCliqueCandidateLengthThreshold(tmp_int);
      // possible choices for selecting the next node in the star clique search
      getIntInPList(pvApiCtx, param_in_addr, "cglclique_rowcliquecandidatelengththreshold", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator5.setRowCliqueCandidateLengthThreshold(tmp_int);
      // possible choices for selecting the next node in the star clique search
      getIntInPList(pvApiCtx, param_in_addr, "cglclique_starcliquereport", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator5.setStarCliqueReport((bool)tmp_int);
      // possible choices for selecting the next node in the star clique search
      getIntInPList(pvApiCtx, param_in_addr, "cglclique_rowcliquereport", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator5.setRowCliqueReport((bool)tmp_int);
      // possible choices for selecting the next node in the star clique search
      getIntInPList(pvApiCtx, param_in_addr, "cglclique_dostarclique", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator5.setDoStarClique((bool)tmp_int);
      // possible choices for selecting the next node in the star clique search
      getIntInPList(pvApiCtx, param_in_addr, "cglclique_dorowclique", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator5.setDoRowClique((bool)tmp_int);
      // possible choices for selecting the next node in the star clique search
      getDoubleInPList(pvApiCtx, param_in_addr, "cglclique_minviolation", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator5.setMinViolation(tmp_double);

      getIntInPList(pvApiCtx, param_in_addr, "cglclique_priority", &tmp_int, &tmp_res, -1, Log, CHECK_NONE);
      model.addCutGenerator(&generator5,tmp_int,"Clique");
    }

  // Generates mixed integer rounding cuts.
  //
  // Another mixed integer rounding cut generator, CglMixedIntegerRounding2, is very similar. 
  // CglMixedIntegerRounding uses CoinPackedVector whereas CglMixedIntegerRounding2 uses CoinIndexedVector. 
  // Depending on the size of the problem, one generator might be faster than the other. 
  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cglmixedintegerrounding");
  if (tmp_int!=-1)
    {
      CglMixedIntegerRounding mixedGen1;

      getIntInPList(pvApiCtx, param_in_addr, "cglmixedintegerrounding_maxaggr", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) mixedGen1.setMAXAGGR_(tmp_int);
      getIntInPList(pvApiCtx, param_in_addr, "cglmixedintegerrounding_multiply", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) mixedGen1.setMULTIPLY_((bool)tmp_int);
      getIntInPList(pvApiCtx, param_in_addr, "cglmixedintegerrounding_criterion", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) mixedGen1.setCRITERION_(tmp_int);
      getIntInPList(pvApiCtx, param_in_addr, "cglmixedintegerrounding_dopreproc", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) mixedGen1.setDoPreproc(tmp_int);

      getIntInPList(pvApiCtx, param_in_addr, "cglmixedintegerrounding_priority", &tmp_int, &tmp_res, -1, Log, CHECK_NONE);
      model.addCutGenerator(&mixedGen1,tmp_int,"MixedIntegerRounding");
    }

  // Generates mixed integer rounding cuts.
  // Another mixed integer rounding cut generator, CglMixedIntegerRounding, is very similar. 
  // CglMixedIntegerRounding uses CoinPackedVector whereas CglMixedIntegerRounding2 uses CoinIndexedVector. 
  // Depending on the size of the problem, one generator might be faster than the other.
  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cglmixedintegerrounding2");
  if (tmp_int!=-1)
    {
      CglMixedIntegerRounding2 mixedGen2;

      getIntInPList(pvApiCtx, param_in_addr, "cglmixedintegerrounding2_maxaggr", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) mixedGen2.setMAXAGGR_(tmp_int);
      getIntInPList(pvApiCtx, param_in_addr, "cglmixedintegerrounding2_multiply", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) mixedGen2.setMULTIPLY_((bool)tmp_int);
      getIntInPList(pvApiCtx, param_in_addr, "cglmixedintegerrounding2_criterion", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) mixedGen2.setCRITERION_(tmp_int);
      getIntInPList(pvApiCtx, param_in_addr, "cglmixedintegerrounding2_dopreproc", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) mixedGen2.setDoPreproc(tmp_int);

      getIntInPList(pvApiCtx, param_in_addr, "cglmixedintegerrounding2_priority", &tmp_int, &tmp_res, -1, Log, CHECK_NONE);
      model.addCutGenerator(&mixedGen2,tmp_int,"MixedIntegerRounding2");
    }

  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cglflowcover");
  if (tmp_int!=-1)
    {
      CglFlowCover flowGen;

      // Functions to query and set the number of cuts have been generated.
      getIntInPList(pvApiCtx, param_in_addr, "cglflowcover_numflowcuts", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) flowGen.setNumFlowCuts(tmp_int);

      getIntInPList(pvApiCtx, param_in_addr, "cglflowcover_priority", &tmp_int, &tmp_res, -1, Log, CHECK_NONE);
      model.addCutGenerator(&flowGen,tmp_int,"FlowCover");
    }

  // Generates odd holes cuts. 
  // This looks at all rows of type sum x(i) <= 1 (or == 1) (with x binary) and sees if there is an odd cycle cut.
  // See Grotschel, Lovasz and Schrijver (1988) for the method. 
  // This is then lifted by using the corresponding Chvatal cut i.e. by summing up all rows in the cycle.
  // The right hand side will be odd and all odd coefficients can be reduced by one. The constraint is
  //
  // sum even(j)*x(j) <= odd
  //
  // which can be replaced by
  //
  // sum (even(j)/2)*x(j) <= (odd-1.0)/2.
  //
  // A similar cut can be generated for sum x(i) >= 1. 
  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cgloddhole");
  if (tmp_int!=-1)
    {
      CglOddHole generator6;

      getDoubleInPList(pvApiCtx, param_in_addr, "cgloddhole_minimumviolation", &tmp_double, &tmp_res, 0.005, Log, CHECK_NONE);
      if (tmp_res!=-1) generator6.setMinimumViolation(tmp_double);
      getDoubleInPList(pvApiCtx, param_in_addr, "cgloddhole_minimumviolationper", &tmp_double, &tmp_res, 0.00002, Log, CHECK_NONE);
      if (tmp_res!=-1) generator6.setMinimumViolationPer(tmp_double);
      // - void setMaximumEntries(int value)         - Minimum violation.
      getIntInPList(pvApiCtx, param_in_addr, "cgloddhole_maximumentries", &tmp_int, &tmp_res, 200, Log, CHECK_NONE);
      if (tmp_res!=-1) generator6.setMaximumEntries(tmp_int);

      getIntInPList(pvApiCtx, param_in_addr, "cgloddhole_priority", &tmp_int, &tmp_res, -1, Log, CHECK_NONE);
      model.addCutGenerator(&generator6,tmp_int,"OddHold");
    }

  // Generates two step mixed integer rounding cuts either from the tableau rows or from the formulation rows. 
  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cgltwomir");
  if (tmp_int!=-1)
    {
      int tmp_int_1, tmp_int_2, tmp_int_3, tmp_int_4;

      CglTwomir generator7;

      // Change criterion on which scalings to use (default = 1,1,1,1)
      // - void setMirScale(int tmin, int tmax)
      getIntInPList(pvApiCtx, param_in_addr, "cgltwomir_mirscal_tmin", &tmp_int_1, &tmp_res, 1, Log, CHECK_NONE);
      getIntInPList(pvApiCtx, param_in_addr, "cgltwomir_mirscal_tmax", &tmp_int_2, &tmp_res, tmp_int_1, Log, CHECK_NONE);
      if (tmp_res!=-1) generator7.setMirScale(tmp_int_1, tmp_int_2);
      // - void setTwomirScale(int qmin, int qmax)
      getIntInPList(pvApiCtx, param_in_addr, "cgltwomir_mirscal_qmin", &tmp_int_1, &tmp_res, 1, Log, CHECK_NONE);
      getIntInPList(pvApiCtx, param_in_addr, "cgltwomir_mirscal_qmax", &tmp_int_2, &tmp_res, tmp_int_1, Log, CHECK_NONE);
      if (tmp_res!=-1) generator7.setTwomirScale(tmp_int_1, tmp_int_2);
      // - void setAMax(int a)
      getIntInPList(pvApiCtx, param_in_addr, "cgltwomir_amax", &tmp_int_1, &tmp_res, 1, Log, CHECK_NONE);
      if (tmp_res!=-1) generator7.setAMax(tmp_int_1);
      // - void setMaxElements(int n)
      getIntInPList(pvApiCtx, param_in_addr, "cgltwomir_maxelements", &tmp_int_1, &tmp_res, 200, Log, CHECK_NONE);
      if (tmp_res!=-1) generator7.setMaxElements(tmp_int_1);
      // - void setCutTypes(bool mir, bool twomir, bool tab, bool form)
      getIntInPList(pvApiCtx, param_in_addr, "cgltwomir_cuttype_mir",    &tmp_int_1, &tmp_res, 0, Log, CHECK_NONE);
      getIntInPList(pvApiCtx, param_in_addr, "cgltwomir_cuttype_twomir", &tmp_int_2, &tmp_res, 0, Log, CHECK_NONE);
      getIntInPList(pvApiCtx, param_in_addr, "cgltwomir_cuttype_tab",    &tmp_int_3, &tmp_res, 0, Log, CHECK_NONE);
      getIntInPList(pvApiCtx, param_in_addr, "cgltwomir_cuttype_form",   &tmp_int_4, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator7.setCutTypes((bool)tmp_int_1, (bool)tmp_int_2, (bool)tmp_int_3, (bool)tmp_int_4);
      // - void setFormulationRows(int n)
      getIntInPList(pvApiCtx, param_in_addr, "cgltwomir_formulationrows", &tmp_int, &tmp_res, 10, Log, CHECK_NONE);
      if (tmp_res!=-1) generator7.setFormulationRows(tmp_int);

      getIntInPList(pvApiCtx, param_in_addr, "cgltwomir_priority", &tmp_int, &tmp_res, -1, Log, CHECK_NONE);
      model.addCutGenerator(&generator7,tmp_int,"Twomir");
    }

  // Generate cuts for the situation where a set of general integer variables must have distinct values. 
  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cglalldifferent");
  if (tmp_int!=-1)
    {
      CglAllDifferent generator8;

      getIntInPList(pvApiCtx, param_in_addr, "cglalldifferent_maxl", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator8.setMaxLook(tmp_int);

      getIntInPList(pvApiCtx, param_in_addr, "cglalldifferent_priority", &tmp_int, &tmp_res, -1, Log, CHECK_NONE);
      model.addCutGenerator(&generator8,tmp_int,"AllDifferent");
    }

  // Fix variables and find duplicate/dominated rows. 
  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cglduplicaterow");
  if (tmp_int!=-1)
    {
      CglDuplicateRow generator9;

      // We only check for dominated amongst groups of columns whose size <= this
      getIntInPList(pvApiCtx, param_in_addr, "cglduplicaterow_maximumrhs", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator9.setMaximumRhs(tmp_int);
      // We only check for dominated amongst groups of columns whose size <= this
      getIntInPList(pvApiCtx, param_in_addr, "cglduplicaterow_maximumdominated", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator9.setMaximumDominated(tmp_int);
      // Set mode. 
      getIntInPList(pvApiCtx, param_in_addr, "cglduplicaterow_mode", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator9.setMode(tmp_int);

      getIntInPList(pvApiCtx, param_in_addr, "cglduplicaterow_priority", &tmp_int, &tmp_res, -1, Log, CHECK_NONE);
      model.addCutGenerator(&generator9, tmp_int, "DuplicateRow");
    }

  // CglLiftAndProject is a very limited implementation of "Lift-and-Project" cuts.
  //
  // The main purpose of this cut generator is to encourage someone to improve it (or better yet, to get them
  // to CONTRIBUTE their own, better, implementation of lift-and-project cuts.
  // CglLiftAndProject is NOT a general purpose cut generators. It makes strong assumptions (i.e., the User knows
  // if the requirments for this cg are satisfied). This implementation uses Normalization 1, and does not include lifting.
  //
  // Assumes the mixed 0-1 problem
  //
  // min {cx: A x >= b}
  //
  // is in canonical form with all bounds, including x_t>=0, -x_t>=-1 for x_t binary, explicitly stated in the constraint matrix.
  // Given canonical problem and the lp-relaxation solution, x, CglLiftAndProject attempts to construct a cut 
  // for every x_j such that 0 < x_j < 1. 
  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cglliftandproject");
  if (tmp_int!=-1)
    {
      CglLiftAndProject generator10;

      // Set the normalization : Either beta=+1 or beta=-1.
      getIntInPList(pvApiCtx, param_in_addr, "cglliftandproject_beta", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
      if (tmp_res!=-1) generator10.setBeta(tmp_int);

      getIntInPList(pvApiCtx, param_in_addr, "cglliftandproject_priority", &tmp_int, &tmp_res, -1, Log, CHECK_NONE);
      model.addCutGenerator(&generator10, tmp_int, "LiftAndProject");
    }
      
  // CglSimpleRounding generates "simple rounding cuts" (described below). 
  // The main purpose of this cut is to provide a simple example of how to write a cut generator to would-be developers.
  //
  // CglSimpleRounding generators simple rounding cuts via the following method:
  //
  // * For each contraint, attempt to derive a <= inequality in all integer variables by netting out any continuous variables. 
  // * Divide the resulting integer inequality through by the greatest common denomimator (gcd) of the lhs coefficients. 
  // * Round down the rhs. 
  //
  // See Nemhauser and Wolsey, Integer and Combinatorial Optimization, 1988, p. 211.
  // Notes:
  //
  // 1. Meant as an example. Not expected to be computationally useful.
  // 2. Warning: Use with careful attention to data precision. 
  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cglsimplerounding");
  if (tmp_int!=-1)
    {
      CglSimpleRounding generator11;

      getIntInPList(pvApiCtx, param_in_addr, "cglsimplerounding_priority", &tmp_int, &tmp_res, -1, Log, CHECK_NONE);
      model.addCutGenerator(&generator11, tmp_int, "SimpleRounding");
    }

  // This is an implementation of a separation algorithm for Residual Capacity Inequalities. 
  // They have been introduced in T. L. Magnanti, P. Mirchandani, and R. Vachani, "The convex hull of two core capacitated
  // network design problems," Math. Programming, 60 (1993), pp. 233-250. 
  // The separation algorithm was given in A. Atamturk and D. Rajan, "On splittable and unsplittable flow capacitated network 
  // design arc-set polyhedra," Math. Program., 92 (2002), pp. 315-333.
  //
  // These inequalities are particularly useful for Network Design and Capacity Planning models.
  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cglresidualcapacity");
  if (tmp_int!=-1)
    {
      CglResidualCapacity generator12;
	  
      getDoubleInPList(pvApiCtx, param_in_addr, "cglresidualcapacity_epsilon", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator12.setEpsilon(tmp_double);
      getDoubleInPList(pvApiCtx, param_in_addr, "cglresidualcapacity_tolerance", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator12.setTolerance(tmp_double);
      getIntInPList(pvApiCtx, param_in_addr, "cglresidualcapacity_dopreproc", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) generator12.setDoPreproc(tmp_int);

      getIntInPList(pvApiCtx, param_in_addr, "cglresidualcapacity_priority", &tmp_int, &tmp_res, -1, Log, CHECK_NONE);
      model.addCutGenerator(&generator12, tmp_int, "ResidualCapacity");
    }

  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cglimplication");
  if (tmp_int!=-1)
    {
      CglImplication generator13;

      getIntInPList(pvApiCtx, param_in_addr, "cglimplication_priority", &tmp_int, &tmp_res, -1, Log, CHECK_NONE);
      model.addCutGenerator(&generator13, tmp_int, "Implication");
    }

  // This procedure implements different variants of the lift-and-project procedure executed in the LP simplex tableau
  // described by Balas and Perregaard.
  // Acknoweldgements:
  // The code has been developped in a joint work with Egon Balas.
  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cgllandp");
  if (tmp_int!=-1)
    {
      CglLandP generator14;

      getIntInPList(pvApiCtx, param_in_addr, "cgllandp_priority", &tmp_int, &tmp_res, -1, Log, CHECK_NONE);
      model.addCutGenerator(&generator14, tmp_int, "LandP");
    }

  /////////////////////////
  // Set some heuristics //
  /////////////////////////

  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cbcrounding");
  if (tmp_int!=-1)
    {
      CbcRounding heuristic1(model);
      getIntInPList(pvApiCtx, param_in_addr, "cbcrounding_seed", &tmp_int, &tmp_res, 12345, Log, CHECK_NONE);
      if (tmp_res!=-1) heuristic1.setSeed(tmp_int);
      model.addHeuristic(&heuristic1,"CbcRounding");
    }

  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cbcheuristicdivecoefficient");
  if (tmp_int!=-1)
    {
      CbcHeuristicDiveCoefficient heuristic3(model);
      model.addHeuristic(&heuristic3,"CbcHeuristicDiveCoefficient");
    }

  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cbcheuristicdivefractional");
  if (tmp_int!=-1)
    {
      CbcHeuristicDiveFractional heuristic4(model);
      model.addHeuristic(&heuristic4,"CbcHeuristicFractional");
    }

  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cbcheuristicdiveguided");
  if (tmp_int!=-1)
    {
      CbcHeuristicDiveGuided heuristic5(model);
      model.addHeuristic(&heuristic5,"CbcHeuristicDiveGuided");
    }

  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cbcheuristicdivevectorlength");
  if (tmp_int!=-1)
    {
      CbcHeuristicDiveVectorLength heuristic6(model);
      model.addHeuristic(&heuristic6,"CbcHeuristicDiveVectorLength");
    }

  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cbcheuristicfpump");
  if (tmp_int!=-1)
    {
      CbcHeuristicFPump heuristic8(model);
      // Set maximum Time (default off) - also sets starttime to current.
      getDoubleInPList(pvApiCtx, param_in_addr, "cbcheuristicfpump_maximumtime", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) heuristic8.setMaximumTime(tmp_double);
      // Set fake cutoff (default COIN_DBL_MAX == off).
      getDoubleInPList(pvApiCtx, param_in_addr, "cbcheuristicfpump_fakecutoff", &tmp_double, &tmp_res, COIN_DBL_MAX, Log, CHECK_NONE);
      if (tmp_res!=-1) heuristic8.setFakeCutoff(tmp_double);
      // Set absolute increment (default 0.0 == off).
      getDoubleInPList(pvApiCtx, param_in_addr, "cbcheuristicfpump_absoluteincrement", &tmp_double, &tmp_res, 0.0, Log, CHECK_NONE);
      if (tmp_res!=-1) heuristic8.setAbsoluteIncrement(tmp_double);
      // Set relative increment (default 0.0 == off).
      getDoubleInPList(pvApiCtx, param_in_addr, "cbcheuristicfpump_relativeincrement", &tmp_double, &tmp_res, 0.0, Log, CHECK_NONE);
      if (tmp_res!=-1) heuristic8.setRelativeIncrement(tmp_double);
      // Set default rounding (default 0.5).
      getDoubleInPList(pvApiCtx, param_in_addr, "cbcheuristicfpump_defaultrounding", &tmp_double, &tmp_res, 0.5, Log, CHECK_NONE);
      if (tmp_res!=-1) heuristic8.setDefaultRounding(tmp_double);
      // Set initial weight (default 0.0 == off).
      getDoubleInPList(pvApiCtx, param_in_addr, "cbcheuristicfpump_initialweight", &tmp_double, &tmp_res, 0.0, Log, CHECK_NONE);
      if (tmp_res!=-1) heuristic8.setInitialWeight(tmp_double);
      // Set weight factor (default 0.1).
      getDoubleInPList(pvApiCtx, param_in_addr, "cbcheuristicfpump_weightfactor", &tmp_double, &tmp_res, 0.1, Log, CHECK_NONE);
      if (tmp_res!=-1) heuristic8.setWeightFactor(tmp_double);
      // Set threshold cost for using original cost - even on continuous (default infinity).
      getDoubleInPList(pvApiCtx, param_in_addr, "cbcheuristicfpump_artificialcost", &tmp_double, &tmp_res, COIN_DBL_MAX, Log, CHECK_NONE);
      if (tmp_res!=-1) heuristic8.setArtificialCost(tmp_double);
      // Set maximum passes (default 100).
      getIntInPList(pvApiCtx, param_in_addr, "cbcheuristicfpump_maximumpasses", &tmp_int, &tmp_res, 100, Log, CHECK_NONE);
      if (tmp_res!=-1) heuristic8.setMaximumPasses(tmp_int);
      // Set maximum retries (default 1).
      getIntInPList(pvApiCtx, param_in_addr, "cbcheuristicfpump_maximumretries", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
      if (tmp_res!=-1) heuristic8.setMaximumRetries(tmp_int);
      // Set use of multiple solutions and solves:
      // - 0 - do not reuse solves, do not accumulate integer solutions for local search 
      // - 1 - do not reuse solves, accumulate integer solutions for local search 
      // - 2 - reuse solves, do not accumulate integer solutions for local search 
      // - 3 - reuse solves, accumulate integer solutions for local search.
      getIntInPList(pvApiCtx, param_in_addr, "cbcheuristicfpump_accumulate", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) heuristic8.setAccumulate(tmp_int);
      //   - void setFixOnReducedCosts(int value)     - Set whether to fix variables on known solution:
      //                                                - 0 - do not fix
      //                                                - 1 - fix integers on reduced costs
      //                                                - 2 - fix integers on reduced costs but only on entry.
      getIntInPList(pvApiCtx, param_in_addr, "cbcheuristicfpump_fixonreducedcosts", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) heuristic8.setFixOnReducedCosts(tmp_int);

      model.addHeuristic(&heuristic8,"CbcHeuristicFPump");
    }

  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cbcheuristicgreedycover");
  if (tmp_int!=-1)
    {
      CbcHeuristicGreedyCover heuristic9(model);

      getIntInPList(pvApiCtx, param_in_addr, "cbcheuristicgreedycover_algorithm", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) heuristic9.setAlgorithm(tmp_int);
      getIntInPList(pvApiCtx, param_in_addr, "cbcheuristicgreedycover_numbertimes", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) heuristic9.setNumberTimes(tmp_int);
	  
      model.addHeuristic(&heuristic9,"CbcHeuristicGreedyCover");
    }

  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cbcheuristicgreedyequality");
  if (tmp_int!=-1)
    {
      CbcHeuristicGreedyEquality heuristic10(model);

      getIntInPList(pvApiCtx, param_in_addr, "cbcheuristicgreedyequality_algorithm", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) heuristic10.setAlgorithm(tmp_int);
      getDoubleInPList(pvApiCtx, param_in_addr, "cbcheuristicgreedyequality_fraction", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) heuristic10.setFraction(tmp_double);
      getIntInPList(pvApiCtx, param_in_addr, "cbcheuristicgreedyequality_numbertimes", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) heuristic10.setNumberTimes(tmp_int);

      model.addHeuristic(&heuristic10,"CbcHeuristicGreedyEquality");
    }

  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cbcheuristiclocal");
  if (tmp_int!=-1)
    {
      CbcHeuristicLocal heuristic11(model);

      getIntInPList(pvApiCtx, param_in_addr, "cbcheuristiclocal_searchtype", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) heuristic11.setSearchType(tmp_int);

      model.addHeuristic(&heuristic11,"CbcHeuristicLocal");
    }

  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cbcheuristicpartial");
  if (tmp_int!=-1)
    {
      CbcHeuristicPartial heuristic12(model);

      getIntInPList(pvApiCtx, param_in_addr, "cbcheuristicpartial_fixpriority", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) heuristic12.setFixPriority(tmp_int);

      model.addHeuristic(&heuristic12,"CbcHeuristicPartial");
    }

  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cbcheuristicrens");
  if (tmp_int!=-1)
    {
      CbcHeuristicRENS heuristic13(model);
      model.addHeuristic(&heuristic13,"CbcHeuristicRENS");
    }

  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cbcheuristicrins");
  if (tmp_int!=-1)
    {
      CbcHeuristicRINS heuristic14(model);

      getIntInPList(pvApiCtx, param_in_addr, "cbcheuristicrins_howoften", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) heuristic14.setHowOften(tmp_int);
      getDoubleInPList(pvApiCtx, param_in_addr, "cbcheuristicrins_decayfactor", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) heuristic14.setDecayFactor(tmp_double);

      model.addHeuristic(&heuristic14,"CbcHeuristicRINS");
    }

  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cbcserendipity");
  if (tmp_int!=-1)
    {
      CbcSerendipity heuristic15(model);
      model.addHeuristic(&heuristic15,"CbcSerendipity");
    }

  //////////////////////////////////
  // Set Branching decision class //
  //////////////////////////////////

  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cbcbranchdefaultdecision");
  if (tmp_int!=-1)
    {
      CbcBranchDefaultDecision branch1;

      getDoubleInPList(pvApiCtx, param_in_addr, "cbcbranchdefaultdecision_bestcriterion", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) branch1.setBestCriterion(tmp_double);

      model.setBranchingMethod(&branch1);
    }

  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cbcbranchdynamicdecision");
  if (tmp_int!=-1)
    {
      CbcBranchDynamicDecision branch2;

      getDoubleInPList(pvApiCtx, param_in_addr, "cbcbranchdynamicdecision_bestcriterion", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
      if (tmp_res!=-1) branch2.setBestCriterion(tmp_double);

      model.setBranchingMethod(&branch2);
    }

  ////////////////////////
  // Set Strategy class //
  ////////////////////////

  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cbcstrategydefaultsubtree");
  if (tmp_int!=-1)
    {
      int tmp_int_1, tmp_int_2, tmp_int_3, tmp_int_4;

      getIntInPList(pvApiCtx, param_in_addr, "cbcstrategydefaultsubtree_cutsonlyatroot",    &tmp_int_1, &tmp_res, 0, Log, CHECK_NONE);
      getIntInPList(pvApiCtx, param_in_addr, "cbcstrategydefaultsubtree_numberstrong",      &tmp_int_2, &tmp_res, 5, Log, CHECK_NONE);
      getIntInPList(pvApiCtx, param_in_addr, "cbcstrategydefaultsubtree_numberbeforetrust", &tmp_int_3, &tmp_res, 0, Log, CHECK_NONE);
      getIntInPList(pvApiCtx, param_in_addr, "cbcstrategydefaultsubtree_printlevel",        &tmp_int_4, &tmp_res, 0, Log, CHECK_NONE);

      CbcStrategyDefaultSubTree strategy2(NULL, (bool)tmp_int_1, tmp_int_2, tmp_int_3, tmp_int_4);
      model.setStrategy(strategy2);
    }

  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cbcstrategydefault");
  if (tmp_int!=-1)
    {
      int tmp_int_1, tmp_int_2, tmp_int_3, tmp_int_4;

      getIntInPList(pvApiCtx, param_in_addr, "cbcstrategydefault_cutsonlyatroot",        &tmp_int_1, &tmp_res, 0, Log, CHECK_NONE);
      getIntInPList(pvApiCtx, param_in_addr, "cbcstrategydefault_numberstrong",          &tmp_int_2, &tmp_res, 0, Log, CHECK_NONE);
      getIntInPList(pvApiCtx, param_in_addr, "cbcstrategydefault_numberbeforetrust",     &tmp_int_3, &tmp_res, 0, Log, CHECK_NONE);
      getIntInPList(pvApiCtx, param_in_addr, "cbcstrategydefault_printlevel",            &tmp_int_4, &tmp_res, 0, Log, CHECK_NONE);

      CbcStrategyDefault strategy1((bool)tmp_int_1, tmp_int_2, tmp_int_3, tmp_int_4);
      model.setStrategy(strategy1);
    }

  tmp_int = hasPartialLabelInPList(pvApiCtx, param_in_addr, "cbcstrategynull");
  if (tmp_int!=-1)
    {
      CbcStrategyNull strategy3;
      model.setStrategy(strategy3);
    }

  //////////////////////////////////////
  // Some Cbc parameters are settable //
  //////////////////////////////////////

  bool tmp_bool;

  // Set cutoff bound on the objective function.
  getDoubleInPList(pvApiCtx, param_in_addr, "cbc_cutoff", &tmp_double, &tmp_res, 1e100, Log, CHECK_NONE);
  if (tmp_res!=-1) model.setCutoff(tmp_double);

  // Set the maximum node limit .
  getIntInPList(pvApiCtx, param_in_addr, "cbc_maximumnodes", &tmp_int, &tmp_res, 2147483647, Log, CHECK_NONE);
  if (tmp_res!=-1) model.setMaximumNodes(tmp_int);

  // Set the maximum number of solutions desired.
  getIntInPList(pvApiCtx, param_in_addr, "cbc_maximumsolutions", &tmp_int, &tmp_res, 9999999, Log, CHECK_NONE);
  if (tmp_res!=-1) model.setMaximumSolutions(tmp_int);

  // Set the printing mode.
  // Adjusts printout 1 does different node message with number unsatisfied on last branch. 
  getIntInPList(pvApiCtx, param_in_addr, "cbc_printingmode", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) model.setPrintingMode(tmp_int);

  // Set the maximum number of seconds desired.
  getDoubleInPList(pvApiCtx, param_in_addr, "cbc_maximumseconds", &tmp_double, &tmp_res, 1e100, Log, CHECK_NONE);
  if (tmp_res!=-1) tmp_bool = model.setMaximumSeconds(tmp_double);

  // Set the integrality tolerance .
  getDoubleInPList(pvApiCtx, param_in_addr, "cbc_integertolerance", &tmp_double, &tmp_res, 1e-6, Log, CHECK_NONE);
  if (tmp_res!=-1) tmp_bool = model.setIntegerTolerance(tmp_double);

  // Set the weight per integer infeasibility .
  getDoubleInPList(pvApiCtx, param_in_addr, "cbc_infeasibilityweight", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) tmp_bool = model.setInfeasibilityWeight(tmp_double);

  // Set the allowable gap between the best known solution and the best possible solution.
  getDoubleInPList(pvApiCtx, param_in_addr, "cbc_allowablegap", &tmp_double, &tmp_res, 1e-10, Log, CHECK_NONE);
  if (tmp_res!=-1) tmp_bool = model.setAllowableGap(tmp_double);

  // Set the fraction allowable gap between the best known solution and the best possible solution.
  getDoubleInPList(pvApiCtx, param_in_addr, "cbc_allowablefractiongap", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) tmp_bool = model.setAllowableFractionGap(tmp_double);

  // Set the percentage allowable gap between the best known solution and the best possible solution.
  getDoubleInPList(pvApiCtx, param_in_addr, "cbc_allowablepercentagegap", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) tmp_bool = model.setAllowablePercentageGap(tmp_double);

  // Set the CbcModel::CbcCutoffIncrement desired.
  getDoubleInPList(pvApiCtx, param_in_addr, "cbc_cutoffincrement", &tmp_double, &tmp_res, 1e-5, Log, CHECK_NONE);
  if (tmp_res!=-1) tmp_bool = model.setCutoffIncrement(tmp_double);

  // Set the minimum drop to continue cuts.
  getDoubleInPList(pvApiCtx, param_in_addr, "cbc_minimumdrop", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) model.setMinimumDrop(tmp_double);

  // Set the maximum number of cut passes at root node (default 20) 
  // Minimum drop can also be used for fine tuning.
  getDoubleInPList(pvApiCtx, param_in_addr, "cbc_maxcutpassesatroot", &tmp_double, &tmp_res, 100, Log, CHECK_NONE);
  if (tmp_res!=-1) model.setMaximumCutPassesAtRoot(tmp_int);

  // Set the maximum number of cut passes at other nodes (default 10)
  // Minimum drop can also be used for fine tuning.
  getIntInPList(pvApiCtx, param_in_addr, "cbc_maxcutpasses", &tmp_int, &tmp_res, 5, Log, CHECK_NONE);
  if (tmp_res!=-1) model.setMaximumCutPasses(tmp_int);

  // Set the maximum number of candidates to be evaluated for strong branching.
  // Maximum number of candidates to consider for strong branching.
  // To disable strong branching, set this to 0.
  getIntInPList(pvApiCtx, param_in_addr, "cbc_numberstrong", &tmp_int, &tmp_res, 5, Log, CHECK_NONE);
  if (tmp_res!=-1) model.setNumberStrong(tmp_int);

  // Set global preferred way to branch -1 down, +1 up, 0 no preference.
  getIntInPList(pvApiCtx, param_in_addr, "cbc_preferredway", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) model.setPreferredWay(tmp_int);

  // NOT AVAILABLE ANYMORE UNDER CBC-2.3
  // // Set size of mini - tree.
  // GET_PARAM_INT("cbc_sizeminitree",  tmp_int, 2, tmp_res);
  // if (tmp_res!=-1) model.setSizeMiniTree(tmp_int);

  // Set the number of branches before pseudo costs believed in dynamic strong branching.
  getIntInPList(pvApiCtx, param_in_addr, "cbc_numberbeforetrust", &tmp_int, &tmp_res, 10, Log, CHECK_NONE);
  if (tmp_res!=-1) model.setNumberBeforeTrust(tmp_int);

  // Set the number of variables for which to compute penalties in dynamic strong branching.
  getIntInPList(pvApiCtx, param_in_addr, "cbc_numberpenalties", &tmp_int, &tmp_res, 20, Log, CHECK_NONE);
  if (tmp_res!=-1) model.setNumberPenalties(tmp_int);

  // Number of analyze iterations to do.
  getIntInPList(pvApiCtx, param_in_addr, "cbc_numberanalyzeiterations", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) model.setNumberAnalyzeIterations(tmp_int);

  // Set scale factor to make penalties match strong.
  getDoubleInPList(pvApiCtx, param_in_addr, "cbc_penaltyscalefactor", &tmp_double, &tmp_res, 3.0, Log, CHECK_NONE);
  if (tmp_res!=-1) model.setPenaltyScaleFactor(tmp_double);

  // Problem type as set by user or found by analysis.
  // This will be extended:
  // - 0 - not known
  // - 1 - Set partitioning <= 
  // - 2 - Set partitioning == 
  // - 3 - Set covering 
  // - 4 - all +- 1 or all +1 and odd
  getIntInPList(pvApiCtx, param_in_addr, "cbc_problemtype", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) model.setProblemType(tmp_int);

  // Set how often to scan global cuts.
  getIntInPList(pvApiCtx, param_in_addr, "cbc_howoftenglobalscan", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) model.setHowOftenGlobalScan(tmp_int);

  // Set the print frequency.
  getIntInPList(pvApiCtx, param_in_addr, "cbc_printfrequency", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) model.setPrintFrequency(tmp_int);

#ifdef THREAD_SUPPORT
  // Set number of threads
  getIntInPList(pvApiCtx, param_in_addr, "cbc_numberthreads", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) model.setNumberThreads(tmp_int);
#endif

  // writemps option
  getStringInPList(pvApiCtx, param_in_addr, "writemps", &tmp_char, &tmp_res, "test.mps", Log, CHECK_NONE);

  if (tmp_res!=-1) 
    {
      solver4->writeMps(tmp_char);
      if (loglevel) 
	{
	  sciprint("%s: writing %s mps file\n",fname, tmp_char);
	}
      FREE(tmp_char);
    }

  ////////////////////////
  // Do complete search //
  ////////////////////////

  // The method assumes that initialSolve() has been called to solve the LP relaxation. 
  // It processes the root node, then proceeds to explore the branch & cut search tree.
  // The search ends when the tree is exhausted or one of several execution limits is reached.
  // If doStatistics is 1 summary statistics are printed
  // if 2 then also the path to best solution (if found by branching)
  // if 3 then also one line per node
  getIntInPList(pvApiCtx, param_in_addr, "cbc_dobranchandbound", &tmp_int, &tmp_res, 0, Log, CHECK_NONE); // Perform branch and bound ?

  if (tmp_res!=-1)
    {
      if (loglevel) 
	{
	  sciprint("%s: do branch and bound\n", fname);
	}

      printer->setLogLevel(loglevel);
      model.passInMessageHandler(printer);

      TRYCATCH(model.branchAndBound(tmp_int))

      status = model.status();
    }
  else
    {
      // Then we get the status from the osi solver
      status  = (int)pow(2.0,1.0)*model.solver()->isAbandoned();                   // Are there numerical difficulties?
      status += (int)pow(2.0,2.0)*model.solver()->isProvenOptimal();               // Is optimality proven?
      status += (int)pow(2.0,3.0)*model.solver()->isProvenPrimalInfeasible();      // Is primal infeasiblity proven?
      status += (int)pow(2.0,4.0)*model.solver()->isProvenDualInfeasible();        // Is dual infeasiblity proven?
      status += (int)pow(2.0,5.0)*model.solver()->isPrimalObjectiveLimitReached(); // Is the given primal objective limit reached?
      status += (int)pow(2.0,6.0)*model.solver()->isDualObjectiveLimitReached();   // Is the given dual objective limit reached?
      status += (int)pow(2.0,7.0)*model.solver()->isIterationLimitReached();       // status = solver2->status();
    }

#ifdef DEBUG
  sciprint("compute outputs\n");
#endif

  // Methods returning info on how the solution process terminated

  double cbc_status = 0, node_count, iteration_count, secondary_status, initsolve_status = 0, res_time;

  cbc_status  = (int)(pow(2.0,0.0)*model.isAbandoned());            // Are there a numerical difficulties?
  cbc_status += (int)(pow(2.0,1.0)*model.isProvenOptimal());        // Is optimality proven?
  cbc_status += (int)(pow(2.0,2.0)*model.isProvenInfeasible());     // Is infeasiblity proven (or none better than cutoff)?
  cbc_status += (int)(pow(2.0,3.0)*model.isContinuousUnbounded());  // Was continuous solution unbounded.
  cbc_status += (int)(pow(2.0,4.0)*model.isProvenDualInfeasible()); // Was continuous solution unbounded.
  cbc_status += (int)(pow(2.0,5.0)*model.isNodeLimitReached());     // Node limit reached?
  cbc_status += (int)(pow(2.0,6.0)*model.isSecondsLimitReached());  // Time limit reached?
  cbc_status += (int)(pow(2.0,7.0)*model.isSolutionLimitReached()); // Solution limit reached?
  
  iteration_count  = model.getIterationCount(); // Get how many iterations it took to solve the problem.
  node_count       = model.getNodeCount();      // Get how many Nodes it took to solve the problem.
  secondary_status = model.secondaryStatus();   // Secondary status of problem:
  //                                               - -1 unset (status_ will also be -1) 
  //                                               -  0 search completed with solution 
  //                                               -  1 linear relaxation not feasible (or worse than cutoff) 
  //                                               -  2 stopped on gap 
  //                                               -  3 stopped on nodes 
  //                                               -  4 stopped on time 
  //                                               -  5 stopped on user event 
  //                                               -  6 stopped on solutions 
  //                                               -  7 linear relaxation unbounded.
  
  initsolve_status  = (int)(pow(2.0,0.0)*model.isInitialSolveAbandoned());              // Are there numerical difficulties (for initialSolve) ?
  initsolve_status += (int)(pow(2.0,1.0)*model.isInitialSolveProvenOptimal());          // Is optimality proven (for initialSolve) ?
  initsolve_status += (int)(pow(2.0,2.0)*model.isInitialSolveProvenPrimalInfeasible()); // Is primal infeasiblity proven (for initialSolve) ?
  initsolve_status += (int)(pow(2.0,3.0)*model.isInitialSolveProvenDualInfeasible());   // Is dual infeasiblity proven (for initialSolve) ? 

  //////////////////////////////
  // Allocate for return data //
  //////////////////////////////

#ifdef DEBUG
  sciprint("DEBUG: allocating data\n");
#endif

  int m_xmin             = ncols, n_xmin             = 1;
  int m_fmin             = 1,     n_fmin             = 1;
  int m_status           = 1,     n_status           = 1;
  int m_lambda           = 1,     n_lambda           = nrows;
  int m_redcosts         = 1,     n_redcosts         = ncols;
  int m_time             = 1,     n_time             = 1;
  int m_cbc_status       = 1,     n_cbc_status       = 1;
  int m_iteration_count  = 1,     n_iteration_count  = 1;
  int m_node_count       = 1,     n_node_count       = 1;
  int m_secondary_status = 1,     n_secondary_status = 1;
  int m_initsolve_status = 1,     n_initsolve_status = 1;

  int m_list_labels,  n_list_labels;
  double * xmin = NULL, * redcosts = NULL, * lambda = NULL, fmin;
  
  _SciErr = allocMatrixOfDouble(pvApiCtx, XMIN_OUT, m_xmin, n_xmin, &xmin);

#ifdef DEBUG
  sciprint("DEBUG: solving\n");
#endif

  status           = (double)model.status();         // Status of problem:
                                                     // - -1 - unknown e.g. before solve or if postSolve says not optimal
                                                     // -  0 - optimal
                                                     // -  1 - primal infeasible
                                                     // -  2 - dual infeasible
                                                     // -  3 - stopped on iterations or time
                                                     // -  4 - stopped due to errors
                                                     // -  5 - stopped by event handler (virtual int ClpEventHandler::event())
  fmin              = model.getObjValue(); // Get best objective function value as minimization.
  res_time          = cpuTime() - time0;
  cbc_status        = (double)cbc_status;
  iteration_count   = (double)iteration_count;
  node_count        = (double)node_count;
  secondary_status  = (double)secondary_status;
  initsolve_status  = (double)initsolve_status;

  /////////////////////////////////
  // Copy solutions if available //
  /////////////////////////////////

#ifdef DEBUG
  sciprint("DEBUG: returning data\n");
#endif

  lambda   = (double *)MALLOC(sizeof(double)*nrows);
  redcosts = (double *)MALLOC(sizeof(double)*ncols);

  if (model.status()>=0)
    {
      // Get pointer to array[getNumCols()] of primal solution vector.
      memcpy(xmin,     model.getColSolution(), sizeof(double)*ncols);
      // Get pointer to array[getNumRows()] of dual prices.
      memcpy(lambda,   model.getRowPrice(),    sizeof(double)*nrows);
      // Get a pointer to array[getNumCols()] of reduced costs.
      memcpy(redcosts, model.getReducedCost(), sizeof(double)*ncols);
    }
  else
    {
      memset(xmin,     0, sizeof(double)*ncols);
      memset(lambda,   0, sizeof(double)*nrows);
      memset(redcosts, 0, sizeof(double)*ncols);
    }

  // Create the 'extra' structure of type plist
  int * extra_addr = NULL;
  char * ListLabels [] = {"lambda","redcosts","time","cbc_status","iteration_count",
                          "node_count","secondary_status","initsolve_status"};

  m_list_labels = 9; n_list_labels = 1;

  _SciErr = createMatrixOfDouble(pvApiCtx, XMIN_OUT,   m_xmin,   n_xmin,   xmin);
  _SciErr = createMatrixOfDouble(pvApiCtx, FMIN_OUT,   m_fmin,   n_fmin,   &fmin);
  _SciErr = createMatrixOfDouble(pvApiCtx, STATUS_OUT, m_status, n_status, &status);

  _SciErr = createPList(pvApiCtx, EXTRA_OUT, &extra_addr, (char **)ListLabels, 8); SCICOINOR_ERROR;

  _SciErr = createColVectorOfDoubleInPList(pvApiCtx, EXTRA_OUT, extra_addr, "lambda",   m_lambda*n_lambda,     lambda); SCICOINOR_ERROR;
  _SciErr = createColVectorOfDoubleInPList(pvApiCtx, EXTRA_OUT, extra_addr, "redcosts", m_redcosts*n_redcosts, redcosts); SCICOINOR_ERROR;
  _SciErr = createDoubleInPList(pvApiCtx, EXTRA_OUT,            extra_addr, "time",             res_time); SCICOINOR_ERROR;
  _SciErr = createIntInPList(pvApiCtx, EXTRA_OUT,               extra_addr, "cbc_status",       cbc_status); SCICOINOR_ERROR;
  _SciErr = createIntInPList(pvApiCtx, EXTRA_OUT,               extra_addr, "iteration_count",  iteration_count); SCICOINOR_ERROR;
  _SciErr = createIntInPList(pvApiCtx, EXTRA_OUT,               extra_addr, "node_count",       node_count); SCICOINOR_ERROR;
  _SciErr = createIntInPList(pvApiCtx, EXTRA_OUT,               extra_addr, "secondary_status", secondary_status); SCICOINOR_ERROR;
  _SciErr = createIntInPList(pvApiCtx, EXTRA_OUT,               extra_addr, "initsolve_status", initsolve_status); SCICOINOR_ERROR;

  LhsVar(1) = XMIN_OUT;
  LhsVar(2) = FMIN_OUT;
  LhsVar(3) = STATUS_OUT;
  LhsVar(4) = EXTRA_OUT;

  //////////////////////////////
  // Delete allocated objects //
  //////////////////////////////

  //if (modelSimplexCopy) delete modelSimplexCopy;
  if (printer)    delete printer;
  if (clpprinter) delete clpprinter;
  if (A_matrix)   delete A_matrix;
  if (Q_matrix)   delete Q_matrix;

  if (lambda)   FREE(lambda);
  if (redcosts) FREE(redcosts);

  return 0;
}
