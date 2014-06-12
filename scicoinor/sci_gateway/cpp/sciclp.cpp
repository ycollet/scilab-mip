////////////////////////////////////////////////////////////////////////////////////////
// sciclp: A scilab interface to the CLP library for linear and quadratic programming //
////////////////////////////////////////////////////////////////////////////////////////

//  Copyright (C) 2008-2010 Yann Collette.
//  Copyright (C) 2001-2007 Johan Loefberg, ETH Zurich, loefberg@control.ee.ethz.ch
//
//  SCICLP is free software; you can redistribute it and/or modify it
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

#include <cstdio>
#include <cstring>
#include <exception>

#undef min
#undef max

#include <Coin_C_defines.h>
#include <CoinPackedMatrix.hpp>
#include <CoinMessageHandler.hpp>
#include <CoinError.hpp>
#include <ClpSimplex.hpp>
#include <ClpInterior.hpp>
#include <ClpPresolve.hpp>
#include <ClpSolve.hpp>
#include <ClpCholeskyBase.hpp>
#include <ClpPdcoBase.hpp>
#include <ClpEventHandler.hpp>

#include <OsiSolverInterface.hpp>
#include <OsiClpSolverInterface.hpp>

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
#define Q_IN        9
#define PARAM_IN    10
#define LASTPARAM   10

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
//       sciprint("sciclp: error %s line %d %s::%s - %s\n", e.fileName().c_str(), \
// 	       e.lineNumber(),						\
// 	       e.className().c_str(),					\
// 	       e.methodName().c_str(),					\
// 	       e.message());						\
//       									\
//       if (modelSimplex)  delete modelSimplex;				\
//       if (modelInterior) delete modelInterior;				\
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
//       sciprint("sciclp: exception raised - %s\n", e.what());		\
//       									\
//       if (modelSimplex)  delete modelSimplex;				\
//       if (modelInterior) delete modelInterior;				\
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
//       sciprint("sciclp: exception raise - %s\n", e.what());		\
//       									\
//       if (modelSimplex)  delete modelSimplex;				\
//       if (modelInterior) delete modelInterior;				\
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
//       sciprint("sciclp: unknow exception raised\n");			\
//       throw;								\
//     }

// To define the function sciclp as a C function and allow this function to be loaded easily in scilab
extern "C" int sciclp(char * fname)
{
  ClpSimplex      * modelSimplex  = NULL;
  ClpInterior     * modelInterior = NULL;
  ClpModel        * modelBase     = NULL;
  ClpCholeskyBase * cholesky      = NULL; 
  DerivedHandler  * printer       = NULL;
  ClpSolve          modelSolve; // YC: we can set some parameters via CbcSolve
  CoinPackedMatrix  A_matrix, Q_matrix;
#ifdef HANDLE_CTRLC
  MyClpEventHandler SciClpEventHandler;
  SciClpEventHandler.setInCbc(false);
#endif
  int Log = 0;
  int ncols = 0, nrows = 0, i, j, nz = 0, writemps = 0;
  int count = 0, status = 0, type;
  char  * writemps_filename = NULL;
  SciSparse S_A, S_Q;
  
  if (Rhs<LASTPARAM) 
    {
      Scierror(999,"%s: 12 inputs required in call to %s. Bug in clp.sci ?...\n", fname, fname);
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
  int n_q,     m_q,     * q_addr     = NULL;
  double * c = NULL, * lhs = NULL, * rhs = NULL, * lower = NULL, * upper = NULL, * a = NULL;
  double * q = NULL;
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
  int     solverchoice = 1, loglevel = 0, fact_freq = 200, perturb = 0;
  int     presolve = 0, red_grad = 1, maximumbarrieriterations = 200;
  int     tmp_int, tmp_res;
  double  tmp_double;
  char *  tmp_char;

  initPList(pvApiCtx, PARAM_IN, &param_in_addr);
  if (!checkPList(pvApiCtx, param_in_addr))
    {
      Scierror(999, "%s: argument n° %d is not a plist\n", fname, PARAM_IN);

      return 0;
    }

  // solver option
  getIntInPList(pvApiCtx, param_in_addr, "solver", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) solverchoice = tmp_int;

  // Load the matrix in the solver
  switch(solverchoice)
    {
    case 6:
      modelInterior = new ClpInterior;
      modelBase = (ClpModel *)modelInterior;
      break;
    default:
      modelSimplex = new ClpSimplex;
      modelBase = (ClpModel *)modelSimplex;
      break;
    }

  // maxnumiterations option
  getIntInPList(pvApiCtx, param_in_addr, "maxnumiterations", &tmp_int, &tmp_res, 1000, Log, CHECK_NONE);
  if (tmp_res!=-1) modelBase->setMaximumIterations(tmp_int);
      
  // maxnumseconds option
  getDoubleInPList(pvApiCtx, param_in_addr, "maxnumseconds", &tmp_double, &tmp_res, 3600, Log, CHECK_NONE);
  if (tmp_res!=-1) modelBase->setMaximumSeconds(tmp_double);
      
  // primaltolerance option
  getDoubleInPList(pvApiCtx, param_in_addr, "primaltolerance", &tmp_double, &tmp_res, 1e-7, Log, CHECK_NONE);
  if (tmp_res!=-1) modelBase->setPrimalTolerance(tmp_double);

  // dualtolerance option
  getDoubleInPList(pvApiCtx, param_in_addr, "dualtolerance", &tmp_double, &tmp_res, 1e-7, Log, CHECK_NONE);
  if (tmp_res!=-1)   modelBase->setDualTolerance(tmp_double);				

  // verbose option
  getIntInPList(pvApiCtx, param_in_addr, "verbose", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) 
    {
      ///////////////////////////////
      // Enable printing in Scilab //
      // and set parameters of clp //
      ///////////////////////////////
	  
      loglevel = tmp_int;

      printer = new DerivedHandler(); // assumed open	
      printer->setLogLevel(loglevel);
    }

  // optim_dir option
  getIntInPList(pvApiCtx, param_in_addr, "optim_dir", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) modelBase->setOptimizationDirection(tmp_int);

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

  // perturb option
  getIntInPList(pvApiCtx, param_in_addr, "perturb", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) perturb = tmp_int;

  // scaling option
  // Sets or unsets scaling:
  // - 0 -off
  // - 1 equilibrium (default)
  // - 2 geometric
  // - 3 auto
  // - 4 auto-but-as-initialSolve-in-bab. 
  getIntInPList(pvApiCtx, param_in_addr, "scaling", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) modelBase->scaling(tmp_int);

  // factorization frequency option
  getIntInPList(pvApiCtx, param_in_addr, "fact_freq", &tmp_int, &tmp_res, 200, Log, CHECK_NONE);
  if (tmp_res!=-1) fact_freq = tmp_int;

  // presolve option
  getIntInPList(pvApiCtx, param_in_addr, "presolve", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) presolve = tmp_int;

  // reduced gradient phase option
  getIntInPList(pvApiCtx, param_in_addr, "red_grad", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) red_grad = tmp_int;

  // maxnumiterationshotstart
  getIntInPList(pvApiCtx, param_in_addr, "maxnumiterationshotstart", &tmp_int, &tmp_res, 1000, Log, CHECK_NONE);
  if (tmp_res!=-1) modelBase->setIntParam(ClpMaxNumIterationHotStart, tmp_int);

  // dualobjectivelimit
  getDoubleInPList(pvApiCtx, param_in_addr, "dualobjectivelimit", &tmp_double, &tmp_res, 1e-7, Log, CHECK_NONE);
  if (tmp_res!=-1) modelBase->setDblParam(ClpDualObjectiveLimit, tmp_double);

  // primalobjectivelimit
  getDoubleInPList(pvApiCtx, param_in_addr, "primalobjectivelimit", &tmp_double, &tmp_res, 1e-7, Log, CHECK_NONE);
  if (tmp_res!=-1) modelBase->setDblParam(ClpPrimalObjectiveLimit, tmp_double);

  // obj offset
  getDoubleInPList(pvApiCtx, param_in_addr, "objoffset", &tmp_double, &tmp_res, 0.0, Log, CHECK_NONE);
  if (tmp_res!=-1) modelBase->setDblParam(ClpObjOffset, tmp_double);

  // presolvetolerance
  getDoubleInPList(pvApiCtx, param_in_addr, "presolvetolerance", &tmp_double, &tmp_res, 1e-7, Log, CHECK_NONE);
  if (tmp_res!=-1) modelBase->setDblParam(ClpPresolveTolerance, tmp_double);
      
  // maximumBarrierIterations
  getIntInPList(pvApiCtx, param_in_addr, "maxnumbarrieriterations", &tmp_int, &tmp_res, 200, Log, CHECK_NONE);
  if (tmp_res!=-1) maximumbarrieriterations = tmp_int;

  /////////////////////////////////////
  // Set the bounds on the variables //
  /////////////////////////////////////

  modelBase->loadProblem(A_matrix,lhs,rhs,c,lower,upper);

  for(i=0;i<ncols; i++)
    {
      modelBase->setColUpper(i, *(upper+i));
      modelBase->setColLower(i, *(lower+i));
      modelBase->setObjectiveCoefficient(i,*(c+i));
      if ((*(vtype+i)=='I')||(*(vtype+i)=='i'))
	{
	  modelBase->setInteger(i);
	}
      else
	{
	  modelBase->setContinuous(i);
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
	  modelBase->setRowUpper(i, *(rhs+i));
	  modelBase->setRowLower(i, -COIN_DBL_MAX);
	  break;
	case 'e':
	case 'E':
	  modelBase->setRowUpper(i, *(rhs+i));
	  modelBase->setRowLower(i, *(rhs+i));
	  break;
	case 'n':
	case 'N':
	  modelBase->setRowUpper(i,  COIN_DBL_MAX);
	  modelBase->setRowLower(i, -COIN_DBL_MAX);
	  break;
	case 'r':
	case 'R':
	  modelBase->setRowUpper(i, *(rhs+i));
	  modelBase->setRowLower(i, *(lhs+i));
	  break;
	case 'g':
	case 'G':
	default:
	  modelBase->setRowUpper(i, COIN_DBL_MAX);
	  modelBase->setRowLower(i, *(lhs+i));
	  break;
	}
#ifdef DEBUG
      sciprint("row lower[%d] = %f row upper[%d] = %f\n", i, modelBase->getRowLower()[i], i, modelBase->getRowUpper()[i]);
#endif
    }
  
  ///////////////////////////////////
  // Affect names to rows and cols //
  ///////////////////////////////////

  modelBase->setIntParam(ClpNameDiscipline,0);

  ////////////////////////
  // Any quadratic part //
  ////////////////////////
  
#ifdef DEBUG
  sciprint("DEBUG: dealing with Q\n");
#endif
  _SciErr = getVarAddressFromPosition(pvApiCtx, Q_IN, &q_addr); SCICOINOR_ERROR;
  _SciErr = getVarType(pvApiCtx, q_addr, &type); SCICOINOR_ERROR;

  if(type!=sci_sparse)
    {
      _SciErr = getMatrixOfDouble(pvApiCtx, q_addr, &m_q, &n_q, &q); SCICOINOR_ERROR;
      
      if (n_q * m_q !=0)
	{
#ifdef DEBUG
	  sciprint("DEBUG: Q_IN is not sparse\n");
	  sciprint("DEBUG: Q size: n = %d m= %d\n", n_q, m_q);
#endif
	  
	  Q_matrix.setDimensions(nrows,ncols);
	  
	  nz = n_q * m_q;
	  
	  if (q==NULL) 
	    {
	      Scierror(999,"%s: invalid value of matrix q\n",fname);
	  
	      freeAllocatedSingleString(ctype);
	      freeAllocatedSingleString(vtype);

	      if (writemps_filename) FREE(writemps_filename);

	      return 0;
	    }
	  
	  for(i=0; i<m_q; i++)
	    {
	      for(j=0; j<n_q; j++)
		{
		  if (*(q+i+j*m_q)!=0) Q_matrix.modifyCoefficient(i,j,*(q+i+j*m_q));
		}
	    }
	  modelBase->loadQuadraticObjective(Q_matrix);
	}
    }
  else
    {
      getAllocatedSparseMatrix(pvApiCtx, q_addr, &S_Q.m, &S_Q.n, &S_Q.nel, &S_Q.mnel, &S_Q.icol, &S_Q.R);
      
      if (S_Q.n * S_Q.m!=0)
	{
#ifdef DEBUG
	  sciprint("DEBUG: Q_IN is sparse\n");
	  sciprint("DEBUG: Q size: n = %d m= %d\n", S_Q.n, S_Q.m);
#endif
	  Q_matrix.setDimensions(nrows,ncols);
	  
	  nz = S_Q.nel;
	  
	  count = 0;
	  for(i=0;i<S_Q.m;i++)
	    {
	      if (S_Q.mnel[i]!=0) 
		{
		  for(j=0;j<S_Q.mnel[i];j++)
		    {
		      count++;
		      Q_matrix.modifyCoefficient(i,S_Q.icol[count-1]-1,S_Q.R[count-1]);
		    }
		}
	    }
	  modelBase->loadQuadraticObjective(Q_matrix);
	}

      freeAllocatedSparseMatrix(S_Q.mnel, S_Q.icol, S_Q.R);
    }
 
  // Solver specific part
  switch(solverchoice)
    {
    case 6:
      modelInterior->setMaximumBarrierIterations(maximumbarrieriterations);
      break;
    default:
      modelSimplex->setFactorizationFrequency(fact_freq);
      modelSimplex->setPerturbation(perturb);
      break;
    }

  ///////////////////////////
  // Pre / Post resolution //
  ///////////////////////////

#ifdef HANDLE_CTRLC
  ///////////////////////////////////////////////////////////
  // Pass in the event handler, to handle ctrl+c in scilab //
  ///////////////////////////////////////////////////////////

  if (modelSimplex) 
    {
      modelSimplex->passInEventHandler(&SciClpEventHandler);
      modelSimplex->passInMessageHandler(printer);
    }
#endif

#ifdef DEBUG
  sciprint("DEBUG: presolve\n");
#endif
  // If the solver is not ClpInterior we can set the presolve option
  if (presolve && solverchoice!=6) 
    {
      TRYCATCH(modelSimplex->initialSolve(modelSolve));
      switch(presolve)
	{
	case 2:
	  TRYCATCH(modelSimplex->initialDualSolve());
	  break;
	case 3:
	  TRYCATCH(modelSimplex->initialPrimalSolve());
	  break;
	case 4:
	  TRYCATCH(modelSimplex->initialBarrierSolve());
	  break;
	case 5:
	  TRYCATCH(modelSimplex->initialBarrierNoCrossSolve());
	  break;
	default:
	  TRYCATCH(modelSimplex->initialSolve());
	  break;
	}
    }

  ////////////////////////////////////////////
  // If needed, write the problem in a file //
  ////////////////////////////////////////////

#ifdef DEBUG
  sciprint("DEBUG: dealing with writemps\n");
#endif
  if (writemps)
    {
      modelBase->writeMps(writemps_filename);
      if (loglevel) sciprint("sciclp: writing %s mps file\n",writemps_filename);
    }

  ////////////////
  // Resolution //
  ////////////////

#ifdef DEBUG
  sciprint("DEBUG: resolution\n");
  sciprint("model number columns = %d\n", modelBase->numberColumns());
  sciprint("model number rows    = %d\n", modelBase->numberRows());
#endif

  if ((solverchoice>=1)&&(solverchoice<=5))
    {
      TRYCATCH(modelSimplex->initialSolve())
    }

#ifdef HANDLE_CTRLC
  ///////////////////////////////////////////////////////////
  // Pass in the event handler, to handle ctrl+c in scilab //
  ///////////////////////////////////////////////////////////

  if (modelSimplex)  
    {
      modelSimplex->passInEventHandler(&SciClpEventHandler);
      modelSimplex->passInMessageHandler(printer);
    }
  if (modelInterior) 
    {
      modelInterior->passInEventHandler(&SciClpEventHandler);
      modelInterior->passInMessageHandler(printer);
    }
#endif

  switch (solverchoice)
    {
    case 2:
      TRYCATCH(status = modelSimplex->dual());
      break;
    case 3:
      TRYCATCH(status = modelSimplex->barrier(false));
      break;
    case 4:
      TRYCATCH(status = modelSimplex->barrier(true));
      break;
    case 5:
      TRYCATCH(status = modelSimplex->reducedGradient(red_grad));
      break;
    case 6:
      cholesky = new ClpCholeskyBase();
      cholesky->setKKT(true);
      modelInterior->setCholesky(cholesky);
      TRYCATCH(status = modelInterior->primalDual());
      break;
    case 7:
      TRYCATCH(status = modelInterior->pdco());
      break;
    default:
      TRYCATCH(status = modelSimplex->primal());
      break;
    }

#ifdef DEBUG
  if (status && loglevel) sciprint("Optimization failed\n");
#endif

  int clp_status = 0;

  clp_status  = (int)(pow(2.0,0.0)*modelBase->isAbandoned());
  clp_status += (int)(pow(2.0,1.0)*modelBase->isProvenOptimal());
  clp_status += (int)(pow(2.0,2.0)*modelBase->isProvenPrimalInfeasible());
  clp_status += (int)(pow(2.0,3.0)*modelBase->isProvenDualInfeasible());
  clp_status += (int)(pow(2.0,4.0)*modelBase->isPrimalObjectiveLimitReached());
  clp_status += (int)(pow(2.0,5.0)*modelBase->isDualObjectiveLimitReached());  
  clp_status += (int)(pow(2.0,6.0)*modelBase->isIterationLimitReached());      

  //////////////////////////////
  // Allocate for return data //
  //////////////////////////////

#ifdef DEBUG
  sciprint("DEBUG: allocating data\n");
#endif

  int m_xmin   = ncols, n_xmin   = 1;
  int m_lambda   = 1, n_lambda   = nrows;
  int * extra_addr = NULL;

  char * ListLabels [] = {"lambda","secondary_status","clp_status"};

  _SciErr = createMatrixOfDouble(pvApiCtx, XMIN_OUT, m_xmin, n_xmin, (double *)modelBase->primalColumnSolution()); SCICOINOR_ERROR;
  createScalarDouble(pvApiCtx, FMIN_OUT, modelBase->getObjValue());
  createScalarDouble(pvApiCtx, STATUS_OUT, (double)modelBase->status());

  _SciErr = createPList(pvApiCtx, EXTRA_OUT, &extra_addr, (char **)ListLabels, 3); SCICOINOR_ERROR;

  _SciErr = createColVectorOfDoubleInPList(pvApiCtx, EXTRA_OUT, extra_addr, "lambda", m_lambda*n_lambda, (double *)modelBase->dualRowSolution()); SCICOINOR_ERROR;
  _SciErr = createIntInPList(pvApiCtx, EXTRA_OUT, extra_addr, "secondary_status", modelBase->secondaryStatus()); SCICOINOR_ERROR;
  _SciErr = createIntInPList(pvApiCtx, EXTRA_OUT, extra_addr, "clp_status",       clp_status); SCICOINOR_ERROR;

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

  if (modelSimplex)  delete modelSimplex;
  if (modelInterior) delete modelInterior;
  if (printer)       delete printer;

  freeAllocatedSingleString(ctype);
  freeAllocatedSingleString(vtype);

  if (writemps_filename) FREE(writemps_filename);

  return 0;
}
