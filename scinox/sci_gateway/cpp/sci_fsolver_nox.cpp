// configure options for trilinos-9.0.3:
// ./configure --prefix=/opt/stow/trilinos-9.0.3 --enable-shared --enable-moocho --enable-nox --enable-anasazi --enable-nox-epetra --enable-nox-epetraext --enable-amesos

#define DEBUG 1

// The Amesos solver seems to be still in devel
#define USE_AZTECOO 1

#include <iostream>
#include <streambuf>
#include <string>

// NOX Objects
#include "NOX.H"
#include "NOX_Epetra.H"
#include "Epetra_Operator.h"
#include "NOX_Epetra_Interface_Required.H" // base class
#include "NOX_Epetra_Interface_Jacobian.H" // base class
#include "NOX_Epetra_Interface_Preconditioner.H" // base class
#ifndef USE_AZTECOO
#include "NOX_Epetra_LinearSystem_Amesos.H"
#else
#include "NOX_Epetra_LinearSystem_AztecOO.H"
#endif

// Trilinos Objects
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"

#include <string.h>
#include <setjmp.h>

extern "C" {
#include <stack-c.h>
#include <sciprint.h>
#include <Scierror.h>
#include <MALLOC.h>
}

#include <api_scilab.h>

#include <parameters.hpp>

#define _(A) (A)
#define MAX(A,B) ((A<B)?B:A)

#define F_IN           1
#define X_IN           2 
#define PARAM_IN       3
#define X_OUT          4
#define INFO_NOX_OUT   5
#define ITERATIONS_OUT 6
#define EXTRA_OUT      7
#define FUNC_TMP_1     Rhs+1
#define FUNC_TMP_2     Rhs+2

////////////////////////////////////
// Some required static variables //
////////////////////////////////////

static jmp_buf call_f_env; 
static int sci_obj, lhs_obj, rhs_obj;

////////////////////////////////////////////////////
// The class which interfaces the problem and NOX //
////////////////////////////////////////////////////

class  Problem_Interface : public NOX::Epetra::Interface::Required,
			   public NOX::Epetra::Interface::Jacobian,
			   public NOX::Epetra::Interface::Preconditioner
{
public:
  ~Problem_Interface();

  //! Compute and return F
  bool computeF(const Epetra_Vector& x, Epetra_Vector& FVec, FillType flag = Residual);

  //! Compute an explicit Jacobian
  bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);

  //! Compute the Epetra_RowMatrix M, that will be used by the Aztec preconditioner instead of the Jacobian.  
  /*! This is used when there is no explicit Jacobian present (i.e. Matrix-Free Newton-Krylov).  
   *  This MUST BE an Epetra_RowMatrix since the Aztec preconditioners need to know the sparsity pattern of the matrix.
   *  Returns true if computation was successful.
   */
  bool computePrecMatrix(const Epetra_Vector& x, Epetra_RowMatrix& M);
  
  //! Computes a user supplied preconditioner based on input vector x.  Returns true if computation was successful.
  bool computePreconditioner(const Epetra_Vector& x, Epetra_Operator& M, Teuchos::ParameterList* precParams = 0);
};

/////////////////////////////////////////
// Implementation of Problem_Interface //
/////////////////////////////////////////

Problem_Interface::~Problem_Interface() { }

bool Problem_Interface::computeF(const Epetra_Vector& x, Epetra_Vector& FVec, FillType flag)
{
  int m_x   = x.GlobalLength(), n_x   = 1, l_x;
  int m_jac = 1,                n_jac = 1, l_jac;
  int i, size_x = x.GlobalLength(), size_f = FVec.GlobalLength();
  int rhs_old = Rhs, nbvars_old = Nbvars;
  int my_lhs_obj = 1;

  Rhs    = FUNC_TMP_1 + MAX(rhs_obj,lhs_obj);
  Nbvars = FUNC_TMP_1 + MAX(rhs_obj,lhs_obj);

  CreateVar(FUNC_TMP_1, MATRIX_OF_DOUBLE_DATATYPE, &m_x, &n_x, &l_x);
  for(i=0;i<size_x;i++) 
    {
      *stk(l_x+i) = x[i];
    }

  CreateVar(FUNC_TMP_2, MATRIX_OF_DOUBLE_DATATYPE, &m_jac, &n_jac, &l_jac);

  PExecSciFunction(FUNC_TMP_1, &sci_obj, &my_lhs_obj, &rhs_obj, "FormFunction", call_f_env);

  GetRhsVar(FUNC_TMP_1, MATRIX_OF_DOUBLE_DATATYPE, &m_x, &n_x, &l_x);

  if (m_x*n_x!=size_f)
    {
      Scierror(999,_("fsolver_nox: error, you must return %d values\n"),size_f);
      Rhs    = rhs_old;
      Nbvars = nbvars_old;
      return 0;
    }

  for(i=0;i<size_f;i++) 
    {
      FVec[i] = *stk(l_x+i);
#ifdef DEBUG
      sciprint("DEBUG: fsolver_nox computeF: FVec[%d] = %f stk[%d] = %f\n", i, FVec[i],i,*stk(l_x+i));
#endif
    }

  Rhs    = rhs_old;
  Nbvars = nbvars_old;

  return true;
}

bool Problem_Interface::computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac)
{
  int m_x   = x.GlobalLength(), n_x   = 1, l_x;
  int m_jac = 1,                n_jac = 1, l_jac;
  int i, j, irow, jcol;
  int rhs_old = Rhs, nbvars_old = Nbvars;
  int size_x = x.GlobalLength(), size_f;
  double A;

  Rhs    = FUNC_TMP_1 + MAX(rhs_obj,lhs_obj);
  Nbvars = FUNC_TMP_1 + MAX(rhs_obj,lhs_obj);

  Epetra_CrsMatrix* Jacobian = dynamic_cast<Epetra_CrsMatrix*>(&Jac);
  if (Jacobian == NULL) 
    {
      sciprint("Error: Problem_Interface::computeJacobian() - The supplied");
      sciprint(" Epetra_Operator is NOT an Epetra_CrsMatrix!\n");
      Rhs    = rhs_old;
      Nbvars = nbvars_old;
      throw;
    }

  size_f = Jacobian->NumGlobalRows();

#ifdef DEBUG
  sciprint("DEBUG: before - computeJacobian\n");
  sciprint("DEBUG: the Jacobian\n");
  Jacobian->Print(cout);
  int * Indices   = NULL;
  double * Values = NULL;
  int num_entries = 0;
  Indices = (int *)MALLOC(size_x*sizeof(int));
  Values  = (double *)MALLOC(size_x*sizeof(double));

  for(i=0;i<size_f;i++)
    {
      Jacobian->ExtractGlobalRowCopy(i, size_f, num_entries, Values, Indices);
      sciprint("Line %d:",i);
      for(j=0;j<num_entries;j++)
	{
	  sciprint(" A[%d] = %f",Indices[j],Values[j]);
	}
      sciprint("\n");
    }

  FREE(Indices);
  FREE(Values);

  sciprint("DEBUG: x\n");
  x.Print(cout);
#endif

  CreateVar(FUNC_TMP_1, MATRIX_OF_DOUBLE_DATATYPE, &m_x, &n_x, &l_x);
  for(i=0;i<size_x;i++) 
    {
      *stk(l_x+i) = x[i];
    }

  CreateVar(FUNC_TMP_2, MATRIX_OF_DOUBLE_DATATYPE, &m_jac, &n_jac, &l_jac);

  PExecSciFunction(FUNC_TMP_1, &sci_obj, &lhs_obj, &rhs_obj, "FormJacobian", call_f_env);

  if (GetType(FUNC_TMP_2)==5) // Sparse matrix
    {
      SciSparse SpResult;
      int i, j, Index;

      GetRhsVar(FUNC_TMP_2, SPARSE_MATRIX_DATATYPE, &m_jac, &n_jac, &SpResult);

      if ((m_jac!=size_f)&&(n_jac!=size_x))
	{
	  Scierror(999,_("fsolver_nox: error, you must return a %d * %d sparse matrix\n"),size_f,size_x);
	  Rhs    = rhs_old;
	  Nbvars = nbvars_old;
	}

      Index = 0;
      for(i=0;i<SpResult.m;i++)
	{
	  for(j=0;j<SpResult.mnel[i];j++)
	    {
	      A = SpResult.R[Index];
	      irow = i;
	      jcol = SpResult.icol[Index]-1;
	      Index++;
#ifdef DEBUG
	      sciprint("DEBUG: fsolver_nox computeJacobian: J[%d][%d] = %f\n", irow, jcol, A);
#endif
	      Jacobian->ReplaceGlobalValues(irow, 1, &A, &jcol);
	    }
	}
    }
  else
    {
      GetRhsVar(FUNC_TMP_2, MATRIX_OF_DOUBLE_DATATYPE, &m_jac, &n_jac, &l_jac);
      
      if ((m_jac!=size_f)&&(n_jac!=size_x))
	{
	  Scierror(999,_("fsolver_nox: error, you must return a %d * %d full matrix\n"),size_f,size_x);
	  Rhs    = rhs_old;
	  Nbvars = nbvars_old;
	}
      
      for(i=0;i<m_jac;i++)
	{
	  for(j=0;j<n_jac;j++)
	    {
	      A = *stk(l_jac+i+j*m_jac);
	      irow = i;
	      jcol = j;
#ifdef DEBUG
	      sciprint("DEBUG: fsolver_nox computeJacobian: J[%d][%d] = %f\n", irow, jcol, A);
#endif
	      Jacobian->ReplaceGlobalValues(irow, 1, &A, &jcol);
	    }
	}
    }

  Jacobian->FillComplete(false); // If true, a storage optimization is performed. The sparsity structure will not be conserved !
  Jacobian->OptimizeStorage();
  Jacobian->MakeDataContiguous();

#ifdef DEBUG
  sciprint("DEBUG: after - computeJacobian\n");
  Jacobian->Print(cout);
  sciprint("DEBUG: the Jacobian\n");
  Indices = (int *)MALLOC(size_x*sizeof(int));
  Values  = (double *)MALLOC(size_x*sizeof(double));

  for(i=0;i<size_f;i++)
    {
      Jacobian->ExtractGlobalRowCopy(i, size_f, num_entries, Values, Indices);
      sciprint("Line %d:",i);
      for(j=0;j<num_entries;j++)
	{
	  sciprint(" A[%d] = %f",Indices[j],Values[j]);
	}
      sciprint("\n");
    }

  Rhs    = rhs_old;
  Nbvars = nbvars_old;

  FREE(Indices);
  FREE(Values);
#endif

  return true;
}

bool Problem_Interface::computePrecMatrix(const Epetra_Vector& x, Epetra_RowMatrix& M)
{
  Epetra_CrsMatrix* precMatrix = dynamic_cast<Epetra_CrsMatrix*>(&M);

#ifdef DEBUG
  sciprint("DEBUG: computePrecMatrix\n");
#endif

  if (precMatrix == NULL) 
    {
      sciprint("Error: Problem_Interface::computePreconditioner() - The supplied");
      sciprint(" Epetra_Operator is NOT an Epetra_CrsMatrix!\n");
      throw;
    }

  return computeJacobian(x, M);
}
bool Problem_Interface::computePreconditioner(const Epetra_Vector& x, Epetra_Operator& M, Teuchos::ParameterList* precParams)
{
#ifdef DEBUG
  sciprint("DEBUG: computePreconditioner\n");
#endif

  sciprint("Error: Problem_Interface::preconditionVector() - Use Explicit Jacobian only for this test problem!\n");

  throw 1;
}

//////////////////////////
// For cout redirection //
//////////////////////////

class ScilabStream : public std::basic_streambuf<char>
{
public:
  ScilabStream(std::ostream &stream) : m_stream(stream)
  {
    m_old_buf = stream.rdbuf();
    stream.rdbuf(this);
  }
  ~ScilabStream()
  {
    // output anything that is left
    if (!m_string.empty())
      sciprint("fsolve_nox: %s\n",m_string.c_str());

    m_stream.rdbuf(m_old_buf);
  }

protected:
  virtual int_type overflow(int_type v)
  {
    if (v == '\n')
      {
	sciprint("fsolve_nox: %s\n",m_string.c_str());
	m_string.clear();
      }
    else
      m_string.push_back(v);
    
    return v;
  }
  
  virtual std::streamsize xsputn(const char *p, std::streamsize n) 
  {
    m_string.append(p, p + n);
    
    int pos = 0;
    while (pos != std::string::npos)
      {
	pos = m_string.find('\n');
	if (pos != std::string::npos)
	  {
	    std::string tmp(m_string.begin(), m_string.begin() + pos);
	    sciprint("fsolve_nox: %s\n",tmp.c_str());
	    m_string.erase(m_string.begin(), m_string.begin() + pos + 1);
	  }
      }
    
    return n;
  }
  
private:
  std::ostream   &m_stream;
  std::streambuf *m_old_buf;
  std::string     m_string;
};

using namespace std;

/////////////////////////////////
// The Scilab interface to NOX //
/////////////////////////////////

extern "C" int sci_fsolver_nox(char * fname)
{
  int m_x_in,           n_x_in,           l_x_in;
  int m_x_out      = 1, n_x_out      = 1;
  int m_x_tmp      = 1, n_x_tmp      = 1, l_x_tmp;
  int m_jac_tmp    = 1, n_jac_tmp    = 1, l_jac_tmp;
  int n, i, j, size_x, size_f;
  int irow, jcol;
  int * param_addr = NULL;
  int * extra_addr = NULL;
  // variables for GET macros
  int      tmp_int, tmp_int_2, tmp_res = -1, tmp_res_2 = -1;
  double   tmp_double, tmp_double_2, A;
  double   tmp_dbl[1], * x_out = NULL;
  // some return variables
  int m_list_labels = 3, n_list_labels = 1;
  static const char * ListLabels [3] = {"plist","info_nox","iterations"};
  Teuchos::RCP<ScilabStream> scicout = Teuchos::rcp(new ScilabStream(std::cout));
  Teuchos::RCP<ScilabStream> scicerr = Teuchos::rcp(new ScilabStream(std::cerr));
  SciErr _SciErr;

  //
  // Input 1 (objective function)
  //

  // Here, we get a Scilab script function

  if (GetType(F_IN) != sci_c_function)
    {
      Scierror(999,_("%s: Wrong type for input argument #%d: A scilab function expected.\n"), fname, 1);
      return 0;
    }
  
  GetRhsVar(F_IN, EXTERNAL_DATATYPE, &lhs_obj, &rhs_obj, &sci_obj);

  //
  // Input 2 (x vector) *
  //

  // Here, we get the vector which will be sent as a parameter to the function
  // The size of the vector is stored in a static variable: size_x to be used later in the solver function
  GetRhsVar(X_IN, "d", &m_x_in, &n_x_in, &l_x_in);
  size_x = m_x_in * n_x_in;

  if (n_x_in != 1)
    {
      sciprint(_("%s: x must be vector of dimension n x 1"),fname);
      return 0;
    }

  n       = size_x;  // number of variables
  m_x_out = m_x_in;
  n_x_out = n_x_in;

  //
  // Input 3 parameter list
  //

  init_parameters(PARAM_IN, &param_addr);

  if (check_parameters(param_addr)) 
    {
      Scierror(999,_("%s: Argument %d is not a plist\n"),fname,PARAM_IN);
      return 0;
    }

#ifdef HAVE_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
#endif

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif
  Epetra_Map MapCol(size_x, 0, Comm);

  // Get the process ID and the total number of processors
  int MyPID   = Comm.MyPID();
  //int NumProc = Comm.NumProc();

  // Get the vector from the Problem
  Teuchos::RCP<Epetra_Vector> soln = Teuchos::rcp(new Epetra_Vector(Copy,MapCol,stk(l_x_in)));

  NOX::Epetra::Vector noxSoln(soln, NOX::Epetra::Vector::CreateView);

  ////////////////////////////
  // Begin Nonlinear Solver //
  ////////////////////////////

  // Create the top level parameter list
  Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr = Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList& nlParams = *nlParamsPtr.get();

  /////////////////////////////////////
  // Set the nonlinear solver method //
  /////////////////////////////////////
  // 'non_linear_solver':
  // - 1 "Line Search Based"
  // - 2 "Trust Region Based"
  // - 3 "Inexact Trust Region Based"
  // - 4 "Tensor Based"
  get_int_parameter(param_addr, "non_linear_solver", &tmp_int, &tmp_res, 1, CHECK_VALUES, 4, 1, 2, 3, 4);
  if (tmp_res!=-1)
    {
      switch(tmp_int)
	{
	case 1:
	  nlParams.set("Nonlinear Solver", "Line Search Based");
	  break;
	case 2:
	  nlParams.set("Nonlinear Solver", "Trust Region Based");
	  // "Minimum Trust Region Radius" ($\Delta_{\min}$) - Minimum allowable trust region radius. Defaults to 1.0e-6.
	  get_double_parameter(param_addr, "ns_trb_minimum_trust_region_radius", &tmp_double,&tmp_res, 1.0e-6, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("Trust Region Based").set("Minimum Trust Region Radius",tmp_double);

	  // "Maximum Trust Region Radius" ($\Delta_{\max}$) - Maximum allowable trust region radius. Defaults to 1.0e+10.
	  get_double_parameter(param_addr, "ns_trb_maximum_trust_region_radius", &tmp_double, &tmp_res, 1.0e10, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("Trust Region Based").set("Maximum Trust Region Radius",tmp_double);

	  // "Minimum Improvement Ratio" ($\rho_{\min}$) - Minimum improvement ratio to accept the step. Defaults to 1.0e-4.
	  get_double_parameter(param_addr, "ns_trb_minimum_improvement_ratio", &tmp_double, &tmp_res, 1.0e-4, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("Trust Region Based").set("Minimum Improvement Ratio",tmp_double);

	  // "Contraction Trigger Ratio" ($\rho_{\rm s}$) - If the improvement ratio is less than this value, then the trust region
	  // is contracted by the amount specified by the "Contraction Factor". Must be larger than "Minimum Improvement Ratio". 
	  // Defaults to 0.1.
	  get_double_parameter(param_addr, "ns_trb_contraction_trigger_ratio", &tmp_double, &tmp_res, 0.1, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("Trust Region Based").set("Contraction Trigger Ratio",tmp_double);

	  // "Contraction Factor" ($\beta_{\rm s}$) - See above. Defaults to 0.25.
	  get_double_parameter(param_addr, "ns_trb_contraction_factor", &tmp_double, &tmp_res, 0.25, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("Trust Region Based").set("Contraction Factor",tmp_double);

	  // "Expansion Trigger Ratio" ($\rho_{\rm e}$) - If the improvement ratio is greater than this value, then the trust region is 
	  // contracted by the amount specified by the "Expansion Factor". Defaults to 0.75.
	  get_double_parameter(param_addr, "ns_trb_expansion_trigger_ratio", &tmp_double, &tmp_res, 0.75, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("Trust Region Based").set("Expansion Trigger Ratio",tmp_double);

	  // "Expansion Factor" ($\beta_{\rm e}$) - See above. Defaults to 4.0.
	  get_double_parameter(param_addr, "ns_trb_expansion_factor", &tmp_double, &tmp_res, 4.0, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("Trust Region Based").set("Expansion Factor",tmp_double);

	  // "Recovery Step" - Defaults to 1.0.
	  get_double_parameter(param_addr, "ns_trb_recovery_step", &tmp_double, &tmp_res, 1.0, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("Trust Region Based").set("Recovery Step",tmp_double);

	  // "Use Ared/Pred Ratio Calculation" (boolean) - Defaults to false. 
	  // If set to true, this option replaces the algorithm used to compute 
	  // the improvement ratio, $ \rho $, as described above. The improvement ratio is replaced by an "Ared/Pred" sufficient decrease
	  // criteria similar to that used in line search algorithms 
	  // (see Eisenstat and Walker, SIAM Journal on Optimization V4 no. 2 (1994) pp 393-422):
	  get_int_parameter(param_addr, "ns_trb_use_ared_pred_ratio_calculation", &tmp_int, &tmp_res, 0, CHECK_VALUES, 2, 0, 1);
	  if (tmp_res!=-1)
	    {
	      if (tmp_int) nlParams.sublist("Trust Region Based").set("Use Ared/Pred Ratio Calculation",true);
	      else         nlParams.sublist("Trust Region Based").set("Use Ared/Pred Ratio Calculation",false);
	    }

	  break;
	case 3:
	  nlParams.set("Nonlinear Solver", "Inexact Trust Region Based");
	  // "Inner Iteration Method" - Choice of trust region algorithm to use. Choices are:
	  // o 1 - "Standard Trust Region"
          // o 2 - "Inexact Trust Region"
	  get_int_parameter(param_addr, "ns_itrb_inner_iteration_method", &tmp_int, &tmp_res, 1, CHECK_VALUES, 2, 1, 2);
	  if (tmp_res!=-1)
	    {
	      switch(tmp_int)
		{
		case 1:
		  nlParams.sublist("Inexact Trust Region Based").set("Inner Iteration Method","Standard Trust Region");
		  break;
		case 2:
		  nlParams.sublist("Inexact Trust Region Based").set("Inner Iteration Method","Inexact Trust Region");
		  break;
		}
	    }

	  // "Minimum Trust Region Radius" ($\Delta_{\min}$) - Minimum allowable trust region radius. Defaults to 1.0e-6.
	  get_double_parameter(param_addr, "ns_itrb_minimum_trust_region_radius", &tmp_double, &tmp_res, 1.0e-6, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("Inexact Trust Region Based").set("Minimum Trust Region Radius",tmp_double);

	  // "Maximum Trust Region Radius" ($\Delta_{\max}$) - Minimum allowable trust region radius. Defaults to 1.0e+10.
	  get_double_parameter(param_addr, "ns_itrb_maximum_trust_region_radius", &tmp_double, &tmp_res, 1.0e10, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("Inexact Trust Region Based").set("Maximum Trust Region Radius",tmp_double);

	  // "Minimum Improvement Ratio" ($\rho_{\min}$) - Minimum improvement ratio to accept the step. Defaults to 1.0e-4.
	  get_double_parameter(param_addr, "ns_itrb_minimum_improvement_ratio", &tmp_double, &tmp_res, 1.0e-4, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("Inexact Trust Region Based").set("Minimum Improvement Ratio",tmp_double);

	  // "Contraction Trigger Ratio" ($\rho_{\rm s}$) - If the improvement ratio is less than this value, then the trust region is
	  // contracted by the amount specified by the "Contraction Factor". Must be larger than "Minimum Improvement Ratio".
	  // Defaults to 0.1.
	  get_double_parameter(param_addr, "ns_itrb_contraction_trigger_ratio", &tmp_double, &tmp_res, 0.1, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("Inexact Trust Region Based").set("Contraction Trigger Ratio",tmp_double);

	  // "Contraction Factor" ($\beta_{\rm s}$) - See above. Defaults to 0.25.
	  get_double_parameter(param_addr, "ns_itrb_contraction_factor", &tmp_double, &tmp_res, 0.25, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("Inexact Trust Region Based").set("Contraction Factor",tmp_double);

	  // "Expansion Trigger Ratio" ($\rho_{\rm e}$) - If the improvement ratio is greater than this value, then the trust region is 
	  // contracted by the amount specified by the "Expansion Factor". Defaults to 0.75.
	  get_double_parameter(param_addr, "ns_itrb_expansion_trigger_ratio", &tmp_double, &tmp_res, 0.75, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("Inexact Trust Region Based").set("Expansion Trigger Ratio",tmp_double);

	  // "Expansion Factor" ($\beta_{\rm e}$) - See above. Defaults to 4.0.
	  get_double_parameter(param_addr, "ns_itrb_expansion_factor", &tmp_double, &tmp_res, 4.0, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("Inexact Trust Region Based").set("Expansion Factor",tmp_double);

	  // "Recovery Step" - Defaults to 1.0.
	  get_double_parameter(param_addr, "ns_itrb_recovery_step", &tmp_double, &tmp_res, 1.0, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("Inexact Trust Region Based").set("Recovery Step",tmp_double);

	  // "Use Ared/Pred Ratio Calculation" (boolean) - Defaults to false. If set to true, this option replaces the algorithm used to
	  // compute the improvement ratio, $ \rho $, as described above. 
	  // The improvement ratio is replaced by an "Ared/Pred" sufficient decrease
	  // criteria similar to that used in line search algorithms (see Eisenstat and Walker, SIAM Journal on 
	  // Optimization V4 no. 2 (1994) pp 393-422):
	  get_int_parameter(param_addr, "ns_itrb_use_ared_pred_ratio_calculation", &tmp_int, &tmp_res, 0, CHECK_VALUES, 2, 0, 1);
	  if (tmp_res!=-1)
	    {
	      if (tmp_int) nlParams.sublist("Inexact Trust Region Based").set("Use Ared/Pred Ratio Calculation",true);
	      else         nlParams.sublist("Inexact Trust Region Based").set("Use Ared/Pred Ratio Calculation",false);
	    }

	  // "Use Cauchy in Newton Direction" - Boolean. Used only by the "Inexact Trust Region" algorithm. 
	  // If set to true, the initial guess
	  // for the Newton direction computation will use the Cauchy direction as the initial guess. Defaults to false.
	  get_int_parameter(param_addr, "ns_itrb_use_cauchy_in_newton_direction", &tmp_int, &tmp_res, 0, CHECK_VALUES, 2, 0, 1);
	  if (tmp_res!=-1)
	    {
	      if (tmp_int) nlParams.sublist("Inexact Trust Region Based").set("Use Cauchy in Newton Direction",true);
	      else         nlParams.sublist("Inexact Trust Region Based").set("Use Cauchy in Newton Direction",false);
	    }

	  // "Use Dogleg Segment Minimization" - Boolean. Used only by the "Inexact Trust Region" algorithm. If set to true, 
	  // the $ \tau $ parameter is minimized over the dogleg line segments instead of being computed at the trust regioin radius. 
	  // Used only by the "Inexact Trust Region" algorithm. Defaults to false.
	  get_int_parameter(param_addr, "ns_itrb_use_dogleg_segment_minimization", &tmp_int, &tmp_res, 0, CHECK_VALUES, 2, 0, 1);
	  if (tmp_res!=-1)
	    {
	      if (tmp_int) nlParams.sublist("Inexact Trust Region Based").set("Use Dogleg Segment Minimization",true);
	      else         nlParams.sublist("Inexact Trust Region Based").set("Use Dogleg Segment Minimization",false);
	    }

	  break;
	case 4:
	  nlParams.set("Nonlinear Solver", "Tensor Based");
	  // "Direction" - Sublist of the direction parameters, passed to the NOX::Direction::Factory constructor. Defaults to an empty list.
	  // * "Method" - Name of the direction to be computed in this solver. "Tensor" and "Newton" are the only two valid choices. 
	  //   A sublist by this name specifies all of the parameters to be passed to the linear solver. See below under "Linear Solver".
	  // - 1 - "Tensor"
	  // - 2 - "Newton"
	  get_int_parameter(param_addr, "ns_tb_direction_method", &tmp_int, &tmp_res, 2, CHECK_VALUES, 2, 1, 2);
	  if (tmp_res!=-1)
	    {
	      switch(tmp_int)
		{
		case 1:
		  nlParams.sublist("Tensor Based").sublist("Direction").set("Method","Tensor");
		  break;
		case 2:
		  nlParams.sublist("Tensor Based").sublist("Direction").set("Method","Newton");
		  break;
		}
	    }

	  // * "Rescue Bad Newton Solve" (Boolean) - If the linear solve does not meet the tolerance specified by the forcing term, 
	  //   then use the step anyway. Defaults to true.
	  get_int_parameter(param_addr, "ns_tb_rescue_bad_newton_solve", &tmp_int, &tmp_res, 1, CHECK_VALUES, 2, 0, 1);
	  if (tmp_res!=-1)
	    {
	      if (tmp_int) nlParams.sublist("Tensor Based").set("Rescue Bad Newton Solve",true);
	      else         nlParams.sublist("Tensor Based").set("Rescue Bad Newton Solve",false);
	    }

	  // * "Linear Solver" - Sublist for the specific linear solver parameters that are passed to NOX::Abstract::Group::computeNewton() 
	  //   and NOX::Abstract::Group::applyJacobianInverse(). "Linear Solver" is itself a sublist of the list specified in "Method" 
	  //   above (i.e., "Tensor" or "Newton"). Below is a partial list of standard parameters usually available in common linear solvers. 
	  //   Check with the specific linear solver being used for other parameters.
	  //   o "Max Iterations" - Maximum number of Arnoldi iterations (also max Krylov space dimension)
	  //   o "Tolerance" - Relative tolerance for solving local model [default = 1e-4]
	  //   o "Output Frequency" - Print output at every number of iterations [default = 20]
	  get_int_parameter(param_addr, "ns_tb_linear_solver_max_iterations", &tmp_int, &tmp_res, 300, CHECK_MIN, 0);
	  if (tmp_res!=-1) nlParams.sublist("Tensor Based").sublist("Linear Solver").set("Max Iterations",tmp_int);

	  get_double_parameter(param_addr, "ns_tb_linear_solver_tolerance", &tmp_double, &tmp_res, 1e-4, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("Tensor Based").sublist("Linear Solver").set("Tolerance",tmp_double);

	  get_int_parameter(param_addr, "ns_tb_linear_solver_output_frequency", &tmp_int, &tmp_res, 20, CHECK_MIN, 0);
	  if (tmp_res!=-1) nlParams.sublist("Tensor Based").sublist("Linear Solver").set("Output Frequency",tmp_int);

	  // "Line Search" - Sublist of the line search parameters. Because the tensor step is not guaranteed to be a descent direction on 
	  // the function, not all "basic" line search approaches would be appropriate. Thus, the LineSearch classes available to Newton's
	  // method (e.g., Polynomial, More-Thuente) are not used here. Instead, this solver class approriately handles technical 
	  // considerations
	  // for tensor methods with its own set of global strategies. The following parameters specify the specific options
	  // for this line search:
	  // * "Method" - Name of the line search available to tensor methods Valid choices are:
	  //   o 1 - "Curvilinear" - Backtrack along the "curvilinear" path that spans the tensor direction and the Newton direction and
	  //          that maintains monotonicity on the tensor model. Recommended because it tends to be more robust and efficient
	  //          than the other choices. [Default]
	  //   o 2 - "Standard" - Backtrack along tensor direction unless it is not a descent direction, in which case backtrack 
	  //          along Newton direction.
	  //   o 3 - "Dual" - Backtrack along both the Newton and tensor directions and choose the better of the two.
	  //   o 4 - "Full Step" - Only use the full step and do not backtrack along both the Newton and tensor directions and
	  //          choose the better of the two.
	  get_int_parameter(param_addr, "ns_tb_line_search_method", &tmp_int, &tmp_res, 1, CHECK_VALUES, 4, 1, 2, 3, 4);
	  if (tmp_res!=-1)
	    {
	      switch(tmp_int)
		{
		case 1:
		  nlParams.sublist("Tensor Based").sublist("Line Search").set("Method","Curvilinear"); // Default
		  break;
		case 2:
		  nlParams.sublist("Tensor Based").sublist("Line Search").set("Method","Standard");
		  break;
		case 3:
		  nlParams.sublist("Tensor Based").sublist("Line Search").set("Method","Dual");
		  break;
		case 4:
		  nlParams.sublist("Tensor Based").sublist("Line Search").set("Method","Full Step");
		  break;
		}
	    }

	  // * "Lambda selection" - Flag for how to calculate the next linesearch parameter lambda. Valid choices are "Quadratic" and 
	  //   "Halving" (default). Quadratic constructs a quadratic interpolating polynomial from the last trial point and uses
	  //   the minimum of this function as the next trial lambda (bounded by 0.1). Halving divides the linesearch parameter
	  //   by 2 before each trial, which is simpler but tends to generate longer steps than quadratic.
	  get_int_parameter(param_addr, "ns_tb_line_search_lambda_selection", &tmp_int, &tmp_res, 2, CHECK_VALUES, 2, 1, 2);
	  if (tmp_res!=-1)
	    {
	      switch(tmp_int)
		{
		case 1:
		  nlParams.sublist("Tensor Based").sublist("Line Search").set("Lambda Selection","Quadratic");
		  break;
		case 2:
		  nlParams.sublist("Tensor Based").sublist("Line Search").set("Lambda Selection","Halving"); // default
		  break;
		}
	    }

	  // * "Default Step" - Starting value of the linesearch parameter (defaults to 1.0)
	  get_double_parameter(param_addr, "ns_tb_default_step", &tmp_double, &tmp_res, 1.0, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("Tensor Based").set("Default Step",tmp_double);

	  // * "Minimum Step" - Minimum acceptable linesearch parameter before the linesearch terminates (defaults to 1.0e-12). 
	  //   If there are many linesearch failures, then lowering this value is one thing to try.
	  get_double_parameter(param_addr, "ns_tb_minimum_step", &tmp_double, &tmp_res, 1.0e-12, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("Tensor Based").set("Minimum Step",tmp_double);

	  // * "Recovery Step Type" - Determines the step size to take when the line search fails. Choices are:
	  //   o 1 - "Constant" [default] - Uses a constant value set in "Recovery Step".
	  //   o 2 - "Last Computed Step" - Uses the last value computed by the line search algorithm.
	  get_int_parameter(param_addr, "ns_tb_recovery_step_type", &tmp_int, &tmp_res, 1, CHECK_VALUES, 2, 1, 2);
	  if (tmp_res!=-1)
	    {
	      switch(tmp_int)
		{
		case 1:
		  nlParams.sublist("Tensor Based").set("Recovery Step Type","Constant"); // Default
		  break;
		case 2:
		  nlParams.sublist("Tensor Based").set("Recovery Step Type","Lag Computed Step");
		  break;
		}
	    }

	  // * "Recovery Step" - Step parameter to take when the line search fails (defaults to value for "Default Step")
	  // - 1 - Default Step
	  // - 2 - Minimum Step
	  // - 3 - Recovery Step
	  get_int_parameter(param_addr, "ns_tb_recovery_step", &tmp_int, &tmp_res, 1, CHECK_VALUES, 3, 1, 2, 3);
	  if (tmp_res!=-1)
	    {
	      switch(tmp_int)
		{
		case 1:
		  nlParams.sublist("Tensor Based").set("Recovery Step","Default Step");
		  break;
		case 2:
		  nlParams.sublist("Tensor Based").set("Recovery Step","Minimum Step");
		  break;
		case 3:
		  nlParams.sublist("Tensor Based").set("Recovery Step","Recovery Step");
		  break;
		}
	    }
       
	  // * "Max Iters" - Maximum number of iterations (i.e., backtracks)
	  get_int_parameter(param_addr, "ns_tb_max_iters", &tmp_int, &tmp_res, 300, CHECK_MIN, 0);
	  if (tmp_res!=-1) nlParams.sublist("Tensor Based").set("Max Iters",tmp_int);

	  break;
	}
    }
  // Set the printing parameters in the "Printing" sublist
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", MyPID); 
  printParams.set("Output Precision", 3);
  printParams.set("Output Processor", 0);

  NOX::Utils::MsgType printingDetails = NOX::Utils::Error;
  // Error 	Errors are always printed.
  // Warning                  2^0
  // OuterIteration           2^1
  // InnerIteration           2^2
  // Parameters               2^3
  // Details                  2^4
  // OuterIterationStatusTest 2^5
  // LinearSolverDetails      2^6
  // TestDetails              2^7
  // Debug                    2^12

  get_int_parameter(param_addr, "printing_warning", &tmp_int, &tmp_res, 0, CHECK_VALUES, 2, 0, 1);
  if (tmp_res!=-1) printingDetails = (NOX::Utils::MsgType)(printingDetails + NOX::Utils::Warning);

  get_int_parameter(param_addr, "printing_outeriterations", &tmp_int, &tmp_res, 0, CHECK_VALUES, 2, 0, 1);
  if (tmp_res!=-1) printingDetails = (NOX::Utils::MsgType)(printingDetails + NOX::Utils::OuterIteration);

  get_int_parameter(param_addr, "printing_inneriterations", &tmp_int, &tmp_res, 0, CHECK_VALUES, 2, 0, 1);
  if (tmp_res!=-1) printingDetails = (NOX::Utils::MsgType)(printingDetails + NOX::Utils::InnerIteration);

  get_int_parameter(param_addr, "printing_parameters", &tmp_int, &tmp_res, 0, CHECK_VALUES, 2, 0, 1);
  if (tmp_res!=-1) printingDetails = (NOX::Utils::MsgType)(printingDetails + NOX::Utils::Parameters);

  get_int_parameter(param_addr, "printing_details", &tmp_int, &tmp_res, 0, CHECK_VALUES, 2, 0, 1);
  if (tmp_res!=-1) printingDetails = (NOX::Utils::MsgType)(printingDetails + NOX::Utils::Details);

  get_int_parameter(param_addr, "printing_outeriterationstatustest", &tmp_int, &tmp_res, 0, CHECK_VALUES, 2, 0, 1);
  if (tmp_res!=-1) printingDetails = (NOX::Utils::MsgType)(printingDetails + NOX::Utils::OuterIterationStatusTest);

  get_int_parameter(param_addr, "printing_linearsolverdetails", &tmp_int, &tmp_res, 0, CHECK_VALUES, 2, 0, 1);
  if (tmp_res!=-1) printingDetails = (NOX::Utils::MsgType)(printingDetails + NOX::Utils::LinearSolverDetails);

  get_int_parameter(param_addr, "printing_testdetails", &tmp_int, &tmp_res, 0, CHECK_VALUES, 2, 0, 1);
  if (tmp_res!=-1) printingDetails = (NOX::Utils::MsgType)(printingDetails + NOX::Utils::TestDetails);

  get_int_parameter(param_addr, "printing_debug", &tmp_int, &tmp_res, 0, CHECK_VALUES, 2, 0, 1);
  if (tmp_res!=-1) printingDetails = (NOX::Utils::MsgType)(printingDetails + NOX::Utils::Debug);

  printParams.set("Output Information", printingDetails);

  // Redirect std::cerr and std::cout to Scilab
  printParams.set("outputStream", scicout); 
  printParams.set("errorStream",  scicerr); 

  // Create printing utilities
  NOX::Utils utils(printParams);
  
  // Sublist for line search 
  Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
  // 'ls_method':
  // - 1 "BackTrack"
  // - 2 "Full Step"
  // - 3 "Polynomial"
  // - 4 "NonLinearCG"
  // - 5 "More'-Thuente"

  get_int_parameter(param_addr, "ls_method", &tmp_int, &tmp_res, 5, CHECK_VALUES, 5, 1, 2, 3, 4, 5);
  if (tmp_res!=-1)
    {
      switch(tmp_int)
	{
	case 1:
	  searchParams.set("Method", "BackTrack");

	  // starting step length (defaults to 1.0) 
	  get_double_parameter(param_addr, "ls_bt_default_step", &tmp_double, &tmp_res, 1.0, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("BackTrack").set("Default Step",tmp_double);

	  // minimum acceptable step length (defaults to 1.0e-12) 
	  get_double_parameter(param_addr, "ls_bt_minimum_step", &tmp_double, &tmp_res, 1.0e-12, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("BackTrack").set("Minimum Step",tmp_double);

	  // step to take when the line search fails (defaults to value for "Default Step") 
	  get_double_parameter(param_addr, "ls_bt_recovery_step", &tmp_double, &tmp_res, 1.0, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("BackTrack").set("Recovery Step",tmp_double);

	  // maximum number of iterations (i.e., RHS computations) 
	  get_int_parameter(param_addr, "ls_bt_max_iters", &tmp_int, &tmp_res, 20, CHECK_MIN, 0);
	  if (tmp_res!=-1) nlParams.sublist("BackTrack").set("Max Iters",tmp_int);

	  get_int_parameter(param_addr, "ls_bt_decrease_condition", &tmp_int, &tmp_res, 1, CHECK_VALUES, 2, 1, 2);
	  if (tmp_res!=-1)
	    {
	      switch(tmp_int)
		{
		case 1:
		  // simple decrease using the infinity norm of the RHS (default) 
		  nlParams.sublist("BackTrack").set("Decrease Condition","Max Norm");
		  break;
		case 2:
		  // simple decrease using the Euclidean norm of the RHS 
		  nlParams.sublist("BackTrack").set("Decrease Condition","Two Norm");
		  break;
		}
	    }

	  // A multiplier between zero and one that reduces the step size between line search iterations
	  get_double_parameter(param_addr, "ls_bt_reduction_factor", &tmp_double, &tmp_res, 0.5, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("BackTrack").set("Reduction Factor",tmp_double);

	  break;
	case 2:
	  searchParams.set("Method", "Full Step");

	  // length of a full step (defaults to 1.0) 
	  get_double_parameter(param_addr, "ls_fs_full_step", &tmp_double, &tmp_res, 1.0, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("Full Step").set("Full Step",tmp_double);
	  break;
	case 3:
	  searchParams.set("Method", "Polynomial");

	  // Starting step length, i.e., $\lambda_0$. Defaults to 1.0.
	  get_double_parameter(param_addr, "ls_pol_default_step", &tmp_double, &tmp_res, 1.0, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("Polynomial").set("Default Step",tmp_double);

	  // Maximum number of line search iterations. The search fails if the number of iterations exceeds this value. Defaults to 100.
	  get_int_parameter(param_addr, "ls_pol_max_iters", &tmp_int, &tmp_res, 100, CHECK_MIN, 0);
	  if (tmp_res!=-1) nlParams.sublist("Polynomial").set("Max Iters",tmp_int);

	  // Minimum acceptable step length. The search fails if the computed $\lambda_k$ is less than this value. Defaults to 1.0e-12.
	  get_double_parameter(param_addr, "ls_pol_minimum_step", &tmp_double, &tmp_res, 1.0e-12, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("Polynomial").set("Minimum Step",tmp_double);

	  // 'pol_recovery_type': Determines the step size to take when the line search fails. Choices are:
	  // - 1 "Constant" - Uses a constant value set in "Recovery Step". 
	  // - 2 "Last Computed Step" - Uses the last value computed by the line search algorithm. 
	  get_int_parameter(param_addr, "ls_pol_recovery_step_type", &tmp_int, &tmp_res, 1, CHECK_VALUES, 2, 1, 2);
	  if (tmp_res!=-1) 
	    {
	      switch(tmp_int)
		{
		case 1:
		  nlParams.sublist("Polynomial").set("Recovery Step Type","Constant");
		  break;
		case 2:
		  nlParams.sublist("Polynomial").set("Recovery Step Type","Last Computed Step");
		  break;
		}
	    }

	  // The value of the step to take when the line search fails. Only used if the "Recovery Step Type" is set to "Constant". 
	  // Defaults to value for "Default Step".
	  get_double_parameter(param_addr, "ls_pol_recovery_step", &tmp_double, &tmp_res, 1.0, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("Polynomial").set("Recovery Step",tmp_double);

	  // Type of interpolation that should be used.
	  get_int_parameter(param_addr, "ls_pol_interpolation_type", &tmp_int, &tmp_res, 1, CHECK_VALUES, 3, 1, 2, 3);
	  if (tmp_res!=-1) 
	    {
	      switch(tmp_int)
		{
		case 1:
		  nlParams.sublist("Polynomial").set("Interpolation Type","Cubic");
		  break;
		case 2:
		  nlParams.sublist("Polynomial").set("Interpolation Type","Quadratic");
		  break;
		case 3:
		  nlParams.sublist("Polynomial").set("Interpolation Type","Quadratic3");
		  break;
		}
	    }

	  // Choice for $ \gamma_{min} $, i.e., the factor that limits the minimum size of the new step based on the previous step. 
	  // Defaults to 0.1.
	  get_double_parameter(param_addr, "ls_pol_min_bounds_factor", &tmp_double, &tmp_res, 0.1, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("Polynomial").set("Min Bounds Factor",tmp_double);

	  // Choice for $ \gamma_{max} $, i.e., the factor that limits the maximum size of the new step based on the previous step.
	  // Defaults to 0.5.
	  get_double_parameter(param_addr, "ls_pol_max_bounds_factor", &tmp_double, &tmp_res, 0.5, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("Polynomial").set("Max Bounds Factor",0.5);

	  // The decrease condition
	  get_int_parameter(param_addr, "ls_pol_sufficient_decrease_condition", &tmp_int, &tmp_res, 2, CHECK_VALUES, 3, 1, 2, 3);
	  if (tmp_res!=-1)
	    {
	      switch(tmp_int)
		{
		case 1:
		  nlParams.sublist("Polynomial").set("Sufficient Decrease Condition","None");
		  break;
		case 2:
		  nlParams.sublist("Polynomial").set("Sufficient Decrease Condition","Armijo-Goldstein");
		  break;
		case 3:
		  nlParams.sublist("Polynomial").set("Sufficient Decrease Condition","Ared/Pred");
		  break;
		}
	    }

	  // Parameter choice for sufficient decrease condition. See checkConvergence() for details. Defaults to 1.0e-4.
	  get_double_parameter(param_addr, "ls_pol_alpha_factor", &tmp_double, &tmp_res, 1.0e-4, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("Polynomial").set("Alpha Factor",tmp_double);

	  // Set to true if at least one interpolation step should be used. The default is false which means that the line search will stop 
	  // if the default step length satisfies the convergence criteria. Defaults to false.
	  get_int_parameter(param_addr, "ls_pol_force_interpolation", &tmp_int, &tmp_res, 0, CHECK_VALUES, 2, 0, 1);
	  if (tmp_res!=-1) 
	    {
	      if (tmp_int) nlParams.sublist("Polynomial").set("Force Interpolation",true);
	      else         nlParams.sublist("Polynomial").set("Force Interpolation",false);
	    }

	  // Maximum index of the nonlinear iteration for which we allow a relative increase. See checkConvergence() 
	  // for further details. Defaults to 0 (zero).
 	  get_int_parameter(param_addr, "ls_pol_maximum_iteration_for_increase", &tmp_int, &tmp_res, 0, CHECK_MIN, 0);
	  if (tmp_res!=-1) nlParams.sublist("Polynomial").set("Maximum Iteration for Increase",tmp_int);

	  // See checkConvergence() for details. Defaults to 100.
	  get_int_parameter(param_addr, "ls_pol_allowed_relative_increase", &tmp_int, &tmp_res, 100, CHECK_MIN, 0);
	  if (tmp_res!=-1) nlParams.sublist("Polynomial").set("Allowed Relative Increase",tmp_int);

	  break;
	case 4:
	  searchParams.set("Method", "NonlinearCG");

	  break;
	case 5:
	  searchParams.set("Method", "More'-Thuente");

	  // Choice to use for the sufficient decrease condition. Options are "Ared/Pred" or "Armijo-Goldstein" 
	  // (defaults to "Armijo-Goldstein").
	  // 1. "Armijo-Goldstein" conditions: $ f(x_{n-1}+ \lambda s) \le f(x_{n-1}) +\alpha \lambda f'(x_{n-1}) $
	  // 2. "Ared/Pred" conditions: 
	  //    $ \| F(x_{n-1}+ \lambda s) \| \le \| F(x_{n-1}) \| (1-\alpha(1-\eta)) 
	  //    $ where $ \eta $ is the linear solve tolerance in the inexact Newton method. 
	  get_int_parameter(param_addr, "ls_mt_sufficient_decrease_condition", &tmp_int, &tmp_res, 1, CHECK_VALUES, 2, 1, 2);
	  if (tmp_res!=-1)
	    {
	      switch(tmp_int)
		{
		case 1:
		  nlParams.sublist("More'-Thuente").set("Sufficient Decrease Condition","Armijo-Goldstein");
		  break;
		case 2:
		  nlParams.sublist("More'-Thuente").set("Sufficient Decrease Condition","Ared/Pred");
		  break;
		}
	    }

	  // The ftol in the sufficient decrease condition (defaults to 1.0e-4)
	  get_double_parameter(param_addr, "ls_mt_sufficient_decrease", &tmp_double, &tmp_res, 1.0e-4, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("More'-Thuente").set("Sufficient Decrease",tmp_double);

	  // The gtol in the curvature condition (defaults to 0.9999)
	  get_double_parameter(param_addr, "ls_mt_curvature_condition", &tmp_double, &tmp_res, 0.9999, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("More'-Thuente").set("Curvature Condition",tmp_double);

	  // If set to true the value of $ s^TJ^TF $ is estimated using a directional derivative in a call to 
	  // NOX::LineSearch::Common::computeSlopeWithOutJac. If false the slope computation is computed with the 
	  // NOX::LineSearch::Common::computeSlope method. Setting this to true eliminates having to compute the Jacobian at each
	  // inner iteration of the More'-Thuente line search (defaults to false).
 	  get_int_parameter(param_addr, "ls_mt_optimize_slope_computation", &tmp_int, &tmp_res, 0, CHECK_VALUES, 2, 0, 1);
	  if (tmp_res!=-1)
	    {
	      if (tmp_int) nlParams.sublist("More'-Thuente").set("Optimize Slope Calculation",true);
	      else         nlParams.sublist("More'-Thuente").set("Optimize Slope Calculation",false);
	    }

	  // The maximum width of the interval containing the minimum of the modified function (defaults to 1.0e-15)
	  get_double_parameter(param_addr, "ls_mt_interval_width", &tmp_double, &tmp_res, 1.0e-15, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("More'-Thuente").set("Interval Width",tmp_double);

	  // maximum allowable step length (defaults to 1.0e6)
	  get_double_parameter(param_addr, "ls_mt_maximum_step", &tmp_double, &tmp_res, 1.0e6, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("More'-Thuente").set("Maximum Step",tmp_double);

	  // minimum allowable step length (defaults to 1.0e-12)
	  get_double_parameter(param_addr, "ls_mt_minimum_step", &tmp_double, &tmp_res, 1.0e-12, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("More'-Thuente").set("Minimum Step",tmp_double);

	  // maximum number of right-hand-side and corresponding Jacobian evaluations (defaults to 20)
 	  get_int_parameter(param_addr, "ls_mt_max_iters", &tmp_int, &tmp_res, 20, CHECK_MIN, 0);
	  if (tmp_res!=-1) nlParams.sublist("More'-Thuente").set("Max Iters",tmp_int);

	  // starting step length (defaults to 1.0)
	  get_double_parameter(param_addr, "ls_mt_default_step", &tmp_double, &tmp_res, 1.0, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("More'-Thuente").set("Default Step",tmp_double);

	  // Determines the step size to take when the line search fails. Choices are:
	  // * "Constant" [default] - Uses a constant value set in "Recovery Step".
	  // * "Last Computed Step" - Uses the last value computed by the line search algorithm.
	  get_int_parameter(param_addr, "ls_mt_recovery_step_type", &tmp_int, &tmp_res, 1, CHECK_VALUES, 2, 1, 2);
	  if (tmp_res!=-1)
	    {
	      switch(tmp_int)
		{
		case 1:
		  nlParams.sublist("More'-Thuente").set("Recovery Step Type","Constant");
		  break;
		case 2:
		  nlParams.sublist("More'-Thuente").set("Recovery Step Type","Last Computed Step");
		  break;
		}
	    }

	  // The value of the step to take when the line search fails. Only used if the "Recovery Step Type" is set to "Constant". 
	  // Defaults to value for "Default Step".
	  get_double_parameter(param_addr, "ls_mt_recovery_step", &tmp_double, &tmp_res, 1.0, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) nlParams.sublist("More'-Thuente").set("Recovery Step",tmp_double);

	  break;
	}
    }

  // Sublist for direction
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  
  // 'dir_method':
  // 1 - Broyden
  // 2 - Newton
  // 3 - Steepest Descent
  // 4 - Nonlinear CG
  get_int_parameter(param_addr, "dir_method", &tmp_int, &tmp_res, 1, CHECK_VALUES, 4, 1, 2, 3, 4);
  if (tmp_res!=-1)
    {
      switch(tmp_int)
	{
	case 1:
	  // Broyden
	  dirParams.set("Method", "Broyden");

	  // How often the Jacobian should be refreshed. A value of 5, for example, means that the Jacobian should
	  // be updated every 5 iterations. Defaults to 10.
 	  get_int_parameter(param_addr, "dir_broyden_restart_frequency", &tmp_int, &tmp_res, 10, CHECK_MIN, 0);
	  if (tmp_res!=-1) dirParams.sublist("Broyden").set("Restart Frequency", tmp_int);

	  // Maximum convergence rate allowed when reusing the Jacobian. The Jacobian will be refreshed if the convergence 
	  // rate, $ \alpha $, is larger than this value. The convergence rate is calculated by 
	  // $ \alpha = \frac{\| F_k \| }{\| F_{k-1} \|} 
	  // $ where F is the nonlinear residual and $ k $ is the nonlinear iteration. Defaults to 1.0.
 	  get_double_parameter(param_addr, "dir_broyden_max_convergence_rate", &tmp_double, &tmp_res, 1.0, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) dirParams.sublist("Broyden").set("Max Convergence Rate", tmp_double);

	  // The maximum number of past updates that can be saved in memory. Defaults to the value of "Restart Frequency".
 	  get_int_parameter(param_addr, "dir_broyden_memory", &tmp_int, &tmp_res, 10, CHECK_MIN, 0);
	  if (tmp_res!=-1) dirParams.sublist("Broyden").set("Memory", tmp_int);

	  break;
	case 2:
	  // Newton
	  dirParams.set("Method", "Newton");

	  // Method to compute the forcing term, i.e., the tolerance for the linear solver.
	  // see http://trilinos.sandia.gov/packages/docs/r4.0/packages/nox/doc/html/classNOX_1_1Direction_1_1Newton.html
 	  get_int_parameter(param_addr, "dir_newton_forcing_term_method", &tmp_int, &tmp_res, 1, CHECK_VALUES, 3, 1, 2, 3);
	  if (tmp_res!=-1)
	    {
	      switch(tmp_int)
		{
		case 1:
		  dirParams.sublist("Newton").set("Forcing Term Method", "Constant");
		  break;
		case 2:
		  dirParams.sublist("Newton").set("Forcing Term Method", "Type 1");
		  break;
		case 3:
		  dirParams.sublist("Newton").set("Forcing Term Method", "Type 2");
		  break;
		}
	    }

	  // $\eta_0$ (initial linear solver tolerance). Defaults to 0.1
 	  get_double_parameter(param_addr, "dir_newton_forcing_term_initial_tolerance", &tmp_double, &tmp_res, 0.1, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) dirParams.sublist("Newton").set("Forcing Term Initial Tolerance", tmp_double);

	  // $\eta_{\min}$. Defaults to 1.0e-6.
 	  get_double_parameter(param_addr, "dir_newton_forcing_term_minimum_tolerance", &tmp_double, &tmp_res, 1.0e-6, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) dirParams.sublist("Newton").set("Forcing Term Minimum Tolerance", tmp_double);

	  // $\eta_{\max}$. Defaults to 0.01.
 	  get_double_parameter(param_addr, "dir_newton_forcing_term_maximum_tolerance", &tmp_double, &tmp_res, 0.01, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) dirParams.sublist("Newton").set("Forcing Term Maximum Tolerance", tmp_double);

	  // $\alpha$ (used only by "Type 2"). Defaults to 1.5.
 	  get_double_parameter(param_addr, "dir_newton_forcing_term_alpha", &tmp_double, &tmp_res, 1.5, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) dirParams.sublist("Newton").set("Forcing Term Alpha", tmp_double);

	  // $\gamma$ (used only by "Type 2"). Defaults to 0.9
 	  get_double_parameter(param_addr, "dir_newton_forcing_term_gamma", &tmp_double, &tmp_res, 1.5, CHECK_MIN, 0.0);
	  if (tmp_res!=-1) dirParams.sublist("Newton").set("Forcing Term Gamma", tmp_double);

	  break;
	case 3:
	  // Steepest Descent
	  dirParams.set("Method", "Steepest Descent");

	  // see http://trilinos.sandia.gov/packages/docs/r6.0/packages/nox/doc/html/classNOX_1_1Direction_1_1SteepestDescent.html
 	  get_int_parameter(param_addr, "dir_newton_forcing_term_method", &tmp_int, &tmp_res, 1, CHECK_VALUES, 4, 1, 2, 3, 4);
	  if (tmp_res!=-1)
	    {
	      switch(tmp_int)
		{
		case 1:
		  dirParams.sublist("Steepest Descent").set("Scaling Type", "None");
		  break;
		case 2:
		  dirParams.sublist("Steepest Descent").set("Scaling Type", "2-Norm");
		  break;
		case 3:
		  dirParams.sublist("Steepest Descent").set("Scaling Type", "F 2-Norm");
		  break;
		case 4:
		  dirParams.sublist("Steepest Descent").set("Scaling Type", "Quadratic Model Min");
		  break;
		}
	    }

	  break;
	case 4:
	  // Nonlinear CG
	  dirParams.set("Method", "NonlinearCG");

	  // An integer specification of the number of nonlinear iterations between restarts [default = 10]. 
	  // Restart corresponds to setting $\beta = 0$. A good heuristic is to limit this value to the number of problem 
	  // degrees of freedom. 
	  // Setting this value to 1 forces $ \beta = 0 $ for every nonlinear iteration which corresponds to suppressing orthogonalization 
	  // against the previous search direction.
 	  get_int_parameter(param_addr, "dir_nonlinear_cg_restart_frequency", &tmp_int, &tmp_res, 10, CHECK_MIN, 0);
	  if (tmp_res!=-1) dirParams.sublist("Nonlinear CG").set("Restart Frequency", tmp_int);

	  // can be either "On" or "Off" [default]: determines whether or not to compute and apply preconditioner $ M $. 
	  // If "Off" is selected, no preconditioner is computed and the behavior is equivalent to $ M = I $ where $ I $ is 
	  // the identity matrix. 
	  // If "On", $ M $ is computed and applied as determined by the underlying implementation of the
	  // "applyRightPreconditioning" method in the Group.
 	  get_int_parameter(param_addr, "dir_nonlinear_cg_precondition", &tmp_int, &tmp_res, 1, CHECK_VALUES, 2, 0, 1);
	  if (tmp_res!=-1)
	    {
	      if (tmp_int) dirParams.sublist("Nonlinear CG").set("Precondition", "On");
	      else         dirParams.sublist("Nonlinear CG").set("Precondition", "Off");
	    }

	  // see http://trilinos.sandia.gov/packages/docs/r7.0/packages/nox/doc/html/classNOX_1_1Direction_1_1NonlinearCG.html
 	  get_int_parameter(param_addr, "dir_nonlinear_cg_orthogonalize", &tmp_int, &tmp_res, 1, CHECK_VALUES, 2, 1, 2);
	  if (tmp_res!=-1)
	    {
	      switch(tmp_int)
		{
		case 1:
		  dirParams.sublist("Nonlinear CG").set("Orthogonalize", "Polak-Ribiere");
		  break;
		case 2:
		  dirParams.sublist("Nonlinear CG").set("Orthogonalize", "Fletcher-Reeves");
		  break;
		}
	    }
	  break;
	}
    }

  // "Aztec Solver" - Determine the iterative technique used in the solve. The following options are valid:
  // * 1 - "GMRES"    - Restarted generalized minimal residual (default).
  // * 2 - "CG"       - Conjugate gradient.
  // * 3 - "CGS"      - Conjugate gradient squared.
  // * 4 - "TFQMR"    - Transpose-free quasi-minimal reasidual.
  // * 5 - "BiCGStab" - Bi-conjugate gradient with stabilization.
  // * 6 - "LU"       - Sparse direct solve (single processor only).
  get_int_parameter(param_addr, "newt_ls_aztec_solver", &tmp_int, &tmp_res, 1, CHECK_VALUES, 6, 1, 2, 3, 4, 5, 6);
  if (tmp_res!=-1)
    {
      switch(tmp_int)
	{
	case 1:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Aztec Solver", "GMRES");  
	  break;
	case 2:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Aztec Solver", "CG");  
	  break;
	case 3:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Aztec Solver", "CGS");  
	  break;
	case 4:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Aztec Solver", "TFQMR");  
	  break;
	case 5:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Aztec Solver", "BiCGStab");  
	  break;
	case 6:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Aztec Solver", "LU");  
	  break;
	}
    }

  // "Size of Krylov Subspace" - When using restarted GMRES this sets the maximum size of the Krylov subspace (defaults to 300).
  get_int_parameter(param_addr, "newt_ls_size_of_krylov_subspace", &tmp_int, &tmp_res, 300, CHECK_MIN, 0);
  if (tmp_res!=-1) dirParams.sublist("Newton").sublist("Linear Solver").set("Size of Krylov Subspace", tmp_int);  

  // "Orthogonalization" - The orthogonalization routine used for the Gram-Schmidt orthogonalization procedure in Aztec. 
  // The following options are valid:
  // * "Classical" - (default).
  // * "Modified"
  get_int_parameter(param_addr, "newt_ls_orthogonalization", &tmp_int, &tmp_res, 1, CHECK_VALUES, 2, 1, 2);
  if (tmp_res!=-1) 
    {
      if (tmp_int==1) dirParams.sublist("Newton").sublist("Linear Solver").set("Orthogonalization", "Classical");  
      else            dirParams.sublist("Newton").sublist("Linear Solver").set("Orthogonalization", "Modified");  
    }

  // "Convergence Test" - Algorithm used to calculate the residual that is used for determining the convergence of the linear solver. 
  // See the Aztec 2.1 manual for more information. The following options are valid:
  // * 1 - "r0" - (default)
  // * 2 - "rhs"
  // * 3 - "norm"
  // * 4 - "no scaling"
  // * 5 - "sol"
  get_int_parameter(param_addr, "newt_ls_convergence_test", &tmp_int, &tmp_res, 1, CHECK_VALUES, 5, 1, 2, 3, 4, 5);
  if (tmp_res!=-1) 
    {
      switch(tmp_int)
	{
	case 1:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Convergence Test", "r0");  
	  break;
	case 2:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Convergence Test", "rhs");  
	  break;
	case 3:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Convergence Test", "norm");  
	  break;
	case 4:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Convergence Test", "no scaling");  
	  break;
	case 5:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Convergence Test", "sol");  
	  break;
	}
    }

  // "Tolerance" - Tolerance used by AztecOO to determine if an iterative linear solve has converged.
  get_double_parameter(param_addr, "newt_ls_tolerance", &tmp_double, &tmp_res, 1.0e-4, CHECK_MIN, 0.0);
  if (tmp_res!=-1) dirParams.sublist("Newton").sublist("Linear Solver").set("Tolerance", tmp_double);  

  // "Ill-Conditioning Threshold" - If the upper hessenberg matrix during GMRES generates a condition number greater
  // than this parameter value, aztec will exit the linear solve returning the it's current solution. The default is 1.0e11.
  get_double_parameter(param_addr, "newt_ls_ill_conditioning_threshold", &tmp_double, &tmp_res, 1.0e11, CHECK_MIN, 0.0);
  if (tmp_res!=-1) dirParams.sublist("Newton").sublist("Linear Solver").set("Ill-Conditioning Threshold", tmp_double);  

  // "Preconditioner Iterations" - Number of iterations an AztecOO_Operator should take when solving the preconditioner. 
  // This is only used if an AztecOO preconditioner is used and the solver makes a call to NOX::Epetra::Group::applyRightPreconditioning(). 
  // This is NOT a recomended approach.
  get_int_parameter(param_addr, "newt_ls_preconditioner_iterations", &tmp_int, &tmp_res, 400, CHECK_MIN, 0);
  if (tmp_res!=-1) dirParams.sublist("Newton").sublist("Linear Solver").set("Preconditioner Iterations", tmp_int);  

  // "Max Iterations" - maximum number of iterations in the linear solve. Default is 400.
  get_int_parameter(param_addr, "newt_ls_max_iterations", &tmp_int, &tmp_res, 400, CHECK_MIN, 0);
  if (tmp_res!=-1) dirParams.sublist("Newton").sublist("Linear Solver").set("Max Iterations", tmp_int);  

  // "Zero Initial Guess" - boolean. Zero out the initial guess for linear solves performed through applyJacobianInverse calls 
  // (i.e. zero out the result vector before the linear solve). Defaults to false.
  get_int_parameter(param_addr, "newt_ls_zero_initial_guess", &tmp_int, &tmp_res, 0, CHECK_VALUES, 2, 0, 1);
  if (tmp_res!=-1) 
    {
      if (tmp_int) dirParams.sublist("Newton").sublist("Linear Solver").set("Zero Initial Guess", true);  
      else         dirParams.sublist("Newton").sublist("Linear Solver").set("Zero Initial Guess", false);  
    }

  // "Throw Error on Prec Failure" - boolean. If set to true, an exception will be thrown if the preconditioner fails to
  //  initialize or recompute/refactor. If set to false, a warning will br printed if the NOX::Utils::Warning is enabled 
  // in the printing utilities (NOX::Utils). Defaults to true.
  get_int_parameter(param_addr, "newt_ls_throw_error_on_prec_failure", &tmp_int, &tmp_res, 1, CHECK_VALUES, 2, 0, 1);
  if (tmp_res!=-1) 
    {
      if (tmp_int) dirParams.sublist("Newton").sublist("Linear Solver").set("Throw Error on Prec Failure", true);  
      else         dirParams.sublist("Newton").sublist("Linear Solver").set("Throw Error on Prec Failure", false);  
    }

  // "Output Frequency" - number of linear solve iterations between output of the linear solve residual. Takes an integer, or one of 
  // the AztecOO flags: AZ_none, AZ_last, or AZ_all as a value. Defaults to AZ_last.
  get_int_parameter(param_addr, "newt_ls_output_frequency", &tmp_int, &tmp_res, 2, CHECK_VALUES, 3, 1, 2, 3);
  if (tmp_res!=-1) 
    {
      switch(tmp_int)
	{
	case 1:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Output Frequency", AZ_none);  
	  break;
	case 2:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Output Frequency", AZ_last);  
	  break;
	case 3:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Output Frequency", AZ_all);  
	  break;
	}
    }

  // "Jacobian Operator" - If a constructor is used that does not supply a Jacobian operator, nox will create an 
  // internal Jacobian operator. 
  // This flag is ONLY valid in such cases. This will determine which Operator is used:
  // * "Matrix-Free" - Create a NOX::Epetra::MatrixFree object.
  // * "Finite Difference" - Create a NOX::Epetra::FiniteDifference object.
  get_int_parameter(param_addr, "newt_ls_jacobian_operator", &tmp_int, &tmp_res, 1, CHECK_VALUES, 2, 1, 2);
  if (tmp_res!=-1) 
    {
      switch(tmp_int)
	{
	case 1:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Jacobian Operator", "Matrix-Free");  
	  break;
	case 2:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Jacobian Operator", "Finite Difference");  
	  break;
	}
    }

  // "Preconditioner" - Sets the choice of the preconditioner to use during linear solves. The validity of the choice of preconditioner 
  // will depend on the types of operators that are available for the Jacobian and preconditioner. 
  // NOTE: This flag will override any constructor details. For example, if you supply a preconditioner operator in the constructor, 
  // it will not be used if this flag is set to "None". If you supply an Epetra_Operator for the preconditioner but the "Preconditioner" 
  // flag is set to "AztecOO" (this requires an Epetra_RowMatrix for the preconditioner operator), this object will exit with a failure. 
  // The valid options and any requirements on the operator type are listed below:
  // * "None" - No preconditioning. (default)
  // * "AztecOO" - AztecOO internal preconditioner. This requires a preconditioner operator that derives from the Epetra_RowMatrix class.
  // * "Ifpack" - Ifpack internal preconditioner. This requires a preconditioner object that derives from the Epetra_RowMatrix class
  //              or it can use a Jacobian if the Jacobian derives from an Epetra_RowMatrix. This option is deprecated. 
  //              Please use "New Ifpack".
  // * "New Ifpack" - Ifpack internal preconditioner. This requires a preconditioner object that derives from the Epetra_RowMatrix 
  //                  class or it can use a Jacobian if the Jacobian derives from an Epetra_RowMatrix.
  get_int_parameter(param_addr, "newt_ls_preconditioner", &tmp_int, &tmp_res, 1, CHECK_VALUES, 5, 1, 2, 3, 4, 5);
  if (tmp_res!=-1) 
    {
      switch(tmp_int)
	{
	case 1:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Preconditioner", "None");  
	  break;
	case 2:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Preconditioner", "AztecOO");  
	  break;
	case 3:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Preconditioner", "Ifpack");  
	  break;
	case 4:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Preconditioner", "New Ifpack");  
	  break;
	case 5:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Preconditioner", "ML");  
	  break;
	}
    }

  // "Preconditioner Operator" - If a constructor is used that does not supply a preconditioner operator, nox will create an internal
  //  preconditioner operator. This flag is ONLY valid in such cases. This will determine which Operator is used:
  // * "Use Jacobian" - Use the Jacobian Operator (it must be an Epetra_RowMatrix derived object).
  // * "Finite Difference" - Create a NOX::Epetra::FiniteDifference object.
  get_int_parameter(param_addr, "newt_ls_preconditioner_operator", &tmp_int, &tmp_res, 1, CHECK_VALUES, 2, 1, 2);
  if (tmp_res!=-1) 
    {
      switch(tmp_int)
	{
	case 1:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Preconditioner Operator", "Use Jacobian");  
	  break;
	case 2:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Preconditioner Operator", "Finite Difference");  
	  break;
	}
    }

  // "Aztec Preconditioner" - If the "Preconditioner" flag is set to "AztecOO" then the specific AztecOO preconditioner
  // is specified with this flag. 
  // Currently supported preconditioners and their corresponding parameters that can be set are shown below
  // (See the Aztec 2.1 manual for more information):
  // * "ilu" - ilu preconditioning. This choice allows the following additional parameters to be specified:
  //   o "Overlap" - defaults to 0
  //   o "Graph Fill" - defaults to 0
  // * "ilut" - ilut preconditioning. This choice allows the following additional parameters to be specified:
  //   o "Overlap" - defaults to 0
  //   o "Fill Factor" - defaults to 1.0
  //   o "Drop Tolerance" - defaults to 1.0e-12
  // * "Jacobi" - k step Jacobi where k is set by the "Steps" flag:
  //   o "Steps" - defaults to 3.
  // * "Symmetric Gauss-Siedel" - Non-overlapping domain decomposition k step symmetric Gauss-Siedel where k is set by the "Steps" flag:
  //   o "Steps" - defaults to 3.
  // * "Polynomial" - Neumann polynomial with order set by the parameter:
  //   o "Polynomial Order" - defaults to 3.
  // * "Least-squares Polynomial" - Least-squares polynomial with order set by the parameter:
  //   o "Polynomial Order" - defaults to 3.
  get_int_parameter(param_addr, "newt_ls_aztec_preconditioner", &tmp_int, &tmp_res, 1, CHECK_VALUES, 6, 1, 2, 3, 4, 5, 6);
  if (tmp_res!=-1) 
    {
      switch(tmp_int)
	{
	case 1:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Aztec Preconditioner", "ilu");  
	  break;
	case 2:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Aztec Preconditioner", "ilut");  
	  break;
	case 3:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Aztec Preconditioner", "Jacobi");  
	  break;
	case 4:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Aztec Preconditioner", "Symmetric");  
	  break;
	case 5:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Aztec Preconditioner", "Polynomial");  
	  break;
	case 6:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Aztec Preconditioner", "Least-Square Polynomial");  
	  break;
	}
    }

  // "Ifpack" - If the "Preconditioner" flag is set to "New Ifpack" then any of the options supported by the Ifpack Create factory can
  //  be specified using a Teuchos::ParameterList containing the Ifpack options and then setting this as a parameter named "Ifpack" in the
  //  "Linear Solver" sublist.

  // "ML" - If the "Preconditioner" flag is set to "ML" then any of the options supported by the ML factory can be specified using a
  //  Teuchos::ParameterList containing the ML options and then setting this as a parameter named "ML" in the "Linear Solver" sublist.

  // "Preconditioner Reuse Policy" - (string) Allows the user to set how and when the preconditioner should be computed. 
  // This flag supports native Aztec, Ifpack and ML preconditioners. There are three options:
  // * "Rebuild" - The "Rebuild" option always completely destroys and then rebuilds the preconditioner each time a linear
  //             solve is requested.
  // * "Reuse" - The group/linear solver will not recompute the preconditioner even if the group's solution vector changes. 
  //             It just blindly reuses what has been constructed. This turns off control of preconditioner recalculation. 
  //             This is a dangerous condition but can really speed up the computations if the user knows what they are doing. 
  //             We don't recommend users trying this.
  // * "Recompute" - Recomputes the preconditioner, but will try to efficiently reuse any objects that don't need to be destroyed. 
  //                 How efficient the "Recompute" option is depends on the type of preconditioner. For example if we are using ILU from the
  //                 Ifpack library, we would like to not destroy and reallocate the graph each solve. With this option, we tell
  //                 Ifpack to reuse the graph from last time - e.g the sparisty pattern has not changed between applications of 
  //                 the preconditioner.
  get_int_parameter(param_addr, "newt_ls_preconditioner_reuse_policy", &tmp_int, &tmp_res, 1, CHECK_VALUES, 3, 1, 2, 3);
  if (tmp_res!=-1) 
    {
      switch(tmp_int)
	{
	case 1:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Preconditioner Reuse Policy", "Rebuild");  
	  break;
	case 2:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Preconditioner Reuse Policy", "Reuse");  
	  break;
	case 3:
	  dirParams.sublist("Newton").sublist("Linear Solver").set("Preconditioner Reuse Policy", "Recompute");  
	  break;
	}
    }

  // "Max Age Of Prec" - (int) If the "Preconditioner Reuse Policy" is set to "Reuse", this integer tells the linear system how many 
  //                     times to reuse the preconditioner before rebuilding it. Defaults to 1.
  get_int_parameter(param_addr, "newt_ls_max_age_of_prec", &tmp_int, &tmp_res, 1, CHECK_MIN, 0);
  if (tmp_res!=-1) dirParams.sublist("Newton").sublist("Linear Solver").set("Max Age Of Prec", tmp_int); 

  // "RCM Reordering" - Enables RCM reordering in conjunction with domain decomp incomplete factorization preconditioning. 
  // The following options are valid:
  // * "Disabled" - (default).
  // * "Enabled"
  get_int_parameter(param_addr, "newt_ls_rcm_reordering", &tmp_int, &tmp_res, 0, CHECK_VALUES, 2, 0, 1);
  if (tmp_res!=-1) 
    {
      if (tmp_int) dirParams.sublist("Newton").sublist("Linear Solver").set("RCM Reordering", "Enabled"); 
      else         dirParams.sublist("Newton").sublist("Linear Solver").set("RCM Reordering", "Disabled"); 
    }

#ifdef USE_AZTECOO
  // "Use Adaptive Linear Solve" - Enables the use of AztecOO's AdaptiveIterate() method instead of calling the Iterate() method. 
  // This causes the preconditioning matrix to be modified to make the linear solves easier. AztecOO will attempt to solve the linear system 
  // multiple times now and if the solves are failing it will modify the preconditioner and try again. 
  // Boolean value, defaults to false. 
  // NOTE: This only works for internal Aztec preconditioners! 
  // The "Preconditioning" parameter must be set to "AztecOO: Jacobian Matrix" or "AztecOO: User RowMatrix". 
  // (NOTE: This parameter is currently NOT supported)
  get_int_parameter(param_addr, "newt_ls_use_adaptive_linear_solve", &tmp_int, &tmp_res, 0, CHECK_VALUES, 2, 0, 1);
  if (tmp_res!=-1) 
    {
      if (tmp_int) dirParams.sublist("Newton").sublist("Linear Solver").set("RCM Reordering", "Enabled"); 
      else         dirParams.sublist("Newton").sublist("Linear Solver").set("RCM Reordering", "Disabled"); 
    }

  // "Compute Scaling Manually" - (boolean) The linear system can be scaled if a NOX::Epetra::Scaling object is supplied to LinearSystemAztecOO. 
  // When to compute the scaling can be handled either manually by the user, or this object can automatically compute the scaling prior
  // to a linear solve. By setting this flag to true, the user will call NOX::Epetra::Scaling::computeScaling() manually - on their own!
  // Setting this to false means the LinearSystemAztecOO object will call the computeScaling function right before it applies the scaling
  //  to the matrix in the applyJacobianInverse function. Default is true (user will call compute scaling).
  get_int_parameter(param_addr, "newt_ls_compute_scaling_manually", &tmp_int, &tmp_res, 0, CHECK_VALUES, 2, 0, 1);
  if (tmp_res!=-1) 
    {
      if (tmp_int) dirParams.sublist("Newton").sublist("Linear Solver").set("Compute Scaling Manually", true); 
      else         dirParams.sublist("Newton").sublist("Linear Solver").set("Compute Scaling Manually", false); 
    }
#endif

  dirParams.sublist("Newton").sublist("Linear Solver").set("Compute Scaling Manually", false); 

  // "Output Solver Details" - (boolean) Write the output sublist below to the parameter list after each linear solve. default is true.
  get_int_parameter(param_addr, "newt_ls_output_solver_details", &tmp_int, &tmp_res, 1, CHECK_VALUES, 2, 0, 1);
  if (tmp_res!=-1) 
    {
      if (tmp_int) dirParams.sublist("Newton").sublist("Linear Solver").set("Output Solver Details", true); 
      else         dirParams.sublist("Newton").sublist("Linear Solver").set("Output Solver Details", false); 
    }

  // Create the interface between the test problem and the nonlinear solver
  // This is created by the user using inheritance of the abstract base class:
  // NOX_Epetra_Interface
  Teuchos::RCP<Problem_Interface> interface = Teuchos::rcp(new Problem_Interface());

  Teuchos::RCP<Epetra_Operator>                  Oper;
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac;
  Teuchos::RCP<Epetra_CrsMatrix>                 Analytic;
  Teuchos::RCP<NOX::Epetra::MatrixFree>          MF;
  bool use_analytic_selected    = false;
  bool use_matrix_free_selected = false;

  // Create the Epetra_RowMatrix.  Uncomment one or more of the following:
  // 1. User supplied (Epetra_RowMatrix)

  get_int_parameter(param_addr, "use_analytic", &tmp_int, &tmp_res, 1, CHECK_VALUES, 2, 0, 1);
  if (tmp_res!=-1)
    {
      m_x_tmp = m_x_in;
      n_x_tmp = n_x_in;
      CreateVar(FUNC_TMP_1, MATRIX_OF_DOUBLE_DATATYPE, &m_x_tmp, &n_x_tmp, &l_x_tmp);
      for(i=0;i<size_x;i++) 
	{
	  *stk(l_x_tmp+i) = (*soln)[i];
	}
      
      m_jac_tmp = 1;
      n_jac_tmp = 1;
      CreateVar(FUNC_TMP_2, MATRIX_OF_DOUBLE_DATATYPE, &m_jac_tmp, &n_jac_tmp, &l_jac_tmp);
      
      PExecSciFunction(FUNC_TMP_1, &sci_obj, &lhs_obj, &rhs_obj, "InitialJacobianComputation", call_f_env);
      
      if (GetType(FUNC_TMP_2)==5) // Sparse matrix
	{
	  SciSparse SpResult;
	  int i, j, Index;
	  
	  GetRhsVar(FUNC_TMP_2, SPARSE_MATRIX_DATATYPE, &m_jac_tmp, &n_jac_tmp, &SpResult);
	  
#ifdef DEBUG
	  sciprint("DEBUG: Sparse matrix m = %d n = %d\n", m_jac_tmp, n_jac_tmp);
#endif
	  
	  if (n_jac_tmp!=size_x)
	    {
	      Scierror(999,_("fsolver_nox: error, you must return a %d * %d sparse matrix\n"),m_jac_tmp,size_x);
	    }
	  
	  size_f = m_jac_tmp;
	  Epetra_Map MapRow(size_f, 0, Comm);
	  // YC: Analytic = Teuchos::rcp(new Epetra_CrsMatrix(View,MapRow,MapCol,size_x,false));
	  Analytic = Teuchos::rcp(new Epetra_CrsMatrix(Copy,MapRow,MapCol,size_x,false));
	  
 	  Index = 0;
 	  for(i=0;i<SpResult.m;i++)
 	    {
 	      for(j=0;j<SpResult.mnel[i];j++)
 		{
 		  A = SpResult.R[Index];
 		  irow = i;
 		  jcol = SpResult.icol[Index]-1;
 		  Index++;
#ifdef DEBUG
 		  sciprint("DEBUG: i = %d, j = %d, val = %f\n", irow, jcol, A);
#endif
 		  Analytic->InsertGlobalValues(irow, 1, &A, &jcol);
 		}
 	    }
	}
      else
	{
	  GetRhsVar(FUNC_TMP_2, MATRIX_OF_DOUBLE_DATATYPE, &m_jac_tmp, &n_jac_tmp, &l_jac_tmp);
	  
#ifdef DEBUG
	  sciprint("DEBUG: Dense matrix m = %d n = %d\n", m_jac_tmp, n_jac_tmp);
#endif
	  
	  if (n_jac_tmp!=size_x)
	    {
	      Scierror(999,_("fsolver_nox: error, you must return a %d * %d full matrix\n"),m_jac_tmp,size_x);
	    }
	  
	  size_f = m_jac_tmp;
	  Epetra_Map MapRow(size_f, 0, Comm);
	  //Analytic = Teuchos::rcp(new Epetra_CrsMatrix(View,MapRow,MapCol,size_x,false));
	  Analytic = Teuchos::rcp(new Epetra_CrsMatrix(Copy,MapRow,MapCol,size_x,false));

 	  for(i=0;i<m_jac_tmp;i++)
 	    {
 	      for(j=0;j<n_jac_tmp;j++)
 		{
 		  A = *stk(l_jac_tmp+i+j*m_jac_tmp);
 		  irow = i;
 		  jcol = j;
#ifdef DEBUG
 		  sciprint("DEBUG: i = %d, j = %d, val = %f\n", irow, jcol, A);
#endif
 		  Analytic->InsertGlobalValues(irow, 1, &A, &jcol);
 		}
 	    }
	}
      
      Analytic->FillComplete(false); // If true, a storage optimization is performed. The sparsity structure will not be conserved !
      Analytic->OptimizeStorage();
      Analytic->MakeDataContiguous();

      // Now, we can set the size of the matrix
      Oper = Analytic;
      //iJac = Teuchos::rcp(new NOX::Epetra::MatrixFree(printParams, interface, noxSoln)); // YC: BEFORE: Analytic;
      iJac = Teuchos::rcp(new Problem_Interface());
      use_analytic_selected = true;
    }

  // 2. Matrix-Free (Epetra_Operator)
  get_int_parameter(param_addr, "use_matrix_free", &tmp_int, &tmp_res, 1, CHECK_VALUES, 2, 0, 1);
  if (tmp_res!=-1)
    {
      Teuchos::RCP<NOX::Epetra::MatrixFree> MF = Teuchos::rcp(new NOX::Epetra::MatrixFree(printParams, interface, noxSoln));
      Oper = MF;
      iJac = MF;
      use_matrix_free_selected = true;
    }

  //////////////////////////////////////////
  // Finite Difference (Epetra_RowMatrix) //
  //////////////////////////////////////////
  Teuchos::RCP<NOX::Epetra::FiniteDifference> FD;
  get_int_parameter(param_addr, "use_finite_difference", &tmp_int, &tmp_res, 1, CHECK_VALUES, 2, 0, 1);
  // By default, is matrix_free and analytic have not been selected, then we will use finite differences
  if ((tmp_res!=-1) || (!use_analytic_selected && !use_matrix_free_selected))
    {
      get_double_parameter(param_addr, "use_finite_difference_alpha", &tmp_double, &tmp_res, 1.0e-4, CHECK_MIN, 0.0);
      get_double_parameter(param_addr, "use_finite_difference_beta", &tmp_double_2, &tmp_res, 1.0e-6, CHECK_MIN, 0.0);
      FD = Teuchos::rcp(new NOX::Epetra::FiniteDifference(printParams, interface, noxSoln, tmp_double_2, tmp_double));
      get_int_parameter(param_addr, "use_finite_difference_type", &tmp_int_2, &tmp_res_2, 1, CHECK_VALUES, 3, 1, 2, 3);
      if (tmp_res_2!=-1)
	{
	  switch(tmp_int_2)
	    {
	    case 1:
	      FD->setDifferenceMethod(NOX::Epetra::FiniteDifference::Forward);
	      break;
	    case 2:
	      FD->setDifferenceMethod(NOX::Epetra::FiniteDifference::Backward);
	      break;
	    case 3:
	      FD->setDifferenceMethod(NOX::Epetra::FiniteDifference::Centered);
	      break;
	    }
	}
      Oper = FD;
      iJac = FD;
    }

  // Create the linear system
  Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = interface;

#ifdef USE_AZTECOO
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys = 
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(nlParams.sublist("Printing"), 
						      dirParams.sublist("Newton").sublist("Linear Solver"),
						      iReq,
						      iJac,
						      Oper, // Analytic, MF, FD
						      noxSoln));
#else
  Teuchos::RCP<NOX::Epetra::LinearSystemAmesos> linSys = 
    Teuchos::rcp(new NOX::Epetra::LinearSystemAmesos(nlParams.sublist("Printing"), 
						     dirParams.sublist("Newton").sublist("Linear Solver"),
						     iReq,
						     iJac,
						     Oper, // Analytic, MF, FD
						     noxSoln));
#endif
  
  // Create the Group
  Teuchos::RCP<NOX::Epetra::Group> grp = Teuchos::rcp(new NOX::Epetra::Group(printParams, iReq, noxSoln, linSys)); 

  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::Combo> combo     = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  Teuchos::RCP<NOX::StatusTest::Combo> converged = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
  combo->addStatusTest(converged);


  // Various convergence tests based on the norm of the residual.
  // the Absolute tolerance type. 
  get_double_parameter(param_addr, "test_absresid", &tmp_double, &tmp_res, 1.0e-8, CHECK_MIN, 0.0);
  if (tmp_res!=-1)
    {
      Teuchos::RCP<NOX::StatusTest::NormF> absresid = Teuchos::rcp(new NOX::StatusTest::NormF(tmp_double));
      converged->addStatusTest(absresid);
    }

  // the Relative tolerance type. 
  get_double_parameter(param_addr, "test_relresid", &tmp_double, &tmp_res, 1.0e-2, CHECK_MIN, 0.0);
  if (tmp_res!=-1)
    {
      Teuchos::RCP<NOX::StatusTest::NormF> relresid = Teuchos::rcp(new NOX::StatusTest::NormF(*grp.get(), tmp_double));
      converged->addStatusTest(relresid);
    }

  // Various convergence tests based on the norm of the change in the solution vector, x, between outer iteration
  // the Absolute ToleranceType and TWO NormType. 
  get_double_parameter(param_addr, "test_update", &tmp_double, &tmp_res, 1.0e-5, CHECK_MIN, 0.0);
  if (tmp_res!=-1)
    {
      Teuchos::RCP<NOX::StatusTest::NormUpdate> update = Teuchos::rcp(new NOX::StatusTest::NormUpdate(tmp_double));
      converged->addStatusTest(update);
    }

  // Convergence test based on the weighted root mean square norm fo the solution update between iterations.
  get_double_parameter(param_addr, "test_wrms_rtol", &tmp_double, &tmp_res, 1.0e-5, CHECK_MIN, 0.0);
  get_double_parameter(param_addr, "test_wrms_atol", &tmp_double_2, &tmp_res_2, 1.0e-5, CHECK_MIN, 0.0);
  if ((tmp_res!=-1) || (tmp_res_2!=-1))
    {
      Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms = Teuchos::rcp(new NOX::StatusTest::NormWRMS(tmp_double, tmp_double_2));
      converged->addStatusTest(wrms);
    }

  // Failure test based on the maximum number of nonlinear solver iterations. 
  get_int_parameter(param_addr, "test_maxiters", &tmp_int, &tmp_res, 20, CHECK_MIN, 0);
  if (tmp_res!=-1)
    {
      Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = Teuchos::rcp(new NOX::StatusTest::MaxIters(tmp_int));
      combo->addStatusTest(maxiters);
    }

  // Failure test based on whether the norm of a vector has a finite value. 
  get_int_parameter(param_addr, "test_finite_value", &tmp_int, &tmp_res, 0, CHECK_VALUES, 2, 0, 1);
  if (tmp_res!=-1)
    {
      Teuchos::RCP<NOX::StatusTest::FiniteValue> fv = Teuchos::rcp(new NOX::StatusTest::FiniteValue);
      combo->addStatusTest(fv);
    }

  // Failure test based on a threshold value of the norm of F. 
  get_double_parameter(param_addr, "test_divergence", &tmp_double, &tmp_res, 1.0e6, CHECK_MIN, 0.0);
  if (tmp_res!=-1)
    {
      Teuchos::RCP<NOX::StatusTest::Divergence> dv = Teuchos::rcp(new NOX::StatusTest::Divergence(tmp_double));
      combo->addStatusTest(dv);
    }

  // Failure test based on the convergence rate between nonlinear iterations. 
  get_double_parameter(param_addr, "test_stagnation_threshold", &tmp_double, &tmp_res, 1.0, CHECK_MIN, 0.0);
  get_int_parameter(param_addr, "test_stagnation_iterations", &tmp_int, &tmp_res_2, 50,CHECK_MIN, 0);
  if ((tmp_res!=-1)||(tmp_res_2!=-1))
    {
      Teuchos::RCP<NOX::StatusTest::Stagnation> stag = Teuchos::rcp(new NOX::StatusTest::Stagnation(tmp_int,tmp_double));
      combo->addStatusTest(stag);
    }

  // Create the method
  Teuchos::RCP<Teuchos::ParameterList> finalParamsPtr = nlParamsPtr;

  Teuchos::RCP<NOX::Solver::Generic> solver = NOX::Solver::buildSolver(grp, combo, finalParamsPtr);

  // Status of the solver:
  // -2 Unevaluated Unevaluated.
  //  0 Unconverged Neither Converged nor Failed.
  //  1 Converged   Converged.
  // -1 Failed      Failed. 

  NOX::StatusTest::StatusType status;

  try
    {
      status = solver->solve();
    }
  catch(const char * msg)
    {
      sciprint("fsolve_nox: Error - %s\n",msg);
      status = NOX::StatusTest::Failed;
    }

  if (status == NOX::StatusTest::Converged)
    utils.out() << "Test Passed!" << endl;
  else 
    {
      if (MyPID==0) utils.out() << "Nonlinear solver failed to converge!" << endl;
    }
  
  // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group& finalGroup    = dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
  const Epetra_Vector&      finalSolution = (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

  m_x_out = size_x; n_x_out = 1;
  _SciErr = allocMatrixOfDouble(pvApiCtx, X_OUT, m_x_out, n_x_out, &x_out);
  for(i=0;i<size_x;i++) x_out[i] = finalSolution[i];

  // End Nonlinear Solver **************************************

#ifdef DEBUG
  // Output the parameter list
  if (utils.isPrintType(NOX::Utils::Parameters)) 
    {
      utils.out() << endl << "Final Parameters" << endl
		  << "****************" << endl;
      solver->getList().print(utils.out());
      utils.out() << endl;
    }
#endif

  _SciErr = createMList(pvApiCtx, EXTRA_OUT, 3, &extra_addr);
  _SciErr = createMatrixOfStringInList(pvApiCtx, EXTRA_OUT, extra_addr, 1, m_list_labels, n_list_labels, (char **)ListLabels);

  // the status of the nox solver
  tmp_dbl[0] = status;
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, EXTRA_OUT, extra_addr, 2, 1, 1, tmp_dbl);
  
  // the iterations spend by the nox solver
  tmp_dbl[0] = solver->getNumIterations();
  _SciErr = createMatrixOfDoubleInList(pvApiCtx, EXTRA_OUT, extra_addr, 3, 1, 1, tmp_dbl);

  LhsVar(1) = X_OUT;
  LhsVar(2) = EXTRA_OUT;

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif
  
  return 0;
}
