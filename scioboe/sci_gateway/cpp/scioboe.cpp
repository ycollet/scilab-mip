/////////////////////////////////////////////////
// scioboe: A scilab interface to OBOE library //
/////////////////////////////////////////////////

//  Copyright (C) 2009-2010 Yann Collette.
//
//  scioboe is free software; you can redistribute it and/or modify it
//  under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 2, or (at your option) any
//  later version.
//
//  This part of code is distributed with the FURTHER condition that it 
//  can be compiled and linked with the Matlab libraries and it can be 
//  used within the Matlab environment.
//
//  scioboe is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with SciCoinOR; see the file COPYING.  If not, write to the Free
//  Software Foundation, 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#include <iomanip>
#include <string>
#include <new>
#include <exception>

#include <lafnames.h>
#include <symd.h>
#include <bmd.h>
#include <lavd.h>
#include <gmd.h>
#include <ltgmd.h>

#include <AccpmBlasInterface.hpp>
#include <Manager.hpp>
#include <Parameters.hpp>
#include <Oracle.hpp>
#include <QpGenerator.hpp>

using namespace Accpm;
using namespace std;

extern "C"
{
#include <stack-c.h>
#include <sciprint.h>
#include <Scierror.h>
#include <MALLOC.h>
#include <api_parameters.h>
}

#include <api_scilab.h>

#include <helper.hpp>

//#define DEBUG 1
//#define DBGPRINTF printf
#define DBGPRINTF sciprint

#include <string.h>

#define X_IN                   1
#define FOBJ_IN                2
#define X_LOWER_IN             3
#define X_UPPER_IN             4
#define D_EQ_IN                5
#define RHS_EQ_IN              6
#define CENTER_BALL_IN         7
#define PARAM_IN               8
#define LAST_PARAMS            PARAM_IN
#define X_OUT                  9
#define STATUS_OUT             10

// A structure which handles informations related to objective functions
struct sci_oboe_info
{
  int fobj_lhs,  fobj_rhs,  l_fobj;
  int ibegin;            // the position of the top of the stack
  int nb_var;
  int nb_constr;
};

class SciOracleFunction : public OracleFunction {

private:
  struct sci_oboe_info * sci_parameters;
public:
  SciOracleFunction() : OracleFunction() {}

  void set_scilab_parameters(struct sci_oboe_info * param) {sci_parameters = param;}
  struct sci_oboe_info * get_scilab_parameters() {return sci_parameters;}

  virtual int eval(const AccpmVector &y, AccpmVector &functionValue, AccpmGenMatrix &subGradients, AccpmGenMatrix *info) 
  {
    int n_function_value, m_function_value, * function_value_addr = NULL;
    int n_subgradient,    m_subgradient,    * subgradient_addr = NULL;
    int n_info,           m_info,           * info_addr = NULL;
    int n_y,              m_y;
    int n_tmp,            m_tmp;
    int rhs_old = Rhs, nbvars_old = Nbvars;
    int i, j, Log = 0, nb_cuts = 0;
    double * int_y = NULL, * int_function_value = NULL, * int_subgradient = NULL, * int_info = NULL;
    SciErr _SciErr;

    //////////////////////////
    // Call to the function //
    //////////////////////////
    
    Rhs    = sci_parameters->ibegin + max(sci_parameters->fobj_rhs,sci_parameters->fobj_lhs);
    Nbvars = sci_parameters->ibegin + max(sci_parameters->fobj_rhs,sci_parameters->fobj_lhs);
    
    // Store the current point into y

    n_y = y.size();
    m_y = 1;

    _SciErr = allocMatrixOfDouble(pvApiCtx, sci_parameters->ibegin+0, n_y, m_y, &int_y); SCICOINOR_ERROR;

    for(i=0;i<n_y;i++) int_y[i] = y(i);

#ifdef DEBUG
    sciprint("eval: n_y = %d m_y = %d\n", n_y, m_y);
    sciprint("eval: Top = %d Rhs+1 = %d ibegin = %d\n", Top, Rhs+1, sci_parameters->ibegin);
#endif

    n_tmp = 1; m_tmp = 1;
    _SciErr = allocMatrixOfDouble(pvApiCtx, sci_parameters->ibegin+1, n_tmp, m_tmp, &int_y); SCICOINOR_ERROR;
    _SciErr = allocMatrixOfDouble(pvApiCtx, sci_parameters->ibegin+2, n_tmp, m_tmp, &int_y); SCICOINOR_ERROR;

    // Call to the scilab function 
    try
    {
      SciFunction(&sci_parameters->ibegin,&sci_parameters->l_fobj,&sci_parameters->fobj_lhs,&sci_parameters->fobj_rhs);
    }
    catch(...)
      {
	Scierror(999,"scioboe: error when calling objective function\n");
	Rhs    = rhs_old;
	Nbvars = nbvars_old;
	return 0;
      }
    
    if (Err>0) 
      {
	Scierror(999,"scioboe: error when calling objective function\n");
	Rhs    = rhs_old;
	Nbvars = nbvars_old;
	return 0;
      } /* End If */
    
    // Get function_value

    // functionValue which is an vector of size NumCuts x 1
    // This value represents the following:
    // * For the non-smooth function f1(.) :
    //  * For optimality cuts, i.e *info = 1, functionValue is typically the function evaluation at the current query point, y.
    //  * For feasiblity cuts, i.e *info = 0, however, this is a value which would make the cut described by the subGradients matrix valid in the 
    //    form :
    //    functionValue + subGradients^T*(y' - y) <= 0, for all feasible points y'
    // * For the smooth function f2(.) :
    //   It again represents the function value f2(y). 

    _SciErr = getVarAddressFromPosition(pvApiCtx, sci_parameters->ibegin+0, &function_value_addr); SCICOINOR_ERROR;
    _SciErr = getMatrixOfDouble(pvApiCtx, function_value_addr, &n_function_value, &m_function_value, &int_function_value); SCICOINOR_ERROR;

    nb_cuts = n_function_value;

    sciprint("eval: n_function_value = %d, m_function_value = %d\n", n_function_value, m_function_value);

    functionValue.resize(n_function_value,m_function_value);
    for(i=0;i<n_function_value*m_function_value;i++) functionValue(i) = int_function_value[i];
    
    // Get subgradient

    // Sub-gradient at y in matrix subGradients which of size NumVariables x NumCuts. 
    // The user can provide more than one cut for each query point, which usually is required when the NumSubProblems is more than 1.
    // The subgradient vector(with abuse of name) is also used to provide a valid cuts incase the given point, y, is not feasible.
    // Hence the eval function is responsible for providing both optimality and feasibility cuts.

    _SciErr = getVarAddressFromPosition(pvApiCtx, sci_parameters->ibegin+1, &subgradient_addr); SCICOINOR_ERROR;
    _SciErr = getMatrixOfDouble(pvApiCtx, subgradient_addr, &n_subgradient, &m_subgradient, &int_subgradient); SCICOINOR_ERROR;

    sciprint("eval: n_subgradient = %d, m_subgradient = %d\n", n_subgradient, m_subgradient);

    if ((n_subgradient != sci_parameters->nb_var) && (m_subgradient != nb_cuts))
      {
	Scierror(999,"%s: subgrad_val must be of size %d x %d\n", "oboe", sci_parameters->nb_var, nb_cuts);
	Rhs    = rhs_old;
	Nbvars = nbvars_old;
	return 0;
      }

    subGradients.resize(n_subgradient, m_subgradient);
    for(i=0;i<n_subgradient;i++)
      {
	for(j=0;j<m_subgradient;j++)
	  {
	    subGradients(i,j) = *(int_subgradient + i + j*n_subgradient);
	  }
      }

    // Get info

    // AccpmGenMatrix *info:
    // * For non-smooth function f1 it is a vector of dimension NumCuts x 1.
    //   It specifies the cut type for each cut:
    //   * (i) For feasibility cut,   *info(i,0) = 0 
    //   * (ii) For optimality cut i, *info(i,0) = i, the index of subproblem.
    // * For the smooth function f2, this vector has the Hessian information. 
    //   It is a vector of dimension NumVariables x NumVariables (if parameter diagHessian is false) otherwise its has dimension NumVariables x 1

    if (info)
      {
	_SciErr = getVarAddressFromPosition(pvApiCtx, sci_parameters->ibegin+2, &info_addr); SCICOINOR_ERROR;
	_SciErr = getMatrixOfDouble(pvApiCtx, info_addr, &n_info, &m_info, &int_info); SCICOINOR_ERROR;
	
	if ((n_info != nb_cuts)&&(m_info != 1))
	  {
	    Scierror(999,"%s: info must be of size %d x 1\n", "oboe", nb_cuts);
	    Rhs    = rhs_old;
	    Nbvars = nbvars_old;
	    return 0;
	  }
	
	sciprint("eval: n_info = %d, m_info = %d\n", n_info, m_info);

	// YC: check the dimension of info
	
	info->resize(n_info,m_info);
	for(i=0;i<n_info;i++)
	  {
	    for(j=0;j<m_info;j++)
	      {
	    	(*info)(i,j) = *(int_info + i + j*n_info);
	      }
	  }
      }

    // The return value is 0 on success and 1 to terminate the Outer Iterations and hence the Query point generation process.

    Rhs    = rhs_old;
    Nbvars = nbvars_old;

    return 0;
  }
};

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
      sciprint("scioboe: %s\n",m_string.c_str());

    m_stream.rdbuf(m_old_buf);
  }

protected:
  virtual int_type overflow(int_type v)
  {
    if (v == '\n')
      {
	sciprint("scioboe: %s\n",m_string.c_str());
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
	    sciprint("scioboe: %s\n",tmp.c_str());
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

/////////////////////////
// QpGenerator changes //
/////////////////////////

class SciQpGenerator : public QpGenerator
{
public:
  virtual void output(std::ostream &os)
  {
    ScilabStream scicout(os);

    _timer.stop();
    if (_manager) 
      {
    	os << "\n-------------------------------------------------------------------------" << endl;
    	os << "OBOE using ACCPM - Vital Statistics" << endl;
    	_manager->output(os);
    	os << "Number of Outer Iterations:\t" << _manager->getNumOuterIteration() << endl;
    	os << "Inner / Outer Iterations:\t" << 
    	  _manager->getNumInnerIterations()*1.0/_manager->getNumOuterIteration() << endl;
    	double totalTime = _timer.getCpuTime();
    	double oraclePercent = 0;
    	if (totalTime > 0)
    	  {
    	    oraclePercent = (_oracleTime * 100) / totalTime;
    	  }
    	os << "Total Elapsed Time:\t\t" << _timer.getRealTime() << " sec" << endl;
    	os << "Total Cpu Time:\t\t\t" << totalTime << " sec" << endl;
    	os << "Oracle Time:\t\t\t" << _oracleTime << " sec (" 
    	   << oraclePercent << " %)" 
    	   << endl;
	
    	if (_param->getIntParameter("Verbosity") > 5) 
    	  {
    	    os << "\nSolution y:\n" << _manager->getBestY() << std::endl;
    	    os << "Input Parameters:" << *_param << std::endl;
    	  }
    	os << "\n-------------------------------------------------------------------------" << endl;
      }
  }
  virtual void printHeader(std::ostream &s = std::cout) const
  {
    ScilabStream scicout(s);

    if (_manager->getNumOuterIteration() == 0) 
      {
    	if (_param->getIntParameter("Verbosity") == 1)
    	  {
    	    s << "Iteration\tObjective" << std::endl;
    	  } 
    	else if (_param->getIntParameter("Verbosity") == 2)
    	  {
    	    s << "Outer\tInner\tFeasible\t\tObjective\t\tUBound\t\tLBound\tRelativeGap" << std::endl;	
    	  }
      }
  }
#ifdef WIN32
  virtual void printIteration(const AccpmVector &_val, const AccpmVector &_val2, const AccpmGenMatrix &_subGrad,
  			      const AccpmGenMatrix &_subProblemIndex, std::ostream &os = std::cout) const
#else
  virtual void printIteration(const AccpmVector &val, const AccpmVector &val2, const AccpmGenMatrix &subGrad,
  			      const AccpmGenMatrix &subProblemIndex, std::ostream &os = std::cout) const
#endif
  {
#ifdef WIN32
    AccpmVector val(_val);
    AccpmVector val2(_val2);
    AccpmGenMatrix subGrad(_subGrad);
    AccpmGenMatrix subProblemIndex(_subProblemIndex);
#endif
    ScilabStream scicout(std::cout);

    const AccpmVector &y = _manager->getCurrentY();
    double fVal = val2(0); 
    if (_manager->isCurrentPointFeasible()) 
      { 
  	fVal += AccpmLADotProd(val, *_param->getPi()); 
  	if (_param->getB()) 
  	  {
  	    fVal +=  AccpmLADotProd(y, *_param->getB());
  	  }
      } 
    else 
      {
  	fVal += val.sum();  /* For now we display the sum of cut values for feasibility cut */ 
      }
    if (_param->getIntParameter("Verbosity") > 2) 
      {
  	std::cout << "Iteration " << _manager->getNumOuterIteration() << ":" << std::endl;
  	std::cout << "y:\n " << y << std::endl;
  	std::cout << "Objective:" << fVal << "\n" << std::endl;
  	if (_param->getIntParameter("Verbosity") > 3) 
  	  {
  	    std::cout << "Subgradient:\n" << subGrad << std::endl;
  	    std::cout << "SubProblemIndex:\n" << subProblemIndex << std::endl;
  	  }
      } 
    else 
      {
  	if (_param->getIntParameter("Verbosity") == 1) 
  	  {
  	    std::cout << std::setw(4) << _manager->getNumOuterIteration() 
  		      << std::setw(18) << fVal << "\n" << std::endl;
  	  } 
  	else if (_param->getIntParameter("Verbosity") == 2) 
  	  {
  	    std::cout << std::setw(4) << _manager->getNumOuterIteration() 
  		      << std::setw(8) << _manager->getNumInnerIterations() 
  		      << std::setw(8) << _manager->isCurrentPointFeasible()
  		      << " " << std::setw(16) << fVal 
  		      << " " << std::setw(16) << getObjUB() 
  		      << " " << std::setw(16) << getObjLB() 
  		      << "\t" << _manager->getRelativeGap() << std::endl;  
  	  }
      }
   }
};
  
extern "C" int scioboe(char * fname)
{
  int m_x,              n_x,              * x_addr = NULL;
  int m_x_out,          n_x_out;
  int m_d_eq_in,        n_d_eq_in,        * d_eq_in_addr = NULL;
  int m_rhs_eq_in,      n_rhs_eq_in,      * rhs_eq_in_addr = NULL;
  int m_center_ball_in, n_center_ball_in, * center_ball_in_addr = NULL;
  int m_fobj,           n_fobj,           l_fobj;
  int m_x_lower,        n_x_lower,        * x_lower_addr = NULL;
  int m_x_upper,        n_x_upper,        * x_upper_addr = NULL;
  int nb_var, i, j, nb_constr;
  int tmp_res, tmp_int, Log = 0;
  int  * param_in_addr = NULL;
  char *  tmp_char = NULL;
  double  tmp_double;
  struct sci_oboe_info * sci_parameters = new struct sci_oboe_info;
  vector<double> start, varLB, varUB;
  StdRealVector CenterBall;
  //Parameters param("param.oboe");
  Parameters oboe_param;
  ScilabStream scicout(std::cout);
  AccpmGenMatrix D_Eq;
  AccpmVector D_Rhs;
  SciOracleFunction f1;
  Oracle oracle(&f1);
  SciQpGenerator qpGen;
  double * x = NULL, * x_lower = NULL, * x_upper = NULL, * d_eq_in = NULL;
  double * rhs_eq_in = NULL, * center_ball_in = NULL;
  SciErr _SciErr;

  if (Rhs<LAST_PARAMS) 
    {
      Scierror(999,"%s: %d inputs required in call to %s.\n", fname, LAST_PARAMS, fname);
      return 0;
    }

  _SciErr = getVarAddressFromPosition(pvApiCtx, X_IN, &x_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, x_addr, &n_x, &m_x, &x); SCICOINOR_ERROR;

  GetRhsVar(FOBJ_IN, EXTERNAL_DATATYPE, &m_fobj, &n_fobj, &l_fobj);

  _SciErr = getVarAddressFromPosition(pvApiCtx, X_LOWER_IN, &x_lower_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, x_lower_addr, &n_x_lower, &m_x_lower, &x_lower); SCICOINOR_ERROR;
  _SciErr = getVarAddressFromPosition(pvApiCtx, X_UPPER_IN, &x_upper_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, x_upper_addr, &n_x_upper, &m_x_upper, &x_upper); SCICOINOR_ERROR;
  _SciErr = getVarAddressFromPosition(pvApiCtx, D_EQ_IN, &d_eq_in_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, d_eq_in_addr, &n_d_eq_in, &m_d_eq_in, &d_eq_in); SCICOINOR_ERROR;
  // YC: check that size(D_EQ_IN,1) = length(X_IN) and that size(D_EQ_IN,2) = length(RHS_EQ_IN)
  _SciErr = getVarAddressFromPosition(pvApiCtx, RHS_EQ_IN, &rhs_eq_in_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, rhs_eq_in_addr, &n_rhs_eq_in, &m_rhs_eq_in, &rhs_eq_in); SCICOINOR_ERROR;
  _SciErr = getVarAddressFromPosition(pvApiCtx, CENTER_BALL_IN, &center_ball_in_addr); SCICOINOR_ERROR;
  _SciErr = getMatrixOfDouble(pvApiCtx, center_ball_in_addr, &n_center_ball_in, &m_center_ball_in, &center_ball_in); SCICOINOR_ERROR;

  if ((m_d_eq_in!=0)&&(n_d_eq_in!=0))
    {
      if ((m_d_eq_in!=m_x*n_x)||(n_d_eq_in!=m_rhs_eq_in*n_rhs_eq_in))
  	{
  	  Scierror(999,"%s: size of D_EQ must be %d x %d\n", fname, m_rhs_eq_in*n_rhs_eq_in, m_x*n_x);
  	  return 0;
  	}
    }
  nb_var    = m_x*n_x;
  nb_constr = n_d_eq_in;

  start.resize(nb_var);
  varLB.resize(nb_var);
  varUB.resize(nb_var);

  for(i=0;i<nb_var;i++) 
    {
      start[i] = x[i];
      varLB[i] = x_lower[i];
      varUB[i] = x_upper[i];
    }

  sci_parameters->ibegin    = Rhs+1; //Top;
  sci_parameters->fobj_lhs  = m_fobj;
  sci_parameters->fobj_rhs  = n_fobj;
  sci_parameters->l_fobj    = l_fobj;
  sci_parameters->nb_var    = nb_var;
  sci_parameters->nb_constr = nb_constr;

  f1.set_scilab_parameters(sci_parameters);

  // Get the parameters stored in the plist
  initPList(pvApiCtx, PARAM_IN, &param_in_addr);
  if (!checkPList(pvApiCtx, param_in_addr))
    {
      Scierror(999, "%s: argument n° %d is not a plist\n", fname, PARAM_IN);

      return 0;
    }

  // The integer valued functions
  // This parameters is automatically set via the variable X_IN.
  oboe_param.setIntParameter("NumVariables",nb_var);

  getIntInPList(pvApiCtx, param_in_addr, "NumSubProblems", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setIntParameter("NumSubProblems",tmp_int);
  //oboe_param.setIntParameter("NumSubProblems",nb_constr);

  getIntInPList(pvApiCtx, param_in_addr, "MaxOuterIterations", &tmp_int, &tmp_res, 1000, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setIntParameter("MaxOuterIterations",tmp_int);
  getIntInPList(pvApiCtx, param_in_addr, "MaxInnerIterations", &tmp_int, &tmp_res, 50, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setIntParameter("MaxInnerIterations",tmp_int);
  getIntInPList(pvApiCtx, param_in_addr, "Verbosity", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setIntParameter("Verbosity",tmp_int);

  // The boolean functions
  getIntInPList(pvApiCtx, param_in_addr, "Filter", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setIntParameter("Filter",tmp_int);
  getIntInPList(pvApiCtx, param_in_addr, "ConvexityCheck", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setIntParameter("ConvexityCheck",tmp_int);
  getIntInPList(pvApiCtx, param_in_addr, "ConvexityFix", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setIntParameter("ConvexityFix",tmp_int);
  getIntInPList(pvApiCtx, param_in_addr, "FixedProximalCenter", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setIntParameter("FixedProximalCenter",tmp_int);
  getIntInPList(pvApiCtx, param_in_addr, "Proximal", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setIntParameter("Proximal",tmp_int);
  getIntInPList(pvApiCtx, param_in_addr, "DynamicRho", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setIntParameter("DynamicRho",tmp_int);
  getIntInPList(pvApiCtx, param_in_addr, "Ball", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setIntParameter("Ball",tmp_int);
  getIntInPList(pvApiCtx, param_in_addr, "DiagHessian", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setIntParameter("DiagHessian",tmp_int);
  getIntInPList(pvApiCtx, param_in_addr, "ComputeLowerBound", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setIntParameter("ComputeLowerBound",tmp_int);
  getIntInPList(pvApiCtx, param_in_addr, "CheckLocSetInterior", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setIntParameter("CheckLocSetInterior",tmp_int);
 
  // The real valued functions
  getDoubleInPList(pvApiCtx, param_in_addr, "ObjectiveLB", &tmp_double, &tmp_res, ACCPM_MINUS_INF, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setRealParameter("ObjectiveLB", tmp_double);
  getDoubleInPList(pvApiCtx, param_in_addr, "ObjectiveUB", &tmp_double, &tmp_res, ACCPM_PLUS_INF, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setRealParameter("ObjectiveUB", tmp_double);
  getDoubleInPList(pvApiCtx, param_in_addr, "Delta", &tmp_double, &tmp_res, 5, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setRealParameter("Delta", tmp_double);
  getDoubleInPList(pvApiCtx, param_in_addr, "Eta", &tmp_double, &tmp_res, 0.99, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setRealParameter("Eta", tmp_double);
  getDoubleInPList(pvApiCtx, param_in_addr, "Gamma", &tmp_double, &tmp_res, 0.99, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setRealParameter("Gamma", tmp_double);
  getDoubleInPList(pvApiCtx, param_in_addr, "Tolerance", &tmp_double, &tmp_res, 1e-6, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setRealParameter("Tolerance", tmp_double);
  getDoubleInPList(pvApiCtx, param_in_addr, "Rho", &tmp_double, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setRealParameter("Rho", tmp_double);
  getDoubleInPList(pvApiCtx, param_in_addr, "RhoMax", &tmp_double, &tmp_res, 100, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setRealParameter("RhoMax", tmp_double);
  getDoubleInPList(pvApiCtx, param_in_addr, "RhoMin", &tmp_double, &tmp_res, 1e-6, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setRealParameter("RhoMin", tmp_double);
  getDoubleInPList(pvApiCtx, param_in_addr, "WeightEpigraphCutInit", &tmp_double, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setRealParameter("WeightEpigraphCutInit", tmp_double);
  getDoubleInPList(pvApiCtx, param_in_addr, "RadiusBall", &tmp_double, &tmp_res, 1e5, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setRealParameter("RadiusBall", tmp_double);
  getDoubleInPList(pvApiCtx, param_in_addr, "EpsilonReal", &tmp_double, &tmp_res, 1e-10, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setRealParameter("EpsilonReal", tmp_double);
  getDoubleInPList(pvApiCtx, param_in_addr, "EpsilonTol", &tmp_double, &tmp_res, 1e-10, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setRealParameter("EpsilonTol", tmp_double);
  getDoubleInPList(pvApiCtx, param_in_addr, "WeightEpigraphCutInt", &tmp_double, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) oboe_param.setRealParameter("WeightEpigraphCutInc", tmp_double);

  // The string values functions
  getStringInPList(pvApiCtx, param_in_addr, "ProblemName", &tmp_char, &tmp_res, "OBOE General Problem", Log, CHECK_NONE);
  if (tmp_res!=1) oboe_param.setStringParameter("ProblemName", string(tmp_char));
  FREE(tmp_char);
  getStringInPList(pvApiCtx, param_in_addr, "ProblemAcronym", &tmp_char, &tmp_res, "OBOE_Problem", Log, CHECK_NONE);
  if (tmp_res!=1) oboe_param.setStringParameter("ProblemAcronym", string(tmp_char));
  FREE(tmp_char);
  getStringInPList(pvApiCtx, param_in_addr, "MethodName", &tmp_char, &tmp_res, "ProximalACDual", Log, CHECK_NONE);
  if (tmp_res!=1) oboe_param.setStringParameter("MethodName", string(tmp_char));
  FREE(tmp_char);
  getStringInPList(pvApiCtx, param_in_addr, "OptimisationType", &tmp_char, &tmp_res, "Min", Log, CHECK_NONE);
  if (tmp_res!=1) oboe_param.setStringParameter("OptimizationType", string(tmp_char));
  FREE(tmp_char);
  getStringInPList(pvApiCtx, param_in_addr, "LPSolverName", &tmp_char, &tmp_res, "GLPK", Log, CHECK_NONE);
  if (tmp_res!=1) oboe_param.setStringParameter("LPSolverName", string(tmp_char));
  FREE(tmp_char);
  
  oboe_param.setStartingPoint(start);
  oboe_param.setVariableLB(varLB);
  oboe_param.setVariableUB(varUB);

  if (m_d_eq_in*n_d_eq_in!=0)
    {
      D_Eq.resize(nb_var,nb_constr);
      for(i=0;i<nb_constr;i++)
  	{
  	  for(j=0;j<nb_var;j++)
  	    {
  	      D_Eq(j,i) = *(d_eq_in + i + j*nb_constr);
  	      //*(D_Eq.addr() + i + j*m_d_eq_in) = *stk(l_d_eq_in + i + j*m_d_eq_in);
  	    }
  	}

      if (m_rhs_eq_in*n_rhs_eq_in!=0)
  	{
  	  D_Rhs.resize(nb_constr);
  	  for(i=0;i<nb_constr;i++) 
  	    {
  	      D_Rhs(i) = rhs_eq_in[i];
  	      //*(D_Rhs.addr()+i) = *stk(l_rhs_eq_in + i);
  	    }
  	}
      else
  	{
  	  D_Rhs.resize(nb_constr);
  	  D_Rhs = 0;
  	}
      oboe_param.addEqualityConstraints(D_Eq, D_Rhs);
    }

  if (m_center_ball_in*n_center_ball_in!=0)
    {
      CenterBall.resize(nb_var);
      for(i=0;i<nb_var;i++) 
  	{
  	  CenterBall[i] = center_ball_in[i];
  	}

      oboe_param.setCenterBall(CenterBall);
    }

  qpGen.init(&oboe_param, &oracle);

  while (!qpGen.run()) 
    {
    }
					   
  qpGen.output(cout);
  qpGen.terminate();
  
  // LOCSET_EMPTY          -5
  // CONVEXITY_FAILURE     -4
  // LA_ERROR              -3
  // CHOLESKY_FAILURE      -2
  // UNKNOWN               -1
  // ITERATING              0
  // RELATIVE_GAP_REACHED   2
  // USER_STOP              3
  // MAX_OUTER_ITERATIONS   4

  double * x_out = NULL;

  m_x_out = qpGen.getCurrentX()->size(); n_x_out = 1;
  x_out = (double *)MALLOC(m_x_out*sizeof(double));

  for(i=0;i<m_x_out;i++) x_out[i] = ((LaVectorDouble *)qpGen.getCurrentX())->operator()(i);

  _SciErr = createMatrixOfDouble(pvApiCtx, X_OUT, m_x_out, n_x_out, x_out); SCICOINOR_ERROR;
  FREE(x_out);

  createScalarDouble(pvApiCtx, STATUS_OUT, (double)qpGen.getExitCode());

  LhsVar(1) = X_OUT;
  LhsVar(2) = STATUS_OUT;

  if (sci_parameters) delete sci_parameters;

  return 0;
}

