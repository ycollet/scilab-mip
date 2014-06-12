/*
 * Scilab interface for the FANN library
 * Author: Dirk Gorissen <dirk.gorissen@ua.ac.be>
 * Author: Yann COLLETTE <ycollet@freesurf.fr>
 * Licence: GPL version 2 or later
 */

#include <helperFann.h>

#include <stdio.h>

#include <stack-c.h>
#include <api_scilab.h>
#include <Scierror.h>

//Calling syntax: [values] = fann_test_matrix(ann,samples);
int sci_fann_test_matrix(char * fname)
{
  int m_xData, n_xData, * pi_xData_addr = NULL;
  int m_fData, n_fData;
  double * xData = NULL, * fData = NULL;
  int res;
  unsigned int numInputs;
  unsigned int numOutputs;
  int sRowLen;
  int sColLen;
  struct fann* ann = NULL;
  SciErr _sciErr;

  if (Rhs!=2)
    {
      Scierror(999,"%s usage: values = %s(ann, samples)\n", fname, fname);
      return 0;
    }

  res = detect_fannlist(1);
  if (res==-1) return 0;

  ann = createCFannStructFromScilabFannStruct(1,&res);
  if (res==-1) return 0;

  if (ann==NULL)
    {
      Scierror(999,"%s: Problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  numInputs  = ann->num_input;
  numOutputs = ann->num_output;

  //Get the samples
  _sciErr = getVarAddressFromPosition(pvApiCtx, 2, &pi_xData_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  
  _sciErr = getMatrixOfDouble(pvApiCtx, pi_xData_addr, &m_xData, &n_xData, &xData);
  sRowLen = m_xData;
  sColLen = n_xData;
  
  if (numInputs!=sColLen)
    {
      Scierror(999,"%s: The network input dimension does not match the dimension of the passed samples!",fname);
      return 0;
    }

  m_fData = 1; n_fData = numOutputs;
  _sciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, m_fData, n_fData, &fData);
  
  //evaluate the network on the given samples
  evaluateNetwork(ann, xData, fData, sRowLen);
  
  LhsVar(1) = Rhs + 1;
  
  return 0;
}
