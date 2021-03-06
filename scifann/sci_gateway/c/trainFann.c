/*
 * Scilab interface for the FANN library
 * Author: Dirk Gorissen <dirk.gorissen@ua.ac.be>
 * Author: Yann COLLETTE <ycollet@freesurf.fr>
 * Licence: GPL version 2 or later
 */

#include <helperFann.h>

#include <stdio.h>
#include <string.h>

#include <stack-c.h>
#include <MALLOC.h>
#include <api_scilab.h>
#include <Scierror.h>

//Calling syntax: [ann] = fann_train_matrix(ann,samples,values,[desired error],[max epochs]);
int sci_fann_train_matrix(char * fname)
{
  int m_xData, n_xData, * pi_xData_addr = NULL;
  int m_fData, n_fData, * pi_fData_addr = NULL;
  int * pi_desiredError_addr = NULL;
  int * pi_epochs_addr = NULL;
  int sRowLen, sColLen, vRowLen, vColLen, res = 0;
  double * xData = NULL, * fData = NULL, tmp_dbl;
  unsigned int numInputs, numOutputs;
  struct fann* ann = NULL;
  struct fann_train_data * data = NULL;
  double desiredError    = 1e-5;
  unsigned int maxEpochs = 5000;
  SciErr _sciErr;

  if ((Rhs<3)||(Rhs>5))
    {
      Scierror(999,"%s: usage fann_out = %s(fann_in, samples, values, [desired_error],[max_epochs])\n", fname, fname);
      return 0;
    }

  res = detect_fannlist(1);
  if (res==-1) return 0;

  //Create the network
  ann  = createCFannStructFromScilabFannStruct(1,&res);
  if (res==-1) return 0;

  if (ann==NULL)
    {
      Scierror(999,"%s: Problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  //Get the samples
  _sciErr = getVarAddressFromPosition(pvApiCtx, 2, &pi_xData_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }

  _sciErr = getMatrixOfDouble(pvApiCtx, pi_xData_addr, &m_xData, &n_xData, &xData);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }

  sRowLen = m_xData;
  sColLen = n_xData;
  
  //Get the values
  _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_fData_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }

  _sciErr = getMatrixOfDouble(pvApiCtx, pi_fData_addr, &m_fData, &n_fData, &fData);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }

  vRowLen = m_fData;
  vColLen = n_fData;
  
  if (sRowLen != vRowLen)
    {
      Scierror(999, "%s: The number of samples and values must be equal",fname);
      return 0;
    }

  if (Rhs >= 4)
    {
      //desired error passed
      _sciErr = getVarAddressFromPosition(pvApiCtx, 4, &pi_desiredError_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      getScalarDouble(pvApiCtx, pi_desiredError_addr, &desiredError);
    }
  
  if (Rhs >= 5)
    {
      //epochs passed
      _sciErr = getVarAddressFromPosition(pvApiCtx, 4, &pi_epochs_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      getScalarDouble(pvApiCtx, pi_desiredError_addr, &tmp_dbl);
      maxEpochs = tmp_dbl;
    }

  numInputs  = sColLen;
  numOutputs = vColLen;
  
  //if training data was passed
  if (sColLen > 0)
    {
      if (numInputs != sColLen)
	{
	  Scierror(999, "%s: The dimension of the passed samples does not match the input dimension of the network\n",fname);
	  return 0;
	}
    
      if (numOutputs != vColLen)
	{
	  Scierror(999, "%s: The dimension of the passed values does not match the output dimension of the network\n",fname);
	  return 0;
	}
    
      //Create the training data structure
      data = read_from_array(xData, fData, sRowLen, numInputs, numOutputs);
    
      if (data==NULL)
	{
	  Scierror(999,"%s: Problem while allocating memory for the training data set\n",fname);
	  return 0;
	}

      //train the network
      fann_train_on_data(ann,data,maxEpochs,maxEpochs,desiredError);
    }
  
  //Create the struct representing this ann in Scilab
  res = createScilabFannStructFromCFannStruct(ann, Rhs + 1);
  if (res==-1) return 0;

  LhsVar(1) = Rhs + 1;

  return 0;
}

// fann_create_train                Creates an empty training data struct.
int sci_fann_create_train(char * fname)
{
  int * pi_numdata_addr = NULL;
  int * pi_numinput_addr = NULL;
  int * pi_numoutput_addr = NULL;
  double numdata, numinput, numoutput;
  int res;
  struct fann_train_data * result_ann_train = NULL;
  SciErr _sciErr;

  // YC: for all the options, check the number of parameters

  if (Rhs!=3)
    {
      Scierror(999,"%s usage: fann_train = %s(num_data, num_input, num_output)\n", fname, fname);
      return 0;
    }

  _sciErr = getVarAddressFromPosition(pvApiCtx, 1, &pi_numdata_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  getScalarDouble(pvApiCtx, pi_numdata_addr, &numdata);

  _sciErr = getVarAddressFromPosition(pvApiCtx, 2, &pi_numinput_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  getScalarDouble(pvApiCtx, pi_numinput_addr, &numinput);

  _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_numoutput_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  getScalarDouble(pvApiCtx, pi_numoutput_addr, &numoutput);
  
  result_ann_train = fann_create_train(numdata,numinput,numoutput);

  res = createScilabFannTrainDataStructFromCFannTrainDataStruct(result_ann_train,Rhs + 1);
  if (res==-1) return 0;

  LhsVar(1) = Rhs + 1;

  return 0;
}

// fann_destroy_train               Destructs the training data and properly deallocates all of the associated data.
int sci_fann_destroy_train(char * fname)
{
  // Get the fann_train_data structure
  int res;
  struct fann_train_data * result_ann_train = NULL;

  if (Rhs!=1)
    {
      Scierror(999,"%s: usage %s(fann_data_in)\n", fname);
      return 0;
    }

  res = detect_fanntraindatalist(1);
  if (res==-1) return 0;

  result_ann_train = createCFannTrainDataStructFromScilabFannTrainDataStruct(1,&res);
  if (res==-1) return 0;

  fann_destroy_train(result_ann_train);

  return 0;
}

// fann_read_train_from_file        Reads a file that stores training data.
int sci_fann_read_train_from_file(char * fname)
{
  int * pi_name_addr = NULL;
  int res;
  char * Name = NULL;
  struct fann_train_data * result_ann_train = NULL;
  SciErr _sciErr;

  if ((Rhs!=1)&&(Lhs!=1))
    {
      Scierror(999,"%s: usage fann_data_out = %s(filename)\n", fname, fname);
      return 0;
    }

  _sciErr = getVarAddressFromPosition(pvApiCtx, 1, &pi_name_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  getAllocatedSingleString(pvApiCtx, pi_name_addr, &Name);

  result_ann_train = fann_read_train_from_file(Name);

  if (result_ann_train==NULL)
    {
      Scierror(999,"%s: problem while reading a training data file\n",fname);
      freeAllocatedSingleString(Name);
      return 0;
    }

  res = createScilabFannTrainDataStructFromCFannTrainDataStruct(result_ann_train, Rhs + 1);
  freeAllocatedSingleString(Name);
  if (res==-1) return 0;

  LhsVar(1) = Rhs + 1;

  return 0;
}

// fann_train                       Train one iteration with a set of inputs, and a set of desired outputs.
int sci_fann_train(char * fname)
{
  int m_inputs,  n_inputs,  * pi_inputs_addr = NULL;
  int m_outputs, n_outputs, * pi_outputs_addr = NULL;
  double * inputs = NULL, * outputs = NULL;
  int res;
  struct fann * result_ann = NULL;
  SciErr _sciErr;

  if ((Rhs!=3)&&(Lhs!=1))
    {
      Scierror(999,"%s: usage fann_out = %s(fann_in, inputs, outputs)\n", fname, fname);
      return 0;
    }

  // Get the fann structure
  res = detect_fannlist(1);
  if (res==-1) return 0;

  result_ann = createCFannStructFromScilabFannStruct(1,&res);
  if (res==-1) return 0;

  if (result_ann==NULL)
    {
      Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  // YC: size verification
  _sciErr = getVarAddressFromPosition(pvApiCtx, 2, &pi_inputs_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  _sciErr = getMatrixOfDouble(pvApiCtx, pi_inputs_addr, &m_inputs, &n_inputs, &inputs);

  _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_outputs_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  _sciErr = getMatrixOfDouble(pvApiCtx, pi_outputs_addr, &m_outputs, &n_outputs, &outputs);

  fann_train(result_ann, (fann_type *)inputs, (fann_type *)outputs);

  res = createScilabFannStructFromCFannStruct(result_ann, Rhs + 1);
  if (res==-1) return 0;

  LhsVar(1) = Rhs + 1;

  return 0;
}

// fann_test                        Test with a set of inputs, and a set of desired outputs.
int sci_fann_test(char * fname)
{
  int m_inputs,  n_inputs,  * pi_inputs_addr  = NULL;
  int m_outputs, n_outputs, * pi_outputs_addr = NULL;
  int m_results, n_results;
  double * inputs = NULL, * outputs = NULL;
  int res, i;
  double * tmp = NULL;
  fann_type * result = NULL;
  struct fann * result_ann = NULL;
  SciErr _sciErr;

  if ((Rhs!=3)&&(Lhs!=2))
    {
      Scierror(999,"%s: usage [fann_out, outputs] = %s(fann_in, inputs, outputs)\n", fname);
      return 0;
    }

  // Get the fann structure
  res = detect_fannlist(1);
  if (res==-1) return 0;

  result_ann = createCFannStructFromScilabFannStruct(1,&res);
  if (res==-1) return 0;

  if (result_ann==NULL)
    {
      Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  // YC: size verification
  _sciErr = getVarAddressFromPosition(pvApiCtx, 2, &pi_inputs_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  _sciErr = getMatrixOfDouble(pvApiCtx, pi_inputs_addr, &m_inputs, &n_inputs, &inputs);

  _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_outputs_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  _sciErr = getMatrixOfDouble(pvApiCtx, pi_outputs_addr, &m_outputs, &n_outputs, &outputs);

  // YC: check size

  result = fann_test(result_ann, (fann_type *)inputs, (fann_type *)outputs);
  tmp = (double *)MALLOC(m_outputs*n_outputs*sizeof(double));
  for(i=0;i<m_outputs*n_outputs;i++) tmp[i] = (double)result[i];

  m_results = m_outputs; n_results = n_outputs;
  _sciErr = createMatrixOfDouble(pvApiCtx, Rhs + 1, m_results, n_results, tmp);
  FREE(tmp);

  res = createScilabFannStructFromCFannStruct(result_ann,Rhs + 2);
  if (res==-1) return 0;
  
  LhsVar(1) = Rhs + 2;
  LhsVar(2) = Rhs + 1;

  return 0;
}

// fann_get_MSE                     Reads the mean square error from the network.
int sci_fann_get_MSE(char * fname)
{
  int res;
  struct fann * result_ann = NULL;

  if (Rhs!=1)
    {
      Scierror(999,"%s: usage MSE = %s(fann_in)\n", fname, fname);
      return 0;
    }

  // Get the fann structure
  res = detect_fannlist(1);
  if (res==-1) return 0;

  result_ann = createCFannStructFromScilabFannStruct(1,&res);
  if (res==-1) return 0;

  if (result_ann==NULL)
    {
      Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  createScalarDouble(pvApiCtx, Rhs + 1, fann_get_MSE(result_ann));

  LhsVar(1) = Rhs + 1;

  return 0;
}

// fann_get_bit_fail                The number of fail bits; means the number of output neurons which differ more than the bit 
//                                  fail limit (see fann_get_bit_fail_limit, fann_set_bit_fail_limit).
int sci_fann_get_bit_fail(char * fname)
{
  int res;
  struct fann * result_ann = NULL;

  if (Rhs!=1)
    {
      Scierror(999,"%s: usage bit_fail = %s(fann_in)\n", fname, fname);
      return 0;
    }

  // Get the fann structure
  res = detect_fannlist(1);
  if (res==-1) return 0;

  result_ann = createCFannStructFromScilabFannStruct(1,&res);
  if (res==-1) return 0;

  if (result_ann==NULL)
    {
      Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  res = (double)fann_get_bit_fail(result_ann);
  createScalarDouble(pvApiCtx, Rhs + 1, res);

  LhsVar(1) = Rhs + 1;

  return 0;
}

// fann_reset_MSE                   Resets the mean square error from the network.
int sci_fann_reset_MSE(char * fname)
{
  int res;
  struct fann * result_ann = NULL;

  if (Rhs!=1)
    {
      Scierror(999,"%s: usage fann_out = %s(fann_in)\n", fname, fname);
      return 0;
    }

  // Get the fann structure
  res = detect_fannlist(1);
  if (res==-1) return 0;

  result_ann = createCFannStructFromScilabFannStruct(1,&res);
  if (res==-1) return 0;

  if (result_ann==NULL)
    {
      Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  fann_reset_MSE(result_ann);

  // YC: is the output really necessary ?
  res = createScilabFannStructFromCFannStruct(result_ann,Rhs + 1);
  if (res==-1) return 0;

  LhsVar(1) = Rhs + 1;

  return 0;
}

// fann_get_output                   get the output i of given train data set
int sci_fann_get_output(char * fname)
{
  int res, i;
  int * pi_index_addr = NULL;
  double _tmp = 0.0, * output = NULL;
  unsigned int index = 1, num_output = 1;
  struct fann_train_data * result_ann_train = NULL;
  SciErr _sciErr;

  if (Rhs!=2)
    {
      Scierror(999,"%s: usage output = %s(fann_data_in, index)\n", fname, fname);
      return 0;
    }

  // Get the fann_train_data structure
  res = detect_fanntraindatalist(1);
  if (res==-1) return 0;

  result_ann_train = createCFannTrainDataStructFromScilabFannTrainDataStruct(1,&res);
  if (res==-1) return 0;

  if (result_ann_train==NULL)
    {
      Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  // Get the index
  _sciErr = getVarAddressFromPosition(pvApiCtx, 2, &pi_index_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  getScalarDouble(pvApiCtx, pi_index_addr, &_tmp);
  index = (unsigned int)_tmp;

  if ((index>result_ann_train->num_data) || (index<1))
    {
      Scierror(999,"%s: index must be between 1 and %d.\n", fname, result_ann_train->num_data);
      return 0;
    }

  _sciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, result_ann_train->num_output, 1, &output);

  for(i=0;i<result_ann_train->num_output; i++) output[i] = result_ann_train->output[index-1][i];

  LhsVar(1) = Rhs + 1;

  return 0;
}

// fann_get_input                   get the input i of given train data set
int sci_fann_get_input(char * fname)
{
  int res, i;
  int * pi_index_addr = NULL;
  int * pi_desired_error_addr = NULL;
  double _tmp = 0.0, * output = NULL;
  unsigned int index = 0;
  struct fann_train_data * result_ann_train = NULL;
  SciErr _sciErr;

  if (Rhs!=2)
    {
      Scierror(999,"%s: usage output = %s(fann_data_in, index)\n", fname, fname);
      return 0;
    }

  // Get the fann_train_data structure
  res = detect_fanntraindatalist(1);
  if (res==-1) return 0;

  result_ann_train = createCFannTrainDataStructFromScilabFannTrainDataStruct(1,&res);
  if (res==-1) return 0;

  if (result_ann_train==NULL)
    {
      Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  // Get the index
  _sciErr = getVarAddressFromPosition(pvApiCtx, 2, &pi_index_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  getScalarDouble(pvApiCtx, pi_index_addr, &_tmp);
  index = (unsigned int)_tmp;

  if ((index>result_ann_train->num_data) || (index<1))
    {
      Scierror(999,"%s: index must be between 1 and %d.\n", fname, result_ann_train->num_data);
      return 0;
    }

  _sciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, result_ann_train->num_input, 1, &output);
  for(i=0; i<result_ann_train->num_input; i++) output[i] = result_ann_train->input[index-1][i];

  LhsVar(1) = Rhs + 1;

  return 0;
}

// fann_set_output                   set the output i of given train data set
int sci_fann_set_output(char * fname)
{
  int res, i;
  int m_value, n_value, * pi_value_addr = NULL;
  int * pi_index_addr = NULL;
  double * value = NULL, _tmp = 0.0;
  unsigned int index = 0;
  struct fann_train_data * result_ann_train = NULL;
  SciErr _sciErr;

  if (Rhs!=3)
    {
      Scierror(999,"%s: usage fann_data_out = %s(fann_data_in, index, value)\n", fname, fname);
      return 0;
    }

  // Get the fann_train_data structure
  res = detect_fanntraindatalist(1);
  if (res==-1) return 0;

  result_ann_train = createCFannTrainDataStructFromScilabFannTrainDataStruct(1,&res);
  if (res==-1) return 0;

  if (result_ann_train==NULL)
    {
      Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  // Get the index
  _sciErr = getVarAddressFromPosition(pvApiCtx, 2, &pi_index_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  getScalarDouble(pvApiCtx, pi_index_addr, &_tmp);
  index = (unsigned int)_tmp;

  // Get the values
  _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_value_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  _sciErr = getMatrixOfDouble(pvApiCtx, pi_value_addr, &m_value, &n_value, &value);

  if (m_value*n_value!=result_ann_train->num_output)
    {
      Scierror(999,"%s: dimension of the output (%d) is not equal to the dimension of the output in the data set (%d)\n",fname, m_value*n_value, result_ann_train->num_output);
      return 0;
    }

  for(i=0; i<result_ann_train->num_output; i++) result_ann_train->output[index-1][i] = (fann_type)value[i];

  res = createScilabFannTrainDataStructFromCFannTrainDataStruct(result_ann_train, Rhs + 1);
  if (res==-1) return 0;

  LhsVar(1) = Rhs + 1;

  return 0;
}

// fann_set_input                   set the input i of given train data set
int sci_fann_set_input(char * fname)
{
  int res, i;
  int m_value, n_value, * pi_value_addr = NULL;
  int * pi_index_addr = NULL;
  double * value = NULL, _tmp = 0.0;
  unsigned int index = 0;
  struct fann_train_data * result_ann_train = NULL;
  SciErr _sciErr;

  if (Rhs!=3)
    {
      Scierror(999,"%s: usage fann_data_out = %s(fann_data_in, index, value)\n", fname, fname);
      return 0;
    }

  // Get the fann_train_data structure
  res = detect_fanntraindatalist(1);
  if (res==-1) return 0;

  result_ann_train = createCFannTrainDataStructFromScilabFannTrainDataStruct(1,&res);
  if (res==-1) return 0;

  if (result_ann_train==NULL)
    {
      Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  // Get the index
  _sciErr = getVarAddressFromPosition(pvApiCtx, 2, &pi_index_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  getScalarDouble(pvApiCtx, pi_index_addr, &_tmp);
  index = (unsigned int)_tmp;

  // Get the values
  _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_value_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  _sciErr = getMatrixOfDouble(pvApiCtx, pi_value_addr, &m_value, &n_value, &value);

  if (m_value*n_value!=result_ann_train->num_input)
    {
      Scierror(999,"%s: dimension of the input (%d) is not equal to the dimension of the input in the data set (%d)\n",fname, m_value*n_value, result_ann_train->num_input);
      return 0;
    }

  for(i=0; i<result_ann_train->num_input; i++) result_ann_train->input[index-1][i] = (fann_type)value[i];

  res = createScilabFannTrainDataStructFromCFannTrainDataStruct(result_ann_train, Rhs + 1);
  if (res==-1) return 0;

  LhsVar(1) = Rhs + 1;

  return 0;
}

// fann_train_on_data               Trains on an entire dataset, for a period of time.
int sci_fann_train_on_data(char * fname)
{
  int res;
  int * pi_epochs_addr = NULL;
  int * pi_desired_error_addr = NULL;
  double epochs = 0.0, desired_error = 0.0;
  struct fann * result_ann = NULL;
  struct fann_train_data * result_ann_train = NULL;
  SciErr _sciErr;

  if ((Rhs!=3)&&(Lhs!=1))
    {
      Scierror(999,"%s: fann_out = %s(fann_in, fann_data_in, epochs, desired_error)\n", fname, fname);
      return 0;
    }

  // Get the fann structure
  res = detect_fannlist(1);
  if (res==-1) return 0;

  result_ann = createCFannStructFromScilabFannStruct(1,&res);
  if (res==-1) return 0;

  if (result_ann==NULL)
    {
      Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  // Get the fann_train_data structure
  res = detect_fanntraindatalist(2);
  if (res==-1) return 0;

  result_ann_train = createCFannTrainDataStructFromScilabFannTrainDataStruct(2,&res);

  if (result_ann_train==NULL)
    {
      Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  // Test on the number of parameters
  _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_epochs_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  getScalarDouble(pvApiCtx, pi_epochs_addr, &epochs);

  _sciErr = getVarAddressFromPosition(pvApiCtx, 4, &pi_desired_error_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  getScalarDouble(pvApiCtx, pi_desired_error_addr, &desired_error);

  fann_train_on_data(result_ann, result_ann_train, (unsigned int)epochs, (unsigned int)(epochs/10), (float)desired_error);

  res = createScilabFannStructFromCFannStruct(result_ann,Rhs + 1);
  if (res==-1) return 0;

  LhsVar(1) = Rhs + 1;

  return 0;
}

// fann_train_on_file               Does the same as fann_train_on_data, but reads the training data directly from a file.
int sci_fann_train_on_file(char * fname)
{
  int res;
  int * pi_filename_addr = NULL;
  int * pi_epochs_addr = NULL;
  int * pi_desired_error_addr = NULL;
  double epochs = 0.0, desired_error = 0.0;
  char * Name = NULL;
  struct fann * result_ann = NULL;
  SciErr _sciErr;

  if ((Rhs!=4)&&(Lhs!=1))
    {
      Scierror(999,"%s: fann_out = %s(fann_in, filename, epochs, desired_error)\n", fname, fname);
      return 0;
    }

  // Get the fann structure
  res = detect_fannlist(1);
  if (res==-1) return 0;

  result_ann = createCFannStructFromScilabFannStruct(1,&res);
  if (res==-1) return 0;

  if (result_ann==NULL)
    {
      Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  // Test on the number of parameters
  // Get the fann_train_data structure
  _sciErr = getVarAddressFromPosition(pvApiCtx, 2, &pi_filename_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  getAllocatedSingleString(pvApiCtx, pi_filename_addr, &Name);

  _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_epochs_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  getScalarDouble(pvApiCtx, pi_epochs_addr, &epochs);

  _sciErr = getVarAddressFromPosition(pvApiCtx, 4, &pi_desired_error_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  getScalarDouble(pvApiCtx, pi_desired_error_addr, &desired_error);

  fann_train_on_file(result_ann, Name, (unsigned int)epochs, (unsigned int)(epochs/10), (float)desired_error);

  if (result_ann==NULL)
    {
      Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  res = createScilabFannStructFromCFannStruct(result_ann,Rhs + 1);
  if (res==-1) return 0;

  LhsVar(1) = Rhs + 1;

  return 0;
}

// fann_train_epoch                 Train one epoch with a set of training data.
int sci_fann_train_epoch(char * fname)
{
  int res;
  struct fann * result_ann = NULL;
  struct fann_train_data * result_ann_train = NULL;

  if ((Rhs!=2)&&(Lhs!=1))
    {
      Scierror(999,"%s: usage fann_out = %s(fann_in, fann_data_in)\n", fname, fname);
      return 0;
    }

  // Get the fann structure
  res = detect_fannlist(1);
  if (res==-1) return 0;

  result_ann = createCFannStructFromScilabFannStruct(1,&res);
  if (res==-1) return 0;

  if (result_ann==NULL)
    {
      Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  // Get the fann_train_data structure
  res = detect_fanntraindatalist(2);
  if (res==-1) return 0;

  result_ann_train = createCFannTrainDataStructFromScilabFannTrainDataStruct(2,&res);

  if (result_ann_train==NULL)
    {
      Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  fann_train_epoch(result_ann, result_ann_train);

  res = createScilabFannStructFromCFannStruct(result_ann, Rhs + 1);
  if (res==-1) return 0;

  LhsVar(1) = Rhs + 1;

  return 0;
}

// fann_test_data                   Test a set of training data and calculates the MSE for the training data.
int sci_fann_test_data(char * fname)
{
  int res;
  struct fann * result_ann = NULL;
  struct fann_train_data * result_ann_train = NULL;

  if ((Rhs!=2)&&(Lhs!=1))
    {
      Scierror(999,"%s: usage fann_out = %s(fann_in, fann_data_in)\n", fname, fname);
      return 0;
    }

  // Get the fann structure
  res = detect_fannlist(1);
  if (res==-1) return 0;

  result_ann = createCFannStructFromScilabFannStruct(1,&res);
  if (res==-1) return 0;

  if (result_ann==NULL)
    {
      Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  // Get the fann_train_data structure
  res = detect_fanntraindatalist(2);
  if (res==-1) return 0;

  result_ann_train = createCFannTrainDataStructFromScilabFannTrainDataStruct(2,&res);
  if (res==-1) return 0;

  if (result_ann_train==NULL)
    {
      Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  fann_test_data(result_ann, result_ann_train);

  res = createScilabFannStructFromCFannStruct(result_ann, Rhs + 1);
  if (res==-1) return 0;

  LhsVar(1) = Rhs + 1;

  return 0;
}

int sci_fann_setup_train_data(char * fname)
{
  int * pi_param_addr  = NULL;
  int * pi_newmin_addr = NULL;
  int * pi_newmax_addr = NULL;
  int * pi_new_input_min_addr  = NULL;
  int * pi_new_input_max_addr  = NULL;
  int * pi_new_output_min_addr = NULL;
  int * pi_new_output_max_addr = NULL;
  int m_weight, n_weight, * pi_weight_addr = NULL;
  int * pi_pos_addr           = NULL;
  int * pi_length_addr        = NULL;
  int * pi_filename_addr      = NULL;
  int * pi_decimal_point_addr = NULL;
  char * Param = NULL;
  int res, managed = 0;
  struct fann_train_data * result_ann_train   = NULL;
  struct fann_train_data * result_ann_train_2 = NULL;
  struct fann_train_data * result_ann_train_3 = NULL;
  struct fann * result_ann = NULL;
  SciErr _sciErr;

  if (Rhs<2)
    {
      Scierror(999,"%s usage: ann_train_data = %s(ann_train_data,'parameter_name')\n", fname, fname);
      return 0;
    }

  // Get the fann_train_data structure
  res = detect_fanntraindatalist(1);
  if (res==-1) return 0;

  result_ann_train = createCFannTrainDataStructFromScilabFannTrainDataStruct(1,&res);
  if (res==-1) return 0;

  if (result_ann_train==NULL)
    {
      Scierror(999,"%s: problem while creating the fann_train_data scilab structure\n",fname);
      return 0;
    }

  // Get the parameter name
  _sciErr = getVarAddressFromPosition(pvApiCtx, 2, &pi_param_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }

  getAllocatedSingleString(pvApiCtx, pi_param_addr, &Param);

  // fann_shuffle_train_data          Shuffles training data, randomizing the order.
  // usage: ann_train_data = setup_train_data(ann_train_data,'shuffle_train_data');
  if (strcmp(Param,"shuffle_train_data")==0)
    {
      freeAllocatedSingleString(Param);
      managed = 1;
      fann_shuffle_train_data(result_ann_train);
      res = createScilabFannTrainDataStructFromCFannTrainDataStruct(result_ann_train, Rhs + 1);
      if (res==-1) return 0;

      LhsVar(1) = Rhs + 1;

      return 0;
    }

  // fann_scale_train                 Scale input and output data based on previously calculated parameters.
  // usage: ann = setup_train_data(ann_train_data,'scale_train',ann);
  if (strcmp(Param,"scale_train")==0)
    {
      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage fann_data_out = %s(fann_data_in, '%s', fann_in)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}
      
      freeAllocatedSingleString(Param);
  
      // Get the ann
      res = detect_fannlist(3);
      if (res==-1) return 0;

      managed = 1;

      result_ann = createCFannStructFromScilabFannStruct(3,&res);
      if (res==-1) return 0;

      if (result_ann==NULL)
	{
	  Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
	  return 0;
	}

      fann_scale_train(result_ann, result_ann_train);

      res = createScilabFannStructFromCFannStruct(result_ann, Rhs + 1);
      if (res==-1) return 0;

      LhsVar(1) = Rhs + 1;

      return 0;
    }

  // fann_descale_train               Descale input and output data based on previously calculated parameters.
  // usage: ann = setup_train_data(ann_train_data,'descale_train',ann);
  if (strcmp(Param,"descale_train")==0)
    {
      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage fann_data_out = %s(fann_data_in, '%s', fann_in)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      freeAllocatedSingleString(Param);

      // Get the ann
      res = detect_fannlist(3);
      if (res==-1) return 0;

      managed = 1;

      result_ann = createCFannStructFromScilabFannStruct(3,&res);
      if (res==-1) return 0;

      if (result_ann==NULL)
	{
	  Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
	  return 0;
	}

      fann_descale_train(result_ann, result_ann_train);

      res = createScilabFannStructFromCFannStruct(result_ann,Rhs + 1);
      if (res==-1) return 0;

      LhsVar(1) = Rhs + 1;

      return 0;
    }

  // fann_set_input_scaling_params    Calculate input scaling parameters for future use based on training data.
  // usage: ann_train_data = setup_train_data(ann_train_data,'set_input_scaling_params',ann,new_input_min,new_input_max);
  if (strcmp(Param,"set_input_scaling_params")==0)
    {
      double newmin=0.0, newmax=0.0;
      
      if (Rhs!=5)
	{
	  Scierror(999,"%s: usage: ann = %s(ann_train_data,'%s',ann,new_input_min,new_input_max)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      freeAllocatedSingleString(Param);

      // Get the ann
      res = detect_fannlist(3);
      if (res==-1) return 0;

      managed = 1;

      result_ann = createCFannStructFromScilabFannStruct(3,&res);
      if (res==-1) return 0;

      if (result_ann==NULL)
	{
	  Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
	  return 0;
	}

      // YC: test min<max
      _sciErr = getVarAddressFromPosition(pvApiCtx, 4, &pi_newmin_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      getScalarDouble(pvApiCtx, pi_newmin_addr, &newmin);

      _sciErr = getVarAddressFromPosition(pvApiCtx, 5, &pi_newmax_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      getScalarDouble(pvApiCtx, pi_newmax_addr, &newmax);

      fann_set_input_scaling_params(result_ann, result_ann_train, (float)newmin, (float)newmax);

      res = createScilabFannStructFromCFannStruct(result_ann, Rhs + 1);
      if (res==-1) return 0;

      LhsVar(1) = Rhs + 1;

      return 0;
    }

  // fann_set_output_scaling_params   Calculate output scaling parameters for future use based on training data.
  // usage: ann = setup_train_data(ann_train_data,'set_output_scaling_params',ann,new_output_min,new_output_max);
  if (strcmp(Param,"set_output_scaling_params")==0)
    {
      double newmin = 0.0, newmax = 0.0;
      
      if (Rhs!=5)
	{
	  Scierror(999,"%s: usage: ann = %s(ann_train_data,'%s',ann,new_output_min,new_output_max)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      freeAllocatedSingleString(Param);

      // Get the ann
      res = detect_fannlist(3);
      if (res==-1) return 0;
      
      managed = 1;

      result_ann = createCFannStructFromScilabFannStruct(3,&res);
      if (res==-1) return 0;

      if (result_ann==NULL)
	{
	  Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 4, &pi_newmin_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      getScalarDouble(pvApiCtx, pi_newmin_addr, &newmin);

      _sciErr = getVarAddressFromPosition(pvApiCtx, 5, &pi_newmax_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      getScalarDouble(pvApiCtx, pi_newmax_addr, &newmax);

      fann_set_output_scaling_params(result_ann, result_ann_train, (float)newmin, (float)newmax);

      res = createScilabFannStructFromCFannStruct(result_ann, Rhs + 1);
      if (res==-1) return 0;

      LhsVar(1) = Rhs + 1;

      return 0;
    }

  // fann_set_scaling_params          Calculate input and output scaling parameters for future use based on training data.
  // usage: ann = setup_train_data(ann_train_data,'set_output_scaling_params',ann,new_input_min,new_input_max,new_output_min,new_output_max);
  if (strcmp(Param,"set_scaling_params")==0)
    {
      double input_min = 0.0, input_max = 0.0, output_min = 0.0, output_max = 0.0;
      
      if (Rhs!=7)
	{
	  Scierror(999,"%s: usage: ann = %s(ann_train_data,'%s',ann,new_input_min,new_input_max,new_output_min,new_output_max)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      freeAllocatedSingleString(Param);

      // Get the ann
      res = detect_fannlist(3);
      if (res==-1) return 0;
      
      result_ann = createCFannStructFromScilabFannStruct(3,&res);
      if (res==-1) return 0;

      if (result_ann==NULL)
	{
	  Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
	  return 0;
	}

      managed = 1;

      _sciErr = getVarAddressFromPosition(pvApiCtx, 4, &pi_new_input_min_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      getScalarDouble(pvApiCtx, pi_new_input_min_addr, &input_min);

      _sciErr = getVarAddressFromPosition(pvApiCtx, 5, &pi_new_input_max_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      getScalarDouble(pvApiCtx, pi_new_input_max_addr, &input_max);

      _sciErr = getVarAddressFromPosition(pvApiCtx, 6, &pi_new_output_min_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      getScalarDouble(pvApiCtx, pi_new_output_min_addr, &output_min);

      _sciErr = getVarAddressFromPosition(pvApiCtx, 7, &pi_new_output_max_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      getScalarDouble(pvApiCtx, pi_new_output_max_addr, &output_max);

      fann_set_scaling_params(result_ann, result_ann_train, (float)input_min, (float)input_max, (float)output_min, (float)output_max);

      res = createScilabFannStructFromCFannStruct(result_ann, Rhs + 1);
      if (res==-1) return 0;

      LhsVar(1) = Rhs + 1;

      return 0;
    }

  // fann_clear_scaling_params        Clears scaling parameters.
  // usage: ann = setup_train_data(ann_train_data,'clear_scaling_params',ann);
  if (strcmp(Param,"clear_scaling_params")==0)
    {
      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage fann_out = %s(fann_data_in, '%s', fann_in)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}
      
      freeAllocatedSingleString(Param);
      
      // Get the ann
      res = detect_fannlist(3);
      if (res==-1) return 0;
      
      result_ann = createCFannStructFromScilabFannStruct(3,&res);
      if (res==-1) return 0;
      
      if (result_ann==NULL)
	{
	  Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
	  return 0;
	}

      managed = 1;

      fann_clear_scaling_params(result_ann);

      res = createScilabFannStructFromCFannStruct(result_ann,Rhs + 1);
      if (res==-1) return 0;

      LhsVar(1) = Rhs + 1;

      return 0;
    }

  // fann_scale_input                 Scale data in input vector before feed it to ann based on previously calculated parameters.
  // usage: ann = setup_train_data(ann_train_data,'scale_input',ann,weight);
  if (strcmp(Param,"scale_input")==0)
    {
      double * weight = NULL;

      if (Rhs!=4)
	{
	  Scierror(999,"%s: usage fann_out = %s(fann_data_in, '%s', fann_in, weight)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      freeAllocatedSingleString(Param);

      // Get the ann
      res = detect_fannlist(3);
      if (res==-1) return 0;
      
      result_ann = createCFannStructFromScilabFannStruct(3,&res);
      if (res==-1) return 0;

      if (result_ann==NULL)
	{
	  Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
	  return 0;
	}

      managed = 1;
      
      _sciErr = getVarAddressFromPosition(pvApiCtx, 4, &pi_weight_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      _sciErr = getMatrixOfDouble(pvApiCtx, pi_weight_addr, &m_weight, &n_weight, &weight);
      
      // YC: il faut verifier la taille de weight == num_input
      fann_scale_input(result_ann,(fann_type *)weight);
      
      res = createScilabFannStructFromCFannStruct(result_ann, Rhs + 1);
      if (res==-1) return 0;

      LhsVar(1) = Rhs + 1;

      return 0;
    }
  
  // fann_scale_output                Scale data in output vector before feed it to ann based on previously calculated parameters.
  // usage: ann = setup_train_data(ann_train_data,'scale_output',ann,weight);
  if (strcmp(Param,"scale_output")==0)
    {
      double * weight = NULL;

      if (Rhs!=4)
	{
	  Scierror(999,"%s: usage fann_out = %s(fann_data_in, '%s', fann_in, weight)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      freeAllocatedSingleString(Param);
      
      // Get the ann
      res = detect_fannlist(3);
      if (res==-1) return 0;
      
      result_ann = createCFannStructFromScilabFannStruct(3,&res);
      if (res==-1) return 0;
      
      if (result_ann==NULL)
	{
	  Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
	  return 0;
	}

      managed = 1;
      
      _sciErr = getVarAddressFromPosition(pvApiCtx, 4, &pi_weight_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      _sciErr = getMatrixOfDouble(pvApiCtx, pi_weight_addr, &m_weight, &n_weight, &weight);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      // YC: il faut verifier la taille de weight == num_output
      fann_scale_output(result_ann,(fann_type *)weight);
      
      res = createScilabFannStructFromCFannStruct(result_ann, Rhs + 1);
      if (res==-1) return 0;

      LhsVar(1) = Rhs + 1;

      return 0;
    }
  
  // fann_descale_input               Scale data in input vector after get it from ann based on previously calculated parameters.
  // usage: ann = setup_train_data(ann_train_data,'descale_input',ann,weight);
  if (strcmp(Param,"descale_input")==0)
    {
      double * weight = NULL;

      if (Rhs!=4)
	{
	  Scierror(999,"%s: usage fann_out = %s(fann_data_in, '%s', fann_in, weight)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      freeAllocatedSingleString(Param);
      
      // Get the ann
      res = detect_fannlist(3);
      if (res==-1) return 0;
      
      result_ann = createCFannStructFromScilabFannStruct(3,&res);
      if (res==-1) return 0;

      if (result_ann==NULL)
	{
	  Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
	  return 0;
	}

      managed = 1;
      
      _sciErr = getVarAddressFromPosition(pvApiCtx, 4, &pi_weight_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      _sciErr = getMatrixOfDouble(pvApiCtx, pi_weight_addr, &m_weight, &n_weight, &weight);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      
      // YC: il faut verifier la taille de weight == num_input
      fann_descale_input(result_ann,(fann_type *)weight);
      
      res = createScilabFannStructFromCFannStruct(result_ann, Rhs + 1);
      if (res==-1) return 0;

      LhsVar(1) = Rhs + 1;

      return 0;
    }
  
  // YC: put descale output and some other functions in a separate function: setparameterstrain and getparameterstrain
  //     we can get one input, one output, modify one input, one output, scale, descale

  // fann_descale_output              Scale data in output vector after get it from ann based on previously calculated parameters.
  // usage: ann_train_data = setup_train_data(ann_train_data,'descale_output',ann,weight);
  if (strcmp(Param,"descale_output")==0)
    {
      double * weight = NULL;

      if (Rhs!=4)
	{
	  Scierror(999,"%s: usage fann_data_out = %s(fann_data_in, '%s', fann_in, weight)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      freeAllocatedSingleString(Param);
      
      // Get the ann
      res = detect_fannlist(3);
      if (res==-1) return 0;
      
      result_ann = createCFannStructFromScilabFannStruct(3,&res);
      if (res==-1) return 0;
      
      if (result_ann==NULL)
	{
	  Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
	  return 0;
	}

      managed = 1;
      
      _sciErr = getVarAddressFromPosition(pvApiCtx, 4, &pi_weight_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      _sciErr = getMatrixOfDouble(pvApiCtx, pi_weight_addr, &m_weight, &n_weight, &weight);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      
      // YC: il faut verifier la taille de weight == num_output
      fann_descale_output(result_ann,(fann_type *)weight);
      
      res = createScilabFannStructFromCFannStruct(result_ann, Rhs + 1);
      if (res==-1) return 0;

      LhsVar(1) = Rhs + 1;

      return 0;
    }
  
  // fann_scale_input_train_data      Scales the inputs in the training data to the specified range.
  // usage: ann_train_data = setup_train_data(ann_train_data,'scale_input_train_data',new_input_min,new_input_max);
  if (strcmp(Param,"scale_input_train_data")==0)
    {
      double input_min = 0.0, input_max = 0.0;

      if (Rhs!=4)
	{
	  Scierror(999,"%s: usage fann_data_out = %s(fann_data_in, '%s', new_input_min, new_input_max)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      freeAllocatedSingleString(Param);
      
      managed = 1;
      
      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_new_input_min_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      getScalarDouble(pvApiCtx, pi_new_input_min_addr, &input_min);

      _sciErr = getVarAddressFromPosition(pvApiCtx, 4, &pi_new_input_max_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      getScalarDouble(pvApiCtx, pi_new_input_max_addr, &input_max);
      
      fann_scale_input_train_data(result_ann_train,(fann_type)input_min,(fann_type )input_max);
      
      res = createScilabFannTrainDataStructFromCFannTrainDataStruct(result_ann_train, Rhs + 1);
      if (res==-1) return 0;

      LhsVar(1) = Rhs + 1;

      return 0;
    }
  
  // fann_scale_output_train_data     Scales the outputs in the training data to the specified range.
  // usage: ann_train_data = setup_train_data(ann_train_data,'scale_output_train_data',new_output_min,new_output_max);
  if (strcmp(Param,"scale_output_train_data")==0)
    {
      double output_min = 0.0, output_max = 0.0;

      if (Rhs!=4)
	{
	  Scierror(999,"%s: usage fann_data_out = %s(fann_data_in, '%s', new_output_min, new_output_max)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      freeAllocatedSingleString(Param);
      
      managed = 1;
      
      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_new_output_min_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      getScalarDouble(pvApiCtx, pi_new_output_min_addr, &output_min);

      _sciErr = getVarAddressFromPosition(pvApiCtx, 4, &pi_new_output_max_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      getScalarDouble(pvApiCtx, pi_new_output_max_addr, &output_max);

      fann_scale_output_train_data(result_ann_train,(fann_type)output_min,(fann_type )output_max);
      
      res = createScilabFannTrainDataStructFromCFannTrainDataStruct(result_ann_train, Rhs + 1);
      if (res==-1) return 0;

      LhsVar(1) = Rhs + 1;

      return 0;
    }
  
  // fann_scale_train_data            Scales the inputs and outputs in the training data to the specified range.
  // usage: ann_train_data = setup_train_data(ann_train_data,'scale_train_data',new_min,new_max);
  if (strcmp(Param,"scale_train_data")==0)
    {
      double newmin = 0.0, newmax = 0.0;

      if (Rhs!=4)
	{
	  Scierror(999,"%s: usage fann_data_out = %s(fann_data_in, '%s', new_min, new_max)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      freeAllocatedSingleString(Param);
      
      managed = 1;
      
      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_newmin_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      getScalarDouble(pvApiCtx, pi_newmin_addr, &newmin);

      _sciErr = getVarAddressFromPosition(pvApiCtx, 4, &pi_newmax_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      getScalarDouble(pvApiCtx, pi_newmax_addr, &newmax);

      fann_scale_train_data(result_ann_train,(fann_type)newmin,(fann_type )newmax);
      
      res = createScilabFannTrainDataStructFromCFannTrainDataStruct(result_ann_train, Rhs + 1);
      if (res==-1) return 0;

      LhsVar(1) = Rhs + 1;

      return 0;
    }
  
  // fann_merge_train_data            Merges the data from data1 and data2 into a new struct fann_train_data.
  // usage: ann_train_data_3 = setup_train_data(ann_train_data1,'merge_train_data',ann_train_data_2);
  if (strcmp(Param,"merge_train_data")==0)
    {
      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage fann_data_out = %s(fann_data_in, '%s', fann_data_in_2)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      freeAllocatedSingleString(Param);

      // Get the ann_train_data_2 structure
      res = detect_fanntraindatalist(3);
      if (res==-1) return 0;

      result_ann_train_2 = createCFannTrainDataStructFromScilabFannTrainDataStruct(3,&res);
      if (res==-1) return 0;
      
      if (result_ann_train==NULL)
	{
	  Scierror(999,"%s: problem while creating the fann_train_data scilab structure\n",fname);
	  return 0;
	}

      managed = 1;

      result_ann_train_3 = fann_merge_train_data(result_ann_train,result_ann_train_2);
      
      res = createScilabFannTrainDataStructFromCFannTrainDataStruct(result_ann_train_3, Rhs + 1);
      if (res==-1) return 0;

      LhsVar(1) = Rhs + 1;

      return 0;
    }
  
  // fann_duplicate_train_data        Returns an exact copy of a struct fann_train_data.
  // usage: ann_train_data_2 = setup_train_data(ann_train_data,'duplicata_train_data');
  if (strcmp(Param,"duplicate_train_data")==0)
    {
      freeAllocatedSingleString(Param);

      managed = 1;

      result_ann_train_2 = fann_duplicate_train_data(result_ann_train);
      
      res = createScilabFannTrainDataStructFromCFannTrainDataStruct(result_ann_train_2, Rhs + 1);
      if (res==-1) return 0;

      LhsVar(1) = Rhs + 1;

      return 0;
    }
  
  // fann_subset_train_data           Returns an copy of a subset of the struct fann_train_data, starting at position pos and length elements forward.
  // usage: ann_train_data = setup_train_data(ann_train_data,'subset_train_data',pos,length);
  if (strcmp(Param,"subset_train_data")==0)
    {
      double pos = 0.0, length = 0.0;

      if (Rhs!=4)
	{
	  Scierror(999,"%s: usage fann_data_out = %s(fann_data_in, '%s', pos, length)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      freeAllocatedSingleString(Param);

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_pos_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      getScalarDouble(pvApiCtx, pi_pos_addr, &pos);

      _sciErr = getVarAddressFromPosition(pvApiCtx, 4, &pi_length_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      getScalarDouble(pvApiCtx, pi_length_addr, &length);

      managed = 1;

      result_ann_train_2 = fann_subset_train_data(result_ann_train,pos,length);
      
      res = createScilabFannTrainDataStructFromCFannTrainDataStruct(result_ann_train_2, Rhs + 1);
      if (res==-1) return 0;

      LhsVar(1) = Rhs + 1;

      return 0;
    }
  
  // fann_length_train_data           Returns the number of training patterns in the struct fann_train_data.
  // usage: length = setup_train_data(ann_train_data,'length_train_data');
  if (strcmp(Param,"length_train_data")==0)
    {
      double tmp = 0.0;

      freeAllocatedSingleString(Param);

      managed = 1;

      tmp = fann_length_train_data(result_ann_train);
      createScalarDouble(pvApiCtx, Rhs + 1, tmp);
      
      LhsVar(1) = Rhs + 1;

      return 0;
    }
  
  // fann_num_input_train_data        Returns the number of inputs in each of the training patterns in the struct fann_train_data.
  // usage: length = setup_train_data(ann_train_data,'num_input_train_data');
  if (strcmp(Param,"num_input_train_data")==0)
    {
      double tmp = 0.0;

      freeAllocatedSingleString(Param);

      managed = 1;

      tmp = fann_num_input_train_data(result_ann_train);
      createScalarDouble(pvApiCtx, Rhs + 1, tmp);

      LhsVar(1) = Rhs + 1;

      return 0;
    }
  
  // fann_num_output_train_data       Returns the number of outputs in each of the training patterns in the struct fann_train_data.
  // usage: length = setup_train_data(ann_train_data,'num_output_train_data');
  if (strcmp(Param,"num_output_train_data")==0)
    {
      double tmp = 0.0;

      freeAllocatedSingleString(Param);

      managed = 1;

      tmp = fann_num_output_train_data(result_ann_train);
      createScalarDouble(pvApiCtx, Rhs + 1, tmp);

      LhsVar(1) = Rhs + 1;

      return 0;
    }
  
  // fann_save_train                  Save the training structure to a file, with the format as specified in fann_read_train_from_file
  // usage: setup_train_data(ann_train_data,'save_train',filename);
  if (strcmp(Param,"save_train")==0)
    {
      char * Filename = NULL;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: %s(fann_data_in, '%s', filename)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      freeAllocatedSingleString(Param);

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_filename_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      
      getAllocatedSingleString(pvApiCtx,  pi_filename_addr, &Filename);

      managed = 1;

      res = fann_save_train(result_ann_train,Filename);
      
      if (res==-1)
	{
	  Scierror(999,"%s: error while writing the file %s\n",fname,Filename);
	  freeAllocatedSingleString(Filename);
	  return 0;
	}

      freeAllocatedSingleString(Filename);
      return 0;
    }
  
  // fann_save_train_to_fixed         Saves the training structure to a fixed point data file.
  // usage: setup_train_data(ann_train_data,'save_train_to_fixed',filename,decimal_point);
  if (strcmp(Param,"save_train_to_fixed")==0)
    {
      char * Filename = NULL;
      double decimal_point = 0.0;

      if (Rhs!=4)
	{
	  Scierror(999,"%s: %s(fann_data_in, '%s', filename, decimal_point)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      freeAllocatedSingleString(Param);

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_filename_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      
      getAllocatedSingleString(pvApiCtx,  pi_filename_addr, &Filename);

      _sciErr = getVarAddressFromPosition(pvApiCtx, 4, &pi_decimal_point_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      getScalarDouble(pvApiCtx, pi_decimal_point_addr, &decimal_point);

      managed = 1;

      res = fann_save_train_to_fixed(result_ann_train,Filename,(unsigned int)decimal_point);
      
      if (res==-1)
	{
	  Scierror(999,"%s: error while writing the file %s\n",fname,Filename);
	  freeAllocatedSingleString(Filename);
	  return 0;
	}

      freeAllocatedSingleString(Filename);
      return 0;
    }
  
  if (!managed)
    {
      Scierror(999,"%s: error, wrong command: %s\n",fname, Param);
      freeAllocatedSingleString(Param);
      return 0;
    }
  
  return 0;
}

// fann_cascadetrain_on_data   Trains on an entire dataset, for a period of time using the Cascade2 training algorithm.
int sci_fann_cascadetrain_on_data(char * fname)
{
  int res;
  int * pi_neurons_addr       = NULL;
  int * pi_desired_error_addr = NULL;
  double neurons = 0.0, desired_error = 0.0;
  struct fann * result_ann = NULL;
  struct fann_train_data * result_ann_train = NULL;
  SciErr _sciErr;

  if (Rhs!=4)
    {
      Scierror(999,"%s: %s(fann_in, fann_data_in, neurons, desired_error)\n", fname, fname);
      return 0;
    }
  
  // Get the fann structure
  res = detect_fannlist(1);
  if (res==-1) return 0;

  result_ann = createCFannStructFromScilabFannStruct(1,&res);
  if (res==-1) return 0;

  if (result_ann==NULL)
    {
      Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  if (fann_get_network_type(result_ann)!=FANN_NETTYPE_SHORTCUT)
    {
      Scierror(999,"%s: error, the network type must be FANN_NETTYPE_SHORTCUT\n",fname);
      return 0;
    }

  if (fann_get_num_layers(result_ann)!=2)
    {
      Scierror(999,"%s: error, the network type must have no hidden layer\n",fname);
      return 0;
    }

  // Get the fann_train_data structure
  res = detect_fanntraindatalist(2);
  if (res==-1) return 0;

  result_ann_train = createCFannTrainDataStructFromScilabFannTrainDataStruct(2,&res);
  if (res==-1) return 0;

  if (result_ann_train==NULL)
    {
      Scierror(999,"%s: problem while creating the fann_train_data scilab structure\n",fname);
      return 0;
    }

  _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_neurons_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  getScalarDouble(pvApiCtx, pi_neurons_addr, &neurons);

  _sciErr = getVarAddressFromPosition(pvApiCtx, 4, &pi_desired_error_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  getScalarDouble(pvApiCtx, pi_desired_error_addr, &desired_error);

  fann_cascadetrain_on_data(result_ann, result_ann_train, (unsigned int)neurons, (unsigned int)(neurons/10), (float)desired_error);

  res = createScilabFannStructFromCFannStruct(result_ann, Rhs + 1);
  if (res==-1) return 0;

  LhsVar(1) = Rhs + 1;

  return 0;
}

// fann_cascadetrain_on_file   Does the same as fann_cascadetrain_on_data, but reads the training data directly from a file.
int sci_fann_cascadetrain_on_file(char * fname)
{
  int res;
  int * pi_neurons_addr       = NULL;
  int * pi_desired_error_addr = NULL;
  int * pi_filename_addr      = NULL;
  char * Filename             = NULL;
  double neurons = 0.0, desired_error = 0.0;
  struct fann * result_ann = NULL;
  SciErr _sciErr;

  if (Rhs!=4)
    {
      Scierror(999,"%s: %s(fann_in, filename, neurons, desired_error)\n", fname, fname);
      return 0;
    }

  // Get the fann structure
  res = detect_fannlist(1);
  if (res==-1) return 0;

  result_ann = createCFannStructFromScilabFannStruct(1,&res);
  if (res==-1) return 0;

  if (result_ann==NULL)
    {
      Scierror(999,"%s: problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  if (fann_get_network_type(result_ann)!=FANN_NETTYPE_SHORTCUT)
    {
      Scierror(999,"%s: error, the network type must be FANN_NETTYPE_SHORTCUT\n",fname);
      return 0;
    }

  if (fann_get_num_layers(result_ann)!=2)
    {
      Scierror(999,"%s: error, the network type must have no hidden layer\n",fname);
      return 0;
    }

  // Get the fann_train_data structure

  _sciErr = getVarAddressFromPosition(pvApiCtx, 2, &pi_filename_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  getAllocatedSingleString(pvApiCtx,  pi_filename_addr, &Filename);

  _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_neurons_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  getScalarDouble(pvApiCtx, pi_neurons_addr, &neurons);

  _sciErr = getVarAddressFromPosition(pvApiCtx, 4, &pi_desired_error_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  getScalarDouble(pvApiCtx, pi_desired_error_addr, &desired_error);

  fann_cascadetrain_on_file(result_ann, Filename, (unsigned int)neurons, (unsigned int)(neurons/10), (float)desired_error);

  freeAllocatedSingleString(Filename);

  res = createScilabFannStructFromCFannStruct(result_ann, Rhs + 1);
  if (res==-1) return 0;

  LhsVar(1) = Rhs + 1;

  return 0;
}

//Calling syntax: ann_train_data = fann_train_matrix(samples,values);
int sci_fann_traindata_from_matrix(char * fname)
{
  int m_xData, n_xData, * pi_xData_addr = NULL;
  int m_fData, n_fData, * pi_fData_addr = NULL;
  int sRowLen, sColLen, vRowLen, vColLen, res = 0;
  double * xData = NULL, * fData = NULL;
  unsigned int numInputs, numOutputs;
  struct fann_train_data* ann_train_data = NULL;
  SciErr _sciErr;

  if (Rhs!=2)
    {
      Scierror(999, "%s usage: fann_train_data = %s(samples, values)\n", fname, fname);
      return 0;
    }

  // Get the samples
  _sciErr = getVarAddressFromPosition(pvApiCtx, 1, &pi_xData_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  _sciErr = getMatrixOfDouble(pvApiCtx, pi_xData_addr, &m_xData, &n_xData, &xData);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }

  sRowLen = m_xData;
  sColLen = n_xData;
  
  // Get the values
  _sciErr = getVarAddressFromPosition(pvApiCtx, 2, &pi_fData_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  _sciErr = getMatrixOfDouble(pvApiCtx, pi_fData_addr, &m_fData, &n_fData, &fData);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }

  vRowLen = m_fData;
  vColLen = n_fData;
  
  if (sRowLen != vRowLen)
    {
      Scierror(999, "%s: The number of samples and values must be equal.\n",fname);
      return 0;
    }

  numInputs  = sColLen;
  numOutputs = vColLen;

  // if training data was passed
  if (sColLen > 0)
    {
      // YC: test to be checked. Above we have numInputs = sColLen and numOutputs = vColLen
      // a similar test exists above
      if (numInputs != sColLen)
	{
	  Scierror(999, "%s: The dimension of the passed samples does not match the input dimension of the network\n",fname);
	  return 0;
	}
    
      if (numOutputs != vColLen)
	{
	  Scierror(999, "%s: The dimension of the passed values does not match the output dimension of the network\n",fname);
	  return 0;
	}
    
      // Create the training data structure
      ann_train_data = read_from_array(xData, fData, sRowLen, numInputs, numOutputs);
    
      if (ann_train_data==NULL)
	{
	  Scierror(999,"%s: Problem while allocating memory for the training data set\n",fname);
	  return 0;
	}
    }
  
  // Create the struct representing this ann in Scilab
  res = createScilabFannTrainDataStructFromCFannTrainDataStruct(ann_train_data, Rhs + 1);
  if (res==-1) return 0;

  LhsVar(1) = Rhs + 1;

  return 0;
}
