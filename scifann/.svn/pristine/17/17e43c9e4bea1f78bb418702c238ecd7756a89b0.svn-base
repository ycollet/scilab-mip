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
#include <api_scilab.h>
#include <MALLOC.h>
#include <Scierror.h>

int sci_fann_create(char * fname)
{
  int * pi_command_addr = NULL;
  int m_layers,  n_layers,  * pi_layers_addr = NULL;
  int * pi_conn_addr = NULL;
  char * Command = NULL;
  double * layers = NULL, conn = 0.0;
  unsigned int * ui_layers = NULL;
  int res, numLayers, i;
  struct fann * result_ann = NULL;
  SciErr _sciErr;

  if (Rhs<2)
    {
      Scierror(999,"%s usage: ann = %s(command,[layers ...])", fname, fname);
      return 0;
    }

  _sciErr = getVarAddressFromPosition(pvApiCtx, 1, &pi_command_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  getAllocatedSingleString(pvApiCtx,  pi_command_addr, &Command);

  _sciErr = getVarAddressFromPosition(pvApiCtx, 2, &pi_layers_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  _sciErr = getMatrixOfDouble(pvApiCtx, pi_layers_addr, &m_layers, &n_layers, &layers);

  if ((n_layers != 1) & (m_layers !=1))
    {
      Scierror(999,"%s: Layers must be a vector!",fname);
      return 0;
    }
  
  numLayers = m_layers * n_layers;
  ui_layers = (unsigned int *)MALLOC(numLayers*sizeof(unsigned int));
  for(i=0; i<numLayers; i++) ui_layers[i] = layers[i];

  if (strcmp(Command,"standard") == 0)
    {
      freeAllocatedSingleString(Command);

      // fann_create_standard_array  Just like fann_create_standard, but with an array of layer sizes instead of individual parameters.
      result_ann = fann_create_standard_array(numLayers,ui_layers);
      FREE(ui_layers);
      if (result_ann==NULL)
	{
	  Scierror(999,"%s: not able to create standard network\n",fname);
	  return 0;
	}
    }
  
  if (strcmp(Command,"sparse") == 0)
    {
      freeAllocatedSingleString(Command);

      // fann_create_sparse_array    Just like fann_create_sparse, but with an array of layer sizes instead of individual parameters.
      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_conn_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      getScalarDouble(pvApiCtx, pi_conn_addr, &conn);

      result_ann = fann_create_sparse_array(conn,numLayers,ui_layers);
      FREE(ui_layers);
      if (result_ann==NULL)
	{
	  Scierror(999,"%s: not able to create sparse network\n",fname);
	  return 0;
	}
    }

  if (strcmp(Command,"shortcut") == 0)
    {
      freeAllocatedSingleString(Command);

      // fann_create_shortcut_array  Just like fann_create_shortcut, but with an array of layer sizes instead of individual parameters.
      result_ann = fann_create_shortcut_array(numLayers,ui_layers);
      FREE(ui_layers);
      if (result_ann==NULL)
	{
	  Scierror(999,"%s: not able to create shortcut network\n",fname);
	  return 0;
	}
    }

  //Create the struct representing this ann in scilab
  res = createScilabFannStructFromCFannStruct(result_ann, Rhs + 1);
  if (res==-1) return 0;

  LhsVar(1) = Rhs + 1;

  return 0;
}

int sci_fann_destroy(char * fname)
{
  int res;
  struct fann * result_fann = NULL;

  if (Rhs!=1)
    {
      Scierror(999,"%s usage: %s(ann_in)", fname, fname);
      return 0;
    }

  res = detect_fannlist(1);
  if (res==-1) return 0;

  result_fann = createCFannStructFromScilabFannStruct(1,&res);
  if (res==-1) return 0;

  if (result_fann==NULL)
    {
      Scierror(999,"%s: Problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  fann_destroy(result_fann);

  return 0;
}

// fann_copy                   Creates a copy of a fann structure.
int sci_fann_copy(char * fname)
{
  int res;
  struct fann * result_ann = NULL;
  struct fann * result_ann_copy = NULL;

  if ((Rhs!=1)&&(Lhs!=1))
    {
      Scierror(999,"%s usage: ann_out = %s(ann_in)", fname, fname);
      return 0;
    }

  // Get the ann
  res = detect_fannlist(1);
  if (res==-1) return 0;
  
  result_ann = createCFannStructFromScilabFannStruct(1,&res);
  if (res==-1) return 0;
  
  if (result_ann==NULL)
    {
      Scierror(999,"%s: Problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  result_ann_copy = fann_copy(result_ann);
  
  res = createScilabFannStructFromCFannStruct(result_ann_copy,Rhs + 1);
  if (res==-1) return 0;

  LhsVar(1) = Rhs + 1;

  return 0;
}

// fann_run                    Will run input through the neural network, returning an array of outputs, the number of which being equal to
//                             the number of neurons in the output layer.
int sci_fann_run(char * fname)
{
  int m_inputs,  n_inputs,  * pi_inputs_addr = NULL;
  int m_outputs, n_outputs;
  int res, i;
  double * inputs = NULL, * outputs = NULL;
  struct fann * result_ann = NULL;
  fann_type * out_tmp = NULL;
  SciErr _sciErr;

  if ((Rhs!=2)&&(Lhs!=1))
    {
      Scierror(999,"%s usage: data_out = %s(ann_in, input)", fname, fname);
      return 0;
    }

  // Get the ann
  res = detect_fannlist(1);
  if (res==-1) return 0;
  
  result_ann = createCFannStructFromScilabFannStruct(1,&res);
  if (res==-1) return 0;

  if (result_ann==NULL)
    {
      Scierror(999,"%s: Problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  _sciErr = getVarAddressFromPosition(pvApiCtx, 2, &pi_inputs_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  _sciErr = getMatrixOfDouble(pvApiCtx, pi_inputs_addr, &m_inputs, &n_inputs, &inputs);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  
  if (result_ann->num_input!=m_inputs*n_inputs)
    {
      Scierror(999,"%s: wrong input dimension\n",fname);
      return 0;
    }

  out_tmp = fann_run(result_ann,(fann_type *)inputs);

  m_outputs = result_ann->num_output; n_outputs = 1;

  _sciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, m_outputs, n_outputs, &outputs);

  for(i=0; i<m_outputs; i++) outputs[i] = out_tmp[i];

  LhsVar(1) = Rhs + 1;

  return 0;
}

// fann_randomize_weights      Give each connection a random weight between min_weight and max_weight
int sci_fann_randomize_weights(char * fname)
{
  int * pi_minweight_addr = NULL;
  int * pi_maxweight_addr = NULL;
  int res;
  double minweight = 0.0, maxweight = 0.0;
  struct fann * result_ann = NULL;
  SciErr _sciErr;

  if ((Rhs!=3)&&(Lhs!=1))
    {
      Scierror(999,"%s usage: ann_out = %s(ann_in, weight_min, weight_max)", fname, fname);
      return 0;
    }

  // Get the ann
  res = detect_fannlist(1);
  if (res==-1) return 0;
  
  result_ann = createCFannStructFromScilabFannStruct(1,&res);
  if (res==-1) return 0;

  if (result_ann==NULL)
    {
      Scierror(999,"%s: Problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  _sciErr = getVarAddressFromPosition(pvApiCtx, 2, &pi_minweight_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  getScalarDouble(pvApiCtx, pi_minweight_addr, &minweight);

  _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_maxweight_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }
  getScalarDouble(pvApiCtx, pi_maxweight_addr, &maxweight);

  fann_randomize_weights(result_ann,(fann_type)minweight,(fann_type)maxweight);
  
  res = createScilabFannStructFromCFannStruct(result_ann, Rhs + 1);
  if (res==-1) return 0;

  LhsVar(1) = Rhs + 1;

  return 0;
}

// fann_init_weights           Initialize the weights using Widrow + Nguyen’s algorithm.
int sci_fann_init_weights(char * fname)
{
  int res;
  struct fann * result_ann = NULL;
  struct fann_train_data * result_ann_train = NULL;

  if ((Rhs!=2)&&(Lhs!=1))
    {
      Scierror(999,"%s usage: ann_out = %s(ann_in, data_train_in)", fname, fname);
      return 0;
    }

  // Get the ann
  res = detect_fannlist(1);
  if (res==-1) return 0;

  result_ann = createCFannStructFromScilabFannStruct(1,&res);
  if (res==-1) return 0;

  if (result_ann==NULL)
    {
      Scierror(999,"%s: Problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  // Get the train data
  res = detect_fanntraindatalist(2);
  if (res==-1) return 0;

  result_ann_train = createCFannTrainDataStructFromScilabFannTrainDataStruct(2,&res);

  fann_init_weights(result_ann, result_ann_train);
  
  res = createScilabFannStructFromCFannStruct(result_ann, Rhs + 1);
  if (res==-1) return 0;

  LhsVar(1) = Rhs + 1;

  return 0;
}

// fann_print_parameters           Display the ANN parameters
int sci_fann_print_parameters(char * fname)
{
  int res;
  struct fann * result_ann = NULL;

  if (Rhs!=1)
    {
      Scierror(999,"%s usage: %s(ann_in)", fname, fname);
      return 0;
    }

  // Get the ann
  res = detect_fannlist(1);
  if (res==-1) return 0;

  result_ann = createCFannStructFromScilabFannStruct(1,&res);

  if (result_ann==NULL)
    {
      Scierror(999,"%s: Problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  fann_print_parameters(result_ann);
  
  LhsVar(1) = 0;

  return 0;
}

// fann_print_connections           Display the ANN connections
int sci_fann_print_connections(char * fname)
{
  int res;
  struct fann * result_ann = NULL;

  if (Rhs!=1)
    {
      Scierror(999,"%s usage: %s(ann_in)", fname, fname);
      return 0;
    }

  // Get the ann
  res = detect_fannlist(1);
  if (res==-1) return 0;

  result_ann = createCFannStructFromScilabFannStruct(1,&res);

  if (result_ann==NULL)
    {
      Scierror(999,"%s: Problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  fann_print_connections(result_ann);
  
  LhsVar(1) = 0;

  return 0;
}

int sci_fann_print_debug(char * fname)
{
  int res;
  struct fann * result_ann = NULL;

  if (Rhs!=1)
    {
      Scierror(999,"%s usage: %s(ann_in)", fname, fname);
      return 0;
    }

  // Get the ann
  res = detect_fannlist(1);
  if (res==-1) return 0;

  result_ann = createCFannStructFromScilabFannStruct(1,&res);

  if (result_ann==NULL)
    {
      Scierror(999,"%s: Problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  printf("DEBUG: learning_rate = %f\n", result_ann->learning_rate);
  printf("DEBUG: num_input     = %d\n", result_ann->num_input);
  printf("DEBUG: num_output    = %d\n", result_ann->num_output);
  printf("DEBUG: total_neurons = %d\n", result_ann->total_neurons);
  
  LhsVar(1) = 0;

  return 0;
}
