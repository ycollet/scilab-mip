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

//--------------------------------------------------------------------------------------------------------
//Calling syntax: value = fann_get_parameters(ann,parameter_name,...);
int sci_fann_get_parameters(char * fname)
{
  int * pi_layer_addr  = NULL;
  int * pi_neuron_addr = NULL;
  int * pi_param_addr  = NULL;
  int m_value,  n_value;
  int res, i, managed = 0;
  unsigned int value_uint, nb_layers, nb_conn, afunc_count;
  double value_double;
  unsigned int * value_uint_ptr = NULL;
  double * value_dbl_ptr        = NULL;
  const char * value_char       = NULL;
  double * step                 = NULL;
  char * Param                  = NULL;
  enum fann_activationfunc_enum * afunc = NULL;
  struct fann_connection * ann_conn     = NULL;
  struct fann* ann          = NULL;
  static char ** afunc_list = NULL;
  SciErr _sciErr;

  if (Rhs<2)
    {
      Scierror(999,"%s usage: value = %s(ann,parameter_name)\n", fname, fname);
      return 0;
    }

  // Get the parameter name
  _sciErr = getVarAddressFromPosition(pvApiCtx, 2, &pi_param_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }

  getAllocatedSingleString(pvApiCtx,  pi_param_addr, &Param);

  // Get the Fann structure
  res = detect_fannlist(1);
  if (res==-1) return 0;

  ann = createCFannStructFromScilabFannStruct(1,&res);
  if (res==-1) return 0;

  if (ann==NULL)
    {
      Scierror(999,"%s: Problem while creating the fann scilab structure\n",fname);
      return 0;
    }

  // fann_get_num_input          Get the number of input neurons.
  if (strcmp(Param,"num_input")==0)
    {
      managed = 1;
      value_uint = fann_get_num_input(ann);
      createScalarDouble(pvApiCtx, Rhs + 1, (double)value_uint);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_num_output         Get the number of output neurons.
  if (strcmp(Param,"num_output")==0)
    {
      managed = 1;
      value_uint = fann_get_num_output(ann);
      createScalarDouble(pvApiCtx, Rhs + 1, (double)value_uint);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_total_neurons      Get the total number of neurons in the entire network.
  if (strcmp(Param,"total_neurons")==0)
    {
      managed = 1;
      value_uint = fann_get_total_neurons(ann);
      createScalarDouble(pvApiCtx, Rhs + 1, (double)value_uint);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_total_connections  Get the total number of connections in the entire network.
  if (strcmp(Param,"total_connections")==0)
    {
      managed = 1;
      value_uint = fann_get_total_connections(ann);
      createScalarDouble(pvApiCtx, Rhs + 1, (double)value_uint);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_network_type       Get the type of neural network it was created as.
  if (strcmp(Param,"network_type")==0)
    {
      managed = 1;
      value_char = FANN_NETTYPE_NAMES[fann_get_network_type(ann)];
      createSingleString(pvApiCtx, Rhs + 1, value_char);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_connection_rate    Get the connection rate used when the network was created
  if (strcmp(Param,"connection_rate")==0)
    {
      managed = 1;
      value_double = (double)fann_get_connection_rate(ann);
      createScalarDouble(pvApiCtx, Rhs + 1, value_double);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_num_layers         Get the number of layers in the network
  if (strcmp(Param,"num_layers")==0)
    {
      managed = 1;
      value_uint = fann_get_num_layers(ann);
      createScalarDouble(pvApiCtx, Rhs + 1, (double)value_uint);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_layer_array        Get the number of neurons in each layer in the network.
  if (strcmp(Param,"layer_array")==0)
    {
      managed = 1;
      nb_layers = fann_get_num_layers(ann);
      value_uint_ptr = (unsigned int *)MALLOC(nb_layers*sizeof(unsigned int));
      value_dbl_ptr  = (double *)MALLOC(nb_layers*sizeof(double));
      fann_get_layer_array(ann,value_uint_ptr);
      for(i=0; i<nb_layers;i++) value_dbl_ptr[i] = (double)value_uint_ptr[i];
      
      _sciErr = createMatrixOfDouble(pvApiCtx, Rhs + 1, nb_layers, 1, value_dbl_ptr);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);
      FREE(value_uint_ptr);
      FREE(value_dbl_ptr);

      return 0;
    }

  // fann_get_bias_array         Get the number of bias in each layer in the network.
  if (strcmp(Param,"bias_array")==0)
    {
      managed = 1;
      nb_layers = fann_get_num_layers(ann);
      value_uint_ptr = (unsigned int *)MALLOC(nb_layers*sizeof(unsigned int));
      value_dbl_ptr = (double *)MALLOC(nb_layers*sizeof(double));
      fann_get_bias_array(ann,value_uint_ptr);
      for(i=0; i<nb_layers;i++) value_dbl_ptr[i] = (double)value_uint_ptr[i];

      _sciErr = createMatrixOfDouble(pvApiCtx, Rhs + 1, nb_layers, 1, value_dbl_ptr);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);
      FREE(value_uint_ptr);
      FREE(value_dbl_ptr);

      return 0;
    }

  // fann_get_connection_array   Get the connections in the network.
  if (strcmp(Param,"connection_array")==0)
    {
      managed = 1;
      nb_conn = fann_get_total_connections(ann);
      ann_conn = (struct fann_connection *)MALLOC(nb_conn*sizeof(struct fann_connection));
      fann_get_connection_array(ann,ann_conn);

      _sciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, nb_conn, 3, &value_dbl_ptr);
      for(i=0;i<(int)nb_conn;i++) 
	{
	  value_dbl_ptr[i*3+0] = (double)ann_conn[i].from_neuron;
	  value_dbl_ptr[i*3+1] = (double)ann_conn[i].to_neuron;
	  value_dbl_ptr[i*3+2] = (double)ann_conn[i].weight;
	}

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);
      FREE(ann_conn);

      return 0;
    }

#ifdef FIXEDFANN
  // fann_get_decimal_point      Returns the position of the decimal point in the ann.
  if (strcmp(Param,"decimal_point")==0)
    {
      managed = 1;
      value_uint = fann_get_decimal_point(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, (double)value_uint);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_multiplier         returns the multiplier that fix point data is multiplied with.
  if (strcmp(Param,"multiplier")==0)
    {
      managed = 1;
      value_uint = fann_get_multiplier(ann);
      
      createScalarDouble(pvApiCtx, Rhs + 1, (double)value_uint);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }
#endif

  // fann_get_training_algorithm                   Return the training algorithm as described by fann_train_enum.
  if (strcmp(Param,"training_algorithm")==0)
    {
      managed = 1;
      value_char = FANN_TRAIN_NAMES[fann_get_training_algorithm(ann)];

      createSingleString(pvApiCtx, Rhs + 1, value_char);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_learning_rate                        Return the learning rate.
  if (strcmp(Param,"learning_rate")==0)
    {
      managed = 1;
      value_double = (double)fann_get_learning_rate(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, value_double);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_learning_momentum                    Get the learning momentum.
  if (strcmp(Param,"learning_momentum")==0)
    {
      managed = 1;
      value_double = (double)fann_get_learning_momentum(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, value_double);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_activation_function                  Get the activation function for neuron number neuron in layer number layer, counting the input layer as layer 0.
  if (strcmp(Param,"activation_function")==0)
    {
      double layer = 0.0, neuron = 0.0;
      managed = 1;
      
      if (Rhs!=4)
	{
	  Scierror(999,"%s: usage activation_name = %s(fann_in, %s, index_layer, index_neuron)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      // Test on the number of parameters 
      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_layer_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      getScalarDouble(pvApiCtx, pi_layer_addr, &layer);

      _sciErr = getVarAddressFromPosition(pvApiCtx, 4, &pi_neuron_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      getScalarDouble(pvApiCtx, pi_neuron_addr, &neuron);

      createSingleString(pvApiCtx, Rhs + 1, FANN_ACTIVATIONFUNC_NAMES[fann_get_activation_function(ann,layer,neuron)]);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_activation_steepness                 Get the activation steepness for neuron number neuron in layer number layer, counting the input layer as layer 0.
  if (strcmp(Param,"activation_steepness")==0)
    {
      double layer = 0.0, neuron = 0.0;
      managed = 1;

      if (Rhs!=4)
	{
	  Scierror(999,"%s: usage activation_steepness = %s(fann_in, %s, index_layer, index_neuron)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      // Test on the number of parameters 
      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_layer_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      getScalarDouble(pvApiCtx, pi_layer_addr, &layer);

      _sciErr = getVarAddressFromPosition(pvApiCtx, 4, &pi_neuron_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      getScalarDouble(pvApiCtx, pi_neuron_addr, &neuron);

      value_double = fann_get_activation_steepness(ann,layer,neuron);

      createScalarDouble(pvApiCtx, Rhs + 1, value_double);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_train_error_function                 Returns the error function used during training.
  if (strcmp(Param,"train_error_function")==0)
    {
      managed = 1;

      createSingleString(pvApiCtx, Rhs + 1, FANN_ERRORFUNC_NAMES[fann_get_train_error_function(ann)]);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_train_stop_function                  Returns the the stop function used during training.
  if (strcmp(Param,"train_stop_function")==0)
    {
      managed = 1;

      createSingleString(pvApiCtx, Rhs + 1, FANN_STOPFUNC_NAMES[fann_get_train_stop_function(ann)]);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_bit_fail_limit                       Returns the bit fail limit used during training.
  if (strcmp(Param,"bit_fail_limit")==0)
    {
      managed = 1;
      value_double = fann_get_bit_fail_limit(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, value_double);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_quickprop_decay                      The decay is a small negative valued number which is the factor that the weights should become smaller
  //                                               in each iteration during quickprop training.
  if (strcmp(Param,"quickprop_decay")==0)
    {
      managed = 1;
      value_double = (double)fann_get_bit_fail_limit(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, value_double);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_rprop_increase_factor                The increase factor is a value larger than 1, which is used to increase the step-size during RPROP training.
  if (strcmp(Param,"rprop_increase_factor")==0)
    {
      managed = 1;
      value_double = (double)fann_get_rprop_increase_factor(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, value_double);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_quickprop_mu                         The mu factor is used to increase and decrease the step-size during quickprop training.
  if (strcmp(Param,"quickprop_mu")==0)
    {
      managed = 1;
      value_double = (double)fann_get_quickprop_mu(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, value_double);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_rprop_decrease_factor                The decrease factor is a value smaller than 1, which is used to decrease the step-size during RPROP training.
  if (strcmp(Param,"rprop_decrease_factor")==0)
    {
      managed = 1;
      value_double = (double)fann_get_rprop_decrease_factor(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, value_double);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_rprop_delta_min                      The minimum step-size is a small positive number determining how small the minimum step-size may be.
  if (strcmp(Param,"rprop_delta_min")==0)
    {
      managed = 1;
      value_double = (double)fann_get_rprop_delta_min(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, value_double);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_rprop_delta_max                      The maximum step-size is a positive number determining how large the maximum step-size may be.
  if (strcmp(Param,"rprop_delta_max")==0)
    {
      managed = 1;
      value_double = (double)fann_get_rprop_delta_max(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, value_double);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_rprop_delta_zero                     The initial step-size is a positive number determining the initial step size.
  if (strcmp(Param,"rprop_delta_zero")==0)
    {
      managed = 1;
      value_double = (double)fann_get_rprop_delta_zero(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, value_double);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_sarprop_weight_decay_shift           The sarprop weight decay shift.
  if (strcmp(Param,"sarprop_weight_decay_shift")==0)
    {
      managed = 1;
      value_double = (double)fann_get_sarprop_weight_decay_shift(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, value_double);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_sarprop_step_error_threshold_factor  The sarprop step error threshold factor.
  if (strcmp(Param,"sarprop_step_error_threshold_factor")==0)
    {
      managed = 1;
      value_double = (double)fann_get_sarprop_step_error_threshold_factor(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, value_double);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_sarprop_step_error_shift             The get sarprop step error shift.
  if (strcmp(Param,"sarprop_step_error_shift")==0)
    {
      managed = 1;
      value_double = (double)fann_get_sarprop_step_error_shift(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, value_double);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_sarprop_temperature                  The sarprop weight decay shift.
  if (strcmp(Param,"sarprop_step_error_shift")==0)
    {
      managed = 1;
      value_double = (double)fann_get_sarprop_temperature(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, value_double);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_cascade_output_change_fraction       The cascade output change fraction is a number between 0 and 1 determining how large a fraction 
  //                                               the fann_get_MSE value should change within fann_get_cascade_output_stagnation_epochs during training 
  //                                               of the output connections, in order for the training not to stagnate.
  if (strcmp(Param,"cascade_output_change_fraction")==0)
    {
      managed = 1;
      value_double = (double)fann_get_cascade_output_change_fraction(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, value_double);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_cascade_output_stagnation_epochs     The number of cascade output stagnation epochs determines the number of epochs training is allowed to 
  //                                               continue without changing the MSE by a fraction of fann_get_cascade_output_change_fraction.
  if (strcmp(Param,"cascade_output_stagnation_epochs")==0)
    {
      managed = 1;
      value_uint = fann_get_cascade_output_stagnation_epochs(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, (double)value_uint);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_cascade_candidate_change_fraction    The cascade candidate change fraction is a number between 0 and 1 determining how large a fraction
  //                                               the fann_get_MSE value should change within fann_get_cascade_candidate_stagnation_epochs during training
  //                                               of the candidate neurons, in order for the training not to stagnate.
  if (strcmp(Param,"cascade_candidate_change_fraction")==0)
    {
      managed = 1;
      value_double = (double)fann_get_cascade_candidate_change_fraction(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, value_double);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_cascade_candidate_stagnation_epochs  The number of cascade candidate stagnation epochs determines the number of epochs training is
  //                                               allowed to continue without changing the MSE by a fraction of fann_get_cascade_candidate_change_fraction.
  if (strcmp(Param,"cascade_candidate_stagnation_epochs")==0)
    {
      managed = 1;
      value_uint = fann_get_cascade_candidate_stagnation_epochs(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, (double)value_uint);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_cascade_weight_multiplier            The weight multiplier is a parameter which is used to multiply the weights from the candidate 
  //                                               neuron before adding the neuron to the neural network.
  if (strcmp(Param,"cascade_weight_multiplier")==0)
    {
      managed = 1;
      value_double = fann_get_cascade_weight_multiplier(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, value_double);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_cascade_candidate_limit              The candidate limit is a limit for how much the candidate neuron may be trained.
  if (strcmp(Param,"cascade_candidate_limit")==0)
    {
      managed = 1;
      value_double = fann_get_cascade_candidate_limit(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, value_double);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_cascade_max_out_epochs               The maximum out epochs determines the maximum number of epochs the output connections may be
  //                                               trained after adding a new candidate neuron.
  if (strcmp(Param,"cascade_max_out_epochs")==0)
    {
      managed = 1;
      value_uint = fann_get_cascade_max_out_epochs(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, (double)value_uint);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_cascade_min_out_epochs               The minimum out epochs determines the minimum number of epochs the output connections must be
  //                                               trained after adding a new candidate neuron.
  if (strcmp(Param,"cascade_min_out_epochs")==0)
    {
      managed = 1;
      value_uint = fann_get_cascade_min_out_epochs(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, (double)value_uint);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_cascade_max_cand_epochs              The maximum candidate epochs determines the maximum number of epochs the input connections to
  //                                               the candidates may be trained before adding a new candidate neuron.
  if (strcmp(Param,"cascade_max_cand_epochs")==0)
    {
      managed = 1;
      value_uint = fann_get_cascade_max_cand_epochs(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, (double)value_uint);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_cascade_min_cand_epochs              The minimum candidate epochs determines the minimum number of epochs the input connections to
  //                                               the candidates may be trained before adding a new candidate neuron.
  if (strcmp(Param,"cascade_min_cand_epochs")==0)
    {
      managed = 1;
      value_uint = fann_get_cascade_min_cand_epochs(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, (double)value_uint);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_cascade_num_candidates               The number of candidates used during training (calculated by multiplying 
  //                                               fann_get_cascade_activation_functions_count, fann_get_cascade_activation_steepnesses_count and 
  //                                               fann_get_cascade_num_candidate_groups).
  if (strcmp(Param,"cascade_num_candidates")==0)
    {
      managed = 1;
      value_uint = fann_get_cascade_num_candidates(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, (double)value_uint);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_cascade_activation_functions_count   The number of activation functions in the fann_get_cascade_activation_functions array.
  if (strcmp(Param,"cascade_activation_functions_count")==0)
    {
      managed = 1;
      value_uint = fann_get_cascade_activation_functions_count(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, (double)value_uint);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_cascade_activation_functions         The cascade activation functions array is an array of the different activation functions used by the candidates.
  if (strcmp(Param,"cascade_activation_functions")==0)
    {
      managed = 1;
      afunc_count = fann_get_cascade_activation_functions_count(ann);
      afunc = fann_get_cascade_activation_functions(ann);

      m_value = afunc_count; n_value = 1;
      afunc_list = NULL;
      afunc_list = (char **)MALLOC(afunc_count*sizeof(char *));

      for(i=0; i<(int)afunc_count; i++)
	{
	  afunc_list[i] = (char *)MALLOC(strlen(FANN_ACTIVATIONFUNC_NAMES[afunc[i]])*sizeof(char));
	  strcpy(afunc_list[i],FANN_ACTIVATIONFUNC_NAMES[afunc[i]]);
	}

      _sciErr = createMatrixOfString(pvApiCtx, Rhs + 1, m_value, n_value, afunc_list);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);
      freeAllocatedMatrixOfString(m_value, n_value, afunc_list);

      return 0;
    }

  // fann_get_cascade_activation_steepnesses_count The number of activation steepnesses in the fann_get_cascade_activation_functions array.
  if (strcmp(Param,"cascade_activation_steepnesses_count")==0)
    {
      managed = 1;
      value_uint = fann_get_cascade_activation_steepnesses_count(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, (double)value_uint);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_cascade_activation_steepnesses       The cascade activation steepnesses array is an array of the different activation functions used by the candidates.
  if (strcmp(Param,"cascade_activation_steepnesses")==0)
    {
      managed = 1;
      value_uint = fann_get_cascade_activation_steepnesses_count(ann);
      step = (double *)fann_get_cascade_activation_steepnesses(ann);

      m_value = value_uint; n_value = 1;
      _sciErr = allocMatrixOfDouble(pvApiCtx, Rhs + 1, m_value, n_value, &value_dbl_ptr);
      for(i=0; i<(int)value_uint; i++) value_dbl_ptr[i] = step[i];

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  // fann_get_cascade_num_candidate_groups         The number of candidate groups is the number of groups of identical candidates which will be used during training.
  if (strcmp(Param,"cascade_num_candidate_groups")==0)
    {
      managed = 1;
      value_uint = fann_get_cascade_num_candidate_groups(ann);

      createScalarDouble(pvApiCtx, Rhs + 1, (double)value_uint);

      LhsVar(1) = Rhs + 1;

      freeAllocatedSingleString(Param);

      return 0;
    }

  if (!managed)
    {
      Scierror(999,"%s: error, wrong command: %s\n",fname,Param);
      freeAllocatedSingleString(Param);
      return 0;
    }

  freeAllocatedSingleString(Param);

  return 0;
}

//--------------------------------------------------------------------------------------------------------

#define FANN_SEARCH_NAME(LIST,LEN_LIST,MESSAGE) \
  Index = -1;				        \
  for(i=0;i<LEN_LIST;i++)			\
    {						\
      if (strcmp(Name,LIST[i])==0)		\
	{					\
	  Index = i;				\
	  break;				\
	}					\
    }						\
  if (Index==-1)				\
    {						\
      Scierror(999,MESSAGE);			\
      return 0;					\
    }

//Calling syntax: ann = fann_set_parameters(ann,parameter_name);
int sci_fann_set_parameters(char * fname)
{
  int * pi_layer_addr  = NULL;
  int * pi_neuron_addr = NULL;
  int * pi_rate_addr   = NULL;
  int * pi_name_addr   = NULL;
  int * pi_momentum_addr  = NULL;
  int * pi_steepness_addr = NULL;
  int * pi_bit_addr       = NULL;
  int * pi_quickprop_addr = NULL;
  int * pi_rprop_addr     = NULL;
  int * pi_sarprop_addr   = NULL;
  int m_from_vect,   n_from_vect,   * pi_from_vect_addr   = NULL;
  int m_to_vect,     n_to_vect,     * pi_to_vect_addr     = NULL;
  int m_weight_vect, n_weight_vect, * pi_weight_vect_addr = NULL;
  int m_cascade,     n_cascade,     * pi_cascade_addr     = NULL;
  int res, i, j, managed = 0;
  int Index;
  char * Name = NULL, * Param = NULL;
  unsigned int actstep_count, numConnections;
  struct fann* ann   = NULL;
  char ** afunc_list = NULL;
  struct fann_connection * connections       = NULL;
  enum fann_activationfunc_enum * afunc_enum = NULL;
  SciErr _sciErr;

  if (Rhs < 2)
    {
      Scierror(999,"%s usage: value = %s(fann_in,'parameter_name')\n", fname, fname);
      return 0;
    }

  // Get the parameter name
  _sciErr = getVarAddressFromPosition(pvApiCtx, 2, &pi_name_addr);
  if (_sciErr.iErr)
    {
      printError(&_sciErr, 0);
      return 0;
    }

  getAllocatedSingleString(pvApiCtx,  pi_name_addr, &Param);

  // Get the Fann structure
  ann = createCFannStructFromScilabFannStruct(1,&res);
  if (res==-1) 
    {
      freeAllocatedSingleString(Param);
      return 0;
    }

  if (ann==NULL)
    {
      Scierror(999,"%s: Problem while creating the fann scilab structure\n",fname);
      freeAllocatedSingleString(Param);
      return 0;
    }

  // Verify that wa have a fannlist struct at pos 1
  res = detect_fannlist(1);
  if (res==-1)
    {
      freeAllocatedSingleString(Param);
      return 0;
    }

  // fann_set_training_algorithm                   Set the training algorithm.
  if (strcmp(Param,"training_algorithm")==0)
    {
      char * _tmp_mess = NULL;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', algorithm_name)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_name_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      
      getAllocatedSingleString(pvApiCtx,  pi_name_addr, &Name);

      _tmp_mess = (char *)MALLOC((strlen("%s: wrong training method - %s\n") + strlen(fname) + strlen(Name) + 10)*sizeof(char));
      sprintf(_tmp_mess,"%s: wrong training method - %s\n",fname, Name);
      FANN_SEARCH_NAME(FANN_TRAIN_NAMES,4,_tmp_mess);
      FREE(_tmp_mess);

      managed = 1;
      fann_set_training_algorithm(ann,(enum fann_train_enum)Index);
    }

  // fann_set_learning_rate                        Set the learning rate.
  if (strcmp(Param,"learning_rate")==0)
    {
      double rate = 0.0;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', learning_rate)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_rate_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_rate_addr, &rate);

      managed = 1;
      fann_set_learning_rate(ann,(float)rate);
    }

  // fann_set_learning_momentum                    Set the learning momentum.
  if (strcmp(Param,"learning_momentum")==0)
    {
      double momentum = 0.0;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', learning_momentum)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_momentum_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_momentum_addr, &momentum);

      managed = 1;
      fann_set_learning_momentum(ann,(float)momentum);
    }

  // fann_set_activation_function                  Set the activation function for neuron number neuron in layer number layer, counting the input layer as layer 0.
  if (strcmp(Param,"activation_function")==0)
    {
      char * _tmp_mess = NULL;
      double layer = 0.0, neuron = 0.0;

      if (Rhs!=5)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', index_layer, index_neuron, activation_function_name)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_layer_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_layer_addr, &layer);

      _sciErr = getVarAddressFromPosition(pvApiCtx, 4, &pi_neuron_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_neuron_addr, &neuron);

      _sciErr = getVarAddressFromPosition(pvApiCtx, 5, &pi_name_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getAllocatedSingleString(pvApiCtx, pi_name_addr, &Name);

      _tmp_mess = (char *)MALLOC((strlen("%s: wrong activation function - %s\n") + strlen(fname) + strlen(Name) + 10)*sizeof(char));
      sprintf(_tmp_mess,"%s: wrong activation function - %s\n",fname,Name);
      FANN_SEARCH_NAME(FANN_ACTIVATIONFUNC_NAMES,18,_tmp_mess)
      FREE(_tmp_mess);

      managed = 1;
      fann_set_activation_function(ann,(enum fann_activationfunc_enum)Index,layer,neuron);
      freeAllocatedSingleString(Name);
    }

  // fann_set_activation_function_layer            Set the activation function for all the neurons in the layer number layer, counting the input layer as layer 0.
  if (strcmp(Param,"activation_function_layer")==0)
    {
      char * _tmp_mess = NULL;
      double layer = 0.0;

      if (Rhs!=4)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', index_layer, function_name)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_layer_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_layer_addr, &layer);

      _sciErr = getVarAddressFromPosition(pvApiCtx, 4, &pi_name_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getAllocatedSingleString(pvApiCtx, pi_name_addr, &Name);

      _tmp_mess = (char *)MALLOC((strlen("%s: wrong activation function - %s\n") + strlen(fname) + strlen(Name) + 10)*sizeof(char));
      sprintf(_tmp_mess,"%s: wrong activation function - %s\n",fname,Name);
      FANN_SEARCH_NAME(FANN_ACTIVATIONFUNC_NAMES,18,_tmp_mess)
      FREE(_tmp_mess);

      managed = 1;
      fann_set_activation_function_layer(ann,(enum fann_activationfunc_enum)Index,layer);
      freeAllocatedSingleString(Name);
    }

  // fann_set_activation_function_hidden           Set the activation function for all of the hidden layers.
  if (strcmp(Param,"activation_function_hidden")==0)
    {
      char * _tmp_mess = NULL;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', function_name)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_name_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getAllocatedSingleString(pvApiCtx, pi_name_addr, &Name);

      _tmp_mess = (char *)MALLOC((strlen("%s: wrong activation function - %s\n") + strlen(fname) + strlen(Name) + 10)*sizeof(char));
      sprintf(_tmp_mess,"%s: wrong activation function - %s\n",fname, Name);
      FANN_SEARCH_NAME(FANN_ACTIVATIONFUNC_NAMES,18,_tmp_mess)
      FREE(_tmp_mess);

      managed = 1;
      fann_set_activation_function_hidden(ann,(enum fann_activationfunc_enum)Index);
      freeAllocatedSingleString(Name);
    }

  // fann_set_activation_function_output           Set the activation function for the output layer.
  if (strcmp(Param,"activation_function_output")==0)
    {
      char * _tmp_mess = NULL;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', function_name)\n", fname ,fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_name_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getAllocatedSingleString(pvApiCtx, pi_name_addr, &Name);

      _tmp_mess = (char *)MALLOC((strlen("%s: wrong activation function - %s\n") + strlen(fname) + strlen(Name) + 10)*sizeof(char));
      sprintf(_tmp_mess,"%s: wrong activation function\n",fname,Name);
      FANN_SEARCH_NAME(FANN_ACTIVATIONFUNC_NAMES,18,_tmp_mess)
      FREE(_tmp_mess);

      managed = 1;
      fann_set_activation_function_output(ann,(enum fann_activationfunc_enum)Index);
    }

  // fann_set_activation_steepness                 Set the activation steepness for neuron number neuron in layer number layer, counting the input layer as layer 0.
  if (strcmp(Param,"activation_steepness")==0)
    {
      double layer = 0.0, neuron = 0.0, steepness = 0.0;

      if (Rhs!=5)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', index_layer, index_neuron, steepness)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_layer_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_layer_addr, &layer);

      _sciErr = getVarAddressFromPosition(pvApiCtx, 4, &pi_neuron_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_neuron_addr, &neuron);

      _sciErr = getVarAddressFromPosition(pvApiCtx, 5, &pi_steepness_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_steepness_addr, &steepness);

      managed = 1;
      fann_set_activation_steepness(ann,(fann_type)steepness,layer,neuron);
    }

  // fann_set_activation_steepness_layer           Set the activation steepness all of the neurons in layer number layer, counting the input layer as layer 0.
  if (strcmp(Param,"activation_steepness_layer")==0)
    {
      double layer = 0.0, steepness = 0.0;

      if (Rhs!=4)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', index_layer, steepness)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_layer_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_layer_addr, &layer);

      _sciErr = getVarAddressFromPosition(pvApiCtx, 4, &pi_steepness_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_steepness_addr, &steepness);

      managed = 1;
      fann_set_activation_steepness_layer(ann,(fann_type)steepness,layer);
    }

  // fann_set_activation_steepness_hidden          Set the steepness of the activation steepness in all of the hidden layers.
  if (strcmp(Param,"activation_steepness_hidden")==0)
    {
      double steepness = 0.0;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', steepness)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_steepness_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_steepness_addr, &steepness);

      managed = 1;
      fann_set_activation_steepness_hidden(ann,(fann_type)steepness);
    }

  // fann_set_activation_steepness_output          Set the steepness of the activation steepness in the output layer.
  if (strcmp(Param,"activation_steepness_output")==0)
    {
      double steepness = 0.0;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', steepness)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_steepness_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_steepness_addr, &steepness);

      managed = 1;
      fann_set_activation_steepness_output(ann,(fann_type)steepness);
    }

  // fann_set_train_error_function                 Set the error function used during training.
  if (strcmp(Param,"train_error_function")==0)
    {
      char * _tmp_mess = NULL;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', function_name)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_name_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      
      getAllocatedSingleString(pvApiCtx,  pi_name_addr, &Name);

      _tmp_mess = (char *)MALLOC((strlen("%s: wrong training error function - %s\n") + strlen(fname) + strlen(Name) + 10)*sizeof(char));
      sprintf(_tmp_mess,"%s: wrong training error function\n",fname,Name);
      FANN_SEARCH_NAME(FANN_ERRORFUNC_NAMES,2,_tmp_mess)
      FREE(_tmp_mess);

      managed = 1;
      fann_set_train_error_function(ann,(enum fann_activationfunc_enum)Index);
      freeAllocatedSingleString(Name);
    }

  // fann_set_train_stop_function                  Set the stop function used during training.
  if (strcmp(Param,"train_stop_function")==0)
    {
      char * _tmp_mess = NULL;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', function_name)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_name_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      
      getAllocatedSingleString(pvApiCtx,  pi_name_addr, &Name);

      _tmp_mess = (char *)MALLOC((strlen("%s: wrong stop training function - %s\n") + strlen(fname) + strlen(Name) + 10)*sizeof(char));
      sprintf(_tmp_mess,"%s: wrong stop training function - %s\n",fname,Name);
      FANN_SEARCH_NAME(FANN_STOPFUNC_NAMES,2,_tmp_mess)
      FREE(_tmp_mess);

      managed = 1;
      fann_set_train_error_function(ann,(enum fann_stopfunc_enum)Index);
      freeAllocatedSingleString(Name);
    }

  // fann_set_bit_fail_limit                       Set the bit fail limit used during training.
  if (strcmp(Param,"bit_fail_limit")==0)
    {
      double bit = 0.0;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', bit_fail_limit)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_bit_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_bit_addr, &bit);

      managed = 1;
      fann_set_bit_fail_limit(ann,(fann_type)bit);
    }

  // fann_set_quickprop_decay                      Sets the quickprop decay factor.
  if (strcmp(Param,"quickprop_decay")==0)
    {
      double quickprop = 0.0;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', quickprop_decay)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_quickprop_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_quickprop_addr, &quickprop);

      managed = 1;
      fann_set_quickprop_decay(ann,(float)quickprop);
    }

  // fann_set_quickprop_mu                         Sets the quickprop mu factor.
  if (strcmp(Param,"quickprop_mu")==0)
    {
      double quickprop = 0.0;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', quickprop_mu)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_quickprop_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_quickprop_addr, &quickprop);

      managed = 1;
      fann_set_quickprop_mu(ann,(float)quickprop);
    }

  // fann_set_rprop_increase_factor                The increase factor used during RPROP training.
  if (strcmp(Param,"rprop_increase_factor")==0)
    {
      double rprop = 0.0;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', rprop_increase_factor)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_rprop_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_rprop_addr, &rprop);

      managed = 1;
      fann_set_rprop_increase_factor(ann,(float)rprop);
    }

  // fann_set_rprop_decrease_factor                The decrease factor is a value smaller than 1, which is used to decrease the step-size during RPROP training.
  if (strcmp(Param,"rprop_decrease_factor")==0)
    {
      double rprop = 0.0;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', rprop_decrease_factor)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_rprop_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_rprop_addr, &rprop);

      managed = 1;
      fann_set_rprop_decrease_factor(ann,(float)rprop);
    }

  // fann_set_rprop_delta_min                      The minimum step-size is a small positive number determining how small the minimum step-size may be.
  if (strcmp(Param,"rprop_delta_min")==0)
    {
      double rprop = 0.0;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', rprop_delta_min)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_rprop_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_rprop_addr, &rprop);

      managed = 1;
      fann_set_rprop_delta_min(ann,(float)rprop);
    }

  // fann_set_rprop_delta_max                      The maximum step-size is a positive number determining how large the maximum step-size may be.
  if (strcmp(Param,"rprop_delta_max")==0)
    {
      double rprop = 0.0;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', rprop_delta_max)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_rprop_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_rprop_addr, &rprop);

      managed = 1;
      fann_set_rprop_delta_max(ann,(float)rprop);
    }

  // fann_set_rprop_delta_zero                     The initial step-size is a positive number determining the initial step size.
  if (strcmp(Param,"rprop_delta_zero")==0)
    {
      double rprop = 0.0;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', rprop_delta_zero)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_rprop_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_rprop_addr, &rprop);

      managed = 1;
      fann_set_rprop_delta_zero(ann,(float)rprop);
    }

  // fann_set_sarprop_weight_decay_shift           Set the sarprop weight decay shift.
  if (strcmp(Param,"sarprop_weight_decay_shift")==0)
    {
      double sarprop = 0.0;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', sarprop_weight_decay_shift)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_sarprop_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_sarprop_addr, &sarprop);

      managed = 1;
      fann_set_sarprop_weight_decay_shift(ann,(float)sarprop);
    }

  // fann_set_sarprop_step_error_threshold_factor  Set the sarprop step error threshold factor.
  if (strcmp(Param,"sarprop_step_error_threshold_factor")==0)
    {
      double sarprop = 0.0;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', sarprop_step_error_threshold_factor)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_sarprop_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_sarprop_addr, &sarprop);

      managed = 1;
      fann_set_sarprop_step_error_threshold_factor(ann,(float)sarprop);
    }

  // fann_set_sarprop_step_error_shift             Set the sarprop step error shift.
  if (strcmp(Param,"sarprop_step_error_shift")==0)
    {
      double sarprop = 0.0;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', sarprop_step_error_shift)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_sarprop_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_sarprop_addr, &sarprop);

      managed = 1;
      fann_set_sarprop_step_error_shift(ann,(float)sarprop);
    }

  // fann_set_sarprop_temperature                  Set the sarprop_temperature.
  if (strcmp(Param,"sarprop_temperature")==0)
    {
      double sarprop = 0.0;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', sarprop_temperature)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_sarprop_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_sarprop_addr, &sarprop);

      managed = 1;
      fann_set_sarprop_temperature(ann,(float)sarprop);
    }
	
  // fann_set_weight_array Set connections in the network.
  if (strcmp(Param,"weight_array")==0)
    {
      double * from_vect = NULL, * to_vect = NULL, * weight_vect = NULL;

      if (Rhs!=5)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', from_vect, to_vect, weight_vect)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_from_vect_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      _sciErr = getMatrixOfDouble(pvApiCtx, pi_from_vect_addr, &m_from_vect, &n_from_vect, &from_vect);

      _sciErr = getVarAddressFromPosition(pvApiCtx, 4, &pi_to_vect_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      _sciErr = getMatrixOfDouble(pvApiCtx, pi_to_vect_addr, &m_to_vect, &n_to_vect, &to_vect);

      _sciErr = getVarAddressFromPosition(pvApiCtx, 5, &pi_weight_vect_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      _sciErr = getMatrixOfDouble(pvApiCtx, pi_weight_vect_addr, &m_weight_vect, &n_weight_vect, &weight_vect);

      if ((m_from_vect!=m_to_vect) | (m_from_vect!=m_weight_vect) |
	  (n_from_vect!=n_to_vect) | (n_from_vect!=n_weight_vect))
	{
	  Scierror(999,"%s: weight_array - wrong size for arguments\n",fname);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      numConnections = m_from_vect * n_from_vect;

      connections = (struct fann_connection *)MALLOC(numConnections*sizeof(struct fann_connection));
      for(i=0; i<(int)numConnections;i++)
	{
	  connections[i].from_neuron = from_vect[i];
	  connections[i].to_neuron   = to_vect[i];
	  connections[i].weight      = (fann_type)weight_vect[i];
	}

      managed = 1;
      fann_set_weight_array(ann,connections,numConnections);

      FREE(connections);
    }

  // fann_set_weight       Set a connection in the network.
  if (strcmp(Param,"weight")==0)
    {
      double from_vect = 0.0, to_vect = 0.0;
      double weight_vect = 0.0;

      if (Rhs!=5)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', from, to, weight)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_from_vect_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      getScalarDouble(pvApiCtx, pi_from_vect_addr, &from_vect);

      _sciErr = getVarAddressFromPosition(pvApiCtx, 4, &pi_to_vect_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      getScalarDouble(pvApiCtx, pi_to_vect_addr, &to_vect);

      _sciErr = getVarAddressFromPosition(pvApiCtx, 5, &pi_weight_vect_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}
      getScalarDouble(pvApiCtx, pi_weight_vect_addr, &weight_vect);

      managed = 1;
      fann_set_weight(ann,(unsigned int)from_vect,(unsigned int)to_vect,(fann_type)weight_vect);
    }

  // fann_set_cascade_output_change_fraction       Sets the cascade output change fraction.
  if (strcmp(Param,"cascade_output_change_fraction")==0)
    {
      double cascade = 0.0;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', cascade_output_change_fraction)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_cascade_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_cascade_addr, &cascade);

      managed = 1;
      fann_set_cascade_output_change_fraction(ann,(float)cascade);
    }

  // fann_set_cascade_output_stagnation_epochs     Sets the number of cascade output stagnation epochs.
  if (strcmp(Param,"cascade_output_stagnation_epochs")==0)
    {
      double cascade = 0.0;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', cascade_output_stagnation_epochs)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_cascade_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_cascade_addr, &cascade);

      managed = 1;
      fann_set_cascade_output_stagnation_epochs(ann,(unsigned int)cascade);
    }

  // fann_set_cascade_candidate_change_fraction    Sets the cascade candidate change fraction.
  if (strcmp(Param,"cascade_candidate_change_fraction")==0)
    {
      double cascade = 0.0;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', cascade_candidate_change_fraction)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_cascade_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_cascade_addr, &cascade);

      managed = 1;
      fann_set_cascade_candidate_change_fraction(ann,(float)cascade);
    }

  // fann_set_cascade_candidate_stagnation_epochs  Sets the number of cascade candidate stagnation epochs.
  if (strcmp(Param,"cascade_candidate_stagnation_epochs")==0)
    {
      double cascade = 0.0;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', cascade_candidate_stagnation_epochs)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_cascade_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_cascade_addr, &cascade);

      managed = 1;
      fann_set_cascade_candidate_stagnation_epochs(ann,(float)cascade);
    }

  // fann_set_cascade_weight_multiplier            Sets the weight multiplier.
  if (strcmp(Param,"cascade_weight_multiplier")==0)
    {
      double cascade = 0.0;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', cascade_weight_multiplier)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_cascade_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_cascade_addr, &cascade);

      managed = 1;
      fann_set_cascade_weight_multiplier(ann,(fann_type)cascade);
    }

  // fann_set_cascade_candidate_limit              Sets the candidate limit.
  if (strcmp(Param,"cascade_candidate_limit")==0)
    {
      double cascade = 0.0;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', cascade_candidate_limit)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_cascade_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_cascade_addr, &cascade);

      managed = 1;
      fann_set_cascade_candidate_limit(ann,(fann_type)cascade);
    }

  // fann_set_cascade_max_out_epochs               Sets the maximum out epochs.
  if (strcmp(Param,"cascade_max_out_epochs")==0)
    {
      double cascade = 0.0;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', cascade_max_out_epochs)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_cascade_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_cascade_addr, &cascade);

      managed = 1;
      fann_set_cascade_max_out_epochs(ann,(unsigned int)cascade);
    }

  // fann_set_cascade_min_out_epochs               Sets the minimum out epochs.
  if (strcmp(Param,"cascade_min_out_epochs")==0)
    {
      double cascade = 0.0;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', cascade_min_out_epochs)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_cascade_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_cascade_addr, &cascade);

      managed = 1;
      fann_set_cascade_min_out_epochs(ann,(unsigned int)cascade);
    }

  // fann_set_cascade_max_cand_epochs              Sets the max candidate epochs.
  if (strcmp(Param,"cascade_max_cand_epochs")==0)
    {
      double cascade = 0.0;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', cascade_max_cand_epochs)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_cascade_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_cascade_addr, &cascade);

      managed = 1;
      fann_set_cascade_max_cand_epochs(ann,(unsigned int)cascade);
    }

  // fann_set_cascade_min_cand_epochs              Sets the min candidate epochs.
  if (strcmp(Param,"cascade_min_cand_epochs")==0)
    {
      double cascade = 0.0;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', cascade_min_cand_epochs)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_cascade_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_cascade_addr, &cascade);

      managed = 1;
      fann_set_cascade_min_cand_epochs(ann,(unsigned int)cascade);
    }

  // fann_set_cascade_activation_functions         Sets the array of cascade candidate activation functions.
  if (strcmp(Param,"cascade_activation_functions")==0)
    {
      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', function_names)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_name_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getAllocatedMatrixOfString(pvApiCtx, pi_name_addr, &m_cascade, &n_cascade, &afunc_list);

      afunc_enum = (enum fann_activationfunc_enum *)MALLOC(m_cascade*sizeof(enum fann_activationfunc_enum));

      for(i=0; i<m_cascade; i++)
	{
	  Index = -1;
	  for(j=0; j<18; j++)
	    {
	      if (strcmp(afunc_list[i],FANN_ACTIVATIONFUNC_NAMES[j])==0)
		{
		  Index = j;
		  break;
		}
	    }
	  if (Index==-1)
	    {
	      Scierror(999,"%s: cascade_activation_functions - wrong activation function name: %s\n", fname, Param);
	      freeAllocatedSingleString(Param);
	      freeAllocatedMatrixOfString(m_cascade, n_cascade, afunc_list);
	      return 0;
	    }
	  else
	    {
	      afunc_enum[i] = Index;
	    }
	}

      managed = 1;
      fann_set_cascade_activation_functions(ann,afunc_enum,m_cascade);
      freeAllocatedMatrixOfString(m_cascade, n_cascade, afunc_list);
    }

  // fann_set_cascade_activation_steepnesses       Sets the array of cascade candidate activation steepnesses.
  if (strcmp(Param,"cascade_activation_steepnesses")==0)
    {
      double * cascade = NULL;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', cascade_activation_steepnesses)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_cascade_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      _sciErr = getMatrixOfDouble(pvApiCtx, pi_cascade_addr, &m_cascade, &n_cascade, &cascade);

      actstep_count = m_cascade * n_cascade;
      
      managed = 1;
      fann_set_cascade_activation_steepnesses(ann,(fann_type *)cascade,actstep_count);
    }

  // fann_set_cascade_num_candidate_groups         Sets the number of candidate groups.
  if (strcmp(Param,"cascade_num_candidate_groups")==0)
    {
      double cascade = 0.0;

      if (Rhs!=3)
	{
	  Scierror(999,"%s: usage %s(fann_in, '%s', cascade_num_candidate_groups)\n", fname, fname, Param);
	  freeAllocatedSingleString(Param);
	  return 0;
	}

      _sciErr = getVarAddressFromPosition(pvApiCtx, 3, &pi_cascade_addr);
      if (_sciErr.iErr)
	{
	  printError(&_sciErr, 0);
	  return 0;
	}

      getScalarDouble(pvApiCtx, pi_cascade_addr, &cascade);

      managed = 1;
      fann_set_cascade_num_candidate_groups(ann,(unsigned int)cascade);
    }

  if (!managed)
    {
      Scierror(999,"%s: error, wrong command: %s\n", fname, Param);
      freeAllocatedSingleString(Param);
      return 0;
    }
  else
    {
      //Create the struct representing this ann in Scilab
      res = createScilabFannStructFromCFannStruct(ann, Rhs + 1);
      freeAllocatedSingleString(Param);
      if (res==-1) return 0;
      
      LhsVar(1) = Rhs + 1;
    }

  return 0;
}
