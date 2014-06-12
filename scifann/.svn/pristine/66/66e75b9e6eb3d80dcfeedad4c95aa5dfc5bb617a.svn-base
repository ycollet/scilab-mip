lines(0);

path = get_absolute_file_path('unit_test.sce');

NbPointsLearn = 1000;
NbPointsValid = 100;
desired_error = 1e-6;
max_epochs    = 5000;

Xmax = [2 2];
Xmin = [-2 -2];

// Generate the learning data set
for i=1:NbPointsLearn
  xData_learn(i,:) = (Xmax - Xmin).*rand(size(Xmax,1),size(Xmax,2)) + Xmin;
  fData_learn(i) = sum(xData_learn(i,:).^2);
end
// Normalisation of the learning data set
x_min_norm = min(xData_learn,'r');
x_max_norm = max(xData_learn,'r');
f_min_norm = min(fData_learn,'r');
f_max_norm = max(fData_learn,'r');
// Now we normalize so as the data spread between -1 and 1.
for i=1:NbPointsLearn
  xData_learn(i,:) = 2 * (xData_learn(i,:) - x_min_norm) ./ (x_max_norm - x_min_norm) - 1;
  fData_learn(i)   = 2 * (fData_learn(i)   - f_min_norm)  / (f_max_norm - f_min_norm) - 1;
end

// Generate the validation data set
for i=1:NbPointsValid
  xData_valid(i,:) = (Xmax - Xmin).*rand(size(Xmax,1),size(Xmax,2)) + Xmin;
  fData_valid(i) = sum(xData_valid(i,:).^2);
end
// Now we also normalize the validation data set using the normalization parameters from
// the learning so as the data spread between -1 and 1.
for i=1:NbPointsValid
  xData_valid(i,:) = 2 * (xData_valid(i,:) - x_min_norm) ./ (x_max_norm - x_min_norm) - 1;
  fData_valid(i)   = 2 * (fData_valid(i)   - f_min_norm)  / (f_max_norm - f_min_norm) - 1;
end

printf('call to createfann\n');
ann = fann_create('sparse',[2 2 1], 0.8);

// get the default parameters
printf('num_input                           : %d\n',fann_get_parameters(ann,'num_input'));
printf('num_output                          : %d\n',fann_get_parameters(ann,'num_output'));
printf('total_neurons                       : %d\n',fann_get_parameters(ann,'total_neurons'));
printf('total_connections                   : %d\n',fann_get_parameters(ann,'total_connections'));
printf('network_type                        : %s\n',fann_get_parameters(ann,'network_type'));
printf('connection_rate                     : %f\n',fann_get_parameters(ann,'connection_rate'));
printf('num_layers                          : %d\n',fann_get_parameters(ann,'num_layers'));
printf('layers_array                        :'); disp(fann_get_parameters(ann,'layer_array'));
printf('bias_array                          :'); disp(fann_get_parameters(ann,'bias_array'));
printf('connection_array                    :'); disp(fann_get_parameters(ann,'connection_array'));
printf('training algorithm                  : %s\n',fann_get_parameters(ann,'training_algorithm'));
printf('learning_rate                       : %f\n',fann_get_parameters(ann,'learning_rate'));
printf('learning_momentum                   : %f\n',fann_get_parameters(ann,'learning_momentum'));
printf('activation_function                 : %s\n',fann_get_parameters(ann,'activation_function',1,1));
printf('activation_steepness                : %f\n',fann_get_parameters(ann,'activation_steepness',1,1));
printf('train_error_function                : %s\n',fann_get_parameters(ann,'train_error_function'));
printf('train_stop_function                 : %s\n',fann_get_parameters(ann,'train_stop_function'));
printf('bit_fail_limit                      : %f\n',fann_get_parameters(ann,'bit_fail_limit'));
printf('quickprop_decay                     : %f\n',fann_get_parameters(ann,'quickprop_decay'));
printf('rprop_increase_factor               : %f\n',fann_get_parameters(ann,'rprop_increase_factor'));
printf('quickprop_mu                        : %f\n',fann_get_parameters(ann,'quickprop_mu'));
printf('rprop_decrease_factor               : %f\n',fann_get_parameters(ann,'rprop_decrease_factor'));
printf('rprop_delta_min                     : %f\n',fann_get_parameters(ann,'rprop_delta_min'));
printf('rprop_delta_max                     : %f\n',fann_get_parameters(ann,'rprop_delta_max'));
printf('rprop_delta_zero                    : %f\n',fann_get_parameters(ann,'rprop_delta_zero'));
printf('cascade_output_change_fraction      : %f\n',fann_get_parameters(ann,'cascade_output_change_fraction'));
printf('cascade_output_stagnation_epochs    : %d\n',fann_get_parameters(ann,'cascade_output_stagnation_epochs'));
printf('cascade_candidate_change_fraction   : %f\n',fann_get_parameters(ann,'cascade_candidate_change_fraction'));
printf('cascade_candidate_stagnation_epochs : %d\n',fann_get_parameters(ann,'cascade_candidate_stagnation_epochs'));
printf('cascade_weight_multiplier           : %f\n',fann_get_parameters(ann,'cascade_weight_multiplier'));
printf('cascade_candidate_limit             : %f\n',fann_get_parameters(ann,'cascade_candidate_limit'));
printf('cascade_max_out_epochs              : %d\n',fann_get_parameters(ann,'cascade_max_out_epochs'));
printf('cascade_max_cand_epochs             : %d\n',fann_get_parameters(ann,'cascade_max_cand_epochs'));
printf('cascade_num_candidates              : %d\n',fann_get_parameters(ann,'cascade_num_candidates'));
printf('cascade_activation_functions_count  : %d\n',fann_get_parameters(ann,'cascade_activation_functions_count'));
printf('cascade_activation_functions        : '); disp(fann_get_parameters(ann,'cascade_activation_functions'));
printf('cascade_activation_steepnesses_count: %d\n',fann_get_parameters(ann,'cascade_activation_steepnesses_count'));
printf('cascade_activation_steepnesses      : '); disp(fann_get_parameters(ann,'cascade_activation_steepnesses'));
printf('cascade_num_candidate_groups        : %d\n',fann_get_parameters(ann,'cascade_num_candidate_groups'));


// set the parameters and display the values
//
// The available learning methods:
//
//	'FANN_TRAIN_INCREMENTAL'
//	'FANN_TRAIN_BATCH'
//	'FANN_TRAIN_RPROP'
//	'FANN_TRAIN_QUICKPROP'
//
// The available activation functions:
//
//	'FANN_LINEAR'
//	'FANN_THRESHOLD'
//	'FANN_THRESHOLD_SYMMETRIC'
//	'FANN_SIGMOID'
//	'FANN_SIGMOID_STEPWISE'
//	'FANN_SIGMOID_SYMMETRIC'
//	'FANN_SIGMOID_SYMMETRIC_STEPWISE'
//	'FANN_GAUSSIAN'
//	'FANN_GAUSSIAN_SYMMETRIC'
//	'FANN_GAUSSIAN_STEPWISE'
//	'FANN_ELLIOT'
//	'FANN_ELLIOT_SYMMETRIC'
//	'FANN_LINEAR_PIECE'
//	'FANN_LINEAR_PIECE_SYMMETRIC'
//	'FANN_SIN_SYMMETRIC'
//	'FANN_COS_SYMMETRIC'
//	'FANN_SIN'
//	'FANN_COS'
//
// The available error functions:
//
//	'FANN_ERRORFUNC_LINEAR'
//	'FANN_ERRORFUNC_TANH'
//
// The available stopping functions:
//
//	'FANN_STOPFUNC_MSE'
//	'FANN_STOPFUNC_BIT'
//
// The available network types:
//
//	'FANN_NETTYPE_LAYER'
//	'FANN_NETTYPE_SHORTCUT'

printf('num_input                           : %d\n',fann_get_parameters(ann,'num_input'));
printf('num_output                          : %d\n',fann_get_parameters(ann,'num_output'));
printf('total_neurons                       : %d\n',fann_get_parameters(ann,'total_neurons'));
printf('total_connections                   : %d\n',fann_get_parameters(ann,'total_connections'));
printf('network_type                        : %s\n',fann_get_parameters(ann,'network_type'));
printf('connection_rate                     : %f\n',fann_get_parameters(ann,'connection_rate'));
printf('num_layers                          : %d\n',fann_get_parameters(ann,'num_layers'));
printf('layer_array                         :'); disp(fann_get_parameters(ann,'layer_array'));
printf('bias_array                          :'); disp(fann_get_parameters(ann,'bias_array'));
printf('connection_array                    :'); disp(fann_get_parameters(ann,'connection_array'));
fann_set_parameters(ann,'training_algorithm','FANN_TRAIN_INCREMENTAL');
printf('training algorithm                  : %s\n',fann_get_parameters(ann,'training_algorithm'));
fann_set_parameters(ann,'learning_rate',0.9);
printf('learning_rate                       : %f\n',fann_get_parameters(ann,'learning_rate'));
fann_set_parameters(ann,'learning_momentum',0.9);
printf('learning_momentum                   : %f\n',fann_get_parameters(ann,'learning_momentum'));
fann_set_parameters(ann,'activation_function',1,1,'FANN_SIGMOID');
printf('activation_function                 : %s\n',fann_get_parameters(ann,'activation_function',1,1));
fann_set_parameters(ann,'activation_function_layer',1,'FANN_SIGMOID');
printf('activation_function                 : %s\n',fann_get_parameters(ann,'activation_function',1,1));
fann_set_parameters(ann,'activation_function_hidden','FANN_SIGMOID');
printf('activation_function                 : %s\n',fann_get_parameters(ann,'activation_function',1,1));
fann_set_parameters(ann,'activation_function_output','FANN_SIGMOID');
printf('activation_function                 : %s\n',fann_get_parameters(ann,'activation_function',1,1));
fann_set_parameters(ann,'activation_steepness',1,1,0.1);
printf('activation_steepness                : %f\n',fann_get_parameters(ann,'activation_steepness',1,1));
fann_set_parameters(ann,'activation_steepness_layer',1,0.1);
printf('activation_steepness                : %f\n',fann_get_parameters(ann,'activation_steepness',1,1));
fann_set_parameters(ann,'activation_steepness_hidden',0.1);
printf('activation_steepness                : %f\n',fann_get_parameters(ann,'activation_steepness',1,1));
fann_set_parameters(ann,'activation_steepness_output',0.1);
printf('activation_steepness                : %f\n',fann_get_parameters(ann,'activation_steepness',1,1));
fann_set_parameters(ann,'train_error_function','FANN_ERRORFUNC_LINEAR');
printf('train_error_function                : %s\n',fann_get_parameters(ann,'train_error_function'));
fann_set_parameters(ann,'train_stop_function','FANN_STOPFUNC_MSE');
printf('train_stop_function                 : %s\n',fann_get_parameters(ann,'train_stop_function'));
fann_set_parameters(ann,'bit_fail_limit',0.1);
printf('bit_fail_limit                      : %f\n',fann_get_parameters(ann,'bit_fail_limit'));
fann_set_parameters(ann,'quickprop_decay',0.1);
printf('quickprop_decay                     : %f\n',fann_get_parameters(ann,'quickprop_decay'));
fann_set_parameters(ann,'quickprop_mu',0.1);
printf('quickprop_mu                        : %f\n',fann_get_parameters(ann,'quickprop_mu'));
fann_set_parameters(ann,'rprop_increase_factor',1.1);
printf('rprop_increase_factor               : %f\n',fann_get_parameters(ann,'rprop_increase_factor'));
fann_set_parameters(ann,'rprop_decrease_factor',0.8);
printf('rprop_decrease_factor               : %f\n',fann_get_parameters(ann,'rprop_decrease_factor'));
fann_set_parameters(ann,'rprop_delta_min',0.2);
printf('rprop_delta_min                     : %f\n',fann_get_parameters(ann,'rprop_delta_min'));
fann_set_parameters(ann,'rprop_delta_max',0.3);
printf('rprop_delta_max                     : %f\n',fann_get_parameters(ann,'rprop_delta_max'));
fann_set_parameters(ann,'rprop_delta_zero',0.1);
printf('rprop_delta_zero                    : %f\n',fann_get_parameters(ann,'rprop_delta_zero'));

from   = [1 2 3];
to     = [2 3 1];
weight = [0.1 0.2 0.3];
fann_set_parameters(ann,'weight_array',from,to,weight);

fann_set_parameters(ann,'weight_array',1,2,0.5);

fann_set_parameters(ann,'cascade_output_change_fraction',0.1);
printf('cascade_output_change_fraction      : %f\n',fann_get_parameters(ann,'cascade_output_change_fraction'));
fann_set_parameters(ann,'cascade_output_stagnation_epochs',10);
printf('cascade_output_stagnation_epochs    : %d\n',fann_get_parameters(ann,'cascade_output_stagnation_epochs'));
fann_set_parameters(ann,'cascade_candidate_change_fraction',0.1);
printf('cascade_candidate_change_fraction   : %f\n',fann_get_parameters(ann,'cascade_candidate_change_fraction'));
fann_set_parameters(ann,'cascade_candidate_stagnation_epochs',10);
printf('cascade_candidate_stagnation_epochs : %d\n',fann_get_parameters(ann,'cascade_candidate_stagnation_epochs'));
fann_set_parameters(ann,'cascade_weight_multiplier',1.1);
printf('cascade_weight_multiplier           : %f\n',fann_get_parameters(ann,'cascade_weight_multiplier'));
fann_set_parameters(ann,'cascade_candidate_limit',10);
printf('cascade_candidate_limit             : %f\n',fann_get_parameters(ann,'cascade_candidate_limit'));
fann_set_parameters(ann,'cascade_max_out_epochs',10);
printf('cascade_max_out_epochs              : %d\n',fann_get_parameters(ann,'cascade_max_out_epochs'));
fann_set_parameters(ann,'cascade_max_cand_epochs',10);
printf('cascade_max_cand_epochs             : %d\n',fann_get_parameters(ann,'cascade_max_cand_epochs'));
fann_set_parameters(ann,'cascade_num_candidate_groups',3);
printf('cascade_num_candidates              : %d\n',fann_get_parameters(ann,'cascade_num_candidate_groups'));
fann_set_parameters(ann,'cascade_activation_functions',['FANN_SIGMOID','FANN_SIGMOID_STEPWISE','FANN_SIGMOID_SYMMETRIC']);
printf('cascade_activation_functions        : '); disp(fann_get_parameters(ann,'cascade_activation_functions'));
printf('cascade_activation_functions_count  : %d\n',fann_get_parameters(ann,'cascade_activation_functions_count'));
fann_set_parameters(ann,'cascade_activation_steepnesses',[0.1,0.1,0.2]);
printf('cascade_activation_steepnesses_count: %d\n',fann_get_parameters(ann,'cascade_activation_steepnesses_count'));
printf('cascade_activation_steepnesses      : '); disp(fann_get_parameters(ann,'cascade_activation_steepnesses'));

// error functions
//ann_error = fann_set_error_log('test_error.log');
ann_error = fann_set_error_log();

printf('fann_get_errno = %d\n',fann_get_errno(ann_error));
printf('fann_get_errstr = %s\n',fann_get_errstr(ann_error));
fann_reset_errno(ann_error);
fann_reset_errstr(ann_error);

// save / read functions

fann_save(ann,path + 'test_save.fann');
fann_savetofixed(ann,path + 'test_savetofixed.fann');
new_ann = fann_read(path + 'test_save.fann');
