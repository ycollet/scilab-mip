//
//  Fast Artificial Neural Network Library (fann)
//  Copyright (C) 2003-2012 Steffen Nissen (sn@leenissen.dk)
//  
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//  
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//  
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//

desired_error = 0.0;
max_neurons   = 30;
neurons_between_reports = 1;
multi = %F;
training_algorithm = 'FANN_TRAIN_RPROP';

printf('Reading data.\n');

traindata = fann_readtrain_from_file("data/parity8.train");
testdata = fann_readtrain_from_file("data/parity8.test");

traindata = fann_setup_train_data(traindata, 'scale_train_data', -1, 1);
testdata  = fann_setup_train_data(testdata,  'scale_train_data', -1, 1);

printf('Creating network.\n');

ann = fann_create('shortcut', [fann_setup_train_data(traindata, 'num_input_train_data'), fann_setup_train_data(traindata, 'num_output_train_data')]);
		
ann = fann_set_parameters(ann, 'training_algorithm', training_algorithm);
ann = fann_set_parameters(ann, 'activation_function_hidden', 'FANN_SIGMOID_SYMMETRIC');
ann = fann_set_parameters(ann, 'activation_function_output', 'FANN_LINEAR');
ann = fann_set_parameters(ann, 'train_error_function', 'FANN_ERRORFUNC_LINEAR');

if (~multi) then
  steepness = 0.5;
  //steepness = 1;
  ann = fann_set_parameters(ann, 'cascade_activation_steepnesses', steepness);
  activation = 'FANN_SIN_SYMMETRIC';
  //activation = 'FANN_SIGMOID_SYMMETRIC';
  
  ann = fann_set_parameters(ann, 'cascade_activation_functions', activation);
  ann = fann_set_parameters(ann, 'cascade_num_candidate_groups', 8);
end	
		
if (training_algorithm=='FANN_TRAIN_QUICKPROP') then
  ann = fann_set_parameters(ann, 'learning_rate', 0.35);
  ann = fann_set_parameters(ann, 'randomize_weights', -2.0, 2.0);
end
	
ann = fann_set_parameters(ann, 'bit_fail_limit', 0.9);
ann = fann_set_parameters(ann, 'train_stop_function', 'FANN_STOPFUNC_BIT');
// fann_print_parameters(ann);
		
fann_save(ann, 'cascade_train2.net');
	
printf('Training network.\n');

fann_casc_train_on_data(ann, traindata, max_neurons, desired_error);
	
fann_print_connections(ann);
	
ann = fann_test_data(ann, traindata);
mse_train = fann_get_MSE(ann);
bit_fail_train = fann_get_bit_fail(ann);
ann = fann_test_data(ann, testdata);
mse_test = fann_get_MSE(ann);
bit_fail_test = fann_get_bit_fail(ann);

printf('\nTrain error: %f, Train bit-fail: %d, Test error: %f, Test bit-fail: %d\n\n', mse_train, bit_fail_train, mse_test, bit_fail_test);

for i=1:fann_setup_train_data(traindata, 'length_train_data')
  data_in  = fann_get_input(traindata, i);
  data_out = fann_get_output(traindata, i);
  output = fann_run(ann, data_in);
  if ((data_out(1) >= 0 & output(1) <= 0) | (data_out(1) <= 0 & output(1) >= 0)) then
    printf("ERROR: %f does not match %f\n", data_out(1), output(1));
  end
end
	
printf('Saving network.\n');

fann_save(ann, 'cascade_train.net');

printf('Cleaning up.\n');

fann_destroy_train(traindata); clear traindata;
fann_destroy_train(testdata); clear testdata;
fann_destroy(ann); clear ann;
