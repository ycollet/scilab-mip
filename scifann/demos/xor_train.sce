// 
// Fast Artificial Neural Network Library (fann)
// Copyright (C) 2003-2012 Steffen Nissen (sn@leenissen.dk)
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// 
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// 

num_input          = 2;
num_output         = 1;
num_neurons_hidden = 3;
desired_error      = 0;
max_epochs         = 1000;

printf('Creating network.\n');

ann = fann_create('standard', [num_input, num_neurons_hidden, num_output]);

data = fann_readtrain_from_file('data/xor.data');

ann = fann_set_parameters(ann, 'activation_steepness_hidden', 1);
ann = fann_set_parameters(ann, 'activation_steepness_output', 1);

ann = fann_set_parameters(ann, 'activation_function_hidden', 'FANN_SIGMOID_SYMMETRIC');
ann = fann_set_parameters(ann, 'activation_function_output', 'FANN_SIGMOID_SYMMETRIC');

ann = fann_set_parameters(ann, 'train_stop_function', 'FANN_STOPFUNC_BIT');
ann = fann_set_parameters(ann, 'bit_fail_limit', 0.01);

ann = fann_set_parameters(ann, 'training_algorithm', 'FANN_TRAIN_RPROP');

ann = fann_init_weights(ann, data);
	
printf('Training network.\n');

ann = fann_train_on_data(ann, data, max_epochs, desired_error, epochs_between_reports);

ann = fann_test_data(ann, data);

printf('Testing network. %f\n', fann_get_MSE(ann));

for i=1:fann_setup_train_data(data, 'length_train_data')
  data_in  = fann_get_input(data, i);
  data_out = fann_get_output(data, i);
  calc_out = fann_run(ann, data_in);
  printf('XOR test (%f,%f) -> %f, should be %f, difference=%f\n', data_in(1), data_in(2), calc_out(1), data_out(1), abs(calc_out(1) - data_out(1)));
end

printf('Saving network.\n');

fann_save(ann, 'xor_float.net');

// fann_savetofixed(ann, 'xor_fixed.net');
// decimal_point = fann_get_parameters(ann, 'decimal_point');
// fann_train_on_data(data, 'save_train_to_fixed', 'xor_fixed.data', decimal_point);

printf("Cleaning up.\n");
fann_destroy_train(data); clear data;
fann_destroy(ann); clear fann;
