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

num_neurons_hidden = 32;
desired_error      = 0.0001;
max_epochs         = 300;

printf('Creating network.\n');

train_data = fann_readtrain_from_file('data/mushroom.train');

ann = fann_create('standard', [fann_setup_train_data(train_data, 'num_input_train_data'), num_neurons_hidden, fann_setup_train_data(train_data, 'num_output_train_data')]);

printf('Training network.\n');

ann = fann_set_parameters(ann, 'activation_function_hidden', 'FANN_SIGMOID_SYMMETRIC_STEPWISE');
ann = fann_set_parameters(ann, 'activation_function_output', 'FANN_SIGMOID_STEPWISE');

// ann = fann_set_parameters(ann, 'training_algorithm', 'FANN_TRAIN_INCREMENTAL');

ann = fann_train_on_data(ann, train_data, desired_error, max_epochs);

printf('Testing network.\n');

test_data = fann_readtrain_from_file('data/mushroom.test');

ann = fann_reset_MSE(ann);

for i=1:fann_setup_train_data(test_data, 'length_train_data')
 data_in  = fann_get_input(test_data, i);
 data_out = fann_get_output(test_data, i);
 fann_test(ann, data_in, data_out);
end	

printf('MSE error on test data: %f\n', fann_get_MSE(ann));

printf('Saving network.\n');

fann_save(ann, 'mushroom_float.net');

printf('Cleaning up.\n');
fann_destroy_train(train_data); clear train_data;
fann_destroy_train(test_data); clear test_data;
fann_destroy(ann); clear ann;
