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

num_layers         = 3;
num_neurons_hidden = 96;
desired_error      = 0.001;

printf('Creating network.\n');

train_data = fann_readtrain_from_file('data/robot.train");

num_inputs  = fann_setup_train_data(train_data, 'num_input_train_data');
num_outputs = fann_setup_train_data(train_data, 'num_output_train_data');

ann = fann_create('standard', [num_inputs, num_neurons_hidden, num_outputs]);

printf('Training network.\n');

ann = fann_set_parameters(ann, 'training_algorithm', 'FANN_TRAIN_INCREMENTAL');
ann = fann_set_parameters(ann, 'learning_momentum', 0.4);

ann = fann_train_on_data(ann, train_data, 3000, desired_error);

printf('Testing network.\n');

test_data = fann_readtrain_from_file('data/robot.test');

fann_reset_MSE(ann);

for i=1:fann_setup_train_data(test_data, 'length_train_data')
  data_in = fann_get_input(test_data,i);
  data_out = fann_get_output(test_data,i);
  fann_test(ann, data_in, data_out);
end

printf('MSE error on test data: %f\n', fann_get_MSE(ann));

printf('Saving network.\n');

fann_save(ann, 'robot_float.net');

printf('Cleaning up.\n');
fann_destroy_train(train_data); clear train_data;
fann_destroy_train(test_data); clear test_data;
fann_destroy(ann); clear ann;
