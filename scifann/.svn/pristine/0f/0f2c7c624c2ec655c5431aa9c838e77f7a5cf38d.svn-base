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
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA02111-1307USA
// 


num_neurons_hidden = 96;
desired_error = 0.001;

traindata = fann_readtrain_from_file('data/robot.train');
testdata = fann_readtrain_from_file('data/robot.test');

for momentum=0.0:0.1:0.7
  printf('============= momentum = %f =============\n', momentum);

  num_input = fann_setup_train_data(traindata, 'num_input_train_data');
  num_output = fann_setup_train_data(traindata, 'num_output_train_data');

  ann = fann_create('standard', [num_input, num_neurons_hidden, num_output]);
  
  ann = fann_set_parameters(ann, 'training_algorithm', 'FANN_TRAIN_INCREMENTAL');
  
  ann = fann_set_parameters(ann, 'learning_momentum', momentum);
  
  ann = fann_train_on_data(ann, traindata, 2000, desired_error);
  
  ann = fann_reset_MSE(ann);
  ann = fann_test_data(ann, traindata);
  printf('MSE error on train data: %f\n', fann_get_MSE(ann));
  ann = fann_reset_MSE(ann);
  ann = fann_test_data(ann, testdata);
  printf('MSE error on test data : %f\n', fann_get_MSE(ann));
  
  fann_print_connections(ann);
  fann_print_parameters(ann);

  fann_destroy(ann); clear ann;
end

fann_destroy_train(traindata); clear traindata;
fann_destroy_train(testdata); clear testdata;
