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

function ann_out = train_on_steepness_file(ann_in, filename, max_epochs, desired_error, steepness_start, steepness_step, steepness_end)
  
  data = fann_readtrain_from_file(filename);
  
  ann_in = fann_set_parameters(ann_in, 'activation_steepness_hidden', steepness_start);
  ann_in = fann_set_parameters(ann_in, 'activation_steepness_output', steepness_start);
  
  for i=1:max_epochs
    // train
    ann_in = fann_train_epoch(ann_in, data);
    _error = fann_get_MSE(ann_in);
    
    if (_error < desired_error) then
      steepness_start = steepness_start + steepness_step;
      if (steepness_start <= steepness_end) then
        printf('Steepness: %f\n', steepness_start);
        ann_in = fann_set_parameters(ann_in, 'activation_steepness_hidden', steepness_start);
        ann_in = fann_set_parameters(ann_in, 'activation_steepness_output', steepness_start);
      else
        break;
      end
    end
  end
  fann_destroy_train(data); clear data;
  ann_out = ann_in;
endfunction

num_input          = 2;
num_output         = 1;
num_neurons_hidden = 3;
desired_error      = 0.001;
max_epochs         = 500000;

ann = fann_create('standard', [num_input, num_neurons_hidden, num_output]);

data = fann_readtrain_from_file('data/xor.data');

ann = fann_set_parameters(ann, 'activation_function_hidden', 'FANN_SIGMOID_SYMMETRIC');
ann = fann_set_parameters(ann, 'activation_function_output', 'FANN_SIGMOID_SYMMETRIC');

ann = fann_set_parameters(ann, 'training_algorithm', 'FANN_TRAIN_QUICKPROP');

ann = train_on_steepness_file(ann, 'data/xor.data', max_epochs, desired_error, 1.0, 0.1, 20.0);

ann = fann_set_parameters(ann, 'activation_function_hidden', 'FANN_THRESHOLD_SYMMETRIC');
ann = fann_set_parameters(ann, 'activation_function_output', 'FANN_THRESHOLD_SYMMETRIC');

for i=1:fann_setup_train_data(data, 'length_train_data')
  data_in  = fann_get_input(data, i);
  data_out = fann_get_output(data, i);
  calc_out = fann_run(ann, data_in);
  printf('XOR test (%f, %f) -> %f, should be %f, difference=%f\n', data_in(1), data_in(2), calc_out(1), data_out(1), abs(calc_out(1) - data_out(1)));
end

fann_save(ann, 'xor_float.net');

fann_destroy(ann); clear ann;
fann_destroy_train(data); clear data;
