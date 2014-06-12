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

printf('Creating network.\n');

[res, ierr] = fileinfo('scaling.net')
if (ierr) then
  printf('ERROR: you must execute scaling_train.sce first !\n');
  abort();
end

ann = fann_read('scaling.net');

if (isempty(ann)) then
   printf('Error creating ann --- ABORTING.\n');
   abort();
end

fann_print_connections(ann);
fann_print_parameters(ann);

printf('Testing network.\n');

data = fann_readtrain_from_file('data/scaling.data');
data_in = fann_get_input(data, 1);
data_out = fann_get_output(data, 1);

ann = fann_setup_train_data(data, 'scale_input', ann, 0.5*ones(length(data_in),1));

for i=1:fann_setup_train_data(data, 'length_train_data')
  ann      = fann_reset_MSE(ann);
  data_in  = fann_get_input(data, i);
  data_out = fann_get_output(data, i);
  calc_out = fann_run(ann, data_in);
  data_out = fann_get_output(data,i);
  data_in  = fann_get_input(data,i);
  // fann_descale_output(ann, calc_out);
  printf('Result %f original %f error %f\n', calc_out(1), data_out(1), abs(calc_out(1) - data_out(1)));
end

printf('Cleaning up.\n');
fann_destroy_train(data); clear data;
fann_destroy(ann); clear ann;
