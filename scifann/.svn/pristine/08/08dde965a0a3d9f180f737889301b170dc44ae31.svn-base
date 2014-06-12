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
desired_error      = 0.001;
max_epochs         = 500000;

ann = fann_create('standard', [num_input, num_neurons_hidden, num_output]);

ann = fann_set_parameters(ann, 'activation_function_hidden', 'FANN_SIGMOID_SYMMETRIC');
ann = fann_set_parameters(ann, 'activation_function_output', 'FANN_SIGMOID_SYMMETRIC');

fann_train_on_file(ann, 'data/xor.data', max_epochs, desired_error);
	
fann_save(ann, 'xor_float.net');

fann_destroy(ann); clear ann;
