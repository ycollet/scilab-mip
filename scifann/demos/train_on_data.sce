path = get_absolute_file_path('train_on_data.sce');

filename = path + 'data/abelone.train';
//filename = path + 'data/census-house.train';
//filename = path + 'data/mushroom.train';
//filename = path + 'data/robot.train';
//filename = path + 'data/bank32fm.train';
//filename = path + 'data/diabetes.train';
//filename = path + 'data/parity13.train';
//filename = path + 'data/soybean.train';
//filename = path + 'data/bank32nh.train';
//filename = path + 'data/gene.train';
//filename = path + 'data/parity8.train';
//filename = path + 'data/thyroid.train';
//filename = path + 'data/building.train';
//filename = path + 'data/kin32fm.train';
//filename = path + 'data/pumadyn-32fm.train';
//filename = path + 'data/two-spiral.train';

max_epochs      = 1000;
desired_error   = 1e-5;
cascade_neurons = 100;
UseStandardLearning = %T;
UseCascadeLearning  = %F;

ann_train_data  = fann_readtrain_from_file(filename);

num_input  = fann_setup_train_data(ann_train_data,'num_input_train_data');
num_output = fann_setup_train_data(ann_train_data,'num_output_train_data');
td_length  = fann_setup_train_data(ann_train_data,'length_train_data');

//ann_train_data = fann_setup_train_data(ann_train_data,'subset_train_data',1,ceil(td_length/2));

//ann_train_data = fann_setup_train_data(ann_train_data,'scale_train_data',-1,1);

printf('Filename          = %s\n',filename);
printf('Number of inputs  = %d\n',num_input);
printf('Number of outputs = %d\n',num_output);
printf('Length of the data set = %d\n',td_length);

if UseStandardLearning then
  ann = fann_create('sparse',[num_input 2 num_output], 0.8);
  //ann = fann_create('standard',[num_input 2 num_output]);
  //ann = fann_create('shortcut',[num_input 2 num_output]);
  
  //ann = fann_setup_train_data(ann_train_data,'set_input_scaling_params',ann,-1,1);
  //ann = fann_setup_train_data(ann_train_data,'set_output_scaling_params',ann,-1,1);

  //ann = fann_randomize_weights(ann, -1, 1);
  //ann = fann_init_weights(ann, ann_train_data);

  ann = fann_set_parameters(ann,'activation_function_hidden','FANN_SIGMOID_SYMMETRIC');
  ann = fann_set_parameters(ann,'activation_function_output','FANN_LINEAR');

  ann = fann_set_parameters(ann,'training_algorithm','FANN_TRAIN_RPROP');
  ann = fann_set_parameters(ann,'train_error_function','FANN_ERRORFUNC_LINEAR');
  ann = fann_set_parameters(ann,'train_stop_function','FANN_STOPFUNC_MSE');

  ann = fann_reset_MSE(ann);
  t_start = getdate();
  ann = fann_train_on_data(ann, ann_train_data, max_epochs, desired_error);
  printf('End of the training phase after %d sec.\n',etime(getdate(),t_start));
end

if UseCascadeLearning then
  ann = fann_reset_MSE(ann);
  t_start = getdate();
  ann = fann_create('shortcut',[num_input num_output]);

  //ann = fann_randomize_weights(ann, -1, 1);
  ann = fann_init_weights(ann, ann_train_data);

  ann = fann_setup_train_data(ann_train_data,'set_input_scaling_params',ann,-1,1);
  ann = fann_setup_train_data(ann_train_data,'set_output_scaling_params',ann,-1,1);
  
  ann = fann_casc_train_on_data(ann,ann_train_data,cascade_neurons,desired_error);
  printf('End of the training phase after %d sec.\n',etime(getdate(),t_start));
end

ann = fann_setup_train_data(ann_train_data,'clear_scaling_params',ann);

ann = fann_test_data(ann, ann_train_data);
printf('Error after training: MSE      = %f\n', fann_get_MSE(ann));
printf('                      bit_fail = %d\n', fann_get_bit_fail(ann));

fann_destroy(ann); clear ann;
fann_destroy_train(ann_train_data); clear ann_train_data;
