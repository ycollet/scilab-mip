NbPointsLearn = 1000;
NbPointsValid = 100;
desired_error = 1e-6;
max_epochs    = 5000;

Xmax = [2 2];
Xmin = [-2 -2];

function plot_netres(xData_valid,fData_valid,fData_estim)
  scf();
  xtitle('Result of the training','x1','x2','f');
  drawlater();
  for i=1:NbPointsValid
    param3d1(xData_valid(i,1), xData_valid(i,2), list(fData_valid(i),-2));
    param3d([xData_valid(i,1) xData_valid(i,1)], [xData_valid(i,2) xData_valid(i,2)], [fData_valid(i) fData_estim(i)]);
  end
  drawnow();
endfunction

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

printf('call to fann_create\n');
ann = fann_create('sparse',[2 2 1], 0.8);
// ann = fann_create('standard',[2 2 1]);
// ann = fann_create('shortcut',[2 2 1]);

ann = fann_set_parameters(ann,'activation_function_hidden','FANN_SIGMOID_SYMMETRIC');
ann = fann_set_parameters(ann,'activation_function_output','FANN_LINEAR');
ann = fann_randomize_weights(ann, -1, 1);

ann = fann_set_parameters(ann,'training_algorithm','FANN_TRAIN_RPROP');
ann = fann_set_parameters(ann,'train_error_function','FANN_ERRORFUNC_LINEAR');
ann = fann_set_parameters(ann,'train_stop_function','FANN_STOPFUNC_MSE');

printf('call to fann_train\n');
ann = fann_train_matrix(ann, xData_learn, fData_learn, desired_error, max_epochs);

valid_err = [];
for i=1:NbPointsValid
  fData_estim(i) = fann_test_matrix(ann, xData_valid(i,:));
  valid_err(i) = fData_estim(i) - fData_valid(i);
  printf('call to testfann %d / %d - f_estim = %f - f_meas = %f\n',i,NbPointsValid,fData_estim(i),fData_valid(i));
end

printf('variance of the modelization error: %f\n', stdev(valid_err));
printf('R2 for the validation data set = %f\n', 1 - (stdev(valid_err)/stdev(fData_valid))^2);

plot_netres(xData_valid,fData_valid,fData_estim);

printf('testing some functions:\n');
disp(fann_test_matrix(ann, xData_valid(i,:)));
disp(fann_run(ann, xData_valid(i,:)));

ann_train_data = fann_traindata_from_matrix(xData_learn,fData_learn);
ann = init_weights(ann, ann_train_data);
// fann_train                       Train one iteration with a set of inputs, and a set of desired outputs.
ann = fann_train(ann,xData_learn, fData_learn);
// fann_test                        Test with a set of inputs, and a set of desired outputs.
[ann,result] = test(ann,xData_learn(1,:), fData_learn(1));
printf("value of the neural network output after a learning phase on 1 point: %f / %f\n", result, fData_learn(1));
// fann_get_MSE                     Reads the mean square error from the network.
result = fann_get_MSE(ann);
printf("MSE = %f\n", result);
// fann_get_bit_fail                The number of fail bits; means the number of output neurons which differ more than the bit 
//                                  fail limit (see fann_get_bit_fail_limit, fann_set_bit_fail_limit).
result = fann_get_bit_fail(ann);
printf("bit_fail = %f\n", result);
// fann_reset_MSE                   Resets the mean square error from the network.
ann = fann_reset_MSE(ann);
// fann_train_on_data               Trains on an entire dataset, for a period of time.
ann = fann_train_on_data(ann, ann_train_data, desired_error, max_epochs);
result = fann_get_MSE(ann);
printf("train_on_data: MSE = %f\n", result);
// fann_train_epoch                 Train one epoch with a set of training data.
ann = fann_reset_MSE(ann);
ann = fann_train_epoch(ann, ann_train_data);
result = fann_get_MSE(ann);
printf("train_epoch: MSE = %f\n", result);
// fann_test_data                   Test a set of training data and calculates the MSE for the training data.
ann = fann_reset_MSE(ann);
ann = fann_test_data(ann, ann_train_data);
result = fann_get_MSE(ann);
printf("test_data: MSE = %f\n", result);

ann_train_data = fann_setup_train_data(ann_train_data,'shuffle_train_data');
ann = fann_setup_train_data(ann_train_data,'scale_train',ann);
ann = fann_setup_train_data(ann_train_data,'descale_train',ann);
ann = fann_setup_train_data(ann_train_data,'set_input_scaling_params',ann,-1,1);
ann = fann_setup_train_data(ann_train_data,'set_output_scaling_params',ann,-1,1);
ann = fann_setup_train_data(ann_train_data,'set_scaling_params',ann,-1,1,-1,1);
ann = fann_setup_train_data(ann_train_data,'clear_scaling_params',ann);
ann = fann_setup_train_data(ann_train_data,'scale_input',ann,[4 4]);
ann = fann_setup_train_data(ann_train_data,'scale_output',ann,4);
ann = fann_setup_train_data(ann_train_data,'descale_input',ann,[4 4]);
ann = fann_setup_train_data(ann_train_data,'descale_output',ann,4);
ann_train_data = fann_setup_train_data(ann_train_data,'scale_input_train_data',-1,1);
ann_train_data = fann_setup_train_data(ann_train_data,'scale_output_train_data',-1,1);
ann_train_data = fann_setup_train_data(ann_train_data,'scale_train_data',-1,1);
ann_train_data_2 = fann_setup_train_data(ann_train_data,'duplicate_train_data');
ann_train_data_3 = fann_setup_train_data(ann_train_data,'merge_train_data',ann_train_data);
ann_train_2 = fann_setup_train_data(ann_train_data_2,'subset_train_data',2,100);
result = fann_setup_train_data(ann_train_data,'length_train_data');
printf('length of ann_train_data: %d\n', result);

result = fann_setup_train_data(ann_train_data,'num_input_train_data');
printf('num_input_train_data: %d\n', result);

result = fann_setup_train_data(ann_train_data,'num_output_train_data');
printf('num_output_train_data: %d\n', result);

fann_setup_train_data(ann_train_data,'save_train','traindata_test.dat');
fann_setup_train_data(ann_train_data,'save_train','traindata_fixed_test.dat',3);

/////////////////////////////////////////////////////
// Generate a new learning and validation data set //
/////////////////////////////////////////////////////

////////////////////////////////////
// Generate the learning data set //
////////////////////////////////////

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

//////////////////////////////////////
// Generate the validation data set //
//////////////////////////////////////

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

ann_train_data = fann_traindata_from_matrix(xData_learn,fData_learn);
fann_setup_train_data(ann_train_data,'save_train','./traindata_test.dat');
fann_setup_train_data(ann_train_data,'save_train','./traindata_fixed_test.dat',3);

// Warning, works only with an initial net of "shortcut" type without hidden layers

// This test must emit an error message: wrong net type
//ann = fann_create('standard',[2 1]);
//ann = fann_casc_train_on_data(ann,ann_train_data,100,1e-4);

// This test must emit an error message: a hidden layer exists
//ann = fann_create('shortcut',[2 2 1]);
//ann = fann_casc_train_on_data(ann,ann_train_data,100,1e-4);
// This one works
ann = fann_create('shortcut',[2 1]);
ann = fann_casc_train_on_data(ann,ann_train_data,100,1e-4);

valid_err = [];
for i=1:NbPointsValid
  fData_estim(i) = test_matrix(ann, xData_valid(i,:));
  valid_err(i) = fData_estim(i) - fData_valid(i);
  printf('call to testfann %d / %d - f_estim = %f - f_meas = %f\n',i,NbPointsValid,fData_estim(i),fData_valid(i));
end

printf('variance of the modelization error: %f\n', stdev(valid_err));
printf('R2 for the validation data set = %f\n', 1 - (stdev(valid_err)/stdev(fData_valid))^2);

plot_netres(xData_valid,fData_valid,fData_estim);

// Warning, works only with an initial net of "shortcut" type without hidden layers
ann = fann_create('shortcut',[2 1]);
ann = fann_casc_train_on_file(ann,'./traindata_test.dat',100,1e-4);

valid_err = [];
for i=1:NbPointsValid
  fData_estim(i) = test_matrix(ann, xData_valid(i,:));
  valid_err(i) = fData_estim(i) - fData_valid(i);
  printf('call to testfann %d / %d - f_estim = %f - f_meas = %f\n',i,NbPointsValid,fData_estim(i),fData_valid(i));
end

printf('variance of the modelization error: %f\n', stdev(valid_err));
printf('R2 for the validation data set = %f\n', 1 - (stdev(valid_err)/stdev(fData_valid))^2);

plot_netres(xData_valid,fData_valid,fData_estim);

// fann_train_on_file               Does the same as fann_train_on_data, but reads the training data directly from a file.
ann = fann_train_on_file(ann, './traindata_test.dat', desired_error, max_epochs);
