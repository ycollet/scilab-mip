lines(0);

X = (0:0.01:3*2*%pi)';
Y = sin(X);
Noise = 0.1*rand(Y,'normal');
Data_In = [X';(Y+Noise)']';
NoiseVar = stdev(Noise);

N = 10;
P = 0.1;
mtest_min = 4;
mtest_max = 100;

subplot(3,1,3);
xtitle('The signal','t','y');
plot(X,Y+Noise,'b-');

params = init_param();

//params = add_param(params, 'select', 1);
//params = add_param(params, 'greedy', 1);
//params = add_param(params, 'greedy_fast', 1);
//params = add_param(params, 'rand', 1);
params = add_param(params, 'starting_n', 2);
//params = add_param(params, 'mtest_min', -1);
//params = add_param(params, 'mtest_max', -1);
params = add_param(params, 'seed', 12345);
params = add_param(params, 'norm', 0);
params = add_param(params, 'bucket_size', 45);
//params = add_param(params, 'start_point', 0);
//params = add_param(params, 'end_point', 0);
params = add_param(params, 'userversion2', 1);
params = add_param(params, 'userup', 0);

[gt_a, gt_b, reg_d, reg_g, scat_d, scat_g] = gammatest(Data_In,N,P,params);

subplot(3,1,1);
drawlater;
xtitle('Gamma Test','Delta', 'Gamma');
plot(scat_d, scat_g, 'r.');
plot(reg_d, reg_g, 'g.');
drawnow;

printf('The R2 value = %f\n", gammatest_r2(Data_In,gt_a));
printf('The variance ration value = %f\n", gammatest_vr(Data_In,gt_a));

params = add_param(params, 'mtest_min', mtest_min);
params = add_param(params, 'mtest_max', mtest_max);

[gt_a, gt_b, reg_d, reg_g, scat_d, scat_g] = gammatest(Data_In,N,P,params);

subplot(3,1,2);
drawlater;
xtitle('MTest','neighbors', 'A');
t=(mtest_min+1:mtest_max)';
plot(t, gammatest_r2(Data_In,gt_a), 'k-');
drawnow;

printf('The R2 value = "); disp(gammatest_r2(Data_In,gt_a));
printf('The variance ratio value = "); disp(gammatest_vr(Data_In,gt_a));
