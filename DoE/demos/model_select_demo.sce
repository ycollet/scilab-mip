my_handle = scf(100001);
clf(my_handle,'reset');
demo_viewCode('model_select_demo.sce');

lines(0);
old_funcprot = funcprot();
funcprot(0);

nb_var = 2;
Log    = %T;

model_2nd  = '1 x1 x2 x1*x2 x1^2 x2^2';

min_levels = [-4 -5];
max_levels = [5 6];

// Test of the factorial DoE on a linear model
H1  = doe_factorial(2);
H1_unnorm = unnorm_doe_matrix(H1, min_levels, max_levels);

RM1 = build_regression_matrix(H1_unnorm,model_2nd);

Y = sum(H1_unnorm,'c') + 0.1*rand(size(H1_unnorm,1),1);

printf('model selection using backward selection\n');
[model_res,coeff_res] = doe_model_bselect(nb_var,model_2nd,Y,Log);
printf('the model = %s\n', model_res);
printf('the coefficients:'); disp(coeff_res);

printf('model selection using forward selection\n');
[model_res,coeff_res] = doe_model_fselect(nb_var,model_2nd,Y,Log);
printf('the model = %s\n', model_res);
printf('the coefficients:'); disp(coeff_res);

funcprot(old_funcprot);

