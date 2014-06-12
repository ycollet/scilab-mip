my_handle = scf(100001);
clf(my_handle,'reset');
demo_viewCode('classic_doe_demo.sce');

lines(0);
old_funcprot = funcprot();
funcprot(0);

model_lin   = '1 x1 x2';
model_2nd   = '1 x1 x2 x1^2 x2^2';
model_inter = '1 x1 x2 x1*x2';

x = -1:0.1:1;
y = -1:0.1:1;
[X,Y] = meshgrid(x,y);

////////////////////////////
// Test of classicals DoE //
////////////////////////////

// Test of the factorial DoE on a linear model
H1  = doe_factorial(2);
RM1 = build_regression_matrix(H1,model_lin);
// Test of the star DoE
H2  = doe_star(2);
RM2 = build_regression_matrix(H2,model_lin);
// Test of the factorial face centered
H3 = doe_union(H1, sqrt(2)*H2);
RM3 = build_regression_matrix(H3,model_2nd);
// Test of the factorial DoE on an interaction model
H4  = doe_factorial(2);
RM4 = build_regression_matrix(H1,model_inter);
// Computation of the variance
for i=1:size(X,1)
  for j=1:size(X,2)
    x_aux   = [X(i,j) Y(i,j)];
    Z1(i,j) = var_regression_matrix(RM1, x_aux, model_lin, 1);
    Z2(i,j) = var_regression_matrix(RM2, x_aux, model_lin, 1);
    Z3(i,j) = var_regression_matrix(RM3, x_aux, model_2nd, 1);
    Z4(i,j) = var_regression_matrix(RM4, x_aux, model_inter, 1);
  end
end

drawlater;

subplot(2,2,1);
f = gcf();
f.color_map = graycolormap(128);
surf(X,Y,Z1);
xtitle('Factorial DoE on a linear model','x1','x2','Var');

subplot(2,2,2);
f = gcf();
f.color_map = graycolormap(128);
surf(X,Y,Z2);
xtitle('Star DoE','x1','x2','Var');

subplot(2,2,3);
f = gcf();
f.color_map = graycolormap(128);
surf(X,Y,Z3);
xtitle('Factorial face centered DoE','x1','x2','Var');

subplot(2,2,4);
f = gcf();
f.color_map = graycolormap(128);
surf(X,Y,Z4);
xtitle('Factorial DoE on a model with interaction','x1','x2','Var');
drawnow;

funcprot(old_funcprot);

