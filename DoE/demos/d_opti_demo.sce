my_handle = scf(100001);
clf(my_handle,'reset');
demo_viewCode('d_opti_demo.sce');

lines(0);
old_funcprot = funcprot();
funcprot(0);

Doe_size   = 25;
Doe_ItMX   = 3;
Doe_plevel = 0.05;
Log        = %T;
size_tabu_list = 6

model_lin        = doe_poly_model('lin',2);
model_poly_inter = doe_poly_model('inter',2,1);
model_poly_2nd   = doe_poly_model('poly',2,2);

model = model_poly_inter;

printf('Model selected:\n');
disp(model);

Doe_Upper = ones(1,2);
Doe_Lower = -1*ones(1,2);

[X, Y] = meshgrid(-1:0.1:1, -1:0.1:1);
Index = 1;
for i=1:size(X,1)
  for j=1:size(X,2)
    Doe_Db(Index,:) = [X(i,j), Y(i,j)];
    Index = Index + 1;
  end
end

Doe_Init = [];

[Doe_Comp, history] = doe_d_opti(Doe_Init, Doe_Db, Doe_size, model, Doe_Lower, Doe_Upper, Doe_ItMX, Doe_plevel, Log, size_tabu_list);

RM = build_regression_matrix(Doe_Comp, model);
// Computation of the variance
[X, Y] = meshgrid(-1:0.1:1, -1:0.1:1);
for i=1:size(X,1)
  for j=1:size(X,2)
    x_aux = [X(i,j) Y(i,j)]';
    Z(i,j) = var_regression_matrix(RM, x_aux, model, 1);
  end
end

drawlater;
if ~isempty(Doe_Init) then
  subplot(3,2,1);
  for i=1:size(Doe_Init,1)
    plot(Doe_Init(i,1), Doe_Init(i,2), 'ro');
  end
  xtitle('Plan initial','x1','x2');
end
if ~isempty(Doe_Init) then
  subplot(3,2,2);
  for i=1:size(Doe_Db,1)
    plot(Doe_Db(i,1), Doe_Db(i,2), 'ro');
  end
  xtitle('Base de points','x1','x2');
end

subplot(3,2,3)
surf(X,Y,Z);
xtitle('D-Opti DoE','x1','x2', 'Var');

subplot(3,2,4)
for i=1:size(Doe_Comp,1)
  plot(Doe_Comp(i,1), Doe_Comp(i,2), 'ro');
end

subplot(3,1,5)
t=1:length(history);
plot(t, history, 'k');
xtitle('History of the criterion','Iteration','Criterion');
drawnow;

printf('D criterion = %f\n', comp_d_opti_crit(Doe_Comp, model));
printf('A criterion = %f\n', comp_a_opti_crit(Doe_Comp, model));
printf('G criterion = %f\n', comp_g_opti_crit(Doe_Comp, model));
printf('O criterion = %f\n', comp_o_opti_crit(Doe_Comp, model));

funcprot(old_funcprot);

