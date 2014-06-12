krig_path= get_absolute_file_path('plot_2d_krig.sce');

[X,Y] = meshgrid(-2:0.1:2,-2:0.1:2);
Model = readKrig(krig_path + '/krig_name.dat');
Total = prod(size(X));
Step  = ceil(Total/10);
wId   = waitbar(0);
Index = 0;
for i=1:size(X,1)
  for j=1:size(X,2)
    Index = Index + 1;
    if (~modulo(Index,Step)) then 
      waitbar(ceil(1000*Index/Total)/1000, wId);
    end;
    [Z_estim(i,j) Var_estim(i,j)] = computeKrig(Model, [X(i,j) Y(i,j)]);
  end
end
winclose(wId);

scf;
f = gcf();
Color = graycolormap(32);
f.color_map = Color(16:$,:);
drawlater;
subplot(2,1,1);
surf(X,Y,Z_estim);
xtitle('Estimation','x1','x2','F');
subplot(2,1,2);
surf(X,Y,abs(Var_estim));
xtitle('Variance','x1','x2','F');
drawnow;

