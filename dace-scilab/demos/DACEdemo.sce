DACEPath = get_absolute_file_path('DACEdemo.sce');

loadmatfile(DACEPath + 'data1.mat');
theta = [10 10];
lob   = [0.1 0.1];
upb   = [20 20];

[dmodel,perf] = dacefit(S, Y, regpoly0, corrgauss, theta, lob, upb);

X = gridsamp([0 0; 100 100], [40 40]);
[YX] = predictor(X, dmodel);

X1 = matrix(X(:,1), 40, 40);
X2 = matrix(X(:,2), 40, 40);
YX = matrix(YX, size(X1,1), size(X1,2));

h = scf();
drawlater;
param3d1(S(:,1), S(:,2), list(Y,-1));
Color = 127*((YX - min(YX))/(max(YX) - min(YX))) + 1;
xset('colormap',hotcolormap(128));
mesh(X1,X2,YX,YX);
xtitle('Kriging interpolation','x1','x2','yx');
drawnow;

[emodel,perf] = dacefit(S, Y, regpoly0, corrgauss, 2);
[tmp,MSE] = predictor(X, dmodel);
printf('max MSE = %f - min MSE = %f\n', max(MSE), min(MSE));

h = scf();
drawlater;
mesh(X1,X2,matrix(MSE,size(X1,1),size(X1,2)));
xtitle('Variance of the kriging model','x1','x2','mse');
drawnow;

[y,dy] = predictor(S(1,:),dmodel);
printf('y = %f - should be 34.1\n',y);
printf('dy = '); disp(dy'); printf('should be [0.2526 0.1774]''\n');

[y,dy,mse,dmse] = predictor([50 50],dmodel);
printf('y = %f - should be 38.0610\n',y);
printf('dy = '); disp(dy'); printf('should be [-0.0141 -0.0431]''\n');
printf('mse = %f - should be 0.7526\n',mse);
printf('dmse = '); disp(dmse'); printf('should be [0.0004 0.0187]''\n');

