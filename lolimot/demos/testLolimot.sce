//////////////////
// Test lolimot //
//////////////////

// Parameters
NbPts      = 100;
NbPtsValid = 10;
sigma      = 0.33;
nbpart     = 10;
maximp     = 0.05;
Log        = %T;
Vec        = %T; // Use vectorial estimation - faster - but uses a lot more memory

/////////////////////////////////////////
// Generation of the learning data set //
/////////////////////////////////////////

x1_learn = -5:10/sqrt(NbPts):5;
x2_learn = -5:10/sqrt(NbPts):5;

Index = 1;
for i=1:length(x1_learn)
  for j=1:length(x2_learn)
    y(Index) = x1_learn(i)^2 + x2_learn(j)^2;
    data_learn(Index,:) = [x1_learn(i) x2_learn(j) y(Index)];
    Index = Index + 1;
  end
end

///////////////////////////////////////////
// Generation of the validation data set //
///////////////////////////////////////////

data_valid = data_learn(1:NbPtsValid,:);

printf('Learning a %d data set\n',size(data_learn,1));

/////////////////////
// Learning models //
/////////////////////

t_global_start = getdate();

tic();
[lolModel,stat] = learn_lolimot(data_learn,sigma,nbpart,maximp,2,Vec,Log);
t = toc();

printf('%f s required for learning - final residual = %f\n',t,lolModel('residual'));

tic();
[lolModel_lin,stat] = learn_lolimot(data_learn,sigma,nbpart,maximp,2,Vec,Log,pLinear=%T);
t = toc();

printf('%f s required for learning a non exponential model - final residual = %f\n',t,lolModel_lin('residual'));

Z_est          = [];
Z_est_lin      = [];
Z_est_vect     = [];
Z_est_vect_lin = [];
Z_mes          = [];

X = [];
Y = [];

Index = 1;
for i=1:length(x1_learn)
  for j=1:length(x2_learn)
    X(i,j) = x1_learn(i);
    Y(i,j) = x2_learn(j);
    Z_mes(i,j) = data_learn(Index,$);
    Index = Index + 1;
  end
end

tic();
Index = 1;
for i=1:length(x1_learn)
  for j=1:length(x2_learn)
    Z_est(i,j)     = estim_lolimot(data_learn(Index,1:$-1),lolModel);
    Z_est_lin(i,j) = estim_lolimot(data_learn(Index,1:$-1),lolModel_lin);
    Index = Index + 1;
  end
end
printf('Lolimot: %.2f s needed to compute sequentially the estimations\n', toc());

tic();
Index = 1;
for i=1:length(x1_learn)
  for j=1:length(x2_learn)
    Z_est_vec(i,j)     = estim_vec_lolimot(data_learn(Index,1:$-1),lolModel);
    Z_est_vec_lin(i,j) = estim_vec_lolimot(data_learn(Index,1:$-1),lolModel_lin);
    Index = Index + 1;
  end
end
printf('Lolimot: %.2f s needed to compute vectorially the estimations\n', toc());

//////////////////////
// Plot the results //
//////////////////////

cmap = graycolormap(128);

h = scf();
h.color_map = cmap;

drawlater;
subplot(2,2,1);
surf(X,Y,Z_mes);
xtitle('Measurement data set','X','Y','Z_mes');

subplot(2,2,2);
surf(X,Y,Z_est);
xtitle('Estimation data set','X','Y','Z_est');

subplot(2,2,3);
surf(X,Y,Z_est_vec);
xtitle('Estimation data set - Vect','X','Y','Z_est_vec');
drawnow;

printf('Writing and reading a learnt lolimot model\n');

///////////////////////////////////
// Save and read again the model //
///////////////////////////////////

err = write_lolimot('lolimot.txt',lolModel);
[err, modelOut] = read_lolimot('lolimot.txt');

Z_read = [];

Index = 1;
for i=1:length(x1_learn)
  for j=1:length(x2_learn)
    Z_read(i,j) = estim_lolimot(data_learn(Index,1:$-1),modelOut);
    Index = Index + 1;
  end
end

drawlater;
subplot(2,2,4);
surf(X,Y,Z_read);
xtitle('Estimation data set from saved lolimot','X','Y','Z_read');
drawnow;

//////////////////////
// Plot the results //
//////////////////////

h0 = scf();
h0.color_map = cmap;

drawlater;
subplot(2,2,1);
surf(X,Y,Z_mes);
xtitle('Measurement data set','X','Y','Z_mes');

subplot(2,2,2);
surf(X,Y,Z_est_lin);
xtitle('Estimation data set - piecewise linear membership function','X','Y','Z_est');

subplot(2,2,3);
surf(X,Y,Z_est_vec_lin);
xtitle('Estimation data set - piecewise linear membership function - Vect','X','Y','Z_est_vec');

printf('Writing and reading a learnt lolimot model\n');

subplot(2,2,4);
surf(X,Y,Z_read);
xtitle('Estimation data set from saved lolimot - piecewise linear membership function','X','Y','Z_read');
drawnow;

h1 = scf();
drawlater;
plot_lolimot_part(lolModel,'First model: Cuts');
drawnow;

h2 = scf();
h2.color_map = cmap;

for i=1:length(stat)
  NbPart(i)   = i;
  Residual(i) = stat(i)(2);
  Time(i)     = stat(i)(3);
  DimCut(i)   = stat(i)(5);
end

drawlater;
subplot(2,2,1);
plot2d(NbPart,Residual);
xtitle('Residual vs Nb partition','Nb partition','Residual');
subplot(2,2,2);
plot2d(NbPart,Time);
xtitle('Time vs Nb partition','Nb partition','Time');
subplot(2,1,2);
plot2d(NbPart,DimCut);
xtitle('Dimension cut vs Nb partition','Nb partition','No dimension');
drawnow;

/////////////////////////////////////////////////
// Perform other learning phase on a new model //
/////////////////////////////////////////////////

nbpart = 5;

printf('learning a 5 partitions / 3 cuts per partition lolimot model\n');

tic();
[lolModel,stat] = learn_lolimot(data_learn,sigma,nbpart,maximp,3,Vec,Log);
t = toc();

printf('%f s required for learning - final residual = %f\n',t,lolModel('residual'));

h3 = scf();

plot_lolimot_part(lolModel,'Second model - phase 1: Cuts');

nbpart = 5;

printf('adding 5 partitions / 2 cuts per partition to the current lolimot model\n');

tic();
[lolModel,stat] = learn_lolimot(data_learn,sigma,nbpart,maximp,2,Vec,Log,lolModel);
t = toc();

printf('%f s required for learning - final residual = %f\n',t,lolModel('residual'));

h4 = scf();

drawlater;
plot_lolimot_part(lolModel,'Second model - phase 2: Cuts');
drawnow;

/////////////////////////////////////////////////
// Learning a model with a validation data set //
/////////////////////////////////////////////////

nbpart = 10;

printf('Learning a 10 partition / 2 cuts per partition lolimot model with a validation data set\n');

tic();
[modelOut,stat] = learn_valid_lolimot(data_learn,data_valid,sigma,nbpart,maximp,2,Vec,Log);
t = toc();

printf('%f s required for learning - final residual = %f\n',t,lolModel('residual'));

h5 = scf();

drawlater;
plot_lolimot_part(lolModel,'Third model: Cuts');
drawnow;

printf('computation of a derivative point\n');

df = estim_der_lolimot([0;0],lolModel);

printf('the derivative value:'); disp(df');

//////////////////////////////////////////////////
// Performing optimization on the lolimot model //
//////////////////////////////////////////////////

// With an exponential membership function

printf('performing the optimization of a lolimot model - exponential membership function\n');

nbpart = 20;

printf('learning a 20 partitions / 2 cuts per partition lolimot model for optimization\n');

[lolModel,stat] = learn_lolimot(data_learn,sigma,nbpart,maximp,2,Vec,Log,pLinear=%F);

deff('[f, df, ind_out] = f_optim(x,ind_in)','f=estim_lolimot(x,lolModel); df=estim_der_lolimot(x,lolModel); ind_out = ind_in;');

[f_opt,x_opt] = optim(f_optim,[4;4]);

printf('Final cost: %f - x_opt = \n',f_opt); disp(x_opt');

// With a piecewise linear exponential membership function

printf('performing the optimization of a lolimot model - piecewise linear membership function\n');

nbpart = 20;

printf('learning a 20 partitions / 2 cuts per partition lolimot model for optimization\n');

[lolModel,stat] = learn_lolimot(data_learn,sigma,nbpart,maximp,2,Vec,Log,pLinear=%T);

deff('[f, df, ind_out] = f_optim(x,ind_in)','f=estim_lolimot(x,lolModel); df=estim_der_lolimot(x,lolModel); ind_out = ind_in;');

[f_opt,x_opt] = optim(f_optim,[4;4]);

printf('Final cost: %f - x_opt = \n',f_opt); disp(x_opt');

///////////////////////////////////
// Testing the derivative model: //
///////////////////////////////////

nbpart = 10;
maximp = 0.05;
NbPts  = 200;

X = -%pi/2:%pi/NbPts:%pi/2;
Y = sin(X);
dY = cos(X);

[lolModelSin,stat] = learn_lolimot([X',Y'],sigma,nbpart,maximp,2,Vec,Log,pLinear=%F);
[lolModelSin_lin,stat] = learn_lolimot([X',Y'],sigma,nbpart,maximp,2,Vec,Log,pLinear=%T);

Y_estim     = []; Y_estim_lin     = [];
Y_estim_der = []; Y_estim_der_lin = [];

for i=1:length(X)
  Y_estim(i)     = estim_lolimot(X(i),lolModelSin);
  Y_estim_lin(i) = estim_lolimot(X(i),lolModelSin_lin);
  Y_estim_der(i) = estim_der_lolimot(X(i),lolModelSin);
  Y_estim_der_lin(i) = estim_der_lolimot(X(i),lolModelSin_lin);
end

h = scf();

drawlater;
subplot(2,1,1);
plot(X,Y,'k-');
plot(X,Y_estim+0.1,'r-');
plot(X,Y_estim_lin-0.1,'g-');
legend(['Measure','Estimations','Estimations lin']);
xtitle('The estimation model','X','Y');
subplot(2,1,2);
plot(X(1:$-1)+0.5*(X(2)-X(1)),-(Y_estim(1:$-1)-Y_estim(2:$))/(X(2)-X(1)),'r-');
plot(X,Y_estim_der,'g-');
plot(X,dY,'k-');
plot(X,Y_estim_der_lin,'b-');
legend(['Finite diff of estim','lolimot der.','anal. der','lolimot lin der.']);
xtitle('Test the derivative model','X','Y');
drawnow;

///////////////////////////////////////////
// Test of the partial derivatives model //
///////////////////////////////////////////

if MSDOS then
  unix('del test_sin_f.sci');
  unix('del test_sin_df.sci');
  unix('del test_sin_lin_f.sci');
  unix('del test_sin_lin_df.sci');
else
  unix('rm test_sin_f.sci');
  unix('rm test_sin_df.sci');
  unix('rm test_sin_lin_f.sci');
  unix('rm test_sin_lin_df.sci');
end

res = export_model('test_sin_f',lolModelSin,'scilab');
res = export_der_model('test_sin_df',lolModelSin,'scilab');
res = export_model('test_sin_lin_f',lolModelSin_lin,'scilab');
res = export_der_model('test_sin_lin_df',lolModelSin_lin,'scilab');

clear test_sin_f; getf('test_sin_f.sci');
clear test_sin_df; getf('test_sin_df.sci');
clear test_sin_lin_f; getf('test_sin_lin_f.sci');
clear test_sin_lin_df; getf('test_sin_lin_df.sci');

Y_estim     = []; Y_estim_lin     = [];
Y_estim_der = []; Y_estim_der_lin = [];

for i=1:length(X)
  Y_estim(i)     = test_sin_f(X(i));
  Y_estim_lin(i) = test_sin_lin_f(X(i));
  Y_estim_der(i) = test_sin_df(X(i));
  Y_estim_der_lin(i) = test_sin_lin_df(X(i));
end

h = scf();
drawlater;
subplot(2,1,1);
plot(X,Y,'k-');
plot(X,Y_estim+0.1,'r-');
plot(X,Y_estim_lin-0.1,'g-');
legend(['Measure','Estimations','Estimations lin']);
xtitle('Generated function - The estimation model','X','Y');
subplot(2,1,2);
plot(X(1:$-1)+0.5*(X(2)-X(1)),-(Y_estim(1:$-1)-Y_estim(2:$))/(X(2)-X(1)),'r-');
plot(X,Y_estim_der+0.1,'g-');
plot(X,dY,'k-');
plot(X,Y_estim_der_lin-0.1,'b-');
legend(['Finite diff of estim','lolimot der.','anal. der','lolimot lin der.']);
xtitle('Generated function - Test the derivative model','X','Y');
drawnow;

t_global_end = getdate();

printf('Global elapsed time = %.2f s\n', etime(t_global_end,t_global_start));
