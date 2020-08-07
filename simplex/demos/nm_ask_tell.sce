getf('step_nelder_mead.sci');

function this = nmomd()
  this = mlist(['nmomd', 'ItMX', 'x0', 'x_init', 'f_init', 'upper','lower', 'kelley_restart', 'kelley_alpha', ...
                'simplex_relsize', 'log', 'stop', 'data_next', 'init', 'f_hist', 'x_hist']);

  this.ItMX = 100;
  this.x0   = [];
  this.kelley_restart  = %F;
  this.kelley_alpha    = 1e-4;
  this.simplex_relsize = 0.1;
  this.log             = %F;
  this.stop            = %F;
  this.upper           = 1e6;
  this.lower           = - 1e6;
  this.data_next       = [];
  this.init            = %T;
  this.x_init          = [];
  this.f_init          = [];
  this.x_hist          = [];
  this.f_hist          = [];
endfunction

function x = %nmomd_ask(this)
  if this.init then
    // We set the initial simplex
    for i=1:length(this.x0)+1
      this.x_init(:,i) = this.x0 + this.simplex_relsize*0.5* ((this.upper - this.lower) .* rand(size(this.x0,1),size(this.x0,2)) + this.lower);
    end
  end
  x = this.x_init;
endfunction

function this = %nmomd_tell(this, x, y)
  this.f_hist = [];
  this.x_hist = list();

  if this.init then
    [this.x_init, this.data_next, eval_func, f_hist, x_hist] = step_nelder_mead(y, x, [], 'init', this.log, this.kelley_restart, this.kelley_alpha);
    this.init = %F;
    this.f_hist = f_hist;
    this.x_hist = x_hist;
  else
    [this.x_init, this.data_next, eval_func, f_hist, x_hist] = step_nelder_mead(y, x, this.data_next, 'run', this.log, this.kelley_restart, this.kelley_alpha);
    this.f_hist = f_hist;
    this.x_hist = x_hist;
  end

  this.ItMX = this.ItMX - 1;
  this.stop = this.stop | (this.ItMX <= 0);
endfunction 

function [yopt, xopt] = %nmomd_best(this) 
  [xopt, yopt] = step_nelder_mead(this.f_init, this.x_init, this.data_next, 'exit', this.log, this.kelley_restart, this.kelley_alpha);
endfunction

////////////////////////////////////
// Definition of the test problem //
////////////////////////////////////

function Res = min_bd_branin()
Res = [-5 0]';
endfunction

function Res = max_bd_branin()
Res = [15 10]';
endfunction

function Res = opti_branin()
Res = [-%pi     12.275; ...
        %pi     12.275; ...
        9.42478  2.475]';
endfunction

function y = branin(x)
y = (x(2)-(5.1/(4*%pi^2))*x(1)^2+5/%pi*x(1)-6)^2+10*(1-1.0/(8*%pi))*cos(x(1))+10;
endfunction

ItMX = 100;

nmopt      = nmomd();
nmopt.ItMX = ItMX;
nmopt.kelley_restart  = %F;
nmopt.kelley_alpha    = 1e-4;
nmopt.simplex_relsize = 0.1;
nmopt.log             = %F;

nmopt.upper = max_bd_branin();
nmopt.lower = min_bd_branin();

nmopt.x0   = (Max - Min).*rand(size(Max,1),size(Max,2)) + Min;

y_min = [];
while ~nmopt.stop
  printf('nmopt running: iteration %d / %d - ', ItMX - nmopt.ItMX + 1, ItMX);
  x = ask(nmopt);
  y = [];
  for i=1:size(x,2)
    y(i) = branin(x(:,i));
  end

  y_min($+1) = min(y);
  printf(' fmin = %f\n', y_min($));
  nmopt = tell(nmopt, x, y);
end 

[f_opt, x_opt] = best(nmopt);

scf;
drawlater;
xgrid(2);

X = Min(1):(Max(1)-Min(1))/10:Max(1);
Y = Min(2):(Max(2)-Min(2))/10:Max(2);
Z = [];

for i=1:length(X)
  for j=1:length(Y)
    Z(i,j) = branin([X(i),Y(j)]);
  end
end

xset('fpf',' ');
contour(X,Y,Z,10);

xtitle('Evolution of the simplex','x1','x2');

for i=4:length(nmopt.x_hist)
  plot(nmopt.x_hist(i)(1)(1), nmopt.x_hist(i)(1)(2), 'ko');
  plot(nmopt.x_hist(i)(2)(1), nmopt.x_hist(i)(2)(2), 'ko');
  plot(nmopt.x_hist(i)(3)(1), nmopt.x_hist(i)(3)(2), 'ko');
  plot([nmopt.x_hist(i)(1)(1) nmopt.x_hist(i)(2)(1)], [nmopt.x_hist(i)(1)(2) nmopt.x_hist(i)(2)(2)], 'k-');
  plot([nmopt.x_hist(i)(2)(1) nmopt.x_hist(i)(3)(1)], [nmopt.x_hist(i)(2)(2) nmopt.x_hist(i)(3)(2)], 'k-');
  plot([nmopt.x_hist(i)(3)(1) nmopt.x_hist(i)(1)(1)], [nmopt.x_hist(i)(3)(2) nmopt.x_hist(i)(1)(2)], 'k-');
end
drawnow;

