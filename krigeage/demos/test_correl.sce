//
//funcname = 'correl_cubic';
//p(1) = 0.1;
//p(2) = 1;

//
//funcname = 'correl_exp';
//p(1) = 1;
//p(2) = 0.1; // p(1)*exp(p(2)*x)

//
//funcname = 'correl_expg';
//p(1) = 1;
//p(2) = -0.1;
//p(3) = 2; // p(1)*exp(p(2)*x^p(3));

//
//funcname = 'correl_gauss';
//p(1) = 1;
//p(2) = 0.1; // p(1)*exp(-p(2)*x^2.0);

//
//funcname = 'correl_lin';
//p(1) = 0.2; // 1.0-p(1).*x;
//p(2) = 1;
//
//funcname = 'correl_materm';
//p(1) = 1; // Sigma > 0
//p(2) = 1; // nu > 0
//p(3) = 1; // rho > 0

//
//funcname = 'correl_sinus';
//p(1) = 1;
//p(2) = 1;
//p(3) = 1;
//p(4) = 1; // p(4)*sin(p(1)*x+p(2))+p(3);

//
//funcname = 'correl_spherical';
//p(1) = 0.1;

//
funcname = 'correl_spline';
p(1) = 0.1;
p(2) = 2;

x=0:0.1:10;

y = [];
for i=1:length(x)
  y(i) = eval(funcname+'(x(i),p)');
end

plot(x,y,'k-');
xtitle(funcname,'delta x','correl');
