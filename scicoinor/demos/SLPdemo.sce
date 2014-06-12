lines(0);

global NbFobj;
NbFobj = 0;

clear f;
clear df;
clear ineqconstraint;
clear df_ineqconstraint;
clear eqconstraint;
clear df_eqconstraint;

f = [];
df = [];
ineqconstraint = [];
df_ineqconstraint = [];
eqconstraint = [];
df_eqconstraint = [];
Integer_Variables = [];
x0 = [];

UseProblem = 1;

// MINLP Test Problems

if UseProblem==1 then
  function y = f(x)
    y = MINLPTP1obj(x);
  endfunction
  
  function y = df(x)
    y = derivative(f,x);
  endfunction
  
  function y = ineqconstraint(x)
    y = MINLPTP1ineq(x);
  endfunction
  
  function y = df_ineqconstraint(x)
    y = derivative(ineqconstraint,x);
  endfunction
  
  if typeof(MINLPTP1eq)=='function' then
    function y = eqconstraint(x)
      y = MINLPTP1eq(x);
    endfunction
    
    function y = df_eqconstraint(x)
      y = derivative(eqconstraint,x);
    endfunction
  end
  
  upper = max_bd_MINLPTP1();
  lower = min_bd_MINLPTP1();
  //x0    = MINLPTP1_x_init();
  x0    = (upper - lower).*rand(size(upper,1),size(upper,2)) + lower;
  Integer_Variables = MINLPTP1_mip_var();
  tmp = ascii('C')*ones(size(x0,1),size(x0,2));
  tmp(Integer_Variables) = ascii('I')*ones(length(Integer_Variables),1);
  Integer_Variables = ascii(tmp);
end

if UseProblem==2 then
  function y = f(x)
    y = MINLPTP2obj(x);
  endfunction
  
  function y = df(x)
    y = derivative(f,x);
  endfunction
  
  function y = ineqconstraint(x)
    y = MINLPTP2ineq(x);
  endfunction
  
  function y = df_ineqconstraint(x)
    y = derivative(ineqconstraint,x);
  endfunction
  
  if typeof(MINLPTP2eq)=='function' then
    function y = eqconstraint(x)
      y = MINLPTP2eq(x);
    endfunction
    
    function y=df_eqconstraint(x)
      y = derivative(eqconstraint,x);
    endfunction
  end
  upper = max_bd_MINLPTP2();
  lower = min_bd_MINLPTP2();
  x0    = MINLPTP2_x_init();
  Index = find(upper~=%inf);
  x0(Index) = (upper(Index) - lower(Index)).*rand(length(Index),1) + lower(Index);
  Integer_Variables = MINLPTP2_mip_var();
  tmp = ascii('C')*ones(size(x0,1),size(x0,2));
  tmp(Integer_Variables) = ascii('I')*ones(length(Integer_Variables),1);
  Integer_Variables = ascii(tmp);
end

if UseProblem==3 then
  function y = f(x)
    y = MINLPTP3obj(x);
  endfunction
  
  function y = df(x)
    y = derivative(f,x);
  endfunction
  
  function y = ineqconstraint(x)
    y = MINLPTP3ineq(x);
  endfunction
  
  function y = df_ineqconstraint(x)
    y = derivative(ineqconstraint,x);
  endfunction
  
  if typeof(MINLPTP3eq)=='function' then
    function y = eqconstraint(x)
      y = MINLPTP3eq(x);
    endfunction
    
    function y = df_eqconstraint(x)
      y = derivative(eqconstraint,x);
    endfunction
  end
  upper = max_bd_MINLPTP3();
  lower = min_bd_MINLPTP3();
  x0    = MINLPTP3_x_init();
  Integer_Variables = MINLPTP3_mip_var();
  tmp = ascii('C')*ones(size(x0,1),size(x0,2));
  tmp(Integer_Variables) = ascii('I')*ones(length(Integer_Variables),1);
  Integer_Variables = ascii(tmp);
end

if UseProblem==4 then
  // NOK - pb convergence - ne pas utiliser - que des variables entieres?
  function y = f(x)
    y = MINLPASAADI1obj(x);
  endfunction
  
  function y = df(x)
    y = derivative(f,x);
  endfunction
  
  function y = ineqconstraint(x)
    y = MINLPASAADI1ineq(x);
  endfunction
  
  function y = df_ineqconstraint(x)
    y = derivative(ineqconstraint,x);
  endfunction
  
  if typeof(MINLPASAADI1eq)=='function' then
    function y = eqconstraint(x)
      y = MINLPASSADI1eq(x);
    endfunction
    
    function y = df_eqconstraint(x)
      y = derivative(eqconstraint,x);
    endfunction
  end
  upper = max_bd_MINLPASAADI1();
  lower = min_bd_MINLPASAADI1();
  x0    = MINLPASAADI1_x_init();
  Integer_Variables = MINLPASAADI1_mip_var();
  tmp = ascii('C')*ones(size(x0,1),size(x0,2));
  tmp(Integer_Variables) = ascii('I')*ones(length(Integer_Variables),1);
  Integer_Variables = ascii(tmp);
end

if UseProblem==5 then
  // NOK - pb convergence ?
  function y = f(x)
    y = MINLPASAADI2obj(x);
  endfunction
  
  function y = df(x)
    y = derivative(f,x);
  endfunction
  
  function y = ineqconstraint(x)
    y = MINLPASAADI2ineq(x);
  endfunction
  
  function y = df_ineqconstraint(x)
    y = derivative(ineqconstraint,x);
  endfunction
  
  if typeof(MINLPASAADI2eq)=='function' then
    function y = eqconstraint(x)
      y = MINLPASSADI2eq(x);
    endfunction
    
    function y = df_eqconstraint(x)
      y = derivative(eqconstraint,x);
    endfunction
  end
  upper = max_bd_MINLPASAADI2();
  lower = min_bd_MINLPASAADI2();
  x0    = MINLPASAADI2_x_init();
  x0    = (upper - lower).*rand(size(upper,1),size(upper,2)) + lower;
  Integer_Variables = MINLPASAADI2_mip_var();
  tmp = ascii('C')*ones(size(x0,1),size(x0,2));
  tmp(Integer_Variables) = ascii('I')*ones(length(Integer_Variables),1);
  Integer_Variables = ascii(tmp);
end

if UseProblem==6 then
  //Presque OK - convergence fine tuning ?
  function y = f(x)
    y = MINLPASAADI3obj(x);
  endfunction
  
  function y = df(x)
    y = derivative(f,x);
  endfunction
  
  function y = ineqconstraint(x)
    y = MINLPASAADI3ineq(x);
  endfunction
  
  function y = df_ineqconstraint(x)
    y = derivative(ineqconstraint,x);
  endfunction
  
  if typeof(MINLPASAADI3eq)=='function' then
    function y = eqconstraint(x)
      y = MINLPASSADI3eq(x);
    endfunction
    
    function y = df_eqconstraint(x)
      y = derivative(eqconstraint,x);
    endfunction
  end
  upper = max_bd_MINLPASAADI3();
  lower = min_bd_MINLPASAADI3();
  x0    = MINLPASAADI3_x_init();
  Integer_Variables = MINLPASAADI3_mip_var();
  tmp = ascii('C')*ones(size(x0,1),size(x0,2));
  tmp(Integer_Variables) = ascii('I')*ones(length(Integer_Variables),1);
  Integer_Variables = ascii(tmp);
end

if UseProblem==7 then
  // Presque OK - Fine tune convergence ?
  function y = f(x)
    y = MINLP2DEXobj(x);
  endfunction
  
  function y = df(x)
    y = derivative(f,x);
  endfunction
  
  function y = ineqconstraint(x)
    y = MINLP2DEXineq(x);
  endfunction
  
  function y = df_ineqconstraint(x)
    y = derivative(ineqconstraint,x);
  endfunction
  
  if typeof(MINLP2DEXeq)=='function' then
    function y = eqconstraint(x)
      y = MINLP2DEXeq(x);
    endfunction
    
    function y = df_eqconstraint(x)
      y = derivative(eqconstraint,x);
    endfunction
  end
  upper = max_bd_MINLP2DEX();
  lower = min_bd_MINLP2DEX();
  x0    = MINLP2DEX_x_init();
  Integer_Variables = MINLP2DEX_mip_var();
  tmp = ascii('C')*ones(size(x0,1),size(x0,2));
  tmp(Integer_Variables) = ascii('I')*ones(length(Integer_Variables),1);
  Integer_Variables = ascii(tmp);
end

if UseProblem==8 then
  // NOK - no constraints, so doesn't work with SLP
  function y = f(x)
    y=MINLPGTDobj(x);
  endfunction
  
  function y = df(x)
    y = derivative(f,x);
  endfunction
  
  if typeof(MINLPGTDeq)=='function' then
    function y = ineqconstraint(x)
      y = MINLPGTDineq(x);
    endfunction
    
    function y = df_ineqconstraint(x)
      y = derivative(ineqconstraint,x);
    endfunction
  end
  if ~isempty(MINLPGTDeq) then
    function y = eqconstraint(x)
      y = MINLPGTDeq(x);
    endfunction
    
    function y = df_eqconstraint(x)
      y = derivative(eqconstraint,x);
    endfunction
  end
  upper = max_bd_MINLPGTD();
  lower = min_bd_MINLPGTD();
  x0    = MINLPGTD_x_init();
  Integer_Variables = MINLPGTD_mip_var();
  tmp = ascii('C')*ones(size(x0,1),size(x0,2));
  tmp(Integer_Variables) = ascii('I')*ones(length(Integer_Variables),1);
  Integer_Variables = ascii(tmp);
end

if UseProblem==9 then
  // NOK
  function y = f(x)
    y = MINLPAVGAS1obj(x); 
  endfunction  
  
  function y = df(x)
    y = derivative(f,x);
  endfunction

  function y = ineqconstraint(x)
    y = MINLPAVGAS1ineq(x);
  endfunction
  
  function y = df_ineqconstraint(x)
    y = derivative(ineqconstraint,x);
  endfunction
  
  if typeof(MINLPAVGAS1eq)=='function' then
    function y = eqconstraint(x)
      y = MINLPAVGAS1eq(x);
    endfunction
    
    function y = df_eqconstraint(x)
      y = derivative(eqconstraint,x);
    endfunction
  end
  upper = max_bd_MINLPAVGAS1();
  lower = min_bd_MINLPAVGAS1();
  x0    = MINLPAVGAS1_x_init();
  Integer_Variables = MINLPAVGAS1_mip_var();
  tmp = ascii('C')*ones(size(x0,1),size(x0,2));
  tmp(Integer_Variables) = ascii('I')*ones(length(Integer_Variables),1);
  Integer_Variables = ascii(tmp);
end

if UseProblem==10 then
  // NOK
  function y = f(x)
    y = MINLPAVGAS2obj(x);
  endfunction
  
  function y = df(x)
    y = derivative(f,x);
  endfunction
  
  function y = ineqconstraint(x)
    y = MINLPAVGAS2ineq(x);
  endfunction
  
  function y = df_ineqconstraint(x)
    y = derivative(ineqconstraint,x);
  endfunction
  
  if typeof(MINLPAVGAS2eq)=='function' then
    function y = eqconstraint(x)
      y = MINLPAVGAS2eq(x);
    endfunction
    
    function y = df_eqconstraint(x)
      y = derivative(eqconstraint,x);
    endfunction
  end
  upper = max_bd_MINLPAVGAS2();
  lower = min_bd_MINLPAVGAS2();
  x0    = MINLPAVGAS2_x_init();
  Integer_Variables = MINLPAVGAS2_mip_var();
  tmp = ascii('C')*ones(size(x0,1),size(x0,2));
  tmp(Integer_Variables) = ascii('I')*ones(length(Integer_Variables),1);
  Integer_Variables = ascii(tmp);
end

// LP Test Problems

if UseProblem==11 then
  function y = f(x)
    y = x(1)^2 + x(2)^2 - 16*x(1) - 10*x(2);
    global NbFobj; ...
    NbFobj = NbFobj + 1;
  endfunction
  
  function y=df(x)
    y=[2*x(1) - 16; ...
       2*x(2) - 10];
  endfunction
  
  function y = ineqconstraint(x)
    y(1,1) = - 11 + x(1)^2 - 6*x(1) + 4*x(2);
    y(2,1) = - x(1)*x(2) + 3*x(2) + exp(x(1) - 3) - 1;
    y(3,1) = - x(1);
    y(4,1) = - x(2);
  endfunction
  
  function y = df_ineqconstraint(x)
    y(1,1) = 2*x(1) - 6;
    y(2,1) = 4;
    y(1,2) = -x(2) + exp(x(1) - 3);
    y(2,2) = -x(1) + 3;
    y(1,3) = -1;
    y(2,3) = 0;
    y(1,4) = 0;
    y(2,4) = -1;
  endfunction
  
  eqconstraint    = [];
  df_eqconstraint = [];
  upper       = [15;9];
  lower       = [0; 0];
  x0          = [4; 3];
end

if UseProblem==12 then
  function y = f(x)
    y = 5*x(1)^2-3*x(2)^2;
    global NbFobj; ...
    NbFobj = NbFobj + 1;
  endfunction
  
  function y = df(x)
    y(1) = 10*x(1);
    y(2) = -6*x(2);
  endfunction
  
  function y = ineqconstraint(x)
    y(1) = -x(1);
    y(2) = -x(2);
  endfunction
  
  function y = df_ineqconstraint(x)
    y(1,1) = -1;
    y(2,1) = 0;
    y(1,2) = 0;
    y(2,2) = -1;
  endfunction
  
  eqconstraint    = [];
  df_eqconstraint = [];
  upper = [4;4];
  lower = [-4;-4];
  x0 = [1;1];
end

if UseProblem==13 then
  function y = f(x)
    y = x(1)^2+x(2)^2;
    global NbFobj;
    NbFobj = NbFobj + 1;
  endfunction
  
  //function y = df(x)
  //  y = derivative(f,x);
  //endfunction
  function y = df(x)
    y(1,1) = 2*x(1);
    y(2,1) = 2*x(2);
  endfunction
  
  function y = eqconstraint(x)
    y(1) = - x(1)^2 - x(2)^2 + 9*x(2) - 4.25;
  endfunction
  
  function y = df_eqconstraint(x)
    y = derivative(eqconstraint,x);
  endfunction
  
  function y = df_eqconstraint(x)
    y(1,1) = -2*x(1);
    y(2,1) = -2*x(2) + 9;
  endfunction
  
  ineqconstraint = [];
  df_ineqconstraint = [];
  upper = [4;4];
  lower = [-4;-4];
  x0    = [2;3.9]; // Non feasible starting point
end

if UseProblem==14 then
  // constr_pb_2
  function y = f(x)
    y = x(1)^2 + x(2)^2 - 16*x(1) - 10*x(2);
    global NbFobj;
    NbFobj = NbFobj + 1;
  endfunction
  
  function y = df(x)
    y = derivative(f,x);
  endfunction
  
  function y = ineqconstraint(x)
    y(1) = - 11 + x(1)^2 - 6*x(1) + 4*x(2);
    y(2) = - x(1)*x(2) + 3*x(2) + exp(x(1) - 3) - 1;
  endfunction
  
  function y = df_ineqconstraint(x)
    y = derivative(ineqconstraint,x);
  endfunction
  
  eqconstraint = [];
  df_eqconstraint = [];

  upper       = [15;9];
  lower       = [0; 0];
  x0          = [4; 3]; // Feasible solution
  // x0          = upper; // Infeasible solution
end

if UseProblem==15 then
  //constr_pb_1
  function y = f(x)
    y = 4*x(1) - x(2)^2 - 12;
    global NbFobj;
    NbFobj = NbFobj + 1;
  endfunction
  
  //function y = df(x)
  //  y = derivative(f,x);
  //endfunction
  
  function y = df(x)
    y(1,1) = 4; ...
    y(2,1) = -2*x(2);
  endfunction

  function y = ineqconstraint(x)
    y(1,1) = - 10*x(1) + x(1)^2 - 10*x(2) + x(2)^2 + 34;
  endfunction
  
  //function y = df_ineqconstraint(x)
  //  y = derivative(ineqconstraint,x);
  //endfunction
  
  function y = df_ineqconstraint(x)
    y(1,1) = -10 + 2*x(1);
    y(2,1) = -10 + 2*x(2);
  endfunction

  function y = eqconstraint(x)
    y(1) = 20 - x(1)^2 - x(2)^2;
  endfunction
  
  //function y = df_eqconstraint(x)
  //  y = derivative(eqconstraint,x);
  //endfunction

  function y = df_eqconstraint(x)
    y(1,1) = -2*x(1);
    y(2,1) = -2*x(2);
  endfunction

  upper = [15;15];
  lower = [-15;-15];
  //x0    = [-12;-12]; // Infeasible starting point
  x0    = [1.5;1.5]; // Feasible starting point
end

if UseProblem==16 then
  function y = f(x)
    y = sum(x.^2);
    global NbFobj;
    NbFobj = NbFobj + 1;
  endfunction
  
  function y = df(x)
    y = derivative(f,x);
  endfunction
  
  function y = ineqconstraint(x)
    x = x(:);
    y(1) = -[1 2]*x-[3 4]*x.^2+1;
    y(2) = -[1 4]*x-[4 5]*x.^2+2;
  endfunction
  
  function y = df_ineqconstraint(x)
    y = derivative(ineqconstraint,x);
  endfunction
  
  eqconstraint    = [];
  df_eqconstraint = [];
  upper       = [2;2];
  lower       = [-2; -2];
  x0          = [1.99; 1.99];
  //x0          = [0; 0];
end

if UseProblem==17 then
  //funcname = 'G1'; // A revoir (min or max) 
  //funcname = 'G2'; // A revoir (min or max) 
  //funcname = 'G3'; // A revoir (min or max) 
  //funcname = 'G4'; // A revoir (min or max) 
  //funcname = 'G5'; // A revoir (min or max) 
  //funcname = 'G6'; // A revoir (min or max) 
  //funcname = 'G7'; // A revoir (min or max) 
  //funcname = 'G8'; // A revoir (min or max) 
  //funcname = 'G9'; // A revoir (vx undefined)
  //funcname = 'G10'; // A revoir (min or max) 
  //funcname = 'G11'; // A revoir (min or max) 
  //funcname = 'G12'; // A revoir
  //funcname = 'G13'; // A revoir (min or max) 
  //funcname = 'weldedbeam';
  //funcname = 'pressurevessel';
  //funcname = 'tensioncompr'; // ??
  //funcname = 'vv_open'; // ??
  //funcname = 'himmelblau_1';
  //funcname = 'himmelblau_3';
  //funcname = 'himmelblau_4';
  //funcname = 'himmelblau_4a';
  //funcname = 'himmelblau_5';
  //funcname = 'himmelblau_8';
  //funcname = 'himmelblau_9';
  //funcname = 'himmelblau_10';
  //funcname = 'himmelblau_11';
  //funcname = 'himmelblau_12';
  //funcname = 'himmelblau_13';
  //funcname = 'himmelblau_14';
  //funcname = 'himmelblau_15';
  //funcname = 'himmelblau_16';
  funcname = 'himmelblau_24';

  x0 = eval(funcname+'_x_init()');
  upper = eval('min_bd_'+funcname+'()');
  lower = eval('max_bd_'+funcname+'()');

  function y = f(x)
    y = eval(funcname+'_obj(x)');
    global NbFobj;
    NbFobj = NbFobj + 1;
  endfunction
  
  function y = df(x)
    y = derivative(f,x);
  endfunction
  
  function y = ineqconstraint(x)
    y = eval(funcname+'_ineq(x)');
  endfunction
  
  function y = df_ineqconstraint(x)
    y = derivative(ineqconstraint,x);
  endfunction
  
  function y = eqconstraint(x)
    y = eval(funcname+'_eq(x)');
  endfunction
  
  function y = df_eqconstraint(x)
    y = derivative(eqconstraint,x);
  endfunction
end

if UseProblem==18 then
  //function [value] = elec_on_sph_obj(x)
  //function [diff_bound] = elec_on_sph_constr(x)
  N = 20;
  upper = ones(3*N,1);
  lower = -ones(3*N,1);
  x0 = (upper - lower) .* rand(size(upper,1),size(upper,2)) + lower;

  function y = f(x)
    y = elec_on_sph_obj(matrix(x,N,3));
    global NbFobj;
    NbFobj = NbFobj + 1;
  endfunction
  
  function y = df(x)
    y = derivative(f,x);
  endfunction
  
  function y = eqconstraint(x)
    y = elec_on_sph_eq(matrix(x,N,3));
  endfunction
  
  function y = df_eqconstraint(x)
    y = derivative(eqconstraint,x);
  endfunction
  
  //function y = df_eqconstraint(x)
  //  y = df_elec_on_sph_eq(matrix(x,N,3));
  //endfunction
  
  ineqconstraint = [];
  df_ineqconstraint = [];

  function plot_elec_on_sph(x_opt)
    N = length(x_opt)/3;
    x = matrix(x_opt,N,3);
    param3d1(x(:,1),x(:,2),list(x(:,3),-1));
    xtitle('plot_elec_on_sph,'x','y','z');
  endfunction
end

if UseProblem==19 then
  //function [value] = larg_sm_poly_obj(r, theta)
  //function [diff_bound] = larg_sm_poly_constr(r, theta)
  // 20 variables r, 20 variables theta
  N = 20;
  upper = max_bd_larg_sm_poly(N);
  lower = min_bd_larg_sm_poly(N);
  x0 = larg_sm_poly_x_init(N);

  function y = f(x)
    y = larg_sm_poly_obj(x(1:N),x(N+1:$));
    global NbFobj;
    NbFobj = NbFobj + 1;
  endfunction
  
  function y = df(x)
    y = derivative(f,x);
  endfunction
  
  function y = ineqconstraint(x)
    y = larg_sm_poly_ineq(x(1:N),x(N+1:$));
  endfunction
  
  function y = df_ineqconstraint(x)
    y = derivative(ineqconstraint,x);
  endfunction
  
  eqconstraint    = [];
  df_eqconstraint = [];

  function plot_larg_sm_poly_obj(x_plot)
    x = [];
    y = [];
    N = length(x_plot)/2;
    for i=1:N
      x(i) = x_plot(i)*cos(x_plot(N+i));
      y(i) = x_plot(i)*sin(x_plot(N+i));
    end
    x($+1) = x(1);
    y($+1) = y(1);
    plot(x,y,'k-');
    xtitle('larg_sm_poly_obj solution','x','y');
  endfunction
end

if UseProblem==20 then
  // MMA problem 1
  mma_dim = 1000;
  function y = f(x)
    y = mma_pb_1_obj(x);
    global NbFobj; ...
    NbFobj = NbFobj + 1;
  endfunction
  
  function y = df(x)
    y = mma_pb_1_df_obj(x);
  endfunction
  
  function y = ineqconstraint(x)
    y = mma_pb_1_ineq(x);
  endfunction
  
  function y = df_ineqconstraint(x)
    y = mma_pb_1_df_ineq(x);
  endfunction
  
  eqconstraint = [];
  df_eqconstraint = [];
  upper       = max_bd_mma_pb_1(mma_dim);
  lower       = min_bd_mma_pb_1(mma_dim);
  x0          = mma_pb_1_x_init(mma_dim);
  //x0          = (upper - lower).*rand(size(lower,1),size(lower,2)) + lower;
end

if UseProblem==21 then
  // C - MMA problem 1
  mma_dim = 1000; //10000;
  function y = f(x)  
    [sv_f,sv_df] = sv1_f_df(x);
    global NbFobj;
    NbFobj = NbFobj + 1;
    y = sv_f;
  endfunction
  
  function y = df(x)
    [sv_f,sv_df] = sv1_f_df(x);
    y = sv_df;
  endfunction
  
  function y = ineqconstraint(x)
    [sv_g1,sv_g2,sv_dg1,sv_dg2] = sv1_g_dg(x);
    y(1,1) = sv_g1;
    y(2,1) = sv_g2;
  endfunction
  
  function y = df_ineqconstraint(x)
    [sv_g1,sv_g2,sv_dg1,sv_dg2] = sv1_g_dg(x);
    y(:,1) = sv_dg1;
    y(:,2) = sv_dg2;
  endfunction
  
  eqconstraint = [];
  df_eqconstraint = [];
  [x0, lower, upper] = sv1_init(mma_dim);
end

if UseProblem==22 then
  // MMA problem 2
  mma_dim = 1000;
  function y = f(x)
    y = mma_pb_2_obj(x);
  endfunction
  
  function y = df(x)
    y = mma_pb_2_df_obj(x);
  endfunction
  
  function y = ineqconstraint(x)
    y = mma_pb_2_ineq(x);
  endfunction
  
  function y = df_ineqconstraint(x)
    y = mma_pb_2_df_ineq(x);
  endfunction
  
  eqconstraint = [];
  df_eqconstraint = [];

  upper       = max_bd_mma_pb_2(mma_dim);
  lower       = min_bd_mma_pb_2(mma_dim);
  x0          = mma_pb_2_x_init(mma_dim);
  //x0          = (upper - lower).*rand(size(lower,1),size(lower,2)) + lower;
end

if UseProblem==23 then
  // C - MMA problem 2
  mma_dim = 1000;

  function y = f(x)
    [sv_f,sv_df] = sv1_f_df(x);
    global NbFobj;
    NbFobj = NbFobj + 1;
    y = -sv_f;
  endfunction
  
  function y = df(x)
    [sv_f,sv_df] = sv1_f_df(x);
    y = -sv_df;
  endfunction
  
  function y = ineqconstraint(x)
    [sv_g1,sv_g2,sv_dg1,sv_dg2] = sv1_g_dg(x);
    y(1,1) = -sv_g1;
    y(2,1) = -sv_g2;
  endfunction
  
  function y = df_ineqconstraint(x)
    [sv_g1,sv_g2,sv_dg1,sv_dg2] = sv1_g_dg(x);
    y(:,1) = -sv_dg1;
    y(:,2) = -sv_dg2;
  endfunction

  eqconstraint = [];
  df_eqconstraint = [];
  [x0, lower, upper] = sv2_init(mma_dim);
end

/////////////////////////////////////////////////////////////////////////

// n_dim    NbOuterIter    MaxEvalFunc per iter
// 1000      177            177*209
// 2000      190            190*224
// 5000      221            221*263
// 10000     251            251*296
// 20000     286            286*316

/////////////////////////////////////////////////////////////////////

ItMX     = 200; // 40
Log      = %T;
param    = [];
Plot     = %T;
Plot1    = [];
Plot2    = [];
StepPlot = 1;
DispNum  = %F;

MaxMinStepCount = 1;
MaxEvalFunc     = ItMX;
CTOL            = 0; // if 0 then Desactivation of the KKT tolerance 
ETOL            = 1e-6;
ITOL            = 1e-6;
XTOL            = 1e-6;
STOL            = 0; // if 0 then Desactivation of the KKT stagnation

if (Plot & max(size(x0))==2)     then Plot2 = %T;
elseif (Plot & max(size(x0))==1) then Plot1 = %T;
else
  Plot1 = %F;
  Plot2 = %F;
end

// Parameters for optim_slp
param = init_param();
param = add_param(param,'itol',ITOL);
param = add_param(param,'etol',ETOL);
param = add_param(param,'ctol',CTOL);
param = add_param(param,'xtol',XTOL);
param = add_param(param,'stol',STOL);
param = add_param(param,'maxevalfunc',MaxEvalFunc);
param = add_param(param,'maxminstepcount',MaxMinStepCount);
param = add_param(param,'debug', %F);
param = add_param(param,'movelimitmax',2.0);  // 0.35
param = add_param(param,'movelimit',1.0);  // 0.20
param = add_param(param,'movelimitmin',0.0);  // 0.001
//param = add_param(param,'reducecoeff', []);
//param = add_param(param,'increasecoeff', []);
param = add_param(param,'reducecoeff', 0.8); // 0.5
param = add_param(param,'increasecoeff', 1/0.8); // 0.7
param = add_param(param,'increase_count',1); 
param = add_param(param,'decrease_count',1); 
param = add_param(param,'randomizemovelimit', %F);
param = add_param(param,'randomizemlfactor', 0.1);
param = add_param(param,'offset_ineq',0.0);
//param = add_param(param,'nu',[]);
param = add_param(param,'nu',100); // 100
//param = add_param(param,'restart',ceil(ItMX/2)+5);
//param = add_param(param,'restart_ml',ceil(ItMX/2)+5);
// Parameters for the clp
param = add_param(param,'maxnumiterations',10000);
param = add_param(param,'maxnumseconds',1000);
param = add_param(param,'primaltolerance',1e-10);
param = add_param(param,'dualtolerance',1e-10);
param = add_param(param,'presolve',0);
param = add_param(param,'solver',1);  // 6 interior - 7 pdco - other simplex
//param = add_param(param,'perturb',10);
//param = add_param(param,'optim_dir', -1);
param = add_param(param,'var_type',Integer_Variables);
param = add_param(param,'verbose',1);
param = add_param(param,'clpverbose',1);
param = add_param(param,'cbc_printfrequency', 2);         // YC: 0
param = add_param(param,'cbc_printingmode', 0); // 0 or 1
param = add_param(param,'cbc_dobranchandbound',1); // 0, 1, 2 or 3
param = add_param(param,'cbcmaininit',1);
// linear solver selection
param = add_param(param,'opt_mip_method','cbc'); // 'clp', 'cbc', 'symphony', 'glpk'
param = add_param(param,'opt_lp_method','symphony'); // 'clp', 'cbc', 'symphony', 'glpk'

// Add some CBC parameters for the MINLP problems:
//YC: bug avec ces deux options. Mauvais parametrage. Voir ce qui est fait dans cbcmaininit
//param = init_param_cbc(x0, param);
//param = add_param(param,'cbc_cutoff', 1);

clear x_opt;
clear x_history;
clear ml_history;

[x_opt, x_history, ml_history] = optim_slp(f, df, ...
                                           ineqconstraint, df_ineqconstraint, ...
                                           eqconstraint, df_eqconstraint, ...
                                           x0, ItMX, upper, lower, Log, param);

printf('Initial point\n');
if (ineqconstraint~=[]) then
  printf('Value of the inequality constraints (must be all negatives or null) :'); disp(ineqconstraint(x0)');
end
if (eqconstraint~=[]) then
  printf('Value of the equality constraints (must be all null) :'); disp(eqconstraint(x0)');
end
printf('Value of the objective function : %f\n', f(x0));

printf('Solution:\n');
if (ineqconstraint~=[]) then
  printf('Value of the inequality constraints (must be all negatives or null) :'); disp(ineqconstraint(x_opt)');
end
if (eqconstraint~=[]) then
  printf('Value of the equality constraints (must be all null) :'); disp(eqconstraint(x_opt)');
end
printf('Value of objective function = %f\n', f(x_opt));
//printf('Solution found:'); disp(x_opt');

clear ListF;
clear ListEC;
clear ListIC;
clear X;
clear Y;
clear Z_fobj;
clear Z_ineq;
clear Z_eq;
clear Rect;
clear Color;

///////////////////////////////////////////////////
// Plots:                                        //
// - the objective function = f(Iterations)      //
// - the equality constraints = f(Iterations)    //
// - the inequality constraints = f(Iterations)  //
// - the size of the move limits = f(Iterations) //
///////////////////////////////////////////////////

if (Plot) then
  scf();
  T = 1:length(x_history);
  ListF = [];
  for i=1:length(x_history);
    ListF(i) = f(x_history(i));
  end
  // We plot the objective function
  plot(T,ListF,'k');
  xtitle('Evolution of the objective function','Iteration','F');
  // We plot the equality and inequality constraints function
  if (eqconstraint~=[]) then
    scf();
    ListEC = [];
    for i=1:length(x_history)
      ListEC(:,i) = eqconstraint(x_history(i));
    end
    for i=1:size(ListEC,1)
      plot(1:length(x_history),ListEC(i,:),'k');
    end
    xtitle('Evolution of the equality constraints functions','Iteration','H');
  end
  if (ineqconstraint~=[]) then
    scf();
    ListIC = [];
    for i=1:length(x_history)
      ListIC(:,i) = ineqconstraint(x_history(i));
    end
    for i=1:size(ListIC,1)
      plot(1:length(x_history),ListIC(i,:),'k');
    end
    xtitle('Evolution of the inequality constraints functions','Iteration','H');
  end
  scf();
  ListML = [];
  for i=1:length(ml_history);
    ListML(i) = sqrt(sum(ml_history(i).^2));
  end
  // We plot the size of the move limits
  plot(T,ListML,'k');
  xtitle('Evolution of the move limits','Iteration','ML');
end

////////////////////////////////////////////////////////////////////
// Plot points computed iterations after iterations in a 2D plane //
////////////////////////////////////////////////////////////////////

if (Plot1) then
  scf();
  drawlater;
  x = lower(1):(upper(1) - lower(1))/20:upper(1);
  for i=1:length(x)
    Z(i) = f(x(i));
  end
  plot(x,Z);

  wId = waitbar(0,'Drawing results');
  for i=1:StepPlot:length(x_history)
    if (modulo(i/StepPlot, ceil((length(x_history)/StepPlot) / 10))==0) then
      waitbar(floor(1000*i/length(x_history))/1000,wId);
    end
    plot(x_history(i)(1), f(x_history(i)(1)), 'ro');
    if (i~=length(x_history)) then
      FrameColor = 'g-';
    else
      FrameColor = 'b-';
    end
    plot([x_history(i)(1) + ml_history(i)(1) x_history(i)(1) - ml_history(i)(1)],[f(x_history(i)(1)) f(x_history(i)(1))], FrameColor);
    if (DispNum) then
      xstring(x_history(i)(1) - ml_history(i)(1),f(x_history(i)(1)), string(i));
    end
  end
  xtitle('SLP','x');
  legends(['Move limits','Solution found'],[3,5],1);
  drawnow;
  winclose(wId);
end

//////////////////////////
// Plot the move limits //
// and the level curves //
//////////////////////////

if (Plot2) then
  scf();
  drawlater;
  
  // The level curves of the objective function
  x = lower(1):(upper(1) - lower(1))/20:upper(1);
  y = lower(2):(upper(2) - lower(2))/20:upper(2);
  Z_fobj = [];
  for i=1:size(x,2)
    for j=1:size(y,2)
      Z_fobj(i,j) = f([x(i) y(j)]);
    end
  end
  xset('fpf',' ');
  contour(x,y,Z_fobj, 10);

  // The inequality constraints
  if (ineqconstraint~=[]) then
    Z_ineq = [];
    tmp = ineqconstraint([x(1) y(1)]);
    nb_constr = length(tmp);
    for i=1:size(x,2)
      for j=1:size(y,2)
        Z_ineq(i,j,1:nb_constr) = ineqconstraint([x(i) y(j)])';
      end
    end
    for i=1:nb_constr
      xset('fpf',' ');
      contour2d(x,y,Z_ineq(:,:,i), [0 0], 21);
    end
  end

  // The equality constraints
  if (eqconstraint~=[]) then
    Z_eq = [];
    tmp = eqconstraint([x(1) y(1)]);
    nb_constr = length(tmp);
    for i=1:size(x,2)
      for j=1:size(y,2)
        Z_eq(i,j,1:nb_constr) = eqconstraint([x(i) y(j)])';
      end
    end
    for i=1:nb_constr
      xset('fpf',' ');
      contour2d(x,y,Z_eq(:,:,i), [0 0], 22);
    end
  end

  // The move limits
  for i=1:StepPlot:length(x_history)
    Rect(1,i) = x_history(i)(1) - ml_history(i)(1);
    Rect(2,i) = x_history(i)(2) + ml_history(i)(2);
    Rect(3,i) = 2*ml_history(i)(1);
    Rect(4,i) = 2*ml_history(i)(2);
    X(i) = x_history(i)(1);
    Y(i) = x_history(i)(2);
    Color(i) = -3;
    if i==length(x_history) then Color(i) = -11; end
    if (DispNum) then
      xstring(x_history(i)(1) - ml_history(i)(1), x_history(i)(2) - ml_history(i)(2), string(i));
    end
  end
  xrects(Rect,Color);
  plot(X, Y, 'ro');
  xtitle('SLP','x1','x2');
  legends(['Move limits','Solution found','inequality','equality'],[3,5,21,22],1);
  drawnow;
end

