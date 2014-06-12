function [x_opt, x_history, ml_history] = optim_slp(slp_f, slp_df, slp_g, slp_dg, slp_h, slp_dh, x0, ItMX, upper, lower, Log, param)

// Verification of the parameters
if ~isdef('x0','local') then
  error('optim_slp: error, x0 is mandatory');
else
  if size(x0,1)<=size(x0,2) then
    error('optim_slp: x0 must be a column vector');
  end
end
if ~isdef('upper','local') then
  upper = 1e6*ones(size(x0,1), size(x0,2));
end
if ~isdef('lower','local') then
  lower = -1e6*ones(size(x0,1), size(x0,2));
end
if ~isdef('param','local') then
  param = init_param();
end
if ~isdef('Log','local') then
  Log = %F;
end

if ~isdef('slp_g','local') | ~isdef('slp_dg','local') then
  NoIneqConstr = %T;
else
  NoIneqConstr = (slp_g==[]);
end

if ~isdef('slp_h','local') | ~isdef('slp_dh','local') then
  NoEqConstr = %T;
else
  NoEqConstr = (slp_h==[]);
end

if ~isdef('slp_f','local') then
  error('optim_slp: slp_f is mandatory');
else
  if typeof(slp_f)=='list' then
    deff('y=_slp_f(x)','y=slp_f(1)(x, slp_f(2:$))');
  else
    deff('y=_slp_f(x)','y=slp_f(x)');
  end
end
if ~isdef('slp_df','local') then
  error('optim_slp: slp_df is mandatory');
else
  if typeof(slp_df)=='list' then
    deff('y=_slp_df(x)','y=slp_df(1)(x, slp_df(2:$))');
  else
    deff('y=_slp_df(x)','y=slp_df(x)');
  end
end

if ~NoEqConstr then
  if typeof(slp_h)=='list' then
    deff('y=_slp_h(x)','y=slp_h(1)(x, slp_h(2:$))');
  else
    deff('y=_slp_h(x)','y=slp_h(x)');
  end

  if typeof(slp_dh)=='list' then
    deff('y=_slp_dh(x)','y=slp_dh(1)(x, slp_dh(2:$))');
  else
    deff('y=_slp_dh(x)','y=slp_dh(x)');
  end
end

if ~NoIneqConstr then
  if typeof(slp_g)=='list' then
    deff('y=_slp_g(x)','y=slp_g(1)(x, slp_g(2:$))');
  else
    deff('y=_slp_g(x)','y=slp_g(x)');
  end

  if typeof(slp_dg)=='list' then
    deff('y=_slp_dg(x)','y=slp_dg(1)(x, slp_vdg(2:$))');
  else
    deff('y=_slp_dg(x)','y=slp_dg(x)');
  end
end

f_value  = _slp_f(x0);
df_value = _slp_df(x0);

if ~NoIneqConstr then
  g_value  = _slp_g(x0);
  dg_value = _slp_dg(x0);
else
  g_value  = [];
  dg_value = [];
end

if ~NoEqConstr then
  h_value  = _slp_h(x0);
  dh_value = _slp_dh(x0);
else
  h_value  = [];
  dh_value = [];
end

NbIneqConstr = max(size(g_value));
NbEqConstr   = max(size(h_value));

[nargout, nargin] = argn();

[ETOL, err]                = get_param(param,'etol',1e-6); // threshold of the equality constraints
[ITOL, err]                = get_param(param,'itol',1e-6); // threshold of the inequality constraints
[CTOL, err]                = get_param(param,'ctol',1e-6); // convergence of the norm of the Lagrangian
[XTOL, err]                = get_param(param,'xtol',1e-6); // stagnation of the x_opt
[STOL, err]                = get_param(param,'stol',1e-6); // stagnation of the kkt condition
[MaxEvalFunc, err]         = get_param(param,'maxevalfunc',ItMX*10);
[delta_ml, err]            = get_param(param,'movelimit',0.01); 
[DeltaMLMax, err]          = get_param(param,'movelimitmax',0.1); 
[DeltaMLMin, err]          = get_param(param,'movelimitmin',0); 
[Inc_Iterations, err]      = get_param(param,'increase_count',1); 
[Dec_Iterations, err]      = get_param(param,'decrease_count',1); 
[MaxMinStepCount, err]     = get_param(param,'maxminstepcount',10);
[Debug, err]               = get_param(param,'debug', %F);
[Reduce_coeff, err]        = get_param(param,'reducecoeff', []);
[Increase_coeff, err]      = get_param(param,'increasecoeff', []);
[RandomizeML, err]         = get_param(param,'randomizemovelimit', %F);
[RandomizeFactor, err]     = get_param(param,'randomizemlfactor', 0.1);
[offset_ineq, err]         = get_param(param,'offset_ineq', 0.0); // offset on constraints
[Restart_Iterations, err]  = get_param(param,'restart',ItMX+1);
[Restart_ML, err]          = get_param(param,'restart_ml',ItMX+1);
[Nu, err]                  = get_param(param,'nu',[]);
[var_type, err]            = get_param(param,'var_type',[]);
[write_mps_it, err]        = get_param(param,'write_mps_it',%F);
[restart_ip_on_fault, err] = get_param(param,'restart_ip_on_fault',%F);
[opt_mip_method, err]      = get_param(param,'opt_mip_method','cbc'); // 'cbc', 'symphony', 'glpk', 'lpsolve'
[opt_lp_method, err]       = get_param(param,'opt_lp_method','clp'); // 'clp', 'cbc', 'symphony', 'glpk'

delta_ml_max = DeltaMLMax.*(upper - lower);
if ~isempty(delta_ml) then delta_ml = delta_ml .*(upper - lower); end
delta_ml_min = DeltaMLMin.*(upper - lower);

Inc_Iterations(find(Inc_Iterations<=0)) = 0;
Dec_Iterations(find(Dec_Iterations<=0)) = 0;

if isempty(Nu) then
  tmp = 0;
  if ~NoIneqConstr then
    for i=1:size(dg_value,2)
      tmp = max(tmp, norm(dg_value(:,i)));
    end
  end
  if ~NoEqConstr then
    for i=1:size(dh_value,2)
      tmp = max(tmp, norm(dh_value(:,i)));
    end
  end
  if ~NoIneqConstr | ~NoEqConstr then
    Nu = 10 * norm(df_value) / max(tmp,1);
  else
    Nu = 10;
  end
  
  if Log then
    printf('optim_slp: Nu = %f\n', Nu);
  end
end

x_history_defined = (nargout>=2);
if (x_history_defined) then
  x_history      = list();
  x_history($+1) = x0;
end

ml_history_defined = (nargout>=3);

// dg must return a matrix NbVar x NbConstr
// x0 must be a column vector of size NbVar
// g must be a column vector of size NbConstr
// df must be a column vector of size NbVar
// delta_ml can either be a scalar or a column vector of size NbVar
// lower and upper must be column vectors of size NbVar

NbConstr = NbEqConstr + NbIneqConstr;

if isempty(Reduce_coeff) then
  [tmp, Index] = max(abs(delta_ml_max - delta_ml_min));
  tmp_min = abs(delta_ml_min(Index));
  tmp_max = abs(delta_ml_max(Index));
  if tmp_min<=100*%eps then tmp_min = tmp_max / 1000; end
  Reduce_coeff = exp(3/ItMX * log(tmp_min / tmp_max));
  if Log then
    printf('optim_slp: Reduce_coeff = %f\n', Reduce_coeff);
  end
end

if isempty(Increase_coeff) then
  if Reduce_coeff<0.7 then
    Increase_coeff = 1 / (Reduce_coeff + 0.1);
  else
    Increase_coeff = 1 / Reduce_coeff;
  end
  if Log then
    printf('optim_slp: Increase_coeff = %f\n', Increase_coeff);
  end
end

EqConstrMatr   = [];
EqBoundMatr    = [];
IneqConstrMatr = [];
IneqBoundMatr  = [];

NbVar = max(size(x0)); 

// Move limits parameters
if RandomizeML then
  RndReduceCoeff   = RandomizeFactor*Reduce_coeff;
  RndIncreaseCoeff = RandomizeFactor*Increase_coeff;
else
  RndReduceCoeff   = 0.0;
  RndIncreaseCoeff = 0.0;
end

x_k_m_1 = x0;
x_k_m_2 = x0;
x       = x0;

f_value_old = f_value;
if ~NoIneqConstr then g_value_old = g_value; end
if ~NoEqConstr   then h_value_old = h_value; end

Iteration  = 0;
EvalFunc   = 0;
XTol_Count = 0;

cum_elapsed_time = 0;

// Verification of the size of all the parameters
if (size(x0,1)~=NbVar) then
  error('optim_slp: x0 must be a column vector of size nbvar');
end
if ~NoIneqConstr then
  if size(g_value,1)~=NbIneqConstr then
    error('optim_slp: g must return a column vector of size nbconstr');
  end
  if ((size(dg_value,1)~=NbVar)&(size(dg_value,2)~=NbIneqConstr)) then
    error('optim_slp: dg must return a matrix of size nbvar x nbconstr');
  end
end
if ~NoEqConstr then
  if size(h_value,1)~=NbEqConstr then
    error('optim_slp: h must return a column vector of size nbconstr');
  end
  if ((size(dh_value,1)~=NbVar)&(size(dh_value,2)~=NbEqConstr)) then
    error('optim_slp: dh must return a matrix of size nbvar x nbconstr');
  end
end
if (size(df_value,1)~=NbVar) then
  error('optim_slp: df must return a column vector of size nbvar');
end

if isempty(delta_ml) then
  // Automatic computation of the Move limits
  // From:
  // K. Y. Chan, S. J. Skerlos, P. Papalambros, "An Adaptative Sequential Linear Programming Algorithm for 
  //                                             Optimal Design Problems With Probabilistic Constraints"
  // Trans. Of The ASME, Vol. 129, February 2007
  
  delta_ml = ones(NbVar,1);
  if ~NoIneqConstr then
    tmp = max(abs(g_value.^2) ./ sum(abs(dg_value).^2),'r');
    delta_ml = max([delta_ml tmp*ones(NbVar,1)], 'c');
  end
  if ~NoEqConstr then
    tmp = max(abs(h_value.^2) ./ sum(abs(dh_value).^2),'r');
    delta_ml = max([delta_ml tmp*ones(NbVar,1)], 'c');
  end
elseif (size(delta_ml)~=[1 1]) then
  error('optim_slp: delta_ml must be a scalar');
end
if (size(lower,1)~=NbVar) then
  error('optim_slp: lower must be a column vector of size NbVar');
end
if (size(upper,1)~=NbVar) then
  error('optim_slp: upper must be a column vector of size NbVar');
end

if (ml_history_defined) then
  ml_history      = list();
  ml_history($+1) = delta_ml;
end

Inc_Count = zeros(size(delta_ml,1),size(delta_ml,2));
Dec_Count = zeros(size(delta_ml,1),size(delta_ml,2));
Iteration = 0;
kkt_value     = %inf;
kkt_value_old = 0;
lambda_old    = [];

btype = ascii(ascii('E')*ones(1,NbEqConstr)) + ascii(ascii('L')*ones(1,NbIneqConstr));
if ~isempty(var_type) then
  if length(var_type)~=NbVar then
    error('optim_slp: var_type must have the same size as x0\n');
  end
  var_type_extended = var_type + ascii(ascii('C')*ones(1,2*NbEqConstr+NbIneqConstr));
else
  var_type_extended = ascii(ascii('C')*ones(1,NbVar+2*NbEqConstr+NbIneqConstr));
  var_type          = ascii(ascii('C')*ones(1,NbVar));
end

is_int_var = (ascii(var_type)==ascii('I')) | (ascii(var_type)==ascii('i')) | ...
             (ascii(var_type)==ascii('B')) | (ascii(var_type)==ascii('b'));

is_not_int_var = (ascii(var_type)==ascii('C')) | (ascii(var_type)==ascii('c'));
             
is_mip = or(is_int_var);
Index_not_int = find(is_int_var==%F);
Index_int     = find(is_int_var==%T);

// We set a big move limits for integer variables to allow them ... to vary
delta_ml_max(Index_int) = 2*(upper(Index_int) - lower(Index_int));
delta_ml(Index_int)     = 2*(upper(Index_int) - lower(Index_int));
delta_ml_min(Index_int) = 2*(upper(Index_int) - lower(Index_int));

x_extended_old = [x' zeros(1,NbIneqConstr) zeros(1,2*NbEqConstr)]'; // We add slack variables to the original vector of variables

while((Iteration<ItMX)&(EvalFunc<MaxEvalFunc))
  Iteration = Iteration + 1;
  EvalFunc  = EvalFunc + 1;
  
  // We allow a random restart in the vicinity of the current point
  if (modulo(Iteration,Restart_Iterations)==0) then
    if Log then printf('optim_slp: random restart triggered\n'); end
    x(Index_not_int) = x(Index_not_int) + 2*delta_ml(Index_not_int).*rand(length(Index_not_int),1) - delta_ml(Index_not_int);
    x(Index_not_int) = max([min([x(Index_not_int) upper(Index_not_int)],'c') lower(Index_not_int)],'c');
  end

  // We allow a random restart of the move limits (they are reinitialized which the max move limits)
  if (modulo(Iteration,Restart_ML)==0) then
    if Log then printf('optim_slp: random restart of move limits triggered\n'); end
    delta_ml = delta_ml_max;
  end
  
  x_k_m_2 = x_k_m_1;
  x_k_m_1 = x;

  f_value_old = f_value;
  if ~NoIneqConstr then g_value_old = g_value; end
  if ~NoEqConstr   then h_value_old = h_value; end
  
  if (Iteration>1) & ~isempty(Index_not_int) then
    // Computation of the KKT condition (we use the KKT coefficients given by clp)
    kkt_value_old = kkt_value;
    
    kkt_value = df_value;
    if ~NoIneqConstr then
      kkt_value = kkt_value - sum(extra('lambda')(1:length(g_value)).*.ones(length(x),1) .* dg_value,'c');
    end
    if ~NoEqConstr then
      kkt_value = kkt_value - sum(extra('lambda')(length(g_value)+1:$).*.ones(length(x),1) .* dh_value,'c');
    end
        
    // Computation of the projected gradient of the Lagrangian
    Index = find(abs(x-upper)<=10*%eps*zeros(size(x,1),size(x,2)));
        
    // Sometimes kkt_value was considered as a sparse value and the following test hangs.
    // We can't mix in [] sparse + full
    kkt_value(Index) = max([full(kkt_value(Index)) zeros(length(Index),1)], 'c');
    if Log then
      if length(Index)~=0 then 
        printf('optim_slp: KKT: gradient of the Laplacian has been upper-projected for %d variables\n',length(Index)); 
      end
    end

    Index = find(abs(x(Index_not_int)-lower(Index_not_int))<=10*%eps*zeros(length(Index_not_int),1));
    kkt_value(Index) = min([full(kkt_value(Index)) zeros(length(Index),1)], 'c');
    if Log then
      if length(Index)~=0 then
        printf('optim_slp: KKT: gradient of the Laplacian has been lower-projected for %d variables\n',length(Index));
      end
    end
    
    kkt_value = norm(kkt_value);
    
    if Log then
      printf('optim_slp: KKT condition = %f', kkt_value);
      if ~NoIneqConstr then
        printf(' max(g, 0) = %f', max(max(g_value), 0));
      end
      if ~NoEqConstr then
        printf(' max(abs(h)) = %f', max(abs(h_value)));
      end
      printf('\n');
    end
    
    constr_test = (kkt_value<CTOL);
    if ~NoIneqConstr then
      constr_test = constr_test & (max(g_value, 0)<ITOL);
    end
    if ~NoEqConstr then
      constr_test = constr_test & (max(abs(h_value))<ETOL);
    end

    if (constr_test) then
      if Log then
        printf('optim_slp: break on tolerance reached on KKT conditions\n');
        printf('optim_slp: cumulated elapsed time used  for linear programming: %.2f sec\n', cum_elapsed_time);
      end    
      x_opt = x;
      return;
    end
    
    if (norm(kkt_value - kkt_value_old) / max(1, norm(kkt_value_old)))<STOL then
      if Log then
        printf('optim_slp: break on KKT conditions stagnation\n');
        printf('optim_slp: cumulated elapsed time used  for linear programming: %.2f sec\n', cum_elapsed_time);
      end    
      x_opt = x;
      return;
    end
  end

  // Building local linear models
  // For the objective function
  if (Iteration~=1) then
    f_value  = _slp_f(x);
    df_value = _slp_df(x);
  end
  // For the constraints
  if (Iteration~=1) then
    if ~NoIneqConstr then
      g_value  = _slp_g(x);
      dg_value = _slp_dg(x);
    end
    if ~NoEqConstr then
      h_value  = _slp_h(x);
      dh_value = _slp_dh(x);
    end
  end

  NbIneqConstr = max(size(g_value));
  NbEqConstr   = max(size(h_value));

  NbConstr = NbEqConstr + NbIneqConstr;

  if Log then
    printf('Iteration %d: f = %f\n', Iteration, f_value);
    if ~NoIneqConstr then
      printf('max overshoot for inequality constraints (must be negative or null): %f\n', max(g_value));
    end
    if ~NoEqConstr then
      printf('max overshoot for equality constraints (must be null): %f\n', max(abs(h_value)));
    end
  end
  
  // Initialisation of the slack variables for equality constraints
  // The part related to w_eq_plus
  w_eq_plus  = zeros(size(h_value,1), size(h_value,2));
  Index      = find(h_value>=0);
  w_eq_plus(Index) = abs(h_value(Index));
  // The part related to w_eq_minus
  w_eq_minus = zeros(size(h_value,1), size(h_value,2));
  Index      = find(h_value<0);
  w_eq_minus(Index) = abs(h_value(Index));
  // Initialisation of the slack variables for inequality constraints
  w_ineq = max([g_value(:) zeros(length(g_value), 1)],'c');
 
  // Building the matrix constraint
  // Preparation of the slack variables for equality and inequality constraints
  if ~NoEqConstr then
    w_mat_eq_plus  =   speye(NbEqConstr, NbEqConstr);
    w_mat_eq_minus = - speye(NbEqConstr, NbEqConstr);
  end
  if ~NoIneqConstr then
    w_mat_ineq = - speye(NbIneqConstr, NbIneqConstr);
  end
  // The Inequality constraints
  if ~NoIneqConstr then
    IneqConstrMatr = sparse([dg_value' w_mat_ineq spzeros(NbIneqConstr, NbEqConstr) spzeros(NbIneqConstr, NbEqConstr)]);
    IneqBoundMatr = [];
    for i=1:NbIneqConstr
      IneqBoundMatr(i) = - g_value(i) + dg_value(:,i)'*x;
    end
    if offset_ineq~=0 then
      IneqBoundMatr = (1 + sign(IneqBoundMatr)*offset_ineq).*IneqBoundMatr;
    end
  else
    w_mat_ineq     = [];
    IneqConstrMatr = [];
    IneqBoundMatr  = [];
  end
  // The equality constraints
  if ~NoEqConstr then
    EqConstrMatr = sparse([dh_value' spzeros(NbEqConstr, NbIneqConstr) w_mat_eq_plus w_mat_eq_minus]);
    EqBoundMatr = [];
    for i=1:NbEqConstr
      EqBoundMatr(i) = - h_value(i) + dh_value(:,i)'*x;
    end
  else
    w_mat_eq_plus  = [];
    w_mat_eq_minus = [];
    EqConstrMatr   = [];
    EqBoundMatr    = [];
  end
  // Optimization by linear programming
  x_extended = [x' w_ineq' w_eq_plus' w_eq_minus']'; // We add slack variables to the original vector of variables
  // First part of the lower / upper bounds are related to move limits and
  // second part are related to slack variables for inequality constraints
  // Sometimes, x+delta is below lower !! This formula prevents this phenomenon
  upper_bound = zeros(size(x_extended,1),size(x_extended,2));
  upper_bound(Index_not_int)  = max([min([x(Index_not_int)+delta_ml(Index_not_int) upper(Index_not_int)],'c') lower(Index_not_int)], 'c');
  upper_bound(Index_int)      = upper(Index_int);
  upper_bound(length(x0)+1:$) = %inf;
  // Sometimes, x-delta is above upper !! This formula prevents this phenomenon
  lower_bound = zeros(size(x_extended,1),size(x_extended,2));
  lower_bound(Index_not_int)  = min([max([x(Index_not_int)-delta_ml(Index_not_int) lower(Index_not_int)],'c') upper(Index_not_int)], 'c');
  lower_bound(Index_int)      = lower(Index_int);
  lower_bound(length(x0)+1:$) = 0;

  df_value_extended = [df_value' Nu*ones(1,length(w_ineq)) Nu*ones(1,length(w_eq_plus)) Nu*ones(1,length(w_eq_minus))]';
  
  if (Debug)
    printf('DEBUG: %d / %d\n', Iteration, ItMX);
    if ~NoIneqConstr then
      printf('optim_slp: inequality constraints\n');
      disp(_slp_g(x_extended(1:length(x0))));
      printf('w_ineq\n');
      disp(w_ineq')
    end
    if ~NoEqConstr then
      printf('optim_slp: equality constraints\n');
      disp(_slp_h(x_extended(1:length(x0))));
      printf('w_eq_plus\n');
      disp(w_eq_plus');
      printf('w_eq_minus\n');
      disp(w_eq_minus');
    end
    printf('Nb Eq Constr = %d Nb Ineq Constr = %d\n', NbEqConstr, NbIneqConstr);
  end

  if write_mps_it then
    [res, err] = is_param(param,'writemps');
    if ~res then
      [param, res] = add_param(param,'writemps','optim_slp_'+string(Iteration)+'.mps');
    else
      [param, res] = set_param(param,'writemps','optim_slp_'+string(Iteration)+'.mps');
    end
    
    if Log then
      [txt, err] = get_param(param,'writemps');
      printf('optim_slp: problem to be saved as %s\n', txt);
    end
  end

  t_start = getdate();
  
  if is_mip then
    if opt_mip_method=='cbc' then
      if Log then printf('optim_slp: launching CBC\n'); end
      [x_extended,fmin,status,extra] = cbc([],df_value_extended,[EqConstrMatr' IneqConstrMatr']',[EqBoundMatr' IneqBoundMatr']', ...
                                           [EqBoundMatr' IneqBoundMatr']',lower_bound,upper_bound,btype,var_type_extended,param);
    end

    if opt_mip_method=='symphony' then 
      if Log then printf('optim_slp: launching SYMPHONY\n'); end
      [x_extended,fmin,status,extra] = symphony(df_value_extended,[],[EqConstrMatr' IneqConstrMatr']',[EqBoundMatr' IneqBoundMatr']', ...
                                                [EqBoundMatr' IneqBoundMatr']',lower_bound,upper_bound,btype,var_type_extended,param);
      status = status - 227; // 227 - optimal solution found
    end
                                              
    if opt_mip_method=='glpk' then
      if Log then printf('optim_slp: launching GLPK\n'); end
      [x_extended,fmin,status,extra] = glpk(df_value_extended,[EqConstrMatr' IneqConstrMatr']',[EqBoundMatr' IneqBoundMatr']', ...
                                            [EqBoundMatr' IneqBoundMatr']',lower_bound,upper_bound,btype,var_type_extended,param);
      status = status - 5; // 5 - optimal solution found
    end

    if opt_mip_method=='lpsolve' then
      if Log then printf('optim_slp: launching LPSOLVE\n'); end
      [x_extended,fmin,status,extra] = lpsolve(df_value_extended,[EqConstrMatr' IneqConstrMatr']',[EqBoundMatr' IneqBoundMatr']', ...
                                               [EqBoundMatr' IneqBoundMatr']',lower_bound,upper_bound,btype,var_type_extended,param);
      status = status - 0; // 0 - optimal solution found
    end
  else
    if opt_lp_method=='clp' then 
      if Log then printf('optim_slp: launching CLP\n'); end
      [x_extended,fmin,status,extra] = clp([],df_value_extended,[EqConstrMatr' IneqConstrMatr']',[EqBoundMatr' IneqBoundMatr']', ...
                                           [EqBoundMatr' IneqBoundMatr']',lower_bound,upper_bound,btype,var_type_extended,param);
      if status & restart_ip_on_fault then
        // If the problem was not solvable using simplex, then we restart with interior point. Sometimes it helps.
        [res, err] = is_param(param,'solver');
        old_method = 1;
        if ~res then
          old_method = 1;
          [param, err] = add_param(param,'solver',6);
        else
          [old_method, err] = get_param(param,'solver');
          [param, err] = set_param(param,'solver',6);
        end

        if Log then printf('optim_slp: status is not null - restart with interior point\n'); end

        [x_extended,fmin,status,extra] = clp([],df_value_extended,[EqConstrMatr' IneqConstrMatr']',[EqBoundMatr' IneqBoundMatr']', ...
                                             [EqBoundMatr' IneqBoundMatr']',lower_bound,upper_bound,btype,var_type_extended,param);
                                           
        // Once solved, switch back to the current method
        [param, err] = set_param(param,'solver',old_method);
      end
    end 
  
    if opt_lp_method=='symphony' then 
      if Log then printf('optim_slp: launching SYMPHONY\n'); end
      [x_extended,fmin,status,extra] = symphony(df_value_extended,[],[EqConstrMatr' IneqConstrMatr']',[EqBoundMatr' IneqBoundMatr']', ...
                                                [EqBoundMatr' IneqBoundMatr']',lower_bound,upper_bound,btype,var_type_extended,param);
      status = status - 227; // 227 - optimal solution found
    end

    if opt_lp_method=='glpk' then 
      if Log then printf('optim_slp: launching GLPK\n'); end
      [x_extended,fmin,status,extra] = glpk(df_value_extended,[EqConstrMatr' IneqConstrMatr']',[EqBoundMatr' IneqBoundMatr']', ...
                                            [EqBoundMatr' IneqBoundMatr']',lower_bound,upper_bound,btype,var_type_extended,param);

      status = status - 5; // 5 - optimal solution found
    end

    if opt_lp_method=='lpsolve' then
      if Log then printf('optim_slp: launching LPSOLVE\n'); end
      [x_extended,fmin,status,extra] = lpsolve(df_value_extended,[EqConstrMatr' IneqConstrMatr']',[EqBoundMatr' IneqBoundMatr']', ...
                                               [EqBoundMatr' IneqBoundMatr']',lower_bound,upper_bound,btype,var_type_extended,param);
      status = status - 0; // 0 - optimal solution found
    end
  end

  if (Debug) then
    printf('optim_slp: x_extended after\n');
    disp(x_extended)
  end

  x          = x_extended(1:length(x0));
  w_ineq     = x_extended(length(x0)+1:length(x0)+1+NbIneqConstr-1);
  w_eq_plus  = x_extended(length(x0)+1+NbIneqConstr:length(x0)+1+NbIneqConstr+NbEqConstr-1);
  w_eq_minus = x_extended(length(x0)+1+NbIneqConstr+NbEqConstr:length(x0)+1+NbIneqConstr+NbEqConstr+NbEqConstr-1);

  lambda_old = extra('lambda');
  
  if Log then
    printf('results of the linear program\n');
    if (opt_mip_method=='symphony') | (opt_lp_method=='symphony') then 
      printf('status = %d\n', status + 227); // 227 - optimal solution found
    end
    if (opt_mip_method=='glpk') | (opt_lp_method=='glpk') then 
      printf('status = %d\n', status + 5); // 5 - optimal solution found
    end
    
    printf('time needed to solve linear program: %.2f sec\n',etime(getdate(),t_start));
    cum_elapsed_time = cum_elapsed_time + etime(getdate(),t_start);
    
    printf('Iteration %d: fmin = %f\n', Iteration, fmin);
    if ~NoIneqConstr then
      tmp = IneqConstrMatr * x_extended_old - IneqBoundMatr;
      printf('max overshoot for the linear model of inequality constraints (must be negative or null): %f\n', max(tmp));
      printf('max value of the slack variable for the linear model of inequality constraints (must be positive or null): %f\n', max(w_ineq));
    end
    if ~NoEqConstr then
      tmp = EqConstrMatr * x_extended_old - EqBoundMatr;
      printf('max overshoot for the linear model of equality constraints (must be null): %f\n', max(abs(tmp)));
      printf('max value of the positive slack variable for the linear model of equality constraints (must be null): %f\n', max(w_eq_plus));
      printf('max value of the negative slack variable for the linear model of equality constraints (must be null): %f\n', max(w_eq_minus));
    end
  end

  if (x_history_defined) then
    x_history($+1) = x;
  end
  
  // Update move limits
  if (Debug) then
    printf('optim_slp: delta_ml before'); disp(delta_ml')
  end
  // If the problem was not solvable, we increase the move limits
  
  if status then
    delta_ml  = delta_ml .* (Increase_coeff + (2*rand(1,1)-1)*RndReduceCoeff);
    delta_ml  = max([min([delta_ml delta_ml_max],'c') delta_ml_min],'c');
    Inc_Count = zeros(size(delta_ml,1),size(delta_ml,2));
    Dec_Count = zeros(size(delta_ml,1),size(delta_ml,2));
  else
    Index_Dec = find(((x(Index_not_int) - x_k_m_1(Index_not_int)).*(x_k_m_1(Index_not_int) - x_k_m_2(Index_not_int))<=0));
    Dec_Count(Index_Dec) = Dec_Count(Index_Dec) + 1;
    Index_Dec = [];

    Index_Inc = find(((x(Index_not_int) - x_k_m_1(Index_not_int)).*(x_k_m_1(Index_not_int) - x_k_m_2(Index_not_int))>0));
    Inc_Count(Index_Inc) = Inc_Count(Index_Inc) + 1;
    Index_Inc = [];

    Index = find(Inc_Count>=Inc_Iterations);
    if ~isempty(Index) then
      delta_ml(Index) = delta_ml(Index) .* (Increase_coeff*ones(length(Index),1) + (2*rand(length(Index),1)-1)*RndIncreaseCoeff);
      Inc_Count(Index) = 0;
    end

    Index = find(Dec_Count>=Dec_Iterations);
    if ~isempty(Index) then
      delta_ml(Index) = delta_ml(Index) .* (Reduce_coeff*ones(length(Index),1) + (2*rand(length(Index),1)-1)*RndReduceCoeff);
      Dec_Count(Index) = 0;
    end
    Index = [];
  end

  x_extended = x_extended_old;

  delta_ml(Index_not_int) = min([max([delta_ml(Index_not_int) delta_ml_min(Index_not_int)],'c') delta_ml_max(Index_not_int)],'c');
  
  if ~isempty(Index_not_int) then 
    printf('delta_ml: min = %f max = %f inc = %f dec = %f\n', min(delta_ml(Index_not_int)), max(delta_ml(Index_not_int)), Increase_coeff, Reduce_coeff);
  end

  if (Debug) then
    printf('optim_slp: delta_ml after'); disp(delta_ml')
  end

  if (ml_history_defined) then
    ml_history($+1) = delta_ml;
  end

  // If the x step is too small, we stop
  if (norm(x - x_k_m_1)<XTOL) then 
    XTol_Count = XTol_Count + 1;
    if Debug then
      printf('optim_slp: norm(x - xk) = %f XTol = %f XTol_Count = %d\n', norm(x - x_k_m_1), XTOL, XTol_Count);
      printf('optim_slp: cumulated elapsed time used  for linear programming: %.2f sec\n', cum_elapsed_time);
    end    
  else
    XTol_Count = 0;
  end

  if (XTol_Count==MaxMinStepCount)&(Iteration>1) then 
    if (Log) then
      printf('optim_slp: stop on x tolerance count\n');
      printf('optim_slp: cumulated elapsed time used  for linear programming: %.2f sec\n', cum_elapsed_time);
    end
    x_opt = x;
    return;
  end
end

if Log then
  printf('optim_slp: cumulated elapsed time used  for linear programming: %.2f sec\n', cum_elapsed_time);
end

x_opt = x;
endfunction
