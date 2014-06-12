function [x_sol, f_sol, extra] = bonmin(x0, _f, _df, _g, _dg, sparse_dg, _dh, sparse_dh, ...
    	                                var_type, var_lin_type, constr_lin_type, ...
					constr_rhs, constr_lhs, ...
					lower, upper, ...
					params)

if ~isdef('x0','local') | ~isdef('_f','local') | ~isdef('_df','local') | ...
   ~isdef('_g','local') | ~isdef('_dg','local') | ...
   ~isdef('sparse_dg','local') then 
  error('bonmin: x0, f, df, g, dg, sparse_dg are mandatory parameters');
end

if isempty(x0) then
  error('bonmin: x0 must not be empty');
end

if ~isdef('_f','local') then
  error('ipopt: f is mandatory');
end

if ~isdef('_df','local') then
  error('ipopt: df is mandatory');
end

if typeof(_f)=='list' then
  deff('y=__f(x,x_new)','y=_f(1)(x,x_new,_f(2:$));');
else
  __f = _f;
end

if typeof(_df)=='list' then
  deff('y=__df(x,x_new)','y=_df(1)(x,x_new,_df(2:$));');
else
  __df = _df;
end
  
if typeof(_g)=='list' then
  deff('y=__g(x,x_new)','y=_g(1)(x,x_new,_g(2:$));');
else
  __g = _g;
end

if typeof(_dg)=='list' then
  deff('y=__dg(x,x_new)','y=_dg(1)(x,x_new,_dg(2:$));');
else
  __dg = _dg;
end

if typeof(_dh)=='list' then
  deff('y=__dh(x,lambda,obj_weight,x_new,lambda_new)','y=_dh(1)(x,lambda,obj_weight,x_new,lambda_new,_dh(2:$));');
else
  __dh = _dh;
end

// Is 0 == CONTINUOUS ?
if ~isdef('var_type','local') then
  var_type = zeros(size(x0,1),size(x0,2));
end
if isempty(var_type) then
  var_type = zeros(size(x0,1),size(x0,2));
end

if ~isdef('var_lin_type','local') then
  var_lin_type = zeros(size(x0,1),size(x0,2));
end
if isempty(var_lin_type) then
  var_lin_type = zeros(size(x0,1),size(x0,2));
end

if ~isdef('constr_lin_type','local') then
  constr_lin_type = ones(length(constr_rhs),1);
end
if isempty(constr_lin_type) then
  constr_lin_type = ones(length(constr_rhs),1);
end

if ~isdef('constr_rhs','local') then
  constr_rhs = -%inf*ones(length(constr_rhs),1);
end
//if isempty(constr_rhs) then
//  error('bonmin: constr_rhs must not be empty');
//end

if ~isdef('constr_lhs','local') then
  constr_lhs = %inf*ones(length(constr_rhs),1);
end
//if isempty(constr_lhs) then
//  error('bonmin: constr_lhs must not be empty');
//end

if ~isdef('lower','local') then
  lower = -%inf*ones(size(x0,1),size(x0,2));
end
if isempty(lower) then
  lower = -%inf*ones(size(x0,1),size(x0,2));
end

if ~isdef('upper','local') then
  upper = %inf*ones(size(x0,1),size(x0,2));
end
if isempty(upper) then
  upper = %inf*ones(size(x0,1),size(x0,2));
end

if ~isdef('params','local') then
  params = init_param();
end

[x_sol, f_sol, extra] = scibonmin(x0, __f, __df, __g, __dg, sparse_dg, __dh, sparse_dh, ...
	                          var_type, var_lin_type, constr_lin_type, ...
				  constr_rhs, constr_lhs, ...
				  lower, upper, ...
				  params);
endfunction
