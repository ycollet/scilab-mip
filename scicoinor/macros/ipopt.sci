// Copyright (C) 2010 - DIGITEO - Yann Collette
// Copyright (C) 2010 - DIGITEO - Michael Baudin
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

//
// Calling sequence
//   Only params is optionnal.
//   [x_sol, f_sol, extra] = ipopt (x0, _f, _df, _g, _dg, sparse_dg, _dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper)
//   [x_sol, f_sol, extra] = ipopt (x0, _f, _df, _g, _dg, sparse_dg, _dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params )
// TODO : make clearer the input arguments which are optionnal or mandatory
// TODO : add unit tests which correspond to each case
//

function [x_sol, f_sol, extra] = ipopt (x0, _f, _df, _g, _dg, sparse_dg, _dh, sparse_dh, ...
  var_lin_type, constr_lin_type, ...
  constr_rhs, constr_lhs, ...
  lower, upper, params)
  
  [lhs,rhs]=argn();
  if ( rhs <> 14 & rhs <> 15 ) then
    errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while 14 or 15 are expected."), "ipopt", rhs);
    error(errmsg)
  end

  
  if ~isdef('x0','local') | ~isdef('_f','local') | ~isdef('_df','local') then 
    error('ipopt: x0, f, df are mandatory parameters');
  end
  
  if isempty(x0) then
    error('ipopt: x0 must not be empty');
  end
  
  if ~isdef('_f','local') then
    error('ipopt: f is mandatory');
  end
  
  if ~isdef('_df','local') then
    error('ipopt: df is mandatory');
  else
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
  //  error('ipopt: constr_rhs must not be empty');
  //end
  if ~isdef('constr_lhs','local') then
    constr_lhs = %inf*ones(length(constr_rhs),1);
  end
  //if isempty(constr_lhs) then
  //  error('ipopt: constr_lhs must not be empty');
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

  [int_cb, _err] = get_param(params,'intermediate_callback',[]);

  apifun_typereal ( x0 , "x0" , 1 )
  apifun_typecallable ( _f , "f" , 2 )
  apifun_typecallable ( _df , "df" , 3 )
  if _g~=[] then
    apifun_typecallable ( _g , "g" , 4 )
    if _dg~=[] then
      apifun_typecallable ( _dg , "dg" , 5 )
    end
    apifun_typereal ( sparse_dg , "sparse_dg" , 6 )
  end
  if _dh~=[] then
    apifun_typecallable ( _dh , "dh" , 7 )
    apifun_typereal ( sparse_dh , "sparse_dh" , 8 )
  end
  apifun_typereal ( var_lin_type , "var_lin_type" , 9 )
  apifun_typereal ( constr_lin_type , "constr_lin_type" , 10 )
  apifun_typereal ( constr_rhs , "constr_rhs" , 11 )
  apifun_typereal ( constr_lhs , "constr_lhs" , 12 )
  apifun_typereal ( lower , "lower" , 13 )
  apifun_typereal ( upper , "upper" , 14 )
  if int_cb~=[] then
    apifun_typecallable ( int_cb , "int_cb" , 15 )
  end
  apifun_typeplist ( params , "params" , 16 )
  
  [x_sol, f_sol, extra] = sciipopt(x0, __f, __df, __g, __dg, sparse_dg, __dh, sparse_dh, ...
                                   var_lin_type, constr_lin_type, ...
                                   constr_rhs, constr_lhs, ...
                                   lower, upper, int_cb, ...
                                   params);
endfunction

////////////////////////
// Internal functions //
////////////////////////

// Generates an error if the given variable is not of type boolean
function apifun_typeboolean ( var , varname , ivar )
  if ( type ( var ) <> 4 ) then
    errmsg = msprintf(gettext("%s: Expected boolean but for variable %s at input #%d, got %s instead."),"apifun_typeboolean", varname , ivar , typeof(var) );
    error(errmsg);
  end
endfunction
// Generates an error if the given variable is not of type real
function apifun_typereal ( var , varname , ivar )
  if ( type ( var ) <> 1 ) then
    errmsg = msprintf(gettext("%s: Expected real variable for variable %s at input #%d, but got %s instead."),"apifun_typereal", varname , ivar , typeof(var) );
    error(errmsg);
  end
endfunction

// Generates an error if the given variable is not of type boolean
function apifun_typestring ( var , varname , ivar )
  if ( type ( var ) <> 10 ) then
    errmsg = msprintf(gettext("%s: Expected string but for variable %s at input #%d, got %s instead."),"apifun_typestring", varname , ivar , typeof(var) );
    error(errmsg);
  end
endfunction
// Generates an error if the given variable is not of type hypermat
function apifun_typehypermat ( var , varname , ivar )
  if ( typeof ( var ) <> "hypermat" ) then
    errmsg = msprintf(gettext("%s: Expected hypermat variable for variable %s at input #%d, but got %s instead."),"apifun_typehypermat", varname , ivar , typeof(var) );
    error(errmsg);
  end
endfunction
// Generates an error if the given variable is not of type cell
function apifun_typecell ( var , varname , ivar )
  if ( typeof ( var ) <> "ce" ) then
    errmsg = msprintf(gettext("%s: Expected cell variable for variable %s at input #%d, but got %s instead."),"apifun_typecell", varname , ivar , typeof(var) );
    error(errmsg);
  end
endfunction

// Generates an error if the given variable is not of type "callable function"
function apifun_typecallable ( var , varname , ivar )
  if ( ~( or ( type ( var ) == [10 11 13 14 15] ) ) ) then
    errmsg = msprintf(gettext("%s: Expected string, macro or list but for variable %s at input #%d, got %s instead."),"apifun_typecallable", varname , ivar , typeof(var) );
    error(errmsg);
  end
endfunction
// Generates an error if the given variable is not of type plist
function apifun_typeplist ( var , varname , ivar )
  if ( typeof ( var ) <> "plist" ) then
    errmsg = msprintf(gettext("%s: Expected plist but for variable %s at input #%d, got %s instead."),"apifun_typeplist", varname , ivar , typeof(var) );
    error(errmsg);
  end
endfunction

