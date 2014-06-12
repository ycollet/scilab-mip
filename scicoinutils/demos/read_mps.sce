mps_filename = 'data/misc07.mps';
mps_type = 0; // 0 -> mps
              // 1 -> gms (GAMS)
              // 2 -> gmpl (AMPL)
Log = 1;

mps_file    = read_mps_file(mps_filename, mps_type, Log);

printf('parameters of problem %s\n', mps_filename);
printf('number of constraints = %d\n',mps_file('nb_constr'));
printf('number of variables   = %d\n',mps_file('nb_obj_var'));
printf('number of coefficients in the constraint matrix = %d\n', mps_file('nb_val_constr_mat'));
printf('number of integer variables = %d\n', mps_file('nb_int_var'));

mps_file_mp = read_mps_file_mp(mps_filename, mps_type, Log);

var_type = string(zeros(1,length(mps_file_mp('obj_var_is_int'))));
var_type(find(mps_file_mp('obj_var_is_int')==0)) = 'C';
var_type(find(mps_file_mp('obj_var_is_int')==1)) = 'I';
var_type = strcat(var_type);

printf('problem name: %s\n', mps_file_mp('pb_name'));
printf('size of constraint matrix      = [%d %d]\n', size(mps_file_mp('constr_mat'),1), size(mps_file_mp('constr_mat'),2));
printf('constraints sense              = %s\n', mps_file_mp('constr_sense'));
printf('size of objective coefficients = %d\n', length(mps_file_mp('obj_coeff')));
printf('size of rhs            = %d\n', length(mps_file_mp('rhs')));
printf('size of lhs            = %d\n', length(mps_file_mp('lhs')));
printf('size of bounds_lower   = %d\n', length(mps_file_mp('bounds_lower')));
printf('size of bounds_upper   = %d\n', length(mps_file_mp('bounds_upper')));
printf('size of obj_var_is_int = %d\n', length(mps_file_mp('obj_var_is_int')));
printf('variable type          = %s\n', var_type);
printf('column name:'); disp(mps_file_mp('col_name'));
printf('row name:');    disp(mps_file_mp('row_name'));

// Input and output parameters of write_mps_file
// C_IN        1   d
// A_IN        2   d
// LHS_IN      3   d
// RHS_IN      4   d
// UB_IN       5   d
// LB_IN       6   d
// BTYPE_IN    7   c
// VARTYPE_IN  8   c
// PB_NAME_IN  9   c
// COL_NAME_IN 10  S
// ROW_NAME_IN 11  S
// FILENAME_IN 12  c
// VERBOSE_IN  13  i
// STATUS_OUT  14  i

printf('Writing MPS file test.mps\n');

status = write_mps_file(mps_file_mp('obj_coeff'), ...
                        mps_file_mp('constr_mat'),...
                        mps_file_mp('lhs'), ...
                        mps_file_mp('rhs'), ...
                        mps_file_mp('bounds_upper'), ...
                        mps_file_mp('bounds_lower'), ...
                        mps_file_mp('constr_sense'), ...
                        var_type, ...
                        mps_file_mp('pb_name'), ...
                        mps_file_mp('col_name'), ...
                        mps_file_mp('row_name'), ...
                        'test.mps', ...
                        0);
                        
printf('Writing LP file test.mps\n');

status = write_lp_file(mps_file_mp('obj_coeff'), ...
                       mps_file_mp('constr_mat'),...
                       mps_file_mp('lhs'), ...
                       mps_file_mp('rhs'), ...
                       mps_file_mp('bounds_upper'), ...
                       mps_file_mp('bounds_lower'), ...
                       mps_file_mp('constr_sense'), ...
                       var_type, ...
                       mps_file_mp('pb_name'), ...
                       mps_file_mp('col_name'), ...
                       mps_file_mp('row_name'), ...
                       'test.lp');

printf('End of demo\n');