gms_filename = 'data/alan.gms';

gms_type = 1; // 0 -> mps
              // 1 -> gms (GAMS)
              // 2 -> gmpl (AMPL)
Log = 1;

gms_file = read_mps_file(gms_filename, gms_type, Log);

printf('parameters of problem %s\n', gms_filename);
printf('number of constraints = %d\n',gms_file('nb_constr'));
printf('number of variables   = %d\n',gms_file('nb_obj_var'));
printf('number of coefficients in the constraint matrix = %d\n', gms_file('nb_val_constr_mat'));
printf('number of integer variables = %d\n', gms_file('nb_int_var'));

gms_file_mp = read_mps_file_mp(gms_filename, gms_type, Log);

var_type = string(zeros(1,length(gms_file_mp('obj_var_is_int'))));
var_type(find(gms_file_mp('obj_var_is_int')==0)) = 'C';
var_type(find(gms_file_mp('obj_var_is_int')==1)) = 'I';
var_type = strcat(var_type);

printf('problem name: %s\n', gms_file_mp('pb_name'));
printf('size of constraint matrix      = [%d %d]\n', size(gms_file_mp('constr_mat'),1), size(gms_file_mp('constr_mat'),2));
printf('constraints sense              = %s\n', gms_file_mp('constr_sense'));
printf('size of objective coefficients = %d\n', length(gms_file_mp('obj_coeff')));
printf('size of rhs            = %d\n', length(gms_file_mp('rhs')));
printf('size of lhs            = %d\n', length(gms_file_mp('lhs')));
printf('size of bounds_lower   = %d\n', length(gms_file_mp('bounds_lower')));
printf('size of bounds_upper   = %d\n', length(gms_file_mp('bounds_upper')));
printf('size of obj_var_is_int = %d\n', length(gms_file_mp('obj_var_is_int')));
printf('variable type          = %s\n', var_type);
printf('column name:'); disp(gms_file_mp('col_name'));
printf('row name:');    disp(gms_file_mp('row_name'));

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

status = write_mps_file(gms_file_mp('obj_coeff'), ...
                        gms_file_mp('constr_mat'),...
                        gms_file_mp('lhs'), ...
                        gms_file_mp('rhs'), ...
                        gms_file_mp('bounds_upper'), ...
                        gms_file_mp('bounds_lower'), ...
                        gms_file_mp('constr_sense'), ...
                        var_type, ...
                        gms_file_mp('pb_name'), ...
                        gms_file_mp('col_name'), ...
                        gms_file_mp('row_name'), ...
                        'test.mps', ...
                        1);
