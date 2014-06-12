lines(0);

stacksize('max');

// Define a list of test files:
// - mps_miplib
// - mps_lplib
// - mps_mitt_c
// - mps_bi_sc
// - mps_Coral
// - mps_coinor_sample
// - mps_coinor_infeas
// - mps_coinor_big
// - mps_coinor_netlib
// - mps_coinor_miplib3
// - mps_coinor_data

path_new = get_absolute_file_path('mps_glpk_test.sce');
path_old = pwd();

cd(path_new)

exec('coinor_data.sce');

/////////////////////////////////
// Choose the problem to solve //
/////////////////////////////////

// See listoffiles.html for informations related to these problems

Solve_miplib   = %F; // 92 files
Solve_lplib    = %F; // 127 files
Solve_mitt_c   = %F; // 10 files
Solve_bi_sc    = %F; // 21 files
Solve_Coral    = %F; // 372 files
Solve_miplib3  = %F; // 164 files
Solve_sample   = %F; // 17 files
Solve_infeas   = %F; // 29 files
Solve_big      = %F; // 1 file
Solve_netlib   = %F; // 90 files
Solve_data     = %T; // 2 files
Solve_OptimSLP = %F; // 3 files

///////////////////////////
// Set global paremeters //
///////////////////////////

Log       = 0;
Debug     = 0;
mps_index = 2; // 1 LP - 2 MIP // Pb avec 1 et 8 de lplib pb with pp08aCUTs.mps
// pb with lplib + 3 et 4
// pb de nan avec sample + 4 + simplex et sans relaxation
// pb de nan dans le status avec sample + 14 avec relaxation
mps_type  = 0;
UseLinpro = %F;
PrintSol  = %F;

if Solve_miplib   then mps_filename = mps_miplib(mps_index);          end
if Solve_lplib    then mps_filename = mps_lplib(mps_index);           end
if Solve_mitt_c   then mps_filename = mps_mitt_c(mps_index);          end
if Solve_bi_sc    then mps_filename = mps_bi_sc(mps_index);           end
if Solve_Coral    then mps_filename = mps_Coral(mps_index);           end
if Solve_miplib3  then mps_filename = mps_coinor_miplib3(mps_index);  end
if Solve_sample   then mps_filename = mps_coinor_sample(mps_index);   end
if Solve_infeas   then mps_filename = mps_coinor_infeas(mps_index);   end
if Solve_big      then mps_filename = mps_coinor_big(mps_index);      end
if Solve_netlib   then mps_filename = mps_coinor_netlib(mps_index);   end
if Solve_data     then mps_filename = mps_coinor_data(mps_index);     end
if Solve_OptimSLP then mps_filename = mps_optimslp_sample(mps_index); end

//  INT: msglev   - must be 0 (no output [default]) or 1 (error messages only) or 2 (normal output) or 3 (full output)
//  INT: scale    - must be 0 (no scaling) or 1 (equilibration scaling [default]) or 2 (geometric mean scaling)
//  INT: dual     - must be 0 (do NOT use dual simplex [default]) or 1 (use dual simplex)
//  INT: price    - must be 0 (textbook pricing) or 1 (steepest edge pricing [default])
//  INT: round    - must be 0 (report all primal and dual values [default]) or 1 (replace tiny primal and dual values by exact zero)");
//  INT: itlim    - Simplex iterations limit
//  INT: outfrq   - Output frequency, in iterations
//  INT: branch   - must be (MIP only) 0 (branch on first variable) or 1 (branch on last variable) 
//                  or 2 (branch using a heuristic by Driebeck and Tomlin [default]
//  INT: btrack   - must be (MIP only) 0 (depth first search) or 1 (breadth first search) 
//                  or 2 (backtrack using the best projection heuristic [default]
//  INT: presol   - must be 0 (do NOT use LP presolver) or 1 (use LP presolver [default])
//  INT: usecuts  - must be 0 (do NOT generate cuts) or 1 (generate Gomory's cuts [default])");
//  INT: lpsolver - must be 1 (simplex method) or 2 (interior point method)");
//  REAL: relax   - Ratio test option
//  REAL: tolbnd  - Relative tolerance used to check if the current basic solution is primal feasible
//  REAL: toldj   - Absolute tolerance used to check if the current basic solution is dual feasible
//  REAL: tolpiv  - Relative tolerance used to choose eligible pivotal elements of the simplex table in the ratio test
//  REAL: objll
//  REAL: objul
//  REAL: tmlim
//  REAL: outdly
//  REAL: tolint
//  REAL: tolobj

param = init_param();
param = add_param(param,'msglev',   3); // msglev: max 3
param = add_param(param,'scale',    0);
param = add_param(param,'dual',     0);
param = add_param(param,'price',    0);
param = add_param(param,'round',    0);
param = add_param(param,'itlim',    100000); // In iterations
param = add_param(param,'outfrq',   1000);   // In iterations
param = add_param(param,'branch',   2);
param = add_param(param,'btrack',   2);
param = add_param(param,'presol',   0);
param = add_param(param,'usecuts',  1);
param = add_param(param,'lpsolver', 1);   // 1: simplex 0: interior 2: exact
param = add_param(param,'basis_type', 0); // 0: standard, 1: advanced, 2: bixby
param = add_param(param,'scale_flag', 0); //
param = add_param(param,'mip_presolve', 1); // 0: presolve method of glp_intopt, 1: simplex, 2: interior
param = add_param(param,'solve_relaxed', 0);
// branching technique:
// GLP_BR_FFV 1 - first fractional variable
// GLP_BR_LFV 2 - last fractional variable
// GLP_BR_MFV 3 - most fractional variable
// GLP_BR_DTH 4 - heuristic by Driebeck and Tomlin
param = add_param(param,'br_tech', 4);
// backtracking technique:
// GLP_BT_DFS 1 - depth first search
// GLP_BT_BFS 2 - breadth first search
// GLP_BT_BLB 3 - best local bound
// GLP_BT_BPH 4 - best projection heuristic
param = add_param(param,'bt_tech', 4);
// preprocessing technique:
// GLP_PP_NONE 0 - disable preprocessing
// GLP_PP_ROOT 1 - preprocessing only on root level
// GLP_PP_ALL  2 - preprocessing on all levels
param = add_param(param,'pp_tech', 1); // It's effective to activate the preprocessing (2 is stronger than 1 but 1 is a good tradeoff)
// relative MIP gap tolerance
param = add_param(param,'mip_gap', 1e-6);
// MIR cuts (GLP_ON/GLP_OFF)
param = add_param(param,'mir_cuts',1);
// GOMORY cuts (GLP_ON/GLP_OFF)
param = add_param(param,'gmi_cuts', 1);
// COVER cuts (GLP_ON/GLP_OFF)
param = add_param(param,'cov_cuts', 0);
// CLIQUE cuts (GLP_ON/GLP_OFF)
param = add_param(param,'clq_cuts', 0);
// try to binarize integer variables 
param = add_param(param,'binarize',0);
//param = add_param(param,'relax',    1);
//param = add_param(param,'tolbnd',   1);
//param = add_param(param,'toldj',    1);
//param = add_param(param,'tolpiv',   1);
//param = add_param(param,'objll',    1);
//param = add_param(param,'objul',    1);
param = add_param(param,'tmlim',    10000*1000); // In milliseconds
//param = add_param(param,'outdly',   1);
param = add_param(param,'tolint',   1e-5);
param = add_param(param,'tolobj',   1e-6);
//param = add_param(param,'sense',   1);
param = add_param(param,'writemps','test.mps');

///////////////////////////////
// Gunzip the chosen problem //
///////////////////////////////

is_gzipped = %F;
if ~isempty(grep(mps_filename,'.gz')) then
  if getos() == 'Windows' then
    unix_w('.' + filesep() + 'data' + filesep() + 'gzip' + filesep() + 'gzip.exe -d ' + mps_filename);
  else
    unix('gunzip ' + mps_filename);
  end
  mps_filename = strsubst(mps_filename,'.gz','');
  is_gzipped = %T;
end

///////////////////
// Read the data //
///////////////////

t_start = getdate();

printf('Reading informations of file |%s|\n', mps_filename);
mps_file = read_mps_file(mps_filename, mps_type);

printf('Reading content of file |%s|\n', mps_filename);
mps_file_mp = read_mps_file_mp(mps_filename, mps_type);

t_end = getdate();

/////////////////////////////
// Gzip the chosen problem //
/////////////////////////////

if is_gzipped then
  if getos() == 'Windows' then
    unix_w('.' + filesep() + 'data' + filesep() + 'gzip' + filesep() + 'gzip.exe ' + mps_filename);
  else
    unix('gzip ' + mps_filename);
  end
end

printf('elapsed time for reading file = %f secondes\n', etime(t_end, t_start));

t_start = t_end;

///////////////////////////
// Set the variable type //
///////////////////////////

// 'I' -> integer

var_type = string(zeros(1,length(mps_file_mp('obj_var_is_int'))));
var_type(find(mps_file_mp('obj_var_is_int')==0)) = 'C';
var_type(find(mps_file_mp('obj_var_is_int')==1)) = 'I';
var_type = strcat(var_type);

/////////////////////////////////////
// Verification of the constraints //
/////////////////////////////////////

// Glpk doesn't seem to support when lower bounds == upper bounds
// So, we locate these case and add a little offset
Index = find((mps_file_mp('bounds_lower')==mps_file_mp('bounds_upper'))==%T);
mps_file_mp('bounds_lower')(Index) = mps_file_mp('bounds_lower')(Index) - 1e-4;

/////////////////////////////
// Set the constraint type //
/////////////////////////////

// 'L' - smaller than - <=
// 'E' - equality     - =
// 'G' - greater than - >=
// 'R' - Range        - <= + >=
// 'N' - Free         - no constraints

btype = mps_file_mp('constr_sense');

printf('nb of constr  = %d\n', size(mps_file_mp('constr_mat'),1));
printf('nb of var     = %d\n', size(mps_file_mp('constr_mat'),2));
printf('nb of int var = %d\n', mps_file('nb_int_var'));

[xmin,fmin,status,extra] = glpk(mps_file_mp('obj_coeff'),mps_file_mp('constr_mat'),mps_file_mp('lhs'),mps_file_mp('rhs'), ...
                                mps_file_mp('bounds_lower'),mps_file_mp('bounds_upper'),btype,var_type,param);

printf('status of problem:\n');

printf('solution status = %d\n', status);
printf('GLP_UNDEF       1  solution is undefined\n');
printf('GLP_FEAS        2  solution is feasible\n');
printf('GLP_INFEAS      3  solution is infeasible\n');
printf('GLP_NOFEAS      4  no feasible solution exists\n');
printf('GLP_OPT         5  solution is optimal\n');
printf('GLP_UNBND       6  solution is unbounded\n');

printf('error status of the solver: %d\n',extra('errnum'));
printf('1 - invalid basis \n');
printf('2 - singular matrix\n');
printf('3 - ill-conditioned matrix\n');
printf('4 - invalid bounds\n');
printf('5 - solver failed\n');
printf('6 - objective lower limit reached\n');
printf('7 - objective upper limit reached\n');
printf('8 - iteration limit exceeded\n');
printf('9 - time limit exceeded\n');
printf('10 - no primal feasible solution\n');
printf('11 - no dual feasible solution\n');
printf('12 - LP optimum not provided\n');
printf('13 - search terminated by application\n');
printf('14 - relative mip gap tolerance reached\n');
printf('15 - no primal/dual feasible solution\n');
printf('16 - no convergence\n');
printf('17 - numerical instability\n');

printf('time    = %f\n', extra('time'));
printf('mem     = %f\n', extra('mem'));

t_end = getdate();

printf('elapsed time for solving problem = %f secondes\n', etime(t_end, t_start));

if UseLinpro then
  index_F = find(mps_file_mp('dir_of_constr')==1);
  index_L = find(mps_file_mp('dir_of_constr')==2);
  index_U = find(mps_file_mp('dir_of_constr')==3);
  index_E = find(mps_file_mp('dir_of_constr')==5);
  constr_mat = mps_file_mp('constr_mat')(index_E,:);
  constr_mat = [constr_mat; -mps_file_mp('constr_mat')(index_L,:)];
  constr_mat = [constr_mat; mps_file_mp('constr_mat')(index_U,:)];
  constr_mat = [constr_mat; mps_file_mp('constr_mat')(index_F,:)];
  bound_mat  = mps_file_mp('rhs')(index_E)';
  bound_mat  = [bound_mat; -mps_file_mp('rhs')(index_L)'];
  bound_mat  = [bound_mat; mps_file_mp('rhs')(index_U)'];
  bound_mat  = [bound_mat; mps_file_mp('rhs')(index_F)'];

  t_start = t_end;

  [xmin_linpro, lambda_linpro, fmin_linpro] = linpro(mps_file_mp('obj_coeff')',full(constr_mat), bound_mat, ...
       		     		               mps_file_mp('bounds_lower')', mps_file_mp('bounds_upper')',length(index_E));

  t_end = getdate();

  printf('elapsed time for solving problem = %f secondes\n', etime(t_end, t_start));
  printf('objective function value = %f\n', fmin_linpro);
end

if PrintSol then
  printf('The solution found:\n');
  for i=1:length(xmin)
    printf('variable %d: %s - %f\n', i, part(var_type,i), xmin(i));
  end
end

cd(path_old)
