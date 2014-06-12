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
// - mps_optimslp_sample

path_new = get_absolute_file_path('mps_lpsolve_test.sce');
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
UseLinpro = %F;
PrintSol  = %F;
mps_type  = 0;

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

param = init_param();

param = add_param(param,'writemps','test.mps');

// 'anti_degen' - int - Specifies if special handling must be done to reduce degeneracy/cycling while solving.    
// Strategy codes to avoid or recover from degenerate pivots, infeasibility or numeric errors via randomized bound relaxation 
// ANTIDEGEN_NONE           0
// ANTIDEGEN_FIXEDVARS      1
// ANTIDEGEN_COLUMNCHECK    2
// ANTIDEGEN_STALLING       4
// ANTIDEGEN_NUMFAILURE     8
// ANTIDEGEN_LOSTFEAS      16
// ANTIDEGEN_INFEASIBLE    32
// ANTIDEGEN_DYNAMIC       64
// ANTIDEGEN_DURINGBB     128
// ANTIDEGEN_RHSPERTURB   256
// ANTIDEGEN_BOUNDFLIP    512
// ANTIDEGEN_DEFAULT      (ANTIDEGEN_FIXEDVARS | ANTIDEGEN_STALLING | ANTIDEGEN_INFEASIBLE)

param = add_param(param,'anti_degen', bitor(1,bitor(2,bitor(3,bitor(4,bitor(8,32))))));
//param = add_param(param,'anti_degen', 0);

// 'verbose' - int - message level
//
// MSG_NONE                 0
// MSG_PRESOLVE             1
// MSG_ITERATION            2
// MSG_INVERT               4
// MSG_LPFEASIBLE           8
// MSG_LPOPTIMAL           16
// MSG_LPEQUAL             32
// MSG_LPBETTER            64
// MSG_MILPFEASIBLE       128
// MSG_MILPEQUAL          256
// MSG_MILPBETTER         512
// MSG_MILPSTRATEGY      1024
// MSG_MILPOPTIMAL       2048
// MSG_PERFORMANCE       4096
// MSG_INITPSEUDOCOST    8192
param = add_param(param,'verbose', 64);

// 'pivoting' -> int -> Sets the pivot rule and mode. PRICER_* and PRICE_* can be ORed
// 
// Pricing methods
//   
// PRICER_FIRSTINDEX        0
// PRICER_DANTZIG           1
// PRICER_DEVEX             2
// PRICER_STEEPESTEDGE      3
// PRICER_LASTOPTION        PRICER_STEEPESTEDGE
//
// Pricing strategies 
//   
// PRICE_PRIMALFALLBACK     4    // In case of Steepest Edge, fall back to DEVEX in primal
// PRICE_MULTIPLE           8    // Enable multiple pricing (primal simplex)
// PRICE_PARTIAL           16    // Enable partial pricing
// PRICE_ADAPTIVE          32    // Temporarily use alternative strategy if cycling is detected 
// PRICE_RANDOMIZE        128    // Adds a small randomization effect to the selected pricer 
// PRICE_AUTOPARTIAL      256    // Detect and use data on the block structure of the model (primal) 
// PRICE_AUTOMULTIPLE     512    // Automatically select multiple pricing (primal simplex) 
// PRICE_LOOPLEFT        1024    // Scan entering/leaving columns left rather than right 
// PRICE_LOOPALTERNATE   2048    // Scan entering/leaving columns alternatingly left/right 
// PRICE_HARRISTWOPASS   4096    // Use Harris' primal pivot logic rather than the default 
// PRICE_FORCEFULL       8192    // Non-user option to force full pricing 
// PRICE_TRUENORMINIT   16384    // Use true norms for Devex and Steepest Edge initializations 
// default value: PRICER_DEVEX | PRICE_ADAPTIVE (34)
param = add_param(param,'pivoting', 34);

// 'epsb' -> double -> the value that is used as a tolerance for the Right Hand Side (RHS) to determine 
//                     whether a value should be considered as 0
param = add_param(param,'epsb', 1e-10);

// 'epsd' -> double -> the value that is used as a tolerance for the reduced costs to determine 
//                     whether a value should be considered as 0.
param = add_param(param,'epsd', 1e-9);

// 'epspivot' -> double -> the value that is used as a tolerance for the pivot element to determine
//                         whether a value should be considered as 0.
param = add_param(param,'epspivot', 2e-7);

// 'epsel' -> double -> the value that is used as a tolerance for rounding values to zero.
param = add_param(param,'epsel', 1e-12);

// 'epsint' -> double -> the tolerance that is used to determine whether a floating-point number is in fact an integer.
param = add_param(param,'epsint', 1e-7);

// 'epsperturb' -> double -> the value that is used as perturbation scalar for degenerative problems.
param = add_param(param,'epsperturb', 1e-5);

// 'infinite' -> double -> Specifies the practical value for "infinite".
param = add_param(param,'infinite', 1e30);

// 'break_at_first' -> boolean -> Specifies if the branch-and-bound algorithm stops at first found solution.
param = add_param(param,'break_at_first', 0);

// 'break_at_value' -> double -> Specifies if the branch-and-bound algorithm stops when the object value is better than a given value.
// default value: (-) infinity
param = add_param(param,'break_at_value', -1e30);

// 'basiscrash' -> int -> Determines a starting base.
//
// Basis crash options
//
// CRASH_NONE               0
// CRASH_NONBASICBOUNDS     1
// CRASH_MOSTFEASIBLE       2
// CRASH_LEASTDEGENERATE    3
param = add_param(param,'basiscrash', 0);

// 'bb_depthlimit' -> int -> Sets the maximum branch-and-bound depth.
param = add_param(param,'bb_depthlimit',-50);

// 'bb_floorfirst' -> int -> Specifies which branch to take first in branch-and-bound algorithm.
//
// BRANCH_CEILING           0
// BRANCH_FLOOR             1
// BRANCH_AUTOMATIC         2
// BRANCH_DEFAULT           3
param = add_param(param,'bb_floorfirst', 2);

//
// B&B strategies 
//
// NODE_FIRSTSELECT         0
// NODE_GAPSELECT           1
// NODE_RANGESELECT         2
// NODE_FRACTIONSELECT      3
// NODE_PSEUDOCOSTSELECT    4
// NODE_PSEUDONONINTSELECT  5    // Kjell Eikland #1 - Minimize B&B depth
// NODE_PSEUDOFEASSELECT   (NODE_PSEUDONONINTSELECT+NODE_WEIGHTREVERSEMODE)
// NODE_PSEUDORATIOSELECT   6    // Kjell Eikland #2 - Minimize a "cost/benefit" ratio 
// NODE_USERSELECT          7
// NODE_STRATEGYMASK        (NODE_WEIGHTREVERSEMODE-1) // Mask for B&B strategies
// NODE_WEIGHTREVERSEMODE   8
// NODE_BRANCHREVERSEMODE  16
// NODE_GREEDYMODE         32
// NODE_PSEUDOCOSTMODE     64
// NODE_DEPTHFIRSTMODE    128
// NODE_RANDOMIZEMODE     256
// NODE_GUBMODE           512
// NODE_DYNAMICMODE      1024
// NODE_RESTARTMODE      2048
// NODE_BREADTHFIRSTMODE 4096
// NODE_AUTOORDER        8192
// NODE_RCOSTFIXING     16384
// NODE_STRONGINIT      32768
param = add_param(param,'bb_rule', 17445); // NODE_PSEUDONONINTSELECT + NODE_GREEDYMODE + NODE_DYNAMICMODE + NODE_RCOSTFIXING 

// 'debug' -> boolean -> Sets a flag if all intermediate results and the branch-and-bound decisions must be printed while solving.
param = add_param(param,'debug', 0);

// 'lag_trace' -> boolean -> Sets a flag if Lagrangian progression must be printed while solving.
param = add_param(param,'lag_trace', 0);

// 'maxpivot' -> int -> Sets the maximum number of pivots between a re-inversion of the matrix.
param = add_param(param,'maxpivot', 250); // The default is 250 for the LUSOL bfp and 42 for the other BFPs. 

// 'mip_gap_abs' -> boolean -> If TRUE then the absolute MIP gap is set, else the relative MIP gap
// 'mip_gap_gap' -> double  -> The MIP gap.
param = add_param(param,'mip_gap_abs', 0);
param = add_param(param,'mip_gap_gap', 1e-11);


// 'preferdual' -> int -> Sets the desired combination of primal and dual simplex algorithms.
//
// SIMPLEX_UNDEFINED        0
// SIMPLEX_Phase1_PRIMAL    1
// SIMPLEX_Phase1_DUAL      2
// SIMPLEX_Phase2_PRIMAL    4
// SIMPLEX_Phase2_DUAL      8
// SIMPLEX_DYNAMIC         16
// SIMPLEX_AUTODUALIZE     32
//
// SIMPLEX_PRIMAL_PRIMAL   (SIMPLEX_Phase1_PRIMAL + SIMPLEX_Phase2_PRIMAL)
// SIMPLEX_DUAL_PRIMAL     (SIMPLEX_Phase1_DUAL   + SIMPLEX_Phase2_PRIMAL)
// SIMPLEX_PRIMAL_DUAL     (SIMPLEX_Phase1_PRIMAL + SIMPLEX_Phase2_DUAL)
// SIMPLEX_DUAL_DUAL       (SIMPLEX_Phase1_DUAL   + SIMPLEX_Phase2_DUAL)
// SIMPLEX_DEFAULT         (SIMPLEX_DUAL_PRIMAL)
param = add_param(param,'preferdual', 6); // The default is SIMPLEX_DUAL_PRIMAL (6). 

// 'simplextype' -> int -> Sets the desired combination of primal and dual simplex algorithms.
// SIMPLEX_PRIMAL_PRIMAL (5)  Phase1 Primal, Phase2 Primal
// SIMPLEX_DUAL_PRIMAL   (6)  Phase1 Dual, Phase2 Primal (default value)
// SIMPLEX_PRIMAL_DUAL   (9)  Phase1 Primal, Phase2 Dual
// SIMPLEX_DUAL_DUAL     (10) Phase1 Dual, Phase2 Dual  
param = add_param(param,'simplextype', 6);

// 'presolve' -> int -> Do presolve in 1
//
// PRESOLVE_NONE            0
// PRESOLVE_ROWS            1
// PRESOLVE_COLS            2
// PRESOLVE_LINDEP          4
// PRESOLVE_SOS            32
// PRESOLVE_REDUCEMIP      64
// PRESOLVE_KNAPSACK      128  // Implementation not tested completely
// PRESOLVE_ELIMEQ2       256
// PRESOLVE_IMPLIEDFREE   512
// PRESOLVE_REDUCEGCD    1024
// PRESOLVE_PROBEFIX     2048
// PRESOLVE_PROBEREDUCE  4096
// PRESOLVE_ROWDOMINATE  8192
// PRESOLVE_COLDOMINATE 16384  // Reduced functionality, should be expanded 
// PRESOLVE_MERGEROWS   32768
// PRESOLVE_IMPLIEDSLK  65536
// PRESOLVE_COLFIXDUAL 131072
// PRESOLVE_BOUNDS     262144
// PRESOLVE_LASTMASKMODE    (PRESOLVE_DUALS - 1)
// PRESOLVE_DUALS      524288
// PRESOLVE_SENSDUALS 1048576
//
// Can be ORed
param = add_param(param,'presolve_do', 0);
//param = add_param(param,'presolve_do', 1048576);
//param = add_param(param,'presolve_do', bitor(1,bitor(2,bitor(3,32768))));

// 'presolve_maxloops' -> int -> maxloops - Specifies if a presolve must be done before solving.
param = add_param(param,'presolve_maxloops', 100);

// 'scalelimit'      -> double -> Sets the relative scaling convergence criterion for the active scaling mode; 
//                                the integer part specifies the maximum number of iterations.
param = add_param(param,'scalelimit',5);

// 'scaling' -> int -> Specifies which scaling algorithm must be used.
//
// SCALE_NONE               0
// SCALE_EXTREME            1
// SCALE_RANGE              2
// SCALE_MEAN               3
// SCALE_GEOMETRIC          4
// SCALE_FUTURE1            5
// SCALE_FUTURE2            6
// SCALE_CURTISREID         7   // Override to Curtis-Reid "optimal" scaling
//  
// Alternative scaling weights 
//  
// SCALE_LINEAR             0
// SCALE_QUADRATIC          8
// SCALE_LOGARITHMIC       16
// SCALE_USERWEIGHT        31
// SCALE_MAXTYPE            (SCALE_QUADRATIC-1)
//  
// Scaling modes 
//  
// SCALE_POWER2            32   // As is or rounded to power of 2 
// SCALE_EQUILIBRATE       64   // Make sure that no scaled number is above 1 
// SCALE_INTEGERS         128   // Apply to integer columns/variables 
// SCALE_DYNUPDATE        256   // Apply incrementally every solve() 
// SCALE_ROWSONLY         512   // Override any scaling to only scale the rows 
// SCALE_COLSONLY        1024   // Override any scaling to only scale the rows 
//
// Can be ORed
param = add_param(param,'scaling', 196); // SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS (196). 

// 'solutionlimit' -> int -> Sets the solution number that must be returned (for problem with binary, integer or semicontinuous variables.
param = add_param(param,'solutionlimit', 1);

// 'timeout' -> int -> set a timeout in second (0: no timeout). 
param = add_param(param,'timeout', 10000);

// 'trace' -> boolean -> Sets a flag if pivot selection must be printed while solving.
param = add_param(param,'trace', 0);

// 'negrange' -> double -> Set negative value below which variables are split into a negative and a positive part.
param = add_param(param,'negrange', -1e6);

// 'epslevel' -> int -> This is a simplified way of specifying multiple eps thresholds that are "logically" consistent.
// EPS_TIGHT  (0) Very tight epsilon values (default)
// EPS_MEDIUM (1) Medium epsilon values
// EPS_LOOSE  (2) Loose epsilon values
// EPS_BAGGY  (3) Very loose epsilon values
param = add_param(param,'epslevel', 0);

// 'improve' -> int -> Specifies the iterative improvement level.
// IMPROVE_NONE      (0) improve none
// IMPROVE_SOLUTION  (1) Running accuracy measurement of solved equations based on Bx=r (primal simplex), remedy is refactorization.
// IMPROVE_DUALFEAS  (2) Improve initial dual feasibility by bound flips (highly recommended, and default)
// IMPROVE_THETAGAP  (4) Low-cost accuracy monitoring in the dual, remedy is refactorization
// IMPROVE_BBSIMPLEX (8) By default there is a check for primal/dual feasibility at optimum only for the relaxed problem, this also activates the test at the node level
param = add_param(param,'improve', 6); // IMPROVE_DUALFEAS + IMPROVE_THETAGAP (6). 

// 'bounds_tighter' -> boolean -> Specifies if set bounds may only be tighter or also less restrictive.
param = add_param(param,'bounds_tighter', 0);

// 'sense' -> int -> optimization direction (-1 minimization, 1 maximization)
param = add_param(param,'sense', -1);


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

[xmin,fmin,status,extra] = lpsolve(mps_file_mp('obj_coeff'),mps_file_mp('constr_mat'),mps_file_mp('lhs'),mps_file_mp('rhs'), ...
                                   mps_file_mp('bounds_lower'),mps_file_mp('bounds_upper'),btype,var_type,param);

// A negative status means that an error occured and no results will be returned
printf('status of problem:\n');
printf('UNKNOWNERROR            -5\n');
printf('DATAIGNORED             -4\n');
printf('NOBFP                   -3\n');
printf('NOMEMORY                -2\n');
printf('NOTRUN                  -1\n');
printf('OPTIMAL                  0\n');
printf('SUBOPTIMAL               1\n');
printf('INFEASIBLE               2\n');
printf('UNBOUNDED                3\n');
printf('DEGENERATE               4\n');
printf('NUMFAILURE               5\n');
printf('USERABORT                6\n');
printf('TIMEOUT                  7\n');
printf('RUNNING                  8\n');
printf('PRESOLVED                9\n');
printf('status  = %d\n', status);
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

  [xmin_linpro, lambda_linpro, fmin_linpro] = linpro(mps_file_mp('obj_coeff'),full(constr_mat), bound_mat, ...
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
