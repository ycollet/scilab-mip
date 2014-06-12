// A LP example which shows all the potentials of LPSOLVE SCILAB

TestLP  = %T;
TestSOS = %T;

param = init_param();

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

if TestLP then
  c = [-1 -1 -3 -2 -2]; // from lpsolve sos documentation
  a = [-1 -1  1  1  0; ...
        1  0  1 -3  0];
  b = [30 30]';
  lb = -1000*[1, 1, 1, 1, 1];
  ub =  [40, 1, 1000, 1000, 1];

  vartype = 'IIIIC';
  constrtype = 'NLL';

  printf('test without special constraints\n');

  constraints = [];

  [xmin,fmin,status,extra] = lpsolve(c,a,b,b,lb,ub,constrtype,vartype,param,constraints);

  printf('solution found: \n');disp(xmin');
  printf('status = %d\n',status);
end


if TestSOS then
  c = -[1 1 3 2 2]; // from lpsolve sos documentation

  a = [-1 -1  1  1  0; ...
        1  0  1 -3  0];
        
  b_lo = [0   0];
  b_up = [30 30];
  
  lb = [ 0, 0, 0, 0, 0];
  ub = [40, 40, 40, 40, 40];

  constrtype = 'LL';

  which  = [1 2 3 4 5];
  weight = [0 1 2 3 4];
  
  vartype = 'IIIII';

  printf('test of SOS order 1\n');

  constraints = [];
  constraints = init_constraint();
  constraints = add_constraint_sos(constraints,1,which,weight,2); // the SOS constraint is to be applied to constraint no 2

  [xmin,fmin,status,extra] = lpsolve(c,a,b_lo,b_up,lb,ub,constrtype,vartype,param,constraints);

  printf('solution found: \n');disp(xmin');
  printf('status = %d\n',status);
  printf('value of objective function = %f\n', c*xmin);
  printf('value of constraint no 1 = %f <= %f\n', a(1,:)*xmin, b_up(1));
  printf('value of constraint no 2 = %f <= %f\n', a(2,:)*xmin, b_up(2));
  
  printf('test of SOS order 2\n');

  constraints = [];
  constraints = init_constraint();
  constraints = add_constraint_sos(constraints,2,which,weight,2); // the SOS constraint is to be applied to constraint no 2

  [xmin,fmin,status,extra] = lpsolve(c,a,b_lo,b_up,lb,ub,constrtype,vartype,param,constraints);

  printf('solution found: \n');disp(xmin');
  printf('status = %d\n',status);
  printf('value of objective function = %f\n', c*xmin);
  printf('value of constraint no 1 = %f <= %f\n', a(1,:)*xmin, b_up(1));
  printf('value of constraint no 2 = %f <= %f\n', a(2,:)*xmin, b_up(2));

  // Warning: CBC doesn't accept SOS constraints with an order equal to 0 or an order strictly above 2.
end

