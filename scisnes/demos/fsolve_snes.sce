// the pipe test problem on the scilab fsolve_snes function

lines(0);
funcprot_old = funcprot();
funcprot(0);

// Load the necessary functions (plot_pipe_graph, etc ...)
getd('.');

Plot = %F;

global jac_ev;
global fun_ev;

use_pipe_pb = %F; // requires metanet

jac_ev = 0;
fun_ev = 0;

params = init_param();

params = add_param(params, 'Qmax', 250); // m^3/s
params = add_param(params, 'Pmax', 30);  // Bars
params = add_param(params, 'Pmin', 0);   // Bars
params = add_param(params, 'k_pipe', 5/100); 
params = add_param(params, 'Debug', %F); 

// Objective function for fsolve
function res = fsolve_obj_pipe(x)
  res = compute_pipe_eq(x,g_in,I_sources,params);
endfunction

// Analytical gradient for fsolve
function res = fsolve_dobj_pipe(x)
  global jac_ev;
  jac_ev = jac_ev + 1;
  res = full(compute_pipe_deq(x,g_in,I_sources,params)');
endfunction

// Finite differences gradient for fsolve
function res = fsolve_dobj_diff_pipe(x)
  res = derivative(fsolve_obj_pipe,x);
endfunction

function [f,J] = fsolver_dobj_pipe(x)
  global fun_ev;
  fun_ev = fun_ev + 1;
  f = compute_pipe_eq(x,g_in,I_sources,params);
  [nargout,nargin] = argn();
  if nargout > 1
    global jac_ev;
    jac_ev = jac_ev + 1;
    //J = fsolve_dobj_diff_pipe(x);
    J = compute_pipe_deq(x,g_in,I_sources,params)';
    //J = full(J);
  end
endfunction

if use_pipe_pb then
  //////////////////////////////////////
  // Initialization of the pipe graph //
  //////////////////////////////////////
  
  // Available tests:
  // - 'simple'
  // - 'simple_loop_1'
  // - 'simple_loop_2'
  // - 'simple_loop_3'
  // - 'test_hc_1'
  // - 'test_hc_1_bis'
  // - 'test_hc_2'
  // - 'test_hc_3'
  // - 'test_pr1'
  // - 'test_pr2'
  // - 'test_pr3'
  
  [g_in, I_sources, I_leaves] = build_test('simple_loop_1',%F);
  
  ///////////////////////////////////////////////////
  // Initialization of the flow and pressures (x0) //
  ///////////////////////////////////////////////////
  
  [x0, lower, upper] = pipe_init(g_in, I_sources,params);
  x0 = (upper - lower).*rand(lower) + lower;
end

/////////////////////
// Test from Petsc //
/////////////////////

function [F,J] = petsc_pb_1(x)
  global jac_ev;
  jac_ev = jac_ev + 1;
  global fun_ev;
  fun_ev = fun_ev + 1;
  F(1,1) = sin(3.0*x(1)) + x(1);
  F(2,1) = x(2);
  J(1,1) = 3*cos(3*x(1)) + 1;
  J(1,2) = 0;
  J(2,1) = 0;
  J(2,2) = 1;
endfunction

////////////////////////
// Rosenbrock problem //
////////////////////////

function [F,J] = genrosenbrock(x)
  disp(size(x))
  global jac_ev;
  jac_ev = jac_ev + 1;
  global fun_ev;
  fun_ev = fun_ev + 1;
  // Derived from the generalized n-dimensional 
  // Rosenbrock function, returns a dense Jacobian 
  n = length(x);         
  if n == 0, error('Input vector, x, is empty.'); end
  if modulo(n,2) ~= 0, 
    error('Input vector, x ,must have an even number of components.');
  end
  // Evaluate the vector function
  odds  = 1:2:n;
  evens = 2:2:n;
  F = zeros(n,1);
  F(odds,1)  = 1-x(odds);
  F(evens,1) = 10.*(x(evens)-x(odds).^2); 
  // Evaluate the Jacobian matrix if nargout > 1
  [nargout,nargin] = argn();
  if nargout > 1
    c = -ones(n/2,1);    C = sparse([odds; odds]',c,[n,n]);
    d = 10*ones(n/2,1);  D = sparse([evens; evens]',d,[n,n]);
    e = -20.*x(odds);    E = sparse([evens; odds]',e,[n,n]);
    J = C + D + E;
  end
  disp(size(F))
  disp(size(J))  
endfunction


function [f,J] = genrosenbrockdense(x)
  [nargout,nargin] = argn();
   
  if nargout == 1 
    f = genrosenbrock(x);
  end
  
  if nargout > 1 
    [f,J] = genrosenbrock(x);
    J = full(J);
  end
endfunction

////////////////////////////
// Set the solver options //
////////////////////////////

// PETSC_DEFAULT == -2
// If you want to get the default parameter for a function, send PETSC_DEFAULT

params_snes = init_param();

// snes_type: Sets the method for the nonlinear solver. 
// SNESLS - 'ls' - Newton's method with line search (systems of nonlinear equations)
// SNESTR - 'tr' - Newton's method with trust region (systems of nonlinear equations) 
params_snes = add_param(params_snes,'snes_type', 'ls');

// snes_tr_tol: Sets the trust region parameter tolerance. 
params_snes = add_param(params_snes,'snes_tr_tol', 0);

// snes_max_nl_step_fail: Sets the maximum number of unsuccessful steps attempted by the nonlinear solver before it gives up. 
params_snes = add_param(params_snes,'snes_max_nl_step_fail', 2);

// snes_max_lin_solve_fail: the number of failed linear solve attempts allowed before SNES returns with a diverged reason of SNES_DIVERGED_LINEAR_SOLVE 
params_snes = add_param(params_snes,'snes_max_lin_solve_fail', 2);

// snes_lag_precond: Determines when the preconditioner is rebuilt in the nonlinear solve. 
// lag: 	
//   - -1 indicates NEVER rebuild, 
//   -  1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 
//   -  2 means every second time the Jacobian is built etc. 
//   - -2 indicates rebuild preconditioner at next chance but then never rebuild after that 
params_snes = add_param(params_snes,'snes_lag_precond',1);

// snes_lag_jac: Determines when the Jacobian is rebuilt in the nonlinear solve. 
// lag: 	
//   - -1 indicates NEVER rebuild, 
//   -  1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 
//   -  2 means every second time the Jacobian is built etc. 
//   - -2 means rebuild at next chance but then never again 
params_snes = add_param(params_snes,'snes_lag_jac', 1);

// snes_pc_type: set the preconditionner type. Can be PCNONE, PCJACOBI, PCILU, or PCBJACOB
// PCNONE            'none'           - This is used when you wish to employ a nonpreconditioned Krylov method. 
// PCJACOBI          'jacobi'         - Jacobi (i.e. diagonal scaling preconditioning) 
// PCSOR             'sor'            - (S)SOR (successive over relaxation, Gauss-Seidel) preconditioning 
// PCLU              'lu'             - Uses a direct solver, based on LU factorization, as a preconditioner 
// PCBJACOBI         'bjacobi'        - Use block Jacobi preconditioning, each block is (approximately) solved with its own KSP object. 
// PCEISENSTAT       'eisenstat'      - An implementation of SSOR (symmetric successive over relaxation, symmetric Gauss-Seidel) preconditioning 
//                                      that incorporates Eisenstat's trick to reduce the amount of computation needed.
// PCILU             'ilu'            - Incomplete factorization preconditioners. 
// PCICC             'icc'            - Incomplete Cholesky factorization preconditioners. 
// PCASM             'asm'            - Use the (restricted) additive Schwarz method, each block is (approximately) solved with its own KSP object. 
// PCSPAI            'spai'           - Use the Sparse Approximate Inverse method of Grote and Barnard as a preconditioner 
//                                      (SIAM J. Sci. Comput.; vol 18, nr 3) 
// PCNN              'nn'             - Balancing Neumann-Neumann for scalar elliptic PDEs. 
// PCCHOLESKY        'cholesky'       - Uses a direct solver, based on Cholesky factorization, as a preconditioner 
// PCPROMETHEUS      'prometheus'     - Prometheus (i.e. diagonal scaling preconditioning) 
// PCOPENMP          'openmp'         - Runs a preconditioner for a single process matrix across several MPI processes 
// PCSUPPORTGRAPH    'supportgraph'   - SupportGraph (i.e. diagonal scaling preconditioning) 
params_snes = add_param(params_snes, 'snes_pc_type', 'lu');

// snes_ksp_rtol, snes_ksp_abstol, snes_ksp_dtol, snes_ksp_maxits:
// rtol   - the relative convergence tolerance (relative decrease in the residual norm)
// abstol - the absolute convergence tolerance (absolute size of the residual norm)
// dtol   - the divergence tolerance (amount residual can increase before KSPDefaultConverged() concludes that the method is diverging)
// maxits - maximum number of iterations to use 
params_snes = add_param(params_snes,'snes_ksp_rtol',   1e-40);
params_snes = add_param(params_snes,'snes_ksp_abstol', 1e-40);
params_snes = add_param(params_snes,'snes_ksp_dtol',   1e6);
params_snes = add_param(params_snes,'snes_ksp_maxits', 100000);

// snes_ksp_comp_sing_val: Sets a flag so that the extreme singular values will be calculated via a Lanczos or Arnoldi process as 
// the linear system is solved.
params_snes = add_param(params_snes,'snes_ksp_comp_sing_val', 1);

// snes_rtol, snes_abstol, snes_stol, snes_maxit, snes_maxf:
// abstol - absolute convergence tolerance
// rtol   - relative convergence tolerance
// stol   - convergence tolerance in terms of the norm of the change in the solution between steps
// maxit  - maximum number of iterations
// maxf   - maximum number of function evaluations
params_snes = add_param(params_snes,'snes_rtol',   1e-40);
params_snes = add_param(params_snes,'snes_abstol', 1e-40);
params_snes = add_param(params_snes,'snes_stol',   1e-40);
params_snes = add_param(params_snes,'snes_maxit',  100);
params_snes = add_param(params_snes,'snes_maxf',   1000);

// Builds KSP for a particular solver. 
// KSPRICHARDSON 'richardson'
// KSPCHEBYCHEV  'chebychev'
// KSPCG         'cg'
//   KSPCGNE       'cgne'
//   KSPNASH       'nash' petsc-3
//   KSPSTCG       'stcg'
//   KSPGLTR       'gltr'
// KSPGMRES      'gmres'
//   KSPFGMRES     'fgmres'
//   KSPLGMRES     'lgmres' // Be careful: seems to have a bug
// KSPTCQMR      'tcqmr'
// KSPBCGS       'bcgs'
// KSPIBCGS        'ibcgs' petsc-3
// KSPBCGSL        'bcgsl'
// KSPCGS        'cgs'
// KSPTFQMR      'tfqmr'
// KSPCR         'cr'
// KSPLSQR       'lsqr'
// KSPPREONLY    'preonly'
// KSPQCG        'qcg' petsc-3
// KSPBICG       'bicg'
// KSPMINRES     'minres'
// KSPSYMMLQ     'symmlq'
// KSPLCD        'lcd'
params_snes = add_param(params_snes,'snes_ksp_type','gmres');

// snes_ls_type Sets the line search routine to be used by the method SNESLS. 
// 1 - SNESLineSearchCubic     - default line search
// 2 - SNESLineSearchQuadratic - quadratic line search
// 3 - SNESLineSearchNo        - the full Newton step (actually not a line search)
// 4 - SNESLineSearchNoNorms   - the full Newton step (calculating no norms; faster in parallel) 
params_snes = add_param(params_snes,'snes_ls_type', 1);

// snes_ls_set_alpha, snes_ls_set_maxstep:
// Sets the parameters associated with the line search routine in the Newton-based method SNESLS.
// Input Parameters:
// alpha   - The scalar such that .5*f_{n+1} . f_{n+1} <= .5*f_n . f_n - alpha |f_n . J . f_n|
// maxstep - The maximum norm of the update vector
// steptol - The minimum norm fraction of the original step after scaling
// Note: Pass in PETSC_DEFAULT for any parameter you do not wish to change.
// We are finding the zero of f() so the one dimensional minimization problem we are solving in the line search is minimize
// .5*f(x_n + lambda*step_direction) . f(x_n + lambda*step_direction)
params_snes = add_param(params_snes,'snes_ls_set_alpha', -2);
params_snes = add_param(params_snes,'snes_ls_set_maxstep', -2);
params_snes = add_param(params_snes,'snes_ls_set_steptol', -2); // Only in petsc-2. This option has been removed in petsc-3

////////////////
// Resolution //
////////////////

NbVar = 4;
x0 = 4*ones(NbVar,1).*rand(NbVar,1) - 2*ones(NbVar,1);

// With analytical Jacobian and with fsolve
tic();
//[x_opt, status] = fsolver_snes(petsc_pb_1,x0,params_snes);
[x_opt, status] = fsolver_snes(genrosenbrockdense,x0,params_snes);
//[x_opt, status] = fsolver_snes(genrosenbrock,x0,params_snes);
//[x_opt, status] = fsolver_snes(fsolver_dobj_pipe,x0,params_snes);
t = toc();

//g_opt = pipe_update_graph(g_in,x_opt);
//if Plot then 
//  plot_pipe_graph(g_opt,'test for fsolve - with analytical Jacobian','X','Y'); 
//end

printf('with Jacobian - result: \n');
printf('info_snes                      = %d\n', status('info_snes'));
printf('info_ksp                       = %d\n', status('info_ksp'));
printf('non linear iterations          = %d\n', status('its'));
printf('linear iterations              = %d\n', status('lits'));
printf('failed non linear iterations   = %d\n', status('fail_its'));
printf('failed linear iterations       = %d\n', status('fail_lits'));
printf('number of function evaluations = %d\n', status('fev'));
printf('elapsed time = %f\n' ,t);

///////////////////////
// SNES return codes //
///////////////////////
// converged
// SNES_CONVERGED_FNORM_ABS         =  2, ||F|| < atol
// SNES_CONVERGED_FNORM_RELATIVE    =  3, ||F|| < rtol*||F_initial||
// SNES_CONVERGED_PNORM_RELATIVE    =  4, Newton computed step size small; || delta x || < tol
// SNES_CONVERGED_ITS               =  5, maximum iterations reached
// SNES_CONVERGED_TR_DELTA          =  7,
// diverged
// SNES_DIVERGED_FUNCTION_DOMAIN    = -1,  
// SNES_DIVERGED_FUNCTION_COUNT     = -2,  
// SNES_DIVERGED_LINEAR_SOLVE       = -3, 
// SNES_DIVERGED_FNORM_NAN          = -4, 
// SNES_DIVERGED_MAX_IT             = -5, means that the solver reached the maximum number of iterations without satisfying any convergence criteria. 
//                                        SNES_CONVERGED_ITS means that SNESSkipConverged() was chosen as the convergence test; 
//                                        thus the usual convergence criteria have not been checked and may or may not be satisfied. 
// SNES_DIVERGED_LS_FAILURE         = -6,
// SNES_DIVERGED_LOCAL_MIN          = -8  || J^T b || is small, implies converged to local minimum of F()
//                                        This can only occur when using the line-search variant of SNES. 
//                                        The line search wants to minimize Q(alpha) = 1/2 || F(x + alpha s) ||^2_2 this occurs at 
//                                        Q'(alpha) = s^T F'(x+alpha s)^T F(x+alpha s) = 0. If s is the Newton direction - F'(x)^(-1)F(x) 
//                                        then you get Q'(alpha) = -F(x)^T F'(x)^(-1)^T F'(x+alpha s)F(x+alpha s); 
//                                        when alpha = 0 Q'(0) = - ||F(x)||^2_2 which is always NEGATIVE if F'(x) is invertible. 
//                                        This means the Newton direction is a descent direction and the line search should succeed 
//                                        if alpha is small enough.
//                                        If F'(x) is NOT invertible AND F'(x)^T F(x) = 0 then Q'(0) = 0 and the Newton direction is NOT a
//                                        descent direction so the line search will fail. 
//                                        All one can do at this point is change the initial guess and try again.
//                                        An alternative explanation: Newton's method can be regarded as replacing the function with its
//                                        linear approximation and minimizing the 2-norm of that. 
//                                        That is F(x+s) approx F(x) + F'(x)s so we minimize || F(x) + F'(x) s ||^2_2; do this using Least Squares. 
//                                        If F'(x) is invertible then s = - F'(x)^(-1)F(x) otherwise F'(x)^T F'(x) s = -F'(x)^T F(x). 
//                                        If F'(x)^T F(x) is NOT zero then there exists a nontrival (that is F'(x)s != 0) solution to the
//                                        equation and this direction is s = - [F'(x)^T F'(x)]^(-1) F'(x)^T F(x) 
//                                        so Q'(0) = - F(x)^T F'(x) [F'(x)^T F'(x)]^(-T) F'(x)^T 
//                                        F(x) = - (F'(x)^T F(x)) [F'(x)^T F'(x)]^(-T) (F'(x)^T F(x)). 
//                                        Since we are assuming (F'(x)^T F(x)) != 0 and F'(x)^T F'(x) has no negative eigenvalues Q'(0) < 0 
//                                        so s is a descent direction and the line search should succeed for small enough alpha.
//                                        Note that this RARELY happens in practice. Far more likely the linear system is not being 
//                                        solved (well enough?) or the Jacobian is wrong. 
// SNES_CONVERGED_ITERATING         =  0

//////////////////////
// KSP return codes //
//////////////////////
// Converged
// KSP_CONVERGED_RTOL               =  2,
// KSP_CONVERGED_ATOL               =  3,
// KSP_CONVERGED_ITS                =  4,
// KSP_CONVERGED_CG_NEG_CURVE       =  5,
// KSP_CONVERGED_CG_CONSTRAINED     =  6,
// KSP_CONVERGED_STEP_LENGTH        =  7,
// KSP_CONVERGED_HAPPY_BREAKDOWN    =  8,
// Diverged
// KSP_DIVERGED_NULL                = -2,
// KSP_DIVERGED_ITS                 = -3,
// KSP_DIVERGED_DTOL                = -4,
// KSP_DIVERGED_BREAKDOWN           = -5,
// KSP_DIVERGED_BREAKDOWN_BICG      = -6,
// KSP_DIVERGED_NONSYMMETRIC        = -7,
// KSP_DIVERGED_INDEFINITE_PC       = -8,
// KSP_DIVERGED_NAN                 = -9,
// KSP_DIVERGED_INDEFINITE_MAT      = -10,
// KSP_CONVERGED_ITERATING          =  0  ierr = VecDestroy(x);CHKERRQ(ierr); 

funcprot(funcprot_old);
