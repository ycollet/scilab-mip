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

path_new = get_absolute_file_path('mps_symphony_test.sce');
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
Solve_miplib3  = %F; // 65 files
Solve_sample   = %F; // 17 files
Solve_infeas   = %F; // 29 files
Solve_big      = %F; // 1 file
Solve_netlib   = %F; // 91 files
Solve_data     = %T; // 2 files
Solve_optimslp = %F; // 3 files

///////////////////////////
// Set global parameters //
///////////////////////////

Log       = 1;
Debug     = 0;
mps_index = 2; // 1 LP - 2 MIP
// netlib index 5: fobj = nan apres 1000 iterations avec interior
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
if Solve_optimslp then mps_filename = mps_optimslp_sample(mps_index); end

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

param = init_param();
//param = add_param(param, 'verbosity', 3);

/////////////////////////
//// Global parameters //
/////////////////////////
//// 'verbosity'  integer (0). Sets the verbosity of all modules to the given value. In general, the greater this number the more verbose each module is.
////                           Experiment to find out what this means.
//param = add_param(param, 'verbosity', 0);
//
//// 'random_seed'  integer (17). A random seed. 
//param = add_param(param, 'random_seed', 17);
//
//// 'granularity'  double (1e-6). Should be set to 'the minimum difference between two distinct objective function values' less the epsilon tolerance. 
////                               E.g., if every variable is integral and the objective coefficients are integral then for any feasible solution the objective
////                               value is integer, so granularity could be correctly set to .99999.
//param = add_param(param, 'granularity', 1e-6);    
//
//// 'upper_bound'  double (none) . The value of the best known upper bound.
//param = add_param(param, 'upper_bound', 1e100);    
//
//// 'lower_bound'  double (no lower bound). This parameter is used if the user wants to artificially impose a lower bound.
//param = add_param(param, 'lower_bound', -1e100);    
//
//// 'probname'  string (empty string). The name of the problem name.
//param = add_param(param, 'probname', 'test');    
//
//// 'infile_name'  string (empty string). The name of the input file that was read by '-F' or the '-L' flag.
//param = add_param(param, 'infile_name', 'test.mps');    
//
//// 'writemps'  string ('test.mps'). The name of the output file generated (if the option is present).
//param = add_param(param, 'writemps', 'test.mps');    
//
////////////////////////////////
//// Master module parameters //
////////////////////////////////
//// 'M_verbosity'  integer (0).
param = add_param(param, 'M_verbosity', 3);
//
//// 'M_random_seed'  integer (17). A random seed just for the Master module.
//param = add_param(param, 'M_random_seed', 17);
//
//// 'upper_bound_estimate'  double (no estimate). This parameter is used if the user wants to provide an estimate of the optimal value which will help guide the
////                                 search. This is used in conjunction with the diving strategy BEST ESTIMATE.
//param = add_param(param, 'upper_bound_estimate', 1e100);    
//
//// 'tm_exe', 'dg_exe'  strings ('tm', 'dg'). The name of the executable files of the TM and DG modules. Note that the TM executable name may have extensions 
////                             that depend on the configuration of the modules, but the default is always set to the file name produced by the makefile.
////                             If you change the name of the treemanager executable from the default, you must set this parameter to the new name.
//param = add_param(param, 'tm_exe', 'tm');    
//param = add_param(param, 'dg_exe', 'dg');    
//
//// 'tm_debug', 'dg_debug'  boolean (both FALSE). Whether these modules should be started under a debugger or not (see 5.6.2 for more details on this).
//param = add_param(param, 'tm_debug', 0);
//param = add_param(param, 'dg_debug', 0);
//
//// 'tm_machine'  string (empty string). On which processor of the virtual machine the TM should be run. Leaving this parameter as an empty string means arbitrary selection.
////param = add_param(param, 'tm_machine', '');    
//
////// 'do_draw_graph'  boolean (FALSE). Whether to start up the DG module or not (see Section 5.6.4 for an introduction to this).
//// 'do_branch_and_cut'  boolean (TRUE). Whether to run the branch and cut algorithm or not. (Set this to FALSE to run the user's heuristics only.)
//param = add_param(param, 'do_branch_and_cut', 0);
//
//// 'mc_search_order'  integer (MC_FIFO). Use the fifo (MC_FIFO) or lifo (MC_LIFO) searh order during the multi criteria solution procedure.
//param = add_param(param, 'mc_search_order', 0);
//
//// mc_warm_start  boolean (FALSE). Whether to solve the corresponding problem of each iteration from a warm start loaded from a base iteration 
////                        (which is the first iteration where gamma = 1.0 and tau = 0.0) or from scratch. Currently, this option is supported if only the supported
////                        solutions are desired to be found.
//param = add_param(param, 'mc_warm_start', 0);
//
//// 'trim_warm_tree'  boolean (FALSE). Whether to trim the warm start tree before re-solving. This consists of locating nodes whose descendants are all likely to be pruned 
////                           in the resolve and eliminating those descendants in favor of processing the parent node itself.
//param = add_param(param, 'trim_warm_tree', 0);
//
//// 'mc_compare_solution_tolerance'  double (0.001). If the difference between the objective values of two solutions to be compared, during the bicriteria solution 
////                                         procedure, are less than this tolerance, then assume them to be equal.
//param = add_param(param, 'mc_compare_solution_tolerance', 0.001);    
//
//// 'mc_binary_search_tolerance'  double (0). The tolerance to be used to differentiate the gamma values if binary search is used during the bicriteria solution procedure. 
////                                       A value greater than zero will cause the binary search to be activated.
//param = add_param(param, 'mc_binary_search_tolerance', 0);    
//
///////////////////////////////
//// Tree manager parameters //
///////////////////////////////
//// 'TM_verbosity'  integer (0). The verbosity of the TM module.
//param = add_param(param, 'TM_verbosity', 0);
//
//// 'lp_exe', 'cg_exe', 'cp_exe'  strings ('lp', 'cg', 'cp'). The name of the LP, CG, and CP module binaries. Note: when running in parallel using PVM, these executables 
////                                       (or links to them) must reside in the PVM ROOT/bin/PVM ARCH/ directory. Also, be sure to note that the executable names may have
////                                        extensions that depend on the configuration of the modules, but the defaults will always be set to the name that the makefile 
////                                        produces.
//param = add_param(param, 'lp_exe', 'lp');    
//param = add_param(param, 'cg_exe', 'cg');    
//param = add_param(param, 'cp_exe', 'cp');    
//
//// 'lp_debug', 'cg_debug', 'cp_debug'  boolean (all FALSE). Whether the modules should be started under a debugger or not.
//param = add_param(param, 'lp_debug', 0);
//param = add_param(param, 'cg_debug', 0);
//param = add_param(param, 'cp_debug', 0);
//
//// 'max_active_nodes'  integer (1). The maximum number of active search tree nodesequal to the number of LP and CG tandems to be started up.
// for misc07.mps:
// - if max_active_nodes == 4 then time needed to solve the problem = 247 sec on my amd64 quadcore
// - if max_active_nodes == 1 then time needed to solve the problem = 760 sec on my amd64 quadcore
////////////////////param = add_param(param, 'max_active_nodes', 1);
//
//// 'max_cp_num'  integer (0). The maximum number of cut pools to be used. BE CAREFUL: this option must be leaved unset in mono processor optimization
//param = add_param(param, 'max_cp_num', 0); 
//
//// 'lp_mach_num', 'cg_mach_num', 'cp_mach_num'  integers (all 0). The number of processors in the virtual machine to run LP (CG, CP) processes. 
////                                                       If this value is 0 then the processes will be assigned to processors in round-robin order. Otherwise
////                                                       the next xx mach num lines describe the processors where the LP (CG, CP) modules must run. The keyword  
////                                                       value pairs on these lines must be TM xx machine and the name or IP address of a processor (the processor
////                                                       names need not be distinct). In this case the actual processes are assigned in a round robin fashion to
////                                                       the processors on this list. This feature is useful if a specific software package is needed for some module, but
////                                                       that software is not licensed for every node of the virtual machine or if a certain process must run on a
////                                                       certain type of machine due to resource requirements.
//param = add_param(param, 'lp_mach_num', 0);
//param = add_param(param, 'cg_mach_num', 0);
//param = add_param(param, 'cp_mach_num', 0);
//
//// 'use_cg'  boolean (FALSE). Whether to use a cut generator or not.
param = add_param(param, 'use_cg', 1);
//
//// 'TM_random_seed'  integer (17). The random seed used in the TM.
//param = add_param(param, 'TM_random_seed', 17);
//
//// 'unconditional_dive_frac'  double (0.0). The fraction of the nodes on which SYMPHONY randomly dives unconditionally into one of the children.
//param = add_param(param, 'unconditional_dive_frac', 0.0);    
//
//// 'diving_strategy'  integer (BEST_ESTIMATE{0}). The strategy employed when deciding whether to dive or not. The BEST_ESTIMATE{0} strategy continues to dive until 
////                            the lower bound in the child to be dived into exceeds the parameter upper_bound_estimate, which is given by the user.
////                            The COMP_BEST_K{1} strategy computes the average lower bound on the best diving_k search tree nodes and decides to dive if
////                            the lower bound of the child to be dived into does not exceed this average by more than the fraction diving_threshold.
////                            The COMP_BEST_K_GAP{2} strategy takes the size of the gap into account when deciding whether to dive. After the average lower bound of
////                            the best diving_k nodes is computed, the gap between this average lower bound and the current upper bound is computed. Diving only occurs
////                            if the difference between the computed average lower bound and the lower bound of the child to be dived into is at most the fraction 
////                            diving_threshold of the gap.
////                            Note that fractional diving settings can override these strategies. See below.
//param = add_param(param, 'diving_strategy', 0);
//
//// 'diving_k', 'diving_threshold'  integer, double (1, 0.05). See above.
//param = add_param(param, 'diving_k', 1);
//param = add_param(param, 'diving_threshold', 0.05);    
//
//// 'fractional_diving_ratio', 'fractional_diving_num'  integer (0.02, 0). Diving occurs automatically if the number of fractional variables in the child to be dived into is 
////                                                         less than fractional diving num or the fraction of total variables that are fractional is less than 
////                                                         fractional_diving_ratio. This overrides the other diving rules. Note that in order for this option to work,
////                                                         the code must be compiled with FRACTIONAL_BRANCHING defined. This is the default. See the makefile
////                                                         for more details.
//param = add_param(param, 'fractional_diving_ratio', 0.02);    
//param = add_param(param, 'fractional_diving_num', 0);
//
//// 'node_selection_rule'  integer (LOWEST_LP_FIRST{0}). The rule for selecting the next search tree node to be processed. This rule selects the one with lowest lower bound. 
////                                                    Other possible values are: HIGHEST_LP_FIRST{1}, BREADTH_FIRST_SEARCH{2} and DEPTH_FIRST_SEARCH{3}.
//param = add_param(param, 'node_selection_rule', 0);
//
//// 'load_balance_level'  integer (-1) A naive attempt at load balancing on problems where significant time is spent in the root node, contributing to a lack of parallel 
////                             speed-up. Only a prescribed number of iterations (load balance iter) are performed in the root node (and in each subsequent node on a 
////                             level less than or equal to load_balance_level) before branching is forced in order to provide additional subproblems for the idle 
////                             processors to work on. This doesn't work well in general. 
//param = add_param(param, 'load_balance_level', -1);
//
//// 'load_balance_iter'  integer (-1) Works in tandem with the load_balance_level to attempt some simple load balancing. See the above description.
//param = add_param(param, 'load_balance_iter', -1);
//
//// 'keep_description_of_pruned'  integer (DISCARD{0}). Whether to keep the description of pruned search tree nodes or not. The reasons to do this are (1) if the user
////                                      wants to write out a proof of optimality using the logging function, (2) for debugging, or (3) to get a visual picture
////                                      of the tree using the software VBCTOOL. Otherwise, keeping the pruned nodes around just takes up memory.
////                                      There are three options if it is desired to keep some description of the pruned nodes around. First, their full description 
////                                      can be written out to disk and freed from memory (KEEP_ON_DISK_FULL{1}). There is not really too much you can do with this 
////                                      kind of file, but theoretically, it contains a full record of the solution process and could be used to provide a certificate 
////                                      of optimality (if we were using exact arithmetic) using an independent verifier.
////                                      In this case, the line following keep_description_of_pruned should be a line containing the keyword pruned_node_file_name
////                                      with its corresponding value being the name of a file to which a description of the pruned nodes can be written. The file 
////                                      does not need to exist and will be over-written if it does exist.
////                                      If you have the software VBCTOOL, then you can alternatively just write out the information VBCTOOL needs to display the tree
////                                      (KEEP_ON_DISK_VBC_TOOL{2}).
////                                      Finally, the user can set the value to of this parameter to KEEP_IN_MEMORY{2}, in which case all pruned nodes will be kept
////                                      in memory and written out to the regular log file if that option is chosen. This is really only useful for debugging. 
////                                      Otherwise, pruned nodes should be flushed.
//param = add_param(param, 'keep_description_of_pruned', 0);
//
//// 'keep_warm_start'  boolean (FALSE). Turning this parameter on will have exactly the same impact with setting the keep_description_of_pruned to KEEP_IN_MEMORY{2}. 
////                          This will allow SYMPHONY to keep all the necessary information obtained from the branching tree of the original problem to be able to 
////                          warm start after a parameter or problem data modification. Thus, if it is intended to warm start later, the user should set this parameter 
////                          before solving the original problem.
//param = add_param(param, 'keep_warm_start', 0);
//
//// 'warm_start_node_limit'  integer (SYM_INFINITY). Setting this parameter will start the warm start routine using only the first warm_start_node_limit nodes generated 
////                                during the previous solve procedure. The rest of the tree will be trimmed.
//param = add_param(param, 'warm_start_node_limit', 0);
//
//// 'warm_start_node_ratio'  double (0.0). Setting this parameter will start the warm start routine using only the first warm_start_node_ratio% of the nodes generated 
////                               during the previous solve procedure.
//param = add_param(param, 'warm_start_node_ratio', 0.0);    
//
//// 'warm_start_node_level'  integer (SYM_INFINITY = 100000000). Setting this parameter will start the warm start routine using all the nodes above the 
////                                  level warm_start_node_level of the tree generated during the previous solve procedure. The rest of the tree will be trimmed.
//param = add_param(param, 'warm_start_node_level', 100000000);
//
//// 'warm_start_node_level_ratio'  double (0.0). Setting this parameter will start the warm start routine using all the nodes above the level warm_start_node_level% 
////                                     of the warm start tree depth. The rest of the tree will be trimmed 
//param = add_param(param, 'warm_start_node_level_ratio', 0.0);    
//
//// 'logging'  integer (NO_LOGGING{0}). Whether or not to write out the state of the search tree and all other necessary data to disk periodically in order to allow 
////                  a warm start in the case of a system crash or to allow periodic viewing with VBCTOOL.
////                  If the value of this parameter is set to FULL_LOGGING{1}, then all information needed to warm start the calculation will written out periodically. 
////                  The next two lines of the parameter file following should contain the keywords tree_log_file_name and cut_log_file_name along with corresponding 
////                  file names as values. These will be the files used to record the search tree and related data and the list of cuts needed to reconstruct the tree.
////                  If the value of the parameter is set to VBC_TOOL{2}, then only the information VBCTOOL needs to display the tree will be logged. 
////                  This is not really a very useful option since a 'live' picture of the tree can be obtained using the vbc_emulation parameter described below.
//param = add_param(param, 'logging', 0);
//
//// 'logging_interval'  integer (1800). Interval (in seconds) between writing out the above log files.
//param = add_param(param, 'logging_interval', 1800);
//
//// 'warm_start' boolean (0). Used to allow the tree manager to make a warm start by reading in previously written log files. If this option is set, then the two 
////                    line following must start with the keywords warm_start_tree_file_name and warm_start_cut_file_name and include the appropriate file names
////                    as the corresponding values.
//param = add_param(param, 'warm_start', 0);
//
//// 'vbc_emulation'  integer (NO_VBC_EMULATION{0}) Determines whether or not to employ the VBCTOOL emulation mode. If one of these modes is chosen, then the tree will 
////                         be displayed in 'real time' using the VBCTOOL Software. When using the option VBC_EMULATION_LIVE{2} and piping the output directly to 
////                         VBCTOOL, the tree will be displayed as it is constructed, with color coding indicating the status of each node. 
////                         With VBC_EMULATION_FILE{1} selected, a log file will be produced which can later be read into VBCTOOL to produce an emulation of the
////                         solution process at any desired speed. If VBC_EMULATION_FILE is selected, the the following line should contain the keyword 
////                         vbc_emulation_file_name along with the corresponding file name for a value.
//param = add_param(param, 'vbc_emulation', 0);
//
//// 'price_in_root'  boolean (FALSE). Whether to price out variables in the root node before the second phase starts (called repricing the root).
//param = add_param(param, 'price_in_root', 0);
//
//// 'trim_search_tree'  boolean (FALSE). Whether to trim the search tree before the second phase starts or not. Useful only if there are two phases. (It is very useful then.)
//param = add_param(param, 'trim_search_tree', 0);
//
//// 'colgen_in_first_phase', 'colgen_in_second_phase'  integers (both 4). These parameters determine if and when to do column generation in the first and second phase of 
////                                                         the algorithm. The value of each parameter is obtained by setting the last four bits. The last two bits refer
////                                                         to what to do when attempting to prune a node. If neither of the last two bits are set, then we don't do
////                                                         anything we just prune it. If only the last bit is set, then we simply save the node for the second phase without 
////                                                         doing any column generation (yet). If only the second to last bit is set, then we do column generation 
////                                                         immediately and resolve if any new columns are found.
////                                                         The next two higher bits determine whether or not to do column generation before branching. If only the third
////                                                         lowest bit is set, then no column generation occurs before branching. If only the fourth lowest bit is set,
////                                                         then column generation is attempted before branching. The default is not to generate columns before branching 
////                                                         or fathoming, which corresponds to only the third lowest bit being set, resulting in a default value of 4.
//param = add_param(param, 'colgen_in_first_phase', 4);
//param = add_param(param, 'colgen_in_second_phase', 4);
//
//// 'time_limit'  double (-1.0). Number of seconds of wall-clock time allowed for solution. When this time limit is reached, the solution process will stop and the 
////                    best solution found to that point, along with other relevant data, will be output. A time limit less than 0.0 means there is no limit.
param = add_param(param, 'time_limit', 10000);
//
//// 'node_limit'  integer (-1). Number of nodes allowed to be analyzed during the solution. When this node limit is reached, the solution process will stop and the 
////                     best solution found to that point, along with other relevant data, will be output. A node limit less than 0 means there is no limit.
param = add_param(param, 'node_limit', 10000000);
//
//// 'gap_limit'  double (-1.0). Target gap limit allowed for solution. When the gap between the lower and the upper bound reaches this point, the solution process will 
////                   stop and the best solution found to that point, along with other relevant data, will be output. A gap limit less than 0 means there is no limit.
//param = add_param(param, 'gap_limit', -1.0);    
//
//// 'find_first_feasible'  boolean (FALSE). Whether to stop after finding the first feasible solution or not.
//param = add_param(param, 'find_first_feasible', 1);
//
//// 'sensitivity_analysis'  boolean (FALSE). If the user wants to do the rudimentary sensitivity analysis, which will give a lower bound for the problem modified by the
////                               right hand side, then, this parameter has to be set before solving the original problem. If it is set, SYMPHONY will keep the 
////                               necessary information from the solution processes of the original problem to be able to do the sensitivity analysis later.
//param = add_param(param, 'sensitivity_analysis', 0);
//
/////////////////////
//// LP Parameters //
/////////////////////
//// 'LP_verbosity'  integer (0). Verbosity level of the LP module.
param = add_param(param, 'LP_verbosity', 10);
//
//// 'set_obj_upper_lim'  boolean (FALSE). Whether to stop solving the LP relaxation when it's optimal value is provably higher than the global upper bound. 
////                              There are some advantages to continuing the solution process anyway. For instance, this results in the highest possible lower
////                              bound. On the other hand, if the matrix is full, this node will be pruned anyway and the rest of the computation is pointless. 
////                              This option should be set at FALSE for column generation since the LP dual values may not be reliable otherwise.
param = add_param(param, 'set_obj_upper_lim', 0);
//
//// 'try_to_recover_from_error'  boolean (TRUE). Indicates what should be done in case the LP solver is unable to solve a particular LP relaxation because of numerical 
////                                    problems. It is possible to recover from this situation but further results may be suspect. On the other hand, the entire solution
////                                    process can be abandoned.
//param = add_param(param, 'try_to_recover_from_error', 1);
//
//// 'problem_type'  integer (ZERO_ONE_PROBLEM{0}). The type of problem being solved. Other values are INTEGER_PROBLEM{1} or MIXED_INTEGER_PROBLEM{2}.
////                       (Caution: The mixed-integer option is not well tested.)
param = add_param(param, 'problem_type', 2);
//
//// 'cut_pool_check_frequency'  integer (10). The number of iterations between sending LP solutions to the cut pool to find violated cuts. It is not advisable to check 
////                                   the cut pool too frequently as the cut pool module can get bogged down and the LP solution generally do not change that 
////                                   drastically from one iteration to the next anyway.
//param = add_param(param, 'cut_pool_check_frequency', 10);
//
//// 'not_fixed_storage_size'  integer (2048). The not fixed list is a partial list of indices of variables not in the matrix that have not been fixed by reduced cost. 
////                                 Keeping this list allows SYMPHONY to avoid repricing variables (an expensive operation) that are not in the matrix because they
////                                 have already been permanently fixed. When this array reaches its maximum size, no more variable indices can be stored. It 
////                                 is therefore advisable to keep the maximum size of this array as large as possible, given memory limitations.
//param = add_param(param, 'not_fixed_storage_size', 2048);
//
//// 'max_non_dual_feas_to_add_min',
//// 'max_non_dual_feas_to_add_max',
//// 'max_non_dual_feas_to_add_frac'  integer, integer, double (20, 200, .05). These three parameters determine the maximum number of non-dual-feasible columns that can be 
////                                added in any one iteration after pricing. This maximum is set to the indicated fraction of the current number of active columns 
////                                unless this numbers exceeds the given maximum or is less than the given minimum, in which case, it is set to the max or min, respectively.
//param = add_param(param, 'max_non_dual_feas_to_add_min', 20);
//param = add_param(param, 'max_non_dual_feas_to_add_max', 200);
//param = add_param(param, 'max_non_dual_feas_to_add_frac', 0.05);    
//
//// 'max_not_fixable_to_add_min',
//// 'max_not_fixable_to_add_max',
//// 'max_not_fixable_to_add_frac'  integer, integer, double (100, 500, .1) As above, these three parameters determine the maximum number of new columns to be added to the
////                              problem because they cannot be priced out. These variables are only added when trying to restore infeasibility and usually, this does 
////                              not require many variables anyway.
//param = add_param(param, 'max_not_fixable_to_add_min', 100);
//param = add_param(param, 'max_not_fixable_to_add_max', 500);
//param = add_param(param, 'max_not_fixable_to_add_frac', 0.1);    
//
//// 'mat_col_compress_num', 'mat_col_compress_ratio'  integer, double (50, .05). Determines when the matrix should be physically compressed. This only happens when the 
////                                               number of columns is high enough to make it "worthwhile". The matrix is physically compressed when the number 
////                                               of deleted columns exceeds either an absolute number and a specified fraction of the current number of active columns.
//param = add_param(param, 'mat_col_compress_num', 50);
//param = add_param(param, 'mat_col_compress_ratio', 0.05);    
//
//// 'mat_row_compress_num', 'mat_row_compress_ratio'  integer, double (20, .05). Same as above except for rows.
//param = add_param(param, 'mat_row_compress_num', 20);
//param = add_param(param, 'mat_row_compress_ratio', 0.05);    
//
//// 'tailoff_gap_backsteps', 'tailoff_gap_frac'  integer, double (2, .99). Determines when tailoff is detected in the LP module. Tailoff is reported if the average 
////                                          ratio of the current gap to the previous iteration's gap over the last tailoff_gap_backsteps iterations wasn't at least
////                                          tailoff_gap_frac.
//param = add_param(param, 'tailoff_gap_backsteps', 2);
//param = add_param(param, 'tailoff_gap_frac', 0.99);    
//
//// 'tailoff_obj_backsteps', 'tailoff_obj_frac'  integer, double (2, .99). Same as above, only the ratio is taken with respect to the change in objective function values 
////                                          instead of the change in the gap.
//param = add_param(param, 'tailoff_obj_backsteps', 2);
//
//// 'ineff_cnt_to_delete'  integer (0). Determines after how many iterations of being deemed ineffective a constraint is removed from the current relaxation.
//param = add_param(param, 'ineff_cnt_to_delete', 0);
//
//// 'eff_cnt_before_cutpool'  integer (3). Determines after how many iterations of being deemed effective each cut will be sent to the global pool.
//param = add_param(param, 'eff_cnt_before_cutpool', 3);
//
//// 'ineffective_constraints'  integer (BASIC_SLACKS_ARE_INEFFECTIVE{2}). Determines under what condition a constraint is deemed ineffective in the current relaxation. 
////                          Other possible values are NO_CONSTRAINT_IS_INEFFECTIVE{0}, NONZERO_SLACKS_ARE_INEFFECTIVE{1}, and ZERO_DUAL_VALUES_ARE_INEFFECTIVE{3}.
//param = add_param(param, 'ineffective_constraints', 2);
//
//// 'base_constraints_always_effective'  boolean (TRUE). Determines whether the base constraints can ever be removed from the relaxation. In some case, removing the 
////                                    base constraints from the problem can be disastrous depending on the assumptions made by the cut generator. 
//param = add_param(param, 'base_constraints_always_effective', 1);
//
//// 'branch_on_cuts'  boolean (FALSE). This informs the framework whether the user plans on branching on cuts or not. If so, there is additional bookkeeping to be done, 
////                         such as maintaining a pool of slack cuts to be used for branching. Therefore, the user should not set this flag unless he actually plans
////                         on using this feature. 
param = add_param(param, 'branch_on_cuts', 1);
//
//// 'discard_slack_cuts'  integer (DISCARD_SLACKS_BEFORE_NEW_ITERATION{0}). Determines when the pool of slack cuts is discarded. The other option is 
////                             DISCARD_SLACKS_WHEN_STARTING_NEW_NODE{1}.
//param = add_param(param, 'discard_slack_cuts', 0);
//
//// 'first_lp_first_cut_time_out',
//// 'first_lp_all_cuts_time_out',
//// 'later_lp_first_cut_time_out',
//// 'later_lp_all_cuts_time_out'  double (0, 0, 5, 1). The next group of parameters determines when the LP should give up waiting for cuts from the cut generator and start 
////                                    to solve the relaxation in its current form or possibly branch if necessary. There are two factors that contribute to determining 
////                                    this timeout. First is whether this is the first LP in the search node of whether it is a later LP. Second is whether any cuts 
////                                    have been added already in this iteration. The four timeout parameters correspond to the four possible combinations of these
////                                    two variables.
//param = add_param(param, 'first_lp_first_cut_time_out', 0.0);    
//param = add_param(param, 'first_lp_all_cuts_time_out', 0.0); 
//param = add_param(param, 'later_lp_first_cut_time_out', 5.0);    
//param = add_param(param, 'later_lp_all_cuts_time_out', 1.0);    
//
//// 'no_cut_timeout'  This keyword does not have an associated value. If this keyword appears on a line by itself or with a value, this tells the framework not to time 
////                   out while waiting for cuts. This is useful for debugging since it enables runs with a single LP module to be duplicated.
//param = add_param(param, 'no_cut_timeout', 1);
//
//// 'all_cut_timeout'  double (no default). This keyword tells the framework to set all of the above timeout parameters to the value indicated.
//param = add_param(param, 'all_cut_timeout', 1.0);    
//
//// 'max_cut_num_per_iter'  integer (20). The maximum number of cuts that can be added to the LP in an iteration. The remaining cuts stay in the local pool to be 
////                                  added in subsequent iterations, if they are strong enough.
//param = add_param(param, 'max_cut_num_per_iter', 20);
//
//// 'do_reduced_cost_fixing'  boolean (FALSE). Whether or not to attempt to fix variables by reduced cost. This option is highly recommended
//param = add_param(param, 'do_reduce_cost_fixing', 0);
//
//// 'gap_as_ub_frac', 'gap_as_last_gap_frac'  double (.1, .7). Determines when reduced cost fixing should be attempted. It is only done when the gap is within the fraction 
////                                              gap as ub frac of the upper bound or when the gap has decreased by the fraction gap as last gap frac since the last
////                                              time variables were fixed.
//param = add_param(param, 'gap_as_ub_frac', 0.1);    
//param = add_param(param, 'gap_as_last_gap_frac', 0.7);    
//
//// 'do_logical_fixing'  boolean (FALSE). Determines whether the user's logical fixing routine should be used.
//param = add_param(param, 'do_logical_fixing', 0);
//
//// 'fixed_to_ub_before_logical_fixing',
//// 'fixed_to_ub_frac_before_logical_fixing'  integer, double (1, .01). Determines when logical fixing should be attempted. It will be called only when a certain
////                                         absolute number and a certain number of variables have been fixed to their upper bounds by reduced cost. This is because
////                                         it is typically only after fixing variables to their upper bound that other variables can be logically fixed.
//param = add_param(param, 'fixed_to_ub_before_logical_fixing', 1);
//param = add_param(param, 'fixed_to_ub_frac_before_logical_fixing', 0.01);    
//
//// 'max_presolve_iter'  integer (10). Number of simplex iterations to be performed in the presolve for strong branching.
//param = add_param(param, 'max_presolve_iter', 10);
//
//// 'strong_branching_cand_num_max',
//// 'strong_branching_cand_num_min',
//// 'strong_branching_red_ratio'  integer (10, 5, 1). These three parameters together determine the number of strong branching candidates to be used by default. 
////                             In the root node, strong_branching_cand_num_max candidates are used. On each succeeding level, this number is reduced by the
////                             number strong_branching_red_ratio multiplied by the square of the level. This continues until the number of candidates is reduced
////                             to strong_branching_cand_num_min and then that number of candidates is used in all lower levels of the tree.
//param = add_param(param, 'strong_branching_cand_num_max', 10);
//param = add_param(param, 'strong_branching_cand_num_min', 5);
//param = add_param(param, 'strong_branching_red_ratio', 1);
//
//// 'is_feasible_default'  integer (TEST_INTEGRALITY{1}). Determines the default test to be used to determine feasibility. This parameter is provided so that the
////                              user can change the default behavior without recompiling. The only other option is TEST_ZERO_ONE{0}. 
//param = add_param(param, 'is_feasible_default', 1);
//
//// 'send_feasible_solution_default'  integer (SEND_NONZEROS{0}). Determines the form in which to send the feasible solution. This parameter is provided so that the 
////                                         user can change the default behavior without recompiling. This is currently the only option.
//param = add_param(param, 'send_feasible_solution_default', 0);
//
//// 'send_lp_solution_default'  integer (SEND_NONZEROS{0}). Determines the default form in which to send the LP solution to the cut generator and cut pool. 
////                                   This parameter is provided so that the user can change the default behavior without recompiling. The other option is
////                                   SEND_FRACTIONS{1}.
//param = add_param(param, 'send_lp_solution_default', 0);
//
//// 'display_solution_default'  integer (DISP_NOTHING{0}). Determines how to display the current LP solution if desired. See the description of user display solution() 
////                                   for other possible values. This parameter is provided so that the user can change the default behavior without recompiling.
//param = add_param(param, 'display_solution_default', 0);
//
//// 'shall_we_branch_default'  integer (USER__BRANCH_IF_MUST{2}). Determines the default branching behavior. Other values are USER__DO_NOT_BRANCH{0} (not recommended
////                                  as a default), USER__DO_BRANCH{1} (also not recommended as a default), and USER__BRANCH_IF_TAILOFF{3}. This parameter 
////                                  is provided so that the user can change the default behavior without recompiling.
//param = add_param(param, 'shall_we_branch_default', 2);
//
//// 'select_candidates_default'  integer (USER__CLOSE_TO_HALF_AND_EXPENSIVE{10}). Determines the default rule for selecting strong branching candidates. Other values
////                                    are USER__CLOSE_TO_HALF{10} and USER__CLOSE_TO_ONE_AND_CHEAP{12}. This parameter is provided so that the user can change the
////                                    default behavior without recompiling.
//param = add_param(param, 'select_candidate_default', 10);
//
//// 'compare_candidates_default'  integer (HIGHEST_LOW_OBJ{2}). Determines the default rule for comparing candidates. See the description of user_compare_candidates() 
////                                     for other values. This parameter is provided so that the user can change the default behavior without recompiling.
//param = add_param(param, 'compare_candidates_default', 2);
//
//// 'select_child_default'  integer (PREFER_LOWER_OBJ_VALUE{0}). Determines the default rule for selecting the child to be processed next. For other possible values, 
////                               see the description user_select_child(). This parameter is provided so that the user can change the default behavior without recompiling.
//param = add_param(param, 'select_child_default', 0);
//
//// 'mc_find_supported_solutions'  boolean (FALSE). By default, sym_mc_solve routine will find all the non-dominated solutions if the problem to be solved is a 
////                                      bicriteria problem. However, if the user plans to find only the supported solutions, then, this parameter has to be set before
////                                      calling sym_mc_solve routine.
//param = add_param(param, 'mc_find_supported_solutions', 0);
//
//// 'mc_rho'  double (0.00001). The value used in augmented Chebyshev norm during the bicriteria solution procedure.
//param = add_param(param, 'mc_rho', 0.00001);    
//
//// 'generate_cgl_cuts'  boolean (TRUE). Whether or not to generate cuts using COIN's cut generation library. Note that, to use CGL cuts, OSI interface has to 
////                            be used and moreover the corresponding flags have to be set during installation. See the makefile for more details.
param = add_param(param, 'generate_cgl_cuts', 1);
//
//// 'generate_cgl_gomory_cuts'  boolean (TRUE). Whether or not to generate Gomory cuts using COIN's cut generation library.
param = add_param(param, 'generate_cgl_gomory_cuts', 1);
//
//// 'generate_cgl_knapsack_cuts'  boolean (TRUE). Whether or not to generate knapsack cover cuts using COIN's cut generation library.
param = add_param(param, 'generate_cgl_knapsack_cuts', 1);
//
//// 'generate_cgl_oddhole_cuts'  boolean (TRUE). Whether or not to generate generalized odd hole cuts using COIN's cut generation library.
param = add_param(param, 'generate_cgl_oddhole_cuts', 1);
//
//// 'generate_cgl_probing_cuts'  boolean (TRUE). Whether or not to generate probing cuts using COIN's cut generation library.
param = add_param(param, 'generate_cgl_probing_cuts', 1);
//
//// 'generate_cgl_clique_cuts'  boolean (TRUE). Whether or not to generate clique cuts using COIN's cut generation library.
param = add_param(param, 'generate_cgl_clique_cuts', 1);
//
//// 'generate_cgl_flow_and_cover_cuts'  boolean (FALSE). Whether or not to generate flow and cover cuts using COIN's cut generation library.
//param = add_param(param, 'generate_cgl_flow_and_cover_cuts', 0);
//
//// 'generate_cgl_rounding_cuts'  boolean (FALSE). Whether or not to generate simple rounding cuts using COIN's cut generation library.
param = add_param(param, 'generate_cgl_rounding_cuts', 1);
//
//// 'generate_cgl_lift_and_project_cuts'  boolean (FALSE). Whether or not to generate lift-andproject cuts using COIN's cut generation library.
param = add_param(param, 'generate_cgl_lift_and_project_cuts', 1);
//
////////////////////////////////
//// Cut generator parameters //
////////////////////////////////
//// 'CG_verbosity'  integer (0). Verbosity level for the cut generator module.
param = add_param(param, 'CG_verbosity', 0);
//
///////////////////////////
//// Cut pool parameters //
///////////////////////////
//// 'CP_verbosity'  integer (0). Verbosity of the cut pool module.
//param = add_param(param, 'CP_verbosity', 0);
//
//// 'cp_logging'  boolean (0). Determines whether the logging option is enabled. In this case, the entire contents of the cut pool are written out periodically to disk 
////                     (at the same interval as the tree manager log files are written). If this option is set, then the line following must start with the keyword 
////                     cp_log_file_name and include the appropriate file name as the value.
//param = add_param(param, 'cp_logging', 0);
//
//// 'cp_warm_start'  boolean (0). Used to allow the cut pool to make a warm start by reading in a previously written log file. If this option is set, then the line
////                        following must start with the keyword cp_warm_start_file_name and include the appropriate file name as the value.
//param = add_param(param, 'cp_warm_start', 0);
//
//// 'block_size'  integer (5000). Indicates the size of the blocks to allocate when more space is needed in the cut list.
//param = add_param(param, 'block_size', 5000);
//
//// 'max_size'  integer (2000000). Indicates the maximum size of the cut pool in bytes. This is the total memory taken up by the cut list, including all data
////                   structures and the array of pointers itself.
//param = add_param(param, 'max_size', 2000000);
//
//// 'max_number_of_cuts'  integer (10000). Indicates the maximum number of cuts allowed to be stored. When this max is reached, cuts are forcibly purged, 
////                             starting with duplicates and then those indicated by the parameter delete_which (see below), until the list is below the allowable size.
//param = add_param(param, 'max_number_of_cuts', 10000);
//
//// 'min_to_delete'  integer (1000). Indicates the number of cuts required to be deleted when the pool reaches it's maximum size.
//param = add_param(param, 'min_to_delete', 1000);
//
//// 'touches_until_deletion'  integer (10). When using the number of touches a cut has as a measure of its quality, this parameter indicates the number of touches a cut 
////                                 can have before being deleted from the pool. The number of touches is the number of times in a row that a cut has been checked
////                                 without being found to be violated. It is a measure of a cut's relevance or effectiveness.
//param = add_param(param, 'touches_until_deletion', 10);
//
//// 'delete_which'  integer (DELETE_BY_TOUCHES{2}). Indicates which cuts to delete when purging the pool. DELETE_BY_TOUCHES indicates that cuts whose number of touches 
////                       is above the threshold (see touches_until_deletion above) should be purged if the pool gets too large. DELETE_BY_QUALITY{1} indicates that
////                       a user-defined measure of quality should be used (see the function user_check_cuts in Section6.3.4).
//param = add_param(param, 'delete_which', 2);
//
//// 'check_which'  integer (CHECK_ALL_CUTS{0}). Indicates which cuts should be checked for violation. The choices are to check all cuts (CHECK_ALL_CUTS{0}); 
////                      only those that have number of touches below the threshold (CHECK_TOUCHES{2}); only those that were generated at a level higher in
////                      the tree than the current one (CHECK_LEVEL{1}); or both (CHECK_LEVEL_AND_TOUCHES{3}). Note that with CHECK_ALL_CUTS set, SYMPHONY will
////                      still only check the first cuts to check cuts in the list ordered by quality (see the function user_check_cut).
//param = add_param(param, 'check_which', 0);
//
//// 'cuts_to_check'  integer (1000). Indicates how many cuts in the pool to actually check. The list is ordered by quality and the first cuts to check cuts are 
////                        checked for violation.
//param = add_param(param, 'cuts_to_check', 1000);

///////////////////////////
// Set the variable type //
///////////////////////////

// 'I' -> integer

var_type = string(zeros(1,length(mps_file_mp('obj_var_is_int'))));
var_type(find(mps_file_mp('obj_var_is_int')==0)) = 'C';
var_type(find(mps_file_mp('obj_var_is_int')==1)) = 'I';
var_type = strcat(var_type);

/////////////////////////////
// Set the constraint type //
/////////////////////////////

// mps_file_mp('constr_sense')
// 'L' - smaller than - <=
// 'E' - equality     - =
// 'G' - greater than - >=
// 'R' - Range        - <= + >=
// 'N' - Free         - no constraints

[xmin,fmin,status,extra] = symphony(mps_file_mp('obj_coeff'),[],mps_file_mp('constr_mat'),mps_file_mp('lhs'),mps_file_mp('rhs'), ...
                                    mps_file_mp('bounds_lower'),mps_file_mp('bounds_upper'),mps_file_mp('constr_sense'),var_type,param);

printf('nb of constr  = %d\n', size(mps_file_mp('constr_mat'),1));
printf('nb of var     = %d\n', size(mps_file_mp('constr_mat'),2));
printf('nb of int var = %d\n', mps_file('nb_int_var'));
printf('status of problem: %d\n', status);
printf('Secondary status of problem - may get extended: %d \n', extra('sym_status'));
printf('iterations: %d\n', extra('iterations'));

t_end = getdate();

printf('elapsed time for solving problem = %f secondes\n', etime(t_end, t_start));

if UseLinpro then
  // 'L' - smaller than - <=
  // 'E' - equality     - =
  // 'G' - greater than - >=
  // 'R' - Range        - <= + >=
  // 'N' - Free         - no constraints
  index_F = find(ascii(mps_file_mp('constr_sense'))==ascii('N'));
  index_L = find(ascii(mps_file_mp('constr_sense'))==ascii('L'));
  index_U = find(ascii(mps_file_mp('constr_sense'))==ascii('U'));
  index_E = find(ascii(mps_file_mp('constr_sense'))==ascii('E'));
  constr_mat = mps_file_mp('constr_mat')(index_E,:);
  constr_mat = [constr_mat; -mps_file_mp('constr_mat')(index_L,:)];
  constr_mat = [constr_mat; mps_file_mp('constr_mat')(index_U,:)];
  constr_mat = [constr_mat; mps_file_mp('constr_mat')(index_F,:)];
  bound_mat  = mps_file_mp('rhs')(index_E)';
  bound_mat  = [bound_mat; -mps_file_mp('rhs')(index_L)'];
  bound_mat  = [bound_mat; mps_file_mp('rhs')(index_U)'];
  bound_mat  = [bound_mat; mps_file_mp('rhs')(index_F)'];

  t_start = t_end;

  n = length(mps_file_mp('obj_coeff'));
  
//  [xmin_linpro] = qpsolve(zeros(n,n),mps_file_mp('obj_coeff')',full(constr_mat), bound_mat, ...
//       		     		          mps_file_mp('bounds_lower')', mps_file_mp('bounds_upper')',length(index_E));
  [xmin_linpro,lagr,fmin] = linpro(mps_file_mp('obj_coeff')',full(constr_mat), bound_mat, ...
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
