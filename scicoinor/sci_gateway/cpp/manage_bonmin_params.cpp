///////////////////////////////////////////////////////////////////////////////////////////////////////////
// scibonmin: a scilab interface to Bonmin, a tool for non linear constrained mixed integer optimization //
///////////////////////////////////////////////////////////////////////////////////////////////////////////

//  Copyright (C) 2008-2010 Yann Collette.
//
//  SCIBONMIN is free software; you can redistribute it and/or modify it
//  under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 2, or (at your option) any
//  later version.
//
//  This part of code is distributed with the FURTHER condition that it 
//  can be compiled and linked with the Matlab libraries and it can be 
//  used within the Matlab environment.
//
//  SCIBONMIN is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with SciCoinOR; see the file COPYING.  If not, write to the Free
//  Software Foundation, 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <IpIpoptApplication.hpp>

#include <string>
#include <new>
#include <string>

extern "C"
{
#include <stack-c.h>
#include <sciprint.h>
#include <Scierror.h>
#include <MALLOC.h>
#include <api_parameters.h>
}

#include <api_scilab.h>

#include <def_solvers.hpp>

int manage_bonmin_params(Ipopt::SmartPtr<Ipopt::OptionsList> options, int * param_in_addr, int Log)
{
  // Get the parameters stored in the plist
  int     tmp_res, tmp_int;
  char *  tmp_char;
  double  tmp_double;

#ifdef DEBUG
  DBGPRINTF("scibonmin: processing options\n");
#endif

  // Some more parameters settable via :
  // - bonmin->setIntParameter(const IntParameter &p, const int v)
  // - bonmin->setDoubleParameter(const DoubleParameter &p, const double v)

  // IntParameters
  // BabLogLevel 	
  // BabLogInterval  Display information every logIntervval nodes.
  // MaxFailures     Max number of failures in a branch.
  // FailureBehavior Behavior of the algorithm in the case of a failure.
  // MaxInfeasible   Max number of consecutive infeasible problem in a branch before fathoming.
  // NumberStrong    Number of candidates for strong branching.
  // MinReliability  Minimum reliability before trust pseudo-costs.
  // MaxNodes        Global node limit.
  // MaxSolutions    limit on number of integer feasible solution.
  // MaxIterations   Global iteration limit.
  // SpecialOption   Spetial option in particular for Cbc.
  // DisableSos      Consider or not SOS constraints.
  // NumCutPasses    Number of cut passes at nodes.
  // NumCutPassesAtRoot Number of cut passes at nodes.

  // DoubleParameters
  // CutoffDecr 	
  // Cutoff               cutoff value
  // AllowableGap         Stop if absolute gap is less than this.
  // AllowableFractionGap Stop if relative gap is less than this.
  // IntTol               Integer tolerance.
  // MaxTime              Global time limit.

  ///////////////////
  // Set algorithm //
  ///////////////////
  // Possible choices:
  // B-BB  simple branch-and-bound algorithm
  // B-OA  OA Decomposition algorithm
  // B-QG  Quesada and Grossmann branch-and-cut algorithm
  // B-Hyb hybrid outer approximation based branch-and-cut
  // B-Ecp ecp cuts based branch-and-cut a la FilMINT

  getStringInPList(pvApiCtx, param_in_addr, "algorithm", &tmp_char, &tmp_res, "B-BB", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.algorithm", tmp_char);
  FREE(tmp_char);

  ////////////////////////////////////////
  //  Bonmin ecp based strong branching //
  ////////////////////////////////////////
  // Set the relative termination tolerance for ECP rounds in strong branching.
  getDoubleInPList(pvApiCtx, param_in_addr, "ecp_abs_tol_strong", &tmp_double, &tmp_res, 1e-6, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("bonmin.ecp_abs_tol_strong", tmp_double);
  // Set the absolute termination tolerance for ECP rounds in strong branching.
  getIntInPList(pvApiCtx, param_in_addr, "ecp_max_rounds_strong", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.ecp_max_rounds_strong", tmp_int);
  // Set the relative termination tolerance for ECP rounds in strong branching.
  getDoubleInPList(pvApiCtx, param_in_addr, "ecp_rel_tol_strong", &tmp_double, &tmp_res, 0.1, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("bonmin.ecp_rel_tol_strong", tmp_double);
  // Choose method to use for warm starting lp in strong branching:
  // - Basis: Use optimal basis of node
  // - Clone: Clone optimal problem of node (Advanced stuff)
  getStringInPList(pvApiCtx, param_in_addr, "lp_strong_warmstart_method", &tmp_char, &tmp_res, "Basis", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.lp_strong_warmstart_method", tmp_char);
  FREE(tmp_char);

  ///////////////////////////////
  //  Branch-and-bound options //
  ///////////////////////////////
  // Specify the value of relative gap under which the algorithm stops.
  // Stop the tree search when the gap between the objective value of the best known solution and the best bound on the objective of any solution is less than this
  // fraction of the absolute value of the best known solution value.
  getDoubleInPList(pvApiCtx, param_in_addr, "allowable_fraction_gap", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("bonmin.allowable_fraction_gap", tmp_double);
  // Specify the value of absolute gap under which the algorithm stops.
  // Stop the tree search when the gap between the objective value of the best known solution and the best bound on the objective of any solution is less than this.
  getDoubleInPList(pvApiCtx, param_in_addr, "allowable_gap", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("bonmin.allowable_gap", tmp_double);
  // Specify cutoff value.
  // cutoff should be the value of a feasible solution known by the user (if any). The algorithm will only look for solutions better than cutoof.
  getDoubleInPList(pvApiCtx, param_in_addr, "cutoff", &tmp_double, &tmp_res, 1e100, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("bonmin.cutoff", tmp_double);
  // Specify cutoff decrement.
  // Specify the amount by which cutoff is decremented below a new best upper-bound (usually a small positive value but in non-convex problems it may be a negative value).
  getDoubleInPList(pvApiCtx, param_in_addr, "cutoff_decr", &tmp_double, &tmp_res, 1e-5, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("bonmin.cutoff_decr", tmp_double);
  // Set integer tolerance
  getDoubleInPList(pvApiCtx, param_in_addr, "integer_tolerance", &tmp_double, &tmp_res, 1e-6, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("bonmin.integer_tolerance", tmp_double);
  // Set the cumulated maximum number of iteration in the algorithm used to process nodes continuous relaxations in the branch-and-bound
  getIntInPList(pvApiCtx, param_in_addr, "iteration_limit", &tmp_int, &tmp_res, INT_MAX, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.iteration_limit", tmp_int);
  // Set the behavior when an NLP or a series of NLP are unsolved by Ipopt (we call unsolved an NLP for which Ipopt is not able to guarantee
  // optimality within the specified tolerances).
  // - stop: Stop when failure happens.
  // - fathom: Continue when failure happens.
  // If set to "fathom", the algorithm will fathom the node when Ipopt fails to find a solution to the nlp at that node whithin the specified tolerances.
  // The algorithm then becomes a heuristic, and the user will be warned that the solution might not be optimal.
  getStringInPList(pvApiCtx, param_in_addr, "nlp_failure_behavior", &tmp_char, &tmp_res, "Stop", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.nlp_failure_behavior", tmp_char);
  FREE(tmp_char);
  // Choose the node selection strategy.
  // - best-bound: choose node with the smallest bound
  // - depth-first: Perform depth first search
  // - breadth-first: Perform breadth first search
  // - dynamic: Cbc dynamic strategy (starts with a depth first search and turn to best bound after 3 integer feasible solutions have been found)
  // - best-guess: choose node with smallest guessed integer solution
  getStringInPList(pvApiCtx, param_in_addr, "node_comparison", &tmp_char, &tmp_res, "dynamic", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.node_comparison", tmp_char);
  FREE(tmp_char);
  // Set the maximum number of nodes explored in the branch-and-bound search
  getIntInPList(pvApiCtx, param_in_addr, "node_limit", &tmp_int, &tmp_res, INT_MAX, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.node_limit", tmp_int);
  // Set the maximum number of cut passes at regular nodes of the branch-and-cut.
  getIntInPList(pvApiCtx, param_in_addr, "num_cut_passes", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.num_cut_passes", tmp_int);
  // Set the maximum number of cut passes at regular nodes of the branch-and-cut.
  getIntInPList(pvApiCtx, param_in_addr, "num_cut_passes_at_root", &tmp_int, &tmp_res, 20, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.num_cut_passes_at_root", tmp_int);
  // Set the number of branches on a variable before its pseudo costs are to be believed in dynamic strong branching.
  getIntInPList(pvApiCtx, param_in_addr, "number_before_trust", &tmp_int, &tmp_res, 8, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.number_before_trust", tmp_int);
  // Choose the maximum number of variables considered for strong branching.
  getIntInPList(pvApiCtx, param_in_addr, "number_strong_branch", &tmp_int, &tmp_res, 20, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.number_strong_branch", tmp_int);
  // Abort after that much integer feasible solution have been found by algorithm
  getIntInPList(pvApiCtx, param_in_addr, "solution_limit", &tmp_int, &tmp_res, INT_MAX, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.solution_limit", tmp_int);
  // Wether or not to activate SOS constraints.
  // - enable
  // - disable
  getStringInPList(pvApiCtx, param_in_addr, "sos_constraints", &tmp_char, &tmp_res, "enable", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.sos_constraints", tmp_char);
  FREE(tmp_char);
  // Set the global maximum computation time (in secs) for the algorithm.
  getDoubleInPList(pvApiCtx, param_in_addr, "time_limit", &tmp_double, &tmp_res, 1e10, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("bonmin.time_limit", tmp_double);
  // Pick a strategy for traversing the tree
  // - probed-dive:
  // - top-node: Always pick the top node as sorted by the node comparison function
  // - dive: Dive in the tree if possible, otherwise pick top node as sorted by the tree comparison function.
  // - probed-dive: Dive in the tree exploring two childs before continuing the dive at each level.
  // - dfs-dive: Dive in the tree if possible doing a depth first search.
  // Backtrack on leaves or when a prescribed depth is attained or when estimate of best possible integer feasible solution in subtree
  // is worst than cutoff. Once a prescribed limit of backtracks is attained pick top node as sorted by the tree comparison function
  // - dfs-dive-dynamic: Same as dfs-dive but once enough solution are found switch to best-bound and if too many nodes switch to depth-first.
  // All strategies can be used in conjunction with any of the node comparison functions. Options which affect dfs-dive are max-backtracks-in-dive and max-dive-depth.
  // The dfs-dive won't work in a non-convex problem where objective does not decrease down branches.
  getStringInPList(pvApiCtx, param_in_addr, "tree_search_strategy", &tmp_char, &tmp_res, "top-node", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.tree_search_strategy", tmp_char);
  FREE(tmp_char);
  // Chooses variable selection strategy
  // - most-fractional: Choose most fractional variable
  // - strong-branching: Perform strong branching
  // - reliability-branching: Use reliability branching
  // - curvature-estimator: Use curvature estimation to select branching variable
  // - qp-strong-branching: Perform strong branching with QP approximation
  // - lp-strong-branching: Perform strong branching with LP approximation
  // - nlp-strong-branching: Perform strong branching with NLP approximation
  // - osi-simple: Osi method to do simple branching
  // - osi-strong: Osi method to do strong branching
  // - random: Method to choose branching variable randomly
  getStringInPList(pvApiCtx, param_in_addr, "variable_selection", &tmp_char, &tmp_res, "strong-branching", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.variable_selection", tmp_char);
  FREE(tmp_char);

  /////////////////////
  //  Diving options //
  /////////////////////
  // Set the number of backtracks in a dive when using dfs-dive tree search strategy.
  getIntInPList(pvApiCtx, param_in_addr, "max_backtracks_in_dive", &tmp_int, &tmp_res, 5, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.max_backtracks_in_dive", tmp_int);
  // When using dfs-dive search. Maximum depth to go to from the diving board (node where the diving started.
  getIntInPList(pvApiCtx, param_in_addr, "max_dive_depth", &tmp_int, &tmp_res, INT_MAX, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.max_dive_depth", tmp_int);
  // Flag indicating whether we stop diving based on guessed feasible objective and the current cutoff
  // - no
  // - yes
  getStringInPList(pvApiCtx, param_in_addr, "stop_diving_on_cutoff", &tmp_char, &tmp_res, "no", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.stop_diving_on_cutoff", tmp_char);
  FREE(tmp_char);

  /////////////////////////////////
  //  Feasibility pump heuristic //
  /////////////////////////////////
  getIntInPList(pvApiCtx, param_in_addr, "feasibility_pump_objective_norm", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.feasibility_pump_objective_norm", tmp_int);
  getStringInPList(pvApiCtx, param_in_addr, "heuristic_feasibility_pump", &tmp_char, &tmp_res, "no", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.heuristic_feasibility_pump", tmp_char);
  FREE(tmp_char);

  //////////////////////////////////////
  //  Fractional diving MIP heuristic //
  //////////////////////////////////////
  getStringInPList(pvApiCtx, param_in_addr, "heuristic_dive_MIP_fractional", &tmp_char, &tmp_res, "no", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.heuristic_dive_MIP_fractional", tmp_char);
  FREE(tmp_char);

  //////////////////////////////////
  //  Fractional diving heuristic //
  //////////////////////////////////
  getStringInPList(pvApiCtx, param_in_addr, "heuristic_dive_fractional", &tmp_char, &tmp_res, "no", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.heuristic_dive_fractional", tmp_char);
  FREE(tmp_char);

  ////////////////////////////////////
  //  Local search based heuristics //
  ////////////////////////////////////
  getStringInPList(pvApiCtx, param_in_addr, "dummy_pump_heuristic", &tmp_char, &tmp_res, "no", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.dummy_pump_heuristic", tmp_char);
  FREE(tmp_char);
  getStringInPList(pvApiCtx, param_in_addr, "fix_and_solve_heuristic", &tmp_char, &tmp_res, "no", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.fix_and_solve_heuristic", tmp_char);
  FREE(tmp_char);
  getStringInPList(pvApiCtx, param_in_addr, "heuristic_RINS", &tmp_char, &tmp_res, "no", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.heuristic_RINS", tmp_char);
  FREE(tmp_char);
  getStringInPList(pvApiCtx, param_in_addr, "heuristic_local_branching", &tmp_char, &tmp_res, "no", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.heuristic_local_branching", tmp_char);
  FREE(tmp_char);

  ////////////////////////////////////
  //  MILP cutting planes in hybrid //
  ////////////////////////////////////
  // Frequency k (in terms of nodes) for generating Gomory cuts in branch-and-cut.
  // If k > 0, cuts are generated every k nodes, if -99 < k < 0 cuts are generated every -k nodes but 
  // Cbc may decide to stop generating cuts, if not enough are generated at the root node
  // if k=-99 generate cuts only at the root node, if k=0 or 100 do not generate cuts
  getIntInPList(pvApiCtx, param_in_addr, "2mir_cuts", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.2mir_cuts", tmp_int);
  // Frequency k (in terms of nodes) for generating Gomory cuts in branch-and-cut.
  getIntInPList(pvApiCtx, param_in_addr, "Gomory_cuts", &tmp_int, &tmp_res, -5, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.Gomory_cuts", tmp_int);
  // Frequency (in terms of nodes) for generating clique cuts in branch-and-cut
  getIntInPList(pvApiCtx, param_in_addr, "clique_cuts", &tmp_int, &tmp_res, -5, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.clique_cuts", tmp_int);
  // Frequency (in terms of nodes) for generating cover cuts in branch-and-cut
  getIntInPList(pvApiCtx, param_in_addr, "cover_cuts", &tmp_int, &tmp_res, -5, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.cover_cuts", tmp_int);
  // Frequency (in terms of nodes) for generating flow cover cuts in branch-and-cut
  getIntInPList(pvApiCtx, param_in_addr, "flow_covers_cuts", &tmp_int, &tmp_res, -5, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.flow_covers_cuts", tmp_int);
  // Frequency (in terms of nodes) for generating lift-and-project cuts in branch-and-cut
  getIntInPList(pvApiCtx, param_in_addr, "lift_and_project_cuts", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.lift_and_project_cuts", tmp_int);
  // Frequency (in terms of nodes) for generating MIR cuts in branch-and-cut
  getIntInPList(pvApiCtx, param_in_addr, "mir_cuts", &tmp_int, &tmp_res, -5, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.mir_cuts", tmp_int);
  // Frequency (in terms of nodes) for generating probing cuts in branch-and-cut
  getIntInPList(pvApiCtx, param_in_addr, "probing_cuts", &tmp_int, &tmp_res, -5, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.probing_cuts", tmp_int);
  // Frequency (in terms of nodes) for generating reduce-and-split cuts in branch-and-cut
  getIntInPList(pvApiCtx, param_in_addr, "reduce_and_split_cuts", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.reduce_and_split_cuts", tmp_int);

  //////////////////////////////
  //  Nlp solution robustness //
  //////////////////////////////
  // (temporarily removed) Number $n$ of consecutive unsolved problems before aborting a branch of the tree.
  getIntInPList(pvApiCtx, param_in_addr, "max_consecutive_failures", &tmp_int, &tmp_res, 10, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.max_consecutive_failures", tmp_int);
  // Set max value r for coordinate of a random point.
  getDoubleInPList(pvApiCtx, param_in_addr, "max_random_point_radius", &tmp_double, &tmp_res, 100000, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("bonmin.max_random_point_radius", tmp_double);
  // Number of iterations over which a node is considered "suspect" (for debugging purposes only, see detailed documentation).
  getIntInPList(pvApiCtx, param_in_addr, "num_iterations_suspect", &tmp_int, &tmp_res, -1, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.num_iterations_suspect", tmp_int);
  // Number $k$ of times that the algorithm will try to resolve an unsolved NLP with a random starting point (we call unsolved an NLP for which
  // Ipopt is not able to guarantee optimality within the specified tolerances).
  getIntInPList(pvApiCtx, param_in_addr, "num_retry_unsolved_random_point", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.num_retry_unsolved_random_point", tmp_int);
  // Amount by which starting point is perturbed when choosing to pick random point by perturbating starting point
  getIntInPList(pvApiCtx, param_in_addr, "random_point_perturbation_interval", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("bonmin.random_point_perturbation_interval", tmp_int);
  // method to choose a random starting point
  getStringInPList(pvApiCtx, param_in_addr, "random_point_type", &tmp_char, &tmp_res, "Jon", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.random_point_type", tmp_char);
  FREE(tmp_char);

  /////////////////////////////////
  //  Nlp solve options in B-Hyb //
  /////////////////////////////////
  // Specify the frequency (in terms of nodes) at which NLP relaxations are solved in B-Hyb.
  getIntInPList(pvApiCtx, param_in_addr, "nlp_solve_frequency", &tmp_int, &tmp_res, 10, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.nlp_solve_frequency", tmp_int);
  // Set maximum depth in the tree at which NLP relaxations are solved in B-Hyb.
  getIntInPList(pvApiCtx, param_in_addr, "nlp_solve_max_depth", &tmp_int, &tmp_res, 10, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.nlp_solve_max_depth", tmp_int);
  // Set average number of nodes in the tree at which NLP relaxations are solved in B-Hyb for each depth.
  getDoubleInPList(pvApiCtx, param_in_addr, "nlp_solves_per_depth", &tmp_double, &tmp_res, 1e30, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("bonmin.nlp_solves_per_depth", tmp_double);

  /////////////////////////////////////////////////////
  //  Options for MILP subsolver in OA decomposition //
  /////////////////////////////////////////////////////
  // specify MILP subsolver log level.
  getIntInPList(pvApiCtx, param_in_addr, "milp_log_level", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.milp_log_level", tmp_int);
  // Choose the solver to solve MILP sub-problems in OA decompositions.
  getStringInPList(pvApiCtx, param_in_addr, "milp_solver", &tmp_char, &tmp_res, "Cbc_D", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.milp_solver", tmp_char);
  FREE(tmp_char);

  ///////////////////////////////////
  //  Options for OA decomposition //
  ///////////////////////////////////
  // display an update on lower and upper bounds in OA every n seconds
  getDoubleInPList(pvApiCtx, param_in_addr, "oa_log_frequency", &tmp_double, &tmp_res, 100, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("bonmin.oa_log_frequency", tmp_double);
  // specify OA iterations log level.
  getIntInPList(pvApiCtx, param_in_addr, "oa_log_level", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.oa_log_level", tmp_int);

  //////////////////////////////////////
  //  Options for ecp cuts generation //
  //////////////////////////////////////
  // Set the absolute termination tolerance for ECP rounds.
  getDoubleInPList(pvApiCtx, param_in_addr, "ecp_abs_tol", &tmp_double, &tmp_res, 1e-6, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("bonmin.ecp_abs_tol", tmp_double);
  // Set the maximal number of rounds of ECP cuts.
  getIntInPList(pvApiCtx, param_in_addr, "ecp_max_rounds", &tmp_int, &tmp_res, 5, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.ecp_max_rounds", tmp_int);
  // Factor appearing in formula for skipping ECP cuts.
  getDoubleInPList(pvApiCtx, param_in_addr, "ecp_propability_factor", &tmp_double, &tmp_res, 1000, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("bonmin.ecp_propability_factor", tmp_double);
  // Set the relative termination tolerance for ECP rounds.
  getDoubleInPList(pvApiCtx, param_in_addr, "ecp_rel_tol", &tmp_double, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("bonmin.ecp_rel_tol", tmp_double);
  // Specify the frequency (in terms of nodes) at which some a la filmint ecp cuts are generated.
  getIntInPList(pvApiCtx, param_in_addr, "filmint_ecp_cuts", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.filmint_ecp_cuts", tmp_int);

  //////////////////////////////////////
  //  Options for non-convex problems //
  //////////////////////////////////////
  // Number of consecutive infeasible subproblems before aborting a branch.
  getIntInPList(pvApiCtx, param_in_addr, "max_consecutive_infeasible", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.max_consecutive_infeasible", tmp_int);
  // Number $k$ of tries to resolve an infeasible node (other than the root) of the tree with different starting point.
  getIntInPList(pvApiCtx, param_in_addr, "num_resolve_at_infeasibles", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.num_resolve_at_infeasibles", tmp_int);
  // Number $k$ of tries to resolve a node (other than the root) of the tree with different starting point.
  getIntInPList(pvApiCtx, param_in_addr, "num_resolve_at_node", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.num_resolve_at_node", tmp_int);
  // Number $k$ of tries to resolve the root node with different starting points.
  getIntInPList(pvApiCtx, param_in_addr, "num_resolve_at_root", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.num_resolve_at_root", tmp_int);

  //////////////////////////////////////////
  //  Outer Approximation cuts generation //
  //////////////////////////////////////////
  // Do we add all OA cuts or only the ones violated by current point?
  getStringInPList(pvApiCtx, param_in_addr, "add_only_violated_oa", &tmp_char, &tmp_res, "no", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.add_only_violated_oa", tmp_char);
  FREE(tmp_char);
  // Determines if and what kind of cut strengthening should be performed.
  getStringInPList(pvApiCtx, param_in_addr, "cut_strengthening_type", &tmp_char, &tmp_res, "none", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.cut_strengthening_type", tmp_char);
  FREE(tmp_char);
  // Determine if and what kind of disjunctive cuts should be computed.
  getStringInPList(pvApiCtx, param_in_addr, "disjunctive_cut_type", &tmp_char, &tmp_res, "none", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.disjunctive_cut_type", tmp_char);
  FREE(tmp_char);
  // level of log when generating OA cuts.
  getIntInPList(pvApiCtx, param_in_addr, "oa_cuts_log_level", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.oa_cuts_log_level", tmp_int);
  // Specify if OA cuts added are to be set globally or locally valid
  getStringInPList(pvApiCtx, param_in_addr, "oa_cuts_scope", &tmp_char, &tmp_res, "global", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.oa_cuts_scope", tmp_char);
  FREE(tmp_char);
  // Value for tiny element in OA cut
  getDoubleInPList(pvApiCtx, param_in_addr, "tiny_element", &tmp_double, &tmp_res, 1e-8, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("bonmin.tiny_element", tmp_double);
  // Value for very tiny element in OA cut
  getDoubleInPList(pvApiCtx, param_in_addr, "very_tiny_element", &tmp_double, &tmp_res, 1e-17, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("bonmin.very_tiny_element", tmp_double);

  ////////////////////////////////////
  //  Output ond log-levels options //
  ////////////////////////////////////
  // bb_log_interval description:
  // Set the interval (in terms of number of nodes) at which a log on node resolutions (consisting of lower and upper bounds) is given.
  getIntInPList(pvApiCtx, param_in_addr, "bb_log_interval", &tmp_int, &tmp_res, 100, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.bb_log_interval",tmp_int );
  // bb_log_level description:
  // Set the level of output of the branch-and-bound :
  // 0 - none, 1 - minimal, 2 - normal low, 3 - normal high
  getIntInPList(pvApiCtx, param_in_addr, "bb_log_level", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.bb_log_level", tmp_int);
  // lp_log_level description:
  // Set the level of output of the linear programming sub-solver in B-Hyb or B-QG :
  // 0 - none, 1 - minimal, 2 - normal low, 3 - normal high, 4 - verbose
  getIntInPList(pvApiCtx, param_in_addr, "lp_log_level", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.lp_log_level", tmp_int);

  /////////////////////////////
  //  Strong branching setup //
  /////////////////////////////
  // Choice of the criterion to choose candidates in strong-branching
  getStringInPList(pvApiCtx, param_in_addr, "candidate_sort_criterion", &tmp_char, &tmp_res, "best-ps-cost", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.candidate_sort_criterion", tmp_char);
  FREE(tmp_char);
  //  Weight towards minimum in of lower and upper branching estimates when a solution has been found.
  getDoubleInPList(pvApiCtx, param_in_addr, "maxmin_crit_have_sol", &tmp_double, &tmp_res, 0.1, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("bonmin.maxmin_crit_have_sol", tmp_double);
  // Weight towards minimum in of lower and upper branching estimates when no solution has been found yet.
  getDoubleInPList(pvApiCtx, param_in_addr, "maxmin_crit_no_sol", &tmp_double, &tmp_res, 0.7, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("bonmin.maxmin_crit_no_sol", tmp_double);
  // Sets minimum number of variables for strong branching (overriding trust)
  getIntInPList(pvApiCtx, param_in_addr, "min_number_strong_branch", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.min_number_strong_branch", tmp_int);
  // Set the number of branches on a variable before its pseudo costs are to be believed during setup of strong branching candidate list.
  getIntInPList(pvApiCtx, param_in_addr, "number_before_trust_list", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.number_before_trust_list", tmp_int);
  // Sets limit of look-ahead strong-branching trials
  getIntInPList(pvApiCtx, param_in_addr, "number_look_ahead", &tmp_int, &tmp_res, 0, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.number_look_ahead", tmp_int);
  // Maximum number of variables considered for strong branching in root node.
  getIntInPList(pvApiCtx, param_in_addr, "number_strong_branch_root", &tmp_int, &tmp_res, INT_MAX, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.number_strong_branch_root", tmp_int);
  // Proportion of strong branching list that has to be taken from most-integer-infeasible list.
  getDoubleInPList(pvApiCtx, param_in_addr, "setup_pseudo_frac", &tmp_double, &tmp_res, 0.5, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("bonmin.setup_pseudo_frac", tmp_double);
  // Wether or not to trust strong branching results for updating pseudo costs.
  getStringInPList(pvApiCtx, param_in_addr, "trust_strong_branching_for_pseudo_cost", &tmp_char, &tmp_res, "yes", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.trust_strong_branching_for_pseudo_cost", tmp_char);
  FREE(tmp_char);

  ////////////////////////////////////////
  //  VectorLength diving MIP heuristic //
  ////////////////////////////////////////
  getStringInPList(pvApiCtx, param_in_addr, "heuristic_dive_MIP_vectorLength", &tmp_char, &tmp_res, "no", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.heuristic_dive_MIP_vectorLength", tmp_char);
  FREE(tmp_char);

  ////////////////////////////////////
  //  VectorLength diving heuristic //
  ////////////////////////////////////
  getStringInPList(pvApiCtx, param_in_addr, "heuristic_dive_vectorLength", &tmp_char, &tmp_res, "no", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.heuristic_dive_vectorLength", tmp_char);
  FREE(tmp_char);

  ///////////////////////////
  //  nlp interface option //
  ///////////////////////////
  getStringInPList(pvApiCtx, param_in_addr, "file_solution", &tmp_char, &tmp_res, "no", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.file_solution", tmp_char);
  FREE(tmp_char);
  // specify NLP solver interface log level (independent from ipopt print_level).
  getIntInPList(pvApiCtx, param_in_addr, "nlp_log_level", &tmp_int, &tmp_res, 1, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.nlp_log_level", tmp_int);
  // Choice of the solver for local optima of continuous nlp's
  getStringInPList(pvApiCtx, param_in_addr, "nlp_solver", &tmp_char, &tmp_res, "Ipopt", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.nlp_solver", tmp_char);
  FREE(tmp_char);
  // Select the warm start method
  getStringInPList(pvApiCtx, param_in_addr, "warm_start", &tmp_char, &tmp_res, "none", Log, CHECK_NONE);
  if (tmp_res!=1) options->SetStringValue("bonmin.warm_start", tmp_char);
  FREE(tmp_char);

  ///////////////////
  // Other options //
  ///////////////////

  // Coefficient of variation threshold (for dynamic definition of cutoff_decr).
  getDoubleInPList(pvApiCtx, param_in_addr, "coeff_var_threshold", &tmp_double, &tmp_res, 0.1, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("bonmin.coeff_var_threshold", tmp_double);
  // Do you want to define the parameter cutoff_decr dynamically?
  getStringInPList(pvApiCtx, param_in_addr, "dynamic_def_cutoff_decr", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=1) options->SetStringValue("bonmin.dynamic_cutoff_decr", tmp_char);
  FREE(tmp_char);
  // Enable dynamic linear and quadratic rows addition in nlp
  getStringInPList(pvApiCtx, param_in_addr, "enable_dynamic_nlp", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=1) options->SetStringValue("bonmin.enable_dynamic_nlp", tmp_char);
  FREE(tmp_char);
  // Choose the type of cuts generated when an integer feasible solution is found
  getStringInPList(pvApiCtx, param_in_addr, "feas_check_cut_types", &tmp_char, &tmp_res, "outer-approximations", Log, 
		       CHECK_VALUES, 2, 
		       "outer-approximations", 
		       "Benders");
  if (tmp_res!=1) options->SetStringValue("bonmin.feas_check_cut_type", tmp_char);
  FREE(tmp_char);
  // How cuts from feasibility checker are discarded
  getStringInPList(pvApiCtx, param_in_addr, "feas_check_discard_policy", &tmp_char, &tmp_res, "detect-cycles", Log, 
		       CHECK_VALUES, 3, 
		       "detect-cycles", 
		       "keep-all",
		       "treated-as-normal");
  if (tmp_res!=1) options->SetStringValue("bonmin.feas_check_discard_policy", tmp_char);
  FREE(tmp_char);
  // The percentage used when, the coeff of variance is smaller than the threshold, to compute the cutoff_decr dynamically.
  getDoubleInPList(pvApiCtx, param_in_addr, "first_perc_for_cutoff_decr", &tmp_double, &tmp_res, -0.02, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("bonmin.first_perc_for_cutoff_decr", tmp_double);
  // display an update on lower and upper bounds in FP every n seconds
  getDoubleInPList(pvApiCtx, param_in_addr, "fp_log_frequency", &tmp_double, &tmp_res, 1.0, Log, 
		    CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetNumericValue("bonmin.fp_log_frequency", tmp_double);
  // specify FP iterations log level.
  getIntInPList(pvApiCtx, param_in_addr, "fp_log_level", &tmp_int, &tmp_res, 0, Log, 
		    CHECK_VALUES, 3, 0, 1, 2);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.fp_log_level", tmp_int);
  // Specify that after so many oa cuts have been generated Benders cuts should be generated instead.
  getIntInPList(pvApiCtx, param_in_addr, "generate_benders_after_so_many_oa", &tmp_int, &tmp_res, 5000, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.generate_benders_after_so_many_oa", tmp_int);
  // Specify a different log level for root relaxtion.
  getIntInPList(pvApiCtx, param_in_addr, "nlp_log_at_root", &tmp_int, &tmp_res, 0, Log, CHECK_MIN, 0);
  if (tmp_res!=-1) options->SetIntegerValue("bonmin.nlp_log_at_root", tmp_int);
  // If yes do initial OA decomposition
  getStringInPList(pvApiCtx, param_in_addr, "oa_decomposition", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=1) options->SetStringValue("bonmin.oa_decomposition", tmp_char);
  FREE(tmp_char);
  // if yes runs FP for MINLP
  getStringInPList(pvApiCtx, param_in_addr, "pump_for_minlp", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");
  if (tmp_res!=1) options->SetStringValue("bonmin.pump_for_minlp", tmp_char);
  FREE(tmp_char);
  // The percentage used when, the coeff of variance is greater than the threshold, to compute the cutoff_decr dynamically.
  getDoubleInPList(pvApiCtx, param_in_addr, "second_perc_for_cutoff_decr", &tmp_double, &tmp_res, -0.05, Log, CHECK_NONE);
  if (tmp_res!=-1) options->SetNumericValue("bonmin.second_perc_for_cutoff_decr", tmp_double);

#ifdef DEBUG
  options->SetStringValue("print_options_documentation", "yes");
  bonmin->mayPrintDoc();
#endif

  // print_user_options: If selected, the algorithm will print the list of all options set by
  //                     the user including their values and whether they have been used.  In
  //                     some cases this information might be incorrect, due to the internal program flow.
  // - no, don't print options (default value)
  // - yes, print options
  getStringInPList(pvApiCtx, param_in_addr, "print_user_options", &tmp_char, &tmp_res, "no", Log, 
		       CHECK_VALUES, 2, "no", "yes");

  if (tmp_res!=-1) options->SetStringValue("print_user_options", tmp_char);
  FREE(tmp_char);

  return 0;
}
