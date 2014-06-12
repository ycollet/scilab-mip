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

path_new = get_absolute_file_path('mps_cbc_test.sce');
path_old = pwd();

cd(path_new)

exec('coinor_data.sce');

/////////////////////////////
// Choose problem to solve //
/////////////////////////////

// See listoffiles.html for informations related to these problems

Solve_miplib  = %F; // 92 files - pb 6,7,10 too large to fit in stack
Solve_lplib   = %F; // 127 files
Solve_mitt_c  = %F; // 10 files
Solve_bi_sc   = %F; // 21 files
Solve_Coral   = %F; // 372 files
Solve_miplib3 = %F; // 164 files
Solve_sample  = %F; // 18 files
Solve_infeas  = %F; // 29 files
Solve_big     = %F; // 1 file
Solve_netlib  = %F; // 90 files
Solve_data    = %T; // 2 files

///////////////////////////
// Set global parameters //
///////////////////////////

Log              = 1; // 6 - maximum informations // 3 - sum-up after each branching // 1 - minimum sum-up, fast with this parameter
clpLog           = 0; // 6 - maximum informations // 3 - sum-up after each branching // 1 - minimum sum-up, fast with this parameter
mps_index        = 2; // 61, 65, 33 et 8 // Pb with 33 (too many nodes ?) and 8 and miplib - miplib 70 for fast tests
mps_type         = 0;
UseLinpro        = %F;
PrintSol         = %F;
DoBranchAndBound = %T;
ManualInit       = %F;

if Solve_miplib  then mps_filename = mps_miplib(mps_index);         end
if Solve_lplib   then mps_filename = mps_lplib(mps_index);          end
if Solve_mitt_c  then mps_filename = mps_mitt_c(mps_index);         end
if Solve_bi_sc   then mps_filename = mps_bi_sc(mps_index);          end
if Solve_Coral   then mps_filename = mps_Coral(mps_index);          end
if Solve_miplib3 then mps_filename = mps_coinor_miplib3(mps_index); end
if Solve_sample  then mps_filename = mps_coinor_sample(mps_index);  end
if Solve_infeas  then mps_filename = mps_coinor_infeas(mps_index);  end
if Solve_big     then mps_filename = mps_coinor_big(mps_index);     end
if Solve_netlib  then mps_filename = mps_coinor_netlib(mps_index);  end
if Solve_data    then mps_filename = mps_coinor_data(mps_index);    end

/////////////////////////////
// Gunzip the gzipped file //
/////////////////////////////

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

////////////////////////////////////////////
// Gzip the file before processing to CBC //
////////////////////////////////////////////

if is_gzipped then
  if getos() == 'Windows' then
    unix_w('.' + filesep() + 'data' + filesep() + 'gzip' + filesep() + 'gzip.exe ' + mps_filename);
  else
    unix('gzip ' + mps_filename);
  end
end

t_start = t_end;

printf('elapsed time for reading file = %f secondes\n', etime(t_end, t_start));

param = init_param();
param = add_param(param,'maxnumiterations',10000000);
param = add_param(param,'maxnumseconds',10000);
param = add_param(param,'primaltolerance',1e-7);
param = add_param(param,'dualtolerance',1e-7);
param = add_param(param,'dualobjectivelimit',1e100);
param = add_param(param,'primalobjectivelimit',1e100);
param = add_param(param,'verbose',Log);
param = add_param(param,'clpverbose',clpLog);
param = add_param(param,'cbcmaininit',1);
param = add_param(param,'writemps','test.mps');
param = add_param(param,'optim_dir', 1); // optimisation direction: 1 - minimize, -1 - maximize, 0 - ignore
param = add_param(param,'stoponfirstsol',0); // we stop once we found a feasible solution

if ManualInit then
  ///////////////////////////
  // Set the compare class //
  ///////////////////////////

  ///////////////////////
  // cbccomparedefault //
  ///////////////////////
  // Description of the method:
  // Description of the options:
  //      - DOUBLE: cbccomparedefault_weight       - 0
  //      - INT:    cbccomparedefault_breadthdepth - 0
  // Uncomment these add_param functions to test:
  param = add_param(param,'cbccomparedefault',1); // YC: 6.46 on miplib 70
  param = add_param(param,'cbccomparedefault_weight',1);
  param = add_param(param,'cbccomparedefault_breadthdepth',10);

  /////////////////////
  // cbccomparedepth //
  /////////////////////
  // Description of the method:
  // Description of the options:
  // Uncomment these add_param functions to test:
  //param = add_param(param,'cbccomparedepth',1); // YC: 7.41 on miplib 70

  ////////////////////////
  // cbccompareestimate //
  ////////////////////////
  // Description of the method:
  // Description of the options:
  // Uncomment these add_param functions to test:
  //param = add_param(param,'cbccompareestimate',1); // YC: 5.06 on miplib 70

  /////////////////////////
  // cbccompareobjective //
  /////////////////////////
  // Description of the method:
  // Description of the options:
  // Uncomment these add_param functions to test:
  //param = add_param(param,'cbccompareobjective',0); // 9.26 on miplib 70

  ////////////////////
  // cbccompareuser //
  ////////////////////
  // Description of the method:
  // Description of the options:
  //      - DOUBLE: ccompareuser_weight          - 0.0
  // Uncomment these add_param functions to test:
  //param = add_param(param,'cbccompareuser',0.0); // 9.26 on miplib 70
  //param = add_param(param,'cbccompareuser_weight',1.0);

  ///////////////////////////////
  // Set the Cgl cut generator //
  ///////////////////////////////

  ///////////////////
  // cglpreprocess //
  ///////////////////
  // Description of the method:
  // Description of the options:
  //param = add_param(param,'cglpreprocess',1); // YC: peut etre pb avec cette option
  // 5 is the number of passes
  param = add_param(param,'cglpreprocess_other',5); // YC: ultra fast on miplib 70
  //param = add_param(param,'cglpreprocess_nondefault',5);
  //// 1 -> clique, 2 -> SOS
  //param = add_param(param,'cglpreprocess_nondefault_type',2);

  ////////////////
  // cglprobing //
  ////////////////
  // Description of the method:
  // Description of the options:
  //      - INT: cglprobing_mode            - 0
  //      - INT: cglprobing_maxpass         - 0
  //      - INT: cglprobing_loglevel        - 0
  //      - INT: cglprobing_maxprobe        - 0
  //      - INT: cglprobing_maxlook         - 0
  //      - INT: cglprobing_maxelements     - 0
  //      - INT: cglprobing_maxpassroot     - 0
  //      - INT: cglprobing_maxproberoot    - 0
  //      - INT: cglprobing_maxlookroot     - 0
  //      - INT: cglprobing_maxelementsroot - 0
  //      - INT: cglprobing_rowcuts         - 0
  //      - INT: cglprobing_usingobjective  - 0
  // Uncomment these add_param functions to test:
  param = add_param(param,'cglprobing',-1);
  param = add_param(param,'cglprobing_mode',0); // 0, 1 or 2
  param = add_param(param,'cglprobing_loglevel',0); // 0, 1 or 2
  //param = add_param(param,'cglprobing_maxelement',200);
  //param = add_param(param,'cglprobing_maxpassroot',5);
  //param = add_param(param,'cglprobing_maxproberoot',1000);
  //param = add_param(param,'cglprobing_maxlookroot',50);
  //param = add_param(param,'cglprobing_maxelementroot',500);
  param = add_param(param,'cglprobing_usingobjective',1); // from sample1.cpp
  param = add_param(param,'cglprobing_maxpass',3); // from sample1.cpp
  param = add_param(param,'cglprobing_maxprobe',100); // from sample1.cpp
  param = add_param(param,'cglprobing_maxlook',50); // from sample1.cpp
  param = add_param(param,'cglprobing_rowcuts',3); // from sample1.cpp

  ///////////////
  // cglgomory //
  ///////////////
  // Description of the method:
  // Description of the options:
  //      - INT: cglgomory_limit       - 50
  //      - INT: cglgomory_limitatroot - 50
  //      - DOUBLE: cglgomory_away     - 0.05
  // Uncomment these add_param functions to test:
  param = add_param(param,'cglgomory',-1);
  //By default, does not generate cuts with more than 50 non zero coefficients. 
  // To get more dense cuts, modify the parameter Limit. See the parameter cglgomory_limit. 
  param = add_param(param,'cglgomory_limit',mps_file('nb_int_var'));
  //param = add_param(param,'cglgomory_limitatroot',50);
  //param = add_param(param,'cglgomory_away',0.05);

  //////////////////////
  // cglknapsackcover //
  //////////////////////
  // Description of the method:
  // Description of the options:
  //      - INT: cglknapsackcover_maxinknapsack      - 0
  //      - INT: cglknapsackcover_switchoffexpensive - 0
  //      - INT: cglknapsackcover_switchonexpensive  - 0
  // Uncomment these add_param functions to test:
  //param = add_param(param,'cglknapsackcover',-1);
  //param = add_param(param,'cglknapsackcover_maxinknapsack',0);
  //param = add_param(param,'cglknapsackcover_switchoffexpensive',0);
  //param = add_param(param,'cglknapsackcover_switchonexpensive',0);

  /////////////////
  // cglredsplit //
  /////////////////
  // Description of the method:
  // Description of the options:
  //      - INT:    cglredsplit_limit         - 50
  //      - DOUBLE: cglredsplit_away          - 0.05
  //      - DOUBLE: cglredsplit_lub           - 1000
  //      - DOUBLE: cglredsplit_eps           - 1e-7
  //      - DOUBLE: cglredsplit_eps_coeff     - 1e-8
  //      - DOUBLE: cglredsplit_eps_coeff_lub - 1e-13
  //      - DOUBLE: cglredsplit_eps_relax     - 1e-8
  //      - DOUBLE: cglredsplit_normiszero    - 1e-5
  //      - DOUBLE: cglredsplit_minreduc      - 0.05
  //      - DOUBLE: cglredsplit_maxtab        - 100
  // Uncomment these add_param functions to test:
  param = add_param(param,'cglredsplit',-1);
  //param = add_param(param,'cglredsplit_limit',50);
  //param = add_param(param,'cglredsplit_away',0.05);
  //param = add_param(param,'cglredsplit_lub',1000);
  //param = add_param(param,'cglredsplit_eps',1e-7);
  //param = add_param(param,'cglredsplit_eps_coeff',1e-8);
  //param = add_param(param,'cglredsplit_eps_coeff_lub',1e-13);
  //param = add_param(param,'cglredsplit_eps_relax',1e-8);
  //param = add_param(param,'cglredsplit_normiszero',1e-5);
  //param = add_param(param,'cglredsplit_minreduc',0.05);
  //param = add_param(param,'cglredsplit_maxtab',100);

  ///////////////
  // cglclique //
  ///////////////
  // Description of the method:
  // Description of the options:
  //      - INT:    cglclique_starcliquecandidatelengththreshold - 0
  //      - INT:    cglclique_rowcliquecandidatelengththreshold  - 0
  //      - INT:    cglclique_StarCliqueReport - 0
  //      - INT:    cglclique_rowcliquereport  - 0
  //      - INT:    cglclique_dostarclique     - 0
  //      - INT:    cglclique_dorowclique      - 0
  //      - DOUBLE: cglclique_minviolation     - 0
  // Uncomment these add_param functions to test:
  //param = add_param(param,'cglclique',-1);
  //param = add_param(param,'cglclique_starcliquecandidatelengththreshold',0);
  //param = add_param(param,'cglclique_rowcliquecandidatelengththreshold',0);
  //param = add_param(param,'cglclique_starcliquereport',0);
  //param = add_param(param,'cglclique_rowcliquereport',0);
  //param = add_param(param,'cglclique_dostartclique',0);
  //param = add_param(param,'cglclique_dorowclisuqe',0);
  //param = add_param(param,'cglclique_minviolation',0);

  /////////////////////////////
  // cglmixedintegerrounding //
  /////////////////////////////
  // Description of the method:
  // Description of the options:
  //      - INT: cglmixedintegerrounding_maxaggr   - 0
  //      - INT: cglmixedintegerrounding_multiply  - 0
  //      - INT: cglmixedintegerrounding_criterion - 0
  //      - INT: cglmixedintegerrounding_dopreproc - 0
  // YC: pb avec cette option
  // Uncomment these add_param functions to test:
  param = add_param(param,'cglmixedintegerrounding',-1);
  //param = add_param(param,'cglmixedintegerrounding_maxaggr',0);
  //param = add_param(param,'cglmixedintegerrounding_multiply',0);
  //param = add_param(param,'cglmixedintegerrounding_criterion',0);
  //param = add_param(param,'cglmixedintegerrounding_dopreproc',0);

  //////////////////////////////
  // cglmixedintegerrounding2 //
  //////////////////////////////
  // Description of the method:
  // Description of the options:
  //      - INT: cglmixedintegerrounding2_maxaggr   - 0
  //      - INT: cglmixedintegerrounding2_multiply  - 0
  //      - INT: cglmixedintegerrounding2_criterion - 0
  //      - INT: cglmixedintegerrounding2_dopreproc - 0
  // Uncomment these add_param functions to test:
  // YC: pb avec cette option
  //param = add_param(param,'cglmixedintegerrounding2',-1);
  //param = add_param(param,'cglmixedintegerrounding2_maxaggr',0);
  //param = add_param(param,'cglmixedintegerrounding2_multiply',0);
  //param = add_param(param,'cglmixedintegerrounding2_criterion',0);
  //param = add_param(param,'cglmixedintegerrounding2_dopreproc',0);

  //////////////////
  // cglflowcover //
  //////////////////
  // Description of the method:
  // Description of the options:
  //      - INT: cglflowcover_numflowcuts - 0
  // Uncomment these add_param functions to test:
  param = add_param(param,'cglflowcover',-1);
  //param = add_param(param,'cglflowcover_numflowcuts',0);

  ////////////////
  // cgloddhole //
  ////////////////
  // Description of the method:
  // Description of the options:
  //      - DOUBLE: cgloddhole_minimumviolation    - 0.005
  //      - DOUBLE: cgloddhole_minimumviolationper - 0.00002
  //      - INT:    cgloddhole_maximumentries      - 200
  // Uncomment these add_param functions to test:
  //param = add_param(param,'cgloddhole',-1);
  //param = add_param(param,'cgloddhole_minimumviolation',0.005);
  //param = add_param(param,'cgloddhole_minimumviolationper',0.00002);
  //param = add_param(param,'cgloddhole_maximumentries',200);

  ///////////////
  // cgltwomir //
  ///////////////
  // Description of the method:
  // Description of the options:
  //      - INT: cgltwomir_mirscal_tmin    - tmp_int_1, 1);
  //      - INT: cgltwomir_mirscal_tmax    - tmp_int_2, tmp_int_1);
  //      - INT: cgltwomir_mirscal_qmin    - tmp_int_1, 1);
  //      - INT: cgltwomir_mirscal_qmax    - tmp_int_2, tmp_int_1);
  //      - INT: cgltwomir_amax            - 1
  //      - INT: cgltwomir_maxelements     - 200
  //      - INT: cgltwomir_cuttype_mir     - 0
  //      - INT: cgltwomir_cuttype_twomir  - 0
  //      - INT: cgltwomir_cuttype_tab     - 0
  //      - INT: cgltwomir_cuttype_form    - 0
  //      - INT: cgltwomir_formulationrows - 10
  // Uncomment these add_param functions to test:
  param = add_param(param,'cgltwomir',-1);
  //param = add_param(param,'cgltwomir_mirscal_tmin',1);
  //param = add_param(param,'cgltwomir_mirscal_tmax',1.1);
  //param = add_param(param,'cgltwomir_mirscal_qmin',1);
  //param = add_param(param,'cgltwomir_mirscal_qmax',1.1);
  //param = add_param(param,'cgltwomir_amax',1);
  //param = add_param(param,'cgltwomir_maxelements',200);
  param = add_param(param,'cgltwomir_cuttype_mir',0);
  param = add_param(param,'cgltwomir_cuttype_twomir',1);
  param = add_param(param,'cgltwomir_cuttype_tab',0);
  param = add_param(param,'cgltwomir_cuttype_form',0);
  param = add_param(param,'cgltwomir_formulation_rows',10);

  /////////////////////
  // cglalldifferent //
  /////////////////////
  // Description of the method:
  // Description of the options:
  //      - INT: cglalldifferent_maxlook - 0
  // Uncomment these add_param functions to test:
  //param = add_param(param,'cglalldifferent',-1);
  //param = add_param(param,'cglalldifferent_maxlook',0);

  /////////////////////
  // cglduplicaterow //
  /////////////////////
  // Description of the method:
  // Description of the options:
  //      - INT: cglduplicaterow_maximumrhs       - 0
  //      - INT: cglduplicaterow_maximumdominated - 0
  //      - INT: cglduplicaterow_mode             - 0
  // Uncomment these add_param functions to test:
  //param = add_param(param,'cglduplicaterow',-1);
  //param = add_param(param,'cglduplicaterow_maximumrhs',0);
  //param = add_param(param,'cglduplicaterow_maximumdominated',0);
  //param = add_param(param,'cglduplicaterow_mode',0);

  ///////////////////////
  // cglliftandproject //
  ///////////////////////
  // Description of the method:
  // Description of the options:
  //      - INT: cglliftandproject_beta - 1
  // Uncomment these add_param functions to test:
  //param = add_param(param,'cglliftandproject',-1);
  //param = add_param(param,'cglliftandproject_beta',1);

  /////////////////////////
  // cglresidualcapacity //
  /////////////////////////
  // Description of the method:
  // Description of the options:
  //      - DOUBLE: cglresidualcapacity_epsilon   - 0
  //      - DOUBLE: cglresidualcapacity_tolerance - 0
  //      - INT:    cglresidualcapacity_dopreproc - 0
  // Uncomment these add_param functions to test:
  //param = add_param(param,'cglresidualcapacity',-1);
  //param = add_param(param,'cglresidualcapacity_epsilon',1e-6);
  //param = add_param(param,'cglresidualcapacity_tolerance',1e-6);
  //param = add_param(param,'cglresidualcapacity_dopreproc',1);

  ////////////////////
  // cglimplication //
  ////////////////////
  // Description of the method:
  // Description of the options:
  // Uncomment these add_param functions to test:
  //param = add_param(param,'cglimplication',-1);

  ///////////////////////
  // cglsimplerounding //
  ///////////////////////
  // Description of the method:
  // Description of the options:
  // Uncomment these add_param functions to test:
  param = add_param(param,'cglsimplerounding',-1);

  //////////////
  // cgllandp //
  //////////////
  // Description of the method:
  // Description of the options:
  // Uncomment these add_param functions to test:
  //param = add_param(param,'cgllandp',-1); // YC: gives a "proven infeasible" on miplib 70

  /////////////////////////
  // Set some heuristics //
  /////////////////////////

  /////////////////
  // cbcrounding //
  /////////////////
  // Description of the method:
  // Description of the options:
  //      - INT: cbcrounding_seed - 12345
  // Uncomment these add_param functions to test:
  param = add_param(param,'cbcrounding',1);
  //param = add_param(param,'cbcrounding_seed',12345);

  /////////////////////////////////
  // cbcheuristicdivecoefficient //
  /////////////////////////////////
  // Description of the method:
  // Description of the options:
  // Uncomment these add_param functions to test:
  //param = add_param(param,'cbcheuristicdivecoefficient',0);

  ////////////////////////////////
  // cbcheuristicdivefractional //
  ////////////////////////////////
  // Description of the method:
  // Description of the options:
  // Uncomment these add_param functions to test:
  //param = add_param(param,'cbcheuristicdivefractional',0);

  ////////////////////////////
  // cbcheuristicdiveguided //
  ////////////////////////////
  // Description of the method:
  // Description of the options:
  // Uncomment these add_param functions to test:
  //param = add_param(param,'cbcheuristicdiveguided',0); // YC: fast on miplib 70

  //////////////////////////////////
  // cbcheuristicdivevectorlength //
  //////////////////////////////////
  // Description of the method:
  // Description of the options:
  // Uncomment these add_param functions to test:
  //param = add_param(param,'cbcheuristicdivevectorlength',0);

  ///////////////////////
  // cbcheuristicfpump //
  ///////////////////////
  // Description of the method:
  // Description of the options:
  //      - DOUBLE: cbcheuristicfpump_maximumtime       - 0
  //      - DOUBLE: cbcheuristicfpump_fakecutoff        - COIN_DBL_MAX
  //      - DOUBLE: cbcheuristicfpump_absoluteincrement - 0.0
  //      - DOUBLE: cbcheuristicfpump_relativeincrement - 0.0
  //      - DOUBLE: cbcheuristicfpump_defaultrounding   - 0.5
  //      - DOUBLE: cbcheuristicfpump_initialweight     - 0.0
  //      - DOUBLE: cbcheuristicfpump_weightfactor      - 0.1
  //      - DOUBLE: cbcheuristicfpump_artificialcost    - COIN_DBL_MAX
  //      - INT:    cbcheuristicfpump_maximumpasses     - 100
  //      - INT:    cbcheuristicfpump_maximumretries    - 1
  //      - INT:    cbcheuristicfpump_accumulate        - 0
  //      - INT:    cbcheuristicfpump_fixonreducedcosts - 0
  // Uncomment these add_param functions to test:
  //param = add_param(param,'cbcheuristicfpump',0); // YC: good choice for miplib 70
  //param = add_param(param,'cbcheuristicfpump_maximumtime',0);
  //param = add_param(param,'cbcheuristicfpump_fakecutoff',1e32);
  //param = add_param(param,'cbcheuristicfpump_absoluteincrement',0.0);
  //param = add_param(param,'cbcheuristicfpump_relativeincrement',0.0);
  //param = add_param(param,'cbcheuristicfpump_defaultrounding',0.5);
  //param = add_param(param,'cbcheuristicfpump_initialweight',0);
  //param = add_param(param,'cbcheuristicfpump_weightfactor',0.1);
  //param = add_param(param,'cbcheuristicfpump_artificialcost',1e32);
  //param = add_param(param,'cbcheuristicfpump_maximumpasses',100);
  //param = add_param(param,'cbcheuristicfpump_maximumretries',1);
  //param = add_param(param,'cbcheuristicfpump_accumulate',0);
  //param = add_param(param,'cbcheuristicfpump_fixonreducedcosts',0);

  /////////////////////////////
  // cbcheuristicgreedycover //
  /////////////////////////////
  // Description of the method:
  // Description of the options:
  //      - INT: cbcheuristicgreedycover_algorithm   - 0
  //      - INT: cbcheuristicgreedycover_numbertimes - 0
  // Uncomment these add_param functions to test:
  //param = add_param(param,'cbcheuristicgreedycover',0);
  //param = add_param(param,'cbcheuristicgreedycover_algorithm',0);
  //param = add_param(param,'cbcheuristicgreedycover_numbertimes',0);

  ////////////////////////////////
  // cbcheuristicgreedyequality //
  ////////////////////////////////
  // Description of the method:
  // Description of the options:
  //      - INT:    cbcheuristicgreedyequality_algorithm   - 0
  //      - DOUBLE: cbcheuristicgreedyequality_fraction    - 0
  //      - INT:    cbcheuristicgreedyequality_numbertimes - 0
  // Uncomment these add_param functions to test:
  //param = add_param(param,'cbcheuristicgreedyequality',0); // YC: good choice for miplib 70
  //param = add_param(param,'cbcheuristicgreedyequality_algorithm',0);
  //param = add_param(param,'cbcheuristicgreedyequality_fraction',0);
  //param = add_param(param,'cbcheuristicgreedyequality_numbertimes',0);

  ///////////////////////
  // cbcheuristiclocal //
  ///////////////////////
  // Description of the method:
  // Description of the options:
  //      - INT: cbcheuristiclocal_searchtype - 0
  // Uncomment these add_param functions to test:
  //param = add_param(param,'cbcheuristiclocal',1);
  //param = add_param(param,'cbcheuristiclocal_searchtype',1);

  /////////////////////////
  // cbcheuristicpartial //
  /////////////////////////
  // Description of the method:
  // Description of the options:
  //      - INT: cbcheuristicpartial_fixpriority - 0
  // Uncomment these add_param functions to test:
  //param = add_param(param,'cbcheuristicpartial',0);
  //param = add_param(param,'cbcheuristicpartial_fixpriority',0);

  //////////////////////
  // cbcheuristicrens //
  //////////////////////
  // Description of the method:
  // Description of the options:
  // Uncomment these add_param functions to test:
  //param = add_param(param,'cbcheuristicrens',0);

  //////////////////////
  // cbcheuristicrins //
  //////////////////////
  // Description of the method:
  // Description of the options:
  //      - INT:    cbcheuristicrins_howoften    - 0
  //      - DOUBLE: cbcheuristicrins_decayfactor - 0
  // Uncomment these add_param functions to test:
  //param = add_param(param,'cbcheuristicrins',0); // YC: fast on miplib 70
  //param = add_param(param,'cbcheuristicrins_howoften',0);
  //param = add_param(param,'cbcheuristicrins_decayfactor',0);

  ////////////////////
  // cbcserendipity //
  ////////////////////
  // Description of the method:
  // Description of the options:
  // Uncomment these add_param functions to test:
  //param = add_param(param,'cbcserendipity',0);

  //////////////////////////////////
  // Set Branching decision class //
  //////////////////////////////////

  //////////////////////////////
  // cbcbranchdefaultdecision //
  //////////////////////////////
  // Description of the method:
  // Description of the options:
  //      - DOUBLE: cbcbranchdefaultdecision_bestcriterion - 0
  // Uncomment these add_param functions to test:
  param = add_param(param,'cbcbranchdefaultdecision',1);
  //param = add_param(param,'cbcbranchdefaultdecision_bestcriterion',1);

  //////////////////////////////
  // cbcbranchdynamicdecision //
  //////////////////////////////
  // Description of the method:
  // Description of the options:
  //      - DOUBLE: cbcbranchdynamicdecision_bestcriterion - 0
  // Uncomment these add_param functions to test:
  //param = add_param(param,'cbcbranchdynamicdecision',1);
  //param = add_param(param,'cbcbranchdynamicdecision_bestcriterion',0);

  ////////////////////////
  // Set Strategy class //
  ////////////////////////

  ///////////////////////////////
  // cbcstrategydefaultsubtree //
  ///////////////////////////////
  // Description of the method:
  // Description of the options:
  //      - INT: cbcstrategydefaultsubtree_cutsonlyatroot    - false
  //      - INT: cbcstrategydefaultsubtree_numberstrong      - 5
  //      - INT: cbcstrategydefaultsubtree_numberbeforetrust - 0
  //      - INT: cbcstrategydefaultsubtree_printlevel        - 0
  // Uncomment these add_param functions to test:1
  //param = add_param(param,'cbcstrategydefaultsubtree',1);
  //YC: if just the option above, an exception is raised. See the default parameters of this option.
  //YC: if just the option below, an exception is raised. See the default parameters of this option.
  //param = add_param(param,'cbcstrategydefaultsubtree_cutsonlyatroot',1);
  //param = add_param(param,'cbcstrategydefaultsubtree_numberstrong',5);
  //param = add_param(param,'cbcstrategydefaultsubtree_numberbeforetrust',10);
  //param = add_param(param,'cbcstrategydefaultsubtree_printlevel',1);

  ////////////////////////
  // cbcstrategydefault //
  ////////////////////////
  // Description of the method:
  // Description of the options:
  //      - INT: cbcstrategydefault_cutsonlyatroot    - false
  //      - INT: cbcstrategydefault_numberstrong      - 5
  //      - INT: cbcstrategydefault_numberbeforetrust - 0
  //      - INT: cbcstrategydefault_printlevel        - 0
  // Uncomment these add_param functions to test:
  //param = add_param(param,'cbcstrategydefault',1);
  //YC: if just the option above, an exception is raised. See the default parameters of this option.
  //param = add_param(param,'cbcstrategydefault_cutsonlyatroot',1);
  //param = add_param(param,'cbcstrategydefault_numberstrong',5);
  //param = add_param(param,'cbcstrategydefault_numberbeforetrust',10);
  //param = add_param(param,'cbcstrategydefault_printlevel',1);

  /////////////////////
  // cbcstrategynull //
  /////////////////////
  // Description of the method:
  // Description of the options:
  // Uncomment these add_param functions to test:
  param = add_param(param,'cbcstrategynull',1);

  //////////////////////////////////////
  // Some Cbc parameters are settable //
  //////////////////////////////////////

  // Description of the options:
  //  DOUBLE: cbc_cutoff                  - 0
  //  INT:    cbc_maximumnodes            - 200
  //  INT:    cbc_maximumsolutions        - 0
  //  INT:    cbc_printingmode            - 0
  //  DOUBLE: cbc_maximumseconds          - 1000
  //  DOUBLE: cbc_integertolerance        - 0
  //  DOUBLE: cbc_infeasibilityweight     - 0
  //  DOUBLE: cbc_allowablegap            - 0
  //  DOUBLE: cbc_allowablefractiongap    - 0
  //  DOUBLE: cbc_allowablepercentagegap  - 0
  //  DOUBLE: cbc_cutoffincrement         - 0
  //  DOUBLE: cbc_minimumdrop             - 0
  //  INT:    cbc_maxcutpassesatroot      - 20
  //  INT:    cbc_maxcutpasses            - 10
  //  INT:    cbc_numberstrong            - 5
  //  INT:    cbc_preferredway            - 0
  //  INT:    cbc_numberbeforetrust       - 0
  //  INT:    cbc_numberpenalties         - 0
  //  INT:    cbc_numberanalyzeiterations - 0
  //  DOUBLE: cbc_penaltyscalefactor      - 0
  //  INT:    cbc_problemtype             - 0
  //  INT:    cbc_howoftenglobalscan      - 0
  //  INT:    cbc_printfrequency          - 0

  // Uncomment these add_param functions to test:
  // YC: trouver de meilleures valeurs par defaut.
  param = add_param(param,'cbc_cutoff', 1e100);               // YC: basee sur la fonction objectif du pb relaxe ?
                                                              // YC: 1e100 (CurrentCutoff constructeur de CbcModel)
  //param = add_param(param,'cbc_cutoffincrement', 1e-5);       // YC: 1e-5
  param = add_param(param,'cbc_printingmode', 1); // 0 or 1
  //param = add_param(param,'cbc_maximumsolutions', 9999999);   // YC: 9999999
  //param = add_param(param,'cbc_integertolerance', 1e-6);      // YC: 1e-6
  //param = add_param(param,'cbc_infeasibilityweight', 0);    // YC: 0.0
  //param = add_param(param,'cbc_allowablegap', 1e-10);       // YC: 1e-10
  //param = add_param(param,'cbc_allowablefractiongap', 0);   // YC: 0.0
  //param = add_param(param,'cbc_allowablepercentagegap', 0);
  //param = add_param(param,'cbc_minimumdrop', 0);            // YC: min(1,fabs(fobj)*1e-3+1e-4)
  param = add_param(param,'cbc_maxcutpassesatroot', 30); // YC: 100 30
  // plus c'est grand, plus le nombre de noeuds explose lorsque l'on teste un noeud.
  // L'occupation memoire augmente rapidement et on peut obtenir facilement un plantage sur de petites instances
  param = add_param(param,'cbc_maxcutpasses', 5);         // YC: 5
  //param = add_param(param,'cbc_preferredway', 0);           // YC: 0
  param = add_param(param,'cbc_numberbeforetrust', 10);     // YC: 10
  //param = add_param(param,'cbc_numberpenalties', 20);       // YC: 20
  //param = add_param(param,'cbc_numberanalyzeiterations', 0);  // YC: 0
  //param = add_param(param,'cbc_penaltyscalefactor', 3.0);   // YC: 3.0
  //param = add_param(param,'cbc_problemtype', 0);            // YC: 0
  param = add_param(param,'cbc_howoftenglobalscan', 10);       // YC: 1 - this option produce a lot of output ?
  param = add_param(param,'cbc_numberstrong', 5);             // YC: entre 5
  param = add_param(param,'cbc_printfrequency', 10);         // YC: 0
  param = add_param(param,'cbc_maximumnodes', 100000);    // YC: 2147483647
  param = add_param(param,'cbc_maximumseconds', 10000);       // YC: 1e100
  //param = add_param(param,'cbc_numberthreads',2);
end

if DoBranchAndBound then
  param = add_param(param,'cbc_dobranchandbound',1); // 0, 1, 2 or 3
end

///////////////////////////
// Set the variable type //
///////////////////////////

// 'I' -> integer
// 'C' -> continuous

var_type = string(zeros(1,length(mps_file_mp('obj_var_is_int'))));
var_type(find(mps_file_mp('obj_var_is_int')==0)) = 'C';
var_type(find(mps_file_mp('obj_var_is_int')==1)) = 'I';
var_type = strcat(var_type);

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

[xmin,fmin,status,extra] = cbc([],mps_file_mp('obj_coeff'),mps_file_mp('constr_mat'),mps_file_mp('lhs'),mps_file_mp('rhs'), ...
                               mps_file_mp('bounds_lower'),mps_file_mp('bounds_upper'),btype,var_type,param);

printf(' status:\n');
printf(' - -1 before branchAndBound\n');
printf(' - 0 finished - check isProvenOptimal or isProvenInfeasible to see if solution found (or check value of best solution)\n');
printf(' - 1 stopped - on maxnodes, maxsols, maxtime \n');
printf(' - 2 difficulties so run was abandoned \n');
printf(' - 5 event user programmed event occurred\n');
printf('\n');
printf(' fields of the extra mlist:\n');
printf(' extra(''lambda'')\n');
printf(' extra(''redcosts'')\n');
printf(' extra(''time'')\n');
printf(' extra(''cbc_status'')\n');
printf('   - bit 1: isAbandoned ?            Are there a numerical difficulties ?\n');
printf('   - bit 2: isProvenOptimal ?        Is optimality proven ?\n');
printf('   - bit 3: isProvenInfeasible ?     Is infeasiblity proven (or none better than cutoff) ?\n');
printf('   - bit 4: isContinuousUnbounded ?  Was continuous solution unbounded ?\n');
printf('   - bit 5: isProvenDualInfeasible ? Was continuous solution unbounded ?\n');
printf('   - bit 6: isNodeLimitReached ?     Node limit reached ?\n');
printf('   - bit 7: isSecondsLimitReached ?  Time limit reached ?\n');
printf('   - bit 8: isSolutionLimitReached ? Solution limit reached ?\n');
printf(' extra(''iteration_count'')\n');
printf('   - Get how many iterations it took to solve the problem.\n');
printf(' extra(''node_count'')\n');
printf('   - Get how many Nodes it took to solve the problem.\n');
printf(' extra(''secondary_status'')\n');
printf('   - -1 unset (status_ will also be -1) \n');
printf('   -  0 search completed with solution \n');
printf('   -  1 linear relaxation not feasible (or worse than cutoff) \n');
printf('   -  2 stopped on gap \n');
printf('   -  3 stopped on nodes \n');
printf('   -  4 stopped on time \n');
printf('   -  5 stopped on user event\n'); 
printf('   -  6 stopped on solutions \n');
printf('   -  7 linear relaxation unbounded.\n');
printf(' extra(''initsolve_status'')\n');
printf('   - bit 1: isInitialSolveAbandoned ?              Are there numerical difficulties (for initialSolve) ?\n');
printf('   - bit 2: isInitialSolveProvenOptimal ?          Is optimality proven (for initialSolve) ?\n');
printf('   - bit 3: isInitialSolveProvenPrimalInfeasible ? Is primal infeasiblity proven (for initialSolve) ?\n');
printf('   - bit 4: isInitialSolveProvenDualInfeasible ?   Is dual infeasiblity proven (for initialSolve) ? \n');
printf('\n');
t_end = getdate();

printf('nb of constr         = %d\n', size(mps_file_mp('constr_mat'),1));
printf('nb of var            = %d\n', size(mps_file_mp('constr_mat'),2));
printf('nb of int var        = %d\n', mps_file('nb_int_var'));
printf('status               = %d\n', status);
printf('cbc status           = %s\n', dec2bin(extra('cbc_status')));
printf('cbc secondary status = %d\n', extra('secondary_status'));
printf('cbc init status      = %s\n', dec2bin(extra('initsolve_status')));
printf('cbc node count       = %d\n', extra('node_count'));
printf('cbc iterations count = %d\n', extra('iteration_count'));
printf('cbc time             = %f\n', extra('time'));

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
