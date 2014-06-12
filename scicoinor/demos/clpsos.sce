// A LP example which shows all the potentials of CBC Scilab

TestLP      = %T;
TestSOS     = %T;
TestClique  = %T; // YC: pb with this one. See CbcClique constructor. Some parameters are wrong in scicbc
TestLotsize = %T;
TestNWay    = %T;

param = init_param();
param = add_param(param,'maxnumiterations',10000000);
param = add_param(param,'maxnumseconds',10000);
param = add_param(param,'primaltolerance',1e-7);
param = add_param(param,'dualtolerance',1e-7);
param = add_param(param,'dualobjectivelimit',1e100);
param = add_param(param,'primalobjectivelimit',1e100);
param = add_param(param,'verbose',2);
param = add_param(param,'clpverbose',1);
param = add_param(param,'optim_dir', 1); // optimisation direction: 1 - minimize, -1 - maximize, 0 - ignore
param = add_param(param,'stoponfirstsol',0); // we stop once we found a feasible solution

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
//param = add_param(param,'cglgomory_limit',mps_file('nb_int_var'));
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

param = add_param(param,'cbc_dobranchandbound',1); // 0, 1, 2 or 3

c = -[1 1 3 2 2]; // from lpsolve sos documentation

a = [-1 -1  1  1  0; ...
      1  0  1 -3  0];

b_lo = [0   0];
b_up = [30 30];
      
constrtype = 'LL';

lb = [0,  0, 0,   0,   0];
ub = [40, 1, 100, 100, 1];

// The last constraint is a SOS constraint:
// - all the variables involved in the SOS constraint must be integer variables
// - the lower bound for the involved variable must be equal to 0
// - the constraint must be composed of 0 and 1
// - the constraint must be an equality constraint == n, the order of the sos constraint

// SOS constraint order 1
// One variable is used to solve the problem, the others are set to lower bounds
// SOS constraints order 2
// Two adjacent variables are used to solve the problem, the others are set to lower bounds
// SOS constraints order n
// N adjacent variables are used to solve the problem, the others are set to lower bounds

if TestLP then
  printf('test without special constraints\n');

  vartype = 'CCCCC';

  constraints = [];

  [xmin,lambda,status] = clp([],c,a,b_lo,b_up,lb,ub,constrtype,vartype,param);

  printf('solution found: \n');disp(xmin');
  printf('status = %d\n',status);
  printf('value of objective function = %f\n', c*xmin);
  printf('value of constraint no 1 = %f <= %f\n', a(1,:)*xmin, b_up(1));
  printf('value of constraint no 2 = %f <= %f\n', a(2,:)*xmin, b_up(2));
end


if TestSOS then
  printf('\n\n');
  printf('test of SOS order 1\n');

  which  = [0 1 2 3 4];
  weight = [0 1 2 3 4];
  
  vartype = 'CCCCC';

  constraints = [];
  constraints = init_constraint();
  constraints = add_constraint_sos(constraints,1,which,weight,2); // the SOS constraint is to be applied to constraint no 2

  [xmin,lambda,status,extra] = cbc([],c,a,b_lo,b_up,lb,ub,constrtype,vartype,param,constraints);

  printf('solution found: \n');disp(xmin');
  printf('status = %d\n',status);
  printf('value of objective function = %f\n', c*xmin);
  printf('value of constraint no 1 = %f <= %f\n', a(1,:)*xmin, b_up(1));
  printf('value of constraint no 2 = %f <= %f\n', a(2,:)*xmin, b_up(2));
  
  printf('\n\n');
  printf('test of SOS order 2\n');

  constraints = [];
  constraints = init_constraint();
  constraints = add_constraint_sos(constraints,2,which,weight,2); // the SOS constraint is to be applied to constraint no 2

  [xmin,lambda,status,extra] = cbc([],c,a,b_lo,b_up,lb,ub,constrtype,vartype,param,constraints);

  printf('solution found: \n');disp(xmin');
  printf('status = %d\n',status);
  printf('value of objective function = %f\n', c*xmin);
  printf('value of constraint no 1 = %f <= %f\n', a(1,:)*xmin, b_up(1));
  printf('value of constraint no 2 = %f <= %f\n', a(2,:)*xmin, b_up(2));

  // Warning: CBC doesn't accept SOS constraints with an order equal to 0 or an order strictly above 2.
end

if TestClique then
  printf('\n\n');
  printf('test of Clique\n');

  vartype = 'IIIII';

  constraints = [];
  constraints = init_constraint();
  constraints = add_constraint_clique(constraints,1,0,[0 1 2 3 4]);

  [xmin,lambda,status,extra] = cbc([],c,a,b_lo,b_up,lb,ub,constrtype,vartype,param,constraints);

  printf('solution found: \n');disp(xmin');
  printf('status = %d\n',status);
  printf('value of objective function = %f\n', c*xmin);
  printf('value of constraint no 1 = %f <= %f\n', a(1,:)*xmin, b_up(1));
  printf('value of constraint no 2 = %f <= %f\n', a(2,:)*xmin, b_up(2));
end

if TestLotsize then
  printf('\n\n');
  printf('test of Lotsize\n');

  lb = [0,  0,   0,   0, 0];
  ub = [50, 1, 100, 100, 1];
  
  vartype = 'IIIII';

  constraints = [];
  constraints = init_constraint();
  //            add_constraint_lotsize(param,column,_range,_min,_max,lower,upper)
  constraints = add_constraint_lotsize(constraints,1,0,[0 10 20 30 40 50],[],lb,ub);
  // Buggy case: we need to add some tests to check the consistency of the ranges
  //constraints = add_constraint_lotsize(constraints,1,1,[0 10 20 30 40 50],[],lb,ub);

  [xmin,lambda,status,extra] = cbc([],c,a,b_lo,b_up,lb,ub,constrtype,vartype,param,constraints);

  printf('solution found: \n');disp(xmin');
  printf('status = %d\n',status);
  printf('value of objective function = %f\n', c*xmin);
  printf('value of constraint no 1 = %f <= %f\n', a(1,:)*xmin, b_up(1));
  printf('value of constraint no 2 = %f <= %f\n', a(2,:)*xmin, b_up(2));

  lb = [0,  0, 0,   0,   0];
  ub = [40, 1, 100, 100, 1];
end

if TestNWay then
  printf('\n\n');
  printf('test of NWay\n');

  vartype = 'IIIII';
  
  constraints = [];
  constraints = init_constraint();
  constraints = add_constraint_nway(constraints,[0 1 2 3 4]);

  // At the solution point, on variable must be at upper bound and the others must be at lower bounds
  
  [xmin,lambda,status,extra] = cbc([],c,a,b_lo,b_up,lb,ub,constrtype,vartype,param,constraints);

  printf('solution found: \n');disp(xmin');
  printf('status = %d\n',status);
  printf('value of objective function = %f\n', c*xmin);
  printf('value of constraint no 1 = %f <= %f\n', a(1,:)*xmin, b_up(1));
  printf('value of constraint no 2 = %f <= %f\n', a(2,:)*xmin, b_up(2));
end

