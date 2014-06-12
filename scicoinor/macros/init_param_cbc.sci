function param = init_param_cbc(x0, param)
if ~isdef('param','local') then
  param = init_param();
end
if typeof(param)~='plist' then
  param = init_param();
end

param = add_param(param,'maxnumiterations',10000000);
param = add_param(param,'maxnumseconds',1000);
param = add_param(param,'primaltolerance',1e-7);
param = add_param(param,'dualtolerance',1e-7);
//param = add_param(param,'cbcmaininit',1);

////////////////////////
// cbccompareestimate //
////////////////////////
// Description of the method:
// Description of the options:
// Uncomment these add_param functions to test:
param = add_param(param,'cbccompareestimate',1); // YC: 5.06 on miplib 70

///////////////////////////////
// Set the Cgl cut generator //
///////////////////////////////

///////////////////
// cglpreprocess //
///////////////////
// Description of the method:
// Description of the options:
param = add_param(param,'cglpreprocess_other',5); // YC: ultra fast on miplib 70

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
param = add_param(param,'cglgomory_limit',length(x0));

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
param = add_param(param,'cgltwomir_cuttype_twomir',1);

///////////////////////
// cglsimplerounding //
///////////////////////
// Description of the method:
// Description of the options:
// Uncomment these add_param functions to test:
param = add_param(param,'cglsimplerounding',-1);

/////////////////////////
// Set some heuristics //
/////////////////////////

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
param = add_param(param,'cbcheuristicfpump',0); // YC: good choice for miplib 70

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

////////////////////////
// Set Strategy class //
////////////////////////

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

// Uncomment these add_param functions to test:
// YC: trouver de meilleures valeurs par defaut.
param = add_param(param,'cbc_cutoff', 1e100);               // YC: basee sur la fonction objectif du pb relaxe ?
                                                            // YC: 1e100 (CurrentCutoff constructeur de CbcModel)
endfunction

