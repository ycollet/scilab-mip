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

path_new = get_absolute_file_path('mps_osi_test.sce');
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
Solve_netlib   = %T; // 91 files
Solve_data     = %F; // 2 files
Solve_optimslp = %F; // 3 files

///////////////////////////
// Set global parameters //
///////////////////////////

Log       = 1;
Debug     = 0;
mps_index = 32; // 1 // 1 LP - 2 MIP
// netlib index 5: fobj = nan apres 1000 iterations avec interior
mps_type  = 0;
UseLinpro = %F;
PrintSol  = %F;
PreSolve  = 1;
Scale     = 1;

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
param = add_param(param,'maxnumiteration',10000);
param = add_param(param,'primaltolerance',1e-7);
param = add_param(param,'dualtolerance',1e-7);
param = add_param(param,'verbose',Log);
param = add_param(param,'solvername','clp');
param = add_param(param,'writemps','test.mps');
param = add_param(param,'optim_dir', 1); // optimisation direction: 1 - minimize, -1 - maximize

// scaling option
// Sets or unsets scaling:
// - 0 -off
// - 1 equilibrium (default)
// - 2 geometric
// - 3 auto
// - 4 auto-but-as-initialSolve-in-bab. 
param = add_param(param,'doscale',0); // Integer

param = add_param(param,'maxnumiterationshotstart',1000);
param = add_param(param,'dualobjectivelimit',1e100);
param = add_param(param,'primalobjectivelimit',1e100);
//param = add_param(param,'objoffset',0.0);

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

[xmin,fmin,status,extra] = osi(mps_file_mp('obj_coeff'),mps_file_mp('constr_mat'),mps_file_mp('lhs'),mps_file_mp('rhs'), ...
                               mps_file_mp('bounds_lower'),mps_file_mp('bounds_upper'),mps_file_mp('constr_sense'),var_type,param);

printf('nb of constr  = %d\n', size(mps_file_mp('constr_mat'),1));
printf('nb of var     = %d\n', size(mps_file_mp('constr_mat'),2));
printf('nb of int var = %d\n', mps_file('nb_int_var'));
printf('status: %s\n', dec2bin(status));
printf(' - bit 1: isAbandoned\n');
printf(' - bit 2: isProvenOptimal\n');
printf(' - bit 3: isProvenPrimalInfeasible\n');
printf(' - bit 4: isProvenDualInfeasible\n');
printf(' - bit 5: isPrimalObjectiveLimitReached\n');
printf(' - bit 6: isDualObjectiveLimitReached\n');
printf(' - bit 7: isIterationLimitReached\n');

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
