stacksize('max');

// These files doesn't seems to work fine ...
// 'arch0.dat-s'; 'arch2.dat-s'; 'arch4.dat-s'; 'arch8.dat-s'; ...
         
Files = ['control10.dat-s'; 'control11.dat-s'; 'control1.dat-s'; 'control2.dat-s'; ...
	 'control3.dat-s'; 'control4.dat-s'; 'control5.dat-s'; 'control6.dat-s'; ...
	 'control7.dat-s'; 'control8.dat-s'; 'control9.dat-s'; ...
	 'equalG11.dat-s'; 'equalG51.dat-s'; ...
	 'gpp100.dat-s'; 'gpp124-1.dat-s'; 'gpp124-2.dat-s'; 'gpp124-3.dat-s'; ...
 	 'gpp124-4.dat-s'; 'gpp250-1.dat-s'; 'gpp250-2.dat-s'; 'gpp250-3.dat-s'; ...
	 'gpp250-4.dat-s'; 'gpp500-1.dat-s'; 'gpp500-2.dat-s'; 'gpp500-3.dat-s'; ...
	 'gpp500-4.dat-s'; ...
	 'hinf10.dat-s'; 'hinf11.dat-s'; 'hinf12.dat-s'; 'hinf13.dat-s'; ...
	 'hinf14.dat-s'; 'hinf15.dat-s'; 'hinf1.dat-s'; 'hinf2.dat-s'; ...
	 'hinf3.dat-s'; 'hinf4.dat-s'; 'hinf5.dat-s'; 'hinf6.dat-s'; ...
	 'hinf7.dat-s'; 'hinf8.dat-s'; 'hinf9.dat-s'; ...
	 'infd1.dat-s'; 'infd2.dat-s'; 'infp1.dat-s'; 'infp2.dat-s'; ...
	 'maxG11.dat-s'; 'maxG32.dat-s'; 'maxG51.dat-s'; 'maxG55.dat-s'; ...
	 'maxG60.dat-s'; ...
	 'mcp100.dat-s'; 'mcp124-1.dat-s'; 'mcp124-2.dat-s'; 'mcp124-3.dat-s'; ...
	 'mcp124-4.dat-s'; 'mcp250-1.dat-s'; 'mcp250-2.dat-s'; 'mcp250-3.dat-s'; ...
	 'mcp250-4.dat-s'; 'mcp500-1.dat-s'; 'mcp500-2.dat-s'; 'mcp500-3.dat-s'; ...
	 'mcp500-4.dat-s'; ...
	 'qap10.dat-s'; 'qap5.dat-s'; 'qap6.dat-s'; 'qap7.dat-s'; 'qap8.dat-s'; ...
	 'qap9.dat-s'; 'qpG11.dat-s'; 'qpG51.dat-s'; ...
	 'ss30.dat-s'; ...
	 'theta1.dat-s'; 'theta2.dat-s'; 'theta3.dat-s'; 'theta4.dat-s'; 'theta5.dat-s'; ...
	 'theta6.dat-s'; 'thetaG11.dat-s'; 'thetaG51.dat-s'; ...
	 'truss1.dat-s'; 'truss2.dat-s'; 'truss3.dat-s'; 'truss4.dat-s'; 'truss5.dat-s'; ...
	 'truss6.dat-s'; 'truss7.dat-s'; 'truss8.dat-s'];

n = x_choose(Files,['Choose the SDPA file'],'Cancel');

if (n==0) then abort(); end

filename = 'data/sdpa/' + Files(n);

[C_out, A_out, b_out, status_out] = sdpa_read_prob(filename,0);

printf('File %s read - status = %d\n', filename, status_out);

printf('C_out - nb blocks = %d\n', length(C_out));
printf('A_out - nb blocks = %d\n', length(A_out));
printf('b_out - nb blocks = %d\n', length(b_out));

res = messagebox("Would you like to solve this problem ?", "modal", "info", ["Yes", "No"]);

// Yes == 1, No == 2

if res==1 then
  params = init_param();
  
  // printlevel: from 0 to 4
  params = add_param(params, 'printlevel',  4);
  params = add_param(params, 'axtol',       1e-6); // 1e-8
  params = add_param(params, 'atytol',      1e-6); // 1e-8
  params = add_param(params, 'objtol',      1e-6); // 1e-8
  params = add_param(params, 'pinftol',     1e10); // 1e8
  params = add_param(params, 'dinftol',     1e10); // 1e8
  params = add_param(params, 'maxiter',     100);
  params = add_param(params, 'minstepfrac', 0.0); // 0.90
  params = add_param(params, 'maxstepfrac', 0.97);
  params = add_param(params, 'minstepp',    1e-6); // 1e-8
  params = add_param(params, 'minstepd',    1e-6); // 1e-8
  params = add_param(params, 'usexzgap',    1);
  params = add_param(params, 'tweakgap',    0);
  params = add_param(params, 'affine',      0);
  params = add_param(params, 'perturbobj',  1);
  params = add_param(params, 'fastmode',    0);
  //params = add_param(params, 'writesdpa',   'test.sdpa');
  
  [x_opt, y_opt, z_opt, f_opt, status, extra] = csdp(C_out, A_out, b_out, params);
  
  printf('First solution:'); disp(x_opt);

  printf('After solving problem %s\n', Files(n));
  
  printf('Status = %d\n', status);
  printf('0 Success\n');
  printf('1 Success: The problem is primal infeasibile, and we have a certificate.\n');
  printf('2 Success: The problem is dual infeasible, and we have a certificate.\n');
  printf('3 Partial Success: Didn''t reach full accuracy.\n');
  printf('4 Failure: Maximum iterations reached.\n');
  printf('5 Failure: Stuck at edge of primal feasibility.\n');
  printf('6 Failure: Stuck at edge of dual feasibility.\n');
  printf('7 Failure: Lack of progress\n');
  printf('8 Failure: X, Z, or O was singular.\n');
  printf('9 Failure: Detected NaN or Inf values.\n');
end
