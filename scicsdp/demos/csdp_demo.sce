test_full   = %t;
test_sparse = %f;

//   An example showing how to call the easy_sdp() interface to CSDP.  In
//   this example, we solve the problem
// 
//      max tr(C*X)
//          tr(A1*X)=1
//          tr(A2*X)=2
//          X >= 0       (X is PSD)
// 
//   where 
// 
//    C=[2 1
//       1 2
//           3 0 1
//           0 2 0
//           1 0 3
//                 0
//                   0]
//
//   A1=[3 1
//       1 3
//           0 0 0
//           0 0 0
//           0 0 0
//                 1
//                   0] 
//
//   A2=[0 0
//       0 0
//           3 0 1
//           0 4 0
//           1 0 5
//                 0
//                   1] 
//
//  Notice that all of the matrices have block diagonal structure.  The first
//  block is of size 2x2.  The second block is of size 3x3.  The third block
//  is a diagonal block of size 2.  

if (test_full) then
  C = list();
  C(1) = [2 1; ...
          1 2];
  C(2) = [3 0 1; ...
          0 2 0; ...
          1 0 3];
  C(3) = [0 0; ...
          0 0];
  
  A = list();
  // Constraint 1
  A(1) = list();
  A(1)(1) = [3 1; ... 
             1 3];
  A(1)(2) = [0 0 0; ...
             0 0 0; ...
             0 0 0];
  A(1)(3) = [1 0; ...
             0 0];
  
  // Constraint 2
  A(2) = list();
  A(2)(1) = [0 0; ...
             0 0];
  A(2)(2) = [3 0 1; ...
             0 4 0; ...
             1 0 5];
  A(2)(3) = [0 0; ...
             0 1];
  
  b = [1 ...
       2];
       
  params = init_param();
  
  // printlevel: from 0 to 4
  params = add_param(params,'printlevel', 4);
  params = add_param(params,'axtol', 1e-8);
  params = add_param(params,'atytol', 1e-8); 
  params = add_param(params,'objtol', 1e-8);
  params = add_param(params,'pinftol', 1e8);
  params = add_param(params,'dinftol', 1e8);
  params = add_param(params,'maxiter', 100);
  params = add_param(params,'minstepfrac', 0.90);
  params = add_param(params,'maxstepfrac', 0.97);
  params = add_param(params,'minstepp', 1e-8);
  params = add_param(params,'minstepd', 1e-8);
  params = add_param(params,'usexzgap', 1);
  params = add_param(params,'tweakgap', 0);
  params = add_param(params,'affine', 0);
  params = add_param(params,'perturbobj', 1);
  params = add_param(params,'fastmode', 0);
  params = add_param(params,'writesdpa', 'test.sdpa');
  
  [x_opt, y_opt, z_opt, f_opt, status, extra] = csdp(C,A,b,params);
  
  printf('First solution:'); disp(x_opt);
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

if (test_sparse) then
  C = list();
  C(1) = [2 1; ...
          1 2];
  C(2) = [3 0 1; ...
          0 2 0; ...
          1 0 3];
  C(3) = [0 0; ...
          0 0];
  
  A = list();
  
  // Constraint 1
  A(1) = list();
  
  A(1)(1) = sparse([3 1; ...
                    1 3]);
  
  A(1)(2) = sparse([0 0 0; ...
                    0 0 0; ...
                    0 0 0]);
  
  A(1)(3) = sparse([1 0; ...
                    0 0]);
  
  // Constraint 2
  A(2) = list();
  
  A(2)(1) = sparse([0 0; ...
                    0 0]);
  
  A(2)(2) = sparse([3 0 1; ...
                    0 4 0; ...
                    1 0 5]);
  
  A(2)(3) = sparse([0 0; ...
                    0 1]);
  
  b = [1 ...
       2];
       
  params = init_param();
  
  // printlevel: from 0 to 4
  params = add_param(params,'printlevel', 4);
  params = add_param(params,'axtol', 1e-8);
  params = add_param(params,'atytol', 1e-8); 
  params = add_param(params,'objtol', 1e-8);
  params = add_param(params,'pinftol', 1e8);
  params = add_param(params,'dinftol', 1e8);
  params = add_param(params,'maxiter', 100);
  params = add_param(params,'minstepfrac', 0.90);
  params = add_param(params,'maxstepfrac', 0.97);
  params = add_param(params,'minstepp', 1e-8);
  params = add_param(params,'minstepd', 1e-8);
  params = add_param(params,'usexzgap', 1);
  params = add_param(params,'tweakgap', 0);
  params = add_param(params,'affine', 0);
  params = add_param(params,'perturbobj', 1);
  params = add_param(params,'fastmode', 0);
  params = add_param(params,'writesdpa', 'test.sdpa');

  [x_opt, y_opt, z_opt, f_opt, status, extra] = csdp(C,A,b,params);
  
  printf('Second solution:'); disp(x_opt);
  
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
