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
       
status = sdpa_write_prob('test_1.sdpa',C,A,b);

printf('sdpa_write - status = %d\n', status);
  
[C_out, A_out, b_out, status_out] = sdpa_read_prob('test_1.sdpa',0);

printf('sdpa_read - status = %d\n', status_out);

printf('C_out ='); disp(C_out);
printf('A_out ='); disp(A_out);
printf('b_out ='); disp(b_out);
