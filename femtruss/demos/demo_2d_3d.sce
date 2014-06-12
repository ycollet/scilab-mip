test_name = '2bars';
//test_name = 'truss';
//test_name = 'dev5';
//test_name = 'bridge2d';
//test_name = 'bridge3d';
//test_name = 'bar';
//test_name = 'pylon3d';
//test_name = 'pylon2d';
//test_name = 'building3d';
//test_name = 'building2d';
//test_name = 'dome3d';

[t,p,e,A,E,rho,F] = build_fem_test(test_name);

[U,P,R,K,M]= femtruss(build_fem_test, %F, test_name);

scf();
plotdeforme(U,p,t,1);
title('Deformation');

scf();
plot_truss(p,t);
title(sprintf('Simple truss: %s', test_name));

printf('Is it a rank efficient structure: %d\n', 1.0*rank_efficient_struct(K,p,e));

[Umod,T_period,Phi]= femmode(build_fem_test, K, M, %F, 3, test_name);

scf();
for i=1:size(Phi,2)
  printf('printing mode %d/%d\n',i,size(Phi,2));
  clf();
  plotdeforme(real(matrix(Phi(:,i),size(p,1),size(p,2))),p,t,1);
  sleep(2000);
end


