// example of use of the simulated annealing method

old_funcprot = funcprot();
funcprot(0);

// filenames = ['a280.tsp','a280.opt.tour'];
// filenames = ['berlin52.tsp','berlin52.opt.tour'];
// filenames = ['bier127.tsp',''];
// filenames = ['brd14051.tsp',''];
filenames = ['ch130.tsp','ch130.opt.tour'];
// filenames = ['ch150.tsp','ch150.opt.tour'];
// filenames = ['d1291.tsp',''];
// filenames = ['d15112.tsp',''];
// filenames = ['d1655.tsp',''];
// filenames = ['d18512.tsp',''];
// filenames = ['d198.tsp',''];
// filenames = ['d2103.tsp',''];
// filenames = ['d493.tsp',''];
// filenames = ['d657.tsp',''];
// filenames = ['eil101.tsp','eil101.opt.tour'];
// filenames = ['eil51.tsp','eil51.opt.tour'];
// filenames = ['eil76.tsp','eil76.opt.tour'];
// filenames = ['fl1400.tsp',''];
// filenames = ['fl1577.tsp',''];
// filenames = ['fl3795.tsp',''];
// filenames = ['fl417.tsp',''];
// filenames = ['fnl4461.tsp',''];
// filenames = ['gil262.tsp',''];
// filenames = ['kroA100.tsp','kroA100.opt.tour'];
// filenames = ['kroA150.tsp',''];
// filenames = ['kroA200.tsp',''];
// filenames = ['kroB100.tsp',''];
// filenames = ['kroB150.tsp',''];
// filenames = ['kroB200.tsp',''];
// filenames = ['kroC100.tsp','kroC100.opt.tour'];
// filenames = ['kroD100.tsp','kroD100.opt.tour'];
// filenames = ['kroE100.tsp',''];
// filenames = ['lin105.tsp','lin105.opt.tour'];
// filenames = ['lin318.tsp',''];
// filenames = ['linhp318.tsp',''];
// filenames = ['nrw1379.tsp',''];
// filenames = ['p654.tsp',''];
// filenames = ['pcb1173.tsp',''];
// filenames = ['pcb3038.tsp',''];
// filenames = ['pcb442.tsp','pcb442.opt.tour'];
// filenames = ['pr1002.tsp',''];
// filenames = ['pr107.tsp',''];
// filenames = ['pr124.tsp',''];
// filenames = ['pr136.tsp',''];
// filenames = ['pr144.tsp',''];
// filenames = ['pr152.tsp',''];
// filenames = ['pr226.tsp',''];
// filenames = ['pr2392.tsp','pr2392.opt.tour'];
// filenames = ['pr264.tsp',''];
// filenames = ['pr299.tsp',''];
// filenames = ['pr439.tsp',''];
// filenames = ['pr76.tsp','pr76.opt.tour'];
// filenames = ['rat195.tsp',''];
// filenames = ['rat575.tsp',''];
// filenames = ['rat783.tsp',''];
// filenames = ['rat99.tsp',''];
// filenames = ['rd100.tsp','rd100.opt.tour'];
// filenames = ['rd400.tsp',''];
// filenames = ['rl11849.tsp',''];
// filenames = ['rl1304.tsp',''];
// filenames = ['rl1323.tsp',''];
// filenames = ['rl1889.tsp',''];
// filenames = ['rl5915.tsp',''];
// filenames = ['rl5934.tsp',''];
// filenames = ['st70.tsp','st70.opt.tour'];
// filenames = ['ts225.tsp',''];
// filenames = ['tsp225.tsp','tsp225.opt.tour'];
// filenames = ['u1060.tsp',''];
// filenames = ['u1432.tsp',''];
// filenames = ['u159.tsp',''];
// filenames = ['u1817.tsp',''];
// filenames = ['u2152.tsp',''];
// filenames = ['u2319.tsp',''];
// filenames = ['u574.tsp',''];
// filenames = ['u724.tsp',''];
// filenames = ['usa13509.tsp',''];
// filenames = ['vm1084.tsp',''];
// filenames = ['vm1748.tsp',''];

Proba_start = 0.8;
It_intern   = 1000;
It_extern   = 30;
It_Pre      = 100;
It_Scramble = 1000;
DoSA    = %F;
DoHuang = %T;

//////////////////////////////////////////

path = get_absolute_file_path('SAdemo.sce');

[TSP_Cycle,TSPSize,TSPComment]  = read_tsp(path + '../tours/'+filenames(1));
TSP0   = init_tsp_var(TSP_Cycle);
TSPMat = compute_tsp_dist(TSP_Cycle,%F);
TSP0   = scramble_tsp(TSP0,It_Scramble);

deff('y=f(x)','y=compute_tsp_fobj(TSPMat,x)');

/////////////////////////
// Simulated Annealing //
/////////////////////////

if DoSA then
  printf('SA: geometrical decrease temperature law\n');

  sa_params = init_param();
  sa_params = add_param(sa_params,'accept_func', accept_func_default); // Optional
  sa_params = add_param(sa_params,'temp_law', temp_law_default); // Optional
  //sa_params = add_param(sa_params,'neigh_func', neigh_func_tsp_swap); 
  sa_params = add_param(sa_params,'neigh_func', neigh_func_tsp_2opt); 
  sa_params = add_param(sa_params,'alpha', 0.8); 

  T0 = compute_initial_temp(TSP0, f, Proba_start, It_Pre, sa_params);
  printf('Initial temperatore T0 = %f\n', T0);

  [tsp_opt, f_opt, sa_mean_list, sa_var_list, temp_list] = optim_sa(TSP0, f, It_extern, It_intern, T0, Log = %T, sa_params);

  plot_tsp(TSP_Cycle,TSP0,'Initial cycle');
  plot_tsp(TSP_Cycle,tsp_opt,'Final cycle');
  
  printf('value of the objective function for the initial cycle = %f\n', compute_tsp_fobj(TSPMat,TSP0));
  printf('value of the objective function for the final cycle = %f\n', f_opt);

  scf();
  drawlater;
  subplot(4,1,1);
  xtitle('Huang annealing','Iteration','Mean');
  t = 1:length(sa_mean_list);
  plot(t,sa_mean_list,'r');
  subplot(4,1,2);
  xtitle('','Iteration','Variance');
  t = 1:length(sa_mean_list);
  plot(t,sa_var_list,'g');
  subplot(2,1,2);
  xtitle('Temperature evolution','Iteration','Temperature');
  for i=1:length(t)-1
    plot([t(i) t(i+1)], [temp_list(i) temp_list(i)],'k-');
  end
  drawnow;
  
  if filenames(2)~='' then
    TSP_Opt = read_tsp_opt(path + '../tours/'+filenames(2));
    plot_tsp(TSP_Cycle,TSP_Opt,'Optimal Cycle');
    printf('value of the objective function for the optimal cycle = %f\n', compute_tsp_fobj(TSPMat,TSP_Opt));
  end
end

///////////
// Huang //
///////////

if DoHuang then
  printf('SA: the Huang annealing\n');

  sa_params = init_param();
  sa_params = add_param(sa_params,'temp_law', temp_law_huang);
  sa_params = add_param(sa_params,'lambda', 0.5);
  //sa_params = add_param(sa_params,'neigh_func', neigh_func_tsp_swap); 
  sa_params = add_param(sa_params,'neigh_func', neigh_func_tsp_2opt); 

  T0 = compute_initial_temp(TSP0, f, Proba_start, It_Pre, sa_params);
  printf('Initial temperatore T0 = %f\n', T0);

  [tsp_opt, f_opt, sa_mean_list, sa_var_list, temp_list] = optim_sa(TSP0, f, It_extern, It_intern, T0, Log = %T, sa_params);

  plot_tsp(TSP_Cycle,TSP0,'Initial cycle');
  plot_tsp(TSP_Cycle,tsp_opt,'Final cycle');
  
  printf('value of the objective function for the initial cycle = %f\n', compute_tsp_fobj(TSPMat,TSP0));
  printf('value of the objective function for the final cycle = %f\n', f_opt);

  scf();
  drawlater;
  subplot(4,1,1);
  xtitle('Huang annealing','Iteration','Mean');
  t = 1:length(sa_mean_list);
  plot(t,sa_mean_list,'r');
  subplot(4,1,2);
  xtitle('','Iteration','Variance');
  t = 1:length(sa_mean_list);
  plot(t,sa_var_list,'g');
  subplot(2,1,2);
  xtitle('Temperature evolution','Iteration','Temperature');
  for i=1:length(t)-1
    plot([t(i) t(i+1)], [temp_list(i) temp_list(i)],'k-');
  end
  drawnow;

  if filenames(2)~='' then
    TSP_Opt = read_tsp_opt(path + '../tours/'+filenames(2));
    plot_tsp(TSP_Cycle,TSP_Opt,'Optimal Cycle');
    printf('value of the objective function for the optimal cycle = %f\n', compute_tsp_fobj(TSPMat,TSP_Opt));
  end
end

funcprot(old_funcprot);

