function cvplot(s_opt, res_mean, res_std)
  // CVPLOT  Simple plotting function for cross validation results. 
  //    CVPLOT(S_OPT, RES_MEAN, RES_STD) plots the mean reconstruction error
  //    with error bars resulting from the function CROSSVALIDATE. The optimal
  //    model choice is marked with a dashed red line.
  //
  // Author: Karl Skoglund, IMM, DTU, kas@imm.dtu.dk

  scf();
  s_sub = linspace(0, 1, 17);
  s_sub = s_sub(2:size(s_sub,2)-1);
  t_sub = round(s_sub*length(res_mean));
  
  errbar(s_sub, res_mean(t_sub), res_std(t_sub), res_std(t_sub));
  s = linspace(0,1,length(res_mean));
  plot(s, res_mean, 'b');
  ax = gca();
  plot([s_opt ; s_opt], [ax.data_bounds(1,2) ; ax.data_bounds(2,2)], 'r-.');
endfunction
