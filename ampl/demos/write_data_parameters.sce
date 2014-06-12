lines(0);

stacksize('max');

exec nl_data.sce;

write_html = %T;

if write_html then
  fid = mopen('listoffiles.html','w+');

  mfprintf(fid, '\nProblem informations:<P>\n');
  mfprintf(fid, 'number of linear binary variables:     nbv<P>\n');
  mfprintf(fid, 'number of linear integer variables:    niv<P>\n');
  mfprintf(fid, 'total number of nonlinear constraints: nlc<P>\n');
  mfprintf(fid, 'number of equality constraints or -1 if unknown (ampl prior to 19970627) : n_eqn<P>\n');
  mfprintf(fid, 'total complementarity conditions:        n_cc<P>\n');
  mfprintf(fid, 'nonlinear complementarity conditions:    nlcc<P>\n');
  mfprintf(fid, 'number of nonlinear network constraints: nlnc<P>\n');
  mfprintf(fid, 'number of nonlinear objectives:          nlo<P>\n');
  mfprintf(fid, 'number of nonlinear variables in both constraints and objectives: nlvb<P>\n');
  mfprintf(fid, 'number of nonlinear variables in constraints:                     nlvc<P>\n');
  mfprintf(fid, 'number of nonlinear variables in objectives nlvc and nlvo include nlvb: nlvo<P>\n'); 
  mfprintf(fid, 'integer nonlinear variables in both constraints and objectives :        nlvbi<P>\n');
  mfprintf(fid, 'integer nonlinear vars just in constraints :    nlvci<P>\n');
  mfprintf(fid, 'integer nonlinear vars just in objectives:      nlvoi<P>\n');
  mfprintf(fid, 'number of (linear) network variables (arcs):    nwv<P>\n');
  mfprintf(fid, 'number of nonzeros in constraints Jacobian:     nzc<P>\n');
  mfprintf(fid, 'number of nonzeros in all objective gradients : nzo<P>\n');
  mfprintf(fid, 'total number of variables:     n_var<P>\n');
  mfprintf(fid, 'total number of constraints:   n_con<P>\n');
  mfprintf(fid, 'total number of objectives:    n_obj<P>\n');
  mfprintf(fid, 'number of logical constraints: n_lcon<P>\n'); 
  
  mfprintf(fid, '<TABLE BORDER=""1"">\n');
  mfprintf(fid, '<CAPTION> Characteritics of NL files </CAPTION>\n');
  mfprintf(fid, '<TR>\n');
  mfprintf(fid, ' <TH> filename <H>\n');
  mfprintf(fid, ' <TH> index <H>\n');
  mfprintf(fid, ' <TH> nbv <H>\n');
  mfprintf(fid, ' <TH> niv <H>\n');
  mfprintf(fid, ' <TH> nlc <H>\n');
  mfprintf(fid, ' <TH> n_eqn <H>\n');
  mfprintf(fid, ' <TH> n_cc <H>\n');
  mfprintf(fid, ' <TH> nlcc <H>\n');
  mfprintf(fid, ' <TH> nlnc <H>\n');
  mfprintf(fid, ' <TH> nlo <H>\n');
  mfprintf(fid, ' <TH> nlvb <H>\n');
  mfprintf(fid, ' <TH> nlvc <H>\n');
  mfprintf(fid, ' <TH> nlvo <H>\n'); 
  mfprintf(fid, ' <TH> nlvbi <H>\n');
  mfprintf(fid, ' <TH> nlvci <H>\n');
  mfprintf(fid, ' <TH> nlvoi <H>\n');
  mfprintf(fid, ' <TH> nwv <H>\n');
  mfprintf(fid, ' <TH> nzc <H>\n');
  mfprintf(fid, ' <TH> nzo <H>\n');
  mfprintf(fid, ' <TH> n_var <H>\n');
  mfprintf(fid, ' <TH> n_con <H>\n');
  mfprintf(fid, ' <TH> n_obj <H>\n');
  mfprintf(fid, ' <TH> n_lcon <H>\n'); 
  mfprintf(fid, '<R>\n');

else
  fid = mopen('listoffiles.txt','w+');

  mfprintf(fid, '\nProblem informations:\n');
  mfprintf(fid, 'number of linear binary variables:     nbv\n');
  mfprintf(fid, 'number of linear integer variables:    niv\n');
  mfprintf(fid, 'total number of nonlinear constraints: nlc\n');
  mfprintf(fid, 'number of equality constraints or -1 if unknown (ampl prior to 19970627) : n_eqn\n');
  mfprintf(fid, 'total complementarity conditions:        n_cc\n');
  mfprintf(fid, 'nonlinear complementarity conditions:    nlcc\n');
  mfprintf(fid, 'number of nonlinear network constraints: nlnc\n');
  mfprintf(fid, 'number of nonlinear objectives:          nlo\n');
  mfprintf(fid, 'number of nonlinear variables in both constraints and objectives: nlvb\n');
  mfprintf(fid, 'number of nonlinear variables in constraints:                     nlvc\n');
  mfprintf(fid, 'number of nonlinear variables in objectives nlvc and nlvo include nlvb: nlvo\n'); 
  mfprintf(fid, 'integer nonlinear variables in both constraints and objectives :        nlvbi\n');
  mfprintf(fid, 'integer nonlinear vars just in constraints :    nlvci\n');
  mfprintf(fid, 'integer nonlinear vars just in objectives:      nlvoi\n');
  mfprintf(fid, 'number of (linear) network variables (arcs):    nwv\n');
  mfprintf(fid, 'number of nonzeros in constraints Jacobian:     nzc\n');
  mfprintf(fid, 'number of nonzeros in all objective gradients : nzo\n');
  mfprintf(fid, 'total number of variables:     n_var\n');
  mfprintf(fid, 'total number of constraints:   n_con\n');
  mfprintf(fid, 'total number of objectives:    n_obj\n');
  mfprintf(fid, 'number of logical constraints: n_lcon\n'); 

  mfprintf(fid, 'Characteritics of NL files\n');
  mfprintf(fid, '\n');
  mfprintf(fid, '|| filename || index || nbv || niv || nlc || n_eqn || n_cc || nlcc || nlnc || nlo || nlvb || nlvc || nlvo'); 
  mfprintf(fid, ' || nlvbi || nlvci || nlvoi || nwv || nzc || nzo || n_var || n_con || n_obj || n_lcon ||\n'); 
end

//////////////
// MacMINLP //
//////////////

for i=1:length(MacMINLP)
  nl_filename = MacMINLP(i);
  
  printf('processing file %s\n', nl_filename);
  
  ///////////////////
  // Read the data //
  ///////////////////

  [asl, x0, bl, bu, v, cl, cu] = ampl_init(nl_filename);

  info = ampl_get_size(asl);

  if write_html then
    mfprintf(fid, '<TR>\n');
    mfprintf(fid, ' <TH> %s <H>\n',nl_filename);
    mfprintf(fid, ' <TH> %d <H>\n',i);
    mfprintf(fid, ' <TH> %d <H>\n',info('nbv'));
    mfprintf(fid, ' <TH> %d <H>\n',info('niv'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlc'));
    mfprintf(fid, ' <TH> %d <H>\n',info('n_eqn'));
    mfprintf(fid, ' <TH> %d <H>\n',info('n_cc'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlcc'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlnc'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlo'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlvb'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlvc'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlvo')); 
    mfprintf(fid, ' <TH> %d <H>\n',info('nlvbi'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlvci'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlvoi'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nwv'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nzc'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nzo'));
    mfprintf(fid, ' <TH> %d <H>\n',info('n_var'));
    mfprintf(fid, ' <TH> %d <H>\n',info('n_con'));
    mfprintf(fid, ' <TH> %d <H>\n',info('n_obj'));
    mfprintf(fid, ' <TH> %d <H>\n',info('n_lcon')); 
    mfprintf(fid, '<R>\n');
  else
    mfprintf(fid, '|| %s ||',nl_filename);
    mfprintf(fid, ' %d ||',i);
    mfprintf(fid, ' %d ||',info('nbv'));
    mfprintf(fid, ' %d ||',info('niv'));
    mfprintf(fid, ' %d ||',info('nlc'));
    mfprintf(fid, ' %d ||',info('n_eqn'));
    mfprintf(fid, ' %d ||',info('n_cc'));
    mfprintf(fid, ' %d ||',info('nlcc'));
    mfprintf(fid, ' %d ||',info('nlnc'));
    mfprintf(fid, ' %d ||',info('nlo'));
    mfprintf(fid, ' %d ||',info('nlvb'));
    mfprintf(fid, ' %d ||',info('nlvc'));
    mfprintf(fid, ' %d ||',info('nlvo')); 
    mfprintf(fid, ' %d ||',info('nlvbi'));
    mfprintf(fid, ' %d ||',info('nlvci'));
    mfprintf(fid, ' %d ||',info('nlvoi'));
    mfprintf(fid, ' %d ||',info('nwv'));
    mfprintf(fid, ' %d ||',info('nzc'));
    mfprintf(fid, ' %d ||',info('nzo'));
    mfprintf(fid, ' %d ||',info('n_var'));
    mfprintf(fid, ' %d ||',info('n_con'));
    mfprintf(fid, ' %d ||',info('n_obj'));
    mfprintf(fid, ' %d ||\n',info('n_lcon')); 
  end

  ampl_free(asl);
end

////////////
// CoinOR //
////////////

for i=1:length(CoinOR)
  nl_filename = CoinOR(i);
  
  printf('processing file %s\n', nl_filename);
  
  ///////////////////
  // Read the data //
  ///////////////////

  [asl, x0, bl, bu, v, cl, cu] = ampl_init(nl_filename);

  info = ampl_get_size(asl);

  if write_html then
    mfprintf(fid, '<TR>\n');
    mfprintf(fid, ' <TH> %s <H>\n',nl_filename);
    mfprintf(fid, ' <TH> %d <H>\n',i);
    mfprintf(fid, ' <TH> %d <H>\n',info('nbv'));
    mfprintf(fid, ' <TH> %d <H>\n',info('niv'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlc'));
    mfprintf(fid, ' <TH> %d <H>\n',info('n_eqn'));
    mfprintf(fid, ' <TH> %d <H>\n',info('n_cc'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlcc'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlnc'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlo'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlvb'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlvc'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlvo')); 
    mfprintf(fid, ' <TH> %d <H>\n',info('nlvbi'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlvci'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlvoi'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nwv'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nzc'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nzo'));
    mfprintf(fid, ' <TH> %d <H>\n',info('n_var'));
    mfprintf(fid, ' <TH> %d <H>\n',info('n_con'));
    mfprintf(fid, ' <TH> %d <H>\n',info('n_obj'));
    mfprintf(fid, ' <TH> %d <H>\n',info('n_lcon')); 
    mfprintf(fid, '<R>\n');
  else
    mfprintf(fid, '|| %s ||',nl_filename);
    mfprintf(fid, ' %d ||',i);
    mfprintf(fid, ' %d ||',info('nbv'));
    mfprintf(fid, ' %d ||',info('niv'));
    mfprintf(fid, ' %d ||',info('nlc'));
    mfprintf(fid, ' %d ||',info('n_eqn'));
    mfprintf(fid, ' %d ||',info('n_cc'));
    mfprintf(fid, ' %d ||',info('nlcc'));
    mfprintf(fid, ' %d ||',info('nlnc'));
    mfprintf(fid, ' %d ||',info('nlo'));
    mfprintf(fid, ' %d ||',info('nlvb'));
    mfprintf(fid, ' %d ||',info('nlvc'));
    mfprintf(fid, ' %d ||',info('nlvo')); 
    mfprintf(fid, ' %d ||',info('nlvbi'));
    mfprintf(fid, ' %d ||',info('nlvci'));
    mfprintf(fid, ' %d ||',info('nlvoi'));
    mfprintf(fid, ' %d ||',info('nwv'));
    mfprintf(fid, ' %d ||',info('nzc'));
    mfprintf(fid, ' %d ||',info('nzo'));
    mfprintf(fid, ' %d ||',info('n_var'));
    mfprintf(fid, ' %d ||',info('n_con'));
    mfprintf(fid, ' %d ||',info('n_obj'));
    mfprintf(fid, ' %d ||\n',info('n_lcon')); 
  end
  
  ampl_free(asl);
end

/////////
// ASL //
/////////

for i=1:length(ASL)
  nl_filename = ASL(i);
  
  printf('processing file %s\n', nl_filename);
  
  ///////////////////
  // Read the data //
  ///////////////////

  [asl, x0, bl, bu, v, cl, cu] = ampl_init(nl_filename);

  info = ampl_get_size(asl);

  if write_html then
    mfprintf(fid, '<TR>\n');
    mfprintf(fid, ' <TH> %s <H>\n',nl_filename);
    mfprintf(fid, ' <TH> %d <H>\n',i);
    mfprintf(fid, ' <TH> %d <H>\n',info('nbv'));
    mfprintf(fid, ' <TH> %d <H>\n',info('niv'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlc'));
    mfprintf(fid, ' <TH> %d <H>\n',info('n_eqn'));
    mfprintf(fid, ' <TH> %d <H>\n',info('n_cc'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlcc'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlnc'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlo'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlvb'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlvc'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlvo')); 
    mfprintf(fid, ' <TH> %d <H>\n',info('nlvbi'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlvci'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlvoi'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nwv'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nzc'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nzo'));
    mfprintf(fid, ' <TH> %d <H>\n',info('n_var'));
    mfprintf(fid, ' <TH> %d <H>\n',info('n_con'));
    mfprintf(fid, ' <TH> %d <H>\n',info('n_obj'));
    mfprintf(fid, ' <TH> %d <H>\n',info('n_lcon')); 
    mfprintf(fid, '<R>\n');
  else
    mfprintf(fid, '|| %s ||',nl_filename);
    mfprintf(fid, ' %d ||',i);
    mfprintf(fid, ' %d ||',info('nbv'));
    mfprintf(fid, ' %d ||',info('niv'));
    mfprintf(fid, ' %d ||',info('nlc'));
    mfprintf(fid, ' %d ||',info('n_eqn'));
    mfprintf(fid, ' %d ||',info('n_cc'));
    mfprintf(fid, ' %d ||',info('nlcc'));
    mfprintf(fid, ' %d ||',info('nlnc'));
    mfprintf(fid, ' %d ||',info('nlo'));
    mfprintf(fid, ' %d ||',info('nlvb'));
    mfprintf(fid, ' %d ||',info('nlvc'));
    mfprintf(fid, ' %d ||',info('nlvo')); 
    mfprintf(fid, ' %d ||',info('nlvbi'));
    mfprintf(fid, ' %d ||',info('nlvci'));
    mfprintf(fid, ' %d ||',info('nlvoi'));
    mfprintf(fid, ' %d ||',info('nwv'));
    mfprintf(fid, ' %d ||',info('nzc'));
    mfprintf(fid, ' %d ||',info('nzo'));
    mfprintf(fid, ' %d ||',info('n_var'));
    mfprintf(fid, ' %d ||',info('n_con'));
    mfprintf(fid, ' %d ||',info('n_obj'));
    mfprintf(fid, ' %d ||\n',info('n_lcon')); 
  end
  
  ampl_free(asl);
end

///////////
// ModNL //
///////////

for i=1:length(ModNL)
  nl_filename = ModNL(i);
  
  printf('processing file %s\n', nl_filename);
  
  ///////////////////
  // Read the data //
  ///////////////////

  [asl, x0, bl, bu, v, cl, cu] = ampl_init(nl_filename);

  info = ampl_get_size(asl);

  if write_html then
    mfprintf(fid, '<TR>\n');
    mfprintf(fid, ' <TH> %s <H>\n',nl_filename);
    mfprintf(fid, ' <TH> %d <H>\n',i);
    mfprintf(fid, ' <TH> %d <H>\n',info('nbv'));
    mfprintf(fid, ' <TH> %d <H>\n',info('niv'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlc'));
    mfprintf(fid, ' <TH> %d <H>\n',info('n_eqn'));
    mfprintf(fid, ' <TH> %d <H>\n',info('n_cc'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlcc'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlnc'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlo'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlvb'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlvc'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlvo')); 
    mfprintf(fid, ' <TH> %d <H>\n',info('nlvbi'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlvci'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nlvoi'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nwv'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nzc'));
    mfprintf(fid, ' <TH> %d <H>\n',info('nzo'));
    mfprintf(fid, ' <TH> %d <H>\n',info('n_var'));
    mfprintf(fid, ' <TH> %d <H>\n',info('n_con'));
    mfprintf(fid, ' <TH> %d <H>\n',info('n_obj'));
    mfprintf(fid, ' <TH> %d <H>\n',info('n_lcon')); 
    mfprintf(fid, '<R>\n');
  else
    mfprintf(fid, '|| %s ||',nl_filename);
    mfprintf(fid, ' %d ||',i);
    mfprintf(fid, ' %d ||',info('nbv'));
    mfprintf(fid, ' %d ||',info('niv'));
    mfprintf(fid, ' %d ||',info('nlc'));
    mfprintf(fid, ' %d ||',info('n_eqn'));
    mfprintf(fid, ' %d ||',info('n_cc'));
    mfprintf(fid, ' %d ||',info('nlcc'));
    mfprintf(fid, ' %d ||',info('nlnc'));
    mfprintf(fid, ' %d ||',info('nlo'));
    mfprintf(fid, ' %d ||',info('nlvb'));
    mfprintf(fid, ' %d ||',info('nlvc'));
    mfprintf(fid, ' %d ||',info('nlvo')); 
    mfprintf(fid, ' %d ||',info('nlvbi'));
    mfprintf(fid, ' %d ||',info('nlvci'));
    mfprintf(fid, ' %d ||',info('nlvoi'));
    mfprintf(fid, ' %d ||',info('nwv'));
    mfprintf(fid, ' %d ||',info('nzc'));
    mfprintf(fid, ' %d ||',info('nzo'));
    mfprintf(fid, ' %d ||',info('n_var'));
    mfprintf(fid, ' %d ||',info('n_con'));
    mfprintf(fid, ' %d ||',info('n_obj'));
    mfprintf(fid, ' %d ||\n',info('n_lcon')); 
  end
  
  ampl_free(asl);
end

/////////////////////
// End of the file //
/////////////////////

if write_html then
  mfprintf(fid,'<TABLE>\n');
end

mclose(fid);

