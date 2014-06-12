/////////////////////////////////////////////////////////////////////////////
// estim_vec_lolimot                                                       //
// For a learnt Lolimot model, computes the estimation of a several points //
/////////////////////////////////////////////////////////////////////////////

function f = estim_vec_lolimot(x,lolModel)
f = 0;
vect_phi = _phi_norm_data(x,1:size(lolModel('listofcutinf'),1),lolModel); // NbPts x NbParams
f = vect_phi' * (lolModel('listofmod')(:,1) + sum(x(:)' .*. ones(size(lolModel('listofmod'),1),1) .* lolModel('listofmod')(:,2:$),'c'));
endfunction

