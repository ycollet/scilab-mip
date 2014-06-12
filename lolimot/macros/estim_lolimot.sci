//////////////////////////////////////////////////////////////////////////
// estim_lolimot                                                        //
// For a learnt Lolimot model, computes the estimation of a given point //
//////////////////////////////////////////////////////////////////////////

function f = estim_lolimot(x,lolModel)
f = 0;
for i=1:size(lolModel('listofcutinf'),1)
  vect_phi = _phi_norm_data(x,i,lolModel); 
  f = f + vect_phi * (lolModel('listofmod')(i,1) + sum(x(:)' .* lolModel('listofmod')(i,2:$)));
end
endfunction

