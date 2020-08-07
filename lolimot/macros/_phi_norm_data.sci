/////////////////////////////////////////////////////////////////////
// _phi_norm_data                                                  //
//                                                                 //
// For a partition, computes the membership level of a given point //
/////////////////////////////////////////////////////////////////////

function vect_phi = _phi_norm_data(x,nopart,lolModel)
// data: data(:,1:$-1) les entrees, data(:,$) la sortie
for i=1:size(lolModel('listofcutinf'),1)
  xc = (lolModel('listofcutinf')(i,:) + lolModel('listofcutplus')(i,:))/2;
  dx = lolModel('listofcutplus')(i,:) - lolModel('listofcutinf')(i,:);
  phi_tmp(i) = exp(-sum(((x(:)-xc(:)).^2)./(2*lolModel('sigma')^2*dx(:).^2)))
end
vect_phi = phi_tmp(nopart) / max([sum(phi_tmp) %eps]);
endfunction

