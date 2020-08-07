/////////////////////////////////////////////////////////////////////
// _phi_norm_data_linear                                           //
//                                                                 //
// For a partition, computes a piecewise linear membership level   //
// of a given point                                                //
/////////////////////////////////////////////////////////////////////

function vect_phi = _phi_norm_data_linear(x,nopart,lolModel)
k_overlap = 1.0/lolModel('sigma'); // 0.5 avant
// data: data(:,1:$-1) the inputs, data(:,$) the output
// Use a piecewise linear membership function instead of a gaussian one
for i=1:size(lolModel('listofcutinf'),1)
  xc = (lolModel('listofcutinf')(i,:) + lolModel('listofcutplus')(i,:))/2;
  dx = lolModel('listofcutplus')(i,:) - lolModel('listofcutinf')(i,:);
  
  phi_tmp(i) = 1;
  
  for j=1:length(xc)
    if ((xc(j) - lolModel('sigma')*dx(j)/2 <= x(j))&(x(j) <= xc(j) + lolModel('sigma')*dx(j)/2)) then
      phi_aux = 1.0;
    elseif ((lolModel('listofcutinf')(i,j) - k_overlap * lolModel('sigma')*dx(j)/2 <= x(j))&(x(j) <= xc(j) - lolModel('sigma')*dx(j)/2)) then
      phi_aux = (x(j) - (lolModel('listofcutinf')(i,j) - k_overlap * lolModel('sigma')*dx(j)/2)) / ...
                ((xc(j) - lolModel('sigma')*dx(j)/2) - (lolModel('listofcutinf')(i,j) - k_overlap * lolModel('sigma')*dx(j)/2));
    elseif ((xc(j) + lolModel('sigma')*dx(j)/2 <= x(j))&(x(j) <= lolModel('listofcutplus')(i,j) + k_overlap * lolModel('sigma')*dx(j)/2)) then
      phi_aux = (x(j) - (lolModel('listofcutplus')(i,j) + k_overlap * lolModel('sigma')*dx(j)/2)) / ...
                ((xc(j) + lolModel('sigma')*dx(j) / 2) - (lolModel('listofcutplus')(i,j) + k_overlap * lolModel('sigma')*dx(j)/2));
    else
      phi_aux = 0;
    end
    phi_tmp(i) = min(phi_tmp(i), phi_aux); // Certainly simpler most the analytical derivative ...
  end
end
vect_phi = phi_tmp(nopart) / max([sum(phi_tmp) %eps]);
endfunction

