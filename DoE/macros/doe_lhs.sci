function lhs_matrix = doe_lhs(nb_dims, x_min, x_max, nb_div, nb_iterations, nb_points, random)

if (~isdef('nb_points','local')) then 
  nb_points = nb_dims;
end
if (~isdef('nb_iterations','local')) then
  nb_iterations = 3*nb_points;
end
if (~(isdef('x_min','local'))&isdef('x_max','local')) then
  error('x_min and x_max are mandatory parameters');
end
if (~isdef('nb_div','local')) then
  nb_div = nb_dims;
end
if (~isdef('random','local')) then
  random = %F;
end
sequence = 1:max(nb_div);
Aux = sequence;
for i=1:ceil(nb_points/nb_div)
  Aux = [Aux sequence];
end
sequence = Aux(1:nb_points);
lhs_matrix = zeros(nb_points, nb_dims);
for i=1:nb_dims
  lhs_matrix(:,i) = sequence';
end
for i=1:nb_iterations
  for j=1:nb_dims
    Index1 = ceil(rand(1,1)*nb_points);
    Index2 = ceil(rand(1,1)*nb_points);
    Aux                  = lhs_matrix(Index1,j);
    lhs_matrix(Index1,j) = lhs_matrix(Index2,j);
    lhs_matrix(Index2,j) = Aux;
  end
end
if (random) then
  for i=1:nb_points
    for j=1:nb_dims
      if random then offset = 1*rand(1,1) - 0.5;
      else           offset = 0; end
      lhs_matrix(i,j) = lhs_matrix(i,j) + offset - 1;
    end
  end
  
  for i=1:size(lhs_matrix,2)
    Min = min(lhs_matrix(:,i));
    Max = max(lhs_matrix(:,i));
    for j=1:size(lhs_matrix,1)
      lhs_matrix(j,i) = (lhs_matrix(j,i) - Min)/(Max - Min);
      lhs_matrix(j,i) = (x_max(i) - x_min(i))*lhs_matrix(j,i) + x_min(i);
    end
  end
end
endfunction
