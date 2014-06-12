function R = build_regression_matrix(H,model,build)
// build_regression_matrix : a function which builds a regression matrix using a DoE matrix and a list of monomials
// H      : a doe matrix or a vector
// model  : a list of monomials like '1 x1 x2, x1*x2' where xi corresponds to the ith column
// R      : the regression matrix

ListOfTokens = tokens(model, ' ');
if (size(H,2)==1) then // Vector mode
  size_index = length(string(size(H,1)));
else
  size_index = length(string(size(H,2)));
end

if (~isdef('build','local')) then
  build = %T*ones(1,size(ListOfTokens,1));
end

// Test: if the vector has wrong direction (lines instead of columns)
if (size(H,1)==1) then
  H = H';
end

// We build the list of monomials to add
Monom_Index = [];
for i=1:size(ListOfTokens,1)
  if (build(i)) then
    Monom_Index = [Monom_Index grep(ListOfTokens, 'x' + repchar(string(0),size_index - length(string(i))) + string(i))];
  end
end

Monom_Index = -gsort(-Monom_Index);
Monom_Index = unique(Monom_Index);

if (size(H,2)==1) then
  nb_var     = size(H,1); // Vector "mode": The number of vars is equal to the number of lines of H
  VectorMode = %T;

  for i=1:nb_var
    for j=1:size(ListOfTokens,1)
      ListOfTokens(j) = strsubst(ListOfTokens(j), 'x' + repchar(string(0),size_index - length(string(i))) + string(i), 'H('+string(i)+')');
    end
  end
else
  nb_var     = size(H,1); // Matrix "mode": The number of vars is equal to the number of columns of H
  VectorMode = %F;

  for i=1:nb_var
    for j=1:size(ListOfTokens,1)
      ListOfTokens(j) = strsubst(ListOfTokens(j), 'x' + repchar(string(0),size_index - length(string(i))) + string(i), 'H(i,'+string(i)+')');
    end
  end
end

// We build the regression matrix
if (VectorMode) then
//  R = zeros(size(Monom_Index,2),1);
  R = zeros(size(ListOfTokens,1),1);

//  for j=1:size(Monom_Index,2)
//    R(j,1) = eval(ListOfTokens(Monom_Index(j)));
  for j=1:size(ListOfTokens,1)
    R(j,1) = eval(ListOfTokens(j));
  end
else
  //R = zeros(size(H,1), size(Monom_Index,2));
  R = zeros(size(H,1), size(ListOfTokens,1));

  for i=1:size(H,1)
    //for j=1:size(Monom_Index,2)
    for j=1:size(ListOfTokens,1)
      //R(i,j) = eval(ListOfTokens(Monom_Index(j)));
      R(i,j) = eval(ListOfTokens(j));
    end
  end
end

return R;
endfunction
