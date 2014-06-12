function [model_new, coeff_new] = doe_model_bselect(nb_var, model_old, measures, Log)
// Last column of measures corresponds to the output. The others columns are inputs measurement
// list_var : a string which contains all the var
if (~isdef('Log','local')) then
  Log = %F;
end

n            = 2^nb_var-1;
Count        = ones(1,nb_var);
Max_Residual = %inf;
Max_Pos      = -1;
ListOfTokens = tokens(model_old,' ');
build        = ones(1,size(ListOfTokens,1));
size_index   = length(string(nb_var));

for l=1:nb_var
  for i=1:nb_var
    if (Count(i)==0) then
      continue;
    end
    
    Count(i) = 0;
    // Preparation of the data set
    build(:) = %T;
    for j=1:nb_var
      if (Count(j)==1) then
        Monom_Index = grep(ListOfTokens,'x' + repchar(string(0),size_index - length(string(j))) + string(j));
        disp(Monom_Index);
        build(Monom_Index) = %F;
      end
    end
    // Construction of the regression matrix
    R = build_regression_matrix(measures(:,1:$-1), model_old, build);
    // Compute the model coefficient
    Coeff_model = R \ measures(:,$);
    // Compute the residual
    Residual = 0;
    for j=1:size(measures,1)
      Residual = Residual + (Coeff_model'*R(j,:)' - measures(j,$))^2;
    end
    if (Residual<Max_Residual) then
      Max_Residual = Residual;
      Max_Count    = Count;
      Max_Pos      = i;
      coeff_new    = Coeff_model;
    end
    if (Log) then
      printf('doe_model_bselect: Residual = %f / %f ', Residual, Max_Residual);
      printf('- Pass %d / %d ', l, nb_var);
      printf('- Variable %d / %d \n', i, nb_var);
    end
        
    Count(i) = 1;
  end

  if ((Log)&(Max_Pos~=-1)) then
    printf('doe_model_bselect: Variable %d removed. Residual = %f\n', Max_Pos, Max_Residual);
  end

  if (Max_Pos==-1) then
    break;
  else
    Count(Max_Pos) = 0;
    Max_Pos = -1;
  end
end

// Build a new model
for i=1:nb_var
  if (Max_Count(i)==0) then
    list_index_op = strindex(model_old,'*');
    list_index_op = [list_index_op strindex(model_old,'/')];
    list_index_op = [list_index_op strindex(model_old,'+')];
    list_index_op = [list_index_op strindex(model_old,'-')];
    list_index_op = [list_index_op strindex(model_old,' ')];
    list_index_op = [list_index_op size(list_index_op,1)];

    list_index_op = gsort(list_index_op)
    list_index_x  = strindex(model_old,'x' + repchar(string(0),size_index - length(string(i))) + string(i));
    for j=1:size(list_index_x,1)
      Index = 1;
      while (list_index_x(j)<list_index_op(Index)) Index = Index + 1; end
      model_new = part(model_old,1:list_index_op(Index)) + part(model_old,list_index_op(Index-1):list_index_op(1));
    end
  end
end
endfunction
