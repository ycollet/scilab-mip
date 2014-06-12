function model = doe_poly_model(mod_type, nb_var, order)
// mod_type : 'lin' (linear model), 'poly' (pure polynomial model), 'inter' (linear model with interaction)
// nb_var   : number of variables
// order    : order of the polynomial
if (~isdef('mod_type','local')) then
  mod_type = 'lin';
end
if (~isdef('nb_var','local')) then
  error('number of variables must be defined');
end
if (~isdef('order','local')) then
  order = 1;
end

if (mod_type=='lin') then
  order = 1;
end

if (mod_type=='inter') then
  order = nb_var;
end

size_index = length(string(nb_var));

model = '1 ';

if (mod_type=='lin') then
  for i=1:nb_var
    model = model + 'x' + repchar(string(0),size_index - length(string(i))) + string(i) + ' ';
  end
else
  // monomial production
  n = (order+1)^nb_var;
  Count = zeros(1,nb_var);
  for i=1:n
    Count(1) = Count(1) + 1;
    if (Count(1)>order) then
      Carry    = 1;
      Count(1) = 0;
    else
      Carry = 0;
    end
    for j=2:nb_var
      Count(j) = Count(j) + Carry;
      if (Count(j)>order) then
        Carry    = 1;
        Count(j) = 0;
      else
        Carry = 0;
      end
    end
    
    if ((sum(Count)<=order)&~((max(Count)>1)&(mod_type=='inter'))) then
      motif = [];
      for j=1:nb_var
        Operator = [];
        if (j~=nb_var) then
          if (sum(Count(j+1:$))==0) then
            Operator = ' ';
          else
            Operator = '*';
          end
        else
          Operator = ' ';
        end

        if (Count(j)~=0) then
          if (Count(j)==1) then
            motif = motif + 'x' + repchar(string(0),size_index - length(string(j))) + string(j) + Operator;
          end
          if ((Count(j)>1)&(mod_type=='poly')) then
            motif = motif + 'x' + repchar(string(0),size_index - length(string(j))) + string(j) + '^' + string(Count(j)) + Operator;
          end
        end // if
      end // for
      model = model + motif;
    end // if
  end // for 
end // else
endfunction

