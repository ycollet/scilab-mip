function [M_doe, history] = doe_d_opti(M_init, M_cand, doe_size, model, l_bounds, u_bounds, ItMX, p_level, Log, size_tabu_list)

[nargout, nargin] = argn();

if (size(M_cand,1)<=doe_size) then
  error("doe_d_opti: candidate set too small");
end

if (isdef('M_init','local')) then
  M_doe = M_init;
  M_init_defined = %T;
else
  M_doe = [];
  M_init_defined = %F;
end

if (~isdef('model','local')) then
  error("doe_d_opti: error, model is mandatory");
end

if (~isdef('l_bounds','local')) then
  l_bounds = zeros(1,size(M_cand,2));
end

if (~isdef('u_bounds','local')) then
  u_bounds = ones(1,size(M_cand,2));
end

if (~isdef('ItMX','local')) then
  ItMX = 10;
end

if (~isdef('p_level','local')) then
  p_level = 0.05;
end

if (~isdef('Log','local')) then
  Log = %F;
end

if (~isdef('size_tabu_list')) then
  size_tabu_list = 6;
end

if (nargout==2) then
  history = [];
  history_defined = %T;
else
  history_defined = %F;
end

M_augm = [];
Tabu_List = zeros(1,size_tabu_list);

for i=1:doe_size-size(M_doe,1)
  M_augm(i,:) = M_cand(i,:);
end
for i=1:size(M_augm,1)
  Index1 = ceil(size(M_augm,1)*rand(1,1));
  Index2 = ceil(size(M_cand,1)*rand(1,1));
  M_augm(Index1,:) = M_cand(Index2,:);
end

M_doe = [M_doe; M_augm];

Progress = 1.0;

// Computation of the initial value of the d-optimality criterion
Value       = comp_d_opti_crit(M_doe,model);
MaxValue    = Value;
OldValue    = Value;
Permutation = -1;

if (history_defined) then
  history = [history Value];
end

Iter = 0;

while((Progress>=p_level)&(Iter<ItMX))
  Iter = Iter + 1;
  
  Tabu_List = zeros(1,size_tabu_list);
  
  for i=max([size(M_init,1),1]):doe_size
    for j=1:size(M_cand,1)
      if (~isempty(find(Tabu_List==j))) then
        if (Log) then
          printf('doe_d_opti: Tabu list used\n');
        end
        continue;
      end
      // Swap ith point from doe set and jth point from candidate set
      Aux         = M_doe(i,:);
      M_doe(i,:)  = M_cand(j,:);
      M_cand(j,:) = Aux;

      // Test the conditionnement of the matrix
      H = build_regression_matrix(M_doe,model);
      if (abs(cond(H'*H))<1e-4) then
        // Undo permutation
        Aux         = M_doe(i,:);
        M_doe(i,:)  = M_cand(j,:);
        M_cand(j,:) = Aux;
        if (Log) then
          printf('doe_d_opti: singular permutation skipped\n');
        end
        continue;
      end
      
      // Compute D-optimality criterion
      Value = comp_d_opti_crit(M_doe,model);
      // Record a good candidate
      if (Value>MaxValue) then
        if (Log) then
          printf('doe_d_opti: D crit = %f / MaxValue = %f / Permutation = %d / Point = %d / %d\n', Value, MaxValue,j, i, doe_size);
        end
        MaxValue    = Value;
        Permutation = j;
      end
      // Undo permutation
      Aux         = M_doe(i,:);
      M_doe(i,:)  = M_cand(j,:);
      M_cand(j,:) = Aux;
    end
    
    // Redo the best permutation found so far
    if (Permutation~=-1) then
      Tabu_List = [Tabu_List(1:$-1), Permutation];

      if (Log) then
        printf('doe_d_opti: Good permutation found between points %d and %d\n', i, Permutation);
      end

      // Swap optimal candidate from candidate set
      Aux                   = M_doe(i,:);
      M_doe(i,:)            = M_cand(Permutation,:);
      M_cand(Permutation,:) = Aux;
      Permutation           = -1;
      Value                 = MaxValue;

      if (history_defined) then
        history = [history Value];
      end
    else
      if (Log) then
        printf('doe_d_opti: No good permutation found for point %d\n', i);
      end
    end
  end

  Progress = abs(MaxValue - OldValue)/max([abs(OldValue) %eps]);

  if (Log) then
    printf('doe_d_opti: Iteration %d / %d - progress = %f / MaxValue = %f / OldValue = %f\n', Iter, ItMX, Progress, MaxValue, OldValue);
  end

  OldValue = MaxValue;
end
endfunction
