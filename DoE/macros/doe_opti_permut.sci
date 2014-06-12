function [M_doe, history] = doe_opti_permut(M_init, M_cand, doe_size, model, criterion, l_bounds, u_bounds, ItMX, p_level, Log)

[nargout, nargin] = argn();

if (size(M_cand,1)<=doe_size) then
  error('doe_opti_permut: candidate set too small');
end

if (isdef('M_init','local')) then
  M_doe = M_init;
  M_init_defined = %T;
else
  M_doe = [];
  M_init_defined = %F;
end

if (~isdef('model','local')) then
  error('doe_opti_permut: error, model is mandatory');
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

if (~isdef('criterion','local')) then
  // The criterion is to be maximized
  error('doe_opti_permut: error, criterion is mandatory');
end

if (nargout==2) then
  history = [];
  history_defined = %T;
else
  history_defined = %F;
end

for i=1:doe_size-size(M_doe,1)
  M_augm(i,:) = rand(1,size(M_cand,2)).*(u_bounds - l_bounds) + l_bounds;
end

M_doe = [M_doe; M_augm];

Progress = 1.0;

// Computation of the initial value of the d-optimality criterion
Value       = criterion(M_doe,model);
MaxValue    = Value;
OldValue    = Value;
Permutation = -1;

if (history_defined) then
  history = [history Value];
end

Iter = 0;

while((Progress>p_level)&(Iter<ItMX))
  Iter = Iter + 1;
  
  for i=max([size(M_init,1),1]):doe_size
    for j=1:size(M_cand,1)
      // Swap ith point from doe set and jth point from candidate set
      Aux         = M_doe(i,:);
      M_doe(i,:)  = M_cand(j,:);
      M_cand(j,:) = Aux;
      // Compute criterion
      Value = criterion(M_doe, model);
      // Record a good candidate
      if (Value>MaxValue) then
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
      // Swap optimal candidate from candidate set
      Aux                   = M_doe(i,:);
      M_doe(i,:)            = M_cand(Permutation,:);
      M_cand(Permutation,:) = Aux;
      Permutation           = -1;
      Value                 = MaxValue;

      if (history_defined) then
        history = [history Value];
      end
    end
  end

  Progress = (MaxValue - OldValue)/max([OldValue %eps]);
  OldValue = MaxValue;

  if (Log) then
    printf('doe_opti_permut: Iteration %d / %d - Progress %3.2f\n', Iter, ItMX, Progress);
  end
end
endfunction
