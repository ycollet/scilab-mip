function  S = gridsamp(_range, q)
//GRIDSAMP  n-dimensional grid over given range
//
// Call:    S = gridsamp(range, q)
//
// range :  2*n matrix with lower and upper limits
// q     :  n-vector, q(j) is the number of points
//          in the j'th direction.
//          If q is a scalar, then all q(j) = q
// S     :  m*n array with points, m = prod(q)

// hbn@imm.dtu.dk  
// Last update June 25, 2002

[mr,n] = size(_range);    dr = diff(_range,1,'r');
if  mr ~= 2 | or(dr < 0) then
  error('range must be an array with two rows and range(1,:) <= range(2,:)')
end 
sq = size(q);
if  min(sq) > 1 | or(q <= 0) then
  error('q must be a vector with non-negative elements')
end
p = length(q);   
if  p == 1 then q = ones(1,n) .*. q; 
elseif p ~= n then
  error(sprintf('length of q must be either 1 or %d',n))
end 

// Check for degenerate intervals
i = find(dr == 0);
if  ~isempty(i) then q(i) = 0*q(i); end

// Recursive computation
if n > 1 then
  A = gridsamp(_range(:,2:$), q(2:$));  // Recursive call
  [m p] = size(A);   q = q(1);
  S = [zeros(m*q,1) ones(q,1) .*. A];
  y = linspace(_range(1,1),_range(2,1), q);
  k = 1:m;
  for  i = 1 : q
    S(k,1) = ones(m,1) .*. y(i);  k = k + m;
  end
else    
  S = linspace(_range(1,1),_range(2,1), q).';
end
endfunction

