function r = doe_halton(dim_num, n, step, seed, leap, base)
// I_TO_HALTON_SEQUENCE: next N elements of an DIM_NUM-dimensional Halton sequence.
//
//  Discussion:
//
//    The DIM_NUM-dimensional Halton sequence is really DIM_NUM separate
//    sequences, each generated by a particular base.
//
//    This routine selects elements of a "leaped" subsequence of the
//    Halton sequence.  The subsequence elements are indexed by a
//    quantity called STEP, which starts at 0.  The STEP-th subsequence
//    element is simply element
//
//      SEED(1:DIM_NUM) + STEP * LEAP(1:DIM_NUM)
//
//    of the original Halton sequence.
//
//  Modified:
//
//    21 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    J H Halton,
//    On the efficiency of certain quasi-random sequences of points
//    in evaluating multi-dimensional integrals,
//    Numerische Mathematik,
//    Volume 2, 1960, pages 84-90.
// 
//    J H Halton and G B Smith,
//    Algorithm 247: Radical-Inverse Quasi-Random Point Sequence,
//    Communications of the ACM,
//    Volume 7, 1964, pages 701-702.
//
//    Ladislav Kocis and William Whiten,
//    Computational Investigations of Low-Discrepancy Sequences,
//    ACM Transactions on Mathematical Software,
//    Volume 23, Number 2, 1997, pages 266-294.
//
//  Parameters:
//
//    Input, integer DIM_NUM, the spatial dimension.
//    1 <= DIM_NUM is required.
//
//    Input, integer N, the number of elements of the sequence.
//
//    Input, integer STEP, the index of the subsequence element.
//    0 <= STEP is required.
//
//    Input, integer SEED(DIM_NUM), the Halton sequence index corresponding
//    to STEP = 0.
//
//    Input, integer LEAP(DIM_NUM), the succesive jumps in the Halton sequence.
//
//    Input, integer BASE(DIM_NUM), the Halton bases.
//
//    Output, double precision R(DIM_NUM,N), the next N elements of the
//    leaped Halton subsequence, beginning with element STEP.
//

//
//  Check the input.
//

if (~isdef('dim_num','local')) then
  error('doe_halton: dim_num is mandatory');
end

if (~isdef('n','local')) then
  error('doe_halton: n is mandatory');
end

if (~isdef('step','local')) then
  step = 0;
end

if (~isdef('seed','local')) then
  seed(1:dim_num) = 0;
end

if (~isdef('leap','local')) then
  leap(1:dim_num) = 1;
end

if (~isdef('base','local')) then
  for i=1:dim_num
    base(i) = doe_prime(i);
  end
end

//
// Setting internal parameters
//

dim_num         = floor(dim_num);
n               = floor(n);
step            = floor(step);
seed(1:dim_num) = floor(seed(1:dim_num));
leap(1:dim_num) = floor(leap(1:dim_num));
base(1:dim_num) = floor(base(1:dim_num));

//
//  Check the input.
//

if (dim_num<1) then
  error('doe_halton: error, dim_num must be greater than 1');
end

if (n<1) then
  error('doe_halton: error, n must be greater than 1');
end

if (step<0) then
  error('doe_halton: error, step must be greater than 0');
end

if (or(seed(1:dim_num)<0)) then
  error('doe_halton: error, seed must be greater than 0');
end

if (or(leap(1:dim_num)<0)) then
  error('doe_halton: error, leap must be greater than 1');
end

if (or(base(1:dim_num)==0)|or(base(1:dim_num)==1)) then
  error('doe_halton: error, base must be equal to 0 or to 1 for some base(I)');
end

//
//  Calculate the data.
//

r(1:dim_num,1:n) = 0.0;

for i=1:dim_num
  seed2(1:n) = seed(i) + step * leap(i) : leap(i) : seed(i) + (step + n - 1) * leap(i);

  base_inv = 1.0 / base(i);

  while (or(seed2 ~= 0 ))
    digit(1:n) = modulo(seed2(1:n), base(i));
    r(i,1:n)   = r(i,1:n) + digit(1:n) * base_inv;
    base_inv   = base_inv / base(i);
    seed2(1:n) = floor(seed2(1:n) / base(i));
  end
end

r = r';

for i=1:size(r,2)
  Min = min(r(:,i));
  Max = max(r(:,i));
  for j=1:size(r,1)
    r(j,i) = (r(j,i) - Min) / (Max - Min);
  end
end
endfunction
