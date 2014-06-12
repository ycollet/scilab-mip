function [quasi, seed] = doe_faure(dim_num, n, seed)

//// FAURE generates a new quasirandom Faure vector with each call.
//
//  Discussion:
//
//    This routine implements a method of H. Faure for computing
//    quasirandom numbers.  It is a merging and adaptation of
//    the routines INFAUR and GOFAUR from ACM TOMS 647.
//
//    Thanks to Ernst Kloppenburg for suggesting the use of persistent
//    variables to improve the MATLAB implementation.
//
//  Modified:
//
//    16 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    H Faure,
//    Discrepance de suites associees a un systeme de numeration
//    (en dimension s),
//    Acta Arithmetica,
//    Volume XLI, 1982, pages 337-351, especially page 342.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//  Parameters:
//
//    Input, integer DIM_NUM, the spatial dimension, which should be
//    at least 2.
//
//    Input, integer SEED, the seed, which indicates the index of the
//    element of the sequence to be calculated.  If SEED is negative,
//    it is effectively replaced by a more suitable value.
//
//    Output, real QUASI(DIM_NUM), the next quasirandom vector.
//
//    Output, integer SEED, the appropriate value of SEED to be
//    used on the next call, if the next element of the sequence is desired.
//

global hisum_save;
global coef;
global qs;

if ~isdef('seed','local') then
  seed = -1;
end

//
//  If the internal variable QS has never been set, then its SIZE is zero.
//
if (size(qs) == 0 | seed <= 0) then
  qs = prime_ge(dim_num);
  hisum_save = -1;
end

for ii=1:n
  //
  //  If SEED < 0, reset for recommended initial skip.
  //
  if (seed < 0) then
    hisum = 3;
    seed = qs^(hisum+1) - 1;
  elseif (seed == 0) then
    hisum = 0;
  else
    hisum = 0;
    i = seed;
    while (qs <= i)
      i = i / qs;
      hisum = hisum + 1;
    end
  end
  //
  //  Is it necessary to recompute the coefficient table?
  //
  if (hisum_save < hisum) then
    hisum_save = hisum;
    coef(1:hisum+1,1:hisum+1) = 0;
    coef(2:hisum+1,1) = 1;
    for j=1:hisum
      coef(j+1,j+1) = 1;
    end
    for j=1:hisum
      for i=j+1:hisum
        coef(i+1,j+1) = modulo(coef(i-1+1,j+1) + coef(i-1+1,j-1+1), qs);
      end
    end
  end
  //
  //  Find QUASI(1) using the method of Faure.
  //
  //  SEED has a representation in base QS of the form:
  //
  //    Sum ( 0 <= J <= HISUM ) YTEMP(J) * QS^J
  //
  //  We now compute the YTEMP(J)'s.
  //
  ktemp = qs^( hisum + 1 );
  ltemp = seed;
  for i=hisum:-1:0
    ktemp = ktemp / qs;
    mtemp = modulo(ltemp,ktemp);
    ytemp(i+1) = (ltemp - mtemp) / ktemp;
    ltemp = mtemp;
  end
  //
  //  QUASI(K) has the form
  //
  //    Sum ( 0 <= J <= HISUM ) YTEMP(J) / QS**(J+1)
  //
  //  Compute QUASI(1) using nested multiplication.
  //
  r = ytemp(hisum+1);
  for i=hisum-1:-1:0
    r = ytemp(i+1) + r / qs;
  end
  quasi(ii,1) = r / qs;
  //
  //  Find components QUASI(2:DIM_NUM) using the Faure method.
  //
  for k=2:dim_num
    quasi(ii,k) = 0.0;
    r = 1.0 / qs;
    for j=0:hisum
      ztemp = ytemp(j+1:hisum+1)' * coef(j+1:hisum+1,j+1);
  //
  //  New YTEMP(J) is:
  //
  //    Sum ( J <= I <= HISUM ) ( old ytemp(i) * binom(i,j) ) mod QS.
  //
      ytemp(j+1) = modulo(ztemp, qs);
      quasi(ii,k) = quasi(ii,k) + ytemp(j+1) * r;
      r = r / qs;
    end
  end
  //
  //  Update SEED.
  //
  seed = seed + 1;
end
endfunction

function p = prime_ge(n)
// PRIME_GE returns the smallest prime greater than or equal to N.
//
//  Discussion:
//
//    The MATLAB version of this program is made much simpler
//    because of the availability of the IS_PRIME logical function.
//
//  Examples:
//
//    N     PRIME_GE
//
//    -10    2
//      1    2
//      2    2
//      3    3
//      4    5
//      5    5
//      6    7
//      7    7
//      8   11
//      9   11
//     10   11
//
//  Modified:
//
//    15 March 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number to be bounded.
//
//    Output, integer P, the smallest prime number that is greater
//    than or equal to N.  
//

p = max(ceil(n),2);
  
while (~isprime(p))
  p = p + 1;
end
endfunction  
