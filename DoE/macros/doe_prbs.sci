
// For a given size of binary code, we give the recommended feedbacks so as
// to produce the maximum sequence of random binary codes
// n  length     feedback
// 2  3          [0,1]
// 3  7          [0,2]
// 4  15         [0,3]
// 5  31         [1,4]
// 6  63         [0,5]
// 7  127        [0,6]
// 8  255        [1,2,3,7]
// 9  511        [3,8]
// 10 1023       [2,9]
// 11 2047       [1,10]
// 12 4095       [0,3,5,11]
// 13 8191       [0,2,3,12]
// 14 16383      [0,2,4,13]
// 15 32767      [0,14]
// 16 65535      [1,2,4,15]
// 17 131071     [2,16]
// 18 262143     [6,17]
// 19 524287     [0,1,4,18]
// 20 1048575    [2,19]
// 21 2097151    [1,20]
// 22 4194303    [0,21]
// 23 8388607    [4,22]
// 24 16777215   [0,2,3,23]
// 25 33554431   [2,24]
// 26 67108863   [0,1,5,25]
// 27 134217727  [0,1,4,26]
// 28 268435455  [2,27]
// 29 536870911  [1,28]
// 30 1073741823 [0.3,5,29]
// 31 2147483647 [2,30]
// 32 4294967295 [1,5,6,31]

function Result = doe_prbs(init, feedback)
if (type(init)==1) then
  init = (init~=0);
end
Aux = xor(init(feedback));
Result = zeros(size(init,1), size(init,2));
Result(2:$) = init(1:$-1);
Result(1)   = Aux;
endfunction

function xRes = xor(A,B)
if (size(A,1)==1) then A = A'; end;
if (isdef('B','local')) then
  if (size(B,1)==1) then B = B'; end;
  Both = %T;
else
  Both = %F;
end

if (Both) then
  xRes = (sum([A B], 'c')==ones(size(A,1),1));
else
  xRes = (sum(A)==1);
end
endfunction
