function [H] = hadamard(n)
// This function computes a Hadamard matrix of size 2^n
if (n<1) then
  error('hadamard: n must be a strictly positive value');
end

// Initial Hadamard matrix (H2)

H = [1, 1; -1, 1];

for i=1:n-1
  H = [H,H;-H,H];
end
endfunction
