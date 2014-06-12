function R = AddToVect(U,i,v)
// R = AddToVect(U,i,V)
// U : input vector
// R : resulting vector
// i : inserting position
// v : value to add
//
// A. Seghir, le 01/08/04

n = length(U);
R([1:i-1]) = U([1:i-1]);
R(i) = v;
R(i+1:n+1) = U([i:n]);
endfunction
