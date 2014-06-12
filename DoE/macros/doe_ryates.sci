function [y, id] = doe_ryates(ef)
// RYATES Reverse Yates' algorithm to give estimated responses
//
// CALL:  y = ryates(ef);
//
//  y  = Estimated response given the effects.
//  ef = vector of average response, main effects and interaction effects.
//
// RYATES applies the reverse Yates' algorithm to the effect EF to obtain 
// the estimated response. EF is assumed to
// be arranged in what is called standard order. (The order of the actual
// running should, of course, be random).  EF(1,:) is the
// average response and EF(2:end,:) contain the main effects and
// interaction effects.
//
// Example:
//   D = ffd(3);                    // complete 2^3 design in standard order.
//   y = [60 72 54 68 52 83 45 80]; // Responses to design D.
//   [ef,id] = yates(y);
//   y1 = ryates(ef);               // gives the same as Y
//
// See also  ffd



// Reference 
// Box, G.E.P, Hunter, W.G. and Hunter, J.S. (1978)
// Statistics for experimenters, John Wiley & Sons, pp 342

// Tested on: Matlab 5.3
// History:
// By Per A. Brodtkorb 16.03.2001

if (size(ef,1)==1) then 
  ef = ef';       // Make sure it is a column vector
end
n = size(ef,1);   // Number of runs

k = log2(n);      // Number of variables.
if (round(k)~=k) then
  error('The length of EF must be in power of two');
end

// Reverse yates algorithm:
y      = ef*(n/2);
y(1,:) = y(1,:)*2;
if (nargout>1) then
  [y,id] = doe_yates(flipud(y));
else
  y = doe_yates(flipud(y));
end
y = flipud(y)/2;
y($,:) = y($,:)*2;
endfunction




