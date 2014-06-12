//////////////////////////////////////////
// Definition of the objective function //
//////////////////////////////////////////

function [f]=calfun(n, x)
// Example of a paraboloide of dimension 4 for which the minimum is 1 1 1 1
f = sum((x-1).^2);
endfunction

//////////////////
// Main program //
//////////////////

n      = 4; // Paraboloide dimension 4
iprint = 2;
maxfun = 5000;
rhoend = 1e-6;
x = 10*ones(n,1); // Initialization far from the minimum

npt = 2*n + 1;
rhobeg = x(1) * 0.01;
	
xopt = newuoa(npt, x, rhobeg, rhoend, iprint, maxfun, calfun);
	
