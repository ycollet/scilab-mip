// Test problems extracted from the following document:
//
// Deterministic Methods for Mixed Integer Nonlinear Programming
// Sven Leyffer
// PhD Thesis
// Department of Mathematics & Computer Science
// University of Dundee
// Dundee
// December 1993

///////////////////////
// TP1 MINLP problem //
///////////////////////

// The integer variables are x(4) x(5) x(6) in {0,1}
// f(x_sol) = 6.010
// f(x_sol_relax) = 0.759

function y = MINLPTP1_mip_var()
y = [4 5 6];
endfunction

function y = MINLPTP1_x_init()
y = ones(6,1);
endfunction

function y = MINLPTP1_x_sol()
y = [1.301 0 1.0 1 0 1]';
endfunction

function y = MINLPTP1_x_sol_relax()
y = [1.147 0.547 1.0 0.273 0.3 0.0]';
endfunction

function Res = min_bd_MINLPTP1()
Res = [0 0 0 0 0 0]';
endfunction

function Res = max_bd_MINLPTP1()
Res = [2 2 1 1 1 1]';
endfunction

function y = MINLPTP1obj(x)
y = 5*x(3+1) + 6*x(3+2) + 8*x(3+3) + 10*x(1) - 7*x(3) - 18*log(x(2) + 1) - 19.2*log(x(1) - x(2) + 1) + 10;
endfunction

function y = MINLPTP1ineq(x)
y(1) = -(0.8*log(x(2)+1) + 0.96*log(x(1) - x(2) + 1) - 0.8*x(3));
y(2) = -(log(x(2) + 1) + 1.2*log(x(1) - x(2) + 1) - x(3) -2*x(3+3) + 2);
y(3) = x(2) - x(1);
y(4) = x(2) - 2*x(3+1);
y(5) = x(1) - x(2) - 2*x(3+2);
y(6) = x(3+1) + x(3+2) - 1;
endfunction

MINLPTP1eq = [];

///////////////////////
// TP2 MINLP problem //
///////////////////////

// The integer variables are x(7) x(8) x(9) x(10) x(11) in {0,1}
// f(x_sol) = 73.035
// f(x_sol_relax) = -0.554

function y = MINLPTP2_mip_var()
y = [7 8 9 10 11];
endfunction

function y = MINLPTP2_x_init()
y = ones(11,1);
endfunction

function y = MINLPTP2_x_sol()
y = [0 2 1.078 0.652 0.326 1.078 0 1 1 1 0]';
endfunction

function y = MINLPTP2_x_sol_relax()
y = [1.903 2 2 1.403 0.701 2 0.571 0.429 0.25 0.21 0]';
endfunction

function Res = min_bd_MINLPTP2()
Res = [0 0 0 0 0 0 0 0 0 0 0]';
endfunction

function Res = max_bd_MINLPTP2()
Res = [2 2 2 %inf %inf 3 1 1 1 1 1]';
endfunction

function y = MINLPTP2obj(x)
y = 5*x(6+1) + 8*x(6+2) + 6*x(6+3) + 10*x(6+4) + 6*x(6+5) - 10*x(1) - 15*x(2) - 15*x(3) + 15*x(4) + 5*x(5) - 20*x(6) + ...
    exp(x(1)) + exp(x(2)/1.2) - 60*log(x(4) + x(5) + 1) + 140;
endfunction

function y = MINLPTP2ineq(x)
y(1)  = -log(x(4) + x(5) + 1);
y(2)  = exp(x(1)) - 10*x(6+1) - 1;
y(3)  = exp(x(2)/1.2) - 10*x(6+2) - 1;
y(4)  = 1.25*x(3) - 10*x(6+3);
y(5)  = x(4) + x(5) - 10*x(6+4);
y(6)  = -2*x(3) + 2*x(6) - 10*x(6+5);
y(7)  = -x(1) - x(2) - 2*x(3) + x(4) +2*x(6);
y(8)  = -x(1) - x(2) - 0.75*x(3) + x(4) + 2*x(6);
y(9)  = x(3) - x(6);
y(10) = 2*x(3) - x(4) - 2*x(6);
y(11) = - 0.5*x(4) + x(5);
y(12) = 0.2*x(4) - x(5);
y(13) = x(6+4) + x(6+5) - 1;
endfunction

function y = MINLPTP2eq(x)
y(1) = x(7) + x(8) - 1;
endfunction

///////////////////////
// TP3 MINLP problem //
///////////////////////

// The integer variables are x(10) x(11) x(12) x(13) x(14) x(15) x(16) x(17) in {0,1}
// f(x_sol) = 6.010
// f(x_sol_relax) = 0.759

function y = MINLPTP3_mip_var()
y = [10 11 12 13 14 15 16 17];
endfunction

function y = MINLPTP3_x_init()
y = ones(17,1);
endfunction

function y = MINLPTP3_x_sol()
y = [0 2 0.468 0.585 2 0 0 0.267 0.585 0 1 0 1 0 1 0 1]';
endfunction

function y = MINLPTP3_x_sol_relax()
y = [1.903 2 0.528 0.659 2 1.083 0.659 0.411 0 0.571 0.430 0.066 0.308 0 0.2 0.108 0.119]';
endfunction

function Res = min_bd_MINLPTP3()
Res = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
endfunction

function Res = max_bd_MINLPTP3()
Res = [2 2 1 2 2 2 2 1 3 1 1 1 1 1 1 1 1]';
endfunction

function y = MINLPTP3obj(x)
y = 5*x(9+1) + 8*x(9+2) + 6*x(9+3) + 10*x(9+4) + 6*x(9+5) + 7*x(9+6) + 4*x(9+7) + 5*x(9+8) - ...
    10*x(1) - 15*x(2) + 15*x(3) + 80*x(4) + 25*x(5) + 35*x(6) - 40*x(7) + ...
    15*x(8) - 35*x(9) + exp(x(1)) + exp(x(2)/1.2) - 65*log(x(3) + x(4) + 1) - ...
    90*log(x(5) + 1) - 80*log(x(6) + 1) + 120;
endfunction

function y = MINLPTP3ineq(x)
y(1)  = -1.5*log(x(5) + 1) - log(x(6) + 1) - x(8);
y(2)  = -log(x(3) + x(4) + 1);
y(3)  = -x(1) - x(2) + x(3) + 2*x(4) + 0.8*x(5) + 0.8*x(6) - 0.5*x(7) - x(8) - 2*x(9);
y(4)  = -x(1) - x(2) + 2*x(4) + 0.8*x(5) + 0.8*x(6) - 2*x(7) - x(8) - 2*x(9);
y(5)  = -2*x(4) - 0.8*x(5) - 0.8*x(6) + 2*x(7) + x(8) + 2*x(9);
y(6)  = -0.8*x(5) - 0.8*x(6) + x(8);
y(7)  = -x(4) + x(7) + x(9);
y(8)  = -0.4*x(5) - 0.4*x(6) + 1.5*x(8);
y(9)  = 0.16*x(5) + 0.16*x(6) - 1.2*x(8);
y(10) =  x(3) - 0.8*x(4);
y(11) = -x(3) + 0.4*x(4);
y(12) = exp(x(1)) - 10*x(9+1) - 1;
y(13) = exp(x(2)/1.2) - 10*x(9+2) - 1;
y(14) = x(7) - 10*x(9+3);
y(15) = 0.8*x(5) + 0.8*x(6) - 10*x(9+4);
y(16) = 2*x(4) - 2*x(7) - 2*x(9) - 10*x(9+5);
y(17) = x(5) - 10*x(9+6);
y(18) = x(6) - 10*x(9+7);
y(19) = x(3) + x(4) - 10*x(9+8);
y(20) = x(9+4) + x(9+5) - 1;
y(21) = x(9+3) - x(9+8);
endfunction

function y = MINLPTP3eq(x)
y(1) = x(9+1) + x(9+2) - 1;
y(2) = - x(9+4) + x(9+6) + x(9+7);
endfunction

///////////////////////////
// ASAADI1 MINLP problem //
///////////////////////////

// The integer variables are x(1) x(2) x(3) x(4) in {0,1}
// f(x_sol) = -40.957
// f(x_sol_relax) = -40.963

function y = MINLPASAADI1_mip_var()
y = [1 2 4];
endfunction

function y = MINLPASAADI1_x_init()
y = ones(4,1);
endfunction

function y = MINLPASAADI1_x_sol()
y = [0 1 2.236 0]';
endfunction

function y = MINLPASAADI1_x_sol_relax()
y = [0 1.038 2.227 0.000]';
endfunction

function Res = min_bd_MINLPASAADI1()
Res = -100*ones(4,1);
endfunction

function Res = max_bd_MINLPASAADI1()
Res = 100*ones(4,1);
endfunction

function y = MINLPASAADI1obj(x)
y = x(1)^2 + x(2)^2 + 2*x(3)^2 + x(4)^2 - 5*x(1) - 5*x(2) - 21*x(3) + 7*x(4);
endfunction

function y = MINLPASAADI1ineq(x)
y(1) = -(-x(1)^2 - x(2)^2 - x(3)^2 - x(4)^2 - x(1) + x(2) - x(3) + x(4) + 8);
y(2) = -(-x(1)^2 - 2*x(2)^2 - x(3)^2 - 2*x(4)^2 + x(1) + x(4) + 10);
y(3) = -(-2*x(1)^2 - x(2)^2 - x(3)^2 - 2*x(1) + x(2) + x(4) + 5);
endfunction

MINLPASAADI1eq = [];

///////////////////////////
// ASAADI2 MINLP problem //
///////////////////////////

// The integer variables are x(1) to x(7)
// f(x_sol) = 694.90
// f(x_sol_relax) = 683.981

function y = MINLPASAADI2_mip_var()
y = [1 2 3 4];
endfunction

function y = MINLPASAADI2_x_init()
y = ones(7,1);
endfunction

function y = MINLPASAADI2_x_sol()
y = [2 2 0 4 -4.63e-4 1.132 1.463]';
endfunction

function y = MINLPASAADI2_x_sol_relax()
y = [2.348 1.935 0 4.298 0 1.048 1.582]';
endfunction

function Res = min_bd_MINLPASAADI2()
Res = [-1 -1 -1 -1 -1 -1 -1]';
endfunction

function Res = max_bd_MINLPASAADI2()
Res = [5 5 5 5 5 5 5]';
endfunction

function y = MINLPASAADI2obj(x)
y = (x(1) - 10)^2 + 5*(x(2) - 12)^2 + x(3)^4 + 3*(x(4) - 11)^2 + 10*x(5)^6 + 7*x(6)^2 + x(7)^4 - 4*x(6)*x(7) - 10*x(6) - 8*x(7);
endfunction

function y = MINLPASAADI2ineq(x)
y(1) = -(-2*x(1)^2 - 3*x(2)^4 - x(3) - 4*x(4)^2 - 5*x(5) + 127);
y(2) = -(7*x(1) - 3*x(2) - 10*x(3)^2 - x(4) + x(5) + 282);
y(3) = -(23*x(1) - x(2)^2 - 6*x(6)^2 + 8*x(7) + 196);
y(4) = -(-4*x(1)^2 - x(2)^2 + 3*x(1)*x(2) - 2*x(3)^2 - 5*x(6) + 11*x(7));
endfunction

MINLPASAADI2eq = [];

///////////////////////////
// ASAADI3 MINLP problem //
///////////////////////////

// The integer variables are x(1) to x(10)
// f(x_sol) = 37.219
// f(x_sol_relax) = 24.306

function y = MINLPASAADI3_mip_var()
y = [1 3 5 7 8 9];
endfunction

function y = MINLPASAADI3_x_init()
y = ones(10,1);
endfunction

function y = MINLPASAADI3_x_sol()
y = [2 2.600 8 5.000 1 .387 2 10 8 8.6]';
endfunction

function y = MINLPASAADI3_x_sol_relax()
y = [2.172 2.364 8.774 5.096 0.991 1.431 1.322 9.829 8.280 8.376]';
endfunction

function Res = min_bd_MINLPASAADI3()
Res = -100*ones(10,1);
endfunction

function Res = max_bd_MINLPASAADI3()
Res = 100*ones(10,1);
endfunction

function y = MINLPASAADI3obj(x)
y = x(1)^2 + x(2)^2 + x(1)*x(2) - 14*x(1) - 16*x(2) + (x(3) - 10)^2 + ...
    4*(x(4) - 5)^2 + (x(5) - 3)^2 + 2*(x(6) - 1)^2 + 5*x(7)^2 + ...
    7*(x(8) - 11)^2 + 2*(x(9) - 10)^2 + (x(10) - 7)^2 + 45;
endfunction

function y = MINLPASAADI3ineq(x)
y(1) = -(-3*(x(1) - 2)^2 - 4*(x(2) - 3)^2 - 2*x(3)^2 + 7*x(4) + 120);
y(2) = -(-5*x(1)^2 - 8*x(2) - (x(3) - 6)^2 + 2*x(4) + 40);
y(3) = -(-1/2*(x(1) - 8)^2 - 2*(x(2) - 4)^2 - 3*x(5)^2 + x(6) + 30);
y(4) = -(-x(1)^2 - 2*(x(2) - 2)^2 + 2*x(1)*x(2) - 14*x(5) + 6*x(6));
y(5) = -(3*x(1) - 6*x(2) - 12*(x(9) - 8)^2 + 7*x(10));
y(6) = -(-4*x(1) - 5*x(2) + 3*x(7) - 9*x(8) + 105);
y(7) = -(-10*x(1) + 8*x(2) + 17*x(7) - 2*x(8));
y(8) = -(8*x(1) - 2*x(2) - 5*x(9) + 2*x(10) + 12);
endfunction

MINLPASAADI3eq = [];

////////////////////////
// 2DEx MINLP problem //
////////////////////////

// The integer variables are x(1) x(2) integers
// f(x_sol) = -56.938
// f(x_sol_relax) = -56.944

function y = MINLP2DEX_mip_var()
y = [1 2];
endfunction

function y = MINLP2DEX_x_init()
y = ones(2,1);
endfunction

function y = MINLP2DEX_x_sol()
y = 4*[4 4.75]';
endfunction

function y = MINLP2DEX_x_sol_relax()
y = 4*[3.947 4.776]';
endfunction

function Res = min_bd_MINLP2DEX()
Res = 4*[0 3]';
endfunction

function Res = max_bd_MINLP2DEX()
Res = 4*[5 5]';
endfunction

function y = MINLP2DEXobj(x)
x = x / 4;
y = 2*x(1)^2 + x(2)^2 - 16*x(1) - 10*x(2);
endfunction

function y = MINLP2DEXineq(x)
x = x / 4;
y(1) = x(1)^2 - 6*x(1) + x(2) - 11;
y(2) = -x(1)*x(2) + 3*x(2)^2 + exp(x(1) - 3) - 1;
endfunction

MINLP2DEXeq = [];

///////////////////////
// GTD MINLP problem //
///////////////////////

// The integer variables are x(1) x(2) integers
// f(x_sol) = 7.779e-5
// f(x_sol_relax) = 9.855e-20

function y = MINLPGTD_mip_var()
y = [1 2 3 4];
endfunction

function y = MINLPGTD_x_init()
y = ones(4,1);
endfunction

function y = MINLPGTD_x_sol()
y = [31 12 12 31]';
endfunction

function y = MINLPGTD_x_sol_relax()
y = [31.592 12 12 31.592]';
endfunction

function Res = min_bd_MINLPGTD()
Res = [12 12 12 12]';
endfunction

function Res = max_bd_MINLPGTD()
Res = [60 60 60 60]';
endfunction

function y = MINLPGTDobj(x)
y = (1/6.931 - x(3)*x(2) / (x(1)*x(4)));
endfunction

MINLPGTDineq = [];

MINLPGTDeq = [];

//////////////////////////
// AVGAS1 MINLP problem //
//////////////////////////

// The integer variables are x(1) to x(8) integers
// f(x_sol) = -4
// f(x_sol_relax) = -8.114

function y = MINLPAVGAS1_mip_var()
y = [1 2 3 4 5 6 7 8];
endfunction

function y = MINLPAVGAS1_x_init()
y = ones(8,1);
endfunction

function y = MINLPAVGAS1_x_sol()
y = [1 0 0 0 0 0 1 0]';
endfunction

function y = MINLPAVGAS1_x_sol_relax()
y = [0.241 0.759 0.389 0.473 0.5 0.095 0.870 0]';
endfunction

function Res = min_bd_MINLPAVGAS1()
Res = -100*ones(8,1);
endfunction

function Res = max_bd_MINLPAVGAS1()
Res = 100*ones(8,1);
endfunction

function y = MINLPAVGAS1obj(x)
G          = 4*ones(8,8);
G(1:7,1:7) = -1*ones(7,7);
G(2:8,2:8) = -1*ones(7,7);
y = -2*x(1) - x(2) - 2*x(3) - 3*x(4) - 4*x(5) - 5*x(6) - 6*x(7) + 8*x(8) + 1/2*x'*G*x;
endfunction

function y = MINLPAVGAS1ineq(x)
y(1)  = -(-x(1) - x(2) + 1);
y(2)  = -(-x(3) - x(4) + 1);
y(3)  = -(-x(5) - x(6) + 1);
y(4)  = -(-x(7) - x(8) + 1);
y(5)  = -(-x(1) - x(3) - x(5) - x(7) + 2);
y(6)  = -(-x(2) - x(4) - x(6) - x(8) + 2);
y(7)  = -(2*x(1) + x(3) - x(7));
y(8)  = -(5*x(1) + 3*x(3) - 3*x(5) - x(7));
y(9)  = -(x(2) - x(4) - 3*x(6) - 5*x(8));
y(10) = -(x(2) - 3*x(6) - 2*x(8));
endfunction

MINLPAVGAS1eq = [];

//////////////////////////
// AVGAS2 MINLP problem //
//////////////////////////

// The integer variables are x(1) to x(8) integers
// f(x_sol) = -4
// f(x_sol_relax) = -6.631

function y = MINLPAVGAS2_mip_var()
y = [1 2 3 4 5 6 7 8];
endfunction

function y = MINLPAVGAS2_x_init()
y = ones(8,1);
endfunction

function y = MINLPAVGAS2_x_sol()
y = [0 0 1 0 0 0 1 0]';
endfunction

function y = MINLPAVGAS2_x_sol_relax()
y = [0.403 0.398 0.146 0.303 0.5 0.032 0.951 0]';
endfunction

function Res = min_bd_MINLPAVGAS2()
Res = -100*ones(8,1);
endfunction

function Res = max_bd_MINLPAVGAS2()
Res = 100*ones(8,1);
endfunction

function y = MINLPAVGAS2obj(x)
G          = 4*ones(8,8);
G(1:7,1:7) = 1*ones(7,7);
G(2:8,2:8) = 1*ones(7,7);
y = -2*x(1) - x(2) - 2*x(3) - 3*x(4) - 4*x(5) - 5*x(6) - 6*x(7) + 8*x(8) + 1/2*x'*G*x;
endfunction

function y = MINLPAVGAS2ineq(x)
y(1)  = -(-x(1) - x(2) + 1);
y(2)  = -(-x(3) - x(4) + 1);
y(3)  = -(-x(5) - x(6) + 1);
y(4)  = -(-x(7) - x(8) + 1);
y(5)  = -(-x(1) - x(3) - x(5) - x(7) + 2);
y(6)  = -(-x(2) - x(4) - x(6) - x(8) + 2);
y(7)  = -(2*x(1) + x(3) - x(7));
y(8)  = -(5*x(1) + 3*x(3) - 3*x(5) - x(7));
y(9)  = -(x(2) - x(4) - 3*x(6) - 5*x(8));
y(10) = -(x(2) - 3*x(6) - 2*x(8));
endfunction

MINLPAVGAS2eq = [];

