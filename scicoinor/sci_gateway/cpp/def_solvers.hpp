//////////////////////////////
// Available linear solvers //
//////////////////////////////

// Because MUMPSis free and opensource, we use only MUMPS
// But if you have access to another solve, you can uncomment
// one or several solvers:

// * ma27: use the Harwell routine MA27
// * ma28: use the Harwell routine MA28
// * ma57: use the Harwell routine MA57
// * pardiso: use the Pardiso package
// * wsmp: use WSMP package
// * mumps: use MUMPS package

//#define USE_MA27    1
//#define USE_MA57    1
//#define USE_PARDISO 1
//#define USE_WSMP    1
//#define USE_MA28    1
#define USE_MUMPS   1
