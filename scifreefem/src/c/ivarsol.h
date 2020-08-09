#ifndef __IVARSOLVE_H
#define __IVARSOLVE_H

#include "fem.h"
// for erreur() routine
#include "rgraph.h"
//
//const int nMaxEDP = 5;  // max number of PDEs in a solve block

template <class M, class V, int N> class Isolve: public Instr
{ // solves a PDE written in strong form
  Analvar* an;				
  Iden *idmesh, *idmatrix; 	// matrix name and mesh name
  Expr* lu;					// with  LU or not LU expression
  EDP* edp;					// the PDE bloc
  Instr* l;				 	// all pdes def and bdy cond
 public:
  Profilmatrix<M,V>* bbn;
 Isolve(Iden* id0,Iden* id1, Expr* ff, EDP* eedp, Instr* ll, Analvar* aan)
   : idmesh(id0), idmatrix(id1), lu(ff), edp(eedp), l(ll), an(aan) { }
  float addmulop(int k);
  void execute ();
};


template <int N> class Ivarsolve: public Instr
{ // solves a PDE written in weak form
 public:
  typedef VectN<float,N> TN;
  typedef MatN<float,N> MN;
 private:
  Analvar* an;
  Iden *id0, *idmesh;
  CTab** f1;
  CTab** f2;
  Expr *e;
  Instr* l;
  int factorize, nedp;
  Profilmatrix<MN,TN>*  aa;
  Vector<TN>* b;
 public:
 Ivarsolve(Iden* iid0, Iden* id1, int nnedp, CTab** ff1, CTab** ff2, Expr* ee,  Instr* ll, 
	   Analvar* aan): id0(iid0),idmesh(id1), nedp(nnedp), f1(ff1), f2(ff2), e(ee), l(ll), an(aan)
  { }
  void execute ();
  void edpdoit ();
};


#ifdef INCLUDE_TEMPLATE_DEFINITION
#define INCLUDE_TEMPLATE
#include "ivarsol.cpp"
#undef INCLUDE_TEMPLATE
#endif
#endif
