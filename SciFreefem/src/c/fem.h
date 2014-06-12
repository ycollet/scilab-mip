#ifndef FEM_H
#define FEM_H

#include "config.h"
#include <iostream>

using namespace std;

#include "vect.h"

//extern const int next[4];
#define mmax(a,b)(a>b?a:b)
#define mmin(a,b)(a<b?a:b)
#define abss(a)(a >= 0 ? a : -(a))
#include "analyse.h"

extern int getMatProfil;

/*
  template <class T> 
  Vector<T> getrhs(Grid& g, float* f);
  template <class T> 
  Vector<T>  intgamma(Grid& g, float* f, int where);
  template <class T, class R> 
  void intmatgamma(Grid& g, Profilmatrix<T,R>& aa, float* f, int where);
  template <class T, class R> 
  void buildmatlaplace(Grid& g, Profilmatrix<T,R>& aa, float* dis
  , float* dif, float* pdx, float* pdy, float* asym, float* pdxy, float* pdyx);
  template <class T, class R> 
  void buildmatstokes (Grid& g, float nu, float alpha, float beta, Vector<R>&,
  Profilmatrix<T,R>& aa );
  template <class T, class R> 
  void buildmatlame (Grid& g, float lam, float mu, float alpha,
  Profilmatrix<T,R>& aa );
  template <class T, class R> 
  float gaussband (Profilmatrix<T,R>& a, Vector<R>& x, int first);

  class EDP;
  template <class M, class V, int N> 
  void solveprofil(Profilmatrix<M,V>* bb, EDP* edp, int old);
  template <class M, class V> 
  float gaussprofil (Profilmatrix<M,V>& a, Vector<V>& x, int first);
*/
class EDP;
template <class M, class V, int N> 
  void solveprofil(Profilmatrix<M,V>* bb, EDP* edp, int old);

static const int next[4] = {1,2,0,1};

class GetMatrixProfil
{ 
 public:
  double * Matrix;
  int * jlow, * jhigh;
  int SizeMatrix, SizeVect,SizeBlock;
 GetMatrixProfil(int pSizeMatrix, int pSizeVect,
		 int pSizeBlock) : SizeMatrix(pSizeMatrix),
    SizeVect(pSizeVect),
    SizeBlock(pSizeBlock)
    {
      Matrix = new double [SizeMatrix];
      jlow = new int [pSizeVect];
      jhigh = new int [pSizeVect];
    }

  ~GetMatrixProfil() 
    { 
      if (jhigh) delete [] jhigh; 
      if (jlow) delete [] jlow;
      if (Matrix) delete [] Matrix;
    } 
};

static GetMatrixProfil * __gmp = NULL;
#ifdef INCLUDE_TEMPLATE_DEFINITION
#define INCLUDE_TEMPLATE
#include "fem.cpp"
#undef INCLUDE_TEMPLATE
#endif
#endif
