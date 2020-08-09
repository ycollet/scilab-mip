#ifndef __CONVEC_H_
#define __CONVEC_H_
//#include "config.h"
//#include "vect.h"
//#include "fem.h"
//#include "graph.h"
//#include "rgraph.h"

typedef VectN<float,2> TVN;
TVN convect(Grid& g, const Vector<TVN>& f, const A<float>& u, const A<float>& v, double dt, const int i1);

int searchTriangle (Grid& g, const A<float>& u, const A<float>& v, int i1,int& k,int& iloc);

int xtoX (Grid& g, const A<float>& u, const A<float>& v, double* xl, double *dt, int *k);

/*
template <class T> T convect(Grid& g, const Vector<T>& f, 
			    const A<float>& u, const A<float>& v, double dt, const int i1);	
*/
/*
#ifdef INCLUDE_TEMPLATE_DEFINITION
#define  INCLUDE_TEMPLATE
#include "convec.cpp"
#undef  INCLUDE_TEMPLATE
#endif
*/
#endif
