#ifndef  __VECT_CPP_
#ifdef INCLUDE_TEMPLATE

#ifdef INCLUDE_TEMPLATE_DEFINITION
#define __VECT_CPP_
#endif

#else 

// no include in vect.h
// -- no include in vect --
//  so reverve the sens of 
//  after vect
//  if INCLUDE_TEMPLATE_DEFINITION => no indeed of vect.cpp
#include "config.h" 
#ifndef INCLUDE_TEMPLATE_DEFINITION
#define __VECT_CPP_
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vect.h>

using namespace std;
#endif



#endif

#ifdef   __VECT_CPP_

//-----------------------------------------------------
template <class T>  A<T>:: A(int csize)
{
  cc = 0;
  size = csize;
  if (size > 0 ) { 
    cc = new T[size]; 
    assert(cc);
  }
} 
//-----------------------------------------------------
template <class T>  A<T>:: A(int csize, T& b)
{
  cc = 0;
  size = csize;
  if (size > 0 ) { 
    cc = new T[size]; assert(cc);
    for(int i=0; i<size; i++)
      cc[i] = b;
  }
} 
//-----------------------------------------------------
template <class T>  A<T>::A(const A<T>& a)
{
  if( a.cc && a.size )
    { cc = 0; size = a.size; cc = new T[size]; assert(cc);
      for(int i=0; i<size;i++) 
	cc[i] = a.cc[i];
    }
}
//-----------------------------------------------------
template <class T> void A<T>::init(int ssize)
{
  assert( !cc && ssize );
  size=ssize; 
  cc= new T[size];              
  assert(cc != 0);
#ifndef NDEBUG  
  //      cout << size <<" created: "<<cc << endl;     
  inspec();   
#endif
} 
//-----------------------------------------------------
template <class T> void A<T>::resize(int ssize)
{
  assert( ssize );
  if(!cc) init(ssize);
  else if(ssize != size) 
    {
      destroy();
      size=ssize; 
      cc= new T[size];              
      assert(cc != 0);
    }
} 
//-----------------------------------------------------
template <class T> T&  A<T>::operator [] (int i) const
{
  assert ( cc&&(i >= 0) && (i < size) );
  return cc[i];
}    
//-----------------------------------------------------
template <class T> A<T>& A<T>::operator=(const A<T>& a)
{
  assert( cc && a.cc && (a.size==size) );
  for(int i=0; i<size;i++) cc[i] = a.cc[i];
  return *this;
}
//-----------------------------------------------------
template <class T> void A<T>::destroy() {if(cc){ 
    // cout << " destroy " << cc << endl;
    if (cc) delete [] cc;size = 0; cc=0;} 
}


/*************************************************
 *                                                *
 *            IMPLEMENTATIONS                     *
 *                                                *
 *************************************************/
const float eps = (float)1.0e-20;

//---IMPLEMENTATION for VectN ---

template <class T, int N> VectN<T,N>::VectN()	
{
  for(int i=0;i<N;i++) val[i]=0;
} 

template <class T, int N> VectN<T,N>::VectN(const VectN<T,N>& r)
{ 	for(int i=0;i<N;i++) 
    val[i]=r.val[i];
} 

template <class T, int N> VectN<T,N>::VectN(T r)	
{
  for(int i=0;i<N;i++) val[i]=r;
} 

template <class T, int N> VectN<T,N> 	VectN<T,N>::operator* ( const T& a)
{
  VectN<T,N> c;
  for(int i=0;i<N;i++) 
    c.val[i] = val[i] * a ; 
  return c; 
}   

template <class T, int N> VectN<T,N> 	VectN<T,N>::operator/ (const T& a)
{
  VectN<T,N> c;
  for(int i=0;i<N;i++) 
    c.val[i] = val[i] / a ; 
  return c; 
}   

template <class T, int N> VectN<T,N>& 	VectN<T,N>::operator+= (const VectN<T,N>& a)
{ 
  for(int i=0;i<N;i++) 
    val[i]+= a.val[i]; 
  return *this; 
}   

template <class T, int N> VectN<T,N> 	VectN<T,N>::operator+ (const VectN<T,N>& a)
{ 
  VectN<T,N> c(a);
  for(int i=0;i<N;i++) 
    c.val[i]+= val[i]; 
  return c; 
}   

template <class T, int N> VectN<T,N>& 	VectN<T,N>::operator-= (const VectN<T,N>& a)
{ 
  for(int i=0;i<N;i++) 
    val[i]-= a.val[i]; 
  return *this; 
}   

template <class T, int N> VectN<T,N> 	VectN<T,N>::operator- (const VectN<T,N>& a)
{ 
  VectN<T,N> c(*this);
  for(int i=0;i<N;i++) 
    c.val[i]-= a.val[i]; 
  return c; 
}   

template <class T, int N> VectN<T,N>& VectN<T,N>::operator=	(const VectN<T,N>& a)
{
  for(int i=0;i<N;i++) 
    val[i] = a.val[i]; 
  return *this; 
}     

template <class T, int N> VectN<T,N>& VectN<T,N>::gauss(const MatN<T,N>& b)
{
  int i,j,k;	 
  T s, s1, s2;
  float saux, smin = float(1e20), eps = 1/smin;
  MatN<T,N> a(b);
  for (i=0;i<N;i++) 
    {
      for(j=0;j<=i;j++)
	{
	  s=0; for (k=0; k<j;k++)  s += a(i,k)*a(k,j) ;
	  a(i,j) -= s ;
	}
      for(j=i+1;j<N;j++)
	{
	  s=0; for (k=0;k<i;k++)  s += a(i,k)*a(k,j) ;
	  s1 = a(i,i);
	  if(saux=norm2(s1), saux< smin) smin=saux;
	  if(saux < eps) {s1 = eps;cout<<"small Gauss subpivot"<<s1<<endl;}
	  a(i,j) = ( a(i,j) - s)/s1;
	}
    }
  for (i=0;i<N;i++)							
    {	
      s2=0; for (k=0;k<i;k++)  s2 += a(i,k) * val[k];
      val[i] = (val[i] - s2) /  a(i,i) ;
    }
  for (i=N-1;i>=0;i--)
    {
      s2=0; for (k=i+1; k<=N-1;k++)  s2 += a(i,k) * val[k];
      val[i] -= s2 ;
    }
  return *this;	
}

template <class T,  int N> VectN<T,N> 	VectN<T,N>::operator/ (const MatN<T,N>& a)
{
  VectN<T,N> c(0);
  if (N==1) { c.val[0] = val[0] / a.val[0][0]; return c;}
  else if(N==2)
    {
      T s = a.val[0][0] * a.val[1][1] - a.val[0][1]* a.val[1][0];
      if(norm2(s)<eps){ s = eps;}
      c.val[0] = (val[0] * a.val[1][1] - a.val[0][1] * val[1]) / s;
      c.val[1] = (a.val[0][0] * val[1] - val[0] * a.val[1][0]) / s;
      return c;
    }
  else
    {
      for(int i =0;i<N;i++) c[i] = val[i];
      return c.gauss(a);
    }
} 

/**********************************************************/


//---IMPLEMENTATION for MatN ---
template <class T,  int N> MatN<T,N>::MatN()	
{ 	
  for(int i=0;i<N;i++) 
    for(int j=0;j<N;j++) 
      val[i][j] = 0;
} 

template <class T, int N> MatN<T,N>::MatN(const MatN<T,N>& r)
{ 	
  for(int i=0;i<N;i++) 
    for(int j=0;j<N;j++) 
      val[i][j]=r.val[i][j];
} 

template <class T, int N> MatN<T,N>::MatN(const T & r)
{ 	
  for(int i=0;i<N;i++) 
    for(int j=0;j<N;j++) 
      val[i][j]=r;
} 

template <class T, int N> MatN<T,N>::MatN( const VectN<T,N>& r)
{
  for(int i=0;i<N;i++) 
    for(int j=0;j<N;j++) 
      if(i==j) val[i][j] = r[i]; else val[i][j]=0;
} 

template <class T,  int N> MatN<T,N>& 	MatN<T,N>::operator+= (const MatN<T,N>& a)
{ 
  for(int i=0;i<N;i++) 
    for(int j=0;j<N;j++) 
      val[i][j] += a.val[i][j];
  return *this; 
}   

template <class T,  int N> MatN<T,N>& 	MatN<T,N>::operator-= (const MatN<T,N>& a)
{
  for(int i=0;i<N;i++) 
    for(int j=0;j<N;j++) 
      val[i][j] -= a.val[i][j];
  return *this; 
}   

template <class T,  int N> MatN<T,N> 	MatN<T,N>::operator+ (const MatN<T,N>& a)
{
  MatN<T,N> c(a);
  for(int i=0;i<N;i++) 
    for(int j=0;j<N;j++) 
      c.val[i][j] += val[i][j];
  return c; 
}   

template <class T,  int N> MatN<T,N> 	MatN<T,N>::operator- (const MatN<T,N>& a)
{
  MatN<T,N> c(*this);
  for(int i=0;i<N;i++) 
    for(int j=0;j<N;j++) 
      c.val[i][j] -= a.val[i][j];
  return c; 
}   

template <class T,  int N> MatN<T,N> 	MatN<T,N>::operator* (const MatN<T,N>& a)
{
  MatN<T,N> c ;
  for(int i=0;i<N;i++) 
    for(int j=0;j<N;j++) 
      for(int k=0;k<N;k++) 
	c.val[i][j] += val[i][k] * a.val[k][j];
  return c; 
}   

template <class T,  int N> MatN<T,N> 	MatN<T,N>::operator* (const T& a)
{ 
  MatN<T,N> c;
  for(int i=0;i<N;i++) 
    for(int j=0;j<N;j++) 
      c.val[i][j] =  val[i][j] * a;
  return c; 
}   

template <class T,  int N> MatN<T,N> 	MatN<T,N>::operator/ (const T& a)
{
  MatN<T,N> c;
  for(int i=0;i<N;i++) 
    for(int j=0;j<N;j++) 
      c.val[i][j] =  val[i][j] / a;
  return c; 
}   

template <class T,  int N> MatN<T,N> 	MatN<T,N>::operator/ (const MatN<T,N>& a)
{ // left inverse:  C = A / B means  BC = A !!
  MatN<T,N> c;
  
  if (N==1){ c.val[0][0] = val[0][0] / a.val[0][0]; return c;}
  else if (N==2)
    {
      T s = a.val[0][0] * a.val[1][1] - a.val[0][1]* a.val[1][0];
      if(norm2(s)<eps){ s = eps;}
      c.val[0][0] = (val[0][0] * a.val[1][1] - a.val[0][1] * val[1][0]) / s;
      c.val[1][0] = (a.val[0][0] * val[1][0] - val[0][0] * a.val[1][0]) / s;
      c.val[0][1] = (val[0][1] * a.val[1][1] - a.val[0][1] * val[1][1]) / s;
      c.val[1][1] = (a.val[0][0] * val[1][1] - val[0][1] * a.val[1][0]) / s;
      return c;
    }
  else
    {
      int i,j;	
      VectN<T,N> f;
      for(j=0; j<N; j++)
 	{
	  for(i=0;i<N;i++) f[i] = val[i][j];
	  f = f / a;
	  for(i=0;i<N;i++) c(i,j) = f[i];
 	}
      return c;
    }
  
} 

template <class T,  int N> VectN<T,N>   MatN<T,N>::operator*( const VectN<T,N>& x)
{
  VectN<T,N> b;
  for(int i=0;i<N;i++) 
    for(int j=0;j<N;j++) 
      b.val[i] += val[i][j] * x.val[j] ;
  return b; 
}     

template <class T,  int N> MatN<T,N>& MatN<T,N>::operator=(const MatN<T,N>& a)
{
  for(int i=0;i<N;i++) 
    for(int j=0;j<N;j++) 
      val[i][j]=a.val[i][j];
  return *this; 
}     

template <class T,  int N> MatN<T,N>& MatN<T,N>::operator=(const T& a)
{
  for(int i=0;i<N;i++) 
    for(int j=0;j<N;j++) 
      if(i==j) val[i][j]=a; else val[i][j]= T(0);
  return *this; 
}


template <class T,  int N> MatN<T,N>& MatN<T,N>::operator=(const VectN<T,N>& a)
{
  for(int i=0;i<N;i++) 
    for(int j=0;j<N;j++) 
      if(i==j) val[i][j]=a.val[i]; else val[i][j]=0;
  return *this; 
}     

template <class T,  int N> MatN<T,N>& 	MatN<T,N>::operator+= (const VectN<T,N>& a)
{ 
  for(int i=0;i<N;i++) 
    val[i][i] += a.val[i];
  return *this; 
}   

template <class T,  int N> MatN<T,N>& 	MatN<T,N>::operator-= (const VectN<T,N>& a)
{
  for(int i=0;i<N;i++) 
    val[i][i] -= a.val[i];
  return *this; 
}   
     
template <class T,  int N>float MatN<T,N>::modul()
{
  float c=0.F, caux; 
  for(int i=0;i<N;i++) if(caux = norm2(val[i][i]), caux>c) c = caux; 
  return (float)sqrt(c);
}

/**********************************************************/


//---IMPLEMENTATION for Vector ---
/*  define in the class A<T> 
    template <class T> Vector<T>::Vector(int csize): A<T>(csize)
    { 	const T zero(0);
    for(int i=0;i<size;i++) cc[i] = zero;
    } 
    template <class T> Vector<T>::Vector(int csize, T& b): A<T>(csize)
    { 	const T zero(0);
    for(int i=0;i<size;i++) cc[i] = b;
    } 

    template <class T> Vector<T>::Vector( const Vector<T>& v )
    {	size = v.size; 
    if (size > 0 ) cc = new T[size];  else cc = NULL;
    for(int i = 0; i<size;i++) cc[i] = v.cc[i];
    }
*/
template <class T> const Vector<T>&  Vector<T>::operator = ( const T& r)
{	
  for (int i=0; i < Vector<T>::size;i++) Vector<T>::cc[i] = r; 
  return *this;  
}

template <class T> const Vector<T>& Vector<T>::operator = ( const Vector<T>& v)
{	
  assert(v.size==Vector<T>::size); 
  for (int i=0; i < Vector<T>::size;i++) Vector<T>::cc[i] = v.cc[i]; 
  return *this;  
}

template <class T> Vector<T>& Vector<T>::operator +=(const Vector<T>& v1)
{	 
  assert(v1.size == Vector<T>::size);
  for(int i=0;i < v1.size;i++) Vector<T>::cc[i] += v1.cc[i]; 
  return *this;
}

template <class T> Vector<T> Vector<T>::operator +(const Vector<T>& v1)
{	 
  assert(v1.size == Vector<T>::size);
  Vector<T> v(v1);
  for(int i=0;i < v1.size;i++) v.cc[i] += Vector<T>::cc[i]; 
  return v;
}

template <class T> Vector<T>& Vector<T>::operator -=(const Vector<T>& v1)
{	 
  assert(v1.size == Vector<T>::size);
  for(int i=0;i < v1.size;i++) Vector<T>::cc[i] -= v1.cc[i]; 
  return *this;
}

template <class T> Vector<T> Vector<T>::operator -(const Vector<T>& v1)
{	 
  assert(v1.size == Vector<T>::size);
  Vector<T> v(*this);
  for(int i=0;i< v1.size;i++) v.cc[i] -= v1.cc[i]; 
  return v;
}

template <class T> Vector<T> Vector<T>::operator *(const float& a)
{	 
  Vector<T> v(Vector<T>::size);
  for(int i=0;i< Vector<T>::size;i++) v.cc[i] = Vector<T>::cc[i] * a; 
  return v;
}

template <class T> Vector<T> Vector<T>::operator /(const float& a)
{	 
  Vector<T> v(Vector<T>::size);
  for(int i=0;i< Vector<T>::size;i++) v.cc[i] =  Vector<T>::cc[i] / a; 
  return v;
}

template <class T> void Vector<T>::show()
{
  for(int i=0; i< Vector<T>::size;i++){ cout<<i<<":"<<Vector<T>::cc[i]<<endl;}
  cout<<"**"<<endl;
}
//-----------------------------------------------------
template <class T> T Vector<T>::min() const
{
  float tmin=+1e30;
  for (int i=0; i< Vector<T>::size; i++)
    if (Vector<T>::cc[i]<tmin) tmin=Vector<T>::cc[i];
  return tmin;
}
//-----------------------------------------------------
template <class T> T Vector<T>::max() const
{
  float tmax=-1e30;
  for (int i=0; i< Vector<T>::size; i++)
    if (Vector<T>::cc[i]>tmax) tmax=Vector<T>::cc[i];
  return tmax;
}
//-----------------------------------------------------
template <class T> T Vector<T>::scal (const Vector<T>& v1) const
{
  T s(0.0);
  for(int i=0; i< Vector<T>::size;i++) s += v1.cc[i] * Vector<T>::cc[i];
  return s;
}

/**********************************************************/


//---IMPLEMENTATION for Profilmatrix<T,R> ---
template <class T, class R> Profilmatrix<T,R>::Profilmatrix(int n,A<int>& alow, A<int>& ajlow, A<int>& ajhigh)
  : low(alow), jlow(ajlow), jhigh(ajhigh)
{ 	size=n; csize = alow[n]; 
  if (csize > 0) 
    cc = new T[csize];  else cc = NULL;
}

template <class T, class R> Profilmatrix<T,R>::Profilmatrix(Profilmatrix<T,R>& b)
  : low(b.low), jlow(b.jlow), jhigh(b.jhigh)
{ 	size=b.size; csize = b.csize;
  if (csize > 0 )
    { 	cc = new T[csize];  for(int i=0;i<csize;i++)cc[i]=b.cc[i]; }
  else cc = NULL;
}


template <class T, class R> int Profilmatrix<T,R>::resize(Grid* g)
{ 	
  if( low[size] != g->low[g->nv] || size != g->nv )
    { 	delete [] cc;
      csize = g->low[g->nv]; cc = new T[csize];  for(int i=0;i<csize;i++)cc[i]=0;
      low.destroy(); low.init(g->low.size); low = g->low;
      jlow.destroy(); jlow.init(g->jlow.size); jlow = g->jlow;
      jhigh.destroy(); jhigh.init(g->jhigh.size); jhigh = g->jhigh;
      return 1;
    }
  else return 0;
}

template <class T, class R> void Profilmatrix<T,R>::destroy()
{
  delete [] cc;
}

template <class T, class R> 
const Profilmatrix<T,R>& Profilmatrix<T,R>::operator = ( const Profilmatrix<T,R>& b)
{	
  assert(b.csize==csize); 
  for (int i=0; i<csize;i++) cc[i] = b.cc[i]; 
  return *this;  
}

template <class T, class R> 
Vector<R> Profilmatrix<T,R>::operator * (const Vector<R>& x )
{
  assert(x.size==size);
  Vector<R> b(size); 
  for (int i=0; i<size;i++) 
    {	
      b[i] = R(0);
      for(int j=low[i];j<=jhigh[i];j++)
	b[i] += cc[low[i] + j - jlow[i]]*x[j];
    }
  return b;
}

template <class T, class R> void Profilmatrix<T,R>::show()
{	
  for (int i=0;i<size;i++) 
    {  	for(int j=jlow[i];j<=jhigh[i];j++)
	cout<< i<<","<<j<<":"<<cc[low[i] + j - jlow[i]] <<"	";
      cout<<endl;
    }
  cout <<"*********"<<endl;
}

template <class T, class R> void Profilmatrix<T,R>::zero()
{	T zero(0.0);
  for (int i=0;i<csize;i++) 
    cc[i] = zero;
}

#endif
#endif
