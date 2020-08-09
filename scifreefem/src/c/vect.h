#ifndef VECT_H
#define VECT_H

#include "config.h"

#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

#include <assert.h>
//#define NDEBUG
#include <stdlib.h>

#ifndef NDEBUG 
//#ifdef __MWERKS1__
#undef assert
void MyAssert(int i,char * ex,char * file,long line);
#define assert(i) MyAssert((i)==0,#i,__FILE__,__LINE__)
//#endif 
#endif

#undef min
#undef max

/***************************************************/
float norm2(float a);
void inspec();  // to paliate the deficiency of the debugger

/***************************************************/
template <class T> class A{		// array with checked bounds
  //------------------------------------------------
 public: 
  T *cc;      
  int size;  
  void init(int ssize);
  void resize(int ssize);
  T&  operator [] (int i) const;
  A(int csize = 0);
  A(int csize, T& b);
  A(const A<T>& a); 
  A<T>& operator=(const A<T>& a); 
  void destroy();
  ~A() { destroy(); }  
  int no( T* t) const { return t - cc;}  // return place in array
};


/***************************************************/
class Complex
//--------------------------------------------------
{
 public: 
  float re,im;      
  
  Complex(int i)	    			{re = (float)i; im = 0.F;}
  Complex(float r=0, float i =0) 		{re =r; im =i;}  
  Complex(const Complex& a) 		{re =a.re; im =a.im;}  
  Complex& operator=(const Complex& a) 	{ re =a.re; im = a.im; return *this; }   
  Complex& operator+=(const Complex& a) 	{ re +=a.re; im += a.im;  return *this; }   
  friend Complex operator+(const Complex& a,const Complex& b) 
  { Complex c(a); c += b; return c; }   
  Complex& operator-=(const Complex& a) 	{ re -=a.re; im -= a.im;  return *this; }   
  friend ostream & operator << (ostream & os, const Complex & a)
  {
    os <<"(" <<a.re << "," << a.im << ")";
    return os;
  }
  friend Complex operator-(const Complex& a,const Complex& b) 
  { Complex c(a); c -= b; return c; }   
  friend Complex operator*(const Complex& a,const Complex& b) 
    { Complex c(a.re*b.re-a.im*b.im,a.re*b.im + a.im*b.re); return c; }   
  friend Complex operator /(const Complex& a, const Complex& b)  
  { 
    Complex c; float xl= b.re*b.re + b.im*b.im; 
    c.re =(a.re*b.re+a.im*b.im)/xl; c.im =(b.re*a.im-b.im*a.re)/xl; return c;
  }
  friend float norm2(const Complex& b) { return b.re*b.re + b.im*b.im;}
  friend float real(const Complex& b) { return b.re;}
  friend float imag(const Complex& b) { return b.im;}
};

template <class T,  int N> class MatN;


/***************************************************/
template <class T, int N> class VectN
  //--------------------------------------------------
{
 public:
  
  T val[N];     
  VectN(const VectN<T,N>& r) ; 
  VectN() ;
  VectN(T r) ;
  T& 			operator[] 	( int i)   		{ assert((i< N)&&(i>=0));  return val[i];}
  const T& 		operator[] 	( int i) const  { assert((i< N)&&(i>=0));  return val[i];}
  VectN<T,N>  	operator*	(const T& a);  
  VectN<T,N>  	operator/	(const T& a);  
  VectN<T,N>& 	operator+=	(const VectN<T,N>& a);
  VectN<T,N>   	operator+	(const VectN<T,N>& a);  
  VectN<T,N>& 	operator-=	(const VectN<T,N>& a);
  VectN<T,N>   	operator-	(const VectN<T,N>& a);  
  VectN<T,N>& 	operator=	(const VectN<T,N>& a);  
  VectN<T,N>&   gauss		(const MatN<T,N>& a);  
  VectN<T,N>   	operator/	(const MatN<T,N>& a);
  
  friend float norm2(VectN<T,N>& a)
  { 
    float s=0; 
    for(int i=0;i<N;i++) s+=norm2(a.val[i]); 
    return s;
  }  
  friend ostream & operator << (ostream & os, const VectN<T,N>& a)
  {
    os << "[";
    for(int i=0;i<N-1;i++) os<<a.val[i]<<" ";
    os<<a.val[N-1] << "]";
    return os;
  }
};

/***************************************************/
template <class T,  int N> class MatN
  //--------------------------------------------------
{
 public:   
  
  T val[N][N];     
  MatN(const MatN<T,N>& r) ; 
  MatN();
  MatN(const T & r) ;  
  MatN( const VectN<T,N>& r);
  
  T& operator() ( int i, int j) 
    { 	assert((i< N)&&(i>=0)&&(j< N)&&(j>=0));  
      return val[i][j];
    }
 
  const T& operator() ( int i, int j) const  
  { 	assert((i< N)&&(i>=0)&&(j< N)&&(j>=0));  
    return val[i][j];
  }
  
  MatN<T,N>& 	operator+=	(const MatN<T,N>& a);
  MatN<T,N>& 	operator-=	(const MatN<T,N>& a);
  MatN<T,N> operator+	(const MatN<T,N>& a);  
  MatN<T,N> operator-	(const MatN<T,N>& a);  
  MatN<T,N> operator*	(const MatN<T,N>& a);  
  MatN<T,N> operator/	(const MatN<T,N>& a);  
  MatN<T,N>& 	operator=	(const MatN<T,N>& a);  
  MatN<T,N>& 	operator=	( const T& a);  
  VectN<T,N> operator*(const VectN<T,N>& x);
  MatN<T,N> operator*	(const T& a);  
  MatN<T,N> operator/	(const T& a);  
  MatN<T,N>& 	operator=	(const VectN<T,N>& a);  
  MatN<T,N>& 	operator+=	(const VectN<T,N>& a);
  MatN<T,N>& 	operator-=	(const VectN<T,N>& a);
  float modul();
  //  friend  float norm2( MatN<T,N>& a);
  /*      { 	float c=0.; 
	  for(int i=0;i<N;i++) 
	  for(int j=0; j<N;j++) c+= norm2(a.val[i][j]);
	  return c;
	  }*/
  friend ostream & operator << (ostream & os, const MatN<T,N>& a)
  {
    os << "[";
    for(int i=0;i<N;i++) for(int j=0;j<N;j++) os<<a.val[i][j]<<" ";
    os << "]";
    return os;
  }
};

template<class T,int N>
  float norm2( MatN<T,N>& a)
{ 	float c=0.; 
  for(int i=0;i<N;i++) 
    for(int j=0; j<N;j++) c+= norm2(a.val[i][j]);
  return c;
}
/***************************************************/
template <class T> class Vector : public A<T>
				  //-------------------------------
{
 public:
 Vector(int csize):A<T>(csize){}
  ~Vector(){}
  Vector ( const Vector<T>& v ): A<T>(v){};
  Vector (int csize, T& b) : A<T>(csize,b){};
  const Vector<T>& operator = (const T& r);
  const Vector<T>& operator= ( const Vector<T>& v);
  Vector<T>& operator += (const Vector<T>& v1);
  Vector<T> operator + (const Vector<T>& v1);
  Vector<T>& operator -= (const Vector<T>& v1);
  Vector<T> operator - (const Vector<T>& v1);
  Vector<T> operator* (const float& a);
  Vector<T> operator/ (const float& a);
  void show();
  T scal(const Vector<T>& v1) const;
  T min() const;
  T max() const;
};

/***************************************************/
class Grid;
//--------------------------------------------------------------------------------
template <class T, class R> class Profilmatrix 
  //--------------------------------------------------------------------------------
{ 
 public: 

  A<int> jlow,low,jhigh;
  int size,csize;		// profil and size
  T* cc;

  Profilmatrix(){}
  Profilmatrix(int n,A<int>& alow, A<int>& ajlow, A<int>& ajhigh);
  
  Profilmatrix(Profilmatrix<T,R>& b);
  
  const Profilmatrix<T,R>& operator = ( const Profilmatrix<T,R>& b);
  int resize(Grid* g);
  void destroy();
  ~Profilmatrix(){ destroy(); }
  
  T& operator()( int i,  int j) const  // access function
    {
      int k = low[i] + j - jlow[i];
      assert(cc && (i>=0)&&(i<size) && (j>=jlow[i]) && (j<= jhigh[i]));
      return cc[k];
    }
  Vector<R> operator * (const Vector<R>& ); 
  void show();
  void zero();
};





const int nbholesmax=100;
class bVertex;
class bTriangle;
class bEdge;

class float3 { public:	float v[3];};
class bPoint { public:	float x,y;};

/***************************************************/
class bVertex {
 public:
  float x, y;       // cordinates          
  int where;        // on which boundary      
  int nsupp;        // nb of triangles which contains this vertex
  A<bTriangle*> supp;// all triangles which contain the vertex
  int nmate;        // number of neighbor vertices
  A<bVertex*> mate;  // all neighbors
  void fill(float xx, float yy, int ngg) ;
  bVertex& operator=(bVertex& v){ fill(v.x,v.y,v.where); return *this;}
};

/***************************************************/
class bTriangle {
 public:
  bVertex* v[3];  // the 3 vertices of the triangle      
  bEdge* e[3];    // pointer to the edges opposite each vertex
  int where;     //  in which region                  
  float area;
  void fill(bVertex& p0, bVertex& p1, bVertex& p2, int wwhere);
  bTriangle& operator=(bTriangle& t)
    { 	v[0]=t.v[0]; v[1] = t.v[1]; v[2] = t.v[2]; where = t.where; return *this;}
};


/***************************************************/
class bEdge { public:
  bVertex *in, *out;      // oriented by first to last vertex
  bTriangle *left, *right;// triangles on each side of edge
  int where;     // on which curve (or boundary)
  float length;          
  void fill(bVertex* iin, bVertex* oout, bTriangle*  lleft, bTriangle* rright);
  bEdge& operator=(bEdge& b)
    { where = b.where; in = b.in; out = b.out; left = b.left; right = b.right; return* this;}
};

/***************************************************/
class Triangles;
class Geometry; // cf Mesh.h (Hecht)
class frontiere;
class Iden;
/***************************************************/

class Grid {
 public:
  int NbRef ; // nb de Ref 
  int nt, nv, ne; // nb of triangles, vertices and edges
  int nbholes;    // nb of holes in the domain
  int bdth;       // bandwidth
  A<bVertex> v;    // all vertices                      
  A<bTriangle> t;  // all triangles          
  A<bEdge> e;      // all edges
  A<int> low, jlow, jhigh;	// profil (see getprofil for comments)
  A<float3> dxhat, dyhat;		// derivatives of hat functions
  A<bPoint> norml;
  
  // Communication between 2 mesh for integrals
  int**	quad;		// for each quadrature point: the triangle nb and the iloc
  Grid* gridFriend;		// the grid with respect to which quad is defined
  void initquad(Grid* t);// in grid.cpp
  
  // Communication with F Hecht data structure
  Triangles* Th;  	
  Geometry *Gh;
  int * NumThinGrid;
  void buildit(frontiere* t, float & waitm);
  void th2t(Triangles* tTh);
  void copyev(const Grid* g); // used only to recopy grids
  
 Grid():NbRef(0),nt(0),nv(0),ne(0),nbholes(nbholesmax),gridFriend(0),NumThinGrid(0) {}// default constructor
  void readgrid(const char *path );     // reads a triangulation
 Grid(const char *path ):NbRef(0),v(),t(),e(),gridFriend(0){ readgrid(path);}
    
 Grid(const Grid* g) :NbRef(0), nt(g->nt), nv(g->nv), ne(g->ne), 
    v(g->v),t(g->nt),e(g->e),low(), jlow(), jhigh(),	
    dxhat(),dyhat(), gridFriend(0), Th(g->Th), Gh(g->Gh) { copyev(g);}
  void destroy();
  ~Grid()  { destroy();NbRef=-1;}
  void save(const char* path, int debugformat) const; //save mesh in Gfem format
  void gnusave(const char* path) const; //save mesh in GNUPlot format
  void dump(const char* path) const;
  int no(bTriangle* tt) const { return t.no(tt);}
  int no(bVertex* tt) const { return v.no(tt);}
  int no(bEdge* tt) const { return e.no(tt);}//place in e of  tt
  void square(int nx, int ny); 	// triangulate a square
  void prepgrid(int dontTouchEdges);	// calls the 5 below
  void computegeom(int dontTouchEdges);	// computes the triangle areas and edges length
  void getbdth();					// computes bandwidth
  void getprofil();				// computes profil
  void fillvsupp();				// these 3 should be called in that order 
  void getnmate() const;			// computes v[].nmate
  void fillmate(int dontTouchEdges);	// computes v[].mate and e[]
  void derivhat();
  int gibbsv (long* ptvoi, long* vois,long* lvois,long* w,long* v);
  int renum();
  void reGrid(const Grid* g);
  void AddRef() {assert(this); NbRef++;}
  void DelRef() {assert(this);if (!NbRef--) delete this;}
  void check(){assert(NbRef>=0);}
  void draw(float & wait);
 private:
  void show();
 
};


/***************************************************/
class EFSpace: public Vector<float>
{ 
 public:
  Grid* g;
 EFSpace(Grid* gg,int n): Vector<float>(n), g(gg) {if (g) g->AddRef();}
  virtual ~EFSpace(){ if(g) g->DelRef();}
  void load(const char* filename);
  void save(const char* filename) const;
  virtual void gnusave(const char* filename) const=0;
  virtual float F(float x, float y) const =0;
  virtual float F(int v) const =0;
  virtual float F(int t,float x[3]) const=0;
  virtual void Vresize()  =0;
};
/***************************************************/
long FindTriangle(Triangles &Th, double x, double y, double* a,int & inside);

class P1: public EFSpace
{ 
 public:
  //  Grid* g;
 P1(Grid* gg): EFSpace(gg,gg?gg->nv:0){ }

  //  ~P1(){ destroy();}
  //  void load(const char* filename);
  //  void save(const char* filename) const;
  void gnusave(const char* filename) const;
  float F(float x, float y) const { 
    double xl[3];
    int inside =1;
    int kt=FindTriangle(*g->Th,x,y,xl,inside);
    return   cc[g->no(g->t[kt].v[0])]*xl[0] 
      + cc[g->no(g->t[kt].v[1])]*xl[1]
      + cc[g->no(g->t[kt].v[2])]*xl[2];}
       
  float F(int v) const { return cc[v];}
  float F(int kt,float xl[3]) const {      
    return   cc[g->no(g->t[kt].v[0])]*xl[0] 
      + cc[g->no(g->t[kt].v[1])]*xl[1]
      + cc[g->no(g->t[kt].v[2])]*xl[2];}

  void Vresize(){ Vector<float>::resize(g->nv);}
  //  float intgrad2(P1& a, P1& b, P1& c) const; //int a.gradf.gradf dx + b.f.f + int_gamma c.f.f
  virtual ~P1() {}
};

class P0: public EFSpace
{ 
 public:
  //  Grid* g;
 P0(Grid* gg): EFSpace(gg,gg?gg->nt:0){}
  //  ~P0(){ destroy();}
  //  void load(const char* filename);
  //  void save(const char* filename) const;
  void gnusave(const char* filename) const;
  float F(float x, float y) const { 
    double xl[3];
    int inside =1;
    int kt=FindTriangle(*g->Th,x,y,xl,inside);
    return   cc[kt];}       
  float F(int v) const ;//{ return cc[g->Th->Number(g->Th->vertices[v].t)]; }
  float F(int kt,float xl[3]) const {return   cc[kt];}

  void Vresize(){ Vector<float>::resize(g->nt);}
  //  float intgrad2(P1& a, P1& b, P1& c) const; //int a.gradf.gradf dx + b.f.f + int_gamma c.f.f
  virtual ~P0() {}
};

/***************************************************/
/*class Laplace{ public:
  Grid* g;
  P1 sol; // contains bdy cond on entry
  P1 neu; // contains neumann conditions
  P1 rhs; // right hand side
  P1 dif; // the diffusion coefficient
  P1 vis; // the viscosity coef
  P1 rob; // the robin coef
  P1 pdx; // the dx() coef
  P1 pdy; // the dy() coef
  P1 asym,pdxy,pdyx;// the second order coef (dif = (dxx+dyy)/2, (d2asym = dxx-dyy)/2)
  
  Laplace(Grid* gg): g(gg), rhs(g), sol(g), neu(g), 
  dif(g), vis(g), rob(g), pdx(g), pdy(g), asym(g),pdxy(g),pdyx(g){}
  void solvegradconj(const int niter, const float precise); 
  //  void solveprofil(Profilmatrix<float,float>* bb, int old);
  void derivener(Vector<float>& gradE) const;
  void erreur_max();
  void erreur_mid();
  };
*/
/***************************************************/
/*
  template <class T> class P1T: public Vector<T> //typically= float, complex or VectN<T,N>
  { public:
  Grid* g;
  P1T(Grid* gg): g(gg), Vector<T>(gg?gg->nv:0){}
  ~P1T(){ destroy();}

  //  void load(const char* filename);
  //  void save(const char* filename) const;
  //  void gnusave(const char* filename) const;
  //  float intgrad2(P1& a, P1& b, P1& c) const; //int a.gradf.gradf dx + b.f.f + int_gamma c.f.f
  };
*/
/***************************************************/
/*
  template <class T> class LaplaceN  //typically T= float, complex or VectN<T,N>
  { public:
  Grid* g;
  P1T<T> sol; // contains bdy cond on entry
  P1T<T> neu; // contains neumann conditions
  P1T<T> rhs; // right hand side
  P1T<T> dif; // the diffusion coefficient
  P1T<T> vis; // the viscosity coef
  P1T<T> rob; // the robin coef
  P1T<T> pdx; // the dx() coef
  P1T<T> pdy; // the dy() coef
  P1T<T> asym,pdxy,pdyx;// the second order coef (dif = (dxx+dyy)/2, (d2asym = dxx-dyy)/2)
  
  LaplaceN(Grid* gg): g(gg), rhs(g), sol(g), neu(g), 
  dif(g), vis(g), rob(g), pdx(g), pdy(g), asym(g),pdxy(g),pdyx(g){}
    
  //  void solvegradconj(const int niter, const float precise); 
  //  void solveprofil(Profilmatrix<M,V>* bb, int old);
  //  void derivener(Vector<float>& gradE) const;
  //  void erreur_max();
  //  void erreur_mid();
  };
*/
#ifdef INCLUDE_TEMPLATE_DEFINITION
#define INCLUDE_TEMPLATE
#include "vect.cpp"
#undef INCLUDE_TEMPLATE
#endif
#endif
