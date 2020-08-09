#ifndef __TNFEM_CPP_
#ifdef INCLUDE_TEMPLATE
#  ifdef INCLUDE_TEMPLATE_DEFINITION
#    define __TNFEM_CPP_
#  endif
#else
#  include "config.h" 
#  ifndef INCLUDE_TEMPLATE_DEFINITION
#    define __TNFEM_CPP_
#    include <assert.h>
//#    define NDEBUG
#    include <cmath>
#    include <stdlib.h>
#    include <iostream>
#    include <fstream>
#    include "fem.h"
#  endif
#endif

#ifdef __TNFEM_CPP_

#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

#include <assert.h>
#include <stdlib.h>

#include "fem2.h"

class EDP;

template <class T, int N> 
void cast(MatN<T,N>& b, const float* a, int k) 
{
  for(int i=0;i<N;i++)
    for(int j= 0; j<N; j++)  b.val[i][j] = a[(k*N+i)*N + j]; 
}

template <class T, int N> void cast(VectN<T,N>& b, const float* a, int k) 
{
  for(int i=0;i<N;i++)
    b.val[i] = a[k*N+i]; 
}


inline void cast(float& b, const float* a, int j) 	{ b= a[j];}

template <int N>  float cast2float(VectN<float,N> a,int j) {return a[j];}
inline float cast2float(float a, int j) {return a;}

template <class T, int N> void 	addcast(int kvar, MatN<T,N>& a, float x) 
{ a.val[kvar][kvar] += x;}

template <class T, int N> void 	addcast(int kvar, VectN<T,N>& b, float x) 
{ b[kvar] += x;}

inline void addcast(int kvar, float& b, float x) { b = x;}

template <class T> 	Vector<T> getrhs(const Grid& g, const float* f)
{
  int iloc;
  T zero(0), dum[3];
  Vector<T> r(g.nv, zero);  // put to zero at creation
  for (int k=0; k<g.nt; k++)
    {
      bTriangle& tk = g.t[k];
      for ( iloc=0; iloc<3; iloc++)
	cast(dum[iloc],f,g.no(tk.e[iloc]));
      for ( iloc=0; iloc<3; iloc++)
	r[g.no(tk.v[iloc])] += (dum[::next[iloc]]+dum[::next[iloc+1]]) *  (tk.area/6) ;
    }
  return r;
}

template <class T> Vector<T> 
intgamma(const Grid& g, const float* fin, const float* fout, int where)
{ 
  const double sq3 = 1./sqrt(3.0);
  T zero(0), dum0,dum1;// computes r_i = int_where(F w_i) on the entire bdy if where<0 
  Vector<T> r(g.nv, zero); // computed by edges using where on edges and values of F in f and g
  for (int k=0; k<g.ne; k++)
    {
      bEdge& ek = g.e[k];
      int i = g.no(ek.in);
      int ip = g.no(ek.out);
      if((ek.where == where)||(ek.where && (where<0)))
	{	
	  // function evaluated at two quadra points per segment. 
	  cast(dum0,fin,k);// first quadrature point is near e.in
	  cast(dum1,fout,k);
	  r[i] += (dum0*(1+sq3) + dum1*(1-sq3))*ek.length/4 ;
	  r[ip]+= (dum0*(1-sq3) + dum1*(1+sq3))*ek.length/4;
	  /*			  cast(dum0,f,i);
				  cast(dum1,f,ip);
				  r[i] += (dum0*2 + dum1)*ek.length/6 ;
				  r[ip]+= (dum1*2 + dum0)*ek.length/6;
	  */			}
    }
  return r;
}

// computed by triangles using where on vertices
template <class M, class V> 
void intmatgamma(const Grid& g, Profilmatrix<M,V>& aa, const float* f, int where)
{ 									// computes a(i,j)+=int_Gamma_where(f w^iw^j) 	
  for (int k=0; k<g.ne; k++)
    {
      bEdge& ek = g.e[k];
      bVertex& vi = *ek.in;
      bVertex& vip = *ek.out;
      int i = g.no(&vi);
      int ip = g.no(&vip);
      if((ek.where == where)||(ek.where && (where<0)))
	{	
	  M tr, dum0;
	  cast(dum0,f,k) ; 	
	  tr = dum0*ek.length/4 ;
	  aa(i,i) += tr ;
	  aa(ip,ip) += tr  ;
	  aa(i,ip) += tr;
	  aa(ip,i) +=tr;
	}
    }
}

template <class M, class V> 
void buildmatlaplace(const Grid& g, Profilmatrix<M,V>& aa, const float* dis
		     ,  const float* dif,  const float* pdx,  const float* pdy,
		     const float* asym, 
		     const float* pdxy,  const float* pdyx)
{  // builds a(i,j) = int(dis w_i w_j + dif grad(w_i) grad(w_j) + trans grad(w_) w_j +...
  
  M alph, beta,betaxy,betayx,betaasym, pdx3,pdy3, dum0,dum1,dum2;
  int i,j,k,ip,jp,ipp,jpp,iloc,jloc;
  float dwidxa,dwjdxa,dwidya,dwjdya;
  M aaloc;
  for ( k=0; k<aa.csize; k++) aa.cc[k] = 0.F;
  
  for ( k=0; k<g.nt; k++)
    {
      bTriangle& tk = g.t[k];
      i = g.no(tk.v[0]);
      ip = g.no(tk.v[::next[1]]);
      ipp = g.no(tk.v[::next[2]]);
      cast(alph,dis,k);  	
      cast(beta,dif,k);
      cast(betaxy,pdxy,k);
      cast(betayx,pdyx,k);  
      cast(betaasym,asym,k); 
      cast(pdx3,pdx,k);
      cast(pdy3,pdy,k); 
      for ( iloc=0; iloc<3; iloc++)
	{
	  i = g.no(tk.v[iloc]);
	  ip = g.no(tk.v[::next[iloc]]);
	  ipp = g.no(tk.v[::next[iloc+1]]);
	  dwidxa = (g.v[ip].y - g.v[ipp].y)/(tk.area * 4); 
	  dwidya = -(g.v[ip].x - g.v[ipp].x)/(tk.area * 4);
	  for ( jloc=0; jloc<3; jloc++)
	    {
	      j = g.no(tk.v[jloc]);			 
	      jp = g.no(tk.v[ ::next[jloc]]);			 
	      jpp = g.no(tk.v[::next[jloc+1]]);
	      dwjdxa = g.v[jp].y - g.v[jpp].y;
	      dwjdya = -(g.v[jp].x - g.v[jpp].x);
	      aaloc = (pdx3 * dwidxa  + pdy3 * dwidya + alph/8 )* tk.area/1.5
		- betaxy * (dwidya * dwjdxa) - betayx * (dwidxa * dwjdya)
		- betaasym * ( dwidxa * dwjdxa - dwidya * dwjdya) 
		- beta * ( dwidxa * dwjdxa + dwidya * dwjdya) ;
	      if(i != j) 
		aa(j,i) += aaloc; 
	      else
		aa(i,i) += aaloc + alph * tk.area/12;
	    }
	}
    }	
}			

template <class M, class V> 
void buildmatstokes (Grid& g, float nu, float alpha, float beta, // Vector<R>& u,
		     Profilmatrix<M,V>& aa )
{
  
  const float alph = alpha/12.0, bet = beta*0.5/nu;
  const float xla = 0., txla = 2 + xla , rxla = 1+xla; // xla = second viscosity
  int i,j,k,kp,ip,jp,ipp,jpp,iloc,jloc;
  float dwidxa,dwjdxa,dwidya,dwjdya, ugradwjwi=0;
  M aaloc;
  
  aa.zero();
  
  for ( k=0; k<g.nt; k++)
    {
      bTriangle& tk = g.t[k];
      float tka = tk.area * 2;
      for ( iloc=0; iloc<3; iloc++)
	{
	  i = g.no(tk.v[iloc]);
	  ip = g.no(tk.v[::next[iloc]]);
	  ipp = g.no(tk.v[::next[iloc+1]]);
	  for ( jloc=0; jloc<3; jloc++)
	    {
	      j = g.no(tk.v[jloc]);			 
	      jp = g.no(tk.v[ ::next[jloc]]);			 
	      jpp = g.no(tk.v[::next[jloc+1]]);
	      dwidxa = (g.v[ip].y - g.v[ipp].y)/tka;
	      dwjdxa = (g.v[jp].y - g.v[jpp].y)/tka;
	      dwidya = -(g.v[ip].x - g.v[ipp].x)/tka;
	      dwjdya = -(g.v[jp].x - g.v[jpp].x)/tka;
	      /* //other formulation for Stokes: nu(grad u + grad u^T)(grad v + grad v^T) -(div u,q) - (div v,p)
		 aaloc(0,0) = nu*(txla * dwidxa * dwjdxa +  dwidya * dwjdya) + alph ;
		 aaloc(1,1) = nu*(dwidxa * dwjdxa + txla * dwidya * dwjdya) + alph ;
		 aaloc(1,0) = nu * (xla *dwidya * dwjdxa +  dwidxa * dwjdya);
		 aaloc(0,1) =  nu * (xla *dwidxa * dwjdya +  dwidya * dwjdxa);
		 aaloc(2,2) =  bet*tk.area*(dwidxa * dwjdxa + dwidya * dwjdya) ;
		 aaloc(2,0) = - dwjdxa/3;
		 aaloc(0,2) = - dwidxa/3;
		 aaloc(1,2) = - dwidya/3;
		 aaloc(2,1) = - dwjdya/3;
		 ugradwjwi  = ( (2 * u[j][0] + u[jp][0] + u[jpp][0]) * dwidxa // convection: ugrad u
		 +(2 * u[j][1] + u[jp][1] + u[jpp][1]) * dwidya )/ 12;
	      */ 			 	aaloc(0,0) =  nu*(dwidxa * dwjdxa + rxla * dwidya * dwjdya) + alph ;//+ ugradwjwi;
	      aaloc(1,1) =  nu*(rxla * dwidxa * dwjdxa + dwidya * dwjdya) + alph ;//+ ugradwjwi;
	      aaloc(2,2) =  bet*tk.area*(dwidxa * dwjdxa + dwidya * dwjdya) + 1.e-5 ;
	      aaloc(1,0) =  nu * xla *dwidya * dwjdxa;
	      aaloc(0,1) =  nu * xla *dwidxa * dwjdya;
	      aaloc(2,0) =  dwidxa/3;
	      aaloc(2,1) =  dwidya/3;
	      aaloc(0,2) =  dwidxa/3;
	      aaloc(1,2) =  dwidya/3;
	      if(i==j) { aaloc(0,0) += alph; aaloc(1,1) += alph; aaloc(2,2) += 1.e-5;} 
	      aa(j,i) += aaloc * tk.area;
	    }
	}
    }	
}			

template <class T, class R> 
void buildmatlame (Grid& g, float lam, float mu, float alpha,
		   Profilmatrix<T,R>& aa )
{	
  int i,j,k,kp,ip,jp,ipp,jpp,iloc,jloc;
  float dwidxa,dwjdxa,dwidya,dwjdya;
  const float lammu = (lam+mu)*2, mu2 = mu*2, lam2 = lam*2, alph = alpha/12.0;;
  T aaloc;
  
  aa.zero();
  for ( k=0; k<g.nt; k++)
    {
      bTriangle& tk = g.t[k];
      float tka = tk.area * 2;
      for ( iloc=0; iloc<3; iloc++)
	{
	  i = g.no(tk.v[iloc]);
	  ip = g.no(tk.v[::next[iloc]]);
	  ipp = g.no(tk.v[::next[iloc+1]]);
	  for ( jloc=0; jloc<3; jloc++)
	    {
	      j = g.no(tk.v[jloc]);			 
	      jp = g.no(tk.v[ ::next[jloc]]);			 
	      jpp = g.no(tk.v[::next[jloc+1]]);
	      dwidxa = (g.v[ip].y - g.v[ipp].y)/tka;
	      dwjdxa = (g.v[jp].y - g.v[jpp].y)/tka;
	      dwidya = -(g.v[ip].x - g.v[ipp].x)/tka;
	      dwjdya = -(g.v[jp].x - g.v[jpp].x)/tka;
	      aaloc(0,0) = lammu * dwidxa * dwjdxa + mu2 * dwidya * dwjdya + alph ;
	      aaloc(1,1) = mu2 * dwidxa * dwjdxa + lammu * dwidya * dwjdya + alph ;
	      aaloc(1,0) = lam2 * dwidya * dwjdxa + mu2 * dwidxa * dwjdya;
	      aaloc(0,1) = lam2 * dwidxa * dwjdya + mu2 * dwidya * dwjdxa;
	      if(i==j) for(kp=0;kp<=1;kp++) aaloc(kp,kp) += alph; 
	      aa(i,j) += aaloc * tk.area;
	    }
	}
    }	
}			


template <class T, class R> float gaussband (Profilmatrix<T,R>& a, Vector<R>& x, int first)
/*----------------------------------------------------
  Factorise (first=1) and/or solve 	Ay = x  with result in x 
  LU is stored in A ; returns the value of the smallest pivot  
  all pivots less than eps are put to eps 
  a[i][j] is stored in a[i-j+bdthl][j]=a[n*(i-j+bdthl)+j] 
  where -bdwth <= i-j <= bdthl 
*/
{
  int i,j,k;	 
  T s, s1;
  R s2, rzero(0);
  float saux, smin = float(1.e9), eps = 1/smin;
  int n = a.size;
  int bdthl = a.bdth;
  
  if (first)			
    for (i=0;i<n;i++) 
      {
	for(j=mmax(i-bdthl,0);j<=i;j++)
	  {
	    s=0.F; for (k=mmax(i-bdthl,0); k<j;k++)  s += a(i,k)*a(k,j) ;
	    a(i,j) -= s ;
	  }
	for(j=i+1;j<=mmin(n-1,i+bdthl);j++)
	  {
	    s= a(i,j); for (k=mmax(j-bdthl,0);k<i;k++)  s -= a(i,k)*a(k,j) ;
	    s1 = a(i,i);
	    if(saux=norm2(s1), saux < smin) smin=saux;
	    if(saux < eps) s1 = eps;
	    a(i,j) = s / s1;
	  }
      }
  x.show();
  for (i=0;i<n;i++)							
    {	
      s2 = x[i];  for (k=mmax(i-bdthl,0);k<i;k++)  s2 -= a(i,k) * x[k];
      x[i] = s2 /  a(i,i) ;
    }
  x.show();
  for (i=n-1;i>=0;i--)
    {
      s2=rzero; for (k=i+1; k<=mmin(n-1,i+bdthl);k++)  s2 += a(i,k) * x[k];
      x[i] -= s2 ;
    }
  x.show();
  return smin;	
}

template <class M, class V> float gaussprofil (Profilmatrix<M,V>& a, Vector<V>& x, int first)
/*----------------------------------------------------
  Factorise (first=1) and/or solve 	Ay = x  with result in x 
  LU is stored in A ; returns the value of the smallest pivot  
  all pivots less than eps are put to eps 
*/
{
  int i,j,k;	 
  M s, s1;
  V s2, rzero(0);
  float saux, smin = float(1.e9), eps = 1/smin;
  int n = a.size;
  
  if (first)			
    for (i=0;i<n;i++) 
      {
	for(j=a.jlow[i];j<=i;j++)
	  {
	    s=0.F; for (k=a.jlow[i]; k<j;k++)
		     if((a.jlow[k]<=j)&&(j<=a.jhigh[k]))  s += a(i,k)*a(k,j) ;
	    a(i,j) -= s ;
	  }
	for(j=i+1;j<=a.jhigh[i];j++)
	  {
	    s= a(i,j); for (k=a.jlow[i];k<i;k++)
			 if((a.jlow[k]<=j)&&(j<=a.jhigh[k]))  s -= a(i,k)*a(k,j) ;
	    s1 = a(i,i);
	    if(saux=norm2(s1), saux < smin) smin=saux;
	    if(saux < eps) s1 = eps;
	    a(i,j) = s / s1;
	  }
      }

  for (i=0;i<n;i++)							
    {	
      s2 = x[i];  for (k=a.jlow[i];k<i;k++)  s2 -= a(i,k) * x[k];
      x[i] = s2 /  a(i,i) ;
    }
  for (i=n-1;i>=0;i--)
    {
      s2=rzero; for (k=i+1; k<=a.jhigh[i];k++)  s2 += a(i,k) * x[k];
      x[i] -= s2 ;
    }
  return smin;	
}

template <class T,int N>
void ffcopy(int & i, double ** m, const MatN<T,N> & mat)
{
  int k,l;
  for (k = 0; k < N; k++)
    for (l = 0; l < N; l++)
      (*m)[i++] = mat(k,l);
}

template <class M, class V, int N> 
void solveprofil(Profilmatrix<M,V>* bb, EDP* edp, int old)
{
  //	Vector<VectN<float,1>> f(nv),sol(nv);
  //	Bandmatrix<MatN<float,1>,VectN<float,1>> aa(nv,g.bdth);
  int i;
  Grid& t=*(edp->g);
  float penal=1e10;
  Profilmatrix<M,V>& aa = *bb;
  
  Vector<V> aux(edp->g->nv);	
  
  if(!old){
    buildmatlaplace<M,V>(t,aa, edp->dis,edp->dif,edp->pdx,edp->pdy,edp->asym,edp->pdxy,edp->pdyx);
    intmatgamma<M,V>(t,aa,edp->rob,-1);
    for(int  i=0;i<t.nv;i++)
      for(int kvar =0; kvar<edp->n;kvar++)
      	if(fabs(edp->sol[i*edp->n +kvar])>0)
	  addcast(kvar,aa(i,i),penal); //Dirichlet
  }
  aux = getrhs<V>(t,edp->rhs) + intgamma<V>(t,edp->neuin,edp->neuout,-1);
  //  aux.show();
  for( i=0;i<t.nv;i++)
    for(int kvar =0; kvar<edp->n;kvar++)
      if(fabs(edp->sol[i*edp->n +kvar])>0)
	addcast(kvar,aux[i],penal*edp->sol[i*edp->n +kvar]);
  
  float pivot = gaussprofil<M,V>(aa,aux,!old) ;
  if(!old) cout <<"\t\t"  << " pivot=" << pivot << endl;
  
  for(int j=0; j<N;j++)
    {
      float maxsol = cast2float(aux[0],j);
      float minsol = maxsol;
      for(  i=0;i<t.nv;i++) 
  	{ 
	  int k = i*edp->n+j;
	  edp->sol[k] = cast2float(aux[i],j);
	  if(maxsol < edp->sol[k] ) maxsol=edp->sol[k];
	  if(minsol > edp->sol[k] ) minsol=edp->sol[k];
  	}
      cout<< "\t min= "<< minsol << "\t max= " << maxsol <<endl;
    }

  if (getMatProfil) 
    {
      int ind = 0;

      __gmp = new GetMatrixProfil(aa.csize*N*N,aa.size,N);
      for (i = 0; i < aa.size; i++)
	__gmp -> jlow[i] = aa.jlow[i];
      for (i = 0; i < aa.size; i++)
	__gmp -> jhigh[i] = aa.jhigh[i];
      for (i = 0; i < aa.csize; i++)
        ffcopy(ind,&(__gmp -> Matrix),(aa.cc)[i]);
    }
}

/*
  template <class T> void deriv(Grid& g, const Vector<T>& f, Vector<T>& rx, Vector<T>& ry)
  { // computes derivatives of f
  for (int k=0; k<g.nt; k++)
  {
  bTriangle& tk = g.t[k];
  float tka = tk.area * 2;
  rx[k] = 0;	ry[k] = 0;
  for (int iloc=0; iloc<3; iloc++)
  {
  int i = g.no(tk.v[iloc]);
  int ip = g.no(tk.v[::next[iloc]]);
  int ipp = g.no(tk.v[::next[iloc+1]]);
  float dwidx = (g.v[ip].y - g.v[ipp].y)/tka;
  float dwidy  = -(g.v[ip].x - g.v[ipp].x)/tka;
  rx[k] += f[i] * dwidx;
  ry[k] += f[i] * dwidy;
  }
  }
  }
*/

#endif
#endif
