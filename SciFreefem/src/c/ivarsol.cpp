#ifndef __IVARSOLVE_CPP
#ifdef INCLUDE_TEMPLATE

#ifdef INCLUDE_TEMPLATE_DEFINITION
# define __IVARSOLVE_CPP
#endif

#else 

// no include in vect.h
// -- no include in vect --
//  so reverve the sens of 
//  after vect
//  if INCLUDE_TEMPLATE_DEFINITION => no indeed of vect.cpp
# include "config.h" 
# ifndef INCLUDE_TEMPLATE_DEFINITION
#  define __IVARSOLVE_CPP
#  include <assert.h>
#  include <math.h>
#  include <stdlib.h>
#  include <iostream>
#  include <fstream>
#  include "ivarsol.h"

using namespace std;
# endif

#include "fem2.h"

#endif
; // empty instruction for metrowerk compiler
#ifdef __IVARSOLVE_CPP

#include "clerror.h"

float addmulop(Expr* e, int k);

extern void ffcopy(int & i, double ** m, float & value);

template <class M, class V, int N> void Isolve<M,V,N>::execute ()
{
  int i, isLU = 0;
  Grid* g = (Grid*)(idmesh->fn);
  Grid* oldActiveMesh = an->activeMesh;

  an->activeMesh = g;
  edp->g = g;
  if(lu) isLU = (int)lu->eval(); // 0 for factorization !=0 if already factorized

  for(i=0; i<N;i++)
    edp->f[i]->resize(g);
  
  int ncoef = edp->g->nt * edp->n * edp->n;
  edp->dif = new float[ncoef];
  edp->dis = new float[ncoef];
  edp->pdx = new float[ncoef];
  edp->pdy = new float[ncoef];
  edp->asym = new float[ncoef];
  edp->pdxy = new float[ncoef];
  edp->pdyx = new float[ncoef];
  {
    for( int j=0; j< ncoef; j++)
      {  edp->dif[j] = 0;
	edp->dis[j] = 0;
	edp->pdx[j] = 0;
	edp->pdy[j] = 0;
	edp->asym[j] = 0;
	edp->pdxy[j] = 0;
	edp->pdyx[j] = 0;
      }}
  
  ncoef = edp->g->ne * edp->n * edp->n;
  edp->rob = new float[ncoef];
  {
    for( int j=0; j< ncoef; j++)
      edp->rob[j] = 0;
  }
  ncoef = edp->g->nv * edp->n ;
  edp->sol = new float[ncoef];
  {
    for( int j=0; j< ncoef; j++)
      edp->sol[j] = 0;
  }
  ncoef = edp->g->ne * edp->n ;
  edp->rhs = new float[ncoef];  
  edp->neuin = new float[ncoef];
  edp->neuout = new float[ncoef];
  {
    for( int j=0; j<ncoef; j++)
      { 	
        edp->rhs[j]    = 0;	
        edp->neuin[j]  = 0;	// rhs of Neumann condition values at quadra pt here and line below
        edp->neuout[j] = 0;
      }
  }
  l->execute(); // fills the coefs and so on
  
  if(idmatrix)
    {
      if(!isLU) 
	if(idmatrix->fn) 
	  bbn = (Profilmatrix<M,V>*)idmatrix->fn;
	else 	bbn = new Profilmatrix<M,V>(g->nv,g->low,g->jlow,g->jhigh);
      else if(idmatrix->fn) 
	bbn = (Profilmatrix<M,V>*)idmatrix->fn;
      else erreur(" Factorized matrix does not exist ");
      idmatrix->fn = bbn;
    }
  else 
    bbn = new Profilmatrix<M,V>(g->nv,g->low,g->jlow,g->jhigh);
  
  if( bbn->low[bbn->size] != g->low[g->nv] ||  bbn->size != g->nv )// grid may have changed at exec time
    {	
      isLU = 0; 
      bbn->destroy();
      bbn = new Profilmatrix<M,V>(g->nv,g->low,g->jlow,g->jhigh);
      if(idmatrix)
	idmatrix->fn = bbn;
    }
  solveprofil<M,V,N>(bbn,edp,isLU);
  //  if (getMatProfil) SaveMatProfil = *bbn;
  
  if(!idmatrix) bbn->destroy();
  
  for(int j=0;j<edp->n;j++)
    for (  i = 0; i < g->nv; i++) 
      edp->f[j]->cc[i] = edp->sol[i*edp->n + j];
  an->activeMesh = oldActiveMesh;
  delete [] edp->dif;	edp->dif=0;
  delete [] edp->dis;	edp->dis=0;
  delete [] edp->rob;	edp->rob=0;
  delete [] edp->pdx;	edp->pdx=0;
  delete [] edp->pdy;	edp->pdy=0;
  delete [] edp->asym;	edp->asym=0;
  delete [] edp->pdxy;	edp->pdyx=0;
  delete [] edp->pdyx;	edp->pdyx=0;
  delete [] edp->sol;	edp->sol=0;
  delete [] edp->rhs;	edp->rhs=0;
  delete [] edp->neuin;	edp->neuin=0;
  delete [] edp->neuout;	edp->neuout=0;
}

template <int N> void Ivarsolve<N>::edpdoit ()
{
  Vector<TN>& bb = *b;
  Grid* oldActiveMesh = an->activeMesh;
  an->activeMesh = (Grid*)idmesh->fn;
  an->gridxyng = an->activeMesh;
  Grid& g= *an->activeMesh;
  Profilmatrix<MN,TN>& aaa = *aa;
  int oldlocal = an->local;
  int i,j,k,n,m;
  float xl[3];
  for(n=0;n<nedp;n++)
    {
      CTab& v = *f2[n];
      CTab& u = *f1[n];
      for(i=0;i<g.nv;i++) { u[i] = 0; v[i] = 0; bb[i][n] = 0;}
    }
  for(k=0; k<g.nt; k++)
    {
      an->trloc = k;
      for(int iloc=0;iloc<3;iloc++)
	{
	  i = g.no(g.t[k].v[iloc]);
	  xl[iloc] = 1;  xl[::next[iloc]] = 0; xl[::next[iloc+1]] = 0;
	  //#if __profile__ //	ProfilerInit(collectSummary,bestTimeBase,10,10);//#endif
	  for(n=0;n<nedp;n++)
	    {
	      CTab& v = *f2[n];
	      v[i] = 1;
	      an->setAn(1, g.t[k].v[iloc]->x,g.t[k].v[iloc]->y,g.t[k].v[iloc]->where,xl,i,iloc,k);
	      l->execute();
	      //#if __profile__ //  ProfilerDump("\pbamg.prof"); ProfilerTerm(); //#undef __profile__ //#endif
	      float bbloc = *id0->storage;
	      bb[i][n] -= bbloc;
	      if(factorize>=0)
		for(m=0;m<nedp;m++)
		  {
		    CTab& u = *f1[m];
		    for(int jloc=0;jloc<3;jloc++)
		      {
			// for(int il=0;il<3;il++) {	int ii=g.no(g.t[k].v[il]); dxu[ii] = g.dxhat[k].v[jloc]; dyu[ii] = g.dyhat[k].v[jloc]; }
			j = g.no(g.t[k].v[jloc]);
			u[j] = 1;
			an->setAn(1, g.t[k].v[iloc]->x,g.t[k].v[iloc]->y,g.t[k].v[iloc]->where,xl,i,iloc,k);
			l->execute();
			float aaloc = *id0->storage;
			MatN<float,N>& amn = aaa(i,j);
			amn(n,m) += aaloc - bbloc;
			u[j] = 0;
			// for(int il=0;il<3;il++) {	int ii=g.no(g.t[k].v[il]); dxu[ii] = 0; dyu[ii] = 0; }
		      } 
		  }
	      v[i] = 0;
	      // for(int il=0;il<3;il++) {	int ii=g.no(g.t[k].v[il]); dxv[ii] = 0; dyv[ii] = 0; }
	    }  	  
	}
    }
  an->activeMesh = oldActiveMesh;
  an->local  = oldlocal;
}


////////////////////////////////////////////
template <int N> void Ivarsolve<N>::execute() 
{
  if(idmesh) an->activeMesh = (Grid*)idmesh->fn;
  Grid& g= *an->activeMesh;
  int  m,n,i;
  for(i=0;i<nedp;i++)
    {
      f1[i]->resize(&g);
      f2[i]->resize(&g);
    }
  factorize = 0;								// nil => 0 create matrix use and destroy
  if(e) 
    {	factorize = (int)e->eval();  	
      if(factorize==0) factorize=1;   			// =0 => 1 create use and keep
      else if(factorize>0) factorize = -1;		// >0 => -1 = reuse and keep
      else if(factorize<0) factorize = -2;		// <0 => -2 reuse and destroy
    }		 
  if(factorize<0)
    {
      if(id0)
	{	
	  if(id0->fn)
	    if(int(*id0->storage)==nedp)
	      aa = (Profilmatrix<MN,TN>*)(id0->fn); 
	    else 
	      throw ErrorExec("Can't re-use a matrix with different size ");
	  else 
	    throw ErrorExec(" Matrix does not exist ");
	}
      else 
	throw ErrorExec("No name to get the matrix from ");
    } 
  else
    {
      aa = new Profilmatrix<MN,TN>(g.nv,g.low,g.jlow,g.jhigh); 
      for(n=0;n<nedp;n++)
	for(m=0;m<nedp;m++)
	  for(i=0;i<aa->csize;i++) 
	    aa->cc[i](m,n) = 0.;
    }
  //  b = &bb;
  // 	if((aa->csize != g.nv)||(aa->jlow != g.jlow)||(aa->jhigh != g.jhigh))
  //  	{	cout <<  "Matrix has been used with a different mesh "<<endl; exit(0);}
  
  b = new Vector<TN>(an->activeMesh->nv);
  edpdoit();
  
  int afactorize = factorize >= 0 ? 1 : 0;
  float pivot = gaussprofil(*aa,*b,afactorize) ;

  if (getMatProfil) 
    {
      Profilmatrix<MN,TN> & a = *aa;
      int ind = 0;
      __gmp = new GetMatrixProfil(a.csize*N*N,a.size,N);
      for (i = 0; i < a.size; i++)
	__gmp -> jlow[i] = (a.jlow)[i];
      for (i = 0; i < a.size; i++)
	__gmp -> jhigh[i] = (a.jhigh)[i];
      for (i = 0; i < a.csize; i++)
	ffcopy(ind,&(__gmp -> Matrix),(a.cc)[i]);
    }
      
  if(afactorize) cout << "\t\t pivot= " << pivot << endl;
  
  for(n=0;n<nedp;n++)
    {	
      CTab& u = *f1[n];
      for(i=0;i<g.nv;i++)
	u[i] = (*b)[i][n];
    } 
  b->destroy();
  if(factorize==0 || factorize==-2) { aa->destroy(); }
  else { id0->fn = aa; *id0->storage = nedp;}
}


#endif
#endif
