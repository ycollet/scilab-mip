/*
  #ifndef  __CONVECT_CPP_
  #ifdef INCLUDE_TEMPLATE
  # ifdef INCLUDE_TEMPLATE_DEFINITION
  #  define __CONVECT_CPP_
  # endif
  #else
  # include "config.h"
  # ifndef INCLUDE_TEMPLATE_DEFINITION
  #  define __CONVECT_CPP_
*/
#  include <assert.h>
#  include <stdlib.h>
//#  define NDEBUG
#  include <cmath>
#  include <iostream>
#  include <fstream>

using namespace std;

#  include "vect.h"
#  include "convec.h"
/*
  # endif
  #endif
  #ifdef __CONVECT_CPP_
*/
const double epsilon = 1e-4;	//- zero for tests on bary coor
const double epsvmin = 1e-10; 	// minimum vlocity to compute characteristics
const double oneeps=1+epsilon;	// 1 for tests on bary coor
static const int next[4] = {1,2,0,1};	

int searchTriangle (Grid& g, const A<float>& u,
		    const A<float>& v, int i1, int& k, int& iloc)
{
  double xu[3];
  int i=0, fit;
	  
  do{ // search in support triangle with best direction
    bTriangle* tk = g.v[i1].supp[i];
    k = g.no(tk); 
    xu[2] = (tk->v[1]->x -tk->v[0]->x)*v[i1] - (tk->v[1]->y - tk->v[0]->y)*u[i1];
    xu[0] = (tk->v[2]->x -tk->v[1]->x)*v[i1] - (tk->v[2]->y - tk->v[1]->y)*u[i1];
    xu[1] = - xu[0] - xu[2];
    iloc = -1;
    do iloc++; while( (iloc<3) && (&g.v[i1]!= tk->v[iloc]));
    if(iloc>=3)
      return 1; // bug in supports of hat functions 
    fit =  xu[iloc] >0 && xu[::next[iloc]] <=1e-30 && xu[::next[iloc+1]] <= 1e-30;
  } while((++i< g.v[i1].nsupp)&&(!fit));
  if(!fit) 
    return 2;	// could not find a good triangle
  else 
    return 0; // OK  	 
}

int xtoX (Grid& g, const A<float>& u, const A<float>& v, double* xl, double *dt, int *k)
// Solves X' = -u, X(0)=x from t=0 to dt
//   x is triangle k and on output X is in new k 
//   ERROR code: <=0: noErr, 

{
  int       iloc,jloc, err=0, kold = *k, count = 0;
  double    mu, xp[3], xu[3], xlold[3], u1k,u2k;
  
  do{
    for(jloc =0; jloc <3;jloc++) //if too near a vertex project on it
      if(xl[jloc] > 1-epsilon) 
	{
	  searchTriangle (g, u,v, g.no(g.t[*k].v[jloc]), kold, iloc);
	  xl[iloc] = 1; xl[::next[iloc]] = 0 ; xl[::next[iloc+1]] = 0;
	  *k = kold; break;
			
	}
				
    int rocked = 0;
  label1:;
    double u1k = 0, u2k = 0;
    for( jloc=0;jloc<3;jloc++) 
      {	int i = g.no(g.t[*k].v[jloc]);
    	u1k += xl[jloc]*u[i];
	u2k += xl[jloc]*v[i];
      } 
    if (fabs(u1k)+fabs(u2k) < epsvmin)
      return -2;		// velocity too small (not an error)
    bTriangle* tk = &g.t[*k];
    double det = 2*tk->area;
    for( jloc=0;jloc<3;jloc++) 
      xu[::next[jloc+1]] = 	 ((tk->v[::next[jloc]]->x - tk->v[jloc]->x)*u2k 
				  - (tk->v[::next[jloc]]->y - tk->v[jloc]->y)*u1k)/det;

    int j=0, fit;
    do{				// get mu & nu with x + mu * u = q^i + nu * (q^{i+}-q^i) 
      int jp = ::next[j], jpp = ::next[jp]; fit =0;
      if(xu[j]>epsvmin && xl[j] >epsilon)
	{	mu = -xl[j] / xu[j]; 
	  xp[j] = 0;
	  xp[jp] = xl[jp]+ mu * xu[jp];
	  xp[jpp] = xl[jpp]+ mu * xu[jpp];
	  fit = (xp[jp] > - epsilon) && (xp[jpp] > - epsilon);
	} 
    } while ((++j<3)&&(!fit));
   
    if(!fit)
      if(!rocked && count++)
	{ // characteristic is near tangent to an edge, use old triangle
	  for(int kloc=0;kloc<3;kloc++) xl[kloc] = xlold[kloc];
	  *k = kold; rocked = 1;
	  goto label1;
	} else
	return 1;		// could not compute a valid  with an edge of this triangle
	  		 	
    if (-mu > *dt)  
      {	mu = -*dt;	*dt = 0; 
	for(int jloc = 0; jloc < 3; jloc++) xl[jloc] += mu * xu[jloc];
	return 0;		// end point found inside the triangle *k
      }

    *dt += mu;
    int kl = (tk->e[--j]->left==&g.t[*k])?  g.no(tk->e[j]->right) 
      : g.no(tk->e[j]->left);
    if(kl<0)
      return -1;	// hit the boundary
    kold = *k;
    *k = kl;
    jloc = -1;
    do jloc++; while(jloc<3 && g.t[kl].v[jloc] != g.t[kold].v[::next[j]]);
    if(jloc == 3)
      return 3; //bug in edge structure
    for(int kloc = 0; kloc <3; kloc++)
      xlold[kloc] = xp[kloc];		
    xl[jloc] = xp[::next[j]];
    if( g.t[kold].v[::next[j+1]] == g.t[kl].v[::next[jloc]])
      {	
	xl[::next[jloc]] = xp[::next[j+1]];
	xl[::next[jloc+1]] = 0;
      } else {	
      xl[::next[jloc+1]] = xp[::next[j+1]];
      xl[::next[jloc]] = 0;
    }
  } while ((*dt > epsvmin) && (count++ <= 50));
  *k = kold;
  return 2;			// trapped in a black hole and playing tennis there 
}


//typedef VectN<float,2> TVN;
//template <class T> 
TVN convect(Grid& g, const Vector<TVN>& f, 
	    const A<float>& u, const A<float>& v, double dt, const int i1)	
// computes foX(i1)  = f(X(q^i1)) where Y'=u(Y), Y(0)= q^i1, X(q^i)= Y(-dt) 
{
  int             fit, k, iloc,err;
  double           dt1, xl[3];

  if(fabs(u[i1])+fabs(v[i1]) < epsvmin)
    return f[i1]; // velocity is too small
  if(err = searchTriangle (g, u, v, i1, k, iloc), err)
    if(err==2 && g.v[i1].where)
      return f[i1];  // inflow part of the boundary
    else 
      return f[i1];	// bug
    
  xl[iloc]= 1; xl[::next[iloc]]=0;xl[::next[iloc+1]]=0;
  //		    rmoveto(x,y);
  dt1 = dt;
  err = xtoX (g,u,v, xl, &dt1, &k); // k is end triangle and xl bary coor of X(x)
  //		    rlineto(x,y);
  return f[g.no(g.t[k].v[0])] * xl[0] 
    + f[g.no(g.t[k].v[1])] * xl[1] 
    + f[g.no(g.t[k].v[2])] * xl[2];
}

//#endif
//#endif
