#include <assert.h>
#include <stdlib.h>

//#define NDEBUG
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

//#include "vect.h"
#include "fem.h"

inline float deter(const float a,const float b,const float c,
		   const float d,const float e,const float f) 
{ float g = (a-c)*(e-f) - (b-c)*(d-f);
  return g;
} 



//-----------------------------------------------------
/* 
   float P1::intgrad2(P1& a, P1& b, P1& c) const
   {// ee = int_vol(a.grad cc . grad cc + b.cc.cc) + int_bdy(c.cc.cc)
   float x[3], y[3], phi[3],gradFx, gradFy, ee=0.;
   for(int k=0;k<g->nt;k++) 
   {
   float a1 = g->t[k].area;
   int n0 = g->no(g->t[k].v[0]), n1 = g->no(g->t[k].v[1]), n2=g->no(g->t[k].v[2]);
   float xa = (a[n0] + a[n1] + a[n2])/3;
   float xb = (b[n0] + b[n1] + b[n2])/3;;
   for(int jloc=0;jloc<3;jloc++) 
   {
   int j = g->no(g->t[k].v[jloc]);
   int jp = g->no(g->t[k].v[(jloc+1)%3]);
   bVertex* vj = g->t[k].v[jloc];
   bVertex* vjp = g->t[k].v[(jloc +1)%3];
   if (vj->where && vjp->where) 
   {	float length = sqrt(pow(vj->x - vjp->x,2) + pow(vj->y - vjp->y,2));
   ee += (c[jp]+c[j])*pow(cc[j]+cc[jp], 2)*length/8;
   }
   x[jloc] = g->v[j].x;
   y[jloc]= g->v[j].y;                                           
   phi[jloc] = cc[j];
   ee += pow(cc[j]+cc[jp], 2) * xa * a1 / 12;
   }
   gradFx = deter(phi[0],phi[1],phi[2],y[0],y[1],y[2]) * 0.5/a1;
   gradFy = deter(x[0],x[1],x[2],phi[0],phi[1],phi[2]) * 0.5/a1;
   ee += xb * a1 * (gradFx*gradFx + gradFy*gradFy);
   }
   return ee;
   }
   /*
   void Laplace::erreur_mid()
   {
   int i;
   float tot=0;
   for (i=0; i< g->nv; i++) 
   tot += fabs(sol[i] - pow(g->v[i].x,2) - pow(g->v[i].y, 2));
   cout <<"\t\t"  << "mean error " << tot / (i+1) << endl;
   }

   void Laplace::erreur_max()
   {
   float errmax = 0;
   float test, vtx, vty;
   for (int i=0; i < g->nv; i++) {
   test = fabs(sol[i] - pow(g->v[i].x,2) - pow(g->v[i].y, 2));
   if ((test > errmax)) 
   errmax = test;
   vtx = g->v[i].x;
   vty = g->v[i].y;
   }
   cout <<"\t\t"  << "maximal absolute error: " << errmax << " en x= " << vtx << " et y= " << vty << endl;
   }

   //-----------------------------------------------------
   void Laplace::solvegradconj(const int niter, const float eps)
   {
   int i,m;
   float g2hm=1, g2h, eps1, ro, gamma, E1, precise=1e-20; 
   P1 hconj(g), u0(g);
   Vector<float> hconjm(g->nv), gradE(g->nv);
   for(i=0;i<g->nv;i++) { hconjm[i] = 0; hconj[i] = 0; u0[i] = sol[i];}
   for(m=0; m<niter; m++) {
   derivener(gradE);
		
   for(i=0;i<g->nv;i++) 
   if(fabs(u0[i])>precise) gradE[i] = 0;
			
   g2h = gradE.scal(gradE);
        	
   if(m==0) { 
   eps1 =  eps*g2h;
   gamma = 0;
   }
   else  gamma = g2h / g2hm;
		
   for(i=0;i<g->nv;i++) 
   hconj[i] = -gradE[i] + gamma * hconjm[i];
    			
   ro = - gradE.scal(hconj) / hconj.intgrad2(dif,vis,rob);
   cout << "iter " <<m<<"\t ro="<<ro<<"\t\t grad= "<<g2h<<endl;
  		
   if(g2h<eps1) 
   break;
  			
   for(i=0;i<g->nv;i++) {
   sol[i] +=  ro* hconj[i];
   hconjm[i] = hconj[i];
   }
   g2hm = g2h;
   }
   hconj.destroy(); u0.destroy(); hconjm.destroy(); gradE.destroy();
   }

   void Laplace::derivener(Vector<float>& r) const
   {
   int i,j,k,jloc;
   float x[3], y[3], phi[3], gradFx,gradFy,gradGx, gradGy;
  
   for(i=0;i<g->nv;i++) r[i]=0; 
  
   for (int l = 0; l < g->nt; l++)
   for(int iloc = 0; iloc < 3 ; iloc++) 
   {
   bVertex* vi = g->t[l].v[iloc];
   bVertex* vj = g->t[l].v[(iloc +1)%3];
   if (vi->where && vj->where) 
   {	int i = g->no(vi), j = g->no(vj);
   float length6 = (sqrt(pow(vi->x - vj->x,2) + pow(vi->y - vj->y,2)))/6;
   r[i] -= (2*neu[i] + neu[j]  - 0.5*(rob[i] +rob[j])*(2*sol[i]+sol[j]))*length6;
   r[j] -= (neu[i] + 2*neu[j]  - 0.5*(rob[i] +rob[j])*(2*sol[j]+sol[i]))*length6;
   }
   }  
   for (k=0; k<g->nt; k++) 
   {
   bTriangle& tk = g->t[k];
   float a1 = 0.5/ tk.area;
   int n0 = g->no(g->t[k].v[0]), n1 = g->no(g->t[k].v[1]), n2 = g->no(g->t[k].v[2]);
   float f0 = rhs[n0] + rhs[n1] + rhs[n2];
   float f1 = sol[n0] + sol[n1] + sol[n2];
   float xvis = (vis[n0] + vis[n1] + vis[n2])/3;
   float xdif = (dif[n0] + dif[n1] + dif[n2])/3;
      
   for (jloc=0; jloc<3; jloc++) 
   {
   j = g->no(tk.v[jloc]);
   r[j] -= ( rhs[j] + f0 - xdif*(sol[j] + f1 )) * (tk.area/12);
   x[jloc] = g->v[j].x;  y[jloc]= g->v[j].y;                                           
   phi[jloc] = sol[j];
   }
   gradFx = a1 * deter(phi[0],phi[1],phi[2],y[0],y[1],y[2]);
   gradFy = a1 * deter(x[0],x[1],x[2],phi[0],phi[1],phi[2]);                 
      
   for(i=0;i<3;i++) 
   {
   gradGy = a1 * deter(x[0],x[1],x[2],i==0,i==1,i==2);
   gradGx = a1 * deter(i==0,i==1,i==2,y[0],y[1],y[2]);
   r[g->no(tk.v[i])] += tk.area * xvis *(gradFx*gradGx + gradFy*gradGy);  
   }
   }
   }

   void Laplace::solveprofil(Profilmatrix<float,float>* bb,int old)
   {
   //	Vector<VectN<float,1>> f(nv),sol(nv);
   //	Bandmatrix<MatN<float,1>,VectN<float,1>> aa(nv,g.bdth);
   int i;
   Grid& t=*g;
   float alpha=0., penal=1e10;
   Profilmatrix<float,float>& aa = *bb;
   Vector<float> aux(sol.size);	
	
   if(!old){
   buildmatlaplace(t,aa, dif, vis,pdx,pdy,asym,pdxy,pdyx);
   intmatgamma(t,aa,rob,-1);
   for(int  i=0;i<t.nv;i++)
   if(fabs(sol[i])>0)
   aa(i,i) = penal; //Dirichlet
   }
   aux = getrhs(t,rhs) + intgamma(t,neu,-1);
   for( i=0;i<t.nv;i++)
   if(fabs(sol[i])>0)
   aux[i] = penal*sol[i];
  
   float pivot = gaussprofil(aa,aux,!old) ;
   float maxsol = aux[0];
   float minsol = aux[0];
   for(  i=0;i<t.nv;i++) { 
   sol[i] = aux[i];
   if(maxsol < sol[i] ) maxsol=sol[i];
   if(minsol > sol[i] ) minsol=sol[i];
   }
   cout<< "min="<< minsol << " max=" << maxsol ;
   if(!old) cout <<"\t\t"  << " pivot=" << pivot;
   cout << endl;
   }

*/
