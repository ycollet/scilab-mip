#include <assert.h>
#include <stdlib.h>

//#define NDEBUG
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

#include "vect.h"
#include "Mesh2.h"
#define abss(a)(a >= 0 ? a : -(a))

float norm2(float a) {return fabs(a);}
void inspec()
{
  float a=0;
  a =1.0;
}

//---------------------------------------------------------
void bVertex::fill(float xx, float yy, int ngg) 
{ x = xx; y = yy; where = ngg; nsupp = 0; nmate=0; supp.cc=0; mate.cc=0; }

//---------------------------------------------------------
void bTriangle::fill(bVertex& p0, bVertex& p1, bVertex& p2, int wwhere) 
{
  v[0] = &p0; v[1] = &p1; v[2] = &p2; where = wwhere; 
}

void bEdge::fill(bVertex* iin, bVertex* oout, 
		 bTriangle*  lleft, bTriangle* rright)
{ 
  in = iin; out=oout; left=lleft; right=rright;
}

/*//---------------------------------------------------------
  Grid::Grid(const char *path ):v(),t(),e()
  {// reads a triangulation in MacGfem format
  int i,where, i0,i1,i2;
  float x,y;
  ifstream file(path);

  assert(!file.fail());
  assert(file >> nv >> nt);
  nbholes =100;  // not in macgfem format
  v.init(nv);
  t.init(nt);
  e.init(nv+nt+nbholes-1);
  for( i=0; i<nv; i++ ) 
  {
  assert(file >> x >> y >> where);
  v[i].fill(x,y,where);
  }
  for( i=0; i<nt; i++ ) 
  {
  assert(file >> i0 >> i1 >> i2 >> t[i].where);
  t[i].fill(v[i0-1],v[i1-1],v[i2-1],where); // fortran convention requires -1
  }
  }

  //---------------------------------------------------------
  void Grid::save(const char *path )
  {// write a triangulation in MacGfem format
  int i,j=0;
  ofstream file(path);

  file << nv <<"	"<< nt<<endl;
  for( i=0; i<nv; i++ ) 
  file << v[i].x <<"	"<< v[i].y <<"	"<< v[i].where<<endl;
  for( i=0; i<nt; i++ ) 
  file << no(t[i].v[0])+1 <<"	"<< no(t[i].v[1])+1 <<"	"<< no(t[i].v[2])+1<<"	"<< j<<endl;
  file.close();
  }
*/

//---------------------------------------------------------
void Grid::computegeom(int dontTouchEdges)
{ 
  for (int i = 0; i < nt; i++) 
    {
      bTriangle& ti = t[i];
      t[i].area = ((ti.v[1]->x - ti.v[0]->x) * (ti.v[2]->y - ti.v[0]->y)
		   - (ti.v[2]->x - ti.v[0]->x) * (ti.v[1]->y - ti.v[0]->y))/2;
    }
  for (int k = 0; k < ne; k++) 
    {
      bEdge& ek = e[k];
      float aux = (ek.out->x - ek.in->x) * (ek.out->x - ek.in->x)
	+(ek.out->y - ek.in->y) * (ek.out->y - ek.in->y);
      e[k].length = (float)sqrt(aux);
    }
}

//---------------------------------------------------------
void Grid::prepgrid(int dontTouchEdges) 
{ 
  bEdge *e1;
  int neold = ne;

  getbdth(); 
  fillvsupp(); 	
  getnmate();	

  if(!dontTouchEdges)
    {	e1 =  new bEdge[ne];
      for(int i=0; i< neold; i++)
	{
	  e1[i].in = e[i].in;
	  e1[i].out = e[i].out;
	  e1[i].where = e[i].where;
	}
    }
	
  fillmate(dontTouchEdges);
  computegeom(dontTouchEdges);
	
  norml.destroy();
  norml.init(nv);  // preparing to compute the normal
  int i;
  for( i=0;i<nv;i++)
    {	norml[i].x = 0; norml[i].y = 0;}
  for( i=0; i<ne; i++)
    {	bEdge& ei = e[i];
      bVertex& vin = *ei.in;
      bVertex& vout = *ei.out;
      if((vin.where !=0)&&(vout.where !=0))
	{			// compute the normal
	  float t1 = vout.x - vin.x, t2 = vout.y - vin.y;
	  float nnn = sqrt(t1*t1+t2*t2)/2;
	  if(fabs(nnn)<1e-30)
	    cout <<" Boundary points are too close: can't get the normal. "<< endl;
	  norml[no(&vin)].x -= t2/nnn;
	  norml[no(&vin)].y += t1/nnn;
	  norml[no(&vout)].x -= t2/nnn;
	  norml[no(&vout)].y += t1/nnn;
	}		 
    }
		
  if(!dontTouchEdges)
    {
      for( i=0; i<ne; i++)
	{	bEdge& ei = e[i];
	  bVertex& vin = *ei.in;
	  bVertex& vout = *ei.out;
	  ei.where = 0;
	  if((vin.where !=0)&&(vout.where !=0))
	    {
	      if( vin.where == vout.where )
		ei.where =  vin.where;
	      else for(int j = 0; j< neold; j++)
		     if(((&vin == e1[j].in ) && (&vout == e1[j].out))||
			((&vin == e1[j].out ) && (&vout == e1[j].in))) ei.where = e1[j].where;
	    }
	}
      delete [] e1;
    }
  getprofil();
  derivhat();
  // compute the normal
}


void Grid::derivhat()
{
  dxhat.resize(nt);
  dyhat.resize(nt);
  float u[3];

  { for(int k=0; k<3;k++) u[k] = 0;}
     
  for(int k=0; k<nt;k++)
    {
      for(int i0=0; i0<3;i0++)
    	{
	  int i1 = (i0+1)%3;
	  int i2 = (i1+1)%3;
	  u[i0]=1;
	  dxhat[k].v[i0] = ((u[i1] - u[i0]) * (t[k].v[i2]->y - t[k].v[i0]->y)
			    - (u[i2] - u[i0]) * (t[k].v[i1]->y - t[k].v[i0]->y))/(2*t[k].area);
	  dyhat[k].v[i0] = ((t[k].v[i1]->x - t[k].v[i0]->x) * (u[i2] - u[i0])
			    - (t[k].v[i2]->x - t[k].v[i0]->x) * (u[i1] - u[i0]))/(2*t[k].area);
	  u[i0] = 0;
        }
    }
}
//---------------------------------------------------------
void Grid::square(int nx, int ny) 
{ // triangulate the unit square and assign bdy id = 4 (left), 3 (top), 2 (right), 1 (bottom)
  int i,j;
  nbholes =0;
  v.init(nv = (nx + 1) * (ny + 1));
  t.init(nt = 2 * nx * ny);
  e.init(nv+nt+nbholes-1);
  float dx = 1.F/nx, dy = 1.F/ny;
  for (i = 0; i <= nx; i++)
    for (j = 0; j <= ny; j++)
      {
	int where = 4*(i==0) + 2*(i==nx);
	if(j==ny) where = 3;
	else if(j==0 ) where = 1;
	v[i + j * (nx + 1)].fill(i * dx, j * dy, where);
      }
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      {
	int k = 2 * (i + j * nx);
	int l = i + j * (nx + 1);
	t[k].fill(v[l], v[l + 1], v[l + nx + 1],0); 
	t[k + 1].fill(v[l + 1], v[l + nx + 2], v[l + nx + 1],0); }
}

//---------------------------------------------------------
void Grid::getbdth() 
{
  int k,iloc, b1;
  bdth = 0;
  for(k=0;k<nt;k++)
    for(iloc=0;iloc<3;iloc++) 
      if(b1 = abss(no(t[k].v[iloc]) 
		   - no(t[k].v[(iloc==2?0:iloc+1)])), b1 > bdth) 
	bdth = b1;
}

//---------------------------------------------------------
void Grid::getprofil() 
{
  int i, k;
  low.resize(nv+1), jlow.resize(nv), jhigh.resize(nv);
  low[0] = 0;
  for(i=0; i<nv;i++)
    { 	jlow[i] = i; jhigh[i] = i;
      for(int j=0; j<v[i].nmate;j++) 
	{	k = no(v[i].mate[j]);
	  if (k < jlow[i]) jlow[i]=k;
	  else if (k > jhigh[i]) jhigh[i]=k;
	}
    }
  for(i=0; i<nv;i++)
    for(k=jlow[i];k<i;k++) if(jhigh[k]<i) jhigh[k]=i;
  for(i=0; i<nv;i++)
    low[i+1] = low[i] + jhigh[i] - jlow[i] +1;
}

//---------------------------------------------------------
//---------------------------------------------------------
void Grid::getnmate() const
{
  int k,iloc, ilocp, ilocpp, ngp, ngpp;
  for(k=0;k<nv;k++) v[k].nmate =0;
  for( k=0;k<nt;k++)
    for( iloc=0; iloc<3; iloc++) 
      {
	ilocp = (iloc+1)%3;
	ilocpp = (ilocp+1)%3;
	if(t[k].v[iloc]->where)
	  {
	    ngp = t[k].v[ilocp]->where != 0;
	    ngpp = t[k].v[ilocpp]->where != 0;
	    if(ngp+ngpp != 2) t[k].v[iloc]->nmate += 2 + ngp + ngpp;
	  } else t[k].v[iloc]->nmate += 2;
      }
  for(k=0;k<nv;k++)
    v[k].nmate /= 2;
  for( k=0;k<nt;k++)  // triangles with 3 bdy points
    for( iloc=0; iloc<3; iloc++) 
      if(t[k].v[0]->where && t[k].v[1]->where && t[k].v[2]->where)
	{	if(!t[k].v[iloc]->nmate) t[k].v[iloc]->nmate = 2;
	  else t[k].v[iloc]->nmate += 1;
	}
}

//---------------------------------------------------------
void Grid::fillmate(int dontTouchEdges) 
{
  int i,j,k,iloc, ilocp, ilocpp;
  bTriangle* tk;
  bVertex* vip;
	
  for(i=0;i<nv;i++) 
    {
      v[i].mate.resize(v[i].nmate+1); 
      for(j=0;j<v[i].nmate;j++) v[i].mate[j] = 0;
      v[i].nmate = 0;
    }
	
  if(!dontTouchEdges) 	ne = 0; 
	
  for( k=0;k<nt;k++)
    for( iloc=0; iloc<3; iloc++) 
      {
	ilocp = (iloc==2)? 0 : iloc+1;
	i = no(t[k].v[iloc]);
	vip = t[k].v[ilocp];
	for(j=0;j<v[i].nmate;j++) if ( v[i].mate[j] == vip) break;
	if(j==v[i].nmate)
	  {
	    v[i].mate[v[i].nmate++] = vip;
	    if(!dontTouchEdges) 
	      {	e[ne].fill(&v[i], vip, &t[k], 0);
		t[k].e[(iloc==0)?2:iloc-1] = &e[ne++];	
	      }
	  }
	ilocpp = (ilocp==2)? 0 : ilocp+1;
	vip = t[k].v[ilocpp];
	for(j=0;j<v[i].nmate;j++) if ( v[i].mate[j] == vip) break;
	if(j==v[i].nmate)	
	  v[i].mate[v[i].nmate++] = vip;
      }
  if(!dontTouchEdges) 
    for(j=0;j<ne;j++) // build e[j].right 
      {
	i = no(e[j].in);
	for( int kloc =0;kloc<v[i].nsupp;kloc++)
	  {
	    tk = v[i].supp[kloc]; 
	    if(tk != e[j].left)
	      for(iloc=0;iloc<3;iloc++) 
		if( tk->v[iloc] == e[j].out)
		  {
		    e[j].right = tk;
		    tk->e[(iloc==0)?2:iloc-1] = &e[j]; //1-12
		    goto l3;
		  }
	  }
      l3: ;
      }
}

//---------------------------------------------------------
void Grid::fillvsupp()
{ // fills v[].nsupp then v[].supp[]
  int i,k,iloc;
  for(i=0;i<nv;i++)  v[i].nsupp = 0;
  for( k=0;k<nt;k++) 
    for( iloc=0; iloc<3; iloc++) 
      v[no(t[k].v[iloc])].nsupp++;
  for(i=0;i<nv;i++)  { v[i].supp.resize(v[i].nsupp+1); v[i].nsupp = 0;}
  for( k=0;k<nt;k++) 
    for( iloc=0; iloc<3; iloc++) 
      {
	i = no(t[k].v[iloc]);
	v[i].supp[v[i].nsupp++] = &t[k];
      }
}

void Grid::readgrid(const char *path )
{// reads a triangulation in MacGfem format
  int i,j, ii;
  ifstream file(path);
  
  file >> nv;
  file >> nt;
  nbholes =nbholesmax;  // not in macgfem format
  //  cout <<"\t\t"  << "reading mesh with nv = " << nv << "nt = "<< nt << endl;
  v.init(nv);
  t.init(nt);
  e.init(ne=nv+nt+nbholes-1);
  for( i=0; i<nv; i++ ) 
    {	
      file >> v[i].x >> v[i].y;
      file >> v[i].where;
    }
  for( i=0; i<nt; i++ ) 
    {  
      for(j=0;j<3;j++){ file >> ii; t[i].v[j] = &v[ii-1];}
      file >> t[i].where;
    }
  for( i=0; i<nt; i++ ) 
    t[i].area = ((t[i].v[1]->x - t[i].v[0]->x) 
		 * (t[i].v[2]->y - t[i].v[0]->y)
		 - (t[i].v[2]->x - t[i].v[0]->x) 
		 * (t[i].v[1]->y - t[i].v[0]->y))/2;
}

//---------------------------------------------------------
void Grid::save(const char *path, int debugformat ) const
{// write a triangulation in MacGfem format
  int i,j=0;
  
  if(debugformat) {   dump(path); return;}
  
  ofstream file(path);
  
  file << nv ; file<<"	"<< nt<<endl;
  for( i=0; i<nv; i++ ) 
    {
      file << v[i].x <<"	"<< v[i].y <<"	";
      file << v[i].where<<endl;
    }
  for( i=0; i<nt; i++ ) 
    {	
      for(int k=0;k<3;k++)
	file << no(t[i].v[k])+1 <<"	";
      file << j<<endl;
    }
}

//----------------------------------------------------------
void Grid::gnusave(const char* path) const
{
  int j=0;
  ofstream file(path);
  assert(!file.fail());
  for (int k=0; k<nt; k++)
    file << v[j=no(t[k].v[0])].x << "\t" << v[j].y << endl
	 << v[j=no(t[k].v[1])].x << "\t" << v[j].y << endl
	 << v[j=no(t[k].v[2])].x << "\t" << v[j].y << endl
	 << v[j=no(t[k].v[0])].x << "\t" << v[j].y << endl
	 << endl;
}

void Grid::dump(const char* path) const
{
  ofstream file(path);
  
  file << "Nb of vertices "<< nv << "	Nb of Trianges " << nt <<"	Nb of Edges " << ne <<endl;
  int i,j;
  for(i=0; i<nv;i++)
    {
      file << "Vertex " << i << "	nsupp =" << v[i].nsupp << "	nmate = " <<v[i].nmate
	   <<"	x= " << v[i].x <<"	y= "<< v[i].y << "	where= " << v[i].where<< endl;
      for( j=0; j < v[i].nsupp; j++)
	file << j <<"		supp is " << no(v[i].supp[j]) << endl;
      for( j=0; j < v[i].nmate; j++)
	file << j <<"		mate is " << no(v[i].mate[j]) << endl;
    }
  file << endl;
  for( i=0; i<ne;i++)
    {
      file << "Edge " << i << "	in =" << no(e[i].in) << "	out = " <<no(e[i].out) << 
	"	left =" << no(e[i].left) << "	right = " <<no(e[i].right) <<endl;
    }
  file << endl;
  for( i=0; i<nt;i++)
    {
      file << "Triangle " << i << "	v[0] = " << no(t[i].v[0]) << "	v[1] = " <<no(t[i].v[1]) << 
	"	v[2] = " << no(t[i].v[2]) <<endl;
      file <<"	\t e[0]=" << no(t[i].e[0]) << "	e[1] = " <<no(t[i].e[1]) << 
	"	e[2] = " << no(t[i].e[2]) << "	where= " << t[i].where<< endl;
    }
}
//----------------------------------------------------------
void Grid::destroy()
{ 
  cout <<"\t\t"  << "Grid::destroy()" << this << Th->identity 
       << " Nbref Th " <<  Th->NbRef << " Nbref Gh " << Gh->NbRef << endl;
  v.destroy(); 
  t.destroy(); 
  e.destroy();
  low.destroy(); 
  jhigh.destroy(); 
  jlow.destroy();
  if (Th && (Th->NbRef)-- ==0)
    delete Th;
  if (Gh && (Gh->NbRef)-- ==0)
    delete Gh;
}
	
//---------------------------------------------------------
void Grid::reGrid(const Grid* g) 
{ 
  int i;
  nv=g->nv;
  nt=g->nt;
  ne = g->ne;
  v.resize(g->nv);
  t.resize(g->nt);
  e.resize(g->ne);
  Th = g->Th; 
  Gh = g->Gh; 
  NumThinGrid = g->NumThinGrid;
  for( i=0; i<nt; i++) 
    for(int j=0;j<3;j++) 
      t[i].v[j] = &v[g->no(g->t[i].v[j])] ;
  v = g->v;
  for( i=0; i<ne; i++) 
    e[i] = g->e[i];
}

//--------------------------------------------------------
void Grid::copyev(const Grid* g) // used only to recopy grids within constructor
{    int i,j;
  for( i=0; i<nt; i++) 
    for( j=0;j<3;j++) 
      {  t[i].v[j] = &v[g->no(g->t[i].v[j])] ;
	t[i].e[j] = &e[g->no(g->t[i].e[j])] ;
      }
  for (i=0;i<ne;i++)
    {
      e[i].in = &v[g->no(g->e[i].in)] ;
      e[i].out = &v[g->no(g->e[i].out)] ;
      if (e[i].left)  e[i].left = &t[g->no(g->e[i].left)] ;
      if (e[i].right) e[i].right = &t[g->no(g->e[i].right)] ;

    }    
  if(g->NumThinGrid)
    {	NumThinGrid = new int[g->nv];
      for(i=0; i<nv;i++) NumThinGrid[i]=g->NumThinGrid[i];    
    } else NumThinGrid = 0;
	    
}

//--------------------------------------------------------
void Grid::initquad(Grid* t)
{	
  int k1,k, kt,inside;
  double a[3]; 

  if(gridFriend==t) return;
  quad = new int*[nt];
  gridFriend = t;

  int* sizequad = new int[nt];
  for(k=0; k<nt;k++) sizequad[k] =0;

	
  for(kt=0; kt<t->nt;kt++)
    for(int iloc=0;iloc<3;iloc++)
      {
    	int i = t->no(t->t[kt].v[iloc]);
    	int inext = t->no(t->t[kt].v[(iloc+1)%3]);
    	int inextnext =  t->no(t->t[kt].v[(iloc+2)%3]);// slightly inside to handle discontinuous functions
    	float x = (t->v[i].x+t->v[inext].x + 0.001*t->v[inextnext].x)/2.001;
    	float y = (t->v[i].y+t->v[inext].y + 0.001*t->v[inextnext].y)/2.001;
	k1=(int)FindTriangle(*Th,x,y,a,inside);
	sizequad[k1]+=1; // this quadrature point will be processed in triangle kt
      }

  for(int  kkk=0; kkk<nt;kkk++) 
    { 		quad[kkk] = new int[1+sizequad[kkk]]; 
      quad[kkk][0] = 0;	// will contain the size of quad[k]
    }
	 
  for(kt=0; kt<t->nt;kt++)
    for(int iloc=0;iloc<3;iloc++)
      {
    	int i = t->no(t->t[kt].v[iloc]);
    	int inext = t->no(t->t[kt].v[(iloc+1)%3]);
    	int inextnext =  t->no(t->t[kt].v[(iloc+2)%3]);// slightly inside to handle discontinuous functions
    	float x = (t->v[i].x+t->v[inext].x + 0.001*t->v[inextnext].x)/2.001;
    	float y = (t->v[i].y+t->v[inext].y + 0.001*t->v[inextnext].y)/2.001;
	k1=(int)FindTriangle(*Th,x,y,a,inside);
	int kk1 = quad[k1][0]+1;
	quad[k1][0] = kk1;
	quad[k1][kk1] = 3*kt+iloc; // this quadrature point will be processed in triangle kt
	k1+=1;
      }
  /*     for( kt=0; kt<nt; kt++)
	 {	for(int j=1;j<=quad[kt][0];j++) cout <<"\t\t"  << quad[kt][j] <<" ";
	 cout <<"\t\t"  << endl;
	 }
	 exit(0);
  */      
}

//---------------------------------------------------------
void EFSpace::save(const char* path) const
{  // writes a P1 continuous function in freefem format
  ofstream file(path);
  int nv = size;
  assert(!file.fail());
  file << nv << endl;
  for(int i=0; i<nv; i++) 
    file << cc[i] << endl;
}

//-----------------------------------------------------------
void P1::gnusave(const char* path) const
{ // writes a P1 continuous function in GNUPlot format 
  ofstream file(path);
  assert(!file.fail());
  
  for (int ln=0; ln<20; ln++) { // 10 lignes de niveaux
    float lambda = min()+ln*(max()-min())/20;
    bVertex p[3];
    
    for (int k=0; k<g->nt; k++) {
      int l = 0;
      for (int n=0; n<3; n++) {
	int i, j;
	i = g->no(g->t[k].v[n]);
	j = g->no(g->t[k].v[(n+1)%3]);
	if (cc[i] != cc[j]) {
	  float mu;
	  mu = (lambda - cc[i]) / (cc[j] - cc[i]);
	  if ((mu >= 0) && (mu <= 1)) {
	    p[l].x = mu * g->v[j].x + (1 - mu) * g->v[i].x;
	    p[l++].y = mu * g->v[j].y + (1 - mu) * g->v[i].y;
	  } //if
	} // if
      } //for n
      
      if (l >= 2)
	file << p[0].x << "\t" << p[0].y << endl
	     <<	p[1].x << "\t" << p[1].y << endl
	     << endl;
    } //for k
  } //for m
}	
//-----------------------------------------------------------
void P0::gnusave(const char* path) const
{ 
  cout << " ################################################# " << endl;
  cout << " Pas de gnu save P0 d»sole sur " << path <<  endl;
  cout << " ################################################# " << endl;
}	

//---------------------------------------------------------
void EFSpace::load(const char* path)
{  // read a P1 continuous function in freefem format
  ifstream file(path);
  int nv ;
  assert(!file.fail());
  file >> nv;
  assert(nv==size);
  for(int i=0; i<nv; i++) file >> cc[i];
}

