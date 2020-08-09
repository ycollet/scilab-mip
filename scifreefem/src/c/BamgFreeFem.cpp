// ********** DO NOT REMOVE THIS BANNER **********
//
// SUMMARY: Bamg: Bidimensional Anisotrope Mesh Generator
// RELEASE: 0 
// USAGE  : You may copy freely these files and use it for    
//          teaching or research. These or part of these may   
//          not be sold or used for a commercial purpose with- 
//          out our consent : fax (33) 1 39 63 55 14       
//
// AUTHOR:   F. Hecht,    
// ORG    :  INRIA
// E-MAIL :   Frederic.Hecht@Inria.fr   
//
// ORIG-DATE:     Dec 97

// #define TRACETRIANGLE 3
#undef NDEBUG
extern int verbosity ;
//#define strcasecmp strcmp
#include <string>
#include <cmath>

using namespace std;

#include <stdio.h>
#include <time.h>
#include <string.h>

#include "Meshio.h"
#include "Mesh2.h"
#include "QuadTree.h"
#include "SetOfE4.h"


#include "vect.h"
#include "BamgFreeFem.h"

/*
  void Grid::MakeTriangles()
  {
  cerr << "õ faire / to do" << endl;
  exit(3);
  }
*/
/*
  Int4  i,j,k;
  Int4 nbv = nv;
  Int4 nbt = nt;
  Int4 nbe = 0; 
  assert(!Th);
  assert(!Gh);
  // generation of the Geometry 
  // Construction of the edges 

  Gh = new Geometry(); 
  Th = new Triangles(nbv,&Gh);
  Triangles &T(*Th);
  Geometry &G(*Gh);
  assert(nbt<=T.nbtx);
  assert(nbv<=T.nbvx);
  T.nbt=nbt;
  T.nbv=nbv;
  for (i=0;i<nv;i++)
  {
  T.vertices[i].r.x = v[i].x;
  T.vertices[i].r.y = v[i].y;
  T.vertices[i].ref = v[i].where;
  }
  for (i=0;i<nt;i++)
  for (j=0;j<3;j++)
  {
  T[i]=trianglesTriangle(this, no(t[i][0]),no(t[i][1]),no(t[i][2]));
  T[i].color = t[i].where;
  }
  // generation of the adjacence of the triangles
  Triangle * triangles = T.triangles;
  Vertex * vertices = T.vertices;

  SetOfEdges4 * edge4= new SetOfEdges4(nbt*3,nbv);
  Int4 * st = new Int4[nbt*3];
  for (i=0;i<nbt*3;i++)
  st[i]=-1;
  Int4 kk =0;
  nbe =0;
  Int4 nbei =0;
  for (i=0;i<nv;i++)
  vetices[i].color=0;
  for (i=0;i<nbt;i++)
  for (int j=0;j<3;j++)
  {
  // Int4 i0,i1;
  Int4 i0 =T.Number(triangles[i][VerticesOfTriangularEdge[j][0]]);
  Int4 i1 = T.Number(triangles[i][VerticesOfTriangularEdge[j][1]]);
  Int4 k =edge4->addtrie(i0,i1);,
  Int4 invisible = triangles[i].Hidden(j);
  if(st[k]==-1)
  st[k]=3*i+j;
  else if(st[k]>=0) 
  {
  nbei++;
  Int4 ii =  st[k] / 3;
  int jj = (int) (st[k]%3);
             
  assert( ! triangles[i].TriangleAdj(j) && 
  !triangles[st[k] / 3].TriangleAdj((int) (st[k]%3)));
             
  triangles[i].SetAdj2(j,triangles + ii,jj);
  if (triangles[i].color != triangles[ii].color)
  {
  triangles[i].SetLocked(j);
  vertices[i0].color=1;
  vertices[i1].color=1;
  nbe++;
  }
  if (invisible)  triangles[i].SetHidden(j);
  st[k]=-2-st[k]; 
  }
  else
  {
  cerr << " The edge (" 
  << T.Number(triangles[i][VerticesOfTriangularEdge[j][0]])
  << " , " 
  << T.Number(triangles[i][VerticesOfTriangularEdge[j][1]])
  << " ) is in more than 2 triangles " <<k <<endl;
  cerr << " Edge " << j << " Of Triangle " << i << endl;
  cerr << " Edge " << (-st[k]+2)%3 << " Of Triangle " << (-st[k]+2)/3  << endl;
  cerr << " Edge " << triangles[(-st[k]+2)/3].NuEdgeTriangleAdj((int)((-st[k]+2)%3))
  << " Of Triangle "
  <<  T.Number(triangles[(-st[k]+2)/3].TriangleAdj((int)((-st[k]+2)%3))) << endl;
  MeshError(9999);}  
  }
  // NbOfEdge
  nbe += edge4->NbOfEdges-nbei;

  // --- 
  kk =0;
  for (i=0;i<nbt;i++)
  for (int j=0;j<3;j++)
  {
  Int4 i0 =T.Number(triangles[i][VerticesOfTriangularEdge[j][0]]);
  Int4 i1 = T.Number(triangles[i][VerticesOfTriangularEdge[j][1]]);
  Int4 k =edge4->addtrie(i0,i1);
  Int4 ii =  st[k] / 3;
  if(st[k]>=0 || (i<ii && triangles[i].Locked(j)) )
  {  
  kk++;
  vertices[i0].color=1;
  vertices[i1].color=1;
  }
         
  }
  assert(kk== nbe);
  // compuation of the geom vertices 
  Gh.nbv =0;
  for (i=0;i<nbv;i++)
  if(vertices[i0].color)
  vertices[i0].color= Gh.nbv++;
  else 
  vertices[i0].color=-1;
  Gh.vertices = new GeometricalVertex[Gh.nbv];
  Gh->nbe = nbe;
  Gh->edges = new GeometricalEdges[nbe];
  kk =0;
  Real4 *len = new Real4[ Gh->nbv];

  for (i=0;i<nbv;i++)
  if(vertices[i0].color)
  {
  Gh.vertices[kk].r = vertices[i0];
  Gh.vertices[kk].color=0;
  len[kk]=0;
  kk++;
  }
  assert(kk==Gh.nbv);
  kk =0;
  Real8 Hmin = HUGE_VAL;// the infinie value 
  R2 zero2(0,0);
  for (i=0;i<nbt;i++)
  for (int j=0;j<3;j++)
  {
  Int4 i0 =T.Number(triangles[i][VerticesOfTriangularEdge[j][0]]);
  Int4 i1 = T.Number(triangles[i][VerticesOfTriangularEdge[j][1]]);
  Int4 k =edge4->addtrie(i0,i1);
  Int4 ii =  st[k] / 3;
  if(st[k]>=0 || (i<ii && triangles[i].Locked(j)) )
  {  
  i0 = vertices[i0].color;
  i1 = vertices[i1].color;
  G.edges[kk].v[0]=  G.vertices + i0;
  G.edges[kk].v[1]=  G.vertices + i1;
  G.edges[kk].tg[0]=zero2;
  G.edges[kk].tg[1]=zero2;
  G.edges[kk].SensAdj[0] =  G.edges[i].SensAdj[1] = -1;
  G.edges[kk].Adj[0] =  G.edges[i].Adj[1] = 0;
  G.edges[kk].flag = 0;
  G.vertices[i0].color++;
  G.vertices[i2].color++;
  R2 x12 = G.vertices[i0].r-G.vertices[i1].r;
  Real8 l12=sqrt(x12*x12);
             
  len[i1] += l12;
  len[i2] += l12;
  Hmin = Min(Hmin,l12);
  kk++;                
  }
         
  }
  assert(kk== nbe);
  for (i=0;i<G.nbv;i++) 
  if (G.vertices[i].color > 0) 
  G.vertices[i].m=  Metric(len[i] /(Real4) vertices[i].color);
  else 
  G.vertices[i].m=  Metric(Hmin);
  delete [] len;
  //  find sub domaine --
  for (i=0;i<nbt;i++)
  triangles[i].color=-1;
  Int4 nbsd =0;
  for (i=0;i<nbt;i++)
  if ( triangles[i].color==-1)
  {
  triangles[i].color=nbsd++;
  Int4 k=0;
  st[k++]=i;
  st[k]=0;
  while(1)
  if((j=st[k]++)<3 && !triangles[st[k-1]].Locked(j))
  {
  Triangle * tt = triangles[st[k-1]].TriangleAdj(j);
               
  }
  else 
  {
  k -= 2;
  if (k<0) break;
  Triangle *t = triangles + st[k-1];
  }
  }
  delete edge4;
  delete st;   
  } */

void Grid::th2t(Triangles* tTh)
{ 
  tTh->ReNumberingTheTriangleBySubDomain();
  //tTh->NbRef++;
  Int4  i,j,k=0;
  nv  =  tTh->nbv;
  nt  =   tTh->nbt - tTh->NbOutT;
  v.init(nv);
  t.init(nt);
  // construction of the edges --

  SetOfEdges4 * edge4= new SetOfEdges4(nt*3,nv);
  // 1 compuation of the nb of edges 
  for (i=0;i<tTh->nbe;i++)
    edge4->addtrie(tTh->Number(tTh->edges[i][0]),tTh->Number(tTh->edges[i][1]));
  for (i=0;i<nt;i++)
    for (j=0;j<3;j++)
      {
        Int4 i0 =tTh->Number(tTh->triangles[i][VerticesOfTriangularEdge[j][0]]);
        Int4 i1 = tTh->Number(tTh->triangles[i][VerticesOfTriangularEdge[j][1]]);
        Int4 k =edge4->addtrie(i0,i1);
      }
  ne = edge4->nb();
  delete edge4;
  if(verbosity>1)
    cout <<"\t\t"  << " ne = " << ne << " nv + nt - ne = " <<  nv + nt - ne 
         << " == Nb of hole " <<   endl;
  e.init(ne);
  for (i=0;i<ne;i++)
    {
      e[i].left=e[i].right=0;
      e[i].where=0;
    }
  // now 
  edge4= new SetOfEdges4(ne,nv);
  for (i=0;i<tTh->nbe;i++)
    {
      Int4 k =edge4->addtrie(tTh->Number(tTh->edges[i][0]),tTh->Number(tTh->edges[i][1]));
      e[k].where=tTh->edges[i].ref;
    }

  for (i=0;i<nt;i++)
    for (j=0;j<3;j++)
      {
        Int4 i0 = tTh->Number(tTh->triangles[i][VerticesOfTriangularEdge[j][0]]);
        Int4 i1 = tTh->Number(tTh->triangles[i][VerticesOfTriangularEdge[j][1]]);
        Int4 k  = edge4->addtrie(i0,i1);
        Int4 ii = edge4->i(k);
        Int4 jj = edge4->j(k);
        e[k].in  = v.cc + ii;
        e[k].out = v.cc + jj;
        t[i].e[j]=e.cc+k;
        
        if ( i0 == ii) 
          assert( !e[k].left), e[k].left = t.cc+i;
        else
          assert( !e[k].right), e[k].right= t.cc+i;
      }
  assert( ne == edge4->nb());
  delete edge4;

  for( i=0;i<tTh->nbv;i++)
    {
      v[i].x =  tTh->vertices[i].r.x;
      v[i].y  =  tTh->vertices[i].r.y;
      v[i].where    =  tTh->vertices[i].ref;
    }
  Int4 *reft = new Int4[tTh->nbt];
  Int4 nbref = tTh->ConsRefTriangle(reft);
  for( i=0,k=0;i<tTh->nbt;i++)
    if(tTh->triangles[i].link)
      { 
        t[k].where    = reft[i];  // a faire
        for(int j=0;j<3;j++)
          t[k].v[j]=  &v[tTh->Number(tTh->triangles[i][j])];
        //renu[i]=k;
        assert(k == i);
        k++;
      }
  //else  renu[i]=-1;
  delete [] reft;
  assert ( nt == k);
  tTh->ReMakeTriangleContainingTheVertex();
  // some verification
  for (k=0;k<ne;k++)
    {
      assert( e[k].left ||  e[k].right);
      if (e[k].left)
        {
          bTriangle & tt=  *e[k].left;
          //  cout <<k <<  no(e[k].in) << " " <<  no(e[k].out) << " " 
          //   << no(tt.v[0]) << " " 
          //    << no(tt.v[1]) << " " << no(tt.v[2]) << " " << endl;
          assert((tt.e[0] == &e[k]) ||(tt.e[1] == &e[k]) ||(tt.e[2] == &e[k]) );
          assert((tt.v[0] == e[k].in) ||(tt.v[1] == e[k].in) ||(tt.v[2] == e[k].in));
          assert((tt.v[0] == e[k].out) ||(tt.v[1] == e[k].out) ||(tt.v[2] == e[k].out));
        }
      if (e[k].right)
        {
          bTriangle & tt=  *e[k].right;
          assert((tt.e[0] == &e[k]) ||(tt.e[1] == &e[k]) ||(tt.e[2] == &e[k]) );
          assert((tt.v[0] == e[k].in) ||(tt.v[1] == e[k].in) ||(tt.v[2] == e[k].in));
          assert((tt.v[0] == e[k].out) ||(tt.v[1] == e[k].out) ||(tt.v[2] == e[k].out));
        }
          
    }
  Th = tTh;     
  Gh = &tTh->Gh;
  cout << " th2t " << Th->NbRef<< " " << Gh->NbRef << endl;
}

Geometry * frontiere2Geometry(const frontiere & fr)
{
  Geometry * Gh =  new Geometry;
  if(verbosity>2)
    cout <<"\t\t"  << "  Begin: ConstGeometry from frontiere"  << endl;
  //assert(empty());
  const char * filename = "FREEFEM.gh";
  Gh->name=new char [strlen(filename)+1];
  strcpy(Gh->name,filename);
  Real8 Hmin = HUGE_VAL;// the infinie value 
  Int4 hvertices =0;
  Int4 i;
  Int4 Version,dim=0;
  Gh->MaximalAngleOfCorner =30.00*Pi/180.0;
  Gh->nbv = fr.nbp;
  Gh->nbvx = Gh->nbv;
          
  Gh->vertices = new GeometricalVertex[Gh->nbvx];
  assert(Gh->nbvx >= Gh->nbv);
  Gh->nbiv = Gh->nbv;
  Int4 k=0;
          
  for (i=0;i<Gh->nbv;i++) 
    {
      Gh->vertices[i].r.x = fr.xy[k++];
      Gh->vertices[i].r.y = fr.xy[k++];
      Gh->vertices[i].ref = fr.ng[i];
      Gh->vertices[i].color =0;
      Gh->vertices[i].Set();
      //  vertices[i].SetCorner();
      //  vertices[i].SetRequired();
    }

  Gh->pmin =  Gh->vertices[0].r;
  Gh->pmax =  Gh->vertices[0].r;
  // recherche des extrema des vertices pmin,pmax
  for (i=0;i<Gh->nbv;i++) 
    {
      Gh->pmin.x = Min(Gh->pmin.x,Gh->vertices[i].r.x);
      Gh->pmin.y = Min(Gh->pmin.y,Gh->vertices[i].r.y);
      Gh->pmax.x = Max(Gh->pmax.x,Gh->vertices[i].r.x);
      Gh->pmax.y = Max(Gh->pmax.y,Gh->vertices[i].r.y);
    }
  Gh->coefIcoor= (MaxICoor)/(Max(Gh->pmax.x-Gh->pmin.x,Gh->pmax.y-Gh->pmin.y));
  assert(Gh->coefIcoor >0);
  if (verbosity>2) 
    {
      cout <<"\t\t"  << "     Geom: min="<< Gh->pmin << "max ="<< Gh->pmax 
           << " hmin = " << Gh->MinimalHmin() <<  endl;
    }
       
  R2 zero2(0,0);
  Gh->nbe = fr.nbs;
  Gh->edges = new GeometricalEdge[Gh->nbe];
  if(verbosity>5) 
    cout <<"\t\t"  << "     Record Edges: Nb of Edge " << Gh->nbe <<endl;
  assert(Gh->edges);
  assert (Gh->nbv >0); 
  Real4 *len =0;
  if (!hvertices) 
    {
      len = new Real4[Gh->nbv];
      for(i=0;i<Gh->nbv;i++)
        len[i]=0;
    }
  k=0;
  for (i=0;i<Gh->nbe;i++)               
    {
      Int4 i1 =  fr.s[k++], i2 =  fr.s[k++];
      Gh->edges[i].ref = fr.ngf[i];
      Gh->edges[i].v[0]=  Gh->vertices + i1;
      Gh->edges[i].v[1]=  Gh->vertices + i2;
      R2 x12 = Gh->vertices[i2].r-Gh->vertices[i1].r;
      Real8 l12=sqrt(x12*x12);
      Gh->edges[i].tg[0]=zero2;
      Gh->edges[i].tg[1]=zero2;
      Gh->edges[i].SensAdj[0] = Gh->edges[i].SensAdj[1] = -1;
      Gh->edges[i].Adj[0] = Gh->edges[i].Adj[1] = 0;
      Gh->edges[i].flag = 0;
      if (!hvertices) 
        {
          Gh->vertices[i1].color++;
          Gh->vertices[i2].color++;
          len[i1] += l12;
          len[i2] += l12;
        }
                  
      Hmin = Min(Hmin,l12);
    }
  // definition  the default of the given mesh size 
  if (!hvertices) 
    {
      for (i=0;i<Gh->nbv;i++) 
        if (Gh->vertices[i].color > 0) 
          Gh->vertices[i].m=  Metric(len[i] /(Real4) Gh->vertices[i].color);
        else 
          Gh->vertices[i].m=  Metric(Hmin);
      delete [] len;
                  
      if(verbosity>3) 
        cout <<"\t\t"  << "     Geom Hmin " << Hmin << endl;
    }

  Gh->NbSubDomains=fr.nbsd;
  if (Gh->NbSubDomains>0) 
    {
      Gh->subdomains = new GeometricalSubDomain[  Gh->NbSubDomains];
      Int4 i0,i1,i2,i3;
      for (i=0;i<Gh->NbSubDomains;i++) 
        {
          i1 = fr.sd[2*i];
          Gh->subdomains[i].sens = 1;
          if(i1<0){i1=-i1; Gh->subdomains[i].sens = -1;}
          assert(i1<Gh->nbe && i1>=0);
          Gh->subdomains[i].edge=Gh->edges + (i1);
          Gh->subdomains[i].ref = i+1;
        }
    }
  Gh->AfterRead();   
  return Gh;  
}

Int4 FindTriangle(Triangles &Th, Real8 x, Real8 y, double* a)
{
  CurrentTh=&Th;
  assert(&Th);
  I2 I = Th.toI2(R2(Min(Max(Th.pmin.x,x),Th.pmax.x),Min(Max(Th.pmin.y,y),Th.pmax.y))); 
  Icoor2 dete[3];
  Triangle & tb = *Th.FindTriangleContening(I,dete);
   
  if  (tb.link) 
    { // internal point in a true triangles
      a[0]= (Real8) dete[0]/ tb.det;
      a[1]= (Real8) dete[1] / tb.det;
      a[2] = (Real8) dete[2] / tb.det;   
      return Th.Number(tb);
    }
  else return -1;
}
Triangles::Triangles(const Triangles & Tho,const int *flag ,const int *bb)
  : Gh(*(new Geometry())), BTh(*this)
{ // truncature
  // 
  
  char cname[] = "trunc";

  int i,j,k,itadj;
  int kt=0;
  int * kk    = new int [Tho.nbv];
  Int4 * reft = new Int4[Tho.nbt];
  Int4 nbInT =    Tho.ConsRefTriangle(reft);
  Int4 * refv = new Int4[Tho.nbv];

  for (i=0;i<Tho.nbv;i++)
    kk[i]=-1;
  for (i=0;i<Tho.nbv;i++)
    refv[i]=0;
  int nbNewBedge =0;
  int nbOldBedge =0;  
  for (i=0;i<Tho.nbt;i++)
    if(  reft[i] >=0 && flag[i]) 
      {
        const Triangle & t = Tho.triangles[i];
        kt++;
        kk[Tho.Number(t[0])]=1;
        kk[Tho.Number(t[1])]=1;
        kk[Tho.Number(t[2])]=1;
        itadj=Tho.Number(t.TriangleAdj(0));
        if (  reft[itadj] >=0 && !flag[itadj])
          { nbNewBedge++;
	    refv[Tho.Number(t[VerticesOfTriangularEdge[0][0]])]=bb[i];
	    refv[Tho.Number(t[VerticesOfTriangularEdge[0][1]])]=bb[i];
          }
        itadj=Tho.Number(t.TriangleAdj(1));
        if (  reft[itadj] >=0 && !flag[itadj])
          { nbNewBedge++;
	    refv[Tho.Number(t[VerticesOfTriangularEdge[1][0]])]=bb[i];
	    refv[Tho.Number(t[VerticesOfTriangularEdge[1][1]])]=bb[i];}
        itadj=Tho.Number(t.TriangleAdj(2));
        if (  reft[itadj] >=0 && !flag[itadj])
          { nbNewBedge++;
	    refv[Tho.Number(t[VerticesOfTriangularEdge[2][0]])]=bb[i];
	    refv[Tho.Number(t[VerticesOfTriangularEdge[2][1]])]=bb[i];}
      }
  k=0;
  for (i=0;i<Tho.nbv;i++)
    if (kk[i]>=0) 
      kk[i]=k++;
  cout << " number of vertices " << k << " remove = " << Tho.nbv - k << endl;
  cout << " number of triangles " << kt << " remove = " << nbInT-kt << endl;
  cout << " number of New boundary edge " << nbNewBedge << endl;
  Int4 inbvx =k;
  PreInit(inbvx,cname);
  for (i=0;i<Tho.nbv;i++)
    if (kk[i]>=0) 
      {
        vertices[nbv] = Tho.vertices[i];
        if (!vertices[nbv].ref)
          vertices[nbv].ref = refv[i];
        nbv++;
      }
  assert(inbvx == nbv);
  for (i=0;i<Tho.nbt;i++)
    if(  reft[i] >=0 && flag[i]) 
      {
        const Triangle & t = Tho.triangles[i];
        int i0 = Tho.Number(t[0]);
        int i1 = Tho.Number(t[1]);
        int i2 = Tho.Number(t[2]);
        assert(i0>=0 && i1 >= 0 && i2  >= 0);
        assert(i0<Tho.nbv && i1 <Tho.nbv && i2  <Tho.nbv);
        // cout <<i<< " F" <<  flag[i] << " T " << nbt << "   = " <<  kk[i0] << " " << kk[i1] << " " << kk[i2] ;
        // cout << " OT  " <<  i0 << " "  << i1 << " " << i2  << " " << reft[i] << endl;
        triangles[nbt] = Triangle(this,kk[i0],kk[i1],kk[i2]);
        triangles[nbt].color = Tho.subdomains[reft[i]].ref; 
        nbt++;           
      }
  assert(kt==nbt);
  if (nbt ==0 && nbv ==0) {
    cout << "Error all triangles was remove " << endl;
    MeshError(999);
  }
  delete [] kk;
  delete [] reft;
  delete [] refv;
  double cutoffradian = 10.0/180.0*Pi;
  ConsGeometry(cutoffradian);
  Gh.AfterRead(); 
  SetIntCoor();
  FillHoleInMesh();
   
  assert(NbSubDomains);
  assert(subdomains[0].head && subdomains[0].head->link);
             
}

