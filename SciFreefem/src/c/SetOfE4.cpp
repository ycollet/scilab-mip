#include <iostream> 

using namespace std;

#include "meshtype.h"
#include "SetOfE4.h"


SetOfEdges4::SetOfEdges4(Int4 mmx,Int4 nnx)
{nx=nnx;
  nbax=mmx;
  NbOfEdges = 0;
  tete= new Int4 [nx];
  Int4 i=nx;
  while (i--)
    tete[i] = -1;// vide 
  Edges =new Int4Edge[nbax];
}
    
Int4 SetOfEdges4::find(Int4 ii,Int4 jj)
{ 
  if (tete == 0 ) {
    cerr <<"SetOfEdges4::find \nplus de tete de liste\n";
    MeshError(888);}
  Int4 n = tete[ Abs( ii ) % nx ];
  
  while (n >= 0) 
    if (ii == Edges[n].i && jj == Edges[n].j)
      return n;
    else n = Edges[n].next;
  return -1; // n'existe pas
}

Int4 SetOfEdges4::add(Int4 ii,Int4 jj)
{
  if (tete == 0 ) {
    cerr << "SetOfEdges4::add\n plus de tete de liste \n" << endl;
    MeshError(888);}
  
  Int4 h;
  Int4 n = tete[ h = Abs( ii ) % nx ];
  while (n >= 0) 
    if (ii == Edges[n].i && jj == Edges[n].j)
      return n;
    else n = Edges[n].next;
  if (nbax <=NbOfEdges ) {
    cerr << " SetOfEdges4::add\noverflow de la pile "  << nbax << " " << NbOfEdges << endl;
    MeshError(888);}
  
  Edges[NbOfEdges].i=ii;
  Edges[NbOfEdges].j=jj;
  Edges[NbOfEdges].next= tete[h];
  tete[h] = NbOfEdges;
  return NbOfEdges ++;
}

