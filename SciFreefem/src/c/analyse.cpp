//file Analyse.cpp
//#include <profiler.h>

#include <cmath>
#include <iostream>
#include <fstream>

//#include <ccstream.h>
#include <assert.h>
//#define NDEBUG
#include <ctype.h>
#include <stdlib.h>
#include "vect.h"
#include "Mesh2.h"
#include "BamgFreeFem.h"

#include "QuadTree.h"
#include "graph.h"
#include "rgraph.h"
#include "convec.h"
#include "ivarsol.h"
#include "clerror.h"
#include "list.h"

#include <sstream>

#define istrstream istringstream

// -------------------------
//int FH = 0;  // debug messages for Frederic (use !=0 for talkative )

//float penal = 1e5; // default value for the Dirichlet penalization
// same order as typedef enum in analyse.h
// Name of  symboles
static char * SymbolName [] = {
  "(",                  ")",             "{",              "}",                 " constante ", 
  " identifier ",     "function ",     "+",              "-",                 "*",
  "/",                "^",             "<",              "<=",                ">", 
  ">=",               "==",            "!=",             ",",                 ";", 
  ":",                "&&",            "||",             " if ",              " then ", 
  " else ",           " for ",         " to ",           " do ",              " erreur ", 
  " end ",            "=",             " mesh ",         " border ",          " function2 ", 
  " derive ",         " plot ",        " plot3d ",       " print ",           " R2 ",
  " string ",         "buildmesh ",    " savemesh ",     " readmesh ",        " movemesh ",
  " save ",           " read ",        " exit ",         " helmholtz ",       " adaptmesh ",
  " int ",            " Id ",          " on ",           " assemble ",        " varsolve ", 
  " convect ",        " dx ",          " dy ",           " array ",           " laplace ",
  " dxx ",            " dxy ",         " dyx ",          " dyy ",             " dnu ", 
  " pde ",            " solve ",       " with ",         " append",           " number ",
  " subroutine ",     "~",             " femp0 ",        " femp1 ",           " plotp0 "
};

// name of the argument of  AdapMesh : warning the order is important 
const int NbArgAdapMesh = 18;
static char * ArgAdapMesh[NbArgAdapMesh] = {
  "hmin","hmax","err","errg","nbvx","nbsmooth","nbjacoby","ratio","omega","iso",
  "abserror","cutoff","verbosity","inquire","splitpbedge","maxsubdiv","anisomax"
  ,"rescaling"};

// Fonctions auxiliaires  
// FH j'ai supprime les statics PB ordre des parameters dans gcc sur hp


/*static*/ float Lt (float x, float y) {return x < y;}
/*static*/ float Ge (float x, float y) {return x >= y;}
/*static*/ float Gt (float x, float y) {return x > y;}
/*static*/ float Eq (float x, float y) {return x == y;}
/*static*/ float Ne (float x, float y) {return x != y;}
/*static*/ float Le (float x, float y) {return x <= y;}
/*static*/ float Et (float x, float y) {return x && y;}
/*static*/ float Ou (float x, float y) {return x || y;}
/*static*/ float aMin (float x, float y) {return x < y ? x : y;}
/*static*/ float aMax (float x, float y) {return x > y ? y : x;}
/*static*/ float Abs (float x) {return x > 0 ? x : -x;}
/*static*/ float Pow (float x, float y) {return pow(x, y);}
/*static*/ float Sqrt (float x) {return sqrt(x);}
/*static*/ float Sin (float x) {return sin(x);}
/*static*/ float Exp (float x) {return exp(x);}
/*static*/ float Log (float x) {return log(x);}
/*static*/ float Cos (float x) {return cos(x);}
/*static*/ float Tan (float x) {return tan(x);}
/*static*/ float Atan (float x) {return atan(x);}
/*static*/ float Atan2 (float x,float y) {return atan2(x,y);}


// -- some graphical variable --
static int ginit=0;
static float gxctr,gyctr,gray;

//--------------------- data for link to Scilab -----------------------

bool toScilab = false;

class DataToScilab
{
public:
  double * an_eval;
  Grid * g;
  
  DataToScilab(Grid * gr) : g(gr) { an_eval = new double [g -> nv]; }
  ~DataToScilab() { if (an_eval) delete [] an_eval; }
} * dts;

//---------------------------------------------------------------------

// Fonction blip = Fonction (1) / (Fonction (1) + Fonction::monome (1, 2));
// Fonction Atan (atan, blip);

// Gestion de la table d'identificateurs
  
Iden *IdenTable::insert (const char *s, int n)
{
  const int increment = 100;
  if ((nb % increment) == 0)
    {
      Iden **newList = new Iden* [nb + increment];
      for (int i = 0; i < nb; i++) newList[i] = list[i];
      delete [] list; list = newList;
    }
  for (int i = nb++; i > n; i--) list[i] = list[i - 1];
  return list[n] = new Iden(s);
}

Iden* IdenTable::find (const char *s)
{
  if ((nb == 0) || (strcmp(list[0]->nom, s) > 0)) return insert(s, 0);
  int l = 0, r = nb;
  while ((r - l) > 1)
    {
      int m = (l + r) / 2;
      if (strcmp(list[m]->nom, s) > 0) r = m; else l = m;
    }
  return strcmp(list[l]->nom, s) ? insert (s, l + 1) : list[l];
}

int IdenTable::remove (const char *s)
{
  /*      Iden **newList = new Iden* [--nb];  
          int j=0;
          for (int i = 0; i <= nb; i++)
          if(j==nb) return 1;
          else if(strcmp(list[i]->nom, s)) newList[j++] = list[i];
          else delete list[i]->storage;
          delete [] list; list = newList;
  */      return 0;
}
void IdenTable::enregistre(const char *s, Symbol sym)
{
  Iden *i = find(s); assert(i->type == Iden::inconnu);
  i->type = Iden::reserve; i->sym = sym;
}
void IdenTable::enregistre(const char *s, float * f)
{
  Iden *i = find(s); assert(i->type == Iden::inconnu);
  i->type = Iden::pfloat; i->pf = f;
}
void IdenTable::enregistre(const char *s, int * l)
{
  Iden *i = find(s); assert(i->type == Iden::inconnu);
  i->type = Iden::pint; i->pi = l;
}
void IdenTable::enregistre(const char *s, const Fonction &f)
{
  Iden *i = find(s); if(i->type != Iden::inconnu)  {
    cout << "PB  IdenTable::enregistre"  << i->type << i->nom << endl;}
  i->type = Iden::fonction; i->f = f;
}
void IdenTable::enregistre(const char *s, const Fonction2 &f)
{
  Iden *i = find(s); assert(i->type == Iden::inconnu);
  i->type = Iden::fonction2; i->f2 = f;
}



float CTab::operator ()(float x, float y) 
{       
  int kkk=0;
  int kt = an->trloc,inside = 1;
  double a[3],xx,yy;
#ifndef NDEBUG1   
 toto:
#endif  

  //#undef NDEBUG
  //  if( an->ivertex == 0) cout <<  "cc = " << cc << " this=" << this << endl;
  assert(cc && g );
  if (g== an->gridxyng) // same grid  
    if(kt>=0)
      return  cc[g->no(g->t[kt].v[0])]*an->xl[0] +cc[g->no(g->t[kt].v[1])]*an->xl[1]
        +cc[g->no(g->t[kt].v[2])]*an->xl[2];
    else if( an->ivertex >= 0) 
      { //if (an->ivertex < 10) cout << cc[an->ivertex] << " " << an->ivertex ;
	return cc[an->ivertex];
      }
  assert(g->Th && g->nv <= size);
  kt=FindTriangle(*g->Th,x,y,a,inside);
  assert(kt>=0 && kt < g->nt);
  int i0= g->no(g->t[kt].v[0]); 
  int i1= g->no(g->t[kt].v[1]); 
  int i2= g->no(g->t[kt].v[2]); 
  /*    if(an->local) 
        {
        a[0] = ((x - g->v[i2].x)*(g->v[i1].y - g->v[i2].y) 
        - (y - g->v[i2].y)*(g->v[i1].x - g->v[i2].x))/(2*g->t[kt].area);
        a[1] = ((x - g->v[i0].x)*(g->v[i2].y - g->v[i0].y) 
        - (y - g->v[i0].y)*(g->v[i2].x - g->v[i0].x))/(2*g->t[kt].area);
        a[2] = ((x - g->v[i1].x)*(g->v[i0].y - g->v[i1].y) 
        - (y - g->v[i1].y)*(g->v[i0].x - g->v[i1].x))/(2*g->t[kt].area);
        }
  */
#ifndef NDEBUG1    
  // verification
  xx = g->v[i0].x*a[0] + g->v[i1].x*a[1] + g->v[i2].x*a[2]; 
  yy = g->v[i0].y*a[0] + g->v[i1].y*a[1] + g->v[i2].y*a[2];
  double err =  Max(Abs(xx-x),Abs(yy-y));
  if (err > 1e-3 && inside && !kkk) 
    {
      //  cout << "\t\t" << x << ", " << y << " localisation  err:" << err << " bary  = " 
      //      <<  a[0] << " " << a[1]<< " " <<  a[2]<< endl;
      kkk++;
      goto toto;
    }
#endif
  assert(cc);
  //if (an->ivertex < 10) 
  // cout << an->ivertex << " " << cc[i0]*a[0] + cc[i1]*a[1] +cc[i2]*a[2] <<  endl;
  return cc[i0]*a[0] + cc[i1]*a[1] +cc[i2]*a[2];
}

float  P0::F(int iv) const {
  // on prend la moyenne 
  double f=0,a=0;
  Triangles &Th(*g->Th);
  int nt=g->nt;
  Vertex &v = Th.vertices[iv];
  Triangle *t=v.t, *tt=t, *ttc;
  int j = EdgesVertexTriangle[v.vint][0], jc;
  int k=0;
  do {
    int kt = Th.Number(tt);
    if (kt<nt) {
      a += tt->det;
      f += tt->det*cc[Th.Number(tt)];
    }
    ttc =  tt->TriangleAdj(j);
    jc = NextEdge[tt->NuEdgeTriangleAdj(j)];
    tt = ttc;
    j = NextEdge[jc];
    assert(k++<2000);        
  } while (t!=tt);
  //  cout << " k = " << k << " " << f/ a << " " << v <<endl;
  return f/a;
}

float CTabP0::operator ()(float x, float y) 
{       
  int kkk=0;
  int kt = an->trloc,inside = 1;
  double a[3],xx,yy;
  assert(g->Th && g->nv <= size);
  //  cout << an->local << " " << g << " " << an->gridxyng << " " << kt << endl;
  if (g== an->gridxyng) // same grid  
    if(kt>=0) {return cc[kt];}
    else if( an->ivertex >= 0) 
      return F(an->ivertex);;
  kt=FindTriangle(*g->Th,x,y,a);
  if (kt>=0) {return cc[kt];}
  else       {return 0;}
}

CTab* Iden::iden2CTab(Iden * grididen, Analvar* an)
{
  if (type ==inconnu)   
    {
      newVar();
      type = ftableau;
      CTab *fff = new CTab (grididen, an);
      return fff;
    } else if (type == Iden::ftableau) return ft;
  else throw(ErrorCompile(" Function must be an array-function "));
  return 0;
}
CTabP0* Iden::iden2CTabP0(Iden * grididen, Analvar* an)
{
  if (type ==inconnu)   
    {
      newVar();
      type = ftableauP0;
      CTabP0 *fff = new CTabP0 (grididen, an);
      return fff;
    } else if (type == Iden::ftableauP0) return fP0;
  else throw(ErrorCompile(" Function must be an array-function "));
  return 0;
}
// eval   string expression
char *EChaine::eval(char *extention)
{
  EChaine * c;
  int l= extention ? strlen(extention) : 0;
  for  (c=this;c;c=c->next) // for all sub chaine 
    {
      if(c->s) l+=strlen(c->s)+1;
      if(c->e) l+=32;
    }
  char *result = new char[++l];
  result[0]='\0';
  char * e = result;
  for  (c=this;c;c=c->next)
    {
      if(c->s)  strcpy(e, c->s), e+=strlen(c->s);
      if(c->e)  e+=sprintf(e,"%g",c->e->eval());
    }
  if (extention)  strcpy(e,extention);
  return result; 
}

typedef struct ArgNode {
  Iden * x;
  ArgNode *n;
  ArgNode (Iden *xx, ArgNode *nn) : x (xx), n(nn) {}
} *ArgList;


struct OnNode {
  Iden * x;
  Expr * e;
  OnNode *n;
  OnNode (Iden *xx, OnNode *nn) : e(0) ,x (xx), n(nn) {}
  OnNode (Expr *ee, OnNode *nn) : e(ee), x (0), n(nn) {} 
  int  ref();// { if (x) return ((IB*)x->fn)->label; OnNode
  //  else  return (int) (0.5+e->eval());}
} ;

typedef struct ArgENode {
  Expr * x;
  ArgENode *n;
  ArgENode (Expr *xx, ArgENode *nn) : x (xx), n(nn) {}
} *ArgEList;

struct CFN { //declaration of functions with an indefinite number of parameters
  Iden *nom;
  ArgList larg;
  Instr *e;
  CFN (Iden *n, ArgList l, Instr *ee) : nom(n), larg (l), e(ee) {}
};

class EFN : public Expr {//definition of functions with an indefinite number of parameters
  CFN *f;
  ArgEList learg;
public:
  EFN (CFN *ff, ArgEList lea) : f(ff), learg (lea) {}
  float eval ();
  Expr *d (Iden *x) {return this;} // can't work...
};

float EFN::eval ()
{
  ArgList a; ArgEList ae;
  for (a = f->larg, ae = learg; a; a = a->n, ae = ae->n)
    *a->x->storage = ae->x->eval ();
  f->e->execute ();
  return *f->nom->storage;
}

struct IB 
{ // Definition de segment de frontiere
  float *x, *y, *ng, *t;
  Expr *start, *stop;
  Instr *l1;
  int label;
public:
  IB( float *xx, float *yy, float *nng, float *tt, Expr *sstart,
      Expr *sstop, Instr *ll1, int llabel) : x(xx), y(yy), ng(nng),
                                             t(tt), start(sstart), stop(sstop), l1(ll1), label(llabel) {}
  //  void execute();
};


int  OnNode::ref()
{ 
  if (x)
    return ((IB*)x->fn)->label;
  else  
    return (int) (e->eval());
}

class EB : public Expr
{ // discretisation of a border
  IB* ib;
  Expr* npas;
  Analvar* an;
public:
  EB(Analvar* aan, IB* iib, Expr* nnpas) : an(aan), ib(iib), npas(nnpas){}
  float eval();
  Expr *d (Iden *x) {cout<<"curves shouldn't be differentiated"<<endl; return 0;} 
} ;

float EB::eval()
{ // execute a border construction
  frontiere& front = *an->front;
  int num=0, oldnum, j=1;
  float a = ib->start->eval();
  float b = ib->stop->eval();
  int n = (int)npas->eval(); 
  if(n<0){j=-1; n=-n;}  // used to leave the domain on the right of the bdy
  if (n < 2) n = 2;
  float delta = (b - a) / (n - 1);
  *ib->ng=ib->label;  // set default value 
  
  for (int i = 0; i < n; i++)
    {
      *ib->t = a + i * delta;
      ib->l1->execute();
      oldnum = num;
      if (front.step) 
        {
          num = front.addPoint(*ib->x, *ib->y, (int)(*ib->ng));
          *ib->t = a + i * delta + delta/2; // add  modif FH 290497   
          ib->l1->execute();//          
          if (i) front.addSegment(oldnum-1, num-1, (int)(*ib->ng)); // fin add  modif FH 290497  
        }
      else if(front.initboundingbox)
        {  
          front.xmin= Min(*ib->x,front.xmin);
          front.xmax= Max(*ib->x,front.xmax);
          front.ymin= Min(*ib->y,front.ymin);
          front.ymax= Max(*ib->y,front.ymax);
        }
      else
        {  
          front.initboundingbox =1;
          front.xmin= *ib->x;
          front.xmax= *ib->x;
          front.ymin= *ib->y;
          front.ymax= *ib->y;
        }
    }
  if (front.step) 
    {
      front.hh[num-1] = front.hh[oldnum-1]; 
      front.sd[2*front.nbsd]= j*(front.nbs-1); 
      front.sd[2*front.nbsd+1]=front.nbsd+1; 
      front.nbsd++;
    }
  return 0.;
}

class Eintegral:public Expr
{ // calcul d'integral
  Expr* e;
  OnList al;
  Analvar* an;
public:
  Eintegral (OnList ial,Expr *ee,Analvar *aan): e(ee), al(ial),an(aan) {}
  float eval();
  void initquad(Grid& coarse, Grid& fine); // in ivarsolve.cpp
  Expr *d (Iden *x) {cout<<"Integrals cannot be differentiated"<<endl;return 0;} 
};


float Eintegral::eval()
{       
  float integr,integral=0;
  Analvar save(*an);
  int ref,i,iloc,inext, inextnext, 
    antrlocold =  an->trloc, localold = an->local, ilocold = an->iloc, ivertexold = an->ivertex;
  int k = an->local ? an->trloc : 0;  // set k to first triangle if not local
  // float xold = *(an->x->storage), yold = *(an->y->storage);
  // float ngold = *(an->ng->storage), xlold[3];
  // xlold[0] = an->xl[0]; xlold[1] = an->xl[1]; xlold[2] = an->xl[2];
  Grid* oldg = an->activeMesh;
  
  if(al->n==0) // it is a volume integral (only Grid is parameter)
    {
      Grid& t= *((Grid*)al->x->fn);
      an->gridxyng=&t;
      assert(t.nv && t.nt);
      an->activeMesh = &t;//(Grid*)al->x->fn;
      an->gridxyng = &t;
      if((&t != oldg)&& localold) 
        { 
          if(&t!=oldg->gridFriend)
            (*oldg).initquad(&t);  // finds all quadrature pt of one mesh in triangle of other
          integral = 0;
          for( int j=0;j<oldg->quad[k][0];j++)
            {
              int kk = oldg->quad[k][j+1] / 3;
              iloc = oldg->quad[k][j+1] -3*kk;
              bTriangle& tk = t.t[kk];  // slightly inside to handle discontinuous functions
              int jloc=::next[iloc], kloc=::next[iloc+1];
              float xx = (tk.v[iloc]->x + tk.v[jloc]->x)/2;// + 0.001*tk.v[kloc]->x)/2.001;
              float yy = (tk.v[iloc]->y + tk.v[jloc]->y)/2;// + 0.001*tk.v[kloc]->y)/2.001;
              float xl[3];
              xl[iloc] = 1/2.0;//01; 
              xl[jloc] = 1/2.0;//01; 
              xl[kloc] = 0.0;//01/2.001;
	      //              int oldlocal = an->local;
              an->setAn(1,xx,yy, tk.e[kloc]->where, xl, -1, 10*(iloc+1)+jloc, kk);
              float z =  e->eval();
	      //              an->local = oldlocal;
              integral += t.t[kk].area * z/3 ;
            }
        } else
	do {
	  integr = 0;
	  for( iloc=0;iloc<3;iloc++)
	    {
	      bTriangle& tk = t.t[k]; // slightly inside to handle discontinuous functions
	      int jloc=::next[iloc], kloc=::next[iloc+1];
	      float xx = (tk.v[iloc]->x + tk.v[jloc]->x)/2;// + 0.001*tk.v[kloc]->x)/2.001;
	      float yy = (tk.v[iloc]->y + tk.v[jloc]->y)/2;// + 0.001*tk.v[kloc]->y)/2.001;
	      float xl[3];
	      xl[iloc] = 1/2.0;//01; 
	      xl[jloc] = 1/2.0;//01; 
	      xl[kloc] = 0.0;//01/2.001;
	      //   int oldlocal = an->local;
	      an->setAn(1,xx,yy, tk.e[kloc]->where, xl,-1, 10*(iloc+1)+jloc, k);
	      integr+= e->eval();
	      // an->local = oldlocal;
	    }
	  integral += integr * t.t[k].area / 3;
	}while((++k<t.nt) && !localold);
    } else {                                                    // it's a boundary integral
    OnList a0 = al; 
    while(a0->n) { a0 = a0->n;}
    Grid& g= *((Grid*)a0->x->fn);
    an->gridxyng = &g;
    if((&g != oldg)&& localold) 
      { cout << endl <<"Can't use different meshes in an integral inside a varPDE " << endl; 
        exit(0);
      }
    do{
      bTriangle& tk = g.t[k];
      for(int iloc=0;iloc<3;iloc++)
	{
	  i = (iloc+1)%3; // i = iloc if check v.where
	  inext = (i+1)%3; 
	  int isonit = 0, refi = tk.e[iloc]->where; // tk.v[i]->where, refinext = tk.v[inext]->where;
	  for (OnList a = al; a!=a0; a = a->n)
	    { ref = a->ref(); //      ((IB*)a->x->fn)->label;
              isonit += refi== ref;
              // isonit += refinext== ref;
	    }
	  if(isonit)
	    { // slightly inside to handle discontinuous functions
	      float xx = (tk.v[i]->x + tk.v[inext]->x)/2;
	      float yy = (tk.v[i]->y + tk.v[inext]->y)/2;
	      float xl[3];
	      xl[i] = 1/2.00; xl[inext] = 1/2.00; xl[iloc] = 0.00;
	      int oldlocal = an->local;
	      an->setAn(1,xx,yy, tk.v[i]->where, xl, -1, 10*(i+1)+inext, k);
	      an->local = oldlocal;
	      integral += e->eval() * sqrt( 
					   (tk.v[i]->x - tk.v[inext]->x) * (tk.v[i]->x - tk.v[inext]->x)
					   + (tk.v[i]->y - tk.v[inext]->y) * (tk.v[i]->y - tk.v[inext]->y) );
	    }
	}
    } while((++k<g.nt) && !localold);
  }                   
  //  an->setAn(localold,xold,yold,ngold,xlold,ivertexold, ilocold, antrlocold);
  //  an->activeMesh = oldg;   
  *an=save;     
  return integral; 
}

class Eop:public Expr 
{ // calcul d'operateur
  Iden* id;
  Symbol s;
  Analvar* an;
public:
  Eop (Symbol ss,Iden *iid, Analvar* aan): id(iid), s(ss),an(aan) {}
  float eval();
  Expr *d (Iden *x) {cout<<"Operator cannot be differentiated"<<endl;return 0;} 
};

float Eop::eval()
{
  if(id->type!=Iden::ftableau){ cerr << "Symbolic derivatives Unimplemented" << endl; exit(1);}
  CTab& u = *(id->ft);
  Grid& t =*(id->ft->g);
  int kt,inside; 
  float deri;
  if(an->local && an->activeMesh == &t) kt = an->trloc; 
  else {
    double a[3];
    assert(t.Th && u.size >= t.nv);
    kt=FindTriangle(*t.Th,*an->x->storage,*an->y->storage,a,inside);
    assert(kt>=0 && kt < t.nt);
  }
  int i0= t.no(t.t[kt].v[0]);   
  int i1= t.no(t.t[kt].v[1]);   
  int i2= t.no(t.t[kt].v[2]);   
  if(s==dx) deri = ((u[i1] - u[i0]) * (t.v[i2].y - t.v[i0].y)
                    - (u[i2] - u[i0]) * (t.v[i1].y - t.v[i0].y))/(2*t.t[kt].area);
  else deri = ((t.v[i1].x - t.v[i0].x) * (u[i2] - u[i0])
               - (t.v[i2].x - t.v[i0].x) * (u[i1] - u[i0]))/(2*t.t[kt].area);
  return deri;
}

class EHeavyside:public Expr
{ // calcul d'une somme de penalisation
  OnList id;
  Analvar *an;
public:
  EHeavyside (OnList iid, Analvar* aan): id(iid), an(aan) {}
  float eval();
  Expr *d (Iden *x) {cout<<"Id cannot be differentiated"<<endl;return 0;} 
};

float EHeavyside::eval()
{
  // une seule frontiere: a faire
  Grid& g = *(an->activeMesh);
  int ref = id->ref();
        
  if(an->local && an->ivertex>=0)       
    return g.v[an->ivertex].where==ref;

  float ruff, x=*(an->x->storage),      y=*(an->y->storage);    
  int kt,inside = 1;
  double a[3],xx,yy;
  assert(g.Th);
  kt=FindTriangle(*g.Th,x ,y ,a,inside);
  assert(kt>=0 && kt < g.nt);
  int i0= g.no(g.t[kt].v[0]);   
  int i1= g.no(g.t[kt].v[1]);   
  int i2= g.no(g.t[kt].v[2]);   
  a[0] = ((x - g.v[i2].x)*(g.v[i1].y - g.v[i2].y) 
          - (y - g.v[i2].y)*(g.v[i1].x - g.v[i2].x))/(2*g.t[kt].area);
  a[1] = ((x - g.v[i0].x)*(g.v[i2].y - g.v[i0].y) 
          - (y - g.v[i0].y)*(g.v[i2].x - g.v[i0].x))/(2*g.t[kt].area);
  a[2] = ((x - g.v[i1].x)*(g.v[i0].y - g.v[i1].y) 
          - (y - g.v[i1].y)*(g.v[i0].x - g.v[i1].x))/(2*g.t[kt].area);
#ifndef NDEBUG    
  // verification
  xx = g.v[i0].x*a[0] + g.v[i1].x*a[1] + g.v[i2].x*a[2]; 
  yy = g.v[i0].y*a[0] + g.v[i1].y*a[1] + g.v[i2].y*a[2];
  double err =  Max(Abs(xx-x),Abs(yy-y));
  if (inside && (err > 1e-5) )
    cout << " \t\tErreur localisation " << x << " , " << y << " = " << err << endl;
#endif
  ruff = (g.v[i0].where==ref)*a[0] + (g.v[i1].where==ref)*a[1] +(g.v[i2].where==ref)*a[2];
  ruff = int(ruff+0.5);
  return ruff;
}

class Epenalty:public Expr
{ // calcul d'une somme de penalisation
  Expr *e0, *e1, *e2; // on(<border>,<exp2>) (<exp0>=<exp1>)
  OnList al;
  Analvar *an;
public:
  Epenalty (OnList iid,Expr *ee0, Expr *ee2, Expr *ee1, Analvar* aan): e0(ee0), e1(ee1), 
								       e2(ee2), al(iid), an(aan) {}
  float eval();
  Expr *d (Iden *x) {cout<<"Id cannot be differentiated"<<endl;return 0;} 
};

float Epenalty::eval() 
{
  Grid& g = *(an->activeMesh);
  int ref = al->ref();
  float ruff, eee, etest; 
  if(an->local && an->ivertex >= 0)
    {   for(OnList a=al;a;a=a->n)
	{
	  ref = a->ref();
	  if(g.v[an->ivertex].where==ref) 
	    {     etest = e2->eval();
	      eee = e0->eval();
	      eee -= e1->eval();
	      return eee * etest *penal;
	    }
	}
      return 0;
    }
  float x=*(an->x->storage),y=*(an->y->storage);        
  int kt,inside = 1;
  double a[3],xx,yy;
  assert(g.Th);
  kt=FindTriangle(*g.Th,x ,y ,a,inside);
  assert(kt>=0 && kt < g.nt);
  int i0= g.no(g.t[kt].v[0]);   
  int i1= g.no(g.t[kt].v[1]);   
  int i2= g.no(g.t[kt].v[2]);   
  a[0] = ((x - g.v[i2].x)*(g.v[i1].y - g.v[i2].y) 
          - (y - g.v[i2].y)*(g.v[i1].x - g.v[i2].x))/(2*g.t[kt].area);
  a[1] = ((x - g.v[i0].x)*(g.v[i2].y - g.v[i0].y) 
          - (y - g.v[i0].y)*(g.v[i2].x - g.v[i0].x))/(2*g.t[kt].area);
  a[2] = ((x - g.v[i1].x)*(g.v[i0].y - g.v[i1].y) 
          - (y - g.v[i1].y)*(g.v[i0].x - g.v[i1].x))/(2*g.t[kt].area);
#ifndef NDEBUG    
  // verification
  xx = g.v[i0].x*a[0] + g.v[i1].x*a[1] + g.v[i2].x*a[2]; 
  yy = g.v[i0].y*a[0] + g.v[i1].y*a[1] + g.v[i2].y*a[2];
  double err =  Max(Abs(xx-x),Abs(yy-y));
  if (inside && (err > 1e-5) )
    cout << " \t\tErreur localisation " << x << " , " << y << " = " << err << endl;
#endif
  ruff = 0;
  for(OnList ai=al;ai;ai=ai->n)
    {   int ref=ai->ref();//((IB*)ai->x->fn)->label ;
      ruff += (g.v[i0].where==ref)*a[0] + (g.v[i1].where==ref)*a[1] +(g.v[i2].where==ref)*a[2];
    }
  ruff = int(ruff+0.5);
  return ruff;
}

/*************************/
class Econvec2: public Expr
{ // convect a function; appears within an expression
  Analvar* an;
  Iden *idmesh, *idu, *idv, *idf;
  Expr *e;
public:
  Econvec2(Iden* iid, Iden* iidf, Iden* iidu, Iden* iidv, Expr* ee, Analvar* ann)
    : idmesh(iid),idf(iidf), idu(iidu), idv(iidv), e(ee),  an(ann) {}
  float eval ();
  Expr *d (Iden *x) {cout<<"convect cannot be differentiated"<<endl;return 0;} 
};

float Econvec2::eval()// convect from within an expression
{
  Grid& g = *(idmesh->fg);
  CTab& u = *(idu->ft);
  CTab& v = *(idv->ft);
  CTab& f = *(idf->ft);
  if((u.g != v.g)||(u.g != f.g) || (u.g !=idmesh->fg)) 
    throw(ErrorExec( " Cannot convect on different meshes from within an expression "));
  double foX;
  double dt = e->eval();  
  double a[3];
  assert(g.Th);
  int inside;

  int k = FindTriangle(*g.Th,*an->x->storage,*an->y->storage,a,inside);
  if(fabs(a[0]*a[1]*a[2]) < 1e-10) // if on border of k, shift it inside
    for(int i=0;i<3;i++) a[i] = (a[i] + 1e-10)/(1+3e-10);
  int err = xtoX (g,u, v, a, &dt, &k);
  return f[g.no(g.t[k].v[0])] * a[0] 
    + f[g.no(g.t[k].v[1])] * a[1] 
    + f[g.no(g.t[k].v[2])] * a[2];
}


class Ebuildmesh: public MeshExpr
{ // mesh construction
  Analvar* an;
  Expr *e;
        
public:
  Ebuildmesh (Expr *ee, Analvar* ann,EChaine * ffile) : e(ee), an(ann),MeshExpr(ffile){}
  Grid * eval ()
  {      
    an->front->init();
    e->eval(); // 
    an->front->step=1; 
    frontiere & f(*an->front);
    assert( f.initboundingbox);
    f.epsilon = Max(f.xmax-f.xmin,f.ymax-f.ymin)*1.0e-5;                
    e->eval(); 
    Grid * g = new Grid();              
    g->buildit(an->front, *an->wait->storage);
    return g;
  }
};

class Edaptmesh: public MeshExpr
{ // adaptation de maillage
  Analvar* an;
  Expr *arg[NbArgAdapMesh];
  Expr *sol[32];
  int nbsol;
  Grid* tloc;
  Iden * bthid;// --
public:
  Edaptmesh(Analyseur & a,Iden * tt,Analvar* ann,EChaine *ffile)
    : bthid(tt),tloc(tt->fg),an(ann),MeshExpr(ffile)
  {nbsol=a.GetExprs(sol,32);a.match(rpar); // get the solution 
    a.FindArgs(NameArg(NbArgAdapMesh,(const char **) ArgAdapMesh,(Expr **)arg));// get the arg
  }
  Grid* eval ();
  //  Expr *d (Iden *x) {cout<<"adaptmesh cannot be differentiated"<<endl;return 0;} 
};

class Emovemesh:public MeshExpr
{// move the x and y of mesh vertices by expr ex and ey
  Expr *ex,*ey;
  Iden* idmoved;
  Analvar* an;
public:
  Emovemesh (Iden* iid,Expr *eex, Expr *eey, Analvar* ann,EChaine *ffile)
    : idmoved(iid), ex(eex), ey(eey), an(ann),MeshExpr(ffile) {}
  Grid* eval ();
};
class Etruncmesh:public MeshExpr
{// move the x and y of mesh vertes by expr ex and ey
  Expr *e,*b;
  Iden* idgrid;
  Analvar* an;
public:
  Etruncmesh (Iden* iid,Expr *ee, Expr *bb, Analvar* ann,EChaine *ffile)
    : idgrid(iid), e(ee), b(bb), an(ann),MeshExpr(ffile) {}
  Grid* eval ();
};


///////////////////////////////////
// Diverses classes d'instruction//
///////////////////////////////////

class IVide : public Instr
{
public:
  IVide () {}
  void execute () {}
};

Instr *Analyseur::rien = new IVide();
        
class Ibecomes : public Instr
{ // Affection
  Iden *v;
  Expr *e;
  Analvar *an;
public:
  Ibecomes (Iden *vv, Expr *ee,Analvar *ann) : v(vv), e(ee), an(ann) {}
  void execute (){
    if(v->type == Iden::pint)         *v->pi       =(int)e->eval(); 
    else if (v->type == Iden::pfloat)  *v->pf       =e->eval()  ;
    else                               *v->storage  = e->eval();}

  Instr *d (Iden *x) {assert (v != x); return new Ibecomes (v, e->d (x),an);}
};

class Iarray : public Instr
{ // affectation for an array function
  Iden *v;
  Expr *e;
  Analvar *an;
public:
  Iarray (Iden *vv, Expr *ee,Analvar *ann) : v(vv), e(ee), an(ann) {}
  void execute ();
  Instr *d (Iden *x) {assert (v != x); return new Iarray (v, e->d (x),an);}
};

class IarrayP0 : public Instr
{ // affectation for an array function
  Iden *Th;
  Iden *v;
  Expr *e;
  Analvar *an;
public:
  IarrayP0 (Iden *vv, Expr *ee,Analvar *ann,Iden *Thh=0) : Th(Thh),v(vv), e(ee), an(ann) {}
  void execute ();
  Instr *d (Iden *x) {assert (v != x); return new IarrayP0 (v, e->d (x),an);}
};


class IIf : public Instr
{ // Instruction if then else
  Expr *e; Instr *I1, *I2;
public:
  IIf (Expr *ee, Instr *i1, Instr *i2) : e(ee), I1(i1), I2(i2) {};
  void execute ()
  {if (e->eval ()) I1->execute (); else if (I2) I2->execute ();}
  Instr *d (Iden * x) {return new IIf (e, I1->d (x), I2->d (x));}
};

class Ifor : public Instr
{ // Instruction for .. to .. do
  Expr *e1, * e2; Instr *i1;
  Iden* id;
public:
  Ifor ( Iden* idd, Expr *ee1, Expr *ee2 , Instr *ii1) : id(idd), e1(ee1), e2(ee2), i1(ii1) {};
  void execute ()
  {     for(*id->storage =e1->eval(); *id->storage <= e2->eval();*id->storage += 1) 
      {  cout << "##########  iteration " << *id->storage << " ------------"<<endl;

	i1->execute();
      }

  }
  Instr *d (Iden * x) {return new Ifor (id, e1, e2, i1->d (x));}
};

class Isavemesh: public Instr
{ // Ecrit la Grid dans un fichier
  EChaine * fname;
  Iden* id;
  Analvar* an;
public:
  Isavemesh(EChaine* ffname,Iden* iid, Analvar* ann)
    :id(iid), an(ann),fname(ffname){ }
  void execute ();
};

class Iscilabmesh : public MeshExpr
{ // read grid "id" from scilab data
  EChaine * ec;
  Analvar * an;
  double * vt;  // vertex tab
  int nbvt;     // nb vertex
  int * tr;     // triangle tab
  int nbtr;     // nb triangle
public:
  Iscilabmesh(EChaine* ec_fp, Analvar * an_fp,
	      double * vt_fp,int nbvt_fp,int * tr_fp, int nbtr_fp)
    : ec(ec_fp),an(an_fp),vt(vt_fp),nbvt(nbvt_fp),tr(tr_fp),nbtr(nbtr_fp) { }
  Grid * eval ();
};

class Ireadmesh: public MeshExpr //Instr
{ // read grid "id" from file
  EChaine * fname;
  Analvar* an;
  Expr *renu;

public:
  Ireadmesh(EChaine* ffname, Analvar* ann,EChaine * ffile,Expr *rr)
    : an(ann),MeshExpr(ffile),fname(ffname),renu(rr){ }
  Grid * eval ();
};

class Isave:public Instr
{// Ecrit une expression dans un fichier
  Analvar* an;
  Expr *e;
  EChaine * fname;
  Iden* id;
public:
  Isave (Iden* iid, EChaine* ffname,Expr *ee, Analvar* ann)
    :id(iid), e(ee), an(ann) ,fname(ffname){}
  void execute ();
};
                
class Iread:public Instr
{// reads a function from a file
  CTab *f;
  EChaine * fname;
  Iden* id;
public:
  Iread (Iden* iid, EChaine* ffname, CTab *ff)
    :id(iid), f(ff),fname(ffname){}
  void execute ();
};

/*
  class Ihelmholtz: public Instr
  { // Solves a Helmholtz eq
  CTab *f2;
  Analvar* an;
  Iden *id1, *id2;   //mesh name and matrix
  Expr *l1, *l2,*l3,*l4,*l5,*l6,*l7;
  int factorize;
  public:
  Ihelmholtz(Iden* iid1, Iden* iid2, CTab *vv, Expr* ll1, Expr* ll2, Expr* ll3, Expr* ll4, 
  Expr* ll5, Expr* ll6, Expr* ll7, Analvar* ann, int dd)
  : id1(iid1),id2(iid2), f2(vv), 
  l1(ll1), l2(ll2), l3(ll3), l4(ll4), l5(ll5), l6(ll6), l7(ll7),
  an(ann), factorize(dd) {}
  void execute ();
  };

*/
class Iprint:public Instr
{// imprime une expression
  float *x, *y, *ng;
  EChaine* e;
  EChaine * fname;
  int ondisk;
public:
  Iprint (int oondisk, EChaine* ffname, EChaine* ee, float *xx, float *yy, float *nng)
    : ondisk(oondisk), e(ee), x(xx), y(yy), ng(nng) ,fname(ffname)  { }
  void execute();
};
// Correction F.H 11/19/1998
void Iprint::execute ()
{ 
  if(e){
    char * xx = e->eval();
    cout << xx ;
    if(ondisk && fname) 
      { 
	char * ff= fname ? fname->eval():0 ; 
        ofstream mfile(ff,ios::app);
        if(mfile.bad()) mfile.open(ff);
        if (mfile) 
          mfile<<xx<<endl;
        else  
	  throw(ErrorExec("Error in Iprint in opening output file"));
        delete ff;
      }  
    delete xx;
  }
  cout << endl;
  
  
}

class Iplot:public Instr
{// Trace une expression liste  
  Analvar* an;
  OnNode * list;
  EChaine * fileps;
public:
  Iplot (/*Iden* iid,*/OnNode * ll, Analvar* ann,EChaine * ffile)
    : /*id(iid),*/ list(ll), an(ann),fileps(ffile) {}
  void execute ();
};

class IplotP0:public Instr
{// Trace une expression liste  
  Analvar* an;
  OnNode * list;
  EChaine * fileps;
public:
  IplotP0 (/*Iden* iid,*/OnNode * ll, Analvar* ann,EChaine * ffile)
    : /*id(iid),*/ list(ll), an(ann),fileps(ffile) {}
  void execute ();
};

class Iplot3d:public Instr
{// plot an expression in 3d
  Analvar* an;  //  information from freefem program. ex x,y,
  Expr* e;      // function to plot
  Iden* id;     // which mesh
  EChaine * fileps;
public:
  Iplot3d (Iden* iid,Expr *ee, Analvar* ann, EChaine  *ffile)
    : id(iid), e(ee), an(ann),fileps(ffile) {}
  void execute ();
};

class Iconvec: public Instr
{ // convect one or 2 functions (global operator)
  Analvar* an;
  CTab *f1,*f2;
  Iden* id;
  Expr *e1,*e2,*e3,*e4,*e5;
public:
  Iconvec(Iden* iid, CTab *ff1, CTab *ff2, Expr* ee1, Expr* ee2, Expr* ee3, 
          Expr* ee4, Expr* ee5, Analvar* ann)
    : id(iid),f1(ff1),f2(ff2), e1(ee1), e2(ee2), e3(ee3), e4(ee4), e5(ee5), an(ann) {}
  void execute ();
};

class IS : public Instr
{ // Suite d'instructions
  Instr *I1, *I2;
public:
  IS (Instr *i1, Instr *i2) : I1(i1), I2(i2) {};
  void execute () {I1->execute(); I2->execute();}
  Instr * d (Iden * x) {return new IS (I1->d(x), I2->d(x));}
};


class MeshCode : public Instr
{ // called whenever "mesh" keyword is used
  MeshExpr *e;
  Iden* id;
  Analvar* an;
public:
  MeshCode (Iden* iid, MeshExpr  *ee,Analvar* ann) 
    : id(iid), e(ee), an(ann) {}
  void execute()
  {
      
    if (e) 
      { // create the new grid
        if (e->fileps ) {char *ff= e->fileps->eval(".ps");  openPS(ff);delete ff;}
        Grid * g = e->eval() ; 
        if(id->fg) 
          { // remove after because  the back ground can be the old 
            cout << "\t\t" << " DELREF " << id->nom << " fg "  
                 << id->fg << " " <<  id->fg->Th->identity 
                 << " NbRef=" << id->fg->NbRef <<  endl;
            id->fg->DelRef(); 
          }
 
        id->fg = g;
        if(e->fileps) closePS();
      }
    cout << "\t\t" << "set the active mesh " << id->nom << " g = " << id->fg  << endl;
    assert( id->fg); 
    an->activeMesh =id->fg ;// set the active mesh 
    assert(an->activeMesh);
  }
};

class Iassemble: public Instr
{ // Resout une EDP
  CTab *f1,*f2;
  Iden* id;
  Expr* e;
  Analvar* an;
public:
  Iassemble(Iden* iid, CTab *ff1, CTab *ff2, Expr* ee, Analvar* aan)
    : id(iid),f1(ff1), f2(ff2), e(ee), an(aan) {}
  void execute ();
};

void Iassemble::execute() 
{
  Analvar save(*an);
  an->activeMesh = (Grid*)id->fn;
  Grid& t= *((Grid*)id->fn);
  CTab& rf1 = *f1;
  CTab& rf2 = *f2;
  rf1.resize(&t);
  rf2.resize(&t);
  float xl[3];
  int i;
  for(i=0;i<3;i++) xl[i] = 0;
  for( i=0;i<t.nv;i++) { rf2[i] = 0;rf1[i] = 0;}
  an->gridxyng = &t;    
  for(int k=0; k<t.nt;k++)
    for(int iloc = 0; iloc<3;iloc++)
      {
        int i = t.no(t.t[k].v[iloc]);
        xl[iloc] = 1;
        int oldlocal = an->local;
        an->setAn(1,t.v[i].x, t.v[i].y, t.v[i].where, xl, i, iloc,k);
        rf1[i] = 1; // basis function
        rf2[i] += e->eval();
        an->local = oldlocal;
        rf1[i] = 0;
        xl[iloc] = 0;
      }
  *an = save;  
} 

class Ipde: public Instr
{       // defines a pde within solve
  EDP* edp;             // parent block of PDEs
  int jedp;             // which PDE within edp
  Expr** exp;
  int* addmul;
  Instr* l;
  Analvar* an;
public:
  Ipde(EDP* eedp, int nn, Expr** eexp, int* aaddmul, Instr* ll, Analvar* aan):
    edp(eedp), jedp(nn), exp(eexp), addmul(aaddmul), l(ll), an(aan){}
  void execute ();
};

float addmulop(Expr* e, int k)
{
  switch (k) {
  case -1:return 0;                             break; // operator to used
  case 0:       return 1;                               break; // + op
  case 1: return -1;                            break; // - op
  case 2: return 1./e->eval();  break; // +op/exp
  case 3: return -1./e->eval();         break; // -op/exp
  case 4: return  e->eval();            break; // +op*exp
  case 5: return -e->eval();            break; // -op*exp
  default: throw(ErrorExec(" Internal bug please fill bug report ")); //return 0;
  }
}

void Ipde::execute()
{
  Analvar save(*an);
  an->activeMesh->check();
  Grid* g = an->activeMesh;
        
  l->execute();
  an->gridxyng = g;
  float addm6, addm9, xx, yy;
  float xl[3], xxl[3] = {1./3.,1./3.,1./3.};
          
  for ( int i = 0; i < g->nt; i++)
    { 
      xx = 0;   yy = 0; bTriangle& tk = g->t[i];
      for(int iloc = 0; iloc<3;iloc++)    
        {   
          int jloc = ::next[iloc], kloc = ::next[iloc+1];
          xx += tk.v[iloc]->x/3; 
          yy += tk.v[iloc]->y/3;
          xl[iloc]=0.5;
          xl[jloc]= 0.5;
          xl[kloc]=0;
          int oldlocal = an->local;
          an->setAn(1, (tk.e[kloc]->in->x + tk.e[kloc]->out->x)/2, 
		    (tk.e[kloc]->in->y + tk.e[kloc]->out->y)/2,
		    tk.e[kloc]->where, xl, -1, 10*(iloc+1)+jloc, i);
          float aux = tk.e[kloc]->where ? 1 : 0.5;
          edp->rhs[g->no(tk.e[kloc])*edp->n+jedp]  += aux * addmulop(exp[0], addmul[0]);
        }
      int oldlocal = an->local;
      an->setAn(0,xx,yy, tk.where, xxl, -1,100, i);
      for(int kvar=0; kvar<edp->n;kvar++) // loop on all unknowns
        {
          addm6 = addmulop(exp[10*kvar+6], addmul[10*kvar+6]);
          addm9 = addmulop(exp[10*kvar+9], addmul[10*kvar+9]);
          int j = edp->n * (g->no(&tk)*edp->n+jedp) + kvar;
          edp->dis[j]  = addmulop(exp[10*kvar+1], addmul[10*kvar+1]);
          edp->pdx[j]  = addmulop(exp[10*kvar+2], addmul[10*kvar+2]);
          edp->pdy[j]  = addmulop(exp[10*kvar+3], addmul[10*kvar+3]);
          edp->dif[j]  = addmulop(exp[10*kvar+5], addmul[10*kvar+5]) + (addm6 + addm9)/2;
          edp->asym[j]  = (addm6 - addm9)/2;
          edp->pdxy[j]  = addmulop(exp[10*kvar+7], addmul[10*kvar+7]);
          edp->pdyx[j]  = addmulop(exp[10*kvar+8], addmul[10*kvar+8]);
        }
    }
  *an=save;
}

class Isur: public Instr
{       // defines a boundary condition within pde
  EDP* edp;
  int jedp;     // which PDE of edp
  OnList larg;
  Expr** exp;
  int* addmul;
  Analvar* an;
public:
  Isur(EDP* eedp, int nn, OnList llarg,  Expr** eexp, int* aaddmul, Analvar* aan): 
    edp(eedp), jedp(nn), exp(eexp), larg(llarg), addmul(aaddmul),  an(aan) {}
  void execute ();
};

void Isur::execute()
{
  Analvar save(*an);
  Grid* g = an->activeMesh;
  an->gridxyng = g;
  float xl[3];
  const double sq3 = 0.5/sqrt(3.0);
  
  for ( int k = 0; k < g->nt; k++)
    for(int iloc=0; iloc<3;iloc++)
      { 
        bTriangle& tk = g->t[k];
        int jloc = ::next[iloc], kloc = ::next[iloc+1];
        int i = g->no(tk.e[kloc]);
        int isonit=0;
        for(OnList a=larg;a;a=a->n)
          if(tk.e[kloc]->where == a->ref())                     //((IB*)a->x->fn)->label)
	    isonit = 1;
	
        if(isonit)
          if(addmul[5*jedp+4]!=-1)                              // there is  neumann cond
            {
              xl[iloc] = 0.5;   xl[jloc] = 0.5;         xl[kloc]=0;
              float xx = (tk.v[iloc]->x+tk.v[jloc]->x)/2;
              float yy = (tk.v[iloc]->y+tk.v[jloc]->y)/2;
              int oldlocal = an->local;
              an->setAn(1,xx , yy, tk.e[kloc]->where,xl, -1,10*(iloc+1)+jloc, k);
              float r = addmulop(exp[0], addmul[0]);    // right hand side  Neumann cond
              float rr = addmulop(exp[5*jedp+4], addmul[5*jedp+4]);
              for(int kvar=0;kvar<edp->n;kvar++)
                edp->rob[edp->n * (i *edp->n +jedp) + kvar] 
                  = addmulop(exp[5*kvar+1], addmul[5*kvar+1])/rr ;
              // function evaluated at two quadra points per segment. 
              int isin,isout;
              if(tk.v[iloc] == tk.e[kloc]->in)
                {               isin = iloc; isout = jloc;}
              else{     isin = jloc; isout = iloc;}
              for(int j=0;j<2;j++)
                {
                  xl[isin] = 0.5 - (j+j-1) * sq3 ; //first quadrature point (j=0) is near e.in
                  xl[isout] = 0.5 + (j+j-1) * sq3 ;  
                  float xx =  xl[isin] * tk.v[isin]->x + xl[isout] * tk.v[isout]->x;
                  float yy =  xl[isin] * tk.v[isin]->y + xl[isout] * tk.v[isout]->y;
                  int oldlocal = an->local;
                  an->setAn(1,xx , yy, tk.e[kloc]->where,xl,-1, 10*(iloc+1)+jloc, k);
                  float r = addmulop(exp[0], addmul[0]);        // right hand side  Neumann cond
                  float rr = addmulop(exp[5*jedp+4], addmul[5*jedp+4]);// coef in front of dnu()
                  if(j==0) edp->neuin[i*edp->n +jedp] = r/rr;
                  else edp->neuout[i*edp->n +jedp] = r/rr;
                }

            }
      }
  xl[0]=0.; xl[1] = 0.; xl[2]=0;
  for(int i=0; i<g->nv;i++)
    { 
      bVertex& v = g->v[i];
      int oldlocal = an->local;
      an->setAn(0,v.x, v.y, v.where, xl, i);
      int isonit=0;
      for(OnList a=larg;a;a=a->n)
        if(v.where == a->ref())                 //((IB*)a->x->fn)->label)
          isonit = 1; 
      if(isonit)
        if(addmul[5*jedp+4]==-1)                                // there is no neumann cond
          {     
            float r = addmulop(exp[0], addmul[0]);      // right hand side f Dirichlet
            if(r==0) edp->sol[i*edp->n +jedp] = 1.0e-20; //Dirichlet homogeneous
            else edp->sol[i*edp->n +jedp] = r/addmulop(exp[5*jedp+1], addmul[5*jedp+1]); 
          }
      /*          else{ // remove this if you go back to the two quadra formula per edge
                  float r = addmulop(exp[0], addmul[0]);        // right hand side  Neumann cond
                  float rr = addmulop(exp[5*jedp+4], addmul[5*jedp+4]);// coef in front of dnu()
                  edp->neu[i*edp->n +jedp] = r/rr;
                  }
      */        }
  *an=save;
}


//////////////////////////////
//  L'analyseur syntaxic   //
/////////////////////////////


void Analyseur::erreur (const char *s, const char *t)
{
  cout <<endl << s;
  if (t) cout << t;
  cout << endl;
  exit(1);
}

// Lexical analyser
inline int isal_num(int c){ return isalnum(c) || c=='_';} 
void Analyseur::nextSym()
{ 
  char c,caux;
  const  char CR(char(13)),LF(char(10));// to avoid PB of  file  
  // comming from other system mac, pc, or unix 
  int incomment =0;
  curIsAlphaNum =0; // not a alphanum global  for use in NameArg
  
  do{
    incomment = 0;
    c = source -> get();

    while (isspace(c)){ c=source->get();} 
    if(c=='/')
      { caux=source->get(); 
	//   caux = source->peek(); 
	if(caux =='/') incomment = 1;
	else if (caux == '*' ) incomment = 2;
	else source->putback(caux);
      }
    if(incomment==1) 
      { do c=source->get(); 
        while( c!= '\n' && c!= CR  && c != LF && c != EOF );
      }
    else if(incomment==2) 
      { do {    
	  c=source->get(); 
	  caux = source->peek();
	} 
        while(c != EOF && !(c=='*' && caux=='/') );
	if(c != EOF)  
	  {       c = source->get(); c = source->get(); }
	else throw(ErrorCompile( " Unterminated comment"));
      }
  } while (incomment);

  bool EndCompile = false;
  if (c == (char)EOF) { curSym = ::end; EndCompile = true; }
  else if (isdigit(c)) {source->putback(c);*source >> curVal; curSym = cste;}
  else if (isalpha(c))
    {
      int i; char buf[256];
      buf[0]=c;
      for (i = 1; i < 256 && isal_num(source->peek()); i++) *source >> buf[i];
      if (i == 256) throw(ErrorCompile("Identifier too long"));
      buf[i] = 0;
      curIsAlphaNum = 1;
      curIden = table.find(buf);
      curSym = curIden->type == Iden::reserve ? curIden->sym : iden;
    }
  else switch (c)
	 {
	   int i;
	   char c;
	 case '"': curSym = chaine; 
	   for (     i = 0,c=source->peek(); 
		     i < 256 &&  (isprint(c) && c !='\n' && c !='"');
		     c=source->peek(), i++
		     ) 
	     curChaine[i]=source->get();
	   if (i == 256) throw(ErrorCompile("String too long"));
	   curChaine[i] = 0;
	   if(source->get() != '"') //cerr <<endl<<"String='" << curChaine<< "'" ,
	     throw(ErrorCompile("End of String could not be found"));
	   break;
	 case '{': curSym = lbrace; break;
	 case '}': curSym = rbrace; break;
	 case '(': curSym = lpar; break;
	 case ')': curSym = rpar; break;
	 case ',': curSym = comma; break;
	 case '~': 
	   curSym = tilde; break;
	 case ';': curSym = semicolon; break;
	 case ':': curSym = colon; break;
	 case '+': curSym = _plus; break;
	 case '-': curSym = _minus; break;
	 case '*': curSym = star; break;
	 case '/': curSym = slash; break;
	 case '^': curSym = power; break;
	 case '<':
	   if (source->peek() == '=') {source->get(); curSym = le;} else curSym = lt; break;
	 case '>':
	   if (source->peek() == '=') {source->get(); curSym = ge;} else curSym = gt; break;
	 case '=':
	   if (source->peek() == '=') {source->get(); curSym = eq;} else curSym = becomes; break;
	 default: throw(ErrorCompile(" Unexpected character"));
	 }
  /* Prints the program on the screen */
  coutmode(1);
  switch (curSym)
    {
    case cste: cout << curVal; break;
    case iden: cout << curIden->nom; break;
    case chaine: cout << '"' << curChaine << '"' ;break;
    default:   if (!EndCompile) cout << SymbolName[curSym];
    }
  coutmode(0);
  
  if(curSym==semicolon) 
    {  cout << endl; coutmode(0); }
}

int Analyseur::IsSym(Symbol s)
{
  int r = (s == curSym);
  if (r) nextSym();
  return r;
}

void Analyseur::match(Symbol s)
{
  if (s == curSym) nextSym();
  else 
    {
      GestChar Text("Unexpected symbol: ");
      Text = Text + SymbolName[curSym];
      throw ErrorCompile(Text.Data());
    }
}


Expr * Analyseur::facteur()
{
  Expr *res;
  int haspar=0;
  switch(curSym)
    {
    case lpar:
      nextSym(); res = expression(); match(rpar); break;
    case cste: 
      res = new EC(curVal); nextSym(); break;
    case iden:
      { 
        Iden *id = curIden;
        switch(id->type)
          {
          case Iden::pint:
            res = new EVpint(id); nextSym(); break;
          case Iden::pfloat:
            res = new EVpfloat(id); nextSym(); break;          
          case Iden::variable:
          case Iden::matrix:
            res = new EV(id); nextSym(); break;
          case Iden::fonction:
            nextSym(); match(lpar);
            res = new EF(id, expression()); match(rpar); break;
                  
          case Iden::fonction2:
            { nextSym(); 
	      haspar=0;
	      if(curSym==lpar) { haspar=1; match(lpar); }
	      Expr *arg1 = expression(); match(comma); 
	      res = new EF2(id, arg1, expression()); 
	      if(haspar) match(rpar);
            } break;

          case Iden::ftableau:
          case Iden::ftableauP0:
            { nextSym(); 
	      if(curSym==lpar) { haspar=1; match(lpar); }
	      if( haspar && curSym != rpar) 
		{
		  Expr *arg1 = expression(); match(comma); 
		  res = new EF2(id, arg1, expression()); 
		} 
	      else 
		{ // to allow f() instead of f(x,y)
		  Expr* arg1 = new EV(an.x);
		  Expr* arg2 = new EV(an.y);
		  res = new EF2(id, arg1,arg2);
		}
	      if(haspar) match(rpar);
            } break;
                  
          case Iden::fonctionN:
            { nextSym(); 
	      if(curSym==lpar) { haspar=1; match(lpar); }
	      CFN *f = (CFN *)id->fn;
	      ArgList l = f->larg;
	      ArgEList le = 0;
	      while (l != 0)
		{
		  le = new ArgENode (expression(), le);
		  l = l->n;
		  if (l) match(comma);
		}
	      if(haspar) match (rpar);
	      res = new EFN(f, le);} break;
                
          case Iden::courbe:
            {
              nextSym(); match(lpar);
              an.front = front;
              res = new EB(&an,(IB*) id->fn, expression()); 
              match(rpar); break;
            }
                          
          default:
            throw(ErrorCompile(" Unknown variable"));
          }
      }
      break;  // end case ident
            
    case dx:
    case dy:
      {
        Symbol what = curSym;
        nextSym(); match(lpar);
        res = new Eop(what, curIden, &an); 
        match(iden);match(rpar); break;
      }
                          
    case intt:          // syntax: intt()(exp) or int(a,b)(exp)
    case Heavyside: // syntax : Id(a,b)
    case sur:   // syntax:  on(a,b)(u1)(exp = exp)
      {         
	Symbol theSym = curSym;
	nextSym(); match(lpar);
	Iden* id;
	if(curIden->type==Iden::maillage)
	  {    
	    id = curIden; match(iden); 
	    if(theSym!=intt|| curSym==comma)
	      {
		match(comma);
		if(curIden->type == Iden::courbe) 
		  throw(ErrorCompile(" Border name expected "));
	      }
	  }
	else id = curMesh;
	OnList larg = 0;
	larg = new OnNode (id, larg);
	if (curSym == rpar && theSym != intt) 
	  throw(ErrorCompile("Ref or Border execpted"));

	if (curSym != rpar )
	  do 
	    {     
	      if (curSym==iden && curIden->type == Iden::courbe)          
		larg = new OnNode (curIden, larg),match(iden);
	      else 
		larg = new OnNode (expression(), larg);               
	    }     while (IsSym(comma)); 
	match(rpar);
	if(theSym==intt) 
	  {       
	    match(lpar); 
	    res = new Eintegral(larg, expression(),&an); match(rpar);
	  }
	else if(theSym==Heavyside) 
	  res = new EHeavyside(larg, &an);
	else if(theSym==sur) 
	  {
	    match(lpar);
	    Expr *arg2 = expression();
	    match(rpar);
	    match(lpar);
	    Expr *arg1 = expression(); 
	    match(becomes); 
	    res = new Epenalty(larg, arg1, arg2, expression(), &an);
	    match(rpar);
	  }
	break;
      }
        
    case convec:  // syntax foX = convec(<mesh,>f,u,v,dt), f,u,v are array functions
      {
        nextSym();
        match(lpar); // optional mesh name
        Iden* idm;
        if(curIden->type==Iden::maillage)
          {     idm = curIden; match(iden); match(comma);}
        else idm = curMesh;
        Iden* idu; 
        if (curIden->type == Iden::ftableau)  
          idu = curIden;
        else throw(ErrorCompile(" Array function expected "));
        match(iden); match(comma); 
        Iden* idv; 
        if (curIden->type == Iden::ftableau)  
          idv = curIden;
        else throw(ErrorCompile(" Array function expected "));
        match(iden); match(comma); 
        Expr* e = expression(); // time step
        match(comma); 
        Iden* idf; 
        if (curIden->type == Iden::ftableau)  
          idf = curIden;
        else throw(ErrorCompile(" Array function expected "));
        match(iden);
        match(rpar);
        res = new Econvec2(idm,idf,idu,idv,e,&an);
        break;
      }

    default: 
      char Text [100];
      strcpy(Text," Expected symbol -> ");
      strcat(Text,SymbolName[curSym]);
      throw(ErrorCompile(Text));
    }
  if (curSym == power)
    {nextSym(); res = new EF2(Pow, "power", res, facteur());}
  return res;
}

Expr * Analyseur::terme()
{
  Expr *res = facteur(); Symbol t;
  while (((t = curSym) == star) || (t == slash))
    {
      nextSym();
      if (t == star)
        res = new EF2(Mul, "produit", res, facteur());
      else
        res = new EF2(Div, "quotient", res, facteur());
    }
  return res;
}

Expr* Analyseur::exprarith()
{
  Expr* res; Symbol t;
  switch(curSym)
    {
    case _plus: nextSym(); break;
    case _minus:
      nextSym();
      res = new EF (Chs, "moins unaire", terme());
      break;
    default: res = terme();
    }
  while (((t = curSym) == _plus) || (t == _minus))
    {
      nextSym();
      if (t == _plus)
        res = new EF2(Add, "somme", res, terme());
      else
        res = new EF2(Sub, "difference", res, terme());
    }
  return res;
}

Expr* Analyseur::exprcomp()
{
  Expr* res = exprarith(); Symbol t;
  while(((t = curSym) >= lt) && (t <= eq))
    {
      nextSym();
      switch(t)
        {
        case lt: res = new EF2(Lt, "inferieur", res, exprarith()); break;
        case le: res = new EF2(Le, "inferieur ou egal", res, exprarith()); break;
        case gt: res = new EF2(Gt, "superieur", res, exprarith()); break;
        case ge: res = new EF2(Ge, "superieur ou egal", res, exprarith()); break;
        case eq: res = new EF2(Eq, "egal", res, exprarith()); break;
        case ne: res = new EF2(Ne, "different", res, exprarith()); break;
        }
    }
  return res;
}
EChaine * Analyseur::expchaine(char * errmsg)
{ 
  EChaine * r=0;
  EChaine **e=&r;
  if (curSym==chaine) // a expchaine  begin with a string
    do {
      if (curSym==chaine)
        *e =new EChaine(curChaine),nextSym();
      else
        *e =new EChaine(expression());
      e=& (**e).next;     
    } while (IsSym(tilde));
  if (r==0 && errmsg) throw(ErrorCompile(errmsg));
  return r;
}
Expr* Analyseur::expression()
{
  Expr *res, *l=0, *ll=0, *lll=0,*l4=0,*l5=0,*l6=0,*l7=0;
  Iden *id;
  Symbol theSym;
  theSym=curSym;
  {
    res = exprcomp(); 
    Symbol t;
    while(((t = curSym) == et) || (t == ou))
      {
        nextSym();
        if (t == et)
          res = new EF2(Et, "et logique", res, terme());
        else
          res = new EF2(Ou, "ou logique", res, terme());
      }
    return res;
  }
}


void Analyseur::readCoef(Expr **exp, int *addmul, int nOp)
{                                       // fills exp[nOp*i+ op] and the addmul = 1,2,3,4,5,6 corresponding to
  int addsym, k;        //  1 if op not used, 2 if no exp but op used, 3=+*,4=-*,5=+/,6=-/
  // exp[nOp*i+ op] is the exp in from of op(k) applied to unknown i
  for(int i=0;i<edp->n*nOp;i++) 
    {   
      exp[i]=0; 
      addmul[i]=-1; // operator is not used in expression
    }
  do{   
    if(curSym==becomes)
      {         nextSym(); 
	exp[0] = expression();                            //rhs of PDE or bdy cond 
	addmul[0] = 4;                                            // it is the rhs: no sign on exp
      } else{
      addsym=0;                                                       // sign in front of operator
      if((curSym==_plus)||(curSym==_minus))     
	{ addsym = (curSym==_minus); nextSym();} // addsym = 0 if + or nothing and 1 if -
      switch(curSym)
	{
	case iden:    k=1;    break;  
	case dx:              k=2;    break;  
	case dy:              k=3;    break;  
	case dnu:             k=4;    break;  
	case laplace: k=5;    break;
	case dxx:     k=6;    break;
	case dxy:     k=7;    break;
	case dyx:     k=8;    break;
	case dyy:     k=9;    break;
	default: throw(ErrorCompile(" Operator expected"));
	}                      
      if(k>1) { nextSym(); match(lpar);}
      int i;
      for(i=0;i<edp->n;i++) 
	if( edp->id[i]==curIden) break;
                
      if(i==edp->n) 
	throw(ErrorCompile(" operator on a VAR which is not in LIST of solve(LIST) "));
      match(iden);
      if(k>1) match(rpar);
      if(curSym==star||curSym==slash)
	{
	  addmul[nOp*i+k] = addsym+2*(1+(curSym==star)); // -/:1+2, +/:2, -*:1+4, +*:4 
	  nextSym(); 
	  exp[nOp*i+k] = terme();
	} else  addmul[nOp*i+k] = addsym;
    }
  }while(curSym!=semicolon && curSym!=rbrace);
  IsSym(semicolon);
}

MeshExpr * Analyseur::genmesh(void)
{
  MeshExpr *res=0;
  Expr *l=0;
  Iden *id = ((curSym==iden) ? curIden : 0);
  Symbol theSym=curSym;
  nextSym();
  match(lpar);
  EChaine * ffile=  (theSym!=readmesh)?expchaine():0;
  if (ffile) match(comma);
  if (theSym==buildmesh)
    { //example mesh th = buildmesh(aa(n) + bb(m));      
      l = expression();
      res = new Ebuildmesh(l,&an,ffile);
      match(rpar);
    }
  else if (theSym==readmesh)
    {  
      Expr *rr=0;
      EChaine * ffilemesh = expchaine("file name expected");
      if (IsSym(comma)) rr = expression();               
      res = new Ireadmesh(ffilemesh,&an,ffile,rr);   
      match(rpar);
    }
  else if (id && id->type == Iden::maillage) 
    {
      Expr  *e = expression();
      match(comma);
      Expr  *b = expression();
      match(rpar);
      res =  new Etruncmesh(id,e,b,&an,ffile);            
    }
  else 
    {
      // set de background mesh 
      if(curIden->type==Iden::maillage)
        { id = curIden; 
	  cout << " - " << id -> fg << "  " << id -> fg -> Th << endl;
	  match(iden); match(comma);}
      else id = curMesh;
      
      if (theSym==adaptmesh)
        return  new Edaptmesh(*this,id,&an,ffile);
      else if(theSym==movemesh)
        {                       
          l = expression();
          match(comma);
          res = new Emovemesh(id,l,expression(),&an,ffile);
          match(rpar);
        }
      else
        throw(ErrorCompile("movemesh, adaptmesh or buildmesh keywork expected"));
    }
  return res;
}

Instr* Analyseur::instruction(void)
{
  Expr *l=0, *l1=0,*l2=0,*l3=0,*l4=0,*l5=0,*l6=0,*l7=0;
  Instr *i1=0, *i2 = 0, *res = rien;
  an.local = 0;
  switch(curSym)
    {
    case fdecl2: lisFonction2 (); break;
    case mesh: 
      { 
        nextSym();
        Iden* id = curIden;
        nextSym();
        if (curIden->type == Iden::inconnu || curIden->type == Iden::maillage ) 
          {     
            if (curIden->type == Iden::inconnu)
              { // init --
                curIden->newVar();
                curIden->type = Iden::maillage;
                //id->fg = new Grid();
                id->fg = 0; // the grid is no construct modif FH                  
              }
            if (IsSym(semicolon))
              res = new MeshCode(id,0, &an);
            else if(IsSym(becomes))
              res = new MeshCode(id,genmesh(), &an);
            else 
              throw(ErrorCompile(" = or ; expected "));
            
            curMesh = id;       // curMesh is the name after mesh=
          } 
        else throw(ErrorCompile(" New name or mesh name expected "));
      } break;
      
    case border: 
      lisBorder(); break;
      
    case si:
      nextSym(); l = expression(); match(alors); i1 = instruction();
      if (curSym == autrement) {nextSym(); i2 = instruction();}
      res = new IIf(l, i1, i2); break;
      
    case pour:
      { 
        nextSym();
        Iden* id = curIden;
        if (curIden->type == Iden::inconnu) 
          {
            curIden->newVar(); 
            curIden->type = Iden::variable;
          }
        match(iden); match(becomes);
        l = expression(); match(jusqua); l1 = expression();
        match(faire); i1 = instruction();
        res = new Ifor(id,l,l1,i1);
        break;
      }
      
    case iden:
      {
        Iden *id = curIden;
        nextSym(); 
        if (id->type == Iden::maillage)
          {
            match(iden); 
            Iden *idp0 = curIden;
            if(idp0->type != Iden::ftableauP0) 
	      throw(ErrorCompile(" expected variable name of type femp0"));
	    IsSym(colon); //heat colon 
            match(becomes); 
            res = new IarrayP0(idp0,expression(),&an,id);                 
          }
        else if (id->type == Iden::ftableau)// it is an old array variable
	  {
	    IsSym(colon); //heat colon 
	    match(becomes) ;
	    res = new Iarray(id, expression(),&an);
	  } 
	else if (id->type == Iden::ftableauP0 )  
	  {
	    IsSym(colon); //heat colon 
	    match(becomes) ;
	    res = new IarrayP0(id, expression(),&an);
	  } 
	else if (id->type == Iden::pint)  
	  {
	    IsSym(colon); //heat colon 
	    match(becomes) ;
	    res = new Ibecomes(id, expression(),&an);
	  } 
	else if (id->type == Iden::pfloat)  
	  {
	    IsSym(colon); //heat colon 
	    match(becomes) ;
	    res = new Ibecomes(id, expression(),&an);
	  } 
	else 
	  {
	    if(IsSym(colon)) // it is a scalar variable
	      {
		if (id->type == Iden::inconnu) 
		  id->newVar();
		id->type = Iden::variable;
	      }
	    if (id != curFunc && id->type != Iden::variable 
		&& id->type != Iden::maillage && id->type != Iden::matrix) // it is a function
	      {     
		lisFonctionN (id); 
		break;
	      }
	    if (id->type == Iden::inconnu) 
	      id->newVar();
	    if (id != curFunc && id->type != Iden::variable 
		&& id->type != Iden::maillage && id->type != Iden::matrix)
	      throw(ErrorCompile(" Unexpected variable name"));
	    match(becomes); 
	    res = new Ibecomes(id, expression(),&an);
	  }
      } 
      break;
     
    case subroutine:  // subroutines' syntax  subroutine sub(arg1,arg2...) <instruction>
      {
        nextSym(); 
        Iden *id = curIden;
        match(iden); 
        if (id->type != Iden::inconnu) throw(ErrorCompile(" New name expected "));
        curFunc = id; 
        if (curFunc->type != Iden::inconnu ) 
          throw(ErrorCompile(" New name expected"));
        curFunc->newVar();
        Iden * arg;
        ArgList larg = 0;
        int haspar=0;
        if(curSym==lpar)
          {     haspar = 1; nextSym();
	    while (curSym != rpar) 
	      { 
		arg = curIden; match(iden); 
		if (curFunc == arg) throw(ErrorCompile(" Name conflict"));
		if (arg->type == Iden::inconnu) arg->newVar();
		if (arg->type != Iden::variable) 
		  throw(ErrorCompile(" Variable name expected "));
		if(curSym!=rpar) match(comma);
		larg = new ArgNode (arg, larg);
	      }
	    match(rpar);
          }
        id->type = Iden::fonctionN;
        id->fn = (void *) new CFN(curFunc, larg, instruction());

        break;
      }
    case fdecl:
      {
        nextSym(); 
        Iden *id = curIden;
        match(iden); 
        if (id->type != Iden::inconnu) throw(ErrorCompile(" New name expected "));
        lisFonctionN (id); 
        break;
       
      }
     
    case number:
      {
        nextSym(); 
        Iden *id = curIden;
        match(iden);
        if (id->type == Iden::inconnu)
          {
            id->newVar();
            id->type = Iden::variable;
          }
        else if(id->type!=Iden::variable) 
          throw(ErrorCompile(" Cannot redeclare this as a number "));
        match(becomes); 
        res = new Ibecomes(id, expression(),&an);
        break;
      }
    case array:
    case arrayP0:
    case arrayP1:
      {   
        Iden::TypeIden typearray = (curSym == arrayP0) ? Iden::ftableauP0 : Iden::ftableau; 
        Iden* idm;
        nextSym();
        if(curSym==lpar)
          {
            nextSym(); 
            idm = curIden; match(iden);
            match(rpar);
          }
        else idm = curMesh;
        if(!idm || idm->type!=Iden::maillage) 
	  throw(ErrorCompile(" Mesh does not exists; can't create array "));
        if (curIden->type != Iden::inconnu && (curIden->type != typearray))
          throw(ErrorCompile(" Only arrays can be redefined as same arrays  "));
        Iden *id = curIden;
        match(iden); 
        match(becomes) ;
        if (id->type == Iden::inconnu)
          {     id->type = typearray;
	    if (typearray== Iden::ftableau) 
	      id->ft =new CTab (idm,&an);
	    else if (typearray== Iden::ftableauP0) 
	      id->fP0 =new CTabP0 (idm,&an);
	    else throw(ErrorCompile(" Error compiler typearray unknow   "));   
          }
        

        
        if (typearray== Iden::ftableau) {
          id->ft->gridid = idm; 
          res = new Iarray(id, expression(),&an);}
        else {
          id->fP0->gridid = idm; 
          res = new IarrayP0(id, expression(),&an);}
        break;
      }
      
    case lbrace:
      nextSym ();
      while (curSym != rbrace)
        {
          i1 = instruction();
          if (res == rien) res = i1; else if (i1 != rien) res = new IS(res, i1);
        }
      match(rbrace); break;
      
    case RR2:
      an.x->type == Iden::inconnu; //because "table.remove(an.x->nom)" doesn' work
      an.y->type == Iden::inconnu;
      nextSym(); match(lpar);
      an.x = curIden;       
      if (an.x->type == Iden::inconnu) an.x->newVar();
      if (an.x->type != Iden::variable) throw(ErrorCompile(" Variable name expected "));
      match(iden); match(comma);
      an.y = curIden;  
      if (an.y->type == Iden::inconnu) an.y->newVar();
      if (an.y->type != Iden::variable) throw(ErrorCompile(" Variable name expected "));
      if (an.x == an.y) throw(ErrorCompile(" Name conflict "));
      match(iden);
      if(curSym!=rpar)
        {
          table.remove(an.ng->nom);
          match(comma);
          an.ng = curIden; 
          if (an.ng->type == Iden::inconnu)  an.ng->newVar(); 
          if (an.ng->type != Iden::variable) 
	    throw(ErrorCompile("  Variable name expected "));
          if ((an.x == an.ng) || (an.y == an.ng)) 
	    throw(ErrorCompile(" Name conflict "));
          match(iden);
        }
      match(rpar);
      break;
      
    case plot:
    case plotP0:
      {
        Symbol sym =  curSym;
        nextSym();
        match(lpar);
        EChaine * ffile= expchaine();
        if (ffile) match(comma);
        OnNode * larg=0,*lend=0;
        // 
        if (curSym!=iden || curIden->type != Iden::maillage)
          lend=larg = new OnNode (curMesh, larg);       
        // new code             
        do //(curSym != rpar )//curSym==iden && curIden->type == Iden::courbe)
          {     
            OnNode * ll ;
            if (curSym==iden && curIden->type == Iden::maillage)                
              ll = new OnNode (curIden, 0),match(iden);
            else 
              ll = new OnNode (expression(), 0);
            // link and end of list 
            if (larg) lend->n  = ll;
            else        larg=ll;
            lend=ll;        
          }     while (IsSym(comma)); 
        if (sym == plot)
          res = new Iplot(larg,&an,ffile);  
        else
          res =  new IplotP0(larg,&an,ffile);       
        match(rpar);
        break;
      }   
    case plot3d:
      {
        nextSym();
        match(lpar);
        EChaine * ffile= expchaine();
        if (ffile) match(comma);
        Iden* id;
        if(curIden->type==Iden::maillage)
          {  id = curIden; match(iden); match(comma);}
        else id = curMesh;
        l = expression();
        res = new Iplot3d(id,l,&an,ffile);
        match(rpar);
        break;
      }
    case print:
    case append:
      {
        int ondisk = curSym==append;
        nextSym();
        match(lpar);
        EChaine *fname;
        
        if(ondisk) 
          {fname = expchaine(" Missing file name ");match(comma);}
        else fname = 0;
        /*      if(fname)       match(comma);
                if(curSym != rpar)
                l = expression();
                else l=0;
                res = new Iprint(ondisk, fname,l,an.x->storage, an.y->storage, an.ng->storage);*/
        { // new code   to  can do print sqdqs,dsqdsqd,dsqdsqd,sqd,sqd,sq,"sdfds"
          EChaine * r=0;
          EChaine **e=&r;
          do {
            if (curSym==chaine)
              *e =new EChaine(curChaine),nextSym();
            else
              *e =new EChaine(expression());
            e=& (**e).next; 
            if( curSym ==  comma) // put a space
              {
                *e =new EChaine(" ");
                e=& (**e).next; 
              }  
          } while (IsSym(comma) || IsSym(tilde) );      
          match(rpar);
          res = new Iprint(ondisk, fname,r,an.x->storage, an.y->storage, an.ng->storage);

        }
        // end new code 
        break;
      }   
      /*    case helmholtz:
            {   
            nextSym();
            match(lpar);
            Iden *id1, *id2;
            if(curIden->type==Iden::maillage)
            {id1 = curIden; match(iden); match(comma);}
            else id1 = curMesh;
            CTab *f ;
            if (curIden->type == Iden::inconnu)
            {
            curIden->type = Iden::ftableau; 
            f = new CTab (id1, &an);
            } 
            else
            if (curIden->type == Iden::ftableau) 
            {
            f = curIden->ft; 
            f->g = (Grid*)id1->fn; 
            }
            else erreur("Can't redefine that function");
            curIden->ft = f;
            match(iden); match(comma);
            l1 = expression(); match(comma);
            l2 = expression(); match(comma);
            l3 = expression(); match(comma);
            l4 = expression(); match(comma);
            l5 = expression(); match(comma);
            l6 = expression(); match(comma);
            l7 = expression();
            int factorize = 0; 
            if(curSym==comma)
            { 
            nextSym(); 
            if(curSym==_plus) 
            {
            nextSym(); factorize = 1;
            if ((curIden->type != Iden::inconnu)&&(curIden->type != Iden::matrix))
            erreur(" New or matrix identifier expected ");
            curIden->type = Iden::matrix; 
            } 
            else if(curSym==_minus) 
            {
            nextSym(); factorize = -1;
            if (curIden->type != Iden::matrix)
            erreur(" This identifier does not refer to a known matrix ");
            curIden->type = Iden::matrix; 
            } 
            else erreur(" + or - expected ");
            id2 = curIden; 
            match(iden);
            }  
            else id2 = 0;
            match(rpar);
            res = new Ihelmholtz(id1,id2,f,l1,l2,l3,l4,l5,l6,l7,&an, factorize);
            break;
            }
      */     
    case solve:         // syntax: solve(<mesh,> u1 <,u2...,uN>) <with A(0)> instr where instr is
      {
        nextSym();      // { pde(u1) u1 * exp + laplace(u1) * exp + dx(u1) * exp... = exp;
        match(lpar);//   on(a1 <,a2,..,aM>) u1 * exp + dnu(u1) * exp + dx(u1) * exp... = exp};
        Iden *idmatrix, *idmesh ;                               // matrix and mesh name
        Expr *factorize;
        if(curIden->type==Iden::maillage)               // finds if mesh name is specified
          {
            idmesh = curIden; 
            match(iden); 
          } 
        else idmesh = curMesh;
        Iden* oldCurMesh = curMesh;
        curMesh = idmesh;                                               // while in solve use this mesh
       
        edp = new EDP;
        edp->f = new CTabpt[edpmax];                    
        edp->id = new Idenpt[edpmax];
        edp->n = 0;
        edp->g = curMesh->fg;
       
        do
          {
            IsSym(comma);
            match(iden);                                                // has found u1,u2...
            edp->id[edp->n] = curIden;
            edp->f[edp->n] = curIden->iden2CTab(curMesh, &an);
            curIden->f2 = edp->f[edp->n++];
          } while(curSym!=rpar);
        match(rpar);
       
        if(curSym == avec)                                              // specifies the matrix to use
          {
            nextSym();
            match(iden);
            if ((curIden->type != Iden::inconnu)&&(curIden->type != Iden::matrix))
              throw(ErrorCompile(" New or matrix identifier expected "));
            curIden->type = Iden::matrix; 
            idmatrix = curIden;
            match(lpar);
            factorize = expression();
            match(rpar);
          }
        else
          { factorize = 0; idmatrix=0;}
       
        switch(edp->n)
          {
          case 1: res = new Isolve<float,float,1>(idmesh,idmatrix,factorize, edp, instruction(),&an); break;
          case 2: res = new Isolve<MatN<float,2>, VectN<float,2>,2>(idmesh,idmatrix,factorize, edp, instruction(),&an); break;
          case 3: res = new Isolve<MatN<float,3>, VectN<float,3>,3>(idmesh,idmatrix,factorize, edp, instruction(),&an); break;
          case 4: res = new Isolve<MatN<float,4>, VectN<float,4>,4>(idmesh,idmatrix,factorize, edp, instruction(),&an); break;
          case 5: res = new Isolve<MatN<float,5>, VectN<float,5>,5>(idmesh,idmatrix,factorize, edp, instruction(),&an); break;
          default: throw(ErrorCompile(" Sorry no more than 5 PDEs"));
          }
        curMesh = oldCurMesh;
        break;
      }         
     
    case pde:           //syntax: pde(u1) readcoef; < on(u1) readcoef...>
      {
        if(edp->n==0) throw(ErrorCompile(" 'pde' should be called with a 'solve' block"));
        nextSym(); 
        match(lpar);
        match(iden); 
        for(edp->j=0;edp->j<edp->n;edp->j++) 
          {  if( edp->id[edp->j]==curIden) break;}
        if(edp->j == edp->n) 
	  throw(ErrorCompile(" VAR is not part of LIST in 'solve(LIST) pde(VAR)' "));
        match(rpar); 
        Expr* *exp;     exp = new Exprpt[10*edp->n];
        int* addmul;    addmul = new int[10*edp->n];
        readCoef(exp,addmul,10);
        res = new Ipde(edp,edp->j, exp, addmul, instruction(), &an);
        break;                  
      }
      
    case sur:
      {
        if(edp->n==0 || edp->j>=edp->n || edp->j <0) 
          throw(ErrorCompile(" 'on' should be called with a 'pde' block"));
        nextSym();
        match(lpar);
        OnList larg = 0;
        do 
          { 
            if (curSym==iden && curIden->type == Iden::courbe)          
              larg = new OnNode (curIden, larg),match(iden);
            else 
              larg = new OnNode (expression(), larg);               
          } while (IsSym(comma));
        match(rpar);
        Expr* *exp; exp = new Exprpt[5*edp->n];
        int* addmul;  addmul = new int[5*edp->n];
        readCoef(exp,addmul,5);
        res = new Isur(edp,edp->j, larg, exp,addmul, &an); 
        break;
      }

    case varsolve:                                              //varsolve(th,oldnew) with aa(u1,v1,u2,v2) = expression
      {
        int nedp = 0;
        Iden *idmesh = curMesh, *oldCurMesh = curMesh;
        l=0;
        nextSym();
        if(curSym==lpar)
          {     
            nextSym();
            if(curIden->type==Iden::maillage)           // finds if th is present
              {
                idmesh = curIden; 
                match(iden); 
                if(curSym!=rpar) match(comma);
              }
            if(curSym!=rpar)
              l = expression();                                         // oldnew                       
            match(rpar);
          }
        curMesh = idmesh;                                                       // while in varsolve
        match(iden);                                                    // finds aa
        if (curIden->type == Iden::inconnu) 
          {  
            curIden->newVar();
            curIden->type = Iden::matrix; 
          }
        else if( curIden->type != Iden::matrix )  
          throw(ErrorCompile(" This name has been used before in a different context "));
        Iden *id0 = curIden;
        match(lpar);
        typedef CTab* CTabpt;
        CTab** f1; f1 = new CTabpt[edpmax];
        CTab** f2; f2= new CTabpt[edpmax];
                 
        do
          {
            IsSym(comma);
            match(iden);                                   // has found u
            f1[nedp] = curIden->iden2CTab(curMesh, &an);
            curIden->f2 = f1[nedp];
            match(comma);
            match(iden);                                // has found v
            f2[nedp] = curIden->iden2CTab(curMesh, &an);
            curIden->f2 = f2[nedp++];
          } while(curSym!=rpar);
                 
        match(rpar);  match(avec);
        switch(nedp)
          {
          case 1: res = new Ivarsolve<1>(id0,idmesh,nedp,f1, f2,l, instruction(), &an); break;
          case 2: res = new Ivarsolve<2>(id0,idmesh,nedp,f1, f2,l, instruction(), &an); break;
          case 3: res = new Ivarsolve<3>(id0,idmesh,nedp,f1, f2,l, instruction(), &an); break;
          }
        break;
      }         
    case assemble:
      {
        nextSym();
        match(lpar);
        Iden* id;
        if(curIden->type==Iden::maillage)
          {     id = curIden; match(iden); match(comma);}
        else id = curMesh;
        CTab *f1, *f2;
        if (curIden->type == Iden::inconnu) 
          {
            curIden->newVar();
            curIden->type = Iden::ftableau;
            f1 = new CTab (curMesh,&an);
            curIden->f2 = f1;
          } else if (curIden->type == Iden::ftableau)  f1 = curIden->ft;
        else throw(ErrorCompile(" Function of space variables expected "));
        match(iden); match(rpar);
        if (curIden->type == Iden::inconnu) 
          {
            curIden->newVar();
            curIden->type = Iden::ftableau;
            f2 = new CTab (curMesh,&an);
            curIden->f2 = f2;
          } else if (curIden->type == Iden::ftableau)  f2 = curIden->ft;
        else throw(ErrorCompile(" Function of space variables expected "));
        match(iden); match(becomes);
        res = new Iassemble(id,f1,f2,expression(), &an);
        break;
      }
  
    case convec:
      {
        nextSym();
        match(lpar); // optional mesh name
        Iden* id;
        CTab *f1, *f2; 
        if(curIden->type==Iden::maillage)
          {     id = curIden; match(iden); match(comma);}
        else id = curMesh;
        l1 = expression(); // u-velocity
        match(comma);
        l2 = expression(); // v-velocity
        match(comma);
        l3 = expression(); // time step
        match(comma);
        if (curIden->type == Iden::inconnu) 
          {
            curIden->newVar();
            curIden->type = Iden::ftableau;
            f1 = new CTab (curMesh,&an);
            curIden->f2 = f1;
          } else if (curIden->type == Iden::ftableau)  f1 = curIden->ft;
        else throw(ErrorCompile(" Function of space variables expected "));
        match(iden); match(comma); 
        l4 = expression();
        if(curSym==comma) 
          {     nextSym();
	    if (curIden->type == Iden::inconnu) 
	      {
		curIden->newVar();
		curIden->type = Iden::ftableau;
		f2 = new CTab (curMesh,&an);
		curIden->f2 = f2;
	      } else if (curIden->type == Iden::ftableau)  f2 = curIden->ft;
	    else throw(ErrorCompile(" Function of space variables expected "));
	    match(iden);match(comma); 
	    l5 = expression();
          } else f2 = 0;
        match(rpar);
        res = new Iconvec(id,f1,f2,l1,l2,l3,l4,l5,&an);
        break;
      }
      
    case savemesh:
      //case readmesh:
      { Symbol theSym = curSym;
	Iden* id;
	nextSym();
	match(lpar);
	EChaine * fname= expchaine("Missing mesh file name");
	if(curSym==comma)
	  {
	    nextSym();
	    if ((theSym==savemesh)&&(curIden->type != Iden::maillage)) 
	      throw(ErrorCompile(" Mesh name expected "));
	    if((theSym==readmesh)&&(curIden->type != Iden::inconnu)) 
	      throw(ErrorCompile(" New Mesh name expected "));
	    if  (curIden->type == Iden::inconnu)
	      {   
		curIden->newVar(); 
		curIden->type = Iden::maillage; 
	      }
	    id = curIden;
	    match(iden);
	    if(theSym==readmesh) 
	      {
		id->fg = 0; // new Grid();modif FH
		curMesh = id;
	      }
	  }
	else
	  {
	    if(theSym==readmesh) throw(ErrorCompile(" Mesh name expected "));
	    else if(curMesh) id=curMesh; else throw(ErrorCompile("Mesh does not exist"));
	  }
	match(rpar);
	if (theSym==savemesh) res =  new Isavemesh(fname,id,&an);
	// else  res = new Ireadmesh(fname,id,&an); 
      
      }
      break;
      
    case save:
    case lire:
      {
        Symbol theSym = curSym;
        nextSym();
        match(lpar);
        EChaine * fname=expchaine("Missing file name");
        match(comma);
        Iden* id;
        if(curIden->type==Iden::maillage)
          {     id = curIden; match(iden); match(comma);}
        else id = curMesh;
        if (theSym==save)
          {
            l = expression();
            res = new Isave(id,fname,l,&an);
          }
        else 
          {
            CTab *f ;
            if (curIden->type == Iden::inconnu)
              {
                curIden->type = Iden::ftableau; 
                f = new CTab (id,&an);
                curIden->ft = f;
              }
            else if (curIden->type == Iden::ftableau) 
              {
                f = curIden->ft; 
                f->gridid = id; // change the id (bof?? modif FH)
              }
            else throw(ErrorCompile(" Can't redefine that function"));
            res = new Iread(id,fname,f);
            curIden->f2 = f;
            nextSym();
          }
        match(rpar);
      }
      break;
      
    case semicolon: nextSym();
    case ::end: break;
    case progexit: exit(0);
    default:
      throw(ErrorCompile(" Unexpected symbol beginning an instruction "));
    }
  return res;
}

class CF : public CVirt
{ // Fonction d'une variable locale
  Iden *f, *a; Instr *i;
public:
  CF (Iden *ff, Iden *aa, Instr *ii) : f(ff), a(aa), i(ii) {}
  float operator() (float x)
  {
    float o = *a->storage;
    *a->storage = x; i->execute(); *a->storage = o;
    return *f->storage;
  }
  CVirt *de (Iden *x) {return new CF (f, a, i->d(x));}
  CVirt *de () {return de (a);}
};

void Analyseur::lisFonctionN(Iden* id)
{
  curFunc = id; 
  if (curFunc->type != Iden::inconnu ) 
    throw(ErrorCompile(" New name expected"));
  curFunc->newVar();
  //  switch (curSym)
  //    {
  //    case lpar:{
  Iden * arg;
  ArgList larg = 0;
  int haspar=0;
  if(curSym==lpar)
    { haspar = 1; nextSym();
      while (curSym != rpar) 
	{ 
	  arg = curIden; match(iden); 
	  if (curFunc == arg) throw(ErrorCompile(" Name conflict"));
	  if (arg->type == Iden::inconnu) arg->newVar();
	  if (arg->type != Iden::variable) throw(ErrorCompile(" Variable name expected "));
	  if(curSym!=rpar) match(comma);
	  larg = new ArgNode (arg, larg);
	}
      match(rpar);
    }
  match(becomes); 
  Instr* res = new Ibecomes(id, expression(),&an);
  id->type = Iden::fonctionN;
  id->fn = (void *) new CFN(curFunc, larg, res);
  /*      break;
          }
          case becomes:
          {match (becomes); match (der); match (lpar);
          Iden * prim = curIden; match (iden); match (comma);
          if (prim->type != Iden::fonction) erreur(" Function name expected ");
          Iden * vari = curIden; match (iden); match (rpar);
          if (vari->type != Iden::variable) erreur(" Variable name expected ");
          curFunc->type = Iden::fonction;
          curFunc->f = ((CF *)prim->f)->de (vari);
          break;}
          default:
          erreur (" Function definition expected ");
          }
  */}

class CF2 : public CVirt2
{ // Fonction de deux variables locales
  Iden *f, *a1, *a2; Instr *i;
public:
  CF2 (Iden *ff, Iden *aa1, Iden *aa2, Instr *ii)
    : f(ff), a1(aa1), a2(aa2), i(ii) {}
  void setInstruction(Instr *ii) {i = ii;}
  float operator() (float x, float y)
  {
    float o1 = *a1->storage, o2 = *a2->storage;
    *a1->storage = x; *a2->storage = y; i->execute();
    *a1->storage = o1; *a2->storage = o2;
    return *f->storage;
  }
  CVirt2 *de (Iden *x) {return new CF2 (f, a1, a2, i->d (x));}
  CVirt2 *de1 () {return de (a1);}
  CVirt2 *de2 () {return de (a2);}  
};

void Analyseur::lisFonction2()
{
  match(fdecl2);
  curFunc = curIden;
  match(iden);
  if (curFunc->type != Iden::inconnu) throw(ErrorCompile(" Function name expected "));
  curFunc->newVar();
  switch (curSym)
    {
    case lpar:
      {match(lpar);
	Iden *arg1 = curIden; match(iden); match(comma);
	Iden *arg2 = curIden; match(iden); match(rpar);
	if ((curFunc == arg1) || (curFunc == arg2) || (arg1 == arg2)) 
	  throw(ErrorCompile(" Name conflict"));
	if (arg1->type == Iden::inconnu) arg1->newVar();
	if (arg1->type != Iden::variable) throw(ErrorCompile(" Variable name expected "));
	if (arg2->type == Iden::inconnu) arg2->newVar();
	if (arg2->type != Iden::variable) throw(ErrorCompile(" Variable name expected "));
	curFunc->type = Iden::fonction2;
	curFunc->f2 = new CF2(curFunc, arg1, arg2, instruction());
	break;}
    case becomes:
      {match (becomes); match (der); match (lpar);
	Iden * prim = curIden; match (iden); match (comma);
	if (prim->type != Iden::fonction) throw(ErrorCompile(" Function name expected "));
	Iden * vari = curIden; match (iden); match (rpar);
	if (vari->type != Iden::variable) throw(ErrorCompile(" Variable name expected "));
	curFunc->type = Iden::fonction2;
	curFunc->f2 = ((CF2 *)prim->f2)->de (vari);
	break;}
    default:
      throw(ErrorCompile(" Function definition expected "));
    }
}

void Analyseur::FindArgs(const NameArg & na)
{  // F Hecht ofr writing  toto=exp tyty=exp ... 
  // where toto , tyty are defined in na.ListOfArg[i]
  int k;
  for ( int i=0;i<na.NbOfArg;i++)
    na.arg[i]=0;// MISE A       zero for all arg
  //  cout << " curIsAlphaNum " << curIsAlphaNum << " " << curIden->nom <<  endl;
  while (curIsAlphaNum) 
    {
      k=-1;    
      cout << "\t\t" << "  Args = " << curIden->nom ;
      for ( int i=0;i<na.NbOfArg;i++) 
        {
          if(!strcmp(curIden->nom,na.ListOfArg[i])) 
            { k=i;break;}
        }
      if(k<0) {
	//cerr << " Argument inconnu " << curIden->nom << endl;
	char * s = " Argument inconnu ";
	strcat(s,curIden->nom);
	throw(ErrorCompile("Argument inconnu"));    
      }
      if (na.arg[k])  throw(ErrorCompile("FindArgs: Arg allready set"));
      nextSym();
      match(becomes);
      na.arg[k] = expression();
      if(!IsSym(comma)) break;
    }
  //  cout << " ------------------- +++++++++++++++ " <<endl;
}
int  Analyseur::GetExprs (Expr ** e,int  n)
{  // F Hecht
  // to store a list of exp0,exp1, ... , expk  in the array e 
  // with k < n
  int i=0;
  
  while (e[i++]=expression(),IsSym(comma) )
    if (i>=n) throw(ErrorCompile(" GetExprs Too Much arg"));
  return i;
}  

void Analyseur::lisBorder()
{ // border id( id1=exp , exp)
  match(border);
  Iden* id = curIden;
  id = curIden;
  match(iden);
  if (id->type != Iden::inconnu) throw(ErrorCompile(" New name expected"));
  id->newVar();
  id->type = Iden::courbe;
  match(lpar);
  Iden *t = curIden; match(iden); match(becomes);
  if (t->type == Iden::inconnu) t->newVar();
  if (t->type != Iden::variable) throw(ErrorCompile(" Variable name expected "));
  Expr *start = expression(); match(comma);
  Expr *stop = expression(); 
  match(rpar);
  *an.ng->storage = (float)(an.bdyLabel);  //most likely useless
  id -> fn = (void *) new IB(an.x -> storage,an.y -> storage,
			     an.ng -> storage, t -> storage,start, stop, 
			     instruction(),an.bdyLabel);
  an.bdyLabel++;
}
        
// inline Fonction GetFonction (IdenTable& t, const char *s) {return t.find (s)->f;}

void Analyseur::programme()
{
  if(source)  nextSym();

  while(curSym != ::end)
    instruction()->execute();
}

MeshExpr * gen_scilabmesh(char * s,double * vt, int nbvt, int * tr, int nbtr,
			  Analvar & an)
{
  MeshExpr * Me = NULL;
  EChaine * ec = new EChaine(s);
  Me = new Iscilabmesh(ec,&an,vt,nbvt,tr,nbtr);

  return Me;
}

Instr *  Analyseur::compile() // for link to scilab 
{
  nextSym();
  Instr * res = rien;
  
  while (curSym != ::end)
    {
      Instr * i1 = instruction();
      if (res == rien) res = i1;
      else if (i1 != rien) 
        res = new IS(res, i1);
    }

  return res;
}

#include "Meshio.h"

void  MeshErrorIO(ios& )
{
  throw(ErrorExec("Reading Error"));
}

Analyseur::Analyseur(istream * s) : source (s)
{
  MeshIstreamErrorHandler = MeshErrorIO;

  curMesh = 0; 
  verbosity=0; // for Mesh generation 
  table.enregistre("function", fdecl);
  table.enregistre("function2", fdecl2);
  table.enregistre("derive", der);
  table.enregistre("if", si);
  table.enregistre("then", alors);
  table.enregistre("else", autrement);
  table.enregistre("for", pour);
  table.enregistre("to", jusqua);
  table.enregistre("do", faire);
  table.enregistre("mesh", mesh);
  table.enregistre("border", border);
  table.enregistre("plot", plot); 
  table.enregistre("plot3d", plot3d); 
  table.enregistre("print", print); 
  table.enregistre("append", append); 
  table.enregistre("R2", RR2); 
  table.enregistre("buildmesh", buildmesh);
  table.enregistre("savemesh", savemesh);
  table.enregistre("readmesh", readmesh);
  table.enregistre("movemesh", movemesh);
  table.enregistre("save", save);
  table.enregistre("read", lire);
  table.enregistre("exit", progexit);
  table.enregistre("helmholtz", helmholtz); 
  table.enregistre("adaptmesh", adaptmesh);
  table.enregistre("int", intt); 
  table.enregistre("Id", Heavyside);
  table.enregistre("on", sur);
  table.enregistre("assemble", assemble);
  table.enregistre("varsolve", varsolve);
  table.enregistre("convect", convec);
  table.enregistre("dx", dx);
  table.enregistre("dy", dy);
  table.enregistre("array", array);
  table.enregistre("femp0", arrayP0);
  table.enregistre("femp1", array);
  table.enregistre("plotp0", plotP0);
  table.enregistre("min", aMin);
  table.enregistre("max", aMax);
  table.enregistre("abs", Abs);
  table.enregistre("pow", Pow);
  table.enregistre("exp", Exp);
  table.enregistre("log", Log);
  table.enregistre("sqrt", Sqrt);
  table.enregistre("sin", Sin);
  table.enregistre("cos", Cos);
  table.enregistre("tan", Tan);
  table.enregistre("atan", Atan);
  table.enregistre("atan2", Atan2);
  table.enregistre("one", Fonction::monome(1,1));// x -> x
  table.enregistre("laplace", laplace);
  table.enregistre("dxx",dxx);
  table.enregistre("dxy",dxy);
  table.enregistre("dyx",dyx);
  table.enregistre("dyy", dyy);
  table.enregistre("dnu", dnu);
  table.enregistre("pde",pde);
  table.enregistre("solve",solve);
  table.enregistre("with",avec);
  table.enregistre("number",number);
  table.enregistre("subroutine", subroutine);
  table.enregistre("verbosity", &verbosity);
  // set the derivate of classic function  
  Fonction2 POW  = table["pow"];
  Fonction2 ATAN2  = table["atan2"];
  Fonction  LOG  = table["log"];
  Fonction  EXP  = table["exp"];
  Fonction  SIN  = table["sin"];
  Fonction  COS  = table["cos"];
  Fonction  TAN  = table["tan"];
  Fonction  ATAN = table["atan"];
  Fonction  SQRT = table["sqrt"];   
  POW.setd(Ordonnee*POW(Abscisse,Ordonnee-Fonction2(1)),LOG(Abscisse)*POW);

  Fonction2 Rayon2(Fonction2::monome(1,0, 2) + Fonction2::monome(1,2,0));   
  ATAN2.setd(Ordonnee/Rayon2,-Abscisse/Rayon2);

  SIN.setd(COS);
  COS.setd(-SIN);
  TAN.setd(1./(COS)/(COS));
  EXP.setd(EXP);
  LOG.setd (Fonction::monome (1, -1));
  ATAN.setd (1. / (Fonction (1) + Fonction::monome (1, 2)));
  SQRT.setd(Fonction::monome(0.5,-1)(SQRT));
  an.x = table.find("x");               an.x->newVar();                         // use to label boundaries
  an.y = table.find("y");               an.y->newVar();                         // use to label boundaries
  an.ng = table.find("label");  an.ng->newVar();                        // use to label boundaries
  an.wait = table.find("wait"); an.wait->newVar();                      // wait variable for graphics prompt
  an.nx = table.find("nrmlx");  an.nx->newVar();                        // first component of normal to curve
  an.ny = table.find("nrmly");  an.ny->newVar();                        // second component of normal to curve
  Iden* ppi = table.find("pi"); ppi->newVar();
  *ppi->storage = 4*atan(1.0);                                                          // puts pi in the table
  *an.wait->storage = 1;                                                                        // press mouse or key between graphs
  an.gridxyng=0; // no grid 
  an.ivertex=-1;  
  an.iloc = -1;
  an.trloc = -1;
  an.front = new frontiere();
  an.bdyLabel = 1;
  front = an.front;
}

void frontiere::addSegment(int n1, int n2, int label)
{
  const int increment = 100;
  if (nbs % increment == 0)
    {
      long *ss = new long[2*(nbs + increment)];
      int *nngf = new int[nbs + increment];
      for (int i = 0; i < nbs; i++)
        {       ss[2*i] = s[2*i]; 
	  ss[2*i+1] = s[2*i+1]; 
	  nngf[i] = ngf[i];
        }
      if(s) delete [] s; s = ss;
      if(ngf) delete [] ngf; ngf = nngf;
    }
  ngf[nbs] = label; s[2*nbs] = n1; s[1+2*nbs] = n2; 
  nbs++;
  hh[n1] = 0.66*sqrt( (xy[2*n2] - xy[2*n1])*(xy[2*n2] - xy[2*n1])
                      + (xy[2*n2+1]-xy[2*n1+1])*(xy[2*n2+1]-xy[2*n1+1]));
}

int frontiere::addPoint(float x, float y, int nng)
{
  //  const float epsilon = (float)1e-6;  // precision for 2 identical bdy points
  assert(step && initboundingbox);
  const int increment = 100;
  for (int i = 0; i < nbp; i++)
    if ((fabs(x - xy[2*i]) < epsilon) && (fabs(y - xy[2*i+1]) < epsilon))
      return i+1;
  if (nbp % increment == 0)
    {
      float *xxy = new float[2*(nbp + increment)];
      float *hhh = new float[nbp + increment];
      int *nnng = new int[nbp + increment];
      for (int i = 0; i < nbp; i++)
        {
          xxy[2*i] = xy[2*i];
          xxy[1+2*i] = xy[1+2*i];
          hhh[i] = hh[i];
          nnng[i] = ng[i];
        }
      if(xy) delete[] xy; xy = xxy;
      if(ng) delete[] ng; ng = nnng;
      if(hh) delete[] hh; hh = hhh;
    }
  xy[2*nbp] = x; xy[1+2*nbp] = y; ng[nbp++] = nng;
  return nbp;
}

void erreurt (const char *s)
{
  cerr << s << endl;
  exit(1);
}

void Iarray::execute () 
{
  Analvar save(*an); 
  assert(v->type == Iden::ftableau);
  CTab & ft = *v->ft; 
  Grid * gn = ft.gridid->fg;// the new grid
  gn->check();    
  Grid * go = ft.g; // the old grid 
  assert(gn); // the new grid must exist exect time
  CTab  val(&ft); // creation of the array of value new size 
  if(verbosity>10) 
    cout << " old g " << go << " new g " << gn << " cc=" << val.cc<< endl;
  float xl[3] = {0.,0.,0.};    
  an->gridxyng = gn;    
  for (int i = 0; i < gn->nv; i++)
    {
      int oldlocal = an->local;
      an->setAn(0,gn->v[i].x, gn->v[i].y, gn->v[i].where, xl, i);     
      val[i] = e->eval(); // with  check
      //cout << "x = " << gn->v[i].x << "\t" << "y = " << gn->v[i].y << endl;
      //cout << "val[" << i << "] = " << val[i] << endl;
      an->local  = oldlocal;
    }
 
  if( verbosity > 10)
    cout << "after ft=" << &val << " " <<  v->nom << "=  (min = " << val.min() << " max = " << val.max() << ")  cc=" << val.cc 
	 << "dl " << v->ft <<endl;
  //  delete v->ft; 
  //  v->ft = &val;  
  v->ft->Moveto(&val); // warning the pointer v->ft must be unchange 
  *an=save;
}       

void IarrayP0::execute () 
{
  // 
  Analvar save(*an);

  assert(v->type == Iden::ftableauP0);
  
  CTabP0 & fP0 = *v->fP0; 
  Grid * gn = fP0.gridid->fg;// the new grid
  gn->check();    
  Grid * go = fP0.g; // the old grid 
  Grid * ggggg=gn;
  int simple = 1;
  if (Th) 
    {
      ggggg = Th->fg;
      assert(ggggg);
      simple = ggggg == gn;
    }
  assert(gn); // the new grid must exist exect time
  // cout << " old g " << go << " new g " << gn << endl;
  CTabP0  val(&fP0); // creation of the array of value new size 
  if (!simple)
    for (int i = 0; i < gn->nt; i++) 
      val[i] =0;
  an->gridxyng = ggggg;
  
  float xl[]={ 1./3.,1./3.,1./3.};      
  if (simple)
    for (int i = 0; i < ggggg->nt; i++)
      {
        int oldlocal = an->local;
        const bTriangle & T = ggggg->t[i];
        const bVertex & v0 = *T.v[0];
        const bVertex & v1 = *T.v[1];
        const bVertex & v2 = *T.v[2];
        float x = v0.x*xl[0] + v1.x*xl[1] + v2.x*xl[2];
        float y = v0.y*xl[0] + v1.y*xl[1] + v2.y*xl[2];                        
        an->setAn(0,x, y, T.where, xl,-1,-1,i);
        val[i] = e->eval(); // with  check
        an->local  = oldlocal;
      }
  else
    for (int i = 0; i < ggggg->nt; i++)
      {  
        int oldlocal = an->local;
        const bTriangle & T = ggggg->t[i];
        const bVertex & v0 = *T.v[0];
        const bVertex & v1 = *T.v[1];
        const bVertex & v2 = *T.v[2];
        int kk=10;
        float dkk= 1./kk;
        float  coef = T.area*(kk*(kk+1))/2;
        for (int ii=0;ii<kk;ii++)
          for (int jj=0;jj<=ii;jj++)
            { 
              xl[0] = dkk*(1./3.+ii);
              xl[1] = dkk*(1./3.+kk-jj-1);
              xl[2] = 1-xl[0]-xl[1];;
              assert (xl[0] >0 && xl[1]>0 && xl[2] >0);
              float x = v0.x*xl[0] + v1.x*xl[1] + v2.x*xl[2];
              float y = v0.y*xl[0] + v1.y*xl[1] + v2.y*xl[2];                        
              an->setAn(0,x, y, T.where, xl,-1,-1,i);
              double a[3];
              int kt=FindTriangle(*gn->Th,x,y,a);
              if(kt>0)  val[kt] += e->eval()*coef;
            }
        an->local  = oldlocal;
        for (int i = 0; i < gn->nt; i++) 
          val[i] /= gn->t[i].area;
        
      }
  an->gridxyng =0;
  v->fP0->Moveto(&val); // warning the pointer v->ft must be unchange not the value
  //   v->ft->detroy(); 
  //  v->fP0 = &val;
  //  for (int i=0;i<gn->nt;i++) 
  //   cout << val[i] << ((i%19)==9 ? "\n" :" "); 
  // cout << endl << endl;
  *an=save;
}       

void CTab::Moveto(CTab *v)
{
  destroy();
  if(g) g->DelRef();
  size = v->size;
  cc = v->cc;
  g = v->g;
  v->size=0; 
  v->cc=0;
  v->g=0;
}
void CTabP0::Moveto(CTabP0 *v)
{
  destroy();
  if(g) g->DelRef();
  size = v->size;
  cc = v->cc;
  g = v->g;
  v->size=0; 
  v->cc=0;
  v->g=0;
}

void CTab::resize(Grid *gg)
{
  if (g != gg )
    { //cout << "\t\t CTab::resize old " << g << " new  " <<  gg  << " " << cc << endl;
      if (cc) g->DelRef();
      gridid->fg = gg;
      g = gg;
      g->AddRef();
      Vresize();
    }
}
void CTabP0::resize(Grid *gg)
{
  if (g != gg )
    { //cout << "\t\t CTab::resize old " << g << " new  " <<  gg  << " " << cc << endl;
      if (cc) g->DelRef();
      gridid->fg = gg;
      g = gg;
      g->AddRef();
      Vresize();
    }
}
void Isavemesh::execute()
{
  int i;
  Grid* g = id->fg;
  //cout << "\t\t save  mesh " << id->nom << " g = " << g << " dans " << fname << endl;
  g->check();

  assert(g);
  int famfmt=0;
  int famdba=0;
  int fmsh=0;
  int fdbg = 0;
  char * ff = fname->eval();
  int l = strlen(ff);  
  famfmt = !strcmp(ff+l-7,".am_fmt");
  int fam = !strcmp(ff+l-3,".am");
  int fnopo = !strcmp(ff+l-5,".nopo");
  famdba = !strcmp(ff+l-6,".amdba");
  fmsh = !strcmp(ff+l-4,".msh");
  fdbg = !strcmp(ff+l-4,".dbg");
  Triangles * Th = g->Th;
  // attention PB renu
  // ordre   = renu 
  // Int4 *renu = new Int4[Th.nbv];
  if (famfmt) 
    {
      if(verbosity >1)
        cout << "\t\t -- write am_fmt file " << ff <<endl;
      ofstream fout(ff);
      if(fout) Th->Write_am_fmt(fout);
    }
  else    if (fam) 
    {
      if(verbosity >1)
        cout << "\t\t -- write am file " << ff <<endl;
      ofstream fout(ff);
      if(fout) Th->Write_am(fout);
    }
  else    if (fnopo) 
    {
      if(verbosity >1)
        cout << "\t\t -- write nopo file " << ff <<endl;
      ofstream fout(ff);
      if(fout) Th->Write_nopo(fout);
    }
  else if (famdba) 
    {
      if(verbosity >1)
        cout << "\t\t -- write amdba file " << ff <<endl;
      ofstream fout(ff);
      if(fout) Th->Write_amdba(fout);
    }
  else if (fmsh)
    {
      if(verbosity >1)
        cout << "\t\t -- write msh file " << ff <<endl;
      g->save(ff,0);
    }
  else if (fdbg)
    {
      if(verbosity >1)
        cout << "\t\t -- write dbg file " << ff <<endl;
      g->save(ff,1); // save as a debug file
    }
  else         
    {
      if(verbosity >1)
        cout << "\t\t -- write mesh out file " << ff  << endl;
      ofstream fout(ff);
      if(fout)  fout << *Th << endl;
    }

  //   g->save(fname);
  delete ff;
 
}

Grid * Iscilabmesh::eval()
{
  Triangles * Th = new Triangles(vt,nbvt,tr,nbtr);
 
  if( ! Th) throw ErrorExec("Create scilab mesh");

  Grid * g = new Grid();
  g -> th2t(Th);
  g -> prepgrid(0);

  return g;
}

Grid * Ireadmesh::eval()
{
  Grid* g = new Grid();
  
  // g->readgrid(fname);
  char * ff = fname->eval();
  Triangles * Th = new Triangles(ff);
  delete ff;
  ff=0;
  if( ! Th) throw(ErrorExec("Read mesh"));
  double hmax = Th -> MaximalHmax();
  //  cout << " hmax = " << hmax << " ------- " << endl;
  Metric M(hmax);
  for (int iv=0;iv < Th->nbv;iv++)
    (*Th)[iv].m = M;

#ifdef DRAWING1
  reffecran();
  Th->InitDraw();
  Th->inquire();
#endif 
  
  g->th2t(Th);
  float rr = 0;
  if  (renu) rr = renu->eval();
  if(rr) {
    cout << " Warning renumbering " << endl;
    g->renum();
  }
  g->prepgrid(0);
   
  g->draw(*an->wait->storage);
  // an->activeMesh=g; // set the activegrid
  return g;
}



Grid * Etruncmesh::eval()
{
  Analvar save(*an);
  assert(idgrid && idgrid->type ==Iden::maillage );
  Grid* go = idgrid->fg;
  assert(go);
  int * flag= new int[go->nt];
  int * bb  = new int[go->nt];
  float xl[]={ 1./3.,1./3.,1./3.};      
  for (int i=0;i<go->nt;i++)
    {
      int oldlocal = an->local;
      const bTriangle & T = go->t[i];
      const bVertex & v0 = *T.v[0];
      const bVertex & v1 = *T.v[1];
      const bVertex & v2 = *T.v[2];
      float x = v0.x*xl[0] + v1.x*xl[1] + v2.x*xl[2];
      float y = v0.y*xl[0] + v1.y*xl[1] + v2.y*xl[2];                        
      an->setAn(0,x, y, T.where, xl,-1,-1,i);
      float ee = e->eval();
      flag[i] = (int) Max(-32000.0F,Min(e->eval(),32000.0F));
      if (b) 
        bb[i]  = (int) Max(-32000.0F,Min(b->eval(),32000.0F));
      else 
        bb[i] = 1;
      //   cout << ee  << " " << flag[i] <<  " " << bb[i] << endl;

      an->local  = oldlocal;
    }
  Grid* g = new Grid();
  // for (int i=0;i<go->nt;i++)
  //  cout << flag[i] << (i%10 == 9 ? '\n' : ' ');
  cout << endl;
  Triangles * Th = new Triangles(*go->Th,flag,bb);
  delete [] flag;
  delete [] bb;
  if( ! Th) throw(ErrorExec("trunc triangulation"));
  double hmax = Th->MaximalHmax();
  //  cout << " hmax = " << hmax << " ------- " << endl;
  Metric M(hmax);
  for (int iv=0;iv < Th->nbv;iv++)
    (*Th)[iv].m = M;
  
#ifdef DRAWING1
  reffecran();
  Th->InitDraw();
  Th->inquire();
#endif 
  
  g->th2t(Th);
  g->renum();
  g->prepgrid(0);
   
  g->draw(*an->wait->storage);
  // an->activeMesh=g; // set the activegrid
  *an=save;
  return g; 
}

Grid * Emovemesh::eval()
{
  Analvar save(*an); 
  int i;
  Grid* go = idmoved->fg;
  Grid* gn = new Grid(go);
  assert(go);
  
  Grid& t =*gn;
  an->gridxyng = go;
  float xl[3] = {0.,0.,0.};
  for ( i = 0; i < go->nv; i++)
    {
      int oldlocal = an->local;
      an->setAn(0, go->v[i].x, go->v[i].y, go->v[i].where,xl, i); 
      t.v[i].x = ex->eval();
      t.v[i].y = ey->eval();
      an->local  = oldlocal;
    }
  an->gridxyng =0;
  t.prepgrid(1);
  //  for(i=0;i<t.nt;i++)
  //  for(int iloc=0;iloc<3;iloc++)
  //    t.t[i].e[iloc] = go->t[i].e[iloc];
  if (!toScilab) t.draw(*an->wait->storage);
  Geometry * tGh = new Geometry(go->Th->Gh);
  Triangles* tTh = new Triangles(*go->Th,tGh);// copy the Triangles
  cout << "\t\t MoveMesh Grid * " << gn << " Gh = " <<  tGh << " Th = " << tTh  << endl;
  //tTh->Write("movemesh.Th");
  
  Triangles & Th = *tTh;
  Geometry & Gh = *tGh;
  Geometry & GhO = go->Th->Gh;
  int * renu = go->NumThinGrid;
  
  cout << "\t\t renu = " << renu << " " << (renu ? renu[0] : 0) << endl;
  // move of the geometry 
  for ( i = 0; i <GhO.nbv;i++)
    {
      int oldlocal = an->local;
      an->setAn(0,GhO.vertices[i].r.x, GhO.vertices[i].r.y, GhO.vertices[i].ref,xl);
      Gh.vertices[i].r.x  = ex->eval();
      Gh.vertices[i].r.y  = ey->eval();
      an->local  = oldlocal;
    }
  // change the tangente 

  for (i=0;i<Gh.nbe;i++) 
    {
      R2 AB = Gh.edges[i].v[1]->r - Gh.edges[i].v[0]->r;        
      Real8 lAB = sqrt(AB*AB); // length of current edge AB
        
      for (int jj=0;jj<2;jj++) 
        if( ! Gh.edges[i].v[jj]->Corner() &&
            (    (jj==0 && Gh.edges[i].TgA()) 
                 || (jj==1 && Gh.edges[i].TgB()) ) )
          {       
                
            // recompute the tangent
            R2 tg =  Gh.edges[i].v[1-jj]->r 
              - Gh.edges[i].Adj[jj]->v[1-Gh.edges[i].SensAdj[jj]]->r;
            Real8 ltg =  sqrt(tg*tg);
            tg =  tg *(lAB/ltg);
            if ( tg*AB < 0) 
              tg = -tg;
            Gh.edges[i].tg[jj] = tg;                 
          }          
    } // for (i=0;i<nbe;i++)
#ifdef DRAWING
  Gh.InitDraw();
#endif
#ifdef DRAWING2    
  reffecran();
  Gh.InitDraw();
  Gh.Draw();
#endif      
  //  move of the bamg mesh 
  Th.pmin =  R2(t.v[0].x,t.v[0].y);
  Th.pmax =  Th.pmin;

  for ( i = 0; i < Th.nbv; i++)
    {   // Be carefull we do a renumbering 
      int j = renu ? renu[i] : i;
      Th.vertices[i].r =   R2(t.v[j].x,t.v[j].y);
      Th.pmin.x = Min( Th.pmin.x, Th.vertices[i].r.x);
      Th.pmin.y = Min( Th.pmin.y, Th.vertices[i].r.y);
      Th.pmax.x = Max( Th.pmax.x, Th.vertices[i].r.x);
      Th.pmax.y = Max( Th.pmax.y, Th.vertices[i].r.y);          
    }
  {R2 P10 = (Th.pmax-Th.pmin)*0.1;
    Th.pmin = Th.pmin - P10;Th.pmax = Th.pmax + P10;}

  Gh.pmin =  Gh.vertices[0].r;
  Gh.pmax =  Gh.pmin;

  for ( i = 0; i < Gh.nbv; i++)
    {
      Gh.pmin.x = Min( Gh.pmin.x, Gh.vertices[i].r.x);
      Gh.pmin.y = Min( Gh.pmin.y, Gh.vertices[i].r.y);
      Gh.pmax.x = Max( Gh.pmax.x, Gh.vertices[i].r.x);
      Gh.pmax.y = Max( Gh.pmax.y, Gh.vertices[i].r.y);          
    }
  {R2 P10 = (Gh.pmax-Gh.pmin)*0.1;
    Gh.pmin = Gh.pmin - P10;Gh.pmax = Gh.pmax + P10;}
    
  if (renu) delete [] renu;
  renu=0;
  
  Gh.coefIcoor = (MaxICoor)/(Max( Gh.pmax.x- Gh.pmin.x, Gh.pmax.y- Gh.pmin.y));
  Th.coefIcoor = (MaxICoor)/(Max( Th.pmax.x- Th.pmin.x, Th.pmax.y- Th.pmin.y));

  //  for ( i = 0; i < Gh.nbv; i++)
  //    Gh.vertices[i].i = Gh.toI2(Gh.vertices[i].r);
  
  for ( i = 0; i < Th.nbv; i++)
    Th.vertices[i].i = Th.toI2(Th.vertices[i].r);
  // remove all the adj  and save the flag
  for ( i= 0; i < Th.nbt; i++)
    { 
      Triangle & t = Th(i);
      for ( int j = 0 ; j<3; j++) 
        t.SetAdj2(j,0,t.GetAllflag(j));
    }
   

  //  Th.quadtree=new QuadTree(&Th);
  Th.nbt = Th.nbt - Th.NbOutT; // remove all the  the ouside triangles 
  Th.SetIntCoor("In movemesh"); 
  Th.FillHoleInMesh(); 
  // delete Th.quadtree; // delete old 
  if (!Th.quadtree)
    Th.quadtree=new QuadTree(&Th);
  Th.ReMakeTriangleContainingTheVertex();
  

#ifdef DRAWING2 
  Th.inquire();
#endif  
  gn->Th = &Th;
  gn->Gh = &Gh;
  Gh.NbRef++;
  gn->nbholes = go->nbholes;
  *an=save;  
  return gn;
}



Grid * Edaptmesh::eval()
{       // tloc is old mesh, g is new mesh
  //Grid* go = bthid->fg->Th;
  Grid & go = *bthid->fg; // 
  Triangles* BTh = go.Th; // the back ground mesh 
  Triangles& Th = *BTh;
  
  assert(BTh);
  
  Real8 hmax = 0.3*Th.MaximalHmax(); // final largest edge 
  Real8 hmin = Th.MinimalHmin();        // final smallest edge
  
  Real8 coef =1;                // a priori don't touch
  // gestion des arguments 
  if(arg[0]) hmin = Max(hmin,(double) arg[0]->eval());
  if(arg[1]) hmax = Min(hmax,(double) arg[1]->eval());
  Real8 err = arg[2] ? (Real8)arg[2]->eval() :0.01 ;    // coef in the metric
  Real8  errg = Min(0.1,arg[3] ? (Real8) arg[3]->eval() : err);
  int nbsx = arg[4] ? Max(100,(int) arg[4]->eval()): 9000;
  int nbsmooth = arg[5] ? (int) arg[5]->eval() : 3;
  int nbjacobi = arg[6] ? Max(0, (int) arg[6]->eval()) :1;              // if increased will be more smooth
  const Real8 raison = arg[7] ? arg[7]->eval() : 1.8;
  const Real8 omega =  arg[8] ? (Real8) arg[8]->eval() : 1; 
  int iso = arg[9]? (int) arg[9]->eval(): 0;
  int AbsError =  arg[10]? (int) arg[10]->eval(): 1;
  assert(11<NbArgAdapMesh);
  Real8 CutOff = arg[11]? (Real8) arg[11]->eval(): 1.0e-6;
  verbosity=  arg[12]? (int) arg[12]->eval(): verbosity;   
  int inq = arg[13] ? (int) arg[13]->eval(): 0;
  int SplitEdgeWith2Boundary = arg[14] ? (int) arg[14]->eval(): 1;
  double maxsubdiv = arg[15] ? Max(Min((double) arg[15]->eval(),10.0),0.1): 10;
  double anisomax =  arg[16] ? Max((double) arg[16]->eval(),1.0) : 1.0e6;
  int rescaling = arg[17] ? (int) arg[17]->eval() : 1;
  if (iso)  anisomax=1;
  if (verbosity) 
    {
      cout << endl  << endl; 
      cout << " \t\t ## adapt : nbsol= " << nbsol << ", nbsx = " << nbsx << ", err = " << err ;
      cout << ", hmin = " << hmin << ", hmax = " << hmax <<endl;
      cout << " \t\t    ratio  = " << raison << ", nbsmooth = " << nbsmooth ;
      cout << ", omega = " << omega <<  ", coef = " << coef << ", iso = " << iso << endl;
      cout << "     AbsError =" << AbsError << ", CutOff = " << CutOff << ", nbjacobi = " << nbjacobi <<endl;
      cout << " \t\t    maxsubdiv = " << maxsubdiv << " splitpbedge = " << SplitEdgeWith2Boundary ;
      cout << "\t\t, anisomax = " << anisomax << ", rescaling = " << rescaling << endl << endl ; 
    }
  Int4 i,iv; 
  double * lessol = new double [Th.nbv*nbsol];
  double *ss = lessol;
  // be careful because renum --
  // the triangle was no renum 
  an->gridxyng =  &go;
  for ( iv=0;iv<Th.nbv;iv++) 
    Th[iv].color=1; // color 
  float xl[3]= {0.,0.,0.};
  for (Int4  it = 0; it < go.nt; it++)
    for (Int4  jt = 0; jt < 3; jt++)
      { 
        Vertex & v= Th(it)[jt];
        if (&v && v.color)
          {
            v.color =0; // uncolor
            xl[jt]=1;
            int oldlocal = an->local;
            an->setAn(1,v.r.x, v.r.y, v.ref,xl,Th.Number(v),jt,it);
            xl[jt]=0;
            ss = lessol + nbsol* Th.Number(v);
            for (int j =0; j < nbsol; j++)
              *ss++= sol[j]->eval() ;
            an->local  = oldlocal;
          }
      }
  an->gridxyng=0;
  // computation of the metric --- 
  // better thing -> create keyword in the language 
  //    a faire F Hecht .
  Metric Mhmax(hmax);
  for ( iv=0;iv<Th.nbv;iv++) 
    Th[iv].m = Mhmax;
 
  Th.IntersectConsMetric(lessol,nbsol,0,hmin,hmax,sqrt(err)*coef,anisomax,AbsError?0.0:CutOff,nbjacobi,rescaling,0);
  delete [] lessol;
  Th.IntersectGeomMetric(errg,iso);
  Th.SmoothMetric(raison);
  Th.MaxSubDivision(maxsubdiv);
  Th.BoundAnisotropy(anisomax);
  // end of metric's computation 
  
  Triangles* nTh = new Triangles(nbsx,Th); // Adaption is here
  
 
  if(SplitEdgeWith2Boundary)
    nTh->SplitInternalEdgeWithBorderVertices();
  
  if(verbosity>2) 
    nTh->ShowHistogram();
  if (nbsmooth)
    nTh->SmoothingVertex(nbsmooth,omega);
  if(verbosity>2 && nbsmooth) 
    nTh->ShowHistogram();
  Metric M(hmax);
  for (iv=0;iv < BTh->nbv;iv++)
    Th[iv].m = M;
   
#ifdef DRAWING
  if (inq) {
    reffecran();
    nTh->InitDraw();
    nTh->Draw();
    nTh->inquire();}
#else
  inq=0;
#endif   
  Grid* g = new Grid();
  
  tloc = g; 
  g->Th = new Triangles (*nTh); // copy to remove link to background mesh and compress--
  delete nTh; // delete old  
  g->Th->MakeQuadTree();
  g->Th->ReMakeTriangleContainingTheVertex();
  //  g->Th = nTh;
  // g->Gh = &g->Th->Gh;
  
  g->th2t(g->Th);
  g->Gh->NbRef++;
  g->renum();
  g->prepgrid(1);
  if (!inq && !toScilab)
    g->draw(*an->wait->storage);
  if (verbosity) 
    cout << "\t\t ### End of  adapt Mesh ### "<<  g->Th->identity  <<"  NbRef Geom " << g->Gh->NbRef << endl ;
  //  an->activeMesh = g;
  return g;
}

void initdraw ( Grid *g,int &init,float  wait)
{
  static   float xmn, xmx, ymn, ymx,x,y;
  if (ginit && wait >1.5) {cadreortho(gxctr,gyctr,gray);return;}
  if(!init) xmn=g->v[0].x, xmx=xmn, ymn=g->v[0].y, ymx=ymn;
  init++;
  for(int  i=1; i<g->nv; i++ )
    {
      x=g->v[i].x; y=g->v[i].y;
      xmx=x>xmx?x:xmx;
      xmn=xmn>x?x:xmn;
      ymx=y>ymx?y:ymx;
      ymn=ymn>y?y:ymn;
    } 
  gxctr=(xmx+xmn)/2;
  gyctr=(ymx+ymn)/2;
  gray=((xmx-gxctr)>(ymx-gyctr)?xmx-gxctr:ymx-gyctr)*1.1;
  ginit = 1;
  cadreortho(gxctr,gyctr,gray);   
}
void initdraw ( frontiere *f,float  wait)
{
  static   float xmn, xmx, ymn, ymx,x,y;
  if (ginit && wait >1.5) {cadreortho(gxctr,gyctr,gray);return;}
  int nbs = f->nbp;
  float * cr = f->xy;
  int i;
  if (nbs<=0) return;
  xmx=xmn=cr[0];
  ymx=ymn=cr[1];  
  for( i=0; i<nbs; i++ )
    {
      x=cr[2*i];
      y=cr[2*i+1];
      xmx=x>xmx?x:xmx;
      xmn=xmn>x?x:xmn;
      ymx=y>ymx?y:ymx;
      ymn=ymn>y?y:ymn;
    }
  gxctr=(xmx+xmn)/2;
  gyctr=(ymx+ymn)/2;
  gray=((xmx-gxctr)>(ymx-gyctr)?xmx-gxctr:ymx-gyctr)*1.1;
  ginit = 1;
  cadreortho(gxctr,gyctr,gray);   
}

void initdraw ( OnList list, float  wait)
{
  int init=0;
  if (ginit && wait >1.5) {cadreortho(gxctr,gyctr,gray);return;}
  OnList a=list;
  for ( a=list;a;a=a->n)
    if(a->x)
      initdraw(a->x->fg,init,wait);
 
}

int loopdraw (float & wait,OnList list)
{
  if (toScilab) return 0; // do not loop if we are in Scilab
  if ( !(((int) wait) % 2) ) return 0;
  int r=1;
  
  float x,y;
  char  c = Getxyc(x,y);
  if (c=='=') initdraw(list);
  else if (c=='r') ; // redraw
  else if (c=='+'|| c==(char)253) gray /= 2.,gxctr=x,gyctr=y;
  else if (c=='-'|| c==(char)250) gray *= 2.,gxctr=x,gyctr=y;
  else if (c==3) erreur("Graphical stop. You enter ^c");
  else if (c=='0') wait=0,r=0;
  else if (c=='1') wait=1,r=0;
  else if (c=='2') wait=2,r=0;
  else if (c=='3') wait=3,r=0;
  else r = 0;
  cadreortho(gxctr,gyctr,gray);
  
  return r;
}

int loopdraw (float & wait,Grid *g)
{
  if ( !(((int) wait) % 2) ) return 0;
  int r=1;
  float x,y;
  int init=0;
  char  c = Getxyc(x,y);
  if (c=='=') initdraw(g,init);
  else if (c=='r') ; // redraw
  else if (c=='+'|| c==(char) 253) gray /= 2.,gxctr=x,gyctr=y;
  else if (c=='-'|| c==(char) 250) gray *= 2.,gxctr=x,gyctr=y;
  else if (c==3) erreur("Graphical stop. You enter ^c");
  else if (c=='0') wait=0,r=0;
  else if (c=='1') wait=1,r=0;
  else if (c=='2') wait=2,r=0;
  else if (c=='3') wait=3,r=0;
  else r = 0;
  cadreortho(gxctr,gyctr,gray);   
  return r;
}

int loopdraw (float & wait,frontiere *g)
{
  if ( !(((int) wait) % 2) ) return 0;
  int r=1;
  float x,y;
  int init=0;
  char  c = Getxyc(x,y);
  if (c=='=') initdraw(g,init);
  else if (c=='r') ; // redraw
  else if (c=='+'|| c==(char) 253) gray /= 2.,gxctr=x,gyctr=y;
  else if (c=='-'|| c==(char) 250) gray *= 2.,gxctr=x,gyctr=y;
  else if (c==3) erreur("Graphical stop. You enter ^c");
  else if (c=='0') wait=0,r=0;
  else if (c=='1') wait=1,r=0;
  else if (c=='2') wait=2,r=0;
  else if (c=='3') wait=3,r=0;
  else r = 0;
  cadreortho(gxctr,gyctr,gray);   
  return r;
} 
void Iplot::execute()
{
  //   
  Analvar save(*an);
  an->activeMesh->check();

  float xmn, xmx, ymn, ymx, xctr, yctr, ray,x,y,fmn=1,fmx=0;
  Grid* g = 0;
  // do not write graphics data if we are in Scilab
  if (!toScilab) initdraw(list,*an->wait->storage);
    
  do {
    if (fileps)  {char *ff= fileps->eval(".ps");  openPS(ff);delete ff;}
    if (!toScilab)  // do not write graphics data if we are in Scilab
      {
	reffecran();
	couleur(1);
	assert(ginit); 
      }

    // COMPUTATION OF THE fmn and fmx  
    int fm =0; // fmx,fmx un set
    int kkk=0;
    int kkkk=0;
    for (OnList a=list;a;a=a->n)
      {
        float xl[3]= {0.,0.,0.};
        Grid* oldActiveMesh = an->activeMesh;
        if(a->x) g = a->x->fg;
        else 
          { 
            an->activeMesh = g;
            an->activeMesh->check();
            an->gridxyng=g;
            int* tmp = new int[g->nv];
            for (int i=0;i<g->nv;i++) 
              tmp[i]=0;
            for (int it=0;it<g->nt;it++)
              {
                kkkk++;
                bTriangle & t = g->t[it];
                if (InRecScreen(Min3(t.v[0]->x,t.v[1]->x,t.v[2]->x),
                                Min3(t.v[0]->y,t.v[1]->y,t.v[2]->y),
                                Max3(t.v[0]->x,t.v[1]->x,t.v[2]->x),
                                Max3(t.v[0]->y,t.v[1]->y,t.v[2]->y))) 
                  for (int j=0;j<3;j++)
                    {
                      bVertex & v = *t.v[j];
                      int iv = g->no(&v);
                      if (tmp[iv]==0) // optimisation 
                        {
                          tmp[iv]=1;
                          int oldlocal = an->local;
                          an->setAn(0,v.x, v.y, v.where,xl,iv);
                          float f = a->e->eval();
                          an->local  = oldlocal;
                          if (fm)
                            fmn=Min(f,fmn),fmx=Max(f,fmx);            
                          else 
                            fmn=fmx=f,fm=1;
                          //  cout << iv << " " << fm << " " << fmn << " " << fmx << endl;
                        }
                   
                   
                    } 
                else kkk++;        
                         
              }
            delete [] tmp; 
          }
      }
    // end COMPUTATION OF THE fmn and fmx  

    if (fm && verbosity>4) 
      cout << "    --  local extrema  " << fmn << " " << fmx << " nb triangle outside " << kkk << " over " << kkkk <<  endl; 
    if (fm)     
      for (OnList a=list;a;a=a->n)
        {
          Grid* oldActiveMesh = an->activeMesh;
          if (a->x) g = a->x->fg;
          else 
            { 
              an->activeMesh = g;
              Grid& t = *g;
              float* temp = new float[t.nv];
              an->gridxyng=&t;
              float xl[3]= {0.,0.,0.};
	      if (toScilab) dts = new DataToScilab(g); // save data to Scilab
              for (int i = 0; i < t.nv; i++)
                {
                  int oldlocal = an->local;
                  an->setAn(0,t.v[i].x, t.v[i].y, t.v[i].where,xl,i);
                  temp[i] = a->e->eval();
		  if (toScilab) dts -> an_eval[i] = temp[i]; // save node values
                  an->local  = oldlocal;
                }
              if (!toScilab) equpot(t,temp,20,0,fmn,fmx); // no plot if we are in Scilab
              delete [] temp;        
              an->activeMesh = oldActiveMesh;
              an->gridxyng=0;
            }
        }
    else
      cerr << " nothing to plot " << endl;    
  } while (loopdraw(*(an->wait->storage),list));

  if (fileps) closePS();
  *an=save;
}

void IplotP0::execute()
{
  //  
  Analvar save(*an); 
  const float xl[3]= {1./3.,1./3.,1./3.};
  
  an->activeMesh->check();
  float xmn, xmx, ymn, ymx, xctr, yctr, ray,x,y,fmn=1,fmx=0;
  Grid* g = 0;//id->fg;
  initdraw(list,*an->wait->storage);
    
  do {
    if (fileps )  {char *ff= fileps->eval(".ps");  openPS(ff);delete ff;}
    reffecran();
    couleur(1);
    
    assert( ginit); 
    // {int i=1; cout << " fordebugging " << i << endl;}
    // COMPUTATION OF THE fmn and fmx  
    int fm =0; // fmx,fmx un set
    int kkk=0;
    int kkkk=0;
    for (OnList a=list;a;a=a->n)
      {
        Grid* oldActiveMesh = an->activeMesh;
        if(a->x) g = a->x->fg;
        else 
          { 
            an->activeMesh = g;
            Expr *e= a->e;

            an->activeMesh->check();
            an->gridxyng=g;
            for (int it=0;it<g->nt;it++)
              {
                kkkk++;
                bTriangle & t = g->t[it];
                if (InRecScreen(Min3(t.v[0]->x,t.v[1]->x,t.v[2]->x),
                                Min3(t.v[0]->y,t.v[1]->y,t.v[2]->y),
                                Max3(t.v[0]->x,t.v[1]->x,t.v[2]->x),
                                Max3(t.v[0]->y,t.v[1]->y,t.v[2]->y))) 
                  {
                    int oldlocal = an->local;
                    const bTriangle & T = g->t[it];
                    const bVertex & v0 = *T.v[0];
                    const bVertex & v1 = *T.v[1];
                    const bVertex & v2 = *T.v[2];
                    float x = (v0.x + v1.x + v2.x)/3.;
                    float y = (v0.y + v1.y + v2.y)/3.;                        
                    an->setAn(0,x, y, T.where, xl,-1,-1,it);     
                    float f = e->eval(); // with  check
                    an->local  = oldlocal;
                    if (fm)
                      fmn=Min(f,fmn),fmx=Max(f,fmx);          
                    else 
                      fmn=fmx=f,fm=1;
                    //  cout << iv << " " << fm << " " << fmn << " " << fmx << endl;
                  }
                else kkk++;                              
              }
          }
      }
    // end COMPUTATION OF THE fmn and fmx  
    // {int i=1; cout << " fordebugging " << i << endl;}
    // if (fm) 
    if(verbosity>3)
      cout << "    --  local extrema  " << fmn << " " << fmx << " nb traingle outside " << kkk << " over " << kkkk <<  endl; 
    if (fm)     
      for (OnList a=list;a;a=a->n)
        {
          Grid* oldActiveMesh = an->activeMesh;
          if(a->x) g = a->x->fg;
          else 
            { 
              an->activeMesh = g;
              Grid& t = *g;
              float* temp = new float[t.nt];
              an->gridxyng=&t;
              for (int i = 0; i < t.nt; i++)
                {
                  int oldlocal = an->local;
                  const bTriangle & T = g->t[i];
                  const bVertex & v0 = *T.v[0];
                  const bVertex & v1 = *T.v[1];
                  const bVertex & v2 = *T.v[2];
                  float x = (v0.x + v1.x + v2.x)/3.;
                  float y = (v0.y + v1.y + v2.y)/3.;                        
                  an->setAn(0,x, y, T.where, xl,-1,-1,i);     
                  temp[i]= a->e->eval(); // with  check
                  an->local  = oldlocal;
                }
              equpotP0(t,temp,20,0,fmn,fmx);
              delete [] temp;        
              an->activeMesh = oldActiveMesh;
              an->gridxyng=0;
            }
        }
    else
      cerr << " nothing to plot " << endl;    
  } while (loopdraw(*(an->wait->storage),list));
  if (fileps) closePS();
  *an=save;
}

void Iplot3d::execute()
{
  Analvar save(*an);
  Grid* g = (Grid*)(id->fn);
  Grid* oldActiveMesh = an->activeMesh;
  if (fileps )   {char *ff= fileps->eval(".ps");  openPS(ff);delete ff;}
  an->activeMesh = g;
  an->activeMesh->check();
    
  Grid& t = *g;
  float* temp = new float[t.nv];
  an->gridxyng=g;
  float xl[3] = {0.,0.,0.};
  for (int i = 0; i < t.nv; i++)
    {
      int oldlocal = an->local;
      an->setAn(0,t.v[i].x, t.v[i].y, t.v[i].where, xl,i);
      temp[i] = e->eval();
      an->local  = oldlocal;
    }
  an->gridxyng=0;
  graph3d(t,temp, (int)(*(an->wait->storage)));
  delete [] temp;
  an->activeMesh = oldActiveMesh;
  if (fileps) closePS();
  *an=save;
}
/*
  void Ihelmholtz::execute()
  { // solve a PDE. Assume files out.msh, f.dta,g.dta,h.dta exist for data
  int i,imax,niter; 
  Grid* g = (Grid*)(id1->fn);
  Grid* oldActiveMesh = an->activeMesh;
  an->activeMesh = g;
  Grid& t = *g;
  Laplace p(g);
  an->gridxyng=g;
  for ( i = 0; i < t.nv; i++)
  {    
  an->ivertex=i;
  *(an->x->storage) = t.v[i].x;
  *(an->y->storage) = t.v[i].y;
  *(an->ng->storage) =t.v[i].where;
  p.dif[i] = l1->eval();
  p.sol[i] = l2->eval();
  p.neu[i] = l3->eval();
  p.rhs[i] = l4->eval();
  p.rob[i] = l5->eval();
  p.vis[i] = l6->eval();
  niter = (int)l7->eval();
  }
  an->gridxyng=0;
    
  if(niter>0) p.solvegradconj(niter, 1e-4);
  else
  {
  Profilmatrix<float,float>* bb;
  if(factorize==-1)
  {
  if(id2)
  { if(id2->fn)
  bb = (Profilmatrix<float,float>*)id2->fn;
  else cout << " Matrix does not exist " << endl;
  } else cout << " No name to get the matrix from " << endl; 
  } else bb = new Profilmatrix<float,float>(t.nv,t.low,t.jlow,t.jhigh);
  p.solveprofil(bb, !factorize);
  if(factorize==1)
  {             if(id2) id2->fn = bb; 
  else cout << " No name to store the matrix into " << endl; 
  } else if(factorize==0) bb->destroy();
  }
  f2->resize(g);
  for ( i = 0; i < t.nv; i++) f2->cc[i] = p.sol[i];
  //  equpot(t,f2->cc,20, 2);
  an->activeMesh = oldActiveMesh;
  }
*/

void Iconvec::execute()
{
  Analvar save(*an);
  Grid* oldActiveMesh = an->activeMesh;
  Grid* g = (Grid*)(id->fn);
  Grid& t = *g;
  typedef VectN<float,2> TVN;
  Vector<TVN> f0(t.nv); // function to be convected
  Vector<float> u(t.nv), v(t.nv); // velocity
  TVN foX;
        
  an->activeMesh = g;
  an->local = 0;
  double dt = e3->eval();  
  f1->resize(g);
  if(f2)
    f2->resize(g);
  an->gridxyng=g;
  float xl[3]= {0.,0.,0.};
  int i;
  for ( i = 0; i < t.nv; i++)
    {
      int oldlocal = an->local;
      an->setAn(0,t.v[i].x, t.v[i].y, t.v[i].where, xl,i);
      u[i] = e1->eval(); 
      v[i] = e2->eval();
      f0[i][0] = e4->eval();
      if(f2) f0[i][1] = e5->eval(); else f0[i][1] = 0;
      an->local  = oldlocal;
    }  
  an->gridxyng=0;
  for ( i = 0; i < t.nv; i++)
    {
      foX = convect(t, f0, u,v, dt,i);
      f1->cc[i] = foX[0];
      if(f2) f2->cc[i] = foX[1];
    }
  an->activeMesh = oldActiveMesh;
  *an=save;
}

void Isave::execute ()
{
  Analvar save(*an); 
  int famgnu=0;
  char * ff=fname->eval();
  
     
  int l = strlen(ff);  
  famgnu = !strcmp(ff+l-4,".gnu");
  ofstream file(ff);
  Grid* oldActiveMesh = an->activeMesh;
  Grid* g = (Grid*)(id->fn);
  an->activeMesh = g;
  Grid& t = *g;
  an->gridxyng=g;
  float xl[3]= {0.,0.,0.};
       
  cout << "\t\t save  function " << id->nom << " grid = " << g << " in " << ff << endl;
  delete ff;ff=0;

  if (famgnu)  // store in gnuplot surface of triangles format
    for(int k=0;k<t.nt; k++)
      {
        for(int jloc=0;jloc<4;jloc++)
          {      int iloc = jloc%3;
	    xl[iloc] = 1;
	    int oldlocal = an->local;
	    an->setAn(1, t.t[k].v[iloc]->x, t.t[k].v[iloc]->y, t.t[k].v[iloc]->where,
		      xl, t.no(t.t[k].v[iloc]), iloc,k);
	    file << t.t[k].v[iloc]->x <<'\t'<< t.t[k].v[iloc]->y <<'\t';
	    float r = e->eval();
	    file << r << endl;
	    xl[iloc] = 0;
	    an->local  = oldlocal;
          }
        file << endl << endl;
      }
  else // store in freefem format
    { 
      file<<t.nv <<endl;
      for (int i = 0; i < t.nv; i++)
        {
          int oldlocal = an->local;
          an->setAn(0, t.v[i].x, t.v[i].y, t.v[i].where,xl, i);
          file<<e->eval()<<endl;
          an->local  = oldlocal;
        }  
    }
  an->gridxyng=0;
  an->activeMesh = oldActiveMesh;
  *an=save;
}

void Iread::execute ()
{
  int np;
  char * ff=fname->eval();
  ifstream file(ff);
  delete ff;ff=0;  
  Grid* g = (Grid*)(id->fn);
  file >>np;
  if(np == g->nv)
    {
      f->resize(g);
      for (int i = 0; i < np; i++)
        file >> f->cc[i];
    } else throw(ErrorExec(" Wrong dimension"));
}

void Grid::buildit(frontiere* tfront, float & waitm)
{
  assert(this);
  frontiere& front = *tfront;
  if (!toScilab) reffecran();  
  if(front.nbp <= 0) 
    throw(ErrorExec("Error no point on the frontiere"));
  if (!toScilab)   // do not write graphics data if we are in Scilab
    {
      initdraw(&front,waitm);
      do 
	{
	  reffecran();
	  showbdy(front.nbp,front.xy,front.nbs,front.s, front.hh, front.ng, front.ngf);
	} while (loopdraw(waitm,&front));
    }
  hinterpole=1; // by def interpolation a h 
  Gh = frontiere2Geometry(front);
  Triangles * nTh = new Triangles(50*front.nbp +1000,*Gh);                              //Achtung!!
  nTh->SplitInternalEdgeWithBorderVertices();
  Th = nTh;
  Th = new Triangles(*nTh);
  delete nTh;
  Th->MakeQuadTree();
  double hmax = Th->MaximalHmax();
  Metric M(hmax);
  for (int iv=0;iv < Th->nbv;iv++)
    (*Th)[iv].m = M;

#ifdef DRAWING1
  reffecran();
  Th->InitDraw();
  Th->inquire();
#endif   
  th2t(Th); 
  renum();
  prepgrid(1);
  if (!toScilab) draw(waitm);  // do not write graphics data if we are in Scilab
} 

void frontiere::save(const char * filename) const
{       int i;
  ofstream file(filename);
  assert(!file.fail());
  file << "MeshVersion  0 \r Dimension 2 \r MaximalAngleOfCorner 360"<< endl;
  file<<endl;
        
  file << "Vertices " << nbp << endl;
  for(i = 0; i< nbp; i++)
    file << xy[2*i] <<'\t'<< xy[2*i+1] <<'\t'<< ng[i] << endl;
  file<<endl;

  file << "Edges " << nbs << endl;
  for(i = 0; i< nbs; i++)
    file << s[2*i]+1 <<'\t'<< s[2*i+1]+1 <<'\t'<< ng[i] << endl;
  file<<endl;

  file << "SubDomain " << nbsd << endl;
  for(i = 0; i<nbsd ; i++)
    file << "2\t" << sd[2*i]+1 << "\t1\t"<< i+1 << endl;
  file<<endl;

  file << "Corners " << nbp << endl;
  for(i = 0; i< nbp; i++)
    file << i+1 << endl;
  file<<endl;
        
  file << "End" << endl;
}

void getMetric(Triangles& BTh, // triangukation
               double* f, // solution for adaption
               double err, // level of error desired
               int iso,         // 1 or 0
               int AbsError, //1 or 0
               double CutOff // if relative error less than cutoff then absolute error is used
               );

void getMetric(Triangles& BTh, // triangukation
               double* f, // solution for adaption
               double err, // level of error desired
               int iso,         // 1 or 0
               int AbsError, //1 or 0
               double CutOff // if relative error less than cutoff then absolute error is used
               )
{
  const int nbjacoby = 1; 
  const int nsol = 1;  // if more than one sol then put them at the end of previous one in array sol
  double hmin = BTh.MinimalHmin();
  double hmax = BTh.MaximalHmax();
  Metric Mhmax(hmax);
  for (Int4 iv=0;iv<BTh.nbv;iv++)
    BTh[iv].m = Mhmax;
  BTh.IntersectConsMetric(f,nsol,0,hmin,hmax,sqrt(err),iso, AbsError?0.0:CutOff,nbjacoby,1);
}

void compile(char *fname)
{
  ifstream f(fname);
  if (!f) throw ErrorMemory("Flux d'entree de fichier (analyse.cpp : ligne 4106)");
    
  Analyseur * a = new Analyseur(&f);
  a -> programme();
  delete a;
}

// for link to Scilab

#include "scilink.h"

// ------------------------------------------------
