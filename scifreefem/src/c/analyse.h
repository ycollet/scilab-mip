#ifndef ANALYSE_H
#define ANALYSE_H

#include <iostream>
#include <string>

using namespace std;

#include <string.h>

#include "fonction.h"

void compile(char *fname);

const float penal = 1e5; // default value for the Dirichlet penalization
const int edpmax = 3;   // max number of unknow functions in edp systems


// same order as static char *SymbolName[] = in analyse.cpp

typedef enum {
  lpar,                 rpar,           lbrace,         rbrace,                 cste,
  iden,                 fdecl,          _plus,           _minus,          star,
  slash,                power,          lt,             le,                     gt, 
  ge,                   eq,             ne,             comma,          semicolon, 
  colon,                et,             ou,             si,                     alors, 
  autrement,            pour,           jusqua,         faire,          error, 
  end,                  becomes,        mesh,           border,                 fdecl2, 
  der,                  plot,           plot3d,         print,          RR2,
  chaine,               buildmesh,      savemesh,       readmesh,               movemesh,
  save,                 lire,           progexit,       helmholtz,              adaptmesh, 
  intt,                 Heavyside,      sur,            assemble,               varsolve, 
  convec,               dx,             dy,             array,          laplace , 
  dxx ,                 dxy ,           dyx ,           dyy,            dnu,  
  pde,                  solve ,         avec,           append,         number,
  subroutine,           tilde ,         arrayP0,        arrayP1,        plotP0        
} Symbol;

class CTab;
class CTabP0;
class Grid;
class Analvar;

class Iden {
 public:
  char *nom;
  enum TypeIden {
    reserve, inconnu, variable, fonction, fonction2, 
    fonctionN, ftableau, maillage,courbe, matrix , ftableauP0 ,
    pint , pfloat
  } type;
  union {float *storage; long *istorage;};
  union {Symbol sym; CVirt *f; CVirt2 *f2; CTab *ft; Grid* fg; void *fn; CTabP0 *fP0;int *pi;float * pf;};
  Iden (const char *s) : type(inconnu), storage(0), fn(0) {
    nom = new char[strlen(s) + 1]; strcpy(nom, s);  }
  void newVar (int size=1) {
    type = variable; 
    if(storage) delete [] storage;
    storage = new float[size];
    fn = 0; }
  
  CTab* iden2CTab(Iden * grididen, Analvar *an);
  CTabP0* iden2CTabP0(Iden * grididen, Analvar *an);
  // operator for cast
  operator Fonction  () const {assert(type==fonction);return f;}
  operator Fonction2  () const {assert(type==fonction2);return f2;}
};

class IdenTable {
  int nb;
  Iden** list;
  Iden *insert(const char *s, int n);
 public:
  void enregistre(const char *s, Symbol sym);
  void enregistre(const char *s, int* );
  void enregistre(const char *s, float * );
  void enregistre(const char *s, const Fonction&);
  void enregistre(const char *s, const Fonction2&);
  void enregistre(const char *s, float (*f) (float,float))
  {enregistre(s,Fonction2(f));}
  void enregistre(const char *s, float (*f) (float))
  {enregistre(s,Fonction(f));}
    
  IdenTable () : nb(0), list(0) {}
  Iden *find(const char *);
  int remove(const char* s);
  Iden & operator[](const char * s){return *find(s);} // modif F Hecht
};

struct Expr {                                         // Expression
  virtual float eval () = 0;
  virtual Expr *d (Iden *) = 0;
};



struct Instr {                                        // Instruction
  virtual void execute () = 0;
  virtual Instr *d (Iden *) {return this;}
};


// ADD FH pour avoir de expression de chaine 

class EChaine{
  friend class Analyseur ;
  EChaine *next;
  char * s;
  Expr * e;
 public:
  char * eval(char *extention=0);
 EChaine(const char * ss) : next(0),s(strcpy(new char[strlen(ss)],ss)),e(0) {}
 EChaine(Expr * ee) : next(0),s(0),e(ee) {}
};

struct MeshExpr {                      // MeshExpr FH
  virtual Grid * eval () = 0;
  EChaine * fileps;
MeshExpr(EChaine * f=0) : fileps(f) {}
};

class frontiere;
//class Grid;
//class Laplace;
//typedef  Laplace* Laplacept;

class EDP
{ public:
  int n;                                // number of equations
  int j;                                // position in id[] of the current pde block variable
  Grid* g;
  CTab** f;                     // arrays for the values of the unknown
  Iden* *id;                    // name of unknown functions
  float* sol; // contains bdy cond on entry
  float* neuin; // contains neumann conditions
  float* neuout; // contains neumann conditions
  float* rhs; // right hand side
  float* dis; // the dissipation coefficient
  float* dif; // the diffusion coef (viscosity)
  float* rob; // the robin coef
  float* pdx; // the dx() coef
  float* pdy; // the dy() coef
  float *asym,*pdxy,*pdyx;// the second order coef (dif = (dxx+dyy)/2, (d2asym = dxx-dyy)/2)
};

class Analvar
{ 
 public: 
  // Global variables used during execution
  Iden *x, *y,*ng;              // current x,y,triangle nb and local vertex nb
  Grid* gridxyng;       // for optimisation, the grid on which x,y,ng are defined
  Iden *wait;                   // wait / nowait for graphics
  Iden *nx, *ny;                // normal to boundaries
  Grid* activeMesh;     // activeMesh is at execution time 
  int local;                    //if true then local assembling is done and antrloc and vloc are defined
  int trloc;                    //if vlocnext != vloc it means the function is to b evaluted at mid-point
  float xl[3];                  // barycentric coor of x,y if local
  int ivertex,iloc;             // number of current vertex and local nb if local. iloc = 10*(in+1)+out if edge
  frontiere* front;             // active boundary front
  int bdyLabel;                 // total number of borders
    
  void setAn(int llocal, float xx, float yy, int nng,const  float* xxl,
             int iivertex=-1, int iiloc=-1, int ttrloc=-1)
  {  
    *x->storage = xx; 
    *y->storage = yy; 
    *ng->storage = nng;
    ivertex = iivertex;
    local = llocal;
    iloc = iiloc;
    trloc = ttrloc;
    if (xxl)
      xl[0] = xxl[0], xl[1] = xxl[1], xl[2] = xxl[2];
    if(ivertex>=0)
      {      
	*(nx->storage) = gridxyng->norml[ivertex].x;  //OP-> gridxyng?
	*(ny->storage) = gridxyng->norml[ivertex].y;      
      }  
    if(!local || iloc > 99) return; // not in local mode or it is a point interior to the triangle
    if( iloc >9 && nng) // it is on an edge
      {  
	int kloc = iloc /10 -1;
	int jloc = iloc - 10*kloc - 10;
	bVertex & v1 = *(gridxyng->t[trloc].v[kloc]);
	bVertex & v2 = *(gridxyng->t[trloc].v[jloc]);
	float t1 = v2.x - v1.x, t2 = v2.y - v1.y;
	float nnn = sqrt(t1*t1+t2*t2);
	if(fabs(nnn)<1e-30)
	  {  
	    cout <<" Boundary points are too close: can't get the normal. "<< endl;
	    exit(0);
	  }
	*(nx->storage) = t2/nnn;
	*(ny->storage) = t1/nnn;
      }
  }
};

class NameArg // Modif F. Hecht pour mesh adapa 
{ public:
  Expr **arg;
  const char **ListOfArg;
  int NbOfArg;
 NameArg(int n,const char **  a,Expr **  e): NbOfArg(n),ListOfArg(a),arg(e){};
};

class Analyseur 
{
  Analvar an;
  istream *source;
  Iden *curIden, *curFunc, *curMesh;    // curMesh is at compile time
  frontiere* front;
  EDP* edp;                                                             // current pde block
  float curVal;
  char curChaine[256];
  Symbol curSym;
  int  curIsAlphaNum;                                   // to avoid keywords
  IdenTable table;
  void erreur (const char *s, const char *t = 0);
  Expr* facteur();
  Expr* terme();
  Expr* exprarith();
  Expr* exprcomp();
  Expr* expression();
  EChaine * expchaine(char*errmesg=0);
  MeshExpr* genmesh(); 
  void  readCoef(Expr** exp, int* addmul, int size);
  Instr* instruction();
  void nextSym();
  void lisFonctionN(Iden*);
  void lisFonction2();
  void lisFonction();
  void lisBorder();
  static Instr *rien;
 public:

  Analyseur(istream* s=0);
  Instr * compile(); // for link to scilab 
  const Analvar & GetAnalvarData() const { return an; } // for link to Scilab
  void InitBorder_1(char * , double * , int ,int); // for link to Scilab
  void InitBorder_2(int); // for link to Scilab
  Instr * InitScilabMesh(char *,double *,int,int *,int); // for link to Scilab
  Instr * BuildScilabMesh(char*); // for link to Scilab
  bool NotInitBuffer() const { return (!source); }
  void setBuffer(istream * is) { source = is; }
  void FindArgs(const NameArg & na);
  void match(Symbol s);
  int IsSym(Symbol s);
  int GetExprs(Expr **,int); 
  void programme(); //DB
};

///////////////////////////////////
// Diverses classes d'expressions//
///////////////////////////////////

class EF2 : public Expr
{  // Appel de fonction a deux variables
  Fonction2 f2; const char *nom; Expr *a, *b;
 public:
 EF2(Fonction2 f, const char *n, Expr * aa, Expr *bb) : f2(f), nom(n), a(aa), b(bb) {}
 EF2(float (*f) (float,float) , const char *n, Expr * aa, Expr *bb) : f2(f), nom(n), a(aa), b(bb) {}
 EF2(Iden *f, Expr *aa, Expr *bb) : f2(f->f2), nom(f->nom), a(aa), b(bb) {}
  float eval() {
    //   cout << "xxx(" << a->eval() << "," <<  b->eval() << ")= " <<  f2(a->eval(), b->eval()) << endl;
    float aa(a->eval()),bb(b->eval());
    return f2(aa,bb);}
  Expr *d (Iden * x) {return new EF2 (Add, "der1", 
                                      new EF2 (Mul, "der2", a->d (x), new EF2 (f2.d1 (), "der3", a, b)),
                                      new EF2 (Mul, "der4", b->d (x), new EF2 (f2.d2 (), "der5", a, b)));}
};

class EF : public Expr
{ // Appel de fonction d'une variable
  Fonction f;  const char *nom; Expr * arg;
 public:
 EF(Fonction ff, const char *n, Expr * a) : f(ff), nom(n), arg(a) {}
 EF(Iden *id, Expr *a) : f(id->f), nom(id->nom), arg(a) {}
  float eval() {return f(arg->eval());}
  Expr *d (Iden * x)
  {return  new EF2(Mul, "der6", arg->d (x), new EF(f.d (), "der7", arg));}
};

class EC : public Expr
{ // Constante
  float value;
 public:
  EC (float a) : value(a) {}
  float eval () {return value;}
  Expr *d (Iden *) {return new EC (0);}
};

class EV : public Expr
{ // Variable
  Iden *v;
 public:
 EV(Iden *vv) : v(vv) {}
  float eval() {return *v->storage;}
  Expr * d (Iden *x) {return new EC (x == v);}
};
class EVpint : public Expr
{ // Variable
  Iden *v;
 public:
 EVpint(Iden *vv) : v(vv) {}
  float eval() { return *v->pi;}
  Expr * d (Iden *x) {return new EC (x == v);}
};
class EVpfloat : public Expr
{ // Variable
  Iden *v;
 public:
 EVpfloat(Iden *vv) : v(vv) {}
  float eval() { return *v->pf;}
  Expr * d (Iden *x) {return new EC (x == v);}
};

void interpol(Triangles* Th, float x, float y, int& kt, double* a);

class CTab: public CVirt2 , public P1 {// Class array function
 public: 
  Iden * gridid; 
  Analvar *an;
 CTab(Iden* tt, Analvar *ann):gridid(tt), an(ann),P1(0) {};
 CTab(CTab * ct):gridid(ct->gridid), an(ct->an),P1(ct->gridid->fg) {};
  float operator ()(float x, float y);
  void resize(Grid * gg);
  void Moveto(CTab *v);
         
};

class CTabP0: public CVirt2 , public P0 {// Class array function
 public: 
  Iden * gridid; 
  Analvar *an;
 CTabP0(Iden* tt, Analvar *ann):gridid(tt), an(ann),P0(0) {};
 CTabP0(CTabP0 * ct):gridid(ct->gridid), an(ct->an),P0(ct->gridid->fg) {};
  float operator ()(float x, float y);
  void resize(Grid * gg);
  void Moveto(CTabP0 *v);
};

typedef CTab* CTabpt;
typedef CTabP0* CTabP0pt;
typedef Iden* Idenpt;
typedef Expr* Exprpt;

#endif
