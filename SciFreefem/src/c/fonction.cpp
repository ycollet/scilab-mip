//file fonction.cpp
#include "fonction.h"
#include <iostream>

using namespace std;

class CMonome : public CVirt {
  float c; int n;
public:
  CMonome (float cc = 0, int k = 0) : c (cc), n (k) {}
  float operator () (float x)
  {float z = c; int k = n;
    while (k-- > 0) z *= x; while (++k < 0) z /= x; return z;}
  CVirt *de () {return n ? new CMonome (c * n, n - 1) : zero;}
};

Fonction Fonction::monome (float x, int k) {return new CMonome (x, k);}

class CFunc : public CVirt {
  float (*f) (float);
public:
  CFunc (float (*ff) (float), CVirt *dd = 0) : f(ff) {md = dd;}
  float operator() (float x) {return (*f)(x);}
};

class CComp : public CVirt {
  CVirt *f; CVirt *g;
public:
  CComp (CVirt *ff, CVirt *gg) : f(ff), g(gg) {}
  float operator() (float x) {return (*f)((*g)(x));}
  CVirt *de () {return Fonction (g->d ()) * Fonction (f->d())(g);} 
};

class CComb : public CVirt {
  CVirt *a1; CVirt *a2; CVirt2(*op);
public:
  CComb (CVirt *aa1, CVirt *aa2, CVirt2 *oop) : a1(aa1), a2(aa2), op(oop) {}
  float operator() (float x) { float a=(*a1)(x),b= (*a2)(x);return (*op) (a,b);}
  CVirt *de () {return Fonction (a1->d ()) * Fonction2 (op->d1 ()) (a1, a2)
      + Fonction (a2->d ()) * Fonction2 (op->d2 ()) (a1, a2);}
};

class CMonome2 : public CVirt2 {
  float c; int n1, n2;
public:
  CMonome2 (float cc = 0, int k1 = 0, int k2 = 0) : c (cc), n1 (k1), n2 (k2) {}
  float operator () (float x, float y)
  {float z = c; int n = n1; while (n-- > 0) z *= x; while (++n < 0) z /= x;
    n = n2; while (n-- > 0) z *= y; while (++n < 0) z /= y; return z;}
  CVirt2 *de1 () {return n1 ? new CMonome2 (c * n1, n1 - 1, n2) : zero2;}
  CVirt2 *de2 () {return n2 ? new CMonome2 (c * n2, n1, n2 - 1) : zero2;}
};

Fonction2 Fonction2::monome (float x, int k, int l) {return new CMonome2 (x, k, l);}

class CFunc2 : public CVirt2 {
  float (*f) (float, float);
public:
  CFunc2 (float (*ff) (float, float), CVirt2 *dd1 = 0, CVirt2 *dd2 = 0) : f(ff)
  {md1 = dd1; md2 = dd2; }
  float operator() (float x, float y) { return  (*f)(x, y);}
};

class CComp2 : public CVirt2 {
  CVirt *f; CVirt2 *g;
public:
  CComp2 (CVirt *ff, CVirt2 *gg) : f(ff), g(gg) {}
  float operator() (float x, float y) {return (*f)((*g)(x, y));}
  CVirt2 *de1 () {return Fonction2 (g->d1 ()) * Fonction (f->d ()) (g);}
  CVirt2 *de2 () {return Fonction2 (g->d2 ()) * Fonction (f->d ()) (g);}
};
// ADD F Hecht pour les fonctions avec des test  operator ? :
class Ctest3: public CVirt2 {
  CVirt2 *ft, *f, *g;
public:
  Ctest3 (CVirt2 *ftt,CVirt2 *aa1, CVirt2 *aa2) : ft(ftt), f(aa1), g(aa2) {}
  float operator() (float x, float y) {return (*ft)(x,y) ? (*f)(x, y): (*g)(x, y);}
  CVirt2 *de1 () {return Fonction2(ft) ? Fonction2 (f->d1 ()): Fonction2 (g->d1 ());}
  CVirt2 *de2 () {return Fonction2(ft) ? Fonction2 (f->d2 ()): Fonction2 (g->d2 ());}
};

class CComb2 : public CVirt2 {
  CVirt2 *a1; CVirt2 *a2; CVirt2(*op);
public:
  CComb2 (CVirt2 *aa1, CVirt2 *aa2, CVirt2 *oop) : a1(aa1), a2(aa2), op(oop) {}
  float operator() (float x, float y) {return (*op) ((*a1)(x, y), (*a2)(x, y));}
  CVirt2 *de1 ()
  {return a1->d1 () * Fonction2 (op->d1 ()) (a1, a2) + a2->d1 () * Fonction2 (op->d2 ()) (a1, a2);}
  CVirt2 *de2 ()
  {return a1->d2 () * Fonction2 (op->d1 ()) (a1, a2) + a2->d2 () * Fonction2 (op->d2 ()) (a1, a2);}
};

CVirt *CVirt::zero = new CMonome;
CVirt2 *CVirt2::zero2 = new CMonome2;

Fonction::Fonction (float (*f) (float)) : p(new CFunc(f)) {}
Fonction::Fonction (float x) : p(new CMonome(x)) {}
Fonction Fonction::operator() (Fonction f) {return new CComp (p, f.p);}
Fonction2 Fonction::operator() (Fonction2 f) {return new CComp2 (p, f.p);}

Fonction2::Fonction2 (float (*f) (float, float)) : p(new CFunc2(f)) {}
Fonction2::Fonction2 (float x) : p(new CMonome2(x)) {}
Fonction Fonction2::operator() (Fonction f, Fonction g) {return new CComb (f.p, g.p, p);}
Fonction2 Fonction2::operator() (Fonction2 f, Fonction2 g) {return new CComb2 (f.p, g.p, p);}

float add (float x, float y) {return x + y;}
float sub (float x, float y) {return x - y;}

Fonction2 Add = new CFunc2 (add, Fonction2 (1), Fonction2 (1));
Fonction2 Sub = new CFunc2 (sub, Fonction2 (1), Fonction2 (-1));
Fonction2 Mul = new CMonome2 (1, 1, 1);
Fonction2 Div = new CMonome2 (1, 1, -1);
Fonction Chs = new CMonome (-1, 1);
Fonction Identity = new CMonome (-1, 1);

Fonction2 Abscisse = new CMonome2 (1, 1, 0);// x
Fonction2 Ordonnee = new CMonome2 (1, 0, 1);// y

