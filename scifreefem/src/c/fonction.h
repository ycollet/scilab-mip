//file fonction.h
#ifndef __FONCTION__
#define __FONCTION__

struct CVirt {
  CVirt *md;
  CVirt () : md (0) {}
  virtual float operator () (float) = 0;
  virtual CVirt *de () {return zero;}
  CVirt *d () {if (md == 0) md = de (); return md;}
  static CVirt *zero;
};

class Fonction2;

class Fonction {                              // Fonction d'une variable
  CVirt *p;
 public:
  operator CVirt *() const {return p;}
  Fonction (CVirt *pp) : p(pp) {}
  Fonction (float (*) (float));               // Creation a partir d'une fonction C
  Fonction (float x);                         // Fonction constante
  Fonction (const Fonction& f) : p(f.p) {}    // Copy constructor
  Fonction d () {return p->d ();}             // Derivee
  void setd (Fonction f) {p->md = f;}
  float operator() (float x) {return (*p)(x);}// Valeur en un point
  Fonction operator() (Fonction);             // Composition de fonctions
  friend class Fonction2;
  Fonction2 operator() (Fonction2);
  static Fonction monome (float, int);
};

struct CVirt2 {
  CVirt2 ():  md1 (0), md2 (0) {}
  virtual float operator () (float, float) = 0;
  virtual CVirt2 *de1 () {return zero2;}
  virtual CVirt2 *de2 () {return zero2;}
  CVirt2 *md1, *md2;
  CVirt2 *d1 () {if (md1 == 0) md1 = de1 (); return md1;}
  CVirt2 *d2 () {if (md2 == 0) md2 = de2 (); return md2;}
  static CVirt2 *zero2;
};

class Fonction2 {                                // Fonction de deux variables
  CVirt2 *p;
 public:
  operator CVirt2 *() const {return p;}
  Fonction2 (CVirt2 *pp) : p(pp) {}
  Fonction2 (float (*) (float, float));          // Creation a partir d'une fonction C
  Fonction2 (float x);                           // Fonction constante
  Fonction2 (const Fonction2& f) : p(f.p) {}     // Copy constructor
  Fonction2 d1 () {return p->d1 ();}
  Fonction2 d2 () {return p->d2 ();} 
  void setd (Fonction2 f1, Fonction2 f2) {p->md1 = f1; p->md2 = f2;}
  float operator() (float x, float y) {return (*p)(x, y);}
  friend class Fonction;                
  Fonction operator() (Fonction, Fonction);      // Composition de fonctions
  Fonction2 operator() (Fonction2, Fonction2);
  static Fonction2 monome (float, int, int);
};

extern Fonction Chs,Identity;
extern Fonction2 Add, Sub, Mul, Div, Abscisse, Ordonnee;

inline Fonction operator+ (Fonction f, Fonction g) {return Add(f, g);}
inline Fonction operator- (Fonction f, Fonction g) {return Sub(f, g);}
inline Fonction operator* (Fonction f, Fonction g) {return Mul(f, g);}
inline Fonction operator/ (Fonction f, Fonction g) {return Div(f, g);}
inline Fonction operator- (Fonction f) {return Chs(f);}

inline Fonction2 operator+ (Fonction2 f, Fonction2 g) {return Add(f, g);}
inline Fonction2 operator- (Fonction2 f, Fonction2 g) {return Sub(f, g);}
inline Fonction2 operator* (Fonction2 f, Fonction2 g) {return Mul(f, g);}
inline Fonction2 operator/ (Fonction2 f, Fonction2 g) {return Div(f, g);}
inline Fonction2 operator- (Fonction2 f) {return Chs(f);}

#endif

