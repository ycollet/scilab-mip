#ifndef MESHIO_H
#define MESHIO_H

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

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <ctype.h>

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

extern  void (*MeshIstreamErrorHandler)(ios & );
extern  int verbosity ;

void WriteStr(ostream & out,char * str);

double * ReadbbFile(const char * file,long & nbsol,long & lsol,const int dim=2,const int typesol=2);
double * ReadBBFile(const char * file,long & nbsol,long & lsol,int *& typesols,const int dim=2,const int typesol=2);
// solution at vertex (P1)

union Char4orLong {  char c[4];    long l;} ;

class  MeshIstream {
 public:
  istream  & in ;
  //  ifstream  fin;
  int LineNumber,LineError,opened;
  const char * CurrentFile;

  istream & cm ()//  mange les blancs et les commentaire 
    { 
      char c;
      int cmm=0;
      while ( in.get(c) &&  
	      ( isspace(c) ? (( ( c=='\n'|| c==char(12) || c==char(15)) && (LineNumber++,cmm=0)),1) : (cmm || (c=='#' && (cmm=1) )) ) );
      if (in.good()) in.putback(c);
      return in;
    }

  // void rewind(){ fin.clear();fin.seekg(0);}
    
  void eol()// go to end of line
  { 
    char c;
    while ( in.get(c) &&  ( c!='\n') );
  }
  void ShowIoErr(int );
  MeshIstream  & err () 
    { 
      if ( ! in.good() ) ShowIoErr(in.rdstate());
      return *this;
    }
  //  MeshIstream(istream & i): in(i),CurrentFile(0),LineNumber(1),LineError(0) {}

 MeshIstream(const char * file_name)
   : in(*new ifstream(file_name)), LineNumber(1),LineError(0),CurrentFile(file_name)
  {    if (!in) {cerr << " Error Opening file " << file_name,CurrentFile=0;ShowIoErr(1);}
    if(verbosity>4) cout << "    Openfile : " << file_name << endl;err();  }

  /*  //  void close()
      {
      if (CurrentFile) {
      if(verbosity>5) cout << "    Closefile: " <<  CurrentFile << endl;
      CurrentFile=0;in.close();}
      } */
  int eof(){return in.eof();}
  ~MeshIstream(){delete &in;}
  int IsString(char* s);
  char * ReadStr();
  MeshIstream&   operator>>(short& i)   { cm() >> i ;return err();}
  MeshIstream&   operator>>(long& i)   { cm() >> i ;return err();}
  MeshIstream&   operator>>(int& i)   { cm() >> i ;return err();}
  MeshIstream&   operator>>(float& i)   { cm() >> i ;return err();}
  MeshIstream&   operator>>(double& i)   { cm() >> i ;return err();}
  MeshIstream&   operator>>(char * & i ) { i=ReadStr();return err();}

};
// Fortran unformatted file  interface ----------

class IFortranUnFormattedFile {
  //  template<class T> friend IFortranUnFormattedFile & operator>>(IFortranUnFormattedFile &f,T & l);
  long i,l,n,j,nb_rec;
  istream * f;
  int to_close;
  const char * file_name;
 public:
 IFortranUnFormattedFile(char *name)
   : i(0),l(0),n((long)-sizeof(long)),
    nb_rec(0),f(new ifstream(name)), to_close(1),file_name(name)
  { if(!*f) Error(0);}

 IFortranUnFormattedFile(MeshIstream & ff)
   : i(0),l(0),n((long)-sizeof(long)),nb_rec(0),
    f(&ff.in), to_close(0),file_name(ff.CurrentFile)
    { if(! *f)  Error(0);}
  
  ~IFortranUnFormattedFile();
  long Record();
  long where(){return j-i;}
  void  read4(char *c,int );// for the fortran 77 char4
  void  read(char * p,const size_t lg);
  void Error(int);
};

class OFortranUnFormattedFile {
  //  template<class T> friend OFortranUnFormattedFile & operator<<(OFortranUnFormattedFile &f,const T & l);
  long i,l,n,j,nb_rec;
  ostream * f;
  const static char * unkown;
  int to_close;
  const char * file_name;
 public:
  
  //YC: : i(0),l(0),n((long) -sizeof(long)),nb_rec(0),f(new ofstream(name,(std::_Ios_Openmode)mm)),to_close(1),file_name(name)
 OFortranUnFormattedFile(char *name,int  mm=ios::trunc)
   : i(0),l(0),n((long) -sizeof(long)),nb_rec(0),f(new ofstream(name,(std::ios_base::openmode)mm)),to_close(1),file_name(name)
  { if(!*f) Error(0);}
 OFortranUnFormattedFile(ostream &ff)
   : i(0),l(0),n((long) -sizeof(long)),nb_rec(0),f(&ff), to_close(0),file_name(unkown)
  { if(!*f) Error(0);}
  
  ~OFortranUnFormattedFile();

  long Record(long ll=0);
  long where(){return j-i;} 
  void write4(const char *c,int );// for the fortran 77 char4
  void write(const char * p,const size_t lg);
  void Error(int );
};

/// ---------- inline -------------------------

inline void  IFortranUnFormattedFile::read(char * p,const size_t lg){  
  f->read(p,lg);
  j+=lg;
  if (j>n) Error(1);
  else if (!f->good()) Error(2) ;
}

inline void  OFortranUnFormattedFile::write(const char * p,const size_t lg){  
  f->write(p,lg);
  j+=lg;
  if (l && j>n)  Error(1);
  else if (!f->good()) Error(2);
}

template<class T> inline
IFortranUnFormattedFile & operator>>(IFortranUnFormattedFile &f,T & l)
{  
  f.read((char *) &l,sizeof(l));return f;
}
/*  bug sur sun  
    template inline 
    OFortranUnFormattedFile & operator<<(OFortranUnFormattedFile &f,const T & l)
    { 
    f.write((char *) &l,sizeof(l));return f;
    }
    on ex les template  */ 

inline 
OFortranUnFormattedFile & operator<<(OFortranUnFormattedFile &f,const int  & l)
{ 
  f.write((char *) &l,sizeof(l));return f;
}
inline 
OFortranUnFormattedFile & operator<<(OFortranUnFormattedFile &f,const long  & l)
{ 
  f.write((char *) &l,sizeof(l));return f;
}
inline 
OFortranUnFormattedFile & operator<<(OFortranUnFormattedFile &f,const double  & l)
{ 
  f.write((char *) &l,sizeof(l));return f;
}
inline 
OFortranUnFormattedFile & operator<<(OFortranUnFormattedFile &f,const float & l)
{ 
  f.write((char *) &l,sizeof(l));return f;
}

inline void OFortranUnFormattedFile::write4(const char *c,int ll)
{
  int i,j;
  Char4orLong ch4;
  for ( i=0;i<ll;i++)
    {
      ch4.l=0;
      for (j=0;j<4;j++)
	ch4.c[j] = (*c) ? (*c++):(' ');
      *this << ch4.l;
    }
}
inline void IFortranUnFormattedFile::read4(char *c,int ll)
{
  int i,j;
  Char4orLong ch4;

  for ( i=0;i<ll;i++)
    {
      *this >> ch4.l;
      for (j=0;j<4;j++)
	*c++= ch4.c[j];
    }
  *c=0;// end of string 
}
#endif
