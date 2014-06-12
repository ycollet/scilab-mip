// ********** DO NOT REMOVE THIS BANNER **********
//
// SUMMARY: Bamg: Bidimensional Anisotrope Mesh Generator
// RELEASE: 0 
// USAGE  : You may copy freely these files and use it for    
//          teaching or research. These or part of these may   
//          not be sold or used for a commercial purpose with- 
//          out our consent : fax (33) 1 39 63 55 14       
//
// AUTHOR:   D. Bernardi, F. Hecht,  O. Pironneau ,    Y. Darmaillac                      
// ORG    :  INRIA
// E-MAIL :   Frederic.Hecht@Inria.fr   
//
// ORIG-DATE:     Dec 97

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include <iostream>

using namespace std;

#define F_OK 0
#define reel float

template<class T> inline T Min (const T &a,const T &b){return a < b ? a : b;}
template<class T> inline T Max (const T &a,const T & b){return a > b ? a : b;}

extern "C"
{
  static FILE *psfile = 0;
  static int  width, height;
  static reel echx,echy,rxmin,rxmax,rymin,rymax;

  void initgraphique() {}
  void closegraphique() {}
  void showgraphic() {}
  void rattente(int waitm) {}
  void cadreortho(float centrex, float centrey, float rayon) {}
  void couleur(int c) {}
  int  LaCouleur() {return 0;}
  void pointe(float x, float y) {}
  void plotstring(const char *s) {}
  void rmoveto(float x, float y) {}
  void rlineto(float x, float y) {}
  void penthickness(int ) {}
  void cercle(float centrex, float centrey, float rayon) {}
  void reffecran() {}
  void raffpoly(int n, float *poly) {}
  void fillpoly(int n, float *poly) {}
  char Getxyc(float &x,float &y) {return '0';}
  void SetColorTable(int nb) {}
  void GetScreenSize(int &ix,int &iy) {}
  void coutmode(short i) {}
  float  GetHeigthFont() { return 1.0;}

  int getprog(char* fn,int argc, char **argv)
  {
    strcpy(fn,argv[1]);
    printf(" file : %s\n",fn);
    return argc;
  }
  
  void execute (char * str)
  { 
    system(str);
  }

  int InRecScreen(reel x1, reel y1,reel x2, reel y2)
  {  
    return (Max(x1,x2)>= rxmin) && (Min(x1,x2) <= rxmax) && (Max(y1,y2) >= rymin) && (Min(y1,y2)  <= rymax);
  }

  int InPtScreen( reel x, reel y)
  {
    return (x >= rxmin) && (x <= rxmax) && (y >= rymin) && (y <= rymax);
  }

  void cadre(reel xmin,reel xmax,reel ymin,reel ymax)
  {
    rxmin = xmin;
    rxmax = xmax;
    rymin = ymin;
    rymax = ymax;
    
    echx = width / (xmax - xmin);
    echy = height / (ymax - ymin);
  }
  
  void getcadre(reel &xmin,reel &xmax,reel &ymin,reel &ymax)
  {
    xmin = rxmin;
    xmax = rxmax;
    ymin = rymin;
    ymax = rymax;
  }

  void closePS(void)
  {
    if(psfile) {
      fprintf(psfile,"showpage\n");
      fclose(psfile);
    }
    psfile=0;
  }

  void openPS(const char *filename )
  { 
    char ffff[32];
    int count=0;
    if(psfile) closePS();
    time_t *timer,t_loc;
    float s=0.5;
    char  username[10];
    time(&t_loc);
    FILE * fileid = NULL;
    if( !filename) 
      {
	sprintf(ffff,"txtgraph_%.3d.ps",count++);
	fileid = fopen(ffff,"w");
	while (!fileid && count<1000)
	  {
	    sprintf(ffff,"txtgraph_%.3d.ps",count++);
	    fileid = fopen(ffff,"w");
	  }

	fclose(fileid);
	cerr << " The postscript file is : " << ffff << endl;
      }
    
    const char *fps (filename?filename:ffff);
    
    
    psfile=fopen(fps,"w");
    if(psfile) 
      {
	fprintf(psfile,"%%!\n%%%%Creator: %s\n%%%%Title: FremFem+\n",username);
	fprintf(psfile,"%%%%CreationDate: %s",ctime(&t_loc));
	fprintf(psfile,"%%%%Pages: 1\n");
	fprintf(psfile,"%%%%BoundingBox:       0 0 %d %d\n",int(width*s),int(height*s));
	fprintf(psfile,"%%%%EndComments\n");
	fprintf(psfile," /L {newpath moveto lineto stroke} def\n");
	fprintf(psfile," /C {setrgbcolor} def\n");
	fprintf(psfile," /rec {newpath 4 copy 8 1 roll moveto 3 -1 roll lineto 4 2 roll exch lineto lineto closepath} def\n");
	fprintf(psfile," %f %f  scale \n",s,s);
	fprintf(psfile," 0 %d 0 %d rec clip\n",int(width),int(height));
	fprintf(psfile," /Helvetica findfont 16 scalefont setfont\n");
	fprintf(psfile," /S {moveto show} def\n");
      }
    else 
      cerr << " Err openning postscript file " << fps << endl;
  }
  
  void *safecalloc(size_t nb, size_t  size)
  {
    void* p=NULL;
    p = calloc(nb, size);
    if (p == NULL) printf("Run out of Memory!\n");
    return p;
  }
  
  void safefree(void** f)
  {
    if(*f){ free((char*) *f); *f=NULL;}
  }
  void myexit(int err=0) {}

  void message(char *s)
  { 
    printf("%s\n",s); 
  }

  void erreur(char *s)
  { 
    printf("%s\n",s); 
    exit(0);
  }
}
