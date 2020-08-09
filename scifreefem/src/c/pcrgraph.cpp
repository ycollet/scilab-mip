//#include <Win32Headers++.mch>
#include <sioux.h>
#include <windows.h>
#include <commdlg.h>
#include <direct.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <setjmp.h>

#include <new>
#include <iostream>
#include <fstream>
//#include <memory.h>

using namespace std; //introduces namespace std
__stdcall int EnumPrintersA(unsigned long, char *, unsigned long, unsigned char *, unsigned long, unsigned long *, unsigned long *)
{ ;}

#include "rgraph.h"

void fillpoly(int n, float *poly){/* a faire */};

template<class T> inline T Min (const T &a,const T &b){return a < b ? a : b;}
template<class T> inline T Max (const T &a,const T & b){return a > b ? a : b;}

void 	out_of_memory ();
void 	NEW_HANDLER (void);
float 	scali(int i);
float 	scalj(int j);
void 	thisexit();
void 	execute(char* what);
char 	Getijc(int & x,int & y);
char 	Getxyc(float &x,float &y);
int 	getcolor();
void	putpixel(int ix,int iy, int couleur);
int 	scalx(float x);
int 	scaly(float y);
void 	compile (char *);
void 	rattente (int);
void 	inittext();
BOOL 	ShowOpenDialogBox(char *fileName);
void 	SetColorTable(int i);
LRESULT CALLBACK (*WndProcfg)( HWND hWnd, UINT messg, 
			       WPARAM wParam, LPARAM lParam );
LRESULT CALLBACK MyWndProc( HWND hWnd, UINT messg,	
			    WPARAM wParam, LPARAM lParam );

char 	szProgName[] = "freefem+"; 
char 	errbuf[255], myiobuf[255];
float 	rayon;
static int width,height;
static FILE *psfile = 0;
static 	float 	aspx, aspy, echx,echy,ech,rxmin,rxmax,rymin,rymax;  
static 	int 	Ax,Ay,Lx,Ly,currx,curry, carre;
int 	NbPlotTotal, jh=20,jv=10;
char 	shortName[256],fullName[256];                     

WNDCLASS rClass;
HDC 	hdc,hdcfg;
HWND    hWnd;  // graph window
HWND 	hWndfg; // text window 

static int  INITGRAPH = 0;
static void *StdCHandle=0;
static int 	cube6[7][3] ={ {65534,0,0},{65534,65534,0},{0,65534,0},
			       {0,65534,65534},{0,0,65534}, {65534,0,65534},{65534,0,0} }; 
static int	nhpen=0;
jmp_buf environ;
static int  myenviron = 0;
static HPEN* 	hpen;

void out_of_memory ()
{
  cout << " error: operator new failed; not enough memory" << endl;
  myexit (1);
}

void erreur(char *s)
{
  char se[256],prj[256];
  strcpy(se,shortName);
  strcat(se," : Syntax Error\n");
  strcat(se,s);
  cout << se << endl;
  myexit(1);
}

void SetColorTable(int n);
void raffpoly(int n, float *poly){}
int  getprog(char* fn,int argc, char** argvptr){  return 0;}
void execute(char* what){}
void penthickness(int pepais){}
void showgraphic(){ShowWindow(hWnd, SW_SHOW );}

void myexit(int i)
{  
  if(INITGRAPH) rattente(i); 
  if(myenviron)
    longjmp(environ,1);
  else   myenviron=0; 
}

void thisexit(){ myexit(0);}

void NEW_HANDLER (void){  set_new_handler (&out_of_memory);}
								

void cadre(float xmin,float xmax,float ymin,float ymax)
{
  RECT rc;
  int a,b,le,Ex,Ey;
  le=20;
  float px,py;
  GetClientRect(hWnd, &rc);
  rxmin = xmin;
  rxmax = xmax;
  rymin = ymin;
  rymax = ymax;
  px=xmax-xmin;
  py=ymax-ymin; 
  if ( px > py)
    {/* Our Object larger than long */
      Ax=le; Lx=((rc.right-rc.left)-2*le);
      Ex=Lx; Ey=(int)(Ex*(py/px));
      b=((rc.bottom-rc.top)-Ey-2*le)/2;
      Ay=b+le; Ly=Ey;
      if (Ey>(rc.bottom-rc.top)-2*le){Ay=le; Ly=((rc.bottom-rc.top)-2*le);
	Ey=Ly; Ex=(int)(Ey*(px/py));
	a=((rc.right-rc.left)-Ex-2*le)/2;
	Ax=a+le; Lx=Ex;
      }
    }
  else
    {/* Our Object longer than large */
      Ay=le; Ly=((rc.bottom-rc.top)-2*le);
      Ey=Ly; Ex=(int)(Ey*(px/py));
      a=((rc.right-rc.left)-Ex-2*le)/2;
      Ax=a+le; Lx=Ex;                     
      if (Ex>(rc.right-rc.left)-2*le){Ax=le; Lx=((rc.right-rc.left)-2*le);
	Ex=Lx; Ey=(int)(Ex*(py/px));
	b=((rc.bottom-rc.top)-Ey-2*le)/2;
	Ay=b+le; Ly=Ey;
      }
     
    }
  echx=1/(xmax-xmin);
  echy=1/(ymax-ymin); 
    
}

int scalx(float x)
{
  int test;
  test=(int)(Ax+Lx*(x-rxmin)*echx);
  return test;
}

int scaly(float y)
{
  int test;
  test=(int)(Ay+Ly*(rymax-y)*echy);
  return test;                                                  
}
void plotstring(const char *s)
{
  if(psfile) fprintf(psfile,"(%s) S\n",s);
}

void rmoveto(float x, float y)
{
  currx = scalx(x);
  curry = scaly(y);
  //MoveTo(hdc,scalx(x), scaly(y));
  if (psfile) 
    fprintf(psfile,"%d %d M\n", currx, height-curry);

}				  

void rlineto(float x, float y)
{
  int newx = scalx(x), newy = scaly(y);
  MoveToEx(hdc,currx,curry,NULL);
  LineTo(hdc,newx,newy);
  if (psfile) 
    fprintf(psfile,"%d %d L\n", newx,height-newy);
 
  currx = newx; curry = newy;
}

void cadreortho(float centrex, float centrey, float ray)
{								  ///
  cadre(centrex-ray, centrex+ray, centrey-ray, centrey+ray);
}

int InRecScreen(float x1, float y1,float x2, float y2)
{  
  float xi = Min(x1,x2),xa=Max(x1,x2);
  float yi = Min(y1,y2),ya=Max(y1,y2);
  return (xa >= rxmin) && (xi <= rxmax) && (ya >= rymin) && (yi <= rymax);
}

int InPtScreen( float x, float y)
{
  return (x >= rxmin) && (x <= rxmax) && (y >= rymin) && (y <= rymax);
}

float scali(int i){  return i/echx  + rxmin;}
float scalj(int j) {  return -j/echy  + rymax;}



LRESULT CALLBACK MyWndProc( HWND hWndd, UINT messg,				/*callback procedure */
			    WPARAM wParam, LPARAM lParam )
{
  HDC hdc; 				/* handle to the device context */
  if (hWndd == hWndfg) 
    return (*WndProcfg)( hWndd, messg, wParam, lParam );
  else 
    switch(messg)
      {
      case WM_PAINT:
	break;
	//case WN_CLOSE:  		  DestroyWindow(hWnd); 		  break;
      case WM_DESTROY: 
	PostQuitMessage(0);
	closegraphique();
	break;	
      default:
	return 0L;//( DefWindowProc( hWndd, messg, wParam, lParam ) );
      }

  return( 0L );
}


void closegraphique(void)
{	
  if(  INITGRAPH) 
    {
      // rattente(1);
      INITGRAPH =0; // before DestroyWindow to avoid loop 
      ReleaseDC(hWnd,hdc); 
      DestroyWindow(hWnd);
      SetWindowLong(hWndfg,GWL_WNDPROC,(long) WndProcfg);	  
    }
}


void initgraphique(void)       
{ 
  if (INITGRAPH) return;
  static char AppName[] = "WinSIOUX" ;
  cout << "Start freefem v0.9.31" <<endl;
  hWndfg = GetForegroundWindow(); //hopefully the text window
  HINSTANCE hInst = (HINSTANCE) GetWindowLong(hWndfg,GWL_HINSTANCE);

  hWnd = CreateWindow( AppName,				/* now create the window */
		       "PC rgraph v.1.0",
		       WS_OVERLAPPEDWINDOW,390,30,405,405,/*
							    CW_USEDEFAULT,
							    CW_USEDEFAULT,
							    CW_USEDEFAULT,
							    CW_USEDEFAULT,*/
		       (HWND)NULL,
		       (HMENU)NULL,
		       hInst,
		       NULL		);
  if (!hWnd) {  cerr << " Cant open a graph window!! " << endl;
    cerr << " HINSTANCE " << hInst << endl;
    cerr << " call back sioux routine " << WndProcfg << endl;
    exit(1); 
  }
  //cout << "initgraphique" << endl;
  /*
    hWndfg = GetForegroundWindow(); //hopefully the text window
    hdcfg=GetDC(hWndfg);
    HINSTANCE hInst = (HINSTANCE) GetWindowLong(hWndfg,GWL_HINSTANCE);
  */
  WndProcfg=(LRESULT CALLBACK (*)( HWND hWnd1, UINT messg,
				   WPARAM wParam, LPARAM lParam )) GetWindowLong(hWndfg,GWL_WNDPROC);
  SetWindowLong(hWndfg,GWL_WNDPROC,(long) &MyWndProc);
  // ShowWindow(hWndfg, SW_SHOW );

  ShowWindow(hWnd, SW_SHOW );
  hdc=GetDC(hWnd);
  hpen=0;
  SetColorTable(2+6);

  RECT rc;
  GetClientRect(hWnd, &rc);
  aspx = (float)(rc.right - rc.left);
  aspy = (float)(rc.bottom - rc.top);
  width = rc.right - rc.left;
  height = rc.bottom - rc.top;
  carre = aspx == aspy;
}

void SetColorTable(int nb)
{
  if (nhpen == nb) return;
  if (nb>2) 
    { 
      if(hpen) delete [] hpen;
      hpen = new HPEN[nb+1];
      nhpen = nb;
      int k=0;
      hpen[k++]= CreatePen(PS_INSIDEFRAME, 1,RGB(255, 255, 255));
      hpen[k++]= CreatePen(PS_INSIDEFRAME, 1,RGB(0, 0, 0));
      nb -= 2;
      for (long i0=0;i0<nb;i0++,k++)
	{  
	  long  i6 = i0*6;
	  long  j0 = i6/nb;// in 0..6
	  long  j1 = j0+1;// in 1..6
	  long  k0 = i0 - (nb*j0)/6L;
	  long  k1 = (nb*j1)/6L-i0;
	  long  kk = (k0+k1);
	  hpen[k]= CreatePen(PS_INSIDEFRAME, 1,RGB(
						   ((cube6[j1][0]*k0+cube6[j0][0]*k1)/kk)%256
						   , ((cube6[j1][1]*k0+cube6[j0][1]*k1)/kk)%256
						   , ((cube6[j1][2]*k0+cube6[j0][2]*k1)/kk)%256
						   ));
	}         
      SelectObject(hdc, hpen[0]);
    }
  else 
    nhpen  =0;
}


void couleur(int c)
{ 
  int c0 = c>nhpen ? nhpen: c;
  SelectObject(hdc, hpen[c0]);
  /*   if (psfile)
       {
       float r=1,g=1,b=1;
       if (colortable) {
       if (c>0 && c < ncolortable)
       {
       r =  (float) colortable[c].red /65535.;
       g =  (float) colortable[c].green /65535.;
       b =  (float) colortable[c].blue /65535.;
       }
       }
       else if (c!=0)
       r=g=b=0;
    
       fprintf(psfile,"%.3f %.3f %.3f C\n",r,g,b);
       }
  */
}



void rattente(int waitm)
{
  //	 char c=' ';
  //    you may prefer to use carriage return to move to the next graph
  //	 getc(stdin);  
  //    ShowWindow(hWnd,  SW_SHOW );
  //    UpdateWindow( hWnd);        
  if(waitm) 
    {
      TextOut(hdc,0,0,"                              ",30);
      TextOut(hdc,0,0,"Click to continue...",19);
      MSG msg;      
      do
	{
	  GetMessage(&msg,hWnd,WM_MOUSEFIRST,WM_MOUSELAST);
	  if (msg.message==WM_RBUTTONDOWN)
	    myexit(2); 
	}
      while (msg.message!=WM_LBUTTONDOWN); 
	
      TextOut(hdc,0,0,"                             ",30);
      TextOut(hdc,0,0,"FreeFem works...",16);
    }
  //    ShowWindow(hWndfg,  SW_SHOW );
  //    UpdateWindow( hWndfg);        
}




/* char Getijc(int & x,int & y){

   char char1=' ';
   if(!INITGRAPH) 
   {
   x = 0;
   y = 0;   
   return char1;
   }
   int  cont=1;
   POINT xy;
   xy.x =0;
   xy.y =0;
   MSG msg;      
   TextOut(hdc,0,0,"                                  ",30);
   TextOut(hdc,0,0,"Click mouse to continue           ",30);
   do
   {   
   GetMessage(&msg,hWnd,0,0);// all message
   GetCursorPos(&xy);
   switch (msg.message)
   {
   case WM_LBUTTONDOWN:char1=char(251), cont=0;break; 
   // with shift 248			   			   
   case WM_RBUTTONDOWN:char1=char(253), cont=0;break;  
   // with shit 250
   // if the 2 buttom, 252, et shith 249;
   case WM_CLOSE: 
   case WM_DESTROY: exit(2);
   case WM_KEYDOWN: char1 = msg.wParam;  cont=0; break;

   }
   }
   while (cont); 
   TextOut(hdc,0,0,"                                  ",33);
   TextOut(hdc,0,0,"freefem works...",16);
   RECT rc;
   GetClientRect(hWnd, &rc);
 
   x = xy.x-rc.left;
   y = xy.x-rc.top;
   cout << " x = " << x << " y = " << (int) char1 << endl;
   return char1;

   }
*/
char Getijc(int & x,int & y)
{
  char char1=' ';
  if(!INITGRAPH)
    {
      x = 0;
      y = 0;
      return char1;
    }
  int  cont=1;
  POINT xy;
  xy.x =0;
  xy.y =0;
  MSG msg;
  TextOut(hdc,0,0,"                                  ",30);
  TextOut(hdc,0,0,"Click mouse to continue           ",30);

  do


    GetMessage(&msg,hWnd,0,0);// all message
  GetCursorPos(&xy);
  switch (msg.message)
    {
    case WM_LBUTTONDOWN:char1=char(251), cont=0;break;
      // with shift 248
    case WM_RBUTTONDOWN:char1=char(253), cont=0;break;
      // with shit 250
      // if the 2 buttom, 252, et shith 249;
    case WM_CLOSE:
    case WM_DESTROY: exit(2);
    case WM_CHAR: char1 = (TCHAR)msg.wParam; cont = 0; break;
      //case WM_KEYDOWN: char1 = (TCHAR)msg.wParam;  cont=0; break;
    default:
      TranslateMessage(&msg);
      DispatchMessage(&msg);
      break;
    }
}
while (cont);

TextOut(hdc,0,0,"                                  ",33);
TextOut(hdc,0,0,"freefem works...",16);
RECT rc;
GetClientRect(hWnd, &rc);

x = xy.x-rc.left;
y = xy.x-rc.top;
cout << " x = " << x << " y = " << y << "char = " << char1 << endl;
return char1;
}
char Getxyc(float &x,float &y)
{ 
  char c=' ';
  int i=0,j=0;
  c = Getijc( i,j);
  x = scali(i);
  y = scalj(j);
  //rattente(1);
  return c;
}


void reffecran(void)
{     
  HBRUSH hbr;
  RECT rc;

  GetClientRect(hWnd, &rc);
  hbr = CreateSolidBrush(RGB(255, 255, 255));
  FillRect(hdc,&rc,hbr);
  DeleteObject(hbr);
}

BOOL ShowOpenDialogBox(char *fileName)
{
  OPENFILENAME ofn; 
  char szDirName[256];   
  char *strFilter="PCgFEM Files (*.edp)\0*.edp\0All Files (*.*)\0*.*\0\0"; 
	
  memset(&ofn, 0, sizeof(OPENFILENAME));
  getcwd(szDirName,sizeof(szDirName));
  ofn.lStructSize = sizeof(OPENFILENAME);
  ofn.hwndOwner = NULL;
  ofn.lpstrFilter = strFilter;
  ofn.lpstrFileTitle = fileName;
  ofn.nMaxFileTitle = 80;
  ofn.lpstrInitialDir=szDirName;
  ofn.lpstrTitle ="Choose you freefem '*.edp' File";
  ofn.Flags=OFN_SHOWHELP|OFN_PATHMUSTEXIST|OFN_FILEMUSTEXIST;
	
  return GetOpenFileName(&ofn);
}    

void coutmode(short r) { ;}// will be done later

int main ( void )
{
  INITGRAPH =0;
  atexit(thisexit);
  /*	hWnd = GetForegroundWindow();
	hdc=GetDC(hWnd);
	hpen=0;
	SetColorTable(2+6); */
  initgraphique();  
  if (ShowOpenDialogBox(shortName)==FALSE) 
    {
      cout << " Annuler " << endl;
      return 1;
    }
  strcpy(fullName,shortName);
  INITGRAPH =1;
  ::UpdateWindow(hWnd);		// move to the top window
  ::SetWindowPos(hWnd, HWND_TOPMOST, 0, 0, 0, 0, SWP_NOMOVE | SWP_NOSIZE);
  if(0==setjmp(environ))
    {
      myenviron = 1;
      compile(fullName);
    }
  myenviron = 0;
  TextOut(hdc,0,0,"                                 ",30);
  cout << "freefem has ended..." << endl;
  TextOut(hdc,0,0,"freefem has ended.              ",30);
  //	if(INITGRAPH)	 closegraphique();
  // RestoreDC(hdc,-1);
  // to test the interface just comment the line above
  // uncomment the line below and make a project with just this file
  //LineTo(hdc,100,100);
  //	rattente(1);
  return 0;
}
void GetSizeScreen(int & ix,int &iy)
{
  RECT rc;
  GetClientRect(hWnd, &rc);
  width = rc.right - rc.left;
  height = rc.bottom - rc.top;

  ix = width ;
  iy = height;
}



void openPS(const char *filename )
{ 
  RECT rc;
  GetClientRect(hWnd, &rc);
  width = rc.right - rc.left;
  height = rc.bottom - rc.top;

  if(psfile) closePS();
  time_t *timer,t_loc;
  float s=0.5;
  char  username[10];
  time(&t_loc);
  printf(" Save Postscript in file '%s'\n",filename?filename:"freefem.ps"),
    psfile=fopen(filename?filename:"freefem.ps","w");
  if(psfile==0) {printf("Erreur %s errno %d\d",filename,errno);exit(1);}
  if(psfile) {
    fprintf(psfile,"%%!\n%%%%Creator: %s\n%%%%Title: FremFem+\n","user");
    fprintf(psfile,"%%%%CreationDate: %s",ctime(&t_loc));
    fprintf(psfile,"%%%%Pages: 1\n");
    fprintf(psfile,"%%%%BoundingBox:       0 0 %d %d\n",int(50+width*s),int(50+height*s));
    fprintf(psfile,"%%%%EndComments\n");
    fprintf(psfile," /L {  lineto stroke} def\n");
    fprintf(psfile," /M {  moveto } def\n");
    fprintf(psfile," /C {setrgbcolor} def\n");
    fprintf(psfile," /rec {newpath 4 copy 8 1 roll moveto 3 -1 roll lineto 4 2 roll exch lineto lineto closepath} def\n");
    fprintf(psfile," 50 50  translate \n");
    fprintf(psfile," %f %f  scale \n",s,s);
    fprintf(psfile," 0 %d 0 %d rec clip newpath\n",int(width),int(height));
    fprintf(psfile," /Helvetica findfont 16 scalefont setfont\n");
    fprintf(psfile," /S { show} def\n");
  }
}
void closePS(void)
{
  if(psfile)   {
    fprintf(psfile,"showpage\n");
    fclose(psfile);
  }
  psfile=0;
}
