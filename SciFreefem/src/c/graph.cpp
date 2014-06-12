//file ugraph2.cpp
/**************DO NOR REMOVE THIS BANNER***************/
/*  FreeFEM : Language for a Finite Element Method    */
/*  -------    Release 1.0:  June 1994.               */
/*  Authors: D. Bernardi, Y. Darmaillac F. Hecht,     */
/*           O. Pironneau                             */
/*  You may copy freely these files and use it for    */
/* teaching or research. These or part of these may   */
/* not be sold or used for a commercial purpose with- */
/* out our consent : fax (33)1 44 27 44 11            */
/* (e-mail)    Olivier.Pironneau@ann.jussieu.fr       */
/******************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <setjmp.h>
#include "vect.h"
#include "rgraph.h"
#include "graph.h"
#include "meshtype.h"
int   *ordre;
float *table;
extern int verbosity;
void Init(Grid&, char *s);
void quicksort(float *tb, int *o, int n);
void initt(void);
void projection(float *f);
void     FillRect(float x0,float y0, float x1, float y1);


void Init(Grid& t, char *s)
{
  int i;
  float xmn=t.v[0].x, xmx=xmn, ymn=t.v[0].y, ymx=ymn;
  float xctr, yctr, ray, x,y;

  for( i=1; i<t.nv; i++ )
    {
      x=t.v[i].x; y=t.v[i].y;
      xmx=x>xmx?x:xmx;
      xmn=xmn>x?x:xmn;
      ymx=y>ymx?y:ymx;
      ymn=ymn>y?y:ymn;
    }
  xctr=(xmx+xmn)/2;
  yctr=(ymx+ymn)/2;
  ray=(xmx-xctr)>(ymx-yctr)?xmx-xctr:ymx-yctr;
  float dx = (xmx-xmn)*0.1;
  float dy = (ymx-ymn)*0.1;
  reffecran();
  if(*s=='o')
    cadreortho(xctr,yctr,ray*1.1);
  else
    cadre(xmn-dx, xmx-dx, ymn+dy, ymx+dy);
  showgraphic();// modif FH 270408
    
}


void quicksort(float *tb, int *o, int n)
{
  int i,j;
  int b;
  float x, y;

  while(n>1)
    {
      x=tb[n/2];
      for(i=0, j=n-1;i<=j;i++, j--)
	{
	  while(tb[i]>x) i++;
	  while(tb[j]<x) j--;
	  if(i>j) break;
	  y=tb[i]; tb[i]=tb[j]; tb[j]=y;
	  b=o[i]; o[i]=o[j]; o[j]=b;
	}
      n-=i;
      if(n>j+1)
	{
	  quicksort(tb, o, j+1);
	  tb+=i;
	  o+=i;
	}
      else
	{
	  quicksort(tb+i, o+i, n);
	  n=j+1;
	}
    }
}
void Grid::draw(float & wait)
{
  int init=0;
  initdraw (this,init,wait);
  
  do {
    reffecran();
    show();
  } while (loopdraw(wait,this));  
}
void Grid::show()
{
  showgraphic();// modif FH 270408
  couleur(1);
  SetColorTable(2+12);
  int i, j;
  couleur(1);
  for(i=0; i<nt; i++)
    {
      rmoveto(v[no(t[i].v[2])].x, v[no(t[i].v[2])].y);
      for( j=0; j<3; j++)
	rlineto(v[no(t[i].v[j])].x,v[no(t[i].v[j])].y);
    }
  for(i=0; i<ne; i++)
    if(e[i].where)
      {  couleur(e[i].where+1);
	rmoveto(e[i].in->x,e[i].in->y);
	rlineto(e[i].out->x,e[i].out->y);
      }
   
}


void contour(Grid& t, int coul)
{
  showgraphic();// modif FH 270408
  couleur(1);
  
  int i = 0;
  for(;i<t.ne;i++)
    if(t.e[i].where) 
      if (InRecScreen(t.e[i].in->x ,t.e[i].in->y ,
		      t.e[i].out->x,t.e[i].out->y) )
        {  
          couleur(t.e[i].where+1);
          rmoveto(t.e[i].in->x ,t.e[i].in->y );
          rlineto(t.e[i].out->x,t.e[i].out->y);
        }
}

void     FillRect(float x0,float y0, float x1, float y1)
{
  float r[8];
  r[0]=x0;r[1]=y0;
  r[2]=x1;r[3]=y0;
  r[4]=x1;r[5]=y1;
  r[6]=x0;r[7]=y1;
  fillpoly(4,r);
}


void equpot(Grid& t,float* f,int nl, int waitm,float fm,float xfm)
{
  showgraphic();// modif FH 270408`
  SetColorTable(nl+2);
  float qp0[4], qp1[4];
  int ns=t.nv, nt=t.nt;
  float xln, xf,fi,fj,xlam;
  int im,i,k,l,ik,jk;
  if (fm>xfm) 
    {
      if(verbosity>1)
	cout <<"\t\t"  << " init plot " << endl;
      Init(t, "o");

      fm= f[0];				/*   search fmin and fmax */
      xfm=fm;
      for(i=1;i<=ns;i++)
	if((f[i-1]) > fm) fm=f[i-1];
	else if((f[i-1]) < xfm) xfm=f[i-1];
    }
  if(fabs(fm-xfm)<1.0e-15) nl = 1;
  float xmin,xmax,ymin,ymax;
  getcadre(xmin,xmax,ymin,ymax);
  float xleft = xmax - (xmax-xmin)*0.1;
  float ytop  = ymax;
  float ydelta = (ymax-ymin)/40;
  ydelta=GetHeigthFont();
  xleft = xmax - 6*ydelta;  
  ytop -= ydelta;
  ytop -= ydelta;
  
  for(l=1;l<=nl;l++)    		/*    loop on the level curves */
    {
      if(nl == 1) xln=0.5F;
      else xln=(l-1.F)/(nl-1.F);
      xf=xfm+(fm-xfm)*xln;

      couleur(l+1);
      FillRect(xleft+ydelta/8.,ytop+ydelta/8.,xleft+ydelta*7./8.,ytop+ydelta*7./8.);
      rmoveto(xleft+ydelta*1.4,ytop+ydelta/4);
      // -- 
      char buf[30];
      sprintf(buf,"%g",xf);
      couleur(1);
      plotstring(buf);
      ytop -= ydelta;

    
      for(k=0;k<nt;k++)			/*   loop on each t.tiangle */
	{
	  bTriangle &tk = t.t[k];
	  if (InRecScreen( Min3(tk.v[0]->x,tk.v[1]->x,tk.v[2]->x),
			   Min3(tk.v[0]->y,tk.v[1]->y,tk.v[2]->y),
			   Max3(tk.v[0]->x,tk.v[1]->x,tk.v[2]->x),
			   Max3(tk.v[0]->y,tk.v[1]->y,tk.v[2]->y))) 
	    {
	      im=0;
	      for(i=0;i<=2;i++)
		{
		  int  j = (i+1)%3;
		  ik = t.no(tk.v[i]);
		  jk = t.no(tk.v[j]);
		  bVertex &vi = *tk.v[i];
		  bVertex &vj = *tk.v[j];
		  fi=(f[ik]);
		  fj=(f[jk]);
		  if(((fi<=xf)&&(fj>=xf))||((fi>=xf)&&(fj<=xf)))
		    {
		      if (fabs(fi-fj)<=0.1e-10) 	/* one side must be drawn */
			{
			  couleur(l+1);
			  rmoveto(vi.x,vi.y);
			  rlineto(vj.x,vj.y);
			}
		      else
			{
			  xlam=(fi-xf)/(fi-fj);
			  qp0[im]  =vi.x*(1.F-xlam)+vj.x*xlam;
			  qp1[im++]=vi.y*(1.F-xlam)+vj.y*xlam;
			}
		    }
		}    
	      if (im>=2) 				/*    draw one segment */
		{
		  couleur(l+1);
		  rmoveto(qp0[0],qp1[0]);
		  rlineto(qp0[1],qp1[1]);
		}
	    }
	} 
    }
  
  contour(t,11);
  if(verbosity>1)
    cout << "\t min="<<fm<<"\n\t max="<<xfm<<endl;
  rattente(waitm);
}

void equpotP0(Grid& t,float* f,int nl, int waitm,float fmin,float fmax)
{
  showgraphic();// modif FH 270408`
  SetColorTable(nl+2);
  int nt=t.nt;
  int i,k;
  if (fmax < fmin) {
    if(verbosity>2)
      cout <<"\t\t"  << " init plot " << endl;
    Init(t, "o");

    fmin= f[0];				/*   search fmin and fmax */
    fmax=fmin;
    for(i=1;i<=nt;i++)
      if((f[i-1]) > fmax) fmax=f[i-1];
      else if((f[i-1]) < fmin) fmin=f[i-1];
  }
    
  float coef = nl/Max((fmax-fmin),1e-20F);
  //    cout << " coef " << coef << " " << fmin << " " << fmax <<  endl;
  for(k=0;k<nt;k++)			/*   loop on each t.tiangle */
    {
      bTriangle &tk = t.t[k];
       
      if (InRecScreen( Min3(tk.v[0]->x,tk.v[1]->x,tk.v[2]->x),
		       Min3(tk.v[0]->y,tk.v[1]->y,tk.v[2]->y),
		       Max3(tk.v[0]->x,tk.v[1]->x,tk.v[2]->x),
		       Max3(tk.v[0]->y,tk.v[1]->y,tk.v[2]->y))) 
	{
	  float v = f[k];
	  //        cout << "  k " << k << " " << v << " c=" << (int) Min(nl,Max(0,(int) ((v-fmin)*coef)))  <<endl;
	  couleur(2+ Min(nl-1,Max(0,(int) ((v-fmin)*coef))));
	  float P[]={tk.v[0]->x,tk.v[0]->y,tk.v[1]->x,tk.v[1]->y,tk.v[2]->x,tk.v[2]->y};        
	  fillpoly(3,P);
	  couleur(0);
	  rmoveto(P[0],P[1]);
	  rlineto(P[2],P[3]);
	  rlineto(P[4],P[5]);
	  rlineto(P[0],P[1]);
        }     
    }
  
  contour(t,11);
  if(verbosity>1)
    cout << "\t min="<<fmin<<"\n\t max="<<fmax<<endl;
  rattente(waitm);
}

void showbdy(long nbs,float* cr, long nba, long* arete, float* hh,int * ng,int * ngf)
{
  showgraphic();// modif FH 270408
  couleur(1);
  int i;
  float d;
  for(i= 0; i<nbs; i++)
    {
      couleur(1+ng[i]);

      d = 0.1F*hh[i];
      rmoveto( cr[2*i]-d, cr[2*i+1]-d);
      rlineto(  d+cr[2*i], cr[2*i+1]-d);
      rlineto(  d+cr[2*i], d+cr[2*i+1]);
      rlineto(  cr[2*i]-d, d+cr[2*i+1]);
      rlineto(  cr[2*i]-d, cr[2*i+1]-d);
    }	
  for(i= 0; i<nba; i++)
    {
      couleur(1+ngf[i]);
      rmoveto( cr[2*arete[2*i]], cr[2*arete[2*i]+1]);
      rlineto( cr[2*arete[2*i+1]], cr[2*arete[2*i+1]+1]);
    }	
}

void graph3d(Grid& t,float* f, int waitm)
{
  showgraphic();// modif FH 270408
  /*
    float xmn=t.v[0].x, xmx=xmn, ymn=t.v[0].y, ymx=ymn, x, y;
    int ns=t.nv, nt=t.nt;
    float fm, xfm ,xlam;
    GRAPH  p1, p2, p3, n;
    extern GLfloat  wAngleY;
    extern GLfloat  wAngleX;
    extern GLfloat  wAngleZ;
    for( i=1; i<t.nv; i++ ) {
    x=t.v[i].x; y=t.v[i].y;
    xmx=x>xmx?x:xmx;
    xmn=xmn>x?x:xmn;
    ymx=y>ymx?y:ymx;
    ymn=ymn>y?y:ymn;
    }
    fm= f[0];            
    xfm=fm;
    for(i=1;i<=ns;i++)
    if((f[i-1]) > fm) fm=f[i-1];
    else if((f[i-1]) < xfm) xfm=f[i-1];
    oglm::InitView(xmx,xmn,ymx,ymn,fm,xfm);
    ::glRotatef(wAngleY, 0.0f, 1.0f, 0.0f);
    ::glRotatef(wAngleX, 1.0f, 0.0f, 0.0f);
    ::glRotatef(wAngleZ, 0.0f, 0.0f, 1.0f);
    ::glScalef(wScale,wScale,wScale);
    oglm::DrawBox(xmx,xmn,ymx,ymn,fm,xfm);
    oglm::SetMaterial();
    for(k=1;k<=nt;k++) {        
    ::glBegin(GL_POLYGON);
    p1.x = t.v[t.no(t.t[k-1].v[0])].x;
    p2.x = t.v[t.no(t.t[k-1].v[1])].x;
    p3.x = t.v[t.no(t.t[k-1].v[2])].x;
    p1.y = t.v[t.no(t.t[k-1].v[0])].y;
    p2.y = t.v[t.no(t.t[k-1].v[1])].y;
    p3.y = t.v[t.no(t.t[k-1].v[2])].y;
    p1.f = f[t.no(t.t[k-1].v[0])];
    p2.f = f[t.no(t.t[k-1].v[1])];
    p3.f = f[t.no(t.t[k-1].v[2])];
    oglm::CalcNormal(&p1, &p2, &p3, &n);
    ::glNormal3f((GLfloat)n.x, (GLfloat)n.y, (GLfloat)n.f);
    ::glVertex3f(p1.x, p1.y, p1.f);
    ::glVertex3f(p2.x, p2.y, p2.f);
    ::glVertex3f(p3.x, p3.y, p3.f);
    ::glEnd();
    }
    oglm::EndView();
  */}


