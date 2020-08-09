#ifndef GRAPH_H
#define GRAPH_H

//file graph.h
class Grid;
struct OnNode;
typedef OnNode * OnList ;
void graph3d(float* f, int waitm);
void showtriangulation(Grid& t, int waitm);
//void equpot(Grid& t, float* f, int nl, int waitm);
void equpot(Grid& t,float* f,int nl, int waitm,float fmax,float fmin);
void equpotP0(Grid& t,float* f,int nl, int waitm,float fmax,float fmin);
void showbdy(long nbs,float* cr, long nba, long* arete, float* hh, int * ng,int * ngf);
void graph3d(Grid& t, float* f, int waitm);  // ohtsuka
void contour(Grid& t, int coul);
int loopdraw (float & wait,OnList list);
int loopdraw (float & wait,Grid *g);
void initdraw ( Grid *g,int &init,float  wait=-1);
void initdraw ( OnList list, float  wait=-1);
void initdraw ( frontiere *f,float  wait);
#endif
