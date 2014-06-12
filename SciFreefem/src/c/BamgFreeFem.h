#ifndef FREEFEMBAMG_H_
#define FREEFEMBAMG_H_
Geometry * frontiere2Geometry(const frontiere & fr);
class frontiere {
 public:
  int nbp;			// number of boundary points
  long nbs;	        // number of edges
  long nbsd;		// number of each connected component
  float *xy;		// cordinates (for point i, x in 2*i , y in 2*i+1)
  int *ng;			// boundary id number for each point
  int *ngf;		// boundary id number for each segment
  long *s;			// edges (edge i: beginning in 2*i, end in 2*i +1)
  long *sd;		// an edge for each connected component
  float *hh;		// weight of points (to insert inner points)
  float xmin,xmax,ymin,ymax; // add FH to compute a related epsilon 
  int initboundingbox;
  int step;
  float epsilon;
  frontiere() 
    {
      step=initboundingbox=nbp =0;
      nbs = nbsd = 0;
      sd = new long[90];
      xy=0;
      ng=0;
      ngf=0;
      s=0;
      hh=0;
    }
  void init()
  { 
    if(xy) delete [] xy;   xy = 0;
    if(ng) delete [] ng;   ng = 0;
    if(s) delete [] s;     s = 0;
    if(hh) delete [] hh;   hh = 0;
    if(ngf) delete [] ngf; ngf = 0;
    step=initboundingbox=nbp = nbs = nbsd = 0;
    xy=0;
    ng=0;
    s=0;
    hh=0;
  } 
  int addPoint(float x, float y, int ng);
  void addSegment(int n1, int n2, int label);
  void save(const char* filename) const; //OP 97
};
#endif
Int4 FindTriangle(Triangles &Th, Real8 x, Real8 y, double* a);
