#include <kdtree_static.h>

bool comp_less(int i, const vector<double> & k1, const vector<double> & k2)
{
  return k1[i] < k2[i];
}

void comp_assign(int i, vector<double> & k1, const vector<double> & k2)
{
  k1[i] = k2[i];
}

double comp_difference(const vector<double> & k1, const vector<double> & k2) 
{
  unsigned int i = 0;
  double dres = 0.0;
  
  for(i=0; i<k1.size()-1; i++)
    {
      dres += (k1[i] - k2[i]) * (k1[i] - k2[i]);
    } /* End For */
  
  return dres;
}
