#ifndef KDTREE_STATIC_H
#define KDTREE_STATIC_H

#include <vector>

using namespace std;

bool comp_less(int i, const vector<double> & k1, const vector<double> & k2);
void comp_assign(int i, vector<double> & k1, const vector<double> & k2);
double comp_difference(const vector<double> & k1, const vector<double> & k2);

class COMP2
{
private:
  class less_impl;
  
public:
  COMP2() : i(0) {};
  static less_impl less(int i)
    {
      return less_impl(i);
    }
  
 private:
  class less_impl 
  {
  public:
    less_impl( int index ) : i(index) {}
      bool operator()(const vector<double> k1, const vector<double> k2) const { return k1[i] < k2[i]; }
      
  private:
      const int i;
  };
  
 private:
  const int i;
};
#endif
