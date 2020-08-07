#include "kdtree.h"

#include <cstdlib>  // rand
#include <iostream>
#include <vector>

using namespace std;

// A key with three random doubles.
class Key 
{
public:
  vector<double> key;
};

// Key adapter for the kdtree.
template < class KEY >
class kdtree_key_adapter 
{
private:
  class less_impl;
  class assign_impl;
  
public:
  static less_impl less( int i ) { return less_impl(i); }
  
  static assign_impl assign( int i ) { return assign_impl(i); }
  
  static float difference( const KEY& k1, const KEY& k2 )
  { return difference_impl( k1, k2 ); }
  
private:
  class less_impl 
  {
  public:
    less_impl( int index ) : i(index) {}
    bool operator()( const KEY& k1, const KEY& k2 ) const
    { return k1.key[i] < k2.key[i]; }
    
  private:
    const int i;
  };

  class assign_impl 
  {
  public:
    assign_impl( int index ) : i(index) {}
    void operator()( KEY& k1, const KEY& k2 ) const
    { k1.key[i] = k2.key[i]; }

  private:
    const int i;
  };

  static float difference_impl( const KEY& k1, const KEY& k2 )
  {
    float dres = 0.0;

    for(unsigned int i=0; i<k1.key.size(); i++)
      {
	dres += (k1.key[i] - k2.key[i]) * (k1.key[i] - k2.key[i]);
      } /* End For */

    return dres;
  }
};

// Example
int main (int , char* [] )
{
  unsigned int i=0, j=0;

  cout << "Generation of the vecKey list" << endl;

  vector<Key> vecKey(100);

  for(i=0; i<vecKey.size(); i++)
    {
      vecKey[i].key.resize(3);

      for(j=0; j<vecKey[i].key.size(); j++)
	{
	  vecKey[i].key[j] = rand()/(double)(RAND_MAX);
	} /* End For */
    } /* End For */

  cout << "Construction of the tree" << endl;

  typedef kdtree<Key, kdtree_key_adapter<Key> > kdtree_type;

  kdtree_type kd_tree;

  kd_tree.setDimension(3);

  kd_tree.run(vecKey.begin(), vecKey.end());

  cout << "Searching for the 10 nearest neighbors of the first point" << endl;

  vector<Key> vecAns = kd_tree.search_n( 10, vecKey[0]);

  cout << "Displaying the nearest neighbors" << endl;

  for(i=0; i<10; i++)
    {
      for(j=0; j<vecAns[i].key.size(); j++)
	{
	  cout << vecAns[i].key[j] << ", ";
	} /* End For */
      cout << endl;
    } /* End For */

  return 0;
}
