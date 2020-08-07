#include "kdtree_2.h"

#include <cstdlib>  // rand
#include <iostream>
#include <vector>

using namespace std;

// Example
int main (int , char* [] )
{
  unsigned int i=0, j=0;

  cout << "Generation of the vecKey list" << endl;

  vector<vector<double> > vecKey(100);

  for(i=0; i<vecKey.size(); i++)
    {
      vecKey[i].resize(3);

      for(j=0; j<vecKey[i].size(); j++)
	{
	  vecKey[i][j] = rand()/(double)(RAND_MAX);
	} /* End For */
    } /* End For */

  cout << "Construction of the tree" << endl;

  typedef kdtree kdtree_type;

  kdtree_type kd_tree;

  kd_tree.setDimension(3);

  kd_tree.run(vecKey.begin(), vecKey.end());

  cout << "Searching for the 10 nearest neighbors of the first point" << endl;

  vector<vector<double> > vecAns = kd_tree.search_n( 10, vecKey[0]);

  cout << "Displaying the nearest neighbors" << endl;

  for(i=0; i<10; i++)
    {
      for(j=0; j<vecAns[i].size(); j++)
	{
	  cout << vecAns[i][j] << ", ";
	} /* End For */
      cout << endl;
    } /* End For */

  return 0;
}
