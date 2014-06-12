// kdtree.h
// Author: ??
// Downloaded from internet
// Modified by Y. Collette
// Date: 19/11/2005
// Renault - DR - 64240

#include <algorithm>	// nth_element, max_element, min_element, swap, reverse
#include <memory>	// auto_ptr
#include <queue>	// priority_queue
#include <vector>

using namespace std;

#include <kdtree_2.h>
#include <kdtree_static.h>

/*!
  \brief A k-d tree is a binary search tree for k-dimensional data.

  I created this as an example of how to use the C++ Standard Template
  Library (STL). <b>I do not recommend that you use this k-d tree for
  anything other than a learning example.</b>

  This implementation is based on another C++ implementation of a
  k-d tree that I found. This code and the code upon which it is
  based are functionaly very similar. However, this implementation
  is only about one third the size of the original, as measured in
  Source Lines of Code (SLOC). Even more dramatic is that this
  implementation has a much lower cyclomatic complexity. The
  methods in this implementation have a maximum complexity around 4
  with an average close to 2; the original had a maximum of 24 with an
  average around 7. Finally, this code constructs the k-d tree much
  more quickly than did the original.

  For an introduction to k-d trees see: <i>CACM</i> <b>18</b> (1975),
  509-517 or <i>IEEE Transactions</i> <b>SE-5</b> (1979), 333-340.

  This k-d tree does not support many of the operations you would
  expect from a search tree. For example, you can not insert or
  delete keys. This class only allows you to construct a tree from
  an existing collection of keys and then search for the nearest N keys
  in the k-d tree.

  To use this k-d tree, the keys it stores require a certain interface.
  The following example shows what is required. The example achieves
  the interface through an Adapter, but it could be built-in.
*/


/*!
  \brief As keys are found during a search, they are tracked
  with this data structure.
*/

/*!
  Constructs a key to be stored in the list of close keys.
*/

queue_item::queue_item(const vector<double> & answer, double distance) : m_answer(answer), m_distance(distance)
{
}

/*!
  Comparison operator required by the priority_queue.
  Sorts keys based on distance from query value. This
  will result in a queue such that queue.top() is
  furthest from the query key.
*/

bool queue_item::operator<(const queue_item & rhs) const
{ 
  return m_distance < rhs.m_distance;
}
  
/*!
  \brief Interface for search tree nodes.
  
  This class, and the children interior_node and leaf_node
  implement the <i>Composite</i> pattern as described in
  <i>Design Patterns: Elements of Reusable Object-Oriented
  Software</i>. They play the roles of component, composite
  and leaf, respectively.
*/

/*!
  Destroys a node.
*/

node::~node()
{
}

/*!
  \brief Interior nodes in the search tree.
  
  Interior nodes store the partition value that divides the left
  and right sub-trees.  Keys in the left sub-tree will be less
  than or equal to \ref m_partition in the \ref m_nsort_dim
  dimension.  Similarly, keys in the right sub-tree will be
  greater than or equal.
  
  This class plays the role of composite in the Composite
  pattern.
  
  \sa node
*/

interior_node::interior_node(int dim,
			     const vector<double> & partition,
			     auto_ptr<node> l,
			     auto_ptr<node> r)
  : m_nsort_dim(dim), m_partition(partition), m_pLChild(l), m_pRChild(r)
{
}

void interior_node::search(int count,
			   const vector<double> & query,
			   const vector<double> & upper,
			   const vector<double> & lower,
			   priority_queue<queue_item> * results) const
{
  if (comp_less(m_nsort_dim, m_partition, query))
    {
      // Search right.
      vector<double> new_lower(lower);
      
      comp_assign(m_nsort_dim, new_lower, m_partition);
      m_pRChild->search(count, query, upper, new_lower, results);
      
      // If query too close to partition value, search other side.
      if (bounds_overlap(query, results->top(), m_nsort_dim, new_lower))
	{
	  vector<double> new_upper(upper);

	  comp_assign(m_nsort_dim, new_upper, m_partition);
	  m_pLChild->search(count, query, new_upper, lower, results);
	}
    }
  else
    {
      vector<double> new_upper(upper);

      comp_assign(m_nsort_dim, new_upper, m_partition);
      m_pLChild->search(count, query, new_upper, lower, results);

      if (bounds_overlap(query, results->top(), m_nsort_dim, new_upper))
	{
	  vector<double> new_lower(lower);

	  comp_assign(m_nsort_dim, new_lower, m_partition);
	  m_pRChild->search(count, query, upper, new_lower, results);
	}
    }
}

/*!
  Pruning heuristic. If true then search should continue.
  
  Returns true if \a query is closer to \a bound than \a top.
*/

bool interior_node::bounds_overlap(const vector<double> & query,
				   const queue_item & top,
				   int dim,
				   const vector<double> & bound) const
{
  // Translate query to partition boundary.
  vector<double> qt(query);

  comp_assign(dim, qt, bound);
  
  /*	If query closer to boundary than to worst nearest answer
    then keep searching
  */
  if (comp_difference(qt, query) <= top.m_distance)
    {
      return true;
    }

  return false;
}

/*!
  \brief Leaf nodes in the search tree.
  
  Leaf nodes contain the key data.
  
  This class plays the role of leaf in the Composite pattern.
  
  \sa node
*/

/*!
  Construct a leaf node by copying the data from \a first
  to \a last.
*/
leaf_node::leaf_node(const_iterator first, const_iterator last) : m_vec(first, last)
{
}
  
/*!
  \todo Remove the static_cast in this method.
*/
void leaf_node::search(int matches,
		       const vector<double> & query,
		       const vector<double> & /*upper*/,
		       const vector<double> & /*lower*/,
		       priority_queue<queue_item> * results) const
{
  leaf_node::const_iterator last = m_vec.end();
  
  for(leaf_node::const_iterator iter=m_vec.begin(); iter!=last; ++iter)
    {
      results->push(queue_item(*iter, comp_difference(*iter, query)));
      
      while(results->size()>static_cast<priority_queue<queue_item>::size_type>(matches))
	{
	  results->pop();
	}
    }
}

kdtree::kdtree() : _Dimension(1), _Bucket(5), m_pRoot(NULL)
{
}

kdtree::kdtree(vector<vector<double> >::iterator first, vector<vector<double> >::iterator last) : _Dimension(1), _Bucket(5), m_pRoot(NULL)
{
}
    
void kdtree::setDimension(int Dimension)
{
  _Dimension = Dimension;
}

int  kdtree::getDimension() 
{
  return _Dimension;
}

void kdtree::setBucket(int Bucket)
{
  _Bucket = Bucket;
}

int  kdtree::getBucket()
{
  return _Bucket;
}

void kdtree::run(vector<vector<double> >::iterator first, vector<vector<double> >::iterator last) 
{
  m_pRoot       = build_tree( first, last, 0);
  m_upper_bound = upper_bound(first, last);
  m_lower_bound = lower_bound(first, last);
}

/*!
  Recursively build the search tree.
	
  Elements from \a first to \a last are copied into the tree.

  Returns the root of the tree.
*/
auto_ptr<node> kdtree::build_tree(vector<vector<double> >::iterator first,
				  vector<vector<double> >::iterator last,
				  int level)
{
  // Base case to terminate recursion.
  if (distance(first, last) <= _Bucket)
    {
      return auto_ptr<node>(new leaf_node(first, last));
    }
  
  /*	Partition the data into a left and right side. Organise the
    data to maintain the search tree invariant:
    left[i] <= key && right[i] >= key
  */
  const int nSortDim = (level + 1) % _Dimension;
  
  vector<vector<double> >::iterator mid = first - 1 + (distance(first, last) / _Bucket + 1) / 2 * _Bucket;
  nth_element(first, mid, last, COMP2::less(nSortDim));
  vector<double> partition(*mid);
  
  // Recursively build sub-trees from partitions.
  return auto_ptr<node>(new interior_node(nSortDim,	
					  partition,
					  build_tree(first, mid + 1, level + 1),
					  build_tree(mid + 1, last, level + 1)));
}

/*!
  For keys in the range \a first to \a last, find the upper bound for
  each dimension.
  
  Returns a key containing the maximum values in each dimension.
*/
vector<double> kdtree::upper_bound(vector<vector<double> >::const_iterator first,
				   vector<vector<double> >::const_iterator last) const
{
  vector<double> u_bound(*first);
  
  for(int i=0; i<_Dimension; ++i)
    {
      vector<vector<double> >::const_iterator dim_max = max_element(first, last, COMP2::less(i));
      comp_assign(i, u_bound, *dim_max);
    }
  
  return u_bound;
}

/*!
  For keys in the range \a first to \a last, find the lower bound for
  each dimension.
	
  Returns a key containing the minimum values in each dimension.
*/
vector<double> kdtree::lower_bound(vector<vector<double> >::const_iterator first,
				   vector<vector<double> >::const_iterator last) const
{
  vector<double> l_bound(*first);
  
  for(int i=0; i<_Dimension;++i)
    {
      vector<vector<double> >::const_iterator dim_min = min_element(first, last, COMP2::less(i));
      comp_assign(i, l_bound, *dim_min);
    }
  
  return l_bound;
}

/*!
  Find \a count closest matches to \a query.
	
  Returns a vector of copies of the closest keys, ordered best-first.

  \todo Describe effect of looking for too many results.
  \todo What if \a count is negative.
*/
vector<vector<double> > kdtree::search_n(int count, const vector<double> & query) const
{
  // Matches will be stored in the queue with the worst on top.
  priority_queue<queue_item> results;
  
  m_pRoot->search(count, query, m_upper_bound, m_lower_bound, &results);
  
  // Copy the results from the queue to a vector, ordered best-first.
  vector<vector<double> > vecResults;
  vecResults.reserve(count);

  while(!results.empty())
    {
      vecResults.push_back(results.top().m_answer);
      results.pop();
    }

  reverse(vecResults.begin(), vecResults.end());

  return vecResults;
}
