// kdtree.h
// Author: ??
// Downloaded from internet
// Modified by Y. Collette
// Date: 19/11/2005
// Renault - DR - 64240

#ifndef KDTREE_2_H
#define KDTREE_2_H

#include <algorithm>	// nth_element, max_element, min_element, swap, reverse
#include <memory>	// auto_ptr
#include <queue>	// priority_queue
#include <vector>

using namespace std;

struct queue_item
{
  queue_item(const vector<double> & answer, double distance);
  bool operator<(const queue_item & rhs) const;
  vector<double> m_answer;
  double m_distance;
};

class node
{
 public:
  virtual void search(int count,
		      const vector<double> & query,
		      const vector<double> & upper_bound,
		      const vector<double> & lower_bound,
		      priority_queue<queue_item> * results) const = 0;
  virtual ~node();
};

class interior_node : public node
{
 public:
  interior_node(int dim, const vector<double> & partition, auto_ptr<node> l, auto_ptr<node> r);
  virtual void search(int count,
		      const vector<double> & query,
		      const vector<double> & upper_bound,
		      const vector<double> & lower_bound,
		      priority_queue<queue_item> * ) const;
  bool bounds_overlap(const vector<double> & query,
		      const queue_item & top,
		      int dim,
		      const vector<double> & bound ) const;
 private:
  interior_node(const interior_node &);
  interior_node & operator=(const interior_node &);
  const int m_nsort_dim;
  const vector<double> m_partition;
  const auto_ptr<node> m_pLChild;
  const auto_ptr<node> m_pRChild;
};

class leaf_node : public node
{
 public:
  typedef vector<vector<double> >::const_iterator const_iterator;
  typedef vector<vector<double> >::iterator iterator;
  leaf_node(const_iterator first, const_iterator last);
  virtual void search(int count,
		      const vector<double> & query,
		      const vector<double> & upper_bound,
		      const vector<double> & lower_bound,
		      priority_queue<queue_item> *) const;
 private:
  leaf_node(const leaf_node &);
  leaf_node & operator=(const leaf_node &);
  const vector<vector<double> > m_vec;
};

class kdtree 
{
public:
  kdtree();
  kdtree(vector<vector<double> >::iterator first, vector<vector<double> >::iterator last);
  vector<vector<double> > search_n(int count, const vector<double> & query) const;
  void setDimension(int Dimension);
  int  getDimension();
  void setBucket(int Bucket);
  int  getBucket();
  void run(vector<vector<double> >::iterator first, vector<vector<double> >::iterator last);
private:
  kdtree(const kdtree &);
  kdtree& operator=(const kdtree &);
  auto_ptr<node> build_tree(vector<vector<double> >::iterator first,
			    vector<vector<double> >::iterator last,
			    int level);
  vector<double> upper_bound(vector<vector<double> >::const_iterator first,
			     vector<vector<double> >::const_iterator last) const;
  vector<double> lower_bound(vector<vector<double> >::const_iterator first,
			     vector<vector<double> >::const_iterator last) const;

  int _Dimension;
  int _Bucket;
  auto_ptr<node> m_pRoot;
  vector<double> m_upper_bound;
  vector<double> m_lower_bound;
};

#endif
