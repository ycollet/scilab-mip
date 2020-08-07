// GammaTest.cpp
// Author: Y. Collette
// Date: 19/11/2005
// Renault - DR - 64240

#include <vector>
#include <limits>
#include <utility>
#include <iostream>
#include <fstream>
#include <iterator>
#include <memory>
#include <queue>
#include <math.h>

//#include <kdtree_2.h>

#include <GammaTest.h>

#ifndef UP
#define UP "\033[A"
#endif

using namespace std;

#include <kdtree_2.h>

// struct queue_item
// {
//   queue_item(const vector<double> & answer, double distance);
//   bool operator<(const queue_item& rhs);
//   vector<double> m_answer;
//   double m_distance;
// };

// class node
// {
//  public:
//   virtual void search(int count,
// 		      const vector<double> & query,
// 		      const vector<double> & upper_bound,
// 		      const vector<double> & lower_bound,
// 		      priority_queue<queue_item> * results);
//   virtual ~node();
// };

// class interior_node : public node
// {
//  public:
//   interior_node(int dim, const vector<double> & partition, auto_ptr<node> l, auto_ptr<node> r);
//   virtual void search(int count,
// 		      const vector<double> & query,
// 		      const vector<double> & upper_bound,
// 		      const vector<double> & lower_bound,
// 		      priority_queue<queue_item> * ) const;
//   bool bounds_overlap(const vector<double> & query,
// 		      const queue_item & top,
// 		      int dim,
// 		      const vector<double> & bound ) const;
//  private:
//   interior_node(const interior_node &);
//   interior_node & operator=(const interior_node &);
//   const int m_nsort_dim;
//   const vector<double> m_partition;
//   const auto_ptr<node> m_pLChild;
//   const auto_ptr<node> m_pRChild;
// };

// class leaf_node : public node
// {
//  public:
//   typedef vector<vector<double> >::const_iterator const_iterator;
//   leaf_node(const_iterator first, const_iterator last);
//   virtual void search(int count,
// 		      const vector<double> & query,
// 		      const vector<double> & upper_bound,
// 		      const vector<double> & lower_bound,
// 		      priority_queue<queue_item> *) const;
//  private:
//   leaf_node(const leaf_node &);
//   leaf_node & operator=(const leaf_node &);
//   const vector<vector<double> > m_vec;
// };

// class kdtree 
// {
// public:
//   kdtree();
//   kdtree(vector<vector<double> >::iterator first, vector<vector<double> >::iterator last);
//   vector<vector<double> > search_n(int count, const vector<double> & query) const;
//   void setDimension(int Dimension);
//   int  getDimension();
//   void setBucket(int Bucket);
//   int  getBucket();
//   void run(vector<vector<double> >::iterator first, vector<vector<double> >::iterator last);
// private:
//   kdtree(const kdtree &);
//   kdtree& operator=(const kdtree &);
//   auto_ptr<node> build_tree(vector<vector<double> >::iterator first,
// 			    vector<vector<double> >::iterator last,
// 			    int level);
//   vector<double> upper_bound(vector<vector<double> >::const_iterator first,
// 			     vector<vector<double> >::const_iterator last) const;
//   vector<double> lower_bound(vector<vector<double> >::const_iterator first,
// 			     vector<vector<double> >::const_iterator last) const;

//   int _Dimension;
//   int _Bucket;
//   auto_ptr<node> m_pRoot;
//   vector<double> m_upper_bound;
//   vector<double> m_lower_bound;
// };


void GammaTest_Estimator(vector<vector<double> > & InputData,
			 vector<vector<double> > & OutputData,
			 vector<double> & Y_Output,
			 unsigned int p,
			 unsigned int Measure,
			 bool useVersion2,
			 unsigned int Starting_n,
			 int BucketSize)
{
  typedef kdtree kdtree_type;

  vector<vector<double> > AuxList;
  double                  AuxDeltaValue, AuxGammaValue;
  vector<vector<double> > AuxPair;
  vector<double>          ListAuxDelta, ListAuxGamma;
  unsigned int            i, j, k;
  bool                    PointFound = false;
  double                  m_x = 0.0, m_y = 0.0, m_dx2 = 0.0, m_dxdy = 0.0;
  double                  m = 0.0, Delta = 0.0;
  kdtree_type           * kd_tree;
  vector<double>          qry;
  vector<vector<double> > vecAns;

  vector<vector<double> >       A, B, C, Lambda;
  vector<vector<vector<unsigned int> > > S;
  vector<double>                Mu, Nu, ListDelta, ListGamma;
  double                        Snew_sqr, a = 0.0, b = 0.0, AuxDbl, AuxDbl2;

  vector<double>::iterator IterDbl;
  vector<unsigned int>::iterator IterUI;

  kd_tree = new kdtree_type;

  kd_tree->setDimension(InputData[0].size());
  kd_tree->setBucket(BucketSize);

  cout << "GammaTest: construction of the data file" << endl;

  // Fusion Inputs + Output
  AuxPair.resize(0);
  for(i=0; i<InputData.size(); i++)
    {
      qry.resize(0);
      for(j=0; j<InputData[i].size(); j++)
	{
	  qry.push_back(InputData[i][j]);
	} /* End For */
      qry.push_back(OutputData[i][Measure]);
      AuxPair.push_back(qry);
    } /* End For */

  // Construction of the kd_tree

  cout << "GammaTest_Estimator: construction of the KD tree" << endl;

  kd_tree->run(AuxPair.begin(), AuxPair.end());

  ListDelta.reserve(p);
  ListGamma.reserve(p);
  ListDelta.resize(p, 0.0);
  ListGamma.resize(p, 0.0);

  S.resize(InputData.size(), vector<vector<unsigned int> >(p, vector<unsigned int>(0)));
  A.resize(InputData.size(), vector<double>(p, 0.0));
  B.resize(InputData.size(), vector<double>(p, 0.0));
  C.resize(InputData.size(), vector<double>(p, 0.0));
  Lambda.resize(InputData.size(), vector<double>(p, 0.0));
  Mu.resize(InputData.size(), 0.0);
  Nu.resize(InputData.size(), 0.0);

  Y_Output.resize(InputData.size(), 0.0);

  cout << endl;
  
  for(i=0; i<InputData.size(); i++)
    {
      // Calcul de la liste des distances par rapport au point i
      cout << UP;
      cout << "GammaTest_Estimator: search for " << p << " neighbors of point " << i << endl;
      
      vecAns = kd_tree->search_n(p, InputData[i]);
            
      // Calcul de Delta et Gamma
      for(j=Starting_n; j<vecAns.size(); j++)
	{
	  // Calcul de Delta
	  AuxDeltaValue = 0.0;
	  AuxGammaValue = 0.0;
	  
	  for(k=0; k<vecAns[j].size() - 1; k++)
	    {
	      AuxDeltaValue += (vecAns[j][k] - InputData[i][k]) * (vecAns[j][k] - InputData[i][k]);
	    } /* End For */
	  ListDelta[j] += AuxDeltaValue;
	  
	  // Calcul de Gamma
	  AuxGammaValue = (vecAns[j][vecAns[j].size() - 1] - OutputData[i][Measure]) *
	    (vecAns[j][vecAns[j].size() - 1] - OutputData[i][Measure]);
	  ListGamma[j] += AuxGammaValue;
	} /* End For */
    } /* End For */

  cout << "GammaTest_Estimator: normalisation of the lists" << endl;
  
  // Normalisation de Delta et Gamma
  for(i=Starting_n; i<p; i++)
    {
      ListDelta[i] = ListDelta[i] / (double)(InputData.size());
      ListGamma[i] = ListGamma[i] / (double)(2*InputData.size());
    } /* End For */
  
  if (useVersion2)
    {
      ListAuxDelta.resize(ListDelta.size(), 0.0);
      ListAuxGamma.resize(ListGamma.size(), 0.0);
      
      for(i=Starting_n; i<p; i++)
	{
	  for(j=Starting_n; j<i; j++)
	    {
	      ListAuxDelta[i] += ListDelta[j];
	      ListAuxGamma[i] += ListGamma[j];
	    } /* End For */
	  
	  ListAuxDelta[i] /= (i+1.0);
	  ListAuxGamma[i] /= (i+1.0);
	} /* End For */
      
      for(i=Starting_n; i<p; i++)
	{
	  ListDelta[i] = ListAuxDelta[i];
	  ListGamma[i] = ListAuxGamma[i];
	} /* End For */
      
      ListAuxDelta.resize(0);
      ListAuxGamma.resize(0);
  
      // Calcul des coeff de regression entre Delta et Gamma: Gamma = A.Delta + B
      
      cout << "GammaTest_Estimator: computation of the regression line" << endl;
      
      m_x    = 0.0;
      m_y    = 0.0;
      m_dx2  = 0.0;
      m_dxdy = 0.0;
      m      = 0.0;
      
      for(i=Starting_n; i<ListDelta.size(); i++)
	{
	  m      += 1;
	  m_x    += ListDelta[i];
	  m_y    += ListGamma[i];
	  m_dx2  += ListDelta[i]*ListDelta[i];
	  m_dxdy += ListDelta[i]*ListGamma[i];
	} /* End For */
      
      Delta = m*m_dx2 - m_x*m_x;
      
      /* In terms of y = a + b x */
      
      a = (m_dx2*m_y - m_x*m_dxdy)/Delta;
      b = (m*m_dxdy - m_x*m_y)/Delta;
    } /* End If */

  Snew_sqr = Delta / (double)p;

  // Calcul of the estimations
  
  // Computation of S(p)

  cout << "GammaTest_Estimator:: Computation of S(p)" << endl << endl;

  for(i=0; i<InputData.size(); i++)
    {
      cout << UP << "Point " << i << " / " << InputData.size() << endl;

      for(j=0; j<InputData.size(); j++)
	{
	  if (i==j) continue;
	  vecAns = kd_tree->search_n(p, InputData[j]);

	  PointFound = false;

	  for(k=Starting_n; k<vecAns.size(); k++)
	    {
	      if ((vecAns[k]==InputData[i])&&!PointFound)
		{
		  S[i][k].push_back(i);
		  PointFound = true;
		} /* End If */
	      else if (PointFound)
		{
		  S[i][k].push_back(i);
		} /* End Else If */
	    } /* End For */
	} /* End For */
    } /* End For */

  // Computation of Lambda

  cout << "GammaTest_Estimator:: Computation of Lambda" << endl;

  for(i=0; i<InputData.size(); i++)
    {
      cout << UP << "Point " << i << " / " << InputData.size() << endl;

      vecAns = kd_tree->search_n(p, InputData[i]);

      for(j=Starting_n; j<vecAns.size(); j++)
	{
	  Lambda[i][j] += vecAns[j][vecAns.size()-1];
	  
	  IterUI = find(S[i][j].begin(), S[i][j].end(), i);
	  
	  if (IterUI!=S[i][j].end()) Lambda[i][j] += vecAns[(*IterUI)][vecAns.size()-1];
	  
	  Lambda[i][j] /= (double)(p*InputData.size());
	} /* End For */
    } /* End For */

  // Computation of Mu

  cout << "GammaTest_Estimator:: Computation of Mu" << endl;

  for(i=Starting_n; i<p; i++)
    {
      cout << UP << "Point " << i << " / " << p << endl;

      Mu[i] = 0.0;

      for(j=0; j<InputData.size(); j++)
	{
	  Mu[i] += Lambda[j][i];
	} /* End For */
    } /* End For */

  // Computation of B(p)

  cout << "GammaTest_Estimator:: Computation of B(p)" << endl;

  for(i=0; i<InputData.size(); i++)
    {
      cout << UP << "Point " << i << " / " << InputData.size() << endl;

      for(j=Starting_n; j<p; j++)
	{
	  B[i][j] = Lambda[i][j] - Mu[j];
	} /* End For */
    } /* End For */

  // Computation of Nu

  cout << "GammaTest_Estimator:: Computation of Nu" << endl;

  for(k=0; k<InputData.size(); k++)
    {
      cout << UP << "Point " << k << " / " << InputData.size() << endl;

      Nu[k] = 0.0;

      for(i=Starting_n; i<p; i++)
	{
	  AuxDbl = 0.0;
	  for(j=Starting_n; j<i; j++)
	    {
	      AuxDbl += S[j].size();
	    } /* End For */
	  AuxDbl /= 1.0/(double)i;
	  Nu[k] += AuxDbl;
	} /* End For */
      Nu[k] /= 1.0/(double)p;
      Nu[k] += 1.0;
      Nu[k] *= 1.0/(double)InputData.size();
    } /* End For */

  // Computation of C(p)

  cout << "GammaTest_Estimator:: Computation of C(p)" << endl;

  for(i=0; i<InputData.size(); i++)
    {
      cout << UP << "Point " << i << " / " << InputData.size() << endl;

      for(j=Starting_n; j<p; j++)
	{
	  C[i][j] = ListDelta[j] - b;
	} /* End For */
    } /* End For */

  // Computation of A(p)

  cout << "GammaTest_Estimator:: Computation of A(p)" << endl;

  for(k=0; k<InputData.size(); k++)
    {
      cout << UP << "Point " << k << " / " << InputData.size() << endl;
	  
      AuxDbl2 = 0.0;
      for(i=Starting_n; i<p; i++)
	{
	  AuxDbl = 0.0;
	  for(j=Starting_n; j<i; j++)
	    {
	      AuxDbl += S[k][j].size();
	    } /* End For */
	  AuxDbl /= (1.0/(double)i);
	  AuxDbl2 += AuxDbl;
	} /* End For */
      AuxDbl2 /= (1.0/(double)p);
      
      for(i=Starting_n; i<p; i++)
	{
	  A[k][i] = 0.0;
	  for(j=Starting_n; j<i; j++)
	    {
	      A[k][i] += S[k][j].size();
	    } /* End For */
	  A[k][i] /= (1.0/(double)i);
	  A[k][i] -= AuxDbl2;
	  A[k][i] /= (1.0/(double)InputData.size());
	} /* End For */
    } /* Ennd For */

  // Computation of y

  cout << "GammaTest_Estimator:: Computation of y" << endl;

  for(i=0; i<InputData.size(); i++)
    {
      cout << UP << "Point " << i << " / " << InputData.size() << endl;

      AuxDbl  = 0.0;
      AuxDbl2 = 0.0;
      for(j=Starting_n; j<p; j++)
	{
	  AuxDbl += C[i][j]*B[i][j];
	  AuxDbl2 += C[i][j]*A[i][j];
	} /* End For */
      AuxDbl /= (double)p;
      AuxDbl2 /= (double)p;

      Y_Output[i] = (Mu[i] - b/Snew_sqr*AuxDbl)/(Nu[i] - b/Snew_sqr*AuxDbl2);
    } /* End For */

  if (kd_tree) delete kd_tree;

  A.clear();
  B.clear();
  C.clear();
  S.clear();
  Nu.clear();
  Mu.clear();
  Lambda.clear();
}
