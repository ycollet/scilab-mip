// GammaTest.cpp
// Author: Y. Collette
// Date: 19/11/2005
// Renault - DR - 64240

#include <vector>
#include <limits>
#include <utility>
#include <iostream>
#include <fstream>
#include <math.h>

#include <kdtree_2.h>

#include <GammaTest.h>

#ifndef UP
#define UP "\033[A"
#endif

using namespace std;

void GammaTest(vector<vector<double> > & InputData,
	       vector<vector<double> > & OutputData, 
	       vector<double> & ListDelta,
	       vector<double> & ListGamma,
	       vector<double> & ListAllDelta, 
	       vector<double> & ListAllGamma,
	       unsigned int p,
	       double & a,
	       double & b,
	       double Proportion,
	       unsigned int Measure,
	       bool useVersion2,
	       bool useUniquePoint,
	       unsigned int IndexUP,
	       unsigned int Starting_n,
	       bool Use_MTest,
	       int MTest_Min,
	       int MTest_Max,
	       vector<double> & MTest_a,
	       vector<double> & MTest_b,
	       int BucketSize)
{
  typedef kdtree kdtree_type;

  vector<vector<double> > AuxList;
  double                  AuxDeltaValue, AuxGammaValue;
  vector<vector<double> > AuxPair;
  vector<double>          ListAuxDelta, ListAuxGamma;
  vector<double>          ListMTestDelta, ListMTestGamma;
  unsigned int            i, j, k, l;
  bool                    Select = false;
  double                  m_x = 0.0, m_y = 0.0, m_dx2 = 0.0, m_dxdy = 0.0;
  double                  dx = 0.0, dy = 0.0, m = 0.0, Delta = 0.0;
  kdtree_type           * kd_tree;
  vector<double>          qry;
  vector<vector<double> > vecAns;

  kd_tree = new kdtree_type;

  kd_tree->setDimension(InputData[0].size());
  kd_tree->setBucket(BucketSize);

#ifdef DEBUG
  cout << "GammaTest: Input nb points  = " << InputData.size() << endl;
  cout << "GammaTest: Output nb points = " << OutputData.size() << endl;
  cout << "GammaTest: Input Dimension  = " << InputData[0].size() << endl;
  cout << "GammaTest: Output Dimension = " << OutputData[0].size() << endl;
  cout << "GammaTest: Bucket size      = " << kd_tree->getBucket() << endl;
#endif

  // Construction of the data file
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

  kd_tree->run(AuxPair.begin(), AuxPair.end());

  ListDelta.reserve(p);
  ListGamma.reserve(p);
  ListDelta.resize(p, 0.0);
  ListGamma.resize(p, 0.0);
  ListAllDelta.resize(0);
  ListAllGamma.resize(0);

  MTest_a.resize(0);
  MTest_b.resize(0);

  if ((useUniquePoint)&&(!Use_MTest))
    {
      // Calcul de la liste des distances par rapport au point i

      i = IndexUP;

      cout << "GammaTest: search for " << p << " neighbors of point " << i + 1 << endl;
      
      vecAns = kd_tree->search_n(p, InputData[i]);
      
#ifdef DEBUG
      cout << "Initial point : ";
      for(j=0; j<InputData[i].size(); j++)
	{
	  cout << InputData[i][j] << " ";
	} /* End For */
      cout << OutputData[i][0] << endl;

      cout << "vecAns list :" << endl;
      for(j=0; j<vecAns.size(); j++)
	{
	  cout << "Point " << j << " : ";
	  for(k=0; k<vecAns[j].size(); k++)
	    {
	      cout << vecAns[j][k] << " ";
	    } /* End For */
	  cout << endl;
	} /* End For */
#endif

      // Calcul de Delta et Gamma
      for(j=Starting_n; j<vecAns.size(); j++)
	{
	  Select = ((rand()/(double)RAND_MAX)<Proportion);
	  
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
	  
	  if (Select)
	    {
	      ListAllDelta.push_back(AuxDeltaValue);
	      ListAllGamma.push_back(AuxGammaValue);
	    } /* End If */
	} /* End For */
    } /* End If */
  else if (!Use_MTest)
    {
      cout << endl;

      for(i=0; i<InputData.size(); i++)
	{
	  // Calcul de la liste des distances par rapport au point i
	  cout << UP;
	  cout << "GammaTest: search for " << p << " neighbors of point " << i + 1 << " / " << InputData.size() << endl;
	  
	  vecAns = kd_tree->search_n(p, InputData[i]);
	  
#ifdef DEBUG
	  cout << "Initial point : ";
	  for(j=0; j<InputData[i].size(); j++)
	    {
	      cout << InputData[i][j] << " ";
	    } /* End For */
	  cout << OutputData[i][0] << endl;

	  cout << "vecAns list :" << endl;
	  for(j=0; j<vecAns.size(); j++)
	    {
	      cout << "Point " << j << " : ";
	      for(k=0; k<vecAns[j].size(); k++)
		{
		  cout << vecAns[j][k] << " ";
		} /* End For */
	      cout << endl;
	    } /* End For */
#endif

	  // Calcul de Delta et Gamma
	  for(j=Starting_n; j<vecAns.size(); j++)
	    {
	      Select = ((rand()/(double)RAND_MAX)<Proportion);
	      
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
	      
	      if (Select)
		{
		  ListAllDelta.push_back(AuxDeltaValue);
		  ListAllGamma.push_back(0.5*AuxGammaValue);
		} /* End If */
	    } /* End For */
	} /* End For */
    } /* End Else */

  if (!Use_MTest)
    {
      cout << "GammaTest: normalisation of the lists" << endl;
      
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
	      
	      ListAuxDelta[i] /= (double)(i+1.0);
	      ListAuxGamma[i] /= (double)(i+1.0);
	    } /* End For */
	  
	  for(i=Starting_n; i<p; i++)
	    {
	      ListDelta[i] = ListAuxDelta[i];
	      ListGamma[i] = ListAuxGamma[i];
	    } /* End For */
	  
	  ListAuxDelta.resize(0);
	  ListAuxGamma.resize(0);
	} /* End If */
      
      // Calcul des coeff de regression entre Delta et Gamma: Gamma = A.Delta + B
      
      cout << "GammaTest: computation of the regression line" << endl;
      
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
      
//       b  = m_dxdy / m_dx2;
//       a  = m_y - m_x * b;
      a = (m_dx2*m_y - m_x*m_dxdy)/(double)Delta;
      b = (m*m_dxdy - m_x*m_y)/(double)Delta;
    } /* End If */

  if (Use_MTest)
    {
      ListDelta.reserve(MTest_Max);
      ListGamma.reserve(MTest_Max);
      ListDelta.resize(MTest_Max, 0.0);
      ListGamma.resize(MTest_Max, 0.0);

      cout << endl;

      for(i=0; i<InputData.size(); i++)
	{
	  // Calcul de la liste des distances par rapport au point i
	  
	  cout << UP;
	  cout << "GammaTest_MTest: searching for " << MTest_Max << " neighbors of point " << i << endl;
	  
	  vecAns = kd_tree->search_n(MTest_Max, InputData[i]);
	  
#ifdef DEBUG
	  cout << "Initial point : ";
	  for(j=0; j<InputData[i].size(); j++)
	    {
	      cout << InputData[i][j] << " ";
	    } /* End For */
	  cout << OutputData[i][0] << endl;
	  
	  cout << "vecAns list :" << endl;
	  for(j=0; j<vecAns.size(); j++)
	    {
	      cout << "Point " << j << " : ";
	      for(k=0; k<vecAns[j].size(); k++)
		{
		  cout << vecAns[j][k] << " ";
		} /* End For */
	      cout << endl;
	    } /* End For */
#endif
	      
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

      // Normalisation de Delta et Gamma
      for(i=Starting_n; i<(unsigned int)MTest_Max; i++)
	{
	  ListDelta[i] /= (double)(InputData.size());
	  ListGamma[i] /= (double)(2*InputData.size());
	} /* End For */
      
      ListMTestDelta.resize(ListDelta.size(), 0.0);
      ListMTestGamma.resize(ListGamma.size(), 0.0);

      for(l=(unsigned int)MTest_Min; l<(unsigned int)MTest_Max; l++)
	{
	  cout << "GammaTest_MTest: normalisation of the lists" << endl;
	  
	  for(i=Starting_n; i<l; i++)
	    {
	      ListMTestDelta[i] = ListDelta[i];
	      ListMTestGamma[i] = ListGamma[i];
	    } /* End For */
	  
	  if (useVersion2)
	    {
	      ListAuxDelta.resize(ListDelta.size(), 0.0);
	      ListAuxGamma.resize(ListGamma.size(), 0.0);
	      
	      for(i=Starting_n; i<l; i++)
		{
		  for(j=Starting_n; j<i; j++)
		    {
		      ListAuxDelta[i] += ListDelta[j];
		      ListAuxGamma[i] += ListGamma[j];
		    } /* End For */
		  
		  ListAuxDelta[i] /= (double)(i+1.0);
		  ListAuxGamma[i] /= (double)(i+1.0);
		} /* End For */
	      
	      for(i=Starting_n; i<l; i++)
		{
		  ListMTestDelta[i] = ListAuxDelta[i];
		  ListMTestGamma[i] = ListAuxGamma[i];
		} /* End For */
	      
	      ListAuxDelta.resize(0);
	      ListAuxGamma.resize(0);
	    } /* End If */
	  
	  // Calcul des coeff de regression entre Delta et Gamma: Gamma = A.Delta + B
	  
	  cout << "GammaTest_MTest: computation of the regression line" << endl;
	  
	  m_x    = 0.0;
	  m_y    = 0.0;
	  m_dx2  = 0.0;
	  m_dxdy = 0.0;
	  
	  for(i=Starting_n; i<l; i++)
	    {
	      m_x += ListMTestDelta[i];
	      m_y += ListMTestGamma[i];
	    } /* End For */
	  
	  m_x = m_x / (double)(l - Starting_n);
	  m_y = m_y / (double)(l - Starting_n);
	  
	  for(i=Starting_n; i<l; i++)
	    {
	      dx = ListMTestDelta[i] - m_x;
	      dy = ListMTestGamma[i] - m_y;
	      
	      m_dx2  += (dx * dx);
	      m_dxdy += (dx * dy);
	    } /* End For */
	  
	  m_dx2  = m_dx2 / (double)(l  - Starting_n);
	  m_dxdy = m_dxdy / (double)(l - Starting_n);
	  
	  /* In terms of y = a + b x */
	  
	  b  = m_dxdy / m_dx2;
	  a  = m_y - m_x * b;
	  
	  MTest_a.push_back(a);
	  MTest_b.push_back(b);
	} /* End For */
    } /* End If */

  if (kd_tree) delete kd_tree;
}
