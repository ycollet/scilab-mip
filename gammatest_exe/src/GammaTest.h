// GammaTest.h
// Author: Y. Collette
// Date: 19/11/2005
// Renault - DR - 64240

#ifndef GAMMATEST_H
#define GAMMATEST_H

#include <vector>

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
	       int BucketSize);
#endif
