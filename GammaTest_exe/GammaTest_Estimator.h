// GammaTest.h
// Author: Y. Collette
// Date: 19/11/2005
// Renault - DR - 64240

#ifndef GAMMATEST_ESTIMATOR_H
#define GAMMATEST_ESTIMATOR_H

#include <vector>

using namespace std;

void GammaTest_Estimator(vector<vector<double> > & InputData,
			 vector<vector<double> > & OutputData,
			 vector<double> & y_Output,
			 unsigned int p,
			 unsigned int Measure,
			 bool useVersion2,
			 unsigned int Starting_n,
			 int BucketSize);
#endif
