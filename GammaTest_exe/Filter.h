// Filter.h
// Author: Y. Collette
// Date: 19/11/2005
// Renault - DR - 64240

#ifndef FILTER_H
#define FILTER_H

#include <vector>

using namespace std;

void Filter_AddMinMax(vector<vector<double> > & InputData,
		      vector<vector<double> > & MeasureData,
		      vector<vector<double> > & Selected_InputData,
		      vector<vector<double> > & Selected_MeasureData);

void Filter_Select(vector<vector<double> > & InputData,
		   vector<vector<double> > & MeasureData,
		   unsigned int Select_Step);

void Filter_EchantStat(vector<vector<double> > & InputData,
		       vector<vector<double> > & MeasureData,
		       unsigned int NbBar,
		       unsigned int NbPoints);

void Filter_EchantRand(vector<vector<double> > & InputData,
		       vector<vector<double> > & MeasureData,
		       unsigned int NbPoints);

void Filter_EchantGreedy(vector<vector<double> > & InputData,
			 vector<vector<double> > & MeasureData,
			 vector<double> & InputMin,
			 vector<double> & InputMax,
			 unsigned int NbPtsToSelect);

void Filter_EchantGreedy_Fast(vector<vector<double> > & InputData,
			      vector<vector<double> > & MeasureData,
			      vector<double> & InputMin,
			      vector<double> & InputMax,
			      unsigned int NbPtsToSelect,
			      unsigned int NbPtsAtTheEnd);

void Filter_Echant_Output(vector<vector<double> > & InputData,
			  vector<vector<double> > & MeasureData,
			  vector<double> & MeasureMin,
			  vector<double> & MeasureMax,
			  unsigned int Measure,
			  unsigned int NbPtsToSelect,
			  bool UseMin);
#endif
