// Author: Y. Collette
// EMail: ycollet@freesurf.fr
// Date: 02/04/2007

/*
This file is part of Lolimot.

    Foobar is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Foobar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <stdlib.h>
#include <math.h>

#include <vector>
#include <limits>
#include <iostream>
#include <limits>
#include <utility>
#include <algorithm>
#include <fstream>

#include "TrainLolimotStruct.h"

using namespace std;

#define UP "\033[A"

class Compare_Pair
{
public:
  int operator()(const pair<float, unsigned int> & P1,
		 const pair<float, unsigned int> & P2) const
  {
    return (P1.first > P2.first);
  }
};

void Filter_AddMinMax(vector<vector<float> > & InputData,
		      vector<vector<float> > & MeasureData,
		      vector<vector<float> > & Selected_InputData,
		      vector<vector<float> > & Selected_MeasureData)
{
  unsigned int i, j, k, IndexMin = 0, IndexMax = 0;
  float        Min, Max;
  bool         Same = true;

  for(i=0; i<InputData[0].size(); i++)
    {
      // We look for the min and the max for dimension i

      cout << "Filter_AddMinMax:: adding points corresponding to min and max value of dimension " << i << endl;

      Min = numeric_limits<float>::max();
      Max = numeric_limits<float>::min();

      for(j=0; j<InputData.size(); j++)
	{
	  if (Min>InputData[j][i])
	    {
	      Min = InputData[j][i];
	      IndexMin = j;
	    } /* End If */
	  if (Max<InputData[j][i])
	    {
	      Max = InputData[j][i];
	      IndexMax = j;
	    } /* End If */
	} /* End For */

      // We verify that these two points are not already in Selected_InputData
      // We verify the min point
      for(j=0; j<Selected_InputData.size(); j++)
	{
	  Same = true;

	  for(k=0; k<Selected_InputData[j].size(); k++)
	    {
	      Same = Same && (Selected_InputData[j][k]==InputData[IndexMin][k]);
	    } /* End For */

	  if (Same) break;
	} /* End For */

      if (!Same) 
	{
	  Selected_InputData.push_back(InputData[IndexMin]);
	  Selected_MeasureData.push_back(MeasureData[IndexMin]);
	} /* End If */

      // We verify the max point
      for(j=0; j<Selected_InputData.size(); j++)
	{
	  Same = true;

	  for(k=0; k<Selected_InputData[j].size(); k++)
	    {
	      Same = Same && (Selected_InputData[j][k]==InputData[IndexMax][k]);
	    } /* End For */

	  if (Same) break;
	} /* End For */

      if (!Same) 
	{
	  Selected_InputData.push_back(InputData[IndexMax]);
	  Selected_MeasureData.push_back(MeasureData[IndexMax]);
	} /* End If */
    } /* End For */


  for(i=0; i<MeasureData[0].size(); i++)
    {
      // We look for the min and the max for dimension i

      cout << "Filter_AddMinMax:: adding points corresponding to min and max value of measure " << i << endl;

      Min = numeric_limits<float>::max();
      Max = numeric_limits<float>::min();

      for(j=0; j<MeasureData.size(); j++)
	{
	  if (Min>MeasureData[j][i])
	    {
	      Min = MeasureData[j][i];
	      IndexMin = j;
	    } /* End If */
	  if (Max<MeasureData[j][i])
	    {
	      Max = MeasureData[j][i];
	      IndexMax = j;
	    } /* End If */
	} /* End For */

      // We verify that these two points are not already in Selected_MeasureData
      // Verification of the min point
      for(j=0; j<Selected_MeasureData.size(); j++)
	{
	  Same = true;

	  for(k=0; k<Selected_MeasureData[j].size(); k++)
	    {
	      Same = Same && (Selected_MeasureData[j][k]==MeasureData[IndexMin][k]);
	    } /* End For */

	  if (Same) break;
	} /* End For */

      if (!Same) 
	{
	  Selected_InputData.push_back(InputData[IndexMin]);
	  Selected_MeasureData.push_back(MeasureData[IndexMin]);
	} /* End If */

      // Verification of the max point
      for(j=0; j<Selected_MeasureData.size(); j++)
	{
	  Same = true;

	  for(k=0; k<Selected_MeasureData[j].size(); k++)
	    {
	      Same = Same && (Selected_MeasureData[j][k]==MeasureData[IndexMax][k]);
	    } /* End For */

	  if (Same) break;
	} /* End For */

      if (!Same) 
	{
	  Selected_InputData.push_back(InputData[IndexMax]);
	  Selected_MeasureData.push_back(MeasureData[IndexMax]);
	} /* End If */
    } /* End For */

  if (InputData.size()==0)
    {
      cerr << "Filter_AddMinMax: Error - Learning data set empty" << endl;
      exit(1);
    } /* End If */
}

void Filter_Filter(vector<vector<float> > & InputData,
		   vector<vector<float> > & MeasureData,
		   unsigned int Filter_Step)
{ 
  vector<vector<float> > Filtered_InputData;
  vector<vector<float> > Filtered_MeasureData;
  vector<float>          OneInput, OneMeasure;
  unsigned int            IndexStep;
  unsigned int            i, j, k;

  OneInput.resize(InputData[0].size());
  OneMeasure.resize(MeasureData[0].size());
  
  Filtered_MeasureData.resize(0);
  Filtered_InputData.resize(0);
  
  for(i=0; i<InputData.size(); i+=Filter_Step)
    {
      for(j=0; j<InputData[i].size(); j++)
	{
	  OneInput[j] = InputData[i][j];
	  
	  // We compute the mean on the Filter_Step samples of the input parameters.
	  if ((i+Filter_Step)<InputData.size())
	    {
	      IndexStep = Filter_Step;
	    } /* End Else */
	  else
	    {
	      IndexStep = InputData.size() - i;
	    } /* End Else */
	  
	  for(k=1; k<IndexStep; k++)
	    {
	      OneInput[j] += InputData[i+k][j];
	    } /* End For */
	  OneInput[j] = OneInput[j] / (float)IndexStep;
	} /* End For */
      Filtered_InputData.push_back(OneInput);
      
      for(j=0; j<MeasureData[i].size(); j++)
	{
	  OneMeasure[j] = MeasureData[i][j];
	  
	  // We compute the mean on the Filter_Step output samples.
	  if ((i+Filter_Step)<MeasureData.size())
	    {
	      IndexStep = Filter_Step;
	    } /* End Else */
	  else
	    {
	      IndexStep = MeasureData.size() - i;
	    } /* End Else */
	  
	  for(k=1; k<IndexStep; k++)
	    {
	      OneMeasure[j] += MeasureData[i+k][j];
	    } /* End For */
	  OneMeasure[j] = OneMeasure[j] / (float)IndexStep;
	} /* End For */
      Filtered_MeasureData.push_back(OneMeasure);
    } /* End For */

  InputData   = Filtered_InputData;
  MeasureData = Filtered_MeasureData;

  Filtered_InputData.clear();
  Filtered_MeasureData.clear();
  OneInput.clear();
  OneMeasure.clear();

  if (InputData.size()==0)
    {
      cerr << "Filter_Filter: Error - Learning data set empty" << endl;
      exit(1);
    } /* End If */
}

void Filter_MovingAverage(vector<vector<float> > & InputData,
			  vector<vector<float> > & MeasureData,
			  unsigned int MAverage)
{
  vector<vector<float> > Filtered_InputData;
  vector<vector<float> > Filtered_MeasureData;
  vector<float>          OneInput, OneMeasure;
  float                  MATemp = 0.0;
  unsigned int            IndexStep;
  unsigned int            i, j, k;

  OneInput.resize(InputData[0].size());
  OneMeasure.resize(MeasureData[0].size());
  
  Filtered_MeasureData.resize(MeasureData.size(), vector<float>(MeasureData[0].size(), 0.0));
  Filtered_InputData.resize(InputData.size(), vector<float>(InputData[0].size(), 0.0));

  for(i=0; i<Filtered_InputData.size(); i++)
    {
      for(j=0; j<Filtered_InputData[i].size(); j++)
	{
	  if ((i+MAverage)<Filtered_InputData.size())
	    {
	      IndexStep = MAverage;
	    } /* End Else */
	  else
	    {
	      IndexStep = Filtered_InputData.size() - i;
	    } /* End Else */
	  
	  MATemp = 0.0;
	  
	  for(k=0; k<IndexStep; k++)
	    {
	      MATemp += Filtered_InputData[i+k][j];
	    } /* End For */
	  
	  if (IndexStep!=0)
	    {
	      MATemp = MATemp / (float)IndexStep;
	    } /* End If */
	  Filtered_InputData[i][j] = MATemp;
	} /* End For */
      
      for(j=0; j<MeasureData[i].size(); j++)
	{
	  if ((i+MAverage)<Filtered_MeasureData.size())
	    {
	      IndexStep = MAverage;
	    } /* End Else */
	  else
	    {
	      IndexStep = MeasureData.size() - i;
	    } /* End Else */
	  
	  MATemp = 0.0;
	  
	  for(k=0; k<IndexStep; k++)
	    {
	      MATemp += MeasureData[i+k][j];
	    } /* End For */
	  
	  if (IndexStep!=0)
	    {
	      MATemp = MATemp / (float)IndexStep;
	    } /* End If */
	  Filtered_MeasureData[i][j] = MATemp;
	} /* End For */
    } /* End For */

  InputData   = Filtered_InputData;
  MeasureData = Filtered_MeasureData;

  Filtered_InputData.clear();
  Filtered_MeasureData.clear();
  OneInput.clear();
  OneMeasure.clear();

  if (InputData.size()==0)
    {
      cerr << "Filter_MovingAverage: Error - Learning data set empty" << endl;
      exit(1);
    } /* End If */
}

void Filter_Select(vector<vector<float> > & InputData,
		   vector<vector<float> > & MeasureData,
		   unsigned int Select_Step)
{ 
  vector<vector<float> > Selected_InputData;
  vector<vector<float> > Selected_MeasureData;
  vector<float>          OneInput, OneMeasure;
  unsigned int            i, j;

  OneInput.resize(InputData[0].size());
  OneMeasure.resize(MeasureData[0].size());
  
  Selected_MeasureData.resize(0);
  Selected_InputData.resize(0);
  
  for(i=0; i<InputData.size(); i+=Select_Step)
    {
      for(j=0; j<InputData[i].size(); j++)
	{
	  OneInput[j] = InputData[i][j];
	} /* End For j */
	  
      Selected_InputData.push_back(OneInput);
      
      for(j=0; j<MeasureData[i].size(); j++)
	{
	  OneMeasure[j] = MeasureData[i][j];
	} /* End For */

      Selected_MeasureData.push_back(OneMeasure);
    } /* End For */

  InputData   = Selected_InputData;
  MeasureData = Selected_MeasureData;

  Selected_InputData.clear();
  Selected_MeasureData.clear();
  OneInput.clear();
  OneMeasure.clear();

  if (InputData.size()==0)
    {
      cerr << "Filter_Select: Error - Learning data set empty" << endl;
      exit(1);
    } /* End If */
}

void Filter_EchantStat(vector<vector<float> > & InputData,
		       vector<vector<float> > & MeasureData,		       
		       unsigned int NbBar,
		       unsigned int NbPoints,
		       bool addMinMax)
{ 
  vector<vector<float> > Selected_InputData, TempList_InputData;
  vector<vector<float> > Selected_MeasureData, TempList_MeasureData;
  vector<unsigned int>   TempList_Index;
  vector<float>          Min, Max;
  vector<int>            Histogram, CumSum;
  unsigned int           NbPointsPerVar, Index, countNbPointsPerVar;
  double                 RandValue;
  vector<bool>           Selected;
  unsigned int           i, j;
  unsigned int           Index_Min, Index_Max;

  Selected_MeasureData.resize(0);
  Selected_InputData.resize(0);

  Min.resize(InputData[0].size(), numeric_limits<float>::max());
  Max.resize(InputData[0].size(), numeric_limits<float>::min());

  Selected.resize(InputData.size(), false);

  Histogram.resize(NbBar, 0);
  CumSum.resize(NbBar, 0);

  NbPointsPerVar = (int)(NbPoints / (float)InputData[0].size());

  cout << "Filter_EchantDyn:: préparation des listes min et max" << endl;

  for(i=0; i<InputData.size(); i++)
    {
      for(j=0; j<InputData[i].size(); j++)
	{
	  if (Min[j]>InputData[i][j]) Min[j] = InputData[i][j];
	  if (Max[j]<InputData[i][j]) Max[j] = InputData[i][j];
	} /* End For j */
    } /* End For i */

  cout << endl << endl << endl << endl;

  for(i=0; i<InputData[0].size(); i++)
    {
      // We build the histogram with respect to variable i
      cout << UP << UP << UP << UP;
      cout << "Filter_EchantStat:: construction of the histogram to select with respect to variable " << i << endl;
      
      for(j=0; j<Histogram.size(); j++)
	{
	  Histogram[j] = 0;
	} /* End For j */
      
      for(j=0; j<InputData.size(); j++)
	{
	  Index = (unsigned int)((InputData[j][i] - Min[i])/(Max[i] - Min[i]) * (Histogram.size()-1));
	  Histogram[Index]++;
	} /* End For j */
      
      // Computation of the cumulated sum
      cout << "Filter_EchantStat:: computation of the cumulated sum" << endl;
      
      CumSum[0] = Histogram[0];
      for(j=1; j<Histogram.size(); j++)
	{
	  CumSum[j] = CumSum[j-1] + Histogram[j];
	} /* End For j */
	        
      countNbPointsPerVar = 0;
      
      cout << endl;

      while(countNbPointsPerVar<NbPointsPerVar)
	{
	  cout << UP << "Filter_EchantStat:: sampling of the point " << countNbPointsPerVar << " out of " << NbPointsPerVar << endl;

	  // Sampling with respect to the probability distribution defined by histogram
	  RandValue = (CumSum[CumSum.size()-1] - CumSum[0])*(rand()/(double)RAND_MAX) + CumSum[0];

	  Index = 0;
	  while ((CumSum[Index]<=RandValue)&&(Index<CumSum.size())) Index++;

	  // We add an offset on the Index so as to stop being bothered by tests on the Index.
	  // It's less rigorous but easier to implement.
	  if (Index>=Histogram.size()-1) Index = Histogram.size()-2;
	  
	  // We get the value contained in the bar
	  TempList_InputData.resize(0);
	  TempList_MeasureData.resize(0);
	  TempList_Index.resize(0);

	  Index_Min = Index;
	  Index_Max = Index+1;

	  double Bound_0 = (Index_Min/(double)(Histogram.size()-1)*(Max[i] - Min[i]) + Min[i]);
	  double Bound_1 = (Index_Max/(double)(Histogram.size()-1)*(Max[i] - Min[i]) + Min[i]);

	  // If, after the sampling with respect to the distribution Histogram, the set TempList_Index
	  // is empty, we widen the bounds either by decreasing Index_Min, or by increasing Index_Max.
	  // If we reach 0 or Histogram.size()-1 as limits, then, we give up and go to the next variable
	  // We do the same with EchantDyn.

	  while((TempList_Index.size()==0)&&(Index_Min>=0)&&(Index_Max<=Histogram.size()-1))
	    {
	      for(j=0; j<InputData.size(); j++)
		{
		  if ((InputData[j][i]<=Bound_1)&&(InputData[j][i]>=Bound_0)&&(!Selected[j]))
		    {
		      TempList_InputData.push_back(InputData[j]);
		      TempList_MeasureData.push_back(MeasureData[j]);
		      TempList_Index.push_back(j);
		    } /* End If */
		} /* End For j */
	
	      if (rand()/(double)RAND_MAX>0.5) Index_Min--;
	      else                             Index_Max++;
	      
	      Bound_0 = (Index_Min/(double)(CumSum.size()-1)*(Max[i] - Min[i]) + Min[i]);
	      Bound_1 = (Index_Max/(double)(CumSum.size()-1)*(Max[i] - Min[i]) + Min[i]);	   
	    } /* End While */

	  // Random sampling of a value in this list
	  if (TempList_InputData.size()>0)
	    {
	      Index = (int)(rand()/(double)RAND_MAX * (TempList_InputData.size()-1));
	      Selected_InputData.push_back(TempList_InputData[Index]);
	      Selected_MeasureData.push_back(TempList_MeasureData[Index]);
	      Selected[TempList_Index[Index]] = true;
	      countNbPointsPerVar++;
	    } /* End If */
	} /* End While */

      cout << "Filter_EchantStat:: variable " << i << " done" << endl;
    } /* End For i */

  // We add the value corresponding to the min and max values for each dimension et and each measure
  if (addMinMax)
    {
      Filter_AddMinMax(InputData, MeasureData, Selected_InputData, Selected_MeasureData);
    } /* End If */
  
  InputData   = Selected_InputData;
  MeasureData = Selected_MeasureData;
  
  Selected_InputData.clear();
  Selected_MeasureData.clear();

  if (InputData.size()==0)
    {
      cerr << "Filter_EchantStat: Error - Learning data set empty" << endl;
      exit(1);
    } /* End If */
}

void Filter_EchantDyn(vector<vector<float> > & InputData,
		      vector<vector<float> > & MeasureData,
		      unsigned int NbBar,
		      unsigned int NbPoints,
		      bool addMinMax)
{ 
  vector<vector<float> > Selected_InputData,   TempList_InputData,  Dyn_InputData;
  vector<vector<float> > Selected_MeasureData, TempList_MeasureData;
  vector<unsigned int>    TempList_Index;
  vector<float>          Min, Max;
  vector<bool>            Selected;
  vector<int>             Histogram, CumSum;
  unsigned int            NbPointsPerVar, Index, countNbPointsPerVar;
  float                  RandValue;
  unsigned int            i, j;
  unsigned int            Index_Min, Index_Max;

  Selected_MeasureData.resize(0);
  Selected_InputData.resize(0);
  
  Selected.resize(InputData.size(), false);

  Histogram.resize(NbBar);
  CumSum.resize(NbBar);

  Max.resize(InputData[0].size(), numeric_limits<float>::min());
  Min.resize(InputData[0].size(), numeric_limits<float>::max());

  Dyn_InputData.resize(InputData.size());
  for(i=0; i<InputData.size(); i++)
    {
      Dyn_InputData[i].resize(InputData[0].size());
    } /* End For */

  NbPointsPerVar = (int)(NbPoints / (float)InputData[0].size());

  // Construction of the lists which contains the inputs dynamics
  for(j=0; j<InputData[0].size(); j++)
    {
      Dyn_InputData[0][j] = 0.0;
    } /* End For j */

  for(i=1; i<InputData.size(); i++)
    {
      for(j=0; j<InputData[i].size(); j++)
	{
	  Dyn_InputData[i][j] = InputData[i][j] - InputData[i-1][j];
	} /* End For j */
    } /* End For i */

  cout << "Filter_EchantDyn:: preparation of the min and max lists" << endl;

  for(i=0; i<InputData.size(); i++)
    {
      for(j=0; j<InputData[i].size(); j++)
	{
	  if (Min[j]>Dyn_InputData[i][j]) Min[j] = Dyn_InputData[i][j];
	  if (Max[j]<Dyn_InputData[i][j]) Max[j] = Dyn_InputData[i][j];
	} /* End For j */
    } /* End For i */

  cout << endl << endl << endl << endl;

  for(i=0; i<Dyn_InputData[0].size(); i++)
    {
      // Construction of the histogram with repect to variable i
      cout << UP << UP << UP << UP;
      cout << "Filter_EchantDyn:: preparation of the selectio histogram with repect to variable " << i << endl;
      
      for(j=0; j<Histogram.size(); j++)
	{
	  Histogram[j] = 0;
	} /* End For j */
      
      for(j=0; j<Dyn_InputData.size(); j++)
	{
	  Index = (unsigned int)((Dyn_InputData[j][i] - Min[i])/(Max[i] - Min[i]) * (Histogram.size()-1));
	  Histogram[Index]++;
	} /* End For j */
      
      cout << "Filter_EchantDyn:: preparation of the cumulated sum" << endl;
      
      CumSum[0] = Histogram[0];
      for(j=1; j<Histogram.size(); j++)
	{
	  CumSum[j] = CumSum[j-1] + Histogram[j];
	} /* End For j */
	        
      countNbPointsPerVar = 0;
      
      cout << endl;

      while(countNbPointsPerVar<NbPointsPerVar)
	{
	  cout << UP << "Filter_EchantDyn:: sampling of the point " << countNbPointsPerVar << " out of " << NbPointsPerVar << endl;

	  // Sampling with respect to the probability distribution defined by histogram
	  RandValue = (CumSum[CumSum.size()-1] - CumSum[0])*rand()/(double)RAND_MAX + CumSum[0];

	  Index = 0;
	  while ((CumSum[Index]<=RandValue)&&(Index<Histogram.size())) Index++;

	  // We shift Index so as to stop being bothered with the tests on Index.
	  // It's a less rigorous sampling, but it's easier to implement.
	  if (Index>=Histogram.size()-1) Index = Histogram.size()-2;
	  
	  // We get the valuse contained in this bar
	  TempList_InputData.resize(0);
	  TempList_MeasureData.resize(0);
	  TempList_Index.resize(0);

	  Index_Min = Index;
	  Index_Max = Index+1;

	  double Bound_0 = (Index_Min/(double)(Histogram.size()-1)*(Max[i] - Min[i]) + Min[i]);
	  double Bound_1 = (Index_Max/(double)(Histogram.size()-1)*(Max[i] - Min[i]) + Min[i]);

	  while((TempList_Index.size()==0)&&(Index_Min!=0)&&(Index_Max!=Histogram.size()-1))
	    {
	      for(j=0; j<Dyn_InputData.size(); j++)
		{
		  if ((Dyn_InputData[j][i]<Bound_1)&&(Dyn_InputData[j][i]>=Bound_0)&&(!Selected[j]))
		    {
		      TempList_InputData.push_back(InputData[j]);
		      TempList_MeasureData.push_back(MeasureData[j]);
		      TempList_Index.push_back(j);
		    } /* End If */
		} /* End For j */

	      if (rand()/(double)RAND_MAX>0.5) Index_Min--;
	      else                             Index_Max++;
	      
	      Bound_0 = (Index_Min/(double)(Histogram.size()-1)*(Max[i] - Min[i]) + Min[i]);
	      Bound_1 = (Index_Max/(double)(Histogram.size()-1)*(Max[i] - Min[i]) + Min[i]);
	    } /* End While */
	  
	  // Random sampling of a value in this list
	  if (TempList_InputData.size()>0)
	    {
	      Index = (unsigned int)(rand()/(double)RAND_MAX * (TempList_InputData.size()-1));
	      Selected_InputData.push_back(TempList_InputData[Index]);
	      Selected_MeasureData.push_back(TempList_MeasureData[Index]);
	      Selected[TempList_Index[Index]] = true;
	      countNbPointsPerVar++;
	    } /* End If */
	} /* End While */
      cout << "Filter_EchantDyn:: variable " << i << " done" << endl;
    } /* End For i */
  
  // We add the points which correspond to the min and max values for each dimension and each measure
  if (addMinMax)
    {
      Filter_AddMinMax(InputData, MeasureData, Selected_InputData, Selected_MeasureData);
    } /* End If */

  InputData   = Selected_InputData;
  MeasureData = Selected_MeasureData;
  
  Selected_InputData.clear();
  Selected_MeasureData.clear();
  Dyn_InputData.clear();

  if (InputData.size()==0)
    {
      cerr << "Filter_EchantDyn: Error - Learning data set empty" << endl;
      exit(1);
    } /* End If */
}

void Filter_EchantDelta(vector<vector<float> > & InputData,
			vector<vector<float> > & MeasureData,
			unsigned int NbPoints,
			unsigned int NbTimesToSelect,
			bool addMinMax,
			unsigned int Measure)
{ 
  vector<vector<float> > Selected_InputData, TempList_InputData;
  vector<vector<float> > Selected_MeasureData, TempList_MeasureData;
  vector<unsigned int>   TempList_Index;
  vector<float>          DeltaF_List, CumSum;
  unsigned int           Index, Index2, countNbPoints, NbPointsPerSelection;
  double                 RandValue;
  vector<bool>           Selected;
  unsigned int           i, j;

  NbPointsPerSelection = (unsigned int)(NbPoints/(double)NbTimesToSelect);

  Selected_MeasureData.resize(0);
  Selected_InputData.resize(0);

  Selected.resize(InputData.size(), false);

  CumSum.resize(InputData.size(), 0);
  DeltaF_List.resize(InputData.size(), 0.0);

  cout << endl << endl << endl << endl;

  for(i=0; i<NbTimesToSelect; i++)
    {
      Index = (unsigned int)((InputData.size()-1) * (rand()/(double)RAND_MAX));
      while(Selected[Index]) Index = (unsigned int)((InputData.size()-1) * (rand()/(double)RAND_MAX));

      cout << UP << UP << UP << UP;
      cout << "Filter_EchantDensity:: preparation of the histogram for selection with respect to the point " << Index << "    " << endl;

      for(j=0; j<InputData.size(); j++)
	{
	  DeltaF_List[j] = MeasureData[j][Measure] - MeasureData[Index][Measure];
	} /* End For */
      
      // Construction de la somme cumulée
      cout << "Filter_EchantDensity:: preparation of the cumulated sum" << endl;
      
      CumSum[0] = DeltaF_List[0];
      for(j=1; j<DeltaF_List.size(); j++)
	{
	  CumSum[j] = CumSum[j-1] + DeltaF_List[j];
	} /* End For j */
	        
      countNbPoints = 0;
      
      cout << endl;

      while(countNbPoints<NbPointsPerSelection)
	{
	  cout << UP << "Filter_EchantDensity:: sampling of point " << countNbPoints + 1 << " out of " << NbPointsPerSelection << "    " << endl;

	  // Sampling with respect to the probability distribution defined by histogram
	  RandValue = (CumSum[CumSum.size()-1] - CumSum[0])*(rand()/(double)RAND_MAX) + CumSum[0];

	  Index2 = 0;
	  while ((CumSum[Index2]<=RandValue)&&(Index2<CumSum.size())) Index2++;

	  // We shift Index so as to stop being bothered with the test on Index.
	  // A less rigorous sampling, but it's easier to implement.
	  if (Index2>=DeltaF_List.size()-1) Index2 = DeltaF_List.size() - 2;
	  
	  Selected_InputData.push_back(InputData[Index2]);
	  Selected_MeasureData.push_back(MeasureData[Index2]);
	  Selected[Index2] = true;
	  countNbPoints++;
	} /* End While */

      cout << "Filter_EchantDensity:: point " << Index << " done" << endl;
    } /* End For i */

  // We add the points which correspond to the min and max values of each dimension and each measure
  if (addMinMax)
    {
      Filter_AddMinMax(InputData, MeasureData, Selected_InputData, Selected_MeasureData);
    } /* End If */
  
  InputData   = Selected_InputData;
  MeasureData = Selected_MeasureData;
  
  Selected_InputData.clear();
  Selected_MeasureData.clear();

  if (InputData.size()==0)
    {
      cerr << "Filter_EchantDensity: Error - Learning data set empty" << endl;
      exit(1);
    } /* End If */
}

void Filter_EchantRand(vector<vector<float> > & InputData,
		       vector<vector<float> > & MeasureData,
		       unsigned int NbPoints,
		       bool addMinMax)
{ 
  vector<vector<float> > Selected_InputData;
  vector<vector<float> > Selected_MeasureData;
  vector<bool>            Selected;
  unsigned int            RandIndex, countNbPoints = 0;

  Selected_MeasureData.resize(0);
  Selected_InputData.resize(0);

  Selected.resize(InputData.size(), false);

  cout << endl;

  while(countNbPoints<NbPoints)
    {
      RandIndex = (unsigned int)(rand()/(double)RAND_MAX * InputData.size());

      while ((Selected[RandIndex])&&(RandIndex<=InputData.size()-1))
	RandIndex = (unsigned int)(rand()/(double)RAND_MAX * (InputData.size()-1));

      Selected_InputData.push_back(InputData[RandIndex]);
      Selected_MeasureData.push_back(MeasureData[RandIndex]);
      Selected[RandIndex] = true;
      countNbPoints++;
      cout << UP << "Filter_EchantRand:: point " << RandIndex << " selected - ";
      cout << countNbPoints << " / " << NbPoints << endl;
    } /* End While */
  cout << "Filter_EchantRand:: done" << endl;
  
  // We add the points which correspond to the min and max values of each dimension and each measure
  if (addMinMax)
    {
      Filter_AddMinMax(InputData, MeasureData, Selected_InputData, Selected_MeasureData);
    } /* End If */

  InputData   = Selected_InputData;
  MeasureData = Selected_MeasureData;
  
  Selected_InputData.clear();
  Selected_MeasureData.clear();
  Selected.clear();

  if (InputData.size()==0)
    {
      cerr << "Filter_EchantRand: Error - Learning data set empty" << endl;
      exit(1);
    } /* End If */
}

void Filter_EchantGreedy_AddMeasure(vector<vector<float> > & InputData,
				    vector<vector<float> > & MeasureData,
				    vector<float> & InputMin,
				    vector<float> & InputMax,
				    vector<float> & MeasureMin,
				    vector<float> & MeasureMax,
				    unsigned int NbPtsToSelect,
				    bool addMinMax,
				    unsigned int Measure)
{ 
  vector<vector<float> > Distance_Matrix;
  vector<bool>           Selected;
  vector<vector<float> > Selected_InputData;
  vector<vector<float> > Selected_MeasureData;
  float                  Distance, D_Aux;
  unsigned int           NbPtsToRemove, Index = 0, NbPtsRemoved = 0;
  bool                   ResultMin, ResultMax;
  unsigned int           i, j, k, l;

  Selected_MeasureData.resize(0);
  Selected_InputData.resize(0);
  
  Selected.resize(InputData.size(), false);

  Distance_Matrix.resize(InputData.size());
  for(i=0; i<InputData.size(); i++)
    {
      Distance_Matrix[i].resize(InputData.size(), 0.0);
    } /* End For */

  // Computation of the distance matrix

  cout << "Filter_EchantGreedy_AddMeasure:: Computation of the distance matrix" << endl;

  for(i=0; i<InputData.size()-1; i++)
    {
      for(j=i+1; j<InputData.size(); j++)
	{
	  Distance = 0.0;
	  for(k=0; k<InputData[0].size(); k++)
	    {
	      D_Aux = (InputData[i][k] - InputData[j][k]) / (InputMax[k] - InputMin[k]);
	      Distance += D_Aux * D_Aux;
	    } /* End For k */
	  D_Aux = (MeasureData[i][Measure] - MeasureData[j][Measure]) / (MeasureMax[Measure] - MeasureMin[Measure]);
	  Distance += D_Aux * D_Aux;

	  Distance_Matrix[i][j] = sqrt(Distance);
	} /* End For j */
    } /* End For i */

  // We remove the N closest points

  NbPtsToRemove = InputData.size() - NbPtsToSelect;

  cout << "Filter_EchantGreedy_AddMeasure:: Removing of " << NbPtsToRemove << endl << endl;

  for(i=0; i<NbPtsToRemove; i++)
    {
      Distance = numeric_limits<float>::max();
      for(j=0; j<InputData.size()-1; j++)
	{
	  if (Selected[j]) continue;
	  for(k=j+1; k<InputData.size(); k++)
	    {
	      if (Selected[k]) continue;

	      if (Distance>Distance_Matrix[j][k])
		{
		  ResultMin = false;
		  ResultMax = false;

		  // We look if a coordinate is equal to Min or Max.
		  // If true, we don't select this point because it corresponds to an
		  // extrem value

		  for(l=0; l<InputData[k].size(); l++)
		    {
		      ResultMin = ResultMin || (InputData[k][l]==InputMin[l]);
		      ResultMax = ResultMax || (InputData[k][l]==InputMax[l]);
		    } /* End For */

		  if (!ResultMin && !ResultMax)
		    {
		      Distance = Distance_Matrix[j][k];
		      Index    = k;
		    } /* End If */
		} /* End If */
	    } /* End For */
	} /* End For */
      Selected[Index] = true;
      NbPtsRemoved++;
      cout << UP << "Filter_EchantGreedy_AddMeasure:: Point " << Index << " removed - ";
      cout << NbPtsRemoved << " / " << NbPtsToRemove << "     " << endl;
    } /* End For */

  // We delete the distance matrix

  for(i=0; i<InputData.size(); i++)
    {
      Distance_Matrix[i].clear();
    } /* End For */
  Distance_Matrix.clear();

  // Creation of the list of selected points

  cout << "Filter_EchantGreedy_AddMeasure:: Creation of the list of selected points" << endl;

  for(i=0; i<InputData.size(); i++)
    {
      if (!Selected[i])
	{
	  Selected_MeasureData.push_back(MeasureData[i]);
	  Selected_InputData.push_back(InputData[i]);
	} /* End If */
    } /* End For */

  // We add the points which correspond to the mina dn max values for each dimension and each measure
  if (addMinMax)
    {
      Filter_AddMinMax(InputData, MeasureData, Selected_InputData, Selected_MeasureData);
    } /* End If */

  InputData   = Selected_InputData;
  MeasureData = Selected_MeasureData;
  
  Selected_InputData.clear();
  Selected_MeasureData.clear();
  Selected.clear();

  if (InputData.size()==0)
    {
      cerr << "Filter_EchantGreedy_AddMeasure:: Error - Learning data set empty" << endl;
      exit(1);
    } /* End If */
}

void Filter_EchantGreedy(vector<vector<float> > & InputData,
			 vector<vector<float> > & MeasureData,
			 vector<float> & InputMin,
			 vector<float> & InputMax,
			 unsigned int NbPtsToSelect,
			 bool addMinMax)
{ 
  vector<vector<float> > Distance_Matrix;
  vector<bool>            Selected;
  vector<vector<float> > Selected_InputData;
  vector<vector<float> > Selected_MeasureData;
  float                  Distance, D_Aux;
  unsigned int            NbPtsToRemove, Index = 0, NbPtsRemoved = 0;
  bool                    ResultMin, ResultMax;
  unsigned int            i, j, k, l;

  Selected_MeasureData.resize(0);
  Selected_InputData.resize(0);
  
  Selected.resize(InputData.size(), false);

  Distance_Matrix.resize(InputData.size());
  for(i=0; i<InputData.size(); i++)
    {
      Distance_Matrix[i].resize(InputData.size(), 0.0);
    } /* End For */

  // Computation of the distance matrix

  cout << "Filter_EchantGreedy:: Computation of the distance matrix" << endl;

  for(i=0; i<InputData.size()-1; i++)
    {
      for(j=i+1; j<InputData.size(); j++)
	{
	  Distance = 0.0;
	  for(k=0; k<InputData[0].size(); k++)
	    {
	      D_Aux = (InputData[i][k] - InputData[j][k]) / (InputMax[k] - InputMin[k]);
	      Distance += D_Aux * D_Aux;
	    } /* End For k */

	  Distance_Matrix[i][j] = sqrt(Distance);
	} /* End For j */
    } /* End For i */

  // Removing the N closest points

  NbPtsToRemove = InputData.size() - NbPtsToSelect;

  cout << "Filter_EchantGreedy:: Selection of " << NbPtsToRemove << " points to remove" << endl << endl;

  for(i=0; i<NbPtsToRemove; i++)
    {
      Distance = numeric_limits<float>::max();
      for(j=0; j<InputData.size()-1; j++)
	{
	  if (Selected[j]) continue;
	  for(k=j+1; k<InputData.size(); k++)
	    {
	      if (Selected[k]) continue;

	      if (Distance>Distance_Matrix[j][k])
		{
		  ResultMin = false;
		  ResultMax = false;

		  // We look if a coordinate is equal to Min or Max.
		  // If true, we don't select this point because it correspond to an
		  // extrem value

		  for(l=0; l<InputData[k].size(); l++)
		    {
		      ResultMin = ResultMin || (InputData[k][l]==InputMin[l]);
		      ResultMax = ResultMax || (InputData[k][l]==InputMax[l]);
		    } /* End For */

		  if (!ResultMin && !ResultMax)
		    {
		      Distance = Distance_Matrix[j][k];
		      Index    = k;
		    } /* End If */
		} /* End If */
	    } /* End For */
	} /* End For */
      Selected[Index] = true;
      NbPtsRemoved++;
      cout << UP << "Filter_EchantGreedy:: Point " << Index << " removed - ";
      cout << NbPtsRemoved << " / " << NbPtsToRemove << "     " << endl;
    } /* End For */

  // We delete the distance matrix

  for(i=0; i<InputData.size(); i++)
    {
      Distance_Matrix[i].clear();
    } /* End For */
  Distance_Matrix.clear();

  // Creation of the list of selected points

  cout << "Filter_EchantGreedy:: Creation of the list of selected points" << endl;

  for(i=0; i<InputData.size(); i++)
    {
      if (!Selected[i])
	{
	  Selected_MeasureData.push_back(MeasureData[i]);
	  Selected_InputData.push_back(InputData[i]);
	} /* End If */
    } /* End For */

  // We add the points which correspond to the min and max values for each dimension and each measure
  if (addMinMax)
    {
      Filter_AddMinMax(InputData, MeasureData, Selected_InputData, Selected_MeasureData);
    } /* End If */

  InputData   = Selected_InputData;
  MeasureData = Selected_MeasureData;
  
  Selected_InputData.clear();
  Selected_MeasureData.clear();
  Selected.clear();

  if (InputData.size()==0)
    {
      cerr << "Filter_EchantGreedy:: Error - Learning data set empty" << endl;
      exit(1);
    } /* End If */
}

void Filter_EchantGreedy_Fast_AddMeasure(vector<vector<float> > & InputData,
					 vector<vector<float> > & MeasureData,
					 vector<float> & InputMin,
					 vector<float> & InputMax,
					 vector<float> & MeasureMin,
					 vector<float> & MeasureMax,
					 unsigned int NbPtsToSelect,
					 unsigned int NbPtsAtTheEnd,
					 bool addMinMax,
					 unsigned int Measure)
{ 
  vector<float>          Distance_List;
  vector<bool>           Selected;
  vector<vector<float> > Selected_InputData;
  vector<vector<float> > Selected_MeasureData;
  float                  Distance, D_Aux;
  unsigned int           NbPtsToRemove, Index = 0, NbPtsRemoved = 0;
  unsigned int           i, j, k;

  Selected_MeasureData.resize(0);
  Selected_InputData.resize(0);
  
  Selected.resize(InputData.size(), false);

  Distance_List.resize(InputData.size());

  NbPtsToRemove = (InputData.size() - NbPtsAtTheEnd) / NbPtsToSelect;

  cout << endl << endl << endl << endl;

  for(i=0; i<NbPtsToSelect; i++)
    {
      Index = (unsigned int)(rand()/(double)RAND_MAX * (InputData.size() - 1));
      while (Selected[Index]) Index = (unsigned int)(rand()/(double)RAND_MAX * (InputData.size() - 1));

      cout << UP << UP << UP << UP;
      cout << "Filter_EchantGreedy_Fast_AddMeasure:: Selection of point " << Index << "       " << endl;

      // Computation of the distance list

      cout << "Filter_EchantGreedy_Fast_AddMeasure:: Computation of the list of distances" << endl;

      for(j=0; j<InputData.size(); j++)
	{
	  if (j==Index)
	    {
	      Distance_List[j] = numeric_limits<float>::max();
	      continue;
	    } /* End If */

	  Distance = 0.0;
	  for(k=0; k<InputData[0].size(); k++)
	    {
	      D_Aux = (InputData[Index][k] - InputData[j][k]) / (InputMax[k] - InputMin[k]);
	      Distance += D_Aux * D_Aux;
	    } /* End For k */
	  D_Aux = (MeasureData[Index][Measure] - MeasureData[j][Measure]) / (MeasureMax[Measure] - MeasureMin[Measure]);
	  Distance += D_Aux * D_Aux;

	  Distance_List[j] = sqrt(Distance);
	} /* End For j */

      // We remove the N closest points

      cout << "Filter_EchantGreedy_Fast_AddMeasure:: Selection of " << NbPtsToRemove << " points to be removed" << endl;

      NbPtsRemoved = 0;

      cout << endl;

      for(j=0; j<NbPtsToRemove; j++)
	{
	  Distance = numeric_limits<float>::max();

	  for(k=0; k<InputData.size(); k++)
	    {
	      if (Selected[k]) continue;

	      if (Distance>Distance_List[k])
		{
		  Distance = Distance_List[k];
		  Index    = k;
		} /* End If */
	    } /* End For */

	  Selected[Index] = true;
	  NbPtsRemoved++;
	  cout << UP << "Filter_EchantGreedy_Fast_AddMeasure:: Point " << Index << " removed - ";
	  cout << NbPtsRemoved << " / " << NbPtsToRemove << " - selected " << i << " / " << NbPtsToSelect << "  " << endl;
	} /* End For */
    } /* End For */

  // We remove the list of distances

  Distance_List.clear();

  // Creation of the list of selected points

  cout << "Filter_EchantGreedy_Fast_AddMeasure:: Creation of the list of selected points" << endl;

  for(i=0; i<InputData.size(); i++)
    {
      if (!Selected[i])
	{
	  Selected_MeasureData.push_back(MeasureData[i]);
	  Selected_InputData.push_back(InputData[i]);
	} /* End If */
    } /* End For */

  // We add the points which correspond to the mina dn max values for each dimension and each measure
  if (addMinMax)
    {
      Filter_AddMinMax(InputData, MeasureData, Selected_InputData, Selected_MeasureData);
    } /* End If */

  InputData   = Selected_InputData;
  MeasureData = Selected_MeasureData;
  
  Selected_InputData.clear();
  Selected_MeasureData.clear();
  Selected.clear();

  if (InputData.size()==0)
    {
      cerr << "Filter_EchantGreedy_Fast_AddMeasure:: Error - Learning data set empty" << endl;
      exit(1);
    } /* End If */
}

void Filter_EchantGreedy_Fast(vector<vector<float> > & InputData,
			      vector<vector<float> > & MeasureData,
			      vector<float> & InputMin,
			      vector<float> & InputMax,
			      unsigned int NbPtsToSelect,
			      unsigned int NbPtsAtTheEnd,
			      bool addMinMax)
{ 
  vector<float>          Distance_List;
  vector<bool>           Selected;
  vector<vector<float> > Selected_InputData;
  vector<vector<float> > Selected_MeasureData;
  float                  Distance, D_Aux;
  unsigned int           NbPtsToRemove, Index = 0, NbPtsRemoved = 0;
  unsigned int           i, j, k;

  Selected_MeasureData.resize(0);
  Selected_InputData.resize(0);
  
  Selected.resize(InputData.size(), false);

  Distance_List.resize(InputData.size());

  NbPtsToRemove = (InputData.size() - NbPtsAtTheEnd) / NbPtsToSelect;

  cout << endl << endl << endl << endl;

  for(i=0; i<NbPtsToSelect; i++)
    {
      Index = (unsigned int)(rand()/(double)RAND_MAX * (InputData.size() - 1));
      while (Selected[Index]) Index = (unsigned int)(rand()/(double)RAND_MAX * (InputData.size() - 1));

      cout << UP << UP << UP << UP;
      cout << "Filter_EchantGreedy_Fast:: Selection of point " << Index << "       " << endl;

      // Computation of the list of distances

      cout << "Filter_EchantGreedy_Fast:: Computation of the list of distances" << endl;

      for(j=0; j<InputData.size(); j++)
	{
	  if (j==Index)
	    {
	      Distance_List[j] = numeric_limits<float>::max();
	      continue;
	    } /* End If */

	  Distance = 0.0;
	  for(k=0; k<InputData[0].size(); k++)
	    {
	      D_Aux = (InputData[Index][k] - InputData[j][k]) / (InputMax[k] - InputMin[k]);
	      Distance += D_Aux * D_Aux;
	    } /* End For k */
	  Distance_List[j] = sqrt(Distance);
	} /* End For j */

      // We remove the N closest points

      cout << "Filter_EchantGreedy_Fast:: Selection of " << NbPtsToRemove << " points to be removed" << endl;

      NbPtsRemoved = 0;

      cout << endl;

      for(j=0; j<NbPtsToRemove; j++)
	{
	  Distance = numeric_limits<float>::max();

	  for(k=0; k<InputData.size(); k++)
	    {
	      if (Selected[k]) continue;

	      if (Distance>Distance_List[k])
		{
		  Distance = Distance_List[k];
		  Index    = k;
		} /* End If */
	    } /* End For */

	  Selected[Index] = true;
	  NbPtsRemoved++;
	  cout << UP << "Filter_EchantGreedy_Fast:: Point " << Index << " removed - ";
	  cout << NbPtsRemoved << " / " << NbPtsToRemove << " - selected " << i << " / " << NbPtsToSelect << "  " << endl;
	} /* End For */
    } /* End For */

  // We remove the list of distances

  Distance_List.clear();

  // Creation of the list of selected points

  cout << "Filter_EchantGreedy_Fast:: Creation of the list of selected points" << endl;

  for(i=0; i<InputData.size(); i++)
    {
      if (!Selected[i])
	{
	  Selected_MeasureData.push_back(MeasureData[i]);
	  Selected_InputData.push_back(InputData[i]);
	} /* End If */
    } /* End For */

  // We add the points which correspond to the min and max values of each dimension and each measure
  if (addMinMax)
    {
      Filter_AddMinMax(InputData, MeasureData, Selected_InputData, Selected_MeasureData);
    } /* End If */

  InputData   = Selected_InputData;
  MeasureData = Selected_MeasureData;
  
  Selected_InputData.clear();
  Selected_MeasureData.clear();
  Selected.clear();

  if (InputData.size()==0)
    {
      cerr << "Filter_EchantGreedy_Fast:: Error - Learning data set empty" << endl;
      exit(1);
    } /* End If */
}

void Filter_Echant_Output(vector<vector<float> > & InputData,
			  vector<vector<float> > & MeasureData,
			  vector<float> & MeasureMin,
			  vector<float> & MeasureMax,
			  unsigned int Measure,
			  unsigned int NbPtsToSelect,
			  bool UseMin,
			  bool addMinMax)
{
  vector<pair<float, unsigned int> > Ranked_Measure;
  vector<vector<float> >             Selected_MeasureData;
  vector<vector<float> >             Selected_InputData;
  unsigned int                        i;

  cout << "Filter_Echant_Output:: Creation of the ranked list" << endl;

  Ranked_Measure.resize(MeasureData.size());

  for(i=0; i<MeasureData.size(); i++)
    {
      Ranked_Measure[i].first  = MeasureData[i][Measure];
      Ranked_Measure[i].second = i;
    } /* End For */

  sort(Ranked_Measure.begin(), Ranked_Measure.end(), Compare_Pair());

  cout << "Filter_Echant_Output:: Selection of " << NbPtsToSelect << " points" << endl;

  Selected_InputData.resize(0);
  Selected_MeasureData.resize(0);

  if (UseMin)
    {
      for(i=MeasureData.size()-NbPtsToSelect; i<MeasureData.size(); i++)
	{
	  Selected_MeasureData.push_back(MeasureData[Ranked_Measure[i].second]);
	  Selected_InputData.push_back(InputData[Ranked_Measure[i].second]);
	} /* End For */
    } /* End If */
  else
    {
      for(i=0; i<NbPtsToSelect; i++)
	{
	  Selected_MeasureData.push_back(MeasureData[Ranked_Measure[i].second]);
	  Selected_InputData.push_back(InputData[Ranked_Measure[i].second]);
	} /* End For */
    } /* End Else */

  // We add the points which correspond to the min and max values of each dimension and each measure
  if (addMinMax)
    {
      Filter_AddMinMax(InputData, MeasureData, Selected_InputData, Selected_MeasureData);
    } /* End If */

  InputData   = Selected_InputData;
  MeasureData = Selected_MeasureData;
  
  Selected_InputData.clear();
  Selected_MeasureData.clear();

  if (InputData.size()==0)
    {
      cerr << "Filter_Echant_Output:: Error - Learning data set empty" << endl;
      exit(1);
    } /* End If */
}


void Filter_Echant_Input_Add(vector<vector<float> > & InputData,
			     vector<vector<float> > & MeasureData,
			     vector<float> & InputMin,
			     vector<float> & InputMax,
			     float Min_Bound,
			     float Max_Bound,
			     unsigned int Input_Index,
			     bool addMinMax)
{
  vector<pair<float, unsigned int> > Ranked_Input;
  vector<vector<float> >             Selected_MeasureData;
  vector<vector<float> >             Selected_InputData;
  unsigned int                       i;

  cout << "Filter_Echant_Input_Add:: Creation of the ranked list" << endl;

  Ranked_Input.resize(InputData.size());

  for(i=0; i<InputData.size(); i++)
    {
      Ranked_Input[i].first  = InputData[i][Input_Index];
      Ranked_Input[i].second = i;
    } /* End For */

  sort(Ranked_Input.begin(), Ranked_Input.end(), Compare_Pair());

  cout << "Filter_Echant_Input_Add:: Selection of points with respect to input " << Input_Index << endl;
  cout << "Interval: between " << Min_Bound << " and " << Max_Bound << endl;

  Selected_InputData.resize(0);
  Selected_MeasureData.resize(0);

  for(i=0; i<Ranked_Input.size(); i++)
    {
      if ((Ranked_Input[i].first<=Max_Bound)&&
	  (Ranked_Input[i].first>=Min_Bound))
	{
	  Selected_MeasureData.push_back(MeasureData[Ranked_Input[i].second]);
	  Selected_InputData.push_back(InputData[Ranked_Input[i].second]);
	} /* End If */
    } /* End For */

  // We add the points which correspond to the min and max values for each dimension and each measure
  if (addMinMax)
    {
      Filter_AddMinMax(InputData, MeasureData, Selected_InputData, Selected_MeasureData);
    } /* End If */

  InputData   = Selected_InputData;
  MeasureData = Selected_MeasureData;
  
  Selected_InputData.clear();
  Selected_MeasureData.clear();

  if (InputData.size()==0)
    {
      cerr << "Filter_Echant_Input_Add:: Error - Learning data set empty" << endl;
      exit(1);
    } /* End If */
}

void Filter_Echant_Input_Remove(vector<vector<float> > & InputData,
				vector<vector<float> > & MeasureData,
				vector<float> & InputMin,
				vector<float> & InputMax,
				float Min_Bound,
				float Max_Bound,
				unsigned int Input_Index,
				bool addMinMax)
{
  vector<pair<float, unsigned int> > Ranked_Input;
  vector<vector<float> >             Selected_MeasureData;
  vector<vector<float> >             Selected_InputData;
  unsigned int                       i;

  cout << "Filter_Echant_Input_Remove:: Creation of the ranked list" << endl;

  Ranked_Input.resize(InputData.size());

  for(i=0; i<InputData.size(); i++)
    {
      Ranked_Input[i].first  = InputData[i][Input_Index];
      Ranked_Input[i].second = i;
    } /* End For */

  sort(Ranked_Input.begin(), Ranked_Input.end(), Compare_Pair());

  cout << "Filter_Echant_Input_Remove:: Selection of points with respect to input " << Input_Index << endl;
  cout << "Interval: out of " << Min_Bound << " - " << Max_Bound << endl;

  Selected_InputData.resize(0);
  Selected_MeasureData.resize(0);

  for(i=0; i<Ranked_Input.size(); i++)
    {
      if ((Ranked_Input[i].first>=Max_Bound)||
	  (Ranked_Input[i].first<=Min_Bound))
	{
	  Selected_MeasureData.push_back(MeasureData[Ranked_Input[i].second]);
	  Selected_InputData.push_back(InputData[Ranked_Input[i].second]);
	} /* End If */
    } /* End For */

  // We add the points which correspond to the min and max values for each dimension and each measure
  if (addMinMax)
    {
      Filter_AddMinMax(InputData, MeasureData, Selected_InputData, Selected_MeasureData);
    } /* End If */

  InputData   = Selected_InputData;
  MeasureData = Selected_MeasureData;
  
  Selected_InputData.clear();
  Selected_MeasureData.clear();

  if (InputData.size()==0)
    {
      cerr << "Filter_Echant_Input_Remove:: Error - Learning data set empty" << endl;
      exit(1);
    } /* End If */
}

void Filter_Normalize(vector<vector<float> > & InputData,
		      vector<vector<float> > & MeasureData,
		      vector<float> & InputMin,
		      vector<float> & InputMax,
		      vector<float> & MeasureMin,
		      vector<float> & MeasureMax,
		      string Filename,
		      bool   SaveData)
{
  unsigned int i, j;
  ofstream     OutFile;

  // Normalisation of the inputs

  for(i=0; i<InputData.size(); i++)
    {
      for(j=0; j<InputData[i].size(); j++)
	{
	  InputData[i][j] = (InputData[i][j] - InputMin[j]) / 
	    (InputMax[j] - InputMin[j]);
	} /* End For */
    } /* End For */

  // Normalisation of the outputs

  for(i=0; i<MeasureData.size(); i++)
    {
      for(j=0; j<MeasureData[i].size(); j++)
	{
	  MeasureData[i][j] = (MeasureData[i][j] - MeasureMin[j]) /
	    (MeasureMax[j] - MeasureMin[j]);
	} /* End For */
    } /* End For */

  if (SaveData)
    {
      OutFile.open(Filename.c_str());
      
      // We export the old bounds
      
      for(i=0; i<InputMin.size(); i++)
	{
	  OutFile << "INPUT MIN " << InputMin[i] << " MAX " << InputMax[i] << endl;
	} /* End For */
      
      for(i=0; i<MeasureMin.size(); i++)
	{
	  OutFile << "OUPUT MIN " << MeasureMin[i] << " MAX " << MeasureMax[i] << endl;
	} /* End For */
      
      OutFile << endl << endl;
      
      // We export the new bounds
      
      for(i=0; i<InputMin.size(); i++)
	{
	  OutFile << "INPUT MIN " << 0 << " MAX " << 1 << endl;
	} /* End For */
      
      for(i=0; i<MeasureMin.size(); i++)
	{
	  OutFile << "OUPUT MIN " << 0 << " MAX " << 1 << endl;
	} /* End For */
      
      OutFile.close();
    } /* End If */
}

void Filter_EchantPartialGreedy_AddMeasure(vector<vector<float> > & InputData,
					   vector<vector<float> > & MeasureData,
					   vector<float> & InputMin,
					   vector<float> & InputMax,
					   vector<float> & MeasureMin,
					   vector<float> & MeasureMax,
					   unsigned int NbPtsToSelect,
					   unsigned int NbTimesToSelect,
					   bool addMinMax,
					   unsigned int Measure)
{ 
  vector<vector<float> > Distance_Matrix;
  vector<bool>           Selected;
  vector<unsigned int>   Index_Selected, Index_NotSelected, AuxIndex_NotSelected;
  vector<vector<float> > Selected_InputData, Aux_InputData;
  vector<vector<float> > Selected_MeasureData, Aux_MeasureData;
  float                  Distance, D_Aux;
  unsigned int           NbPtsToRemove, Index = 0, NbPtsRemoved = 0, SizeOfSet = 0;
  unsigned int           i, j, k, l;

  NbPtsToRemove = (unsigned int)((InputData.size() - NbPtsToSelect)/(double)NbTimesToSelect);
  NbPtsToSelect = (unsigned int)(NbPtsToSelect/(double)NbTimesToSelect);
  SizeOfSet     = (unsigned int)(InputData.size()/(double)NbTimesToSelect);

  Selected_MeasureData.resize(0);
  Selected_InputData.resize(0);
  
  Index_Selected.resize(0);
  Index_NotSelected.resize(0);
  for(i=0; i<InputData.size(); i++)
    {
      Index_NotSelected.push_back(i);
    } /* End For */
  Selected.resize(SizeOfSet, false);

  Distance_Matrix.resize(SizeOfSet);
  for(i=0; i<SizeOfSet; i++)
    {
      Distance_Matrix[i].resize(SizeOfSet, 0.0);
    } /* End For */

  // Sampling
  for(i=0; i<NbTimesToSelect; i++)
    {
      cout << "Filter_EchantPartialGreedy:: selection of point " << i + 1 << " / " << NbTimesToSelect << endl;
      // selection of SizeOfSet points which haven't been selected yet
      Selected.clear();
      Aux_InputData.clear();
      Aux_MeasureData.clear();

      Selected.resize(Index_NotSelected.size(), false);
      Aux_InputData.resize(0);
      Aux_MeasureData.resize(0);

      for(j=0; j<SizeOfSet; j++)
	{
	  Index = (unsigned int)(rand()/(double)RAND_MAX * (Index_NotSelected.size() - 1));
	  while (Selected[Index]==true)
	    {
	      Index = (unsigned int)(rand()/(double)RAND_MAX * (Index_NotSelected.size() - 1));
	    } /* End While */

	  Selected[Index] = true;

	  Aux_InputData.push_back(InputData[Index_NotSelected[Index]]);
	  Aux_MeasureData.push_back(MeasureData[Index_NotSelected[Index]]);
	} /* End For */

      // Updating the differences list
      AuxIndex_NotSelected.clear();
      AuxIndex_NotSelected.resize(0);

      for(j=0; j<Selected.size(); j++)
	{
	  if (Selected[j])
	    {
	      Index_Selected.push_back(Index_NotSelected[j]);
	    } /* End If */
	  else
	    {
	      AuxIndex_NotSelected.push_back(Index_NotSelected[j]);
	    } /* End Else */
	} /* End For */

      Index_NotSelected = AuxIndex_NotSelected;

      // Computation of the distance matrix
      cout << "Filter_EchantPartialGreedy_AddMeasure:: Computation of the distance matrix" << endl;

      for(j=0; j<SizeOfSet; j++)
	{
	  for(k=j+1; k<SizeOfSet; k++)
	    {
	      Distance = 0.0;
	      for(l=0; l<InputData[0].size(); l++)
		{
		  D_Aux = (Aux_InputData[j][l] - Aux_InputData[k][l]) / (InputMax[l] - InputMin[l]);
		  Distance += D_Aux * D_Aux;
		} /* End For l */
	      D_Aux = (Aux_MeasureData[j][Measure] - Aux_MeasureData[k][Measure]) / (MeasureMax[Measure] - MeasureMin[Measure]);
	      Distance += D_Aux * D_Aux;
	      
	      Distance_Matrix[j][k] = sqrt(Distance);
	      Distance_Matrix[k][j] = sqrt(Distance);
	    } /* End For k */
	} /* End For j */

      // We remove the N closest points

      Selected.clear();
      Selected.resize(SizeOfSet, false);

      NbPtsRemoved = 0;

      cout << "Filter_EchantPartialGreedy_AddMeasure:: Selection of " << NbPtsToRemove << " points to be removed" << endl << endl;
      
      for(j=0; j<NbPtsToRemove; j++)
	{
	  Distance = numeric_limits<float>::max();
	  for(k=0; k<SizeOfSet; k++)
	    {
	      if (Selected[k]) continue;
	      for(l=k+1; l<SizeOfSet; l++)
		{
		  if (Selected[l]) continue;
		  
		  if (Distance>Distance_Matrix[k][l])
		    {
		      Distance = Distance_Matrix[k][l];
		      Index    = l;
		    } /* End If */
		} /* End For */
	    } /* End For */
	  Selected[Index] = true;
	  NbPtsRemoved++;
	  cout << UP << "Filter_EchantPartialGreedy_AddMeasure:: Point " << Index << " removed - ";
	  cout << NbPtsRemoved << " / " << NbPtsToRemove << "     " << endl;
	} /* End For */

      // We add the selected points
      for(j=0; j<SizeOfSet; j++)
	{
	  if (!Selected[j])
	    {
	      Selected_InputData.push_back(Aux_InputData[j]);
	      Selected_MeasureData.push_back(Aux_MeasureData[j]);
	    } /* End If */
	} /* End For */
    } /* End For */

  // We remove the distance matrix
  for(i=0; i<SizeOfSet; i++)
    {
      Distance_Matrix[i].clear();
    } /* End For */
  Distance_Matrix.clear();


  // We add the points which correspond to the min and max values of each dimension and each measure
  if (addMinMax)
    {
      Filter_AddMinMax(InputData, MeasureData, Selected_InputData, Selected_MeasureData);
    } /* End If */

  InputData   = Selected_InputData;
  MeasureData = Selected_MeasureData;
  
  Selected_InputData.clear();
  Selected_MeasureData.clear();
  Selected.clear();

  if (InputData.size()==0)
    {
      cerr << "Filter_EchantPartialGreedy_AddMeasure:: Error - Learning data set empty" << endl;
      exit(1);
    } /* End If */
}

void Filter_EchantPartialGreedyDelta(vector<vector<float> > & InputData,
				     vector<vector<float> > & MeasureData,
				     vector<float> & InputMin,
				     vector<float> & InputMax,
				     vector<float> & MeasureMin,
				     vector<float> & MeasureMax,
				     unsigned int NbPtsToSelect,
				     unsigned int NbTimesToSelect,
				     bool addMinMax,
				     unsigned int Measure)
{ 
  vector<vector<float> > Distance_Matrix;
  vector<bool>           Selected;
  vector<unsigned int>   Index_Selected, Index_NotSelected, AuxIndex_NotSelected;
  vector<vector<float> > Selected_InputData, Aux_InputData;
  vector<vector<float> > Selected_MeasureData, Aux_MeasureData;
  float                  Distance;
  unsigned int           NbPtsToRemove, Index = 0, NbPtsRemoved = 0, SizeOfSet = 0;
  unsigned int           i, j, k, l;

  // We need that NbTimesToSelect <= NbPtsToSelect (we need enough points in SizeOfSet to
  // sample the right proportion of points in the set
  NbPtsToRemove = (unsigned int)((InputData.size() - NbPtsToSelect)/(double)NbTimesToSelect);
  SizeOfSet     = (unsigned int)(InputData.size()/(double)NbTimesToSelect);

  Selected_MeasureData.resize(0);
  Selected_InputData.resize(0);
  
  Index_Selected.resize(0);
  Index_NotSelected.resize(0);
  for(i=0; i<InputData.size(); i++)
    {
      Index_NotSelected.push_back(i);
    } /* End For */
  Selected.resize(SizeOfSet, false);

  Distance_Matrix.resize(SizeOfSet);
  for(i=0; i<SizeOfSet; i++)
    {
      Distance_Matrix[i].resize(SizeOfSet, 0.0);
    } /* End For */

  // Sampling
  for(i=0; i<NbTimesToSelect; i++)
    {
      cout << "Filter_EchantPartialGreedyDelta:: selection of point " << i + 1 << " / " << NbTimesToSelect << endl;
      // selection of SizeOfSet points which haven't been selected yet
      Selected.clear();
      Aux_InputData.clear();
      Aux_MeasureData.clear();

      Selected.resize(Index_NotSelected.size(), false);
      Aux_InputData.resize(0);
      Aux_MeasureData.resize(0);

      for(j=0; j<SizeOfSet; j++)
	{
	  Index = (unsigned int)(rand()/(double)RAND_MAX * (Index_NotSelected.size() - 1));
	  while (Selected[Index]==true)
	    {
	      Index = (unsigned int)(rand()/(double)RAND_MAX * (Index_NotSelected.size() - 1));
	    } /* End While */

	  Selected[Index] = true;

	  Aux_InputData.push_back(InputData[Index_NotSelected[Index]]);
	  Aux_MeasureData.push_back(MeasureData[Index_NotSelected[Index]]);
	} /* End For */

      // Updating of the differences list
      AuxIndex_NotSelected.clear();
      AuxIndex_NotSelected.resize(0);

      for(j=0; j<Selected.size(); j++)
	{
	  if (Selected[j])
	    {
	      Index_Selected.push_back(Index_NotSelected[j]);
	    } /* End If */
	  else
	    {
	      AuxIndex_NotSelected.push_back(Index_NotSelected[j]);
	    } /* End Else */
	} /* End For */

      Index_NotSelected = AuxIndex_NotSelected;

      // Computation of the distance matrix
      cout << "Filter_EchantPartialGreedyDensity:: Computation of the distance matrix" << endl;

      for(j=0; j<SizeOfSet; j++)
	{
	  Distance_Matrix[j][j] = numeric_limits<float>::max();

	  for(k=j+1; k<SizeOfSet; k++)
	    {
	      Distance_Matrix[j][k] = Aux_MeasureData[j][Measure] - Aux_MeasureData[k][Measure];
	      Distance_Matrix[k][j] = Aux_MeasureData[j][Measure] - Aux_MeasureData[k][Measure];
	    } /* End For k */
	} /* End For j */

      // We remove the N closest points

      Selected.clear();
      Selected.resize(SizeOfSet, false);

      NbPtsRemoved = 0;

      cout << "Filter_EchantPartialGreedyDensity:: Selection of " << NbPtsToRemove << " points to be removed" << endl << endl;
      
      for(j=0; j<NbPtsToRemove; j++)
	{
	  Distance = numeric_limits<float>::max();
	  for(k=0; k<SizeOfSet; k++)
	    {
	      if (Selected[k]) continue;
	      for(l=k+1; l<SizeOfSet; l++)
		{
		  if (Selected[l]) continue;
		  
		  if (Distance>Distance_Matrix[k][l])
		    {
		      Distance = Distance_Matrix[k][l];
		      Index    = l;
		    } /* End If */
		} /* End For */
	    } /* End For */
	  Selected[Index] = true;
	  NbPtsRemoved++;
	  cout << UP << "Filter_EchantPartialGreedyDensity:: Point " << Index << " removed - ";
	  cout << NbPtsRemoved << " / " << NbPtsToRemove << "     " << endl;
	} /* End For */

      // We add the selected points
      for(j=0; j<SizeOfSet; j++)
	{
	  if (!Selected[j])
	    {
	      Selected_InputData.push_back(Aux_InputData[j]);
	      Selected_MeasureData.push_back(Aux_MeasureData[j]);
	    } /* End If */
	} /* End For */
    } /* End For */

  // We remove the distance matrix
  for(i=0; i<SizeOfSet; i++)
    {
      Distance_Matrix[i].clear();
    } /* End For */
  Distance_Matrix.clear();


  // We add the points which correspond to the min and max values for each dimension and each measure
  if (addMinMax)
    {
      Filter_AddMinMax(InputData, MeasureData, Selected_InputData, Selected_MeasureData);
    } /* End If */

  InputData   = Selected_InputData;
  MeasureData = Selected_MeasureData;
  
  Selected_InputData.clear();
  Selected_MeasureData.clear();
  Selected.clear();

  if (InputData.size()==0)
    {
      cerr << "Filter_EchantPartialGreedyDelta:: Error - Learning data set empty" << endl;
      exit(1);
    } /* End If */
}

void Skip(vector<vector<float> > & InputData,
	  vector<vector<float> > & MeasureData,
	  vector<Data_Skip> & List_Skip,
	  unsigned int NumOfFile)
{
  vector<vector<float> >::iterator itInpBegin, itInpEnd, itMeasBegin, itMeasEnd;

  cout << "Skip" << endl;

  for(unsigned int i=0; i<List_Skip.size(); i++)
    {
      if (List_Skip[i].NumOfFile==(int)NumOfFile)
	{
	  cout << "Skipping sequence from index " << List_Skip[i].Begin_Part;
	  cout << " to index " << List_Skip[i].End_Part << endl;
	  
	  itInpBegin = InputData.begin();
	  itInpEnd   = InputData.begin();
	  advance(itInpBegin, List_Skip[i].Begin_Part);
	  advance(itInpEnd,   List_Skip[i].End_Part + 1);
	  InputData.erase(itInpBegin, itInpEnd);
	  
	  itMeasBegin = MeasureData.begin();
	  itMeasEnd   = MeasureData.begin();
	  advance(itMeasBegin, List_Skip[i].Begin_Part);
	  advance(itMeasEnd,   List_Skip[i].End_Part + 1);
	  MeasureData.erase(itMeasBegin, itMeasEnd);
	} /* End If */
    } /* End For */
}
