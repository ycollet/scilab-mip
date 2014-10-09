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

#include <limits>
#include <vector>
#include <fstream>
#include <iostream>

#include "LL_Lolimot.h"
#include "LL_Mesure.h"
#include "LL_Partition.h"

using namespace std;

#define UP "\033[A"

//#define DEBUG 1

enum Model_Input_Type {Input = 0, Output = 1, RetroInput = 2, Weight = 3};

struct Model_Input
{
  string       Name;
  unsigned int Pos;
  unsigned int Tempo;
  unsigned int FileNumber;
  bool         DontCut;
  unsigned int GroupNb;
  Model_Input_Type Type;
};

void Init_Temporize(vector<vector<float> > & Temporize, vector<unsigned int> & Index_Tempo,
		    vector<Model_Input> & List_Model_Input)
{
  vector<float> X_Test;
  unsigned int i=0;

  Temporize.resize(0);
  Index_Tempo.resize(0);

  // Creation of the temporisation for the retro inputs
  for(i=0; i<List_Model_Input.size(); i++)
    {
      if (List_Model_Input[i].Type==RetroInput)
	{
	  X_Test.resize(List_Model_Input[i].Tempo, 0.0);
	  Temporize.push_back(X_Test);
	  Index_Tempo.push_back(i);
	} /* End If */
    } /* End For */
}

void Update_Temporize(vector<vector<float> > & Temporize, float Y_Res)
{
  unsigned int j=0, k=0;

  // Storage of this measure in temporize
  for(j=0; j<Temporize.size(); j++)
    {
      for(k=Temporize[j].size()-1; k>0; k--)
	{
	  Temporize[j][k] = Temporize[j][k-1];
	} /* End For */
      Temporize[j][0] = Y_Res;
    } /* End For */
}

void Learning_NoMeth(LL_Lolimot * Lolimot,
		     vector<vector<float> > & InputData,
		     vector<vector<float> > & MeasureData,
		     vector<float> & WeightData,
		     vector<Model_Input> & List_Model_Input,
		     vector<float> & InputMin,
		     vector<float> & InputMax,
		     vector<float> & MeasureMin,
		     vector<float> & MeasureMax,
		     float Sigma,
		     int NbCut,
		     int Measure,
		     bool Loaded_Lolimot)
{
  LL_Mesure  * Mesure = NULL;
  unsigned int i, j;

  cout << "Learning NoMeth: No learning of the LOLIMOT network" << endl;

  if (!Loaded_Lolimot)
    {
      for(i=0; i<List_Model_Input.size(); i++)
	{
	  Lolimot->addDimension(List_Model_Input[i].Name.c_str(), InputMin[i], InputMax[i], Sigma, NbCut);
	} /* End For */
      
      for(i=0; i<List_Model_Input.size(); i++)
	{
	  if (List_Model_Input[i].DontCut) Lolimot->setDontCutVar(i, true);
	} /* End For */
    } /* End If */
      
  for(i=0; i<MeasureData.size(); i++)
    {
      Mesure = Lolimot->addMesure(MeasureData[i][Measure],  MeasureMin[Measure],  MeasureMax[Measure]);
      
      for(j=0; j<List_Model_Input.size(); j++)
	{
	  Mesure->setDimensionValue(j, InputData[i][j]);
	} /* End For */
    } /* End For */

  if (WeightData.size()!=0)
    {
      for(i=0; i<MeasureData.size(); i++)
	{
	  Lolimot->addWeight(WeightData[i]);
	} /* End For */
    } /* End If */
  else
    {
      Lolimot->useDefaultWeight();
    } /* End Else */

  if (Loaded_Lolimot)
    {
      float          residu;
      LL_Partition * AuxPart = NULL;

      Lolimot->updateCoefficients(residu, AuxPart);
      cout << "DEBUG: residu = " << residu << endl;
    } /* End If */
}

void Learning_Classic(LL_Lolimot * Lolimot,
		      vector<vector<float> > & InputData,
		      vector<vector<float> > & MeasureData,
		      vector<float> & WeightData,
		      vector<Model_Input> & List_Model_Input,
		      vector<float> & InputMin,
		      vector<float> & InputMax,
		      vector<float> & MeasureMin,
		      vector<float> & MeasureMax,
		      float Sigma,
		      int NbCut,
		      int Measure,
		      bool Loaded_Lolimot)
{
  LL_Mesure  * Mesure = NULL;
  unsigned int i, j;

  cout << "Learning Classic: Learning of the LOLIMOT network" << endl;

  if (!Loaded_Lolimot)
    {
      for(i=0; i<List_Model_Input.size(); i++)
	{
	  Lolimot->addDimension(List_Model_Input[i].Name.c_str(), InputMin[i], InputMax[i], Sigma, NbCut);
	} /* End For */
      
      for(i=0; i<List_Model_Input.size(); i++)
	{
	  if (List_Model_Input[i].DontCut) Lolimot->setDontCutVar(i, true);
	} /* End For */
    } /* End If */
      
  for(i=0; i<MeasureData.size(); i++)
    {
      Mesure = Lolimot->addMesure(MeasureData[i][Measure],  MeasureMin[Measure],  MeasureMax[Measure]);
      
      for(j=0; j<List_Model_Input.size(); j++)
	{
	  Mesure->setDimensionValue(j, InputData[i][j]);
	} /* End For */
    } /* End For */

  if (WeightData.size()!=0)
    {
      for(i=0; i<MeasureData.size(); i++)
	{
	  Lolimot->addWeight(WeightData[i]);
	} /* End For */
    } /* End If */
  else
    {
      Lolimot->useDefaultWeight();
    } /* End Else */

  if (Loaded_Lolimot)
    {
      float          residu;
      LL_Partition * AuxPart = NULL;

      Lolimot->updateCoefficients(residu, AuxPart);
    } /* End If */

  // If a Model has been already loaded, we continue the learning
  if (Loaded_Lolimot) Lolimot->setContinue(true);

  if (!Lolimot->optimise(true))
    {
      cerr << "Apprentissage - Problem during the learning of the LOLIMOT network" << endl;
      exit(1);
    } /* End If */
}

void Learning_Bootstrap_Light(LL_Lolimot * Lolimot,
			      vector<vector<float> > & InputData,
			      vector<vector<float> > & MeasureData,
			      vector<float> & WeightData,
			      vector<Model_Input> & List_Model_Input,
			      vector<float> & MeasureMin,
			      vector<float> & MeasureMax, 
			      vector<float> & InputMin,
			      vector<float> & InputMax,
			      float Sigma, int NbCut,
			      float Bootstrap_Proportion,
			      unsigned int NbPartToExplore,
			      int Measure, 
			      string ResultDir,
			      string Filename,
			      bool Loaded_Lolimot)
{
  vector<vector<float> >      Input_TrainingData, Input_ValidationData;
  vector<vector<float> >      Measure_TrainingData, Measure_ValidationData;
  unsigned int                NbValidationPoints, countNbValidationPoints;
  vector<float>               X_Test, Y_LL;
  vector<bool>                PointForValidation;
  float                       Residual, BestResidual;
  int                         Index, Old_NbPartition;
  LL_Mesure                 * Mesure      = NULL;
  LL_Lolimot                * BestLolimot = NULL;
  ofstream                  * OutFile     = NULL;
  unsigned int                i, j, ii, ij;

  OutFile = new ofstream;
  OutFile->open((ResultDir + Filename).c_str());

  BestLolimot = new LL_Lolimot;

  cout << "Learning Bootstrap Light - Apprentissage du réseau LOLIMOT" << endl;

  Input_TrainingData.resize(0);
  Input_ValidationData.resize(0);
  Measure_TrainingData.resize(0);
  Measure_ValidationData.resize(0);

  PointForValidation.resize(InputData.size(), false);

  // We produce 2 data set: Learning and Validation
  NbValidationPoints = (int)(Bootstrap_Proportion*InputData.size());

  countNbValidationPoints = 0;
  BestResidual            = numeric_limits<float>::max();

  cout << "Learning Bootstrap Light - Generation of the validation file" << endl << endl;

  while(countNbValidationPoints<NbValidationPoints)
    {
      cout << UP << "Extraction of the Validation / Learning data - number of computed points : ";
      cout << countNbValidationPoints << " / " << NbValidationPoints << endl;

      do
	{
	  Index = (int)(rand()/(float)RAND_MAX * (InputData.size() - 1));
	}
      while(PointForValidation[Index]);

      PointForValidation[Index] = true;

      countNbValidationPoints++;
    } /* End While */

  for(i=0; i<InputData.size(); i++)
    {
      if (PointForValidation[i])
	{
	  Input_ValidationData.push_back(InputData[i]);
	  Measure_ValidationData.push_back(MeasureData[i]);
	} /* End If */
      else
	{
	  Input_TrainingData.push_back(InputData[i]);
	  Measure_TrainingData.push_back(MeasureData[i]);
	} /* End Else */
    } /* End For */

  cout << "Learning Bootstrap Light - Initialisation of the LOLIMOT network" << endl;

  if (!Loaded_Lolimot)
    {
      for(i=0; i<List_Model_Input.size(); i++)
	{
	  Lolimot->addDimension(List_Model_Input[i].Name.c_str(), InputMin[i], InputMax[i], Sigma, NbCut);
	} /* End For */
      
      for(i=0; i<List_Model_Input.size(); i++)
	{
	  if (List_Model_Input[i].DontCut) Lolimot->setDontCutVar(i, true);
	} /* End For */
    } /* End If */

  for(i=0; i<Measure_TrainingData.size(); i++)
    {
      Mesure = Lolimot->addMesure(Measure_TrainingData[i][Measure],  MeasureMin[Measure],  MeasureMax[Measure]);
      
      for(j=0; j<List_Model_Input.size(); j++)
	{
	  Mesure->setDimensionValue(j, Input_TrainingData[i][j]);
	} /* End For */
    } /* End For */

  if (WeightData.size()!=0)
    {
      for(i=0; i<MeasureData.size(); i++)
	{
	  Lolimot->addWeight(WeightData[i]);
	} /* End For */
    } /* End If */
  else
    {
      Lolimot->useDefaultWeight();
    } /* End Else */

  Old_NbPartition = Lolimot->getPartitionSet().size();

  if (!Loaded_Lolimot) Lolimot->setNbMaxPartitions(0);
  else                 Lolimot->setNbMaxPartitions(Old_NbPartition);

  if (Loaded_Lolimot)
    {
      float          residu;
      LL_Partition * AuxPart = NULL;

      Lolimot->updateCoefficients(residu, AuxPart);
    } /* End If */

  for(i=Old_NbPartition; i<NbPartToExplore; i++)
    {
      if (!Lolimot->step_optimise(true))
	{
	  cerr << "Learning - BootstrapLight - Problem during the learning of the LOLIMOT network" << endl;
	  exit(1);
	} /* End If */

      // Compute the validation residual
      cout << "Learning Bootstrap Light - Validation on a model which has " << i + 2 << " partitions" << endl;

      X_Test.resize(List_Model_Input.size());
      Y_LL.resize(Measure_ValidationData.size());
      
      for(ii=0; ii<Input_ValidationData.size(); ii++)
	{
	  for(ij=0; ij<List_Model_Input.size(); ij++)
	    {
	      X_Test[ij] = Input_ValidationData[ii][ij];
	    } /* End For */
	  
	  Y_LL[ii] = Lolimot->calculeF_Untransformed(X_Test);
	} /* End For */

      Residual = 0.0;

      for(j=0;j<Y_LL.size(); j++)
	{
	  Residual += pow(Y_LL[j] - Measure_ValidationData[j][Measure],(float)2.0);
	} /* End For j */

      Residual = Residual / (float)Y_LL.size();

      (*OutFile) << Residual << endl;

      // We store the lolimot model which has the best validation residual

      if (Residual<BestResidual)
	{
	  cout << "Learning Bootstrap Light - A good model has been found : Residual = ";
	  cout << Residual << "/" << BestResidual << endl;

	  BestResidual   = Residual;
	  (*BestLolimot) = (*Lolimot);
	} /* End If */
      else
	{
	  cout << "Learning Bootstrap Light - No good model found : Residual = ";
	  cout << Residual << "/" << BestResidual << endl;
	} /* End Else */

      if (Lolimot->isThereAnyPartitionAvailable()==false) break;
    } /* End For i */

  // We return the best value

  if (BestLolimot!=NULL)
    {
      (*Lolimot) = (*BestLolimot);
    } /* End If */

  if (BestLolimot) delete BestLolimot;

  OutFile->close();
  if (OutFile) delete OutFile;
}

void Learning_Compute_R2(LL_Lolimot * Lolimot,
			 vector<vector<float> > & InputData,
			 vector<vector<float> > & MeasureData,
			 vector<float> & WeightData,
			 vector<Model_Input> & List_Model_Input,
			 vector<float> & MeasureMin,
			 vector<float> & MeasureMax, 
			 vector<float> & InputMin,
			 vector<float> & InputMax,
			 float Sigma,
			 int NbCut,
			 int Measure,
			 string ResultDir,
			 string Filename)
{
  vector<float>        X_Test, Y_LL, List_R2;
  vector<unsigned int> ListInputsOfGroup;
  vector<int>          ListGroup;
  float                Residual, BestResidual;
  LL_Lolimot         * BestLolimot = NULL;
  LL_Partition       * PartitionWithBiggestResidu = NULL;
  float                StdResidu, StdMeas, MeasMean, R2A, Residu;
  int                  NbInputs = 0;
  ofstream           * OutFile  = NULL;
  unsigned int         MinGroup, MaxGroup;
  unsigned int         i, j, GroupInputs;

  MinGroup = numeric_limits<unsigned int>::max();
  MaxGroup = numeric_limits<unsigned int>::min();

  ListGroup.resize(0);
  List_R2.resize(0);

  // We look for the interval of varation of the group indexes
  for(i=0; i<List_Model_Input.size(); i++)
    {
      if (MinGroup>List_Model_Input[i].GroupNb) MinGroup = List_Model_Input[i].GroupNb;
      if (MaxGroup<List_Model_Input[i].GroupNb) MaxGroup = List_Model_Input[i].GroupNb;
    } /* End For */

  OutFile = new ofstream;
  OutFile->open((ResultDir + Filename).c_str());

  cout << "Learning Compute R2 - Learning of the LOLIMOT network" << endl;

  BestLolimot  = new LL_Lolimot;

  BestResidual = numeric_limits<float>::max();

  // Compute the validation residual
  cout << "Learning Compute R2 - Validation Of the complete model" << endl;
      
  X_Test.resize(List_Model_Input.size());
  Y_LL.resize(MeasureData.size());
      
  for(i=0; i<InputData.size(); i++)
    {
      NbInputs = 0;
      
      for(j=0; j<List_Model_Input.size(); j++)
	{
	  X_Test[j] = InputData[i][j];
	} /* End For */
	  
      Y_LL[i] = Lolimot->calculeF_Untransformed(X_Test);
    } /* End For */

  MeasMean  = 0.0;
  StdResidu = 0.0;
  StdMeas   = 0.0;

  for(i=0; i<InputData.size(); i++)
    {	  
      MeasMean   += MeasureData[i][Measure];
    } /* End For */
      
  MeasMean   = MeasMean / (float)InputData.size();

  for(i=0; i<InputData.size(); i++)
    {	  
      StdResidu  += pow(MeasureData[i][Measure] - Y_LL[i], (float)2.0);
      StdMeas    += pow(MeasureData[i][Measure] - MeasMean, (float)2.0);
    } /* End For */
  
  StdResidu  = StdResidu  / (float)(InputData.size());
  StdMeas    = StdMeas    / (float)(InputData.size());
  
  R2A        = 1 - StdResidu / (float)StdMeas;
  
  List_R2.push_back(R2A);
  ListGroup.push_back(-1);

  (*OutFile) << -1 << " " << R2A << endl;

  cout << "Learning Compute R2 - R2 for the complete model = " << R2A << endl;
  
  // We now remove the inputs of the model one after the other. We just update the model
  // without touching the partitions

  for(GroupInputs=MinGroup; GroupInputs<=MaxGroup; GroupInputs++)
    {
      // We search for the inputs which belongs to the GroupInputs group
      ListInputsOfGroup.resize(0);
      for(i=0; i<List_Model_Input.size(); i++)
	{
	  if (List_Model_Input[i].GroupNb==GroupInputs) ListInputsOfGroup.push_back(i);
	} /* End For */

      // If there are no inputs in this group, we do nothing and continue
      if (ListInputsOfGroup.size()==0) continue;

      // Otherwise, we memorise the number of the group
      ListGroup.push_back(GroupInputs);

      // On froze the corresponding inputs
      for(i=0; i<ListInputsOfGroup.size(); i++)
	{
	  Lolimot->setInhibe(ListInputsOfGroup[i], true);
	  cout << "Learning Compute R2 - Model with input " << List_Model_Input[ListInputsOfGroup[i]].Name;
	  cout << " removed" << endl;
	} /* End For */
      
      Lolimot->updateCoefficients(Residu, PartitionWithBiggestResidu);

      // Compute the validation residual
      cout << "Learning Compute R2 - Validation" << endl;
            
      X_Test.resize(List_Model_Input.size());
      Y_LL.resize(MeasureData.size());

      for(i=0; i<InputData.size(); i++)
	{
	  NbInputs = 0;

	  for(j=0; j<List_Model_Input.size(); j++)
	    {
	      X_Test[j] = InputData[i][j];
	    } /* End For */
	  
	  Y_LL[i] = Lolimot->calculeF_Untransformed(X_Test);
	} /* End For */

      // We "unfroze" the corresponding inputs
      for(i=0; i<ListInputsOfGroup.size(); i++)
	{
	  Lolimot->setInhibe(ListInputsOfGroup[i], false);
	} /* End For */

      MeasMean  = 0.0;
      StdResidu = 0.0;
      StdMeas   = 0.0;

      for(i=0; i<InputData.size(); i++)
	{	  
	  MeasMean  += MeasureData[i][Measure];
	} /* End For */

      MeasMean = MeasMean / (float)InputData.size();

      for(i=0; i<InputData.size(); i++)
	{	  
	  StdResidu  += pow(MeasureData[i][Measure] - Y_LL[i], (float)2.0);
	  StdMeas    += pow(MeasureData[i][Measure] - MeasMean, (float)2.0);
	} /* End For */

      StdResidu  = StdResidu  / (float)(InputData.size());
      StdMeas    = StdMeas    / (float)(InputData.size());
      
      R2A        = 1 - StdResidu / (float)StdMeas;

      List_R2.push_back(R2A);

      (*OutFile) << GroupInputs << " " << R2A << endl;

      cout << "Learning Compute R2 - we remove " << ListInputsOfGroup.size() << " inputs - ";
      cout << " Belong to group : " << GroupInputs << endl;
      cout << "Learning Compute R2 - R2 = " << R2A << endl;

      // Computation of the residual which allow to get the best model
      Residual = 0.0;
      
      for(j=0;j<Y_LL.size(); j++)
	{
	  Residual += pow(Y_LL[j] - MeasureData[j][Measure], (float)2.0);
	} /* End For j */
      
      Residual = Residual / (float)Y_LL.size();
    } /* End For GroupInputs */

  // Display the information related to the R2 results
  for(i=1; i<List_R2.size(); i++)
    {
      cout << "R2 Model with input " << ListGroup[i] << " removed = " << List_R2[i] << endl;
    } /* End For */

  cout << "R2 Complete model = " << List_R2[0] << endl;

  // Updating the model
  Lolimot->updateCoefficients(Residu, PartitionWithBiggestResidu);
}

void Optimize_Sigma(LL_Lolimot * Lolimot,
		    vector<vector<float> > & InputData,
		    vector<vector<float> > & MeasureData,
		    vector<float> & WeightData,
		    vector<Model_Input> & List_Model_Input,
		    int Measure,
		    float SigmaMin,
		    float SigmaMax, 
		    float SigmaStep)
{
  vector<float> X_Test;
  vector<float> Y_LL;
  float         Residu = 0.0;
  LL_Partition * partitionWithBiggestResidu = NULL;
  vector<float> ResiduList, SigmaList;
  float         ResiduMin, Min_Sigma, old_sigma = 0.0;
  unsigned int   i, j;

  old_sigma = Lolimot->getSigma();

  X_Test.resize(List_Model_Input.size());
  Y_LL.resize(MeasureData.size());
  ResiduList.resize(0);
  SigmaList.resize(0);

  // Computation of the initial residual

  for(i=0; i<InputData.size(); i++)
    {
      for(j=0; j<List_Model_Input.size(); j++)
	{
	  X_Test[j] = InputData[i][j];
	} /* End For */
      
      Y_LL[i] = Lolimot->calculeF_Untransformed(X_Test);
    } /* End For */
  
  // We compute the residual
  
  Residu = 0.0;

  for(i=0; i<InputData.size(); i++)
    {
      Residu += pow(MeasureData[i][Measure] - Y_LL[i], (float)2.0);
    } /* End For */
  
  Residu = Residu / (float)InputData.size();

  ResiduMin = Residu;
  Min_Sigma = old_sigma;
  
  for(float sigma_step = SigmaMin; sigma_step<=SigmaMax; sigma_step+=SigmaStep)
    {
      SigmaList.push_back(sigma_step);

      Lolimot->setSigma(sigma_step);
      Lolimot->updatePartitions();
      Lolimot->updateCoefficients(Residu, partitionWithBiggestResidu);

      for(i=0; i<InputData.size(); i++)
	{
	  for(j=0; j<List_Model_Input.size(); j++)
	    {
	      X_Test[j] = InputData[i][j];
	    } /* End For */
	  
	  Y_LL[i] = Lolimot->calculeF_Untransformed(X_Test);
	} /* End For */

      // We compute the residual
      
      Residu = 0.0;

      for(i=0; i<InputData.size(); i++)
	{
	  Residu += pow(MeasureData[i][Measure] - Y_LL[i], (float)2.0);
	} /* End For */
      
      Residu = Residu / (float)InputData.size();

      ResiduList.push_back(Residu);
    } /* End For */

  cout << "Optimize_Sigma: List Of Residual values with respect to sigma values" << endl;

  for(i=0; i<ResiduList.size(); i++)
    {
      if (ResiduMin>ResiduList[i])
	{
	  ResiduMin = ResiduList[i];
	  Min_Sigma = SigmaList[i];
	} /* End If */
      cout << "Sigma = " << SigmaList[i] << " -> Residu = " << ResiduList[i] << endl;
    } /* End For */

  cout << "Best Sigma value found:        Sigma    = " << Min_Sigma << endl;
  cout << "Residual value for this sigma: Residual = " << ResiduMin << endl;
  cout << "Old Sigma value :              Sigma    = " << old_sigma << endl;
  cout << "The Lolimot model has been set up to the optimized value of Sigma" << endl;
}

void Learning_LeavOneOut(LL_Lolimot * Lolimot,
			 vector<vector<float> > & InputData,
			 vector<vector<float> > & MeasureData,
			 vector<float> & WeightData,
			 vector<Model_Input> & List_Model_Input,
			 vector<float> & InputMin,
			 vector<float> & InputMax,
			 vector<float> & MeasureMin,
			 vector<float> & MeasureMax,
			 int Measure,
			 string Filename)
{
  vector<float> ResiduList;
  vector<float> X_Test, Y_LL;
  LL_Partition * partitionWithBiggestResidu = NULL;
  LL_Mesure *    Mesure = NULL;
  ofstream       OutFile;
  float         Residu = 0.0;
  unsigned int   i, j, k;

  ResiduList.resize(0);
  Y_LL.resize(InputData.size());
  X_Test.resize(List_Model_Input.size());

  cout << "Learning Leave One Out" << endl << endl << endl;

  for(i=0; i<InputData.size(); i++)
    {
      Lolimot->cleanMesures();
      
      cout << UP << UP << "Learning Leave One Out - Removing point " << i << endl;
      
      for(j=0; j<InputData.size(); j++)
	{
	  if (i!=j)
	    {
	      Mesure = Lolimot->addMesure(MeasureData[j][Measure],  MeasureMin[Measure],  MeasureMax[Measure]);
	      
	      for(k=0; k<List_Model_Input.size(); k++)
		{
		  Mesure->setDimensionValue(k, InputData[j][k]);
		} /* End For */
	    } /* End If */
	} /* End For */

      if (WeightData.size()!=0)
	{
	  for(j=0; j<MeasureData.size(); j++)
	    {
	      if (j==i) continue;
	      Lolimot->addWeight(WeightData[j]);
	    } /* End For */
	} /* End If */
      else
	{
	  Lolimot->useDefaultWeight();
	} /* End Else */
      
      cout << "Learning Leave One Out - Updating the model - Computation of the residual " << endl;
      
      Lolimot->updateCoefficients(Residu, partitionWithBiggestResidu);
      
      for(i=0; i<InputData.size(); i++)
	{
	  for(j=0; j<List_Model_Input.size(); j++)
	    {
	      X_Test[j] = InputData[i][j];
	    } /* End For */
	  
	  Y_LL[i] = Lolimot->calculeF_Untransformed(X_Test);
	} /* End For */
      
      // We compute the residual
      
      Residu = 0.0;
	  
      for(i=0; i<InputData.size(); i++)
	{
	  Residu += pow(MeasureData[i][Measure] - Y_LL[i], (float)2.0);
	} /* End For */
      
      Residu = Residu / (float)InputData.size();
      
      ResiduList.push_back(Residu);
    } /* End For */
  
  cout << "Learning Leave One Out - We write the results in the file " << Filename << endl;
  
  OutFile.open(Filename.c_str());
  
  for(i=0; i<ResiduList.size(); i++)
    {
      OutFile << ResiduList[i] << endl;
    } /* End For */
  
  OutFile.close();

  Lolimot->cleanMesures();
  
  cout << "Learning Leave One Out - Get back to the original model " << endl;
  
  for(i=0; i<InputData.size(); i++)
    {
      Mesure = Lolimot->addMesure(MeasureData[i][Measure],  MeasureMin[Measure],  MeasureMax[Measure]);
      
      for(j=0; j<List_Model_Input.size(); j++)
	{
	  Mesure->setDimensionValue(j, InputData[i][j]);
	} /* End For */
    } /* End For */
  
  Lolimot->updateCoefficients(Residu, partitionWithBiggestResidu);
}

void Learning_Recursive(LL_Lolimot * Lolimot,
			vector<vector<float> > & InputData,
			vector<vector<float> > & MeasureData,
			vector<float> & WeightData,
			vector<Model_Input> & List_Model_Input,
			vector<float> & MeasureMin,
			vector<float> & MeasureMax, 
			vector<float> & InputMin,
			vector<float> & InputMax,
			float Sigma, int NbCut,
			unsigned int NbPartToExplore_Phase1,
			unsigned int NbPartToExplore_Phase2,
			int Measure, 
			bool Loaded_Lolimot)
{
  vector<vector<float> > Input_TrainingData;
  vector<vector<float> > Measure_TrainingData;
  vector<vector<float> > Temporize;
  vector<unsigned int>   Index_Tempo;
  vector<float>          X_Test;
  vector<vector<float> > Y_LL;
  int                    Old_NbPartition;
  LL_Mesure            * Mesure = NULL;
  unsigned int           i, j, k, ii, ij;
  float                  Y_Res = 0.0, Residu = 0.0;
  LL_Partition         * partitionWithBiggestResidu = NULL;

  cout << "Learning Recursive - Learning of the LOLIMOT network" << endl;

  Input_TrainingData.resize(0);
  Measure_TrainingData.resize(0);

  Init_Temporize(Temporize, Index_Tempo, List_Model_Input);

  // Copy of InputData and MeasureData file
  for(i=0; i<InputData.size(); i++)
    {
      Input_TrainingData.push_back(InputData[i]);
      Measure_TrainingData.push_back(MeasureData[i]);
    } /* End For */
  
  cout << "Learning Recursive - Initialisation of the LOLIMOT network" << endl;

  if (!Loaded_Lolimot)
    {
      for(i=0; i<List_Model_Input.size(); i++)
	{
	  Lolimot->addDimension(List_Model_Input[i].Name.c_str(), InputMin[i], InputMax[i], Sigma, NbCut);
	} /* End For */
      
      for(i=0; i<List_Model_Input.size(); i++)
	{
	  if (List_Model_Input[i].DontCut) Lolimot->setDontCutVar(i, true);
	} /* End For */
    } /* End If */

  for(i=0; i<Measure_TrainingData.size(); i++)
    {
      Mesure = Lolimot->addMesure(Measure_TrainingData[i][Measure],  MeasureMin[Measure],  MeasureMax[Measure]);
      
      for(j=0; j<List_Model_Input.size(); j++)
	{
	  Mesure->setDimensionValue(j, Input_TrainingData[i][j]);
	} /* End For */
    } /* End For */

  if (WeightData.size()!=0)
    {
      for(i=0; i<MeasureData.size(); i++)
	{
	  Lolimot->addWeight(WeightData[i]);
	} /* End For */
    } /* End If */
  else
    {
      Lolimot->useDefaultWeight();
    } /* End Else */

  Old_NbPartition = Lolimot->getPartitionSet().size();

  if (!Loaded_Lolimot)
    {
      Lolimot->cleanPartitions();
      Lolimot->setNbMaxPartitions(0);
    }
  else Lolimot->setNbMaxPartitions(Old_NbPartition);

  Old_NbPartition = Lolimot->getPartitionSet().size();

  cout << "Learning - Recursive - Learning of the LOLIMOT network phase 1 - Using the measures" << endl;
    
  for(i=Old_NbPartition; i<NbPartToExplore_Phase1; i++)
    {
      if (!Lolimot->step_optimise(true))
	{
	  cerr << "Learning - Recursive - Problem during the learning of the LOLIMOT network" << endl;
	  exit(1);
	} /* End If */
    } /* End For */

  cout << "Learning - Recursive - Learning of the LOLIMOT network phase 2 - Updating the partitions" << endl;

  // Initialisation of the temporisation of the output
  X_Test.resize(List_Model_Input.size());

  Y_LL.resize(Input_TrainingData.size());
  for(i=0; i<Input_TrainingData.size(); i++)
    {
      Y_LL[i].resize(Temporize.size(), 0.0);
    } /* End For */

  for(ii=0; ii<Input_TrainingData.size(); ii++)
    {
      for(ij=0; ij<List_Model_Input.size(); ij++)
	{
	  X_Test[ij] = Input_TrainingData[ii][ij];
	} /* End For */
      
      Y_Res = Lolimot->calculeF_Untransformed(X_Test);
      
      Update_Temporize(Temporize, Y_Res);
      
      // Get the temporized values and store them in Y_LL
      for(j=0; j<Temporize.size(); j++)
	{
	  Y_LL[ii][j] = Temporize[j][Temporize[j].size()-1];
	} /* End For */
    } /* End For */

  for(k=0; k<Measure_TrainingData.size(); k++)
    {
      Mesure = Lolimot->getMesure(k);
      
      for(j=0; j<Temporize.size(); j++)
	{
	  Mesure->setDimensionValue(Index_Tempo[j], Y_LL[k][j]);
	  Input_TrainingData[k][Index_Tempo[j]] = Y_LL[k][j];
	} /* End For */
    } /* End For */
  
  Lolimot->updateCoefficients(Residu, partitionWithBiggestResidu);

  cout << "Learning - Recursive - Learning of the LOLIMOT network phase 2 - Using the estimations" << endl;

  for(i=0; i<NbPartToExplore_Phase2; i++)
    {
      // We prepare the estimation column and the looped inputs
      for(ii=0; ii<Input_TrainingData.size(); ii++)
	{
	  for(ij=0; ij<List_Model_Input.size(); ij++)
	    {
	      X_Test[ij] = Input_TrainingData[ii][ij];
	    } /* End For */
	  
	  Y_Res = Lolimot->calculeF_Untransformed(X_Test);

	  Update_Temporize(Temporize, Y_Res);
	    
	  // Get the temporized values and store them in Y_LL
	  for(j=0; j<Temporize.size(); j++)
	    {
	      Y_LL[ii][j] = Temporize[j][Temporize[j].size()-1];
	    } /* End For */
	} /* End For */

      for(k=0; k<Measure_TrainingData.size(); k++)
	{
	  Mesure = Lolimot->getMesure(k);
	  
	  for(j=0; j<Temporize.size(); j++)
	    {
	      Mesure->setDimensionValue(Index_Tempo[j], Y_LL[k][j]);
	      Input_TrainingData[k][Index_Tempo[j]] = Y_LL[k][j];
	    } /* End For */
	} /* End For */

      // Add a new partition
      if (!Lolimot->step_optimise(true))
	{
	  cerr << "Learning - Recursive - Problem during the learning of the LOLIMOT network" << endl;
	  exit(1);
	} /* End If */
    } /* End For */
}

void Learning_Update(LL_Lolimot * Lolimot,
		     vector<vector<float> > & InputData,
		     vector<vector<float> > & MeasureData,
		     vector<float> & WeightData,
		     vector<Model_Input> & List_Model_Input,
		     vector<float> & InputMin,
		     vector<float> & InputMax,
		     vector<float> & MeasureMin,
		     vector<float> & MeasureMax,
		     float Sigma,
		     int NbCut,
		     int Measure,
		     unsigned int NbIter,
		     bool Loaded_Lolimot)

{
  vector<vector<float> > Input_TrainingData;
  vector<vector<float> > Measure_TrainingData;
  vector<vector<float> > Temporize;
  vector<unsigned int>   Index_Tempo;
  vector<float>          X_Test;
  vector<vector<float> > Y_LL;
  LL_Mesure            * Mesure = NULL;
  unsigned int           i, j, k, ii, ij;
  float                  Y_Res = 0.0, Residu = 0.0;
  LL_Partition         * partitionWithBiggestResidu = NULL;

  X_Test.resize(List_Model_Input.size());
  Y_LL.resize(MeasureData.size());
  
  cout << "Learning Update - Learning of the LOLIMOT network" << endl;

  // Loading data
  if (!Loaded_Lolimot)
    {
      for(i=0; i<List_Model_Input.size(); i++)
	{
	  Lolimot->addDimension(List_Model_Input[i].Name.c_str(), InputMin[i], InputMax[i], Sigma, NbCut);
	} /* End For */
      
      for(i=0; i<List_Model_Input.size(); i++)
	{
	  if (List_Model_Input[i].DontCut) Lolimot->setDontCutVar(i, true);
	} /* End For */
    } /* End If */

  for(i=0; i<MeasureData.size(); i++)
    {
      Mesure = Lolimot->addMesure(MeasureData[i][Measure],  MeasureMin[Measure],  MeasureMax[Measure]);
      
      for(j=0; j<List_Model_Input.size(); j++)
	{
	  Mesure->setDimensionValue(j, InputData[i][j]);
	} /* End For */
    } /* End For */

  if (WeightData.size()!=0)
    {
      for(i=0; i<MeasureData.size(); i++)
	{
	  Lolimot->addWeight(WeightData[i]);
	} /* End For */
    } /* End If */
  else
    {
      Lolimot->useDefaultWeight();
    } /* End Else */

  Lolimot->updatePartitions();

  Input_TrainingData.resize(0);
  Measure_TrainingData.resize(0);

  Init_Temporize(Temporize, Index_Tempo, List_Model_Input);

  // Copy of InputData and MeasureData file
  for(i=0; i<InputData.size(); i++)
    {
      X_Test.resize(0);

      for(j=0; j<InputData[i].size(); j++)
	{
	  X_Test.push_back(InputData[i][j]);
	} /* End For */

      Input_TrainingData.push_back(X_Test);
    } /* End For */

  for(i=0; i<MeasureData.size(); i++)
    {
      X_Test.resize(0);

      for(j=0; j<MeasureData[i].size(); j++)
	{
	  X_Test.push_back(MeasureData[i][j]);
	} /* End For */

      Measure_TrainingData.push_back(X_Test);
    } /* End For */

  cout << "Apprentissage - Update - Using the estimations" << endl;

  X_Test.resize(List_Model_Input.size());
  Y_LL.resize(Input_TrainingData.size());

  for(i=0; i<Input_TrainingData.size(); i++)
    {
      Y_LL[i].resize(Temporize.size(), 0.0);
    } /* End For */

  for(i=0; i<NbIter; i++)
    {
      cout << "Apprentissage - Update - Iteration " << i << endl;

      for(ii=0; ii<Input_TrainingData.size(); ii++)
	{
	  for(ij=0; ij<List_Model_Input.size(); ij++)
	    {
	      X_Test[ij] = Input_TrainingData[ii][ij];
	    } /* End For */

	  Y_Res = Lolimot->calculeF_Untransformed(X_Test);

	  Update_Temporize(Temporize, Y_Res);
	    
	  // Get the temporized values and store them in Y_LL
	  for(j=0; j<Temporize.size(); j++)
	    {
	      Y_LL[ii][j] = Temporize[j][Temporize[j].size()-1];
	    } /* End For */
	} /* End For */

      for(k=0; k<Measure_TrainingData.size(); k++)
	{
	  Mesure = Lolimot->getMesure(k);

	  for(j=0; j<Temporize.size(); j++)
	    {
	      Mesure->setDimensionValue(Index_Tempo[j], Y_LL[k][j]);
	      Input_TrainingData[k][Index_Tempo[j]] = Y_LL[k][j];
	    } /* End For */
	} /* End For */

      Lolimot->updateCoefficients(Residu, partitionWithBiggestResidu);
    } /* End For */
}

void Learning_CrossValidation(LL_Lolimot * Lolimot,
			      vector<vector<float> > & InputData,
			      vector<vector<float> > & MeasureData,
			      vector<float> & WeightData,
			      vector<Model_Input> & List_Model_Input,
			      vector<float> & MeasureMin,
			      vector<float> & MeasureMax, 
			      vector<float> & InputMin,
			      vector<float> & InputMax,
			      unsigned int NbDataFiles,
			      int Measure, 
			      string ResultDir,
			      string Filename,
			      float & MeanResidual,
			      float & StdResidual)
{
  vector<vector<vector<float> > > Input_TrainAndValidData;
  vector<vector<vector<float> > > Measure_TrainAndValidData;
  vector<vector<float> >          Weight_TrainAndValidData;
  vector<unsigned int>            Selected_Pts, NotSelected_Pts;
  vector<float>                   X_Test, Y_LL;
  vector<bool>                    PointForValidation;
  float                           Residual;
  unsigned int                    Index, NbValidationPoints;
  ofstream                      * OutFile     = NULL;
  LL_Partition                  * AuxPart     = NULL;
  LL_Mesure                     * Mesure      = NULL;
  unsigned int                    i, j, k, l, ii, ij;

  OutFile = new ofstream;
  OutFile->open((ResultDir + Filename).c_str());

  cout << "Learning Cross-Validation" << endl;

  Input_TrainAndValidData.resize(NbDataFiles);
  Measure_TrainAndValidData.resize(NbDataFiles);
  Weight_TrainAndValidData.resize(NbDataFiles);

  for(i=0; i<NbDataFiles; i++)
    {
      Input_TrainAndValidData[i].resize(0);
      Measure_TrainAndValidData[i].resize(0);
      Weight_TrainAndValidData[i].resize(0);
    } /* End For */

  Selected_Pts.resize(InputData.size());
  NotSelected_Pts.resize(0);

  MeanResidual = 0.0;
  StdResidual  = 0.0;

  for(i=0; i<InputData.size(); i++)
    {
      Selected_Pts[i] = i;
    } /* End For */

  // We build 2 data set: Learning and Validation
  NbValidationPoints = (int)(InputData.size()/(float)NbDataFiles);

  cout << "Learning Cross-Validation - Generation of the validation file" << endl << endl;

  for(i=0; i<NbDataFiles; i++)
    {
      cout << "Extraction of the validation / Learning data - File : ";
      cout << i + 1 << " / " << NbDataFiles << endl;

      PointForValidation.clear();
      PointForValidation.resize(Selected_Pts.size(), false);

      for(j=0; j<(unsigned int)NbValidationPoints; j++)
	{
	  Index = (unsigned int)(rand()/(float)RAND_MAX * (Selected_Pts.size()));
	  if (Index==Selected_Pts.size()) Index--;
	  while (PointForValidation[Index])
	    {
	      Index = (unsigned int)(rand()/(float)RAND_MAX * (Selected_Pts.size()));
	      if (Index==Selected_Pts.size()) Index--;
	    } /* End While */

	  Input_TrainAndValidData[i].push_back(InputData[Selected_Pts[Index]]);
	  Measure_TrainAndValidData[i].push_back(MeasureData[Selected_Pts[Index]]);

	  if (WeightData.size()!=0)
	    {
	      Weight_TrainAndValidData[i].push_back(WeightData[Selected_Pts[Index]]);
	    } /* End If */

	  PointForValidation[Index] = true;
	} /* End For */

      NotSelected_Pts.clear();
      NotSelected_Pts.resize(0);

      for(j=0; j<Selected_Pts.size(); j++)
	{
	  if (!PointForValidation[j])
	    {
	      NotSelected_Pts.push_back(Selected_Pts[j]);
	    } /* End If */
	} /* End For */
      Selected_Pts = NotSelected_Pts;
    } /* End For */

  cout << "Learning Cross-Validation : Performing cross-validation" << endl;

  for(i=0; i<NbDataFiles; i++)
    {
      Lolimot->cleanMesures();

      for(j=0; j<NbDataFiles; j++)
	{
	  if (i!=j)
	    {
	      for(k=0; k<Measure_TrainAndValidData[j].size(); k++)
		{
		  Mesure = Lolimot->addMesure(Measure_TrainAndValidData[j][k][Measure],
					      MeasureMin[Measure],
					      MeasureMax[Measure]);
		  
		  for(l=0; l<List_Model_Input.size(); l++)
		    {
		      Mesure->setDimensionValue(l, Input_TrainAndValidData[j][k][l]);
		    } /* End For */
		} /* End For */
	      
	      if (WeightData.size()!=0)
		{
		  for(k=0; k<MeasureData.size(); k++)
		    {
		      Lolimot->addWeight(Weight_TrainAndValidData[j][k]);
		    } /* End For */
		} /* End If */
	    } /* End If */
	} /* End For */

      if (WeightData.size()==0)
	{
	  Lolimot->useDefaultWeight();
	} /* End Else */

      Lolimot->updateCoefficients(Residual, AuxPart);
	    
      // Compute the validation residual
      cout << "Learning Cross-Validation - Validation on the data set " << i + 1 << " / " << NbDataFiles << endl;

      X_Test.resize(List_Model_Input.size());
      Y_LL.resize(Measure_TrainAndValidData[i].size());
	      
      for(ii=0; ii<Input_TrainAndValidData[i].size(); ii++)
	{
	  for(ij=0; ij<List_Model_Input.size(); ij++)
	    {
	      X_Test[ij] = Input_TrainAndValidData[i][ii][ij];
	    } /* End For */
		  
	  Y_LL[ii] = Lolimot->calculeF_Untransformed(X_Test);
	} /* End For */
	      
      Residual = 0.0;
      
      for(ii=0; ii<Y_LL.size(); ii++)
	{
	  Residual += pow(Y_LL[ii] - Measure_TrainAndValidData[i][ii][Measure],(float)2.0);
	} /* End For j */
	      
      Residual = Residual / (float)Y_LL.size();

      MeanResidual += Residual;
      StdResidual  += Residual*Residual;

      (*OutFile) << Residual << endl;
    } /* End For */

  MeanResidual /= (float)NbDataFiles;
  StdResidual  /= (float)NbDataFiles;
  StdResidual  -= MeanResidual*MeanResidual;

  (*OutFile) << "Mean residual = " << MeanResidual << endl;
  (*OutFile) << "Std residual  = " << StdResidual << endl;

  OutFile->close();
  if (OutFile) delete OutFile;

  // Updating the model

  cout << "Learning Cross-Validation - Let's get back to the original model" << endl;

  Lolimot->cleanMesures();

  for(i=0; i<NbDataFiles; i++)
    {
      for(j=0; j<NbDataFiles; j++)
	{
	  for(k=0; k<Measure_TrainAndValidData[j].size(); k++)
	    {
	      Mesure = Lolimot->addMesure(Measure_TrainAndValidData[j][k][Measure],
					  MeasureMin[Measure],
					  MeasureMax[Measure]);
	      
	      for(l=0; l<List_Model_Input.size(); l++)
		{
		  Mesure->setDimensionValue(l, Input_TrainAndValidData[j][k][l]);
		} /* End For */
	    } /* End For */
	  
	  if (WeightData.size()!=0)
	    {
	      for(k=0; k<MeasureData.size(); k++)
		{
		  Lolimot->addWeight(WeightData[k]);
		} /* End For */
	    } /* End If */
	  else
	    {
	      Lolimot->useDefaultWeight();
	    } /* End Else */
	} /* End If */
    } /* End For */

  Lolimot->updateCoefficients(Residual, AuxPart);
}
