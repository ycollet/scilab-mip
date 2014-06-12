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

#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>

using namespace std;

#include <TrainLolimotStruct.h>
#include <LL_Lolimot.h>
#include <Post.h>
#include <Analysis.h>

void Validate(LL_Lolimot * Lolimot, vector<vector<float> > & InputData, vector<vector<float> > & MeasureData,
	      vector<float> & InputMin, vector<float> & InputMax,
	      vector<float> & MeasureMin, vector<float> & MeasureMax,
	      unsigned int Measure, string ModelName, vector<vector<float> > & TempoInput,
	      vector<Model_Input> & List_Model_Output,vector<Model_Input> & List_Model_Input,
	      vector<unsigned int> & List_Retro_Input, int R2_Start, int R2_End, bool Display,
	      bool CrossValid, unsigned int CrossValid_NbFiles, float CrossValid_MeanResidual, float CrossValid_StdResidual,
	      string SequenceName, float UseTransformLogEps, float UseTransformLogAlpha, bool Validation,
	      unsigned int NumOfValidateFile, bool Compute_C1, bool Compute_C2)
{
  double StdResidu = 0.0, MeasMean = 0.0, StdMeas = 0.0, EstimMean = 0.0, StdEstim = 0.0;
  double CumEstim  = 0.0, CumMeas = 0.0, ResidualMean = 0.0, R2A = 0.0;
  stringstream FileName;
  ofstream * ExportFile = NULL;
  vector<float> X_Test, Y_LL;
  unsigned int L_Start, L_End;

  FileName.str("");
  
  if (!Validation)
    {
      if (List_Retro_Input.size()!=0)
	{
	  FileName << ModelName << "_Data_FB_Meas_" << SequenceName;
	  FileName << "_" << List_Model_Output[Measure].Name << ".dat";
	} /* End If */
      else
	{
	  FileName << ModelName << "_Data_" << SequenceName;
	  FileName << "_" << List_Model_Output[Measure].Name << ".dat";
	} /* End Else */
    } /* End If */
  else
    {
      if (List_Retro_Input.size()!=0)
	{
	  FileName << ModelName << "_Data_FB_Meas_" << SequenceName << "_" << NumOfValidateFile;
	  FileName << "_" << List_Model_Output[Measure].Name << ".dat";
	} /* End If */
      else
	{
	  FileName << ModelName << "_Data_" << SequenceName << "_" << NumOfValidateFile;
	  FileName << "_" << List_Model_Output[Measure].Name << ".dat";
	} /* End Else */
    } /* End Else */
  
  ExportFile = new ofstream;
  
  ExportFile->open(FileName.str().c_str(), ios_base::out);
  
  cout << SequenceName << " - Test of the LOLIMOT network" << endl;
  
  X_Test.resize(List_Model_Input.size());
  Y_LL.resize(MeasureData.size());
  
  cout << endl << endl;
  cout << "Test of the looped network: (if there are loops, then we perform the loops on the measures)" << endl;
  
  for(unsigned int i=0; i<InputData.size(); i++)
    {	      
      for(unsigned int j=0; j<List_Model_Input.size(); j++)
	{
	  X_Test[j] = InputData[i][j];
	} /* End For */
	  
      Y_LL[i] = Lolimot->calculeF_Untransformed(X_Test);
    } /* End For */
  
  StdResidu    = 0.0;
  MeasMean     = 0.0;
  StdMeas      = 0.0;
  EstimMean    = 0.0;
  StdEstim     = 0.0;
  CumEstim     = 0.0;
  CumMeas      = 0.0;
  ResidualMean = 0.0;
  
  FileName.str("");
  
  L_Start = 0;
  L_End   = InputData.size();
  if (R2_Start!=-1) L_Start = R2_Start;
  else              L_Start = 0;
  if (R2_End!=-1)   L_End   = InputData.size() - R2_End;
  else              L_End   = InputData.size();
  
  for(unsigned int i=0; i<InputData.size(); i++)
    {
      if (Display)
	{
	  for(unsigned int j=0; j<InputData[i].size(); j++)
	    {
	      cout << InputData[i][j] << " ";
	    } /* End For */
	  cout << " = (Measure / Estime) - (" << MeasureData[i][Measure] << " / " << Y_LL[i] << ")" << endl;
	} /* End If */
      
      for(unsigned int j=0; j<InputData[i].size(); j++)
	{
	  (*ExportFile) << InputData[i][j] << " ";
	} /* End For */
      
      (*ExportFile) << MeasureData[i][Measure] << " " << Y_LL[i] << endl;
    } /* End For */	      
  
  for(unsigned int i=0; i<InputData.size(); i++)
    {
      if ((i>=L_Start)&&(i<=L_End))
	{
	  CumEstim   += Y_LL[i];
	  CumMeas    += MeasureData[i][Measure];
	  
	  MeasMean     += MeasureData[i][Measure];
	  EstimMean    += Y_LL[i];
	  ResidualMean += MeasureData[i][Measure] - Y_LL[i];
	    } /* End If */
    } /* End For */
  
  MeasMean     = MeasMean   / (float)(L_Start - L_End);
  EstimMean    = EstimMean  / (float)(L_Start - L_End);
  ResidualMean = ResidualMean / (float)(L_Start - L_End);
  
  for(unsigned int i=0; i<InputData.size(); i++)
    {
      if ((i>=L_Start)&&(i<=L_End))
	{
	  StdResidu  += pow(MeasureData[i][Measure] - Y_LL[i], (float)2.0);
	  StdMeas    += pow(MeasureData[i][Measure] - MeasMean, (float)2.0);
	  StdEstim   += pow(Y_LL[i] - EstimMean, (float)2.0);
	} /* End If */
    } /* End For */
  
  StdResidu  = StdResidu  / (float)(L_Start - L_End);
  StdMeas    = StdMeas    / (float)(L_Start - L_End);
  
  R2A = 1 - StdResidu / (float)StdMeas;
  
  ExportFile->close();
  delete ExportFile;
  
  FileName.str("");
  
  if (!Validation)
    {
      if (List_Retro_Input.size()!=0)
	{
	  FileName << ModelName << "_Data_FB_Meas_Sumup_" << SequenceName;
	  FileName << "_" << List_Model_Output[Measure].Name << ".dat";
	} /* End If */
      else
	{
	  FileName << ModelName << "_Data_Sumup_" << SequenceName;
	  FileName << "_" << List_Model_Output[Measure].Name << ".dat";	      
	} /* End Else */
    } /* End If */
  else
    {
      if (List_Retro_Input.size()!=0)
	{
	  FileName << ModelName << "_Data_FB_Meas_Sumup_" << SequenceName << "_" << NumOfValidateFile;
	  FileName << "_" << List_Model_Output[Measure].Name << ".dat";
	} /* End If */
      else
	{
	  FileName << ModelName << "_Data_Sumup_" << SequenceName << "_" << NumOfValidateFile;
	  FileName << "_" << List_Model_Output[Measure].Name << ".dat";
	} /* End Else */
    } /* End Else */
  
  ExportFile = new ofstream;
  
  ExportFile->open(FileName.str().c_str(), ios_base::out);
  
  if (!Validation)
    {
      if (Compute_C1)
	{
	  vector<float> Result_C1;
	  Analysis_Response_NL(Lolimot, Result_C1);
	  for(unsigned int i=0; i<Result_C1.size(); i++)
	    {
	      cout << SequenceName << " - C1(" << i << ") = " << Result_C1[i] << endl;
	      (*ExportFile) << SequenceName << " - C1(" << i << ") = " << Result_C1[i] << endl;
	    } /* End For */
	} /* End If */
      
      if (Compute_C2)
	{
	  vector<float> Result_C2;
	  Analysis_Response(Lolimot, Result_C2, InputMin, InputMax);
	  for(unsigned int i=0; i<Result_C2.size(); i++)
	    {
	      cout << SequenceName << " - C2(" << i << ") = " << Result_C2[i] << endl;
	      (*ExportFile) << SequenceName << " - C2(" << i << ") = " << Result_C2[i] << endl;
	    } /* End For */
	} /* End If */
    } /* End If */
  
  cout << SequenceName << " - Initial residual : " << Lolimot->getInitialResidu()            << endl;
  cout << SequenceName << " - Final residual   : " << Lolimot->getResidu()                   << endl;
  cout << SequenceName << " - Std Residual     : " << sqrt(StdResidu)                        << endl;
  cout << SequenceName << " - Mean Measure     : " << MeasMean                               << endl;
  cout << SequenceName << " - Mean Estim       : " << EstimMean                              << endl;
  cout << SequenceName << " - R2               : " << R2A                                    << endl;
  cout << SequenceName << " - Cumulated error  : " << (CumEstim - CumMeas)                   << endl;
  cout << SequenceName << " - Cumul in %       : " << fabs(CumEstim - CumMeas)/CumMeas*100.0 << endl;
  cout << SequenceName << " - Cumul Estim      : " << CumEstim                               << endl;
  cout << SequenceName << " - Cumul Measure    : " << CumMeas                                << endl;
  
  if (CrossValid)
    {
      cout << SequenceName << " - Nb of divisions : " << CrossValid_NbFiles      << endl;
      cout << SequenceName << " - Mean CrossValid : " << CrossValid_MeanResidual << endl;
      cout << SequenceName << " - Std  CrossValid : " << CrossValid_StdResidual  << endl;
    } /* End If */
  
  (*ExportFile) << SequenceName << " - Initial residual : " << Lolimot->getInitialResidu()            << endl;
  (*ExportFile) << SequenceName << " - Final residual   : " << Lolimot->getResidu()                   << endl;
  (*ExportFile) << SequenceName << " - Std residual     : " << sqrt(StdResidu)                       << endl;
  (*ExportFile) << SequenceName << " - Mean Measure     : " << MeasMean                              << endl;
  (*ExportFile) << SequenceName << " - Mean Estim       : " << EstimMean                             << endl;
  (*ExportFile) << SequenceName << " - R2               : " << R2A                                   << endl;
  (*ExportFile) << SequenceName << " - Cumulated error  : " << (CumEstim - CumMeas)                  << endl;
  (*ExportFile) << SequenceName << " - Cumul in %       : " << fabs(CumEstim - CumMeas)/CumMeas*100.0 << endl;
  (*ExportFile) << SequenceName << " - Cumul Estim      : " << CumEstim                              << endl;
  (*ExportFile) << SequenceName << " - Cumul Measure    : " << CumMeas                               << endl;
  
  if (CrossValid)
    {
      (*ExportFile) << SequenceName << " - Nb of divisions : " << CrossValid_NbFiles      << endl;
      (*ExportFile) << SequenceName << " - Mean CrossValid : " << CrossValid_MeanResidual << endl;
      (*ExportFile) << SequenceName << " - Std  CrossValid : " << CrossValid_StdResidual  << endl;
    } /* End If */
  
  ExportFile->close();
  delete ExportFile;

  // Second validation
  
  // Sizing the input temporizations

  if (List_Retro_Input.size()!=0)
    {
      FileName.str("");
      
      if (!Validation)
	{
	  FileName << ModelName << "_Data_FB_Estim_" << SequenceName;
	  FileName << "_" << List_Model_Output[Measure].Name << ".dat";
	} /* End If */
      else
	{
	  FileName << ModelName << "_Data_FB_Estim_" << SequenceName << "_" << NumOfValidateFile;
	  FileName << "_" << List_Model_Output[Measure].Name << ".dat";
	} /* End Else */
      
      ExportFile = new ofstream;
      
      ExportFile->open(FileName.str().c_str(), ios_base::out);

      cout << "Test of the looped network : loop performed on the estimations" << endl;
  
      for(unsigned int i=0; i<List_Retro_Input.size(); i++)
	{
	  TempoInput[List_Retro_Input[i]].resize(List_Model_Input[List_Retro_Input[i]].Tempo-1, 0.0);
	  // We initialize the feedback with the measure at t=0. This allows to limit the convergence 
	  // time to the equilibrium position
	  for(unsigned int j=0; j<TempoInput[List_Retro_Input[i]].size(); j++)
	    {
	      TempoInput[List_Retro_Input[i]][j] = MeasureData[0][Measure];
	    } /* End For */
	  // We initialize the looped outputs - We store the initial value
	  InputData[0][List_Retro_Input[i]] = Temporize(TempoInput[List_Retro_Input[i]], MeasureData[0][Measure]);
	} /* End For */
      
      for(unsigned int i=0; i<InputData.size(); i++)
	{
	  if (i!=0)
	    {
	      // We memorize the output to produce the looped outputs
	      for(unsigned int j=0; j<List_Retro_Input.size(); j++)
		{
		  InputData[i][List_Retro_Input[j]] = Temporize(TempoInput[List_Retro_Input[j]], Y_LL[i-1]);
		} /* End For */
	    } /* End If */
	  
	  for(unsigned int j=0; j<List_Model_Input.size(); j++)
	    {
	      X_Test[j] = InputData[i][j];
	    } /* End For */
	  
	  Y_LL[i] = Lolimot->calculeF_Untransformed(X_Test);
	} /* End For */
      
      StdResidu    = 0.0;
      MeasMean     = 0.0;
      StdMeas      = 0.0;
      EstimMean    = 0.0;
      StdEstim     = 0.0;
      CumEstim     = 0.0;
      CumMeas      = 0.0;
      ResidualMean = 0.0;
      
      L_Start = 0;
      L_End   = InputData.size();
      if (R2_Start!=-1) L_Start = R2_Start;
      else              L_Start = 0;
      if (R2_End!=-1)   L_End   = InputData.size() - R2_End;
      else              L_End   = InputData.size();
      
      for(unsigned int i=0; i<InputData.size(); i++)
	{
	  if (Display)
	    {
	      for(unsigned int j=0; j<InputData[i].size(); j++)
		{
		  cout << InputData[i][j] << " ";
		} /* End For */
	      cout << " = (Measure / Estime) - (" << MeasureData[i][Measure] << " / " << Y_LL[i] << ")" << endl;
	    } /* End If */
	  
	  for(unsigned int j=0; j<InputData[i].size(); j++)
	    {
	      (*ExportFile) << InputData[i][j] << " ";
	    } /* End For */
	  
	  (*ExportFile) << MeasureData[i][Measure] << " " << Y_LL[i] << endl;
	} /* End For */
      
      for(unsigned int i=0; i<InputData.size(); i++)
	{
	  if ((i>=L_Start)&&(i<=L_End))
	    {
	      CumEstim   += Y_LL[i];
	      CumMeas    += MeasureData[i][Measure];
	      
	      MeasMean     += MeasureData[i][Measure];
	      EstimMean    += Y_LL[i];
	      ResidualMean += MeasureData[i][Measure] - Y_LL[i];
	    } /* End If */
	} /* End For */
      
      MeasMean     = MeasMean   / (float)(L_Start - L_End);
      EstimMean    = EstimMean  / (float)(L_Start - L_End);
      ResidualMean = ResidualMean / (float)(L_Start - L_End);
      
      for(unsigned int i=0; i<InputData.size(); i++)
	{
	  if ((i>=L_Start)&&(i<=L_End))
	    {
	      StdResidu  += pow(MeasureData[i][Measure] - Y_LL[i], (float)2.0);
	      StdMeas    += pow(MeasureData[i][Measure] - MeasMean, (float)2.0);
	      StdEstim   += pow(Y_LL[i] - EstimMean, (float)2.0);
	    } /* End For */
	} /* End For */
      
      StdResidu  = StdResidu  / (float)(L_Start - L_End);
      StdMeas    = StdMeas    / (float)(L_Start - L_End);
      
      R2A = 1 - StdResidu / (float)StdMeas;

      ExportFile->close();
      delete ExportFile;

      FileName.str("");
      
      if (!Validation)
	{
	  FileName << ModelName << "_Data_FB_Estim_Sumup_" << SequenceName;
	  FileName << "_" << List_Model_Output[Measure].Name << ".dat";
	} /* End If */
      else
	{
	  FileName << ModelName << "_Data_FB_Estim_Sumup_" << SequenceName << "_" << NumOfValidateFile;
	  FileName << "_" << List_Model_Output[Measure].Name << ".dat";
	} /* End Else */
      
      ExportFile = new ofstream;
      
      ExportFile->open(FileName.str().c_str(), ios_base::out);

      cout << SequenceName << " - Initial residual : " << Lolimot->getInitialResidu()            << endl;
      cout << SequenceName << " - Final residual   : " << Lolimot->getResidu()                   << endl;
      cout << SequenceName << " - Std residual     : " << sqrt(StdResidu)                        << endl;
      cout << SequenceName << " - Mean Measure     : " << MeasMean                               << endl;
      cout << SequenceName << " - Mean Estim       : " << EstimMean                              << endl;
      cout << SequenceName << " - R2               : " << R2A                                    << endl;
      cout << SequenceName << " - Cumulated error  : " << (CumEstim - CumMeas)                   << endl;
      cout << SequenceName << " - Cumul in %       : " << fabs(CumEstim - CumMeas)/CumMeas*100.0 << endl;
      cout << SequenceName << " - Cumul Estim      : " << CumEstim                               << endl;
      cout << SequenceName << " - Cumul Measure    : " << CumMeas                                << endl;
      
      if (!Validation)
	{
	  if (CrossValid)
	    {
	      (*ExportFile) << SequenceName << " - Nb of divisions : " << CrossValid_NbFiles      << endl;
	      (*ExportFile) << SequenceName << " - Mean CrossValid : " << CrossValid_MeanResidual << endl;
	      (*ExportFile) << SequenceName << " - Std  CrossValid : " << CrossValid_StdResidual  << endl;
	    } /* End If */
	  
	  if (Compute_C1)
	    {
	      vector<float> Result_C1;
	      Analysis_Response_NL(Lolimot, Result_C1);
	      for(unsigned int i=0; i<Result_C1.size(); i++)
		{
		  cout << SequenceName << " - C1(" << i << ") = " << Result_C1[i] << endl;
		  (*ExportFile) << SequenceName << " - C1(" << i << ") = " << Result_C1[i] << endl;
		} /* End For */
	    } /* End If */
	  
	  if (Compute_C2)
	    {
	      vector<float> Result_C2;
	      Analysis_Response(Lolimot, Result_C2, InputMin, InputMax);
	      for(unsigned int i=0; i<Result_C2.size(); i++)
		{
		  cout << SequenceName << " - C2(" << i << ") = " << Result_C2[i] << endl;
		  (*ExportFile) << SequenceName << " - C2(" << i << ") = " << Result_C2[i] << endl;
		} /* End For */
	    } /* End If */
	} /* End If */
      
      ExportFile->close();
      
      delete ExportFile;
    } /* End If */
}
