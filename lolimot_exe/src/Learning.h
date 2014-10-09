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

#ifndef LEARNING_H
#define LEARNING_H

#include <vector>

#include "LL_Lolimot.h"

struct Model_Input;

void Learning_NoMeth(LL_Lolimot * Lolimot,
		     std::vector<std::vector<float> > & InputData,
		     std::vector<std::vector<float> > & MeasureData,
		     std::vector<float> & WeightData,
		     std::vector<Model_Input> & List_Model_Input,
		     std::vector<float> & InputMin,
		     std::vector<float> & InputMax,
		     std::vector<float> & MeasureMin,
		     std::vector<float> & MeasureMax,
		     float Sigma,
		     int NbCut,
		     int Measure,
		     bool Loaded_Lolimot);

void Learning_Classic(LL_Lolimot * Lolimot,
		      std::vector<std::vector<float> > & InputData,
		      std::vector<std::vector<float> > & MeasureData,
		      std::vector<float> & WeightData,
		      std::vector<Model_Input> & List_Model_Input,
		      std::vector<float> & InputMin, 
		      std::vector<float> & InputMax,
		      std::vector<float> & MeasureMin,
		      std::vector<float> & MeasureMax,
		      float Sigma,
		      int NbCut, 
		      int Measure,
		      bool Loaded_Lolimot);

void Learning_Bootstrap_Light(LL_Lolimot * Lolimot,
			      std::vector<std::vector<float> > & InputData,
			      std::vector<std::vector<float> > & MeasureData,
			      std::vector<float> & WeightData,
			      std::vector<Model_Input> & List_Model_Input,
			      std::vector<float> & MeasureMin,
			      std::vector<float> & MeasureMax, 
			      std::vector<float> & InputMin,
			      std::vector<float> & InputMax,
			      float Sigma,
			      int NbCut,
			      float Bootstrap_Proportion,
			      unsigned int NbPartToExplore,
			      int Measure, 
			      string ResultDir,
			      string Filename,
			      bool Loaded_Lolimot);

void Learning_Compute_R2(LL_Lolimot * Lolimot,
			 std::vector<std::vector<float> > & InputData,
			 std::vector<std::vector<float> > & MeasureData,
			 std::vector<float> & WeightData,
			 std::vector<Model_Input> & List_Model_Input,
			 std::vector<float> & MeasureMin, 
			 std::vector<float> & MeasureMax, 
			 std::vector<float> & InputMin,
			 std::vector<float> & InputMax, 
			 float Sigma,
			 int NbCut,
			 int Measure,
			 string ResultDir,
			 string Filename);

void Optimize_Sigma(LL_Lolimot * Lolimot,
		    std::vector<std::vector<float> > & InputData, 
		    std::vector<std::vector<float> > & MeasureData,
		    std::vector<float> & WeightData,
		    std::vector<Model_Input> & List_Model_Input,
		    int Measure, 
		    float SigmaMin,
		    float SigmaMax, 
		    float SigmaStep);

void Learning_LeavOneOut(LL_Lolimot * Lolimot,
			 std::vector<std::vector<float> > & InputData,
			 std::vector<std::vector<float> > & MeasureData,
			 std::vector<float> & WeightData,
			 std::vector<Model_Input> & List_Model_Input,
			 std::vector<float> & InputMin,
			 std::vector<float> & InputMax,
			 std::vector<float> & MeasureMin,
			 std::vector<float> & MeasureMax,
			 int Measure,
			 string Filename);

void Learning_Recursive(LL_Lolimot * Lolimot,
			std::vector<std::vector<float> > & InputData,
			std::vector<std::vector<float> > & MeasureData,
			std::vector<float> & WeightData,
			std::vector<Model_Input> & List_Model_Input,
			std::vector<float> & MeasureMin,
			std::vector<float> & MeasureMax, 
			std::vector<float> & InputMin,
			std::vector<float> & InputMax,
			float Sigma, int NbCut,
			unsigned int NbPartToExplore_Phase1,
			unsigned int NbPartToExplore_Phase2,
			int Measure, 
			bool Loaded_Lolimot);

void Learning_Update(LL_Lolimot * Lolimot,
		     std::vector<std::vector<float> > & InputData,
		     std::vector<std::vector<float> > & MeasureData,
		     std::vector<float> & WeightData,
		     std::vector<Model_Input> & List_Model_Input,
		     std::vector<float> & InputMin,
		     std::vector<float> & InputMax,
		     std::vector<float> & MeasureMin,
		     std::vector<float> & MeasureMax,
		     float Sigma,
		     int NbCut,
		     int Measure,
		     unsigned int NbIter,
		     bool Loaded_Lolimot);

void Learning_CrossValidation(LL_Lolimot * Lolimot,
			      std::vector<std::vector<float> > & InputData,
			      std::vector<std::vector<float> > & MeasureData,
			      std::vector<float> & WeightData,
			      std::vector<Model_Input> & List_Model_Input,
			      std::vector<float> & MeasureMin,
			      std::vector<float> & MeasureMax, 
			      std::vector<float> & InputMin,
			      std::vector<float> & InputMax,
			      unsigned int NbDataFiles,
			      int Measure, 
			      string ResultDir,
			      string Filename,
			      float & MeanResidual,
			      float & StdResidual);
#endif
