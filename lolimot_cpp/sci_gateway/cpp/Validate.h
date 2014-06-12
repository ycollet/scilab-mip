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

#ifndef VALIDATE_H
#define VALIDATE_H

#include <vector>

using namespace std;

#include <TrainLolimotStruct.h>
#include <LL_Lolimot.h>

void Validate(LL_Lolimot * Lolimot, vector<vector<float> > & InputData, vector<vector<float> > & MeasureData,
	      vector<float> & InputMin, vector<float> & InputMax,
	      vector<float> & MeasureMin, vector<float> & MeasureMax,
	      unsigned int Measure, string ModelName, vector<vector<float> > & TempoInput,
	      vector<Model_Input> & List_Model_Output,vector<Model_Input> & List_Model_Input,
	      vector<unsigned int> & List_Retro_Input, int R2_Start, int R2_End, bool Display,
	      bool CrossValid, unsigned int CrossValid_NbFiles, float CrossValid_MeanResidual, float CrossValid_StdResidual,
	      string SequenceName, float UseTransformLogEps, float UseTransformLogAlpha, bool Validation,
	      unsigned int NumOfValidateFile, bool Compute_C1, bool Compute_C2);
#endif
