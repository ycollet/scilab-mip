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
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <limits>

#include "Lasso/lasso.h"

using namespace std;

void ExportMinMaxMap(string Modelname, 
		     const vector<vector<float> > & MinData, 
		     const vector<vector<float> > & MaxData,
		     float SpeedEngMin,
		     float SpeedEngMax,
		     float TorqueMin,
		     float TorqueMax)
{
  ofstream   * OutFile = NULL;
  stringstream SFilename;

  OutFile = new ofstream;

  SFilename.str("");

  SFilename << Modelname << "_MaxMap.cpp";

  OutFile->open(SFilename.str().c_str());
  
  (*OutFile) << "float get_" << Modelname << "_MaxMap(float SpeedEng, float Torque)" << endl;
  (*OutFile) << "{" << endl;
  
  (*OutFile) << "  unsigned int Index_SE = 0;" << endl;
  (*OutFile) << "  unsigned int Index_T  = 0;" << endl;

  (*OutFile) << " float MaxData[" << MaxData.size() << "][" << MaxData[0].size() << "] = {" << endl;

  for(unsigned int i=0; i<MaxData.size(); i++)
    {
      (*OutFile) << "  ";
      for(unsigned int j=0; j<MaxData[i].size(); j++)
	{
	  (*OutFile) << ", " << MaxData[i][j];
	} /* End For */
      (*OutFile) << endl;
    } /* End For */
  (*OutFile) << "  };" << endl;

  (*OutFile) << "  Index_SE = (unsigned int)((SpeedEng - SpeedEngMin) / (SpeedEngMax - SpeedEngMin));" << endl;
  (*OutFile) << "  Index_T  = (unsigned int)((Torque - TorqueMin) / (TorqueMax - TorqueMin));" << endl;

  (*OutFile) << "  return MaxData[Index_SE][Index_T];" << endl;
  (*OutFile) << "}" << endl;

  SFilename.str("");

  SFilename << Modelname << "_MinMap.cpp";

  delete OutFile;

  OutFile = new ofstream;

  OutFile->open(SFilename.str().c_str());
  
  (*OutFile) << "float get_" << Modelname << "_MinMap(float SpeedEng, float Torque)" << endl;
  (*OutFile) << "{" << endl;
  
  (*OutFile) << "  unsigned int Index_SE = 0;" << endl;
  (*OutFile) << "  unsigned int Index_T  = 0;" << endl;

  (*OutFile) << " float MinData[" << MinData.size() << "][" << MinData[0].size() << "] = {" << endl;

  for(unsigned int i=0; i<MaxData.size(); i++)
    {
      for(unsigned int j=0; j<MaxData[i].size(); j++)
	{
	  (*OutFile) << ", " << MaxData[i][j];
	} /* End For */
      (*OutFile) << endl;
    } /* End For */
  (*OutFile) << "  };" << endl;

  (*OutFile) << "  Index_SE = (unsigned int)((SpeedEng - SpeedEngMin) / (SpeedEngMax - SpeedEngMin));" << endl;
  (*OutFile) << "  Index_T  = (unsigned int)((Torque - TorqueMin) / (TorqueMax - TorqueMin));" << endl;

  (*OutFile) << "  return MinData[Index_SE][Index_T];" << endl;
  (*OutFile) << "}" << endl;

  OutFile->close();

  delete OutFile;
}

void ExportMinMaxRegression(string Modelname, 
			    const vector<vector<float> > & MinData, 
			    const vector<vector<float> > & MaxData,
			    float SpeedEngMin,
			    float SpeedEngMax,
			    float TorqueMin,
			    float TorqueMax)
{
  ofstream    * OutFile = NULL;
  stringstream  SFilename;
  vector<float> MaxModel, MinModel;

  double * X = NULL;
  double * Y = NULL, * YHat = NULL;
  double * Residual = NULL;
  double * Param = NULL;
  int      Verbose = 0, PasSub = 0, PSuc, nbDimLasso = 3, NbPoints;
  double   Bound, Lagrangian, SpeedEng, Torque;

  MaxModel.resize(3, 0.0);
  MinModel.resize(3, 0.0);

  NbPoints = MinData.size()*MinData[0].size();

  X        = new double[3*NbPoints];
  Y        = new double[NbPoints];
  YHat     = new double[NbPoints];
  Residual = new double[NbPoints];
  Param    = new double[3];

  // Computation of the regression model
  
  // The Min model

  // Initialisation of the model
  for(unsigned int i=0; i<3; i++) Param[i] = 0.0;

  // We fill the Y matrix
  for(unsigned int i = 0; i < MinData.size(); i++)
    {
      for(unsigned int j = 0; j<MinData[i].size(); j++)
	{
	  Y[i*MinData[i].size() + j] = MinData[i][j];
	} /* End If */
    } /* End For */

  // We fill the X matrix
  for(unsigned int i=0; i<MinData.size(); i++)
    {
      for(unsigned int j=0; j<MinData[i].size(); j++)
	{
	  SpeedEng = (float)(i/(float)MinData.size() * (SpeedEngMax - SpeedEngMin) + SpeedEngMin);
	  Torque   = (float)(i/(float)MinData[i].size() * (TorqueMax - TorqueMin) + TorqueMin);

	  X[2*(i*MinData[i].size() + j) + 0] = 1.0;
	  X[2*(i*MinData[i].size() + j) + 1] = SpeedEng;
	  X[2*(i*MinData[i].size() + j) + 2] = Torque;
	} /* End For */
    } /* End For */

  Bound      = numeric_limits<float>::max();
  Lagrangian = 1.0;

  // Resolution using Lasso (module lasso extracted from R)
  lasso(X, &NbPoints, &nbDimLasso, &Bound, Param, Y, YHat, Residual, &Lagrangian, &PSuc,  &Verbose, &PasSub);

  MinModel[0] = Param[0]; // Coeff for the constant
  MinModel[1] = Param[1]; // Coeff for the speed
  MinModel[2] = Param[2]; // Coeff for the torque

  // The Max model

  // Initialisation of the model
  for(unsigned int i=0; i<3; i++) Param[i] = 0.0;

  // We fill the Y matrix
  for(unsigned int i=0; i<MaxData.size(); i++)
    {
      for(unsigned int j=0; j<MaxData[i].size(); j++)
	{
	  Y[i*MaxData[i].size() + j] = MaxData[i][j];
	} /* End If */
    } /* End For */

  // We fill the X matrix
  for(unsigned int i=0; i<MaxData.size(); i++)
    {
      for(unsigned int j=0; j<MaxData[i].size(); j++)
	{
	  SpeedEng = (float)(i/(float)MaxData.size() * (SpeedEngMax - SpeedEngMin) + SpeedEngMin);
	  Torque   = (float)(i/(float)MaxData[i].size() * (TorqueMax - TorqueMin) + TorqueMin);

	  X[2*(i*MaxData[i].size() + j) + 0] = 1.0;
	  X[2*(i*MaxData[i].size() + j) + 1] = SpeedEng;
	  X[2*(i*MaxData[i].size() + j) + 2] = Torque;
	} /* End For */
    } /* End For */

  Bound      = numeric_limits<float>::max();
  Lagrangian = 1.0;

  // Resolution using Lasso (module lasso extracted from R)
  lasso(X, &NbPoints, &nbDimLasso, &Bound, Param, Y, YHat, Residual, &Lagrangian, &PSuc,  &Verbose, &PasSub);

  MaxModel[0] = Param[0]; // Coeff for the constant
  MaxModel[1] = Param[1]; // Coeff for the speed
  MaxModel[2] = Param[2]; // Coeff for the torque

  // Cleaning the data set

  delete [] X;
  delete [] Y;
  delete [] YHat;
  delete [] Residual;
  delete [] Param;

  // Export of the regression model

  OutFile = new ofstream;

  SFilename.str("");

  SFilename << Modelname << "_MaxRegression.cpp";

  OutFile->open(SFilename.str().c_str());
  
  (*OutFile) << "float get_" << Modelname << "_MaxRegression(float SpeedEng, float Torque)" << endl;
  (*OutFile) << "{" << endl;
  (*OutFile) << "  return " << MaxModel[0] << " + " << MaxModel[1] << " * SpeedEng + " << MaxModel[2] << " * Torque;" << endl;
  (*OutFile) << "}" << endl;

  SFilename.str("");

  SFilename << Modelname << "_MinRegression.cpp";

  delete OutFile;

  OutFile = new ofstream;

  OutFile->open(SFilename.str().c_str());
  
  (*OutFile) << "float get_" << Modelname << "_MinRegression(float SpeedEng, float Torque)" << endl;
  (*OutFile) << "{" << endl;
  (*OutFile) << "  return " << MinModel[0] << " + " << MinModel[1] << " * SpeedEng + " << MinModel[2] << " * Torque;" << endl;
  (*OutFile) << "}" << endl;

  OutFile->close();

  delete OutFile;
}
