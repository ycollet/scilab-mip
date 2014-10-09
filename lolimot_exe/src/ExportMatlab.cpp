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
#include <sstream>
#include <string>
#include <vector>

#include "LL_Lolimot.h"
#include "LL_Partition.h"
#include "LL_Dimension.h"
#include "LL_Mesure.h"
#include "Post.h"
#include "TrainLolimotStruct.h"

using namespace std;

void export_Destroy_Matlab(LL_Lolimot * Lolimot, ofstream & F_Lolimot)
{
  F_Lolimot << "  clear Ei;" << endl;
  F_Lolimot << "  clear Li;" << endl;
}

void export_Init_Data_Matlab(LL_Lolimot * Lolimot, ofstream & F_Lolimot)
{
  unsigned int nbPartition = Lolimot->getPartitionSet().size();
  unsigned int nbDimension = Lolimot->getDimensionSet().size();

  F_Lolimot << " Coeff     = zeros(" << nbPartition << ", " << nbDimension << ");" << endl;
  F_Lolimot << " Centre    = zeros(" << nbPartition << ", " << nbDimension << ");" << endl;
  F_Lolimot << " EcartType = zeros(" << nbPartition << ", " << nbDimension << ");" << endl;
  F_Lolimot << " Coeff0    = zeros(1, " << nbPartition << ");" << endl;
}

void export_Destroy_Data_Matlab(LL_Lolimot * Lolimot, ofstream & F_Lolimot)
{
  F_Lolimot << " clear Coeff;" << endl;
  F_Lolimot << " clear Centre;" << endl;
  F_Lolimot << " clear EcartType;" << endl;
  F_Lolimot << " clear Coeff0;" << endl;
}

void export_Declare_Variable_Matlab(LL_Lolimot * Lolimot, ofstream & F_Lolimot, const string & VarName)
{
  F_Lolimot << " " << VarName << " = 0.0;" << endl;
}

void export_Data_Matlab_Text(LL_Lolimot * Lolimot, ofstream & F_Lolimot)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  unsigned int   noDimension, nbDimension = Lolimot->getDimensionSet().size();
  LL_Partition * curPartition = NULL;

  F_Lolimot << endl;
  F_Lolimot << "    %%----------------------------" << endl;
  F_Lolimot << "    %% List of the partitions data" << endl;
  F_Lolimot << "    %%----------------------------" << endl;
  F_Lolimot << endl;

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];
      F_Lolimot << " Coeff0(" << noPartition + 1 << ") = " << curPartition->getCoeff0() << ";" << endl;
    } /* End For */

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];
      for(noDimension = 0; noDimension < nbDimension; noDimension ++)
	{
	  F_Lolimot << " Coeff(" << noPartition + 1<< ", " << noDimension + 1 << ") = ";
	  F_Lolimot << curPartition->getCoeff(noDimension) << ";" << endl;
	} /* End For */
    } /* End For */

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];
      for(noDimension = 0; noDimension < nbDimension; noDimension ++)
	{
	  F_Lolimot << " Centre(" << noPartition + 1 << ", " << noDimension + 1 << ") = ";
	  F_Lolimot << curPartition->getCentre(noDimension) << ";" << endl;
	} /* End For */
    } /* End For */

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];
      for(noDimension = 0; noDimension < nbDimension; noDimension ++)
	{
	  F_Lolimot << " EcartType(" << noPartition + 1 << ", " << noDimension + 1 << ") = ";
	  F_Lolimot << curPartition->getEcartType(noDimension) << ";" << endl;
	} /* End For */
    } /* End For */
}

void export_dV_x_Matlab_Text(LL_Lolimot * Lolimot, ofstream & F_Lolimot, const string & ParamName, const string & VarName)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();

  F_Lolimot << endl;
  F_Lolimot << "    %%-------------------------------------------------" << endl;
  F_Lolimot << "    %% Computation of the derivative of the denominator" << endl;
  F_Lolimot << "    %%-------------------------------------------------" << endl;
  F_Lolimot << endl;

  F_Lolimot << " " << VarName << " = 0.0;" << endl;

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      F_Lolimot << " " << VarName << " = " << VarName << " - ((x(" << ParamName << ") - Centre(";
      F_Lolimot << noPartition + 1 << ", " << ParamName << ")) / ";
      F_Lolimot << "(EcartType(" << noPartition + 1 << ", " << ParamName << ")*EcartType(";
      F_Lolimot << noPartition + 1 << ", " << ParamName << "))";
      F_Lolimot << " * Ei(" << noPartition + 1 << ");" << endl;
    } /* End For */
}

void export_dU_x_Matlab_Text(LL_Lolimot * Lolimot, ofstream & F_Lolimot, const string & ParamName, const string & VarName)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();

  F_Lolimot << endl;
  F_Lolimot << "    %%-----------------------------------------------" << endl;
  F_Lolimot << "    %% Computation of the derivative of the numerator" << endl;
  F_Lolimot << "    %%-----------------------------------------------" << endl;
  F_Lolimot << endl;

  F_Lolimot << " " << VarName << " = 0.0;" << endl;

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      F_Lolimot << " " << VarName << " = " << VarName << " + (Coeff(" << noPartition + 1 << ", " << ParamName << ") - ";
      F_Lolimot << "(x(" << ParamName << ") - Centre(" << noPartition + 1 << ", " << ParamName << ")) / ";
      F_Lolimot << "(EcartType(" << noPartition + 1 << ", " << ParamName << ")*EcartType(";
      F_Lolimot << noPartition + 1 << ", " << ParamName << "))";
      F_Lolimot << " * Li(" << noPartition + 1 << ")) * Ei(" << noPartition + 1 << ");" << endl;
    } /* End For */
}

void export_Init_Matlab(LL_Lolimot * Lolimot, ofstream & F_Lolimot)
{
  unsigned int nbPartition = Lolimot->getPartitionSet().size();

  F_Lolimot << "Ei = zeros(1, " << nbPartition << ");" << endl;
  F_Lolimot << "Li = zeros(1," << nbPartition << ");" << endl;
}

void export_Ei_x_Matlab(LL_Lolimot * Lolimot, ofstream & F_Lolimot)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  unsigned int   noDimension, nbDimension = Lolimot->getDimensionSet().size();
  LL_Partition * curPartition = NULL;

  F_Lolimot << endl;
  F_Lolimot << "    %%------------------------------------------" << endl;
  F_Lolimot << "    %% Pour chaque partition i, on calcule Ei(x)" << endl;
  F_Lolimot << "    %%------------------------------------------" << endl;
  F_Lolimot << endl;

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];

      F_Lolimot << " Ei(" << noPartition + 1 << ") = exp(-0.5 * (";

      for(noDimension = 0; noDimension < nbDimension; noDimension ++)
	{
	  F_Lolimot << " + ((x(" << noDimension + 1 << ") - " << curPartition->getCentre(noDimension) << ")^2.0) / ";
	  F_Lolimot <<  curPartition->getEcartType(noDimension) * curPartition->getEcartType(noDimension);
	} /* End For */

      F_Lolimot << "));" << endl;
    } /* End For */
}

void export_Li_x_Matlab(LL_Lolimot * Lolimot, ofstream & F_Lolimot)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  unsigned int   noDimension, nbDimension = Lolimot->getDimensionSet().size();
  LL_Partition * curPartition = NULL;

  F_Lolimot << endl;
  F_Lolimot << "    %%------------------------------------------" << endl;
  F_Lolimot << "    %% Pour chaque partition i, on calcule Li(x)" << endl;
  F_Lolimot << "    %%------------------------------------------" << endl;
  F_Lolimot << endl;

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];

      F_Lolimot << "    Li(" << noPartition + 1 << ") = " << curPartition->getCoeff0();

      for(noDimension = 0; noDimension < nbDimension; noDimension ++)
	{
	  F_Lolimot << " + " << curPartition->getCoeff(noDimension) << " * x(" << noDimension + 1 << ")";
	} /* End For */

      F_Lolimot << ";" << endl;
    } /* End For */
}

void export_U_x_Matlab(LL_Lolimot * Lolimot, ofstream & F_Lolimot, const string & VarName)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();

  F_Lolimot << endl;
  F_Lolimot << "    %%-----------------------------" << endl;
  F_Lolimot << "    %% Computation of the numerator" << endl;
  F_Lolimot << "    %%-----------------------------" << endl;
  F_Lolimot << endl;

  F_Lolimot << " " << VarName << " = 0.0;" << endl;

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      F_Lolimot << " " << VarName << " = " << VarName << " + Li(" << noPartition + 1;
      F_Lolimot << ") * Ei(" << noPartition + 1 << ");" << endl;
    } /* End For */
}

void export_V_x_Matlab(LL_Lolimot * Lolimot, ofstream & F_Lolimot, const string & VarName)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();

  F_Lolimot << endl;
  F_Lolimot << "    %%-------------------------------" << endl;
  F_Lolimot << "    %% Computation of the denominator" << endl;
  F_Lolimot << "    %%-------------------------------" << endl;
  F_Lolimot << endl;

  F_Lolimot << " " << VarName << " = 0.0;" << endl;

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      F_Lolimot << " " << VarName << " = " << VarName << " + Ei(" << noPartition + 1 << ");" << endl;
    } /* End For */
}

void export_dU_x_Matlab(LL_Lolimot * Lolimot, ofstream & F_Lolimot, unsigned int VarNb1, const string & VarName)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;

  // Warning: VarNb1 index must start from 1

  F_Lolimot << endl;
  F_Lolimot << "    %%-----------------------------------------------" << endl;
  F_Lolimot << "    %% Computation of the derivative of the numerator" << endl;
  F_Lolimot << "    %%-----------------------------------------------" << endl;
  F_Lolimot << endl;

  F_Lolimot << " " << VarName << " = 0.0;" << endl;

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];

      F_Lolimot << " " << VarName << " = " << VarName << " + (";
      F_Lolimot << curPartition->getCoeff(VarNb1 - 1) << " - ";
      F_Lolimot << "((x(" << VarNb1 << ") - " << curPartition->getCentre(VarNb1 - 1);
      F_Lolimot << ") / " << curPartition->getEcartType(VarNb1 - 1)*curPartition->getEcartType(VarNb1 - 1) << ") * ";
      F_Lolimot << "Li(" << noPartition + 1 << ")) * ";
      F_Lolimot << "Ei(" << noPartition + 1 << ");" << endl;
    } /* End For */
}

void export_dV_x_Matlab(LL_Lolimot * Lolimot, ofstream & F_Lolimot, unsigned int VarNb1, const string & VarName)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;

  // Warning: VarNb1 index must start from 1

  F_Lolimot << endl;
  F_Lolimot << "    %%-------------------------------------------------" << endl;
  F_Lolimot << "    %% Computation of the derivative of the denominator" << endl;
  F_Lolimot << "    %%-------------------------------------------------" << endl;
  F_Lolimot << endl;

  F_Lolimot << " " << VarName << " = 0.0;" << endl;

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];

      F_Lolimot << " " << VarName << " = " << VarName << " - ";
      F_Lolimot << "((x(" << VarNb1 << ") - " << curPartition->getCentre(VarNb1 - 1);
      F_Lolimot << ") / " << curPartition->getEcartType(VarNb1 - 1)*curPartition->getEcartType(VarNb1 - 1) << ")";
      F_Lolimot << " * Ei(" << noPartition + 1 << ");" << endl;
    } /* End For */
}

void export_List_Partitions_Matlab(LL_Lolimot * Lolimot, ofstream & F_Lolimot)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;

  F_Lolimot << endl;
  F_Lolimot << "    %%-------------------" << endl;
  F_Lolimot << "    %% List of partitions" << endl;
  F_Lolimot << "    %%-------------------" << endl;
  F_Lolimot << endl;

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];

      F_Lolimot << "    %% Partition n " << noPartition + 1 << " " << curPartition->printDimensions().c_str() << endl;
    } /* End For */
}

void export_List_Variables_Matlab(LL_Lolimot * Lolimot, ofstream & F_Lolimot)
{
  unsigned int noDimension, nbDimension = Lolimot->getDimensionSet().size();

  F_Lolimot << endl;
  F_Lolimot << "    %%------------------------------------------" << endl;
  F_Lolimot << "    %% Liste of variables name and of boundaries" << endl;
  F_Lolimot << "    %%------------------------------------------" << endl;
  F_Lolimot << endl;

  for(noDimension = 0; noDimension < nbDimension; noDimension ++)
    {
      F_Lolimot << "%% ";
      F_Lolimot << Lolimot->getDimensionSet()[noDimension]->getName() << " - [";
      F_Lolimot << Lolimot->getDimensionSet()[noDimension]->getMin()  << "  ";
      F_Lolimot << Lolimot->getDimensionSet()[noDimension]->getMax()  << "]" << endl;
    } /* End For */

  F_Lolimot << "%% ";
  F_Lolimot << "Measure - [" << Lolimot->getMesureSet()[0]->getMinValueExtrapolation();
  F_Lolimot << " " << Lolimot->getMesureSet()[0]->getMaxValueExtrapolation() << "]" << endl;
}

void export_ddU_x_Matlab(LL_Lolimot * Lolimot, ofstream & F_Lolimot, unsigned int VarNb1, unsigned int VarNb2, const string & VarName)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;

  // Warning: VarNb1 and VarNb2 indexes must start from 1

  F_Lolimot << endl;
  F_Lolimot << "    %%------------------------------------------------------" << endl;
  F_Lolimot << "    %% Computation of the second derivative of the numerator" << endl;
  F_Lolimot << "    %%------------------------------------------------------" << endl;
  F_Lolimot << endl;

  F_Lolimot << " " << VarName << " = 0.0;" << endl;

  if (VarNb1!=VarNb2)
    {
      for(noPartition = 0; noPartition < nbPartition; noPartition ++)
	{
	  curPartition = Lolimot->getPartitionSet()[noPartition];
	  
	  F_Lolimot << " " << VarName << " = " << VarName << " - (";
	  F_Lolimot << curPartition->getCoeff(VarNb2-1) << " * ";
	  F_Lolimot << "((x(" << VarNb1 << ") - " << curPartition->getCentre(VarNb1 - 1) << ") / ";
	  F_Lolimot << curPartition->getEcartType(VarNb1 - 1) * curPartition->getEcartType(VarNb1 - 1) << ") + ";
	  F_Lolimot << "(" << curPartition->getCoeff(VarNb1-1) << " - ";
	  F_Lolimot << "((x(" << VarNb1 << ") - " << curPartition->getCentre(VarNb1 - 1) << ") / ";
	  F_Lolimot << curPartition->getEcartType(VarNb1 - 1) * curPartition->getEcartType(VarNb1 - 1) << ") * ";
	  F_Lolimot << "Li(" << noPartition + 1 << ")) * ";
	  F_Lolimot << "((x(" << VarNb2 << ") - " << curPartition->getCentre(VarNb2 - 1) << ") / ";
	  F_Lolimot << curPartition->getEcartType(VarNb2 - 1) * curPartition->getEcartType(VarNb2 - 1) << ")) * ";
	  F_Lolimot << "Ei(" << noPartition + 1 << ");" << endl;
	} /* End For */
    } /* End If */
  else
    {
      for(noPartition = 0; noPartition < nbPartition; noPartition ++)
	{
	  curPartition = Lolimot->getPartitionSet()[noPartition];
	  
	  F_Lolimot << " " << VarName << " = " << VarName << " - (Li(" << noPartition + 1 << ") / ";
	  F_Lolimot << curPartition->getEcartType(VarNb1 - 1)*curPartition->getEcartType(VarNb1 - 1) << " + 2 * ";
	  F_Lolimot << curPartition->getCoeff(VarNb1 - 1) << " * ";
	  F_Lolimot << "((x(" << VarNb1 << ") - " << curPartition->getCentre(VarNb1 - 1) << ") / ";
	  F_Lolimot << curPartition->getEcartType(VarNb1 - 1) * curPartition->getEcartType(VarNb1 - 1) << ") - ";
	  F_Lolimot << "Li(" << noPartition + 1 << ") * ";
	  F_Lolimot << "((x(" << VarNb1 << ") - " << curPartition->getCentre(VarNb1 - 1) << ") / ";
	  F_Lolimot << curPartition->getEcartType(VarNb1 - 1) * curPartition->getEcartType(VarNb1 - 1) << ")^2.0) * ";
	  F_Lolimot << "Ei(" << noPartition + 1 << ");" << endl;
	} /* End For */
    } /* End Else */
}

void export_ddV_x_Matlab(LL_Lolimot * Lolimot, ofstream & F_Lolimot, unsigned int VarNb1, unsigned int VarNb2, const string & VarName)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;

  // Warning: VarNb1 and VarNb2 indexes must start from 1

  F_Lolimot << endl;
  F_Lolimot << "    %%--------------------------------------------------------" << endl;
  F_Lolimot << "    %% Computation of the second derivative of the denominator" << endl;
  F_Lolimot << "    %%--------------------------------------------------------" << endl;
  F_Lolimot << endl;

  F_Lolimot << " " << VarName << " = 0.0;" << endl;

  if (VarNb1!=VarNb2)
    {
      for(noPartition = 0; noPartition < nbPartition; noPartition ++)
	{
	  curPartition = Lolimot->getPartitionSet()[noPartition];
	  
	  F_Lolimot << " " << VarName << " = " << VarName << " + (";
	  F_Lolimot << "((x(" << VarNb2 << ") - " << curPartition->getCentre(VarNb2 - 1) << ") / ";
	  F_Lolimot << curPartition->getEcartType(VarNb2 - 1) * curPartition->getEcartType(VarNb2 - 1) << ") * ";
	  F_Lolimot << "((x(" << VarNb1 << ") - " << curPartition->getCentre(VarNb1 - 1) << ") / ";
	  F_Lolimot << curPartition->getEcartType(VarNb1 - 1) * curPartition->getEcartType(VarNb1 - 1) << ")) * ";
	  F_Lolimot << "Ei(" << noPartition + 1 << ");" << endl;
	} /* End For */
    } /* End If */
  else
    {
      for(noPartition = 0; noPartition < nbPartition; noPartition ++)
	{
	  curPartition = Lolimot->getPartitionSet()[noPartition];
	  
	  F_Lolimot << " " << VarName << " = " << VarName << " - ";
	  F_Lolimot << "(1.0 / " << curPartition->getEcartType(VarNb1 - 1) * curPartition->getEcartType(VarNb1 - 1) << " - ";
	  F_Lolimot << "((x(" << VarNb1 << ") - " << curPartition->getCentre(VarNb1 - 1) << ") / ";
	  F_Lolimot << curPartition->getEcartType(VarNb1 - 1) * curPartition->getEcartType(VarNb1 - 1) << ")^2.0) * ";
	  F_Lolimot << "Ei(" << noPartition + 1 << ");" << endl;
	} /* End For */
    } /* End Else */
}

void export_ddU_x_Matlab_Text(LL_Lolimot * Lolimot, ofstream & F_Lolimot,
			      const string & ParamName1, const string & ParamName2, const string & VarName)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();

  F_Lolimot << endl;
  F_Lolimot << "    //------------------------------------------------------" << endl;
  F_Lolimot << "    // Computation of the second derivative of the numerator" << endl;
  F_Lolimot << "    //------------------------------------------------------" << endl;
  F_Lolimot << endl;

  F_Lolimot << " " << VarName << " = 0.0;" << endl;

  if (ParamName1!=ParamName2)
    {
      for(noPartition = 0; noPartition < nbPartition; noPartition ++)
	{
	  F_Lolimot << " " << VarName << " -= " << " (";
	  F_Lolimot << "Coeff(" << noPartition + 1<< "," << ParamName2 << ") * ";
	  F_Lolimot << "((x(" << ParamName1 << ") - Centre(" << noPartition + 1 << "," << ParamName1 << ")) / ";
	  F_Lolimot << "(EcartType(" << noPartition + 1 << "," << ParamName1 << ") * EcartType(";
	  F_Lolimot << noPartition + 1 << "," << ParamName1 << "))) + ";
	  F_Lolimot << "(Coeff(" << noPartition + 1 << "," << ParamName1 << ") - ";
	  F_Lolimot << "((x(" << ParamName1 << ") - Centre(" << noPartition + 1 << "," << ParamName1 << ")) / ";
	  F_Lolimot << "(EcartType(" << noPartition + 1 << "," << ParamName1 << ") * EcartType(";
	  F_Lolimot << noPartition + 1 << "," << ParamName1 << "))) * ";
	  F_Lolimot << "Li(" << noPartition + 1 << ")) * ";
	  F_Lolimot << "((x(" << ParamName2 << ") - Centre(" << noPartition + 1 << "," << ParamName2 << ")) / ";
	  F_Lolimot << "(EcartType(" << noPartition + 1 << "," << ParamName2 << ") * EcartType(";
	  F_Lolimot << noPartition + 1 << "," << ParamName2 << ")))) * ";
	  F_Lolimot << "Ei(" << noPartition + 1 << ");" << endl;
	} /* End For */
    } /* End If */
  else
    {
      for(noPartition = 0; noPartition < nbPartition; noPartition ++)
	{
	  F_Lolimot << " " << VarName << " -= (Li(" << noPartition + 1 << ") / ";
	  F_Lolimot << "(EcartType(" << noPartition + 1 << "," << ParamName1 << ") * EcartType(";
	  F_Lolimot << noPartition + 1 << "," << ParamName1 << ")) + 2 * ";
	  F_Lolimot << "Coeff(" << noPartition + 1 << "," << ParamName1 << ") * ";
	  F_Lolimot << "((x(" << ParamName1 << ") - Centre(" << noPartition + 1 << "," << ParamName1 << ")) / ";
	  F_Lolimot << "(EcartType(" << noPartition + 1 << "," << ParamName1 << ") * EcartType(";
	  F_Lolimot << noPartition + 1 << "," << ParamName1 << "))) - ";
	  F_Lolimot << "Li(" << noPartition + 1 << ") * ";
	  F_Lolimot << "pow((x(" << ParamName1 << ") - Centre(" << noPartition + 1 << "," << ParamName1 << ")) / ";
	  F_Lolimot << "(EcartType(" << noPartition + 1 << "," << ParamName1 << ") * EcartType(";
	  F_Lolimot << noPartition + 1 << "," << ParamName1 << ")), 2.0)) * ";
	  F_Lolimot << "Ei(" << noPartition + 1 << "));" << endl;
	} /* End For */
    } /* End Else */
}

void export_ddV_x_Matlab_Text(LL_Lolimot * Lolimot, ofstream & F_Lolimot, 
			      const string & ParamName1, const string & ParamName2, const string & VarName)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();

  F_Lolimot << endl;
  F_Lolimot << "    //--------------------------------------------------------" << endl;
  F_Lolimot << "    // Computation of the second derivative of the denominator" << endl;
  F_Lolimot << "    //--------------------------------------------------------" << endl;
  F_Lolimot << endl;

  F_Lolimot << " " << VarName << " = 0.0;" << endl;

  if (ParamName1!=ParamName2)
    {
      for(noPartition = 0; noPartition < nbPartition; noPartition ++)
	{
	  F_Lolimot << " " << VarName << " += (";
	  F_Lolimot << "((x(" << ParamName2 << ") - Centre(" << noPartition + 1 << "," << ParamName2 << ")) / ";
	  F_Lolimot << "(EcartType(" << noPartition + 1 << "," << ParamName2 << ") * EcartType(";
	  F_Lolimot << noPartition + 1 << "," << ParamName2 << "))) * ";
	  F_Lolimot << "((x(" << ParamName1 << ") - Centre(" << noPartition + 1 << "," << ParamName1 << ")) / ";
	  F_Lolimot << "(EcartType(" << noPartition + 1 << "," << ParamName1 << ") * EcartType(";
	  F_Lolimot << noPartition + 1 << "," << ParamName1 << ")))) * ";
	  F_Lolimot << "Ei(" << noPartition + 1 << ");" << endl;
	} /* End For */
    } /* End If */
  else
    {
      for(noPartition = 0; noPartition < nbPartition; noPartition ++)
	{
	  F_Lolimot << " " << VarName << " -= "; 
	  F_Lolimot << "(1.0 / (EcartType(" << noPartition + 1 << "," << ParamName1 << ") * EcartType(";
	  F_Lolimot << noPartition + 1 << "," << ParamName1 << ")) - ";
	  F_Lolimot << "pow((x(" << ParamName1 << ") - Centre(" << noPartition + 1 << "," << ParamName1 << ")) / ";
	  F_Lolimot << "(EcartType(" << noPartition + 1 << "," << ParamName1 << ") * EcartType(";
	  F_Lolimot << noPartition + 1 << "," << ParamName1 << ")), 2.0)) * ";
	  F_Lolimot << "Ei(" << noPartition + 1 << ");" << endl;
	} /* End For */
    } /* End Else */
}

bool exportFunctionInMatlab(LL_Lolimot * Lolimot, const string & path, const string & nameFunction, 
			    const vector<Model_Input> & ListInput)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;
  ofstream       F_Lolimot;
  stringstream   sTmp;
  string         strTmp;

  strTmp = path + "/" + nameFunction + ".m";
  F_Lolimot.open(strTmp.c_str());

  if (!F_Lolimot.is_open())
    {
      cerr << "** Error exportFunctionInMatlab() : can't create the file ";
      cerr << strTmp << "." << endl;
      return false;
    } /* End If */

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];
      curPartition->updateEcartType(Lolimot);
    } /* End For */

  F_Lolimot << "function res = " << nameFunction << "(x" << Lolimot->getCustomArguments() << ")" << endl;

  export_Init_Matlab(Lolimot, F_Lolimot);

  export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "U");
  export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "V");

  export_List_Partitions_Matlab(Lolimot, F_Lolimot);
  export_List_Variables_Matlab(Lolimot, F_Lolimot);

  if (Lolimot->getUseTransformLog())
    {
      for(unsigned int i=0; i<ListInput.size(); i++)
	{
	  if (ListInput[i].Type==RetroInput)
	    {
	      sTmp.str("");
	      sTmp << "x(" << i + 1 << ")";
	      Export_Direct_Log_Transform(F_Lolimot, sTmp.str(),
					  Lolimot->getTransformLogAlpha(),
					  Lolimot->getTransformLogMin(),
					  Lolimot->getTransformLogMax(),
					  Lolimot->getTransformLogEps());
	    } /* End If */
	} /* End For */
    } /* End For */

  export_Ei_x_Matlab(Lolimot, F_Lolimot);
  export_Li_x_Matlab(Lolimot, F_Lolimot);

  export_U_x_Matlab(Lolimot, F_Lolimot, "U");
  export_V_x_Matlab(Lolimot, F_Lolimot, "V");

  export_Destroy_Matlab(Lolimot, F_Lolimot);

  F_Lolimot << " res = " << Lolimot->getCustomReturn_Before() << "(U/V)" << Lolimot->getCustomReturn_After() << ";" << endl;

  if (Lolimot->getUseTransformLog()) Export_Inverse_Log_Transform(F_Lolimot, string("res"),
								  Lolimot->getTransformLogAlpha(),
								  Lolimot->getTransformLogMin(),
								  Lolimot->getTransformLogMax(),
								  Lolimot->getTransformLogEps());

  F_Lolimot.close();

  return true;
}

bool exportDerivativeFunctionInMatlab(LL_Lolimot * Lolimot, const string & path, const string & nameFunction,
				      unsigned int VarNb, const vector<Model_Input> & ListInput)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;
  ofstream       F_Lolimot;
  stringstream   sTmp;
  string         strTmp;

  strTmp = path + "/" + nameFunction + ".m";
  F_Lolimot.open(strTmp.c_str());

  if (!F_Lolimot.is_open())
    {
      cerr << "** Error exportDerivativeFunctionInMatlab() : can't create the file ";
      cerr << strTmp << "." << endl;
      return false;
    } /* End If */

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];
      curPartition->updateEcartType(Lolimot);
    } /* End For */

  F_Lolimot << "function res = " << nameFunction << "(x" << Lolimot->getCustomArguments() << ")" << endl;

  export_Init_Matlab(Lolimot, F_Lolimot);

  export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "U");
  export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "V");
  export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "dU");
  export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "dV");

  export_List_Partitions_Matlab(Lolimot, F_Lolimot);
  export_List_Variables_Matlab(Lolimot, F_Lolimot);

  if (Lolimot->getUseTransformLog())
    {
      for(unsigned int i=0; i<ListInput.size(); i++)
	{
	  if (ListInput[i].Type==RetroInput)
	    {
	      sTmp.str("");
	      sTmp << "x(" << i + 1 << ")";
	      Export_Direct_Log_Transform(F_Lolimot, sTmp.str(),
					  Lolimot->getTransformLogAlpha(),
					  Lolimot->getTransformLogMin(),
					  Lolimot->getTransformLogMax(),
					  Lolimot->getTransformLogEps());
	    } /* End If */
	} /* End For */
    } /* End For */

  export_Ei_x_Matlab(Lolimot, F_Lolimot);
  export_Li_x_Matlab(Lolimot, F_Lolimot);

  export_U_x_Matlab(Lolimot, F_Lolimot, "U");
  export_V_x_Matlab(Lolimot, F_Lolimot, "V");

  export_dV_x_Matlab(Lolimot, F_Lolimot, VarNb, "dV");
  export_dU_x_Matlab(Lolimot, F_Lolimot, VarNb, "dU");

  export_Destroy_Matlab(Lolimot, F_Lolimot);

  F_Lolimot << "    res = " << Lolimot->getCustomReturn_Before() << "((dU*V-dV*U)/(V*V))";
  F_Lolimot << Lolimot->getCustomReturn_After() << ";" << endl;

  if (Lolimot->getUseTransformLog()) Export_Inverse_Log_Transform_Deriv(F_Lolimot, "res",
									Lolimot->getTransformLogAlpha(),
									Lolimot->getTransformLogMin(),
									Lolimot->getTransformLogMax(),
									Lolimot->getTransformLogEps());

  F_Lolimot.close();
  
  return true;
}

bool exportDerivativeFunctionInMatlab2(LL_Lolimot * Lolimot, const string & path, const string & nameFunction,
				       const vector<Model_Input> & ListInput)
{
  unsigned int   noPartition = 0, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;
  ofstream       F_Lolimot;
  stringstream   sTmp;
  string         strTmp;

  strTmp = path + "/" + nameFunction + ".m";
  F_Lolimot.open(strTmp.c_str());

  if (!F_Lolimot.is_open())
    {
      cerr << "** Error exportDerivativeFunctionInMatlab2() : can't create the file ";
      cerr << strTmp << "." << endl;
      return false;
    } /* End If */

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];
      curPartition->updateEcartType(Lolimot);
    } /* End For */

  F_Lolimot << "function res = " << nameFunction << "(x, VarNb" << Lolimot->getCustomArguments() << ")" << endl;

  export_Init_Matlab(Lolimot, F_Lolimot);

  export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "U");
  export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "V");
  export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "dU");
  export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "dV");

  export_Init_Data_Matlab(Lolimot, F_Lolimot);

  export_List_Partitions_Matlab(Lolimot, F_Lolimot);
  export_List_Variables_Matlab(Lolimot, F_Lolimot);

  if (Lolimot->getUseTransformLog())
    {
      for(unsigned int i=0; i<ListInput.size(); i++)
	{
	  if (ListInput[i].Type==RetroInput)
	    {
	      sTmp.str("");
	      sTmp << "x(" << i + 1 << ")";
	      Export_Direct_Log_Transform(F_Lolimot, sTmp.str(),
					  Lolimot->getTransformLogAlpha(),
					  Lolimot->getTransformLogMin(),
					  Lolimot->getTransformLogMax(),
					  Lolimot->getTransformLogEps());
	    } /* End If */
	} /* End For */
    } /* End For */

  export_Ei_x_Matlab(Lolimot, F_Lolimot);
  export_Li_x_Matlab(Lolimot, F_Lolimot);

  export_U_x_Matlab(Lolimot, F_Lolimot, "U");
  export_V_x_Matlab(Lolimot, F_Lolimot, "V");

  export_dV_x_Matlab_Text(Lolimot, F_Lolimot, "VarNb", "dV");
  export_dU_x_Matlab_Text(Lolimot, F_Lolimot, "VarNb", "dU");

  export_Destroy_Data_Matlab(Lolimot, F_Lolimot);
  export_Destroy_Matlab(Lolimot, F_Lolimot);

  F_Lolimot << "    res = " << Lolimot->getCustomReturn_Before() << "((dU*V-dV*U)/(V*V))";
  F_Lolimot << Lolimot->getCustomReturn_After() << ";" << endl;

  if (Lolimot->getUseTransformLog()) Export_Inverse_Log_Transform_Deriv(F_Lolimot, "res",
									Lolimot->getTransformLogAlpha(),
									Lolimot->getTransformLogMin(),
									Lolimot->getTransformLogMax(),
									Lolimot->getTransformLogEps());

  F_Lolimot.close();
  
  return true;
}

bool exportSecondDerivativeFunctionInMatlab(LL_Lolimot * Lolimot, const string & path, const string & nameFunction, 
					    unsigned int VarNb1, unsigned int VarNb2, const vector<Model_Input> & ListInput)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;
  ofstream       F_Lolimot;
  stringstream   sTmp;
  string         strTmp;

  strTmp = path + "/" + nameFunction + ".m";
  F_Lolimot.open(strTmp.c_str());

  if (!F_Lolimot.is_open())
    {
      cerr << "** Error exportSecondDerivativeFunctionInMatlab() : can't create the file ";
      cerr << strTmp << "." << endl;
      return false;
    } /* End If */

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];
      curPartition->updateEcartType(Lolimot);
    } /* End For */

  F_Lolimot << "function res = " << nameFunction << "(x" << Lolimot->getCustomArguments() << ")" << endl;
  
  export_Init_Matlab(Lolimot, F_Lolimot);

  if (Lolimot->getUseTransformLog())
    {
      for(unsigned int i=0; i<ListInput.size(); i++)
	{
	  if (ListInput[i].Type==RetroInput)
	    {
	      sTmp.str("");
	      sTmp << "x(" << i + 1 << ")";
	      Export_Direct_Log_Transform(F_Lolimot, sTmp.str(), 
					  Lolimot->getTransformLogAlpha(),
					  Lolimot->getTransformLogMin(),
					  Lolimot->getTransformLogMax(),
					  Lolimot->getTransformLogEps());

	    } /* End If */
	} /* End For */
    } /* End For */

  if (VarNb1!=VarNb2)
    {
      export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "U");
      export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "V");
      export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "dU1");
      export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "dV1");
      export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "dU2");
      export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "dV2");
      export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "ddU");
      export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "ddV");
      
      export_List_Partitions_Matlab(Lolimot, F_Lolimot);
      export_List_Variables_Matlab(Lolimot, F_Lolimot);
      
      export_Ei_x_Matlab(Lolimot, F_Lolimot);
      export_Li_x_Matlab(Lolimot, F_Lolimot);
      
      export_U_x_Matlab(Lolimot, F_Lolimot, "U");
      export_V_x_Matlab(Lolimot, F_Lolimot, "V");
      
      export_dV_x_Matlab(Lolimot, F_Lolimot, VarNb1, "dV1");
      export_dU_x_Matlab(Lolimot, F_Lolimot, VarNb1, "dU1");
      export_dV_x_Matlab(Lolimot, F_Lolimot, VarNb2, "dV2");
      export_dU_x_Matlab(Lolimot, F_Lolimot, VarNb2, "dU2");
      
      export_ddV_x_Matlab(Lolimot, F_Lolimot, VarNb1, VarNb2, "ddV");
      export_ddU_x_Matlab(Lolimot, F_Lolimot, VarNb1, VarNb2, "ddU");
    } /* End If */
  else
    {
      export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "U");
      export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "V");
      export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "dU1");
      export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "dV1");
      export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "ddU");
      export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "ddV");
      
      export_List_Partitions_Matlab(Lolimot, F_Lolimot);
      export_List_Variables_Matlab(Lolimot, F_Lolimot);
      
      export_Ei_x_Matlab(Lolimot, F_Lolimot);
      export_Li_x_Matlab(Lolimot, F_Lolimot);
      
      export_U_x_Matlab(Lolimot, F_Lolimot, "U");
      export_V_x_Matlab(Lolimot, F_Lolimot, "V");
      
      export_dV_x_Matlab(Lolimot, F_Lolimot, VarNb1, "dV1");
      export_dU_x_Matlab(Lolimot, F_Lolimot, VarNb1, "dU1");
      
      export_ddV_x_Matlab(Lolimot, F_Lolimot, VarNb1, VarNb1, "ddV");
      export_ddU_x_Matlab(Lolimot, F_Lolimot, VarNb1, VarNb1, "ddU");
    } /* End Else */

  export_Destroy_Matlab(Lolimot, F_Lolimot);

  if (VarNb1!=VarNb2)
    {
      F_Lolimot << "    res =  " << Lolimot->getCustomReturn_Before();
      F_Lolimot << "(((ddU*V + dU1*dV2 - ddV*U - dV1*dU2)*V*V - (dU1*V-dV1*U)*2*V*dV2)/(V*V*V*V))";
      F_Lolimot << Lolimot->getCustomReturn_After() << ";" << endl;
    } /* End If */
  else
    {
      F_Lolimot << "    res = " << Lolimot->getCustomReturn_Before();
      F_Lolimot << "(((ddU*V - ddV*U)*V*V - (dU1*V-dV1*U)*2*V*dV1)/(V*V*V*V))";
      F_Lolimot << Lolimot->getCustomReturn_After() << ";" << endl;
    } /* End Else */

  if (Lolimot->getUseTransformLog()) Export_Inverse_Log_Transform_Deriv_2(F_Lolimot, "res",
									  Lolimot->getTransformLogAlpha(),
									  Lolimot->getTransformLogMin(),
									  Lolimot->getTransformLogMax(),
									  Lolimot->getTransformLogEps());
  
  F_Lolimot.close();
  
  return true;
}

bool exportSecondDerivativeFunctionInMatlab2(LL_Lolimot * Lolimot, const string & path, const string & nameFunction, 
					     const vector<Model_Input> & ListInput)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;
  ofstream       F_Lolimot;
  stringstream   sTmp;
  string         strTmp;

  strTmp = path + "/" + nameFunction + ".m";
  F_Lolimot.open(strTmp.c_str());

  if (!F_Lolimot.is_open())
    {
      cerr << "** Error exportSecondDerivativeFunctionInMatlab2() : can't create the file ";
      cerr << strTmp << "." << endl;
      return false;
    } /* End If */

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];
      curPartition->updateEcartType(Lolimot);
    } /* End For */

  F_Lolimot << "function res = " << nameFunction << "(x" << Lolimot->getCustomArguments() << ", VarNb1, VarNb2)" << endl;
  
  export_Init_Matlab(Lolimot, F_Lolimot);

  if (Lolimot->getUseTransformLog())
    {
      for(unsigned int i=0; i<ListInput.size(); i++)
	{
	  if (ListInput[i].Type==RetroInput)
	    {
	      sTmp.str("");
	      sTmp << "x(" << i + 1 << ")";
	      Export_Direct_Log_Transform(F_Lolimot, sTmp.str(), 
					  Lolimot->getTransformLogAlpha(),
					  Lolimot->getTransformLogMin(),
					  Lolimot->getTransformLogMax(),
					  Lolimot->getTransformLogEps());

	    } /* End If */
	} /* End For */
    } /* End For */

  export_Init_Data_Matlab(Lolimot, F_Lolimot);

  export_List_Partitions_Matlab(Lolimot, F_Lolimot);
  export_List_Variables_Matlab(Lolimot, F_Lolimot);

  export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "U");
  export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "V");
  export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "dU1");
  export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "dV1");
  export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "dU2");
  export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "dV2");
  export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "ddU");
  export_Declare_Variable_Matlab(Lolimot, F_Lolimot, "ddV");
  
  export_List_Partitions_Matlab(Lolimot, F_Lolimot);
  export_List_Variables_Matlab(Lolimot, F_Lolimot);
  
  export_Ei_x_Matlab(Lolimot, F_Lolimot);
  export_Li_x_Matlab(Lolimot, F_Lolimot);
  
  export_U_x_Matlab(Lolimot, F_Lolimot, "U");
  export_V_x_Matlab(Lolimot, F_Lolimot, "V");
  
  export_dV_x_Matlab_Text(Lolimot, F_Lolimot, "VarNb1", "dV1");
  export_dU_x_Matlab_Text(Lolimot, F_Lolimot, "VarNb1", "dU1");
  export_dV_x_Matlab_Text(Lolimot, F_Lolimot, "VarNb2", "dV2");
  export_dU_x_Matlab_Text(Lolimot, F_Lolimot, "VarNb2", "dU2");
  
  export_ddV_x_Matlab_Text(Lolimot, F_Lolimot, "VarNb1", "VarNb2", "ddV");
  export_ddU_x_Matlab_Text(Lolimot, F_Lolimot, "VarNb1", "VarNb2", "ddU");

  export_Destroy_Matlab(Lolimot, F_Lolimot);

  F_Lolimot << "    res =  " << Lolimot->getCustomReturn_Before();
  F_Lolimot << "(((ddU*V + dU1*dV2 - ddV*U - dV1*dU2)*V*V - (dU1*V-dV1*U)*2*V*dV2)/(V*V*V*V))";
  F_Lolimot << Lolimot->getCustomReturn_After() << ";" << endl;

  if (Lolimot->getUseTransformLog()) Export_Inverse_Log_Transform_Deriv_2(F_Lolimot, "res",
									  Lolimot->getTransformLogAlpha(),
									  Lolimot->getTransformLogMin(),
									  Lolimot->getTransformLogMax(),
									  Lolimot->getTransformLogEps());

  F_Lolimot.close();
  
  return true;
}
