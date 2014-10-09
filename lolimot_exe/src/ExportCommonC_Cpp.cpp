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

#include <fstream>
#include <string>
#include <vector>

#include "LL_Lolimot.h"
#include "LL_Partition.h"
#include "LL_Dimension.h"
#include "LL_Mesure.h"

using namespace std;

void export_Declare_Variable_C(LL_Lolimot * Lolimot, ofstream & F_Lolimot, const string & VarType, const string & VarName)
{
  F_Lolimot << VarType << " " << VarName << " = 0.0;" << endl;
}

void export_Ei_x_C(LL_Lolimot * Lolimot, ofstream & F_Lolimot, const string & VarType)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;
  unsigned int   noDimension, nbDimension = Lolimot->getDimensionSet().size();

  F_Lolimot << endl;
  F_Lolimot << "    //---------------------------------------------------------------" << endl;
  F_Lolimot << "    // For each partition i, we compute the membership function Ei(x)" << endl;
  F_Lolimot << "    //---------------------------------------------------------------" << endl;
  F_Lolimot << endl;

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];

      F_Lolimot << " Ei[" << noPartition << "] = exp(-0.5 * (";

      for(noDimension = 0; noDimension < nbDimension; noDimension ++)
	{
	  F_Lolimot << " + pow(x[" << noDimension << "] - " << curPartition->getCentre(noDimension) << ", 2.0) / ";
	  F_Lolimot << "(" << VarType << ")";
	  F_Lolimot << curPartition->getEcartType(noDimension) * curPartition->getEcartType(noDimension);
	} /* End For */

      F_Lolimot << "));" << endl;
    } /* End For */
}

void export_Li_x_C(LL_Lolimot * Lolimot, ofstream & F_Lolimot, const string & VarType)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;
  unsigned int   noDimension, nbDimension = Lolimot->getDimensionSet().size();

  F_Lolimot << endl;
  F_Lolimot << "    //--------------------------------------------------------" << endl;
  F_Lolimot << "    // For each partition i, we compute the linear model Li(x)" << endl;
  F_Lolimot << "    //--------------------------------------------------------" << endl;
  F_Lolimot << endl;

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];

      F_Lolimot << "    Li[" << noPartition << "] = " << curPartition->getCoeff0();

      for(noDimension = 0; noDimension < nbDimension; noDimension ++)
	{
	  F_Lolimot << " + " << curPartition->getCoeff(noDimension) << " * x[" << noDimension << "]";
	} /* End For */

      F_Lolimot << ";" << endl;
    } /* End For */
}

void export_U_x_C(LL_Lolimot * Lolimot, ofstream & F_Lolimot, const string & VarType, const string & VarName)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;

  F_Lolimot << endl;
  F_Lolimot << "    //-----------------------------" << endl;
  F_Lolimot << "    // Computation of the numerator" << endl;
  F_Lolimot << "    //-----------------------------" << endl;
  F_Lolimot << endl;

  F_Lolimot << " " << VarName << " = 0.0;" << endl;

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];

      F_Lolimot << " " << VarName << " += Li[" << noPartition << "] * Ei[" << noPartition << "];" << endl;
    } /* End For */
}

void export_V_x_C(LL_Lolimot * Lolimot, ofstream & F_Lolimot, const string & VarType, const string & VarName)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;

  F_Lolimot << endl;
  F_Lolimot << "    //-------------------------------" << endl;
  F_Lolimot << "    // Computation of the denominator" << endl;
  F_Lolimot << "    //-------------------------------" << endl;
  F_Lolimot << endl;

  F_Lolimot << " " << VarName << " = 0.0;" << endl;

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];

      F_Lolimot << " " << VarName << " += Ei[" << noPartition << "];" << endl;
    } /* End For */
}

void export_dU_x_C(LL_Lolimot * Lolimot, ofstream & F_Lolimot, unsigned int VarNb1, const string & VarType, const string & VarName)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;

  F_Lolimot << endl;
  F_Lolimot << "    //-----------------------------------------------" << endl;
  F_Lolimot << "    // Computation of the derivative of the numerator" << endl;
  F_Lolimot << "    //-----------------------------------------------" << endl;
  F_Lolimot << endl;

  F_Lolimot << " " << VarName << " = 0.0;" << endl;

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];

      F_Lolimot << " " << VarName << " += (" << curPartition->getCoeff(VarNb1) << " - ";
      F_Lolimot << "(x[" << VarNb1 << "] - " << curPartition->getCentre(VarNb1) << ") / ";
      F_Lolimot << "(" << VarType << ")";
      F_Lolimot << curPartition->getEcartType(VarNb1)*curPartition->getEcartType(VarNb1);
      F_Lolimot << " * Li[" << noPartition << "]) * Ei[" << noPartition << "];" << endl;
    } /* End For */
}

void export_dV_x_C(LL_Lolimot * Lolimot, ofstream & F_Lolimot, unsigned int VarNb1, const string & VarType, const string & VarName)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;

  F_Lolimot << endl;
  F_Lolimot << "    //-------------------------------------------------" << endl;
  F_Lolimot << "    // Computation of the derivative of the denominator" << endl;
  F_Lolimot << "    //-------------------------------------------------" << endl;
  F_Lolimot << endl;

  F_Lolimot << " " << VarName << " = 0.0;" << endl;

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];

      F_Lolimot << " " << VarName << " += - ((x[" << VarNb1 << "] - " << curPartition->getCentre(VarNb1) << ") / ";
      F_Lolimot << "(" << VarType << ")";
      F_Lolimot << curPartition->getEcartType(VarNb1)*curPartition->getEcartType(VarNb1) << ")";
      F_Lolimot << " * Ei[" << noPartition << "];" << endl;
    } /* End For */
}

void export_Data_C_Text(LL_Lolimot * Lolimot, ofstream & F_Lolimot, const string & VarType)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  unsigned int   noDimension, nbDimension = Lolimot->getDimensionSet().size();
  LL_Partition * curPartition = NULL;

  F_Lolimot << endl;
  F_Lolimot << "    //----------------------------" << endl;
  F_Lolimot << "    // List of the partitions data" << endl;
  F_Lolimot << "    //----------------------------" << endl;
  F_Lolimot << endl;

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];
      F_Lolimot << " Coeff0[" << noPartition << "] = " << curPartition->getCoeff0() << ";" << endl;
    } /* End For */

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];
      for(noDimension = 0; noDimension < nbDimension; noDimension ++)
	{
	  F_Lolimot << " Coeff[" << noPartition << "][" << noDimension << "] = ";
	  F_Lolimot << curPartition->getCoeff(noDimension) << ";" << endl;
	} /* End For */
    } /* End For */

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];
      for(noDimension = 0; noDimension < nbDimension; noDimension ++)
	{
	  F_Lolimot << " Centre[" << noPartition << "][" << noDimension << "] = ";
	  F_Lolimot << curPartition->getCentre(noDimension) << ";" << endl;
	} /* End For */
    } /* End For */

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];
      for(noDimension = 0; noDimension < nbDimension; noDimension ++)
	{
	  F_Lolimot << " EcartType[" << noPartition << "][" << noDimension << "] = ";
	  F_Lolimot << curPartition->getEcartType(noDimension) << ";" << endl;
	} /* End For */
    } /* End For */
}

void export_dV_x_C_Text(LL_Lolimot * Lolimot, ofstream & F_Lolimot, const string & ParamName, 
			const string & VarType, const string & VarName)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;

  F_Lolimot << endl;
  F_Lolimot << "    //-------------------------------------------------" << endl;
  F_Lolimot << "    // Computation of the derivative of the denominator" << endl;
  F_Lolimot << "    //-------------------------------------------------" << endl;
  F_Lolimot << endl;

  F_Lolimot << " " << VarName << " = 0.0;" << endl;

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];

      F_Lolimot << " " << VarName << " += - ((x[" << ParamName << "] - Centre[" << noPartition << "][";
      F_Lolimot << ParamName << "]) / ";
      F_Lolimot << "(" << VarType << ")";
      F_Lolimot << "(EcartType[" << noPartition << "][" << ParamName << "]*EcartType[" << noPartition << "][";
      F_Lolimot << ParamName << "])";
      F_Lolimot << " * Ei[" << noPartition << "];" << endl;
    } /* End For */
}

void export_dU_x_C_Text(LL_Lolimot * Lolimot, ofstream & F_Lolimot, const string & ParamName, 
			const string & VarType, const string & VarName)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;

  F_Lolimot << endl;
  F_Lolimot << "    //-----------------------------------------------" << endl;
  F_Lolimot << "    // Computation of the derivative of the numerator" << endl;
  F_Lolimot << "    //-----------------------------------------------" << endl;
  F_Lolimot << endl;

  F_Lolimot << " " << VarName << " = 0.0;" << endl;

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];

      F_Lolimot << " " << VarName << " += (Coeff[" << noPartition << "][" << ParamName << "] - ";
      F_Lolimot << "(x[" << ParamName << "] - Centre[" << noPartition << "][" << ParamName << "]) / ";
      F_Lolimot << "(" << VarType << ")";
      F_Lolimot << "(EcartType[" << noPartition << "][" << ParamName << "]*EcartType[";
      F_Lolimot << noPartition << "][" << ParamName << "])";
      F_Lolimot << " * Li[" << noPartition << "]) * Ei[" << noPartition << "];" << endl;
    } /* End For */
}

void export_List_Partitions_C(LL_Lolimot * Lolimot, ofstream & F_Lolimot)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;

  F_Lolimot << endl;
  F_Lolimot << "    //-------------------" << endl;
  F_Lolimot << "    // List of partitions" << endl;
  F_Lolimot << "    //-------------------" << endl;
  F_Lolimot << endl;

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];

      F_Lolimot << "    // Partition n " << noPartition+1 << " " << curPartition->printDimensions().c_str() << endl;
    } /* End For */
}

void export_List_Variables_C(LL_Lolimot * Lolimot, ofstream & F_Lolimot)
{
  unsigned int noDimension, nbDimension = Lolimot->getDimensionSet().size();

  F_Lolimot << endl;
  F_Lolimot << "    //------------------------------------------" << endl;
  F_Lolimot << "    // Liste of variables name and of boundaries" << endl;
  F_Lolimot << "    //------------------------------------------" << endl;
  F_Lolimot << endl;

  for(noDimension = 0; noDimension < nbDimension; noDimension ++)
    {
      F_Lolimot << "// ";
      F_Lolimot << Lolimot->getDimensionSet()[noDimension]->getName() << " - [";
      F_Lolimot << Lolimot->getDimensionSet()[noDimension]->getMin()  << "  ";
      F_Lolimot << Lolimot->getDimensionSet()[noDimension]->getMax()  << "]" << endl;
    } /* End For */

  F_Lolimot << "// ";
  F_Lolimot << "Measure - [" << Lolimot->getMesureSet()[0]->getMinValueExtrapolation();
  F_Lolimot << " " << Lolimot->getMesureSet()[0]->getMaxValueExtrapolation() << "]" << endl;
}

void export_ddU_x_C(LL_Lolimot * Lolimot, ofstream & F_Lolimot, unsigned int VarNb1, unsigned int VarNb2, 
		    const string & VarType, const string & VarName)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;

  F_Lolimot << endl;
  F_Lolimot << "    //------------------------------------------------------" << endl;
  F_Lolimot << "    // Computation of the second derivative of the numerator" << endl;
  F_Lolimot << "    //------------------------------------------------------" << endl;
  F_Lolimot << endl;

  F_Lolimot << " " << VarName << " = 0.0;" << endl;

  if (VarNb1!=VarNb2)
    {
      for(noPartition = 0; noPartition < nbPartition; noPartition ++)
	{
	  curPartition = Lolimot->getPartitionSet()[noPartition];
	  
	  F_Lolimot << " " << VarName << " -= " << " (";
	  F_Lolimot << curPartition->getCoeff(VarNb2) << " * ";
	  F_Lolimot << "((x[" << VarNb1 << "] - " << curPartition->getCentre(VarNb1) << ") / ";
	  F_Lolimot << curPartition->getEcartType(VarNb1) * curPartition->getEcartType(VarNb1) << ") + ";
	  F_Lolimot << "(" << curPartition->getCoeff(VarNb1) << " - ";
	  F_Lolimot << "((x[" << VarNb1 << "] - " << curPartition->getCentre(VarNb1) << ") / ";
	  F_Lolimot << curPartition->getEcartType(VarNb1) * curPartition->getEcartType(VarNb1) << ") * ";
	  F_Lolimot << "Li[" << noPartition  << "]) * ";
	  F_Lolimot << "((x[" << VarNb2 << "] - " << curPartition->getCentre(VarNb2) << ") / ";
	  F_Lolimot << curPartition->getEcartType(VarNb2) * curPartition->getEcartType(VarNb2) << ")) * ";
	  F_Lolimot << "Ei[" << noPartition << "];" << endl;
	} /* End For */
    } /* End If */
  else
    {
      for(noPartition = 0; noPartition < nbPartition; noPartition ++)
	{
	  curPartition = Lolimot->getPartitionSet()[noPartition];
	  
	  F_Lolimot << " " << VarName << " -= (Li[" << noPartition << "] / ";
	  F_Lolimot << curPartition->getEcartType(VarNb1)*curPartition->getEcartType(VarNb1) << " + 2 * ";
	  F_Lolimot << curPartition->getCoeff(VarNb1) << " * ";
	  F_Lolimot << "((x[" << VarNb1 << "] - " << curPartition->getCentre(VarNb1) << ") / ";
	  F_Lolimot << curPartition->getEcartType(VarNb1) * curPartition->getEcartType(VarNb1) << ") - ";
	  F_Lolimot << "Li[" << noPartition << "] * ";
	  F_Lolimot << "pow((x[" << VarNb1 << "] - " << curPartition->getCentre(VarNb1) << ") / ";
	  F_Lolimot << curPartition->getEcartType(VarNb1) * curPartition->getEcartType(VarNb1) << ", 2.0)) * ";
	  F_Lolimot << "Ei[" << noPartition << "];" << endl;
	} /* End For */
    } /* End Else */
}

void export_ddV_x_C(LL_Lolimot * Lolimot, ofstream & F_Lolimot, unsigned int VarNb1, unsigned int VarNb2, 
		    const string & VarType, const string & VarName)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;

  F_Lolimot << endl;
  F_Lolimot << "    //--------------------------------------------------------" << endl;
  F_Lolimot << "    // Computation of the second derivative of the denominator" << endl;
  F_Lolimot << "    //--------------------------------------------------------" << endl;
  F_Lolimot << endl;

  F_Lolimot << " " << VarName << " = 0.0;" << endl;

  if (VarNb1!=VarNb2)
    {
      for(noPartition = 0; noPartition < nbPartition; noPartition ++)
	{
	  curPartition = Lolimot->getPartitionSet()[noPartition];
	  
	  F_Lolimot << " " << VarName << " += (";
	  F_Lolimot << "((x[" << VarNb2 << "] - " << curPartition->getCentre(VarNb2) << ") / ";
	  F_Lolimot << curPartition->getEcartType(VarNb2) * curPartition->getEcartType(VarNb2) << ") * ";
	  F_Lolimot << "((x[" << VarNb1 << "] - " << curPartition->getCentre(VarNb1) << ") / ";
	  F_Lolimot << curPartition->getEcartType(VarNb1) * curPartition->getEcartType(VarNb1) << ")) * ";
	  F_Lolimot << "Ei[" << noPartition << "];" << endl;
	} /* End For */
    } /* End If */
  else
    {
      for(noPartition = 0; noPartition < nbPartition; noPartition ++)
	{
	  curPartition = Lolimot->getPartitionSet()[noPartition];
	  
	  F_Lolimot << " " << VarName << " -= "; 
	  F_Lolimot << "(1.0 / " << curPartition->getEcartType(VarNb1) * curPartition->getEcartType(VarNb1) << " - ";
	  F_Lolimot << "pow((x[" << VarNb1 << "] - " << curPartition->getCentre(VarNb1) << ") / ";
	  F_Lolimot << curPartition->getEcartType(VarNb1) * curPartition->getEcartType(VarNb1) << ", 2.0)) * ";
	  F_Lolimot << "Ei[" << noPartition << "];" << endl;
	} /* End For */
    } /* End Else */
}

void export_ddU_x_C_Text(LL_Lolimot * Lolimot, ofstream & F_Lolimot, const string & ParamName1, 
			 const string & ParamName2, const string & VarType, const string & VarName)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;

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
	  curPartition = Lolimot->getPartitionSet()[noPartition];
	  
	  F_Lolimot << " " << VarName << " -= " << " (";
	  F_Lolimot << "Coeff[" << noPartition << "][" << ParamName2 << "] * ";
	  F_Lolimot << "((x[" << ParamName1 << "] - Centre[" << noPartition << "][" << ParamName1 << "]) / ";
	  F_Lolimot << "(EcartType[" << noPartition << "][" << ParamName1 << "] * EcartType[";
	  F_Lolimot << noPartition << "][" << ParamName1 << "])) + ";
	  F_Lolimot << "(Coeff[" << noPartition << "][" << ParamName1 << "] - ";
	  F_Lolimot << "((x[" << ParamName1 << "] - Centre[" << noPartition << "][" << ParamName1 << "]) / ";
	  F_Lolimot << "(EcartType[" << noPartition << "][" << ParamName1 << "] * EcartType[";
	  F_Lolimot << noPartition << "][" << ParamName1 << "])) * ";
	  F_Lolimot << "Li[" << noPartition  << "]) * ";
	  F_Lolimot << "((x[" << ParamName2 << "] - Centre[" << noPartition << "][" << ParamName2 << "]) / ";
	  F_Lolimot << "(EcartType[" << noPartition << "][" << ParamName2 << "] * EcartType[";
	  F_Lolimot << noPartition << "][" << ParamName2 << "]))) * ";
	  F_Lolimot << "Ei[" << noPartition << "];" << endl;
	} /* End For */
    } /* End If */
  else
    {
      for(noPartition = 0; noPartition < nbPartition; noPartition ++)
	{
	  curPartition = Lolimot->getPartitionSet()[noPartition];
	  
	  F_Lolimot << " " << VarName << " -= (Li[" << noPartition << "] / ";
	  F_Lolimot << "(EcartType[" << noPartition << "][" << ParamName1 << "] * EcartType[";
	  F_Lolimot << noPartition << "][" << ParamName1 << "]) + 2 * ";
	  F_Lolimot << "Coeff[" << noPartition << "][" << ParamName1 << "] * ";
	  F_Lolimot << "((x[" << ParamName1 << "] - Centre[" << noPartition << "][" << ParamName1 << "]) / ";
	  F_Lolimot << "(EcartType[" << noPartition << "][" << ParamName1 << "] * EcartType[";
	  F_Lolimot << noPartition << "][" << ParamName1 << "])) - ";
	  F_Lolimot << "Li[" << noPartition << "] * ";
	  F_Lolimot << "pow((x[" << ParamName1 << "] - Centre[" << noPartition << "][" << ParamName1 << "]) / ";
	  F_Lolimot << "(EcartType[" << noPartition << "][" << ParamName1 << "] * EcartType[";
	  F_Lolimot << noPartition << "][" << ParamName1 << "]), 2.0)) * ";
	  F_Lolimot << "Ei[" << noPartition << "]);" << endl;
	} /* End For */
    } /* End Else */
}

void export_ddV_x_C_Text(LL_Lolimot * Lolimot, ofstream & F_Lolimot, const string & ParamName1, 
			 const string & ParamName2, const string & VarType, const string & VarName)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;

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
	  curPartition = Lolimot->getPartitionSet()[noPartition];
	  
	  F_Lolimot << " " << VarName << " += (";
	  F_Lolimot << "((x[" << ParamName2 << "] - Centre[" << noPartition << "][" << ParamName2 << "]) / ";
	  F_Lolimot << "(EcartType[" << noPartition << "][" << ParamName2 << "] * EcartType[";
	  F_Lolimot << noPartition << "][" << ParamName2 << "])) * ";
	  F_Lolimot << "((x[" << ParamName1 << "] - Centre[" << noPartition << "][" << ParamName1 << "]) / ";
	  F_Lolimot << "(EcartType[" << noPartition << "][" << ParamName1 << "] * EcartType[";
	  F_Lolimot << noPartition << "][" << ParamName1 << "]))) * ";
	  F_Lolimot << "Ei[" << noPartition << "];" << endl;
	} /* End For */
    } /* End If */
  else
    {
      for(noPartition = 0; noPartition < nbPartition; noPartition ++)
	{
	  curPartition = Lolimot->getPartitionSet()[noPartition];
	  
	  F_Lolimot << " " << VarName << " -= "; 
	  F_Lolimot << "(1.0 / (EcartType[" << noPartition << "][" << ParamName1 << "] * EcartType[";
	  F_Lolimot << noPartition << "][" << ParamName1 << "]) - ";
	  F_Lolimot << "pow((x[" << ParamName1 << "] - Centre[" << noPartition << "][" << ParamName1 << "]) / ";
	  F_Lolimot << "(EcartType[" << noPartition << "][" << ParamName1 << "] * EcartType[";
	  F_Lolimot << noPartition << "][" << ParamName1 << "]), 2.0)) * ";
	  F_Lolimot << "Ei[" << noPartition << "];" << endl;
	} /* End For */
    } /* End Else */
}
