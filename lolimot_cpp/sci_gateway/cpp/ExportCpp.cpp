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

#include <ExportCommonC_Cpp.h>
#include <LL_Lolimot.h>
#include <LL_Partition.h>
#include <Post.h>
#include <TrainLolimotStruct.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

void export_Init_Cpp(LL_Lolimot * Lolimot, ofstream & F_Lolimot, const string & VarType)
{
  unsigned int  nbPartition = Lolimot->getPartitionSet().size();

  F_Lolimot << "vector<" << VarType << "> Ei(" << nbPartition << ");" << endl;
  F_Lolimot << "vector<" << VarType << "> Li(" << nbPartition << ");" << endl;
}

void export_Destroy_Cpp(LL_Lolimot * Lolimot, ofstream & F_Lolimot, const string & VarType)
{
  F_Lolimot << "  Ei.clear();" << endl;
  F_Lolimot << "  Li.clear();" << endl;
}

void export_Init_Data_Cpp(LL_Lolimot * Lolimot, ofstream & F_Lolimot, const string & VarType)
{
  unsigned int nbPartition = Lolimot->getPartitionSet().size();
  unsigned int nbDimension = Lolimot->getDimensionSet().size();

  F_Lolimot << "vector<vector<" << VarType << "> > Coeff;" << endl;
  F_Lolimot << "vector<vector<" << VarType << "> > Centre;" << endl;
  F_Lolimot << "vector<vector<" << VarType << "> > EcartType;" << endl;
  F_Lolimot << "vector<" << VarType << "> Coeff0;" << endl;
  F_Lolimot << "unsigned int i = 0;" << endl;
  
  F_Lolimot << " Coeff.resize(" << nbPartition << ", vector<" << VarType << ">(" << nbDimension << ", 0.0));" << endl;
  F_Lolimot << " Centre.resize(" << nbPartition << ", vector<" << VarType << ">(" << nbDimension << ", 0.0));" << endl;
  F_Lolimot << " EcartType.resize(" << nbPartition << ", vector<" << VarType << ">(" << nbDimension << ", 0.0));" << endl;
  F_Lolimot << " Coeff0.resize(" << nbPartition << ", 0.0);" << endl;
}

void export_Destroy_Data_Cpp(LL_Lolimot * Lolimot, 
			     ofstream & F_Lolimot, 
			     const string & VarType)
{
  F_Lolimot << "Coeff.clear();" << endl;
  F_Lolimot << "Centre.clear();" << endl;
  F_Lolimot << "EcartType.clear();" << endl;
  F_Lolimot << "Coeff0.clear();" << endl;
}

bool exportHeaderInCpp2(LL_Lolimot * Lolimot,
			const string & nameFunction, 
			const vector<Model_Input> & ListInput)
{
  ofstream       F_Lolimot;
  string         strTmp;
  string         VarType;

  strTmp = nameFunction + ".h";
  F_Lolimot.open(strTmp.c_str());

  if (!F_Lolimot.is_open())
    {
      cerr << "** Error exportHeaderInC() : can't create the file ";
      cerr << strTmp << "." << endl;
      return false;
    } /* End If */

  if (Lolimot->getExportType_double())
    {
      VarType = "double";
    } /* End If */
  else
    {
      VarType = "float";
    } /* End Else */

  F_Lolimot << "extern \"C\" ";
  F_Lolimot << VarType << " " << nameFunction << "(" << VarType << " * x, unsigned int VarNb";
  F_Lolimot << Lolimot->getCustomArguments() << ");" << endl;

  F_Lolimot.close();

  return true;
}

bool exportHeaderInCpp(LL_Lolimot * Lolimot,
		       const string & nameFunction, 
		       const vector<Model_Input> & ListInput)
{
  ofstream       F_Lolimot;
  string         strTmp;
  string         VarType;

  strTmp = nameFunction + ".h";
  F_Lolimot.open(strTmp.c_str());

  if (!F_Lolimot.is_open())
    {
      cerr << "** Error exportHeaderInCpp() : can't create the file ";
      cerr << strTmp << "." << endl;
      return false;
    } /* End If */

  if (Lolimot->getExportType_double())
    {
      VarType = "double";
    } /* End If */
  else
    {
      VarType = "float";
    } /* End Else */

  F_Lolimot << "extern \"C\" ";
  F_Lolimot << VarType << " " << nameFunction << "(" << VarType << " * x" << Lolimot->getCustomArguments() << ");" << endl;

  F_Lolimot.close();

  return true;
}

bool exportFunctionInCpp(LL_Lolimot * Lolimot, 
			 const string & nameFunction, 
			 const vector<Model_Input> & ListInput)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;
  ofstream       F_Lolimot;
  stringstream   sTmp;
  string         strTmp;
  string         VarType;

  strTmp = nameFunction + ".cpp";
  F_Lolimot.open(strTmp.c_str());

  if (!F_Lolimot.is_open())
    {
      cerr << "** Error exportFunctionInCpp() : can't create the file ";
      cerr << strTmp << "." << endl;
      return false;
    } /* End If */

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];
      curPartition->updateEcartType(Lolimot);
    } /* End For */

  if (Lolimot->getExportType_double())
    {
      VarType = "double";
    } /* End If */
  else
    {
      VarType = "float";
    } /* End Else */

  F_Lolimot << "#include <vector>" << endl;
  F_Lolimot << "#include <cmath>" << endl << endl;
  F_Lolimot << endl << "using namespace std;" << endl;
  F_Lolimot << "extern \"C\" ";
  F_Lolimot << VarType << " " << nameFunction << "(vector<" << VarType << "> & x";
  F_Lolimot << Lolimot->getCustomArguments() << ")" << endl;
  F_Lolimot << "{" << endl;
  F_Lolimot << VarType << " Result;" << endl;

  export_Init_Cpp(Lolimot, F_Lolimot, VarType);

  export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "U");
  export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "V");

  export_List_Partitions_C(Lolimot, F_Lolimot);
  export_List_Variables_C(Lolimot, F_Lolimot);

  if (Lolimot->getUseTransformLog())
    {
      for(unsigned int i=0; i<ListInput.size(); i++)
	{
	  if (ListInput[i].Type==RetroInput)
	    {
	      sTmp.str("");
	      sTmp << "x[" << i << "]";
	      Export_Direct_Log_Transform(F_Lolimot, sTmp.str(),
					  Lolimot->getTransformLogAlpha(),
					  Lolimot->getTransformLogMin(),
					  Lolimot->getTransformLogMax(),
					  Lolimot->getTransformLogEps());
	    } /* End If */
	} /* End For */
    } /* End For */

  export_Ei_x_C(Lolimot, F_Lolimot, VarType);
  export_Li_x_C(Lolimot, F_Lolimot, VarType);

  export_U_x_C(Lolimot, F_Lolimot, VarType, "U");
  export_V_x_C(Lolimot, F_Lolimot, VarType, "V");

  export_Destroy_Cpp(Lolimot, F_Lolimot, VarType);

  F_Lolimot << " Result = " << Lolimot->getCustomReturn_Before() << "(U/(" << VarType << ")V)";
  F_Lolimot << Lolimot->getCustomReturn_After() << ";" << endl;

  if (Lolimot->getUseTransformLog()) Export_Inverse_Log_Transform(F_Lolimot, "Result",
								  Lolimot->getTransformLogAlpha(),
								  Lolimot->getTransformLogMin(),
								  Lolimot->getTransformLogMax(),
								  Lolimot->getTransformLogEps());

  F_Lolimot << " return Result;" << endl;
  F_Lolimot << "}" << endl;

  F_Lolimot.close();

  return true;
}

bool exportDerivativeFunctionInCpp(LL_Lolimot * Lolimot, 
				   const string & nameFunction, 
				   unsigned int VarNb, 
				   const vector<Model_Input> & ListInput)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;
  ofstream       F_Lolimot;
  stringstream   sTmp;
  string         strTmp;
  string         VarType;

  strTmp = nameFunction + ".cpp";
  F_Lolimot.open(strTmp.c_str());

  if (!F_Lolimot.is_open())
    {
      cerr << "** Error exportDerivativeFunctionInCpp() : can't create the file ";
      cerr << strTmp << "." << endl;
      return false;
    } /* End If */

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];
      curPartition->updateEcartType(Lolimot);
    } /* End For */

  if (Lolimot->getExportType_double())
    {
      VarType = "double";
    } /* End If */
  else
    {
      VarType = "float";
    } /* End Else */

  F_Lolimot << "#include <vector>" << endl;
  F_Lolimot << "#include <cmath>" << endl << endl;
  F_Lolimot << endl << "using namespace std;" << endl;
  F_Lolimot << "extern \"C\" ";
  F_Lolimot << VarType << " " << nameFunction << "(vector<" << VarType << "> & x";
  F_Lolimot << Lolimot->getCustomArguments() << ")" << endl;
  F_Lolimot << "{" << endl;
  F_Lolimot << VarType << " Result;" << endl;

  export_Init_Cpp(Lolimot, F_Lolimot, VarType);

  export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "U");
  export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "V");
  export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "dU");
  export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "dV");

  export_List_Partitions_C(Lolimot, F_Lolimot);
  export_List_Variables_C(Lolimot, F_Lolimot);

  if (Lolimot->getUseTransformLog())
    {
      for(unsigned int i=0; i<ListInput.size(); i++)
	{
	  if (ListInput[i].Type==RetroInput)
	    {
	      sTmp.str("");
	      sTmp << "x[" << i << "]";
	      Export_Direct_Log_Transform(F_Lolimot, sTmp.str(),
					  Lolimot->getTransformLogAlpha(),
					  Lolimot->getTransformLogMin(),
					  Lolimot->getTransformLogMax(),
					  Lolimot->getTransformLogEps());
	    } /* End If */
	} /* End For */
    } /* End For */

  export_Ei_x_C(Lolimot, F_Lolimot, VarType);
  export_Li_x_C(Lolimot, F_Lolimot, VarType);

  export_U_x_C(Lolimot, F_Lolimot, VarType, "U");
  export_V_x_C(Lolimot, F_Lolimot, VarType, "V");

  export_dV_x_C(Lolimot, F_Lolimot, VarNb, VarType, "dV");
  export_dU_x_C(Lolimot, F_Lolimot, VarNb, VarType, "dU");

  export_Destroy_Cpp(Lolimot, F_Lolimot, VarType);

  F_Lolimot << "    Result = " << Lolimot->getCustomReturn_Before() << "((dU*V-dV*U)/(" << VarType << ")(V*V))";
  F_Lolimot << Lolimot->getCustomReturn_After() << ";" << endl;
  
  if (Lolimot->getUseTransformLog()) Export_Inverse_Log_Transform_Deriv(F_Lolimot, "Result",
									Lolimot->getTransformLogAlpha(),
									Lolimot->getTransformLogMin(),
									Lolimot->getTransformLogMax(),
									Lolimot->getTransformLogEps());


  F_Lolimot << " return Result;" << endl;

  F_Lolimot << "}" << endl;
  F_Lolimot.close();

  return true;
}

bool exportDerivativeFunctionInCpp2(LL_Lolimot * Lolimot, 
				    const string & nameFunction,
				    const vector<Model_Input> & ListInput)
{
  unsigned int   noPartition = 0, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;
  ofstream       F_Lolimot;
  stringstream   sTmp;
  string         strTmp;
  string         VarType;

  strTmp = nameFunction + ".cpp";
  F_Lolimot.open(strTmp.c_str());

  if (!F_Lolimot.is_open())
    {
      cerr << "** Error exportDerivativeFunctionInCpp2() : can't create the file ";
      cerr << strTmp << "." << endl;
      return false;
    } /* End If */

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];
      curPartition->updateEcartType(Lolimot);
    } /* End For */

  if (Lolimot->getExportType_double())
    {
      VarType = "double";
    } /* End If */
  else
    {
      VarType = "float";
    } /* End Else */

  F_Lolimot << "#include <vector>" << endl;
  F_Lolimot << "#include <cmath>" << endl << endl;
  F_Lolimot << endl << "using namespace std;" << endl;
  F_Lolimot << "extern \"C\" ";
  F_Lolimot << VarType << " " << nameFunction << "(vector<" << VarType << "> & x, unsigned int VarNb";
  F_Lolimot << Lolimot->getCustomArguments() << ")" << endl;
  F_Lolimot << "{" << endl;
  F_Lolimot << VarType << " Result;" << endl;

  export_Init_Cpp(Lolimot, F_Lolimot, VarType);

  export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "U");
  export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "V");
  export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "dU");
  export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "dV");

  export_Init_Data_Cpp(Lolimot, F_Lolimot, VarType);

  export_List_Partitions_C(Lolimot, F_Lolimot);
  export_List_Variables_C(Lolimot, F_Lolimot);

  if (Lolimot->getUseTransformLog())
    {
      for(unsigned int i=0; i<ListInput.size(); i++)
	{
	  if (ListInput[i].Type==RetroInput)
	    {
	      sTmp.str("");
	      sTmp << "x[" << i << "]";
	      Export_Direct_Log_Transform(F_Lolimot, sTmp.str(), 
					  Lolimot->getTransformLogAlpha(),
					  Lolimot->getTransformLogMin(),
					  Lolimot->getTransformLogMax(),
					  Lolimot->getTransformLogEps());
	    } /* End If */
	} /* End For */
    } /* End For */

  export_Ei_x_C(Lolimot, F_Lolimot, VarType);
  export_Li_x_C(Lolimot, F_Lolimot, VarType);

  export_U_x_C(Lolimot, F_Lolimot, VarType, "U");
  export_V_x_C(Lolimot, F_Lolimot, VarType, "V");

  export_dV_x_C_Text(Lolimot, F_Lolimot, "VarNb", VarType, "dV");
  export_dU_x_C_Text(Lolimot, F_Lolimot, "VarNb", VarType, "dU");

  export_Destroy_Data_Cpp(Lolimot, F_Lolimot, VarType);
  export_Destroy_Cpp(Lolimot, F_Lolimot, VarType);

  F_Lolimot << "    Result = " << Lolimot->getCustomReturn_Before() << "((dU*V-dV*U)/(" << VarType << ")(V*V))";
  F_Lolimot << Lolimot->getCustomReturn_After() << ";" << endl;

  if (Lolimot->getUseTransformLog()) Export_Inverse_Log_Transform_Deriv(F_Lolimot, "Result",
									Lolimot->getTransformLogAlpha(),
									Lolimot->getTransformLogMin(),
									Lolimot->getTransformLogMax(),
									Lolimot->getTransformLogEps());

  F_Lolimot << " return Result;" << endl;
  F_Lolimot << "}" << endl;

  F_Lolimot.close();

  return true;
}

bool exportSecondDerivativeFunctionInCpp(LL_Lolimot * Lolimot, 
					 const string & nameFunction,
					 unsigned int VarNb1, 
					 unsigned int VarNb2, 
					 const vector<Model_Input> & ListInput)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;
  ofstream       F_Lolimot;
  stringstream   sTmp;
  string         strTmp;
  string         VarType;

  strTmp = nameFunction + ".cpp";
  F_Lolimot.open(strTmp.c_str());

  if (!F_Lolimot.is_open())
    {
      cerr << "** Error exportSecondDerivativeFunctionInCpp() : can't create the file ";
      cerr << strTmp << "." << endl;
      return false;
    } /* End If */

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];
      curPartition->updateEcartType(Lolimot);
    } /* End For */

  if (Lolimot->getExportType_double())
    {
      VarType = "double";
    } /* End If */
  else
    {
      VarType = "float";
    } /* End Else */

  F_Lolimot << "#include <vector>" << endl;
  F_Lolimot << "#include <cmath>" << endl << endl;
  F_Lolimot << endl << "using namespace std;" << endl;
  F_Lolimot << "extern \"C\" ";
  F_Lolimot << VarType << " " << nameFunction << "(vector<" << VarType << "> & x";
  F_Lolimot << Lolimot->getCustomArguments() << ")" << endl;
  F_Lolimot << "{" << endl;
  F_Lolimot << " " << VarType << " Result = 0.0;" << endl << endl;

  export_Init_Cpp(Lolimot, F_Lolimot, VarType);

  if (Lolimot->getUseTransformLog())
    {
      for(unsigned int i=0; i<ListInput.size(); i++)
	{
	  if (ListInput[i].Type==RetroInput)
	    {
	      sTmp.str("");
	      sTmp << "x[" << i << "]";
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
      export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "U");
      export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "V");
      export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "dU1");
      export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "dV1");
      export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "dU2");
      export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "dv2");
      export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "ddU");
      export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "ddV");
      
      export_List_Partitions_C(Lolimot, F_Lolimot);
      export_List_Variables_C(Lolimot, F_Lolimot);
      
      export_Ei_x_C(Lolimot, F_Lolimot, VarType);
      export_Li_x_C(Lolimot, F_Lolimot, VarType);
      
      export_U_x_C(Lolimot, F_Lolimot, VarType, "U");
      export_V_x_C(Lolimot, F_Lolimot, VarType, "V");
      
      export_dV_x_C(Lolimot, F_Lolimot, VarNb1, VarType, "dV1");
      export_dU_x_C(Lolimot, F_Lolimot, VarNb1, VarType, "dU1");
      export_dV_x_C(Lolimot, F_Lolimot, VarNb2, VarType, "dV2");
      export_dU_x_C(Lolimot, F_Lolimot, VarNb2, VarType, "dU2");
      
      export_ddV_x_C(Lolimot, F_Lolimot, VarNb1, VarNb2, VarType, "ddV");
      export_ddU_x_C(Lolimot, F_Lolimot, VarNb1, VarNb2, VarType, "ddU");
    } /* End If */
  else
    {
      export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "U");
      export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "V");
      export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "dU1");
      export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "dV1");
      export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "ddU");
      export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "ddV");
      
      export_List_Partitions_C(Lolimot, F_Lolimot);
      export_List_Variables_C(Lolimot, F_Lolimot);
      
      export_Ei_x_C(Lolimot, F_Lolimot, VarType);
      export_Li_x_C(Lolimot, F_Lolimot, VarType);
      
      export_U_x_C(Lolimot, F_Lolimot, VarType, "U");
      export_V_x_C(Lolimot, F_Lolimot, VarType, "V");
      
      export_dV_x_C(Lolimot, F_Lolimot, VarNb1, VarType, "dV1");
      export_dU_x_C(Lolimot, F_Lolimot, VarNb1, VarType, "dU1");
      
      export_ddV_x_C(Lolimot, F_Lolimot, VarNb1, VarNb1, VarType, "ddV");
      export_ddU_x_C(Lolimot, F_Lolimot, VarNb1, VarNb1, VarType, "ddU");
    } /* End Else */

  export_Destroy_Cpp(Lolimot, F_Lolimot, VarType);

  if (VarNb1!=VarNb2)
    {
      F_Lolimot << "    Result = " << Lolimot->getCustomReturn_Before();
      F_Lolimot << "(((ddU*V + dU1*dV2 - ddV*U - dV1*dU2)*V*V - (dU1*V-dV1*U)*2*V*dV2)/(" << VarType << ")(V*V*V*V))";
      F_Lolimot << Lolimot->getCustomReturn_After() << ";" << endl;
    } /* End If */
  else
    {
      F_Lolimot << "    Result = " << Lolimot->getCustomReturn_Before();
      F_Lolimot << "(((ddU*V - ddV*U)*V*V - (dU1*V-dV1*U)*2*V*dV1)/(" << VarType << ")(V*V*V*V))";
      F_Lolimot << Lolimot->getCustomReturn_After() << ";" << endl;
    } /* End Else */

  if (Lolimot->getUseTransformLog()) Export_Inverse_Log_Transform_Deriv_2(F_Lolimot, "Result",
									  Lolimot->getTransformLogAlpha(),
									  Lolimot->getTransformLogMin(),
									  Lolimot->getTransformLogMax(),
									  Lolimot->getTransformLogEps());
  
  F_Lolimot << " return Result;" << endl;
  F_Lolimot << "}" << endl;

  F_Lolimot.close();

  return true;
}

bool exportSecondDerivativeFunctionInCpp2(LL_Lolimot * Lolimot, 
					  const string & nameFunction,
					  const vector<Model_Input> & ListInput)
{
  unsigned int   noPartition, nbPartition = Lolimot->getPartitionSet().size();
  LL_Partition * curPartition = NULL;
  ofstream       F_Lolimot;
  stringstream   sTmp;
  string         strTmp;
  string         VarType;

  strTmp = nameFunction + ".cpp";
  F_Lolimot.open(strTmp.c_str());

  if (!F_Lolimot.is_open())
    {
      cerr << "** Error exportSecondDerivativeFunctionInCpp2() : can't create the file ";
      cerr << strTmp << "." << endl;
      return false;
    } /* End If */

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = Lolimot->getPartitionSet()[noPartition];
      curPartition->updateEcartType(Lolimot);
    } /* End For */

  if (Lolimot->getExportType_double())
    {
      VarType = "double";
    } /* End If */
  else
    {
      VarType = "float";
    } /* End Else */

  F_Lolimot << "#include <vector>" << endl;
  F_Lolimot << "#include <cmath>" << endl << endl;
  F_Lolimot << endl << "using namespace std;" << endl;
  F_Lolimot << "extern \"C\" ";
  F_Lolimot << VarType << " " << nameFunction << "(vector<" << VarType << "> & x, unsigned int VarNb1, unsigned int VarNb2";
  F_Lolimot << Lolimot->getCustomArguments() << ")" << endl;
  F_Lolimot << "{" << endl;
  F_Lolimot << " " << VarType << " Result = 0.0;" << endl << endl;

  export_Init_Cpp(Lolimot, F_Lolimot, VarType);

  if (Lolimot->getUseTransformLog())
    {
      for(unsigned int i=0; i<ListInput.size(); i++)
	{
	  if (ListInput[i].Type==RetroInput)
	    {
	      sTmp.str("");
	      sTmp << "x[" << i << "]";
	      Export_Direct_Log_Transform(F_Lolimot, sTmp.str(),
					  Lolimot->getTransformLogAlpha(),
					  Lolimot->getTransformLogMin(),
					  Lolimot->getTransformLogMax(),
					  Lolimot->getTransformLogEps());
	    } /* End If */
	} /* End For */
    } /* End For */

  export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "U");
  export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "V");
  export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "dU1");
  export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "dV1");
  export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "dU2");
  export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "dv2");
  export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "ddU");
  export_Declare_Variable_C(Lolimot, F_Lolimot, VarType, "ddV");
  
  export_List_Partitions_C(Lolimot, F_Lolimot);
  export_List_Variables_C(Lolimot, F_Lolimot);
  
  export_Ei_x_C(Lolimot, F_Lolimot, VarType);
  export_Li_x_C(Lolimot, F_Lolimot, VarType);
  
  export_U_x_C(Lolimot, F_Lolimot, VarType, "U");
  export_V_x_C(Lolimot, F_Lolimot, VarType, "V");
  
  export_dV_x_C_Text(Lolimot, F_Lolimot, "VarNb1", VarType, "dV1");
  export_dU_x_C_Text(Lolimot, F_Lolimot, "VarNb1", VarType, "dU1");
  export_dV_x_C_Text(Lolimot, F_Lolimot, "VarNb2", VarType, "dV2");
  export_dU_x_C_Text(Lolimot, F_Lolimot, "VarNb2", VarType, "dU2");
      
  export_ddV_x_C_Text(Lolimot, F_Lolimot, "VarNb1", "VarNb2", VarType, "ddV");
  export_ddU_x_C_Text(Lolimot, F_Lolimot, "VarNb1", "VarNb2", VarType, "ddU");

  export_Destroy_Cpp(Lolimot, F_Lolimot, VarType);

  F_Lolimot << "    Result = " << Lolimot->getCustomReturn_Before();
  F_Lolimot << "(((ddU*V + dU1*dV2 - ddV*U - dV1*dU2)*V*V - (dU1*V-dV1*U)*2*V*dV2)/(" << VarType << ")(V*V*V*V))";
  F_Lolimot << Lolimot->getCustomReturn_After() << ";" << endl;

  if (Lolimot->getUseTransformLog()) Export_Inverse_Log_Transform_Deriv_2(F_Lolimot, "Result",
									  Lolimot->getTransformLogAlpha(),
									  Lolimot->getTransformLogMin(),
									  Lolimot->getTransformLogMax(),
									  Lolimot->getTransformLogEps());

  F_Lolimot << " return Result;" << endl;

  F_Lolimot << "}" << endl;

  F_Lolimot.close();

  return true;
}
