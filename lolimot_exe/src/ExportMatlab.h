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

#ifndef EXPORTMATLAB_H
#define EXPORTMATLAB_H

#include <fstream>
#include <string>
#include <vector>

#include "LL_Lolimot.h"
#include "TrainLolimotStruct.h"

void export_Destroy_Matlab(LL_Lolimot * Lolimot, std::ofstream & F_Lolimot);
void export_Init_Data_Matlab(LL_Lolimot * Lolimot, std::ofstream & F_Lolimot);
void export_Destroy_Data_Matlab(LL_Lolimot * Lolimot, std::ofstream & F_Lolimot);
void export_Declare_Variable_Matlab(LL_Lolimot * Lolimot, std::ofstream & F_Lolimot, const std::string & VarName);
void export_Data_Matlab_Text(LL_Lolimot * Lolimot, std::ofstream & F_Lolimot);
void export_dV_x_Matlab_Text(LL_Lolimot * Lolimot, std::ofstream & F_Lolimot, const std::string & ParamName, const std::string & VarName);
void export_dU_x_Matlab_Text(LL_Lolimot * Lolimot, std::ofstream & F_Lolimot, const std::string & ParamName, const std::string & VarName);
void export_Init_Matlab(LL_Lolimot * Lolimot, std::ofstream & F_Lolimot);
void export_Ei_x_Matlab(LL_Lolimot * Lolimot, std::ofstream & F_Lolimot);
void export_Li_x_Matlab(LL_Lolimot * Lolimot, std::ofstream & F_Lolimot);
void export_U_x_Matlab(LL_Lolimot * Lolimot, std::ofstream & F_Lolimot, const std::string & VarName);
void export_V_x_Matlab(LL_Lolimot * Lolimot, std::ofstream & F_Lolimot, const std::string & VarName);
void export_dU_x_Matlab(LL_Lolimot * Lolimot, std::ofstream & F_Lolimot, unsigned int VarNb1, const std::string & VarName);
void export_dV_x_Matlab(LL_Lolimot * Lolimot, std::ofstream & F_Lolimot, unsigned int VarNb1, const std::string & VarName);
void export_List_Partitions_Matlab(LL_Lolimot * Lolimot, std::ofstream & F_Lolimot);
void export_List_Variables_Matlab(LL_Lolimot * Lolimot, std::ofstream & F_Lolimot);
void export_ddU_x_Matlab(LL_Lolimot * Lolimot, std::ofstream & F_Lolimot,
			 unsigned int VarNb1, unsigned int VarNb2, const std::string & VarName);
void export_ddV_x_Matlab(LL_Lolimot * Lolimot, std::ofstream & F_Lolimot,
			 unsigned int VarNb1, unsigned int VarNb2, const std::string & VarName);
void export_ddV_x_Matlab_Text(LL_Lolimot * Lolimot, std::ofstream & F_Lolimot, const std::string & ParamName1, 
			      const std::string & ParamName2, const std::string & VarType, const std::string & VarName);
void export_ddU_x_Matlab_Text(LL_Lolimot * Lolimot, std::ofstream & F_Lolimot, const std::string & ParamName1, 
			      const std::string & ParamName2, const std::string & VarType, const std::string & VarName);
bool exportFunctionInMatlab(LL_Lolimot * Lolimot, const std::string & path, const std::string & nameFunction, 
			    const std::vector<Model_Input> & ListInput);
bool exportDerivativeFunctionInMatlab(LL_Lolimot * Lolimot, const std::string & path, const std::string & nameFunction,
				      unsigned int VarNb, const std::vector<Model_Input> & ListInput);
bool exportDerivativeFunctionInMatlab2(LL_Lolimot * Lolimot, const std::string & path, const std::string & nameFunction,
				       const std::vector<Model_Input> & ListInput);
bool exportSecondDerivativeFunctionInMatlab(LL_Lolimot * Lolimot, const std::string & path, const std::string & nameFunction,
					    unsigned int VarNb1, unsigned int VarNb2,
					    const std::vector<Model_Input> & ListInput);
bool exportSecondDerivativeFunctionInMatlab2(LL_Lolimot * Lolimot, const std::string & path, const std::string & nameFunction, 
					     const std::vector<Model_Input> & ListInput);
#endif
