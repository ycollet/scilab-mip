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

#ifndef EXPORT_CPP_H
#define EXPORT_CPP_H

#include <LL_Lolimot.h>
#include <TrainLolimotStruct.h>

#include <fstream>
#include <string>
#include <vector>

using namespace std;

bool exportHeaderInCpp2(LL_Lolimot * Lolimot,
			const string & nameFunction, 
			const vector<Model_Input> & ListInput);
bool exportHeaderInCpp(LL_Lolimot * Lolimot,
		       const string & nameFunction, 
		       const vector<Model_Input> & ListInput);
void export_Init_Cpp(LL_Lolimot * Lolimot, 
		     ofstream & F_Lolimot, 
		     const string & VarType);
void export_Destroy_Cpp(LL_Lolimot * Lolimot, 
			ofstream & F_Lolimot, 
			const string & VarType);
void export_Init_Data_Cpp(LL_Lolimot * Lolimot, 
			  ofstream & F_Lolimot, 
			  const string & VarType);
void export_Destroy_Data_Cpp(LL_Lolimot * Lolimot, 
			     ofstream & F_Lolimot, 
			     const string & VarType);
bool exportFunctionInCpp(LL_Lolimot * Lolimot, 
			 const string & nameFunction, 
			 const vector<Model_Input> & ListInput);
bool exportDerivativeFunctionInCpp(LL_Lolimot * Lolimot, 
				   const string & nameFunction, 
				   unsigned int VarNb, 
				   const vector<Model_Input> & ListInput);
bool exportDerivativeFunctionInCpp2(LL_Lolimot * Lolimot, 
				    const string & nameFunction,
				    const vector<Model_Input> & ListInput);
bool exportSecondDerivativeFunctionInCpp(LL_Lolimot * Lolimot, 
					 const string & nameFunction,
					 unsigned int VarNb1, 
					 unsigned int VarNb2, 
					 const vector<Model_Input> & ListInput);
bool exportSecondDerivativeFunctionInCpp2(LL_Lolimot * Lolimot, 
					  const string & nameFunction,
					  const vector<Model_Input> & ListInput);
#endif
