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

#ifndef POST_H
#define POST_H

#include <fstream>
#include <string>

using namespace std;

double Direct_Log_Transform(double Var, double Alpha, double Min, double Max, double eps);
double Inverse_Log_Transform(double Var, double Alpha, double Min, double Max, double eps);
void Export_Direct_Log_Transform(ofstream & OutFile, const string & VarName,
				 double Alpha, double Min, double Max, double eps);
void Export_Inverse_Log_Transform(ofstream & OutFile, const string & VarName,
				  double Alpha, double Min, double Max, double eps);
void Export_Direct_Log_Transform_Deriv(ofstream & OutFile, const string & VarName, 
				       double Alpha, double Min, double Max, double eps);
void Export_Inverse_Log_Transform_Deriv(ofstream & OutFile, const string & VarName, 
					double Alpha, double Min, double Max, double eps);
void Export_Direct_Log_Transform_Deriv_2(ofstream & OutFile, const string & VarName, 
					 double Alpha, double Min, double Max, double eps);
void Export_Inverse_Log_Transform_Deriv_2(ofstream & OutFile, const string & VarName, 
					  double Alpha, double Min, double Max, double eps);
#endif
