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
#include <cmath>

using namespace std;

double Inverse_Log_Transform(double Var, double Alpha, double Min, double Max, double eps)
{
  double Result = 0.0;

  if (Var<=Alpha) Result = Alpha*exp(Var/Alpha - 1.0);
  else            Result = Var;

  Result = (Result - eps)*(Max - Min) + Min;

  return Result;
}

double Direct_Log_Transform(double Var, double Alpha, double Min, double Max, double eps)
{
  double Result = 0.0;

  Result = (Var - Min)/(Max - Min) + eps;

  if (Result<=Alpha) Result = Alpha*(1.0 + log(Result/Alpha));

  return Result;
}

void Export_Inverse_Log_Transform(ofstream & OutFile, const string & VarName, 
				  double Alpha, double Min, double Max, double eps)
{
  OutFile << "if (" << VarName << "<=" << Alpha << ")";
  OutFile << "  " << VarName << " = " << Alpha << " * exp(" << VarName << "/" << Alpha << " - 1.0);" << endl;
  OutFile << "  " << VarName << " = (" << VarName << " - " << eps << ") * ((";
  OutFile << Max << ") - (" << Min << ")) + (" << Min << ");" << endl;
}

void Export_Direct_Log_Transform(ofstream & OutFile, const string & VarName,
				 double Alpha, double Min, double Max, double eps)
{
  OutFile << "  " << VarName << " = (" << VarName << " - (" << Min << ")) / ((";
  OutFile << Max << ") - (" << Min << ")) + (" << eps << ");" << endl;

  OutFile << "if (" << VarName << "<=" << Alpha << ")";
  OutFile << "  " << VarName << " = " << Alpha << "*(1.0 + log(" << VarName << "/" << Alpha << "));" << endl;
}

void Export_Inverse_Log_Transform_Deriv(ofstream & OutFile, const string & VarName, 
					double Alpha, double Min, double Max, double eps)
{
  OutFile << "if (" << VarName << "<=" << Alpha << ")" << endl;
  OutFile << "  {" << endl;
  OutFile << VarName << " = " << (Max - Min) << " * exp(" << VarName << "/" << Alpha << " - 1.0);" << endl;
  OutFile << "  }" << endl;
  OutFile << "else" << endl;
  OutFile << "  {" << endl;
  OutFile << "  " << VarName << " = (" << VarName << " * ("<< Max - Min << "));" << endl;
  OutFile << "  }" << endl;
}

void Export_Direct_Log_Transform_Deriv(ofstream & OutFile, const string & VarName,
				       double Alpha, double Min, double Max, double eps)
{
  OutFile << "if (" << VarName << "<=" << Alpha << ")" << endl;
  OutFile << "  {" << endl;
  OutFile << VarName << " = " << Alpha*Alpha*(Max-Min) << " / (";
  OutFile << VarName << " - (" << Min - eps*(Max-Min) << "));" << endl;
  OutFile << "  }" << endl;
  OutFile << "else" << endl;
  OutFile << "  {" << endl;
  OutFile << "  " << VarName << " = " << VarName << " / (" << Max << " - " << Min<< ");" << endl;
  OutFile << "  }" << endl;
}

void Export_Inverse_Log_Transform_Deriv_2(ofstream & OutFile, const string & VarName, 
					  double Alpha, double Min, double Max, double eps)
{
  OutFile << "if (" << VarName << "<=" << Alpha << ")" << endl;
  OutFile << "  {" << endl;
  OutFile << VarName << " = " << (Max - Min)/Alpha << " * exp(" << VarName << "/" << Alpha << " - 1.0);" << endl;
  OutFile << "  }" << endl;
  OutFile << "else" << endl;
  OutFile << "  {" << endl;
  OutFile << "  " << VarName << " = (" << VarName << " * (" << Max - Min << ");" << endl;
  OutFile << "  }" << endl;
}

void Export_Direct_Log_Transform_Deriv_2(ofstream & OutFile, const string & VarName,
					 double Alpha, double Min, double Max, double eps)
{
  OutFile << "if (" << VarName << "<=" << Alpha << ")" << endl;
  OutFile << "  {" << endl;
  OutFile << VarName << " = - (" << Alpha*Alpha*(Max-Min) << ") / pow(";
  OutFile << VarName << " - (" << Min - eps*(Max-Min) << "), 2.0);" << endl;
  OutFile << "  }" << endl;
  OutFile << "else" << endl;
  OutFile << "  {" << endl;
  OutFile << "  " << VarName << " = " << VarName << " / (" << Max << " - " << Min<< ");" << endl;
  OutFile << "  }" << endl;
}
