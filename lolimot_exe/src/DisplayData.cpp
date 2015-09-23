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

#include <vector>
#include <iostream>

using namespace std;

void Display_DisplayData(vector<vector<float> > & InputData,
			 vector<vector<float> > & MeasureData)
{
  unsigned int i, j;

  for(i = 0; i<InputData.size(); i++)
    {
      cout << "Input " << i << " = ";
      for(j = 0; j<InputData[i].size(); j++)
	{
	  cout << InputData[i][j] << " ";
	} /* End For */
      cout << endl;
    } /* End For */
  
  for(i = 0; i<MeasureData.size(); i++)
    {
      cout << "Measure " << i << " = ";
      for(j = 0; j<MeasureData[i].size(); j++)
	{
	  cout << MeasureData[i][j] << " ";
	} /* End For */
      cout << endl;
    } /* End For */
}

void Display_DisplayData(vector<vector<float> > & InputData)
{
  unsigned int i, j;

  for(i = 0; i<InputData.size(); i++)
    {
      cout << "Input " << i << " = ";
      for(j = 0; j<InputData[i].size(); j++)
	{
	  cout << InputData[i][j] << " ";
	} /* End For */
      cout << endl;
    } /* End For */
}