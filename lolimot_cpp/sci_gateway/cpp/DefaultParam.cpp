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
#include <limits>
#include <math.h>

using namespace std;

float DefFctResidu(vector<float> & Measure,
		   vector<float> & Estim,
		   vector<float> & ParamList)
{
  float Res = 0;
  
  for(unsigned int i=0; i<Measure.size(); i++)
    {
      Res += fabs(Measure[i] - Estim[i]);
    } /* End For */

  return Res/(float)Measure.size();
}

float DefFctResidu_L1(vector<float> & Measure,
		      vector<float> & Estim,
		      vector<float> & ParamList)
{
  float       Res = 0;
  unsigned int i;

  for(i=0; i<Measure.size(); i++)
    {
      Res += pow(Measure[i] - Estim[i], 2.0);
    } /* End For */

  return Res/(float)Measure.size();
}

float DefFctResidu_L1_Rel(vector<float> & Measure,
			  vector<float> & Estim,
			  vector<float> & ParamList)
{
  float       Res  = 0.0;
  unsigned int i;

  for(i=0; i<Measure.size(); i++)
    {
      Res += fabs((Measure[i] - Estim[i])/Measure[i]);
    } /* End For */


  return Res/(float)Measure.size();
}

float DefFctResidu_L2(vector<float> & Measure,
		      vector<float> & Estim,
		      vector<float> & ParamList)
{
  float       Res = 0;
  unsigned int i;

  for(i=0; i<Measure.size(); i++)
    {
      Res += pow(Measure[i] - Estim[i], 2.0);
    } /* End For */

  return sqrt(Res/(float)Measure.size());
}

float DefFctResidu_L2_Rel(vector<float> & Measure,
			  vector<float> & Estim,
			  vector<float> & ParamList)
{
  float       Res = numeric_limits<float>::min();
  float       Tmp = 0.0;
  unsigned int i;

  for(i=0; i<Measure.size(); i++)
    {
      Tmp = fabs(Measure[i] - Estim[i]);

      if (Res<Tmp) Res = Tmp;
    } /* End For */

  return Res;
}

float DefFctResidu_NoNeg(vector<float> & Measure,
			 vector<float> & Estim,
			 vector<float> & ParamList)
{
  float       Res = 0, Aux;
  float       Min = numeric_limits<float>::max();
  unsigned int i;

  for(i=0; i<Measure.size(); i++)
    {
      Aux = pow(Measure[i] - Estim[i], 2.0);
      if (Aux>Min) Aux = Min;
      Res += Aux;
    } /* End For */

  if (Min>0) Min = 0;

  Res += fabs(Min*Measure.size());

  return Res/(float)Measure.size();
}
