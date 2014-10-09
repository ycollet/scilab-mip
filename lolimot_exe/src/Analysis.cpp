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
#include <algorithm>

using namespace std;

#include "LL_Lolimot.h"
#include "LL_Partition.h"
#include "LL_Dimension.h"

// C1

void Analysis_Response_NL(LL_Lolimot * Lolimot, vector<float> & Result)
{
  float       Mean  = 0.0;
  float       Total = 0.0;
  unsigned int i, j;

  Result.resize(0);

  // Computation of the mean distances

  for(i=0; i<Lolimot->getDimensionSet().size(); i++)
    {
      Mean = 0.0;
      for(j=0; j<Lolimot->getPartitionSet().size(); j++)
	{
	  Mean += Lolimot->getPartitionSet()[j]->getEcartType(i);
	} /* End For */

      Mean = Mean / (float)Lolimot->getPartitionSet().size();

      Result.push_back(Mean);
    } /* End For */

  // Computation of C1(xi)
  for(i=0; i<Result.size(); i++)
    {
      Total += (1/(float)Result[i] - 1.0);
    } /* End For */

  for(i=0; i<Result.size(); i++)
    {
      Result[i] = 100*(1/(float)Result[i] - 1.0)/(float)Total;
    } /* End For */
}

// C2

void Analysis_Response(LL_Lolimot * Lolimot, vector<float> & Result,
		       vector<float> & InputMin, vector<float> & InputMax)
{
  float       Mean  = 0.0;
  float       Total = 0.0;
  unsigned int i, j;

  Result.resize(0);

  // Computation of the mean coefficients (the +1 is here to take into account the constant)

  for(i=0; i<Lolimot->getDimensionSet().size(); i++)
    {
      Mean = 0.0;
      for(j=0; j<Lolimot->getPartitionSet().size(); j++)
	{
	  if (i==0)
	    {
	      Mean += fabs(Lolimot->getPartitionSet()[j]->getCoeff0()/(float)(InputMax[0]-InputMin[0]));
	    } /* End If */
	  else
	    {
	      Mean += fabs(Lolimot->getPartitionSet()[j]->getCoeff(i)/(float)(InputMax[j]-InputMin[j]));
	    } /* End Else */
	} /* End For */

      Mean = Mean / (float)Lolimot->getPartitionSet().size();

      Result.push_back(Mean);
    } /* End For */

  // Computation of C2(xi)
  for(i=0; i<Result.size(); i++)
    {
      Total += Result[i];
    } /* End For */

  for(i=0; i<Result.size(); i++)
    {
      Result[i] = 100.0*Result[i]/(float)Total;
    } /* End For */
}

void Analysis_List_DistribOfCut(LL_Lolimot * Lolimot, string Filename)
{
  ofstream       Outfile;
  vector<double> ListOfCut, ListOfCutDest;
  vector<double>::iterator itListOfCut;
  unsigned int   i, j;

  Outfile.open(Filename.c_str());
  if (!Outfile.is_open())
    {
      cerr << "Analysis_List_DistribOfCut: Can't open file " << Filename << endl;
      exit(1);
    } /* End If */

  for(i=0; i<Lolimot->getDimensionSet().size(); i++)
    {
      ListOfCut.resize(0);
      ListOfCutDest.resize(0);
      // We list all the value of the cut of a dimension
      for(j=0; j<Lolimot->getPartitionSet().size(); j++)
	{
	  ListOfCut.push_back(Lolimot->getPartitionSet()[j]->getDimensionMin(i));
	  ListOfCut.push_back(Lolimot->getPartitionSet()[j]->getDimensionMax(i));
	} /* End For */
      // We remove all duplicated values
      // First we sort the list
      sort(ListOfCut.begin(), ListOfCut.end());
      // Next we remove equal adjacent points
      itListOfCut = unique(ListOfCut.begin(), ListOfCut.end());
      // And we copy this ordered list into a vector
      copy(ListOfCut.begin(), itListOfCut, back_inserter(ListOfCutDest));
      ListOfCut = ListOfCutDest;

      Outfile << Lolimot->getDimensionSet()[i]->getName() << " : ";
      for(j=0; j<ListOfCut.size(); j++) Outfile << ListOfCut[j] << " ";
      Outfile << endl;
    } /* End For */

  Outfile.close();
}
