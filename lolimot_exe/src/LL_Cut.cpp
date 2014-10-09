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
#include <algorithm>

#include "LL_Cut.h"
#include "LL_Mesure.h"
#include "LL_Partition.h"

using namespace std;

LL_Cut::LL_Cut()
{
  _separationSet.resize(0);
  _mesureSet.resize(0);
}

void LL_Cut::Uniform_Cut(LL_Partition * curPartition, unsigned int noDimension, unsigned int NbCut)
{
  float dimensionMinPartition = curPartition->getDimensionMin(noDimension);
  float dimensionMaxPartition = curPartition->getDimensionMax(noDimension);
  float stepSeparation        = (dimensionMaxPartition - dimensionMinPartition) / (float)(NbCut+1);

  _separationSet.clear();
  _separationSet.resize(0);

  for(unsigned int i=1; i<NbCut; i++)
    {
      _separationSet.push_back(dimensionMinPartition + i * stepSeparation);
    } /* End For */
}

void LL_Cut::Distributed_Cut(LL_Partition * curPartition, unsigned int noDimension, unsigned int NbCut)
{
  vector<float> rankedSet;
  float         rsStep = 0.0;
  unsigned int  Index, Step, i;

  rankedSet.clear();
  rankedSet.resize(0);

  // selection of the measures which belong to the partition
  for(i=0; i<_mesureSet.size(); i++)
    {
      if (curPartition->isMesureInPartition(_mesureSet[i]))
	{
	  rankedSet.push_back(_mesureSet[i]->getDimensionValue(noDimension));
	} /* End If */
    } /* End For */

  // We sort the measure which belong to the partition
  sort(rankedSet.begin(), rankedSet.end());

  _separationSet.clear();
  _separationSet.resize(0, 0.0);

  // We can have a problem when a lot of points have the same value :
  // the first bound can have the same value as the min bound. In that case, we are sure to put nothing in the partition
  // which induce a computation problem

  // We compute the step index of rankedSet
  rsStep = rankedSet.size() / (float)NbCut;

  float LastStep = rankedSet[0];

  Index = 1;
  Step  = 0;
  while(Index<NbCut)
    {
      Step += (unsigned int)(rsStep);

      if (rankedSet[Step]!=LastStep)
	{
	  _separationSet.push_back(rankedSet[Step]);
	  LastStep = rankedSet[Step];
	} /* End If */
      else
	{
	  while((rankedSet[Step]==LastStep)&&(Step<rankedSet.size())) Step ++;

	  _separationSet.push_back(rankedSet[Step]);
	  LastStep = rankedSet[Step];

	  // We recompute the step index of rankedSet
	  rsStep = (rankedSet.size() - Step) / (float)NbCut;
	} /* End Else */
      Index++;
    } /* End While */
}

float LL_Cut::getStep(unsigned int noStep)
{
  return _separationSet[noStep];
}

unsigned int LL_Cut::getNbSeparations()
{
  return _separationSet.size();
}

void LL_Cut::addMesures(vector<LL_Mesure *> & Mesures)
{
  _mesureSet = Mesures;
}

vector<LL_Mesure *> & LL_Cut::getMesure()
{
  return _mesureSet;
}
