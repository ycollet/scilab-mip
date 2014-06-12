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

#include "LL_Partition.h"
#include "LL_Mesure.h"
#include "LL_Lolimot.h"
#include "LL_Dimension.h"

#include <stdlib.h>
#include <math.h>

#include <sstream>

using namespace std;

LL_Partition::LL_Partition()
{
  _coeff0 = 0.0;
  _freezeCoeff0 = false;
  _residu = -1;
  _coeffSet.resize(2, 0.0);
  _freezeCoeffSet.resize(2, false);
  _dimensionMinSet.resize(2, 0.0);
  _dimensionMaxSet.resize(2, 0.0);
  _DoNotCut = false;
  _ecartTypeSet.resize(2, 0.0);
  _Name = "none";
}

LL_Partition::LL_Partition(int nbDimensions)
{
  _coeff0 = 0.0;
  _freezeCoeff0 = false;
  _residu = -1;
  _coeffSet.resize(nbDimensions, 0.0);
  _freezeCoeffSet.resize(nbDimensions, false);
  _dimensionMinSet.resize(nbDimensions, 0.0);
  _dimensionMaxSet.resize(nbDimensions, 0.0);
  _DoNotCut = false;
  _inhibeSet.resize(nbDimensions, false);
  _ecartTypeSet.resize(nbDimensions, 0.0);
  _Name = "none";
}

LL_Partition::LL_Partition(const LL_Partition & Var)
{
  _dimensionMinSet = Var._dimensionMinSet;
  _dimensionMaxSet = Var._dimensionMaxSet;
  _coeffSet        = Var._coeffSet;
  _coeff0          = Var._coeff0; 
  _freezeCoeffSet  = Var._freezeCoeffSet;
  _freezeCoeff0    = Var._freezeCoeff0;
  _residu          = Var._residu;
  _DoNotCut        = Var._DoNotCut;
  _ecartTypeSet    = Var._ecartTypeSet;
  strtmp           = Var.strtmp;
  _inhibeSet       = Var._inhibeSet;
  _Name            = Var._Name;
}

LL_Partition::~LL_Partition()
{
}


void LL_Partition::setCoeff(int i, float value)
{
  _coeffSet[i] = value;
}

void LL_Partition::setNbDimension(int nbDimensions)
{
  _coeff0 = 0.0;
  _freezeCoeff0 = false;
  _residu = -1;
  _coeffSet.resize(nbDimensions, 0.0);
  _freezeCoeffSet.resize(nbDimensions, false);
  _dimensionMinSet.resize(nbDimensions, 0.0);
  _dimensionMaxSet.resize(nbDimensions, 0.0);
  _DoNotCut = false;
  _inhibeSet.resize(nbDimensions, false);
  _ecartTypeSet.resize(nbDimensions, 0.0);
}


void LL_Partition::setDimensionMin(int i, float dimensionMin)
{
  _dimensionMinSet[i] = dimensionMin;
}


void LL_Partition::setDimensionMax(int i, float dimensionMax)
{
  _dimensionMaxSet[i] = dimensionMax;
}

string LL_Partition::getRandName(string Prefix, unsigned int NbRandChar)
{
  unsigned int i;
  stringstream Result;

  Result.str("");
  Result << Prefix;

  for(i=0; i<NbRandChar; i++)
    {
      Result << (char)(rand()/(float)RAND_MAX * ('a'-'z') + 'z');
    } /* End For */

  return Result.str();
}

// this method copies the dimensions of curPartition in this partition
bool LL_Partition::setDimensions(const LL_Partition * curPartition, LL_Lolimot * Lolimot)
{
  int    noDimension, nbDimension = curPartition->getDimensionMinSet().size();
  bool   Same = false;
  string PartName;
  unsigned int   i;

  // Computation of a random name

  PartName = getRandName("part_", 5);

  for(i=0; i<Lolimot->getPartitionSet().size(); i++)
    {
      Same = Same || (PartName==Lolimot->getPartitionSet()[i]->getName());
    } /* End For */

  while(Same)
    {
      PartName = getRandName("part_", 5);
      
      Same = false;
      
      for(i=0; i<Lolimot->getPartitionSet().size(); i++)
	{
	  Same = Same || (PartName==Lolimot->getPartitionSet()[i]->getName());
	} /* End For */
    } /* End While */

  _Name = PartName;

  _coeff0 = curPartition->getCoeff0();
  //_freezeCoeff0 = false;
  _freezeCoeff0 = curPartition->getFreezeCoeff0();
  _coeffSet.resize(nbDimension, 0.0);
  _freezeCoeffSet.resize(nbDimension, false);
  _dimensionMinSet.resize(nbDimension, 0.0);
  _dimensionMaxSet.resize(nbDimension, 0.0);
  _ecartTypeSet.resize(nbDimension, 0.0);
  _DoNotCut = curPartition->DoNotCut();
  _inhibeSet.resize(nbDimension, false);

  for(noDimension = 0; noDimension < nbDimension; noDimension ++)
    {
      setDimensionMin(noDimension, curPartition->getDimensionMin(noDimension));
      setDimensionMax(noDimension, curPartition->getDimensionMax(noDimension));
      setCoeff(noDimension, curPartition->getCoeff(noDimension));
      setFreezeCoeff(noDimension, curPartition->getFreezeCoeff(noDimension));
      setInhibe(noDimension, curPartition->getInhibe(noDimension));
    } /* End For */

  return true;
}

void LL_Partition::setInhibe(int i, bool Value)
{
  _inhibeSet[i] = Value;
}

string LL_Partition::printDimensions()
{
  stringstream streamtmp;

  int noDimension, nbDimension = _dimensionMinSet.size();

  // the name of the partition
  streamtmp << "partition name : " << getName() << " - ";

  for(noDimension = 0; noDimension < nbDimension; noDimension ++)
    {
      streamtmp << "[" << getDimensionMin(noDimension) << ", " << getDimensionMax(noDimension) << "]";
    }

  // coefficients
  streamtmp << "(" << getCoeff0();

  for(noDimension = 0; noDimension < nbDimension; noDimension ++)
    {
      streamtmp << " , " << getCoeff(noDimension);
    }

  streamtmp << ")";

  strtmp = streamtmp.str();

  return strtmp;
}


bool LL_Partition::isMesureInPartition(const LL_Mesure * curMesure) const
{
  unsigned int noDimension, nbDimension = _dimensionMinSet.size();

  for(noDimension = 0; noDimension < nbDimension; noDimension ++)
    {
      if ((curMesure->getDimensionValue(noDimension) < this->getDimensionMin(noDimension)) ||
	  (curMesure->getDimensionValue(noDimension) > this->getDimensionMax(noDimension)))
	return false;
    } /* End For */

  // if the mesure repects all the dimensions, then the measure is in the partition
  return true;
}

void LL_Partition::updateEcartType(const LL_Lolimot * curLolimot)
{
  int d;
  int nbDimensions = curLolimot->getDimensionSet().size();
  const LL_Dimension * curDimension;

  for(d=0; d < nbDimensions; d++)
    {
      curDimension     = curLolimot->getDimension(d);
      _ecartTypeSet[d] = curDimension->getEcartType()
	                 * (getDimensionMax(d) - getDimensionMin(d));
    } /* End For */
}


// compute z for given values for each dimension
// z = e(-1/2 (sum[d]   (x1 - centre(d))^2 / std(d)^2 ))
// d belongs to the set of dimensions
float LL_Partition::calculeZ(const vector<float> & x) const
{
  int d;
  int nbDimensions = _ecartTypeSet.size();
  float res = 0.0;

  // we perform the sum
  for(d=0; d < nbDimensions; d++)
    {
      if (!_inhibeSet[d])
	res += pow(x[d] - getCentre(d), (float)2.0) / (float)pow(getEcartType(d), (float)2.0);
    } /* End For */

  res = exp(-0.5 * res);

  return res;
}


// we compute the resulting value for given values for each dimension
// using the regression coefficients
// val = coeff0 + sum[d] (coeff(d) * x(d))
// d belongs to the set of dimensions

float LL_Partition::calculeValeurSelonRegression(const vector<float> & x) const
{
  int d;
  int nbDimensions = _ecartTypeSet.size();
  float res = getCoeff0();

  // we perform the sum
  for(d=0; d < nbDimensions; d++)
    {
      if (!_inhibeSet[d])
	res += getCoeff(d) * x[d];
    } /* End For */

  return res;
}

LL_Partition & LL_Partition::operator=(LL_Partition & right)
{
  _dimensionMinSet = right._dimensionMinSet;
  _dimensionMaxSet = right._dimensionMaxSet;
  _coeffSet        = right._coeffSet;
  _coeff0          = right._coeff0; 
  _freezeCoeffSet  = right._freezeCoeffSet;
  _freezeCoeff0    = right._freezeCoeff0;
  _residu          = right._residu;
  _DoNotCut        = right._DoNotCut;
  _ecartTypeSet    = right._ecartTypeSet;
  strtmp           = right.strtmp;
  _inhibeSet       = right._inhibeSet;

  return *this;
}
