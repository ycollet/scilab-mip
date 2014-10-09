// Authors: L. Pajou, Y. Collette
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

/////////////////////////////////////////////////////////////////////
// Description        : Modelisation using lolimot functions       //
//                                                                 //
// Date of creation   :      11/12/2002                            //
//                                                                 //
// Nom du developpeur : Laurent PAJOU                              //
// Comment            :                           Ste EURODECISION //
// Modified by        : Yann COLLETTE (Renault)                    //
/////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "LL_Lolimot.h"
#include "LL_Dimension.h"
#include "LL_Mesure.h"
#include "LL_Partition.h"
#include "LL_Cut.h"
#include "Post.h"

extern "C" {
#include <lasso.h>
}

#ifdef WIN32
#define UP "\x1b[A"
#else
#define UP "\033[A"
#endif
#define NANVALUE 0.0

#ifdef WIN32
#define isnan(A) (false)
#endif

#if WIN32
extern "C" 
{
#include <f2c.h>
#include <blaswrap.h>
} 
#endif

// we remove 2 #define from f2c which interfere with the method numeric_limits<float>::max() and
// numeric_limits<float>::min()

#undef max
#undef min

//#define DEBUG 1

using namespace std;

#if !defined(WIN32)
#define USEDOUBLE 1
typedef long int integer;
#ifdef USEDOUBLE
typedef double doublereal;
#else
typedef float doublereal; 
#endif
#undef USEDOUBLE
#endif

// declaration of the least square function from Lapack.
extern "C" int dgelss_(integer *, integer *, integer *,
		       doublereal *, integer *, doublereal *, integer *, 
		       doublereal *,  doublereal *, integer *,
		       doublereal *, integer *, integer *);

void funcs(float * x, float * afunc, int ma)
{
  int i;

  for(i=2; i<ma+1; i++)
    {
      afunc[i] = x[i];
    } /* End For */

  afunc[1]=1;
}

int compareDouble( const void *arg1, const void *arg2 )
{
  float d1 = *((float*)arg1);
  float d2 = *((float*)arg2);
  if (d1 < d2)
    return -1;
  else
    return 1;
}

int LL_Lolimot::getNbMaxPartitions() const
{
  return _nbMaxPartitions;
}

float LL_Lolimot::getResiduGapPercentage() const
{
  return _residuGapPercentage;
}

void LL_Lolimot::setNbMaxPartitions(int nb)
{
  _nbMaxPartitions = nb;
}

void LL_Lolimot::setResiduGapPercentage(float res)
{
  _residuGapPercentage = res;
}

float LL_Lolimot::getResidu() const
{
  return _residu;
}

void LL_Lolimot::setCustomArguments(string Arguments)
{
  _customArguments = Arguments;
}

string LL_Lolimot::getCustomArguments() const
{
  return _customArguments;
}

void LL_Lolimot::setCustomReturn_Before(string Arguments)
{
  _customReturn_Before = Arguments;
}

string LL_Lolimot::getCustomReturn_Before() const
{
  return _customReturn_Before;
}

void LL_Lolimot::setCustomReturn_After(string Arguments)
{
  _customReturn_After = Arguments;
}

string LL_Lolimot::getCustomReturn_After() const
{
  return _customReturn_After;
}

float LL_Lolimot::getInitialResidu() const
{
  return _initialResidu;
}

void LL_Lolimot::setInitialResidu(float InitRes)
{
  _initialResidu = InitRes;
}

vector<LL_Dimension*> LL_Lolimot::getDimensionSet() const
{
  return _dimensionSet;
}

LL_Dimension* LL_Lolimot::getDimension(int i) const
{
  return _dimensionSet[i];
}

vector<LL_Mesure*> LL_Lolimot::getMesureSet() const
{
  return _mesureSet;
}

LL_Mesure* LL_Lolimot::getMesure(int i) const
{
  return _mesureSet[i];
}

vector<LL_Partition*> LL_Lolimot::getPartitionSet() const
{
  return _partitionSet;
}

LL_Partition* LL_Lolimot::getPartition(int i) const
{
  return _partitionSet[i];
}

void LL_Lolimot::setNbParam(unsigned int nbParam)
{
  _paramList.resize(nbParam, 0.0);
}

unsigned int LL_Lolimot::getNbParam() const
{
  return _paramList.size();
}

void LL_Lolimot::setParam(unsigned int nbParam, float ParamValue)
{
  _paramList[nbParam] = ParamValue;
}

float LL_Lolimot::getParam(unsigned int nbParam) const
{
  return _paramList[nbParam];
}

void LL_Lolimot::addFctResidu(PtrToFctResidu FctResidu)
{
  _FctResidu = FctResidu;
}

void LL_Lolimot::setMaxNbPointsPerPartitions(int MaxPtPerPart)
{
  _MaxNbPointsPerPartitions = MaxPtPerPart;
}

int LL_Lolimot::getMaxNbPointsPerPartitions() const
{
  return _MaxNbPointsPerPartitions;
}

bool LL_Lolimot::isThereAnyPartitionAvailable() const
{
  return _anyPartitionAvailable;
}

void LL_Lolimot::clearAnyPartitionAvailable()
{
  _anyPartitionAvailable = true;
}

void LL_Lolimot::setUseUniformCutting(bool Value)
{
  _useuniformcutting     = Value;
  _usedistributedcutting = !Value;
}

bool LL_Lolimot::getUseUniformCutting() const
{
  return _useuniformcutting;
}

void LL_Lolimot::setUseLasso(bool Value)
{
  _useLasso = Value;
}

bool LL_Lolimot::getUseLasso() const
{
  return _useLasso;
}

void LL_Lolimot::setUseDistributedCutting(bool Value)
{
  _usedistributedcutting = Value;
  _useuniformcutting     = !Value;
}

bool LL_Lolimot::getUseDistributedCutting() const
{
  return _usedistributedcutting;
}

void LL_Lolimot::setMembershipThreshold(float Value)
{
  _MembershipThreshold = Value;
}

float LL_Lolimot::getMembershipThreshold() const
{
  return _MembershipThreshold;
}

LL_Lolimot::LL_Lolimot()
{
  _opt_continue              = false;
  _MaxNbPointsPerPartitions  = 3;
  _nbMaxPartitions           = 0;
  _residuGapPercentage       = 0.0;
  _residu                    = 0.0;
  _initialResidu             = 0.0;
  _anyPartitionAvailable     = true;

  _dimensionSet.resize(0);
  _mesureSet.resize(0);
  _partitionSet.resize(0);

  _FctResidu = DefFctResidu;
  _cutFcn = new LL_Cut;

  _useuniformcutting            = true;
  _usedistributedcutting        = false;
  _exportTypeDouble             = true;
  _Weight.resize(0);
  _customArguments     = "";
  _customReturn_Before = "";
  _customReturn_After  = "";
  _useTransformLog     = false;
  _TransformLogAlpha   = 1.0;
  _TransformLogMin     = 0.0;
  _TransformLogMax     = 1.0;
  _TransformLogEps     = 0.0;
  _MembershipThreshold = 0.0;
  _useLasso            = false;
}

LL_Lolimot::LL_Lolimot(string name, int nbMaxPartitions, float residuGapPercentage)
{
  _opt_continue              = false;
  _MaxNbPointsPerPartitions  = 3;
  _nbMaxPartitions           = nbMaxPartitions;
  _residuGapPercentage       = residuGapPercentage;
  _residu                    = 0.0;
  _initialResidu             = 0.0;
  _FctResidu                 = DefFctResidu;
  _anyPartitionAvailable     = true;

  _dimensionSet.resize(0);
  _mesureSet.resize(0);
  _partitionSet.resize(0);

  _FctResidu = DefFctResidu;
  _cutFcn = new LL_Cut;

  _useuniformcutting     = true;
  _usedistributedcutting = false;
  _exportTypeDouble      = true;
  _Weight.resize(0);
  _customArguments     = "";
  _customReturn_Before = "";
  _customReturn_After  = "";
  _useTransformLog     = false;
  _TransformLogAlpha   = 1.0;
  _TransformLogMin     = 0.0;
  _TransformLogMax     = 1.0;
  _TransformLogEps     = 0.0;
  _MembershipThreshold = 0.0;
  _useLasso            = false;
}

LL_Lolimot::~LL_Lolimot()
{
  int            noDimension, nbDimension = getDimensionSet().size();
  int            noMesure,    nbMesure    = _mesureSet.size();
  int            noPartition, nbPartition = _partitionSet.size();
  LL_Dimension * curDimension = NULL;
  LL_Mesure    * curMesure    = NULL;
  LL_Partition * curPartition = NULL;

  for(noDimension = 0; noDimension < nbDimension; noDimension ++)
    {
      curDimension = _dimensionSet[noDimension];
      delete curDimension;
    } /* End For */

  for(noMesure = 0; noMesure < nbMesure; noMesure ++)
    {
      curMesure = _mesureSet[noMesure];
      delete curMesure;
    } /* End For */

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = _partitionSet[noPartition];
      delete curPartition;
    } /* End For */

  delete _cutFcn;
}

LL_Lolimot & LL_Lolimot::operator=(const LL_Lolimot & right)
{
  unsigned int i;

  _MaxNbPointsPerPartitions  = right._MaxNbPointsPerPartitions;
  _nbMaxPartitions           = right.getNbMaxPartitions();
  _residuGapPercentage       = right.getResiduGapPercentage();
  _residu                    = right.getResidu();
  _initialResidu             = right.getInitialResidu();
  _anyPartitionAvailable     = right._anyPartitionAvailable;

  _useuniformcutting            = right._useuniformcutting;
  _usedistributedcutting        = right._usedistributedcutting;
  _exportTypeDouble             = right._exportTypeDouble;

  for(i=0; i<_dimensionSet.size(); i++)
    {
      if (_dimensionSet[i]) delete _dimensionSet[i];
    } /* End For i */
  
  _dimensionSet.resize(right.getDimensionSet().size());

  for(i=0; i<right.getDimensionSet().size(); i++)
    {
      _dimensionSet[i] = new LL_Dimension;
      (*_dimensionSet[i]) = (*right.getDimensionSet()[i]);
    } /* End For i */

  for(i=0; i<_mesureSet.size(); i++)
    {
      if (_mesureSet[i]) delete _mesureSet[i];
    } /* End For i */
  
  _mesureSet.resize(right.getMesureSet().size());

  for(i=0; i<right.getMesureSet().size(); i++)
    {
      _mesureSet[i] = new LL_Mesure;
      (*_mesureSet[i]) = (*right.getMesureSet()[i]);
    } /* End For i */

  for(i=0; i<_partitionSet.size(); i++)
    {
      if (_partitionSet[i]) delete _partitionSet[i];
    } /* End For i */
  
  _partitionSet.resize(right.getPartitionSet().size());

  for(i=0; i<right.getPartitionSet().size(); i++)
    {
      _partitionSet[i] = new LL_Partition();
      (*_partitionSet[i]) = (*right.getPartitionSet()[i]);
    } /* End For i */

  _Weight   = right._Weight;

  _customArguments     = right._customArguments;
  _customReturn_Before = right._customReturn_Before;
  _customReturn_After  = right._customReturn_After;

  _MembershipThreshold = right._MembershipThreshold;
  _useLasso            = right._useLasso;

  return *this;
}

void LL_Lolimot::setSigma(float SigmaVar)
{
  unsigned int noDimension;

  for(noDimension=0; noDimension<_dimensionSet.size(); noDimension++)
    {
      _dimensionSet[noDimension]->setEcartType(SigmaVar);
    } /* End For */
}

float LL_Lolimot::getSigma() const
{
  float Aux;

  if (_dimensionSet.size()==0) Aux = 0;
  else                         Aux = _dimensionSet[0]->getEcartType();

  return Aux;
}

LL_Dimension * LL_Lolimot::addDimension(string name, float min, float max, float ecartType, int nbDiscretisations)
{
  LL_Dimension * curDimension = new LL_Dimension(name, min, max, ecartType, nbDiscretisations);

  _dimensionSet.push_back(curDimension);

  _MaxNbPointsPerPartitions = _dimensionSet.size() + 1;

  ListOfVarToDontCut.resize(_dimensionSet.size(), false);

  return curDimension;
}


LL_Mesure * LL_Lolimot::addMesure(float value, float minValueExtrapolation, float maxValueExtrapolation)
{
  float AuxMesure;

  if (_useTransformLog) AuxMesure = Direct_Log_Transform(value,
							 _TransformLogAlpha,
							 _TransformLogMin,
							 _TransformLogMax,
							 _TransformLogEps);
  else                  AuxMesure = value;

  LL_Mesure * curMesure = new LL_Mesure(AuxMesure, _dimensionSet.size(), minValueExtrapolation, maxValueExtrapolation);

  _mesureSet.push_back(curMesure);

  return curMesure;
}

string LL_Lolimot::getRandName(string Prefix, unsigned int NbRandChar)
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

LL_Partition * LL_Lolimot::addPartition()
{
  LL_Partition * curPartition = new LL_Partition(_dimensionSet.size());
  bool           Same = false;
  string         PartName;
  unsigned int   i;

  _partitionSet.push_back(curPartition);

  PartName = getRandName("part_", 5);

  for(i=0; i<_partitionSet.size(); i++)
    {
      Same = Same || (PartName==_partitionSet[i]->getName());
    } /* End For */

  while(Same)
    {
      PartName = getRandName("part_", 5);
      
      Same = false;
      
      for(i=0; i<_partitionSet.size(); i++)
	{
	  Same = Same || (PartName==_partitionSet[i]->getName());
	} /* End For */
    } /* End While */
  
  curPartition->setName(PartName);
  
  return curPartition;
}

int LL_Lolimot::cleanMesures()
{
  _mesureSet.clear();
  return 0;
}

int LL_Lolimot::cleanPartitions()
{
  _partitionSet.clear();
  return 0;
}

int LL_Lolimot::cleanDimensions()
{
  _dimensionSet.clear();

  _MaxNbPointsPerPartitions = 3;

  return 0;
}

vector<float> LL_Lolimot::getWeightSet() const
{
  return _Weight;
}

float LL_Lolimot::getWeight(int i) const
{
  return _Weight[i];
}

void LL_Lolimot::addWeight(float Weight)
{
  _Weight.push_back(Weight);
}

void LL_Lolimot::cleanWeight()
{
  _Weight.resize(0);
}

void LL_Lolimot::useDefaultWeight()
{
  _Weight.resize(_mesureSet.size(), 1);
}

bool LL_Lolimot::calculePhi(const vector<float> & x, vector<float> & phi) const
{
  int            p;
  int            nbPartitions = _partitionSet.size();
  LL_Partition * curPartition = NULL;
  float         sommeZ = 0.0;

  // on calcule les z
  for(p=0; p < nbPartitions; p++)
    {
      curPartition = getPartition(p);
      phi[p]       = curPartition->calculeZ(x);
      sommeZ      += phi[p];
    } /* End For */

  // on en deduit phi
  for(p=0; p < nbPartitions; p++)
    phi[p] = phi[p] / (float)sommeZ;

  return true;
}

float LL_Lolimot::calculeF(const vector<float> & x, const vector<float> & phi) const
{
  int            p;
  int            nbPartitions = _partitionSet.size();
  LL_Partition * curPartition = NULL;
  float          res = 0.0;

  for(p=0; p < nbPartitions; p++)
    {
      curPartition = getPartition(p);
      res         += curPartition->calculeValeurSelonRegression(x) * phi[p];
    } /* End For */

  return res;
}

float LL_Lolimot::calculeF_Untransformed(const vector<float> & x, const vector<float> & phi) const
{
  float Result;
  if (_useTransformLog) Result = Inverse_Log_Transform(calculeF(x, phi),
						       _TransformLogAlpha,
						       _TransformLogMin,
						       _TransformLogMax,
						       _TransformLogEps);
  else                  Result = calculeF(x, phi);

  return Result;
}

// compute the estimation for values given for each dimension
// f(x) = sum[p] (ValeurSelonRegression(p) * phi(p))
// p belongs to the set of partitions
// phi(p) = z(p) / sum[p2] z(p2)   (phi(p) is a probability distribution)

float LL_Lolimot::calculeF(const vector<float> & x) const
{
  vector<float> phi(_partitionSet.size(), 0.0);
  calculePhi(x, phi);
  return calculeF(x, phi);
}

float LL_Lolimot::calculeF_Untransformed(const vector<float> & x) const
{
  float Result;

  if (_useTransformLog) Result = Inverse_Log_Transform(calculeF(x),
						       _TransformLogAlpha,
						       _TransformLogMin,
						       _TransformLogMax,
						       _TransformLogEps);
  else                  Result = calculeF(x);

  return Result;
}

bool LL_Lolimot::updatePartitions()
{
  int            noPartition, nbPartition = _partitionSet.size();
  LL_Partition * curPartition = NULL;

  // we update the steady deviation of each partition-dimension
  for(noPartition=0; noPartition < nbPartition; noPartition++)
    {
      curPartition = getPartition(noPartition);
      curPartition->updateEcartType(this);
    } /* End For */

  return true;
}

bool LL_Lolimot::updateCoefficients(float & residu, LL_Partition * & partitionWithBiggestResidu)
{
  bool Result = true;

  if (!_useLasso)
    Result = calculeCoefficients(residu, partitionWithBiggestResidu);
  else
    Result = calculeCoefficients_Lasso(residu, partitionWithBiggestResidu);
    

  return Result;
}

// computation of the coefficients (and the corresponding residual) with respect to the partitions
// (we call the lapack dgelss function for each partition)

bool LL_Lolimot::calculeCoefficients(float & residu, LL_Partition * & partitionWithBiggestResidu)
{
  cerr << "LL_Lolimot::calculeCoefficients " << endl;
  float          residuPartition = 0.0;
  int            noPartition, nbPartition = _partitionSet.size();
  LL_Partition * curPartition = NULL;
  LL_Mesure *    curMesure = NULL;
  int            noDimension, nbDimension = _dimensionSet.size();
  int            noMesure, nbMesure       = _mesureSet.size();
  int            i;
  double         biggestResiduPartition   = 0.0;
  // variables for lapack
  long int     M    = nbMesure;             // nb of lines for the matrix A  
  long int     N    = nbDimension + 1;      // nb of columns for the matrix A  
  long int     NRHS = 1;                    // nb of columns of matrix B and X
  doublereal * ptrA = new doublereal[M * N];    // Matrix A
  long int     LDA  = M;                    // Leading Dimension for A : max (1, M) = M
  long int     LDB;                         // Leading Dimension for B = max(N, M)
  if (M > N)
    LDB = M;
  else
    LDB = N;
  doublereal * ptrB = new doublereal[LDB];      // MAtrix B - will hold the solution

  int     sizeS;
  if (M < N) sizeS = M;
  else       sizeS = N;

  doublereal * ptrS  = new doublereal[sizeS];   // Sum-up : Table of singular values of dimension Min(M,N)
  doublereal   RCOND = 1E-10;                   // float precision
  long int     RANK;                            // Sum-up : rank of A (we don't use this)

  long int     LWORK   = 5*LDB;                 // Storage Table (WORK) of dimension 5*Max(M,N)=5*LDB
  doublereal * ptrWORK = new doublereal[LWORK];
  long int     INFOMC;                          // = 0 : successful
                                                // < 0 : illegal argument
                                                // > 0 : no convergence


  partitionWithBiggestResidu = NULL;

  // we update the steady deviate of each partition-dimension
  for(noPartition=0; noPartition < nbPartition; noPartition++)
    {
      curPartition = getPartition(noPartition);
      curPartition->updateEcartType(this);
    } /* End For */

  // computation of phi

  // phi[e][p] : e : experiment of measure - p : partition
  vector<vector<float> > phiEP;
  phiEP.resize(nbMesure, vector<float>(nbPartition, 0.0));

  // for each experiment, we compute phi[e]
  for(noMesure = 0; noMesure < nbMesure; noMesure ++)
    {
      curMesure = _mesureSet[noMesure];
      calculePhi(curMesure->getDimensionValueSet(), phiEP[noMesure]);
    } /* End For */

  // we build the least square problem and solve it for each partition

  cerr << endl;

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      cerr << UP << "LL_Lolimot: computing partition " << noPartition + 1 << " / " << nbPartition << endl;
      curPartition = _partitionSet[noPartition];

      // we fill the matrix B
      for(noMesure = 0; noMesure < nbMesure; noMesure ++)
	{
	  if (phiEP[noMesure][noPartition]>=_MembershipThreshold)
	    {
	      curMesure = _mesureSet[noMesure];
	      ptrB[noMesure] = phiEP[noMesure][noPartition] * _Weight[noMesure] * curMesure->getValue();
	    } /* End If */
	} /* End For */

      // we fill the matrix A
      for(noMesure = 0; noMesure < nbMesure; noMesure ++)
	{
	  if (phiEP[noMesure][noPartition]>=_MembershipThreshold)
	    {
	      curMesure = _mesureSet[noMesure];
	      
	      for(noDimension = 0; noDimension < nbDimension; noDimension ++)
		{
		  // M = nbMesure
		  ptrA[M * (noDimension+1) + noMesure] = phiEP[noMesure][noPartition] * _Weight[noMesure] 
		    * curMesure->getDimensionValue(noDimension);
		} /* End For */
	      
	      // We fill the first column of A
	      ptrA[noMesure] = phiEP[noMesure][noPartition] * _Weight[noMesure];
	    } /* End If */
	} /* End For */

      // Resolution using the least square method
      // We use LAPACK for the least square problem
      // Parameters (see the definition of the function CLAPACK/dgelss on http://www.netlib.org/):
      // Call to the SVD least square computation (Clapack procedure)
      dgelss_(&M, &N, &NRHS, ptrA, &LDA, ptrB, &LDB, ptrS , &RCOND , &RANK , ptrWORK, &LWORK, &INFOMC);

      // test in case of degeneracy
      //  1) M < N (too few experiments): in that case, the problem has an infinite number of solutions with ||Ax-b||=0.
      // Detection of the situation 1) ==> Warning
      if (M < N)
        {
	  cerr << "** ERROR LL_Lolimot::calculeCoefficients() : Too few experiments." << endl;
        } /* End If */

      //  2) M >=N and rank(A)<N: In that case, the problem has also an infinite number of solutions.
      // Detection of the situation 2) ==> Warning
      // In the case where M >=N : Detection of the case where the matrix is of deficient rank (badly chosen experiments)
      if ((M >= N) && (ptrS[N-1] < 1e-10))
        {
	  // display some informations related to the cases of rank deficient problems
	  cerr << "** ERROR LL_Lolimot::calculeCoefficients() : M >= N and rank(A) < N." << endl;

	  cerr << "   Partition involved : n " << noPartition+1 << " : " << curPartition->printDimensions() << endl;

	  for(i = 0; i < N; i++)
	    {
	      cerr << "    Singular value " << i << " = " << ptrS[i] << endl;
	    } /* End For */

	  delete[] ptrA;
	  delete[] ptrB;
	  delete[] ptrS;
	  delete[] ptrWORK;

	  residu = numeric_limits<float>::max();
	  return true;
        } /* End If */

      // we get the coefficients
      if (isnan(ptrB[0])) ptrB[0] = NANVALUE;
      if (!curPartition->getFreezeCoeff0()) curPartition->setCoeff0(ptrB[0]);

      for(noDimension = 0; noDimension < nbDimension; noDimension ++)
	{
	  if (isnan(ptrB[noDimension + 1])) ptrB[noDimension + 1] = NANVALUE;
	  if (!curPartition->getFreezeCoeff(noDimension)) curPartition->setCoeff(noDimension, ptrB[noDimension + 1]);
	} /* End For */

      residuPartition = 0.0;

      if (M > N)
	{
	  int i;
	  for(i = N; i < M; i++)
	    {
	      residuPartition += pow(ptrB[i], 2.0);
	    } /* End For */
	} /* End For */

      // We memorize the performances of the partition
      if ((residuPartition > biggestResiduPartition))
	{
	  biggestResiduPartition     = residuPartition;
	  partitionWithBiggestResidu = curPartition;
	} /* End If */
    } /* End For */

  // computation of the residual
  _residu = computeGlobalResidual();
  residu  = biggestResiduPartition;
  
  // Destruction of the intermediate arrays
  delete [] ptrA;
  delete [] ptrB;
  delete [] ptrS;
  delete [] ptrWORK; 

  return true;
}

// computation of the coefficients (and the associated residuals) with respect to the partitions

bool LL_Lolimot::calculeCoefficients_Lasso(float & residu, LL_Partition * & partitionWithBiggestResidu)
{
  cerr << "LL_Lolimot::calculeCoefficients_Lasso " << endl;
  float          residuPartition = 0.0;
  int            noPartition, nbPartition = _partitionSet.size();
  LL_Partition * curPartition = NULL;
  LL_Mesure *    curMesure = NULL;
  int            noDimension, nbDimension = _dimensionSet.size();
  int            noMesure, nbMesure       = _mesureSet.size();
  int            i;
  double         biggestResiduPartition   = 0.0;
  // variables for the Lasso method
  double * X = NULL;
  double * Y = NULL, * YHat = NULL;
  double * Residual = NULL;
  double * Param = NULL;
  int      Verbose = 0, PasSub = 0, PSuc, nbDimLasso = nbDimension + 1;
  double   Bound, Lagrangian;

  X          = new double[nbDimLasso*nbMesure];
  Y          = new double[nbMesure];
  YHat       = new double[nbMesure];
  Residual   = new double[nbMesure];
  Param      = new double[nbDimLasso];

  // Initialisation of the model
  for(i=0; i<nbDimension+1; i++) Param[i] = 0.0;

  partitionWithBiggestResidu = NULL;

  // we update the steady deviation of each partition-dimension
  for(noPartition=0; noPartition < nbPartition; noPartition++)
    {
      curPartition = getPartition(noPartition);
      curPartition->updateEcartType(this);
    } /* End For */

  // computation of phi

  // phi[e][p] : e : experiment or measure - p : partition
  vector<vector<float> > phiEP;
  phiEP.resize(nbMesure, vector<float>(nbPartition, 0.0));

  // for each experiment, we compute phi[e]
  for(noMesure = 0; noMesure < nbMesure; noMesure ++)
    {
      curMesure = _mesureSet[noMesure];
      calculePhi(curMesure->getDimensionValueSet(), phiEP[noMesure]);
    } /* End For */

  // we build and solve the least square problem for each partition

  cerr << endl;

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      cerr << UP << "LL_Lolimot: computing partition " << noPartition + 1 << " / " << nbPartition << endl;
      curPartition = _partitionSet[noPartition];

      // we fill the Y matrix
      for(noMesure = 0; noMesure < nbMesure; noMesure ++)
	{
	  if (phiEP[noMesure][noPartition]>=_MembershipThreshold)
	    {
	      curMesure = _mesureSet[noMesure];
	      Y[noMesure] = phiEP[noMesure][noPartition] * _Weight[noMesure] * curMesure->getValue();
	    } /* End If */
	} /* End For */

      // we fill the X matrix
      for(noMesure = 0; noMesure < nbMesure; noMesure ++)
	{
	  if (phiEP[noMesure][noPartition]>=_MembershipThreshold)
	    {
	      curMesure = _mesureSet[noMesure];
	      
	      for(noDimension = 1; noDimension < nbDimLasso; noDimension ++)
		{
		  X[noDimension*nbMesure + noMesure] = phiEP[noMesure][noPartition] * _Weight[noMesure] 
		    * curMesure->getDimensionValue(noDimension);
		} /* End For */
	      
	      // we fill the first column of matrix A
	      X[0*nbMesure + noMesure] = phiEP[noMesure][noPartition] * _Weight[noMesure];
	    } /* End If */
	} /* End For */

      Bound      = numeric_limits<float>::max();
      Lagrangian = 1.0;

      // We solve the problem using the Lasso method (lasso module extracted from the Lasso R package)
      lasso(X, &nbMesure, &nbDimLasso, &Bound, Param, Y, YHat, Residual, &Lagrangian, &PSuc,  &Verbose, &PasSub);

      // Test of degeneracy cases
      //  1) nbMesures < nbDimension (too few experiments): in that case, the problem has an infinite number of solutions
      // with ||Ax-b||=0.
      // Detection of the situation 1) ==> Warning
      if (nbMesure < nbDimension)
        {
	  cerr << "** ERROR LL_Lolimot::calculeCoefficients() : Too few experiments." << endl;
        } /* End If */

      // we get the coefficients
      if (!curPartition->getFreezeCoeff0()) curPartition->setCoeff0(Param[0]);

      for(noDimension = 1; noDimension < nbDimLasso; noDimension ++)
	{
	  if (!curPartition->getFreezeCoeff(noDimension)) curPartition->setCoeff(noDimension, Param[noDimension]);
	} /* End For */

      residuPartition = 0.0;
      for(i=0; i<nbMesure; i++) residuPartition += Residual[i]*Residual[i];

      // We memorize the performances of the partition
      if ((residuPartition > biggestResiduPartition))
	{
	  biggestResiduPartition     = residuPartition;
	  partitionWithBiggestResidu = curPartition;
	} /* End If */
    } /* End For */

  // computation of the residual
  _residu = computeGlobalResidual();
  residu  = biggestResiduPartition;
  
  // destruction of the intermediate storage
  delete [] X;
  delete [] Y;
  delete [] YHat;
  delete [] Residual;
  delete [] Param;

  return true;
}

// method which add the best partition (i.e. which reduces the most the global residual)
// if partitionToSeparate is equal to NULL, we try to separate all the partition and then, 
// we separate the most interesting one

int LL_Lolimot::addTheBestPartition(LL_Partition * partitionToSeparate, LL_Partition * & partitionWithBiggestResidu)
{
  unsigned int   noPartition, nbPartition = _partitionSet.size();
  LL_Partition * curPartition = NULL;
  LL_Dimension * curDimension = NULL;
  unsigned int   noDimension, nbDimension = _dimensionSet.size();
  unsigned int   noSeparation; // nb of possible separations for a dimension
  float          valeurSeparation;
  float          dimensionMinPartition, dimensionMaxPartition;

  // variables to store the best found partition
  float          bestResidu           = numeric_limits<float>::max(); // best residual found
  float          bestResiduPourcent   = bestResidu/100.0; // value of the relative residual
  LL_Partition * bestPartition        = NULL; // partition to cut
  int            bestDimension        = -1; // index of the dimension to cut
  float          bestValeurSeparation = -1; // value of the separation
  float          residu = 0, residuPourcent = 0;
  float          precResidu = 0.0, bestPrecResidu = 0.0;
  LL_Partition * locPartitionWithBiggestResidu = NULL;

  cout << "LL_Lolimot::addTheBestPartition" << endl;

  // In precResidu, we store the value of the preceding residual before cutting.
  // This, so as to compare the improvement of the adequacy between the model and the data
  
  if (updateCoefficients(residu, locPartitionWithBiggestResidu)==false)
    {
      cerr << "** ERROR LL_Lolimot::addTheBestPartition() : during the computation of the models." << endl;
      return 1;
    } /* End If */
  
  precResidu     = residu;
  residuPourcent = residu/100.0;

  // we create a new partition
  LL_Partition * newPartition = addPartition();
  partitionWithBiggestResidu  = NULL;

  // for each already existing partition
  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = _partitionSet[noPartition];

      // if we have given a precise partition to separate
      if (partitionToSeparate) curPartition = partitionToSeparate;

      // we initialize the new partition : the new partition is the same hyper-rectangle
      // than curPartition (we will modify in dimension at a time next)
      if (newPartition->setDimensions(curPartition, this) == false)
	{
	  cerr << "** ERROR LL_Lolimot::addTheBestPartition() : in LL_Partition::setDimensions." << endl;
	  return 1;
	} /* End If */

      for(noDimension = 0; noDimension < nbDimension; noDimension ++)
	{
	  cerr << "LL_Lolimot::addTheBestPartition() : " << _dimensionSet[noDimension]->getName();
	  cerr << " - dimension " << noDimension + 1 << " / " << nbDimension << endl;

	  // Memorisation of the model - be careful, this model has a partition in double
	  push_model(precResidu);

	  curDimension          = _dimensionSet[noDimension];
	  dimensionMinPartition = curPartition->getDimensionMin(noDimension);
	  dimensionMaxPartition = curPartition->getDimensionMax(noDimension);

	  if (_useuniformcutting)     _cutFcn->Uniform_Cut(curPartition, noDimension, curDimension->getNbDiscretisations());
	  if (_usedistributedcutting) _cutFcn->Distributed_Cut(curPartition, noDimension, curDimension->getNbDiscretisations());

	  cout << endl;

 	  for(noSeparation = 0; noSeparation < _cutFcn->getNbSeparations(); noSeparation ++)
 	    {
	      cout << UP << "LL_Lolimot_addTheBestPartition() : test of the cutting ";
	      cout << noSeparation + 1 << " / " << _cutFcn->getNbSeparations() << endl;

 	      valeurSeparation = _cutFcn->getStep(noSeparation);

	      // We perform the cutting in two parts at the position given by valeurSeparation
	      curPartition->setDimensionMax(noDimension, valeurSeparation);
	      newPartition->setDimensionMin(noDimension, valeurSeparation);
	      
	      locPartitionWithBiggestResidu = NULL;

	      // we evaluate the cutting because it's valid

	      if (updateCoefficients(residu, locPartitionWithBiggestResidu)==false)
		{
		  cerr << "** ERROR LL_Lolimot::addTheBestPartition() : during the computation of the models." << endl;
		  return 1;
		} /* End If */

	      residuPourcent = residu/100.0;

	      if (residu < bestResidu) 
		{
		  bestResidu                 = residu;
		  bestResiduPourcent         = bestResidu/100.0;
		  bestPartition              = curPartition;
		  bestPrecResidu             = precResidu;
		  bestDimension              = noDimension;
		  bestValeurSeparation       = valeurSeparation; // value of the separation
		  partitionWithBiggestResidu = locPartitionWithBiggestResidu;

		  push_good_model(bestResidu);
		} /* End If */ // if this residual is more interesting than the best solution found so far ...
	    } /* End For noSeparation */

	  // We get back the old model
	  pop_model(residu);

	  // We update curPartition and newPartition
 	  curPartition->setDimensionMax(noDimension, dimensionMaxPartition);
 	  newPartition->setDimensionMin(noDimension, dimensionMinPartition);
	} /* End For noDimension */

      // If we have given a partition to separate, we quit
      if (partitionToSeparate) 	break;
    } /* End For noPartition */

  if (partitionWithBiggestResidu==NULL)
    {
      cerr << "LL_Lolimot::addTheBestPartition() : no cutting possible" << endl;

      // We get back the old model
      pop_model(residu);
            
      // We delete the partition which just has been added
      newPartition = _partitionSet[_partitionSet.size()-1];

      _partitionSet.pop_back();

      delete newPartition;

      return -1;
    } /* End If */

  _relativeGap = fabs((precResidu - bestResidu) / precResidu);

  // display the current solution 

  cerr << endl << "----------------------------------------" << endl;
  cerr << "LL : Nb of partitions : " << _partitionSet.size() << endl;
  cerr << "LL : precResidu              = " << precResidu << endl;
  cerr << "LL : bestResidu              = " << bestResidu << endl;
  cerr << "LL : _initialResidu (global) = " << _initialResidu << endl;
  cerr << "LL : relative gap            = " << _relativeGap * 100.0 << "  % (- ";
  cerr << precResidu - bestResidu << ")" << endl;
  
  // If we are not in continuous mode, we stop when the relative improvement is below a certain threshold
  if (!_opt_continue)
    {
      // stopping test on the percentage of improvement with respect to the initial residual
      if (_relativeGap < getResiduGapPercentage())
	{
	  cerr << "** LL_Lolimot::addTheBestPartition() : The test has stopped : percentage of improvement with respect to the initial residual" << endl;
	  
	  // We get back the old model
	  pop_good_model(residu);

	  return -1;
	} /* End If */
    } /* End If */

  // If the partition to cut has too few data, we haven(t been able to perform a cutting
  // So, partitionWithBiggestResidu is equal to NULL. We have no right to do the cutting and go to the next partition
  // If we mark the partition as "uncuttable", we can do the first part of the if
  // Else, we don't have to consider the worsening cuts because they can induce the loop

  cerr << "LL_Lolimot::addTheBestPartition() : good cutting - ";
  cerr << "bestPrecResidu = " << bestPrecResidu << " / bestResidu = " << bestResidu << endl;

  // We get back the old model
  pop_good_model(residu);
	        
  // Control test : we must have residu == bestResidu !!
  if (fabs((residu - bestResidu)) > 1e-4)
    {
      cerr << "** ERROR LL_Lolimot::addTheBestPartition() : incoherent value of the residual (";
      cerr << residu << ", " << bestResidu << ") !!." << endl;

      return 1;
    } /* End If */

  // Display the list of partitions
  
  cout << "LL_Lolimot::addTheBestPartition - List of partitions" << endl;
  
  nbPartition = _partitionSet.size();
  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = _partitionSet[noPartition];
      
      cerr << "n " << noPartition+1 << " " << curPartition->printDimensions() << endl;
    } /* End For */
  
  return 0;
}
  
void LL_Lolimot::clear_stack_model()
{
  unsigned int i;

  // Clearing the stack
  if (stack_dimensionSet.size()!=0)
    {
      for(i=0; i<stack_dimensionSet.size(); i++)
	{
	  delete stack_dimensionSet[i];
	} /* End For */
      
      for(i=0; i<stack_mesureSet.size(); i++)
	{
	  delete stack_mesureSet[i];
	} /* End For */
      
      for(i=0; i<stack_partitionSet.size(); i++)
	{
	  delete stack_partitionSet[i];
	} /* End For */
      
      stack_dimensionSet.resize(0);
      stack_mesureSet.resize(0);
      stack_partitionSet.resize(0);
    } /* End If */
}

void LL_Lolimot::push_model(float residuval)
{
  unsigned int i;

  clear_stack_model();

  stack_dimensionSet.resize(_dimensionSet.size(), NULL);
  stack_mesureSet.resize(_mesureSet.size(), NULL);
  stack_partitionSet.resize(_partitionSet.size(), NULL);

  // Pushing the new model
  for(i=0; i<_dimensionSet.size(); i++)
    {
      stack_dimensionSet[i]    = new LL_Dimension;
      (*stack_dimensionSet[i]) = (*_dimensionSet[i]);
    } /* End For */

  for(i=0; i<_mesureSet.size(); i++)
    {
      stack_mesureSet[i]    = new LL_Mesure;
      (*stack_mesureSet[i]) = (*_mesureSet[i]);
    } /* End For */

  for(i=0; i<_partitionSet.size(); i++)
    {
      stack_partitionSet[i]    = new LL_Partition;
      (*stack_partitionSet[i]) = (*_partitionSet[i]);
    } /* End For */

  stack_residu = residuval;
}

void LL_Lolimot::pop_model(float & residuval)
{
  unsigned int i;

  if (stack_dimensionSet.size()!=0)
    {
      for(i=0; i<stack_dimensionSet.size(); i++)
	{
	  (*_dimensionSet[i]) = (*stack_dimensionSet[i]);
	} /* End For */
      
      for(i=0; i<stack_mesureSet.size(); i++)
	{
	  (*_mesureSet[i]) = (*stack_mesureSet[i]);
	} /* End For */
      
      for(i=0; i<stack_partitionSet.size(); i++)
	{
	  (*_partitionSet[i]) = (*stack_partitionSet[i]);
	} /* End For */
      
      residuval = stack_residu;
    } /* End If */
  else
    {
      cout << "LL_Lolimot:: pop an empty stack - pop_model" << endl;
      cout << "LL_Lolimot:: do nothing" << endl;
      residuval = numeric_limits<float>::max();
    } /* End Else */
}

void LL_Lolimot::clear_good_stack()
{
  unsigned int i;

  // Clearing the stack
  if (stack_good_dimensionSet.size()!=0)
    {
      for(i=0; i<stack_good_dimensionSet.size(); i++)
	{
	  delete stack_good_dimensionSet[i];
	} /* End For */
      
      for(i=0; i<stack_good_mesureSet.size(); i++)
	{
	  delete stack_good_mesureSet[i];
	} /* End For */
      
      for(i=0; i<stack_good_partitionSet.size(); i++)
	{
	  delete stack_good_partitionSet[i];
	} /* End For */
    } /* End Else */

  stack_good_dimensionSet.resize(0);
  stack_good_mesureSet.resize(0);
  stack_good_partitionSet.resize(0);
}

void LL_Lolimot::push_good_model(float residuval)
{
  unsigned int i;

  clear_good_stack();

  // Pushing the best model

  stack_good_dimensionSet.resize(_dimensionSet.size(), NULL);
  stack_good_mesureSet.resize(_mesureSet.size(), NULL);
  stack_good_partitionSet.resize(_partitionSet.size(), NULL);

  for(i=0; i<_dimensionSet.size(); i++)
    {
      stack_good_dimensionSet[i]    = new LL_Dimension;
      (*stack_good_dimensionSet[i]) = (*_dimensionSet[i]);
    } /* End For */

  for(i=0; i<_mesureSet.size(); i++)
    {
      stack_good_mesureSet[i]    = new LL_Mesure;
      (*stack_good_mesureSet[i]) = (*_mesureSet[i]);
    } /* End For */

  for(i=0; i<_partitionSet.size(); i++)
    {
      stack_good_partitionSet[i]    = new LL_Partition;
      (*stack_good_partitionSet[i]) = (*_partitionSet[i]);
    } /* End For */

  stack_good_residu = residuval;
}

void LL_Lolimot::pop_good_model(float & residuval)
{
  unsigned int i;

  if (stack_good_dimensionSet.size()==0)
    {
      cout << "LL_Lolimot:: pop an empty stack - pop_good_model()" << endl;

      pop_model(residuval);
    } /* End If */
  else
    {
      for(i=0; i<stack_good_dimensionSet.size(); i++)
	{
	  (*_dimensionSet[i]) = (*stack_good_dimensionSet[i]);
	} /* End For */
      
      for(i=0; i<stack_good_mesureSet.size(); i++)
	{
	  (*_mesureSet[i]) = (*stack_good_mesureSet[i]);
	} /* End For */
      
      for(i=0; i<stack_good_partitionSet.size(); i++)
	{
	  (*_partitionSet[i]) = (*stack_good_partitionSet[i]);
	} /* End For */
      
      residuval = stack_good_residu;
    } /* End Else */
}

bool LL_Lolimot::step_optimise(bool separatePartitionWithBiggestResidu)
{
  float OldGapPercentage = 0.0;
  bool  Result = false;

  OldGapPercentage = getResiduGapPercentage();

  setNbMaxPartitions(getNbMaxPartitions()+1);
  setResiduGapPercentage(0.0);
  setContinue(true);

  Result = optimise(separatePartitionWithBiggestResidu);

  setContinue(false);
  setResiduGapPercentage(OldGapPercentage);
  
  return Result;
}

void LL_Lolimot::setInhibe(unsigned int Dimension, bool Value)
{
  unsigned int i;

  for(i=0; i<getPartitionSet().size(); i++)
    {
      getPartition(i)->setInhibe(Dimension, Value);
    } /* End For */
}

bool LL_Lolimot::getInhibe(unsigned int Dimension) const
{
  return getPartition(0)->getInhibe(Dimension);
}

void LL_Lolimot::setContinue(bool opt_continue)
{
  _opt_continue = opt_continue;
}

bool LL_Lolimot::getContinue() const
{
  return _opt_continue;
}

bool LL_Lolimot::optimise(bool separatePartitionWithBiggestResidu)
{
  LL_Partition * partitionToSeparate        = NULL;
  LL_Partition * partitionWithBiggestResidu = NULL;
  LL_Partition * newPartition = NULL;
  LL_Dimension * curDimension = NULL;
  unsigned int   noDimension;
  float          residu         = 0.0;
  float          residuPourcent = 0.0;

  clearAnyPartitionAvailable();

  _cutFcn->addMesures(_mesureSet);

  // When we are not in continuous mode, a partition already exists
  // We don't have to create a new partition except if the max number of partition is equal to 1
  
  if ((!_opt_continue)||(_nbMaxPartitions==1))
    {
      // If there are already partition, we stop
      if (getPartitionSet().size() > 0)
	{
	  cerr << "** ERROR LL_Lolimot::optimise() : an optimization has already be done." << endl;
	  return false;
	} /* End If */

      // we create a first partition for initialization
      newPartition = addPartition();

      // we initialize this partition with the bounds of the dimensions
      for(noDimension = 0; noDimension < _dimensionSet.size(); noDimension ++)
	{
	  curDimension = _dimensionSet[noDimension];
	  newPartition->setDimensionMin(noDimension, curDimension->getMin());
	  newPartition->setDimensionMax(noDimension, curDimension->getMax());
	} /* End For */

      newPartition->updateEcartType(this);
    } /* End If */

  // we solve the problem with 1 partition to see what is the value of the residual
  // We recompute all the residuals for all the partitions
  // In partitionWithBiggestResidu, we highligth the partition which has the greatest value of residual
  // It's the partition we have to cut
  // residu contains the value of the greatest residual

  if (updateCoefficients(residu, partitionWithBiggestResidu)==false)
    {
      cerr << "** ERROR LL_Lolimot::addTheBestPartition() : during the computation of the models." << endl;
      return 1;
    } /* End If */

  _residu = computeGlobalResidual();
  
  residuPourcent = _residu/100.0;
    
  _initialResidu = _residu;
  
  cerr << "LL_Lolimot::optimise - Initial residual before optimization : " << getInitialResidu() << endl;

  if (_nbMaxPartitions>1)
    {
      // while the stopping conditions are not satisfied, we add a partition
      bool swCreated = true;
      do
	{
	  if (separatePartitionWithBiggestResidu && (partitionWithBiggestResidu!=NULL))
	    {
	      partitionToSeparate = partitionWithBiggestResidu;
	      cerr << "LL_Lolimot::optimise - partition to cut : ";
	      cerr << partitionWithBiggestResidu->printDimensions() << endl;
	    } /* End If */
	  else
	    partitionToSeparate = NULL;
	  
	  int code;
	  
	  code = addTheBestPartition(partitionToSeparate, partitionWithBiggestResidu);
	  
	  _residu = computeGlobalResidual();
	  
	  if (code > 0)
	    {
	      cerr << "** ERROR LL_Lolimot::optimise() : in addTheBestPartition()." << endl;
	      return false;
	    } /* End If */
	  
	  if (code < 0)
	    {
	      swCreated = false;
	    } /* End If */
	  
	  partitionToSeparate = partitionWithBiggestResidu;
	} /* End Do */
      while(((_partitionSet.size() < (unsigned int)getNbMaxPartitions()) && (swCreated == true)));
    } /* End If */

  clear_stack_model();
  clear_good_stack();
  
  return true;
}

float LL_Lolimot::computeGlobalResidual()
{
  LL_Mesure             * curMesure = NULL;
  LL_Partition          * curPartition = NULL;
  int                     noMesure,    nbMesures   = _mesureSet.size();
  int                     noPartition, nbPartition = _partitionSet.size();
  vector<vector<float> > phiEP;
  float                  Residu = 0.0;
  vector<float>          MeasureList, EstimList;

  // we update the steady deviation of each partition-dimension
  for(noPartition=0; noPartition < nbPartition; noPartition++)
    {
      curPartition = _partitionSet[noPartition];
      curPartition->updateEcartType(this);
    } /* End For */

  phiEP.resize(nbMesures, vector<float>(nbPartition, 0.0));
  MeasureList.resize(nbMesures, 0.0);
  EstimList.resize(nbMesures, 0.0);

  // for each experiment, we compute phi[e]
  for(noMesure = 0; noMesure < nbMesures; noMesure ++)
    {
      curMesure = _mesureSet[noMesure];
      calculePhi(curMesure->getDimensionValueSet(), phiEP[noMesure]);
    } /* End For */

  for(noMesure = 0; noMesure < nbMesures; noMesure ++)
    {
      curMesure             = _mesureSet[noMesure];
      EstimList[noMesure]   = calculeF(curMesure->getDimensionValueSet(), phiEP[noMesure]);
      MeasureList[noMesure] = curMesure->getValue();
    } /* End For */

  // Updating of the residual. We store in residu the value of the residual computed on the whole domain

  Residu = (*_FctResidu)(MeasureList, EstimList, _paramList);

  return Residu;
}

void LL_Lolimot::setUseTransformLog(bool UseLogTransform)
{
  _useTransformLog = UseLogTransform;
}

bool LL_Lolimot::getUseTransformLog() const
{
  return _useTransformLog;
}

void LL_Lolimot::setTransformLogAlpha(float alpha)
{
  _TransformLogAlpha = alpha;
}

float LL_Lolimot::getTransformLogAlpha() const
{
  return _TransformLogAlpha;
}

void LL_Lolimot::setTransformLogMin(float Min)
{
  _TransformLogMin = Min;
}

float LL_Lolimot::getTransformLogMin() const
{
  return _TransformLogMin;
}

void LL_Lolimot::setTransformLogMax(float Max)
{
  _TransformLogMax = Max;
}

float LL_Lolimot::getTransformLogMax() const
{
  return _TransformLogMax;
}

void LL_Lolimot::setTransformLogEps(float Eps)
{
  _TransformLogEps = Eps;
}

float LL_Lolimot::getTransformLogEps() const
{
  return _TransformLogEps;
}

bool LL_Lolimot::exportAllPartitions(string path, string Filename)
{
  int            noPartition, nbPartition = getPartitionSet().size();
  LL_Partition * curPartition = NULL;
  LL_Dimension * curDimension = NULL;
  int            noDimension, nbDimension = getDimensionSet().size();
  ofstream       F_Lolimot;
  string         strEcartType;
  string         strTmp;
  unsigned int   i;

  strTmp = path + "/" + Filename;
  F_Lolimot.open(strTmp.c_str());

  if (!F_Lolimot.is_open())
    {
      cerr << "** Error LL_Lolimot::exportAllPartitions() : can't create the file ";
      cerr << strTmp << "." << endl;
      return false;
    } /* End If */

  for(noDimension = 0; noDimension < nbDimension; noDimension ++)
    {
      curDimension = _dimensionSet[noDimension];
      F_Lolimot << "List of Partitions for Dimension " << noDimension << " (" << curDimension->getName();
      F_Lolimot << ")" << endl;

      for(noPartition = 0; noPartition < nbPartition; noPartition ++)
	{
	  curPartition = _partitionSet[noPartition];

	  F_Lolimot << "Partition " << noPartition << " :";
	  F_Lolimot << " Min = " << curPartition->getDimensionMin(noDimension) << " -";
	  F_Lolimot << " Max = " << curPartition->getDimensionMax(noDimension) << endl;
	} /* End For */
    } /* End For */ 

  F_Lolimot << endl << "List of Models for each Partition: " << endl;

  for(noPartition = 0; noPartition < nbPartition; noPartition ++)
    {
      curPartition = _partitionSet[noPartition];
      
      F_Lolimot << "Partition " << noPartition << " : " << endl;
      F_Lolimot << curPartition->getCoeff0() << " ";
      for(i=0; i<curPartition->getCoeffSet().size(); i++)
	{
	  F_Lolimot << curPartition->getCoeff(i) << " ";
	} /* End For */
      F_Lolimot << endl;
    } /* End For */

 F_Lolimot.close();

  return true;
}

void LL_Lolimot::setDontCutVar(string Name, bool Value)
{
  int Index = -1;

  do
    {
      Index++;
    }
  while((Name!=_dimensionSet[(unsigned int)Index]->getName())&&((unsigned int)Index<_dimensionSet.size()));

  if ((unsigned int)Index<_dimensionSet.size()) ListOfVarToDontCut[(unsigned int)Index] = Value;
}

bool LL_Lolimot::getDontCutVar(string Name) const
{
  int  Index = -1;
  bool Result;

  do
    {
      Index++;
    }
  while((Name!=_dimensionSet[(unsigned int)Index]->getName())&&((unsigned int)Index<_dimensionSet.size()));

  if ((unsigned int)Index<_dimensionSet.size()) Result = ListOfVarToDontCut[Index];
  else                                          Result = false;

  return Result;
}

void LL_Lolimot::setDontCutVar(unsigned int Index, bool Value)
{
  ListOfVarToDontCut[Index] = Value;
}

bool LL_Lolimot::getDontCutVar(unsigned int Index) const
{
  return ListOfVarToDontCut[Index];
}

void LL_Lolimot::clearDontCutList()
{
  ListOfVarToDontCut.resize(_dimensionSet.size(), false);
}

void LL_Lolimot::setExportType_double(bool Var)
{
  _exportTypeDouble = Var;
}

bool LL_Lolimot::getExportType_double() const
{
  return _exportTypeDouble;
}

void LL_Lolimot::setExportType_float(bool Var)
{
  _exportTypeDouble = !Var;
}

bool LL_Lolimot::getExportType_float() const
{
  return !_exportTypeDouble;
}
