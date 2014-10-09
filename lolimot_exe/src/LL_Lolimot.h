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

///////////////////////////////////////////////////////////////////
// Description      : Modelisation using Lolimot functions       //
//                                                               //
// Date of creation :      11/12/2002                            //
//                                                               //
// Author           : Laurent PAJOU                              //
// Comment          :                           Ste EURODECISION //
// Modified by      : Yann COLLETTE (Renault)                    //
///////////////////////////////////////////////////////////////////

#ifndef LL_LOLIMOT_H
#define LL_LOLIMOT_H

#include <math.h>

#include <vector>
#include <string>

#include "DefaultParam.h"

class LL_Dimension;
class LL_Mesure;
class LL_Partition;
class LL_Cut;

typedef float (*PtrToFctResidu)(std::vector<float> &, std::vector<float> &, std::vector<float> &);

///////////////////
// Lolimot Class //
///////////////////

class LL_Lolimot
{
 public:
  LL_Lolimot();
  LL_Lolimot(std::string name, int nbMaxPartitions, float residuGapPercentage);
  virtual ~LL_Lolimot();
  
  LL_Lolimot & operator=(const LL_Lolimot & right);
  
  int    getNbMaxPartitions() const;
  float  getResiduGapPercentage() const;
  
  void   setNbMaxPartitions(int nb);
  void   setResiduGapPercentage(float res);
  
  float  getResidu() const;
  float  getInitialResidu() const;
  void   setInitialResidu(float InitRes);

  void   setSigma(float SigmaVar);
  float  getSigma() const;

  // managing the dimensions
  std::vector<LL_Dimension*> getDimensionSet() const;
  LL_Dimension *             getDimension(int i) const;
  LL_Dimension *             addDimension(std::string name, float min, float max, float ecartType, int nbDiscretisations);
  int                        cleanDimensions();
  
  // managing the measurements
  std::vector<LL_Mesure*> getMesureSet() const;
  LL_Mesure *             getMesure(int i) const;
  LL_Mesure *             addMesure(float value, float minValueExtrapolation = -1e39, float maxValueExtrapolation = 1e39);
  int                     cleanMesures();
  
  // managing the partitions
  std::vector<LL_Partition*> getPartitionSet() const;
  LL_Partition *             getPartition(int i) const;
  LL_Partition *             addPartition();
  std::string                getRandName(std::string Prefix, unsigned int NbRandChar);
  int                        cleanPartitions();

  // managing the weights
  std::vector<float>    getWeightSet() const;
  float                 getWeight(int i) const;
  void                  addWeight(float Weight);
  void                  cleanWeight();
  void                  useDefaultWeight();

  void                  setMaxNbPointsPerPartitions(int MaxPtPerPart);
  int                   getMaxNbPointsPerPartitions() const;

  // managing the function which allows to computed the global residual
  void                  setNbParam(unsigned int nbParam);
  unsigned int          getNbParam() const;
  void                  setParam(unsigned int nbParam, float ParamValue);
  float                 getParam(unsigned int nbParam) const;
  void                  addFctResidu(PtrToFctResidu FctResidu);
  
  bool                  exportAllPartitions(std::string path, std::string Filename);

  // managing the function which allows to normalize the input / output
  void                  setUseTransformLog(bool UseLogTransform);
  bool                  getUseTransformLog() const;
  void                  setTransformLogAlpha(float alpha);
  float                 getTransformLogAlpha() const;
  void                  setTransformLogMin(float min);
  float                 getTransformLogMin() const;
  void                  setTransformLogMax(float max);
  float                 getTransformLogMax() const;
  void                  setTransformLogEps(float eps);
  float                 getTransformLogEps() const;

  // Computation of the estimation using data given for each dimension
  // f(x) = sum[p] (ValeurSelonRegression(p) * phi(p))
  // p belongs to the set of partitions
  // phi(p) = z(p) / sum[p2] z(p2)   (phi(p) is a probability distribution)
  float                 calculeF(const std::vector<float> & x) const;
  float                 calculeF_Untransformed(const std::vector<float> & x) const;
  // same thing, but we don't recompute phi(x)
  float                 calculeF(const std::vector<float> & x, const std::vector<float> & phi) const;
  float                 calculeF_Untransformed(const std::vector<float> & x, const std::vector<float> & phi) const;
  // compute phi(x)
  bool                  calculePhi(const std::vector<float> & x, std::vector<float> & phi) const;

  bool                  optimise(bool separatePartitionWithBiggestResidu);
  bool                  step_optimise(bool separatePartitionWithBiggestResidu);

  void                  push_model(float);
  void                  pop_model(float &);
  void                  clear_stack_model();

  void                  push_good_model(float);
  void                  pop_good_model(float &);
  void                  clear_good_stack();

  // compute the coefficients (and the associated residual) with respect to the partitions
  // (we use the lapack function dgelss for each partition)
  // (this method modifies the value of the coefficients for each partition)
  bool                  updatePartitions();
  bool                  updateCoefficients(float & residu, LL_Partition * & partitionWithBiggestResidu);
  bool                  calculeCoefficients(float & residu, LL_Partition * & partitionWithBiggestResidu);
  bool                  calculeCoefficients_Lasso(float & residu, LL_Partition * & partitionWithBiggestResidu);

  // this method computes the lolimot coefficients so as to minimize the gap between estimation and measurement
  // and so as to verify the bounds on these measures
  void                  setContinue(bool opt_continue);
  bool                  getContinue() const;

  void                  setMembershipThreshold(float Value);
  float                 getMembershipThreshold() const;

  void                  setInhibe(unsigned int Dimension, bool Value);
  bool                  getInhibe(unsigned int Dimension) const;

  void                  setDontCutVar(std::string Name, bool Value);
  bool                  getDontCutVar(std::string Name) const;
  void                  setDontCutVar(unsigned int Index, bool Value);
  bool                  getDontCutVar(unsigned int Index) const;
  void                  clearDontCutList();

  bool                  isThereAnyPartitionAvailable() const;
  void                  clearAnyPartitionAvailable();

  void                  setUseUniformCutting(bool Value);
  bool                  getUseUniformCutting() const;

  void                  setUseLasso(bool Value);
  bool                  getUseLasso() const;

  void                  setUseDistributedCutting(bool Value);
  bool                  getUseDistributedCutting() const;

  void                  setExportType_double(bool Var);
  bool                  getExportType_double() const;
  void                  setExportType_float(bool Var);
  bool                  getExportType_float() const;

  float                 computeGlobalResidual();

  void                  setCustomArguments(std::string Arguments);
  std::string           getCustomArguments() const;
  void                  setCustomReturn_Before(std::string Arguments);
  std::string           getCustomReturn_Before() const;
  void                  setCustomReturn_After(std::string Arguments);
  std::string           getCustomReturn_After() const;
 private:
  int                        _nbMaxPartitions;
  float                      _residuGapPercentage;

  std::vector<LL_Dimension*> _dimensionSet;
  std::vector<LL_Mesure*>    _mesureSet;
  std::vector<LL_Partition*> _partitionSet;
  std::vector<float>         _Weight;

  std::vector<LL_Dimension*> stack_dimensionSet;
  std::vector<LL_Mesure*>    stack_mesureSet;
  std::vector<LL_Partition*> stack_partitionSet;
  float                      stack_residu;

  std::vector<LL_Dimension*> stack_good_dimensionSet;
  std::vector<LL_Mesure*>    stack_good_mesureSet;
  std::vector<LL_Partition*> stack_good_partitionSet;
  float                      stack_good_residu;

  std::vector<float>     _paramList;
  PtrToFctResidu         _FctResidu;
  LL_Cut *               _cutFcn;  
  float                  _residu;
  float                  _initialResidu;
  float                  _relativeGap;
  // method which add the best partition (i.e. which reduce the most the global residual)
  int                    addTheBestPartition(LL_Partition * partitionToSeparate, LL_Partition * & partitionWithBiggestResidu);
  // returns 0 if OK
  //        -1 if no partition has been added (stopping criterion)
  //         1 if an error occured
  bool              _opt_continue;
  bool              _useuniformcutting;
  bool              _usedistributedcutting;
  bool              _anyPartitionAvailable;
  bool              _useTransformLog;
  bool              _useLasso;
  float             _TransformLogAlpha;
  float             _TransformLogMin;
  float             _TransformLogMax;
  float             _TransformLogEps;
  float             _MembershipThreshold;
  int               _MaxNbPointsPerPartitions;
  std::vector<bool> ListOfVarToDontCut;
  bool              _exportTypeDouble;
  std::string       _customArguments;
  std::string       _customReturn_Before;
  std::string       _customReturn_After;
};
#endif
