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

#ifndef LL_PARTITION_H
#define LL_PARTITION_H

///////////////////////
//LL_Partition class //
///////////////////////

// A partition corresponds to an hyper-rectangle in the domain
// It's the result of the optimization
// It's defined by the min and max with respect to each dimension
// It's also defined by the coefficients (w) of the regression

#include <vector>
#include <string>

class LL_Lolimot;
class LL_Mesure;

class LL_Partition
{
public:
   LL_Partition();
   LL_Partition(int nbDimensions);
   LL_Partition(const LL_Partition &);
  ~LL_Partition();
  
  void          setNbDimension(int nbDimensions);

  float         getCoeff0() const {return _coeff0;}
  void          setCoeff0(float coeff0) {_coeff0 = coeff0;}

  void          setFreezeCoeff0(bool Value) {_freezeCoeff0 = Value;}
  bool          getFreezeCoeff0() const {return _freezeCoeff0;}
  
  std::vector<float> getCoeffSet()   const {return _coeffSet;}
  float              getCoeff(int i) const {return _coeffSet[i];}
  void               setCoeff(int i, float value);

  std::vector<bool>  getFreezeCoeffSet() const {return _freezeCoeffSet;}
  void               setFreezeCoeff(int i, bool Value) {_freezeCoeffSet[i] = Value;}
  bool               getFreezeCoeff(int i) const {return _freezeCoeffSet[i];}

  std::vector<float> getDimensionMinSet()   const {return _dimensionMinSet;}
  float              getDimensionMin(int i) const {return _dimensionMinSet[i];}
  void               setDimensionMin(int i, float value);
  
  std::vector<float> getDimensionMaxSet()   const {return _dimensionMaxSet;}
  float              getDimensionMax(int i) const {return _dimensionMaxSet[i];}
  void               setDimensionMax(int i, float value);
  
  float         getResidu() const {return _residu;}
  void          setResidu(float residu) {_residu = residu;}

  bool          getInhibe(int i) const {return _inhibeSet[i];}
  void          setInhibe(int i, bool value);

  // this method copies the dimensions of curPartition in this partition
  bool          setDimensions(const LL_Partition * curPartition, LL_Lolimot * Lolimot);
  std::string   printDimensions();

  // managing the steady deviation for each partition
  std::vector<float> getEcartTypeSet()   const {return _ecartTypeSet;}
  float              getEcartType(int i) const {return _ecartTypeSet[i];}
  void               updateEcartType(const LL_Lolimot * curLolimot);
  
  // centre(i) = (min(i) + max(i)) / 2
  float         getCentre(int i) const {return (getDimensionMin(i) + getDimensionMax(i)) / 2.0;}
  
  // compute z for given values for each dimension
  // z = e(-1/2 (sum[d] (x(d) - centre(d))^2 / std(d)^2 ))
  // d belongs to the set of dimensions
  float         calculeZ(const std::vector<float> & x) const;
  
  // compute the resulting value for given values for each dimension
  // using the coefficients of the regression
  // val = coeff0 + sum[d] (coeff(d) * x(d))
  // d belong to the set of dimension
  float         calculeValeurSelonRegression(const std::vector<float> & x) const;
  bool          isMesureInPartition(const LL_Mesure * curMesure) const;
  bool          DoNotCut() const {return _DoNotCut;} 
  void          SetDoNotCut(bool DoNotCut) {_DoNotCut = DoNotCut;}
  void          setName(std::string Name) {_Name = Name;}
  std::string   getName() const {return _Name;}
  std::string   getRandName(std::string Prefix, unsigned int NbRandChar);
  LL_Partition & operator=(LL_Partition & right);
 private:
  std::vector<float> _dimensionMinSet;
  std::vector<float> _dimensionMaxSet;
  std::vector<float> _coeffSet;
  std::vector<bool>  _freezeCoeffSet;
  std::vector<bool>  _inhibeSet;
  float              _coeff0;
  bool               _freezeCoeff0;
  float              _residu;
  bool               _DoNotCut;
  std::string        _Name;
  // steady deviation for each dimension
  // EcartType(d) = std(d) * (max(d) - min(d))
  // ecartTypeDefault(d) is defined in LL_Dimension
  std::vector<float> _ecartTypeSet;
  std::string        strtmp;
}; // LL_Partition
#endif
