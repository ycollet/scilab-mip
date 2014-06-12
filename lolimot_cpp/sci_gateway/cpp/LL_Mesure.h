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

#ifndef LL_MESURE_H
#define LL_MESURE_H

#include "LL_Lolimot.h"
#include "LL_Partition.h"

/////////////////////
// LL_Mesure class //
/////////////////////

// A measure (experiment) est constituted by a value (for the measuremebt) and the value of each dimension
// (vector which has the same order than the list of dimensions in LL_Lolimot)

class LL_Mesure
{
public:
  LL_Mesure(float value, int nbDimensions, float minValueExtrapolation, float maxValueExtrapolation);
  LL_Mesure();
  LL_Mesure(const LL_Mesure &);
  virtual ~LL_Mesure();

  void  setValue(float value);
  float getValue() const;
  float getMinValueExtrapolation() const {return _minValueExtrapolation;}
  float getMaxValueExtrapolation() const {return _maxValueExtrapolation;}
  float getValueExtrapolation()    const {return _valueExtrapolation;}

  vector<float> getDimensionValueSet() const;
  float getDimensionValue(int i) const;
  void  setDimensionValue(int i, float value);

  LL_Mesure & operator=(LL_Mesure & right);
private:
  float _value;
  float _minValueExtrapolation;
  float _maxValueExtrapolation;
  float _valueExtrapolation;
  vector<float> _dimensionValueSet;
}; // LL_Mesure
#endif
