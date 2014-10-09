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

#include "LL_Mesure.h"

///////////////
// LL_Mesure //
///////////////

void  LL_Mesure::setValue(float value)
{
  _value = value;
}

float LL_Mesure::getValue() const
{
  return _value;
}

std::vector<float> LL_Mesure::getDimensionValueSet() const
{
  return _dimensionValueSet;
}

float LL_Mesure::getDimensionValue(int i) const
{
  return _dimensionValueSet[i];
}

LL_Mesure::LL_Mesure(float value, int nbDimensions, float minValueExtrapolation, float maxValueExtrapolation)
  : _value(value), _minValueExtrapolation(minValueExtrapolation), _maxValueExtrapolation(maxValueExtrapolation)
{
  _dimensionValueSet.resize(nbDimensions, 0.0);

  _valueExtrapolation = 0.0;
}

LL_Mesure::LL_Mesure() : _value(-1), _minValueExtrapolation(-1), _maxValueExtrapolation(-1)
{
}

LL_Mesure::LL_Mesure(const LL_Mesure & Var)
{
  _value                 = Var._value;
  _minValueExtrapolation = Var._minValueExtrapolation;
  _maxValueExtrapolation = Var._maxValueExtrapolation;
  _valueExtrapolation    = Var._valueExtrapolation;
  _dimensionValueSet     = Var._dimensionValueSet;
}

LL_Mesure::~LL_Mesure()
{
}


void LL_Mesure::setDimensionValue(int i, float value)
{
  _dimensionValueSet[i] = value;
}

LL_Mesure & LL_Mesure::operator=(LL_Mesure & right)
{
  _value                 = right._value;
  _minValueExtrapolation = right._minValueExtrapolation;
  _maxValueExtrapolation = right._maxValueExtrapolation;
  _valueExtrapolation    = right._valueExtrapolation;
  _dimensionValueSet     = right._dimensionValueSet;

  return *this;
}
