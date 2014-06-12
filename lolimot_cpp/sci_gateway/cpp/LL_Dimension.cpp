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

#include "LL_Dimension.h"

//////////////////
// LL_Dimension //
//////////////////

float LL_Dimension::getMin() const
{
  return _min;
}

float LL_Dimension::getMax() const
{
  return _max;
}

float LL_Dimension::getEcartType() const
{
  return _ecartType;
}

int LL_Dimension::getNbDiscretisations() const
{
  return _nbDiscretisations;
}

void LL_Dimension::setName(string name)
{
  _name = name;
}

string LL_Dimension::getName() const
{
  return _name;
}

LL_Dimension::LL_Dimension(string name, float min, float max, float ecartType, int nbDiscretisations)
  : _name(name), _min(min), _max(max), _ecartType(ecartType), _nbDiscretisations(nbDiscretisations) 
{
}

LL_Dimension::LL_Dimension() : _name("Not Initialized"), _min(-1), _max(-1), _ecartType(-1), _nbDiscretisations(-1)
{
}

LL_Dimension::LL_Dimension(const LL_Dimension & Var)
{
  _name              = Var._name;
  _min               = Var._min;
  _max               = Var._max;
  _ecartType         = Var._ecartType;
  _nbDiscretisations = Var._nbDiscretisations;
}

LL_Dimension::~LL_Dimension()
{
}

LL_Dimension & LL_Dimension::operator=(LL_Dimension & right)
{
  _name              = right._name;
  _min               = right._min;
  _max               = right._max;
  _ecartType         = right._ecartType;
  _nbDiscretisations = right._nbDiscretisations;

  return *this;
}

void LL_Dimension::setRecursive(bool RecursiveVar)
{
  _isRecursive = RecursiveVar;
}

bool LL_Dimension::getRecursive() const
{
  return _isRecursive;
}
