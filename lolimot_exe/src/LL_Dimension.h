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

#ifndef LL_DIMENSION_H
#define LL_DIMENSION_H

#include <string>

/////////////////////////
// classe LL_Dimension //
/////////////////////////

// A dimension corresponds to one of the axis in the space

class LL_Dimension
{
 public:
   LL_Dimension();
   LL_Dimension(std::string name, float min, float max, float ecartType, int nbDiscretisations);
   LL_Dimension(const LL_Dimension &);
  ~LL_Dimension();
  void   setMin(float min) {_min = min;}
  float  getMin() const;
  void   setMax(float max) {_max = max;}
  float  getMax() const;
  void   setEcartType(float ecartType) {_ecartType = ecartType;}
  float  getEcartType() const;
  void   setNbDiscretisations(int nbDiscretisations) {_nbDiscretisations = nbDiscretisations;}
  int    getNbDiscretisations() const;
  void   setName(std::string name);
  std::string getName() const;
  void   setRecursive(bool RecursiveVar);
  bool   getRecursive() const;
  LL_Dimension & operator=(LL_Dimension & right);
 private:
  std::string _name;
  float       _min;
  float       _max;
  float       _ecartType;
  int         _nbDiscretisations;
  bool        _isRecursive;
}; // LL_Dimension
#endif
