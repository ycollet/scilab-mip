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

#ifndef LL_CUT_H
#define LL_CUT_H

#include <vector>

#include "LL_Mesure.h"
#include "LL_Partition.h"

class LL_Cut
{
 public:
  LL_Cut();
  void         Uniform_Cut(LL_Partition * curPartition, unsigned int noDimension, unsigned int NbCut);
  void         Distributed_Cut(LL_Partition * curPartition, unsigned int noDimension, unsigned int NbCut);
  float        getStep(unsigned int noStep);
  unsigned int getNbSeparations();
  void         addMesures(std::vector<LL_Mesure *> & Mesures);
  std::vector<LL_Mesure *> & getMesure();
 protected:
  std::vector<float>       _separationSet;
  std::vector<LL_Mesure *> _mesureSet;
};
#endif
