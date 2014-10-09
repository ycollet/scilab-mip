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

#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <vector>
#include <string>

#include "LL_Lolimot.h"

void Analysis_Response_NL(LL_Lolimot * Lolimot, std::vector<float> & Result);
void Analysis_Response(LL_Lolimot * Lolimot, std::vector<float> & Result,
		       std::vector<float> & InputMin, std::vector<float> & InputMax);
void Analysis_List_DistribOfCut(LL_Lolimot * Lolimot, std::string Filename);
#endif
