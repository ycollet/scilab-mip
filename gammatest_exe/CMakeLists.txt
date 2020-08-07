#=============================================================================
#  GammaTest
#  $Id:$
#
#  Copyright (C) 2020 by Yann Collette
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License version 2.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#=============================================================================

cmake_minimum_required(VERSION 2.8)

project(GammaTest CXX)

include_directories(.)

set(GAMMATEST_SRCS Arguments.cpp
                   Filter.cpp
                   GammaTest.cpp
		   GammaTest_Estimator.cpp
		   kdtree_2.cpp
		   kdtree_static.cpp
		   MainGT.cpp)

add_executable(gammatest ${GAMMATEST_SRCS})
