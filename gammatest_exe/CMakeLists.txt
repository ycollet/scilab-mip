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

option(GT_ENABLE_WARNINGS "Enable more warnings from the compiler" OFF)

include(CheckCXXCompilerFlag)

if (GT_ENABLE_WARNINGS)
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter -Wduplicated-cond -Wduplicated-branches -Wlogical-op") # -Wunused
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wrestrict -Wnull-dereference -Wuseless-cast") # -Wold-style-cast
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wdelete-non-virtual-dtor -Wdouble-promotion -Wshadow -Wformat=2") 
    check_cxx_compiler_flag(-Wextra HAVE_EXTRA_FLAG)
    if (HAVE_EXTRA_FLAG)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
    endif ()
    check_cxx_compiler_flag(-Wc++11-compat HAVE_CXX11_COMPAT_FLAG)
    if (HAVE_CXX11_COMPAT_FLAG)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wc++11-compat")
    endif ()
  elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall") # Trigger all warnings for clang
    if (KL_ENABLE_WARNINGS)
      check_cxx_compiler_flag(-Weverything HAVE_EVERYTHING_FLAG)
      if (HAVE_EVERYTHING_FLAG)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weverything")
      endif ()
    endif ()
  endif ()
endif ()

include_directories(src)

set(GAMMATEST_SRCS src/Arguments.cpp
                   src/Filter.cpp
                   src/GammaTest.cpp
		   src/GammaTest_Estimator.cpp
		   src/kdtree_2.cpp
		   src/kdtree_static.cpp
		   src/MainGT.cpp)

add_executable(gammatest ${GAMMATEST_SRCS})
