#=============================================================================
#  Lolimot
#  $Id:$
#
#  Copyright (C) 2013 by Yann Collette and others
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

project(Lolimot C CXX)

cmake_minimum_required(VERSION 2.8)

option(ENABLE_DEBUG_MESSAGES "Enable debug messages" OFF)
option(ENABLE_WARNINGS "Enable more warnings from the compiler" OFF)

Enable_Testing()

if (ENABLE_DEBUG_MESSAGES)
  add_definitions("-DDEBUG -D_DEBUG")
endif (ENABLE_DEBUG_MESSAGES)

include(CheckCXXCompilerFlag)

if (ENABLE_WARNINGS)
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

if (NOT WIN32)
  find_package(LAPACK REQUIRED)
endif ()

include_directories("${CMAKE_SOURCE_DIR}/src"
                    "${CMAKE_SOURCE_DIR}/src/Lasso"
                    "${CMAKE_BINARY_DIR}")

set(LOLIMOT_SRCS src/Analysis.cpp
                 src/DefaultParam.cpp
		 src/DisplayData.cpp
		 src/ExportC.cpp
		 src/ExportCommonC_Cpp.cpp
		 src/ExportCpp.cpp
		 src/ExportMatlab.cpp
		 src/ExportMinMaxMap.cpp
		 src/Filter.cpp
		 src/Learning.cpp
		 src/LL_Cut.cpp
		 src/LL_Dimension.cpp
		 src/LL_Lolimot.cpp
		 src/LL_Mesure.cpp
		 src/LL_Partition.cpp
		 src/Post.cpp
		 src/TrainLolimot2.cpp
		 src/TrainLolimotStruct.cpp
		 src/Validate.cpp)

set(LASSO_SRCS src/Lasso/lasso.c)

set(SMALL_LAPACK_SRCS src/SmallLapack/d_abs.c
                      src/SmallLapack/dbdsqr.c
		      src/SmallLapack/dcopy.c
		      src/SmallLapack/dgebd2.c
		      src/SmallLapack/dgebrd.c
		      src/SmallLapack/dgelq2.c
		      src/SmallLapack/dgelqf.c
		      src/SmallLapack/dgelss.c
		      src/SmallLapack/dgemm.c
		      src/SmallLapack/dgemv.c
		      src/SmallLapack/dgeqr2.c
		      src/SmallLapack/dgeqrf.c
		      src/SmallLapack/dger.c
		      src/SmallLapack/dlabad.c
		      src/SmallLapack/dlabrd.c
		      src/SmallLapack/dlacpy.c
		      src/SmallLapack/dlamch.c
		      src/SmallLapack/dlange.c
		      src/SmallLapack/dlapy2.c
		      src/SmallLapack/dlarfb.c
		      src/SmallLapack/dlarf.c
		      src/SmallLapack/dlarfg.c
		      src/SmallLapack/dlarft.c
		      src/SmallLapack/dlartg.c
		      src/SmallLapack/dlas2.c
		      src/SmallLapack/dlascl.c
		      src/SmallLapack/dlaset.c
		      src/SmallLapack/dlasq1.c
		      src/SmallLapack/dlasq2.c
		      src/SmallLapack/dlasq3.c
		      src/SmallLapack/dlasq4.c
		      src/SmallLapack/dlasq5.c
		      src/SmallLapack/dlasq6.c
		      src/SmallLapack/dlasr.c
		      src/SmallLapack/dlasrt.c
		      src/SmallLapack/dlassq.c
		      src/SmallLapack/dlasv2.c
		      src/SmallLapack/d_lg10.c
		      src/SmallLapack/d_lg.c
		      src/SmallLapack/dnrm2.c
		      src/SmallLapack/dorg2r.c
		      src/SmallLapack/dorgbr.c
		      src/SmallLapack/dorgl2.c
		      src/SmallLapack/dorglq.c
		      src/SmallLapack/dorgqr.c
		      src/SmallLapack/dorm2r.c
		      src/SmallLapack/dormbr.c
		      src/SmallLapack/dorml2.c
		      src/SmallLapack/dormlq.c
		      src/SmallLapack/dormqr.c
		      src/SmallLapack/drot.c
		      src/SmallLapack/drscl.c
		      src/SmallLapack/dscal.c
		      src/SmallLapack/d_sign.c
		      src/SmallLapack/d_sqrt.c
		      src/SmallLapack/dswap.c
		      src/SmallLapack/dtrmm.c
		      src/SmallLapack/dtrmv.c
		      src/SmallLapack/exit_.c
		      src/SmallLapack/F77_aloc.c
		      src/SmallLapack/ieeeck.c
		      src/SmallLapack/ilaenv.c
		      src/SmallLapack/lsame.c
		      src/SmallLapack/pow_dd.c
		      src/SmallLapack/pow_di.c
		      src/SmallLapack/s_cat.c
		      src/SmallLapack/s_cmp.c
		      src/SmallLapack/s_copy.c
		      src/SmallLapack/xerbla.c)

add_library(lasso STATIC ${LASSO_SRCS})

if (WIN32)
  add_library(small_lapack STATIC ${SMALLLAPACK_SRCS})
endif ()

add_executable(lolimot ${LOLIMOT_SRCS})
target_link_libraries(lolimot lasso)
if (WIN32)
  target_link_libraries(lolimot small_lapack)
else ()
  target_link_libraries(lolimot ${LAPACK_LIBRARIES})
endif ()

add_test(NAME lolimot_cone_noise 
         COMMAND $<TARGET_FILE:lolimot> ${CMAKE_SOURCE_DIR}/tests/Modele/Modele_Test_Cone_Noise.txt
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests)
set_tests_properties(lolimot_cone_noise PROPERTIES TIMEOUT 60)

add_test(NAME lolimot_sinus_noise 
         COMMAND $<TARGET_FILE:lolimot> ${CMAKE_SOURCE_DIR}/tests/Modele/Modele_Test_Sinus_Noise.txt
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests)
set_tests_properties(lolimot_sinus_noise PROPERTIES TIMEOUT 30)

add_test(NAME lolimot_temporel_noise 
         COMMAND $<TARGET_FILE:lolimot> ${CMAKE_SOURCE_DIR}/tests/Modele/Modele_Test_Temporel_Noise.txt
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests)
set_tests_properties(lolimot_temporel_noise PROPERTIES TIMEOUT 30)

add_test(NAME lolimot_cone_wo_noise 
         COMMAND $<TARGET_FILE:lolimot> ${CMAKE_SOURCE_DIR}/tests/Modele/Modele_Test_Cone_woNoise.txt
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests)
set_tests_properties(lolimot_cone_wo_noise PROPERTIES TIMEOUT 30)

add_test(NAME lolimot_sinus_wo_noise 
         COMMAND $<TARGET_FILE:lolimot> ${CMAKE_SOURCE_DIR}/tests/Modele/Modele_Test_Sinus_woNoise.txt
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests)
set_tests_properties(lolimot_sinus_wo_noise PROPERTIES TIMEOUT 30)

add_test(NAME lolimot_temporel_wo_noise 
         COMMAND $<TARGET_FILE:lolimot> ${CMAKE_SOURCE_DIR}/tests/Modele/Modele_Test_Temporel_woNoise.txt
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests)
set_tests_properties(lolimot_temporel_wo_noise PROPERTIES TIMEOUT 30)

add_test(NAME lolimot_plan_noise 
         COMMAND $<TARGET_FILE:lolimot> ${CMAKE_SOURCE_DIR}/tests/Modele/Modele_Test_Plan_Noise.txt
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests)
set_tests_properties(lolimot_plan_noise PROPERTIES TIMEOUT 30)

add_test(NAME lolimot_temporel_gap_noise 
         COMMAND $<TARGET_FILE:lolimot> ${CMAKE_SOURCE_DIR}/tests/Modele/Modele_Test_Temporel_Gap_Noise.txt
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests)
set_tests_properties(lolimot_temporel_gap_noise PROPERTIES TIMEOUT 30)

add_test(NAME lolimot_plan_wo_noise 
         COMMAND $<TARGET_FILE:lolimot> ${CMAKE_SOURCE_DIR}/tests/Modele/Modele_Test_Plan_woNoise.txt
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests)
set_tests_properties(lolimot_plan_wo_noise PROPERTIES TIMEOUT 30)

add_test(NAME lolimot_temporel_gap_wo_noise 
         COMMAND $<TARGET_FILE:lolimot> ${CMAKE_SOURCE_DIR}/tests/Modele/Modele_Test_Temporel_Gap_woNoise.txt
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests)
set_tests_properties(lolimot_temporel_gap_wo_noise PROPERTIES TIMEOUT 30)
