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

set(SMALL_LAPACK_SRCS Thirdparty/SmallLapack/d_abs.c
                      Thirdparty/SmallLapack/dbdsqr.c
		      Thirdparty/SmallLapack/dcopy.c
		      Thirdparty/SmallLapack/dgebd2.c
		      Thirdparty/SmallLapack/dgebrd.c
		      Thirdparty/SmallLapack/dgelq2.c
		      Thirdparty/SmallLapack/dgelqf.c
		      Thirdparty/SmallLapack/dgelss.c
		      Thirdparty/SmallLapack/dgemm.c
		      Thirdparty/SmallLapack/dgemv.c
		      Thirdparty/SmallLapack/dgeqr2.c
		      Thirdparty/SmallLapack/dgeqrf.c
		      Thirdparty/SmallLapack/dger.c
		      Thirdparty/SmallLapack/dlabad.c
		      Thirdparty/SmallLapack/dlabrd.c
		      Thirdparty/SmallLapack/dlacpy.c
		      Thirdparty/SmallLapack/dlamch.c
		      Thirdparty/SmallLapack/dlange.c
		      Thirdparty/SmallLapack/dlapy2.c
		      Thirdparty/SmallLapack/dlarfb.c
		      Thirdparty/SmallLapack/dlarf.c
		      Thirdparty/SmallLapack/dlarfg.c
		      Thirdparty/SmallLapack/dlarft.c
		      Thirdparty/SmallLapack/dlartg.c
		      Thirdparty/SmallLapack/dlas2.c
		      Thirdparty/SmallLapack/dlascl.c
		      Thirdparty/SmallLapack/dlaset.c
		      Thirdparty/SmallLapack/dlasq1.c
		      Thirdparty/SmallLapack/dlasq2.c
		      Thirdparty/SmallLapack/dlasq3.c
		      Thirdparty/SmallLapack/dlasq4.c
		      Thirdparty/SmallLapack/dlasq5.c
		      Thirdparty/SmallLapack/dlasq6.c
		      Thirdparty/SmallLapack/dlasr.c
		      Thirdparty/SmallLapack/dlasrt.c
		      Thirdparty/SmallLapack/dlassq.c
		      Thirdparty/SmallLapack/dlasv2.c
		      Thirdparty/SmallLapack/d_lg10.c
		      Thirdparty/SmallLapack/d_lg.c
		      Thirdparty/SmallLapack/dnrm2.c
		      Thirdparty/SmallLapack/dorg2r.c
		      Thirdparty/SmallLapack/dorgbr.c
		      Thirdparty/SmallLapack/dorgl2.c
		      Thirdparty/SmallLapack/dorglq.c
		      Thirdparty/SmallLapack/dorgqr.c
		      Thirdparty/SmallLapack/dorm2r.c
		      Thirdparty/SmallLapack/dormbr.c
		      Thirdparty/SmallLapack/dorml2.c
		      Thirdparty/SmallLapack/dormlq.c
		      Thirdparty/SmallLapack/dormqr.c
		      Thirdparty/SmallLapack/drot.c
		      Thirdparty/SmallLapack/drscl.c
		      Thirdparty/SmallLapack/dscal.c
		      Thirdparty/SmallLapack/d_sign.c
		      Thirdparty/SmallLapack/d_sqrt.c
		      Thirdparty/SmallLapack/dswap.c
		      Thirdparty/SmallLapack/dtrmm.c
		      Thirdparty/SmallLapack/dtrmv.c
		      Thirdparty/SmallLapack/exit_.c
		      Thirdparty/SmallLapack/F77_aloc.c
		      Thirdparty/SmallLapack/ieeeck.c
		      Thirdparty/SmallLapack/ilaenv.c
		      Thirdparty/SmallLapack/lsame.c
		      Thirdparty/SmallLapack/pow_dd.c
		      Thirdparty/SmallLapack/pow_di.c
		      Thirdparty/SmallLapack/s_cat.c
		      Thirdparty/SmallLapack/s_cmp.c
		      Thirdparty/SmallLapack/s_copy.c
		      Thirdparty/SmallLapack/xerbla.c)

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
set_tests_properties(lolimot_cone_noise PROPERTIES TIMEOUT 30)

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
